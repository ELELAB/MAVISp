# MAVISp - classes for data ingest and manipulation
# Copyright (C) 2022 Matteo Tiberti, Danish Cancer Society
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
from functools import reduce
from collections import defaultdict
import pandas as pd
import numpy as np
import logging as log
from datetime import date
from mavisp.error import *
from mavisp.methods import *
from mavisp.modules import *
from mavisp.modes import *
import yaml
from termcolor import colored
from tabulate import tabulate

class MAVISpFileSystem:

    supported_modes = { 'simple_mode'   : MAVISpSimpleMode(),
                        'ensemble_mode' : MAVISpEnsembleMode()}

    def __init__(self, modes=None, include_proteins=None, exclude_proteins=None, data_dir="database", verbose=True):

        self.log = log.getLogger(__name__)

        if verbose:
            level = log.INFO
        else:
            level = log.WARNING

        self.log.setLevel(level)

        if modes is None:
            modes = self.supported_modes

        modes_diff = set(modes).difference(set(self.supported_modes.keys()))
        if len(modes_diff) > 0:
            raise TypeError(f"the following modes are not supported: {modes_diff}")

        self.data_dir = data_dir

        self._tree = self._traverse(self.data_dir)

        self.dataset_tables = {}
        for mode_name, mode in self.supported_modes.items():
            self.dataset_tables[mode_name] = self._gen_dataset_table(mode, include=include_proteins, exclude=exclude_proteins)

    def _gen_dataset_table(self, mode, include, exclude):

        self.log.info("generating dataset table")

        df_list = []

        for system in self._dir_list(self._tree):
            if (include is not None and system not in include) or (exclude is not None and system in exclude) and (not (exclude is None and include is None)):
                self.log.warning(f"ignoring not selected protein {system}")
                continue

            self.log.info(f"adding {[system, mode]} to dataset")
            mavisp_criticals = []

            try:
                mutation_list = self._parse_mutation_list(system, mode.name)
            except Exception as e:
                self.log.error(e)
                mutation_list = None
                mavisp_criticals.append(MAVISpCriticalError("the mutation list was not available, readable or in the expected format"))

            try:
                metadata, metadata_mavisp_criticals = mode.parse_metadata(self.data_dir, system)
            except IOError:
                self.log.error("Couldn't parse metadata file")
                metadata = {k : None for k in mode.supported_metadata}
                metadata_mavisp_criticals = []

            mavisp_criticals.extend(metadata_mavisp_criticals)

            metadata['system'] = system
            metadata['mutations'] = mutation_list
            metadata['mode'] = mode.name
            metadata['mavisp_criticals'] = mavisp_criticals

            df_list.append(metadata)

        main_df = pd.DataFrame.from_records(df_list)
        self.log.debug(f"identified datasets:\n{main_df}")

        return main_df

    def _dir_list(self, d):

        return [ k for k,v in d.items() if v is not None ]

    def _file_list(self, d):

        return [ k for k,v in d.items() if v is None ]

    def _traverse(self, rootdir):

        self.log.info(f"building directory tree for {rootdir}")

        rootdir = str(rootdir)

        tree = {}
        rootdir = rootdir.rstrip(os.sep)
        start = rootdir.rfind(os.sep) + 1
        for path, dirs, files in os.walk(rootdir, followlinks=True):
            folders = path[start:].split(os.sep)
            subdir = dict.fromkeys(files)
            parent = reduce(dict.get, folders[:-1], tree)
            parent[folders[-1]] = subdir

        return tree[os.path.basename(os.path.normpath(rootdir))]

    def _parse_mutation_list(self, system, mode):

        self.log.info(f"Gathering mutation list for {system} {mode}")

        mutation_files = self._file_list(self._tree[system][mode]['mutation_list'])
        most_recent_mut_file = self._select_most_recent_file(mutation_files)
        self.log.info(f"selected {most_recent_mut_file} as mutation file")

        mut_path = os.path.join(self.data_dir, system, mode, 'mutation_list', most_recent_mut_file)

        try:
            mutations = pd.read_csv(mut_path, delim_whitespace=True)
        except Exception as e:
            raise TypeError(f"Couldn't parse mutation list {mut_path} with error {e}")

        # remove duplicates, sort
        mutations = mutations.drop_duplicates(ignore_index=True)
        if len(mutations['mutation']) != mutations['mutation'].unique().shape[0]:
            raise TypeError(f"{system}, {mode}: Duplicate mutations with inconsistent references in the mutation list")

        mutations['res_num'] = mutations['mutation'].str[1:-1].astype(int)
        mutations['alt'] = mutations['mutation'].str[-1]
        mutations = mutations.sort_values(['res_num', 'alt'])
        mutations = mutations.drop(columns=['res_num', 'alt'])

        self.log.debug(f"found mutations: {mutations['mutation'].tolist()}")

        return mutations

    def _select_most_recent_file(self, fnames):

        dates = {}

        for fname in fnames:

            basename = os.path.splitext(fname)[0]

            try:
                dates[date(int(basename[-4:]),
                           int(basename[-6:-4]),
                           int(basename[-8:-6]))] = fname
            except (ValueError, TypeError):
                log.error("file {fname} doesn't contain a valid date string at the end of the file name")
                raise TypeError

        selected_file = dates[max(dates.keys())]

        self.log.debug(f"file names and their dates {dates}")
        self.log.debug(f"selected most recent file {selected_file} among {fnames}")

        return selected_file

    def ingest(self, stop_at='critical'):

        for mode_name, mode in self.supported_modes.items():

            # for every protein (and supported mode) in the dataset
            mavisp_warnings_column = []
            mavisp_errors_column = []
            mavisp_criticals_column = []
            mavisp_dataset_column = []

            for _, r in self.dataset_tables[mode_name].iterrows():

                mavisp_modules = defaultdict(lambda: None)
                mavisp_warnings = defaultdict(list)
                mavisp_errors = defaultdict(list)
                mavisp_criticals = r['mavisp_criticals']

                system = r['system']
                mutations = r['mutations']

                if len(mavisp_criticals) > 0:
                    mavisp_dataset_column.append(mavisp_modules)
                    mavisp_errors_column.append(mavisp_errors)
                    mavisp_warnings_column.append(mavisp_warnings)
                    mavisp_criticals_column.append(mavisp_criticals)

                    continue

                self.log.info(f"Gathering data for {r['system']} {mode_name}")

                analysis_basepath = os.path.join(self.data_dir, system, mode_name)

                # for every available module:
                for mod in mode.supported_modules:

                    # check if the dataset is available
                    if mod.module_dir in self._dir_list(self._tree[system][mode_name]):
                        try:
                            self.log.info(f"processing module {mod.name} for {system}, {mode_name}")
                            this_module = mod(analysis_basepath)
                            this_module.ingest(mutations['mutation'].tolist())

                        except MAVISpMultipleError as e:
                            mavisp_errors[mod.name].extend(e.critical)
                            mavisp_warnings[mod.name].extend(e.warning)
                            if len(e.critical) != 0:
                                mavisp_modules[mod.name] = None
                                continue
                            if len(e.warning) != 0 and stop_at == 'warning':
                                mavisp_modules[mod.name] = None
                                continue

                        mavisp_modules[mod.name] = this_module
                mavisp_dataset_column.append(mavisp_modules)
                mavisp_errors_column.append(mavisp_errors)
                mavisp_warnings_column.append(mavisp_warnings)
                mavisp_criticals_column.append(mavisp_criticals)

            self.dataset_tables[mode_name]['modules'] = mavisp_dataset_column
            self.dataset_tables[mode_name]['criticals'] = mavisp_criticals_column
            self.dataset_tables[mode_name]['errors'] = mavisp_errors_column
            self.dataset_tables[mode_name]['warnings'] = mavisp_warnings_column

    def get_datasets_table_view(self, mode_name):
        return self.dataset_tables[mode_name][self.supported_modes[mode_name].index_cols]

    def get_annotation_tables_view(self, mode_name):

        all_tables = {}

        for _, r in self.dataset_tables[mode_name].iterrows():

            annotation_table = pd.DataFrame(index=r['mutations'])
            for mod in r['modules']:
                annotation_table = annotation_table.join(mod.get_dataset_view())

                all_tables[f"{r['system']}_{mode_name}"] = annotation_table

        return all_tables

    def get_datasets_table_summary(self, mode_name):

        data = defaultdict(list)

        for _, r in self.dataset_tables[mode_name].iterrows():
            data['System'].append(r['system'])
            data['Mode'].append(r['mode'])
            if len(r['criticals']) > 0:
                data['Status'].append(colored("CRITICAL", 'magenta'))
            elif sum( [ len(x) for x in list(r['errors'].values())]) > 0:
                data['Status'].append(colored("ERROR", 'red'))
            elif sum([ len(x) for x in list(r['warnings'].values())]) > 0:
                data['Status'].append(colored("WARNING", 'yellow'))
            else:
                data['Status'].append(colored("OK", 'green'))

        return pd.DataFrame(data)

    def get_datasets_table_details(self, mode_name):

        data = defaultdict(list)

        for _, r in self.dataset_tables[mode_name].iterrows():

            if len(r['criticals']) > 0:
                data['system'].append(r['system'])
                data['mode'].append(r['mode'])
                data['module'].append(pd.NA)
                data['details_crit'].append(r['criticals'])
                data['details_warn'].append(list())
                data['details_err'].append(list())
                data['status'].append('critical')
                continue

            for this_m in self.supported_modes[r['mode']].supported_modules:
                data['system'].append(r['system'])
                data['mode'].append(r['mode'])
                data['module'].append(this_m.name)

                # module not ran/available
                if this_m.name not in r['modules'].keys():
                    data['status'].append("not_available")
                    data['details_crit'].append(list())
                    data['details_warn'].append(list())
                    data['details_err'].append(list())
                    continue

                # all good
                if len(r['errors'][this_m.name]) == 0 and len(r['warnings'][this_m.name]) == 0:
                    data['status'].append("ok")
                    data['details_crit'].append(list())
                    data['details_warn'].append(list())
                    data['details_err'].append(list())
                    continue

                # errors
                if len(r['errors'][this_m.name]) > 0:
                    data['status'].append("error")
                    data['details_crit'].append(list())
                    data['details_err'].append([str(x).strip() for x in r['errors'][this_m.name]])
                    data['details_warn'].append([str(x).strip() for x in r['warnings'][this_m.name]])
                    continue

                # warnings
                if len(r['warnings'][this_m.name]) > 0:
                    data['status'].append("warning")
                    data['details_crit'].append(list())
                    data['details_err'].append(list())
                    data['details_warn'].append([str(x).strip() for x in r['warnings'][this_m.name]])

        return pd.DataFrame(data)
