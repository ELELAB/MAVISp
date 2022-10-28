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
import yaml
from termcolor import colored
from tabulate import tabulate

class DataType(object):
    def __init__(self, data_dir=None, stop_at='critical'):

        self.data_dir = data_dir
        
        if data_dir is None:
            return

    def ingest(self, stop_at='critical'):
        pass

    def process(self):
        pass

    @property
    def data_dir(self):
        return self._data_dir

    @data_dir.setter
    def data_dir(self, val):
        if val is None:
            self._data_dir = val
            return
        if os.path.exists(val) and os.path.isdir(val):
            self._data_dir = val
        else:
            raise MAVISpMultipleError(critical=[MAVISpCriticalError("the input directory pathway doesn't exist or is not a directory")],
                                      warning=[])

    @data_dir.getter
    def data_dir(self):
        return self._data_dir


class CancermutsTable(DataType):

    module_dir = "cancermuts"
    name = "cancermuts"

    def __init__(self, data_dir=None):

        super().__init__(data_dir)

    def ingest(self, mutations):
        warnings = []

        # check that we have only one file
        cancermuts_files = os.listdir(os.path.join(self.data_dir, self.module_dir))
        if len(cancermuts_files) != 1:
            this_error = f"multiple files found in {cancermuts_files}; only one expected"
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[MAVISpCriticalError(this_error)])
        
        cancermuts_file = cancermuts_files[0]

        log.info(f"parsing Cancermuts file {cancermuts_file}")

        # parse cancermuts table
        try:
            cancermuts = pd.read_csv(os.path.join(self.data_dir, self.module_dir, cancermuts_file))
        except:
            this_error = f"Failed parsing Cancermuts table {cancermuts_file}"
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[MAVISpCriticalError(this_error)])

        # check if all the required columns are present
        required_columns = ['aa_position', 'ref_aa', 'alt_aa', 'gnomad_genome_af', 'gnomad_exome_af', 'REVEL_score', 'sources']
        if not set(required_columns).issubset(cancermuts.columns):
            this_error = f"input table doesn't have all the required columns"
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[MAVISpCriticalError(this_error)])

        # process table
        cancermuts = cancermuts[ ~ pd.isna(cancermuts.alt_aa)]
        cancermuts['mutation_index'] = cancermuts.ref_aa + cancermuts.aa_position.astype(str) + cancermuts.alt_aa
        cancermuts = cancermuts.set_index('mutation_index')

        # check if all mutations are present in the cancermuts table
        available_mutations = cancermuts.index.intersection(mutations)
        if len(available_mutations) != len(mutations):
            warnings.append(MAVISpWarningError("Not all the annotated mutations were found in the Cancermuts table"))

        # keep rows with mutations only
        cancermuts = cancermuts.loc[available_mutations]

        # check if all the mutations have a defined REVEL score
        if np.any(pd.isna(cancermuts['REVEL_score'])):
            warnings.append(MAVISpWarningError("One or more mutations in the Cancermuts table don't have an associated REVEL score"))

        # filter by column and pretty rename column names
        cancermuts = cancermuts[['gnomad_genome_af', 'gnomad_exome_af', 'REVEL_score', 'sources']]

        self.data = cancermuts.rename(columns={ 'gnomad_genome_af' : 'gnomAD genome allele frequency',
                                                'gnomad_exome_af'  : 'gnomAD exome allele frequency',
                                                'REVEL_score'      : 'REVEL score',
                                                'sources'          : 'Mutation sources' })

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[])

    def get_dataset_view(self):
        return self.data


class MAVISpFileSystem:

    supported_mdoules = ['cancermuts']
    supported_modes = ['simple_mode']
    supported_stability_methods = ['foldx5', 'rosetta_ref2015', 'rosetta_cartddg2020_ref2015']
    supported_interaction_methods = ['foldx5']

    supported_modules = [ CancermutsTable ]

    def __init__(self, modes=None, proteins=None, data_dir="database", verbose=True):

        self.log = log.getLogger(__name__)

        if verbose:
            level = log.INFO
        else:
            level = log.WARNING

        self.log.setLevel(level)

        self.proteins = proteins
        if modes is None:
            modes = self.supported_modes            
        
        modes_diff = set(modes).difference(set(self.supported_modes))
        if len(modes_diff) > 0: 
            raise TypeError(f"the following modes are not supported: {modes_diff}")

        self.data_dir=data_dir

        self._tree = self._traverse(self.data_dir)
        self.dataset_table = self._dataset_table()

    def _dataset_table(self):

        warnings = []

        self.log.info("generating dataset table")

        df_list = []

        for system in  self._dir_list(self._tree):
            if self.proteins is not None and system not in self.proteins:
                self.log.warning(f"ignoring unsupported protein {system}")
                warnings.append(f"ignoring unsupported protein {system}")
                continue
            for mode in self._dir_list(self._tree[system]):
                if mode not in self.supported_modes:
                    self.log.warning(f"ignoring unsupported mode {mode} in {[system, mode]}")
                    warnings.append(f"ignoring unsupported mode {mode} in {[system, mode]}")
                    continue
                self.log.info(f"adding {[system, mode]} to dataset")

                try: 
                    mutation_list = self._parse_mutation_list(system, mode)
                except:
                    self.log.error(f"Couldn't parse mutation list table")

                try:
                    metadata = self._parse_metadata(system, mode)
                except IOError:
                    log.error("Couldn't parse metadata file")
                    raise IOError

                curators= ', '.join(
                    [ f"{curator} ({', '.join(metadata['curators'][curator]['affiliation'])})" for curator in metadata['curators'].keys() ]
                    )

                df_list.append((system, mode, mutation_list, curators))

        main_df = pd.DataFrame.from_records(df_list, columns=['system', 'mode', 'mutations', 'curators'])
        log.debug(f"identified datasets:\n{main_df}")

        return main_df

    def _dir_list(self, d):

        return [ k for k,v in d.items() if v is not None ]

    def _file_list(self, d):

        return [ k for k,v in d.items() if v is None ]

    def _traverse(self, rootdir):

        log.info(f"building directory tree for {rootdir}")

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
            
        log.info(f"Gathering mutation list for {system} {mode}")

        mutation_files = self._file_list(self._tree[system][mode]['mutation_list'])
        most_recent_mut_file = self._select_most_recent_file(mutation_files)
        log.info(f"selected {most_recent_mut_file} as mutation file")

        mut_path = os.path.join(self.data_dir, system, mode, 'mutation_list', most_recent_mut_file)

        try:
            with open(mut_path) as fh:
                lines = fh.read().splitlines()
        except IOError:
            log.error("Couldn't parse mutation list {mut_path}")
            raise IOError

        # remove duplicates, sort, remove empty lines
        lines = sorted(list(set(lines)), key=lambda x: int(x[1:-1]))
        mutations = list(filter(lambda x: len(x) != 0, lines))

        log.debug(f"found mutations: {mutations}")
        
        return mutations

    def _parse_metadata(self, system, mode):
        log.info("parsing metadata file")

        metadata_path = os.path.join(self.data_dir, system, mode, 'metadata.yaml')

        with open(metadata_path) as fh:
            return yaml.safe_load(fh)

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
        
        log.debug(f"file names and their dates {dates}")
        log.debug(f"selected most recent file {selected_file} among {fnames}")

        return selected_file

    def ingest(self, stop_at='critical'):

        # for every protein (and supported mode) in the dataset
        mavisp_warnings_column = []
        mavisp_errors_column = []
        mavisp_dataset_column = []

        for _, r in self.dataset_table.iterrows():

            system = r['system']
            mode = r['mode']
            mutations = r['mutations']

            mavisp_modules = []
            mavisp_warnings = defaultdict(list)
            mavisp_errors = defaultdict(list)
            
            log.info(f"Gathering data for {r['system']} {r['mode']}")
            
            this_df = pd.DataFrame({'Mutation': mutations})
            this_df = this_df.set_index('Mutation')
            
            analysis_basepath = os.path.join(self.data_dir, system, mode)

            # for every available module:
            for mod in self.supported_modules:
                
                # check if the dataset is available
                if mod.module_dir in self._dir_list(self._tree[system][mode]):
                    try:
                        this_module = mod(analysis_basepath)
                        this_module.ingest(mutations)
                        
                    except MAVISpMultipleError as e:
                        mavisp_errors[mod.name].extend(e.critical)
                        mavisp_warnings[mod.name].extend(e.warning)
                        if len(e.critical) != 0:
                            continue
                        if len(e.warning) != 0 and stop_at == 'warning':
                            continue
                    
                    mavisp_modules.append(this_module)
            mavisp_dataset_column.append(mavisp_modules)
            mavisp_errors_column.append(mavisp_errors)
            mavisp_warnings_column.append(mavisp_warnings)
        
        self.dataset_table['modules'] = mavisp_dataset_column
        self.dataset_table['errors'] = mavisp_errors_column
        self.dataset_table['warnings'] = mavisp_warnings_column

    def get_datasets_table_view(self):
        return self.dataset_table[['system', 'mode', 'curators']]

    def get_annotation_tables_view(self):
        
        all_tables = {}

        for _, r in self.dataset_table.iterrows():

            annotation_table = pd.DataFrame(index=r['mutations'])
            for mod in r['modules']:
                annotation_table = annotation_table.join(mod.get_dataset_view())

                all_tables[f"{r['system']}_{r['mode']}"] = annotation_table

        return all_tables

    def get_datasets_table_summary(self):

        text = ""

        data = defaultdict(list)

        for _, r in self.dataset_table.iterrows():
            data['System'].append(r['system'])
            data['Mode'].append(r['mode'])
            if sum( [ len(x) for x in list(r['errors'].values())]) > 0:
                data['Status'].append(colored("ERROR", 'red'))
            elif sum([ len(x) for x in list(r['warnings'].values())]) > 0:
                data['Status'].append(colored("WARNING", 'yellow'))
            else: 
                data['Status'].append(colored("OK", 'green'))

        return pd.DataFrame(data)
