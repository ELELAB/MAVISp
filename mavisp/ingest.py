# MAVISp - classes for data ingest
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
import pandas as pd
import numpy as np
import logging as log
from datetime import date

class MAVISFileSystem:

    excluded_proteins = []
    supported_modes = ['basic_mode']
    supported_stability_methods = ['foldx5', 'rosetta_ref2015', 'rosetta_cartddg2020_ref2015']
    supported_interaction_methods = ['foldx5']

    def __init__(self, data_dir="/data/raw_data/computational_data/mavisp_data/"):
        
        log.info(f"initializing MAVISFileSystem from {data_dir}")

        self.data_dir = data_dir
        self.dataset_table = None
        self.mutation_datasets = None

        self.ingest_data()

    def ingest_data(self):

        self._tree = self._traverse(self.data_dir)
        self.dataset_table = self._dataset_table()
        self.mutation_datasets = self._mutation_tables()

    def _dir_list(self, d):

        return [ k for k,v in d.items() if v is not None ]
        """
        out = []
        for k,v in d.items():
            if v is not None:
                out.append(k)

        return out
        """
    def _file_list(self, d):

        out = []
        for k,v in d.items():
            if v is None:
                out.append(k)
        return out

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

    def _parse_mutation_list(self, fname):
        log.info(f"parsing mutation file {fname}")

        try:
            with open(fname) as fh:
                lines = fh.read().splitlines()
        except IOError:
            log.error("Couldn't parse mutation list {fname}")
            raise

        # remove duplicates, sort, remove empty lines
        lines = sorted(list(set(lines)), key=lambda x: int(x[1:-1]))
        lines = list(filter(lambda x: len(x) != 0, lines))

        log.debug(f"found mutations: {lines}")

        return lines

    def _parse_foldx_summary(self, fname, type='STABILITY', version='FoldX5', unit='kcal/mol'):

        log.info(f"parsing FoldX summary file {fname}")

        try:
            data = pd.read_csv(fname, delim_whitespace=True, header=None)
        except IOError:
            log.error("Couldn't parse FoldX summary file {fname}")
            raise

        # remove chain ID from identifier, keep mutation and DDG average, rename cols
        data[0] = data[0].apply(lambda x: x[0] + x[2:])

        data = data[[0,1]].rename(columns={0:'mutations',
                                           1:f"{type} ({version}, {unit})"}).set_index('mutations')

        log.debug(f"collected data: {data}")
        return(data)

    def _parse_rosetta_aggregate(self, fname):

        log.info(f"parsing Rosetta aggregate file {fname}")

        try:
            mutation_data =  pd.read_csv(fname)
        except IOError:
            log.error("Couldn't parse Rosetta energy file {fname}")
            raise

        energy = mutation_data[mutation_data['state'] == 'ddg']['total_score'].to_numpy()[0]

        log.debug("collected Rosetta energy: {energy}")

        return energy

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
                raise

        selected_file = dates[max(dates.keys())]
        
        log.debug(f"file names and their dates {dates}")
        log.debug(f"selected most recent file {selected_file} among {fnames}")

        return selected_file

    def _dataset_table(self):
        
        log.info("generating dataset table")

        df_list = {}

        k = 0

        for system in  self._dir_list(self._tree):
            if system in self.excluded_proteins:
                log.warning(f"ignoring unsupported protein {system}")
                continue
            for mode in self._dir_list(self._tree[system]):
                if mode not in self.supported_modes:
                    log.warning(f"ignoring unsupported mode {mode} in {[system, mode]}")
                    continue
                log.info(f"adding {[system, mode]} to dataset")

                try: 
                    mutation_list = self._mutation_list(system, mode)
                except:
                    exit(1)
                    
                df_list[k] = [system, mode, mutation_list]

                k += 1

        main_df = pd.DataFrame.from_dict(df_list, orient='index', columns=['system', 'mode', 'mutations'])
        log.debug(f"identified datasets:\n{main_df}")

        return main_df

    def _mutation_list(self, system, mode):
            
        log.info(f"Gathering mutation list for {system} {mode}")

        mutation_files = self._file_list(self._tree[system][mode]['mutation_list'])
        most_recent_mut_file = self._select_most_recent_file(mutation_files)
        log.info(f"selected {most_recent_mut_file} as mutation file")

        mut_path = os.path.join(self.data_dir, system, mode, 'mutation_list', most_recent_mut_file)

        try:
            mutations = self._parse_mutation_list(mut_path)
        except IOError:
            log.error(f"Couldn't parse mutation list {mut_path}")
            exit(1)

        log.debug(f"final mutation lists: {mutations}")
        
        return mutations

    def _mutation_tables(self):
        if self.dataset_table is None:
            return None

        data_dfs = {}

        for idx, r in self.dataset_table.iterrows():

            system = r['system']
            mode = r['mode']
            mutations = r['mutations']
            
            log.info(f"Gathering data for {r['system']} {r['mode']}")
            
            this_df = pd.DataFrame({'Mutation': mutations})
            this_df = this_df.set_index('Mutation')
            
            analysis_basepath = os.path.join(self.data_dir, system, mode)

            if 'stability' in self._dir_list(self._tree[system][mode]):

                tmp = self._dir_list(self._tree[system][mode]['stability'])
                if len(tmp) != 1:
                    log.error("stability folder has to contain only 1 dataset. It will be skipped")
                    continue
                structure_ID, residue_range = tmp[0].split("_")

                tmp = self._dir_list(self._tree[system][mode]['stability'][f'{structure_ID}_{residue_range}'])
                if len(tmp) != 1:
                    log.error("stability folder has to contain only 1 dataset. It will be skipped")
                    continue
                method = tmp[0]

                tmp = self._dir_list(self._tree[system][mode]['stability'][f'{structure_ID}_{residue_range}'][method])
                if len(tmp) != 1:
                    log.error("stability folder has to contain only 1 dataset. It will be skipped")
                    continue
                model = tmp[0]

                stability_methods = self._dir_list(self._tree[system][mode]['stability'][f'{structure_ID}_{residue_range}'][method][model])

                sm_basepath = os.path.join(analysis_basepath, 'stability', f'{structure_ID}_{residue_range}', method, model)

                log.info(f"found methods for stability: {stability_methods}")

                for sm in stability_methods:

                    if sm not in self.supported_stability_methods:
                        log.warning(f"WARNING: stability method {sm} not supported, will be skipped")
                        continue

                    if sm == 'foldx5':

                        log.info("parsing data for foldx5")

                        foldx_dir = os.path.join(sm_basepath, 'foldx5')

                        foldx_files = os.listdir(foldx_dir)
                        if len(foldx_files) != 1:
                            log.error(f"multiple files found in {foldx_dir}; only one expected")
                            exit(1)
                        foldx_file = foldx_files[0]

                        try:
                            data = self._parse_foldx_summary(os.path.join(foldx_dir, foldx_file))
                        except IOError:
                            exit(1)

                        this_df = this_df.join(data)

                        log.debug(f"adding foldx5 data {this_df}")

                    if sm == 'rosetta_ref2015' or sm == 'rosetta_cartddg2020_ref2015':

                        log.info(f"parsing data for {sm}")

                        data = []

                        rosetta_dir = os.path.join(sm_basepath, sm)

                        for mutation in this_df.index:
                            rosetta_file = os.path.join(rosetta_dir, f"{mutation}_aggregate.csv")
                            try:
                                data.append(self._parse_rosetta_aggregate(rosetta_file))
                            except IOError:
                                log.error("couldn't open expected energy file {rosetta_file}")
                                exit(1)

                        log.debug(f"adding {sm} data {data}")
                        this_df['STABILITY (Rosetta, ref2015, kcal/mol)'] = data

            if 'local_interactions' in self._dir_list(self._tree[system][mode]):

                interaction_methods = self._dir_list(self._tree[system][mode]['local_interactions'])

                log.info(f"found methods for interaction: {stability_methods}")

                for method in interaction_methods:

                    if method not in self.supported_interaction_methods:
                        log.warning(f"Method {method} for interaction is not supported and it will be skipped")
                        continue
                    if method == 'foldx5':

                        log.info("parsing data for foldx5")

                        int_basepath = os.path.join(analysis_basepath, 'local_interactions', 'foldx5')

                        interactors = self._dir_list(self._tree[system][mode]['local_interactions']['foldx5'])
                        if len(interactors) == 0:
                            log.error("zero interactors found for FoldX local interactions")
                            exit(1)

                        for interactor in interactors:

                            interactor_dir = os.path.join(int_basepath, interactor)

                            foldx_files = os.listdir(interactor_dir)

                            if len(foldx_files) != 1:
                                log.error(f"zero or multiple files found in {foldx_dir}; exactly one expected")
                                exit(1)
                            foldx_file = foldx_files[0]

                            try:
                                data = self._parse_foldx_summary(os.path.join(interactor_dir, foldx_file), type="LOCAL INT", version=f"Binding with {interactor}, Foldx5")
                            except IOError:
                                exit(1)

                            this_df = this_df.join(data)

                            log.debug(f"adding foldx5 data {this_df}")

            this_df = this_df.reset_index()
            data_dfs[(system, mode)] = this_df

        return data_dfs

