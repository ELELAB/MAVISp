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
    supported_interaction_methods = ['foldx5', 'rosetta_flexddg_talaris2014']

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

    def _parse_mutation_list(self, fname):
        log.info(f"parsing mutation file {fname}")

        try:
            with open(fname) as fh:
                lines = fh.read().splitlines()
        except IOError:
            log.error("Couldn't parse mutation list {fname}")
            raise IOError

        # remove duplicates, sort, remove empty lines
        lines = sorted(list(set(lines)), key=lambda x: int(x[1:-1]))
        lines = list(filter(lambda x: len(x) != 0, lines))

        log.debug(f"found mutations: {lines}")

        return lines

    def _parse_foldx_csv(self, fname, type='STABILITY', version='FoldX5', unit='kcal/mol'):

        log.info(f"parsing FoldX csv file {fname}")

        try:
            df = pd.read_csv(fname)
        except IOError:
            log.error("Couldn't parse FoldX summary file {fname}")
            raise IOError

        # create residue column
        df['residue'] = df['WT residue type'] + df['Residue #'].astype(str)

        # drop column we don't need
        df = df.drop(['WT residue type', 'Residue #', 'chain ID'], axis=1)

        # stack remaining columns
        df = df.set_index('residue')
        df = df.stack()
        df = df.reset_index()

        # create mutation column
        df['mutations'] = df['residue'] + df['level_1']
        df = df.set_index('mutations')

        # drop now useless columns, rename
        df = df.drop(['residue', 'level_1'], axis=1)
        df = df.rename(columns={0 : f"{type} ({version}, {unit})"})

        log.debug(f"collected data: {df}")
        return(df)

    def _parse_rosetta_aggregate(self, fname, type='STABILITY', version='Rosetta Flexddg2020', unit='kcal/mol'):

        log.info(f"parsing Rosetta aggregate file {fname}")

        try:
            mutation_data = pd.read_csv(fname)
        except IOError:
            log.error(f"Couldn't parse Rosetta energy file {fname}")
            raise IOError

        mutation_data = mutation_data[mutation_data['state'] == 'ddg']
        mutation_data = mutation_data.set_index('mutation_label')
        mutation_data = mutation_data[['total_score']]
        mutation_data = mutation_data.rename(columns={'total_score':f'{type} ({version}, {unit})'})

        return mutation_data

    def _parse_cancermuts(self, fname):

        log.info(f"parsing Cancermuts file {fname}")

        try:
            cancermuts = pd.read_csv(fname)
        except IOError:
            log.error(f"Couldn't parse Cancermuts file {fname}")
            raise

        cancermuts = cancermuts[ ~ pd.isna(cancermuts.alt_aa)]
        cancermuts['mutation_index'] = cancermuts.ref_aa + cancermuts.aa_position.astype(str) + cancermuts.alt_aa
        cancermuts = cancermuts.set_index('mutation_index')
        cancermuts = cancermuts[['gnomad_genome_af', 'gnomad_exome_af', 'REVEL_score', 'sources']]

        return       cancermuts.rename(columns={ 'gnomad_genome_af' : 'gnomAD genome allele frequency',
                                                 'gnomad_exome_af'  : 'gnomAD exome allele frequency',
                                                 'REVEL_score'      : 'REVEL score',
                                                 'sources'          : 'Mutation sources' })

    def _parse_pmid(self, fname):

        log.info(f"parsing PMID file {fname}")

        try:
            pmid = pd.read_csv(fname, delim_whitespace=True)
        except IOError:
            log.error(f"Couldn't parse PMID file {fname}")
            raise

        return pmid.set_index('mutation')

    def _parse_ptm(self, fname):

        try:
            ptms = pd.read_csv(fname, delim_whitespace=True)
        except:
            log.error(f"file {fname} not in the right format")
            raise TypeError

        ptms['mutation'] = os.path.splitext(os.path.basename(fname))[0]

        if len(ptms) != 1:
            tmp_data = {}
            for c in ptms.columns:
                tmp_data[c] = ", ".join(map(str, ptms[c].to_list()))
            ptms = pd.DataFrame(tmp_data)

        ptms = ptms.set_index('mutation')

        return ptms

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

    def _dataset_table(self):
        
        log.info("generating dataset table")

        df_list = []

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
                    
                df_list.append((system, mode, mutation_list))

        main_df = pd.DataFrame.from_records(df_list, columns=['system', 'mode', 'mutations'])
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
            raise IOError

        log.debug(f"final mutation lists: {mutations}")
        
        return mutations

    def _process_stability(self, row):
        keys = [ k for k in row.keys() if k.startswith('STABILITY') ]

        if len(keys) == 2:
            if   'Rosetta' in keys[0] and 'FoldX' in keys[1]:
                rosetta_header, foldx_header    = keys
            elif 'Rosetta' in keys[1] and 'FoldX' in keys[0]:
                foldx_header,   rosetta_header  = keys
            else:
                raise TypeError

        stab_co = 3.0
        neut_co = 2.0

        if rosetta_header not in row.index or foldx_header not in row.index:
            return pd.NA

        if row[foldx_header] > stab_co and row[rosetta_header] > stab_co:
            return 'Destabilizing'
        if row[foldx_header] < (- stab_co) and row[rosetta_header] < (- stab_co):
            return 'Stabilizing'
        if (- neut_co) < row[foldx_header] < neut_co and (- neut_co) < row[rosetta_header] < neut_co:
            return('Neutral')
        return 'Uncertain'

    def _process_table(self, table, which='all'):

        log.info("Processing metatable")

        functions = {'Stability classification' : self._process_stability}

        if which == 'all':
            this_run = functions.keys()
        else:
            this_run = list( set(which).intersection(set(list(functions.keys()))) )

        for colname in this_run:
            log.info(f"Processing {colname}")
            table[colname] = table.apply(functions[colname], axis=1)

        return table

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
                structure_ID, residue_range = tmp[0].split("_", maxsplit=1)

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
                            data = self._parse_foldx_csv(os.path.join(foldx_dir, foldx_file), version='FoldX5')
                        except IOError:
                            exit(1)

                        this_df = this_df.join(data)

                        log.debug(f"adding foldx5 data {this_df}")

                    if sm == 'rosetta_ref2015' or sm == 'rosetta_cartddg2020_ref2015':

                        log.info(f"parsing data for {sm}")

                        data = []

                        rosetta_dir = os.path.join(sm_basepath, sm)

                        rosetta_files = os.listdir(rosetta_dir)
                        if len(rosetta_files) != 1:
                            log.error(f"multiple files found in {rosetta_dir}; only one expected")
                            exit(1)
                        rosetta_file = rosetta_files[0]

                        try:
                            data = self._parse_rosetta_aggregate(os.path.join(rosetta_dir, rosetta_file))
                        except IOError:
                            log.error("couldn't open expected energy file {rosetta_file}")
                            exit(1)

                        this_df = this_df.join(data)

                        log.debug(f"adding {sm} data")

            this_df = self._process_table(this_df, which=['Stability classification'])

            if 'local_interactions' in self._dir_list(self._tree[system][mode]):

                interaction_methods = self._dir_list(self._tree[system][mode]['local_interactions'])

                log.info(f"found methods for interaction: {interaction_methods}")

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
                                data = self._parse_foldx_csv(os.path.join(interactor_dir, foldx_file), type="LOCAL INT", version=f"Binding with {interactor}, Foldx5")
                            except IOError:
                                exit(1)

                            this_df = this_df.join(data)

                            log.debug(f"adding foldx5 data {this_df}")

                    if method == 'rosetta_flexddg_talaris2014':

                        log.info("parsing data for rosetta_flexddg_talaris2014")

                        int_basepath = os.path.join(analysis_basepath, 'local_interactions', 'rosetta_flexddg_talaris2014')

                        interactors = self._dir_list(self._tree[system][mode]['local_interactions']['rosetta_flexddg_talaris2014'])
                        if len(interactors) == 0:
                            log.error("zero interactors found for Rosetta FlexDDG local interactions")
                            exit(1)

                        for interactor in interactors:

                            interactor_dir = os.path.join(int_basepath, interactor)

                            rosetta_files = os.listdir(interactor_dir)

                            if len(rosetta_files) != 1:
                                log.error(f"zero or multiple files found in {interactor_dir}; exactly one expected")
                                exit(1)
                            rosetta_file = rosetta_files[0]

                            try:
                                data = self._parse_rosetta_aggregate(os.path.join(interactor_dir, rosetta_file), type="LOCAL INT", version=f"Binding with {interactor}, Rosetta Talaris 2014")
                            except IOError:
                                exit(1)

                            this_df = this_df.join(data)

                            log.info(f"Adding local interaction information for {method}, {interactor}")

                            log.debug(f"adding rosetta_flexddg_talaris2014 data {this_df}")

            if 'cancermuts' in self._dir_list(self._tree[system][mode]):

                cancermuts_dir = os.path.join(analysis_basepath, 'cancermuts')

                cancermuts_files = os.listdir(cancermuts_dir)
                if len(cancermuts_files) != 1:
                    log.error(f"multiple files found in {cancermuts_files}; only one expected")
                    exit(1)
                cancermuts_file = cancermuts_files[0]

                try:
                    cancermuts_data = self._parse_cancermuts(os.path.join(analysis_basepath, 'cancermuts', cancermuts_file))
                except IOError:
                    exit(1)

                this_df = this_df.join(cancermuts_data)

            if 'pmid_list' in self._dir_list(self._tree[system][mode]):

                pmid_dir = os.path.join(analysis_basepath, 'pmid_list')

                pmid_files = os.listdir(pmid_dir)
                if len(pmid_files) != 1:
                    log.error(f"multiple or no files found in {pmid_files}; only one expected")
                    exit(1)
                pmid_file = pmid_files[0]

                try:
                    pmid_data = self._parse_pmid(os.path.join(analysis_basepath, 'pmid_list', pmid_file))
                except IOError:
                    exit(1)

                this_df = this_df.join(pmid_data)

            if 'ptm' in self._dir_list(self._tree[system][mode]):

                ptm_dir = os.path.join(analysis_basepath, 'ptm')

                ptm_files = os.listdir(ptm_dir)
                if len(ptm_files) == 0:
                    log.error(f"no files found in {ptm_files}")

                ptm_data = []
                for fname in ptm_files:
                    try:
                        ptm_data.append(self._parse_ptm(os.path.join(analysis_basepath, 'ptm', fname)))
                    except (IOError, TypeError):
                        exit(1)

                ptm_data = pd.concat(ptm_data)
                ptm_data = ptm_data.rename(columns={ '#ptm'                 : "PTMs",
                                                     'SASA(%)'              : "PTM residue SASA (%)" ,
                                                     'ddG(foldX5,kcal/mol)' : "Change in stability with PTM (FoldX5, kcal/mol)",
                                                     'effect_regulation'    : "PTM effect in regulation",
                                                     'effect_stability'     : "PTM effect in stability" ,
                                                     'effect_function'      : "PTM effect in function",
                                                     'notes'                : "PTM notes"})

                this_df = this_df.join(ptm_data)

            this_df = this_df.reset_index()
            this_df = this_df.fillna(pd.NA)
            data_dfs[(system, mode)] = this_df

        return data_dfs

