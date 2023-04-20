# MAVISp - classes for handling different methods
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


import pandas as pd
import numpy as np
from mavisp.error import *
import os

class Method(object):
    name = "name"

    def __init__(self, version):
        self.version = version

        self.data = None

class MutateXStability(Method):

    unit = "kcal/mol"
    type = "Stability"

    def parse(self, dir_path):

        warnings = []

        mutatex_files = os.listdir(dir_path)

        if len(mutatex_files) != 1:
            this_error = f"zero or multiple files found in {dir_path}; only one expected"
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[MAVISpCriticalError(this_error)])

        mutatex_file = mutatex_files[0]

        try:
            df = pd.read_csv(os.path.join(dir_path, mutatex_file))
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the MutateX csv file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[MAVISpCriticalError(this_error)])

        # create residue column
        df['residue'] = df['WT residue type'] + df['Residue #'].astype(str)

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
        df = df.rename(columns={0 : f"{self.type} ({self.version}, {self.unit})"})

        self.data = df

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[])

class MutateXBinding(Method):

    unit = "kcal/mol"
    type = "Local Int."
    chain = "A"
    measure = "Binding with"

    def parse(self, dir_path):

        warnings = []

        interactors = os.listdir(dir_path)

        if len(interactors) == 0:
            raise MAVISpMultipleError(critical=[MAVISpCriticalError("no interactor folders found")],
                                      warning=[])

        all_data = None

        for interactor in interactors:

            interactor_dir = os.path.join(dir_path, interactor)

            mutatex_files = os.listdir(interactor_dir)

            if len(mutatex_files) != 1:
                raise MAVISpMultipleError(critical=[MAVISpCriticalError("zero or multiple files found in {interactor_dir}; exactly one expected")],
                                        warning=[])

            mutatex_file = mutatex_files[0]
        
            try:
                df = pd.read_csv(os.path.join(interactor_dir, mutatex_file))
            except Exception as e:
                this_error = f"Exception {type(e).__name__} occurred when parsing the MutateX csv file. Arguments:{e.args}"
                raise MAVISpMultipleError(warning=warnings, 
                                        critical=[MAVISpCriticalError(this_error)])

            # create residue column
            df['residue'] = df['WT residue type'] + df['Residue #'].astype(str)

            if self.chain is not None:
                df = df[ df['chain ID'] == self.chain ]

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

            # handle space around measure
            if self.measure == "":
                measure = ""
            else:
                measure = "{self.measure} "

            df = df.rename(columns={0 : f"{self.type} ({self.measure}{interactor}, {self.version}, {self.unit})"})

            if all_data is None:
                all_data = df
            else:
                all_data = all_data.join(df)

        self.data = all_data

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[])

class MutateXDNABinding(Method):

    unit = "kcal/mol"
    type = "Local Int. With DNA"
    chain = "A"
    measure = ""

    parse = MutateXBinding.parse

class RosettaDDGPredictionStability(Method):

    unit = "kcal/mol"
    type = "Stability"

    def parse(self, dir_path):

        warnings = []

        rosetta_files = os.listdir(dir_path)

        if len(rosetta_files) != 1:
            this_error = f"multiple files found in {dir_path}; only one expected"
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[MAVISpCriticalError(this_error)])

        rosetta_file = rosetta_files[0]

        try:
            mutation_data = pd.read_csv(os.path.join(dir_path, rosetta_file))
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the Rosetta csv file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[MAVISpCriticalError(this_error)])

        mutation_data = mutation_data[mutation_data['state'] == 'ddg']
        mutation_data = mutation_data.set_index('mutation_label')
        mutation_data = mutation_data[['total_score']]
        mutation_data = mutation_data.rename(columns={'total_score':f'{self.type} ({self.version}, {self.unit})'})

        self.data = mutation_data

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[])

class RosettaDDGPredictionBinding(Method):

    unit = "kcal/mol"
    type = "Local Int."
    chain = 'A'

    def parse(self, dir_path):

        warnings = []

        interactors = os.listdir(dir_path)

        if len(interactors) == 0:
            raise MAVISpMultipleError(critical=[MAVISpCriticalError("no interactor folders found")],
                                        warning=[])

        all_data = None

        for interactor in interactors:

            interactor_dir = os.path.join(dir_path, interactor)

            rosetta_files = os.listdir(interactor_dir)

            if len(rosetta_files) != 1:
                raise MAVISpMultipleError(critical=[MAVISpCriticalError("zero or multiple files found in {interactor_dir}; exactly one expected")],
                                        warning=[])

            rosetta_file = rosetta_files[0]
        
            try:
                mutation_data = pd.read_csv(os.path.join(interactor_dir, rosetta_file))
            except Exception as e:
                this_error = f"Exception {type(e).__name__} occurred when parsing the Rosetta csv file. Arguments:{e.args}"
                raise MAVISpMultipleError(warning=warnings, 
                                        critical=[MAVISpCriticalError(this_error)])

            mutation_data = mutation_data[mutation_data['state'] == 'ddg']
            if self.chain is not None:
                mutation_data = mutation_data[ mutation_data['mutation'].str[0] == self.chain ]
            mutation_data = mutation_data.set_index('mutation_label')
            mutation_data = mutation_data[['total_score']]
            mutation_data = mutation_data.rename(columns={'total_score':f'{self.type} (Binding with {interactor}, {self.version}, {self.unit})'})

            if all_data is None:
                all_data = mutation_data
            else:
                all_data = all_data.join(mutation_data)

        self.data = all_data

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[])

class AlloSigma(Method):
    name = "AlloSigma"

    allowed_fnames = set(['allosigma_mut.txt', 'filtered_down_mutations.tsv', 'filtered_up_mutations.tsv'])

    def _process_allosigma2_tables(self, row, filt_up, filt_down, cutoff):
        if pd.isna(row['allosigma-mode']):
            return '-'

        if row['allosigma-mode'] == 'UP':
            filt = filt_up
        elif row['allosigma-mode'] == 'DOWN':
            filt = filt_down

        if filt is None:
            return 'neutral'

        if row['mutations'] not in filt.index:
            return 'neutral'

        entry = filt.loc[row['mutations']]
        entry = entry.drop(['avg_dG', 'n_mutations'])
        values = entry[~ pd.isna(entry)].values

        if np.all(values > cutoff):
            return 'destabilizing'
        elif np.all(values < (-cutoff)):
            return 'stabilizing'
        else:
            return 'mixed_effects'


    def parse(self, allosigma2_dir):

        warnings = []

        # check what files are available
        allosigma2_files = set(os.listdir(allosigma2_dir))
        if   np.all([ os.path.isfile(f"{allosigma2_dir}/{f}") for f in allosigma2_files ]):
            allosigma2_domains = {'.' : allosigma2_files}
        elif np.all([ os.path.isdir(f"{allosigma2_dir}/{d}")  for d in allosigma2_files ]):
            allosigma2_domains = {}
            for d in allosigma2_files:
                allosigma2_domains[d] = set(os.listdir(f"{allosigma2_dir}/{d}"))
        else:
            raise MAVISpMultipleError(critical=[MAVISpCriticalError(f"the AlloSigma2 directory should contain either only files or only directories")],
                        warning=[])

        allosigma2_data = []

        for dirname, allosigma2_files in allosigma2_domains.items():

            if len(allosigma2_files) not in [1, 2, 3]:
                raise MAVISpMultipleError(critical=[MAVISpCriticalError(f"the AlloSigma2 directory {allosigma2_dir}/{dirname} should contain only 2 or 3 files")],
                                        warning=[])
                                        
            if 'allosigma_mut.txt' not in allosigma2_files:
                raise MAVISpMultipleError(critical=[MAVISpCriticalError(f"the allosigma_mut.txt file must be present in the AlloSigma2 directory {allosigma2_dir}/{dirname}")],
                                        warning=[])

            if not allosigma2_files.issubset(self.allowed_fnames):
                raise MAVISpMultipleError(critical=[MAVISpCriticalError(f"the only allowed file names in the allosigma2 directory {allosigma2_dir}/{dirname} are {', '.join(list(self.allowed_fnames))}")],
                                        warning=[])

            try:
                all_mut = pd.read_csv(os.path.join(allosigma2_dir, dirname, 'allosigma_mut.txt'), sep='\t')
            except Exception as e:
                this_error = f"Exception {type(e).__name__} occurred when parsing the csv files. Arguments:{e.args}"
                raise MAVISpMultipleError(warning=warnings, 
                                        critical=[MAVISpCriticalError(this_error)])

            all_mut = all_mut.drop_duplicates()

            all_mut['mutations'] = all_mut.wt_residue + all_mut.position.astype(str) + all_mut.mutated_residue

            if 'filtered_down_mutations.tsv' in allosigma2_files:
                try:
                    filt_down = pd.read_csv(os.path.join(allosigma2_dir, dirname, 'filtered_down_mutations.tsv'), sep='\t', index_col=0)
                except Exception as e:
                    this_error = f"Exception {type(e).__name__} occurred when parsing filtered_down_mutations.tsv. Arguments:{e.args}"
                    raise MAVISpMultipleError(warning=warnings, 
                                            critical=[MAVISpCriticalError(this_error)])

                filt_down['mutations'] = filt_down['mutations'].str.split()
                filt_down = filt_down.explode('mutations')
                filt_down = filt_down.set_index('mutations')
            else:
                filt_down = None

            if 'filtered_up_mutations.tsv' in allosigma2_files:
                try:
                    filt_up   = pd.read_csv(os.path.join(allosigma2_dir, dirname, 'filtered_up_mutations.tsv'), sep='\t', index_col=0)
                except Exception as e:
                    this_error = f"Exception {type(e).__name__} occurred when parsing filtered_up_mutations.tsv. Arguments:{e.args}"
                    raise MAVISpMultipleError(warning=warnings, 
                                            critical=[MAVISpCriticalError(this_error)])

                filt_up['mutations'] = filt_up['mutations'].str.split()
                filt_up = filt_up.explode('mutations')
                filt_up = filt_up.set_index('mutations')
            else:
                filt_up = None

            all_mut['allosigma-consequence'] = all_mut.apply(self._process_allosigma2_tables, filt_up=filt_up, filt_down=filt_down, cutoff=1, axis=1)

            all_mut = all_mut[['mutations', 'allosigma-mode', 'allosigma-consequence']].set_index('mutations')

            all_mut = all_mut.rename(columns={'allosigma-mode': f'AlloSigma{self.version} mutation type',
                                            'allosigma-consequence' : f'AlloSigma{self.version} predicted consequence'}).fillna('-')

            allosigma2_data.append(all_mut)
        
        self.data = pd.concat(allosigma2_data)

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[])

class RaSP(Method):

    unit = "kcal/mol"
    type = "Stability"

    def parse(self, dir_path):

        warnings = []

        rasp_files = os.listdir(dir_path)

        if len(rasp_files) != 1:
            this_error = f"zero or multiple files found in {dir_path}; only one expected"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        rasp_file = rasp_files[0]

        try:
            mutation_data = pd.read_csv(os.path.join(dir_path, rasp_file),
                                        usecols=['variant', 'RaSP_ddG'], index_col='variant')

        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the RaSP csv file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        mutation_data = mutation_data.rename(columns={'RaSP_ddG' : f"{self.type} ({self.version}, {self.unit})"})

        self.data = mutation_data

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])
