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

        # drop empty column if it exists
        cols_to_drop = [ c for c in df.columns if c.startswith('Unnamed') ]
        df = df.drop(columns=cols_to_drop)

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

        return df, warnings

class MutateXBinding(Method):

    unit = "kcal/mol"
    type = "Local Int."
    heterodimer_chains = set(['A'])
    homodimer_chains   = set(['AB'])
    target_chain       = 'A'
    measure = "Binding with"
    complex_status = "heterodimer"

    def __init__(self, version, complex_status=None):

        super().__init__(version)

        if complex_status is not None:
            self.complex_status = complex_status

        self.interactors = []

    def parse(self, dir_path):

        warnings = []

        interactors = os.listdir(dir_path)
        self.interactors = interactors

        if len(interactors) == 0:
            raise MAVISpMultipleError(critical=[MAVISpCriticalError("no interactor folders found")],
                                      warning=warnings)

        all_data = None

        for interactor in interactors:

            interactor_dir = os.path.join(dir_path, interactor)

            mutatex_files = os.listdir(interactor_dir)

            if len(mutatex_files) != 1:
                raise MAVISpMultipleError(critical=[MAVISpCriticalError(f"zero or multiple files found in {interactor_dir}; exactly one expected")],
                                        warning=warnings)

            mutatex_file = mutatex_files[0]

            try:
                df = pd.read_csv(os.path.join(interactor_dir, mutatex_file))
            except Exception as e:
                this_error = f"Exception {type(e).__name__} occurred when parsing the MutateX csv file. Arguments:{e.args}"
                raise MAVISpMultipleError(warning=warnings,
                                        critical=[MAVISpCriticalError(this_error)])

            # create residue column
            df['residue'] = df['WT residue type'] + df['Residue #'].astype(str)

            # detect and handle homodimer case
            chains = set(df['chain ID'].unique())

            if self.target_chain in chains:
                df = df[ df['chain ID'] == self.target_chain ]

            elif set(df['chain ID'].unique()) != self.homodimer_chains:
                message = "chain ID in FoldX energy file must be either A or B (heterodimer case) or AB (homodimer case)"
                raise MAVISpMultipleError(critical=[MAVISpCriticalError(message)],
                                          warning=[])


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
                measure = f"{self.measure} "

            df = df.rename(columns={0 : f"{self.type} ({measure}{interactor}, {self.complex_status}, {self.version}, {self.unit})"})

            if all_data is None:
                all_data = df
            else:
                all_data = all_data.join(df, how='outer')

        return all_data, warnings

class MutateXDNABinding(Method):

    unit = "kcal/mol"
    type = "Local Int. With DNA"
    heterodimer_chains = set(['A'])
    homodimer_chains   = set(['AB'])
    target_chain       = 'A'
    measure = ""
    complex_status = "heterodimer"

    parse = MutateXBinding.parse

class RosettaDDGPredictionStability(Method):

    unit = "kcal/mol"
    type = "Stability"

    def _parse_aggregate_csv(self, csvf, warnings):
        try:
            df = pd.read_csv(csvf)
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the Rosetta csv file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        df = df[df['state'] == 'ddg']
        df = df[['total_score', 'mutation_label']]

        # homodimer mode and heterodimer mode can be treated similarly
        df['mutation_label'] = df['mutation_label'].str.split('_')

        # check if all mutation labels have the same number of residues
        if len(set(df['mutation_label'].apply(len))) != 1:
            this_error = f"RosettaDDG aggregate file contains values for both homodimer and heterodimer modes"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        # check if all homodimer modes refer to the same residue (e.g. there are no pairs with different mutations)
        if not np.all(df['mutation_label'].apply(lambda r: len(set(r))) == 1):
            this_error = f"RosettaDDG aggregate file contains entries that refer to residues with different numbering"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        df['mutation_label'] = df['mutation_label'].str[0]
        df = df.set_index('mutation_label')
        df = df[['total_score']]

        return df

    def parse(self, dir_path):

        warnings = []

        rosetta_files = os.listdir(dir_path)

        if len(rosetta_files) == 1 and os.path.isfile(os.path.join(dir_path, rosetta_files[0])):
            rosetta_file = rosetta_files[0]

            mutation_data = self._parse_aggregate_csv(os.path.join(dir_path, rosetta_file), warnings)

        else:
            csv_files = []
            rosetta_folder = os.listdir(dir_path)

            # Check if all available files are directories
            for folder in rosetta_folder:
                if not os.path.isdir(os.path.join(dir_path, folder)):
                    this_error = f"{folder} in {dir_path} is not a directory"
                    raise MAVISpMultipleError(warning=warnings,
                                                critical=[MAVISpCriticalError(this_error)])

                ddg_files = os.listdir(os.path.join(dir_path, folder))
                # check one file per directory is available
                if len(ddg_files) != 1:
                    this_error = f"zero or multiples files found in {dir_path}; only one expected"
                    raise MAVISpMultipleError(warning=warnings,
                                              critical=[MAVISpCriticalError(this_error)])

                ddg_file = ddg_files[0]

                if not os.path.splitext(ddg_file)[-1] == '.csv':
                    this_error = f"file {ddg_file} in {dir_path}/{folder} doesn't have csv as extension"
                    raise MAVISpMultipleError(warning=warnings,
                                              critical=[MAVISpCriticalError(this_error)])

                csv_files.append(os.path.join(dir_path, folder, ddg_file))

            list_mutation_label = None
            mutation_data = None

            for fname in csv_files:
                tmp = self._parse_aggregate_csv(fname, warnings)

                # Check if the mutation labels are the same in the different csv files
                if list_mutation_label is None:
                    list_mutation_label = set(tmp.index)
                elif list_mutation_label != set(tmp.index):
                    this_error = f"the mutation labels are not the same in the different csv files"
                    raise MAVISpMultipleError(warning=warnings,
                                            critical=[MAVISpCriticalError(this_error)])

                # Allow to merge the data from the different cl folders
                if mutation_data is None:
                    mutation_data = tmp
                else:
                    mutation_data = mutation_data.join(tmp, rsuffix="_")

            # merge the data from the different cl folders and keep average
            ddg_colname = f'{self.type} ({self.version}, {self.unit})'
            mutation_data[ddg_colname] = mutation_data.mean(axis=1)
            mutation_data = mutation_data[[ddg_colname]]

            # Sort the data by mutation_label
            mutation_data = mutation_data.sort_index()

        return mutation_data, warnings

class RosettaDDGPredictionBinding(Method):

    unit = "kcal/mol"
    type = "Local Int."
    chain = 'A'
    complex_status = 'heterodimer'

    _parse_aggregate_csv = RosettaDDGPredictionStability._parse_aggregate_csv

    def __init__(self, version, complex_status=None):

        super().__init__(version)

        if complex_status is not None:
            self.complex_status = complex_status

        self.interactors = []

    def parse(self, dir_path):

        warnings = []

        interactors = os.listdir(dir_path)
        self.interactors = interactors

        if len(interactors) == 0:
            raise MAVISpMultipleError(critical=[MAVISpCriticalError("no interactor folders found")],
                                        warning=[])

        all_data = None

        for interactor in interactors:

            interactor_dir = os.path.join(dir_path, interactor)
            rosetta_files = os.listdir(interactor_dir)

            if len(rosetta_files) == 1 and os.path.isfile(os.path.join(interactor_dir, rosetta_files[0])):

                rosetta_file = rosetta_files[0]
                mutation_data = self._parse_aggregate_csv(os.path.join(interactor_dir, rosetta_file), warnings)

            elif len(rosetta_files) > 1 and all( [ os.path.isdir(os.path.join(interactor_dir, f)) for f in rosetta_files] ):
                mutation_data = None
                for c, conformer_dir in enumerate(rosetta_files):

                    conformer_files = os.listdir(os.path.join(interactor_dir, conformer_dir))

                    if len(conformer_files) != 1:
                        text = "only one file per conformer is supported for RosettaDDGPrediction"
                        raise MAVISpMultipleError(critical=[MAVISpCriticalError(text)],
                                                  warning=warnings)

                    conformer_data = self._parse_aggregate_csv(os.path.join(interactor_dir, conformer_dir, conformer_files[0]), warnings)

                    conformer_data = conformer_data.rename(columns={'total_score' : f'total_score_{c}'})

                    if mutation_data is None:
                        mutation_data = conformer_data
                    else:
                        mutation_data = mutation_data.join(conformer_data)

                mutation_data = pd.DataFrame(mutation_data.mean(axis=1), columns=['total_score'])

            else:
                text = f"dataset {interactor_dir} was not either a single files, or multiple directories containing one file"
                raise MAVISpMultipleError(critical=[MAVISpCriticalError(text)],
                                            warning=warnings)

            mutation_data = mutation_data.rename(columns={'total_score':f'{self.type} (Binding with {interactor}, {self.complex_status}, {self.version}, {self.unit})'})

            if all_data is None:
                all_data = mutation_data
            else:
                all_data = all_data.join(mutation_data, how='outer')

        return all_data, warnings

class AlloSigma(Method):

    name = "AlloSigma"

    datasets = {'sub_cat'  : 'AlloSigma{version} predicted consequence - active sites',
                'cofactor' : 'AlloSigma{version} predicted consequence - cofactor sites',
                'pockets'  : 'AlloSigma{version} predicted consequence - pockets and interfaces'}

    allowed_fnames = [f'filtered_up_{x}.tsv' for x in datasets.keys()] +\
                     [f'filtered_down_{x}.tsv' for x in datasets.keys()] +\
                     ['allosigma_mut.txt']

    def _process_allosigma2_tables(self, row, filt_up, filt_down, cutoff):

        if pd.isna(row['allosigma-mode']):
            return 'uncertain'

        if row['allosigma-mode'] == 'UP':
            filt = filt_up
        elif row['allosigma-mode'] == 'DOWN':
            filt = filt_down

        if filt is None:
            return 'neutral'

        if row['mutations'] not in filt.index:
            return 'neutral'

        entry = filt.loc[row['mutations']]
        values = entry[~ pd.isna(entry)].values

        if np.all(values > cutoff):
            return 'destabilizing'
        elif np.all(values < (-cutoff)):
            return 'stabilizing'
        else:
            return 'mixed_effects'

    def _parse_allosigma2_energy_table(self, fname):

        filt = pd.read_csv(fname, sep='\t', index_col=0)

        filt = filt.reset_index()
        filt['mutations'] = filt['mutations'].str.split()
        filt = filt.explode('mutations')
        filt = filt.set_index('mutations')

        return filt

    def _parse_allosigma2_mutation_file(self, fname):

        all_mut = pd.read_csv(fname, sep='\t')

        all_mut['mutations'] = all_mut.wt_residue + all_mut.position.astype(str) + all_mut.mutated_residue

        return all_mut

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

            # check if the files are the expected ones
            if 'allosigma_mut.txt' not in allosigma2_files:
                raise MAVISpMultipleError(critical=[MAVISpCriticalError(f"the allosigma_mut.txt file must be present in the AlloSigma2 directory {allosigma2_dir}/{dirname}")],
                                        warning=[])

            if not allosigma2_files.issubset(self.allowed_fnames):
                raise MAVISpMultipleError(critical=[MAVISpCriticalError(f"the only allowed file names in the allosigma2 directory {allosigma2_dir}/{dirname} are {', '.join(list(self.allowed_fnames))}")],
                                        warning=[])

            # parse the mutation file
            try:
                all_mut = self._parse_allosigma2_mutation_file(os.path.join(allosigma2_dir, dirname, 'allosigma_mut.txt'))
            except Exception as e:
                this_error = f"Exception {type(e).__name__} occurred when parsing the csv files. Arguments:{e.args}"
                raise MAVISpMultipleError(warning=warnings,
                                          critical=[MAVISpCriticalError(this_error)])

            # for every data type
            for suffix, colname in self.datasets.items():

                data = {}

                # for either up or down mutations for the type
                for direction in ['up', 'down']:

                    fname = f'filtered_{direction}_{suffix}.tsv'

                    if fname in allosigma2_files:
                        try:
                            data[direction] = self._parse_allosigma2_energy_table(os.path.join(allosigma2_dir, dirname, fname))
                        except Exception as e:
                            this_error = f"Exception {type(e).__name__} occurred when parsing {fname}. Arguments:{e.args}"
                            raise MAVISpMultipleError(warning=warnings,
                                                    critical=[MAVISpCriticalError(this_error)])
                    else:
                        data[direction] = None

                all_mut[suffix] = all_mut.apply(self._process_allosigma2_tables,
                                                                          filt_up=data['up'],
                                                                          filt_down=data['down'],
                                                                          cutoff=2,
                                                                          axis=1)

            all_mut = all_mut[['mutations', 'allosigma-mode'] + list(self.datasets.keys())]
            all_mut = all_mut.set_index('mutations')

            allosigma2_data.append(all_mut)

        out_data = pd.concat(allosigma2_data, axis=0)

        rename_scheme = {'allosigma-mode' : f'AlloSigma{self.version} mutation type'}
        for k, v in self.datasets.items():
            rename_scheme[k] = v.format(version=self.version)

        return out_data.rename(columns=rename_scheme).fillna('-'), warnings

class RaSP(Method):

    unit = "kcal/mol"
    type = "Stability"

    def _parse_postprocessed_csv(self, fname, warnings):

        warnings = []

        try:
            mutation_data = pd.read_csv(fname,
                                        usecols=['variant', 'RaSP_ddG'], index_col='variant')

        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the RaSP csv file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        return mutation_data

    def parse(self, dir_path):

        warnings = []

        rasp_files = os.listdir(dir_path)

        if len(rasp_files) == 1 and os.path.isfile(os.path.join(dir_path, rasp_files[0])):

            rasp_file = rasp_files[0]

            mutation_data = self._parse_postprocessed_csv(os.path.join(dir_path, rasp_file), warnings)

        else:
            csv_files = []
            rasp_folder = os.listdir(dir_path)

            # Check if all available files are directories
            for folder in rasp_folder:
                if not os.path.isdir(os.path.join(dir_path, folder)):
                    this_error = f"{folder} in {dir_path} is not a directory"
                    raise MAVISpMultipleError(warning=warnings,
                                                critical=[MAVISpCriticalError(this_error)])

                ddg_files = os.listdir(os.path.join(dir_path, folder))
                # check one file per directory is available
                if len(ddg_files) != 1:
                    this_error = f"zero or multiples files found in {dir_path}; only one expected"
                    raise MAVISpMultipleError(warning=warnings,
                                              critical=[MAVISpCriticalError(this_error)])

                ddg_file = ddg_files[0]

                if not os.path.splitext(ddg_file)[-1] == '.csv':
                    this_error = f"file {ddg_file} in {dir_path}/{folder} doesn't have csv as extension"
                    raise MAVISpMultipleError(warning=warnings,
                                              critical=[MAVISpCriticalError(this_error)])

                csv_files.append(os.path.join(dir_path, folder, ddg_file))

            list_mutation_label = None
            mutation_data = None

            for fname in csv_files:
                tmp = self._parse_postprocessed_csv(fname, warnings)

                # Check if the mutation labels are the same in the different csv files
                if list_mutation_label is None:
                    list_mutation_label = set(tmp.index)
                elif list_mutation_label != set(tmp.index):
                    this_error = f"the mutation labels are not the same in the different csv files"
                    raise MAVISpMultipleError(warning=warnings,
                                            critical=[MAVISpCriticalError(this_error)])

                # Allow to merge the data from the different cl folders
                if mutation_data is None:
                    mutation_data = tmp
                else:
                    mutation_data = mutation_data.join(tmp, rsuffix="_")

            # merge the data from the different cl folders and keep average
        ddg_colname = f'{self.type} ({self.unit})'
        mutation_data[ddg_colname] = mutation_data.mean(axis=1)
        mutation_data = mutation_data[[ddg_colname]]

        mutation_data = mutation_data.sort_index()

        return mutation_data, warnings
