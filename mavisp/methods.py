# MAVISp - classes for handling different methods
# Copyright (C) 2022 Matteo Tiberti, Danish Cancer Society
#           (C) 2026 Eszter Toldi, Technical University of Denmark (DTU)
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
import re

class Method(object):
    name = "name"

    def __init__(self, version):
        self.version = version

class MutateXStability(Method):

    unit = "kcal/mol"
    type = "Stability"
    averages_filename = 'energies.csv'
    stds_filename = 'energies_std.csv'
    conf_filename = 'damaging_proportion.csv'

    def _parse_mutatex_energy_file(self, fname, data_type):
        try:
            df = pd.read_csv(fname)
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

        if data_type is None or data_type == '':
            colname =  f"{self.type} ({self.version}, {self.unit})"
        else:
            colname =  f"{self.type} ({self.version}, {self.unit}, {data_type})"

        return df.rename(columns={0 : colname})

    def _parse_conformations_file(self, fname):
        try:
            df = pd.read_csv(fname, usecols=['file', 'proportion_damaging'])
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the MutateX conformation-dependentfile. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        pattern = r'^([A-Za-z])(\d+)_heatmap_ranking_([A-Za-z])\.csv$'

        df['mutations'] = df['file'].str.extract(pattern).apply(lambda x: f"{x[0]}{x[1]}{x[2]}", axis=1)
        df = df.set_index('mutations')
        df = df[['proportion_damaging']]
        df['proportion_damaging'] = df['proportion_damaging'].round(2)
        df.columns = [f"{self.version} proportion of damaging conformations"]

        return df

    def parse(self, dir_path):

        warnings = []

        mutatex_files = os.listdir(dir_path)

        if self.averages_filename not in mutatex_files:
            this_error = f"energies.csv file not found in {dir_path}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        averages_df = self._parse_mutatex_energy_file(os.path.join(dir_path, self.averages_filename), '')

        if self.stds_filename in mutatex_files:
            stds_df = self._parse_mutatex_energy_file(os.path.join(dir_path, self.stds_filename), 'st. dev.')
        else:
            warnings.append(MAVISpWarningError("standard deviation file not found for MutateX data"))
            stds_df = None

        if self.conf_filename in mutatex_files:
            conf_mutation_data = self._parse_conformations_file(os.path.join(dir_path, self.conf_filename))
        else:
            warnings.append(MAVISpWarningError("Conformation-dependent stability population file not found for MutateX data"))
            conf_mutation_data = None

        return averages_df, stds_df, conf_mutation_data, warnings

class MutateXBinding(Method):
    unit = "kcal/mol"
    type = "Local Int."
    heterodimer_chains = set(['A'])
    homodimer_chains   = set(['AB'])
    target_chain       = 'A'
    measure = "Binding with"
    averages_filename = 'energies.csv'
    stds_filename = 'energies_std.csv'
    complex_status = "heterodimer"

    def __init__(self, version, complex_status=None):

        super().__init__(version)

        if complex_status is not None:
            self.complex_status = complex_status

        self.interactors = []

    # data_type is either '' or 'st. dev.'
    def _parse_mutatex_binding_energy_file(self, fname, interactor, data_type):
        """Parse a single MutateX binding energy file (average or std)."""

        try:
            df = pd.read_csv(fname)
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the MutateX csv file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        # Create residue column
        df['residue'] = df['WT residue type'] + df['Residue #'].astype(str)

        # Detect and handle homodimer case
        chains = set(df['chain ID'].unique())

        # If chain A exists, KEEP ONLY rows where chain ID == 'A'.
        if self.target_chain in chains:
            df = df[ df['chain ID'] == self.target_chain ]
        # If chain A does NOT exist, the ONLY valid alternative is homodimer AB.
        elif set(df['chain ID'].unique()) != self.homodimer_chains:
            message = "chain ID in FoldX energy file must be either A or B (heterodimer case) or AB (homodimer case)"
            raise MAVISpMultipleError(critical=[MAVISpCriticalError(message)],
                                      warning=[])

        # Drop unnecessary columns
        df = df.drop(['WT residue type', 'Residue #', 'chain ID'], axis=1)

        # Stack remaining columns
        df = df.set_index('residue') # set 'residue' column as index
        df = df.stack() # rotates columns downward and makes the dataframe long-format (level_1 contains the original column names and 0 contains the values)
        df = df.reset_index() # reset index to turn the index into a column

        # Create mutation column
        df['mutations'] = df['residue'] + df['level_1'] # concatenate 'residue' and 'level_1' columns to create 'mutations' column
        df = df.set_index('mutations') # set 'mutations' column as index

        # Drop now useless columns, rename
        df = df.drop(['residue', 'level_1'], axis=1)

        # handle space around measure
        if self.measure == "":
            measure = ""
        else:
            measure = f"{self.measure} "

        # rename column Local Int. (Binding with B, heterodimer, FoldX5, kcal/mol)
        if data_type is None or data_type == '':
            data_type = ""
        else:
            data_type = f", {data_type}"

        colname = f"{self.type} ({measure}{interactor}, {self.complex_status}, {self.version}, {self.unit}{data_type})"

        return df.rename(columns={0 : colname}) # rename the sinle column named 0 to the formatted name

    def parse(self, dir_path):
        """
        reads the MutateX output files (energies.csv + energies_std.csv) for each interactor,
        converts them into mutation-indexed dataframes, and returns them to the Local interaction module.
        """

        warnings = []
        all_data = None

        interactors = os.listdir(dir_path) # list of subfolders in dir_path
        self.interactors = interactors # store interactors in the instance variable

        if len(interactors) == 0:
            raise MAVISpMultipleError(critical=[MAVISpCriticalError("no interactor folders found")],
                                      warning=warnings)

        for interactor in interactors:

            interactor_dir = os.path.join(dir_path, interactor)
            mutatex_files = os.listdir(interactor_dir)

            # expect energies.csv file per interactor
            if self.averages_filename not in mutatex_files:
                this_error = f"energies.csv file not found in {interactor_dir}"
                raise MAVISpMultipleError(warning=warnings,
                                          critical=[MAVISpCriticalError(this_error)])

            # Parse averages file
            averages_df = self._parse_mutatex_binding_energy_file(os.path.join(interactor_dir, self.averages_filename), interactor, '')

            # Parse stds file if it exists
            if self.stds_filename in mutatex_files:
                stds_df = self._parse_mutatex_binding_energy_file(os.path.join(interactor_dir, self.stds_filename), interactor, 'st. dev.')
            else:
                warnings.append(MAVISpWarningError("standard deviation file not found for MutateX binding energy data"))
                stds_df = None

            # Combine averages and stds data
            if stds_df is not None:
                interactor_data = averages_df.join(stds_df, how='outer')
            else:
                interactor_data = averages_df

            # Combine data across interactors
            if all_data is None:
                all_data = interactor_data
            else:
                all_data = all_data.join(interactor_data, how='outer')

        return all_data, warnings

class MutateXDNABinding(MutateXBinding):

    unit = "kcal/mol"
    type = "Local Int. With DNA"
    heterodimer_chains = set(['A'])
    homodimer_chains   = set(['AB'])
    target_chain       = 'A'
    measure = ""
    complex_status = "heterodimer"

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

            avg_mutation_data = self._parse_aggregate_csv(os.path.join(dir_path, rosetta_file), warnings)
            std_mutation_data = None

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

            avg_ddg_colname = f'{self.type} ({self.version}, {self.unit})'
            std_ddg_colname = f'{self.type} ({self.version}, {self.unit}, st. dev.)'

            mutation_data_mean = mutation_data.mean(axis=1)
            mutation_data_std = mutation_data.std(axis=1)

            mutation_data[avg_ddg_colname] = mutation_data_mean
            mutation_data[std_ddg_colname] = mutation_data_std

            mutation_data = mutation_data[[avg_ddg_colname, std_ddg_colname]]

            # Sort the data by mutation_label
            mutation_data = mutation_data.sort_index()

            avg_mutation_data = mutation_data[[ c for c in mutation_data.columns if not 'st. dev.' in c ]]
            std_mutation_data = mutation_data[[ c for c in mutation_data.columns if ', st. dev.)' in c ]]

        return avg_mutation_data, std_mutation_data, None, warnings

class RosettaDDGPredictionBinding(RosettaDDGPredictionStability):

    unit = "kcal/mol"
    type = "Local Int."
    chain = 'A'
    complex_status = 'heterodimer'
    aggregate_fname = 'ddg_mutations_aggregate.csv'
    structures_fname = 'ddg_mutations_structures.csv'

    def __init__(self, version, complex_status=None):

        super().__init__(version)

        if complex_status is not None:
            self.complex_status = complex_status

        self.interactors = []

    def _parse_structure_csv(self, csvf, warnings):
        """Parse the RosettaDDGPrediction binding structure CSV file."""
        try:
            df = pd.read_csv(csvf)
        except Exception as e:
            this_error = f"Exception {type(e).__name__} while reading structure CSV: {e.args}"
            raise MAVISpMultipleError(
                warning=warnings,
                critical=[MAVISpCriticalError(this_error)]
            )

        #keep only ddg rows
        df = df[df["state"] == "ddg"]

        #group by mutation and compute stdev of total_score
        std_series = df.groupby("mutation_label")["total_score"].std()
        average_series = df.groupby("mutation_label")["total_score"].mean()

        #turn into DataFrame
        std_df = std_series.to_frame(name ="total_score")
        average_df = average_series.to_frame(name ="total_score")

        return average_df, std_df

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

            # Identify the correct files
            agg_file = None
            struct_file = None

            # Expect either a single file or multiple directories containing one file each
            if len(rosetta_files) == 1 and os.path.isfile(os.path.join(interactor_dir, rosetta_files[0])) and rosetta_files[0] == self.aggregate_fname:
                rosetta_file = rosetta_files[0]
                # Parse single aggregate CSV file
                mutation_data = self._parse_aggregate_csv(os.path.join(interactor_dir, rosetta_file), warnings)
                mutation_data = mutation_data.rename(columns={'total_score':f'{self.type} (Binding with {interactor}, {self.complex_status}, {self.version}, {self.unit})'})

            # or structures file (which we can use to calculate average and stdev)
            elif (len(rosetta_files) == 1 and\
                    os.path.isfile(os.path.join(interactor_dir, rosetta_files[0])) and \
                    rosetta_files[0] == self.structures_fname)    or\
                 (set(rosetta_files) == set([self.aggregate_fname, self.structures_fname]) and\
                    os.path.isfile(os.path.join(interactor_dir, rosetta_files[0])) and\
                    os.path.isfile(os.path.join(interactor_dir, rosetta_files[1]))):

                if len(rosetta_files) == 2:
                    warnings.append(MAVISpWarningError(f"for {interactor}, both Rosetta aggregate and structures file were found; the aggregate file will be ignored"))

                mutation_data, std_df = self._parse_structure_csv(os.path.join(interactor_dir, self.structures_fname), warnings)

                mutation_data = mutation_data.rename(columns={'total_score':f'{self.type} (Binding with {interactor}, {self.complex_status}, {self.version}, {self.unit})'})
                std_df = std_df.rename(columns={'total_score' : f"{self.type} (Binding with {interactor}, {self.complex_status}, {self.version}, {self.unit}, st. dev.)"})

                mutation_data = mutation_data.join(std_df, how="outer")

            # Multiple directories containing one file each
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

                # Average total_score across conformers

                mutation_data = pd.DataFrame({f'{self.type} (Binding with {interactor}, {self.complex_status}, {self.version}, {self.unit})' : mutation_data.mean(axis=1),
                                              f'{self.type} (Binding with {interactor}, {self.complex_status}, {self.version}, {self.unit}, st. dev.)' : mutation_data.std(axis=1)})

            else:
                text = f"dataset {interactor_dir} did not contain an expected folder structure"
                raise MAVISpMultipleError(critical=[MAVISpCriticalError(text)],
                                            warning=warnings)

            if all_data is None:
                all_data = mutation_data
            else:
                all_data = all_data.join(mutation_data, how='outer')

        # return the combined data for all interactors
        return all_data, warnings

class AlloSigma(Method):

    name = "AlloSigma"

    dataset_aliases = {"sub_cat": "active sites",
        "cofactor": "cofactor sites",
        "pockets": "pockets and interfaces"}

    exp_files = ['allosigma_mut.txt']

    directions = ['up', 'down']

    mutation_re = re.compile('[ACDEFGHIKLMNPQRSTVWY][0-9]+')
    filtered_re = re.compile(r"^filtered_(up|down)_(.+)\.tsv$")

    def _get_dataset_label(self, suffix):
        label = self.dataset_aliases.get(suffix, suffix.replace('_', ' '))
        return f"AlloSigMA {self.version} predicted consequence - {label}"

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

        # process dataframe
        filt = filt.reset_index()
        filt['mutations'] = filt['mutations'].str.split()
        filt = filt.explode('mutations')
        filt = filt.set_index('mutations')

        # remove columns not be considered, if present
        for c in ['n_mutations', 'avg_dG', 'index']:
            if c in filt.columns:
                filt = filt.drop(columns=c)

        # raise error if columns don't match the residue format
        forbidden_cols = []
        for c in filt.columns:
            if not re.fullmatch(self.mutation_re.pattern, c):
                forbidden_cols.append(c)
        if len(forbidden_cols) > 0:
            raise KeyError(f"column names {', '.join(forbidden_cols)} are not in the expected residue format")

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
            if not set(self.exp_files).issubset(set(allosigma2_files)):
                this_error = (f"the allosigma_mut.txt file must be present in the AlloSigma2 directory {allosigma2_dir}/{dirname}")
                raise MAVISpMultipleError(warning=warnings,
                                        critical=[MAVISpCriticalError(this_error)])

            found_data_files = [f for f in allosigma2_files if f.startswith("filtered_") and f.endswith(".tsv")]

            if len(found_data_files) < 1:
                this_error = (f"no filtered_*.tsv files were found in {dirname}; "
                f"found {', '.join(allosigma2_files)}")
                raise MAVISpMultipleError(warning=warnings,
                                            critical=[MAVISpCriticalError(this_error)])

            # parse the mutation file
            try:
                all_mut = self._parse_allosigma2_mutation_file(os.path.join(allosigma2_dir, dirname, 'allosigma_mut.txt'))
            except Exception as e:
                this_error = f"Exception {type(e).__name__} occurred when parsing the csv files. Arguments:{e.args}"
                raise MAVISpMultipleError(warning=warnings,
                                          critical=[MAVISpCriticalError(this_error)])

            available_suffixes = []

            # identify suffixes in the directory
            for fname in sorted(found_data_files):
                m = self.filtered_re.fullmatch(fname)
                if m is None:
                    continue
                suffix = m.group(2)
                if suffix not in available_suffixes:
                    available_suffixes.append(suffix)

            # parse each discovered suffix
            for suffix in available_suffixes:
                data = {}

                for direction in self.directions:
                    fname = f"filtered_{direction}_{suffix}.tsv"
                    if fname in allosigma2_files:
                        try:
                            data[direction] = self._parse_allosigma2_energy_table(
                                os.path.join(allosigma2_dir, dirname, fname))
                        except Exception as e:
                            this_error = (
                                f"Exception {type(e).__name__} occurred when parsing {fname}. "
                                f"Arguments:{e.args}")
                            raise MAVISpMultipleError(
                                warning=warnings,
                                critical=[MAVISpCriticalError(this_error)])
                    else:
                        data[direction] = None

                all_mut[suffix] = all_mut.apply(
                    self._process_allosigma2_tables,
                    filt_up=data["up"],
                    filt_down=data["down"],
                    cutoff=2,
                    axis=1)

            all_mut = all_mut[["mutations", "allosigma-mode"] + available_suffixes]
            all_mut = all_mut.set_index("mutations")
            
            allosigma2_data.append(all_mut)

        out_data = pd.concat(allosigma2_data, axis=0)

        rename_scheme = {suffix: self._get_dataset_label(suffix) for suffix in available_suffixes}
        rename_scheme["allosigma-mode"] = f"AlloSigMA {self.version} mutation type"

        return out_data.rename(columns=rename_scheme).fillna("-"), warnings

class RaSP(Method):

    unit = "kcal/mol"
    type = "Stability"
    conf_filename = 'damaging_proportion.csv'

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

    def _parse_conformation_csv(self, fname, warnings):

        warnings = []

        try:
            df = pd.read_csv(fname, usecols= ['mutation_interest', 'fraction_destabilizing'])
            df.rename(columns={'mutation_interest': 'variant'}, inplace=True)
            df = df.set_index('variant')
            df.columns = [f"{self.version} proportion of damaging conformations"]
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the RaSP conformation csv file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        return df

    def parse(self, dir_path):

        warnings = []

        rasp_files = os.listdir(dir_path)

        if len(rasp_files) == 1 and os.path.isfile(os.path.join(dir_path, rasp_files[0])):

            rasp_file = rasp_files[0]

            avg_mutation_data = self._parse_postprocessed_csv(os.path.join(dir_path, rasp_file), warnings)
            std_mutation_data = None
            conf_mutation_data = None 


        else:
            csv_files = []
            rasp_folder = os.listdir(dir_path)

            if self.conf_filename not in rasp_folder:
                                warnings.append(MAVISpWarningError("Conformation-dependent stability population file not found for RaSP data"))
                                conf_mutation_data = None

            # Check if all available files are directories
            for item in rasp_folder:
                if item == self.conf_filename:
                    conf_mutation_data = self._parse_conformation_csv(os.path.join(dir_path, self.conf_filename), warnings)
                    continue

                if not os.path.isdir(os.path.join(dir_path, item)):
                    this_error = f"{item} in {dir_path} is not a directory"
                    raise MAVISpMultipleError(warning=warnings,
                                                critical=[MAVISpCriticalError(this_error)])

                ddg_files = os.listdir(os.path.join(dir_path, item))
                # check one file per directory is available
                if len(ddg_files) != 1:
                    this_error = f"zero or multiples files found in {dir_path}; only one expected"
                    raise MAVISpMultipleError(warning=warnings,
                                              critical=[MAVISpCriticalError(this_error)])

                ddg_file = ddg_files[0]

                if not os.path.splitext(ddg_file)[-1] == '.csv':
                    this_error = f"file {ddg_file} in {dir_path}/{item} doesn't have csv as extension"
                    raise MAVISpMultipleError(warning=warnings,
                                              critical=[MAVISpCriticalError(this_error)])

                csv_files.append(os.path.join(dir_path, item, ddg_file))

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

            avg_ddg_colname = f'{self.type} ({self.unit}'
            std_ddg_colname = f'{self.type} ({self.unit}, st. dev.)'

            mutation_data_mean = mutation_data.mean(axis=1)
            mutation_data_std = mutation_data.std(axis=1)

            mutation_data[avg_ddg_colname] = mutation_data_mean
            mutation_data[std_ddg_colname] = mutation_data_std

            mutation_data = mutation_data[[avg_ddg_colname, std_ddg_colname]]

            mutation_data = mutation_data.sort_index()

            avg_mutation_data = mutation_data[[ c for c in mutation_data.columns if not 'st. dev.' in c ]]
            std_mutation_data = mutation_data[[ c for c in mutation_data.columns if ', st. dev.)' in c ]]

        return avg_mutation_data, std_mutation_data, conf_mutation_data, warnings


class ThermoMPNN(Method):

    unit = "kcal/mol"
    type = "Stability"

    def _parse_postprocessed_csv(self, fname, warnings):

        try:
            mutation_data = pd.read_csv(fname,
                                        usecols=['variant', 'ThermoMPNN_ddG'],
                                        index_col='variant')
        except Exception as e:
            this_error = (
                f"Exception {type(e).__name__} occurred when parsing the "
                f"ThermoMPNN csv file. Arguments: {e.args}"
            )
            raise MAVISpMultipleError(
                warning=warnings,
                critical=[MAVISpCriticalError(this_error)]
            )

        return mutation_data

    def parse(self, dir_path):

        warnings = []

        thermompnn_files = os.listdir(dir_path)
        conf_mutation_data = None 
        
        if len(thermompnn_files) == 1 and os.path.isfile(os.path.join(dir_path, thermompnn_files[0])):

            thermompnn_file = thermompnn_files[0]

            avg_mutation_data = self._parse_postprocessed_csv(os.path.join(dir_path, thermompnn_file), warnings)
            std_mutation_data = None

        else:
            csv_files = []
            thermompnn_folder = os.listdir(dir_path)

            for folder in thermompnn_folder:
                if not os.path.isdir(os.path.join(dir_path, folder)):
                    this_error = f"{folder} in {dir_path} is not a directory"
                    raise MAVISpMultipleError(warning=warnings,
                                              critical=[MAVISpCriticalError(this_error)])

                ddg_files = os.listdir(os.path.join(dir_path, folder))
                if len(ddg_files) != 1:
                    this_error = f"Zero or multiple files found in {os.path.join(dir_path, folder)}; only one expected"
                    raise MAVISpMultipleError(warning=warnings,
                                              critical=[MAVISpCriticalError(this_error)])

                ddg_file = ddg_files[0]

                if not os.path.splitext(ddg_file)[-1] == '.csv':
                    this_error = f"File {ddg_file} in {dir_path}/{folder} is not a CSV file."
                    raise MAVISpMultipleError(warning=warnings,
                                              critical=[MAVISpCriticalError(this_error)])

                csv_files.append(os.path.join(dir_path, folder, ddg_file))

            list_mutation_label = None
            mutation_data = None

            for fname in csv_files:
                tmp = self._parse_postprocessed_csv(fname, warnings)

                if list_mutation_label is None:
                    list_mutation_label = set(tmp.index)
                elif list_mutation_label != set(tmp.index):
                    this_error = "The mutation labels are not the same in the different csv files."
                    raise MAVISpMultipleError(warning=warnings,
                                              critical=[MAVISpCriticalError(this_error)])

                if mutation_data is None:
                    mutation_data = tmp
                else:
                    mutation_data = mutation_data.join(tmp, rsuffix="_")

            std_ddg_colname = f'{self.type} ({self.unit}, st. dev.)'
            avg_ddg_colname = f'{self.type} ({self.unit})'

            mutation_data_mean = mutation_data.mean(axis=1)
            mutation_data_std  = mutation_data.std(axis=1)

            mutation_data[avg_ddg_colname] = mutation_data_mean
            mutation_data[std_ddg_colname] = mutation_data_std

            mutation_data = mutation_data[[avg_ddg_colname, std_ddg_colname]]

            mutation_data = mutation_data.sort_index()

            avg_mutation_data = mutation_data[[c for c in mutation_data.columns if 'st. dev.' not in c]]
            std_mutation_data = mutation_data[[c for c in mutation_data.columns if ', st. dev.)' in c]]

        return avg_mutation_data, std_mutation_data, conf_mutation_data, warnings
