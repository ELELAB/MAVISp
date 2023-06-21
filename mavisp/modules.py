# MAVISp - classes for handling different MAVISp modules
# Copyright (C) 2022 Matteo Tiberti, Danish Cancer Society
#           (C) 2023 Jérémy Vinhas, Danish Cancer Society
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import os
import pandas as pd
from mavisp.methods import *
from mavisp.utils import three_to_one
import logging as log
import re
import pkg_resources

class DataType(object):
    def __init__(self, data_dir=None, stop_at='critical'):

        self.data_dir = data_dir
        self.data = None

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

    def get_dataset_view(self):
        return self.data

class MultiMethodDataType(DataType):
    def __init__(self, data_dir=None):

        super().__init__(data_dir)

    def ingest(self, mutations):

        self.data = pd.DataFrame({'mutations' : mutations}).set_index('mutations')

        warnings = []

        method_dirs = os.listdir(os.path.join(self.data_dir, self.module_dir))
        method_dirs = [ d for d in method_dirs if os.path.isdir(os.path.join(self.data_dir, self.module_dir, d)) ]

        if not set(method_dirs).issubset(set(self.methods.keys())):
            this_error = f"One or more {self.name} methods are not supported"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        for method_dir in method_dirs:
            data, this_warnings = self.methods[method_dir].parse(os.path.join(self.data_dir, self.module_dir, method_dir))
            self.data = self.data.join(data)
            warnings += this_warnings

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])

class Stability(MultiMethodDataType):

    module_dir = "stability"
    name = "stability"
    methods = {'foldx5'                      : MutateXStability(version="FoldX5"),
               'rosetta_cartddg2020_ref2015' : RosettaDDGPredictionStability(version='Rosetta Cartddg2020'),
               'rosetta_ref2015'             : RosettaDDGPredictionStability(version='Rosetta Cartddg'),
               'rasp'                        : RaSP(version='RaSP')}

    def ingest(self, mutations):

        self.data = pd.DataFrame({'mutations' : mutations}).set_index('mutations')

        warnings = []

        this_error = "Stability folder has to contain only 1 dataset"
        tmp = os.listdir(os.path.join(self.data_dir, self.module_dir))
        if len(tmp) != 1:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        structure_ID, residue_range = tmp[0].split("_", maxsplit=1)

        # tmp is equal to the list of methods
        tmp = os.listdir(os.path.join(self.data_dir, self.module_dir, f'{structure_ID}_{residue_range}'))

        # We loop over the methods
        for method in tmp:
            # all_methods is a list of all the methods (Rosetta,FoldX,RaSP) in the folder, it is defined for every method (nmr, cabsflex...), so it is reset to an empty list and we can check if there are no duplicated methods (FoldX,Rosetta).
            all_methods = []
            models = os.listdir(os.path.join(self.data_dir, self.module_dir, f'{structure_ID}_{residue_range}', method))
            model_data = pd.DataFrame({'mutations' : mutations}).set_index('mutations')
            model_data_list = []
        # check that:
            # all methods are supported
            # there are no duplicated methods (possible since they are in different model dirs)
            for model in models:
                method_dirs = os.listdir(os.path.join(self.data_dir, self.module_dir, f'{structure_ID}_{residue_range}', method, model))
                if not set(method_dirs).issubset(set(self.methods.keys())):
                    this_error = f"One or more {self.name} methods are not supported"
                    raise MAVISpMultipleError(warning=warnings,
                                            critical=[MAVISpCriticalError(this_error)])
                all_methods.extend(method_dirs)

                if len(all_methods) != len(set(all_methods)):
                                this_error = f"Only using one single instance of any given method is supported"
                                raise MAVISpMultipleError(warning=warnings,
                                                        critical=[MAVISpCriticalError(this_error)])
            for model in models:
                method_dirs = os.listdir(os.path.join(self.data_dir, self.module_dir, f'{structure_ID}_{residue_range}', method, model))

                for method_dir in method_dirs:
                    model_data, this_warnings = self.methods[method_dir].parse(os.path.join(self.data_dir, self.module_dir, f'{structure_ID}_{residue_range}', method, model, method_dir))
                    warnings += this_warnings
                    model_data.columns = [ f"Stability ({self.methods[method_dir].version}, {method}, {self.methods[method_dir].unit})" ]
                    model_data_list.append(model_data)
                    model_data = pd.concat(model_data_list, axis=1)

            keys = [ k for k in model_data.columns if k.startswith('Stability') ]

            if any(['FoldX' in k for k in keys]):
                foldx_col = [k for k in keys if 'FoldX' in k]
                assert foldx_col is not None
                foldx_header = foldx_col[0]
            else:
                foldx_header = None

            if any(['Rosetta' in k for k in keys]):
                rosetta_col = [k for k in keys if 'Rosetta' in k]
                assert rosetta_col is not None
                rosetta_header = rosetta_col[0]
            else:
                rosetta_header = None

            if any(['RaSP' in k for k in keys]):
                rasp_col = [k for k in keys if 'RaSP' in k]
                assert rasp_col is not None
                rasp_header = rasp_col[0]
            else:
                rasp_header = None

            # check if we have both FoldX and Rosetta/RaSP col
            if rosetta_header is not None and foldx_header is not None:
                model_data[f'Stability classification, {method}, (Rosetta, FoldX)'] = model_data.apply(self._generate_stability_classification, foldx_header=foldx_header, rosetta_header=rosetta_header, axis=1)
            else:
                warnings.append(MAVISpWarningError(f"Stability classification (Rosetta, FoldX) for {method} method can only be calculated if exactly one Rosetta and one MutateX datasets are available"))

            if rasp_header is not None and foldx_header is not None:
                model_data[f'Stability classification, {method}, (RaSP, FoldX)'] = model_data.apply(self._generate_stability_classification, foldx_header=foldx_header, rosetta_header=rasp_header, axis=1)
            else:
                warnings.append(MAVISpWarningError(f"Stability classification (RaSP, FoldX) for {method} method can only be calculated if exactly one RaSP and one MutateX datasets are available"))

            self.data = self.data.join(model_data)

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])

    def _generate_stability_classification(self, row, foldx_header, rosetta_header):

        stab_co = 3.0
        neut_co = 2.0

        if pd.isna(row[foldx_header]) or pd.isna(row[rosetta_header]):
            return pd.NA

        if row[foldx_header] >= stab_co and row[rosetta_header] >= stab_co:
            return 'Destabilizing'
        if row[foldx_header] <= (- stab_co) and row[rosetta_header] <= (- stab_co):
            return 'Stabilizing'
        if (- neut_co) < row[foldx_header] < neut_co and (- neut_co) < row[rosetta_header] < neut_co:
            return 'Neutral'
        return 'Uncertain'

class LocalInteractions(MultiMethodDataType):
    module_dir = "local_interactions"
    name = "local_interactions"
    methods = {'foldx5'                      : MutateXBinding(version="FoldX5",
                                                              complex_status='heterodimer'),
               'rosetta_flexddg_talaris2014' : RosettaDDGPredictionBinding(version='Rosetta Talaris 2014',
                                                                           complex_status='heterodimer')}

    def ingest(self, mutations):

        warnings = []

        try:
            super().ingest(mutations)
        except MAVISpMultipleError as e:
            if len(e.critical) > 0:
                raise
        else:
            e = None

        module_dir_files = os.listdir(os.path.join(self.data_dir, self.module_dir))
        if 'sasa.rsa' not in module_dir_files:
            this_error = f"required sasa.rsa file not found in {self.module_dir}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        try:
            rsa = pd.read_fwf(os.path.join(self.data_dir, self.module_dir, 'sasa.rsa'),
                skiprows=4, skipfooter=4, header=None, widths=[4,4,1,4,9,6,7,6,7,6,7,6,7,6],
                names = ['entry', 'rest', 'chain', 'resn', 'all_abs', 'sas_all_rel', 'sas_sc_abs',
                'sas_sc_rel', 'sas_mc_abs', 'sas_mc_rel', 'sas_np_abs', 'sas_np_rel', 'sas_ap_abs',
                'sas_ap_rel'],
                usecols = ['resn', 'sas_sc_rel'],
                )
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the sasa.rsa file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                        critical=[MAVISpCriticalError(this_error)])

        rsa['resn'] = rsa['resn'].astype("string")
        rsa = rsa.set_index('resn')
        self.data['res_num'] = self.data.index.str[1:-1]
        self.data = self.data.join(rsa, on='res_num')

        common_interactors = set.intersection(*[ set(m.interactors) for k, m in self.methods.items() ])
        if np.any([ len(m.interactors) != len(common_interactors) for k, m in self.methods.items()]):
            warnings.append(MAVISpWarningError("All supported methods must be available to generate the classification for a specific interactor"))

        for ci in common_interactors:
            self.data[f'Local Int. classification ({ci})'] = self.data.apply(self._generate_local_interactions_classification, axis=1, ci=ci)

        self.data = self.data.drop(columns=['res_num', 'sas_sc_rel'])

        if e is None and len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])
        elif len(warnings) > 0:
            e.warning.extend(warnings)
            raise e

    def _generate_local_interactions_classification(self, row, ci, stab_co=1.0):

        colnames = [ f"{m.type} (Binding with {ci}, {m.complex_status}, {m.version}, {m.unit})" for k, m in self.methods.items() ]

        if np.any( [ pd.isna(row[h]) for h in colnames ] ):
            if row['sas_sc_rel'] >= 20:
                return 'Uncertain'
            return pd.NA
        if np.all( [ row[h] > stab_co for h in colnames ] ):
            return 'Destabilizing'
        if np.all( [ row[h] < (-stab_co) for h in colnames ] ):
            return 'Stabilizing'
        if np.all( [ (- stab_co) <= row[h] <= stab_co for h in colnames ] ):
            return 'Neutral'
        return 'Uncertain'

class LocalInteractionsDNA(MultiMethodDataType):
    module_dir = "local_interactions_DNA"
    name = "local_interactions_DNA"
    methods = {'foldx5' : MutateXDNABinding(version="FoldX5")}

    def ingest(self, mutations):

        warnings = []

        try:
            super().ingest(mutations)
        except MAVISpMultipleError as e:
            if len(e.critical) > 0:
                raise
        else:
            e = None

        module_dir_files = os.listdir(os.path.join(self.data_dir, self.module_dir))
        if 'sasa.rsa' not in module_dir_files:
            this_error = f"required sasa.rsa file not found in {self.module_dir}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        try:
            rsa = pd.read_fwf(os.path.join(self.data_dir, self.module_dir, 'sasa.rsa'),
                skiprows=4, skipfooter=4, header=None, widths=[4,4,1,4,9,6,7,6,7,6,7,6,7,6],
                names = ['entry', 'rest', 'chain', 'resn', 'all_abs', 'sas_all_rel', 'sas_sc_abs',
                'sas_sc_rel', 'sas_mc_abs', 'sas_mc_rel', 'sas_np_abs', 'sas_np_rel', 'sas_ap_abs',
                'sas_ap_rel'],
                usecols = ['resn', 'sas_sc_rel'],
                )
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the sasa.rsa file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                        critical=[MAVISpCriticalError(this_error)])

        rsa['resn'] = rsa['resn'].astype("string")
        rsa = rsa.set_index('resn')
        self.data['res_num'] = self.data.index.str[1:-1]
        self.data = self.data.join(rsa, on='res_num')

        keys = [ k for k in self.data.columns if k.startswith('Local Int. With DNA') ]

        if len(keys) != 1:
            warnings.append(MAVISpWarningError("Exactly one data column expected to calculate classification"))

        self.data['Local Int. classification With DNA'] = self.data.apply(self._generate_local_interactions_DNA_classification, axis=1)

        self.data = self.data.drop(columns=['res_num', 'sas_sc_rel'])

        if e is None and len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])
        elif len(warnings) > 0:
            e.warning.extend(warnings)
            raise e

    def _generate_local_interactions_DNA_classification(self, row):

        keys = [ k for k in row.keys() if k.startswith('Local Int. With DNA') ]

        if len(keys) != 1:
            return pd.NA

        stab_co =  1.0

        header = keys[0]

        if pd.isna(row[header]):
            if row['sas_sc_rel'] >= 20:
                return 'Uncertain'
            return pd.NA

        if row[header] > stab_co:
            return 'Destabilizing'
        elif row[header] < (- stab_co):
            return 'Stabilizing'
        else:
            return 'Neutral'

class LocalInteractionsHomodimer(LocalInteractions):
    module_dir = "local_interactions_homodimers"
    name = "local_interactions_homodimers"
    methods = {'foldx5'                      : MutateXBinding(version="FoldX5",
                                                              complex_status='homodimer'),
               'rosetta_flexddg_talaris2014' : RosettaDDGPredictionBinding(version='Rosetta Talaris 2014',
                                                                           complex_status='homodimer')}

class LongRange(MultiMethodDataType):

    module_dir = "long_range"
    name = "long_range"
    methods = {'allosigma2' : AlloSigma(version=2)}

class SAS(DataType):

    module_dir = "sas"
    name = "sas"

    def ingest(self, mutations):
        warnings = []

        sas_file = os.listdir(os.path.join(self.data_dir, self.module_dir))
        if len(sas_file) != 1:
            this_error = f"multiple or no files found in {sas_file}; only one expected"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        sas_file = sas_file[0]

        log.info(f"parsing sas file {sas_file}")

        try:
            rsa = pd.read_fwf(os.path.join(self.data_dir, self.module_dir, sas_file),
                skiprows=4, skipfooter=4, header=None, widths=[4,4,1,4,9,6,7,6,7,6,7,6,7,6],
                names = ['entry', 'rest', 'chain', 'resn', 'all_abs', 'sas_all_rel', 'sas_sc_abs',
                'sas_sc_rel', 'sas_mc_abs', 'sas_mc_rel', 'sas_np_abs', 'sas_np_rel', 'sas_ap_abs',
                'sas_ap_rel'],
                usecols = ['resn', 'sas_sc_rel'],
                )
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the sasa.rsa file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                        critical=[MAVISpCriticalError(this_error)])

        mut_resn = [ mut[1:-1] for mut in mutations ]
        df = pd.DataFrame({'mutation' : mutations, 'position_mutation' : mut_resn})

        rsa["resn"]= rsa["resn"].astype(str)
        rsa = rsa.set_index("resn")

        result = pd.merge(df, rsa, left_on="position_mutation", right_on="resn", how="left")
        result = result[['mutation', 'sas_sc_rel']]
        self.data = result.rename(columns={'mutation' : 'mutation',
                                           'sas_sc_rel' : 'Relative Side Chain Solvent Accessibility in wild-type'})
        self.data = self.data.set_index('mutation')

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])

class PTMs(DataType):

    module_dir = "ptm"
    name = "ptms"
    allowed_ptms = ['s', 'p', 'y']
    allowed_ptm_muts = {'S' : 's',
                        'T' : 'p',
                        'Y' : 'y'}
    protein_chain = 'A'
    slim_pattern = re.compile('.* \((CLV|DEG|DOC|LIG|MOD|TRG)_[A-Za-z0-9_-]+\), [0-9]+-[0-9]+, ELM')

    def _assign_regulation_class(self, row):
        ref = row.name[0]
        alt = row.name[-1]
        r_a = set([ref, alt]) # set of WT and mutation residue types
        S_T = set(['S', 'T']) # set of S and T residue types
        allowed_wt_res = self.allowed_ptm_muts.keys()

        # if site is not a known PTM
        if row['phosphorylation_site'] != 'P':
            # if wt residue is S,T,Y or mutation is NOT S->T or T->S
            # return uncertain, otherwise neutral
            if ref in allowed_wt_res and not r_a == S_T:
                return 'uncertain'
            else:
                return 'neutral'

        # if site is a known PTM,
        else:
            # if S to T or T to S, return neutral
            if r_a == S_T:
                return 'neutral'

            # otherwise, cases:
                # any mutation sas < 20% or
                # any T/S to Y or
                # any Y to T/S
                # then return uncertain
            elif row['sas_sc_rel'] < 20.0 or\
            (ref in S_T and alt == 'Y') or\
            (ref == 'Y' and (alt in S_T)):
                return 'uncertain'

            # otherwise, if mutation sas >= 20% return damaging
            elif row['sas_sc_rel'] >= 20.0:
                return 'damaging'

        return '???'

    def _assign_ddg_class(self, row, ddg_col_name, stab_co=3.0, neut_co=2.0):

        if pd.isna(row[ddg_col_name]):
            return pd.NA

        if row[ddg_col_name] > stab_co:
            return 'Destabilizing'
        elif row[ddg_col_name] < (- stab_co):
            return 'Stabilizing'
        elif (- neut_co) <= row[ddg_col_name] <= neut_co:
            return 'Neutral'
        else:
            return 'Uncertain'

    def _assign_stability_class(self, row, mut_col_name, ptm_col_name):

        # if WT residue is not PTMable, mutation is neutral
        if not row.name[0] in self.allowed_ptm_muts.keys():
            return 'neutral'

        # otherwise if site is not a known PTM, return NA
        elif row['phosphorylation_site'] != 'P':
            return pd.NA

        # otherwise if either classification is NA, return NA
        elif pd.isna(row[ptm_col_name]) or pd.isna(row[mut_col_name]):
            return pd.NA

        # otherwise if PTM OR mut are classified as Uncertain, return Uncertain
        elif row[ptm_col_name] == 'Uncertain' or row[mut_col_name] == 'Uncertain':
            return 'uncertain'

        # otherwise, if PTM and mut they are classified the same, return Neutral.
        # if they are not classified the same, return damaging
        elif row[mut_col_name] == row[ptm_col_name]:
            return 'neutral'
        else:
            return 'damaging'

        return '???'

    def _validate_elms(self, row):
        if pd.isna(row['linear_motif']):
            return True

        motifs = row['linear_motif'].split('|')
        for motif in motifs:
            if self.slim_pattern.match(motif) is None:
                return False

        return True

    def _mut_in_phospho_slim(self, row, phospho_slims):

        if pd.isna(row['linear_motif']):
            return False

        pslim_identifier = phospho_slims['ELM_Identifier'].tolist()

        slims = row['linear_motif'].split('|')
        descriptions = [ d.split(',')[0] for d in slims ]
        for desc in descriptions:
            if any(ident in desc for ident in pslim_identifier):
                return True
        return False

    def _assign_function_class(self, row, mut_col_name, ptm_col_name):

        # if WT residue is not PTMable, mutation is neutral
        if not row.name[0] in self.allowed_ptm_muts.keys():
            return 'neutral'

        # otherwise if site is not a known PTM, return NA
        elif row['phosphorylation_site'] != 'P':
            return pd.NA

        # otherwise if either classification is NA, if it is phospho-motif
        # and solvent exposed, return potentially_damaging
        elif pd.isna(row[ptm_col_name]) or pd.isna(row[mut_col_name]):
            if row['site_in_slim'] and row['sas_sc_rel'] >= 20:
                return 'potentially_damaging'
            else:
                return 'uncertain'

        # otherwise if PTM OR mut are classified as Uncertain, return Uncertain
        elif row[ptm_col_name] == 'Uncertain' or row[mut_col_name] == 'Uncertain':
            return 'uncertain'

        # otherwise, if PTM and mut they are classified the same, return Neutral.
        # if they are not classified the same, return damaging
        elif row[mut_col_name] == row[ptm_col_name]:
            return 'neutral'
        else:
            return 'damaging'

        return '???'

    def ingest(self, mutations):

        warnings = []

        expected_files = ['summary_stability.txt',
                          'sasa.rsa',
                          'metatable.csv']

        ptm_files = os.listdir(os.path.join(self.data_dir, self.module_dir))

        if not set(expected_files).issubset(ptm_files):
            this_error = f"required file(s) not found"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        try:
            ddg_stability = pd.read_csv(os.path.join(self.data_dir, self.module_dir, 'summary_stability.txt'),
                delim_whitespace=True,
                header=None,
                names=['mutation', 'ddg_avg', 'ddg_std', 'ddg_min', 'ddg_max', 'idx'])
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the summary.txt file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                        critical=[MAVISpCriticalError(this_error)])

        try:
            ddg_binding = pd.read_csv(os.path.join(self.data_dir, self.module_dir, 'summary_binding.txt'),
                delim_whitespace=True,
                header=None,
                names=['mutation', 'ddg_avg', 'ddg_std', 'ddg_min', 'ddg_max', 'idx'])
        except FileNotFoundError as e:
            ddg_binding = pd.DataFrame(columns=['mutation', 'ddg_avg', 'ddg_std', 'ddg_min', 'ddg_max', 'idx'])
            warnings.append(MAVISpWarningError(f"summary_binding.txt not found - changes in free energy will not be used to classify function"))
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the summary_binding.txt file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                        critical=[MAVISpCriticalError(this_error)])

        try:
            rsa_monomer = pd.read_fwf(os.path.join(self.data_dir, self.module_dir, 'sasa.rsa'),
                skiprows=4, skipfooter=4, header=None, widths=[4,4,1,4,9,6,7,6,7,6,7,6,7,6],
                names = ['entry', 'rest', 'chain', 'resn', 'all_abs', 'sas_all_rel', 'sas_sc_abs',
                'sas_sc_rel', 'sas_mc_abs', 'sas_mc_rel', 'sas_np_abs', 'sas_np_rel', 'sas_ap_abs',
                'sas_ap_rel'],
                usecols = ['rest', 'resn', 'sas_sc_rel'],
                index_col = 'resn').fillna(pd.NA)
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the sasa.rsa file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                        critical=[MAVISpCriticalError(this_error)])

        try:
            cancermuts = pd.read_csv(os.path.join(self.data_dir, self.module_dir, 'metatable.csv'))
        except:
            this_error = f"Exception {type(e).__name__} occurred when parsing the cancermuts.csv file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                        critical=[MAVISpCriticalError(this_error)])

        # check if all the required columns are present
        required_columns = ['aa_position', 'ref_aa', 'alt_aa', 'linear_motif']
        if not set(required_columns).issubset(cancermuts.columns):
            this_error = f"input Cancermuts table doesn't have all the required columns"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        # load phospho-SLiM data file
        fname = pkg_resources.resource_filename(__name__, 'data/phosphoSLiMs_15062023.csv')
        phospho_slims = pd.read_csv(fname)

        # create final table
        final_table = pd.DataFrame({'mutation' : mutations})
        final_table['position'] = final_table['mutation'].str[1:-1].astype(int)
        final_table = final_table.set_index('mutation')

        # join cancermuts info
        cancermuts['mutation']  = cancermuts['ref_aa']\
                                + cancermuts['aa_position'].astype(str)\
                                + cancermuts['alt_aa']

        # remove rows with no mutation defined
        cancermuts = cancermuts[~ pd.isna(cancermuts['mutation'])]
        cancermuts = cancermuts.set_index('mutation')

        # validate found ELMs (i.e. they must have a ELM ID)
        are_elms_valid = cancermuts.apply(self._validate_elms, axis=1).all()
        if not are_elms_valid:
            warnings.append(MAVISpWarningError(f"one or more ELM entry does not contain ELM ID; PTM function will be set to NA"))


        # define whether each site is in a phospho-SLiM
        cancermuts['site_in_slim'] = cancermuts.apply(self._mut_in_phospho_slim, phospho_slims=phospho_slims, axis=1)

        # join cancermuts table to final table
        final_table = final_table.join(cancermuts)

        # join RSA data
        rsa_monomer = rsa_monomer.drop(columns=['rest'])
        final_table = final_table.join(rsa_monomer, on='position')

        # process DDG information
        ddg_stability['ref'] = ddg_stability['mutation'].str[0]
        ddg_stability['alt'] = ddg_stability['mutation'].str[-1]
        ddg_stability['number'] = ddg_stability['mutation'].str[2:-1].astype(int)
        ddg_stability['mutation'] = ddg_stability['ref'] + ddg_stability['number'].astype(str) + ddg_stability['alt']

        ddg_binding['ref'] = ddg_binding['mutation'].str[0]
        ddg_binding['alt'] = ddg_binding['mutation'].str[-1]
        ddg_binding['number'] = ddg_binding['mutation'].str[2:-1].astype(int)
        ddg_binding['mutation'] = ddg_binding['ref'] + ddg_binding['number'].astype(str) + ddg_binding['alt']

        if ddg_stability.shape[0] != 0 and ddg_binding.shape[0] != 0 and not (ddg_binding['mutation'] == ddg_stability['mutation']).all():
            this_error = f"stability DDG summary has different residues or a different order of residues than binding DDG summary"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        is_ptm_stability = ddg_stability.apply(lambda r: r['alt'].islower(), axis=1)
        stability_ptm_muts    = ddg_stability[  is_ptm_stability]
        stability_cancer_muts = ddg_stability[~ is_ptm_stability]

        is_ptm_binding = ddg_binding.apply(lambda r: r['alt'].islower(), axis=1)
        binding_ptm_muts    = ddg_binding[  is_ptm_binding]
        binding_cancer_muts = ddg_binding[~ is_ptm_binding]


        # process DDG PTM information. We do only one check since we already
        # checked that binding and stability have the same mutations
        if not set(stability_ptm_muts['alt']).issubset(set(self.allowed_ptms)):
            warnings.append(MAVISpWarningError(f"in DDG summary files, some of the mutations were not to {', '.join(self.allowed_ptms)}. These will be filtered out."))
            cancer_muts = cancer_muts[cancer_muts.apply(lambda r: r['alt'] in self.allowed_ptms)]

        if not set(stability_ptm_muts.number.unique()).issubset(set(stability_cancer_muts.number.unique())):
            this_error = f"in DDG summary files, some cancer mutations had no corresponding PTM mutation"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        if not (stability_ptm_muts.apply(lambda r: self.allowed_ptm_muts[r['ref']] == r['alt'], axis=1)).all():
            this_error = f"in DDG summary files, the wild type residue and the wild type residue upon PTM didn't correspond"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        # further process DDG dataframes and join them to final table
        stability_ptm_muts = stability_ptm_muts.set_index('number')
        stability_ptm_muts = stability_ptm_muts[['ddg_avg']].rename(columns={'ddg_avg' : 'stability_ddg_ptm'})
        final_table = final_table.join(stability_ptm_muts, on='position')

        stability_cancer_muts = stability_cancer_muts[['mutation', 'ddg_avg']].rename(columns={'ddg_avg' : 'stability_ddg_mut'}).set_index('mutation')
        final_table = final_table.join(stability_cancer_muts, on='mutation')

        binding_ptm_muts = binding_ptm_muts.set_index('number')
        binding_ptm_muts = binding_ptm_muts[['ddg_avg']].rename(columns={'ddg_avg' : 'binding_ddg_ptm'})
        final_table = final_table.join(binding_ptm_muts, on='position')

        binding_cancer_muts = binding_cancer_muts[['mutation', 'ddg_avg']].rename(columns={'ddg_avg' : 'binding_ddg_mut'}).set_index('mutation')
        final_table = final_table.join(binding_cancer_muts, on='mutation')

        # calculate class for DDG values, stability and binding
        final_table['cancer_mut_stab_class'] = final_table.apply(self._assign_ddg_class, ddg_col_name='stability_ddg_mut', axis=1)
        final_table['ptm_stab_class']        = final_table.apply(self._assign_ddg_class, ddg_col_name='stability_ddg_ptm', axis=1)

        final_table['cancer_mut_bind_class'] = final_table.apply(self._assign_ddg_class, ddg_col_name='binding_ddg_mut', axis=1, stab_co=1.0, neut_co=1.0)
        final_table['ptm_bind_class']        = final_table.apply(self._assign_ddg_class, ddg_col_name='binding_ddg_ptm', axis=1, stab_co=1.0, neut_co=1.0)

        # generate final classification for each type
        final_table['regulation'] = final_table.apply(self._assign_regulation_class, axis=1)
        final_table['stability'] = final_table.apply(self._assign_stability_class,
                                                    mut_col_name='cancer_mut_stab_class',
                                                    ptm_col_name='ptm_stab_class',
                                                    axis=1)
        if are_elms_valid:
            final_table['function'] = final_table.apply(self._assign_function_class,
                                                        mut_col_name='cancer_mut_bind_class',
                                                        ptm_col_name='ptm_bind_class',
                                                        axis=1)
        else:
            final_table['function'] = pd.NA

        # final processing for output table
        final_table = final_table[[ 'phosphorylation_site', 'site_in_slim', 'sas_sc_rel',
                                    'stability_ddg_ptm', 'binding_ddg_ptm', 'binding_ddg_mut',
                                    'regulation', 'stability', 'function' ]]

        self.data = final_table.rename(columns={'phosphorylation_site' : "PTMs",
                                                'site_in_slim'         : "is site part of phospho-SLiM",
                                                'sas_sc_rel'           : "PTM residue SASA (%)" ,
                                                'stability_ddg_ptm'    : "Change in stability with PTM (FoldX5, kcal/mol)",
                                                'binding_ddg_mut'      : "Change in binding with mutation (FoldX5, kcal/mol)",
                                                'binding_ddg_ptm'      : "Change in binding with PTM (FoldX5, kcal/mol)",
                                                'regulation'           : "PTM effect in regulation",
                                                'stability'            : "PTM effect in stability" ,
                                                'function'             : "PTM effect in function"
                                                })

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                        critical=[])


class CancermutsTable(DataType):

    module_dir = "cancermuts"
    name = "cancermuts"

    def __init__(self, data_dir=None):

        super().__init__(data_dir)

    def _process_sources(self, row):

        manual_pattern = 'Manual annotations from (.+)\..+'

        sources = row['sources'].split(',')
        for i,s in enumerate(sources):
            match = re.search(manual_pattern, s)

            if match is None:
                continue
            else:
                sources[i] = match.groups()[0]

        return ",".join(sources)

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
        cancermuts['sources'] = cancermuts.apply(self._process_sources, axis=1)

        # check if all mutations are present in the cancermuts table
        available_mutations = cancermuts.index.intersection(mutations)
        if len(available_mutations) != len(mutations):
            warnings.append(MAVISpWarningError("Not all the annotated mutations were found in the Cancermuts table"))

        # keep rows with mutations only
        cancermuts = cancermuts.loc[available_mutations]

        # check if all the mutations have a defined REVEL score
        # if np.any(pd.isna(cancermuts['REVEL_score'])):
        #     warnings.append(MAVISpWarningError("One or more mutations in the Cancermuts table don't have an associated REVEL score"))

        # filter by column and pretty rename column names
        cancermuts = cancermuts[['gnomad_genome_af', 'gnomad_exome_af', 'REVEL_score', 'sources']]

        self.data = cancermuts.rename(columns={ 'gnomad_genome_af' : 'gnomAD genome allele frequency',
                                                'gnomad_exome_af'  : 'gnomAD exome allele frequency',
                                                'REVEL_score'      : 'REVEL score',
                                                'sources'          : 'Mutation sources' })

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])

class ClinVar(DataType):

    module_dir = "clinvar"
    name = "clinvar"
    found_fname = "variants_output.csv"
    missing_fname = "entry_not_found.csv"

    def ingest(self, mutations):
        warnings = []

        clinvar_files = os.listdir(os.path.join(self.data_dir, self.module_dir))

        if len(clinvar_files) not in [1, 2] or not set(clinvar_files).issubset(set([self.found_fname, self.missing_fname])):
            this_error = f"One or two files expected for ClinVar, they must be {self.found_fname} and/or {self.missing_fname}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])


        if self.found_fname not in clinvar_files:
            this_error = f"variants_output.csv expected for ClinVar"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        log.info(f"parsing ClinVar files")

        try:
            clinvar_found = pd.read_csv(os.path.join(self.data_dir, self.module_dir, self.found_fname), sep=';')
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the csv files. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        missing_entries_path = os.path.join(self.data_dir, self.module_dir, self.missing_fname)
        if os.path.isfile(missing_entries_path):
            try:
                clinvar_missing = pd.read_csv(missing_entries_path, sep=';')
            except Exception as e:
                this_error = f"Exception {type(e).__name__} occurred when parsing the csv files. Arguments:{e.args}"
                raise MAVISpMultipleError(warning=warnings,
                                        critical=[MAVISpCriticalError(this_error)])
            missing_muts = clinvar_missing['variant_name'].to_list()
        else:
            warnings.append(MAVISpWarningError(f"file {self.missing} not found in ClinVar module"))
            missing_muts = []

        if len(set(missing_muts).intersection(set(clinvar_found['variant_name']))) > 0:
            warnings.append(MAVISpWarningError(f"some mutations are both in the ClinVar annotation and not found file"))

        if not ('clinvar_code' in clinvar_found.columns) ^ ('variant_id' in clinvar_found.columns):
            this_error = f"The variants_output.csv file must contain either the clinvar_code or the variant_id column"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        if 'clinvar_code' in clinvar_found.columns:
            warnings.append(MAVISpWarningError(f"the input file has the old style clinvar_code column"))
            id_col = 'clinvar_code'
        else:
            id_col = 'variant_id'

        if 'interpretation' not in clinvar_found.columns:
            this_error = f"The variants_output.csv file must contain the interpretation column"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])
        clinvar_found[id_col] = clinvar_found[id_col].astype(str)
        if "number_of_stars" in clinvar_found.columns:
            clinvar_found['number_of_stars'] = clinvar_found['number_of_stars'].astype(str)
            clinvar_found = clinvar_found.groupby('variant_name').agg(lambda x: ", ".join(list(x)))[[id_col, 'interpretation',"number_of_stars"]]
            self.data = clinvar_found.rename({ id_col           : 'ClinVar Variation ID',
                                               'interpretation' : 'ClinVar Interpretation',
                                               'number_of_stars': 'ClinVar Review Status'}, axis=1)
        else:
            warnings.append(MAVISpWarningError(f"the variant_output.csv file doesn't contain the number_of_stars column (ClinVar review status)"))
            clinvar_found[id_col] = clinvar_found[id_col].astype(str)
            clinvar_found = clinvar_found.groupby('variant_name').agg(lambda x: ", ".join(list(x)))[[id_col, 'interpretation']]
            self.data = clinvar_found.rename({ id_col        : 'ClinVar Variation ID',
                                            'interpretation' : 'ClinVar Interpretation',}, axis=1)

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])

class EVE(DataType):

    module_dir = "eve"
    name = "eve"

    def __init__(self, data_dir=None):

        super().__init__(data_dir)

    def _process_sources(self, row):

        manual_pattern = 'Manual annotations from (.+)\..+'

        sources = row['sources'].split(',')
        for i,s in enumerate(sources):
            match = re.search(manual_pattern, s)

            if match is None:
                continue
            else:
                sources[i] = match.groups()[0]

        return ",".join(sources)

    def ingest(self, mutations):
        warnings = []

        # check that we have only one file
        eve_files = os.listdir(os.path.join(self.data_dir, self.module_dir))
        if len(eve_files) != 1:
            this_error = f"multiple files found in {eve_files}; only one expected"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        eve_file = eve_files[0]

        log.info(f"parsing EVE file {eve_file}")

        # parse EVE table
        try:
            eve = pd.read_csv(os.path.join(self.data_dir, self.module_dir, eve_file))
        except:
            this_error = f"Failed parsing EVE csv file {eve_file}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        # check if all the required columns are present
        required_columns = ['protein_name', 'mutations', 'evol_indices', 'EVE_scores', 'EVE_classes_100_pct_retained', 'uncertainty']
        if not set(required_columns).issubset(eve.columns):
            this_error = f"input table doesn't have all the required columns"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        # process table
        eve = eve[['mutations', 'EVE_scores', 'EVE_classes_75_pct_retained']]
        eve = eve.set_index('mutations')

        self.data = eve.rename(columns={ 'EVE_scores' : 'EVE score',
                                         'EVE_classes_75_pct_retained' : 'EVE classification (25% Uncertain)'})

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])


class AlphaFoldMetadata(DataType):

    module_dir = "alphafold"
    name = "alphafold"

    def ingest(self, mutations):
        warnings = []

        af_files = os.listdir(os.path.join(self.data_dir, self.module_dir))
        if len(af_files) != 1:
            this_error = f"multiple or no files found in {af_files}; only one expected"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        af_file = af_files[0]

        log.info(f"parsing AlphaFold metadata file {af_file}")

        try:
            afmd = pd.read_csv(os.path.join(self.data_dir, self.module_dir, af_file))
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the csv files. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        if not set(['resnum', 'resname', 'pLDDT', 'secstruc']).issubset(set(afmd.columns)):
            this_error = f"The CSV file must contain the following columns: resnum, resname, pLDDT, secstruc"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        afmd['resname'] = afmd.apply(lambda r: three_to_one[r['resname']], axis=1)
        nres = afmd.shape[0]
        afmd = afmd.iloc[np.arange(nres).repeat(len(three_to_one))]
        afmd['alt_aa'] = list(three_to_one.values()) * nres
        afmd['mutations'] = afmd['resname'] + afmd['resnum'].astype(str) + afmd['alt_aa']
        afmd = afmd.set_index('mutations')
        afmd = afmd[['pLDDT', 'secstruc']]
        afmd = afmd.rename(columns={'pLDDT'    : 'AlphaFold2 model pLDDT score',
                                    'secstruc' : 'AlphaFold2 model secondary structure'})

        if not set(mutations).issubset(set(afmd.index)):
            diff = set(mutations).difference(set(afmd.index))
            warnings.append(MAVISpWarningError(f"The following mutations had a wild-type residue that did not correspond to the wild-type residue in the AF model: {diff}"))

        self.data = afmd

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])

class DeMaSk(DataType):

    module_dir = "demask"
    name = "demask"

    def _classify(self, row):
        if row['score'] > 0:
            return 'gain_of_function'
        elif row['score'] < 0:
            return 'loss_of_function'
        else:
            return 'neutral'

    def ingest(self, mutations):

        warnings = []

        demask_files = os.listdir(os.path.join(self.data_dir, self.module_dir))
        if len(demask_files) != 1:
            this_error = f"multiple or no files found in {demask_files}; only one expected"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        demask_file = demask_files[0]

        log.info(f"parsing DeMaSk data file {demask_file}")

        try:
            demask = pd.read_csv(os.path.join(self.data_dir, self.module_dir, demask_file), delim_whitespace=True)
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the csv files. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        if not set(['pos', 'WT', 'var', 'score', 'entropy', 'log2f_var', 'matrix']).issubset(set(demask.columns)):
            this_error = f"The input file doesn't have the columns expected for a DeMaSk output file"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        demask['mutations'] = demask['WT'] + demask['pos'].astype(str) + demask['var']
        demask = demask[['mutations', 'score', 'entropy', 'log2f_var']]
        demask = demask.set_index('mutations')
        demask['DeMaSk predicted consequence'] = demask.apply(self._classify, axis=1)

        self.data = demask.rename(columns = {'score'     : 'DeMaSk delta fitness',
                                             'entropy'   : 'DeMaSk Shannon entropy',
                                             'log2f_var' : 'DeMaSk log2 variant frequency'})
        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])

class GEMME(DataType):

    module_dir = "gemme"
    name = "gemme"
    accepted_filenames = ['normPred_evolCombi.txt']

    def ingest(self, mutations):

        warnings = []

        gemme_files = os.listdir(os.path.join(self.data_dir, self.module_dir))
        if len(gemme_files) != 1:
            this_error = f"multiple or no files found in {gemme_files}; only one expected"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        gemme_file = gemme_files[0]

        if gemme_file not in self.accepted_filenames:
            this_error = f"the input file for GEMME must be named {gemme_file}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])


        log.info(f"parsing GEMME data file {gemme_file}")

        try:
            gemme = pd.read_csv(os.path.join(self.data_dir, self.module_dir, gemme_file),
                                sep=' ', index_col=None)
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the csv files. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        # rename columns and rows as better names
        gemme = gemme.rename(columns={h:h[1:] for h in gemme.columns},
                             index={r:r.upper() for r in gemme.index})

        # melt matrix to one mutation per row
        gemme = gemme.reset_index().melt('index').rename(columns={'index':'mut', 'variable':'res', 'value':'score'})

        # calculate which residue is WT for each row (for every residue, this would
        # be tone with score 'None')
        wts = gemme.groupby('res').apply(lambda x: x[pd.isna(x['score'])]['mut'].to_list()[0])
        wts.name = 'wt'

        # join WT definition on main dataframe
        gemme = gemme.join(wts, on='res')

        # reconstruct mutations in the usual format
        gemme['mutations'] = gemme['wt'] + gemme['res'] + gemme['mut']

        # drop unnecessary columns and rename for pretty
        gemme = gemme.drop(columns=['wt', 'res', 'mut']).rename(columns={'score': 'GEMME Score'})
        self.data = gemme.set_index('mutations')

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])
