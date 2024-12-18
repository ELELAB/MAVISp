# MAVISp - classes for handling different MAVISp modules
# Copyright (C) 2022 Matteo Tiberti, Danish Cancer Society
#           (C) 2023 Jérémy Vinhas, Danish Cancer Society
#           (C) 2024 Pablo Sánchez-Izquierdo, Danish Cancer Society
#           (C) 2024 Eleni Kiahaki, Danish Cancer Society
#           (C) 2024 Karolina Krzesińska, Danish Cancer Society & DTU
#
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
from mavisp.utils import three_to_one, three_to_one_hgvsp
from collections import defaultdict
import logging as log
import re
import pkg_resources
import numbers
import yaml

class MavispModule(object):
    def __init__(self, data_dir=None, module_dir=None, stop_at='critical'):

        self.data_dir = data_dir
        self.data = None

        if module_dir is not None:
            self.module_dir = module_dir

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

class MavispMultiEnsembleModule(MavispModule):
    def __init_subclass__(cls, module_class, **kwargs):
        super().__init_subclass__(**kwargs)

        cls.base_module = module_class

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.ensembles = {}

        ensembles = os.listdir(os.path.join(self.data_dir, self.module_dir))

        for ensemble in ensembles:
            base_module_dir = os.path.join(self.module_dir, ensemble)

            self.ensembles[ensemble] = self.base_module(self.data_dir,
                                                        module_dir=base_module_dir)

    def ingest(self, mutations):

        self.data = pd.DataFrame({'mutations' : mutations}).set_index('mutations')

        warnings = []

        if len(self.ensembles) == 0:
            message = "module contained no ensembles"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(message)])

        for name, obj in self.ensembles.items():
            try:
                obj.ingest(mutations)
            except MAVISpMultipleError as e:
                ensemble_warnings = [MAVISpWarningError(f"[{name}] {exc.args[0]}") for exc in e.warning]
                warnings += ensemble_warnings
                if len(e.critical) != 0:
                    raise MAVISpMultipleError(warning=warnings,
                                              critical=e.critical)

            new_colnames = {c : f"{c} [{name}]" for c in obj.data.columns }
            obj.data = obj.data.rename(columns=new_colnames)

            self.data = self.data.join(obj.data)

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                    critical=list())

class MultiMethodMavispModule(MavispModule):
    def __init__(self, data_dir=None, module_dir=None):

        super().__init__(data_dir, module_dir)

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

class Stability(MultiMethodMavispModule):

    module_dir = "stability"
    name = "stability"
    methods = {'foldx5'                      : MutateXStability(version="FoldX5"),
               'rosetta_cartddg2020_ref2015' : RosettaDDGPredictionStability(version='Rosetta Cartddg2020'),
               'rosetta_ref2015'             : RosettaDDGPredictionStability(version='Rosetta Cartddg'),
               'rasp'                        : RaSP(version='RaSP')}

    def ingest(self, mutations):

        self.data = pd.DataFrame({'mutations' : mutations}).set_index('mutations')

        warnings = []

        tmp = os.listdir(os.path.join(self.data_dir, self.module_dir))

        if len(tmp) != 1:
            this_error = "Stability folder has to contain only 1 dataset"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        structure_ID, residue_range = tmp[0].split("_", maxsplit=1)

        # tmp is equal to the list of methods
        tmp = os.listdir(os.path.join(self.data_dir, self.module_dir, f'{structure_ID}_{residue_range}'))

        # all_methods is a list of all the methods (Rosetta,FoldX,RaSP) in the folder, it is defined for every method (nmr, cabsflex...), so it is reset to an empty list and we can check if there are no duplicated methods (FoldX,Rosetta).
        all_methods = []
        models = os.listdir(os.path.join(self.data_dir, self.module_dir, f'{structure_ID}_{residue_range}'))

        model_data_stds = None
        model_data_list = []
        model_data_stds_list = []
    # check that:
        # all methods are supported
        # there are no duplicated methods (possible since they are in different model dirs)
        for model in models:
            method_dirs = os.listdir(os.path.join(self.data_dir, self.module_dir, f'{structure_ID}_{residue_range}', model))
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
            method_dirs = os.listdir(os.path.join(self.data_dir, self.module_dir, f'{structure_ID}_{residue_range}', model))

            for method_dir in method_dirs:
                model_averages, model_stds, this_warnings = self.methods[method_dir].parse(os.path.join(self.data_dir, self.module_dir, f'{structure_ID}_{residue_range}', model, method_dir))

                warnings += this_warnings
                model_averages.columns = [ f"Stability ({self.methods[method_dir].version}, {self.methods[method_dir].unit})" ]

                model_data_list.append(model_averages)
                model_data = pd.concat(model_data_list, axis=1)

                if model_stds is not None:
                    model_stds.columns = [ f"Stability ({self.methods[method_dir].version}, {self.methods[method_dir].unit}, st. dev.)" ]

                    model_data_stds_list.append(model_stds)
                    model_data_stds = pd.concat(model_data_stds_list, axis=1)

        keys = [ k for k in model_data.columns if k.startswith('Stability') and not 'st. dev.' in k ]

        if any(['FoldX' in k for k in keys]):
            foldx_col = [k for k in keys if 'FoldX' in k and 'st. dev.' not in k]
            assert foldx_col is not None
            foldx_header = foldx_col[0]
        else:
            foldx_header = None

        if any(['Rosetta' in k for k in keys]):
            rosetta_col = [k for k in keys if 'Rosetta' in k and 'st. dev.' not in k]
            assert rosetta_col is not None
            rosetta_header = rosetta_col[0]
        else:
            rosetta_header = None

        if any(['RaSP' in k for k in keys]):
            rasp_col = [k for k in keys if 'RaSP' in k and 'st. dev.' not in k]
            assert rasp_col is not None
            rasp_header = rasp_col[0]
        else:
            rasp_header = None

        # check if we have both FoldX and Rosetta/RaSP col
        if rosetta_header is not None and foldx_header is not None:
            model_data[f'Stability classification, (Rosetta, FoldX)'] = model_data.apply(self._generate_stability_classification, foldx_header=foldx_header, rosetta_header=rosetta_header, axis=1)
        else:
            warnings.append(MAVISpWarningError(f"Stability classification (Rosetta, FoldX) can only be calculated if exactly one Rosetta and one MutateX datasets are available"))

        if rasp_header is not None and foldx_header is not None:
            model_data[f'Stability classification, (RaSP, FoldX)'] = model_data.apply(self._generate_stability_classification, foldx_header=foldx_header, rosetta_header=rasp_header, axis=1)
        else:
            warnings.append(MAVISpWarningError(f"Stability classification (RaSP, FoldX) can only be calculated if exactly one RaSP and one MutateX datasets are available"))

        self.data = self.data.join(model_data)
        if model_data_stds is not None:
            self.data = self.data.join(model_data_stds)

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

class SimpleStability(Stability):
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

                    model_data, model_stds, this_warnings = self.methods[method_dir].parse(os.path.join(self.data_dir, self.module_dir, f'{structure_ID}_{residue_range}', method, model, method_dir))
                    warnings += this_warnings

                    model_data.columns = [ f"Stability ({self.methods[method_dir].version}, {method}, {self.methods[method_dir].unit})" ]
                    model_data_list.append(model_data)

                    if model_stds is not None:
                        model_stds.columns = [ f"Stability ({self.methods[method_dir].version}, {method}, {self.methods[method_dir].unit}, st. dev.)" ]
                        model_data_list.append(model_stds)

            model_data = pd.concat(model_data_list, axis=1)

            keys = [ k for k in model_data.columns if k.startswith('Stability') and not 'st. dev.' in k]

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


class EnsembleStability(MavispMultiEnsembleModule, module_class=Stability):
    module_dir = "stability"
    name = "stability"

class LocalInteractions(MultiMethodMavispModule):
    module_dir = "local_interactions"
    name = "local_interactions"
    methods = {'foldx5'                      : MutateXBinding(version="FoldX5",
                                                              complex_status='heterodimer'),
               'rosetta_flexddg_talaris2014' : RosettaDDGPredictionBinding(version='Rosetta Talaris 2014',
                                                                           complex_status='heterodimer')}
    sas_filename = 'sasa.rsa'

    def _parse_sas(self, fname, warnings):

        try:
            rsa = pd.read_fwf(fname,
                skiprows=4, skipfooter=4, header=None, widths=[4,4,1,4,9,6,7,6,7,6,7,6,7,6],
                names = ['entry', 'rest', 'chain', 'resn', 'all_abs', 'sas_all_rel', 'sas_sc_abs',
                'sas_sc_rel', 'sas_mc_abs', 'sas_mc_rel', 'sas_np_abs', 'sas_np_rel', 'sas_ap_abs',
                'sas_ap_rel'],
                usecols = ['resn', 'sas_sc_rel']).fillna(pd.NA)
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the sasa.rsa file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                        critical=[MAVISpCriticalError(this_error)])

        rsa['resn'] = rsa['resn'].astype("string")

        return rsa.set_index('resn')

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

        rsa = self._parse_sas(os.path.join(self.data_dir, self.module_dir, self.sas_filename), warnings)

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

class TaccLocalInteractions(LocalInteractions):

    sas_filename = 'acc_REL.csv'

    def _parse_sas(self, sas_file, warnings):

        try:
            rsa = pd.read_csv(sas_file)
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the acc_REL.csv file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        rsa = rsa.drop(columns=['acc_std'])

        rsa = rsa.rename(columns={'residue'     : 'resn',
                                  'acc_average' : 'sas_sc_rel'})

        rsa['resn'] = rsa['resn'].astype(str)

        return rsa.set_index('resn')

class EnsembleLocalInteractions(MavispMultiEnsembleModule, module_class=TaccLocalInteractions):
    module_dir = "local_interactions"
    name = "local_interactions"

class LocalInteractionsDNA(MultiMethodMavispModule):
    module_dir = "local_interactions_DNA"
    name = "local_interactions_DNA"
    methods = {'foldx5' : MutateXDNABinding(version="FoldX5")}
    sas_filename = 'sasa.rsa'

    def _parse_sas(self, fname, warnings):

        try:
            rsa = pd.read_fwf(fname,
                skiprows=4, skipfooter=4, header=None, widths=[4,4,1,4,9,6,7,6,7,6,7,6,7,6],
                names = ['entry', 'rest', 'chain', 'resn', 'all_abs', 'sas_all_rel', 'sas_sc_abs',
                'sas_sc_rel', 'sas_mc_abs', 'sas_mc_rel', 'sas_np_abs', 'sas_np_rel', 'sas_ap_abs',
                'sas_ap_rel'],
                usecols = ['resn', 'sas_sc_rel']).fillna(pd.NA)
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the sasa.rsa file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                        critical=[MAVISpCriticalError(this_error)])

        rsa['resn'] = rsa['resn'].astype(str)

        return rsa.set_index('resn')

    def ingest(self, mutations):

        warnings = []

        try:
            super().ingest(mutations)
        except MAVISpMultipleError as e:
            if len(e.critical) > 0:
                raise
        else:
            e = None

        rsa = self._parse_sas(os.path.join(self.data_dir, self.module_dir, self.sas_filename), warnings)

        module_dir_files = os.listdir(os.path.join(self.data_dir, self.module_dir))

        self.data['res_num'] = self.data.index.str[1:-1]

        self.data = self.data.join(rsa, on='res_num')

        common_interactors = set.intersection(*[ set(m.interactors) for k, m in self.methods.items() ])

        for ci in common_interactors:
            self.data[f'Local Int. With DNA classification ({ci})'] = self.data.apply(self._generate_local_interactions_DNA_classification, axis=1, ci=ci)

        self.data = self.data.drop(columns=['res_num', 'sas_sc_rel'])

        if e is None and len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])
        elif len(warnings) > 0:
            e.warning.extend(warnings)
            raise e

    def _generate_local_interactions_DNA_classification(self, row, ci, stab_co=1.0):

        colnames = [ f"{m.type} ({ci}, {m.complex_status}, {m.version}, {m.unit})" for k, m in self.methods.items() ]

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

class TaccLocalInteractionsDNA(LocalInteractionsDNA):

    sas_filename = 'acc_REL.csv'

    def _parse_sas(self, sas_file, warnings):

        try:
            rsa = pd.read_csv(sas_file)
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the acc_REL.csv file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        rsa = rsa.drop(columns=['acc_std'])

        rsa = rsa.rename(columns={'residue'     : 'resn',
                                  'acc_average' : 'sas_sc_rel'})

        rsa['resn'] = rsa['resn'].astype(str)

        return rsa.set_index('resn')

class EnsembleLocalInteractionsDNA(MavispMultiEnsembleModule, module_class=TaccLocalInteractionsDNA):
    module_dir = "local_interactions_DNA"
    name = "local_interactions_DNA"

class LocalInteractionsHomodimer(LocalInteractions):
    module_dir = "local_interactions_homodimers"
    name = "local_interactions_homodimers"
    methods = {'foldx5'                      : MutateXBinding(version="FoldX5",
                                                              complex_status='homodimer'),
               'rosetta_flexddg_talaris2014' : RosettaDDGPredictionBinding(version='Rosetta Talaris 2014',
                                                                           complex_status='homodimer')}

class TaccLocalInteractionsHomodimer(TaccLocalInteractions):
    module_dir = "local_interactions_homodimers"
    name = "local_interactions_homodimers"
    methods = {'foldx5'                      : MutateXBinding(version="FoldX5",
                                                              complex_status='homodimer'),
               'rosetta_flexddg_talaris2014' : RosettaDDGPredictionBinding(version='Rosetta Talaris 2014',
                                                                           complex_status='homodimer')}

class EnsembleLocalInteractionsHomodimer(MavispMultiEnsembleModule, module_class=TaccLocalInteractionsHomodimer):
    module_dir = "local_interactions_homodimers"
    name = "local_interactions_homodimers"

class LongRange(MultiMethodMavispModule):

    module_dir = "long_range"
    name = "long_range"
    methods = {'allosigma2' : AlloSigma(version=2)}

class SAS(MavispModule):

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

class TaccSAS(MavispModule):

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
            rsa = pd.read_csv(os.path.join(self.data_dir, self.module_dir, sas_file))
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the sasa.rsa file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                        critical=[MAVISpCriticalError(this_error)])

        mut_resn = [ mut[1:-1] for mut in mutations ]
        df = pd.DataFrame({'mutation' : mutations, 'position_mutation' : mut_resn})

        rsa["residue"]= rsa["residue"].astype(str)
        rsa = rsa.set_index("residue")

        result = pd.merge(df, rsa, left_on="position_mutation", right_on="residue", how="left")
        result = result[['mutation', 'acc_average', 'acc_std']]
        self.data = result.rename(columns={'mutation' : 'mutation',
                                           'acc_average' : 'Relative Side Chain Solvent Accessibility in wild-type (average)',
                                           'acc_std'     : 'Relative Side Chain Solvent Accessibility in wild-type (standard deviation)'})

        self.data = self.data.set_index('mutation')

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])

class EnsembleSAS(MavispMultiEnsembleModule, module_class=TaccSAS):
    module_dir = "sas"
    name = "sas"

class EFoldMine(MavispModule):

    module_dir = "efoldmine"
    name = "efoldmine"

    def ingest(self, mutations):
        warnings = []
        efoldmine_file = os.listdir(os.path.join(self.data_dir, self.module_dir))
        if len(efoldmine_file) != 1:
            this_error = f"multiple or no files found in {efoldmine_file}; only one expected"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        efoldmine_file = efoldmine_file[0]

        log.info(f"parsing efoldmine file {efoldmine_file}")

        try:
            ef_res = pd.read_csv(os.path.join(self.data_dir, self.module_dir, efoldmine_file))
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the efoldmine file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                        critical=[MAVISpCriticalError(this_error)])

        # Filter required columns:
        efoldmine_parsed = ef_res[['residue_index', 'earlyFolding']]
        if 'residue_index' not in ef_res.columns or 'earlyFolding' not in ef_res.columns:
            this_error = "Required columns 'residue_index' or 'earlyFolding' are missing from the file."
            raise MAVISpMultipleError(warning=warnings, critical=[MAVISpCriticalError(this_error)])

        # Inspect early folding scores:
        try:
            pd.to_numeric(efoldmine_parsed['earlyFolding'], errors='raise')
        except ValueError:
            this_error = "Invalid or missing early folding scores in the file."
            raise MAVISpMultipleError(warning=warnings, critical=[MAVISpCriticalError(this_error)])

        # Precompute early folding regions:
        efoldmine_scores = efoldmine_parsed['earlyFolding'].tolist()
        seq_length = len(efoldmine_scores)
        window_size = 3
        thres = 0.169

        # Initialize early_foding_regions list:
        early_folding_regions = [False] * seq_length

        # Find the residues in early folding regions:
        for i in range(seq_length - window_size + 1):
            window = efoldmine_scores[i:i + window_size]
            if all(score > thres for score in window):
                early_folding_regions[i:i + window_size] = [True] * window_size

        efoldmine_parsed['is_early_folding'] = early_folding_regions

        # Compare with mutations:
        result = []

        for mut in mutations:
            mut_resn = int(mut[1:-1])
            row = efoldmine_parsed[efoldmine_parsed['residue_index']+1 == mut_resn]
            if len(row) != 1:
                this_error = f"Expected exactly one row for residue index {mut_resn}, but found {len(row)} rows."
                raise MAVISpMultipleError(warning=warnings, critical=[MAVISpCriticalError(this_error)])
            is_early_folding = row['is_early_folding'].iloc[0]
            efoldmine_score = row['earlyFolding'].iloc[0]
            result.append((is_early_folding, efoldmine_score))

        # Create DataFrame:
        df = pd.DataFrame(result, columns=['efoldmine_is_early_folding', 'efoldmine_score'], index=mutations)
        self.data = df.rename(columns={'efoldmine_is_early_folding' : 'EFoldMine - part of early folding region',
                                       'efoldmine_score'            : 'EFoldMine score'})

class DenovoPhospho(MavispModule):
    module_dir = "denovo_phospho"
    name = "denovo_phospho"
    expected_files = ['aggregated_filtered_output.csv', 'sasa.rsa']
    sasa_fname = expected_files[1]

    def _parse_sas(self, fname):

        return pd.read_fwf(fname,
            skiprows=4, skipfooter=4, header=None, widths=[4,4,1,4,9,6,7,6,7,6,7,6,7,6],
            names = ['entry', 'rest', 'chain', 'resn', 'all_abs', 'sas_all_rel', 'sas_sc_abs',
            'sas_sc_rel', 'sas_mc_abs', 'sas_mc_rel', 'sas_np_abs', 'sas_np_rel', 'sas_ap_abs',
            'sas_ap_rel'],
            usecols = ['resn', 'sas_sc_rel'],
            index_col = 'resn')

    def ingest(self, mutations):
        warnings = []
        mf_files = os.listdir(os.path.join(self.data_dir, self.module_dir))
        if not set(self.expected_files).issubset(set(mf_files)):
            this_error = f"the input files for Muts On Phospho must be named {', '.join(self.expected_files)}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        # Load and process aggregated data
        try:
            aggregated_df = pd.read_csv(os.path.join(self.data_dir, self.module_dir, 'aggregated_filtered_output.csv'))
        except Exception as e:
            this_error = f"Failed to load 'aggregated_filtered_output.csv': {e}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        # Parse SAS data
        try:
            sas_data = self._parse_sas(os.path.join(self.data_dir, self.module_dir, self.sasa_fname))
        except Exception as e:
            this_error = f"Failed to parse solvent accessibility data ('{self.sasa_fname}'): {e}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        try:
            # Data type alignment and merge
            sas_data = sas_data.reset_index()
            aggregated_df = pd.merge(aggregated_df, sas_data, left_on='resnum', right_on='resn', how='left')
            aggregated_df['restype_resnum_kinase'] = aggregated_df['restype'] + aggregated_df['resnum'].astype(str) + '_' + aggregated_df['kinase']
        except Exception as e:
            this_error = f"Error during data preparation: {e}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        gain_of_function = defaultdict(list)
        loss_of_function = defaultdict(list)

        # Main analysis loop
        try:
            available_mutations = list(set(mutations).intersection(set(aggregated_df.columns)))

            for mutation in available_mutations:
                for index, row in aggregated_df.iterrows():
                    restype_resnum_kinase = row['restype_resnum_kinase']
                    if pd.isna(row['WT']) and not pd.isna(row[mutation]) and row['sas_sc_rel'] > 20:
                        gain_of_function[mutation].append(restype_resnum_kinase)
                    elif not pd.isna(row['WT']) and pd.isna(row[mutation]):
                        loss_of_function[mutation].append(restype_resnum_kinase)

        except Exception as e:
            this_error = f"Error during mutation analysis: {e}"
            raise MAVISpMultipleError(warning=warnings,
                              critical=[MAVISpCriticalError(this_error)])

        # Convert the dictionaries into DataFrames
        try:
            gain_of_function_df = pd.DataFrame(
                [(k, ','.join(v)) for k, v in gain_of_function.items()],
                columns=['mutation', 'Phosphorylation - gain of function']).set_index('mutation')
            loss_of_function_df = pd.DataFrame(
                [(k, ','.join(v)) for k, v in loss_of_function.items()],
                columns=['mutation', 'Phosphorylation - loss of function']).set_index('mutation')
            final_table = pd.DataFrame({'mutation': mutations}).set_index('mutation')
            final_table = final_table.join(loss_of_function_df, how='left').join(gain_of_function_df, how='left')
            # Add 'PTM effect by mutation' based on gain or loss of function
            final_table['Mutation predicted to add new phosphorylation site'] = final_table['Phosphorylation - gain of function'].apply(lambda x: False if pd.isna(x) else True)
            self.data = final_table
        except Exception as e:
            this_error = f"Error compiling results: {e}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings, critical=[])

class TaccDenovoPhospho(DenovoPhospho):

    expected_files = ['aggregated_filtered_output.csv', 'acc_REL.csv']
    sasa_fname = 'acc_REL.csv'

    def _parse_sas(self, fname):

        sas_data = pd.read_csv(fname, usecols=['residue', 'acc_average'])
        sas_data.rename(columns={'residue': 'resn','acc_average': 'sas_sc_rel'}, inplace=True)
        sas_data.set_index('resn', inplace=True)
        return sas_data

class EnsembleDenovoPhospho(MavispMultiEnsembleModule, module_class=TaccDenovoPhospho):
    module_dir = "denovo_phospho"
    name = "denovo_phospho"

class PTMs(MavispModule):

    module_dir = "ptm"
    name = "ptms"
    expected_files = ['summary_stability.txt',
                      'sasa.rsa',
                      'metatable.csv']

    sasa_fname = expected_files[1]

    allowed_ptms = ['s', 'p', 'y']
    allowed_ptm_muts = {'S' : 's',
                        'T' : 'p',
                        'Y' : 'y'}
    protein_chain = 'A'
    slim_pattern = re.compile('.* \((CLV|DEG|DOC|LIG|MOD|TRG)_[A-Za-z0-9_-]+\), [0-9]+-[0-9]+, (ggetELM|ELM)')

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

    def _parse_sas(self, fname):

        try:
            return pd.read_fwf(fname,
                skiprows=4, skipfooter=4, header=None, widths=[4,4,1,4,9,6,7,6,7,6,7,6,7,6],
                names = ['entry', 'rest', 'chain', 'resn', 'all_abs', 'sas_all_rel', 'sas_sc_abs',
                'sas_sc_rel', 'sas_mc_abs', 'sas_mc_rel', 'sas_np_abs', 'sas_np_rel', 'sas_ap_abs',
                'sas_ap_rel'],
                usecols = ['resn', 'sas_sc_rel'],
                index_col = 'resn').fillna(pd.NA)

        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the sasa.rsa file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                        critical=[MAVISpCriticalError(this_error)])

    def ingest(self, mutations):

        warnings = []

        ptm_files = os.listdir(os.path.join(self.data_dir, self.module_dir))

        if not set(self.expected_files).issubset(ptm_files):
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
            binding_energies_available = True
        except FileNotFoundError as e:
            ddg_binding = pd.DataFrame(columns=['mutation', 'ddg_avg', 'ddg_std', 'ddg_min', 'ddg_max', 'idx'])
            warnings.append(MAVISpWarningError(f"summary_binding.txt not found - changes in free energy will not be used to classify function"))
            binding_energies_available = False
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the summary_binding.txt file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                        critical=[MAVISpCriticalError(this_error)])

        try:
            rsa_monomer = self._parse_sas(os.path.join(self.data_dir, self.module_dir, self.sasa_fname))
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the {self.sasa_fname} file. Arguments:{e.args}"
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
        cols = [ 'phosphorylation_site', 'site_in_slim', 'sas_sc_rel',
                 'stability_ddg_ptm', 'binding_ddg_ptm', 'binding_ddg_mut',
                 'regulation', 'stability', 'function', 'acc_rel' ]

        final_table = final_table[final_table.columns.intersection(cols)]

        if not binding_energies_available:
            final_table = final_table.drop(columns=['binding_ddg_mut', 'binding_ddg_ptm'])

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

class TaccPTMs(PTMs):

    expected_files = ['summary_stability.txt',
                      'acc_REL.csv',
                      'metatable.csv']

    sasa_fname = expected_files[1]

    def _parse_sas(self, fname):

        return pd.read_csv(fname).set_index('residue').rename(columns={'acc_average' : 'sas_sc_rel'})

    def ingest(self, mutations):

        super().ingest(mutations)

        self.data.rename(columns = {"PTM residue SASA (%)" : "PTM residue SASA (%), average",
                                    "acc_std" : "PTM residue SASA (%), standard deviation"})

class EnsemblePTMs(MavispMultiEnsembleModule, module_class=TaccPTMs):
    module_dir = "ptm"
    name = "ptms"

class CancermutsTable(MavispModule):
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
        required_columns = ['aa_position', 'ref_aa', 'alt_aa', 'gnomad_genome_af', 'gnomad_exome_af', 'REVEL_score', 'sources', 'genomic_mutation']
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
            this_error = f"Not all the annotated mutations were found in the Cancermuts table"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        # keep rows with mutations only
        cancermuts = cancermuts.loc[available_mutations]

        # check if all the mutations have a defined REVEL score
        # if np.any(pd.isna(cancermuts['REVEL_score'])):
        #     warnings.append(MAVISpWarningError("One or more mutations in the Cancermuts table don't have an associated REVEL score"))

        # filter by column and pretty rename column names
        cancermuts = cancermuts[['genomic_mutation', 'gnomad_genome_af', 'gnomad_exome_af', 'REVEL_score', 'sources']]

        self.data = cancermuts.rename(columns={ 'genomic_mutation' : 'HGVSg',
                                                'gnomad_genome_af' : 'gnomAD genome allele frequency',
                                                'gnomad_exome_af'  : 'gnomAD exome allele frequency',
                                                'REVEL_score'      : 'REVEL score',
                                                'sources'          : 'Mutation sources' })

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])

class ClinVar(MavispModule):

    module_dir = "clinvar"
    name = "clinvar"
    variants_fname = "variants_output.csv"
    mutation_pattern = re.compile('\(p\.([A-Z][a-z]{2}[0-9]+[A-Z][a-z]{2})\)')

    def _get_mutation_string(self, row):
        matches = self.mutation_pattern.findall(row['variant_name'])
        if len(matches) != 1:
            raise TypeError(f"processing ot variant {row['variant_id']} name failed")
        mutation = matches[0]
        try:
            mutation = three_to_one[mutation[0:3].upper()] + \
                       mutation[3:-3] + \
                       three_to_one[mutation[-3:].upper()]
        except KeyError:
            raise TypeError(f"non-standard residue name in variant {row['variant_id']}")

        return mutation

    def ingest(self, mutations):
        warnings = []

        clinvar_files = os.listdir(os.path.join(self.data_dir, self.module_dir))

        if self.variants_fname not in clinvar_files:
            this_error = f"Required file {self.variants_fname} was not found"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        try:
            clinvar_found = pd.read_csv(os.path.join(self.data_dir, self.module_dir, self.variants_fname), sep=';')
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the csv files. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        if 'interpretation' not in clinvar_found.columns:
            this_error = f"The variants_output.csv file must contain the interpretation column"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        if not ('clinvar_code' in clinvar_found.columns) ^ ('variant_id' in clinvar_found.columns):
            this_error = f"The variants_output.csv file must contain either the clinvar_code or the variant_id column"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        if 'clinvar_code' in clinvar_found.columns:
            warnings.append(MAVISpWarningError(f"the input file has the old style clinvar_code column"))
            clinvar_found.rename({'clinvar_code' : 'variant_id'})

        clinvar_found['variant_id'] = clinvar_found['variant_id'].astype(str)

        try:
            clinvar_found['mutations'] = clinvar_found.apply(self._get_mutation_string, axis=1)
        except TypeError as e:
            this_error = f"Error processing the variants file: {e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        if "number_of_stars" in clinvar_found.columns:
            clinvar_found['number_of_stars'] = clinvar_found['number_of_stars'].astype(str)
            clinvar_found = clinvar_found.groupby('mutations').agg(lambda x: ", ".join(list(x)))[['variant_id', 'interpretation', 'number_of_stars']]
            self.data = clinvar_found.rename({ 'variant_id'     : 'ClinVar Variation ID',
                                               'interpretation' : 'ClinVar Interpretation',
                                               'number_of_stars': 'ClinVar Review Status'}, axis=1)
        else:
            warnings.append(MAVISpWarningError(f"the variant_output.csv file doesn't contain the number_of_stars column (ClinVar review status)"))
            clinvar_found = clinvar_found.groupby('mutations').agg(lambda x: ", ".join(list(x)))[['variant_id', 'interpretation']]
            self.data = clinvar_found.rename({ 'variant_id'     : 'ClinVar Variation ID',
                                               'interpretation' : 'ClinVar Interpretation',}, axis=1)

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])

class EVE(MavispModule):

    module_dir = "eve"
    name = "eve"

    def __init__(self, data_dir=None):

        super().__init__(data_dir)

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
        required_columns = ['mutations', 'EVE_scores', 'EVE_classes_75_pct_retained']
        if not set(required_columns).issubset(eve.columns):
            this_error = f"input table doesn't have all the required columns"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        # process table
        eve = eve[required_columns]
        eve = eve.set_index('mutations')

        self.data = eve.rename(columns={ 'EVE_scores' : 'EVE score',
                                         'EVE_classes_75_pct_retained' : 'EVE classification (25% Uncertain)'})

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])

class AlphaFoldMetadata(MavispModule):

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

class DeMaSk(MavispModule):

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

class AlphaMissense(MavispModule):

    module_dir = "alphamissense"
    name = "alphamissense"

    def ingest(self, mutations):

        warnings = []

        am_files = os.listdir(os.path.join(self.data_dir, self.module_dir))
        if len(am_files) != 1:
            this_error = f"multiple or no files found in {am_files}; only one expected"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        am_file = am_files[0]

        log.info(f"parsing AlphaMissense data file {am_file}")

        try:

            afm = pd.read_csv(os.path.join(self.data_dir, self.module_dir, am_file),
                             usecols=['protein_variant', 'am_pathogenicity', 'am_class'],
                             dtype={ 'protein_variant' :  'string',
                                     'am_pathogenicity' : 'float32',
                                     'am_class' :         'category' },
                             sep='\t',
                             index_col='protein_variant')
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the csv files. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        self.data = afm.rename(columns = {'am_pathogenicity'     : 'AlphaMissense pathogenicity score',
                                          'am_class'   : 'AlphaMissense classification',})

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])

class GEMME(MavispModule):

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
        # be the one with score 'None')
        wts = gemme.groupby('res').apply(lambda x: x[pd.isna(x['score'])]['mut'].to_list()[0])
        wts.name = 'wt'

        # join WT definition on main dataframe
        gemme = gemme.join(wts, on='res')

        # reconstruct mutations in the usual format
        gemme['mutations'] = gemme['wt'] + gemme['res'] + gemme['mut']

        # drop unnecessary columns and rename for pretty
        gemme = gemme.drop(columns=['wt', 'res', 'mut'])

        # rank-normalize
        max_s = gemme['score'].max()
        min_s = gemme['score'].min()
        gemme['score_rn'] = (gemme['score'] - min_s) / (max_s - min_s)

        # rename for pretty
        gemme = gemme.rename(columns={'score'    : 'GEMME Score',
                                      'score_rn' : 'GEMME Score (rank-normalized)'})

        self.data = gemme.set_index('mutations')

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])

class FunctionalDynamics(MavispModule):

    module_dir = "functional_dynamics"
    name = "functional_dynamics"
    allowed_values = set(['uncertain', 'neutral', 'destabilizing', 'stabilizing'])

    def _parse_table(self, fname):

        df = pd.read_csv(fname, sep='\t', index_col=False)

        if len(df.columns) != 2:
            this_error = f"Functional Dynamics table must have two columns (mutations and classification)"
            raise TypeError(this_error)

        if 'mutation' not in df.columns:
            this_error = f"Functional Dynamics table must have a column named mutation"
            raise TypeError(this_error)

        df = df.set_index('mutation')

        other_col = df.columns.to_list()[0]

        if not set(df[other_col]).issubset(self.allowed_values):
            this_error = f"Functional Dynamics table must have only the following values: {', '.join(self.allowed_values)}"
            raise TypeError(this_error)

        df = df.rename(columns={other_col : f"Functional dynamics ({other_col})"})

        return df

    def ingest(self, mutations):

        warnings = []

        self.data = pd.DataFrame({'mutations' : mutations}).set_index('mutations')

        # check that we have one or more directories
        fd_dirs = os.listdir(os.path.join(self.data_dir, self.module_dir))

        if len(fd_dirs) == 0 or not all([os.path.isdir(os.path.join(self.data_dir, self.module_dir, d)) for d in fd_dirs]):
            this_error = f"functional_dynamics directory must contain one or more subdirectory"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])
        for fd_dir in fd_dirs:

            # check that we have only one file
            fd_files = os.listdir(os.path.join(self.data_dir, self.module_dir, fd_dir))
            if len(fd_files) != 1:
                this_error = f"multiple or no files found in {fd_files}; only one expected"
                raise MAVISpMultipleError(warning=warnings,
                                          critical=[MAVISpCriticalError(this_error)])

            fd_file = fd_files[0]

            log.info(f"parsing Functional Dynamics data file {fd_file}")

            # parse functional dynamics table
            try:
                fd = self._parse_table(os.path.join(self.data_dir, self.module_dir, fd_dir, fd_file))
            except Exception as e:
                this_error = f"Exception {type(e).__name__} occurred when parsing the csv files. Arguments:{e.args}"
                raise MAVISpMultipleError(warning=warnings,
                                          critical=[MAVISpCriticalError(this_error)])

            # join this functional dynamics table to main table
            self.data = self.data.join(fd)

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])

class EnsembleFunctionalDynamics(MavispMultiEnsembleModule, module_class=FunctionalDynamics):
    module_dir = "functional_dynamics"
    name = "functional_dynamics"

class FunctionalSites(MavispModule):

    module_dir = "functional_sites"
    name = "functional_sites"
    accepted_filenames = ['cofactor_local_aggregate.txt',
                            'active_site_local_aggregate.txt']
    column_headers = {'cofactor_local_aggregate.txt'    : 'Functional sites (cofactor)',
                      'active_site_local_aggregate.txt' : 'Functional sites (active site)'}
    _site_re = re.compile('[ACDEFGHIKLMNPQRSTVWY][0-9]+')

    def _parse_table(self, fname):

        df = pd.read_csv(fname, index_col=False, header=None)

        if len(df.columns) != 1:
            this_error = f"Functional Sites table {fname} must have one column (sites)"
            raise TypeError(this_error)

        df[0] = df[0].str.strip()

        # check if each element of the column is in the format Xn, where X is a
        # residue type and N is an arbitrary number
        if df[0].apply(lambda x: self._site_re.match(x) is None).any():
            this_error = f"Functional Sites table {fname} must have sites in the format Xn, where X is a residue type and N is an arbitrary number"
            raise TypeError(this_error)

        df = df.drop_duplicates()
        df['sites_type_table'] = df[0].str[0]
        df['sites'] = pd.Series(df[0].str[1:].astype(int))
        df = df.drop(columns=[0])

        if df['sites'].duplicated().any():
            this_error = f"Functional Sites table {fname} has duplicated sites with different residue types"
            raise TypeError(this_error)

        df = df.set_index('sites')
        df['classification'] = pd.Series(index=df.index, data='damaging')

        return df

    def ingest(self, mutations):

        warnings = []

        fs_files = os.listdir(os.path.join(self.data_dir, self.module_dir))

        if not set(fs_files).issubset(set(self.accepted_filenames)):
            this_error = f"the input files for Functional Sites must be named {', '.join(self.accepted_filenames)}"
            raise MAVISpMultipleError(warning=warnings,
                                        critical=[MAVISpCriticalError(this_error)])

        out_df = pd.DataFrame({'mutations' : mutations})
        out_df['sites'] = out_df['mutations'].str[1:-1].astype(int)
        out_df['sites_type_mut'] = out_df['mutations'].str[0]
        out_df = out_df.set_index('mutations')

        for fs_file in fs_files:

            log.info(f"parsing Functional Sites data file {fs_file}")

            try:
                df = self._parse_table(os.path.join(self.data_dir, self.module_dir, fs_file))
            except Exception as e:
                this_error = f"Exception {type(e).__name__} occurred when parsing the csv files. Arguments:{e.args}"
                raise MAVISpMultipleError(warning=warnings,
                                            critical=[MAVISpCriticalError(this_error)])
            out_df = out_df.join(df, on='sites')

            if not out_df.apply(lambda r: r['sites_type_table'] == r['sites_type_mut'] if not pd.isna(r['classification']) else True, axis=1).all():
                this_error = f"Functional Sites table {fs_file} has sites with residue type different than in mutations"
                raise MAVISpMultipleError(warning=warnings,
                                            critical=[MAVISpCriticalError(this_error)])

            out_df = out_df.rename(columns={'classification' : self.column_headers[fs_file]})
            out_df = out_df.drop(columns=['sites_type_table'])

        self.data = out_df.drop(columns=['sites_type_mut', 'sites']).fillna('neutral')

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])

class AllosigmaPSNLongRange(MavispModule):

    module_dir = "long_range"
    name = "long_range"
    method_dir = "allosigma2_psn"
    exp_files = ['allosigma_mut.txt', 'results_summary.txt']

    def _generate_allosigma_psn_classification(self, row, ensemble_data):
        res_num = row['res_num']
        allosigma_mode = row['allosigma-mode']
        variant_sites = ensemble_data[ensemble_data['Variant_Site'].astype(str) == res_num]

        # Mutation not classified by Allosigma
        if not allosigma_mode in ['UP', 'DOWN']:
            return 'uncertain'

        # No predicted Allosigma effects
        if variant_sites.empty:
            return 'neutral'

        # Predicted Allosigma effects
        total_paths = variant_sites['Total_Paths'].max()
        if total_paths >= 1:
            # Predictions validated
            return 'damaging'
        else:
            # Predictions could not be validated
            return 'uncertain'

    def _read_file(self, file_path, columns, warnings):
        try:
            # Read the file with the specified columns
            df = pd.read_csv(file_path, sep='\t', usecols=columns)
        except Exception as e:
            this_error = f"Error while parsing file {file_path}: {e}"
            raise MAVISpMultipleError(warning=warnings,
                                            critical=[MAVISpCriticalError(this_error)])
        return df

    def ingest(self, mutations):

        warnings = []

        base_path = os.path.join(self.data_dir, self.module_dir, self.method_dir)
        psn_files = os.listdir(base_path)

        if not set(psn_files).issubset(set(self.exp_files)):
            this_error = (f"the input files for AlloSigma2-PSN must be named {', '.join(self.exp_files)}, "
            f"found {', '.join(psn_files)}")
            raise MAVISpMultipleError(warning=warnings,
                                        critical=[MAVISpCriticalError(this_error)])

        # Parsing input files + required columns
        df_simple_data = self._read_file(
            file_path=os.path.join(base_path, self.exp_files[0]),
            columns=['wt_residue', 'position', 'mutated_residue', 'allosigma-mode'],
            warnings=warnings)
        df_ensemble_data = self._read_file(
            file_path=os.path.join(base_path, self.exp_files[1]),
            columns=['Variant_Site', 'Total_Paths'],
            warnings=warnings)

        # Build mutations column + order columns
        df_simple_data['mutations'] = (df_simple_data['wt_residue'] +
            df_simple_data['position'].astype(str) +
            df_simple_data['mutated_residue'])
        df_simple_data = df_simple_data[['mutations', 'allosigma-mode']]

        #Define working copy of data
        psn = pd.DataFrame({'mutations' : mutations}).set_index('mutations')

        # Add residue number column
        psn['res_num'] = psn.index.str[1:-1].astype(str)

        # Merge allosigma data with psn
        psn = psn.merge(df_simple_data, on='mutations', how='left')

        try:
            # Perform classification
            psn['AlloSigma2-PSN classification'] = psn.apply(self._generate_allosigma_psn_classification, ensemble_data=df_ensemble_data, axis=1)
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred while performing allosigma-psn classifcation."
            raise MAVISpMultipleError(warning=warnings,
                                            critical=[MAVISpCriticalError(this_error)])
        # Drop unnecessary columns
        psn = psn.drop(columns=['res_num', 'allosigma-mode'])

        # Add new Allosigma-PSN column to data
        self.data = psn.set_index('mutations')

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])

class EnsembleAllosigmaPSNLongRange(MavispMultiEnsembleModule, module_class=AllosigmaPSNLongRange):
    module_dir = "long_range"
    name = "long_range"

class ExperimentalData(MavispModule):

    module_dir = "experimental_data"
    name = "experimental_data"

    hgvsp_regexp = "p\.[A-Z][a-z][a-z][0-9]+[A-Z][a-z][a-z]"

    def _hgvs_to_mavisp(self, hgvsp, offset):
        return f"{three_to_one_hgvsp[hgvsp[2:5]]}{int(hgvsp[5:-3]) + int(offset)}{three_to_one_hgvsp[hgvsp[-3:]]}"

    def _get_classification(self, series, thresholds, threshold_type):
        masks = []
        mask_descriptions = []

        # create a mask for every threshold we have
        for desc, thres in thresholds.items():
            if threshold_type == 'values':
                if isinstance(thres, numbers.Number):
                    mask = series == thres
                    masks.append(mask)
                    mask_descriptions.append(desc)
                else:
                    for t in thres:
                        mask = series == t
                        masks.append(mask)
                        mask_descriptions.append(desc)
            else:
                if len(thres) != 2 or not thres[0] < thres[1]:
                    raise RuntimeError("When using threshold type ranges, classes need to be made of a list of two values (min, max)")
                mask = (series >= thres[0]) & (series < thres[1])
                masks.append(mask)
                mask_descriptions.append(desc)

        all_masks = pd.concat(masks, axis=1)

        # check if any threshold overlap
        if any(all_masks.sum(axis=1) > 1):
            raise RuntimeError("One or more mutations belong to multiple classes; are your definitions overlapping?")

        # generate classification
        out_series = series.copy()

        for mask, desc in zip(masks, mask_descriptions):
            out_series[mask] = desc

        return out_series

    def ingest(self, mutations):
        warnings = []

        module_path = os.path.join(self.data_dir, self.module_dir)

        yaml_files = os.listdir(module_path)
        yaml_files = [ f for f in yaml_files if f.endswith('.yaml') ]

        if len(yaml_files) == 0:
            this_error = f"At least one metadata yaml file is required"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        all_data = pd.DataFrame({'mutations' : mutations}).set_index('mutations')

        for yaml_file in yaml_files:
            try:
                with open(os.path.join(module_path, yaml_file)) as fh:
                    metadata = yaml.safe_load(fh)
            except Exception as e:
                this_error = f"Error parsing YAML metadata {yaml_file}: {e}"
                raise MAVISpMultipleError(warning=warnings,
                                        critical=[MAVISpCriticalError(this_error)])

            for col in metadata['columns'].keys():

                col_metadata = metadata['columns'][col]

                if col_metadata['threshold_type'] not in ['values', 'ranges']:
                    this_error = f"Error in {yaml_file}, {col}: threshold_type can either be values or ranges"
                    raise MAVISpMultipleError(warning=warnings,
                                            critical=[MAVISpCriticalError(this_error)])


                try:
                    data = pd.read_csv(os.path.join(module_path, col_metadata['data_file']))
                except Exception as e:
                    this_error = f"Error parsing data file {col_metadata['data_file']}: {e}"
                    raise MAVISpMultipleError(warning=warnings,
                                            critical=[MAVISpCriticalError(this_error)])

                if col_metadata['mutation_format'] == 'hgvsp':
                    data = data[ ~ data[col_metadata['mutation_column']].str.endswith('=') ]
                    data = data[ ~ data[col_metadata['mutation_column']].str.endswith('Ter') ]
                    full_data_len = data.shape[0]
                    data = data[   data[col_metadata['mutation_column']].str.contains(self.hgvsp_regexp, regex=True, na=False) ]
                    if data.shape[0] != full_data_len:
                        warnings.append(MAVISpWarningError("rows with inconsistent HGVSp notation in mutation column were removed from the dataset"))

                    data['mutations'] = data[col_metadata['mutation_column']].apply(self._hgvs_to_mavisp, offset=col_metadata['offset'])

                elif col_metadata['mutation_format'] == 'mavisp':
                    data['mutations'] = data[col_metadata['mutation_column']]

                else:
                    this_error = f"in file {yaml_file}, mutations format should be either hgvsp or mavisp"
                    raise MAVISpMultipleError(warning=warnings,
                                              critical=[MAVISpCriticalError(this_error)])


                data = data[['mutations', col]].set_index('mutations')
                try:
                    data[f"{col} classification"] = self._get_classification(data[col], col_metadata['thresholds'], col_metadata['threshold_type'])
                except Exception as e:
                    this_error = f"Error while generating classification for {col} in {yaml_file}: {e}"
                    raise MAVISpMultipleError(warning=warnings,
                                              critical=[MAVISpCriticalError(this_error)])

                data = data.rename(columns={col                     : f"Experimental data ({metadata['assay']}, {col_metadata['header']})",
                                            f"{col} classification" : f"Experimental data classification ({metadata['assay']}, {col_metadata['header']})"})

                all_data = all_data.join(data)

        self.data = all_data

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])
