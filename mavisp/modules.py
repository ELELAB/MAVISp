# MAVISp - classes for handling different MAVISp modules
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
import pandas as pd
from mavisp.methods import *
import logging as log

class DataType(object):
    def __init__(self, data_dir=None, stop_at='critical'):

        self.data_dir = data_dir
        self.data = None

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

    def get_dataset_view(self):
        return self.data

class MultiMethodDataType(DataType):
    def __init__(self, data_dir=None):

        super().__init__(data_dir)

    def ingest(self, mutations):

        self.data = pd.DataFrame({'mutations' : mutations}).set_index('mutations')

        warnings = []

        method_dirs = os.listdir(os.path.join(self.data_dir, self.module_dir))
        
        if not set(method_dirs).issubset(set(self.methods.keys())):
            this_error = f"One or more {self.name} methods are not supported"
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[MAVISpCriticalError(this_error)])
        
        for method_dir in method_dirs:
            self.methods[method_dir].parse(os.path.join(self.data_dir, self.module_dir, method_dir))
            self.data = self.data.join(self.methods[method_dir].data)

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[])

class Stability(MultiMethodDataType):

    module_dir = "stability"
    name = "stability"
    methods = {'foldx5'                      : MutateXStability(version="FoldX5"),
               'rosetta_cartddg2020_ref2015' : RosettaDDGPredictionStability(version='Rosetta Flexddg2020'),
               'rosetta_ref2015'             : RosettaDDGPredictionStability(version='Rosetta Flexddg')}

    def ingest(self, mutations):

        self.data = pd.DataFrame({'mutations' : mutations}).set_index('mutations')

        warnings = []

        this_error = "Stability folder has to contain only 1 dataset"
        tmp = os.listdir(os.path.join(self.data_dir, self.module_dir))
        if len(tmp) != 1:
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[MAVISpCriticalError(this_error)])

        structure_ID, residue_range = tmp[0].split("_", maxsplit=1)

        tmp = os.listdir(os.path.join(self.data_dir, self.module_dir, f'{structure_ID}_{residue_range}'))
        if len(tmp) != 1:
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[MAVISpCriticalError(this_error)])
        method = tmp[0]

        tmp = os.listdir(os.path.join(self.data_dir, self.module_dir, f'{structure_ID}_{residue_range}', method))
        if len(tmp) != 1:
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[MAVISpCriticalError(this_error)])
        model = tmp[0]

        method_dirs = os.listdir(os.path.join(self.data_dir, self.module_dir, f'{structure_ID}_{residue_range}', method, model))

        if not set(method_dirs).issubset(set(self.methods.keys())):
            this_error = f"One or more {self.name} methods are not supported"
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[MAVISpCriticalError(this_error)])

        for method_dir in method_dirs:
            self.methods[method_dir].parse(os.path.join(self.data_dir, self.module_dir, f'{structure_ID}_{residue_range}', method, model, method_dir))
            self.data = self.data.join(self.methods[method_dir].data)

        keys = self.data.keys().to_list()
        if len(keys) != 2 or not ( 'Rosetta' in keys[0] and 'FoldX' in keys[1] or 'Rosetta' in keys[1] and 'FoldX' in keys[0]):
            warnings.append(MAVISpWarningError("Stability classification can only be calculated if exactly one Rosetta and one MutateX datasets are available"))
        
        self.data['Stability classification'] = self.data.apply(self._generate_stability_classification, axis=1)

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[])

    def _generate_stability_classification(self, row):

        keys = [ k for k in row.keys() if k.startswith('Stability') ]

        if len(keys) == 2:
            if   'Rosetta' in keys[0] and 'FoldX' in keys[1]:
                rosetta_header, foldx_header    = keys
            elif 'Rosetta' in keys[1] and 'FoldX' in keys[0]:
                foldx_header,   rosetta_header  = keys
            else:
                return pd.NA
        else:
            return pd.NA

        stab_co = 3.0
        neut_co = 2.0

        if rosetta_header not in row.index or foldx_header not in row.index:
            return pd.NA

        if row[foldx_header] > stab_co and row[rosetta_header] > stab_co:
            return 'Destabilizing'
        if row[foldx_header] < (- stab_co) and row[rosetta_header] < (- stab_co):
            return 'Stabilizing'
        if (- neut_co) < row[foldx_header] < neut_co and (- neut_co) < row[rosetta_header] < neut_co:
            return 'Neutral'
        return 'Uncertain'


class LocalInteractions(MultiMethodDataType):
    module_dir = "local_interactions"
    name = "local_interactions"
    methods = {'foldx5'                      : MutateXBinding(version="FoldX5"),
                'rosetta_flexddg_talaris2014' : RosettaDDGPredictionBinding(version='Rosetta Talaris 2014')}

    def ingest(self, mutations):

        warnings = []

        try:
            super().ingest(mutations)
        except MAVISpMultipleError as e:
            if len(e.critical) > 0:
                raise
        else:
            e = None
        keys = self.data.keys().to_list()
        if len(keys) != 2 or not ( 'Rosetta' in keys[0] and 'FoldX' in keys[1] or 'Rosetta' in keys[1] and 'FoldX' in keys[0]):
            warnings.append(MAVISpWarningError("Stability classification can only be calculated if exactly one Rosetta and one MutateX datasets are available"))

        self.data['Local Int. classification'] = self.data.apply(self._generate_local_interactions_classification, axis=1)

        if e is None and len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[])
        elif len(warnings) > 0:
            e.warning.extend(warnings)
            raise e

    def _generate_local_interactions_classification(self, row):

        keys = [ k for k in row.keys() if k.startswith('Local Int.') ]

        if len(keys) == 2:
            if   'Rosetta' in keys[0] and 'FoldX' in keys[1]:
                rosetta_header, foldx_header    = keys
            elif 'Rosetta' in keys[1] and 'FoldX' in keys[0]:
                foldx_header,   rosetta_header  = keys
            else:
                return pd.NA
        else:
            return pd.NA

        stab_co =  1.0

        if rosetta_header not in row.index or foldx_header not in row.index:
            return pd.NA
        if row[foldx_header] > stab_co and row[rosetta_header] > stab_co:
            return 'Destabilizing'
        if row[foldx_header] < (- stab_co) and row[rosetta_header] < (- stab_co):
            return 'Stabilizing'
        if (- stab_co) <= row[foldx_header] <= stab_co and (- stab_co) <= row[rosetta_header] <= stab_co:
            return 'Neutral'
        return 'Uncertain'

class LongRange(MultiMethodDataType):

    module_dir = "long_range"
    name = "long_range"
    methods = {'allosigma2' : AlloSigma(version=2)}

class References(DataType):

    module_dir = "pmid_list"
    name = "references"

    def ingest(self, mutations):
        warnings = []

        pmid_files = os.listdir(os.path.join(self.data_dir, self.module_dir))
        if len(pmid_files) != 1:
            this_error = f"multiple or no files found in {pmid_files}; only one expected"
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[MAVISpCriticalError(this_error)])

        pmid_file = pmid_files[0]

        log.info(f"parsing PMID file {pmid_file}")

        try:
            pmid = pd.read_csv(os.path.join(self.data_dir, self.module_dir, pmid_file), delim_whitespace=True)
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the csv files. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[MAVISpCriticalError(this_error)])

        if 'mutation' not in pmid.columns or 'PMID' not in pmid.columns:
            this_error = f"The CSV file must contain the following columns: mutation, PMID"
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[MAVISpCriticalError(this_error)])

        if len(pmid.columns) != 2:
            warnings.append(MAVISpWarningError("The CSV file has more than two columns"))

        pmid = pmid[['mutation', 'PMID']]
        pmid = pmid.set_index('mutation')
        pmid = pmid[~ pmid.index.duplicated(keep='first')]
        self.data = pmid.rename(columns={'PMID':'PMID / DOI'})

        if len(warnings) > 0:
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[])

class PTMs(DataType):

    module_dir = "ptm"
    name = "ptms"

    def ingest(self, mutations):
        warnings = []

        ptm_files = os.listdir(os.path.join(self.data_dir, self.module_dir))
        if len(ptm_files) == 0:
            this_error = f"no files found in {self.module_dir}; one or more expected"
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=[MAVISpCriticalError(this_error)])

        ptm_data = []
        for fname in ptm_files:
            try:
                ptms = pd.read_csv(os.path.join(self.data_dir, self.module_dir, fname), delim_whitespace=True)
            except Exception as e:
                this_error = f"Exception {type(e).__name__} occurred when parsing the csv files. Arguments:{e.args}"
                raise MAVISpMultipleError(warning=warnings, 
                                          critical=[MAVISpCriticalError(this_error)])

            ptms['mutation'] = os.path.splitext(os.path.basename(fname))[0]

            if len(ptms) != 1:
                tmp_data = {}
                for c in ptms.columns:
                    tmp_data[c] = ", ".join(map(str, ptms[c].to_list()))
                ptms = pd.DataFrame(tmp_data)

            ptms = ptms.set_index('mutation')
            ptms = ptms[~ ptms.index.duplicated(keep='first')]
            ptm_data.append(ptms)

        ptm_data = pd.concat(ptm_data)
        self.data = ptm_data.rename(columns={   '#ptm'                 : "PTMs",
                                                'SASA(%)'              : "PTM residue SASA (%)" ,
                                                'ddG(foldX5,kcal/mol)' : "Change in stability with PTM (FoldX5, kcal/mol)",
                                                'effect_regulation'    : "PTM effect in regulation",
                                                'effect_stability'     : "PTM effect in stability" ,
                                                'effect_function'      : "PTM effect in function",
                                                'notes'                : "PTM notes"})

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
