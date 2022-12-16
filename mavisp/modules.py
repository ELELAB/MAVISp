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
from mavisp.utils import three_to_one
import logging as log
import re

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
    allowed_ptms = ['s', 'p', 'y']
    allowed_ptm_muts = {'S' : 's',
                        'T' : 'p',
                        'Y' : 'y'}

    def _assign_regulation_class(self, row):
        # case S to T or T to S
        if set([row['alt'], row['ref']]) == set(['S', 'T']):
            return 'neutral'
        # cases:
            # any mutation sas < 20% (excluding neutral cases)
            # any T/S to Y
            # any Y to T/S
        elif row['sas_sc_rel'] < 20.0 or row['alt'] == 'Y' or (row['ref'] == 'Y' and (row['alt'] in ['S', 'T'])):
            return 'unkown'

        # cases:
            # any mutation that is not T to S, S to T, T/S to Y or Y to T/S and SASA >=20%:
        elif row['sas_sc_rel'] >= 20.0:
            return 'damaging'

        return '???'

    def _assign_ddg_stability_class(self, row, ddg_col_name):

        stab_co = 3.0
        neut_co = 2.0

        if row[ddg_col_name] > stab_co:
            return 'Destabilizing'
        elif row[ddg_col_name] < (- stab_co):
            return 'Stabilizing'
        elif (- neut_co) <= row[ddg_col_name] <= neut_co:
            return 'Neutral'
        else:
            return 'Uncertain'

    def _assign_stability_class(self, row, mut_col_name, ptm_col_name):

        if pd.isna(row[ptm_col_name]) or pd.isna(row[mut_col_name]):
            return pd.NA

        if row[ptm_col_name] == 'Uncertain' or row[mut_col_name] == 'Uncertain':
            return 'Uncertain'
        elif row[mut_col_name] == row[ptm_col_name]:
            return 'Neutral'
        else:
            return 'Damaging'

    def _assign_function_class(self, row):
        if row['site_in_slim']:
            return 'potentially_damaging'
        else:
            return 'unknown'

    def ingest(self, mutations):
        warnings = []

        expected_files = ['summary_stability.txt', 'sasa.rsa', 'metatable.csv']

        ptm_files = os.listdir(os.path.join(self.data_dir, self.module_dir))

        if not set(expected_files).issubset(ptm_files):
            this_error = f"required summary_stability.txt, sasa.rsa or metatable.csv file not found in {self.module_dir}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        try:
            ptms = pd.read_csv(os.path.join(self.data_dir, self.module_dir, 'summary_stability.txt'),
                delim_whitespace=True,
                header=None,
                names=['mutation', 'ddg_avg', 'ddg_std', 'ddg_min', 'ddg_max', 'idx'])
        except Exception as e:
            this_error = f"Exception {type(e).__name__} occurred when parsing the summary.txt file. Arguments:{e.args}"
            raise MAVISpMultipleError(warning=warnings,
                                        critical=[MAVISpCriticalError(this_error)])

        try:
            rsa = pd.read_fwf(os.path.join(self.data_dir, self.module_dir, 'sasa.rsa'),
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
            this_error = f"input table doesn't have all the required columns"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])

        # process table
        in_slims = cancermuts[['aa_position', 'linear_motif']].drop_duplicates().set_index('aa_position')
        in_slims['site_in_slim'] = ~ pd.isna(in_slims['linear_motif'])

        ptms['ref'] = ptms['mutation'].str[0]
        ptms['alt'] = ptms['mutation'].str[-1]
        ptms['number'] = ptms['mutation'].str[2:-1].astype(int)
        ptms['mutation'] = ptms['ref'] + ptms['number'].astype(str) + ptms['alt']

        ptms = ptms.set_index('number')
        ptms = ptms.join(rsa, on='number')

        ptm_muts = ptms[ptms.apply(lambda r: r['alt'].islower(), axis=1)]
        cancer_muts = ptms[~ ptms.apply(lambda r: r['alt'].islower(), axis=1)]

        if not set(ptm_muts['alt']).issubset(set(self.allowed_ptms)):
            warnings.append(MAVISpWarningError(f"in summary.txt, some of the mutations were not to {', '.join(self.allowed_ptms)}. These will be filtered out."))

            cancer_muts = cancer_muts[cancer_muts.apply(lambda r: r['alt'] in self.allowed_ptms)]

        if not set(ptm_muts.index).issubset(set(cancer_muts.index)):
            this_error = f"in summary.txt, some cancer mutations had no corresponding PTM mutation"
            raise MAVISpMultipleError(warning=warnings,
                                        critical=[MAVISpCriticalError(this_error)])

        if not np.all(ptm_muts.apply(lambda r: self.allowed_ptm_muts[r['ref']] == r['alt'], axis=1)):
            this_error = f"in summary.txt, the wild type residue and the wild type residue upon PTM didn't correspond"
            raise MAVISpMultipleError(warning=warnings,
                                        critical=[MAVISpCriticalError(this_error)])

        cancer_muts = cancer_muts.join(in_slims)

        cancer_muts = cancer_muts.rename(columns={'ddg_avg' : 'ddg_mut'})
        cancer_muts = cancer_muts.join(ptm_muts[['ddg_avg']]).rename(columns={'ddg_avg' : 'ddg_ptm'})

        cancer_muts = cancer_muts.set_index('mutation')

        cancer_muts['regulation'] = cancer_muts.apply(self._assign_regulation_class, axis=1)

        cancer_muts['cancer_mut_stab_class'] = cancer_muts.apply(self._assign_ddg_stability_class, ddg_col_name='ddg_mut', axis=1)
        cancer_muts['ptm_stab_class'] = cancer_muts.apply(self._assign_ddg_stability_class, ddg_col_name='ddg_ptm', axis=1)
        cancer_muts['stability'] = cancer_muts.apply(self._assign_stability_class, mut_col_name='cancer_mut_stab_class', ptm_col_name='ptm_stab_class', axis=1)
        cancer_muts['function'] = cancer_muts.apply(self._assign_function_class, axis=1)

        cancer_muts['PTMs'] = 'P'
        pd.set_option('display.max_columns', None)

        cancer_muts = cancer_muts[['PTMs', 'sas_sc_rel', 'ddg_ptm', 'regulation', 'stability', 'function']]

        self.data = cancer_muts.rename(columns={'sas_sc_rel'      : "PTM residue SASA (%)" ,
                                                'ddg_ptm'         : "Change in stability with PTM (FoldX5, kcal/mol)",
                                                'regulation'      : "PTM effect in regulation",
                                                'stability'       : "PTM effect in stability" ,
                                                'function'        : "PTM effect in function",
                                                'site in slim'    : "is PTM site part of SLiM",
                                                'effect_function' : "PTM effect in function"})

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
            this_error = f"One or two files expected for Clinvar, they must be {self.found_fname} and/or {self.missing_fname}"
            raise MAVISpMultipleError(warning=warnings,
                                      critical=[MAVISpCriticalError(this_error)])


        if self.found_fname not in clinvar_files:
            this_error = f"variants_output.csv expected for Clinvar"
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
            warnings.append(MAVISpWarningError(f"file {self.missing} not found in Clinvar module"))
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
        clinvar_found = clinvar_found.groupby('variant_name').agg(lambda x: ", ".join(list(x)))[[id_col, 'interpretation']]
        self.data = clinvar_found.rename({ id_col          : 'Clinvar Variation ID',
                                          'interpretation' : 'ClinVar interpretation'}, axis=1)

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
