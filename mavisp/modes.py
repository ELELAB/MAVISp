# MAVISp - classes for handling different MAVISp modules
# Copyright (C) 2022 Matteo Tiberti, Danish Cancer Society
#           (C) 2023 Jérémy Vinhas, Danish Cancer Society
#           (C) 2024 Pablo Sánchez-Izquierdo, Danish Cancer Society
#           (C) 2024 Eleni Kiahaki, Danish Cancer Society
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

from mavisp.modules import *
import yaml

class MAVISpMode:

    name = None
    supported_modules = None
    supported_metadata = None

    def _parse_metadata_file(self, data_dir, system, mode):
        log.info("parsing metadata file")

        metadata_path = os.path.join(data_dir, system, self.name, 'metadata.yaml')

        with open(metadata_path) as fh:
            return yaml.safe_load(fh)

class MAVISpSimpleMode(MAVISpMode):

    name = 'simple_mode'
    supported_modules = [ CancermutsTable,
                          Pfam,
                          TED,
                          PTMs,
                          DenovoPhospho,
                          LongRange,
                          SimpleStability,
                          LocalInteractions,
                          SAS,
                          LocalInteractionsDNA,
                          LocalInteractionsHomodimer,
                          FunctionalSites,
                          ClinVar,
                          AlphaFoldMetadata,
                          DeMaSk,
                          GEMME,
                          EVE,
                          AlphaMissense,
                          EFoldMine,
                          ExperimentalData ]
    module_order = ['cancermuts', 'pfam', 'ted', 'stability', 'efoldmine', 'local_interactions',
    'local_interactions_DNA', 'local_interactions_homodimers', 'sas', 'ptms',
    'denovo_phospho', 'long_range', 'functional_sites', 'clinvar', 'alphafold',
    'demask', 'gemme', 'eve', 'alphamissense', 'experimental_data']
    supported_metadata = ['uniprot_ac', 'refseq_id', 'review_status', 'curators', 'gitbook_entry',
                          'allosigma_distance_cutoff', 'allosigma_distance_mode', 
                          'structure_source', 'structure_description', 'linker_design', 'pdb_id', 'cofactors_in_structure']
    index_cols = ['system', 'uniprot_ac', 'refseq_id', 'review_status', 'curators', 'gitbook_entry',
                  'allosigma_distance_cutoff', 'allosigma_distance_mode', 'structure_source',
                  'structure_description', 'linker_design', 'pdb_id', 'cofactors_in_structure']
    index_col_labels = {'system': "Protein",
                        'uniprot_ac': 'Uniprot AC',
                        'refseq_id': "RefSeq ID",
                        'review_status': 'Review status',
                        'curators': 'Curators',
                        'gitbook_entry': 'GitBook report',
                        'allosigma_distance_cutoff': 'Distance cut-off used for AlloSigma2',
                        'allosigma_distance_mode' : 'Contact calculation mode for AlloSigma2 filtering',
                        'structure_source': 'Structure source',
                        'structure_description': 'Description of structure source',
                        'linker_design': 'Linker design included',
                        'pdb_id': 'PDB ID',
                        'cofactors_in_structure' : 'Cofactors in starting structure'}
    
    structure_sources = {
        "AFDB": "AlphaFold database",
        "AF3": "AlphaFold3 webserver",
        "AF2": "AlphaFold2 standalone",
        "PDB": "Experimental PDB structure",
        "Mod": "Homology model (PDB template, reconstruction)"}

    allosigma_modes = { "CA-CA", "atomic_contacts"}

    supported_cofactors = {"Zn2+", "Mg2+", "ADP", " ATP", "GDP", "GTP", "NADH", "NAD+", "FADH", "FAD+", "Ca2+", "Mn2+", "Fe2+", "Fe3+"}

    def parse_metadata(self, data_dir, system):
        out_metadata = {k: None for k in self.supported_metadata}
        mavisp_criticals = []
        structure_source = None
        metadata = self._parse_metadata_file(data_dir, system, self.name)

        try:
            curators = ', '.join(
                [f"{curator} ({', '.join(metadata['curators'][curator]['affiliation'])})" for curator in metadata['curators'].keys()]
            )
        except KeyError:
            log.debug("There is no curators field in metadata file")
            curators = None
            mavisp_criticals.append(MAVISpCriticalError("curators field not found in metadata file"))
        out_metadata['curators'] = curators

        for k in ['uniprot_ac', 'refseq_id', 'review_status']:
            try:
                out_metadata[k] = str(metadata[k])
            except KeyError:
                log.debug(f"There is no {k} field in metadata file")
                out_metadata[k] = None
                mavisp_criticals.append(MAVISpCriticalError(f"{k} was not found in the metadata file"))

        if 'allosigma_distance_cutoff' in metadata.keys():
            out_metadata['allosigma_distance_cutoff'] = ', '.join(map(str, metadata['allosigma_distance_cutoff']))
        else:
            out_metadata['allosigma_distance_cutoff'] = ''

        if 'allosigma_distance_mode' in metadata.keys():

            if not 'allosigma_distance_cutoff' in metadata.keys():
                mavisp_criticals.append(MAVISpCriticalError(f"in metadata, allosigma_distance_cutoff must be used when allosigma_distance_mode is used"))
            elif not len(metadata['allosigma_distance_cutoff']) == len(metadata['allosigma_distance_mode']):
                mavisp_criticals.append(MAVISpCriticalError(f"in metadata, allosigma_distance_mode must be have the same number of entries as allosigma_distance_cutoff"))

            if not set(metadata['allosigma_distance_mode']).issubset(self.allosigma_modes):
                mavisp_criticals.append(MAVISpCriticalError(f"only these modes are allowed for allosigma_distance_mode: {', '.join(list(self.allosigma_modes))}"))

            out_metadata['allosigma_distance_mode'] = ', '.join(map(str, metadata['allosigma_distance_mode']))
        else:
            out_metadata['allosigma_distance_mode'] = ''

        if 'gitbook_entry' in metadata.keys():
            out_metadata['gitbook_entry'] = metadata['gitbook_entry']
        else:
            out_metadata['gitbook_entry'] = ''

        try:
            structure_source = metadata["structure_source"]
        except KeyError as e:
            log.debug(f"Missing optional field in metadata: {e}")
        else:
            if structure_source in self.structure_sources:
                out_metadata["structure_source"] = structure_source
                out_metadata["structure_description"] = self.structure_sources[structure_source]
            else:
                mavisp_criticals.append(MAVISpCriticalError(f"Invalid structure_source '{structure_source}' in metadata file. "
                        f"Expected one of: {', '.join(self.structure_sources.keys())}"))

        try:
            linker_design = metadata["linker_design"]
        except KeyError as e:
            log.debug(f"Missing optional field in metadata: {e}")
        else:
            if isinstance(linker_design, bool):
                out_metadata["linker_design"] = linker_design
            else:
                mavisp_criticals.append(MAVISpCriticalError(f"Invalid value for linker_design: {linker_design}. Must be True or False."))

        if structure_source in ["PDB", "Mod"]:
            try:
                pdb_id = metadata["pdb_id"]
            except KeyError:
                log.debug(f"'pdb_id' key missing for structure_source '{structure_source}'")
                pdb_id = None
            else:
                if not pdb_id or str(pdb_id).strip().lower() == 'none':
                    mavisp_criticals.append(MAVISpCriticalError(f"structure_source '{structure_source}' requires a non-empty 'pdb_id' field."))
                    pdb_id = None

            out_metadata["pdb_id"] = pdb_id

        if 'cofactors_in_structure' in metadata.keys():
            cofactors = metadata['cofactors_in_structure']
            if cofactors is None:
                cofactors = []
            if not isinstance(cofactors, list):
                cofactors = [ cofactors ]
            if not set(cofactors).issubset(self.supported_cofactors):
                mavisp_criticals.append(MAVISpCriticalError(f"only these cofactors are allowed for cofactors_in_structure: {', '.join(list(self.supported_cofactors))}"))
            else:
                out_metadata['cofactors_in_structure'] = ', '.join(list(cofactors))
        else:
            out_metadata['cofactors_in_structure'] = ''

        return out_metadata, mavisp_criticals

class MAVISpEnsembleMode(MAVISpMode):

    supported_modules = [ CancermutsTable,
                          Pfam,
                          TED,
                          EnsemblePTMs,
                          EnsembleAllosigmaPSNLongRange,
                          EnsembleStability,
                          EnsembleLocalInteractions,
                          EnsembleSAS,
                          EnsembleDenovoPhospho,
                          EnsembleLocalInteractionsDNA,
                          EnsembleLocalInteractionsHomodimer,
                          EnsembleFunctionalDynamics,
                          FunctionalSites,
                          ClinVar,
                          AlphaFoldMetadata,
                          DeMaSk,
                          GEMME,
                          EVE,
                          AlphaMissense,
                          EFoldMine,
                          ExperimentalData ]
    module_order = ['cancermuts', 'pfam', 'ted', 'stability', 'efoldmine', 'local_interactions', 'local_interactions_DNA',
    'local_interactions_homodimers', 'sas', 'ptms', 'denovo_phospho', 'long_range',
    'functional_dynamics', 'functional_sites', 'clinvar', 'alphafold', 'demask',
    'gemme', 'eve', 'alphamissense', 'experimental_data']
    name = 'ensemble_mode'
    supported_metadata = ['uniprot_ac', 'refseq_id', 'ensemble_sources', 'ensemble_size_foldx',
    'ensemble_size_rosetta', 'sampling_functional_dynamics', 'interfaces_functional_dynamics',
    'review_status', 'curators', 'gitbook_entry', 'ensemble_files_osf', 'structure_source',
    'structure_description', 'linker_design', 'pdb_id', 'cofactors_in_structure']
    index_cols = ['system', 'uniprot_ac', 'refseq_id', 'ensemble_sources', 'ensemble_size_foldx',
    'ensemble_size_rosetta',  'sampling_functional_dynamics', 'interfaces_functional_dynamics', 'simulation_length', 'simulation_force_field',
    'review_status', 'curators', 'gitbook_entry', 'ensemble_files_osf', 'structure_source',
    'structure_description', 'linker_design', 'pdb_id', 'cofactors_in_structure']
    index_col_labels = {'system' : "Protein",
                        'uniprot_ac' : 'Uniprot AC',
                        'refseq_id' : "RefSeq ID",
                        'ensemble_sources' : "Ensemble sources",
                        'ensemble_size_foldx' : 'Ensemble sizes (FoldX)',
                        'ensemble_size_rosetta' : 'Ensemble sizes (Rosetta)',
                        'ensemble_files_osf' : 'OSF repository for ensemble data',
                        'sampling_functional_dynamics' : "Sampling methods for functional dynamics",
                        'interfaces_functional_dynamics' : "Regions of interest for functional dynamics",
                        'simulation_length' : 'Simulation length (ns)',
                        'simulation_force_field' : 'Simulation force field',
                        'review_status' : 'Review status',
                        'curators' : 'Curators',
                        'gitbook_entry' : 'GitBook report',
                        'structure_source': 'Structure source',
                        'structure_description': 'Description of structure source',
                        'linker_design': 'Linker design included',
                        'pdb_id' : 'PDB ID',
                        'cofactors_in_structure': 'Cofactors in starting structures'}
    structure_sources = {
        "AFDB": "AlphaFold database",
        "AF3": "AlphaFold3 webserver",
        "AF2": "AlphaFold2 standalone",
        "PDB": "Experimental PDB structure",
        "Mod": "Homology model (PDB template, reconstruction)"}

    supported_cofactors = {"Zn2+", "Mg2+", "ADP", " ATP", "GDP", "GTP", "NADH", "NAD+", "FADH", "FAD+"}

    ensemble_names = {"md", "cabsflex", "bioemu", "nmr", "metad"}

    def parse_metadata(self, data_dir, system):

        out_metadata = {k : None for k in self.supported_metadata}
        mavisp_criticals = []
        structure_source = None
        metadata = self._parse_metadata_file(data_dir, system, self.name)

        try:
            curators = ', '.join(
                [ f"{curator} ({', '.join(metadata['curators'][curator]['affiliation'])})" for curator in metadata['curators'].keys() ]
            )
        except KeyError:
            log.debug("There is no curators field in metadata file")
            curators = None
            mavisp_criticals.append(MAVISpCriticalError("curators field not found in metadata file"))
        out_metadata['curators'] = curators

        for k in ['uniprot_ac', 'refseq_id', 'review_status']:
            try:
                out_metadata[k] = str(metadata[k])
            except KeyError:
                log.debug(f"There is no {k} field in metadata file")
                out_metadata[k] = None
                mavisp_criticals.append(MAVISpCriticalError(f"{k} was not found in the metadata file"))

        for k in ['ensemble_sources', 'ensemble_size_foldx', 'ensemble_size_rosetta', 'simulation_length',
        'simulation_force_field']:
            try:
                if isinstance(metadata[k], list):
                    out_metadata[k] = ", ".join([str(value) for value in metadata[k]])
                else:
                    out_metadata[k] = None
                    mavisp_criticals.append(MAVISpCriticalError(f"{k} was not a list of values in the metadata file"))
            except KeyError:
                out_metadata[k] = None
                mavisp_criticals.append(MAVISpCriticalError(f"{k} was not found in the metadata file"))
        try:
            if not len(set( len(metadata[m]) for m in ["ensemble_size_foldx", "ensemble_size_rosetta", "ensemble_sources"])) == 1:
                mavisp_criticals.append(MAVISpCriticalError(f"the ensemble metadata must all have the same length"))

        except (KeyError, TypeError): #TypeError is raised when the Key is present in metadata, but the value is None
            pass

        for k in ['sampling_functional_dynamics', 'interfaces_functional_dynamics']:
            try:
                out_metadata[k] = str(metadata[k])
            except KeyError:
                out_metadata[k] = ""

        try:
            curators = ', '.join(
                [ f"{curator} ({', '.join(metadata['curators'][curator]['affiliation'])})" for curator in metadata['curators'].keys() ]
            )
        except KeyError:
            log.debug("There is no curators field in metadata file")
            curators = None
            mavisp_criticals.append(MAVISpCriticalError("curators field not found in metadata file"))

        for k in ['ensemble_files_osf', 'gitbook_entry']:
            if k not in metadata.keys():
                metadata[k] = ''
            out_metadata[k] = metadata[k]

        try:
            structure_source = metadata["structure_source"]
        except KeyError as e:
            log.debug(f"Missing optional field in metadata: {e}")
        else:
            if structure_source in self.structure_sources:
                out_metadata["structure_source"] = structure_source
                out_metadata["structure_description"] = self.structure_sources[structure_source]
            else:
                mavisp_criticals.append(MAVISpCriticalError(f"Invalid structure_source '{structure_source}' in metadata file. "
                        f"Expected one of: {', '.join(self.structure_sources.keys())}"))

        try:
            linker_design = metadata["linker_design"]
        except KeyError as e:
            log.debug(f"Missing optional field in metadata: {e}")
        else:
            if isinstance(linker_design, bool):
                out_metadata["linker_design"] = linker_design
            else:
                mavisp_criticals.append(MAVISpCriticalError(f"Invalid value for linker_design: {linker_design}. Must be True or False."))

        if structure_source in ["PDB", "Mod"]:
            try:
                pdb_id = metadata["pdb_id"]
            except KeyError:
                log.debug(f"'pdb_id' key missing for structure_source '{structure_source}'")
                pdb_id = None
            else:
                if not pdb_id or str(pdb_id).strip().lower() == 'none':
                    mavisp_criticals.append(MAVISpCriticalError(f"structure_source '{structure_source}' requires a non-empty 'pdb_id' field."))
                    pdb_id = None

            out_metadata["pdb_id"] = pdb_id
        
        if 'cofactors_in_structure' in metadata.keys():
            cofactors = metadata['cofactors_in_structure']
            out_string = []
            if cofactors is None:
                cofactors = {}
            if not isinstance(cofactors, dict):
                mavisp_criticals.append(MAVISpCriticalError(f"cofactors_in_structure must be a dictionary, with ensemble names as keys and either a single cofactor or a list of cofactors as value"))
            else:
                for k in cofactors:
                    if cofactors[k] is None:
                        cofactors[k] = []
                    if not isinstance(cofactors[k], list):
                        cofactors[k] = [ cofactors[k] ]
                    if k not in self.ensemble_names:
                        mavisp_criticals.append(MAVISpCriticalError(f"only these ensembles are allowed for cofactors_in_structure: {', '.join(list(self.ensemble_names))}"))
                    elif not set(cofactors[k]).issubset(self.supported_cofactors):
                        mavisp_criticals.append(MAVISpCriticalError(f"only these cofactors are allowed for cofactors_in_structure: {', '.join(list(self.supported_cofactors))}"))
                    else:
                        out_string.append(f"{k}: {', '.join(list(cofactors[k]))}")
                out_metadata['cofactors_in_structure'] = '; '.join(out_string)
        else:
            out_metadata['cofactors_in_structure'] = ''

        return out_metadata, mavisp_criticals


