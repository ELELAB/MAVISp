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
                                PTMs,
                          LongRange,
                          Stability,
                          LocalInteractions,
                          SAS,
                          LocalInteractionsDNA,
                          LocalInteractionsHomodimer,
                          ClinVar,
                          AlphaFoldMetadata,
                          DeMaSk,
                          GEMME,
                          EVE ]
    module_order = ['stability', 'local_interactions', 'local_interactions_DNA', 'local_interactions_homodimers', 'sas', 'cancermuts', 'ptms', 'long_range', 'clinvar', 'alphafold', 'demask', 'gemme', 'eve']
    supported_metadata = ['uniprot_ac', 'refseq_id', 'review_status', 'curators']
    index_cols = ['system', 'uniprot_ac', 'refseq_id', 'mode', 'review_status', 'curators']
    index_col_labels = {'system' : "Protein",
                        'mode'  : "Mode",
                        'uniprot_ac' : 'Uniprot AC',
                        'refseq_id' : "RefSeq ID",
                        'review_status' : 'Review status',
                        'curators' : 'Curators',
                        }

    def parse_metadata(self, data_dir, system):

        out_metadata = {k : None for k in self.supported_metadata}
        mavisp_criticals = []

        metadata = self._parse_metadata_file(data_dir, system, self.name)

        try:
            curators= ', '.join(
                [ f"{curator} ({', '.join(metadata['curators'][curator]['affiliation'])})" for curator in metadata['curators'].keys() ]
            )
        except KeyError:
            log.debug("There is no curators field in metadata file")
            out_metadata['curators'] = None
            mavisp_criticals.append(MAVISpCriticalError("curators field not found in metadata file"))

        for k in ['uniprot_ac', 'refseq_id', 'review_status']:
            try:
                out_metadata[k] = str(metadata[k])
            except KeyError:
                log.debug(f"There is no {k} field in metadata file")
                out_metadata[k] = None
                mavisp_criticals.append(MAVISpCriticalError(f"{k} was not found in the metadata file"))

        return out_metadata, mavisp_criticals

class MAVISpEnsembleMode(MAVISpMode):

    supported_modules = [ CancermutsTable,
                                PTMs,
                          LongRange,
                          Stability,
                          LocalInteractions,
                          SAS,
                          LocalInteractionsDNA,
                          LocalInteractionsHomodimer,
                          ClinVar,
                          AlphaFoldMetadata,
                          DeMaSk,
                          GEMME,
                          EVE ]
    module_order = ['stability', 'local_interactions', 'local_interactions_DNA', 'local_interactions_homodimers', 'sas', 'cancermuts', 'ptms', 'long_range', 'clinvar', 'alphafold', 'demask', 'gemme', 'eve']
    name = 'ensemble_mode'
    supported_modules = None
    supported_metadata = ['uniprot_ac', 'refseq_id', 'ensemble_sources', 'ensemble_size_foldx', 'ensemble_size_rosetta', 'review_status', 'curators']
    index_cols = ['system', 'uniprot_ac', 'refseq_id', 'mode','ensemble_sources','ensemble_size_foldx','ensemble_size_rosetta','review_status', 'curators']
    index_col_labels = {'system' : "Protein",
                        'mode'  : "Mode",
                        'uniprot_ac' : 'Uniprot AC',
                        'refseq_id' : "RefSeq ID",
                        'ensemble_sources' : "Ensemble sources",
                        'ensemble_size_foldx' : 'Ensemble sizes (FoldX)',
                        'ensemble_size_rosetta' : 'Ensemble sizes (Rosetta)',
                        'review_status' : 'Review status',
                        'curators' : 'Curators',
                        }


    def parse_metadata(self, data_dir, system):

        out_metadata = {k : None for k in self.supported_metadata}
        mavisp_criticals = []

        metadata = self._parse_metadata_file(data_dir, system, self.name)

        try:
            curators= ', '.join(
                [ f"{curator} ({', '.join(metadata['curators'][curator]['affiliation'])})" for curator in metadata['curators'].keys() ]
            )
        except KeyError:
            log.debug("There is no curators field in metadata file")
            out_metadata['curators'] = None
            mavisp_criticals.append(MAVISpCriticalError("curators field not found in metadata file"))

        for k in ['uniprot_ac', 'refseq_id', 'review_status']:
            try:
                out_metadata[k] = str(metadata[k])
            except KeyError:
                log.debug(f"There is no {k} field in metadata file")
                out_metadata[k] = None
                mavisp_criticals.append(MAVISpCriticalError(f"{k} was not found in the metadata file"))

        for k in ['ensemble_sources', 'ensemble_size_foldx', 'ensemble_size_rosetta']:
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

        return out_metadata, mavisp_criticals