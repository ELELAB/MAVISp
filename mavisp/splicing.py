"""Attach spliceAI_lookup delta scores to the final genomic annotations."""
from pathlib import Path
import re
import ast
import yaml
from functools import lru_cache
from urllib.parse import urlencode
from urllib.request import urlopen
import xml.etree.ElementTree as ET

import pandas as pd


SPLICING_OUTPUTS = {
    'pangolin_output.csv': (
        ('splice_loss', 'Pangolin Splice_Loss Δ_score'),
        ('splice_gain', 'Pangolin Splice_Gain Δ_score'),
    ),
    'spliceai_output.csv': (
        ('acceptor_loss', 'SpliceAI Acceptor_Loss Δ_score'),
        ('acceptor_gain', 'SpliceAI Acceptor_Gain Δ_score'),
        ('donor_loss', 'SpliceAI Donor_Loss Δ_score'),
        ('donor_gain', 'SpliceAI Donor_Gain Δ_score'),
    ),
}


@lru_cache(maxsize=None)
def _transcript_from_protein(protein_id):
    """Resolve the RefSeq transcript from the protein's NCBI coded_by field.

    Cache within this process so simple/ensemble modes reuse the lookup.
    Do not guess when the record has no unique NM transcript.
    """
    if not re.fullmatch(r'NP_\d+(?:\.\d+)?', protein_id):
        raise ValueError(f'Expected a RefSeq NP accession, got {protein_id!r}')
    query = urlencode({
        'db': 'protein', 'id': protein_id,
        'rettype': 'gp', 'retmode': 'xml', 'tool': 'mavisp_splicing',
    })
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + query
    try:
        with urlopen(url, timeout=30) as response:
            root = ET.fromstring(response.read())
    except Exception as error:
        raise ValueError(
            f'Cannot retrieve NCBI protein record for {protein_id}: {error}'
        ) from error
    transcripts = set()
    for qualifier in root.findall('.//GBQualifier'):
        if qualifier.findtext('GBQualifier_name') == 'coded_by':
            coded_by = qualifier.findtext('GBQualifier_value', default='')
            transcripts.update(re.findall(r'NM_\d+(?:\.\d+)?', coded_by))
    if len(transcripts) != 1:
        raise ValueError(
            f'{protein_id}: expected one NM transcript in NCBI coded_by; '
            f'found {sorted(transcripts)}'
        )
    transcript = transcripts.pop()
    return transcript


def _coordinate(value):
    # Preserve genome assembly: hg19 and hg38 are distinct lookup keys.
    return re.sub(r'\s+', '', str(value))


def _annotations(value):
    if pd.isna(value) or not str(value).strip():
        return []
    # The comma within "hg38,17:g..." belongs to a single annotation.
    parts = [part.strip() for part in str(value).split(',')]
    result = []
    i = 0
    while i < len(parts):
        if re.fullmatch(r'hg\d+', parts[i], flags=re.IGNORECASE):
            if i + 1 >= len(parts) or not parts[i + 1]:
                raise ValueError(f'Incomplete genomic annotation: {value!r}')
            result.append(_coordinate(parts[i] + ',' + parts[i + 1]))
            i += 2
        else:
            result.append(_coordinate(parts[i]))
            i += 1
    return result


def add_splicing_scores(dataset, mode_path):
    """Read optional splicing files and match Mutation + exact HGVSg.

    Multiple genomic annotations receive comma-separated scores in HGVSg
    order, with literal NA placeholders for missing matches. Repeated input
    rows are accepted only if their delta scores agree; conflicting transcript
    predictions require curation rather than an implicit aggregation rule.
    """
    folder = Path(mode_path) / 'splicing'
    if not folder.is_dir():
        return dataset
    metadata_path = Path(mode_path) / 'metadata.yaml'
    with metadata_path.open() as handle:
        metadata = yaml.safe_load(handle) or {}
    protein_id = metadata.get('refseq_id')
    if not protein_id:
        raise ValueError(f'{metadata_path}: missing refseq_id')
    transcript_id = _transcript_from_protein(str(protein_id).strip())

    dataset = dataset.copy()
    hgvs_column = next((c for c in ('HGVSg', 'HGVS') if c in dataset), None)
    for filename, output_columns in SPLICING_OUTPUTS.items():
        lookup = {}
        source = folder / filename
        if source.is_file():
            data = pd.read_csv(source, dtype=str)
            required = {'Mutation', 'variant_coordinate', 'Δ_type', 'Δ_score', 'ref_seq_id'}
            missing = required.difference(data.columns)
            if missing:
                raise ValueError(f'{source}: missing columns {sorted(missing)}')
            def matches_transcript(value):
                if pd.isna(value):
                    return False
                try:
                    refs = ast.literal_eval(value)
                except (ValueError, SyntaxError) as error:
                    raise ValueError(f'{source}: invalid ref_seq_id {value!r}') from error
                if not isinstance(refs, list):
                    raise ValueError(f'{source}: ref_seq_id must contain a list: {value!r}')
                return transcript_id in refs

            data = data.loc[data['ref_seq_id'].map(matches_transcript)]
            if data.empty:
                raise ValueError(f'{source}: no rows for transcript {transcript_id}')
            accepted = {kind for kind, _ in output_columns}
            for _, row in data.iterrows():
                kind = str(row['Δ_type']).strip().casefold()
                if kind not in accepted or pd.isna(row['Δ_score']):
                    continue
                if pd.isna(row['Mutation']) or pd.isna(row['variant_coordinate']):
                    raise ValueError(f'{source}: missing Mutation or variant_coordinate')
                score = str(row['Δ_score']).strip()
                numeric_score = float(score)
                if not pd.notna(numeric_score) or abs(numeric_score) == float('inf'):
                    raise ValueError(f'{source}: non-finite Δ_score {score!r}')
                key = (str(row['Mutation']).strip(), _coordinate(row['variant_coordinate']), kind)
                if key in lookup and float(lookup[key]) != numeric_score:
                    raise ValueError(f'{source}: conflicting Δ_scores for {key}; curate transcript predictions')
                lookup.setdefault(key, score)
        for kind, column in output_columns:
            values = []
            for mutation, row in dataset.iterrows():
                annotations = _annotations(row[hgvs_column]) if hgvs_column else []
                scores = [lookup.get((str(mutation).strip(), annotation, kind)) for annotation in annotations]
                if not scores or all(score is None for score in scores):
                    values.append(pd.NA)
                elif len(scores) == 1:
                    values.append(scores[0])
                else:
                    values.append(', '.join(score if score is not None else 'NA' for score in scores))
            dataset[column] = values
    return dataset


class Splicing:
    """Splicing module participating in MAVISp validation and reporting."""
    name = 'splicing'
    module_dir = 'splicing'

    def __init__(self, data_dir):
        self.data_dir = data_dir
        self.data = None
        self.metadata = None

    def get_dataset_view(self):
        return self.data

    def get_metadata_view(self):
        return self.metadata

    def ingest(self, mutations, genomic_annotations=None):
        from mavisp.error import (
            MAVISpMultipleError, MAVISpCriticalError, MAVISpWarningError,
        )
        warnings = []
        try:
            if genomic_annotations is None or 'HGVSg' not in genomic_annotations:
                raise ValueError('Splicing requires HGVSg from the cancermuts module')
            folder = Path(self.data_dir) / self.module_dir
            available = [name for name in SPLICING_OUTPUTS if (folder / name).is_file()]
            if not available:
                raise ValueError(f'{folder}: neither expected splicing CSV is present')
            for name in SPLICING_OUTPUTS:
                if name not in available:
                    warnings.append(MAVISpWarningError(f'{name} not found; corresponding scores are NA'))
            base = genomic_annotations[['HGVSg']].reindex(mutations)
            annotated = add_splicing_scores(base, self.data_dir)
            columns = [column for outputs in SPLICING_OUTPUTS.values() for _, column in outputs]
            self.data = annotated[columns]
            with (Path(self.data_dir) / 'metadata.yaml').open() as handle:
                metadata = yaml.safe_load(handle) or {}
            protein = str(metadata['refseq_id']).strip()
            self.metadata = {
                'refseq_id': protein,
                'refseq_transcript_id': _transcript_from_protein(protein),
            }
            has_hgvs = base['HGVSg'].notna() & base['HGVSg'].astype(str).str.strip().ne('')
            for filename, outputs in SPLICING_OUTPUTS.items():
                if filename not in available:
                    continue
                for _, column in outputs:
                    values = self.data.loc[has_hgvs, column]
                    missing = values.isna() | values.fillna('').astype(str).str.contains(r'(?:^|,\s*)NA(?:,|$)', regex=True)
                    count = int(missing.sum())
                    if count:
                        warnings.append(MAVISpWarningError(
                            f'{column}: {count} rows with HGVSg have missing or partial scores'
                        ))
        except Exception as error:
            raise MAVISpMultipleError(
                warning=warnings, critical=[MAVISpCriticalError(str(error))]
            ) from error
        if warnings:
            raise MAVISpMultipleError(warning=warnings, critical=[])
