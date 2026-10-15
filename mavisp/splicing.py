"""Attach spliceAI_lookup delta scores to the final genomic annotations."""
from pathlib import Path
import re

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
    dataset = dataset.copy()
    hgvs_column = next((c for c in ('HGVSg', 'HGVS') if c in dataset), None)
    for filename, output_columns in SPLICING_OUTPUTS.items():
        lookup = {}
        source = folder / filename
        if source.is_file():
            data = pd.read_csv(source, dtype=str)
            required = {'Mutation', 'variant_coordinate', 'Δ_type', 'Δ_score'}
            missing = required.difference(data.columns)
            if missing:
                raise ValueError(f'{source}: missing columns {sorted(missing)}')
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

    def ingest(self, mutations):
        from mavisp.error import (
            MAVISpMultipleError, MAVISpCriticalError, MAVISpWarningError,
        )
        warnings = []
        try:
            table_folder = Path(self.data_dir) / 'cancermuts'
            tables = list(table_folder.iterdir())
            if len(tables) != 1 or not tables[0].is_file():
                raise ValueError(f'{table_folder}: expected exactly one Cancermuts metatable')
            table = pd.read_csv(tables[0])
            required = {'ref_aa', 'aa_position', 'alt_aa', 'genomic_mutation'}
            if not required.issubset(table.columns):
                raise ValueError(f'{tables[0]}: missing columns {sorted(required - set(table.columns))}')
            table = table.loc[table['alt_aa'].notna()].copy()
            table.index = table['ref_aa'] + table['aa_position'].astype(str) + table['alt_aa']
            if table.index.has_duplicates:
                raise ValueError(f'{tables[0]}: duplicate protein mutations')
            genomic_annotations = table[['genomic_mutation']].rename(columns={'genomic_mutation': 'HGVSg'})
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
            self.metadata = {'input': 'pre-filtered splicing outputs'}
            has_hgvs = base['HGVSg'].notna() & base['HGVSg'].astype(str).str.strip().ne('')
            for filename, outputs in SPLICING_OUTPUTS.items():
                if filename not in available:
                    continue
                values = self.data.loc[has_hgvs, [column for _, column in outputs]]
                # Count once per predictor; all-NA rows and partial rows are distinct.
                absent = values.isna().all(axis=1)
                incomplete = values.isna().any(axis=1)
                for column in values.columns:
                    incomplete |= values[column].fillna('').astype(str).str.contains(
                        r'(?:^|,\s*)NA(?:,|$)', regex=True
                    )
                missing_count = int(absent.sum())
                partial_count = int((incomplete & ~absent).sum())
                if missing_count or partial_count:
                    predictor = 'Pangolin' if filename == 'pangolin_output.csv' else 'SpliceAI'
                    warnings.append(MAVISpWarningError(
                        f'{predictor}: splicing annotations were missing or incomplete '
                        f'for {missing_count + partial_count} MAVISp mutations in {filename}.'
                    ))
        except Exception as error:
            raise MAVISpMultipleError(
                warning=warnings, critical=[MAVISpCriticalError(str(error))]
            ) from error
        if warnings:
            raise MAVISpMultipleError(warning=warnings, critical=[])
