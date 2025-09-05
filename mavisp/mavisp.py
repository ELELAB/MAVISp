#!/usr/bin/env python3
# MAVISp - main data ingest and conversion script
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

from pathlib import Path
import argparse
from tabulate import tabulate
import pandas as pd
import yaml
from mavisp.core import MAVISpFileSystem
from mavisp.utils import mutation_to_HGVSp
import logging as log
from termcolor import colored
from time import strftime
from time import gmtime
header = """
                        .__
  _____  _____   ___  __|__|  ____________
 /     \ \__  \  \  \/ /|  | /  ___/\____ \\
|  Y Y  \ / __ \_ \   / |  | \___ \ |  |_> >
|__|_|  /(____  /  \_/  |__|/____  >|   __/
      \/      \/                 \/ |__|


If you use MAVISp for your research, please cite:

    Matteo Arnaudi, Mattia Utichi, Kristine Degn, Ludovica Beltrame et al.
    MAVISp: A Modular Structure-Based Framework for Protein Variant Effects
    bioRxiv, https://doi.org/10.1101/2022.10.22.513328

for more information about MAVISp see:
    https://www.github.com/ELELAB/MAVISp
    https://services.healthtech.dtu.dk/services/MAVISp-1.0/


============================================


"""

legend = f"""legend:
{colored("    ! error", 'magenta')}          the entry has errors in key components and will not be annotated at all
{colored("    ✓ module_name", 'green')}    all good - no warnings or errors for this module
{colored("    X module_name", 'red')}    errors detected in this module. There might be additional warnings as well
{colored("    ~ module_name", 'yellow')}    warnings detected in this module, with no errors
    = module_name    this module is not available for this protein"""

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--data-dir",
                        dest="data_dir",
                        default="./data",
                        help="directory where the MAVISp data files are located (default: ./data)")
    parser.add_argument("-o", "--database-dir",
                        dest="database_dir",
                        default="./database",
                        help="output directory where the csv database is written (default: ./database")
    parser.add_argument("-m", "--mode",
                        dest="modes",
                        choices=['all', 'simple_mode', 'ensemble_mode'],
                        default="all",
                        help="which mode should considered (all, simple_mode, ensemble_mode; default: all)")
    parser.add_argument("-p", "--included-proteins",
                        dest="included_proteins",
                        default=None,
                        help="list of proteins to be included (comma-separated). The remaining ones won't be considered. Default: all of them")
    parser.add_argument("-e", "--excluded-proteins",
                        dest="excluded_proteins",
                        default=None,
                        help="list of proteins to be excluded (comma-separated). All the remaining ones will be considered. Default: none of them")
    parser.add_argument("-n", "--dry-run",
                        dest="dry_run",
                        default=False,
                        action="store_true",
                        help="perform processing as requested, but do not write output files")
    parser.add_argument("-f", "--force",
                        dest="force_write",
                        default=False,
                        action="store_true",
                        help="do not stop if output directory exists and overwrite files if necessary")
    parser.add_argument("-v", "--verbose",
                        dest="verbose",
                        default=False,
                        action="store_true",
                        help="toggle verbose mode")



    args = parser.parse_args()

    if args.verbose:
        log.basicConfig(level=log.DEBUG)
    else:
        log.basicConfig(level=log.WARNING)

    if args.excluded_proteins is not None and args.included_proteins is not None:
        log.error("options -p and -e cannot be specified at the same time; exiting...")

    if args.excluded_proteins is not None:
        excluded_proteins = args.excluded_proteins.split(',')
    else:
        excluded_proteins = None

    if args.included_proteins is not None:
        included_proteins = args.included_proteins.split(',')
    else:
        included_proteins = None

    in_path = Path(args.data_dir)

    if not in_path.is_dir():
        log.error("Input directory isn't a directory or doesn't exist; exiting...")
        exit(1)

    if not args.dry_run:
        out_path = Path(args.database_dir)

        if out_path.exists():
            if out_path.is_file():
                log.error("Specified output directory is currently a file; exiting...")
                exit(1)

            if out_path.is_dir() and not args.force_write:
                log.error("Specified output directory already exists; exiting...")
                exit(1)

    if args.modes == 'all':
        args.modes = None

    print(header)

    mfs = MAVISpFileSystem( data_dir=in_path,
                            exclude_proteins=excluded_proteins,
                            include_proteins=included_proteins,
                            modes=args.modes)

    mfs.ingest()

    all_modes_error_count = 0
    all_modes_critical_count = 0

    for mode_name in mfs.supported_modes.keys():
        summary = mfs.get_datasets_table_summary(mode_name)

        print(f"\n\n*** SUMMARY - {mode_name}***\n")

        if len(summary) == 0:
            print(colored(f"No entry found for {mode_name}\n", 'magenta'))
            continue

        print(tabulate(summary, headers=summary.columns, showindex=False, tablefmt="double_outline"))

        details_text = "\n\n*** DETAILED REPORT ***\n\n"
        details_text += legend + "\n\n"

        error_count = 0
        warning_count = 0
        critical_count = 0

        details = mfs.get_datasets_table_details(mode_name)
        systems = details['system'].unique()

        for s in systems:
            this_system = details[details['system'] == s].reset_index()
            if 'critical' in this_system['status'].values:
                details_text += colored(f"{this_system['system'].unique()[0]} - {mode_name}\n", 'magenta', attrs=['bold'])
                for crit in this_system.iloc[0]['details_crit']:
                    details_text += colored(f"    ! {str(crit)}\n", 'magenta')
                    critical_count += 1
                details_text += '\n'
                continue

            details_text += colored(f"{this_system['system'].unique()[0]} - {mode_name}\n", 'cyan', attrs=['bold'])
            for _, r in this_system.iterrows():
                if r['status'] == "ok":
                    details_text += colored(f"    ✓ {r['module']}\n", 'green')
                if r['status'] == "error":
                    details_text += colored(f"    X {r['module']}\n", 'red')
                    for err in r['details_err']:
                        details_text += f"        {colored(err, 'red')}\n"
                        error_count += 1
                    for warn in r['details_warn']:
                        details_text += f"        {colored(warn, 'yellow')}\n"
                        warning_count += 1
                if r['status'] == "warning":
                    details_text += colored(f"    ~ {r['module']}\n", 'yellow')
                    for warn in r['details_warn']:
                        details_text += f"        {colored(warn, 'yellow')}\n"
                        warning_count += 1
                if r['status'] == "not_available":
                    details_text += f"    = {r['module']}\n"
            details_text += '\n'

        if error_count == 0 and warning_count == 0 and critical_count == 0:
            details_text += colored("*** ALL GOOD! 👏👏👏 ***\n", "green")
        else:
            if critical_count > 0:
                critical_text = colored(f"{critical_count} critical error(s)", 'magenta')
            else:
                critical_text = ""

            if error_count > 0:
                error_text = colored(f"{error_count} error(s)", 'red')
            else:
                error_text = ""

            if warning_count > 0:
                warning_text = colored(f"{warning_count} warning(s)", 'yellow')
            else:
                warning_text = ""

            details_text += f"*** {', '.join( [ x for x in [ critical_text, error_text, warning_text ] if x != ''] ) } ***\n"

        all_modes_error_count += error_count
        all_modes_critical_count += critical_count

        print(details_text)

    if args.dry_run:
        log.info("Exiting without writing any file, as dry-run mode is active")
        exit()

    if all_modes_error_count > 0 or all_modes_critical_count > 0:
        log.error("One or more error detected. Will not proceed to generate the database. Exiting...")
        exit(1)

    try:
        out_path.mkdir(exist_ok=True)
    except FileNotFoundError:
        log.error("Couldn't create the specified output directory; parent folder(s) doesn't exist; exiting...")
        exit(1)
    except PermissionError:
        log.error("Couldn't create the specified output directory; exiting...")
        exit(1)

    all_indexes = []

    for mode_name, mode in mfs.supported_modes.items():

        if len(mfs.dataset_tables[mode_name]) == 0:
            continue

        out_table = mfs.dataset_tables[mode_name].copy(deep=True)

        out_index_table = out_table.copy(deep=True)
        out_index_table['mode'] = mode_name
        all_indexes.append(out_index_table)

        mode_path = out_path / Path(mode_name)
        mode_path.mkdir(exist_ok=True)

        out_table = out_table[mode.index_cols]
        out_table = out_table.rename(columns=mode.index_col_labels)

        out_table.to_csv(mode_path / 'index.csv', index=False)

        dataset_tables_path = mode_path / 'dataset_tables'
        metadata_path = mode_path / 'metadata'

        dataset_tables_path.mkdir(exist_ok=True)
        metadata_path.mkdir(exist_ok=True)

        for _, r in mfs.dataset_tables[mode_name].iterrows():

            rows_n = {}

            this_refseq_id = out_table[out_table['Protein'] == r['system']]
            assert(this_refseq_id.shape[0]) == 1
            this_refseq_id = this_refseq_id.iloc[0]['RefSeq ID']

            this_df = r['mutations']
            this_df = this_df.rename(columns={'mutation' : 'Mutation',
                                              'PMID'     : 'References'})
            this_df['HGVSp'] = this_df.apply(lambda r: f"{this_refseq_id}:{mutation_to_HGVSp(r['Mutation'])}", axis=1)
            this_df = this_df.set_index('Mutation')

            module_metadata = {}

            for mod_name in mode.module_order:
                mod = r['modules'][mod_name]
                if mod is None:
                    continue
                this_df = this_df.join(mod.get_dataset_view())
                rows_n[mod_name] = len(this_df)
                module_metadata[mod.name] = mod.get_metadata_view()

            # move Reference column to last
            this_df = this_df[[c for c in this_df.columns if c != 'References'] + ['References']]

            # fill NAs
            this_df = this_df.fillna(pd.NA)

            # save final dataframe
            this_df.to_csv(dataset_tables_path / f"{r['system']}-{mode.name}.csv")

            if any(this_df.index.duplicated()):
                log.warning(f"duplicated mutations found in {r['system']}, {mode.name}")
                log.debug("number of rows after each module:")
                log.debug(rows_n)

            # save metadata
            with open(metadata_path / f"{r['system']}.yaml", 'w') as fh:
                yaml.dump(module_metadata, fh, sort_keys=False)

    all_indexes = pd.concat(all_indexes, ignore_index=True)

    # Generate a csv file that contains the number of mutations and the date of the run
    time = strftime("%Y-%m-%d", gmtime())

    # turn mutations to list
    #print(all_indexes.iloc[0])
    #print(all_indexes.iloc[0]['mutations'])

    all_indexes['mutations'] = all_indexes.apply(lambda r: r.mutations['mutation'].tolist(), axis=1)

    # Count number of unique mutations
    nb_mutations = all_indexes.explode('mutations').drop_duplicates(['system', 'mutations']).shape[0]

    # Group the rows by the "system" column and count the number of unique modes for each group
    grouped = all_indexes.groupby('system')['mode'].nunique()

    # Group and count instances in which a protein is found to have two modes
    nb_both_modes = grouped[grouped == 2].shape[0]

    nb_simple_mode =   len(all_indexes[all_indexes['mode'] == 'simple_mode']) - nb_both_modes
    nb_ensemble_mode = len(all_indexes[all_indexes['mode'] == 'ensemble_mode']) - nb_both_modes

    nb_proteins = nb_simple_mode + nb_ensemble_mode + nb_both_modes

    mutation_table = pd.DataFrame({ 'Date of run': time,
                                    'Number of mutations': nb_mutations,
                                    'Number of proteins': nb_proteins,
                                    'Number of proteins in simple mode only' : nb_simple_mode,
                                    'Number of proteins in ensemble mode only': nb_ensemble_mode,
                                    'Number of proteins in both modes': nb_both_modes,
                                  }, index=[0])
    mutation_table.to_csv(out_path / 'dataset_info.csv', index=False)
