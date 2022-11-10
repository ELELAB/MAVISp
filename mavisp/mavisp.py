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
from mavisp.core import MAVISpFileSystem
import logging as log
from termcolor import colored

header = """
                        .__                 
  _____  _____   ___  __|__|  ____________  
 /     \ \__  \  \  \/ /|  | /  ___/\____ \ 
|  Y Y  \ / __ \_ \   / |  | \___ \ |  |_> >
|__|_|  /(____  /  \_/  |__|/____  >|   __/ 
      \/      \/                 \/ |__|    


============================================

If you use MAVISp for your research, please cite:

    Matteo Arnaudi, Ludovica Beltrame, Kristine 
    Degn, Mattia Utichi, Alberto Pettenella, 
    Simone Scrima, Peter Wad Sackett, Matteo 
    Lambrughi, Matteo Tiberti, Elena Papaleo. 
    MAVISp: Multi-layered Assessment of VarIants
    by Structure for proteins. bioRxiv 2022.10.22.513328;
    doi: https://doi.org/10.1101/2022.10.22.513328

============================================


"""

legend = f"""legend:
{colored("    ✓ module_name", 'green')}    all good - no warnings or errors for this module
{colored("    X module_name", 'red')}    errors detected in this module. There might be additional warnings as well
{colored("    ~ module_name", 'yellow')}    warnings detected in this module, with no errors
    = module_name    this module is not available for this protein"""

module_order = ['stability', 'local_interactions', 'cancermuts', 'references', 'ptms', 'long_range', 'clinvar', 'alphafold']

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
    parser.add_argument("-w", "--stop-on-warnings", 
                        dest="stop_on_warnings",
                        default=False,
                        help="do not write output if any warning is found (default: false)")
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
                        help="only perform check, do not write output files")
    parser.add_argument("-f", "--force",
                        dest="force_write",
                        default=False,
                        action="store_true",
                        help="do not stop if output directory exists and overwrite files if necessary")

    args = parser.parse_args()

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
        try:
            out_path.mkdir(exist_ok=True)
        except FileNotFoundError:
            log.error("Coudln't create the specified output directory; parent folder(s) doesn't exist; exiting...")
            exit(1)
        except PermissionError:
            log.error("Coudln't create the specified output directory; exiting...")
            exit(1)
    mfs = MAVISpFileSystem(data_dir=in_path, exclude_proteins=excluded_proteins, include_proteins=included_proteins)

    mfs.ingest()

    summary = mfs.get_datasets_table_summary()
    print("\n\n*** SUMMARY ***\n")
    print(tabulate(summary, headers=summary.columns, showindex=False, tablefmt="double_outline"))

    details = mfs.get_datasets_table_details()

    details_text = "\n\n*** DETAILED REPORT ***\n\n"
    details_text += legend + "\n\n"

    error_count = 0
    warning_count = 0

    details['index'] = list(zip(details['system'], details['mode']))
    systems = details['index'].unique()

    for s in systems:
        this_system = details[(details['system'] == s[0]) & (details['mode'] == s[1])]
        details_text += colored(f"{this_system['system'].unique()[0]} - {this_system['mode'].unique()[0]}\n", 'cyan', attrs=['bold'])
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
    
    if error_count == 0 and warning_count == 0:
        details_text += colored("ALL GOOD! 👏👏👏\n", "green")
    else:
        if error_count > 0:
            error_text = colored(f"{error_count} error(s)", 'red')
        else:
            error_text = ""
        if warning_count > 0:
            warning_text = colored(f"{warning_count} warning(s)", 'yellow')
        else:
            warning_text = ""

        if error_text != "" and warning_text != "":
            separator = ", "
        else:
            separator = ""

        details_text += f"*** {error_text}{separator}{warning_text} ***\n"

    print(details_text)

    if args.dry_run:
        log.info("Exiting without writing any file, as dry-run mode is active")
        exit()

    if error_count > 0:
        log.warning("One or more critical error detected. Will not proceed to generate the database. Exiting...")
        exit(1)

    if args.stop_on_warnings and warning_count > 0:
        log.warning("One or more warnings detected, and you asked to stop on warnings. Will not proceed to generate the database. Exiting...")

    out_table = mfs.dataset_table[['system', 'mode', 'curators']]
    out_table = out_table.rename(columns={'system' : "Protein",
                                          'mode'  : "Mode",
                                          'curators' : 'Curators'})
    out_table.to_csv(out_path / 'index.csv', index=False)

    dataset_tables_path = out_path / 'dataset_tables'
    dataset_tables_path.mkdir(exist_ok=True)
    for _, r in mfs.dataset_table.iterrows():

        this_df = pd.DataFrame({'Mutation': r['mutations']})
        this_df = this_df.set_index('Mutation')

        for mod_name in module_order:
            mod = r['modules'][mod_name]
            if mod is None:
                continue
            this_df = this_df.join(mod.get_dataset_view())

        this_df = this_df.fillna(pd.NA)
        this_df.to_csv(dataset_tables_path / f"{r['system']}-{r['mode']}.csv")
