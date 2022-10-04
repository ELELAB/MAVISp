# MAVISp - classes for data ingest and manipulation
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
from functools import reduce
import pandas as pd
import numpy as np
import logging as log
from datetime import date
from error import *

class DataType(object):
    def __init__(self, data_dir=None):

        self.data_dir = data_dir
        
        if data_dir is None:
            return

        self.ingest()
        self.process()

    def ingest(self):
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
            raise MAVISpCriticalError("the input directory pathway doesn't exist or is not a directory")

    @data_dir.getter
    def data_dir(self):
        return self._data_dir


class CancermutsTable(DataType):
    def __init__(self, data_dir=None):

        super().__init__(data_dir)

    def ingest(self):
        warnings = []

        cancermuts_files = os.listdir(self.data_dir)
        if len(cancermuts_files) != 1:
            this_error = f"multiple files found in {cancermuts_files}; only one expected"
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=MAVISpCriticalError(this_error))
        
        cancermuts_file = cancermuts_files[0]

        try:
            self.data = self._parse(os.path.join(self.data_dir, cancermuts_file))
        except IOError:
            this_error = "failed to parse Cancermuts table file"
            raise MAVISpMultipleError(warning=warnings, 
                                      critical=MAVISpCriticalError(this_error))

    def _parse(self, fname):

        log.info(f"parsing Cancermuts file {fname}")

        cancermuts = pd.read_csv(fname)

        cancermuts = cancermuts[ ~ pd.isna(cancermuts.alt_aa)]
        cancermuts['mutation_index'] = cancermuts.ref_aa + cancermuts.aa_position.astype(str) + cancermuts.alt_aa
        cancermuts = cancermuts.set_index('mutation_index')
        cancermuts = cancermuts[['gnomad_genome_af', 'gnomad_exome_af', 'REVEL_score', 'sources']]

        return       cancermuts.rename(columns={ 'gnomad_genome_af' : 'gnomAD genome allele frequency',
                                                 'gnomad_exome_af'  : 'gnomAD exome allele frequency',
                                                 'REVEL_score'      : 'REVEL score',
                                                 'sources'          : 'Mutation sources' })

    def main_table_view(self):
        return self.data

    def dataset_view(self):
        return self.data


class MAVISpFileSystem:

    supported_modes = ['basic_mode']
    supported_stability_methods = ['foldx5', 'rosetta_ref2015', 'rosetta_cartddg2020_ref2015']
    supported_interaction_methods = ['foldx5']

    def __init__(self, modes=None, proteins=None, data_dir="/data/raw_data/computational_data/mavisp_data/", verbose=True):

        self.log = log.getLogger(__name__)

        if verbose:
            level = log.info
            message = "%(message)"
        else:
            level = log.warning
            message = "%(message)"

        self.level.basicConfig(level=level, format=message)

        self.proteins = proteins
        
        modes_diff = set(modes).difference(set(self.supported_modes))
        if len(modes_diff) > 0: 
            raise TypeError(f"the following modes are not supported: {modes_diff}")

        self._tree = self._traverse(self.data_dir)
        self.dataset_table = self._dataset_table()

    def _dataset_table(self):

        warnings = []

        self.log.info("generating dataset table")

        df_list = []

        for system in  self._dir_list(self._tree):
            if self.protein is not None and system not in self.proteins:
                self.log.warning(f"ignoring unsupported protein {system}")
                warnings.append(f"ignoring unsupported protein {system}")
                continue
            for mode in self._dir_list(self._tree[system]):
                if mode not in self.supported_modes:
                    self.log.warning(f"ignoring unsupported mode {mode} in {[system, mode]}")
                    warnings.append(f"ignoring unsupported mode {mode} in {[system, mode]}")
                    continue
                self.log.info(f"adding {[system, mode]} to dataset")

                try: 
                    mutation_list = self._mutation_list(system, mode)
                except:
                    exit(1)
                    
                df_list.append((system, mode, mutation_list))

        main_df = pd.DataFrame.from_records(df_list, columns=['system', 'mode', 'mutations'])
        log.debug(f"identified datasets:\n{main_df}")

        return main_df

    def _dir_list(self, d):

        return [ k for k,v in d.items() if v is not None ]

    def _file_list(self, d):

        return [ k for k,v in d.items() if v is None ]

    def _traverse(self, rootdir):

        log.info(f"building directory tree for {rootdir}")

        tree = {}
        rootdir = rootdir.rstrip(os.sep)
        start = rootdir.rfind(os.sep) + 1
        for path, dirs, files in os.walk(rootdir, followlinks=True):
            folders = path[start:].split(os.sep)
            subdir = dict.fromkeys(files)
            parent = reduce(dict.get, folders[:-1], tree)
            parent[folders[-1]] = subdir

        return tree[os.path.basename(os.path.normpath(rootdir))]

    def _parse_mutation_list(self, system, mode):
            
        log.info(f"Gathering mutation list for {system} {mode}")

        mutation_files = self._file_list(self._tree[system][mode]['mutation_list'])
        most_recent_mut_file = self._select_most_recent_file(mutation_files)
        log.info(f"selected {most_recent_mut_file} as mutation file")

        mut_path = os.path.join(self.data_dir, system, mode, 'mutation_list', most_recent_mut_file)

        try:
            with open(fname) as fh:
                lines = fh.read().splitlines()
        except IOError:
            log.error("Couldn't parse mutation list {fname}")
            raise IOError

        # remove duplicates, sort, remove empty lines
        lines = sorted(list(set(lines)), key=lambda x: int(x[1:-1]))
        mutations = list(filter(lambda x: len(x) != 0, lines))

        log.debug(f"found mutations: {mutations}")
        
        return mutations

