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
        errors = []

        cancermuts_files = os.listdir(self.data_dir)
        if len(cancermuts_files) != 1:
            this_error = f"multiple files found in {cancermuts_files}; only one expected"
            raise MAVISpMultipleError(warning=errors, 
                                      critical=MAVISpCriticalError(this_error))
        
        cancermuts_file = cancermuts_files[0]

        try:
            self.data = self._parse(os.path.join(self.data_dir, cancermuts_file))
        except IOError:
            this_error = "failed to parse Cancermuts table file"
            raise MAVISpMultipleError(warning=errors, 
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
