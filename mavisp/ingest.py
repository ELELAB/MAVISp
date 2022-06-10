# MAVISp - classes for data ingest
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

class MAVISFileSystem:

    allowed_methods=['nmr', 'xray']#, 'alphafold']
    
    def __init__(self, data_dir="/data/raw_data/computational_data/mavisp_data/"):
        self.data_dir = data_dir
        self.dataset_table = None
        self.mutation_lists = None
        self.mutation_datasets = None

        self.ingest_data()

    def ingest_data(self):
        self._tree = self._traverse(self.data_dir)
        self.dataset_table = self._dataset_table()
        self.mutation_lists = self._mutation_lists()
        self.mutation_tables = self._mutation_tables()

    def _dir_list(self, d):
        out = []
        for k,v in d.items():
            if v is not None:
                out.append(k)

        return out
  
    def _file_list(self, d):
        out = []
        for k,v in d.items():
            if v is None:
                out.append(k)
        return out

    def _traverse(self, rootdir):
        tree = {}
        rootdir = rootdir.rstrip(os.sep)
        start = rootdir.rfind(os.sep) + 1
        print(rootdir, start)
        for path, dirs, files in os.walk(rootdir):
            folders = path[start:].split(os.sep)
            subdir = dict.fromkeys(files)
            parent = reduce(dict.get, folders[:-1], tree)
            parent[folders[-1]] = subdir
        return tree[os.path.basename(os.path.normpath(rootdir))]

    def _parse_mutation_list(self, fname):
        out = []
        with open(fname) as fh:
            for line in fh:
                tmp = line.strip()
                #out.append((tmp[0], tmp[1:-1], tmp[-1]))
                out.append(tmp)

        return out

    def _dataset_table(self):
        df_list = {}
        mut_lists = {}
        data_dfs = {}

        k = 0

        for system in  self._dir_list(self._tree):
            for mode in self._dir_list(self._tree[system]):
                if mode == 'full_mode':
                    continue
                for structure in self._dir_list(self._tree[system][mode]):
                    st_id, residues = structure.split('_')
                    for method in self._dir_list(self._tree[system][mode][structure]):
                        if method in self.allowed_methods:
                            for model in self._dir_list(self._tree[system][mode][structure][method]):
                                df_list[k] = [system, mode, st_id, residues, method, model]
                                k += 1

        main_df = pd.DataFrame.from_dict(df_list, orient='index', columns=['system', 'mode', 'structure ID', 'residue range', 'method', 'model'])
        return main_df


    def _mutation_lists(self):
        if self.dataset_table is None:
            return None

        mut_lists = {}
        for idx, r in self.dataset_table.iterrows():
            mut_fname = self._file_list(self._tree[r['system']][r['mode']][f"{r['structure ID']}_{r['residue range']}"][r['method']][r['model']]['mutation_list'])[0]
            mut_path = os.path.join(self.data_dir, r['system'], r['mode'], f"{r['structure ID']}_{r['residue range']}", r['method'], r['model'], 'mutation_list', mut_fname)
            mut_lists[(r['system'], r['mode'], r['structure ID'], r['residue range'], r['method'], r['model'])] = self._parse_mutation_list(mut_path)

        return mut_lists


    def _mutation_tables(self):
        if self.dataset_table is None:
            return None

        data_dfs = {}

        for idx, r in self.dataset_table.iterrows():
            mutations = self.mutation_lists[(r['system'], r['mode'], r['structure ID'], r['residue range'], r['method'], r['model'])]
            this_df = pd.DataFrame({'mutations': mutations})
            this_df = this_df.set_index('mutations')
            analysis_basepath = os.path.join(self.data_dir, r['system'], r['mode'], f"{r['structure ID']}_{r['residue range']}", r['method'], r['model'])
            stability_methods = self._dir_list(self._tree[r['system']][r['mode']][f"{r['structure ID']}_{r['residue range']}"][r['method']][r['model']]['stability'])
            for method in stability_methods:
                if method == 'foldx5':
                    foldx_dir = os.path.join(analysis_basepath, 'stability', 'foldx5')
                    foldx_file = self._file_list(self._tree[r['system']][r['mode']][f"{r['structure ID']}_{r['residue range']}"][r['method']][r['model']]['stability']['foldx5'])[0]
                    data = pd.read_csv(os.path.join(foldx_dir, foldx_file), sep='\t', header=None)
                    data[0] = data[0].apply(lambda x: x[0] + x[2:])
                    data = data[[0,1]].set_index(0).rename(columns={1:'STABILITY (FoldX5, kcal/mol)'})
                    this_df = this_df.join(data)
                if method == 'rosetta_ref2015':
                    data = []
                    rosetta_dir = os.path.join(analysis_basepath, 'stability', 'rosetta_ref2015')
                    for mutation in this_df.index:
                        rosetta_file = os.path.join(rosetta_dir, f"{mutation}_aggregate.csv")
                        try:
                            mutation_data =  pd.read_csv(rosetta_file)
                            data.append(mutation_data[mutation_data['state'] == 'ddg']['total_score'][0])
                        except IOError:
                            data.append(np.nan)
                    this_df['STABILITY (Rosetta, ref2015, kcal/mol)'] = data
            data_dfs[(r['system'], r['mode'], r['structure ID'], r['residue range'], r['method'], r['model'])] = this_df

        return data_dfs


    


"""

print(main_df)
print(mut_lists)
"""
