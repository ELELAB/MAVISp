# MAVISp - Streamlit application
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

import streamlit as st
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode
from mavisp.ingest import MAVISFileSystem
import pandas as pd

mfs = MAVISFileSystem()

st.image('assets/logo.png')

st.write('Welcome to MAVISp!')

gb_datasets_grid = GridOptionsBuilder.from_dataframe(mfs.dataset_table)
gb_datasets_grid.configure_selection(selection_mode='single',
                       use_checkbox=True)

datasets_grid = AgGrid(mfs.dataset_table,
                      gridOptions=gb_datasets_grid.build(),
                      update_mode=GridUpdateMode.SELECTION_CHANGED,
                      fit_columns_on_grid_load = True)

if len(datasets_grid["selected_rows"]) == 1:
    print(datasets_grid["selected_rows"][0].values())
    print(tuple(datasets_grid["selected_rows"][0].values()))
    this_dataset = mfs.mutation_datasets[tuple(datasets_grid["selected_rows"][0].values())]

    this_gb = GridOptionsBuilder.from_dataframe(this_dataset)
    this_gb.configure_grid_options(alwaysShowHorizontalScroll=True)

    mutations_grid = AgGrid(this_dataset,
                            gridOptions=this_gb.build(),
                            width=3000)

