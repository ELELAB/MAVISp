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
import os
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode
import pandas as pd
from collections import defaultdict
from streamlit_utils import *

st.set_page_config(layout="wide",
    page_title="Database",
    page_icon="👋")

@st.cache_data
def load_dataset(data_dir, protein, mode):
    return pd.read_csv(os.path.join(data_dir, 'dataset_tables', f'{protein}-{mode}.csv'))

@st.cache_data
def load_main_table(data_dir):
    return pd.read_csv(os.path.join(data_dir, 'index.csv'))

database_dir = get_database_dir()

show_table = load_main_table(database_dir)

add_logo("assets/logo_small.png")

st.write('Please choose a dataset in the table below and click on the "View dataset"'
'button.  The corresponding MAVISp results table will be displayed underneath. '
'Click on the "Download dataset" button to download the corresponding CSV file.')

gb_datasets_grid = GridOptionsBuilder.from_dataframe(show_table)

gb_datasets_grid.configure_selection(selection_mode='single',
                       use_checkbox=True)

datasets_grid = AgGrid(show_table,
                      gridOptions=gb_datasets_grid.build(),
                      update_mode=GridUpdateMode.SELECTION_CHANGED,
                      fit_columns_on_grid_load = True,
                      reload_data=False)

if len(datasets_grid["selected_rows"]) != 1:
    button_disabled = True
else:
    button_disabled = False

if st.button('View dataset',
            disabled=button_disabled):

    protein, mode = ( datasets_grid["selected_rows"][0]['Protein'],
                    datasets_grid["selected_rows"][0]['Mode']    )

    st.write(f"Currently viewing: {protein}, {mode}")

    this_dataset = load_dataset(database_dir, protein, mode)

    with open(os.path.join(database_dir, 'dataset_tables', f'{protein}-{mode}.csv')) as data:
        st.download_button(label="Download dataset",
                            data=data,
                            file_name=f'{protein}-{mode}.csv',
                            mime="text/csv",
                            key='download-csv')

    this_gb = GridOptionsBuilder.from_dataframe(this_dataset)
    this_gb.configure_grid_options(alwaysShowHorizontalScroll=True)

    this_dataset = this_dataset.fillna(pd.NA)
    mutations_grid = AgGrid(this_dataset,
                            gridOptions=this_gb.build(),
                            reload_data=False)

