# MAVISp - Streamlit application
# Copyright (C) 2022 Matteo Tiberti, Danish Cancer Society
#
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

import streamlit as st
import os
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode
import pandas as pd
from collections import defaultdict
from streamlit_utils import *

mode = 'simple_mode'

st.set_page_config(layout="wide",
    page_title="MAVISp simple mode",
    page_icon="📖")

database_dir = get_database_dir()

add_mavisp_logo("static/logo_small.png")

add_affiliation_logo()

st.header('MAVISp simple mode')

st.write('''Please choose a dataset in the table below and click on the "View 
dataset" button. The corresponding MAVISp results table will be displayed underneath.''')

st.write('''Click on the "Download dataset" button to download the corresponding CSV file.''')

st.write('''Please see the Acknowledgement and data usage page for information on our data
sources, licensing term, and data reuse permissions''')

try:
    original_show_table = load_main_table(database_dir, mode)
except FileNotFoundError:
    st.write('No entries are currently available for simple mode.')
    st.stop()

filter_string = st.text_input("Search for protein name, UniProt AC or RefSeq ID")

if filter_string == "":
    show_table = original_show_table
else:
    filter_string = filter_string.upper()
    show_table = original_show_table.copy()
    show_table = show_table[ (show_table['Protein'].str.contains(filter_string)) | (show_table['Uniprot AC'].str.contains(filter_string)) | (show_table['RefSeq ID'].str.contains(filter_string))]

gb_datasets_grid = GridOptionsBuilder.from_dataframe(show_table)

if "id_row" not in st.session_state:
    st.session_state["id_row"] = ''
    st.session_state.selected_row = '0'
else:
    st.session_state.selected_row = st.session_state["id_row"].get('selectedItems')[0]['_selectedRowNodeInfo'][
        'nodeRowIndex']


gb_datasets_grid.configure_selection(selection_mode='single',
                                    pre_selected_rows=[str(st.session_state.selected_row)],
                                       use_checkbox=True)

datasets_grid = AgGrid(show_table,
                      gridOptions=gb_datasets_grid.build(),
                      update_mode=GridUpdateMode.SELECTION_CHANGED,
                      fit_columns_on_grid_load = True,
                      reload_data = True,
                      height=200,
                      key="id_row")

if len(datasets_grid["selected_rows"]) != 1:
    button_disabled = True
else:
    button_disabled = False
    #st.session_stat['id_row'] = datasets_grid['selected_rows']

if st.button('View dataset',
            disabled=button_disabled):

    protein = datasets_grid["selected_rows"][0]['Protein']
    upac = datasets_grid["selected_rows"][0]['Uniprot AC']

    st.write(f"Currently viewing: {protein}")

    this_dataset = load_dataset(database_dir, protein, mode)
    print(this_dataset.columns)

    with open(os.path.join(database_dir, mode, 'dataset_tables', f'{protein}-{mode}.csv')) as data:
        st.download_button(label="Download dataset",
                            data=data,
                            file_name=f'{protein}-{mode}.csv',
                            mime="text/csv",
                            key='download-csv')

    this_dataset_table = this_dataset.copy()
    this_dataset_table = this_dataset_table.fillna(pd.NA)
    this_dataset_table['UniProtAC'] = upac

    this_gb = GridOptionsBuilder.from_dataframe(this_dataset_table)
    this_gb.configure_grid_options(alwaysShowHorizontalScroll=True)
    this_gb.configure_column('UniProt AC', hide=True)
    for col in ['Mutation sources', 'PTMs']:
        this_gb.configure_column(col, cellRenderer=cell_renderers[col])

    mutations_grid = AgGrid(this_dataset_table,
                            gridOptions=this_gb.build(),
                            reload_data=False,
                            allow_unsafe_jscode=True,
                            key = st.session_state.grid_key)

    plots = plot_dotplot(this_dataset, demask_co=0.3, revel_co=0.5)
    st.write(f"N plots: {len(plots)}")
    for plot in plots:
        st.pyplot(plot)
