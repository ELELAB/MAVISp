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
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from io import BytesIO

mode = 'ensemble_mode'

st.set_page_config(layout="wide",
    page_title="MAVISp ensemble mode",
    page_icon="📖")

database_dir = get_database_dir()

add_mavisp_logo("static/logo_small.png")

add_affiliation_logo()

st.header("MAVISp ensemble mode")

st.write('''Please choose a dataset in the table below by clicking on the
checkbox on the left. The corresponding MAVISp results table will be displayed underneath.
Alternatively, check out the Classification tab for a graphical representation
of the final classification performed by MAVISp.''')

st.write('''Please see the Acknowledgement and data usage page for information on our data sources, licensing term, and data reuse permissions''')

try:
    show_table = load_main_table(database_dir, mode)
except FileNotFoundError:
    st.write('No entries are currently available for ensemble mode.')
    st.stop()


gb_datasets_grid = GridOptionsBuilder.from_dataframe(show_table)

if "id_row" not in st.session_state:
    st.session_state["id_row"] = ''
    st.session_state.selected_row = ''
else:
    try:
        st.session_state.selected_row = st.session_state["id_row"].get('selectedItems')[0]['_selectedRowNodeInfo']['nodeRowIndex']
    except:
        pass

gb_datasets_grid.configure_selection(selection_mode='single',
                                     use_checkbox=True)
gb_datasets_grid.configure_column('OSF repository for ensemble data', cellRenderer=cell_renderers['OSF repository for ensemble data'])
gb_datasets_grid.configure_column('GitBook report', cellRenderer=cell_renderers['GitBook report'])

datasets_grid = AgGrid(show_table,
                      gridOptions=gb_datasets_grid.build(),
                      update_mode=GridUpdateMode.SELECTION_CHANGED,
                      fit_columns_on_grid_load = True,
                      reload_data=False,
                      allow_unsafe_jscode=True,
                      height=200,
                      key="id_row")

if len(datasets_grid["selected_rows"]) == 1:

    protein = datasets_grid["selected_rows"][0]['Protein']
    upac = datasets_grid["selected_rows"][0]['Uniprot AC']

    st.write(f"Currently viewing: {protein}")

    this_dataset = load_dataset(database_dir, protein, mode)
    dataset, dotplots, lolliplots = st.tabs(["Dataset", "Classification", "Damaging mutations"])

    ptm_columns = [c for c in this_dataset.columns if c.startswith('PTMs')]

    with dataset:

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
        col = 'Mutation sources'
        this_gb.configure_column(col, cellRenderer=cell_renderers[col])
        for col in ptm_columns:
            this_gb.configure_column(col, cellRenderer=cell_renderers['PTMs'])

        mutations_grid = AgGrid(this_dataset,
                                gridOptions=this_gb.build(),
                                reload_data=False,
                                allow_unsafe_jscode=True,
                                height=400)

    with dotplots:

        this_dataset_table = this_dataset.copy()
        this_dataset_table = this_dataset_table.set_index('Mutation')

        col1, col2 = st.columns(2)

        with col1:
            do_revel = st.checkbox('Show available REVEL classification', )
            revel_co = st.number_input("Cutoff for REVEL score (between 0 and 1)", value=0.5, min_value=0.0, max_value=1.0)
            demask_co = st.number_input("Cutoff for DeMaSk score (absolute value)", value=0.3, min_value=0.0)
            gemme_co = st.number_input("Cutoff for GEMME", value=0.3)
        with col2:
            n_muts = st.number_input("Number of mutations per plot", value=50, min_value=0)
            fig_width = st.number_input("Plot width (inches)", value=14, min_value=0)
            fig_height = st.number_input("Plot height (inches)", value=6, min_value=0)

        plots = plot_dotplot(this_dataset_table,
                             demask_co=demask_co,
                             revel_co=revel_co,
                             gemme_co=gemme_co,
                             n_muts=n_muts,
                             fig_width=fig_width,
                             fig_height=fig_height,
                             do_revel=do_revel,
                             do_demask=do_demask)

        with BytesIO() as pdf_stream:
            with PdfPages(pdf_stream) as pdf:
                for fig in plots:
                    fig.savefig(pdf, format='pdf', dpi=300)

            st.download_button(label="Download as PDF",
                               data=pdf_stream,
                               mime="application/pdf",
                               file_name=f'{protein}-{mode}_dotplots.pdf')

        for fig in plots:
            st.pyplot(fig, use_container_width=False)

    with lolliplots:

        this_dataset_table = this_dataset.copy()
        this_dataset_table = this_dataset_table.set_index('Mutation')

        plots = plot_lolliplots(this_dataset_table)

        with BytesIO() as pdf_stream:
            with PdfPages(pdf_stream) as pdf:
                for fig in plots:
                    fig.savefig(pdf, format='pdf', dpi=300)

            st.download_button(label="Download as PDF",
                               data=pdf_stream,
                               mime="application/pdf",
                               file_name=f'{protein}-{mode}_lolliplots.pdf')

        for fig in plots:
            st.pyplot(fig, use_container_width=False)