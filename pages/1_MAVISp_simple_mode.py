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
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from io import BytesIO

mode = 'simple_mode'

st.set_page_config(layout="wide",
    page_title="MAVISp simple mode",
    page_icon="📊")

database_dir = get_database_dir()

add_mavisp_logo("static/logo_small.png")

add_affiliation_logo()

st.header('MAVISp simple mode')

st.write('''Please choose a dataset in the table below by clicking on the
checkbox on the left. The corresponding MAVISp results table will be displayed underneath.
Alternatively, check out the Classification tab for a graphical representation
of the final classification performed by MAVISp.''')

st.write('''Please see the Acknowledgement and data usage page for information on our data sources, licensing term, and data reuse permissions''')

try:
    show_table = load_main_table(database_dir, mode)
except FileNotFoundError:
    st.write('No entries are currently available for simple mode.')
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
gb_datasets_grid.configure_column('GitBook report', cellRenderer=cell_renderers['GitBook report'])

datasets_grid = AgGrid(show_table,
                      gridOptions=gb_datasets_grid.build(),
                      update_mode=GridUpdateMode.SELECTION_CHANGED,
                      fit_columns_on_grid_load = True,
                      reload_data = False,
                      height=200,
                      allow_unsafe_jscode=True,
                      key="id_row",
                      custom_css={"#gridToolBar" : {
                                      "padding-bottom": "0px !important",
                                      }
                                  })

if len(datasets_grid["selected_rows"]) == 1:

    protein = datasets_grid["selected_rows"][0]['Protein']
    upac = datasets_grid["selected_rows"][0]['Uniprot AC']

    st.write(f"Currently viewing: {protein}")

    this_dataset = load_dataset(database_dir, protein, mode)

    dataset, dotplots, lolliplots = st.tabs(["Dataset", "Classification", "Damaging mutations"])

    with dataset:

        data_type = st.radio("Select the data representation to be shown or downloaded",
                            options=["Compact dataset", "Full dataset"],
                            index=0)

        this_dataset_table = this_dataset.copy()
        this_dataset_table = this_dataset_table.fillna(pd.NA)
        this_dataset_table['UniProtAC'] = upac

        if data_type == 'Compact dataset':
            this_dataset_table = get_compact_dataset(this_dataset_table)
            st.download_button(label="Download dataset",
                                data=this_dataset_table.to_csv(),
                                file_name=f'{protein}-{mode}-compact.csv',
                                mime="text/csv",
                                key='download-csv-compact')

        else:
            with open(os.path.join(database_dir, mode, 'dataset_tables', f'{protein}-{mode}.csv')) as data:
                st.download_button(label="Download dataset",
                                    data=data,
                                    file_name=f'{protein}-{mode}.csv',
                                    mime="text/csv",
                                    key='download-csv')

        this_gb = GridOptionsBuilder.from_dataframe(this_dataset_table)
        this_gb.configure_grid_options(alwaysShowHorizontalScroll=True)
        this_gb.configure_column('UniProtAC', hide=True)
        for col in ['Mutation sources', 'PTMs']:
            if col in this_dataset_table.columns:
                this_gb.configure_column(col, cellRenderer=cell_renderers[col])

        mutations_grid = AgGrid(this_dataset_table,
                                gridOptions=this_gb.build(),
                                reload_data=False,
                                allow_unsafe_jscode=True,
                                height=400,
                                custom_css={"#gridToolBar": {
                                                "padding-bottom": "0px !important",
                                                }
                                            })

    with dotplots:

        this_dataset_table = this_dataset.copy()
        this_dataset_table = this_dataset_table.set_index('Mutation')

        col1, col2 = st.columns(2)

        with col1:
            do_revel = st.checkbox('Show available REVEL classification', )
            revel_co = st.number_input("Cutoff for REVEL score (between 0 and 1)", value=0.5, min_value=0.0, max_value=1.0)
            demask_co = st.number_input("Cutoff for DeMaSk score (absolute value)", value=0.3, min_value=0.0)
            gemme_co = st.number_input("Curoff for GEMME", value=0.3)
        with col2:
            do_demask = st.checkbox('Show available DeMaSk classification', value=True)
            n_muts = st.number_input("Number of mutations per plot", value=50, min_value=0)
            fig_width = st.number_input("Plot width (inches)", value=14, min_value=0)
            fig_height = st.number_input("Plot height (inches)", value=4, min_value=0)

        st.write("""Select up to 50 mutations to be shown in the dot plot below. They should be expressed in plain `[reference amino acid][residue number][mutant amino acid]` format, for example `A49G`.""")

        selected_muts = st.multiselect(label="Mutations to be displayed",
                                       options=this_dataset_table.index,
                                       default=None,
                                       max_selections=50,
                                       placeholder="Type or select one mutation or more")

        if st.button('Generate plot',
                      disabled=len(selected_muts) == 0):
            this_dataset_table = this_dataset_table.loc[selected_muts]
            this_dataset_table.to_csv('does_error.csv')

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
                        fig.savefig(pdf, format='pdf', dpi=300, bbox_inches='tight')

                st.download_button(label="Download as PDF",
                                data=pdf_stream,
                                mime="application/pdf",
                                file_name=f'{protein}-{mode}_dotplots.pdf')

            for fig in plots:
                st.pyplot(fig, use_container_width=False)

    with lolliplots:

        this_dataset_table = this_dataset.copy()
        this_dataset_table = this_dataset_table.set_index('Mutation')
        this_dataset_table = process_df_for_lolliplot(this_dataset_table)

        st.write(f"""Select one or more mutations below, up to 50, to be included
        in the plot. These are only those mutations that are classified as pathogenic
        for AlphaMissense and have an explanation for MAVISp. They are
        {this_dataset_table.shape[0]} in this daataset.""")

        selected_muts = st.multiselect(label="Mutations to be displayed",
                                       options=this_dataset_table.index,
                                       default=None,
                                       placeholder="Type or select one mutation or more",
                                       max_selections=50,
                                       key='sj17h39')

        if st.button('Generate plot',
                     disabled=len(selected_muts) == 0,
                     key='qwe123'):

            this_dataset_table = this_dataset_table.loc[selected_muts]

            plots = plot_lolliplots(this_dataset_table)

            with BytesIO() as pdf_stream:
                with PdfPages(pdf_stream) as pdf:
                    for fig in plots:
                        fig.savefig(pdf, format='pdf', dpi=300, bbox_inches='tight')

                st.download_button(label="Download as PDF",
                                data=pdf_stream,
                                mime="application/pdf",
                                file_name=f'{protein}-{mode}_lolliplots.pdf',
                                key='1231')

            for fig in plots:
                st.pyplot(fig, use_container_width=False)
