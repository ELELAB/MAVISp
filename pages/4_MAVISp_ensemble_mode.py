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
import py3Dmol
from stmol import showmol
import requests as rq
from requests import HTTPError
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from io import BytesIO

mode = 'ensemble_mode'

st.set_page_config(layout="wide",
    page_title="MAVISp ensemble mode",
    page_icon="📊")

add_mavisp_logo("static/logo_small.png", image_width='50%')

add_affiliation_logo()

database_dir = get_database_dir()

st.title("MAVISp ensemble mode")

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

datasets_grid = AgGrid(show_table, enable_enterprise_modules=False,
                      gridOptions=gb_datasets_grid.build(),
                      update_mode=GridUpdateMode.SELECTION_CHANGED,
                      fit_columns_on_grid_load = True,
                      reload_data=False,
                      allow_unsafe_jscode=True,
                      height=200,
                      key="id_row",
                      custom_css={"#gridToolBar": {
                                      "padding-bottom": "0px !important",
                                      }
                                  })

if datasets_grid["selected_rows"] is not None and len(datasets_grid["selected_rows"]) == 1:

    protein = datasets_grid["selected_rows"].iloc[0]['Protein']
    upac = datasets_grid["selected_rows"].iloc[0]['Uniprot AC']

    st.write(f"Currently viewing: {protein}")

    this_dataset = load_dataset(database_dir, protein, mode)
    dataset, dotplots, lolliplots, structure = st.tabs(["Dataset", "Classification", "Damaging mutations", "Damaging mutations on structure"])

    ptm_columns = [c for c in this_dataset.columns if c.startswith('PTMs')]

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
        if 'Mutation sources' in this_dataset_table.columns:
            this_gb.configure_column('Mutation sources', cellRenderer=cell_renderers['Mutation sources'])
        for col in ptm_columns:
            if col in this_dataset_table.columns:
                this_gb.configure_column(col, cellRenderer=cell_renderers['PTMs'])

        mutations_grid = AgGrid(this_dataset, enable_enterprise_modules=False,
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
            demask_co = st.number_input("Cutoff for DeMaSk score (absolute value)", value=0.25, min_value=0.0)
            gemme_co = st.number_input("Cutoff for GEMME", value=3.0)
        with col2:
            do_demask = st.checkbox('Show available DeMaSk classification', value=True)
            n_muts = st.number_input("Number of mutations per plot", value=50, min_value=0)
            fig_width = st.number_input("Plot width (inches)", value=14, min_value=0)
            fig_height = st.number_input("Plot height (inches)", value=6, min_value=0)

        st.write("""Select up to 50 mutations to be shown in the dot plot below. They should be expressed in plain `[reference amino acid][residue number][mutant amino acid]` format, for example `A49G`.""")

        selected_muts = st.multiselect(label="Mutations to be displayed",
                                       options=this_dataset_table.index,
                                       default=None,
                                       max_selections=50,
                                       placeholder="Type or select one mutation or more")

        if st.button('Generate plot',
                      disabled=len(selected_muts) == 0):
            this_dataset_table = this_dataset_table.loc[selected_muts]

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
        in the plot. These are only those mutations that are at the same time
        i) classified as pathogenic for AlphaMissense, ii) classified as loss
        of function or gain of function for either GEMME or DeMaSk and
        iii) damaging for the respective module in MAVISp. They are
        {this_dataset_table.shape[0]} in this dataset.""")

        disable_lolliplot = False
        if this_dataset_table.shape[0] == 0:
            st.write("""There are no suitable mutations in this dataset, therefore
            this section has been disabled""")
            disable_lolliplot = True

        selected_muts = st.multiselect(label="Mutations to be displayed",
                                       options=this_dataset_table.index,
                                       default=None,
                                       placeholder="Type or select one mutation or more",
                                       max_selections=50,
                                       key='sj17h39')

        if st.button('Generate plot',
                     disabled=len(selected_muts) == 0 or disable_lolliplot,
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
                                file_name=f'{protein}-{mode}_lolliplots.pdf')

            for fig in plots:
                st.pyplot(fig, use_container_width=False)

    with structure:

        structure_colors = {'Stability'  : 'redCarbon',
                            'Local Int.' : 'orangeCarbon',
                            'PTM'        : 'greenCarbon',
                            'Long Range' : 'blueCarbon',
                            'Multiple'   : 'purpleCarbon'}

        this_dataset_table = this_dataset.copy()
        this_dataset_table = this_dataset_table.set_index('Mutation')
        this_dataset_table = process_df_for_lolliplot(this_dataset_table)

        st.write("""This tab displays the AlphaFold model for the selected protein,
        if available. The checkbox below will activate colouring of mutations that
        are at the same time i) classified as pathogenic for AlphaMissense, ii) classified as loss
        of function or gain of function for either GEMME or DeMaSk and
        iii) damaging for the respective module in MAVISp. You
        can have multiple checkboxes active at the same time; residues with mutations
        that have multiple effects for MAVISp will be coloured in purple.""")

        st.write("""We offer two types of analysis: one that colors only those residues
        for which the number of damaging mutations is higher than a user-selected
        threshold. This can be useful to spot mutational hotspots. In the second,
        the user can choose to color selected residues of interest which will be
        displayed on the structure.""")

        disable_structure = False
        if this_dataset_table.shape[0] == 0:
            st.write("""There are no suitable mutations in this dataset, therefore
            this tab has been disabled""")
            disable_structure = True

        if not disable_structure:

            # download model and stop if it can't be found
            try:
                response = rq.get(f"https://alphafold.ebi.ac.uk/files/AF-{upac}-F1-model_v4.pdb")
                response.raise_for_status()
            except ConnectionError:
                st.write("Failed connecting to the AlphaFold Protein Structure Database")
                st.stop()
            except HTTPError:
                st.write("Could not fetch protein structure model from the AlphaFold Protein Structure Database")
                st.stop()
            else:
                model = response.text

            # set up viewer
            viewer = py3Dmol.view(width=900, height=600)
            viewer.addModel(model, 'pdb')
            viewer.setStyle({ "cartoon": { "color" : "lightgray", "style" : "parabola" } })

            # decide which classification terms to consider
            st.write("Classification terms to be considered (choose one or more):")

            col_structure1, col_structure2 = st.columns(2)

            with col_structure1:
                structure_stability = st.checkbox('Stability')
                structure_li        = st.checkbox('Local Int.')
            with col_structure2:
                structure_ptm       = st.checkbox('PTM')
                structure_lr        = st.checkbox('Long range')

            interesting_cols = []
            for col, checkbox in [('Stability',  structure_stability),
                                  ('Local Int.', structure_li),
                                  ('PTM',        structure_ptm),
                                  ('Long Range', structure_lr)]:
                if checkbox:
                    interesting_cols.append(col)

            # pre-process the dataframe
            this_dataset_table = this_dataset_table[interesting_cols]
            this_dataset_table = this_dataset_table.loc[this_dataset_table[interesting_cols].sum(axis=1) > 0]

            this_dataset_table['residue'] = this_dataset_table.index.str[1:-1]
            tmp_df1 = this_dataset_table.reset_index().groupby('residue').agg({'Mutation':lambda x: " ".join(x.tolist())})
            tmp_df2 = this_dataset_table.groupby('residue').agg(sum)
            this_dataset_table = tmp_df1.join(tmp_df2)

            # select analysis type and act accordingly
            analysis_type = st.radio("Type of analysis", options=['Hotspots', 'Custom sites'])

            if analysis_type == 'Hotspots':
                min_muts = st.slider(label="Minimum number of damaging mutations",
                                     min_value=1, max_value=19, step=1,
                                     value=5)
                this_dataset_table = this_dataset_table.loc[this_dataset_table['Mutation'].str.split(' ').apply(len) >= min_muts]
            elif analysis_type == 'Custom sites':
                selected_muts = st.multiselect(label="Sites of interest",
                                            options=this_dataset_table.index,
                                            default=None,
                                            max_selections=50,
                                            placeholder="Type or select one site or more")
                this_dataset_table = this_dataset_table.loc[selected_muts]

            labels = st.radio("Show residue labels:", options=['none', 'for mutations', 'for sites'])
            labels_in_front_checkbox = st.checkbox("Labels are always in front of the structure")

            if labels_in_front_checkbox:
                front_labels = 'true'
            else:
                front_labels = 'false'

            # stop unless at least one classification has been selected
            if len(interesting_cols) == 0:
                st.stop()

            st.markdown('''
            Color legend:

            🔴 Stability
            🟠 Local interactions
            🟢 PTM
            🔵 Long Range
            🟣 Mutations with multiple classifications
            ''')

            # prepare colors
            this_dataset_table['color'] = ''
            for col in interesting_cols:
                this_dataset_table.loc[this_dataset_table[col] > 0, 'color'] = structure_colors[col]

            this_dataset_table.loc[(this_dataset_table[interesting_cols] > 0).sum(axis=1) > 1, 'color'] = structure_colors['Multiple']
            for color in this_dataset_table['color'].unique():
                residues = this_dataset_table.loc[this_dataset_table['color'] == color].index.tolist()
                viewer.addStyle({'resi': residues}, {'cartoon': {'colorscheme': color}})

            if labels == 'for mutations':
                for idx, row in this_dataset_table.iterrows():
                    viewer.addLabel(row['Mutation'], {'fontColor':'black', 'backgroundColor':'lightgray'}, {'resi':idx})

            elif labels == 'for sites':
                for idx, row in this_dataset_table.iterrows():
                    viewer.addLabel(idx, {'fontColor':'black', 'backgroundColor':'lightgray', 'inFront':front_labels}, {'resi':idx})

            viewer.zoomTo()
            showmol(viewer, width=900, height=600)
