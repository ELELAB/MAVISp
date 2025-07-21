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
import pandas as pd
from streamlit_utils import *
import py3Dmol
from stmol import showmol
import requests as rq
from requests import HTTPError
from matplotlib.backends.backend_pdf import PdfPages
from io import BytesIO

mode = 'simple_mode'

st.set_page_config(layout="wide",
    page_title="MAVISp simple mode",
    page_icon="📊")

add_mavisp_logo("static/logo_small.png", image_width='50%')

add_affiliation_logo()

database_dir = get_database_dir()

# disable download button for dataframes
# st.markdown(
#     """
#     <style>
#     [data-testid="stElementToolbar"] {
#         display: none;
#     }
#     </style>
#     """,
#     unsafe_allow_html=True
# )

st.title('MAVISp simple mode')

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

show_table = replace_boolean_col(show_table, 'Linker design included')

filter_text = st.text_input("Filter by UniProt AC, gene name or protein RefSeq ID. Press Enter to apply the filter")

if filter_text:
    mask = (
        show_table['Protein'].str.contains(filter_text, case=False, na=False) |
        show_table['Uniprot AC'].str.contains(filter_text, case=False, na=False) |
        show_table['RefSeq ID'].str.contains(filter_text, case=False, na=False)
    )
    filtered_show_table = show_table[mask]
else:
    filtered_show_table = show_table

protein_table = st.dataframe(filtered_show_table,
                                 hide_index=True,
                                 use_container_width=True,
                                 on_select='rerun',
                                 selection_mode='single-row',
                                 column_config = { 'OSF repository for ensemble data' : st.column_config.LinkColumn(display_text='OSF'),
                                                   'GitBook report' : st.column_config.LinkColumn(display_text='report')})

if len(protein_table.selection['rows']) == 0:
    st.stop()

protein_row = protein_table.selection['rows'][0]

protein = show_table.iloc[protein_row]['Protein']
upac = show_table.iloc[protein_row]['Uniprot AC']

st.write(f"Currently viewing: {protein}")

this_dataset = load_dataset(database_dir, protein, mode)

dataset, dotplots, lolliplots, structure = st.tabs(["Dataset", "Classification", "Damaging mutations", "Damaging mutations on structure"])

with dataset:

    data_type = st.radio("Select the data representation to be shown or downloaded",
                        options=["Compact dataset", "Full dataset"],
                        index=0)

    this_dataset_table = this_dataset.copy()
    this_dataset_table = this_dataset_table.fillna(pd.NA)

    if data_type == 'Compact dataset':
        this_dataset_table = get_compact_dataset(this_dataset_table)

    display_dataset_table = this_dataset_table.copy()

    print(display_dataset_table.columns)

    binary_cols = ['EFoldMine - part of early folding region',
                   'is site part of phospho-SLiM',
                   'Mutation predicted to add new phosphorylation site']

    for bcol in binary_cols:
        if bcol in display_dataset_table.columns:
            replace_boolean_col(display_dataset_table, bcol)

    if 'PTMs' in display_dataset_table.columns:
        ptm_link = f"http://www.phosphosite.org/uniprotAccAction?id={upac}"
        display_dataset_table['PTMs'] = display_dataset_table['PTMs'].replace(to_replace='P', value=f'http://www.phosphosite.org/uniprotAccAction?id=P{upac}')

    available_data_sources = list(set(",".join(display_dataset_table['Mutation sources'].tolist()).split(",")))

    selected_sources = st.multiselect(label="Data sources to be considered",
                                        options=available_data_sources,
                                    default=available_data_sources,
                                    placeholder="Type or select a data source or more")

    mask = display_dataset_table['Mutation sources'].apply(lambda x: any(term in x for term in selected_sources))

    filtered_display_dataset_table = display_dataset_table[mask]

    mutation_format = st.radio("Mutation column to filter on", options=['Mutation', 'HGVSp', 'HGVSg'])

    filtering_input_type = st.radio("Input method for mutations to filter on", options=['Select from list', 'Free text format'])

    if filtering_input_type == 'Select from list':

        selected_mutations = st.multiselect(label="Mutations to be considered (all if empty)",
                                            options=filtered_display_dataset_table[mutation_format].dropna().unique().tolist(),
                                            default=None,
                                            placeholder="Type or select a mutation")
        
    else:
        placeholder_text = """Insert one mutation per row according to the selected column, e.g.\n"""
        placeholder_text += f"{'\n'.join(filtered_display_dataset_table[mutation_format].dropna().unique().tolist()[0:3])}\n..."

        selected_mutations_input = st.text_area("Insert list of mutations, one per line, according to the format selected previously",
                     height=100,
                     placeholder=placeholder_text)
        
        print("SMI", selected_mutations_input)
        if selected_mutations_input is None or selected_mutations_input == "":
            selected_mutations = None
        else:
            selected_mutations_input = selected_mutations_input.split('\n')
            selected_mutations_input = list(filter(lambda x: x, selected_mutations_input))
            if len(selected_mutations_input) == len(filtered_display_dataset_table[filtered_display_dataset_table[mutation_format].isin(selected_mutations_input)]):
                selected_mutations = selected_mutations_input
            else:
                st.error("One or more mutations not present in the specified column, please double check their format")
                selected_mutations = None

    if selected_mutations:
        filtered_display_dataset_table = filtered_display_dataset_table[filtered_display_dataset_table[mutation_format].isin(selected_mutations)]

    st.dataframe(filtered_display_dataset_table,
                    hide_index=True,
                    use_container_width=True,
                    column_config = { 'Mutation sources' : st.column_config.ListColumn(),
                                      'PTMs' : st.column_config.LinkColumn(display_text='P')})

    st.download_button(label="Download current dataset view",
                             data=filtered_display_dataset_table.to_csv(),
                             file_name=f'{protein}-{mode}-view.csv',
                             mime="text/csv",
                             key='download-csv-compact')

    with open(os.path.join(database_dir, mode, 'dataset_tables', f'{protein}-{mode}.csv')) as data:
        st.download_button(label="Download original full dataset",
                                 data=data,
                                 file_name=f'{protein}-{mode}.csv',
                                 mime="text/csv",
                                 key='download-csv')

with dotplots:

    this_dataset_table = this_dataset.copy()

    col1, col2 = st.columns(2)

    with col1:
        do_revel = st.checkbox('Show available REVEL classification', )
        revel_co = st.number_input("Cutoff for REVEL score (between 0 and 1)", value=0.5, min_value=0.0, max_value=1.0)
        demask_co = st.number_input("Cutoff for DeMaSk score (absolute value)", value=0.25, min_value=0.0)
        gemme_co = st.number_input("Curoff for GEMME", value=3.0)
    with col2:
        do_demask = st.checkbox('Show available DeMaSk classification', value=True)
        n_muts = st.number_input("Number of mutations per plot", value=50, min_value=0)
        fig_width = st.number_input("Plot width (inches)", value=14, min_value=0)
        fig_height = st.number_input("Plot height (inches)", value=4, min_value=0)

    st.write("""Select up to 50 mutations to be shown in the dot plot below.""")

    mutation_format_dotplot = st.radio("Mutation column to select on", options=['Mutation', 'HGVSp', 'HGVSg'], key="sel_mut_dotplots")

    filtering_input_type_dotplot = st.radio("Input method for mutations to select on", options=['Select from list', 'Free text format'], key="sel_col_dotplot")

    if filtering_input_type_dotplot == 'Select from list':
        selected_mutations_dotplot = st.multiselect(label="Mutations to be considered (all if empty)",
                                                    options=this_dataset_table[mutation_format_dotplot].dropna().unique().tolist(),
                                                    default=None,
                                                    placeholder="Type or select a mutation",
                                                    key="mut_select_dotplot")
        
    else:
        placeholder_text = """Insert one mutation per row according to the selected column, e.g.\n"""
        placeholder_text += f"{'\n'.join(this_dataset_table[mutation_format_dotplot].dropna().unique().tolist()[0:3])}\n..."

        selected_mutations_input_dotplot = st.text_area("Insert list of mutations, one per line, according to the format selected previously",
                     height=100,
                     placeholder=placeholder_text,
                     key="text_area_dotplot")
        
        if selected_mutations_input_dotplot is None or selected_mutations_input_dotplot == "":
            selected_mutations_dotplot = None
        else:
            selected_mutations_input_dotplot = selected_mutations_input_dotplot.split('\n')
            selected_mutations_input_dotplot = list(filter(lambda x: x, selected_mutations_input_dotplot))
            if len(selected_mutations_input_dotplot) == len(this_dataset_table[this_dataset_table[mutation_format_dotplot].isin(selected_mutations_input_dotplot)]):
                selected_mutations_dotplot = selected_mutations_input_dotplot
            else:
                st.error("One or more mutations not present in the specified column, please double check their format")
                selected_mutations_dotplot = None

    if st.button('Generate plot',
                 disabled=not bool(selected_mutations_dotplot)):
        
        this_dataset_table_dotplot = this_dataset_table[this_dataset_table[mutation_format_dotplot].isin(selected_mutations_dotplot)].set_index('Mutation')

        plots = plot_dotplot(this_dataset_table_dotplot,
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

    this_dataset_table_lolliplot = this_dataset.copy()
    this_dataset_table_lolliplot = this_dataset_table_lolliplot.set_index('Mutation')
    this_dataset_table_lolliplot = process_df_for_lolliplot(this_dataset_table_lolliplot)
    this_dataset_table_lolliplot = this_dataset_table_lolliplot.join(this_dataset.set_index('Mutation')[['HGVSp', 'HGVSg']])
    this_dataset_table_lolliplot = this_dataset_table_lolliplot.reset_index()

    if this_dataset_table_lolliplot.shape[0] == 0:
        st.write("""There are no suitable mutations in this dataset, therefore
        this section has been disabled""")
        st.stop()

    st.write(f"""Select one or more mutations below, up to 50, to be included
    in the plot. These are only those mutations that are at the same time
    i) classified as pathogenic for AlphaMissense, ii) classified as loss
    of function or gain of function for either GEMME or DeMaSk and
    iii) damaging for the respective module in MAVISp. They are
    {this_dataset_table.shape[0]} in this dataset.""")

    mutation_format_lolliplot = st.radio("Mutation column to select on", options=['Mutation', 'HGVSp', 'HGVSg'], key="sel_mut_lolliplots")

    filtering_input_type_lolliplot = st.radio("Input method for mutations to select on", options=['Select from list', 'Free text format'], key="sel_col_lolliplots")

    if filtering_input_type_lolliplot == 'Select from list':
        selected_mutations_lolliplot = st.multiselect(label="Mutations to be considered (all if empty)",
                                                    options=this_dataset_table[mutation_format_lolliplot].dropna().unique().tolist(),
                                                    default=None,
                                                    placeholder="Type or select a mutation",
                                                    key="mut_select_lolliplots")
        
    else:
        placeholder_text = """Insert one mutation per row according to the selected column, e.g.\n"""
        placeholder_text += f"{'\n'.join(this_dataset_table[mutation_format_lolliplot].dropna().unique().tolist()[0:3])}\n..."

        selected_mutations_input_lolliplot = st.text_area("Insert list of mutations, one per line, according to the format selected previously",
                     height=100,
                     placeholder=placeholder_text,
                     key="text_area_lolliplot")
        
        if selected_mutations_input_lolliplot is None or selected_mutations_input_lolliplot == "":
            selected_mutations_lolliplot = None
        else:
            selected_mutations_input_lolliplot = selected_mutations_input_lolliplot.split('\n')
            selected_mutations_input_lolliplot = list(filter(lambda x: x, selected_mutations_input_lolliplot))
            if len(selected_mutations_input_lolliplot) == len(this_dataset_table_lolliplot[this_dataset_table_lolliplot[mutation_format_lolliplot].isin(selected_mutations_input_lolliplot)]):
                selected_mutations_lolliplot = selected_mutations_input_lolliplot
            else:
                st.error("One or more mutations not present in the specified column, please double check their format")
                selected_mutations_lolliplot = None

    if st.button('Generate plot',
                    disabled=not bool(selected_mutations_lolliplot),
                    key='qwe123'):
        print(this_dataset_table_lolliplot.columns)
        print(this_dataset_table_lolliplot.index)

        this_dataset_table_lolliplot = this_dataset_table_lolliplot[this_dataset_table_lolliplot[mutation_format_lolliplot].isin(selected_mutations_lolliplot)].set_index('Mutation')

        plots = plot_lolliplots(this_dataset_table_lolliplot)

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

    # download model and stop if it can't be found
    if not disable_structure:
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
