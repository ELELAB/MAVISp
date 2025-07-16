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
from streamlit_utils import *
import yaml

st.set_page_config(layout="wide",
    page_title="Datasets and metadata",
    page_icon="📊")

add_mavisp_logo("static/logo_small.png", image_width='50%')
add_affiliation_logo()

database_dir = get_database_dir()

st.title('Datasets and metadata')

st.write('''This page is a repository for the metadata about the calculation
and data included in MAVISp. Metadata is meant to display the methodology behind
each module as well as possible data sources and their references.

Please select a mode, gene, and module below - the corresponding metadata will
be displayed.

Currently, we only include metadata for the `EXPERIMENTAL_DATA` module, and we plan
on adding more in following releases''')

st.write('''Please see the Acknowledgement and data usage page for information on our data sources, licensing term, and data reuse permissions''')

# Mode selection
mode = st.selectbox(
    "Select MAVISp mode", 
    options=['simple_mode', 'ensemble_mode'], 
    index=None, 
    key='mode',
)

# Gene selection
gene_options = []
if mode is not None:
    try:
        main_table = load_main_table(database_dir, mode)
        gene_options = main_table['Protein'].to_list()
    except FileNotFoundError:
        st.write(f"No entries are currently available for {mode}.")
        st.stop()

gene = st.selectbox(
    "Select gene", 
    options=gene_options, 
    index=None, 
    key='gene',
)

# Module selection
modules_options = []
if gene is not None:
    try:
        with open(f'{database_dir}/{mode}/metadata/{gene}.yaml') as fh:
            st.session_state.data = yaml.safe_load(fh)
        modules_options = [k for k in st.session_state.data.keys() if st.session_state.data[k] is not None]
    except FileNotFoundError:
        st.write(f"Metadata file for {mode}, {gene} was not found on the database. Please contact the website administrators")
        st.stop()
    except Exception:
        st.write(f"There was a problem loading metadata for {mode}, {gene}. Please contact the website administrators")
        st.stop()

module = st.selectbox(
    "Select module", 
    options=modules_options, 
    index=None,
    key='module'
)

# Display data only if all three selections are made
if mode is not None and gene is not None and module is not None and st.session_state.data is not None:
    st.json(st.session_state.data[module])
