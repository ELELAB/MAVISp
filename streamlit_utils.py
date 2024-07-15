# MAVISp - various utilities for MAVISp web server
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
import base64
import os
import pandas as pd
from st_aggrid import JsCode
from dot_plot import plot as do_dotplots
from dot_plot import process_input as process_input_for_dotplot
from dot_plot import generate_summary
from lolliplot import process_input as process_input_for_lolliplot
from lolliplot import plot as do_lolliplot

@st.cache_data
def get_base64_of_bin_file(png_file):
    with open(png_file, "rb") as f:
        data = f.read()
    return base64.b64encode(data).decode()

def build_markup_for_logo(
    png_file,
    background_position="50% 10%",
    margin_top="0%",
    margin_bottom="10%",
    image_width="60%",
    image_height="",
):
    binary_string = get_base64_of_bin_file(png_file)
    return """
            <style>

                [data-testid="stSidebarNav"] {
                    background-image: url("data:image/png;base64,%s");
                    background-repeat: no-repeat;
                    background-position: %s;
                    margin-top: %s;
                    margin-bottom: %s;
                    background-size: %s %s;
                }
            [data-testid="stSidebarNav"]::before {
                content: "  ";
                margin-left: 20px;
                margin-top: 50px;
                margin-bottom: 50px;
                font-size: 30px;
                position: relative;
                top: 100px;
            }
        </style>
            """ % (
        binary_string,
        background_position,
        margin_top,
        margin_bottom,
        image_width,
        image_height,
    )

def add_mavisp_logo(png_file, *args, **kwargs):
    logo_markup = build_markup_for_logo(png_file, *args, **kwargs)
    st.markdown(
        logo_markup,
        unsafe_allow_html=True,
    )

@st.cache_data
def get_database_dir(var_name='MAVISP_DATABASE', default_dir='./database'):
    dir_name = os.getenv(var_name)
    if dir_name is None:
        return default_dir
    return dir_name

def add_affiliation_logo():
    columns = st.sidebar.columns(2)

    with columns[0]:
        st.write("""<div style="width:100%;text-align:center;"><a href="https://www.cancer.dk" style="float:center"><img src="app/static/dcs_logo.png" width="60px"></img></a></div>""", unsafe_allow_html=True)

    with columns[1]:
        st.write("""<div style="width:100%;text-align:center;"><a href="https://www.dtu.dk" style="float:center"><img src="app/static/dtu_logo.png" width="60px"></img></a></div>""", unsafe_allow_html=True)

@st.cache_data
def load_dataset(data_dir, protein, mode):
    return pd.read_csv(os.path.join(data_dir, mode, 'dataset_tables', f'{protein}-{mode}.csv'))

@st.cache_data
def load_main_table(data_dir, mode):
    return pd.read_csv(os.path.join(data_dir, mode, 'index.csv')).sort_values('Protein')

@st.cache_data
def plot_dotplot(df, demask_co, revel_co, gemme_co, fig_width=14, fig_height=4, n_muts=50, do_revel=False, do_demask=True):
    df = df.copy()
    processed_df, _ = process_input_for_dotplot(df, d_cutoff=demask_co, r_cutoff=revel_co, g_cutoff=gemme_co)

    if not do_revel:
        processed_df = processed_df.drop(columns=['REVEL'])

    my_plots = do_dotplots(processed_df, fig_width, fig_height, n_muts, do_demask)
    return my_plots


@st.cache_data
def process_df_for_lolliplot(df):
    _, processed_df = process_input_for_dotplot(df, d_cutoff=0.3, r_cutoff=0.5, g_cutoff=0.3)
    summary_df = generate_summary(processed_df)
    return process_input_for_lolliplot(summary_df)

@st.cache_data
def plot_lolliplots(df, muts_per_plot=50):
    return do_lolliplot(df, muts_per_plot)

@st.cache_data
def get_compact_dataset(this_dataset_table):
    default_cols = ['Mutation', 'HGVSp', 'HGVSg', 'Mutation sources']
    selected_cols = [c for c in this_dataset_table.columns if "classification" in c ]
    return this_dataset_table[default_cols + selected_cols + ['References']]

# JavaScript column renderers, to dynamically add web links

cell_renderers = {}

cell_renderers['GitBook report'] = JsCode('''
class PTMCellRenderer {
    init(params) {
        this.eGui = document.createElement('span');
        if (params.value == null) {
            this.eGui.innerHTML = '';
        } else {
            this.eGui.innerHTML = `<a target="_parent" href="${params.value}">report</a>`;
        }
    }
  getGui() {
    return this.eGui;
  }
}
''')

cell_renderers['OSF repository for ensemble data'] = JsCode('''
class PTMCellRenderer {
    init(params) {
        this.eGui = document.createElement('span');
        if (params.value == null) {
            this.eGui.innerHTML = '';
        } else {
            this.eGui.innerHTML = `<a target="_parent" href="${params.value}">OSF</a>`;
        }
    }
  getGui() {
    return this.eGui;
  }
}
''')

cell_renderers['Mutation sources'] = JsCode('''
class SourceCellRenderer {
  init(params) {
    this.eGui = document.createElement('span');
    // Split the string into an array using comma as separator
    var array = params.value.split(',');

    // Loop through array and create a link for each item
    for (var i = 0; i < array.length; i++) {
        var item = array[i].trim();
        var link = "";

        // Check the item and assign the corresponding link
        if (item === "COSMIC") {
            link = "https://cancer.sanger.ac.uk/cosmic";
        } else if (item === "cBioPortal") {
            link = "https://www.cbioportal.org";
        } else if (item === "clinvar") {
            item = "ClinVar";
            link = "https://www.ncbi.nlm.nih.gov/clinvar/";
        }

        // If a link is assigned, create the anchor tag
        if (link !== "") {
            array[i] = `<a target="_parent" href=${link}>${item}</a>`;
        } else {
            array[i] = item;
        }
    }

    this.eGui.innerHTML = array.join(', ');

  }
  getGui() {
    return this.eGui;
  }
}''')

cell_renderers['PTMs'] = JsCode('''
class PTMCellRenderer {
    init(params) {
        this.eGui = document.createElement('span');
        if (params.value == null) {
            this.eGui.innerHTML = '';
        } else {
            this.eGui.innerHTML = `<a target="_parent" href="http://www.phosphosite.org/uniprotAccAction?id=${params.data.UniProtAC}">${params.value}</a>`;
        }
    }
  getGui() {
    return this.eGui;
  }
}''')
