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
from fsspec.implementations.dirfs import DirFileSystem
from fsspec.implementations.zip   import ZipFileSystem
from pathlib import Path
from dot_plot import plot as do_dotplots
from dot_plot import process_input as process_input_for_dotplot
from dot_plot import generate_summary, filter_am_summary
from lolliplot import process_input as process_input_for_lolliplot
from lolliplot import plot as do_lolliplot

@st.cache_data
def get_base64_of_bin_file(png_file):
    with open(png_file, "rb") as f:
        data = f.read()
    return base64.b64encode(data).decode()

def build_markup_for_logo(
    png_file,
    background_position="center top",
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
def get_database_filesystem(dir_var_name='MAVISP_DATABASE_PATH',
                            db_var_name='MAVISP_DATABASE_NAME',
                            default_dir_name='.',
                            default_db_name='database'):
    
    dir_name = os.getenv(dir_var_name)
    db_name = os.getenv(db_var_name)

    if dir_name is None:
        dir_name = default_dir_name

    if db_name is None:
        db_name = default_db_name

    db_path = Path(dir_name) / Path(db_name)

    print('MAVISP_DATABASE_PATH', os.getenv(dir_var_name))
    print('MAVISP_DATABASE_NAME', os.getenv(db_var_name))
    print(dir_name, '|', db_name, '|', db_path)

    if not db_path.exists():
        raise FileExistsError(f"provided database path {db_path} does not exist")

    if db_path.is_file() and db_path.suffix == ".zip":
        fs = ZipFileSystem(db_path)
    elif db_path.is_dir():
        fs = DirFileSystem(db_path)
    else:
        raise TypeError(f"database must be either a directory or a zip file with .zip extensions. Current database is: {db_path}")

    return fs

def add_affiliation_logo():
    columns = st.sidebar.columns(2)

    with columns[0]:
        st.write("""<div style="width:100%;text-align:center;"><a href="https://www.cancer.dk" style="float:center"><img src="app/static/dcs_logo.png" width="60px"></img></a></div>""", unsafe_allow_html=True)

    with columns[1]:
        st.write("""<div style="width:100%;text-align:center;"><a href="https://www.dtu.dk" style="float:center"><img src="app/static/dtu_logo.png" width="60px"></img></a></div>""", unsafe_allow_html=True)

@st.cache_data
def load_dataset(_data_fs, protein, mode):
    with _data_fs.open(os.path.join(mode, 'dataset_tables', f'{protein}-{mode}.csv')) as fh:
        return pd.read_csv(fh)


@st.cache_data
def load_main_table(_data_fs, mode):
    with _data_fs.open(os.path.join(mode, 'index.csv')) as fh:
        return pd.read_csv(fh).sort_values('Protein')

@st.cache_data
def load_clinvar_dict(tsv_file):
    clinvar_dict = pd.read_csv(tsv_file,
                                sep='\t',
                                header=None,
                                names=['clinvar', 'internal_category'])
    return clinvar_dict.set_index('clinvar')['internal_category'].to_dict()

@st.cache_data
def plot_dotplot(df, demask_co, revel_co, gemme_co, fig_width=14, fig_height=4, n_muts=50, do_revel=False, do_demask=True):

    df = df.copy()

    if 'ClinVar Interpretation' not in df.columns:
        df['ClinVar Interpretation'] = None

    clinvar_dict = load_clinvar_dict('mavisp/data/clinvar_interpretation_internal_dictionary.txt')

    processed_df, full_df, clinvar_mapped_df = process_input_for_dotplot(df,
                                                            d_cutoff=demask_co,
                                                            r_cutoff=revel_co,
                                                            g_cutoff=gemme_co,
                                                            residues=None,
                                                            mutations=None,
                                                            clinvar_dict=clinvar_dict,
                                                            plot_Revel=True,
                                                            plot_Demask=True,
                                                            plot_Source=None,
                                                            plot_Clinvar=None,
                                                            color_Clinvar=True)

    if not do_revel:
        processed_df = processed_df.drop(columns=['REVEL'])

    my_plots = do_dotplots(processed_df, clinvar_mapped_df, fig_width, fig_height, n_muts, False, True)

    return my_plots

@st.cache_data
def process_df_for_lolliplot(df):
    df = df.copy()

    clinvar_dict = load_clinvar_dict('mavisp/data/clinvar_interpretation_internal_dictionary.txt')

    processed_df, full_df, clinvar_mapped_df = process_input_for_dotplot(df,
                                                            r_cutoff=0.5,
                                                            d_cutoff=0.25,
                                                            g_cutoff=3.0,
                                                            residues=None,
                                                            mutations=None,
                                                            clinvar_dict=clinvar_dict,
                                                            plot_Revel=False,
                                                            plot_Demask=True,
                                                            plot_Source=None,
                                                            plot_Clinvar=None,
                                                            color_Clinvar=False)

    text, summary_df = generate_summary(full_df, d_cutoff=0.25, r_cutoff=0.5)

    filtered_summary_df = filter_am_summary(summary_df, processed_df, True)

    return process_input_for_lolliplot(filtered_summary_df)

@st.cache_data
def plot_lolliplots(df, muts_per_plot=50):
    return do_lolliplot(df, muts_per_plot)

@st.cache_data
def get_compact_dataset(this_dataset_table):
    default_cols = ['Mutation', 'HGVSp', 'HGVSg', 'Mutation sources']
    selected_cols = [c for c in this_dataset_table.columns if "classification" in c ]
    return this_dataset_table[default_cols + selected_cols + ['References']]

def replace_boolean_col(df, col, dictionary={True : 'Yes', False : 'No'}):

    df[col] = df[col].astype(str)

    for k,v in dictionary.items():
        k, v = str(k), str(v)
        df[col] = df[col].replace(to_replace=k, value=v)
    return df
