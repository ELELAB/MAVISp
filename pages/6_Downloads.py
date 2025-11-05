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
from io import BytesIO
import yaml

st.set_page_config(layout="wide",
    page_title="MAVISp - Datasets and metadata",
    page_icon="🔬")

add_mavisp_logo("static/logo_small.png", image_width='50%')
add_affiliation_logo()

database_dir = get_database_dir()

db_files = find_database_files(database_dir)
db_files_display = db_files.drop(columns=['File name'])

st.title('Download MAVISp dataset')

st.write("""In this page, you are able to download complete MAVISp datasets in bulk. Please select the
dataset of interest by selecting a row of interest in the table below, finally click on the Download button to
start your download.""")

download = st.dataframe(db_files,
                        hide_index=True,
                        use_container_width=True,
                        on_select='rerun',
                        selection_mode='single-row')

selected_row = download['selection']['rows']
if len(selected_row) == 1:
    selected_row = download['selection']['rows'][0]
    download_button_disabled = False
else:
    selected_row = None
    download_button_disabled = True

if not download_button_disabled:
    fh = open(db_files['File name'].loc[selected_row], "rb")
    file_name = f'mavisp_database_{db_files['Date of generation'].loc[selected_row]}.zip'
else:
    fh = BytesIO()
    file_name = 'none.zip'

st.download_button(label="Download dataset",
                    data=fh,
                    mime="application/zip",
                    file_name=file_name,
                    disabled = download_button_disabled)