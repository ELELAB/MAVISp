# MAVISp - various utilities for MAVISp web server
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
import base64
import os

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

def add_logo(png_file, *args, **kwargs):
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
    return var_name
