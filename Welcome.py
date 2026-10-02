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

from streamlit_utils import *
import streamlit as st

st.set_page_config(layout="wide",
    page_title="MAVISp - Welcome",
    page_icon="👋")

add_mavisp_logo("static/logo_small.png", image_width='50%')

add_affiliation_logo()

st.header("Welcome to MAVISp!")

st.write("""MAVISp includes different types of predictions for the effect of mutations
on protein function, structure or more. Its final goal is to predict the effect of
relevant variants of unknown significance and identify their mechanism of action.

MAVISp has been designed by the Cancer Structural Biology (Danish Cancer Institute, Denmark)
and the Cancer Systems Biology (Department of Health and Technology at the Technical University of Demark)
groups, headed by Dr. Elena Papaleo

Please use the menu on the left to navigate the website.

If you use data from MAVISp from your research, please cite
[our publication](https://onlinelibrary.wiley.com/doi/10.1002/pro.70548):""")

st.code("""Arnaudi M, Utichi M, Degn K, Tiberti M, Beltrame L, Krzesińska K, et al.
MAVISp: A modular structure-based framework for protein variant effects.
Protein Science. 2026;35(5):e70548. https://doi.org/10.1002/pro.70548""", language=None)

st.write("""Please see the "Acknowledgement and data usage" section for information about our data
sources, data license terms and data reuse""")

