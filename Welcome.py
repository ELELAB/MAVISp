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
#from streamlit_extras.app_logo import add_logo
import pandas as pd
from streamlit_utils import *
import streamlit as st



st.set_page_config(layout="wide",
    page_title="Hello",
    page_icon="👋")

add_logo("assets/logo_small.png", image_width='80%')

st.header("Welcome to MAVISp!")

st.write("""MAVISp includes different types of predictions for the effect of mutations on protein function, structure or more. Its final goal is to predict the effect of relevant variants of unknown significance and identify their mechanism of action.

MAVISp has been designed by the Cancer Structural/Systems Biology group, headed by Elena Papaleo, with shared affiliation at the Danish Cancer Research Center and the Danish Technical University, Department of Health and Technology. 

Please use the menu on the left to navigate the website.

If you use data from MAVISp from your research, please cite our preprint:""")

st.code("""MAVISp: Multi-layered Assessment of VarIants by Structure for proteins
Matteo Arnaudi, Ludovica Beltrame, Kristine Degn, Mattia Utichi, Pablo Sánchez-Izquierdo,
Simone Scrima, Francesca Maselli, Karolina Krzesińska, Terézia Dorčaková, Jordan Safer,
Katrine Meldgård, Julie Bruun Brockhoff, Amalie Drud Nielsen, Alberto Pettenella, Jérémy Vinhas,
Peter Wad Sackett, Claudia Cava, Sumaiya Iqbal, Matteo Lambrughi, Matteo Tiberti, Elena Papaleo
bioRxiv 2022.10.22.51332; doi: https://doi.org/10.1101/2022.10.22.513328""", language=None)

