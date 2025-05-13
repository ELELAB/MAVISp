# MAVISp - News page
# Copyright (C) 2024 Matteo Tiberti, Danish Cancer Society
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
import streamlit as st

st.set_page_config(layout="wide",
    page_title="Events",
    page_icon="🤝")

add_mavisp_logo("static/logo_small.png", image_width='50%')

add_affiliation_logo()

st.title("Events")

st.header("Planned events")

st.subheader("MAVISp BioCurator Training Workshop 2025")

st.markdown("""
- When: 3rd to 5th of September 2025
- Where: fully online
- Participation fee: none (free of charge)
- Deadline for registering: June 27th 2025
- Requirements:
  - Basic competence with Linux and its terminal
  - Basic knowledge of structural and computational biology
  - Reading: MAVISp preprint and MAVISp GitLab site
- Registration [here](https://docs.google.com/forms/d/e/1FAIpQLScQmPatyYt43JwA6wyQ2V4Pyh7nLVo0uWa9kAAk3kyZawvSlg/viewform?usp=dialog)

We are organizing the first **MAVISp BioCurator Training Workshop**! 

During the event, you will learn the basics behind the MAVISp methodology and work towards becoming a curator for the MAVISp server and website, on your own proteins of interest. The program includes:
- An overview of the methodological framework behind MAVISp
- Hands on sessions, to learn how to use the tools supported by MAVISp to curate your own protein of choice
- Access to a computing server to perform the analyses

The methodology you will learn can be used to generate data for novel proteins to be added to the MAVISp dataset, as well as standalone for your own research.

Read more about MAVISp on [our preprint](https://www.biorxiv.org/content/10.1101/2022.10.22.513328v6) and the [MAVISp gitbook](https://elelab.gitbook.io/mavisp/documentation/how-to-contribute-as-a-biocurator)

Please register by filling out [this form](https://docs.google.com/forms/d/e/1FAIpQLScQmPatyYt43JwA6wyQ2V4Pyh7nLVo0uWa9kAAk3kyZawvSlg/viewform?usp=dialog)
""" 
)