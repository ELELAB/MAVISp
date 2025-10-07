# MAVISp - News page
# Copyright (C) 2025 Matteo Tiberti, Danish Cancer Society
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

from streamlit_utils import *
import streamlit as st

st.set_page_config(layout="wide",
    page_title="MAVISp - News",
    page_icon="🗞️")

add_mavisp_logo("static/logo_small.png", image_width='50%')

add_affiliation_logo()

st.title("News and updates")

st.header("In evidence")

st.subheader("26/05/2025 - MAVISp BioCurator Training Workshop 2025 - deadline extended and programme available")
st.badge("Training Event", icon="📖", color="green")
st.text("""We have extended the deadline for joining the MAVISp BioCurator Training Workshop 2025 to July 4th. We have also made available a tentative programme for the workshop. For more information, see the Events page""")
st.divider()

st.subheader("15/05/2025 - MAVISp BioCurator Training Workshop 2025 announced")
st.badge("Training Event", icon="📖", color="green")
st.text("""We are organizing the first MAVISp BioCurator Training Workshop, to be held in early September 2025, fully online. For more information, see the Events page""")
st.divider()

st.header("Others")

st.subheader("26/05/2025 - New implementation of the LONG_RANGE module")
st.badge("Training Event", icon="📖", color="green")
st.markdown("""Following a benchmarking effort that we recently did to improve the prediction performance of the LONG_RANGE module, we have devised a new and improved LONG_RANGE workflow to be used for new MAVISp entries. At the time being, the MAVISp database contains both entries done using this new protocol as well as entries made using the previous one. Which of them was used for any entry can be discerned in the main proteins table of the MAVISp website, as the new workflow uses simple interatomics contacts for filtering, while the older one used inter-carbon alpha distances. For more details about the new protocol see our [GitHub repository](https://github.com/ELELAB/MAVISp_allosigma2_workflow/), and you'll be able to read all the relevant details in our soon to be published paper.""")
st.divider()
