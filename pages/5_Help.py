# MAVISp - Streamlit application
# Copyright (C) 2023 Matteo Tiberti, Danish Cancer Society
#               2023 Elena Papaleo, Danish Cancer Society 
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
#
# Text in this page was reviewed by Elena Papaleo (Danish Cancer Society,
# Danish Technical University) on 2023-11-20

import streamlit as st

st.set_page_config(layout="wide",
    page_title="Help",
    page_icon="❓")

st.title('Usage of the web app')

st.write("""This section contains instructions and directions on how to make
the most of the MAVISp web app. For lengthier documentation on the collected
data and its interpretation please see the Documentation section""")

st.header("Introduction")

st.write("""The MAVISp web app has been designed to make the MAVISp database
easier to access and explore. It was devised using Streamlit, a web app 
framework for Python, and its source code is available on the MAVISp GitHub 
repository. If you have any suggetion or complaint, please feel free to open an 
issue in the Issues section of the repository.

Navigating the website is as simple as selecting one of the links on the 
panel on the left. These are:

  - Welcome: a welcome page
  - MAVISp simple mode: this page contains all the data that we collected in MAVISp using the simple mode in MAVISp, as detailed in the MAVISp paper. Briefly, this includes data collected working on a single protein structure, most often a model from the AlphaFold Protein Structure Database
  - MAVISp ensemble mode: similarly as the previous page, this contains the data collected using the ensemble mode of MAVISp, i.e. data collected on structural ensembles, most often (but not exclusively) from atomistic simulations
  - Acknowledgements and data usage: this page contains acknowledgements on our data sources, used software, and instructions, information on data reuse and how to cite our work.
  - Documentation: (to be added)
  - Help: this page""")

st.header("MAVISp simple mode and MAVISp ensemble mode")

st.subheader("First look at the dataset")

st.write("""These two pages have very similar layout and content and they are designed for users to download the content of the database and explore its content.

Users are first shown a list of proteins to choose from, by clicking on the 
checkbox at the beginning of each row in the table. Each row also contains a 
link to the corresponding report, hosted on GitBook. In the case of ensemble 
mode, the table also includes a link to the Open Science Foundation repository 
where the structural ensemble is stored.

Once a protein of interest has been selected, users can choose betwen three tabs to explore and visualize the data: XX, YY and ZZ.
""")

st.subheader("First exploration and download of the dataset")

st.write("""The Dataset tab contains, a table with the MAVISp dataset for the selected protein. Users can select between two visualization formats for the dataset: Compact and full. Both dataset contains one protein-level mutation per row. The compact dataset is a more streamlined version of the full dataset, containing only the most important outcomes of the MAVISp framework (mostly, the effect classification for each module), while the Full dataset contains all the availbile data in MAVISp. These datasets only differ in the number of columns and contain the same mutations.

By using the button on each column header, which becomes visibile on when hovering the mouse point over the header itself, users can further explore and refine the dataset, by for instance adding or removing columns, hiding or showing specific rows. Clicking on any table header allows to sort the table by that column, in either ascending or descending order.""")