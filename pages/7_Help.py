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
from streamlit_utils import add_mavisp_logo, add_affiliation_logo

st.set_page_config(layout="wide",
    page_title="Help",
    page_icon="❓")

add_mavisp_logo("static/logo_small.png", image_width='50%')

add_affiliation_logo()

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
  - MAVISp simple mode: this page contains all the data that we collected
  in MAVISp using the simple mode in MAVISp, as detailed in the MAVISp paper.
  Briefly, this includes data collected working on a single protein structure,
  most often a model from the AlphaFold Protein Structure Database
  - MAVISp ensemble mode: similarly as the previous page, this contains the
  data collected using the ensemble mode of MAVISp, i.e. data collected on
  structural ensembles, most often (but not exclusively) from atomistic
  simulations
  - Acknowledgements and data usage: this page contains acknowledgements
  on our data sources, used software, and instructions, information on
  data reuse and how to cite our work.
  - Documentation: (to be added)
  - Help: this page""")

st.header("MAVISp simple mode and MAVISp ensemble mode")

st.subheader("First look at the dataset")

st.write("""These two pages have very similar layout and content and they
are designed for users to download the content of the database and explore
its content.

Users are first shown a list of proteins to choose from, by clicking on the 
checkbox at the beginning of each row in the table. Each row also contains a 
link to the corresponding report, hosted on GitBook. In the case of ensemble 
mode, the table also includes a link to the Open Science Foundation repository 
where the structural ensemble is stored.

Once a protein of interest has been selected, users can choose betwen three
tabs to explore and visualize the data: XX, YY and ZZ.
""")

st.subheader("First exploration and download of the dataset")

st.write("""The Dataset tab contains, a table with the MAVISp dataset for
the selected protein. Users can select between two visualization formats
for the dataset: Compact and full. Both dataset contains one protein-level
mutation per row. The compact dataset is a more streamlined version of the
full dataset, containing only the most important outcomes of the MAVISp
framework (mostly, the effect classification for each module), while the
Full dataset contains all the availbile data in MAVISp. These datasets
only differ in the number of columns and contain the same mutations.

By using the button on each column header, which becomes visibile on when
hovering the mouse point over the header itself, users can further explore
and refine the dataset, by for instance adding or removing columns, hiding
or showing specific rows. Clicking on any table header allows to sort the
table by that column, in either ascending or descending order.""")

st.subheader("Plotting MAVISp classification for a custom set of mutations")

st.write("""This is performed by selecting the second tab, Classification. This
tab aims at creating a dot plot that shows the classification of each mutation
according to each of the available MAVISp module for this protein, as well as
for reference pathogenicity score. The user interface allows to customize whether
to show or hide the classification for the DeMaSk and REVEL predictors (if available) 
and configure the cut-off for the classification of the mutations done using such
predictors. It also allows to customize the size of the plot(s) and how many mutations
should be plotted per plot. If the total value of mutations plotted exceeds this
value, more than one plot will be generated.

The plots contain two color definitions: the top one refers to the dots in the table,
while the bottom one refers to the ClinVar classification of the mutations, which is
displayed as the color in the mutation label below the X axis. If the color is black,
no classification was available.

Finally the user is asked to select a number of mutations to plot in the dataset.
There is no hard limit to how many these should be, nonetheless, generating many
plot of the same time is considerably slower than just one at the time.""")

st.subheader("Plotting summary of MAVISp classification using a reduced set of classifications")

st.write("""The MAVISp website allows to download an even more succint representation
of the classification performed by MAVISp. This is carried out in the Damaging mutations tab.
This visualization only considers mutations that are at the same time i) classified as pathogenic for AlphaMissense,
ii) Predicted to cause loss of function or gain of function for either GEMME or DeMaSk, iii) are
found to be damaging in at least one MAVISp module. Consequences are grouped in the following
broad categories: stability, PTM, long range effects, local interactions. The user is asked
to select a number of mutations as done for the previous tab. It should be noticed that not
all the mutations available in the dataset will be available in this panel, as not all mutations
will satisfy the aforamentioned criteria.

The final representation is a lollipop plot, in which mutations are found on the X axis and
one or more vertical bars are present in corresponde of the mutation, displaying which of the
broad categories are affected by the mutation.""")

st.subheader("Visualizing a summary of MAVISp classification on the 3D structure")

st.write("""The MAVISp website also include an interactive 3D visualization of the
protein structure on which it is possible to visualize the same mutation sites
and mutations as per the previous section. It also allows to customize which
type of MAVISp classification should be considered. It supports two analysis modes:
one in which only the mutation sites that have at least a specified number of mutations
that are classified as damaging are displayed; and another in which the user
can freely select a set of mutation sites to be shown on the protein structure.""")
