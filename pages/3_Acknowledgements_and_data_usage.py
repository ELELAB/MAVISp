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

st.header('Acknowledgements and data usage')

st.write('''MAVISp is developed by the Cancer Structural Biology lab at the
Danish Cancer Institute, Copenhagen, Denmark, and the Cancer Systems Biology lab
at the Danish Technical University, Lyngby, Denmark. The project is led by
Elena Papaleo, head of both labs.

We thank all the people that have contributed to either the MAVISp codebase or
as MAVISp data curators over time. Their contributions can be found,
respectively, on our [GitHub repository](https://www.github.com/ELELAB/MAVISp)
and in the dataset tables in the simple or ensemble mode sections.

We also would like to acknowledge our data sources and software dependencies,
without which MAVISp would not be possible. Each dataset we have used has
different licensing terms, which are listed below for reuse. Unless stated
otherwise, MAVISp data is released under the [Creative Commons Attribution 4.0
International (CC BY 4.0) license](https://creativecommons.org/licenses/by/4.0/)''')

st.subheader('Mutations')

st.markdown('''MAVISp uses the following datasets for the mutation data:

  - [**COSMIC**](https://cancer.sanger.ac.uk/cosmic): The Catalogue Of Somatic Mutations In Cancer. We use version 96 in our current releases. The data in MAVISp is released according to the COSMIC Non-Commercial Terms and Conditions, and available for non-commercial use.
  - [**cBioPortal**](https://www.cbioportal.org/) for cancer genomics. Data from
  cBioPortal present in MAVISp is released under the [Open Data Commons Open Database License (ODbL) v1.0](https://opendatacommons.org/licenses/odbl/1-0/)
  - [**ClinVar**](https://www.ncbi.nlm.nih.gov/clinvar/)

''')