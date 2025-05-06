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
    page_title="Acknowledgements and data usage",
    page_icon="🙏")

add_mavisp_logo("static/logo_small.png", image_width='50%')

add_affiliation_logo()

st.title('Acknowledgements and data usage')

st.subheader("Citing us")

st.write("""MAVISp is developed by the Cancer Structural Biology lab at the
Danish Cancer Institute, Copenhagen, Denmark, and the Cancer Systems Biology lab
at the Danish Technical University, Lyngby, Denmark. The project is led by
Elena Papaleo, head of both labs.""")

st.write("""If you use MAVISp in your research, please cite
  [our preprint](https://www.biorxiv.org/content/10.1101/2022.10.22.513328v4):""")

st.code("""MAVISp: A Modular Structure-Based Framework for Genomic Variant Interpretation
Matteo Arnaudi, Ludovica Beltrame, Kristine Degn, Mattia Utichi, Simone Scrima,
Pablo Sanchez Izquierdo, Karolina Krzesinska, Francesca Maselli, Terezia Dorcakova,
Jordan Safer, Alberte Heering Estad, Katrine Meldgard, Philipp Becker, Julie Bruun Brockhoff,
Amalie Drud Nielsen,  View ORCID ProfileValentina Sora, Alberto Pettenella, Jeremy Vinhas,
Peter Wad Sackett, Claudia Cava, Anna Rohlin, Mef Nilbert, Sumaiya Iqbal, Matteo Lambrughi,
Matteo Tiberti, Elena Papaleo. bioRxiv https://doi.org/10.1101/2022.10.22.513328""", language=None)

st.subheader("Data and software availability")

st.markdown('''We thank all the people that have contributed to either the MAVISp codebase or
as MAVISp data curators over time. Their contributions can be found,
respectively, on our [GitHub repository](https://www.github.com/ELELAB/MAVISp)
and in the dataset tables in the simple or ensemble mode sections.

We also would like to acknowledge our data sources and software dependencies,
without which MAVISp would not have been possible. The datasets we have used are
available under a variety of licensing terms, which are listed below for reuse.

Unless stated otherwise in the following text, MAVISp data is released under
the [Creative Commons Attribution 4.0 International (CC BY 4.0) license](https://creativecommons.org/licenses/by/4.0/).

Bulk download of MAVISp data tables is available at [our OSF repository](https://osf.io/ufpzm/)

The MAVISp software and source code for this website is available at our
[MAVISp GitHub repository](https://github.com/ELELAB/MAVISp)
and released under open source license.''')

st.subheader('Mutations')

st.markdown('''MAVISp uses the following datasets for the mutation data:

  - [**COSMIC**](https://cancer.sanger.ac.uk/cosmic): The Catalogue Of Somatic
    Mutations In Cancer. We use version 96 in our current releases. The data in
    MAVISp is released according to the COSMIC Non-Commercial Terms and
    Conditions, and released in MAVISp for non-commercial use only.
  - [**cBioPortal**](https://www.cbioportal.org/) for cancer genomics. Data from
    cBioPortal present in MAVISp is released under the [Open Data Commons Open
    Database License (ODbL) v1.0](https://opendatacommons.org/licenses/odbl/1-0/)
  - [**ClinVar**](https://www.ncbi.nlm.nih.gov/clinvar/), the NIH database of
    clinically relevant genetic variants. It is released under public domain.''')

st.image("static/NCBI_powered.png")

st.subheader('Mutations metadata')

st.markdown('''MAVISp uses the following datasets for the mutation metadata:

   - [**REVEL**](https://sites.google.com/site/revelgenomics/) scores were
   downloaded by [myvariants.info](https://myvariant.info/), whose license
   terms are available on their website. REVEL scores are available as
   integrated by dbSNP which is released under public domain.
   - [**gnomAD**](https://gnomad.broadinstitute.org/), database of human genetic
    variation. It is released under the public domain.''')

st.subheader('Protein metadata and predictions')

st.markdown('''MAVISp uses the following datasets or software for the prediction
or annotation of protein-level features, or structure:

   - [**UniProt**](https://www.uniprot.org). It is released under the  Creative
  Commons Attribution 4.0 International (CC BY 4.0) License
  - the [**Protein Data Bank**](https://www.rcsb.org). It is released under
  public domain.
  - the [**AlphaFold Protein Structure Database**](https://alphafold.ebi.ac.uk). It
  is released under the Creative Commons Attribution 4.0 International (CC BY 4.
  0 DEED) License.
  - the [**AlphaFill**](https://alphafill.eu) database of protein ligands for
  AlphaFold structures.
  - [**MobiDB**](https://mobidb.org/) for the prediction of protein disorder.
  It is released under the Creative Commons Attribution 4.0 International (CC
  BY 4.0 DEED) License.
  - [**ELM**](http://elm.eu.org), the Eukaryotic Linear Motif resource for functional sites in proteins.
  - [**PhosphoSitePlus**](https://www.phosphosite.org), a database of
  post-translational modifications. Data from PhosphoSitePlus in MAVISp, (i.e.
  the PTMs column in our database files) is released according to the
  [PhosphoSitePlus' terms and conditions](https://www.phosphosite.org/staticDownloads),
  and it is not available for commercial use.''')

st.subheader('Prediction of the effect of mutations')

st.markdown('''MAVISp uses the following software for the prediction of the effect of mutations:

  - [**AlloSigma**](http://allosigma.bii.a-star.edu.sg/home/), a webserver for the prediction of allosteric effects.
  - [**FoldX**](http://foldxsuite.crg.eu/), a software suite for the prediction of changes of protein stability upon mutation. FoldX was licensed to us according to the Academic license terms.
  - [**Rosetta**](https://www.rosettacommons.org/), a software suite for the prediction of protein structure and stability. We have used the Rosetta suite licensed to us according to the Academic license terms, meaning that data produced with Rosetta is released in MAVISp for non-profit research
    purposes only, in accordance to such license terms.
  - [**RaSP**](https://github.com/KULL-Centre/_2022_ML-ddG-Blaabjerg/), a machine learning method for the prediction of changes of protein stability upon mutation. We have used RaSP by means of [our pipeline](https://github.com/ELELAB/RaSP_workflow), which is derived from and similar to the original implementation
  - [**Demask**](https://demask.princeton.edu), released under GPL v.3 license
  - [**EVE**](https://evemodel.org), released under the MIT license
  - [**GEMME**](http://www.lcqb.upmc.fr/GEMME/download.html)
  - [**AlphaMissense**](https://github.com/google-deepmind/alphamissense). In
    particular, we rely on the [publicly available predictions v.2](https://zenodo.org/records/8360242),
    released under the [CC BY-NC-SA 4.0 License](https://creativecommons.org/licenses/by-nc-sa/4.0/).
    This means that AlphaMissense predictions available in MAVISp are strictly for non-commercial use only.''')

st.subheader('Protein-protein interaction databases')

st.markdown('''MAVISp uses the following protein-protein interaction databases:

  - [**mentha**](https://mentha.uniroma2.it), a database of protein-protein
  interactions.
  - [**HPC-Atlas**](http://www.yulpan.top/HPC-Atlas/), a database of protein
  complexes.''')
