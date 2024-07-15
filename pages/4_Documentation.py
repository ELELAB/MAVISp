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
import pandas as pd

st.set_page_config(layout="wide",
    page_title="Documentation",
    page_icon="📖")

st.header('Documentation')

This page contains an introduction to the MAVISp framework and a guide on the
available data and on its interpretation. If you're looking on a guide on how to
use the MAVISp website, please look into the Help section on the left.

st.subheader("Introduction to MAVISp")

A complete description of the MAVISp framework is available in the MAVISp manuscript:

ref

This documentation assumes the reader to be familiar with the main concepts exposed
in the paper, which contains in-depth information on how the framework is structured.

Briefly, MAVISp is a framework designed for the prediction on the effect of mutations
on proteins, using mostly structure-based methods. It has a modular structure, with each
module predicting a different possible consequence of a mutation. For instance, the `stability`
module predicts the effect of mutations on protein stability; the `long range` module predicts
long-range effects of mutations and so on. Modules use custom or already available
methodologies - see the acknowledgements section or the MAVISp paper to know more. It
should be noted that we also include as modules some that perform some preparatory steps
for the subsequent analyses, such as downloading protein structures. We will not refer to
them in this documentation as they are typically not user-facing.

We will now follow with a short description of each module, what data it delivers to the final
MAVISp dataset and how to interpret their results. Please notice that modules are described more
in-depth in the MAVISp paper.

st.subheader("Structure of a typical MAVISp dataset")
The dataset for a protein is organized as a table, where each row is a protein mutation, and each
column contains data that has been calculated for that specific mutation or for the wild-type
residue in that position, depending on the context. The dataset also contains specific classification
column that summarize the final classification outcome for a specific module.

Furthermore, MAVISp operates on either a single protein structure (simple mode) or an
ensemble of structures (ensemble mode). In ensemble mode, we occasionally include more than
one ensemble (for instance, from molecular dynamics simulations or NMR experiments) in the
dataset; for this reason, some of the columns will have a `[tag]` that defines which ensemble
the column refers to. For instance, the `Stability classification, (Rosetta, FoldX) [md]` column
refers to a MD ensemble.

Finally, it should be noted that not all data will be available for all the rows. This can depend
on the trimming we perform on the structural models we have available (e.g. because the stability
module would not be reliable if ran on disordered regions) or other factors.

    
st.subheader("Stability")

st.write("""The Stability module predicts changes of folding free energy upon mutation, i.e. how the stability
of the protein changes respect to the reference protein sequence (the wild-type). In this module,
we use both FoldX and Rosetta or RAsP to calculate changes of free energy of folding associated
with the mutation and build a consensus from their results as explained in the MAVISp paper.

A MAVISp dataset will typically have one or more columns named:

Stability (FoldX5, alphafold, kcal/mol)
Stability (Rosetta Cartddg2020, alphafold, kcal/mol)
Stability (RaSP, alphafold, kcal/mol)
Stability classification, alphafold, (Rosetta, FoldX)
Stability classification, alphafold, (RaSP, FoldX)

The row name includes the method with which the calculation has been performed, the
source of the structural model used for it and the unit the free energy changes are expressed in.

The `Stability` columns contain free energy changes values in kcal/mol associated to the
specific mutations, with positive values indicating that the mutation destabilizies the protein
structure and negative values indicating that the mutation improves the stability of the protein
structure. The final classification is derived from these values, as described in the MAVISp paper.
The classification column is always generated by considering FoldX and Rosetta, or FoldX and RaSP.
This means that the classification will only be available if the respective datasets are also available.

The possible classification values are:""")

data = [ ( 'Destabilizing', 'The mutation is destabilizing for the protein structure'),
         ( 'Stabilizing'  , 'The mutation is destabilizing for the protein structure'),
         ( 'Neutral'      , 'The mutation has no significant effect on stability'),
         ( 'Uncertain'    , 'The mutation has a border line effect on stability, or the two methods are not in agreement',
         ( 'N.A'          , 'No data available to perform the classification') ]
st.table(pd.DataFrame(data, columns=['Value', 'Meaning']))

st.subheader("Local interactions")

st.write("""The `Local interactions` module calculates the change of free energy of binding upon mutation
between our protein of interest and one target protein (heterodimer) or the same protein (homodimer),
starting from a protein complex structure. These are calculated either using FoldX and Rosetta;
a consensus approach is used to build up a final classification.

A MAVISp dataset may have one or more columns related to Local interactions, even though they are
currently available for a small fraction of the dataset. For example:

Local Int. (Binding with MAP1LC3B_AFmulti, heterodimer, FoldX5, kcal/mol)"
Local Int. (Binding with MAP1LC3B_AFmulti, heterodimer, Rosetta Talaris 2014, kcal/mol)"
Local Int. classification (MAP1LC3B_AFmulti),
Local Int. (Binding with OPTN_AFmulti, homodimer, FoldX5, kcal/mol)"
Local Int. (Binding with OPTN_AFmulti, homodimer, Rosetta Talaris 2014, kcal/mol)",
Local Int. classification (OPTN_AFmulti)

The column names will have information on the protein that forms a complex with our
protein of interest and on the origin of such complex structure, whether they are
in a hetero- or homodimer, the method used to perform the calculation and the unit in which
the changes in binding free energy are expressed. Positive value indicate a less strong
binding upon mutation, while negative values indicate a stronger binding.

Similarly as for the `Stability` module, MAVISp uses a consensus approach between FoldX
and Rosetta to calculate the final classification. Additionally, the module also consider
the relative solvent accessible surface area for the side chain of the wild-type residue (SAS)
in the considered structure of ensemble to classify cases for which information is not
available.

The possible classification values are:""")

data = [ ( 'Destabilizing', 'The mutation is destabilizes the binding between the two proteins'),
         ( 'Stabilizing'  , 'The mutation is stabilizes the binding between the two proteins'),
         ( 'Neutral'      , 'The mutation has no significant effect on the binding between the two proteins'),
         ( 'Uncertain'    , 'the two methods are not in agreement, or free energy values are not available and SAS >= 20%',
         ( 'N.A'          , 'free energy values are not available and SAS < 20%',') ]
st.table(pd.DataFrame(data, columns=['Value', 'Meaning']))

st.subheader("Local interactions with DNA")

st.write("""This module has similar purpose as `Local interactions`, just considering interactions
between our protein of interest and DNA. It uses FoldX only to calculate the change of binding free energy
upon mutation for our set of mutations of interest, so the classification is performed using FoldX energy
values only, together with SAS as described for `Local interactions`.

Columns generated by the `Local interactions with DNA` module look like:

Local Int. With DNA (modeller_2KK0_1KQQ, heterodimer, FoldX5, kcal/mol)
Local Int. classification With DNA

they follow a similar structure as those for `Local interactions`. The classification column che have the
following values:""")

data = [ ( 'Destabilizing', 'The mutation is destabilizes the binding between our protein and DNA'),
         ( 'Stabilizing'  , 'The mutation is stabilizes the binding between our protein and DNA'),
         ( 'Neutral'      , 'The mutation has no significant effect on the binding between our protein and DNA'),
         ( 'Uncertain'    , 'free energy values are not available and SAS >= 20%',
         ( 'N.A'          , 'free energy values are not available and SAS < 20%',') ]
st.table(pd.DataFrame(data, columns=['Value', 'Meaning']))
