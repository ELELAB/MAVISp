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

st.write("""This page contains an introduction to the MAVISp framework and a guide on the
available data and on its interpretation. If you're looking on a guide on how to
use the MAVISp website, please look into the Help section on the left.""")

st.subheader("Introduction to MAVISp")

st.write("""A complete description of the MAVISp framework is available in the MAVISp manuscript:

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
module would not be reliable if ran on disordered regions) or other factors.""")

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
         ( 'Uncertain'    , 'The mutation has a border line effect on stability, or the two methods are not in agreement'),
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
         ( 'Uncertain'    , 'the two methods are not in agreement, or free energy values are not available and SAS >= 20%'),
         ( 'N.A'          , 'free energy values are not available and SAS < 20%',) ]
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
         ( 'Uncertain'    , 'free energy values are not available and SAS >= 20%'),
         ( 'N.A'          , 'free energy values are not available and SAS < 20%',) ]
st.table(pd.DataFrame(data, columns=['Value', 'Meaning']))

st.subheader("Cancermuts table")

st.write("""MAVISp includes data from a table which is the output of Cancermuts.
Cancermuts is our Python package for gathering and annotating known cancer-associated
mutations from cancer genomics and mutations database, as well as to annotate the
mutations with a wealth of information. A pipeline running cancermuts is part of our
standard data collection strategy. From the Cancermuts output, we collect the following
information:""")

data = [ ( 'HGVSg'                           , 'Genomic mutation(s) cause of the protein mutation'),
         ( 'gnomAD genome allele frequency'  , 'gnomAD allele frequency on genomes dataset'),
         ( 'gnomAD exome allele frequency'   , 'gnomAD allelel frequency on exome dataset'),
         ( 'REVEL score'                     , 'REVEL score(s) associated with the genomic mutation(s)'),
         ( 'Mutation sources'                , 'which database or dataset the mutations were gathered from',) ]
st.table(pd.DataFrame(data, columns=['Column', 'Description']))

st.subheader("PTMs")

st.write("""The `PTMs` module tries to predict the effect of mutations on residues
that are affected by post-translational modifications. It currently supports 
phosphorylation only as well as phosphorylatable residues. The module predicts the
effect of the mutation on three different aspects: regulation, function and stability,
as detailed in the MAVISp paper. The final columns produced by this module are:""")

data = [ ( "PTMs", "whether this position was found to be phosphorylatable or not in the protein", "P for phosphorylatable residues, nothing otherwise"),
         ( "is site part of phospho-SLiM", "Whether this residue is included in a short linear motif that is known to be phosphorylatable", "True if it is, False if it's not"),
         ( "PTM residue SASA (%)" , "SAS for wild-type, unmodified residue in the selected structure or structural ensemble", "number (%)"),
         ( "Change in stability with PTM (FoldX5, kcal/mol)", "Change in folding free energy upon phosphorylation", "number (kcal/mol)"),
         ( "Change in binding with mutation (FoldX5, kcal/mol)", "Change in binding free energy upon mutation", "number (kcal/mol)"),
         ( "Change in binding with PTM (FoldX5, kcal/mol)", "Change in binding free energy upon phosphorylation", "number (kcal/mol)"),
         ( "PTM effect in regulation", "Final classification of the effect of PTM in terms of regulation", "see below"),
         ( "PTM effect in stability", "Final classification of the effect of PTM in terms of stability", "see below"),
         ( "PTM effect in function", "Final classification of the effect of PTM in terms of function", "see below") ]
st.table(pd.DataFrame(data, columns=['Column', 'Description', 'Possible values']))

st.write("""The PTM regulation classification predicts on whether the mutation
will have consequences on the functional regulation of the protein. The classification
follows the following flowchart:

The PTM stability classification predicts whether the presence of the mutation is
likely to have an effect on stability, by removing the possibility of a residue
to be phosphorylated. This is important because phosphorylation itself can have
an effect on stability. The classification follows the following flowchart:

The PTM function classification predicts whether the presence of the mutation is
likely to have consequences on function. In this context, we consider binding with
other protein as the function that we test. The classification follows the
following flowchart:""")

st.subheader("""Denovo Phospho""")

st.write("""This module predicts ...""")

st.subheader("""Long range""")

st.write("""This module predicts the capability of mutations to affect other residues
long range, through allosteric mechanisms. This module is based on the output of
the AlloSigma2 web server. The module processes the output of the web server to
understand whether the mutations are likely to have a long range effect.

Briefly, AlloSigma2 provides with an allosteric free energy value which is a measure
of how much a mutation at a certain position is likely to have a long-range effect on
any other residue of the protein (allosteric signaling map). AlloSigma2 considers two
possible classes of mutations: smaller to larger residues (UP mutations) or larger to
smaller residue (DOWN mutations). In MAVISp, we use residue a size cut-off to classify our
mutations in UP, DOWN or neither; for UP and DOWN mutations, we then identify if
other residues in the protein are affected by the mutation by using AlloSigma2.
Of these, we only consider residues that are likely to have a functional meaning,
by filtering them for their belonging to a pocket in the protein surface, identified
using fPocket. Alternatively, we also occasionally consider residues that are in
functional sites or active sites.

The dataset columns generated by this module are:""")

data = [ ( "AlloSigma2 mutation type", "Mutation type for Allosigma2", "DOWN, UP or empty"),
         ( "AlloSigma2 predicted consequence - pockets and interfaces", "Classification for long range effects", "see below") ]
st.table(pd.DataFrame(data, columns=['Column', 'Description', 'Possible values']))

st.write("""The classification column can have the following values:""")

data = [ ( 'destabilizing'  , 'The mutation has destabilizing long-range effects on the protein structure'),
         ( 'stabilizing'    , 'The mutation has stabilizing long-range effects on the protein structure'),
         ( 'mixed_effects'  , 'The mutation has both stabilizing an destabilizing long-range effects on the protein structure'),
         ( 'neutral'        , 'The mutation has neither stabilizing or destabilizing long-range effects on the protein structure'),
         ( 'uncertain'      , 'The mutation results in a too small change of side-chain volume to be considered either UP or DOWN'),
         ( 'N.A.'           , 'The mutation site could not be predicted (for instance, because it is located in an unstructured region)') ]
st.table(pd.DataFrame(data, columns=['Value', 'Meaning']))

st.subheader("")

damaging = if there is a path and allosigma2 has a prediction
uncertain = if allosigma predicted an effect but effect not validated by path analysis or mutation is not classified UP or DOWN
neutral = if the mutation is not found in the file BUT it was assigned UP or DOWN
NA = the mutation was not in allosigma_mut.Text


