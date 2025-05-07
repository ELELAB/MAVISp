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

import streamlit as st
import pandas as pd
from streamlit_utils import add_mavisp_logo, add_affiliation_logo
st.set_page_config(layout="wide",

    page_title="Documentation",
    page_icon="📖")

add_mavisp_logo("static/logo_small.png", image_width='50%')

add_affiliation_logo()

st.header('Documentation')

st.write("""This page contains an introduction to the MAVISp framework and a guide on the
available data and on its interpretation. If you're looking on a guide on how to
use the MAVISp website, please look into the Help section on the left.""")

st.subheader("Introduction to MAVISp")

st.write("""A complete description of the MAVISp framework is available in the
MAVISp manuscript (see Acknowledgements section).

This documentation assumes the reader to be familiar with the main concepts exposed
in the paper, which contains in-depth information on how the framework is structured.

Briefly, MAVISp is a framework designed for the prediction on the effect of mutations
on proteins, using mostly structure-based methods. It has a modular structure, with each
module predicting a different possible consequence of a mutation. For instance, the Stability
module predicts the effect of mutations on protein stability; the Long Range module predicts
long range effects of mutations and so on. Alternatively, modules can annotate features
of the structural model in use or of the protein unders study. Modules use custom or already available
methodologies - see the acknowledgements section or the MAVISp paper to know more.

It should be noted that MAVISp also includes modules that perform preparatory steps
to perform analyses for other modules, for instance identifying and obtaining
protein structures. We will not refer to them in this documentation as they are
typically not user-facing.""")

st.subheader("Structure of a typical MAVISp dataset")

st.write("""The dataset for a protein is organized as a table, where each row is a protein mutation, and each
column contains data that has been generated for that specific mutation or for the wild-type
residue in that position, depending on the context. The dataset also contains classification
columns that summarize the final classification outcome for a specific module.

Furthermore, MAVISp operates on either a single protein structure (simple mode) or an
ensemble of structures (ensemble mode). In ensemble mode, we occasionally include more than
one ensemble (for instance, from molecular dynamics simulations or NMR experiments) in the
dataset; for this reason, some of the columns will have a `[tag]` that defines which ensemble
the column refers to. For instance, the `Stability classification, (Rosetta, FoldX) [md]` column
refers to a MD ensemble. If this is the case, some columns might be repeated, once
per supported ensemble.

Finally, it should be noted that not all data will be available for all the rows. This can depend
on the trimming we perform on the structural models we have available (e.g. because the stability
module would not be reliable if ran on disordered regions) or other factors.

We will now follow with a short description of each module, what data it delivers
to the final MAVISp dataset and how to interpret their results.""")

st.subheader("SAS")

st.write("""The SAS (Solvent-Accessible Surface) module is designed to annotate
each residue (and therefore each mutation) with a solvent-accessible surface
area value calculated on the protein structure. This corresponds to the
relative side-chain solvent-accessible surface area in the wild-type protein
structure for simple mode, and the average over the ensemble of the same value
over the ensemble for ensemble mode. The SAS is calculated using the
NACCESS program and expressed in percentage.

Since the value is calculated by normalizing the absolute SAS by the SAS of the
residue in a Ala-X-Ala tripeptide, the SAS value might in some cases exceed 100%.
The SAS value is written in the \"Relative Side Chain Solvent Accessibility
in wild-type\" column.

It should be noted that the very same SAS values are used in other modules as
well - see below for more details""")

st.subheader("Stability")

st.write("""The Stability module predicts changes of folding free energy upon mutation, i.e. how the stability
of the protein changes respect to the reference protein sequence (the wild-type). In this module,
we use both FoldX and Rosetta or RAsP to calculate changes of free energy of folding associated
with the mutation and build a consensus from their results as explained in the MAVISp paper.

A MAVISp dataset will typically have one or more columns named:""")

data = [ ("Stability (FoldX5, alphafold, kcal/mol)",
          "change of folding free energy upon mutation calculated using FoldX",
          "value (kcal/mol)"),
         ("Stability (Rosetta Cartddg2020, alphafold, kcal/mol)",
         "change of folding free energy upon mutation calculated using Rosetta",
         "value (kcal/mol)"),
         ("Stability (RaSP, alphafold, kcal/mol)",
         "change of folding free energy upon mutation calculated using RaSP",
         "value (kcal/mol)"),
         ("Stability classification, alphafold, (Rosetta, FoldX)",
         "Consensus classification using Rosetta and FoldX data",
         "see below"),
         ("Stability classification, alphafold, (RaSP, FoldX)",
         "Consensus classification using Rosetta and FoldX data",
         "see below") ]
st.table(pd.DataFrame(data, columns=['Column', 'Description', 'Possible values']))

st.write("""The row name includes the method with which the calculation has been performed, the
source of the structural model used for it and the unit the free energy changes are expressed in.

The `Stability` columns contain free energy changes values in kcal/mol associated to the
specific mutations, with positive values indicating that the mutation destabilizies the protein
structure and negative values indicating that the mutation makes the protein more stable.
The final classification is derived from these values, as described in the MAVISp paper.
The classification column is always generated by considering FoldX and Rosetta, or FoldX and RaSP.
This means that the classification will only be available if the respective datasets are also available.

The possible classification values are:""")

data = [ ( 'Destabilizing', 'The mutation is destabilizing for the protein structure'),
         ( 'Stabilizing'  , 'The mutation is destabilizing for the protein structure'),
         ( 'Neutral'      , 'The mutation has no significant effect on stability'),
         ( 'Uncertain'    , 'The mutation has a borderline effect on stability, or the two methods are not in agreement'),
         ( 'N.A'          , 'No data available to perform the classification') ]
st.table(pd.DataFrame(data, columns=['Value', 'Meaning']))

st.subheader("Local interactions")

st.write("""This module calculates the change of free energy of binding upon mutation
between our protein of interest and one target protein (heterodimer) or the same protein (homodimer),
starting from a protein complex structure. These are calculated either using FoldX and Rosetta;
a consensus approach is used to build up a final classification. Currently, results
for the Local interactions module are only available for a small fraction of the
dataset.

The module generates multiple columns, and several of them might be present
depending on the system and on other factors. Typical columns look like:""")

data = [ ("Local Int. (Binding with MAP1LC3B_AFmulti, heterodimer, FoldX5, kcal/mol)",
          "change of binding free energy upon mutation calculated using FoldX",
          "value (kcal/mol)"),
("Local Int. (Binding with MAP1LC3B_AFmulti, heterodimer, Rosetta Talaris 2014, kcal/mol)",
          "change of binding free energy upon mutation calculated using FoldX",
          "value (kcal/mol)"),
("Local Int. classification (MAP1LC3B_AFmulti)",
          "change of binding free energy upon mutation calculated using FoldX",
          "see below"),
("Local Int. (Binding with OPTN_AFmulti, homodimer, FoldX5, kcal/mol)",
          "change of binding free energy upon mutation calculated using FoldX",
          "value (kcal/mol)"),
("Local Int. (Binding with OPTN_AFmulti, homodimer, Rosetta Talaris 2014, kcal/mol)",
          "change of binding free energy upon mutation calculated using FoldX",
          "value (kcal/mol)"),
("Local Int. classification (OPTN_AFmulti)",
          "change of binding free energy upon mutation calculated using FoldX",
          "see below")]
st.table(pd.DataFrame(data, columns=['Column', 'Description', 'Possible values']))

st.write("""The column names will have information on the protein that forms a complex with our
protein of interest and on the origin of such complex structure, whether they are
in a hetero- or homodimer, the method used to perform the calculation and the unit in which
the changes in binding free energy are expressed. Positive value indicate a less strong
binding upon mutation, while negative values indicate a stronger binding.

Similarly as for the `Stability` module, MAVISp uses a consensus approach between FoldX
and Rosetta to calculate the final classification. Additionally, the module also consider
the relative solvent accessible surface area for the side chain of the wild-type residue
(see SAS module) in the considered structure or ensemble to classify cases for
which free energy information is not available.

The possible classification values are:""")

data = [ ( 'Destabilizing', 'The mutation is destabilizes the binding between the two proteins'),
         ( 'Stabilizing'  , 'The mutation is stabilizes the binding between the two proteins'),
         ( 'Neutral'      , 'The mutation has no significant effect on the binding between the two proteins'),
         ( 'Uncertain'    , 'the two methods are not in agreement, or free energy values are not available and SAS >= 20%'),
         ( 'N.A'          , 'free energy values are not available and SAS < 20%',) ]
st.table(pd.DataFrame(data, columns=['Value', 'Meaning']))

st.write("""Notice that we typically calculate changes of free energy exclusively for
residues located at the binding interface, meaning that most residues will likely not feature an
associated free energy value. If such residues are solvent-exposed, they will still
be classified as Uncertain, even though they are not located on the
binding interface for the specific interactor of interest, according to the previously
outlined classification rules.""")

st.subheader("Local interactions with DNA")

st.write("""This module has similar purpose as `Local interactions`, just considering
interactions between our protein of interest and DNA. It uses FoldX only to calculate
the change of binding free energy upon mutation for our set of mutations of interest,
so the classification is performed using FoldX energy values only, together with
SAS as described for `Local interactions`.

Columns generated by the `Local interactions with DNA` module look like:""")

data = [ ("Local Int. With DNA (modeller_2KK0_1KQQ, heterodimer, FoldX5, kcal/mol)",
          "change of binding free energy upon mutation between protein and DNA",
          "value (kcal/mol)"),
         ("Local Int. classification With DNA",
          "classification of the mutation according to binding with DNA",
          "see below") ]
st.table(pd.DataFrame(data, columns=['Column', 'Description', 'Possible values']))

st.write("""they follow a similar structure as those for `Local interactions`.
The classification column che have the following values:""")

data = [ ( 'Destabilizing', 'The mutation destabilizes the binding between our protein and DNA'),
         ( 'Stabilizing'  , 'The mutation stabilizes the binding between our protein and DNA'),
         ( 'Neutral'      , 'The mutation has no significant effect on the binding between our protein and DNA'),
         ( 'Uncertain'    , 'free energy values are not available and SAS >= 20%'),
         ( 'N.A.'          , 'free energy values are not available and SAS < 20%',) ]
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

data = [ ( "PTMs",
           "whether this position was found to be phosphorylatable or not in the protein",
           "P for phosphorylatable residues, nothing otherwise"),
         ( "is site part of phospho-SLiM",
            "whether this residue is included in a short linear motif that is known to be phosphorylatable",
            "True if it is, False if it's not"),
         ( "PTM residue SASA (%)",
           "SAS for wild-type, unmodified residue in the selected structure or structural ensemble",
           "value (%)"),
         ( "Change in stability with PTM (FoldX5, kcal/mol)",
           "Change in folding free energy upon phosphorylation",
           "value (kcal/mol)"),
         ( "Change in binding with mutation (FoldX5, kcal/mol)",
           "Change in binding free energy upon mutation",
           "value (kcal/mol)"),
         ( "Change in binding with PTM (FoldX5, kcal/mol)",
           "Change in binding free energy upon phosphorylation",
           "value (kcal/mol)"),
         ( "PTM effect in regulation",
           "Final classification of the effect of PTM in terms of regulation",
           "see below"),
         ( "PTM effect in stability",
           "Final classification of the effect of PTM in terms of stability",
           "see below"),
         ( "PTM effect in function",
           "Final classification of the effect of PTM in terms of function",
           "see below") ]
st.table(pd.DataFrame(data, columns=['Column', 'Description', 'Possible values']))

st.write("""The PTM regulation classification predicts on whether the mutation
will have consequences on the functional regulation of the protein.

The PTM stability classification predicts whether the presence of the mutation is
likely to have an effect on stability, by removing the possibility of a residue
to be phosphorylated. This is important because phosphorylation itself can have
an effect on stability.

The PTM function classification predicts whether the presence of the mutation is
likely to have consequences on function. In this context, we consider binding with
other protein as the function that we test.

The three types of classification are performed by using a custom logic for each.
The classifications themselves and a flowchart on how they are derived is available
in the supplementary material of the MAVISp paper.""")

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

data = [ ( "AlloSigma2 mutation type",
           "Mutation type for Allosigma2",
           "DOWN, UP or empty"),
         ( "AlloSigma2 predicted consequence - pockets and interfaces",
           "Classification for long range effects",
           "see below") ]
st.table(pd.DataFrame(data, columns=['Column', 'Description', 'Possible values']))

st.write("""The classification column can have the following values:""")

data = [ ( 'destabilizing'  , 'The mutation has destabilizing long-range effects on the protein structure'),
         ( 'stabilizing'    , 'The mutation has stabilizing long-range effects on the protein structure'),
         ( 'mixed_effects'  , 'The mutation has both stabilizing an destabilizing long-range effects on the protein structure'),
         ( 'neutral'        , 'The mutation has neither stabilizing nor destabilizing long-range effects on the protein structure'),
         ( 'uncertain'      , 'The mutation results in a too small change of side-chain volume to be considered either UP or DOWN'),
         ( 'N.A.'           , 'The mutation site could not be predicted (for instance, because it is located in an unstructured region)') ]
st.table(pd.DataFrame(data, columns=['Value', 'Meaning']))

st.subheader("""Long range: AlloSigma2-PSN""")

st.write("""This module is available only for ensemble mode.

The module aims to validate the simple mode Long Range module predictions, performed on a single 
structure using AlloSigma2, by analysing an ensemble of structures using a Protein Structure Network
(PSN)-based method.

The module takes as input the predictions of the Long Range module, specifically the mutations 
predicted to affect pocket residues, that satisfy the thresholds described in the previous section. 

Briefly, these predictions are then validated by identifying potential 
communication paths between mutations and their respective response pocket sites using
[PyInteraph2](https://doi.org/10.1021/acs.jcim.3c00574).

PyInteraph2 derives an inter-residue network from the contacts between pairs of residues
in an ensemble of structures. Such network can then be analyzed to identify potential
pathways through which long-range happens.

As input for constructing the network, an ensemble structure is used (e.g., trajectory of an MD simulation). 
Subsequently, the path_analysis tool of PyInteraph2 is used to identify
shortest paths of communication between the mutation and pocket sites. 
We retain only paths that are 4 residues in length or more. 
If such a path is identified, it is considered a validation of the AlloSigma2 predicted allosteric 
effect of a mutation on the pocket residue(s).

The module aims to leverage these two methods to validate the predictions and apply a consensus 
approach to build the final classification of a mutation's long-range effects.

The dataset columns generated by this module are:""")

data = [ ( "AlloSigma2-PSN classification",
           "Classification of long range AlloSigma2-PSN effects",
           "see below") ]
st.table(pd.DataFrame(data, columns=['Column', 'Description', 'Possible values']))

st.write("""The classification column can have the following values:""")

data = [ ( 'damaging'       , 'The mutation effect on pocket site(s) was identified by both methods'),
         ( 'neutral'        , 'The mutation was not predicted to have any effect on pocket site(s)'),
         ( 'uncertain'      , 'The the two methods are not in agreement, or the mutation results in a too small change of side-chain volume to be considered either UP or DOWN in AlloSigma2') ]
st.table(pd.DataFrame(data, columns=['Value', 'Meaning']))

st.subheader("AlphaFold Metadata")

st.write("""This module is designed to add some information about the used
AlphaFold model (if any) to the databae. It generates two columns: one
containing the AlphaFold2 pLDDT score for each residue (\"AlphaFold2 model pLDDT
score\"), and another one including its secondary structure definition as
calculated by DSSP ("AlphaFold2 model secondary structure"), and therefore
expressed using the secondary structure dictionary used by the DSSP software.""")

st.subheader("ClinVar")

st.write("""This module adds to the table data regarding available ClinVar
classification of our mutations, when available. It should be noted that any
protein mutation can be associated with more than one ClinVar ID, because ClinVar
IDs typically refer to their genomic mutation, and multiple genomic mutations can
have the same protein consequence. The ClinVar module produces the following
columns:""")

data = [ ( "ClinVar Variation ID",
           "one or more ClinVar ID associated to the protein mutation",
           "comma-separated list of IDs"),
         ( "ClinVar Interpretation",
           "Clinvar interpretation associated with the variant IDs - one per ID",
           "comma-separated list of interpretations"),
         ( "ClinVar Review Status",
           "Review status associated with the variant",
           "comma-separated number, which is the number of stars associated with the variant") ]
st.table(pd.DataFrame(data, columns=['Column', 'Description', 'Possible values']))

st.subheader("DeMaSk, GEMME, EVE and AlphaMissense modules")

st.write("""Each of this modules is named after a different predictor of
pathogenicity for mutations, and it adds the results of the predictions of the
respective method to the MAVISp dataset. For some of them, we also report
their own or our own variant classification calculated from the available
scores. Each of these module will have its own columns:

  - DeMaSk module:""")

data = [ ( "DeMaSk delta fitness", "Delta fitness value from DeMaSk", "value"),
         ( "DeMaSk Shannon entropy", "Shannon entropy value for DeMaSk", "value"),
         ( "DeMaSk log2 variant frequency" , "log2 of DeMaSk variant frequency", "value"),
         ( "DeMaSk predicted consequence", "MAVISp-predicted consequence for DeMaSk",
            "gain_of_function, when Delta fitness > 0\n\
            loss_of_function, when Delta fitness < 0\n\
            neutral, when Delta fitness = 0" ) ]
st.table(pd.DataFrame(data, columns=['Column', 'Description', 'Possible values']))

st.write("""  - GEMME module:""")

data = [ ( "GEMME Score", "score from GEMME", "value"),
         ( "GEMME Score (rank-normalized)", "Rank-normalized gemme score", "value") ]
st.table(pd.DataFrame(data, columns=['Column', 'Description', 'Possible values']))

st.write("""  - EVE module:""")

data = [ ( "EVE score", "score from EVE", "value"),
         ( "EVE classification (25% Uncertain)", "Classification performed by EVE at 25% uncertainty", "Benign, Uncertain or Pathogenic") ]
st.table(pd.DataFrame(data, columns=['Column', 'Description', 'Possible values']))

st.write("""  - AlphaMissense module:""")

data = [ ( "'AlphaMissense pathogenicity score", "pathogenicity sore from AlphaMissense", "value"),
         ( "AlphaMissense classification", "Classification of mutation by AlphaMissense", "benign, pathogenic or ambiguous") ]
st.table(pd.DataFrame(data, columns=['Column', 'Description', 'Possible values']))

st.subheader('EFoldMine module')

st.write("""This module uses results from the EFoldMine software to annotate early
folding regions in the protein. EFoldMine expresses the likelihood of any residue
to be part of an early folding region as a score. MAVISp identifies continuous
stretches of residues that are above a pre-determined threshold of a minimum
length of 3 residues to assign whether each residue is part of a early folding
region. This module produces the following columns:""")

data = [ ( "'EFoldMine score", "Score from EFoldMine", "value"),
         ( "EFoldMine - part of early folding region",
           "Whether the residue was predicted to be part of an early folding region",
           "`True` if it is, `False` otherwise") ]
st.table(pd.DataFrame(data, columns=['Column', 'Description', 'Possible values']))

st.subheader("De novo Phosphosites module")

st.write("""this module was designed to identify whether each mutation has the
potential of adding a new phosphorylation site to the protein. It uses the Netphos
method to predict phosphorylation sites in both the wild-type and each mutant, and
compares the results to identify novel phosphorylation sites. Information about
solvent accessibility is also included in the prediction. It generates the
following columns:""")

data = [ ( "Phosphorylation - gain of function",
           "New phosphorylation sites that are predicted to be available upon mutation",
           "see below"),
         ( "Phosphorylation - loss of function",
           "Phosphorylation sites that are predicted for the wild-type but not predicted to exist upon mutation",
           "see below"),
         ("Mutation predicted to add new phosphorylation site",
          "Whether a mutation is predicted to add a new phosphorylation site",
          "`True` if it is, `False` if it isn't") ]
st.table(pd.DataFrame(data, columns=['Column', 'Description', 'Possible values']))

st.write("""The gain and loss of function columns express whether the mutation has
caused a new phosphosite to be predicted, or if the mutation has caused
a phosphosite to not be predicted anymore respect to the wild-type, respectively.
If no value is present in such column, than the mutation had no effect in this
regard; if themutation had an effect, the column will contain the phosphorylation
site(s) together with the kinase that would phosphorylate the site.""")

st.subheader("Functional Sites module")

st.write("""This module annotates whether mutations have a known effect on functional
aspects of the protein, starting from manualy annotated data. The module produces
the following columns:""")

data = [ ( "Functional sites (cofactor)",
           "Consequence that the mutation can have on the binding of a cofactor",
           "`neutral` or `damaging`"),
         ( "Functional sites (active site)",
           "Consequence that the mutation can have on the active site of the protein",
           "`neutral` or `damaging`") ]
st.table(pd.DataFrame(data, columns=['Column', 'Description', 'Possible values']))

st.subheader("Functional Dynamics module")

st.write("""this module is only available for ensemble mode. It doesn't correspond
to any single specific analysis or method - it is used to annotate known effects
of any mutation to protein dynamics, when this is supposed to affect function.
The columns produced by this module are simply called "Functional dynamics (evidence)",
where "evidence" is the type of evidence that column expresses. This column contains
whatever classification makes sense for the type of evidence. One or more columns
with a similar format can be specified for the same protein.""")

st.subheader('Pfam module')

st.write("""This module uses results retrieved from Pfam, reporting the identified
         domains for the given protein and the respective residue ranges for each domain.
         This module produces the following columns:""")

data = [ ( "Pfam domain classification",
           "Pfam domains associated with the given residue in the protein",
           "single string reporting the Pfam domain description and accession code") ]
st.table(pd.DataFrame(data, columns=['Column', 'Description', 'Possible values']))

st.subheader('TED module')

st.write("""This module uses results retrieved from The Encyclopedia of Domains (TED), reporting the identified
         domain annotation for the given protein and the respective residue ranges for each domain.
         This module produces the following columns:""")

data = [ ( "TED-CATH domain classification",
           "TED domains and their associated CATH-assigned labels for the given residue in the protein",
           "single string reporting the CATH-label code") ]
st.table(pd.DataFrame(data, columns=['Column', 'Description', 'Possible values']))