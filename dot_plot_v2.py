#!/usr/bin/env python

# Copyright (C) 2023 Ludovica Beltrame <beltrameludo@gmail.com>,
# Simone Scrima <simonescrima@gmail.com>, Karolina Krzesińska <kzokr@dtu.dk>,
# Matteo Tiberti
# Danish Cancer Society & Technical University of Denmark,
# 2026 Eszter Toldi, Technical University of Denmark

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
import math
import argparse
import argcomplete
import logging as log
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
import re


def adjust_figsize_w(x_mut, original_w, original_m):
    ''' Adjust width of figuresize if the last plot has less mutations
    on the x axis.

    The function takes as input the number of remaining mutations and
    the original width and number of mutations per axis, and returns
    the new width adapted to remaining mutations.
    It also ensures that the figure width does not exceed certain
    or become too small for proper display of legends.

    Parameters
    ----------
    x_mut: int
        number of mutations of the last plot
    original_w: float
        original figuresize width
    original_m: int
        original number of mutations per x-axis

    Returns
    ----------
    new_w: float
        adapted width

    '''
    # Define the minimum and maximum figure widths
    min_width = 10.5
    max_width = 20

    # Define scale factor
    scale_factor = (original_w - min_width) / np.sqrt(original_m)

    # Get the new width
    new_w = min_width + scale_factor * np.sqrt(x_mut)

    # Ensure the new width stays within the defined bounds
    new_w = max(min_width, min(new_w, max_width))

    return new_w

def color_markers(data):
    ''' Assign a color to the dots.

    The function takes as input the array and returns a list of colors
    to be plotted on the final matrix.

    Parameters
    ----------
    data: np.array
        array containing the effect information from the dataframe

    Returns
    ----------
    colors: list
        list of colors associated to the effectt

    '''

    # Dictionary of the colors associated to the effect
    # Damaging : 1 : "black"
    # Neutral : 2 : "beige"
    # Not available : 3 : "white"
    # Uncertain : 4 : "lightgray"
    # Loss of function : 'chinarose'
    # Gain of function : 'payne's gray'
    # Early folding region True : 'penn blue'
    mapping_colors = {1 : 'black',
                      5 : 'black',
                      6 : 'black',
                      2 : '#d5bdaf',
                      3 : 'white',
                      4 : 'lightgray',
                      7 : '#8E5466',
                      8 : '#636D87',
                      9 : '#9BB295'}

    # Make np.zeros array with the shape of the plotted grid
    zeros_array = np.zeros([data.shape[0], data.shape[1]], dtype = int)

    # Make color array mapping the effect to the colors in the
    # mapping_colors dictionary
    c_array = np.vectorize(mapping_colors.get)(data)

    # Combine the zeros_array and the c_array in a 3D array
    colors_array = np.stack([zeros_array,c_array])

    # Define the color list to be used for plotting
    colors = []

    # Iterate over the vertical dimension of the array
    for y in range(data.shape[1]):

        # Iterate over the horizontal dimension of the array
        for x in range(data.shape[0]):

            # Append the color corresponding to the position
            # in the array
            colors.append(colors_array[1,x,y])

    return colors

def color_yticklabels(labels):
    ''' Assign a color to the ylabel according to the type of effect.

    The function takes as input a list of y lables and returns a list
    of colors associated to the label depending on the kind of effect.

    Parameters
    ----------
    labels: list
        list of ylabels

    Returns
    ----------
    tick_colors: list
        list of colors associated to the ylabels

    '''
    # Make empty out list
    tick_colors = []

    # Define the prefix to classify the labels as having
    # an effect on stability, function or other
    stability = ['Stability classification',
                 'Phosphorylation effect in stability',
                 "Classification of change in folding free energy with phosphorylation",
                 'EFoldMine - part of early folding region',
                 'Loss of disulfide bridge',
                 'Predicted de-novo disulfide bridge']


    function = ['Local Int.',
                "Classification of change in binding free energy with phosphorylation",
                'Functional sites (cofactor)',
                'Functional sites (active site)',
                'Phosphorylation',
                'AlloSigMA']
    # Iterate over the labels and append to the color list
    # the color based on the effect
    for l in labels:

        # Stability
        if l.startswith(tuple(stability)):
            tick_colors.append('#6a4c93')

        # Function
        elif l.startswith(tuple(function)):
            tick_colors.append('#1982c4')

        # Other scores
        else:
            tick_colors.append('#AACC72')
            # og color 8ac926

    return tick_colors


def color_xticklabels(mutations, clinvar_col):
    ''' Assign a color to the xtick label according to ClinVar Interpretation.

    Parameters
    ----------
    mutations: list
        List of xlabels i.e. mutation indices
    clinvar_col: pandas.Series
        ClinVar Interpretation values associated with the mutations

    Returns
    ----------
    tick_colors: list
        List of colors associated with the xlabels
    '''
    # Define ClinVar color map
    clinvar_colors = {
        'conflicting interpretation': '#EBD479',
        'pathogenic': '#A71B32',
        'likely pathogenic': '#A71B32',
        'benign': '#3A8752',
        'likely benign': '#3A8752',
        'uncertain': '#878E99',
        "oncogenic": "#A71B32",
        "likely oncogenic": "#A71B32",
        "strong": "#A71B32",
        "potential": "#A71B32",
        "unknown": "#878E99",
        "benign/likely benign": "#3A8752"}

    # Create a list to store the colors
    tick_colors = []

    # Iterate through mutations and assign colors
    for mutation in mutations:
        clinvar_value = clinvar_col.get(mutation, None)
        if pd.isna(clinvar_value):
            tick_colors.append('black')
            continue

        # Find the color based on the presence of patterns
        clinvar_value_lower = clinvar_value.strip().lower()
        color_assigned = clinvar_colors.get(clinvar_value_lower, 'black')

        tick_colors.append(color_assigned)

    return tick_colors

def convert_to_float(value):
    ''' Read REVEL score and convert it to float.

    The function takes as input a a REVEL score as a string
    and converts it to float for subsequent analyses.

    Parameters
    ----------
    value: str
        REVEL score as it is provided in the MAVISp table
    Returns
    ----------
    value: float
        converted REVEL score. If mutliple REVEL scores are found,
        their average is returned

    '''
    if isinstance(value, str):

        # If multiple REVEL scores
        if ',' in value:

            # Split the string and convert each value to float
            values = [float(x.strip()) for x in value.split(',')]

            # Calculate the average of the float values
            return np.mean(values)

        # If single REVEL score
        else:
            return float(value)
    else:
        return value

def get_clinvar_columns(df, clinvar_class_type: str):
    split = { "aggregated": ("ClinVar Interpretation", "ClinVar Review Status"),
              "germline": ("ClinVar Germline Interpretation", "ClinVar Germline Review Status"),
              "oncogenicity": ("ClinVar Oncogenicity Interpretation", "ClinVar Oncogenicity Review Status"),
              "clinical_impact": ("ClinVar Clinical Impact Interpretation", "ClinVar Clinical Impact Review Status")}
    pair = split[clinvar_class_type]
    if all(col in df.columns for col in pair):
        return pair
    if clinvar_class_type == "aggregated":
        raise ValueError( "ClinVar mode 'aggregated' requested (default), but aggregated columns are missing. "
                          "This looks like a 'new style' entry csv. "
                          "Use -cct germline|oncogenicity|clinical_impact.")
    raise ValueError(f"ClinVar mode '{clinvar_class_type}' requested, but the corresponding columns are missing.")


def process_input(full_df, r_cutoff, p_cutoff, d_cutoff, g_cutoff, residues, mutations,
    clinvar_dict, plot_Revel, plot_popEVE, plot_Demask, plot_Source, plot_Clinvar, color_Clinvar, clinvar_cols):
    ''' Read MAVISp aggregated table.

    The function takes as input a MAVISp csv file and returns
    a dataframe formatted in a plot-compatible way.

    Parameters
    ----------
    full_df: dataframe
	   input MAVISp dataframe
    r_cutoff: float
        REVEL score cutoff
    d_cutoff: float
        DeMaSk cutoff
    g_cutoff: float
        GEMME cutoff
    residues: list or None
        list of residues to be considered
    mutations: list or None
        list of mutations to be considered
    clinvar_dict: pandas.DataFrame
        dataframe containing the ClinVar dictionary
    plot_Revel: bool
        whether to keep Revel results
    plot_Demask: bool
        whether to keep DeMaSk results
    plot_Source: list
        list of sources to be included (if any)
    plot_Clinvar: list
        list of ClinVar
    color_Clinvar: bool
        whether to color output for ClinVar

    Returns
    ----------
    plot_df: dataframe
        copy of the processed dataframe trimmed according to plotting flags.
    df: dataframe
        processed dataframe retaining all classification columns for downstream filters.
    full_df: dataframe
              mavisp csv without altered columns

    '''

    # Assert the required column is present


    if (plot_Clinvar or color_Clinvar) and clinvar_cols[0] not in full_df.columns:
        log.error(f"{clinvar_cols[0]} column is missing from the input data.")
        raise TypeError(f"{clinvar_cols[0]} column is missing from the input data.")

    f = lambda x: '(ThermoMPNN)' in x or \
                    '(Foldetta from FoldX and Rosetta)' in x or \
                    '(Foldetta from FoldX and RaSP)' in x or \
                    '(Rosetta, FoldX)' in x or \
                    '(RaSP, FoldX)' in x or \
                    "Classification of change in folding free energy with phosphorylation" in x or \
                    'Local Int. classification' in x or \
                    'Local Int. With DNA classification' in x or \
                    "Classification of change in binding free energy with phosphorylation" in x or \
                    'Functional sites (cofactor)' in x or \
                    'Functional sites (active site)' in x or \
                    'AlloSigMA 2 predicted consequence - active sites' in x or \
                    'AlloSigMA 2 predicted consequence - cofactor sites' in x or \
                    'AlloSigMA 2 predicted consequence - pockets and interfaces' in x or \
                    ('AlloSigMA2-PSN classification' in x and not 'AlloSigMA 2 mutation type' in x) or\
                    'Phosphorylation effect in ' in x or 'REVEL score' in x or \
                    'popEVE score' in x or \
                    'EVE classification (25% Uncertain)' in x or \
                    'DeMaSk delta fitness' in x or \
                    'DeMaSk predicted consequence' in x or \
                    ('GEMME Score' in x and not 'GEMME Score (rank-normalized)' in x) or \
                    'AlphaMissense classification' in x or \
                    ('Mutation' in x and not 'Mutation sources' in x and not 'Mutation predicted to' in x) or \
                    'EFoldMine - part of early folding region' in x or \
                    'Experimental data classification' in x or \
                    'Loss of disulfide bridge' in x or \
                    'Predicted de-novo disulfide bridge' in x

    df = full_df.copy()
    filter_columns = [col for col in full_df.columns if f(col)]
    df = df[filter_columns]

    # if both these columns are not present or absent, exit immediately
    if ('Loss of disulfide bridge' in filter_columns) ^ ('Predicted de-novo disulfide bridge' in filter_columns):
        raise TypeError ('Columns produced by the disulfide bridges module should either be both '
                         'available or both unavailable')

    # If pltRevel is True then convert REVEL score column to float
    # In case multiple REVEL scores are present it returns the average
    # of the values

    for d in [df, full_df]:
        d['REVEL score'] = d['REVEL score'].apply(convert_to_float)

        # Convert popEVE score column to numeric values.
        # Invalid/non-numeric values are converted to NaN.
        if 'popEVE score' in d.columns:
            d['popEVE score'] = pd.to_numeric(d['popEVE score'], errors='coerce')


    # Add REVEL score interpretation column
    df['REVEL'] = np.where(df['REVEL score'].isna(), None,
                                np.where(df['REVEL score'] >= r_cutoff,
                                    'Damaging', 'Neutral'))

    try:
        # Add popEVE score interpretation column
        df['popEVE'] = np.where(df['popEVE score'].isna(), None,
                                np.where(df['popEVE score'] < p_cutoff,
                                         'Damaging', 'Neutral'))
    except:
        log.warning(f'- no popEVE found in MAVISp csv.')

    try:
        # Convert GEMME score into absolute value
        df['GEMME predicted consequence'] = np.where(
        df['GEMME Score'].isna(), None,
        np.where(df['GEMME Score'] <= g_cutoff, 'Damaging', 'Neutral')
    )
    except:
        log.warning(f'- no GEMME found in MAVISp csv.')

    # Convert Demask delta fitness into absolute value
    df['DeMaSk predicted consequence'] = np.where(df['DeMaSk delta fitness'].isna(), None,
                np.where(df['DeMaSk delta fitness'] >= d_cutoff,'GoF',
                np.where(df['DeMaSk delta fitness'] <= -d_cutoff,'LoF',
                'Neutral')))

    # Drop score columns
    score_cols_to_drop = ['REVEL score',
                        'DeMaSk delta fitness']

    if 'GEMME Score' in df.columns:
        score_cols_to_drop.append('GEMME Score')

    if 'popEVE score' in df.columns:
        score_cols_to_drop.append('popEVE score')

    df.drop(columns=score_cols_to_drop, inplace=True)


    # Sort columns based on broad effect categories
    functional_cols = [col for col in df.columns if ('functional' in col.lower()
                                                   and 'experimental data classification' not in col.lower())
                       or "Classification of change in binding free energy with phosphorylation" in col]
    disulfide_cols = [col for col in df.columns if 'disulfide bridge' in col.lower()]
    stability_cols = [col for col in df.columns if ('stability' in col.lower()
                                                  and 'experimental data classification' not in col.lower())
                      or "Classification of change in folding free energy with phosphorylation" in col]
    stability_cols += disulfide_cols
    efold_col = [col for col in df.columns if 'efoldmine' in col.lower()]
    other_cols = [col for col in df.columns if 'functional' not in col.lower()
                  and 'stability' not in col.lower()
                  and "Classification of change in folding free energy with phosphorylation" not in col
                  and "Classification of change in binding free energy with phosphorylation" not in col
                  and 'disulfide bridge' not in col.lower()
                  and 'experimental data classification' not in col.lower()
                  and 'efoldmine' not in col.lower()]
    experimental_cols = []
    for col in df.columns:
        if 'Experimental data classification' in col and col not in experimental_cols:
            experimental_cols.append(col)

    # Get ensemble labels only from stability columns
    ensembles = []
    for col in stability_cols:
        m = re.search(r"\[(.*?)\]", col)
        ensemble = m.group(1).lower() if m else "simple"
        if ensemble not in ensembles:
            ensembles.append(ensemble)

    # Reorder stability columns by ensemble label, Foldetta last
    ordered_stability_cols = []

    for ensemble in ensembles:
        if ensemble == "simple":
            ensemble_cols = [
                col for col in stability_cols
                if "[" not in col
            ]
        else:
            ensemble_cols = [
                col for col in stability_cols
                if f"[{ensemble}]" in col.lower()
            ]

        foldetta_cols = [
            col for col in ensemble_cols
            if 'foldetta' in col.lower()
        ]
        other_stability_cols = [
            col for col in ensemble_cols
            if 'foldetta' not in col.lower()
        ]

        ordered_stability_cols += foldetta_cols + other_stability_cols

    stability_cols = ordered_stability_cols

    # Combine lists in order
    df = df[stability_cols + efold_col + functional_cols + other_cols + experimental_cols]

    # Reorder rows so VEP classifiers appear at the top of the plot (Alpha, EVE, GEMME, REVEL),
    # followed by DeMaSk and Experimental rows, while keeping the remaining blocks unchanged.
    preferred_vep_order = [
        'AlphaMissense classification',
        'EVE classification (25% Uncertain)',
        'GEMME predicted consequence',
        'REVEL',
        'popEVE'
    ]
    demask_pred_col = 'DeMaSk predicted consequence'
    experimental_present = [col for col in experimental_cols if col in df.columns]

    remaining_cols = [
        col for col in df.columns
        if col not in preferred_vep_order
        and col != demask_pred_col
        and col not in experimental_present
    ]

    ordered_cols = remaining_cols + experimental_present
    if demask_pred_col in df.columns:
        ordered_cols.append(demask_pred_col)

    for col in reversed(preferred_vep_order):
        if col in df.columns:
            ordered_cols.append(col)

    df = df[ordered_cols]

    # Define a dictionary of effect: code_number
    # Damaging : 1, 5, 6
    # Neutral : 2
    # Not available : 3
    # Uncertain : 4
    # Added support for EVE
    # -> pathogenic = 1
    # -> benign = 2
    # -> uncertain = 4
    # -> loss of function = 7
    # -> gain of function = 8
    # -> early folding region = 9
    effect_code = {'(?i)destabilizing': 1,
                   '(?i)pathogenic': 1,
                   '(?i)damaging': 1,
                   'mixed_effects': 5,
                   '(?i)stabilizing': 6,
                   '(?i)neutral': 2,
                   '(?i)benign': 2,
                   '(?i)uncertain': 4,
                   'ambiguous': 4,
                   'loss_of_function' : 7,
                   '(?i)loss[ _]?of[ _]?function' : 7,
                   '(?i)LoF' : 7,
                   'gain_of_function' : 8,
                   '(?i)gain[ _]?of[ _]?function' : 8,
                   '(?i)GoF' : 8,
                    True : 9,
                    False : 2}

    # First fill None values with code
    df.fillna(value=3, inplace=True)

    # Replace the string nomenclature with the
    # effect code
    df.replace(effect_code,
               regex = True,
               inplace = True)

    # Ensure 'DeMaSk predicted consequence' is the last column

    last_columns = []
    if 'DeMaSk' in df.columns:
        last_columns.append('DeMaSk')

    columns_to_keep = [col for col in df.columns if col not in last_columns]
    df = df[columns_to_keep + last_columns]

    # Keep only mutations associated to user-defined positions
    if residues:

        # Extract and store in a temporary 'pos' column the
        # positions associated to the mutations
        df['pos'] = df.index.astype(str).str[1:-1]
        log.info(f'Using residues: {residues[0]}')

        # Raise warning if position not in dataframe
        found_residues = False
        for r in residues[0]:
            if r not in df['pos'].to_list():
                log.warning(f'{r} not found in table. Continuing...')
            else:
                found_residues = True

        if not found_residues:
            log.error('ERROR: none of the residues specified with --residues was foundin the input dataframe')
            raise TypeError('ERROR: none of the residues specified with --residues was foundin the input dataframe')

        # Filter dataframe according to the user-defined residues
        df = df[df['pos'].isin(residues[0])]

        # Remove temporary 'pos' column
        df.drop(columns= ['pos'],
                inplace=True)

    # Keep only user-defined mutations
    if mutations:

        # Raise warning if mutation not in dataframe
        found_mutations = False
        for m in mutations[0]:
            if m not in df.index.to_list():
                log.warning(f'{m} not found in table. Continuing...')
            else:
                found_mutations = True


        if not found_mutations:
            log.error('ERROR: none of the mutations specified with --mutations was found in the input dataframe')
            raise TypeError('ERROR: none of the mutations specified with --mutations was found in the input dataframe')

        # Filter dataframe according to the user-defined mutations
        df = df[df.index.isin(mutations[0])]

    # Initialize indexes for additive filtering
    source_idx = pd.Index([])
    clinvar_idx = pd.Index([])

    if plot_Source:
        # Accommodate multiple user options
        pattern = '|'.join(plot_Source)

        # Filter df based on desired mutation source
        source_idx = full_df[full_df['Mutation sources'].str.contains(pattern, case=False, na=False)].index

        # Raise error if filtering leaves df empty
        if source_idx.empty:
            log.warning(f"No mutations found matching the mutation source: {plot_Source}")
        else:
            log.info(f"Found {len(source_idx)} mutations matching the source filter: {plot_Source}")

    # Define an empty df if flags -pltC/-colC are not used
    clinvar_mapped_df = None

    # Assert the required column is present
    if plot_Clinvar or color_Clinvar:

        # Filter df according to poss. previous arguments
        clinvar_mapped_df = full_df.loc[full_df.index.isin(df.index)]
        # Map dictionary to df
        clinvar_mapped_df = map_clinvar_categories(clinvar_mapped_df, clinvar_dict, clinvar_cols)

    # Keep only user-defined mutations according to ClinVar annotation
    if plot_Clinvar:
        # Define a mapping of user-input to categories
        filter_terms = {
            'benign': ['benign', 'benign/likely benign'],
            'likely_benign' : ['likely benign'],
            'pathogenic': ['pathogenic'],
            'likely_pathogenic': ['likely pathogenic'],
            'uncertain': ['uncertain'],
            'conflicting': ['conflicting interpretation'],
            'oncogenic': ['oncogenic'],
            'likely_oncogenic': ['likely oncogenic'],
            'strong': ['strong'],
            'potential': ['potential'],
            'unknown': ['unknown']}

        # Filter dataframe according to flag option
        if 'all' in plot_Clinvar:
            clinvar_idx = clinvar_mapped_df[clinvar_mapped_df['ClinVar Category'].notna()].index
        else:
            # Accommodate multiple user options
            for clinvar_option in plot_Clinvar:
                if clinvar_option in filter_terms:
                    category_values = filter_terms[clinvar_option]
                    # Match categories
                    pattern = '|'.join(category_values)
                    #  Get indexes of mutations matching current ClinVar category
                    indexes_to_add = clinvar_mapped_df[clinvar_mapped_df['ClinVar Category'].str.contains(pattern, case=False, na=False)].index
                    clinvar_idx = clinvar_idx.union(indexes_to_add)
                else:
                    log.warning(f"Unrecognized ClinVar filter option: {clinvar_option}")
                    raise TypeError(f"Unrecognized ClinVar filter option: {clinvar_option}")

        # Raise error if filtering leaves df empty
        if len(clinvar_idx) == 0:
            log.warning(f"No mutations found matching the ClinVar filter {plot_Clinvar}")
        else:
            log.info(f"Found {len(clinvar_idx)} mutations matching the ClinVar filter {plot_Clinvar}")


    # Additively filter if -pltS -pltC
    if plot_Source or plot_Clinvar:
        #Save individual pngs
        combined_indexes = source_idx.union(clinvar_idx)
        df = df[df.index.isin(combined_indexes)]

        log.info(f"Found {len(df)} mutations matching filters.")
        # Raise error if filtering leaves df empty
        if df.empty:
            log.warning("No mutations remain after filtering")

    # Create plotting dataframe and drop columns only used for visualization
    plot_df = df.copy()

    if not plot_Revel and 'REVEL' in plot_df.columns:
        plot_df = plot_df.drop(columns=['REVEL'])

    if not plot_popEVE and 'popEVE' in plot_df.columns:
        plot_df = plot_df.drop(columns=['popEVE'])

    if not plot_Demask and 'DeMaSk predicted consequence' in plot_df.columns:
        plot_df = plot_df.drop(columns=['DeMaSk predicted consequence'])

    if plot_Demask and 'DeMaSk' in plot_df.columns:
        plot_df = plot_df.drop(columns=['DeMaSk'])

    return plot_df, df, full_df, clinvar_mapped_df

def plot(df, full_df, width, height, xlim, clinvar_flag, clinvar_class_type, clinvar_col):
    ''' Plot

    The function is aimed at plotting a matrix showing the effect of
    selected mutations.

    Parameters
    ----------
    df: dataframe
        input MAVISp dataframe
     width: int
        Width of each figure
    height: int
        Height of each figure
    xlim: int
        Number of mutations to include on the x-axis of each plot
    output: str
        Base name for output files
    save_png: bool
        Flag indicating whether to save individual plots as PNG files.
    clinvar_flag: bool
        Flag indicating whether to color xticks according to ClinVar Int. values
    clinvar_col: bool
        Flag indicating whether to color xticks according to ClinVar Int. values

    Returns
    ----------
    figures: list
        list of Figure object, one per plot
    '''
    # Define font of the plot
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = [ 'Arial' ]

    # List axes label colors
    yticklabel_color = color_yticklabels(labels = df.columns.values)

    # Override display label for GEMME row
    display_yticklabels = [
                                'GEMME classification'
                                if lbl == 'GEMME predicted consequence'
                                else 'De novo phosphosite stability'
                                if lbl.startswith('Classification of change in folding free energy with phosphorylation')
                                else 'De novo phosphosite binding'
                                if lbl.startswith('Classification of change in binding free energy with phosphorylation')
                                else lbl
                                for lbl in df.columns.values
                          ]
    # Lower and upper limits
    l = 0
    u = xlim

    expected_plots = math.ceil(df.shape[0]/u)

    figures = []

    if clinvar_flag or clinvar_col:
        height += 0.3

    # Iterate over the number of resulting plots
    for p in range(0, expected_plots):

        # Filter dataframe for the range of rows
        # defined by the upper and lower limits
        filtered_df = df.iloc[l:u, :]

        # Adjust width in figsize for last plot
        if p == range(0, expected_plots)[-1]:

            # If the remainig mutations are less than
            # the chosen xlimit
            if filtered_df.shape[0] != xlim:

                # Adjust width
                width = adjust_figsize_w(x_mut = filtered_df.shape[0],
                                         original_w = width,
                                         original_m = xlim)

        # Adjust height if the number of y-axis labels is large
        if filtered_df.shape[1] >= 15 and height <= 5.9:
            height = 6

        # Convert filtered dataframe to array
        data = filtered_df.to_numpy()

        # List dot colors
        grid_colors = color_markers(data = data)

        # Define grid with the shape of the dataframe
        x = np.linspace(0, data.shape[0], data.shape[0])
        y = np.linspace(0, data.shape[1], data.shape[1])
        grid = np.meshgrid(x, y)

        figsize = tuple((width, height))

        # Define subplots
        fig, ax = plt.subplots(figsize=figsize)

        # Plot scatter
        _ = ax.scatter(x = grid[0],
                    y = grid[1],
                    c = grid_colors,
                    s = 150)

        # Customize plot
        _ = ax.set_xticks(x)
        _ = ax.set_yticks(y)
        _ = ax.set_xticklabels(filtered_df.index, rotation=90)
        _ = ax.set_yticklabels(display_yticklabels, fontweight = 'bold')
        _ = ax.set_xlabel('Mutations')

        # Color ytick labels
        for ticklabel, tickcolor in zip(plt.gca().get_yticklabels(), yticklabel_color):
            ticklabel.set_color(tickcolor)

        # Color xtick labels if flags
        if clinvar_flag or clinvar_col:
            xtick_colors = color_xticklabels(filtered_df.index, full_df['ClinVar Category'])
            for xtick, color in zip(ax.get_xticklabels(), xtick_colors):
                xtick.set_color(color)

        # Adjust padding if less than 10 datapoints
        if data.shape[0] <= 10:
            plt.margins(x = 0.5)

        # Put a legend below current axis
        legend_dict = {'black': 'Damaging',
                    '#d5bdaf': 'Neutral',
                    'lightgray': 'Uncertain',
                    'white': 'Not_available',
                    '#9BB295': 'Early Folding Region',
                    '#8E5466' : 'LoF',
                    '#636D87' : 'GoF'}

        legend_list = [Line2D([0], [0],
                        color='w',
                        marker='o',
                        markersize=10,
                        markerfacecolor=k,
                        label=legend_dict[k],
                        markeredgecolor='k') for k in legend_dict.keys()]

        first_legend = ax.legend(handles = legend_list,
                loc='upper right',
                bbox_to_anchor=(1, 1.2),
                ncol = 8,
                fancybox=True,
                shadow=True)

        # Build second axis + second legend if flags
        if clinvar_flag or clinvar_col:
            ax2 = ax.twinx()
            ax2.get_yaxis().set_visible(False)
            if clinvar_class_type in ("aggregated", "germline"):
                xtick_legend_dict = {
                        '#000000': 'Not available',
                        '#EBD479': 'Conflicting',
                        '#878E99': 'Uncertain ',
                        '#3A8752': 'Benign/Likely Benign',
                        '#A71B32': 'Pathogenic/Likely Pathogenic'}
                if clinvar_class_type == "aggregated":
                    clinvar_label = "ClinVar classification"
                else:
                    clinvar_label = "ClinVar germline classification"
            elif clinvar_class_type == "oncogenicity":
                clinvar_label = "ClinVar Oncogenicity classification"
                xtick_legend_dict = {
                        '#000000': 'Not available',
                        '#878E99': 'Uncertain significance',
                        '#EBD479': 'Conflicting',
                        '#3A8752': 'Benign/Likely benign',
                        '#A71B32': 'Oncogenic/Likely oncogenic'}
            elif clinvar_class_type == "clinical_impact":
                clinvar_label = "ClinVar Clinical Impact"
                xtick_legend_dict = {
                        '#000000': 'Not available',
                        '#EBD479': 'Conflicting',
                        '#A71B32': 'Strong',
                        '#C05A2A': 'Potential',
                        '#878E99': 'Unknown',
                        '#3A8752': 'Benign/Likely benign'}

            xtick_legend_list = [Line2D([0], [0],
                                color='w',
                                marker='o',
                                markersize=10,
                                markerfacecolor=k,
                                label=xtick_legend_dict[k],
                                markeredgecolor='k') for k in xtick_legend_dict.keys()]

            fake = Line2D([], [], color="none")
            labels = [clinvar_label] + list(xtick_legend_dict.values())
            handles = [fake] + xtick_legend_list
            second_legend = ax2.legend(handles=handles,
                                       labels=labels,
                                       loc='lower right',
                                       bbox_to_anchor=(1, -0.4),
                                       ncol=6,
                                       fancybox=True,
                                       shadow=True)

        # Save plot as png and pdf
        fig.tight_layout()

        figures.append(fig)

        # Update upper and lower limits
        l += xlim
        u += xlim

    return figures


def get_damaging_mutations(data):
    ''' Select the mutations with a damaging effect in at least one MAVISp module. '''
    data = data.apply(lambda col: col.map(lambda x: str(x) if isinstance(x, bool) else x))
    data = data.apply(lambda x: x.str.lower() if x.dtype == object else x)

    # Identify the columns of each module

    filter_col_stability = [col for col in data if 'ThermoMPNN' in col
                            or '(Foldetta from FoldX and Rosetta)' in col
                            or '(Foldetta from FoldX and RaSP)' in col
                            or '(Rosetta, FoldX)' in col
                            or '(RaSP, FoldX)' in col
                            or "Classification of change in folding free energy with phosphorylation" in col]
    filter_col_local = [col for col in data if 'Local Int. classification' in col
                        or 'Local Int. With DNA classification' in col
                        or "Classification of change in binding free energy with phosphorylation" in col]
    ptm_stab_clmn = [col for col in data if 'Phosphorylation effect in stability' in col]
    ptm_reg_clmn = [col for col in data if 'Phosphorylation effect in regulation' in col]
    ptm_funct_clmn = [col for col in data if 'Phosphorylation effect in function' in col]
    allosigma_clmn = [col for col in data if 'AlloSigMA 2 predicted consequence' in col]
    long_range_psn_clmn = [col for col in data if 'AlloSigMA2-PSN classification' in col]
    cofactor_active = [col for col in data if 'Functional sites' in col]
    disulfide_clmn = [col for col in data if 'disulfide bridge' in col]

    # Damaging stability
    d_s = ['destabilizing', 'stabilizing']
    all_stability = []
    for col in filter_col_stability:
        all_stability_list = data.index[data[col].isin(d_s)].to_list()
        all_stability.extend(all_stability_list)

    # Damaging local interactions
    all_local = []
    for col in filter_col_local:
        loc_d_list = data.index[data[col].isin(d_s)].to_list()
        all_local.extend(loc_d_list)
    
    # Damaging PTM_stability
    ptm_s_d_list = []
    for col in ptm_stab_clmn:
        ptm_s_d_list.extend(data.index[data[col] == 'damaging'].to_list())

    # Damaging PTM_regulation
    ptm_r_d_list = []
    for col in ptm_reg_clmn:
        ptm_r_d_list.extend(data.index[data[col] == 'damaging'].to_list())

    # Damaging PTM_function
    dmg_str = ['damaging', 'potentially_damaging']
    ptm_f_d_list = []
    for col in ptm_funct_clmn:
        ptm_f_d_list.extend(data.index[data[col].isin(dmg_str)].to_list())

    all_ptm = ptm_s_d_list + ptm_r_d_list + ptm_f_d_list

    # Damaging Long-Range AlloSigma2
    multiple = []
    for col in allosigma_clmn:
        long_d_list = data.index[data[col] == 'destabilizing'].to_list()
        long_me_list = data.index[data[col] == 'mixed_effects'].to_list()
        long_s_list = data.index[data[col] == 'stabilizing'].to_list()
        multiple.extend(long_d_list + long_me_list + long_s_list)

    # Damaging disulfide bridges
    all_disulfide = []
    for col in disulfide_clmn:
        disulfide_d_list = data.index[data[col] == 'damaging'].to_list()
        all_disulfide.extend(disulfide_d_list)

    # Damaging functional sites
    multiple_1 = []
    for col in cofactor_active:
        funct_sites_d_list = data.index[data[col] == 'damaging'].to_list()
        multiple_1.extend(funct_sites_d_list)

    # Damaging Long-Range AlloSigma2-PSN
    psn_d_list = []
    for col in long_range_psn_clmn:
        psn_d_list.extend(data.index[data[col] == 'damaging'].to_list())

    all_long_range = multiple + psn_d_list

    # Combine the damaging mutations of all modules
    damaging_full = set(all_stability + all_local + all_ptm + multiple_1 + all_long_range + all_disulfide)
    full_data_d = data.query('index in @damaging_full')

    return full_data_d

def effect_summary(df):
    '''Add effect column to AM output csv.

    The function takes the df with AM +
    damaging mutations and summarizes the identified
    MAVISp effects in one column as a list, extending
    the df.
    Effects considered: stability, long range,
    local int., functional and ptm.

    Parameters
    ----------
    df: dataframe
        input df
    Returns
    ----------
    df: dataframe
        df with added column of mechanistic indicators
    '''
    # Define code to replace descriptive values
    mechanism_code = {
        '(?i)destabilizing': 1,
        '(?i)pathogenic': 1,
        '(?i)damaging': 1,
        'mixed_effects': 1,
        '(?i)stabilizing': 1,
        '(?i)neutral': 0,
        '(?i)benign': 0,
        '(?i)uncertain': 0,
        'ambiguous': 0}

    # Rename Denovo phospho columns
    df.columns = df.columns.str.replace("Classification of change in folding free energy with phosphorylation", "De novo phosphorylation stability", regex=False)
    df.columns = df.columns.str.replace("Classification of change in binding free energy with phosphorylation", "De novo phosphorylation binding", regex=False)

    # Replace values with 0/1s
    temp_df = df.replace(mechanism_code, regex=True)

    # Wrangle disulfide module columns to make just one
    if 'Loss of disulfide bridge' and 'Predicted de-novo disulfide bridge' in temp_df.columns:
        temp_df['Disulfide bridges'] = temp_df['Loss of disulfide bridge'] + temp_df['Predicted de-novo disulfide bridge']
        temp_df['Disulfide bridges'] = temp_df['Disulfide bridges'].clip(upper=1)
        temp_df = temp_df.drop(columns=['Loss of disulfide bridge', 'Predicted de-novo disulfide bridge'])

    # Define broad effect categories + corresponding regex patterns
    effect_categories = {
        'Stability': 'Stability classification|De novo phosphorylation stability',
        'Local Int.': 'Local Int.|De novo phosphorylation binding',
        'PTM': 'Phosphorylation effect',
        'Long Range': 'AlloSigMA',
        'Functional': 'Functional sites',
        'Disulfide bridges': 'Disulfide bridges'}

    # Initialise series for storing results
    effects_summary = pd.Series(index=df.index, dtype='object')

    # Iterate over each mutation in df
    for idx, row in temp_df.iterrows():
        effects = []
        # For each effect cattegory check if effect found
        # i.e. value is 1
        for category, pattern in effect_categories.items():
            if pd.to_numeric(row.filter(regex=pattern), errors="coerce").max() == 1:
                effects.append(category)
        effects_summary.at[idx] = ', '.join(effects)

    # Add the new column to the original DataFrame
    df['MAVISp Effects'] = effects_summary

    return df

def load_clinvar_dict(tsv_file):
    """Load the ClinVar dictionary.

    The function aims to read a tab-delimited
    file containing two columns, delineating
    ClinVar annotations and their respective internal
    categories.

    Parameters
    ----------
    tsv_file : str
        input clinvar dictionary file
    Returns
    -------
    clinvar_dict : dict
        dictionary of ClinVar annotations + internal categories
    """
    clinvar_dict = pd.read_csv(tsv_file,
                                sep='\t',
                                header=None,
                                names=['clinvar', 'internal_category'])
    return clinvar_dict.set_index('clinvar')['internal_category'].to_dict()

def map_clinvar_categories(dataframe, clinvar_dict, clinvar_cols):
    """Translate ClinVar to internal categories.

    The function takes the full dataframe of the input
    CSV and checks whether a given mutation was sourced
    from ClinVar, if so the function translates the
    ClinVar annotation to an internal nomenclature
    in a new column.

    Parameters
    ----------
    dataframe : dataframe
        Full dataframe containing 'ClinVar Interpretation' column.
    clinvar_dict : dict
        Dictionary mapping ClinVar annotations to categories.
    Returns
    -------
    dataframe : dataframe
        Updated dataframe with added 'ClinVar Category' column.
    """
    dataframe['ClinVar Category'] = dataframe.apply(
    lambda row: clinvar_dict.get(str(row[clinvar_cols[0]]).replace(", ", ","))
    if (pd.notna(row['Mutation sources'])
        and 'clinvar' in row['Mutation sources'].lower()
        and pd.notna(row[clinvar_cols[0]]))
    else None,
    axis=1)
    return dataframe

def filter_vep_summary(summary, df, vep_filter, glof_filter):
    """Filter the AlphafoldMissense summary dataframe

    This function filters the input dataframe to keep only columns of interest
    with a classification, and that are present in the original input.

    Parameters
    ----------
    summary : dataframe
        summary data frame
    df : dataframe
        full data frame
    vep_filter: str
        which VEP to use for filtering for pathogenic variants. "none" implies no filter
    glof_filter: bool
        whether to filter for cases that have DeMaSk LoF/GoF assignment

    Returns
    -------
    dataframe : dataframe
        filtered dataframe for alphamissense summary
    """

    f = lambda x: 'Stability classification' in x or \
              'Local Int. classification' in x or \
              'Local Int. With DNA classification' in x or \
              "Classification of change in folding free energy with phosphorylation" in x or \
              "Classification of change in binding free energy with phosphorylation" in x or \
              'AlloSigMA 2 predicted consequence' in x or \
              'AlloSigMA2-PSN classification' in x or \
              'Phosphorylation effect' in x or \
              'Functional sites' in x or \
              'Loss of disulfide bridge' in x or \
              'Predicted de-novo disulfide bridge' in x or \
              ('Mutation' in x and 'Mutation sources' not in x)

    # Filter the columns to keep only classification info
    filter_columns = [col for col in summary.columns if f(col)]

    # Create a new df with the selected columns
    filtered_am = summary[filter_columns]

    # Filter output according to previously applied args
    filtered_am = filtered_am[filtered_am.index.isin(df.index)]

    if vep_filter == 'alphamissense':
        filtered_index_vep = df[df['AlphaMissense classification'] == 1].index
    elif vep_filter == 'gemme':
        filtered_index_vep = df[df['GEMME predicted consequence'] == 1].index
    elif vep_filter == 'revel':
        filtered_index_vep = df[df['REVEL'] == 1].index
    elif vep_filter == 'eve':
        filtered_index_vep = df[df['EVE classification (25% Uncertain)'] == 1].index
    elif vep_filter == 'popeve':
        filtered_index_vep = df[df['popEVE'] == 1].index
    elif vep_filter == 'none':
        filtered_index_vep = df.index

    if glof_filter:
        filtered_index_glof = df[df['DeMaSk predicted consequence'].isin([7, 8])].index
    else:
        filtered_index_glof = df.index

    filtered_index = list(set(filtered_index_vep).intersection(set(filtered_index_glof)))

    filtered_am = filtered_am[filtered_am.index.isin(filtered_index)]
    
    # Add new column with identified broad effects
    filtered_am = effect_summary(filtered_am)

    if filtered_am.empty:
        log.warning(f'No mutations found with AM pathogenic annotation'\
                    f' output file will be empty...')


    return filtered_am


def main():

    # Basic logging configuration
    log.basicConfig(level=log.INFO, \
                    format='%(levelname)s - %(message)s')

    # Add arguments required to the script to a argparse.ArgumentParser instance.
    description = "Plot of all the analyses performed on the target " \
                  "protein to associate a mutation to the predicted " \
                  "damaging/neutral/uncertain effect. "
    parser = argparse.ArgumentParser(description = description)

    i_helpstr = "Input: MAVISp aggregated csv."
    parser.add_argument("-i", "--input",
                        action = "store",
                        type = str,
                        help = i_helpstr,
                        required = True)

    v_helpstr = "Input: MAVISp ClinVar dictionary file. Only required if options -pltC or -colC are used." \
                "for old entries or new entries with germline classification → use the standard  dictionary" \
                "for oncogenicity  → use the onco dictionary" \
                "for clinical impact → use the clinical dictionary"
    parser.add_argument("-v", "--clinvar-dictionary",
                        dest = "clinvar_dict",
                        action = "store",
                        type = str,
                        help = v_helpstr,
                        default = "dictionary.csv")

    o_default = 'dot_plot'
    o_helpstr = f"Output filename. Default = {o_default}"
    parser.add_argument("-o","--output",
                        type = str,
                        default = o_default,
                        help = o_helpstr)

    m_helpstr = "Mutations of interest. Please provide " \
                "comma-separated mutations (e.g., -m M1A,K77R)"
    parser.add_argument("-m","--mutations",
                        nargs = '+',
                        type = lambda arg: arg.strip().split(','),
                        help = m_helpstr)

    r_helpstr = "Residue numbers of interest. Please provide " \
                "comma-separated residues (e.g., -r 1,5,6,18)"
    parser.add_argument("-r","--residues",
                        nargs = '+',
                        type= lambda arg: arg.strip().split(','),
                        help = r_helpstr)

    R_default = 0.5
    R_helpstr = f"Threshold to classify a mutation according to the " \
                f"revel score. (Default = {R_default})"
    parser.add_argument("-R","--revel_threshold",
                        default = R_default,
	                    type = float,
                        help = R_helpstr)

    P_default = -4.617
    P_helpstr = f"Threshold to classify a mutation according to the " \
                f"popEVE score. (Default = {P_default})"
    parser.add_argument("-P", "--popeve_threshold",
                        default=P_default,
                        type=float,
                        help=P_helpstr)

    D_default = 0.25
    D_helpstr = f"Threshold to classify a mutation according to the " \
                f"DeMask score. (Default = {D_default})"
    parser.add_argument("-D","--demask_threshold",
                        default = D_default,
	                    type = float,
                        help = D_helpstr)

    G_default = -3.0
    G_helpstr = f"Threshold to classify a mutation according to the " \
                f"GEMME  score. (Default = {G_default})"
    parser.add_argument("-G","--gemme_threshold",
                        default = G_default,
	                    type = float,
                        help = G_helpstr)

    x_default = 50
    x_helpstr = f"Number of mutations to plot on the x axis. " \
                f"(Default = {x_default})"
    parser.add_argument("-x","--x_lim",
                        default = x_default,
                        type = int,
                        help = x_helpstr)

    f_default = [14, 5]
    f_helpstr = f"Figure size of each plot. It is suggested to use the default " \
                f"with 40/50 mutations on the x-axis and 7/8 labels on the " \
                f"y-axis. (Default = {f_default})"
    parser.add_argument("-f","--figsize",
                        default = f_default,
                        nargs = 2,
                        type = float,
                        help = f_helpstr)

    pltR_helpstr = f"Plotting of Revel scores. " \
                   f"(Default = None)"
    parser.add_argument("-pltR","--plot_Revel",
                        action = 'store_true',
                        help = pltR_helpstr)

    pltP_helpstr = f"Plotting of popEVE classification. (Default = None)"
    parser.add_argument("-pltP", "--plot_popEVE",
                        action='store_true',
                        help=pltP_helpstr)

    pltD_helpstr = f"Plotting of Demask LoF/GoF if" \
                    f" mutation is above demask threshold. " \
                    f"(Default = None)"
    parser.add_argument("-pltD", "--plot_Demask",
                        action = 'store_true',
                        help = pltD_helpstr)

    clinvarclasstype_helpstr = f"ClinVar classification type to use for plotting/coloring. " \
                               f"Choices:aggregated, germline, oncogenicity, clinical_impact."
    parser.add_argument("-cct", "--clinvar_class_type",
                        default = "aggregated",
                        choices = ["aggregated", "germline", "oncogenicity", "clinical_impact"],
                        help = clinvarclasstype_helpstr)

    pltC_helpstr =  f"Plotting of Clinvar variants " \
                    f"choose from: all, uncertain, " \
                    f"benign, likely benign, pathogenic, " \
                    f"likely pathogenic, conflicting, oncogenic, " \
                    f"likely_oncogenic, potential, strong, unknown " \
                    f"tier_IV. For combinations, provide space separated " \
                    f"options (e.g. benign uncertain). Make sure that the " \
                    f"chosen options are aligned with the chosen cct flag." \
                    f"(Default = None)"
    parser.add_argument("-pltC", "--plot_Clinvar",
                        nargs='+',
                        choices = ["all",
                                    "uncertain",
                                    "benign",
                                    "likely_benign",
                                    "pathogenic",
                                    "likely_pathogenic",
                                    "conflicting",
                                    "oncogenic",
                                    "likely_oncogenic",
                                    "potential",
                                    "strong",
                                    "unknown"],
                        help = pltC_helpstr)

    colC_helpstr =  f"Color x-axis of ClinVar mutations" \
                    f"(Default = None)"
    parser.add_argument("-colC", "--color_Clinvar",
                        action = 'store_true',
                        help = colC_helpstr)

    pltS_helpstr =  f"Plotting of mutations based on source " \
                    f"choose from: saturation, cosmic, cbioportal. " \
                    f"For combinations, provide space separated options " \
                    f"(e.g. saturation cosmic). Default = all."
    parser.add_argument("-pltS", "--plot_Source",
                        nargs = "+",
                        choices = ["saturation",
                                   "cosmic",
                                   "cbioportal"],
                        help = pltS_helpstr)

    AMx_helpstr = "Restrict mechanisitc indicators output to pathogenic variants. Choose the VEP to use" \
                  "to detect pathogenic variants between none, alphamissense, revel, gemme, eve, popeve. If this option" \
                  "is used without argument it will default to alphamissense"
    parser.add_argument("-vep", "--vep-filter",
                        choices=["none", "alphamissense", "revel", "gemme", "eve", "popeve"],
                        nargs="?",
                        const="alphamissense",
                        default="none",
                        help = AMx_helpstr,
                        dest='AMx')
    
    lgof_helpstr = "Restrict mechanisitc indicators output to variants that are classified as LoF or GoF for DeMaSk"
    parser.add_argument("-lgof", "--vep-filter-lgof",
                        action = 'store_true',
                        help = lgof_helpstr,
                        dest='adf',
                        default=False)
                        
    # Enable shell autocompletion when available (safe guard)
    try:
        argcomplete.autocomplete(parser)
    except Exception:
        pass
    args = parser.parse_args()

    # Assert the user didn't provide both -r and -m options
    if args.residues and args.mutations:
        log.error("Please choose only one option between user-defined " \
                  "positions (-r) and mutations (-m). Exiting...")
        sys.exit(1)

    # If user provided custom residues,
    if args.residues:
        # Assert length of the argument; must be 1
        # as expecting comma-separated values
        if len(args.residues) > 1:
            log.error('Please provide comma-separated residues. Exiting...')
            sys.exit(1)

    # If user provided custom mutations,
    if args.mutations:
        # Assert length of the argument; must be 1
        # as expecting comma-separated values
        if len(args.mutations) > 1:
            log.error('Please provide comma-separated mutations. Exiting...')
            sys.exit(1)
        log.info(f'Using mutations: {args.mutations[0]}')

    # Assert the user didn't provide both -pltC and -colC options
    if args.plot_Clinvar and args.color_Clinvar:
        log.error("Please choose only one option between " \
                  "plotting specific ClinVar variants (-pltC) and " \
                  "coloring ClinVar variants (-colC) among all mutations. " \
                  "Exiting...")
        sys.exit(1)

    save_png = bool(args.residues) or bool(args.mutations) or bool(args.plot_Source) or bool(args.plot_Clinvar)

    ################################# INPUT #################################

    # Read input

    full_df = pd.read_csv(args.input, index_col = 'Mutation')
    try:
        clinvar_cols = get_clinvar_columns(full_df, args.clinvar_class_type)
    except ValueError as e:
        log.error(f"ERROR: {e}")
        exit(1)

    if args.plot_Clinvar or args.color_Clinvar:
        try:
            # Load the ClinVar dictionary
            clinvar_dict = load_clinvar_dict(args.clinvar_dict)
        except FileNotFoundError as e:
            log.error(f"ClinVar dictionary file not found: {e}. Exiting...")
            exit(1)
        # check if correct dict was provided:
        expected_header = {"aggregated": "#ClinVar",
                           "germline": "#ClinVar",
                           "oncogenicity": "#Oncogenicity",
                           "clinical_impact": "#Clinical Impact"}[args.clinvar_class_type]
        first_cell = pd.read_csv(args.clinvar_dict, sep="\t", header=None).iat[0, 0]
        if str(first_cell).strip() != expected_header:
            log.error(f"ERROR: Wrong ClinVar dictionary for --clinvar_class_type {args.clinvar_class_type}")
            exit(1)
    else:
        clinvar_dict = None

    try:
        plot_df, classification_df, dataframe, clinvar_mapped_df =  process_input(full_df = full_df,
                                                          r_cutoff = args.revel_threshold,
                                                          p_cutoff=args.popeve_threshold,
                                                          d_cutoff = args.demask_threshold,
                                                          g_cutoff= args.gemme_threshold,
                                                          residues = args.residues,
                                                          mutations = args.mutations,
                                                          clinvar_dict = clinvar_dict,
                                                          plot_Revel = args.plot_Revel,
                                                          plot_popEVE=args.plot_popEVE,
                                                          plot_Demask = args.plot_Demask,
                                                          plot_Source = args.plot_Source,
                                                          plot_Clinvar = args.plot_Clinvar,
                                                          color_Clinvar = args.color_Clinvar,
                                                          clinvar_cols = clinvar_cols)
    except TypeError as e:
        print(f"ERROR: {str(e)}")
        exit(1)
    except FileNotFoundError as e:
        print(f"ERROR: {str(e)}")
        exit(1)

    ############################### SUMMARY ###############################

    # Select the mutations with a damaging effect in at least one module
    summary_df = get_damaging_mutations(data = dataframe)
    
    ################################# PLOT #################################
    '''
    filter_col = [col for col in plot_df if col.startswith("Local Int.")]
    try:
        func_sites = df[['Functional sites (cofactor)', 'Functional sites (active site)']]
        tmp_df = df.drop(columns=['Functional sites (cofactor)', 'Functional sites (active site)'])
        df = pd.concat([func_sites, tmp_df], axis=1)
        new_order = ['Functional sites (cofactor)',
                    'Functional sites (active site)'] + [col for col in tmp_df.columns]
        df = df[new_order]
    except:
        pass
    '''

    # Plot dot plot
    figures = plot(df = plot_df,
                   full_df = clinvar_mapped_df,
                   width = args.figsize[0],
                   height = args.figsize[1],
                   xlim = args.x_lim,
                   clinvar_flag = bool(args.plot_Clinvar),
                   clinvar_col = bool(args.color_Clinvar),
                   clinvar_class_type = args.clinvar_class_type)


    if len(figures) > 0:
        with PdfPages(f"{args.output}.pdf") as pdf:
            for i, figure in enumerate(figures):
                figure.savefig(pdf, format='pdf', dpi=300)
                if save_png:
                    figure.savefig(f'{args.output}_{i}.png', dpi=300)
    else:
        log.warning("No figures available to produce; figure output will not be written")

################################# AM CSV #################################
    filtered_am = filter_vep_summary(summary_df, classification_df, args.AMx, args.adf)

    filtered_am.to_csv(f'mechanistic_indicators_out_{args.clinvar_class_type}.csv', index=True)

if __name__ == '__main__':
    main()
