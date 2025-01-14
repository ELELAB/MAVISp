#!/usr/bin/env python

# Copyright (C) 2023 Ludovica Beltrame <beltrameludo@gmail.com>,
# Simone Scrima <simonescrima@gmail.com>, Karolina Krzesińska <kzokr@dtu.dk>,
# Matteo Tiberti
# Danish Cancer Society & Technical University of Denmark

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
from itertools import combinations
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
                 'PTM effect in stability',
                 'EFoldMine - part of early folding region']
    
    function = ['Local Int.',
                'Functional sites (cofactor)',
                'Functional sites (active site)',
                'PTM',
                'AlloSigma2']

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
        'uncertain': '#878E99'}

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

def process_input(full_df, r_cutoff, d_cutoff, g_cutoff):
    ''' Read MAVISp aggregated table.

    The function takes as input a MAVISp csv file and returns
    a dataframe formatted in a plot-compatible way.

    Parameters
    ----------
    data: dataframe
	   input MAVISp dataframe
    r_cutoff: float
        REVEL score cutoff
    Returns
    ----------
    df: dataframe
        output dataframe in which the classification is
        converted to number to facilitate the plotting.
    full_df: dataframe
              mavisp csv without altered columns

    '''
    f = lambda x: '(Rosetta, FoldX)' in x or \
                    '(RaSP, FoldX)' in x or \
                    'Local Int. classification' in x or \
                    'Local Int. With DNA classification' in x or \
                    'Functional sites (cofactor)' in x or \
                    'Functional sites (active site)' in x or \
                    'AlloSigma2 predicted consequence - active sites' in x or \
                    'AlloSigma2 predicted consequence - cofactor sites' in x or \
                    'AlloSigma2 predicted consequence - pockets and interfaces' in x or \
                    ('AlloSigma2-PSN classification' in x and not 'AlloSigma2 mutation type' in x) or\
                    'PTM effect in ' in x or 'REVEL score' in x or \
                    'EVE classification (25% Uncertain)' in x or \
                    'DeMaSk delta fitness' in x or \
                    'DeMaSk predicted consequence' in x or \
                    ('GEMME Score' in x and not 'GEMME Score (rank-normalized)' in x) or \
                    'AlphaMissense classification' in x or \
                    ('Mutation' in x and not 'Mutation sources' in x and not 'Mutation predicted to' in x) or \
                    'EFoldMine - part of early folding region' in x or \
                    'Experimental data classification' in x

    df = full_df.copy()
    filter_columns = [col for col in full_df.columns if f(col)]
    df = df[filter_columns]

    # If pltRevel is True then convert REVEL score column to float
    # In case multiple REVEL scores are present it returns the average
    # of the values

    for d in [df, full_df]:
        d['REVEL score'] = d['REVEL score'].apply(convert_to_float)


    # Add REVEL score interpretation column
    df['REVEL'] = np.where(df['REVEL score'].isna(), None,
                                np.where(df['REVEL score'] >= r_cutoff,
                                    'Damaging', 'Neutral'))

    try:
        # Convert GEMME score into absolute value
        df['GEMME predicted consequence'] = np.where(df['GEMME Score'].isna(), None,
                np.where(df['GEMME Score'] >= g_cutoff,'gain_of_function',
                np.where(df['GEMME Score'] <= -g_cutoff,'loss_of_function',
                'Neutral')))
    except:
        log.warning(f'- no GEMME found in MAVISp csv.')

    # Convert Demask delta fitness into absolute value
    df['DeMaSk'] = np.where(df['DeMaSk delta fitness'].isna(), None,
                np.where(df['DeMaSk delta fitness'] >= d_cutoff,'Damaging',
                np.where(df['DeMaSk delta fitness'] <= -d_cutoff,'Damaging',
                'Neutral')))

    # Only keep Demask consequence for those mutations that satisfy the Demask threshold
    df['DeMaSk predicted consequence'] = np.where(df['DeMaSk'].isna(), None,
            np.where(df['DeMaSk'] != 'Damaging', 'Neutral', df['DeMaSk predicted consequence']))

    # Drop score columns
    if 'GEMME Score' in df.columns:
        df.drop(columns = ['REVEL score','DeMaSk delta fitness', 'GEMME Score'],
            inplace = True)
    else:
        df.drop(columns = ['REVEL score','DeMaSk delta fitness'],
            inplace = True)
   
    # Sort columns based on broad effect categories
    functional_cols = [col for col in df.columns if 'functional' in col.lower() and 'experimental data classification' not in col.lower()]
    stability_cols = [col for col in df.columns if 'stability' in col.lower() and 'experimental data classification' not in col.lower()]
    efold_col = [col for col in df.columns if 'efoldmine' in col.lower()]
    other_cols = [col for col in df.columns if 'functional' not in col.lower() and 'stability' not in col.lower() and 'experimental data classification' not in col.lower() and 'efoldmine' not in col.lower()]
    experimental_cols = list(set([col for col in df.columns if 'Experimental data classification' in col]))
    
    # Combine lists in order
    df = df[stability_cols + efold_col + functional_cols + other_cols + experimental_cols]

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
                   'gain_of_function' : 8,
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
    if 'DeMaSk predicted consequence' in df.columns:
        last_columns.append('DeMaSk predicted consequence')
    if 'GEMME predicted consequence' in df.columns:
        last_columns.append('GEMME predicted consequence')

    # Filter out the last columns from the main columns list
    columns_to_keep = [col for col in df.columns if col not in last_columns]

    # Append the last columns to the final list
    df = df[columns_to_keep + last_columns]
    
    return df, full_df

def plot(df, full_df, width, height, xlim, clinvar_flag, clinvar_col):
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
        _ = ax.set_yticklabels(df.columns.values, fontweight = 'bold')
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
                    '#8E5466' : 'Loss of Function',
                    '#636D87' : 'Gain of Function'}

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

            xtick_legend_dict = {
                    '#EBD479': 'Conflicting',
                    '#878E99': 'Uncertain ',
                    '#3A8752': 'Benign/Likely Benign',
                    '#A71B32': 'Pathogenic/Likely Pathogenic'}

            xtick_legend_list = [Line2D([0], [0],
                                color='w',
                                marker='o',
                                markersize=10,
                                markerfacecolor=k,
                                label=xtick_legend_dict[k],
                                markeredgecolor='k') for k in xtick_legend_dict.keys()]

            second_legend = ax2.legend(handles=xtick_legend_list,
            loc='upper right',
            bbox_to_anchor=(1, 1.11),
            ncol=4,
            fancybox=True,
            shadow=True)

        # Save plot as png and pdf
        fig.tight_layout()

        figures.append(fig)

        # Update upper and lower limits
        l += xlim
        u += xlim

    return figures

def generate_summary(data,d_cutoff,r_cutoff):
    ''' Summary log.txt file.

    The function is aimed at summarizing the number of mutations
    with uncertain or damaging consequences for each type of
    effect (structural and/or functional).
    Additionally, it identified mutations with a mavisp
    effect and with AM pathogenic classification and return
    them as a df.

    Parameters
    ----------
    data: dataframe matrix
    out: output file to be written
    Returns
    ----------
    out: str
        string containing the summary information
    alpha_d_df: df
        df containing damaging mutations with AM pathogenic class.
    '''

    # Convert booleans to strings and then all the strings to lowercase
    data = data.apply(lambda col: col.map(lambda x: str(x) if isinstance(x, bool) else x))
    data = data.apply(lambda x: x.str.lower() if x.dtype == object else x)

    # Colnames ptm and allosigma

    filter_col_stability = [col for col in data if '(Rosetta, FoldX)' in col or '(RaSP, FoldX)' in col]
    filter_col_local = [col for col in data if 'Local Int. classification' in col or 'Local Int. With DNA classification' in col]
    ptm_stab_clmn = [col for col in data if 'PTM effect in stability' in col]
    ptm_reg_clmn = [col for col in data if 'PTM effect in regulation' in col]
    ptm_funct_clmn = [col for col in data if 'PTM effect in function' in col]
    allosigma_clmn = [col for col in data if 'AlloSigma2 predicted consequence' in col]
    long_range_psn_clmn = [col for col in data if 'AlloSigma2-PSN classification' in col]
    cofactor_active = [col for col in data if 'Functional sites' in col]
    # # Check if columns with well-defined names are in the dataframe
    # for col in [ptm_stab_clmn, ptm_reg_clmn, ptm_funct_clmn, allosigma_clmn]:
    #     if col not in data.columns:
    #         log.warning(f'{col} not found in csv. Skipping associated analyses for summary file.')

    out = ''
    out += 'Summary of the effects.\n\n'

    ### UNCERTAIN ###
    out += '>>> UNCERTAIN effects:\n\n'
    pattern = r'\[.*?\]'

    ### Stability ###
    for col in filter_col_stability:
        try:
            ensemble = re.findall(pattern, col)

            # Calculate number of uncertain mutations
            stab_un = str(len(data[data[col] == 'uncertain']))

            # List uncertain mutations
            stab_un_list = data.index[data[col] == 'uncertain'].to_list()

            # Print list of mutations only if they are less than 10
            mutlist = out_list(stab_un_list)

            if ensemble:
                out += f'- {stab_un} variants remain uncertain for {col}. {mutlist}'

                #out += f'- {stab_un} variants remain uncertain for STABILITY {ensemble[0]}. {mutlist}'
            else:
                out += f'- {stab_un} variants remain uncertain for {col}. {mutlist}'
                #out += f'- {stab_un} variants remain uncertain for STABILITY. {mutlist}'

        except KeyError:
            log.warning(f'None of the STABILITY modules were found in MAVISp csv')

    ### Local interactions ###
    interactors = r'\((.*?)\)' # interactor are delimited by ()
    if filter_col_local:
        # Iterate over the interactors
        for col in filter_col_local:
            try:
                ensemble = re.findall(pattern, col)
                inter = re.findall(interactors, col) # search for interactors

                # Calculate number or uncertain mutations
                loc_un = str(len(data[data[col] == 'uncertain']))

                # List uncertain mutations
                loc_un_list = data.index[data[col] == 'uncertain'].to_list()

                # Print list of mutations only if the are less than 10
                mutlist = out_list(loc_un_list)

                if ensemble:
                    if 'DNA' in col:
                        out += f'- {loc_un} variants remain uncertain for LOCAL_INTERACTION with DNA {inter}. {ensemble[0]}. {mutlist}'
                        continue
                    if inter:
                        out += f'- {loc_un} variants remain uncertain for LOCAL_INTERACTION {inter} {ensemble[0]}. {mutlist}'
                else:
                    if 'DNA' in col:
                        out += f'- {loc_un} variants remain uncertain for LOCAL_INTERACTION with DNA {inter}. {mutlist}'
                        continue
                    if inter:
                        out += f'- {loc_un} variants remain uncertain for LOCAL_INTERACTION {inter}. {mutlist}'
            except KeyError:
                log.warning(f'{col} not found in MAVISp csv.')
    else:
        out += '- No classification for LOCAL_INTERACTION available yet.\n'

    ### PTM_stability ###
    for col in ptm_stab_clmn:
        try:
            ensemble = re.findall(pattern, col)

            # Calculate number of uncertain mutations
            ptm_s_un = str(len(data[data[col] == 'uncertain']))

            # List uncertain mutations
            ptm_s_un_list = data.index[data[col] == 'uncertain'].to_list()

            # Print list of mutations only if the are less than 10
            ptm_mutlist_stab =  out_list(ptm_s_un_list)

            if ensemble:
                out += f'- {ptm_s_un} variants remain uncertain for PTM-stability. {ensemble[0]} {ptm_mutlist_stab} '
            else:
                out += f'- {ptm_s_un} variants remain uncertain for PTM-stability. {ptm_mutlist_stab} '
        except KeyError:
            log.warning(f'{ptm_stab_clmn} not found in MAVISp csv.')

    ### PTM_regulation ### 
    for col in ptm_reg_clmn:
        try:
            ensemble = re.findall(pattern, col)

            # Calculate number of uncertain mutations
            ptm_r_un = str(len(data[data[col] == 'uncertain']))

            # List uncertain mutations
            ptm_r_un_list = data.index[data[col] == 'uncertain'].to_list()

            # Print list of mutations only if the are less than 10
            ptm_mutlist_reg = out_list(ptm_r_un_list)
            if ensemble:
                out += f'- {ptm_r_un} variants remain uncertain for PTM-regulation. {ensemble[0]} {ptm_mutlist_reg} '
            else:
                out += f'- {ptm_r_un} variants remain uncertain for PTM-regulation. {ptm_mutlist_reg} '

        except KeyError:
            log.warning(f'{ptm_reg_clmn} not found in MAVISp csv.')

    # Uncertain PTM_function
    for col in ptm_funct_clmn:
        try:
            ensemble = re.findall(pattern, col)

            # Calculate number of uncertain mutations
            ptm_f_un = str(len(data[data[col] == 'uncertain']))

            # List uncertain mutations
            ptm_f_un_list = data.index[data[col] == 'uncertain'].to_list()

            # Print list of mutations only if the are less than 10
            ptm_mutlist_funct = out_list(ptm_f_un_list)
            if ensemble:
                out += f'- {ptm_f_un} variants remain uncertain for PTM-function. {ensemble[0]} {ptm_mutlist_funct} '
            else:
                out += f'- {ptm_f_un} variants remain uncertain for PTM-function. {ptm_mutlist_funct} '
        except KeyError:
            log.warning(f'- {ptm_funct_clmn} not found in MAVISp csv.')

    ### Assess presence of PTM module ####
    # Start counter
    count_ptm_na = 0
    d = {}
    # Iterate over the ptm columns
    ptm_module = ptm_stab_clmn + ptm_reg_clmn + ptm_funct_clmn

    # Create dict for counting ptm per ensemble
    for ptm in ptm_module:
        matches = re.findall(pattern, ptm)
        if matches:
            d[matches[0]] = 0

    for col in ptm_module:
        ensemble = re.findall(pattern, col)
        # If column not in dataframe
        if ensemble:
            for ens in d:
                if col not in data.columns and ensemble[0] == ens:
                    d[ens] += 1
                else:
                    continue
        else:
            if col not in data.columns:
                # Increase counter by 1
                count_ptm_na += 1
    # Per ensemble if none is in the table
    for ensemble in d:
        if d[ensemble] == 3:
             out += f'- No classification for PTM available yet for the {ensemble} ensemble.' \
                     'No phosphorylation sites have been found in correspondence of the mutation sites.\n'

    # If none of them in the table
    if count_ptm_na == 3:
        out += '- No classification for PTM available yet. No phosphorylation sites ' \
               'have been found in correspondence of the mutation sites.\n'

    lr_type = r'-(.*)'
    ### Not evaluated long-range Allosigma2 ###
    for col in allosigma_clmn:
        try:
            lr = re.findall(lr_type, col)

            # Calculate number of uncertain mutations
            na_long = str(len(data[data[col] == 'uncertain']))
            out += f'- {na_long} variants couldn\'t be assessed with the LONG_RANGE {lr} module ' \
                   f'due to small changes in side-chain volume upon mutation. '

            # List uncertain mutations
            allosigma_un_list = data.index[data[col] == 'uncertain'].to_list()

            # Print list of mutations only if the are less than 10
            out += out_list(allosigma_un_list)
        except KeyError:
            log.warning(f'- {col} not found in MAVISp csv.')

    funct_sites = r'\((.*?)\)'
    ### Not evaluated long-range Allosigma2 ###
    for col in cofactor_active:
        try:
            site = re.findall(funct_sites, col)
            # Calculate number of uncertain mutations
            site_unc = str(len(data[data[col] == 'uncertain']))
            out += f'- {site_unc} variants couldn\'t be assessed with the functional site {site} module ' \
            # List uncertain mutations
            func_un_list = data.index[data[col] == 'uncertain'].to_list()
            # Print list of mutations only if the are less than 10
            out += out_list(func_un_list)
        except KeyError:
            log.warning(f'- {col} not found in MAVISp csv.')

    ### Not evaluated for Long-Range AlloSigma2-PSN ###
    if long_range_psn_clmn:
        for col in long_range_psn_clmn:
            try:
                ensemble = re.findall(pattern, col)

                # Calculate number of uncertain mutations
                psn_un = str(len(data[data[col] == 'uncertain']))

                # List uncertain mutations
                psn_un_list = data.index[data[col] == 'uncertain'].to_list()

                mutlist = out_list(psn_un_list)

                # Print list of mutations only if the are less than 10
                if ensemble:
                    out += f'- {psn_un} variants remain uncertain for Long Range-PSN. {ensemble[0]} {mutlist} '
                else:
                    out += f'- {psn_un} variants remain uncertain for Long Range-PSN. {mutlist} '
            except:
                continue
    else:
        log.warning(f'- Long Range-PSN not found in MAVISp csv.')
        out += f'- No classification for Long Range-PSN available yet.\n'

    ### Not classified by any module ###

    # Count mutations that cannot be classified
    data['is_NC'] = data.apply(lambda row: all(value in ['uncertain', 'unknown', np.nan, None] for value in row), axis=1)
    count_NC = data['is_NC'].sum()
    out += f'- {count_NC} variants could not be classified with any of the MAVISp modules.'

    # List mutations that cannot be classified
    NC_list = data.index[data['is_NC'] == True].to_list()

    # Print list of mutations only if the are less than 10
    out += out_list(NC_list)

    ### DAMAGING ###
    out += '\n>>> DAMAGING/DESTABILIZING effects:\n\n'

    ### Damaging stability ###
    all_stability = []

    # Define columns and effects' variables
    clinvar_clmn = 'ClinVar Interpretation'
    vus = '(?i)uncertain'
    conflict = 'conflicting interpretations of pathogenicity'
    d_s = ['destabilizing', 'stabilizing']

    if filter_col_stability:
        for col in filter_col_stability:
            try:
                # Calculate number of damaging mutations
                stab_d = str(len(data[data[col].isin(d_s)]))

                # List damaging mutations
                all_stability_list = data.index[data[col].isin(d_s)].to_list()
                all_stability.extend(all_stability_list)

                # Print list of mutations only if the are less than 10
                dmg_mutlist = out_list(all_stability_list)

                out += f'- {stab_d} variants are predicted with damaging effects for {col}. {dmg_mutlist}'

                try:
                    # Calculate number of VUS
                    stab_d_vus = str(len(data[(data[col].isin(d_s) & \
                                        data[clinvar_clmn].astype(str).str.contains(vus))]))
                    out += f'-- {stab_d_vus} of these are vus: '

                    # List VUS mutations
                    vus_list = data.index[(data[col].isin(d_s) & \
                                        data[clinvar_clmn].astype(str).str.contains(vus))].to_list()

                    # Print list of mutations only if the are less than 10
                    out += out_list(vus_list)

                    # Calculate number of variants with conflicting interpretations
                    stab_d_conf = str(len(data[(data[col].isin(d_s) & \
                                        data[clinvar_clmn].astype(str).str.contains(conflict))]))
                    out += f'-- {stab_d_conf} of these have conflicting interpretations of pathogenicity: '

                    # List mutations with conflicting interpretations
                    conf_list = data.index[(data[col].isin(d_s) & \
                                        data[clinvar_clmn].astype(str).str.contains(conflict))].to_list()

                    # Print list of mutations only if the are less than 10
                    out += out_list(conf_list)

                except KeyError:
                    log.warning('Missing ClinVar Interpretation column, please check the input file. Continuing...')
                    out += f'ClinVar Interpretation column not available yet; not possible to predict the number ' \
                           f'of vus and/or variants with conflicting interpretations with a destabilizing effect.\n'
            except:
                continue
    else:
        all_stability = []
        log.warning(f'- Damaging mutation for (RaSP, FoldX) and (Rosetta, FoldX) not found in MAVISp csv.')

    all_stability = [item for sublist in all_stability if sublist != [] for item in (sublist if isinstance(sublist, list) else [sublist])]

    ### Damaging local interactions ###
    all_local = []

    # If local interactions are found
    if len(filter_col_local) > 0:
        for col in filter_col_local:
            ensemble = re.findall(pattern, col)
            inter = re.findall(interactors, col)

            # Calculate number of damaging mutations
            loc_d = str(len(data[data[col].isin(d_s)]))

            # List damaging mutations
            loc_d_list = data.index[data[col].isin(d_s)].to_list()
            all_local.extend(loc_d_list)

            # Print list of mutations only if the are less than 10
            mutlist = out_list(loc_d_list)
            if ensemble:
                if 'DNA' in col:
                    out += f'- {loc_d} variants are destabilizing for the LOCAL_INTERACTION with DNA {inter}. {ensemble[0]}. {mutlist}'
                    continue
                if inter:
                    out += f'- {loc_d} variants are destabilizing for the LOCAL_INTERACTION {inter} {ensemble[0]}. {mutlist}'
            else:
                if 'DNA' in col:
                    out += f'- {loc_d} variants are destabilizing for the LOCAL_INTERACTION with DNA {inter}. {mutlist}'
                    continue
                if inter:
                    out += f'- {loc_d} variants are destabilizing for the LOCAL_INTERACTION {inter}. {mutlist}'
    else:
        all_local = []
        out += '- No classification for LOCAL_INTERACTION available yet.\n'
    all_local = [item for sublist in all_local if sublist != [] for item in (sublist if isinstance(sublist, list) else [sublist])]

    ### Damaging PTM_stability ###
    if ptm_stab_clmn:
        for col in ptm_stab_clmn:
            try:
                ensemble = re.findall(pattern, col)
                # Calculate number of damaging mutations
                ptm_s_d = str(len(data[data[col] == 'damaging']))

                # List damaging mutations
                ptm_s_d_list = data.index[data[col] == 'damaging'].to_list()

                # Print list of mutations only if the are less than 10
                ptm_mutlist_stab = out_list(ptm_s_d_list)

                if ensemble:
                    out += f'- {ptm_s_d} variants are predicted with damaging effects on stability '\
                        f'after the removal of a PTM. {ensemble[0]} {ptm_mutlist_stab} '
                else:
                    out += f'- {ptm_s_d} variants are predicted with damaging effects on stability ' \
                        f'after the removal of a PTM. {ptm_mutlist_stab} '
            except:
                continue
    else:
        log.warning(f'- PTM effect in stability not found in MAVISp csv.')
        ptm_s_d_list = []


    ### Damaging PTM_regulation ###
    if ptm_reg_clmn:
        for col in ptm_reg_clmn:
            try:
                ensemble = re.findall(pattern, col)
                # Calculate number of damaging mutations
                ptm_r_d = str(len(data[data[col] == 'damaging']))

                # List damaging mutations
                ptm_r_d_list = data.index[data[col] == 'damaging'].to_list()

                # Print list of mutations only if the are less than 10
                ptm_mutlist_reg = out_list(ptm_r_d_list)
                if ensemble:
                    out += f'- {ptm_r_d} variants are likely to alter the regulation provided by ' \
                           f'a phosphorylation. {ensemble[0]} {ptm_mutlist_reg} '
                else:
                    out += f'- {ptm_r_d} variants are likely to alter the regulation provided by ' \
                           f'a phosphorylation. {ptm_mutlist_reg} '
            except:
                continue

    else:
        ptm_r_d_list = []
        log.warning(f'- PTM effect in regulation not found in MAVISp csv.')


    ### Damaging PTM_function ###
    dmg_str = ['damaging', 'potentially_damaging']
    if ptm_funct_clmn:
        for col in ptm_funct_clmn:
            try:
                ensemble = re.findall(pattern, col)

                # Calculate number of damaging mutations
                ptm_f_d = str(len(data[data[col].isin(dmg_str)]))

                # List damaging mutations
                ptm_f_d_list = data.index[data[col].isin(dmg_str)].to_list()

                # Print list of mutations only if the are less than 10
                ptm_mutlist_funct = out_list(ptm_f_d_list)
                if ensemble:
                    out += f'- {ptm_f_d} variants are potentially damaging in terms of removing ' \
                        f'the function related to a phosphorylation. {ensemble[0]} {ptm_mutlist_funct} '
                else:
                    out += f'- {ptm_f_d} variants are potentially damaging in terms of removing ' \
                        f'the function related to a phosphorylation. {ptm_mutlist_funct}'
            except:
                continue
    else:
        log.warning(f'- PTM effect in function not found in MAVISp csv.')
        ptm_f_d_list = []

    ### Assess presence of PTM module ####
    # Start counter
    count_ptm_na = 0
    d = {}
    # Iterate over the ptm columns
    ptm_module = ptm_stab_clmn + ptm_reg_clmn + ptm_funct_clmn

    # Create dict for counting ptm per ensemble
    for ptm in ptm_module:
        matches = re.findall(pattern, ptm)
        if matches:
            d[matches[0]] = 0

    for col in ptm_module:
        ensemble = re.findall(pattern, col)
        # If column not in dataframe
        if ensemble:
            for ens in d:
                if col not in data.columns and ensemble[0] == ens:
                    d[ens] += 1
                else:
                    continue
        else:
            if col not in data.columns:
                # Increase counter by 1
                count_ptm_na += 1
    # Per ensemble if none is in the table
    for ensemble in d:
        if d[ensemble] == 3:
            out += f'- No classification for PTM available yet for the {ensemble[0]} ensemble.'
    # If none of the PTM columns is found in the table
    if count_ptm_na == 3:
        out += '- No classification for PTM available yet.\n'

    all_ptm = ptm_s_d_list + ptm_r_d_list + ptm_f_d_list

    ### Damaging long-range Allosigma2 ###
    multiple = []
    for col in allosigma_clmn:
        try:
            lr = re.findall(lr_type, col)
            # Calculate number of damaging mutations (by class)
            long_d = str(len(data[data[col] == 'destabilizing']))
            long_me = str(len(data[data[col] == 'mixed_effects']))
            long_s = str(len(data[data[col] == 'stabilizing']))
            tot = str(int(long_d) + int(long_me) + int(long_s))

            # List damaging mutations (by class)
            long_d_list = data.index[data[col] == 'destabilizing'].to_list()
            long_me_list = data.index[data[col] == 'mixed_effects'].to_list()
            long_s_list = data.index[data[col] == 'stabilizing'].to_list()
            all_long = long_d_list + long_me_list + long_s_list

            # add long range type to a list to evaluate multiple effect later
            multiple.append(all_long)

            # Print list of mutations only if the are less than 10
            out += f'- {tot} variants are damaging for LONG_RANGE {lr} '
            out += out_list(all_long)
            out += f'-- {long_d} as destabilizing '
            out += out_list(long_d_list)
            out += f'-- {long_s} as stabilizing '
            out += out_list(long_s_list)
            out += f'-- {long_me} with mixed effects '
            out += out_list(long_me_list)

        except KeyError:
            out += f'- {col} not found in MAVISp csv.\n'
            multiple.append([])
    multiple = [item for sublist in multiple if sublist != [] for item in (sublist if isinstance(sublist, list) else [sublist])]

    ### Damaging for functional sites ###
    multiple_1 = []
    for col in cofactor_active:
        try:
            funct = re.findall(funct_sites, col)
            # Calculate number of damaging mutations (by class)
            funct_sites_d = str(len(data[data[col] == 'damaging']))

            # List damaging mutations (by class)
            funct_sites_d_list = data.index[data[col] == 'damaging'].to_list()

            # add long range type to a list to evaluate multiple effect later
            multiple_1.append(funct_sites_d_list)

            # Print list of mutations only if the are less than 10
            out += f'- {funct_sites_d} variants are damaging for Functional site {funct} module '
            out += out_list(funct_sites_d_list)

        except KeyError:
            out += f'- {col} not found in MAVISp csv.\n'
            multiple_1.append([])

    multiple_1= [item for sublist in multiple_1 if sublist != [] for item in (sublist if isinstance(sublist, list) else [sublist])]

    ### Damaging for Long-Range AlloSigma2-PSN ###
    psn_d_list=[]
    if long_range_psn_clmn:
        for col in long_range_psn_clmn:
            try:
                ensemble = re.findall(pattern, col)
                    # Calculate number of damaging mutations
                psn_d = str(len(data[data[col] == 'damaging']))

                # List damaging mutations
                psn_d_list = data.index[data[col] == 'damaging'].to_list()
                # Print list of mutations only if the are less than 10
                psn_d_mutlist = out_list(psn_d_list)

                if ensemble:
                    out += f'- {psn_d} variants are predicted with damaging long range effects on pocket site(s). {ensemble[0]} {psn_d_mutlist} '
                else:
                    out += f'- {psn_d} variants are predicted with damaging long range effects on pocket site(s). {psn_d_mutlist} '
            except:
                continue
    else:
        out += f'- No classification for Long Range-PSN available yet.\n'
        psn_d_list = []

    # Add long-range and long-range psn to single list
    all_long_range = multiple + psn_d_list

    ### Multiple effects ###

    lists = [("stability", all_stability),
             ("local_interactions", all_local),
             ("long_range", all_long_range),
             ("functional site",multiple_1),
             ("ptm", all_ptm)]

    # Find shared elements for each combination size
    shared_elements = {}
    for r in range(2, len(lists) + 1):

        # Make combinations
        combinations_lists = list(combinations(lists, r))

        # Iterate over combinations
        for combination in combinations_lists:

            # Define combined name
            combined_name = ", ".join(name for name, _ in combination)

            # Initialize a set with the mutations from the first sublist
            shared = set(combination[0][1])

            # Iterate over the remaining sublists in the combination
            for _, sublist in combination[1:]:

                # Keep only the mutations that are common to all sublists
                shared.intersection_update(sublist)

            # Create an empty list for the combined name if it doesn't exist
            if combined_name not in shared_elements.keys():
                shared_elements[combined_name] = []

            # Add value to dictionary key
            shared_elements[combined_name].append(shared)

    # Remove duplicate mutations from subsets
    filt_shared_elements = remove_duplicate_mutations(shared_elements)

    # Write out string
    tmp_shared = ''

    # Iterate over filtered dictionary items
    for combination, shared in filt_shared_elements.items():

        # If the shared elements contains something
        if len(shared) > 0:

            # Iterate over each value
            for val in shared:

                # If the value is not an empty set
                if len(val) != 0:

                    # Convert val set to string
                    V = ', '.join(val)

                    # Define out string
                    tmp_shared += f"-- {V} -> [{combination}]\n"
        else:
            continue

    # Write output file
    if len(tmp_shared) > 0:
        out += "- Variants with multiple effects:\n"
        out += tmp_shared
    else:
        out += "- No variants with multiple effects have been found.\n"

    ### DAMAGING ###
    damaging_all = set(all_stability + all_local + ptm_s_d_list + ptm_r_d_list + all_long_range)
    data_d = data.query('index in @damaging_all')

    # For alphamissense output include ptm_f as well in set and df
    damaging_full = set(all_stability + all_local + all_ptm + multiple_1 + all_long_range)
    full_data_d = data.query('index in @damaging_full')

    # REVEL score > 0.5 (default)
    revel_d = data_d.index[data_d['REVEL score'] >= r_cutoff].to_list()

    out += f'- We aggregated all the variants that have at least one of the MAVISp modules ' \
        f'with a predicted damaging effect (except for PTM.function) and retained only ' \
        f'the ones with a REVEL score higher than 0.5 for a total of {len(revel_d)} variants ' \
        f'which could be of interest for further investigation:\n'
    out += f'-- {revel_d}\n'

    # Demask
    demask_d = data_d.index[(data_d['DeMaSk delta fitness'] >= d_cutoff) | (data_d['DeMaSk delta fitness'] <= -d_cutoff)].to_list()
    out += f'- We aggregated all the variants that have at least one of the MAVISp modules ' \
        f'with a predicted damaging effect (except for PTM.function) and retained only ' \
        f'the ones with a Demask delta fitness >= {d_cutoff} or <= -{d_cutoff} ' \
        f'for a total of {len(demask_d)} variants which could be of interest for further investigation:\n'
    out += f'-- {demask_d}\n'

    try:
        # EVE
        eve_d = data_d.index[data_d['EVE classification (25% Uncertain)'] == 'pathogenic'].to_list()
        out += f'- We aggregated all the variants that have at least one of the MAVISp modules ' \
            f'with a predicted damaging effect (except for PTM.function) and retained only ' \
            f'the ones with an EVE classification of Pathogenic for a total of {len(eve_d)} variants ' \
            f'which could be of interest for further investigation:\n '
        out += f'-- {eve_d}\n'
    except KeyError:
        out += '\n- EVE grouped classification not available.\n'

    # Initialise AlphaMissense df
    alpha_d_df = pd.DataFrame()

    try:
        # Alphamissense
        # Subset the df for damaging mutations and with alphamissense pathogenic classification
        alpha_d_df = full_data_d[full_data_d['AlphaMissense classification'] == 'pathogenic']
        alpha_d = data_d.index[data_d['AlphaMissense classification'] == 'pathogenic'].to_list()
        out += f'- We aggregated all the variants that have at least one of the MAVISp modules ' \
               f'with a predicted damaging effect (except for PTM.function) and retained only ' \
               f'the ones with an AlphaMissense classification of Pathogenic for a total of {len(alpha_d)} variants ' \
               f'which could be of interest for further investigation:\n '
        out += f'-- {alpha_d}\n'
    except KeyError:
        log.warning(f'AlphaMissense grouped classification not available. Outfile will be empty...')
        out += '\n- AlphaMissense grouped classification not available.\n'

    # Experimental Classification
    try:
        # ID experimental classification columns in df
        exp_clm = [col for col in data if col.startswith('Experimental data classification')]

        # If no experimental columns are present
        if not exp_clm:
            raise KeyError

        for col in exp_clm:
            exper_d = data_d.index[data_d[col] == 'damaging'].to_list()
            experiment = re.findall(r'\([^)]*\)', col)
            out += f'- We aggregated all the variants that have at least one of the MAVISp modules ' \
                f'with a predicted damaging effect (except for PTM.function) and retained only ' \
                f'the ones with an experimentally predicted damaging effect within the context of the experiment ' \
                f'{experiment[0]} for a total of {len(exper_d)} variants ' \
                f'which could be of interest for further investigation:\n '
            out += f'-- {exper_d}\n'
    except KeyError:
        out += '\n- Experimental classification not available.\n'

    # Early Folding Regions
    try:
        earlyfolding_reg = data.index[data['EFoldMine - part of early folding region'] == 'true'].to_list()
        out += f'\n- We identified continuous stretches of residues that are above a '\
               f'pre-determined threshold of a minimum length of 3 residues to ' \
               f'assign whether each residue is part of a early folding region ' \
               f'totalling to {len(earlyfolding_reg)} variants. {out_list(earlyfolding_reg)}\n'
    except KeyError:
        out += '\n- EfoldMine annotation not available.\n'

    # Clinvar classification
    try:
        clinvar_nan = data[clinvar_clmn].isna().sum()
        out += '\nNot reported in clinvar: ' + str(clinvar_nan) + '\n'
        clinvar_grouped = data.groupby([clinvar_clmn])[clinvar_clmn].count()
        out += str(clinvar_grouped) + '\n'
    except KeyError:
        out += '\nClinvar grouped classification not available.\n'

    return out, alpha_d_df

def remove_duplicate_mutations(dictionary):
    ''' Remove duplucate mutations in mutations with
    multiple effects.

    If a mutation is shared by more than 2 classes,
    it must be removed from the list of the subclasses.

    Parameters
    ----------
    dictionary: dict
        dictionary to be filtered
    Returns
    ----------
    filtered_dictionary: dict
        filtered dictionary

    '''

    # Make empty out dict
    filtered_dictionary = {}

    # Iterate over dict keys and val
    for key, value in dictionary.items():

        # Set subset flag to false
        is_subset = False

        # Iterate over the dict keys
        for other_key in dictionary.keys():

            # If the new key is not the same as the previous one
            # and the key is a subset of the previous one
            if key != other_key and set(key).issubset(set(other_key)):

                # Change flax to true
                is_subset = True
                break

         # If flag is set to true
        if is_subset:

            # Remove the repeated mutations
            filtered_value = [mutation for mutation in value if mutation not in dictionary[other_key]]

        else:

            # Keep the original list of mutations
            filtered_value = value

        # Replace the value associated to a key with the filtered list
        filtered_dictionary[key] = filtered_value

    return filtered_dictionary


def out_list(x):
    ''' Print list or message in summary.

    The function returns the output that will be printed
    in the summary file. If the length of the list containing
    the mutations of a category is 10 or less it will return
    the list, otherwise "see csv file for details".

    Parameters
    ----------
    x: list
        input list of interest
    Returns
    ----------
    o: str
        string containing the summary message

    '''

    # If length of list is shorter than 10
    if len(x) <= 10:

        # Print list
        o = f'{x}\n'

    else:

        # Else print message
        o = 'See csv file for details.\n'

    return o

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

    # Replace values with 0/1s
    temp_df = df.replace(mechanism_code, regex=True)

    # Define broad effect categories + corresponding regex patterns
    effect_categories = {
        'Stability': 'Stability classification',
        'Local Int.': 'Local Int.',
        'PTM': 'PTM effect',
        'Long Range': 'AlloSigma2',
        'Functional': 'Functional sites'}

    # Initialise series for storing results
    effects_summary = pd.Series(index=df.index, dtype='object')

    # Iterate over each mutation in df
    for idx, row in temp_df.iterrows():
        effects = []
        # For each effect cattegory check if effect found
        # i.e. value is 1
        for category, pattern in effect_categories.items():
            if row.filter(regex=pattern).max() == 1:
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

def map_clinvar_categories(dataframe, clinvar_dict):
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
        lambda row: clinvar_dict.get(row['ClinVar Interpretation'])
        if pd.notna(row['Mutation sources']) and 'clinvar' in row['Mutation sources'].lower()
        else None,
        axis=1)
    return dataframe


def filter_am_summary(am_summary, df, amx_flag):
    """Filter the AlphafoldMissense summary dataframe

    This function filters the input dataframe to keep only columns of interest
    with a classification, and that are present in the original input.

    Parameters
    ----------
    am_summary : dataframe
        dataframe for alphamissense summary
    df : dataframe
        full data frame
    amx_flag: bool
        denotes whether to filter df on G/D thresholds
    Returns
    -------
    dataframe : dataframe
        filtered dataframe for alphamissense summary
    """

    f = lambda x: 'Stability classification' in x or \
              'Local Int. classification' in x or \
              'Local Int. With DNA classification' in x or \
              'AlloSigma2 predicted consequence' in x or \
              'AlloSigma2-PSN classification' in x or \
              'PTM effect' in x or \
              'Functional sites' in x or \
              ('Mutation' in x and 'Mutation sources' not in x)

    # Filter the columns to keep only classification info
    filter_columns = [col for col in am_summary.columns if f(col)]

    # Create a new df with the selected columns
    filtered_am = am_summary[filter_columns]

    # Filter output according to previously applied args
    filtered_am = filtered_am[filtered_am.index.isin(df.index)]

    if amx_flag:
        # Filter for LOF/GOF in at least one Gemme/Demask 
        filtered_index = df[
        (df['DeMaSk predicted consequence'].isin([7, 8]) if 'DeMaSk predicted consequence' in df.columns else False) |
        (df['GEMME predicted consequence'].isin([7, 8]) if 'GEMME predicted consequence' in df.columns else False)].index
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

    D_default = 0.25
    D_helpstr = f"Threshold to classify a mutation according to the " \
                f"DeMask score. (Default = {D_default})"
    parser.add_argument("-D","--demask_threshold",
                        default = D_default,
	                    type = float,
                        help = D_helpstr)

    G_default = 3.0
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

    pltD_helpstr = f"Plotting of Demask LOF/GOF if" \
                    f" mutation is above demask threshold. " \
                    f"(Default = None)"
    parser.add_argument("-pltD", "--plot_Demask",
                        action = 'store_true',
                        help = pltD_helpstr)

    pltC_helpstr =  f"Plotting of Clinvar variants " \
                    f"choose from: all, uncertain, " \
                    f"benign, likely benign, pathogenic, " \
                    f"likely pathogenic, conflicting. " \
                    f"For combinations, provide space separated " \
                    f"options (e.g. benign uncertain)" \
                    f"(Default = None)"
    parser.add_argument("-pltC", "--plot_Clinvar",
                        nargs='+',
                        choices = ["all",
                                    "uncertain",
                                    "benign",
                                    "likely_benign",
                                    "pathogenic",
                                    "likely_pathogenic",
                                    "conflicting"],
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

    AMx_helpstr =   f"Include AM pathogenic mutations that " \
                    f"also satisfy at least one of either GEMME " \
                    f"or DeMaSk thresholds for " \
                    f"LOF/GOF, default or specificied using respective " \
                    f"flags. "
    parser.add_argument("-amx", "--AMx",
                        action = 'store_true',
                        help = AMx_helpstr)
    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    ################################# INPUT #################################

    # Read input

    full_df = pd.read_csv(args.input, index_col = 'Mutation')

    df , dataframe =  process_input(full_df = full_df,
                                    r_cutoff = args.revel_threshold,
                                    d_cutoff = args.demask_threshold,
                                    g_cutoff= args.gemme_threshold)

    # Only save individual pngs of plots if flags r/m used
    save_png = False
    # Define an empty df if flags -pltC/-colC are not used
    clinvar_mapped_df = None

    # Assert the user didn't provide both -r and -m options
    if args.residues and args.mutations:
        log.error("Please choose only one option between user-defined " \
                  "positions (-r) and mutations (-m). Exiting...")
        sys.exit(1)

    # Keep only mutations associated to user-defined positions
    if args.residues:

        # Assert length of the argument; must be 1
        # as expecting comma-separated values
        if len(args.residues) > 1:
            log.error('Please provide comma-separated residues. Exiting...')
            sys.exit(1)

        # Extract and store in a temporary 'pos' column the
        # positions associated to the mutations
        df['pos'] = df.index.astype(str).str[1:-1]
        log.info(f'Using residues: {args.residues[0]}')

        # Raise warning if position not in dataframe
        for r in args.residues[0]:
            if r not in df['pos'].to_list():
                log.warning(f'{r} not found in table. Continuing...')

        # Filter dataframe according to the user-defined residues
        df = df[df['pos'].isin(args.residues[0])]

        # Remove temporary 'pos' column
        df.drop(columns= ['pos'],
                inplace=True)

        # Save individual pngs
        save_png = True

    # Keep only user-defined mutations
    if args.mutations:

        # Assert length of the argument; must be 1
        # as expecting comma-separated values
        if len(args.mutations) > 1:
            log.error('Please provide comma-separated mutations. Exiting...')
            sys.exit(1)
        log.info(f'Using mutations: {args.mutations[0]}')

        # Raise warning if mutation not in dataframe
        for m in args.mutations[0]:
            if m not in df.index.to_list():
                log.warning(f'{m} not found in table. Continuing...')

        # Filter dataframe according to the user-defined mutations
        df = df[df.index.isin(args.mutations[0])]

        # Save individual pngs
        save_png = True

    # If the flag is not selected drop REVEL
    if not args.plot_Revel:
        df.drop(columns=['REVEL'],
                inplace=True)

    # If the flag is not selected drop Demask consequence
    if not args.plot_Demask:
        df.drop(columns=['DeMaSk predicted consequence'],
                    inplace=True)

    # If the flag is selected drop DeMaSk fitness for plotting
    if args.plot_Demask:
        df = df.drop(columns=['DeMaSk'])
    
    # Initialize indexes for additive filtering
    source_idx = pd.Index([])
    clinvar_idx = pd.Index([])
    
    if args.plot_Source:
        try:
            # Accommodate multiple user options
            pattern = '|'.join(args.plot_Source)
            
            # Filter df based on desired mutation source
            source_idx = dataframe[dataframe['Mutation sources'].str.contains(pattern, case=False, na=False)].index
            
            # Raise error if filtering leaves df empty
            if source_idx.empty:
                log.warning(f"No mutations found matching the mutation source: {args.plot_Source}. Exiting...")
                sys.exit(1)
            else:
                log.info(f"Found {len(source_idx)} mutations matching the source filter: {args.plot_Source}")

        except Exception as e:
            log.error(f"Error occurred: {e}. Exiting...")
            sys.exit(1)

    # Assert the user didn't provide both -pltC and -colC options
    if args.plot_Clinvar and args.color_Clinvar:
        log.error("Please choose only one option between " \
                  "plotting specific ClinVar variants (-pltC) and " \
                  "coloring ClinVar variants (-colC) among all mutations. " \
                  "Exiting...")
        sys.exit(1)

    # Assert the required column is present
    if args.plot_Clinvar or args.color_Clinvar:
        
        if 'ClinVar Interpretation' not in dataframe.columns:
            log.error("ClinVar Interpretation column is missing from the input data.")
            sys.exit(1)
        try:
            # Load the ClinVar dictionary
            clinvar_dict = load_clinvar_dict('dictionary.csv')
            # Filter df according to poss. previous arguments
            clinvar_mapped_df = dataframe.loc[dataframe.index.isin(df.index)]
            # Map dictionary to df
            clinvar_mapped_df = map_clinvar_categories(clinvar_mapped_df, clinvar_dict)

        except FileNotFoundError as e:
            log.error(f"ClinVar dictionary file not found: {e}. Exiting...")
            sys.exit(1)
        except Exception as e:
            log.error(f"Error occurred: {e}. Exiting...")
            sys.exit(1)

    # Keep only user-defined mutations according to ClinVar annotation
    if args.plot_Clinvar:
        
        try:
            # Define a mapping of user-input to categories
            filter_terms = {
                'benign': ['benign'],
                'likely_benign' : ['likely benign'],
                'pathogenic': ['pathogenic'],
                'likely_pathogenic': ['likely pathogenic'],
                'uncertain': ['uncertain'],
                'conflicting': ['conflicting']}

            # Filter dataframe according to flag option
            if 'all' in args.plot_Clinvar:
                clinvar_idx = clinvar_mapped_df[clinvar_mapped_df['ClinVar Category'].notna()].index
            else:
                # Accommodate multiple user options
                for clinvar_option in args.plot_Clinvar:
                    if clinvar_option in filter_terms:
                        category_values = filter_terms[clinvar_option]
                        # Match categories
                        pattern = '|'.join(category_values)
                        #  Get indexes of mutations matching current ClinVar category
                        indexes_to_add = clinvar_mapped_df[clinvar_mapped_df['ClinVar Category'].str.contains(pattern, case=False, na=False)].index
                        clinvar_idx = clinvar_idx.union(indexes_to_add)
                    else:
                        log.warning(f"Unrecognized ClinVar filter option: {clinvar_option}")
                        sys.exit(1)
            
            # Raise error if filtering leaves df empty
            if len(clinvar_idx) == 0:
                log.warning(f"No mutations found matching the ClinVar filter {args.plot_Clinvar}. Exiting...")
                sys.exit(1)
            else:
                # Log the number of mutations remaining after filtering
                log.info(f"Found {len(clinvar_idx)} mutations matching the ClinVar filter {args.plot_Clinvar}")

        except Exception as e:
            log.error(f"Error occurred: {e}. Exiting...")
            sys.exit(1)
    
    # Additively filter if -pltS -pltC
    if args.plot_Source or args.plot_Clinvar:
        #Save individual pngs
        save_png = True
        combined_indexes = source_idx.union(clinvar_idx)
        df = df[df.index.isin(combined_indexes)]
        
        log.info(f"Found {len(df)} mutations matching filters.")
        # Raise error if filtering leaves df empty
        if df.empty:
                log.warning("No mutations remain after filtering. Exiting...")
                sys.exit(1)

    ############################### SUMMARY ###############################

    # Write summary output file
    with open('log.txt', 'w') as out:
        summary, am_summary = generate_summary(data = dataframe,
                                d_cutoff = args.demask_threshold,
                                r_cutoff= args.revel_threshold)
        out.write(summary)

    ################################# PLOT #################################
    '''
    filter_col = [col for col in df if col.startswith("Local Int.")]
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
    figures = plot(df = df,
                   full_df = clinvar_mapped_df,
                   width = args.figsize[0],
                   height = args.figsize[1],
                   xlim = args.x_lim,
                   clinvar_flag = bool(args.plot_Clinvar),
                   clinvar_col = bool(args.color_Clinvar))

    with PdfPages(f"{args.output}.pdf") as pdf:
        for i, figure in enumerate(figures):
            figure.savefig(pdf, format='pdf', dpi=300)

            if save_png:
                figure.savefig(f'{args.output}_{i}.png', dpi=300)


################################# AM CSV #################################

    filtered_am = filter_am_summary(am_summary, df, args.AMx)

    filtered_am.to_csv('alphamissense_out.csv', index=True)

if __name__ == '__main__':
    main()
