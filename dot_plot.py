#!/usr/bin/env python

# Copyright (C) 2023 Ludovica Beltrame, Simone Scrima, Matteo Tiberti
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
    # Define a score as the original width - 3 (default space for y labels)
    # and divided by the number of mutations
    score = (original_w - 3) / original_m

    # Get the new width by adding 3 (default space for y lables)
    # to the produt of the score and the remaining mutations
    new_w = 3 + score * x_mut

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
    mapping_colors = {1 : 'black',
                      5 : 'black',
                      6 : 'black',
                      2 : '#d5bdaf',
                      3 : 'white',
                      4 : 'lightgray',
                      7 : '#8E5466',
                      8 : '#636D87'}

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
                 'Functional sites (cofactor)',
                 'Functional sites (active site)']
    function = ['Local Int.',
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
            tick_colors.append('#8ac926')

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
                    x == 'Local Int. With DNA' or \
                    'Functional sites (cofactor)' in x or \
                    'Functional sites (active site)' in x or \
                    'AlloSigma2 predicted consequence - active sites' in x or \
                    'AlloSigma2 predicted consequence - cofactor sites' in x or \
                    'AlloSigma2 predicted consequence - pockets and interfaces' in x or \
                    'PTM effect in ' in x or 'REVEL score' in x or \
                    'EVE classification (25% Uncertain)' in x or \
                    'DeMaSk delta fitness' in x or \
                    'DeMaSk predicted consequence' in x or \
                    'GEMME Score (rank-normalized)' in x or \
                    'AlphaMissense classification' in x or \
                    'Mutation' in x and not 'Mutation sources' in x

    df = full_df.copy()
    df = df[df.columns[list(map(f, df.columns))]]

    # If pltRevel is True then convert REVEL score column to float
    # In case multiple REVEL scores are present it returns the average
    # of the values

    for d in [df, full_df]:
        d['REVEL score'] = d['REVEL score'].apply(convert_to_float)


    # Add REVEL score interpretation column
    df['REVEL'] = np.where(df['REVEL score'].isna(), None,
                                np.where(df['REVEL score'] >= r_cutoff,
                                    'Damaging', 'Neutral'))
                                        # Add Demask score interpretation column

    try:
        # Add GEMME score interpretation column
        df['GEMME Score (rank-normalized)'] = np.where(df['GEMME Score (rank-normalized)'].isna(), None,
                                    np.where(df['GEMME Score (rank-normalized)'] >= g_cutoff,
                                        'Damaging', 'Neutral'))
    except:
        log.warning(f'- no GEMME found in MAVISp csv.')


    # Convert Demask delta fitness into absolute value
    df['DeMaSk'] = np.where(df['DeMaSk delta fitness'].isna(), None,
                np.where(df['DeMaSk delta fitness'] >= d_cutoff,'Damaging',
                np.where(df['DeMaSk delta fitness'] <= -d_cutoff,'Damaging',
                'Neutral')))

    # Only keep Demask consequence for those mutations that satisfy the Demask threshold
    df['DeMaSk predicted consequence'] = np.where(df['DeMaSk'] == 'Damaging', df['DeMaSk predicted consequence'], None)

    # Drop REVEL score column
    df.drop(columns = ['REVEL score','DeMaSk delta fitness'],
            inplace = True)

    # Sort df columns by group: stability effects, functional effects, others
    df = df[
    [col for col in df.columns if 'functional' in col.lower()] +
    [col for col in df.columns if 'stability' in col.lower()] +
    [col for col in df.columns if 'functional' not in col.lower() and 'stability' not in col.lower()]
    ]

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
    effect_code = {'(?i)destabilizing': 1,
                   '(?i)pathogenic': 1,
                   '(?i)damaging': 1,
                   'mixed_effects': 5,
                   '(?i)stabilizing': 6,
                   '(?i)neutral': 2,
                   '(?i)benign': 2,
                   None : 3,
                   '(?i)uncertain': 4,
                   'ambiguous': 4,
                   'loss_of_function' : 7,
                   'gain_of_function' : 8}
    # Replace the string nomenclature with the
    # effect code
    df.replace(effect_code,
               regex = True,
               inplace = True)

    # Ensure 'DeMaSk predicted consequence' is the last column
    if 'DeMaSk predicted consequence' in df.columns:
        columns_to_keep = [col for col in df.columns if col != 'DeMaSk predicted consequence']
        columns_to_keep.append('DeMaSk predicted consequence')
        df = df[columns_to_keep]

    return df, full_df

def plot(df, width, height, xlim, use_demask_classification):
    ''' Plot

    The function is aimed at plotting a matrix showing the effect of
    selected mutations

    Parameters
    ----------
    data: dataframe
        input MAVISp dataframe
    Returns
    ----------
    data: dataframe
        output dataframe in which the classification is
        converted to number to facilitate the plotti.ng

    '''
    # Define font of the plot
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = [ 'Arial' ]

    # List y label colors
    yticklabel_color = color_yticklabels(labels = df.columns.values)

    # Lower and upper limits
    l = 0
    u = xlim

    expected_plots = math.ceil(df.shape[0]/u)

    figures = []

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
        for ticklabel, tickcolor in zip(plt.gca().get_yticklabels(), yticklabel_color):
            ticklabel.set_color(tickcolor)

        if data.shape[0] <= 10:
            plt.margins(x = 0.4)

        # Put a legend below current axis
        legend_dict = {'black': 'Damaging',
                    '#d5bdaf': 'Neutral',
                    'lightgray': 'Uncertain',
                    'white': 'Not_available',
                    }
        if use_demask_classification is True:
            legend_dict.update({'#8E5466' : 'Loss of Function',
                    '#636D87' : 'Gain of Function'})

        legend_list = []

        for k in legend_dict.keys():
            legend_list.append(Line2D([0], [0],
                                    color = 'w',
                                    marker = 'o',
                                    markersize = 10,
                                    markerfacecolor=k,
                                    label=legend_dict[k],
                                    markeredgecolor='k'))

        ax.legend(handles = legend_list,
                loc='upper right',
                bbox_to_anchor=(1, 1.2),
                ncol = 6,
                fancybox=True,
                shadow=True)

        fig.tight_layout()

        figures.append(fig)

        # Update upper and lower limits
        l += xlim
        u += xlim
    # close merged pdf file

    return figures

def generate_summary(data):
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

    # Convert all the strings to lowercase
    data = data.apply(lambda x: x.str.lower() if x.dtype == object else x)


    # Colnames ptm and allosigma

    filter_col_stability = [col for col in data if '(Rosetta, FoldX)' in col or '(RaSP, FoldX)' in col]
    filter_col_local = [col for col in data if ('Local Int. classification') in col]
    ptm_stab_clmn = [col for col in data if 'PTM effect in stability' in col]
    ptm_reg_clmn = [col for col in data if 'PTM effect in regulation' in col]
    ptm_funct_clmn = [col for col in data if 'PTM effect in function' in col]
    allosigma_clmn = [col for col in data if 'AlloSigma2 predicted consequence' in col]
    cofactor_active = [col for col in data if 'Functional sites' in col]
    # # Check if columns with well-defined names are in the dataframe
    # for col in [ptm_stab_clmn, ptm_reg_clmn, ptm_funct_clmn, allosigma_clmn]:
    #     if col not in data.columns:
    #         log.warning(f'{col} not found in csv. Skipping associated analyses for summary file.')

    ### Damaging stability ###
    all_stability = []

    # Define columns and effects' variables
    clinvar_clmn = 'ClinVar Interpretation'
    vus = '(?i)uncertain'
    conflict = 'conflicting interpretations of pathogenicity'
    d_s = ['destabilizing', 'stabilizing']

    pattern = r'\[.*?\]'
    interactors = r'\((.*?)\)'
    lr_type = r'-(.*)'

    if filter_col_stability:
        for col in filter_col_stability:
            try:
                ensemble = re.findall(pattern, col)
                # Calculate number of damaging mutations
                stab_d = str(len(data[data[col].isin(d_s)]))

                # List damaging mutations
                all_stability_list = data.index[data[col].isin(d_s)].to_list()
                all_stability.extend(all_stability_list)
            except:
                continue

    ### Damaging local interactions ###
    all_local = []

    # If local interactions are found
    if len(filter_col_local) > 0:
        for col in filter_col_local:
            ensemble = re.findall(pattern, col)
            inter = re.findall(interactors, col) # search for interactors
            # Calculate number of damaging mutations
            loc_d = str(len(data[data[col].isin(d_s)]))

            # List damaging mutations
            loc_d_list = data.index[data[col].isin(d_s)].to_list()
            all_local.extend(loc_d_list)

    ### Damaging PTM_stability ###
    if ptm_stab_clmn:
        for col in ptm_stab_clmn:
            try:
                ensemble = re.findall(pattern, col)
                # Calculate number of damaging mutations
                ptm_s_d = str(len(data[data[col] == 'damaging']))

                # List damaging mutations
                ptm_s_d_list = data.index[data[col] == 'damaging'].to_list()
            except:
                continue
    else:
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
            except:
                continue
    else:
        log.warning(f'- PTM effect in regulation not found in MAVISp csv.')
        ptm_r_d_list = []



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
            except:
                continue
    else:
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

        except KeyError:
            multiple.append([])

    multiple = [item for sublist in multiple if sublist != [] for item in (sublist if isinstance(sublist, list) else [sublist])]


    multiple_1 = []
    for col in cofactor_active:
        try:
            funct = re.findall(funct_sites, col)
            # Calculate number of damaging mutations (by class)
            funct_sites_d = str(len(data[data[col] == 'damaging']))
            #tot = str(int(funct_sites_d) + int(long_me) + int(long_s))

            # List damaging mutations (by class)
            funct_sites_d_list = data.index[data[col] == 'damaging'].to_list()
            #all_long = long_d_list + long_me_list + long_s_list

            # add long range type to a list to evaluate multiple effect later
            multiple_1.append(funct_sites_d_list)
        except KeyError:
            multiple_1.append([])

    multiple_1= [item for sublist in multiple_1 if sublist != [] for item in (sublist if isinstance(sublist, list) else [sublist])]


    ### DAMAGING ###
    damaging_full = set(all_stability + all_local + all_ptm + multiple + multiple_1)
    full_data_d = data.query('index in @damaging_full')


    # Initialise AlphaMissense df
    alpha_d_df = pd.DataFrame()

    try:
        # Alphamissense
        # Subset the df for damaging mutations and with alphamissense pathogenic classification
        alpha_d_df = full_data_d[full_data_d['AlphaMissense classification'] == 'pathogenic']
        #alpha_d = data_d.index[data_d['AlphaMissense classification'] == 'pathogenic'].to_list()
    except KeyError:
        pass

    f = lambda x: 'Stability classification' in x or \
              'Local Int. classification' in x or \
              'AlloSigma2 predicted consequence' in x or \
              'PTM effect' in x or \
              'Functional sites' in x or \
              ('Mutation' in x and 'Mutation sources' not in x)

    # Filter the columns to keep only classification info
    filter_columns = [col for col in alpha_d_df.columns if f(col)]

    # Create a new df with the selected columns
    filtered_am = alpha_d_df[filter_columns]

    return filtered_am
