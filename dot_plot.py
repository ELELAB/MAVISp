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
                      4 : 'lightgray'}

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
    stability = ['Stability classification', 'PTM effect in stability']
    function = ['Local Int.', 'PTM', 'AlloSigma2']

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

def process_input(df, d_cutoff, r_cutoff, all):
    '''Process MAVISp aggregated table.

    The function takes as input a MAVISp csv file and returns
    a dataframe formatted in a plot-compatible way.

    Parameters
    ----------
    df: dataframe
	   input MAVISp dataframe
    r_cutoff: float
        REVEL score cutoff
    d_cutoff:
        DeMaSk score cut-off
    all:
        whether to keep both Rosetta and Rasp results
    Returns
    ----------
    df: dataframe
        output dataframe in which the classification is
        converted to number to facilitate the plotting.
    '''

    # Define columns of interest

    df['REVEL score'] = df['REVEL score'].apply(convert_to_float)
    df_full = df.copy()

    likes = ['(Rosetta, FoldX)',
              'Local Int. classification',
              'Local Int. With DNA',
              'AlloSigma2 predicted consequence - active sites',
              'AlloSigma2 predicted consequence - cofactor sites',
              'AlloSigma2 predicted consequence - pockets and interfaces',
              'Functional sites (cofactor)',
              'Functional sites (active site)',
              'PTM effect in ',
              'REVEL score',
              'EVE classification (25% Uncertain)',
              'DeMaSk delta fitness',
              'AlphaMissense classification',
              'Mutation']
    drops = ['Mutation sources']

    if all:
        likes.insert(1, '(RaSP, FoldX)')

    selected_cols = []
    for c in df.columns:
        for l in likes:
            if l in c:
                selected_cols.append(c)
                break
    df = df[selected_cols]
    df = df.drop(columns=drops)

    # Add REVEL score interpretation column
    df['REVEL'] = np.where(df['REVEL score'].isna(), None,
                                np.where(df['REVEL score'] >= r_cutoff,
                                    'Damaging', 'Neutral'))
                                        # Add Demask score interpretation column

    # Convert Demask delta fitness into absolute value
    df['DeMaSk'] = np.where(df['DeMaSk delta fitness'].isna(), None,
                   np.where(df['DeMaSk delta fitness'] >= d_cutoff,'Damaging',
                   np.where(df['DeMaSk delta fitness'] <= -d_cutoff,'Damaging',
                   'Neutral')))

    # Drop REVEL score column
    df.drop(columns = ['REVEL score','DeMaSk delta fitness'],
            inplace = True)

    # Sort df columns by group: stability effects, functional effects, others
    df = df[[col for col in df.columns if 'stability' in col.lower()] +
            [col for col in df.columns if not 'stability' in col.lower()]]

    # Define a dictionary of effect: code_number
    # Damaging : 1, 5, 6
    # Neutral : 2
    # Not available : 3
    # Uncertain : 4
    # Added support for EVE
    # -> pathogenic = 1
    # -> benign = 2
    # -> uncertain = 4
    effect_code = {'(?i)destabilizing': 1,
                   '(?i)pathogenic': 1,
                   '(?i)damaging': 1,
                   'mixed_effects': 5,
                   '(?i)stabilizing': 6,
                   '(?i)neutral': 2,
                   '(?i)benign': 2,
                   None : 3,
                   '(?i)uncertain': 4,
                   'ambiguous': 4}
    # Replace the string nomenclature with the
    # effect code
    df.replace(effect_code,
               regex = True,
               inplace = True)
    return df, df_full

def plot(df, width, height, xlim, reshape_last=True):
    '''
    The function is aimed at plotting a matrix showing the effect of
    selected mutations. Plots several dot plots, depending on the
    requested number of mutations per dot-plot

    Parameters
    ----------
    data: pandas.DataFrame
        input MAVISp dataframe
    width: float
        width of the plot in inches
    height: float
        height of the plot in inches
    xlim: int
        number of mutations per dot plot
    reshape_last: bool
        whether to change the figure width of a panel
        if the number of mutations is lower than xlim

    Returns
    ----------
    data: list of matplotlib.pyplot.Figure objects
        list of matplotlib Figure objects, one per plot
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
            if filtered_df.shape[0] != xlim and reshape_last:

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
                    'white': 'Not_available'}

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
                ncol = 4,
                fancybox=True,
                shadow=True)

        # Save plot as png and pdf
        fig.tight_layout()

        # add figure to output list
        figures.append(fig)

        # Update upper and lower limits
        l += xlim
        u += xlim

    return figures

