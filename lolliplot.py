"""
# Copyright (C) 2024 Karolina Krzesińska <kzokr@dtu.dk> 
# Danish Cancer Institute 

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
"""
# Imports
import math
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import argparse
import logging as log

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

def process_input(df):
    '''Processes AlphaMissense CSV file. 
    
    The funtion takes as input a dataframe from a MAVISp CSV file and returns a
    dataframe with converted values for subsequent plotting.

    Parameters 
    ----------
    data: dataframe
        input AlphaMissense dataframe from CSV.
    Returns
    ----------
    df: dataframe
        filtered dataframe for classification columns where values 
        have been converted to 0/1s for plotting.
    '''
    # Define columns with effect classification
    f = lambda x: 'Stability classification' in x or \
                        'Local Int. ' in x or \
                        'AlloSigma2 predicted consequence' in x or \
                        'PTM effect' in x or \
                        'Mutation' in x and not 'Mutation sources' in x
    
    # Check at least one column for each broad effect
    # has been identified
    regex_patterns = {
        'Stability classification': 'Stability',
        'AlloSigma2 predicted consequence': 'Long Range',
        'Local Int. ': 'Local Int.',
        'PTM effect': 'PTM'
    }
    # Print warning if not one column identified
    for pattern, name in regex_patterns.items():
        if not any(df.filter(regex=pattern).columns):
            log.warning(f" No columns matching '{pattern}' found for '{name}'.")

    # Define a dictionary of effect nomenclature 
    # to numerical value
    mechanism_code = {'(?i)destabilizing': 1,
                   '(?i)pathogenic': 1,
                   '(?i)damaging': 1,
                   'mixed_effects': 1,
                   '(?i)stabilizing': 1,
                   '(?i)neutral': 0,
                   '(?i)benign': 0,
                   None : 0,
                   '(?i)uncertain': 0,
                   'ambiguous': 0}
    # Replace the string values with the effect code
    df.replace(mechanism_code,
               regex = True,
               inplace = True)
    
    # Group columns by broad effect 
    df['Stability'] = df.filter(regex='Stability classification').max(axis=1)
    df['Long Range'] = df.filter(regex='AlloSigma2 predicted consequence').max(axis=1)
    df['Local Int.'] = df.filter(regex='Local Int.').max(axis=1)
    df['PTM'] = df.filter(regex='PTM effect').max(axis=1)
    
    # Filter df for relevant columns
    df = df[['Stability', 'Long Range', 'Local Int.', 'PTM']]
    
    return df

def plot(df, xlim):
    '''Plot stacked lollis per identified mutational effect. 

    The function takes the filtered df and plots a lolliplot
    per identified effect for each mutation.

    
    Parameters
    ----------
    data: dataframe
        filtered dataframe with replaced values to 0/1s.
    xlim: int
        number of mutations to plot on the x axis.

    Returns
    ----------
        plots saved as pngs
    '''
    
    # Set the font
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = [ 'Arial' ]
    mpl.rcParams['font.size'] = 10

    # Define the colors for each effect
    effect_colors = {
            'Stability': '#E34A6F',
            'Local Int.': '#3E5965',
            'PTM': '#FF8B1F',
            'Long Range': '#F6CE3C'}
    
    # Define range of df to be plotted
    lower = 0 
    upper = xlim 

    # Define default width and height of fig
    width = 15
    height = 6

    # Define expected number of plots 
    expected_plots = math.ceil(df.shape[0]/xlim)

    out_figures = []

    for p in range(expected_plots):
        # Filter df for range of mutations to be plotted
        filtered_df = df.iloc[lower:upper, :]

        # If last plot adjust width 
        if p == range(0, expected_plots)[-1]:
            
            # If the remainig mutations are less than xlim
            if filtered_df.shape[0] != xlim:
                # Adjust width of plot 
                width = adjust_figsize_w(x_mut = filtered_df.shape[0],
                                         original_w = width,
                                         original_m = xlim)                
        # Define figure and axis
        fig, ax = plt.subplots(figsize=(width, height))

        # Initialize the x positions for each mutation
        x_positions = np.arange(len(filtered_df))

        # Iterate over mutations + repsecitve effects
        for x_pos, (mutation, row) in enumerate(filtered_df.iterrows()):
            y_pos = 0
            
            for effect, color in effect_colors.items():
                value = row[effect]
                # If effect identified
                if value == 1:
                    # Increment y-axis position for each effect plotted
                    y_pos += 1
                    # Define the sticks of the lolli
                    ax.vlines(x_pos, 
                            y_pos - 0.9, 
                            y_pos, 
                            color=color, 
                            alpha=0.7)
                    # Plot the head of the lolli
                    ax.scatter(x_pos, 
                            y_pos, 
                            color=color, 
                            s=100)
        
        # Define legend elements
        legend_elements = [plt.Line2D([0], [0], 
                marker='o', 
                color=color, 
                label=effect, 
                markersize=10,
                linestyle='')
                for effect, color in effect_colors.items()]

        ax.legend(handles=legend_elements, 
                title='Effects', 
                loc='upper left', 
                bbox_to_anchor=(1, 1))

        # Customizing the plot
        ax.set_ylim(0, len(effect_colors)+0.5)
        ax.set_xticks(x_positions)      
        ax.set_xticklabels(filtered_df.index, rotation=90)       
        ax.set_xlabel('Mutations')
        ax.set_title('Identified Mutational Effects')
        ax.spines['bottom'].set_visible(True)
        ax.yaxis.set_visible(False)
                
        # Update ranges for next plot
        lower += xlim
        upper += xlim

        out_figures.append(fig)

    return out_figures

def main(data, xlim):
    # Logging configuration
    log.basicConfig(level=log.INFO, \
            format='%(levelname)s - %(message)s')

    # Read input file
    df = read_input(data)
    
    # Plot data 
    plot(df, xlim)

if __name__ == "__main__":
    description="Plotting AM pathogenic mutations and " \
        "identified respective module effects in categories:" \
        "stability, ptm, long range, local int."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', 
                        dest= "input_file", 
                        type=str, 
                        required=True, 
                        help="AlphaMissense output csv of dotplot.py")
    xlim_d = 15
    parser.add_argument('-x',
                        dest= "xlim",
                        default = xlim_d, 
                        type = int,
                        help = "Number of mutations to plot on the x axis." \
                            " (default 15)")

    args = parser.parse_args()

    main(args.input_file, args.xlim)
