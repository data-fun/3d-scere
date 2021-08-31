import math
import matplotlib
import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import plotly
import plotly.graph_objects as go
import plotly.express as px
import sqlite3


def display_module_version():
    """Display dependencies versions.
    """
    print("sqlite3 version:", sqlite3.version)
    print("pandas version:", pd.__version__)
    print("matplotlib version:", matplotlib.__version__)
    print("numpy version:", np.__version__)
    print("plotly version:", plotly.__version__)
    print("ipywidgets version:", widgets.__version__)




def format_coordinates(coordinates, space_between_chromosomes):
    """Format the locus coordinates for Plotly visualization. 
    
    Each locus is represented by three rows: 
    x1, x2 (the two values are in the column x) and none.
    The third row allow the separation between lines.

    Parameters
    ----------
    coordinates : Pandas dataframe
        Dataframe created from SQL query (contains locus coordinates).
    space_between_chromosomes : int
        Graphical space to leave between chromosomes.

    Returns
    -------
    Pandas dataframe
    """
    
    genome_data = pd.DataFrame(columns=coordinates.columns)
    row_null = {"Start_coordinate": "none", "Stop_coordinate": "none"}
    
    for chromosome_id in range(1, coordinates["Chromosome"].max() + 1):
        chrom = coordinates[coordinates["Chromosome"] == chromosome_id]
        row_one = chrom.copy()
        row_one.index = range(0, len(chrom)*3, 3)
        row_one = row_one.drop("Stop_coordinate", axis = 1)
        row_one = row_one.transpose()
        
        row_two = chrom.copy()
        row_two.index = range(1, len(chrom)*3, 3)
        row_two["Start_coordinate"] = row_two["Stop_coordinate"]
        row_two = row_two.drop("Stop_coordinate", axis = 1)
        row_two = row_two.transpose()
        
        row_three = chrom.assign(Start_coordinate = "none")
        row_three.index = range(2, len(chrom)*3, 3)
        row_three = row_three.drop("Stop_coordinate", axis = 1)
        row_three = row_three.transpose()
        
        chrom_data = pd.merge(row_one, row_two, left_index = True, right_index = True)
        chrom_data = pd.merge(chrom_data, row_three, left_index = True, right_index = True)
        chrom_data = chrom_data.transpose()
        chrom_data = chrom_data.sort_index()

        # Add y-coordinates   
        chrom_data["Stop_coordinate"] = (chromosome_id - 1) * space_between_chromosomes
        chrom_data["Stop_coordinate"] = chrom_data.apply(lambda x: x["Stop_coordinate"] + 0.2 if x["Strand"] == "C" else x["Stop_coordinate"] - 0.2, axis=1)
        
        chrom_data = chrom_data.append(row_null, ignore_index = True)
        genome_data = genome_data.append(chrom_data)
        
    genome_data = genome_data.rename(columns={"Start_coordinate": "x", "Stop_coordinate": "y"})
        
    return genome_data



# Chromosome shapes.

def get_chromosome_lenght(chrom_number):
    #SQL request
    db_connexion = sqlite3.connect('../SCERE.db')

    cursor = db_connexion.cursor()

    chromosome_length = cursor.execute("""
    SELECT length
    FROM chromosome_length
    """)
    
    chromosome_length = chromosome_length.fetchall()
    chromosome_length = pd.DataFrame(chromosome_length, columns = ["length"], index = list(range(1,18)))
    
    return chromosome_length.loc[chrom_number][0]

def format_chromosomes(y1, y2):
    
    chromosomes = pd.DataFrame(columns = ["x", "y", "Chromosome"])
    
    for c in range(1,18):
        chrom_lenght = get_chromosome_lenght(c)
        chromosomes = chromosomes.append({"x": 0, 
                                          "y": y1[c-1], 
                                          "Chromosome": 0, 
                                          "Feature_type": "0"}, ignore_index = True)
        chromosomes = chromosomes.append({"x": chrom_lenght, 
                                          "y": y1[c-1], 
                                          "Chromosome": 0, 
                                          "Feature_type": "0"}, ignore_index = True)
        chromosomes = chromosomes.append({"x": "none", 
                                          "y": "none", 
                                          "Chromosome": 0, 
                                          "Feature_type": "0"}, ignore_index = True)
        
        chromosomes = chromosomes.append({"x": 0, 
                                          "y": y2[c-1], 
                                          "Chromosome": 0, 
                                          "Feature_type": "0"}, ignore_index = True)
        chromosomes = chromosomes.append({"x": chrom_lenght, 
                                          "y": y2[c-1],
                                          "Chromosome": 0, 
                                          "Feature_type": "0"}, ignore_index = True)
        chromosomes = chromosomes.append({"x": "none", 
                                          "y": "none", 
                                          "Chromosome": 0, 
                                          "Feature_type": "0"}, ignore_index = True)
    
    return chromosomes

# Genome drawing.

def genome_drawing(genome_data, mode, parameter, 
                   values = "null", values_colors = "null", threshold = 10**40, hover = []):
    
    chromosomes = format_chromosomes(list(i + 0.2 for i in range(0,108,6)), list(i - 0.2 for i in range(0,108,6)))
    
    genome_data = chromosomes.append(genome_data)
    genome_data.index = range(1, len(genome_data) + 1)
    
    
    if mode == "continuous":
        colors = get_color_continuous(genome_data[parameter])
        colors.index = range(1, len(colors) + 1)
        
        genome_data["colors"] = colors
        
        fig = px.line(genome_data,
                      x = "x",
                      y = "y",
                      color = "colors",
                      color_discrete_map = "identity", 
                      hover_name = "Primary_SGDID")
    
    if mode == "semi_continuous":
        
        colors_and_intervals = get_color_semi_continuous(genome_data[parameter], threshold)
        colors = colors_and_intervals[0]
        colors.index = range(1, len(colors) + 1)
        genome_data[parameter] = colors
        
        intervals = colors_and_intervals[1]
        color_discrete_map = zip(intervals, px.colors.sequential.Viridis_r)
        color_discrete_map = dict(color_discrete_map)
        color_discrete_map = { "null": "lightgrey", **color_discrete_map}
        
        hover_formating = [True] * len(hover)
        hover_data = dict(zip(hover, hover_formating))
        
        fig = px.line(genome_data,
                      x = "x",
                      y = "y",
                      color = parameter,
                      color_discrete_map = color_discrete_map, 
                      hover_name = "Primary_SGDID", 
                      hover_data = {**hover_data, "y": False})
                    #no order in legend because locus are not drawed when there is an order
    
    if mode == "discreet":
        
        genome_data = get_color_discreet(genome_data, parameter, values)
        
        color_discrete_map = dict(zip(values, values_colors))
        
        fig = px.line(genome_data,
                      x = "x",
                      y = "y",
                      color = "colors", 
                      color_discrete_map = {"Other": "darkgrey", "Background": "lightgrey", **color_discrete_map},
                      hover_name = "Primary_SGDID")

    fig.update_traces(line = dict(width = 9))
    
    fig.update_layout(plot_bgcolor = "white", 
                      xaxis_showgrid = False, 
                      yaxis_showgrid = False, 
                      showlegend = True)
    
    fig.update_yaxes(tickmode = "array",
                     tickvals = list(range(0,102,6)),
                     ticktext = ["1", "2", "3", "4", "5", "6", "7", "8", "9",
                                 "10", "11", "12", "13", "14", "15", "16", "mitochondrial"],
                     title = "Chromosomes number")
    fig.update_xaxes(title = "Coordinates (bp)")
    
    fig.update_layout(hoverlabel = dict(bgcolor="white",
                                        font_size=16))
    
    return fig

# Adding color.

def get_color_discreet(genome_data, parameter, values):
    
    genome_data.loc[genome_data[parameter] != values[0], "colors"] = "Other"
    
    for v in values :
        genome_data.loc[genome_data[parameter] == v, "colors"] = v
    
    genome_data.loc[genome_data["Chromosome"] == 0, "colors"] = "Background"
    
    return genome_data

def get_color_semi_continuous(parameter, threshold):
    
    parameter = parameter.apply(float)
    parameter.index = range(1, len(parameter) + 1)
    limit = min(parameter.max(), threshold)
    STEP = (limit/9)
    
    conditions = [(parameter <= STEP), (parameter <= STEP * 2), (parameter <= STEP * 3),
                  (parameter <= STEP * 4), (parameter <= STEP * 5), (parameter <= STEP * 6), 
                  (parameter <= STEP * 7), (parameter <= STEP * 8), (parameter <= STEP * 9), 
                  (parameter > STEP * 9)]
    
    choices = ["0-" + str(round(STEP)), 
               str(round(STEP)) + "-" + str(round(STEP * 2)), 
               str(round(STEP * 2)) + "-" + str(round(STEP * 3)), 
               str(round(STEP * 3)) + "-" + str(round(STEP * 4)), 
               str(round(STEP * 4)) + "-" + str(round(STEP * 5)), 
               str(round(STEP * 5)) + "-" + str(round(STEP * 6)), 
               str(round(STEP * 6)) + "-" + str(round(STEP * 7)), 
               str(round(STEP * 7)) + "-" + str(round(STEP * 8)), 
               str(round(STEP * 8)) + "-" + str(round(STEP * 9)), 
               str(round(STEP * 9)) + "<"]
    
    
    right_parameter = np.select(conditions, choices, default = "null")
    
    return [pd.Series(right_parameter), choices]

def get_color_continuous(parameter):
    cmap = matplotlib.cm.get_cmap('viridis')
    parameter = parameter.apply(float)
    MIN = min(parameter)
    MAX = max(parameter)
    colors = []

    for i in range(1, len(parameter)+1):
        
        if parameter[i] == 0 or parameter[i] == "" or parameter[i] == "NaN" :
            color = "lightgrey"
            colors = colors + [color]
        
        else :
            color = cmap((parameter[i] - MIN) / (MAX - MIN))
            color = "rgb" + str(color[:3])
            
            colors = colors + [color]
    
    return pd.Series(colors)






