import matplotlib
import numpy as np
import pandas as pd
import plotly
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
    """Get the chromosome length from the SQL database.

    Parameters
    ----------
    chrom_number : int

    Returns
    -------
    int
    """

    #SQL request
    db_connexion = sqlite3.connect('./static/SCERE.db')

    cursor = db_connexion.cursor()

    chromosome_length = cursor.execute("""
    SELECT length
    FROM chromosome_length
    """)

    chromosome_length = chromosome_length.fetchall()
    chromosome_length = pd.DataFrame(chromosome_length, columns = ["length"], index = list(range(1,18)))

    return chromosome_length.loc[chrom_number][0]

def format_chromosomes(y1, y2):
    """Format the chromosomes coordinates for Plotly visualization.

    Each chromosome is represented by three rows:
    x1, x2 (the two values are in the column x) and none.
    The third row allow the separation between lines.

    Parameters
    ----------
    y1 : list
        List containing chromosomes y coordinates (+ strand).
    y2 : list
        List containing chromosomes y coordinates (- strand).

    Returns
    -------
    Pandas dataframe
    """

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

def genome_drawing(genome_data, parameter, values = "null", values_colors = "null", hover = []):
    """Draw the 2D plotly figure, representing the 16 chromosomes (+ mitochondrial plasmid) in lightgrey and all the loci in darkgrey.

    Parameters
    ----------
    genome_data : Pandas dataframe
        2D coordinates of all the loci for Plotly visualization.
    parameter : str
        The name of the genome_data column containing the coloring parameter.
    values : list
        The values in genome_data[parameter] that will be colored according to values_colors.
    values_colors : str
        A list of color names. Loci with values[0] in the column genome_data[parameter] will be colored in values_colors[0].
    hover : list

    Returns
    -------
    Plotly figure
    """

    chromosomes = format_chromosomes(list(i + 0.2 for i in range(0,108,6)), list(i - 0.2 for i in range(0,108,6)))

    genome_data = chromosomes.append(genome_data)
    genome_data.index = range(1, len(genome_data) + 1)

    genome_data = get_color_discreet(genome_data, parameter, values)

    color_discrete_map = dict(zip(values, values_colors))

    fig = px.line(genome_data,
                    x = "x",
                    y = "y",
                    color = "legend",
                    color_discrete_map = {"Other": "darkgrey", "Background": "lightgrey", **color_discrete_map},
                    hover_name = "Feature_name")

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
    """Create a column "colors" used by Plotly as a reference for the color_discrete_map.

    Parameters
    ----------
    genome_data : Pandas dataframe
        2D coordinates of all the loci for Plotly visualization.
    parameter : str
        The name of the genome_data column containing the coloring parameter.
    values : list
        The values in genome_data[parameter] that will be colored according to values_colors.

    Returns
    -------
    Pandas dataframe
    """
    genome_data.loc[genome_data[parameter] != values[0], "legend"] = "Other"

    for v in values :
        genome_data.loc[genome_data[parameter] == v, "legend"] = v

    genome_data.loc[genome_data["Chromosome"] == 0, "legend"] = "Background"

    return genome_data