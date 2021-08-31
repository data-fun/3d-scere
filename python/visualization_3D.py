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

def format_segments_dataframe(chrom_number, atoms_coordinates):
    
    segments_coordinates = calcul_segments_coordinates(chrom_number, atoms_coordinates)
    
    segments_sizes = segments_coordinates.apply(calcul_segment_size, axis = 1)
    sum_segments_sizes = sum(segments_sizes)
    
    chrom_lenght = get_chromosome_lenght(chrom_number)
    
    sizes_on_chromosome = segments_sizes.apply(calcul_size_on_chrom, args = [chrom_lenght, sum_segments_sizes])
    
    start_bp = [0] + list(sizes_on_chromosome.cumsum().iloc[:len(sizes_on_chromosome) - 1])
    
    segments = pd.DataFrame(data = {"Segment_ID": sizes_on_chromosome.index,
                                              "start_bp": start_bp,
                                              "stop_bp": sizes_on_chromosome.cumsum()})
    
    return segments

def calcul_segments_coordinates(chrom_number, atoms_coordinates):
    segment_start = atoms_coordinates[atoms_coordinates["chrom"] == chrom_number]
    segment_start.index = range(1, len(segment_start) + 1)
    segment_start = segment_start.drop([len(segment_start)], axis = 0)
    segment_start.index = range(1, len(segment_start) + 1)

    segment_stop = atoms_coordinates[atoms_coordinates['chrom'] == chrom_number]
    segment_start.index = range(1, len(segment_start) + 1)
    segment_stop = segment_stop.iloc[1:]
    segment_stop.index = range(1, len(segment_stop) + 1)
    
    segments_coordinates = segment_start.merge(segment_stop, left_index = True,
                                               right_index = True,
                                               suffixes = ["_start", "_stop"] )
    segments_coordinates = segments_coordinates.drop(["chrom_start","chrom_stop"], axis = 1)

    return segments_coordinates

def calcul_segment_size(segments_coordinates):
    
    return math.sqrt((segments_coordinates.x_stop - segments_coordinates.x_start) ** 2 +
                     (segments_coordinates.y_stop - segments_coordinates.y_start) ** 2 +
                     (segments_coordinates.z_stop - segments_coordinates.z_start) ** 2)

def get_chromosome_lenght(chrom_number):
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

def calcul_size_on_chrom(segments_sizes, chrom_lenght, sum_segments_sizes):
    
    return (segments_sizes * chrom_lenght) / sum_segments_sizes

# ###Segments and locus connexions.

def connect(segments, loci):
    loci = loci[loci.Feature_type != 'CDS']
    loci = loci[loci.Start_coordinate != ""]
    loci.index = range(1, len(loci) + 1)
    
    intersections = segments.apply(get_intersection, args = [loci], axis = "columns")
    
    return intersections

def get_intersection(segment, loci):
    
    on = loci.assign(start_bp = segment.start_bp)
    on = on.assign(stop_bp = segment.stop_bp)
    
    conditions = [(on.Start_coordinate <= on.start_bp) & (on.start_bp <= on.Stop_coordinate),
                  (on.start_bp <= on.Start_coordinate) & (on.Start_coordinate <= on.stop_bp)]
    
    choices = [True, True]
    
    on = on.assign(on_segment = np.select(conditions, choices, default = False))
    
    locus_on_segment = on.Primary_SGDID[on.on_segment == True]
    locus_on_segment.index = range(1, len(locus_on_segment) + 1)
    
    return locus_on_segment

def string_to_list(string):
    return list(string.split(", "))

def get_chrom_info(database, chrom_number, info, source, condition, group_by, order_by):
    #SQL request
    db_connexion = sqlite3.connect(database)
    
    cursor = db_connexion.cursor()
    
    chrom = cursor.execute("""
    SELECT """ + info + """
    FROM """ + source + """ 
    WHERE Chromosome == """ + chrom_number + condition + group_by + """
    ORDER BY """ + order_by)
    
    
    list_info = string_to_list(info)
    
    #pandas dataframe formatting
    chrom = chrom.fetchall()
    chrom = pd.DataFrame(chrom, columns = list_info)
    
    return chrom

# ###Format atoms coordinates with segments

def format_atoms_coordinates_V2(atoms_coordinates, segments_loci):
    
    plotly_segments = pd.DataFrame(columns = ["chrom", "x", "y", "z"])
    row_null = {"chrom": 0, "x": "none", "y": "none", "z": "none"}
    start_index_minus = 0
    stop_index_minus = 1

    for c in range(1, 17):
        chrom_coordinates = atoms_coordinates[atoms_coordinates["chrom"] == c]
        chrom_segments_loci = segments_loci.loc[chrom_coordinates.index[0] - start_index_minus:chrom_coordinates.index.max() - stop_index_minus]
        start_index_minus = start_index_minus + 1
        stop_index_minus = stop_index_minus + 1
        
        chrom_coordinates = chrom_coordinates.merge(chrom_segments_loci, left_index = True, right_index = True)
        chrom_coordinates = chrom_coordinates.rename(columns = {1: "Primary_SGDID"})
        
        row_one = chrom_coordinates.copy()
        
        row_two = chrom_coordinates.copy()
        row_two = row_two[1:]
        row_one.index = range(2, len(row_one) + 2)
        row_two["Primary_SGDID"] = row_one["Primary_SGDID"]
        
        row_three = chrom_coordinates.assign(chrom = 0, x = "none", y = "none", z = "none")
        
        row_one.index = range(0, len(chrom_coordinates) * 3, 3)
        row_two.index = range(1, (len(chrom_coordinates)-1) * 3, 3)
        row_three.index = range(2, len(chrom_coordinates) * 3, 3)
        
        row_one = row_one.transpose()
        row_two = row_two.transpose()
        row_three = row_three.transpose()
        
        plotly_chrom = pd.merge(row_one, row_two, how = "left", left_index = True, right_index = True)
        plotly_chrom = pd.merge(plotly_chrom, row_three, how = "left", left_index = True, right_index = True)
        plotly_chrom = plotly_chrom.transpose()
        plotly_chrom = plotly_chrom.sort_index()
        
        plotly_chrom = plotly_chrom.append(row_null, ignore_index = True)
        
        plotly_segments = plotly_segments.append(plotly_chrom)
        plotly_segments.index = range(1, len(plotly_segments) + 1)

    return plotly_segments

# ###Adding colors in 3D

def get_color_discreet_3D(genome_data, parameter, values, values_colors):
    
    genome_data.loc[genome_data[parameter] != values[0], "colors"] = "darkgrey"
    
    for v, c in zip(values, values_colors):
        genome_data.loc[genome_data[parameter] == v, "colors"] = c
        
    genome_data.loc[genome_data["Primary_SGDID"].isna() == True, "colors"] = "whitesmoke"

    return genome_data

# ###3D Genome drawing.

def genome_drawing(whole_genome_segments):

    fig = go.Figure(data=[go.Scatter3d(x = whole_genome_segments.x,
                                       y = whole_genome_segments.y,
                                       z = whole_genome_segments.z,
                                       mode = "lines",
                                       name = "",
                                       line = {"color": whole_genome_segments["colors"],
                                               "width": 12},
                                       customdata = whole_genome_segments.Primary_SGDID,
                                       hovertemplate = ("<b>SGDID :</b> %{customdata} <br>"
                                                        "<b>x :</b> %{x} <br>"),
                                       hoverlabel = dict(bgcolor = "white", font_size = 16))])

    fig.update_layout(scene=dict(xaxis = dict(showgrid = False, backgroundcolor = "white"),
                                 yaxis = dict(showgrid = False, backgroundcolor = "white"),
                                 zaxis = dict(showgrid = False, backgroundcolor = "white")))
    fig.update_layout(height=800)
    
    return fig





