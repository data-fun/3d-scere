import plotly.graph_objects as go


#3D Genome drawing.

def genome_drawing(whole_genome_segments):
    """Draw the 3D plotly figure, representing the 16 chromosomes in lightgrey and all the loci in darkgrey.

    Parameters
    ----------
    whole_genome_segments : Pandas dataframe
        3D segments coordinates for Plotly visualization and their associated locus.

    Returns
    -------
    Plotly figure
    """
    fig = go.Figure(data=[go.Scatter3d(x = whole_genome_segments.x,
                                       y = whole_genome_segments.y,
                                       z = whole_genome_segments.z,
                                       mode = "lines",
                                       name = "",
                                       line = {"color": whole_genome_segments["colors"],
                                               "width": 12},
                                       customdata = whole_genome_segments["Feature_name"],
                                       hovertemplate = ("<b>YORF :</b> %{customdata} <br>"),
                                       hoverlabel = dict(bgcolor = "white", font_size = 16))])

    fig.update_layout(scene=dict(xaxis = dict(showgrid = False, backgroundcolor = "white"),
                                 yaxis = dict(showgrid = False, backgroundcolor = "white"),
                                 zaxis = dict(showgrid = False, backgroundcolor = "white")))
    fig.update_layout(height=800)

    return fig

#Adding colors in 3D

def get_color_discreet_3D(genome_data, parameter, values, values_colors):
    """Create a column "colors" used by Plotly to color each 3D segment.

    Parameters
    ----------
    genome_data : Pandas dataframe
        2D coordinates of all the loci for Plotly visualization.
    parameter : str
        The name of the genome_data column containing the coloring parameter.
    values : list
        The values in genome_data[parameter] that will be colored according to values_colors.
    values_colors : list
        A list of color names.

    Returns
    -------
    Pandas dataframe
    """

    genome_data.loc[genome_data[parameter] != values[0], "colors"] = "darkgrey"

    for v, c in zip(values, values_colors):
        genome_data.loc[genome_data[parameter] == v, "colors"] = c

    genome_data.loc[genome_data["Primary_SGDID"].isna() == True, "colors"] = "whitesmoke"

    return genome_data