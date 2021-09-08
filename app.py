import sys

sys.path.insert(0, "./lib/")

import tools
import visualization_2D as vis2D
import visualization_3D as vis3D
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import networkx as nx

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
from dash.dependencies import Input, Output, State

########################
############APP_INITIALIZATION############
########################

NAME = "3D-Scere"
FONTAWESOME = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css"
LITERA = "https://cdn.jsdelivr.net/npm/bootswatch@4.5.2/dist/litera/bootstrap.min.css"

GO_terms = pd.read_csv("./static/GO_terms.csv")
GO_terms_options = [{"label": GO, "value": GO} for GO in GO_terms["GO_terms"]]

plotly_segments = pd.read_csv("./static/plotly_segments.csv")
edges_list = pd.read_parquet("./static/3D_distances.parquet.gzip", engine="pyarrow")

# Get all features for all gene
sql_query = \
"""SELECT Primary_SGDID, Standard_gene_name, Chromosome, Feature_name, Strand, Stop_coordinate, Start_coordinate, Description
FROM SGD_features
"""
all_feature_name = tools.get_locus_info("./static/SCERE.db", sql_query)

#3D distance histogram constants
bin_number = 50
all_x = list(edges_list["3D_distances"])
H2, X = np.histogram(all_x, bins=bin_number, range=(0, 200))
H2 = H2/len(all_x)
F2 = np.cumsum(H2)/sum(H2)

basic_stylesheet = [{"selector": "node", "style": {"background-color": "#BFD7B5"}},
                    {"selector": "node", "style": {"label": "data(label)"}}]

colors = ["darkred", "red", "darkorange", "orange", "gold", "green",
"mediumseagreen", "turquoise", "deepskyblue", "dodgerblue",
"blueviolet", "purple", "magenta", "deeppink", "crimson", "black"]

app = dash.Dash(name=NAME, assets_folder="./assets", external_stylesheets=[dbc.themes.LUX, LITERA])
app.title = NAME
app.config.suppress_callback_exceptions = True

########################
############DASHBOARD_LAYOUT############
########################

############APP_HEADER############

header = html.Div(
        [dbc.Row(
            [
                html.Img(src="./static/yeast_icon.png", height="70px"),
                html.H1("3D-Scere", style={"padding-left": "2%", "padding-top": "1%"})
            ])
        ],
        style = {"padding-down": "4%", "padding-top": "2%"})

summary = html.Details([html.Summary([html.H4("Introduction")]),
                        html.Div("""3D-Scere is an open-source online tool for interactive visualization and exploration.
                                    This tool allows the visualization of any list of genes in the context of the 3D model of S. cerevisiae genome.
                                    Further information can easily be added like functional annotations (GO terms) or gene expression measurements.
                                    Qualitative or quantitative functional properties are thus highlighted in the large-scale 3D context of the genome
                                    with only a few mouse clicks.""")
                        ], open = True)

############APP_INPUTS_COMPONENTS############

input_tab1 = html.Div(
        [   dbc.Row(
            [
                dbc.Col(
                [
                    dbc.Row([html.H4("csv file upload", style={"padding-right" : "2%", "padding-left" : "2%"}), html.Abbr("\u003f\u20dd", title="Upload a one column .csv file with YORF")]),
                    dcc.Upload(id="upload_data_tab1", children=html.Div(
                    ["Drag and Drop or ",
                     html.A("Select Files")
                    ]),
                    style={"height": "60px",
                    "lineHeight": "60px",
                    "borderWidth": "1px",
                    "borderStyle": "dashed",
                    "borderRadius": "5px",
                    "textAlign": "center",
                    "margin": "10px"},
                    multiple=True),
                    dcc.Loading(children=[html.Div(id="output_data_upload_tab1")]),
                ]),
                dbc.Col(
                [
                ])
            ]),
            dbc.Row(style={"height" : 35}),
            dbc.Row(
            [
                dbc.Col(
                [
                    dbc.Row([html.H4("GO terms", style={"padding-right" : "2%", "padding-left" : "2%"}), html.Abbr("\u003f\u20dd", title="Choose a GO term to tag")]),
                    dcc.Dropdown(
                        id="GoTerm-dropdown",
                        options=GO_terms_options,
                        placeholder="select a GO term"),
                ]),
                dbc.Col(
                [
                    dbc.Row([html.H4("Color", style={"padding-right" : "2%", "padding-left" : "2%"}), html.Abbr("\u003f\u20dd", title="Choose the tagging color of the GO term")]),
                    dcc.Dropdown(
                        id="color-dropdown",
                        options=[
                            {"label": "Blue", "value": "blue"},
                            {"label": "Red", "value": "red"},
                            {"label": "Green", "value": "green"},
                            {"label": "Yellow", "value": "yellow"}],
                        placeholder="select a color"),
                ]),
            ]),
            dbc.Row(style={"height" : 25}),
            dbc.Row(
            [
                dbc.Button("Submit", id="Submit_tab1", outline=True, color="primary", className="mr-1", style={"vertical-align": "middle"})
            ],
            justify="end"
            )
        ],
        className="shadow p-3 mb-5 bg-body rounded", style={"padding-top" : "1%"})

input_tab2 = html.Div(
        [   dbc.Row(
            [
                dbc.Col(
                [
                    dbc.Row([html.H4("csv file upload", style={"padding-right" : "2%", "padding-left" : "2%"}),
                    html.Abbr("\u003f\u20dd", title="Upload a .csv file with YORF in the first column")]),
                    dcc.Upload(id="upload_data_tab2", children=html.Div(
                    ["Drag and Drop or ",
                     html.A("Select Files")
                    ]),
                    style={"height": "60px",
                    "lineHeight": "60px",
                    "borderWidth": "1px",
                    "borderStyle": "dashed",
                    "borderRadius": "5px",
                    "textAlign": "center",
                    "margin": "10px"},
                    multiple=True),
                    dcc.Loading(children=[html.Div(id="output_data_upload_tab2")]),
                ]),
                dbc.Col(
                [
                    dbc.Row([html.H4("Color scale", style={"padding-right" : "2%", "padding-left" : "2%"}), html.Abbr("\u003f\u20dd", title="Choose a colorscale")]),
                    dcc.Dropdown(
                        id="color_scale_dropdown",
                        options=[
                            {"label": "rainbow (diverging scale)", "value": "Rainbow"},
                            {"label": "picnic (diverging scale)", "value": "Picnic"},
                            {"label": "viridis", "value": "Viridis"},
                            {"label": "plasma", "value": "Plasma"},
                            {"label": "thermal", "value": "thermal"}],
                            placeholder="select a color scale"),
                ])
            ]),
            dbc.Row(style={"height" : 25}),
            dbc.Row(
            [
                dbc.Button("Submit", id="Submit_tab2", outline=True, color="primary", className="mr-1", style={"vertical-align": "middle"})
            ],
            justify="end"
            )
        ],
        className="shadow p-3 mb-5 bg-body rounded", style={"padding-top" : "1%"})

input_tab3 = html.Div(
        [   dbc.Row(
            [
                dbc.Col(
                [
                    dbc.Row([html.H4("csv file upload", style={"padding-right" : "2%", "padding-left" : "2%"}), html.Abbr("\u003f\u20dd", title="Upload a one column .csv file with YORF")]),
                    dcc.Upload(id="upload_data_tab3", children=html.Div(
                    ["Drag and Drop or ",
                     html.A("Select Files")
                    ]),
                    style={"height": "60px",
                    "lineHeight": "60px",
                    "borderWidth": "1px",
                    "borderStyle": "dashed",
                    "borderRadius": "5px",
                    "textAlign": "center",
                    "margin": "10px"},
                    multiple=True),
                    dcc.Loading(children=[html.Div(id="output_data_upload_tab3")]),
                ]),
                dbc.Col(
                [
                ]),
            ]),
            dbc.Row(style={"height" : 25}),
            dbc.Row(
            [
                dbc.Button("Submit", id="Submit_tab3", outline=True, color="primary", className="mr-1", style={"vertical-align": "middle"})
            ],
            justify="end"
            )
        ],
        className="shadow p-3 mb-5 bg-body rounded", style={"padding-top" : "1%"})

slider_tab3 = html.Div(
        [   dbc.Row(
            [
                dbc.Col(
                [
                    html.H4("3D distance treshold"),
                    dcc.Slider(id="treshold_slider",
                               min=0,
                               max=10,
                               step=1,
                               value=5
                              )
                ]),
                dbc.Col(
                [
                ]),
            ])
        ],
        className="shadow p-3 mb-5 bg-body rounded", style={"padding-top" : "1%"})

############APP_VISUALIZATIONS_COMPONENTS############

visualization_tab1 = html.Div(
        [   dbc.Row(
            [
                dbc.Col(
                [
                    html.H3("2D Visualization"),
                    dcc.Loading(children=[dcc.Graph(id="2D_representation")]),
                ])
            ]),
            dbc.Row(
            [
                dbc.Col(
                [
                    html.H3("Chromosomes repartition"),
                    dcc.Loading(children=[dcc.Graph(id="Chromosomes_repartition")]),
                ])
            ]),
            dbc.Row(
            [
                dbc.Col(
                [
                    html.H3("3D Visualization"),
                    dcc.Loading(children=[dcc.Graph(id="3D_representation")]),
                ])
            ]),
            dbc.Row(
            [
                dbc.Col(
                [
                    dcc.Graph(id="3D_representation_chrom"),
                ])
            ])
        ],
        className="shadow p-3 mb-5 bg-body rounded", style={"padding-top" : "1%"})

visualization_tab2 = html.Div(
        [   dbc.Row(
            [
                dbc.Col(
                [
                    html.H3("3D Visualization"),
                    dcc.Loading(children=[dcc.Graph(id="3D_representation_tab2")]),
                ])
            ])
        ],
        className="shadow p-3 mb-5 bg-body rounded", style={"padding-top" : "1%"})

visualization_tab3_hist = html.Div(
        [   dbc.Row(
            [
                dbc.Col(
                [
                html.H3("3D distances histogram"),
                dbc.Row(style={"height" : 10}),
                dcc.Loading(children=[html.Img(id="hist", src="")])
                ]),
            ])
        ],
        className="shadow p-3 mb-5 bg-body rounded", style={"padding-top" : "1%"})

visualization_tab3_network = html.Div(
        [   dbc.Row(
            [
                dbc.Col(
                [
                    html.H3("Network visualization"),
                    dcc.Loading(children=[cyto.Cytoscape(id="network",
                                                            stylesheet=basic_stylesheet,
                                                            elements=[],
                                                            style={"width": "100%", "height": "400px"},
                                                            layout={"name": "random"})])
                ])
            ])
        ],
        className="shadow p-3 mb-5 bg-body rounded", style={"padding-top" : "1%"})

visualization_tab3_metrics = html.Div(
        [   dbc.Row(
            [
                dbc.Col(
                [
                    html.H3("Network metrics"),
                    html.Div(id="output_edges_number_tab3"),
                    html.Div(id="output_nodes_number_tab3"),
                    dcc.Loading(children=[dcc.Graph(id="Degrees_hist")])
                ])
            ])
        ],
        className="shadow p-3 mb-5 bg-body rounded", style={"padding-top" : "1%"})


############APP_LAYOUT############

app.layout = dbc.Container(
      [ header,
        dbc.Row(style={"height" : 25}),
        summary,
        dbc.Row(style={"height" : 25}),
        dcc.Tabs([
        dcc.Tab(label="GO term projection", children=[
            dbc.Row(style={"height" : 45}),
            html.Div("""The projected list of genes can be colored uniformly or according to a selected GO term.
                        Upload the genes list as a one column .csv file containing YORF, then click the submit button.
                        To simultaneously color genes in the list that are associated to a given GO term, select it with an associated color before clicking on submit."""),
            dbc.Row(style={"height" : 45}),
            input_tab1,
            visualization_tab1
        ]),
        dcc.Tab(label="Quantitative variable projection", children=[
            dbc.Row(style={"height" : 45}),
            html.Div("""The projected list of genes can be colored according to a given quantitative variable.
                        Upload the genes list as a .csv file, with YORF in the first column. Then select the column corresponding the quantitative variable and a color scale before
                        clicking on submit."""),
            dbc.Row(style={"height" : 45}),
            input_tab2,
            visualization_tab2
        ]),
        dcc.Tab(label="3D distances histogram and network", children=[
            dbc.Row(style={"height" : 45}),
            html.Div("""All the 3D distances between genes in the list are summarized into a histogram and a network.
                        Upload the genes list as a one column .csv file containing YORF, then click the submit button.
                        The slider determines the threshold under witch 3D distances are used to construct the network."""),
            dbc.Row(style={"height" : 45}),
            input_tab3,
            slider_tab3,
            visualization_tab3_hist,
            visualization_tab3_network,
            visualization_tab3_metrics
        ]),
        ])
      ])

########################
############CALLBACKS############
########################

############TAB1_UPLOAD############
@app.callback(Output("output_data_upload_tab1", "children"),
              Input("upload_data_tab1", "contents"),
              State("upload_data_tab1", "filename"))
def update_output(list_of_contents, list_of_names):
    if list_of_contents is not None:
        children=[
            tools.parse_contents(c, n, "datatable_tab1") for c, n in
            zip(list_of_contents, list_of_names)]
        return children

############TAB1_UPLOAD_STYLE############
@app.callback(
    Output("datatable_tab1", "style_data_conditional"),
    Input("datatable_tab1", "selected_columns"))
def update_styles_tab1(selected_columns):
    return [{
        "if": { "column_id": i },
        "background_color": "#D2F3FF"
    } for i in selected_columns]

############TAB1_2D_GRAPH############
@app.callback(Output("2D_representation", "figure"),
              Input("Submit_tab1", "n_clicks"),
              State("GoTerm-dropdown", "value"),
              State("color-dropdown", "value"),
              State("datatable_tab1", "derived_virtual_data"),
              State("datatable_tab1", "selected_columns"))
def update_2D_graphs_tab1(n_clicks, GoTerm, color, data, column):

    sql_query_gobal = \
"""SELECT Primary_SGDID, count(SGDID), Feature_name, Start_coordinate, Stop_coordinate, Chromosome, Strand, GO_slim_term
FROM SGD_features, go_slim_mapping
WHERE SGDID == Primary_SGDID
GROUP BY SGDID
ORDER BY Start_coordinate
"""
    sql_query_specific = \
"""SELECT Primary_SGDID, count(SGDID), Feature_name, Start_coordinate, Stop_coordinate, Chromosome, Strand, GO_slim_term
FROM SGD_features, go_slim_mapping
WHERE SGDID == Primary_SGDID
AND (GO_slim_term == """ + "'" + str(GoTerm) + "'" + """)
GROUP BY SGDID
ORDER BY Start_coordinate
"""
    all_loci = tools.get_locus_info("./static/SCERE.db", sql_query_gobal)
    selected_loci = tools.get_locus_info("./static/SCERE.db", sql_query_specific)

    loci = pd.concat([all_loci, selected_loci]).drop_duplicates(subset=["Primary_SGDID"], keep="last")

    if column != []:
        unfiltered_data = pd.DataFrame(data)
        filtered_data = unfiltered_data[str(column[0])]
        loci = loci.assign(FT_target=loci.Feature_name.isin(filtered_data))

        loci.loc[loci.FT_target == True, "colors_parameters"]="Targets"
        loci.loc[(loci.GO_slim_term == str(GoTerm)) & (loci.FT_target == True), "colors_parameters"]=str(GoTerm)

        loci = vis2D.format_coordinates(loci, 6)
        fig = vis2D.genome_drawing(loci, "colors_parameters", [str(GoTerm), "Targets"], [str(color), "Black"])

    else :
        loci = vis2D.format_coordinates(loci, 6)
        fig = vis2D.genome_drawing(loci, "GO_slim_term", [str(GoTerm)], [str(color)])

    return fig

############TAB1_CHROMOSOME_REPARTITION############
@app.callback(Output("Chromosomes_repartition", "figure"),
              Input("Submit_tab1", "n_clicks"),
              State("datatable_tab1", "derived_virtual_data"),
              State("datatable_tab1", "selected_columns"))
def update_chrom_repartition_tab1(n_clicks, data, column):

    sql_query_2 = \
"""SELECT Primary_SGDID, Feature_name, Start_coordinate, Stop_coordinate, Chromosome, Strand
FROM SGD_features
ORDER BY Start_coordinate
"""
    if column != []:
        unfiltered_data = pd.DataFrame(data)
        filtered_data = unfiltered_data[str(column[0])]

        loci = tools.get_locus_info("./static/SCERE.db", sql_query_2)
        loci = loci.assign(FT_target=loci.Feature_name.isin(filtered_data))

        loci = loci[loci.FT_target == True].drop(["FT_target"], axis=1)

        fig = px.histogram(loci, x="Chromosome", nbins=30, range_x=[1, 17], color_discrete_sequence=["#A0E8AF"])
        fig.update_layout(plot_bgcolor="white",
                        xaxis_showgrid=False,
                        yaxis_showgrid=False,
                        showlegend=True)

        return fig

############TAB1_3D_GRAPH_FEATURE############
@app.callback(Output("3D_representation", "figure"),
              Input("Submit_tab1", "n_clicks"),
              State("GoTerm-dropdown", "value"),
              State("color-dropdown", "value"),
              State("datatable_tab1", "derived_virtual_data"),
              State("datatable_tab1", "selected_columns"))
def update_3D_graph_tab1(n_clicks, GoTerm, color, data, column):

    sql_query_gobal = \
"""SELECT Primary_SGDID, count(SGDID), Feature_name, Start_coordinate, Stop_coordinate, Chromosome, Strand, GO_slim_term
FROM SGD_features, go_slim_mapping
WHERE SGDID == Primary_SGDID
GROUP BY SGDID
ORDER BY Start_coordinate
"""

    sql_query_3 = \
"""SELECT Primary_SGDID, Feature_name, Start_coordinate, Stop_coordinate, Chromosome, Strand, GO_slim_term
FROM SGD_features, go_slim_mapping
WHERE SGDID == Primary_SGDID
AND (GO_slim_term == """ + "'" + str(GoTerm) + "'" + """)
GROUP BY SGDID
ORDER BY Start_coordinate
"""
    all_loci = tools.get_locus_info("./static/SCERE.db", sql_query_gobal)
    selected_loci = tools.get_locus_info("./static/SCERE.db", sql_query_3)

    if column != []:
        unfiltered_data = pd.DataFrame(data)
        filtered_data = unfiltered_data[str(column[0])]
        loci = all_loci.assign(FT_target=all_loci.Feature_name.isin(filtered_data))
        loci = loci.assign(GoTerm=loci.Primary_SGDID.isin(selected_loci.Primary_SGDID))

        loci.loc[loci.FT_target == True, "colors_parameters"]="Targets"
        loci.loc[(loci.GoTerm == True) & (loci.FT_target == True), "colors_parameters"]=str(GoTerm)

        loci_segments = plotly_segments.merge(loci, on="Primary_SGDID", how="left", copy=False)
        loci_segments.index = range(1, len(loci_segments) + 1)
        loci_segments = vis3D.get_color_discreet_3D(loci_segments, "colors_parameters", [str(GoTerm), "Targets"], [str(color), "blue"])
        fig = vis3D.genome_drawing(loci_segments)

    else :
        selected_loci_segments = plotly_segments.merge(selected_loci, on="Primary_SGDID", how="left", copy=False)
        selected_loci_segments.index = range(1, len(selected_loci_segments) + 1)
        selected_loci_segments = vis3D.get_color_discreet_3D(selected_loci_segments, "GO_slim_term", [str(GoTerm)], [str(color)])

        fig = vis3D.genome_drawing(selected_loci_segments)

    return fig

############TAB1_3D_GRAPH_CHROMOSOMES############
@app.callback(Output("3D_representation_chrom", "figure"),
              Input("Submit_tab1", "n_clicks"))
def update_3D_graph_chrom_tab1(n_clicks):

    sql_query_4 = \
"""SELECT Primary_SGDID, Start_coordinate, Stop_coordinate, Chromosome, Strand
FROM SGD_features
ORDER BY Start_coordinate
"""
    selected_loci = tools.get_locus_info("./static/SCERE.db", sql_query_4)

    selected_loci_segments = plotly_segments.merge(selected_loci, on="Primary_SGDID", how="left", copy=False)
    selected_loci_segments.index = range(1, len(selected_loci_segments) + 1)

    selected_loci_segments = vis3D.get_color_discreet_3D(selected_loci_segments, "Chromosome", list(range(1, 17)), colors)

    return vis3D.genome_drawing(selected_loci_segments)

############TAB2_UPLOAD############
@app.callback(Output("output_data_upload_tab2", "children"),
              Input("upload_data_tab2", "contents"),
              State("upload_data_tab2", "filename"))
def update_output_tab2(list_of_contents, list_of_names):
    if list_of_contents is not None:
        children=[
            tools.parse_contents(c, n, "datatable") for c, n in
            zip(list_of_contents, list_of_names)]
        return children

############TAB2_COLUMN_SELECTION_UPLOAD############
@app.callback(
    Output("datatable", "style_data_conditional"),
    Input("datatable", "selected_columns"))
def update_styles_tab2(selected_columns):
    return [{
        "if": { "column_id": i },
        "background_color": "#D2F3FF"
    } for i in selected_columns]

############TAB2_3D_GRAPH############
@app.callback(Output("3D_representation_tab2", "figure"),
              Input("Submit_tab2", "n_clicks"),
              State("datatable", "derived_virtual_data"),
              State("datatable", "selected_columns"),
              State("color_scale_dropdown", "value"))
def update_3D_graphs_tab2(n_clicks, input1, input2, input3):

    unfiltered_data = pd.DataFrame(input1)
    filtered_data = unfiltered_data[[str(input2[0]), str(input2[1])]]

    sql_query_5 = \
"""SELECT Primary_SGDID, Start_coordinate, Stop_coordinate, Chromosome, Feature_name, Strand
FROM gene_literature, SGD_features
WHERE SGDID == Primary_SGDID
GROUP BY SGDID
ORDER BY Start_coordinate
"""

    whole_genome = tools.get_locus_info("./static/SCERE.db", sql_query_5)

    whole_genome_segments = plotly_segments.merge(whole_genome, on="Primary_SGDID", how="left", copy=False)
    whole_genome_segments.index = range(1, len(whole_genome_segments) + 1)

    whole_genome_segments = whole_genome_segments.merge(filtered_data, left_on="Feature_name", right_on="YORF", how="left", copy=False)
    whole_genome_segments.iloc[: , -1].fillna("whitesmoke", inplace=True)

    fig = go.Figure(data=[go.Scatter3d(x=whole_genome_segments.x,
                                   y=whole_genome_segments.y,
                                   z=whole_genome_segments.z,
                                   mode="lines",
                                   name="",
                                   line={"color": whole_genome_segments.iloc[: , -1],
                                           "colorscale": input3,
                                           "showscale": True,
                                           "width": 12},
                                   customdata=whole_genome_segments.Primary_SGDID,
                                   hovertemplate=("<b>SGDID :</b> %{customdata} <br>"
                                                    "<b>x :</b> %{x} <br>"),
                                   hoverlabel=dict(bgcolor="white", font_size=16))])

    fig.update_layout(scene=dict(xaxis=dict(showgrid=False, backgroundcolor="white"),
                             yaxis=dict(showgrid=False, backgroundcolor="white"),
                             zaxis=dict(showgrid=False, backgroundcolor="white")))
    fig.update_layout(height=800)

    return fig

############TAB3_UPLOAD############
@app.callback(Output("output_data_upload_tab3", "children"),
              Input("upload_data_tab3", "contents"),
              State("upload_data_tab3", "filename"))
def update_output_tab3(list_of_contents, list_of_names):
    if list_of_contents is not None:
        children=[
            tools.parse_contents(c, n, "datatable_tab3") for c, n in
            zip(list_of_contents, list_of_names)]
        return children

############TAB3_UPLOAD_STYLE############
@app.callback(
    Output("datatable_tab3", "style_data_conditional"),
    Input("datatable_tab3", "selected_columns"))
def update_styles_tab3(selected_columns):
    return [{
        "if": { "column_id": i },
        "background_color": "#D2F3FF"
    } for i in selected_columns]

############TAB3_SLIDER_AND_NETWORK############
@app.callback(Output("network", "elements"),
              Output("treshold_slider", "max"),
              Output("treshold_slider", "min"),
              Input("Submit_tab3", "n_clicks"),
              State("datatable_tab3", "derived_virtual_data"))
def update_network(n_clicks, input1):

    genes_list = pd.DataFrame(input1)

    sql_query_6 = \
"""SELECT Primary_SGDID, Chromosome, Feature_name, Strand, Stop_coordinate, Start_coordinate
FROM SGD_features
"""

    Feature_name = tools.get_locus_info("./static/SCERE.db", sql_query_6)
    Feature_name = Feature_name.merge(genes_list, left_on="Feature_name", right_on=genes_list.columns[0])

    nodes = [{"data": {"id": Primary_SGDID, "label": Feature_name}}
         for Primary_SGDID, Feature_name in zip(Feature_name["Primary_SGDID"], Feature_name["Feature_name"])
        ]

    edges_list_select = tools.get_edges_list(genes_list, edges_list, all_feature_name)

    edges = [{"data": {"source": source, "target": target, "weight": float(weight)}}
             for source, target, weight in zip(edges_list_select["Primary_SGDID_bis"], edges_list_select["Primary_SGDID"], edges_list_select["3D_distances"])
            ]

    elements = nodes + edges
    slider_max = max(edges_list_select["3D_distances"])
    slider_min = min(edges_list_select["3D_distances"])

    return elements, slider_max, slider_min

############TAB3_HIST############
@app.callback(Output("hist", component_property="src"),
              Input("Submit_tab3", "n_clicks"),
              Input("treshold_slider", "value"),
              State("datatable_tab3", "derived_virtual_data"))
def update_hist(n_clicks, input1, input2):

    genes_list = pd.DataFrame(input2)

    fig = tools.distri(genes_list, edges_list, all_feature_name, H2, F2, bin_number, input1)

    out_url = tools.fig_to_uri(fig)

    return out_url

############TAB3_NETWORK_TRESHOLD############
@app.callback(Output("network", "stylesheet"),
              Input("treshold_slider", "value"))
def update_stylesheet_(treshold):
    new_styles = [{"selector": "[weight >" + str(treshold) + "]", "style": {"opacity": 0}}]
    stylesheet = basic_stylesheet + new_styles

    return stylesheet

############TAB3_NETWORK_METRICS############
@app.callback(Output("output_nodes_number_tab3", "children"),
              Input("treshold_slider", "value"),
              Input("network", "elements"))
def update_metrics_1(treshold, elements):

    subgraph_edges = pd.DataFrame(elements)
    subgraph_edges = pd.json_normalize(subgraph_edges["data"])
    subgraph_edges = subgraph_edges[subgraph_edges["weight"] < treshold]

    G = nx.from_pandas_edgelist(subgraph_edges, source="source", target="target")

    return "number of nodes : " + str(G.number_of_nodes())

@app.callback(Output("output_edges_number_tab3", "children"),
              Input("treshold_slider", "value"),
              Input("network", "elements"))
def update_metrics_2(treshold, elements):

    subgraph_edges = pd.DataFrame(elements)
    subgraph_edges = pd.json_normalize(subgraph_edges["data"])
    subgraph_edges = subgraph_edges[subgraph_edges["weight"] < treshold]

    G = nx.from_pandas_edgelist(subgraph_edges, source="source", target="target")

    return "number of edges : " + str(G.number_of_edges())

@app.callback(Output("Degrees_hist", "figure"),
              Input("treshold_slider", "value"),
              Input("network", "elements"))
def update_metrics_3(treshold, elements):

    subgraph_edges = pd.DataFrame(elements)
    subgraph_edges = pd.json_normalize(subgraph_edges["data"])
    subgraph_edges = subgraph_edges[subgraph_edges["weight"] < treshold]

    G = nx.from_pandas_edgelist(subgraph_edges, source="source", target="target")

    degrees = [val for (node, val) in G.degree()]
    fig = px.histogram(degrees, nbins= 70, color_discrete_sequence=["#A0E8AF"])
    fig.update_layout(plot_bgcolor="white",
                      xaxis_showgrid=False,
                      yaxis_showgrid=False,
                      showlegend=True)

    return fig


if __name__ == "__main__":
    app.run_server(debug=False)
