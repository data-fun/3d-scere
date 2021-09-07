import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sqlite3
import dash_html_components as html
import dash_table
import base64
from io import BytesIO
import io


def display_module_version():
    """Display dependencies versions.
    """
    print("sqlite3 version:", sqlite3.version)
    print("pandas version:", pd.__version__)

def get_locus_info(database, query):
    """Query the SQLite database.
    
    Parameters
    ----------
    database : str
        Path to the SQLite database.
    query : str
        SQL query.

    Returns
    -------
    Pandas Dataframe
    """  
    # Connect to database.
    db_connexion = sqlite3.connect(database)
    cursor = db_connexion.cursor()
    
    # Query database.
    chrom_info = cursor.execute(query)
    
    # Convert to Pandas dataframe
    column_names = [column[0] for column in chrom_info.description]
    chrom_info_df = pd.DataFrame(chrom_info.fetchall(), columns=column_names)
    
    # Select only strands + and -
    chrom_info_df = chrom_info_df[ (chrom_info_df["Strand"] == "C") | (chrom_info_df["Strand"] == "W") ]
    # Remove "2-micron" plasmid
    chrom_info_df = chrom_info_df[ chrom_info_df["Chromosome"] != "2-micron" ]
    # Convert chromosome id to int
    chrom_info_df["Chromosome"] = chrom_info_df["Chromosome"].astype(int)

    return chrom_info_df

############UPLOAD_PARSING############

def parse_contents(contents, filename, datatable_id):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])
    
    return html.Div([
        html.H5(filename),
        dash_table.DataTable(
            id=datatable_id,
            data=df.to_dict('records'),
            columns=[{'name': i, 'id': i, "selectable": True} for i in df.columns],
            page_size=10,
            column_selectable="multi",
            selected_columns=[df.columns[0]],
            style_cell={'textAlign': 'left'},
            style_data_conditional=[{'if': {'row_index': 'odd'},
                                     'backgroundColor': 'rgb(248, 248, 248)'}],
            style_header={'backgroundColor': 'rgb(230, 230, 230)',
                          'fontWeight': 'bold'}),

        html.Hr(),  # horizontal line

        # For debugging, display the raw contents provided by the web browser
        #html.Div('Raw Content'),
        #html.Pre(contents[0:200] + '...', style={'whiteSpace': 'pre-wrap','wordBreak': 'break-all'})
    ])

def get_edges_list(gene_list, edges_list, feature_name):
    
    # Add SGDID
    feature_name = feature_name.merge(gene_list, left_on = "Feature_name", right_on = gene_list.columns[0])
    
    # Extract distances for selected genes list
    edges_list_select = edges_list[edges_list["Primary_SGDID"].isin(feature_name["Primary_SGDID"])]
    edges_list_select = edges_list_select[edges_list_select["Primary_SGDID_bis"].isin(feature_name["Primary_SGDID"])]
    edges_list_select.index = range(1, len(edges_list_select) + 1)
    
    return edges_list_select

def distri(genes_list, edges_list, feature_name, H2, F2, bin_number, input1):
    
    edges_list_select = get_edges_list(genes_list, edges_list, feature_name)
    x = list(edges_list_select["3D_distances"])
    H, X1 = np.histogram(x, bins = bin_number, range = (0, 200))
    F1 = np.cumsum(H)/len(x)

    H = H/len(x)

    fig, ax = plt.subplots()
    ax.hist(X1[:-1], X1, weights=H2, color="#5767FF", alpha=0.3, label="All distances")
    ax.set_xlim(0, 200)
    plt.xlabel("3D distances", size = 16)
    ax.hist(X1[:-1], X1, weights=H, color="#FA3824", alpha=0.3, label="TF's targets")
    plt.ylabel("Density", size = 16)
    ax2=ax.twinx()
    plt.ylabel("CDF", size = 16)
    ax2.plot(X1[:-1], F1, label="CDF (TF's targets)", color = "#FA3824")
    ax2.plot(X1[:-1], F2, label="CDF (all)", color = "#5767FF")
    ax.legend(bbox_to_anchor = (0.6, 0.9), loc="upper left")
    ax2.legend(bbox_to_anchor = (0.6, 0.7), loc="upper left")

    plt.axvline(x=input1, color='black', linestyle='--')
    
    return fig

#From https://github.com/4QuantOSS/DashIntro/blob/master/notebooks/Tutorial.ipynb
def fig_to_uri(in_fig, close_all=True, **save_args):
    """
    Save a figure as a URI
    """
    out_img = BytesIO()
    in_fig.savefig(out_img, format='png', **save_args)
    if close_all:
        in_fig.clf()
        plt.close('all')
    out_img.seek(0)  # rewind file
    encoded = base64.b64encode(out_img.read()).decode("ascii").replace("\n", "")
    return "data:image/png;base64,{}".format(encoded)