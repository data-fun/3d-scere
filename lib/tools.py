import pandas as pd
import sqlite3
import dash_html_components as html
import dash_table
import base64
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