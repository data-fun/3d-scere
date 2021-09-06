import numpy as np
import pandas as pd
import sqlite3


def display_module_version():
    """Display dependencies versions.
    """
    print("sqlite3 version:", sqlite3.version)
    print("pandas version:", pd.__version__)
    print("numpy version:", np.__version__)


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