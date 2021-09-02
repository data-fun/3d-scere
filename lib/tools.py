import math

import numpy as np
import pandas as pd
import sqlite3


def display_module_version():
    """Display dependencies versions.
    """
    print("sqlite3 version:", sqlite3.version)
    print("pandas version:", pd.__version__)
    print("numpy version:", np.__version__)


def string_to_list(string):
    return list(string.split(", "))


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



