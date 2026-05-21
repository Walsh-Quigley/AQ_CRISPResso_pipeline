import pandas as pd
from pathlib import Path


# Reads CRISPResso ouput files: allele frequency tables, mapping statistics

def read_mapping_stats(path: Path) -> tuple[int, int]:
    """Opens the CRISPRessso_mapping_statistics file and collects the total and
    aligned read number
    Args: 
        path: Path to the CRISPResso_mapping_statstics file
    Returns:
        tuple[int, int]: returns a tuple containing total reads and aligned reads for a given sample
    Raises:
        ValueError: total reads value is 0 in the file
        ValueError: aligned reads exceeds total reads.
        FileNotFoundError: if the 'open' function fails
    """
    with open(path) as f:
        headers = f.readline().strip().split("\t")
        values = f.readline().strip().split("\t")
    
    row = dict(zip(headers,values))
    # row is now {"READS AFTER PREPROCESSING": "######",
    #             "READS ALIGNED": "##### }

    reads_total = int(row["READS AFTER PREPROCESSING"])
    reads_aligned = int(row["READS ALIGNED"])

    if reads_total == 0:
        raise ValueError(f"No reads found in {path} — file may be corrupt or empty")

    if reads_aligned > reads_total:
        raise ValueError(f"Reads aligned ({reads_aligned}) exceeds total reads ({reads_total}) in {path}")


    return reads_total, reads_aligned

def read_allele_table(path: Path) -> pd.DataFrame:
    """Creates a dataframe for the allele_frequency_table
    Args:
        path: file path to the allele_frequency_table
    Returns:
        pd.DataFrame: a pandas dataframe of the allele_frequency_table_data
    Raises:
        FileNotFoundError: if read_csv's open() call fails
    """
    df = pd.read_csv(path, sep="\t")
    return df

def read_quant_window(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df = df.astype(float)
    return df

