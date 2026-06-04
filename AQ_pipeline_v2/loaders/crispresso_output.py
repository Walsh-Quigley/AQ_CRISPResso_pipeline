import pandas as pd
from pathlib import Path


# Reads CRISPResso output files: allele frequency tables, mapping statistics

def read_mapping_stats(path: Path) -> tuple[int, int]:
    """Opens the CRISPResso_mapping_statistics file and collects the total and
    aligned read number
    Args: 
        path: Path to the CRISPResso_mapping_statistics file
    Returns:
        tuple[int, int]: returns a tuple containing total reads and aligned reads for a given sample
    Raises:
        ValueError: total reads value is 0 in the file
        ValueError: aligned reads exceeds total reads.
        ValueError: aligned reads below 10% of total reads
        FileNotFoundError: if the 'open' function fails
    """
    with open(path, encoding="utf-8") as f:
        headers = f.readline().strip().split("\t")
        values = f.readline().strip().split("\t")
    
    row = dict(zip(headers,values))
    # row contains all columns from the stats file:
    # {
    #     "READS IN INPUTS":         "######",
    #     "READS AFTER PREPROCESSING":  "######",
    #     "READS ALIGNED":           "######",
    #     "N_COMPUTED_ALN":          "######",
    #     "N_CACHED_ALN":            "######",
    #     "N_COMPUTED_NOTALN":       "######",
    #     "N_CACHED_NOTALN":         "######",
    # }
    # We only use READS AFTER PREPROCESSING and READS ALIGNED.

    reads_total = int(row["READS AFTER PREPROCESSING"])
    reads_aligned = int(row["READS ALIGNED"])

    if reads_total == 0:
        raise ValueError(f"No reads found in {path} — file may be corrupt or empty")

    if reads_aligned > reads_total:
        raise ValueError(f"Reads aligned ({reads_aligned}) exceeds total reads ({reads_total}) in {path}")

    pct_of_reads_aligned = (reads_aligned / reads_total) * 100

    if pct_of_reads_aligned < 10:
        raise ValueError(
            f"Reads aligned unusually low ({pct_of_reads_aligned:.1f}%): "
            f"{reads_aligned} out of {reads_total} in {path}"
        )


    return reads_total, reads_aligned

def read_allele_table(path: Path) -> pd.DataFrame:
    """Creates a dataframe for the allele_frequency_table
    Args:
        path: file path to the allele_frequency_table
    Returns:
        pd.DataFrame: a pandas dataframe of the allele_frequency_table_data
    Raises:
        FileNotFoundError: if read_csv's open() call fails
        ValueError: df.astype(float) interacts with non-numeric value
    """
    df = pd.read_csv(path, sep="\t")
    return df

def read_quant_window(path: Path) -> pd.DataFrame:
    """Reads the CRISPResso quantification to create a dataframe for downstream use
    Args:
        path: file path to the quantification table
    Returns:
        pd.DataFrame: a pandas dataframe of the quantification window
    Raises:
        FileNotFoundError: if read_csv's open() call fails"""
    df = pd.read_csv(path, sep="\t", index_col=0)
    df = df.astype(float)
    return df

def read_editing_frequency(path: Path) -> dict:
    """Reads the CRISPResso_quantification_of_editing_frequency file and returns the
    modification breakdown for the (single) reference amplicon.
    Args:
        path: Path to the CRISPResso_quantification_of_editing_frequency.txt file
    Returns:
        dict: modified/unmodified percentages and insertion/deletion/substitution read counts
    Raises:
        FileNotFoundError: if the file does not exist
        ValueError: if the file does not contain exactly one data row
    """
    with open(path, encoding="utf-8") as f:
        lines = []
        for line in f:
            stripped = line.strip()
            if stripped:
                lines.append(stripped)

    if len(lines[1:]) != 1:
        raise ValueError(
            f"Expected an editing frequency table with 1 data row, found {len(lines[1:])} in {path}"
        )

    headers = lines[0].split("\t")
    values = lines[1].split("\t")
    row = dict(zip(headers, values))

    return {
        "modified_pct": float(row["Modified%"]),
        "unmodified_pct": float(row["Unmodified%"]),
        "insertions": int(row["Insertions"]),
        "deletions": int(row["Deletions"]),
        "substitutions": int(row["Substitutions"]),
    }