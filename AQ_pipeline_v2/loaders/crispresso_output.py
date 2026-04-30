import pandas as pd
from pathlib import Path


# Reads CRISPResso ouput files: allele frequency tables, mapping statistics

def read_mapping_stats(path: Path) -> tuple[int, int]:
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
    df = pd.read_csv(path, sep="\t")
    return df