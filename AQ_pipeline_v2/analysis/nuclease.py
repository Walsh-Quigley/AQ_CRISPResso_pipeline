import pandas as pd

def calculate_frameshift(allele_table_df: pd.DataFrame) -> dict:
    """% of aligned reads that are frameshifted (net indel not divisible by 3)
    vs in-frame indels, from the allele frequency table.
    net indel per allele = n_inserted - n_deleted; frameshift if net % 3 != 0.
    Args:
        allele_table_df: dataframe containing read data from the sample's allele frequency table
    Returns:
        dict: a dictionary matching percentage of in frame and frameshift indels to their values"""
    if allele_table_df.empty:
        return {"pct_frameshift_indels": 0.0, "pct_inframe_indels": 0.0}

    net = allele_table_df["n_inserted"] - allele_table_df["n_deleted"]
    has_indel = (allele_table_df["n_inserted"] != 0) | (allele_table_df["n_deleted"] != 0)
    is_frameshift = (net % 3) != 0                  # net%3!=0 implies an indel exists
    is_inframe_indel = has_indel & ~is_frameshift   # in-frame but still has an indel (e.g. 3bp del)

    return {
        "pct_frameshift_indels": round(allele_table_df.loc[is_frameshift, "%Reads"].sum(), 2),
        "pct_inframe_indels": round(allele_table_df.loc[is_inframe_indel, "%Reads"].sum(), 2),
    }
