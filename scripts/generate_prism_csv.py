# scripts/generate_prism_csv.py
import re
import pandas as pd

def generate_prism_csv(df):
    
    df = df.dropna(how="all")

    if "directory" not in df.columns:
        raise ValueError("DataFrame must contain a 'directory' column.")
    
    df["rep"] = df["directory"].str.extract(r'[-_](\d+)(?=(?:_L\d{3}|-ds|$))')
    df["rep"] = df["rep"].fillna("1")

    df["base_sample"] = df["directory"].str.replace(r'_(\d+)(?:-ds.*)?$', '', regex=True)

    numeric_cols = [
        "reads_aligned",
        "reads_total",
        "correction_without_bystanders",
        "w_bystanders_less_wo_bystanders",
        "indep_less_w_bystanders",
    ]
    
    numeric_df = df.pivot(index="base_sample", columns="rep", values=numeric_cols)
    numeric_df.columns = [f"{col[0]}_rep{col[1]}" for col in numeric_df.columns]
    numeric_df = numeric_df.reset_index()

    dir_df = df.pivot(index="base_sample", columns="rep", values="directory")
    dir_df.columns = [f"directory_rep{col}" for col in dir_df.columns]
    dir_df = dir_df.reset_index()

    prism_df = pd.merge(numeric_df, dir_df, on="base_sample", how="left")

    cols_to_keep = [
        "base_sample",
        "directory_rep1", "directory_rep2", "directory_rep3",
        "reads_aligned_rep1", "reads_aligned_rep2", "reads_aligned_rep3",
        "reads_total_rep1", "reads_total_rep2", "reads_total_rep3",
        "correction_without_bystanders_rep1", "correction_without_bystanders_rep2", "correction_without_bystanders_rep3",
        "w_bystanders_less_wo_bystanders_rep1", "w_bystanders_less_wo_bystanders_rep2", "w_bystanders_less_wo_bystanders_rep3",
        "indep_less_w_bystanders_rep1", "indep_less_w_bystanders_rep2", "indep_less_w_bystanders_rep3",
    ]

    for col in cols_to_keep:
        if col not in prism_df.columns:
            prism_df[col] = pd.NA

    prism_df = prism_df[cols_to_keep]

    return prism_df
