# scripts/generate_prism_csv.py
import pandas as pd

def generate_prism_csv(df):

    # Drop empty rows
    df = df.dropna(how="all")

    if "directory" not in df.columns:
        raise ValueError("Missing 'directory' column — cannot extract base_sample or rep info.")

    expected_cols = [
        "directory", "reads_aligned", "reads_total", "correction_without_bystanders",
        "w_bystanders_less_wo_bystanders", "indep_less_w_bystanders"
    ]
    available_cols = [col for col in expected_cols if col in df.columns]
    reduced_df = df[available_cols].copy()

    # Skip empty directory rows
    reduced_df = reduced_df[reduced_df["directory"].notna() & (reduced_df["directory"] != "")]

    # Extract replicate (matches _1_L001 or -1_L001)
    reduced_df["rep"] = reduced_df["directory"].str.extract(r'[-_]([0-9]+)_L\d{3}')

    # Extract base sample name (remove the replicate + lane section)
    reduced_df["base_sample"] = reduced_df["directory"].str.replace(r'[-_]([0-9]+)_L\d{3}.*', '', regex=True)

    # Fill missing rep with "1" so single-replicate samples don't get dropped
    reduced_df["rep"] = reduced_df["rep"].fillna("1")

    # Pivot to wide format
    prism_df = reduced_df.pivot(index="base_sample", columns="rep")

    # Flatten MultiIndex
    prism_df.columns = [f"{col[0]}_rep{col[1]}" for col in prism_df.columns]
    prism_df = prism_df.reset_index()

    # Define desired output columns
    cols_to_keep = [
        "base_sample",
        "directory_rep1", "directory_rep2", "directory_rep3",
        "reads_aligned_rep1", "reads_aligned_rep2", "reads_aligned_rep3",
        "reads_total_rep1", "reads_total_rep2", "reads_total_rep3",
        "correction_without_bystanders_rep1", "correction_without_bystanders_rep2", "correction_without_bystanders_rep3",
        "w_bystanders_less_wo_bystanders_rep1", "w_bystanders_less_wo_bystanders_rep2", "w_bystanders_less_wo_bystanders_rep3",
        "indep_less_w_bystanders_rep1", "indep_less_w_bystanders_rep2", "indep_less_w_bystanders_rep3",
    ]

    # Add missing columns as NaN
    for col in cols_to_keep:
        if col not in prism_df.columns:
            prism_df[col] = pd.NA

    # Reorder
    prism_df = prism_df[cols_to_keep]


    return prism_df
