# scripts/generate_prism_csv.py
import re
import logging
import pandas as pd

def generate_prism_csv(df):
    
    df = df.dropna(how="all")

    if "sample" not in df.columns:
        raise ValueError("DataFrame must contain a 'sample' column.")
    
    # Separate out samples with notes (ONE-seq or unanalyzed)
    if "note" in df.columns:
        noted_samples = df[df['note'].notna()].copy()
        # Only process samples without notes for Prism analysis
        df_to_process = df[df['note'].isna()].copy()
    else:
        noted_samples = pd.DataFrame()  # Empty dataframe
        df_to_process = df.copy()

    if df_to_process.empty:
        logging.warning("No samples available for Prism analysis after filtering")
        return pd.DataFrame(columns=["base_sample", "status"])
    
    # Extract replicate number
    df_to_process["rep"] = df_to_process["sample"].str.extract(r'[-_](\d+)$')
    df_to_process["rep"] = df_to_process["rep"].fillna("1")

    # Extract base sample name (remove the trailing _#)
    df_to_process["base_sample"] = df_to_process["sample"].str.replace(r'[-_](\d+)$', '', regex=True)

    # Debug logging
    logging.info("Sample name extraction:")
    for idx, row in df_to_process[['sample', 'base_sample', 'rep']].head().iterrows():
        logging.info(f"  Sample: {row['sample']} → Base: {row['base_sample']}, Rep: {row['rep']}")

    #Numeric columns to pivot
    numeric_cols = [
        "reads_aligned",
        "reads_total",
        "correction_without_bystanders",
        "w_bystanders_less_wo_bystanders",
        "indep_less_w_bystanders",
    ]
    
    #Pivot numeric data
    numeric_df = df_to_process.pivot(index="base_sample", columns="rep", values=numeric_cols)
    numeric_df.columns = [f"{col[0]}_rep{col[1]}" for col in numeric_df.columns]
    numeric_df = numeric_df.reset_index()

    #Pivot directory/sample names
    dir_df = df_to_process.pivot(index="base_sample", columns="rep", values="sample")
    dir_df.columns = [f"sample_rep{col}" for col in dir_df.columns]
    dir_df = dir_df.reset_index()

    prism_df = pd.merge(numeric_df, dir_df, on="base_sample", how="left")

    #define expected columns
    cols_to_keep = [
        "base_sample",
        "sample_rep1", "sample_rep2", "sample_rep3",
        "reads_aligned_rep1", "reads_aligned_rep2", "reads_aligned_rep3",
        "reads_total_rep1", "reads_total_rep2", "reads_total_rep3",
        "correction_without_bystanders_rep1", "correction_without_bystanders_rep2", "correction_without_bystanders_rep3",
        "w_bystanders_less_wo_bystanders_rep1", "w_bystanders_less_wo_bystanders_rep2", "w_bystanders_less_wo_bystanders_rep3",
        "indep_less_w_bystanders_rep1", "indep_less_w_bystanders_rep2", "indep_less_w_bystanders_rep3",
    ]

    #add missing columns
    for col in cols_to_keep:
        if col not in prism_df.columns:
            prism_df[col] = pd.NA

    prism_df = prism_df[cols_to_keep]

    # Now add the noted samples at the bottom
    if not noted_samples.empty:
        logging.info(f"Excluding {len(noted_samples)} samples from Prism analysis (ONE-seq or unanalyzed)")
        
        # Create rows for noted samples
        excluded_rows = []
        for _, row in noted_samples.iterrows():
            excluded_row = {col: pd.NA for col in cols_to_keep}
            excluded_row["base_sample"] = row["sample"]
            excluded_row["sample_rep1"] = f"Not included in Prism analysis - {row['note']}"
            excluded_rows.append(excluded_row)
        
        # Append to prism_df
        excluded_df = pd.DataFrame(excluded_rows)
        prism_df = pd.concat([prism_df, excluded_df], ignore_index=True)

    return prism_df
