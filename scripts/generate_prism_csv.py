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

    # Map old column names to prism output names
    column_mapping = {
        "column E minus column D": "w_bystanders_minus_wo_bystanders",
        "column F minus column E": "any_AtoG_minus_w_bystanders",
        "column G minus column F": "any_protospacer_change_minus_any_AtoG",
    }


    # Rename columns if old names are present
    df_to_process = df_to_process.rename(columns=column_mapping)

    #Numeric columns to pivot
    numeric_cols = [
        "reads_aligned",
        "reads_total",
        "correction_without_bystanders",
        "w_bystanders_minus_wo_bystanders",
        "any_AtoG_minus_w_bystanders",
        "any_protospacer_change_minus_any_AtoG",
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
        "w_bystanders_minus_wo_bystanders_rep1", "w_bystanders_minus_wo_bystanders_rep2", "w_bystanders_minus_wo_bystanders_rep3",
        "any_AtoG_minus_w_bystanders_rep1", "any_AtoG_minus_w_bystanders_rep2", "any_AtoG_minus_w_bystanders_rep3",
        "any_protospacer_change_minus_any_AtoG_rep1", "any_protospacer_change_minus_any_AtoG_rep2", "any_protospacer_change_minus_any_AtoG_rep3",
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

def generate_prism_csv_het(df):
    """Generate Prism-formatted CSVs for heterozygous samples.
    
    Returns a tuple of two DataFrames (allele1_df, allele2_df), each with the 
    same structure as the non-het prism output.
    """
    
    df = df.dropna(how="all")

    if "sample" not in df.columns:
        raise ValueError("DataFrame must contain a 'sample' column.")
    
    # Separate out samples with notes
    if "note" in df.columns:
        noted_samples = df[df['note'].notna()].copy()
        df_to_process = df[df['note'].isna()].copy()
    else:
        noted_samples = pd.DataFrame()
        df_to_process = df.copy()

    if df_to_process.empty:
        logging.warning("No het samples available for Prism analysis after filtering")
        empty_df = pd.DataFrame(columns=["base_sample", "status"])
        return (empty_df, empty_df)
    
    # Extract replicate number
    df_to_process["rep"] = df_to_process["sample"].str.extract(r'[-_](\d+)$')
    df_to_process["rep"] = df_to_process["rep"].fillna("1")

    # Extract base sample name
    df_to_process["base_sample"] = df_to_process["sample"].str.replace(r'[-_](\d+)$', '', regex=True)

    # Map old column names for allele1
    column_mapping_allele1 = {
        "reads_aligned_allele1": "reads_aligned",
        "correction_without_bystanders_allele1": "correction_without_bystanders",
        "column L minus column K": "w_bystanders_minus_wo_bystanders",
        "column M minus column L": "any_AtoG_minus_w_bystanders",
        "column N minus column M": "any_protospacer_change_minus_any_AtoG",
    }

    # Map old column names for allele2
    column_mapping_allele2 = {
        "reads_aligned_allele2": "reads_aligned",
        "correction_without_bystanders_allele2": "correction_without_bystanders",
        "column S minus column R": "w_bystanders_minus_wo_bystanders",
        "column T minus column S": "any_AtoG_minus_w_bystanders",
        "column U minus column T": "any_protospacer_change_minus_any_AtoG",
    }

    # Numeric columns to pivot (no reads_total for het)
    numeric_cols = [
        "reads_aligned",
        "correction_without_bystanders",
        "w_bystanders_minus_wo_bystanders",
        "any_AtoG_minus_w_bystanders",
        "any_protospacer_change_minus_any_AtoG",
    ]

    # Define expected columns (no reads_total for het)
    cols_to_keep = [
        "base_sample",
        "sample_rep1", "sample_rep2", "sample_rep3",
        "reads_aligned_rep1", "reads_aligned_rep2", "reads_aligned_rep3",
        "correction_without_bystanders_rep1", "correction_without_bystanders_rep2", "correction_without_bystanders_rep3",
        "w_bystanders_minus_wo_bystanders_rep1", "w_bystanders_minus_wo_bystanders_rep2", "w_bystanders_minus_wo_bystanders_rep3",
        "any_AtoG_minus_w_bystanders_rep1", "any_AtoG_minus_w_bystanders_rep2", "any_AtoG_minus_w_bystanders_rep3",
        "any_protospacer_change_minus_any_AtoG_rep1", "any_protospacer_change_minus_any_AtoG_rep2", "any_protospacer_change_minus_any_AtoG_rep3",
        "het_position", "het_alleles",
    ]

    # Drop non-allele-specific columns that would conflict with renamed allele columns
    cols_to_drop = [
        "reads_aligned",
        "reads_total",
        "correction_without_bystanders",
        "correction_with_bystanders",
        "correction_with_any_AtoG_change",
        "correction_with_any_change_in_protospacer",
        "column E minus column D",
        "column F minus column E",
        "column G minus column F",
    ]
    df_to_process = df_to_process.drop(columns=[c for c in cols_to_drop if c in df_to_process.columns])

    allele1_df = _process_allele_for_prism(df_to_process, column_mapping_allele1, numeric_cols, cols_to_keep, noted_samples)
    allele2_df = _process_allele_for_prism(df_to_process, column_mapping_allele2, numeric_cols, cols_to_keep, noted_samples)

    return (allele1_df, allele2_df)


def _process_allele_for_prism(df_input, column_mapping, numeric_cols, cols_to_keep, noted_samples):
    """Helper function to process a single allele for Prism output."""
    df_allele = df_input.copy()
    df_allele = df_allele.rename(columns=column_mapping)
    
    # Pivot numeric data
    numeric_df = df_allele.pivot(index="base_sample", columns="rep", values=numeric_cols)
    numeric_df.columns = [f"{col[0]}_rep{col[1]}" for col in numeric_df.columns]
    numeric_df = numeric_df.reset_index()

    # Pivot directory/sample names
    dir_df = df_allele.pivot(index="base_sample", columns="rep", values="sample")
    dir_df.columns = [f"sample_rep{col}" for col in dir_df.columns]
    dir_df = dir_df.reset_index()

    # Get het info (take first value since should be same across reps)
    het_info_df = df_allele.groupby("base_sample")[["het_position", "het_alleles"]].first().reset_index()

    prism_df = pd.merge(numeric_df, dir_df, on="base_sample", how="left")
    prism_df = pd.merge(prism_df, het_info_df, on="base_sample", how="left")

    # Add missing columns
    for col in cols_to_keep:
        if col not in prism_df.columns:
            prism_df[col] = pd.NA

    prism_df = prism_df[cols_to_keep]

    # Add noted samples at the bottom
    if not noted_samples.empty:
        excluded_rows = []
        for _, row in noted_samples.iterrows():
            excluded_row = {col: pd.NA for col in cols_to_keep}
            excluded_row["base_sample"] = row["sample"]
            excluded_row["sample_rep1"] = f"Not included in Prism analysis - {row['note']}"
            excluded_rows.append(excluded_row)
        
        excluded_df = pd.DataFrame(excluded_rows)
        prism_df = pd.concat([prism_df, excluded_df], ignore_index=True)

    return prism_df
