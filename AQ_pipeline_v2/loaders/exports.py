import pandas as pd
import logging

#translates Quantification_Summary into wide format for GraphPad Prism

def generate_prism_csv(df: pd.DataFrame) -> pd.DataFrame:
    """Translates the Quantification_Summary dataframe into a format that is acceptable by
    Prism.
    Args:
        df: the dataframe that contains the Quantification_Summary's information
    Returns:
        pd.DataFrame: a newly formatted dataframe
    Raises:
        ValueError: if the sample column is missing from the input DataFrame
        ValueError: if more than 3 replicates are detected

    """
    df = df.dropna(how="all")
    if df.empty:
        logging.warning("Passed through DataFrame empty, returning empty output")
        return pd.DataFrame()
    if "sample" not in df.columns:
        error_msg = "Sample column missing in quantification summary"
        raise ValueError(error_msg)
  
    dfcopy = df.copy()
    dfcopy["rep"] = df["sample"].str.extract(r'[-_](?:rep)?(\d+)$')
    if dfcopy["rep"].isna().any():
        missing = dfcopy.loc[dfcopy["rep"].isna(), "sample"].tolist()
        logging.warning(f"Could not extract rep number from samples: {missing}")
    dfcopy["rep"] = dfcopy["rep"].fillna("1")
    
    unique_reps = dfcopy["rep"].unique()
    if len(unique_reps) > 3:
        raise ValueError(f"Prism export currently supports up to 3 replicates, got {len(unique_reps)}")


    dfcopy["base_sample"] = dfcopy["sample"].str.replace(r'[-_](?:rep)?\d+$', '', regex=True)

    numeric_cols = ["reads_total", "reads_aligned", "correction_without_bystanders", 
                    "w_bystanders_minus_wo_bystanders", "any_AtoG_minus_w_bystanders",
                    "any_change_minus_any_AtoG"]
    
    numeric_df = dfcopy.pivot(index="base_sample", columns="rep", values=numeric_cols)

    numeric_df.columns = [f"{col[0]}_rep{col[1]}" for col in numeric_df.columns]
    numeric_df = numeric_df.reset_index()

    dir_df = dfcopy.pivot(index="base_sample", columns="rep", values="sample")
    dir_df.columns = [f"sample_rep{col}" for col in dir_df.columns]
    dir_df = dir_df.reset_index()

    prism_df = pd.merge(numeric_df, dir_df, on="base_sample", how="left")

    cols_to_keep = [
        "base_sample",
        "sample_rep1", "sample_rep2", "sample_rep3",
        "reads_aligned_rep1", "reads_aligned_rep2", "reads_aligned_rep3",
        "reads_total_rep1", "reads_total_rep2", "reads_total_rep3",
        "correction_without_bystanders_rep1", "correction_without_bystanders_rep2", "correction_without_bystanders_rep3",
        "w_bystanders_minus_wo_bystanders_rep1", "w_bystanders_minus_wo_bystanders_rep2", "w_bystanders_minus_wo_bystanders_rep3",
        "any_AtoG_minus_w_bystanders_rep1", "any_AtoG_minus_w_bystanders_rep2", "any_AtoG_minus_w_bystanders_rep3",
        "any_change_minus_any_AtoG_rep1", "any_change_minus_any_AtoG_rep2", "any_change_minus_any_AtoG_rep3",
    ]

    for col in cols_to_keep:
        if col not in prism_df.columns:
            prism_df[col] = pd.NA

    return prism_df[cols_to_keep]

def generate_prism_csv_het(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Translates a heterozygous-samples-only Quantification_Summary dataframe into
    two Prism-formatted dataframes — one per allele — suitable for direct import
    into GraphPad Prism.
    
    Args:
        df: dataframe containing ONLY heterozygous samples — every row must have
            non-null values in the allele1/allele2 correction columns. Callers
            should pre-filter by `df["het_position"].notna()`.
    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: (allele1_df, allele2_df) — two wide-format
            DataFrames, one per allele, each with up to 3 replicates per base_sample.
    Raises:
        ValueError: if the sample column is missing from the input DataFrame
        ValueError: if more than 3 replicates are detected
    Note:
        reads_total is intentionally excluded from the het Prism output — only
        reads_aligned is included, since the per-allele metrics are computed
        against per-allele totals, not the sample-wide total.
    """
    df = df.dropna(how="all")
    if df.empty:
        logging.warning("Passed through DataFrame empty, returning empty output")
        return pd.DataFrame()
    if "sample" not in df.columns:
        raise ValueError("Sample column missing in quantification summary")

    dfcopy = df.copy()

    # drop global columns that would collide with the renamed allele columns
    conflicting_globals = [
        "reads_aligned",   # NEW — replaced by per-allele reads in het case
        "w_bystanders_minus_wo_bystanders",
        "any_AtoG_minus_w_bystanders",
        "any_change_minus_any_AtoG",
    ]
    dfcopy = dfcopy.drop(columns=[c for c in conflicting_globals if c in dfcopy.columns])

    # rep extraction
    dfcopy["rep"] = df["sample"].str.extract(r'[-_](?:rep)?(\d+)$')
    if dfcopy["rep"].isna().any():
        missing = dfcopy.loc[dfcopy["rep"].isna(), "sample"].tolist()
        logging.warning(f"Could not extract rep number from samples: {missing}")
    dfcopy["rep"] = dfcopy["rep"].fillna("1")

    dfcopy["base_sample"] = df["sample"].str.replace(r'[-_](?:rep)?\d+$', '', regex=True)

    # rep count guard
    unique_reps = dfcopy["rep"].unique()
    if len(unique_reps) > 3:
        raise ValueError(f"Prism export currently supports up to 3 replicates, got {len(unique_reps)}")

    allele1_df = _process_het_allele(dfcopy, "allele1")
    allele2_df = _process_het_allele(dfcopy, "allele2")
    return (allele1_df, allele2_df)

def _process_het_allele(df_input: pd.DataFrame, allele_suffix: str) -> pd.DataFrame:
    """Internal helper: process one allele for het Prism output. Called twice 
    (once per allele) by generate_prism_csv_het.
    Args:
        df_input: a het samples dataframe with rep and base_sample columns 
            already extracted by the caller
        allele_suffix: either "allele1" or "allele2" — determines which 
            allele-specific source columns are pulled
    Returns:
        pd.DataFrame: wide-format Prism output for the specified allele
    """
    
    column_mapping = {
        f"correction_wo_bystanders_{allele_suffix}": "correction_wo_bystanders",
        f"w_bystanders_minus_wo_bystanders_{allele_suffix}": "w_bystanders_minus_wo_bystanders",
        f"any_AtoG_minus_w_bystanders_{allele_suffix}": "any_AtoG_minus_w_bystanders",
        f"any_change_minus_any_AtoG_{allele_suffix}": "any_change_minus_any_AtoG",
        f"reads_aligned_{allele_suffix}": "reads_aligned",
    }
    renamed = df_input.rename(columns=column_mapping)

    numeric_cols = [
        "reads_aligned",
        "correction_wo_bystanders",
        "w_bystanders_minus_wo_bystanders",
        "any_AtoG_minus_w_bystanders",
        "any_change_minus_any_AtoG",
    ]

    numeric_df = renamed.pivot(index="base_sample", columns="rep", values=numeric_cols)
    numeric_df.columns = [f"{col[0]}_rep{col[1]}" for col in numeric_df.columns]
    numeric_df = numeric_df.reset_index()

    dir_df = renamed.pivot(index="base_sample", columns="rep", values="sample")
    dir_df.columns = [f"sample_rep{col}" for col in dir_df.columns]
    dir_df = dir_df.reset_index()

    # het_position and het_alleles are identical across reps for a given base_sample
    # group and take first so we can merge them back in
    het_info_df = renamed.groupby("base_sample")[["het_position", "het_alleles"]].first().reset_index()

    wide = pd.merge(numeric_df, dir_df, on="base_sample", how="left")
    wide = pd.merge(wide, het_info_df, on="base_sample", how="left")

    cols_to_keep = [
        "base_sample",
        "sample_rep1", "sample_rep2", "sample_rep3",
        "reads_aligned_rep1", "reads_aligned_rep2", "reads_aligned_rep3",
        "correction_wo_bystanders_rep1", "correction_wo_bystanders_rep2", "correction_wo_bystanders_rep3",
        "w_bystanders_minus_wo_bystanders_rep1", "w_bystanders_minus_wo_bystanders_rep2", "w_bystanders_minus_wo_bystanders_rep3",
        "any_AtoG_minus_w_bystanders_rep1", "any_AtoG_minus_w_bystanders_rep2", "any_AtoG_minus_w_bystanders_rep3",
        "any_change_minus_any_AtoG_rep1", "any_change_minus_any_AtoG_rep2", "any_change_minus_any_AtoG_rep3",
        "het_position",
        "het_alleles",
    ]

    for col in cols_to_keep:
        if col not in wide.columns:
            wide[col] = pd.NA

    return wide[cols_to_keep]
