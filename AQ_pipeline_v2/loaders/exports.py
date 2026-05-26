import pandas as pd
import logging

#Prism input generator function

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
    dfcopy["rep"] = df["sample"].str.extract(r'[-_](\d+)$')
    if dfcopy["rep"].isna().any():
        missing = dfcopy.loc[dfcopy["rep"].isna(), "sample"].tolist()
        logging.warning(f"Could not extract rep number from samples: {missing}")
    dfcopy["rep"] = dfcopy["rep"].fillna("1")
    
    unique_reps = dfcopy["rep"].unique()
    if len(unique_reps) > 3:
        raise ValueError(f"Prism export currently supports up to 3 replicates, got {len(unique_reps)}")


    dfcopy["base_sample"] = dfcopy["sample"].str.replace(r'[-_]\d+$', '', regex=True)

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



