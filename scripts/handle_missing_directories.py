# scripts/handle_missing_directories.py

import os
import pandas as pd

def add_unanalyzed_directories(df: pd.DataFrame, skip_dirs=None, note="Directory not analyzed") -> pd.DataFrame:
    """
    Append directories that exist on disk but were not analyzed in the DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing analyzed directories in the "directory" column.
    skip_dirs : list, optional
        List of directories to ignore (e.g., ["scripts", "unprocessed_data"]).
    note : str, optional
        Note to write in the 'corrected_locus_with_bystanders' column for unanalyzed directories.

    Returns
    -------
    pd.DataFrame
        Updated DataFrame with placeholders for unanalyzed directories.
    """
    if skip_dirs is None:
        skip_dirs = ["scripts", "unprocessed_data"]

    # Get all directories in the current folder, excluding skip_dirs
    all_dirs = [d for d in os.listdir() if os.path.isdir(d) and d not in skip_dirs]

    # Directories already in the DataFrame
    analyzed_dirs = df["directory"].tolist()

    # Find missing directories
    missing_dirs = [d for d in all_dirs if d not in analyzed_dirs]

    # Create placeholder rows
    placeholder_rows = []
    for d in missing_dirs:
        placeholder_rows.append({
            "directory": d,
            "reads_aligned": "NA",
            "reads_total": "NA",
            "correction_with_bystanders": "NA",
            "correction_without_bystanders": "NA",
            "independent_correction": "NA",
            "indep_less_w_bystanders": "NA",
            "w_bystanders_less_wo_bystanders": "NA",
            "target_locus": "",
            "perfect_correction": "",
            "corrected_locus_with_bystanders": note
        })

    # Append placeholders
    if placeholder_rows:
        df = pd.concat([df, pd.DataFrame(placeholder_rows)], ignore_index=True)

    # Sort by directory
    df = df.sort_values(by="directory").reset_index(drop=True)

    return df