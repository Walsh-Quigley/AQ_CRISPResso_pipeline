# scripts/handle_missing_directories.py

import re
import os
import pandas as pd

def add_unanalyzed_directories(df: pd.DataFrame, skip_dirs, convert_to_sample_name=False, note="Directory not analyzed") -> pd.DataFrame:
    """
    add rows for directories that weren't analyzed

    Args:
        df: DataFramge with analyzed results
        skip_dirs: List of directory names to skip
        convert_to_sample_name: Whether to convert directory names to sample names by removing suffixes
        note: Note to add for unanalyzed directories

    Returns:
        pd.DataFrame
            Updated DataFrame with placeholders for unanalyzed directories.
    """

    if skip_dirs is None:
        skip_dirs = ["scripts", "unprocessed_data"]

    # Get all directories in the current folder, excluding skip_dirs
    all_dirs = [d for d in os.listdir() if os.path.isdir(d) and d not in skip_dirs]

    # Get analyzed sample names from dataframe
    analyzed_samples = set(df['sample'].dropna())

    cleaned_all_dirs = [re.sub(r'(_L\d{3})?-ds\..*', '', d) for d in all_dirs]

    # Find unanalyzed directories
    unanalyzed_dirs = []
    for directory in all_dirs:
        # Convert to sample name using same logic as in processing
        if convert_to_sample_name:
            sample_name = re.sub(r'(_L\d{3})?-ds\..*', '', directory)
        else:
            sample_name = directory
        
        if sample_name not in analyzed_samples:
            unanalyzed_dirs.append(sample_name)

    # Directories already in the DataFrame
    analyzed_dirs = df["sample"].tolist()

    # Find missing directories
    missing_dirs = [
        d for d, cleaned in zip(all_dirs, cleaned_all_dirs)
        if cleaned not in analyzed_dirs
    ]

    # Create rows for unanalyzed directories
    if unanalyzed_dirs:
        unanalyzed_rows = []
        for sample in unanalyzed_dirs:
            new_row = {col: None for col in df.columns}
            new_row['sample'] = sample
            new_row['note'] = note
            unanalyzed_rows.append(new_row)
        
        # Append to dataframe
        df_unanalyzed = pd.DataFrame(unanalyzed_rows)
        df = pd.concat([df, df_unanalyzed], ignore_index=True)
    
    return df