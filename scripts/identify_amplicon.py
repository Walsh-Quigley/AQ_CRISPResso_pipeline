# scripts/identify_amplicon.py
import logging

def identify_amplicon(directory, amplicon_names):
    """
    Identify which amplicon a directory corresponds to based on its name.
    
    Args:
        directory: Directory name to match
        amplicon_names: List of valid amplicon names
        
    Returns:
        str: Matched amplicon name
        
    Raises:
        ValueError: If no matching amplicon is found
    """

    matched_name = None
    directory_upper = directory.upper()

    # Sort names by length, longest first to prevent substring collisions
    for name in sorted(amplicon_names, key=len, reverse=True):
        if name.upper() in directory_upper:
            matched_name = name
            break

    if not matched_name:
        error_msg = f"No valid amplicon match found for directory: {directory}"
        logging.error(error_msg)
        raise ValueError(error_msg)
    
    return matched_name