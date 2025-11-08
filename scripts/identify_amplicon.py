# scripts/identify_amplicon.py

def identify_amplicon(directory, amplicon_names):
    matched_name = None
    directory_upper = directory.upper()

    # Sort names by length, longest first to prevent substring collisions
    for name in sorted(amplicon_names, key=len, reverse=True):
        if name.upper() in directory_upper:
            matched_name = name
            break

    if not matched_name:
        print(f"No valid match found in amplicon list; {directory}")
 
    return matched_name