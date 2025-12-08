# scripts/identify_independent_correction.py
import os
import glob
import csv


def identify_independent_correction(orientation, intended_edit, directory_path):

    if orientation == "F":
        edit = "G"
    elif orientation == "R":
        edit = "C"

    crispr_dirs = [d for d in glob.glob(os.path.join(directory_path, "CRISPResso_on_*")) if os.path.isdir(d)]

    if not crispr_dirs:
        print(f"No CRISPResso_on* subdirectory found in {directory_path}")
        return "NA"
    
    crispresso_subfolder = crispr_dirs[0]  # first matching CRISPResso subdirectory
    quantification_window_file_path = glob.glob(os.path.join(crispresso_subfolder, "Quantification_window_nucleotide_percentage_table.txt"))
    if not quantification_window_file_path:
        print(f"No quantification_window_file file found in {crispresso_subfolder}")
        return "NA"
    
    with open(quantification_window_file_path[0], newline='') as f:
        reader = csv.reader(f, delimiter="\t")
        table = list(reader)
    
    if orientation == "F":
        col_index = intended_edit  
    elif orientation == "R":
        col_index = 20 - intended_edit + 1 

    for row in table:
        if row[0] == edit:
            try:
                correction_percentage = float(row[col_index])*100
                print(f"the non-read-based correction percentage is {correction_percentage}")
                return float(row[col_index])*100
            except IndexError:
                print(f"Position {intended_edit} out of range, this likely means the CRISPResso function call had to be run without a quantification window.")
                print(f"this happens when the quantification window is too close to the end of the read/amplicon")
                return "NA"

    print(f"No {edit} row found in table")
    return "NA"