# scripts/read_extraction.py
import os
import glob
import csv


def read_extraction(directory_path):
    crispr_dirs = [d for d in glob.glob(os.path.join(directory_path, "CRISPResso_on_*")) if os.path.isdir(d)]

    if not crispr_dirs:
        print(f"No CRISPResso_on* subdirectory found in {directory_path}")
        return "NA","NA"
    
    crispresso_subfolder = crispr_dirs[0]  # first matching CRISPResso subdirectory
    statistics_file_path = glob.glob(os.path.join(crispresso_subfolder, "CRISPResso_mapping_statistics.txt"))
    if not statistics_file_path:
        print(f"No stats file found in {crispresso_subfolder}")
        return "NA","NA"
    
    with open(statistics_file_path[0], newline='') as f:
        reader = csv.reader(f, delimiter="\t")
        table = list(reader)

    if len(table) < 2:
        return "NA","NA"

    data = table[1]

    reads_aligned = data[2]
    reads_total = data[1]


    return reads_aligned, reads_total