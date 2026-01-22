import os
import glob
import csv
import logging
from .reverse_complement import reverse_complement
from scripts.logging_setup import setup_logging




def total_A_to_G(orientation, intended_edit, guide_seq, directory_path, tolerated_edits):
    #log_file = setup_logging()
    
    if orientation == "F":
        target_base = "G"
        original_base = "A"
    elif orientation == "R":
        intended_edit = len(guide_seq) - intended_edit + 1
        target_base = "C"
        original_base = "T"
        guide_seq = reverse_complement(guide_seq)
    
    crispr_dirs = []
    for d in glob.glob(os.path.join(directory_path, "CRISPResso_on_*")):
        if os.path.isdir(d):
            crispr_dirs.append(d)
    
    #if not crispr_dirs:
        #logging.error(f"no CRISPResso directories found")

    crispresso_subfolder = crispr_dirs[0]
    alleles_freqency_table_path = glob.glob(os.path.join(crispresso_subfolder, "Alleles_frequency_table_around_*"))
    #if not alleles_freqency_table_path:
        #logging.error(f"could not find allele frequency table")
    
    f = open(alleles_freqency_table_path[0], newline='')
    reader = csv.reader(f, delimiter="\t")
    table = []
    for row in reader:
        table.append(row)
    f.close()

    total_reads_any = 0
    total_reads_correct_index = 0
    total_reads_A_to_G = 0
    A_to_G_sequences = []

    for row in table[1:]:
        cur_sequence = row[0]
        reads = int(row[6])
        total_reads_any += reads

        if cur_sequence[intended_edit - 1] == target_base:
            total_reads_correct_index += reads
            
            for i in range(len(guide_seq)):
                if i == intended_edit - 1:
                    continue
                if guide_seq[i] == original_base and cur_sequence[i] == target_base:
                    total_reads_A_to_G += reads
                    A_to_G_sequences.append(cur_sequence)
                    break


    correction_with_any_change_in_protospacer = (total_reads_correct_index/total_reads_any) * 100
    total_A_to_G = (total_reads_A_to_G/total_reads_any) * 100

    return (correction_with_any_change_in_protospacer, total_A_to_G)

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
