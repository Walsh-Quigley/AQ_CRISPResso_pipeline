import os
import glob
import csv
import logging
from .reverse_complement import reverse_complement
from scripts.logging_setup import setup_logging
from .filter_alleles_file import find_het_position



def total_A_to_G_hetero(orientation, intended_edit, guide_seq, directory_path, tolerated_edits):
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

    het_pos, base1, base2 = find_het_position(crispresso_subfolder)
    if het_pos is None:
        print(f"No heterozygous position found in {crispresso_subfolder}")
        return "NA", "NA", "NA", "NA"
    
    alleles_freqency_table_path = glob.glob(os.path.join(crispresso_subfolder, "Alleles_frequency_table_around_*"))
    #if not alleles_freqency_table_path:
        #logging.error(f"could not find allele frequency table")
    
    f = open(alleles_freqency_table_path[0], newline='')
    reader = csv.reader(f, delimiter="\t")
    table = []
    for row in reader:
        table.append(row)
    f.close()

    #allele 1:
    total_reads_base1 = 0
    correct_index_base1 = 0
    A_to_G_base1 = 0

    #allele 2:
    total_reads_base2 = 0
    correct_index_base2 = 0
    A_to_G_base2 = 0

    for row in table[1:]:
        cur_sequence = row[0]
        reads = int(row[6])
        
        if cur_sequence[het_pos] == base1:
            total_reads_base1 += reads
            if cur_sequence[intended_edit - 1] == target_base:
                correct_index_base1 += reads
                has_non_AtoG_change = False
                for i in range(len(guide_seq)):
                    if cur_sequence[i] != guide_seq[i]:
                        if not (guide_seq[i] == original_base and cur_sequence[i] == target_base):
                            has_non_AtoG_change = True
                            break
                if not has_non_AtoG_change:
                    A_to_G_base1 += reads


        elif cur_sequence[het_pos] == base2:
            total_reads_base2 += reads
            if cur_sequence[intended_edit - 1] == target_base:
                correct_index_base2 += reads
                has_non_AtoG_change = False
                for i in range(len(guide_seq)):
                    if cur_sequence[i] != guide_seq[i]:
                        if not (guide_seq[i] == original_base and cur_sequence[i] == target_base):
                            has_non_AtoG_change = True
                            break
                if not has_non_AtoG_change:
                    A_to_G_base2 += reads

    correction_with_any_change_in_protospacer_allele1 = ((correct_index_base1/total_reads_base1) * 100) if total_reads_base1 > 0 else 0
    correction_with_any_change_in_protospacer_allele2 = ((correct_index_base2/total_reads_base2) * 100) if total_reads_base2 > 0 else 0
   
    pct_A_to_G_base1 = ((A_to_G_base1/total_reads_base1) * 100) if total_reads_base1 > 0 else 0
    pct_A_to_G_base2 = ((A_to_G_base2/total_reads_base2) * 100) if total_reads_base2 > 0 else 0


    return (correction_with_any_change_in_protospacer_allele1, 
            correction_with_any_change_in_protospacer_allele2,
            pct_A_to_G_base1,
            pct_A_to_G_base2)


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
    #A_to_G_sequences = []

    for row in table[1:]:
        cur_sequence = row[0]
        reads = int(row[6])
        total_reads_any += reads

        if cur_sequence[intended_edit - 1] == target_base:
            total_reads_correct_index += reads
            
            has_non_AtoG_change = False
            for i in range(len(guide_seq)):
                if cur_sequence[i] != guide_seq[i]:
                    if not (guide_seq[i] == original_base and cur_sequence[i] == target_base):
                        has_non_AtoG_change = True
                        break

            if not has_non_AtoG_change:
                total_reads_A_to_G += reads


    correction_with_any_change_in_protospacer = (total_reads_correct_index/total_reads_any) * 100
    pct_A_to_G = (total_reads_A_to_G/total_reads_any) * 100

    return (correction_with_any_change_in_protospacer, pct_A_to_G)

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
