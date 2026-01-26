#scripts/process_ABE_case.py
"""
import os
import re
import logging
from .generate_search_sequences import generate_search_sequences
from .read_extraction import read_extraction
from .filter_alleles_file import filter_alleles_file
from .filter_alleles_file import filter_alleles_file_hetero
#from .identify_independent_correction import identify_independent_correction
from .identify_independent_correction import total_A_to_G
from .identify_independent_correction import total_A_to_G_hetero


def process_ABE_case(directory_path, guide_seq, orientation, editor, intended_edit, tolerated_edits):
    
    Process non-ONE-seq target amplicon sequencing samples
    
    Args:
        directory_path: Full path to the directory to process
        guide_seq: the off target locus
        orientation: Orientation of the guide relative to amplicon
        editor: The base editor used
        intended_edit: The position of the intended edit
        tolerated_edits: List of positions of tolerated bystander edits
        
    Returns:
        dict: Result dictionary with metrics from CRISPResso output
    

    #generate search sequences
    search_strings = generate_search_sequences(guide_seq, orientation, editor, intended_edit, tolerated_edits, directory_path)

    #filter CRISPResso_ouput*/alleles_frequency_table file for search sequences
    correction_with_bystander, correction_without_bystanders = filter_alleles_file(search_strings, directory_path)
    pct_w_base1, pct_w_base2, pct_wo_base1, pct_wo_base2, het_pos, base1, base2 = filter_alleles_file_hetero(search_strings, directory_path, orientation, guide_seq)


    #identify independent correction from CRISPResso_output*/Quantification_window_nucleotide_percentage_table.txt
    #independent_correction = identify_independent_correction(orientation, intended_edit, directory_path)
    correction_with_any_change_in_protospacer, total_A_to_G_value = total_A_to_G(orientation, intended_edit, guide_seq, directory_path, tolerated_edits)    
    correction_protospacer_allele1, correction_protospacer_allele2, total_A_to_G_allele1, total_A_to_G_allele2 = total_A_to_G_hetero(orientation, intended_edit, guide_seq, directory_path, tolerated_edits)

    #Extract sample name
    directory_name = os.path.basename(directory_path.rstrip('/'))
    sample_name = re.sub(r'(_L\d{3})?-ds\..*', '', directory_name)
    
    if correction_with_bystander == "NA" or correction_without_bystanders == "NA" or correction_with_any_change_in_protospacer == "NA" or total_A_to_G_value == "NA":
        error_msg = f"Missing correction data for {sample_name}: " \
                    f"with_bystanders={correction_with_bystander}, " \
                    f"without_bystanders={correction_without_bystanders},"\
                    f"independent_correction={correction_with_any_change_in_protospacer}"
        logging.error(error_msg)
        raise ValueError(error_msg)
    

    #log the outputs for the user
    logging.info(f"The directory {directory_path} has the following metrics")
    logging.info(f"correction_with_bystander: {correction_with_bystander}")
    logging.info(f"correction_without_bystanders: {correction_without_bystanders}")
    logging.info(f"frequency of correction independent of read: {correction_with_any_change_in_protospacer}")

    #get the read counts:
    reads_aligned, reads_total = read_extraction(directory_path)
    indep_less_w_bystanders = correction_with_any_change_in_protospacer - correction_with_bystander
    w_bystanders_less_wo_bystanders = correction_with_bystander - correction_without_bystanders

    # Allele-specific difference calculations
    indep_less_w_bystanders_allele1 = correction_protospacer_allele1 - pct_w_base1 if (correction_protospacer_allele1 != "NA" and pct_w_base1 != "NA") else "NA"
    indep_less_w_bystanders_allele2 = correction_protospacer_allele2 - pct_w_base2 if (correction_protospacer_allele2 != "NA" and pct_w_base2 != "NA") else "NA"
    w_bystanders_less_wo_bystanders_allele1 = pct_w_base1 - pct_wo_base1 if (pct_w_base1 != "NA" and pct_wo_base1 != "NA") else "NA"
    w_bystanders_less_wo_bystanders_allele2 = pct_w_base2 - pct_wo_base2 if (pct_w_base2 != "NA" and pct_wo_base2 != "NA") else "NA"
        
    # Warn about unexpected negative values but keep the numbers
    if indep_less_w_bystanders < 0:
        logging.warning(f"Negative independent correction rate for {sample_name}: {indep_less_w_bystanders:.2f}")
    
    if w_bystanders_less_wo_bystanders < 0:
        logging.warning(f"Negative read-based correction rate for {sample_name}: {w_bystanders_less_wo_bystanders:.2f}")

    result = {
        "sample":sample_name,
        "reads_aligned": reads_aligned,
        "reads_total": reads_total,

        "correction_without_bystanders":correction_without_bystanders,
        "correction_with_bystanders":correction_with_bystander,
        "correction_with_any_AtoG_change": total_A_to_G_value,
        "correction_with_any_change_in_protospacer": correction_with_any_change_in_protospacer,
        "column F minus column E": w_bystanders_less_wo_bystanders,
        "column": indep_less_w_bystanders,
        
        "correction_w_bystanders_allele1": pct_w_base1,
        "correction_wo_bystanders_allele1": pct_wo_base1,
        "correction_w_bystanders_allele2": pct_w_base2,
        "correction_wo_bystanders_allele2": pct_wo_base2,
        "correction_protospacer_allele1": correction_protospacer_allele1,
        "total_A_to_G_allele1": total_A_to_G_allele1,
        "correction_protospacer_allele2": correction_protospacer_allele2,
        "total_A_to_G_allele2": total_A_to_G_allele2,
        "indep_less_w_bystanders_allele1": indep_less_w_bystanders_allele1,
        "indep_less_w_bystanders_allele2": indep_less_w_bystanders_allele2,
        "w_bystanders_less_wo_bystanders_allele1": w_bystanders_less_wo_bystanders_allele1,
        "w_bystanders_less_wo_bystanders_allele2": w_bystanders_less_wo_bystanders_allele2,
        "het_position": het_pos,
        "het_alleles": f"{base1}/{base2}" if het_pos != "NA" else "NA",
        "target_locus":guide_seq,
        "perfect_correction":search_strings[0],
        "corrected_locus_with_bystanders": ";".join(search_strings)
        }
    
    return result
"""