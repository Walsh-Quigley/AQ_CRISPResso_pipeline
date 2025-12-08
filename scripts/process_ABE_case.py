#scripts/process_ABE_case.py

import os
import re
import logging
from .generate_search_sequences import generate_search_sequences
from .read_extraction import read_extraction
from .filter_alleles_file import filter_alleles_file
from .identify_independent_correction import identify_independent_correction


def process_ABE_case(directory_path, guide_seq, orientation, editor, intended_edit, tolerated_edits):
    """
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
    """

    #generate search sequences
    search_strings = generate_search_sequences(guide_seq, orientation, editor, intended_edit, tolerated_edits, directory_path)

    #filter CRISPResso_ouput*/alleles_frequency_table file for search sequences
    correction_with_bystander, correction_without_bystanders = filter_alleles_file(search_strings, directory_path)

    #identify independent correction from CRISPResso_output*/Quantification_window_nucleotide_percentage_table.txt
    independent_correction = identify_independent_correction(orientation, intended_edit, directory_path)

    #Extract sample name
    directory_name = os.path.basename(directory_path.rstrip('/'))
    sample_name = re.sub(r'(_L\d{3})?-ds\..*', '', directory_name)
    
    if  correction_with_bystander == "NA" or correction_without_bystanders == "NA" or independent_correction == "NA":
        error_msg = f"Missing correction data for {sample_name}: " \
                    f"with_bystanders={correction_with_bystander}, " \
                    f"without_bystanders={correction_without_bystanders},"\
                    f"independent_correction={independent_correction}"
        logging.error(error_msg)
        raise ValueError(error_msg)

    #log the outputs for the user
    logging.info(f"The directory {directory_path} has the following metrics")
    logging.info(f"correction_with_bystander: {correction_with_bystander}")
    logging.info(f"correction_without_bystanders: {correction_without_bystanders}")
    logging.info(f"frequency of correction independent of read: {independent_correction}")


    #get the read counts:
    reads_aligned, reads_total = read_extraction(directory_path)
    indep_less_w_bystanders = independent_correction - correction_with_bystander
    w_bystanders_less_wo_bystanders = correction_with_bystander - correction_without_bystanders
    
    # Warn about unexpected negative values but keep the numbers
    if indep_less_w_bystanders < 0:
        logging.warning(f"Negative independent correction rate for {sample_name}: {indep_less_w_bystanders:.2f}")
    
    if w_bystanders_less_wo_bystanders < 0:
        logging.warning(f"Negative read-based correction rate for {sample_name}: {w_bystanders_less_wo_bystanders:.2f}")

    result = {
        "sample":sample_name,
        "reads_aligned": reads_aligned,
        "reads_total": reads_total,
        "correction_with_bystanders":correction_with_bystander,
        "correction_without_bystanders":correction_without_bystanders,
        "independent_correction": independent_correction,
        "indep_less_w_bystanders": indep_less_w_bystanders,
        "w_bystanders_less_wo_bystanders": w_bystanders_less_wo_bystanders,
        "target_locus":guide_seq,
        "perfect_correction":search_strings[0],
        "corrected_locus_with_bystanders": ";".join(search_strings)
        }
    
    return result