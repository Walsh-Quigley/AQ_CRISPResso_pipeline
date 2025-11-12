#scripts/process_oneseq.py

import os
import re
import logging
from .generate_search_sequences import first_10bp_A_to_G_sequences
from .read_extraction import read_extraction
from .filter_alleles_file import filter_alleles_file



def process_oneseq(directory_path, guide_seq, orientation):
    """
    Process ONE-seq target amplicon sequencing samples
    
    Args:
        directory_path: Full path to the directory to process
        guide_seq: the off target locus
        orientation: Orientation of the guide relative to amplicon
        
    Returns:
        dict: Result dictionary with ONESEQ metrics
    """

    logging.info(f"Processing {directory_path} for ONE-seq analysis.")

    #generate search sequences
    search_strings = first_10bp_A_to_G_sequences(guide_seq, orientation)
    logging.info(f"The search strings generated for ONESEQ are: {search_strings}")

    #filter alleles file
    correction_with_bystander, correction_without_bystanders = filter_alleles_file(search_strings, directory_path)
    logging.info(f"editing at this site in the first 10 bases: {correction_with_bystander}")

    # Extract sample name and read counts
    directory_name = os.path.basename(directory_path.rstrip('/'))
    sample_name = re.sub(r'(_L\d{3})?-ds\..*', '', directory_name)

    # get read counts
    reads_aligned, reads_total = read_extraction(directory_path)

    #build results dictionary
    result = ({
        "sample":sample_name,
        "reads_aligned": reads_aligned,
        "reads_total": reads_total,
        "Percent_reads_with_edit_in_edit_window":correction_with_bystander,
        "guide_seq": guide_seq,
        "search_sequences": ";".join(search_strings)
        })
    
    return result