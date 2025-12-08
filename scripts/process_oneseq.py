#scripts/process_oneseq.py

import os
import re
import logging
from .generate_search_sequences import A_to_G_sequences
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

    #generate search sequences for the editing window
    search_strings_first_10bp, search_strings_anyA = A_to_G_sequences(guide_seq, orientation)
    logging.info(f"The search strings generated for ONESEQ are: {search_strings_first_10bp}")

    #filter alleles file for A to G changes in the first 10bp
    A_to_G_in_first_ten, meaningless = filter_alleles_file(search_strings_first_10bp, directory_path)
    logging.info(f"editing at this site in the first 10 bases: {A_to_G_in_first_ten}")

    #filter alleles file for any A to G changes
    A_to_G_anywhere, meaningless = filter_alleles_file(search_strings_anyA, directory_path)
    


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
        "Percent_of_reads_with_A>G_in_first_10bp":A_to_G_in_first_ten,
        "Percent_of_reads_with_A>G_in_protospacer":A_to_G_anywhere,
        "guide_seq": guide_seq,
        "A>G_10bp_search_sequences": ";".join(search_strings_first_10bp),
        "A>G_any_search_sequences": ";".join(search_strings_anyA)
        })
    
    return result