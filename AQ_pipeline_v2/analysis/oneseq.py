import pandas as pd

# ONESEQ analysis: A-to-G combinations across the first 10bp and full protospacer

def calculate_oneseq(allele_table: pd.DataFrame,
                     first_10bp_seqs: list[str],
                     any_bp_seqs: list[str]) -> tuple[float, float]:
    """calculates the percentage of correction in subsections of the full sequence
    Args:
        allele_table: the allele frequency table from the relevant CRISPResso sample folder.
        first_10bp_seqs: a list of sequences with correction in the first 10 base pairs.
        any_bp_seqs: a list of sequences with correction anywhere in the protospacer.
    Returns:
        tuple[float, float]: returns a tuple of floats containing the percentage of editing in
        the first 10 bp anywhere in the protospacer respectively
    """
        
    first_10_mask = allele_table["Aligned_Sequence"].isin(first_10bp_seqs)
    pct_first_10_bp_editing = allele_table[first_10_mask]["%Reads"].sum()

    any_bp_mask = allele_table["Aligned_Sequence"].isin(any_bp_seqs)
    pct_any_bp_editing = allele_table[any_bp_mask]["%Reads"].sum()

    return(pct_first_10_bp_editing, pct_any_bp_editing)
