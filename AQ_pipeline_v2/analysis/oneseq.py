from utils.sequences import reverse_complement
import pandas as pd

# ONE_seq analysis: A-to-G combinations across the first 10bp and full protospacer

def calculate_oneseq(allele_table: pd.DataFrame,
                     first_10bp_seqs: list[str],
                     any_bp_seqs: list[str]) -> tuple[float, float]:
        
    first_10_mask = allele_table["Aligned_Sequence"].isin(first_10bp_seqs)
    pct_first_10_bp_editing = allele_table[first_10_mask]["%Reads"].sum()

    any_bp_mask = allele_table["Aligned_Sequence"].isin(any_bp_seqs)
    pct_any_bp_editing = allele_table[any_bp_mask]["%Reads"].sum()

    return(pct_first_10_bp_editing, pct_any_bp_editing)
