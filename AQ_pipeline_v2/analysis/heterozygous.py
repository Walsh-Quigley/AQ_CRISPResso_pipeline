import pandas as pd
from utils.sequences import reverse_complement
import logging

# Detects heterozygous positions and splits allele tables per allele


def find_het_position(quant_window_df: pd.DataFrame) -> tuple[list[int], str | None, str | None]:
    """finds all positions in protospacer that differ from one allele vs the other
    Args:
        quant_window_df: dataframe for the given sample's quantification window
    Returns:
        tuple[list[int], str, str]: returns a tuple containing the list of positions of het nt, and the 
            two bases that rest at the primary het index."""
    het_positions = []

    primary_base1 = None
    primary_base2 = None
    for idx, col in enumerate(quant_window_df.columns):
        col_data = quant_window_df[col]
        bases_in_range = []
        for base in col_data.index:
            if base not in ("A", "C", "G", "T"):
                continue
            if 0.40 <= col_data[base] <= 0.60:
                bases_in_range.append(base)
        # A/G and C/T pairs at ~50/50 are indistinguishable from ABE editing artifacts,
        # so they're filtered. Filtered positions are NOT added to het_positions,
        # meaning they get counted as "edits" in any_AtoG / any_change protospacer walks —
        # the conservative choice when we can't tell het from editing.
        if len(bases_in_range) == 2:                        
            if set(bases_in_range) in [{"A","G"}, {"C","T"}]:
                continue                                     
            het_positions.append(idx)
            if primary_base1 is None and primary_base2 is None:
                primary_base1 = bases_in_range[0]
                primary_base2 = bases_in_range[1]
    return (het_positions, primary_base1, primary_base2)

def calculate_het_correction(allele_table_df: pd.DataFrame,
                            search_seqs: list[str],
                            het_pos: list[int],
                            base1: str,
                            base2: str) -> dict:
    """calculates correction with bystanders and correction without bystanders for heterozygous samples.
    Args:
        allele_table_df: the dataframe for a given sample's allele frequency table
        search_seqs: the list of strings to search for in the allele frequency table
        het_pos: the heterozygous positions in the protospacer
        base1: nt at primary het_pos for allele 1
        base2: nt at primary het_pos for allele 2
    Returns:
        A dictionary containing corrections with bystanders and without bystanders broken up by allele
    Note:
        it is of note that reads_w_baseX are "reads with tolerated bystanders for base X". It is NOT "read with base X"
    """
    unsorted_reads_pct = 0.0
    primary_het_pos = het_pos[0]
    total_reads_base1, total_reads_base2 = 0, 0
    reads_wo_base1, reads_wo_base2 = 0, 0
    reads_w_base1, reads_w_base2 = 0, 0
    alignment_shift_reads_pct = 0.0

    for _, row in allele_table_df.iterrows():
        if len(row["Aligned_Sequence"]) != len(search_seqs[0]):
            alignment_shift_reads_pct += row["%Reads"]
            continue
        if row["Aligned_Sequence"][primary_het_pos] == base1:
            total_reads_base1 += row["%Reads"]
            if row["Aligned_Sequence"] in search_seqs:
                reads_w_base1 += row["%Reads"]
            if row["Aligned_Sequence"] == search_seqs[0]:
                reads_wo_base1 += row["%Reads"]
        elif row["Aligned_Sequence"][primary_het_pos] == base2:
            total_reads_base2 += row["%Reads"]
            if row["Aligned_Sequence"] in search_seqs:
                reads_w_base2 += row["%Reads"]
            if row["Aligned_Sequence"] == search_seqs[0]:
                reads_wo_base2 += row["%Reads"]
        else:
            unsorted_reads_pct += row["%Reads"]
            continue
    
    if unsorted_reads_pct > 3:
        logging.warning(
            f"{unsorted_reads_pct:.2f}% of aligned reads had a base other than "
            f"'{base1}' or '{base2}' at het position {primary_het_pos + 1} — "
            f"these reads were excluded from per-allele metrics. "
            f"High values may indicate sequencing errors or a third allele."
        )

    if alignment_shift_reads_pct > 3:
        logging.warning(
            f"Skipped {alignment_shift_reads_pct:.2f}% of reads with alignment shifts "
            f"(likely insertions in/near the protospacer region). "
            f"These reads cannot be reliably analyzed for per-protospacer metrics."
        )
            

    pct_wo_base1 = ((reads_wo_base1 / total_reads_base1) * 100)  if total_reads_base1 > 0 else 0.0
    pct_wo_base2 = ((reads_wo_base2 / total_reads_base2) * 100)  if total_reads_base2 > 0 else 0.0
    pct_w_base1 = ((reads_w_base1 / total_reads_base1) * 100)  if total_reads_base1 > 0 else 0.0
    pct_w_base2 = ((reads_w_base2 / total_reads_base2) * 100)  if total_reads_base2 > 0 else 0.0

    results_dict = {
        "correction_wo_bystanders_allele1": pct_wo_base1,
        "correction_w_bystanders_allele1": pct_w_base1,
        "correction_wo_bystanders_allele2": pct_wo_base2,
        "correction_w_bystanders_allele2": pct_w_base2,    
        "total_pct_allele1": total_reads_base1,
        "total_pct_allele2": total_reads_base2,
        "total_pct_allele1": total_reads_base1,
        "total_pct_allele2": total_reads_base2,
    }

    return results_dict

def calculate_het_protospacer_metrics(allele_table: pd.DataFrame,   
                                      protospacer: str,
                                      intended_edit: int,
                                      orientation: str,
                                      het_pos: list[int],
                                      base1: str,
                                      base2: str) -> dict:
    """Calculates correction with A to G change and correction with any change in protospacer for
        heterozygous samples.
    Args:
        allele_table: dataframe containing the allele frequency table for a given sample
        protospacer: protospacer for the sample
        intended_edit: the intended edit index
        orientation: orientation of the protospacer relative to the corresponding amplicon
        het_pos: the heterozygous positions in the protospacer
        base1: nt at primary het_pos for allele 1
        base2: nt at primary het_pos for allele 2
    Returns:
        returns a dictionary containing correction with A to G change and correction with any
            change in protospacer, broken out by allele.
    Raises:
        ValueError: raises value error if orientation is neither forward nor reverse
    """
    
    any_change_in_protospacer_base1 = 0
    any_change_in_protospacer_base2 = 0
    any_AtoG_change_in_protospacer_base1 = 0
    any_AtoG_change_in_protospacer_base2 = 0

    total_reads_base1 = 0
    total_reads_base2 = 0

    alignment_shift_reads_pct = 0.0


    rc_protospacer = reverse_complement(protospacer)

    primary_het_pos = het_pos[0]
    for _, row in allele_table.iterrows():
        if len(row["Aligned_Sequence"]) != len(protospacer):
            alignment_shift_reads_pct += row["%Reads"]
            continue
        if orientation == "F":
            if row["Aligned_Sequence"][primary_het_pos] == base1:
                total_reads_base1 += row["%Reads"]
                if row["Aligned_Sequence"][intended_edit-1] == "G":
                    any_change_in_protospacer_base1 += row["%Reads"]
                    only_AtoG = True
                    for idx, c in enumerate(row["Aligned_Sequence"]):
                        if idx in het_pos:
                            continue
                        if c != protospacer[idx]:
                            if not (protospacer[idx] == "A" and c == "G"):
                                only_AtoG = False
                                break
                    if only_AtoG:
                        any_AtoG_change_in_protospacer_base1 += row["%Reads"]
            elif row["Aligned_Sequence"][primary_het_pos] == base2:
                total_reads_base2 += row["%Reads"]
                if row["Aligned_Sequence"][intended_edit-1] == "G":
                    any_change_in_protospacer_base2 += row["%Reads"]
                    only_AtoG = True
                    for idx, c in enumerate(row["Aligned_Sequence"]):
                        if idx in het_pos:
                            continue
                        if c != protospacer[idx]:
                            if not (protospacer[idx] == "A" and c == "G"):
                                only_AtoG = False
                                break
                    if only_AtoG:
                        any_AtoG_change_in_protospacer_base2 += row["%Reads"]
            else:
                continue
        elif orientation == "R":
            if row["Aligned_Sequence"][primary_het_pos] == base1:
                total_reads_base1 += row["%Reads"]
                if row["Aligned_Sequence"][len(rc_protospacer) - intended_edit] == "C":
                    any_change_in_protospacer_base1 += row["%Reads"]
                    only_AtoG = True
                    for idx, c in enumerate(row["Aligned_Sequence"]):
                        if idx in het_pos:
                            continue
                        if c != rc_protospacer[idx]:
                            if not (rc_protospacer[idx] == "T" and c == "C"):
                                only_AtoG = False
                                break
                    if only_AtoG:
                        any_AtoG_change_in_protospacer_base1 += row["%Reads"]
            elif row["Aligned_Sequence"][primary_het_pos] == base2:
                total_reads_base2 += row["%Reads"]
                if row["Aligned_Sequence"][len(rc_protospacer) - intended_edit] == "C":
                    any_change_in_protospacer_base2 += row["%Reads"]
                    only_AtoG = True
                    for idx, c in enumerate(row["Aligned_Sequence"]):
                        if idx in het_pos:
                            continue
                        if c != rc_protospacer[idx]:
                            if not (rc_protospacer[idx] == "T" and c == "C"):
                                only_AtoG = False
                                break
                    if only_AtoG:
                        any_AtoG_change_in_protospacer_base2 += row["%Reads"]
            else:
                continue
        else:
            raise ValueError(f"orientation must be 'F' or 'R', got '{orientation}'")
    
    if alignment_shift_reads_pct > 3:
        logging.warning(
            f"Skipped {alignment_shift_reads_pct:.2f}% of reads with alignment shifts "
            f"(likely insertions in/near the protospacer region). "
            f"These reads cannot be reliably analyzed for per-protospacer metrics."
        )

    AtoG_base1 = ((any_AtoG_change_in_protospacer_base1 / total_reads_base1) * 100) if total_reads_base1 > 0 else 0.0
    AtoG_base2 = ((any_AtoG_change_in_protospacer_base2 / total_reads_base2) * 100) if total_reads_base2 > 0 else 0.0
    any_change_base_1 = ((any_change_in_protospacer_base1 / total_reads_base1) * 100) if total_reads_base1 > 0 else 0.0
    any_change_base_2 = ((any_change_in_protospacer_base2 / total_reads_base2) * 100) if total_reads_base2 > 0 else 0.0
    
    results_dict = {
        "correction_with_any_AtoG_change_allele1": AtoG_base1,
        "correction_with_any_change_in_protospacer_allele1": any_change_base_1,
        "correction_with_any_AtoG_change_allele2": AtoG_base2,
        "correction_with_any_change_in_protospacer_allele2": any_change_base_2,
    }
    
    
    return results_dict