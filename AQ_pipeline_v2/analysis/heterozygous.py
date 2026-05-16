import pandas as pd
from utils.sequences import reverse_complement

# Detects heterozygous positions and splits allele tables per allele


def find_het_position(quant_window_df: pd.DataFrame) -> tuple[list[int], str | None, str | None]:
    het_positions = []
    primary_base1 = None
    primary_base2 = None
    for idx, col in enumerate(quant_window_df.columns):
        col_data = quant_window_df[col]
        bases_in_range = []
        for base in col_data.index:
            if 0.40 <= col_data[base] <= 0.60:
                bases_in_range.append(base)
        if len(bases_in_range) == 2:                        
            if set(bases_in_range) in [{"A","G"}, {"C","T"}]:
                continue                                     
            het_positions.append(idx)
            if primary_base1 == None and primary_base2 == None:
                primary_base1 = bases_in_range[0]
                primary_base2 = bases_in_range[1]
    return (het_positions, primary_base1, primary_base2)

def calculate_het_correction(allele_table_df: pd.DataFrame,
                            search_seqs: list[str],
                            het_pos: list[int],
                            base1: str,
                            base2: str) -> dict:
    
    primary_het_pos = het_pos[0]
    base1_mask = allele_table_df["Aligned_Sequence"].str[primary_het_pos] == base1
    base2_mask = allele_table_df["Aligned_Sequence"].str[primary_het_pos] == base2

    total_reads_base1 = allele_table_df[base1_mask]["%Reads"].sum()
    total_reads_base2 = allele_table_df[base2_mask]["%Reads"].sum()

    any_match_mask = allele_table_df["Aligned_Sequence"].isin(search_seqs)
    reads_w_base1 = allele_table_df[base1_mask & any_match_mask]["%Reads"].sum()
    reads_w_base2 = allele_table_df[base2_mask & any_match_mask]["%Reads"].sum()

    without_bystanders_match = allele_table_df["Aligned_Sequence"] == search_seqs[0]
    reads_wo_base1 = allele_table_df[base1_mask & without_bystanders_match]["%Reads"].sum()
    reads_wo_base2 = allele_table_df[base2_mask & without_bystanders_match]["%Reads"].sum()

    pct_wo_base1 = ((reads_wo_base1 / total_reads_base1) * 100)  if total_reads_base1 > 0 else 0.0
    pct_wo_base2 = ((reads_wo_base2 / total_reads_base2) * 100)  if total_reads_base2 > 0 else 0.0
    pct_w_base1 = ((reads_w_base1 / total_reads_base1) * 100)  if total_reads_base1 > 0 else 0.0
    pct_w_base2 = ((reads_w_base2 / total_reads_base2) * 100)  if total_reads_base2 > 0 else 0.0

    resultsDict = {
    "correction_w_bystanders_allele1": pct_w_base1,
    "correction_wo_bystanders_allele1": pct_wo_base1,
    "correction_w_bystanders_allele2": pct_w_base2,
    "correction_wo_bystanders_allele2": pct_wo_base2,
    }

    return resultsDict

def calculate_het_protospacer_metrics(allele_table: pd.DataFrame, 
                                      protospacer: str,
                                      intended_edit: int,
                                      orientation: str,
                                      het_pos: list[int],
                                      base1: str,
                                      base2: str) -> dict:
    
    any_change_in_protospacer_base1 = 0
    any_change_in_protospacer_base2 = 0
    any_AtoG_change_in_protospacer_base1 = 0
    any_AtoG_change_in_protospacer_base2 = 0

    total_reads_base1 = 0
    total_reads_base2 = 0

    rc_protospacer = reverse_complement(protospacer)

    primary_het_position = het_pos[0]
    for _, row in allele_table.iterrows():
        if orientation == "F":
            if row["Aligned_Sequence"][primary_het_position] == base1:
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
            elif row["Aligned_Sequence"][primary_het_position] == base2:
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
            if row["Aligned_Sequence"][primary_het_position] == base1:
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
            elif row["Aligned_Sequence"][primary_het_position] == base2:
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
    
    AtoG_base1 = ((any_AtoG_change_in_protospacer_base1 / total_reads_base1) * 100) if total_reads_base1 > 0 else 0.0
    AtoG_base2 = ((any_AtoG_change_in_protospacer_base2 / total_reads_base2) * 100) if total_reads_base2 > 0 else 0.0
    Any_change_base_1 = ((any_change_in_protospacer_base1 / total_reads_base1) * 100) if total_reads_base1 > 0 else 0.0
    Any_change_base_2 = ((any_change_in_protospacer_base2 / total_reads_base2) * 100) if total_reads_base2 > 0 else 0.0
    
    resultsDict = {
        "correction_with_any_AtoG_change_allele1": AtoG_base1,
        "correction_with_any_change_in_protospacer_allele1": Any_change_base_1,
        "correction_with_any_AtoG_change_allele2": AtoG_base2,
        "correction_with_any_change_in_protospacer_allele2": Any_change_base_2,
    }
    
    
    return resultsDict