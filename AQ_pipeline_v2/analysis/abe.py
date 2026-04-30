import pandas as pd
from utils.sequences import reverse_complement

# ABE-specific metric: correction rate with/without bystanders, any A-to-G, any change

def calculate_correction(allele_table: pd.DataFrame,
                        search_sequences: list[str]) -> tuple[float, float]:

    exact_match_mask = allele_table["Aligned_Sequence"] == search_sequences[0]
    pct_without_bystanders = allele_table[exact_match_mask]["%Reads"].sum()

    any_match_mask = allele_table["Aligned_Sequence"].isin(search_sequences)
    pct_with_bystanders = allele_table[any_match_mask]["%Reads"].sum()
        
    return (pct_without_bystanders, pct_with_bystanders)

def calculate_protospacer_metrics(allele_table: pd.DataFrame,
                                  protospacer: str,
                                  intended_edit: int,
                                  orientation: str) -> tuple[float, float]:
    any_change_in_protospacer = 0
    any_AtoG_change_in_protospacer = 0
    rc_protospacer = reverse_complement(protospacer)

    for _, row in allele_table.iterrows():
        if orientation == "F":
            if row["Aligned_Sequence"][intended_edit-1] == "G":
                any_change_in_protospacer += row["%Reads"]
                only_AtoG = True
                for idx, c in enumerate(row["Aligned_Sequence"]):
                    if c != protospacer[idx]:
                        if not (protospacer[idx] == "A" and c == "G"):
                            only_AtoG = False
                            break
                if only_AtoG:
                    any_AtoG_change_in_protospacer += row["%Reads"]
        elif orientation == "R":
            if row["Aligned_Sequence"][len(rc_protospacer) - intended_edit] == "C":
                any_change_in_protospacer += row["%Reads"]
                only_AtoG = True
                for idx, c in enumerate(row["Aligned_Sequence"]):
                    if c != rc_protospacer[idx]:
                        if not (rc_protospacer[idx] == "T" and c == "C"):
                            only_AtoG = False
                            break
                if only_AtoG:
                    any_AtoG_change_in_protospacer += row["%Reads"]
        else:
            raise ValueError(f"orientation must be 'F' or 'R', got '{orientation}'")

    return(any_AtoG_change_in_protospacer, any_change_in_protospacer,)