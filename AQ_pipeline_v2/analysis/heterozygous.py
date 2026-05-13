import pandas as pd

# Detects heterozygous positions and splits allele tables per allele


def find_het_position(quant_window_df: pd.DataFrame) -> tuple[int | None, str | None, str | None]:
    for idx, col in enumerate(quant_window_df.columns):
        col_data = quant_window_df[col]
        bases_in_range = []
        for base in col_data.index:
            if 0.40 <= col_data[base] <= 0.60:
                bases_in_range.append(base)
        if len(bases_in_range) == 2:                        
            if set(bases_in_range) in [{"A","G"}, {"C","T"}]:
                continue                                     
            return (idx, bases_in_range[0], bases_in_range[1])
    return (None, None, None)

def calculate_het_correction(allele_table_df, search_seqs, het_pos, base1, base2) -> dict:
    base1_mask = allele_table_df["Aligned_Sequence"].str[het_pos] == base1
    base2_mask = allele_table_df["Aligned_Sequence"].str[het_pos] == base2

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