from pathlib import Path
from glob import glob
import pandas as pd
import re
from config import AmpliconConfig
from loaders.crispresso_output import read_mapping_stats, read_allele_table, read_quant_window
from utils.sequences import generate_search_sequences, generate_oneseq_search_sequences, reverse_complement
from analysis.abe import calculate_correction, calculate_protospacer_metrics
from analysis.oneseq import calculate_oneseq
from analysis.heterozygous import calculate_het_correction, calculate_het_protospacer_metrics, find_het_position


# Stage 2 -> parses CRISPResso outputs, calls analysis modules,
# assembles final result
def quantify_abe_sample(amplicon_row: AmpliconConfig, 
                        sample_name: str,
                        allele_table_df: pd.DataFrame,
                        reads_total: int,
                        reads_aligned: int) -> dict:
    """creates the dictionary for samples with ABE editing
    Args:
        amplicon_row: the AmpliconConfig object for the given CRISPResso sample
        sample_name: name of the CRISPResso sample
        allele_table_df: the dataframe containing the allele frequency table information
        reads_total: the number of total reads for a sample
        reads_aligned: the number of reads aligned to the given amplicon
    Returns:
        dict: the dictionary containing all the information relevant to the ABE analysis of a crispresso sample
        """

    search_seqs = generate_search_sequences(
        protospacer=amplicon_row.protospacer,
        intended_edit=amplicon_row.intended_edit,
        tolerated_edits=amplicon_row.tolerated_edits,
        orientation=amplicon_row.orientation
    )

    without_bystanders, with_bystanders = calculate_correction(allele_table_df, search_seqs)

    any_AtoG, any_change = calculate_protospacer_metrics(
        allele_table_df,
        amplicon_row.protospacer,
        amplicon_row.intended_edit,
        amplicon_row.orientation
    )

    return {
        "sample": re.sub(r'(_L\d{3})?-ds\..*', '', sample_name),
        "reads_total": reads_total, #B
        "reads_aligned": reads_aligned, #C
        "correction_without_bystanders": without_bystanders, #D
        "correction_with_tolerated_bystanders": with_bystanders, #E
        "correction_with_any_AtoG_change": any_AtoG, #F
        "correction_with_any_change_in_protospacer": any_change, #G
        "w_bystanders_minus_wo_bystanders": round(with_bystanders - without_bystanders, 2),
        "any_AtoG_minus_w_bystanders": round(any_AtoG - with_bystanders, 2),
        "any_change_minus_any_AtoG": round(any_change - any_AtoG, 2),
        "target_locus":amplicon_row.protospacer,
        "perfect_correction": search_seqs[0],
        "corrected_locus_with_bystanders": ";".join(search_seqs)
    }    

def quantify_het_sample(amplicon_row: AmpliconConfig,
                        sample_name: str,
                        allele_table_df: pd.DataFrame,
                        reads_total: int,
                        reads_aligned: int,
                        het_pos: list[int],
                        het_base1: str,
                        het_base2: str) -> dict:
    """Creates the result dictionary for heterozygous ABE samples.
    Args:
        amplicon_row: the AmpliconConfig object for a given sample
        sample_name: the name of the current sample
        allele_table_df: the dataframe for the current sample's allele frequency table
        reads_total: total number of reads in the current sample's fastq
        reads_aligned: total number of ALIGNED reads in the current sample's fastq
        het_pos: positions of the heterozygous differences, het_pos[0] is the primary het_pos,
            which determines allele sorting downstream.
        het_base1: nt at het_pos for allele 1
        het_base2: nt at het_pos for allele 2
    Returns:
        A dictionary containing all sample metrics, broken out by allele for relevant samples.
    Note:
        het_base1 will always correspond to the protospacer allele
    """
    base = None
    if amplicon_row.orientation == "F":
        base = amplicon_row.protospacer[het_pos[0]]
    elif amplicon_row.orientation == "R":
        base = reverse_complement(amplicon_row.protospacer)[het_pos[0]]
    if base == het_base2:
        het_base1, het_base2 = het_base2, het_base1

    search_seqs = generate_search_sequences(
        protospacer=amplicon_row.protospacer,
        intended_edit=amplicon_row.intended_edit,
        tolerated_edits=amplicon_row.tolerated_edits,
        orientation=amplicon_row.orientation
    )

    het_correction_dict = calculate_het_correction(allele_table_df,
                                                    search_seqs,
                                                    het_pos,
                                                    het_base1,
                                                    het_base2)
    het_protospacer_metrics_dict = calculate_het_protospacer_metrics(allele_table_df,
                                                                     amplicon_row.protospacer,
                                                                     amplicon_row.intended_edit,
                                                                     amplicon_row.orientation,
                                                                     het_pos,
                                                                     het_base1,
                                                                     het_base2)

    without_bystanders, with_bystanders = calculate_correction(allele_table_df, search_seqs)

    any_AtoG, any_change = calculate_protospacer_metrics(
        allele_table_df,
        amplicon_row.protospacer,
        amplicon_row.intended_edit,
        amplicon_row.orientation
    )

    return {
        "sample": re.sub(r'(_L\d{3})?-ds\..*', '', sample_name),
        "reads_total": reads_total, #B
        "reads_aligned": reads_aligned, #C
        "correction_without_bystanders": without_bystanders, #D
        "correction_with_tolerated_bystanders": with_bystanders, #E
        "correction_with_any_AtoG_change": any_AtoG, #F
        "correction_with_any_change_in_protospacer": any_change, #G
        "w_bystanders_minus_wo_bystanders": round(with_bystanders - without_bystanders, 2),
        "any_AtoG_minus_w_bystanders": round(any_AtoG - with_bystanders, 2),
        "any_change_minus_any_AtoG": round(any_change - any_AtoG, 2),
        **het_correction_dict,
        **het_protospacer_metrics_dict,
        "w_bystanders_minus_wo_bystanders_allele1": round(het_correction_dict["correction_w_bystanders_allele1"] - het_correction_dict["correction_wo_bystanders_allele1"], 2),
        "w_bystanders_minus_wo_bystanders_allele2": round(het_correction_dict["correction_w_bystanders_allele2"] - het_correction_dict["correction_wo_bystanders_allele2"], 2),
        "any_AtoG_minus_w_bystanders_allele1": round(het_protospacer_metrics_dict["correction_with_any_AtoG_change_allele1"] - het_correction_dict["correction_w_bystanders_allele1"], 2),
        "any_AtoG_minus_w_bystanders_allele2": round(het_protospacer_metrics_dict["correction_with_any_AtoG_change_allele2"] - het_correction_dict["correction_w_bystanders_allele2"], 2),
        "any_change_minus_any_AtoG_allele1": round(het_protospacer_metrics_dict["correction_with_any_change_in_protospacer_allele1"] - het_protospacer_metrics_dict["correction_with_any_AtoG_change_allele1"], 2),
        "any_change_minus_any_AtoG_allele2": round(het_protospacer_metrics_dict["correction_with_any_change_in_protospacer_allele2"] - het_protospacer_metrics_dict["correction_with_any_AtoG_change_allele2"], 2),
        "het_position": het_pos[0] + 1,
        "het_alleles": f"{het_base1}/{het_base2}",
        "target_locus":amplicon_row.protospacer,
        "perfect_correction": search_seqs[0],
        "corrected_locus_with_bystanders": ";".join(search_seqs)
    }    

def quantify_oneseq_sample(amplicon_row: AmpliconConfig, 
                            sample_name: str,
                            allele_table_df: pd.DataFrame,
                            reads_total: int,
                            reads_aligned: int) -> dict:
    """creates a dictionary containing relevant analytics about a given CRISPResso sample.
    Args:
        amplicon_row: the AmpliconConfig object for the specific sample
        sample_name: the name of the sample analysis is being performed on
        allele_table_df: the dataframe containing the information from the Allele Frequency Table from a sample
        reads_total: the number of total reads in a fastq
        reads_aligned: the number of reads that aligned in the fastq
    Returns:
        dict: returns a dictionary of all information from a given ONESEQ sample"""
    first_10_seqs, full_seqs = generate_oneseq_search_sequences(
        protospacer=amplicon_row.protospacer,
        orientation=amplicon_row.orientation
    )

    first_10_pct, full_sequence_pct = calculate_oneseq(allele_table_df, first_10_seqs, full_seqs)
    
    return {
        "sample": re.sub(r'(_L\d{3})?-ds\..*', '', sample_name),
        "reads_total": reads_total,
        "reads_aligned": reads_aligned,
        "pct_AtoG_first_10bp": first_10_pct,
        "pct_AtoG_anywhere": full_sequence_pct,
        "guide_seq": amplicon_row.protospacer,
        "search_sequences_first_10bp": ";".join(first_10_seqs),
        "search_sequences_any": ";".join(full_seqs)
    }

def quantify_sample(amplicon_row: AmpliconConfig, crispresso_dir: Path) -> dict:
    """the guiding path for the sample quantification, determined by what kind of editor
    Args:
        amplicon_row: the AmpliconConfig object for the given CRISPResso sample
        crispresso_dir: the path to the crispresso directory
    Returns:
        dict: the result dictionary passed through from the relevent analysis branch
    Raises:
        FileNotFoundError: no CRISPResso output folder found in the directory
        ValueError: multiple allele tables found in the CRISPResso subfolder
        ValueError: unknown editor type"""
    
    matches = glob(str(crispresso_dir / "CRISPResso_on_*"))
    if not matches:
        raise FileNotFoundError(f"No CRISPResso output folder found in {crispresso_dir}")
    crispresso_subfolder = Path(matches[0])

    allele_file_matches = glob(str(crispresso_subfolder / "Alleles_frequency_table_around_sgRNA_*.txt"))
    if len(allele_file_matches) != 1:
        raise ValueError(f"Expected exactly one allele table in {crispresso_subfolder}, found {len(allele_file_matches)}")

    allele_file = Path(allele_file_matches[0])    

    stats_file = crispresso_subfolder / "CRISPResso_mapping_statistics.txt"
    
    reads_total, reads_aligned = read_mapping_stats(stats_file)

    allele_table_df = read_allele_table(allele_file)

    if amplicon_row.intended_edit == "ONESEQ":
        results_dict = quantify_oneseq_sample(amplicon_row, crispresso_dir.name, allele_table_df, reads_total, reads_aligned)
    else:
        quant_window = crispresso_subfolder / "Quantification_window_nucleotide_percentage_table.txt"
        quant_window_df = read_quant_window(quant_window)
        het_pos, base1, base2 = find_het_position(quant_window_df)
        if amplicon_row.editor == "ABE" and het_pos:
            results_dict = quantify_het_sample(amplicon_row, crispresso_dir.name, allele_table_df, reads_total, reads_aligned, het_pos, base1, base2)
        elif amplicon_row.editor == "ABE":
            results_dict = quantify_abe_sample(amplicon_row, crispresso_dir.name, allele_table_df, reads_total, reads_aligned)
        else:
            raise ValueError(f"Unknown editor type: {amplicon_row.editor}")
    
    return results_dict
