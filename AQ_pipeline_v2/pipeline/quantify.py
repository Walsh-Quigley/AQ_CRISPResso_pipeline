from pathlib import Path
from glob import glob
import re
from config import AmpliconConfig
from loaders.crispresso_output import read_mapping_stats, read_allele_table
from utils.sequences import generate_search_sequences
from analysis.abe import calculate_correction, calculate_protospacer_metrics

# Stage 2 -> parses CRISPResso outputs, calls analysis modules,
# assembles final reasult

def quantify_sample(amplicon_row: AmpliconConfig, crispresso_dir: Path) -> dict:

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
        "sample": re.sub(r'(_L\d{3})?-ds\..*', '', crispresso_dir.name),
        "reads_total": reads_total, #B
        "reads_aligned": reads_aligned, #C
        "correction_without_bystanders": without_bystanders, #D
        "correction_with_tolerated_bystanders": with_bystanders, #E
        "correction_with_any_AtoG_change": any_AtoG, #F
        "correction_with_any_change_in_protospacer": any_change, #G
        "column E minus column D": round(with_bystanders - without_bystanders, 2),
        "column F minus column E": round(any_AtoG - with_bystanders, 2),
        "column G minus column F": round(any_change - any_AtoG, 2),
        "target_locus":amplicon_row.protospacer,
        "perfect_correction": search_seqs[0],
        "corrected_locus_with_bystanders": ";".join(search_seqs)
    }    
