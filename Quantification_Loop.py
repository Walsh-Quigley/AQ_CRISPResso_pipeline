import logging
import pandas as pd
from pathlib import Path
from loaders.amplicon_list import load_amplicon_list, find_amplicon_list
from pipeline.crispresso import identify_amplicon
from pipeline.quantify import quantify_sample
from loaders.exports import generate_prism_csv, generate_prism_csv_het

#Entry point for the Quantification stage of the pipeline

SKIP_DIRS = {"__pycache__", ".git", ".DS_Store"}


CANONICAL_ABE_COLUMNS = [
    "sample",
    "reads_total",
    "reads_aligned",
    "correction_without_bystanders",
    "correction_with_tolerated_bystanders",
    "correction_with_any_AtoG_change",
    "correction_with_any_change_in_protospacer",
    "w_bystanders_minus_wo_bystanders",
    "any_AtoG_minus_w_bystanders",
    "any_change_minus_any_AtoG",
    "correction_wo_bystanders_allele1",
    "correction_w_bystanders_allele1",
    "correction_wo_bystanders_allele2",
    "correction_w_bystanders_allele2",
    "correction_with_any_AtoG_change_allele1",
    "correction_with_any_change_in_protospacer_allele1",
    "correction_with_any_AtoG_change_allele2",
    "correction_with_any_change_in_protospacer_allele2",
    "w_bystanders_minus_wo_bystanders_allele1",
    "w_bystanders_minus_wo_bystanders_allele2",
    "any_AtoG_minus_w_bystanders_allele1",
    "any_AtoG_minus_w_bystanders_allele2",
    "any_change_minus_any_AtoG_allele1",
    "any_change_minus_any_AtoG_allele2",
    "het_position",
    "het_alleles",
    "reads_aligned_allele1",
    "reads_aligned_allele2",
    "target_locus",
    "perfect_correction",
    "corrected_locus_with_bystanders",
]


log_dir = Path("logs")
log_dir.mkdir(exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(message)s",
    handlers=[
        logging.StreamHandler(),                        #log toterminal
        logging.FileHandler(log_dir / "quantification_loop.log"),
    ]
)

def main():
    """Entry point for Stage 2. Iterates over all fastq subdirectories, matches each
    to an AmpliconConfig object, and performs data analysis on each sample"""
    failed_samples = []
    error_count = 0
    completed_count = 0
    results_by_type = {"ABE": [], "ONESEQ": []}
    fastqs_dir = Path("fastqs")
    amplicon_configs = load_amplicon_list(find_amplicon_list())
    for sample_dir in fastqs_dir.iterdir():
        if not sample_dir.is_dir():
            continue
        if sample_dir.name in SKIP_DIRS:
            continue
        try:
            config = identify_amplicon(sample_dir.name, amplicon_configs)
            logging.info(f"Processing {sample_dir.name}")
            result = quantify_sample(config, sample_dir)
            if config.intended_edit == "ONESEQ":
                results_by_type["ONESEQ"].append(result)
            else:
                results_by_type["ABE"].append(result)
            logging.info(f"Done: {sample_dir.name}")
            completed_count += 1
        except Exception as e:
            logging.error(f"Error processing {sample_dir.name}: {e}")
            error_count += 1
            failed_samples.append({
                "sample": sample_dir.name,
                "error_type": type(e).__name__,
                "error_message": str(e),
            })

    abe_df = None
    for each in results_by_type:
        if results_by_type[each]:
            if each == "ABE":
                abe_df = pd.DataFrame(results_by_type["ABE"])
                abe_df = abe_df.sort_values(by="sample")
                known = [c for c in CANONICAL_ABE_COLUMNS if c in abe_df.columns]
                unknown = [c for c in abe_df.columns if c not in CANONICAL_ABE_COLUMNS]
                if unknown:
                    logging.warning(f"Unexpected columns not in canonical order, appended at end: {unknown}")
                abe_df = abe_df.reindex(columns=known + unknown)
                abe_df.to_csv("ABE_Quantification_Summary.csv", index=False)
            elif each == "ONESEQ":
                oneseq_df = pd.DataFrame(results_by_type["ONESEQ"])
                oneseq_df = oneseq_df.sort_values(by="sample")
                oneseq_df.to_csv("ONESEQ_Quantification_Summary.csv", index=False)
            else:
                raise ValueError(f"Unknown editor type")
    if results_by_type["ABE"]:
        while True:
            make_prism = input("Would you like to generate a Prism input file? (y/n)").strip().lower()
            if make_prism in ("yes", "no", "y", "n"):
                break
            print("Please enter (y/n)")
        if make_prism == "y" or make_prism == "yes":
            if "het_position" in abe_df.columns:
                het_df = abe_df[abe_df["het_position"].notna()].copy()
                non_het_df = abe_df[abe_df["het_position"].isna()].copy()
            else:
                het_df = pd.DataFrame()
                non_het_df = abe_df.copy()

            if not non_het_df.empty:
                prism_df = generate_prism_csv(non_het_df)
                prism_df.to_csv("Prism_Input.csv", index=False)

            if not het_df.empty:
                allele1_df, allele2_df = generate_prism_csv_het(het_df)
                allele1_df.to_csv("Prism_Input_het_allele1.csv", index=False)
                allele2_df.to_csv("Prism_Input_het_allele2.csv", index=False)


    if failed_samples:
        failed_df = pd.DataFrame(failed_samples).sort_values("sample")
        failed_df.to_csv("failed_samples.csv", index=False)
        logging.warning(f"{len(failed_samples)} sample(s) failed — see failed_samples.csv")
        for f in failed_samples:
            logging.info(f"  FAILED: {f['sample']} — {f['error_type']}: {f['error_message']}")

    logging.info(f"Samples processed correctly: {completed_count}")
    logging.info(f"Samples encountered with errors: {error_count}")

            


if __name__ == "__main__":
    main()