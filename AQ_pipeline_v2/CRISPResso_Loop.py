import logging
from pathlib import Path
from loaders.amplicon_list import load_amplicon_list, find_amplicon_list
from pipeline.crispresso import identify_amplicon, run_crispresso

#Entry point for the CRISPResso loop


log_dir = Path("logs")
log_dir.mkdir(exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(message)s",
    handlers=[
        logging.StreamHandler(),                        #  log to terminal
        logging.FileHandler(log_dir / "crispresso_loop.log"),
    ]
)


def main():
    """Entry point for Stage 1. Iterates over fastq subdirectories, matches each
    to an AmpliconConfig object, and runs CRISPResso on each sample
    """
    error_count = 0
    completed_count = 0
    fastqs_dir = Path("fastqs")
    amplicon_configs = load_amplicon_list(find_amplicon_list())

    for sample_dir in fastqs_dir.iterdir():
        if not sample_dir.is_dir():
            continue
        try:
            logging.info(f"Processing {sample_dir.name}")
            config = identify_amplicon(sample_dir.name, amplicon_configs)
            run_crispresso(config, sample_dir)
            logging.info(f"Done: {sample_dir.name}")
            completed_count += 1
        except Exception as e:
            logging.error(f"Error processing {sample_dir.name}: {e}") 
            error_count += 1
    
    logging.info(f"Samples processed correctly: {completed_count}")
    logging.info(f"Samples encountered with errors: {error_count}")


    
if __name__ == "__main__":
    main()