import logging
import pandas as pd
from pathlib import Path
from loaders.amplicon_list import load_amplicon_list, find_amplicon_list
from pipeline.crispresso import identify_amplicon
from pipeline.quantify import quantify_sample

#Entry point for the Quantification stage of the pipeline

log_dir = Path("logs")
log_dir.mkdir(exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(message)s",
    handlers=[
        logging.StreamHandler(),                        # terminal
        logging.FileHandler(log_dir / "quantification_loop.log"),
    ]
)

def main():
    error_count = 0
    completed_count = 0
    results = []
    fastqs_dir = Path("fastqs")
    amplicon_configs = load_amplicon_list(find_amplicon_list())
    for sample_dir in fastqs_dir.iterdir():
        if not sample_dir.is_dir():
            continue
        try:
            logging.info(f"Processing {sample_dir.name}")
            config = identify_amplicon(sample_dir.name, amplicon_configs)
            results.append(quantify_sample(config, sample_dir))
            logging.info(f"Done: {sample_dir.name}")
            completed_count += 1
        except Exception as e:
            logging.error(f"Error processing {sample_dir.name}: {e}")  # something went wrong
            error_count += 1

    df = pd.DataFrame(results)
    df = df.sort_values(by="sample")
    df.to_csv("quantification_summary.csv", index=False)
    logging.info(f"Samples processed correctly: {completed_count}")
    logging.info(f"Samples encountered with errors: {error_count}")

            


if __name__ == "__main__":
    main()