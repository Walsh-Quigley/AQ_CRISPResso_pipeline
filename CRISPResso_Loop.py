import os
import logging
import sys
from scripts.logging_setup import setup_logging
from scripts.banner import print_banner
from scripts.verify_amplicon_list import verify_amplicon_list
from scripts.identify_amplicon import identify_amplicon
from scripts.CRISPResso_inputs import CRISPResso_inputs
from scripts.run_CRISPResso import run_CRISPResso



def main():
    # Set up logging
    log_file = setup_logging()

    #print the super cool banner
    print_banner()

    #Get a list of all the names in the amplicon_list
    try:
        amplicon_names = verify_amplicon_list()  # Auto-search for amplicon_list file
        logging.info(f"Amplicon names: {amplicon_names}")
    except ValueError as e:
        logging.error(f"Amplicon list verification failed: {e}")
        sys.exit(1)
    except FileNotFoundError as e:
        logging.error(f"Amplicon list file error: {e}")
        sys.exit(1)

    # Move into the /fastqs directory
    fastqs_dir = os.path.join(os.getcwd(), "fastqs")
    os.chdir(fastqs_dir)
    logging.info(f"Changed to directory: {fastqs_dir}")

    # Track processing statistics
    processed_count = 0
    skipped_count = 0
    error_count = 0

    # Track specific error types
    error_types = {
        'empty_fastq': 0,
        'spaces_in_path': 0,
        'no_fastq': 0,
        'too_many_fastq': 0,
        'amplicon_mismatch': 0,
        'missing_inputs': 0,
        'other': 0
    }


    for directory in os.listdir():
        logging.info("-----------")

        if directory in ["scripts", "unprocessed_data"]:
            logging.info(f"Skipping: {directory}")
            skipped_count += 1

        elif os.path.isdir(directory):
            logging.info(f"Processing directory: {directory}")
            directory_path = os.path.join(os.getcwd(), directory)
            
            try:
                # Identify the amplicon name from the directory name
                matched_name = identify_amplicon(directory.upper(), amplicon_names)

                # Get inputs for CRISPResso call from amplicon list
                guide_seq, amplicon_seq, orientation, editor, intended_edit, tolerated_edits = CRISPResso_inputs(matched_name)
                
                if not guide_seq or not amplicon_seq or not orientation or not editor:
                    logging.error(f"Missing required inputs for {directory}. Skipping this directory.")
                    error_count += 1
                    error_types['missing_inputs'] += 1
                    continue

                logging.info(f"guide_seq: {guide_seq}, amplicon_seq: {amplicon_seq}, orientation: {orientation}, editor: {editor}")

                
                # Handle PE cases
                if editor == "PE":
                    continue
                    logging.info(f"Running CRISPResso for PE editor on {directory}")
                    success, error_type = run_CRISPResso(guide_seq, amplicon_seq, orientation, editor, directory_path)
                    if success:
                        processed_count += 1
                        logging.info(f"Successfully processed {directory}")
                    else:
                        error_count += 1
                        if error_type in error_types:
                            error_types[error_type] += 1
                # Handle CBE cases
                if editor == "CBE":
                    logging.warning("CBE editor case not yet implemented")
                    skipped_count += 1
                    continue

                # Handle ABE cases
                if editor == "ABE":
                    logging.info(f"Running CRISPResso for ABE editor on {directory}")
                    success, error_type = run_CRISPResso(guide_seq, amplicon_seq, orientation, editor, directory_path)
                    if success:
                        processed_count += 1
                        logging.info(f"Successfully processed {directory}")
                    else:
                        error_count += 1
                        if error_type in error_types:
                            error_types[error_type] += 1

            except ValueError as e:
                logging.error(f"Amplicon identification failed for {directory}: {str(e)}")
                error_count += 1
                error_types['amplicon_mismatch'] += 1
                continue
            except Exception as e:
                logging.error(f"Unexpected error processing {directory}: {str(e)}", exc_info=True)
                error_count += 1
                error_types['other'] += 1
                continue
        else:
            logging.info(f"Skipping non-directory item: {directory}")
            skipped_count += 1
    
    # Summary of processing
    logging.info("="*50)
    logging.info("PIPELINE SUMMARY")
    logging.info(f"Directories processed: {processed_count}")
    logging.info(f"Directories skipped: {skipped_count}")
    logging.info(f"Errors encountered: {error_count}")
    logging.info(f"Log file saved to: {log_file}")
    logging.info("="*50)

    # If errors occurred, log detailed breakdown to file only
    if error_count > 0:
        # Create a logger that only writes to the file (not console)
        file_logger = logging.getLogger('file_only')
        file_logger.setLevel(logging.INFO)
        file_logger.propagate = False  # Don't pass messages to root logger
        # Only add the file handler (not the console handler)
        if not file_logger.handlers:  # Avoid adding duplicate handlers on reruns
            for handler in logging.getLogger().handlers:
                if isinstance(handler, logging.FileHandler):
                    file_logger.addHandler(handler)
                    break

        file_logger.info("")
        file_logger.info("ERROR BREAKDOWN:")
        for err_type, count in error_types.items():
            if count > 0:
                file_logger.info(f"  - {err_type}: {count}")

        file_logger.info("")
        file_logger.info("LIKELY ISSUES:")

        # Find the most common error
        most_common_error = max(error_types, key=error_types.get)
        most_common_count = error_types[most_common_error]

        if most_common_count > 0:
            if most_common_error == 'empty_fastq':
                file_logger.info("  -> Most errors are from empty FASTQ files (0 reads).")
                file_logger.info("     Check that your FASTQ files contain data and are not corrupted.")
            elif most_common_error == 'spaces_in_path':
                file_logger.info("  -> Most errors are from spaces in file paths.")
                file_logger.info("     CRISPResso does not handle spaces in paths well.")
                file_logger.info("     Move your data to a path without spaces or rename directories.")
            elif most_common_error == 'no_fastq':
                file_logger.info("  -> Most errors are from directories with no FASTQ files.")
                file_logger.info("     Ensure each sample directory contains .fastq or .fastq.gz files.")
            elif most_common_error == 'too_many_fastq':
                file_logger.info("  -> Most errors are from directories with more than 2 FASTQ files.")
                file_logger.info("     Each directory should have 1 (single-end) or 2 (paired-end) FASTQ files.")
            elif most_common_error == 'amplicon_mismatch':
                file_logger.info("  -> Most errors are from amplicon name mismatches.")
                file_logger.info("     Check that directory names match amplicon names in your amplicon list file.")
            elif most_common_error == 'missing_inputs':
                file_logger.info("  -> Most errors are from missing inputs in your amplicon list file.")
                file_logger.info("     Ensure all required columns are filled for each amplicon.")
            else:
                file_logger.info("  -> Review the error messages above for details.")

        file_logger.info("="*50)

        # Simple message to terminal
        print(f"\n\033[1;33mSee {log_file} for error breakdown and likely cause.\033[0m")



if __name__ == "__main__":
    main()