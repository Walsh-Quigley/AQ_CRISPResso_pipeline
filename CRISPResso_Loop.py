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
        amplicon_names = verify_amplicon_list("amplicon_list.csv")
        logging.info(f"Amplicon names: {amplicon_names}")
    except ValueError as e:
        logging.error(f"Amplicon list verification failed: {e}")
        sys.exit(1)
    except FileNotFoundError:
        logging.error("amplicon_list.csv not found")
        sys.exit(1)

    # Move into the /fastqs directory
    fastqs_dir = os.path.join(os.getcwd(), "fastqs")
    os.chdir(fastqs_dir)
    logging.info(f"Changed to directory: {fastqs_dir}")

    # Track processing statistics
    processed_count = 0
    skipped_count = 0
    error_count = 0


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
                    continue

                logging.info(f"guide_seq: {guide_seq}, amplicon_seq: {amplicon_seq}, orientation: {orientation}, editor: {editor}")

                
                # Handle PE cases
                if editor == "PE":
                    logging.warning("PE editor case not yet implemented")
                    skipped_count += 1
                    continue
                
                # Handle CBE cases
                if editor == "CBE":
                    logging.warning("CBE editor case not yet implemented")
                    skipped_count += 1
                    continue

                # Handle ABE cases
                if editor == "ABE":
                    logging.info(f"Running CRISPResso for ABE editor on {directory}")
                    run_CRISPResso(guide_seq, amplicon_seq, orientation, editor, directory_path)
                    processed_count += 1
                    logging.info(f"Successfully processed {directory}")

            except ValueError as e:
                logging.error(f"Amplicon identification failed for {directory}: {str(e)}")
                error_count += 1
                continue
            except Exception as e:
                logging.error(f"Unexpected error processing {directory}: {str(e)}", exc_info=True)
                error_count += 1
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



if __name__ == "__main__":
    main()