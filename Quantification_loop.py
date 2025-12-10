import os
import pandas as pd
import re
import logging
from scripts.banner import print_banner
from scripts.verify_amplicon_list import verify_amplicon_list
from scripts.identify_amplicon import identify_amplicon
from scripts.CRISPResso_inputs import CRISPResso_inputs
from scripts.handle_missing_directories import add_unanalyzed_directories
from scripts.generate_prism_csv import generate_prism_csv
from scripts.logging_setup import setup_logging
from scripts.process_oneseq import process_oneseq
from scripts.process_ABE_case import process_ABE_case
from scripts.yes_no import yes_no




def main():

    results = []
    one_seq_results = []
    
    #set up logging
    log_file = setup_logging()

    print_banner()

    #Get a list of all the names in the amplicon_list
    amplicon_names = verify_amplicon_list("amplicon_list.csv")
    logging.info(f"Amplicon names: {amplicon_names}")

    # Move into the /fastqs directory
    fastqs_dir = os.path.join(os.getcwd(), "fastqs")
    os.chdir(fastqs_dir)

    # Track processing statistics
    processed_count = 0
    skipped_count = 0
    error_count = 0

    # Track specific error types
    error_types = {
        'missing_crispresso_output': 0,
        'missing_correction_data': 0,
        'amplicon_mismatch': 0,
        'file_not_found': 0,
        'other': 0
    }



    for directory in os.listdir():
        logging.info("-----------")
        if directory in ["scripts", "unprocessed_data", "unadultered_fastqs"]:
            logging.info(f"Skipping: {directory}")
            skipped_count += 1

        elif os.path.isdir(directory):
            logging.info(f"Processing directory: {directory}")
            directory_path = os.path.join(os.getcwd(), directory)

            try:
                # Identify the amplicon name from the directory name
                matched_name = identify_amplicon(directory.upper(), amplicon_names)

                # Get inputs for CRISPResso call from amplicon list
                guide_seq, amplicon_seq, orientation, editor, intended_edit, tolerated_edits  = CRISPResso_inputs(matched_name)

                if intended_edit == "ONESEQ":
                    result = process_oneseq(directory_path, guide_seq, orientation)
                    one_seq_results.append(result)
                    processed_count += 1
                    continue

                # Handle PE cases
                if editor == "PE":
                    logging.warning("Still need to build out PE case")
                    skipped_count += 1
                    continue
                
                # Handle CBE cases
                if editor == "CBE":
                    logging.warning("Still need to build out CBE case")
                    skipped_count += 1
                    continue

                # Handle ABE cases
                if editor == "ABE":
                    result = process_ABE_case(directory_path, guide_seq, orientation, editor, intended_edit, tolerated_edits)
                    results.append(result)
                    processed_count += 1
                    continue


            except ValueError as e:
                error_str = str(e)
                logging.error(f"Skipping: {directory}: {e}")
                error_count += 1
                # Categorize the error
                if "Missing correction data" in error_str:
                    error_types['missing_correction_data'] += 1
                elif "No amplicon" in error_str or "not found in amplicon" in error_str:
                    error_types['amplicon_mismatch'] += 1
                else:
                    error_types['other'] += 1
                continue
            except FileNotFoundError as e:
                logging.error(f"File not found for {directory}: {e}")
                error_count += 1
                if "CRISPResso" in str(e):
                    error_types['missing_crispresso_output'] += 1
                else:
                    error_types['file_not_found'] += 1
                continue
            except Exception as e:
                logging.error(f"Unexpected error processing {directory}: {e}")
                error_count += 1
                error_types['other'] += 1
                continue

        else:
            logging.warning(f"Skipping non-directory item: {directory}")

    #Saving the ONE-seq results if any
    if one_seq_results:
        one_seq_csv_file = os.path.join(os.getcwd(), "quantification_summary_ONE-seq.csv")   
        df_one_seq = pd.DataFrame(one_seq_results, columns=[
            "sample",
            "reads_aligned",
            "reads_total",
            "Percent_of_reads_with_A>G_in_first_10bp",
            "Percent_of_reads_with_A>G_in_protospacer",
            "guide_seq",
            "A>G_10bp_search_sequences",
            "A>G_any_search_sequences"
            ])
        df_one_seq.to_csv(one_seq_csv_file, index=False)
        logging.info(f"Saved ONE-seq summary to {one_seq_csv_file}")
    else:
        logging.info("No ONE-seq results to save.")


    #Saving non ONE-seq results if any
    if results:
        df = pd.DataFrame(results, columns=["sample",
                                            "reads_aligned",
                                            "reads_total",
                                            "correction_with_bystanders",
                                            "correction_without_bystanders",
                                            "independent_correction",
                                            "indep_less_w_bystanders",
                                            "w_bystanders_less_wo_bystanders",
                                            "target_locus",
                                            "perfect_correction",
                                            "corrected_locus_with_bystanders"])
        
        #adding unanalyzed directories to the main results dataframe
        df = add_unanalyzed_directories(
            df, 
            skip_dirs=["scripts", "unprocessed_data"], 
            convert_to_sample_name=True,
            note="Directory not analyzed"
        )

        if one_seq_results:
            oneseq_directories = [result['sample'] for result in one_seq_results]
            df.loc[df['sample'].isin(oneseq_directories), 'note'] = "Analyzed in quantification_loop_ONE-seq.csv"

        # sort by directory (alphabetical by default)
        df = df.sort_values(by="sample")

        # Prompt user to generate Prism formatted CSV
        if yes_no("Would you like to generate a csv output formatted for prism?"):
            logging.info("="*50)
            logging.info("Performing Prism csv generation")
            prism_csv_file = os.path.join(os.getcwd(), "prism_formatted_output.csv")   
            prism_df = generate_prism_csv(df)
            prism_df.to_csv(prism_csv_file, index=False)
            logging.info(f"Prism formatted csv saved to: {prism_csv_file}")
        else:
            logging.info("Skipping Prism CSV generation")

        
        # save to CSV in the current working directory
        out_file = os.path.join(os.getcwd(), "quantification_summary.csv")
        df.to_csv(out_file, index=False)

        logging.info(f"Saved results to {out_file}")

    # At the very end of main(), before the function closes:
    logging.info("="*50)
    logging.info("PIPELINE SUMMARY")
    logging.info(f"Directories processed: {processed_count}")
    logging.info(f"Directories skipped: {skipped_count}")
    logging.info(f"Errors encountered: {error_count}")
    logging.info(f"ONE-seq results: {len(one_seq_results)}")
    logging.info(f"Standard results: {len(results)}")
    logging.info(f"Log file: {log_file}")
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
            if most_common_error == 'missing_crispresso_output':
                file_logger.info("  -> Most errors are from missing CRISPResso output folders.")
                file_logger.info("     Run CRISPResso_Loop.py first to generate the required output files.")
            elif most_common_error == 'missing_correction_data':
                file_logger.info("  -> Most errors are from missing correction data in CRISPResso output.")
                file_logger.info("     Check that CRISPResso ran successfully and generated all output files.")
            elif most_common_error == 'amplicon_mismatch':
                file_logger.info("  -> Most errors are from amplicon name mismatches.")
                file_logger.info("     Check that directory names match amplicon names in amplicon_list.csv.")
            elif most_common_error == 'file_not_found':
                file_logger.info("  -> Most errors are from missing files.")
                file_logger.info("     Check that all required input files exist.")
            else:
                file_logger.info("  -> Review the error messages above for details.")

        file_logger.info("="*50)

        # Simple message to terminal
        print(f"\n\033[1;33mSee {log_file} for error breakdown and likely cause.\033[0m")




if __name__ == "__main__":
    main()