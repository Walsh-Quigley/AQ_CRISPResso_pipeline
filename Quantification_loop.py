import os
import pandas as pd
from scripts.banner import print_banner
from scripts.gather_amplicon_names import gather_amplicon_names
from scripts.identify_amplicon import identify_amplicon
from scripts.CRISPResso_inputs import CRISPResso_inputs
from scripts.generate_search_sequences import generate_search_sequences
from scripts.filter_alleles_file import filter_alleles_file
from scripts.identify_independent_correction import identify_independent_correction
from scripts.handle_missing_directories import add_unanalyzed_directories
from scripts.read_extraction import read_extraction 
from scripts.generate_prism_csv import generate_prism_csv



results = []
def main():
    print_banner()

    #Get a list of all the names in the amplicon_list
    amplicon_names = gather_amplicon_names("amplicon_list.csv")
    print("Amplicon names:", amplicon_names)


    for directory in os.listdir():
        print("-----------")
        if directory in ["scripts", "unprocessed_data", "unadultered_fastqs"]:
            print(f"Skipping: {directory}")

        elif os.path.isdir(directory):
            print(f"Processing directory: {directory}")
            directory_path = os.path.join(os.getcwd(), directory)

            # Identify the amplicon name from the directory name
            matched_name = identify_amplicon(directory.upper(), amplicon_names)

                        # Skip if no valid match found
            if not matched_name:
                print(f"No valid match found in amplicon list; skipping directory: {directory}")
                continue

            # Get inputs for CRISPResso call from amplicon list
            guide_seq, amplicon_seq, orientation, editor, intended_edit, tolerated_edits  = CRISPResso_inputs(matched_name)

            print(f"guide_seq: {guide_seq}, amplicon_seq: {amplicon_seq}, orientation: {orientation}, editor: {editor}, intended_edit: {intended_edit}, tolerated_edits: {tolerated_edits}")


            # Handle PE cases
            if editor == "PE":
                print("Still need to build out PE case")
                continue
            
            # Handle CBE cases
            if editor == "CBE":
                print("Still need to build out CBE case")
                continue

            # Handle ABE cases
            if editor == "ABE":
                search_strings = generate_search_sequences(guide_seq, orientation, editor, intended_edit, tolerated_edits, directory_path)
                correction_with_bystander, correction_without_bystanders = filter_alleles_file(search_strings, directory_path)
                print(f"correction_with_bystander: {correction_with_bystander}, correction_without_bystanders: {correction_without_bystanders}")
                directory_name = os.path.basename(directory_path.rstrip('/'))
                independent_correction = identify_independent_correction(orientation, intended_edit, directory_path)
                if  correction_with_bystander == "NA" or correction_without_bystanders == "NA":
                    print(f"Skipping directory {directory_name} due to missing data.")
                    continue
                if independent_correction == "NA":
                    reads_aligned, reads_total = read_extraction(directory_path)
                    w_bystanders_less_wo_bystanders = correction_with_bystander - correction_without_bystanders
                    results.append({"directory":directory_name,
                                    "reads_aligned": reads_aligned,
                                    "reads_total": reads_total,
                                    "correction_with_bystanders":correction_with_bystander,
                                    "correction_without_bystanders":correction_without_bystanders,
                                    "independent_correction": "NA",
                                    "indep_less_w_bystanders": "NA",
                                    "w_bystanders_less_wo_bystanders": w_bystanders_less_wo_bystanders,
                                    "target_locus":guide_seq,
                                    "perfect_correction":search_strings[0],
                                    "corrected_locus_with_bystanders": ";".join(search_strings)
                                    })
                    continue

                indep_less_w_bystanders = independent_correction - correction_with_bystander
                w_bystanders_less_wo_bystanders = correction_with_bystander - correction_without_bystanders
                print(f"indep_less_w_bystanders: {indep_less_w_bystanders}, w_bystanders_less_wo_bystanders: {w_bystanders_less_wo_bystanders}")
                if indep_less_w_bystanders < 0:
                    print("Warning: Unexpected values detected in correction rates (independent). Please check the input data and calculations.")
                    indep_less_w_bystanders = "error"
                if w_bystanders_less_wo_bystanders < 0:
                    print("Warning: Unexpected values detected in correction rates (Read-based). Please check the input data and calculations.")
                    w_bystanders_less_wo_bystanders = "error"

                reads_aligned, reads_total = read_extraction(directory_path)
                results.append({"directory":directory_name,
                                "reads_aligned": reads_aligned,
                                "reads_total": reads_total,
                                "correction_with_bystanders":correction_with_bystander,
                                "correction_without_bystanders":correction_without_bystanders,
                                "independent_correction": independent_correction,
                                "indep_less_w_bystanders": indep_less_w_bystanders,
                                "w_bystanders_less_wo_bystanders": w_bystanders_less_wo_bystanders,
                                "target_locus":guide_seq,
                                "perfect_correction":search_strings[0],
                                "corrected_locus_with_bystanders": ";".join(search_strings)
                                 })
            
        else:
            print(f"Skipping non-directory item: {directory}")

    df = pd.DataFrame(results, columns=["directory",
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
    
    df = add_unanalyzed_directories(df, skip_dirs=["scripts", "unprocessed_data"], note="Directory not analyzed")

    # sort by directory (alphabetical by default)
    df = df.sort_values(by="directory")

    while True:
        prism_choice = input(f"Would you like to generate a csv output formatted for prism? (y/n)").strip().lower()
        if prism_choice in ("y", "yes"):
            prism_csv_file = os.path.join(os.getcwd(), "prism_formatted_output.csv")    
            prism_df = generate_prism_csv(df)
            prism_df.to_csv(prism_csv_file, index=False)
            print(f"Prism formatted csv saved to: {prism_csv_file}")
            break
        elif prism_choice in ("n", "no"):
            print("Skipping Prism CSV generation")
            break
        else:
            print("Invalid, input. please type y or n")
    
    # save to CSV in the current working directory
    out_file = os.path.join(os.getcwd(), "quantification_summary.csv")
    df.to_csv(out_file, index=False)

    print(f"Saved results to {out_file}")




if __name__ == "__main__":
    main()