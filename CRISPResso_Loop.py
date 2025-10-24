import os
from scripts.banner import print_banner
from scripts.gather_amplicon_names import gather_amplicon_names
from scripts.identify_amplicon import identify_amplicon
from scripts.CRISPResso_inputs import CRISPResso_inputs
from scripts.run_CRISPResso import run_CRISPResso



def main():
    print_banner()

    #Get a list of all the names in the amplicon_list
    amplicon_names = gather_amplicon_names("amplicon_list.csv")
    print("Amplicon names:", amplicon_names)


    for directory in os.listdir():
        print("-----------")
        if directory in ["scripts", "unprocessed_data"]:
            print(f"Skipping: {directory}")

        elif os.path.isdir(directory):
            print(f"Processing directory: {directory}")
            directory_path = os.path.join(os.getcwd(), directory)

            # Identify the amplicon name from the directory name
            matched_name = identify_amplicon(directory.upper(), amplicon_names)

            # Get inputs for CRISPResso call from amplicon list
            guide_seq, amplicon_seq, orientation, editor, intended_edit, tolerated_edits = CRISPResso_inputs(matched_name)
            print(f"guide_seq: {guide_seq}, amplicon_seq: {amplicon_seq}, orientation: {orientation}, editor: {editor}")

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
                run_CRISPResso(guide_seq, amplicon_seq, orientation, editor, directory_path)

            
        else:
            print(f"Skipping non-directory item: {directory}")
        




if __name__ == "__main__":
    main()