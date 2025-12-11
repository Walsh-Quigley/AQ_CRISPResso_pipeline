# scripts/CRISPResso_inputs.py
import csv
import os
from scripts.verify_amplicon_list import find_amplicon_list_file

def CRISPResso_inputs(matched_name):
    """
    Extract CRISPResso parameters for a given amplicon name.
    
    Args:
        matched_name: Name of the amplicon to look up
        
    Returns:
        tuple: (guide_seq, amplicon_seq, orientation, editor, intended_edit, tolerated_edits)
    """

    guide_seq = ""
    amplicon_seq = ""
    orientation = ""
    editor = ""
    intended_edit = None
    tolerated_edits = []

    # Check if matched_name is None before proceeding
    if matched_name is None:
        print("ERROR: No amplicon name provided (matched_name is None)")
        print("This usually means the directory name didn't match any amplicon in the amplicon list file")
        return guide_seq, amplicon_seq, orientation, editor, intended_edit, tolerated_edits
    
    
    print(f"Looking for CRISPResso inputs for amplicon named: {matched_name}")


    try:
        # Search for amplicon list file in parent directory
        amplicon_file = find_amplicon_list_file("..")
        amplicon_file_path = os.path.join("..", amplicon_file)

        with open(amplicon_file_path, newline='', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            found = False

            for row in reader:
                if row["name"].strip().upper() == matched_name.strip().upper():
                    guide_seq = row["protospacer_or_PEG"].upper().strip()

                    amplicon_seq = row["amplicon"].upper().strip()
                    orientation = row["guide_orientation_relative_to_amplicon"].upper().strip()
                    editor = row["editor"].upper().strip()
                    
                    # Handle intended_edit column flexibly for ONE-seq
                    intended_edit_raw = row["intended_edit"].strip().upper().replace("-", "")
                    if intended_edit_raw == "ONESEQ":
                        intended_edit = "ONESEQ"
                        print("Intended edit set to ONESEQ")
                    elif intended_edit_raw.isdigit():
                        intended_edit = int(intended_edit_raw)
                    else:
                        print(f"Warning: intended_edit value '{row['intended_edit']}' is not recognized. Setting to None.")
                    

                    tolerated_edits = [int(x) for x in row["tolerated_edits"].split(",") if x.strip().isdigit()]
                    found = True
                    break
            if not found:
                print(f"ERROR: Amplicon name '{matched_name}' not found in amplicon list file")

    except FileNotFoundError as e:
        print(f"ERROR: Amplicon list file error: {e}")
    except KeyError as e:
        print(f"ERROR: Missing expected column in amplicon list file: {e}")
    except Exception as e:
        print(f"ERROR: Unexpected error reading amplicon list file: {e}")
     
    return guide_seq, amplicon_seq, orientation, editor, intended_edit, tolerated_edits
