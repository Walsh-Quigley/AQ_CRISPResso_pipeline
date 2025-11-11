# scripts/CRISPResso_inputs.py
import csv

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
        print("This usually means the directory name didn't match any amplicon in the amplicon_list.csv")
        return guide_seq, amplicon_seq, orientation, editor, intended_edit, tolerated_edits
    
    
    print(f"Looking for CRISPResso inputs for amplicon named: {matched_name}")


    try:
        with open("../amplicon_list.csv", newline='', encoding='utf-8-sig') as f:
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
                print(f"ERROR: Amplicon name '{matched_name}' not found in amplicon_list.csv")
    
    except FileNotFoundError:
        print("ERROR: amplicon_list.csv file not found in current directory")
    except KeyError as e:
        print(f"ERROR: Missing expected column in amplicon_list.csv: {e}")
    except Exception as e:
        print(f"ERROR: Unexpected error reading amplicon_list.csv: {e}")
     
    return guide_seq, amplicon_seq, orientation, editor, intended_edit, tolerated_edits
