# scripts/CRISPResso_inputs.py
import csv

def CRISPResso_inputs(matched_name):
    guide_seq = ""
    amplicon_seq = ""
    orientation = ""
    editor = ""
    intended_edit = None
    tolerated_edits = []
    
    print(f"Looking for CRISPResso inputs for amplicon named: {matched_name}")

    with open("./amplicon_list.csv", newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
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
                break

    return guide_seq, amplicon_seq, orientation, editor, intended_edit, tolerated_edits
