# scripts/CRISPResso_inputs.py
import csv

def CRISPResso_inputs(matched_name):
    guide_seq = ""
    amplicon_seq = ""
    orientation = ""
    editor = ""
    
    print(f"Looking for CRISPResso inputs for amplicon named: {matched_name}")
    with open("./amplicon_list.csv", newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["name"] == matched_name:
                guide_seq = row["protospacer_or_PEG"].upper().strip()
                amplicon_seq = row["amplicon"].upper().strip()
                orientation = row["guide_orientation_relative_to_amplicon"].upper().strip()
                editor = row["editor"].upper().strip()
                intended_edit = int(row["intended_edit"].upper().strip())
                tolerated_edits = [int(x) for x in row["tolerated_edits"].split(",") if x.strip().isdigit()]
                break

    return guide_seq, amplicon_seq, orientation, editor, intended_edit, tolerated_edits
