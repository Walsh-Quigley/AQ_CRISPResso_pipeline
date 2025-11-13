# scripts/verify_amplicon_list.py
import csv
import logging
import shutil

def verify_amplicon_list(file):


    required_headers = [
        'name',
        'protospacer_or_PEG',
        'editor',
        'guide_orientation_relative_to_amplicon',
        'amplicon',
        'tolerated_edits',
        'intended_edit'
    ]
    with open(file, newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        
        # Verify headers
        actual_headers = reader.fieldnames
        missing = set(required_headers) - set(actual_headers)
        if missing:
            error_msg = "Header mismatch in amplicon list file:\n"
            error_msg += f"  Missing headers: {', '.join(missing)}\n"
            error_msg += f"  Expected: {', '.join(required_headers)}"
            
            raise ValueError(error_msg)

        
        #verify that amplicons only contain letters
        amplicon_names = []
        non_standard_guides = []
        all_rows = []
        for row_num, row in enumerate(reader, start=2):  # start=2 because row 1 is headers
            # Validate amplicon sequence contains only letters
            amplicon = row['amplicon'].strip()
            if not amplicon.isalpha():
                raise ValueError(
                    f"Invalid amplicon sequence at row {row_num} (name: '{row['name']}'):\n"
                    f"  Amplicon must contain only letters (A-Z, a-z).\n"
                    f"  Found: '{amplicon}'"
                )
                        # Check protospacer length
            protospacer = row['protospacer_or_PEG'].strip()
            if len(protospacer) != 20:
                non_standard_guides.append({
                    'name': row['name'].strip(),
                    'length': len(protospacer),
                    'sequence': protospacer
                })
            #compile the list of amplicon names for the user   
            amplicon_names.append(row['name'].strip().upper())
            all_rows.append(row)
        
        if non_standard_guides:
            logging.warning("="*60)
            logging.warning("Non-standard guide lengths detected")
            logging.warning("="*60)
            for guide in non_standard_guides:
                logging.warning(f"  Name: {guide['name']}")
                logging.warning(f"  Length: {guide['length']} (expected 20)")
                logging.warning(f"  Sequence: {guide['sequence']}")
                logging.warning("")

        #check if the user wants to do anything about the guides that are nonstandard length
        while True:
            response = input("Would you like to truncate these guides to 20bp? (y/n) (probably not if these are for ONE-seq): ").strip().lower()
            
            if response == 'y':
                # Save original file as backup
                backup_file = file.replace('.csv', '_untruncated.csv')
                shutil.copy2(file, backup_file)
                logging.info(f"Original file backed up to: {backup_file}")
                
                # Truncate the protospacer sequences and overwrite the file
                for row in all_rows:
                    protospacer = row['protospacer_or_PEG'].strip()
                    if len(protospacer) > 20:
                        row['protospacer_or_PEG'] = protospacer[:20]
                
                # Write back to the file
                with open(file, 'w', newline='', encoding='utf-8-sig') as f:
                    writer = csv.DictWriter(f, fieldnames=required_headers)
                    writer.writeheader()
                    writer.writerows(all_rows)
                
                logging.info("Truncated non-standard guides to 20bp and updated amplicon_list.csv")
                break
            elif response == 'n':
                logging.info("User declined truncation. Continuing with non-standard guide lengths.")
                break
            else:
                print("Invalid input. Please enter 'y' or 'n'.")
    
    return amplicon_names
