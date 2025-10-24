# scripts/amplicon_names.py
import csv

def identify_amplicon_names_list(file):
    names = []
    with open(file, newline='', encoding ='utf-8-sig') as f:
        reader = csv.DictReader(f)
        
        for row in reader:
            names.append(row['name'].strip().upper())
    return names
