# scripts/gather_amplicon_names.py
import csv

def gather_amplicon_names(file):
    names = []
    with open(file, newline='', encoding ='utf-8-sig') as f:
        reader = csv.DictReader(f)
        
        for row in reader:
            names.append(row['name'].strip().upper())
    return names
