# scripts/filter_alleles_file.py
import os
import glob
import csv

def filter_alleles_file(search_strings, directory_path):
    crispr_dirs = [d for d in glob.glob(os.path.join(directory_path, "CRISPResso_on_*")) if os.path.isdir(d)]

    if not crispr_dirs:
        print(f"No CRISPResso_on* subdirectory found in {directory_path}")
        return "NA", "NA"
    
    crispresso_subfolder = crispr_dirs[0]  # first matching CRISPResso subdirectory
    allele_files = glob.glob(os.path.join(crispresso_subfolder, "Alleles_frequency_table*.txt"))
        
    if not allele_files:
        print(f"No Alleles_frequency_table file found in {crispresso_subfolder}")
        return "NA", "NA"

    allele_file = allele_files[0]  # first matching Alleles_frequency_table file

    output_file_w = os.path.join(directory_path, "AQ_read_quant_w_bystanders.csv")
    output_file_wo = os.path.join(directory_path, "AQ_read_quant_wo_bystanders.csv")

    filtered_rows_w = []   # rows matching any search string
    filtered_rows_wo = []  # rows matching only search_strings[0]

    with open(allele_file, newline='', encoding='utf-8-sig') as infile:
        reader = csv.reader(infile, delimiter='\t')
        header = next(reader)

        for row in reader:
            if any(s in row[0] for s in search_strings):
                filtered_rows_w.append(row)
            if search_strings[0] in row[0]:
                filtered_rows_wo.append(row)

    # Write "with bystanders" file
    with open(output_file_w, 'w', newline='', encoding='utf-8-sig') as outfile_w:
        writer = csv.writer(outfile_w)
        writer.writerow(header)
        writer.writerows(filtered_rows_w)

    # Write "without bystanders" file
    with open(output_file_wo, 'w', newline='', encoding='utf-8-sig') as outfile_wo:
        writer = csv.writer(outfile_wo)
        writer.writerow(header)
        writer.writerows(filtered_rows_wo)

    print(f"Filtered rows written to {output_file_w} (with bystanders)")
    print(f"Filtered rows written to {output_file_wo} (without bystanders)")

    # Sum last column
    sum_w = sum(float(row[-1]) for row in filtered_rows_w if row[-1].replace('.','',1).isdigit())
    sum_wo = sum(float(row[-1]) for row in filtered_rows_wo if row[-1].replace('.','',1).isdigit())

    return sum_w, sum_wo