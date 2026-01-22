# scripts/filter_alleles_file.py
import os
import glob
import csv

def find_het_position(crispresso_subfolder):
    quant_files = glob.glob(os.path.join(crispresso_subfolder, "Quantification_window_nucleotide_percentage_table.txt"))
    if not quant_files:
        return None, None, None
    
    with open(quant_files[0], newline='', encoding='utf-8-sig') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)  # positions
        
        # Read percentages for each base
        percentages = {}
        for row in reader:
            if row[0] in ['A', 'C', 'G', 'T']:
                percentages[row[0]] = [float(x) for x in row[1:]]
        
        # Find position where two bases are both between 40-60%
        for pos_idx in range(len(percentages['A'])):
            bases_in_range = []
            for base in ['A', 'C', 'G', 'T']:
                pct = percentages[base][pos_idx]
                if 0.40 <= pct <= 0.60:
                    bases_in_range.append(base)
            
            if len(bases_in_range) == 2:
                # Skip A/G pairs - these are likely editing, not heterozygosity
                if set(bases_in_range) == {'A', 'G'} or set(bases_in_range) == {"C","T"}:
                    continue
                return pos_idx, bases_in_range[0], bases_in_range[1]
    
    return None, None, None  # No het position found



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

def filter_alleles_file_hetero(search_strings, directory_path, orientation, guide_seq):
    crispr_dirs = [d for d in glob.glob(os.path.join(directory_path, "CRISPResso_on_*")) if os.path.isdir(d)]

    if not crispr_dirs:
        print(f"No CRISPResso_on* subdirectory found in {directory_path}")
        return "NA", "NA", "NA", "NA", "NA", "NA", "NA"
    
    crispresso_subfolder = crispr_dirs[0]
    # Find heterozygous position dynamically
    het_pos, base1, base2 = find_het_position(crispresso_subfolder)
    
    if het_pos is None:
        print(f"No heterozygous position found in {crispresso_subfolder}")
        return "NA", "NA", "NA", "NA", "NA", "NA", "NA"
        
    print(f"Found het position: {het_pos}, bases: {base1}/{base2}")

    allele_files = glob.glob(os.path.join(crispresso_subfolder, "Alleles_frequency_table*.txt"))

    if not allele_files:
        print(f"No Alleles_frequency_table file found in {crispresso_subfolder}")
        return "NA", "NA", "NA", "NA", "NA", "NA", "NA"

    allele_file = allele_files[0]

    filtered_rows_w_base1 = []
    filtered_rows_w_base2 = []
    filtered_rows_wo_base1 = []
    filtered_rows_wo_base2 = []     

    with open(allele_file, newline='', encoding='utf-8-sig') as infile:
        reader = csv.reader(infile, delimiter='\t')
        header = next(reader)

        total_reads_base1 = 0
        total_reads_base2 = 0

        for row in reader:
            sequence = row[0]
            base_at_pos = sequence[het_pos]
                
            read_count = int(row[6])
            if base_at_pos == base1:
                total_reads_base1 += read_count
            elif base_at_pos == base2:
                total_reads_base2 += read_count

            # Check with bystanders (any search string)
            if any(s in sequence for s in search_strings):
                if base_at_pos == base1:
                    filtered_rows_w_base1.append(row)
                elif base_at_pos == base2:
                    filtered_rows_w_base2.append(row)
            
            # Check without bystanders (only first search string)
            if search_strings[0] in sequence:
                if base_at_pos == base1:
                    filtered_rows_wo_base1.append(row)
                elif base_at_pos == base2:
                    filtered_rows_wo_base2.append(row)

    reads_w_base1 = sum(int(row[6]) for row in filtered_rows_w_base1)
    reads_w_base2 = sum(int(row[6]) for row in filtered_rows_w_base2)
    reads_wo_base1 = sum(int(row[6]) for row in filtered_rows_wo_base1)
    reads_wo_base2 = sum(int(row[6]) for row in filtered_rows_wo_base2)

    pct_w_base1 = (reads_w_base1 / total_reads_base1 * 100) if total_reads_base1 > 0 else 0
    pct_w_base2 = (reads_w_base2 / total_reads_base2 * 100) if total_reads_base2 > 0 else 0
    pct_wo_base1 = (reads_wo_base1 / total_reads_base1 * 100) if total_reads_base1 > 0 else 0
    pct_wo_base2 = (reads_wo_base2 / total_reads_base2 * 100) if total_reads_base2 > 0 else 0

    return pct_w_base1, pct_w_base2, pct_wo_base1, pct_wo_base2, het_pos + 1, base1, base2
