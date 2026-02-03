# scripts/verify_allele_results.py
"""
Verification script to compare allele-specific results in quantification_summary.csv
against the AQ_allele_frequency_table_pos*_*.csv files.
"""
import os
import glob
import csv
import pandas as pd

def load_quantification_summary(summary_path):
    """Load the quantification summary CSV."""
    return pd.read_csv(summary_path)

def find_sample_directory(fastqs_dir, sample_name):
    """Find the full directory path for a sample."""
    pattern = os.path.join(fastqs_dir, f"{sample_name}*")
    matches = glob.glob(pattern)
    dirs = [m for m in matches if os.path.isdir(m)]
    return dirs[0] if dirs else None

def load_allele_frequency_table(file_path):
    """Load an AQ_allele_frequency_table CSV and return total reads and sequence->reads dict."""
    if not os.path.exists(file_path):
        return None, None, None

    sequences = {}
    total_reads = 0

    with open(file_path, newline='', encoding='utf-8-sig') as f:
        reader = csv.reader(f)
        header = next(reader)

        for row in reader:
            if len(row) >= 7:
                seq = row[0]
                reads = int(row[6])
                pct = float(row[7]) if len(row) > 7 else 0
                sequences[seq] = {'reads': reads, 'pct': pct}
                total_reads += reads

    return sequences, total_reads, header

def verify_sample(sample_name, summary_row, sample_dir, perfect_correction):
    """Verify allele-specific results for a single sample."""
    results = {
        'sample': sample_name,
        'issues': [],
        'details': {}
    }

    het_position = summary_row.get('het_position')
    het_alleles = summary_row.get('het_alleles')

    if pd.isna(het_position) or pd.isna(het_alleles):
        results['details']['status'] = 'Not heterozygous (no allele-specific data expected)'
        return results

    het_pos = int(het_position)
    allele1, allele2 = het_alleles.split('/')

    # Load allele frequency tables
    file1 = os.path.join(sample_dir, f"AQ_allele_frequency_table_pos{het_pos}_{allele1}.csv")
    file2 = os.path.join(sample_dir, f"AQ_allele_frequency_table_pos{het_pos}_{allele2}.csv")

    seqs1, total1, _ = load_allele_frequency_table(file1)
    seqs2, total2, _ = load_allele_frequency_table(file2)

    if seqs1 is None:
        results['issues'].append(f"Missing file: {os.path.basename(file1)}")
    if seqs2 is None:
        results['issues'].append(f"Missing file: {os.path.basename(file2)}")

    if seqs1 is None or seqs2 is None:
        return results

    results['details']['het_position'] = het_pos
    results['details']['het_alleles'] = het_alleles
    results['details']['allele1_total_reads'] = total1
    results['details']['allele2_total_reads'] = total2

    # Find reads matching perfect correction in each allele
    reads_corrected_allele1 = 0
    reads_corrected_allele2 = 0

    for seq, data in seqs1.items():
        if perfect_correction in seq:
            reads_corrected_allele1 += data['reads']

    for seq, data in seqs2.items():
        if perfect_correction in seq:
            reads_corrected_allele2 += data['reads']

    # Calculate percentages
    pct_corrected_allele1 = (reads_corrected_allele1 / total1 * 100) if total1 > 0 else 0
    pct_corrected_allele2 = (reads_corrected_allele2 / total2 * 100) if total2 > 0 else 0

    results['details']['calculated_correction_allele1'] = pct_corrected_allele1
    results['details']['calculated_correction_allele2'] = pct_corrected_allele2
    results['details']['reads_corrected_allele1'] = reads_corrected_allele1
    results['details']['reads_corrected_allele2'] = reads_corrected_allele2

    # Compare with summary values
    summary_allele1 = summary_row.get('correction_without_bystanders_allele1', 0)
    summary_allele2 = summary_row.get('correction_without_bystanders_allele2', 0)

    if pd.isna(summary_allele1):
        summary_allele1 = 0
    if pd.isna(summary_allele2):
        summary_allele2 = 0

    results['details']['summary_correction_allele1'] = summary_allele1
    results['details']['summary_correction_allele2'] = summary_allele2

    # Check for discrepancies (allow 0.01% tolerance for floating point)
    tolerance = 0.01
    diff1 = abs(pct_corrected_allele1 - summary_allele1)
    diff2 = abs(pct_corrected_allele2 - summary_allele2)

    if diff1 > tolerance:
        results['issues'].append(
            f"Allele1 mismatch: calculated={pct_corrected_allele1:.4f}%, summary={summary_allele1:.4f}%, diff={diff1:.4f}%"
        )
    if diff2 > tolerance:
        results['issues'].append(
            f"Allele2 mismatch: calculated={pct_corrected_allele2:.4f}%, summary={summary_allele2:.4f}%, diff={diff2:.4f}%"
        )

    if not results['issues']:
        results['details']['status'] = 'VERIFIED - values match'

    return results

def verify_all_samples(fastqs_dir):
    """Run verification on all samples in the quantification summary."""
    summary_path = os.path.join(fastqs_dir, "quantification_summary.csv")

    if not os.path.exists(summary_path):
        print(f"Error: quantification_summary.csv not found at {summary_path}")
        return

    summary_df = load_quantification_summary(summary_path)

    print("=" * 80)
    print("ALLELE-SPECIFIC RESULTS VERIFICATION")
    print("=" * 80)

    all_verified = True

    for idx, row in summary_df.iterrows():
        sample_name = row['sample']
        perfect_correction = row.get('perfect_correction', '')

        sample_dir = find_sample_directory(fastqs_dir, sample_name)

        if sample_dir is None:
            print(f"\n{sample_name}: SKIPPED - directory not found")
            continue

        result = verify_sample(sample_name, row, sample_dir, perfect_correction)

        print(f"\n{sample_name}:")
        print("-" * 40)

        if result['issues']:
            all_verified = False
            print("  STATUS: ISSUES FOUND")
            for issue in result['issues']:
                print(f"    - {issue}")
        else:
            status = result['details'].get('status', 'OK')
            print(f"  STATUS: {status}")

        if 'het_position' in result['details']:
            print(f"  Het position: {result['details']['het_position']}")
            print(f"  Het alleles: {result['details']['het_alleles']}")
            print(f"  Allele1 total reads: {result['details'].get('allele1_total_reads', 'N/A')}")
            print(f"  Allele2 total reads: {result['details'].get('allele2_total_reads', 'N/A')}")
            print(f"  Correction (allele1): calculated={result['details'].get('calculated_correction_allele1', 0):.4f}%, "
                  f"summary={result['details'].get('summary_correction_allele1', 0):.4f}%")
            print(f"  Correction (allele2): calculated={result['details'].get('calculated_correction_allele2', 0):.4f}%, "
                  f"summary={result['details'].get('summary_correction_allele2', 0):.4f}%")

    print("\n" + "=" * 80)
    if all_verified:
        print("OVERALL: All samples verified successfully!")
    else:
        print("OVERALL: Some samples have discrepancies - please review above")
    print("=" * 80)

if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        fastqs_dir = sys.argv[1]
    else:
        # Default path
        fastqs_dir = r"c:\Users\senti\Desktop\WorkStuff\Hanbing\fastqs"

    verify_all_samples(fastqs_dir)
