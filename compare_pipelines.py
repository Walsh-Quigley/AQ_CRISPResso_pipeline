"""
Compare ABE summary outputs from the old and new pipelines.

Usage:
    python compare_pipelines.py

Inner-joins on sample name, so samples missing from either side are skipped
(only those that succeeded in BOTH pipelines are compared). Numeric columns
flagged when |old - new| exceeds TOLERANCE.
"""

import pandas as pd

OLD_CSV = "FileToCompare/quantification_summary.csv"
NEW_CSV = "FileToCompare/ABE_Quantification_Summary.csv"
OUT_FULL = "comparison_full.csv"
OUT_DIFFS = "comparison_diffs.csv"
TOLERANCE = 0.5  # percentage points

# old name -> new name
OLD_TO_NEW = {
    "correction_with_bystanders":      "correction_with_tolerated_bystanders",
    "column E minus column D":          "w_bystanders_minus_wo_bystanders",
    "column F minus column E":          "any_AtoG_minus_w_bystanders",
    "column G minus column F":          "any_change_minus_any_AtoG",
    # het-only renames
    "correction_without_bystanders_allele1": "correction_wo_bystanders_allele1",
    "correction_with_bystanders_allele1":    "correction_w_bystanders_allele1",
    "column L minus column K":               "w_bystanders_minus_wo_bystanders_allele1",
    "column M minus column L":               "any_AtoG_minus_w_bystanders_allele1",
    "column N minus column M":               "any_change_minus_any_AtoG_allele1",
    "correction_without_bystanders_allele2": "correction_wo_bystanders_allele2",
    "correction_with_bystanders_allele2":    "correction_w_bystanders_allele2",
    "column S minus column R":               "w_bystanders_minus_wo_bystanders_allele2",
    "column T minus column S":               "any_AtoG_minus_w_bystanders_allele2",
    "column U minus column T":               "any_change_minus_any_AtoG_allele2",
}

# numeric columns to compare (using NEW names)
NUMERIC_COLS = [
    "reads_total",
    "reads_aligned",
    "correction_without_bystanders",
    "correction_with_tolerated_bystanders",
    "correction_with_any_AtoG_change",
    "correction_with_any_change_in_protospacer",
    "w_bystanders_minus_wo_bystanders",
    "any_AtoG_minus_w_bystanders",
    "any_change_minus_any_AtoG",
    # het-only
    "correction_wo_bystanders_allele1",
    "correction_w_bystanders_allele1",
    "correction_with_any_AtoG_change_allele1",
    "correction_with_any_change_in_protospacer_allele1",
    "w_bystanders_minus_wo_bystanders_allele1",
    "any_AtoG_minus_w_bystanders_allele1",
    "any_change_minus_any_AtoG_allele1",
    "correction_wo_bystanders_allele2",
    "correction_w_bystanders_allele2",
    "correction_with_any_AtoG_change_allele2",
    "correction_with_any_change_in_protospacer_allele2",
    "w_bystanders_minus_wo_bystanders_allele2",
    "any_AtoG_minus_w_bystanders_allele2",
    "any_change_minus_any_AtoG_allele2",
]


def main():
    old = pd.read_csv(OLD_CSV).rename(columns=OLD_TO_NEW)
    new = pd.read_csv(NEW_CSV)

    old_samples = set(old["sample"])
    new_samples = set(new["sample"])

    print(f"Old pipeline samples: {len(old_samples)}")
    print(f"New pipeline samples: {len(new_samples)}")
    print(f"Samples in BOTH:      {len(old_samples & new_samples)}")
    print(f"Only in old:          {len(old_samples - new_samples)}")
    print(f"Only in new:          {len(new_samples - old_samples)}")
    print()

    only_old = sorted(old_samples - new_samples)
    only_new = sorted(new_samples - old_samples)
    if only_old:
        print(f"Samples only in old ({len(only_old)}):")
        for s in only_old[:20]:
            print(f"  {s}")
        if len(only_old) > 20:
            print(f"  ...and {len(only_old) - 20} more")
        print()
    if only_new:
        print(f"Samples only in new ({len(only_new)}):")
        for s in only_new[:20]:
            print(f"  {s}")
        if len(only_new) > 20:
            print(f"  ...and {len(only_new) - 20} more")
        print()

    merged = old.merge(new, on="sample", suffixes=("_old", "_new"), how="inner")
    print(f"Comparing {len(merged)} samples...")
    print()

    diff_rows = []  # long-format records of mismatches
    delta_cols = []  # column names of computed deltas (added to merged for full CSV)
    any_mismatches = False

    for col in NUMERIC_COLS:
        old_col, new_col = f"{col}_old", f"{col}_new"
        if old_col not in merged.columns or new_col not in merged.columns:
            print(f"[SKIP] {col}: not in both files")
            continue

        # Coerce to numeric (old has 'NA' strings in some columns)
        merged[old_col] = pd.to_numeric(merged[old_col], errors="coerce")
        merged[new_col] = pd.to_numeric(merged[new_col], errors="coerce")

        delta_col = f"{col}_delta"
        merged[delta_col] = merged[new_col] - merged[old_col]
        delta_cols.append(delta_col)

        # only flag where BOTH values are real numbers
        comparable = merged[merged[old_col].notna() & merged[new_col].notna()]
        mismatches = comparable[comparable[delta_col].abs() > TOLERANCE]

        if mismatches.empty:
            print(f"[OK]   {col}")
        else:
            any_mismatches = True
            print(f"[DIFF] {col} - {len(mismatches)} sample(s) exceed tolerance")
            # Add to long-format diff records
            for _, row in mismatches.iterrows():
                diff_rows.append({
                    "sample": row["sample"],
                    "metric": col,
                    "old_value": row[old_col],
                    "new_value": row[new_col],
                    "delta": row[delta_col],
                })

    if not any_mismatches:
        print("\nAll numeric columns within tolerance.")

    # Write full merged CSV (sample + all old/new/delta columns)
    merged.to_csv(OUT_FULL, index=False)
    print(f"\nFull comparison saved to: {OUT_FULL}")

    # Write diffs-only CSV (long format)
    if diff_rows:
        pd.DataFrame(diff_rows).to_csv(OUT_DIFFS, index=False)
        print(f"Diffs-only saved to:      {OUT_DIFFS}  ({len(diff_rows)} mismatches)")
    else:
        print(f"No diffs to save.")


if __name__ == "__main__":
    main()
