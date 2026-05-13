# Changelog

All notable changes to this project will be documented in this file.

## [1.0.3] - 2026-01-28

### Added

- **Prism CSV generation for heterozygous samples**:
  - New function `generate_prism_csv_het()` in `scripts/generate_prism_csv.py`
  - Outputs separate files: `prism_formatted_output_het_allele1.csv` and `prism_formatted_output_het_allele2.csv`
  - Allele-specific files contain the same column structure as non-het prism output

- **Allele read counts** in heterozygous output:
  - New columns `reads_aligned_allele1` and `reads_aligned_allele2` in quantification summary
  - `total_A_to_G_hetero()` now returns read counts per allele

- **Test suite** for core functionality:
  - `tests/test_het_cases.py` - heterozygous detection and processing tests
  - `tests/test_identify_independent_correction.py` - independent correction calculation tests
  - `tests/verify_allele_results.py` - allele filtering verification

- **README documentation** for heterozygous output format and column descriptions

### Changed

- **Re-enabled Prism CSV generation** in `Quantification_loop.py`:
  - Now splits output by het vs non-het samples
  - Het samples generate allele-specific Prism files

- **Column naming in Prism output** - standardized to "minus" convention:
  - `w_bystanders_less_wo_bystanders` → `w_bystanders_minus_wo_bystanders`
  - `indep_less_w_bystanders` → `any_AtoG_minus_w_bystanders`
  - Added `any_protospacer_change_minus_any_AtoG` column

- **Cleaned up** commented-out code blocks in `scripts/process_ABE_case.py`

## [1.0.2] - 2026-01-26

### Changed

- **Refactored het detection logic** in `scripts/process_ABE_case.py`:
  - Het position now detected once via `find_het_position()` at the start
  - Het-specific functions only called when heterozygosity is detected
  - Non-het samples no longer trigger redundant het function calls

- **Dynamic output columns** in `Quantification_loop.py`:
  - Het-specific columns (K-Z) only included if ANY sample has het data
  - Non-het datasets now produce cleaner output without NA-filled het columns

- **Variable naming consistency**:
  - `total_A_to_G_value` → `pct_A_to_G_value` (reflects that it's a percentage)
  - `total_A_to_G_allele1/2` → `pct_A_to_G_allele1/2`
  - `het_pos` vs `het_pos_display` distinction (0-indexed internal vs 1-indexed output)

### Fixed

- Typo: `corection` → `correction` in variable names

## [1.0.1] - 2026-01-22

### Added

- **Heterozygous sample support** - New functionality to analyze samples with heterozygous positions
  - New function `find_het_position()` in `scripts/filter_alleles_file.py` that automatically detects heterozygous positions by finding bases at 40-60% frequency (excludes A/G and C/T pairs to avoid confusing editing with heterozygosity)
  - New function `filter_alleles_file_hetero()` to calculate correction metrics separately for each allele
  - New output columns in `quantification_summary.csv`:
    - `correction_w_bystanders_allele1` / `correction_w_bystanders_allele2`
    - `correction_wo_bystanders_allele1` / `correction_wo_bystanders_allele2`
    - `het_position` - The detected heterozygous position (1-indexed)
    - `het_alleles` - The two bases found at the het position (e.g., "A/C")

### Changed

- **Renamed metrics for clarity**:
  - `independent_correction` → `correction_with_any_change_in_protospacer`
  - `non_tolerated_A_to_G` → `total_A_to_G_changes`
  - Function `non_tol_A_to_G()` renamed to `total_A_to_G()` in `scripts/identify_independent_correction.py`
- Simplified `total_A_to_G()` function - removed tolerated edits filtering from A-to-G counting logic

## [1.0.0] - 2026-01-20

### Added

- **New metric: `non_tolerated_A_to_G`** - Percentage of reads containing the intended edit AND at least one non-tolerated A-to-G bystander edit
  - New function `non_tol_A_to_G()` in `scripts/identify_independent_correction.py`
  - New column in `quantification_summary.csv` output
  - Properly handles both Forward (F) and Reverse (R) guide orientations
  - Excludes tolerated bystander positions from the count

### Changed

- `scripts/process_ABE_case.py` now calls `non_tol_A_to_G()` instead of `identify_independent_correction()` for read-based correction metrics
- `Quantification_loop.py` DataFrame now includes `non_tolerated_A_to_G` column

### Technical Details

The `non_tol_A_to_G()` function:
- Returns a tuple: `(independent_correction_percentage, non_tolerated_A_to_G_percentage)`
- For Forward orientation: detects A→G changes
- For Reverse orientation: reverse complements the guide sequence and detects T→C changes (equivalent to A→G on the edited strand)
- Uses 0-indexed positions internally while accepting 1-indexed input positions
- Reads from `Alleles_frequency_table_around_*.txt` in CRISPResso output

### Movement to Distribution
The updated loop will be moved to distribution after perfmormace testing with hardcoded data to verify calculation logic.
