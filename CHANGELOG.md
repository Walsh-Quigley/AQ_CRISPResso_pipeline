# Changelog

All notable changes to this project will be documented in this file.

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
