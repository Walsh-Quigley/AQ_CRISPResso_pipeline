"""
Tests the filter_alleles_file.py functions.

Tests:
    - find_het_position: Detecting heterozygous positions from base percentages
    - filter_alleles_file: Filtering alleles for correction quantification
    - filter_alleles_file_hetero: Allele-specific filtering for het samples
    - write_hetero_frequency_tables: Writing split allele frequency files
"""

import os
import pytest

from scripts.filter_alleles_file import (
    find_het_position,
    filter_alleles_file,
    filter_alleles_file_hetero,
    write_hetero_frequency_tables
)

class TestFindHetPosition:
    """Tests for find_het_position() function."""

    def test_detects_AC_het_at_50_50(self, crispresso_output_dir, file_factory):
        """Test: A/C heterozygous detection at 50%/50% split

        Expected return: (index=2, 'A', 'C')

        Reasoning: Two bases in 40-60% range at same position, assume this means heterozygous.
        Test to make sure we grab the first instance of this, and that is correct
        """
        temp_dir, crispresso_subdir = crispresso_output_dir

        #pos 3 (index 2) has 50%-A 50%-C
        base_pcts ={
            'A': [0.95, 0.02, 0.50, 0.02, 0.95],
            'C': [0.02, 0.95, 0.50, 0.95, 0.02],
            'G': [0.02, 0.02, 0.00, 0.02, 0.02],
            'T': [0.01, 0.01, 0.00, 0.01, 0.01],
        }
        file_factory.create_quant_window_file(crispresso_subdir, base_pcts)

        het_pos, base1, base2 = find_het_position(crispresso_subdir)

        assert het_pos == 2 #0 index
        assert set([base1, base2]) == {'A', 'C'}
        
    def test_detects_het_at_boundary_40_60(self, crispresso_output_dir, file_factory):
        """
        Test: Het detection at the boundary of the 40-60% range.
        
        Expected return: (index=0, 'A', 'C')
        
        Reasoning: The function uses 0.40 <= pct <= 0.60 as the range. This tests
        that the boundary values (exactly 40% and 60%) are correctly included.
        Edge cases at boundaries reveal off-by-one errors.
        """
        temp_dir, crispresso_subdir = crispresso_output_dir
        
        # Position 1 (index 0) has A=40%, C=60% - right at boundaries
        base_pcts = {
            'A': [0.40, 0.95, 0.02],
            'C': [0.60, 0.02, 0.95],
            'G': [0.00, 0.02, 0.02],
            'T': [0.00, 0.01, 0.01],
        }
        file_factory.create_quant_window_file(crispresso_subdir, base_pcts)
        
        het_pos, base1, base2 = find_het_position(crispresso_subdir)
        
        assert het_pos == 0 #0 index
        assert set([base1, base2]) == {'A', 'C'}

    def test_no_het_when_one_base_dominant(self, crispresso_output_dir, file_factory):
        """
        Test: No heterozygous position when each position has one dominant base.
        
        Expected return: (None, None, None)
        
        Reasoning: If no position has two bases both in 40-60% range,
        the sample is homozygous and should return None for all values.
        """
        temp_dir, crispresso_subdir = crispresso_output_dir
        
        # All positions have one dominant base (>90%)
        base_pcts = {
            'A': [0.95, 0.02, 0.02],
            'C': [0.02, 0.95, 0.02],
            'G': [0.02, 0.02, 0.95],
            'T': [0.01, 0.01, 0.01],
        }
        file_factory.create_quant_window_file(crispresso_subdir, base_pcts)
        
        het_pos, base1, base2 = find_het_position(crispresso_subdir)
        
        assert het_pos is None
        assert base1 is None
        assert base2 is None

    def test_excludes_AG_pairs_as_editing_artifact(self, crispresso_output_dir, file_factory):
        """
        Test: A/G at 50%/50% should NOT be detected as heterozygous.
        
        Expected return: (None, None, None)
        
        Reasoning: A/G pairs are excluded because they likely represent 
        ABE editing (A->G conversion), not true heterozygosity.
        """
        temp_dir, crispresso_subdir = crispresso_output_dir
        
        # Position 1 has A=50%, G=50% - should be excluded
        base_pcts = {
            'A': [0.50, 0.95, 0.02],
            'C': [0.00, 0.02, 0.95],
            'G': [0.50, 0.02, 0.02],
            'T': [0.00, 0.01, 0.01],
        }
        file_factory.create_quant_window_file(crispresso_subdir, base_pcts)
        
        het_pos, base1, base2 = find_het_position(crispresso_subdir)
        
        assert het_pos is None

    def test_excludes_CT_pairs_as_editing_artifact(self, crispresso_output_dir, file_factory):
        """
        Test: C/T at 50%/50% should NOT be detected as heterozygous.
        
        Expected return: (None, None, None)
        
        Reasoning: C/T pairs are excluded because they likely represent
        CBE editing (C->T conversion), not true heterozygosity.
        """
        temp_dir, crispresso_subdir = crispresso_output_dir
        
        # Position 1 has C=50%, T=50% - should be excluded
        base_pcts = {
            'A': [0.00, 0.95, 0.02],
            'C': [0.50, 0.02, 0.95],
            'G': [0.00, 0.02, 0.02],
            'T': [0.50, 0.01, 0.01],
        }
        file_factory.create_quant_window_file(crispresso_subdir, base_pcts)
        
        het_pos, base1, base2 = find_het_position(crispresso_subdir)
        
        assert het_pos is None

class TestFilterAllelesFile:
    """Tests for filter_alleles_file() function."""

    def test_calculates_correction_percentages(self, crispresso_output_dir, file_factory):
        """
        Test: Correctly sums percentages for matching sequences.
        
        Expected return: (sum_w_bystanders, sum_wo_bystanders)
        
        Reasoning: The function searches for corrected sequences in the allele table.
        'with bystanders' matches ANY search string, 'without bystanders' matches 
        only the first (perfect correction). Percentages should sum correctly.
        """
        temp_dir, crispresso_subdir = crispresso_output_dir
        
        # Create alleles file with known sequences and percentages
        rows = [
            ("ACGTGCGTACGT", 50),   # Perfect correction (search_strings[0])
            ("GCGTGCGTACGT", 30),   # With bystander edit (search_strings[1])
            ("ACGTACGTACGT", 20),   # No correction (original)
        ]
        file_factory.create_alleles_frequency_table(crispresso_subdir, rows)
        
        search_strings = ["ACGTGCGTACGT", "GCGTGCGTACGT"]
        
        sum_w, sum_wo = filter_alleles_file(search_strings, temp_dir)
        
        # With bystanders: 50% + 30% = 80%
        # Without bystanders: 50% only
        assert sum_w == 80.0
        assert sum_wo == 50.0

    def test_returns_NA_when_no_crispresso_dir(self, temp_dir):
        """
        Test: Returns 'NA' when no CRISPResso directory exists.
        
        Expected return: ("NA", "NA")
        
        Reasoning: If CRISPResso hasn't run or directory is missing,
        function should gracefully return NA rather than crash.
        """
        search_strings = ["ACGTGCGT"]
        
        sum_w, sum_wo = filter_alleles_file(search_strings, temp_dir)
        
        assert sum_w == "NA"
        assert sum_wo == "NA"

    def test_returns_NA_when_no_alleles_file(self, crispresso_output_dir):
        """
        Test: Returns 'NA' when alleles file doesn't exist.
        
        Expected return: ("NA", "NA")
        
        Reasoning: CRISPResso directory exists but alleles file is missing.
        Should return NA gracefully.
        """
        temp_dir, crispresso_subdir = crispresso_output_dir
        # Don't create any files
        
        search_strings = ["ACGTGCGT"]
        
        sum_w, sum_wo = filter_alleles_file(search_strings, temp_dir)
        
        assert sum_w == "NA"
        assert sum_wo == "NA"

class TestFilterAllelesFileHetero:
    """Tests for filter_alleles_file_hetero() function."""

    def test_splits_reads_by_allele_and_calculates_percentages(self, crispresso_output_dir, file_factory):
        """
        Test: Correctly splits reads by het allele and calculates per-allele percentages.
        
        Expected return: (pct_w_base1, pct_w_base2, pct_wo_base1, pct_wo_base2, het_pos, base1, base2)
        
        Reasoning: For het samples, we need separate correction rates for each allele.
        Reads are split based on which base they have at the het position, then
        correction percentages are calculated within each allele group.
        """
        temp_dir, crispresso_subdir = crispresso_output_dir
        
        # Het at position 3 (index 2): A vs C
        base_pcts = {
            'A': [0.95, 0.02, 0.50, 0.02],
            'C': [0.02, 0.95, 0.50, 0.95],
            'G': [0.02, 0.02, 0.00, 0.02],
            'T': [0.01, 0.01, 0.00, 0.01],
        }
        file_factory.create_quant_window_file(crispresso_subdir, base_pcts)
        
        # Sequences: het position is index 2
        # Allele A sequences have 'A' at index 2, Allele C have 'C' at index 2
        rows = [
            ("GCAC", 30),  # Allele A, has correction (G at pos 0)
            ("GCCC", 20),  # Allele C, has correction (G at pos 0)
            ("ACAC", 25),  # Allele A, no correction
            ("ACCC", 25),  # Allele C, no correction
        ]
        file_factory.create_alleles_frequency_table(crispresso_subdir, rows)
        
        guide_seq = "ACAC"
        search_strings = ["GCAC", "GCCC"]  # Corrected sequences
        
        result = filter_alleles_file_hetero(search_strings, temp_dir, "F", guide_seq)
        pct_w_base1, pct_w_base2, pct_wo_base1, pct_wo_base2, het_pos, base1, base2 = result
        
        assert het_pos == 3  # 1-indexed for display
        assert set([base1, base2]) == {'A', 'C'}

    def test_returns_NA_when_no_het_position(self, crispresso_output_dir, file_factory):
        """
        Test: Returns all NA values when no het position is found.
        
        Expected return: ("NA", "NA", "NA", "NA", "NA", "NA", "NA")
        
        Reasoning: If sample is homozygous, het-specific analysis isn't applicable.
        """
        temp_dir, crispresso_subdir = crispresso_output_dir
        
        # No het - all positions have dominant base
        base_pcts = {
            'A': [0.95, 0.02, 0.02],
            'C': [0.02, 0.95, 0.02],
            'G': [0.02, 0.02, 0.95],
            'T': [0.01, 0.01, 0.01],
        }
        file_factory.create_quant_window_file(crispresso_subdir, base_pcts)
        
        rows = [("ACGT", 100)]
        file_factory.create_alleles_frequency_table(crispresso_subdir, rows)
        
        result = filter_alleles_file_hetero(["GCGT"], temp_dir, "F", "ACGT")
        
        assert result == ("NA", "NA", "NA", "NA", "NA", "NA", "NA")

class TestWriteHeteroFrequencyTables:
    """Tests for write_hetero_frequency_tables() function."""

    def test_creates_two_allele_files(self, crispresso_output_dir, file_factory):
        """
        Test: Creates separate allele frequency files for each het allele.
        
        Expected: Two files created - one for each allele at het position.
        
        Reasoning: For downstream analysis, we need the allele frequency data
        split by het allele so each can be analyzed independently.
        """
        temp_dir, crispresso_subdir = crispresso_output_dir
        
        rows = [
            ("ACAC", 50),  # Has 'A' at index 2
            ("ACCC", 50),  # Has 'C' at index 2
        ]
        alleles_file = file_factory.create_alleles_frequency_table(crispresso_subdir, rows)
        
        write_hetero_frequency_tables(alleles_file, het_pos=2, base1='A', base2='C', directory_path=temp_dir)
        
        # Check files were created (het_pos + 1 for 1-indexed filename)
        file_a = os.path.join(temp_dir, "AQ_allele_frequency_table_pos3_A.csv")
        file_c = os.path.join(temp_dir, "AQ_allele_frequency_table_pos3_C.csv")
        
        assert os.path.exists(file_a)
        assert os.path.exists(file_c)

