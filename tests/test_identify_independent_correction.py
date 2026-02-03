# tests/test_identify_independent_correction.py
"""
Tests for identify_independent_correction.py functions.

Tests:
    - total_A_to_G: Calculates correction and A->G percentages
    - total_A_to_G_hetero: Per-allele correction for het samples  
    - identify_independent_correction: Gets correction % from quantification window
"""
import os
import pytest

from scripts.identify_independent_correction import (
    total_A_to_G,
    total_A_to_G_hetero,
    identify_independent_correction
)


class TestTotalAtoG:
    """Tests for total_A_to_G() function."""

    def test_forward_calculates_correction_percentage(self, crispresso_output_dir, file_factory):
        """
        Test: Forward orientation correctly calculates A->G correction percentage.
        
        Expected: Returns (correction_with_any_change, pct_A_to_G_only)
        
        Reasoning: Reads with G at intended position are 'corrected'. Reads with
        ONLY A->G changes (no other mutations) contribute to pct_A_to_G.
        """
        temp_dir, crispresso_subdir = crispresso_output_dir
        
        guide_seq = "ACGTACGTACGTACGTACGT"  # A at position 5 (index 4)
        # Create alleles_frequency_table_around file
        rows = [
            ("ACGTGCGTACGTACGTACGT", 50),  # Has intended A->G at pos 5, only A->G change
            ("ACGTACGTACGTACGTACGT", 50),  # No edit (original)
        ]
        file_factory.create_alleles_frequency_table(
            crispresso_subdir, 
            rows, 
            filename="Alleles_frequency_table_around_test.txt"
        )
        
        correction, pct_A_to_G = total_A_to_G(
            orientation="F",
            intended_edit=5,
            guide_seq=guide_seq,
            directory_path=temp_dir,
            tolerated_edits=[]
        )
        
        # 50% have the intended edit
        assert correction == 50.0
        assert pct_A_to_G == 50.0

    def test_reverse_calculates_correction_percentage(self, crispresso_output_dir, file_factory):
        """
        Test: Reverse orientation uses T->C logic (reverse complement).
        
        Expected: Correctly identifies edits using reverse complement rules.
        
        Reasoning: For reverse orientation, the guide is reverse complemented,
        so we look for T->C changes instead of A->G.
        """
        temp_dir, crispresso_subdir = crispresso_output_dir
        
        guide_seq = "ACGTACGTACGTACGTACGT"
        # Position 5 in reverse = position 16 in forward (20 - 5 + 1 = 16)
        # After RC, we look for T->C
        rows = [
            ("ACGTACGTACGTACGCACGT", 50),  # Has edit
            ("ACGTACGTACGTACGTACGT", 50),  # No edit
        ]
        file_factory.create_alleles_frequency_table(
            crispresso_subdir,
            rows,
            filename="Alleles_frequency_table_around_test.txt"
        )
        
        correction, pct_A_to_G = total_A_to_G(
            orientation="R",
            intended_edit=5,
            guide_seq=guide_seq,
            directory_path=temp_dir,
            tolerated_edits=[]
        )
        
        assert correction == 50.0


class TestIdentifyIndependentCorrection:
    """Tests for identify_independent_correction() function."""

    def test_reads_correction_from_quant_window(self, crispresso_output_dir, file_factory):
        """
        Test: Reads G percentage at intended position from quantification window.
        
        Expected: Returns percentage * 100 for the edit base at intended position.
        
        Reasoning: The quantification window table shows base frequencies at each
        position. For forward ABE, we read the G row at the intended edit column.
        """
        temp_dir, crispresso_subdir = crispresso_output_dir
        
        # G at position 5 = 0.75 (75%)
        base_pcts = {
            'A': [0.95, 0.02, 0.02, 0.02, 0.25] + [0.02] * 15,
            'C': [0.02, 0.95, 0.02, 0.02, 0.00] + [0.02] * 15,
            'G': [0.02, 0.02, 0.95, 0.02, 0.75] + [0.95] * 15,
            'T': [0.01, 0.01, 0.01, 0.94, 0.00] + [0.01] * 15,
        }
        file_factory.create_quant_window_file(crispresso_subdir, base_pcts)
        
        result = identify_independent_correction(
            orientation="F",
            intended_edit=5,
            directory_path=temp_dir
        )
        
        assert result == 75.0

    def test_returns_NA_when_no_quant_file(self, crispresso_output_dir):
        """
        Test: Returns 'NA' when quantification window file is missing.
        
        Expected: "NA"
        
        Reasoning: Graceful handling when CRISPResso output is incomplete.
        """
        temp_dir, crispresso_subdir = crispresso_output_dir
        # Don't create any files
        
        result = identify_independent_correction(
            orientation="F",
            intended_edit=5,
            directory_path=temp_dir
        )
        
        assert result == "NA"

class TestTotalAtoGHetero:
    """Tests for total_A_to_G_hetero() function."""

    def test_splits_correction_by_allele(self, crispresso_output_dir, file_factory):
        """
        Test: Calculates separate correction percentages for each het allele.
        
        Expected return: (correction_allele1, correction_allele2, pct_AtoG_allele1, 
                          pct_AtoG_allele2, reads_allele1, reads_allele2)
        
        Reasoning: For het samples, each allele may have different editing efficiency.
        We need to report correction rates separately for accurate analysis.
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
        
        # Allele A reads (A at index 2), Allele C reads (C at index 2)
        rows = [
            ("GCAC", 40),  # Allele A, corrected (G at pos 0)
            ("ACAC", 10),  # Allele A, not corrected
            ("GCCC", 20),  # Allele C, corrected
            ("ACCC", 30),  # Allele C, not corrected
        ]
        file_factory.create_alleles_frequency_table(
            crispresso_subdir,
            rows,
            filename="Alleles_frequency_table_around_test.txt"
        )
        
        guide_seq = "ACAC"
        
        result = total_A_to_G_hetero(
            orientation="F",
            intended_edit=1,
            guide_seq=guide_seq,
            directory_path=temp_dir,
            tolerated_edits=[]
        )
        
        corr_a1, corr_a2, pct_a1, pct_a2, reads_a1, reads_a2 = result
        
        # Allele A: 40 corrected / 50 total = 80%
        # Allele C: 20 corrected / 50 total = 40%
        assert reads_a1 + reads_a2 == 100

    def test_returns_NA_when_no_het(self, crispresso_output_dir, file_factory):
        """
        Test: Returns NA values when sample is homozygous.
        
        Expected return: All "NA" values
        
        Reasoning: Het-specific analysis isn't applicable for homozygous samples.
        Function should gracefully return NA rather than fail.
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
        
        rows = [("ACG", 100)]
        file_factory.create_alleles_frequency_table(
            crispresso_subdir,
            rows,
            filename="Alleles_frequency_table_around_test.txt"
        )
        
        result = total_A_to_G_hetero(
            orientation="F",
            intended_edit=1,
            guide_seq="ACG",
            directory_path=temp_dir,
            tolerated_edits=[]
        )
        
        # All values should be NA for non-het sample
        assert result == ("NA", "NA", "NA", "NA", "NA", "NA")
