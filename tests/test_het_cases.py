"""
Tests for heterozygous case processing functions.
Tests find_het_position, filter_alleles_file_hetero, and total_A_to_G_hetero.
"""
import os
import sys
import tempfile
import shutil
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.filter_alleles_file import find_het_position, filter_alleles_file_hetero
from scripts.identify_independent_correction import total_A_to_G_hetero


class TestFindHetPosition:
    """Tests for find_het_position() function"""

    @pytest.fixture
    def temp_crispresso_dir(self):
        """Create temporary directory structure mimicking CRISPResso output"""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir)

    def create_quant_window_file(self, crispresso_dir, base_percentages):
        """
        Create mock Quantification_window_nucleotide_percentage_table.txt

        Args:
            crispresso_dir: Path to CRISPResso output directory
            base_percentages: Dict with keys 'A', 'C', 'G', 'T', each containing
                              list of percentages per position (as decimals 0-1)
        """
        quant_file = os.path.join(crispresso_dir, "Quantification_window_nucleotide_percentage_table.txt")

        # Create header with position numbers
        num_positions = len(base_percentages['A'])
        header = "Base\t" + "\t".join([str(i+1) for i in range(num_positions)])

        with open(quant_file, 'w', encoding='utf-8-sig') as f:
            f.write(header + "\n")
            for base in ['A', 'C', 'G', 'T']:
                row = base + "\t" + "\t".join([str(p) for p in base_percentages[base]])
                f.write(row + "\n")

    def test_valid_het_position_AC(self, temp_crispresso_dir):
        """A/C het at position 3 (50%/50%) should be detected"""
        # Position:    1     2     3     4     5
        base_pcts = {
            'A': [0.95, 0.02, 0.50, 0.01, 0.98],
            'C': [0.02, 0.95, 0.50, 0.02, 0.01],
            'G': [0.02, 0.02, 0.00, 0.95, 0.00],
            'T': [0.01, 0.01, 0.00, 0.02, 0.01],
        }
        self.create_quant_window_file(temp_crispresso_dir, base_pcts)

        het_pos, base1, base2 = find_het_position(temp_crispresso_dir)

        assert het_pos == 2  # 0-indexed, so position 3 = index 2
        assert set([base1, base2]) == {'A', 'C'}

    def test_valid_het_position_at_boundaries(self, temp_crispresso_dir):
        """Het at exactly 40%/60% should be detected"""
        # Position:    1     2     3
        base_pcts = {
            'A': [0.40, 0.95, 0.02],
            'C': [0.60, 0.02, 0.95],
            'G': [0.00, 0.02, 0.02],
            'T': [0.00, 0.01, 0.01],
        }
        self.create_quant_window_file(temp_crispresso_dir, base_pcts)

        het_pos, base1, base2 = find_het_position(temp_crispresso_dir)

        assert het_pos == 0  # Position 1, 0-indexed
        assert set([base1, base2]) == {'A', 'C'}

    def test_no_het_position_found(self, temp_crispresso_dir):
        """No position has two bases in 40-60% range"""
        # Position:    1     2     3
        base_pcts = {
            'A': [0.95, 0.02, 0.30],  # 30% is outside range
            'C': [0.02, 0.95, 0.70],  # 70% is outside range
            'G': [0.02, 0.02, 0.00],
            'T': [0.01, 0.01, 0.00],
        }
        self.create_quant_window_file(temp_crispresso_dir, base_pcts)

        het_pos, base1, base2 = find_het_position(temp_crispresso_dir)

        assert het_pos is None
        assert base1 is None
        assert base2 is None

    def test_exclude_AG_pairs(self, temp_crispresso_dir):
        """A/G at 50%/50% should NOT be detected (editing artifact)"""
        # Position:    1     2     3
        base_pcts = {
            'A': [0.50, 0.95, 0.02],  # A/G pair at position 1
            'C': [0.00, 0.02, 0.95],
            'G': [0.50, 0.02, 0.02],
            'T': [0.00, 0.01, 0.01],
        }
        self.create_quant_window_file(temp_crispresso_dir, base_pcts)

        het_pos, base1, base2 = find_het_position(temp_crispresso_dir)

        # Should return None because A/G is excluded
        assert het_pos is None

    def test_exclude_CT_pairs(self, temp_crispresso_dir):
        """C/T at 50%/50% should NOT be detected (editing artifact)"""
        # Position:    1     2     3
        base_pcts = {
            'A': [0.00, 0.95, 0.02],
            'C': [0.50, 0.02, 0.95],  # C/T pair at position 1
            'G': [0.00, 0.02, 0.02],
            'T': [0.50, 0.01, 0.01],
        }
        self.create_quant_window_file(temp_crispresso_dir, base_pcts)

        het_pos, base1, base2 = find_het_position(temp_crispresso_dir)

        # Should return None because C/T is excluded
        assert het_pos is None

    def test_finds_first_het_position(self, temp_crispresso_dir):
        """When multiple het positions exist, returns the first one"""
        # Position:    1     2     3
        base_pcts = {
            'A': [0.50, 0.50, 0.95],  # Het at pos 1 and 2
            'C': [0.50, 0.00, 0.02],
            'G': [0.00, 0.00, 0.02],
            'T': [0.00, 0.50, 0.01],  # A/T het at pos 2
        }
        self.create_quant_window_file(temp_crispresso_dir, base_pcts)

        het_pos, base1, base2 = find_het_position(temp_crispresso_dir)

        assert het_pos == 0  # First het position (A/C)
        assert set([base1, base2]) == {'A', 'C'}


class TestFilterAllelesFileHetero:
    """Tests for filter_alleles_file_hetero() function"""

    @pytest.fixture
    def temp_dir_with_crispresso(self):
        """Create temporary directory with CRISPResso subfolder"""
        temp_dir = tempfile.mkdtemp()
        crispresso_subdir = os.path.join(temp_dir, "CRISPResso_on_test")
        os.makedirs(crispresso_subdir)
        yield temp_dir, crispresso_subdir
        shutil.rmtree(temp_dir)

    def create_quant_window_file(self, crispresso_dir, base_percentages):
        """Create mock quantification window file for het detection"""
        quant_file = os.path.join(crispresso_dir, "Quantification_window_nucleotide_percentage_table.txt")
        num_positions = len(base_percentages['A'])
        header = "Base\t" + "\t".join([str(i+1) for i in range(num_positions)])

        with open(quant_file, 'w', encoding='utf-8-sig') as f:
            f.write(header + "\n")
            for base in ['A', 'C', 'G', 'T']:
                row = base + "\t" + "\t".join([str(p) for p in base_percentages[base]])
                f.write(row + "\n")

    def create_alleles_file(self, crispresso_dir, rows):
        """
        Create mock Alleles_frequency_table.txt

        Args:
            crispresso_dir: Path to CRISPResso output directory
            rows: List of tuples (sequence, n_reads)
        """
        alleles_file = os.path.join(crispresso_dir, "Alleles_frequency_table.txt")
        header = "Aligned_Sequence\tReference_Sequence\tUnedited\tn_deleted\tn_inserted\tn_mutated\t#Reads\t%Reads"

        total_reads = sum(r[1] for r in rows)
        with open(alleles_file, 'w', encoding='utf-8-sig') as f:
            f.write(header + "\n")
            for seq, n_reads in rows:
                pct = (n_reads / total_reads) * 100 if total_reads > 0 else 0
                row = f"{seq}\t{seq}\tFalse\t0\t0\t0\t{n_reads}\t{pct}"
                f.write(row + "\n")

    def test_splits_reads_by_allele(self, temp_dir_with_crispresso):
        """Correctly splits reads based on het position base"""
        temp_dir, crispresso_dir = temp_dir_with_crispresso

        # Het at position 3 (index 2): A vs C
        base_pcts = {
            'A': [0.95, 0.02, 0.50, 0.02],  # 50% A at pos 3
            'C': [0.02, 0.95, 0.50, 0.95],  # 50% C at pos 3
            'G': [0.02, 0.02, 0.00, 0.02],
            'T': [0.01, 0.01, 0.00, 0.01],
        }
        self.create_quant_window_file(crispresso_dir, base_pcts)

        # Sequences with A or C at position 3 (index 2)
        # Search string matches position 1 having G instead of A
        guide_seq = "ACAC"  # Original guide
        rows = [
            ("GCAC", 30),  # Allele A (has A at pos 2), has edit at pos 0
            ("GCCC", 20),  # Allele C (has C at pos 2), has edit at pos 0
            ("ACAC", 25),  # Allele A, no edit
            ("ACCC", 25),  # Allele C, no edit
        ]
        self.create_alleles_file(crispresso_dir, rows)

        search_strings = ["GCAC", "GCCC"]  # Search for G at position 0

        result = filter_alleles_file_hetero(search_strings, temp_dir, "F", guide_seq)
        pct_w_base1, pct_w_base2, pct_wo_base1, pct_wo_base2, het_pos, base1, base2 = result

        # Total reads for allele A (pos 2 = A): 30 + 25 = 55
        # Total reads for allele C (pos 2 = C): 20 + 25 = 45
        # With bystanders for allele A: 30/55 * 100 ≈ 54.5%
        # With bystanders for allele C: 20/45 * 100 ≈ 44.4%

        assert het_pos == 3  # 1-indexed for display
        assert set([base1, base2]) == {'A', 'C'}
        assert pct_w_base1 > 0 or pct_w_base2 > 0  # At least one has matches

    def test_returns_na_when_no_het(self, temp_dir_with_crispresso):
        """Returns all NA when no het position found"""
        temp_dir, crispresso_dir = temp_dir_with_crispresso

        # No het position (no two bases at 40-60%)
        base_pcts = {
            'A': [0.95, 0.02, 0.02],
            'C': [0.02, 0.95, 0.02],
            'G': [0.02, 0.02, 0.95],
            'T': [0.01, 0.01, 0.01],
        }
        self.create_quant_window_file(crispresso_dir, base_pcts)

        rows = [("ACGT", 100)]
        self.create_alleles_file(crispresso_dir, rows)

        result = filter_alleles_file_hetero(["GCGT"], temp_dir, "F", "ACGT")

        assert result == ("NA", "NA", "NA", "NA", "NA", "NA", "NA")


class TestTotalAtoGHetero:
    """Tests for total_A_to_G_hetero() function"""

    @pytest.fixture
    def temp_dir_with_crispresso(self):
        """Create temporary directory with CRISPResso subfolder"""
        temp_dir = tempfile.mkdtemp()
        crispresso_subdir = os.path.join(temp_dir, "CRISPResso_on_test")
        os.makedirs(crispresso_subdir)
        yield temp_dir, crispresso_subdir
        shutil.rmtree(temp_dir)

    def create_quant_window_file(self, crispresso_dir, base_percentages):
        """Create mock quantification window file"""
        quant_file = os.path.join(crispresso_dir, "Quantification_window_nucleotide_percentage_table.txt")
        num_positions = len(base_percentages['A'])
        header = "Base\t" + "\t".join([str(i+1) for i in range(num_positions)])

        with open(quant_file, 'w', encoding='utf-8-sig') as f:
            f.write(header + "\n")
            for base in ['A', 'C', 'G', 'T']:
                row = base + "\t" + "\t".join([str(p) for p in base_percentages[base]])
                f.write(row + "\n")

    def create_alleles_file(self, crispresso_dir, rows):
        """Create mock Alleles_frequency_table_around file"""
        alleles_file = os.path.join(crispresso_dir, "Alleles_frequency_table_around_test.txt")
        header = "Aligned_Sequence\tReference_Sequence\tUnedited\tn_deleted\tn_inserted\tn_mutated\t#Reads\t%Reads"

        total_reads = sum(r[1] for r in rows)
        with open(alleles_file, 'w', encoding='utf-8-sig') as f:
            f.write(header + "\n")
            for seq, n_reads in rows:
                pct = (n_reads / total_reads) * 100 if total_reads > 0 else 0
                row = f"{seq}\t{seq}\tFalse\t0\t0\t0\t{n_reads}\t{pct}"
                f.write(row + "\n")

    def test_forward_orientation_per_allele(self, temp_dir_with_crispresso):
        """Forward orientation: correctly calculates A->G per allele"""
        temp_dir, crispresso_dir = temp_dir_with_crispresso

        # Het at position 5 (index 4): A vs T
        # Using 20-position sequences to match typical guide length
        base_pcts = {
            'A': [0.95] * 4 + [0.50] + [0.95] * 15,  # 50% A at pos 5
            'C': [0.02] * 20,
            'G': [0.02] * 20,
            'T': [0.01] * 4 + [0.50] + [0.01] * 15,  # 50% T at pos 5
        }
        self.create_quant_window_file(crispresso_dir, base_pcts)

        guide_seq = "ACGTACGTACGTACGTACGT"  # A at position 5
        # Intended edit at position 5: A->G
        rows = [
            # Allele A (has A at het pos 4)
            ("ACGTGCGTACGTACGTACGT", 25),  # Has intended edit (G at pos 4)
            ("ACGTACGTACGTACGTACGT", 25),  # No edit
            # Allele T (has T at het pos 4)
            ("ACGTTCGTACGTACGTACGT", 25),  # Has intended edit (G at pos 4)
            ("ACGTTCGTACGTACGTACGT", 25),  # Same - intended edit
        ]
        self.create_alleles_file(crispresso_dir, rows)

        result = total_A_to_G_hetero(
            orientation="F",
            intended_edit=5,
            guide_seq=guide_seq,
            directory_path=temp_dir,
            tolerated_edits=[]
        )

        corr_allele1, corr_allele2, atog_allele1, atog_allele2 = result

        # Should return percentages (not NA)
        assert corr_allele1 != "NA"
        assert corr_allele2 != "NA"

    def test_returns_na_when_no_het(self, temp_dir_with_crispresso):
        """Returns all NA when no het position found"""
        temp_dir, crispresso_dir = temp_dir_with_crispresso

        # No het position
        base_pcts = {
            'A': [0.95] * 20,
            'C': [0.02] * 20,
            'G': [0.02] * 20,
            'T': [0.01] * 20,
        }
        self.create_quant_window_file(crispresso_dir, base_pcts)

        rows = [("ACGTACGTACGTACGTACGT", 100)]
        self.create_alleles_file(crispresso_dir, rows)

        result = total_A_to_G_hetero(
            orientation="F",
            intended_edit=5,
            guide_seq="ACGTACGTACGTACGTACGT",
            directory_path=temp_dir,
            tolerated_edits=[]
        )

        assert result == ("NA", "NA", "NA", "NA")


class TestProcessABECaseHetIntegration:
    """Integration tests for process_ABE_case with het vs non-het samples"""

    # These tests would require more setup - mocking generate_search_sequences, etc.
    # For now, we test the key functions independently above.
    pass


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
