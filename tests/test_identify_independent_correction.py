"""
Tests for identify_independent_correction.py
"""
import os
import sys
import tempfile
import shutil
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.identify_independent_correction import total_A_to_G

class TestTotalAtoG:
    """test case for total_A_to_G"""
    @pytest.fixture
    def temp_crispresso_dir(self):
        """Create tempororary directory structure mimicking CRISPResso output"""
        temp_dir = tempfile.mkdtemp()
        crispresso_subdir = os.path.join(temp_dir, "CRISPResso_on_test")
        os.makedirs(crispresso_subdir)
        yield temp_dir, crispresso_subdir
        shutil.rmtree(temp_dir)
    
    def create_alleles_file(self, crispresso_subdir, rows):
        """mock alleles frequency table"""
        alleles_file = os.path.join(crispresso_subdir, "Alleles_frequency_table_around_test.txt")
        header = "Aligned_Sequence\tReference_Sequence\tUnedited\tn_deleted\tn_inserted\tn_mutated\t#Reads\t%Reads"

        with open(alleles_file, 'w') as f:
            f.write(header + "\n")
            for seq, n_reads in rows:
                 row = f"{seq}\t{seq}\tFalse\t0\t0\t0\t{n_reads}\t0.0"
                 f.write(row + "\n")

    def test_forward_intended_edit_only(self, temp_crispresso_dir):
        """forward orientation, 50% have intended edit, 50% dont.
        Since the only change is the intended A->G edit, pct_A_to_G should equal indep.
        expected: indep = 50%, pct_A_to_G = 50%"""
        temp_dir, crispresso_subdir = temp_crispresso_dir
        guide_seq = "ACGTACGTACGTACGTACGT"
        rows = [
            ("ACGTGCGTACGTACGTACGT", 50),  # Has intended A->G at pos 5, only A->G change
            ("ACGTACGTACGTACGTACGT", 50),  # No edit
        ]
        self.create_alleles_file(crispresso_subdir, rows)
        indep, pct_A_to_G = total_A_to_G(
            orientation="F",
            intended_edit=5,
            guide_seq=guide_seq,
            directory_path=temp_dir,
            tolerated_edits=[]
        )

        assert indep == 50.0
        assert pct_A_to_G == 50.0  # Only A->G changes, so equals indep

    def test_reverse_intended_edit_only(self, temp_crispresso_dir):
        """reverse orientation, 50% have intended edit, 50% dont.
        expected: indep = 50%, non_tol_a_to_g = 0%"""
        temp_dir, crispresso_subdir = temp_crispresso_dir
        guide_seq = "ACGTACGTACGTACGTACGT"
        rows = [
            ("ACGTACGTACGTACGCACGT", 50),
            ("ACGTACGTACGTACGTACGT", 50),
        ]
        self.create_alleles_file(crispresso_subdir, rows)
        indep, pct_A_to_G = total_A_to_G(
            orientation="R",
            intended_edit=5,
            guide_seq=guide_seq,
            directory_path=temp_dir,
            tolerated_edits=[]
        )

        assert indep == 50.0
        assert pct_A_to_G == 50.0

    def test_forward_intended_edit_and_A_to_G(self, temp_crispresso_dir):
        """forward orientation, 50% have intended edit, 25% have intendeded edit as well 
        as non-tolerated A to G changes, 25% of them have .
        expected: indep = 50%, non_tol_a_to_g = 0%"""
        temp_dir, crispresso_subdir = temp_crispresso_dir
        guide_seq = "ACGTACGTACGTACGTACGT"
        rows = [
            ("ACGTACGTACGTACGTACGT", 50),
            ("GCGTGCGTACGTACGTACGT", 25),
            ("GCGTACGTACGTACGTACGT", 25),
        ]
        self.create_alleles_file(crispresso_subdir, rows)
        indep, pct_A_to_G = total_A_to_G(
            orientation="F",
            intended_edit=5,
            guide_seq=guide_seq,
            directory_path=temp_dir,
            tolerated_edits=[]
        )

        assert indep == 25.0
        assert pct_A_to_G == 25.0

##need a case with tolerated edits


    def test_reverse_intended_edit_and_A_to_G(self, temp_crispresso_dir):
        """reverse oreintation, 50% have intended edtit"""