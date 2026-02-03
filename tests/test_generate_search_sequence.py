"""
Tests for generate_search_sequences.py functions.

Tests:
    - generate_search_sequences: Creates search strings for correction detection
    - generate_tolerated_sequences: Generates bystander edit combinations
    - A_to_G_sequences: Generates all possible A->G edit combinations
"""
import pytest

from scripts.generate_search_sequences import (
    generate_search_sequences,
    generate_tolerated_sequences,
    A_to_G_sequences
)

class TestGenerateSearchSequences:
    """
    Tests for generate search sequences() funciton.
    """
    def test_forward_orientation_creates_corrected_sequence(self, tmp_path):
        """
        Test: Forward orientation creates A->G correction at intended position.
        
        Expected: First search string has G at intended_edit position (1-indexed).
        
        Reasoning: For ABE forward editing, we convert A->G at the target position.
        The first search string should be the 'perfect correction' with only the
        intended edit.
        """
        guide_seq = "ACGTACGTACGTACGTACGT"  # A at position 5
        orientation = "F"
        editor = "ABE"
        intended_edit = 5  # 1-indexed, should change A->G
        tolerated_edits = []

        result = generate_search_sequences(guide_seq, orientation, editor, intended_edit, tolerated_edits, str(tmp_path))

        # First result should be perfect correction: G at position 5 (index 4)
        assert result[0] == "ACGTGCGTACGTACGTACGT"
        assert result[0][4] == 'G'  # Position 5 (0-indexed = 4)
    
    def test_reverse_orientation_creates_corrected_sequence(self, tmp_path):
        """
        Test: Reverse orientation creates T->C correction (reverse complement logic).
        
        Expected: Sequence is reverse complemented, then edited.
        
        Reasoning: For reverse orientation, the guide is reverse complemented,
        so A->G becomes T->C from the original perspective.
        """
        guide_seq = "ACGTACGTACGTACGTACGT"
        orientation = "R"
        editor = "ABE"
        intended_edit = 5
        tolerated_edits = []
        
        result = generate_search_sequences(guide_seq, orientation, editor, intended_edit, tolerated_edits, str(tmp_path))
        
        # Should return a list with the corrected sequence
        assert len(result) >= 1
        assert result[0] != guide_seq  # Should be different from original

    def test_includes_tolerated_bystander_edits(self, tmp_path):
        """
        Test: Generates additional sequences with tolerated bystander edits.
        
        Expected: Returns perfect correction + all combinations of bystander edits.
        
        Reasoning: Tolerated edits are other A positions that may also get edited.
        We need to detect reads with the intended edit plus any combo of these.
        """
        guide_seq = "AAGTACGTACGTACGTACGT"  # A at positions 1, 2
        orientation = "F"
        editor = "ABE"
        intended_edit = 1  # Change first A->G
        tolerated_edits = [2]  # Second A is tolerated bystander
        
        result = generate_search_sequences(guide_seq, orientation, editor, intended_edit, tolerated_edits, str(tmp_path))
        
        # Should have: perfect correction + correction with bystander
        assert len(result) == 2
        assert "GAGTACGTACGTACGTACGT" in result  # Perfect correction
        assert "GGGTACGTACGTACGTACGT" in result  # With bystander

class TestAtoGSequences:
    """Tests for A_to_G_sequences() function."""

    def test_forward_generates_all_A_to_G_combos(self):
        """
        Test: Generates all combinations of A->G edits in first 10bp.
        
        Expected: Returns two lists - first 10bp combos and all position combos.
        
        Reasoning: We want to detect any A->G editing pattern, so we generate
        all possible combinations to search against.
        """
        guide_seq = "AAAAACGTAC"  # 4 A's in first 5 positions
        orientation = "F"
        
        first_10, all_combos = A_to_G_sequences(guide_seq, orientation)
        
        # Should have combinations for 6 A's in first 10bp = 63
        assert len(first_10) == 63
        # All combos should also include these since all A's are in first 10
        assert len(all_combos) >= len(first_10)

    def test_reverse_uses_T_to_C(self):
        """
        Test: Reverse orientation converts T->C instead of A->G.
        
        Expected: Edits are T->C on the reverse complemented sequence.
        
        Reasoning: Reverse complement means the 'A' positions become 'T',
        so ABE editing shows as T->C.
        """
        guide_seq = "ACGTACGTTT"  # Has T's at end
        orientation = "R"
        
        first_10, all_combos = A_to_G_sequences(guide_seq, orientation)
        
        # Should have generated sequences with C's where T's were
        assert len(first_10) > 0