import pytest
from utils.sequences import reverse_complement, generate_search_sequences

#tests/test_sequences.py
# Tests for reverse complement and search sequence generation

def test_reverse_complement_basic():
    assert reverse_complement("ATCG") == "CGAT"

def test_reverse_complement_all_same():
    assert reverse_complement("AAAA") == "TTTT"

def test_reverse_complement_palindrome():
    assert reverse_complement("AATT") == "AATT"

def test_reverse_complement_lowercase():
    assert reverse_complement("atcg") == "CGAT"

def test_reverse_complement_mixed_case():
    assert reverse_complement("AtCg") == "CGAT"

def test_reverse_complement_invalid_base():
    with pytest.raises(KeyError):
        reverse_complement("ATXG")

def test_generate_search_sequences_no_bystanders():
    result = generate_search_sequences(
        protospacer="aACGAACGAA",
        intended_edit=1,
        tolerated_edits=[],
        orientation="F",
    )
    assert result == ["GACGAACGAA"]

def test_generate_search_sequences_no_bystander_reverse():
    result = generate_search_sequences(
        protospacer="aTTCGGTCG",
        intended_edit=1,
        tolerated_edits=[],
        orientation="R",
    )
    assert result == ["CGACCGAAC"]

def test_generate_search_sequences_one_bystander():
    # position 3 (t) is intended edit, position 7 (a) is the tolerated bystander
    result = generate_search_sequences(
        protospacer="aATACGGAaCGT",
        intended_edit=1,
        tolerated_edits=[9],
        orientation="F",
    )
    assert result == ["GATACGGAACGT", "GATACGGAGCGT"]

def test_generate_search_sequences_two_bystanders():
    # pos 3 intended, pos 1 and pos 7 tolerated
    result = generate_search_sequences(
        protospacer="aTaCGAaCGT",
        intended_edit=3,
        tolerated_edits=[1, 7],
        orientation="F",
    )
    assert result == [
        "ATGCGAACGT",   # intended only
        "GTGCGAACGT",   # intended + pos 1
        "ATGCGAGCGT",   # intended + pos 7
        "GTGCGAGCGT",   # intended + pos 1 + pos 7
    ]

def test_generate_search_sequences_edit_at_position_1():
    result = generate_search_sequences(
        protospacer="AATCGAACGT",
        intended_edit=1,
        tolerated_edits=[],
        orientation="F",
    )
    assert result == ["GATCGAACGT"]

def test_generate_search_sequences_edit_at_last_position():
    result = generate_search_sequences(
        protospacer="AATCGAACGA",
        intended_edit=10,
        tolerated_edits=[],
        orientation="F",
    )
    assert result == ["AATCGAACGG"]

def test_generate_search_sequences_one_bystander_reverse():
    # same input as one_bystander test but orientation R
    result = generate_search_sequences(
        protospacer="AAaCGAaCGT",
        intended_edit=3,
        tolerated_edits=[7],
        orientation="R",
    )
    assert result == ["ACGTTCGCTT", "ACGCTCGCTT"]

def test_generate_search_sequences_two_bystanders_reverse():
    # pos 3 intended, pos 1 and 7 tolerated, reverse orientation
    result = generate_search_sequences(
        protospacer="aAaCGAaCGT",
        intended_edit=3,
        tolerated_edits=[1, 7],
        orientation="R",
    )
    assert result == [
        "ACGTTCGCTT",   # intended only
        "ACGTTCGCTC",   # intended + pos 1
        "ACGCTCGCTT",   # intended + pos 7
        "ACGCTCGCTC",   # intended + pos 1 + pos 7
    ]

def test_generate_search_sequences_invalid_orientation():
    with pytest.raises(ValueError):
        generate_search_sequences(
            protospacer="AATCGAACGT",
            intended_edit=3,
            tolerated_edits=[],
            orientation="X",
        )

def test_generate_search_sequences_position_out_of_range():
    with pytest.raises(IndexError):
        generate_search_sequences(
            protospacer="AATCGAACGT",
            intended_edit=99,
            tolerated_edits=[],
            orientation="F",
        )

def test_generate_search_sequences_wrong_intended_base_forward():
    with pytest.raises(ValueError):
        generate_search_sequences(
            protospacer="GGGGAGGGGGGGGGGGGGGG",
            intended_edit=2,
            tolerated_edits=[],
            orientation="F"
        )

def test_generate_search_sequences_wrong_intended_base_reverse():
    with pytest.raises(ValueError):
        generate_search_sequences(
            protospacer="GGGGAGGGGGGGGGGGGGGG",
            intended_edit=2,
            tolerated_edits=[],
            orientation="R"
        )

def test_generate_search_sequences_wrong_tolerated_base_forward():
    with pytest.raises(ValueError):
        generate_search_sequences(
            protospacer="GGGGAGGGGGGGGGGGGGGG",
            intended_edit=5,
            tolerated_edits=[7, 8, 9],
            orientation="F"
        )

def test_generate_search_sequences_wrong_tolerated_base_reverse():
    with pytest.raises(ValueError):
        generate_search_sequences(
            protospacer="GGGGAGGGGGGGGGGGGGGG",
            intended_edit=5,
            tolerated_edits=[8, 9, 10],
            orientation="R"
        )