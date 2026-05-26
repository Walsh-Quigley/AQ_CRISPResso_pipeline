import pandas as pd
import pytest
from analysis.abe import calculate_correction, calculate_protospacer_metrics
import logging

# Tests for ABE metric calcualtions
"""Tests for analysis/abe.py - covers perfect and tolerated analysis in the 
forward and reverse orientation"""


def test_perfect_correction():
    table = pd.DataFrame({
        "Aligned_Sequence": ["GATCGAACGT", "AATCGAACGT"],
        "%Reads": [60.0, 40.0]
    })
    search_sequences = ["GATCGAACGT"]

    without, with_ = calculate_correction(table, search_sequences)

    assert without == 60.0
    assert with_ == 60.0

def test_perfect_and_tolerated_correction():
    table = pd.DataFrame({
        "Aligned_Sequence": ["GATCGAACGT", "GATCGGACGT" ,"AATCGAACGT"],
        "%Reads": [35.0, 25.0, 40.0]
    })
    search_sequences = ["GATCGAACGT", "GATCGGACGT"]

    without, with_ = calculate_correction(table, search_sequences)

    assert without == 35.0
    assert with_ == 60.0

def test_calculate_protospacer_metrics_forward():
    table = pd.DataFrame({
        "Aligned_Sequence": ["GTTTTTTT",   # intended A→G only → both metrics
                            "GTTTCTTT",   # intended A→G + T→C elsewhere → any_change only
                            "ATTTTTTT"],  # no intended edit → neither
        "%Reads": [40.0, 30.0, 30.0]
    })
    protospacer = "ATTTTTTT"   # A at pos 1

    any_AtoG, any_change = calculate_protospacer_metrics(table, protospacer, intended_edit=1, orientation="F")

    assert any_AtoG == 40.0
    assert any_change == 70.0

def test_calcualte_protospacer_metrics_reverse():
    table = pd.DataFrame({
        "Aligned_Sequence": ["AAAAAAAC",   # intended A→G only → both metrics
                            "AAAAGAAC",   # intended A→G + T→C elsewhere → any_change only
                            "AAAAAAAT"],  # no intended edit → neither
        "%Reads": [40.0, 30.0, 30.0]
    })
    protospacer = "ATTTTTTT"

    any_AtoG, any_change = calculate_protospacer_metrics(table, protospacer, intended_edit=1, orientation="R")

    assert any_AtoG == 40.0
    assert any_change == 70.0

def test_protospacer_metrics_invalid_orientation_FORCED_FAIL():
    with pytest.raises(ValueError):
        table = pd.DataFrame({
            "Aligned_Sequence": ["AAAAAAAC",   # intended A→G only → both metrics
                                "AAAAGAAC",   # intended A→G + T→C elsewhere → any_change only
                                "AAAAAAAT"],  # no intended edit → neither
            "%Reads": [40.0, 30.0, 30.0]
        })
        protospacer = "ATTTTTTT"
        calculate_protospacer_metrics(table, protospacer, intended_edit=1, orientation="X")

def test_calculate_correction_no_matches():
    table = pd.DataFrame({
        "Aligned_Sequence": ["GATCGAACGT", "GATCGGACGT" ,"AATCGAACGT"],
        "%Reads": [35.0, 25.0, 40.0]
    })
    search_sequences = ["AGTCAGTCAG", "AGTCAGTCAG"]

    without, with_ = calculate_correction(table, search_sequences)

    assert without == 0
    assert with_ == 0

def test_calcualte_correction_deletion_at_intended_edit_position():
    table = pd.DataFrame({
        "Aligned_Sequence": ["GATCG-ACGT", "GATCGGACGT" ,"AATCGAACGT"],
        "%Reads": [35.0, 25.0, 40.0]
    })
    search_sequences = ["GATCGAACGT", "GATCGGACGT"]

    without, with_ = calculate_correction(table, search_sequences)
    
    assert without == 0.0
    assert with_ == 25.0

def test_protospacer_metrics_deletion_at_intended_position():
    table = pd.DataFrame({
        "Aligned_Sequence": ["-TTTTTTT",   # deletion at intended position → neither metric
                             "GTTTTTTT",   # A→G at intended position → both metrics
                             "ATTTTTTT"],  # no edit → neither
        "%Reads": [30.0, 40.0, 30.0]
    })
    protospacer = "ATTTTTTT"

    any_AtoG, any_change = calculate_protospacer_metrics(table, protospacer, intended_edit=1, orientation="F")

    assert any_AtoG == 40.0
    assert any_change == 40.0

def test_single_insertion_F_skipped_and_warned(caplog):
    table = pd.DataFrame({
        "Aligned_Sequence": ["GTTTTTTT",   # intended A→G only → both metrics
                            "GTTTCTTT",   # intended A→G + T→C elsewhere → any_change only
                            "ATTGTTTTT"],  # insertion
        "%Reads": [40.0, 30.0, 30.0]
    })
    protospacer = "ATTTTTTT"   # A at pos 1

    with caplog.at_level(logging.WARNING):
        any_AtoG, any_change = calculate_protospacer_metrics(
            table, protospacer, intended_edit=1, orientation="F"
        )

    assert any_AtoG == 40.0
    assert any_change == 70.0
    assert "alignment shift" in caplog.text

def test_single_insertion_R_skipped_and_warned(caplog):
    table = pd.DataFrame({
        "Aligned_Sequence": ["AAAAAAAC",   # intended A→G only → both metrics
                            "AAAAGAAC",   # intended A→G + T→C elsewhere → any_change only
                            "AACAAAAAT"],  # insertion
        "%Reads": [40.0, 30.0, 30.0]
    })
    protospacer = "ATTTTTTT"

    with caplog.at_level(logging.WARNING):
        any_AtoG, any_change = calculate_protospacer_metrics(
            table, protospacer, intended_edit=1, orientation="R"
        )

    assert any_AtoG == 40.0
    assert any_change == 70.0
    assert "alignment shift" in caplog.text

def test_multiple_insertions_cumulative_warning(caplog):
    table = pd.DataFrame({
        "Aligned_Sequence": ["GTTTTTTT",   # intended A→G only → both metrics
                            "GTTTCTTT",   # intended A→G + T→C elsewhere → any_change only
                            "ATTGTTTTT", #insertion
                            "ACCCCCCCC"],  # insertion
        "%Reads": [40.0, 30.0, 15.0, 15.0]
    })
    protospacer = "ATTTTTTT"   # A at pos 1

    with caplog.at_level(logging.WARNING):
        any_AtoG, any_change = calculate_protospacer_metrics(
            table, protospacer, intended_edit=1, orientation="F"
        )

    assert any_AtoG == 40.0
    assert any_change == 70.0
    assert "alignment shift" in caplog.text
    assert "30" in caplog.text

def test_deletion_NOT_skipped(caplog):
    table = pd.DataFrame({
        "Aligned_Sequence": ["GTTTTTTT",   # intended A→G only → both metrics
                            "GTTTCT-T",
                            "TTTTTTTT"],   # intended A→G + T→C elsewhere → any_change only

        "%Reads": [40.0, 30.0, 30.0]
    })
    protospacer = "ATTTTTTT"   # A at pos 1

    with caplog.at_level(logging.WARNING):
        any_AtoG, any_change = calculate_protospacer_metrics(
            table, protospacer, intended_edit=1, orientation="F"
        )

    assert any_AtoG == 40.0
    assert any_change == 70.0
    assert "alignment shift" not in caplog.text

def test_all_rows_are_insertions(caplog):
    table = pd.DataFrame({
        "Aligned_Sequence": ["GTTTTTTCT",   # intended A→G only → both metrics
                            "GTTTCTTTC",   # intended A→G + T→C elsewhere → any_change only
                            "ATTGTTTTT", #insertion
                            "ACCCCCCCC"],  # insertion
        "%Reads": [40.0, 30.0, 15.0, 15.0]
    })
    protospacer = "ATTTTTTT"   # A at pos 1

    with caplog.at_level(logging.WARNING):
        any_AtoG, any_change = calculate_protospacer_metrics(
            table, protospacer, intended_edit=1, orientation="F"
        )

    assert any_AtoG == 0
    assert any_change == 0
    assert "alignment shift" in caplog.text
    assert "100" in caplog.text

def test_aligned_shorter_than_protospacer_also_skipped(caplog):
    table = pd.DataFrame({
        "Aligned_Sequence": ["GTTTTTTTTT",   # intended A→G only → both metrics
                            "GTTTCTTTTT",   # intended A→G + T→C elsewhere → any_change only
                            "GTTGTTTTT", #shorter
                            "ACCCCCCCCCCC"],  # insertion
        "%Reads": [40.0, 30.0, 15.0, 15.0]
    })
    protospacer = "ATTTTTTTTT"   # A at pos 1

    with caplog.at_level(logging.WARNING):
        any_AtoG, any_change = calculate_protospacer_metrics(
            table, protospacer, intended_edit=1, orientation="F"
        )

    assert any_AtoG == 40
    assert any_change == 70
    assert "alignment shift" in caplog.text
    assert "30" in caplog.text

def test_warning_message_content(caplog):
    table = pd.DataFrame({
        "Aligned_Sequence": ["GTTTTTTT",   # intended A→G only → both metrics
                            "GTTTCTTT",   # intended A→G + T→C elsewhere → any_change only
                            "ATTTTTTTT"],  # no intended edit → neither
        "%Reads": [40.0, 52.5, 7.5]
    })
    protospacer = "ATTTTTTT"   # A at pos 1

    any_AtoG, any_change = calculate_protospacer_metrics(table, protospacer, intended_edit=1, orientation="F")

    assert any_AtoG == 40.0
    assert any_change == 92.5
    assert "alignment shift" in caplog.text
    assert "7.5" in caplog.text

