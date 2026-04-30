import pandas as pd
import pytest
from analysis.abe import calculate_correction, calculate_protospacer_metrics

# Tests for ABE metric calcualtions


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

def test_protospacer_metrics_invalid_orientation():
    with pytest.raises(ValueError):
        table = pd.DataFrame({
            "Aligned_Sequence": ["AAAAAAAC",   # intended A→G only → both metrics
                                "AAAAGAAC",   # intended A→G + T→C elsewhere → any_change only
                                "AAAAAAAT"],  # no intended edit → neither
            "%Reads": [40.0, 30.0, 30.0]
        })
        protospacer = "ATTTTTTT"
        calculate_protospacer_metrics(table, protospacer, intended_edit=1, orientation="X")
