import pytest
import pandas as pd
from loaders.exports import generate_prism_csv
import logging

"""Tests for loaders/exports - covers generate prism formated output, empty dataframe, 
only two replicates, no trailing replicate numebr warning, and a missing sample column 
forced fail"""

def test_generate_prism_csv_happy_path():
    df = pd.DataFrame({
        "sample": ["Icosa03_PAH1_UTD_1", "Icosa03_PAH1_UTD_2", "Icosa03_PAH1_UTD_3"],
        "reads_total": [6500, 6400, 6300],
        "reads_aligned": [6000, 5900, 5800],
        "correction_without_bystanders": [15.0, 14.0, 13.0],
        "correction_with_tolerated_bystanders": [28.33, 27.33, 26.33],
        "correction_with_any_AtoG_change": [40.0, 39.0, 38.0],
        "correction_with_any_change_in_protospacer": [50.0, 49.0, 48.0],
        "column E minus column D": [13.33, 13.33, 13.33],
        "column F minus column E": [11.67, 11.67, 11.67],
        "column G minus column F": [10.0, 10.0, 10.0],
        "target_locus": ["X", "X", "X"],
        "perfect_correction": ["X", "X", "X"],
        "corrected_locus_with_bystanders": ["X", "X", "X"],
    })

    result = generate_prism_csv(df)
    assert result.iloc[0]["base_sample"] == "Icosa03_PAH1_UTD"
    assert result.iloc[0]["sample_rep1"] == "Icosa03_PAH1_UTD_1"
    assert result.iloc[0]["reads_aligned_rep1"] == 6000
    assert result.iloc[0]["w_bystanders_minus_wo_bystanders_rep1"] == 13.33

def test_missing_sample_column_FORCED_FAIL():
    df = pd.DataFrame({
        "reads_total": [6500, 6400, 6300],
        "reads_aligned": [6000, 5900, 5800],
        "correction_without_bystanders": [15.0, 14.0, 13.0],
        "correction_with_tolerated_bystanders": [28.33, 27.33, 26.33],
        "correction_with_any_AtoG_change": [40.0, 39.0, 38.0],
        "correction_with_any_change_in_protospacer": [50.0, 49.0, 48.0],
        "column E minus column D": [13.33, 13.33, 13.33],
        "column F minus column E": [11.67, 11.67, 11.67],
        "column G minus column F": [10.0, 10.0, 10.0],
        "target_locus": ["X", "X", "X"],
        "perfect_correction": ["X", "X", "X"],
        "corrected_locus_with_bystanders": ["X", "X", "X"],
    })
    with pytest.raises(ValueError):
        generate_prism_csv(df)

def test_empty_dataframe_return_empty():
    df = pd.DataFrame({})
    result = generate_prism_csv(df)
    assert result.empty

def test_only_two_replicates():
    df = pd.DataFrame({
        "sample": ["Icosa03_PAH1_UTD_1", "Icosa03_PAH1_UTD_2"],
        "reads_total": [6500, 6400],
        "reads_aligned": [6000, 5900],
        "correction_without_bystanders": [15.0, 14.0],
        "correction_with_tolerated_bystanders": [28.33, 27.33],
        "correction_with_any_AtoG_change": [40.0, 39.0],
        "correction_with_any_change_in_protospacer": [50.0, 49.0],
        "column E minus column D": [13.33, 13.33],
        "column F minus column E": [11.67, 11.67],
        "column G minus column F": [10.0, 10.0],
        "target_locus": ["X", "X"],
        "perfect_correction": ["X", "X"],
        "corrected_locus_with_bystanders": ["X", "X"],
    })
    result = generate_prism_csv(df)
    assert pd.isna(result.iloc[0]["reads_aligned_rep3"])

def test_no_trailing_rep_number_warns(caplog):
    df = pd.DataFrame({
        "sample": ["Icosa03_PAH1_UTD"],
        "reads_total": [6500],
        "reads_aligned": [6000],
        "correction_without_bystanders": [15.0],
        "correction_with_tolerated_bystanders": [28.33],
        "correction_with_any_AtoG_change": [40.0],
        "correction_with_any_change_in_protospacer": [50.0],
        "column E minus column D": [13.33],
        "column F minus column E": [11.67],
        "column G minus column F": [10.0],
        "target_locus": ["X"],
        "perfect_correction": ["X"],
        "corrected_locus_with_bystanders": ["X"],
    })
    with caplog.at_level(logging.WARNING):
        result = generate_prism_csv(df)
    assert "Icosa03_PAH1_UTD" in caplog.text
    assert not result.empty

