import pytest
import pandas as pd
from loaders.exports import generate_prism_csv, generate_prism_csv_het
import logging

"""Tests for loaders/exports - covers generate prism formated output, empty dataframe, 
only two replicates, no trailing replicate number warning, and a missing sample column 
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
        "w_bystanders_minus_wo_bystanders": [13.33, 13.33, 13.33],
        "any_AtoG_minus_w_bystanders": [11.67, 11.67, 11.67],
        "any_change_minus_any_AtoG": [10.0, 10.0, 10.0],
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
        "w_bystanders_minus_wo_bystanders": [13.33, 13.33, 13.33],
        "any_AtoG_minus_w_bystanders": [11.67, 11.67, 11.67],
        "any_change_minus_any_AtoG": [10.0, 10.0, 10.0],
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
        "w_bystanders_minus_wo_bystanders": [13.33, 13.33],
        "any_AtoG_minus_w_bystanders": [11.67, 11.67],
        "any_change_minus_any_AtoG": [10.0, 10.0],
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
        "w_bystanders_minus_wo_bystanders": [13.33],
        "any_AtoG_minus_w_bystanders": [11.67],
        "any_change_minus_any_AtoG": [10.0],
        "target_locus": ["X"],
        "perfect_correction": ["X"],
        "corrected_locus_with_bystanders": ["X"],
    })
    with caplog.at_level(logging.WARNING):
        result = generate_prism_csv(df)
    assert "Icosa03_PAH1_UTD" in caplog.text
    assert not result.empty

def test_prism_rep_with_rep_prefix(tmp_path):
    df = pd.DataFrame({
        "sample": ["Sample_rep1", "Sample_rep2", "Sample_rep3"],
        "reads_total": [1000, 1000, 1000],
        "reads_aligned": [900, 900, 900],
        "correction_without_bystanders": [10.0, 11.0, 12.0],
        "correction_with_tolerated_bystanders": [20.0, 21.0, 22.0],
        "correction_with_any_AtoG_change": [30.0, 31.0, 32.0],
        "correction_with_any_change_in_protospacer": [40.0, 41.0, 42.0],
        "w_bystanders_minus_wo_bystanders": [10.0, 10.0, 10.0],
        "any_AtoG_minus_w_bystanders": [10.0, 10.0, 10.0],
        "any_change_minus_any_AtoG": [10.0, 10.0, 10.0],
        "target_locus": ["X", "X", "X"],
        "perfect_correction": ["X", "X", "X"],
        "corrected_locus_with_bystanders": ["X", "X", "X"],
    })
    result = generate_prism_csv(df)
    assert result.iloc[0]["base_sample"] == "Sample"
    assert result.iloc[0]["sample_rep1"] == "Sample_rep1"

def test_generate_prism_csv_het_happy_path():
    df = pd.DataFrame({
        "sample": ["Sample_PCSK9_1", "Sample_PCSK9_2", "Sample_PCSK9_3"],
        "reads_total": [1000, 1000, 1000],   # intentionally NOT in het output
        "reads_aligned": [900, 900, 900],
        "correction_wo_bystanders_allele1": [10.0, 11.0, 12.0],
        "correction_w_bystanders_allele1":  [20.0, 21.0, 22.0],
        "correction_wo_bystanders_allele2": [5.0, 6.0, 7.0],
        "correction_w_bystanders_allele2":  [15.0, 16.0, 17.0],
        "w_bystanders_minus_wo_bystanders_allele1": [10.0, 10.0, 10.0],
        "w_bystanders_minus_wo_bystanders_allele2": [10.0, 10.0, 10.0],
        "any_AtoG_minus_w_bystanders_allele1":      [5.0, 5.0, 5.0],
        "any_AtoG_minus_w_bystanders_allele2":      [5.0, 5.0, 5.0],
        "any_change_minus_any_AtoG_allele1":        [3.0, 3.0, 3.0],
        "any_change_minus_any_AtoG_allele2":        [3.0, 3.0, 3.0],
        "het_position": [10, 10, 10],
        "het_alleles":  ["G/C", "G/C", "G/C"],
        "reads_aligned_allele1": [540, 540, 540],
        "reads_aligned_allele2": [360, 360, 360],
    })

    allele1_df, allele2_df = generate_prism_csv_het(df)

    # one row each (3 reps collapse to one base_sample)
    assert len(allele1_df) == 1
    assert len(allele2_df) == 1

    # base_sample and rep names extracted correctly
    assert allele1_df.iloc[0]["base_sample"] == "Sample_PCSK9"
    assert allele1_df.iloc[0]["sample_rep1"] == "Sample_PCSK9_1"
    assert allele1_df.iloc[0]["sample_rep3"] == "Sample_PCSK9_3"

    # allele1 file has allele1 values
    assert allele1_df.iloc[0]["correction_wo_bystanders_rep1"] == 10.0

    # allele2 file has allele2 values (different from allele1)
    assert allele2_df.iloc[0]["correction_wo_bystanders_rep1"] == 5.0

    # het metadata preserved
    assert allele1_df.iloc[0]["het_position"] == 10
    assert allele1_df.iloc[0]["het_alleles"] == "G/C"

    # reads_total intentionally NOT in either output
    assert "reads_total_rep1" not in allele1_df.columns
    assert "reads_total_rep1" not in allele2_df.columns

    # reads_aligned in the allele files should be per-allele counts, NOT the total
    assert allele1_df.iloc[0]["reads_aligned_rep1"] == 540
    assert allele2_df.iloc[0]["reads_aligned_rep1"] == 360

    # verify the total is NOT what's in the allele files
    assert allele1_df.iloc[0]["reads_aligned_rep1"] != 900


def test_generate_prism_csv_het_returns_tuple():
    df = pd.DataFrame({
        "sample": ["Sample_1"],
        "reads_aligned": [900],
        "reads_aligned_allele1": [540],
        "reads_aligned_allele2": [360],
        "correction_wo_bystanders_allele1": [10.0],
        "correction_w_bystanders_allele1":  [20.0],
        "correction_wo_bystanders_allele2": [5.0],
        "correction_w_bystanders_allele2":  [15.0],
        "w_bystanders_minus_wo_bystanders_allele1": [10.0],
        "w_bystanders_minus_wo_bystanders_allele2": [10.0],
        "any_AtoG_minus_w_bystanders_allele1": [5.0],
        "any_AtoG_minus_w_bystanders_allele2": [5.0],
        "any_change_minus_any_AtoG_allele1": [3.0],
        "any_change_minus_any_AtoG_allele2": [3.0],
        "het_position": [10],
        "het_alleles": ["G/C"],
    })
    result = generate_prism_csv_het(df)
    assert isinstance(result, tuple)
    assert len(result) == 2
    assert isinstance(result[0], pd.DataFrame)
    assert isinstance(result[1], pd.DataFrame)

def test_generate_prism_csv_het_missing_sample_column_FORCED_FAIL():
    df = pd.DataFrame({
        "reads_aligned": [900],
        "correction_wo_bystanders_allele1": [10.0],
        "het_position": [10],
        "het_alleles": ["G/C"],
    })
    with pytest.raises(ValueError):
        generate_prism_csv_het(df)

def test_generate_prism_csv_het_too_many_reps_FORCED_FAIL():
    df = pd.DataFrame({
        "sample": [f"Sample_{i}" for i in range(1, 5)],   # 4 reps
        "reads_aligned": [900] * 4,
        "correction_wo_bystanders_allele1": [10.0] * 4,
        "correction_w_bystanders_allele1":  [20.0] * 4,
        "correction_wo_bystanders_allele2": [5.0] * 4,
        "correction_w_bystanders_allele2":  [15.0] * 4,
        "w_bystanders_minus_wo_bystanders_allele1": [10.0] * 4,
        "w_bystanders_minus_wo_bystanders_allele2": [10.0] * 4,
        "any_AtoG_minus_w_bystanders_allele1": [5.0] * 4,
        "any_AtoG_minus_w_bystanders_allele2": [5.0] * 4,
        "any_change_minus_any_AtoG_allele1": [3.0] * 4,
        "any_change_minus_any_AtoG_allele2": [3.0] * 4,
        "het_position": [10] * 4,
        "het_alleles": ["G/C"] * 4,
    })
    with pytest.raises(ValueError):
        generate_prism_csv_het(df)
