import pandas as pd
from analysis.heterozygous import find_het_position
from pipeline.quantify import calculate_het_correction, calculate_het_protospacer_metrics
import logging

"""Tests for analysis/heterozygous.py - covers find_het_position (basic detection, 
threshold edges, A/G and C/T edit skip, deletion/non-canonical base skip, multiple 
het positions), calculate_het_correction (basic, multiple positions, insertion 
skip with and without warning), and calculate_het_protospacer_metrics (basic, 
multiple positions, insertion skip F/R, allele2-only insertion, deletion 
non-regression, all-insertion edge case, shorter-than-protospacer skip, warning 
message content)."""



def test_find_primary_het_position_basic():
    data = {
        "0": [0.95, 0.02, 0.02, 0.01],  # position 0: A dominant
        "1": [0.50, 0.50, 0.00, 0.00],  # position 1: A/C het
        "2": [0.02, 0.95, 0.02, 0.01],  # position 2: C dominant
    }
    df = pd.DataFrame(data, index=["A","C","G","T"])
    indexes, base1, base2 = find_het_position(df)
    assert indexes == [1]
    assert base1 == "A"
    assert base2 == "C"

def test_find_primary_het_position_edge():
    data = {
        "0": [0.95, 0.02, 0.02, 0.01],  # position 0: A dominant
        "1": [0.60, 0.40, 0.00, 0.00],  # position 1: A/C het
        "2": [0.02, 0.95, 0.02, 0.01],  # position 2: C dominant
    }
    df = pd.DataFrame(data, index=["A","C","G","T"])
    indexes, base1, base2 = find_het_position(df)
    assert indexes == [1]
    assert base1 == "A"
    assert base2 == "C"

def test_find_primary_het_position_AG_edit_skipped():
    data = {
        "0": [0.95, 0.02, 0.02, 0.01],  # position 0: A dominant
        "1": [0.60, 0.00, 0.40, 0.00],  # position 1: A/C het
        "2": [0.02, 0.95, 0.02, 0.01],  # position 2: C dominant
    }
    df = pd.DataFrame(data, index=["A","C","G","T"])
    indexes, base1, base2 = find_het_position(df)
    assert indexes == []
    assert base1 is None
    assert base2 is None   

def test_find_primary_het_position_CT_edit_skipped():
    data = {
        "0": [0.95, 0.02, 0.02, 0.01],  # position 0: A dominant
        "1": [0.00, 0.60, 0.00, 0.40],  # position 1: A/C het
        "2": [0.02, 0.95, 0.02, 0.01],  # position 2: C dominant
    }
    df = pd.DataFrame(data, index=["A","C","G","T"])
    indexes, base1, base2 = find_het_position(df)
    assert indexes == []
    assert base1 is None
    assert base2 is None   

def test_find_het_position_no_split_skipped():
    data = {
        "0": [0.95, 0.02, 0.02, 0.01],  # position 0: A dominant
        "1": [0.00, 0.90, 0.00, 0.10],  # position 1: A/C het
        "2": [0.02, 0.95, 0.02, 0.01],  # position 2: C dominant
    }
    df = pd.DataFrame(data, index=["A","C","G","T"])
    indexes, base1, base2 = find_het_position(df)
    assert indexes == []
    assert base1 is None
    assert base2 is None

def test_find_het_position_returns_first():
    data = {
        "0": [0.50, 0.50, 0.00, 0.00],  # position 0: A dominant
        "1": [0.50, 0.50, 0.00, 0.10],  # position 1: A/C het
        "2": [0.02, 0.95, 0.02, 0.01],  # position 2: C dominant
    }
    df = pd.DataFrame(data, index=["A","C","G","T"])
    indexes, base1, base2 = find_het_position(df)
    assert indexes == [0, 1]
    assert base1 == "A"
    assert base2 == "C"

def test_find_het_position_returns_first_valid():
    data = {
        "0": [0.50, 0.00, 0.50, 0.00],  # position 0: A dominant
        "1": [0.50, 0.50, 0.00, 0.10],  # position 1: A/C het
        "2": [0.02, 0.95, 0.02, 0.01],  # position 2: C dominant
    }
    df = pd.DataFrame(data, index=["A","C","G","T"])
    indexes, base1, base2 = find_het_position(df)
    assert indexes == [1]
    assert base1 == "A"
    assert base2 == "C"

def test_find_ALL_het_positions():
    data = {
        "0": [0.40, 0.02, 0.02, 0.60],  # position 0: A dominant
        "1": [0.50, 0.50, 0.00, 0.00],  # position 1: A/C het
        "2": [0.02, 0.95, 0.02, 0.01],  # position 2: C dominant
    }
    df = pd.DataFrame(data, index=["A","C","G","T"])
    indexes, base1, base2 = find_het_position(df)
    assert len(indexes) == 2
    assert indexes == [0, 1]
    assert base1 == "A"
    assert base2 == "T"

def test_calculate_het_correction_basic():
    table = pd.DataFrame({
        "Aligned_Sequence": [
            "CCACCHTTTTCCCCCTTTTT",  # allele1, matches search_seqs[0] (wo bystanders)
            "CCGCCHTTTTCCCCCTTTTT",  # allele1, no match
            "CCACCXTATTCCCCCTTTTT",  # allele2, matches search_seqs[1] (w bystanders only)
            "CCGCCXTATTCCCCCTTTTT",  # allele2, no match
        ],
        "%Reads": [40.0, 20.0, 30.0, 10.0]
    })
    search_seqs = ["CCACCHTTTTCCCCCTTTTT", "CCACCXTATTCCCCCTTTTT"]
    het_pos = [5]
    base1 = "H"
    base2 = "X"

    result = calculate_het_correction(table, search_seqs, het_pos, base1, base2)

    assert round(result["correction_w_bystanders_allele1"], 2) == 66.67
    assert round(result["correction_wo_bystanders_allele1"], 2) == 66.67
    assert round(result["correction_w_bystanders_allele2"], 2) == 75.0
    assert round(result["correction_wo_bystanders_allele2"], 2) == 0.0

def test_calculate_het_protospacer_metrics_basic():
    table = pd.DataFrame({
        "Aligned_Sequence": [
            "CCGCCCTATTCCCCCTTTTT",  # allele1, matches search_seqs[0] (wo bystanders)
            "CCACCCTATTCCCCCTTTTT",  # allele1, no match
            "CCGCCTTGTTCCCCCTTTTT",  # allele2, matches search_seqs[1] (w bystanders only)
            "CCACCTTGTTCCCCCTTTTT",  # allele2, no match
        ],
        "%Reads": [40.0, 20.0, 30.0, 10.0]
    })
    protospacer = "CCACCCTATTCCCCCTTTTT"
    intended_edit = 3
    orientation = "F"
    het_pos = [5]
    base1 = "C"
    base2 = "T"
    result = calculate_het_protospacer_metrics(table, 
                                      protospacer, 
                                      intended_edit, 
                                      orientation,
                                      het_pos,
                                      base1,
                                      base2)
    assert round(result["correction_with_any_AtoG_change_allele1"], 2) == 66.67
    assert round(result["correction_with_any_change_in_protospacer_allele1"], 2) == 66.67
    assert round(result["correction_with_any_AtoG_change_allele2"], 2) == 75.00
    assert round(result["correction_with_any_change_in_protospacer_allele2"], 2) == 75.0

def test_calculate_het_correction_multiple_pos():
    table = pd.DataFrame({
        "Aligned_Sequence": [
            "CCACCHTTTTCCHCCTTTTT",  # allele1, matches search_seqs[0] (wo bystanders)
            "CCGCCHTTTTCCHCCTTTTT",  # allele1, no match
            "CCACCXTATTCCXCCTTTTT",  # allele2, matches search_seqs[1] (w bystanders only)
            "CCGCCXTATTCCXCCTTTTT",  # allele2, no match
        ],
        "%Reads": [40.0, 20.0, 30.0, 10.0]
    })
    search_seqs = ["CCACCHTTTTCCHCCTTTTT", "CCACCXTATTCCXCCTTTTT"]
    het_pos = [5,12]
    base1 = "H"
    base2 = "X"
    result = calculate_het_correction(table, search_seqs, het_pos, base1, base2)

    assert round(result["correction_w_bystanders_allele1"], 2) == 66.67
    assert round(result["correction_wo_bystanders_allele1"], 2) == 66.67
    assert round(result["correction_w_bystanders_allele2"], 2) == 75.0
    assert round(result["correction_wo_bystanders_allele2"], 2) == 0.0

def test_calculate_het_protospacer_metrics_multiple_pos():
    table = pd.DataFrame({
        "Aligned_Sequence": [
            "CCGCCCTATTCCCCCTTTTT",  # allele1, matches search_seqs[0] (wo bystanders)
            "CCACCCTATTCCCCCTTTTT",  # allele1, no match
            "CCGCCTTGTTCTCCCTTTTT",  # allele2, matches search_seqs[1] (w bystanders only)
            "CCACCTTGTTCTCCCTTTTT",  # allele2, no match
        ],
        "%Reads": [40.0, 20.0, 30.0, 10.0]
    })
    protospacer = "CCGCCCTATTCCCCCTTTTT"
    intended_edit = 3
    orientation = "F"
    het_pos = [5,11]
    base1 = "C"
    base2 = "T"
    result = calculate_het_protospacer_metrics(table, 
                                      protospacer, 
                                      intended_edit, 
                                      orientation,
                                      het_pos,
                                      base1,
                                      base2)
    assert round(result["correction_with_any_AtoG_change_allele1"], 2) == 66.67
    assert round(result["correction_with_any_change_in_protospacer_allele1"], 2) == 66.67
    assert round(result["correction_with_any_AtoG_change_allele2"], 2) == 75.00
    assert round(result["correction_with_any_change_in_protospacer_allele2"], 2) == 75.0

def test_find_primary_het_position_skip_deletion():
    data = {
        "0": [0.50, 0.02, 0.02, 0.01, 0.00, 0.46],  # position 0: A dominant
        "1": [0.50, 0.50, 0.00, 0.00, 0.00, 0.00],  # position 1: A/C het
        "2": [0.02, 0.95, 0.02, 0.01, 0.00, 0.00],  # position 2: C dominant
    }
    df = pd.DataFrame(data, index=["A","C","G","T","N","-"])
    indexes, base1, base2 = find_het_position(df)
    assert indexes == [1]
    assert base1 == "A"
    assert base2 == "C"

def test_het_single_insertion_F_skipped_and_warned(caplog):
    table = pd.DataFrame({
        "Aligned_Sequence": [
            "CCGCCCTATTCCCCCTTTTT",  # allele1, clean A→G at intended_edit pos 3 — both metrics for allele1
            "CCACCCTATTCCCCCTTTTT",  # allele1, no edit at intended position — counts toward total only
            "CCGCCTTGTTCTCCCTTTTT",  # allele2, A→G at intended + A→G bystander — both metrics for allele2
            "CCACCTTGTTCTCCCTTTTT",  # allele2, no edit at intended position — counts toward total only
            "CCACCTTGTTCTCCCTTTTTT", # length 21 — insertion, should be skipped
        ],
        "%Reads": [40.0, 20.0, 30.0, 5.0, 5.0]
    })
    protospacer = "CCGCCCTATTCCCCCTTTTT"
    intended_edit = 3
    orientation = "F"
    het_pos = [5,11]
    base1 = "C"
    base2 = "T"

    with caplog.at_level(logging.WARNING):
        result = calculate_het_protospacer_metrics(
            table, protospacer, intended_edit, orientation, het_pos, base1, base2
        )

    assert round(result["correction_with_any_AtoG_change_allele1"], 2) == 66.67
    assert round(result["correction_with_any_change_in_protospacer_allele1"], 2) == 66.67
    assert round(result["correction_with_any_AtoG_change_allele2"], 2) == 85.71
    assert round(result["correction_with_any_change_in_protospacer_allele2"], 2) == 85.71
    assert "alignment shift" in caplog.text
    assert "5" in caplog.text

def test_het_single_insertion_R_skipped_and_warned(caplog):
    table = pd.DataFrame({
        "Aligned_Sequence": [
            "AAAAAAAAAC",   # allele1, clean T→C edit at intended position (= A→G on guide strand)
            "AAAAAAAAAT",   # allele1, no edit at intended position
            "AAACAAAAAC",   # allele2, clean T→C edit + het base C at pos 3
            "AAACAAAAAT",   # allele2, no edit at intended position
            "AAAAAAAAACT",  # length 11 — insertion, should be skipped
        ],
        "%Reads": [40.0, 20.0, 30.0, 5.0, 5.0]
    })
    protospacer = "ATTTTTTTTT"
    intended_edit = 1
    orientation = "R"
    het_pos = [3]
    base1 = "A"
    base2 = "C"

    with caplog.at_level(logging.WARNING):
        result = calculate_het_protospacer_metrics(
            table, protospacer, intended_edit, orientation, het_pos, base1, base2
        )

    assert round(result["correction_with_any_AtoG_change_allele1"], 2) == 66.67
    assert round(result["correction_with_any_change_in_protospacer_allele1"], 2) == 66.67
    assert round(result["correction_with_any_AtoG_change_allele2"], 2) == 85.71
    assert round(result["correction_with_any_change_in_protospacer_allele2"], 2) == 85.71
    assert "alignment shift" in caplog.text
    assert "5" in caplog.text

def test_het_single_insertion_allele2_only(caplog):
    table = pd.DataFrame({
        "Aligned_Sequence": [
            "CCGCCCTATTCCCCCTTTTT",   # allele1, clean A→G — both metrics for allele1
            "CCACCCTATTCCCCCTTTTT",   # allele1, no edit at intended position — total only
            "CCGCCTTGTTCTCCCTTTTT",   # allele2, clean A→G + bystander — both metrics for allele2
            "CCACCTTGTTCTCCCTTTTTT",  # length 21 — insertion, T at pos 5 (would be allele2 if not skipped)
        ],
        "%Reads": [40.0, 20.0, 30.0, 10.0]
    })
    protospacer = "CCGCCCTATTCCCCCTTTTT"
    intended_edit = 3
    orientation = "F"
    het_pos = [5, 11]
    base1 = "C"
    base2 = "T"

    with caplog.at_level(logging.WARNING):
        result = calculate_het_protospacer_metrics(
            table, protospacer, intended_edit, orientation, het_pos, base1, base2
        )

    assert round(result["correction_with_any_AtoG_change_allele1"], 2) == 66.67
    assert round(result["correction_with_any_change_in_protospacer_allele1"], 2) == 66.67
    assert round(result["correction_with_any_AtoG_change_allele2"], 2) == 100.0
    assert round(result["correction_with_any_change_in_protospacer_allele2"], 2) == 100.0
    assert "alignment shift" in caplog.text
    assert "10" in caplog.text

def test_deletion_not_skipped(caplog):
    table = pd.DataFrame({
        "Aligned_Sequence": [
            "CCGCCCTATTCCCCCTTTTT",   # allele1, clean A→G — both metrics for allele1
            "CCACCCTATTCCCCCTTTTT",   # allele1, no edit at intended — total only
            "CCGCCTTGTTCTCCCTTTTT",   # allele2, clean A→G + bystander — both metrics for allele2
            "CCGCCTTGTT-TCCCTTTTT",   # allele2, A→G at intended + deletion at non-het position — any_change only
        ],
        "%Reads": [40.0, 20.0, 20.0, 10.0]
    })
    protospacer = "CCGCCCTATTCCCCCTTTTT"
    intended_edit = 3
    orientation = "F"
    het_pos = [5, 11]
    base1 = "C"
    base2 = "T"


    with caplog.at_level(logging.WARNING):
        result = calculate_het_protospacer_metrics(
            table, protospacer, intended_edit, orientation, het_pos, base1, base2
        )

    assert round(result["correction_with_any_AtoG_change_allele1"], 2) == 66.67
    assert round(result["correction_with_any_change_in_protospacer_allele1"], 2) == 66.67
    assert round(result["correction_with_any_AtoG_change_allele2"], 2) == 66.67
    assert round(result["correction_with_any_change_in_protospacer_allele2"], 2) == 100.0
    assert "alignment shift" not in caplog.text

def test_all_rows_are_insertions(caplog):
    table = pd.DataFrame({
        "Aligned_Sequence": [
            "CCGCCCTATTCCCCCTTTTTT",  # 21 — would be allele1 clean A→G
            "CCACCCTATTCCCCCTTTTTT",  # 21 — would be allele1 no edit
            "CCGCCTTGTTCTCCCTTTTTT",  # 21 — would be allele2 clean A→G
            "CCACCTTGTTCTCCCTTTTTT",  # 21 — would be allele2 no edit
        ],
        "%Reads": [40.0, 20.0, 30.0, 10.0]
    })
    protospacer = "CCGCCCTATTCCCCCTTTTT"
    intended_edit = 3
    orientation = "F"
    het_pos = [5, 11]
    base1 = "C"
    base2 = "T"

    with caplog.at_level(logging.WARNING):
        result = calculate_het_protospacer_metrics(
            table, protospacer, intended_edit, orientation, het_pos, base1, base2
        )

    assert result["correction_with_any_AtoG_change_allele1"] == 0.0
    assert result["correction_with_any_change_in_protospacer_allele1"] == 0.0
    assert result["correction_with_any_AtoG_change_allele2"] == 0.0
    assert result["correction_with_any_change_in_protospacer_allele2"] == 0.0
    assert "alignment shift" in caplog.text
    assert "100" in caplog.text

def test_het_shorter_than_protospacer_also_skipped(caplog):
    table = pd.DataFrame({
        "Aligned_Sequence": [
            "CCGCCCTATTCCCCCTTTTT",   # length 20, allele1 clean A→G
            "CCACCCTATTCCCCCTTTTT",   # length 20, allele1 no edit
            "CCGCCTTGTTCTCCCTTTTT",   # length 20, allele2 clean A→G
            "CCGCCCTATTCCCCCTTTT",    # length 19 — shorter, would be allele1 if not skipped
            "CCGCCTTGTTCTCCCTTTTTT",  # length 21 — longer, would be allele2 if not skipped
        ],
        "%Reads": [40.0, 20.0, 20.0, 15.0, 5.0]
    })
    protospacer = "CCGCCCTATTCCCCCTTTTT"
    intended_edit = 3
    orientation = "F"
    het_pos = [5, 11]
    base1 = "C"
    base2 = "T"

    with caplog.at_level(logging.WARNING):
        result = calculate_het_protospacer_metrics(
            table, protospacer, intended_edit, orientation, het_pos, base1, base2
        )

    assert round(result["correction_with_any_AtoG_change_allele1"], 2) == 66.67
    assert round(result["correction_with_any_change_in_protospacer_allele1"], 2) == 66.67
    assert round(result["correction_with_any_AtoG_change_allele2"], 2) == 100.0
    assert round(result["correction_with_any_change_in_protospacer_allele2"], 2) == 100.0
    assert "alignment shift" in caplog.text
    assert "20" in caplog.text

def test_warning_message_content(caplog):
    table = pd.DataFrame({    
        "Aligned_Sequence": [
            "CCGCCCTATTCCCCCTTTTT",   # allele1 clean
            "CCACCCTATTCCCCCTTTTT",   # allele1 no edit
            "CCGCCTTGTTCTCCCTTTTT",   # allele2 clean
            "CCGCCTTGTTCTCCCTTTTTT",  # length 21 — insertion at 7.5%
        ],
        "%Reads": [40.0, 22.5, 30.0, 7.5]
    })
    protospacer = "CCGCCCTATTCCCCCTTTTT"
    intended_edit = 3
    orientation = "F"
    het_pos = [5, 11]
    base1 = "C"
    base2 = "T"

    with caplog.at_level(logging.WARNING):
        result = calculate_het_protospacer_metrics(
            table, protospacer, intended_edit, orientation, het_pos, base1, base2
        )
    assert round(result["correction_with_any_AtoG_change_allele1"], 2) == 64.0
    assert round(result["correction_with_any_change_in_protospacer_allele1"], 2) == 64.0
    assert round(result["correction_with_any_AtoG_change_allele2"], 2) == 100.0
    assert round(result["correction_with_any_change_in_protospacer_allele2"], 2) == 100.0

    # warning content checks
    assert "alignment shift" in caplog.text
    assert "7.50" in caplog.text                    # specific formatted percentage
    assert "insertions" in caplog.text              # explanatory text

def test_het_correction_single_insertion_skipped_and_warned(caplog):
    table = pd.DataFrame({
        "Aligned_Sequence": [
            "CCACCHTTTTCCCCCTTTTT",  # allele1, matches search_seqs[0] (perfect correction)
            "CCGCCHTTTTCCCCCTTTTT",  # allele1, no match (no edit at intended)
            "CCACCXTATTCCCCCTTTTT",  # allele2, matches search_seqs[1] (with bystanders)
            "CCGCCXTATTCCCCCTTTTT",  # allele2, no match
            "CCACCHTTTTCCCCCTTTTTT", # length 21 - insertion, should be skipped
        ],
        "%Reads": [40.0, 20.0, 30.0, 10.0, 5.0]
    })
    search_seqs = ["CCACCHTTTTCCCCCTTTTT", "CCACCXTATTCCCCCTTTTT"]
    het_pos = [5]
    base1 = "H"
    base2 = "X"

    with caplog.at_level(logging.WARNING):
        result = calculate_het_correction(table, search_seqs, het_pos, base1, base2)

    assert round(result["correction_w_bystanders_allele1"], 2) == 66.67
    assert round(result["correction_wo_bystanders_allele1"], 2) == 66.67
    assert round(result["correction_w_bystanders_allele2"], 2) == 75.0
    assert round(result["correction_wo_bystanders_allele2"], 2) == 0.0
    assert "alignment shift" in caplog.text
    assert "5" in caplog.text

def test_het_correction_single_insertion_skipped_and_not_warned(caplog):
    table = pd.DataFrame({
        "Aligned_Sequence": [
            "CCACCHTTTTCCCCCTTTTT",  # allele1, matches search_seqs[0] (perfect correction)
            "CCGCCHTTTTCCCCCTTTTT",  # allele1, no match (no edit at intended)
            "CCACCXTATTCCCCCTTTTT",  # allele2, matches search_seqs[1] (with bystanders)
            "CCGCCXTATTCCCCCTTTTT",  # allele2, no match
            "CCACCHTTTTCCCCCTTTTTT", # length 21 - insertion, should be skipped
        ],
        "%Reads": [40.0, 20.0, 30.0, 10.0, 2.5]
    })
    search_seqs = ["CCACCHTTTTCCCCCTTTTT", "CCACCXTATTCCCCCTTTTT"]
    het_pos = [5]
    base1 = "H"
    base2 = "X"

    with caplog.at_level(logging.WARNING):
        result = calculate_het_correction(table, search_seqs, het_pos, base1, base2)

    assert round(result["correction_w_bystanders_allele1"], 2) == 66.67
    assert round(result["correction_wo_bystanders_allele1"], 2) == 66.67
    assert round(result["correction_w_bystanders_allele2"], 2) == 75.0
    assert round(result["correction_wo_bystanders_allele2"], 2) == 0.0
    assert "alignment shift" not in caplog.text
