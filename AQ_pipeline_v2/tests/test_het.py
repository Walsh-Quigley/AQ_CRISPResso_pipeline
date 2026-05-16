import pandas as pd
from analysis.heterozygous import find_het_position
from pipeline.quantify import calculate_het_correction, calculate_het_protospacer_metrics

# Tests for het position detection and allele splitting


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

def test_find_primary_het_position_AG_edit_FORCED_FAIL():
    data = {
        "0": [0.95, 0.02, 0.02, 0.01],  # position 0: A dominant
        "1": [0.60, 0.00, 0.40, 0.00],  # position 1: A/C het
        "2": [0.02, 0.95, 0.02, 0.01],  # position 2: C dominant
    }
    df = pd.DataFrame(data, index=["A","C","G","T"])
    indexes, base1, base2 = find_het_position(df)
    assert indexes == []
    assert base1 == None
    assert base2 == None   

def test_find_primary_het_position_CT_edit_FORCED_FAIL():
    data = {
        "0": [0.95, 0.02, 0.02, 0.01],  # position 0: A dominant
        "1": [0.00, 0.60, 0.00, 0.40],  # position 1: A/C het
        "2": [0.02, 0.95, 0.02, 0.01],  # position 2: C dominant
    }
    df = pd.DataFrame(data, index=["A","C","G","T"])
    indexes, base1, base2 = find_het_position(df)
    assert indexes == []
    assert base1 == None
    assert base2 == None   

def test_find_het_position_no_split_FORCED_FAIL():
    data = {
        "0": [0.95, 0.02, 0.02, 0.01],  # position 0: A dominant
        "1": [0.00, 0.90, 0.00, 0.10],  # position 1: A/C het
        "2": [0.02, 0.95, 0.02, 0.01],  # position 2: C dominant
    }
    df = pd.DataFrame(data, index=["A","C","G","T"])
    indexes, base1, base2 = find_het_position(df)
    assert indexes == []
    assert base1 == None
    assert base2 == None

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
    indexs, base1, base2 = find_het_position(df)
    assert len(indexs) == 2
    assert indexs == [0, 1]
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

def test_caclulate_het_protospacer_metrics_basic():
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

def test_calcualte_het_protospacer_metrics_multiple_pos():
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

