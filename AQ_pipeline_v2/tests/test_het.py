import pandas as pd
from analysis.heterozygous import find_het_position

# Tests for het position detection and allele splitting


def test_find_het_position_basic():
    data = {
        "0": [0.95, 0.02, 0.02, 0.01],  # position 0: A dominant
        "1": [0.50, 0.50, 0.00, 0.00],  # position 1: A/C het
        "2": [0.02, 0.95, 0.02, 0.01],  # position 2: C dominant
    }
    df = pd.DataFrame(data, index=["A","C","G","T"])
    index, base1, base2 = find_het_position(df)
    assert index == 1
    assert base1 == "A"
    assert base2 == "C"

def test_find_het_position_edge():
    data = {
        "0": [0.95, 0.02, 0.02, 0.01],  # position 0: A dominant
        "1": [0.60, 0.40, 0.00, 0.00],  # position 1: A/C het
        "2": [0.02, 0.95, 0.02, 0.01],  # position 2: C dominant
    }
    df = pd.DataFrame(data, index=["A","C","G","T"])
    index, base1, base2 = find_het_position(df)
    assert index == 1
    assert base1 == "A"
    assert base2 == "C"

def test_find_het_position_AG_edit_FORCED_FAIL():
    data = {
        "0": [0.95, 0.02, 0.02, 0.01],  # position 0: A dominant
        "1": [0.60, 0.00, 0.40, 0.00],  # position 1: A/C het
        "2": [0.02, 0.95, 0.02, 0.01],  # position 2: C dominant
    }
    df = pd.DataFrame(data, index=["A","C","G","T"])
    index, base1, base2 = find_het_position(df)
    assert index == None
    assert base1 == None
    assert base2 == None   

def test_find_het_position_CT_edit_FORCED_FAIL():
    data = {
        "0": [0.95, 0.02, 0.02, 0.01],  # position 0: A dominant
        "1": [0.00, 0.60, 0.00, 0.40],  # position 1: A/C het
        "2": [0.02, 0.95, 0.02, 0.01],  # position 2: C dominant
    }
    df = pd.DataFrame(data, index=["A","C","G","T"])
    index, base1, base2 = find_het_position(df)
    assert index == None
    assert base1 == None
    assert base2 == None   

def test_find_het_position_no_split_FORCED_FAIL():
    data = {
        "0": [0.95, 0.02, 0.02, 0.01],  # position 0: A dominant
        "1": [0.00, 0.90, 0.00, 0.10],  # position 1: A/C het
        "2": [0.02, 0.95, 0.02, 0.01],  # position 2: C dominant
    }
    df = pd.DataFrame(data, index=["A","C","G","T"])
    index, base1, base2 = find_het_position(df)
    assert index == None
    assert base1 == None
    assert base2 == None

def test_find_het_position_returns_first():
    data = {
        "0": [0.50, 0.50, 0.00, 0.00],  # position 0: A dominant
        "1": [0.50, 0.50, 0.00, 0.10],  # position 1: A/C het
        "2": [0.02, 0.95, 0.02, 0.01],  # position 2: C dominant
    }
    df = pd.DataFrame(data, index=["A","C","G","T"])
    index, base1, base2 = find_het_position(df)
    assert index == 0
    assert base1 == "A"
    assert base2 == "C"

def test_find_het_position_returns_first_valid():
    data = {
        "0": [0.50, 0.00, 0.50, 0.00],  # position 0: A dominant
        "1": [0.50, 0.50, 0.00, 0.10],  # position 1: A/C het
        "2": [0.02, 0.95, 0.02, 0.01],  # position 2: C dominant
    }
    df = pd.DataFrame(data, index=["A","C","G","T"])
    index, base1, base2 = find_het_position(df)
    assert index == 1
    assert base1 == "A"
    assert base2 == "C"

