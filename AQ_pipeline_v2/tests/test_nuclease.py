import pandas as pd
from analysis.nuclease import calculate_frameshift
from tests.helper import make_nuclease_fastq_gz

def test_calcualte_frameshift_basic_mix():
    df = pd.DataFrame({
        "n_inserted": [0, 0, 0, 1],
        "n_deleted":  [0, 1, 3, 0],
        "%Reads":     [40.0, 20.0, 20.0, 20.0],
    })
    result = calculate_frameshift(df)
    assert result["pct_frameshift_indels"] == 40
    assert result["pct_inframe_indels"] == 20

def test_calcualte_frameshift_empty():
    df = pd.DataFrame({

    })
    result = calculate_frameshift(df)
    assert result["pct_frameshift_indels"] == 0
    assert result["pct_inframe_indels"] == 0

def test_calculate_frameshift_insertion_deletion_cancel_inframe():
    df = pd.DataFrame({
        "n_inserted": [0, 0, 1, 1],
        "n_deleted":  [0, 1, 4, 0],
        "%Reads":     [40.0, 10.0, 20.0, 30.0],
    })
    result = calculate_frameshift(df)
    assert result["pct_frameshift_indels"] == 40
    assert result["pct_inframe_indels"] == 20

def test_calcuate_frameshift_all_unmodified():
    df = pd.DataFrame({
        "n_inserted": [0, 0, 0, 0],
        "n_deleted":  [0, 0, 0, 0],
        "%Reads":     [40.0, 20.0, 20.0, 20.0],
    })
    result = calculate_frameshift(df)
    assert result["pct_frameshift_indels"] == 0
    assert result["pct_inframe_indels"] == 0

def test_calcualte_frameshift_insertion_only_frameshift():
    df = pd.DataFrame({
        "n_inserted": [2],
        "n_deleted":  [0],
        "%Reads":     [100],
    })
    result = calculate_frameshift(df)
    assert result["pct_frameshift_indels"] == 100
    assert result["pct_inframe_indels"] == 0