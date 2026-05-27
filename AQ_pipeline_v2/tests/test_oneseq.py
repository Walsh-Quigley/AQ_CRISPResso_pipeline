import pytest
import pandas as pd
from analysis.oneseq import calculate_oneseq

"""Tests for analysis/oneseq.py - covers calculate_oneseq with matches in both 
first 10bp and full protospacer, no matches, all-match, matches only in full, 
and the 10bp boundary edge cases (position 10 vs position 11)."""

def test_calculate_oneseq_basic():
    table = pd.DataFrame({
        "Aligned_Sequence": ["CCXCCTTTTTCCCCCTTTTT", "CCCCCTTTTTCCXCCTTTTT" ,"CCCCCTTTTTCCCCCTTGTT"],
        "%Reads": [35.0, 25.0, 40.0]
    })
    first_10_seqs = ["CCXCCTTTTTCCCCCTTTTT"]
    full_seqs = ["CCXCCTTTTTCCCCCTTTTT", "CCCCCTTTTTCCXCCTTTTT"]

    first_10, full_proto = calculate_oneseq(table, first_10_seqs, full_seqs)
    assert first_10 == 35.0
    assert full_proto == 60.0

def test_calculate_oneseq_no_matches():
    table = pd.DataFrame({
        "Aligned_Sequence": ["CCXCCTTTTTCCCCCTTTTT", "CCCCCTTTTTCCXCCTTTTT" ,"CCCCCTTTTTCCCCCTTGTT"],
        "%Reads": [35.0, 25.0, 40.0]
    })
    first_10_seqs = ["CCCCCTTTTTCCCCCTTTTT"]
    full_seqs = ["CCCCCTTTTTCCCCCTTTTT", "CCCCCTTTTTCCCCCTTTTT"]

    first_10, full_proto = calculate_oneseq(table, first_10_seqs, full_seqs)
    assert first_10 == 0.0
    assert full_proto == 0.0

def test_calculate_oneseq_all_reads_match():
    table = pd.DataFrame({
        "Aligned_Sequence": ["CXCCCTTTTTCCCCCTTTTT", "CCCCCTTTTTCCXCCTTTTT"],
        "%Reads": [35.0, 65.0]
    })
    first_10_seqs = ["CXCCCTTTTTCCCCCTTTTT"]
    full_seqs = ["CXCCCTTTTTCCCCCTTTTT", "CCCCCTTTTTCCXCCTTTTT"]

    first_10, full_proto = calculate_oneseq(table, first_10_seqs, full_seqs)
    assert first_10 == 35.0
    assert full_proto == 100.0

def test_calculate_oneseq_only_in_full():
    table = pd.DataFrame({
        "Aligned_Sequence": ["CCCCCTTTTTCCXCCTTTTT"],
        "%Reads": [35.0]
    })
    first_10_seqs = []
    full_seqs = ["CCCCCTTTTTCCXCCTTTTT"]

    first_10, full_proto = calculate_oneseq(table,first_10_seqs,full_seqs)
    assert first_10 == 0.0
    assert full_proto == 35.0

def test_calculate_oneseq_pos_10():
    table = pd.DataFrame({
        "Aligned_Sequence": ["CCCCCTTTTXCCCCCTTTTT"],
        "%Reads": [35.0]
    })
    first_10_seqs = ["CCCCCTTTTXCCCCCTTTTT"]
    full_seqs = ["CCCCCTTTTXCCCCCTTTTT"]

    first_10, full_proto = calculate_oneseq(table,first_10_seqs,full_seqs)
    assert first_10 == 35.0
    assert full_proto == 35.0

def test_calculate_oneseq_pos_11():
    table = pd.DataFrame({
        "Aligned_Sequence": ["CCCCCTTTTTXCCCCTTTTT"],
        "%Reads": [35.0]
    })
    first_10_seqs = []
    full_seqs = ["CCCCCTTTTTXCCCCTTTTT"]

    first_10, full_proto = calculate_oneseq(table,first_10_seqs,full_seqs)
    assert first_10 == 0.0
    assert full_proto == 35.0