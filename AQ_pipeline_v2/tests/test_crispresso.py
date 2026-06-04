import pytest
from pipeline.crispresso import identify_amplicon, pair_fastq_files, build_window_args
from config import AmpliconConfig


"""Tests for pipeline/crispresso - covers amplicon matching, longest first priority,
stripping random suffixes, paired fastq files, and no matches being found"""

def test_basic_amplicon_match():
    configs = [
        AmpliconConfig(name="PAH1", protospacer="AAAAAAAAAAAAAAAAAAAA", editor="ABE",
                       orientation="F", amplicon="AAAA", intended_edit=5,
                       tolerated_edits=[], note="")
    ]
    result = identify_amplicon("Icosa03_PAH1_UTD_1", configs)
    assert result == configs[0]

def test_case_insensitive_match():
    configs = [
        AmpliconConfig(name="PAH1", protospacer="AAAAAAAAAAAAAAAAAAAA", editor="ABE",
                       orientation="F", amplicon="AAAA", intended_edit=5,
                       tolerated_edits=[], note="")
    ]
    result = identify_amplicon("Icosa03_pah1_UTD_1", configs)
    assert result == configs[0]

def test_longest_first_priority():
    configs = [
        AmpliconConfig(name="PAH", protospacer="TTTTTTTTTTTTTTTTTTTT", editor="ABE",
                       orientation="F", amplicon="TTTT", intended_edit=4,
                       tolerated_edits=[], note=""),
        AmpliconConfig(name="PAH1", protospacer="AAAAAAAAAAAAAAAAAAAA", editor="ABE",
                       orientation="F", amplicon="AAAA", intended_edit=5,
                       tolerated_edits=[], note="")
    ]
    result = identify_amplicon("Icosa03_pah_UTD_1", configs)
    assert result == configs[0]

def test_stripping_random_suffix():
    configs = [
        AmpliconConfig(name="PAH", protospacer="TTTTTTTTTTTTTTTTTTTT", editor="ABE",
                       orientation="F", amplicon="TTTT", intended_edit=4,
                       tolerated_edits=[], note=""),
        AmpliconConfig(name="PAH1", protospacer="AAAAAAAAAAAAAAAAAAAA", editor="ABE",
                       orientation="F", amplicon="AAAA", intended_edit=5,
                       tolerated_edits=[], note="")
    ]
    result = identify_amplicon("Icosa03_pah_UTD_1.PAH1PAH1PAH1", configs)
    assert result == configs[0]

def test_no_match_FORCED_FAIL():
    configs = [
        AmpliconConfig(name="PAH", protospacer="TTTTTTTTTTTTTTTTTTTT", editor="ABE",
                       orientation="F", amplicon="TTTT", intended_edit=4,
                       tolerated_edits=[], note=""),
        AmpliconConfig(name="PAH1", protospacer="AAAAAAAAAAAAAAAAAAAA", editor="ABE",
                       orientation="F", amplicon="AAAA", intended_edit=5,
                       tolerated_edits=[], note="")
    ]
    with pytest.raises(ValueError):
        identify_amplicon("Icosa03_R408W_UTD_1.PAH1PAH1PAH1", configs)
    
def test_pair_standard_R1_R2_order():
    read1, read2 = pair_fastq_files(["Sample_R1.fastq.gz", "Sample_R2.fastq.gz"])
    assert read1 == "Sample_R1.fastq.gz"
    assert read2 == "Sample_R2.fastq.gz"

def test_pair_reversed_order():
    read1, read2 = pair_fastq_files(["Sample_R2.fastq.gz", "Sample_R1.fastq.gz"])
    assert read1 == "Sample_R1.fastq.gz"
    assert read2 == "Sample_R2.fastq.gz"

def test_pair_case_insensitive():
    read1, read2 = pair_fastq_files(["Sample_r1.fastq.gz", "Sample_r2.fastq.gz"])
    assert read1 == "Sample_r1.fastq.gz"
    assert read2 == "Sample_r2.fastq.gz"

def test_pair_two_R1s_FORCED_FAIL():
    with pytest.raises(ValueError):
        pair_fastq_files(["Sample_R1.fastq.gz", "Sample_R1_resequenced.fastq.gz"])

def test_pair_two_R2s_FORCED_FAIL():
    with pytest.raises(ValueError):
        pair_fastq_files(["Sample_R2.fastq.gz", "Sample_R2_resequenced.fastq.gz"])

def test_pair_no_markers_FORCED_FAIL():
    with pytest.raises(ValueError):
        pair_fastq_files(["Sample_a.fastq.gz", "Sample_b.fastq.gz"])

def test_window_args_nuclease():
    cfg = AmpliconConfig(name="x", protospacer="A"*20, editor="NUCLEASE",
                        orientation="F", amplicon="A"*40,
                        intended_edit=None, tolerated_edits=[], min_alignment_score=75)
    assert build_window_args(cfg) == [
        '--quantification_window_center', '-3',
        '--quantification_window_size', '15',
        '--plot_window_size', '15',
        '--default_min_aln_score', '75',
    ]

def test_window_args_abe():
    cfg = AmpliconConfig(name="x", protospacer="A"*20, editor="ABE",
                     orientation="F", amplicon="A"*40,
                     intended_edit=5, tolerated_edits=[])   # min_alignment_score defaults 60
    assert build_window_args(cfg) == [
        '--plot_window_size', '10',
        '--quantification_window_center', '-10',
        '--quantification_window_size', '10',
    ]

def test_window_args_oneqseq():
    cfg = AmpliconConfig(name="x", protospacer="A"*20, editor="ONESEQ",
                     orientation="F", amplicon="A"*40,
                     intended_edit="ONESEQ", tolerated_edits=[])
    assert build_window_args(cfg) == [
        '--plot_window_size', '10', '--quantification_window_center', '-10',
        '--quantification_window_size', '10',
    ]
