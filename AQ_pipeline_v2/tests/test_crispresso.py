import pytest
from pipeline.crispresso import identify_amplicon, pair_fastq_files
from config import AmpliconConfig


"""Tests for pipeline/crispresso - covers amplicon matching, longest first priority,
stripping random suffixs, and no matches being found"""

def test_basic_amplicon_match(tmp_path):
    configs = [
        AmpliconConfig(name="PAH1", protospacer="AAAAAAAAAAAAAAAAAAAA", editor="ABE",
                       orientation="F", amplicon="AAAA", intended_edit=5,
                       tolerated_edits=[], note="")
    ]
    result = identify_amplicon("Icosa03_PAH1_UTD_1", configs)
    assert result == configs[0]

def test_case_insensitive_match(tmp_path):
    configs = [
        AmpliconConfig(name="PAH1", protospacer="AAAAAAAAAAAAAAAAAAAA", editor="ABE",
                       orientation="F", amplicon="AAAA", intended_edit=5,
                       tolerated_edits=[], note="")
    ]
    result = identify_amplicon("Icosa03_pah1_UTD_1", configs)
    assert result == configs[0]

def test_longest_first_priority(tmp_path):
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

def test_striping_random_suffix(tmp_path):
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

def test_no_match_forced_fail(tmp_path):
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

    