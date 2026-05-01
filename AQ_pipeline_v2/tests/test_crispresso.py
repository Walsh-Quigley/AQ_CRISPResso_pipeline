import pytest
from pipeline.crispresso import identify_amplicon
from config import AmpliconConfig

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
    

    