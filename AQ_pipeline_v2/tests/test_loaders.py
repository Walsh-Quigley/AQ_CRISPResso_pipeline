import pytest
from loaders.amplicon_list import load_amplicon_list
from config import AmpliconConfig

#Tests for reading amplicon_list.csv

def test_parse_single_row(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "R186W,CGCTGCATTTCTGCTGGGCC,ABE,R,TAGAGCAACAGT,min humanize,15,8\n"
    )
    result = load_amplicon_list(csv_file)
    assert result[0].name == "R186W"
    assert result[0].protospacer == "CGCTGCATTTCTGCTGGGCC"
    assert result[0].editor == "ABE"
    assert result[0].orientation == "R"
    assert result[0].intended_edit == 8
    assert result[0].tolerated_edits == [15]
    assert len(result) == 1

def test_parse_multiple_rows(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "R186W,CGCTGCATTTCTGCTGGGCC,ABE,R,TAGAGCAACAGT,note1,15,8\n"
        "G542X,ATCGATCGATCGATCGATCG,ABE,F,GCTAGCTAGCTA,note2,3,5\n"
    )
    result = load_amplicon_list(csv_file)
    assert len(result) == 2
    assert result[1].name == "G542X"
    assert result[1].intended_edit == 5

def test_parse_multiple_tolerated_edits(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "R186W,CGCTGCATTTCTGCTGGGCC,ABE,R,TAGAGCAACAGT,note,\"4,7,15\",8\n"
    )
    result = load_amplicon_list(csv_file)
    assert result[0].tolerated_edits == [4, 7, 15]

def test_parse_empty_tolerated_edits(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "R186W,CGCTGCATTTCTGCTGGGCC,ABE,R,TAGAGCAACAGT,note,,8\n"
    )
    result = load_amplicon_list(csv_file)
    assert result[0].tolerated_edits == []

def test_parse_invalid_intended_edit(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "R186W,CGCTGCATTTCTGCTGGGCC,ABE,R,TAGAGCAACAGT,note,15,INVALID\n"
    )
    with pytest.raises(ValueError):
        load_amplicon_list(csv_file)