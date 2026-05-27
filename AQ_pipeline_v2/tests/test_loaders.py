import pytest
from loaders.amplicon_list import load_amplicon_list, find_amplicon_list
from config import AmpliconConfig
from loaders.crispresso_output import read_allele_table, read_mapping_stats


"""Tests for loaders/amplicon_list.py - covers CSV parsing, tolerated edit formats,
and amplicon list file discovery."""

def test_parse_single_row(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "R186W,CGCTGCATTTCTGCTGGGCC,ABE,R,AAAAAGGCCCAGCAGAAATGCAGCGTTTTT,min humanize,15,8\n"
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
        "R186W,CGCTGCATTTCTGCTGGGCC,ABE,R,AAAAAGGCCCAGCAGAAATGCAGCGTTTTT,note1,15,8\n"
        "G542X,ATCGATCGATCGATCGATCG,ABE,F,AAAAAATCGATCGATCGATCGATCGTTTTT,note2,3,5\n"
    )
    result = load_amplicon_list(csv_file)
    assert len(result) == 2
    assert result[1].name == "G542X"
    assert result[1].intended_edit == 5

def test_parse_multiple_tolerated_edits_spaces(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "R186W,CGCTGCATTTCTGCTGGGCC,ABE,R,AAAAAGGCCCAGCAGAAATGCAGCGTTTTT,note,4 7 15,8\n"
    )
    result = load_amplicon_list(csv_file)
    assert result[0].tolerated_edits == [4, 7, 15]

def test_parse_multiple_tolerated_edits_commas(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "R186W,CGCTGCATTTCTGCTGGGCC,ABE,R,AAAAAGGCCCAGCAGAAATGCAGCGTTTTT,note,\"4, 7, 15\",8\n"
    )
    result = load_amplicon_list(csv_file)
    assert result[0].tolerated_edits == [4, 7, 15]

def test_parse_multiple_tolerated_edits_mix(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "R186W,CGCTGCATTTCTGCTGGGCC,ABE,R,AAAAAGGCCCAGCAGAAATGCAGCGTTTTT,note,\"4, 7 15\",8\n"
    )
    result = load_amplicon_list(csv_file)
    assert result[0].tolerated_edits == [4, 7, 15]

def test_parse_multiple_tolerated_edits_misentered_FORCED_FAIL(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "R186W,CGCTGCATTTCTGCTGGGCC,ABE,R,AAAAAGGCCCAGCAGAAATGCAGCGTTTTT,note,4, 7, 15,8\n"
    )
    with pytest.raises(ValueError):
        result = load_amplicon_list(csv_file)

def test_parse_empty_tolerated_edits(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "R186W,CGCTGCATTTCTGCTGGGCC,ABE,R,AAAAAGGCCCAGCAGAAATGCAGCGTTTTT,note,,8\n"
    )
    result = load_amplicon_list(csv_file)
    assert result[0].tolerated_edits == []

def test_parse_invalid_intended_edit_FORCED_FAIL(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "R186W,CGCTGCATTTCTGCTGGGCC,ABE,R,AAAAAGGCCCAGCAGAAATGCAGCGTTTTT,note,15,INVALID\n"
    )
    with pytest.raises(ValueError):
        load_amplicon_list(csv_file)

def test_duplicate_name_FORCED_FAIL(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "R186W,CGCTGCATTTCTGCTGGGCC,ABE,R,AAAAAGGCCCAGCAGAAATGCAGCGTTTTT,note1,15,8\n"
        "R186W,ATCGATCGATCGATCGATCG,ABE,F,AAAAAATCGATCGATCGATCGATCGTTTTT,note2,3,5\n"
    )
    with pytest.raises(ValueError):
        load_amplicon_list(csv_file)

def test_non_nt_base_in_amplicon_FORCED_FAIL(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "R186W,CGCTGCATTTCTGCTGGGCC,ABE,R,TAGAGXAACAGT,note1,15,8\n"
        "R186W,ATCGATCGATCGATCGATCG,ABE,F,GCTAGXTAGCTA,note2,3,5\n"
    )
    with pytest.raises(ValueError):
        load_amplicon_list(csv_file)

def test_non_nt_base_in_proto_FORCED_FAIL(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "R186W,CGCTGCATTTXTGCTGGGCC,ABE,R,TAGAGAAACAGT,note1,15,8\n"
        "R186W,ATCGATCGATXGATCGATCG,ABE,F,GCTAGATAGCTA,note2,3,5\n"
    )
    with pytest.raises(ValueError):
        load_amplicon_list(csv_file)

def test_invalid_orientation_FORCED_FAIL(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "R186W,CGCTGCATTTTTGCTGGGCC,ABE,X,TAGAGGAACAGT,note1,15,8\n"
        "R186W,ATCGATCGATTGATCGATCG,ABE,fr,GCTAGGTAGCTA,note2,3,5\n"
    )
    with pytest.raises(ValueError):
        load_amplicon_list(csv_file)

def test_protospacer_not_in_amplicon_FORCED_FAIL(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "R186W,CGCTGCATTTCTGCTGGGCC,ABE,F,TAGAGGAACAGT,note1,15,8\n"
        "R186W,ATCGATCGATCGATCGATCG,ABE,R,GCTAGGTAGCTA,note2,3,5\n"
    )
    with pytest.raises(ValueError):
        load_amplicon_list(csv_file)
        
def test_read_mapping_stats(tmp_path):
    stats_file = tmp_path / "CRISPResso_mapping_statistics.txt"
    stats_file.write_text(
        "READS IN INPUTS\tREADS AFTER PREPROCESSING\tREADS ALIGNED\tN_COMPUTED_ALN\tN_CACHED_ALN\tN_COMPUTED_NOTALN\tN_CACHED_NOTALN\n"
        "61046\t51046\t24357\t1861\t22496\t7338\t19351\n"
    )
    result = read_mapping_stats(stats_file)
    assert result[0] == 51046
    assert result[1] == 24357 

def test_read_mapping_stats_file_not_found_FORCED_FAIL(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_mapping_stats(tmp_path / "nonexistent.txt")

def test_read_mapping_stats_zero_total_reads_FORCED_FAIL(tmp_path):
    stats_file = tmp_path / "CRISPResso_mapping_statistics.txt"
    stats_file.write_text(
        "READS IN INPUTS\tREADS AFTER PREPROCESSING\tREADS ALIGNED\tN_COMPUTED_ALN\tN_CACHED_ALN\tN_COMPUTED_NOTALN\tN_CACHED_NOTALN\n"
        "0\t0\t24357\t1861\t22496\t7338\t19351\n"
    )
    with pytest.raises(ValueError):
        read_mapping_stats(stats_file)

def test_read_mapping_stats_greater_aligned_than_total_FORCED_FAIL(tmp_path):
    stats_file = tmp_path / "CRISPResso_mapping_statistics.txt"
    stats_file.write_text(
        "READS IN INPUTS\tREADS AFTER PREPROCESSING\tREADS ALIGNED\tN_COMPUTED_ALN\tN_CACHED_ALN\tN_COMPUTED_NOTALN\tN_CACHED_NOTALN\n"
        "10\t10\t24357\t1861\t22496\t7338\t19351\n"
    )
    with pytest.raises(ValueError):
        read_mapping_stats(stats_file)

def test_read_allele_table(tmp_path):
    allele_file = tmp_path / "Alleles_frequency_table_around.txt"
    allele_file.write_text(
        "Aligned_Sequence\tReference_Sequence\tUnedited\tn_deleted\tn_inserted\tn_mutated\t#Reads\t%Reads\n"
        "AAGCGAACGT\tAATCGAACGT\tFalse\t0\t0\t1\t11185\t45.92\n"
        "AATCGAACGT\tAATCGAACGT\tTrue\t0\t0\t0\t10150\t41.67\n"
    )
    result = read_allele_table(allele_file)
    assert len(result) == 2
    assert list(result.columns) == ["Aligned_Sequence","Reference_Sequence","Unedited","n_deleted","n_inserted","n_mutated","#Reads","%Reads"] 
    assert result["#Reads"][0] == 11185

def test_read_allele_table_file_not_found_FORCED_FAIL(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_allele_table(tmp_path / "nonexistent.txt")

def test_read_allele_table_empty(tmp_path):
    allele_file = tmp_path / "Alleles_frequency_table_around.txt"
    allele_file.write_text(
        "Aligned_Sequence\tReference_Sequence\tUnedited\tn_deleted\tn_inserted\tn_mutated\t#Reads\t%Reads\n"
    )
    result = read_allele_table(allele_file)
    assert len(result) == 0

def test_find_amplicon_list_one_file(tmp_path):
    amplicon_list_path = tmp_path / "common_amplicon_list.csv"
    amplicon_list_path.touch()
    result = find_amplicon_list(tmp_path)
    assert result == amplicon_list_path

def test_find_amplicon_list_no_file(tmp_path):
    amplicon_list_path = tmp_path / "common_amplicon_list.csv"
    with pytest.raises(FileNotFoundError):
        find_amplicon_list(tmp_path)

def test_find_amplicon_list_multiple_files_FORCED_FAIL(tmp_path):
    amplicon_list_path = tmp_path / "common_amplicon_list1.csv"
    amplicon_list_path.touch()
    amplicon_list_path = tmp_path / "common_amplicon_list2.csv"
    amplicon_list_path.touch()
    with pytest.raises(ValueError):
        find_amplicon_list(tmp_path)
