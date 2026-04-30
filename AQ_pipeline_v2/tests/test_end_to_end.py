import pytest
import pandas as pd
from loaders.amplicon_list import load_amplicon_list
from utils.sequences import generate_search_sequences
from analysis.abe import calculate_correction, calculate_protospacer_metrics
from config import AmpliconConfig
from loaders.crispresso_output import read_allele_table, read_mapping_stats

# Full pipeline test using sample data

def test_load_amplicon_list_UNTIL_generate_search_sequences(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "PAH1,TCACAGTTCGGGGGTATACA,ABE,F,XXXXX,P281L Hexa,\"3,16\",5\n"
        "R186W,CAGCAGCCACTCAGAGTCTC,ABE,R,XXXXX,note,13,9\n"
    )
    amplicon_list_object = load_amplicon_list(csv_file)
    
    result1 = generate_search_sequences(protospacer=amplicon_list_object[0].protospacer,
                                        intended_edit=amplicon_list_object[0].intended_edit,
                                        tolerated_edits=amplicon_list_object[0].tolerated_edits,
                                        orientation=amplicon_list_object[0].orientation)

    result2 = generate_search_sequences(protospacer=amplicon_list_object[1].protospacer,
                                        intended_edit=amplicon_list_object[1].intended_edit,
                                        tolerated_edits=amplicon_list_object[1].tolerated_edits,
                                        orientation=amplicon_list_object[1].orientation)

    assert result1 == ["TCACGGTTCGGGGGTATACA", "TCGCGGTTCGGGGGTATACA", "TCACGGTTCGGGGGTGTACA", "TCGCGGTTCGGGGGTGTACA"]
    assert result2 == ["GAGACTCTGAGCGGCTGCTG", "GAGACTCCGAGCGGCTGCTG"]

def test_load_amplicon_list_UNTIL_calcualte_correction(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "PAH1,TCACAGTTCGGGGGTATACA,ABE,F,XXXXX,P281L Hexa,\"3,16\",5\n"
        "R186W,CAGCAGCCACTCAGAGTCTC,ABE,R,XXXXX,note,13,9\n"
    )

    table1 = pd.DataFrame({
        "Aligned_Sequence": ["TCACGGTTCGGGGGTATACA", "TCGCGGTTCGGGGGTATACA" ,"TCACGGTTCGGGGGTGTACA", "TCGCGGTTCGGGGGTGTACA", "CCCCCCCCCCCCCCCCCCCC"],
        "%Reads": [10.0, 20.0, 25.0, 5.0,40.0]
    })

    table2 = pd.DataFrame({
        "Aligned_Sequence": ["GAGACTCTGAGCGGCTGCTG", "GAGACTCCGAGCGGCTGCTG" , "GGGGGGGGGGGGGGGGGGGG"],
        "%Reads": [10.0, 20.0, 70.0]
    })

    amplicon_list_object = load_amplicon_list(csv_file)

    amplicon1_search_seq = generate_search_sequences(protospacer=amplicon_list_object[0].protospacer,
                                        intended_edit=amplicon_list_object[0].intended_edit,
                                        tolerated_edits=amplicon_list_object[0].tolerated_edits,
                                        orientation=amplicon_list_object[0].orientation)

    amplicon2_search_seq = generate_search_sequences(protospacer=amplicon_list_object[1].protospacer,
                                        intended_edit=amplicon_list_object[1].intended_edit,
                                        tolerated_edits=amplicon_list_object[1].tolerated_edits,
                                        orientation=amplicon_list_object[1].orientation)

    result1 = calculate_correction(table1, amplicon1_search_seq)
    result2 = calculate_correction(table2, amplicon2_search_seq)

    assert result1 == (10.0, 60.0)
    assert result2 == (10.0, 30.0)


def test_load_amplicon_list_UNTIL_calcualte_correction(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "PAH1,TCACAGTTCGGGGGTATACA,ABE,F,XXXXX,P281L Hexa,\"3,16\",5\n"
        "R186W,CAGCAGCCACTCAGAGTCTC,ABE,R,XXXXX,note,13,9\n"
    )

    table1 = pd.DataFrame({
        "Aligned_Sequence": ["TCACGGTTCGGGGGTATACA", "TCGCGGTTCGGGGGTATACA" ,"TCACGGTTCGGGGGTGTACA", "TCGCGGTTCGGGGGTGTACA", "TCGCGGTTCGGGGGTGTGCG", "CCCCGCCCCCCCCCCCCCCC"],
        "%Reads": [5.0, 10.0, 15.0, 20.0, 22.5, 27.5]
    })

    table2 = pd.DataFrame({
        "Aligned_Sequence": ["GAGACTCTGAGCGGCTGCTG", "GAGACTCCGAGCGGCTGCTG" , "GAGACTCCGAGCGGCTGCCG", "GAGAC-CTGAGCGCCTGCTG", "GGGGGGGGGGGGGGGGGGGG"],
        "%Reads": [10.0, 15.0, 20.0, 25.0, 30.0]
    })

    amplicon_list_object = load_amplicon_list(csv_file)

    amplicon1_search_seq = generate_search_sequences(protospacer=amplicon_list_object[0].protospacer,
                                        intended_edit=amplicon_list_object[0].intended_edit,
                                        tolerated_edits=amplicon_list_object[0].tolerated_edits,
                                        orientation=amplicon_list_object[0].orientation)

    amplicon2_search_seq = generate_search_sequences(protospacer=amplicon_list_object[1].protospacer,
                                        intended_edit=amplicon_list_object[1].intended_edit,
                                        tolerated_edits=amplicon_list_object[1].tolerated_edits,
                                        orientation=amplicon_list_object[1].orientation)

    result1_with_without = calculate_correction(table1, amplicon1_search_seq)
    result2_with_without = calculate_correction(table2, amplicon2_search_seq)

    result1_AtoG_anychange = calculate_protospacer_metrics(
        table1,
        amplicon_list_object[0].protospacer,
        amplicon_list_object[0].intended_edit,
        amplicon_list_object[0].orientation
    )

    result2_AtoG_anychange = calculate_protospacer_metrics(
        table2,
        amplicon_list_object[1].protospacer,
        amplicon_list_object[1].intended_edit,
        amplicon_list_object[1].orientation
    )

    assert result1_with_without == (5.0, 50.0)
    assert result1_AtoG_anychange == (72.5, 100.0)

    assert result2_with_without == (10.0, 25.0)
    assert result2_AtoG_anychange == (45.0, 70.0)