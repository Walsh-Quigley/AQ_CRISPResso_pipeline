import pytest
import pandas as pd
from loaders.amplicon_list import load_amplicon_list
from utils.sequences import generate_search_sequences
from analysis.abe import calculate_correction, calculate_protospacer_metrics
from analysis.heterozygous import calculate_het_correction, calculate_het_protospacer_metrics
from config import AmpliconConfig
from loaders.crispresso_output import read_allele_table, read_mapping_stats
from tests.helper import make_quant_window
from pipeline.quantify import quantify_sample

# Full pipeline test using sample data
"""Tests end to end workflow of the pipeline - from load amplicon list to generate search
sequences, load amplicon list to calcualte correction, load amplicon list to calculate 
protospacer metric"""

def test_load_amplicon_list_UNTIL_generate_search_sequences(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "PAH1,TCACAGTTCGGGGGTATACA,ABE,F,XXXXX,P281L Hexa,3 16,5\n"
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
        "PAH1,TCACAGTTCGGGGGTATACA,ABE,F,XXXXX,P281L Hexa,3 16,5\n"
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

def test_load_amplicon_list_UNTIL_calcualte_protospacer_metrics(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "PAH1,TCACAGTTCGGGGGTATACA,ABE,F,XXXXX,P281L Hexa,3 16,5\n"
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

def test_load_amplicon_list_UNTIL_calcualte_het_metrics_reverse(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    CRISPResso_mapping_statistic1 = tmp_path / "PAH1_dir" / "CRISPResso_on_sample" / "CRISPResso_mapping_statistics.txt" 
    CRISPResso_mapping_statistic2 = tmp_path / "R186W_dir" / "CRISPResso_on_sample" / "CRISPResso_mapping_statistics.txt" 
    Alleles_frequency_table_around_sgRNA_XXXX1 = tmp_path / "PAH1_dir" / "CRISPResso_on_sample" / "Alleles_frequency_table_around_sgRNA_XXXX.txt" 
    Alleles_frequency_table_around_sgRNA_XXXX2 = tmp_path / "R186W_dir" / "CRISPResso_on_sample" / "Alleles_frequency_table_around_sgRNA_XXXX.txt" 

    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "PAH1,TCACAGTTCGGGGGTATACA,ABE,F,XXXXX,P281L Hexa,3,5\n"
        "R186W,CAGCAGCCACTCAGAGTCTC,ABE,R,XXXXX,note,13,9\n"
    )

    Quantification_window1 = tmp_path / "PAH1_dir" / "CRISPResso_on_sample" / "Quantification_window_nucleotide_percentage_table.txt"
    Quantification_window2 = tmp_path / "R186W_dir" / "CRISPResso_on_sample" / "Quantification_window_nucleotide_percentage_table.txt"

    # protospacer, orientation, intended index, tolerated index, het positions
    PAH1_df = make_quant_window("TCACAGTTCGGGGGTATACA", "F", (4, "A", "G"), [(2, "A", "G")], [(0, "T", "A")]) 
    R186W_df = make_quant_window("CAGCAGCCACTCAGAGTCTC", "R", (8, "A", "G"), [(12, "A", "G")],  [(0, "C", "G")])
    
    Quantification_window1.parent.mkdir(parents=True)
    PAH1_df.rename_axis("Base").to_csv(Quantification_window1, sep="\t")
    Quantification_window2.parent.mkdir(parents=True)
    R186W_df.rename_axis("Base").to_csv(Quantification_window2, sep="\t")

    CRISPResso_mapping_statistic1.write_text(
        "READS IN INPUTS\tREADS AFTER PREPROCESSING\tREADS ALIGNED\tN_COMPUTED_ALN\tN_CACHED_ALN\tN_COMPUTED_NOTALN\tN_CACHED_NOTALN\n"
        "7000\t6500\t6000\t543\t4261\t490\t2"
    )
    CRISPResso_mapping_statistic2.write_text(
        "READS IN INPUTS\tREADS AFTER PREPROCESSING\tREADS ALIGNED\tN_COMPUTED_ALN\tN_CACHED_ALN\tN_COMPUTED_NOTALN\tN_CACHED_NOTALN\n"
        "10000\t9000\t8000\t543\t4261\t490\t2"
    )

    Alleles_frequency_table_around_sgRNA_XXXX1.write_text(
        "Aligned_Sequence\tReference_Sequence\tUnedited\tn_deleted\tn_inserted\tn_mutated\t#Reads\t%Reads\n"
        #allele 1 vvvvvvv
        "TCACAGTTCGGGGGTATACA\tTCACAGTTCGGGGGTATACA\tTrue\t0\t0\t0\t600\t10\n"     #no correction
        "TCGCAGTTCGGGGGTATACA\tTCGCAGTTCGGGGGTATACA\tFalse\t0\t0\t1\t300\t5\n" #no correction but bystander correction
        "TCACGGTTCGGGGGTATACA\tTCACAGTTCGGGGGTATACA\tFalse\t0\t0\t1\t600\t10\n"     #perfect correction
        "TCGCGGTTCGGGGGTATACA\tTCACAGTTCGGGGGTATACA\tFalse\t0\t0\t2\t600\t10\n"  #correction with tolerated bystander
        "TCACGGTTCGGGGGTGTACA\tTCACAGTTCGGGGGTATACA\tFalse\t0\t0\t3\t600\t10\n"  #correction with A to G change 
        "TCACGGTTC-------TACA\tTCACAGTTCGGGGGTATACA\tFalse\t8\t0\t1\t300\t5\n"  #correction with other change in protospacer
        #allele 2 vvvvvvv
        "ACACGGTTCGGGGGTATACA\tTCACAGTTCGGGGGTATACA\tTrue\t0\t0\t0\t600\t10\n"     #no correction
        "ACACGGTTCGGGGGTATACA\tTCGCAGTTCGGGGGTATACA\tFalse\t0\t0\t1\t300\t5\n" #no correction but bystander correction
        "ACACGGTTCGGGGGTATACA\tTCACAGTTCGGGGGTATACA\tFalse\t0\t0\t1\t600\t10\n"     #perfect correction
        "ACACGGTTCGGGGGTATACA\tTCACAGTTCGGGGGTATACA\tFalse\t0\t0\t2\t600\t10\n"  #correction with tolerated bystander
        "ACACGGTTCGGGGGTGTACA\tTCACAGTTCGGGGGTATACA\tFalse\t0\t0\t3\t600\t10\n"  #correction with A to G change 
        "ACACGGTTC-------TACA\tTCACAGTTCGGGGGTATACA\tFalse\t8\t0\t1\t300\t5\n"  #correction with other change in protospacer
    )
    Alleles_frequency_table_around_sgRNA_XXXX2.write_text(
        "Aligned_Sequence\tReference_Sequence\tUnedited\tn_deleted\tn_inserted\tn_mutated\t#Reads\t%Reads\n"
        #allele 1 vvvvvvv
        "GAGACTCTGAGTGGCTGCTG\tGAGACTCTGAGTGGCTGCTG\tTrue\t0\t0\t0\t800\t10\n"      #no correction
        "GAGACTCTGAGTGGCCGCTG\tGAGACTCTGAGTGGCTGCTG\tFalse\t0\t0\t1\t400\t5\n"      #no correction but bystander correction
        "GAGACTCTGAGCGGCTGCTG\tGAGACTCTGAGTGGCTGCTG\tFalse\t0\t0\t1\t800\t10\n"     #perfect correction
        "GAGACTCCGAGCGGCTGCTG\tGAGACTCTGAGTGGCTGCTG\tFalse\t0\t0\t2\t800\t10\n"     #correction with tolerated bystander
        "GAGACTCTGAGCGGCCGCCG\tGAGACTCTGAGTGGCTGCTG\tFalse\t0\t0\t3\t800\t10\n"     #correction with A to G change 
        "G-----CTGAGCGGCTGCTG\tGAGACTCTGAGTGGCTGCTG\tFalse\t8\t0\t1\t400\t5\n"      #correction with other change in protospacer
        #allele 2 vvvvvvv
        "GAGACTCTGAGCGGCTGCTC\tGAGACTCTGAGTGGCTGCTG\tTrue\t0\t0\t0\t800\t10\n"      #no correction
        "GAGACTCTGAGCGGCCGCTC\tGAGACTCTGAGTGGCTGCTG\tFalse\t0\t0\t1\t400\t5\n"      #no correction but bystander correction
        "GAGACTCTGAGCGGCTGCTC\tGAGACTCTGAGTGGCTGCTG\tFalse\t0\t0\t1\t800\t10\n"     #perfect correction
        "GAGACTCCGAGCGGCTGCTC\tGAGACTCTGAGTGGCTGCTG\tFalse\t0\t0\t2\t800\t10\n"     #correction with tolerated bystander
        "GAGACTCTGAGCGGCCGCCC\tGAGACTCTGAGTGGCTGCTG\tFalse\t0\t0\t3\t800\t10\n"     #correction with A to G change 
        "G-----CTGAGCGGCTGCTC\tGAGACTCTGAGTGGCTGCTG\tFalse\t8\t0\t1\t400\t5\n"      #correction with other change in protospacer
    )

    amplicon_list = load_amplicon_list(csv_file)
    PAH1_dict = quantify_sample(amplicon_list[0], tmp_path / "PAH1_dir")
    R186W_dict = quantify_sample(amplicon_list[1], tmp_path / "R186W_dir")

    assert PAH1_dict["reads_total"] == 6500
    assert R186W_dict["reads_total"] == 9000

    assert PAH1_dict["reads_aligned"] == 6000
    assert R186W_dict["reads_aligned"] == 8000

    assert PAH1_dict["correction_wo_bystanders_allele1"] == 20
    assert PAH1_dict["correction_wo_bystanders_allele2"] == 0
    assert R186W_dict["correction_wo_bystanders_allele1"] == 20
    assert R186W_dict["correction_wo_bystanders_allele2"] == 0

    assert PAH1_dict["correction_w_bystanders_allele1"] == 40
    assert PAH1_dict["correction_w_bystanders_allele2"] == 0
    assert R186W_dict["correction_w_bystanders_allele1"] == 40
    assert R186W_dict["correction_w_bystanders_allele2"] == 0

    assert PAH1_dict["correction_with_any_AtoG_change_allele1"] == 60
    assert PAH1_dict["correction_with_any_AtoG_change_allele2"] == 90
    assert R186W_dict["correction_with_any_AtoG_change_allele1"] == 60
    assert R186W_dict["correction_with_any_AtoG_change_allele2"] == 90

    assert PAH1_dict["correction_with_any_change_in_protospacer_allele1"] == 70
    assert PAH1_dict["correction_with_any_change_in_protospacer_allele2"] == 100
    assert R186W_dict["correction_with_any_change_in_protospacer_allele1"] == 70
    assert R186W_dict["correction_with_any_change_in_protospacer_allele2"] == 100

    assert PAH1_dict["het_position"] == 1
    assert R186W_dict["het_position"] == 20

    assert PAH1_dict["het_alleles"] == "T/A"
    assert R186W_dict["het_alleles"] == "G/C"

def test_non_het_ABE_end_to_end(tmp_path):
    csv_file = tmp_path / "amplicon_list.csv"
    CRISPResso_mapping_statistic1 = tmp_path / "PAH1_dir" / "CRISPResso_on_sample" / "CRISPResso_mapping_statistics.txt" 
    CRISPResso_mapping_statistic2 = tmp_path / "R186W_dir" / "CRISPResso_on_sample" / "CRISPResso_mapping_statistics.txt" 
    Alleles_frequency_table_around_sgRNA_XXXX1 = tmp_path / "PAH1_dir" / "CRISPResso_on_sample" / "Alleles_frequency_table_around_sgRNA_XXXX.txt" 
    Alleles_frequency_table_around_sgRNA_XXXX2 = tmp_path / "R186W_dir" / "CRISPResso_on_sample" / "Alleles_frequency_table_around_sgRNA_XXXX.txt" 

    csv_file.write_text(
        "name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n"
        "PAH1,TCACAGTTCGGGGGTATACA,ABE,F,XXXXX,P281L Hexa,3,5\n"
        "R186W,CAGCAGCCACTCAGAGTCTC,ABE,R,XXXXX,note,13,9\n"
    )

    
    Quantification_window1 = tmp_path / "PAH1_dir" / "CRISPResso_on_sample" / "Quantification_window_nucleotide_percentage_table.txt"
    Quantification_window2 = tmp_path / "R186W_dir" / "CRISPResso_on_sample" / "Quantification_window_nucleotide_percentage_table.txt"

    # protospacer, orientation, intended index, tolerated index, het positions
    PAH1_df = make_quant_window("TCACAGTTCGGGGGTATACA", "F", (4, "A", "G"), [(2, "A", "G")]) 
    R186W_df = make_quant_window("CAGCAGCCACTCAGAGTCTC", "R", (8, "A", "G"), [(12, "A", "G")])
    

    Quantification_window1.parent.mkdir(parents=True)
    PAH1_df.rename_axis("Base").to_csv(Quantification_window1, sep="\t")
    Quantification_window2.parent.mkdir(parents=True)
    R186W_df.rename_axis("Base").to_csv(Quantification_window2, sep="\t")

    CRISPResso_mapping_statistic1.write_text(
        "READS IN INPUTS\tREADS AFTER PREPROCESSING\tREADS ALIGNED\tN_COMPUTED_ALN\tN_CACHED_ALN\tN_COMPUTED_NOTALN\tN_CACHED_NOTALN\n"
        "7000\t6500\t6000\t543\t4261\t490\t2"
    )
    CRISPResso_mapping_statistic2.write_text(
        "READS IN INPUTS\tREADS AFTER PREPROCESSING\tREADS ALIGNED\tN_COMPUTED_ALN\tN_CACHED_ALN\tN_COMPUTED_NOTALN\tN_CACHED_NOTALN\n"
        "10000\t9000\t8000\t543\t4261\t490\t2"
    )

    Alleles_frequency_table_around_sgRNA_XXXX1.write_text(
        "Aligned_Sequence\tReference_Sequence\tUnedited\tn_deleted\tn_inserted\tn_mutated\t#Reads\t%Reads\n"
        "TCACAGTTCGGGGGTATACA\tTCACAGTTCGGGGGTATACA\tTrue\t0\t0\t0\t600\t10\n"     #no correction
        "TCGCAGTTCGGGGGTATACA\tTCGCAGTTCGGGGGTATACA\tFalse\t0\t0\t1\t300\t5\n" #no correction but bystander correction
        "TCACGGTTCGGGGGTATACA\tTCACAGTTCGGGGGTATACA\tFalse\t0\t0\t1\t600\t10\n"     #perfect correction
        "TCGCGGTTCGGGGGTATACA\tTCACAGTTCGGGGGTATACA\tFalse\t0\t0\t2\t600\t10\n"  #correction with tolerated bystander
        "TCACGGTTCGGGGGTGTACA\tTCACAGTTCGGGGGTATACA\tFalse\t0\t0\t3\t600\t10\n"  #correction with A to G change 
        "TCACGGTTC-------TACA\tTCACAGTTCGGGGGTATACA\tFalse\t8\t0\t1\t300\t5\n"  #correction with other change in protospacer
    )
    Alleles_frequency_table_around_sgRNA_XXXX2.write_text(
        "Aligned_Sequence\tReference_Sequence\tUnedited\tn_deleted\tn_inserted\tn_mutated\t#Reads\t%Reads\n"
        #allele 1 vvvvvvv
        "GAGACTCTGAGTGGCTGCTG\tGAGACTCTGAGTGGCTGCTG\tTrue\t0\t0\t0\t800\t10\n"      #no correction
        "GAGACTCTGAGTGGCCGCTG\tGAGACTCTGAGTGGCTGCTG\tFalse\t0\t0\t1\t400\t5\n"      #no correction but bystander correction
        "GAGACTCTGAGCGGCTGCTG\tGAGACTCTGAGTGGCTGCTG\tFalse\t0\t0\t1\t800\t10\n"     #perfect correction
        "GAGACTCCGAGCGGCTGCTG\tGAGACTCTGAGTGGCTGCTG\tFalse\t0\t0\t2\t800\t10\n"     #correction with tolerated bystander
        "GAGACTCTGAGCGGCCGCCG\tGAGACTCTGAGTGGCTGCTG\tFalse\t0\t0\t3\t800\t10\n"     #correction with A to G change 
        "G-----CTGAGCGGCTGCTG\tGAGACTCTGAGTGGCTGCTG\tFalse\t8\t0\t1\t400\t5\n"      #correction with other change in protospacer
    )

    amplicon_list = load_amplicon_list(csv_file)
    PAH1_dict = quantify_sample(amplicon_list[0], tmp_path / "PAH1_dir")
    R186W_dict = quantify_sample(amplicon_list[1], tmp_path / "R186W_dir")

    assert PAH1_dict["reads_total"] == 6500
    assert R186W_dict["reads_total"] == 9000

    assert PAH1_dict["reads_aligned"] == 6000
    assert R186W_dict["reads_aligned"] == 8000

    assert PAH1_dict["correction_without_bystanders"] == 10
    assert R186W_dict["correction_without_bystanders"] == 10

    assert PAH1_dict["correction_with_tolerated_bystanders"] == 20
    assert R186W_dict["correction_with_tolerated_bystanders"] == 20

    assert PAH1_dict["correction_with_any_AtoG_change"] == 30
    assert R186W_dict["correction_with_any_AtoG_change"] == 30

    assert PAH1_dict["correction_with_any_change_in_protospacer"] == 35
    assert R186W_dict["correction_with_any_change_in_protospacer"] == 35
