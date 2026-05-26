import pytest
from pathlib import Path
from pipeline.quantify import quantify_sample
from config import AmpliconConfig
from unittest.mock import patch
import pandas as pd
from tests.helper import make_dummy_crispresso_output

"""Tests for pipeline/quantify.py - covers quantify sample, missing CRISPResso_subfolder
forced fail, and multiple allele frequency tables forced fail"""



def test_quantify_sample(tmp_path):
    CRISPResso_mapping_statistic = tmp_path / "sample_dir" / "CRISPResso_on_sample" / "CRISPResso_mapping_statistics.txt" 
    Quantification_window = tmp_path / "sample_dir" / "CRISPResso_on_sample" / "Quantification_window_nucleotide_percentage_table.txt"
    Alleles_frequency_table_around_sgRNA_XXXX = tmp_path / "sample_dir" / "CRISPResso_on_sample" / "Alleles_frequency_table_around_sgRNA_XXXX.txt" 
    CRISPResso_mapping_statistic.parent.mkdir(parents=True)
    CRISPResso_mapping_statistic.write_text(
        "READS IN INPUTS\tREADS AFTER PREPROCESSING\tREADS ALIGNED\tN_COMPUTED_ALN\tN_CACHED_ALN\tN_COMPUTED_NOTALN\tN_CACHED_NOTALN\n"
        "7000\t6500\t6000\t543\t4261\t490\t2"
    )
    Alleles_frequency_table_around_sgRNA_XXXX.write_text(
        "Aligned_Sequence\tReference_Sequence\tUnedited\tn_deleted\tn_inserted\tn_mutated\t#Reads\t%Reads\n"
        "TGTATACCCCCGAACTGTGA\tTGTATACCCCCGAACTGTGA\tTrue\t0\t0\t0\t2000\t33\n"     #no correction
        "TGTATACCCCCGAGCTGTGA\tTGTATACCCCCGAACTGTGA\tFalse\t0\t0\t1\t1000\t16.67\n" #no correction but bystander correction
        "TGTATGCCCCCGAACTGTGA\tTGTATACCCCCGAACTGTGA\tFalse\t0\t0\t1\t900\t15\n"     #perfect correction
        "TGTATGCCCCCGAGCTGTGA\tTGTATACCCCCGAACTGTGA\tFalse\t0\t0\t2\t800\t13.33\n"  #correction with tolerated bystander
        "TGTATGCCCCCGGGCTGTGA\tTGTATACCCCCGAACTGTGA\tFalse\t0\t0\t3\t700\t11.67\n"  #correction with A to G change 
        "TGTATGCCCC--------GA\tTGTATACCCCCGAACTGTGA\tFalse\t8\t0\t1\t600\t10.00\n"  #correction with other change in protospacer
    )
    configs = [
        AmpliconConfig(name="TEST1", protospacer="TGTATACCCCCGAACTGTGA", editor="ABE", orientation="F", 
                       amplicon="CCTTTTTTTAGATGGCGCTCATTGTGCCTGGCAACTGGTAGCTGGAGGACAGTACTGTATACCCCCGAACTGTGATGGGCTTGGATCCATGTCTGATGTACTGTGTGCAGCAAGACCTCAATCCTTTGGGTGTATGGGTCG", 
                       intended_edit=6, tolerated_edits=[14], note="")
    ]
    Quantification_window.write_text(
        "Base\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\t19\t20\n"
        "A\t0.95\t0.02\t0.02\t0.95\t0.95\t0.95\t0.02\t0.02\t0.02\t0.02\t0.02\t0.02\t0.95\t0.02\t0.02\t0.02\t0.02\t0.02\t0.02\t0.95\n"
        "C\t0.02\t0.02\t0.95\t0.02\t0.02\t0.02\t0.02\t0.95\t0.95\t0.95\t0.95\t0.95\t0.02\t0.02\t0.02\t0.02\t0.02\t0.02\t0.02\t0.02\n"
        "G\t0.02\t0.02\t0.02\t0.02\t0.02\t0.02\t0.02\t0.02\t0.02\t0.02\t0.02\t0.02\t0.02\t0.95\t0.95\t0.02\t0.02\t0.02\t0.02\t0.02\n"
        "T\t0.01\t0.95\t0.01\t0.01\t0.01\t0.01\t0.95\t0.01\t0.01\t0.01\t0.01\t0.01\t0.01\t0.01\t0.01\t0.95\t0.95\t0.95\t0.95\t0.01\n"
    )
    result = quantify_sample(configs[0], tmp_path/"sample_dir")
    assert result["sample"] == "sample_dir"
    assert result["reads_total"] == 6500
    assert result["reads_aligned"] == 6000
    assert result["correction_without_bystanders"] == 15
    assert result["correction_with_tolerated_bystanders"] == 28.33
    assert result["correction_with_any_AtoG_change"] == 40
    assert result["correction_with_any_change_in_protospacer"] == 50
    assert result["w_bystanders_minus_wo_bystanders"] == 13.33
    assert result["any_AtoG_minus_w_bystanders"] == 11.67
    assert result["any_change_minus_any_AtoG"] == 10
    assert result["target_locus"] == "TGTATACCCCCGAACTGTGA"
    assert result["perfect_correction"] == "TGTATGCCCCCGAACTGTGA"
    assert result["corrected_locus_with_bystanders"] == "TGTATGCCCCCGAACTGTGA;TGTATGCCCCCGAGCTGTGA"

def test_missing_CRISPResso_subfolder_FORCED_FAIL(tmp_path):
    configs = [
        AmpliconConfig(name="TEST1", protospacer="TGTATACCCCCGAACTGTGA", editor="ABE", orientation="F", 
                       amplicon="CCTTTTTTTAGATGGCGCTCATTGTGCCTGGCAACTGGTAGCTGGAGGACAGTACTGTATACCCCCGAACTGTGATGGGCTTGGATCCATGTCTGATGTACTGTGTGCAGCAAGACCTCAATCCTTTGGGTGTATGGGTCG", 
                       intended_edit=6, tolerated_edits=[14], note="")
    ]
    with pytest.raises(FileNotFoundError):
        result = quantify_sample(configs[0], tmp_path/"sample_dir")

def test_multiple_allele_frequency_tables_FORCED_FAIL(tmp_path):
    CRISPResso_mapping_statistic = tmp_path / "sample_dir" / "CRISPResso_on_sample" / "CRISPResso_mapping_statistics.txt" 
    Alleles_frequency_table_around_sgRNA_XXXX = tmp_path / "sample_dir" / "CRISPResso_on_sample" / "Alleles_frequency_table_around_sgRNA_XXXX.txt" 
    Alleles_frequency_table_around_sgRNA_YYYY = tmp_path / "sample_dir" / "CRISPResso_on_sample" / "Alleles_frequency_table_around_sgRNA_YYYY.txt" 
    CRISPResso_mapping_statistic.parent.mkdir(parents=True)
    CRISPResso_mapping_statistic.write_text(
        "READS IN INPUTS\tREADS AFTER PREPROCESSING\tREADS ALIGNED\tN_COMPUTED_ALN\tN_CACHED_ALN\tN_COMPUTED_NOTALN\tN_CACHED_NOTALN\n"
        "7000\t6500\t6000\t543\t4261\t490\t2"
    )
    Alleles_frequency_table_around_sgRNA_XXXX.write_text(
        "Aligned_Sequence\tReference_Sequence\tUnedited\tn_deleted\tn_inserted\tn_mutated\t#Reads\t%Reads\n"
        "TGTATACCCCCGAACTGTGA\tTGTATACCCCCGAACTGTGA\tTrue\t0\t0\t0\t2000\t33\n"     #no correction
        "TGTATACCCCCGAGCTGTGA\tTGTATACCCCCGAACTGTGA\tFalse\t0\t0\t1\t1000\t16.67\n" #no correction but bystander correction
        "TGTATGCCCCCGAACTGTGA\tTGTATACCCCCGAACTGTGA\tFalse\t0\t0\t1\t900\t15\n"     #perfect correction
        "TGTATGCCCCCGAGCTGTGA\tTGTATACCCCCGAACTGTGA\tFalse\t0\t0\t2\t800\t13.33\n"  #correction with tolerated bystander
        "TGTATGCCCCCGGGCTGTGA\tTGTATACCCCCGAACTGTGA\tFalse\t0\t0\t3\t700\t11.67\n"  #correction with A to G change 
        "TGTATGCCCC--------GA\tTGTATACCCCCGAACTGTGA\tFalse\t8\t0\t1\t600\t10.00\n"  #correction with other change in protospacer
    )
    Alleles_frequency_table_around_sgRNA_YYYY.write_text(
        "Aligned_Sequence\tReference_Sequence\tUnedited\tn_deleted\tn_inserted\tn_mutated\t#Reads\t%Reads\n"
        "TGTATACCCCCGAACTGTGA\tTGTATACCCCCGAACTGTGA\tTrue\t0\t0\t0\t2000\t33\n"     #no correction
        "TGTATACCCCCGAGCTGTGA\tTGTATACCCCCGAACTGTGA\tFalse\t0\t0\t1\t1000\t16.67\n" #no correction but bystander correction
        "TGTATGCCCCCGAACTGTGA\tTGTATACCCCCGAACTGTGA\tFalse\t0\t0\t1\t900\t15\n"     #perfect correction
        "TGTATGCCCCCGAGCTGTGA\tTGTATACCCCCGAACTGTGA\tFalse\t0\t0\t2\t800\t13.33\n"  #correction with tolerated bystander
        "TGTATGCCCCCGGGCTGTGA\tTGTATACCCCCGAACTGTGA\tFalse\t0\t0\t3\t700\t11.67\n"  #correction with A to G change 
        "TGTATGCCCC--------GA\tTGTATACCCCCGAACTGTGA\tFalse\t8\t0\t1\t600\t10.00\n"  #correction with other change in protospacer
    )
    configs = [
        AmpliconConfig(name="TEST1", protospacer="TGTATACCCCCGAACTGTGA", editor="ABE", orientation="F", 
                       amplicon="CCTTTTTTTAGATGGCGCTCATTGTGCCTGGCAACTGGTAGCTGGAGGACAGTACTGTATACCCCCGAACTGTGATGGGCTTGGATCCATGTCTGATGTACTGTGTGCAGCAAGACCTCAATCCTTTGGGTGTATGGGTCG", 
                       intended_edit=6, tolerated_edits=[14], note="")
    ]
    with pytest.raises(ValueError):
        result = quantify_sample(configs[0], tmp_path/"sample_dir")

def test_dispatch_routes_ONESEQ_to_oneseq_sample(tmp_path):
    sample_dir = make_dummy_crispresso_output(tmp_path)
    config = AmpliconConfig(
        name="TEST", protospacer="A" * 20, editor="ONESEQ",
        orientation="F", amplicon="A" * 30, intended_edit="ONESEQ",
        tolerated_edits=[], note=""
    )

    with patch("pipeline.quantify.read_mapping_stats", return_value=(100, 80)), \
         patch("pipeline.quantify.read_allele_table", return_value=pd.DataFrame()), \
         patch("pipeline.quantify.read_quant_window", return_value=pd.DataFrame()), \
         patch("pipeline.quantify.find_het_position", return_value=([], None, None)), \
         patch("pipeline.quantify.quantify_oneseq_sample", return_value={"branch": "oneseq"}) as mock_oneseq, \
         patch("pipeline.quantify.quantify_abe_sample", return_value={"branch": "abe"}) as mock_abe, \
         patch("pipeline.quantify.quantify_het_sample", return_value={"branch": "het"}) as mock_het:
        result = quantify_sample(config, sample_dir)

    mock_oneseq.assert_called_once()
    mock_abe.assert_not_called()
    mock_het.assert_not_called()
    assert result == {"branch": "oneseq"}


def test_dispatch_routes_ABE_non_het_to_abe_sample(tmp_path):
    sample_dir = make_dummy_crispresso_output(tmp_path)
    config = AmpliconConfig(
        name="TEST", protospacer="A" * 20, editor="ABE",
        orientation="F", amplicon="A" * 30, intended_edit=5,
        tolerated_edits=[], note=""
    )

    with patch("pipeline.quantify.read_mapping_stats", return_value=(100, 80)), \
         patch("pipeline.quantify.read_allele_table", return_value=pd.DataFrame()), \
         patch("pipeline.quantify.read_quant_window", return_value=pd.DataFrame()), \
         patch("pipeline.quantify.find_het_position", return_value=([], None, None)), \
         patch("pipeline.quantify.quantify_oneseq_sample", return_value={"branch": "oneseq"}) as mock_oneseq, \
         patch("pipeline.quantify.quantify_abe_sample", return_value={"branch": "abe"}) as mock_abe, \
         patch("pipeline.quantify.quantify_het_sample", return_value={"branch": "het"}) as mock_het:
        result = quantify_sample(config, sample_dir)

    mock_abe.assert_called_once()
    mock_oneseq.assert_not_called()
    mock_het.assert_not_called()
    assert result == {"branch": "abe"}


def test_dispatch_routes_ABE_het_to_het_sample(tmp_path):
    sample_dir = make_dummy_crispresso_output(tmp_path)
    config = AmpliconConfig(
        name="TEST", protospacer="A" * 20, editor="ABE",
        orientation="F", amplicon="A" * 30, intended_edit=5,
        tolerated_edits=[], note=""
    )

    with patch("pipeline.quantify.read_mapping_stats", return_value=(100, 80)), \
         patch("pipeline.quantify.read_allele_table", return_value=pd.DataFrame()), \
         patch("pipeline.quantify.read_quant_window", return_value=pd.DataFrame()), \
         patch("pipeline.quantify.find_het_position", return_value=([5], "C", "T")), \
         patch("pipeline.quantify.quantify_oneseq_sample", return_value={"branch": "oneseq"}) as mock_oneseq, \
         patch("pipeline.quantify.quantify_abe_sample", return_value={"branch": "abe"}) as mock_abe, \
         patch("pipeline.quantify.quantify_het_sample", return_value={"branch": "het"}) as mock_het:
        result = quantify_sample(config, sample_dir)

    mock_het.assert_called_once()
    mock_oneseq.assert_not_called()
    mock_abe.assert_not_called()
    assert result == {"branch": "het"}


def test_dispatch_unknown_editor_FORCED_FAIL(tmp_path):
    sample_dir = make_dummy_crispresso_output(tmp_path)
    config = AmpliconConfig(
        name="TEST", protospacer="A" * 20, editor="CBE",
        orientation="F", amplicon="A" * 30, intended_edit=5,
        tolerated_edits=[], note=""
    )

    with patch("pipeline.quantify.read_mapping_stats", return_value=(100, 80)), \
         patch("pipeline.quantify.read_allele_table", return_value=pd.DataFrame()), \
         patch("pipeline.quantify.read_quant_window", return_value=pd.DataFrame()), \
         patch("pipeline.quantify.find_het_position", return_value=([], None, None)), \
         patch("pipeline.quantify.quantify_oneseq_sample"), \
         patch("pipeline.quantify.quantify_abe_sample"), \
         patch("pipeline.quantify.quantify_het_sample"):
        with pytest.raises(ValueError):
            quantify_sample(config, sample_dir)
