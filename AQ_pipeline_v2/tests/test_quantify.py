import pytest
from pathlib import Path
from pipeline.quantify import quantify_sample
from config import AmpliconConfig

"""Tests for pipeline/quantify.py - covers quantify sample, missing CRISPResso_subfolder
forced fail, and multiple allele frequency tables forced fail"""


def test_quantify_sample(tmp_path):
    CRISPResso_mapping_statistic = tmp_path / "sample_dir" / "CRISPResso_on_sample" / "CRISPResso_mapping_statistics.txt" 
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
    result = quantify_sample(configs[0], tmp_path/"sample_dir")
    assert result["sample"] == "sample_dir"
    assert result["reads_total"] == 6500
    assert result["reads_aligned"] == 6000
    assert result["correction_without_bystanders"] == 15
    assert result["correction_with_tolerated_bystanders"] == 28.33
    assert result["correction_with_any_AtoG_change"] == 40
    assert result["correction_with_any_change_in_protospacer"] == 50
    assert result["column E minus column D"] == 13.33
    assert result["column F minus column E"] == 11.67
    assert result["column G minus column F"] == 10
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

