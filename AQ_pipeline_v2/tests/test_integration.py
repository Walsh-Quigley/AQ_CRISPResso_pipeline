import subprocess
from pathlib import Path
from tests.helper import make_fastq_gz, make_nuclease_fastq_gz
from utils.sequences import reverse_complement
import os
import pandas as pd
import shutil

"""End-to-end integration tests for the AQ pipeline - covers ABE (forward + reverse), 
het ABE (heterozygous forward + reverse), and ONESEQ (forward + reverse) full 
workflows. Each test builds a tmp_path directory with fastq files and amplicon_list, 
then invokes CRISPResso_Loop.py and Quantification_Loop.py via subprocess. Output 
CSVs are saved to tests/test_output/ for manual inspection.

Requires: Linux or macOS environment with CRISPResso installed (Windows users 
can use WSL)."""


env = os.environ.copy()
env["PYTHONPATH"] = str(Path(__file__).parent.parent)  # points to AQ_pipeline_v2/


def test_ABE(tmp_path):
    fastq_dir_f = tmp_path / "fastqs" / "Sample_PAH1_1"
    fastq_dir_f.mkdir(parents=True)
    fastq_dir_r = tmp_path / "fastqs" / "Sample_R186W_1"
    fastq_dir_r.mkdir(parents=True)

    PAH1_proto = "TCACAGTTCGGGGGTATACA"
    R186W_proto_rc = reverse_complement("CAGCAGCCACTCAGAGTCTC")

    PAH1_amplicon_seq = "A" * 60 + PAH1_proto + "A" * 71 #numbers come from illumina 151 length
    R186W_amplicon_seq = "A" * 60 + R186W_proto_rc + "A" * 71

    with open(tmp_path / "amplicon_list.csv", "w") as f:
        f.write("name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n")
        f.write(f"PAH1,TCACAGTTCGGGGGTATACA,ABE,F,{PAH1_amplicon_seq},,3,5\n")
        f.write(f"R186W,CAGCAGCCACTCAGAGTCTC,ABE,R,{R186W_amplicon_seq},,13,9")

    make_fastq_gz(
        fastq_dir_f / "reads_R1.fastq.gz",
        protospacer=PAH1_proto,
        n_reads=1000,
        outcomes={
            "unedited": (0.4, []),
            "corrected": (0.2, [(4, "A", "G")]),
            "corrected_with_bystander": (0.2, [(4, "A", "G"), (2, "A", "G")]),
            "corrected_with_atog": (0.1, [(4, "A", "G"), (15, "A", "G")]),
            "corrected_with_other": (0.1, [(4, "A", "G"), (1, "C", "T")])
        }
    )
    make_fastq_gz(
        fastq_dir_r / "reads_R1.fastq.gz",
        protospacer=R186W_proto_rc,
        n_reads=2000,
        outcomes={
            "unedited": (0.4, []),
            "corrected": (0.2, [(11, "T", "C")]),
            "corrected_with_bystander": (0.2, [(11, "T", "C"), (7, "T", "C")]),
            "corrected_with_atog": (0.1, [(11, "T", "C"), (5, "T", "C")]),
            "corrected_with_other": (0.1, [(11, "T", "C"), (0, "G", "A")])
        }
    )

    subprocess.run(
        ["python", str(Path(__file__).parent.parent / "CRISPResso_Loop.py")],
        cwd=tmp_path,
        env=env,
        check=True
    )
    subprocess.run(
        ["python", str(Path(__file__).parent.parent / "Quantification_Loop.py")],
        cwd=tmp_path,
        env=env,
        input="y\n",
        text=True,
        check=True,
    )

    shutil.copy(
        tmp_path / "ABE_Quantification_Summary.csv",
        Path("tests/test_output") / "ABE_Quantification_Summary_test_ABE.csv"
    )
    shutil.copy(
        tmp_path / "Prism_Input.csv",
        Path("tests/test_output") / "Prism_Input_test_ABE.csv"
    )
    df = pd.read_csv(tmp_path / "ABE_Quantification_Summary.csv")
    pah1 = df[df["sample"] == "Sample_PAH1_1"].iloc[0]
    r186w = df[df["sample"] == "Sample_R186W_1"].iloc[0]
    for row in [pah1, r186w]:
        assert row["correction_without_bystanders"] == 20.0
        assert row["correction_with_tolerated_bystanders"] == 40.0
        assert row["correction_with_any_AtoG_change"] == 50.0
        assert row["correction_with_any_change_in_protospacer"] == 60.0
    prism_csv = tmp_path / "Prism_Input.csv"
    assert prism_csv.exists()
    prism_df = pd.read_csv(prism_csv)
    # both base_samples present
    base_samples = set(prism_df["base_sample"])
    assert "Sample_PAH1" in base_samples
    assert "Sample_R186W" in base_samples
    # spot-check the value pipeline (correction_without_bystanders_rep1 == 20.0 from the asserts above)
    pah1_row = prism_df[prism_df["base_sample"] == "Sample_PAH1"].iloc[0]
    assert pah1_row["correction_without_bystanders_rep1"] == 20.0

def test_het_ABE(tmp_path):
    fastq_dir_f = tmp_path / "fastqs" / "Sample_PAH1_1"
    fastq_dir_f.mkdir(parents=True)
    fastq_dir_r = tmp_path / "fastqs" / "Sample_R186W_1"
    fastq_dir_r.mkdir(parents=True)

    PAH1_proto = "TCACAGTTCGGGGGTATACA"
    R186W_proto_rc = reverse_complement("CAGCAGCCACTCAGAGTCTC")

    PAH1_amplicon_seq = "A" * 60 + PAH1_proto + "A" * 71
    R186W_amplicon_seq = "A" * 60 + R186W_proto_rc + "A" * 71

    with open(tmp_path / "amplicon_list.csv", "w") as f:
        f.write("name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n")
        f.write(f"PAH1,TCACAGTTCGGGGGTATACA,ABE,F,{PAH1_amplicon_seq},,3,5\n")
        f.write(f"R186W,CAGCAGCCACTCAGAGTCTC,ABE,R,{R186W_amplicon_seq},,13,9")
    
    make_fastq_gz(
        fastq_dir_f / "reads_R1.fastq.gz",
        protospacer="TCACAGTTCGGGGGTATACA",
        n_reads=1000,
        outcomes={
            "allele1_unedited":             (0.20, []),
            "allele1_corrected":            (0.10, [(4, "A", "G")]),
            "allele1_corrected_w_bystander":(0.10, [(4, "A", "G"), (2, "A", "G")]),
            "allele1_corrected_w_atog":     (0.05, [(4, "A", "G"), (15, "A", "G")]),
            "allele1_corrected_w_other":    (0.05, [(4, "A", "G"), (1, "C", "T")]),
            "allele2":                      (0.50, [(0, "T", "A"), (18, "C", "G"), (4, "A", "G")])
        }
    )
    make_fastq_gz(
        fastq_dir_r / "reads_R1.fastq.gz",
        protospacer=R186W_proto_rc,
        n_reads=2000,
        outcomes={
            "allele1_unedited":              (0.20, []),
            "allele1_corrected":             (0.10, [(11, "T", "C")]),
            "allele1_corrected_w_bystander": (0.10, [(11, "T", "C"), (7, "T", "C")]),
            "allele1_corrected_w_atog":      (0.05, [(11, "T", "C"), (5, "T", "C")]),
            "allele1_corrected_w_other":     (0.05, [(11, "T", "C"), (0, "G", "A")]),
            "allele2":                       (0.50, [(19, "G", "C"), (1, "A", "T"), (11, "T", "C")])
        }
    )
    subprocess.run(
        ["python", str(Path(__file__).parent.parent / "CRISPResso_Loop.py")],
        cwd=tmp_path,
        env=env,
        check=True
    )
    subprocess.run(
        ["python", str(Path(__file__).parent.parent / "Quantification_Loop.py")],
        cwd=tmp_path,
        env=env,
        input="y\n",
        text=True,
        check=True
    )
    shutil.copy(
        tmp_path / "ABE_Quantification_Summary.csv",
        Path("tests/test_output") / "ABE_Quantification_Summary_test_het_ABE.csv"
    )
    shutil.copy(
        tmp_path / "Prism_Input_het_allele1.csv",
        Path("tests/test_output") / "Prism_Input_het_allele1_test_het_ABE.csv"
    )
    shutil.copy(
        tmp_path / "Prism_Input_het_allele2.csv",
        Path("tests/test_output") / "Prism_Input_het_allele2_test_het_ABE.csv"
    )
    df = pd.read_csv(tmp_path / "ABE_Quantification_Summary.csv")
    pah1_het = df[df["sample"] == "Sample_PAH1_1"].iloc[0]
    r186w_het = df[df["sample"] == "Sample_R186W_1"].iloc[0]
    for row in [pah1_het, r186w_het]:
        assert row["correction_wo_bystanders_allele1"] == 20.0
        assert row["correction_w_bystanders_allele1"] == 40.0
        assert row["correction_with_any_AtoG_change_allele1"] == 50.0
        assert row["correction_with_any_change_in_protospacer_allele1"] == 60.0
        assert row["correction_wo_bystanders_allele2"] == 0.0
        assert row["correction_w_bystanders_allele2"] == 0.0
        assert row["correction_with_any_AtoG_change_allele2"] == 100.0
        assert row["correction_with_any_change_in_protospacer_allele2"] == 100.0
    
    allele1_csv = tmp_path / "Prism_Input_het_allele1.csv"
    allele2_csv = tmp_path / "Prism_Input_het_allele2.csv"
    assert allele1_csv.exists()
    assert allele2_csv.exists()

    allele1_df = pd.read_csv(allele1_csv)
    allele2_df = pd.read_csv(allele2_csv)

    # both base_samples present in both files
    for df in [allele1_df, allele2_df]:
        base_samples = set(df["base_sample"])
        assert "Sample_PAH1" in base_samples
        assert "Sample_R186W" in base_samples

    # spot-check values match the per-allele assertions above
    pah1_a1 = allele1_df[allele1_df["base_sample"] == "Sample_PAH1"].iloc[0]
    assert pah1_a1["correction_wo_bystanders_rep1"] == 20.0   # matches allele1 wo_bystanders assert above

    pah1_a2 = allele2_df[allele2_df["base_sample"] == "Sample_PAH1"].iloc[0]
    assert pah1_a2["correction_wo_bystanders_rep1"] == 0.0    # matches allele2 wo_bystanders assert above

    # per-allele reads_aligned should be ~half of total since het is 50/50
    assert pah1_a1["reads_aligned_rep1"] == 500
    assert pah1_a2["reads_aligned_rep1"] == 500

    # R186W
    r186w_a1 = allele1_df[allele1_df["base_sample"] == "Sample_R186W"].iloc[0]
    r186w_a2 = allele2_df[allele2_df["base_sample"] == "Sample_R186W"].iloc[0]
    assert r186w_a1["reads_aligned_rep1"] == 1000
    assert r186w_a2["reads_aligned_rep1"] == 1000

def test_ONESEQ(tmp_path):
    fastq_dir_f= tmp_path / "fastqs" / "Sample_forward_ONESEQ_1"
    fastq_dir_f.mkdir(parents=True)
    fastq_dir_r= tmp_path / "fastqs" / "Sample_reverse_ONESEQ_1"
    fastq_dir_r.mkdir(parents=True)

    forward_proto = "AACAGTTTCGGGGGTATACA"
    reverse_proto_rc = reverse_complement("TTCAGTTTCGAGGGTATACA")

    forward_example_amplicon_seq = "A" * 60 + forward_proto + "A" * 71
    reverse_example_amplicon_seq = "A" * 60 + reverse_proto_rc + "A" * 71

    with open(tmp_path / "amplicon_list.csv", "w") as f:
        f.write("name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n")
        f.write(f"forward,AACAGTTTCGGGGGTATACA,ONESEQ,F,{forward_example_amplicon_seq},,,ONE-SEQ\n")
        f.write(f"reverse,TTCAGTTTCGAGGGTATACA,ONESEQ,R,{reverse_example_amplicon_seq},,,ONE-SEQ")
        
    make_fastq_gz(
        fastq_dir_f / "reads_R1.fastq.gz",
        protospacer="AACAGTTTCGGGGGTATACA",
        n_reads=1000,
        outcomes={
            "unedited":       (0.4, []),
            "first10":        (0.4, [(0, "A", "G")]),
            "anywhere":       (0.2, [(15, "A", "G")])
        }
    )
    make_fastq_gz(
        fastq_dir_r / "reads_R1.fastq.gz",
        protospacer="TGTATACCCTCGAAACTGAA",
        n_reads=1000,
        outcomes={
            "unedited":       (0.4, []),
            "first10":        (0.4, [(16, "T", "C")]),
            "anywhere":       (0.2, [(9, "T", "C")])
        }
    )
    subprocess.run(
        ["python", str(Path(__file__).parent.parent / "CRISPResso_Loop.py")],
        cwd=tmp_path,
        env=env,
        check=True
    )
    subprocess.run(
        ["python", str(Path(__file__).parent.parent / "Quantification_Loop.py")],
        cwd=tmp_path,
        env=env,
        input="n\n",
        text=True,
        check=True
    )
    shutil.copy(
        tmp_path / "ONESEQ_Quantification_Summary.csv",
        Path("tests/test_output") / "ONESEQ_Quantification_Summary_test_ONESEQ.csv"
    )
    df = pd.read_csv(tmp_path / "ONESEQ_Quantification_Summary.csv")
    forward_ONESEQ = df[df["sample"] == "Sample_forward_ONESEQ_1"].iloc[0]
    reverse_ONESEQ = df[df["sample"] == "Sample_reverse_ONESEQ_1"].iloc[0]
    for row in [forward_ONESEQ, reverse_ONESEQ]:
        assert row["pct_AtoG_first_10bp"] == 40.0
        assert row["pct_AtoG_anywhere"] == 60.0

def test_NUCLEASE(tmp_path):
    fastq_dir = tmp_path / "fastqs" / "Sample_KO1_1"
    fastq_dir.mkdir(parents=True)

    KO1_proto = "GCATGACTAGTCGTACGCTG"
    KO1_amplicon_seq = "A" * 60 + KO1_proto + "A" * 71

    with open(tmp_path / "amplicon_list.csv", "w") as f:
        f.write("name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n")
        f.write(f"KO1,{KO1_proto},NUCLEASE,F,{KO1_amplicon_seq},,,\n")

    make_nuclease_fastq_gz(
        fastq_dir / "reads_R1.fastq.gz",
        protospacer=KO1_proto,
        n_reads=1000,
        outcomes={
            "unmodified": (0.40, []),
            "del1":       (0.20, [("del", 17, 1)]),
            "del3":       (0.20, [("del", 16, 3)]),
            "ins2":       (0.20, [("ins", 17, "TT")]),
        },
    )

    subprocess.run(
        ["python", str(Path(__file__).parent.parent / "CRISPResso_Loop.py")],
        cwd=tmp_path,
        env=env,
        check=True,
    )
    subprocess.run(
        ["python", str(Path(__file__).parent.parent / "Quantification_Loop.py")],
        cwd=tmp_path,
        env=env,
        input="n\n",
        text=True,
        check=True,
    )

    shutil.copy(
        tmp_path / "NUCLEASE_Quantification_Summary.csv",
        Path("tests/test_output") / "NUCLEASE_Quantification_Summary_test_NUCLEASE.csv",
    )

    df = pd.read_csv(tmp_path / "NUCLEASE_Quantification_Summary.csv")
    ko1 = df[df["sample"] == "Sample_KO1_1"].iloc[0]

    assert ko1["pct_unmodified"] == 40.0
    assert ko1["pct_modified"] == 60.0
    assert ko1["pct_deletions"] == 40.0
    assert ko1["pct_insertions"] == 20.0
    assert ko1["pct_frameshift_indels"] == 40.0
    assert ko1["pct_inframe_indels"] == 20.0

def test_failed_samples_csv_generated(tmp_path):
    # Create fastqs dir with a sample whose name won't match any amplicon
    fastq_dir = tmp_path / "fastqs" / "UnmatchableName_1"
    fastq_dir.mkdir(parents=True)

    # Create an amplicon_list that doesn't contain "UnmatchableName"
    PAH1_proto = "TCACAGTTCGGGGGTATACA"
    PAH1_amplicon_seq = "A" * 60 + PAH1_proto + "A" * 71
    with open(tmp_path / "amplicon_list.csv", "w") as f:
        f.write("name,protospacer_or_PEG,editor,guide_orientation_relative_to_amplicon,amplicon,note,tolerated_edits,intended_edit\n")
        f.write(f"PAH1,{PAH1_proto},ABE,F,{PAH1_amplicon_seq},,3,5\n")

    # Run Quantification_Loop only — no CRISPResso needed because
    # the failure is at the identify_amplicon step, before CRISPResso output is read
    subprocess.run(
        ["python", str(Path(__file__).parent.parent / "Quantification_Loop.py")],
        cwd=tmp_path,
        env=env,
        input="n\n",
        text=True,
        check=True,
    )

    # Verify failed_samples.csv was created with the expected content
    failed_csv = tmp_path / "failed_samples.csv"
    assert failed_csv.exists()

    df = pd.read_csv(failed_csv)
    assert len(df) == 1
    assert df["sample"].iloc[0] == "UnmatchableName_1"
    assert df["error_type"].iloc[0] == "ValueError"
    assert "amplicon" in df["error_message"].iloc[0].lower()
    assert not (tmp_path / "ABE_Quantification_Summary.csv").exists()
