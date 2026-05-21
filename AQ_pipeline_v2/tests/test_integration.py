import subprocess
from pathlib import Path
from tests.helper import make_fastq_gz
from utils.sequences import reverse_complement
import os
import pandas as pd
import shutil

env = os.environ.copy()
env["PYTHONPATH"] = str(Path(__file__).parent.parent)  # points to AQ_pipeline_v2/


def test_ABE(tmp_path):
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
        input="n\n",
        text=True,
        check=True
    )
    shutil.copy(
        tmp_path / "ABE_Quantification_Summary.csv",
        Path("tests/test_output") / "ABE_Quantification_Summary_test_ABE.csv"
    )
    df = pd.read_csv(tmp_path / "ABE_Quantification_Summary.csv")
    pah1 = df[df["sample"] == "Sample_PAH1_1"].iloc[0]
    r186w = df[df["sample"] == "Sample_R186W_1"].iloc[0]
    for row in [pah1, r186w]:
        assert row["correction_without_bystanders"] == 20.0
        assert row["correction_with_tolerated_bystanders"] == 40.0
        assert row["correction_with_any_AtoG_change"] == 50.0
        assert row["correction_with_any_change_in_protospacer"] == 60.0

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
        input="n\n",
        text=True,
        check=True
    )
    shutil.copy(
        tmp_path / "ABE_Quantification_Summary.csv",
        Path("tests/test_output") / "ABE_Quantification_Summary_test_het_ABE.csv"
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