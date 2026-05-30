import sys
import csv
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "AQ_pipeline_v2"))

from tests.helper import make_fastq_gz
from utils.sequences import reverse_complement

SAMPLES = [
    {
        "name": "hetF_1",
        "guide_seq": "CCACACCGGGTCCCTGCGCT",
        "orientation": "F",
        "intended_edit": 5,
        "tolerated_edits": [3],
        "het_pairs": [(9, "G", "T")],      # (guide 0-idx, allele1 base, allele2 base)
    },
    {
        "name": "hetF_2",
        "guide_seq": "CCTGGCACGTGCCTCACGCT",
        "orientation": "F",
        "intended_edit": 7,
        "tolerated_edits": [],
        "het_pairs": [(3, "G", "C"), (15, "A", "T")],
    },
    {
        "name": "hetF_3",
        "guide_seq": "GCAGACCATCGCGTCGGTCC",
        "orientation": "F",
        "intended_edit": 5,
        "tolerated_edits": [3, 8],
        "het_pairs": [(11, "C", "G"), (16, "G", "T")],
    },
    {
        "name": "hetR_1",
        "guide_seq": "GCTGCAGCGTCGCATGCTGC",
        "orientation": "R",
        "intended_edit": 6,
        "tolerated_edits": [14],
        "het_pairs": [(10, "C", "G")],
    },
    {
        "name": "hetR_2",
        "guide_seq": "GCTGGCTGACTGCGTCGTGC",
        "orientation": "R",
        "intended_edit": 9,
        "tolerated_edits": [],
        "het_pairs": [(3, "G", "T"), (15, "C", "G")],
    },
]

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

def make_outcomes(guide_seq, orientation, intended_edit, tolerated_edits, het_pairs):
    n = len(guide_seq)

    if orientation == "F":
        proto = guide_seq
        intended_tuple = (intended_edit - 1, "A", "G")
        tolerated_tuples = [(t - 1, "A", "G") for t in tolerated_edits]
        allele2_het_edits = [(idx, b1, b2) for (idx, b1, b2) in het_pairs]
    else:
        proto = reverse_complement(guide_seq)
        intended_tuple = ((n - 1) - (intended_edit - 1), "T", "C")
        tolerated_tuples = [((n - 1) - (t - 1), "T", "C") for t in tolerated_edits]
        allele2_het_edits = [
            ((n - 1) - idx, COMP[b1], COMP[b2])
            for (idx, b1, b2) in het_pairs
        ]

    if tolerated_edits:
        outcomes = {
            "allele1_unedited":              (0.20, []),
            "allele1_corrected":             (0.10, [intended_tuple]),
            "allele1_corrected_w_bystander": (0.20, [intended_tuple] + tolerated_tuples),
            "allele2_unedited":              (0.20, allele2_het_edits),
            "allele2_corrected":             (0.10, [intended_tuple] + allele2_het_edits),
            "allele2_corrected_w_bystander": (0.20, [intended_tuple] + tolerated_tuples + allele2_het_edits),
        }
    else:
        outcomes = {
            "allele1_unedited":  (0.25, []),
            "allele1_corrected": (0.25, [intended_tuple]),
            "allele2_unedited":  (0.25, allele2_het_edits),
            "allele2_corrected": (0.25, [intended_tuple] + allele2_het_edits),
        }

    return proto, outcomes

def main():
    FASTQS_DIR = Path("fastqs")
    N_READS = 2000
    amplicon_rows = []

    for sample in SAMPLES:
        proto, outcomes = make_outcomes(
            sample["guide_seq"], sample["orientation"],
            sample["intended_edit"], sample["tolerated_edits"], sample["het_pairs"]
        )

        amplicon = "A" * 60 + proto + "A" * (151 - 60 - len(proto))
        fastq_dir = FASTQS_DIR / f"{sample['name']}_1"
        fastq_dir.mkdir(parents=True, exist_ok=True)

        make_fastq_gz(fastq_dir / "reads_R1.fastq.gz", proto, N_READS, outcomes)

        amplicon_rows.append({
            "name":                                 sample["name"],
            "protospacer_or_PEG":                   sample["guide_seq"],
            "editor":                               "ABE",
            "guide_orientation_relative_to_amplicon": sample["orientation"],
            "amplicon":                             amplicon,
            "note":                                 "",
            "tolerated_edits":                      " ".join(str(t) for t in sample["tolerated_edits"]),
            "intended_edit":                        sample["intended_edit"],
        })
        print(f"Generated {sample['name']}: {N_READS} reads, {len(outcomes)} outcomes")

    with open("het_amplicon_list.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "name", "protospacer_or_PEG", "editor",
            "guide_orientation_relative_to_amplicon",
            "amplicon", "note", "tolerated_edits", "intended_edit"
        ])
        writer.writeheader()
        writer.writerows(amplicon_rows)

    print("\nWrote het_amplicon_list.csv")


if __name__ == "__main__":
    main()