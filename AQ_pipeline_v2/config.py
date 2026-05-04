from dataclasses import dataclass
from pathlib import Path

# Dataclasses that hold all settings: amplicon info, threshholds, 
# file paths, editor type



@dataclass
class AmpliconConfig:
    """data class for an entry(row) in an amplicon_list.csv
    Attributes:
        name: name of the amplicon
        protospacer: guide sequence 
        editior: type of editing being used
        orientation: guides oriention relative to the amplicon
        amplicon: the amplicon sequence itself
        intended_edit: the intended edit position within the protospacer, 1-indexed
        tolerated_edits: list of tolerated edit positions with the protospacer, 1-indexed
        note: optional note for user input 
    """
    name: str
    protospacer: str
    editor: str                 # "ABE", "ONE_seq"
    orientation: str            # "F" or "R"
    amplicon: str
    intended_edit: int          # 1-indexed position of the target base
    tolerated_edits: list[int]  # 1-indexed position of tolerated bystander edits
    note: str = ""


@dataclass
class RunConfig:
    """
    """
    fastqs_dir: Path
    het_min_freq: float = 0.40
    het_max_freq: float = 0.60
    generate_prism_output: bool = False
    output_filename: str = "Quantification_Summary.csv"
