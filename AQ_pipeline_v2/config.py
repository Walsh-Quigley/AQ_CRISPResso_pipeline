from dataclasses import dataclass
from pathlib import Path

# AmpliconConfig dataclass — one instance per row of amplicon_list.csv



@dataclass
class AmpliconConfig:
    """data class for an entry (row) in an amplicon_list.csv
    Attributes:
        name: name of the amplicon
        protospacer: guide sequence 
        editor: type of editing being used
        orientation: guide's orientation relative to the amplicon
        amplicon: the amplicon sequence itself
        intended_edit: the intended edit position within the protospacer, 1-indexed
        tolerated_edits: list of tolerated edit positions within the protospacer, 1-indexed
        note: optional note for user input 
    """
    name: str
    protospacer: str
    editor: str                 # "ABE", "ONESEQ"
    orientation: str            # "F" or "R"
    amplicon: str
    intended_edit: int | str         # 1-indexed position of the target base
    tolerated_edits: list[int]  # 1-indexed position of tolerated bystander edits
    note: str = ""

