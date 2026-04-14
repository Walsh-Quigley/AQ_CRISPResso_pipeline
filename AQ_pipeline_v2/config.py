from dataclasses import dataclass, field
from pathlib import Path

# Dataclasses that hold all settings: amplicon info, threshholds, 
# file paths, editor type


# data class for the amplicon list entries
# note section is set witha  default value of 
# an empty string so that the cell can be left 
# blank and not cause an error
@dataclass
class AmpliconConfig:
    name: str
    protospacer: str
    editor: str                 # "ABE", "ONE_seq"
    orientation: str            # "F" or "R"
    amplicon: str
    intended_edit: int          # 1-indexed position of the target base
    tolerated_edits: list[int]  # 1-indexed position of tolerated bystander edits
    note: str = ""


# data class for configuration. needs a path to the fastqs
# het_min and het_max are set but can be changed by user later. 
# prism ouput does nto generate automatically but a user can elect to
# output file name remains quantification summary but can be changed
# by user
@dataclass
class RunConfig:
    fastqs_dir: Path
    het_min_freq: float = 0.40
    het_max_freq: float = 0.60
    generate_prism_output: bool = False
    output_filename: str = "quantification_summary.csv"
