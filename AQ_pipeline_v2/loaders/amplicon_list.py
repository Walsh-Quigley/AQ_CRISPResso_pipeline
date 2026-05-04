import csv
import logging
from pathlib import Path
from config import AmpliconConfig

# io/amplicon_list.py
# Reads and validates amplicon_list.sv, returns list of AmpliconConfig objects
# One AmpliconConfig is created per row. Raises ValueError on missing or invalid data

def load_amplicon_list(path: Path) -> list[AmpliconConfig]:
    """Parses amplicon list file and returns a list of the custom made 
    AmpliconConfig objects.
    Args:
        path: A path to the amplicon_list.csv
    Returns:
        list[AmpliconConfig]: a list of AmpliconConfig objects
    Raises:
        ValueError: if there are too many columns in the amplicon_list.csv
        ValueError: if there is a non-digit in the intended_edit column (excluding "ONESEQ")
    """
    configs = []

    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"].strip()

            if row.get(None):
                raise ValueError(
                    f"Row for '{name}' has too many columns — this usually means unquoted commas in tolerated_edits. "
                    f"Use spaces instead (e.g. '4 7 15')."
                )
            
            protospacer = row["protospacer_or_PEG"].strip().upper()
            editor = row["editor"].strip().upper()
            orientation = row["guide_orientation_relative_to_amplicon"].strip().upper()
            amplicon = row["amplicon"].strip().upper()
            note = row.get("note","").strip()

            if len(amplicon) < len(protospacer):
                logging.warning(f"Amplicon for '{name}' is shorter than its protospacer — check your amplicon_list.")

            intended_edit_raw = row["intended_edit"].strip().upper().replace("-", "")
            if intended_edit_raw == "ONESEQ":
                intended_edit = "ONESEQ"
            elif intended_edit_raw.isdigit():
                intended_edit = int(intended_edit_raw)
            else:
                raise ValueError(f"Invalid intended_edit value '{row['intended_edit']}' for amplicon '{name}'")

            tolerated_edits_raw = row.get("tolerated_edits", "").strip()
            if tolerated_edits_raw:
                tolerated_edits = [int(x) for x in tolerated_edits_raw.replace(",", " ").split()]
            else:
                tolerated_edits = []


            configs.append(AmpliconConfig(
                name=name,
                protospacer=protospacer,
                editor=editor,
                orientation=orientation,
                amplicon=amplicon,
                intended_edit=intended_edit,
                tolerated_edits=tolerated_edits,
                note=note,
            ))
    return configs

def find_amplicon_list(search_dir: Path = Path(".")) -> Path:
    """Parse through user's project direcotry for potential amplicon_list.csv
    Args:
        search_dir: the directory path to search within, set to the current directory by default
    Returns:
        Path: a path to the amplicon_list.csv itself
    Raises:
        FileNotFoundError: if no file contianing 'amplicon_list' is found.
        ValueError: if multiple files containing 'amplicon_list" are found.
    """
    found_lists = 0
    found_names = []

    for file in search_dir.iterdir():
        if "amplicon_list" in file.name.lower() and file.is_file():
            found_lists += 1
            found_names.append(file.name)
            to_be_returned = file
    if found_lists == 0:
        raise FileNotFoundError(f"No file containing 'amplicon_list' found in {Path('.').resolve()}")
    elif found_lists == 1:
        return to_be_returned
    elif found_lists > 1:
        raise ValueError(f"Multiple amplicon list files found: {', '.join(found_names)}")
    else:
        print("it should be impossible for you to see this message")
