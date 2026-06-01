import csv
import logging
from pathlib import Path
from config import AmpliconConfig
from utils.sequences import reverse_complement


# io/amplicon_list.py
# Reads and validates amplicon_list.csv, returns list of AmpliconConfig objects
# One AmpliconConfig is created per row. Raises ValueError on missing or invalid data

def load_amplicon_list(path: Path) -> list[AmpliconConfig]:
    """Parses amplicon list file and returns a list of the custom made 
    AmpliconConfig objects.
    Args:
        path: A path to the amplicon_list.csv
    Returns:
        list[AmpliconConfig]: a list of AmpliconConfig objects
    Raises:
        ValueError: if editor is ONESEQ but intended_edit is a numeric position
    """
    configs = []

    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        seen_names = set()
        for row in reader:
            name = row["name"].strip()

            if name in seen_names:
                raise ValueError(f"Duplicate amplicon name '{name}' in {path}")
            seen_names.add(name)

            if row.get(None):
                raise ValueError(
                    f"Row for '{name}' has too many columns — this usually means unquoted commas in tolerated_edits. "
                    f"Use spaces instead (e.g. '4 7 15')."
                )
            
            protospacer = row["protospacer_or_PEG"].strip().upper()
            editor = row["editor"].strip().upper()
            orientation = row["guide_orientation_relative_to_amplicon"].strip().upper()
            if orientation not in ("F", "R"):
                raise ValueError(f"Orientation for '{name}' must be 'F' or 'R', got '{orientation}'")
            amplicon = row["amplicon"].strip().upper()
            note = row.get("note","").strip()

            invalid_protospacer = set(protospacer) - {"A","C","G","T"}
            if invalid_protospacer:
                raise ValueError(f"Protospacer for '{name}' contains non-DNA characters: {invalid_protospacer}")

            invalid_amplicon = set(amplicon) - {"A","C","G","T"}
            if invalid_amplicon:
                raise ValueError(f"Amplicon for '{name}' contains non-DNA characters: {invalid_amplicon}")
            
            if protospacer not in amplicon and reverse_complement(protospacer) not in amplicon:
                raise ValueError(f"Protospacer for '{name}' (or its reverse complement) not found in the amplicon sequence — check for swapped rows or typos.")

            if len(amplicon) < len(protospacer):
                logging.warning(f"Amplicon for '{name}' is shorter than its protospacer — check your amplicon_list.")

            intended_edit_raw = row["intended_edit"].strip().upper().replace("-", "")
            if intended_edit_raw == "ONESEQ":
                intended_edit = "ONESEQ"
            elif intended_edit_raw.isdigit():
                intended_edit = int(intended_edit_raw)
            else:
                raise ValueError(f"Invalid intended_edit value '{row['intended_edit']}' for amplicon '{name}'")
            if editor == "ONESEQ" and intended_edit != "ONESEQ":
                raise ValueError(f"Amplicon '{name}' has editor=ONESEQ but intended_edit is '{intended_edit}' — set intended_edit to ONE-SEQ.")


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
    """Parse through user's project directory for potential amplicon_list.csv
    Args:
        search_dir: the directory path to search within, set to the current directory by default
    Returns:
        Path: a path to the amplicon_list.csv itself
    Raises:
        FileNotFoundError: if no file containing 'amplicon_list' is found.
        ValueError: if multiple files containing 'amplicon_list' are found.
    """
    found_lists = 0
    found_names = []

    for file in search_dir.iterdir():
        if "amplicon_list" in file.name.lower() and file.is_file():
            found_lists += 1
            found_names.append(file.name)
            to_be_returned = file
    if found_lists == 0:
        raise FileNotFoundError(f"No file containing 'amplicon_list' found in {search_dir.resolve()}")
    elif found_lists == 1:
        return to_be_returned
    elif found_lists > 1:
        raise ValueError(f"Multiple amplicon list files found: {', '.join(found_names)}")
