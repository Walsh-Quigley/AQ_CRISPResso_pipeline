# utils/sequences.py
from itertools import combinations


# DNA sequence utilities used across the pipeline.
# (intended edit + tolerated bystander combinations) for allele table filtering



# takes in a passed through str sequence and 
# returns the reverse complement of that sequence
# as a string
def reverse_complement(seq:str) -> str:
    complement = {"A": "T", 
                  "T" : "A",
                  "C" : "G",
                  "G" : "C" }
    return "".join(complement[base] for base in reversed(seq.upper()))
    
def generate_search_sequences(
                            protospacer: str,
                            intended_edit: int,
                            tolerated_edits: list[int],
                            orientation: str,
                        ) -> list[str]:
    
    if orientation not in ("F", "R"):
        raise ValueError(f"orientation must be 'F' or 'R', got '{orientation}'")
    bases = list(protospacer.upper())

    target_base = "A"

    if bases[intended_edit - 1] != target_base:
        raise ValueError(
            f"Expected '{target_base}' at position {intended_edit} but found '{bases[intended_edit - 1]}'. "
            f"Check intended_edit in amplicon_list.csv."
        )
    for pos in tolerated_edits:
        if bases[pos - 1] != target_base:
            raise ValueError(
                f"Expected '{target_base}' at tolerated edit position {pos} but found '{bases[pos - 1]}'. "
                f"Check tolerated_edits in amplicon_list.csv."
            )

        
    bases[intended_edit - 1] = "G"

    sequences = []
    for size in range(len(tolerated_edits) + 1):
        for combo in combinations(tolerated_edits, size):
            current = bases.copy()
            for pos in combo:
                current[pos - 1] = "G"
            sequences.append("".join(current))
        
    if orientation == "R":
        sequences = [reverse_complement(s) for s in sequences]

    return sequences
