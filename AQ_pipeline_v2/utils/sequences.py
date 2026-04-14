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
    base = list(protospacer.upper())
        
    base[intended_edit - 1] = "G"

    sequences = []
    for size in range(len(tolerated_edits) + 1):
        for combo in combinations(tolerated_edits, size):
            current = base.copy()
            for pos in combo:
                current[pos - 1] = "G"
            sequences.append("".join(current))
        
    if orientation == "R":
        sequences = [reverse_complement(s) for s in sequences]

    return sequences
