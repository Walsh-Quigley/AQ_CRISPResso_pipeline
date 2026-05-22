# utils/sequences.py
from itertools import combinations


# DNA sequence utilities used across the pipeline.
# (intended edit + tolerated bystander combinations) for allele table filtering



def reverse_complement(seq:str) -> str:
    """Reverse_complements a given string
    Args:
        seq: a sequence of characters in str format
    Returns:
        str: the reverse complement of the inputted string
    Raises:
        KeyError: if the seq contains a base that is not in our dictionary
    """
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
    """Generate all ABE search sequences based on users tolerated and intended edits
    Args:
        protospacer: the user's guide sequence
        intended_edit: the user's intended_edit location
        tolerated_edits: the user's tolerated_edit locations
        orientation: the orientaion of the guide relative to the amplicon
    Returns:
        list[str]: a list of search sequences to be used
    Raises:
        ValueError: if the orientaion is neither forward nor reversed
        ValueError: if the intended edit target base is not an A
        ValueError: if the tolerated edit location is not an A
    Note: 
        Currently ABE specific, this will need to be changed when opening up to other editing
        strategies.
    """
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

def generate_oneseq_search_sequences(protospacer:str,
                                     orientation:str) -> tuple[list[str], list[str]]:
    """Creates all search sequences for ONE_seq case. Includes all A to G changes, and A to G changes
        only in the first 10 bp.
    Args: 
        protospacer: the users guide sequence
        orientation: the orientation of the guide sequence relative to the amplicon
    Returns:
        tuple[list[str], list[str]]: returns a tuple of lists containing sequences with edits in
            the first 10bp and edits anywhere in the protospacer respectivly 
    Raises:
        ValueError: orientation is neither forward nor reverse
    Note:
        For R orientation, edited sequences are reverse-complemented before being returned
        so they match the amplicon-strand sequences reported by CRISPResso's allele table.
        An A->G edit on the guide strand appears as T->C in the allele table.
    """
    
    if orientation not in ("F", "R"):
        raise ValueError(f"Could not determine whether forward or reverse orientation")
    
    working_seq = protospacer
    
    a_pos_full_protospacer = []
    for index, c in enumerate(working_seq):
        if c.upper() == "A":
            a_pos_full_protospacer.append(index)

    a_pos_first_ten_of_protospacer = []

    for index, c in enumerate(working_seq[:10]):
        if c.upper() == "A":
            a_pos_first_ten_of_protospacer.append(index)
    first_10bp_seqs = []
    for i in range(1, len(a_pos_first_ten_of_protospacer)+1):
        for combo in combinations(a_pos_first_ten_of_protospacer, i):
            seq = list(working_seq)
            for pos in combo:
                seq[pos] = "G"
            first_10bp_seqs.append("".join(seq))
    full_bp_seqs = []
    for i in range(1, len(a_pos_full_protospacer) + 1):
        for combo in combinations(a_pos_full_protospacer, i):
            seq = list(working_seq)
            for pos in combo:
                seq[pos] = "G"
            full_bp_seqs.append("".join(seq))
    
    if orientation == "R":
        for index, seq in enumerate(first_10bp_seqs):
            first_10bp_seqs[index] = reverse_complement(seq)
        for index, seq in enumerate(full_bp_seqs):
            full_bp_seqs[index] = reverse_complement(seq)

        
    return (first_10bp_seqs, full_bp_seqs)