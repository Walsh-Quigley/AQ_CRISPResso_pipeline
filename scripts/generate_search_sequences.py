# scripts/generate_search_sequences.py
import itertools
import os
from .console_align import console_align 
from .reverse_complement import reverse_complement



def generate_tolerated_sequences(edit, point_correction, corrected_sequence, guide_seq, tolerated_edits):
    print("\nGenerating tolerated sequences with bystander edits...")
    tolerated_sequences = []

    # Generate all combinations of tolerated edits (from 1 to all)
    for r in range(1, len(tolerated_edits)+1):
        for combo in itertools.combinations(tolerated_edits, r):
            seq = list(corrected_sequence)  # convert to list for mutability
            valid_combo = True

            for edit_loc in combo:
                if seq[edit_loc-1] == edit:
                    seq[edit_loc-1] = point_correction
                else:
                    print(f"\033[4mWARNING:\033[0m Tolerated edit at position {edit_loc} is not an {edit} in the corrected sequence. Skipping this combination.")
                    valid_combo = False
                    break  # skip this combo if any position is invalid

            if valid_combo:
                tolerated_seq = "".join(seq)
                tolerated_sequences.append(tolerated_seq)

    # Align each tolerated sequence to the guide
    for seq in tolerated_sequences:
        match_line_str, mismatches = console_align(guide_seq, seq)
        print(f"                             original locus: {guide_seq}")
        print(f"                                             {match_line_str}")
        print(f"corrected locus with tolerated bystander(s): {seq}")

    return tolerated_sequences

def generate_search_sequences(guide_seq, orientation, editor, intended_edit, tolerated_edits, directory_path):

    if orientation not in ["F", "R"]:
        print("\033[4mERROR:\033[0m Orientation must be 'F' for forward or 'R' for reverse.")
        return
    
    if orientation == "R":
        edit = "T"
        point_correction = "C"
        guide_seq = reverse_complement(guide_seq)
        intended_edit = len(guide_seq) - intended_edit + 1
        tolerated_edits = [len(guide_seq) - i + 1 for i in tolerated_edits]

    if orientation == "F":
        edit = "A"
        point_correction = "G"

    if not guide_seq[intended_edit-1] == edit:
        print("\033[4mERROR:\033[0m CorrectionLocationIndex is empty, not a valid integer, or not in the range 0 to 19. i.e. position 1 to 20")

    for i in tolerated_edits:
        if not guide_seq[i-1] == edit:
            print(f"\033[4mERROR:\033[0m One of the ToleratedEditIndices is not an {edit}, not a valid integer, empty or not in the range 0 to 19. i.e. position 1 to 20")
            return


    corrected_sequence = guide_seq[:intended_edit-1] + point_correction + guide_seq[intended_edit:]


    match_line_str, mismatches = console_align(guide_seq, corrected_sequence)

    print(f"original target locus:  {guide_seq}")
    print(f"                        {match_line_str}")
    print(f"corrected target locus: {corrected_sequence}")

    tolerated_sequences = generate_tolerated_sequences(edit, point_correction, corrected_sequence, guide_seq, tolerated_edits)

    search_strings = [corrected_sequence] + tolerated_sequences

    print(f"Search sequences generated: {search_strings}")

    return search_strings

def A_to_G_sequences(guide_seq, orientation):

    search_strings_first_10bp = []
    search_strings_anyA = []
    base_to_correct = "A"
    guide_seq = guide_seq.upper()
    correction = "G"

    if orientation == "R":
        guide_seq = reverse_complement(guide_seq)
        base_to_correct = "T"
        correction = "C"

    # Find positions of the base to correct
    all_edit_positions = [i for i, base in enumerate(guide_seq) if base == base_to_correct]


    # Restrict to first/last 10 positions
    if orientation == "F":
        subset_seq = guide_seq[:10]
        offset = 0  # positions relative to full guide
    elif orientation == "R":
        subset_seq = guide_seq[-10:]
        offset = len(guide_seq) - 10
    else:
        subset_seq = guide_seq
        offset = 0

    first_10_edit_positions = [i + offset for i, base in enumerate(subset_seq) if base == base_to_correct]

    # Generate all combinations of A-to-G edits in fist 10bp
    for r in range(1, len(first_10_edit_positions) + 1):
        for combo in itertools.combinations(first_10_edit_positions, r):
            seq_list = list(guide_seq)
            for pos in combo:
                seq_list[pos] = correction
            search_strings_first_10bp.append("".join(seq_list))
    
    # Generate all combinations of edits for ANY A position
    for r in range(1, len(all_edit_positions) + 1):
        for combo in itertools.combinations(all_edit_positions, r):
            seq_list = list(guide_seq)
            for pos in combo:
                seq_list[pos] = correction
            search_strings_anyA.append("".join(seq_list))

    return search_strings_first_10bp, search_strings_anyA 

