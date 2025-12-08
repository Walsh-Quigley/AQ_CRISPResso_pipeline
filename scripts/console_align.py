#scripts/console_align.py

def console_align(seq1, seq2):
    """
    Align two sequences in the console and show matches with | and mismatches with space.
    seq1: reference (amplicon region)
    seq2: query (RTT + guide)
    """
    # Make sure sequences are uppercase
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    # Pad the shorter sequence
    if len(seq1) < len(seq2):
        seq1 += "-" * (len(seq2) - len(seq1))
    elif len(seq2) < len(seq1):
        seq2 += "-" * (len(seq1) - len(seq2))

    # Build match line and record mismatches
    match_line = []
    mismatches = []
    for i, (a, b) in enumerate(zip(seq1, seq2)):
        if a == b:
            match_line.append("|")
        else:
            match_line.append("*")
            mismatches.append(i)  # store mismatch position

    match_line_str = "".join(match_line)



    return match_line_str, mismatches