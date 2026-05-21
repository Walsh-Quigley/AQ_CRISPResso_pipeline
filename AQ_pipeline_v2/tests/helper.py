import pandas as pd
from utils.sequences import reverse_complement
from pathlib import Path
import gzip
from utils.sequences import reverse_complement

def make_quant_window(protospacer: str, 
                      orientation: str,
                      intended_edit: tuple = None,
                      tolerated_edits: list[tuple] = None,
                      het_positions: list[tuple] = None) -> pd.DataFrame:
    

    if orientation == "R":
        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        if intended_edit is not None:
            intended_edit = (len(protospacer) -  1 - intended_edit[0], 
                            complement[intended_edit[1]], 
                            complement[intended_edit[2]])
        if tolerated_edits is not None:
            flipped = []
            for t in tolerated_edits:
                flipped.append((len(protospacer) -  1 - t[0],
                                complement[t[1]],
                                complement[t[2]]))
            tolerated_edits = flipped
        if het_positions is not None:
            flipped = []
            for h in het_positions:
                flipped.append((len(protospacer) -  1 - h[0],
                                complement[h[1]],
                                complement[h[2]]))
            het_positions = flipped


        protospacer = reverse_complement(protospacer)
    base_to_idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    data = {}
    for idx, c in enumerate(protospacer):
        tmp_list = [0.02, 0.02, 0.02, 0.02]
        if c == "A":
            tmp_list[0] = .94
        elif c == "C":
            tmp_list[1] = .94
        elif c == "G":
            tmp_list[2] = .94
        elif c == "T":
            tmp_list[3] = .94
        else:
            print("you broke it")
        data[str(idx + 1)] = tmp_list

    if intended_edit is not None:
        col = str(intended_edit[0] + 1)
        data[col] = [0.02, 0.02, 0.02, 0.02]
        data[col][base_to_idx[intended_edit[1]]] = 0.35
        data[col][base_to_idx[intended_edit[2]]] = 0.65
    if tolerated_edits is not None:
        for each in tolerated_edits:
            col = str(each[0] + 1)
            data[col] = [0.02, 0.02, 0.02, 0.02]
            data[col][base_to_idx[each[1]]] = 0.20
            data[col][base_to_idx[each[2]]] = 0.80
    if het_positions is not None:
        for each in het_positions:
            col = str(each[0] + 1)
            data[col] = [0.02, 0.02, 0.02, 0.02]
            data[col][base_to_idx[each[1]]] = 0.50
            data[col][base_to_idx[each[2]]] = 0.50
    return pd.DataFrame(data, index=["A","C","G","T"])

def make_fastq_gz(path: str | Path,
                  protospacer: str,
                  n_reads: int,
                  outcomes: dict[str, tuple[float, list[tuple[int, str, str]]]]) -> None:
    
    protospacer_offset = 60
    right_pad = 151 - protospacer_offset - len(protospacer)

    amplicon = ("A" * protospacer_offset + protospacer + "A" * right_pad)
    read_counter = 0
    with gzip.open(path, "wt") as f:
        for outcome_name, (fraction, edits) in outcomes.items():
            edited = list(amplicon)
            for (pos, from_base, to_base) in edits:
                edited[protospacer_offset + pos] = to_base
            edited_amplicon = "".join(edited)
            for i in range(round(n_reads * fraction)):
                f.write(f"@read_{read_counter}\n")
                f.write(f"{edited_amplicon}\n")
                f.write("+\n")
                f.write(f"{'I' * len(edited_amplicon)}\n")
                read_counter += 1
                
