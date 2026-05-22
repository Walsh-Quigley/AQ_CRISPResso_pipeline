# AQ Pipeline V2 — Architecture

## Call Graph

```mermaid
graph TD
    CRISPResso_Loop.py --> pipeline/crispresso.py
    Quantification_Loop.py --> pipeline/quantify.py

    pipeline/crispresso.py --> load_amplicon_list
    pipeline/crispresso.py --> run_crispresso

    pipeline/quantify.py --> read_mapping_stats
    pipeline/quantify.py --> read_allele_table
    pipeline/quantify.py --> read_quant_window
    pipeline/quantify.py --> generate_search_sequences
    pipeline/quantify.py --> generate_oneseq_search_sequences
    pipeline/quantify.py --> analysis/abe.py
    pipeline/quantify.py --> analysis/oneseq.py
    pipeline/quantify.py --> analysis/heterozygous.py
    pipeline/quantify.py --> export

    load_amplicon_list --> AmpliconConfig
    generate_search_sequences --> reverse_complement

    subgraph loaders/
        load_amplicon_list
        read_mapping_stats
        read_allele_table
        read_quant_window
        export
    end

    subgraph utils/
        reverse_complement
        generate_search_sequences
        generate_oneseq_search_sequences
    end

    subgraph config.py
        AmpliconConfig
    end
```

## Module Summary

| File | Status | Purpose |
|---|---|---|
| `config.py` | Done | AmpliconConfig dataclass |
| `utils/sequences.py` | Done | reverse_complement, generate_search_sequences, generate_oneseq_search_sequences |
| `loaders/amplicon_list.py` | Done | Parse amplicon_list.csv into AmpliconConfig objects |
| `loaders/crispresso_output.py` | Done | Read CRISPResso allele table and mapping stats |
| `loaders/export.py` | Done | Write summary CSV and PRISM output |
| `analysis/abe.py` | Done | ABE correction metric calculations |
| `analysis/heterozygous.py` | Done | Het position detection and allele splitting |
| `analysis/oneseq.py` | Done | ONE-seq A-to-G analysis |
| `pipeline/crispresso.py` | Done | Stage 1 — amplicon matching, run CRISPResso on each sample |
| `pipeline/quantify.py` | Done | Stage 2 — orchestrate analysis, assemble results |
| `CRISPResso_Loop.py` | Done | Entry point for Stage 1 |
| `Quantification_Loop.py` | Done | Entry point for Stage 2 |
