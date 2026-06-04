# Changelog

All notable changes to the AQ CRISPResso Pipeline are documented here.

## [Unreleased]
### Added
- Nuclease (Cas9 KO) editor support — Phase 1, single-guide (in progress)
  - Optional `min_alignment_score` column in amplicon_list.csv (default 60)
  - `editor=NUCLEASE` rows load without `intended_edit`/`tolerated_edits`
  - CRISPResso runs nuclease samples with a cut-site-centered quantification
    window (`-wc -3`, `-w 15`) and the configured min alignment score

