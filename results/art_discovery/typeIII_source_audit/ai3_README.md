# Deposit accompanying the manuscript "Low confidence is not absence: what alignment depth lets structure prediction say about partner proteins of phage reverse transcriptases"

Author: Chulmin Pyo (AI3 Discovery), discovery@ai3.kr. This deposit is deliberately minimal: it contains what is needed to check and re-derive the results reported in the manuscript, and nothing else.

## License
CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/) for all files in this deposit. Third-party tools and data keep their own licenses.

## Contents
- `manuscript_en.pdf`, `manuscript_ko.pdf` — the English (canonical) and Korean manuscripts. `figs/` — the figures of the manuscript (PNG).
- `numbers/NUMBERS.json`, `numbers/NUMBERS_D13.json` — every number quoted in the manuscript, computed from the result files by code (the manuscript text is filled from these two files). Each carries the SHA-256 of its source result files.
- `sequences/` — FASTA files of every analysed protein (targets, controls, tGCN5 reference), `sequences.tsv` (accession, length, source) and `d13_partners/` (the four partner sequences of the last stage with `targets.json`).
- `tables/` — per-stage result tables (per-seed calls, per-unit tables, family and member lists, census, Pfam and motif results, shuffled-alignment and long-non-GNAT analyses). Stage names follow the manuscript (D3-D13). Some JSON files contain file-system paths of the compute machine as plain strings; they are left unchanged so that the hashes stay valid.
- `models_alignments/<stage>.tar.gz` — predicted structures (PDB), the multiple sequence alignments given to the predictors (A3M), AlphaFold2 and Boltz-2 confidence outputs (pLDDT, PAE), Foldseek hit tables, and per-window score tables, with `<stage>.members.sha256` listing the SHA-256 and size of every member. All structures are predictions; no claim about enzymatic activity is made.
- `integrity/` — `SHA256SUMS` (every file of the deposit), `verify_release.py` (checks the sums and every archive member; standard library only), `preregistration_commitments.json`, `withheld_files.tsv`.

## Layout of the Zenodo record
Zenodo records are flat and hold at most 100 files, so the record holds: this README, the two manuscript PDFs, `core.tar.gz` (the directories `figs/`, `numbers/`, `sequences/`, `tables/`, `integrity/` and `models_alignments/<stage>.members.sha256`), one archive per analysis stage, and `CHECKSUMS_TOPLEVEL.sha256` (SHA-256 of every file of the record except itself). Archives larger than 90 MB are split into numbered parts, `<stage>.tar.gz.part001`, `.part002`, and so on; `cat <stage>.tar.gz.part* > <stage>.tar.gz` restores the archive, and the verifier also reads the parts directly. Download all files into one folder, then:

    tar xzf core.tar.gz
    python3 integrity/verify_release.py .

## Verify
`integrity/verify_release.py` checks `integrity/SHA256SUMS` (every file of the deposit, with the layout `models_alignments/<stage>.tar.gz`; flat copies and parts are found automatically) and every archive member.

## Not included (by design)
The analysis, scoring and pipeline code; the calibrated threshold and gate parameter files; the texts of the pre-registration documents; search profiles; vendored helper modules; internal notes and run logs; and the compute-stage work directories. They are not deposited. `integrity/preregistration_commitments.json` records, for each pre-registration document and amendment, its SHA-256, size and the git commits that touched it (commit identifiers and dates). `integrity/withheld_files.tsv` records, per category of withheld files, the number of files and one SHA-256 over the sorted per-file hashes. These records let the withheld material be checked against this deposit if it is disclosed later. The git history behind them is a local record without an external timestamp; the date of this deposit is the first external timestamp.

## Public inputs (not redistributed)
ColabFold server alignments (queries were sent to the public server), Logan protein clusters (s3://logan-pub; the object ETag and size are in `numbers/NUMBERS_D13.json` and in the `*DONE.json` files under `tables/`), PDB entries, Pfam-A 38.2. Third-party tools and data keep their own licenses.
