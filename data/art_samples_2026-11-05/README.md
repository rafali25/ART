# ART candidate sample dataset — 5 November 2026

`art_samples.csv` is the main shareable data file: **46 deduplicated RT–partner
pairs (36 type I, 10 type II)**. Each row includes source identifiers, sequences,
coding completeness, the six evidence categories, and validation limitations.
This collection includes prior references/controls and newly acquired candidates;
`cohort`, `record_role`, and `new_in_20261105_expansion` distinguish them.
The collection contains 31 complete type I pairs, 18 passing the strict ART screen.

These are computational candidates. ART-versus-NART validation is unresolved;
no experimental validation was performed in this run. Evidence indicators are
overlapping checks, not six mutually exclusive classes or calibrated probabilities.
False for profile/array evidence does not by itself establish biological absence.
Missing source metadata is blank or `unresolved` and is not an independence claim.

## Contents

- `art_samples.csv` / `.tsv`: identical sample rows; protein, CDS and upstream-window
  sequences are included directly, so the main table can be used on its own.
- `data_dictionary.tsv`: definitions and units for every column.
- `sequences/`: RT/partner protein and CDS FASTAs plus oriented locus FASTA, keyed
  by the same stable `ARTS0001`–`ARTS0046` identifiers.
- `loci/`: annotated GenBank, locus FASTA and source-coordinate maps per sample.
- `sources/`: original source sequence files; shared sources are stored once.
- `provenance/`: source manifest, three excluded aliases, original evidence,
  previous mapping verification, and input hashes.
- `summary.json`, `verification.json`, `SHA256SUMS.json`: counts and export checks.
- `code/`: the exporter used to build this snapshot.

## Evidence interpretation

1. ART-profile support uses the frozen sequence/domain E<=1e-5 and >=60% HMM
   coverage screen. It is not proof of ART function.
2. Tree affinity is descriptive. Verified NART controls are missing; all family
   validation statuses remain unresolved.
3. RVT_1 and a `[YF]xDD` pattern support RT identity. Motif presence does not prove
   catalytic activity or establish the homologous catalytic position by itself.
4. NTD length is the sequence before the first significant RVT_1 envelope. The
   >=150-aa indicator is descriptive, not a Yoon acceptance threshold.
5. Architecture requires a directly downstream co-oriented selected partner and
   >=500 nt of observed upstream annotation-free sequence. The available upstream
   window is a different quantity and can contain coding genes. Source gene
   annotations and contig boundaries limit this assessment.
6. Supported arrays retain historical cohort-specific assessments, including
   outstanding copy/boundary/family review. No joint expanded-cohort statistical
   retest was performed, and no new accepted arrays were generated in this run.

Coordinate maps use zero-based half-open intervals. CDS DNA is oriented in the
coding direction and preserves raw feature sequence; apply the recorded
`codon_start` and genetic code before translating. Complete CDSs include stop
codons; protein FASTAs omit terminal stops. GenBank coordinates retain their
standard representation. Sample IDs identify sequence-pair representatives,
not independent organisms or experimental replicates.

The three annotation/exact-pair aliases are retained as provenance, not additional
sample rows. Public source identifiers and retrieval links are supplied for
attribution; this export does not change the original data providers' licenses.
The literature and complete database archives are not bundled.
