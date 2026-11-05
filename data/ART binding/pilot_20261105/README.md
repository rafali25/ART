# ART binding development pilot — 5 November 2026

SA1 and MarsHill form a two-reference 2×2 matrix under query-only (S) and
explicitly unpaired (U) conditions. The initial stage uses one fixed seed per
cell: eight predictions. This measures input controls, execution, scoring and
runtime; it is not a confirmatory test of binding or specificity. The two loci
are close relatives; alternative pairings are not experimentally verified
nonbinders. Their close clades are reserved from future confirmation.

U uses shallow local homolog alignments. Native Boltz unpaired flags do not
guarantee that homologs occupy separate rows; a recorded input adapter separates
every non-query homolog into a chain-specific row. Processed matrices and audits
are retained, and the exact adapter must be used during subsequent calibration.
S is query-only. Templates, interaction restraints and affinity prediction are
disabled; no external MSA service or API key is used.

The current catalogue has zero eligible three-locus confirmatory blocks. Family,
independence, calibration and power requirements remain unresolved. Inspect
`manifest.json`, `cells.tsv`, `runs/`, and the eventual `results.tsv` for actual
execution status. Never interpret a one-seed pilot contrast as a binding result.
