# Tolerant ART array discovery with cohort assessment

The published-style replication results remain unchanged in `../art_arrays/`. This assessment adds a fixed-size inferential local-block null, an exact dinucleotide null and a correction across every input locus, including those without candidates. No fixed repeat motif is used.

## Acceptance rule

Main candidates require >=3 PWM-recovered copies and Bonferroni-adjusted p <= 0.05. Candidate p is the maximum of empirical tail probabilities against the two full 6-kb-search nulls: local-block and dinucleotide. Each uses (1 + null best scores >= observed chain score)/(1 + fixed draws). Multiplying the maximum p by all input loci controls the family across loci without an independence assumption. The 6-kb null maximizes over all rediscovered seeds/chains, so candidate scores are compared with a within-window search maximum rather than a fitted-motif null. The separate per-locus `background_tests.tsv` also tests the 3-kb window using its own observed best score as a diagnostic. Candidate 6-kb scores compared with 3-kb shuffle maxima are labelled diagnostic tail fractions, not inferential p-values. A candidate entirely outside the initial 3-kb window is therefore not blocked by absence of a 3-kb signal.

Coding overlap, intervening genes, low-complexity repeats, a stronger <60-nt period, ambiguous array sequence, boundary saturation and initial/final copy-count disagreement prevent automatic acceptance and retain the candidate for review. Contig limitations elsewhere in the locus are reported without automatically rejecting an intact array. Two-copy candidates are explicitly exploratory, have no borrowed main-call p-value, and cannot be accepted.

## Alignment and full-array evaluation

Seed-centred segments are aligned with MAFFT before fitting the PWM. Conserved aligned bounds define a fixed-width PWM, calibrated conditionally against local-block shuffles for copy recovery. Copies are rescanned/chained and their common boundaries are checked again before exporting actual spacers. Gapped alignments are explicitly projected onto columns occupied by DNA in at least 80% of copies; minority-copy gaps are missing observations. The full alignment, projection mask and indel flags are retained. Such a fixed-width PWM may miss copies with internal indels, which can instead appear as missing copies or review flags. The complete-array information score is descriptive, not substituted into the pre-PWM seed-chain null. All coordinates are zero-based, half-open.

## Nulls and limitations

The local-block null preserves mononucleotide counts within successive 50-nt blocks. Exact dinucleotide shuffles preserve endpoints and every adjacent-base count. They use the Euler last-exit-tree construction, enumerating weighted arborescences over at most four DNA vertices. Tree multiplicity weights account for parallel edges. ACGT runs separated by ambiguity are shuffled separately with fixed endpoints. This is based on the Euler construction described by [Jiang et al., uShuffle](https://pmc.ncbi.nlm.nih.gov/articles/PMC2375906/); it does not use an arbitrary random Euler walk. Both nulls rerun the entire seed-discovery/best-chain search.

Inference uses fixed-size Monte Carlo samples, avoiding an empirical p-value from the replication procedure's adaptive 200/2,000 stopping rule. Existing draws are reused only as a fixed prefix with verified input checksums; additional local draws follow the same seed schedule. The paper-style maximum-shuffle decision remains a distinct descriptive replication field. These are computational candidates, not biological validation or guarantees of precise RNA boundaries.

The exploratory file reports at most the configured number of highest-scoring two-copy tracts per locus; it is a prioritization list, not a complete calibrated two-copy census.

## Results

Status counts: `{'accepted_array_candidate': 5, 'statistically_supported_requires_review': 1, 'not_supported_after_locus_correction': 3}`.

| Locus | Copies | Adjusted p | Assessment | Flags |
|---|---:|---:|---|---|
| LC680885.1__BDE75684.1 | 5 | 0.008991 | accepted_array_candidate | none |
| MN091626.1__QDJ97640.1 | 5 | 0.01798 | accepted_array_candidate | none |
| MW218148.1__QPI16926.1 | 5 | 0.01798 | accepted_array_candidate | none |
| MW248466.1__QQM14740.1 | 5 | 0.02697 | accepted_array_candidate | none |
| MW349129.1__QQO92858.1 | 5 | 0.01798 | statistically_supported_requires_review | repeat_overlaps_CDS |
| MZ422438.1__QXN67953.1 | 5 | 0.008991 | accepted_array_candidate | none |
| MZ422438.1__QXN67955.1 | 0 | 0.1888 | not_supported_after_locus_correction | repeat_overlaps_CDS;low_complexity_repeat |
| MZ779063.1__UAJ17040.1 | 5 | 0.06294 | not_supported_after_locus_correction | none |
| OR836606.1__WRM43423.1 | 5 | 0.08991 | not_supported_after_locus_correction | gene_between_array_and_RT |

## Files and reproduction

`locus_assessments.tsv`, `array_assessments.tsv`, `repeats.tsv`, `spacers.tsv`, the corresponding FASTAs, `alignments/`, `details/`, `background_tests.tsv`, `null_best_scores.tsv` and `exploratory_two_copy.tsv` contain the scores, coordinates, sequences, alignments and decisions. `accepted_arrays.fna`, `accepted_repeats.fna` and `accepted_spacers.fna` contain only automatically accepted candidates; other candidates remain in the complete assessment files.

```powershell
.venv/Scripts/python.exe scripts/art_array_assessment.py
# Resume or redo downstream annotation using compatible saved nulls:
.venv/Scripts/python.exe scripts/art_array_assessment.py --reuse-nulls
```
