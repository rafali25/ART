# Dedicated ART array detector

Custom DNA repeat scanning based on Yoon et al., Methods pp. 32-33. This is a documented reimplementation; the authors' executable code and underspecified choices were not available.

Both modes run independently on every RT-oriented upstream locus. A failed 3-kb initial scan never prevents the 6-kb scan. `initial_mode_and_extended` means both modes passed this implementation; `extended_only` means only the extended procedure yielded a delimited candidate. Neither label claims exact numerical reproduction of the authors' analysis.

## Results

Scanned 9 loci; initial-mode passes: 7; delimited arrays: 3; extended-only arrays: 1. Retested 5 borderline loci with 2000 shuffles per window.

| Locus | Initial pass | Delimited arrays | Status |
|---|---|---:|---|
| LC680885.1__BDE75684.1 | True | 1 | array_candidate_detected |
| MN091626.1__QDJ97640.1 | True | 0 | initial_only_not_confirmed_by_delimitation |
| MW218148.1__QPI16926.1 | True | 0 | initial_only_not_confirmed_by_delimitation |
| MW248466.1__QQM14740.1 | True | 1 | array_candidate_detected |
| MW349129.1__QQO92858.1 | True | 0 | initial_only_not_confirmed_by_delimitation |
| MZ422438.1__QXN67953.1 | False | 1 | array_candidate_detected |
| MZ422438.1__QXN67955.1 | False | 0 | no_supported_array_detected_not_evidence_of_absence |
| MZ779063.1__UAJ17040.1 | True | 0 | initial_only_not_confirmed_by_delimitation |
| OR836606.1__WRM43423.1 | True | 0 | initial_only_not_confirmed_by_delimitation |

## Outputs

- `locus_calls.tsv`: both modes, failures, short upstream flanks and null maxima.
- `arrays.tsv`: final boundaries, copy counts, repeat lengths, start spacings, actual spacer lengths and range flags.
- `repeats.tsv`, `spacers.tsv`, and their FASTAs: sequences and original-source coordinate segments.
- `all_seed_chains.tsv`: every scored extended candidate, including coding/non-adjacent/rejected chains.
- `null_scores.tsv`: every initial, refinement and PWM null maximum, with reproducible RNG seeds.
- `alignments/`: ungapped flanks, MAFFT alignments and logs; `details/`: seed, PWM and boundary evidence.
- `loci/`: original oriented annotations plus candidate repeat_region features.

All TSV/JSON positions are zero-based, half-open in the oriented locus; source segments explicitly record the orientation back to the original contig. GenBank files use one-based inclusive locations. Spacer length = next repeat start - current full repeat end. Start spacing is reported separately.

## Implementation choices and limits

- Initial seeds: top 20 exact-frequency 14-mers after the published diversity/homopolymer filters; up to 2 mismatches, non-overlapping matches, gaps restricted to 100-450 nt. The common spacing tolerance is 30% (the initial-method paragraph does not specify a numerical tolerance). The dense veto compares the longest <100-nt run across all tested seeds with the best regular run. The exceptional one-off exact-12-mer rescue in the paper is not applied automatically.
- Extended seeds: every 10-mer with >=2 exact valid occurrences before mismatch expansion, up to 1 mismatch, no top-seed cap or initial-mode low-complexity filter. All recurring seeds are rescanned in every shuffled window. Ambiguous bases never match. Exact recurrence as the seed eligibility rule is explicit because the paper does not specify that detail.
- Chains use a shared feasible spacing interval, a base spacing of 60-600 nt with 30% tolerance, and at most one skipped copy over a whole chain. This is a conservative interpretation of "single skipped copies allowed". Longest chains are retained per seed; their score is (copies-1) times summed column KL information relative to regional mononucleotide background. The best score across all seeds defines each shuffle maximum. The paper does not give the information-content formula or tie-breaking rules.
- Conserved boundaries are the longest block containing the seed midpoint with >=80% consensus per column, allowing one low-conservation column, with conserved endpoints. Preliminary ungapped flanks extend 80 nt either side of each seed. A PWM with 0.5 pseudocount per base is calibrated against 200 new local-block shuffles; above-maximum non-overlapping matches are chained again. The final boundaries are measured in ungapped copies and MAFFT-aligned flanks; aligned bounds determine exported repeats/spacers. Boundary search-limit flags are recorded. This finite flank is a computational search parameter, not an asserted biological maximum repeat length.
- The initial statistic exceeds all 100 mononucleotide shuffle run lengths. Extended candidates must strictly exceed local-block shuffle maxima in BOTH the 3-kb and 6-kb backgrounds, initially 200 each. Provisionally passing chains below twice the initial maximum trigger extension to 2,000 total shuffles per window (the first 200 are retained). Already failed chains cannot be rescued by adding shuffles and are not retested. Independent 50-nt blocks retain their own base composition, including the final short block.
- Empirical p = (1 + null scores >= observed)/(1 + shuffle count). These are finite-null comparisons for each window, not probabilities of ART identity or dataset-wide FDR. Null search maxima account for the within-window seed search. PWM thresholds are conditional on the fitted repeat model and are not an independent biological validation.
- Chains inside one annotated CDS with all start spacings divisible by three are set aside as coding repeats. A CDS >=300 nt wholly between the final repeat and RT prevents adjacency. Shifted seed descriptions of a tract are merged when at least two repeat intervals overlap; the chain with more copies, then higher score, supplies the PWM. Full seed evidence is retained.
- Published 15-49 nt repeats and 120-220 nt spacers are descriptive flags, never hard filters. Array spans are not capped at 4.1 kb. Failure to detect is not evidence of biological absence. No folding, RNA processing, protein binding or functional activity is inferred.

## Reproduction

```powershell
.venv/Scripts/python.exe scripts/art_array_detector.py
.venv/Scripts/python.exe -m unittest discover -s scripts -p test_art_array_detector.py -v
```

## Initial candidates and disagreement review

`initial_arrays.tsv`, `initial_repeats.tsv`, `initial_spacers.tsv` and matching FASTAs retain initial-mode candidates even when extended delimitation fails. Their boundaries are separately estimated by the same ungapped/MAFFT method. These are candidate annotations, and are not promoted to confirmed delimitation calls. `review_flags.tsv` identifies unconfirmed initial detections and copy-count/boundary disagreements.

`arrays.fna` and `initial_arrays.fna` contain complete candidate array tracts, including their spacers. `array_sequence_map.tsv` maps these full tracts to source coordinates.

| Locus | Initial copies | Refined copies at same tract | Review flags |
|---|---:|---|---|
| LC680885.1__BDE75684.1 | 5 | 5 | none |
| MN091626.1__QDJ97640.1 | 5 | none | initial_detection_not_confirmed_by_delimitation |
| MW218148.1__QPI16926.1 | 5 | none | initial_detection_not_confirmed_by_delimitation |
| MW248466.1__QQM14740.1 | 5 | 3 | initial_vs_refined_copy_count_disagreement;initial_vs_refined_repeat_boundary_disagreement |
| MW349129.1__QQO92858.1 | 5 | none | initial_detection_not_confirmed_by_delimitation |
| MZ779063.1__UAJ17040.1 | 5 | none | initial_detection_not_confirmed_by_delimitation |
| OR836606.1__WRM43423.1 | 5 | none | initial_detection_not_confirmed_by_delimitation |

Repeat boundaries and copy counts remain computational estimates. A strict PWM can retain a more conserved subset of an array and expand the inferred boundary into flanks shared by that subset. Disagreement with the initial seed chain requires manual review before using spacer sequences for downstream structural or evolutionary inference.
