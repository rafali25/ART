# Tuned PILER-CR inside the selection-aware shuffle wrapper (5 November 2026)

This experiment asks whether detector D adds anything beyond the shuffle wrapper. It
reuses the frozen D (`refinement_scan`, `config/art_array_search.json`) and the frozen
seeded nulls. No detector, threshold or shuffle code was changed. All computation ran as
Slurm CPU jobs on `cscc-cpu-p`/`cscc-cpu-qos`.

## Design

The inputs are the 31 baseline-comparison windows with copy-level truth: 24 designed
cases with planted repeat intervals and the 7 published reference loci (6-kb upstream
windows). Their stored 1,000 local-block and 1,000 exact-dinucleotide shuffles were
regenerated from the stored per-draw seeds. The regeneration was verified by rescoring
372 shuffles with D: all reproduce the stored scores exactly, and all 31 observed D scores
also match. Tuned PILER-CR 1.06 (`-minrepeat 15 -maxrepeat 45 -minspacer 25
-maxspacer 420`, all else default) ran on each observed window and on all 62,000
shuffles. The same conservative empirical p, max(p_local, p_dinucleotide) with
(1 + exceedances)/(1 + 1000), was computed for both methods.

PILER-CR reports no numeric array score. The statistic is therefore
(copies - 1) x consensus repeat length x mean copy identity of its best array (0 if none),
mirroring S = (R - 1) * I. The best-array copy count is a secondary statistic and gives
identical decisions.

PILER-CR 1.06 segfaulted deterministically on 14 of the 62,000 shuffles, all of them
published-reference shuffles; an unlimited stack does not help. The primary p counts a
crash as an exceedance (conservative); the sensitivity analysis counts it as no array.
Every count below is identical under both.

## Result

| Inputs | D, p <= 0.05 | D, Bonferroni(31) <= 0.05 | PILER-CR wrapper, p <= 0.05 | PILER-CR wrapper, Bonferroni(31) <= 0.05 |
|---|---|---|---|---|
| 24 designed | 24 | 24 | 22 | 22 |
| 7 published references | 7 | 2 | **0** | **0** |
| all 31 | 31 | 26 | 22 | 22 |

- **Designed cases:** the two methods are close. The wrapped PILER-CR misses
  `missing_copy_1` and `missing_copy_2`, where it calls no array.
- **Published references:** the wrapped PILER-CR has no power, because it calls no array
  on any reference window. D's reference p-values range from 0.000999 to 0.00999, and the
  Bonferroni(31) floor is 31/1001 = 0.031.
- **Shuffles:** tuned PILER-CR called an array on none of the 62,000 shuffles. Its wrapped
  p is therefore the floor 1/1001 whenever it calls the observed window and 1 otherwise,
  so the wrapper adds nothing beyond PILER-CR's own call. The calibrated separation on real
  arrays comes from D, not from the wrapper.

**Discrepancy with the draft's Table 2.** This reconstruction of tuned PILER-CR calls an
array on 22/31 observed windows, whereas the draft reports 24/31 cases with a call. The
original baseline harness (`verify_array_finder_benchmark.py`, ART_BENCH repository) was
not available here. The difference may come from input extent (6-kb upstream windows
here) or from the exact tuned options. Reconcile this before quoting both numbers together.

## Files

- `pilercr_shuffle_wrapper/comparison.tsv`: one row per input with D's p-values, PILER-CR's
  observed call, best copies, score, per-background p-values, crash counts, and both
  crash treatments.
- `pilercr_shuffle_wrapper/summary.json`: the counts in the table above.
- `pilercr_shuffle_wrapper/pilercr_null_scores.tsv.gz`: 62,000 rows, one per shuffle, with
  seed, PILER-CR score, copies and crash flag.
- `pilercr_shuffle_wrapper/shuffle_regeneration_checks.tsv`: the 372 stored-versus-regenerated
  D scores.
- Finder installation: `.venv-array-finders/` (gitignored), built by
  `scripts/array_finders_setup.sbatch`, with source URLs and SHA-256 hashes in
  `.venv-array-finders/provenance.json`.

## Reproduction

```bash
sbatch scripts/array_finders_setup.sbatch
sbatch scripts/array_finder_shuffle_wrapper.sbatch
```

Code: `scripts/array_finder_shuffle_wrapper.py` and `scripts/array_finder_tools.py`.
