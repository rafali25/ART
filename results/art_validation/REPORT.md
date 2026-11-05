# ART detection validation — first milestone

**Milestone partially met; scale search remains on hold.** All 35 reference-core copies in the seven available Staphylococcus genomes are recovered. Long spacers are retained. This is not validation of the complete published metagenomic census or all ART partner types.

## Published loci and orientation

The independent reference anchors use the paper’s ATATGAATACGTAT core with at most one mismatch, solely to score reference recovery. This motif never enters discovery. Full source genomes and original annotations are retained. All nine real RT loci have identical oriented DNA, RT CDS sequence and equivalent genomic repeat coordinates after reverse-complementing the complete source genome.

| Available reference | Core copies recovered | Final locus assessment | Reason / limit |
| --- | --- | --- | --- |
| S6 | 5/5 | supported_ART_locus_candidate | Family + partner + accepted array evidence; adjusted p=0.00899 |
| PALS2 | 5/5 | supported_ART_locus_candidate | Family + partner + accepted array evidence; adjusted p=0.01798 |
| SA1 | 5/5 | supported_ART_locus_candidate | Family + partner + accepted array evidence; adjusted p=0.01798 |
| MarsHill | 5/5 | supported_ART_locus_candidate | Family + partner + accepted array evidence; adjusted p=0.02697 |
| Madawaska | 5/5 | ART_locus_requires_review | repeat_overlaps_CDS; adjusted p=0.01798 |
| MZ422438.1__QXN67953.1 | No independent copy-coordinate gold set | supported_ART_locus_candidate | Family + partner + accepted array evidence; adjusted p=0.00899 |
| MZ422438.1__QXN67955.1 | No independent copy-coordinate gold set | ART_locus_requires_review | repeat_overlaps_CDS;low_complexity_repeat; adjusted p=0.18881 |
| UFV_DC4 | 5/5 | ART_locus_requires_review | not_supported_after_locus_correction; adjusted p=0.06294 |
| LY01 | 5/5 | ART_locus_requires_review | gene_between_array_and_RT; adjusted p=0.08991 |

MarsHill repeat-start spacings reproduce Figure 2G: 213, 197, 210 and 173 nt. These are start-to-start spacings. Its conservation-defined 17-nt repeats give actual spacers of 196, 180, 193 and 156 nt. The paper also uses fixed 22-nt core-centred windows for RNAfold; none of the 35 conservation-defined boundaries exactly equals that 22-nt definition. Core recovery therefore does **not** demonstrate exact full-boundary reproduction. The independent anchor table records this distinction.

UFV_DC4 and LY01 arrays are recovered but fail the existing nine-locus corrected statistical rule (p=0.06294 and 0.08991). LY01 also has an intervening annotated gene. Madawaska requires review because repeats overlap two retained CDS annotations. These annotations were not deleted to make the references pass. The second LPJP1 RT has insufficient recovered copies and p=0.18881. LPJP1’s first RT remains supported, including a 264-nt spacer; observed spacer ranges are not hard filters.

## Synthetic and external genomic tests

There are 27 synthetic cases: three independently generated motifs per test, plus three external genomic controls. Synthetic RT/partner support is supplied as known test metadata to exercise the evidence gate; it is not a simulated HMM sensitivity estimate. Synthetic fixtures are regression/development tests, including cases that exposed alignment bugs, and are not claimed as untouched evaluation data.

| Test | Cases passing expected behaviour | Notes |
| --- | --- | --- |
| long_spacer_120 | 3/3 | Planted 28-nt repeats, 120-nt spacers |
| long_spacer_180 | 3/3 | Planted 28-nt repeats, 180-nt spacers |
| long_spacer_220 | 3/3 | Planted 28-nt repeats, 220-nt spacers |
| imperfect_copies | 1/3 | Two copies carry a substitution; copy precision is scored as well as recall |
| missing_copy | 2/3 | One repeat absent; require exactly one skipped-copy link |
| extended_beyond_3kb | 3/3 | 20 copies spanning 3,980 nt; extended search must recover additional copies |
| truncated_upstream | 3/3 | Only 800 nt upstream available; incomplete assessment required |
| unrelated_RT_with_array | 3/3 | Strong long-spacer DNA array must not imply ART family membership |
| repetitive_DNA | 3/3 | CRISPR-like short-spacer array without ART RT or partner support |
| retron | 1/1 | Ec86 M24408.1 / AAA61471.1; source fragment has only 547 nt upstream |
| DGR | 1/1 | BPP-1 NC_005357.1 / NP_958675.1, reverse-strand Brt |
| CRISPR_RT | 1/1 | MMB-1 CP002583.1 / ADZ89953.1, Marme_0669 Cas6–RT–Cas1 |

Overall: **27/30** case behaviours pass. Failures: imperfect_copies_2, imperfect_copies_3, missing_copy_3.

An array can recover all planted copies and still fail copy accuracy by adding a spurious copy. `unmatched_detected_copies` and exact-boundary counts are retained. No thresholds were relaxed to turn these errors into passes. Full repeat intervals, training alignments, PWM models and shuffled-search results accompany each case.

Each of the 30 sequence windows has 1,000 local 50-nt block shuffles and 1,000 exact dinucleotide-preserving shuffles (60,000 full-search scores). Every shuffle rediscovers recurring seeds and maximizes the chain score. The selected real chain is compared with these full-search maxima; p=(1+exceedances)/(1+1000). The larger background p is Bonferroni-corrected across all 30 benchmark loci. Minimum possible adjusted p is 30/1001≈0.02997. The original nine-reference cohort retains its separate nine-locus correction and 36,000 null scores. Initial 3-kb results are diagnostic here; the full 6-kb scan always runs.

## Protein specificity and group exclusion

The comparison FASTA contains 746 protein records (744 distinct sequences): 96 labelled retron, 502 DGR and 130 CRISPR-associated myRT references, plus the three anchor RTs and available translated neighboring CDS. These are not all independent evolutionary examples. None passes the production ART RT, type I partner or type II partner screens at the frozen E≤1e-5/domain coverage≥0.6 settings (HMMER Z=domZ=2076). Raw subthreshold hits are retained.

Pairwise full-protein identity components are recorded, then merged conservatively into Staphylococcus-phage and LPJP1 lineage groups. SA1 stays with its host-genus group even though it separates at 80% identity; LPJP1 paralogs stay together. Entire groups are excluded from training in both evaluation folds. All nine RTs are retrieved. However, the profile trained on the single LPJP1 seed also retrieves retron NP_782489.1_Retrons (coverage 0.6652). Thus RT sensitivity across these two groups does not establish family specificity or broad generalization. No cutoff was retuned using these evaluation results.

The seven Staphylococcus genomes were already used during development; SA1 and MarsHill are sequence holdouts in the production profile, not independent lineage holdouts. The myRT proteins were also screened in the earlier pilot: this is a labelled comparison panel, not a previously untouched negative set. Type I has one training partner and no positive holdout. Type II has only one conservative lineage group, so no independent partner-lineage evaluation is available. Type III first and second partner profiles remain separate and unbuilt.

## Missing comparisons and remaining work

**NART is unavailable in this benchmark.** The paper’s NART examples use internal labels / contig lengths without a retrieved accession mapping. The [source audit](NART_SOURCE_AUDIT.md) did not recover verified NART proteins or contigs. AR9-like phage membership is not sufficient: downloaded KU878088.1 has no annotated RT and was excluded. `comparison_source_manifest.tsv` explicitly records the missing NART set; there is no NART specificity result. Missing NART comparisons and type III references prevent a complete validation claim.

Before scale search: resolve residual copy-precision errors, obtain independently mapped full repeat boundaries, review the three flagged Staphylococcus references, add verified NART loci and type III partner seeds, and evaluate profiles on additional independent ART lineages and competing RT families. These are limitations of the available data and present detector, not grounds to alter the frozen acceptance thresholds.

## Repairs, evidence and rerun

Validation exposed two alignment errors. A chance first seed could veto a genuine majority repeat block; the anchor now uses the median mapped seed midpoint. A chance seed aligned entirely as gaps could discard every supported copy; it is now excluded and its copy index recorded. Original seed discovery, dynamic programming, scoring, thresholds and shuffle procedures are unchanged. Pre-repair code, partial synthetic results and original nine-locus outputs are preserved under `baseline_before_anchor_fix/`. Null-cache migration verifies source hashes and that only the downstream alignment function changed. All nine production loci were reprocessed after the repairs.

Main files: `reference_recovery.tsv`, `reference_gold.tsv`, `reverse_complement_checks.tsv`, `sequence_cases.tsv`, `synthetic_truth.tsv`, `comparison_source_manifest.tsv`, `comparison_profile_hits.tsv`, `relatedness_groups.tsv`, `group_holdout_results.tsv`, `null_scores.tsv`, `summary.json`, `verification.json`, `frozen_inputs.json`, `output_sha256.json`. Source records, genomic coordinate maps, raw HMMER files, alignments and code snapshots are retained.

```powershell
.venv/Scripts/python.exe scripts/validate_art_detection.py
.venv/Scripts/python.exe scripts/verify_art_validation.py
.venv/Scripts/python.exe -m unittest discover -s scripts -p 'test_*art*.py'
```

Sources: supplied Yoon PDF Fig.2G and Methods (PDF pp.7,32–37); [ART preprint](https://doi.org/10.64898/2026.10.05.755809); [myRT reference collection](https://github.com/mgtools/myRT) at fa8e3624adf4c4e57ebf3350c06d626f805e42f7; [Ec86](https://www.ebi.ac.uk/ena/browser/view/M24408.1); [BPP-1](https://www.ncbi.nlm.nih.gov/nuccore/NC_005357.1); [MMB-1 genome](https://www.ebi.ac.uk/ena/browser/view/CP002583.1); [Cas6–RT–Cas1 primary study](https://doi.org/10.1016/j.molcel.2018.09.013).
