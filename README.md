# ART protein profiles and two-route locus discovery

## RECOMB 2027 paper

Code and data for *Selection-Aware Significance Testing for de novo Repeat-Array
Discovery in Single Genomic Loci, with Application to Array-Associated Reverse
Transcriptases* (Rafiq Ali\*, Ebad Shabbir\*, Imran Razzak; MBZUAI;
\*equal contribution).

| Paper component | Code | Outputs |
|---|---|---|
| Discovery procedure *D* and selection-aware null (Sec. 3) | `scripts/art_array_detector.py`, `scripts/art_array_assessment.py`, `config/art_array_search.json` | `results/art_arrays_assessed/` |
| Benchmark, reference loci, specificity panel (Sec. 4) | `scripts/validate_art_detection.py` | `results/art_validation/` |
| Tuned PILER-CR inside the same shuffle wrapper (Sec. 4.1) | `scripts/array_finder_shuffle_wrapper.py`, `scripts/array_finder_tools.py`, `scripts/array_finders_setup.sbatch` | `results/array_calibration_20261105/` |
| Type I cohorts and the 17-locus array testing family (Sec. 5) | `scripts/art_type1_pipeline.py`, `scripts/build_art_type1_catalog.py` | `results/art_type1_public/`, `results/art_type1_mgnify/`, `results/art_type1_catalog/` |
| Five-step acquisition (archive leads, HPGC, JPGD, UHGV, MGnify) | `scripts/art_expansion_*.py` and `.sbatch` | `results/art_expansion_20261105/` (see its `README.md` and `checkpoints/`) |
| Evidence census and released dataset | `scripts/count_art_evidence.py`, `scripts/export_art_sample_dataset.py` | `results/art_evidence_census_20261105/`, `data/art_samples_2026-11-05/` |
| Figures 1-2 / appendix figures | `paper/tools/make_figures.py` / `paper_figures/paper_figures/make_figures.py` | `paper/figures/`, `paper_figures/paper_figures/` |

**Released dataset.** `data/art_samples_2026-11-05/` holds 46 RT-partner pairs
(36 type I, 10 type II) with sequences, annotated loci, original source records,
provenance and `SHA256SUMS.json`; its `README.md` and `data_dictionary.tsv`
define every column. All statuses are computational and provisional.

**Running.** Every analysis was run as a Slurm CPU batch job; each workflow has
an `.sbatch` wrapper next to its script, and Python dependencies are pinned in
`requirements-art*.txt`. Large public inputs (the HPGC and JPGD genome archives)
are not stored here; the expansion workflow downloads them and verifies their
published checksums. The Table 1 baseline harness
(`verify_array_finder_benchmark.py`) is not yet included in this repository.


Start with [HANDOFF.md](HANDOFF.md) for the project rationale, proposal revisions,
completed work, current sequence counts, limitations, and continuation commands
for a collaborator or LLM coding agent. The latest [combined type I catalog](results/art_type1_catalog/README.md)
contains 13 provisional primary pairs (11 public additions plus two LPJP1 controls),
nine complete pair representatives, and four supported candidate array tracts.
The earlier sections below describe the history and individual workflow cohorts.

## Clone this repository

Install Git LFS before cloning to retrieve the large generated comparison table:

```powershell
git lfs install
git clone https://github.com/rafali25/ART-repo.git
cd ART-repo
git lfs pull
```

The repository includes source data, analysis results, reports, and the handoff.
Python environments, downloaded tool installations, and MMseqs scratch directories
are excluded. Recreate the environment and tools using the setup commands in
[HANDOFF.md](HANDOFF.md) before rerunning analyses.

## Type I acquisition method and local FASTA rescue

The [type I method](results/art_type1_acquisition/METHOD.md) specifies a
partner-first environmental search, complementary RT retrieval, source-contig
mapping, array extraction and diversity selection for the proposal.
The [local rescue audit](results/art_type1_acquisition/AUDIT.md) rechecks seven
local genomes using fresh viral gene calls and six-frame translation. It adds
no type I candidate at the current thresholds. All six local FASTAs duplicate
previously searched DNA.

[Current type I sequences](results/art_type1_acquisition/current_typeI/) contain
the two known LPJP1 RT-partner pairs and the one existing candidate array, with
coordinates and source annotations. They are one source-genome group, not new
independent samples. The new public acquisition search is a proposed next step.

```powershell
.venv/Scripts/python.exe scripts/audit_art_type1_fasta.py
.venv/Scripts/python.exe scripts/verify_art_type1_acquisition.py
```

## Staged public search and separate family reconstruction

The [expansion report](results/art_expansion/REPORT.md) and
[combined locus export](results/art_expansion/candidate_loci.tsv) record a small
three-genome public-phage pilot. It adds one mapped Machias RT candidate to the
nine baseline RT loci. Sequence placement and observed architecture are exported
separately. The 101-tip RT-domain tree is provisional because verified NART
references are unavailable; no complete ART-locus claim is made for the new hit.

```powershell
.venv/Scripts/python.exe scripts/setup_art_reconstruction_tools.py
.venv/Scripts/python.exe scripts/art_staged_search.py --stage 1
.venv/Scripts/python.exe scripts/reconstruct_art_family.py
.venv/Scripts/python.exe scripts/art_staged_search.py --stage 2 --prepare-only
.venv/Scripts/python.exe scripts/verify_art_expansion.py
.venv/Scripts/python.exe scripts/report_art_expansion.py
```

Stage 2's fixed two-record viral/metagenomic panel is prepared, including viral
gene calls for its unannotated contig. Its profile search is gated on the current
validation failures and missing references. Stage 3 also requires measured stage
2 cost/yield, a collection manifest and a record budget. Raw timings are retained.

The [published-procedure branch](results/art_family_reconstruction/PUBLISHED_REPLICATION.md)
uses DIAMOND links, MMseqs2 clustering, longest-contig representatives, FAMSA and
FastTree. Missing original study inputs prevent reproducing its 823/230/95 counts.
The architecture-defined clade is confined to that branch. The independent
RT-domain tree uses sequence labels and never uses partner/array evidence to set
family placement. Out-of-range proteins remain available for sensitivity analysis.

## Validation before scale search

Run `.venv/Scripts/python.exe scripts/validate_art_detection.py`, then
`.venv/Scripts/python.exe scripts/verify_art_validation.py`. The benchmark and
its limitations are in [the validation report](results/art_validation/REPORT.md).
It checks published repeat cores and full-genome reverse complements, independently
planted long-spacer/mismatch/missing-copy/extended/truncated cases, unrelated RT
controls, and profiles with entire conservative lineage groups held out.

All 35 reference-core copies are recovered; exact full-boundary reproduction is
not established. Validation exposed and repaired two alignment errors while
retaining the frozen search and statistical thresholds. Residual copy-precision
errors, missing verified NART controls and missing type III partner seeds keep
scale search on hold. A strong DNA array alone cannot receive supported ART-locus
status: `scripts/art_locus_evidence.py` also requires ART RT/partner homology,
compatible context and a complete assessment. The existing seven Staphylococcus
genomes are development references, not independent lineage-level evaluation.

## Tolerant array assessment: current final-call entry point

After generating the replication scan below, run:

```powershell
.venv/Scripts/python.exe scripts/art_array_assessment.py
# Reuse checksum-compatible null samples when repeating assessment:
.venv/Scripts/python.exe scripts/art_array_assessment.py --reuse-nulls
# Verify p-values, source mappings, alignments and sampled null replay:
.venv/Scripts/python.exe scripts/verify_art_array_assessment.py
```

This writes `results/art_arrays_assessed/`. It retains the locus-specific seed
scanner and explicit dynamic-programming chains, aligns repeat-centred segments
before fitting each PWM, rescans for degenerate copies and evaluates the whole
array. Each of 1,000 local-block and 1,000 exact dinucleotide shuffles per window
reruns seed discovery and the full best-chain search. The published-style
replication outputs remain in `results/art_arrays/`.

Final significance uses the full 6-kb search, the larger p-value from its local
and dinucleotide backgrounds, and a Bonferroni correction across all input loci
at familywise alpha 0.05. The 3-kb tests are separate diagnostics; failure there
does not prevent recovery of a distant array. Inference uses fixed sample sizes,
separately from the replication scan's adaptive 200/2,000-shuffle rule.

`array_assessments.tsv` and `locus_assessments.tsv` hold final status and flags;
`background_tests.tsv` and `null_best_scores.tsv` retain the statistical evidence.
Repeat coordinates, alignments, spacer sequences, information scores and source
coordinate maps accompany them. `accepted_*.fna` files contain only accepted
candidates. `exploratory_two_copy.tsv` has a separate evidence label and no
borrowed main-call p-values. Coding overlap, stronger short-period patterns,
low complexity, intervening genes, ambiguity and boundaries are recorded, with
blocking conditions and method choices documented in the result's `REPORT.md`.

## Dedicated ART array detector

```powershell
.venv/Scripts/python.exe scripts/art_array_detector.py
.venv/Scripts/python.exe -m unittest discover -s scripts -p test_art_array_detector.py -v
```

The custom DNA scanner reads the RT-oriented loci and writes `results/art_arrays/`.
Its parameters are in `config/art_array_search.json`. Both modes always execute:
the 3-kb initial scan (top twenty 14-mers, two mismatches, 100-450 nt start spacing,
100 mononucleotide shuffles), and the 6-kb delimitation search (all recurring
10-mers, one mismatch, base spacing 60-600 nt with 30% tolerance and a skipped
copy). A failed initial scan never gates the extended search.

Delimitation scores use conserved-block information content and must exceed
locally shuffled search maxima in both 3-kb and 6-kb windows. Each background
starts with 200 shuffles within consecutive 50-nt blocks; provisionally passing
chains below twice the largest null score trigger 2,000 total shuffles per window.
A calibrated PWM retrieves degenerate copies. Full repeat boundaries are measured
on ungapped copies and MAFFT alignments, and exported spacer lengths use the
**next full repeat start minus the current full repeat end**. Start spacings are
separate fields. Published repeat/spacer ranges are flags and never hard filters.

`locus_calls.tsv` retains initial failures and non-confirmed initial detections.
`arrays.tsv`, `repeats.tsv`, `spacers.tsv` and sequence FASTAs hold delimited
candidates and original-contig coordinate mappings. `all_seed_chains.tsv` retains
scored alternatives and coding/context exclusions; `null_scores.tsv` records the
full shuffle comparisons. Individual locus GenBank files retain original
annotations and add candidate repeat regions. See `REPORT.md` for implementation
choices, including the exact recurrence rule, finite boundary flanks and the
one-skipped-copy interpretation. Labels distinguish candidates supported by both
modes from `extended_only` discoveries. This is a documented reimplementation of
the methods, not the authors' original code or a claim of exact numerical replication.

Initial candidates remain available in `initial_arrays.tsv` and its repeat/spacer
tables and FASTAs, even if strict delimitation fails. `review_flags.tsv` records
copy-count and boundary disagreements between modes. Complete candidate tracts
are in `arrays.fna` and `initial_arrays.fna`, with `array_sequence_map.tsv` for
their original coordinates. Review disagreements before downstream RNA analysis.

## RT-oriented DNA extraction

```powershell
.venv/Scripts/python.exe scripts/extract_art_loci.py
.venv/Scripts/python.exe -m unittest discover -s scripts -p test_extract_art_loci.py -v
```

This reads the current candidate and source-record tables in
`results/art_discovery/` and writes `results/oriented_art_loci/`. Every candidate
RT gets its full annotated CDS, 6 kb upstream and 6 kb downstream where available.
The interval expands when needed to retain whole candidate partners. A negative
strand RT receives upstream DNA from higher source coordinates, and the complete
window is then reverse-complemented. Every resulting RT has positive orientation.

The outputs include combined and individual FASTA/GenBank records, individual
GFF3 annotations, `locus_manifest.tsv`, `gene_annotations.tsv`,
`coordinate_map.tsv`, per-locus `.map.json` files, and
`sequence_and_annotation_issues.tsv`. Original records are copied verbatim into
`original_sources/`, with feature-indexed JSON snapshots. Extracted context
features that cross a window boundary are marked as clipped; their complete
translations remain in the original snapshot rather than on the clipped CDS.

All internal coordinates and TSV/JSON intervals are **zero-based, half-open**.
GenBank and GFF3 coordinates are explicitly **one-based, inclusive**. For a
mapping segment `[s,e)` starting at locus offset `o`, a source interval `[a,b)`
maps to `[o+a-s,o+b-s)` on a positive segment, or `[o+e-b,o+e-a)` on a negative
segment. Segments are ordered by increasing oriented-locus coordinate. Inverse
mapping and source-sequence identity are verified during extraction.

Only an explicit `topology=circular` in the authoritative source record enables
wrapping. Unknown/linear records are clipped at their ends. A circular extraction
never repeats a source base, even when the requested flanks exceed the contig.
The manifest records topology support, wrapping, reached boundaries and flank
shortfalls. Annotation-free intervals, ambiguous DNA, explicit assembly-gap
features and annotations clipped by extraction are recorded separately. Absence
of a CDS annotation can describe a noncoding array and does not imply an assembly
gap or a missing gene.

All eight source genomes in the current candidate panel already have CDS
annotations; no genes were called or replaced in this step. For an unannotated
viral contig, generate candidate CDS annotations with a viral gene caller such
as [Pyrodigal-gv](https://github.com/althonos/pyrodigal-gv), retain caller/model/code
provenance, and run discovery before extraction. The extractor requires each
candidate RT and named partner to resolve to its source CDS and fails explicitly
when that mapping is absent. It does not repair split gene calls.

`RT_review.json` checks annotated translation tables, translation agreement,
fuzzy/joined locations, source boundaries, ambiguous RT bases and explicit gap
overlap. Profile coverage below 90% is an additional review heuristic, not a
classification or truncation threshold. Original CDS boundaries remain intact.
The run's `verification.json`, code snapshots and checksums accompany the data.

## Two-route discovery workflow

Run the independent ART RT-first and partner-first searches with:

```powershell
.venv/Scripts/python.exe scripts/art_discovery.py --config config/art_discovery.json
```

The combined output is `results/art_discovery/candidate_loci.tsv`. It records
both route flags, separate partner-family IDs, ART profile evidence, genomic
coordinates, source accessions, training-sequence status, repeat inspection
hints and unavailable search families. Each partner hit has a gene-centered
plus/minus 6-kb DNA window and an inspection of its first three upstream CDSs,
including opposite-strand and overlapping genes.

The available data support an ART RT profile, a provisional single-sequence
type I partner profile, and a type II partner profile. The two type III genes
remain separate, unbuilt family slots because curated seeds are unavailable.
Their missing profiles are recorded explicitly, not replaced by generic domain
models. Separate `seed_fasta` and optional `heldout_fasta` paths can be added to
each family in `config/art_discovery.json`; document the curation evidence in
that family's `curation` field before activating it.

The follow-up source audit is in
`results/art_discovery/typeIII_source_audit/REPORT.md`. It identifies the paper's
type III labels as ART_30, ART_31 and ART_32/L10 and documents the original
8-sequence and 6-sequence profile inputs. Verified sequences or public accession
mappings were not retrieved from the checked sources. The related public LPJP1
P1/P2 sequence deposit describes type I paralogs and does not resolve these two
type III seed sets.

MarsHill and SA1 are withheld from profile training. Calibration reports their
related-sequence recovery, three reviewed GNAT controls, and annotated phage
background controls. Type I has no independent positive holdout. The small
datasets do not justify a well-validated family-wide cutoff. Candidate type I
associations require an ART-associated partner-profile match and compatible RT
context; a GNAT annotation alone is never sufficient.

The run includes the six local annotated genomes, public GenBank annotation
matching the local SA1 FASTA, and LPJP1 from ENA. Input records, downloads,
profiles, commands, calibration details, code snapshots and hashes accompany
the main table. Some targets are training examples and are flagged accordingly.
See `results/art_discovery/REPORT.md` for interpretation limits.

Tests cover independent partner-first recovery, separate type III gene families,
orientation, overlapping CDSs, contig-end clipping and the GNAT annotation rule:

```powershell
.venv/Scripts/python.exe -m unittest discover -s scripts -p test_art_discovery.py -v
```

## Original MarsHill RT pilot

This workspace implements the first protein-profile step of the revised ART
proposal. The pilot uses local annotated phage genomes and leaves MarsHill out
of its training alignment. It does not yet define or calibrate a family-wide
ART classifier.

## Run

From PowerShell in this directory:

```powershell
python -m venv .venv
.venv/Scripts/python.exe -m pip install -r requirements-art.txt
.venv/Scripts/python.exe scripts/setup_art_tools.py
.venv/Scripts/python.exe scripts/art_profile_pilot.py all
```

The tool setup downloads the official MAFFT 7.526 Windows distribution and
official HMMER 3.4 source, and builds HMMER with project-local Cygwin tools.
The Cygwin index and package checksums, source-download hashes, and build logs
are retained under `.tools/`. The existing WSL installation is not changed.
The reference check uses the public `mgtools/myRT` repository and records its
commit and input hashes. The observed reference version for the initial run is
`fa8e3624adf4c4e57ebf3350c06d626f805e42f7`.

All derived outputs are in `results/art_profile_marshill/`. Original files under
`datasetss/` are preserved. The `prepare` stage can be run without MAFFT/HMMER;
the `run` stage performs alignment, model building, searches and reference checks.

## Inputs and curation

The search target is the annotated proteome of MarsHill, accession MW248466.1,
in `datasetss/SA1 on Marshill/`. Its known ART RT is QQM14740.1. The five training
proteins come from LY01, PALS_2, S6, Madawaska and UFV_DC4. These are phages whose
ART arrays are discussed in the supplied Yoon et al. paper (Methods, PDF p. 36).
Protein accessions, genomic coordinates, sequence hashes, translation checks,
neighbouring genes and repeat-core observations are recorded in
`curated_RT_manifest.tsv`. The seed IDs are selected explicitly; a generic
"reverse transcriptase" annotation or a catalytic motif is not the selection
criterion. S6 and UFV_DC4 RTs are annotated as hypothetical proteins.

The `datasetss/SA1/` directory contains an SA1 FASTA (MW218148.1) and a Madawaska
GenBank record (MW349129.1). They are different genomes. This workflow keys
records by accession and uses the latter as Madawaska. The SA1 genome FASTA is
not used to infer proteins, and there is no SA1 result in this pilot.

The five seeds are all from the Staphylococcus/type-II branch. They are related
to MarsHill, so its exclusion is a target-sequence holdout, not an independent
lineage holdout. Seed/target identities are reported explicitly. A profile built
from these inputs may miss other ART lineages. Type I and III seeds, curated
NART controls and separate lineage-held-out validation are still needed for a
family-wide profile.

## Commands and reporting

The effective biological commands are:

```bash
mafft --auto --thread 4 art_seed_proteins.faa > art_seed_alignment.faa
hmmbuild --amino --cpu 4 --seed 42 -n ART_Staph_typeII_pilot art_family.hmm art_seed_alignment.faa
hmmsearch --cpu 4 --seed 42 -E 1e-5 --domE 1e-5 --tblout art_candidates.tbl --domtblout art_candidates.domtbl art_family.hmm MarsHill_proteins.faa
```

The Windows launcher uses MAFFT's documented `--out` argument instead of shell
redirection. `commands.jsonl` contains the actual executed arguments. Searches
are repeated at E <= 1e-3 for sensitivity. All reported hits are preserved;
350–900 aa is a recorded replication flag, not a deletion rule. Domain endpoints
and profile coverage appear in `art_candidates.tsv`; 60% coverage is a proposed
screening convention, not a calibrated biological boundary. E-values refer to
the searched database. No 207-bit threshold or other learned decision cutoff is
borrowed from the paper.

The profile is also searched against 2,339 published myRT reference sequences.
The candidate proteins are searched with myRT's existing class profiles to
record competing RT-class matches. Those historical labels are not a newly
curated ART/NART negative set. Cross-hits require review; lack of cross-hits is
not an estimate of specificity. Scores from different profiles are not directly
treated as calibrated class probabilities or compared using a fixed margin.

Only annotated CDS translations are searched in this pilot. An unannotated or
incorrectly called gene could be missed. Genome/FASTA agreement is checked when
a matching-accession FASTA is present, and translation QC is retained for every
annotated protein. Repeat-core counts are a context check, not an implementation
of the paper's array-delimitation or shuffled-null procedure. An adjacent
hypothetical protein is a candidate partner, not evidence of physical binding.

## Main outputs

- `art_seed_proteins.faa`, `art_seed_alignment.faa`, `art_family.hmm`: model inputs and profile.
- `MarsHill_proteins.faa`, `protein_annotations.tsv`: target database and accession/coordinate mapping.
- `art_candidates.faa`, `art_candidates.tsv`: retrieved proteins and source-contig links.
- `art_candidates.*`, `art_sensitivity.*`: raw HMMER outputs at both retrieval thresholds.
- `curated_RT_manifest.tsv`, `source_records.tsv`, `seed_and_holdout_identity.tsv`: curation and provenance.
- `validation_summary.json`, `myRT_reference_*`, `competing_RT_profile_matches.tsv`: reference checks.
- `commands.jsonl`, `python_requirements.txt`, `output_sha256.json`: execution record and checksums.

## Method sources

- MAFFT Windows distribution and command-line usage: https://mafft.ddbj.nig.ac.jp/alignment/software/windows_without_cygwin.html
- HMMER source and installation: https://hmmer.org/installation.html
- HMMER guide: https://eddylab.org/software/hmmer/CURRENT/Userguide.pdf
- myRT source: https://github.com/mgtools/myRT
- myRT publication: https://doi.org/10.1093/nar/gkab1207
- ART publication: https://doi.org/10.64898/2026.10.05.755809
