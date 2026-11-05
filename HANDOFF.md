# CRISP ART project handoff

Snapshot date: **21 October 2026**, Asia/Dubai. Written for a collaborator and an LLM coding agent taking over this workspace. All paths below are relative to the project root unless stated otherwise.

This handoff summarizes the files and recorded analyses available in the project. Earlier chat history and decisions not saved in the workspace cannot be reconstructed from those files. Remote job statuses are the last saved statuses, not fresh online checks.

## Start here

We are studying **array-associated reverse transcriptases (ARTs)**, particularly whether their N-terminal region has a specific predicted interaction with the neighbouring partner protein, and how array-derived RNA might affect that interaction. We first need a traceable, sufficiently diverse collection of naturally occurring loci. The implemented work is currently **sequence acquisition, locus reconstruction, repeat-array detection, validation, and provisional evolutionary analysis**.

The latest combined type I collection contains **13 candidate RT–partner pairs: 11 public candidates added to this workspace and two LPJP1 controls**. Ten pairs have complete, translation-matched coding sequences; nine complete pairs are selected across partner clusters. Four loci have statistically supported candidate repeat tracts, including three public additions. These are computational candidates, not experimentally validated ART systems.

Across the staged expansion and both public acquisition exports, there are **25 distinct mapped candidate RT locus IDs** after removing repeated LPJP1/baseline rows. That broader inventory includes type II references and RT candidates without a primary type I partner; it must not be described as 25 validated ART systems or 25 independent samples.

**The structural specificity experiment has not been completed.** No executed Boltz/AlphaFold matrix, RNA-complex analysis, or final mechanistic conclusion was found in the recorded outputs. The current sample collection does not meet the more rigorous v2 proposal's confirmatory design gate.

Recommended reading order:

1. This file for the scientific question, counts, progress, limitations, and continuation steps.
2. [Combined type I catalog](results/art_type1_catalog/README.md) and [sample table](results/art_type1_catalog/typeI_samples.tsv) for the latest type I inventory.
3. [Validation report](results/art_validation/REPORT.md) for what the detector can and cannot establish.
4. [Revised proposal](extra/ART_Research_Proposal_v2.docx) for the detailed intended experimental design.
5. [Main README](README.md) for the earlier pipelines and their commands. Its opening acquisition section predates the public catalog; use the catalog for current type I counts.

## Biological context and why we are doing this

In the locus architecture under study, an upstream noncoding repeat array is followed by an RT gene and a separate partner gene. The **NTD** is the N-terminal portion of the RT protein; it is not the array or the partner. The catalytic RT core is another region of that same protein.

The supplied Yoon et al. study motivates a model involving RT, array-derived RNA, and partner. Our question is whether predicted contacts distinguish native RT–partner pairings from comparable alternative pairings, and whether the evidence localizes those contacts to the NTD. Shared ancestry, broad family compatibility, input alignment pairing, prediction variability, and incomplete locus annotations can all create misleading patterns. Controlled comparisons are needed before making a mechanistic interpretation.

Partner families I, II, and III are distinct search targets. The current type II examples are largely related Staphylococcus phages. Type I provides an opportunity to broaden the collection beyond that group. Type III involves two separate small partner-gene families, both still lacking verified training seeds here. **NART** refers to non-array-associated RT comparative examples; verified NART accession-linked references have not been recovered in this workspace.

An ART-like genomic arrangement does not establish CRISPR-like targeting, adaptive immunity, gene editing, reverse-transcribed product formation, or binding. Likewise, a structure-confidence score is not a measured affinity. Our immediate contribution is a reproducible catalog and evidence assessment; controlled structural comparisons are the proposed next research phase.

## The proposal and its revisions

Several proposals coexist. Their existence does not prove which draft the project owner has formally adopted. For continuation, the **v2 draft is the most detailed methodological revision found**, and its stronger controls should be reviewed before implementing the prediction experiment.

| File | What it contains and how to use it |
|---|---|
| [Original ART NTD proposal](literature/ARTNT%20CLUADE_v1.docx) | Original three aims: reconstruct/census the family, test cognate specificity with within-family co-folding matrices, and examine RNA/ternary alternatives. |
| [ART Research Proposal v2](extra/ART_Research_Proposal_v2.docx) | Revised two-part biological study with matched blocks, explicit input controls, calibration, phylogenetic uncertainty, simulation gates, domain controls, and restrained conclusions. |
| [Alternative computational proposal](extra/claude%20ko%20do.docx) | Another two-part comparative-genomics and predicted-interface formulation. Retain as a draft, not evidence of executed work. |
| [Two-part biology and AI proposal](extra/ART_Two_Part_Research_Proposal_gpt_v1.docx) | Extends the biology question to an evaluation of AI agents' mechanistic inference. That benchmark is not implemented by the current acquisition scripts. |
| [ART Bench research plan](extra/ART_Bench_Research_Plan.docx) and [workflow](extra/ART_Bench_Workflow.docx) | Earlier plans for an evidence-grounded genomic agent benchmark. Proposed dataset sizes and evaluation runs are not completed results. |
| [Testing Cognate Specificity note](extra/Testing%20Cognate%20Specificity%20in%20the%20ART%20N-Terminal%20Domain.docx) | The local file contains only a short title/date note; do not assume it contains the full proposal. |

The original proposal envisaged a 12–15-locus type I all-against-all matrix, three to five seeds per cell, an NTD-only arm, NART/cross-type comparisons, optional type II replication, and RNA-containing complexes. Its approximate 200–300 GPU-hour estimate was a planning estimate, not a measured project cost. Its assumption that type I is already well powered is not supported by the current reviewed sample collection.

The v2 proposal makes the following commitments. These are **planned requirements**, not implemented results:

- **Part 1:** versioned public RT-first and partner-first discovery; source-linked locus annotations; completeness and repeat evidence; separate core, NTD, and partner comparisons; and auditable selection/exclusion records.
- **Part 2:** disjoint blocks of three loci with complete 3 × 3 pairing matrices. Proposed matching requires 50–90% identity in both core and partner, at least 80% reciprocal alignment coverage, matching partner architecture, and partner length ratios of 0.8–1.25.
- Target **8–12 evaluable evolutionary blocks, requiring 24–36 eligible loci**. Eight blocks is a minimum design gate, not proof of adequate power. Neither 12 partner clusters nor nine shortlisted pairs demonstrates this eligibility.
- Run both **U, fixed unpaired per-chain alignments**, and **S, single-sequence inputs**, with templates and restraints disabled. Verify the actual processed features; choosing an engine alone does not guarantee that pairing information is absent.
- Use three fixed seeds per cell and the median score. Compare the diagonal mean with the off-diagonal mean **within each block**. Seeds and cells do not increase the biological sample size.
- Calibrate the scoring procedure on experimentally characterized specificity comparisons from at least two homologous interaction families before confirmatory ART inference.
- Model uncertainty at the block level with precomputed RT-core/partner tree dependence; retain the proposal's prior, relatedness, exclusion, and failure sensitivity analyses. A simple partial correlation is not sufficient phylogenetic correction.
- Before confirmation, simulate the complete decision procedure. The planned resolution target is a mean contrast of 0.10 score units, with at most 5% null positive calls and at least 80% power across registered plausible scenarios.
- For NTD attribution, establish a core-profile boundary independently of interface scores; test baseline and ±20-residue NTD cuts, core-only constructs, monomers, and full-length contact localization. The current PF00078-envelope NTD exports are operational acquisition boundaries, not this completed attribution procedure.
- RNA is a secondary, calibrated arm with fixed candidate boundaries and matched controls. Compare like configurations and contact maps; do not rank protein versus RNA binding by incompatible confidence scores.
- At eight blocks, the illustrative mandatory protein inventory is **1,656 runs** before calibration, failures, or RNA work. Runtime and GPU budget must be measured on the intended hardware.

A positive controlled contrast would support computational discrimination under the specified conditions. A null or underpowered result does not establish absence of binding or prove that the source study's interface is an artefact. Neither outcome alone proves interface co-evolution. These revised limits take precedence over stronger interpretations in early drafts when describing the intended analysis.

## Counts and what each number means

The latest type I source of truth is [summary.json](results/art_type1_catalog/summary.json). Older result directories describe their own cohorts and should not be added without deduplication.

| Quantity | Count | Meaning |
|---|---:|---|
| Broader mapped RT inventory across expansion/public/MGnify exports | 25 | Unique locus IDs, including reference loci and unpaired RT candidates; not a validated ART census. |
| Primary type I candidate pairs in the combined catalog | **13** | 11 public additions plus two LPJP1 paralog controls. |
| Public additions from the archive/Logan pilot | 4 | In `art_type1_public`; new to this workspace. |
| Public primary pairs from the MGnify follow-up | 7 | In `art_type1_mgnify`, in addition to its repeated controls. |
| Exact RT–partner sequence-pair groups in the catalog | 13 | All 13 exported RT sequences and all 13 partners are individually distinct. |
| Partner clusters at 70% identity and 80% mutual coverage | 12 | Sequence groups, not proven independent ancestry or samples. |
| Complete translation-matched primary pairs | 10 | Three other pairs require gene-boundary/completeness review. |
| Complete pairs selected across eligible partner clusters | **9** | Starting shortlist, not a confirmatory structural panel. |
| Primary pairs with less than 6 kb of upstream DNA | 6 | Upstream assessment is incomplete; this is separate from protein completeness. |
| Unique RT candidates in the combined public array-testing family | 17 | Includes four RT candidates without a primary type I pair. Duplicate LPJP1 controls counted once. |
| Supported candidate arrays in the combined type I catalog | **4** | Three public additions plus one LPJP1 locus. |
| Repeat and spacer records exported for those four arrays | 18 repeats; 14 spacers | Components of four tracts, not extra independent loci. |
| Earlier baseline discovery | 9 RT loci in 8 genomes | Seven Staphylococcus/type II reference loci plus two LPJP1/type I loci. |
| Staged expansion | 10 RT loci | The same nine baseline loci plus Machias; not ten further additions. |
| Baseline tolerant assessment | 5 accepted arrays; 1 supported requiring review | A different nine-locus cohort; do not add its LPJP1 result to the catalog again. |
| Verified NART references; built type III partner profiles | 0; 0 | Missing inputs, not evidence those systems do not exist. |

The 25-locus inventory is the set union of `locus_id` in `results/art_expansion/candidate_loci.tsv`, `results/art_type1_public/candidate_loci.tsv`, and `results/art_type1_mgnify/candidate_loci.tsv`. It contains eight type II-associated candidates including Machias, 13 primary type I pairs, and four further RT candidates. It is a count of mapped records, not a deduplication into independent evolutionary samples.

The source paper's 823 proteins / 230 clusters / 95 ART representatives are **published comparison outcomes recorded in the reconstruction report**, not our recovered sequence counts. The original proposal also quotes a much larger partner census. Such literature numbers require checking against the cited source version before reuse; none describes the current local collection.

## Current type I sequence inventory

The RT ID is the final component of each `locus_id`; use the exact IDs in [typeI_samples.tsv](results/art_type1_catalog/typeI_samples.tsv). `Complete` below means both CDS translations match; it does not mean the entire genomic locus is assessable.

| Source contig | RT ID | Primary partner ID | Protein pair | Combined array evidence |
|---|---|---|---|---|
| SRR8925777_19723 | SRR8925777_19723_pgv9 | SRR8925777_19723_pgv10 | Complete | Not supported |
| ERR10896470_6836 | ERR10896470_6836_pgv1 | ERR10896470_6836_pgv2 | Complete | Incomplete upstream, 232 nt |
| ERR8055547_17635 | ERR8055547_17635_pgv12 | ERR8055547_17635_pgv13 | Complete | Supported candidate, 5 copies |
| MZ422438.1, LPJP1 | QXN67953.1 | QXN67952.1 | Complete | Supported candidate, 5 copies |
| SRR26805139_340 | SRR26805139_340_pgv1 | SRR26805139_340_pgv2 | Complete | Incomplete upstream, 233 nt |
| MZ422438.1, LPJP1 | QXN67955.1 | QXN67954.1 | Complete | Not supported |
| MGYC001394629236 | MGYC001394629236_pgv6 | MGYC001394629236_pgv5 | Review required | Incomplete upstream, 1 nt |
| MGYC003138808455 | MGYC003138808455_pgv4 | MGYC003138808455_pgv3 | Review required | Incomplete upstream, 2 nt |
| MGYC003138808455 | MGYC003138808455_pgv2 | MGYC003138808455_pgv1 | Review required | Incomplete upstream, 3,724 nt |
| MGYC002687634053 | MGYC002687634053_pgv52 | MGYC002687634053_pgv53 | Complete | Supported candidate, 5 copies; long-gap review |
| MGYC001394626043 | MGYC001394626043_pgv9 | MGYC001394626043_pgv10 | Complete | Not supported |
| MGYC003526048914 | MGYC003526048914_pgv179 | MGYC003526048914_pgv180 | Complete | Supported candidate, 3 copies |
| MGYC002345864168 | MGYC002345864168_pgv4 | MGYC002345864168_pgv3 | Complete | Incomplete upstream, 5,791 nt |

ERR10896470 and SRR26805139 share a partner cluster; the complete-pair shortlist chooses SRR26805139. The two LPJP1 pairs come from one genome and remain one training-source group, even though they are distinct paralogs. Two primary RT loci also share MGYC003138808455. Neither contig IDs nor cluster IDs establish biological independence.

| Supported candidate tract | Repeat copies | Inferred spacer gaps in nt | Combined corrected p |
|---|---:|---|---:|
| ERR8055547_17635 / pgv12 | 5 | 284, 326, 630, 289 | 0.0084957521 |
| LPJP1 / QXN67953.1 | 5 | 148, 264, 137, 199 | 0.0084957521 |
| MGYC002687634053 / pgv52 | 5 | 219, 226, 210, 689 | 0.0084957521 |
| MGYC003526048914 / pgv179 | 3 | 122, 121 | 0.0169915042 |

The combined p-values correct across all 17 public-cohort RT candidates. MGYC002687634053 has a flagged gap that may span multiple missing repeats; inferred gaps are not necessarily single biological spacers. Statistical significance supports a repetitive tract, not its full boundaries, transcription product, or function.

## Work completed so far

### 1 Original MarsHill profile pilot

`scripts/art_profile_pilot.py` built a five-seed Staphylococcus/type II RT profile from LY01, PALS_2, S6, Madawaska, and UFV_DC4, withholding MarsHill. It searched the annotated MarsHill proteome and compared with the available myRT references. MarsHill RT is QQM14740.1, accession MW248466.1. This demonstrated local retrieval and established provenance/QC outputs, not family-wide classifier validation. See [pilot report](results/art_profile_marshill/REPORT.md).

### 2 Two-route discovery and source audits

`scripts/art_discovery.py` implemented independent RT-first and partner-first searches, per-family profiles, context inspection, and source-linked candidate exports. All nine baseline candidates were found by both routes among 2,076 annotated proteins in eight genomes. MarsHill and SA1 were withheld from production-profile training as related-sequence checks.

The production ART RT model uses six seeds; the original type I model uses one partner seed; the type II model uses five. Type III first and second partners remain separate unbuilt slots. The [type III source audit](results/art_discovery/typeIII_source_audit/REPORT.md) traced internal paper labels and the related archive but did not recover verified seed sets. See [discovery report](results/art_discovery/REPORT.md).

### 3 RT-oriented locus extraction

`scripts/extract_art_loci.py` extracted the nine baseline RT CDSs with up to 6 kb upstream/downstream, expanding to retain full candidate partners. Outputs include oriented FASTA/GenBank/GFF3, original source records, annotations, manifests, coordinate maps, translation checks, and checksums. Negative-strand RT loci are reverse-complemented into positive RT orientation. See [extraction report](results/oriented_art_loci/REPORT.md).

### 4 Dedicated DNA repeat-array detection

`scripts/art_array_detector.py` implemented a locus-specific seed scanner, dynamic-programming repeat chains, information-content scoring, shuffled nulls, PWM extension, and repeat/spacer exports. Both the initial 3-kb and extended 6-kb modes execute independently. The original replication-style run yielded seven initial candidates and three delimited arrays across nine loci. Those historical counts are not the final tolerant-assessment count. See [detector report](results/art_arrays/REPORT.md).

`scripts/art_array_assessment.py` added aligned repeat-centred PWM fitting, tolerant copy recovery, context flags, exact dinucleotide and local-block backgrounds, and locus-wide multiple-testing correction. The latest baseline assessment has five accepted array candidates, one statistically supported candidate requiring review, and three not supported after correction. Two-copy exploratory tracts remain separately labelled. See [assessment report](results/art_arrays_assessed/REPORT.md).

### 5 Validation and detector repairs

`scripts/validate_art_detection.py` checked source-genome reverse complements, the seven Staphylococcus reference arrays, planted synthetic cases, competing RT controls, and conservative lineage-group exclusions. It recovered **35/35 reference repeat cores** and preserved equivalent sequence/coordinate mappings for all nine real loci after reversing the full source genomes.

The biological behaviour benchmark passes **27/30 cases**, with unresolved failures in `imperfect_copies_2`, `imperfect_copies_3`, and `missing_copy_3`. Detecting a planted copy plus an extra false copy still fails precision. Full reference repeat boundaries have not been independently reproduced. A verification JSON can report `passed: true` for internal consistency while recording these scientific validation failures.

Two alignment bugs were repaired: a chance first seed vetoing a genuine repeat block, and an all-gap seed discarding supported copies. Earlier outputs/code are retained under `results/art_validation/baseline_before_anchor_fix/`. Thresholds were not relaxed to manufacture passes. Full-search nulls and source mappings are retained. The family-evidence gate in `scripts/art_locus_evidence.py` prevents a strong repeat array alone from becoming a supported ART-locus call.

The production screens did not accept the labelled competing-RT comparison panel at frozen thresholds. However, the LPJP1-only group-held-out RT profile cross-hit a retron, and no verified NART controls were available. These checks do not establish universal specificity. See [validation report](results/art_validation/REPORT.md) and [NART audit](results/art_validation/NART_SOURCE_AUDIT.md).

### 6 Staged expansion and family reconstruction

`scripts/art_staged_search.py` searched three additional public annotated phages: Machias, AR9, and phiR1-37. Machias contributed RT QQO92593.1, MW349128.1, with downstream type II partner QQO92594.1. No calibrated supported array was established for this addition. AR9 and phiR1-37 had no accepted candidates at the available thresholds; AR9 was not promoted to a NART control.

`scripts/reconstruct_art_family.py` ran two separate analyses:

- A published-procedure branch using DIAMOND connectivity, MMseqs2 90% identity / 80% mutual coverage, longest-contig representatives, FAMSA, and FastTree. The available MarsHill component has ten RT proteins and seven representatives.
- An independent sequence-only RT-domain analysis using PF00078 envelopes ±40 aa, MAFFT L-INS-i, 50% occupancy trimming, and FastTree LG+Gamma. The staged tree has 101 tips and 272 columns. It is provisional, not the full published phylogeny.

Architecture-defined clade selection is confined to the replication branch. Arrays and partners do not assign sequence-family placement. Missing original reference membership, database representatives, and cluster-member mappings prevent exact reproduction of the paper's family counts. Available myRT input has 2,339 entries, of which 2,328 canonical proteins enter this branch; it is not silently substituted for the paper's reported 2,019-reference set.

Stage 2 prepared PX848769.1 and PV020493.1, including 191 Pyrodigal-gv CDS predictions for the latter unannotated partial record. Stage 2 profile search and Stage 3 search were not executed because their gates remain unmet. Later bounded public pilots were separate provisional acquisitions; they do not mean these gates passed. See [expansion report](results/art_expansion/REPORT.md) and [replication requirements](results/art_family_reconstruction/PUBLISHED_REPLICATION.md).

### 7 Local type I rescue audit

`scripts/audit_art_type1_fasta.py` checked seven unique local genomes, with LPJP1 as an additional control, using source annotations, Pyrodigal-gv, and six-frame translation. It searched 15,471 translations and added **zero** local type I candidates at the current thresholds. All six local FASTA files duplicate previously searched DNA. This audit did not demonstrate biological absence. See [method](results/art_type1_acquisition/METHOD.md) and [audit](results/art_type1_acquisition/AUDIT.md).

### 8 Bounded public acquisition

The archive/Logan pilot recovered four new provisional primary type I pairs, plus the two LPJP1 controls. Archived protein leads were treated as discovery inputs, then exactly mapped back to historical and current DNA. Changed contig numbers were not accepted without sequence evidence. One new supported array was found. Historical/current bounded compressed retrieval totalled 59,602,523 bytes. See [public pilot report](results/art_type1_public/REPORT.md).

Four EBI HMMER/MGnify jobs completed. Their collected results contain 418 query-hit rows representing 293 unique MGYP proteins; each RT query was capped at 200 results and each partner query returned nine hits. These are remote hit counts, not mapped locus counts. Query-balanced acquisition retained ten retrievable MGnify contigs; 13 queried proteins mapped exactly. The resulting cohort contains 13 RT candidates, nine primary type I pairs including controls, seven public primary additions, and two new supported array candidates. See [MGnify report](results/art_type1_mgnify/REPORT.md).

Both pilots retain original annotations, predicted genes, six-frame rescue, exact translation mappings, source snapshots, database versions, raw HTTP responses, ETags/checksums, oriented loci, CDS DNA, operational NTDs, arrays, repeats/spacers, and sequence-only placement. The additional two-seed type I retrieval model uses the two LPJP1 paralogs, still one training-genome group. No expanded model was trained on the newly recovered candidates.

### 9 Combined catalog

`scripts/build_art_type1_catalog.py` merges both public cohorts, removes duplicated LPJP1 controls, clusters primary partners, selects complete representatives, and applies a combined array correction over all 17 unique candidate RT loci. This is the current handoff entry point for type I sequence selection. Cohort reports retain their own six- and thirteen-locus corrections; use combined p-values for statements about the merged catalog.

## Main files and data to give the collaborator

| Path | Purpose |
|---|---|
| `datasetss/` | Original supplied phage genomes. Preserve the spelling and original files. |
| `literature/` and `extra/` | Source study, original proposal, revised proposals, and benchmark plans. |
| `scripts/` | Active implementation and tests. Result-local `code/` folders are archived snapshots. |
| `config/` | Search thresholds, array rules, validation settings, bounded source manifests. |
| `results/art_type1_catalog/typeI_samples.tsv` | Full current 13-pair type I inventory and evidence/status columns. |
| `results/art_type1_catalog/all_RT_candidate_loci.tsv` | All 17 unique RT candidates in the combined public testing family. |
| `results/art_type1_catalog/complete_pair_representatives.tsv` | Nine complete pair representatives across partner clusters. |
| `results/art_type1_catalog/RT_proteins.faa` | Thirteen primary type I RT protein sequences. |
| `results/art_type1_catalog/primary_partners.faa` | Thirteen matched primary partner proteins, keyed by locus ID. |
| `results/art_type1_catalog/NTDs.faa` | Thirteen operational NTD extractions; boundaries require further review. |
| `results/art_type1_catalog/supported_array_candidates.tsv` and `.fna` | Four supported candidate DNA tracts with source mapping links. |
| `results/art_type1_catalog/array_repeats.fna` and `array_spacers.fna` | Eighteen inferred repeat copies and fourteen gaps from those tracts. |
| `results/art_type1_public/` and `results/art_type1_mgnify/` | Cohort evidence, original source records, full locus files, arrays, nulls, and provenance. |
| `results/art_validation/` | Benchmark truth, failures, controls, group-held-out results, code/input snapshots. |
| `results/art_expansion/` and `results/art_family_reconstruction/` | Broader mapped RT panel, staged gates, alignments and trees. |
| `.tools/` | Local external executables, myRT/Pfam references, download/build provenance. Needed for reruns. |
| `.venv/` | Existing Windows Python environment. Recreate it on another machine. |
| `tmp/reading/` | Convenience text extractions. Original DOCX/PDF and result tables remain the authorities. |

Copy the **whole project folder** when possible. The handoff itself is a guide; it does not embed the sequence data or replace the source folders. If transfer space is limited, preserve `datasetss/`, `literature/`, `extra/`, `config/`, `scripts/`, requirements, README/handoff, and all result evidence/provenance. `.venv/` is machine-specific; `.tools/` can be rebuilt, but retain tool provenance/reference snapshots. Cached nulls are useful for reproducible, faster reruns. The workspace is now being versioned for `https://github.com/rafali25/ART-repo.git`; the large `results/art_family_reconstruction/replication/all_against_all.tsv` uses Git LFS. Install Git LFS and run `git lfs pull` after cloning. Rebuildable MMseqs scratch directories, local environments, and tool installations are excluded from Git.

### Important naming and coordinate details

- The local `datasetss/SA1/` directory contains an SA1 FASTA **MW218148.1** and a Madawaska GenBank **MW349129.1**. They are different genomes. Key by accession and sequence, not folder name. MarsHill is MW248466.1; Machias is MW349128.1.
- The baseline `art_discovery/candidate_loci.tsv` explicitly labels its coordinate columns as **one-based**. Later internal mapping/assessment TSV and JSON intervals are generally **zero-based, half-open**. GenBank/GFF3 use **one-based, inclusive**. Read the actual field names/schema before converting; do not apply one blanket convention to every historical table.
- An RT on the negative source strand is oriented positively in its extracted locus. Use its `.map.json` segments for original-coordinate mapping, including explicitly circular records that wrap.
- Spacer length is **next full repeat start minus current full repeat end**. Repeat start-to-start spacing is a different quantity. For MarsHill, recorded start spacings are 213, 197, 210, 173 nt; inferred spacer gaps are 196, 180, 193, 156 nt.
- `pgv` IDs refer to predicted gene calls. Operational NTDs end before the first significant PF00078 envelope and are not experimentally established domains.
- Logan SRR/ERR contig identifiers and MGYC/MGYP identifiers are release/source snapshots, not invented INSDC accession versions. Keep release IDs, ETags, exact mappings, and hashes.

## Statistical rules already implemented

The frozen production profile screen uses sequence/domain E ≤ 1e-5 and at least 0.6 profile coverage; relaxed reporting is retained separately. A generic GNAT annotation is not a type I partner assignment. Partner-first searches inspect up to the first three upstream CDSs, including opposite-strand/overlapping genes. Primary pairing respects orientation/context and intervening RTs.

The initial and extended repeat scans both run. The initial 3-kb result is not a gate that can prevent a distant 6-kb candidate from being examined. Repeat seeds are discovered within each locus; the known Staphylococcus core is used for reference scoring, not hardcoded into discovery. Published repeat/spacer ranges are review flags, not universal rejection filters.

Final array support uses full-search maxima from both local-block and dinucleotide-preserving nulls. Each shuffle repeats seed discovery and best-chain selection. The larger empirical p-value is used, then corrected across the declared locus cohort. Baseline tolerant assessment uses 1,000 draws per background per window; public-pilot positive calls require 2,000 draws per background after a 200-draw screen. Keep these procedures distinct from the original replication-style adaptive 200/2,000 rule.

The combined catalog uses Bonferroni correction across 17 unique tested RT loci, not just the 13 successful primary pairs. This is a conditional selected-panel comparison; it does not correct the external lead-selection process or establish a database-wide error rate. Changed DNA, search settings, detector code, or testing-family membership require reassessment of compatible nulls and correction. Do not transfer a p-value to a different tract or exploratory two-copy call.

## Current blockers and unanswered questions

1. **Data sufficiency:** nine complete representatives do not meet the v2 target of at least 24 eligible loci/eight blocks. Matching and phylogenetic independence have not been demonstrated by the 70% partner clustering.
2. **Repeat-copy precision:** three benchmark cases still fail. Repeat cores are recovered, but full boundaries and missed/extra-copy behaviour need independent review.
3. **Incomplete loci:** three type I pairs need gene review; six primary loci lack a full 6-kb upstream interval. Protein completeness and upstream completeness are different decisions.
4. **Family calibration:** the expanded type I model still has one training-source genome and no verified independent positive holdout group. New candidates are not automatically gold-standard positives.
5. **Missing comparative references:** no verified NART panel or verified seeds for either type III partner family. Do not manufacture these labels from architecture alone.
6. **Remote acquisition:** NCBI BLAST batch `DEP6X50X016` was `WAITING` at the last saved poll. Unretrieved MGnify hits remain deferred, not negative. Current remote state is unknown until polled.
7. **Reconstruction scope:** original study database inputs are missing; the small executed branch does not reproduce its exact numerical family size or ART/NART boundary.
8. **Prediction phase:** calibration, engine/input-feature audit, block selection, core-based NTD definitions, simulations, GPU cost measurements, protocol freeze, and structural runs remain to be done.
9. **RNA interpretation:** exported arrays are DNA. Transcript ends, mature RNA processing, and the appropriate molecular composition for a calibration complex have not been established here.

## How to continue

First verify the current catalog and inspect the nine complete representatives, the four supported candidate tracts, the three gene-review pairs, and all selection flags. Review repeat boundaries without tuning to make the known references pass. Fix copy-precision defects against truth cases while retaining the frozen baseline and documenting any new method version.

Next obtain verified independent partner/RT references and NART/type III accession mappings where possible. Continue bounded public acquisition from recorded jobs and manifests. Exact protein-to-DNA mapping and retained provenance are prerequisites for calling a mapped locus. Track all retrieved, deferred, excluded, and ambiguous records.

Then build core/partner alignments and reviewed trees and test whether the v2 matching/block gate is feasible. Keep development groups separate from confirmatory groups. If the required diversity is unavailable, report a catalog/feasibility study or exploratory prediction pilot; do not lower eligibility thresholds after seeing scores.

Only after selecting an adequate panel should the project commit to the mandatory calibration, U/S input controls, simulated analysis checks, measured GPU budget, registered NTD boundary procedure, and frozen prediction manifest. Type II replication and the AI-agent benchmark are separate future arms, not completed or automatically authorized by the existence of their drafts.

## Commands for a coding agent

Current workspace: `D:\MBZUAI\CRISP`, Windows PowerShell. Run from the project root. Active scripts resolve their root from the source tree, but existing environments/binaries and some archived logs are Windows-specific.

For a fresh environment, create a project-local venv and install the pinned requirements. The reconstruction requirements include the base requirements:

```powershell
python -m venv .venv
.venv/Scripts/python.exe -m pip install -r requirements-art-reconstruction.txt
.venv/Scripts/python.exe scripts/setup_art_tools.py
.venv/Scripts/python.exe scripts/setup_art_reconstruction_tools.py
```

These setup scripts download/build tools and references; they are not needed merely to read the handoff or existing sequences. Recorded tools include MAFFT 7.526, HMMER 3.4 with local Cygwin support, DIAMOND 2.1.15, MMseqs2 17-b804f, FAMSA 2.5.2, FastTree, Pyrodigal-gv 0.3.2/Pyrodigal 3.7.1, myRT, and Pfam RVT_1. Consult `.tools/` provenance for exact builds/hashes. On another OS, adapt and verify executable setup before assuming the Windows wrappers will work.

To check existing scientific exports and software, without re-running acquisition:

```powershell
.venv/Scripts/python.exe scripts/verify_art_type1_public.py
.venv/Scripts/python.exe scripts/verify_art_type1_public.py --config config/art_type1_mgnify.json
.venv/Scripts/python.exe scripts/verify_art_validation.py
.venv/Scripts/python.exe scripts/verify_art_expansion.py
.venv/Scripts/python.exe -m unittest discover -s scripts -p 'test_*art*.py'
```

Verifier scripts may refresh their verification artifacts. Passing software/internal-consistency checks does not erase the 27/30 biological benchmark result. These commands are provided for continuation; they were not all rerun during handoff authoring.

To rebuild a cohort, run its stages **sequentially** because they update shared files. Below is the MGnify cohort; omit `--config` to use the original public pilot:

```powershell
.venv/Scripts/python.exe scripts/art_type1_pipeline.py process --config config/art_type1_mgnify.json
.venv/Scripts/python.exe scripts/art_type1_pipeline.py nulls --config config/art_type1_mgnify.json
.venv/Scripts/python.exe scripts/art_type1_tree.py --config config/art_type1_mgnify.json
.venv/Scripts/python.exe scripts/art_type1_pipeline.py diversity --config config/art_type1_mgnify.json
.venv/Scripts/python.exe scripts/finalize_art_type1_public.py --config config/art_type1_mgnify.json
.venv/Scripts/python.exe scripts/verify_art_type1_public.py --config config/art_type1_mgnify.json
```

After finalizing both cohorts, rebuild the combined catalog:

```powershell
.venv/Scripts/python.exe scripts/build_art_type1_catalog.py
```

To resume existing remote work rather than create duplicate submissions:

```powershell
.venv/Scripts/python.exe scripts/art_type1_public_search.py poll
.venv/Scripts/python.exe scripts/art_type1_public_search.py collect
.venv/Scripts/python.exe scripts/art_type1_public_search.py acquire --config config/art_type1_mgnify.json
```

Inspect the recorded source manifest and acquisition limits before fetching another panel. The MGnify config caps acquisition at ten contigs, 24 acquisition proteins, 200 collected hits per query, and 20 MB of new source downloads per invocation. The archived Logan prefix budget is separate. Changing a cap does not automatically authorize calling an expanded panel a representative census.

Earlier baseline reruns are documented in README and the linked reports: `art_profile_pilot.py all`, `art_discovery.py`, `extract_art_loci.py`, `art_array_detector.py`, `art_array_assessment.py`, and `validate_art_detection.py`. Run those in dependency order only when needed. Stage 2/3 expansion gates are still active; do not bypass them to turn incomplete validation into a definitive locus claim.

## LLM coding agent operating notes

Read this handoff and the relevant result report before editing. Treat `scripts/` and `config/` as active code/settings; preserve original inputs and frozen result-local code snapshots. Determine which cohort a proposed change affects before rerunning anything. Do not silently overwrite older evidence to make a new method appear to be the original run.

Keep separate fields for profile evidence, sequence placement, partner context, array evidence, protein completeness, upstream assessability, and biological interpretation. Record missing evidence explicitly. Do not use arrays to train/label a sequence-family boundary and then claim independent array validation of that boundary. Do not replace missing partner families with generic domains.

When sample membership or analysis code changes, refresh affected hashes, null-cache compatibility, testing-family corrections, exports, verifications, reports, and this handoff together. Retain failed jobs and exclusions. New loci are new to the workspace unless novelty has been established by an independent literature/source review.

The handoff creation checked all **18 catalog output checksum entries** and all **six input-manifest hashes**, with no mismatches. It also checked the 13 RT/partner/NTD FASTA records, the four supported array records, the 18 repeat/14 spacer records, and the 25-ID cross-export union. It did not rerun the entire scientific pipeline, update remote jobs, or validate molecular function.

Suggested prompt to begin a future agent session:

```text
Read HANDOFF.md, results/art_type1_catalog/README.md,
results/art_type1_catalog/summary.json, and results/art_validation/REPORT.md.
Inspect the relevant scripts/configs and current verification artifacts.
Summarize the current candidate counts and remaining validation/design gates.
Continue the specific task I assign while preserving source provenance,
coordinate conventions, cohort corrections, and provisional evidence labels.
Use ART_Research_Proposal_v2.docx for the revised design unless I specify
another proposal. Do not assume the structural experiment has already run.
```

## Source references

Scientific background and literature quantities above are attributed to the supplied study and local proposal/report files; this handoff did not conduct a new online literature audit. The primary local study is [Yoon et al. 2026](literature/Yoon%20et%20al.%202026%20Array%20Associated%20Reverse%20Transcriptases%20%282%29.pdf). Its DOI recorded by the project is `10.64898/2026.10.05.755809`. The archive used as discovery input is recorded as `10.5281/zenodo.23086637`; the myRT snapshot commit is `fa8e3624adf4c4e57ebf3350c06d626f805e42f7`. Use the source/download manifests and original documents for accession-level and source-version verification before publication.
