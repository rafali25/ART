# Part 2 — Controlled tests of predicted ART RT–partner specificity

**Research question:** Does an array-associated reverse transcriptase (ART) receive greater predicted interface support with its naturally co-occurring partner than with comparable partners from other loci, and is that discrimination associated with its N-terminal domain (NTD)?

**Author:** Rafiq Ali, Mohamed bin Zayed University of Artificial Intelligence.

**Version:** 5 November 2026.

**Status:** Development experiments started on 5 November 2026. Eight actual local input-feature checks, two local GPU smoke predictions, and four hosted query-only predictions are complete. The hosted endpoint's fixed-seed, processed-feature and PAE limitations are recorded separately; its results are exploratory. Confirmatory binding analyses have not been executed. See the [experiment record](<../data/ART binding/README.md>).

This document develops the protein-binding questions in Aims 2 and 3 of [ARTNT CLUADE_v1.docx](<ARTNT CLUADE_v1.docx>). It incorporates the stronger experimental-design and interpretation safeguards in [ART Research Proposal v2](../extra/ART_Research_Proposal_v2.docx), and uses the current sample inventory to state what is feasible. The sample catalogue is an input to this study; reconstruction and agent-benchmark experiments are outside this document's scope.

## 1. Biological question and competing explanations

An ART-related locus can contain an upstream noncoding array, an RT with an extended N-terminal region, and a separately encoded downstream partner. The NTD is part of the RT protein; it is neither the array nor the partner. The source study reported predicted RT–partner interfaces and array transcription, while leaving RT activity, physiological RNA substrates, RT–partner interaction, and biological function unresolved. See [Yoon et al., especially pp. 6, 13 and 37–38](<Yoon et al. 2026 Array Associated Reverse Transcriptases (2).pdf>).

We will distinguish four explanations that can coexist:

1. **Native-pair compatibility:** naturally co-occurring RTs and partners have features that support a preferential predicted interface.
2. **Broad family compatibility:** homologous partners receive similar support regardless of their locus of origin.
3. **Input or model dependence:** apparent discrimination arises from supplied MSA pairing, learned sequence associations, or prediction instability.
4. **Context-dependent association:** the interface involves the catalytic core, depends on intact RT architecture, or changes when array-derived RNA is included.

The primary outcome is computational discrimination between specified pairings. Natural co-occurrence is an observational pairing label, not an experimentally verified interaction label. A positive result would motivate binding experiments; a flat matrix would not establish absence of binding or prove that an earlier model is an artefact.

## 2. Available samples and present readiness

The versioned [sample dataset](../data/art_samples_2026-11-05/art_samples.csv) contains the following representatives:

| Quantity | Current count | Relevance to the binding study |
|---|---:|---|
| Deduplicated RT–partner candidates | 46 | Includes references, controls, and public acquisitions |
| Type I candidates | 36 | Initial family for panel selection |
| Complete type I pairs | 31 | Both coding sequences pass completeness/translation checks |
| Complete type I pairs passing the strict ART-profile screen | 18 | Starting shortlist; ART-versus-NART classification remains unresolved |
| Type II candidates | 10 | Potential separate replication arm; independence is not established |
| Verified NART references | 0 | Missing comparative controls |
| Compatible pairwise matches among strict complete type I candidates | 1 | Under the sequence-matching rules below |
| Eligible three-locus blocks | 0 | The confirmatory panel is not currently available |

These counts come from the [dataset summary](../data/art_samples_2026-11-05/summary.json) and [block-feasibility assessment](../results/art_expansion_20261105/step5/summary.json). Allowing tree-affinity review candidates still produced no matching triples. Partner clusters, accession counts, and complete CDSs do not establish independent biological samples.

**Present decision:** prepare the protocol and resolve sample eligibility. Any small prediction pilot must be labelled development or feasibility work and excluded from subsequent confirmation. The current inventory does not support a claim that the proposed specificity experiment is adequately powered.

## 3. Define and freeze an eligible protein panel

### 3.1 Locus-level requirements

Each selected locus must have source-linked RT and partner sequences, translation-matched complete CDSs, reviewed gene pairing, and reviewed family/domain assignments. Inspect frameshifts, ambiguous starts, alternative partner annotations, and contig-boundary effects before inclusion. Keep weak ART-family matches in a separate review/sensitivity panel.

Use sequence and genomic evidence to assess ART-family membership, including verified ART, NART, and competing RT references when available. An array, an extended NTD, or a generic RT motif alone cannot establish ART identity. Missing NART controls must remain an explicit family-classification limitation.

A detectable upstream array is not required for the protein-only comparison. RNA eligibility is assessed separately. Type I, II, and III partner families are not interchangeable controls; each primary matrix must contain homologous partners with compatible architecture.

### 3.2 Matched three-locus blocks

Select disjoint blocks of three loci using sequence and metadata only, before examining complex-prediction scores. Within each block, every pair must meet the proposed matching rules:

- 50–90% amino-acid identity in both RT core and partner.
- At least 80% reciprocal alignment coverage in both comparisons.
- Matching partner-domain architecture and a shorter/longer partner-length ratio of at least 0.8.
- Distinct source groups, with identical loci, gene-model aliases, and near-duplicate protein pairs removed.

These are design choices to obtain informative alternatives, not biological definitions of ART. Review supported RT-core and partner trees; do not assume that a 70%-identity cluster is an independent lineage. Avoid reusing RTs, partners, BioSamples, or near-duplicate loci across blocks, and retain broader clade relationships for the uncertainty analysis.

Target **8–12 evaluable blocks, requiring 24–36 eligible loci**. Eight is a minimum design gate, not proof of power. Type II can become a replication arm only if it independently meets the panel, calibration, and precision requirements. Type III remains descriptive unless a separate adequate design is established.

Freeze selected and reserve blocks, selection rules, exclusions, and the development/confirmation split. If no suitable blocks exist, acquire better-matched homologs or report a feasibility study. Do not fill blocks with distant proteins or change thresholds after seeing favourable interface scores.

## 4. Calibrate the prediction workflow before ART confirmation

Evaluate the fixed workflow on experimentally characterized cognate and tested alternative interactions from at least two homologous protein-interaction families. Record the experimental conditions and detection limits. Unmeasured alternative pairings are not negative controls, and weak non-cognate binding must not be relabelled as non-binding.

Apply the same input conditions, template policy, score definition, sampling, and aggregation planned for ART. Require a positive mean native-minus-tested-alternative contrast in each calibration family and input condition, and a combined 95% interval above zero from 10,000 resamples of experimental blocks within each fixed family. Group comparisons sharing endpoints together. Insufficient independent blocks make calibration inconclusive.

This establishes whether the workflow can resolve the intended kind of comparison; it does not establish a transferable ART effect size, eliminate training-data exposure, or convert confidence scores into affinities. Archive calibration failures and floor/ceiling effects.

Select a separate ART development panel by sequence and length strata to measure runtime, verify preprocessing, and test parsers. Exclude its evolutionary groups and close homologues from confirmation. Do not select an engine or score using whichever ART pilot produces the strongest diagonal.

## 5. Control the information supplied to the predictor

| Condition | Input treatment | Role |
|---|---|---|
| **U: fixed unpaired MSAs** | Independently prepared per-chain alignments; no cross-chain homologue pairing | Primary protein comparison |
| **S: single sequence** | Query sequences only | Mandatory sensitivity analysis |
| **P: naturally paired MSAs** | Explicit, documented locus pairing | Optional diagnostic of supplied pairing information |

Use a pinned Boltz-2 implementation as the initial candidate engine, subject to calibration and feature verification. Its documentation supports custom MSA input and an explicit single-sequence mode. It also documents an affinity module intended for protein–small-molecule prediction; that module will not be used to estimate RT–partner affinity. See the [official prediction documentation](https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md).

For U, prepare each chain's alignment once and reuse its sequence content, row order, depth, and preprocessing settings whenever that chain occurs. Remove shared pairing identifiers and verify that non-query homologues do not become cross-chain paired feature rows. Disable automatic multichain MSA generation. For S, verify that the processed features contain only the query sequences.

Disable templates and interaction restraints in all primary cells. Keep chain order, stoichiometry, recycling, and sampling settings constant. Use neutral identifiers rather than source/locus labels in prediction inputs. Archive raw inputs and processed features: choosing a prediction engine alone does not guarantee that pairing information is absent.

An engine that cannot implement or verify the required conditions is ineligible for confirmation. U and S control supplied information; they do not remove evolutionary signal from the sequences or unknown information acquired during model training.

## 6. Predict complete matrices and measure a balanced contrast

Each block contributes a complete 3 × 3 matrix under both U and S:

| | Partner 1 | Partner 2 | Partner 3 |
|---|---|---|---|
| RT 1 | Native | Alternative | Alternative |
| RT 2 | Alternative | Native | Alternative |
| RT 3 | Alternative | Alternative | Native |

Run three fixed independent seeds per cell and use the median score across seeds. Never select the best-scoring run. If the implementation produces additional samples per seed, fix their aggregation before prediction. Seeds measure prediction variability; they do not increase biological sample size.

The primary endpoint is the pinned engine's binary RT–partner interface-confidence score, with the exact output field and code version archived. For the planned Boltz-2 binary comparison, record the `iptm` field after confirming its definition in the pinned implementation. Treat it as a confidence readout, not binding probability or affinity. Raw scores from different engines are not pooled.

Secondary outputs are interchain predicted aligned error (PAE), interface contacts, local confidence, and seed agreement. Define a contact using non-hydrogen atoms within 5 Å and repeat summaries at 4 and 8 Å. Report both directions of interchain PAE and their mean over contacting residue pairs. A model with no contacts has zero contacts and undefined contact-PAE, rather than a favourable imputed value.

For block `g` and condition `c`, let `Y(g,i,j,c)` be the seed-median score. Define:

```text
D(g,c) = mean of the 3 native diagonal cells
         − mean of the 6 alternative off-diagonal cells
```

Each RT and each partner contributes equally to both terms, so additive row and column effects cancel. All remaining cell dependence stays inside the block. The primary estimand is the expected block contrast under U for the declared eligible-block selection scheme; S is the mandatory sensitivity counterpart. Report every block contrast and the equal-block mean.

This comparison evaluates native pairings against specified homologous alternatives. It does not isolate an effect beyond every form of sequence similarity or ancestry. In particular, do not regress off-diagonal scores on partner distance and extrapolate to the native pairing's distance of zero as if comparable alternative observations existed there.

## 7. Statistical analysis, failure handling, and decision gate

### 7.1 Uncertainty at the block and clade levels

Follow the v2 proposal's pre-specified block-level model separately for U and S:

```text
D(c) ~ MultivariateNormal(mu(c) × 1, sigma(c)^2 × R)
mu(c) ~ Normal(0, 0.2^2), restricted to [-1, 1]
sigma(c) ~ HalfNormal(0, 0.2^2)
R = 0.5 I + 0.5 K
```

Construct `K` before predictions from the mean of normalized block kernels derived from RT-core and partner patristic distances. For each tree, use `exp(-d / ell)`, where `ell` is its median nonzero distance among selected loci; average kernel entries over each pair of blocks and normalize to unit diagonal. Archive trees, alignments, branch-length treatment, and covariance code. Unsupported tree distances prevent confirmatory use of this model.

Report posterior medians and central 95% credible intervals, alongside raw block contrasts. Check convergence, effective sample sizes, and posterior predictive adequacy. Repeat with independent blocks (`R = I`), stronger dependence (`R = 0.25 I + 0.75 K`), each tree kernel alone, and prior scales 0.1 and 0.4. More than 5% predictive block contrasts outside the possible range [-1, 1], or material failure to reproduce the observed distribution, makes this model inadequate for confirmation.

Perform leave-one-block-out and leave-one-broader-clade-out analyses. A robust claim requires agreement across the registered dependence/prior sensitivities and no reversal of the observed mean direction on removing one group. Do not replace a failing model after viewing ART results and retain a confirmatory label.

Partner-label permutations within complete blocks provide a secondary reference, preserving endpoint reuse. Enumerate assignments when feasible or use 10,000 recorded permutations. Exchangeability is an assumption; these permutations are not an exact test of physical non-interaction or a cure for phylogenetic dependence.

### 7.2 Missing predictions

Retry a technical failure once using the same configuration and seed; retain all failure records. A primary cell requires all three seeds, and a primary block requires every planned U/S cell. Do not discard low scores or implausible interfaces as technical failures.

For unresolved failures, report complete-block results and worst-case contrast bounds using the score range [0, 1]. Replacements must come from the pre-ranked reserve. If failure handling changes the conclusion, or fewer than eight eligible blocks remain, report exploratory results only.

### 7.3 Precision before confirmation

Use a planning resolution target of **0.10 score units**, explicitly a computational target rather than a biological effect threshold. Before opening confirmatory ART scores, simulate the full decision procedure under plausible null and positive effects, phylogenetic dependence, heterogeneous block effects, seed variation, and missing outputs. Require at most 5% null positive decisions and at least 80% power at the planning target across the registered scenarios.

A strong positive interpretation requires successful calibration, adequate simulated resolution, positive U and S contrasts with 95% intervals above zero, and agreement with the registered robustness/failure checks. If those requirements fail, report the narrower condition-dependent or inconclusive result.

An interval entirely within [-0.10, 0.10] can support a precise computational null only when calibration, resolution, and failure checks succeed. Even then, it does not exclude physical binding. Extra seeds cannot compensate for too few eligible evolutionary groups.

## 8. Attribute predicted contacts to the NTD

### 8.1 Define boundaries before examining interfaces

The dataset's current PF00078-based NTD lengths are operational annotations, not reviewed domain boundaries. Build a conserved-core reference profile with independently supported catalytic-domain annotation. Register the core match state used to assign each sequence's core start; define the NTD as the preceding sequence. Ambiguous or missing boundary states exclude a primary domain assignment.

Inspect full-length monomer models before opening complex scores. Flag boundaries cutting a helix/strand, uncertain core alignments, or extensions without a credible independently folded domain. Test the baseline NTD cut and fixed cuts 20 residues shorter and longer; identify the extended construct's inclusion of assigned core residues. Do not select the boundary giving the strongest native score.

### 8.2 Mandatory localization and domain controls

- Map full-length RT–partner contacts to the pre-defined NTD and core, accounting for domain size and solvent exposure.
- Repeat the complete pairing matrices under U for all three NTD constructs and a complementary core-only construct.
- Predict full RT and all four constructs as monomers under the registered unpaired input treatment.
- Assess construct confidence and retained-domain similarity to the full-length monomer, and report seed agreement.

The strongest computational NTD attribution combines reproducible NTD contacts in full-length models, native-pair discrimination across structurally credible NTD cuts, and weaker discrimination in a credible core-only comparison. This supports NTD-associated predicted contacts; it does not demonstrate autonomous NTD binding.

A distorted or low-confidence isolated NTD cannot refute NTD participation in the full-length protein. If domain controls cannot be completed, retain the full-length comparison and omit NTD-specific conclusions.

## 9. Additional controls and conditional RNA extension

### Protein controls

Use verified, source-linked NART RTs as comparative controls when obtained, matching length, architecture, and prediction quality where possible. Do not designate a sequence NART merely because its array was undetected. Cross-type pairings are broad compatibility diagnostics; they are not the primary specificity negatives. High NART scores or low cross-type scores cannot independently establish or falsify an ART interface.

An AlphaFold-Multimer comparison can be run on a subset selected before scores are seen, using the same input-information controls and its own calibration. Agreement across engines is sensitivity evidence; shared training or modelling biases can remain.

### RNA and ternary comparisons

This is a secondary, hypothesis-generating extension. Choose loci and RNA units by a fixed diversity rule before viewing protein-matrix scores. Use transcript-supported boundaries when available; otherwise register a small set of plausible cuts and label their uncertainty. Do not choose whichever RNA fragment gives the highest confidence.

First assess recovery of a known RT–nucleic-acid interface using the actual composition of an experimentally resolved assembly. A retron complex containing RNA–DNA products or additional components must not be silently simplified into an equivalent RT–RNA reference. If this calibration is unavailable or fails, ART RNA models remain illustrative.

Compare the same protein inputs with native array RNA, non-cognate natural RNA, length/composition-matched unrelated RNA, and multiple dinucleotide-preserving shuffles. Record differences in RNA folding. For ternary models, keep the RT–partner pair fixed while varying RNA, and summarize changes in pair-specific contacts and PAE. Report RT–RNA and partner–RNA interfaces separately.

Do not rank RNA versus protein binding using their confidence scores, or compare global binary and ternary scores as interchangeable measures. An electropositive patch, a conserved repeat stem, or RNA-associated contact rearrangement can guide a hypothesis but cannot establish ligand preference, RNA substrate identity, or a partner-release mechanism.

## 10. Resources, execution stages, and outputs

Measure runtime, peak memory, storage, and failure rates on short, median, and long development complexes before requesting the full compute allocation. The v1 estimate of 200–300 GPU-hours is not treated as a measured requirement.

At eight blocks, the planned mandatory protein inventory is:

| Component | Calculation | Seeded prediction runs |
|---|---|---:|
| Full-length U/S matrices | 8 blocks × 9 cells × 2 conditions × 3 seeds | 432 |
| Three NTD-boundary matrices under U | 8 × 9 × 3 boundaries × 3 seeds | 648 |
| Core-only matrices under U | 8 × 9 × 3 seeds | 216 |
| RT/construct monomers under U | 24 loci × 5 representations × 3 seeds | 360 |
| **Subtotal** | Before calibration, retries, and optional arms | **1,656** |

Estimate GPU-hours per stratum from measured minutes per seeded run, with explicit additional lines for calibration, development, failures, RNA, and other engines. Under resource pressure, defer optional paired-MSA diagnostics, additional engines, and RNA before reducing required U/S, calibration, complete-block, or uncertainty controls. Preserve domain controls whenever making NTD claims.

| Stage | Exit criterion |
|---|---|
| Sample readiness | Reviewed family/architecture assignments, eligible disjoint blocks, and declared exclusions |
| Calibration and development | Score sensitivity demonstrated; U/S features audited; runtime and failure handling measured |
| Analysis freeze | Protocol, selected/reserve blocks, seeds, boundaries, model/checkpoint hashes, simulation report, and budget archived |
| Mandatory protein predictions | All planned full-length and domain-control jobs accounted for, including failures |
| Locked analysis | Block-level effects, robustness checks, contact localization, and permissible conclusions reported |
| Conditional extensions | RNA or additional-engine results reported with their own calibration and limits |

All future execution must follow [ciai cluster.docx](<../ciai cluster.docx>) and [AGENTS.md](../AGENTS.md): computation through Slurm, appropriate CPU/GPU partitions and QoS, explicit resource requests, and compliance with account limits.

Release a versioned protocol; sample/block and exclusion manifests; sequences; alignments and trees; raw and processed prediction inputs; software/checkpoint hashes; all models, scores, contacts and failures; resource records; and the analysis scripts. Maintain separate development, confirmatory, and exploratory outputs. A future experimental collaboration could test direct cognate versus alternative binding, NTD contributions, and RT activity as distinct biological questions.

## 11. Interpretations fixed before prediction

| Outcome | Permissible conclusion |
|---|---|
| Robust positive U/S contrasts after calibration and readiness gates | The workflow distinguishes native pairs from the specified homologous alternatives; physical interaction and interface co-evolution remain unproven |
| Advantage only with naturally paired MSAs | Discrimination depends on supplied pairing information |
| U positive but S is insensitive or fails | The result is alignment-dependent or sensitivity-limited |
| Flat contrast with wide uncertainty or inadequate groups | Pairing discrimination is unresolved |
| Precise computational null meeting the registered bounds | The workflow excludes a contrast of the target size in this sampled setting; binding is not excluded |
| Stable full-length NTD contacts and credible boundary/core controls agree | The predicted interface is consistently associated with the NTD |
| NTD truncations are unstable or inconsistent | Domain attribution is unresolved; full-length participation remains possible |
| Controlled RNA changes reproduce across comparable configurations | An RNA-associated structural hypothesis is supported for further testing |

## 12. Improvements over the original binding plan

The core question is retained, while several stronger claims and design assumptions are revised:

- Replace an assumed well-powered 12–15-locus grid with matched blocks and explicit data/precision gates.
- Replace reliance on engine choice to avoid pairing leakage with verified unpaired and query-only input conditions.
- Replace a zero-distance regression extrapolation with balanced within-block contrasts and declared evolutionary uncertainty.
- Add empirical calibration, fixed seed aggregation, complete failure accounting, and measured compute requirements.
- Add independent NTD boundaries, boundary perturbations, core-only comparisons, and monomer integrity checks.
- Treat flat matrices, NART scores, isolated-domain failures, and RNA scores as bounded computational evidence rather than decisive mechanistic verdicts.

A weaker NTD–partner identity correlation than RT-core–partner correlation does not itself refute an NTD interaction. Likewise, native-pair discrimination does not by itself prove co-evolution. The study is designed to produce an auditable computational result with clear limits and useful priorities for subsequent experiments.
