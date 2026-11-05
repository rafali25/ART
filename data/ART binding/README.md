# ART binding experiments — started 5 November 2026

The active experiment is a **technical development pilot**. Four hosted
predictions are complete in [api_pilot_20261105/README.md](api_pilot_20261105/README.md).
The separate local pilot is documented in
[pilot_20261105/README.md](pilot_20261105/README.md). Both use the existing
[binding framework](../../literature/ARTNT_binding_framework.md) and original
proposal. The linked `Testing Cognate Specificity in the ART N-Terminal Domain`
DOCX contains a title, date, and computational-study note only; it does not
provide additional experiment settings.

The first pilot freezes SA1 and MarsHill references in a 2×2 native/alternative
pairing matrix, under query-only (S) and explicitly unpaired (U) inputs. Eight
one-seed predictions are registered. This is for input validation, runtime,
scoring and pipeline checks. The two related references do not establish
independent biological replication or experimental binding labels.

## Completed preparation

- Isolated Boltz 2.2.1 / PyTorch 2.5.1+cu124 environment and public model assets:
  [environment provenance](env/boltz_provenance.json).
- Source-linked, immutable pilot configuration and matrix:
  [manifest](pilot_20261105/manifest.json), [cells](pilot_20261105/cells.tsv).
- Four meaningful input-control tests passed in the source-analysis and model
  environments. They verify chain separation, preservation of homolog/deletion
  content, query-only inputs, rejection of leakage, and serializable audits.
- All eight actual Boltz parsing/featurization cases passed:
  [feature audit](pilot_20261105/cpu_feature_audit/summary.json). U contains 19
  rows (one query and 18 chain-specific homolog rows); S contains only one query
  row. No non-query cross-chain rows or templates remain. Processed arrays are
  retained alongside this audit.

Boltz's native unpaired row construction can co-populate homolog rows from both
chains. The pilot uses an archived input adapter that explicitly separates
those rows. The same adapter must be used and calibrated in any later study.

## Execution record

| Slurm job | Work | Outcome |
|---|---|---|
| 373437 | Initial GPU hardware inventory | A100 40-GB hardware, driver 570.195.03 |
| 373440 | Isolated environment and public checkpoints | Completed; dependencies checked and asset hashes recorded |
| 373455 | Control tests, alignments, eight frozen inputs | Completed |
| 373461 | First GPU pilot launch | Stopped before prediction: allocated CUDA device unavailable |
| 373467 | CPU preprocessing audit | Exposed missing RDKit all-property pickle setting; corrected to match official CLI |
| 373479 / 373488 | GPU initialization diagnosis | `gpu-05` assigned Slurm device 3, while its driver inventory exposed only indices 0–2; CUDA initialization returned code 100 |
| 373480 | Full CPU input audit | Feature checks completed; NumPy JSON serialization issue corrected |
| 373489 | CPU feature audit and corrected code snapshots | Completed: all eight actual input-feature cases passed |
| 373504 | Short S/U smoke stage on another GPU node | Completed two local models; outputs retained separately |
| 373517 | Remaining local eight-cell pilot | Cancelled after hosted API inference succeeded |
| 373536 | Bounded CPU prediction probe | Reduced-settings probe exhausted its 20-minute timing budget |
| 373559 | Identify supplied credential's service | NVIDIA hosted-inference key recognized; value never logged |
| 373580 | Hosted query-only 2×2 matrix | Four predictions completed with no failures; 21.2–26.6 seconds end-to-end each |
| 373586 | Verify and export hosted results | Four structures and all eight protein chain sequences verified |

Raw logs are in the repository's `tmp/art-binding-*JOBID.log`; per-cell attempts
and status records live under `pilot_20261105/runs/`. Failed attempts are retained.
Current launch code is identified in
[active_code.json](pilot_20261105/active_code.json), with versioned snapshots.

The latest outputs are [results.csv](api_pilot_20261105/results.csv) and the
[downloadable API pilot bundle](api_pilot_20261105.zip). The hosted run uses the
full 3-recycle/200-step prediction settings, with one stochastic sample per cell.
The endpoint exposes neither fixed seed control nor the processed features needed
to verify our strict U adapter, so only query-only S was submitted remotely.
Full PAE matrices are not returned by its hosted JSON response. Local and remote
scores are not pooled. No allocations remain active at this handoff; completed
outputs and failed attempts are preserved.

## Credentials and interpretation

The local workflow required no API key. The user then authorized the key in
`boltz.md` for hosted Boltz-2 inference and explicitly requested that cost checks
be skipped. The key is read only by the request client, excluded from version
control, and absent from the result bundle. Remote predictions run through
NVIDIA's hosted endpoint; no local GPU allocation is needed for those calls.
Neither local nor remote prediction uses templates, interaction restraints, or
the small-molecule affinity head.

The current catalogue has zero eligible confirmatory three-locus blocks and no
verified NART comparison set. Empirical specificity calibration and power gates
remain unmet. No pilot confidence score should be reported as binding affinity,
ART validation, or evidence of cognate specificity. The future full study needs
the sample, calibration and statistical gates recorded in the binding framework.

All computation runs through Slurm under the supplied cluster policy. The GPU
diagnosis did not change visibility settings or use unallocated GPUs.
