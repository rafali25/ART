# Hosted Boltz-2 development pilot — 5 November 2026

The user-provided key was used successfully with NVIDIA's hosted Boltz-2 service.
All **four query-only native/alternative pair predictions** completed, with three
recycling steps, 200 sampling steps, and one diffusion sample per request.
End-to-end latency was **21.2–26.6 seconds per request**;
two requests ran concurrently. Median server-reported runtime was
**10.5 seconds**. This is a measured pilot result,
not a guarantee for larger complexes, throughput, or future service queues.

Use `results.csv` or `results.tsv` for the four-cell table. Each cell directory
contains its exact request, full response, runtime/status records, parsed result,
and downloaded `model_0.cif`. Both chain sequences in every CIF were checked
against the submitted proteins. The credential is not stored in this dataset.

## Scope and interpretation

This is the SA1/MarsHill 2×2 type II development panel, not the full type I study.
MSAs contain only the query for each chain. No template, interaction restraint,
ligand, or affinity prediction was supplied. The four results provide technical
validation of remote execution and output retrieval; the pairs are close
relatives and alternative pairings have no experimental non-binding label.

The API endpoint does not expose fixed-seed control or the processed features
needed to verify our strict unpaired-MSA adapter. Consequently only the S/query
condition was submitted; these are not the registered fixed-seed GPU replicates.
Full PAE matrices are not available in the hosted JSON response. Local GPU/CPU
outputs are retained separately. Backend calibration is required before making
biological conclusions. Low scores, a flat matrix, or a small descriptive contrast
cannot establish absence of interaction or a prediction artefact.

The original slow CPU probe reached its 20-minute reduced-settings budget.
The local GPU smoke test separately completed two models. The queued remaining
local matrix was cancelled once the hosted service worked, preserving completed
outputs. API cost checks were omitted as instructed by the user.

References: [NVIDIA hosted endpoint](https://docs.api.nvidia.com/nim/reference/mit-boltz2-infer),
[model card](https://docs.api.nvidia.com/nim/reference/mit-boltz2), and
[NIM input/output documentation](https://docs.nvidia.com/nim/bionemo/boltz2/inference).
The model card identifies Boltz2 v2.2.1; an exact deployed checkpoint hash is not
exposed by the hosted response and is not claimed here.
