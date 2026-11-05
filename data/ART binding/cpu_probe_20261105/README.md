# CPU fallback timing probe — 5 November 2026

The CPU node can load and begin execution of the local Boltz-2 model. The bounded
timing probe ran as Slurm job **373536** on `cn-11`, using 32 CPU cores and 64 GB RAM.
Its per-job directory contains the frozen input, command, prediction log, model
outputs, and status. It exhausted its 20-minute model-execution budget without
producing a finished prediction. See `373536/status.json` for the recorded outcome.

The input is the same SA1 RT/partner sequence pair as `S_r1_p1` in the GPU pilot
(501 + 265 residues). It uses query-only features, one seed, no templates,
no interaction restraints, no external MSA service, and no API key.

This is a **timing/compatibility probe**, using 1 recycling step and 10 sampling
steps. The registered experiment uses 3 recycling steps and 200 sampling steps.
The probe's scores must not be substituted for full-settings experiment results.

CPU inference uses explicitly recorded FP32 arithmetic; the pinned native GPU
configuration uses BF16 mixed precision. Hardware/precision variants are kept
separate, and all cells in a future CPU matrix must use the same registered CPU
settings. Neither identical inputs nor a shared checkpoint guarantees identical
CPU/GPU prediction scores.

The local GPU smoke job 373504 subsequently completed two models, and the
remaining local job 373517 was cancelled after the hosted API succeeded. Four
full-settings hosted predictions completed in 21.2–26.6 seconds each; their
separate results are in `../api_pilot_20261105/`. CPU preprocessing is supported
and audited, but a full CPU prediction matrix has not been launched.
