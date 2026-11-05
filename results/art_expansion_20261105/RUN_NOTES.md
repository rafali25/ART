# Five-step expansion execution notes

The user authorized completing all five immediate execution steps in the
2026-11-04 PLAN, with a checkpoint after every step, and requires compliance with
`ciai cluster.docx`. `AGENTS.md` persists the applicable execution rules.

All downloads, parsing, gene calling, searches, tests and phylogenetic work run
in Slurm CPU batch jobs (`cscc-cpu-p`, `cscc-cpu-qos`). Login-node work is file
preparation, reading logs and scheduler administration only. No GPUs, salloc,
exclusive allocations, idle allocation jobs or external messages are used.

Execution jobs submitted:

- 371266: checkpoint 1 passed; initial HPGC download passed checksum, then a
  request-log formatting error stopped further preparation. The error was fixed.
- 371284: reused verified HPGC; JPGD downloaded and passed published MD5.
- 371302 and 371327: archive search superseded after preserving cached prefixes
  and mappings, to replace slow per-contig translation prefiltering and serial
  transfers. Only these task-owned jobs were cancelled.
- 371364: resumed 45-lead/44-source-run archive and bounded MGnify occurrence
  acquisition, then native two-route screening and checkpoint 2.
- 371321: full HPGC then nonoverlapping JPGD screens, with optional bounded UHGV
  selection; checkpoints 3 and 4. Runs alongside the network-bound archive job.
- 371381: validation, independent RT-core/type-I-partner trees, matching/block
  feasibility and checkpoint 5, dependent on both 371425 and 371321 succeeding.
- 371425: after 371364, a four-representative MGnify Parquet member/occurrence
  benchmark through a proxy enforcing a 100-MB remote-read cap; updates checkpoint
  2 and screens any additional exactly mapped source contigs.

Source code: `scripts/art_expansion_workflow.py`, `art_expansion_leads.py`,
`art_expansion_bulk.py`, `art_expansion_extra.py`, and
`art_expansion_feasibility.py`, plus corresponding Slurm scripts and mapping tests.
No previous scientific cohort/profile/null distribution is changed.

Durable checkpoints are `checkpoints/01.json` through `05.json` and Markdown
companions. Repeated checkpoint versions are preserved in `checkpoints/history/`.
The presence of a checkpoint records completed bounded work and its scientific
limitations; it does not imply the study is adequately powered or validated.

Slurm logs are `tmp/art-{prepare,download,leads,catalogs,feasibility}-JOBID.log`.
Use `squeue`/`sacct` and these logs to resume; never execute the analysis directly
on ciai-login-1. Analysis drivers reject execution without SLURM_JOB_ID.


Final closure:

- 371381 passed all 12 tests and source-map replay but stopped because MMseqs
  rejected the spelling `.7`/`.8` for numeric thresholds. Arguments were changed
  to `0.7`/`0.8`, with the same biological thresholds.
- 371557 completed the rerun and saved checkpoint 5. No source acquisition or
  bulk gene searches were repeated for this correction.
- The bounded Parquet member attempt in 371425 failed on an official endpoint
  connection refusal, transferring zero Parquet bytes. No member absence is
  inferred. Source-run ENA metadata retrieval and archive screening completed.
- Final audit runs as a separate CPU Slurm job, preserves a source snapshot,
  rechecks all checkpoint artifacts and original baseline hashes, and writes
  `verification.json`. No active analysis jobs remain after closure.
- `sacct` could not reach its accounting database. This is recorded rather than
  substituting invented accounting statistics.
