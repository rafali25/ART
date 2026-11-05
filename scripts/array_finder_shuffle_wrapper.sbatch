#!/bin/bash
#SBATCH --job-name=pilercr-shuffle-wrapper
#SBATCH --partition=cscc-cpu-p
#SBATCH --qos=cscc-cpu-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/nfs-stor/rafiq.ali/projects/CRISP/tmp/pilercr-shuffle-wrapper-%j.log
# Part B: tuned PILER-CR inside the selection-aware shuffle wrapper on the 31 inputs.
# Requires .venv-array-finders/ from scripts/array_finders_setup.sbatch.
set -euo pipefail
cd /nfs-stor/rafiq.ali/projects/CRISP
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export TMPDIR=$(mktemp -d "/tmp/pilercr-wrapper-${SLURM_JOB_ID}-XXXX")
trap 'rm -rf "$TMPDIR"' EXIT
.venv/bin/python -u scripts/array_finder_shuffle_wrapper.py --workers "$SLURM_CPUS_PER_TASK"
