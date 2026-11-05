#!/bin/bash
#SBATCH --job-name=art-binding-cpu-probe
#SBATCH --partition=cscc-cpu-p
#SBATCH --qos=cscc-cpu-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=/nfs-stor/rafiq.ali/projects/CRISP/tmp/art-binding-cpu-probe-%j.log
set -euo pipefail
cd /nfs-stor/rafiq.ali/projects/CRISP
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=32
export MKL_NUM_THREADS=32
export OPENBLAS_NUM_THREADS=32
export WANDB_MODE=disabled
.venv-boltz/bin/python -u scripts/art_binding_cpu_probe.py
