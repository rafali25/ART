#!/bin/bash
#SBATCH --job-name=art-bindings-exploratory
#SBATCH --partition=cscc-cpu-p
#SBATCH --qos=cscc-cpu-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=24:00:00
#SBATCH --output=/nfs-stor/rafiq.ali/projects/CRISP/tmp/art-bindings-exploratory-%j.log
set -euo pipefail
cd /nfs-stor/rafiq.ali/projects/CRISP
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1
.venv/bin/python -u scripts/art_bindings_run.py --phase full --retry-gateway-failures
