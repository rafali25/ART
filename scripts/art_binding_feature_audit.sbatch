#!/bin/bash
#SBATCH --job-name=art-binding-feature-audit
#SBATCH --partition=cscc-cpu-p
#SBATCH --qos=cscc-cpu-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/nfs-stor/rafiq.ali/projects/CRISP/tmp/art-binding-feature-audit-%j.log
set -euo pipefail
cd /nfs-stor/rafiq.ali/projects/CRISP
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1
.venv-boltz/bin/python -m unittest discover -s scripts -p test_art_binding_controls.py -v
.venv-boltz/bin/python -u scripts/art_binding_feature_audit.py
