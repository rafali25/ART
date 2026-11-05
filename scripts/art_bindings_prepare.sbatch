#!/bin/bash
#SBATCH --job-name=art-bindings-prepare
#SBATCH --partition=cscc-cpu-p
#SBATCH --qos=cscc-cpu-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/nfs-stor/rafiq.ali/projects/CRISP/tmp/art-bindings-prepare-%j.log
set -euo pipefail
cd /nfs-stor/rafiq.ali/projects/CRISP
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1
.venv/bin/python -m unittest discover -s scripts -p art_bindings_tests.py -v
.venv/bin/python -u scripts/art_bindings_prepare.py
