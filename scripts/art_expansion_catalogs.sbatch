#!/bin/bash
#SBATCH --job-name=art-catalogs
#SBATCH --partition=cscc-cpu-p
#SBATCH --qos=cscc-cpu-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/nfs-stor/rafiq.ali/projects/CRISP/tmp/art-catalogs-%j.log
set -euo pipefail
cd /nfs-stor/rafiq.ali/projects/CRISP
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1
.venv/bin/python -u scripts/art_expansion_workflow.py HPGC
.venv/bin/python -u scripts/art_expansion_workflow.py JPGD
