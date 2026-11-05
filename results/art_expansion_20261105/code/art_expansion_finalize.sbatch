#!/bin/bash
#SBATCH --job-name=art-feasibility
#SBATCH --partition=cscc-cpu-p
#SBATCH --qos=cscc-cpu-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/nfs-stor/rafiq.ali/projects/CRISP/tmp/art-feasibility-%j.log
set -euo pipefail
cd /nfs-stor/rafiq.ali/projects/CRISP
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1
.venv/bin/python -m unittest discover -s scripts -p test_art_expansion_workflow.py -v
.venv/bin/python -m unittest discover -s scripts -p test_art_type1_public.py -v
.venv/bin/python -u scripts/art_expansion_workflow.py feasibility
