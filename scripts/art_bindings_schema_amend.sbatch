#!/bin/bash
#SBATCH --job-name=art-bindings-schema-fix
#SBATCH --partition=cscc-cpu-p
#SBATCH --qos=cscc-cpu-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=/nfs-stor/rafiq.ali/projects/CRISP/tmp/art-bindings-schema-fix-%j.log
set -euo pipefail
cd /nfs-stor/rafiq.ali/projects/CRISP
export PYTHONDONTWRITEBYTECODE=1
.venv/bin/python -m unittest discover -s scripts -p art_bindings_tests.py -v
.venv/bin/python -u scripts/art_bindings_schema_amend.py
