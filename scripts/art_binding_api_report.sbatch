#!/bin/bash
#SBATCH --job-name=art-binding-api-report
#SBATCH --partition=cscc-cpu-p
#SBATCH --qos=cscc-cpu-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/nfs-stor/rafiq.ali/projects/CRISP/tmp/art-binding-api-report-%j.log
set -euo pipefail
cd /nfs-stor/rafiq.ali/projects/CRISP
export PYTHONDONTWRITEBYTECODE=1
.venv/bin/python -u scripts/art_binding_api_report.py
