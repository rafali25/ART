#!/bin/bash
#SBATCH --job-name=art-bindings-report-env
#SBATCH --partition=cscc-cpu-p
#SBATCH --qos=cscc-cpu-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/nfs-stor/rafiq.ali/projects/CRISP/tmp/art-bindings-report-env-%j.log
set -euo pipefail
cd /nfs-stor/rafiq.ali/projects/CRISP
export PYTHONDONTWRITEBYTECODE=1
export PIP_CACHE_DIR="/tmp/${USER}-art-report-pip-${SLURM_JOB_ID}"
trap 'rm -rf "$PIP_CACHE_DIR"' EXIT
if [ ! -x .venv-report/bin/python ]; then /usr/bin/python3.12 -m venv .venv-report; fi
.venv-report/bin/python -m pip install numpy==1.26.4 biopython==1.84 matplotlib==3.10.7
.venv-report/bin/python -m pip check
.venv-report/bin/python -m pip freeze > data/results_bindings/provenance/report_environment.txt
.venv-report/bin/python -m py_compile scripts/art_bindings_report.py
