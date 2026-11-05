#!/bin/bash
#SBATCH --job-name=art-bindings-report
#SBATCH --partition=cscc-cpu-p
#SBATCH --qos=cscc-cpu-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=01:00:00
#SBATCH --output=/nfs-stor/rafiq.ali/projects/CRISP/tmp/art-bindings-report-%j.log
set -euo pipefail
cd /nfs-stor/rafiq.ali/projects/CRISP
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1
if [ -x .venv-report/bin/python ]; then
    .venv-report/bin/python -u scripts/art_bindings_report.py --final
    .venv-report/bin/python -u scripts/art_bindings_verify.py
else
    .venv/bin/python -u scripts/art_bindings_report.py --final
    .venv/bin/python -u scripts/art_bindings_verify.py
fi
