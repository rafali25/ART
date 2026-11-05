#!/bin/bash
#SBATCH --job-name=art-binding-environment
#SBATCH --partition=cscc-cpu-p
#SBATCH --qos=cscc-cpu-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/nfs-stor/rafiq.ali/projects/CRISP/tmp/art-binding-environment-%j.log
set -euo pipefail
cd /nfs-stor/rafiq.ali/projects/CRISP
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1
export PIP_CACHE_DIR="/tmp/${USER}-art-binding-pip-${SLURM_JOB_ID}"
trap 'rm -rf "$PIP_CACHE_DIR"' EXIT
mkdir -p 'data/ART binding/env' .tools/boltz_cache
if [ ! -x .venv-boltz/bin/python ]; then
    /usr/bin/python3.12 -m venv .venv-boltz
fi
.venv-boltz/bin/python -m pip install --upgrade pip wheel
.venv-boltz/bin/python -m pip install --index-url https://download.pytorch.org/whl/cu124 torch==2.5.1
.venv-boltz/bin/python -m pip install -r requirements-art-binding.txt
.venv-boltz/bin/python -m pip check
.venv-boltz/bin/python -m pip freeze > 'data/ART binding/env/pip_freeze.txt'
.venv-boltz/bin/python -u scripts/art_binding_environment.py
