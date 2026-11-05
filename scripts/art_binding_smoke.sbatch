#!/bin/bash
#SBATCH --job-name=art-binding-smoke
#SBATCH --partition=cscc-gpu-p
#SBATCH --qos=gpu-debug-qos
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu-05
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:25:00
#SBATCH --output=/nfs-stor/rafiq.ali/projects/CRISP/tmp/art-binding-smoke-%j.log
set -euo pipefail
cd /nfs-stor/rafiq.ali/projects/CRISP
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1
export WANDB_MODE=disabled
.venv-boltz/bin/python -u scripts/art_binding_run.py --cell-ids S_r1_p1 U_r1_p1
.venv-boltz/bin/python -u scripts/art_binding_analyse.py
