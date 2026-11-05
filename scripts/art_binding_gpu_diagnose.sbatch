#!/bin/bash
#SBATCH --job-name=art-binding-cuda-diagnose
#SBATCH --partition=cscc-gpu-p
#SBATCH --qos=gpu-debug-qos
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:03:00
#SBATCH --output=/nfs-stor/rafiq.ali/projects/CRISP/tmp/art-binding-cuda-diagnose-%j.log
set -euo pipefail
cd /nfs-stor/rafiq.ali/projects/CRISP
export PYTHONDONTWRITEBYTECODE=1
srun .venv-boltz/bin/python -u scripts/art_binding_gpu_diagnose.py
