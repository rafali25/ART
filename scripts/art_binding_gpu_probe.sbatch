#!/bin/bash
#SBATCH --job-name=art-binding-gpu-probe
#SBATCH --partition=cscc-gpu-p
#SBATCH --qos=gpu-debug-qos
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:03:00
#SBATCH --output=/nfs-stor/rafiq.ali/projects/CRISP/tmp/art-binding-gpu-probe-%j.log
set -euo pipefail
cd /nfs-stor/rafiq.ali/projects/CRISP
export PYTHONDONTWRITEBYTECODE=1
srun .venv-boltz/bin/python -c 'import json,os,torch; print(json.dumps({"torch":torch.__version__,"compiled_CUDA":torch.version.cuda,"cuda_available":torch.cuda.is_available(),"visible_device_count":torch.cuda.device_count(),"CUDA_VISIBLE_DEVICES":os.environ.get("CUDA_VISIBLE_DEVICES"),"SLURM_JOB_GPUS":os.environ.get("SLURM_JOB_GPUS"),"SLURM_STEP_GPUS":os.environ.get("SLURM_STEP_GPUS")})); print(torch.cuda.get_device_name(0) if torch.cuda.is_available() else "CUDA unavailable")'
