"""Read CUDA driver diagnostics within the single GPU assigned by Slurm."""
import ctypes
import json
import os
import socket
import subprocess
import traceback
import torch

print(json.dumps(dict(host=socket.gethostname(),job=os.environ.get('SLURM_JOB_ID'),
    CUDA_VISIBLE_DEVICES=os.environ.get('CUDA_VISIBLE_DEVICES'),SLURM_JOB_GPUS=os.environ.get('SLURM_JOB_GPUS'),
    SLURM_STEP_GPUS=os.environ.get('SLURM_STEP_GPUS'),LD_LIBRARY_PATH=os.environ.get('LD_LIBRARY_PATH'),
    torch=torch.__version__,compiled_CUDA=torch.version.cuda)),flush=True)
subprocess.run(['nvidia-smi','--query-gpu=index,uuid,driver_version','--format=csv'],check=True)
driver=ctypes.CDLL('libcuda.so.1');result=driver.cuInit(0);error=ctypes.c_char_p()
driver.cuGetErrorString(result,ctypes.byref(error))
print('cuInit='+str(result)+' '+str(error.value),flush=True)
try:
    torch.cuda.init();print('torch CUDA initialized: '+torch.cuda.get_device_name(0),flush=True)
except Exception:traceback.print_exc()
print('CUDA available: '+str(torch.cuda.is_available()),flush=True)
