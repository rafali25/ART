#!/bin/bash
#SBATCH --job-name=art-binding-setup
#SBATCH --partition=cscc-cpu-p
#SBATCH --qos=cscc-cpu-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/nfs-stor/rafiq.ali/projects/CRISP/tmp/art-binding-setup-%j.log
# Build the pinned Boltz-2 environment and fetch its public checkpoints.
# Network-bound; runs on a CPU node per AGENTS.md. No GPU is requested here.
set -euo pipefail
cd /nfs-stor/rafiq.ali/projects/CRISP
export PYTHONDONTWRITEBYTECODE=1
ENV=.venv-boltz
CACHE=.tools/boltz_cache
OUT="data/ART binding/env"
mkdir -p "$CACHE" "$OUT"
export PIP_CACHE_DIR="/tmp/${USER}-art-binding-pip-${SLURM_JOB_ID}"
trap 'rm -rf "$PIP_CACHE_DIR"' EXIT

echo "host=$(hostname) job=${SLURM_JOB_ID} start=$(date -Is)"
if [ ! -x "$ENV/bin/python" ]; then
  /usr/bin/python3.12 -m venv "$ENV"
fi
"$ENV/bin/python" -m pip install --upgrade pip wheel
"$ENV/bin/python" -m pip index versions boltz | tee "$OUT/boltz_pip_index.txt"
BOLTZ_VERSION=$("$ENV/bin/python" -m pip index versions boltz | head -1 | sed -E 's/.*\((.*)\).*/\1/')
echo "pinning boltz==${BOLTZ_VERSION}"
"$ENV/bin/python" -m pip install "boltz[cuda]==${BOLTZ_VERSION}" gemmi
"$ENV/bin/python" -m pip freeze > "$OUT/pip_freeze.txt"

"$ENV/bin/python" - "$CACHE" "$OUT" <<'EOF'
import hashlib, json, sys, importlib.metadata as md
from pathlib import Path
from boltz.main import download_boltz2
cache, out = Path(sys.argv[1]).resolve(), Path(sys.argv[2])
download_boltz2(cache)
files = {}
for p in sorted(cache.rglob("*")):
    if p.is_file():
        h = hashlib.sha256()
        with p.open("rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 24), b""):
                h.update(chunk)
        files[str(p.relative_to(cache))] = {"bytes": p.stat().st_size, "sha256": h.hexdigest()}
info = {
    "boltz_version": md.version("boltz"),
    "torch_version": md.version("torch"),
    "cache_dir": str(cache),
    "checkpoint_files": files,
}
(out / "boltz_provenance.json").write_text(json.dumps(info, indent=2) + "\n")
print(json.dumps(info, indent=2))
EOF
echo "end=$(date -Is)"
