#!/bin/bash
#SBATCH --job-name=array-finders-setup
#SBATCH --partition=cscc-cpu-p
#SBATCH --qos=cscc-cpu-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/nfs-stor/rafiq.ali/projects/CRISP/tmp/array-finders-setup-%j.log
# Install the baseline array finders used by the PILER-CR shuffle-wrapper comparison
# into a dedicated environment, .venv-array-finders/: PILER-CR 1.06 (C++ source,
# built here) and MinCED 0.4.2 (release jar plus a bin/minced wrapper).
set -euo pipefail
cd /nfs-stor/rafiq.ali/projects/CRISP
ENV=.venv-array-finders
[ -x "$ENV/bin/python" ] || python3 -m venv "$ENV"
mkdir -p "$ENV/downloads" "$ENV/share/minced"
curl -sSfL -o "$ENV/downloads/pilercr1.06.tar.gz" https://drive5.com/pilercr/pilercr1.06.tar.gz
curl -sSfL -o "$ENV/downloads/minced-0.4.2.jar" https://github.com/ctSkennerton/minced/releases/download/0.4.2/minced.jar
BUILD=$(mktemp -d "/tmp/array-finders-${SLURM_JOB_ID}-XXXX")
trap 'rm -rf "$BUILD"' EXIT
tar -xzf "$ENV/downloads/pilercr1.06.tar.gz" -C "$BUILD"
SRC=$(dirname "$(find "$BUILD" -name Makefile | head -1)")
echo "PILER-CR source: $SRC"; cat "$SRC/Makefile"
( cd "$SRC" && make -j4 ) || ( cd "$SRC" && make clean >/dev/null 2>&1 || true; make -j4 CFLAGS="-O3 -DNDEBUG=1 -fpermissive -w" )
cp "$SRC/pilercr" "$ENV/bin/pilercr"
cp "$ENV/downloads/minced-0.4.2.jar" "$ENV/share/minced/minced.jar"
cat > "$ENV/bin/minced" <<'EOF'
#!/bin/bash
exec java -Xmx2g -jar "$(dirname "$(readlink -f "$0")")/../share/minced/minced.jar" "$@"
EOF
chmod +x "$ENV/bin/minced"
"$ENV/bin/pilercr" 2>&1 | head -5 || true
"$ENV/bin/minced" --version 2>&1 | head -3 || true
# Smoke test on one published reference locus (MarsHill); output kept for the parser.
SMOKE="$ENV/smoke"; mkdir -p "$SMOKE"
cp results/oriented_art_loci/loci/MW248466.1__QQM14740.1.fna "$SMOKE/in.fna"
cp results/art_validation/loci/repetitive_DNA_1.fna "$SMOKE/repetitive.fna"
for name in in repetitive; do
  "$ENV/bin/pilercr" -in "$SMOKE/$name.fna" -out "$SMOKE/$name.pilercr_default.txt" -quiet -noinfo || true
  "$ENV/bin/pilercr" -in "$SMOKE/$name.fna" -out "$SMOKE/$name.pilercr_tuned.txt" -quiet -noinfo -minrepeat 15 -maxrepeat 45 -minspacer 25 -maxspacer 420 || true
  "$ENV/bin/minced" "$SMOKE/$name.fna" "$SMOKE/$name.minced_default.txt" || true
  "$ENV/bin/minced" -minRL 15 -maxRL 45 -minSL 25 -maxSL 420 "$SMOKE/$name.fna" "$SMOKE/$name.minced_tuned.txt" || true
done
grep -H "putative CRISPR" "$SMOKE"/*.pilercr_*.txt || true
grep -c "^CRISPR" "$SMOKE"/*.minced_*.txt || true
"$ENV/bin/python" - <<'EOF'
import hashlib, json, platform, subprocess
from pathlib import Path
env = Path('.venv-array-finders')
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
record = dict(
    pilercr=dict(version='1.06', source_url='https://drive5.com/pilercr/pilercr1.06.tar.gz',
                 source_sha256=digest(env/'downloads/pilercr1.06.tar.gz'), binary_sha256=digest(env/'bin/pilercr'),
                 compiler=subprocess.run(['g++', '--version'], capture_output=True, text=True).stdout.splitlines()[0]),
    minced=dict(version='0.4.2', source_url='https://github.com/ctSkennerton/minced/releases/download/0.4.2/minced.jar',
                jar_sha256=digest(env/'share/minced/minced.jar'),
                java=subprocess.run(['java', '-version'], capture_output=True, text=True).stderr.splitlines()[0]),
    host=platform.node())
(env/'provenance.json').write_text(json.dumps(record, indent=2)+'\n')
print(json.dumps(record, indent=2))
EOF
