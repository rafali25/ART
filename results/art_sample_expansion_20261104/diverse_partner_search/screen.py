"""Native exploratory screening of newly acquired diverse-query contigs.

Uses frozen baseline family profiles and thresholds without retraining.
"""
from pathlib import Path
import importlib.util
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))
import art_type1_pipeline as pipeline
from art_reconstruction_common import table, sha
from extract_art_loci import read_tsv

DRIVER = ROOT / "results/art_sample_expansion_20261104/deferred_acquisition/screen_deferred.py"
spec = importlib.util.spec_from_file_location("native_deferred_driver", DRIVER)
driver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(driver)


def main():
    out = HERE / "processed"
    out.mkdir(exist_ok=True)
    sources = read_tsv(HERE / "acquired_source_records.tsv")
    table(out / "acquired_source_records.tsv", sources)
    cfg = json.loads((ROOT / "config/art_type1_mgnify.json").read_text())
    cfg.update(output=out.relative_to(ROOT).as_posix(), threads=2)
    (out / "run_config.json").write_text(json.dumps(cfg, indent=2) + "\n")
    pipeline.Runner = driver.NativeRunner
    profile = ROOT / "results/art_type1_public/profiles/TypeI_two_seed_exploratory.hmm"
    pipeline.profile_setup = lambda out, cfg, run: profile
    original_run = subprocess.run

    def native_mafft(args, **kwargs):
        if args and str(args[0]) == "cmd.exe":
            args = [str(x) for x in args]
            target = Path(args[args.index("--out") + 1])
            command = [str(ROOT / ".venv/bin/mafft"), "--auto", "--thread", "1", "--inputorder", args[-1]]
            with target.open("w") as aligned:
                result = original_run(command, stdout=aligned, stderr=kwargs.get("stdout"), check=True)
            with (out / "commands.jsonl").open("a") as stream:
                stream.write(json.dumps(dict(args=command, stdout=str(target), returncode=result.returncode)) + "\n")
            return result
        return original_run(args, **kwargs)

    subprocess.run = native_mafft
    try:
        pipeline.process(cfg, out, include_control=False)
    finally:
        subprocess.run = original_run
    inputs = [Path(__file__), DRIVER, ROOT / "scripts/art_type1_pipeline.py", profile,
        ROOT / "results/art_discovery/profiles/ART_RT.hmm",
        ROOT / "results/art_discovery/profiles/Type_I_partner.hmm",
        ROOT / "results/art_discovery/profiles/Type_II_partner.hmm",
        ROOT / ".tools/art-reconstruction/PF00078.hmm",
        ROOT / "config/art_array_search.json", ROOT / "config/art_array_assessment.json"]
    inputs += [ROOT / r["path"] for r in sources]
    (out / "input_sha256.json").write_text(json.dumps({str(p.relative_to(ROOT)):sha(p) for p in inputs}, indent=2) + "\n")


if __name__ == "__main__":
    main()
