"""Screen newly retrieved DNA with frozen profiles and native Linux tools.

This driver adapts executable paths in memory; baseline scripts/results stay intact.
It runs exploratory array discovery without borrowing any old cohort p-values.
"""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))

from Bio import SeqIO
import art_type1_pipeline as pipeline
from art_reconstruction_common import parse_domains, table, sha


class NativeRunner:
    def __init__(self, out):
        self.out = Path(out)
        self.binaries = {"RVT_1": ROOT / ".tools/art-reconstruction/PF00078.hmm"}

    def run(self, args, log, stdout=None, expected=()):
        args = [str(x) for x in args]
        started = time.perf_counter()
        status = None
        try:
            with Path(log).open("w") as error:
                if stdout:
                    with Path(stdout).open("w") as dest:
                        result = subprocess.run(args, stdout=dest, stderr=error, check=True)
                else:
                    result = subprocess.run(args, stdout=error, stderr=subprocess.STDOUT, check=True)
                status = result.returncode
            for p in expected:
                if not Path(p).exists() or not Path(p).stat().st_size:
                    raise ValueError("Missing expected output: " + str(p))
        finally:
            with (self.out / "commands.jsonl").open("a") as stream:
                stream.write(json.dumps(dict(args=args, log=str(log), stdout=str(stdout),
                    elapsed_seconds=time.perf_counter()-started, returncode=status)) + "\n")

    def hmm(self, hmm, fasta, stem, z=2076, report_e=.001):
        stem = Path(stem)
        started = time.perf_counter()
        self.run([ROOT / ".venv/bin/hmmsearch", "--cpu", "2", "--seed", "42",
            "-Z", str(z), "--domZ", str(z), "-E", str(report_e), "--domE", str(report_e),
            "--tblout", str(stem)+".tbl", "--domtblout", str(stem)+".domtbl", hmm, fasta],
            str(stem)+".log", expected=[str(stem)+".domtbl"])
        return parse_domains(str(stem)+".domtbl"), time.perf_counter()-started


def main():
    out = HERE / "processed"
    out.mkdir(exist_ok=True)
    sources = []
    for name, cid, role in [
        ("ERR3946871_166.current_matched.fna", "ERR3946871_166", "new_Logan_archive_lead"),
        ("ERR2847180_159.current_matched.fna", "ERR2847180_159", "new_Logan_archive_lead"),
        ("SRR7532897_75.current_matched.fna", "SRR7532897_75", "new_Logan_archive_lead"),
        ("MGYC001898389771.fna", "MGYC001898389771", "new_MGnify_deferred_partner"),
    ]:
        raw = HERE / name
        record = SeqIO.read(raw, "fasta")
        record.id = cid
        record.name = cid
        record.annotations.update(molecule_type="DNA", topology="linear")
        dest = HERE / (cid+".source.gb")
        SeqIO.write(record, dest, "genbank")
        if cid.startswith("MGYC"):
            request = json.loads((HERE / (cid+".request.json")).read_text())
            url, release = request["url"], "MGnify_live_metadata_20261104"
        else:
            request = json.loads((HERE / (cid.split("_")[0]+".download.json")).read_text())
            url, release = request["url"], "Logan_current_object;"+request["headers"].get("ETag", "")
        sources.append(dict(accession=cid, path=dest.relative_to(ROOT).as_posix(), format="genbank",
            file_sha256=sha(dest), sequence_sha256=hashlib.sha256(str(record.seq).upper().encode()).hexdigest(),
            source_url=url, source_version_kind="retrieval_date_ETag_and_sequence_checksum",
            database_release=release, role=role))
    table(out / "acquired_source_records.tsv", sources)
    cfg = json.loads((ROOT / "config/art_type1_mgnify.json").read_text())
    cfg.update(output=out.relative_to(ROOT).as_posix(), threads=2)
    (out / "run_config.json").write_text(json.dumps(cfg, indent=2)+"\n")
    pipeline.Runner = NativeRunner
    profile = ROOT / "results/art_type1_public/profiles/TypeI_two_seed_exploratory.hmm"
    pipeline.profile_setup = lambda out, cfg, run: profile
    original_run = subprocess.run

    def native_mafft(args, **kwargs):
        if args and str(args[0]) == "cmd.exe":
            args = [str(x) for x in args]
            target = Path(args[args.index("--out")+1])
            input_path = args[-1]
            command = [str(ROOT / ".venv/bin/mafft"), "--auto", "--thread", "1", "--inputorder", input_path]
            # The existing function supplies a log stream as stdout. Native MAFFT
            # writes aligned FASTA on stdout and diagnostic messages on stderr.
            with target.open("w") as aligned:
                result = original_run(command, stdout=aligned, stderr=kwargs.get("stdout"), check=True)
            with (out / "commands.jsonl").open("a") as stream:
                stream.write(json.dumps(dict(args=command, stdout=str(target), returncode=result.returncode))+"\n")
            return result
        return original_run(args, **kwargs)

    subprocess.run = native_mafft
    try:
        pipeline.process(cfg, out, include_control=False)
    finally:
        subprocess.run = original_run
    inputs = [Path(__file__), ROOT / "scripts/art_type1_pipeline.py", profile,
        ROOT / "results/art_discovery/profiles/ART_RT.hmm",
        ROOT / "results/art_discovery/profiles/Type_I_partner.hmm",
        ROOT / "results/art_discovery/profiles/Type_II_partner.hmm",
        ROOT / ".tools/art-reconstruction/PF00078.hmm",
        ROOT / "config/art_array_search.json", ROOT / "config/art_array_assessment.json"]
    (out / "input_sha256.json").write_text(json.dumps({p.relative_to(ROOT).as_posix():sha(p) for p in inputs}, indent=2)+"\n")


if __name__ == "__main__":
    main()
