"""Bounded, recorded public phmmer queries from three complete diverse pairs.

Exploratory lead acquisition only; no active profiles or acceptance rules change.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import io
import json
import re
import sys
import time
import urllib.parse

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from art_type1_public_search import RecordedHTTP, HMMER, MGNIFY, normalize_hmmer, locate_peptide


def rows(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def dump(name, data):
    (OUT / name).write_text(json.dumps(data, indent=2) + "\n")


def table(name, data):
    if not data:
        (OUT / name).write_text("status\n")
        return
    fields = list(dict.fromkeys(k for r in data for k in r))
    with (OUT / name).open("w", newline="") as f:
        w = csv.DictWriter(f, fields, delimiter="\t")
        w.writeheader()
        w.writerows(data)


def prepare():
    catalog = rows(ROOT / "results/art_type1_catalog/complete_pair_representatives.tsv")
    proteins = {r.id: str(r.seq) for r in SeqIO.parse(ROOT / "results/art_type1_catalog/primary_partners.faa", "fasta")}
    selected_sources = ["SRR8925777_19723", "ERR8055547_17635", "MGYC002687634053"]
    selected = [next(r for r in catalog if r["source_genome"] == name) for name in selected_sources]
    assert len({r["partner_cluster70"] for r in selected}) == 3
    queries = []
    for row in selected:
        assert row["pair_completeness"] == "both_CDS_complete"
        assert row["partner_translation_qc"] == "complete_translation_matches"
        folder = ROOT / row["result_folder"]
        source = next(r for r in rows(folder / "source_records.tsv") if r["accession"] == row["source_genome"])
        source_path = ROOT / source["path"]
        if source.get("file_sha256"):
            assert hashlib.sha256(source_path.read_bytes()).hexdigest() == source["file_sha256"]
        rec = SeqIO.read(source_path, source["format"])
        partner = next(x for x in json.loads(row["partner_coordinates"]) if x["protein_id"] == row["partner_id"])
        gene = rec.seq[partner["start0"]:partner["end0"]]
        if partner["strand"] == -1:
            gene = gene.reverse_complement()
        peptide = str(gene.translate(table=11, cds=True))
        assert peptide == proteins[row["locus_id"]]
        rt = json.loads(row["RT_coordinates"])
        dna = rec.seq[rt["start0"]:rt["end0"]]
        if rt["strand"] == -1:
            dna = dna.reverse_complement()
        expected = next(str(r.seq) for r in SeqIO.parse(ROOT / "results/art_type1_catalog/RT_proteins.faa", "fasta") if r.id == row["locus_id"])
        assert str(dna.translate(table=11, cds=True)) == expected
        queries.append(dict(query_id=row["source_genome"], family="Type_I_partner", sequence=peptide,
            source_locus=row["locus_id"], partner_id=row["partner_id"], cluster70=row["partner_cluster70"],
            sequence_sha256=hashlib.sha256(peptide.encode()).hexdigest(), source_path=source["path"],
            source_file_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
            full_partner_and_RT_translation_verified=True, length_aa=len(peptide)))
    dump("queries.json", queries)
    SeqIO.write([SeqRecord(Seq(q["sequence"]), id=q["query_id"], description=q["partner_id"]) for q in queries], OUT / "queries.faa", "fasta")
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    comparisons = []
    anchors = [r for r in catalog if r["source_genome"] == "MZ422438.1"] + selected[:2]
    for row in [r for r in catalog if r["source_genome"].startswith("MGYC")]:
        for other in anchors:
            a, b = proteins[row["locus_id"]], proteins[other["locus_id"]]
            ident = aligner.align(a, b)[0].counts().identities / max(len(a), len(b))
            comparisons.append(dict(candidate=row["source_genome"], comparison=other["source_genome"] + "|" + other["partner_id"],
                global_identity_over_longer_length=ident, selected=row["source_genome"] == selected_sources[2]))
    table("seed_diversity.tsv", comparisons)
    dump("run_config.json", dict(database="mgnify30_c2", expected_release="2026_07", max_jobs=3,
        max_hits_per_query=50, remote_report_evalue=0.001, max_metadata_proteins=5, max_contigs=3,
        max_source_bytes=5000000, purpose="exploratory_public_leads_no_profile_training_or_final_ART_claim",
        seed_selection="Two public archive representatives in distinct 70pct clusters; MGYC representative minimizing maximum global identity to both LPJP1 controls and the first two queries"))
    print(json.dumps([{k:v for k,v in q.items() if k != "sequence"} for q in queries], indent=2))


def submit():
    cfg = json.loads((OUT / "run_config.json").read_text())
    http = RecordedHTTP(OUT / "remote", timeout=30)
    databases = http.json("databases", HMMER + "/search/databases", limit=1000000)
    db = next(d for d in databases if d["id"] == cfg["database"] and d["status"] == "enabled")
    dump("database.json", db)
    if db["version"] != cfg["expected_release"]:
        raise ValueError("Database release differs from frozen configuration")
    jobs = []
    for q in json.loads((OUT / "queries.json").read_text()):
        payload = dict(input=">" + q["query_id"] + "\n" + q["sequence"] + "\n", database=db["id"],
            E=cfg["remote_report_evalue"], domE=cfg["remote_report_evalue"], incE=1e-5, incdomE=1e-5)
        try:
            job = http.json("submit_" + q["query_id"], HMMER + "/search/phmmer", payload, limit=1000000)
            jobs.append(dict(query_id=q["query_id"], family=q["family"], job_id=job["id"], status="submitted"))
        except Exception as exc:
            jobs.append(dict(query_id=q["query_id"], family=q["family"], status="submission_failed", error=str(exc)))
        dump("jobs.json", jobs)
    print(json.dumps(jobs, indent=2))


def poll():
    cfg = json.loads((OUT / "run_config.json").read_text())
    http = RecordedHTTP(OUT / "remote", timeout=30)
    jobs = json.loads((OUT / "jobs.json").read_text())
    for j in jobs:
        if not j.get("job_id") or j["status"] == "SUCCESS":
            continue
        url = HMMER + "/result/" + j["job_id"] + "?" + urllib.parse.urlencode(dict(page=1, page_size=cfg["max_hits_per_query"], with_domains="true"))
        name = "result_" + j["query_id"] + "_" + str(int(time.time()))
        try:
            data = http.json(name, url, cached=False, limit=4000000)
            j["status"] = data.get("status", "unknown")
            j["result_path"] = str((OUT / "remote" / (name + ".response")).relative_to(ROOT))
        except Exception as exc:
            j["last_error"] = str(exc)
        dump("jobs.json", jobs)
    print(json.dumps(jobs, indent=2))


def collect():
    old = set()
    for cohort in ("art_type1_public", "art_type1_mgnify"):
        old.update(r["protein_accession"] for r in rows(ROOT / "results" / cohort / "remote_hits.tsv"))
    hits = []
    summaries = []
    for j in json.loads((OUT / "jobs.json").read_text()):
        if j["status"] != "SUCCESS":
            summaries.append(dict(query=j["query_id"], status=j["status"]))
            continue
        data = json.loads((ROOT / j["result_path"]).read_text())
        found = normalize_hmmer(data, j)
        assert len(found) <= 50, "Server returned more than requested bounded hits"
        for hit in found:
            hit["in_original_saved_hit_lists"] = hit["protein_accession"] in old
            hit["evidence_level"] = "protein_homology_lead_only"
        hits += found
        summaries.append(dict(query=j["query_id"], status=j["status"], returned_hits=len(found),
            incremental_unique_accessions=len({r["protein_accession"] for r in found} - old),
            page_count=data.get("page_count"), result_keys=list((data.get("result") or {}).keys())))
    table("remote_hits.tsv", hits)
    dump("search_summary.json", dict(queries=summaries, returned_rows=len(hits),
        unique_accessions=len({r["protein_accession"] for r in hits}),
        incremental_unique_accessions=len({r["protein_accession"] for r in hits} - old),
        complete_ART_candidates=0, comparison="Both original public and MGnify remote_hits.tsv accession lists; bounded top-50 results, not exhaustive database absence"))
    print((OUT / "search_summary.json").read_text())


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("operation", choices=("prepare", "submit", "poll", "collect"))
    globals()[p.parse_args().operation]()
