"""Verify bounded acquisition, exact translations, and incremental diversity."""
from pathlib import Path
import csv
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices
from art_reconstruction_common import table
from extract_art_loci import read_tsv


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    processed = HERE / "processed"
    oldrows = read_tsv(ROOT / "results/art_type1_catalog/typeI_samples.tsv")
    oldrt = {r.id:str(r.seq) for r in SeqIO.parse(ROOT / "results/art_type1_catalog/RT_proteins.faa", "fasta")}
    oldpt = {r.id:str(r.seq) for r in SeqIO.parse(ROOT / "results/art_type1_catalog/primary_partners.faa", "fasta")}
    newrt = {r.id:str(r.seq) for r in SeqIO.parse(processed / "candidate_RT_proteins.faa", "fasta")}
    newpt = {r.id:str(r.seq) for r in SeqIO.parse(processed / "typeI_partners.faa", "fasta")}
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    comparisons = []
    for pair in read_tsv(processed / "typeI_pairs.tsv"):
        rt, pt = newrt[pair["RT_key"]], newpt[pair["partner_key"]]
        for old in oldrows:
            key = old["locus_id"]
            x,y = oldrt[key],oldpt[key]
            comparisons.append(dict(new_locus=pair["locus_id"], prior_locus=key,
                RT_exact_sequence_match=rt==x, partner_exact_sequence_match=pt==y,
                RT_global_identity_over_longer_length=aligner.align(rt,x)[0].counts().identities/max(len(rt),len(x)),
                partner_global_identity_over_longer_length=aligner.align(pt,y)[0].counts().identities/max(len(pt),len(y))))
    comparisons.sort(key=lambda r:-r["partner_global_identity_over_longer_length"])
    table(HERE / "prior_catalog_similarity.tsv", comparisons)
    jobs = json.loads((HERE / "jobs.json").read_text())
    summary = json.loads((HERE / "acquisition_summary.json").read_text())
    assert len(jobs) == 3 and all(r["status"]=="SUCCESS" for r in jobs)
    assert summary["metadata_proteins_attempted"]<=5 and summary["source_contigs"]<=3 and summary["source_bytes"]<=5000000
    sources = read_tsv(HERE / "acquired_source_records.tsv")
    for source in sources:
        assert digest(ROOT/source["path"]) == source["file_sha256"]
    for path,expected in json.loads((processed / "input_sha256.json").read_text()).items():
        assert digest(ROOT/path) == expected
    for q in json.loads((HERE / "queries.json").read_text()):
        assert digest(ROOT/q["source_path"]) == q["source_file_sha256"]
        assert hashlib.sha256(q["sequence"].encode()).hexdigest() == q["sequence_sha256"]
    pair_checks = []
    for candidate in read_tsv(processed / "candidate_loci.tsv"):
        record = SeqIO.read(ROOT/candidate["source_path"],"fasta")
        rtloc = json.loads(candidate["rt_coordinates_and_strand"])
        nt = record.seq[rtloc["start0"]:rtloc["end0"]]
        if rtloc["strand"] == -1:
            nt = nt.reverse_complement()
        assert str(nt.translate(table=11,cds=True)) == newrt[candidate["protein_key"]]
        for p in json.loads(candidate["partner_family_and_coordinates"]):
            nt = record.seq[p["start0"]:p["end0"]]
            if p["strand"] == -1:
                nt = nt.reverse_complement()
            assert str(nt.translate(table=11,cds=True)) == newpt[p["protein_key"]]
        pair_checks.append(candidate["locus_id"])
    verification = dict(all_three_live_jobs_successful=True, budgets_respected=True,
        source_hashes_verified=True, screening_input_hashes_verified=True,
        seed_source_hashes_verified=True, complete_pair_translations_replayed=pair_checks,
        supported_complete_ART_loci=0, new_independent_lineages_established=0,
        inference="One generic-RT/type-I-partner context pair is new to the source catalog but closely matches an existing pair; ART coverage fails frozen screen")
    (HERE / "verification.json").write_text(json.dumps(verification,indent=2)+"\n")
    inputs = [ROOT/"scripts/art_type1_public_search.py",ROOT/"results/art_type1_catalog/typeI_samples.tsv",
        ROOT/"results/art_type1_catalog/RT_proteins.faa",ROOT/"results/art_type1_catalog/primary_partners.faa",
        ROOT/"results/art_type1_public/remote_hits.tsv",ROOT/"results/art_type1_mgnify/remote_hits.tsv"]
    (HERE / "input_sha256.json").write_text(json.dumps({str(p.relative_to(ROOT)):digest(p) for p in inputs},indent=2)+"\n")
    (HERE / "output_sha256.json").write_text(json.dumps({str(p.relative_to(HERE)):digest(p) for p in sorted(HERE.rglob("*"))
        if p.is_file() and p.name!="output_sha256.json" and "__pycache__" not in p.parts},indent=2)+"\n")
    print(json.dumps(verification,indent=2))


if __name__ == "__main__":
    main()
