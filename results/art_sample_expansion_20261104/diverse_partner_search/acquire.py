"""Bounded metadata/contig acquisition for incremental diverse-query hits."""
from pathlib import Path
import csv
import hashlib
import io
import json
import sys
import urllib.parse

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from art_type1_public_search import RecordedHTTP, MGNIFY, locate_peptide


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def dump(name, data):
    (OUT / name).write_text(json.dumps(data, indent=2) + "\n")


def table(name, data):
    fields = list(dict.fromkeys(k for r in data for k in r)) if data else ["status"]
    with (OUT / name).open("w", newline="") as f:
        writer = csv.DictWriter(f, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(data)


def main():
    cfg = json.loads((OUT / "run_config.json").read_text())
    http = RecordedHTTP(OUT / "remote", timeout=30)
    (OUT / "sources").mkdir(exist_ok=True)
    old_contigs = set()
    for folder in ("art_type1_public", "art_type1_mgnify"):
        old_contigs.update(r["accession"] for r in read(ROOT / "results" / folder / "source_records.tsv"))
    old_contigs.update(p.stem for p in (ROOT / "results/art_sample_expansion_20261104/deferred_acquisition").glob("*.fna"))
    hits = [r for r in read(OUT / "remote_hits.tsv") if r["in_original_saved_hit_lists"] == "False"]
    hits.sort(key=lambda r: (int(r["rank"]), r["query_id"], float(r["evalue"])))
    unique = {}
    for hit in hits:
        unique.setdefault(hit["protein_accession"], hit)
    metadata, peptides, locations = [], [], []
    for accession, hit in list(unique.items())[:cfg["max_metadata_proteins"]]:
        try:
            data = http.json("protein_" + accession, MGNIFY + "/protein/" + accession, limit=500000)
            assert data["mgyp"] == accession
            peptide = data["sequence"]
            peptides.append(SeqRecord(Seq(peptide), id=accession, description="incremental_homology_lead"))
            locs = data.get("study_assembly_contigs", [])
            metadata.append(dict(protein_accession=accession, status="public_metadata_retrieved",
                length_aa=len(peptide), full_length=data.get("full_length"), cluster_size=data.get("cluster_size"),
                locations=len(locs), best_evalue=min(float(r["evalue"]) for r in hits if r["protein_accession"] == accession),
                query_ids=";".join(sorted({r["query_id"] for r in hits if r["protein_accession"] == accession})),
                peptide_sha256=hashlib.sha256(peptide.encode()).hexdigest(),
                evidence_level="protein_homology_lead_only"))
            for loc in locs:
                locations.append(dict(protein_accession=accession, peptide=peptide, full_length=data.get("full_length"),
                    known_contig=loc["contig_accession"] in old_contigs, location=loc))
        except Exception as exc:
            metadata.append(dict(protein_accession=accession, status="metadata_failed", error=str(exc)))
    table("protein_metadata.tsv", metadata)
    SeqIO.write(peptides, OUT / "remote_homologs.faa", "fasta")
    dump("ranked_locations.json", [{k:v for k,v in r.items() if k != "peptide"} for r in locations])
    # Full-length metadata is a ranking aid only. DNA translation verifies maps.
    locations.sort(key=lambda r: (r["known_contig"], not bool(r["full_length"]),
        min(int(r["location"]["start"]), int(r["location"]["end"])) < 1800,
        r["protein_accession"], r["location"]["contig_accession"]))
    sources, mappings, attempts, used = [], [], [], 0
    seen = set()
    for item in locations:
        loc = item["location"]
        contig = loc["contig_accession"]
        if item["known_contig"] or contig in seen:
            continue
        if len(seen) >= cfg["max_contigs"] or used >= cfg["max_source_bytes"]:
            break
        seen.add(contig)
        try:
            assembly = loc["assembly_accession"]
            url = "https://www.ebi.ac.uk/metagenomics/api/v2/assemblies/" + assembly + "/analyses?page_size=1"
            data = http.json("assembly_" + assembly, url, limit=500000)
            analyses = data.get("items", data.get("results", []))
            if not analyses:
                raise ValueError("No public per-contig analysis metadata")
            analysis = analyses[0].get("accession") or analyses[0].get("mgya_accession")
            endpoint = "https://www.ebi.ac.uk/metagenomics/api/v1/analyses/" + analysis + "/contigs/" + urllib.parse.quote(loc["contig_name"], safe="")
            raw = http.request("contig_" + contig, endpoint, limit=cfg["max_source_bytes"]-used)
            used += len(raw)
            record = SeqIO.read(io.StringIO(raw.decode()), "fasta")
            if record.id != loc["contig_name"]:
                raise ValueError("Source contig name differs from protein metadata")
            record.id = contig
            record.description = loc["contig_name"]
            record.annotations["molecule_type"] = "DNA"
            record.annotations["topology"] = "linear"
            path = OUT / "sources" / (contig + ".fna")
            SeqIO.write(record, path, "fasta")
            sources.append(dict(accession=contig, path=str(path.relative_to(ROOT)), format="fasta", source_url=endpoint,
                file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                sequence_sha256=hashlib.sha256(str(record.seq).upper().encode()).hexdigest(), bytes=len(raw),
                source_version_kind="release_snapshot", database_release=cfg["expected_release"], assembly=assembly,
                analysis=analysis, source_contig_name=loc["contig_name"], role="diverse_query_public_homology_lead"))
            matches = locate_peptide(record, item["peptide"], {"+":1, "-":-1}.get(loc.get("strand")))
            status = "exact_translated_peptide_mapping" if len(matches) == 1 else "ambiguous_mapping" if matches else "translation_mapping_failed"
            mappings.append(dict(protein_accession=item["protein_accession"], accession=contig, status=status,
                coordinates_json=json.dumps(matches), original_mapping_json=json.dumps(loc),
                evidence_level="peptide_mapped_RT_partner_family_and_array_not_yet_assessed"))
            attempts.append(dict(accession=contig, status=status, length_nt=len(record)))
        except Exception as exc:
            attempts.append(dict(accession=contig, status="retrieval_failed", error=str(exc)))
        table("acquired_source_records.tsv", sources)
        table("homolog_to_contig.tsv", mappings)
        dump("acquisition_attempts.json", attempts)
    table("acquired_source_records.tsv", sources)
    table("homolog_to_contig.tsv", mappings)
    dump("acquisition_attempts.json", attempts)
    dump("acquisition_summary.json", dict(metadata_proteins_attempted=len(metadata), source_contigs=len(sources),
        source_contigs_attempted=len(seen), exact_peptide_maps=sum(r["status"] == "exact_translated_peptide_mapping" for r in mappings),
        source_bytes=used, complete_ART_candidates=0, limits=cfg,
        interpretation="Additional protein leads; no ART RT or complete locus claim before native profiling and source checks"))
    print((OUT / "acquisition_summary.json").read_text())


if __name__ == "__main__":
    main()
