"""Recorded, resumable public homology acquisition; never invent nucleotide maps."""
import argparse
from datetime import datetime, timezone
import hashlib
import io
import json
import math
from pathlib import Path
import re
import time
import urllib.error
import urllib.parse
import urllib.request

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from art_reconstruction_common import ROOT, sha, table
from extract_art_loci import write_json

HMMER = "https://www.ebi.ac.uk/Tools/hmmer/api/v1"
BLAST = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
MGNIFY = "https://www.ebi.ac.uk/metagenomics/proteins/api/v1"


class RecordedHTTP:
    def __init__(self, out, timeout=45):
        self.out = Path(out)
        self.out.mkdir(parents=True, exist_ok=True)
        self.timeout = timeout

    def request(self, name, url, payload=None, form=False, cached=True, limit=20000000):
        path = self.out / (name + ".response")
        meta = self.out / (name + ".request.json")
        request = dict(url=url, payload=payload, form=form, method="POST" if payload is not None else "GET")
        if cached and path.exists() and meta.exists():
            recorded = json.loads(meta.read_text())
            if recorded["request"] != request:
                raise ValueError("Cache request differs: " + name)
            if recorded["status"] == 200 and not recorded.get("error"):
                assert sha(path) == recorded["response_sha256"]
                if path.stat().st_size > limit:
                    raise ValueError("Cached response exceeds byte budget")
                return path.read_bytes()
        data = None
        headers = {"User-Agent": "ART-locus-research/1.0", "Accept": "application/json"}
        if payload is not None:
            data = (urllib.parse.urlencode(payload) if form else json.dumps(payload)).encode()
            headers["Content-Type"] = "application/x-www-form-urlencoded" if form else "application/json"
        started = time.perf_counter()
        status = None
        body = b""
        error = None
        try:
            with urllib.request.urlopen(urllib.request.Request(url, data=data, headers=headers),
                                        timeout=self.timeout) as response:
                status = response.status
                body = response.read(limit + 1)
                if len(body) > limit:
                    raise ValueError("Response exceeds byte budget")
        except urllib.error.HTTPError as exc:
            status = exc.code
            body = exc.read(200000)
            error = str(exc)
        except Exception as exc:
            error = type(exc).__name__ + ": " + str(exc)
        path.write_bytes(body)
        record = dict(request=request, status=status, error=error,
                      response_sha256=sha(path), bytes=len(body),
                      time_utc=datetime.now(timezone.utc).isoformat(),
                      elapsed_seconds=round(time.perf_counter() - started, 4))
        write_json(meta, record)
        with (self.out / "requests.jsonl").open("a", encoding="utf8") as stream:
            stream.write(json.dumps(dict(name=name, **record)) + "\n")
        if error:
            raise RuntimeError(name + ": " + error)
        return body

    def json(self, *args, **kwargs):
        return json.loads(self.request(*args, **kwargs))


def queries(cfg):
    rows = []
    for family, key in (("Type_I_partner", "query_partner_fasta"), ("ART_RT", "query_RT_fasta")):
        for record in SeqIO.parse(ROOT / cfg[key], "fasta"):
            rows.append(dict(query_id=record.id, family=family, sequence=str(record.seq),
                             sequence_sha256=hashlib.sha256(str(record.seq).encode()).hexdigest()))
    return rows


def submit(cfg, out):
    http = RecordedHTTP(out / "remote")
    qs = queries(cfg)
    write_json(out / "queries.json", qs)
    write_json(out / "run_config.json", cfg)
    databases = http.json("hmmer_databases", HMMER + "/search/databases")
    db = next(x for x in databases if x["id"] == cfg["hmmer_database"] and x["status"] == "enabled")
    if db["version"] != cfg["hmmer_expected_release"]:
        raise ValueError("MGnify release changed; freeze a new explicit configuration")
    write_json(out / "hmmer_database.json", db)
    statuses = []
    for q in qs:
        try:
            payload = dict(input=f">{q['query_id']}\n{q['sequence']}\n", database=db["id"],
                           E=cfg["remote_report_evalue"], domE=cfg["remote_report_evalue"],
                           incE=1e-5, incdomE=1e-5)
            job = http.json("hmmer_submit_" + q["query_id"], HMMER + "/search/phmmer", payload)
            statuses.append(dict(service="HMMER", query_id=q["query_id"], family=q["family"],
                                 job_id=job["id"], status="submitted"))
        except Exception as exc:
            statuses.append(dict(service="HMMER", query_id=q["query_id"], family=q["family"],
                                 status="submission_failed", error=str(exc)))
    # One batch job, four public protein queries. No email address is supplied.
    fasta = "".join(f">{q['query_id']}\n{q['sequence']}\n" for q in qs)
    payload = dict(CMD="Put", PROGRAM="blastp", DATABASE=cfg["blast_database"], QUERY=fasta,
                   EXPECT=cfg["remote_report_evalue"], HITLIST_SIZE=cfg["maximum_hits_per_query"],
                   FILTER="L", COMPOSITION_BASED_STATISTICS=2)
    try:
        body = http.request("blast_submit", BLAST, payload, form=True).decode()
        match = re.search(r"RID\s*=\s*(\S+)", body)
        if not match:
            raise ValueError("NCBI did not return a BLAST job identifier")
        statuses.append(dict(service="NCBI_BLAST", query_id="four_query_batch",
                             job_id=match.group(1), status="submitted"))
    except Exception as exc:
        statuses.append(dict(service="NCBI_BLAST", query_id="four_query_batch",
                             status="submission_failed", error=str(exc)))
    write_json(out / "jobs.json", statuses)
    print(json.dumps(statuses, indent=2))


def poll(cfg, out):
    http = RecordedHTTP(out / "remote")
    jobs = json.loads((out / "jobs.json").read_text())
    for job in jobs:
        if "job_id" not in job:
            continue
        if job["service"] == "HMMER":
            url = HMMER + "/result/" + job["job_id"] + "?" + urllib.parse.urlencode(
                dict(page=1, page_size=cfg["maximum_hits_per_query"], with_domains="true"))
            name = "hmmer_result_" + job["query_id"]
            try:
                data = http.json(name, url, cached=job.get("status") == "SUCCESS")
                job["status"] = data.get("status", "unknown")
                job["result_path"] = (out / "remote" / (name + ".response")).relative_to(ROOT).as_posix()
            except Exception as exc:
                job["last_error"] = str(exc)
        else:
            url = BLAST + "?" + urllib.parse.urlencode(
                dict(CMD="Get", RID=job["job_id"], FORMAT_TYPE="JSON2_S"))
            try:
                raw = http.request("blast_result", url, cached=job.get("status") == "SUCCESS")
                if raw.lstrip().startswith(b"{"):
                    data = json.loads(raw)
                    if "BlastOutput2" in data:
                        job["status"] = "SUCCESS"
                        job["result_path"] = (out / "remote/blast_result.response").relative_to(ROOT).as_posix()
                    else:
                        job["status"] = "unrecognized_JSON"
                else:
                    match = re.search(r"Status=(\w+)", raw.decode(errors="replace"))
                    job["status"] = match.group(1) if match else "UNRECOGNIZED_RESPONSE"
            except Exception as exc:
                job["last_error"] = str(exc)
    write_json(out / "jobs.json", jobs)
    print(json.dumps(jobs, indent=2))


def normalize_blast(data, qs):
    """Keep all accession aliases for an nr hit; retain query/target coverage."""
    rows = []
    entries = data["BlastOutput2"]
    if isinstance(entries, dict):
        entries = [entries]
    for entry in entries:
        search = entry["report"]["results"]["search"]
        qid = search["query_title"].split()[0]
        q = next(q for q in qs if q["query_id"] == qid)
        for rank, hit in enumerate(search.get("hits", []), 1):
            if not hit.get("hsps"):
                continue
            hsp = min(hit["hsps"], key=lambda h: h["evalue"])
            for description in hit.get("description", [])[:10]:
                rows.append(dict(service="NCBI_BLAST", query_id=qid, family=q["family"],
                    rank=rank, protein_accession=description["accession"],
                    description=description.get("title", ""), bits=hsp["bit_score"],
                    evalue=hsp["evalue"], query_coverage=(hsp["query_to"]-hsp["query_from"]+1)/search["query_len"],
                    target_coverage=(abs(hsp["hit_to"]-hsp["hit_from"])+1)/hit["len"],
                    target_length=hit["len"], query_from1=hsp["query_from"],
                    query_to1=hsp["query_to"], raw_hit_aliases=len(hit.get("description", []))))
    return rows


def normalize_hmmer(data, job):
    rows = []
    for rank, hit in enumerate((data.get("result") or {}).get("hits", []), 1):
        matches = re.findall(r"MGYP\d+", " ".join(str(x) for x in (
            hit.get("name",""),hit.get("acc",""),(hit.get("metadata") or {}).get("accession",""))))
        if not matches:
            raise ValueError("MGnify hit lacks an MGYP accession")
        ev = hit.get("evalue", hit.get("e_value"))
        rows.append(dict(service="HMMER", query_id=job["query_id"], family=job["family"],
                         rank=rank, protein_accession=matches[0], bits=hit.get("score", ""),
                         evalue=ev if ev is not None else "",
                         description=hit.get("desc", hit.get("description", "")),
                         reported_by_server=True))
    return rows


def collect(cfg, out):
    qs = queries(cfg)
    jobs = json.loads((out / "jobs.json").read_text())
    rows, status = [], []
    for job in jobs:
        if job.get("status") != "SUCCESS":
            status.append(dict(service=job["service"], query_id=job["query_id"],
                               status=job.get("status"), complete=False))
            continue
        data = json.loads((ROOT / job["result_path"]).read_text())
        found = normalize_blast(data, qs) if job["service"] == "NCBI_BLAST" else normalize_hmmer(data, job)
        rows += found
        status.append(dict(service=job["service"], query_id=job["query_id"], complete=True,
                           hits_returned=len(found), cap=cfg["maximum_hits_per_query"],
                           potentially_capped=bool(data.get("page_count", 1) > 1) if
                           job["service"] == "HMMER" else any(r["rank"] >= cfg["maximum_hits_per_query"] for r in found)))
    table(out / "remote_hits.tsv", rows, fields=None if rows else ["service", "query_id", "protein_accession"])
    write_json(out / "retrieval_status.json", status)
    print(json.dumps(dict(reported_hit_rows=len(rows), search_status=status), indent=2))


def accession_from_contig_name(name):
    # ENA's processed-contig header embeds an actual accession, not an MGYC.
    match = re.search(r"(?:^|[-|])([A-Z]{1,6}_?\d{5,12}\.\d+)(?:[-| ]|$)", name)
    if not match:
        return None
    accession=match.group(1)
    # ENA analysis/assembly IDs followed by a contig number are not sequence
    # accession versions (e.g. ERZ2641806.86530).
    if re.match(r"^(?:ERZ|SRZ|DRZ|MGY|ERR|SRR|DRR|ERS|SRS|DRS)",accession):
        return None
    return accession


def fetch_nucleotide(http, accession, out, remaining):
    if not re.fullmatch(r"[A-Z]{1,6}_?\d{5,12}\.\d+", accession):
        raise ValueError("Expected a nucleotide accession.version: " + accession)
    errors = []
    for fmt, url in (
        ("embl", "https://www.ebi.ac.uk/ena/browser/api/embl/" + accession + "?download=false"),
        ("genbank", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(
            dict(db="nuccore", id=accession, rettype="gbwithparts", retmode="text")))):
        try:
            raw = http.request("nt_" + accession + "_" + fmt, url, limit=remaining)
            record = SeqIO.read(io.StringIO(raw.decode()), fmt)
            if record.id != accession:
                raise ValueError("Returned nucleotide version differs")
            path = out / "sources" / (accession + (".embl" if fmt == "embl" else ".gb"))
            path.write_bytes(raw)
            return record, dict(accession=accession, path=path.relative_to(ROOT).as_posix(),
                format=fmt, source_url=url, file_sha256=sha(path),
                sequence_sha256=hashlib.sha256(str(record.seq).upper().encode()).hexdigest(),
                bytes=len(raw), source_version_kind="accession_version", database_release="")
        except Exception as exc:
            errors.append(str(exc))
    raise ValueError("; ".join(errors))


def mgnify_contig(http, location, cfg, out, remaining):
    contig_name = location["contig_name"]
    acc = accession_from_contig_name(contig_name)
    if acc:
        return fetch_nucleotide(http, acc, out, remaining)
    assembly = location["assembly_accession"]
    url = "https://www.ebi.ac.uk/metagenomics/api/v2/assemblies/" + assembly + "/analyses?page_size=100"
    data = http.json("assembly_" + assembly, url)
    analyses = data.get("items", data.get("results", []))
    errors = []
    for analysis in analyses:
        aid = analysis.get("accession") or analysis.get("mgya_accession")
        if not aid:
            continue
        # v1 is the documented historical per-contig endpoint; v2 lacks it.
        endpoint = "https://www.ebi.ac.uk/metagenomics/api/v1/analyses/" + aid + "/contigs/" + urllib.parse.quote(contig_name, safe="")
        try:
            raw = http.request("contig_" + location["contig_accession"] + "_" + aid, endpoint, limit=remaining)
            rec = SeqIO.read(io.StringIO(raw.decode()), "fasta")
            if rec.id != contig_name:
                raise ValueError("Returned contig name differs")
            rec.id = location["contig_accession"]
            rec.description = contig_name
            rec.annotations["molecule_type"] = "DNA"
            rec.annotations["topology"] = "linear"
            path = out / "sources" / (rec.id + ".fna")
            SeqIO.write(rec, path, "fasta")
            return rec, dict(accession=rec.id, path=path.relative_to(ROOT).as_posix(),
                format="fasta", source_url=endpoint, file_sha256=sha(path),
                sequence_sha256=hashlib.sha256(str(rec.seq).upper().encode()).hexdigest(),
                bytes=len(raw), source_version_kind="release_snapshot",
                database_release=cfg["hmmer_expected_release"], assembly=assembly, analysis=aid,
                source_contig_name=contig_name)
        except Exception as exc:
            errors.append(str(exc))
    raise ValueError("No bounded per-contig retrieval; assembly bulk download deferred: " + "; ".join(errors))


def locate_peptide(record, peptide, strand=None):
    """Require an exact genomic translation; do not guess from header coordinates.

    Record alternative initiator handling separately. Stop-free fragments may be
    partial; completeness is determined by the downstream annotation workflow.
    """
    matches = []
    peptide = peptide.rstrip("*")
    for code in (11, 15):
        for orientation in (1, -1):
            if strand and orientation != strand:
                continue
            seq = record.seq if orientation == 1 else record.seq.reverse_complement()
            for frame in range(3):
                usable = (len(seq)-frame)//3*3
                protein = str(seq[frame:frame+usable].translate(table=code))
                # Search through all occurrences; don't silently choose a paralog.
                probe = peptide[1:] if peptide.startswith("M") else peptide
                offset = 1 if peptide.startswith("M") else 0
                for match in re.finditer("(?=" + re.escape(probe) + ")", protein):
                    aa = match.start() - offset
                    if aa < 0:
                        continue
                    a, b = frame + 3*aa, frame + 3*(aa+len(peptide))
                    coding = seq[a:b]
                    actual = str(coding.translate(table=code))
                    start_ok = actual == peptide or (peptide.startswith("M") and
                        str(coding[:3]) in {"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"} and actual[1:] == peptide[1:])
                    if not start_ok:
                        continue
                    x,y = (a,b) if orientation == 1 else (len(record)-b,len(record)-a)
                    matches.append(dict(start0=x,end0=y,strand=orientation,genetic_code=code,
                                        initiator_normalized=actual != peptide))
    # Same coordinate under both codes is a single mapping; prefer standard code.
    unique = {}
    for m in matches:
        unique.setdefault((m["start0"],m["end0"],m["strand"]),m)
    return list(unique.values())


def acquire(cfg, out):
    """Bounded, route-balanced acquisition, retaining failures and missing maps."""
    from extract_art_loci import read_tsv
    http = RecordedHTTP(out / "remote")
    (out / "sources").mkdir(exist_ok=True)
    hits = read_tsv(out / "remote_hits.tsv")
    # Round-robin rank across both seeds/routes and databases before deduplication.
    hits.sort(key=lambda r:(int(r["rank"]), r["query_id"], r["service"], r["protein_accession"]))
    unique = {}
    for h in hits:
        unique.setdefault((h["service"],h["protein_accession"]),h)
    sources, sequences, mappings, failures = {}, {}, [], []
    if (out/"remote_homologs.faa").exists():
        sequences={r.id:r for r in SeqIO.parse(out/"remote_homologs.faa","fasta")}
    if (out/"homolog_to_contig.tsv").exists():
        mappings=read_tsv(out/"homolog_to_contig.tsv")
    if (out/"unresolved_or_deferred.tsv").exists():
        failures=read_tsv(out/"unresolved_or_deferred.tsv")
    if (out/"acquired_source_records.tsv").exists():
        for meta in read_tsv(out/"acquired_source_records.tsv"):
            if sha(ROOT/meta["path"])!=meta["file_sha256"]:raise ValueError("Existing source checksum mismatch")
            sources[meta["accession"]]=(SeqIO.read(ROOT/meta["path"],meta["format"]),meta)
    used_bytes = 0
    for index, ((service, accession), hit) in enumerate(unique.items()):
        if index >= cfg.get("maximum_acquisition_proteins",40):
            failures.append(dict(protein_accession=accession,service=service,status="deferred_protein_attempt_budget"))
            continue
        if any(m.get("protein_accession")==accession and m.get("service")==service and
               m.get("status")=="exact_translated_mapping" for m in mappings):
            continue
        if len(sources) >= cfg["maximum_contigs"]:
            failures.append(dict(protein_accession=accession,status="deferred_contig_budget"))
            continue
        try:
            if service == "HMMER":
                data = http.json("protein_" + accession, MGNIFY + "/protein/" + accession)
                peptide = data["sequence"]
                locations = data.get("study_assembly_contigs", [])
            else:
                url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(
                    dict(db="protein",id=accession,rettype="gp",retmode="text"))
                raw = http.request("protein_" + accession, url)
                protein = SeqIO.read(io.StringIO(raw.decode()), "genbank")
                peptide = str(protein.seq)
                coded = [x for f in protein.features for x in f.qualifiers.get("coded_by", [])]
                ntids = sorted(set(re.findall(r"([A-Z]{1,6}_?\d{5,12}\.\d+):", ";".join(coded))))
                locations = [dict(nucleotide_accession=a, coded_by=";".join(coded)) for a in ntids]
            sequences[accession] = SeqRecord(Seq(peptide), id=accession, description=service)
            if not locations:
                raise ValueError("Protein-only record has no nucleotide mapping")
            for loc in locations[:3]:
                if len(sources) >= cfg["maximum_contigs"]:
                    failures.append(dict(protein_accession=accession,status="deferred_contig_budget"))
                    break
                remaining = cfg["maximum_source_download_bytes"] - used_bytes
                if remaining <= 0:
                    raise ValueError("Source download budget exhausted")
                existing = loc.get("nucleotide_accession") or accession_from_contig_name(loc.get("contig_name","")) or loc.get("contig_accession")
                if existing in sources:
                    record, meta = sources[existing]
                elif service == "HMMER":
                    record, meta = mgnify_contig(http,loc,cfg,out,remaining)
                else:
                    record, meta = fetch_nucleotide(http,loc["nucleotide_accession"],out,remaining)
                if record.id not in sources:
                    meta.setdefault("role","public_remote_candidate")
                    sources[record.id] = (record,meta)
                    used_bytes += meta["bytes"]
                strand = {"-":-1, "+":1}.get(loc.get("strand"))
                matches = locate_peptide(record,peptide,strand)
                status = "exact_translated_mapping" if len(matches) == 1 else "ambiguous_multiple_mappings" if matches else "translation_mapping_failed"
                mappings.append(dict(protein_accession=accession, service=service,
                    accession=record.id, status=status, coordinates_json=json.dumps(matches),
                    original_mapping_json=json.dumps(loc), peptide_sha256=hashlib.sha256(peptide.encode()).hexdigest(),
                    discovery_queries=";".join(sorted({r["query_id"] for r in hits if r["protein_accession"]==accession})),
                    discovery_families=";".join(sorted({r["family"] for r in hits if r["protein_accession"]==accession}))))
                if status == "exact_translated_mapping":
                    break
        except Exception as exc:
            failures.append(dict(protein_accession=accession,service=service,status="unresolved",reason=str(exc)))
        # Flush after each record so failed/interrupted runs retain progress.
        table(out / "acquired_source_records.tsv", [m for r,m in sources.values()])
        table(out / "homolog_to_contig.tsv", mappings)
        table(out / "unresolved_or_deferred.tsv", failures)
        SeqIO.write(list(sequences.values()),out / "remote_homologs.faa","fasta")
        print(json.dumps(dict(processed=index+1,unique_proteins=len(unique),contigs=len(sources),
                              mapped=sum(m["status"]=="exact_translated_mapping" for m in mappings),
                              source_bytes=used_bytes)),flush=True)
    table(out / "acquired_source_records.tsv", [m for r,m in sources.values()])
    table(out / "homolog_to_contig.tsv", mappings)
    table(out / "unresolved_or_deferred.tsv", failures)
    SeqIO.write(list(sequences.values()),out / "remote_homologs.faa","fasta")
    write_json(out / "acquisition_summary.json", dict(unique_reported_proteins=len(unique),
        downloaded_protein_sequences=len(sequences), source_contigs=len(sources),
        exact_protein_maps=sum(m["status"]=="exact_translated_mapping" for m in mappings),
        unresolved_or_deferred=len(failures), source_bytes=used_bytes,
        source_budget=cfg["maximum_source_download_bytes"], selection="round_robin_query_rank_bounded",
        not_a_prevalence_estimate=True))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("operation", choices=("submit", "poll", "collect", "acquire"))
    ap.add_argument("--config", type=Path, default=ROOT / "config/art_type1_public.json")
    args = ap.parse_args()
    cfg = json.loads(args.config.read_text())
    out = ROOT / cfg["output"]
    out.mkdir(parents=True, exist_ok=True)
    {"submit": submit, "poll": poll, "collect": collect, "acquire": acquire}[args.operation](cfg, out)


if __name__ == "__main__":
    main()
