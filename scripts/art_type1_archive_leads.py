"""Retrieve a bounded Logan panel from archived leads, with exact DNA validation.

The external archive supplies leads, not verified ART classifications. Logan
protein IDs were made from older assemblies: never assume a contig number alone
proves a protein-to-DNA mapping. Save fetched compressed prefixes and object ETags.
"""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import tarfile
import urllib.request

from Bio import SeqIO
import zstandard

from art_reconstruction_common import ROOT, Runner, sha, table
from art_type1_public_search import locate_peptide
from extract_art_loci import read_tsv, write_json

ARCHIVE = ROOT/"results/art_discovery/typeIII_source_audit/ai3_core.tar.gz"
OUT = ROOT/"results/art_type1_public"
PANEL = [
    ("SRR26805139",816,"L1_SRR26805139_SRR26805139_816_2_partner.faa"),
    ("ERR10896470",6887,"L1_ERR10896470_ERR10896470_6887_1_partner.faa"),
    ("ERR8055547",17666,"ART_34_ERR8055547_17666_12_partner.faa"),
    ("SRR8925777",19762,"ART_69_SRR8925777_19762_9_partner.faa"),
]


def scan_prefix(raw, target):
    """Return a complete target FASTA record only; a truncated last record fails."""
    dec=zstandard.ZstdDecompressor().decompressobj()
    pending=b"";lines=[];active=False
    for pos in range(0,len(raw),65536):
        pending+=dec.decompress(raw[pos:pos+65536])
        parts=pending.split(b"\n");pending=parts.pop()
        for line in parts:
            if line.startswith(b">"):
                if active:
                    return SeqIO.read(io.StringIO(b"\n".join(lines).decode()),"fasta")
                active=line[1:].split()[0].decode()==target
            if active:lines.append(line)
    if dec.eof and active:
        if pending:lines.append(pending)
        return SeqIO.read(io.StringIO(b"\n".join(lines).decode()),"fasta")
    return None


def complete_records(raw):
    dec=zstandard.ZstdDecompressor().decompressobj()
    text=dec.decompress(raw).decode()
    if not dec.eof:
        text=text[:text.rfind("\n>")+1] if "\n>" in text else ""
    return SeqIO.parse(io.StringIO(text),"fasta")


def download_current_prefixes():
    """Cache four 5 MB current-assembly prefixes, reproducibly from a fresh run."""
    folder=OUT/"archive_leads";folder.mkdir(parents=True,exist_ok=True)
    def download(item):
        run=item[0];path=folder/(run+".prefix.zst");mp=folder/(run+".download.json")
        if path.exists() and mp.exists():
            if sha(path)!=json.loads(mp.read_text())["compressed_prefix_sha256"]:
                raise ValueError("Current prefix checksum mismatch")
            return
        url=f"https://s3.amazonaws.com/logan-pub/c/{run}/{run}.contigs.fa.zst"
        req=urllib.request.Request(url,headers={"Range":"bytes=0-4999999","User-Agent":"ART-locus-research/1.0"})
        with urllib.request.urlopen(req,timeout=60) as response:
            raw=response.read(5000000)
            meta=dict(run=run,url=url,status=response.status,headers=dict(response.headers),
                retrieved_utc=datetime.now(timezone.utc).isoformat(),bytes=len(raw),
                compressed_prefix_sha256=hashlib.sha256(raw).hexdigest())
        path.write_bytes(raw);write_json(mp,meta)
    with ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(download,PANEL))


def reconcile_current():
    """Search current contigs by peptide, never by obsolete numeric contig ID."""
    folder=OUT/"archive_leads";rows=[];audits=[]
    old=read_tsv(OUT/"acquired_source_records.tsv")
    archive=tarfile.open(ARCHIVE)
    for run,number,filename in PANEL:
        oldrow=next((r for r in old if r["accession"]==f"{run}_{number}"),None)
        if oldrow is None:continue
        member=next(m for m in archive.getmembers() if m.name.endswith("/"+filename))
        peptide=SeqIO.read(io.StringIO(archive.extractfile(member).read().decode()),"fasta")
        raw=(folder/(run+".prefix.zst")).read_bytes()
        download=json.loads((folder/(run+".download.json")).read_text())
        matches=[]
        for rec in complete_records(raw):
            if len(rec)<len(peptide)*3:continue
            # Cheap unique peptide probe before full exact mapping.
            probe=str(peptide.seq[10:40]);possible=False
            for seq in (rec.seq,rec.seq.reverse_complement()):
                for frame in range(3):
                    n=(len(seq)-frame)//3*3
                    if probe in str(seq[frame:frame+n].translate(table=11)):
                        possible=True;break
                if possible:break
            if not possible:continue
            maps=locate_peptide(rec,str(peptide.seq))
            if maps:matches.append((rec,maps))
        audit=dict(historical_contig=oldrow["accession"],matches=[dict(contig_id=r.id,maps=m) for r,m in matches],
            current_prefix_bytes=len(raw),current_source_url=download["url"])
        if len(matches)==1 and len(matches[0][1])==1:
            rec,maps=matches[0];rec.annotations.update(molecule_type="DNA",topology="linear")
            dest=folder/(rec.id+".current_verified.gb");SeqIO.write(rec,dest,"genbank")
            row=dict(oldrow,accession=rec.id,path=dest.relative_to(ROOT).as_posix(),
                file_sha256=sha(dest),sequence_sha256=hashlib.sha256(bytes(rec.seq)).hexdigest(),
                source_url=download["url"],database_release="Logan_current_object;"+download["headers"].get("ETag",""),
                length_nt=len(rec),bytes=len(raw),assembly_status="current_contig_exact_peptide_confirmed")
            oldrec=SeqIO.read(ROOT/oldrow["path"],"genbank")
            audit.update(status="current_mapping_confirmed",old_length_nt=len(oldrec),current_length_nt=len(rec),
                         identical_DNA=str(oldrec.seq)==str(rec.seq) or oldrec.seq.reverse_complement()==rec.seq)
            rows.append(row)
        else:
            oldrow["assembly_status"]="historical_v1.0_looping_artifact_risk_current_confirmation_pending"
            rows.append(oldrow);audit["status"]="current_mapping_unresolved_bounded_prefix"
        audits.append(audit)
    archive.close()
    handled={a["historical_contig"] for a in audits}
    rows += [r for r in old if r["accession"] not in handled]
    table(OUT/"acquired_source_records.tsv",list({r["accession"]:r for r in rows}.values()))
    if audits:write_json(folder/"current_assembly_reconciliation.json",audits)
    print(json.dumps(audits,indent=2),flush=True)


def extend_two_prefixes():
    """Additional bounded 20 MB for the two leads outside the first prefixes."""
    tasks=[(OUT/"archive_leads"/version,run) for version in ("","logan_c1.0")
           for run in ("ERR10896470","ERR8055547")]
    def extend(task):
        folder,run=task;path=folder/(run+".prefix.zst");mp=folder/(run+".download.json")
        meta=json.loads(mp.read_text());raw=path.read_bytes()
        if sha(path)!=meta["compressed_prefix_sha256"]:raise ValueError("Checksum mismatch")
        if len(raw)>=10000000:return
        req=urllib.request.Request(meta["url"],headers={"Range":f"bytes={len(raw)}-9999999",
            "If-Match":meta["headers"]["ETag"],"User-Agent":"ART-locus-research/1.0"})
        with urllib.request.urlopen(req,timeout=60) as response:
            if response.status!=206 or not response.headers.get("Content-Range","").startswith(f"bytes {len(raw)}-"):
                raise ValueError("Range not honored")
            extra=response.read(10000000-len(raw))
            meta.setdefault("extensions",[]).append(dict(range=response.headers["Content-Range"],
                bytes=len(extra),utc=datetime.now(timezone.utc).isoformat()))
        path.write_bytes(raw+extra)
        meta.update(bytes=len(raw)+len(extra),compressed_prefix_sha256=sha(path))
        write_json(mp,meta)
    with ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(extend,tasks))


def fetch(item):
    run,number,filename=item
    folder=OUT/"archive_leads"/"logan_c1.0";folder.mkdir(parents=True,exist_ok=True)
    target=f"{run}_{number}"
    url=f"https://s3.amazonaws.com/logan-pub/c1.0/{run}/{run}.contigs.fa.zst"
    path=folder/(run+".prefix.zst")
    meta_path=folder/(run+".download.json")
    # Fixed 5 MB per run, maximum 20 MB for four runs, no bulk database download.
    if path.exists() and meta_path.exists():
        meta=json.loads(meta_path.read_text())
        if sha(path)!=meta["compressed_prefix_sha256"]:raise ValueError("Prefix checksum mismatch")
        raw=path.read_bytes()
    else:
        req=urllib.request.Request(url,headers={"Range":"bytes=0-4999999","User-Agent":"ART-locus-research/1.0"})
        with urllib.request.urlopen(req,timeout=60) as response:
            raw=response.read(5000000)
            meta=dict(run=run,url=url,status=response.status,headers=dict(response.headers),
                      retrieved_utc=datetime.now(timezone.utc).isoformat(),bytes=len(raw),
                      compressed_prefix_sha256=hashlib.sha256(raw).hexdigest())
        path.write_bytes(raw);write_json(meta_path,meta)
    rec=scan_prefix(raw,target)
    if rec is None:
        return None,dict(contig_id=target,status="target_not_in_bounded_prefix",bytes=len(raw))
    with tarfile.open(ARCHIVE) as archive:
        member=next(m for m in archive.getmembers() if m.name.endswith("/"+filename))
        peptide=SeqIO.read(io.StringIO(archive.extractfile(member).read().decode()),"fasta")
    maps=locate_peptide(rec,str(peptide.seq))
    provenance=dict(contig_id=target,archive_member=member.name,archive_sha256=sha(ARCHIVE),
        archive_doi="10.5281/zenodo.23086637",lead_protein_id=peptide.id,
        lead_protein_sha256=hashlib.sha256(bytes(peptide.seq)).hexdigest(),
        coordinates_json=json.dumps(maps),status="exact_translated_mapping" if len(maps)==1 else "unverified_mapping",
        source_url=url,bytes=len(raw))
    # Retain mismatches for diagnosis, but never promote them to the acquired panel.
    rec.annotations.update(molecule_type="DNA",topology="linear")
    rec.description+="; Logan v1.0; reported assembly looping artifacts require current-assembly confirmation"
    dest=folder/(target+".gb");SeqIO.write(rec,dest,"genbank")
    if len(maps)!=1:return None,provenance
    return dict(accession=rec.id,path=dest.relative_to(ROOT).as_posix(),format="genbank",
        file_sha256=sha(dest),sequence_sha256=hashlib.sha256(bytes(rec.seq)).hexdigest(),
        source_url=url,source_version_kind="Logan_object_ETag_and_content_SHA256",
        database_release="Logan_c1.0;"+meta["headers"].get("ETag",meta["headers"].get("Etag","")),
        role="public_archive_lead",length_nt=len(rec),bytes=len(raw)),provenance


def main():
    folder=OUT/"archive_leads";folder.mkdir(parents=True,exist_ok=True)
    # Cross-screen every archived F2 sequence and known GNAT controls; no automatic
    # profile training on these putative homologs.
    proteins=[];controls=[]
    with tarfile.open(ARCHIVE) as archive:
        for member in archive.getmembers():
            if member.name.endswith("tables/d5d/F2.faa"):
                proteins=list(SeqIO.parse(io.StringIO(archive.extractfile(member).read().decode()),"fasta"))
            if Path(member.name).name in {"Eis_P9WFK7.faa","GlmU_P0ACC7.faa","MshD_P9WJM7.faa","RimI_P0A944.faa","tGCN5_1QSR_A.faa"}:
                controls+=list(SeqIO.parse(io.StringIO(archive.extractfile(member).read().decode()),"fasta"))
    SeqIO.write(proteins+controls,folder/"leads_and_GNAT_controls.faa","fasta")
    run=Runner(folder);hits=[]
    for family,hmm in [
        ("Type_I_partner",ROOT/"results/art_discovery/profiles/Type_I_partner.hmm"),
        ("TypeI_two_seed",OUT/"profiles/TypeI_two_seed_exploratory.hmm"),
        ("ART_RT",ROOT/"results/art_discovery/profiles/ART_RT.hmm")]:
        hs,_=run.hmm(hmm,folder/"leads_and_GNAT_controls.faa",folder/family,z=len(proteins)+len(controls))
        hits.extend(dict(h,family=family) for h in hs)
    table(folder/"profile_cross_screen.tsv",hits)
    table(folder/"control_labels.tsv",[dict(protein_id=r.id,label="unrelated_GNAT_from_external_archive") for r in controls])
    with ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(fetch,PANEL))
    rows=[r for r,p in results if r];provenance=[p for r,p in results]
    table(folder/"logan_c1.0"/"lead_mapping.tsv",provenance)
    previous=read_tsv(OUT/"acquired_source_records.tsv") if (OUT/"acquired_source_records.tsv").exists() else []
    merged={r["accession"]:r for r in previous+rows}
    table(OUT/"acquired_source_records.tsv",list(merged.values()))
    write_json(folder/"logan_c1.0"/"summary.json",dict(archive_doi="10.5281/zenodo.23086637",
        archive_sha256=sha(ARCHIVE),archived_protein_leads=len(proteins),GNAT_controls=len(controls),
        panel_requested=len(PANEL),exact_mapped_contigs=len(rows),compressed_bytes=sum(p["bytes"] for p in provenance),
        status="external_leads_require_independent_family_and_locus_assessment",
        mappings=provenance))
    print((folder/"logan_c1.0"/"summary.json").read_text(),flush=True)


if __name__=="__main__":
    import sys
    if "--reconcile-current" in sys.argv:reconcile_current()
    elif "--extend-two-prefixes" in sys.argv:extend_two_prefixes()
    elif "--download-current" in sys.argv:download_current_prefixes()
    else:main()
