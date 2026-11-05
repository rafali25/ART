"""Bounded type I locus processing, separate from the frozen replication outputs.

All DNA intervals are zero based and half open. Public sequence/coordinate
evidence is preserved. Array and family evidence never silently become a
validated biological classification.
"""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import copy
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import shutil
import time

from Bio import AlignIO, Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation, BeforePosition, AfterPosition
import pyrodigal_gv

from art_reconstruction_common import ROOT, Runner, sha, table, passes, validation_gate
from extract_art_loci import (read_tsv, write_json, pid, plan_locus, extract_record,
                              source_pieces, write_gff)
from art_discovery import neighbors
from audit_art_type1_fasta import six_frame_orfs, matching_annotations
from art_array_detector import (refinement_scan, initial_scan, empirical, stable_seed,
                                BASES, null_scores)
from art_array_assessment import (aligned_pwm, quality_flags, short_period_scan,
                                  sensitivity_batch)


def profile_setup(out, cfg, run):
    folder = out / "profiles"
    folder.mkdir(exist_ok=True)
    source = ROOT / cfg["query_partner_fasta"]
    seeds = folder / "typeI_two_LPJP1_queries.faa"
    shutil.copy2(source, seeds)
    aligned = folder / "typeI_two_LPJP1_queries.aligned.faa"
    run.run(["cmd.exe","/d","/c",ROOT/".tools/mafft-win/mafft.bat","--localpair",
             "--maxiterate","1000","--thread","1","--out",aligned,seeds],
            folder/"MAFFT.log",expected=[aligned])
    hmm = folder / "TypeI_two_seed_exploratory.hmm"
    run.run([ROOT/".tools/hmmer-install/bin/hmmbuild.exe","--amino",hmm,aligned],
            folder/"hmmbuild.log",expected=[hmm])
    write_json(folder/"calibration_status.json",dict(
        seed_proteins=[r.id for r in SeqIO.parse(seeds,"fasta")],source_genome_groups=1,
        positive_independent_holdout_groups=0,status="exploratory_retrieval_only",
        production_single_seed_profile_unchanged=True,
        limitation="Two paralogs from LPJP1 cannot validate broad type I sensitivity or specificity"))
    return hmm


def merged_gene_calls(record, out):
    """Preserve source genes; add predictions only when no same-frame CDS exists."""
    rec = copy.deepcopy(record)
    rec.annotations.setdefault("molecule_type","DNA")
    rec.annotations.setdefault("topology","linear")
    finder = pyrodigal_gv.ViralGeneFinder(meta=True)
    predicted = finder.find_genes(bytes(record.seq))
    alternatives = []
    for i, g in enumerate(predicted,1):
        row = dict(start0=g.begin-1,end0=g.end,strand=g.strand)
        source_matches = matching_annotations(row,record.features)
        alternatives.append(dict(contig_id=rec.id,prediction_id=f"{rec.id}_pgv{i}",
            **row,genetic_code=g.translation_table,source_CDS_overlap=source_matches,
            action="alternative_of_source_gene" if source_matches else "added_candidate_CDS",
            translation=g.translate().rstrip("*")))
        if source_matches:
            continue
        a = BeforePosition(g.begin-1) if g.partial_begin else g.begin-1
        b = AfterPosition(g.end) if g.partial_end else g.end
        rec.features.append(SeqFeature(SimpleLocation(a,b,strand=g.strand),type="CDS",qualifiers={
            "protein_id":[f"{rec.id}_pgv{i}"],"locus_tag":[f"{rec.id}_pgv{i}"],
            "translation":[g.translate().rstrip("*")],"transl_table":[str(g.translation_table)],
            "inference":["ab initio prediction:Pyrodigal-gv:"+pyrodigal_gv.__version__],
            "note":["Candidate gene call; original deposited features retained"]}))
    SeqIO.write(rec,out/"annotated"/(rec.id+".gb"),"genbank")
    return rec, alternatives, dict(contig_id=rec.id,model=str(predicted.metagenomic_bin),
        codes=sorted({g.translation_table for g in predicted}),predictions=len(predicted),
        original_CDS=sum(f.type=="CDS" for f in record.features))


def ingest_records(records, metadata, out):
    genes, features, proteins = {}, {}, []
    by_contig = defaultdict(list)
    for rec in records:
        for i,f in enumerate(rec.features):
            if f.type != "CDS" or not f.qualifiers.get("translation"):
                continue
            text = f.qualifiers["translation"][0].rstrip("*")
            if not text or "*" in text:
                continue
            key = "t"+hashlib.sha256((rec.id+"|"+str(f.location)+"|"+pid(f)).encode()).hexdigest()[:18]
            code = int(f.qualifiers.get("transl_table",["11"])[0])
            try:
                computed = str(f.extract(rec.seq)[int(f.qualifiers.get("codon_start",["1"])[0])-1:].translate(table=code,cds=True))
                qc = "complete_translation_matches" if computed == text else "translation_disagrees"
            except Exception:
                qc = "complete_translation_not_confirmed"
            g = dict(protein_key=key,gene_key=key,accession=rec.id,protein_id=pid(f) or f"CDS{i}",
                     start0=int(f.location.start),end0=int(f.location.end),strand=f.location.strand,
                     length_aa=len(text),genetic_code=code,translation_qc=qc,
                     source_parts_json=json.dumps([[int(p.start),int(p.end),p.strand] for p in f.location.parts]),
                     annotation_method="pyrodigal_gv" if "_pgv" in pid(f) else "source_original",
                     role=metadata[rec.id]["role"])
            genes[key]=g
            features[key]=f
            by_contig[rec.id].append(g)
            proteins.append(SeqRecord(Seq(text),id=key,description=rec.id+" "+g["protein_id"]))
    return genes,features,by_contig,proteins


def best_hits(hits, cfg):
    best = defaultdict(dict)
    for h in hits:
        family = h["family"]
        if family == "RVT_1":
            eligible = h["domain_evalue"] <= .001
        else:
            eligible = passes(h,cfg)
        if eligible:
            previous = best[family].get(h["target"])
            better = (h["protein_envelope_start0"] < previous["protein_envelope_start0"]
                      if previous is not None and family == "RVT_1"
                      else previous is None or h["domain_bits"] > previous["domain_bits"])
            if previous is None or better:
                best[family][h["target"]]=h
    return best


def process(cfg,out,include_control=True):
    started=time.perf_counter()
    for name in ("annotated","searches","loci","alignments","array_details","profiles"):
        (out/name).mkdir(parents=True,exist_ok=True)
    run=Runner(out)
    source_rows=read_tsv(out/"acquired_source_records.tsv") if (out/"acquired_source_records.tsv").exists() else []
    if include_control:
        row=next(r for r in read_tsv(ROOT/"results/art_discovery/source_records.tsv") if r["accession"]=="MZ422438.1")
        source_rows.append(dict(row,role="LPJP1_control",source_url="https://www.ebi.ac.uk/ena/browser/view/MZ422438.1",
                               source_version_kind="accession_version",database_release=""))
    records,metadata,alternatives,models,aliases=[],{},[],[],[]
    seen={}
    for source in source_rows:
        if sha(ROOT/source["path"])!=source["file_sha256"]:
            raise ValueError("Source checksum mismatch")
        original=SeqIO.read(ROOT/source["path"],source["format"])
        digest=hashlib.sha256(str(original.seq).upper().encode()).hexdigest()
        if digest in seen:
            aliases.append(dict(accession=original.id,canonical=seen[digest],status="exact_DNA_duplicate"))
            continue
        seen[digest]=original.id
        source.setdefault("role","public_candidate")
        metadata[original.id]=source
        rec,alt,model=merged_gene_calls(original,out)
        records.append(rec);alternatives+=alt;models.append(model)
    if not records:
        raise ValueError("No retrievable nucleotide records; acquisition remains incomplete")
    table(out/"source_records.tsv",list(metadata.values()))
    table(out/"source_aliases.tsv",aliases,fields=["accession","canonical","status"])
    table(out/"alternative_gene_calls.tsv",alternatives)
    write_json(out/"gene_call_models.json",models)
    genes,features,by_contig,proteins=ingest_records(records,metadata,out)
    record_map={r.id:r for r in records}
    # Search six-frame fragments as retrieval representations; they do not
    # automatically enter the ranked gene neighborhood.
    orfs={}
    for rec in records:
        codes=next(m["codes"] for m in models if m["contig_id"]==rec.id)
        for code in sorted(set([11]+codes)):
            for i,o in enumerate(six_frame_orfs(rec.seq,code,cfg["minimum_six_frame_aa"]),1):
                text=o.pop("peptide")
                key="o"+hashlib.sha256(f"{rec.id}|{code}|{i}|{o['start0']}".encode()).hexdigest()[:18]
                orfs[key]=dict(o,protein_key=key,gene_key=key,accession=rec.id,
                    protein_id=key,length_aa=len(text),annotation_method="six_frame_fragment",
                    translation_qc="CDS_boundaries_unconfirmed",role=metadata[rec.id]["role"],
                    source_CDS_overlap=matching_annotations(o,rec.features))
                proteins.append(SeqRecord(Seq(text),id=key,description=rec.id+" six_frame"))
    SeqIO.write(proteins,out/"proteins.faa","fasta")
    protein_map={r.id:r for r in proteins}
    newprofile=profile_setup(out,cfg,run)
    hits=[]
    for family,hmm in (
        ("ART_RT",ROOT/"results/art_discovery/profiles/ART_RT.hmm"),
        ("Type_I_partner",ROOT/"results/art_discovery/profiles/Type_I_partner.hmm"),
        ("TypeI_two_seed",newprofile),
        ("Type_II_partner",ROOT/"results/art_discovery/profiles/Type_II_partner.hmm"),
        ("RVT_1",run.binaries["RVT_1"])):
        hs,_=run.hmm(hmm,out/"proteins.faa",out/"searches"/family,z=len(proteins))
        hits += [dict(h,family=family) for h in hs]
    table(out/"profile_hits.tsv",hits)
    best=best_hits(hits,cfg)
    rescues=[]
    for key,o in orfs.items():
        families=[f for f in ("ART_RT","Type_I_partner","TypeI_two_seed") if key in best[f]]
        if not families:
            continue
        if o["source_CDS_overlap"]:
            rescues.append(dict(**o,matching_families=";".join(families),status="alternative_of_mapped_CDS"))
            # Copy genuine domain evidence to the canonical annotation only if
            # the canonical sequence was itself a passing match; no borrowed score.
            continue
        genes[key]=o;by_contig[o["accession"]].append(o)
        f=SeqFeature(SimpleLocation(BeforePosition(o["start0"]),AfterPosition(o["end0"]),strand=o["strand"]),
            type="CDS",qualifiers={"protein_id":[key],"locus_tag":[key],
                "translation":[str(protein_map[key].seq)],"transl_table":[str(o["genetic_code"])],
                "note":["Six-frame homology rescue; gene start and completeness unconfirmed"]})
        features[key]=f;record_map[o["accession"]].features.append(f)
        rescues.append(dict(**o,matching_families=";".join(families),status="new_partial_candidate_requires_gene_review"))
    table(out/"six_frame_rescue.tsv",rescues)
    for rec in records:
        SeqIO.write(rec,out/"annotated"/(rec.id+".gb"),"genbank")
    table(out/"protein_to_contig.tsv",list(genes.values()))
    table(out/"six_frame_map.tsv",list(orfs.values()))
    # Union independently retrieved family profiles; preserve which model hit.
    typeI={k for fam in ("Type_I_partner","TypeI_two_seed") for k in best[fam] if k in genes}
    rts={k for k in best["ART_RT"] if k in genes}
    routes=defaultdict(set)
    for k in rts:routes[k].add("ART_protein_first")
    neighbor_rows=[];orphans=[];partner_windows=[]
    for p in sorted(typeI):
        pg=genes[p];rec=record_map[pg["accession"]]
        # Partner +/-6kb window is exported even if no upstream RT is found.
        a,b=max(0,pg["start0"]-6000),min(len(rec),pg["end0"]+6000)
        window=rec[a:b]
        window.id=p
        SeqIO.write(window,out/"loci"/(p+".partner_window.gb"),"genbank")
        partner_windows.append(dict(partner_key=p,accession=rec.id,start0=a,end0=b,
            orientation=1,partner_strand=pg["strand"],left_clipped=a==0,right_clipped=b==len(rec),
            path=(out/"loci"/(p+".partner_window.gb")).relative_to(ROOT).as_posix()))
        found=False
        for rank,g in enumerate(neighbors(pg,by_contig[pg["accession"]],"upstream",3),1):
            has_rt=g["protein_key"] in best["ART_RT"] or g["protein_key"] in best["RVT_1"]
            neighbor_rows.append(dict(partner_key=p,upstream_rank=rank,**g,
                same_orientation=pg["strand"]==g["strand"],RT_profile_match=has_rt))
            if has_rt:
                routes[g["protein_key"]].add("partner_first");found=True
        if not found:orphans.append(dict(partner_key=p,accession=pg["accession"],status="typeI_partner_without_RT_in_first_three_genes"))
    table(out/"partner_first_neighbors.tsv",neighbor_rows)
    table(out/"partner_window_coordinate_maps.tsv",partner_windows)
    table(out/"partner_orphans.tsv",orphans,fields=["partner_key","accession","status"])
    acfg=json.loads((ROOT/"config/art_array_search.json").read_text())["art_array_search"]
    assess=json.loads((ROOT/"config/art_array_assessment.json").read_text())
    rows=[];sequence_rows=[];architecture_rows=[];pair_rows=[]
    rtout=[];partnerout=[];ntdout=[];repeatout=[];spacerout=[];arrayout=[]
    for key in sorted(routes):
        g=genes[key];rec=record_map[g["accession"]];rt=features[key]
        locus=g["accession"]+"__"+g["protein_id"]
        pairs=[]
        for rank,p in enumerate(neighbors(g,by_contig[g["accession"]],"downstream",3),1):
            if p["protein_key"] not in typeI:continue
            # Check rank in partner->RT direction too; rank need not be symmetric.
            pup=neighbors(p,by_contig[p["accession"]],"upstream",3)
            up_rank=next((i for i,v in enumerate(pup,1) if v["protein_key"]==key),None)
            intervening_rts = [v["protein_key"] for v in pup[:up_rank-1] if v["protein_key"] in routes] if up_rank else []
            compatible=p["strand"]==g["strand"] and up_rank is not None and not intervening_rts
            phs=[best[f][p["protein_key"]] for f in ("Type_I_partner","TypeI_two_seed") if p["protein_key"] in best[f]]
            pairs.append(dict(protein_key=p["protein_key"],protein_id=p["protein_id"],
                family="Type_I_partner",start0=p["start0"],end0=p["end0"],strand=p["strand"],
                RT_upstream_gene_rank=up_rank,partner_downstream_gene_rank=rank,
                intervening_RT_keys=intervening_rts,
                compatible_context=compatible,translation_qc=p["translation_qc"],profile_matches=phs))
        plan=plan_locus(rec,rt,[features[p["protein_key"]] for p in pairs])
        oriented=extract_record(rec,rt,plan,locus)
        for fmt,ext in (("genbank",".gb"),("fasta",".fna")):SeqIO.write(oriented,out/"loci"/(locus+ext),fmt)
        write_gff(out/"loci"/(locus+".gff3"),oriented)
        write_json(out/"loci"/(locus+".map.json"),dict(source_accession=rec.id,
            source_path=metadata[rec.id]["path"],source_sha256=metadata[rec.id]["file_sha256"],
            segments=[asdict(s) for s in plan["segments"]],
            metadata={k:v for k,v in plan.items() if k not in ("segments","plus_segments")}))
        rtf=next(f for f in oriented.features if pid(f)==g["protein_id"] and f.type=="CDS")
        rtstart=int(rtf.location.start)
        upstream=str(oriented.seq[max(0,rtstart-6000):rtstart]).upper()
        offset=max(0,rtstart-6000)
        initial=initial_scan(upstream[-3000:],acfg);extended=refinement_scan(upstream,acfg)
        # Use the maximum-score chain whose statistic is tested by shuffled maxima.
        model=aligned_pwm(upstream,extended["candidates"][0],acfg,out,locus) if extended["candidates"] else {}
        intervals=model.get("final_intervals",[])
        locus_intervals=[[a+offset,b+offset] for a,b in intervals]
        short=short_period_scan(upstream,acfg)
        flags=quality_flags(oriented,locus_intervals,rtstart,plan,acfg,assess,short,extended["score"]) if intervals else {"flags":[]}
        if len(upstream)<6000:flags["flags"].append("incomplete_6000nt_upstream")
        if any(c not in "ACGT" for c in upstream):flags["flags"].append("ambiguous_upstream_bases")
        if plan["source_origin_reached"] or plan["source_end_reached"]:flags["flags"].append("extracted_locus_reaches_contig_boundary")
        if "historical_v1.0" in metadata[rec.id].get("assembly_status",""):
            flags["flags"].append("historical_assembly_looping_risk")
        reps=[dict(start0=a,end0=b,sequence=str(oriented.seq[a:b]),
                   source_segments=source_pieces(a,b,plan["segments"])) for a,b in locus_intervals]
        spacers=[dict(start0=b,end0=c,sequence=str(oriented.seq[b:c]),length_nt=c-b,
                      source_segments=source_pieces(b,c,plan["segments"]))
                 for (a,b),(c,d) in zip(locus_intervals,locus_intervals[1:])]
        array_status="exploratory_array_requires_copy_review_and_nulls" if len(reps)>=3 else "no_supported_array_exploratory"
        if len(upstream)<6000 and len(reps)<3:array_status="incomplete_upstream_assessment"
        rvt=best["RVT_1"].get(key)
        ntd={} if not rvt else dict(start0=0,end0=rvt["protein_envelope_start0"],definition="before_first_significant_PF00078_envelope")
        row=dict(locus_id=locus,source_accession_and_version=rec.id,source_version_kind=metadata[rec.id].get("source_version_kind","accession_version"),
            database_release=metadata[rec.id].get("database_release",""),contig_id=rec.id,role=g["role"],
            rt_coordinates_and_strand=json.dumps(dict(start0=g["start0"],end0=g["end0"],strand=g["strand"])),
            RT_protein_id=g["protein_id"],protein_key=key,
            art_profile_score_and_coverage=json.dumps(best["ART_RT"].get(key,{})),
            rt_completeness=g["translation_qc"],ntd_boundary=json.dumps(ntd),
            partner_family_and_coordinates=json.dumps(pairs),array_assessment_status=array_status,
            repeat_coordinates_and_consensus=json.dumps(reps),spacer_lengths=json.dumps([s["length_nt"] for s in spacers]),
            upstream_sequence_available=len(upstream),phylogenetic_assignment="pending_sequence_only_tree",
            discovery_route=";".join(sorted(routes[key])),final_evidence_category="provisional_mapped_RT_locus",
            flags=";".join(flags["flags"]),array_score=extended["score"],repeat_copies=len(reps),
            source_path=metadata[rec.id]["path"],source_file_sha256=metadata[rec.id]["file_sha256"])
        rows.append(row)
        rtout.append(SeqRecord(protein_map[key].seq,id=key,description=locus))
        if ntd:ntdout.append(SeqRecord(protein_map[key].seq[:ntd["end0"]],id=key,description=locus))
        sequence_rows.append(dict(locus_id=locus,protein_key=key,ART_profile=row["art_profile_score_and_coverage"],
                                  RVT_1_profile=json.dumps(rvt),NTD=json.dumps(ntd),uses_architecture=False))
        architecture_rows.append({k:row[k] for k in ("locus_id","partner_family_and_coordinates","array_assessment_status",
            "repeat_coordinates_and_consensus","spacer_lengths","upstream_sequence_available","flags")})
        compatible_pairs=[p for p in pairs if p["compatible_context"]]
        primary_key=min(compatible_pairs,key=lambda p:p["partner_downstream_gene_rank"])["protein_key"] if compatible_pairs else None
        for p in pairs:
            if not p["compatible_context"]:continue
            partnerout.append(SeqRecord(protein_map[p["protein_key"]].seq,id=p["protein_key"],description=locus))
            pair_rows.append(dict(locus_id=locus,RT_key=key,partner_key=p["protein_key"],partner_id=p["protein_id"],
                source_genome=rec.id,role=g["role"],RT_translation_qc=g["translation_qc"],
                association="primary_nearest_compatible" if p["protein_key"]==primary_key else "alternative_nearby_homolog",
                partner_translation_qc=p["translation_qc"],array_status=array_status,
                pair_completeness="both_CDS_complete" if g["translation_qc"]==p["translation_qc"]=="complete_translation_matches" else "gene_review_required"))
        for i,rep in enumerate(reps,1):repeatout.append(SeqRecord(Seq(rep["sequence"]),id=f"{locus}.repeat{i}",description="provisional"))
        for i,sp in enumerate(spacers,1):spacerout.append(SeqRecord(Seq(sp["sequence"]),id=f"{locus}.spacer{i}",description="provisional"))
        if reps:arrayout.append(SeqRecord(oriented.seq[locus_intervals[0][0]:locus_intervals[-1][1]],id=locus,description="exploratory"))
        write_json(out/"array_details"/(locus+".json"),dict(locus_id=locus,upstream_sequence=upstream,
            upstream_offset0=offset,initial=initial,extended=extended,model=model,repeats=reps,spacers=spacers,
            quality=flags,statistical_score="extended_seed_discovery_best_chain_score",array_status=array_status))
    table(out/"candidate_loci.tsv",rows)
    table(out/"sequence_family_placement.tsv",sequence_rows)
    table(out/"observed_locus_architecture.tsv",architecture_rows)
    table(out/"typeI_pairs.tsv",pair_rows)
    for name,seqs in (("candidate_RT_proteins.faa",rtout),("typeI_partners.faa",partnerout),
        ("NTDs.faa",ntdout),("repeats.fna",repeatout),("spacers.fna",spacerout),("arrays_exploratory.fna",arrayout)):
        SeqIO.write(list({r.id:r for r in seqs}.values()),out/name,"fasta")
    write_json(out/"processing_summary.json",dict(source_contigs=len(records),proteins_and_ORFs=len(proteins),
        mapped_RT_loci=len(rows),typeI_pair_associations=len(pair_rows),
        primary_typeI_pairs=sum(r["association"]=="primary_nearest_compatible" for r in pair_rows),
        new_public_primary_pairs=sum(r["role"]!="LPJP1_control" and r["association"]=="primary_nearest_compatible" for r in pair_rows),
        exploratory_arrays=sum(int(r["repeat_copies"])>=3 for r in rows),
        validation_limits=validation_gate(),elapsed_seconds=round(time.perf_counter()-started,3)))
    print((out/"processing_summary.json").read_text())


def null_assessment(cfg,out):
    rows=read_tsv(out/"candidate_loci.tsv")
    acfg=json.loads((ROOT/"config/art_array_search.json").read_text())["art_array_search"]
    n=cfg["null_shuffles"]
    if 1/(n+1)*len(rows)>.05:
        raise ValueError("Insufficient empirical resolution for declared cohort")
    null_dir=out/"nulls";null_dir.mkdir(exist_ok=True)
    tests=[]
    with ProcessPoolExecutor(max_workers=cfg["threads"]) as pool:
        for row in rows:
            locus=row["locus_id"];detail=json.loads((out/"array_details"/(locus+".json")).read_text())
            sequence=detail["upstream_sequence"];scores={};screen={}
            # Screening may stop only for non-support. Any positive decision
            # requires the fixed 2,000-draw test; the screen cannot create a call.
            for stage_n in (cfg["screening_shuffles"],n):
                for kind in ("local_blocks","dinucleotide"):
                    path=null_dir/(locus+"."+kind+"."+str(stage_n)+".json")
                    provenance=dict(sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                        search_config=acfg,seed=cfg["random_seed"],n=stage_n,
                        detector_sha256=sha(ROOT/"scripts/art_array_detector.py"),
                        assessment_sha256=sha(ROOT/"scripts/art_array_assessment.py"))
                    if path.exists():
                        cached=json.loads(path.read_text())
                        if cached["provenance"]!=provenance:raise ValueError("Null cache inputs changed")
                        values=cached["draws"]
                    else:
                        seeds=[stable_seed(cfg["random_seed"],locus,kind,i) for i in range(stage_n)]
                        tasks=[(sequence,kind,acfg,seeds[i:i+20]) for i in range(0,stage_n,20)]
                        values=[v for batch in pool.map(sensitivity_batch,tasks) for v in batch]
                        write_json(path,dict(provenance=provenance,draws=values))
                    scores[kind]=empirical(float(row["array_score"]),[v for s,v in values])["empirical_p"]
                if stage_n==cfg["screening_shuffles"]:
                    screen=dict(scores)
                    # In this bounded pilot, refine every >=3-copy tract,
                    # including borderline screen results. Do not lose a small
                    # three-copy array because 200 draws have coarse resolution.
                    if int(row["repeat_copies"])<3:
                        break
            adjusted=min(1,max(scores.values())*len(rows))
            status="statistically_supported_array_requires_copy_and_family_review" if stage_n==n and adjusted<=.05 and int(row["repeat_copies"])>=3 else "not_supported_by_new_cohort_test"
            if int(row["upstream_sequence_available"])<6000:
                status="incomplete_upstream_assessment"
            row["array_assessment_status"]=status
            row["array_empirical_p_local"]=scores["local_blocks"]
            row["array_empirical_p_dinucleotide"]=scores["dinucleotide"]
            row["array_bonferroni_p"]=adjusted
            row["testing_family_loci"]=len(rows)
            initial_values=null_scores(sequence[-3000:],"initial",acfg,100,locus,"typeI_public_initial",pool)
            initial_p=empirical(detail["initial"]["R"],[v for s,v in initial_values])["empirical_p"]
            write_json(null_dir/(locus+".initial100.json"),dict(draws=initial_values,
                sequence_sha256=hashlib.sha256(sequence[-3000:].encode()).hexdigest(),search_config=acfg,
                detector_sha256=sha(ROOT/"scripts/art_array_detector.py"),empirical_p=initial_p))
            row["initial_scan_empirical_p"]=initial_p
            row["search_mode_label"]="extended_only" if status.startswith("statistically_supported") and initial_p>.05 else "initial_and_extended_examined"
            tests.append(dict(locus_id=locus,local_p=scores["local_blocks"],dinucleotide_p=scores["dinucleotide"],
                corrected_p=adjusted,cohort_loci=len(rows),shuffles_per_background=stage_n,status=status,
                screen_local_p=screen["local_blocks"],screen_dinucleotide_p=screen["dinucleotide"],
                initial_p=initial_p,screen_can_only_stop_for_non_support=True,
                full_discovery_repeated=True,complete_ART_locus_claim=False))
            print(json.dumps(tests[-1]),flush=True)
    table(out/"candidate_loci.tsv",rows)
    table(out/"background_tests.tsv",tests)
    latest={r["locus_id"]:r for r in rows}
    architecture=read_tsv(out/"observed_locus_architecture.tsv")
    for row in architecture:
        row["array_assessment_status"]=latest[row["locus_id"]]["array_assessment_status"]
    table(out/"observed_locus_architecture.tsv",architecture)


def diversity(cfg,out):
    run=Runner(out)
    fasta=out/"typeI_partners.faa"
    if not list(SeqIO.parse(fasta,"fasta")):
        write_json(out/"diversity_summary.json",dict(partner_sequences=0,clusters70=0))
        return
    folder=out/"diversity";folder.mkdir(exist_ok=True)
    prefix=folder/"partners70";scratch=folder/"mmseqs_tmp"
    run.run([run.binaries["mmseqs"],"easy-cluster",fasta,prefix,scratch,
        "--min-seq-id",cfg["partner_cluster_identity"],"-c",cfg["partner_cluster_coverage"],
        "--cov-mode","0","--threads",cfg["threads"],"-k","5","--split-memory-limit","1G"],
        folder/"MMseqs.log",expected=[str(prefix)+"_cluster.tsv"])
    cluster_of={}
    for line in Path(str(prefix)+"_cluster.tsv").read_text().splitlines():
        a,b=line.split("\t");cluster_of[b]=a
    pairs=read_tsv(out/"typeI_pairs.tsv")
    latest={r["locus_id"]:r for r in read_tsv(out/"candidate_loci.tsv")}
    for p in pairs:
        p["partner_cluster70"]=cluster_of[p["partner_key"]]
        p["array_status"]=latest[p["locus_id"]]["array_assessment_status"]
    table(out/"typeI_pairs.tsv",pairs)
    primary=[p for p in pairs if p["association"]=="primary_nearest_compatible"]
    table(out/"primary_typeI_pairs.tsv",primary)
    primary_ids={p["partner_key"] for p in primary}
    SeqIO.write([r for r in SeqIO.parse(fasta,"fasta") if r.id in primary_ids],out/"primary_typeI_partners.faa","fasta")
    write_json(out/"diversity_summary.json",dict(partner_sequences=len(cluster_of),
        clusters70=len(set(cluster_of.values())),source_genome_groups=len({p["source_genome"] for p in pairs}),
        primary_pairs=len(primary),primary_partner_clusters70=len({p["partner_cluster70"] for p in primary}),
        primary_public_pairs=sum(p["role"]!="LPJP1_control" for p in primary),
        independent_sample_count_not_established=True,identity=.7,coverage_both=.8,
        complete_pairs=sum(p["pair_completeness"]=="both_CDS_complete" for p in pairs)))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("operation",choices=("process","nulls","diversity"))
    ap.add_argument("--config",type=Path,default=ROOT/"config/art_type1_public.json")
    args=ap.parse_args();cfg=json.loads(args.config.read_text());out=ROOT/cfg["output"]
    {"process":process,"nulls":null_assessment,"diversity":diversity}[args.operation](cfg,out)


if __name__=="__main__":main()
