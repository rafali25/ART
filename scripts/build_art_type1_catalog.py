"""Combine primary pairs without double-counting controls or hiding cohort tests."""
from collections import defaultdict
import hashlib
import json
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from art_reconstruction_common import ROOT, Runner, sha, table
from extract_art_loci import read_tsv, write_json


def main():
    out=ROOT/"results/art_type1_catalog";out.mkdir(exist_ok=True)
    cohorts=[ROOT/"results/art_type1_public",ROOT/"results/art_type1_mgnify"]
    all_loci={};occurrences=defaultdict(list);samples={};rtseqs={};partners={};ntds={}
    for folder in cohorts:
        loci={r["locus_id"]:r for r in read_tsv(folder/"candidate_loci.tsv")}
        tests={r["locus_id"]:r for r in read_tsv(folder/"background_tests.tsv")}
        rs={r.id:r for r in SeqIO.parse(folder/"candidate_RT_proteins.faa","fasta")}
        ps={r.id:r for r in SeqIO.parse(folder/"primary_typeI_partners.faa","fasta")}
        ns={r.id:r for r in SeqIO.parse(folder/"NTDs.faa","fasta")}
        for lid,r in loci.items():
            occurrences[lid].append(dict(cohort=folder.name,local_test=tests[lid]))
            all_loci.setdefault(lid,dict(r,result_folder=folder.relative_to(ROOT).as_posix(),
                tested_shuffles=int(tests[lid]["shuffles_per_background"])))
        for p in read_tsv(folder/"primary_typeI_pairs.tsv"):
            lid=p["locus_id"];r=loci[lid]
            rt=rs[p["RT_key"]];partner=ps[p["partner_key"]]
            digest=hashlib.sha256((str(rt.seq)+"|"+str(partner.seq)).encode()).hexdigest()
            if lid in samples:
                assert samples[lid]["exact_pair_sequence_sha256"]==digest
                continue
            samples[lid]=dict(p,contig_id=r["contig_id"],result_folder=folder.relative_to(ROOT).as_posix(),
                upstream_nt=int(r["upstream_sequence_available"]),repeat_copies=int(r["repeat_copies"]),
                spacer_lengths=r["spacer_lengths"],flags=r["flags"],phylogenetic_assignment=r["phylogenetic_assignment"],
                rt_length_aa=len(rt),partner_length_aa=len(partner),exact_pair_sequence_sha256=digest,
                rt_replication_length_range_350_900=350<=len(rt)<=900,
                NTD_boundary=r["ntd_boundary"],RT_coordinates=r["rt_coordinates_and_strand"],
                partner_coordinates=r["partner_family_and_coordinates"])
            rtseqs[lid]=SeqRecord(rt.seq,id=lid,description=p["RT_key"]+" provisional_typeI_RT")
            partners[lid]=SeqRecord(partner.seq,id=lid,description=p["partner_id"]+" primary_partner")
            if p["RT_key"] in ns:ntds[lid]=SeqRecord(ns[p["RT_key"]].seq,id=lid,description="operational_NTD")
    # Global correction includes RT loci that lack a type I partner, not only
    # the successful pairs. A duplicated control is a single tested locus.
    m=len(all_loci)
    for lid,r in all_loci.items():
        r["cohort_tests_json"]=json.dumps(occurrences[lid])
        r["combined_testing_family_loci"]=m
        r["combined_array_p"]=min(1,max(float(r["array_empirical_p_local"]),float(r["array_empirical_p_dinucleotide"]))*m)
        statistically_supported=(r["tested_shuffles"]==2000 and int(r["repeat_copies"])>=3 and r["combined_array_p"]<=.05)
        r["combined_array_assessment"]=(
            "supported_tract_in_incomplete_upstream" if statistically_supported and int(r["upstream_sequence_available"])<6000 else
            "statistically_supported_tract_requires_copy_and_family_review" if statistically_supported else
            "incomplete_upstream_assessment" if int(r["upstream_sequence_available"])<6000 else "not_supported")
        if lid in samples:
            samples[lid]["combined_array_p"]=r["combined_array_p"]
            samples[lid]["combined_array_assessment"]=r["combined_array_assessment"]
    table(out/"all_RT_candidate_loci.tsv",list(all_loci.values()))
    SeqIO.write(list(rtseqs.values()),out/"RT_proteins.faa","fasta")
    SeqIO.write(list(partners.values()),out/"primary_partners.faa","fasta")
    SeqIO.write(list(ntds.values()),out/"NTDs.faa","fasta")
    run=Runner(out);prefix=out/"partners70"
    run.run([run.binaries["mmseqs"],"easy-cluster",out/"primary_partners.faa",prefix,out/"mmseqs_tmp",
        "--min-seq-id","0.7","-c","0.8","--cov-mode","0","--threads","4","-k","5","--split-memory-limit","1G"],
        out/"MMseqs.log",expected=[str(prefix)+"_cluster.tsv"])
    groups=defaultdict(list)
    for line in (out/"partners70_cluster.tsv").read_text().splitlines():
        center,member=line.split("\t");samples[member]["partner_cluster70"]=center
        groups[center].append(member)
    representatives=[]
    for center,members in groups.items():
        complete=[k for k in members if samples[k]["pair_completeness"]=="both_CDS_complete"]
        if not complete:continue
        best=min(complete,key=lambda k:(not samples[k]["combined_array_assessment"].startswith("statistically_supported"),
            bool(samples[k]["flags"]),-samples[k]["upstream_nt"],samples[k]["combined_array_p"],k))
        representatives.append(dict(samples[best],selection_reason="one_complete_pair_per_70pct_partner_cluster"))
    table(out/"typeI_samples.tsv",list(samples.values()))
    table(out/"complete_pair_representatives.tsv",representatives)
    supported=[r for r in samples.values() if r["combined_array_assessment"].startswith("statistically_supported")]
    array_sequences=[];repeat_sequences=[];spacer_sequences=[];array_rows=[]
    for r in supported:
        folder=ROOT/r["result_folder"];lid=r["locus_id"]
        d=json.loads((folder/"array_details"/(lid+".json")).read_text())
        locus=SeqIO.read(folder/"loci"/(lid+".fna"),"fasta")
        a=d["repeats"][0]["start0"];b=d["repeats"][-1]["end0"]
        array_sequences.append(SeqRecord(locus.seq[a:b],id=lid,description="statistically_supported_provisional_array_DNA"))
        for label,items,records in (("repeat",d["repeats"],repeat_sequences),("spacer",d["spacers"],spacer_sequences)):
            for i,x in enumerate(items,1):
                assert str(locus.seq[x["start0"]:x["end0"]])==x["sequence"]
                records.append(SeqRecord(locus.seq[x["start0"]:x["end0"]],id=f"{lid}.{label}{i}",
                    description="inferred_boundaries"))
        array_rows.append(dict(locus_id=lid,role=r["role"],locus_start0=a,locus_end0=b,
            length_nt=b-a,copies=r["repeat_copies"],spacer_gaps_nt=r["spacer_lengths"],
            corrected_p=r["combined_array_p"],flags=r["flags"],
            locus_genbank=(folder/"loci"/(lid+".gb")).relative_to(ROOT).as_posix(),
            coordinate_map=(folder/"loci"/(lid+".map.json")).relative_to(ROOT).as_posix(),
            array_details=(folder/"array_details"/(lid+".json")).relative_to(ROOT).as_posix()))
    table(out/"supported_array_candidates.tsv",array_rows)
    SeqIO.write(array_sequences,out/"supported_array_candidates.fna","fasta")
    SeqIO.write(repeat_sequences,out/"array_repeats.fna","fasta")
    SeqIO.write(spacer_sequences,out/"array_spacers.fna","fasta")
    public=[r for r in samples.values() if r["role"]!="LPJP1_control"]
    summary=dict(primary_loci=len(samples),new_public_primary_loci=len(public),
        exact_protein_pair_groups=len({r["exact_pair_sequence_sha256"] for r in samples.values()}),
        partner_clusters70=len(groups),complete_pair_representatives=len(representatives),
        complete_primary_pairs=sum(r["pair_completeness"]=="both_CDS_complete" for r in samples.values()),
        incomplete_gene_pairs=sum(r["pair_completeness"]!="both_CDS_complete" for r in samples.values()),
        supported_array_loci=len(supported),new_public_supported_array_loci=sum(r["role"]!="LPJP1_control" for r in supported),
        testing_family_distinct_RT_loci=m,all_statuses_provisional=True,
        incomplete_upstream_loci=sum(r["upstream_nt"]<6000 for r in samples.values()),
        input_manifest_sha256={str(p.relative_to(ROOT)):sha(p) for folder in cohorts for p in
            (folder/"candidate_loci.tsv",folder/"primary_typeI_pairs.tsv",folder/"background_tests.tsv")})
    write_json(out/"summary.json",summary)
    lines=["# Combined provisional type I ART sample catalog","",
        f"**{len(samples)} primary loci: {len(public)} new public candidates plus the two LPJP1 controls.** "
        f"The primary partners form {len(groups)} clusters at 70% identity and 80% mutual coverage. "
        "Different contigs and clusters are not proof of independent biological samples.","",
        f"{summary['complete_primary_pairs']} pairs have complete, translation-matched CDSs; "
        f"{summary['incomplete_gene_pairs']} require gene-boundary/completeness review. "
        f"`complete_pair_representatives.tsv` selects {len(representatives)} complete pairs, one per eligible partner cluster.","",
        f"Array evidence is corrected across **all {m} distinct RT candidate loci** tested in either cohort, "
        "including candidates without type I partners. This is a conditional comparison within the selected panels, "
        "not a database-wide false discovery rate. Duplicate LPJP1 controls are counted once.","",
        "| Locus | Copies | Spacer gaps, nt | Combined corrected p | Flags |",
        "|---|---:|---|---:|---|"]
    for r in supported:lines.append(f"| {r['locus_id']} | {r['repeat_copies']} | {r['spacer_lengths']} | {r['combined_array_p']:.5g} | {r['flags']} |")
    lines += ["",
        "Statistical support establishes a repeat tract, not ART function or correct full repeat boundaries. "
        "The largest gap in MGYC002687634053 may span more than one missing copy despite the permissive seed-chain fit; "
        "it is explicitly flagged for review. Two-copy and weak tracts remain exploratory.",
        "",
        "The type I profile still lacks verified independent positive holdout groups. NART comparisons and "
        "independent full-boundary validation remain incomplete. Partial genes are excluded from the "
        "complete-pair shortlist. Short upstream windows do not establish array absence; complete protein "
        "pairs may still be retained when their upstream DNA is incomplete.",
        "",
        "Use `typeI_samples.tsv` for the complete inventory, `complete_pair_representatives.tsv` for a "
        "diversity-oriented starting selection, and `RT_proteins.faa`, `primary_partners.faa`, `NTDs.faa` "
        "for sequences. Each row links through its `result_folder` and `locus_id` to oriented DNA, "
        "GenBank/GFF annotations, coordinate maps, CDS DNA, arrays, repeat alignments and spacers.",
        "",
        "`supported_array_candidates.fna` contains the four supported candidate tracts (three new plus LPJP1). "
        "`supported_array_candidates.tsv` links each to its coordinate map and locus annotation. "
        "`array_repeats.fna` and `array_spacers.fna` contain the inferred components. These are DNA "
        "extractions; transcription start/end sites and mature RNA boundaries have not been inferred.",
        "",
        "The two original result folders retain their own cohort corrections. Source manifests, "
        "exact genomic translation checks, raw HTTP responses, database releases, ETags and SHA-256 checksums "
        "are preserved. MGnify retrieval was capped at ten contigs; remaining hits are deferred, not negative.",
        "",
        "Rebuild with `.venv/Scripts/python.exe scripts/build_art_type1_catalog.py` after finalizing both cohorts."]
    (out/"README.md").write_text("\n".join(lines)+"\n",encoding="utf8")
    assert len(rtseqs)==len(partners)==len(samples)
    assert sum(r["role"]=="LPJP1_control" for r in samples.values())==2
    assert all(r["locus_id"] in all_loci for r in representatives)
    write_json(out/"verification.json",dict(passed=True,unique_primary_loci=len(samples),
        RT_and_partner_sequences=len(samples),duplicate_LPJP1_controls_removed=True,
        exported_array_DNA_records=len(array_sequences),repeat_and_spacer_slices_checked=len(repeat_sequences)+len(spacer_sequences),
        statistics_include_all_RT_candidates=True,source_cohort_verifications=[
            json.loads((f/"verification.json").read_text()) for f in cohorts]))
    write_json(out/"output_sha256.json",{p.relative_to(out).as_posix():sha(p) for p in out.rglob("*")
        if "mmseqs_tmp" not in p.parts and p.name!="output_sha256.json" and p.is_file()})
    print(json.dumps(summary,indent=2))


if __name__=="__main__":main()
