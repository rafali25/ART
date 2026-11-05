"""Independent sequence/coordinate verification of the public type I pilot."""
import json
from Bio import SeqIO
from Bio.Seq import Seq
from art_reconstruction_common import ROOT, sha
from extract_art_loci import read_tsv, write_json, pid


def main(cfg=None):
    cfg=cfg or json.loads((ROOT/"config/art_type1_public.json").read_text())
    out=ROOT/cfg["output"];sources={}
    for r in read_tsv(out/"source_records.tsv"):
        assert sha(ROOT/r["path"])==r["file_sha256"]
        sources[r["accession"]]=SeqIO.read(ROOT/r["path"],r["format"])
    proteins={r.id:r.seq for r in SeqIO.parse(out/"proteins.faa","fasta")}
    checked=0
    for g in read_tsv(out/"six_frame_map.tsv"):
        seq=sources[g["accession"]].seq[int(g["start0"]):int(g["end0"])]
        if g["strand"]=="-1":seq=seq.reverse_complement()
        assert str(seq.translate(table=int(g["genetic_code"])))==str(proteins[g["protein_key"]])
        checked+=1
    def source_sequence(parts,source):
        seq=Seq("")
        for s in parts:
            fragment=source.seq[s["source_start0"]:s["source_end0"]]
            if s["orientation"]==-1:fragment=fragment.reverse_complement()
            seq+=fragment
        return seq
    loci=read_tsv(out/"candidate_loci.tsv");repeat_checks=0;cds_checks=0
    for r in loci:
        lid=r["locus_id"];record=SeqIO.read(out/"loci"/(lid+".gb"),"genbank")
        mapping=json.loads((out/"loci"/(lid+".map.json")).read_text())
        source=sources[r["contig_id"]]
        assert source_sequence(mapping["segments"],source)==record.seq
        assert SeqIO.read(out/"loci"/(lid+".fna"),"fasta").seq==record.seq
        rt=next(f for f in record.features if f.type=="CDS" and pid(f)==r["RT_protein_id"])
        assert rt.location.strand==1
        if r["rt_completeness"]=="complete_translation_matches":
            assert rt.extract(record.seq).translate(table=int(rt.qualifiers.get("transl_table",["11"])[0]),cds=True)==proteins[r["protein_key"]]
            cds_checks+=1
        d=json.loads((out/"array_details"/(lid+".json")).read_text())
        assert d["array_status"]==r["array_assessment_status"]
        assert str(record.seq[d["upstream_offset0"]:int(rt.location.start)])==d["upstream_sequence"]
        for x in d["repeats"]+d["spacers"]:
            assert str(record.seq[x["start0"]:x["end0"]])==x["sequence"]
            assert str(source_sequence(x["source_segments"],source))==x["sequence"]
            repeat_checks+=1
        for left,right,spacer in zip(d["repeats"],d["repeats"][1:],d["spacers"]):
            assert spacer["length_nt"]==right["start0"]-left["end0"]
    primary=read_tsv(out/"primary_typeI_pairs.tsv")
    assert len(primary)==len({r["locus_id"] for r in primary})
    control=[r for r in primary if r["role"]=="LPJP1_control"]
    assert {(r["locus_id"].split("__")[-1],r["partner_id"]) for r in control}=={
        ("QXN67953.1","QXN67952.1"),("QXN67955.1","QXN67954.1")}
    placements={r["protein_key"]:r["assignment"] for r in read_tsv(out/"domain_phylogeny/sequence_only_placement.tsv")}
    assert all(r["phylogenetic_assignment"]==placements[r["protein_key"]] for r in loci)
    for r in read_tsv(out/"background_tests.tsv"):
        assert int(r["cohort_loci"])==len(loci)
        assert abs(float(r["corrected_p"])-min(1,max(float(r["local_p"]),float(r["dinucleotide_p"]))*len(loci)))<1e-10
        if r["status"].startswith("statistically_supported"):assert int(r["shuffles_per_background"])==2000
    reconciliation=[]
    if (out/"archive_leads/current_assembly_reconciliation.json").exists():
        reconciliation=json.loads((out/"archive_leads/current_assembly_reconciliation.json").read_text())
        assert len(reconciliation)==4 and all(r["status"]=="current_mapping_confirmed" for r in reconciliation)
    mapped_remote=0
    if (out/"homolog_to_contig.tsv").exists():
        peptides={r.id:str(r.seq) for r in SeqIO.parse(out/"remote_homologs.faa","fasta")}
        for m in read_tsv(out/"homolog_to_contig.tsv"):
            if m["status"]!="exact_translated_mapping":continue
            intervals=json.loads(m["coordinates_json"]);assert len(intervals)==1
            x=intervals[0];seq=sources[m["accession"]].seq[x["start0"]:x["end0"]]
            if x["strand"]==-1:seq=seq.reverse_complement()
            observed=str(seq.translate(table=x["genetic_code"]));expected=peptides[m["protein_accession"]].rstrip("*")
            assert observed==expected or (x["initiator_normalized"] and expected.startswith("M") and observed[1:]==expected[1:])
            mapped_remote+=1
    summary=dict(passed=True,source_records=len(sources),oriented_loci=len(loci),
        six_frame_translations=checked,complete_RT_CDS_translations=cds_checks,
        repeat_and_spacer_source_mappings=repeat_checks,primary_pairs=len(primary),
        current_assembly_exact_peptide_mappings=len(reconciliation),
        remote_protein_DNA_mappings=mapped_remote,
        LPJP1_pairing_regression_passed=True)
    write_json(out/"verification.json",summary)
    inputs=[out/"processing_config.json"]+list((ROOT/"scripts").glob("*art_type1*.py"))
    inputs += [
        ROOT/"config/art_array_search.json",ROOT/"config/art_array_assessment.json",
        ROOT/"scripts/art_array_detector.py",ROOT/"scripts/art_array_assessment.py",
        ROOT/"scripts/extract_art_loci.py",ROOT/"scripts/art_discovery.py",
        ROOT/"scripts/art_reconstruction_common.py",
        ROOT/"results/art_discovery/profiles/ART_RT.hmm",
        ROOT/"results/art_discovery/profiles/Type_I_partner.hmm",
        ROOT/"results/art_discovery/profiles/Type_II_partner.hmm",
        ROOT/".tools/myRT/Models/HMM/RVT-All.hmm",
        ROOT/"results/art_family_reconstruction/domain_phylogeny/domains.faa",
        ROOT/"results/art_family_reconstruction/domain_phylogeny/tree_tip_annotations.tsv",
        ROOT/"literature/Yoon et al. 2026 Array Associated Reverse Transcriptases (2).pdf",
        ROOT/"literature/ARTNT CLUADE_v1.docx",
        ROOT/"results/art_discovery/typeIII_source_audit/ai3_core.tar.gz"]
    write_json(out/"input_sha256.json",{p.relative_to(ROOT).as_posix():sha(p) for p in inputs})
    write_json(out/"output_sha256.json",{p.relative_to(out).as_posix():sha(p) for p in out.rglob("*")
        if "mmseqs_tmp" not in p.parts and p.name!="output_sha256.json" and p.is_file()})
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    import argparse
    from pathlib import Path
    ap=argparse.ArgumentParser()
    ap.add_argument("--config",type=Path,default=ROOT/"config/art_type1_public.json")
    main(json.loads(ap.parse_args().config.read_text()))
