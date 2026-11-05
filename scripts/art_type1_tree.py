"""Independent RT-domain tree using the frozen, labelled comparison panel."""
from collections import Counter
import json
from Bio import AlignIO, Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from art_reconstruction_common import ROOT, Runner, sha, table
from extract_art_loci import read_tsv, write_json


def main(cfg=None):
    cfg=cfg or json.loads((ROOT/"config/art_type1_public.json").read_text())
    out=ROOT/cfg["output"];folder=out/"domain_phylogeny";folder.mkdir(exist_ok=True)
    frozen=ROOT/"results/art_family_reconstruction/domain_phylogeny"
    metadata={r["protein_key"]:r for r in read_tsv(frozen/"tree_tip_annotations.tsv")
              if r["role"] in ("protein_only_labelled_reference","published_ART_reference")}
    fragments=[r for r in SeqIO.parse(frozen/"domains.faa","fasta") if r.id in metadata]
    queries=list(SeqIO.parse(out/"candidate_RT_proteins.faa","fasta"))
    hits=read_tsv(out/"profile_hits.tsv");first={}
    for h in hits:
        if h["family"]!="RVT_1" or float(h["domain_evalue"])>.001:continue
        if h["target"] not in first or int(h["protein_envelope_start0"])<int(first[h["target"]]["protein_envelope_start0"]):
            first[h["target"]]=h
    coordinates=[]
    for r in queries:
        metadata[r.id]=dict(protein_key=r.id,original_id=r.description,broad_class="query",role="mapped_candidate")
        if r.id not in first:continue
        h=first[r.id];a=max(0,int(h["protein_envelope_start0"])-40);b=min(len(r),int(h["protein_envelope_end0"])+40)
        fragments.append(SeqRecord(r.seq[a:b],id=r.id,description=""))
        coordinates.append(dict(protein_key=r.id,start0=a,end0=b))
    SeqIO.write(fragments,folder/"domains.faa","fasta");table(folder/"query_domain_coordinates.tsv",coordinates)
    run=Runner(folder);alignment=folder/"domains.mafft.faa"
    run.run(["cmd.exe","/d","/c",ROOT/".tools/mafft-win/mafft.bat","--localpair","--maxiterate","1000",
        "--thread","4","--inputorder","--out",alignment,folder/"domains.faa"],folder/"MAFFT.log",expected=[alignment])
    msa=AlignIO.read(alignment,"fasta")
    keep=[j for j in range(msa.get_alignment_length()) if sum(r.seq[j]!="-" for r in msa)>=.5*len(msa)]
    if len(keep)<50:raise ValueError("Insufficient domain alignment")
    trimmed=[SeqRecord(Seq("".join(str(r.seq[j]) for j in keep)),id=r.id,description="") for r in msa]
    SeqIO.write(trimmed,folder/"domains.trimmed.faa","fasta")
    write_json(folder/"alignment_columns.json",dict(input_columns0=keep,occupancy=.5))
    treepath=folder/"RT_domain.FastTree.nwk"
    run.run([run.binaries["fasttree"],"-lg","-gamma","-seed","42",folder/"domains.trimmed.faa"],
        folder/"FastTree.log",stdout=treepath,expected=[treepath])
    tree=Phylo.read(treepath,"newick");tips={t.name:t for t in tree.get_terminals()}
    refs=[k for k in tips if metadata[k]["broad_class"]!="query"]
    placement=[]
    for r in queries:
        if r.id not in tips:
            placement.append(dict(protein_key=r.id,assignment="RVT_domain_unavailable",uses_architecture=False))
            continue
        near=sorted([(tree.distance(tips[r.id],tips[k]),k) for k in refs])[:5]
        votes=Counter(metadata[k]["broad_class"] for d,k in near[:3]);label,n=votes.most_common(1)[0]
        placement.append(dict(protein_key=r.id,
            assignment="ART_reference_affinity_NART_unresolved" if label=="ART" else label+"_affinity_provisional",
            nearest_three_agreement=n/3,uses_architecture=False,
            nearest_references_json=json.dumps([dict(protein_key=k,label=metadata[k]["broad_class"],
                original_id=metadata[k]["original_id"],distance=d) for d,k in near])))
    table(folder/"sequence_only_placement.tsv",placement)
    table(folder/"tree_tip_annotations.tsv",[metadata[k] for k in tips])
    lookup={r["protein_key"]:r["assignment"] for r in placement}
    loci=read_tsv(out/"candidate_loci.tsv")
    for r in loci:r["phylogenetic_assignment"]=lookup[r["protein_key"]]
    table(out/"candidate_loci.tsv",loci)
    seq=read_tsv(out/"sequence_family_placement.tsv")
    for r in seq:r["phylogenetic_assignment"]=lookup[r["protein_key"]]
    table(out/"sequence_family_placement.tsv",seq)
    write_json(folder/"summary.json",dict(tips=len(tips),columns=len(keep),query_count=len(queries),
        reference_domains_sha256=sha(frozen/"domains.faa"),
        reference_labels_sha256=sha(frozen/"tree_tip_annotations.tsv"),NART_references=0,
        scope="Descriptive sequence-only affinity; not calibrated family classification",
        method="PF00078 envelope +/-40aa; MAFFT L-INS-i; 50% occupancy; FastTree LG+Gamma"))
    print((folder/"summary.json").read_text())


if __name__=="__main__":
    import argparse
    from pathlib import Path
    ap=argparse.ArgumentParser()
    ap.add_argument("--config",type=Path,default=ROOT/"config/art_type1_public.json")
    args=ap.parse_args()
    main(json.loads(args.config.read_text()))
