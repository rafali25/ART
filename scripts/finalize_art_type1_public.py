"""Consolidate completed analyses, export useful sequences, and write a report."""
import json
import statistics
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from art_reconstruction_common import ROOT, Runner, passes, sha, table
from extract_art_loci import read_tsv, write_json, pid


def main(cfg=None):
    cfg=cfg or json.loads((ROOT/"config/art_type1_public.json").read_text())
    out=ROOT/cfg["output"]
    rows=read_tsv(out/"candidate_loci.tsv")
    placement={r["protein_key"]:r for r in read_tsv(out/"domain_phylogeny/sequence_only_placement.tsv")}
    backgrounds={r["locus_id"]:r for r in read_tsv(out/"background_tests.tsv")}
    run=Runner(out)
    hits,_=run.hmm(ROOT/".tools/myRT/Models/HMM/RVT-All.hmm",out/"candidate_RT_proteins.faa",
        out/"searches/competing_myRT",z=len(rows))
    table(out/"competing_RT_family_matches.tsv",hits)
    families={}
    for h in hits:families.setdefault(h["target"],[]).append(h)
    details={};repeat_rows=[];spacer_rows=[];cds=[];upstream=[]
    primary=read_tsv(out/"primary_typeI_pairs.tsv")
    pairmap={p["locus_id"]:p for p in primary}
    for r in rows:
        r["phylogenetic_assignment"]=placement[r["protein_key"]]["assignment"]
        r["competing_RT_profiles"]=json.dumps(families.get(r["protein_key"],[]))
        supported=r["array_assessment_status"].startswith("statistically_supported")
        if r["locus_id"] in pairmap:
            r["final_evidence_category"]=(
                "provisional_typeI_pair_with_statistically_supported_array" if supported else
                "provisional_typeI_pair_array_assessment_incomplete" if int(r["upstream_sequence_available"])<6000 else
                "provisional_typeI_pair_without_supported_array")
        path=out/"array_details"/(r["locus_id"]+".json")
        d=json.loads(path.read_text());d["array_status"]=r["array_assessment_status"];d["background_test"]=backgrounds[r["locus_id"]]
        base_spacing=d["extended"]["candidates"][0]["base_spacing"] if d["extended"]["candidates"] else 0
        gaps=[b["start0"]-a["start0"] for a,b in zip(d["repeats"],d["repeats"][1:])]
        typical_spacing=statistics.median(gaps) if gaps else 0
        if len(gaps)>=3 and max(gaps)>2.6*typical_spacing:
            flag="gap_may_span_multiple_missing_repeats"
            if flag not in d["quality"]["flags"]:d["quality"]["flags"].append(flag)
            r["flags"]=";".join(dict.fromkeys(r["flags"].split(";")+[flag])).strip(";")
        aligned=d.get("model",{}).get("final_aligned_sequences",[])
        r["repeat_consensus"]=("".join(Counter(column).most_common(1)[0][0] for column in zip(*aligned))
                              if aligned else "")
        r["inferred_repeat_lengths_nt"]=json.dumps([len(x["sequence"]) for x in d["repeats"]])
        for i,x in enumerate(d["spacers"]):
            start_gap=d["repeats"][i+1]["start0"]-d["repeats"][i]["start0"]
            x["possible_missing_repeat"]=bool(base_spacing and 1.4*base_spacing<=start_gap<=2.6*base_spacing)
            x["spacing_based_missing_copy_estimate"]=max(0,round(start_gap/typical_spacing)-1) if typical_spacing else 0
        write_json(path,d);details[r["locus_id"]]=d
        for label,items,target in (("repeat",d["repeats"],repeat_rows),("spacer",d["spacers"],spacer_rows)):
            for i,x in enumerate(items,1):
                target.append(dict(locus_id=r["locus_id"],copy_index=i,locus_start0=x["start0"],locus_end0=x["end0"],
                    sequence=x["sequence"],length_nt=len(x["sequence"]),source_segments_json=json.dumps(x["source_segments"]),
                    possible_missing_repeat=x.get("possible_missing_repeat",False) if label=="spacer" else "",
                    spacing_based_missing_copy_estimate=x.get("spacing_based_missing_copy_estimate",""),
                    status=r["array_assessment_status"],boundary_status="inferred_requires_independent_confirmation"))
        locus=SeqIO.read(out/"loci"/(r["locus_id"]+".gb"),"genbank")
        feature_ids={r["RT_protein_id"]} if r["locus_id"] in pairmap else set()
        if r["locus_id"] in pairmap:feature_ids.add(pairmap[r["locus_id"]]["partner_id"])
        for f in locus.features:
            if f.type=="CDS" and pid(f) in feature_ids:
                cds.append(SeqRecord(f.extract(locus.seq),id=r["locus_id"]+"__"+pid(f),description="source_oriented_CDS"))
        upstream.append(SeqRecord(Seq(d["upstream_sequence"]),id=r["locus_id"],
            description="RT_transcription_orientation available_nt="+r["upstream_sequence_available"]))
    table(out/"candidate_loci.tsv",rows)
    latest={r["locus_id"]:r for r in rows}
    architecture=read_tsv(out/"observed_locus_architecture.tsv")
    for r in architecture:r["flags"]=latest[r["locus_id"]]["flags"]
    table(out/"observed_locus_architecture.tsv",architecture)
    table(out/"repeats.tsv",repeat_rows);table(out/"spacers.tsv",spacer_rows)
    SeqIO.write(cds,out/"primary_pair_CDS.fna","fasta")
    SeqIO.write(upstream,out/"upstream_sequences.fna","fasta")
    sequence=read_tsv(out/"sequence_family_placement.tsv")
    for r in sequence:
        r["phylogenetic_assignment"]=placement[r["protein_key"]]["assignment"]
        r["competing_RT_profiles"]=json.dumps(families.get(r["protein_key"],[]))
    table(out/"sequence_family_placement.tsv",sequence)
    diversity=json.loads((out/"diversity_summary.json").read_text())
    calibration_source=ROOT/"results/art_type1_public/archive_leads"
    controls={r["protein_id"] for r in read_tsv(calibration_source/"control_labels.tsv")}
    cross=read_tsv(calibration_source/"profile_cross_screen.tsv")
    passing={}
    for family in ("Type_I_partner","TypeI_two_seed"):
        passing[family]=sorted({h["target"] for h in cross if h["family"]==family and
            float(h["sequence_evalue"])<=cfg["sequence_evalue"] and
            float(h["domain_evalue"])<=cfg["domain_evalue"] and
            float(h["profile_coverage"])>=cfg["minimum_ART_profile_coverage"]})
    calibration=dict(independent_verified_positive_groups=0,training_source_groups=1,
        unrelated_GNAT_controls=len(controls),
        reported_GNAT_control_matches=[h for h in cross if h["target"] in controls],
        archived_leads_passing={f:len(ids) for f,ids in passing.items()},
        labels_are_external_leads_not_gold_standard=True,status="exploratory_not_broadly_calibrated")
    write_json(out/"profiles/independent_screen.json",calibration)
    write_json(out/"processing_config.json",cfg)
    public=[r for r in rows if r["role"]!="LPJP1_control" and r["locus_id"] in pairmap]
    jobs_path=ROOT/"results/art_type1_public/jobs.json"
    summary=dict(new_public_primary_loci=len(public),control_loci=sum(r["role"]=="LPJP1_control" for r in rows),
        all_RT_candidate_loci=len(rows),
        total_primary_pairs=len(primary),new_statistically_supported_arrays=sum(
            r["array_assessment_status"].startswith("statistically_supported") for r in public),
        new_incomplete_upstream=sum(int(r["upstream_sequence_available"])<6000 for r in public),
        primary_partner_clusters70=diversity["primary_partner_clusters70"],
        final_classification="provisional",software_test_log="software_tests.log",
        remote_jobs=json.loads(jobs_path.read_text()),calibration=calibration)
    write_json(out/"summary.json",summary)
    lines=[
        "# Type I ART acquisition pilot",
        "",
        "Implemented and run on four public environmental contigs plus LPJP1. "
        "Four new provisional primary RT–type I partner pairs were recovered, giving six primary pairs "
        "in this pilot. These are new to this workspace, not newly discovered biological families. "
        "The four public targets were supplied as leads by an external archive, not a representative prevalence sample.",
        "",
        "| Current contig | RT | Partner | Available upstream | Array assessment | Spacer lengths, nt |",
        "|---|---|---|---:|---|---|",
    ]
    for r in rows:
        p=pairmap.get(r["locus_id"])
        if not p:continue
        status=("supported statistically; provisional boundaries" if r["array_assessment_status"].startswith("statistically_supported")
            else "incomplete upstream" if int(r["upstream_sequence_available"])<6000 else "not supported")
        lines.append(f"| {r['contig_id']} | {r['RT_protein_id']} | {p['partner_id']} | "
            f"{r['upstream_sequence_available']} | {status} | {r['spacer_lengths']} |")
    lines += [
        "",
        f"The six primary partners form **{diversity['primary_partner_clusters70']} clusters** at 70% identity "
        "and 80% mutual coverage. Different contigs or protein clusters do not establish independent biological samples. "
        "Two additional nearby partner homolog associations are retained separately; they are not counted as primary pairs.",
        "",
        "## What was implemented",
        "",
        "- Reproducible MAFFT/HMMER profile construction and independent ART-first and partner-first retrieval. "
        "The additional two-seed type I model uses both LPJP1 paralogs; the original profiles remain unchanged.",
        "- Original annotations, Pyrodigal-gv predictions and six-frame ORF rescue; exact protein-to-DNA mapping. "
        "First-three-upstream-gene checks include orientation. An intervening RT prevents assigning a partner across that RT.",
        "- Partner ±6 kb windows, RT-oriented loci, source coordinate maps, GenBank/GFF annotations, "
        "RT/partner/NTD proteins, CDS DNA, upstream DNA, repeat alignments, repeats and spacers.",
        "- Independent 3 kb initial and 6 kb extended repeat searches. Spacers use inferred full repeat ends; "
        "published length ranges are descriptive, not rejection filters.",
        "- Full seed rediscovery and best-chain search on 200 local-block and 200 dinucleotide shuffles per locus. "
        "Potentially supported candidates advance to fixed 2,000-draw tests per background. "
        "The 200-draw screen can stop only for non-support; positive status requires the fixed 2,000-draw test. "
        "The larger background p-value is Bonferroni-adjusted across all six loci. "
        "This controls the declared pilot comparison, not the earlier external lead-selection process or a database-wide discovery rate. "
        "Each locus also has 100 nucleotide shuffles for the initial scan.",
        "- A separate 97-tip RT-domain tree with labelled RT classes, plus competing myRT-family HMM matches. "
        "The tree uses sequence only: PF00078 ±40 aa, MAFFT L-INS-i, 50% occupancy trimming, FastTree LG+Gamma. "
        "All six queries have descriptive ART-reference affinity; the panel cannot resolve ART versus NART.",
        "",
        "## New array",
        "",
        "ERR8055547_17635 has five inferred 16-nt repeats (`CCCTAGGATCTAGGGT`) spanning 1,609 nt. "
        "Source repeat intervals (zero-based, half-open) are [17737,17753), [18037,18053), "
        "[18379,18395), [19025,19041), [19330,19346). "
        "Spacers are 284, 326, 630 and 289 nt; the large gap is consistent with a missed/absent repeat. "
        "There are no detected CDS overlaps or genes between this tract and the RT. "
        "Neither background produced a score at least as high in 2,000 shuffles, so each empirical p is 1/2001 "
        "and the six-locus corrected p is approximately 0.003. This supports the repetitive tract, not its biological function.",
        "",
        "## Acquisition and provenance",
        "",
        "The existing local FASTA rescue audit found no extra type I loci. Public leads came from "
        "[the AI3 supporting archive](https://doi.org/10.5281/zenodo.23086637), "
        "whose 89-sequence homolog list is treated as unverified discovery input. "
        "Thirty leads pass the original type I profile and 35 pass the two-seed profile. "
        "Five archived unrelated GNAT controls had no reported profile matches; this small set does not establish broad specificity.",
        "",
        "All four target peptides were exactly mapped in Logan v1.0 and independently located in the current contigs. "
        "Contig numbers changed; numeric identifiers alone were rejected as mapping evidence. "
        "Three current DNA sequences are identical to their historical contigs; SRR26805139_340 is longer "
        "(16,150 nt versus 5,048 nt), with the partner still exactly encoded. "
        "Current sequences, object ETags, compressed prefixes, download metadata and SHA-256 checksums are retained. "
        "Logan run/contig IDs are not represented as INSDC accession versions. "
        "The bounded historical/current retrieval used 59,602,523 compressed bytes in total.",
        "",
        "The four EBI HMMER jobs have completed; NCBI BLAST remains queued at the latest poll. "
        "The MGnify follow-up is stored separately in `results/art_type1_mgnify/`; "
        "its larger testing cohort does not replace this six-locus pilot. Resume with the commands in README.md.",
        "",
        "## Limits and next selection",
        "",
        "No expanded profile has been trained on the new candidates. The exploratory type I model still has one "
        "training-genome group and no verified independent positive holdout group. "
        "Known synthetic repeat-copy precision failures and incomplete full-boundary validation remain unresolved. "
        "NART references and both verified type III partner families remain unavailable. "
        "Two new loci have only 232/233 nt upstream and cannot establish array absence. "
        "The other array-negative public locus has 6 kb upstream but no supported tract under this detector. "
        "Keep all these labels provisional; prioritize the ERR8055547 locus for independent boundary review and "
        "use primary partner clusters to select distinct protein pairs for the proposal.",
        "",
        "The proposal's DIAMOND/MMseqs2/FAMSA/FastTree replication workflow remains separate in "
        "`results/art_family_reconstruction/`; this bounded acquisition pilot does not claim to reconstruct "
        "the entire published ART family.",
    ]
    if out.name!="art_type1_public":
        sources=read_tsv(out/"source_records.tsv")
        lines=[
            "# MGnify follow-up to the Type I ART pilot","",
            f"Processed {sum(r['role']!='LPJP1_control' for r in sources)} public contigs plus LPJP1. "
            f"Recovered {len(public)} provisional public primary type I RT–partner loci. "
            f"All RT candidates ({len(rows)}, including controls and RTs without type I partners) "
            "belong to the declared array-testing cohort. These are candidate loci, not validated ART systems.","",
            "| Contig | Primary partner | Upstream nt | Copies | Spacer gaps nt | Corrected p | Assessment |",
            "|---|---|---:|---:|---|---:|---|"]
        for r in rows:
            if r["locus_id"] not in pairmap:continue
            p=pairmap[r["locus_id"]]
            lines.append(f"| {r['contig_id']} | {p['partner_id']} | {r['upstream_sequence_available']} | "
                f"{r['repeat_copies']} | {r['spacer_lengths']} | {r['array_bonferroni_p']} | {r['array_assessment_status']} |")
        lines += ["",
            "The profile searches cover both RT-first and type I partner-first routes. Original annotations are retained; "
            "unannotated DNA is called with Pyrodigal-gv and searched in six frames. Every downloaded query homolog "
            "has a recorded exact translated DNA mapping. MGYC identifiers are release snapshots with sequence checksums.",
            "",
            "Both 3 kb and 6 kb array searches run independently. The 200-draw screen can only stop for non-support; "
            "positive statistical status requires 2,000 local-block and 2,000 dinucleotide shuffles with complete "
            "seed rediscovery and chain selection. The larger p is corrected across all candidate loci in this cohort. "
            "Short upstream windows remain explicitly incomplete even when they contain a repeat tract.",
            "",
            "The separate RT-domain tree and competing myRT HMM searches do not use array or partner architecture. "
            "ART/NART separation remains unresolved without verified NART references. Type I profiles have one "
            "training-genome group and no independent verified positive holdout group. The five GNAT negatives and "
            "archived protein cross-screen are inherited from the pilot, not a new validation experiment. "
            "Known repeat-copy precision and full-boundary validation limitations still apply.",
            "",
            "Acquisition is bounded by the recorded config and source manifest: the first ten retrievable contigs "
            "from query-balanced ranked hits were retained. Unretrieved hits are deferred, not negative. "
            "The four MGnify searches produced 418 collected query-hit rows (293 unique MGYP proteins); "
            "both RT queries were capped at 200 results. Each partner query returned nine hits. "
            "NCBI BLAST was still queued at the latest poll. This is not a prevalence estimate.",
            "",
            "Use `primary_typeI_pairs.tsv`, `candidate_loci.tsv`, `background_tests.tsv`, `loci/`, "
            "`primary_typeI_partners.faa`, `candidate_RT_proteins.faa`, `NTDs.faa`, `primary_pair_CDS.fna`, "
            "`upstream_sequences.fna`, `repeats.tsv`, `spacers.tsv`, and `verification.json`. "
            "Coordinates are zero-based/half-open internally. Spacer intervals flagged `possible_missing_repeat` "
            "may span an undetected copy and should not be interpreted as a confirmed single biological spacer.",
            "",
            "Run the same sequential processing commands as the pilot with `--config config/art_type1_mgnify.json`. "
            "The original pilot results and its six-locus statistical correction remain separate."]
    (out/"REPORT.md").write_text("\n".join(lines)+"\n",encoding="utf8")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    import argparse
    from pathlib import Path
    ap=argparse.ArgumentParser()
    ap.add_argument("--config",type=Path,default=ROOT/"config/art_type1_public.json")
    main(json.loads(ap.parse_args().config.read_text()))
