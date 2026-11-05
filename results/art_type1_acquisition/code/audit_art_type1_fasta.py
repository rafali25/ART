"""Bounded annotation-rescue audit of local DNA, with LPJP1 as a retrieval control.

This does not train profiles, classify new loci, or run a new array-significance
cohort. Source annotations and the existing assessments remain unchanged.
All nucleotide coordinates are zero-based, half-open.
"""
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import time

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pyrodigal
import pyrodigal_gv

from art_reconstruction_common import Runner, sha, table
from extract_art_loci import ROOT, pid, read_tsv, write_json

OUT = ROOT / "results/art_type1_acquisition"


def six_frame_orfs(sequence, code=11, minimum_aa=40):
    """Maximal translated segments between stops/X, without requiring a start.

    These are retrieval fragments, not CDS annotations. Break at ambiguous
    translations so an ambiguous codon cannot falsely join two fragments.
    Coordinates exclude the terminating stop, and translation is reversible.
    """
    sequence = Seq(str(sequence).upper())
    length = len(sequence)
    for strand in (1, -1):
        oriented = sequence if strand == 1 else sequence.reverse_complement()
        for frame in range(3):
            usable = (length - frame) // 3 * 3
            peptide = str(oriented[frame:frame + usable].translate(table=code))
            for match in re.finditer(r"[^*X]+", peptide):
                if len(match.group()) < minimum_aa:
                    continue
                a = frame + 3 * match.start()
                b = frame + 3 * match.end()
                start, end = (a, b) if strand == 1 else (length - b, length - a)
                yield dict(start0=start, end0=end, strand=strand,
                           frame_in_oriented_sequence=frame, genetic_code=code,
                           peptide=match.group(),
                           boundary_5="sequence_edge" if match.start() == 0 else
                           "stop" if peptide[match.start() - 1] == "*" else "ambiguous",
                           boundary_3="sequence_edge" if match.end() == len(peptide) else
                           "stop" if peptide[match.end()] == "*" else "ambiguous")


def protein_interval_to_source(row, aa_start, aa_end):
    """Map an alignment in an unspliced translation to source DNA."""
    if not 0 <= aa_start <= aa_end <= row["length_aa"]:
        raise ValueError("Protein interval outside the translated segment")
    if row["strand"] == 1:
        return row["start0"] + 3 * aa_start, row["start0"] + 3 * aa_end
    return row["end0"] - 3 * aa_end, row["end0"] - 3 * aa_start


def matching_annotations(row, annotations):
    """Report same-frame overlap, without turning an ORF fragment into a gene."""
    matches = []
    for feature in annotations:
        if feature.type != "CDS" or len(feature.location.parts) != 1:
            continue
        a, b = int(feature.location.start), int(feature.location.end)
        if feature.location.strand != row["strand"]:
            continue
        offset = int(feature.qualifiers.get("codon_start", ["1"])[0]) - 1
        a, b = (a + offset, b) if row["strand"] == 1 else (a, b - offset)
        delta = row["start0"] - a if row["strand"] == 1 else row["end0"] - b
        overlap = max(0, min(b, row["end0"]) - max(a, row["start0"]))
        shorter = min(b - a, row["end0"] - row["start0"])
        if shorter and delta % 3 == 0 and overlap / shorter >= .8:
            matches.append((overlap, pid(feature)))
    return ";".join(p for _, p in sorted(matches, reverse=True))


def package_current_type1(records):
    """Export the two known LPJP1 pairs and the one existing array candidate."""
    folder = OUT / "current_typeI"
    folder.mkdir(exist_ok=True)
    rec = records["MZ422438.1"]
    features = {pid(f): f for f in rec.features if f.type == "CDS"}
    master = read_tsv(ROOT / "results/art_expansion/candidate_loci.tsv")
    rows = [r for r in master if r["locus_id"] in
            {"MZ422438.1__QXN67953.1", "MZ422438.1__QXN67955.1"}]
    if len(rows) != 2:
        raise ValueError("Expected the two existing LPJP1 locus records")
    rt_proteins, partners, ntds, cds, upstream, arrays, units = [], [], [], [], [], [], []
    pairs = []
    for row in rows:
        locus_id = row["locus_id"]
        rt_id = row["RT_protein_id"]
        ps = [p for p in json.loads(row["partner_family_and_coordinates"])
              if p["compatible_adjacent_context"] and p["family"] == "Type_I_partner"]
        if len(ps) != 1:
            raise ValueError("Expected one immediately adjacent compatible type I partner")
        partner_id = ps[0]["protein_id"]
        rtf, pf = features[rt_id], features[partner_id]
        rtseq = rtf.qualifiers["translation"][0]
        pseq = pf.qualifiers["translation"][0]
        for feature, translation in ((rtf, rtseq), (pf, pseq)):
            code = int(feature.qualifiers.get("transl_table", ["11"])[0])
            assert str(feature.extract(rec.seq).translate(table=code, cds=True)) == translation
            cds.append(SeqRecord(feature.extract(rec.seq), id=pid(feature),
                                 description=f"{locus_id} source_CDS_in_transcription_direction"))
        rt_proteins.append(SeqRecord(Seq(rtseq), id=rt_id, description=locus_id))
        partners.append(SeqRecord(Seq(pseq), id=partner_id, description=locus_id))
        ntd = json.loads(row["ntd_boundary"])
        ntds.append(SeqRecord(Seq(rtseq[ntd["start0"]:ntd["end0"]]), id=rt_id + "_NTD",
                             description=f"aa0={ntd['start0']}:{ntd['end0']} PF00078_proxy"))
        locus_path = ROOT / row["oriented_locus_path"]
        oriented = SeqIO.read(locus_path, "genbank")
        oriented_rt = next(f for f in oriented.features if pid(f) == rt_id and f.type == "CDS")
        rtstart = int(oriented_rt.location.start)
        upstream.append(SeqRecord(oriented.seq[max(0, rtstart - 6000):rtstart],
                                  id=locus_id, description="RT_oriented_upstream"))
        for suffix in (".gb", ".fna", ".gff3", ".map.json"):
            path = locus_path.parent / (locus_id + suffix)
            if path.exists():
                shutil.copy2(path, folder / path.name)
        repeat_info = json.loads(row["repeat_coordinates_and_consensus"])
        repeats = repeat_info["repeats"]
        if repeats and row["array_assessment_status"] == "accepted_array_candidate":
            start, end = repeats[0]["locus_start0"], repeats[-1]["locus_end0"]
            arrays.append(SeqRecord(oriented.seq[start:end], id=locus_id,
                                   description="existing_computational_array_candidate"))
            for i, (left, right) in enumerate(zip(repeats, repeats[1:]), 1):
                units.append(SeqRecord(oriented.seq[left["locus_start0"]:right["locus_end0"]],
                                       id=f"{locus_id}_unit{i}",
                                       description="repeat_spacer_repeat_DNA_not_measured_transcript"))
        pairs.append(dict(locus_id=locus_id, RT=rt_id, partner=partner_id,
                          RT_length_aa=len(rtseq), partner_length_aa=len(pseq),
                          NTD_end0=ntd["end0"], array_status=row["array_assessment_status"],
                          repeat_copies=len(repeats), spacer_lengths=row["spacer_lengths"],
                          independent_genome_group="LPJP1",
                          seed_status="publication_based_reconstruction_not_independent_validation"))
    for name, seqs in (("RT_proteins.faa", rt_proteins), ("partners.faa", partners),
                       ("NTDs.faa", ntds), ("CDS.fna", cds), ("upstream.fna", upstream),
                       ("array_candidates.fna", arrays), ("repeat_spacer_repeat_units.fna", units)):
        SeqIO.write(seqs, folder / name, "fasta")
    table(folder / "pairs.tsv", pairs)
    table(folder / "locus_evidence.tsv", rows)
    for name in ("repeats.tsv", "spacers.tsv"):
        selected = [r for r in read_tsv(ROOT / "results/art_arrays_assessed" / name)
                    if r["locus_id"] in {x["locus_id"] for x in rows}]
        table(folder / name, selected)
    return pairs


def main():
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "searches").mkdir(exist_ok=True)
    sources = read_tsv(ROOT / "results/art_discovery/source_records.tsv")
    records, source_rows, by_hash = {}, [], {}
    for row in sources:
        path = ROOT / row["path"]
        if sha(path) != row["file_sha256"]:
            raise ValueError("Original source changed")
        rec = SeqIO.read(path, row["format"])
        records[rec.id] = rec
        by_hash[hashlib.sha256(str(rec.seq).upper().encode()).hexdigest()] = rec.id
    fasta_rows = []
    for path in sorted((ROOT / "datasetss").rglob("*")):
        if path.suffix.lower() not in {".fasta", ".fa", ".fna"}:
            continue
        for rec in SeqIO.parse(path, "fasta"):
            digest = hashlib.sha256(str(rec.seq).upper().encode()).hexdigest()
            if digest not in by_hash:
                raise ValueError("New unmatched local FASTA requires an explicit source mapping")
            fasta_rows.append(dict(path=path.relative_to(ROOT).as_posix(),
                                   record_id=rec.id, length_nt=len(rec),
                                   file_sha256=sha(path), sequence_sha256=digest,
                                   matching_annotated_accession=by_hash[digest],
                                   previously_searched_identical_DNA=True))
    table(OUT / "local_FASTA_inventory.tsv", fasta_rows)
    mapping, proteins, issues = [], [], []
    finder = pyrodigal_gv.ViralGeneFinder(meta=True)
    for source in sources:
        rec = records[source["accession"]]
        role = "LPJP1_retrieval_control" if rec.id == "MZ422438.1" else "local_genome"
        predicted = finder.find_genes(bytes(rec.seq))
        codes = sorted({11, *(g.translation_table for g in predicted)})
        count_before = len(mapping)

        def add(row, peptide, method, label):
            row = dict(row, accession=rec.id, role=role, method=method,
                       source_label=label, length_aa=len(peptide))
            key = "audit_" + hashlib.sha256(
                f"{rec.id}|{method}|{label}|{row['start0']}|{row['end0']}|{row['strand']}".encode()
            ).hexdigest()[:20]
            row["target"] = key
            row["source_CDS_overlap"] = matching_annotations(row, rec.features)
            mapping.append(row)
            proteins.append(SeqRecord(Seq(peptide), id=key, description=f"{rec.id} {method} {label}"))

        for i, f in enumerate(rec.features):
            if f.type != "CDS":
                continue
            peptide = f.qualifiers.get("translation", [""])[0].rstrip("*")
            if not peptide or len(f.location.parts) != 1:
                issues.append(dict(accession=rec.id, feature=pid(f),
                                   reason="missing_translation_or_compound_location"))
                continue
            if set(peptide) - set("ACDEFGHIKLMNPQRSTVWYXBZJUO"):
                raise ValueError("Unsupported original translation")
            a, b = int(f.location.start), int(f.location.end)
            offset = int(f.qualifiers.get("codon_start", ["1"])[0]) - 1
            a, b = (a + offset, b) if f.location.strand == 1 else (a, b - offset)
            add(dict(start0=a, end0=b, strand=f.location.strand,
                     genetic_code=int(f.qualifiers.get("transl_table", ["11"])[0])),
                peptide, "source_annotation", pid(f) or str(i))
        for i, gene in enumerate(predicted, 1):
            add(dict(start0=gene.begin - 1, end0=gene.end, strand=gene.strand,
                     genetic_code=gene.translation_table,
                     partial_begin=gene.partial_begin, partial_end=gene.partial_end),
                gene.translate().rstrip("*"), "pyrodigal_gv", f"pgv{i}")
        for code in codes:
            for i, orf in enumerate(six_frame_orfs(rec.seq, code), 1):
                peptide = orf.pop("peptide")
                # Independent genomic retranslation verifies every fragment.
                dna = rec.seq[orf["start0"]:orf["end0"]]
                if orf["strand"] == -1:
                    dna = dna.reverse_complement()
                assert str(dna.translate(table=code)) == peptide
                add(orf, peptide, "six_frame", f"code{code}_orf{i}")
        source_rows.append(dict(source, role=role, selected_gene_model=str(predicted.metagenomic_bin),
                                predicted_CDS=len(predicted), six_frame_codes=json.dumps(codes),
                                translations_searched=len(mapping) - count_before,
                                circular_boundary_translation="not_joined_in_this_bounded_audit"))
    keys = {r["target"]: r for r in mapping}
    assert len(keys) == len(mapping)
    table(OUT / "source_records.tsv", source_rows)
    table(OUT / "protein_to_contig.tsv", mapping)
    table(OUT / "annotation_issues.tsv", issues)
    SeqIO.write(proteins, OUT / "proteins.faa", "fasta")
    run = Runner(OUT)
    hits, profile_hashes = [], {}
    for family in ("ART_RT", "Type_I_partner", "Type_II_partner", "RVT_1"):
        path = run.binaries["RVT_1"] if family == "RVT_1" else (
            ROOT / "results/art_discovery/profiles" / (family + ".hmm"))
        profile_hashes[family] = sha(path)
        found, _ = run.hmm(path, OUT / "proteins.faa", OUT / "searches" / family,
                           z=len(proteins), report_e=.001)
        for h in found:
            row = keys[h["target"]]
            a, b = protein_interval_to_source(row, h["protein_alignment_start0"],
                                              h["protein_alignment_end0"])
            hits.append(dict(h, family=family, **{k: row[k] for k in (
                "accession", "role", "method", "source_label", "source_CDS_overlap",
                "start0", "end0", "strand", "genetic_code", "length_aa")},
                alignment_source_start0=a, alignment_source_end0=b,
                passes_provisional_retrieval=h["sequence_evalue"] <= 1e-5
                and h["domain_evalue"] <= 1e-5 and h["profile_coverage"] >= .6))
    table(OUT / "all_reported_profile_hits.tsv", hits)
    passed = [h for h in hits if h["passes_provisional_retrieval"]]
    table(OUT / "passing_profile_hits.tsv", passed)
    counts = []
    for accession in records:
        for family in ("ART_RT", "Type_I_partner", "Type_II_partner", "RVT_1"):
            for method in ("source_annotation", "pyrodigal_gv", "six_frame"):
                selected = [h for h in passed if h["accession"] == accession
                            and h["family"] == family and h["method"] == method]
                counts.append(dict(accession=accession, family=family, method=method,
                                   passing_translations=len({h["target"] for h in selected}),
                                   mapped_existing_CDS=";".join(sorted({p for h in selected
                                       for p in h["source_CDS_overlap"].split(";") if p})),
                                   unmapped_translation_hits=sum(not h["source_CDS_overlap"] for h in selected)))
    table(OUT / "recovery_by_genome.tsv", counts)
    pairs = package_current_type1(records)
    summary = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                   local_FASTA_files=len({r["path"] for r in fasta_rows}),
                   local_FASTA_records=len(fasta_rows),
                   all_local_FASTA_identical_to_previous_DNA=True,
                   local_unique_genomes_searched=sum(s["role"] == "local_genome" for s in source_rows),
                   additional_LPJP1_control_genomes=1,
                   translations_searched=len(proteins),
                   translation_methods=dict(Counter(r["method"] for r in mapping)),
                   typeI_local_reported_translations=len({h["target"] for h in hits
                       if h["family"] == "Type_I_partner" and h["role"] == "local_genome"}),
                   typeI_local_passing_translations=len({h["target"] for h in passed
                       if h["family"] == "Type_I_partner" and h["role"] == "local_genome"}),
                   local_unmapped_ART_or_typeI_hits=[h for h in passed
                       if h["role"] == "local_genome" and h["family"] in {"ART_RT", "Type_I_partner"}
                       and not h["source_CDS_overlap"]],
                   current_typeI_pairs=pairs, new_array_tests_executed=False,
                   profiles_retrained=False, profile_sha256=profile_hashes,
                   pyrodigal_gv_version=pyrodigal_gv.__version__,
                   pyrodigal_version=pyrodigal.__version__,
                   limitations=["single-seed type I profile; no diverse positive calibration",
                       "maximal six-frame fragments are not complete gene annotations",
                       "linear translation does not rescue genes spanning circular sequence boundaries",
                       "recovery is provisional; no biological absence or new sample claim"],
                   elapsed_seconds=round(time.perf_counter() - started, 3))
    write_json(OUT / "summary.json", summary)
    for source in sources:
        assert sha(ROOT / source["path"]) == source["file_sha256"]
    (OUT / "code").mkdir(exist_ok=True)
    shutil.copy2(Path(__file__), OUT / "code" / Path(__file__).name)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
