"""Independently check the small type I package and annotation-rescue output."""
import json
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from art_reconstruction_common import sha
from extract_art_loci import ROOT, pid, read_tsv, write_json


def main():
    out = ROOT / "results/art_type1_acquisition"
    summary = json.loads((out / "summary.json").read_text())
    records = {}
    for row in read_tsv(out / "source_records.tsv"):
        assert sha(ROOT / row["path"]) == row["file_sha256"]
        records[row["accession"]] = SeqIO.read(ROOT / row["path"], row["format"])
    translations = {s.id: str(s.seq) for s in SeqIO.parse(out / "proteins.faa", "fasta")}
    mappings = {r["target"]: r for r in read_tsv(out / "protein_to_contig.tsv")}
    assert set(translations) == set(mappings)
    checked = 0
    for key, row in mappings.items():
        if row["method"] != "six_frame":
            continue
        source = records[row["accession"]].seq[int(row["start0"]):int(row["end0"])]
        if row["strand"] == "-1":
            source = source.reverse_complement()
        assert str(source.translate(table=int(row["genetic_code"]))) == translations[key]
        checked += 1
    profile_coordinates = 0
    for hit in read_tsv(out / "all_reported_profile_hits.tsv"):
        source = records[hit["accession"]].seq[
            int(hit["alignment_source_start0"]):int(hit["alignment_source_end0"])]
        if hit["strand"] == "-1":
            source = source.reverse_complement()
        expected = translations[hit["target"]][
            int(hit["protein_alignment_start0"]):int(hit["protein_alignment_end0"])]
        computed = str(source.translate(table=int(hit["genetic_code"])))
        # CDS translations convert a bacterial alternative start codon to M.
        assert computed == expected or (
            int(hit["protein_alignment_start0"]) == 0
            and expected.startswith("M") and computed[1:] == expected[1:]
            and hit["method"] != "six_frame")
        profile_coordinates += 1
    folder = out / "current_typeI"
    original = records["MZ422438.1"]
    features = {pid(f): f for f in original.features if f.type == "CDS"}
    for name in ("RT_proteins.faa", "partners.faa"):
        seqs = list(SeqIO.parse(folder / name, "fasta"))
        assert len(seqs) == 2
        for seq in seqs:
            assert str(seq.seq) == features[seq.id].qualifiers["translation"][0]
    for seq in SeqIO.parse(folder / "CDS.fna", "fasta"):
        assert seq.seq == features[seq.id].extract(original.seq)
    loci = {}
    for path in folder.glob("*.map.json"):
        mapping = json.loads(path.read_text())
        seq = Seq("")
        for segment in mapping["segments"]:
            fragment = original.seq[segment["source_start0"]:segment["source_end0"]]
            if segment["orientation"] == -1:
                fragment = fragment.reverse_complement()
            seq += fragment
        locus_id = path.name.removesuffix(".map.json")
        loci[locus_id] = SeqIO.read(folder / (locus_id + ".gb"), "genbank")
        assert seq == loci[locus_id].seq
        assert seq == SeqIO.read(folder / (locus_id + ".fna"), "fasta").seq
    repeat_checks = 0
    for filename in ("repeats.tsv", "spacers.tsv"):
        for row in read_tsv(folder / filename):
            locus = loci[row["locus_id"]]
            assert str(locus.seq[int(row["locus_start0"]):int(row["locus_end0"])]) == row["sequence"]
            seq = Seq("")
            for segment in json.loads(row["source_segments_json"]):
                fragment = original.seq[segment["source_start0"]:segment["source_end0"]]
                if segment["orientation"] == -1:
                    fragment = fragment.reverse_complement()
                seq += fragment
            assert str(seq) == row["sequence"]
            repeat_checks += 1
    assert summary["typeI_local_passing_translations"] == 0
    assert not summary["local_unmapped_ART_or_typeI_hits"]
    recovery = read_tsv(out / "recovery_by_genome.tsv")
    for method in ("source_annotation", "pyrodigal_gv", "six_frame"):
        control = next(r for r in recovery if r["accession"] == "MZ422438.1"
                       and r["family"] == "Type_I_partner" and r["method"] == method)
        assert set(control["mapped_existing_CDS"].split(";")) == {"QXN67952.1", "QXN67954.1"}
        local_rt = [r for r in recovery if r["accession"] != "MZ422438.1"
                    and r["family"] == "ART_RT" and r["method"] == method]
        assert len(local_rt) == 7 and all(r["passing_translations"] == "1" for r in local_rt)
    result = dict(passed=True, six_frame_source_translations_checked=checked,
                  profile_hit_coordinate_mappings_checked=profile_coordinates,
                  oriented_loci_checked=len(loci), repeat_and_spacer_source_mappings=repeat_checks,
                  LPJP1_control_recovered_in_all_three_layers=True,
                  seven_local_known_RT_genes_recovered_in_all_three_layers=True,
                  new_typeI_local_candidates=0)
    write_json(out / "verification.json", result)
    inputs = [ROOT / "literature/Yoon et al. 2026 Array Associated Reverse Transcriptases (2).pdf",
              ROOT / "literature/ARTNT CLUADE_v1.docx",
              ROOT / "scripts/audit_art_type1_fasta.py",
              ROOT / "scripts/test_art_type1_fasta.py", Path(__file__)]
    write_json(out / "input_sha256.json", {p.relative_to(ROOT).as_posix(): sha(p) for p in inputs})
    outputs = {p.relative_to(out).as_posix(): sha(p) for p in sorted(out.rglob("*"))
               if p.is_file() and p.name != "output_sha256.json"}
    write_json(out / "output_sha256.json", outputs)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
