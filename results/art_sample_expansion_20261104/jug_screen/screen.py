"""Bounded protein-profile screen of deposited curated Jug phage genomes.

Coordinates are zero-based, half-open; genomic source DNA is never edited.
This script performs gene calling and homology retrieval, not array inference.
"""
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pyrodigal
import pyrodigal_gv

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent


def sha(data):
    return hashlib.sha256(data).hexdigest()


def table(path, rows, fields=None):
    fields = fields or list(dict.fromkeys(k for row in rows for k in row)) or ['status']
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def call_genes(item):
    index, genome_id, dna = item
    finder = pyrodigal_gv.ViralGeneFinder(meta=True)
    called = finder.find_genes(dna.encode())
    rows = []
    for number, gene in enumerate(called, 1):
        a, b, strand = gene.begin - 1, gene.end, gene.strand
        assert 0 <= a < b <= len(dna)
        cds = Seq(dna[a:b])
        if strand == -1:
            cds = cds.reverse_complement()
        text = gene.translate().rstrip('*')
        computed = str(cds.translate(table=gene.translation_table)).rstrip('*')
        # Prodigal begin/end partial flags refer to increasing genomic
        # coordinates. The coding 5' end is the genomic end on a minus gene.
        partial_5prime = gene.partial_begin if strand == 1 else gene.partial_end
        partial_3prime = gene.partial_end if strand == 1 else gene.partial_begin
        if not partial_5prime and computed:
            computed = 'M' + computed[1:]
        compatible = len(computed) == len(text) and all(a == b or b == 'X' for a,b in zip(computed,text))
        assert compatible, (genome_id, number, str(cds[:3]), computed[:30], text[:30])
        rows.append(dict(protein_id=f'jug{index:04d}_pgv{number:04d}', genome_id=genome_id,
                         start0=a, end0=b, strand=strand, genetic_code=gene.translation_table,
                         partial_begin=gene.partial_begin, partial_end=gene.partial_end,
                         partial_5prime=partial_5prime, partial_3prime=partial_3prime,
                         protein_length_aa=len(text), protein_sha256=sha(text.encode()),
                         CDS_sha256=sha(str(cds).encode()), translation_matches=computed == text,
                         translation_compatible=compatible,
                         translation_qc='exact_match' if computed == text else 'caller_X_at_ambiguous_codon',
                         ambiguous_residues=sum(c not in 'ACDEFGHIKLMNPQRSTVWY' for c in text),
                         translation=text))
    return rows, dict(genome_id=genome_id, genome_length_nt=len(dna),
                      sequence_sha256=sha(dna.encode()), genes=len(rows),
                      selected_model=str(called.metagenomic_bin),
                      genetic_codes_json=json.dumps(dict(Counter(g['genetic_code'] for g in rows))),
                      caller='Pyrodigal-gv', caller_version=pyrodigal_gv.__version__,
                      pyrodigal_version=pyrodigal.__version__, mode='ViralGeneFinder(meta=True)')


def search_profile(item):
    family, profile, protein_count = item
    stem = OUT / 'searches' / family
    args = [str(ROOT / '.venv/bin/hmmsearch'), '--cpu', '1', '--seed', '42', '--noali',
            '-Z', str(protein_count), '--domZ', str(protein_count), '-E', '0.001', '--domE', '0.001',
            '--tblout', str(stem) + '.tbl', '--domtblout', str(stem) + '.domtbl',
            str(profile), str(OUT / 'proteins.faa')]
    start = time.perf_counter()
    with Path(str(stem) + '.log').open('w') as log:
        result = subprocess.run(args, stdout=log, stderr=subprocess.STDOUT)
    record = dict(family=family, args=args, returncode=result.returncode,
                  elapsed_seconds=time.perf_counter() - start, profile_path=str(profile.relative_to(ROOT)),
                  profile_sha256=sha(profile.read_bytes()))
    if result.returncode:
        raise RuntimeError(json.dumps(record))
    rows = []
    for line in Path(str(stem) + '.domtbl').read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        c = line.split(maxsplit=22)
        coverage = (int(c[16]) - int(c[15]) + 1) / int(c[5])
        seq_e, domain_e = float(c[6]), float(c[12])
        rows.append(dict(family=family, protein_id=c[0], profile_name=c[3],
                         sequence_evalue=seq_e, domain_evalue=domain_e,
                         sequence_bits=float(c[7]), domain_bits=float(c[13]),
                         profile_coverage=coverage, hmm_start0=int(c[15])-1, hmm_end0=int(c[16]),
                         protein_start0=int(c[17])-1, protein_end0=int(c[18]),
                         envelope_start0=int(c[19])-1, envelope_end0=int(c[20]),
                         passes_frozen_screen=seq_e <= 1e-5 and domain_e <= 1e-5 and coverage >= 0.6,
                         interpretation='generic_RT_domain_only' if family == 'RVT_1' else
                         'exploratory_two_same_genome_paralog_profile' if family == 'TypeI_two_seed_exploratory' else
                         'production_profile_retrieval_not_validated_family_assignment'))
    return rows, record


def main():
    started = time.perf_counter()
    (OUT / 'searches').mkdir(exist_ok=True)
    source = OUT / 'sources/deposited/139.curated.Mac.phage.genomes.fa'
    records = list(SeqIO.parse(source, 'fasta'))
    assert len(records) == 139, len(records)
    assert len({r.id for r in records}) == len(records), 'Duplicate source IDs'
    assert all(set(str(r.seq).upper()) <= set('ACGTRYSWKMBDHVN') for r in records)
    genes, genomes = [], []
    jobs = [(i, r.id, str(r.seq).upper()) for i, r in enumerate(records, 1)]
    with ProcessPoolExecutor(max_workers=4) as pool:
        for i, (rows, genome) in enumerate(pool.map(call_genes, jobs), 1):
            genes.extend(rows)
            genomes.append(genome)
            if i % 20 == 0:
                print(f'Called {i}/139 genomes; {len(genes)} proteins', flush=True)
    gene_seconds = time.perf_counter() - started
    SeqIO.write([SeqRecord(Seq(g['translation']), id=g['protein_id'],
                          description=f"{g['genome_id']} {g['start0']}:{g['end0']}:{g['strand']}")
                 for g in genes], OUT / 'proteins.faa', 'fasta')
    gene_fields = [k for k in genes[0] if k != 'translation']
    table(OUT / 'gene_coordinates.tsv', [{k:g[k] for k in gene_fields} for g in genes], gene_fields)
    table(OUT / 'genome_manifest.tsv', genomes)
    with (OUT / 'predicted_genes.gff3').open('w') as handle:
        handle.write('##gff-version 3\n')
        for r in records:
            handle.write(f'##sequence-region {r.id} 1 {len(r)}\n')
        for g in genes:
            attrs = f"ID={g['protein_id']};protein_id={g['protein_id']};transl_table={g['genetic_code']};partial_begin={int(g['partial_begin'])};partial_end={int(g['partial_end'])}"
            handle.write(f"{g['genome_id']}\tPyrodigal-gv\tCDS\t{g['start0']+1}\t{g['end0']}\t.\t{'+' if g['strand']==1 else '-'}\t0\t{attrs}\n")
    profiles = [('ART_RT', ROOT / 'results/art_discovery/profiles/ART_RT.hmm'),
                ('Type_I_partner', ROOT / 'results/art_discovery/profiles/Type_I_partner.hmm'),
                ('Type_II_partner', ROOT / 'results/art_discovery/profiles/Type_II_partner.hmm'),
                ('TypeI_two_seed_exploratory', ROOT / 'results/art_type1_public/profiles/TypeI_two_seed_exploratory.hmm'),
                ('RVT_1', ROOT / '.tools/art-reconstruction/PF00078.hmm')]
    hits, commands = [], []
    with ThreadPoolExecutor(max_workers=4) as pool:
        for rows, command in pool.map(search_profile, [(family, path, len(genes)) for family, path in profiles]):
            hits.extend(rows)
            commands.append(command)
    lookup = {g['protein_id']:g for g in genes}
    merged = [dict(h, **{k:lookup[h['protein_id']][k] for k in gene_fields if k != 'protein_id'}) for h in hits]
    fields = list(merged[0]) if merged else ['family','protein_id','passes_frozen_screen']
    table(OUT / 'all_profile_domain_hits.tsv', merged, fields)
    passing = [h for h in merged if h['passes_frozen_screen']]
    table(OUT / 'production_ART_candidates.tsv', [h for h in passing if h['family']=='ART_RT'], fields)
    table(OUT / 'partner_candidates.tsv', [h for h in passing if h['family'] in ('Type_I_partner','Type_II_partner','TypeI_two_seed_exploratory')], fields)
    table(OUT / 'generic_RT_candidates.tsv', [h for h in passing if h['family']=='RVT_1'], fields)
    selected = {h['protein_id'] for h in passing}
    SeqIO.write([SeqRecord(Seq(lookup[p]['translation']), id=p, description=lookup[p]['genome_id'])
                 for p in sorted(selected)], OUT / 'passing_profile_proteins.faa', 'fasta')
    # Preserve near-neighbor coordinates for review without inferring an ART locus.
    neighbors = []
    for protein in sorted(selected):
        anchor = lookup[protein]
        same_genome = sorted([g for g in genes if g['genome_id']==anchor['genome_id']], key=lambda g:(g['start0'],g['end0']))
        position = next(i for i,g in enumerate(same_genome) if g['protein_id']==protein)
        for other in same_genome[max(0,position-3):position+4]:
            if other['protein_id']==protein:
                continue
            neighbors.append(dict(anchor_protein=protein, neighbor_protein=other['protein_id'],
                                  genome_id=other['genome_id'], start0=other['start0'], end0=other['end0'],
                                  strand=other['strand'], same_orientation=other['strand']==anchor['strand'],
                                  source_order_difference=next(i for i,g in enumerate(same_genome) if g['protein_id']==other['protein_id'])-position,
                                  passing_profiles=';'.join(sorted({h['family'] for h in passing if h['protein_id']==other['protein_id']}))))
    table(OUT / 'candidate_neighbors.tsv', neighbors,
          ['anchor_protein','neighbor_protein','genome_id','start0','end0','strand','same_orientation','source_order_difference','passing_profiles'])
    write_json(OUT / 'commands.json', commands)
    summary = dict(genomes=len(records), unique_genome_sequences=len({g['sequence_sha256'] for g in genomes}),
                   nucleotide_nt=sum(len(r) for r in records), proteins=len(genes),
                   genetic_codes=dict(Counter(g['genetic_code'] for g in genes)),
                   caller_models=dict(Counter(g['selected_model'] for g in genomes)),
                   partial_proteins=sum(g['partial_begin'] or g['partial_end'] for g in genes),
                   translation_exact_differences=sum(not g['translation_matches'] for g in genes),
                   translation_unexplained_mismatches=sum(not g['translation_compatible'] for g in genes),
                   reported_domains_by_family=dict(Counter(h['family'] for h in hits)),
                   passing_domains_by_family={family:sum(h['family']==family for h in passing) for family,path in profiles},
                   passing_proteins_by_family={family:len({h['protein_id'] for h in passing if h['family']==family}) for family,path in profiles},
                   retrieval_evalue=0.001, sequence_evalue=1e-5, domain_evalue=1e-5,
                   minimum_profile_coverage=0.6, hmm_Z=len(genes), hmm_domZ=len(genes),
                   gene_call_seconds=round(gene_seconds,3), total_seconds=round(time.perf_counter()-started,3),
                   workers=4, pyrodigal_gv_version=pyrodigal_gv.__version__, pyrodigal_version=pyrodigal.__version__,
                   source_fasta_sha256=sha(source.read_bytes()), script_sha256=sha(Path(__file__).read_bytes()),
                   hmmsearch_sha256=sha((ROOT/'.venv/bin/hmmsearch').read_bytes()),
                   array_inference_performed=False, validated_ART_locus_claims=0,
                   interpretation='Bounded additional-dataset protein retrieval; generic RT hits do not establish ART family or locus membership')
    write_json(OUT / 'summary.json', summary)
    write_json(OUT / 'output_sha256.json', {str(p.relative_to(OUT)):sha(p.read_bytes())
                                         for p in OUT.rglob('*') if p.is_file() and p.name!='output_sha256.json'})
    print(json.dumps(summary,indent=2),flush=True)


if __name__ == '__main__':
    main()
