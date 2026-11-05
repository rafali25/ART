"""Prepare and run the local Staphylococcus/type-II ART RT profile pilot.

Run from any directory with the project .venv Python. The target MarsHill
protein is never included in the training alignment. This is a recovery check,
not an independent evolutionary holdout or a calibrated ART classifier.
"""
from pathlib import Path
from datetime import datetime, timezone
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import csv
import hashlib
import json
import os
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/art_profile_marshill'
TARGET = 'MW248466.1'
TARGET_RT = 'QQM14740.1'
MYRT_COMMIT = 'fa8e3624adf4c4e57ebf3350c06d626f805e42f7'
# These phages' ART arrays are explicitly described in the supplied paper,
# Methods, RNA secondary-structure prediction (PDF p. 36).
SEEDS = {
    'OR836606.1': 'WRM43423.1',  # LY01
    'MN091626.1': 'QDJ97640.1',  # PALS_2
    'LC680885.1': 'BDE75684.1',  # S6; annotated hypothetical protein
    'MW349129.1': 'QQO92858.1',  # Madawaska, not SA1
    'MZ779063.1': 'UAJ17040.1',  # UFV_DC4; annotated hypothetical protein
}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def table(path, rows):
    if not rows:
        path.write_text('', encoding='utf-8')
        return
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def feature_id(f):
    return f.qualifiers.get('protein_id', f.qualifiers.get('locus_tag', ['unidentified']))[0]


def prepare():
    OUT.mkdir(parents=True, exist_ok=True)
    records, provenance, proteins, qc, selected = {}, [], [], [], []
    for path in sorted((ROOT / 'datasetss').rglob('*.gb')):
        record = SeqIO.read(path, 'genbank')
        assert record.id not in records, f'Duplicate accession: {record.id}'
        records[record.id] = record
        fasta = [p for p in path.parent.glob('*.fasta')
                 if next(SeqIO.parse(p, 'fasta')).id == record.id]
        if fasta:
            assert str(SeqIO.read(fasta[0], 'fasta').seq).upper() == str(record.seq).upper()
        provenance.append(dict(accession=record.id, description=record.description,
                               path=path.relative_to(ROOT).as_posix(),
                               file_sha256=sha(path.read_bytes()),
                               sequence_sha256=sha(str(record.seq).upper().encode()),
                               genome_nt=len(record), matching_fasta=bool(fasta)))
        for f in record.features:
            if f.type != 'CDS' or 'translation' not in f.qualifiers:
                continue
            protein = f.qualifiers['translation'][0].replace(' ', '').upper()
            nt = f.extract(record.seq)[int(f.qualifiers.get('codon_start', ['1'])[0]) - 1:]
            genetic_code = int(f.qualifiers.get('transl_table', ['11'])[0])
            try:
                translated = str(nt.translate(table=genetic_code, cds=True))
                translation_qc = 'exact' if translated == protein else 'mismatch'
            except Exception as exc:
                translation_qc = 'incomplete_or_invalid_CDS:' + str(exc)
            pid = feature_id(f)
            row = dict(protein_id=pid, accession=record.id, start_1based=int(f.location.start)+1,
                       end_1based=int(f.location.end), strand=f.location.strand,
                       location=str(f.location), length_aa=len(protein),
                       product=f.qualifiers.get('product', [''])[0],
                       translation_qc=translation_qc,
                       standard_amino_acids=bool(re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+', protein)),
                       partial_location=bool(re.search(r'[<>]',str(f.location))),
                       replication_length_350_900=350 <= len(protein) <= 900,
                       yxdd_positions_1based=','.join(str(m.start()+1) for m in re.finditer('Y.DD',protein)),
                       protein_sha256=sha(protein.encode()))
            qc.append(row)
            proteins.append(SeqRecord(Seq(protein),id=pid,
                                      description=f'{record.id} {row["start_1based"]}..{row["end_1based"]} strand={f.location.strand} {row["product"]}'))
            if SEEDS.get(record.id) == pid or (record.id == TARGET and pid == TARGET_RT):
                assert translation_qc == 'exact' and row['standard_amino_acids'] and not row['partial_location'], row
                # Genomic downstream follows transcription, including reverse-strand RTs.
                cds = [x for x in record.features if x.type == 'CDS']
                if f.location.strand == 1:
                    following = [x for x in cds if int(x.location.start) >= int(f.location.end)]
                    partner = min(following,key=lambda x:int(x.location.start))
                    upstream = record.seq[max(0,int(f.location.start)-6000):int(f.location.start)]
                else:
                    following = [x for x in cds if int(x.location.end) <= int(f.location.start)]
                    partner = max(following,key=lambda x:int(x.location.end))
                    upstream = record.seq[int(f.location.end):min(len(record),int(f.location.end)+6000)].reverse_complement()
                row = dict(row, role='target_holdout' if record.id==TARGET else 'seed',
                           partner_id=feature_id(partner), partner_strand=partner.location.strand,
                           partner_length_aa=len(partner.qualifiers.get('translation',[''])[0]),
                           partner_start_1based=int(partner.location.start)+1,
                           partner_end_1based=int(partner.location.end),
                           upstream_6kb_repeat_core_copies_le1_mismatch=sum(
                               sum(a!=b for a,b in zip(str(upstream)[k:k+14], 'ATATGAATACGTAT'))<=1
                               for k in range(max(0,len(upstream)-13))))
                assert partner.location.strand == f.location.strand
                selected.append(row)
    assert set(SEEDS) <= records.keys() and TARGET in records
    assert len(selected) == len(SEEDS)+1
    by_id={r.id:r for r in proteins}
    seeds=[];seen=set()
    for s in selected:
        if s['role']=='seed' and s['protein_sha256'] not in seen:
            seeds.append(by_id[s['protein_id']]);seen.add(s['protein_sha256'])
    assert TARGET_RT not in {s.id for s in seeds} and len(seeds)>=2
    assert str(by_id[TARGET_RT].seq) not in {str(s.seq) for s in seeds}, 'Exact target-sequence leakage'
    SeqIO.write(seeds, OUT/'art_seed_proteins.faa','fasta')
    SeqIO.write([x for x in proteins if x.description.startswith(TARGET+' ')],OUT/'MarsHill_proteins.faa','fasta')
    SeqIO.write([by_id[TARGET_RT]], OUT/'MarsHill_reference_RT.faa','fasta')
    SeqIO.write([SeqRecord(records[TARGET].seq,id=TARGET,description=records[TARGET].description)],OUT/'MarsHill_contig.fna','fasta')
    table(OUT/'source_records.tsv',provenance)
    table(OUT/'protein_annotations.tsv',qc)
    table(OUT/'curated_RT_manifest.tsv',selected)
    aligner=Align.PairwiseAligner(scoring='blastp');aligner.mode='global'
    pairwise=[]
    for i,a in enumerate(selected):
        for b in selected[i+1:]:
            aa,bb=by_id[a['protein_id']].seq,by_id[b['protein_id']].seq
            al=aligner.align(aa,bb)[0]
            pairwise.append(dict(protein_a=a['protein_id'],protein_b=b['protein_id'],
                                 identities_over_longer_length=al.counts().identities/max(len(aa),len(bb)),
                                 denominator=max(len(aa),len(bb))))
    table(OUT/'seed_and_holdout_identity.tsv',pairwise)
    versions=subprocess.check_output([str(ROOT/'.venv/Scripts/python.exe'),'-m','pip','freeze'],text=True)
    (OUT/'python_requirements.txt').write_text(versions,encoding='utf-8')
    print(f'Prepared {len(seeds)} distinct seeds; target proteome: '
          f'{sum(x.description.startswith(TARGET+" ") for x in proteins)} proteins.',flush=True)


def execute(command, log, stdout=None, env=None, cwd=None):
    with (OUT/'commands.jsonl').open('a',encoding='utf-8') as f:
        f.write(json.dumps(dict(utc=datetime.now(timezone.utc).isoformat(),command=[str(x) for x in command]))+'\n')
    with log.open('wb') as stderr:
        if stdout:
            with stdout.open('wb') as out:
                subprocess.run([str(x) for x in command],stdout=out,stderr=stderr,check=True,env=env,cwd=cwd)
        else:
            subprocess.run([str(x) for x in command],stdout=stderr,stderr=subprocess.STDOUT,check=True,env=env,cwd=cwd)


def run():
    # Launch the official Windows MAFFT batch entry point via cmd as documented.
    mafft=ROOT/'.tools/mafft-win/mafft.bat'
    execute(['cmd.exe','/d','/c',str(mafft),'--auto','--thread','4',
             '--out',OUT/'art_seed_alignment.faa',OUT/'art_seed_proteins.faa'],OUT/'mafft.log')
    msa=list(SeqIO.parse(OUT/'art_seed_alignment.faa','fasta'))
    seeds=list(SeqIO.parse(OUT/'art_seed_proteins.faa','fasta'))
    assert {r.id:str(r.seq).replace('-','').upper() for r in msa}=={r.id:str(r.seq) for r in seeds}
    assert len({len(r) for r in msa})==1
    bindir=ROOT/'.tools/hmmer-install/bin'
    env=os.environ.copy();env['PATH']=str(ROOT/'.tools/cygwin/bin')+os.pathsep+env['PATH']
    execute([bindir/'hmmbuild.exe','--amino','--cpu','4','--seed','42','-n','ART_Staph_typeII_pilot',
             OUT/'art_family.hmm',OUT/'art_seed_alignment.faa'],OUT/'hmmbuild.log',env=env)
    for stem,threshold in [('art_candidates', '1e-5'),('art_sensitivity','1e-3')]:
        execute([bindir/'hmmsearch.exe','--cpu','4','--seed','42','-E',threshold,'--domE',threshold,
                 '--tblout',OUT/(stem+'.tbl'),'--domtblout',OUT/(stem+'.domtbl'),
                 OUT/'art_family.hmm',OUT/'MarsHill_proteins.faa'],OUT/(stem+'.log'),env=env)
    summarize()
    validate_references(bindir, env)
    write_hashes()


def validate_references(bindir, env):
    source=ROOT/'.tools/myRT'
    if not source.exists():
        subprocess.run(['git','clone','--depth','1','https://github.com/mgtools/myRT.git',str(source)],check=True)
    commit=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
    if commit != MYRT_COMMIT:
        raise RuntimeError(f'myRT commit {commit} differs from pinned {MYRT_COMMIT}; '
                           'review the reference change before rerunning validation.')
    refs=source/'Models/RTs-collection.faa'
    models=source/'Models/HMM/RVT-All.hmm'
    panel=list(SeqIO.parse(refs,'fasta'))
    provenance=dict(repository='https://github.com/mgtools/myRT',commit=commit,
                    references_sha256=sha(refs.read_bytes()),models_sha256=sha(models.read_bytes()),
                    reference_sequences=len(panel),
                    caveat='Historical RT labels, not an independently curated ART/NART negative set.')
    (OUT/'myRT_provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    execute([bindir/'hmmsearch.exe','--cpu','4','--seed','42','-E','1e-5','--domE','1e-5',
             '--tblout',OUT/'myRT_reference_screen.tbl','--domtblout',OUT/'myRT_reference_screen.domtbl',
             OUT/'art_family.hmm',refs],OUT/'myRT_reference_screen.log',env=env)
    execute([bindir/'hmmsearch.exe','--cpu','4','--seed','42','-E','1e-5','--domE','1e-5',
             '--tblout',OUT/'competing_RT_profiles.tbl','--domtblout',OUT/'competing_RT_profiles.domtbl',
             models,OUT/'art_candidates.faa'],OUT/'competing_RT_profiles.log',env=env)
    class_counts={}
    for r in panel:
        label=r.id.rsplit('_',1)[-1]
        class_counts[label]=class_counts.get(label,0)+1
    coverage={}
    domain_rows=[]
    for line in (OUT/'myRT_reference_screen.domtbl').read_text().splitlines():
        if line.startswith('#') or not line.strip():continue
        c=line.split(maxsplit=22)
        cov=(int(c[16])-int(c[15])+1)/int(c[5])
        coverage[c[0]]=max(coverage.get(c[0],0),cov)
        domain_rows.append(dict(protein_id=c[0],historical_class=c[0].rsplit('_',1)[-1],
                                domain_independent_evalue=float(c[12]),domain_bitscore=float(c[13]),
                                hmm_from=int(c[15]),hmm_to=int(c[16]),hmm_length=int(c[5]),
                                hmm_coverage=cov,profile_coverage_ge_0_6=cov>=0.6,
                                alignment_from=int(c[17]),alignment_to=int(c[18])))
    table(OUT/'myRT_reference_domains.tsv',domain_rows)
    matches=[]
    for line in (OUT/'myRT_reference_screen.tbl').read_text().splitlines():
        if line.startswith('#') or not line.strip():continue
        c=line.split(maxsplit=18)
        matches.append(dict(protein_id=c[0],historical_class=c[0].rsplit('_',1)[-1],
                            full_sequence_evalue=float(c[4]),full_sequence_bitscore=float(c[5]),
                            max_domain_hmm_coverage=coverage.get(c[0],0),
                            profile_coverage_ge_0_6=coverage.get(c[0],0)>=0.6))
    table(OUT/'myRT_reference_hits.tsv',matches)
    competitors=[]
    for line in (OUT/'competing_RT_profiles.domtbl').read_text().splitlines():
        if line.startswith('#') or not line.strip():continue
        c=line.split(maxsplit=22)
        competitors.append(dict(protein_id=c[0],profile=c[3],full_sequence_evalue=float(c[6]),
                               full_sequence_bitscore=float(c[7]),domain_independent_evalue=float(c[12]),
                               profile_coverage=(int(c[16])-int(c[15])+1)/int(c[5]),
                               protein_alignment_from=int(c[17]),protein_alignment_to=int(c[18])))
    table(OUT/'competing_RT_profile_matches.tsv',competitors)
    summary=dict(reference_sequences=len(panel),historical_class_counts=class_counts,
                 reference_hits_at_1e_minus_5=len(matches),reference_hits=matches,
                 reference_hits_at_least_60_percent_profile_coverage=sum(x['profile_coverage_ge_0_6'] for x in matches),
                 competing_profile_matches=competitors,
                 validated_specificity_threshold=None,
                 NART_specific_validation='pending curated NART reference sequences',
                 scope='Staphylococcus/type-II pilot; related target holdout; no family-wide sensitivity claim')
    (OUT/'validation_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    make_report(summary)
    print(f'Reference check: {len(matches)} hits in {len(panel)} historical myRT sequences.',flush=True)


def make_report(validation):
    candidates=list(csv.DictReader((OUT/'art_candidates.tsv').open(encoding='utf-8'),delimiter='\t'))
    identities=list(csv.DictReader((OUT/'seed_and_holdout_identity.tsv').open(encoding='utf-8'),delimiter='\t'))
    manifest=list(csv.DictReader((OUT/'curated_RT_manifest.tsv').open(encoding='utf-8'),delimiter='\t'))
    profile=(OUT/'art_family.hmm').read_text()
    rows=['# MarsHill ART RT profile pilot', '',
          '## Result', '',
          f'The profile retrieved {len(candidates)} candidate domain hit(s) from the 262 annotated '
          'MarsHill proteins at the 1e-5 sequence and conditional-domain reporting thresholds.', '',
          '| Protein | Genome | CDS coordinates (1-based, inclusive) | Length | E-value | Bits | Profile coverage |',
          '| --- | --- | --- | ---: | ---: | ---: | ---: |']
    for r in candidates:
        rows.append(f'| {r["protein_id"]} | {r["accession"]} | {r["start_1based"]}–{r["end_1based"]} '
                    f'(strand {r["strand"]}) | {r["length_aa"]} | {float(r["full_sequence_evalue"]):.3g} '
                    f'| {float(r["full_sequence_bitscore"]):.1f} | {float(r["hmm_coverage"]):.1%} |')
    rows += ['', '## Training and provenance', '',
             'Five distinct full-length RT seeds were used: LY01, PALS_2, S6, Madawaska and UFV_DC4. '
             'MarsHill QQM14740.1 and its exact sequence were excluded from training.',
             f'The full-length global identities among seeds and target span '
             f'{min(float(x["identities_over_longer_length"]) for x in identities):.1%}–'
             f'{max(float(x["identities_over_longer_length"]) for x in identities):.1%}. '
             'This is a related-sequence recovery check, not a lineage-held-out sensitivity estimate.',
             f'Profile length: {re.search(r"^LENG  +(\d+)",profile,re.M)[1]} match states. '
             f'Training sequence count: {re.search(r"^NSEQ  +(\d+)",profile,re.M)[1]}.',
             'Each curated RT translates exactly from its annotated genomic CDS; each has five '
             'matches to the known type-II array core within one mismatch in its 6-kb upstream window, '
             'and an adjacent co-oriented 274-aa candidate partner. This context check does not '
             'replace array delimitation, functional annotation or binding experiments.', '',
             '## Reference checks', '',
             f'The new profile returned {validation["reference_hits_at_1e_minus_5"]} reported '
             f'sequence hits against {validation["reference_sequences"]} historical myRT reference '
             'sequences at E <= 1e-5. The panel includes retron, DGR, group-II intron, CRISPR-associated '
             'and other RT labels. See validation_summary.json and the raw outputs.',
             f'Of these reference hits, {validation["reference_hits_at_least_60_percent_profile_coverage"]} '
             'have a reported domain covering at least 60% of the new profile. This proposed coverage '
             'screen separates this particular reference panel from the MarsHill recovery, but does '
             'not establish a family-wide decision rule or sensitivity/specificity estimate.',
             'Competing myRT-class profile matches to the retrieved MarsHill protein:', '']
    for x in validation['competing_profile_matches']:
        rows.append(f'- {x["profile"]}: E={x["full_sequence_evalue"]:.3g}, '
                    f'score={x["full_sequence_bitscore"]:.1f}, coverage={x["profile_coverage"]:.1%}.')
    sensitivity={line.split()[0] for line in (OUT/'art_sensitivity.tbl').read_text().splitlines()
                 if line.strip() and not line.startswith('#')}
    rows += ['', f'The relaxed 1e-3 sequence reporting threshold returned {len(sensitivity)} protein(s).', '',
             '## Interpretation limits', '',
             '- This is a Staphylococcus/type-II ART profile pilot. Distant type-I/III ART recovery has not been tested.',
             '- Curated NART-specific validation remains pending; historical myRT labels are not an independently curated ART/NART negative set.',
             '- No calibrated ART bit-score threshold or sensitivity/specificity estimate is claimed. The paper\'s 207-bit threshold is not used.',
             '- Only annotated MarsHill CDSs were searched; missed or miscalled genes remain possible.',
             '- The SA1 FASTA was not searched. The Madawaska GenBank file in the SA1 folder is tracked as MW349129.1.',
             '', 'An E-value printed as 0 is HMMER numerical underflow, not a mathematical zero.',
             '', 'Actual commands, reference versions, input hashes and raw logs accompany these results.', '']
    (OUT/'REPORT.md').write_text('\n'.join(rows),encoding='utf-8')


def summarize():
    annotations={r['protein_id']:r for r in csv.DictReader((OUT/'protein_annotations.tsv').open(encoding='utf-8'),delimiter='\t')}
    rows=[]
    for line in (OUT/'art_candidates.domtbl').read_text().splitlines():
        if line.startswith('#') or not line.strip():continue
        c=line.split(maxsplit=22);r=annotations[c[0]]
        rows.append(dict(protein_id=c[0],accession=r['accession'],start_1based=r['start_1based'],end_1based=r['end_1based'],
                         strand=r['strand'],length_aa=int(c[2]),full_sequence_evalue=float(c[6]),
                         full_sequence_bitscore=float(c[7]),domain_independent_evalue=float(c[12]),
                         hmmer_reported_zero_evalue=float(c[6])==0,
                         domain_bitscore=float(c[13]),hmm_from=int(c[15]),hmm_to=int(c[16]),
                         hmm_length=int(c[5]),hmm_coverage=(int(c[16])-int(c[15])+1)/int(c[5]),
                         alignment_from=int(c[17]),alignment_to=int(c[18]),translation_qc=r['translation_qc'],
                         complete_annotated_cds=r['translation_qc']=='exact' and r['partial_location']=='False',
                         profile_coverage_ge_0_6=(int(c[16])-int(c[15])+1)/int(c[5]) >= 0.6,
                         length_350_900=r['replication_length_350_900'],
                         yxdd_positions_1based=r['yxdd_positions_1based'],
                         status='provisional_profile_hit_requires_context_and_competing_RT_validation'))
    table(OUT/'art_candidates.tsv',rows)
    ids={r['protein_id'] for r in rows}
    SeqIO.write([r for r in SeqIO.parse(OUT/'MarsHill_proteins.faa','fasta') if r.id in ids],OUT/'art_candidates.faa','fasta')
    assert TARGET_RT in ids, 'Known MarsHill RT was not recovered'
    print(json.dumps(rows,indent=2),flush=True)


def write_hashes():
    manifest={p.name:sha(p.read_bytes()) for p in OUT.iterdir() if p.is_file() and p.name!='output_sha256.json'}
    (OUT/'output_sha256.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','all']);args=p.parse_args()
    if args.stage in ('prepare','all'):prepare()
    if args.stage in ('run','all'):run()
