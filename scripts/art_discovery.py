"""Independent RT-first and partner-first profile searches with accession-linked context.

Missing seed families are explicitly unavailable, never replaced with generic
domain models. Coordinates in TSV outputs are 1-based inclusive; internal
coordinates are 0-based half-open. Run with the project's .venv Python.
"""
from pathlib import Path
from collections import Counter, defaultdict
from datetime import datetime, timezone
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import csv
import hashlib
import json
import math
import os
import re
import subprocess
import shutil
import sys
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
FAMILIES = ('ART_RT', 'Type_I_partner', 'Type_II_partner',
            'Type_III_first_partner', 'Type_III_second_partner')


def digest(data):
    return hashlib.sha256(data).hexdigest()


def tsv(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or (list(rows[0]) if rows else ['status'])
    with path.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fields, delimiter='\t', extrasaction='raise')
        writer.writeheader()
        writer.writerows(rows)


def fetch(url, path):
    if not path.exists():
        request = urllib.request.Request(url, headers={'User-Agent': 'ART-public-locus-study/1.0'})
        with urllib.request.urlopen(request, timeout=60) as r:
            data = r.read()
        path.write_bytes(data)
    return dict(path=path.relative_to(ROOT).as_posix(), url=url,
                sha256=digest(path.read_bytes()),
                retrieval_recorded_utc=datetime.now(timezone.utc).isoformat())


def protein_id(feature):
    return feature.qualifiers.get('protein_id', feature.qualifiers.get('locus_tag', ['']))[0]


def neighbors(anchor, genes, direction, count=3):
    """Order all CDSs along the anchor's transcription axis, including overlaps.

    On a + anchor the ordering key is genomic start; on a - anchor it is
    negative genomic end. Other-strand neighbors remain in the ranking.
    """
    if anchor['strand'] not in (-1,1):
        return []
    axis = lambda g: g['start0'] if anchor['strand'] == 1 else -g['end0']
    pos = axis(anchor)
    candidates = [g for g in genes if g['gene_key'] != anchor['gene_key']
                  and (axis(g) < pos if direction == 'upstream' else axis(g) > pos)]
    return sorted(candidates, key=lambda g: (axis(g), g['gene_key']),
                  reverse=direction == 'upstream')[:count]


def window_bounds(gene, contig_length, flank):
    return max(0, gene['start0']-flank), min(contig_length, gene['end0']+flank)


def repeat_hint(record, rt, genes, flank=6000):
    """Conservative exact-word context hint, not a validated ART array caller.

    Scan the oriented upstream DNA, masking annotated CDSs >=300 nt. Require
    an exact 10-mer with >=3 copies at roughly regular 100..450-nt spacing,
    allowing a skipped copy. This is not a claim about RNA processing.
    No shuffling-based significance or biological absence claim is made.
    """
    if rt['strand'] not in (-1,1):
        return dict(upstream_nt_available=0,repeat_hint_copies=0,repeat_hint_word='',
                    repeat_hint_positions_1based='',repeat_hint_status='unassessable_unknown_RT_orientation')
    if rt['strand'] == 1:
        start, end = max(0, rt['start0']-flank), rt['start0']
        sequence = str(record.seq[start:end]).upper()
    else:
        start, end = rt['end0'], min(len(record), rt['end0']+flank)
        sequence = str(record.seq[start:end].reverse_complement()).upper()
    mask = [False]*len(sequence)
    for g in genes:
        if g['end0']-g['start0'] < 300:
            continue
        a, b = max(start, g['start0']), min(end, g['end0'])
        if a >= b:
            continue
        lo, hi = (a-start, b-start) if rt['strand'] == 1 else (end-b, end-a)
        mask[lo:hi] = [True]*(hi-lo)
    words = defaultdict(list)
    for i in range(max(0, len(sequence)-9)):
        word = sequence[i:i+10]
        counts = Counter(word)
        entropy = -sum((n/10)*math.log2(n/10) for n in counts.values())
        if (any(mask[i:i+10]) or not set(word) <= set('ACGT') or len(counts)<3
                or entropy<1.2 or re.search(r'(.)\1{5}', word)):
            continue
        words[word].append(i)
    best = ('', [])
    for word, positions in words.items():
        if len(positions)<3 or len(positions)>80:
            continue
        for i, first in enumerate(positions):
            for second in positions[i+1:]:
                step = second-first
                if step<100 or step>450:
                    continue
                chain = [first, second]
                for nxt in positions:
                    if nxt<=second:
                        continue
                    gap = nxt-chain[-1]
                    if gap>=100 and any(abs(gap-m*step)<=0.3*m*step for m in (1,2)):
                        chain.append(nxt)
                if len(chain)>len(best[1]):
                    best = word, chain
    word, chain = best
    if len(chain)<3:
        word, chain = '', []
    absolute = [start+i+1 if rt['strand']==1 else end-i-9 for i in chain]
    return dict(upstream_nt_available=len(sequence), repeat_hint_copies=len(chain),
                repeat_hint_word=word, repeat_hint_positions_1based=','.join(map(str, absolute)),
                repeat_hint_status='exact_word_hint_not_validated_array' if chain else
                ('no_exact_word_hint_unassessable_short_flank' if len(sequence)<1500 else
                 'no_exact_word_hint_not_evidence_of_absence'))


class Pipeline:
    def __init__(self, config_path):
        self.cfg = json.loads(config_path.read_text(encoding='utf-8'))
        assert tuple(self.cfg['families']) == FAMILIES
        self.out = ROOT/self.cfg['output']
        self.out.mkdir(parents=True, exist_ok=True)
        for sub in ('profiles','searches','windows','reference_sources','calibration'):
            (self.out/sub).mkdir(exist_ok=True)
        self.env = os.environ.copy()
        self.env['PATH'] = str(ROOT/'.tools/cygwin/bin')+os.pathsep+self.env['PATH']
        self.bindir = ROOT/'.tools/hmmer-install/bin'
        self.records = {}
        self.genes = {}
        self.by_contig = defaultdict(list)
        self.proteins = {}
        self.registry = []
        self.hits = []
        self.accepted = defaultdict(dict)
        self.raw = defaultdict(dict)
        self.generic = {}
        (self.out/'run_config.json').write_text(json.dumps(self.cfg, indent=2), encoding='utf-8')
        (self.out/'code').mkdir(exist_ok=True)
        shutil.copy2(Path(__file__), self.out/'code/art_discovery.py')
        test_source=Path(__file__).with_name('test_art_discovery.py')
        if test_source.exists():shutil.copy2(test_source,self.out/'code/test_art_discovery.py')

    def command(self, args, log):
        args = [str(x) for x in args]
        with (self.out/'commands.jsonl').open('a', encoding='utf-8') as f:
            f.write(json.dumps(dict(utc=datetime.now(timezone.utc).isoformat(), args=args))+'\n')
        with log.open('wb') as f:
            subprocess.run(args, env=self.env, stdout=f, stderr=subprocess.STDOUT,
                           check=True, cwd=ROOT)

    def prepare(self):
        ref = self.out/'reference_sources'
        downloads = []
        downloads.append(fetch('https://www.ebi.ac.uk/ena/browser/api/embl/MZ422438.1?download=false',
                               ref/'LPJP1_MZ422438.1.gb'))
        downloads.append(fetch('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=MW218148.1&rettype=gbwithparts&retmode=text',
                               ref/'SA1_MW218148.1.gb'))
        query='(gene_exact:rimI OR gene_exact:rimJ OR gene_exact:rimL) AND organism_id:83333 AND reviewed:true'
        url='https://rest.uniprot.org/uniprotkb/search?'+urllib.parse.urlencode({'query':query,'format':'fasta'})
        downloads.append(fetch(url, ref/'unrelated_GNATs.faa'))
        (self.out/'external_source_provenance.json').write_text(json.dumps(downloads,indent=2),encoding='utf-8')
        inputs=[(p,'genbank','local_target') for p in sorted((ROOT/'datasetss').rglob('*.gb'))]
        inputs += [(ref/'SA1_MW218148.1.gb','genbank','public_annotation_for_local_SA1_FASTA'),
                   (ref/'LPJP1_MZ422438.1.gb','embl','public_reference_also_screened')]
        sources=[]
        for path,fmt,role in inputs:
            record=SeqIO.read(path,fmt)
            if record.id in self.records:
                raise ValueError(f'Duplicate accession: {record.id}')
            self.records[record.id]=record
            for fasta in (ROOT/'datasetss').rglob('*.fasta'):
                for f in SeqIO.parse(fasta,'fasta'):
                    if f.id==record.id:
                        assert str(f.seq).upper()==str(record.seq).upper(), f'FASTA mismatch: {f.id}'
            sources.append(dict(accession=record.id,path=path.relative_to(ROOT).as_posix(),
                                format=fmt,data_role=role,length_nt=len(record),
                                file_sha256=digest(path.read_bytes()),
                                sequence_sha256=digest(str(record.seq).upper().encode())))
            for index,feature in enumerate(record.features):
                if feature.type!='CDS':continue
                pid=protein_id(feature) or f'{record.id}_CDS_{index}'
                key=f'{record.id}|{pid}'
                if pid in self.proteins:raise ValueError(f'Nonunique protein accession: {pid}')
                sequence=feature.qualifiers.get('translation',[''])[0].upper()
                qc='no_annotated_translation'
                if sequence:
                    nt=feature.extract(record.seq)[int(feature.qualifiers.get('codon_start',['1'])[0])-1:]
                    try:
                        translated=str(nt.translate(table=int(feature.qualifiers.get('transl_table',['11'])[0]),cds=True))
                        qc='exact' if translated==sequence else 'translation_mismatch'
                    except Exception as exc:qc='partial_or_invalid:'+str(exc)
                    self.proteins[pid]=SeqRecord(Seq(sequence),id=pid,description=record.id)
                gene=dict(gene_key=key,protein_id=pid,accession=record.id,
                          start0=int(feature.location.start),end0=int(feature.location.end),
                          strand=int(feature.location.strand or 0),length_aa=len(sequence),
                          product=feature.qualifiers.get('product',[''])[0],translation_qc=qc,
                          partial_location=bool(re.search('[<>]',str(feature.location))),
                          source_path=path.relative_to(ROOT).as_posix())
                self.genes[pid]=gene;self.by_contig[record.id].append(gene)
        self.sources=sources
        SeqIO.write(self.proteins.values(),self.out/'target_proteins.faa','fasta')
        SeqIO.write(self.records.values(),self.out/'target_contigs.fna','fasta')
        tsv(self.out/'source_records.tsv',sources)
        tsv(self.out/'genes.tsv',[dict(g, start_1based=g['start0']+1,end_1based=g['end0']) for g in self.genes.values()])
        self.train={};self.heldout={};seed_manifest=[]
        for family,cfg in self.cfg['families'].items():
            selected=[]
            for pid in cfg['seed_ids']:
                g=self.genes[pid]
                assert g['translation_qc']=='exact' and not g['partial_location'], (family,pid,g)
                selected.append(self.proteins[pid])
            if cfg.get('seed_fasta'):
                selected += list(SeqIO.parse(ROOT/cfg['seed_fasta'],'fasta'))
            unique={}
            for r in selected:
                assert re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+',str(r.seq)),r.id
                unique.setdefault(str(r.seq),r)
                seed_manifest.append(dict(family=family,protein_id=r.id,length_aa=len(r),
                                          sequence_sha256=digest(str(r.seq).encode()),curation=cfg['curation']))
            self.train[family]=list(unique.values())
            positive=[self.proteins[pid] for pid in cfg['heldout_ids']]
            if cfg.get('heldout_fasta'):positive+=list(SeqIO.parse(ROOT/cfg['heldout_fasta'],'fasta'))
            self.heldout[family]=positive
            assert not ({str(r.seq) for r in positive}&set(unique)), f'Exact-sequence holdout leakage: {family}'
            SeqIO.write(self.train[family],self.out/'profiles'/f'{family}.seeds.faa','fasta')
        tsv(self.out/'seed_manifest.tsv',seed_manifest)
        self.make_calibration()
        print(f'Prepared {len(self.records)} genomes, {len(self.proteins)} annotated proteins.',flush=True)

    def make_calibration(self):
        controls=list(SeqIO.parse(self.out/'reference_sources/unrelated_GNATs.faa','fasta'))
        assert len(controls)==3,'Expected RimI, RimJ and RimL GNAT references'
        annotations=[]
        for r in controls:
            annotations.append(dict(protein_id=r.id,role='reviewed_unrelated_GNAT',family='',description=r.description))
        used=set()
        for family,rs in self.heldout.items():
            for r in rs:
                if r.id not in used:controls.append(r);used.add(r.id)
                annotations.append(dict(protein_id=r.id,role='related_sequence_positive_holdout',family=family,description=r.description))
        alltrain={r.id for rs in self.train.values() for r in rs}
        background=[]
        for pid,g in sorted(self.genes.items()):
            if pid in alltrain or pid in used or pid not in self.proteins:continue
            if re.search(r'capsid|terminase|helicase|ligase|RNA polymerase|ribonuclease',g['product'],re.I):
                background.append(self.proteins[pid])
        for r in background[:50]:
            controls.append(r)
            annotations.append(dict(protein_id=r.id,role='annotated_phage_background_not_verified_nonbinder',family='',description=self.genes[r.id]['product']))
        SeqIO.write(controls,self.out/'calibration/controls.faa','fasta')
        tsv(self.out/'calibration/control_manifest.tsv',annotations)
        self.control_annotations=annotations

    def search(self, hmm, database, stem):
        self.command([self.bindir/'hmmsearch.exe','--cpu',self.cfg['threads'],'--seed','42',
                      '-Z',len(self.proteins),'--domZ',len(self.proteins),
                      '-E',self.cfg['retrieval_evalue'],'--domE',self.cfg['retrieval_evalue'],
                      '--tblout',str(stem)+'.tbl','--domtblout',str(stem)+'.domtbl',hmm,database],
                     Path(str(stem)+'.log'))
        domains=[]
        for line in Path(str(stem)+'.domtbl').read_text().splitlines():
            if line.startswith('#') or not line.strip():continue
            c=line.split(maxsplit=22)
            row=dict(protein_id=c[0],family=c[3],length_aa=int(c[2]),
                     sequence_evalue=float(c[6]),sequence_bitscore=float(c[7]),
                     domain_independent_evalue=float(c[12]),domain_bitscore=float(c[13]),
                     profile_length=int(c[5]),hmm_from=int(c[15]),hmm_to=int(c[16]),
                     protein_from=int(c[17]),protein_to=int(c[18]),
                     profile_coverage=(int(c[16])-int(c[15])+1)/int(c[5]))
            row['passes_screen']=(row['sequence_evalue']<=self.cfg['sequence_evalue'] and
                                  row['domain_independent_evalue']<=self.cfg['domain_independent_evalue'] and
                                  row['profile_coverage']>=self.cfg['minimum_profile_coverage'])
            domains.append(row)
        return domains

    def build_and_search(self):
        calibration=[]
        for family in FAMILIES:
            seeds=self.train[family];cfg=self.cfg['families'][family]
            status=dict(family=family,unique_seed_sequences=len(seeds),seed_loci=len(cfg['seed_ids']),
                        profile_path='',profile_sha256='',heldout_positive_count=len(self.heldout[family]),
                        heldout_positive_recovered=0,reviewed_GNAT_hits=0,annotated_background_hits=0,
                        status='',calibration_status='',curation=cfg['curation'])
            if not seeds:
                status.update(status='unavailable_no_curated_seeds',calibration_status='not_testable')
                self.registry.append(status);continue
            stem=self.out/'profiles'/family
            msa=Path(str(stem)+'.alignment.faa');hmm=Path(str(stem)+'.hmm')
            self.command(['cmd.exe','/d','/c',ROOT/'.tools/mafft-win/mafft.bat','--auto','--thread',
                          self.cfg['threads'],'--out',msa,Path(str(stem)+'.seeds.faa')],Path(str(stem)+'.mafft.log'))
            aligned=list(SeqIO.parse(msa,'fasta'))
            assert {r.id:str(r.seq).replace('-','').upper() for r in aligned}=={r.id:str(r.seq) for r in seeds}
            self.command([self.bindir/'hmmbuild.exe','--amino','--cpu',self.cfg['threads'],
                          '--seed','42','-n',family,hmm,msa],Path(str(stem)+'.hmmbuild.log'))
            hits=self.search(hmm,self.out/'target_proteins.faa',self.out/'searches'/family)
            self.hits+=hits
            for hit in hits:
                pid=hit['protein_id']
                if pid not in self.raw[family] or hit['domain_bitscore']>self.raw[family][pid]['domain_bitscore']:
                    self.raw[family][pid]=hit
                if hit['passes_screen'] and (pid not in self.accepted[family] or hit['domain_bitscore']>self.accepted[family][pid]['domain_bitscore']):
                    self.accepted[family][pid]=hit
            checks=self.search(hmm,self.out/'calibration/controls.faa',self.out/'calibration'/family)
            calibration+=checks
            passed={r['protein_id'] for r in checks if r['passes_screen']}
            positive={r.id for r in self.heldout[family]}
            gnat={r['protein_id'] for r in self.control_annotations if r['role']=='reviewed_unrelated_GNAT'}
            background={r['protein_id'] for r in self.control_annotations if r['role'].startswith('annotated_phage')}
            status.update(profile_path=hmm.relative_to(ROOT).as_posix(),profile_sha256=digest(hmm.read_bytes()),
                          status='provisional_retrieval_profile',heldout_positive_recovered=len(positive&passed),
                          reviewed_GNAT_hits=len(gnat&passed),annotated_background_hits=len(background&passed),
                          calibration_status='insufficient_positive_references_single_seed' if len(seeds)==1 else
                          'limited_related_holdout_check_not_lineage_validation')
            self.registry.append(status)
            print(f'{family}: {len(seeds)} seeds; {len(self.accepted[family])} screened protein hits.',flush=True)
        tsv(self.out/'profile_registry.tsv',self.registry)
        tsv(self.out/'protein_hits.tsv',self.hits)
        tsv(self.out/'calibration/domain_results.tsv',calibration)
        generic=self.search(ROOT/'.tools/myRT/Models/HMM/RVT-All.hmm',self.out/'target_proteins.faa',
                            self.out/'searches/generic_RT_context_only')
        for h in generic:
            if (h['sequence_evalue']<=1e-5 and h['domain_independent_evalue']<=1e-5 and
                    h['profile_coverage']>=0.6):
                pid=h['protein_id']
                if pid not in self.generic or h['domain_bitscore']>self.generic[pid]['domain_bitscore']:self.generic[pid]=h
        tsv(self.out/'generic_RT_context_matches.tsv',generic)

    def context(self, anchor, direction):
        record=self.records[anchor['accession']]
        a,b=window_bounds(anchor,len(record),self.cfg['flank_nt'])
        rows=[]
        for rank,g in enumerate(neighbors(anchor,self.by_contig[anchor['accession']],direction),1):
            art=self.accepted['ART_RT'].get(g['protein_id'])
            weak=self.raw['ART_RT'].get(g['protein_id'])
            generic=self.generic.get(g['protein_id'])
            rows.append(dict(anchor_protein=anchor['protein_id'],accession=anchor['accession'],rank=rank,
                             protein_id=g['protein_id'],start_1based=g['start0']+1,end_1based=g['end0'],
                             strand=g['strand'],same_orientation=g['strand']==anchor['strand'],
                             overlap_nt=max(0,min(anchor['end0'],g['end0'])-max(anchor['start0'],g['start0'])),
                             fully_inside_6kb_window=g['start0']>=a and g['end0']<=b,
                             ART_profile_match=bool(art),weak_ART_profile_match=bool(weak),
                             ART_bitscore=weak['sequence_bitscore'] if weak else '',
                             ART_profile_coverage=weak['profile_coverage'] if weak else '',
                             partner_profile_matches=';'.join(f for f in FAMILIES[1:] if g['protein_id'] in self.accepted[f]),
                             generic_RT_profile=generic['family'] if generic else '',product=g['product']))
        return rows

    def assemble(self):
        loci={};upstream=[];downstream=[];windows={}
        missing=';'.join(r['family'] for r in self.registry if r['status'].startswith('unavailable'))
        training={r.id for rs in self.train.values() for r in rs}

        def locus(accession,rtid='',partnerid=''):
            key=f'{accession}|RT:{rtid}' if rtid else f'{accession}|partner:{partnerid}'
            if key not in loci:
                loci[key]=dict(accession=accession,rt_id=rtid,partners=defaultdict(set),
                               art_first=False,partner_first=False,anchors=set(),context_notes=set())
            return loci[key]

        for pid in self.accepted['ART_RT']:
            g=self.genes[pid];entry=locus(g['accession'],rtid=pid)
            entry['art_first']=True;entry['anchors'].add(pid)
            rows=self.context(g,'downstream');downstream+=rows
            # Partner associations are attached below using the partner's own
            # upstream ranking and extracted window. This prevents asymmetric
            # RT-window associations from adding an out-of-window partner.
            windows[pid]=g

        for family in FAMILIES[1:]:
            for pid in self.accepted[family]:
                g=self.genes[pid];rows=self.context(g,'upstream')
                for row in rows:row['anchor_partner_family']=family
                upstream+=rows
                strong=[r for r in rows if r['ART_profile_match'] and r['same_orientation'] and r['fully_inside_6kb_window']]
                generic=[r for r in rows if r['generic_RT_profile'] and r['same_orientation'] and r['fully_inside_6kb_window']]
                choices=strong if strong else generic
                if choices:
                    # Retain multiple compatible RTs as separate candidate associations,
                    # recording ambiguity instead of silently picking the closest.
                    for row in choices:
                        entry=locus(g['accession'],rtid=row['protein_id'])
                        entry['partner_first']=True;entry['partners'][family].add(pid);entry['anchors'].add(pid)
                        if not strong:entry['context_notes'].add('generic_RT_neighbor_without_screened_ART_match')
                        if len(choices)>1:entry['context_notes'].add('multiple_upstream_RT_candidates')
                else:
                    entry=locus(g['accession'],partnerid=pid)
                    entry['partner_first']=True;entry['partners'][family].add(pid);entry['anchors'].add(pid)
                    entry['context_notes'].add('no_cooriented_RT_match_in_first3_upstream_genes_inside_window')
                windows[pid]=g

        fasta=[];window_rows=[]
        for pid,g in sorted(windows.items()):
            record=self.records[g['accession']];a,b=window_bounds(g,len(record),self.cfg['flank_nt'])
            wid=f'{g["accession"]}__{pid}__{a+1}_{b}'
            fasta.append(SeqRecord(record.seq[a:b],id=wid,description=f'anchor={pid} strand={g["strand"]}; genomic plus strand'))
            path=self.out/'windows'/f'{wid}.fna'
            SeqIO.write([fasta[-1]],path,'fasta')
            window_rows.append(dict(anchor_protein=pid,accession=g['accession'],start_1based=a+1,end_1based=b,
                                    requested_flank_nt=self.cfg['flank_nt'],left_flank_nt=g['start0']-a,
                                    right_flank_nt=b-g['end0'],left_clipped=a==0 and g['start0']<self.cfg['flank_nt'],
                                    right_clipped=b==len(record) and len(record)-g['end0']<self.cfg['flank_nt'],
                                    path=path.relative_to(ROOT).as_posix()))
        SeqIO.write(fasta,self.out/'context_windows.fna','fasta')
        window_paths={r['anchor_protein']:r['path'] for r in window_rows}
        output=[]
        for key,entry in sorted(loci.items()):
            rtid=entry['rt_id'];rt=self.genes.get(rtid)
            pids=set().union(*entry['partners'].values()) if entry['partners'] else set()
            involved=pids|({rtid} if rtid else set())
            genes=[self.genes[p] for p in involved]
            art=self.accepted['ART_RT'].get(rtid)
            hint=repeat_hint(self.records[entry['accession']],rt,self.by_contig[entry['accession']]) if rt else dict(
                upstream_nt_available='',repeat_hint_copies='',repeat_hint_word='',repeat_hint_positions_1based='',repeat_hint_status='no_RT_anchor')
            coor=bool(rt and any(self.genes[p]['strand']==rt['strand'] for p in pids))
            row=dict(locus_id=key,accession=entry['accession'],source_path=genes[0]['source_path'],
                     locus_start_1based=min(g['start0'] for g in genes)+1,locus_end_1based=max(g['end0'] for g in genes),
                     RT_protein_id=rtid,RT_strand=rt['strand'] if rt else '',RT_length_aa=rt['length_aa'] if rt else '',
                     RT_start_1based=rt['start0']+1 if rt else '',RT_end_1based=rt['end0'] if rt else '',
                     RT_profile_bitscore=art['sequence_bitscore'] if art else '',
                     RT_profile_coverage=art['profile_coverage'] if art else '',
                     found_by_ART_protein_first=entry['art_first'],found_by_partner_first=entry['partner_first'],
                     discovery_routes=';'.join((["ART_protein_first"] if entry['art_first'] else [])+(["partner_first"] if entry['partner_first'] else [])),
                     Type_I_partner_ids=';'.join(sorted(entry['partners']['Type_I_partner'])),
                     Type_II_partner_ids=';'.join(sorted(entry['partners']['Type_II_partner'])),
                     Type_III_first_partner_ids=';'.join(sorted(entry['partners']['Type_III_first_partner'])),
                     Type_III_second_partner_ids=';'.join(sorted(entry['partners']['Type_III_second_partner'])),
                     partner_family_profiles=';'.join(f for f in FAMILIES[1:] if entry['partners'][f]),
                     cooriented_RT_partner_context=coor,
                     candidate_status='RT_partner_homology_and_compatible_context_provisional' if coor else
                     ('ART_RT_profile_hit_without_matched_partner' if art else 'partner_profile_hit_context_unresolved'),
                     contains_training_protein=bool(involved&training),training_protein_ids=';'.join(sorted(involved&training)),
                     context_anchor_ids=';'.join(sorted(entry['anchors'])),
                     context_window_paths=';'.join(window_paths[p] for p in sorted(entry['anchors'])),
                     unavailable_partner_search_families=missing,notes=';'.join(sorted(entry['context_notes'])),**hint)
            output.append(row)
        self.loci=output
        tsv(self.out/'candidate_loci.tsv',output)
        tsv(self.out/'partner_upstream_genes.tsv',upstream)
        tsv(self.out/'RT_downstream_genes.tsv',downstream)
        tsv(self.out/'windows.tsv',window_rows)
        candidates={row['RT_protein_id'] for row in output if row['RT_protein_id']}
        SeqIO.write([self.proteins[p] for p in sorted(candidates)],self.out/'candidate_RT_proteins.faa','fasta')
        for family in FAMILIES[1:]:
            ids=set().union(*(entry['partners'][family] for entry in loci.values())) if loci else set()
            SeqIO.write([self.proteins[p] for p in sorted(ids)],self.out/f'{family}.candidates.faa','fasta')
        return output

    def report(self):
        both=sum(r['found_by_ART_protein_first'] and r['found_by_partner_first'] for r in self.loci)
        rt_only=sum(r['found_by_ART_protein_first'] and not r['found_by_partner_first'] for r in self.loci)
        partner_only=sum(r['found_by_partner_first'] and not r['found_by_ART_protein_first'] for r in self.loci)
        summary=dict(genomes=len(self.records),annotated_proteins=len(self.proteins),candidate_loci=len(self.loci),
                     found_by_both=both,ART_first_only=rt_only,partner_first_only=partner_only,
                     profile_registry=self.registry,calibrated_familywide_threshold=False)
        (self.out/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
        rows=['# ART discovery through independent RT and partner searches','',
              f'Screened {len(self.proteins)} annotated proteins from {len(self.records)} genomes. '
              f'Candidate loci: {len(self.loci)}; both routes: {both}; RT-first only: {rt_only}; partner-first only: {partner_only}.','',
              'The search panel includes the six local GenBank genomes, matching public annotation for the local SA1 FASTA, '
              'and public LPJP1. Several target proteins are also training seeds, explicitly flagged in the locus table. '
              'These are recovery results, not an estimate of independent discovery sensitivity.','',
              '| Family | Unique seeds | Held-out recovery | Reviewed GNAT hits | Status |',
              '| --- | ---: | ---: | ---: | --- |']
        for r in self.registry:
            rows.append(f'| {r["family"]} | {r["unique_seed_sequences"]} | {r["heldout_positive_recovered"]}/{r["heldout_positive_count"]} | {r["reviewed_GNAT_hits"]} | {r["calibration_status"]} |')
        rows += ['', '## Discovery rules','',
                 'Each available family profile is searched against the entire protein database independently. '
                 'Partner searches are never restricted to proteins near prior RT hits. Sequence E <= 1e-5, '
                 'independent domain E <= 1e-5 and domain profile coverage >=60% define a provisional retrieval screen. '
              'Raw hits to E <= 1e-3 are retained. Calibration and discovery searches use the same fixed '
              'HMMER -Z and --domZ search space (the target proteome size). These settings are not calibrated '
              'family-membership probabilities or validated bit cutoffs.','',
                 'For every partner hit the complete gene plus 6 kb on each side is extracted, clipped at contig ends. '
                 'The first three upstream annotated CDSs are ranked on the partner transcription axis, including other-strand '
                 'and overlapping CDSs. All are reported with orientation, coordinates, ART-profile evidence and generic RT evidence. '
                 'A compatible association requires a co-oriented ART or generic RT among those three, within the extracted window. '
                 'Generic RT evidence alone never sets the RT-first discovery flag. Multiple compatible RTs remain explicit alternative associations.','',
                 'The type I seed is an accession-level reconstruction of the paper\'s LPJP1 representative: '
                 'array-associated RT QXN67953.1 and downstream overlapping gene QXN67952.1. The seed assignment is provisional, '
                 'not experimentally established. Type I hits require homology to this ART-associated protein, and genomic-context '
                 'compatibility is recorded separately. A GNAT product label never triggers a type I call.','',
                 '## Calibration limits','',
                 'MarsHill and SA1 supply held-out type II partner and RT sequences with no exact training-sequence duplicates; '
                 'they are related Staphylococcus references, not independent lineages. Three reviewed E. coli GNATs (RimI, RimJ, RimL) '
                 'are unrelated-family controls; up to 50 annotated phage housekeeping proteins are an additional background screen, '
                 'not experimentally verified noninteractors. The small control set cannot establish a family-wide false-positive rate.','',
                 'Type I has one provisional seed and no independent verified positive holdout; its single-sequence HMM is an exploratory '
                 'homology retrieval model, not a well-validated family profile. Type III first and second genes have separate slots, '
                 'empty seed FASTAs and explicit unavailable status. No placeholder HMMs or substitute generic domains are used. '
                 'Supply separate curated seed FASTAs and curation evidence in the configuration to activate them. Their absence means '
                 'the corresponding discovery routes were not assessed, not that no type III loci exist.','',
                 '## DNA and annotation limits','',
                 'The upstream exact-word repeat scan masks annotated CDSs >=300 nt and searches for >=3 approximately regularly '
                 'spaced copies of a non-low-complexity 10-mer. It is only an inspection hint: there is no shuffled-null significance '
                 'test or validated array absence call. The windows and full contigs are released for inspection. Only annotated '
                 'CDSs are searched; de novo ORF rescue, revised gene calling and systematic partial-CDS recovery remain future work.','',
                 '## Reproduce','',
                 '```powershell', '.venv/Scripts/python.exe scripts/art_discovery.py --config config/art_discovery.json', '```','',
                 'Main output: candidate_loci.tsv. See partner_upstream_genes.tsv, windows.tsv, profile_registry.tsv, '
                 'calibration/, searches/, seed_manifest.tsv, source_records.tsv and commands.jsonl for evidence.','',
                 'Method sources: supplied Yoon et al. PDF, family description and partner methods (pp. 8, 36-38); '
                 'https://doi.org/10.64898/2026.10.05.755809 ; LPJP1 MZ422438.1 at ENA; SA1 MW218148.1 at NCBI; '
                 'reviewed GNAT controls at UniProt; myRT class profiles https://github.com/mgtools/myRT .','']
        (self.out/'REPORT.md').write_text('\n'.join(rows),encoding='utf-8')
        print(json.dumps({k:v for k,v in summary.items() if k!='profile_registry'},indent=2),flush=True)

    def verify_and_hash(self):
        for source in self.sources:
            assert digest((ROOT/source['path']).read_bytes())==source['file_sha256']
        for row in self.loci:
            assert row['found_by_ART_protein_first'] or row['found_by_partner_first']
        expected={('MW248466.1','QQM14740.1'),('MW218148.1','QPI16926.1')}
        recovered={(r['accession'],r['RT_protein_id']) for r in self.loci if r['found_by_ART_protein_first'] and r['found_by_partner_first']}
        assert expected<=recovered,'MarsHill/SA1 dual-route recovery failed'
        assert len({r['locus_id'] for r in self.loci})==len(self.loci)
        with (self.out/'windows.tsv').open(encoding='utf-8') as f:
            window_anchors={r['anchor_protein'] for r in csv.DictReader(f,delimiter='\t')}
        for g in self.genes.values():
            if g['protein_id'] not in window_anchors:continue
            a,b=window_bounds(g,len(self.records[g['accession']]),self.cfg['flank_nt'])
            assert 0<=a<=g['start0']<g['end0']<=b<=len(self.records[g['accession']])
        (self.out/'verification.json').write_text(json.dumps(dict(source_hashes_unchanged=True,
            exact_holdout_sequence_overlap=False,unique_candidate_locus_ids=True,
            MarsHill_and_SA1_both_routes_recovered=True,source_records=len(self.records)),indent=2),encoding='utf-8')
        myrt_commit=subprocess.check_output(['git','-C',str(ROOT/'.tools/myRT'),'rev-parse','HEAD'],text=True).strip()
        hmmer_help=subprocess.check_output([str(self.bindir/'hmmbuild.exe'),'-h'],env=self.env,text=True)
        mafft_source=(ROOT/'.tools/mafft-win/usr/bin/mafft').read_text(encoding='utf-8')
        tool_versions=dict(python=sys.version,HMMER=re.search(r'HMMER (\S+)',hmmer_help)[1],
                           MAFFT=re.search(r'version="v([\d.]+)',mafft_source)[1],
                           generic_RT_reference_commit=myrt_commit,
                           generic_RT_profiles_sha256=digest((ROOT/'.tools/myRT/Models/HMM/RVT-All.hmm').read_bytes()),
                           hmmbuild_binary_sha256=digest((self.bindir/'hmmbuild.exe').read_bytes()),
                           hmmsearch_binary_sha256=digest((self.bindir/'hmmsearch.exe').read_bytes()))
        (self.out/'software_provenance.json').write_text(json.dumps(tool_versions,indent=2),encoding='utf-8')
        hashes={p.relative_to(self.out).as_posix():digest(p.read_bytes()) for p in self.out.rglob('*')
                if p.is_file() and p.name!='output_sha256.json'}
        (self.out/'output_sha256.json').write_text(json.dumps(hashes,indent=2),encoding='utf-8')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--config',type=Path,default=ROOT/'config/art_discovery.json')
    args=parser.parse_args()
    pipeline=Pipeline(args.config)
    pipeline.prepare();pipeline.build_and_search();pipeline.assemble();pipeline.report();pipeline.verify_and_hash()


if __name__=='__main__':main()
