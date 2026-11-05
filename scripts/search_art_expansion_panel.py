"""Bounded, separate ART retrieval pilot on public related/huge phage genomes.

This does not alter frozen cohorts or assign final array/family status. Native
executables are resolved from the current Python environment before PATH.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import hashlib
from html.parser import HTMLParser
import io
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import urllib.request

from Bio import Align, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pyrodigal_gv

from art_reconstruction_common import parse_domains
from audit_art_type1_fasta import matching_annotations

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/art_sample_expansion_20261104'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def table(path, rows, fields):
    with path.open('w', newline='', encoding='utf8') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


class MetadataTable(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows = []
        self.row = None
        self.cell = None

    def handle_starttag(self, tag, attrs):
        if tag == 'tr':
            self.row = []
        if tag in ('td', 'th') and self.row is not None:
            self.cell = []

    def handle_data(self, text):
        if self.cell is not None:
            self.cell.append(text)

    def handle_endtag(self, tag):
        if tag in ('td', 'th') and self.cell is not None:
            self.row.append(' '.join(''.join(self.cell).split()))
            self.cell = None
        if tag == 'tr' and self.row is not None:
            if self.row:
                self.rows.append(self.row)
            self.row = None


def discover():
    meta = OUT / 'panel_metadata'
    selected = {}
    for r in csv.DictReader(io.StringIO((meta/'huge_phage_project.response').read_text()), delimiter='\t'):
        selected[r['accession']] = dict(accession=r['accession'], description=r['description'],
            expected_nt=int(r['base_count']), stratum='Al_Shayeb_2020_PRJEB35371',
            selection='all_sequence_records_returned_by_project_query')
    for label, cap in [('Listeria', 40), ('Wallmarkvirinae', 20), ('BpSP', 4), ('AR9', 4)]:
        path = meta / ('inphared_' + label + '_selected.html')
        if not path.exists():
            continue
        parser = MetadataTable()
        parser.feed(path.read_text())
        n = 0
        for cells in parser.rows[1:]:
            if len(cells) < 3 or not re.fullmatch(r'[A-Z]+_?\d+(?:\.\d+)?', cells[0]):
                continue
            accession, description, length = cells[:3]
            if accession in selected:
                continue
            selected[accession] = dict(accession=accession, description=description,
                expected_nt=int(length.replace(',', '')), stratum='INPHARED2_'+label,
                selection='length_descending_bounded_metadata_panel')
            n += 1
            if n == cap:
                break
    known = set()
    for folder in ['art_discovery','art_expansion','art_type1_public','art_type1_mgnify']:
        path = ROOT/'results'/folder/'source_records.tsv'
        if path.exists():
            known.update(r['accession'].split('.')[0] for r in csv.DictReader(path.open(),delimiter='\t') if r.get('accession'))
    # Two current ICTV-listed near relatives absent from the previous export.
    for accession, name in [('OR455461','Staphylococcus phage AH12'),('OR770614','vB_StaM_PB50')]:
        selected.setdefault(accession,dict(accession=accession,description=name,expected_nt=0,
            stratum='ICTV_2026_Machiasvirus',selection='named_near_relative_not_in_previous_catalog'))
    rows = [dict(r,previously_searched_source=r['accession'].split('.')[0] in known) for r in selected.values()]
    table(OUT/'panel_manifest.tsv', rows,
        ['accession','description','expected_nt','stratum','selection','previously_searched_source'])
    print(json.dumps(dict(selected=len(rows), expected_nt=sum(r['expected_nt'] for r in rows))), flush=True)


def acquire():
    rows = list(csv.DictReader((OUT/'panel_manifest.tsv').open(), delimiter='\t'))
    if sum(int(r['expected_nt']) for r in rows) > 50_000_000:
        raise ValueError('Panel exceeds the 50-Mb DNA pilot budget')
    source = OUT/'panel_sources'
    source.mkdir(exist_ok=True)

    def fetch(r):
        fmt = 'genbank' if r['accession'].startswith('NC_') else 'embl'
        path = source/(r['accession']+('.gb' if fmt=='genbank' else '.embl'))
        url = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gbwithparts&retmode=text&id='+r['accession']
            if fmt=='genbank' else 'https://www.ebi.ac.uk/ena/browser/api/embl/'+r['accession']+'?download=false')
        try:
            if not path.exists():
                req = urllib.request.Request(url, headers={'User-Agent':'ART-public-discovery/1.0'})
                with urllib.request.urlopen(req, timeout=35) as response:
                    raw = response.read(5_000_001)
                    headers = dict(response.headers)
                if len(raw) > 5_000_000:
                    raise ValueError('Per-genome 5-MB download cap exceeded')
                rec = SeqIO.read(io.StringIO(raw.decode()), fmt)
                if int(r['expected_nt']) and len(rec) != int(r['expected_nt']):
                    raise ValueError('Source length differs from frozen metadata')
                path.write_bytes(raw)
                path.with_suffix('.request.json').write_text(json.dumps(dict(
                    url=url, retrieved_utc=datetime.now(timezone.utc).isoformat(),
                    bytes=len(raw), sha256=digest(raw), response_headers=headers), indent=2))
            rec = SeqIO.read(path, fmt)
            result = dict(r, source_accession_version=rec.id, status='retrieved',
                path=path.relative_to(ROOT).as_posix(), format=fmt, file_sha256=digest(path.read_bytes()),
                sequence_sha256=digest(str(rec.seq).upper().encode()),
                original_CDS=sum(f.type=='CDS' for f in rec.features), length_nt=len(rec))
        except Exception as e:
            result = dict(r, status='retrieval_failed', reason=type(e).__name__+': '+str(e))
        print(json.dumps({k:result[k] for k in ('accession','status')}), flush=True)
        return result

    results = list(ThreadPoolExecutor(max_workers=3).map(fetch, rows))
    fields = list(dict.fromkeys(k for r in results for k in r))
    table(OUT/'panel_acquisition.tsv', results, fields)


def screen():
    out = OUT/'panel_screen'
    out.mkdir(exist_ok=True)
    rows = list(csv.DictReader((OUT/'panel_acquisition.tsv').open(), delimiter='\t'))
    proteins, genes, sources, seen = [], [], [], {}
    finder = pyrodigal_gv.ViralGeneFinder(meta=True)
    for r in rows:
        if r['status'] != 'retrieved':
            continue
        rec = SeqIO.read(ROOT/r['path'], r['format'])
        if digest((ROOT/r['path']).read_bytes()) != r['file_sha256']:
            raise ValueError('Source checksum changed')
        dh = digest(str(rec.seq).upper().encode())
        if dh in seen:
            sources.append(dict(accession=rec.id, stratum=r['stratum'], duplicate_of=seen[dh]))
            continue
        seen[dh] = rec.id
        intervals = set()
        for i, f in enumerate(rec.features):
            if f.type != 'CDS' or not f.qualifiers.get('translation'):
                continue
            peptide = f.qualifiers['translation'][0].replace(' ', '').rstrip('*')
            code = int(f.qualifiers.get('transl_table', ['11'])[0])
            codon_start = int(f.qualifiers.get('codon_start', ['1'])[0])-1
            coding = f.extract(rec.seq)[codon_start:]
            computed = str(coding.translate(table=code)).rstrip('*')
            qc = computed == peptide or (computed[1:] == peptide[1:] and peptide.startswith('M'))
            key = rec.id+'__source'+str(i)
            if set(peptide) <= set('ACDEFGHIKLMNPQRSTVWY'):
                proteins.append(SeqRecord(Seq(peptide), id=key, description='original_source_CDS'))
                genes.append(dict(key=key, accession=rec.id, start0=int(f.location.start),
                    end0=int(f.location.end), strand=f.location.strand, method='source_CDS',
                    genetic_code=code, translation_matches=qc,
                    complete_location=not bool(re.search(r'[<>]', str(f.location))),
                    protein_id=f.qualifiers.get('protein_id', [''])[0],ranked_neighborhood=True))
                intervals.add((int(f.location.start),int(f.location.end),f.location.strand,peptide))
        predicted = finder.find_genes(bytes(rec.seq))
        for i, g in enumerate(predicted, 1):
            peptide = g.translate().rstrip('*')
            if (g.begin-1,g.end,g.strand,peptide) in intervals:
                continue
            key = rec.id+'__pgv'+str(i)
            proteins.append(SeqRecord(Seq(peptide), id=key, description='alternative_viral_gene_call'))
            genes.append(dict(key=key, accession=rec.id, start0=g.begin-1,end0=g.end,
                strand=g.strand,method='Pyrodigal_gv',genetic_code=g.translation_table,
                translation_matches=True,complete_location=not(g.partial_begin or g.partial_end),
                protein_id='',ranked_neighborhood=not bool(matching_annotations(
                    dict(start0=g.begin-1,end0=g.end,strand=g.strand),rec.features))))
        sources.append(dict(accession=rec.id,stratum=r['stratum'],length_nt=len(rec),
            gene_predictions=len(predicted),genetic_codes=sorted({g.translation_table for g in predicted}),
            model=str(predicted.metagenomic_bin),duplicate_of=''))
    fasta = out/'proteins.faa'
    SeqIO.write(proteins, fasta, 'fasta')
    table(out/'genes.tsv', genes, list(genes[0]))
    table(out/'sources.tsv', sources, list(dict.fromkeys(k for r in sources for k in r)))
    profiles = {
        'ART_RT':ROOT/'results/art_discovery/profiles/ART_RT.hmm',
        'Type_I_partner':ROOT/'results/art_discovery/profiles/Type_I_partner.hmm',
        'Type_II_partner':ROOT/'results/art_discovery/profiles/Type_II_partner.hmm',
        'Type_I_two_seed_exploratory':ROOT/'results/art_type1_public/profiles/TypeI_two_seed_exploratory.hmm',
        'RVT_1':ROOT/'.tools/art-reconstruction/PF00078.hmm',
    }
    exe = Path(sys.executable).parent/'hmmsearch'
    if not exe.exists():
        exe = Path(shutil.which('hmmsearch') or '')
    hits, commands = [], []
    bykey = {r['key']:r for r in genes}
    for family, profile in profiles.items():
        args = [str(exe),'--cpu','4','--seed','42','-E','0.001','--domE','0.001',
            '--domtblout',str(out/(family+'.domtbl')),str(profile),str(fasta)]
        with (out/(family+'.log')).open('w') as log:
            subprocess.run(args,stdout=log,stderr=subprocess.STDOUT,check=True)
        commands.append(dict(args=args,profile_sha256=digest(profile.read_bytes()),
            search_space_protein_records=len(proteins)))
        for h in parse_domains(out/(family+'.domtbl')):
            accepted = h['sequence_evalue']<=1e-5 and h['domain_evalue']<=1e-5 and h['profile_coverage']>=.6
            hits.append(dict(family=family,passes_retrieval_threshold=accepted,
                **h,**bykey[h['target']]))
    fields = list(dict.fromkeys(k for r in hits for k in r)) or ['family','target']
    table(out/'profile_hits.tsv',hits,fields)
    # Keep paired gene representations and their evidence explicit. This is a
    # retrieval table: array testing and sequence-placement review follow later.
    accepted = [h for h in hits if h['passes_retrieval_threshold']]
    rt = {h['target']:h for h in accepted if h['family']=='ART_RT'}
    partners = [h for h in accepted if h['family'].startswith('Type_')]
    pairs = []
    for h in partners:
        axis = lambda g: int(g['start0']) if int(h['strand'])==1 else -int(g['end0'])
        neighborhood = sorted([g for g in genes if g['accession']==h['accession'] and
            g['ranked_neighborhood'] and g['key']!=h['target'] and axis(g)<axis(h)],key=axis,reverse=True)[:3]
        for i,g in enumerate(neighborhood,1):
            if g['key'] in rt and g['strand']==h['strand']:
                if any(x['key'] in rt for x in neighborhood[:i-1]):
                    continue
                pairs.append(dict(accession=h['accession'],RT_key=g['key'],partner_key=h['target'],
                    partner_family=h['family'],upstream_gene_rank=i,
                    RT_complete=g['complete_location'],partner_complete=h['complete_location'],
                    RT_source_translation_matches=g['translation_matches'],partner_source_translation_matches=h['translation_matches'],
                    status='provisional_profile_and_context_pair_requires_sequence_and_array_review'))
    table(out/'candidate_pairs.tsv',pairs,list(pairs[0]) if pairs else ['accession','RT_key','partner_key','status'])
    review = [h for h in hits if h['family']=='ART_RT' and not h['passes_retrieval_threshold']]
    table(out/'RT_review_hits.tsv',review,fields)
    (out/'commands.json').write_text(json.dumps(commands,indent=2))
    summary = dict(unique_DNA_sequences=len(seen),genomes_with_source_records=len(sources),
        protein_records=len(proteins),reported_domain_hits=len(hits),
        accepted_by_profile={family:len({h['target'] for h in hits if h['family']==family and h['passes_retrieval_threshold']}) for family in profiles},
        candidate_pair_representations=len(pairs),candidate_pair_source_genomes=len({p['accession'] for p in pairs}),
        interpretation='provisional_retrieval_only_not_final_ART_or_array_classification',
        Pyrodigal_gv_version=pyrodigal_gv.__version__,source_annotated_and_predicted_representations_retained=True,
        caveat='This pilot screens source CDSs and alternative viral calls; six-frame rescue is a follow-up around hits, not a genome-wide completed test.')
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2),flush=True)


def export():
    """Export new source-CDS pairs separately from alternative representations."""
    from dataclasses import asdict
    from art_array_detector import initial_scan, refinement_scan
    from extract_art_loci import plan_locus, extract_record, write_gff
    out = OUT/'panel_screen'
    loci = out/'new_loci'
    loci.mkdir(exist_ok=True)
    genes = {r['key']:r for r in csv.DictReader((out/'genes.tsv').open(),delimiter='\t')}
    proteins = {r.id:r for r in SeqIO.parse(out/'proteins.faa','fasta')}
    pairs = list(csv.DictReader((out/'candidate_pairs.tsv').open(),delimiter='\t'))
    sources = {r.get('source_accession_version'):r for r in csv.DictReader((OUT/'panel_acquisition.tsv').open(),delimiter='\t') if r['status']=='retrieved'}
    refs = {}
    paths = list((ROOT/'results').glob('*/candidate_RT_proteins.faa'))
    paths += list((ROOT/'results/art_expansion').glob('stage*/candidate_RT_proteins.faa'))
    for path in paths:
        for r in SeqIO.parse(path,'fasta'):
            refs.setdefault(r.id,r)
    aligner = Align.PairwiseAligner(scoring='blastp')
    aligner.mode = 'global'
    cfg = json.loads((ROOT/'config/art_array_search.json').read_text())['art_array_search']
    rows, seen = [], set()
    for pair in pairs:
        accession = pair['accession']
        if accession not in {'OR455461.1','OR770614.1'}:
            continue
        rt,partner = genes[pair['RT_key']],genes[pair['partner_key']]
        if partner['method']!='source_CDS' or (accession,rt['key']) in seen:
            continue
        seen.add((accession,rt['key']))
        src = sources[accession]
        record = SeqIO.read(ROOT/src['path'],src['format'])
        def feature(g):
            return next(f for f in record.features if f.type=='CDS' and f.qualifiers.get('protein_id',[''])[0]==g['protein_id'])
        rf,pf = feature(rt),feature(partner)
        lid = accession+'__'+rt['protein_id']
        plan = plan_locus(record,rf,[pf])
        locus = extract_record(record,rf,plan,lid)
        SeqIO.write(locus,loci/(lid+'.gb'),'genbank')
        SeqIO.write(locus,loci/(lid+'.fna'),'fasta')
        write_gff(loci/(lid+'.gff3'),locus)
        serial = dict(plan,segments=[asdict(s) for s in plan['segments']],plus_segments=[asdict(s) for s in plan['plus_segments']])
        (loci/(lid+'.map.json')).write_text(json.dumps(serial,indent=2))
        seqs = [SeqRecord(proteins[g['key']].seq,id=g['protein_id'],description=lid) for g in [rt,partner]]
        SeqIO.write(seqs,loci/(lid+'.proteins.faa'),'fasta')
        SeqIO.write([SeqRecord(f.extract(record.seq),id=g['protein_id'],description=lid) for f,g in [(rf,rt),(pf,partner)]],loci/(lid+'.CDS.fna'),'fasta')
        if int(rt['strand'])==1:
            upstream=record.seq[max(0,int(rt['start0'])-6000):int(rt['start0'])]
        else:
            upstream=record.seq[int(rt['end0']):int(rt['end0'])+6000].reverse_complement()
        initial = initial_scan(str(upstream[-3000:]),cfg)
        extended = refinement_scan(str(upstream),cfg)
        (loci/(lid+'.exploratory_array.json')).write_text(json.dumps(dict(initial=initial,extended=extended,
            status='uncalibrated_diagnostic_only_no_p_value_or_array_acceptance'),indent=2))
        identity,closest=max((aligner.align(seqs[0].seq,r.seq)[0].counts().identities/max(len(seqs[0]),len(r)),name) for name,r in refs.items())
        rows.append(dict(locus_id=lid,source_accession_version=accession,RT_protein_id=rt['protein_id'],
            partner_protein_id=partner['protein_id'],RT_start0=rt['start0'],RT_end0=rt['end0'],
            partner_start0=partner['start0'],partner_end0=partner['end0'],strand=rt['strand'],
            partner_family='Type_II',protein_CDS_translation_matches=True,
            RT_length_aa=len(seqs[0]),partner_length_aa=len(seqs[1]),closest_previous_RT_key=closest,
            closest_RT_identity_over_longer_length=identity,
            exact_RT_sequence_previously_exported=any(seqs[0].seq==r.seq for r in refs.values()),
            upstream_available_nt=len(upstream),exploratory_extended_seed_copies=len(extended['candidates'][0]['positions']) if extended['candidates'] else 0,
            array_status='not_calibrated',status='complete_type_II_profile_context_candidate_new_to_workspace_not_independent_lineage_claim'))
    table(out/'new_complete_pair_candidates.tsv',rows,list(rows[0]) if rows else ['locus_id','status'])
    (out/'review_summary.json').write_text(json.dumps(dict(new_complete_type_II_candidates=len(rows),
        exact_new_RT_sequences=sum(not r['exact_RT_sequence_previously_exported'] for r in rows),
        array_tests_not_performed=True,source_CDS_pairs_deduplicated=True),indent=2))
    print(json.dumps(rows,indent=2),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation',choices=['discover','acquire','screen','export'])
    args = parser.parse_args()
    OUT.mkdir(exist_ok=True)
    dict(discover=discover,acquire=acquire,screen=screen,export=export)[args.operation]()
