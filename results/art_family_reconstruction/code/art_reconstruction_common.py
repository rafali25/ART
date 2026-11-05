"""Shared, recorded command execution and sequence-only reconstruction utilities."""
from collections import defaultdict
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

from extract_art_loci import ROOT,read_tsv,write_tsv,write_json


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def table(path,rows,fields=None):
    write_tsv(path,rows,fields or list(dict.fromkeys(k for r in rows for k in r)) or ['status'])


class Runner:
    def __init__(self,out):
        self.out=Path(out);self.out.mkdir(parents=True,exist_ok=True)
        self.metadata=json.loads((ROOT/'.tools/art-reconstruction/provenance.json').read_text())
        self.binaries={k:(ROOT/v).resolve() for k,v in self.metadata['binaries'].items()}
        self.rvt_accession=next(line.split()[1] for line in self.binaries['RVT_1'].read_text().splitlines() if line.startswith('ACC '))
        # Direct invocation avoids the upstream .bat wrapper's UAC fallback and
        # its unconditional success exit. BusyBox helpers are local hardlinks.
        self.binaries['mmseqs']=self.binaries['mmseqs'].parent/'bin/mmseqs.exe'
        self.env=os.environ.copy()
        self.env['PATH']=str(self.binaries['mmseqs'].parent)+os.pathsep+str(ROOT/'.tools/cygwin/bin')+os.pathsep+self.env['PATH']
    def run(self,args,log,stdout=None,expected=()):
        args=[x.resolve().as_posix() if isinstance(x,Path) else str(x) for x in args]
        start=time.perf_counter();status=None
        try:
            with Path(log).open('w',encoding='utf8') as err:
                if stdout:
                    with Path(stdout).open('w',encoding='utf8') as dest:
                        result=subprocess.run(args,env=self.env,stdout=dest,stderr=err,check=True)
                else:result=subprocess.run(args,env=self.env,stdout=err,stderr=subprocess.STDOUT,check=True)
            status=result.returncode
            for p in expected:
                if not Path(p).exists() or Path(p).stat().st_size==0:raise ValueError('Missing tool output: '+str(p))
        except subprocess.CalledProcessError as exc:
            status=exc.returncode
            raise
        finally:
            with (self.out/'commands.jsonl').open('a',encoding='utf8') as dest:
                dest.write(json.dumps(dict(args=args,log=str(log),stdout=str(stdout) if stdout else None,
                    elapsed_seconds=round(time.perf_counter()-start,4),returncode=status))+'\n')
        return time.perf_counter()-start
    def hmm(self,hmm,fasta,stem,z=2076,report_e=.001):
        stem=Path(stem)
        elapsed=self.run([ROOT/'.tools/hmmer-install/bin/hmmsearch.exe','--cpu','4','--seed','42',
            '-Z',z,'--domZ',z,'-E',report_e,'--domE',report_e,'--tblout',str(stem)+'.tbl',
            '--domtblout',str(stem)+'.domtbl',hmm,fasta],str(stem)+'.log',expected=[str(stem)+'.domtbl'])
        return parse_domains(str(stem)+'.domtbl'),elapsed


def parse_domains(path):
    rows=[]
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith('#'):continue
        c=line.split(maxsplit=22)
        rows.append(dict(target=c[0],profile=c[3],sequence_evalue=float(c[6]),domain_evalue=float(c[12]),
            sequence_bits=float(c[7]),domain_bits=float(c[13]),profile_coverage=(int(c[16])-int(c[15])+1)/int(c[5]),
            hmm_start0=int(c[15])-1,hmm_end0=int(c[16]),protein_alignment_start0=int(c[17])-1,
            protein_alignment_end0=int(c[18]),protein_envelope_start0=int(c[19])-1,protein_envelope_end0=int(c[20])))
    return rows


def passes(row,cfg):
    return row['sequence_evalue']<=cfg['sequence_evalue'] and row['domain_evalue']<=cfg['domain_evalue'] and row['profile_coverage']>=cfg['minimum_ART_profile_coverage']


def connected_component(edges,anchor):
    """Sequence links only: architecture is deliberately absent from this API."""
    graph=defaultdict(set)
    for a,b in edges:graph[a].add(b);graph[b].add(a)
    seen={anchor};todo=[anchor]
    while todo:
        for nxt in graph[todo.pop()]-seen:seen.add(nxt);todo.append(nxt)
    return seen


def longest_contig_representatives(cluster_rows,metadata):
    groups=defaultdict(list)
    for center,member in cluster_rows:groups[center].append(member)
    return {center:sorted(members,key=lambda p:(-int(metadata[p]['contig_length_nt']),p))[0] for center,members in groups.items()}


def validation_gate():
    summary=json.loads((ROOT/'results/art_validation/summary.json').read_text())
    reasons=[]
    if summary['sequence_cases_passing']!=summary['sequence_cases']:reasons.append('repeat_copy_validation_failures')
    if summary['missing_comparison_sets']:reasons.append('missing_comparison_sets:'+','.join(summary['missing_comparison_sets']))
    if summary['missing_partner_families']:reasons.append('missing_partner_families:'+','.join(summary['missing_partner_families']))
    if not summary.get('full_boundary_validation_passed',False):reasons.append('full_repeat_boundaries_not_independently_confirmed')
    return dict(passed=not reasons,reasons=reasons,summary_sha256=sha(ROOT/'results/art_validation/summary.json'),
        scope='small_mapped_protein_pilot_allowed; definitive_locus_calls_and_further_expansion_gated')


def stage_gate(stage,cfg,metrics=None):
    if stage==1:return dict(allowed=True,mode='bounded_provisional_protein_pilot',validation=validation_gate())
    gate=validation_gate();reasons=list(gate['reasons'])
    if not metrics or not metrics.get('records_searched') or 'search_seconds' not in metrics or 'new_RT_candidates' not in metrics:
        reasons.append('previous_stage_cost_and_yield_not_measured')
    if stage==3 and not cfg['stage3_enabled']:reasons.append('large_collection_manifest_and_budget_not_enabled')
    if stage==3 and (not cfg.get('stage3_assembly_manifest') or cfg.get('stage3_maximum_records',0)<=0):
        reasons.append('large_collection_manifest_or_budget_missing')
    return dict(allowed=not reasons,mode='gated',reasons=reasons)
