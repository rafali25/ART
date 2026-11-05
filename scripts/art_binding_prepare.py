"""Freeze a bounded development matrix and independently generated chain MSAs."""
from datetime import datetime,timezone
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

from Bio import AlignIO,SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/ART binding/pilot_20261105'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    if not os.environ.get('SLURM_JOB_ID'):raise SystemExit('Use Slurm.')
    config_path=ROOT/'config/art_binding_pilot_20261105.json'
    cfg=json.loads(config_path.read_text())
    source=ROOT/'data/art_samples_2026-11-05/art_samples.csv'
    with source.open(newline='') as f:catalog=list(csv.DictReader(f))
    selected=[next(r for r in catalog if r['sample_id']==sid) for sid in cfg['sample_ids']]
    assert all(r['pair_complete']=='True' and r['partner_type']=='Type_II' and r['category1_ART_profile_pass']=='True' for r in selected)
    if (OUT/'manifest.json').exists():raise SystemExit('Preserve the existing frozen pilot; use a new version to change inputs.')
    for d in ['alignments','msas','inputs','code','logs']: (OUT/d).mkdir(parents=True,exist_ok=True)
    cohort=[r for r in catalog if r['partner_type']=='Type_II' and r['pair_complete']=='True']
    msa_paths={}
    for kind,column in [('RT','RT_protein_sequence'),('partner','partner_protein_sequence')]:
        raw=OUT/'alignments'/(kind+'.faa');aligned=OUT/'alignments'/(kind+'.aligned.faa')
        SeqIO.write([SeqRecord(Seq(r[column]),id=r['sample_id'],description='') for r in cohort],raw,'fasta')
        with aligned.open('w') as stdout,(OUT/'logs'/(kind+'.mafft.log')).open('w') as stderr:
            subprocess.run([str(ROOT/'.venv/bin/mafft'),'--localpair','--maxiterate','1000','--thread','2','--inputorder',str(raw)],stdout=stdout,stderr=stderr,check=True)
        alignment={r.id:str(r.seq) for r in AlignIO.read(aligned,'fasta')}
        for query in selected:
            sid=query['sample_id'];q=alignment[sid];records=[]
            for other in [sid]+sorted(k for k in alignment if k!=sid):
                a3m=''.join(h.lower() if c=='-' else h for c,h in zip(q,alignment[other]) if c!='-' or h!='-')
                assert len(''.join(c for c in a3m if not c.islower()))==len(query[column])
                records.append(dict(sequence=a3m,key=''))
            assert records[0]['sequence']==query[column]
            path=OUT/'msas'/(sid+'_'+kind+'.csv')
            with path.open('w',newline='') as f:
                writer=csv.DictWriter(f,['sequence','key']);writer.writeheader();writer.writerows(records)
            msa_paths[sid,kind]=path
    cells=[]
    # Native/alternative labels are stored only in the manifest, never in model inputs.
    for condition in cfg['conditions']:
        for i,rt in enumerate(selected):
            for j,partner in enumerate(selected):
                name=f'{condition}_r{i+1}_p{j+1}';path=OUT/'inputs'/(name+'.json')
                content=dict(version=1,sequences=[dict(protein=dict(id='A',sequence=rt['RT_protein_sequence'],msa='empty' if condition=='S' else str(msa_paths[rt['sample_id'],'RT']))),
                    dict(protein=dict(id='B',sequence=partner['partner_protein_sequence'],msa='empty' if condition=='S' else str(msa_paths[partner['sample_id'],'partner'])))])
                path.write_text(json.dumps(content,indent=2)+'\n')
                # The Boltz parser accepts JSON through its YAML loader; use the
                # documented .yaml extension without adding a prep dependency.
                yaml_path=path.with_suffix('.yaml');path.replace(yaml_path)
                for seed in cfg['initial_seeds']:
                    cells.append(dict(cell_id=name,condition=condition,RT_sample_id=rt['sample_id'],partner_sample_id=partner['sample_id'],
                        native_pair=i==j,seed=seed,input_path=str(yaml_path.relative_to(ROOT)),input_sha256=sha(yaml_path),
                        RT_length_aa=len(rt['RT_protein_sequence']),partner_length_aa=len(partner['partner_protein_sequence'])))
    with (OUT/'cells.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,list(cells[0]),delimiter='\t');w.writeheader();w.writerows(cells)
    cfg_copy=OUT/'registered_config.json';shutil.copy2(config_path,cfg_copy)
    for name in ['art_binding_prepare.py','art_binding_msa_control.py','art_binding_boltz_entry.py','art_binding_run.py','art_binding_analyse.py','test_art_binding_controls.py']:
        p=ROOT/'scripts'/name
        if p.exists():shutil.copy2(p,OUT/'code'/name)
    manifest=dict(experiment_id=cfg['experiment_id'],status='inputs_frozen',prepared_utc=datetime.now(timezone.utc).isoformat(),
        slurm_job_id=os.environ['SLURM_JOB_ID'],config=cfg,cells=cells,development_loci=[r['locus_id'] for r in selected],
        homolog_samples=[r['sample_id'] for r in cohort],input_source_sha256=sha(source),
        hashes={str(p.relative_to(OUT)):sha(p) for p in OUT.rglob('*') if p.is_file()},
        API_key_required=False,calibration_complete=False,confirmatory_inference_allowed=False)
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (OUT/'README.md').write_text('''# ART binding development pilot — 5 November 2026

SA1 and MarsHill form a two-reference 2×2 matrix under query-only (S) and
explicitly unpaired (U) conditions. The initial stage uses one fixed seed per
cell: eight predictions. This measures input controls, execution, scoring and
runtime; it is not a confirmatory test of binding or specificity. The two loci
are close relatives; alternative pairings are not experimentally verified
nonbinders. Their close clades are reserved from future confirmation.

U uses shallow local homolog alignments. Native Boltz unpaired flags do not
guarantee that homologs occupy separate rows; a recorded input adapter separates
every non-query homolog into a chain-specific row. Processed matrices and audits
are retained, and the exact adapter must be used during subsequent calibration.
S is query-only. Templates, interaction restraints and affinity prediction are
disabled; no external MSA service or API key is used.

The current catalogue has zero eligible three-locus confirmatory blocks. Family,
independence, calibration and power requirements remain unresolved. Inspect
`manifest.json`, `cells.tsv`, `runs/`, and the eventual `results.tsv` for actual
execution status. Never interpret a one-seed pilot contrast as a binding result.
''')
    print(json.dumps(dict(prepared=True,cells=len(cells),homolog_samples=len(cohort),output=str(OUT)),indent=2),flush=True)


if __name__=='__main__':main()
