"""Collect technical pilot outputs; no confirmatory specificity inference."""
import csv
import json
import os
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser,NeighborSearch

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/ART binding/pilot_20261105'


def main():
    if not os.environ.get('SLURM_JOB_ID'):raise SystemExit('Use Slurm.')
    manifest=json.loads((OUT/'manifest.json').read_text())
    with (ROOT/'data/art_samples_2026-11-05/art_samples.csv').open(newline='') as f:samples={r['sample_id']:r for r in csv.DictReader(f)}
    rows=[]
    for cell in manifest['cells']:
        folder=OUT/'runs'/(cell['cell_id']+'_seed'+str(cell['seed']))
        if not (folder/'status.json').exists():continue
        state=json.loads((folder/'status.json').read_text())
        if state['status']!='complete':continue
        confidence=json.loads((ROOT/state['confidence_file']).read_text())
        cif=next(folder.glob('boltz_results_*/predictions/*/*_model_0.cif'))
        model=MMCIFParser(QUIET=True,auth_chains=False,auth_residues=False).get_structure(cell['cell_id'],cif)[0]
        chain_a=model['A'];chain_b=model['B'];assert len(list(chain_a.get_residues()))==cell['RT_length_aa'];assert len(list(chain_b.get_residues()))==cell['partner_length_aa']
        atoms=[a for a in model.get_atoms() if a.element not in ['H','D']]
        search=NeighborSearch(atoms);pairs=search.search_all(5.0,level='A');contacts=set()
        for first,second in pairs:
            ra,rb=first.get_parent(),second.get_parent()
            ca,cb=ra.get_parent().id,rb.get_parent().id
            if ca==cb:continue
            if ca=='B':ra,rb=rb,ra
            contacts.add((ra.id[1]-1,rb.id[1]-1))
        pae_path=next(folder.glob('boltz_results_*/predictions/*/pae_*_model_0.npz'))
        arrays=np.load(pae_path);pae=arrays['pae'];n=cell['RT_length_aa']
        if pae.ndim==3:pae=pae[0]
        assert pae.shape==(n+cell['partner_length_aa'],)*2
        contact_pae=float(np.mean([(pae[i,n+j]+pae[n+j,i])/2 for i,j in contacts])) if contacts else None
        metrics=json.loads((folder/'audit/execution_metrics.json').read_text())
        audit=json.loads(next((folder/'audit').glob('feature_audit_*.json')).read_text())
        ntd=int(samples[cell['RT_sample_id']]['NTD_length_aa'])
        rows.append(dict(cell_id=cell['cell_id'],condition=cell['condition'],RT_sample_id=cell['RT_sample_id'],partner_sample_id=cell['partner_sample_id'],
            native_pair=cell['native_pair'],seed=cell['seed'],iptm=confidence['iptm'],ptm=confidence.get('ptm'),
            residue_contacts_5A=len(contacts),operational_NTD_contact_pairs=sum(i<ntd for i,j in contacts),
            mean_bidirectional_contact_PAE=contact_pae,seconds=metrics['seconds'],peak_GPU_allocated_bytes=metrics['peak_GPU_allocated_bytes'],
            MSA_rows=audit['rows'],nonquery_copaired_rows=audit['nonquery_multiple_chain_rows'],
            confidence_file=state['confidence_file'],structure_file=str(cif.relative_to(ROOT))))
    if rows:
        with (OUT/'results.tsv').open('w',newline='') as f:
            w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    contrasts={}
    for condition in manifest['config']['conditions']:
        group=[r for r in rows if r['condition']==condition]
        if len(group)==4:
            contrasts[condition]=dict(native_mean=float(np.mean([r['iptm'] for r in group if r['native_pair']])),
                alternative_mean=float(np.mean([r['iptm'] for r in group if not r['native_pair']])))
            contrasts[condition]['descriptive_difference']=contrasts[condition]['native_mean']-contrasts[condition]['alternative_mean']
    result=dict(completed_predictions=len(rows),planned_predictions=len(manifest['cells']),descriptive_pilot_contrasts=contrasts,
        confirmatory_inference=False,experimentally_validated_binding=False,calibration_complete=False,
        limitations=manifest['config']['limitations'],slurm_job_id=os.environ['SLURM_JOB_ID'])
    (OUT/'pilot_summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
