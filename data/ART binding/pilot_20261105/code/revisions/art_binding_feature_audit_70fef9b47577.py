"""Exercise actual Boltz parsing and featurization on CPU before GPU prediction."""
import json
import os
from pathlib import Path
import hashlib
import shutil

import numpy as np
import torch

from art_binding_msa_control import separate_chains,audit_msa


def main():
    if not os.environ.get('SLURM_JOB_ID'):raise SystemExit('Use Slurm.')
    root=Path(__file__).resolve().parents[1];out=root/'data/ART binding/pilot_20261105'
    target=out/'cpu_feature_audit';target.mkdir(exist_ok=True)
    manifest=json.loads((out/'manifest.json').read_text());cache=root/'.tools/boltz_cache'
    from boltz.main import process_inputs
    from rdkit import Chem
    # Match the official CLI: atom names must survive multiprocessing pickle.
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
    from boltz.data.types import Manifest
    from boltz.data import const
    from boltz.data.mol import load_canonicals,load_molecules
    from boltz.data.module.inferencev2 import load_input
    from boltz.data.tokenize.boltz2 import Boltz2Tokenizer
    import boltz.data.feature.featurizerv2 as module
    process_inputs(data=[root/c['input_path'] for c in manifest['cells']],out_dir=target,
        ccd_path=cache/'ccd.pkl',mol_dir=cache/'mols',msa_server_url='',msa_pairing_strategy='greedy',
        max_msa_seqs=64,use_msa_server=False,boltz2=True,preprocessing_threads=2)
    parsed=Manifest.load(target/'processed/manifest.json')
    original=module.construct_paired_msa;current={};audits=[];canonicals=load_canonicals(cache/'mols')
    def controlled(*args,**kwargs):
        data=kwargs.get('data',args[0] if args else None)
        raw=original(*args,**kwargs);arrays=[x.numpy() for x in raw]
        chains=np.asarray(data.tokens['asym_id']);gap=const.token_ids['-']
        if current['condition']=='U':arrays=separate_chains(*arrays,chains,gap)
        audit=audit_msa(arrays[0],arrays[2],chains,gap,current['condition'])
        audit['cell_id']=current['cell_id'];audit['native_rows']=int(raw[0].shape[1]);audits.append(audit)
        return tuple(torch.as_tensor(a,dtype=t.dtype) for a,t in zip(arrays,raw))
    module.construct_paired_msa=controlled
    for record in parsed.records:
        current.update(next(c for c in manifest['cells'] if c['cell_id']==record.id))
        data=load_input(record=record,target_dir=target/'processed/structures',msa_dir=target/'processed/msa',
            constraints_dir=target/'processed/constraints',template_dir=target/'processed/templates',extra_mols_dir=target/'processed/mols')
        tokens=Boltz2Tokenizer().tokenize(data)
        molecules=dict(canonicals);molecules.update(data.extra_mols)
        molecules.update(load_molecules(cache/'mols',set(tokens.tokens['res_name'].tolist())-set(molecules)))
        features=module.Boltz2Featurizer().process(tokens,random=np.random.default_rng(42),molecules=molecules,
            training=False,max_seqs=64,max_atoms=None,max_tokens=None,pad_to_max_seqs=False,
            single_sequence_prop=0.0,compute_frames=True,compute_constraint_features=True)
        assert features['msa'].shape[0]==audits[-1]['rows']
        assert not features['msa_paired'][1:].any()
        assert not features['template_mask'].any()
        selected={k:features[k].numpy() for k in ['msa','msa_paired','profile','asym_id','res_type','template_mask','deletion_value']}
        np.savez_compressed(target/(record.id+'.model_features.npz'),**selected)
    assert len(audits)==8 and all(a['nonquery_multiple_chain_rows']==0 for a in audits)
    result=dict(passed=True,cells_audited=8,audits=audits,actual_processed_features_saved=True,templates_present=False,
        API_key_required=False,slurm_job_id=os.environ['SLURM_JOB_ID'])
    revisions=out/'code/revisions';revisions.mkdir(exist_ok=True)
    active={}
    for name in ['art_binding_run.py','art_binding_boltz_entry.py','art_binding_msa_control.py','art_binding_analyse.py','art_binding_feature_audit.py']:
        source=root/'scripts'/name;digest=hashlib.sha256(source.read_bytes()).hexdigest()
        target_file=revisions/(source.stem+'_'+digest[:12]+'.py');shutil.copy2(source,target_file)
        active[name]=dict(sha256=digest,snapshot=str(target_file.relative_to(out)))
    (out/'active_code.json').write_text(json.dumps(active,indent=2)+'\n')
    (target/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
