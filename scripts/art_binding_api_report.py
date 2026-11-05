"""Verify downloaded API models and export a portable pilot result table."""
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import statistics
import zipfile

from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1

from art_binding_api_credentials import credential


def main():
    if not os.environ.get('SLURM_JOB_ID'):raise SystemExit('Use Slurm.')
    root=Path(__file__).resolve().parents[1];out=root/'data/ART binding/api_pilot_20261105'
    key,provider,mode=credential(root)
    summary=json.loads((out/'summary.json').read_text());assert summary['status']=='complete' and summary['completed']==4
    rows=[];models=0
    for result in sorted(summary['results'],key=lambda r:r['cell_id']):
        folder=out/result['cell_id'];request=json.loads((folder/'request.json').read_text())
        structure_path=root/result['structure_file']
        model=MMCIFParser(QUIET=True,auth_chains=False,auth_residues=False).get_structure('model',structure_path)[0]
        for chain,polymer in zip(model.get_chains(),request['polymers']):
            observed=''.join(seq1(res.resname) for res in chain.get_residues())
            assert observed==polymer['sequence'],result['cell_id']
            align=polymer['msa']['query_only']['a3m']['alignment']
            assert align=='>query\n'+polymer['sequence']+'\n'
        assert request['sampling_steps']==200 and request['recycling_steps']==3
        assert request['diffusion_samples']==1 and request['without_potentials'] is True
        assert not request['constraints'] and not request['ligands']
        assert hashlib.sha256(structure_path.read_bytes()).hexdigest()==result['structure_sha256']
        models+=1
        rows.append(dict(cell_id=result['cell_id'],RT_sample_id=result['RT_sample_id'],partner_sample_id=result['partner_sample_id'],
            native_pair=result['native_pair'],condition=result['condition'],protein_iptm=result['protein_iptm'],
            ptm=result['ptm'][0],confidence=result['confidence'][0],contact_residue_pairs_5A=result['contacts_5A_residue_pairs'],
            client_seconds=result['client_elapsed_seconds'],server_seconds=result['server_runtime_metrics'].get('total_time_seconds'),
            structure_file=result['cell_id']+'/model_0.cif',sampling_steps=200,recycling_steps=3,
            exact_seed_control=False,full_PAE_available=result['full_PAE_available'],confirmatory=False))
    for name,delimiter in [('results.csv',','),('results.tsv','\t')]:
        with (out/name).open('w',newline='') as f:
            writer=csv.DictWriter(f,list(rows[0]),delimiter=delimiter);writer.writeheader();writer.writerows(rows)
    native=[r['protein_iptm'] for r in rows if r['native_pair']]
    alternative=[r['protein_iptm'] for r in rows if not r['native_pair']]
    timing=[r['client_seconds'] for r in rows]
    report=dict(predictions=4,verified_structures=4,protein_CIF_sequences_verified=8,
        client_seconds_min=min(timing),client_seconds_median=statistics.median(timing),client_seconds_max=max(timing),
        server_seconds_median=statistics.median(r['server_seconds'] for r in rows),
        native_protein_iptm_mean=statistics.mean(native),alternative_protein_iptm_mean=statistics.mean(alternative),
        descriptive_native_minus_alternative=statistics.mean(native)-statistics.mean(alternative),
        conclusion='Technical API pilot only; one unseeded sample per cell with closely related references cannot establish binding or specificity.',
        exact_seed_control=False,processed_feature_access=False,full_PAE_available=False,confirmatory=False,
        API_provider=provider,cost_checks_performed=False,slurm_job_id=os.environ['SLURM_JOB_ID'])
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    (out/'README.md').write_text(f'''# Hosted Boltz-2 development pilot — 5 November 2026

The user-provided key was used successfully with NVIDIA's hosted Boltz-2 service.
All **four query-only native/alternative pair predictions** completed, with three
recycling steps, 200 sampling steps, and one diffusion sample per request.
End-to-end latency was **{min(timing):.1f}–{max(timing):.1f} seconds per request**;
two requests ran concurrently. Median server-reported runtime was
**{report['server_seconds_median']:.1f} seconds**. This is a measured pilot result,
not a guarantee for larger complexes, throughput, or future service queues.

Use `results.csv` or `results.tsv` for the four-cell table. Each cell directory
contains its exact request, full response, runtime/status records, parsed result,
and downloaded `model_0.cif`. Both chain sequences in every CIF were checked
against the submitted proteins. The credential is not stored in this dataset.

## Scope and interpretation

This is the SA1/MarsHill 2×2 type II development panel, not the full type I study.
MSAs contain only the query for each chain. No template, interaction restraint,
ligand, or affinity prediction was supplied. The four results provide technical
validation of remote execution and output retrieval; the pairs are close
relatives and alternative pairings have no experimental non-binding label.

The API endpoint does not expose fixed-seed control or the processed features
needed to verify our strict unpaired-MSA adapter. Consequently only the S/query
condition was submitted; these are not the registered fixed-seed GPU replicates.
Full PAE matrices are not available in the hosted JSON response. Local GPU/CPU
outputs are retained separately. Backend calibration is required before making
biological conclusions. Low scores, a flat matrix, or a small descriptive contrast
cannot establish absence of interaction or a prediction artefact.

The original slow CPU probe reached its 20-minute reduced-settings budget.
The local GPU smoke test separately completed two models. The queued remaining
local matrix was cancelled once the hosted service worked, preserving completed
outputs. API cost checks were omitted as instructed by the user.

References: [NVIDIA hosted endpoint](https://docs.api.nvidia.com/nim/reference/mit-boltz2-infer),
[model card](https://docs.api.nvidia.com/nim/reference/mit-boltz2), and
[NIM input/output documentation](https://docs.nvidia.com/nim/bionemo/boltz2/inference).
The model card identifies Boltz2 v2.2.1; an exact deployed checkpoint hash is not
exposed by the hosted response and is not claimed here.
''')
    hashes={str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in out.rglob('*') if p.is_file() and p.name not in ['SHA256SUMS.json','verification.json']}
    for relative in hashes:
        assert key.encode() not in (out/relative).read_bytes(),'Credential unexpectedly present in an export file; content withheld.'
    (out/'SHA256SUMS.json').write_text(json.dumps(hashes,indent=2)+'\n')
    verification=dict(passed=True,structures_verified=models,protein_sequences_verified=8,request_settings_verified=True,
        output_files_hashed=len(hashes),credential_value_absent_from_export=True,slurm_job_id=os.environ['SLURM_JOB_ID'])
    (out/'verification.json').write_text(json.dumps(verification,indent=2)+'\n')
    archive=out.with_suffix('.zip')
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for path in sorted(out.rglob('*')):
            if path.is_file():z.write(path,str(Path(out.name)/path.relative_to(out)))
    with zipfile.ZipFile(archive) as z:assert z.testzip() is None
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
