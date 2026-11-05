"""Explicit non-query chain separation for the registered unpaired pilot input."""
import numpy as np


def separate_chains(msa, deletion, paired, chain_ids, gap_token):
    """Keep query column; split homolog columns into chain-specific columns.

    Inputs use Boltz's pre-transpose [tokens, rows] layout. Residue/deletion
    content for each homolog is unchanged; other chains are explicitly gapped.
    No inferred taxonomy, locus correspondence, or cross-chain pairing survives.
    """
    msa=np.asarray(msa);deletion=np.asarray(deletion);paired=np.asarray(paired)
    chain_ids=np.asarray(chain_ids)
    if msa.ndim!=2 or deletion.shape!=msa.shape or paired.shape!=msa.shape or len(chain_ids)!=msa.shape[0]:
        raise ValueError('Incompatible MSA/deletion/pairing/chain shapes')
    if msa.shape[1]<1:raise ValueError('Query row is required')
    columns=[msa[:,0].copy()];deletions=[deletion[:,0].copy()];flags=[paired[:,0].copy()]
    for chain in np.unique(chain_ids):
        selected=chain_ids==chain
        for i in range(1,msa.shape[1]):
            if not np.any(msa[selected,i]!=gap_token):continue
            row=np.full(msa.shape[0],gap_token,dtype=msa.dtype)
            d=np.zeros(msa.shape[0],dtype=deletion.dtype)
            row[selected]=msa[selected,i];d[selected]=deletion[selected,i]
            columns.append(row);deletions.append(d);flags.append(np.zeros(msa.shape[0],dtype=paired.dtype))
    return np.stack(columns,axis=1),np.stack(deletions,axis=1),np.stack(flags,axis=1)


def audit_msa(msa,paired,chain_ids,gap_token,condition):
    msa=np.asarray(msa);paired=np.asarray(paired);chain_ids=np.asarray(chain_ids)
    occupancy=[sum(np.any(msa[chain_ids==c,i]!=gap_token) for c in np.unique(chain_ids)) for i in range(msa.shape[1])]
    result=dict(condition=condition,rows=int(msa.shape[1]),tokens=int(msa.shape[0]),
        query_occupied_chains=int(occupancy[0]),nonquery_multiple_chain_rows=int(sum(v>1 for v in occupancy[1:])),
        nonquery_paired_feature_values=int(np.count_nonzero(paired[:,1:])),
        nonquery_rows=int(sum(v>0 for v in occupancy[1:])))
    if result['nonquery_multiple_chain_rows'] or result['nonquery_paired_feature_values']:
        raise ValueError('Cross-chain homolog information remains: '+str(result))
    if condition=='S' and result['nonquery_rows']:
        raise ValueError('Single-sequence condition contains homologs')
    if condition=='U' and not result['nonquery_rows']:
        raise ValueError('Unpaired condition silently became single-sequence')
    return result
