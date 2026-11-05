"""Summarize preserved source checks without inferring missing protein identities.

The online responses and the small, checksum-verified related data archive were
retrieved separately. This script records their scope and the unresolved inputs;
it does not substitute internal figure labels for accession identifiers.
"""
from pathlib import Path
import csv
import hashlib
import json
import re
import subprocess

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/art_discovery/typeIII_source_audit'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    metadata=json.loads((OUT/'biorxiv_api.response').read_text())['collection']
    assert metadata and metadata[0]['doi']=='10.64898/2026.10.05.755809'
    online=(OUT/'biorxiv_full.txt').read_text(encoding='utf-8')
    assert 'the 8 and 6 database members of the two type III genes' in online
    related=json.loads((OUT/'ai3_zenodo_metadata.json').read_text())
    expected=next(line.split()[0] for line in (OUT/'ai3_CHECKSUMS_TOPLEVEL.sha256').read_text().splitlines()
                  if line.split()[-1].lstrip('*')=='core.tar.gz')
    assert sha(OUT/'ai3_core.tar.gz')==expected
    with (OUT/'ai3_sequences__sequences.tsv').open(encoding='utf-8') as f:
        catalogue=list(csv.DictReader(f,delimiter='\t'))
    lpjp=[r for r in catalogue if r['name'] in ('P1_QXN67954','P2_QXN67952')]
    assert len(lpjp)==2 and all('Type I' in r['role'] for r in lpjp)
    assert not [r for r in catalogue if re.search(r'type[ _-]*iii',r['role'],re.I)]
    result=dict(
        audit_date='2026-10-19',
        outcome='verified_type_III_seeds_not_retrieved_from_checked_sources',
        original_publication_doi='10.64898/2026.10.05.755809',
        publication_versions_returned=[r['version'] for r in metadata],
        paper_type_III_labels=['ART_30','ART_31','ART_32 (L10)'],
        figure_evidence='Supplementary Figure 3, PDF p.17; dark-purple type III labels visually inspected.',
        original_profile_seed_counts={'Type_III_first_partner':8,'Type_III_second_partner':6},
        original_census_loci_per_family=10,
        method_evidence='Database-wide occurrence of partner families, PDF p.37 and archived bioRxiv full text.',
        original_sequence_accession_mapping_retrieved=False,
        original_type_III_FASTA_or_HMM_retrieved=False,
        related_deposit_doi='10.5281/zenodo.23086637',
        related_deposit_title=related['metadata']['title'],
        related_deposit_core_sha256=sha(OUT/'ai3_core.tar.gz'),
        related_deposit_scope='Metadata, README, checksums, core archive file inventory, sequence catalogue and small text tables. Large structure/alignment archives were not downloaded.',
        related_LPJP1_partners=[dict(accession=r['accession'],role=r['role']) for r in lpjp],
        related_core_explicit_type_III_catalogue_entries=0,
        related_code_repository='https://github.com/agent-next/art-discovery-repro',
        related_code_commit=subprocess.check_output(['git','-C',str(ROOT/'.tools/art-discovery-repro'),'rev-parse','HEAD'],text=True).strip(),
        acquired_verified_type_III_sequences=0,
        profiles_built={'Type_III_first_partner':False,'Type_III_second_partner':False},
        required_inputs=['Public contig/protein accession mapping for ART_30, ART_31 and ART_32/L10',
                         'Separate FASTA alignments/sequences or HMMs for the first and second type III partner families',
                         'Genomic coordinates and order to verify partner-family assignments'],
        limits='This is an audit of the listed sources, not proof that no type III sequence exists elsewhere. An inaccessible source is not evidence of data absence. Neither size nor helical structure identifies these families.'
    )
    (OUT/'findings.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    report='''# Type III partner seed source audit

Checked on 19 October 2026. No verified type III first- or second-partner sequences were retrieved from the sources inspected. Both profiles remain unavailable; no candidate sequence was relabelled as a verified seed.

## What the source paper provides

- Supplementary Figure 3 (PDF p. 17) marks ART_30, ART_31 and ART_32 (also labelled L10) as type III. The figure's dark-purple key was visually checked. These are internal study labels, not nucleotide or protein accessions.
- The database-wide partner census methods (PDF p. 37) state that the two type III profiles used 8 and 6 database-member sequences, respectively. Each family returned 10 genomic loci. These counts do not supply the actual sequences or independent lineages.
- The type III structural example in Figure 3 is labelled “locus A”. The protein is described as approximately 170 aa, with no clear database-domain homolog. Length and predicted helicity cannot establish family membership.
- The supplied PDF has no embedded attachments. Its hyperlinks did not identify a type III sequence deposit or accession table.

## Online sources checked

| Source | Scope and result |
| --- | --- |
| [bioRxiv publication](https://www.biorxiv.org/content/10.64898/2026.10.05.755809v1.full) and [API metadata](https://api.biorxiv.org/details/biorxiv/10.64898/2026.10.05.755809) | Version 1 was returned. Full text was archived. No type III FASTA/HMM download or label-to-accession mapping was located. |
| Requested bioRxiv supplementary-material endpoint | Returned an article page without supplementary sequence-file links. This does not prove that no supplementary data exist elsewhere. |
| JATS source XML link reported by the API | Returned HTTP 403; an alternate versioned source-XML endpoint returned 404. These are access failures, not evidence of data absence. |
| Crossref and Europe PMC records | Metadata checked; no usable type III sequence mapping was recovered. |
| DataCite and Zenodo title-keyword searches | Exact-title keyword queries returned no records. These queries alone are not an exhaustive repository search. |
| [AI3 Discovery's related ART study](https://ai3-discovery.onpod.ai/research/art-partner-structure-alignment-depth) and [Zenodo deposit](https://doi.org/10.5281/zenodo.23086637) | The small core archive, sequence catalogue, README, checksums and tables were inspected. Its LPJP1 P1 (QXN67954.1) and P2 (QXN67952.1) are explicitly type I partners of two RT paralogs. There is no explicitly assigned type III seed set in the inspected catalogue/tables. The core archive checksum matched the deposit. Large prediction archives were not downloaded. |
| [Independent reproduction repository](https://github.com/agent-next/art-discovery-repro) | Source and recorded public-data runs checked at the commit recorded in findings.json. No traceable type III first/second partner seed collection was located. This is a third-party reproduction, not the original authors' data release. |

The related Zenodo record is useful evidence for type I partner sequences, but its P1/P2 naming must not be confused with the two distinct type III genes.

## Data needed to activate the profiles

The unresolved dependency is a mapping from ART_30, ART_31 and ART_32/L10 to public contigs/proteins, or the source study's separate 8-sequence and 6-sequence type III FASTAs/alignments/HMMs. Coordinates and gene order are needed to keep the first and second families correctly assigned. The pipeline already has separate configuration slots and route tracking for both.

This search does not establish that the sequences are globally unavailable or that the source authors have never released them. It establishes that a verified seed set was not obtained from the listed sources. No emails or other messages were sent to authors or third parties.

Raw responses, repository/deposit identifiers, archive inventories and checksums accompany this report. `findings.json` records the structured result. The existing candidate-locus table and trained profiles were not changed by this audit.
'''
    (OUT/'REPORT.md').write_text(report,encoding='utf-8')
    hashes={p.relative_to(OUT).as_posix():sha(p) for p in OUT.rglob('*') if p.is_file() and p.name!='audit_sha256.json'}
    (OUT/'audit_sha256.json').write_text(json.dumps(hashes,indent=2),encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['outcome','paper_type_III_labels','original_profile_seed_counts','acquired_verified_type_III_sequences','profiles_built']},indent=2))


if __name__=='__main__':main()
