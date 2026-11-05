# Type III partner seed source audit

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
