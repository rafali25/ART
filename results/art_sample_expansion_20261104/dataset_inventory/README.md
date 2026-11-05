# HPGC and JPGD acquisition inventory

Public REST metadata and HTTP HEAD checks, recorded for the 2026-11-04 search.
No genome or proteome archive was downloaded by this inventory step. The saved
metadata, the small published JPGD analysis recipe, and HEAD results total 24,088
bytes. Source responses and SHA-256 hashes are in `inventory_manifest.json`;
`files.tsv` contains every file advertised by the two deposits.

| Resource | Records / scope | Archive | Compressed bytes | Access |
| --- | --- | --- | ---: | --- |
| Curated Jug genomes, HPGC project | 139 curated genomes; article v2 | `139.curated.Mac.phage.genomes.fa.zip` | 14,678,419 | HEAD 200; length confirmed |
| HPGC | 7,230 nonredundant huge-phage genomes; article v1 | `HPGC.7230.genomes.final.fasta.gz` | 560,950,249 | HEAD 200; length confirmed |
| JPGD | Paper describes 10,754 jumbo-phage genomes; Zenodo version 2 | `Jumbo-phage.all.fasta.gz` | 871,869,424 | HEAD 200; API supplies length |

The curated Jug subset is the cheapest concrete next search panel. It fits the
existing 20 MB per-source acquisition budget and supplies long genomic context
from a different phage lineage. Its 139 genomes include related genomes and do
not represent 139 independent lineages. The filenames retain the older name
"Mac" although the current deposit and paper use "Jug".

The APIs advertise nucleotide genomes but no protein FASTA, GFF coordinates,
prebuilt sequence-search database, or hosted HMM search endpoint. A bounded next
step is to checksum the 14.7 MB curated archive, retain deposited nucleotide
identifiers, call candidate proteins with the repository's viral gene caller,
retain DNA/CDS/strand/genetic-code mappings, and screen RT plus partner profiles.
Pair evidence and array assessment must remain separate from taxonomy. Confirm
record counts and sequence hashes only after acquisition; this inventory did not
inspect sequence contents.

Before increasing download scope, deduplicate each candidate genome against the
current catalog and previous search inputs. HPGC and JPGD compile material from
other public resources, so overlap between these catalogs and IMG/VR, UHGV,
Logan, or NCBI is expected. New database identifiers alone do not establish new
independent samples. The two full gzip archives provide no separately listed
small representative download in the inspected deposit versions.

The JPGD deposit also supplies `megaphage_analysis.sh` (9,655 bytes). It is a
published analysis recipe with Markdown fences and references to external
helpers/databases, rather than a self-contained executable workflow. It was
saved for inspection and was not executed. Fresh gene calls should use the
repository's recorded tools rather than silently reuse this recipe's defaults.

Sources:

- HPGC project API: https://api.figshare.com/v2/projects/254252/articles
- HPGC genome record: https://figshare.com/articles/dataset/HPGC_genomes/30506864
- Curated Jug record: https://figshare.com/articles/dataset/Jug_phages_curated_genomes/29414588
- HPGC paper: https://doi.org/10.1126/sciadv.aeb6265
- JPGD record/API: https://zenodo.org/records/14974448 and https://zenodo.org/api/records/14974448
- JPGD paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC13037316/
