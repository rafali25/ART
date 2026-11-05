# Curated Jug phage protein-profile pilot

The additional Jug dataset produced **no reported ART RT, type I partner, type II
partner, exploratory two-LPJP1 partner, or generic RVT_1 domain hits** at the
recorded E <= 0.001 retrieval threshold. This result describes this panel and
these models; it does not establish absence of ART across all huge phages.

The panel contains 139 deposited curated genomes, 138 distinct DNA sequences,
52,152,061 total nucleotide positions, and 69,120 predicted proteins. One exact
genome duplication is preserved and identified in `duplicate_genome_groups.json`.
The archive is 14,678,419 bytes and matches the published MD5
`c80580eab632cebf86ba4f04f7be7559`. Original archive and extracted FASTA remain in
`sources/`; `archive_inventory.json` records filenames, CRCs and file hashes.
Extraction rejected absolute/traversing paths, symlinks and excess archive size.

Pyrodigal-gv 0.3.2 / Pyrodigal 3.7.1 ran in viral metagenomic mode with four
workers. All genomes selected `Gut_phage_code_11e`; the selected model and code
are recorded per genome. The caller's 153 partial predictions and ambiguous
translations were retained. All 69,120 proteins are compatible with their
source-coordinate DNA translations. Five exact differences arise because the
caller emits X for ambiguous codons that Biopython can resolve to a single amino
acid; there are zero unexplained translation discrepancies.

The five native HMMER searches used the existing profiles without retraining,
`--seed 42`, E/domain-report E <= 0.001 and explicit Z=domZ=69,120. Candidate
criteria remained sequence E <= 1e-5, domain E <= 1e-5 and profile coverage >=0.6.
The exploratory two-seed partner model remains separately labelled. Generic
RVT_1 evidence cannot supply an ART-family assignment. `commands.json` records
arguments, elapsed time and source-profile hashes; each raw `.tbl`, `.domtbl` and
log is in `searches/`. Empty candidate tables are intentional.

`verification.json` checks all genomic CDS mappings, protein hashes and GFF
coordinates. Separate known-seed controls recover positives with every profile
using the same effective database size. These controls validate executable/model
operation; they do not validate family specificity or independent-lineage recall.

Internal/TSV coordinates are zero-based and half-open. GFF3 coordinates are
one-based and inclusive. `gene_coordinates.tsv` records both the caller's genomic
begin/end partial flags and the strand-aware coding 5-prime/3-prime flags.
`genome_manifest.tsv`, `proteins.faa` and `predicted_genes.gff3` preserve the
DNA-to-protein link. No arrays were inferred or assessed, and no validated ART
locus claim is made.

The 139 curated Jug genomes form a related, predominantly gut-associated clade;
they are an additional dataset, not 139 independent evolutionary tests. This
negative pilot supports broadening the next panel across HPGC lineages rather
than treating genome size alone as evidence of ART. Full HPGC/JPGD download costs
and URLs are separately inventoried in `../dataset_inventory/`.

Rerun locally from the repository root:

```bash
.venv/bin/python results/art_sample_expansion_20261104/jug_screen/screen.py
.venv/bin/python results/art_sample_expansion_20261104/jug_screen/verify.py
```

Sources: [curated genome deposit, version 2](https://doi.org/10.6084/m9.figshare.29414588.v2),
[Chen et al. 2026](https://doi.org/10.1126/sciadv.aeb6265). This inventory and pilot
do not download either complete HPGC or JPGD collection.
