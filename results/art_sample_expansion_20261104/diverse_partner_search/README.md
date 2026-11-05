# Bounded diverse-partner public search, 4 November 2026

Three live public EBI phmmer searches produced **two incremental MGnify peptide
accessions and two newly retrieved DNA contigs**. Native screening found one
complete translated generic-RT/type-I-partner context pair. It closely resembles
an existing catalog pair and its RT fails the frozen ART-profile coverage screen.
This branch establishes **no new supported ART lineage**.

## Queries and search scope

All query partners and their cognate RTs were retranslated exactly from source DNA
after verifying the recorded source hashes. The queries occupy three different
existing 70%-identity partner clusters; no new profile was trained.

| Query source | Partner length | Returned hits | Incremental MGYP accessions |
| --- | ---: | ---: | ---: |
| SRR8925777_19723 | 594 aa | 10 | 1 |
| ERR8055547_17635 | 890 aa | 11 | 2 |
| MGYC002687634053 | 587 aa | 9 | 0 |

The live database registry confirmed `mgnify30_c2`, **release 2026_07**, before
submission. This database is MGnify30-C2 non-singletons. Searches used reporting
E/domain E <= 0.001 and requested at most 50 hits per query. Each completed result
had one page; the union was 30 rows and 11 MGYP accessions. Two were absent from
the union of the original public and MGnify saved hit tables. This is a comparison
of these specific recorded searches, not an exhaustive census of MGnify.

The MGYC seed minimized maximum full-length global identity to both LPJP1
controls and the first two public queries among the available complete MGYC
representatives. Alignment parameters and comparisons are recorded in
`search.py` and `seed_diversity.tsv`.

Job IDs, complete query sequences, database registry, raw requests and responses,
intermediate polls and final results remain in `jobs.json`, `queries.json`,
`database.json` and `remote/`.

## Incremental sequence retrieval

| Incremental peptide | Best reported E | Retrieved contig | DNA length | Exact peptide mapping, 0-based half-open |
| --- | ---: | --- | ---: | --- |
| MGYP011007627867 | 7.69e-8 | MGYC003966598318 | 252,959 nt | [49043,49739), +; 232 aa |
| MGYP003591901103 | 7.34e-8 | MGYC001042694086 | 1,171 nt | [11,914), -; 301 aa |

Both peptides map uniquely by exact genomic translation. Their API full-length
flags describe those peptides, and do not establish a complete ART partner or RT.
Two protein metadata records and two source-contig requests used 254,236 source
bytes, below the caps of five metadata proteins, three contigs and 5 MB.
The short reverse-strand hit has only 257 nt upstream in its transcription
direction, preventing an intact upstream RT assessment.

## Native screening and novelty check

`screen.py` uses native HMMER and MAFFT, viral gene calls and six-frame rescue
through the existing pipeline. It reuses frozen ART RT, type-I single-seed,
type-I two-LPJP1 exploratory, type-II and RVT_1 profiles. No acceptance threshold
or active profile changed. The two contigs yielded 1,729 annotated-protein/ORF
representations.

MGYC003966598318 contains one complete translated, co-oriented adjacent pair:

- RT `MGYC003966598318_pgv45`: **609 aa**, [52908,54738), +.
- Type-I partner `MGYC003966598318_pgv46`: **587 aa**, [54724,56488), +.
- Partner profile: E=2.2e-36, coverage **67.6%**; compatible nearest-gene context.
- ART RT profile: E=8.4e-28, 94.2 bits, coverage **37.2%**, below the unchanged
  required 60%. RVT_1 coverage is 82.5%. Retrieval was partner-first with generic
  RT evidence; ART-family placement remains unresolved.
- Four recovered repeat copies with spacers 229, 219 and 538 nt are exploratory.
  No new shuffled-null assessment or supported-array claim was made.

Full-length global identities to the existing MGYC002687634053 pair are **98.03%
for RT and 97.44% for partner**. Neither protein is an exact duplicate, but this
does not provide a new diverse partner representative or establish an independent
lineage. Comparisons against every prior primary pair are preserved in
`prior_catalog_similarity.tsv`. New MGYC003966598318 and prior MGYC002687634053
belong to the same study, ERP150573.

The short MGYC001042694086 contig produced no mapped RT locus under the frozen
screens. Lack of sufficient context remains explicit.

## Files and verification

- `remote_hits.tsv`, `search_summary.json`: query hits and incremental counts.
- `protein_metadata.tsv`, `ranked_locations.json`: acquisition evidence.
- `acquired_source_records.tsv`, `homolog_to_contig.tsv`, `sources/`: source DNA,
  provenance and exact peptide maps.
- `processed/`: source annotations, gene calls, raw HMMER searches, candidate
  loci, provisional pair sequences, source mappings and exploratory arrays.
- `verification.json`: successful live jobs, budget checks, input/source hashes
  and replay of the complete pair's two CDS translations.
- `input_sha256.json`, `output_sha256.json`: frozen input references and artifact
  checksums. Active profile training and earlier result folders were untouched.

Reproduction commands, from the repository root:

```bash
.venv/bin/python results/art_sample_expansion_20261104/diverse_partner_search/search.py prepare
.venv/bin/python results/art_sample_expansion_20261104/diverse_partner_search/search.py submit
.venv/bin/python results/art_sample_expansion_20261104/diverse_partner_search/search.py poll
.venv/bin/python results/art_sample_expansion_20261104/diverse_partner_search/search.py collect
.venv/bin/python results/art_sample_expansion_20261104/diverse_partner_search/acquire.py
.venv/bin/python results/art_sample_expansion_20261104/diverse_partner_search/screen.py
.venv/bin/python results/art_sample_expansion_20261104/diverse_partner_search/finalize.py
```

Completed submission responses and successful job states are cached. Start a new
output directory and explicitly freeze its database release for a genuinely new
search rather than reusing these historical jobs.
