# RT-oriented locus extraction

Extracted 9 candidate RT loci from 8 source records, including 3 reverse-strand RTs. Every locus retains its full annotated RT and named partner CDSs. Requested flanks: 6,000 nt upstream and 6,000 nt downstream. Discovery classifications are inherited from the input candidate table and remain provisional.

## Files

- `oriented_loci.fna` and `oriented_loci.gb`: all oriented loci.
- `loci/`: individual FASTA, GenBank, GFF3 and JSON coordinate maps.
- `coordinate_map.tsv`: exact source/locus mapping segments, in locus order.
- `locus_manifest.tsv`: sequence paths, RT position, partners, topology, flank availability and boundary flags.
- `gene_annotations.tsv`: gene/CDS coordinates, original feature index, source location, strand and qualifiers.
- `sequence_and_annotation_issues.tsv`: ambiguous sequence, explicit gaps, annotation-free intervals and clipped features.
- `original_sources/`: unchanged source files and indexed annotation snapshots.
- `RT_review.json`, `verification.json`, `summary.json`, `output_sha256.json`: QC and provenance.

## Coordinate convention

All TSV/JSON intervals are zero-based, half-open. GenBank/GFF3 are one-based, inclusive. On a mapping segment [s,e) with locus offset o, [a,b) maps to [o+a-s,o+b-s) for orientation +1, or [o+e-b,o+e-a) for orientation -1. To map a locus interval that crosses a segment boundary, split it at that boundary first. Strand in the map is the RT orientation relative to the original contig. Mapped gene strands are relative to the oriented locus.

## Current results

There are 0 ambiguous nucleotides and 0 wrapped loci. Full requested flanks are available for 9/9 upstream and 9/9 downstream intervals. Circular wrapping is enabled only by explicit circular topology in the source record; unknown/linear topology never enables wrapping.

| Locus | Original RT strand | Length (nt) |
|---|---:|---:|
| LC680885.1__BDE75684.1 | +1 | 13503 |
| MN091626.1__QDJ97640.1 | +1 | 13503 |
| MW218148.1__QPI16926.1 | +1 | 13506 |
| MW248466.1__QQM14740.1 | +1 | 13506 |
| MW349129.1__QQO92858.1 | +1 | 13506 |
| MZ422438.1__QXN67953.1 | -1 | 13365 |
| MZ422438.1__QXN67955.1 | -1 | 13530 |
| MZ779063.1__UAJ17040.1 | +1 | 13503 |
| OR836606.1__WRM43423.1 | -1 | 13503 |

9/9 RT translations pass checks under the annotated genetic code. Review flags:

- QXN67955.1: RT_profile_coverage_below_review_threshold_0.9. Source translation verified: True.

Review flags retain the source CDS for inspection; profile coverage alone does not establish biological truncation. See `RT_review.json` for all sequence/annotation checks.

All source files are unchanged. Every extracted base, fully retained feature sequence, complete RT and partner, GenBank location round trip and GFF coordinate conversion passed verification.

## Annotation limits

All current sources are annotated; no gene caller was used. Unannotated viral contigs require candidate CDS calls and profile discovery before this extractor can map their RTs. Pyrodigal-gv is an option: https://github.com/althonos/pyrodigal-gv . Keep caller/model/genetic-code provenance and review suspected split RTs. This step preserves the provided gene boundaries and does not repair gene calls.

Annotation-free intervals describe annotation coverage; they may contain noncoding ART arrays and are not inferred assembly gaps. Context annotations cut by a locus boundary are marked clipped; full original qualifiers and translations remain in the archived source. Extracted clipped CDSs omit the full-length translation and adjust codon_start for the retained reading frame.

Implementation uses Biopython 1.88 plus custom segment/interval handling. See ../../README.md and ../../scripts/extract_art_loci.py for commands and parameters.
