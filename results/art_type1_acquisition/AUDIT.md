# Local FASTA rescue audit and current type I sequences

**No additional type I candidate was recovered from the local genomes at the available thresholds.** This is a bounded retrieval result, not proof that all divergent type I genes are absent.

Six genomic FASTA files in `datasetss/` match previously searched nucleotide sequences exactly. The audit additionally includes the local GenBank-only Madawaska genome, for seven local genomes in total. LPJP1 is a separate positive retrieval control. The SA1 folder contains a Madawaska GenBank file and a different SA1 FASTA; the audit maps them by sequence, not folder name.

Across these eight genomes the audit searches 15,471 translated records:

| Translation source | Records |
| --- | ---: |
| Original single-interval CDS translations | 2,061 |
| Fresh Pyrodigal-gv predictions | 2,079 |
| Six-frame translated segments of at least 40 aa | 11,331 |

These are overlapping representations of the same DNA, not independent proteins or samples. Fifteen compound/missing-translation source CDS records are logged in `annotation_issues.tsv` and are excluded from this audit's simple-coordinate annotation layer. The existing broader annotated search is unchanged.

The four searched profiles are ART RT, type I partner, type II partner and generic Pfam RVT_1. Search-space size is the actual translated-record count. All reported hits at E <= 1e-3 are retained; the provisional passing rule is sequence and independent-domain E <= 1e-5 and profile coverage >=0.60.

All seven known local ART-like RT genes are recovered independently through fresh gene calls and six-frame translation. Neither layer produces an additional unmapped ART-profile candidate. No local translated record even reaches the type I reporting threshold. In the separate LPJP1 control, both RTs and both type I partners are recovered through all three translation layers. This verifies recovery of the available close controls; the current type I profile still has only one training seed.

The current LPJP1 package contains:

| RT and adjacent partner | Array status | Use |
| --- | --- | --- |
| QXN67953.1 / QXN67952.1 | Existing five-copy candidate; 898-nt array span | Main current type I array example |
| QXN67955.1 / QXN67954.1 | Not supported by the existing corrected assessment | Second homology query pair, array unresolved |

Both pairs are on MZ422438.1, on the reverse strand. They remain one source-genome group. The exported sequences are the original deposited translations, not alternative starts selected by the new gene caller. This task does not add a new type I sample.

Files in `current_typeI/`:

- `RT_proteins.faa`, `partners.faa`, `NTDs.faa`: two sequences each.
- `CDS.fna`: four source-verified CDSs in transcriptional orientation.
- `upstream.fna`: two 6-kb upstream regions.
- `array_candidates.fna`: the one existing candidate array.
- `repeat_spacer_repeat_units.fna`: four DNA units; not measured RNA boundaries.
- `pairs.tsv`, `locus_evidence.tsv`, `repeats.tsv`, `spacers.tsv`.
- Individual oriented `.gb`, `.fna`, `.gff3` and `.map.json` files.

The array status, copy coordinates and old cohort-adjusted probability are inherited from `results/art_arrays_assessed/`; no new array search, null simulation, model training or cohort correction was performed. The exported array retains the original detector's limitations, including uncertain full repeat boundaries. It is a computational candidate, not an experimentally validated type I system.

`protein_to_contig.tsv`, `all_reported_profile_hits.tsv`, `passing_profile_hits.tsv`, and `recovery_by_genome.tsv` preserve the new search evidence. `local_FASTA_inventory.tsv` records duplicate DNA and file hashes. `source_records.tsv` records source hashes, gene-calling models and translation tables. Every six-frame fragment was retranslated from its mapped source interval. Original source files and assessments were not edited.

The linear six-frame audit does not reconstruct genes spanning circular sequence boundaries. Gene calls and stop-to-stop fragments are retrieval hypotheses. Fresh predictions do not overwrite deposited annotations.

The proposed next acquisition and extraction workflow is in [METHOD.md](METHOD.md).
