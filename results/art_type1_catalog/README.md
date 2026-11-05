# Combined provisional type I ART sample catalog

**13 primary loci: 11 new public candidates plus the two LPJP1 controls.** The primary partners form 12 clusters at 70% identity and 80% mutual coverage. Different contigs and clusters are not proof of independent biological samples.

10 pairs have complete, translation-matched CDSs; 3 require gene-boundary/completeness review. `complete_pair_representatives.tsv` selects 9 complete pairs, one per eligible partner cluster.

Array evidence is corrected across **all 17 distinct RT candidate loci** tested in either cohort, including candidates without type I partners. This is a conditional comparison within the selected panels, not a database-wide false discovery rate. Duplicate LPJP1 controls are counted once.

| Locus | Copies | Spacer gaps, nt | Combined corrected p | Flags |
|---|---:|---|---:|---|
| ERR8055547_17635__ERR8055547_17635_pgv12 | 5 | [284, 326, 630, 289] | 0.0084958 |  |
| MZ422438.1__QXN67953.1 | 5 | [148, 264, 137, 199] | 0.0084958 |  |
| MGYC002687634053__MGYC002687634053_pgv52 | 5 | [219, 226, 210, 689] | 0.0084958 | gap_may_span_multiple_missing_repeats |
| MGYC003526048914__MGYC003526048914_pgv179 | 3 | [122, 121] | 0.016992 |  |

Statistical support establishes a repeat tract, not ART function or correct full repeat boundaries. The largest gap in MGYC002687634053 may span more than one missing copy despite the permissive seed-chain fit; it is explicitly flagged for review. Two-copy and weak tracts remain exploratory.

The type I profile still lacks verified independent positive holdout groups. NART comparisons and independent full-boundary validation remain incomplete. Partial genes are excluded from the complete-pair shortlist. Short upstream windows do not establish array absence; complete protein pairs may still be retained when their upstream DNA is incomplete.

Use `typeI_samples.tsv` for the complete inventory, `complete_pair_representatives.tsv` for a diversity-oriented starting selection, and `RT_proteins.faa`, `primary_partners.faa`, `NTDs.faa` for sequences. Each row links through its `result_folder` and `locus_id` to oriented DNA, GenBank/GFF annotations, coordinate maps, CDS DNA, arrays, repeat alignments and spacers.

`supported_array_candidates.fna` contains the four supported candidate tracts (three new plus LPJP1). `supported_array_candidates.tsv` links each to its coordinate map and locus annotation. `array_repeats.fna` and `array_spacers.fna` contain the inferred components. These are DNA extractions; transcription start/end sites and mature RNA boundaries have not been inferred.

The two original result folders retain their own cohort corrections. Source manifests, exact genomic translation checks, raw HTTP responses, database releases, ETags and SHA-256 checksums are preserved. MGnify retrieval was capped at ten contigs; remaining hits are deferred, not negative.

Rebuild with `.venv/Scripts/python.exe scripts/build_art_type1_catalog.py` after finalizing both cohorts.
