# MGnify follow-up to the Type I ART pilot

Processed 10 public contigs plus LPJP1. Recovered 7 provisional public primary type I RT–partner loci. All RT candidates (13, including controls and RTs without type I partners) belong to the declared array-testing cohort. These are candidate loci, not validated ART systems.

| Contig | Primary partner | Upstream nt | Copies | Spacer gaps nt | Corrected p | Assessment |
|---|---|---:|---:|---|---:|---|
| MGYC001394629236 | MGYC001394629236_pgv5 | 1 | 0 | [] | 1 | incomplete_upstream_assessment |
| MGYC003138808455 | MGYC003138808455_pgv3 | 2 | 0 | [] | 1 | incomplete_upstream_assessment |
| MGYC003138808455 | MGYC003138808455_pgv1 | 3724 | 0 | [] | 1 | incomplete_upstream_assessment |
| MGYC002687634053 | MGYC002687634053_pgv53 | 6000 | 5 | [219, 226, 210, 689] | 0.006496751624187906 | statistically_supported_array_requires_copy_and_family_review |
| MZ422438.1 | QXN67952.1 | 6000 | 5 | [148, 264, 137, 199] | 0.006496751624187906 | statistically_supported_array_requires_copy_and_family_review |
| MGYC001394626043 | MGYC001394626043_pgv10 | 6000 | 0 | [] | 1 | not_supported_by_new_cohort_test |
| MGYC003526048914 | MGYC003526048914_pgv180 | 6000 | 3 | [122, 121] | 0.012993503248375811 | statistically_supported_array_requires_copy_and_family_review |
| MGYC002345864168 | MGYC002345864168_pgv3 | 5791 | 0 | [] | 0.06467661691542288 | incomplete_upstream_assessment |
| MZ422438.1 | QXN67954.1 | 6000 | 0 | [] | 0.32338308457711445 | not_supported_by_new_cohort_test |

The profile searches cover both RT-first and type I partner-first routes. Original annotations are retained; unannotated DNA is called with Pyrodigal-gv and searched in six frames. Every downloaded query homolog has a recorded exact translated DNA mapping. MGYC identifiers are release snapshots with sequence checksums.

Both 3 kb and 6 kb array searches run independently. The 200-draw screen can only stop for non-support; positive statistical status requires 2,000 local-block and 2,000 dinucleotide shuffles with complete seed rediscovery and chain selection. The larger p is corrected across all candidate loci in this cohort. Short upstream windows remain explicitly incomplete even when they contain a repeat tract.

The separate RT-domain tree and competing myRT HMM searches do not use array or partner architecture. ART/NART separation remains unresolved without verified NART references. Type I profiles have one training-genome group and no independent verified positive holdout group. The five GNAT negatives and archived protein cross-screen are inherited from the pilot, not a new validation experiment. Known repeat-copy precision and full-boundary validation limitations still apply.

Acquisition is bounded by the recorded config and source manifest: the first ten retrievable contigs from query-balanced ranked hits were retained. Unretrieved hits are deferred, not negative. The four MGnify searches produced 418 collected query-hit rows (293 unique MGYP proteins); both RT queries were capped at 200 results. Each partner query returned nine hits. NCBI BLAST was still queued at the latest poll. This is not a prevalence estimate.

Use `primary_typeI_pairs.tsv`, `candidate_loci.tsv`, `background_tests.tsv`, `loci/`, `primary_typeI_partners.faa`, `candidate_RT_proteins.faa`, `NTDs.faa`, `primary_pair_CDS.fna`, `upstream_sequences.fna`, `repeats.tsv`, `spacers.tsv`, and `verification.json`. Coordinates are zero-based/half-open internally. Spacer intervals flagged `possible_missing_repeat` may span an undetected copy and should not be interpreted as a confirmed single biological spacer.

Run the same sequential processing commands as the pilot with `--config config/art_type1_mgnify.json`. The original pilot results and its six-locus statistical correction remain separate.
