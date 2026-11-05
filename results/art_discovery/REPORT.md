# ART discovery through independent RT and partner searches

Screened 2076 annotated proteins from 8 genomes. Candidate loci: 9; both routes: 9; RT-first only: 0; partner-first only: 0.

The search panel includes the six local GenBank genomes, matching public annotation for the local SA1 FASTA, and public LPJP1. Several target proteins are also training seeds, explicitly flagged in the locus table. These are recovery results, not an estimate of independent discovery sensitivity.

| Family | Unique seeds | Held-out recovery | Reviewed GNAT hits | Status |
| --- | ---: | ---: | ---: | --- |
| ART_RT | 6 | 2/2 | 0 | limited_related_holdout_check_not_lineage_validation |
| Type_I_partner | 1 | 0/0 | 0 | insufficient_positive_references_single_seed |
| Type_II_partner | 5 | 2/2 | 0 | limited_related_holdout_check_not_lineage_validation |
| Type_III_first_partner | 0 | 0/0 | 0 | not_testable |
| Type_III_second_partner | 0 | 0/0 | 0 | not_testable |

## Discovery rules

Each available family profile is searched against the entire protein database independently. Partner searches are never restricted to proteins near prior RT hits. Sequence E <= 1e-5, independent domain E <= 1e-5 and domain profile coverage >=60% define a provisional retrieval screen. Raw hits to E <= 1e-3 are retained. Calibration and discovery searches use the same fixed HMMER -Z and --domZ search space (the target proteome size). These settings are not calibrated family-membership probabilities or validated bit cutoffs.

For every partner hit the complete gene plus 6 kb on each side is extracted, clipped at contig ends. The first three upstream annotated CDSs are ranked on the partner transcription axis, including other-strand and overlapping CDSs. All are reported with orientation, coordinates, ART-profile evidence and generic RT evidence. A compatible association requires a co-oriented ART or generic RT among those three, within the extracted window. Generic RT evidence alone never sets the RT-first discovery flag. Multiple compatible RTs remain explicit alternative associations.

The type I seed is an accession-level reconstruction of the paper's LPJP1 representative: array-associated RT QXN67953.1 and downstream overlapping gene QXN67952.1. The seed assignment is provisional, not experimentally established. Type I hits require homology to this ART-associated protein, and genomic-context compatibility is recorded separately. A GNAT product label never triggers a type I call.

## Calibration limits

MarsHill and SA1 supply held-out type II partner and RT sequences with no exact training-sequence duplicates; they are related Staphylococcus references, not independent lineages. Three reviewed E. coli GNATs (RimI, RimJ, RimL) are unrelated-family controls; up to 50 annotated phage housekeeping proteins are an additional background screen, not experimentally verified noninteractors. The small control set cannot establish a family-wide false-positive rate.

Type I has one provisional seed and no independent verified positive holdout; its single-sequence HMM is an exploratory homology retrieval model, not a well-validated family profile. Type III first and second genes have separate slots, empty seed FASTAs and explicit unavailable status. No placeholder HMMs or substitute generic domains are used. Supply separate curated seed FASTAs and curation evidence in the configuration to activate them. Their absence means the corresponding discovery routes were not assessed, not that no type III loci exist.

## DNA and annotation limits

The upstream exact-word repeat scan masks annotated CDSs >=300 nt and searches for >=3 approximately regularly spaced copies of a non-low-complexity 10-mer. It is only an inspection hint: there is no shuffled-null significance test or validated array absence call. The windows and full contigs are released for inspection. Only annotated CDSs are searched; de novo ORF rescue, revised gene calling and systematic partial-CDS recovery remain future work.

## Reproduce

```powershell
.venv/Scripts/python.exe scripts/art_discovery.py --config config/art_discovery.json
```

Main output: candidate_loci.tsv. See partner_upstream_genes.tsv, windows.tsv, profile_registry.tsv, calibration/, searches/, seed_manifest.tsv, source_records.tsv and commands.jsonl for evidence.

Method sources: supplied Yoon et al. PDF, family description and partner methods (pp. 8, 36-38); https://doi.org/10.64898/2026.10.05.755809 ; LPJP1 MZ422438.1 at ENA; SA1 MW218148.1 at NCBI; reviewed GNAT controls at UniProt; myRT class profiles https://github.com/mgtools/myRT .
