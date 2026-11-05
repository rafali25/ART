# Type I ART pilot

Start with [REPORT.md](REPORT.md), [primary_typeI_pairs.tsv](primary_typeI_pairs.tsv)
and [candidate_loci.tsv](candidate_loci.tsv). Results are provisional.

The pilot adds four mapped public RT–partner loci to the two LPJP1 controls.
One new locus has a statistically supported candidate array. Two new loci have
insufficient upstream DNA to assess array absence. No new type I locus was found
by the preceding local FASTA rescue audit.

## Outputs

- `primary_typeI_pairs.tsv`: six primary pairs; excludes nearby alternative homolog associations.
- `typeI_pairs.tsv`: all primary and alternative associations, explicitly labelled.
- `candidate_RT_proteins.faa`, `primary_typeI_partners.faa`, `NTDs.faa`: protein sequences.
- `primary_pair_CDS.fna`, `upstream_sequences.fna`: extracted DNA.
- `loci/`: oriented GenBank/FASTA/GFF3 and original-coordinate maps.
- `partner_window_coordinate_maps.tsv`: source coordinates for partner ±6 kb windows.
- `repeats.tsv`, `spacers.tsv`, `repeats.fna`, `spacers.fna`, `alignments/`:
  inferred repeat boundaries, alignments, sequences and spacer gaps.
- `arrays_exploratory.fna`: all inferred tracts; use `background_tests.tsv` for status.
- `array_details/`: full seed, chain, PWM and quality evidence, plus background status.
- `sequence_family_placement.tsv`, `domain_phylogeny/`: sequence-only RT evidence.
- `observed_locus_architecture.tsv`: architecture evidence kept separate.
- `competing_RT_family_matches.tsv`: myRT HMM matches; scores from different
  profiles are not interchangeable classification cutoffs.
- `archive_leads/`: external protein leads, unrelated GNAT controls, cached
  Logan prefixes, ETags and historical-to-current mapping evidence.
- `source_records.tsv`, `input_sha256.json`, `output_sha256.json`, `verification.json`:
  provenance and sequence/coordinate verification.

All internal coordinates are **zero-based, half-open**. GFF3 and GenBank writers
perform their explicit external coordinate conversion. Logan IDs are run/contig
identifiers backed by object ETags and checksums, not fabricated accession versions.
NTDs are operationally defined as sequence before the first significant PF00078
envelope; they are not experimentally established structural boundaries.

## Re-run the bounded pilot

From `D:\MBZUAI\CRISP`, using the existing project environment and reference files:

```powershell
.venv/Scripts/python.exe scripts/run_art_type1_pilot.py
```

This reuses checksum-checked downloads and array nulls. The fixed panel is four
archived leads, with historical and current assembly prefixes totalling at most
60 MB compressed. The historical archive and installed tools are existing project
dependencies. The pipeline verifies the target peptide in current DNA; it does
not trust old contig numbering or promote a protein-only hit to a locus.

Individual analysis stages:

```powershell
.venv/Scripts/python.exe scripts/art_type1_pipeline.py process
.venv/Scripts/python.exe scripts/art_type1_pipeline.py nulls
.venv/Scripts/python.exe scripts/art_type1_tree.py
.venv/Scripts/python.exe scripts/art_type1_pipeline.py diversity
.venv/Scripts/python.exe scripts/finalize_art_type1_public.py
.venv/Scripts/python.exe scripts/verify_art_type1_public.py
```

Run these stages sequentially because they update shared result tables. Finalization
synchronizes the tree and array results and exports CDS/repeat/spacer tables.
The null calculation tests all candidate loci, including controls and short loci,
as one declared cohort. Potential positives require 2,000 shuffles per background;
the 200-draw screen cannot independently yield a positive decision.

## Resume separately queued public homology searches

Four HMMER/MGnify jobs and one four-query NCBI BLAST batch are recorded in
`jobs.json`. The four MGnify jobs have completed; their bounded follow-up is in
`results/art_type1_mgnify/`. NCBI BLAST remains pending. Pending searches are
**incomplete**, not negative.

```powershell
.venv/Scripts/python.exe scripts/art_type1_public_search.py poll
.venv/Scripts/python.exe scripts/art_type1_public_search.py collect
```

To acquire a separate follow-up from completed hits, use its config:

```powershell
.venv/Scripts/python.exe scripts/art_type1_public_search.py acquire --config config/art_type1_mgnify.json
```

Then run the six sequential analysis stages above with the same `--config`.
Acquisition preserves the
existing archived-lead sources, retrieves bounded nucleotide records, checks
translation and records unresolved protein-to-contig mappings. Its own limits are
limits in that config (the initial MGnify follow-up uses ten contigs) and 20 MB
of newly acquired source records per invocation. The
archive-prefix budget is separate. Database search caps and unavailable contigs
limit recall; no prevalence estimate is supported.

Do not resubmit jobs merely because polling remains pending. Remote response
normalization and MGnify protein-to-contig retrieval have now been exercised on
completed responses. NCBI result acquisition remains untested while its job is queued.

## Interpretation limits

The exploratory type I profile has two seeds from one genome and no verified
independent positive holdout group. Five unrelated GNAT negatives are a small
specificity check. Published/local validation still has repeat-copy precision
failures and lacks independently confirmed full boundaries and NART references.
Array significance does not establish ART function. The two type III families
remain separate unavailable-reference families; no type III profile is inferred
from these type I sequences.
