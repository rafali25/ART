# Completed five-step ART acquisition run

This bounded run adds **21 distinct computational type I candidate pairs**, of
which **19 have complete CDSs**, to the pre-run collection. The merged collection
has **36 type I pairs, 31 complete pairs, and 18 complete pairs passing the frozen
strict ART screen**. It also retains 10 type II pairs as a separate arm.
These counts do not establish biological function or independent samples.

| Step | Executed work | Checkpoint |
|---|---|---|
| 1 | Reviewed four prior additions; retained strict and review tiers; verified 24 frozen baseline hashes | [01](checkpoints/01.md) |
| 2 | Attempted all 45 leads from 44 runs; exactly mapped 24 source contigs; screened them and queried MGnify occurrences/members | [02](checkpoints/02.md) |
| 3 | Downloaded and checksum-verified all 7,230 HPGC records; called/searched 7,151 records after initial DNA deduplication | [03](checkpoints/03.md) |
| 4 | Downloaded and checksum-verified all 10,754 JPGD records; called/searched 4,853 after initial deduplication; retrieved 107 selected UHGV genomes and screened 66 additional DNA sequences | [04](checkpoints/04.md) |
| 5 | Replayed source maps and CDS translations; built separate RT-core/partner trees; deduplicated pairs; evaluated matched blocks | [05](checkpoints/05.md) |

HPGC contributed one distinct complete type I pair. Its type II hit is the old
Machias locus with a different partner gene boundary. JPGD contributed four
distinct complete type I pairs, including one strict ART-screen pair; two other
JPGD pair records were exact duplicates. The archive screen contributed 16
distinct type I pairs. Selected UHGV screening found no ART/partner pair.
See [cohort counts](audit/cohort_counts.tsv) and [aliases](step5/exact_pair_aliases.tsv).

The 36 type I pairs span 22 partner clusters at 70% identity/80% mutual coverage;
the complete-pair subset has 19 such representatives. Cluster counts are distinct
from complete strict-screen samples and from phylogenetically independent loci.

The strict-screen matching graph has **one compatible pair and zero matching
triples**, hence **zero disjoint three-locus blocks** under the proposed criteria
(50–90% identity in both RT core and partner; >=80% mutual coverage; compatible
lengths; distinct source groups). Allowing tree-affinity review candidates gives
six compatible pairs but still zero triples. The proposed eight-block minimum
is unmet. Architecture, phylogenetic independence, calibration, power, and RNA
eligibility remain additional gates. See [matching comparisons](step5/matching_pair_comparisons.tsv)
and [summary](step5/summary.json).

The RT-core tree contains 137 sequences including references; the partner tree
contains 36. FastTree internal values are local SH-like support, not bootstraps.
Tree affinity does not validate low-coverage ART hits or replace missing NART
controls. No new array-significance tests or accepted arrays were produced.

All computation, downloads, and tests ran in Slurm CPU batch jobs under
`cscc-cpu-p`/`cscc-cpu-qos`, following the supplied cluster document. Twelve
mapping tests passed; all 49 raw pair source maps and 98 CDS translations were
replayed successfully. [Final verification](verification.json),
[artifact checks](audit/verified_artifacts.tsv), and [code snapshot](code/) retain
the closing evidence. Original catalogs, profiles, nulls, and significance
results were preserved.

Retrieval limitations remain explicit: 21 archive leads were unresolved within
the 200-MB successful-transfer cap; the MGnify remote Parquet endpoint refused
connections (zero data bytes), so member expansion remains unresolved; 13 of 120
selected UHGV genomes could not be retrieved. These are not absence calls.
A separate IMG/VR download was not performed. The initial DNA registry omitted
the historical type-II source manifest; [source-overlap audit](audit/prior_source_overlap.tsv)
records that overlap, and final locus deduplication removes the Machias alias.
The 7,151/4,853 counts describe records called in the initial registry, not novel
genomes. Slurm accounting (`sacct`) was unavailable because its database
connection was refused; scheduler output and compute logs retain run evidence.

The final SHA-256 manifest covers every durable file in this run except itself,
`verification.json`, and temporary files. Checkpoint revisions are retained in
`checkpoints/history/`; each step's recorded artifact hashes were rechecked.
