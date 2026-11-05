# Type I ART acquisition pilot

Implemented and run on four public environmental contigs plus LPJP1. Four new provisional primary RT–type I partner pairs were recovered, giving six primary pairs in this pilot. These are new to this workspace, not newly discovered biological families. The four public targets were supplied as leads by an external archive, not a representative prevalence sample.

| Current contig | RT | Partner | Available upstream | Array assessment | Spacer lengths, nt |
|---|---|---|---:|---|---|
| SRR8925777_19723 | SRR8925777_19723_pgv9 | SRR8925777_19723_pgv10 | 6000 | not supported | [] |
| ERR10896470_6836 | ERR10896470_6836_pgv1 | ERR10896470_6836_pgv2 | 232 | incomplete upstream | [] |
| ERR8055547_17635 | ERR8055547_17635_pgv12 | ERR8055547_17635_pgv13 | 6000 | supported statistically; provisional boundaries | [284, 326, 630, 289] |
| MZ422438.1 | QXN67953.1 | QXN67952.1 | 6000 | supported statistically; provisional boundaries | [148, 264, 137, 199] |
| SRR26805139_340 | SRR26805139_340_pgv1 | SRR26805139_340_pgv2 | 233 | incomplete upstream | [] |
| MZ422438.1 | QXN67955.1 | QXN67954.1 | 6000 | not supported | [] |

The six primary partners form **5 clusters** at 70% identity and 80% mutual coverage. Different contigs or protein clusters do not establish independent biological samples. Two additional nearby partner homolog associations are retained separately; they are not counted as primary pairs.

## What was implemented

- Reproducible MAFFT/HMMER profile construction and independent ART-first and partner-first retrieval. The additional two-seed type I model uses both LPJP1 paralogs; the original profiles remain unchanged.
- Original annotations, Pyrodigal-gv predictions and six-frame ORF rescue; exact protein-to-DNA mapping. First-three-upstream-gene checks include orientation. An intervening RT prevents assigning a partner across that RT.
- Partner ±6 kb windows, RT-oriented loci, source coordinate maps, GenBank/GFF annotations, RT/partner/NTD proteins, CDS DNA, upstream DNA, repeat alignments, repeats and spacers.
- Independent 3 kb initial and 6 kb extended repeat searches. Spacers use inferred full repeat ends; published length ranges are descriptive, not rejection filters.
- Full seed rediscovery and best-chain search on 200 local-block and 200 dinucleotide shuffles per locus. Potentially supported candidates advance to fixed 2,000-draw tests per background. The 200-draw screen can stop only for non-support; positive status requires the fixed 2,000-draw test. The larger background p-value is Bonferroni-adjusted across all six loci. This controls the declared pilot comparison, not the earlier external lead-selection process or a database-wide discovery rate. Each locus also has 100 nucleotide shuffles for the initial scan.
- A separate 97-tip RT-domain tree with labelled RT classes, plus competing myRT-family HMM matches. The tree uses sequence only: PF00078 ±40 aa, MAFFT L-INS-i, 50% occupancy trimming, FastTree LG+Gamma. All six queries have descriptive ART-reference affinity; the panel cannot resolve ART versus NART.

## New array

ERR8055547_17635 has five inferred 16-nt repeats (`CCCTAGGATCTAGGGT`) spanning 1,609 nt. Source repeat intervals (zero-based, half-open) are [17737,17753), [18037,18053), [18379,18395), [19025,19041), [19330,19346). Spacers are 284, 326, 630 and 289 nt; the large gap is consistent with a missed/absent repeat. There are no detected CDS overlaps or genes between this tract and the RT. Neither background produced a score at least as high in 2,000 shuffles, so each empirical p is 1/2001 and the six-locus corrected p is approximately 0.003. This supports the repetitive tract, not its biological function.

## Acquisition and provenance

The existing local FASTA rescue audit found no extra type I loci. Public leads came from [the AI3 supporting archive](https://doi.org/10.5281/zenodo.23086637), whose 89-sequence homolog list is treated as unverified discovery input. Thirty leads pass the original type I profile and 35 pass the two-seed profile. Five archived unrelated GNAT controls had no reported profile matches; this small set does not establish broad specificity.

All four target peptides were exactly mapped in Logan v1.0 and independently located in the current contigs. Contig numbers changed; numeric identifiers alone were rejected as mapping evidence. Three current DNA sequences are identical to their historical contigs; SRR26805139_340 is longer (16,150 nt versus 5,048 nt), with the partner still exactly encoded. Current sequences, object ETags, compressed prefixes, download metadata and SHA-256 checksums are retained. Logan run/contig IDs are not represented as INSDC accession versions. The bounded historical/current retrieval used 59,602,523 compressed bytes in total.

The four EBI HMMER jobs have completed; NCBI BLAST remains queued at the latest poll. The MGnify follow-up is stored separately in `results/art_type1_mgnify/`; its larger testing cohort does not replace this six-locus pilot. Resume with the commands in README.md.

## Limits and next selection

No expanded profile has been trained on the new candidates. The exploratory type I model still has one training-genome group and no verified independent positive holdout group. Known synthetic repeat-copy precision failures and incomplete full-boundary validation remain unresolved. NART references and both verified type III partner families remain unavailable. Two new loci have only 232/233 nt upstream and cannot establish array absence. The other array-negative public locus has 6 kb upstream but no supported tract under this detector. Keep all these labels provisional; prioritize the ERR8055547 locus for independent boundary review and use primary partner clusters to select distinct protein pairs for the proposal.

The proposal's DIAMOND/MMseqs2/FAMSA/FastTree replication workflow remains separate in `results/art_family_reconstruction/`; this bounded acquisition pilot does not claim to reconstruct the entire published ART family.
