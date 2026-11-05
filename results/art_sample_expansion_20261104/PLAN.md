# ART sample expansion: evidence, pilot results, and acquisition plan

Prepared 4 November 2026, Asia/Dubai. Both files in `literature/` were read: the
Yoon study and `ARTNT CLUADE_v1.docx`. The original proposal supplies the research
question; the revised `extra/ART_Research_Proposal_v2.docx` supplies the stronger
sampling and confirmation gate already recorded in HANDOFF.md.

## Decision and sample target

Prioritize additional **complete, source-linked type I RT–partner pairs that fill
new sequence-space and source-group gaps**. Use type II as a separate replication
arm. Keep array-bearing/RNA-eligible samples as a linked subset of the protein
pair collection. A full upstream array is not required to retain an intact
protein pair, but an incomplete upstream window cannot establish array absence.

The baseline catalog contains 13 primary type I pairs, 10 complete pairs, nine
complete representatives, and 12 partner clusters. The user's figure of 19 is
not the number of complete independent type I samples in that catalog. Historical
six- and thirteen-RT testing cohorts also overlap. Use the actual sample tables,
protein completeness, source groups, and deduplicated testing families.

The revised proposal targets 24–36 eligible loci arranged in 8–12 disjoint
three-locus blocks, with both core and partner identities of 50–90%, reciprocal
coverage >=80%, compatible architecture, and partner length ratios 0.8–1.25.
This is a design gate, not a power guarantee. More sequence-divergent clusters
alone will not make every proposed block feasible. Reserve close alternatives
within each subfamily as well as representatives across subfamilies. Complete
pairs need reviewed RT-core/partner trees before independence or block
eligibility can be claimed.

For the next acquisition round, use a planning target of **50–100 traceable
candidate pairs** from which to assess the 24–36-locus gate. This is an operational
target with no promised yield. Record candidate, complete-pair, diversity-group,
eligible-block, and array-bearing counts separately after every batch.

## What the literature supports changing

These are limitations/opportunities, not evidence that every reported Yoon result
is erroneous. PDF page numbers refer to the supplied local study.

| Observation | Evidence in Yoon | Consequence for this project |
|---|---|---|
| Family ascertainment starts close to known phages | pp31–32: twelve RT seeds, ten Staphylococcus and two LPJP1; 350–900-aa representative filtering and >=400-aa expansion | Search diverse RT and partner queries in parallel; retain length/coverage outliers in review categories; build independently assessed subfamily profiles. |
| The published family boundary uses locus architecture | p32: smallest clade containing type I partner loci with >=500 nt upstream noncoding DNA | Keep sequence-family placement independent of partner/array labels; inspect neighboring sequence clades without promoting their members automatically. |
| A representative is not its complete cluster | p32: retrieval of cluster members added 706 proteins | Expand promising public protein clusters and genomic occurrences; a handful of representatives is not a locus census. |
| Initial scan and exceptions have ascertainment limits | p32: top-20 14-mers in 3 kb; one R=3 exception; incompleteness categories | Freeze a consistent final inferential procedure and preserve diagnostics/exceptions separately. |
| Binary predictions do not establish cognate specificity | pp37–38: locus-paired MSAs, limited noncognate controls; mechanism remains tentative | Acquire pairs for matched within-type comparisons, audited unpaired/single-sequence inputs, calibration, and source/lineage-level uncertainty. |

Yoon **already performed** database-wide partner-first HMM searches, annotated
protein plus six-frame rescue, 6-kb recurring-10-mer delimitation, skipped-copy
tolerance, shuffled full-search maxima, and PWM degenerate-copy recovery
(pp32–33,37). These are useful foundations, not missing methods we can claim as
new improvements. The study also explicitly withholds a P value for its dependent
identity pairs (p38), and acknowledges that RT activity, RNA substrates, and
RT–partner interaction remain unknown (p13).

The original ARTNT proposal overstates some possible conclusions. A flat
co-folding matrix would not by itself prove an interface artifact; a weaker
whole-NTD identity correlation would not refute specificity; and RNA-versus-protein
confidence scores are not comparable binding energies. These points strengthen
the need for adequate, matched samples and the revised proposal's controls.

## Searches actually executed in this task

All work is isolated in this dated result directory. Original catalogs, profiles,
null samples, significance calls, and historical result snapshots were preserved.
All additions below are new to this workspace; biological novelty is not claimed.

### Related/public-genome pilot

`panel_manifest.tsv` froze 78 records: 32 sequence records returned for ENA
PRJEB35371, 40 length-prioritized Listeria reference records, named BpSP/AR9
comparisons, and AH12/PB50 near relatives. INPHARED2's attempted TSV export
returned HTML, so the actual metadata tables were parsed and retained. RefSeq
records were recovered from NCBI because ENA does not serve those identifiers.

All 78 records were retrieved; exact DNA deduplication yielded **66 genomes**.
Source CDSs and alternative native Pyrodigal-gv calls yielded **22,472 protein
representations**. Original ART, type I, type II, exploratory two-LPJP1, and
generic RVT_1 profiles were searched. Known LPJP1 controls were recovered.

Two additional complete source-CDS type II pairs pass the frozen retrieval
thresholds and have compatible immediate gene context:

| Source | RT | Type II partner | Protein lengths | Upstream DNA |
|---|---|---|---|---:|
| AH12, OR455461.1 | WNX97802.1 | WNX97801.1 | 517 / 291 aa | 6,000 nt |
| PB50, OR770614.1 | WPH64254.1 | WPH64253.1 | 517 / 290 aa | 6,000 nt |

Both source CDS translations match. Their RTs are approximately 96.3% and 98.1%
identical to the existing Machias RT using global identity over the longer
protein. They broaden the mapped type II collection, without establishing new
independent evolutionary groups. PB50 was already discussed by Yoon; it is not a
newly discovered ART system. Eight-copy seed chains are uncalibrated diagnostics,
not accepted arrays. Protein/CDS/oriented DNA/GenBank/GFF/maps are exported under
`panel_screen/new_loci/`; the clean pair table is
`panel_screen/new_complete_pair_candidates.tsv`.

Environmental LR756488.1 has an RT-like hit with ART-profile coverage 55%, below
the unchanged 60% threshold, and no accepted partner. A follow-up myRT scan also
finds a strong retron-profile match. It remains class-unresolved comparative
evidence, not a new primary ART sample or a verified NART control. Raw matches are
in `panel_screen/environmental_RT_myRT.domtbl`; scores across different profiles
are not calibrated class probabilities.

### Deferred/archive leads

The frozen lead manifest contains 280 previously unacquired MGnify accessions and
48 previously unacquired archive source-run leads. Bounded source retrieval used
4,919,302 successful response bytes and recovered three Logan contigs plus one
MGnify contig. Another MGnify contig returned 404.

The new Logan source **SRR7532897_75** contains a complete primary type I pair:
RT 607 aa, ART E=8.2e-51/coverage 87.6%; partner 851 aa, passing the original type I
profile. Both CDSs retranslate exactly and are mutual first neighbors. Historical
contig number 76 was reconciled with current number 75 by exact peptide-to-DNA
mapping.

Only 66 nt upstream are available, so array assessment is incomplete. Combined
MMseqs2 clustering of the baseline 13 pairs plus this pair gives **14 pairs but
still 12 partner clusters**. It adds a distinct sequence pair, not a new diverse
representative. The other retrieved partner contigs did not yield an RT under
the tested gene/six-frame screens. **45 archive source-run leads remain
unacquired.** See `deferred_acquisition/README.md` and
`deferred_acquisition/processed/primary_pair_summary.tsv`.

### Diverse-query MGnify search

Three complete, source-retranslated public partner representatives from different
existing clusters were submitted to live EBI phmmer, instead of reusing only
LPJP1. MGnify30-C2 release 2026_07 was verified before submission. Completed
searches returned 30 rows / 11 unique proteins; only two peptide accessions were
absent from earlier saved lists. Both new source contigs were recovered and
peptide mappings were verified, using 254,236 source bytes.

MGYC003966598318 contains an intact generic-RT/type-I-partner context pair, but
the ART RT profile covers only 37.2%. Its RT/partner are 98.03%/97.44% identical
to an existing pair from the same study, ERP150573. It is a partner-first review
lead, not an additional passing ART-family pair or a new diversity group. Four
repeat copies are exploratory without null significance. The other new source
is too short to assess upstream RT context. See `diverse_partner_search/README.md`.

### A genuinely additional environmental dataset

Downloaded the published **139-genome curated Jug subset** of the HPGC project
(14,678,419 bytes) and verified its published MD5. This yielded 138 exact DNA
sequences and 69,120 predicted proteins. Five profile searches reported no domains
at E<=0.001. DNA translations, coordinate mappings, and known-seed positive
controls at the same search-space size passed. This is a negative result for this
subset and tested representations; it does not establish absence across HPGC or
all huge phages. See `jug_screen/README.md`.

**Executed yield:** one additional passing complete type I candidate, two
additional complete type II candidates, one additional type-I-context review
lead, and **zero newly accepted arrays or demonstrated new eligible blocks**.
The baseline catalog has deliberately not been silently merged with these tiers.

## Ranked next datasets and viral groups

| Priority | Search substrate | Concrete next action | Reason and qualification |
|---|---|---|---|
| 1 | Remaining 45 archive source-run leads and all occurrences/members of promising public protein clusters | Rank source context/completeness first; expand 100–300 DNA-mapped leads across source studies and sequence strata; retrieve bounded contigs/windows | The first archive follow-up yielded a complete pair at low transfer cost. Avoid spending the whole budget on the closest homologs. |
| 2 | Full HPGC: 7,230 nonredundant huge-phage genomes, 560,950,249 compressed bytes | Download once, verify published MD5, screen mapped genes in diverse clades/habitats, excluding already tested exact genomes | Purpose-built huge-phage substrate; the Jug subset's zero yield favors breadth across other clades. It is not 7,230 ART systems or independent samples. |
| 3 | JPGD: 10,754 jumbo-phage genomes, 871,869,424 compressed bytes | Deduplicate against HPGC/baseline before gene calling; retain alternative-code and completeness metadata | Adds multiple habitats and lineages. Published metadata includes genetic code; overlap and assembly quality need review. |
| 4 | INPHARED2 / current NCBI complete isolates | Search Wallmarkvirinae, including Madawaskavirus/Machiasvirus/Lentusvirus, plus accession-linked LPJP1-like relatives and broader Firmicutes phages | Useful for complete contexts/type II replication. Current taxonomy rather than old family labels matters; close Staph genomes alone do not solve type I diversity. |
| 5 | UHGV and high-confidence IMG/VR | Search their representative proteins locally, expand matching members, recover linked DNA/GFF and source studies | Broader environmental/gut coverage; avoid counting overlapping catalogs repeatedly. Check actual JGI download access. |
| 6 | Targeted Logan assemblies, including oral/gut, soil, freshwater and wastewater studies | Retrieve assemblies/contigs supported by specific homologs or related viral core markers; use translated local search near leads | Supports long context and annotation rescue. Prioritize plausible sources rather than a blind planetary-scale download. |

Current NCBI places MarsHill in Wallmarkvirinae/Madawaskavirus. LPJP1 should not
be assigned to Herelleviridae simply because it infects Listeria. Herelleviridae
and AR9/BpSP-like Bacillus phages are useful broader comparison strata; host,
capsid homology, genome size, or taxonomic proximity does not establish ART
homology. Cultured AH12/PB50 retrieval is evidence for testing the exact named
relatives rather than relying solely on a possibly outdated taxonomy filter.

The new MGYC003966598318 and its existing relative come from Papua New Guinean
oral metagenome assemblies (ERP150573/PRJEB65450 derived from PRJEB54966).
This gives a concrete oral-metagenome expansion stratum. Use other studies and
related viral genomes as well, because two assemblies from that study are not
independent cohorts.

## Smarter search and retrieval workflow

1. **Freeze public inputs and divide development from evaluation by source/lineage.**
   Keep full RT, RT-core, and partner sequences linked to DNA. Preserve original
   profiles beside experimental expansions. Never add every new hit to training
   and describe re-recovery as independent validation.

2. **Search complementary sequence representations.**
   Use DIAMOND/MMseqs2 for fast broad protein retrieval, followed by HMMER using
   separate ART RT/core and partner subfamily profiles. A variable NTD can reduce
   whole-protein coverage; retain strong core/partial matches for review without
   lowering the frozen acceptance gate. Use translated DNA near partner hits and
   viral gene calls to resolve missed genes. Assess genetic code explicitly;
   introns, frameshifts, and uncertain starts remain review issues rather than
   fabricated complete genes. Structure-profile tools such as HHsearch/Foldseek
   can prioritize weak candidates after sequence/context checks; a generic GNAT
   fold is not a type I assignment.

3. **Expand members and occurrences rather than repeatedly searching the same representatives.**
   MGnify30-C2 searches only representatives of non-singleton 30%-identity
   clusters. Three diverse queries produced only two incremental accessions here.
   MGnify's newly documented 2026_07 Parquet files expose cluster members and
   protein occurrences with contig length, coordinates, strand, biome, and study.
   Use filtered DuckDB/Parquet queries or equivalent indexed retrieval to obtain
   source contexts for promising clusters, then download selected DNA. This is a
   metadata/member-expansion step, not a replacement for homology search.
   Benchmark a small predicate-restricted query's byte/time cost first: remote
   queries, joins, or sorting can still scan large columns. Singletons excluded
   by the web search need a separate local/public DNA search strategy.

4. **Rank for usable context and new sequence space.**
   Prefer complete partners with an upstream RT and sufficient upstream DNA;
   retain intact protein-only pairs separately for the protein arm. Select a
   mixture of close alternatives needed for 50–90% matched blocks and distinct
   subfamilies. Penalize repeated BioSamples, identical DNA, >95% protein-pair
   redundancy, and overrepresented studies. Do not assume every new accession,
   partner cluster, or contig is an independent sample.

5. **Verify DNA, pairing, and family evidence before expensive array inference.**
   Retranslate both CDSs, record code/starts/partial flags, reconcile old/current
   sequence identifiers, require compatible orientation and actual gene context,
   and retain alternative partners. A non-ART RT plus a type I-like partner is a
   review tier. Keep competing retron/DGR/CRISPR profiles and sequence-only trees
   separate from architectural labels.

6. **Assess arrays only in the declared, mapped candidate family.**
   Retain independent initial/extended scans, full-search local/dinucleotide
   nulls, quality flags, and source-coordinate maps. Resolve the existing three
   copy-precision benchmark failures on development data and test fresh held-out
   motifs before promoting a changed method. Freeze final testing-family
   membership and refresh compatible nulls/corrections when code or membership
   changes. Do not borrow old p-values for new loci.

   Monte Carlo resolution matters: with N shuffles, minimum empirical p is
   1/(N+1). Bonferroni across m loci requires at least
   `ceil(m/0.05)-1` shuffles merely to make significance mathematically possible.
   At N=1,000, m>=51 cannot pass; at N=2,000, m>=101 cannot pass. These are
   resolution limits, not recommended precision targets. Larger inference needs
   a predeclared larger null budget or a validated alternative procedure. Do not
   split one family into convenient small batches to evade multiplicity.

7. **Stop acquisition against eligibility, cost, and incremental diversity.**
   For each frozen batch export: genomes/proteins searched; mapped loci; complete
   primary pairs; exact/near duplicates; source groups; new partner/core strata;
   eligible disjoint blocks; supported arrays; bytes/CPU; deferred/excluded
   records. Redirect later batches toward underfilled clades/blocks when yield
   is dominated by near duplicates. Proceed to confirmatory structural work only
   after the revised design, calibration, simulation, and input-control gates.

## Immediate execution sequence

1. Review the three passing new pairs and the separate MGYC review lead, keeping
   their tiers and incomplete/uncalibrated array status explicit.
2. Run a small public cluster-member/occurrence expansion and rank the remaining
   45 archive runs by context, novelty, and study balance.
3. Download/checksum the full 561-MB HPGC archive and screen the nonduplicated
   genomes with the native, provenance-recording workflow. Preserve the observed
   Jug zero-yield result; use fresh source strata for further sampling.
4. Add the nonoverlapping JPGD collection, then selected UHGV/IMG/VR strata if the
   complete-pair and block yield is still inadequate. Fix each manifest and
   search settings before measuring that batch's yield.
5. Build reviewed separate RT-core/partner trees and the v2 block-feasibility
   table. Array inference and RNA eligibility are additional evidence layers.

No full HPGC/JPGD archive or blanket large collection was downloaded in this
task. Exact file inventories, sizes, URLs, licenses and MD5 values are saved in
`dataset_inventory/files.tsv`, so the next acquisition is concrete and reviewable.

## Sources and reproducible artifacts

- Local Yoon PDF and original proposal in `../../literature/`; revised proposal
  in `../../extra/ART_Research_Proposal_v2.docx`.
- [INPHARED2](https://phagebank.uk/inphared2/), current displayed update 3 October 2026.
- [NCBI MarsHill lineage](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2796499).
- [ICTV 2026 proposal summaries](https://ictv.global/system/files/proposals/pending/2025/All%20Proposals/Taxonomy-Proposal-Summaries-for-Ratification-Vote%2C-2026.pdf).
- [Al-Shayeb huge-phage study](https://www.nature.com/articles/s41586-020-2007-4)
  and ENA PRJEB35371; only 32 sequence records were returned by this task's
  project query, not the entire published 351-phage inventory.
- [HPGC/Jug study](https://doi.org/10.1126/sciadv.aeb6265),
  [HPGC data](https://figshare.com/articles/dataset/HPGC_genomes/30506864),
  [curated Jug data](https://figshare.com/articles/dataset/Jug_phages_curated_genomes/29414588).
- [JPGD study](https://doi.org/10.1186/s42523-026-00534-z),
  [JPGD sequences](https://zenodo.org/records/14974448).
- [UHGV v1.0](https://zenodo.org/records/17402089),
  [IMG/VR release/file documentation](https://genome.jgi.doe.gov/portal/IMG_VR/IMG_VR.home.html).
- [MGnify search scope](https://docs.mgnify.org/src/docs/mgnify-proteins-sequence-search.html),
  [Parquet/remote-query schemas](https://docs.mgnify.org/src/docs/mgnify-proteins-parquet-queries.html),
  [ERP150573 study metadata](https://www.ebi.ac.uk/metagenomics/api/v1/studies/MGYS00006631?format=api).

The native related-panel driver is `../../scripts/search_art_expansion_panel.py`
(`discover`, `acquire`, `screen`, `export`). Each other branch retains its own
scripts, source requests, commands, hashes, verification and README. Rerunning
`acquire` may fetch public data; rerunning `screen`/`export` refreshes only this
dated panel's artifacts. Existing scientific pipelines and catalog were not
rerun or reclassified.
