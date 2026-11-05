# Finding and extracting additional type I ART loci

Prepared 20 October 2026 from the supplied Yoon paper, the ARTNT proposal, and the current local results.

**Recommendation: use a type I partner-first search of public environmental proteins, recover their source contigs, and run an independent ART RT search in parallel.** Reusing the local FASTAs is useful for rescuing missed gene annotations, but the current files contain the same seven Staphylococcus genomes already examined. A new annotation-rescue audit found no additional type I candidate at the available profile thresholds.

This document specifies the next acquisition workflow. The local audit and the existing LPJP1 sequence package are completed. The public acquisition search described below has not been run in this task.

## What the paper implies for sampling

The supplied paper reports type I partners at 59 of the 95 family representatives, principally in Listeria and environmental lineages; type II is associated with the Staphylococcus clade (PDF pp. 8-9). Its broader partner census retrieves 1,011 type I partner loci, of which Supplementary Figure 7 marks 712 as having an adjacent ART RT (p. 21). These are different quantities from independent, intact, array-bearing systems. The paper delimits arrays adjacent to RTs at 28 of the 95 family loci overall, without establishing that all type I pairs have an observable array.

The census used 46 type I partner sequences to construct its profile, compared with our current one-seed profile. The reported scale therefore does not establish the sensitivity of our current search or guarantee a comparable number of accessible samples. Original accession mappings, if obtained, would be valuable starting material; internal labels such as ART_57 or L0050 cannot substitute for a retrievable contig.

Your proposal's Aim 1.2 calls for an independent partner-first census, Aim 1.3 clusters partners at 70% identity, and Aim 2 requires cognate RT-partner pairs. Aim 3 additionally requires the associated array RNA sequence. Maintain two linked cohorts:

| Cohort | Required evidence | Proposed use |
| --- | --- | --- |
| Type I RT-partner candidates | Traceable same-contig proteins; type I family homology; ART-like RT sequence placement; compatible orientation/context; complete proteins | Acquisition pool for the protein-pair analysis, with array status retained |
| Type I array-bearing candidates | The above plus a reviewed upstream array and suitable statistical evidence | Array comparisons and candidate RNA-containing complexes |

An array-unresolved protein pair must remain explicitly labelled. It is not interchangeable with an array-bearing locus. The proposal's statement that type I is already well powered is a feasibility hypothesis: the usable independent sample count still has to be measured.

## 1. Bootstrap with two traceable LPJP1 pairs

Use both LPJP1 RTs and both partners as separate homology queries:

| RT | Type I partner | RT / partner length | Current array evidence |
| --- | --- | --- | --- |
| QXN67953.1 | QXN67952.1 | 454 / 481 aa | Existing five-copy computational candidate |
| QXN67955.1 | QXN67954.1 | 509 / 661 aa | Not supported by the current corrected array assessment |

Both belong to the same source genome, MZ422438.1. They are paralogous examples, not two independent host/phage lineages. Their identities, source CDS translations, reverse-strand orientation and coordinates have been checked against the deposited record.

The [current_typeI](current_typeI/) directory contains both full RTs, partners, NTDs, CDSs, upstream regions, oriented GenBank/GFF loci, coordinate maps and the existing array/repeat/spacer outputs. The NTD boundary is the sequence before the first significant PF00078 envelope; it is an annotation proxy, not a demonstrated structural boundary. Repeat-spacer-repeat DNA units are candidate inputs for RNA modelling, not experimentally measured transcript boundaries.

Use both partners in initial BLASTP/phmmer searches. A two-sequence alignment can provide an additional exploratory profile, but it must not be called a validated type I family model. Keep the current frozen production profile as a comparison. A generic GNAT domain, protein length near 600 aa, or predicted helical/acetyltransferase structure alone is insufficient.

## 2. Acquire new source DNA through targeted homology searches

The best immediate expansion is homology-directed retrieval from environmental sequence collections, alongside public annotated phages. Avoid selecting more genomes solely because they resemble MarsHill.

| Priority | Resource and operation | Required mapping |
| --- | --- | --- |
| First | MGnify protein sequence search with both type I partners and both LPJP1 RTs; union the results | MGYP to study, assembly, contig, coordinates and strand; recover nucleotide sequence |
| Parallel | Public GenBank proteins/annotated phages: BLASTP or phmmer with the same queries; TBLASTN on selected nucleotide assemblies for annotation rescue | Protein accession/version to the actual CDS on a nucleotide accession/version |
| Next | A fixed release/subset of IMG/VR, searching predicted proteins and recovering viral contigs | Protein ID to viral contig ID and genomic coordinates; validate by translation |
| Later | Selected Logan assemblies whose sample/contig identifiers come from promising hits | SRA run, Logan release/object, contig identifier, sequence checksum and gene calls |

MGnify's protein search and metadata API are different services: the metadata `/protein/search` endpoint filters by biome or Pfam and is not a homology search. Its API documentation also says detailed records cover cluster representatives. For members without a direct record, use the published cluster/member and contig-retrieval procedures; retain failures as `protein_only_unmapped`. Never count them as recovered loci. The general MGnify API is currently v2; the separate Proteins API is v1.

IMG/VR distributes both protein and nucleotide FASTAs plus viral metadata. Check that coordinates are present in the selected release; if they are not, obtain source annotations or call genes on the retrieved contigs and map the homolog by its actual translation. A protein header alone is not a coordinate map.

For a first targeted public batch, freeze the queries, database releases, retrieval date and all result identifiers. A proposed manageable budget is at most 200 reported homologs per query and 30 distinct, DNA-mapped contigs selected across sequence clusters, with a 20-MB source-download cap. These are operational budgets, not biological filters. Record truncated hit lists and deferred records. This is a directed discovery panel, not an estimate of prevalence or sensitivity.

Prioritize complete loci and sequence diversity. Do not require a known Listeria host, an annotation containing "ART", a contig longer than 200 kb, a repeat-array annotation, or a particular spacer length. Environmental host labels and fragmented assembly lengths can be misleading.

## 3. Search proteins by both routes and rescue missed genes

For every selected assembly:

1. Retain the original FASTA and annotations, source identifiers, versions and checksums.
2. Search the annotated proteins with the ART RT and type I partner profiles independently. Also compare generic RT and available retron, DGR and CRISPR-associated RT references.
3. For unannotated contigs, make Pyrodigal-gv candidate gene calls, recording model and genetic code. For annotated contigs, keep alternative calls in a separate layer.
4. Around every partner hit, retrieve the partner gene plus approximately 6 kb on each genomic side. Inspect the first three upstream genes in partner transcriptional direction, including opposing genes and overlaps.
5. Search six-frame translated segments of at least 40 aa in the candidate windows, as in the paper's census, to rescue missing/fragmented RT annotations. Use the supported translation table and a separately labelled alternative-code sensitivity analysis when indicated. A stop-to-stop fragment is not a complete gene.
6. Merge overlapping representations of the same genomic gene and record discovery route and annotation source; do not count a deposited CDS, a new prediction and a six-frame fragment as three loci.

For provisional local retrieval, record sequence/domain E-values, score and profile coverage; the current audit uses reporting E <= 1e-3 and provisional sequence/domain E <= 1e-5 with profile coverage >= 0.60. Keep weaker/partial hits in a review tier. These settings retrieve candidates, not calibrated family probabilities. Do not transfer the paper's ART HMM score of 207 bits to a new profile.

A compatible principal pair has a co-oriented RT within the first three upstream genes, with direct adjacency strongest. Record the RT-partner gap and intervening genes. The paper used a <=60-nt RT-partner spacing to assemble a particular paired MSA; it is not a universal discovery filter. RT length 350-900 aa belongs in the separate replication workflow, not a hard acquisition exclusion.

## 4. Improve type I coverage without contaminating validation

Review candidate type I homology across substantial portions of the protein, its domain architecture, and competing unrelated GNAT-containing proteins. In parallel, evaluate RT placement using sequence-only profiles/tree references; array and partner architecture must not define the sequence-placement label.

Add clearly mapped homologs to a separate provisional expansion set. Before rebuilding the type I HMM, reserve entire relatedness groups for evaluation. Align curated positives with MAFFT L-INS-i, inspect and trim poorly supported terminal stretches conservatively, then use HMMER hmmbuild. If distinct type I subfamilies align poorly across full length, use subfamily profiles with provenance and calibration for each, rather than forcing all sequences into one model.

Calibrate retrieval on held-out lineages and unrelated GNATs, other phage proteins and the other partner families. Keep the original seed and expanded-model results side by side. A candidate recovered by a newly trained model is not independent validation of that model. Two LPJP1 partners cannot establish broad sensitivity; NART specificity also remains unresolved until traceable NART references are available.

The paper's three type I-associated NTD structural groups are diversity targets only after those structures/groups can be established. Do not infer those groups from partner length or from an arbitrary count of sequence clusters.

## 5. Extract and assess arrays independently within each locus

For an RT on source interval [a,b), extract at least 6 kb upstream where available and enough downstream DNA to contain the full candidate partner. A simple 6-kb-flank interval is:

```text
source interval = [max(0, a-6000), min(contig_length, b+6000))
plus-strand RT: keep the source orientation
minus-strand RT: reverse-complement the whole extracted interval
```

For a minus-strand RT, its upstream sequence is at larger source coordinates, beyond b, before reverse-complementing. If a mapped partner extends outside the initial window, extend to contain the full partner. Save piecewise coordinate maps. Only wrap a boundary when circularity is supported. Preserve original and alternative annotations and report upstream shortfall, uncertain starts, ambiguity and intervening genes.

Run both scanner modes:

- The published 3-kb/top-20 14-mer mode provides a replication diagnostic.
- The 6-kb/all-recurring-10-mer mode is mandatory even if the first mode fails.

Discover seeds anew within each locus; chain approximate occurrences by dynamic programming at approximately regular spacing with a missing-copy allowance; align repeat-centred regions; estimate boundaries and a PWM; rescan for degenerate copies. Do not use the MarsHill repeat as a universal motif or fix type I repeat length to the current LPJP1 example.

Spacers are `next_repeat_start - current_repeat_end`. The current LPJP1 candidate has spacer lengths 148, 264, 137 and 199 nt; the 264-nt interval demonstrates why 120-220 nt cannot be an absolute filter. Review it for a missing copy or uncertain boundaries, while retaining the observation.

Rerun complete seed discovery and best-chain selection on local 50-nt-block shuffles and dinucleotide-preserving shuffles. Keep paper-replication and inferential rules separate. For a new fixed cohort of L loci using Bonferroni alpha=0.05, choose the null size before testing: at least `ceil(L/0.05)-1` shuffles are needed merely to make the minimum possible empirical p-value small enough. For 30 loci, 599 is the mathematical minimum; 1,999 per null is a proposed starting budget. Larger cohorts require more resolution or a pre-specified alternative multiple-testing design. Do not add newly discovered loci to the old nine-locus corrected probabilities.

The existing detector has three synthetic cases with an extra spurious copy, despite recovering the true repeats. Address that copy-selection problem on development cases and evaluate on fresh held-out motifs before accepting a larger array cohort. A leave-one-copy-out PWM score is a useful proposed safeguard against a copy supporting itself; it must be validated and any altered test statistic requires new null simulations.

Flag coding overlap, stronger short-period repeats, low complexity, intervening genes, contig ends and ambiguous bases. Keep two-copy candidates exploratory. Array-unresolved candidates remain in the protein-pair acquisition pool with that status.

## 6. Select samples for the proposal and export them

Cluster partners at 70% identity with an explicit coverage rule (proposed: 80% of both proteins, MMseqs2 `--cov-mode 0`). Report sensitivity to identity/coverage settings, because partial proteins can distort clustering. Also report RT-core groups, assembly/sample provenance and duplicated genomes. A 70% protein cluster is an operational diversity group, not proof of statistical independence.

Choose complete pairs spanning partner and RT/NTD diversity, prefer a well-supported locus on the longest suitable contig, and avoid counting repeated assemblies or the two LPJP1 paralogs as independent genomes. Ten or more diverse intact type I pairs would be a useful initial feasibility target for a pilot matrix, followed by a power calculation using observed variance; it is not a guaranteed or sufficient sample size.

Export per locus:

```text
locus_id, source_accession_and_version, assembly_and_sample_id, contig_id
RT_coordinates_and_strand, RT_protein_sequence, NTD_boundary_and_sequence
ART_profile_score_coverage_and_competing_matches, RT_completeness
partner_family, partner_coordinates, partner_protein_sequence, partner_coverage
discovery_route, gene_call_method, translation_table
oriented_locus.fna, original_and_alternative_annotations, coordinate_map.json
upstream_sequence_available, array_status, repeats.fna, repeat_alignment.fna
repeat_source_coordinates, consensus, spacers.fna, spacer_lengths
array_score, empirical_p_values, cohort_adjusted_p, flags
protein_sequence_assignment, observed_architecture
partner_cluster_70, RT_group, sample_group, final_evidence_category
```

Stop counting at each evidence level: homologous proteins, DNA-mapped loci, intact compatible pairs, independent diversity groups, and supported array-bearing loci. This makes clear which samples support Aim 2 and which can also support Aim 3.

## Completed local check

See [AUDIT.md](AUDIT.md) and [summary.json](summary.json). Rerun:

```powershell
.venv/Scripts/python.exe scripts/audit_art_type1_fasta.py
.venv/Scripts/python.exe -m unittest discover -s scripts -p test_art_type1_fasta.py
.venv/Scripts/python.exe scripts/verify_art_type1_acquisition.py
```

## Sources

- Supplied Yoon PDF: `literature/Yoon et al. 2026 Array Associated Reverse Transcriptases (2).pdf`, Results pp. 6-9, Supplementary Figure 7 p. 21, array/family methods pp. 31-33, partner census p. 37.
- Supplied proposal: `literature/ARTNT CLUADE_v1.docx`, Aim 1 Tasks 1.1-1.4 and Aims 2-3.
- [LPJP1 source nucleotide record](https://www.ncbi.nlm.nih.gov/nuccore/MZ422438.1).
- [MGnify sequence-search portal](https://www.ebi.ac.uk/metagenomics/proteins/), [Proteins API and coordinate metadata](https://docs.mgnify.org/src/docs/mgnify-proteins-api.html), [general API](https://docs.mgnify.org/src/docs/api.html).
- [IMG/VR sequence release files](https://genome.jgi.doe.gov/portal/IMG_VR/IMG_VR.home.html).
- [Logan official resource](https://github.com/IndexThePlanet/Logan).
- [Pyrodigal-gv](https://github.com/althonos/pyrodigal-gv), [MMseqs2 guide](https://github.com/soedinglab/MMseqs2/wiki), [NCBI BLAST API](https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html).

Online resource documentation checked 20 October 2026. A documented service or mapping field does not guarantee recoverability of every candidate.
