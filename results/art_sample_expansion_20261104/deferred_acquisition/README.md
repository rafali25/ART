# Deferred MGnify and archived partner follow-up, 4 November 2026

This bounded follow-up recovered **one new complete provisional RT–type I partner pair**, `SRR7532897_75__SRR7532897_75_pgv1`. Its RT and primary partner differ from all existing exported sequences, but it belongs to an existing partner cluster. Combining it with the 13-pair baseline gives 14 candidate pairs and still 12 clusters at 70% identity and 80% mutual coverage. This branch does not modify the baseline catalog.

The source is a current Logan contig of 9,827 nt. The old archived partner label `SRR7532897_76_2` was mapped by exact genomic translation to current contig `SRR7532897_75`; the historical contig number was not assumed to remain valid.

| Component | Source interval, zero-based/half-open | Strand | Protein length | Evidence |
|---|---|---|---:|---|
| RT, `SRR7532897_75_pgv1` | [66,1890) | + | 607 aa | ART RT profile E=8.2e-51; coverage 0.876 |
| Primary partner, `SRR7532897_75_pgv2` | [1870,4426) | + | 851 aa | Original type I profile E=5.8e-38, coverage 0.626; two-seed profile E=1.7e-79, coverage 0.700 |

Both Pyrodigal-gv CDS translations match their genomic DNA exactly and satisfy complete-CDS translation checks. They are mutual first neighbors and overlap by 20 nt. A further downstream partner homolog is retained as an alternative association and is not counted as another primary pair.

Only **66 nt upstream of the RT** are present. Array assessment is incomplete; this cannot establish array absence. No array significance test or sequence-family tree placement was performed for this new pair. Existing cohort p-values were not reused.

The nearest existing pair is the ERR10896470/SRR26805139 group. Descriptive global alignments give approximately 97.69% RT identity and 93.07% partner identity to the nearest existing sequences. Native MMseqs2 confirms that the new partner joins the SRR26805139/ERR10896470 cluster. This is an additional source and distinct protein pair, without an additional 70% partner-diversity group or an established count of independent biological samples.

## Other retrieved leads

- `ERR3946871_166` (61,637 nt) and `ERR2847180_159` (13,500 nt) encode exact matches to archived partner leads. Their partner genes pass the exploratory two-seed profile, but neither contig has an ART RT or Pfam RVT_1 match at the E<=0.001 reporting threshold in the searched gene and six-frame representations. They remain partner homolog orphans.
- Deferred MGnify partner `MGYP010182298134` maps exactly to `MGYC001898389771` [0,1536), +, on a 16,365-nt contig. It starts at the source boundary and cannot provide the required upstream RT. Its translated gene/ORF matches also fall below the frozen 0.60 profile-coverage passing rule.
- `MGYP001235072058` metadata maps to a 1,370-nt source, `MGYC000372146758`. The per-contig endpoint returned 404. This remains unresolved; metadata alone is not a recovered pair.

The native screen covered **676 gene and six-frame translated records across four source contigs**. Source records, Pyrodigal-gv models, alternative calls, six-frame maps, all reported HMM hits, gene neighborhoods and command logs are retained under `processed/`.

## Files

- `processed/primary_pair_summary.tsv`, `processed/screening_summary.json`: candidate coordinates, protein completeness, profile evidence, sequence comparisons and limitations.
- `processed/primary_RT_proteins.faa`, `processed/primary_partners.faa`, `processed/primary_NTDs.faa`, `processed/primary_pair_CDS.fna`: primary protein and genomic CDS exports. The NTD is an operational PF00078-based boundary.
- `processed/loci/`: oriented FASTA, GenBank, GFF3 and source coordinate maps.
- `processed/orphan_partner_summary.tsv`, `processed/orphan_partner_proteins.faa`: partner-only results.
- `processed/diversity/summary.json`: actual combined MMseqs2 clustering result.
- `*.download.json`, `*.request.json`, `*.scan.json`, `retrieval_summary.json`: source URLs, dates, ETags, hashes, bounded-prefix scans, exact peptide mappings and retained HTTP failures.

Successful public HTTP response bodies total **4,919,302 bytes**, below the 5,000,000-byte cap. No remote homology jobs were submitted by this branch. The three Logan downloads used bounded Range requests and searched only complete decompressed FASTA records. A missing match in a prefix would not establish biological absence.

`screen_deferred.py` adapts the existing pipeline's executable paths in memory to native Linux HMMER/MAFFT and uses the frozen original and two-LPJP1 profiles. It leaves old scripts and result cohorts unchanged. Rerun from the repository root with `.venv/bin/python results/art_sample_expansion_20261104/deferred_acquisition/screen_deferred.py`. The selected primary exports and diversity tables are derived follow-up exports; the screen itself does not finalize a combined catalog or calibrate the ART family.
