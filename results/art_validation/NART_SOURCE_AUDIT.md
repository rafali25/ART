# NART comparison-set dependency

Checked 20 October 2026. **No verified NART protein/contig mapping was retrieved. No NART specificity estimate is reported.**

The supplied Yoon PDF defines NART as the sister clade of ART, rather than any RT without an array. Supplementary Figure 5 (PDF p.19, visually inspected) shows 29 NART sequences and 17 assessable array-negative loci. Its illustrated examples are labelled only “54-kb contig”, “333-kb contig” and “14-kb contig”. Neither those labels nor the 164-aa example N-terminal domain provide an accession or sufficient sequence to reconstruct a reference.

Sources checked include the supplied PDF and archived bioRxiv full text, the downloaded myRT reference collection, the existing public-data reproduction checkout and the related ART sequence catalogue retained in the [previous source audit](../art_discovery/typeIII_source_audit/REPORT.md). Targeted online searches for NART/ART and AR9-related RT accessions did not supply a traceable mapping. These are bounded checks; they do not establish that an accession mapping is absent everywhere.

The AR9 genome KU878088.1 was downloaded from ENA and retained under `sources/`. Its deposited CDS product annotations contain no annotated reverse transcriptase. Being an AR9-like phage cannot establish NART membership. It was therefore excluded from the NART evaluation panel. The paper notes that most related polymerase phages lack the RT, making such a substitution particularly inappropriate.

The required input is a source-study NART-labelled protein FASTA linked to contigs, or the accession/coordinate mapping for its 29 phylogenetic tips. A generic retron sequence, lack of detected repeats, or RNA-polymerase adjacency alone will not be used to relabel a negative as NART.

The comparison manifest includes a NART row with status `unavailable_verified_sequence_mapping`; it is not an evaluated locus or a counted true negative. The scale-search readiness flag remains false.

Primary sources: [Yoon et al. ART preprint](https://doi.org/10.64898/2026.10.05.755809), [AR9 deposited genome](https://www.ebi.ac.uk/ena/browser/view/KU878088.1), [AR9 genome study](https://pmc.ncbi.nlm.nih.gov/articles/PMC13189113/).
