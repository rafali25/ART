# MarsHill ART RT profile pilot

## Result

The profile retrieved 1 candidate domain hit(s) from the 262 annotated MarsHill proteins at the 1e-5 sequence and conditional-domain reporting thresholds.

| Protein | Genome | CDS coordinates (1-based, inclusive) | Length | E-value | Bits | Profile coverage |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| QQM14740.1 | MW248466.1 | 215089–216594 (strand 1) | 501 | 0 | 1060.7 | 100.0% |

## Training and provenance

Five distinct full-length RT seeds were used: LY01, PALS_2, S6, Madawaska and UFV_DC4. MarsHill QQM14740.1 and its exact sequence were excluded from training.
The full-length global identities among seeds and target span 85.6%–96.4%. This is a related-sequence recovery check, not a lineage-held-out sensitivity estimate.
Profile length: 500 match states. Training sequence count: 5.
Each curated RT translates exactly from its annotated genomic CDS; each has five matches to the known type-II array core within one mismatch in its 6-kb upstream window, and an adjacent co-oriented 274-aa candidate partner. This context check does not replace array delimitation, functional annotation or binding experiments.

## Reference checks

The new profile returned 56 reported sequence hits against 2339 historical myRT reference sequences at E <= 1e-5. The panel includes retron, DGR, group-II intron, CRISPR-associated and other RT labels. See validation_summary.json and the raw outputs.
Of these reference hits, 0 have a reported domain covering at least 60% of the new profile. This proposed coverage screen separates this particular reference panel from the MarsHill recovery, but does not establish a family-wide decision rule or sensitivity/specificity estimate.
Competing myRT-class profile matches to the retrieved MarsHill protein:

- RVT-CRISPR-G1: E=2e-10, score=26.6, coverage=77.9%.
- RVT-CRISPR-G2: E=6.5e-17, score=47.9, coverage=96.6%.
- RVT-CRISPR-G5: E=4.6e-12, score=31.9, coverage=96.6%.
- RVT-DGRs: E=2.9e-14, score=39.3, coverage=67.3%.
- RVT-G2L4: E=3e-07, score=16.2, coverage=26.9%.
- RVT-G2Lb: E=1.4e-07, score=17.6, coverage=55.7%.
- RVT-G2L: E=7.4e-08, score=18.3, coverage=71.4%.
- RVT-GII-I: E=2.9e-11, score=29.3, coverage=96.4%.
- RVT-GII-II: E=1.9e-18, score=53.1, coverage=98.1%.
- RVT-Retrons: E=2.4e-35, score=108.0, coverage=99.0%.
- RVT-UG11: E=8.7e-10, score=24.6, coverage=16.3%.
- RVT-UG19: E=7.3e-13, score=34.9, coverage=66.9%.
- RVT-UG7: E=3.3e-12, score=32.7, coverage=93.7%.

The relaxed 1e-3 sequence reporting threshold returned 1 protein(s).

## Interpretation limits

- This is a Staphylococcus/type-II ART profile pilot. Distant type-I/III ART recovery has not been tested.
- Curated NART-specific validation remains pending; historical myRT labels are not an independently curated ART/NART negative set.
- No calibrated ART bit-score threshold or sensitivity/specificity estimate is claimed. The paper's 207-bit threshold is not used.
- Only annotated MarsHill CDSs were searched; missed or miscalled genes remain possible.
- The SA1 FASTA was not searched. The Madawaska GenBank file in the SA1 folder is tracked as MW349129.1.

An E-value printed as 0 is HMMER numerical underflow, not a mathematical zero.

Actual commands, reference versions, input hashes and raw logs accompany these results.
