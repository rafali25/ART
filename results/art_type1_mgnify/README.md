# MGnify Type I follow-up

See [REPORT.md](REPORT.md), [primary_typeI_pairs.tsv](primary_typeI_pairs.tsv)
and [candidate_loci.tsv](candidate_loci.tsv). The cross-cohort inventory is
in `results/art_type1_catalog/`.

The four completed phmmer searches returned 418 collected query-hit rows
(293 unique MGYP proteins). Each partner query returned nine hits. Both RT
queries were capped at 200 collected hits. The fixed acquisition panel retained
the first ten retrievable contigs under the documented query-balanced ranking.
Thirteen searched proteins mapped exactly to these contigs. Remaining hits
are explicitly deferred, not negative.

This follow-up has its own array-testing cohort. The original six-locus pilot
remains in `results/art_type1_public/`. The combined catalog additionally corrects
array tests across unique RT loci in both cohorts, without counting LPJP1 twice.

From the project root, run these stages **sequentially**:

```powershell
.venv/Scripts/python.exe scripts/art_type1_pipeline.py process --config config/art_type1_mgnify.json
.venv/Scripts/python.exe scripts/art_type1_pipeline.py nulls --config config/art_type1_mgnify.json
.venv/Scripts/python.exe scripts/art_type1_tree.py --config config/art_type1_mgnify.json
.venv/Scripts/python.exe scripts/art_type1_pipeline.py diversity --config config/art_type1_mgnify.json
.venv/Scripts/python.exe scripts/finalize_art_type1_public.py --config config/art_type1_mgnify.json
.venv/Scripts/python.exe scripts/verify_art_type1_public.py --config config/art_type1_mgnify.json
.venv/Scripts/python.exe scripts/build_art_type1_catalog.py
```

All tracts with three or more inferred repeats receive 2,000 shuffles per
background; fewer-copy outputs stop at 200. Both backgrounds rerun the full
seed-discovery and chain-selection search. Array support is not functional ART
validation. Boundary, copy-number, family-calibration and NART limitations from
the pilot still apply.
