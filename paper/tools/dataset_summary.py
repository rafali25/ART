import csv, sys, collections
rows = list(csv.DictReader(open(sys.argv[1], newline=""), delimiter="\t"))
cols = ["sample_id","locus_id","partner_type","cohort","new_in_20261105_expansion","pair_complete",
        "category1_ART_profile_pass","category4_NTD_at_least_150aa","NTD_length_aa",
        "category5_direct_partner_and_500nt","upstream_annotation_free_nt","available_upstream_window_nt",
        "category6_historical_supported_array","array_status","historical_array_scope","partner_cluster70","source_identifier"]
for r in rows:
    print("\t".join(str(r.get(c,""))[:40] for c in cols))
print()
for c in ["cohort","record_role","array_status","historical_array_scope","partner_type","source_file_format"]:
    print(c, dict(collections.Counter(r[c] for r in rows)))
t1 = [r for r in rows if r["partner_type"]=="Type_I"]
short = [r for r in t1 if int(float(r["available_upstream_window_nt"] or 0)) < 6000]
print("typeI with available upstream <6000:", len(short))
src = collections.Counter(r["source_identifier"][:3] for r in rows)
print("source prefixes", dict(src))
print("distinct source groups", len({r["source_group"] for r in rows}), "typeI", len({r["source_group"] for r in t1}))
