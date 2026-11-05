"""Figures for the RECOMB 2027 paper, drawn from the released dataset and the
frozen array-test outputs. Run inside a Slurm job:

  make_figures.py PROJECT_ROOT OUTDIR
"""
import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon, Rectangle

ROOT, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)

# Validated default palette (dataviz references/palette.md), light mode.
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
RT_CORE, RT_NTD = "#2a78d6", "#86b6ef"      # blue 450 / blue 250: one entity, two parts
PARTNER = "#eb6834"                          # categorical slot 2
OTHER_FILL, OTHER_EDGE = "#ecebe7", "#a9a8a2"
BG_DINUC, BG_BLOCK = "#1baf7a", "#6250d6"    # aqua / violet: null backgrounds
WINDOW = "#e3eefb"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7, "axes.labelsize": 7,
    "axes.titlesize": 7.5, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5, "axes.edgecolor": INK2, "axes.linewidth": 0.6,
    "xtick.color": INK2, "ytick.color": INK2, "axes.labelcolor": INK,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "pdf.fonttype": 42,
})

rows = {r["locus_id"]: r for r in csv.DictReader(
    open(ROOT / "data/art_samples_2026-11-05/art_samples.tsv", newline=""), delimiter="\t")}
arrays = list(csv.DictReader(
    open(ROOT / "results/art_type1_catalog/supported_array_candidates.tsv", newline=""), delimiter="\t"))

LABEL = {
    "ERR8055547_17635__ERR8055547_17635_pgv12": "ERR8055547_17635",
    "MZ422438.1__QXN67953.1": "LPJP1 QXN67953.1",
    "MGYC002687634053__MGYC002687634053_pgv52": "MGYC002687634053",
    "MGYC003526048914__MGYC003526048914_pgv179": "MGYC003526048914",
}


def find_cds(record, ident):
    for f in record.features:
        if f.type != "CDS":
            continue
        vals = {v for vs in f.qualifiers.values() for v in vs}
        if ident in vals or any(v.endswith(ident) for v in vals):
            return f
    raise KeyError(ident)


def arrow(ax, x0, x1, strand, y, h, fc, ec, lw=0.5, z=2):
    """Gene arrow from x0 to x1 (kb); strand +1 points right."""
    length = x1 - x0
    head = min(0.18, 0.35 * length)
    if strand < 0:
        pts = [(x1, y - h / 2), (x0 + head, y - h / 2), (x0, y), (x0 + head, y + h / 2), (x1, y + h / 2)]
    else:
        pts = [(x0, y - h / 2), (x1 - head, y - h / 2), (x1, y), (x1 - head, y + h / 2), (x0, y + h / 2)]
    ax.add_patch(Polygon(pts, closed=True, fc=fc, ec=ec, lw=lw, zorder=z, joinstyle="round"))


def load_nulls(details_path, locus):
    cohort_dir = (ROOT / details_path).parent.parent
    out = {}
    for bg in ("dinucleotide", "local_blocks"):
        d = json.load(open(cohort_dir / "nulls" / f"{locus}.{bg}.2000.json"))
        out[bg] = np.array([s for _, s in d["draws"]], dtype=float)
    return out


# ---------------------------------------------------------------- Figure 1
tracks = []
for a in arrays:
    locus = a["locus_id"]
    det = json.load(open(ROOT / a["array_details"]))
    rec = SeqIO.read(ROOT / a["locus_genbank"], "genbank")
    r = rows[locus]
    rt = find_cds(rec, r["RT_feature_id"])
    partner = find_cds(rec, r["partner_feature_id"])
    rt0 = int(rt.location.start)
    rel = lambda x: (x - rt0) / 1000.0
    lo, hi = float(a["locus_start0"]) - 350, int(partner.location.end) + 250
    others = [f for f in rec.features if f.type == "CDS" and f is not rt and f is not partner
              and int(f.location.end) > lo and int(f.location.start) < hi]
    tracks.append(dict(
        locus=locus, rt=(rel(rt0), rel(int(rt.location.end))),
        ntd_end=rel(rt0 + 3 * int(r["NTD_length_aa"])),
        partner=(rel(int(partner.location.start)), rel(int(partner.location.end)), partner.location.strand),
        others=[(rel(int(f.location.start)), rel(int(f.location.end)), f.location.strand) for f in others],
        repeats=[rel((rp["start0"] + rp["end0"]) / 2) for rp in det["repeats"]],
        lo=rel(lo), hi=rel(hi), copies=int(a["copies"]), p=float(a["corrected_p"]),
        observed=float(det["extended"]["score"]), nulls=load_nulls(a["array_details"], locus),
        strict=rows[locus]["category1_ART_profile_pass"] == "True",
    ))

xlo = min(t["lo"] for t in tracks)
xhi = max(t["hi"] for t in tracks)
smax = max(max(t["observed"], max(v.max() for v in t["nulls"].values())) for t in tracks) * 1.06
fig = plt.figure(figsize=(6.5, 2.5))
gs = fig.add_gridspec(len(tracks), 2, width_ratios=[2.15, 1], hspace=0.62, wspace=0.10,
                      left=0.015, right=0.985, top=0.86, bottom=0.14)
B = 2000
grid = np.linspace(0, smax, 600)
for i, t in enumerate(tracks):
    ax = fig.add_subplot(gs[i, 0])
    ax.plot([t["lo"], t["hi"]], [0, 0], color=MUTED, lw=0.8, zorder=1, solid_capstyle="round")
    for x0, x1, s in t["others"]:
        arrow(ax, x0, x1, s, 0, 0.62, OTHER_FILL, OTHER_EDGE)
    arrow(ax, t["rt"][0], t["rt"][1], 1, 0, 0.62, RT_CORE, RT_CORE, z=3)
    ax.add_patch(Rectangle((t["rt"][0], -0.31), t["ntd_end"] - t["rt"][0], 0.62,
                           fc=RT_NTD, ec=RT_CORE, lw=0.5, zorder=4))
    arrow(ax, t["partner"][0], t["partner"][1], t["partner"][2], 0, 0.62, PARTNER, PARTNER, z=3)
    ax.scatter(t["repeats"], [0] * len(t["repeats"]), marker="D", s=9, color=INK,
               edgecolor="white", linewidth=0.5, zorder=5)
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(-0.42, 0.42)
    ax.set_yticks([])
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    if i < len(tracks) - 1:
        ax.spines["bottom"].set_visible(False)
        ax.set_xticks([])
    else:
        ax.set_xlabel("position relative to RT start (kb)", labelpad=1.5)
    ax.text(xlo, 0.5, LABEL[t["locus"]], ha="left", va="bottom", fontsize=6.6, color=INK,
            fontweight="bold", clip_on=False)
    ax.text(xhi, 0.5,
            f"{t['copies']} repeat copies,  corrected $p$ = {t['p']:.4f}",
            ha="right", va="bottom", fontsize=6.0, color=INK2, clip_on=False)

    ax2 = fig.add_subplot(gs[i, 1])
    for bg, color, ls in (("dinucleotide", BG_DINUC, "-"), ("local_blocks", BG_BLOCK, "--")):
        v = t["nulls"][bg]
        surv = np.array([(v >= g).mean() for g in grid])
        keep = surv > 0
        ax2.step(grid[keep], surv[keep], where="post", color=color, lw=1.1, ls=ls)
    ax2.axvline(t["observed"], color=INK, lw=1.0)
    ax2.set_yscale("log")
    ax2.set_ylim(0.5 / B, 1.6)
    ax2.set_xlim(0, smax)
    ax2.set_yticks([1e-3, 1e-1])
    ax2.set_yticklabels(["$10^{-3}$", "$10^{-1}$"])
    ax2.tick_params(axis="y", which="minor", left=False)
    for side in ("right", "top"):
        ax2.spines[side].set_visible(False)
    if i < len(tracks) - 1:
        ax2.set_xticklabels([])
    else:
        ax2.set_xlabel("best-chain score $s$", labelpad=1.5)

fig.text(0.004, 0.965, "a", fontsize=8, fontweight="bold", color=INK)
fig.text(0.676, 0.965, "b", fontsize=8, fontweight="bold", color=INK)
arch = [
    Line2D([], [], marker="D", ls="", color=INK, markersize=3.5, label="repeat copy"),
    Patch(fc=RT_NTD, ec=RT_CORE, lw=0.5, label="RT NTD"),
    Patch(fc=RT_CORE, ec=RT_CORE, label="RT core"),
    Patch(fc=PARTNER, ec=PARTNER, label="type I partner"),
    Patch(fc=OTHER_FILL, ec=OTHER_EDGE, lw=0.5, label="other CDS"),
]
nulls = [
    Line2D([], [], color=BG_DINUC, lw=1.1, label="dinucleotide"),
    Line2D([], [], color=BG_BLOCK, lw=1.1, ls="--", label="local block"),
    Line2D([], [], color=INK, lw=1.0, label="observed"),
]
fig.legend(handles=arch, loc="upper left", bbox_to_anchor=(0.02, 1.005), ncol=5, frameon=False,
           handlelength=1.4, columnspacing=1.0, handletextpad=0.4)
fig.legend(handles=nulls, loc="upper left", bbox_to_anchor=(0.69, 1.005), ncol=3, frameon=False,
           handlelength=1.3, columnspacing=0.55, handletextpad=0.3, fontsize=6.0)
fig.savefig(OUT / "fig_loci_nulls.pdf")
fig.savefig(OUT / "fig_loci_nulls.png", dpi=200)
plt.close(fig)

# ---------------------------------------------------------------- Figure 2 (appendix)
pairs = list(csv.DictReader(open(ROOT / "results/art_expansion_20261105/step5/matching_pair_comparisons.tsv",
                                 newline=""), delimiter="\t"))
x = np.array([float(p["RT_core_identity"]) for p in pairs])
y = np.array([float(p["partner_identity"]) for p in pairs])
ok = np.array([p["passes_sequence_matching"] == "True" for p in pairs])
inwin = (x >= 0.5) & (x <= 0.9) & (y >= 0.5) & (y <= 0.9)

order = ["Type_I", "Type_II"]
samples = sorted(rows.values(), key=lambda r: (order.index(r["partner_type"]), r["sample_id"]))
fields = [
    ("retained in this expansion", lambda r: r["new_in_20261105_expansion"] == "True"),
    ("complete CDS pair", lambda r: r["pair_complete"] == "True"),
    ("ART RT profile screen", lambda r: r["category1_ART_profile_pass"] == "True"),
    ("RVT_1 + [YF]xDD", lambda r: r["category3_RVT1_domain"] == "True" and r["category3_catalytic_motif_pattern"] == "True"),
    ("NTD $\\geq$150 aa", lambda r: r["category4_NTD_at_least_150aa"] == "True"),
    ("partner anchor", lambda r: r["category5_direct_partner_and_500nt"] == "True"),
    ("full 6 kb upstream", lambda r: int(float(r["available_upstream_window_nt"] or 0)) >= 6000),
    ("supported array", lambda r: r["category6_historical_supported_array"] == "True"),
]
M = np.array([[f(r) for r in samples] for _, f in fields], dtype=bool)

fig = plt.figure(figsize=(6.5, 2.35))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 2.75], wspace=0.62, left=0.07, right=0.995, top=0.84, bottom=0.17)
ax = fig.add_subplot(gs[0, 0])
ax.add_patch(Rectangle((0.5, 0.5), 0.4, 0.4, fc=WINDOW, ec="none", zorder=0))
ax.scatter(x[~inwin], y[~inwin], s=9, color=MUTED, edgecolor="white", linewidth=0.4, label="outside window", zorder=2)
if (inwin & ~ok).any():
    ax.scatter(x[inwin & ~ok], y[inwin & ~ok], s=14, facecolor="white", edgecolor=PARTNER, linewidth=1.0,
               label="in window, fails coverage/length", zorder=3)
ax.scatter(x[ok], y[ok], s=22, color=RT_CORE, edgecolor="white", linewidth=0.6, label="eligible edge", zorder=4)
ax.text(0.66, 0.66, "eligible\nwindow", ha="center", va="center", fontsize=5.8, color=INK2)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_aspect("equal")
ax.set_xlabel("RT-core identity")
ax.set_ylabel("partner identity")
ax.set_title(f"a   {len(pairs)} pairs of 18 strict loci", loc="left", fontsize=7, color=INK, pad=4)
for side in ("right", "top"):
    ax.spines[side].set_visible(False)
ax.legend(loc="upper left", frameon=False, ncol=1, handletextpad=0.2, borderaxespad=0.1,
          fontsize=5.4, labelspacing=0.3)

ax = fig.add_subplot(gs[0, 1])
present = np.where(M, 1.0, 0.0)
from matplotlib.colors import ListedColormap
ax.imshow(present, aspect="auto", cmap=ListedColormap(["#ecebe7", "#256abf"]), interpolation="nearest")
ax.set_yticks(range(len(fields)))
ax.set_yticklabels([f for f, _ in fields])
ax.set_xticks(range(len(samples)))
ax.set_xticklabels([r["sample_id"].replace("ARTS00", "") for r in samples], rotation=90, fontsize=4.8)
ax.set_xticks(np.arange(-0.5, len(samples)), minor=True)
ax.set_yticks(np.arange(-0.5, len(fields)), minor=True)
ax.grid(which="minor", color="white", linewidth=1.0)
ax.tick_params(which="minor", length=0)
ax.tick_params(axis="both", length=0)
n1 = sum(r["partner_type"] == "Type_I" for r in samples)
ax.axvline(n1 - 0.5, color=INK, lw=1.2)
ax.text((n1 - 1) / 2, -0.75, f"type I ({n1})", ha="center", fontsize=6.5, color=INK)
ax.text(n1 + (len(samples) - n1 - 1) / 2, -0.75, f"type II ({len(samples) - n1})", ha="center", fontsize=6.5, color=INK)
ax.set_title("b   evidence fields per sample (ARTS identifier; blue = present)", loc="left", fontsize=7,
             color=INK, pad=13)
for side in ax.spines.values():
    side.set_visible(False)
fig.savefig(OUT / "fig_census.pdf")
fig.savefig(OUT / "fig_census.png", dpi=200)
plt.close(fig)

print("in-window pairs:", int(inwin.sum()), "passing:", int(ok.sum()), "of", len(pairs))
print("column totals type I:", M[:, :n1].sum(axis=1).tolist())
for t in tracks:
    print(LABEL[t["locus"]], "obs", round(t["observed"], 1),
          {bg: int((v >= t["observed"]).sum()) for bg, v in t["nulls"].items()},
          "others", len(t["others"]), "ntd_kb", round(t["ntd_end"] - t["rt"][0], 2))
