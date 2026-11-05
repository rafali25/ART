#!/usr/bin/env python3
"""Regenerate every supplementary figure for the ART RECOMB paper from the raw result files.

Usage: python make_figures.py [fig ids...]      e.g.  python make_figures.py calibration binding
Nothing is hard-coded except Table 1 of the paper (its raw harness output is not in the results
folder) and the acquisition counts, which are read from their summary files where one exists.
Binding figures read the run folders directly, so re-running after the batch finishes updates them.
Requires: pip install numpy pandas matplotlib biopython
"""
import glob, gzip, json, math, os, re, sys
from collections import Counter, defaultdict
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle, FancyArrowPatch, PathPatch
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
from matplotlib.transforms import Affine2D
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap

ROOT = "/Users/ebadshabbir/Desktop/ART-BENCH/illustrations"
RES, DATA = f"{ROOT}/results", f"{ROOT}/data"
OUT = os.path.dirname(os.path.abspath(__file__))

# Palette carried over from Figures 1-2 of the paper, so new panels read as the same system.
NTD, CORE, PARTNER = "#8AB8F0", "#2A76D2", "#EE6C36"
OTHER, OTHER_E = "#EDECE6", "#A9A89F"
DINUC, BLOCK = "#14A37F", "#5B4FD6"
INK, MUTED, FAINT, GRID = "#1F2328", "#5F6670", "#A7ADB5", "#E7E9EC"
GOOD, WARN, BAD = "#2A76D2", "#C98A12", "#CF3A32"
SHADE = "#DCEAFB"

mpl.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
    "xtick.labelsize": 6.3, "ytick.labelsize": 6.3, "legend.fontsize": 6.2,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 400, "figure.dpi": 150,
    "mathtext.fontset": "dejavusans", "hatch.linewidth": 0.35,
})

def save(fig, name):
    fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(f"{OUT}/{name}.png", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("wrote", name)

def tag(ax, letter, x=-0.02, y=1.03, **kw):
    ax.text(x, y, letter, transform=ax.transAxes, fontweight="bold", fontsize=8.5, va="bottom", ha="right", **kw)

def tsv(path, **kw): return pd.read_csv(path, sep="\t", **kw)

SHORT = {"LC680885.1__BDE75684.1": "S6", "MN091626.1__QDJ97640.1": "PALS2", "MW218148.1__QPI16926.1": "SA1",
         "MW248466.1__QQM14740.1": "MarsHill", "MW349129.1__QQO92858.1": "Madawaska", "MZ779063.1__UAJ17040.1": "UFV_DC4",
         "OR836606.1__WRM43423.1": "LY01", "MZ422438.1__QXN67953.1": "LPJP1 QXN67953.1", "MZ422438.1__QXN67955.1": "LPJP1 QXN67955.1"}

def null_scores(path):
    d = json.load(open(path))
    return {k: np.array([x[1] for x in v], float) for k, v in d["scores"].items()}

# ------------------------------------------------------------------ calibration
def cross_background():
    """Every stored null draw is re-used as a pseudo-observed score, tested within and across backgrounds."""
    alphas = np.array([0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2])
    recs = []
    for f in sorted(glob.glob(f"{RES}/art_validation/nulls/*.json")):
        s = null_scores(f); a, b = s["dinucleotide"], s["local_blocks"]
        if a.max() == 0 and b.max() == 0: continue          # degenerate truncated / empty windows
        def pvals(x, ref, same):
            # (1 + #{ref >= x}) / (1 + |ref|); the draw itself is left out when it comes from ref
            srt = np.sort(ref); ge = len(ref) - np.searchsorted(srt, x, side="left")
            return (ge if same else 1 + ge) / (len(ref) if same else 1 + len(ref))
        for src_, oth, name in ((a, b, "dinucleotide"), (b, a, "local_blocks")):
            within, across = pvals(src_, src_, True), pvals(src_, oth, False)
            cons = np.maximum(within, across)
            for al in alphas:
                recs.append(dict(case=os.path.basename(f)[:-5], source=name, alpha=al, within=(within <= al).mean(),
                                 across=(across <= al).mean(), conservative=(cons <= al).mean()))
    return pd.DataFrame(recs), alphas

def fig_calibration():
    fig = plt.figure(figsize=(6.5, 4.2))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.3, 1, 1], hspace=0.56, wspace=0.46, left=0.135, right=0.995, top=0.95, bottom=0.095)
    ax = fig.add_subplot(gs[:, 0]); tag(ax, "a", x=-0.26, y=1.0)
    sc = tsv(f"{RES}/art_validation/sequence_cases.tsv")
    rr = tsv(f"{RES}/art_validation/reference_recovery.tsv")
    assess = tsv(f"{RES}/art_arrays_assessed/array_assessments.tsv").set_index("locus_id")
    tests = tsv(f"{RES}/art_arrays_assessed/background_tests.tsv").set_index("locus_id")
    rows = []
    for name, g in sc.groupby("test", sort=False):
        pool = np.concatenate([np.concatenate(list(null_scores(f"{RES}/art_validation/nulls/{c}.json").values())) for c in g.case_id])
        rows.append((name.replace("_", " "), pool, g.best_search_score.values, "b", None))
    for loc in rr.locus_id:
        n = null_scores(f"{RES}/art_arrays_assessed/nulls/{loc}.json")
        pool = np.concatenate([n["local_blocks_6000"], n["dinucleotide_6000"]])        # inferential 6-kb window only
        reported, top = assess.loc[loc, "seed_chain_score"], tests.loc[loc, "local_blocks_6000_best_score"]
        rows.append((SHORT[loc], pool, np.array([reported]), "r", top if top > reported + 1e-6 else None))
    y = 0; ticks = []; labs = []
    for r in [r for r in rows if r[3] == "b"] + [None] + [r for r in rows if r[3] == "r"]:
        if r is None: y += 0.8; continue
        name, pool, obs, kind, top = r
        p = pool[pool > 0]
        if len(p):
            q = np.percentile(p, [5, 50, 95])
            ax.plot([q[0], p.max()], [y, y], color=FAINT, lw=0.8, solid_capstyle="butt", zorder=1)
            ax.plot([q[0], q[2]], [y, y], color=MUTED, lw=2.6, solid_capstyle="butt", zorder=2)
            ax.plot([q[1]], [y], marker="|", ms=3.6, color="white", mew=0.8, zorder=3)
        o = obs[obs > 0]
        if len(o): ax.scatter(o, [y] * len(o), s=11, marker="D", color=CORE if kind == "b" else PARTNER, edgecolor="white", linewidth=0.3, zorder=4)
        if top: ax.scatter([top], [y], s=10, marker="D", facecolor="white", edgecolor=PARTNER, linewidth=0.7, zorder=4)
        ticks.append(y); labs.append(name); y += 1
    ax.set_yticks(ticks); ax.set_yticklabels(labs, fontsize=5.8); ax.set_ylim(y - 0.4, -0.7)
    ax.set_xscale("log"); ax.set_xlim(12, 1500); ax.set_xlabel("best-chain score $S$ (log scale)")
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.grid(axis="x", color=GRID, lw=0.5); ax.set_axisbelow(True); ax.tick_params(axis="y", length=0)
    ax.set_title("observed score vs. search null", loc="left", pad=3)
    ax.legend(handles=[Line2D([], [], color=MUTED, lw=2.6, label="null 5\u201395%"), Line2D([], [], color=FAINT, lw=0.9, label="up to null max"),
                       Line2D([], [], marker="D", ls="", color=CORE, ms=3, label="designed"), Line2D([], [], marker="D", ls="", color=PARTNER, ms=3, label="reference")],
              loc="center right", handlelength=1.3, borderpad=0.1, labelspacing=0.3, fontsize=5.6, bbox_to_anchor=(1.0, 0.31))
    # b -- DGR: the two backgrounds disagree
    ax = fig.add_subplot(gs[0, 1]); tag(ax, "b", x=-0.2)
    s = null_scores(f"{RES}/art_validation/nulls/DGR__NP_958675.1.json")
    obs = float(sc.set_index("case_id").loc["DGR__NP_958675.1", "best_search_score"])
    for k, c, ls in (("dinucleotide", DINUC, "-"), ("local_blocks", BLOCK, "--")):
        x = np.sort(s[k]); yy = 1 - np.arange(len(x)) / len(x)
        ax.step(x, yy, where="post", color=c, ls=ls, lw=1.0)
    ax.axvline(obs, color=INK, lw=0.9); ax.set_yscale("log"); ax.set_ylim(8e-4, 1.3); ax.set_xlim(20, 100)
    pd_ = (1 + (s["dinucleotide"] >= obs).sum()) / 1001; pb = (1 + (s["local_blocks"] >= obs).sum()) / 1001
    ax.text(obs + 2, 0.5, f"observed\n$S$ = {obs:.1f}", fontsize=5.8, va="center")
    ax.text(obs + 2.2, pd_ * 1.05, f"dinucleotide\n$\\hat p$ = {pd_:.3f}", fontsize=5.8, color=DINUC, va="bottom", linespacing=1.05)
    ax.text(obs - 2, pb * 0.55, f"local block\n$\\hat p$ = {pb:.3f}", fontsize=5.8, color=BLOCK, va="center", ha="right", linespacing=1.05)
    ax.set_xlabel("best-chain score $s$"); ax.set_ylabel("fraction of draws $\\geq s$")
    ax.set_title("DGR control", loc="left", pad=3)
    # c -- cross-background calibration, as a ratio so anticonservatism is readable
    ax = fig.add_subplot(gs[0, 2]); tag(ax, "c", x=-0.2)
    cal, alphas = cross_background()
    m = cal.groupby("alpha")[["within", "across", "conservative"]]
    mean, hi = m.mean(), m.quantile(0.975)
    worst = m.max()
    ax.axhline(1, color=FAINT, lw=0.7, zorder=1)
    ax.fill_between(alphas, (mean["across"] / alphas), hi["across"] / alphas, color=BAD, alpha=0.13, lw=0, zorder=0)
    for col, c, ls, lab, mk in (("across", BAD, "--", "other background", "v"), ("within", MUTED, "-", "own background", "o"), ("conservative", CORE, "-", "max of the two", "s")):
        ax.plot(alphas, mean[col] / alphas, color=c, ls=ls, lw=1.0, marker=mk, ms=2.4, label=lab, zorder=3)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(8e-4, 0.25); ax.set_ylim(0.2, 16)
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v:g}\u00d7"))
    ax.set_xticks([1e-3, 1e-2, 1e-1]); ax.set_xlabel("nominal $\\alpha$"); ax.set_ylabel("rejection rate / $\\alpha$")
    ax.text(1.1e-3, hi["across"].iloc[0] / alphas[0] * 1.12, "97.5th percentile of windows", fontsize=5.3, color=BAD, va="bottom")
    ax.set_title("rejections on null draws", loc="left", pad=3); ax.legend(loc="lower right", handlelength=1.8, borderpad=0.1, labelspacing=0.15, fontsize=5.3)
    # d -- corrected floor against family size
    ax = fig.add_subplot(gs[1, 1]); tag(ax, "d", x=-0.2)
    ms = np.array([1.8, 170])
    for B, ls, xl in ((200, ":", 40), (1000, "--", 105), (2000, "-", 62)):
        ax.plot(ms, ms / (B + 1), color=MUTED, ls=ls, lw=0.9)
        ax.text(xl, xl / (B + 1), f"$B$={B:,}", fontsize=5.3, color=MUTED, ha="center", va="center", zorder=3, bbox=dict(fc="white", ec="none", pad=0.5))
    ax.axhline(0.05, color=BAD, lw=0.7); ax.text(2.1, 0.056, "0.05", color=BAD, fontsize=5.6, va="bottom", ha="left")
    for mm, BB, lab, xt, yt, ha in ((9, 1000, "9 refs", 7.8, 0.0125, "right"), (17, 2000, "17 type I", 21, 0.0058, "left"), (30, 1000, "30 cases", 25, 0.041, "right")):
        ax.scatter([mm], [mm / (BB + 1)], s=14, color=CORE, edgecolor="white", linewidth=0.4, zorder=4)
        ax.text(xt, yt, lab, fontsize=5.5, ha=ha, va="center")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1.8, 170); ax.set_ylim(1.5e-3, 0.8)
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_xlabel("declared family size $m$"); ax.set_ylabel("floor $m/(B+1)$")
    ax.set_title("attainable corrected $p$", loc="left", pad=3)
    # e -- references change status with the family
    ax = fig.add_subplot(gs[1, 2]); tag(ax, "e", x=-0.2)
    cmp_ = tsv(f"{RES}/array_calibration_20261105/pilercr_shuffle_wrapper/comparison.tsv").set_index("input_id")
    refs = cmp_[cmp_.group == "published_reference"]
    mm = np.array([1, 50.0])
    for l, r in refs.iterrows():
        ax.plot(mm, np.minimum(1, mm * r.D_conservative_p), color=CORE, lw=0.75, alpha=0.8)
    reported = assess.loc["MW248466.1__QQM14740.1", "max_background_p"]
    ax.plot(mm, np.minimum(1, mm * reported), color=CORE, lw=0.75, ls=(0, (2, 1.5)))
    ax.axhline(0.05, color=BAD, lw=0.7); ax.text(1.15, 0.056, "0.05", color=BAD, fontsize=5.6, va="bottom")
    n9 = (9 * refs.D_conservative_p <= 0.05).sum(); n31 = (31 * refs.D_conservative_p <= 0.05).sum()
    for xx, nn in ((9, n9), (31, n31)):
        ax.axvline(xx, color=FAINT, lw=0.6, ls=":", zorder=0)
        ax.text(xx, 0.0009, f"$m$={xx}: {nn}/7", fontsize=5.5, ha="right", va="bottom", color=INK, rotation=90, bbox=dict(fc="white", ec="none", pad=0.3))
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1, 50); ax.set_ylim(7e-4, 1.1)
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_xlabel("declared family size $m$"); ax.set_ylabel("corrected $\\tilde p$")
    ax.set_title("7 references vs. family size", loc="left", pad=3)
    save(fig, "fig_calibration")
    out = mean.add_prefix("mean_").join(hi.add_prefix("p97.5_")).join(worst.add_prefix("worst_"))
    out.to_csv(f"{OUT}/source_data_calibration.csv"); print(out.round(5).to_string())

# ------------------------------------------------------------------ benchmark ladder
def fig_benchmark():
    sc = tsv(f"{RES}/art_validation/sequence_cases.tsv"); truth = tsv(f"{RES}/art_validation/synthetic_truth.tsv").set_index("case_id")
    fig, ax = plt.subplots(figsize=(6.5, 4.4)); fig.subplots_adjust(left=0.165, right=0.70, top=0.945, bottom=0.075)
    y = 0; last = None; ylab = []; ypos = []
    txt = lambda yy, xa, s, **k: ax.text(xa, yy, s, transform=ax.get_yaxis_transform(), va="center", fontsize=5.9, clip_on=False, **k)
    columns = [(1.05, "planted"), (1.16, "found"), (1.255, "extra"), (1.355, "$\\tilde p$"), (1.475, "behaviour")]
    for c in sc.itertuples():
        if last is not None and c.test != last: y += 0.75
        last = c.test
        planted = json.loads(truth.loc[c.case_id, "repeat_intervals_json"]) if c.case_id in truth.index else []
        det = json.loads(c.detected_intervals_json) if isinstance(c.detected_intervals_json, str) else []
        ok = bool(c.pass_expected_behavior)
        if not ok: ax.axhspan(y - 0.5, y + 0.5, color="#FBE9E7", lw=0, zorder=0)
        ax.plot([0, 6000], [y, y], color=GRID, lw=0.6, zorder=1)
        for s, e in planted: ax.add_patch(Rectangle((s, y - 0.42), max(e - s, 20), 0.42, fc=FAINT, ec="none", zorder=2))
        pset = [tuple(p) for p in planted]; extra = 0
        for s, e in det:
            match = any(s < pe and e > ps for ps, pe in pset)
            ax.add_patch(Rectangle((s, y), max(e - s, 20), 0.42, fc=CORE if match else BAD, ec="none", zorder=3))
            if not match and planted:
                extra += 1; ax.scatter([s + 10], [y + 0.62], marker="^", s=9, color=BAD, zorder=5, clip_on=False, linewidth=0)
        vals = [len(planted) if planted else "\u2013", len(det), extra if planted else "\u2013", f"{c.selected_adjusted_p:.3f}" if c.selected_adjusted_p < 1 else "1.0"]
        for i, ((xa, _), v) in enumerate(zip(columns[:4], vals)):
            txt(y, xa, str(v), ha="center", color=BAD if (i == 2 and v not in ("\u2013", 0)) else INK)
        txt(y, columns[4][0], "pass" if ok else "fail", ha="center", fontweight="normal" if ok else "bold", color=MUTED if ok else BAD)
        ypos.append(y); ylab.append(c.case_id.split("__")[0].replace("_", " ")); y += 1
    for xa, s in columns: ax.text(xa, -1.2, s, transform=ax.get_yaxis_transform(), ha="center", va="bottom", fontsize=5.9, fontweight="bold", clip_on=False)
    ax.set_yticks(ypos); ax.set_yticklabels(ylab, fontsize=5.7); ax.tick_params(axis="y", length=0)
    ax.set_ylim(y - 0.2, -0.9); ax.set_xlim(0, 6000); ax.set_xlabel("position in the assessed window (nt)"); ax.spines["left"].set_visible(False)
    h = [Rectangle((0, 0), 1, 1, fc=FAINT), Rectangle((0, 0), 1, 1, fc=CORE), Line2D([], [], marker="^", ls="", color=BAD, ms=3.5)]
    ax.legend(h, ["planted copy", "detected, matches a planted copy", "detected, spurious"], ncol=3, loc="lower left", bbox_to_anchor=(-0.02, 1.002), handlelength=1.1, columnspacing=1.0, borderpad=0)
    save(fig, "fig_benchmark")

# ------------------------------------------------------------------ baselines + wrapper
# Table 1 of the paper; the raw baseline harness output is not in the results folder.
TABLE1 = {("PILER-CR", "default"): (0.231, 0.952, 0.952), ("PILER-CR", "tuned"): (0.812, 0.876, 0.900),
          ("CRT", "default"): (0.231, 1.0, 1.0), ("CRT", "tuned"): (0.231, 1.0, 1.0),
          ("MinCED", "default"): (0.231, 1.0, 1.0), ("MinCED", "tuned"): (0.231, 1.0, 1.0), ("This work", ""): (0.923, 0.988, 0.988)}

def fig_baselines():
    fig = plt.figure(figsize=(6.5, 2.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.62], wspace=0.34, left=0.075, right=0.885, top=0.9, bottom=0.17)
    ax = fig.add_subplot(gs[0, 0]); tag(ax, "a", x=-0.14)
    ax.set_xlim(0, 1.04); ax.set_ylim(0.855, 1.012); ax.grid(color=GRID, lw=0.5); ax.set_axisbelow(True)
    p = TABLE1
    ax.annotate("", xy=p[("PILER-CR", "tuned")][:2], xytext=p[("PILER-CR", "default")][:2], arrowprops=dict(arrowstyle="-|>", color=PARTNER, lw=0.9, shrinkA=3, shrinkB=3, mutation_scale=7))
    ax.scatter(*p[("PILER-CR", "default")][:2], s=20, facecolor="white", edgecolor=PARTNER, linewidth=1, zorder=4)
    ax.scatter(*p[("PILER-CR", "tuned")][:2], s=20, color=PARTNER, zorder=4)
    ax.scatter(p[("PILER-CR", "tuned")][0], p[("PILER-CR", "tuned")][2], s=16, marker="_", color=PARTNER, linewidth=1.2, zorder=4)
    ax.scatter(*p[("CRT", "default")][:2], s=20, color=MUTED, marker="s", zorder=4)
    ax.scatter(*p[("This work", "")][:2], s=42, color=CORE, marker="D", edgecolor="white", linewidth=0.5, zorder=5)
    ax.text(0.27, 1.0, "CRT, MinCED\n(default = tuned)", fontsize=5.8, va="center", linespacing=1.05)
    ax.text(0.2, 0.944, "PILER-CR\ndefault", fontsize=5.8, va="top", ha="center", linespacing=1.05)
    ax.text(0.785, 0.872, "PILER-CR tuned", fontsize=5.8, va="center", ha="right")
    ax.text(0.84, 0.900, "lenient", fontsize=5.2, va="center", ha="left", color=MUTED)
    ax.text(0.885, 0.984, "this work", fontsize=6.2, va="center", ha="right", color=CORE, fontweight="bold")
    ax.set_xlabel("copy recall"); ax.set_ylabel("copy precision"); ax.set_title("search: 31 inputs, 260 copies", loc="left", pad=3)
    # b -- wrapper, one column per input
    ax = fig.add_subplot(gs[0, 1]); tag(ax, "b", x=-0.03)
    w = tsv(f"{RES}/array_calibration_20261105/pilercr_shuffle_wrapper/comparison.tsv")
    w["fam"] = w.group.str.replace("designed_", "").str.replace("_", " ")
    thr = 0.05 / 31
    x = 0; xs = []; lastf = None
    for r in w.itertuples():
        if lastf is not None and r.fam != lastf: x += 0.6
        lastf = r.fam; xs.append(x); x += 1
    w["x"] = xs
    def mark(a, xv, yv, p_):
        if p_ <= thr: a.scatter([xv], [yv], s=21, marker="s", color=CORE, linewidth=0, zorder=3)
        elif p_ <= 0.05: a.scatter([xv], [yv], s=19, marker="s", facecolor="white", edgecolor=CORE, linewidth=0.8, zorder=3)
        else: a.scatter([xv], [yv], s=12, marker="x", color=BAD, linewidth=0.8, zorder=3)
    for r in w.itertuples():
        mark(ax, r.x, 1, r.D_conservative_p); mark(ax, r.x, 0, r.pilercr_conservative_p)
    ax.set_yticks([1, 0]); ax.set_yticklabels(["selection-aware $D$", "tuned PILER-CR,\nsame wrapper"], fontsize=6.2, linespacing=1.05); ax.set_ylim(-1.75, 1.75); ax.tick_params(axis="y", length=0)
    fams = w.groupby("fam", sort=False).x.agg(["min", "max"])
    labels = {"long spacer 120": "120", "long spacer 180": "180", "long spacer 220": "220", "imperfect copies": "imp.", "missing copy": "miss.", "extended beyond 3kb": ">3 kb", "unrelated RT with array": "unrel.", "repetitive DNA": "rep.", "published reference": "references"}
    for f, r_ in fams.iterrows():
        ax.text((r_["min"] + r_["max"]) / 2, -0.62, labels[f], fontsize=5.3, ha="center", va="center", color=MUTED)
    ax.set_xticks([]); ax.spines["left"].set_visible(False); ax.spines["bottom"].set_visible(False); ax.set_xlim(-0.8, w.x.max() + 0.8)
    pub = fams.loc["published reference"]
    ax.add_patch(Rectangle((pub["min"] - 0.55, -0.85), pub["max"] - pub["min"] + 1.1, 2.5, fc=SHADE, alpha=0.6, lw=0, zorder=0))
    s_ = json.load(open(f"{RES}/array_calibration_20261105/pilercr_shuffle_wrapper/summary.json"))["all"]
    t1, t2 = w.x.max() + 2.9, w.x.max() + 7.5
    ax.text(t1, 1.5, "$\\hat p\\leq$0.05", fontsize=5.3, ha="center", va="center", color=MUTED, clip_on=False); ax.text(t2, 1.5, "Bonf.(31)", fontsize=5.3, ha="center", va="center", color=MUTED, clip_on=False)
    for yy, k in ((1, "D_conservative_p"), (0, "pilercr_conservative_p")):
        ax.text(t1, yy, f"{s_[k]['uncorrected_le_0_05']}/31", fontsize=6.0, ha="center", va="center", clip_on=False, fontweight="bold")
        ax.text(t2, yy, f"{s_[k]['bonferroni_31_le_0_05']}/31", fontsize=6.0, ha="center", va="center", clip_on=False, fontweight="bold")
    ax.set_title("calibration: same shuffles, same 2×1,000 draws", loc="left", pad=3)
    ax.legend(handles=[Line2D([], [], marker="s", ls="", color=CORE, ms=4, label="Bonferroni(31) $\\leq$ 0.05"), Line2D([], [], marker="s", ls="", mfc="white", mec=CORE, ms=3.8, label="$\\hat p \\leq$ 0.05 only"), Line2D([], [], marker="x", ls="", color=BAD, ms=3.5, label="no array, $\\hat p$ = 1")],
              ncol=3, loc="lower left", bbox_to_anchor=(-0.02, 0.135), handletextpad=0.1, columnspacing=1.0, borderpad=0, fontsize=5.8)
    ax.text(-0.015, 0.06, "Tuned PILER-CR called an array on 0 of 62,000 shuffles, so its wrapped $\\hat p$ only", transform=ax.transAxes, fontsize=5.7, va="center", color=INK)
    ax.text(-0.015, -0.015, "restates whether it called the observed window. No reference window was called.", transform=ax.transAxes, fontsize=5.7, va="center", color=INK)
    save(fig, "fig_baselines")

# ------------------------------------------------------------------ gene arrows shared by locus figures
def gene(ax, s, e, y, h, strand, fc, ec, lw=0.5, z=2):
    s, e = float(s), float(e); head = min(110, (e - s) * 0.55)
    if strand >= 0: pts = [(s, y - h), (e - head, y - h), (e, y), (e - head, y + h), (s, y + h)]
    else: pts = [(e, y - h), (s + head, y - h), (s, y), (s + head, y + h), (e, y + h)]
    ax.add_patch(Polygon(pts, closed=True, fc=fc, ec=ec, lw=lw, zorder=z, joinstyle="round"))

def fig_references():
    from Bio import SeqIO
    samples = pd.read_csv(f"{DATA}/art_samples_2026-11-05/art_samples.csv").set_index("locus_id")
    rr = tsv(f"{RES}/art_validation/reference_recovery.tsv").set_index("locus_id")
    repeats = tsv(f"{RES}/art_arrays_assessed/repeats.tsv")
    order = ["LC680885.1__BDE75684.1", "MN091626.1__QDJ97640.1", "MW218148.1__QPI16926.1", "MW248466.1__QQM14740.1", "MZ422438.1__QXN67953.1",
             "MW349129.1__QQO92858.1", "MZ779063.1__UAJ17040.1", "OR836606.1__WRM43423.1", "MZ422438.1__QXN67955.1"]
    why = {"MW349129.1__QQO92858.1": "repeats overlap CDS", "MZ779063.1__UAJ17040.1": "fails locus correction",
           "OR836606.1__WRM43423.1": "fails correction, intervening gene", "MZ422438.1__QXN67955.1": "too few copies"}
    fig = plt.figure(figsize=(6.5, 3.55))
    gs = fig.add_gridspec(1, 2, width_ratios=[3.25, 1], wspace=0.2, left=0.012, right=0.992, top=0.9, bottom=0.115)
    ax = fig.add_subplot(gs[0, 0]); tag(ax, "a", x=0.0, y=1.075)
    h = 0.17
    for i, loc in enumerate(order):
        rec = SeqIO.read(f"{RES}/oriented_art_loci/loci/{loc}.gb", "genbank"); sm = samples.loc[loc]
        cds = {f.qualifiers.get("protein_id", [""])[0]: f for f in rec.features if f.type == "CDS"}
        rt = cds[sm.RT_feature_id]; o = int(rt.location.start); y = i
        ax.plot([-2300, 4400], [y, y], color=FAINT, lw=0.5, zorder=1)
        for pid, f in cds.items():
            s, e = int(f.location.start) - o, int(f.location.end) - o
            if e < -2300 or s > 4400: continue
            s2, e2 = max(s, -2300), min(e, 4400); st = f.location.strand or 1
            if pid == sm.RT_feature_id:
                b = int(sm.NTD_length_aa) * 3
                ax.add_patch(Rectangle((0, y - h), b, 2 * h, fc=NTD, ec=CORE, lw=0.4, zorder=3))
                gene(ax, b, e, y, h, 1, CORE, CORE, z=3)
            elif pid == sm.partner_feature_id: gene(ax, s2, e2, y, h, st, PARTNER, "#B9481B", z=3)
            else: gene(ax, s2, e2, y, h, st, OTHER, OTHER_E)
        rp = repeats[repeats.locus_id == loc]
        xs = (rp.locus_start0 + rp.locus_end0).values / 2 - o
        if len(xs): ax.plot(xs, [y] * len(xs), ls="", marker="D", ms=2.5, color=INK, mew=0, zorder=5); ax.plot([xs.min(), xs.max()], [y, y], color=INK, lw=0.6, zorder=4)
        r = rr.loc[loc]; ok = r.locus_status.startswith("supported")
        ax.text(-2330, y - 0.31, SHORT[loc], fontsize=6.2, fontweight="bold", va="center", ha="left")
        copies = int(r.detected_copies) if r.detected_copies else 0
        st = "supported" if ok else f"review: {why[loc]}"
        ax.text(4400, y - 0.31, f"{st}  \u00b7  {copies} copies  \u00b7  $\\tilde p$ = {r.adjusted_p:.4f}".replace("0 copies", "no call"), fontsize=5.7, va="center", ha="right", color=INK if ok else "#8A5F0A")
    ax.axhline(4.5, color=MUTED, lw=0.5, ls=(0, (3, 2)))
    ax.set_xlim(-2350, 4420); ax.set_ylim(len(order) - 0.45, -0.62)
    ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xticks(range(-2000, 4001, 1000)); ax.set_xticklabels([f"{x / 1000:g}" for x in range(-2000, 4001, 1000)])
    ax.set_xlabel("position relative to RT start (kb)")
    hs = [Line2D([], [], marker="D", ls="", color=INK, ms=3), Rectangle((0, 0), 1, 1, fc=NTD, ec=CORE, lw=0.4), Rectangle((0, 0), 1, 1, fc=CORE), Rectangle((0, 0), 1, 1, fc=PARTNER), Rectangle((0, 0), 1, 1, fc=OTHER, ec=OTHER_E, lw=0.4)]
    ax.legend(hs, ["repeat copy", "RT NTD", "RT core", "partner", "other CDS"], ncol=5, loc="lower left", bbox_to_anchor=(0.03, 1.0), handlelength=1.1, columnspacing=1.0, borderpad=0)
    # b -- core recovery is not boundary reproduction
    ax = fig.add_subplot(gs[0, 1]); tag(ax, "b", x=-0.14, y=1.075)
    for loc in order[:4] + order[5:8]:
        r = rr.loc[loc]; g, me = json.loads(r.gold_spacers_nt_json), json.loads(r.measured_spacers_nt_json)
        sa1 = loc == "MW218148.1__QPI16926.1"
        ax.scatter(g, np.array(me) - np.array(g), s=15, marker="^" if sa1 else "o", facecolor=PARTNER if sa1 else "none", edgecolor=PARTNER if sa1 else CORE, linewidth=0.6)
    ax.axhline(0, color=FAINT, lw=0.7)
    ax.set_ylim(-0.6, 6.6); ax.set_xlim(143, 198); ax.set_xlabel("spacer, 22-nt core windows (nt)"); ax.set_ylabel("our spacer \u2212 published (nt)", labelpad=1)
    ax.text(170.5, 5.55, "17-nt repeats: +5", fontsize=5.6, ha="center", va="bottom")
    ax.text(170.5, 2.55, "20-nt repeats (SA1): +2", fontsize=5.6, ha="center", va="bottom")
    ax.set_title("boundary offset", loc="left", pad=3)
    save(fig, "fig_references")

# ------------------------------------------------------------------ geometry + logos
LOGO_COL = {"A": "#14A37F", "C": "#2A76D2", "G": "#C98A12", "T": "#CF3A32"}
_FP = FontProperties(family="DejaVu Sans", weight="bold")
def draw_logo(ax, seqs):
    L = max(map(len, seqs)); n = len(seqs)
    for c in range(L):
        col = [s[c] for s in seqs if c < len(s) and s[c] in "ACGT"]
        cnt = Counter(col); N = sum(cnt.values())
        if N < 0.5 * n: continue            # a column most copies do not reach would look falsely conserved
        H = -sum(v / N * math.log2(v / N) for v in cnt.values()); R = 2 - H
        y = 0
        for b, v in sorted(cnt.items(), key=lambda kv: kv[1]):
            ht = v / N * R
            if ht <= 0: continue
            tp = TextPath((0, 0), b, size=1, prop=_FP); ext = tp.get_extents()
            tr = Affine2D().translate(-ext.x0, -ext.y0).scale(0.9 / ext.width, ht / ext.height).translate(c + 0.05, y)
            ax.add_patch(PathPatch(tp, transform=tr + ax.transData, fc=LOGO_COL[b], lw=0)); y += ht
    ax.set_xlim(0, 30.2); ax.set_ylim(0, 2.05); ax.set_yticks([0, 2]); ax.set_yticklabels(["0", "2"], fontsize=5.3); ax.set_xticks([])
    ax.spines["bottom"].set_color(FAINT); ax.tick_params(axis="y", length=1.5, pad=1)

def fig_geometry():
    fig = plt.figure(figsize=(6.5, 2.75))
    gs = fig.add_gridspec(5, 2, width_ratios=[1.15, 1], wspace=0.17, hspace=0.42, left=0.07, right=0.995, top=0.93, bottom=0.135)
    ax = fig.add_subplot(gs[:, 0]); tag(ax, "a", x=-0.09)
    rep = tsv(f"{RES}/art_arrays_assessed/repeats.tsv"); sp = tsv(f"{RES}/art_arrays_assessed/spacers.tsv")
    sup = tsv(f"{RES}/art_type1_catalog/supported_array_candidates.tsv")
    ax.axvspan(20, 50, color=OTHER, lw=0, zorder=0)
    ax.text(53, 43.9, "\u2190 finder defaults, \u224820\u201350 nt", fontsize=5.6, ha="left", va="top", color=MUTED)
    ax.add_patch(Rectangle((25, 15), 395, 30, fc="none", ec=PARTNER, lw=0.7, ls=(0, (3, 2)), zorder=1))
    ax.text(27.5, 15.6, "tuned baseline range", fontsize=5.6, ha="left", va="bottom", color="#B9481B")
    groups = []; rng = np.random.default_rng(3)
    staph = [l for l in rep.locus_id.unique() if not l.startswith("MZ422438")]
    groups.append(("7 Staphylococcus references", "o", CORE, sp[sp.locus_id.isin(staph)].length_nt.values, sp[sp.locus_id.isin(staph)].merge(rep.groupby("locus_id").length_nt.first().rename("w"), on="locus_id").w.values))
    mark = {"ERR8055547": ("ERR8055547_17635", "D"), "MZ422438": ("LPJP1 QXN67953.1", "s"), "MGYC002687634053": ("MGYC002687634053", "^"), "MGYC003526048914": ("MGYC003526048914", "v")}
    all_w, all_s = [], []
    sup_seq = {}
    for r in sup.itertuples():
        d = json.load(open(f"{RES}/{r.array_details.replace('results/', '', 1)}"))
        wd = [q["end0"] - q["start0"] for q in d["repeats"]]; gaps = json.loads(r.spacer_gaps_nt)
        key = next(k for k in mark if r.locus_id.startswith(k)); lab, mk = mark[key]
        groups.append((lab, mk, PARTNER, np.array(gaps), np.array([wd[0]] * len(gaps))))
        sup_seq[lab] = [q["sequence"] for q in d["repeats"]]; all_w += wd; all_s += gaps
    truth = tsv(f"{RES}/art_validation/synthetic_truth.tsv")
    dg = truth[truth.test.isin(["long_spacer_120", "long_spacer_180", "long_spacer_220", "repetitive_DNA"])]
    sx, sy = [], []
    for r in dg.itertuples():
        iv = json.loads(r.repeat_intervals_json)
        for g in set(json.loads(r.known_spacers_json)): sx.append(g); sy.append(iv[0][1] - iv[0][0])
    ax.scatter(sx, sy, s=26, marker="+", color=MUTED, linewidth=0.7, zorder=2, label="designed benchmark cases")
    for lab, mk, c, xs, ys in groups:
        j = rng.uniform(-0.28, 0.28, len(xs))
        ax.scatter(xs, ys + j, s=11 if c == CORE else 16, marker=mk, color=c, edgecolor="white", linewidth=0.3, alpha=0.9, zorder=3, label=lab)
    ax.set_xscale("log"); ax.set_xlim(18, 900); ax.set_ylim(13, 46)
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v:g}")); ax.set_xticks([20, 50, 100, 200, 400, 800])
    ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    ax.set_xlabel("gap between consecutive repeats (nt)"); ax.set_ylabel("repeat width (nt)")
    ax.set_title("array geometry vs. what finders search", loc="left", pad=3)
    ax.legend(loc="upper left", bbox_to_anchor=(0.30, 0.90), fontsize=5.2, handletextpad=0.15, labelspacing=0.2, borderpad=0.1)
    print("supported tracts: repeat width", min(all_w), max(all_w), "gap", min(all_s), max(all_s))
    print("references: repeat width", rep.length_nt.min(), rep.length_nt.max(), "spacer", sp.length_nt.min(), sp.length_nt.max())
    # b -- logos
    staph_all = []
    for l in staph:
        for rec in open(f"{RES}/art_arrays_assessed/alignments/{l}__assessed1.final_repeat_alignment.fna").read().split(">")[1:]:
            staph_all.append("".join(rec.split("\n")[1:]))
    # centre every copy on the published core so widths of 17 and 20 nt stack in register
    core = "ATATGAATACGTAT"; aligned = []
    for q in staph_all:
        q = q.replace("-", "")
        k = min(range(len(q) - len(core) + 1), key=lambda i: sum(a != b for a, b in zip(q[i:], core)))
        aligned.append(" " * (7 - k) + q)
    logos = [("7 Staph. references  (%d copies)" % len(staph_all), aligned)] + [(f"{k}  ({len(v)} copies)", v) for k, v in sup_seq.items()]
    axs = []
    for i, (lab, seqs) in enumerate(logos):
        a = fig.add_subplot(gs[i, 1]); draw_logo(a, seqs); axs.append(a)
        a.text(0.998, 1.02, lab, transform=a.transAxes, fontsize=5.5, ha="right", va="bottom")
        if i == 0:
            tag(a, "b", x=-0.045, y=1.02); a.set_title("repeat logos (bits)", loc="left", pad=2, fontsize=7)
            a.plot([7, 21], [-0.16, -0.16], color=INK, lw=0.7, clip_on=False)
            a.text(6.7, -0.16, "published core", fontsize=4.9, ha="right", va="center", color=MUTED)
    save(fig, "fig_geometry")

# ------------------------------------------------------------------ acquisition flow
def fig_acquisition():
    hp = json.load(open(f"{RES}/art_expansion_20261105/HPGC/summary.json")); jp = json.load(open(f"{RES}/art_expansion_20261105/JPGD/summary.json"))
    st = json.load(open(f"{RES}/art_expansion_20261105/step5/summary.json"))
    fig, ax = plt.subplots(figsize=(6.5, 2.95)); fig.subplots_adjust(0, 0, 1, 1); ax.set_xlim(0, 650); ax.set_ylim(0, 295); ax.axis("off")
    def box(x, y, w, h, big, small, fc="white", ec=FAINT, col=INK, bs=8.5):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=3.5", fc=fc, ec=ec, lw=0.7, zorder=2))
        if small:
            ax.text(x + w / 2, y + h * 0.66, big, ha="center", va="center", fontsize=bs, fontweight="bold", color=col, zorder=3)
            ax.text(x + w / 2, y + h * 0.28, small, ha="center", va="center", fontsize=5.3, color=MUTED, zorder=3, linespacing=1.03)
        elif big:
            ax.text(x + w / 2, y + h / 2, big, ha="center", va="center", fontsize=bs, fontweight="bold", color=col, zorder=3)
    def arrow(x0, y0, x1, y1, c=FAINT, lw=0.8):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=6, color=c, lw=lw, shrinkA=0, shrinkB=0, zorder=1))
    def bus(pts, c=FAINT, lw=0.8):
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color=c, lw=lw, zorder=1, solid_joinstyle="round")
    head = lambda x, s: ax.text(x, 284, s, fontsize=6.3, fontweight="bold", va="center")
    head(4, "bounded sources"); head(114, "after dedup and mapping"); head(261, "new pairs"); head(376, "merged"); head(510, "evidence tiers")
    rows = [("HPGC huge phage", f"{hp['records_in_archive']:,}", "records", f"{hp['unique_new_genomes']:,}", f"genomes · {hp['proteins_searched'] / 1e6:.2f} M proteins", 1, 226),
            ("JPGD jumbo phage", f"{jp['records_in_archive']:,}", "records", f"{jp['unique_new_genomes']:,}", f"genomes · {jp['proteins_searched'] / 1e6:.2f} M proteins", 4, 180),
            ("Archive leads", "45", "from 44 runs", "24", "contigs mapped · 21 unresolved", 16, 134),
            ("UHGV gut virome", "120", "selected genomes", "107", "retrieved · 13 not retrieved", 0, 88)]
    for name, a, al, b, bl, c, y in rows:
        box(4, y, 94, 38, a, f"{name}\n{al}", bs=8)
        box(114, y, 128, 38, b, bl, bs=8)
        box(261, y + 4, 62, 30, str(c), "pair" if c == 1 else "pairs", fc=SHADE if c else "white", ec=CORE if c else FAINT, bs=8)
        arrow(98, y + 19, 114, y + 19); arrow(242, y + 19, 261, y + 19); bus([(323, y + 19), (350, y + 19)])
    box(4, 16, 238, 56, "", "", fc="#F7F7F5")
    ax.text(12, 62, "retained from earlier cohorts", fontsize=5.8, fontweight="bold", va="center")
    for i, (n, l) in enumerate((("13", "frozen type I catalog"), ("2", "earlier follow-up pairs"), ("10", "type II references"))):
        ax.text(14, 49 - i * 12.2, n, fontsize=6.8, fontweight="bold", va="center", ha="left"); ax.text(32, 49 - i * 12.2, l, fontsize=5.8, va="center", color=MUTED)
    bus([(242, 44), (350, 44)]); bus([(350, 44), (350, 245)]); arrow(350, 165, 376, 165)
    box(376, 138, 106, 54, "49", "raw pair records\n(all source maps replayed)", bs=10)
    arrow(429, 138, 429, 119)
    ax.text(424, 128.5, "3 aliases collapsed", fontsize=5.3, color=MUTED, va="center", ha="right")
    box(376, 66, 106, 53, f"{st['analysis_locus_representatives']}", "RT–partner pairs\n(checksum-verified)", fc=SHADE, ec=CORE, bs=11)
    tiers = [(226, st["typeI_pairs"], "type I pairs", f"in {st['typeI_partner70_clusters']} partner clusters"), (180, st["complete_typeI_pairs"], "complete type I", "translation-matched CDS"),
             (134, st["complete_strict_screen_typeI_pairs"], "strict complete type I", "pass the ART RT screen"), (88, st["block_feasibility"]["matching_pair_edges"], "matched-block edge", "among 153 strict pairs"), (42, st["block_feasibility"]["matching_triples"], "eligible triples", "8 disjoint blocks required")]
    for y, n, l, sl in tiers:
        box(510, y, 136, 38, "", "", fc="white" if n else "#FBE9E7", ec=FAINT if n else BAD)
        ax.text(519, y + 19, str(n), fontsize=9.5, fontweight="bold", va="center", color=INK if n else BAD, zorder=3)
        ax.text(543, y + 25, l, fontsize=5.9, va="center", zorder=3); ax.text(543, y + 13, sl, fontsize=5.3, va="center", color=MUTED, zorder=3)
    for a, b in zip(tiers[:-1], tiers[1:]): arrow(578, a[0], 578, b[0] + 38)
    bus([(482, 92), (496, 92), (496, 245)]); arrow(496, 245, 510, 245)
    ax.text(429, 44, "10 type II pairs\nkept as a separate arm", fontsize=5.5, va="center", ha="center", color=MUTED, linespacing=1.05)
    arrow(429, 66, 429, 54)
    save(fig, "fig_acquisition")

# ------------------------------------------------------------------ phylogeny + tanglegram
def fig_trees():
    from Bio import Phylo
    samples = pd.read_csv(f"{DATA}/art_samples_2026-11-05/art_samples.csv")
    step5 = f"{RES}/art_expansion_20261105/step5"
    cat = tsv(f"{step5}/pair_catalog.tsv")
    key = {}
    for r in tsv(f"{RES}/art_evidence_census_20261105/sample_evidence.tsv")[["tree_key", "locus_id"]].itertuples(): key[r.tree_key] = r.locus_id
    sid = samples.set_index("locus_id")
    labels = {}
    for line in open(f"{step5}/RT_core.input.faa"):
        if line.startswith(">"): k, _, v = line[1:].strip().partition(" "); labels[k] = v
    def cls(t):
        if t.startswith("L"):
            l = key.get(t)
            if l in sid.index: return "type I" if sid.loc[l, "partner_type"] == "Type_I" else "type II"
            return "other"
        v = labels.get(t, "")
        if t.startswith("a"): return "ART ref"
        for suf, c in (("Retrons", "retron"), ("DGRs", "DGR"), ("CRISPR", "CRISPR"), ("GII", "group II")):
            if v.endswith(suf): return c
        return "other RT"
    tree = Phylo.read(f"{step5}/RT_core.FastTree.nwk", "newick"); tree.root_at_midpoint(); tree.ladderize()
    fig = plt.figure(figsize=(6.5, 3.75))
    ax = fig.add_axes([0.0, 0.0, 0.5, 0.97], projection="polar"); ax.axis("off")
    tips = tree.get_terminals(); n = len(tips); theta = {t: 2 * np.pi * i / n for i, t in enumerate(tips)}
    depths = tree.depths(); maxd = max(depths[t] for t in tips)
    def ang(c):
        if c.is_terminal(): return theta[c]
        return np.mean([ang(k) for k in c.clades])
    tcol = {"type I": PARTNER, "type II": CORE, "ART ref": INK, "retron": "#7A7F87", "DGR": "#A7ADB5", "CRISPR": "#C6B38E", "group II": "#BFC7D1", "other RT": "#D6D8DC", "other": "#D6D8DC"}
    def edges(c):
        a = ang(c); r = depths[c]
        for k in c.clades:
            ak = ang(k); rk = depths[k]
            ax.plot([ak, ak], [r, rk], color="#9AA0A8", lw=0.42, solid_capstyle="round", zorder=2)
            aa = np.linspace(min(a, ak), max(a, ak), max(4, int(abs(a - ak) * 30)))
            ax.plot(aa, [r] * len(aa), color="#9AA0A8", lw=0.42, zorder=2)
            edges(k)
    edges(tree.root)
    rr = maxd * 1.04; counts = Counter()
    for t in tips:
        c = cls(t.name); counts[c] += 1
        ax.plot([theta[t], theta[t]], [depths[t], rr], color=GRID, lw=0.3, zorder=1)
        hi = c in ("type I", "type II", "ART ref")
        ax.scatter([theta[t]], [rr], s=8 if hi else 5, color=tcol[c], marker={"type I": "o", "type II": "s", "ART ref": "D"}.get(c, "."), linewidth=0, zorder=4)
    ax.set_ylim(0, rr * 1.03)
    print("tree classes", counts)
    h = [Line2D([], [], marker=m, ls="", color=tcol[k], ms=ms, label=f"{k} ({counts[k]})") for k, m, ms in (("type I", "o", 3.2), ("type II", "s", 3), ("ART ref", "D", 2.8), ("retron", ".", 4), ("DGR", ".", 4), ("CRISPR", ".", 4), ("group II", ".", 4), ("other RT", ".", 4))]
    fig.legend(handles=h, loc="lower left", bbox_to_anchor=(0.0, -0.005), ncol=4, fontsize=5.3, handletextpad=0.05, columnspacing=0.55, borderpad=0)
    fig.text(0.005, 0.975, "a", fontweight="bold", fontsize=8.5, va="top"); fig.text(0.03, 0.975, "RT-core tree, 137 sequences", fontsize=7.5, va="top")
    # b -- tanglegram over the 36 type I pairs
    partner = Phylo.read(f"{step5}/TypeI_partner.FastTree.nwk", "newick"); partner.root_at_midpoint(); partner.ladderize()
    keep = {t.name for t in partner.get_terminals()}
    rt = Phylo.read(f"{step5}/RT_core.FastTree.nwk", "newick")
    for t in list(rt.get_terminals()):
        if t.name not in keep: rt.prune(t)
    rt.root_at_midpoint(); rt.ladderize()
    def layout(tr, tiporder=None):
        tips_ = tr.get_terminals(); ypos = {t.name: i for i, t in enumerate(tips_)}
        return ypos
    # barycentric untangling of the partner tree against the RT tip order
    left_y = layout(rt)
    def untangle(c):
        for k in c.clades: untangle(k)
        c.clades.sort(key=lambda k: np.mean([left_y[t.name] for t in k.get_terminals()]))
    def untangle_rt(c, ref):
        for k in c.clades: untangle_rt(k, ref)
        c.clades.sort(key=lambda k: np.mean([ref[t.name] for t in k.get_terminals()]))
    for _ in range(3):
        untangle(partner.root); right_y = layout(partner); untangle_rt(rt.root, right_y); left_y = layout(rt)
    untangle(partner.root); right_y = layout(partner)
    ax = fig.add_axes([0.5, 0.045, 0.495, 0.885]); ax.axis("off")
    gap_l, gap_r = 0.40, 0.60
    def draw(tr, ypos, x0, sign, width):
        d = tr.depths(); md = max(d[t] for t in tr.get_terminals()); sc = width / md
        def yy(c): return ypos[c.name] if c.is_terminal() else np.mean([yy(k) for k in c.clades])
        def rec(c):
            y0 = yy(c); x = x0 + sign * (d[c] * sc)
            for k in c.clades:
                xk = x0 + sign * (d[k] * sc); yk = yy(k)
                ax.plot([x, x], [y0, yk], color=MUTED, lw=0.5, solid_capstyle="round"); ax.plot([x, xk], [yk, yk], color=MUTED, lw=0.5, solid_capstyle="round")
                conf = k.confidence
                if not k.is_terminal() and conf is not None and conf >= 0.9: ax.scatter([xk], [yk], s=3.5, color=INK, zorder=3, linewidth=0)
                rec(k)
        rec(tr.root)
    draw(rt, left_y, 0.0, 1, gap_l - 0.07); draw(partner, right_y, 1.0, -1, 1 - gap_r - 0.07)
    for t in keep:
        l = key.get(t); s = sid.loc[l] if l in sid.index else None
        sup = bool(s is not None and s.category6_historical_supported_array)
        strict = bool(s is not None and s.pair_complete and s.category1_ART_profile_pass)
        ya, yb = left_y[t], right_y[t]
        xs = np.linspace(gap_l + 0.015, gap_r - 0.015, 30); u = (xs - xs[0]) / (xs[-1] - xs[0]); ys = ya + (yb - ya) * (3 * u**2 - 2 * u**3)
        ax.plot(xs, ys, color=PARTNER if sup else ("#7DA9DF" if strict else "#D5D8DD"), lw=1.1 if sup else 0.6, zorder=3 if sup else (2 if strict else 1))
        nm = s.sample_id.replace("ARTS00", "") if s is not None else t
        ax.text(gap_l - 0.062, ya, nm, fontsize=4.9, ha="left", va="center", color=INK, fontweight="bold" if sup else "normal")
        ax.text(gap_r + 0.062, yb, nm, fontsize=4.9, ha="right", va="center", color=INK, fontweight="bold" if sup else "normal")
    ax.set_ylim(len(keep) - 0.4, -0.6); ax.set_xlim(-0.01, 1.01)
    fig.text(0.505, 0.975, "b", fontweight="bold", fontsize=8.5, va="top"); fig.text(0.53, 0.975, "36 type I pairs: RT core (left) vs. partner (right)", fontsize=7.5, va="top")
    fig.legend(handles=[Line2D([], [], color=PARTNER, lw=1.1, label="supported array"), Line2D([], [], color="#7DA9DF", lw=0.8, label="strict complete"), Line2D([], [], color="#D5D8DD", lw=0.8, label="other"), Line2D([], [], marker=".", ls="", color=INK, ms=3.2, label="SH-like support \u2265 0.9")],
               loc="lower right", bbox_to_anchor=(0.995, -0.005), ncol=4, fontsize=5.3, handletextpad=0.25, columnspacing=0.8, handlelength=1.2, borderpad=0)
    save(fig, "fig_trees")

# ------------------------------------------------------------------ dataset landscape
def fig_landscape():
    s = pd.read_csv(f"{DATA}/art_samples_2026-11-05/art_samples.csv")
    def st(r):
        a = r.array_status
        if r.category6_historical_supported_array: return "supported"
        if a == "statistically_supported_requires_review": return "supported, held for review"
        if a == "incomplete_upstream_assessment": return "truncated upstream"
        if a in ("not_supported", "not_supported_after_locus_correction", "no_supported_array_exploratory", "no_supported_array_call_exploratory_scan"): return "assessed, no support"
        return "not assessed or calibrated"
    s["cat"] = s.apply(st, axis=1); s["t"] = s.partner_type.map({"Type_I": 0, "Type_II": 1})
    s = s.sort_values(["t", "available_upstream_window_nt", "RT_length_aa"], ascending=[True, False, False]).reset_index(drop=True)
    x = np.arange(len(s), dtype=float); x[s.t == 1] += 1.2
    cols = {"supported": CORE, "supported, held for review": "#7DA9DF", "assessed, no support": MUTED, "not assessed or calibrated": FAINT, "truncated upstream": PARTNER}
    fig = plt.figure(figsize=(6.5, 3.15))
    gs = fig.add_gridspec(2, 2, width_ratios=[3.3, 1], height_ratios=[1, 1.5], hspace=0.14, wspace=0.05, left=0.07, right=0.995, top=0.93, bottom=0.085)
    ax = fig.add_subplot(gs[0, 0]); tag(ax, "a", x=-0.03)
    ax.bar(x, s.available_upstream_window_nt / 1000, width=0.72, color=[cols[c] for c in s.cat], linewidth=0)
    ax.axhline(6, color=INK, lw=0.6, ls=(0, (3, 2))); ax.set_ylim(0, 7.3); ax.set_xlim(-0.8, x.max() + 0.8); ax.set_yticks([0, 2, 4, 6])
    ax.set_ylabel("upstream DNA (kb)"); ax.tick_params(axis="x", length=0, labelbottom=False); ax.spines["bottom"].set_color(FAINT)
    for arr, l in ((x[s.t == 0], "type I (36)"), (x[s.t == 1], "type II (10)")):
        ax.text(arr.mean(), 6.75, l, ha="center", fontsize=6.2, fontweight="bold", va="center")
    ax.set_title("how much of the 6-kb window each locus has, coloured by array status", loc="left", pad=3)
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax); tag(ax2, "b", x=-0.03, y=0.98)
    ntd = s.NTD_length_aa.values; rt = s.RT_length_aa.values; pa = s.partner_length_aa.values
    ax2.bar(x, ntd, width=0.72, color=NTD, linewidth=0); ax2.bar(x, rt - ntd, bottom=ntd, width=0.72, color=CORE, linewidth=0)
    ax2.bar(x, -pa, width=0.72, color=PARTNER, linewidth=0); ax2.axhline(0, color=INK, lw=0.5)
    inc = ~s.pair_complete.values
    ax2.scatter(x[inc], (rt[inc] + 34), marker="v", s=6, color=INK, linewidth=0, zorder=4)
    ax2.set_ylim(-1120, 790); ax2.set_yticks([-1000, -500, 0, 500]); ax2.set_yticklabels(["1,000", "500", "0", "500"]); ax2.set_ylabel("partner  \u00b7  RT (aa)", labelpad=1)
    ax2.set_xticks(x); ax2.set_xticklabels([v.replace("ARTS00", "") for v in s.sample_id], fontsize=4.3, rotation=90); ax2.tick_params(axis="x", length=0, pad=1)
    ax2.spines["bottom"].set_visible(False)
    ax2.legend(handles=[Rectangle((0, 0), 1, 1, fc=NTD), Rectangle((0, 0), 1, 1, fc=CORE), Rectangle((0, 0), 1, 1, fc=PARTNER), Line2D([], [], marker="v", ls="", color=INK, ms=2.8)], labels=["RT NTD", "RT core", "partner", "incomplete CDS"], ncol=4, loc="lower left", bbox_to_anchor=(0.2, -0.01), handlelength=1, columnspacing=0.9, borderpad=0, fontsize=5.6)
    ax3 = fig.add_subplot(gs[:, 1]); ax3.axis("off"); ax3.set_xlim(0, 1); ax3.set_ylim(0, 1)
    ax3.text(0.08, 0.985, "array status", fontsize=6.8, fontweight="bold", va="top")
    y = 0.92
    for cat in ["supported", "supported, held for review", "assessed, no support", "not assessed or calibrated", "truncated upstream"]:
        c1, c2 = ((s.cat == cat) & (s.t == 0)).sum(), ((s.cat == cat) & (s.t == 1)).sum()
        ax3.add_patch(Rectangle((0.08, y - 0.035), 0.07, 0.035, fc=cols[cat], lw=0))
        ax3.text(0.19, y - 0.018, cat, fontsize=5.7, va="center")
        ax3.text(0.19, y - 0.062, f"type I: {c1}   type II: {c2}", fontsize=5.5, va="center", color=MUTED)
        y -= 0.118
    y -= 0.03
    ax3.text(0.08, y, f"{(s[s.t == 0].available_upstream_window_nt < 6000).sum()} of {(s.t == 0).sum()} type I loci", fontsize=6.3, fontweight="bold", va="top")
    ax3.text(0.08, y - 0.055, "lack the full window. Each is\nrecorded as an absence of\nevidence, not a negative.", fontsize=5.6, va="top", color=MUTED, linespacing=1.15)
    y -= 0.235
    tI = s[s.t == 0]; tII = s[s.t == 1]
    ax3.text(0.08, y, "median NTD", fontsize=6.3, fontweight="bold", va="top")
    ax3.text(0.08, y - 0.055, f"type I {tI.NTD_length_aa.median():.1f} aa (range {tI.NTD_length_aa.min()}\u2013{tI.NTD_length_aa.max()})\ntype II {tII.NTD_length_aa.median():.0f} aa (range {tII.NTD_length_aa.min()}\u2013{tII.NTD_length_aa.max()})", fontsize=5.6, va="top", color=MUTED, linespacing=1.15)
    print(s.cat.value_counts(), s.groupby(['t','cat']).size())
    save(fig, "fig_landscape")

# ------------------------------------------------------------------ binding (exploratory)
def load_binding():
    B = f"{DATA}/results_bindings"
    exp = pd.read_csv(f"{B}/experiments.csv", usecols=["case_id", "stage", "family", "condition", "RT_sample_id", "partner_sample_id", "native_pair", "request_hash", "chain_metadata_json"])
    rows = []
    for p in glob.glob(f"{B}/runs/*/normalized.json"):
        h = p.split("/")[-2]; n = json.load(open(p))
        if not n.get("verified"): continue
        meta = exp[exp.request_hash == h].iloc[0]; ch = json.loads(meta.chain_metadata_json); ntd_end = ch[0].get("source_NTD_end0"); a_start = ch[0].get("start0", 0)
        for m in n["models"]:
            r = dict(request_hash=h, model_index=m["model_index"], plddt_A=m["chain_summary"][0]["mean_CA_plddt"], plddt_B=m["chain_summary"][1]["mean_CA_plddt"] if len(m["chain_summary"]) > 1 else np.nan,
                     iptm=m["scores"]["protein_iptm_scores"], ptm=m["scores"]["ptm_scores"], ntd_end=ntd_end, len_A=m["chain_summary"][0]["length"])
            cf = f"{B}/runs/{h}/model_{m['model_index']}.contacts.csv"
            if os.path.exists(cf):
                c = pd.read_csv(cf); c = c[c.min_heavy_atom_distance_A <= 5.0]
                r["n5"] = len(c); r["ntd5"] = (c.chain0_residue0 < ntd_end).sum() if ntd_end is not None else np.nan
            rows.append(r)
    d = pd.DataFrame(rows).merge(exp.drop(columns=["chain_metadata_json"]), on="request_hash")
    return d, exp

def fig_binding():
    B = f"{DATA}/results_bindings"
    d, exp = load_binding(); planned = len(exp); done = d.request_hash.nunique()
    print("binding: models", len(d), "requests", done, "of", planned)
    fig = plt.figure(figsize=(6.5, 3.3))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.45, 1, 1], height_ratios=[1, 1], wspace=0.42, hspace=0.62, left=0.065, right=0.995, top=0.89, bottom=0.11)
    # a -- the partial matrix, missing cells stay missing
    ax = fig.add_subplot(gs[:, 0]); tag(ax, "a", x=-0.06, y=1.012)
    s = d[(d.stage == "full_length_Type_I") & (d.condition == "S_query_only")]
    ids = sorted(exp[exp.stage == "full_length_Type_I"].RT_sample_id.unique())
    M = s.groupby(["RT_sample_id", "partner_sample_id"]).iptm.median().unstack().reindex(index=ids, columns=ids)
    cmap = LinearSegmentedColormap.from_list("seq", ["#F1F6FD", "#8AB8F0", "#2A76D2", "#123E7C"])
    ax.add_patch(Rectangle((-0.5, -0.5), len(ids), len(ids), facecolor="#F5F5F3", hatch="/////", edgecolor="#DADCDF", lw=0, zorder=0))
    im = ax.imshow(np.ma.masked_invalid(M.values.astype(float)), cmap=cmap, vmin=0.0, vmax=0.6, zorder=2, interpolation="nearest")
    for i in range(len(ids)):
        ax.add_patch(Rectangle((i - 0.5, i - 0.5), 1, 1, fc="none", ec=INK, lw=0.6, zorder=3))
    lab = [v.replace("ARTS00", "") for v in ids]
    ax.set_xticks(range(len(ids))); ax.set_xticklabels(lab, fontsize=4.8, rotation=90); ax.set_yticks(range(len(ids))); ax.set_yticklabels(lab, fontsize=4.8)
    ax.tick_params(length=0, pad=1); ax.set_xlabel("partner", labelpad=1); ax.set_ylabel("RT", labelpad=1)
    for sp_ in ax.spines.values(): sp_.set_visible(False)
    filled = int(M.notna().sum().sum())
    ax.set_title(f"median ipTM: {filled} of {len(ids)**2} cells returned", loc="left", pad=3)
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.04, pad=0.105, aspect=34, ticks=[0, 0.2, 0.4, 0.6]); cb.ax.tick_params(labelsize=5.5, length=1.5, pad=1); cb.outline.set_linewidth(0.3)
    # b -- MSA input raises per-residue confidence, not interface confidence
    ax = fig.add_subplot(gs[0, 1]); tag(ax, "b", x=-0.14)
    pair = d[d.RT_sample_id.isin(["ARTS0001", "ARTS0014"]) & (d.RT_sample_id == d.partner_sample_id) & d.stage.str.startswith("full_length")]
    med = pair.groupby(["RT_sample_id", "condition"])[["plddt_A", "plddt_B", "iptm"]].median().reset_index()
    for sid_, mk in (("ARTS0001", "o"), ("ARTS0014", "s")):
        g = med[med.RT_sample_id == sid_].set_index("condition")
        xs = [0, 1]; ys = [(g.loc["S_query_only", "plddt_A"] + g.loc["S_query_only", "plddt_B"]) / 2, (g.loc["U_requested_independent_MSA_order", "plddt_A"] + g.loc["U_requested_independent_MSA_order", "plddt_B"]) / 2]
        ax.plot(xs, ys, color=MUTED, marker=mk, ms=3.4, lw=0.9, label=sid_)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["query\nonly", "with\nMSA"], fontsize=5.8); ax.set_xlim(-0.3, 1.3); ax.set_ylim(0, 100); ax.set_yticks([0, 50, 100]); ax.set_ylabel("mean C$\\alpha$ pLDDT")
    ax.set_title("per-residue confidence", loc="left", pad=3); ax.legend(loc="lower right", fontsize=5.2, handlelength=1.3, borderpad=0.1, labelspacing=0.2)
    ax = fig.add_subplot(gs[0, 2]); tag(ax, "c", x=-0.14)
    for sid_, mk in (("ARTS0001", "o"), ("ARTS0014", "s")):
        g = med[med.RT_sample_id == sid_].set_index("condition")
        ax.plot([0, 1], [g.loc["S_query_only", "iptm"], g.loc["U_requested_independent_MSA_order", "iptm"]], color=MUTED, marker=mk, ms=3.4, lw=0.9)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["query\nonly", "with\nMSA"], fontsize=5.8); ax.set_xlim(-0.3, 1.3); ax.set_ylim(0, 1); ax.set_yticks([0, 0.5, 1]); ax.set_ylabel("ipTM", labelpad=1)
    ax.axhline(0.6, color=FAINT, lw=0.7, ls=(0, (3, 2))); ax.text(1.28, 0.62, "below 0.6: usually read\nas a failed prediction", fontsize=5.2, color=MUTED, ha="right", va="bottom", linespacing=1.0)
    ax.set_title("interface confidence", loc="left", pad=3)
    # d -- every native pair: where do predicted contacts land?
    ax = fig.add_subplot(gs[1, 1]); tag(ax, "d", x=-0.14)
    nat = s[s.native_pair == True].groupby("RT_sample_id").agg(n5=("n5", "median"), ntd5=("ntd5", "median"), ntd=("ntd_end", "first"), L=("len_A", "first"), iptm=("iptm", "median"))
    nat["share_res"] = nat.ntd / nat.L; nat["share_con"] = nat.ntd5 / nat.n5
    ax.plot([0, 1], [0, 1], color=FAINT, lw=0.7)
    ax.scatter(nat.share_res, nat.share_con, s=11, color=NTD, edgecolor=CORE, linewidth=0.4, zorder=3)
    ax.set_xlim(0.2, 0.6); ax.set_ylim(0.0, 0.8); ax.set_xticks([0.2, 0.4, 0.6]); ax.set_yticks([0, 0.4, 0.8])
    ax.set_xlabel("NTD share of RT residues", labelpad=1); ax.set_ylabel("NTD share of contacts", labelpad=1); ax.set_title("18 native pairs", loc="left", pad=3)
    # e -- the highest-confidence MSA-informed model, C-alpha trace
    ax = fig.add_subplot(gs[1, 2]); tag(ax, "e", x=-0.02)
    u = d[(d.RT_sample_id == "ARTS0001") & (d.partner_sample_id == "ARTS0001") & (d.condition.str.startswith("U")) & (d.stage == "full_length_Type_I")].sort_values("iptm")
    row = u.iloc[len(u) // 2]
    z = np.load(f"{B}/runs/{row.request_hash}/model_{row.model_index}.coordinates.npz")
    A, Bc = z["chain0_CA"], z["chain1_CA"]; X = np.vstack([A, Bc]); X = X - X.mean(0)
    U_, S_, Vt = np.linalg.svd(X, full_matrices=False); P = X @ Vt.T
    from matplotlib.collections import LineCollection
    def segs(P_): return np.stack([P_[:-1, :2], P_[1:, :2]], axis=1)
    pa, pb = P[:len(A)], P[len(A):]; ne = int(row.ntd_end)
    smooth = lambda q: np.vstack([q[0], (q[:-2] + q[1:-1] + q[2:]) / 3, q[-1]])       # light 3-residue smoothing, display only
    pa, pb = smooth(pa), smooth(pb)
    sg = np.vstack([segs(pa), segs(pb)]); dep = np.concatenate([(pa[:-1, 2] + pa[1:, 2]) / 2, (pb[:-1, 2] + pb[1:, 2]) / 2])
    colr = np.array([mpl.colors.to_rgba(NTD)] * ne + [mpl.colors.to_rgba(CORE)] * (len(pa) - 1 - ne) + [mpl.colors.to_rgba(PARTNER)] * (len(pb) - 1))
    o = np.argsort(dep); w = 0.55 + 0.75 * (dep - dep.min()) / (dep.max() - dep.min())
    ax.add_collection(LineCollection(sg[o], colors=colr[o], linewidths=w[o], capstyle="round"))
    ax.set_aspect("equal"); ax.autoscale_view(); ax.axis("off")
    ax.set_title("ARTS0001, with MSA", loc="left", pad=3, x=0.0)
    ax.text(0.5, -0.02, f"ipTM {row.iptm:.2f}  \u00b7  pLDDT {(row.plddt_A + row.plddt_B) / 2:.0f}", transform=ax.transAxes, ha="center", va="top", fontsize=5.4, color=MUTED)
    for yy, t_, c_ in ((1.0, "NTD", "#5E97E0"), (0.93, "core", CORE), (0.86, "partner", "#B9481B")):
        ax.text(-0.13, yy, t_, color=c_, fontsize=5.3, transform=ax.transAxes, ha="left", va="top")
    fig.text(0.008, 0.995, f"Exploratory and partial: {done:,} of {planned:,} requests returned; no contrast computed", fontsize=6.2, ha="left", va="top", color=BAD, fontweight="bold")
    save(fig, "fig_binding")
    print(nat.round(3).to_string()); print(med.round(3).to_string())

# ------------------------------------------------------------------ method overview
def fig_method():
    fig, ax = plt.subplots(figsize=(6.5, 1.85)); fig.subplots_adjust(0, 0, 1, 1); ax.set_xlim(0, 650); ax.set_ylim(0, 185); ax.axis("off")
    def rbox(x, y, w, h, fc="white", ec=FAINT, lw=0.7, z=2, r=4):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw, zorder=z))
    def arrow(a, b, c=MUTED, lw=0.8):
        ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=6, color=c, lw=lw, shrinkA=0, shrinkB=0, zorder=3))
    def stripe(x, y, w, reps=(), c=OTHER_E, h=7, shuffled=None):
        ax.add_patch(Rectangle((x, y), w, h, fc="#F3F3F0", ec=c, lw=0.4, zorder=3))
        for f in reps: ax.add_patch(Rectangle((x + f * w, y), 3.4, h, fc=INK, lw=0, zorder=4))
        if shuffled is not None:
            for f in shuffled: ax.add_patch(Rectangle((x + f * w, y + 1.1), 1.6, h - 2.2, fc=c, alpha=0.55, lw=0, zorder=4))
    # inputs
    ax.text(4, 175, "observed", fontsize=6.0, fontweight="bold", va="center")
    stripe(4, 153, 94, reps=(0.40, 0.55, 0.69, 0.83))
    ax.text(4, 144, "locus $L$, 6 kb upstream of RT", fontsize=5.2, va="center", color=MUTED)
    ax.text(4, 121, "null: $B$ shuffles each", fontsize=6.0, fontweight="bold", va="center")
    rng = np.random.default_rng(1)
    for lab, c, yy in (("$\\mathcal{M}_{\\mathrm{dinuc}}$ Euler-trail shuffle", DINUC, 94), ("$\\mathcal{M}_{\\mathrm{block}}$ 50-nt block shuffle", BLOCK, 48)):
        for k in range(3): stripe(4 + k * 2.6, yy + 13 - k * 4.8, 88, c=c, shuffled=rng.uniform(0.03, 0.96, 12))
        ax.text(4, yy - 3.5, lab, fontsize=5.2, va="center", color=c)
    ax.text(4, 20, "two backgrounds, run\nindependently, never pooled", fontsize=4.9, va="center", color=MUTED, style="italic", linespacing=1.15)
    # D
    rbox(112, 5, 138, 175, fc="#F2F6FC", ec=CORE, lw=0.9)
    ax.text(181, 171, "the same deterministic $D$", fontsize=6.0, fontweight="bold", ha="center", va="center", color="#174E93")
    steps = [("1", "seed enumeration", "$k$ = 10; no filter, no cap"), ("2", "tolerant recovery", "$\\leq$ 1 mismatch"), ("3", "chaining, one spacing", "DP on intervals; $\\leq$ 1 skip"),
             ("4", "boundaries, information", "$\\pm$80 nt; 80% conserved"), ("5", "score and argmax", "$S=(R-1)\\cdot I$; total order")]
    for i, (n, a, b) in enumerate(steps):
        y = 149 - i * 28.5
        rbox(118, y - 12.5, 126, 25, fc="white", ec="#BFD5F2", lw=0.5, r=3)
        ax.text(126, y, n, fontsize=6.6, fontweight="bold", color=CORE, ha="center", va="center", zorder=4)
        ax.text(133, y + 5.2, a, fontsize=5.4, va="center", zorder=4, fontweight="bold")
        ax.text(133, y - 5.2, b, fontsize=5.0, va="center", color=MUTED, zorder=4)
    ax.text(181, 15, "re-executed in every replicate", fontsize=5.3, ha="center", va="center", style="italic", color="#174E93")
    arrow((98, 157), (112, 157), c=INK); arrow((98, 104), (112, 104), c=DINUC); arrow((98, 58), (112, 58), c=BLOCK)
    # real null, drawn on inset axes
    n = null_scores(f"{RES}/art_arrays_assessed/nulls/MZ422438.1__QXN67953.1.json")
    obs = tsv(f"{RES}/art_arrays_assessed/array_assessments.tsv").set_index("locus_id").loc["MZ422438.1__QXN67953.1", "seed_chain_score"]
    ix, iy, iw, ih = 264, 50, 104, 95
    ins = fig.add_axes([ix / 650, iy / 185, iw / 650, ih / 185]); ins.patch.set_alpha(0)
    bins = np.linspace(0, 190, 58)
    ins.hist(n["dinucleotide_6000"], bins=bins, histtype="step", color=DINUC, lw=0.8)
    ins.hist(n["local_blocks_6000"], bins=bins, histtype="step", color=BLOCK, lw=0.8, ls=(0, (2.5, 1.2)))
    ins.axvline(obs, color=INK, lw=1.0, ymax=0.92); ins.set_yticks([]); ins.set_xticks([0, 100, 200]); ins.set_xlim(0, 200); ins.spines["left"].set_visible(False)
    ins.tick_params(labelsize=5, length=1.5, pad=1); ins.set_xlabel("best-chain score", fontsize=5.2, labelpad=0)
    ax.text(ix + obs / 200 * iw, iy + ih + 1, "$S(D(L))$", fontsize=5.4, ha="center", va="bottom")
    ax.text(ix + iw / 2, 175, "search-outcome null", fontsize=6.0, fontweight="bold", ha="center", va="center")
    ax.text(ix + iw / 2, 21, f"real draws, LPJP1 QXN67953.1\n{len(n['dinucleotide_6000']):,} per background", fontsize=4.9, ha="center", va="center", color=MUTED, linespacing=1.12)
    arrow((250, 92), (262, 92))
    # p chain
    arrow((370, 145), (383, 145))
    rbox(383, 110, 121, 70, ec=FAINT)
    ax.text(443.5, 164, "$\\hat p_{\\mathcal{M}}=\\dfrac{1+|\\{b:S(D(X_b))\\geq S(D(L))\\}|}{1+B}$", fontsize=4.8, ha="center", va="center", zorder=4)
    ax.text(446, 134, "$\\hat p=\\max(\\hat p_{\\mathrm{dinuc}},\\hat p_{\\mathrm{block}})$", fontsize=6.0, ha="center", va="center", zorder=4)
    ax.text(446, 120, "a tract must survive both", fontsize=5.0, ha="center", va="center", color=MUTED, zorder=4)
    arrow((446, 110), (446, 97))
    rbox(383, 54, 121, 43, ec=FAINT)
    ax.text(446, 84, "$\\tilde p=\\min(1,\\,m\\,\\hat p)$", fontsize=6.3, ha="center", va="center", zorder=4)
    ax.text(446, 67, "declared family of $m$ loci\nrequires $m/(B+1)\\leq 0.05$", fontsize=5.0, ha="center", va="center", color=MUTED, zorder=4, linespacing=1.15)
    ax.text(446, 28, "relative to the two shuffle\nmodels only: not function,\nnot a census", fontsize=5.0, ha="center", va="center", color=MUTED, linespacing=1.15, style="italic")
    # gate
    arrow((504, 76), (514, 76))
    rbox(514, 5, 132, 175, fc="#FFF8F4", ec=PARTNER, lw=0.9)
    ax.text(580, 171, "non-circular evidence gate", fontsize=5.7, fontweight="bold", ha="center", va="center", color="#B9481B")
    def field(y, t, c, ec):
        rbox(520, y - 7.5, 120, 15, fc=c, ec=ec, lw=0.5, r=2.5); ax.text(580, y, t, fontsize=5.3, ha="center", va="center", zorder=4)
    for i, t in enumerate(("profile evidence", "partner gene context", "RT-only placement")): field(153 - i * 17, t, "white", "#CFD3D8")
    ax.text(580, 106, "may license a family call", fontsize=5.2, ha="center", va="center", fontweight="bold")
    field(87, "array: $\\tilde p$, copies, spacers", "white", INK)
    ax.text(580, 72, "may license an array call", fontsize=5.2, ha="center", va="center", fontweight="bold")
    field(53, "protein completeness", "#F4F4F2", "#CFD3D8"); field(36, "upstream assessability", "#F4F4F2", "#CFD3D8")
    ax.text(580, 19, "can only withhold a call", fontsize=5.2, ha="center", va="center", fontweight="bold")
    save(fig, "fig_method")

ALL = ["method", "calibration", "benchmark", "baselines", "references", "geometry", "acquisition", "trees", "landscape", "binding"]
if __name__ == "__main__":
    for w in (sys.argv[1:] or ALL): globals()[f"fig_{w}"]()
