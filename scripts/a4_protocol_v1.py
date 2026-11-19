"""Protocol v1's criteria and exploratory comparisons, from the saved picks and curves (U-44).

Usage: python scripts/a4_protocol_v1.py [--tables] [--seeds 3-9] [--pbmc-seeds 0-4]
Reads results/v1/{picks,curves,picks_chooser,chooser_r}/ and results/v1/regression.csv, and
writes results/a4_protocol_v1_numbers.json. With --tables it also writes the five
latex/tables/a4_*.tex, latex/figures/a4_seed_means.pdf and latex/numbers_a4.tex. The criteria and their thresholds are fixed in
docs/protocol_matched_v1.md. Summary of a design: the per-dataset median ARI over the seeds,
then the mean over the nine small datasets, as in Assignment 3 (scripts/a3_report.py, agg).

Scorers, as named in the protocol: S0 ours at B=25 (the original protocol), S1 ours at B=100,
S2 ours at B=100 on all cells (PBMC14K), S3-V, S3-G and S3-G2 chooseR's own code at B=100.
"""
import json
import sys
from glob import glob
from pathlib import Path

import numpy as np
import pandas as pd

from scminer_clust import stability
from scminer_clust.paths import RESULTS

V1 = RESULTS / "v1"
TABLES = RESULTS.parent / "latex" / "tables"
FIGS = RESULTS.parent / "latex" / "figures"
NAMES = {"yan": "Yan", "goolam": "Goolam", "buettner": "Buettner", "pollen": "Pollen", "chung": "Chung",
         "usoskin": "Usoskin", "kolodziejczyk": "Kolodziejczyk", "klein": "Klein", "zeisel": "Zeisel",
         "pbmc14k": "PBMC14K"}
TYPES = {"yan": 8, "goolam": 5, "buettner": 3, "pollen": 11, "chung": 4, "usoskin": 4, "kolodziejczyk": 3,
         "klein": 4, "zeisel": 7, "pbmc14k": 7}
MNAME = {"mica": "MICA", "leiden_pca": "Leiden", "kmeans_pca": "$k$-means"}
SMALL = ["yan", "goolam", "buettner", "pollen", "chung", "usoskin", "kolodziejczyk", "klein", "zeisel"]
FIVE = ["yan", "goolam", "buettner", "pollen", "chung"]  # at most 317 cells
METHODS = ["mica", "leiden_pca", "kmeans_pca"]
F, M = "stability_finest", "stability_most"


def short(ds):
    return ds.replace("paper_", "")


def scorer_name(r):
    if r.scorer.startswith("chooseR"):
        return "S3-" + r.scorer.split("_")[1]
    if r.scorer == "ours_all":
        return "S2"
    return "S0" if int(r.b) == 25 else "S1"


def load(sub):
    files = sorted(glob(str(V1 / sub / "*.csv")))
    files = [f for f in files if "crosscheck" not in f]
    if not files:
        return pd.DataFrame()
    d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    d["ds"] = d.dataset.map(short)
    d["S"] = [scorer_name(r) for r in d.itertuples()]
    return d


def summary(p, S, design, reading, method, seeds, datasets=SMALL):
    s = p[(p.S == S) & (p.design == design) & (p.reading == reading) & (p.method == method)
          & p.seed.isin(seeds) & p.ds.isin(datasets)]
    if not len(s):
        return None
    per_ds = s.groupby("ds").ari.median()
    if per_ds.size < len(datasets):  # a mean over fewer datasets is not comparable
        return None
    per_seed = s.groupby("seed").agg(mean=("ari", "mean"), n=("ari", "size"))
    per_seed = per_seed[per_seed.n == len(datasets)]["mean"]
    return {"mean": round(float(per_ds.mean()), 4), "n_datasets": int(per_ds.size),
            "n_runs": int(len(s)), "seeds": sorted(int(x) for x in s.seed.unique()),
            "per_seed": {int(k): round(float(v), 4) for k, v in per_seed.items()},
            "per_dataset": {k: round(float(v), 4) for k, v in per_ds.items()},
            "k_median": {k: float(v) for k, v in s.groupby("ds").n_clusters.median().items()}}


def lead(p, S, design, reading, a, b, seeds):
    x, y = summary(p, S, design, reading, a, seeds), summary(p, S, design, reading, b, seeds)
    if not x or not y:
        return None
    common = sorted(set(x["per_seed"]) & set(y["per_seed"]))
    leads = sum(x["per_seed"][s] > y["per_seed"][s] for s in common)
    return {"a": a, "b": b, "mean_a": x["mean"], "mean_b": y["mean"], "gap": round(x["mean"] - y["mean"], 4),
            "a_leads_seeds": f"{leads} of {len(common)}", "leads": leads, "n_seeds": len(common)}


def verdict_c1(r, need):
    if r is None or r["n_seeds"] < need:
        return "not yet computable"
    if r["gap"] > 0 and r["leads"] >= 5:
        return "holds"
    return "weak" if r["gap"] > 0 else "fails"


def reading_change(p, S, design, method, seeds):
    s = p[(p.S == S) & (p.design == design) & (p.method == method) & p.seed.isin(seeds) & p.ds.isin(SMALL)]
    med = s.pivot_table(index="ds", columns="reading", values="ari", aggfunc="median")
    if F not in med or M not in med:
        return None
    return round(float((med[F] - med[M]).abs().mean()), 4)


def flat_top(c, seeds, design="mid", b=25):
    """Share of settings within 0.01 of the best, pooling each dataset's curves over the seeds,
    as in scripts/a3_matched.py (mechanism)."""
    s = c[(c.scorer == "ours") & (c.b == b) & c[f"in_{design}"] & c.seed.isin(seeds) & c.ds.isin(SMALL)]
    rows = []
    for (ds, m), g in s.groupby(["ds", "method"]):
        t = g[g.stability >= g.stability.max() - 0.01]
        rows.append({"ds": ds, "method": m, "share_tied": len(t) / len(g),
                     "span": int(t.n_clusters.max() - t.n_clusters.min()),
                     "share_at_one": float((g.stability >= 0.9999).mean())})
    K = pd.DataFrame(rows)
    if not len(K):
        return {}
    K["size"] = np.where(K.ds.isin(FIVE), "five_smallest", "four_larger")
    out = {f"{sz}/{m}": {"share_tied": round(float(g.share_tied.mean()), 3), "span": round(float(g.span.mean()), 2),
                         "share_at_one": round(float(g.share_at_one.mean()), 3)} for (sz, m), g in K.groupby(["size", "method"])}
    out["all"] = {m: {"share_tied": round(float(g.share_tied.mean()), 3), "span": round(float(g.span.mean()), 2)}
                  for m, g in K.groupby("method")}
    return out


def exploratory_values():
    """Each criterion's value in the exploratory runs, as their own report recorded it
    (results/a3_matched_numbers.json, from scripts/a3_matched.py). The protocol quotes these."""
    f = RESULTS / "a3_matched_numbers.json"
    if not f.exists():
        return {}
    d = json.loads(f.read_text())
    mm, s, mech = d["matched"], d["summary"], d["matched"]["mechanism_by_size"]
    means = mm["means_nine_small"]
    return {"F.mica": means[f"{F}/mica"], "F.leiden": means[f"{F}/leiden_pca"],
            "M.mica": means[f"{M}/mica"], "M.leiden": means[f"{M}/leiden_pca"],
            "C2.mica": mm["mean_abs_reading_change"]["mica"], "C2.leiden": mm["mean_abs_reading_change"]["leiden_pca"],
            "C3.five.mica": mech["at most 317 cells/mica"]["share_tied"],
            "C3.four.mica": mech["622 cells or more/mica"]["share_tied"],
            "C3.five.leiden": mech["at most 317 cells/leiden_pca"]["share_tied"],
            "C5.leiden": s["Most stable, original grids"]["leiden_pca"], "C5.mica": s["Most stable, original grids"]["mica"]}


def label_free_signals(p, c, seeds, design="mid", b=25, gap=5):
    """Warning signs a user can read without labels, run by run, beside the ARI the labels later give.

    A run is one small dataset, seed and method, scored as the original protocol (S0) on the given
    candidates. The signs need no labels: the share of settings within 0.01 of the run's best
    stability, the span of cluster counts among them, the count each reading picks, and whether a
    pick is the smallest count offered. The ARI is reported only to show what the signs foretold.
    """
    cs = c[(c.scorer == "ours") & (c.b == b) & c[f"in_{design}"] & c.seed.isin(seeds) & c.ds.isin(SMALL)]
    rows = []
    for (ds, seed, m), g in cs.groupby(["ds", "seed", "method"]):
        t = g[g.stability >= g.stability.max() - 0.01]
        rows.append({"ds": ds, "seed": int(seed), "method": m, "share_tied": len(t) / len(g),
                     "span": int(t.n_clusters.max() - t.n_clusters.min()), "k_floor": int(g.n_clusters.min())})
    R = pd.DataFrame(rows)
    pk = p[(p.S == "S0") & (p.design == design) & p.seed.isin(seeds) & p.ds.isin(SMALL)]
    w = pk.pivot_table(index=["ds", "seed", "method"], columns="reading", values=["n_clusters", "ari"])
    w.columns = [f"{a}_{'F' if r == F else 'M'}" for a, r in w.columns]
    R = R.merge(w.reset_index(), on=["ds", "seed", "method"]).dropna(subset=["ari_F", "ari_M"])
    R["k_gap"] = (R.n_clusters_F - R.n_clusters_M).abs()
    R["ari_gap"] = (R.ari_F - R.ari_M).abs()
    per = {}
    for (ds, m), g in R.groupby(["ds", "method"]):
        per[f"{ds}/{m}"] = {"share_tied": round(float(g.share_tied.median()), 3), "span": float(g.span.median()),
                            "k_F": kr(int(g.n_clusters_F.min()), int(g.n_clusters_F.max())),
                            "k_M": kr(int(g.n_clusters_M.min()), int(g.n_clusters_M.max())),
                            "ari_F": round(float(g.ari_F.median()), 3), "ari_M": round(float(g.ari_M.median()), 3),
                            "true_k": TYPES[ds], "n_runs": int(len(g))}
    dis = R[R.k_gap >= gap]
    out = {"per_dataset": per, "n_runs": int(len(R)), "gap": gap,
           "disagree": {"n": int(len(dis)), "n_mica": int((dis.method == "mica").sum()),
                        "mean_ari_gap": round(float(dis.ari_gap.mean()), 3)},
           "agree": {"n": int(len(R) - len(dis)), "mean_ari_gap": round(float(R[R.k_gap < gap].ari_gap.mean()), 3)}}
    for tag in ("F", "M"):
        fl = R[f"n_clusters_{tag}"] == R.k_floor
        out[f"floor_{tag}"] = {"n": int(fl.sum()), "mean_ari": round(float(R.loc[fl, f"ari_{tag}"].mean()), 3),
                               "n_interior": int((~fl).sum()), "mean_ari_interior": round(float(R.loc[~fl, f"ari_{tag}"].mean()), 3)}
    return out


def curve_figure(c, p, seeds, show=("buettner", "pollen", "zeisel"), design="mid", b=25):
    """Figure 1: stability and ARI against the number of clusters, from the confirmatory runs.

    Medians over the seeds of each count's stability (solid) and ARI (dashed), for MICA and Leiden on
    matched candidates, scored as the original protocol. The markers show the counts the two readings
    picked (median over the seeds), and the dotted line the true number of types.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    blue, orange, muted, ink = "#2a6f97", "#d1792b", "#8a8a85", "#2b2b2b"
    plt.rcParams.update({"font.size": 7, "font.family": "serif", "pdf.fonttype": 42})
    cs = c[(c.scorer == "ours") & (c.b == b) & c[f"in_{design}"] & c.seed.isin(seeds)]
    fig, axes = plt.subplots(1, len(show), figsize=(7.0, 1.95), sharey=True)
    for ax, ds in zip(axes, show):
        for m, col, lab in [("mica", blue, "MICA"), ("leiden_pca", orange, "Leiden")]:
            s = cs[(cs.ds == ds) & (cs.method == m)].groupby("n_clusters").agg(
                stab=("stability", "median"), ari=("ari", "median")).reset_index()
            if not len(s):
                continue
            ax.plot(s.n_clusters, s.stab, color=col, lw=1.5, label=f"{lab} stability")
            ax.plot(s.n_clusters, s.ari, color=col, lw=1.2, ls=(0, (2, 1.5)), label=f"{lab} ARI")
            for reading, mk in ((F, "^"), (M, "v")):
                k = p[(p.S == "S0") & (p.design == design) & (p.reading == reading) & (p.method == m)
                      & (p.ds == ds) & p.seed.isin(seeds)].n_clusters.median()
                if pd.notna(k):
                    y = float(np.interp(k, s.n_clusters, s.stab))
                    ax.plot([k], [y], marker=mk, color=col, ms=4, mec="white", mew=0.5, ls="none")
        ax.axvline(TYPES[ds], color=muted, lw=0.8, ls=":")
        ax.text(TYPES[ds], 1.04, "true", fontsize=6, color=muted, ha="center")
        ax.set_title(NAMES[ds], fontsize=7.5, color=ink)
        ax.set_xlabel("clusters")
        ax.set_ylim(-0.05, 1.1)
        ax.grid(axis="y", color="#e6e5e1", lw=0.5)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("score")
    axes[0].legend(fontsize=5.5, frameon=False, loc="lower left", ncol=2, columnspacing=0.8,
                   handlelength=1.5, labelspacing=0.2)
    from matplotlib.lines import Line2D
    picks = [Line2D([], [], marker="^", color=muted, ls="none", ms=4, label="finest stable pick"),
             Line2D([], [], marker="v", color=muted, ls="none", ms=4, label="most stable pick")]
    axes[1].legend(handles=picks, fontsize=5.5, frameon=False, loc="lower right", handletextpad=0.2,
                   labelspacing=0.2)
    fig.tight_layout(pad=0.3)
    fig.savefig(FIGS / "a4_stability_curve.pdf", metadata={"CreationDate": None})
    fig.savefig(FIGS / "a4_stability_curve.png", dpi=200)
    plt.close(fig)


def repick(curve_rows, settings):
    """stability.pick on a subset of one run's curve, in the given order."""
    sc = {}
    for s in settings:
        r = curve_rows.loc[s]
        sc[s] = {"stability": r.stability, "ci_low": r.ci_low, "ci_high": r.ci_high, "n_clusters": int(r.n_clusters)}
    return stability.pick(sc)


def decomposition(c, seeds, b=25):
    """E3: Leiden's mean under four candidate sets that change one thing at a time."""
    s = c[(c.scorer == "ours") & (c.b == b) & (c.method == "leiden_pca") & c.seed.isin(seeds) & c.ds.isin(SMALL)]
    out = []
    for (ds, seed), g in s.groupby(["ds", "seed"]):
        g = g.drop_duplicates("setting").set_index("setting").sort_index()
        orig = g[g.in_original]
        old_counts = sorted(orig.n_clusters.unique())
        by = {}
        for st, r in orig.iterrows():
            by.setdefault(int(r.n_clusters), []).append(st)
        orig_one = [v[(len(v) - 1) // 2] for k, v in sorted(by.items())]
        mid = g[g.in_mid].sort_values("n_clusters")
        mid_old = mid[mid.n_clusters.isin(old_counts)]
        sets = {"original": list(orig.index), "original_one_per_count": orig_one,
                "mid_old_counts": list(mid_old.index), "mid": list(mid.index)}
        for name, settings in sets.items():
            if not settings:
                continue
            for reading, info in repick(g, settings).items():
                out.append({"ds": ds, "seed": seed, "set": name, "reading": reading,
                            "ari": float(g.loc[info["chosen_setting"], "ari"]), "k": info["n_clusters"]})
    D = pd.DataFrame(out)
    if not len(D):
        return {}
    res = {}
    for (name, reading), g in D.groupby(["set", "reading"]):
        res[f"{name}/{reading}"] = round(float(g.groupby("ds").ari.median().mean()), 4)
    return res


def chooser_agreement(p, seeds, design="mid"):
    """E4: S1 against S3-G, finest reading: same chosen count, and the ARI change."""
    a = p[(p.S == "S1") & (p.reading == F) & (p.design == design) & p.seed.isin(seeds)]
    b = p[(p.S == "S3-G") & (p.design == design) & p.seed.isin(seeds)]
    j = a.merge(b, on=["ds", "seed", "method"], suffixes=("_s1", "_s3"))
    if not len(j):
        return {}
    return {m: {"runs": int(len(g)), "same_count": int((g.n_clusters_s1 == g.n_clusters_s3).sum()),
                "mean_abs_ari_change": round(float((g.ari_s1 - g.ari_s3).abs().mean()), 4)}
            for m, g in j.groupby("method")}


def fmt(x, nd=3):
    return "n/a" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.{nd}f}"


def kr(lo, hi):
    lo, hi = int(lo), int(hi)
    return f"{lo}" if lo == hi else f"{lo}\\textrm{{-}}{hi}"


def write_tables(p, c, out, conf, pseeds):
    """The manuscript's protocol v1 tables (latex/tables/a4_*.tex) and the per-seed figure."""
    TABLES.mkdir(parents=True, exist_ok=True)
    allseeds = list(range(10))
    # Table: confirmatory runs, matched candidates, S0
    rows = []
    for ds in SMALL:
        cells = []
        for reading in (F, M):
            for m in METHODS:
                s = p[(p.S == "S0") & (p.design == "mid") & (p.reading == reading) & (p.method == m)
                      & (p.ds == ds) & p.seed.isin(conf)]
                cells.append("n/a" if not len(s) else f"{s.ari.median():.3f}\\,({kr(s.n_clusters.min(), s.n_clusters.max())})")
        rows.append(f"{NAMES[ds]} ({TYPES[ds]}) & " + " & ".join(cells) + " \\\\")
    rows.append("\\midrule")
    for label, seeds in (("Mean, seeds 3 to 9", conf), ("Mean, seeds 0 to 2", [0, 1, 2]), ("Mean, all ten seeds", allseeds)):
        vals = [summary(p, "S0", "mid", r, m, seeds) for r in (F, M) for m in METHODS]
        rows.append(f"\\textit{{{label}}} & " + " & ".join("n/a" if v is None else
                    (f"\\textbf{{{v['mean']:.3f}}}" if label.startswith("Mean, seeds 3") else f"{v['mean']:.3f}") for v in vals) + " \\\\")
    (TABLES / "a4_confirm.tex").write_text(
        "\\begin{table*}[t]\n\\centering\n\\caption{The confirmatory runs of protocol v1: the stability rule on matched "
        "candidates at the seven new seeds, 3 to 9, scored as the original protocol (S0, the first 25 of the 100 "
        "subsamples). ARI of the chosen partition, median over the seven seeds, with the smallest and largest chosen "
        "count in parentheses and the number of cell types after each dataset name. No labels are used in the choice. "
        "The mean rows average the per-dataset medians over the nine datasets, for the confirmatory seeds, for the "
        "exploratory seeds 0 to 2, and for all ten.}\n\\label{tab:confirm}\n\\footnotesize\n"
        "\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{lcccccc}\n\\toprule\n"
        " & \\multicolumn{3}{c}{Finest stable} & \\multicolumn{3}{c}{Most stable} \\\\\n"
        "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\nDataset (types) & MICA & Leiden & $k$-means & MICA & Leiden & $k$-means \\\\\n"
        "\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n\\end{table*}\n")

    # Table: every design and scorer, means over the nine datasets at the confirmatory seeds
    spec = [("Original grids", "original", "S0"), ("Matched, middle setting", "mid", "S0"),
            ("Matched, lowest setting", "low", "S0"), ("Matched, highest setting", "high", "S0"),
            ("Matched, every setting", "dense", "S0"), ("Matched, 100 subsamples (S1)", "mid", "S1"),
            ("Original, chooseR's code (S3)", "original", "S3-G"), ("Matched, chooseR's code (S3)", "mid", "S3-G")]
    rows = []
    for label, d, S in spec:
        vals = []
        for r in (F, M):
            for m in METHODS:
                v = None if (S.startswith("S3") and r == M) else summary(p, S, d, r, m, conf)
                vals.append("n/a" if (S.startswith("S3") and r == M) else fmt(v["mean"] if v else None))
        rows.append(f"{label} & " + " & ".join(vals) + " \\\\")
    (TABLES / "a4_designs.tex").write_text(
        "\\begin{table*}[t]\n\\centering\n\\caption{How far each method's mean ARI moves with the candidate design and "
        "the scorer. Means over the nine small datasets of the per-dataset medians, confirmatory seeds 3 to 9. The "
        "lowest, middle and highest setting per count and every setting differ only for Leiden, since $k$ has one "
        "setting per count. S1 scores the same partitions on all 100 subsamples. S3 is chooseR's own R code on all "
        "100, with a point interval where its bootstrap gives none (the protocol's primary repair). chooseR has only "
        "the finest-stable reading, and the most-stable reading uses no interval, so it equals S1.}\n"
        "\\label{tab:designs}\n\\footnotesize\n\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{lcccccc}\n\\toprule\n"
        " & \\multicolumn{3}{c}{Finest stable} & \\multicolumn{3}{c}{Most stable} \\\\\n"
        "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\nDesign and scorer & MICA & Leiden & $k$-means & MICA & Leiden & $k$-means \\\\\n"
        "\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n\\end{table*}\n")

    # Table: PBMC14K, medians over its seeds
    rows = []
    for m in METHODS:
        for S, label in (("S0", "S0"), ("S1", "S1"), ("S2", "S2, all cells"), ("S3-G", "S3, chooseR")):
            cells = []
            for d in ("original", "mid"):
                for r in (F, M):
                    s = p[(p.ds == "pbmc14k") & (p.S == S) & (p.design == d) & (p.reading == r) & (p.method == m)
                          & p.seed.isin(pseeds)]
                    cells.append("n/a" if (S == "S3-G" and r == M) else ("n/a" if not len(s) else
                                 f"{s.ari.median():.3f}\\,({kr(s.n_clusters.min(), s.n_clusters.max())})"))
            rows.append(f"{MNAME[m] if S == 'S0' else ''} & {label} & " + " & ".join(cells) + " \\\\")
        rows.append("\\midrule")
    (TABLES / "a4_pbmc.tex").write_text(
        "\\begin{table*}[t]\n\\centering\n\\caption{PBMC14K, the one large benchmark dataset and the only one in MICA's "
        "graph mode, at five seeds. ARI of the chosen partition, median over the seeds, with the range of chosen counts. "
        "MICA's node2vec and Louvain take no seed, so every MICA run is a new draw. S0 keeps the co-clustering on a "
        "random 3,000 of the 13,605 cells and uses the first 25 subsamples, S1 the same cells and all 100, S2 all cells "
        "and all 100, S3 chooseR's code on all cells and all 100. The true number of types is 7.}\n"
        "\\label{tab:pbmc}\n\\footnotesize\n\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{llcccc}\n\\toprule\n"
        " & & \\multicolumn{2}{c}{Original grids} & \\multicolumn{2}{c}{Matched candidates} \\\\\n"
        "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\nMethod & Scorer & Finest & Most & Finest & Most \\\\\n\\midrule\n"
        + "\n".join(rows[:-1]) + "\n\\bottomrule\n\\end{tabular}\n\\end{table*}\n")

    # Figure: per-seed means, the variation the draw alone causes
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 4, figsize=(7.0, 1.9), sharey=True)
        colors = {"mica": "#1f77b4", "leiden_pca": "#d62728", "kmeans_pca": "#2ca02c"}
        for ax, (d, r, title) in zip(axes, (("original", F, "Original, finest"), ("original", M, "Original, most"),
                                             ("mid", F, "Matched, finest"), ("mid", M, "Matched, most"))):
            for m in METHODS:
                v = summary(p, "S0", d, r, m, allseeds)
                if v and v["per_seed"]:
                    xs = sorted(v["per_seed"])
                    ax.plot(xs, [v["per_seed"][x] for x in xs], "o-", ms=2.5, lw=1, color=colors[m],
                            label=MNAME[m].replace("$", ""))
            ax.axvspan(-0.5, 2.5, color="0.9", lw=0)
            ax.set_title(title, fontsize=7)
            ax.set_xticks(range(0, 10, 3))
            ax.tick_params(labelsize=6)
            ax.set_xlabel("seed", fontsize=6)
        axes[0].set_ylabel("mean ARI, nine datasets", fontsize=6)
        axes[0].legend(fontsize=5.5, frameon=False, loc="lower left")
        fig.tight_layout(pad=0.3)
        FIGS.mkdir(parents=True, exist_ok=True)
        fig.savefig(FIGS / "a4_seed_means.pdf", metadata={"CreationDate": None})
        plt.close(fig)
    except Exception as e:  # the numbers do not depend on the figure
        print("figure skipped:", e)


def write_tables2(p, out, conf):
    """Criteria verdicts, the decomposition of Leiden's change, and chooseR's failures."""
    def v(x, k="verdict"):
        return (out.get(x) or {}).get(k, "n/a")
    c1f, c1m = out["C1"][F], out["C1"][M]
    a3 = exploratory_values()
    rows = [
        ("C1, finest stable: MICA above Leiden", f"{fmt(a3.get('F.mica'))} vs {fmt(a3.get('F.leiden'))}", f"{fmt(c1f.get('mean_a'))} vs {fmt(c1f.get('mean_b'))}, leads at {c1f.get('a_leads_seeds', 'n/a')} seeds", c1f["verdict"]),
        ("C1, most stable: MICA above Leiden", f"{fmt(a3.get('M.mica'))} vs {fmt(a3.get('M.leiden'))}", f"{fmt(c1m.get('mean_a'))} vs {fmt(c1m.get('mean_b'))}, leads at {c1m.get('a_leads_seeds', 'n/a')} seeds", c1m["verdict"]),
        ("C2, the reading moves MICA more", f"{fmt(a3.get('C2.mica'))} vs {fmt(a3.get('C2.leiden'))}",
         f"{fmt(out['C2']['mean_abs_reading_change'].get('mica'))} vs {fmt(out['C2']['mean_abs_reading_change'].get('leiden_pca'))}", out["C2"]["verdict"]),
    ]
    c3 = out.get("C3", {})
    if c3:
        rows.append(("C3, MICA's flat top in the five smallest", f"{fmt(a3.get('C3.five.mica'))} vs {fmt(a3.get('C3.four.mica'))}, Leiden {fmt(a3.get('C3.five.leiden'))}",
                     f"{fmt(c3.get('five_smallest/mica', {}).get('share_tied'))} vs {fmt(c3.get('four_larger/mica', {}).get('share_tied'))}, "
                     f"Leiden {fmt(c3.get('five_smallest/leiden_pca', {}).get('share_tied'))}", c3["verdict"]))
    c4 = out.get("C4", {})
    rows.append(("C4, PBMC14K tie under most stable, seeds 0 to 4", "3 of 3 seeds", c4.get("ties", "n/a"), c4.get("verdict", "n/a")))
    c5 = out["C5"]
    rows.append(("C5, original grids: Leiden above MICA, most stable", f"{fmt(a3.get('C5.leiden'))} vs {fmt(a3.get('C5.mica'))}",
                 f"{fmt(c5.get('mean_a'))} vs {fmt(c5.get('mean_b'))}", c5["verdict"]))
    f1 = out["F1"]
    rows.append(("F1, chooseR's code: MICA above Leiden, finest", "not run", f"{fmt(f1.get('mean_a'))} vs {fmt(f1.get('mean_b'))}", f1["verdict"]))
    fc = out.get("FC", {})
    rows.append(("FC, chooseR's tally equals ours", "not run", f"{fc['identical']:,} of {fc['matrices']:,}" if 'identical' in fc else "n/a", fc.get("verdict", "n/a")))
    body = "\n".join(f"{a} & {b} & {c} & {d} \\\\" for a, b, c, d in rows)
    (TABLES / "a4_criteria.tex").write_text(
        "\\begin{table*}[t]\n\\centering\n\\caption{The criteria fixed in protocol v1 before the confirmatory runs, with "
        "the exploratory value (seeds 0 to 2) and the value at the seven new seeds, 3 to 9, except C4 (PBMC14K, seeds 0 to 4) "
        "and FC (every seed). Mean ARI over the nine small "
        "datasets of the per-dataset medians, on matched candidates scored as the original protocol (S0), unless stated. "
        "A criterion of C1 holds if the mean is higher and MICA leads at 5 or more of the 7 per-seed means, and is weak "
        "if only the mean is higher.}\n\\label{tab:criteria}\n\\footnotesize\n\\setlength{\\tabcolsep}{4pt}\n"
        "\\begin{tabular}{lllc}\n\\toprule\nCriterion & Exploratory & Protocol v1 & Verdict \\\\\n\\midrule\n"
        + body + "\n\\bottomrule\n\\end{tabular}\n\\end{table*}\n")
    e3 = out.get("E3_decomposition_leiden", {})
    steps = [("Original grid, capped at 2 to 20", "original"), ("Original grid, one setting per count", "original_one_per_count"),
             ("Matched grid, only counts the original grid reached", "mid_old_counts"), ("Matched grid, all matched counts", "mid")]
    body = "\n".join(f"{lab} & {fmt(e3.get(f'{k}/{F}'))} & {fmt(e3.get(f'{k}/{M}'))} \\\\" for lab, k in steps)
    (TABLES / "a4_decomp.tex").write_text(
        "\\begin{table}[t]\n\\centering\n\\caption{Why matching changes Leiden's score: its mean ARI over the nine small "
        "datasets as the candidate set changes one step at a time, seeds 3 to 9, original protocol (S0). Each step is "
        "rescored from the same partitions, so only the candidates differ.}\n\\label{tab:decomp}\n\\footnotesize\n"
        "\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{lcc}\n\\toprule\nLeiden's candidates & Finest & Most \\\\\n\\midrule\n"
        + body + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n")


def write_numbers(p, out, conf, pseeds, root):
    """latex/numbers_a4.tex: every protocol v1 figure the text quotes, as \V{key}, from this run's results.

    The text never types a protocol v1 number by hand. \V{key} prints the value, and an unknown
    key prints ??key in bold, so a stale or missing number shows in the PDF.
    """
    vals = {}
    short_m = {"mica": "mica", "leiden_pca": "leiden", "kmeans_pca": "kmeans"}
    short_r = {F: "F", M: "M"}
    for S in ("S0", "S1", "S3-G", "S3-G2"):
        for d in protocol_designs():
            for r in (F, M):
                for m in METHODS:
                    for tag, seeds in (("", conf), (".s02", [0, 1, 2]), (".all", list(range(10)))):
                        v = summary(p, S, d, r, m, seeds)
                        if v:
                            vals[f"{S.replace('-', '')}.{d}.{short_r[r]}.{short_m[m]}{tag}"] = f"{v['mean']:.3f}"
    for r in (F, M):
        c = out["C1"][r]
        if "a_leads_seeds" in c:
            vals[f"C1.{short_r[r]}.leads"] = c["a_leads_seeds"]
            vals[f"C1.{short_r[r]}.verdict"] = c["verdict"]
        for m in ("mica", "leiden_pca"):
            ps = summary(p, "S0", "mid", r, m, conf)
            if ps and ps["per_seed"]:
                vals[f"S0.mid.{short_r[r]}.{short_m[m]}.seedmin"] = f"{min(ps['per_seed'].values()):.3f}"
                vals[f"S0.mid.{short_r[r]}.{short_m[m]}.seedmax"] = f"{max(ps['per_seed'].values()):.3f}"
    for m, v in out["C2"]["mean_abs_reading_change"].items():
        if v is not None:
            vals[f"C2.{short_m[m]}"] = f"{v:.3f}"
    for k, val in exploratory_values().items():
        vals[f"A3.{k}"] = f"{val:.3f}"
    for m in ("mica", "leiden_pca"):  # the same quantity at the exploratory seeds, which R0 reproduces
        v = reading_change(p, "S0", "mid", m, [0, 1, 2])
        if v is not None:
            vals[f"C2.{short_m[m]}.s02"] = f"{v:.3f}"
    for r in (F, M):
        x, y = summary(p, "S0", "mid", r, "mica", [0, 1, 2]), summary(p, "S0", "mid", r, "leiden_pca", [0, 1, 2])
        if x and y:
            vals[f"C1.{short_r[r]}.gap.s02"] = f"{x['mean'] - y['mean']:.3f}"
    for k, v in out.get("C3", {}).items():
        if isinstance(v, dict) and "share_tied" in v:
            sz, m = k.split("/")
            vals[f"C3.{sz}.{short_m[m]}.share"] = f"{v['share_tied']:.3f}"
            vals[f"C3.{sz}.{short_m[m]}.span"] = f"{v['span']:.1f}"
            vals[f"C3.{sz}.{short_m[m]}.atone"] = f"{100 * v['share_at_one']:.0f}"
    lf = out.get("label_free_signals")
    if lf:
        for key, v in lf["per_dataset"].items():
            ds, m = key.split("/")
            base = f"L.{ds}.{short_m[m]}"
            vals[f"{base}.share"] = f"{100 * v['share_tied']:.0f}"
            vals[f"{base}.span"] = f"{v['span']:.0f}"
            vals[f"{base}.kF"], vals[f"{base}.kM"] = v["k_F"], v["k_M"]
            vals[f"{base}.ariF"], vals[f"{base}.ariM"] = f"{v['ari_F']:.3f}", f"{v['ari_M']:.3f}"
            vals[f"{base}.truek"] = str(v["true_k"])
        vals["L.runs"], vals["L.gap"] = str(lf["n_runs"]), str(lf["gap"])
        vals["L.disagree.n"], vals["L.disagree.mica"] = str(lf["disagree"]["n"]), str(lf["disagree"]["n_mica"])
        vals["L.disagree.arigap"] = f"{lf['disagree']['mean_ari_gap']:.3f}"
        vals["L.agree.n"], vals["L.agree.arigap"] = str(lf["agree"]["n"]), f"{lf['agree']['mean_ari_gap']:.3f}"
        for tag in ("F", "M"):
            fl = lf[f"floor_{tag}"]
            vals[f"L.floor.{tag}.n"], vals[f"L.floor.{tag}.ari"] = str(fl["n"]), f"{fl['mean_ari']:.3f}"
            vals[f"L.interior.{tag}.n"], vals[f"L.interior.{tag}.ari"] = str(fl["n_interior"]), f"{fl['mean_ari_interior']:.3f}"
    for k in ("C1", "C2", "C3", "C4", "C5", "F1", "FC"):
        if k in out and isinstance(out[k], dict) and "verdict" in out[k] and k != "C1":
            vals[f"{k}.verdict"] = out[k]["verdict"]
    if "C4" in out:
        vals["C4.ties"] = out["C4"]["ties"]
    if "FC" in out:
        vals["FC.identical"] = f"{out['FC']['identical']:,}"
        vals["FC.matrices"] = f"{out['FC']['matrices']:,}"
    for k, v in out.get("E3_decomposition_leiden", {}).items():
        name, r = k.split("/")
        vals[f"E3.{name.replace('_', '')}.{short_r[r]}"] = f"{v:.3f}"
    for k, v in out.get("E2_menu_vs_representation", {}).items():
        d, r = k.split("/")
        for kk, x in v.items():
            if x is not None:
                vals[f"E2.{d}.{short_r[r]}.{kk.replace('_', '')}"] = f"{x:.3f}"
                if kk == "kmeans_minus_leiden":
                    vals[f"E2.{d}.{short_r[r]}.leidenminuskmeans"] = f"{-x:.3f}"
    # Leiden's range over the candidate variants
    for r in (F, M):
        xs = [summary(p, "S0", d, r, "leiden_pca", conf) for d in ("mid", "low", "high", "dense")]
        xs = [x["mean"] for x in xs if x]
        if xs:
            vals[f"variants.{short_r[r]}.leiden.min"] = f"{min(xs):.3f}"
            vals[f"variants.{short_r[r]}.leiden.max"] = f"{max(xs):.3f}"
    reg = root / "regression.csv"
    if reg.exists():
        rg = pd.read_csv(reg)
        rg = rg[~rg.picks.astype(str).str.startswith("no saved")]
        bad = rg[[c for c in ("picks", "stability_curve", "full_sweep") if c in rg]].map(
            lambda v: str(v) not in ("equal", "absent in saved file", "nan")
            and not str(v).startswith("not expected")).any(axis=1)
        vals["R0.compared"], vals["R0.mismatch"] = str(len(rg)), str(int(bad.sum()))
        small = ~rg.folder.str.contains("pbmc14k")
        vals["R0.small.compared"], vals["R0.small.mismatch"] = str(int(small.sum())), str(int(bad[small].sum()))
        lp = rg.folder.str.contains("pbmc14k") & (rg.note.astype(str) == "Leiden only")
        vals["R0.pbmcleiden.compared"], vals["R0.pbmcleiden.mismatch"] = str(int(lp.sum())), str(int(bad[lp].sum()))
    e3a = out.get("E3_added_counts", {})
    for k, v in e3a.items():
        vals[f"E3.{k.replace('_', '')}"] = str(v)
    f3 = out.get("S3_failures")
    if f3:
        vals["S3.settings"], vals["S3.failed"] = f"{f3['settings']:,}", f"{f3['failed']:,}"
        for m, x in f3["runs_without_choice_V"].items():
            vals[f"S3.nochoice.{short_m[m]}"] = x
        for ds, x in f3["by_dataset"].items():
            vals[f"S3.failed.{ds}"] = x
    for d, block in out.get("E4_chooser_agreement", {}).items():
        for m, x in block.items():
            vals[f"E4.{d}.{short_m[m]}.same"] = f"{x['same_count']} of {x['runs']}"
            vals[f"E4.{d}.{short_m[m]}.change"] = f"{x['mean_abs_ari_change']:.3f}"
    for k, v in out.get("ami_means", {}).items():
        r, m = k.split("/")
        vals[f"AMI.{short_r[r]}.{short_m[m]}"] = f"{v:.3f}"
    for k, v in out.get("seed_range", {}).items():
        r, m = k.split("/")
        vals[f"range.{short_r[r]}.{short_m[m]}.mean"] = f"{v['mean']:.2f}"
        vals[f"range.{short_r[r]}.{short_m[m]}.max"] = f"{v['max']:.2f}"
        vals[f"range.{short_r[r]}.{short_m[m]}.maxds"] = NAMES.get(v["max_dataset"], v["max_dataset"])
    for k, x in out.get("E5_dataset_robustness", {}).items():
        S, d, r = k.split("/")
        key = f"E5.{S.replace('-', '')}.{d}.{short_r[r]}"
        vals[key + ".micawins"] = str(x["mica_wins"])
        vals[key + ".leidenwins"] = str(x["leiden_wins"])
        vals[key + ".level"] = str(x["level"])
        vals[key + ".loo"] = f"{x['loo_gap']:.3f}"
        vals[key + ".loods"] = x["loo_gap_without"]
        vals[key + ".flips"] = str(x["sign_flips_without_one"])
        vals[key + ".seedleads"] = x["mica_leads_seeds"]
    for k, x in out.get("E5_mica_vs_kmeans", {}).items():
        S, d, r = k.split("/")
        key = f"E5k.{S.replace('-', '')}.{d}.{short_r[r]}"
        vals[key + ".micawins"] = str(x["mica_wins"])
        vals[key + ".kmeanswins"] = str(x["kmeans_wins"])
        vals[key + ".level"] = str(x["level"])
        vals[key + ".loo"] = f"{x['loo_gap']:.3f}"
        vals[key + ".loods"] = x["loo_gap_without"]
    # PBMC14K medians over its seeds, every scorer, design and reading
    pb = p[(p.ds == "pbmc14k") & p.seed.isin(pseeds)]
    for (S, d, r, m), g in pb.groupby(["S", "design", "reading", "method"]):
        key = f"pbmc.{S.replace('-', '')}.{d}.{short_r.get(r, 'F')}.{short_m[m]}"
        if S == "S3-V":  # as written, chooseR's code can make no choice, so count those runs instead
            vals[key + ".nochoice"] = f"{int(g.ari.isna().sum())} of {len(g)}"
        g = g.dropna(subset=["ari"])
        if not len(g):
            continue
        vals[key] = f"{g.ari.median():.3f}"
        vals[key + ".k"] = kr(g.n_clusters.min(), g.n_clusters.max())
        vals[key + ".seeds"] = str(int(g.seed.nunique()))
    # per-subsample cost from the provenance files
    prov = [json.loads(Path(f).read_text()) for f in glob(str(root / "provenance" / "paper_*.json"))]
    for m in METHODS:
        xs = [d["seconds_per_subsample"][m] for d in prov if "seconds_per_subsample" in d and m in d["seconds_per_subsample"]
              and not str(d.get("dataset")).endswith("pbmc14k")]
        if xs:
            vals[f"cost.small.{short_m[m]}.median"] = f"{float(np.median(xs)):.1f}"
    lines = ["% Generated by scripts/a4_protocol_v1.py from results/v1. Do not edit by hand.",
             f"% seeds {conf[0]} to {conf[-1]} confirmatory, PBMC14K seeds {pseeds}"]
    for k in sorted(vals):
        lines.append(f"\\expandafter\\def\\csname v@{k}\\endcsname{{{vals[k]}}}")
    (RESULTS.parent / "latex" / "numbers_a4.tex").write_text("\n".join(lines) + "\n")
    return vals


def protocol_designs():
    return ("original", "mid", "low", "high", "dense")


def seeds_of(arg, default):
    if not arg:
        return default
    a, b = arg.split("-")
    return list(range(int(a), int(b) + 1))


def parse_args(argv):
    """--tables is a switch. The old form '--tables 1' still works. Every other option takes a value."""
    opts, i = {}, 0
    while i < len(argv):
        if argv[i] == "--tables":
            opts["--tables"] = True
            i += 2 if i + 1 < len(argv) and not argv[i + 1].startswith("--") else 1
        else:
            opts[argv[i]] = argv[i + 1] if i + 1 < len(argv) else None
            i += 2
    return opts


def main(argv):
    opts = parse_args(argv)
    global V1
    if "--root" in opts:
        V1 = Path(opts["--root"])
    conf = seeds_of(opts.get("--seeds"), list(range(3, 10)))
    pseeds = seeds_of(opts.get("--pbmc-seeds"), list(range(0, 5)))
    allseeds = list(range(0, 10))
    p = pd.concat([load("picks"), load("picks_chooser")], ignore_index=True)
    c = load("curves")
    out = {"seeds_confirmatory": conf, "pbmc_seeds": pseeds,
           "runs_present": {f"{ds}": sorted(int(s) for s in g.seed.unique()) for ds, g in p.groupby("ds")}}

    # R0
    reg = V1 / "regression.csv"
    if reg.exists():
        r = pd.read_csv(reg)
        out["R0"] = r.fillna("").to_dict(orient="records")
    # C1 and C5
    out["C1"] = {}
    for reading in (F, M):
        lr = lead(p, "S0", "mid", reading, "mica", "leiden_pca", conf)
        out["C1"][reading] = {**(lr or {}), "verdict": verdict_c1(lr, len(conf))}
    lr = lead(p, "S0", "original", M, "leiden_pca", "mica", conf)
    out["C5"] = {**(lr or {}), "verdict": "not yet computable" if lr is None else ("holds" if lr["gap"] > 0 else "fails")}
    # C2
    rc = {m: reading_change(p, "S0", "mid", m, conf) for m in ("mica", "leiden_pca", "kmeans_pca")}
    out["C2"] = {"mean_abs_reading_change": rc,
                 "verdict": "not yet computable" if None in (rc["mica"], rc["leiden_pca"])
                 else ("holds" if rc["mica"] > rc["leiden_pca"] else "fails")}
    # C3
    ft = flat_top(c, conf)
    if ft:
        a, b_, cl = (ft.get("five_smallest/mica", {}).get("share_tied"), ft.get("four_larger/mica", {}).get("share_tied"),
                     ft.get("five_smallest/leiden_pca", {}).get("share_tied"))
        ok = None if None in (a, b_, cl) else (a > b_ and a > cl)
        out["C3"] = {**ft, "verdict": "not yet computable" if ok is None else ("holds" if ok else "fails")}
    out["label_free_signals"] = label_free_signals(p, c, conf)
    # C4
    pb = p[(p.ds == "pbmc14k") & (p.S == "S0") & (p.design == "mid") & (p.reading == M) & p.seed.isin(pseeds)]
    if len(pb):
        w = pb.pivot_table(index="seed", columns="method", values="ari")
        if {"mica", "leiden_pca"} <= set(w.columns):
            w = w.dropna(subset=["mica", "leiden_pca"])
            ties = int(((w.mica - w.leiden_pca).abs() <= 0.02).sum())
            out["C4"] = {"per_seed": w.round(4).to_dict(orient="index"), "ties": f"{ties} of {len(w)}",
                         "verdict": "not yet computable" if len(w) < len(pseeds) else ("holds" if ties >= 3 else "fails")}
    # F1
    lr = lead(p, "S3-G", "mid", F, "mica", "leiden_pca", conf)
    out["F1"] = {**(lr or {}), "verdict": "not yet computable" if lr is None else ("holds" if lr["gap"] > 0 else "fails")}
    # FC
    cc = sorted(glob(str(V1 / "chooser_r" / "crosscheck_*.csv")))
    if cc:
        x = pd.concat([pd.read_csv(f) for f in cc])
        out["FC"] = {"matrices": int(len(x)), "identical": int(x.equal.sum()), "max_abs_diff": float(x.max_abs_diff.max()),
                     "verdict": "holds" if bool(x.equal.all()) else "fails"}
    # summaries for the tables: every scorer, design, reading and method, confirmatory and all seeds
    tab = {}
    for (S, design, reading, method), _ in p.groupby(["S", "design", "reading", "method"]):
        for label, seeds in (("conf", conf), ("all", allseeds)):
            s = summary(p, S, design, reading, method, seeds)
            if s:
                tab[f"{S}/{design}/{reading}/{method}/{label}"] = {k: s[k] for k in ("mean", "n_datasets", "n_runs", "per_seed")}
    out["summaries"] = tab
    out["E2_menu_vs_representation"] = {
        f"{d}/{r}": {"kmeans_minus_leiden": (lead(p, "S0", d, r, "kmeans_pca", "leiden_pca", conf) or {}).get("gap"),
                     "mica_minus_kmeans": (lead(p, "S0", d, r, "mica", "kmeans_pca", conf) or {}).get("gap")}
        for d in ("original", "mid") for r in (F, M)}
    out["E3_decomposition_leiden"] = decomposition(c, conf)
    # Exploratory, added Nov 13 after reading Table a4_confirm: how far the mean gap rests on single datasets.
    # Per-dataset medians over the confirmatory seeds, MICA minus Leiden. Wins count gaps above 0.01.
    rob = {}
    for S, d, r in (("S0", "mid", F), ("S0", "mid", M), ("S1", "mid", F), ("S1", "mid", M), ("S3-G", "mid", F),
                    ("S0", "original", F), ("S0", "original", M), ("S1", "original", F), ("S3-G", "original", F)):
        x, y = summary(p, S, d, r, "mica", conf), summary(p, S, d, r, "leiden_pca", conf)
        if not x or not y:
            continue
        g = {ds: x["per_dataset"][ds] - y["per_dataset"][ds] for ds in SMALL}
        loo = {ds: float(np.mean([g[o] for o in SMALL if o != ds])) for ds in SMALL}
        worst = min(loo, key=loo.get) if sum(g.values()) > 0 else max(loo, key=loo.get)
        rob[f"{S}/{d}/{r}"] = {"gap": round(float(np.mean(list(g.values()))), 4),
                               "mica_leads_seeds": lead(p, S, d, r, "mica", "leiden_pca", conf)["a_leads_seeds"],
                               "mica_wins": int(sum(v > 0.01 for v in g.values())),
                               "leiden_wins": int(sum(v < -0.01 for v in g.values())),
                               "level": int(sum(abs(v) <= 0.01 for v in g.values())),
                               "loo_gap_without": NAMES[worst], "loo_gap": round(loo[worst], 4),
                               "sign_flips_without_one": int(sum(np.sign(v) != np.sign(np.mean(list(g.values())))
                                                                 for v in loo.values())),
                               "per_dataset_gap": {k: round(v, 4) for k, v in g.items()}}
    out["E5_dataset_robustness"] = rob
    # The same check for MICA against PCA plus k-means, the representation comparison (E2)
    robk = {}
    for S, d, r in (("S0", "mid", F), ("S0", "mid", M), ("S0", "original", F), ("S0", "original", M)):
        x, y = summary(p, S, d, r, "mica", conf), summary(p, S, d, r, "kmeans_pca", conf)
        if not x or not y:
            continue
        g = {ds: x["per_dataset"][ds] - y["per_dataset"][ds] for ds in SMALL}
        loo = {ds: float(np.mean([g[o] for o in SMALL if o != ds])) for ds in SMALL}
        worst = min(loo, key=loo.get)
        robk[f"{S}/{d}/{r}"] = {"gap": round(float(np.mean(list(g.values()))), 4),
                                "mica_wins": int(sum(v > 0.01 for v in g.values())),
                                "kmeans_wins": int(sum(v < -0.01 for v in g.values())),
                                "level": int(sum(abs(v) <= 0.01 for v in g.values())),
                                "loo_gap_without": NAMES[worst], "loo_gap": round(loo[worst], 4),
                                "per_dataset_gap": {k: round(v, 4) for k, v in g.items()}}
    out["E5_mica_vs_kmeans"] = robk
    # which counts the matched design adds for Leiden, and whether the plain maximum takes them
    lc = c[(c.scorer == "ours") & (c.b == 25) & (c.method == "leiden_pca") & c.seed.isin(conf) & c.ds.isin(SMALL)]
    old_counts = {k: set(g[g.in_original].n_clusters) for k, g in lc.groupby(["ds", "seed"])}
    added = [k for key, g in lc.groupby(["ds", "seed"]) for k in set(g[g.in_mid].n_clusters) - old_counts[key]]
    q = p[(p.S == "S0") & (p.design == "mid") & (p.reading == M) & (p.method == "leiden_pca") & p.seed.isin(conf) & p.ds.isin(SMALL)]
    new_pick = sum(1 for r in q.itertuples() if r.n_clusters not in old_counts.get((r.ds, r.seed), set()))
    med = q.groupby("ds").n_clusters.median()
    out["E3_added_counts"] = {"added_at_most_3": int(sum(k <= 3 for k in added)), "added_above_3": int(sum(k > 3 for k in added)),
                              "plain_max_picks_new_count": f"{new_pick} of {len(q)}",
                              "median_pick_2_or_3": f"{int(((med >= 2) & (med <= 3)).sum())} of {len(med)}"}
    out["E4_chooser_agreement"] = {"mid": chooser_agreement(p, allseeds, "mid"),
                                   "original": chooser_agreement(p, allseeds, "original")}
    st = p[p.scorer.str.startswith("chooseR")]
    if len(st):
        out["S3_failed_settings"] = {m: int(g.n_failed.max()) for m, g in st.groupby("method")}
        out["S3_no_choice_runs_V"] = {m: int(g.chosen_setting.isna().sum()) for m, g in st[st.S == "S3-V"].groupby("method")}
    if "--tables" in opts:
        write_tables(p, c, out, conf, pseeds)
        curve_figure(c, p, conf)
    files = [f for f in sorted(glob(str(V1 / "chooser_r" / "paper_*_seed*.csv")))]
    cr = pd.concat([pd.read_csv(f).assign(ds=short(Path(f).stem.rsplit("_seed", 1)[0]),
                                          seed=int(Path(f).stem.rsplit("_seed", 1)[1])) for f in files],
                   ignore_index=True) if files else pd.DataFrame()
    if len(cr):
        cr["failed"] = cr.status.isin(["error", "no_interval"])
        out["S3_failures"] = {
            "settings": int(len(cr)), "failed": int(cr.failed.sum()),
            "by_status": cr.status.value_counts().to_dict(),
            "by_dataset": {k: f"{int(g.failed.sum())} of {len(g)}" for k, g in cr.groupby("ds")},
            "runs_without_choice_V": {m: f"{int(g.chosen_setting.isna().sum())} of {len(g)}"
                                      for m, g in p[p.S == "S3-V"].groupby("method")}}
    sm = p[(p.S == "S0") & (p.design == "mid") & p.ds.isin(SMALL) & p.seed.isin(conf)]
    if len(sm):
        out["ami_means"] = {f"{r}/{m}": round(float(g.groupby("ds").ami.median().mean()), 4)
                            for (r, m), g in sm.groupby(["reading", "method"])}
        rng_ = sm.groupby(["reading", "method", "ds"]).ari.agg(lambda x: x.max() - x.min())
        out["seed_range"] = {f"{r}/{m}": {"mean": round(float(g.mean()), 3), "max": round(float(g.max()), 3),
                                          "max_dataset": g.idxmax()[2]} for (r, m), g in rng_.groupby(level=[0, 1])}
    if "--tables" in opts:
        write_tables2(p, out, conf)
        out["latex_values"] = write_numbers(p, out, conf, pseeds, V1)
    (RESULTS / "a4_protocol_v1_numbers.json").write_text(json.dumps(out, indent=1, default=str))
    for k in ("C1", "C2", "C3", "C4", "C5", "F1", "FC"):
        if k in out:
            v = out[k]
            print(k, v.get("verdict") if isinstance(v, dict) and "verdict" in v else "",
                  {kk: vv for kk, vv in (v.items() if isinstance(v, dict) else []) if kk not in ("per_seed",)} if k != "C3" else "")
    print("wrote", RESULTS / "a4_protocol_v1_numbers.json")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
