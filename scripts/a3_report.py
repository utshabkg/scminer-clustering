"""Tables and figures for the Assignment 3 manuscript, built from results/.

Usage: python scripts/a3_report.py
Writes latex/tables/a3_*.tex, latex/figures/a3_*.pdf and
results/a3_report_numbers.json, which holds every number the text quotes so it can
be checked. Nothing here reruns an experiment.
"""
import json
import sys
from glob import glob

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scminer_clust.paths import REPO, RESULTS

TABLES, FIGS = REPO / "latex" / "tables", REPO / "latex" / "figures"
ORDER = ["yan", "goolam", "buettner", "pollen", "chung", "usoskin", "kolodziejczyk", "klein", "zeisel", "pbmc14k"]
NAMES = {"yan": "Yan", "goolam": "Goolam", "buettner": "Buettner", "pollen": "Pollen", "chung": "Chung",
         "usoskin": "Usoskin", "kolodziejczyk": "Kolodziejczyk", "klein": "Klein", "zeisel": "Zeisel",
         "pbmc14k": "PBMC14K"}
RULES = [("default", "Default"), ("silhouette_shared", "Silhouette"),
         ("stability_finest", "Finest stable"), ("stability_most", "Most stable")]
BLUE, ORANGE, INK, MUTED = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"


def short(ds):
    return ds.replace("paper_", "")


def load_stability():
    rows = []
    patterns = ["rule_stability_*", "rule_silhouette_*", "rule_default_*"]
    for f in sorted(f for pat in patterns for f in glob(str(RESULTS / pat / "metrics.csv"))):
        d = pd.read_csv(f)
        if "seed" not in d:
            d["seed"] = 0
        rows.append(d)
    if not rows:
        return pd.DataFrame()
    d = pd.concat(rows, ignore_index=True)
    # Leiden's shipped default was scored in the Assignment 2 runs and lives in the tracker
    t = pd.read_csv(REPO / "docs" / "results_tracker.csv", low_memory=False)
    lei = t[(t.k_setting == "default_resolution_1") & (t.gene_setting == "all genes")
            & t.dataset.str.startswith("paper_") & (t.method == "leiden_pca")].copy()
    lei["k_setting"] = "default"
    d = pd.concat([d, lei], ignore_index=True)
    d["ds"] = d.dataset.map(short)
    return d


def agg(d, rule):
    """Median and range of ARI and chosen cluster count over the repeats of a rule."""
    s = d[d.k_setting == rule]
    if not len(s):
        return pd.DataFrame()
    g = s.groupby(["ds", "method"]).agg(
        ari=("ari", "median"), ari_min=("ari", "min"), ari_max=("ari", "max"),
        ami=("ami", "median"), k=("n_clusters", "median"), k_min=("n_clusters", "min"),
        k_max=("n_clusters", "max"), runs=("ari", "size"), seconds=("seconds", "median")).reset_index()
    return g


def cell(r):
    if r is None or not len(r):
        return "n/a"
    v = r.iloc[0]
    if v.runs > 1 and v.ari_max - v.ari_min >= 0.005:
        return f"{v.ari:.3f} [{v.ari_min:.2f}, {v.ari_max:.2f}]"
    return f"{v.ari:.3f}"


def kcell(r, truth):
    if r is None or not len(r):
        return "n/a"
    v = r.iloc[0]
    k = f"{int(v.k)}" if v.k_min == v.k_max else f"{int(v.k_min)} to {int(v.k_max)}"
    return k


def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    nums = {}
    facts = pd.read_csv(REPO / "docs" / "dataset_facts.csv").set_index("dataset")
    d = load_stability()
    if not len(d):
        print("no stability results yet")
        return 1
    present = [x for x in ORDER if x in set(d.ds)]
    truth_k = {x: int(facts.loc[f"paper_{x}", "label_classes"]) for x in present}

    # ---------- dataset table ----------
    lines = []
    for x in ORDER:
        f = facts.loc[f"paper_{x}"]
        lines.append(f"{NAMES[x]} & {int(f.cells_after_filtering):,} & {int(f.genes_after_filtering):,} & "
                     f"{int(f.label_classes)} & {f.smallest_class_percent:.1f} & {f.percent_zeros:.0f} \\\\")
    (TABLES / "a3_datasets.tex").write_text(
        "\\begin{table}[t]\n\\centering\n\\caption{The ten benchmark datasets, as released by the scMINER authors "
        "and used here without further preprocessing. The smallest class column is that class as a share of all "
        "cells, and the last column is the share of the matrix that is zero.}\n\\label{tab:data}\n\\footnotesize\n"
        "\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{lrrrrr}\n\\toprule\nDataset & Cells & Genes & Types & "
        "Smallest (\\%) & Zeros (\\%) \\\\\n\\midrule\n" + "\n".join(lines) +
        "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n")

    # ---------- headline rule table ----------
    tabs = {r: agg(d, r) for r, _ in RULES}
    have = [(r, lab) for r, lab in RULES if len(tabs[r])]
    body, kbody = [], []
    for x in present:
        cells, kcells = [], []
        for r, _ in have:
            t = tabs[r]
            for m in ("mica", "leiden_pca"):
                sel = t[(t.ds == x) & (t.method == m)] if len(t) else None
                cells.append(cell(sel))
                kcells.append(kcell(sel, truth_k[x]))
        body.append(f"{NAMES[x]} & " + " & ".join(cells) + " \\\\")
        kbody.append(f"{NAMES[x]} ({truth_k[x]}) & " + " & ".join(kcells) + " \\\\")
    ncol = 2 * len(have)
    groups = " & ".join(f"\\multicolumn{{2}}{{c}}{{{lab}}}" for _, lab in have)
    cmid = "".join(f"\\cmidrule(lr){{{2 + 2 * i}-{3 + 2 * i}}}" for i in range(len(have)))
    heads = " & ".join(["MICA & Leiden"] * len(have))
    (TABLES / "a3_rules.tex").write_text(
        "\\begin{table*}[t]\n\\centering\n\\caption{ARI of the partition each label-free rule chose, with no "
        "reference labels used in the choice. Median over three seeds, with the range in brackets where the runs "
        "differed by more than 0.005. Both methods saw the same input, the same candidate range of 2 to 20 "
        "clusters and the same subsamples.}\n\\label{tab:rules}\n\\footnotesize\n\\setlength{\\tabcolsep}{5pt}\n"
        f"\\begin{{tabular}}{{l{'cc' * len(have)}}}\n\\toprule\n & {groups} \\\\\n{cmid}\nDataset & {heads} \\\\\n"
        "\\midrule\n" + "\n".join(body) + "\n\\bottomrule\n\\end{tabular}\n\\end{table*}\n")
    (TABLES / "a3_clusters.tex").write_text(
        "\\begin{table}[t]\n\\centering\n\\caption{Number of clusters each rule chose. The true number of types "
        "is in brackets after the dataset name. A rule can reach a fair ARI with a badly wrong cluster count, "
        "which a user would notice.}\n\\label{tab:clusters}\n\\scriptsize\n\\setlength{\\tabcolsep}{3pt}\n"
        f"\\begin{{tabular}}{{l{'cc' * len(have)}}}\n\\toprule\n & {groups} \\\\\n{cmid}\nDataset (types) & "
        f"{heads} \\\\\n\\midrule\n" + "\n".join(kbody) + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n")

    # ---------- numbers for the text ----------
    for r, _ in RULES:
        t = tabs[r]
        if not len(t):
            continue
        piv = t.pivot_table(index="ds", columns="method", values="ari")
        if {"mica", "leiden_pca"} <= set(piv.columns):
            piv = piv.dropna()
            nums[r] = {"per_dataset": piv.round(4).to_dict(),
                       "mica_mean": round(float(piv.mica.mean()), 4),
                       "leiden_mean": round(float(piv.leiden_pca.mean()), 4),
                       "mica_wins": int((piv.mica > piv.leiden_pca).sum()), "n_datasets": int(len(piv)),
                       "median_gap": round(float((piv.mica - piv.leiden_pca).median()), 4)}
        kp = t.pivot_table(index="ds", columns="method", values="k")
        if len(kp):
            nums.setdefault(r, {})["chosen_k"] = kp.round(1).to_dict()
        nums.setdefault(r, {})["seconds_median"] = t.groupby("method").seconds.median().round(1).to_dict()
    # ---------- mechanism: can the rule tell the settings apart at all? ----------
    from scipy.stats import spearmanr
    curve_all = d[d.k_setting == "stability_curve"].dropna(subset=["stability"])
    TIE = 0.01  # two settings this close are not meaningfully different
    mrows = []
    for (ds_, m), g in curve_all.groupby(["ds", "method"]):
        tied = g[g.stability >= g.stability.max() - TIE]
        rho = spearmanr(g.stability, g.ari).correlation if len(g.dropna(subset=["ari"])) > 4 else np.nan
        mrows.append({"ds": ds_, "method": m, "n_settings": len(g),
                      "share_tied": round(len(tied) / len(g), 3),
                      "k_tied_min": int(tied.n_clusters.min()), "k_tied_max": int(tied.n_clusters.max()),
                      "k_tied_span": int(tied.n_clusters.max() - tied.n_clusters.min()),
                      "spread": round(float(g.stability.max() - g.stability.min()), 3),
                      "spearman_stability_ari": round(float(rho), 3) if rho == rho else None})
    mech = pd.DataFrame(mrows)
    chosen = d[d.k_setting == "stability_finest"].groupby(["ds", "method"]).n_clusters.median().rename("chosen_k")
    mech = mech.merge(chosen.reset_index(), on=["ds", "method"], how="left")
    mech["true_k"] = mech.ds.map(truth_k)
    mech["excess_clusters"] = mech.chosen_k - mech.true_k
    nums["mechanism"] = {
        "tie_window": TIE,
        "mean_share_tied": mech.groupby("method").share_tied.mean().round(3).to_dict(),
        "mean_k_span_among_tied": mech.groupby("method").k_tied_span.mean().round(1).to_dict(),
        "mean_spread": mech.groupby("method").spread.mean().round(3).to_dict(),
        "mean_spearman_stability_ari": mech.groupby("method").spearman_stability_ari.mean().round(3).to_dict(),
        "spearman_spread_vs_excess_clusters": {
            m: round(float(spearmanr(g.spread, g.excess_clusters).correlation), 3)
            for m, g in mech.dropna(subset=["excess_clusters"]).groupby("method") if len(g) > 3},
        "n_datasets": int(mech.ds.nunique()),
        "table": mech.to_dict(orient="records")}
    mech.to_csv(RESULTS / "a3_mechanism.csv", index=False)
    print("mechanism: share of settings tied with the best " + str(nums["mechanism"]["mean_share_tied"]))
    print("mechanism: mean span of cluster counts among tied settings " + str(nums["mechanism"]["mean_k_span_among_tied"]))

    # ---------- does MICA's shipped default suit this benchmark in particular? ----------
    dflt = d[(d.k_setting == "default") & (d.method == "mica")]
    if len(dflt):
        g = dflt.groupby("ds").agg(ari=("ari", "median")).reset_index()
        g["true_k"] = g.ds.map(truth_k)
        near = g[g.true_k.between(3, 5)]
        far = g[~g.true_k.between(3, 5)]
        nums["mica_default"] = {
            "shipped_k": 4, "source": "MICA/mica_mds.py, -nck default",
            "median_true_k_of_benchmark": int(np.median(list(truth_k.values()))),
            "datasets_with_true_k_4": [k for k, v in truth_k.items() if v == 4],
            "mean_ari_true_k_3_to_5": round(float(near.ari.mean()), 3),
            "mean_ari_true_k_outside": round(float(far.ari.mean()), 3),
            "n_near": int(len(near)), "n_far": int(len(far))}
        print(f"MICA default k=4: ARI {nums['mica_default']['mean_ari_true_k_3_to_5']} where the true "
              f"count is 3 to 5, {nums['mica_default']['mean_ari_true_k_outside']} outside")

    # ---------- how sensitive is the choice to the threshold itself? ----------
    # Relaxing the threshold admits more settings, so the finest choice can only move up.
    # Only this direction is meaningful: raising it can empty the passing set entirely,
    # which says nothing about the rule. The shifts are the size of the difference the
    # bootstrap flavour could plausibly make to the threshold.
    curve0 = d[(d.k_setting == "stability_curve") & d.ci_low.notna()]
    def finest_at(g, bump):
        u = g[g.n_clusters >= 2]
        passing = u[u.stability >= u.ci_low.max() + bump]
        if not len(passing):
            return None
        return int(passing.sort_values(["n_clusters", "stability"], ascending=[False, False]).iloc[0].n_clusters)
    srows = []
    for (ds_, m), g in curve0.groupby(["ds", "method"]):
        g = g.groupby("n_clusters", as_index=False).agg({"stability": "median", "ci_low": "median"})
        base = finest_at(g, 0.0)
        shifted = {f"minus_{abs(b):.2f}": finest_at(g, b) for b in (-0.01, -0.02)}
        srows.append({"ds": ds_, "method": m, "chosen_k": base, **shifted,
                      "changed": any(v is not None and v != base for v in shifted.values())})
    sens = pd.DataFrame(srows)
    if len(sens):
        nums["threshold_sensitivity"] = {
            "shifts": [-0.01, -0.02],
            "changed_of_total": f"{int(sens.changed.sum())} of {len(sens)}",
            "by_method": sens.groupby("method").changed.sum().to_dict(),
            "table": sens.to_dict(orient="records")}
        sens.to_csv(RESULTS / "a3_threshold_sensitivity.csv", index=False)
        print(f"threshold sensitivity: the chosen cluster count changes in "
              f"{nums['threshold_sensitivity']['changed_of_total']} cases when the threshold is "
              f"relaxed by at most 0.02")

    nums["true_k"] = truth_k
    (RESULTS / "a3_report_numbers.json").write_text(json.dumps(nums, indent=1, default=str))

    # ---------- mechanism figure ----------
    curve = d[d.k_setting == "stability_curve"].copy()
    if len(curve):
        show = [x for x in ["buettner", "pollen", "zeisel"] if x in set(curve.ds)][:3]
        if show:
            fig, axes = plt.subplots(1, len(show), figsize=(7.0, 2.1), sharey=True)
            axes = np.atleast_1d(axes)
            plt.rcParams.update({"font.size": 7, "font.family": "serif", "pdf.fonttype": 42})
            for ax, x in zip(axes, show):
                for m, col, lab in [("mica", BLUE, "MICA"), ("leiden_pca", ORANGE, "Leiden")]:
                    s = curve[(curve.ds == x) & (curve.method == m)].groupby("n_clusters").agg(
                        stab=("stability", "median"), ari=("ari", "median")).reset_index().sort_values("n_clusters")
                    if not len(s):
                        continue
                    ax.plot(s.n_clusters, s.stab, color=col, lw=1.5, label=f"{lab} stability")
                    ax.plot(s.n_clusters, s.ari, color=col, lw=1.2, ls=(0, (2, 1.5)), label=f"{lab} ARI")
                ax.axvline(truth_k[x], color=MUTED, lw=0.8, ls=":")
                ax.text(truth_k[x], 1.04, "true", fontsize=6, color=MUTED, ha="center")
                ax.set_title(NAMES[x], fontsize=7.5, color=INK)
                ax.set_xlabel("clusters")
                ax.set_ylim(-0.05, 1.1)
                ax.grid(axis="y", color="#e6e5e1", lw=0.5)
                ax.set_axisbelow(True)
                for sp in ("top", "right"):
                    ax.spines[sp].set_visible(False)
            axes[0].set_ylabel("score")
            axes[0].legend(fontsize=5.5, frameon=False, loc="lower left", ncol=2, columnspacing=0.8,
                           handlelength=1.5, labelspacing=0.2)
            fig.tight_layout(pad=0.3)
            fig.savefig(FIGS / "a3_stability_curve.pdf")
            fig.savefig(FIGS / "a3_stability_curve.png", dpi=200)
            print("wrote the stability curve figure")

    for k, v in nums.items():
        if isinstance(v, dict) and "mica_mean" in v:
            print(f"{k:20s} MICA {v['mica_mean']:.3f}  Leiden {v['leiden_mean']:.3f}  "
                  f"MICA wins {v['mica_wins']} of {v['n_datasets']}")
    print(f"wrote tables to {TABLES}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
