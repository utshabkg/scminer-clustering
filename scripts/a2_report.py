"""Tables and figure for the Assignment 2 report, built straight from results/.

Usage: python scripts/a2_report.py
Writes latex/tables/a2_reproduction.tex, latex/tables/a2_selection.tex,
latex/figures/a2_neighbor_recall.pdf, and results/a2_report_numbers.json (every
number the report quotes, for checking). Nothing here reruns an experiment.
"""
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scminer_clust.paths import REPO, RESULTS

TABLES, FIGS = REPO / "latex" / "tables", REPO / "latex" / "figures"
ALL10 = ["yan", "goolam", "buettner", "pollen", "chung", "usoskin", "kolodziejczyk", "klein", "zeisel", "pbmc14k"]
NAMES = {"yan": "Yan", "goolam": "Goolam", "buettner": "Buettner", "pollen": "Pollen", "chung": "Chung",
         "usoskin": "Usoskin", "kolodziejczyk": "Kolodziejczyk", "klein": "Klein", "zeisel": "Zeisel", "pbmc14k": "PBMC14K"}
GRAPH_MODE = {"pbmc14k": ["reproduction_paper_pbmc14k"] + [f"repeat_paper_pbmc14k_rep{i}" for i in range(1, 6)],
              "zeisel": ["reproduction_paper_zeisel_ge"] + [f"repeat_paper_zeisel_ge_rep{i}" for i in range(1, 6)]}
LABELFREE_MICA = {"pbmc14k": ["labelfree_mica_paper_pbmc14k"] + [f"repeat_labelfree_paper_pbmc14k_rep{i}" for i in range(1, 6)],
                  "zeisel": ["labelfree_mica_paper_zeisel_ge"] + [f"repeat_labelfree_paper_zeisel_rep{i}" for i in range(1, 6)]}
BLUE, ORANGE, INK, MUTED = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"


def read(exp):
    return pd.read_csv(RESULTS / exp / "metrics.csv")


def spread(values):
    v = np.asarray(values, dtype=float)
    return {"median": float(np.median(v)), "min": float(v.min()), "max": float(v.max()), "n": int(len(v)), "values": [round(x, 4) for x in v]}


def fmt(s):
    if s["n"] == 1 or s["max"] - s["min"] < 0.005:
        return f"{s['median']:.3f}"
    return f"{s['median']:.3f} [{s['min']:.2f}, {s['max']:.2f}]"


def k_range(ks):
    ks = sorted(set(int(k) for k in ks))
    if len(ks) == 1:
        return str(ks[0])
    return f"{ks[0]} or {ks[1]}" if len(ks) == 2 else f"{ks[0]} to {ks[-1]}"


def reproduction(nums):
    pub = pd.read_csv(REPO / "docs" / "published_scores.csv").set_index("dataset")
    lines = []
    for ds in ALL10:
        # The paper ran MICA's small-data mode up to 5,000 cells, so only PBMC14K is graph mode in its headline.
        if ds == "pbmc14k":
            ours = spread([read(e).query("method == 'mica'").ari.iloc[0] for e in GRAPH_MODE[ds]])
        else:
            ours = spread([read(f"reproduction_paper_{ds}").query("method == 'mica'").ari.mean()])
        if ds == "zeisel":
            nums["reproduction"]["zeisel_graph_mode"] = spread([read(e).query("method == 'mica'").ari.iloc[0] for e in GRAPH_MODE[ds]])
        leiden = spread(read(f"baselines_paper_{ds}").query("method == 'leiden_pca' and k_setting == 'paper_protocol'").ari)
        p = pub.loc[ds]
        nums["reproduction"][ds] = {"cells": int(p.n_cells), "published_mica": float(p.MICA_ARI), "ours_mica": ours,
                                    "published_scanpy": float(p.Scanpy_ARI), "ours_leiden_all_genes": leiden}
        mark = {"pbmc14k": "$^\\dagger$", "kolodziejczyk": "$^\\ast$", "zeisel": "$^\\ddagger$"}.get(ds, "")
        lines.append(f"{NAMES[ds]} & {int(p.n_cells):,} & {p.MICA_ARI:.3f} & {ours['median']:.3f}{mark} & {p.Scanpy_ARI:.3f} & {leiden['median']:.3f} \\\\")
    body = "\n".join(lines)
    (TABLES / "a2_reproduction.tex").write_text(
        "\\begin{table}[t]\n\\centering\n\\caption{MICA and PCA + Leiden on the authors' released inputs with all genes, scored by ARI under the paper's rule, "
        "where the true number of types is given and the labels pick the setting.}\n\\label{tab:repro}\n\\footnotesize\n\\setlength{\\tabcolsep}{3.5pt}\n"
        "\\begin{tabular}{lrcccc}\n\\toprule\n & & \\multicolumn{2}{c}{MICA} & \\multicolumn{2}{c}{PCA + Leiden} \\\\\n"
        "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\nDataset & Cells & paper & ours & paper & ours \\\\\n\\midrule\n" + body +
        "\n\\bottomrule\n\\end{tabular}\n\\\\[2pt]\n\\parbox{\\columnwidth}{\\scriptsize As in the paper, sets of up to 5,000 cells use MICA's small-data mode, where we give the final $k$-means 100 restarts and report the mean of three seeds, and PBMC14K uses graph mode. $^\\dagger$Median of six runs, 0.850 to 0.857. $^\\ddagger$Graph mode gives 0.755 to 0.781 over six runs. $^\\ast$Not reproduced from the released files (Sec.~\\ref{sec:results}). The paper's PCA + Leiden is Scanpy on 2,000 variable genes. Ours uses all genes, like MICA, and we report the median of three seeds.}\n\\end{table}\n")


def selection(nums):
    rows = {"default": {}, "silhouette": {}, "paper": {}}
    for ds in ["pbmc14k", "zeisel"]:
        lf = pd.concat([read(e).query("method == 'mica'") for e in LABELFREE_MICA[ds]])
        ld = read(f"labelfree_paper_{ds}").query("method == 'leiden_pca'")
        for rule, mica_key, leiden_key in [("default", "default_resolution_1.8", "default_resolution_1"), ("silhouette", "silhouette", "silhouette")]:
            m, l = lf[lf.k_setting == mica_key], ld[ld.k_setting == leiden_key]
            rows[rule][ds] = {"mica": spread(m.ari), "mica_ami": spread(m.ami), "mica_k": k_range(m.n_clusters),
                              "leiden": spread(l.ari), "leiden_ami": spread(l.ami), "leiden_k": k_range(l.n_clusters)}
        mp = spread([read(e).query("method == 'mica'").ari.iloc[0] for e in GRAPH_MODE[ds]])
        mpa = spread([read(e).query("method == 'mica'").ami.iloc[0] for e in GRAPH_MODE[ds]])
        bl = read(f"baselines_paper_{ds}").query("method == 'leiden_pca' and k_setting == 'paper_protocol'")
        rows["paper"][ds] = {"mica": mp, "mica_ami": mpa, "mica_k": "7", "leiden": spread(bl.ari), "leiden_ami": spread(bl.ami), "leiden_k": "7"}
    nums["selection"] = rows

    label = {"default": "Default", "silhouette": "Silhouette", "paper": "Paper's rule"}
    lines = []
    for block, rules in [("No labels", ["default", "silhouette"]), ("Labels pick the setting", ["paper"])]:
        lines.append(f"\\multicolumn{{6}}{{l}}{{\\textit{{{block}}}}} \\\\")
        for rule in rules:
            for i, (meth, key) in enumerate([("MICA", "mica"), ("Leiden", "leiden")]):
                c = []
                for ds in ["pbmc14k", "zeisel"]:
                    c += [fmt(rows[rule][ds][key]), rows[rule][ds][key + "_k"]]
                lines.append(f"{label[rule] if i == 0 else ''} & {meth} & " + " & ".join(c) + " \\\\")
        if block == "No labels":
            lines.append("\\midrule")
    (TABLES / "a2_selection.tex").write_text(
        "\\begin{table}[t]\n\\centering\n\\caption{ARI of the chosen clustering on the authors' PBMC14K and Zeisel inputs (all genes, "
        "7 types each), as median [min, max] over six MICA graph-mode runs or three Leiden seeds. A single value means all runs agree "
        "within 0.005. $k$ is the number of clusters chosen. In the top block the labels are used only for scoring.}\n"
        "\\label{tab:selection}\n\\scriptsize\n\\setlength{\\tabcolsep}{1.7pt}\n\\begin{tabular}{llcccc}\n\\toprule\n"
        " & & \\multicolumn{2}{c}{PBMC14K} & \\multicolumn{2}{c}{Zeisel} \\\\\n\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\n"
        "Rule & Method & ARI & $k$ & ARI & $k$ \\\\\n\\midrule\n" + "\n".join(lines) + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n")


def figure(nums):
    big = read("rq2_large_paper_pbmc14k")
    al = read("rq2_allen_100k")
    q = big[big.stage == "query sweep"]
    corr = q[q.candidates == "correlation"].sort_values("c")
    pd_ = q[q.candidates == "pca + detected"].sort_values("c")
    full = big[(big.stage == "all cells, round 0") & (big.candidates == "correlation")].set_index("c")
    hn_p = q[q.candidates.str.startswith("MICA HNSW")].iloc[0]
    qa = al[al.stage == "query sweep"]
    ca = qa[qa.candidates == "correlation"].sort_values("c")
    fa = al[al.stage == "all cells, round 0"].set_index("c")
    hn_a = qa[qa.candidates.str.startswith("MICA HNSW")].iloc[0]
    bf_p = q[q.candidates.str.startswith("brute")].iloc[0]
    bf_a = qa[qa.candidates.str.startswith("brute")].iloc[0]
    nums["neighbors"] = {
        "pbmc14k": {"correlation": dict(zip(corr.c.astype(int), corr.recall_at_k.round(4))),
                    "pca_plus_detected": dict(zip(pd_.mean_scored.round(0).astype(int), pd_.recall_at_k.round(4))),
                    "correlation_minutes_all_cells": {int(c): round(s / 60, 1) for c, s in full.seconds_all_cells.items()},
                    "hnsw_recall": round(float(hn_p.recall_at_k), 4), "hnsw_minutes": round(float(hn_p.seconds_all_cells) / 60, 1),
                    "brute_force_minutes_est": round(float(bf_p.seconds_all_cells) / 60, 1)},
        "allen_100k": {"correlation": dict(zip(ca.c.astype(int), ca.recall_at_k.round(4))),
                       "correlation_hours_all_cells": {int(c): round(s / 3600, 2) for c, s in fa.seconds_all_cells.items()},
                       "hnsw_recall": round(float(hn_a.recall_at_k), 4), "hnsw_hours": round(float(hn_a.seconds_all_cells) / 3600, 2),
                       "brute_force_hours_est": round(float(bf_a.seconds_all_cells) / 3600, 1)},
    }
    plt.rcParams.update({"font.size": 7, "font.family": "serif", "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                         "ytick.major.width": 0.6, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                         "axes.labelcolor": INK, "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(3.5, 1.9), sharey=True)
    lab = dict(fontsize=6.3, color=INK, textcoords="offset points")
    ax = axes[0]
    ax.axhline(100 * hn_p.recall_at_k, color=MUTED, lw=0.8, ls=(0, (3, 2)), label="MICA HNSW")
    ax.text(102, 100 * hn_p.recall_at_k + 1.0, "MICA HNSW, 15 min", color=MUTED, fontsize=6.3, va="bottom")
    ax.plot(corr.c, 100 * corr.recall_at_k, color=BLUE, lw=1.4, marker="o", ms=3.2, mec="#fcfcfb", mew=0.6, label="Correlation")
    ax.plot(pd_.mean_scored, 100 * pd_.recall_at_k, color=ORANGE, lw=1.4, marker="s", ms=3.0, mec="#fcfcfb", mew=0.6, label="PCA + detected")
    ax.annotate(f"{full.seconds_all_cells[800] / 60:.1f} min", (800, 100 * full.recall_at_k[800]), xytext=(3, -11), ha="left", **lab)
    ax.annotate("PCA + detected", (pd_.mean_scored.iloc[-1], 100 * pd_.recall_at_k.iloc[-1]), xytext=(4, -2), ha="left", va="center", **lab)
    ax.set_title("PBMC14K, 13,605 cells", fontsize=7, color=INK, pad=3)
    ax.set_ylabel("True neighbors found (%)")
    ax.legend(loc="lower right", fontsize=6.3, frameon=False, handlelength=1.6, borderaxespad=0.1, labelspacing=0.25)
    ax = axes[1]
    ax.axhline(100 * hn_a.recall_at_k, color=MUTED, lw=0.8, ls=(0, (3, 2)))
    ax.text(205, 100 * hn_a.recall_at_k + 1.2, "MICA HNSW, 5.2 h", color=MUTED, fontsize=6.3, va="bottom")
    ax.plot(ca.c, 100 * ca.recall_at_k, color=BLUE, lw=1.4, marker="o", ms=3.2, mec="#fcfcfb", mew=0.6)
    ax.annotate(f"{fa.seconds_all_cells[1600] / 3600:.1f} h", (1600, 100 * fa.recall_at_k[1600]), xytext=(5, -10), ha="left", **lab)
    ax.annotate(f"{fa.seconds_all_cells[3200] / 3600:.1f} h", (3200, 100 * fa.recall_at_k[3200]), xytext=(-5, 2), ha="right", **lab)
    ax.set_title("Allen subset, 100,000 cells", fontsize=7, color=INK, pad=3)
    for ax in axes:
        ax.set_xscale("log", base=2)
        ax.set_ylim(30, 106)
        ax.grid(axis="y", color="#e6e5e1", lw=0.5)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(length=2, pad=1.5)
    axes[0].set_xticks([100, 200, 400, 800, 1600]); axes[0].set_xticklabels(["100", "200", "400", "800", "1.6k"])
    axes[1].set_xticks([200, 400, 800, 1600, 3200]); axes[1].set_xticklabels(["200", "400", "800", "1.6k", "3.2k"])
    axes[0].minorticks_off(); axes[1].minorticks_off()
    fig.supxlabel("Candidates scored exactly per cell (log scale)", fontsize=7, y=0.02, color=INK)
    fig.tight_layout(pad=0.3, w_pad=0.6, rect=(0, 0.04, 1, 1))
    fig.savefig(FIGS / "a2_neighbor_recall.pdf")
    fig.savefig(FIGS / "a2_neighbor_recall.png", dpi=300)


def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    nums = {"reproduction": {}}
    reproduction(nums)
    selection(nums)
    figure(nums)
    (RESULTS / "a2_report_numbers.json").write_text(json.dumps(nums, indent=1, default=str))
    print(json.dumps(nums["selection"], indent=1, default=str)[:3000])
    return 0


if __name__ == "__main__":
    sys.exit(main())
