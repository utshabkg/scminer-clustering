"""Figures for the Dec 2 slides, from the same saved runs as the paper.

Usage: python scripts/tools/slide_figures.py
Writes docs/slides/pollen_curve.png: MICA and Leiden on Pollen, the confirmatory runs on matched candidates
(seeds 3 to 9, 25 subsamples), the stability curve each reading reads and the ARI the labels give, with
each reading's pick marked. It is the middle panel of the paper's Figure 1, drawn large for a slide.
"""
from glob import glob

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

from scminer_clust.paths import REPO, RESULTS

V1 = RESULTS / "v1"
OUT = REPO / "docs" / "slides"
SEEDS = range(3, 10)
BLUE, ORANGE, INK, MUTED, PAPER = "#2A6F97", "#C96A1B", "#14213D", "#6B7280", "#FAF8F3"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    c = pd.concat(pd.read_csv(f) for f in glob(str(V1 / "curves" / "paper_pollen_seed*.csv")))
    c = c[(c.scorer == "ours") & (c.b == 25) & c.in_mid & c.seed.isin(SEEDS)]
    p = pd.concat(pd.read_csv(f) for f in glob(str(V1 / "picks" / "paper_pollen_seed*.csv")))
    p = p[(p.scorer == "ours") & (p.b == 25) & (p.design == "mid") & p.seed.isin(SEEDS)]
    plt.rcParams.update({"font.size": 15, "font.family": "sans-serif", "axes.edgecolor": MUTED,
                         "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK})
    fig, ax = plt.subplots(figsize=(9.6, 5.4), dpi=200)
    fig.patch.set_facecolor(PAPER)
    ax.set_facecolor(PAPER)
    for m, col, lab in [("mica", BLUE, "MICA"), ("leiden_pca", ORANGE, "Leiden")]:
        s = c[c.method == m].groupby("n_clusters").agg(stab=("stability", "median"), ari=("ari", "median")).reset_index()
        ax.plot(s.n_clusters, s.stab, color=col, lw=3.2, label=f"{lab}: stability (no labels)")
        ax.plot(s.n_clusters, s.ari, color=col, lw=2.4, ls=(0, (3, 2)), label=f"{lab}: ARI (needs labels)")
        for reading, mk in (("stability_finest", "^"), ("stability_most", "v")):
            k = p[(p.method == m) & (p.reading == reading)].n_clusters.median()
            y = float(s.set_index("n_clusters").stab.reindex([k]).interpolate().iloc[0]) if k in set(s.n_clusters) else \
                float(pd.Series(s.stab.values, index=s.n_clusters).sort_index().loc[:k].iloc[-1])
            ax.plot([k], [y], marker=mk, color=col, ms=15, mec=PAPER, mew=1.5, ls="none", zorder=5)
    ax.axvline(11, color=MUTED, lw=1.6, ls=":")
    ax.text(11.3, 0.42, "11 true types", color=MUTED, fontsize=14)
    ax.set_xlabel("number of clusters")
    ax.set_ylabel("score")
    ax.set_ylim(0, 1.06)
    ax.set_xlim(1.5, 20.5)
    ax.set_xticks(range(2, 21, 2))
    ax.grid(axis="y", color="#E3E1DA", lw=0.8)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    handles, labels = ax.get_legend_handles_labels()
    handles += [Line2D([], [], marker="^", color=INK, ls="none", ms=11), Line2D([], [], marker="v", color=INK, ls="none", ms=11)]
    labels += ["chooseR’s reading picks", "plain maximum picks"]
    ax.legend(handles, labels, fontsize=12.5, frameon=False, loc="lower right", ncol=2, columnspacing=1.2)
    fig.tight_layout(pad=0.4)
    fig.savefig(OUT / "pollen_curve.png", facecolor=PAPER, metadata={"Software": None})
    plt.close(fig)
    print("wrote", OUT / "pollen_curve.png")


if __name__ == "__main__":
    main()
