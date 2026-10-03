"""RQ1 part one: what happens when nobody tells the method how many clusters to find.

Usage: python scripts/summarize_label_free.py
Reads docs/results_tracker.csv for the ten benchmark sets on the authors' input
and writes results/summary_label_free.csv and .md. For each method and gene
setting: the score under the paper's rule (told the true number), the score when
the method picks by silhouette, the chosen and true cluster counts, and how often
the choice landed on the edge of the search range (2 or 20 clusters).
"""
import sys

import numpy as np
import pandas as pd

from scminer_clust.paths import REPO, RESULTS

ALL10 = ["yan", "goolam", "buettner", "pollen", "chung", "usoskin", "kolodziejczyk", "klein", "zeisel", "pbmc14k"]
METHODS = {"mica": "MICA", "leiden_pca": "Leiden", "kmeans_pca": "k-means", "hierarchical_pca": "hierarchical"}


def main():
    t = pd.read_csv(REPO / "docs" / "results_tracker.csv")
    t = t[t["dataset"].isin([f"paper_{d}" for d in ALL10]) & ~t["experiment"].str.endswith("_ge")]
    rows = []
    for ds in ALL10:
        for g in sorted(t["gene_setting"].dropna().unique()):
            for m, label in METHODS.items():
                sel = t[(t["dataset"] == f"paper_{ds}") & (t["gene_setting"] == g) & (t["method"] == m)]
                ora, sil = sel[sel.k_setting == "paper_protocol"], sel[sel.k_setting == "silhouette"]
                dfl = sel[sel.k_setting.astype(str).str.startswith("default")]
                if not len(ora) or not len(sil):
                    continue
                edge = sil["n_clusters"].isin([2, 20]).mean() if m != "leiden_pca" else np.nan
                rows.append({"dataset": ds, "gene_setting": g, "method": label, "true_k": int(ora["n_true_classes"].iloc[0]),
                             "ari_told_k": round(ora["ari"].mean(), 3), "ari_silhouette": round(sil["ari"].mean(), 3),
                             "loss": round(ora["ari"].mean() - sil["ari"].mean(), 3),
                             "k_chosen": round(sil["n_clusters"].mean(), 1), "share_at_edge": edge,
                             "ari_default": round(dfl["ari"].mean(), 3) if len(dfl) else None,
                             "ami_told_k": round(ora["ami"].mean(), 3), "ami_silhouette": round(sil["ami"].mean(), 3)})
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "summary_label_free.csv", index=False)
    piv = df.pivot_table(index=["dataset", "gene_setting"], columns="method", values="ari_silhouette").round(3)
    loss = df.groupby(["method", "gene_setting"])["loss"].agg(["mean", "median", "max"]).round(3)
    over = df.groupby("method").apply(lambda d: pd.Series({
        "mean ARI told k": d.ari_told_k.mean(), "mean ARI silhouette": d.ari_silhouette.mean(),
        "mean loss": d.loss.mean(), "cases chosen k equals true k": f"{int((d.k_chosen.round() == d.true_k).sum())} of {len(d)}"})).round(3)
    md = ("# Choosing the number of clusters without labels\n\n"
          "Ten benchmark sets on the authors' own input. Told k is the paper's rule. Silhouette means each method "
          "tries k = 2 to 20 (or its resolution sweep) and keeps the best silhouette in its own embedding.\n\n"
          "## Overall\n\n" + over.to_markdown() + "\n\n## Loss from not knowing k, by gene setting\n\n" + loss.to_markdown()
          + "\n\n## ARI with silhouette-chosen k\n\n" + piv.to_markdown() + "\n")
    (RESULTS / "summary_label_free.md").write_text(md)
    print(over.to_string()); print(); print(loss.to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
