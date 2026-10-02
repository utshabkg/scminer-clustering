"""One table for the gene-set question: is it the genes or the method?

Usage: python scripts/summarize_gene_sets.py
Reads docs/results_tracker.csv and writes results/summary_gene_sets.csv and
results/summary_gene_sets.md. Each cell is the mean ARI over seeds, with the
number of seeds behind it. Every number uses the paper's selection rule.
"""
import sys

import pandas as pd

from scminer_clust.paths import REPO, RESULTS

ORDER = ["all genes", "2000 genes (seurat)", "2000 genes (cell_ranger)", "2000 genes (seurat_v3)"]
DATASETS = ["zeisel", "pbmc14k", "usoskin", "kolodziejczyk", "klein"]
PUBLISHED = {"zeisel": 0.744, "pbmc14k": 0.853, "usoskin": 0.847, "kolodziejczyk": 1.0, "klein": 0.98}


def main(metric="ari"):
    t = pd.read_csv(REPO / "docs" / "results_tracker.csv")
    t = t[t["k_setting"] == "paper_protocol"]
    rows = []
    for ds in DATASETS:
        paper = t[(t["dataset"] == f"paper_{ds}") & (t["method"] == "mica") & ~t["experiment"].str.endswith("_ge")]
        row = {"dataset": ds, "published MICA": PUBLISHED[ds],
               "MICA, authors' input": round(paper[metric].mean(), 3) if len(paper) else None,
               "authors' gene count": int(paper["n_genes"].iloc[0]) if len(paper) and pd.notna(paper["n_genes"].iloc[0]) else None}
        ge = t[t["experiment"].isin([f"reproduction_{ds}_ge", f"reproduction_paper_{ds}_ge"])]
        row["MICA graph mode, ours / authors' input"] = (
            " / ".join(f"{ge[ge.experiment == e][metric].mean():.3f}" for e in (f"reproduction_{ds}_ge", f"reproduction_paper_{ds}_ge"))
            if len(ge) else None)
        for method, label in (("mica", "MICA"), ("leiden_pca", "Leiden"), ("kmeans_pca", "k-means"), ("hierarchical_pca", "hierarchical")):
            for g in ORDER:
                sel = t[(t["dataset"] == ds) & (t["method"] == method) & (t["gene_setting"] == g) & ~t["experiment"].str.endswith("_ge")]
                ari = round(sel[metric].mean(), 3) if len(sel) else None
                # a star marks cells where the method missed the true cluster count, so it
                # was scored at the wrong number of clusters (Leiden's grid cannot always hit k)
                missed = len(sel) and (sel["n_clusters"] != sel["n_true_classes"]).any()
                row[f"{label}, {g}"] = (f"{ari:.3f}*" if missed else ari) if ari is not None else None
                row[f"{label}, {g}, seeds"] = len(sel)
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / f"summary_gene_sets_{metric}.csv", index=False)

    show = (["dataset", "published MICA", "MICA, authors' input", "authors' gene count", "MICA graph mode, ours / authors' input"]
            + [f"MICA, {g}" for g in ORDER] + [f"Leiden, {g}" for g in ORDER])
    md = df[show].rename(columns=lambda c: c.replace("2000 genes ", "2k ").replace("all genes", "all"))
    note = (f"\nAll cells are the mean {metric.upper()} over seeds under the paper's selection rule (true number of clusters). "
            "A star means at least one seed could not produce the true number of clusters and was scored at the closest count. "
            "Published values are ARI.\n")
    print(f"=== {metric.upper()}"); print(md.to_string(index=False))
    return f"## {metric.upper()}\n\n" + md.to_markdown(index=False) + "\n" + note


if __name__ == "__main__":
    parts = [main("ari"), main("ami")]
    (RESULTS / "summary_gene_sets.md").write_text("# MICA against Leiden, by gene set\n\n" + "\n".join(parts))
    sys.exit(0)
