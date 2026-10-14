"""Build and run notebooks/assignment2_preliminary.ipynb, which rebuilds the Assignment 2 tables and figure.

Usage: python scripts/make_a2_notebook.py
The notebook reads only saved results (results/*/metrics.csv), so it runs in
under a minute. The experiments themselves are run by the scripts it lists.
"""
import sys

import nbformat
from nbclient import NotebookClient

from scminer_clust.paths import REPO

CELLS = [
    ("md", """# Assignment 2: preliminary results

This notebook rebuilds every table and the figure in our Assignment 2 report from the saved result files. It does not rerun experiments. Each number comes from `results/<experiment>/metrics.csv`, written by the scripts below.

| Result | Script that produced it |
|---|---|
| MICA reproduction on the authors' inputs | `scripts/run_mica_job.sh`, `scripts/mica_seeded_kmeans.py`, `scripts/score_mica.py` |
| Leiden, k-means and hierarchical baselines | `scripts/run_baselines.py` |
| Label-free choice (defaults, silhouette) | `scripts/run_label_free.py`, `scripts/mica_label_free.py` |
| Neighbor search on PBMC14K and Allen 100K | `scripts/run_rq2_large.py`, `scripts/run_rq2_allen.py` |
| MICA's later steps on each graph | `scripts/mica_downstream.py`, `scripts/score_downstream.py` |

The data live outside git (see `README.md`). The environment is `envs/scminer.yml`."""),
    ("code", """import json, sys
from pathlib import Path
import pandas as pd
from IPython.display import Image, display
from scminer_clust.paths import REPO, RESULTS
sys.path.insert(0, str(REPO / "scripts"))
import a2_report
a2_report.main()   # writes latex/tables/a2_*.tex, latex/figures/a2_neighbor_recall.pdf and results/a2_report_numbers.json
nums = json.loads((RESULTS / "a2_report_numbers.json").read_text())"""),
    ("md", "## Table I. Reproduction on the authors' inputs (ARI, the paper's rule)"),
    ("code", """t1 = pd.DataFrame({ds: {"cells": r["cells"], "MICA, paper": r["published_mica"], "MICA, ours": round(r["ours_mica"]["median"], 3),
                         "Scanpy, paper": r["published_scanpy"], "Leiden all genes, ours": round(r["ours_leiden_all_genes"]["median"], 3)}
                    for ds, r in nums["reproduction"].items() if "cells" in r}).T
print("Zeisel in graph mode, six runs:", nums["reproduction"]["zeisel_graph_mode"]["values"])
t1["within 0.05"] = (t1["MICA, ours"] - t1["MICA, paper"]).abs() <= 0.05
t1"""),
    ("md", "## Table II. Choosing the clustering without and with labels (PBMC14K and Zeisel, all genes)\n\nMedian [min, max] ARI over six MICA graph-mode runs or three Leiden seeds, and the number of clusters chosen."),
    ("code", """rows = []
for rule, name in [("default", "No labels: default"), ("silhouette", "No labels: silhouette"), ("paper", "Labels pick: paper's rule")]:
    for meth in ["mica", "leiden"]:
        r = {"rule": name, "method": meth}
        for ds in ["pbmc14k", "zeisel"]:
            s = nums["selection"][rule][ds]
            r[f"{ds} ARI"] = f"{s[meth]['median']:.3f} [{s[meth]['min']:.3f}, {s[meth]['max']:.3f}]"
            r[f"{ds} AMI"] = round(s[meth + '_ami']["median"], 3)
            r[f"{ds} k"] = s[meth + "_k"]
        rows.append(r)
pd.DataFrame(rows)"""),
    ("md", "## Figure 1. Neighbor search against brute force"),
    ("code", """display(Image(filename=str(REPO / "latex" / "figures" / "a2_neighbor_recall.png"), width=700))
pd.DataFrame(nums["neighbors"]).T"""),
    ("md", "## MICA's own later steps (node2vec and Louvain) on each graph\n\nThe same code runs on our fast graph and on MICA's own HNSW graph, so only the graph differs."),
    ("code", """runs = []
for f in sorted(RESULTS.glob("rq2_downstream_*/metrics.csv")):
    d = pd.read_csv(f)
    runs.append({"run": f.parent.name, "dataset": d.dataset[0], "graph": d.graph[0], "ARI": round(d.ari[0], 3), "AMI": round(d.ami[0], 3), "clusters": d.n_clusters[0]})
pd.DataFrame(runs)"""),
]


def main():
    nb = nbformat.v4.new_notebook()
    nb.cells = [nbformat.v4.new_markdown_cell(src) if kind == "md" else nbformat.v4.new_code_cell(src) for kind, src in CELLS]
    nb.metadata["kernelspec"] = {"name": "scminer", "display_name": "scminer", "language": "python"}
    out = REPO / "notebooks"
    out.mkdir(exist_ok=True)
    NotebookClient(nb, timeout=600, kernel_name="scminer", resources={"metadata": {"path": str(out)}}).execute()
    nbformat.write(nb, out / "assignment2_preliminary.ipynb")
    print(f"wrote {out / 'assignment2_preliminary.ipynb'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
