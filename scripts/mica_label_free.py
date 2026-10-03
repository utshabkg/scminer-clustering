"""MICA with the number of clusters chosen without labels (U-12, U-13).

Usage: python scripts/mica_label_free.py <dataset> <mica_output_dir> [tag] [--exp ID]
Small-data mode: k-means on MICA's saved MDS embedding (mi_reduced.pkl) for
k = 2..20, three seeds, keep the best silhouette in that embedding. Graph mode:
the resolution sweep MICA already wrote, scored by silhouette on its saved
node2vec embedding, plus MICA's default resolution (1.822, the nearest grid value
1.8). Same rule as scripts/run_label_free.py for the baselines.
"""
import argparse
import pickle
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

from scminer_clust import metrics, selection
from scminer_clust.paths import PROCESSED


def main(name, out_dir, tag=None, exp=None, seeds=(0, 1, 2)):
    out_dir = Path(out_dir)
    truth = ad.read_h5ad(PROCESSED / f"{name}.h5ad", backed="r").obs["label"]
    m = re.match(r"hvg(\d+)_(.+)", tag or "")
    genes, flavor = (int(m.group(1)), m.group(2)) if m else ("all", "")
    base = {"dataset": name, "method": "mica", "genes": genes, "hvg_flavor": flavor}
    files = sorted(out_dir.glob("clustering_*.txt"))
    ids = pd.read_csv(files[0], sep="\t", index_col=0).index.astype(str)
    assert set(ids) == set(truth.index), "MICA cells differ from the processed file"
    y = truth.reindex(ids).to_numpy()
    rows = []
    if (out_dir / "mi_reduced.pkl").exists() and len(files) == 1:  # small-data mode
        with open(out_dir / "mi_reduced.pkl", "rb") as f:
            emb = pickle.load(f)
        for s in seeds:
            # enough starts to reach the best fit for each k (see scripts/mica_seeded_kmeans.py)
            parts = {k: KMeans(n_clusters=k, n_init=50, random_state=s).fit_predict(emb).astype(str) for k in range(2, selection.K_MAX + 1)}
            key, info = selection.pick_by_silhouette(emb, parts, s)
            r = {**base, "seed": s, "k_setting": "silhouette", "mode": "small-data", **info}
            r.update(metrics.score(y, parts[key])); rows.append(r)
    else:  # graph mode
        emb_file = next(out_dir.glob("knn_node2vec*_emb*.txt"))
        e = pd.read_csv(emb_file, sep=r"\s+", skiprows=1, header=None, index_col=0).sort_index()
        emb = e.to_numpy()
        assert len(emb) == len(ids), "embedding rows differ from the cell count"
        parts = {f.stem.split("_")[-1]: pd.read_csv(f, sep="\t", index_col=0)["label"].astype(str).to_numpy() for f in files}
        key, info = selection.pick_by_silhouette(emb, parts)
        r = {**base, "seed": 0, "k_setting": "silhouette", "mode": "graph", **info}
        r.update(metrics.score(y, parts[key])); rows.append(r)
        if "1.8" in parts:
            r = {**base, "seed": 0, "k_setting": "default_resolution_1.8", "mode": "graph"}
            r.update(metrics.score(y, parts["1.8"])); rows.append(r)
    exp = exp or (f"labelfree_mica_{name}" + (f"_{tag}" if tag else ""))
    d = metrics.write(exp, rows, {"dataset": name, "mica_output": str(out_dir), "tag": tag})
    print(f"{exp}: " + ", ".join(f"{r['k_setting']} k={r['n_clusters']} ARI {r['ari']:.3f}" for r in rows))
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset"); ap.add_argument("out_dir"); ap.add_argument("tag", nargs="?"); ap.add_argument("--exp")
    a = ap.parse_args()
    sys.exit(main(a.dataset, a.out_dir, a.tag or None, a.exp))
