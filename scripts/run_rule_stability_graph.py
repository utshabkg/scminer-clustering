"""The stability rule for MICA's graph mode and for Leiden, on PBMC14K (U-32).

Usage: python scripts/run_rule_stability_graph.py [seed] [b_subsamples]
PBMC14K has 13,605 cells, above the 5,000 cell point where scMINER switches MICA
to graph mode, so the small-data runner cannot be used for it. Each subsample
goes through MICA's real graph pipeline: an 80 nearest neighbor graph under the
mutual information distance, node2vec, then Louvain over the resolution grid.
node2vec lives in the mica-legacy environment, so that step runs as a subprocess
there, through the same scripts/mica_downstream.py we used for the Assignment 2
graph comparison.

As in the small-data runner, the full distance matrix is computed once and
subsampled, which is exact (tests/test_subsample_invariance.py).

This is the most expensive job in the project after the atlas pilot. Budget about
six minutes per subsample for MICA and under a minute for Leiden.
"""
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from scminer_clust import cluster, metrics, mi, mi_knn, stability
from scminer_clust.paths import DATA, PROCESSED, REPO

NAME = "paper_pbmc14k"
K_NEIGHBORS = 79  # MICA's 80 counts the cell itself
RESOLUTIONS = np.round(np.arange(0.1, 4.0 + 1e-9, 0.1), 1)  # MICA's own Louvain grid
LEIDEN_RESOLUTIONS = np.round(np.arange(0.05, 2.0 + 1e-9, 0.05), 2)
WORK = DATA / "runs" / "rules" / "pbmc14k_graph"


def full_distance(X):
    cache = DATA / "runs" / "rules" / f"{NAME}_ge_distance.npy"
    cache.parent.mkdir(parents=True, exist_ok=True)
    if cache.exists():
        return np.load(cache, mmap_mode="r")
    t0 = time.time()
    D = mi.distance_matrix(X, "ge")
    np.save(cache, D)
    print(f"  full MI distance for {len(X)} cells in {time.time() - t0 / 60:.0f} min", flush=True)
    return np.load(cache, mmap_mode="r")


def mica_graph_partitions(D, cells, tag, workers=8):
    """MICA's graph pipeline on a set of cells. Returns a dict of resolution to labels."""
    sub = np.asarray(cells)
    d = np.asarray(D[np.ix_(sub, sub)])
    np.fill_diagonal(d, np.inf)
    k = min(K_NEIGHBORS, len(sub) - 1)
    idx = np.argpartition(d, k, axis=1)[:, :k]
    dst = np.take_along_axis(d, idx, axis=1)
    order = np.argsort(dst, axis=1)
    idx, dst = np.take_along_axis(idx, order, axis=1), np.take_along_axis(dst, order, axis=1)
    gdir = WORK / tag
    if gdir.exists():
        shutil.rmtree(gdir)
    mi_knn.write_mica_edgelist(idx, dst, gdir, [str(c) for c in sub])
    out = WORK / f"{tag}_out"
    if out.exists():
        shutil.rmtree(out)
    cmd = ["conda", "run", "-n", "mica-legacy", "python", "-u", str(REPO / "scripts" / "mica_downstream.py"),
           str(gdir), str(out), str(workers)]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO)
    if r.returncode != 0:
        raise RuntimeError(f"node2vec step failed for {tag}: {r.stderr[-800:]}")
    parts = {}
    for f in sorted(out.glob("clustering_*.txt")):
        s = pd.read_csv(f, sep="\t", index_col=0)["label"]
        parts[float(f.stem.split("_")[-1])] = s.reindex([str(c) for c in sub]).to_numpy().astype(str)
    shutil.rmtree(gdir, ignore_errors=True)
    shutil.rmtree(out, ignore_errors=True)
    return parts


def leiden_partitions(adata, cells, seed=0):
    g, _ = cluster._graph(adata[np.asarray(cells)].copy(), seed=seed)
    out = {}
    for r in LEIDEN_RESOLUTIONS:
        sc.tl.leiden(g, resolution=float(r), random_state=seed, flavor="igraph", n_iterations=2, directed=False)
        out[float(r)] = g.obs["leiden"].to_numpy().astype(str)
    return out


def run(method, maker, n, truth, seed, b):
    t0 = time.time()
    full = maker(np.arange(n), "full")
    settings = [s for s, lab in full.items() if 2 <= len(np.unique(lab)) <= 20]
    print(f"  {method}: {len(settings)} settings give 2 to 20 clusters on the full data", flush=True)
    ref = stability.reference_cells(n, seed=seed)
    co = {s: stability.CoClustering(ref) for s in settings}
    for i, cells in enumerate(stability.subsamples(n, b=b, seed=seed), 1):
        parts = maker(cells, f"sub{i}")
        for s in settings:
            if s in parts:
                co[s].add(cells, parts[s])
        print(f"    subsample {i} of {b} done ({(time.time() - t0) / 60:.0f} min elapsed)", flush=True)
    scores = {s: stability.score_setting(co[s].distance(), np.asarray(full[s])[ref], seed=seed) for s in settings}
    picks = stability.pick(scores)
    rows = []
    for rule, info in picks.items():
        rec = {"dataset": NAME, "method": method, "k_setting": rule, "seed": seed, "genes": "all",
               "mode": "graph" if method == "mica" else "pca", "b_subsamples": b,
               "seconds": round(time.time() - t0, 1), **info}
        rec.update(metrics.score(truth, full[info["chosen_setting"]]))
        rows.append(rec)
        print(f"    {rule}: setting {info['chosen_setting']}, {rec['n_clusters']} clusters, "
              f"ARI {rec['ari']:.3f}", flush=True)
    for s, sc_ in scores.items():
        rows.append({"dataset": NAME, "method": method, "k_setting": "stability_curve", "setting": s,
                     "seed": seed, **sc_, "ari": metrics.score(truth, full[s])["ari"]})
    return rows


def main(seed=0, b=stability.B_SUBSAMPLES):
    seed, b = int(seed), int(b)
    WORK.mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(PROCESSED / f"{NAME}.h5ad")
    X = np.asarray(a.X.todense() if sparse.issparse(a.X) else a.X, dtype=np.float32)
    truth = a.obs["label"].to_numpy()
    n = a.n_obs
    print(f"{NAME}: {n} cells, {X.shape[1]} genes, seed {seed}, {b} subsamples", flush=True)
    D = full_distance(X)
    rows = []
    rows += run("mica", lambda c, t: mica_graph_partitions(D, c, f"s{seed}_{t}"), n, truth, seed, b)
    rows += run("leiden_pca", lambda c, t: leiden_partitions(a, c, seed), n, truth, seed, b)
    d = metrics.write(f"rule_stability_{NAME}_graph_seed{seed}", rows,
                      {"dataset": NAME, "rule": "chooseR-style subsampling stability", "mica_mode": "graph",
                       "b_subsamples": b, "fraction": stability.FRACTION, "seed": seed,
                       "mica_resolutions": RESOLUTIONS.tolist(), "leiden_resolutions": LEIDEN_RESOLUTIONS.tolist()})
    print(f"wrote {d / 'metrics.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:3]))
