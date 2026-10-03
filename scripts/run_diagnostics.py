"""RQ1 diagnostics: co-detection share and depth sensitivity of MICA (U-14).

Usage: python scripts/run_diagnostics.py zeisel
Uses our processed data and our MICA distance, which matches the authors' code
to 1e-6 (tests/test_mi_matches_mica.py), with MICA's small-data pipeline (MDS
then k-means at the true number of clusters, three seeds). Writes
results/rq1_diagnostics_<name>/metrics.csv. Distance matrices are cached in
data/runs/diagnostics/ so reruns are cheap.
"""
import sys
import time

import anndata as ad
import numpy as np
from sklearn.metrics import adjusted_mutual_info_score as ami, adjusted_rand_score as ari

from scminer_clust import cluster, diagnostics as dg, metrics, mi
from scminer_clust.paths import DATA, PROCESSED

FRACTIONS = (0.75, 0.5, 0.25)


def cached_distance(path, compute):
    if path.exists():
        return np.load(path)
    D = compute()
    np.save(path, D)
    return D


def main(name):
    a = ad.read_h5ad(PROCESSED / f"{name}.h5ad")
    y = a.obs["label"].to_numpy()
    k = len(np.unique(y))
    X = np.asarray(a.X.todense() if hasattr(a.X, "todense") else a.X, dtype=np.float32)
    cache = DATA / "runs" / "diagnostics"; cache.mkdir(parents=True, exist_ok=True)
    rows = []

    def score(analysis, method, labels_list, ref_list=None, **extra):
        r = {"dataset": name, "analysis": analysis, "method": method, "n_cells": len(y), "n_genes": X.shape[1],
             "ari": float(np.mean([ari(y, l) for l in labels_list])), "ami": float(np.mean([ami(y, l) for l in labels_list]))}
        if ref_list is not None:
            r["ari_vs_full"] = float(np.mean([ari(rf, l) for rf, l in zip(ref_list, labels_list)]))
        r.update(extra); rows.append(r)
        print(f"  {analysis:16s} {method:14s} ARI {r['ari']:.3f}" + (f"  vs full {r['ari_vs_full']:.3f}" if "ari_vs_full" in r else "")
              + "".join(f"  {k_}={v:.3f}" for k_, v in extra.items() if isinstance(v, float)))

    t0 = time.time()
    D = cached_distance(cache / f"{name}_full.npy", lambda: mi.distance_matrix(X, "mds"))
    _, full = dg.mica_small_data(D, k)
    print(f"{name}: {len(y)} cells, {X.shape[1]} genes, distance in {time.time() - t0:.0f}s")
    score("full", "mica", full)

    # co-detection: random pairs, split by whether the two cells share a label
    rng = np.random.default_rng(0)
    pairs = rng.integers(0, len(y), size=(20000, 2)); pairs = pairs[pairs[:, 0] != pairs[:, 1]]
    share, bin0 = dg.codetection_share(X, pairs)
    same = y[pairs[:, 0]] == y[pairs[:, 1]]
    Ddet = cached_distance(cache / f"{name}_detection.npy", lambda: dg.detection_distance(X))
    _, det = dg.mica_small_data(Ddet, k)
    score("detection_only", "mica", det, full, codetection_share_median=float(np.nanmedian(share)),
          share_same_type=float(np.nanmedian(share[same])), share_diff_type=float(np.nanmedian(share[~same])),
          detected_values_in_bin0=bin0, knn_overlap=dg.knn_overlap(D, Ddet))

    # depth: thin the counts, renormalize the same genes, recluster with both methods
    counts = a.layers["counts"]
    leiden_full, _, _ = cluster.leiden_paper_protocol(a, y, seed=0)
    score("full", "leiden_pca", [leiden_full])
    for p in FRACTIONS:
        Xp = dg.log_cpm(dg.thin(counts, p, seed=0))
        Dp = cached_distance(cache / f"{name}_thin{p}.npy", lambda: mi.distance_matrix(Xp, "mds"))
        _, lab = dg.mica_small_data(Dp, k)
        score(f"thin_{p}", "mica", lab, full, knn_overlap=dg.knn_overlap(D, Dp))
        b = ad.AnnData(X=Xp, obs=a.obs[["label"]].copy())
        b.uns["log1p"] = {"base": 2}
        lp, _, _ = cluster.leiden_paper_protocol(b, y, seed=0)
        score(f"thin_{p}", "leiden_pca", [lp], [leiden_full])

    d = metrics.write(f"rq1_diagnostics_{name}", rows, {"dataset": name, "fractions": FRACTIONS, "pairs": len(pairs),
                                                        "pipeline": "our MICA distance, MDS 19 dims, k-means at true k, seeds 0-2"})
    print(f"wrote {d / 'metrics.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
