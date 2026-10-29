"""Pick each method's clustering by subsampling stability, with no labels (U-32).

Usage: python scripts/run_rule_stability.py <dataset> [n_hvg] [hvg_flavor] [seed] [--matched]
Runs the chooseR-style rule in src/scminer_clust/stability.py for MICA and for
PCA plus Leiden on the same input, over a matched range of cluster counts (2 to
20). Both methods are fully rerun on each of 25 random 80 percent subsamples.
The labels are read only at the end, to score the partition the rule chose.

MICA's distance between two cells depends only on those two cells, because each
cell is binned over its own range. Subsampling therefore cannot change it, which
is checked in tests/test_subsample_invariance.py. So the full matrix is computed
once per dataset and cached, and each subsample takes a submatrix of it. Every
other step, the MDS embedding included, is redone per subsample.

With --matched the baseline uses the widened Leiden grid in src/scminer_clust/grids.py, so
both methods are offered the same cluster counts, and the run writes
results/rule_matched_stability_<exp>/. MICA's k grid already reaches every count.

Writes results/rule_stability_<exp>/. Besides the chosen settings and each setting's
stability, it records every setting's full-data cluster count and ARI ("full_sweep" rows),
which the transfer comparison needs.
"""
import sys
import time

import anndata as ad
import numpy as np
import scanpy as sc
from scipy import sparse

from scminer_clust import cluster, diagnostics, grids, metrics, mi, stability
from scminer_clust.paths import DATA, PROCESSED

K_RANGE = range(2, 21)
N_INIT = 50  # as in diagnostics.KMEANS_RESTARTS, enough to reach MICA's own optimum
RESOLUTIONS = grids.leiden()
SEED = 0


def dense(a):
    return np.asarray(a.X.todense() if sparse.issparse(a.X) else a.X, dtype=np.float32)


def mica_distance(name, X):
    """Full small-data-mode MI distance, cached. Subsampling it is exact (see the docstring)."""
    cache = DATA / "runs" / "rules" / f"{name}_mds_distance.npy"
    cache.parent.mkdir(parents=True, exist_ok=True)
    if cache.exists():
        D = np.load(cache)
        assert D.shape == (len(X), len(X)), f"cached distance is {D.shape}, data is {len(X)} cells"
        return D
    t0 = time.time()
    D = mi.distance_matrix(X, "mds")
    np.save(cache, D)
    print(f"  MI distance for {len(X)} cells in {time.time() - t0:.0f} s")
    return D


def mica_partitions(D, cells, seed=SEED):
    """MICA's small-data pipeline on a set of cells: MDS of their distances, then k-means."""
    emb = diagnostics.classical_mds(D[np.ix_(cells, cells)])
    from sklearn.cluster import KMeans
    return emb, {k: KMeans(n_clusters=k, n_init=N_INIT, random_state=seed).fit_predict(emb) for k in K_RANGE}


def leiden_partitions(adata, cells, seed=SEED, resolutions=None):
    """The Scanpy pipeline on a set of cells: principal components, neighbors, Leiden sweep.

    Each resolution is a separate seeded Leiden run on the same graph, so sweeping only some
    of them gives the same partitions as sweeping all.
    """
    sub = adata[cells].copy()
    g, emb = cluster._graph(sub, seed=seed)
    out = {}
    for r in (RESOLUTIONS if resolutions is None else resolutions):
        sc.tl.leiden(g, resolution=float(r), random_state=seed, flavor="igraph", n_iterations=2, directed=False)
        out[float(r)] = g.obs["leiden"].to_numpy().astype(str)
    return emb, out


def run_method(name, adata, D, make_partitions, n, seed=SEED, full=None, settings=None):
    """The stability rule for one method: full-data partitions, then B subsample reruns.

    The matched comparison passes the full-data partitions and the candidate settings in,
    since they are chosen across both methods.
    """
    t0 = time.time()
    allc = np.arange(n)
    if full is None:
        _, full = make_partitions(allc)
    if settings is None:
        # Only settings whose full-data partition has 2 to 20 clusters, a cap shared by both
        # methods. A cap is not a match, which the matched comparison fixes.
        settings = [s for s, lab in full.items() if 2 <= len(np.unique(lab)) <= 20]
    ref = stability.reference_cells(n, seed=seed)
    co = {s: stability.CoClustering(ref) for s in settings}
    subs = stability.subsamples(n, seed=seed)
    for i, cells in enumerate(subs, 1):
        _, parts = make_partitions(cells, settings) if grids.matched() else make_partitions(cells)
        for s in settings:
            co[s].add(cells, parts[s])
        print(f"    subsample {i} of {len(subs)}", flush=True)
    scores = {}
    for s in settings:
        lab_ref = np.asarray(full[s])[ref]
        scores[s] = stability.score_setting(co[s].distance(), lab_ref, seed=seed)
    picks = stability.pick(scores)
    secs, nc = round(time.time() - t0, 1), len(settings) * (len(subs) + 1)
    for v in picks.values():
        v["seconds"], v["n_clusterings"] = secs, nc
    return picks, full, scores


def main(name, n_hvg=None, flavor="seurat", seed=SEED):
    cluster.N_HVG = int(n_hvg) if n_hvg else None
    cluster.HVG_FLAVOR = flavor
    seed = int(seed)
    tag = f"{name}" + (f"_hvg{n_hvg}_{flavor}" if n_hvg else "")
    exp_tag = tag + (f"_seed{seed}" if seed != SEED else "")
    adata = ad.read_h5ad(PROCESSED / f"{name}.h5ad")
    truth = adata.obs["label"].to_numpy()
    n = adata.n_obs
    genes = cluster.select_hvg(adata, int(n_hvg), flavor) if n_hvg else None
    X = dense(adata[:, genes] if genes is not None else adata)
    print(f"{name}: {n} cells, {X.shape[1]} genes, {len(np.unique(truth))} true types")
    D = mica_distance(tag, X)

    rows = []
    makers = [("mica", lambda c, res=None: mica_partitions(D, c, seed)),
              ("leiden_pca", lambda c, res=None: leiden_partitions(adata, c, seed, res))]
    fulls, chosen = {}, {}
    if grids.matched():
        # Both methods are offered the same cluster counts, one setting each.
        fulls = {m: mk(np.arange(n))[1] for m, mk in makers}
        common = grids.common_counts(*fulls.values())
        chosen = {m: list(grids.one_per_count(fulls[m], common).values()) for m in fulls}
        print(f"  matched counts: {common} ({len(common)} of 19)", flush=True)
    for method, maker in makers:
        print(f"  {method}: stability rule")
        picks, full, scores = run_method(name, adata, D, maker, n, seed, fulls.get(method), chosen.get(method))
        for rule, info in picks.items():
            rec = {"dataset": name, "method": method, "k_setting": rule, "seed": seed,
                   "genes": int(n_hvg) if n_hvg else "all", "hvg_flavor": flavor if n_hvg else "",
                   "b_subsamples": stability.B_SUBSAMPLES, "fraction": stability.FRACTION, **info}
            rec.update(metrics.score(truth, full[info["chosen_setting"]]))
            rows.append(rec)
            print(f"    {rule}: setting {info['chosen_setting']}, {rec['n_clusters']} clusters, "
                  f"ARI {rec['ari']:.3f}, {info['seconds']:.0f} s, {info['n_clusterings']} clusterings", flush=True)
        # every setting's full-data count and ARI, in or out of the 2 to 20 range, for the
        # transfer comparison, which tunes a setting on the other datasets' labels
        for s, lab in full.items():
            rows.append({"dataset": name, "method": method, "k_setting": "full_sweep", "setting": s,
                         "seed": seed, "genes": int(n_hvg) if n_hvg else "all",
                         "n_clusters": len(np.unique(lab)), "ari": metrics.score(truth, lab)["ari"]})
        # every setting's stability score, so the choice can be inspected
        for s, sc_ in scores.items():
            rows.append({"dataset": name, "method": method, "k_setting": "stability_curve", "setting": s,
                         "seed": seed, "genes": int(n_hvg) if n_hvg else "all", **sc_,
                         "ari": metrics.score(truth, full[s])["ari"]})
    d = metrics.write(grids.tag(f"rule_stability_{exp_tag}"), rows,
                      {"dataset": name, "rule": "chooseR-style subsampling stability", "n_hvg": n_hvg,
                       "grid": "matched" if grids.matched() else "original",
                       "b_subsamples": stability.B_SUBSAMPLES, "fraction": stability.FRACTION,
                       "k_range": [2, 20], "resolutions": RESOLUTIONS.tolist(), "seed": seed})
    print(f"wrote {d / 'metrics.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None,
                  sys.argv[3] if len(sys.argv) > 3 else "seurat",
                  sys.argv[4] if len(sys.argv) > 4 else SEED))
