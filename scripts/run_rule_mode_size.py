"""Separate MICA's mode from the dataset size, which every earlier run confounds (U-39).

Usage: python scripts/run_rule_mode_size.py [n_pool] [seed] [b_subsamples] [--matched]

The stability results so far compare nine datasets under 3,000 cells, all in MICA's
small-data mode, against one dataset of 13,605 cells in its graph mode. Mode and size
move together, so a difference between them cannot be attributed to either.

This holds the data and the size fixed and changes only the mode. We draw one pool of
n_pool cells from PBMC14K and run the same stability rule on those same cells three
times: MICA's small-data mode, MICA's graph mode, and the Leiden baseline. Comparing
the first two gives the mode effect at fixed size; comparing the second against the
full 13,605 cell graph run gives the size effect at fixed mode.

The two MICA modes differ in the distance as well as the embedding, because scMINER
normalises the mutual information by the square root of the marginal entropies in
small-data mode and by the joint entropy in graph mode. We keep that difference rather
than removing it, since the question is about the two pipelines scMINER actually runs.

The graph-mode distance comes from the cached full matrix by taking a submatrix, which
is exact (tests/test_subsample_invariance.py). The small-data distance uses the other
normalisation, so it is computed on the pool directly.

With --matched every arm uses the matched grids in src/scminer_clust/grids.py and the run
writes results/rule_matched_mode_size_pbmc14k_n<pool>_seed<seed>/.

Writes results/rule_mode_size_pbmc14k_n<pool>_seed<seed>/.
"""
import sys
import time

import anndata as ad
import numpy as np
from scipy import sparse

from scminer_clust import cluster, grids, metrics, mi, stability
from scminer_clust.paths import DATA, PROCESSED

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import run_rule_stability as small
import run_rule_stability_graph as graph

NAME = "paper_pbmc14k"
N_POOL = 2500


def pool_cells(n, n_pool, seed):
    """A fixed random pool of cells, the shared input for every arm."""
    return np.sort(np.random.default_rng(seed).choice(n, size=n_pool, replace=False))


def mds_distance_for_pool(X, tag):
    """Small-data-mode distance on the pool. Cached, since it is the one real cost here."""
    cache = DATA / "runs" / "rules" / f"{tag}_mds_distance.npy"
    cache.parent.mkdir(parents=True, exist_ok=True)
    if cache.exists():
        D = np.load(cache)
        assert D.shape == (len(X), len(X)), f"cached {D.shape} but pool is {len(X)} cells"
        return D
    t0 = time.time()
    D = mi.distance_matrix(X, "mds")
    np.save(cache, D)
    print(f"  small-data MI distance for {len(X)} cells in {(time.time() - t0) / 60:.1f} min", flush=True)
    return D


def run_arm(method, mode, maker, n, truth, seed, b, full=None, settings=None):
    """The stability rule for one arm. Mirrors the other two runners so results are comparable."""
    t0 = time.time()
    if full is None:
        full = maker(np.arange(n), "full")
    if settings is None:
        settings = [s for s, lab in full.items() if 2 <= len(np.unique(lab)) <= 20]
    print(f"  {method} ({mode}): {len(settings)} settings give 2 to 20 clusters", flush=True)
    ref = stability.reference_cells(n, seed=seed)
    co = {s: stability.CoClustering(ref) for s in settings}
    subs = stability.subsamples(n, b=b, seed=seed)
    for i, cells in enumerate(subs, 1):
        parts = maker(cells, f"sub{i}", settings) if grids.matched() else maker(cells, f"sub{i}")
        for s in settings:
            if s in parts:
                co[s].add(cells, parts[s])
        print(f"    subsample {i} of {b} done ({(time.time() - t0) / 60:.1f} min)", flush=True)
    scores = {s: stability.score_setting(co[s].distance(), np.asarray(full[s])[ref], seed=seed)
              for s in settings}
    picks = stability.pick(scores)
    rows = []
    for rule, info in picks.items():
        rec = {"dataset": NAME, "method": method, "mica_mode": mode, "k_setting": rule, "seed": seed,
               "n_cells": n, "b_subsamples": b, "seconds": round(time.time() - t0, 1),
               "n_settings": len(settings), **info}
        rec.update(metrics.score(truth, full[info["chosen_setting"]]))
        rows.append(rec)
        print(f"    {rule}: setting {info['chosen_setting']}, {rec['n_clusters']} clusters, "
              f"ARI {rec['ari']:.3f}", flush=True)
    # The span of cluster counts among settings that tie with the best is the quantity the
    # paper's mechanism argument rests on, so record it per arm rather than recomputing later.
    best = max(v["stability"] for v in scores.values())
    tied = [s for s, v in scores.items() if v["stability"] >= best - 0.01]
    ks = [len(np.unique(full[s])) for s in tied]
    print(f"    tied within 0.01: {len(tied)} of {len(settings)} settings, "
          f"cluster counts {min(ks)} to {max(ks)}", flush=True)
    for s, sc_ in scores.items():
        rows.append({"dataset": NAME, "method": method, "mica_mode": mode,
                     "k_setting": "stability_curve", "setting": s, "seed": seed, "n_cells": n,
                     **sc_, "ari": metrics.score(truth, full[s])["ari"]})
    return rows, {"method": method, "mode": mode, "n_cells": n, "n_settings": len(settings),
                  "n_tied": len(tied), "share_tied": round(len(tied) / len(settings), 3),
                  "k_tied_min": int(min(ks)), "k_tied_max": int(max(ks)),
                  "k_tied_span": int(max(ks) - min(ks))}


def main(n_pool=N_POOL, seed=0, b=stability.B_SUBSAMPLES):
    n_pool, seed, b = int(n_pool), int(seed), int(b)
    assert graph.LEGACY_PYTHON.exists(), f"mica-legacy interpreter missing at {graph.LEGACY_PYTHON}"
    graph.WORK.mkdir(parents=True, exist_ok=True)
    cluster.N_HVG = None  # all genes, as in every other stability run
    a = ad.read_h5ad(PROCESSED / f"{NAME}.h5ad")
    pool = pool_cells(a.n_obs, n_pool, seed)
    sub = a[pool].copy()
    X = np.asarray(sub.X.todense() if sparse.issparse(sub.X) else sub.X, dtype=np.float32)
    truth = sub.obs["label"].to_numpy()
    print(f"{NAME}: pool of {n_pool} of {a.n_obs} cells, {X.shape[1]} genes, "
          f"{len(np.unique(truth))} true types, seed {seed}, {b} subsamples", flush=True)

    tag = f"{NAME}_pool{n_pool}_seed{seed}"
    Dmds = mds_distance_for_pool(X, tag)
    # Graph-mode distance: an exact submatrix of the full cached matrix.
    Dge = np.asarray(np.load(DATA / "runs" / "rules" / f"{NAME}_ge_distance.npy",
                             mmap_mode="r")[np.ix_(pool, pool)])

    rows, summary = [], []
    arms = [
        ("mica", "small", lambda c, t, res=None: small.mica_partitions(Dmds, np.asarray(c), seed)[1]),
        ("mica", "graph", lambda c, t, res=None: graph.mica_graph_partitions(Dge, c, f"ms{seed}_{t}",
                                                                              resolutions=res)),
        ("leiden_pca", "pca", lambda c, t, res=None: small.leiden_partitions(sub, np.asarray(c), seed, res)[1]),
    ]
    fulls, chosen = {}, {}
    if grids.matched():
        # all three arms are offered the same cluster counts, one setting each
        fulls = {mode: mk(np.arange(n_pool), "full") for _, mode, mk in arms}
        common = grids.common_counts(*fulls.values())
        chosen = {mode: list(grids.one_per_count(fulls[mode], common).values()) for mode in fulls}
        print(f"  matched counts: {common} ({len(common)} of 19)", flush=True)
    for method, mode, maker in arms:
        r, s = run_arm(method, mode, maker, n_pool, truth, seed, b, fulls.get(mode), chosen.get(mode))
        rows += r
        summary.append(s)

    d = metrics.write(grids.tag(f"rule_mode_size_pbmc14k_n{n_pool}_seed{seed}"), rows,
                      {"dataset": NAME, "question": "mode against size at fixed data",
                       "grid": "matched" if grids.matched() else "original",
                       "n_pool": n_pool, "pool_seed": seed, "b_subsamples": b,
                       "fraction": stability.FRACTION, "seed": seed,
                       "k_range": list(small.K_RANGE),
                       "mica_resolutions": graph.RESOLUTIONS.tolist(),
                       "leiden_resolutions": small.RESOLUTIONS.tolist()})
    import pandas as pd
    pd.DataFrame(summary).to_csv(d / "tie_summary.csv", index=False)
    print(f"\nwrote {d / 'metrics.csv'}")
    print(pd.DataFrame(summary).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:4]))
