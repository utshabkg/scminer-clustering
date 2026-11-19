"""The Allen whole mouse brain atlas on every cell (task U-46).

Usage: python scripts/large_scale/allen.py <stage> [--scope full|sample100k|ctxsp] [--workers N]

Run the stages in this order. Each one is its own process (scripts/large_scale/run.sh wraps it in
scripts/tools/peak_mem.sh), so the time and peak memory of every stage are measured separately.
  manifest   the annotated 10Xv3 cells, in the order data.load_allen reads them
  scan       pass 1 over the raw counts: genes detected per gene, per-gene mean and variance of
             the counts, total counts and genes per cell
  hvg        pass 2: the 2,000 seurat_v3 genes, computed as scanpy computes them but streamed
  normalize  pass 3: preprocess.normalize on each block of cells, keeping the 2,000 genes
  pca        50 principal components with cluster.pca's scaling, from the exact covariance
  kmeans     cluster.kmeans_pca with k = 34 classes, seeds 0, 1 and 2
  leiden     cluster._graph's 15-neighbour graph, then Leiden over a few resolutions
  assign     MICA's 100,000-cell sample run, every other cell given the majority cluster of its
             15 nearest sample cells in the full-data principal components
  collect    scores, dimensions, time and memory, into results/large_scale/

Scopes: full is every annotated cell (2,341,350). sample100k is the 100,000 cells of the pilot
(data/processed/allen_100k.h5ad), used to check this pipeline against the pilot's own numbers.
ctxsp is one region, for a quick dry run.

Shared code in src/scminer_clust is imported and used unchanged.
"""
import argparse
import json
import os
import resource
import sys
import time
from multiprocessing import get_context
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from scminer_clust import cluster, metrics, preprocess
from scminer_clust.data import ALLEN_REGIONS
from scminer_clust.paths import PROCESSED, RAW, REPO, RESULTS

N_HVG = 2000
N_PCS = 50
K = 34              # the atlas's classes, the label level every method is set to
SEEDS = (0, 1, 2)
BLOCK = 50_000      # rows read from a region file at a time
LEIDEN_ROUND1 = (0.05, 0.10, 0.17, 0.30, 0.50, 1.00)  # 0.17 was the pilot's pick, 1.0 scanpy's default
MICA_SAMPLE = REPO / "results" / "reproduction" / "mica_allen_100k"
PILOT_RUNS = REPO / "data" / "runs" / "baselines_allen_100k_hvg2000_seurat_v3"


# ---------------------------------------------------------------- bookkeeping
class Scope:
    def __init__(self, name):
        self.name = name
        suffix = "" if name == "full" else f"_{name}"
        self.runs = REPO / "data" / "runs" / f"large_scale{suffix}"
        self.results = RESULTS / f"large_scale{suffix}"
        for d in (self.runs, self.results):
            d.mkdir(parents=True, exist_ok=True)

    def path(self, name):
        return self.runs / name


def record(scope, stage, extra=None):
    """Append this stage's time and memory to stages.jsonl.

    peak_rss_gb_self is this process, peak_rss_gb_child the largest single worker. The sum over
    the process tree is measured from outside by peak_mem.sh, and counts pages that forked workers
    share more than once, so it is an upper bound.
    """
    row = {"stage": stage, "scope": scope.name, "wall_s": round(time.time() - T0, 1),
           "peak_rss_gb_self": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20, 2),
           "peak_rss_gb_child": round(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 2**20, 2),
           "threads": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS")},
           "load_avg": os.getloadavg(), "finished": time.strftime("%Y-%m-%d %H:%M:%S")}
    row.update(extra or {})
    with open(scope.results / "stages.jsonl", "a") as f:
        f.write(json.dumps(row) + "\n")
    (scope.path(f"{stage}.done")).write_text(json.dumps(row))
    print(json.dumps(row), flush=True)


def save_npz(path, **arrays):
    tmp = path.with_name(path.name + ".tmp.npz")
    np.savez(tmp, **arrays)
    tmp.replace(path)


def load_manifest(scope):
    return pd.read_csv(scope.path("manifest.csv.gz"), index_col="cell_label")


def region_blocks(region, rows, block=BLOCK):
    """Yield (positions in rows, CSR of those rows) for one region file, read block by block.

    rows are file row numbers in increasing order. Each block of the file is read once, whole, and
    only the wanted rows are kept, so reading cost does not depend on how many rows are wanted.
    """
    path = RAW / "allen" / f"WMB-10Xv3-{region}-raw.h5ad"
    with h5py.File(path, "r") as f:
        x = f["X"]
        indptr = x["indptr"][:].astype(np.int64)
        n_genes = int(x.attrs["shape"][1])
        for start in range(0, len(indptr) - 1, block):
            stop = min(start + block, len(indptr) - 1)
            lo_i, hi_i = np.searchsorted(rows, [start, stop])
            if hi_i == lo_i:
                continue
            lo, hi = indptr[start], indptr[stop]
            m = sp.csr_matrix((x["data"][lo:hi], x["indices"][lo:hi], indptr[start:stop + 1] - lo),
                              shape=(stop - start, n_genes))
            yield np.arange(lo_i, hi_i), m[rows[lo_i:hi_i] - start]


def _pool(workers):
    return get_context("fork").Pool(workers)


# ---------------------------------------------------------------- S0 manifest
def stage_manifest(scope, args):
    meta = pd.read_csv(RAW / "allen" / "cell_metadata_with_cluster_annotation.csv",
                       usecols=["cell_label", "feature_matrix_label", "library_method", "class",
                                "subclass", "supertype", "cluster"], index_col="cell_label")
    meta = meta[meta["library_method"] == "10Xv3"]
    sample = set(_sample_cells())
    parts = []
    for r in ALLEN_REGIONS:
        want = meta.index[meta["feature_matrix_label"] == f"WMB-10Xv3-{r}"]
        with h5py.File(RAW / "allen" / f"WMB-10Xv3-{r}-raw.h5ad", "r") as f:
            names = pd.Index(f["obs"][f["obs"].attrs["_index"]][:].astype(str))
        rows = np.flatnonzero(names.isin(want))  # file order, as in data.load_allen
        assert len(rows) == len(want), f"{r}: {len(want) - len(rows)} annotated cells not in the file"
        df = meta.loc[names[rows], ["class", "subclass", "supertype", "cluster"]].copy()
        df.insert(0, "row", rows)
        df.insert(0, "region", r)
        df["in_sample"] = df.index.isin(sample)
        parts.append(df)
    man = pd.concat(parts)
    man.index.name = "cell_label"
    assert man["in_sample"].sum() == len(sample), "every pilot cell must be an annotated 10Xv3 cell"
    if scope.name == "sample100k":
        man = man[man["in_sample"]]
    elif scope.name == "ctxsp":
        man = man[man["region"] == "CTXsp"]
    man.to_csv(scope.path("manifest.csv.gz"))
    n_rows = 0
    for r in ALLEN_REGIONS:
        with h5py.File(RAW / "allen" / f"WMB-10Xv3-{r}-raw.h5ad", "r") as f:
            n_rows += int(f["X"].attrs["shape"][0])
    record(scope, "manifest", {"n_cells": int(len(man)), "n_rows_in_files": n_rows,
                               "n_classes": int(man["class"].nunique()),
                               "n_subclasses": int(man["subclass"].nunique()),
                               "n_in_sample": int(man["in_sample"].sum())})


def _sample_cells():
    with h5py.File(PROCESSED / "allen_100k.h5ad", "r") as f:
        return f["obs"][f["obs"].attrs["_index"]][:].astype(str)


# ---------------------------------------------------------------- S1 scan
def _scan_region(job):
    region, rows = job
    n_det = sums = sumsq = None
    n_umi, n_gene = np.zeros(len(rows)), np.zeros(len(rows), dtype=np.int64)
    for pos, m in region_blocks(region, rows):
        g = m.shape[1]
        if n_det is None:
            n_det, sums, sumsq = np.zeros(g, np.int64), np.zeros(g), np.zeros(g)
        pos_data = m.data > 0
        n_det += np.bincount(m.indices[pos_data], minlength=g)
        sums += np.bincount(m.indices, weights=m.data.astype(np.float64), minlength=g)
        # scanpy squares float32 counts in float32 and sums them in float64 (_get_mean_var)
        sumsq += np.bincount(m.indices, weights=np.square(m.data).astype(np.float64), minlength=g)
        n_umi[pos] = np.asarray(m.sum(axis=1)).ravel()
        n_gene[pos] = np.diff(m.indptr) if not (~pos_data).any() else np.asarray((m > 0).sum(axis=1)).ravel()
    return region, n_det, sums, sumsq, n_umi, n_gene


def stage_scan(scope, args):
    man = load_manifest(scope)
    jobs = [(r, man.loc[man["region"] == r, "row"].to_numpy()) for r in ALLEN_REGIONS
            if (man["region"] == r).any()]
    jobs.sort(key=lambda j: -len(j[1]))
    with _pool(args.workers) as p:
        out = {r[0]: r for r in p.imap_unordered(_scan_region, jobs)}
    n = len(man)
    n_det = sum(out[r][1] for r in out)
    sums = sum(out[r][2] for r in out)
    sumsq = sum(out[r][3] for r in out)
    per_cell = {r: (out[r][4], out[r][5]) for r in out}
    n_umi = np.concatenate([per_cell[r][0] for r in ALLEN_REGIONS if r in per_cell])
    n_gene = np.concatenate([per_cell[r][1] for r in ALLEN_REGIONS if r in per_cell])
    with h5py.File(RAW / "allen" / "WMB-10Xv3-CTXsp-raw.h5ad", "r") as f:
        ids = f["var"][f["var"].attrs["_index"]][:].astype(str)
        symbols = f["var"]["gene_symbol"]
        symbols = symbols[:].astype(str) if isinstance(symbols, h5py.Dataset) else \
            np.asarray(symbols["categories"][:].astype(str))[symbols["codes"][:]]
    mean = sums / n
    var = (sumsq / n - mean ** 2) * n / (n - 1)
    genes = pd.DataFrame({"gene_symbol": symbols, "n_cells_detected": n_det, "mean": mean, "var": var},
                         index=pd.Index(ids, name="gene_identifier"))
    genes["kept"] = genes["n_cells_detected"] >= 3  # sc.pp.filter_genes(min_cells=3), preprocess.py
    genes.to_csv(scope.path("genes.csv.gz"))
    save_npz(scope.path("cells_scan.npz"), n_umi=n_umi, n_gene=n_gene)
    dims = {"n_cells": int(n), "n_genes_raw": int(len(genes)), "n_genes_kept": int(genes["kept"].sum()),
            "median_umi_per_cell": float(np.median(n_umi)), "median_genes_per_cell": float(np.median(n_gene)),
            "n_nonzero": int(n_gene.sum()), "share_zero": float(1 - n_gene.sum() / (n * len(genes)))}
    (scope.results / "dims.json").write_text(json.dumps(dims, indent=2))
    record(scope, "scan", dims)


# ---------------------------------------------------------------- S2 hvg
def _clip_region(job):
    region, rows, kept_pos, clip = job
    s = np.zeros(len(clip))
    s2 = np.zeros(len(clip))
    for _, m in region_blocks(region, rows):
        g = kept_pos[m.indices]
        ok = g >= 0
        # scanpy's _sum_and_sum_squares_clipped: min(float64(x), clip), summed with its square
        e = np.minimum(m.data[ok].astype(np.float64), clip[g[ok]])
        s += np.bincount(g[ok], weights=e, minlength=len(clip))
        s2 += np.bincount(g[ok], weights=e * e, minlength=len(clip))
    return s, s2


def seurat_v3_from_sums(mean, var, n, clipped_sum=None, clipped_sumsq=None, span=0.3):
    """The parts of scanpy's seurat_v3 (scanpy 1.11.5, _highly_variable_genes.py l.98-175), one batch.

    Called twice: first without the clipped sums, to get the clipping values from the loess fit of
    log10 variance on log10 mean, then with them, to get the normalized variance that ranks genes.
    """
    from skmisc.loess import loess

    not_const = var > 0
    estimat_var = np.zeros(len(mean))
    model = loess(np.log10(mean[not_const]), np.log10(var[not_const]), span=span, degree=2)
    model.fit()
    estimat_var[not_const] = model.outputs.fitted_values
    reg_std = np.sqrt(10 ** estimat_var)
    clip_val = reg_std * np.sqrt(n) + mean
    if clipped_sum is None:
        return clip_val
    return (1 / ((n - 1) * np.square(reg_std))) * (n * np.square(mean) + clipped_sumsq - 2 * clipped_sum * mean)


def top_genes(norm_var, n_top):
    """Rank as scanpy does for one batch: argsort twice, keep ranks below n_top, in gene order."""
    ranks = np.argsort(np.argsort(-norm_var))
    return ranks < n_top, ranks


def stage_hvg(scope, args):
    man = load_manifest(scope)
    genes = pd.read_csv(scope.path("genes.csv.gz"), index_col=0)
    kept = genes[genes["kept"]]
    n = len(man)
    clip = seurat_v3_from_sums(kept["mean"].to_numpy(), kept["var"].to_numpy(), n)
    kept_pos = np.full(len(genes), -1, dtype=np.int64)
    kept_pos[np.flatnonzero(genes["kept"].to_numpy())] = np.arange(len(kept))
    jobs = [(r, man.loc[man["region"] == r, "row"].to_numpy(), kept_pos, clip) for r in ALLEN_REGIONS
            if (man["region"] == r).any()]
    jobs.sort(key=lambda j: -len(j[1]))
    with _pool(args.workers) as p:
        parts = p.map(_clip_region, jobs)
    s = sum(x[0] for x in parts)
    s2 = sum(x[1] for x in parts)
    norm = seurat_v3_from_sums(kept["mean"].to_numpy(), kept["var"].to_numpy(), n, s, s2)
    hv, ranks = top_genes(norm, N_HVG)
    out = kept[["gene_symbol", "mean", "var"]].copy()
    out["variances_norm"], out["rank"], out["highly_variable"] = norm, ranks, hv
    out.to_csv(scope.path("hvg.csv.gz"))
    record(scope, "hvg", {"n_hvg": int(hv.sum()), "n_genes_ranked": int(len(kept))})


# ---------------------------------------------------------------- S3 normalize
def _normalize_region(job):
    import anndata as ad

    region, rows, kept_idx, hvg_in_kept = job
    blocks, totals = [], np.zeros(len(rows))
    for pos, m in region_blocks(region, rows):
        a = ad.AnnData(X=m[:, kept_idx].astype(np.float32))
        totals[pos] = np.asarray(a.X.sum(axis=1)).ravel()
        a = preprocess.normalize(a)  # CPM over the kept genes, then log2(x + 1), row by row
        blocks.append(sp.csr_matrix(a.X[:, hvg_in_kept]))
    return region, sp.vstack(blocks).tocsr(), totals


def stage_normalize(scope, args):
    man = load_manifest(scope)
    genes = pd.read_csv(scope.path("genes.csv.gz"), index_col=0)
    hvg = pd.read_csv(scope.path("hvg.csv.gz"), index_col=0)
    kept_idx = np.flatnonzero(genes["kept"].to_numpy())
    hvg_in_kept = np.flatnonzero(hvg["highly_variable"].to_numpy())
    jobs = [(r, man.loc[man["region"] == r, "row"].to_numpy(), kept_idx, hvg_in_kept) for r in ALLEN_REGIONS
            if (man["region"] == r).any()]
    jobs.sort(key=lambda j: -len(j[1]))
    (scope.runs / "hvg").mkdir(exist_ok=True)
    with _pool(args.workers) as p:
        for region, m, totals in p.imap_unordered(_normalize_region, jobs):
            save_npz(scope.runs / "hvg" / f"{region}.npz", data=m.data, indices=m.indices,
                     indptr=m.indptr.astype(np.int64), shape=np.array(m.shape), totals=totals)
    record(scope, "normalize", {"n_hvg": int(len(hvg_in_kept))})


def load_hvg_matrix(scope):
    man = load_manifest(scope)
    mats = []
    for r in ALLEN_REGIONS:
        if not (man["region"] == r).any():
            continue
        z = np.load(scope.runs / "hvg" / f"{r}.npz")
        mats.append(sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])))
    return sp.vstack(mats, format="csr"), man


# ---------------------------------------------------------------- S4 pca
def cov_pca(x, n_comps=N_PCS, block=25_000):
    """cluster.pca's result without making the scaled matrix whole.

    cluster.pca runs sc.pp.scale(max_value=10), which centres each gene, divides by its standard
    deviation and clips at +-10 (scanpy 1.11.5 _scale.py), then sklearn's PCA, which centres again
    and keeps the top components, signed so each component's largest entry is positive. Here the
    scaled rows are made one block at a time, the covariance is summed exactly, and its top
    eigenvectors are the same components.
    """
    from scanpy.preprocessing._utils import _get_mean_var

    mean, var = _get_mean_var(x)
    std = np.sqrt(var)
    std[std == 0] = 1
    n, g = x.shape
    gram, total = np.zeros((g, g)), np.zeros(g)

    def scaled(i):
        z = (x[i:i + block].toarray().astype(np.float64) - mean) / std
        return np.clip(z, -10, 10)

    for i in range(0, n, block):
        z = scaled(i)
        gram += z.T @ z
        total += z.sum(axis=0)
    zbar = total / n
    cov = (gram - n * np.outer(zbar, zbar)) / (n - 1)
    evals, evecs = np.linalg.eigh(cov)
    order = np.argsort(evals)[::-1][:n_comps]
    v = evecs[:, order].T                                    # components, one per row
    v *= np.sign(v[np.arange(n_comps), np.abs(v).argmax(axis=1)])[:, None]
    emb = np.empty((n, n_comps), dtype=np.float32)
    for i in range(0, n, block):
        emb[i:i + block] = (scaled(i) - zbar) @ v.T
    return emb, evals[order], float(evals.sum())


def stage_pca(scope, args):
    x, man = load_hvg_matrix(scope)
    emb, evals, total_var = cov_pca(x)
    np.save(scope.path("pca.npy"), emb)
    save_npz(scope.path("pca_variance.npz"), eigenvalues=evals, total_variance=total_var)
    record(scope, "pca", {"n_cells": int(x.shape[0]), "n_genes": int(x.shape[1]), "n_pcs": N_PCS,
                          "variance_explained": float(evals.sum() / total_var)})


def embedding_adata(scope):
    """An AnnData holding only the cells and the principal components, under the key cluster.pca uses."""
    import anndata as ad

    man = load_manifest(scope)
    emb = np.load(scope.path("pca.npy"))
    a = ad.AnnData(obs=pd.DataFrame(index=man.index.astype(str)))
    for s in SEEDS:
        a.obsm[f"X_pca_{N_PCS}_{s}_all_"] = emb  # computed once. arpack's start vector only adds rounding noise
    return a, man, emb


# ---------------------------------------------------------------- S5 kmeans
def stage_kmeans(scope, args):
    a, man, _ = embedding_adata(scope)
    cluster.N_HVG = None
    (scope.runs / "partitions").mkdir(exist_ok=True)
    timing = {}
    for s in SEEDS:
        f = scope.runs / "partitions" / f"kmeans_k{K}_seed{s}.npy"
        if f.exists():
            continue
        t = time.time()
        labels, _ = cluster.kmeans_pca(a, K, N_PCS, seed=s, n_init=10)
        timing[f"kmeans_seed{s}_s"] = round(time.time() - t, 1)
        np.save(f, labels.astype(np.int32))
    record(scope, "kmeans", timing)


# ---------------------------------------------------------------- S6 and S7 leiden
_GRAPH = None


def _leiden_job(job):
    import scanpy as sc

    res, seed = job
    t = time.time()
    sc.tl.leiden(_GRAPH, resolution=res, random_state=seed, flavor="igraph", n_iterations=2, directed=False)
    labels = _GRAPH.obs["leiden"].to_numpy().astype(np.int32)
    return res, seed, labels, round(time.time() - t, 1)


def stage_leiden(scope, args):
    import anndata as ad
    import scanpy as sc

    global _GRAPH
    sc.settings.n_jobs = args.threads
    graph_file = scope.path("graph.h5ad")
    if graph_file.exists():
        _GRAPH = ad.read_h5ad(graph_file)
        graph_s = json.loads(scope.path("graph.json").read_text())["seconds"]
    else:
        a, man, _ = embedding_adata(scope)
        t = time.time()
        _GRAPH, _ = cluster._graph(a, N_PCS, 15, 0)  # sc.pp.neighbors(n_neighbors=15, use_rep="X_pca")
        graph_s = round(time.time() - t, 1)
        _GRAPH.write_h5ad(graph_file)
        scope.path("graph.json").write_text(json.dumps({"seconds": graph_s, "n_jobs": args.threads}))
    man = load_manifest(scope)
    truth = man["class"].to_numpy()
    target = int(man["class"].nunique())  # 34 on the full atlas, fewer in a one-region dry run
    out = scope.runs / "partitions"
    out.mkdir(exist_ok=True)
    search = scope.results / "leiden_search.csv"
    done = pd.read_csv(search) if search.exists() else pd.DataFrame(columns=["resolution", "seed", "n_clusters", "seconds"])

    def run(jobs):
        nonlocal done
        jobs = [j for j in jobs if not (np.isclose(done["resolution"].astype(float), j[0])
                                        & (done["seed"].astype(int) == j[1])).any()]
        if not jobs:
            return
        with _pool(min(args.workers, len(jobs))) as p:
            for res, seed, labels, secs in p.imap_unordered(_leiden_job, jobs):
                np.save(out / f"leiden_res{res:.2f}_seed{seed}.npy", labels)
                done = pd.concat([done, pd.DataFrame([{"resolution": res, "seed": seed,
                                                       "n_clusters": int(len(np.unique(labels))),
                                                       "seconds": secs}])], ignore_index=True)
                done.to_csv(search, index=False)
                print(f"leiden res {res} seed {seed}: {len(np.unique(labels))} clusters, {secs} s", flush=True)

    run([(r, 0) for r in LEIDEN_ROUND1])
    grid = np.round(np.arange(0.01, 2.0 + 1e-9, 0.01), 2)  # the paper's grid, cluster.PAPER_RESOLUTIONS
    for _ in range(4):  # refine around the resolutions that bracket the number of classes
        d0 = done[done["seed"] == 0].sort_values("resolution")
        if (d0["n_clusters"] == target).any():
            r = d0.loc[d0["n_clusters"] == target, "resolution"].iloc[0]
            run([(round(r + dr, 2), 0) for dr in (-0.01, 0.01)])
            break
        below = d0[d0["n_clusters"] < target]["resolution"].max()
        above = d0[d0["n_clusters"] > target]["resolution"].min()
        lo = below if pd.notna(below) else 0.0
        hi = above if pd.notna(above) else 2.0
        inside = grid[(grid > lo + 1e-9) & (grid < hi - 1e-9)]
        if len(inside) == 0:
            break
        run([(float(r), 0) for r in inside[np.unique(np.linspace(0, len(inside) - 1, min(6, len(inside))).round().astype(int))]])
    d0 = done[done["seed"] == 0]
    parts = {float(r): np.load(out / f"leiden_res{r:.2f}_seed0.npy") for r in d0["resolution"]}
    res, info = cluster.pick_paper_protocol(parts, truth)
    run([(res, s) for s in SEEDS if s != 0] + [(1.0, s) for s in SEEDS if s != 0])
    info.update({"graph_seconds": graph_s, "n_evaluated": int(len(d0)), "picked_resolution": res})
    (scope.results / "leiden_pick.json").write_text(json.dumps(info, indent=2))
    record(scope, "leiden", info)


# ---------------------------------------------------------------- S8 assign
def vote(neigh_labels):
    """Majority label of each row's neighbours, nearest first. A tie goes to the tied label met first."""
    out = np.empty(len(neigh_labels), dtype=neigh_labels.dtype)
    for i, row in enumerate(neigh_labels):
        vals, first, counts = np.unique(row, return_index=True, return_counts=True)
        best = counts == counts.max()
        out[i] = vals[best][np.argmin(first[best])]
    return out


def knn(ref, query, k=15, threads=16):
    from sklearn.neighbors import NearestNeighbors

    nn = NearestNeighbors(n_neighbors=k, algorithm="brute", n_jobs=threads).fit(ref)
    idx = np.empty((len(query), k), dtype=np.int32)
    for i in range(0, len(query), 20_000):
        idx[i:i + 20_000] = nn.kneighbors(query[i:i + 20_000], return_distance=False)
    return idx


def sample_partitions(sample_cells, truth):
    """The partitions made on the 100,000-cell sample: MICA's and the pilot's k-means and Leiden."""
    parts, info = {}, {}
    mica = {}
    for f in sorted(MICA_SAMPLE.glob("clustering_UMAP_euclidean_20_*.txt")):
        res = float(f.stem.rsplit("_", 1)[1])
        d = pd.read_csv(f, sep="\t", index_col="ID")["label"]
        mica[res] = d.reindex(sample_cells).to_numpy()
    assert all(not pd.isna(v).any() for v in mica.values()), "MICA output must cover every sample cell"
    res, info["mica"] = cluster.pick_paper_protocol(mica, truth)  # the pilot's rule
    parts["mica"] = mica[res].astype(np.int32)
    parts["mica_default"] = mica[1.8].astype(np.int32)  # MICA's shipped graph-mode default
    for m in ("kmeans", "leiden"):
        d = pd.read_csv(PILOT_RUNS / f"{m}_pca_seed0.csv", index_col="cell_label")["label"]
        parts[m] = d.reindex(sample_cells).to_numpy().astype(np.int32)
    return parts, info


def stage_assign(scope, args):
    from sklearn.metrics import adjusted_rand_score
    from sklearn.model_selection import KFold

    man = load_manifest(scope)
    emb = np.load(scope.path("pca.npy"))
    is_s = man["in_sample"].to_numpy()
    sample_cells = man.index[is_s]
    parts, info = sample_partitions(sample_cells, man.loc[is_s, "class"].to_numpy())
    t = time.time()
    idx = knn(emb[is_s], emb[~is_s], threads=args.threads) if (~is_s).any() else np.empty((0, 15), np.int32)
    knn_s = round(time.time() - t, 1)
    np.save(scope.path("knn_idx.npy"), idx)
    (scope.runs / "partitions").mkdir(exist_ok=True)
    vote_s = {}
    for name, lab in parts.items():
        t = time.time()
        full = np.empty(len(man), dtype=np.int32)
        full[is_s] = lab
        if len(idx):
            full[~is_s] = vote(lab[idx])
        vote_s[name] = round(time.time() - t, 1)
        np.save(scope.runs / "partitions" / f"assigned_{name}.npy", full)
        np.save(scope.runs / "partitions" / f"sample_{name}.npy", lab)
    # The assignment step checked without labels: hold out a tenth of the sample, assign it from the
    # other nine tenths, and compare with the cluster each held-out cell actually got.
    rows = []
    es = emb[is_s]
    for fold, (tr, te) in enumerate(KFold(10, shuffle=True, random_state=0).split(es)):
        nb = knn(es[tr], es[te], threads=args.threads)
        for name, lab in parts.items():
            pred = vote(lab[tr][nb])
            rows.append({"fold": fold, "partition": name, "accuracy": float((pred == lab[te]).mean()),
                         "ari": float(adjusted_rand_score(lab[te], pred))})
    pd.DataFrame(rows).to_csv(scope.results / "assignment_holdout.csv", index=False)
    info.update({"knn_seconds": knn_s, "vote_seconds": vote_s, "n_reference": int(is_s.sum()),
                 "n_assigned": int((~is_s).sum()), "k": 15, "space": "full-data principal components"})
    (scope.results / "assign_info.json").write_text(json.dumps(info, indent=2, default=float))
    record(scope, "assign", {"knn_seconds": knn_s})


# ---------------------------------------------------------------- S9 collect
def stage_collect(scope, args):
    from sklearn.metrics import adjusted_rand_score

    man = load_manifest(scope)
    truth = man["class"].to_numpy()
    is_s = man["in_sample"].to_numpy()
    part = scope.runs / "partitions"
    rows = []

    def add(method, data, setting, rule, labels, clustered):
        r = {"method": method, "data": data, "setting": setting, "selection": rule,
             "n_cells_clustered": clustered, "n_cells_scored": int(len(labels))}
        sc_all = metrics.score(truth if len(labels) == len(truth) else truth[is_s], labels)
        r.update({k: sc_all[k] for k in ("n_clusters", "ari", "ami", "nmi")})
        if len(labels) == len(truth) and is_s.any():
            s = metrics.score(truth[is_s], labels[is_s])
            r.update({"ari_sample_cells": s["ari"], "ami_sample_cells": s["ami"]})
            if (~is_s).any() and data != "full data":
                r["ari_assigned_cells"] = metrics.score(truth[~is_s], labels[~is_s])["ari"]
        rows.append(r)

    km = {s: np.load(part / f"kmeans_k{K}_seed{s}.npy") for s in SEEDS if (part / f"kmeans_k{K}_seed{s}.npy").exists()}
    if km:
        add("PCA + k-means", "full data", f"k = {K}", "k set to the number of classes", km[0], len(man))
    pick = json.loads((scope.results / "leiden_pick.json").read_text()) if (scope.results / "leiden_pick.json").exists() else None
    if pick:
        r = pick["picked_resolution"]
        add("PCA + Leiden", "full data", f"resolution {r:.2f}", "closest count to the classes, then ARI (paper rule)",
            np.load(part / f"leiden_res{r:.2f}_seed0.npy"), len(man))
        add("PCA + Leiden", "full data", "resolution 1.00", "scanpy default, no labels",
            np.load(part / "leiden_res1.00_seed0.npy"), len(man))
    for name, method, setting, rule in [("kmeans", "PCA + k-means", f"k = {K}", "k set to the number of classes"),
                                        ("leiden", "PCA + Leiden", "pilot resolution", "paper rule on the sample"),
                                        ("mica", "MICA (graph mode)", "paper rule", "paper rule on the sample"),
                                        ("mica_default", "MICA (graph mode)", "resolution 1.8", "MICA's default, no labels")]:
        f = part / f"sample_{name}.npy"
        if f.exists() and scope.name == "full":
            add(method, "100k sample", setting, rule, np.load(f), int(is_s.sum()))
            add(method, "sample + assigned", setting, rule, np.load(part / f"assigned_{name}.npy"), int(is_s.sum()))
    df = pd.DataFrame(rows)
    # label-free reproducibility: agreement between seeds of the same method
    agree = []
    for s, t in [(0, 1), (0, 2), (1, 2)]:
        if s in km and t in km:
            agree.append({"method": "PCA + k-means", "setting": f"k = {K}", "seeds": f"{s}-{t}",
                          "ari": float(adjusted_rand_score(km[s], km[t]))})
        if pick:
            for r in (pick["picked_resolution"], 1.0):
                fs, ft = part / f"leiden_res{r:.2f}_seed{s}.npy", part / f"leiden_res{r:.2f}_seed{t}.npy"
                if fs.exists() and ft.exists():
                    agree.append({"method": "PCA + Leiden", "setting": f"resolution {r:.2f}", "seeds": f"{s}-{t}",
                                  "ari": float(adjusted_rand_score(np.load(fs), np.load(ft)))})
    pd.DataFrame(agree).to_csv(scope.results / "seed_agreement.csv", index=False)
    if km:
        metrics.per_type_recall(truth, km[0]).to_csv(scope.results / "per_type_recall_kmeans.csv", index=False)
    stages = [json.loads(l) for l in (scope.results / "stages.jsonl").read_text().splitlines()]
    pd.DataFrame(stages).drop(columns=["threads"]).to_csv(scope.results / "stages.csv", index=False)
    params = {"scope": scope.name, "n_hvg": N_HVG, "n_pcs": N_PCS, "k": K, "seeds": list(SEEDS),
              "hvg_flavor": "seurat_v3, streamed, all cells", "normalization": "preprocess.normalize: CPM then log2(x + 1)",
              "gene_filter": "detected in at least 3 cells", "cell_filter": "none, the atlas authors' QC",
              "labels": "class (34) from cell_metadata_with_cluster_annotation.csv",
              "mica_sample": "uniform random 100,000 cells, seed 0 (data.load_allen)",
              "assignment": "15 nearest sample cells in the full-data principal components, majority vote"}
    df.to_csv(scope.results / "metrics.csv", index=False)
    (scope.results / "params.json").write_text(json.dumps(params, indent=2))
    print(df.to_string(), flush=True)
    record(scope, "collect")


STAGES = {"manifest": stage_manifest, "scan": stage_scan, "hvg": stage_hvg, "normalize": stage_normalize,
          "pca": stage_pca, "kmeans": stage_kmeans, "leiden": stage_leiden, "assign": stage_assign,
          "collect": stage_collect}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=list(STAGES))
    ap.add_argument("--scope", default="full", choices=["full", "sample100k", "ctxsp"])
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--threads", type=int, default=16)
    a = ap.parse_args()
    T0 = time.time()
    STAGES[a.stage](Scope(a.scope), a)
