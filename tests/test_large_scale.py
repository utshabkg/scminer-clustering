"""The Allen pipeline (scripts/large_scale/allen.py) gives what the in-memory code gives.

Each piece that the pipeline streams or rewrites for size is compared with the code it replaces,
on synthetic counts small enough to hold whole.
"""
import importlib.util
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
import scanpy as sc
import scipy.sparse as sp

from scminer_clust import cluster, preprocess

spec = importlib.util.spec_from_file_location(
    "allen", Path(__file__).resolve().parents[1] / "scripts" / "large_scale" / "allen.py")
allen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(allen)


@pytest.fixture(scope="module")
def counts():
    rng = np.random.default_rng(0)
    n, g = 3000, 600
    mu = rng.gamma(0.3, 2.0, size=g)
    groups = rng.integers(0, 4, size=n)
    shift = rng.lognormal(0, 1, size=(4, g))
    lam = mu * shift[groups] * rng.lognormal(0, 0.3, size=(n, 1))
    x = rng.poisson(lam).astype(np.float32)
    x[:, :5] = 0
    x[:2, :5] = 1  # five genes seen in only two cells, so the gene filter drops them
    return sp.csr_matrix(x)


def _chunks(m, size=700):
    return [m[i:i + size] for i in range(0, m.shape[0], size)]


def test_streamed_seurat_v3_matches_scanpy(counts):
    keep = np.asarray((counts > 0).sum(axis=0)).ravel() >= 3
    x = counts[:, keep]
    a = ad.AnnData(X=x.copy())
    a.layers["counts"] = x.copy()
    sc.pp.highly_variable_genes(a, n_top_genes=200, flavor="seurat_v3", layer="counts")
    n = x.shape[0]
    s = sum(np.asarray(c.sum(axis=0)).ravel().astype(np.float64) for c in _chunks(x))
    s2 = sum(np.bincount(c.indices, weights=np.square(c.data).astype(np.float64), minlength=x.shape[1])
             for c in _chunks(x))
    mean = s / n
    var = (s2 / n - mean ** 2) * n / (n - 1)
    np.testing.assert_allclose(mean, a.var["means"], rtol=1e-12)
    np.testing.assert_allclose(var, a.var["variances"], rtol=1e-9)
    clip = allen.seurat_v3_from_sums(mean, var, n)
    cs = np.zeros(x.shape[1])
    cs2 = np.zeros(x.shape[1])
    for c in _chunks(x):
        e = np.minimum(c.data.astype(np.float64), clip[c.indices])
        cs += np.bincount(c.indices, weights=e, minlength=x.shape[1])
        cs2 += np.bincount(c.indices, weights=e * e, minlength=x.shape[1])
    norm = allen.seurat_v3_from_sums(mean, var, n, cs, cs2)
    np.testing.assert_allclose(norm, a.var["variances_norm"], rtol=1e-9)
    hv, _ = allen.top_genes(norm, 200)
    assert (hv == a.var["highly_variable"].to_numpy()).all()


def test_blockwise_normalize_matches_whole(counts):
    whole = preprocess.normalize(ad.AnnData(X=counts.copy())).X
    parts = sp.vstack([preprocess.normalize(ad.AnnData(X=c.copy())).X for c in _chunks(counts)]).tocsr()
    assert (abs(whole - parts)).max() == 0


def test_cov_pca_matches_cluster_pca(counts):
    x = preprocess.normalize(ad.AnnData(X=counts.copy())).X.tocsr()
    a = ad.AnnData(X=x.copy())
    cluster.N_HVG = None
    ref = cluster.pca(a, n_comps=20, seed=0, n_hvg=None)
    emb, evals, total = allen.cov_pca(x, n_comps=20, block=500)
    for j in range(15):  # the trailing components of a small matrix are close in variance
        r = np.corrcoef(ref[:, j], emb[:, j])[0, 1]
        assert r > 0.9999, (j, r)


def test_vote_breaks_ties_by_nearest():
    neigh = np.array([[3, 1, 1, 3, 2],    # 1 and 3 tie, 3 is nearer
                      [2, 2, 5, 5, 5],    # 5 wins outright
                      [7, 7, 7, 7, 7]])
    assert allen.vote(neigh).tolist() == [3, 5, 7]
