"""MICA's distance between two cells must not depend on which other cells are present.

Each cell is binned over its own minimum and maximum, so a pair's distance is a
function of that pair alone. scripts/run_rule_stability.py relies on this: it
computes the full distance matrix once and takes a submatrix per subsample. If
this test ever fails, that shortcut is invalid and the stability runs are wrong.
"""
import anndata as ad
import numpy as np
import pytest
from scipy import sparse

from scminer_clust import mi
from scminer_clust.paths import PROCESSED


@pytest.mark.parametrize("mode", ["mds", "ge"])
def test_submatrix_equals_recomputed(mode):
    a = ad.read_h5ad(PROCESSED / "paper_buettner.h5ad")
    X = np.asarray(a.X.todense() if sparse.issparse(a.X) else a.X, dtype=np.float32)
    full = mi.distance_matrix(X, mode)
    sub = np.sort(np.random.default_rng(0).choice(len(X), len(X) // 2, replace=False))
    recomputed = mi.distance_matrix(X[sub], mode)
    assert np.array_equal(recomputed, full[np.ix_(sub, sub)])
