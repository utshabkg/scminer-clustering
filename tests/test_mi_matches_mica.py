"""Our mutual information distance must match the authors' MICA functions.

Run scripts/make_mi_reference.py in the mica-legacy env first, then
`pytest tests/test_mi_matches_mica.py` in the scminer env. This test gates all
later work that relies on our own implementation (U-14, U-17 onward).
"""
import numpy as np
import pytest
from scipy import sparse

from scminer_clust import mi
from scminer_clust.paths import PROCESSED

DATASETS = ["zeisel"]


@pytest.mark.parametrize("name", DATASETS)
@pytest.mark.parametrize("mode", ["mds", "ge"])
def test_matches_mica(name, mode):
    ref_file = PROCESSED / f"mi_reference_{name}.npz"
    if not ref_file.exists():
        pytest.skip(f"run scripts/make_mi_reference.py {name} in mica-legacy first")
    ref = np.load(ref_file)
    X = np.asarray(sparse.load_npz(PROCESSED / f"{name}_exchange" / "X.npz").todense(), dtype=np.float32)
    ours = mi.pair_distances(X, ref["pairs"], mode=mode, bins=int(ref["bins"]))
    theirs = ref[f"d_{mode}"]
    np.testing.assert_allclose(ours, theirs, rtol=0, atol=1e-6)
