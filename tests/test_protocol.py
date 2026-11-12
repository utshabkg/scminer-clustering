"""Protocol v1 rests on a few exact claims (docs/protocol_matched_v1.md). These tests check them.

1. The first 25 of 100 subsamples are the 25 the original protocol drew.
2. The matrix-product co-clustering counts equal stability.CoClustering.add, so the distance
   and the score are identical, with and without a reference subset.
3. Subsample labels can be stored as codes: only equality matters for co-clustering.
4. The candidate variants: mid is the original one_per_count, low and high bracket it,
   dense holds every setting, and a method reaching every count leaves the matched counts alone.
5. chooseR's frequency is together / seen where the pair was ever together, else 0.
6. pick_chooser follows the example script: highest lower bound as threshold, then the most
   clusters among settings whose median clears it, ties to the highest setting value, and no
   choice when any interval is missing.
"""
import numpy as np
import pytest

from scminer_clust import grids, protocol, stability


@pytest.mark.parametrize("n", [124, 2903, 13605])
@pytest.mark.parametrize("seed", [0, 3, 9])
def test_subsamples_nest(n, seed):
    a = stability.subsamples(n, b=100, seed=seed)[:25]
    b = stability.subsamples(n, b=25, seed=seed)
    assert all(np.array_equal(x, y) for x, y in zip(a, b))


def _random_runs(n, b, seed, k=4):
    rng = np.random.default_rng(seed)
    cells = stability.subsamples(n, b=b, seed=seed)
    labels = [rng.integers(0, k, len(c)).astype(str) for c in cells]  # strings, as Leiden returns
    return cells, labels


@pytest.mark.parametrize("n,ref_n", [(150, None), (3400, 3000)])
def test_counter_equals_coclustering(n, ref_n):
    cells, labels = _random_runs(n, 12, seed=1)
    ref = stability.reference_cells(n, seed=1)
    assert ref_n is None or len(ref) == ref_n
    co = stability.CoClustering(ref)
    for c, lab in zip(cells, labels):
        co.add(c, lab)
    ctr = protocol.Counter(n, ref, cells, 12)
    tog = ctr.together(np.stack([protocol.encode(lab) for lab in labels]))
    assert np.array_equal(tog, co.together) and np.array_equal(ctr.seen, co.seen)
    assert np.array_equal(protocol.distance_ours(tog, ctr.seen), co.distance())
    full = np.random.default_rng(2).integers(0, 3, n).astype(str)[ref]
    assert stability.score_setting(protocol.distance_ours(tog, ctr.seen), full, seed=1) == \
        stability.score_setting(co.distance(), full, seed=1)


def test_encode_keeps_equality():
    lab = np.array(["10", "2", "2", "0", "10"])
    code = protocol.encode(lab)
    assert np.array_equal(lab[:, None] == lab[None, :], code[:, None] == code[None, :])


def test_candidate_variants():
    rng = np.random.default_rng(0)
    full = {r: rng.integers(0, 2 + int(r * 10) % 6, 50) for r in np.round(np.arange(0.1, 3.0, 0.1), 1)}
    counts = sorted({len(np.unique(v)) for v in full.values()})
    assert grids.per_count(full, counts, "mid") == grids.one_per_count(full, counts)
    low, high, dense = grids.per_count(full, counts, "low"), grids.per_count(full, counts, "high"), \
        grids.dense_per_count(full, counts)
    for k in counts:
        assert low[k] <= grids.one_per_count(full, counts)[k] <= high[k]
        assert low[k] == dense[k][0] and high[k] == dense[k][-1]
    kmeans_like = {k: np.arange(50) % k for k in range(2, 21)}
    assert grids.common_counts(full, kmeans_like) == grids.common_counts(full) == [c for c in counts if 2 <= c <= 20]


def test_freq_chooser():
    together = np.array([[3, 0, 2], [0, 4, 0], [2, 0, 5]], dtype=np.float32)
    seen = np.array([[3, 2, 3], [2, 4, 0], [3, 0, 5]], dtype=np.float32)
    f = protocol.freq_chooser(together, seen)
    assert f.dtype == np.float64
    assert np.array_equal(f, np.array([[1, 0, 2 / 3], [0, 1, 0], [2 / 3, 0, 1]]))


def test_pick_chooser():
    rows = [{"setting": 0.5, "n_clusters": 3, "med": 0.90, "low_med": 0.80},
            {"setting": 0.7, "n_clusters": 5, "med": 0.86, "low_med": 0.70},
            {"setting": 0.8, "n_clusters": 5, "med": 0.85, "low_med": 0.60},
            {"setting": 0.9, "n_clusters": 7, "med": 0.70, "low_med": 0.65}]
    assert protocol.pick_chooser(rows) == (0.8, 0.80)  # 5 clusters, the higher setting of the tie
    rows[1]["low_med"] = None
    assert protocol.pick_chooser(rows)[0] is None
