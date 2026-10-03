"""MICA's cell-to-cell mutual information distance, reimplemented.

MICA treats the genes of two cells as paired samples, puts each cell's values
into equal-width bins between that cell's own minimum and maximum (floor(G^(1/3))
bins for G genes), and estimates mutual information from the joint histogram.

Its two modes differ in two ways, and both are reproduced here:

- Binning. Small-data mode uses numpy's histogram (edges from linspace, the
  maximum in the last bin). Graph mode uses floor(b * (x - min) / (max - min)).
  They agree except for values that land exactly on a bin edge.
- Normalization. Small-data mode uses 1 - I / sqrt(H_x * H_y)
  (MICA/lib/mutual_info.py, mi_norm). Graph mode and the paper's Eq. 3 use
  1 - I / H(x, y) (MICA/lib/distance.py), which is a true metric (Vinh et al.
  2010). The square-root version is not guaranteed to be one.

tests/test_mi_matches_mica.py checks both against the authors' own functions.
"""
import numpy as np

MODES = {"mds": ("numpy", "sqrt"), "ge": ("floor", "joint")}


def n_bins(n_genes):
    return int(n_genes ** (1 / 3.0))


def bin_codes(X, bins, rule="numpy"):
    """Bin each cell (row) over its own range. Returns integer codes in [0, bins)."""
    X = np.asarray(X, dtype=np.float64)
    lo = X.min(axis=1, keepdims=True)
    hi = X.max(axis=1, keepdims=True)
    codes = np.empty(X.shape, dtype=np.int16)
    if rule == "floor":
        width = np.where(hi > lo, hi - lo, 1.0)
        c = np.floor(bins * (X - lo) / width).astype(np.int64)
        codes[:] = np.clip(c, 0, bins - 1)
    else:  # replicate numpy.histogram / histogram2d edge handling exactly
        for i in range(X.shape[0]):
            edges = np.linspace(lo[i, 0], hi[i, 0], bins + 1)
            c = np.searchsorted(edges, X[i], side="right") - 1
            c[X[i] == hi[i, 0]] = bins - 1
            codes[i] = c
    return codes


def _entropy_from_counts(counts, total):
    p = counts / total
    with np.errstate(divide="ignore", invalid="ignore"):
        return -np.where(p > 0, p * np.log(p), 0.0).sum(axis=-1)


def marginal_entropy(codes, bins):
    counts = np.stack([np.bincount(row, minlength=bins) for row in codes])
    return _entropy_from_counts(counts, codes.shape[1])


def joint_entropy(ci, cj_block, bins):
    """Joint entropy of one cell's codes against each row of a block of cells."""
    m, g = cj_block.shape
    comb = ci.astype(np.int64)[None, :] * bins + cj_block.astype(np.int64)
    comb += (np.arange(m, dtype=np.int64) * bins * bins)[:, None]
    counts = np.bincount(comb.ravel(), minlength=m * bins * bins).reshape(m, bins * bins)
    return _entropy_from_counts(counts, g)


def distance_matrix(X, mode="mds", bins=None, block=256):
    """Full cell-by-cell MICA distance. Quadratic, so for small datasets only."""
    rule, norm = MODES[mode]
    bins = bins or n_bins(np.asarray(X).shape[1])
    codes = bin_codes(X, bins, rule)
    h = marginal_entropy(codes, bins)
    n = codes.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for s in range(0, i + 1, block):
            e = min(i + 1, s + block)
            hij = joint_entropy(codes[i], codes[s:e], bins)
            mi = h[i] + h[s:e] - hij
            if norm == "sqrt":
                d = 1 - mi / np.sqrt(h[i] * h[s:e])
            else:
                d = 1 - mi / hij
            D[i, s:e] = d
            D[s:e, i] = d
    return D


def pair_distances(X, pairs, mode="mds", bins=None):
    """Distances for a list of (i, j) cell pairs, without the full matrix."""
    rule, norm = MODES[mode]
    bins = bins or n_bins(np.asarray(X).shape[1])
    codes = bin_codes(X, bins, rule)
    h = marginal_entropy(codes, bins)
    out = np.empty(len(pairs))
    for k, (i, j) in enumerate(pairs):
        hij = joint_entropy(codes[i], codes[j][None, :], bins)[0]
        mi = h[i] + h[j] - hij
        out[k] = 1 - mi / np.sqrt(h[i] * h[j]) if norm == "sqrt" else 1 - mi / hij
    return out
