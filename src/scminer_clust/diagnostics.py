"""What does MICA's cell-to-cell mutual information measure? (RQ1 diagnostic, U-14)

Two questions.

Co-detection. Single-cell counts are mostly zeros, and MICA bins each cell between
its own minimum and maximum, so the lowest bin usually holds exactly the genes
that were not detected. Let D be 1 where a cell's bin code is above zero. Because
D is a function of the code, I(D_x; D_y) <= I(X; Y), and the ratio says how much
of the similarity comes from which genes are detected in both cells rather than
from how much they are expressed. We also cluster with a detection-only version
of MICA's distance to see how far detection alone gets.

Depth. Binomial thinning keeps each count with probability p, which mimics
sequencing the same cells less deeply. If MICA's clusters change a lot when only
depth changes, the similarity depends on depth.
"""
import numpy as np
from scipy import sparse
from sklearn.cluster import KMeans

from . import mi

MDS_DIMS = 19  # the width of MICA's saved embeddings (mi_reduced.pkl)


def classical_mds(D, dims=MDS_DIMS):
    """MICA's MDS (MICA/lib/transform.py, mds_trans): double-centred squared distances."""
    n = D.shape[0]
    H = np.eye(n) - 1.0 / n
    B = -H @ (D ** 2) @ H / 2
    evals, evecs = np.linalg.eigh(B)
    order = np.argsort(evals)[::-1]
    evals, evecs = evals[order], evecs[:, order]
    pos = evals > 0
    return (evecs[:, pos] * np.sqrt(evals[pos]))[:, :dims]


KMEANS_RESTARTS = 50  # see below


def mica_small_data(D, k, seeds=(0, 1, 2)):
    """MICA's small-data clustering on a given distance matrix: MDS, then k-means.

    MICA keeps the best of 10 k-means starts (its sklearn 0.24 default). Current
    sklearn defaults to one start, and even 10 starts often miss the best fit: on
    Usoskin, 200 single starts give 41 distinct solutions, and fits within 0.2
    percent of the best objective range from ARI 0.79 to 0.87. We use 50 starts so
    the result reliably reaches MICA's own optimum, which makes full and thinned
    runs comparable.
    """
    emb = classical_mds(D)
    return emb, [KMeans(n_clusters=k, n_init=KMEANS_RESTARTS, random_state=s).fit_predict(emb) for s in seeds]


def detection_distance(X, mode="mds", bins=None):
    """MICA's distance computed on detection only: bin codes collapsed to 0 or 1."""
    rule, norm = mi.MODES[mode]
    bins = bins or mi.n_bins(X.shape[1])
    det = (mi.bin_codes(X, bins, rule) > 0).astype(np.int16)
    h = mi.marginal_entropy(det, 2)
    n = det.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        hij = mi.joint_entropy(det[i], det[: i + 1], 2)
        m = h[i] + h[: i + 1] - hij
        with np.errstate(divide="ignore", invalid="ignore"):
            d = 1 - m / np.sqrt(h[i] * h[: i + 1]) if norm == "sqrt" else 1 - m / hij
        D[i, : i + 1] = np.nan_to_num(d, nan=1.0)
        D[: i + 1, i] = D[i, : i + 1]
    np.fill_diagonal(D, 0)
    return D


def codetection_share(X, pairs, bins=None, rule="numpy"):
    """I(D_x; D_y) / I(X; Y) for cell pairs, plus how often bin 0 holds detected genes."""
    X = np.asarray(X)
    bins = bins or mi.n_bins(X.shape[1])
    codes = mi.bin_codes(X, bins, rule)
    det = (codes > 0).astype(np.int16)
    h, hd = mi.marginal_entropy(codes, bins), mi.marginal_entropy(det, 2)
    share = np.empty(len(pairs))
    for k, (i, j) in enumerate(pairs):
        full = h[i] + h[j] - mi.joint_entropy(codes[i], codes[j][None], bins)[0]
        part = hd[i] + hd[j] - mi.joint_entropy(det[i], det[j][None], 2)[0]
        share[k] = part / full if full > 0 else np.nan
    nonzero_in_bin0 = float(((codes == 0) & (X > 0)).sum() / max((X > 0).sum(), 1))
    return share, nonzero_in_bin0


def thin(counts, p, seed=0):
    """Keep each count with probability p (binomial thinning), cell by cell."""
    rng = np.random.default_rng(seed)
    c = sparse.csr_matrix(counts, copy=True).astype(np.int64)
    c.data = rng.binomial(c.data, p)
    c.eliminate_zeros()
    return c.astype(np.float32)


def log_cpm(counts):
    """log2(CPM + 1), the normalization the rest of the project uses."""
    c = sparse.csr_matrix(counts, dtype=np.float64)
    tot = np.asarray(c.sum(axis=1)).ravel()
    c = sparse.diags(1e6 / np.where(tot > 0, tot, 1)) @ c
    c.data = np.log2(c.data + 1)
    return np.asarray(c.todense(), dtype=np.float32)


def knn_overlap(D1, D2, k=20):
    """Mean share of each cell's k nearest neighbors that two distance matrices agree on."""
    n1 = np.argsort(D1, axis=1)[:, 1 : k + 1]
    n2 = np.argsort(D2, axis=1)[:, 1 : k + 1]
    return float(np.mean([len(set(a) & set(b)) / k for a, b in zip(n1, n2)]))
