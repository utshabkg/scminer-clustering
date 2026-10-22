"""Choosing a clustering by subsampling stability, with no labels (U-32).

This follows chooseR (Patterson-Cross et al., BMC Bioinformatics 2021). For each
candidate setting the method is rerun on B random 80 percent subsamples. How
often two cells land in the same cluster gives a co-clustering frequency, and
1 minus that frequency is used as a distance. Scoring the full-data partition in
that distance says how reproducible the setting is. We read the scores two ways. chooseR keeps the finest setting that is still
reproducible, meaning the one with the most clusters whose score clears a
threshold taken from the settings themselves. We also report the plain maximum,
which has no preference for fine partitions. The two bracket how much the choice
depends on that preference.

Two departures from chooseR, both for cost. The bootstrap is a plain percentile
bootstrap over cells rather than BCa, and on datasets above REFERENCE cells the
co-clustering bookkeeping is kept on a fixed random reference subset. Clustering
always runs on real subsamples of the whole dataset.

Nothing here reads the reference labels.
"""
import numpy as np
from sklearn.metrics import silhouette_samples

B_SUBSAMPLES = 25  # chooseR uses 100 and reports that 20 gives the same choice
FRACTION = 0.8
REFERENCE = 3000
BOOTSTRAP = 2000


def reference_cells(n, seed=0):
    """A fixed subset the co-clustering matrix is kept on. All settings share it."""
    if n <= REFERENCE:
        return np.arange(n)
    return np.sort(np.random.default_rng(seed).choice(n, REFERENCE, replace=False))


def subsamples(n, b=B_SUBSAMPLES, fraction=FRACTION, seed=0):
    rng = np.random.default_rng(seed)
    size = int(round(fraction * n))
    return [np.sort(rng.choice(n, size, replace=False)) for _ in range(b)]


class CoClustering:
    """Counts, over the reference cells, how often a pair is seen together and clustered together."""

    def __init__(self, ref):
        self.ref = np.asarray(ref)
        self.pos = {c: i for i, c in enumerate(self.ref)}
        m = len(self.ref)
        self.together = np.zeros((m, m), dtype=np.float32)
        self.seen = np.zeros((m, m), dtype=np.float32)

    def add(self, cells, labels):
        """cells are indices into the full dataset, labels their cluster for this subsample."""
        cells = np.asarray(cells)
        labels = np.asarray(labels)
        keep = np.array([c in self.pos for c in cells])
        if not keep.any():
            return
        rows = np.array([self.pos[c] for c in cells[keep]])
        lab = labels[keep]
        self.seen[np.ix_(rows, rows)] += 1
        same = lab[:, None] == lab[None, :]
        self.together[np.ix_(rows, rows)] += same

    def distance(self):
        with np.errstate(invalid="ignore", divide="ignore"):
            freq = np.where(self.seen > 0, self.together / self.seen, 0.0)
        d = 1.0 - freq
        np.fill_diagonal(d, 0.0)
        return np.clip((d + d.T) / 2, 0, 1).astype(np.float64)


def score_setting(dist, labels_ref, seed=0):
    """chooseR's score: per-cell silhouette in the co-clustering distance, averaged per
    cluster, then the median across clusters. Returns the median and a bootstrap interval."""
    labels_ref = np.asarray(labels_ref)
    uniq = np.unique(labels_ref)
    if len(uniq) < 2:
        return {"stability": -1.0, "ci_low": -1.0, "ci_high": -1.0, "n_clusters": int(len(uniq))}
    sil = silhouette_samples(dist, labels_ref, metric="precomputed")
    per_cluster = np.array([sil[labels_ref == u].mean() for u in uniq])
    rng = np.random.default_rng(seed)
    draws = np.empty(BOOTSTRAP)
    idx_by_cluster = [np.flatnonzero(labels_ref == u) for u in uniq]
    for b in range(BOOTSTRAP):
        means = [sil[rng.choice(ix, len(ix), replace=True)].mean() for ix in idx_by_cluster]
        draws[b] = np.median(means)
    return {"stability": float(np.median(per_cluster)),
            "ci_low": float(np.percentile(draws, 2.5)), "ci_high": float(np.percentile(draws, 97.5)),
            "n_clusters": int(len(uniq))}


def pick(scores):
    """Two label-free readings of the same stability scores.

    finest: chooseR as published. The threshold is the highest lower interval bound
    over the settings, and the choice is the setting with the most clusters that
    clears it. The authors note this leans toward overclustering.
    most_stable: the highest score outright, ties going to fewer clusters. It has no
    preference for fine partitions, so the pair brackets how much the choice depends
    on that preference rather than on the data.
    """
    usable = {k: v for k, v in scores.items() if v["n_clusters"] >= 2}
    if not usable:
        return {}
    threshold = max(v["ci_low"] for v in usable.values())
    passing = {k: v for k, v in usable.items() if v["stability"] >= threshold} or usable
    finest = max(passing, key=lambda k: (passing[k]["n_clusters"], passing[k]["stability"]))
    best = min(usable, key=lambda k: (-usable[k]["stability"], usable[k]["n_clusters"]))
    common = {"threshold": float(threshold), "n_candidates": len(scores), "n_passing": len(passing)}
    return {"stability_finest": {"chosen_setting": finest, **common, **passing[finest]},
            "stability_most": {"chosen_setting": best, **common, **usable[best]}}
