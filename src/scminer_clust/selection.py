"""Choosing a partition without labels (RQ1, U-12).

The paper picks the number of clusters with the true labels. Here every method
picks with the same label-free rule instead: among its candidate partitions,
keep the one with the highest silhouette in that method's own embedding (MICA's
MDS or node2vec embedding, principal components for the baselines). Candidates
are k = 2..K_MAX for k-based methods and the resolution sweep for graph methods.
The range is fixed for every dataset, so it does not leak the true number.
"""
import numpy as np
from sklearn.metrics import silhouette_score

K_MAX = 20
SAMPLE = 3000  # silhouette on a fixed subsample, for speed on the larger datasets


def silhouette(emb, labels, seed=0):
    labels = np.asarray(labels)
    if len(np.unique(labels)) < 2:
        return -1.0
    n = min(SAMPLE, len(labels))
    return float(silhouette_score(emb, labels, sample_size=n, random_state=seed))


def pick_by_silhouette(emb, partitions, seed=0):
    """partitions maps a setting (k or resolution) to labels. Ties go to fewer clusters."""
    scored = []
    for key, labels in partitions.items():
        scored.append((-silhouette(emb, labels, seed), len(np.unique(labels)), key))
    s, n_clusters, key = min(scored)
    return key, {"silhouette": -s, "n_candidates": len(partitions), "chosen_setting": key}
