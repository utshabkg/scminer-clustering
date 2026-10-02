"""Scoring and run bookkeeping. Every experiment writes the same two files."""
import json
import resource
import time
from contextlib import contextmanager

import numpy as np
import pandas as pd
from sklearn import metrics as skm

from .paths import results_dir


def score(labels_true, labels_pred, embedding=None):
    """Agreement with reference labels, plus silhouette if an embedding is given."""
    out = {
        "n_clusters": int(len(np.unique(labels_pred))),
        "n_true_classes": int(len(np.unique(labels_true))),
        "ari": float(skm.adjusted_rand_score(labels_true, labels_pred)),
        "ami": float(skm.adjusted_mutual_info_score(labels_true, labels_pred)),
        "nmi": float(skm.normalized_mutual_info_score(labels_true, labels_pred)),
    }
    if embedding is not None and out["n_clusters"] > 1:
        out["silhouette"] = float(skm.silhouette_score(embedding, labels_pred))
    return out


def per_type_recall(labels_true, labels_pred):
    """For each reference type, the share of its cells in the cluster that holds most of them.

    This is the rare-population view: a type that gets scattered scores low even
    when the overall agreement looks fine.
    """
    rows = []
    labels_true, labels_pred = np.asarray(labels_true), np.asarray(labels_pred)
    for t in np.unique(labels_true):
        mask = labels_true == t
        clusters, counts = np.unique(labels_pred[mask], return_counts=True)
        rows.append({
            "type": t,
            "n_cells": int(mask.sum()),
            "share_of_data": float(mask.mean()),
            "best_cluster": clusters[counts.argmax()],
            "recall": float(counts.max() / mask.sum()),
            "n_clusters_spanned": int(len(clusters)),
        })
    return pd.DataFrame(rows).sort_values("n_cells")


@contextmanager
def timed(record: dict, key: str = "seconds"):
    start = time.perf_counter()
    yield record
    record[key] = round(time.perf_counter() - start, 2)
    record["peak_rss_gb"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2, 2)


def write(exp_id: str, rows, params: dict):
    """Write metrics.csv and params.json for one experiment."""
    d = results_dir(exp_id)
    df = pd.DataFrame(rows if isinstance(rows, list) else [rows])
    df.to_csv(d / "metrics.csv", index=False)
    (d / "params.json").write_text(json.dumps(params, indent=2, default=str))
    return d
