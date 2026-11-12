"""Protocol v1: run the stability rule once at B=100 and score the saved partitions many ways.

Written for Jian Liu's Assignment 3 feedback (Nov 2026). The design is fixed in
docs/protocol_matched_v1.md and tagged matched-protocol-v1 before the confirmatory runs.

Why one run is enough:
- stability.subsamples draws its subsamples in sequence from one generator, and every
  subsample is clustered with the same fixed seed. The first 25 subsamples of a B=100 run
  are therefore exactly the B=25 run of the original protocol (tests/test_protocol.py).
- A setting's score depends only on that setting's partitions, so any candidate set
  (original, matched middle, low, high, dense) can be scored after the fact from one
  curve over the union of settings. Only stability.pick looks at the whole set.

The saved store keeps, per method, every subsample's labels as int16 codes. Codes only need
to preserve which cells share a cluster. The full-data labels keep their original values,
because their sorted order sets the order of the per-cluster bootstrap in
stability.score_setting.

Nothing here reads the reference labels.
"""
import hashlib
import json
import os
import platform
import subprocess
import sys

import numpy as np

from . import stability
from .paths import DATA, REPO

STORE = DATA / "runs" / "protocol_v1"
TAG = "matched-protocol-v1"
B_FULL = 100
B_OLD = stability.B_SUBSAMPLES  # 25, the original protocol
METHODS = ("mica", "leiden_pca", "kmeans_pca")
DESIGNS = ("original", "mid", "low", "high", "dense")


def encode(labels):
    """Subsample labels as int16 codes. Only equality matters for co-clustering."""
    _, inv = np.unique(np.asarray(labels), return_inverse=True)
    assert inv.max() < np.iinfo(np.int16).max
    return inv.astype(np.int16)


def save_store(folder, cells, methods):
    """cells: B x size int array. methods: {name: {"settings": [...], "full": {s: labels},
    "codes": B x S x size int16}}. Uncompressed npy so the scorers can memory-map it."""
    folder.mkdir(parents=True, exist_ok=True)
    np.save(folder / "cells.npy", np.asarray(cells, dtype=np.int32))
    for m, d in methods.items():
        np.save(folder / f"{m}_codes.npy", d["codes"])
        np.save(folder / f"{m}_settings.npy", np.asarray(d["settings"], dtype=np.float64))
        full = np.stack([np.asarray(d["full"][s]) for s in d["settings"]])
        np.save(folder / f"{m}_full.npy", full)  # original dtype: ints or fixed-width strings


def load_store(folder, method):
    cells = np.load(folder / "cells.npy", mmap_mode="r")
    codes = np.load(folder / f"{method}_codes.npy", mmap_mode="r")
    settings = np.load(folder / f"{method}_settings.npy").tolist()
    full = np.load(folder / f"{method}_full.npy")
    return cells, codes, settings, full


class Counter:
    """together and seen counts over the reference cells, from the first b subsamples.

    One float32 matrix product replaces the per-subsample updates of
    stability.CoClustering.add. Counts are whole numbers far below 2**24, so float32 holds
    them exactly and the result equals what add() accumulates (tests/test_protocol.py).
    """

    def __init__(self, n, ref, cells, b):
        self.ref = np.asarray(ref)
        self.pos = np.full(n, -1, dtype=np.int64)
        self.pos[self.ref] = np.arange(len(self.ref))
        self.b = b
        m = len(self.ref)
        self.rows, self.keep = [], []
        P = np.zeros((b, m), dtype=np.float32)
        for t in range(b):
            idx = self.pos[np.asarray(cells[t])]
            keep = idx >= 0
            self.rows.append(idx[keep])
            self.keep.append(keep)
            P[t, idx[keep]] = 1.0
        self.seen = P.T @ P

    def together(self, codes):
        """codes: b x size array of one setting's subsample codes, aligned with cells."""
        blocks, offset = [], 0
        m = len(self.ref)
        sizes = [int(np.asarray(codes[t])[self.keep[t]].max()) + 1 if self.keep[t].any() else 0 for t in range(self.b)]
        M = np.zeros((m, sum(sizes)), dtype=np.float32)
        for t in range(self.b):
            if not sizes[t]:
                continue
            lab = np.asarray(codes[t])[self.keep[t]].astype(np.int64)
            M[self.rows[t], offset + lab] = 1.0
            offset += sizes[t]
        return M @ M.T


def distance_ours(together, seen):
    """stability.CoClustering.distance on given counts, byte for byte."""
    co = stability.CoClustering.__new__(stability.CoClustering)
    co.together, co.seen = together, seen
    return co.distance()


def freq_chooser(together, seen):
    """chooseR's co-clustering frequency in float64: Re(x) / (n - Im(x)) where the pair was
    ever together, else 0 (examples/1_seurat_pipeline.R, the "Scoring" step)."""
    t = together.astype(np.float64)
    s = seen.astype(np.float64)
    out = np.zeros_like(t)
    np.divide(t, s, out=out, where=t > 0)
    return out


def pick_chooser(rows):
    """chooseR's choice, ported from examples/1_seurat_pipeline.R:
        threshold <- max(meds$low_med)
        meds %>% filter(med >= threshold) %>% arrange(n_clusters) %>% tail(n = 1) %>% pull(res)
    rows: list of dicts with setting, n_clusters, med, low_med. group_by(res) sorts by res and
    arrange() is stable, so ties in n_clusters go to the highest setting value. A missing
    low_med (a failed bootstrap) makes max() NA in R, so no choice is made, returned as None.
    """
    if not rows or any(r["low_med"] is None or np.isnan(r["low_med"]) for r in rows):
        return None, float("nan")
    threshold = max(r["low_med"] for r in rows)
    passing = sorted((r for r in rows if r["med"] >= threshold), key=lambda r: (r["n_clusters"], r["setting"]))
    return (passing[-1]["setting"] if passing else None), threshold


def sha256(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


# The files a protocol v1 job executes. require_tag compares exactly these with the tag, so notes,
# the manuscript and the analysis script can change while the runs go on.
CODE = ["src/scminer_clust", "scripts/run_protocol_small.py", "scripts/run_protocol_pbmc.py",
        "scripts/run_rule_stability.py", "scripts/run_rule_stability_graph.py", "scripts/mica_downstream.py",
        "scripts/run_chooser_r.py", "scripts/score_chooser_r.R", "scripts/run_protocol_v1.sh",
        "scripts/fetch_chooser.sh", "envs/chooser-r.yml", "envs/chooser-r.lock.yml", "docs/protocol_matched_v1.md"]


def git_state():
    def run(*a):
        r = subprocess.run(["git", *a], cwd=REPO, capture_output=True, text=True)
        return r.stdout.strip()
    tag_commit = run("rev-list", "-n", "1", TAG) if run("tag", "-l", TAG) else ""
    changed = run("diff", "--name-only", TAG, "--", *CODE).split() if tag_commit else ["(no tag)"]
    untracked = run("ls-files", "--others", "--exclude-standard", "--", *CODE).split()
    return {"head": run("rev-parse", "HEAD"), "describe": run("describe", "--tags", "--always", "--dirty"),
            "dirty": bool(run("status", "--porcelain", "--untracked-files=no")), "tag_commit": tag_commit,
            "code_changed_since_tag": changed, "code_untracked": untracked,
            "code_matches_tag": bool(tag_commit) and not changed and not untracked}


def require_tag():
    """Refuse to run a confirmatory job unless the code it executes is exactly the tagged code."""
    g = git_state()
    if not g["code_matches_tag"]:
        sys.exit(f"refusing to run: the protocol code differs from {TAG}: changed {g['code_changed_since_tag']}, "
                 f"untracked {g['code_untracked']}")
    return g


def provenance(extra=None):
    import anndata
    import scanpy
    import sklearn
    return {"git": git_state(), "module": __file__, "repo": str(REPO),
            "python": sys.version.split()[0], "platform": platform.platform(),
            "numpy": np.__version__, "scanpy": scanpy.__version__, "sklearn": sklearn.__version__,
            "anndata": anndata.__version__, "pythonpath": os.environ.get("PYTHONPATH", ""),
            "threads": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                                                       "MKL_NUM_THREADS", "NUMBA_NUM_THREADS")},
            **(extra or {})}


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, default=str))
