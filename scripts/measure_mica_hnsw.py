"""Peak memory and time of MICA's own HNSW neighbour search, separated from the rest of MICA.

Usage (mica-legacy environment): python scripts/measure_mica_hnsw.py <dataset> [workers]

The other half of the graph cost comparison, beside scripts/measure_graph_cost.py. It calls
hnswlib_ann from MICA's own mica_ge.py, with the settings `mica ge` uses as we ran it
(80 neighbours, ef 400, M 8, bin power and bin size 0, 8 workers as in run_mica_job.sh), on
the matrix `mica ge` reads. That matrix is rebuilt in memory from the exchange folder exactly
as scripts/export_for_mica.py builds it, so no temporary h5ad is written.

Only the neighbour search is measured. node2vec and Louvain run on both graphs and are timed
separately in results/rq2_downstream_*.

Peak resident memory is read from the kernel and is a high-water mark. hnswlib runs its
threads inside this process, so they are counted.

This script runs under Python 3.9 in the authors' environment, which cannot import
scminer_clust, so it writes its own results folder.
Writes results/graph_cost_<dataset>_mica_hnsw/.
"""
import json
import resource
import sys
import time
from pathlib import Path

import anndata
import hnswlib
import numpy as np
import pandas as pd
from scipy import sparse
from MICA.mica_ge import hnswlib_ann

ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
SETTINGS = {"num_neighbors": 80, "ef": 400, "M": 8, "bin_power": 0, "bin_size": 0}  # mica ge defaults


def peak_gb():
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2, 2)


def main(name, workers=8):
    workers = int(workers)
    t0 = time.time()
    exch = PROCESSED / f"{name}_exchange"
    X = sparse.load_npz(exch / "X.npz")
    cells = (exch / "cells.txt").read_text().splitlines()
    genes = (exch / "genes.txt").read_text().splitlines()
    adata = anndata.AnnData(X=np.asarray(X.todense(), dtype=np.float32),
                            obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))
    del X
    data = adata.to_df().to_numpy()   # what mica_ge passes to hnswlib_ann
    del adata
    n, g = data.shape
    t_load, m_load = time.time() - t0, peak_gb()
    print(f"{name}: {n} cells, {g} genes, loaded in {t_load:.0f} s, peak {m_load} GB", flush=True)

    hnswlib.set_bin_power(SETTINGS["bin_power"])
    hnswlib.set_bin_size(SETTINGS["bin_size"])
    t0 = time.time()
    idx, dst = hnswlib_ann(data, num_neighbors=SETTINGS["num_neighbors"], ef=SETTINGS["ef"],
                           M=SETTINGS["M"], num_jobs=workers)
    t_knn, m_all = time.time() - t0, peak_gb()
    assert idx.shape == (n, SETTINGS["num_neighbors"]), idx.shape
    print(f"  MICA HNSW: {t_knn:.0f} s, peak {m_all} GB", flush=True)

    out = ROOT / "results" / f"graph_cost_{name}_mica_hnsw"
    out.mkdir(parents=True, exist_ok=True)
    row = {"dataset": name, "graph": "mica_hnsw", "k": SETTINGS["num_neighbors"] - 1, "workers": workers,
           "n_cells": n, "n_genes": g, "seconds_load_and_bin": round(t_load, 1), "seconds_candidates": "",
           "seconds_rescoring": "", "seconds_graph": round(t_knn, 1), "peak_rss_gb_after_load": m_load,
           "peak_rss_gb_after_candidates": "", "peak_rss_gb": m_all}
    pd.DataFrame([row]).to_csv(out / "metrics.csv", index=False)
    (out / "params.json").write_text(json.dumps(
        {"dataset": name, "workers": workers, **SETTINGS,
         "call": "MICA.mica_ge.hnswlib_ann on adata.to_df().to_numpy(), as mica_ge does",
         "memory": "resource.getrusage(RUSAGE_SELF).ru_maxrss, a high-water mark"}, indent=2))
    print(f"wrote {out / 'metrics.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:3]))
