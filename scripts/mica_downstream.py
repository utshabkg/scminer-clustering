"""Run MICA's own steps after the neighbor graph (node2vec, then Louvain) on any edge list (U-19).

Run in the mica-legacy env:
  python scripts/mica_downstream.py <graph_dir> <out_dir> [workers]
<graph_dir> holds edgelist.txt and cells.txt, as written by export_fast_graph.py.
Uses MICA's library functions and MICA's defaults (node2vec: 20 dimensions, walk
length 60, 110 walks, window 10, p 2.8, q 0.4; Euclidean 20-neighbor graph on the
embedding; Louvain from 0.1 to 4.0 in steps of 0.1, as in our MICA runs). The UMAP
plots are skipped. Writes one clustering_<resolution>.txt per resolution (ID, label).
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from MICA.lib import clustering as cl
from MICA.lib import dimension_reduction as dr
from MICA.lib import neighbor_graph as ng

DIM = 20


def main(graph_dir, out_dir, workers=8):
    graph_dir, out_dir = Path(graph_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cells = np.array(graph_dir.joinpath("cells.txt").read_text().split())
    t0 = time.time()
    emb_file = out_dir / f"knn_node2vec_graph_emb_on_gene_to_{DIM}.txt"
    dr.dim_reduce_node2vec_pecanpy(str(graph_dir / "edgelist.txt"), str(emb_file), dim=DIM, num_jobs=workers,
                                   walk_len=60, n_walks=110, context_size=10, hyper_p=2.8, hyper_q=0.4)
    t_emb = time.time() - t0
    emb = pd.read_csv(emb_file, delimiter=" ", skiprows=1, index_col=0, names=np.arange(1, DIM + 1)).sort_index()
    assert len(emb) == len(cells), f"embedding has {len(emb)} rows for {len(cells)} cells"
    G = ng.build_graph(emb.to_numpy(), dis_metric="euclidean", num_neighbors=20, num_jobs=workers)
    parts = cl.graph_clustering_parallel(G, min_resolution=0.1, max_resolution=4.0, step_size=0.1, num_workers=workers)
    for par, reso in parts:
        labels = [x + 1 for x in par.values()]
        pd.DataFrame({"label": labels}, index=pd.Index(cells[emb.index], name="ID")).to_csv(
            out_dir / f"clustering_{round(reso, 5)}.txt", sep="\t")
    total = time.time() - t0
    (out_dir / "seconds.txt").write_text(f"node2vec {t_emb:.1f}\ntotal {total:.1f}\n")
    print(f"node2vec {t_emb:.0f} s, clustering done in {total:.0f} s, {len(parts)} resolutions, wrote {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 8))
