"""Rerun MICA's final k-means step with fixed seeds. Run inside the mica-legacy env.

Usage: python scripts/mica_seeded_kmeans.py results/reproduction/mica_zeisel 7

MICA's small-data mode saves its embedding (mi_reduced.pkl) and then calls
sklearn KMeans(n_clusters=k) once with no seed (MICA/mica_mds.py), which keeps the
best of 10 starts. That often stops short of the best fit: on Yan, Chung,
Buettner and Usoskin, the published MICA score equals the ARI of the k-means fit
with the lowest objective, while 10 starts can land on fits scoring 0.1 to 0.26
lower. So we run the same call in the same sklearn version with 100 starts and
seeds 0, 1 and 2. That reaches a near-best fit, though not always the very best
(on Pollen the best fit scores 0.980 against the published 0.942). Cell order
comes from MICA's own clustering file.
"""
import pickle
import sys
from pathlib import Path

import pandas as pd
from sklearn.cluster import KMeans


RESTARTS = 100


def main(out_dir, k, seeds=(0, 1, 2)):
    out_dir = Path(out_dir)
    with open(out_dir / "mi_reduced.pkl", "rb") as f:
        emb = pickle.load(f)
    ids = pd.read_csv(next(out_dir.glob("clustering_*.txt")), sep="\t", index_col=0).index
    assert len(ids) == emb.shape[0], "embedding and clustering file disagree on the cell count"
    (out_dir / "seeded").mkdir(exist_ok=True)
    for s in seeds:
        labels = KMeans(n_clusters=k, n_init=RESTARTS, random_state=s).fit(emb).labels_ + 1
        pd.DataFrame({"label": labels}, index=ids).to_csv(out_dir / "seeded" / f"seed_{s}.txt", sep="\t")
    print(f"{out_dir.name}: {len(seeds)} seeded draws, embedding {emb.shape}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], int(sys.argv[2])))
