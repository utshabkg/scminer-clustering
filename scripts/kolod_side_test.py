"""U-29: what recovers the published Kolodziejczyk score of 1.0?

Usage: python scripts/kolod_side_test.py

On the authors' released input MICA gets 0.727, and it splits the 2i cells by
plate. The authors' unreleased input had small negative values that shift by
plate, which suggests the plate effect was removed per gene. Here we try that on
the released input and rerun MICA's small-data pipeline (our distance, which
matches theirs to 1e-6, then MDS and k-means at k=3):
  1. the released input as is
  2. each gene centered on its plate mean (location only)
  3. each gene centered and scaled within its plate (location and scale)
Plate is not the label, but it is partly confounded with culture (plate 1 holds
only serum cells, plates 4 and 5 only 2i), so using it can leak label
information. Writes results/kolod_side_test/metrics.csv.
"""
import sys

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_mutual_info_score as ami, adjusted_rand_score as ari

from scminer_clust import diagnostics as dg, metrics, mi
from scminer_clust.paths import PROCESSED


def main():
    a = ad.read_h5ad(PROCESSED / "paper_kolodziejczyk.h5ad")
    X = np.asarray(a.X.todense(), dtype=np.float64)
    y = a.obs["label"].to_numpy()
    # plate is in the cell names, e.g. ola_mES_2i_3_14.counts -> plate 3
    plate = a.obs_names.str.extract(r"_(\d+)_\d+")[0].to_numpy()
    assert not pd.isna(plate).any(), "could not read the plate from every cell name"
    variants = {"released input": X}
    Xc = X.copy()
    for p in np.unique(plate):
        Xc[plate == p] -= Xc[plate == p].mean(axis=0)
    variants["plate-centered genes"] = Xc
    Xs = X.copy()
    for p in np.unique(plate):
        m = plate == p
        sd = Xs[m].std(axis=0); sd[sd == 0] = 1
        Xs[m] = (Xs[m] - Xs[m].mean(axis=0)) / sd
    variants["plate-centered and scaled genes"] = Xs
    rows = []
    for name, V in variants.items():
        D = mi.distance_matrix(V.astype(np.float32), "mds")
        _, labs = dg.mica_small_data(D, 3)
        r = {"variant": name, "ari": float(np.mean([ari(y, l) for l in labs])), "ami": float(np.mean([ami(y, l) for l in labs])),
             "min_value": float(V.min()), "share_negative": float((V < 0).mean()),
             "plate_ari": float(np.mean([ari(plate, l) for l in labs]))}
        rows.append(r)
        print(f"{name:32s} ARI {r['ari']:.3f}  AMI {r['ami']:.3f}  agreement with plate {r['plate_ari']:.3f}  negatives {r['share_negative']:.3f}")
    d = metrics.write("kolod_side_test", rows, {"dataset": "paper_kolodziejczyk", "pipeline": "our MICA distance, MDS, k-means k=3, 50 restarts, seeds 0-2",
                                                "note": "plate correction uses batch labels that partly line up with culture"})
    print(f"wrote {d / 'metrics.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
