"""What each method gives when run unmodified, with no setting supplied (U-32 companion).

Usage: python scripts/run_rule_default.py <dataset> [...]
This is the weakest label-free rule and the one a user hits first. Scanpy's Leiden
defaults to resolution 1.0. MICA's graph mode defaults to resolution 1.822, which we
score at 1.8, the nearest point on its own grid. MICA's small-data mode defaults to
four clusters: in MICA/mica_mds.py the -nck argument is optional with default [4],
so a user who does not pass it gets four clusters on any dataset, whatever its size.

Writes results/rule_default_<dataset>/.
"""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
from sklearn.cluster import KMeans

from scminer_clust import cluster, diagnostics, metrics
from scminer_clust.paths import PROCESSED

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_rule_stability import N_INIT, dense, mica_distance  # noqa: E402

MICA_DEFAULT_K = 4  # MICA/mica_mds.py, -nck default=[4]
SEEDS = (0, 1, 2)


def main(names):
    rows = []
    for name in names:
        adata = ad.read_h5ad(PROCESSED / f"{name}.h5ad")
        truth = adata.obs["label"].to_numpy()
        D = mica_distance(name, dense(adata))
        emb = diagnostics.classical_mds(D)
        for seed in SEEDS:
            lab = KMeans(n_clusters=MICA_DEFAULT_K, n_init=N_INIT, random_state=seed).fit_predict(emb)
            rec = {"dataset": name, "method": "mica", "k_setting": "default", "seed": seed,
                   "genes": "all", "mode": "small-data", "chosen_setting": MICA_DEFAULT_K,
                   "n_true_classes": int(len(np.unique(truth)))}
            rec.update(metrics.score(truth, lab))
            rows.append(rec)
        best = [r for r in rows if r["dataset"] == name]
        print(f"  {name}: MICA default k={MICA_DEFAULT_K} against {best[0]['n_true_classes']} true types, "
              f"ARI {np.median([r['ari'] for r in best]):.3f}")
    d = metrics.write("rule_default_mica_small_data", rows,
                      {"rule": "the setting MICA uses when none is given", "mica_small_data_k": MICA_DEFAULT_K,
                       "source": "MICA/mica_mds.py -nck default", "seeds": list(SEEDS)})
    print(f"wrote {d / 'metrics.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or ["paper_yan", "paper_goolam", "paper_buettner", "paper_pollen", "paper_chung",
                                   "paper_usoskin", "paper_kolodziejczyk", "paper_klein", "paper_zeisel"]))
