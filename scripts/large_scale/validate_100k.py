"""Check the Allen pipeline against the 100,000-cell pilot before trusting it on every cell.

Usage: python scripts/large_scale/validate_100k.py   (after run.sh all sample100k)

The pilot loaded the same 100,000 cells whole and used the shared code in src/scminer_clust. Run on
those cells, the streamed pipeline should give back the same genes, values, gene selection,
components and partitions. Writes results/large_scale_sample100k/validation.json.
"""
import json
import sys

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from scminer_clust import cluster
from scminer_clust.paths import PROCESSED, REPO

sys.path.insert(0, str(REPO / "scripts" / "large_scale"))
import allen  # noqa: E402

scope = allen.Scope("sample100k")
man = allen.load_manifest(scope)
genes = pd.read_csv(scope.path("genes.csv.gz"), index_col=0)
hvg = pd.read_csv(scope.path("hvg.csv.gz"), index_col=0)
pilot = ad.read_h5ad(PROCESSED / "allen_100k.h5ad")
out = {}

# 1. the same cells in the same order, and the same genes after the at-least-3-cells filter
out["same_cell_order"] = bool((man.index.astype(str) == pilot.obs_names).all())
kept = genes.index[genes["kept"]]
out["n_genes_kept"] = int(len(kept))
out["same_kept_genes"] = bool(len(kept) == pilot.n_vars and (kept == pilot.var_names).all())

# 2. the same 2,000 genes as cluster.select_hvg(seurat_v3) on the whole pilot matrix
ref_hvg = cluster.select_hvg(pilot, 2000, "seurat_v3")
ours = hvg.index[hvg["highly_variable"]]
out["hvg_overlap"] = int(len(set(ref_hvg) & set(ours)))
out["same_hvg"] = out["hvg_overlap"] == 2000

# 3. the same log2(CPM + 1) values on those genes
x, _ = allen.load_hvg_matrix(scope)
ref_x = pilot[:, list(ours)].X
out["max_abs_diff_log2cpm"] = float(abs(x - ref_x).max())

# 4. the same principal components as cluster.pca with 2,000 seurat_v3 genes, as in the pilot
cluster.HVG_FLAVOR = "seurat_v3"
ref_emb = cluster.pca(pilot, 50, seed=0, n_hvg=2000)
emb = np.load(scope.path("pca.npy"))
corr = [float(np.corrcoef(ref_emb[:, j], emb[:, j])[0, 1]) for j in range(50)]
out["pc_corr_min_first40"] = min(corr[:40])
out["pc_corr"] = corr

# 5. the same k-means and Leiden partitions, scored against the 34 classes
truth = man["class"].to_numpy()
part = scope.runs / "partitions"
km = np.load(part / "kmeans_k34_seed0.npy")
pk = pd.read_csv(allen.PILOT_RUNS / "kmeans_pca_seed0.csv", index_col="cell_label")["label"].reindex(man.index)
out["kmeans_ari_classes"] = float(adjusted_rand_score(truth, km))
out["kmeans_ari_vs_pilot_partition"] = float(adjusted_rand_score(pk, km))
out["pilot_kmeans_ari_classes"] = float(adjusted_rand_score(truth, pk))
le = np.load(part / "leiden_res0.17_seed0.npy")
pl = pd.read_csv(allen.PILOT_RUNS / "leiden_pca_seed0.csv", index_col="cell_label")["label"].reindex(man.index)
out["leiden017_n_clusters"] = int(len(np.unique(le)))
out["leiden017_ari_classes"] = float(adjusted_rand_score(truth, le))
out["leiden017_ari_vs_pilot_partition"] = float(adjusted_rand_score(pl, le))
out["pilot_leiden_ari_classes"] = float(adjusted_rand_score(truth, pl))

out["passed"] = bool(out["same_cell_order"] and out["same_kept_genes"] and out["hvg_overlap"] >= 1998
                     and out["max_abs_diff_log2cpm"] <= 1e-5 and out["pc_corr_min_first40"] >= 0.999
                     and abs(out["kmeans_ari_classes"] - out["pilot_kmeans_ari_classes"]) <= 0.01)
(scope.results / "validation.json").write_text(json.dumps(out, indent=2))
print(json.dumps({k: v for k, v in out.items() if k != "pc_corr"}, indent=2))
