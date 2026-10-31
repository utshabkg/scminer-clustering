"""The silhouette rule, scored for both methods in one shared space (U-35).

Usage: python scripts/run_rule_silhouette.py <dataset> [n_hvg] [hvg_flavor] [--matched]
In Assignment 2 each method's silhouette was computed in its own embedding, so the
score mixed the quality of the partition with the geometry of the space. Here every
partition from both methods is scored in the same space, the 50 principal
components, and the rule keeps the setting with the highest mean silhouette.

To show how much the space matters, each partition is also scored in MICA's MDS
space. The rule's choice under each space is reported, so the two can be compared.
No labels are used to choose. Writes results/rule_silhouette_<exp>/.

With --matched both methods are offered the same cluster counts with one setting per count,
as in the matched stability runs (src/scminer_clust/grids.py), and the run writes
results/rule_matched_silhouette_<exp>/. The Assignment 2 plan promised this rule over the
same range of cluster counts for both methods, which a shared cap does not give.
"""
import sys
from pathlib import Path

import anndata as ad
import numpy as np

from scminer_clust import cluster, diagnostics, grids, metrics, selection
from scminer_clust.paths import PROCESSED

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_rule_stability import dense, leiden_partitions, mica_distance, mica_partitions  # noqa: E402

SEED = 0


def main(name, n_hvg=None, flavor="seurat"):
    cluster.N_HVG = int(n_hvg) if n_hvg else None
    cluster.HVG_FLAVOR = flavor
    tag = f"{name}" + (f"_hvg{n_hvg}_{flavor}" if n_hvg else "")
    adata = ad.read_h5ad(PROCESSED / f"{name}.h5ad")
    truth = adata.obs["label"].to_numpy()
    n = adata.n_obs
    genes = cluster.select_hvg(adata, int(n_hvg), flavor) if n_hvg else None
    X = dense(adata[:, genes] if genes is not None else adata)
    D = mica_distance(tag, X)
    allc = np.arange(n)

    spaces = {"pca": cluster.pca(adata, seed=SEED), "mds": diagnostics.classical_mds(D)}
    _, mica_parts = mica_partitions(D, allc)
    _, leiden_parts = leiden_partitions(adata, allc)
    partitions = {"mica": mica_parts, "leiden_pca": leiden_parts}
    chosen = {}
    if grids.matched():
        common = grids.common_counts(mica_parts, leiden_parts)
        chosen = {m: set(grids.one_per_count(parts, common).values()) for m, parts in partitions.items()}
        print(f"  matched counts: {common} ({len(common)} of 19)", flush=True)

    rows = []
    for method, parts in partitions.items():
        keep = ({s: parts[s] for s in chosen[method]} if chosen else
                {s: lab for s, lab in parts.items() if 2 <= len(np.unique(lab)) <= 20})
        for space, emb in spaces.items():
            key, info = selection.pick_by_silhouette(emb, keep, seed=SEED)
            rec = {"dataset": name, "method": method, "k_setting": f"silhouette_{space}", "space": space,
                   "seed": SEED, "genes": int(n_hvg) if n_hvg else "all", "n_candidates": len(keep), **info}
            rec.update(metrics.score(truth, keep[key]))
            rows.append(rec)
            print(f"  {method} scored in {space}: setting {key}, {rec['n_clusters']} clusters, "
                  f"ARI {rec['ari']:.3f}", flush=True)
    # The shared space is the principal components, so that scoring is the rule itself.
    # Iterate over a copy: appending to the list being iterated would never terminate.
    rows += [{**r, "k_setting": "silhouette_shared"} for r in list(rows) if r["space"] == "pca"]
    d = metrics.write(grids.tag(f"rule_silhouette_{tag}"), rows,
                      {"dataset": name, "rule": "highest mean silhouette, scored in a shared space",
                       "grid": "matched" if grids.matched() else "original",
                       "shared_space": "50 principal components", "also_scored_in": "MICA MDS", "seed": SEED})
    print(f"wrote {d / 'metrics.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None,
                  sys.argv[3] if len(sys.argv) > 3 else "seurat"))
