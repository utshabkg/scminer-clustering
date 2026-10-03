"""How sure are we about each Leiden minus MICA gap? (U-28)

Usage: python scripts/bootstrap_gaps.py [n_boot]

For each of the ten benchmark datasets on the authors' own input and each gene
setting, this resamples cells with replacement and recomputes ARI and AMI for
both methods on the same resampled cells (a paired bootstrap), averaging over the
saved seeds. It writes the observed gap, a 95 percent interval and the share of
resamples where Leiden is ahead to results/uncertainty/bootstrap_gaps.csv, and a
paired Wilcoxon test across datasets for each gene setting.

The bootstrap measures how much the score depends on which cells were sampled.
It does not measure run-to-run randomness of the methods; repeated runs cover
that. Partitions come from the paper's selection rule (true number of clusters).
"""
import json
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.metrics import adjusted_mutual_info_score as ami, adjusted_rand_score as ari

from scminer_clust.paths import DATA, PROCESSED, REPO, RESULTS

DATASETS = ["yan", "goolam", "buettner", "pollen", "chung", "usoskin", "kolodziejczyk", "klein", "zeisel", "pbmc14k"]
SETTINGS = {"all genes": "", "2000 genes (seurat)": "_hvg2000", "2000 genes (cell_ranger)": "_hvg2000_cell_ranger",
            "2000 genes (seurat_v3)": "_hvg2000_seurat_v3"}


def mica_partitions(ds, suffix):
    """MICA's partitions picked by the paper's rule: seeded draws, or the chosen graph-mode setting."""
    tag = suffix.replace("_hvg2000", "hvg2000_seurat") if suffix == "_hvg2000" else suffix.lstrip("_")
    run = REPO / "results" / "reproduction" / (f"mica_paper_{ds}" + (f"_{tag}" if tag else ""))
    seeded = sorted((run / "seeded").glob("seed_*.txt"))
    if seeded:
        return [pd.read_csv(f, sep="\t", index_col=0)["label"] for f in seeded]
    exp = f"reproduction_paper_{ds}" if not tag else f"mica_paper_{ds}_{tag}"
    res = json.loads((RESULTS / exp / "params.json").read_text())["resolution"]
    f = next(run.glob(f"clustering_*_{res:g}.txt"), None) or next(run.glob(f"clustering_*_{res}.txt"))
    return [pd.read_csv(f, sep="\t", index_col=0)["label"]]


def leiden_partitions(ds, suffix):
    runs = DATA / "runs" / f"baselines_paper_{ds}{suffix}"
    return [pd.read_csv(f, index_col=0)["label"] for f in sorted(runs.glob("leiden_pca_seed*.csv"))]


def main(n_boot=1000):
    rng = np.random.default_rng(0)
    rows = []
    for ds in DATASETS:
        truth = ad.read_h5ad(PROCESSED / f"paper_{ds}.h5ad", backed="r").obs["label"]
        cells = truth.index
        y = truth.to_numpy()
        for setting, suffix in SETTINGS.items():
            try:
                M = [p.astype(str).reindex(cells).to_numpy() for p in mica_partitions(ds, suffix)]
                L = [p.astype(str).reindex(cells).to_numpy() for p in leiden_partitions(ds, suffix)]
            except (StopIteration, FileNotFoundError, KeyError):
                print(f"skip {ds} {setting}: missing partitions"); continue
            if not M or not L or any(pd.isna(p).any() for p in M + L):
                print(f"skip {ds} {setting}: missing or misaligned partitions"); continue
            obs = {f: np.mean([s(y, p) for p in L]) - np.mean([s(y, p) for p in M]) for f, s in (("ari", ari), ("ami", ami))}
            boot = {"ari": [], "ami": []}
            for _ in range(n_boot):
                idx = rng.integers(0, len(y), len(y))
                for f, s in (("ari", ari), ("ami", ami)):
                    boot[f].append(np.mean([s(y[idx], p[idx]) for p in L]) - np.mean([s(y[idx], p[idx]) for p in M]))
            r = {"dataset": ds, "gene_setting": setting, "n_cells": len(y), "mica_draws": len(M), "leiden_seeds": len(L)}
            for f in ("ari", "ami"):
                b = np.asarray(boot[f])
                r.update({f"gap_{f}": round(obs[f], 4), f"ci_low_{f}": round(np.percentile(b, 2.5), 4),
                          f"ci_high_{f}": round(np.percentile(b, 97.5), 4), f"share_leiden_ahead_{f}": round(float((b > 0).mean()), 3)})
            rows.append(r)
            print(f"{ds:14s} {setting:26s} Leiden-MICA ARI {obs['ari']:+.3f} [{r['ci_low_ari']:+.3f}, {r['ci_high_ari']:+.3f}]")
    out = RESULTS / "uncertainty"; out.mkdir(exist_ok=True)
    df = pd.DataFrame(rows); df.to_csv(out / "bootstrap_gaps.csv", index=False)
    tests = []
    for setting in SETTINGS:
        d = df[df.gene_setting == setting]
        for f in ("ari", "ami"):
            gaps = d[f"gap_{f}"].to_numpy()
            p = wilcoxon(gaps).pvalue if len(gaps) >= 5 and np.any(gaps != 0) else np.nan
            tests.append({"gene_setting": setting, "metric": f, "n_datasets": len(gaps), "mean_gap": round(gaps.mean(), 4),
                          "median_gap": round(np.median(gaps), 4), "leiden_ahead": int((gaps > 0).sum()),
                          "ci_excludes_zero": int(((d[f"ci_low_{f}"] > 0) | (d[f"ci_high_{f}"] < 0)).sum()), "wilcoxon_p": p})
    pd.DataFrame(tests).to_csv(out / "dataset_level_tests.csv", index=False)
    print(pd.DataFrame(tests).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else 1000))
