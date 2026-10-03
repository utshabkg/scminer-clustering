"""Score the authors' MICA output with the same rules as our baselines.

Usage: python scripts/score_mica.py zeisel results/reproduction/mica_zeisel [tag]

Without a tag this is a reproduction of a published number, written to
results/reproduction_<dataset>/. With a tag such as hvg2000_seurat it is MICA on
a gene subset, written to results/mica_<dataset>_<tag>/.

MICA writes one clustering_*.txt per setting (ID, X, Y, label). With a
resolution sweep there is one file per resolution, and we pick among them with
the paper's protocol, which is what the published numbers used.
"""
import re
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

from scminer_clust import cluster, metrics
from scminer_clust.paths import PROCESSED

# The paper's own per-dataset scores, from the authors' Benchmark_ARI_metrics files
PUBLISHED_ARI = pd.read_csv(Path(__file__).resolve().parents[1] / "docs" / "published_scores.csv",
                            index_col="dataset")["MICA_ARI"].to_dict()


def _fed_genes(name, tag):
    f = PROCESSED / f"{name}{'_' + tag if tag else ''}_exchange" / "genes.txt"
    return f.read_text().splitlines() if f.exists() else None


def main(name, out_dir, tag=None, exp_id=None):
    out_dir = Path(out_dir)
    truth = ad.read_h5ad(PROCESSED / f"{name}.h5ad", backed="r").obs["label"]
    files = sorted(out_dir.glob("clustering_*.txt"))  # seeded draws live in out_dir/seeded/
    if not files:
        print(f"no clustering_*.txt in {out_dir}")
        return 1

    exch = PROCESSED / f"{name}{'_' + tag if tag else ''}_exchange" / "cells.txt"
    fed = exch.read_text().splitlines() if exch.exists() else None

    def aligned(df, where):
        """MICA's labels in the order of our cells, or a clear error if the cells differ."""
        ids = df.index.astype(str)
        if ids.duplicated().any():
            raise ValueError(f"{where}: duplicate cell IDs in MICA output")
        if set(ids) != set(truth.index):
            raise ValueError(f"{where}: MICA clustered {len(ids)} cells, the processed file has {len(truth)}; "
                             f"{len(set(ids) - set(truth.index))} unknown, {len(set(truth.index) - set(ids))} missing")
        if fed is not None and set(ids) != set(fed):
            raise ValueError(f"{where}: MICA output does not match the cells it was given ({exch})")
        return df.set_axis(ids)["label"].reindex(truth.index).astype(str).to_numpy()

    partitions = {}
    for f in files:
        setting = f.stem.split("_")[-1]
        partitions[setting] = aligned(pd.read_csv(f, sep="\t", index_col=0), f.name)

    setting, info = cluster.pick_paper_protocol(partitions, truth.to_numpy())
    rec = {"dataset": name, "method": "mica", "seed": 0, "k_setting": "paper_protocol",
           "mica_setting": setting, "n_settings_tried": len(partitions)}
    m = re.match(r"hvg(\d+)_(.+)", tag or "")
    rec["genes"], rec["hvg_flavor"] = (int(m.group(1)), m.group(2)) if m else ("all", "")
    rec.update(metrics.score(truth.to_numpy(), partitions[setting]))
    log = Path("logs") / f"{out_dir.name}.log"
    if log.exists():
        wall = re.findall(r"WALL_SECONDS (\d+)", log.read_text())
        rec["seconds"] = int(wall[-1]) if wall else None
    if tag is None:
        rec["published_ari"] = PUBLISHED_ARI.get(name.removeprefix("paper_"))
        rec["gap_vs_published"] = round(rec["ari"] - rec["published_ari"], 3) if rec["published_ari"] else None

    rec["n_cells"], rec["n_genes"] = len(truth), (len(fed_genes) if (fed_genes := _fed_genes(name, tag)) else None)
    rows = [rec]
    # MICA's small-data mode runs one unseeded k-means. scripts/mica_seeded_kmeans.py
    # reruns that exact call on MICA's saved embedding with fixed seeds, so MICA gets
    # the same three draws as the baselines.
    seeded = sorted((out_dir / "seeded").glob("seed_*.txt"))
    if seeded:
        rows = []
        for f in seeded:
            r = {**rec, "seed": int(f.stem.split("_")[1]), "mica_setting": f"k={info['n_clusters']} seeded"}
            r.update(metrics.score(truth.to_numpy(), aligned(pd.read_csv(f, sep="\t", index_col=0), f.name)))
            rows.append(r)
        rec = {**rec, "ari": sum(r["ari"] for r in rows) / len(rows)}
        if "published_ari" in rec and rec["published_ari"]:
            for r in rows:
                r["gap_vs_published"] = round(rec["ari"] - rec["published_ari"], 3)

    exp = exp_id or (f"reproduction_{name}" if tag is None else f"mica_{name}_{tag}")
    d = metrics.write(exp, rows, {"dataset": name, "mica_output": str(out_dir), "tag": tag,
                                  "seeded_draws": len(seeded), **info})
    print(f"{exp}: MICA ARI {rec['ari']:.3f}{' (mean of %d seeds)' % len(seeded) if seeded else ''} (published {rec.get('published_ari')}), "
          f"AMI {rec['ami']:.3f}, clusters {rec['n_clusters']}, setting {setting}")
    print(f"wrote {d / 'metrics.csv'}")
    return 0


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset"); ap.add_argument("out_dir"); ap.add_argument("tag", nargs="?")
    ap.add_argument("--exp", help="experiment id to write, for runs like graph mode on a small dataset")
    a = ap.parse_args()
    sys.exit(main(a.dataset, a.out_dir, a.tag or None, a.exp))
