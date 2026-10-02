"""Score the authors' MICA output with the same rules as our baselines.

Usage: python scripts/score_mica.py zeisel results/reproduction/mica_zeisel

MICA writes one clustering_*.txt per setting (ID, X, Y, label). With a
resolution sweep there is one file per resolution, and we pick among them with
the paper's protocol, which is what the published numbers used.
"""
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

from scminer_clust import cluster, metrics
from scminer_clust.paths import PROCESSED

PUBLISHED_ARI = {"zeisel": 0.744, "pbmc14k": 0.853}  # scMINER Source Data, per dataset


def main(name, out_dir):
    out_dir = Path(out_dir)
    truth = ad.read_h5ad(PROCESSED / f"{name}.h5ad", backed="r").obs["label"]
    files = sorted(out_dir.glob("clustering_*.txt"))
    if not files:
        print(f"no clustering_*.txt in {out_dir}")
        return 1

    partitions = {}
    for f in files:
        df = pd.read_csv(f, sep="\t", index_col=0)
        setting = f.stem.split("_")[-1]
        partitions[setting] = df["label"].astype(str).reindex(truth.index).to_numpy()

    setting, info = cluster.pick_paper_protocol(partitions, truth.to_numpy())
    rec = {"dataset": name, "method": "mica", "seed": 0, "k_setting": "paper_protocol",
           "mica_setting": setting, "n_settings_tried": len(partitions)}
    rec.update(metrics.score(truth.to_numpy(), partitions[setting]))
    rec["published_ari"] = PUBLISHED_ARI.get(name)
    rec["gap_vs_published"] = round(rec["ari"] - rec["published_ari"], 3) if rec["published_ari"] else None

    d = metrics.write(f"reproduction_{name}", rec, {"dataset": name, "mica_output": str(out_dir), **info})
    print(f"{name}: MICA ARI {rec['ari']:.3f} (published {rec['published_ari']}), "
          f"AMI {rec['ami']:.3f}, clusters {rec['n_clusters']}, setting {setting}")
    print(f"wrote {d / 'metrics.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
