"""Score MICA's later steps run on a given neighbor graph (U-19).

Usage: python scripts/score_downstream.py paper_pbmc14k <run_dir> <graph_label> <exp_id>
Reads clustering_<resolution>.txt from mica_downstream.py, picks by the paper's
rule, and writes results/<exp_id>/. The method is recorded as
mica_downstream so these rows never mix with the reproduction tables.
"""
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

from scminer_clust import cluster, metrics
from scminer_clust.paths import PROCESSED


def main(name, run_dir, graph_label, exp_id):
    run_dir = Path(run_dir)
    truth = ad.read_h5ad(PROCESSED / f"{name}.h5ad", backed="r").obs["label"]
    parts = {}
    for f in sorted(run_dir.glob("clustering_*.txt")):
        lab = pd.read_csv(f, sep="\t", index_col=0)["label"]
        assert set(lab.index.astype(str)) == set(truth.index), f"{f.name}: cells differ from {name}"
        parts[f.stem.split("_")[-1]] = lab.set_axis(lab.index.astype(str)).reindex(truth.index).astype(str).to_numpy()
    setting, info = cluster.pick_paper_protocol(parts, truth.to_numpy())
    rec = {"dataset": name, "method": "mica_downstream", "graph": graph_label, "k_setting": "paper_protocol",
           "mica_setting": setting, "n_settings_tried": len(parts), "seed": 0}
    rec.update(metrics.score(truth.to_numpy(), parts[setting]))
    secs = run_dir / "seconds.txt"
    if secs.exists():
        rec["seconds_downstream"] = float(secs.read_text().split("total")[1].split()[0])
    d = metrics.write(exp_id, [rec], {"run_dir": str(run_dir), "graph": graph_label, **info})
    print(f"{exp_id}: {graph_label}  ARI {rec['ari']:.3f}  AMI {rec['ami']:.3f}  clusters {rec['n_clusters']} at {setting}")
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:5]))
