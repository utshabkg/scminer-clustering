"""Collect every experiment's metrics into one sheet.

Usage: python scripts/update_tracker.py
Writes docs/results_tracker.csv, which is the file the team reads in meetings.
"""
import sys

import pandas as pd

from scminer_clust.paths import REPO, RESULTS

KEEP = ["experiment", "dataset", "method", "k_setting", "seed", "n_clusters",
        "ari", "ami", "nmi", "silhouette", "seconds", "peak_rss_gb"]


def main():
    frames = []
    for f in sorted(RESULTS.glob("*/metrics.csv")):
        df = pd.read_csv(f)
        df.insert(0, "experiment", f.parent.name)
        frames.append(df)
    if not frames:
        print("no results yet")
        return 0
    all_df = pd.concat(frames, ignore_index=True)
    cols = [c for c in KEEP if c in all_df.columns] + [c for c in all_df.columns if c not in KEEP]
    out = REPO / "docs" / "results_tracker.csv"
    all_df[cols].to_csv(out, index=False)
    print(f"wrote {out} with {len(all_df)} rows from {len(frames)} experiments")
    return 0


if __name__ == "__main__":
    sys.exit(main())
