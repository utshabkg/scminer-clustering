"""One table for the RQ1 diagnostics (U-14): co-detection and depth.

Usage: python scripts/summarize_diagnostics.py
Collects results/rq1_diagnostics_*/metrics.csv into results/summary_diagnostics.md.
"""
import sys

import pandas as pd

from scminer_clust.paths import RESULTS


def main():
    frames = [pd.read_csv(f) for f in sorted(RESULTS.glob("rq1_diagnostics_*/metrics.csv"))]
    if not frames:
        print("no diagnostics yet"); return 0
    d = pd.concat(frames, ignore_index=True)
    co = d[d.analysis == "detection_only"][["dataset", "codetection_share_median", "share_same_type", "share_diff_type",
                                            "detected_values_in_bin0", "ari", "ari_vs_full", "knn_overlap"]]
    co = co.merge(d[(d.analysis == "full") & (d.method == "mica")][["dataset", "ari"]].rename(columns={"ari": "ari_full_mica"}), on="dataset")
    co = co.rename(columns={"ari": "ari_detection_only", "ari_vs_full": "agreement_with_full_mica"}).round(3)
    depth = d[d.analysis.str.startswith("thin") | (d.analysis == "full")].copy()
    depth["depth"] = depth.analysis.map(lambda a: "100%" if a == "full" else f"{int(float(a.split('_')[1]) * 100)}%")
    dp = depth.pivot_table(index=["dataset", "method"], columns="depth", values="ari").round(3)
    dp = dp[[c for c in ("100%", "75%", "50%", "25%") if c in dp.columns]]
    agree = depth[depth.analysis != "full"].pivot_table(index=["dataset", "method"], columns="depth", values="ari_vs_full").round(3)
    md = ("# What MICA's similarity measures\n\n## Co-detection\n\nShare of the mutual information between two cells that "
          "comes only from which genes are detected in both (median over 20,000 random pairs), and how MICA does when it "
          "sees detection only.\n\n" + co.to_markdown(index=False)
          + "\n\n## Depth\n\nARI against the labels after thinning every cell's counts to a fraction of its reads.\n\n"
          + dp.to_markdown() + "\n\nAgreement (ARI) with the clusters found at full depth.\n\n" + agree.to_markdown() + "\n")
    (RESULTS / "summary_diagnostics.md").write_text(md)
    print(md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
