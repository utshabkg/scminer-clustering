"""Candidate settings for the stability comparison, in one place.

The first stability runs capped both methods at 2 to 20 clusters, but a cap is not a match.
MICA's small-data k grid reaches every count from 2 to 20, while the baseline's Leiden grid
of 0.05 to 2.0 never went above 8 clusters on Yan or below 4 on Pollen, and MICA's graph-mode
Louvain grid stopped at 5 clusters on PBMC14K. The matched grids widen the baseline's grid at
both ends and refine MICA's graph-mode grid, so that both methods are offered the same
cluster counts where both can produce them. They were fixed from results/grid_coverage_*, which sweeps the full
data with no labels and no stability scores, before any matched run.

A runner uses the matched grids when --matched is on its command line.
"""
import sys

import numpy as np

LEIDEN_OLD = np.round(np.arange(0.05, 2.0 + 1e-9, 0.05), 2)
MICA_GRAPH_MAX_OLD = 4.0  # MICA's Louvain from 0.1 in steps of 0.1. Higher means fewer clusters.

# 300 log-spaced resolutions, the grid results/grid_coverage_fine probed. Even this skips up to
# 5 of the 19 counts in a run, which common_counts handles.
LEIDEN_MATCHED = np.round(np.geomspace(0.001, 50, 300), 5)
# MICA's graph-mode Louvain, from results/grid_coverage_mica_graph: up to 30 it never gives
# fewer than 4 clusters on PBMC14K or its 2,500 cell pools, and its 0.1 step skips counts
# between 10 and 20, which change fastest at low resolution. So 0.01 steps up to 2.0, then
# MICA's own 0.1 steps to 4.0, then 0.5 steps to 30.
MICA_GRAPH_MATCHED = np.round(np.concatenate([np.arange(0.1, 2.0, 0.01), np.arange(2.0, 4.0, 0.1),
                                              np.arange(4.0, 30.0 + 1e-9, 0.5)]), 2)

_matched = None


def matched():
    """True if --matched was given. Removes the flag so positional arguments still line up."""
    global _matched
    if _matched is None:
        _matched = "--matched" in sys.argv
        while "--matched" in sys.argv:
            sys.argv.remove("--matched")
    return _matched


def leiden():
    return LEIDEN_MATCHED if matched() else LEIDEN_OLD


def mica_graph_grid():
    """The matched Louvain list, or None for MICA's own 0.1 to 4.0 in steps of 0.1."""
    return MICA_GRAPH_MATCHED if matched() else None


def common_counts(*fulls):
    """Cluster counts from 2 to 20 that every method reaches on the full data.

    Each argument maps a setting to its full-data labels. A fixed grid cannot promise every
    count, because Leiden and Louvain jump over some as the resolution moves, so the matched
    comparison offers each method only the counts all of them can produce.
    """
    reached = [{len(np.unique(lab)) for lab in f.values()} for f in fulls]
    return sorted(set(range(2, 21)).intersection(*reached))


def one_per_count(full, counts):
    """One setting per cluster count, chosen on the full data with no labels.

    MICA's k grid has exactly one setting per count, while a dense resolution grid has
    several, and more candidates raise the rule's threshold (the best lower bound among
    them) by chance. So each method gets one: the middle of the settings, in grid order,
    whose full-data partition has that count. Returns {count: setting}.
    """
    by = {}
    for s in sorted(full):
        by.setdefault(len(np.unique(full[s])), []).append(s)
    return {k: by[k][(len(by[k]) - 1) // 2] for k in counts if k in by}


def tag(exp):
    """Result folder name: matched runs get their own prefix, so old globs never pick them up."""
    return exp.replace("rule_", "rule_matched_", 1) if matched() else exp
