"""Timings that exist only in the logs, saved where they can be checked.

Usage: python scripts/log_timings.py
logs/ is gitignored, so the paper's log-derived figures could not be verified from the
repository. This reads them once and writes results/a3_log_timings.csv:

  distance_matrix   MICA's one-off distance matrix per dataset, logs/rules_small.log
  matched_small     a matched small-dataset run, logs/matched_small_seed<s>.log
  matched_pbmc      a matched PBMC14K seed, logs/matched_pbmc_seed<s>.log
  downstream_pbmc   node2vec and the whole downstream step, logs/downstream_pbmc.log

For the two matched kinds, wall_s is the time between log timestamps and recorded_s the
selection seconds saved in that run's metrics.csv for both methods. The matched runs started
their timer after the full-data sweep, so unrecorded_s is that sweep plus data loading.
"""
import re
import sys
from datetime import datetime

import pandas as pd

from scminer_clust.paths import REPO, RESULTS

LOGS = REPO / "logs"


def stamp(text):
    return datetime.strptime(text.strip(), "%a %b %d %I:%M:%S %p %Z %Y")


def recorded(run):
    m = pd.read_csv(RESULTS / f"rule_matched_stability_paper_{run}" / "metrics.csv")
    return float(m[m.k_setting == "stability_finest"].seconds.sum())


def main():
    rows, name = [], None
    for line in (LOGS / "rules_small.log").read_text().splitlines():
        if line.startswith("######## "):
            name = line.split()[1]
        hit = re.search(r"MI distance for (\d+) cells in (\d+) s", line)
        if hit:
            rows.append({"kind": "distance_matrix", "run": name, "n_cells": int(hit[1]), "wall_s": int(hit[2])})
    for s in (0, 1, 2):
        lines = (LOGS / f"matched_small_seed{s}.log").read_text().splitlines()
        heads = [re.search(r"seed \d+ (\S+) (.*)$", x) for x in lines if x.startswith("######## matched stability")]
        times = [stamp(h[2]) for h in heads] + [stamp(re.search(r"DONE (.*)$", lines[-1])[1])]
        for h, t0, t1 in zip(heads, times, times[1:]):
            run = h[1].replace("paper_", "") + ("" if s == 0 else f"_seed{s}")
            rows.append({"kind": "matched_small", "run": run, "wall_s": (t1 - t0).total_seconds(), "recorded_s": recorded(run)})
        lines = (LOGS / f"matched_pbmc_seed{s}.log").read_text().splitlines()
        t0, t1 = stamp(re.search(r"seed \d+ (.*)$", lines[0])[1]), stamp(re.search(r"DONE (.*)$", lines[-1])[1])
        run = f"pbmc14k_graph_seed{s}"
        rows.append({"kind": "matched_pbmc", "run": run, "wall_s": (t1 - t0).total_seconds(), "recorded_s": recorded(run)})
    for line in (LOGS / "downstream_pbmc.log").read_text().splitlines():
        hit = re.search(r"node2vec (\d+) s, clustering done in (\d+) s, (\d+) resolutions, wrote \S*downstream_(\S+)", line)
        if hit:
            rows.append({"kind": "downstream_pbmc", "run": hit[4], "node2vec_s": int(hit[1]), "wall_s": int(hit[2]),
                         "n_resolutions": int(hit[3])})
    d = pd.DataFrame(rows)
    d["unrecorded_s"] = d.wall_s - d.recorded_s
    d.to_csv(RESULTS / "a3_log_timings.csv", index=False)
    print(d.groupby("kind").agg(n=("run", "size"), wall_min=("wall_s", "min"), wall_max=("wall_s", "max"),
                                unrecorded_min=("unrecorded_s", "min"), unrecorded_max=("unrecorded_s", "max")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
