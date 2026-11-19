"""The Allen atlas table and numbers for the paper, from results/large_scale/.

Usage: python scripts/large_scale/make_tables.py   (after run.sh <stage> full for every stage)

Writes results/large_scale/stage_costs.csv (time and peak memory per stage, copied from the gitignored
logs so the numbers stay in the repository), latex/tables/large_scale.tex and latex/numbers_large.tex.
Full-data results, results on the 100,000-cell sample and results on the sample with every other cell
assigned are kept in separate groups of the table.
"""
import json
import re

import pandas as pd

from scminer_clust.paths import REPO, RESULTS

OUT = RESULTS / "large_scale"
TABLES = REPO / "latex" / "tables"
# The 100,000-cell pilot of Oct 3, run whole in memory (docs/allen_pilot.md, logs/allen_pilot_*.log):
# MICA graph mode 25,727 s and 139.8 GB, the baselines' job 1,388 s and 26.6 GB with 2,000 genes.
PILOT = {"mica_s": 25727, "mica_gb": 139.8, "baselines_s": 1388, "baselines_gb": 26.6,
         "mica_extrap_days": "6 to 8", "mica_extrap_tb": 3.2}
STAGES = ["manifest", "scan", "hvg", "normalize", "pca", "kmeans", "leiden", "assign", "collect"]


def stage_costs():
    """Wall time and peak memory per stage. peak_tree_gb is peak_mem.sh's sum over the process and its
    workers, sampled every 10 s, an upper bound when forked workers share pages. peak_process_gb is the
    largest single process, from the stage's own record."""
    rec = {}
    for line in (OUT / "stages.jsonl").read_text().splitlines():
        d = json.loads(line)
        rec[d["stage"]] = d  # the last run of a stage wins
    rows = []
    for s in STAGES:
        if s not in rec:
            continue
        r = {"stage": s, "wall_s": rec[s]["wall_s"],
             "peak_process_gb": max(rec[s]["peak_rss_gb_self"], rec[s]["peak_rss_gb_child"])}
        log = REPO / "logs" / f"large_scale_full_{s}.log"
        if log.exists():
            txt = log.read_text()
            peaks = re.findall(r"PEAK_RSS_GB ([\d.]+)", txt)
            walls = re.findall(r"WALL_SECONDS (\d+)", txt)
            if peaks:
                r["peak_tree_gb"] = float(peaks[-1])
                r["wall_s_outer"] = int(walls[-1])
        rows.append(r)
    df = pd.DataFrame(rows)
    if "peak_tree_gb" not in df:
        df["peak_tree_gb"] = float("nan")
    df["peak_gb"] = df[["peak_process_gb", "peak_tree_gb"]].max(axis=1)
    df.to_csv(OUT / "stage_costs.csv", index=False)
    return df.set_index("stage")


def fmt_time(s):
    return f"{s / 3600:.1f} h" if s >= 5400 else f"{s / 60:.0f} min" if s >= 90 else f"{s:.0f} s"


def main():
    costs = stage_costs()
    m = pd.read_csv(OUT / "metrics.csv")
    dims = json.loads((OUT / "dims.json").read_text())
    pick = json.loads((OUT / "leiden_pick.json").read_text())
    assign = json.loads((OUT / "assign_info.json").read_text())
    search = pd.read_csv(OUT / "leiden_search.csv")
    agree = pd.read_csv(OUT / "seed_agreement.csv")
    hold = pd.read_csv(OUT / "assignment_holdout.csv").groupby("partition")[["accuracy", "ari"]].mean()
    kseconds = json.loads((REPO / "data" / "runs" / "large_scale" / "kmeans.done").read_text()).get("kmeans_seed0_s")
    prep = costs.loc[["scan", "hvg", "normalize", "pca"]]
    res = pick["picked_resolution"]
    leiden_s = float(search[(search.seed == 0) & ((search.resolution - res).abs() < 1e-9)].seconds.iloc[0])
    default_s = float(search[(search.seed == 0) & ((search.resolution - 1.0).abs() < 1e-9)].seconds.iloc[0])
    graph_s = pick["graph_seconds"]
    # the pilot baselines on the same 100,000 cells: time per method (k-means includes the gene selection
    # and PCA, Leiden the graph and its 200-resolution sweep), one peak for the whole job
    pilot = pd.read_csv(RESULTS / "baselines_allen_100k_hvg2000_seurat_v3" / "metrics.csv").set_index("method")

    def row(method, data, setting):
        r = m[(m.method == method) & (m.data == data) & (m.setting == setting)]
        assert len(r) == 1, (method, data, setting)
        return r.iloc[0]

    km = row("PCA + k-means", "full data", "k = 34")
    le = row("PCA + Leiden", "full data", f"resolution {res:.2f}")
    ld = row("PCA + Leiden", "full data", "resolution 1.00")
    groups = {
        "full": [("PCA + $k$-means", "$k = 34$", km, kseconds, costs.loc["kmeans", "peak_gb"]),
                 ("PCA + Leiden", f"resolution {res:.2f}", le, graph_s + leiden_s, costs.loc["leiden", "peak_gb"]),
                 ("PCA + Leiden", "resolution 1.0 (default)", ld, graph_s + default_s, costs.loc["leiden", "peak_gb"])],
        "sample": [], "assigned": []}
    for key, name in [("kmeans", "PCA + $k$-means"), ("leiden", "PCA + Leiden"), ("mica", "MICA, graph mode")]:
        setting = {"kmeans": "k = 34", "leiden": "pilot resolution", "mica": "paper rule"}[key]
        s, a = row(m_name(key), "100k sample", setting), row(m_name(key), "sample + assigned", setting)
        cost_s = PILOT["mica_s"] if key == "mica" else float(pilot.loc[f"{key}_pca", "seconds"])
        cost_gb = PILOT["mica_gb"] if key == "mica" else PILOT["baselines_gb"]
        groups["sample"].append((name, label_setting(key, assign, pilot), s, cost_s, cost_gb))
        groups["assigned"].append((name, label_setting(key, assign, pilot), a,
                                   assign["knn_seconds"] + assign["vote_seconds"][key], costs.loc["assign", "peak_gb"]))

    n_s = assign["n_reference"]
    lines = [r"% Written by scripts/large_scale/make_tables.py from results/large_scale/. Do not edit by hand.",
             r"\begin{table*}[t]", r"\centering",
             r"\caption{The Allen whole mouse brain atlas, " + f"{dims['n_cells']:,}" + r" cells by "
             + f"{dims['n_genes_kept']:,}" + r" genes after preprocessing. Full data: the method ran on every cell. "
             r"Sample: a uniform random " + f"{n_s:,}" + r"-cell sample, the only size at which MICA runs on our "
             r"machine, scored on those cells. Sample and assigned: every other cell takes the majority cluster of "
             r"its 15 nearest sample cells in the full-data principal components, scored on every cell. Each method "
             r"is set to the 34 classes, except Leiden's default. Time is the clustering itself, after the shared "
             r"preprocessing (" + fmt_time(prep.wall_s.sum()) + r", at most " + f"{prep.peak_gb.max():.0f}"
             + r"~GB). Sample times and memory are from the whole-job pilot runs, and the assignment rows add only "
             r"the neighbour search and vote.}",
             r"\label{tab:large}", r"\footnotesize", r"\setlength{\tabcolsep}{4pt}",
             r"\begin{tabular}{llrrrrr}", r"\toprule",
             r"Method & Setting & Clusters & ARI & AMI & Time & Peak memory \\", r"\midrule"]
    titles = {"full": f"\\multicolumn{{7}}{{l}}{{\\textit{{Full data, {dims['n_cells']:,} cells}}}} \\\\",
              "sample": f"\\multicolumn{{7}}{{l}}{{\\textit{{Sample, {n_s:,} cells}}}} \\\\",
              "assigned": f"\\multicolumn{{7}}{{l}}{{\\textit{{Sample and assigned, {dims['n_cells']:,} cells}}}} \\\\"}
    for g in ("full", "sample", "assigned"):
        if g != "full":
            lines.append(r"\midrule")
        lines.append(titles[g])
        for name, setting, r, secs, gb in groups[g]:
            t = "n/a" if secs is None else fmt_time(secs)
            lines.append(f"{name} & {setting} & {int(r.n_clusters)} & {r.ari:.3f} & {r.ami:.3f} & {t} & {gb:.0f}~GB \\\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table*}"]
    (TABLES / "large_scale.tex").write_text("\n".join(lines) + "\n")

    v = {"n_cells": f"{dims['n_cells']:,}", "n_genes_raw": f"{dims['n_genes_raw']:,}",
         "n_genes": f"{dims['n_genes_kept']:,}", "n_hvg": "2,000", "n_pcs": "50",
         "n_sample": f"{n_s:,}", "sample_pct": f"{100 * n_s / dims['n_cells']:.1f}",
         "n_assigned": f"{assign['n_assigned']:,}",
         "prep.minutes": f"{prep.wall_s.sum() / 60:.0f}", "prep.peak_gb": f"{prep.peak_gb.max():.0f}",
         "kmeans.ari": f"{km.ari:.3f}", "kmeans.ami": f"{km.ami:.3f}", "kmeans.minutes": f"{kseconds / 60:.1f}",
         "kmeans.peak_gb": f"{costs.loc['kmeans', 'peak_gb']:.0f}",
         "leiden.res": f"{res:.2f}", "leiden.k": str(int(le.n_clusters)), "leiden.ari": f"{le.ari:.3f}",
         "leiden.ami": f"{le.ami:.3f}", "leiden.minutes": f"{(graph_s + leiden_s) / 60:.0f}",
         "leiden.graph_minutes": f"{graph_s / 60:.0f}", "leiden.peak_gb": f"{costs.loc['leiden', 'peak_gb']:.0f}",
         "leiden.n_res": str(int((search.seed == 0).sum())),
         "leiden.default.k": str(int(ld.n_clusters)), "leiden.default.ari": f"{ld.ari:.3f}",
         "mica.hours": f"{PILOT['mica_s'] / 3600:.1f}", "mica.peak_gb": f"{PILOT['mica_gb']:.0f}",
         "mica.extrap_days": PILOT["mica_extrap_days"], "mica.extrap_tb": f"{PILOT['mica_extrap_tb']:.1f}",
         "assign.minutes": f"{assign['knn_seconds'] / 60:.0f}"}
    full_s = prep.wall_s.sum() + max(kseconds, graph_s + leiden_s)
    v["full.max_hours"] = f"{full_s / 3600:.0f}" if full_s >= 3600 else f"{full_s / 3600:.1f}"
    v["full.max_gb"] = f"{max(prep.peak_gb.max(), costs.loc['kmeans', 'peak_gb'], costs.loc['leiden', 'peak_gb']):.0f}"
    for key in ("kmeans", "leiden", "mica"):
        setting = {"kmeans": "k = 34", "leiden": "pilot resolution", "mica": "paper rule"}[key]
        s, a = row(m_name(key), "100k sample", setting), row(m_name(key), "sample + assigned", setting)
        v[f"sample.{key}.ari"], v[f"sample.{key}.k"] = f"{s.ari:.3f}", str(int(s.n_clusters))
        v[f"assigned.{key}.ari"] = f"{a.ari:.3f}"
        v[f"assigned.{key}.ari_assigned_cells"] = f"{a.ari_assigned_cells:.3f}"
        v[f"holdout.{key}.acc"] = f"{100 * hold.loc[key, 'accuracy']:.1f}"
    v["sample.mica.res"] = f"{assign['mica']['resolution']:.1f}"
    md = row("MICA (graph mode)", "100k sample", "resolution 1.8")
    v["sample.mica_default.ari"], v["sample.mica_default.k"] = f"{md.ari:.3f}", str(int(md.n_clusters))
    for meth, short in (("PCA + k-means", "kmeans"), ("PCA + Leiden", "leiden")):
        a = agree[agree.method == meth]
        if meth == "PCA + Leiden":
            a = a[a.setting == f"resolution {res:.2f}"]
        if len(a):
            v[f"{short}.seed_ari_min"] = f"{a.ari.min():.3f}"
    v["kmeans.vs_sample"] = f"{km.ari_sample_cells:.3f}"
    lines = [r"% Allen atlas numbers, written by scripts/large_scale/make_tables.py. Do not edit by hand."]
    lines += [f"\\expandafter\\def\\csname v@large.{k}\\endcsname{{{val}}}" for k, val in sorted(v.items())]
    (REPO / "latex" / "numbers_large.tex").write_text("\n".join(lines) + "\n")
    (OUT / "numbers_large.json").write_text(json.dumps(v, indent=1))
    print(json.dumps(v, indent=1))


def m_name(key):
    return {"kmeans": "PCA + k-means", "leiden": "PCA + Leiden", "mica": "MICA (graph mode)"}[key]


def label_setting(key, assign, pilot):
    if key == "mica":
        return f"resolution {assign['mica']['resolution']:.1f}"
    return "$k = 34$" if key == "kmeans" else f"resolution {pilot.loc['leiden_pca', 'resolution']:.2f}"


if __name__ == "__main__":
    main()
