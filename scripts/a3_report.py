"""Tables and figures for the Assignment 3 manuscript, built from results/.

Usage: python scripts/a3_report.py
Writes latex/tables/a3_*.tex, latex/figures/a3_*.pdf and
results/a3_report_numbers.json, which holds every number the text quotes so it can
be checked. Nothing here reruns an experiment.
"""
import json
import sys
from glob import glob
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scminer_clust.paths import REPO, RESULTS

TABLES, FIGS = REPO / "latex" / "tables", REPO / "latex" / "figures"
ORDER = ["yan", "goolam", "buettner", "pollen", "chung", "usoskin", "kolodziejczyk", "klein", "zeisel", "pbmc14k"]
NAMES = {"yan": "Yan", "goolam": "Goolam", "buettner": "Buettner", "pollen": "Pollen", "chung": "Chung",
         "usoskin": "Usoskin", "kolodziejczyk": "Kolodziejczyk", "klein": "Klein", "zeisel": "Zeisel",
         "pbmc14k": "PBMC14K"}
RULES = [("default", "Default"), ("silhouette_shared", "Silhouette"),
         ("stability_finest", "Finest stable"), ("stability_most", "Most stable")]
# PBMC14K runs MICA's graph mode and is reported on its own. Every mean, win count and
# mechanism figure in the text is over the nine small-data datasets only, so adding it
# here must not change them.
SMALL = [x for x in ORDER if x != "pbmc14k"]
BLUE, ORANGE, INK, MUTED = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"


def short(ds):
    return ds.replace("paper_", "")


def load_stability():
    rows = []
    patterns = ["rule_stability_*", "rule_silhouette_*", "rule_default_*"]
    for f in sorted(f for pat in patterns for f in glob(str(RESULTS / pat / "metrics.csv"))):
        d = pd.read_csv(f)
        if "seed" not in d:
            d["seed"] = 0
        rows.append(d)
    if not rows:
        return pd.DataFrame()
    d = pd.concat(rows, ignore_index=True)
    # Leiden's shipped default was scored in the Assignment 2 runs and lives in the tracker
    t = pd.read_csv(REPO / "docs" / "results_tracker.csv", low_memory=False)
    lei = t[(t.k_setting == "default_resolution_1") & (t.gene_setting == "all genes")
            & t.dataset.str.startswith("paper_") & (t.method == "leiden_pca")].copy()
    lei["k_setting"] = "default"
    # MICA's graph-mode default (resolution 1.822, scored at 1.8) on PBMC14K, the one dataset
    # where MICA runs in graph mode. Its small-data default lives in results/rule_default_*.
    mge = t[(t.k_setting == "default_resolution_1.8") & (t.gene_setting == "all genes")
            & (t.dataset == "paper_pbmc14k") & (t.method == "mica")].copy()
    mge["k_setting"] = "default"
    d = pd.concat([d, lei, mge], ignore_index=True)
    d["ds"] = d.dataset.map(short)
    return d


def agg(d, rule):
    """Median and range of ARI and chosen cluster count over the repeats of a rule."""
    s = d[d.k_setting == rule]
    if not len(s):
        return pd.DataFrame()
    g = s.groupby(["ds", "method"]).agg(
        ari=("ari", "median"), ari_min=("ari", "min"), ari_max=("ari", "max"),
        ami=("ami", "median"), k=("n_clusters", "median"), k_min=("n_clusters", "min"),
        k_max=("n_clusters", "max"), runs=("ari", "size"), seconds=("seconds", "median")).reset_index()
    return g


def cell(r):
    if r is None or not len(r):
        return "n/a"
    v = r.iloc[0]
    if v.runs > 1 and v.ari_max - v.ari_min >= 0.005:
        return f"{v.ari:.3f} [{v.ari_min:.2f}, {v.ari_max:.2f}]"
    return f"{v.ari:.3f}"


def kcell(r, truth):
    if r is None or not len(r):
        return "n/a"
    v = r.iloc[0]
    k = f"{int(v.k)}" if v.k_min == v.k_max else f"{int(v.k_min)} to {int(v.k_max)}"
    return k


def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    nums = {}
    facts = pd.read_csv(REPO / "docs" / "dataset_facts.csv").set_index("dataset")
    d = load_stability()
    if not len(d):
        print("no stability results yet")
        return 1
    present = [x for x in ORDER if x in set(d.ds)]
    truth_k = {x: int(facts.loc[f"paper_{x}", "label_classes"]) for x in present}

    # ---------- dataset table ----------
    lines = []
    for x in ORDER:
        f = facts.loc[f"paper_{x}"]
        lines.append(f"{NAMES[x]} & {int(f.cells_after_filtering):,} & {int(f.genes_after_filtering):,} & "
                     f"{int(f.label_classes)} & {f.smallest_class_percent:.1f} & {f.percent_zeros:.0f} \\\\")
    (TABLES / "a3_datasets.tex").write_text(
        "\\begin{table}[t]\n\\centering\n\\caption{The ten benchmark datasets, as released by the scMINER authors "
        "and used here without further preprocessing. The smallest class column is that class as a share of all "
        "cells, and the last column is the share of the matrix that is zero.}\n\\label{tab:data}\n\\footnotesize\n"
        "\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{lrrrrr}\n\\toprule\nDataset & Cells & Genes & Types & "
        "Smallest (\\%) & Zeros (\\%) \\\\\n\\midrule\n" + "\n".join(lines) +
        "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n")

    # ---------- rule tables ----------
    # The body carries the two stability readings, which is the argument the paper makes and
    # the only place the chosen cluster count matters. All four rules go to the appendix.
    tabs = {r: agg(d, r) for r, _ in RULES}
    have = [(r, lab) for r, lab in RULES if len(tabs[r])]

    def rule_table(rules, label, caption, star):
        rows = []
        for x in present:
            cells = []
            for r, _ in rules:
                t = tabs[r]
                for m in ("mica", "leiden_pca"):
                    sel = t[(t.ds == x) & (t.method == m)] if len(t) else None
                    if sel is None or not len(sel):
                        cells.append("n/a")
                        continue
                    v = sel.iloc[0]
                    k = f"{int(v.k)}" if v.k_min == v.k_max else f"{int(v.k_min)}\\textrm{{-}}{int(v.k_max)}"
                    cells.append(f"{v.ari:.3f}\\,({k})")
            rows.append(f"{NAMES[x]} ({truth_k[x]}) & " + " & ".join(cells) + " \\\\")
        means = []
        for r, _ in rules:
            t = tabs[r]
            for m in ("mica", "leiden_pca"):
                sel = t[(t.ds.isin(SMALL)) & (t.method == m)]
                means.append(f"\\textbf{{{sel.ari.mean():.3f}}}" if len(sel) else "n/a")
        rows.append("\\midrule\n\\textit{Mean, nine small} & " + " & ".join(means) + " \\\\")
        g = " & ".join(f"\\multicolumn{{2}}{{c}}{{{lab}}}" for _, lab in rules)
        cm = "".join(f"\\cmidrule(lr){{{2 + 2 * i}-{3 + 2 * i}}}" for i in range(len(rules)))
        hd = " & ".join(["MICA & Leiden"] * len(rules))
        env = "table*" if star else "table"
        # A single-column table is scaled to the column width so it can never overflow it.
        open_box, close_box = ("", "") if star else ("\\resizebox{\\columnwidth}{!}{%\n", "}\n")
        return (f"\\begin{{{env}}}[t]\n\\centering\n\\caption{{{caption}}}\n\\label{{{label}}}\n"
                "\\footnotesize\n\\setlength{\\tabcolsep}{4pt}\n" + open_box +
                f"\\begin{{tabular}}{{l{'cc' * len(rules)}}}\n\\toprule\n & {g} \\\\\n{cm}\n"
                f"Dataset (types) & {hd} \\\\\n\\midrule\n" + "\n".join(rows) +
                f"\n\\bottomrule\n\\end{{tabular}}\n" + close_box + f"\\end{{{env}}}\n")

    stab = [(r, lab) for r, lab in have if r.startswith("stability")]
    (TABLES / "a3_rules.tex").write_text(rule_table(
        stab, "tab:rules",
        "ARI of the partition the stability rule chose under each of its two readings, with "
        "the number of clusters chosen in parentheses and the true number of cell types after "
        "each dataset name. No reference labels are used in the choice. Values are medians "
        "over three seeds, and a pair such as (15-18) gives the smallest and largest count "
        "over those seeds. The mean row covers the nine small-data datasets, the set every "
        "mean in the text refers to. PBMC14K runs in MICA\'s graph mode and is discussed in "
        "Section~\\ref{sec:pbmc}. The other two rules are in Table~\\ref{tab:rulesall}, and "
        "Table~\\ref{tab:common} gives this rule on the cluster counts the baseline produces.",
        star=False))
    (TABLES / "a3_rules_all.tex").write_text(rule_table(
        have, "tab:rulesall",
        "ARI under all four label-free rules, with the chosen cluster count in parentheses. "
        "The two stability columns repeat Table~\\ref{tab:rules} and are included so the four "
        "rules can be compared in one place. The default rule was run at three seeds and the "
        "silhouette rule at one. On PBMC14K MICA runs in graph mode, so its default there is "
        "the graph-mode resolution, 1.822 scored at 1.8, from one Assignment 2 run, and the "
        "silhouette rule was not run.",
        star=True))

    # ---------- AMI, promised in the Measures subsection ----------
    arows = []
    for x in present:
        cells = []
        for r, _ in have:
            t = tabs[r]
            for m in ("mica", "leiden_pca"):
                sel = t[(t.ds == x) & (t.method == m)]
                cells.append(f"{sel.iloc[0].ami:.3f}" if len(sel) and sel.iloc[0].ami == sel.iloc[0].ami
                             else "n/a")
        arows.append(f"{NAMES[x]} ({truth_k[x]}) & " + " & ".join(cells) + " \\\\")
    ameans = []
    for r, _ in have:
        t = tabs[r]
        for m in ("mica", "leiden_pca"):
            sel = t[(t.ds.isin(SMALL)) & (t.method == m)].dropna(subset=["ami"])
            ameans.append(f"\\textbf{{{sel.ami.mean():.3f}}}" if len(sel) else "n/a")
    arows.append("\\midrule\n\\textit{Mean, nine small} & " + " & ".join(ameans) + " \\\\")
    groups = " & ".join(f"\\multicolumn{{2}}{{c}}{{{lab}}}" for _, lab in have)
    cmid = "".join(f"\\cmidrule(lr){{{2 + 2 * i}-{3 + 2 * i}}}" for i in range(len(have)))
    heads = " & ".join(["MICA & Leiden"] * len(have))
    (TABLES / "a3_ami.tex").write_text(
        "\\begin{table*}[t]\n\\centering\n\\caption{Adjusted mutual information for the same chosen partitions "
        "as Table~\\ref{tab:rulesall}, on the same runs. The AMI and the ARI agree on which method leads under "
        "every rule. The AMI sits above the ARI whenever the chosen count is wrong, in either direction, and it is "
        "much more forgiving of a partition that merges cell types, such as MICA\'s 2-cluster pick on "
        "Pollen.}\n\\label{tab:ami}\n\\footnotesize\n\\setlength{\\tabcolsep}{4pt}\n"
        f"\\begin{{tabular}}{{l{'cc' * len(have)}}}\n\\toprule\n & {groups} \\\\\n{cmid}\nDataset (types) & "
        f"{heads} \\\\\n\\midrule\n" + "\n".join(arows) + "\n\\bottomrule\n\\end{tabular}\n\\end{table*}\n")

    # ---------- per-dataset seed variation, which the main table has no room for ----------
    sv = d[d.k_setting.isin(["stability_finest", "stability_most"])].groupby(
        ["ds", "method", "k_setting"]).agg(
        runs=("ari", "size"), ari_median=("ari", "median"), ari_min=("ari", "min"),
        ari_max=("ari", "max"), k_median=("n_clusters", "median"), k_min=("n_clusters", "min"),
        k_max=("n_clusters", "max")).round(4).reset_index()
    sv["ari_range"] = (sv.ari_max - sv.ari_min).round(4)
    sv.to_csv(RESULTS / "a3_seed_variation.csv", index=False)
    print(f"seed variation: largest ARI range {sv.ari_range.max():.3f} "
          f"({sv.loc[sv.ari_range.idxmax(), 'ds']}, {sv.loc[sv.ari_range.idxmax(), 'method']})")

    # ---------- the same rule on the cluster counts both methods can reach ----------
    # The 2 to 20 filter caps the candidates but does not match them: MICA's k covers 2 to 20
    # on every dataset, while the baseline's resolution grid reaches only part of that range.
    # Re-applying the rule's own selection (stability.pick) to MICA's saved curve, restricted to
    # the set of cluster counts the baseline's partitions produced in the same run, separates the
    # method from the range. It is the set, not the interval between the fewest and most: on Yan
    # the baseline makes 2 and 4 to 8 clusters but never 3. Restricting uses no labels, since it
    # depends only on which counts each grid produced. It is done after the results were seen.
    def repick(c):
        c = c[c.n_clusters >= 2]
        thr = c.ci_low.max()
        pas = c[c.stability >= thr]
        pas = pas if len(pas) else c
        fin = pas.sort_values(["n_clusters", "stability"]).iloc[-1]
        most = c.sort_values(["stability", "n_clusters"], ascending=[False, True]).iloc[0]
        tied = c[c.stability >= c.stability.max() - TIE_WINDOW]
        return fin, most, int(tied.n_clusters.max() - tied.n_clusters.min())

    TIE_WINDOW = 0.01
    crows = []
    for f in sorted(glob(str(RESULTS / "rule_stability_paper_*" / "metrics.csv"))):
        if "graph" in f:
            continue
        run = Path(f).parent.name.replace("rule_stability_paper_", "")
        x = run.split("_seed")[0]
        cur = pd.read_csv(f)
        cur = cur[cur.k_setting == "stability_curve"]
        mi, le = cur[cur.method == "mica"], cur[cur.method == "leiden_pca"]
        lo, hi = int(le.n_clusters.min()), int(le.n_clusters.max())
        made = set(le.n_clusters.astype(int))
        for arm, c in [("mica_own", mi), ("mica_common", mi[mi.n_clusters.isin(made)]), ("leiden", le)]:
            fin, most, span = repick(c)
            avail = set(c.n_clusters.astype(int))
            crows.append({"ds": x, "run": run, "arm": arm, "reach_lo": lo, "reach_hi": hi,
                          "baseline_counts": " ".join(map(str, sorted(made))),
                          "finest_ari": fin.ari, "finest_k": int(fin.n_clusters),
                          "most_ari": most.ari, "most_k": int(most.n_clusters), "tied_span": span,
                          # is the finest-stable pick at an end of what this arm could choose from?
                          "finest_at_edge": int(fin.n_clusters) in (min(avail), max(avail))})
    C = pd.DataFrame(crows)
    C.to_csv(RESULTS / "a3_common_range.csv", index=False)
    cm = C.groupby(["ds", "arm"])[["finest_ari", "most_ari", "tied_span"]].median().reset_index()
    cmean = cm.groupby("arm")[["finest_ari", "most_ari", "tied_span"]].mean()
    reach = C[C.arm == "leiden"].groupby("ds").agg(lo=("reach_lo", "min"), hi=("reach_hi", "max"))
    cm["abs_gap"] = (cm.finest_ari - cm.most_ari).abs()
    abs_gap = cm.groupby("arm").abs_gap.mean()
    edges = C.groupby("arm").finest_at_edge.agg(["sum", "size"])
    nums["common_range"] = {
        "restriction": "set of cluster counts the baseline produced in the same run",
        "mean_abs_reading_change": {a: round(float(v), 4) for a, v in abs_gap.items()},
        "finest_pick_at_edge": {a: f"{int(r['sum'])} of {int(r['size'])} runs" for a, r in edges.iterrows()},
        "means": {a: {"finest": round(float(r.finest_ari), 4), "most": round(float(r.most_ari), 4),
                      "gap": round(float(r.finest_ari - r.most_ari), 4), "tied_span": round(float(r.tied_span), 2)}
                  for a, r in cmean.iterrows()},
        "baseline_reach": {k: [int(v.lo), int(v.hi)] for k, v in reach.iterrows()}}
    # The span quoted in the text pools the three seeds before taking the window (see the
    # mechanism block), so the common-range span is computed the same way to be comparable.
    made_any = {x: set().union(*[set(map(int, b.split())) for b in g_.baseline_counts])
                for x, g_ in C[C.arm == "leiden"].groupby("ds")}
    pooled = {}
    for (x, m_), g in d[(d.k_setting == "stability_curve") & d.ds.isin(SMALL)].dropna(subset=["stability"]).groupby(["ds", "method"]):
        if m_ == "mica":
            g = g[g.n_clusters.isin(made_any[x])]
        t_ = g[g.stability >= g.stability.max() - TIE_WINDOW]
        pooled.setdefault(m_, []).append(int(t_.n_clusters.max() - t_.n_clusters.min()))
    nums["common_range"]["pooled_tied_span"] = {m_: round(float(np.mean(v)), 2) for m_, v in pooled.items()}
    crow = []
    for x in [d_ for d_ in ORDER if d_ in set(cm.ds)]:
        g = cm[cm.ds == x].set_index("arm")
        r_ = reach.loc[x]
        crow.append(f"{NAMES[x]} ({truth_k[x]}) & {int(r_.lo)}\\textrm{{-}}{int(r_.hi)} & "
                    f"{g.loc['mica_own'].finest_ari:.3f} & {g.loc['mica_own'].most_ari:.3f} & "
                    f"{g.loc['mica_common'].finest_ari:.3f} & {g.loc['mica_common'].most_ari:.3f} & "
                    f"{g.loc['leiden'].finest_ari:.3f} & {g.loc['leiden'].most_ari:.3f} \\\\")
    mrow = " & ".join(f"\\textbf{{{cmean.loc[a][col]:.3f}}}" for a in ("mica_own", "mica_common", "leiden")
                      for col in ("finest_ari", "most_ari"))
    crow.append("\\midrule\n\\textit{Mean} & & " + mrow + " \\\\")
    (TABLES / "a3_common_range.tex").write_text(
        "\\begin{table*}[t]\n\\centering\n\\caption{The stability rule on its own grids and on a common "
        "range. The second column is the span of cluster counts the baseline's resolution grid produced over "
        "the three seeds. MICA's $k$ covers 2 to 20 on every dataset. The middle pair re-applies the rule's own "
        "selection to MICA's saved scores restricted to the set of counts the baseline produced in the same run, "
        "which can have gaps within that span. The restriction uses no labels but was chosen after the results "
        "were seen and applies to MICA only. Values are medians over three seeds, F the finest-stable reading and "
        "M the plain maximum.}\n\\label{tab:common}\n\\footnotesize\n\\setlength{\\tabcolsep}{4pt}\n"
        "\\begin{tabular}{lccccccc}\n\\toprule\n & & \\multicolumn{2}{c}{MICA, own grid} & "
        "\\multicolumn{2}{c}{MICA, common range} & \\multicolumn{2}{c}{Leiden} \\\\\n"
        "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}\n"
        "Dataset (types) & Baseline span & F & M & F & M & F & M \\\\\n\\midrule\n"
        + "\n".join(crow) + "\n\\bottomrule\n\\end{tabular}\n\\end{table*}\n")
    print("common range, pooled tied span: " + str(nums["common_range"]["pooled_tied_span"]))
    print("common range: " + ", ".join(f"{a} F {v['finest']:.3f} M {v['most']:.3f} span {v['tied_span']}"
                                       for a, v in nums["common_range"]["means"].items()))

    # the mode test on a common range: the small-data arm held to the counts the graph arm reaches
    mrange = []
    for run in sorted(glob(str(RESULTS / "rule_mode_size_pbmc14k_n*_seed*"))):
        m = pd.read_csv(Path(run) / "metrics.csv")
        cur = m[m.k_setting == "stability_curve"]
        sm = cur[(cur.method == "mica") & (cur.mica_mode == "small")]
        gr = cur[(cur.method == "mica") & (cur.mica_mode == "graph")]
        lo, hi = int(gr.n_clusters.min()), int(gr.n_clusters.max())
        for arm, c in [("small_own", sm), ("small_common", sm[sm.n_clusters.between(lo, hi)]), ("graph", gr)]:
            fin, most, span = repick(c)
            mrange.append({"seed": int(run.split("seed")[-1]), "arm": arm, "graph_lo": lo, "graph_hi": hi,
                           "gap": round(float(fin.ari - most.ari), 4), "tied_span": span})
    M = pd.DataFrame(mrange)
    nums["mode_size_common_range"] = {
        a: {"readings_disagree": f"{int((g.gap.abs() > 0.02).sum())} of {len(g)}",
            "gaps": g.gap.tolist(), "spans": g.tied_span.tolist(),
            "graph_reach": [[int(r.graph_lo), int(r.graph_hi)] for _, r in g.iterrows()]}
        for a, g in M.groupby("arm")}
    print("mode test on a common range: " + ", ".join(
        f"{a} disagree {v['readings_disagree']} spans {v['spans']}" for a, v in nums["mode_size_common_range"].items()))

    # ---------- mode against size, one fixed pool of cells run through both pipelines ----------
    pool = sorted(glob(str(RESULTS / "rule_mode_size_pbmc14k_n*_seed*")))
    if pool:
        trows, prows = [], []
        for run in pool:
            seed = int(run.split("seed")[-1])
            t = pd.read_csv(Path(run) / "tie_summary.csv"); t["seed"] = seed
            trows.append(t)
            m = pd.read_csv(Path(run) / "metrics.csv")
            m = m[m.k_setting.isin(["stability_finest", "stability_most"])]
            for _, r in m.iterrows():
                prows.append({"seed": seed, "arm": f"{r.method}/{r.mica_mode}",
                              "rule": r.k_setting, "k": int(r.n_clusters), "ari": r.ari})
        T, P = pd.concat(trows), pd.DataFrame(prows)
        w = P.pivot_table(index=["arm", "seed"], columns="rule", values="ari")
        w["gap"] = w.stability_finest - w.stability_most
        ARMS = [("mica", "small", "MICA, small-data mode"), ("mica", "graph", "MICA, graph mode"),
                ("leiden_pca", "pca", "Leiden, as a control")]
        body, ms = [], {}
        for meth, mode, label in ARMS:
            g = T[(T.method == meth) & (T["mode"] == mode)]
            if not len(g):
                continue
            arm = f"{meth}/{mode}"
            gaps = w.loc[arm, "gap"] if arm in w.index.get_level_values(0) else pd.Series(dtype=float)
            span = f"{int(g.k_tied_span.median())} ({int(g.k_tied_span.min())}\\textrm{{-}}{int(g.k_tied_span.max())})"
            dis = int((gaps > 0.001).sum())
            body.append(f"{label} & {g.share_tied.median():.2f} & {span} & "
                        f"{dis} of {len(gaps)} & {gaps.median():.3f} \\\\")
            ms[arm] = {"share_tied_median": round(float(g.share_tied.median()), 3),
                       "k_tied_span_median": int(g.k_tied_span.median()),
                       "k_tied_span_range": [int(g.k_tied_span.min()), int(g.k_tied_span.max())],
                       "readings_disagree": f"{dis} of {len(gaps)}",
                       "gap_median": round(float(gaps.median()), 3),
                       "gaps": [round(float(x), 3) for x in gaps]}
        (TABLES / "a3_mode_size.tex").write_text(
            "\\begin{table}[t]\n\\centering\n\\caption{Mode against size. One fixed pool of 2,500 PBMC14K "
            "cells, drawn afresh at each of three seeds, run through both of MICA\'s pipelines and the "
            "baseline under the same rule. The data, the size, the subsamples and the labels are identical "
            "across the three rows. The pipelines differ, and so do the cluster counts their grids reach: "
            "graph mode cannot make fewer than 4 or 5 clusters here, and Section~\\ref{sec:modesize} shows "
            "that this, not the mode, accounts for the difference. Tied settings are those within 0.01 of the "
            "best. The span is the range of cluster counts among them, reported as the median over the three "
            "pools with the observed range. The last two columns say how often the two readings of the rule "
            "chose differently and by how much in ARI.}\n\\label{tab:modesize}\n\\footnotesize\n"
            "\\setlength{\\tabcolsep}{4pt}\n\\resizebox{\\columnwidth}{!}{%\n\\begin{tabular}{lcccc}\n\\toprule\n"
            "Pipeline & Share tied & Span & Readings differ & Median gap \\\\\n\\midrule\n"
            + "\n".join(body) + "\n\\bottomrule\n\\end{tabular}}\n\\end{table}\n")
        nums["mode_size"] = {"n_pool": 2500, "seeds": sorted(P.seed.unique().tolist()), "arms": ms}
        print("mode against size: " + ", ".join(
            f"{k} span {v['k_tied_span_median']} disagree {v['readings_disagree']}" for k, v in ms.items()))

    # ---------- graph quality and cost: correlation shortlist against MICA's own HNSW ----------
    # Graph build time and memory come from scripts/measure_graph_cost.py and measure_mica_hnsw.py,
    # recall from the brute-force audits, and the final clustering from MICA's own node2vec and
    # Louvain run on each graph (results/rq2_downstream_*). The downstream ARIs let the labels pick
    # the resolution, because this compares graphs, not selection rules.
    GRAPHS = [("paper_pbmc14k", "PBMC14K", 800, "rq2_large_paper_pbmc14k"),
              ("allen_100k", "Allen 100k", 3200, "rq2_allen_100k")]
    grows, gnums = [], {}
    for ds, label, c, audit in GRAPHS:
        aud = pd.read_csv(RESULTS / audit / "metrics.csv")
        mica_recall = float(aud[aud.candidates.astype(str).str.startswith("MICA HNSW")].recall_at_k.iloc[0])
        for kind, gtag, cost_dir in [("MICA HNSW", "mica_hnsw", f"graph_cost_{ds}_mica_hnsw"),
                                     (f"Shortlist, $c={c}$", f"corr{c}", f"graph_cost_{ds}_corr{c}")]:
            f = RESULTS / cost_dir / "metrics.csv"
            cost = pd.read_csv(f).iloc[0] if f.exists() else None
            down = []
            for rep in (1, 2):
                df = RESULTS / f"rq2_downstream_{ds}_{gtag}_rep{rep}" / "metrics.csv"
                if df.exists():
                    down.append(pd.read_csv(df).iloc[0])
            # MICA's recall comes from the Assignment 2 audit and does not need the cost run.
            if gtag == "mica_hnsw":
                recall = mica_recall
            elif cost is not None:
                recall = float(cost.recall_500_queries)
            else:
                recall = None
            # MICA's search takes the raw matrix and bins inside its distance, while the shortlist bins
            # once up front, so both builds count loading (and binning) to cover the same work.
            g_s = float(cost.seconds_load_and_bin) + float(cost.seconds_graph) if cost is not None else None
            mem = float(cost.peak_rss_gb) if cost is not None else None
            d_s = float(np.median([x.seconds_downstream for x in down])) if down else None
            aris = [x.ari for x in down]
            amis = [x.ami for x in down]

            def fmt_s(v):
                if v is None:
                    return "pending"
                return f"{v / 3600:.1f} h" if v >= 3600 else (f"{v / 60:.1f} min" if v >= 120 else f"{v:.0f} s")

            cand = (f" ({fmt_s(float(cost.seconds_candidates))})"
                    if cost is not None and str(cost.seconds_candidates) not in ("", "nan") else "")
            grows.append(f"{label} & {kind} & {fmt_s(g_s)}{cand} & "
                         f"{'pending' if mem is None else f'{mem:.1f}'} & "
                         f"{'pending' if recall is None else f'{recall:.3f}'} & "
                         f"{fmt_s(None if g_s is None or d_s is None else g_s + d_s)} & "
                         f"{' / '.join(f'{a:.3f}' for a in aris) or 'n/a'} & "
                         f"{' / '.join(f'{a:.3f}' for a in amis) or 'n/a'} \\\\")
            gnums[f"{ds}/{gtag}"] = {"seconds_build_with_load": g_s, "peak_rss_gb": mem, "recall": recall,
                                     "seconds_downstream_median": d_s, "ari": aris, "ami": amis,
                                     "seconds_candidates": (None if cost is None or
                                                            str(cost.seconds_candidates) in ("", "nan")
                                                            else float(cost.seconds_candidates))}
        if label == "PBMC14K":
            grows.append("\\midrule")
    (TABLES / "a3_graph_cost.tex").write_text(
        "\\begin{table*}[t]\n\\centering\n\\caption{Graph quality and cost. MICA's own HNSW search against "
        "a correlation shortlist of $c$ candidates per cell rescored with the exact MICA distance, which keeps "
        "each cell's 79 closest. Nothing guarantees the shortlist holds every true neighbour, so recall is "
        "measured against brute force on 500 query cells. Graph build times are on 8 workers and include loading "
        "the data and, for the shortlist, binning it once, which MICA's search does inside its distance. The cost "
        "of generating the shortlist is in parentheses. Peak memory covers the same steps. For MICA, recall, ARI "
        "and AMI come from the graph it built in Assignment 2, and time and memory from a rebuild with the same "
        "settings. "
        "Total adds MICA's node2vec and Louvain, which run unchanged on either graph. ARI and AMI are from two "
        "runs of that downstream pipeline on each graph. They let the labels pick the resolution, as in the "
        "published protocol, because this table compares graphs and not selection rules.}\n"
        "\\label{tab:graph}\n\\footnotesize\n\\setlength{\\tabcolsep}{4pt}\n"
        "\\begin{tabular}{llcccccc}\n\\toprule\n"
        "Data & Graph & Build (shortlist) & Peak GB & Recall@79 & Total & ARI, two runs & AMI, two runs "
        "\\\\\n\\midrule\n" + "\n".join(grows) + "\n\\bottomrule\n\\end{tabular}\n\\end{table*}\n")
    nums["graph_cost"] = gnums
    print("graph cost: " + "; ".join(f"{k} build {v['seconds_build_with_load']} s, {v['peak_rss_gb']} GB, recall "
                                      f"{v['recall']}" for k, v in gnums.items()))

    # ---------- numbers for the text ----------
    for r, _ in RULES:
        t = tabs[r]
        if not len(t):
            continue
        piv = t.pivot_table(index="ds", columns="method", values="ari")
        if {"mica", "leiden_pca"} <= set(piv.columns):
            piv = piv.dropna()
            piv = piv[piv.index.isin(SMALL)]
            nums[r] = {"per_dataset": piv.round(4).to_dict(),
                       "mica_mean": round(float(piv.mica.mean()), 4),
                       "leiden_mean": round(float(piv.leiden_pca.mean()), 4),
                       "mica_wins": int((piv.mica > piv.leiden_pca).sum()), "n_datasets": int(len(piv)),
                       "median_gap": round(float((piv.mica - piv.leiden_pca).median()), 4)}
        kp = t.pivot_table(index="ds", columns="method", values="k")
        if len(kp):
            nums.setdefault(r, {})["chosen_k"] = kp.round(1).to_dict()
        nums.setdefault(r, {})["seconds_median"] = t.groupby("method").seconds.median().round(1).to_dict()
    # ---------- mechanism: can the rule tell the settings apart at all? ----------
    from scipy.stats import spearmanr
    curve_all = d[d.k_setting == "stability_curve"].dropna(subset=["stability"])
    TIE = 0.01  # two settings this close are not meaningfully different
    mrows = []
    for (ds_, m), g in curve_all.groupby(["ds", "method"]):
        tied = g[g.stability >= g.stability.max() - TIE]
        rho = spearmanr(g.stability, g.ari).correlation if len(g.dropna(subset=["ari"])) > 4 else np.nan
        mrows.append({"ds": ds_, "method": m, "n_settings": len(g),
                      "share_tied": round(len(tied) / len(g), 3),
                      "k_tied_min": int(tied.n_clusters.min()), "k_tied_max": int(tied.n_clusters.max()),
                      "k_tied_span": int(tied.n_clusters.max() - tied.n_clusters.min()),
                      "spread": round(float(g.stability.max() - g.stability.min()), 3),
                      "spearman_stability_ari": round(float(rho), 3) if rho == rho else None})
    mech_all = pd.DataFrame(mrows)
    chosen = d[d.k_setting == "stability_finest"].groupby(["ds", "method"]).n_clusters.median().rename("chosen_k")
    mech_all = mech_all.merge(chosen.reset_index(), on=["ds", "method"], how="left")
    mech_all["true_k"] = mech_all.ds.map(truth_k)
    mech_all["excess_clusters"] = mech_all.chosen_k - mech_all.true_k
    # Aggregates are over the nine small-data datasets, the set the text quotes.
    mech = mech_all[mech_all.ds.isin(SMALL)].reset_index(drop=True)
    nums["mechanism"] = {
        "tie_window": TIE,
        "mean_share_tied": mech.groupby("method").share_tied.mean().round(3).to_dict(),
        "mean_k_span_among_tied": mech.groupby("method").k_tied_span.mean().round(1).to_dict(),
        "mean_spread": mech.groupby("method").spread.mean().round(3).to_dict(),
        "mean_spearman_stability_ari": mech.groupby("method").spearman_stability_ari.mean().round(3).to_dict(),
        "spearman_spread_vs_excess_clusters": {
            m: round(float(spearmanr(g.spread, g.excess_clusters).correlation), 3)
            for m, g in mech.dropna(subset=["excess_clusters"]).groupby("method") if len(g) > 3},
        "n_datasets": int(mech.ds.nunique()),
        "table": mech.to_dict(orient="records")}
    mech_all.to_csv(RESULTS / "a3_mechanism.csv", index=False)
    print("mechanism: share of settings tied with the best " + str(nums["mechanism"]["mean_share_tied"]))
    print("mechanism: mean span of cluster counts among tied settings " + str(nums["mechanism"]["mean_k_span_among_tied"]))

    # ---------- does MICA's shipped default suit this benchmark in particular? ----------
    # The k=4 analysis is about MICA's small-data default, so PBMC14K's graph-mode default stays out.
    dflt = d[(d.k_setting == "default") & (d.method == "mica") & d.ds.isin(SMALL)]
    if len(dflt):
        g = dflt.groupby("ds").agg(ari=("ari", "median")).reset_index()
        g["true_k"] = g.ds.map(truth_k)
        near = g[g.true_k.between(3, 5)]
        far = g[~g.true_k.between(3, 5)]
        nums["mica_default"] = {
            "shipped_k": 4, "source": "MICA/mica_mds.py, -nck default",
            "median_true_k_of_benchmark": int(np.median(list(truth_k.values()))),
            "datasets_with_true_k_4": [k for k, v in truth_k.items() if v == 4],
            "mean_ari_true_k_3_to_5": round(float(near.ari.mean()), 3),
            "mean_ari_true_k_outside": round(float(far.ari.mean()), 3),
            "n_near": int(len(near)), "n_far": int(len(far))}
        print(f"MICA default k=4: ARI {nums['mica_default']['mean_ari_true_k_3_to_5']} where the true "
              f"count is 3 to 5, {nums['mica_default']['mean_ari_true_k_outside']} outside")

    # ---------- how sensitive is the choice to the threshold itself? ----------
    # Relaxing the threshold admits more settings, so the finest choice can only move up.
    # Only this direction is meaningful: raising it can empty the passing set entirely,
    # which says nothing about the rule. The shifts are the size of the difference the
    # bootstrap flavour could plausibly make to the threshold.
    curve0 = d[(d.k_setting == "stability_curve") & d.ci_low.notna()]
    def finest_at(g, bump):
        u = g[g.n_clusters >= 2]
        passing = u[u.stability >= u.ci_low.max() + bump]
        if not len(passing):
            return None
        return int(passing.sort_values(["n_clusters", "stability"], ascending=[False, False]).iloc[0].n_clusters)
    # The check is applied to each run's own curve, setting by setting, exactly as the rule ran,
    # so its starting pick is the chosen count in Table I. An earlier version pooled the seeds
    # into a median curve per cluster count, which is not the rule and disagreed with Table I.
    srows = []
    for f in sorted(glob(str(RESULTS / "rule_stability_paper_*" / "metrics.csv"))):
        run = Path(f).parent.name.replace("rule_stability_paper_", "")
        ds_ = run.split("_graph")[0].split("_seed")[0]
        g_all = pd.read_csv(f)
        g_all = g_all[g_all.k_setting == "stability_curve"]
        for m, g in g_all.groupby("method"):
            base = finest_at(g, 0.0)
            shifted = {f"minus_{abs(b):.2f}": finest_at(g, b) for b in (-0.01, -0.02)}
            srows.append({"ds": ds_, "run": run, "method": m, "chosen_k": base, **shifted,
                          "changed": any(v is not None and v != base for v in shifted.values())})
    sens = pd.DataFrame(srows)
    if len(sens):
        # Headline over the nine small-data datasets, to match every other aggregate.
        sm = sens[sens.ds.isin(SMALL)]
        pb = sens[sens.ds == "pbmc14k"]
        nums["threshold_sensitivity"] = {
            "shifts": [-0.01, -0.02], "unit": "runs, each with its own curve",
            "changed_of_total": f"{int(sm.changed.sum())} of {len(sm)}",
            "by_method": {k: f"{int(v.sum())} of {len(v)}" for k, v in sm.groupby("method").changed},
            "pbmc14k": {k: f"{int(v.sum())} of {len(v)}" for k, v in pb.groupby("method").changed},
            "pbmc14k_shifted_counts": {k: v["minus_0.02"].tolist() for k, v in pb.groupby("method")},
            "table": sens.to_dict(orient="records")}
        sens.to_csv(RESULTS / "a3_threshold_sensitivity.csv", index=False)
        print(f"threshold sensitivity: the chosen cluster count changes in "
              f"{nums['threshold_sensitivity']['changed_of_total']} cases when the threshold is "
              f"relaxed by at most 0.02")

    nums["true_k"] = truth_k
    (RESULTS / "a3_report_numbers.json").write_text(json.dumps(nums, indent=1, default=str))

    # ---------- mechanism figure ----------
    curve = d[d.k_setting == "stability_curve"].copy()
    if len(curve):
        show = [x for x in ["buettner", "pollen", "zeisel"] if x in set(curve.ds)][:3]
        if show:
            fig, axes = plt.subplots(1, len(show), figsize=(7.0, 1.85), sharey=True)
            axes = np.atleast_1d(axes)
            plt.rcParams.update({"font.size": 7, "font.family": "serif", "pdf.fonttype": 42})
            for ax, x in zip(axes, show):
                for m, col, lab in [("mica", BLUE, "MICA"), ("leiden_pca", ORANGE, "Leiden")]:
                    s = curve[(curve.ds == x) & (curve.method == m)].groupby("n_clusters").agg(
                        stab=("stability", "median"), ari=("ari", "median")).reset_index().sort_values("n_clusters")
                    if not len(s):
                        continue
                    ax.plot(s.n_clusters, s.stab, color=col, lw=1.5, label=f"{lab} stability")
                    ax.plot(s.n_clusters, s.ari, color=col, lw=1.2, ls=(0, (2, 1.5)), label=f"{lab} ARI")
                ax.axvline(truth_k[x], color=MUTED, lw=0.8, ls=":")
                ax.text(truth_k[x], 1.04, "true", fontsize=6, color=MUTED, ha="center")
                ax.set_title(NAMES[x], fontsize=7.5, color=INK)
                ax.set_xlabel("clusters")
                ax.set_ylim(-0.05, 1.1)
                ax.grid(axis="y", color="#e6e5e1", lw=0.5)
                ax.set_axisbelow(True)
                for sp in ("top", "right"):
                    ax.spines[sp].set_visible(False)
            axes[0].set_ylabel("score")
            axes[0].legend(fontsize=5.5, frameon=False, loc="lower left", ncol=2, columnspacing=0.8,
                           handlelength=1.5, labelspacing=0.2)
            fig.tight_layout(pad=0.3)
            # No creation date, so regenerating an unchanged figure leaves the file
            # byte-identical and the repository clean.
            fig.savefig(FIGS / "a3_stability_curve.pdf", metadata={"CreationDate": None})
            fig.savefig(FIGS / "a3_stability_curve.png", dpi=200)
            print("wrote the stability curve figure")

    for k, v in nums.items():
        if isinstance(v, dict) and "mica_mean" in v:
            print(f"{k:20s} MICA {v['mica_mean']:.3f}  Leiden {v['leiden_mean']:.3f}  "
                  f"MICA wins {v['mica_wins']} of {v['n_datasets']}")
    print(f"wrote tables to {TABLES}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
