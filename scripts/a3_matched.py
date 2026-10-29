"""Matched-grid results (U-40), the transfer comparison (U-34) and the summary table (U-37).

Usage: python scripts/a3_matched.py
Reads results/rule_matched_* and the original runs through a3_report.py's loader. Nothing
here reruns an experiment. Writes latex/tables/a3_matched.tex, a3_transfer.tex and
a3_summary.tex, results/a3_transfer.csv and results/a3_matched_numbers.json, which holds every
number the text quotes from these runs.

The design was fixed before the matched runs (PROGRESS.md, Oct 28, "Matched grids, declared
before any matched run"): every method is offered the cluster counts from 2 to 20 that all
methods reach on the full data, with one setting per count.
"""
import json
import sys
from glob import glob
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import a3_report as rep
from scminer_clust.paths import REPO, RESULTS

TABLES = REPO / "latex" / "tables"
TIE_WINDOW = 0.01  # as in a3_report.py
READINGS = [("stability_finest", "F"), ("stability_most", "M")]


def load_matched():
    """Picks, stability curves and full-data sweeps of the matched stability runs."""
    rows = []
    for f in sorted(glob(str(RESULTS / "rule_matched_stability_paper_*" / "metrics.csv"))):
        d = pd.read_csv(f)
        d["run"] = Path(f).parent.name.replace("rule_matched_stability_paper_", "")
        rows.append(d)
    if not rows:
        return pd.DataFrame()
    d = pd.concat(rows, ignore_index=True)
    d["ds"] = d.dataset.map(rep.short)
    d["seed"] = d.seed.fillna(0).astype(int)
    return d


def seed_means(d, rule, method, datasets):
    """Mean ARI over the datasets at each seed, for the range of the mean across seeds."""
    s = d[(d.k_setting == rule) & (d.method == method) & d.ds.isin(datasets)]
    per = s.groupby("seed").agg(mean=("ari", "mean"), n=("ari", "size"))
    return per[per.n == len(datasets)]["mean"]


def main():
    nums = {}
    facts = pd.read_csv(REPO / "docs" / "dataset_facts.csv").set_index("dataset")
    m = load_matched()
    if not len(m):
        print("no matched runs yet")
        return 1
    present = [x for x in rep.ORDER if x in set(m.ds)]
    small = [x for x in rep.SMALL if x in present]
    truth_k = {x: int(facts.loc[f"paper_{x}", "label_classes"]) for x in present}
    picks = m[m.k_setting.isin([r for r, _ in READINGS])]
    curves = m[m.k_setting == "stability_curve"]

    # ---------- the matched counts each run offered ----------
    offered = (curves[curves.method == "mica"].groupby(["ds", "seed"]).n_clusters
               .agg(lambda s: sorted(set(s.astype(int)))))
    nums["offered_counts"] = {f"{x}/seed{s}": v for (x, s), v in offered.items()}
    n_off = offered.map(len)

    # ---------- the rule table on matched grids ----------
    tabs = {r: rep.agg(m, r) for r, _ in READINGS}
    body = []
    for x in present:
        cells = []
        for r, _ in READINGS:
            for meth in ("mica", "leiden_pca"):
                sel = tabs[r][(tabs[r].ds == x) & (tabs[r].method == meth)]
                if not len(sel):
                    cells.append("n/a")
                    continue
                v = sel.iloc[0]
                k = f"{int(v.k)}" if v.k_min == v.k_max else f"{int(v.k_min)}\\textrm{{-}}{int(v.k_max)}"
                cells.append(f"{v.ari:.3f}\\,({k})")
        no = n_off.loc[x] if x in n_off.index.get_level_values(0) else pd.Series(dtype=int)
        offer = f"{int(no.min())}" if no.min() == no.max() else f"{int(no.min())}\\textrm{{-}}{int(no.max())}"
        body.append(f"{rep.NAMES[x]} ({truth_k[x]}) & {offer} & " + " & ".join(cells) + " \\\\")
    means = {}
    for r, _ in READINGS:
        for meth in ("mica", "leiden_pca"):
            sel = tabs[r][tabs[r].ds.isin(small) & (tabs[r].method == meth)]
            means[(r, meth)] = float(sel.ari.mean()) if len(sel) == len(small) else None
    mrow = " & ".join("n/a" if means[(r, meth)] is None else f"\\textbf{{{means[(r, meth)]:.3f}}}"
                      for r, _ in READINGS for meth in ("mica", "leiden_pca"))
    body.append("\\midrule\n\\textit{Mean, nine small} & & " + mrow + " \\\\")
    (TABLES / "a3_matched.tex").write_text(
        "\\begin{table}[t]\n\\centering\n\\caption{The stability rule on matched candidates. In each run both "
        "methods are offered the same cluster counts, those from 2 to 20 that both reach on the full data, with one "
        "setting per count. The second column is how many of the 19 counts that is, over the three seeds. ARI with "
        "the chosen count in parentheses, medians over three seeds, F the finest-stable reading and M the plain "
        "maximum. The design was fixed before these runs and after the post hoc check in "
        "Table~\\ref{tab:common}.}\n\\label{tab:matched}\n\\footnotesize\n\\setlength{\\tabcolsep}{3pt}\n"
        "\\resizebox{\\columnwidth}{!}{%\n\\begin{tabular}{lccccc}\n\\toprule\n & & \\multicolumn{2}{c}{F} & "
        "\\multicolumn{2}{c}{M} \\\\\n\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\n"
        "Dataset (types) & Counts & MICA & Leiden & MICA & Leiden \\\\\n\\midrule\n"
        + "\n".join(body) + "\n\\bottomrule\n\\end{tabular}}\n\\end{table}\n")

    # ---------- numbers for the text ----------
    out = {"means_nine_small": {f"{r}/{meth}": None if v is None else round(v, 4) for (r, meth), v in means.items()}}
    for r, _ in READINGS:
        per = {meth: seed_means(m, r, meth, small) for meth in ("mica", "leiden_pca")}
        if all(len(v) for v in per.values()):
            out[f"{r}_per_seed"] = {meth: {int(s): round(float(v), 4) for s, v in per[meth].items()} for meth in per}
            lead = (per["mica"] - per["leiden_pca"]).dropna()
            out[f"{r}_mica_leads_at"] = f"{int((lead > 0).sum())} of {len(lead)} seeds"
    med = {r: tabs[r][tabs[r].ds.isin(small)].pivot(index="ds", columns="method", values="ari") for r, _ in READINGS}
    gap = (med["stability_finest"] - med["stability_most"]).abs()
    out["mean_abs_reading_change"] = {meth: round(float(gap[meth].mean()), 4) for meth in gap.columns}
    out["sign_mean_reading_change"] = {meth: round(float((med["stability_finest"] - med["stability_most"])[meth].mean()), 4)
                                       for meth in gap.columns}
    out["datasets_where_mica_leads"] = {r: f"{int((med[r]['mica'] > med[r]['leiden_pca']).sum())} of {len(med[r])}"
                                        for r, _ in READINGS}
    # is the finest-stable pick at an end of the counts the run offered?
    edge = []
    for _, p in picks[(picks.k_setting == "stability_finest") & picks.ds.isin(small)].iterrows():
        offs = offered.loc[(p.ds, p.seed)]
        edge.append({"method": p.method, "edge": int(p.n_clusters) in (min(offs), max(offs)),
                     "top": int(p.n_clusters) == max(offs)})
    E = pd.DataFrame(edge)
    out["finest_pick_at_edge"] = {meth: f"{int(g.edge.sum())} of {len(g)} runs, {int(g.top.sum())} at the top"
                                  for meth, g in E.groupby("method")}
    # mechanism on matched candidates: how much of the curve ties with its best
    mech = []
    for (x, meth, s), g in curves[curves.ds.isin(small)].groupby(["ds", "method", "seed"]):
        t = g[g.stability >= g.stability.max() - TIE_WINDOW]
        mech.append({"ds": x, "method": meth, "seed": s, "share_tied": len(t) / len(g),
                     "span": int(t.n_clusters.max() - t.n_clusters.min()),
                     "share_at_one": float((g.stability >= 0.9999).mean())})
    K = pd.DataFrame(mech)
    out["mechanism"] = {meth: {"share_tied": round(float(g.share_tied.mean()), 3),
                               "tied_span": round(float(g.span.mean()), 2),
                               "share_at_one": round(float(g.share_at_one.mean()), 3)}
                        for meth, g in K.groupby("method")}
    out["offered_counts_per_run"] = {"min": int(n_off[n_off.index.get_level_values(0).isin(small)].min()),
                                     "max": int(n_off[n_off.index.get_level_values(0).isin(small)].max()),
                                     "mean": round(float(n_off[n_off.index.get_level_values(0).isin(small)].mean()), 1)}
    if "pbmc14k" in present:
        pb = picks[picks.ds == "pbmc14k"][["method", "k_setting", "seed", "n_clusters", "ari"]]
        out["pbmc14k"] = {f"{r.method}/{r.k_setting}/seed{r.seed}": {"k": int(r.n_clusters), "ari": round(float(r.ari), 4)}
                          for _, r in pb.iterrows()}
        out["pbmc14k_offered"] = {f"seed{s}": v for (x, s), v in offered.items() if x == "pbmc14k"}
    # MICA's small-data k grid was already complete, so on runs where every count is offered
    # its picks must equal the original runs. A check that the matched runner changes nothing else.
    orig = rep.load_stability()
    chk = []
    for (x, s), offs in offered.items():
        if x not in small or len(offs) != 19:
            continue
        for r, _ in READINGS:
            a = picks[(picks.ds == x) & (picks.seed == s) & (picks.method == "mica") & (picks.k_setting == r)]
            b = orig[(orig.ds == x) & (orig.seed == s) & (orig.method == "mica") & (orig.k_setting == r)]
            if len(a) and len(b):
                chk.append(int(a.n_clusters.iloc[0]) == int(b.n_clusters.iloc[0]))
    out["mica_unchanged_where_all_counts_offered"] = f"{sum(chk)} of {len(chk)} picks identical"
    nums["matched"] = out
    print("matched means: " + ", ".join(f"{k} {v}" for k, v in out["means_nine_small"].items()))
    print("matched abs reading change: " + str(out["mean_abs_reading_change"]) + ", edges " + str(out["finest_pick_at_edge"]))
    print("matched mechanism: " + str(out["mechanism"]) + ", " + out["mica_unchanged_where_all_counts_offered"])

    # ---------- the mode test on matched candidates ----------
    mode_rows = []
    for run in sorted(glob(str(RESULTS / "rule_matched_mode_size_pbmc14k_n*_seed*"))):
        seed = int(run.split("seed")[-1])
        mm = pd.read_csv(Path(run) / "metrics.csv")
        ts = pd.read_csv(Path(run) / "tie_summary.csv")
        for _, t in ts.iterrows():
            p = mm[(mm.method == t.method) & (mm.mica_mode == t["mode"]) & mm.k_setting.isin([r for r, _ in READINGS])]
            w = p.set_index("k_setting")
            mode_rows.append({"seed": seed, "arm": f"{t.method}/{t['mode']}", "n_settings": int(t.n_settings),
                              "share_tied": float(t.share_tied), "span": int(t.k_tied_span),
                              "finest_k": int(w.loc["stability_finest", "n_clusters"]),
                              "most_k": int(w.loc["stability_most", "n_clusters"]),
                              "gap": float(w.loc["stability_finest", "ari"] - w.loc["stability_most", "ari"])})
    if mode_rows:
        Mo = pd.DataFrame(mode_rows)
        nums["matched_mode_size"] = {
            a: {"counts_offered": g.n_settings.tolist(), "spans": g.span.tolist(),
                "share_tied_median": round(float(g.share_tied.median()), 3),
                "readings_disagree": f"{int((g.gap.abs() > 0.001).sum())} of {len(g)}",
                "gaps": [round(x, 3) for x in g.gap], "finest_k": g.finest_k.tolist(), "most_k": g.most_k.tolist()}
            for a, g in Mo.groupby("arm")}
        print("matched mode test: " + ", ".join(f"{a} spans {v['spans']} disagree {v['readings_disagree']}"
                                               for a, v in nums["matched_mode_size"].items()))

    # ---------- transfer: the setting that did best on the other datasets' labels ----------
    sw = m[m.k_setting == "full_sweep"].copy()
    trows = []
    for meth in ("mica", "leiden_pca"):
        for s in sorted(sw.seed.unique()):
            g = sw[(sw.method == meth) & (sw.seed == s) & sw.ds.isin(small)]
            P = g.pivot_table(index="setting", columns="ds", values="ari")
            K_ = g.pivot_table(index="setting", columns="ds", values="n_clusters")
            if P.shape[1] != len(small):
                continue
            P = P.dropna()
            for x in small:
                train = P.drop(columns=x).mean(axis=1)
                best = train[train == train.max()].index.min()  # ties go to the smaller setting
                trows.append({"ds": x, "method": meth, "seed": int(s), "setting": float(best),
                              "n_clusters": int(K_.loc[best, x]), "ari": float(P.loc[best, x])})
            if meth == "leiden_pca":
                pb = sw[(sw.method == meth) & (sw.seed == s) & (sw.ds == "pbmc14k")].set_index("setting")
                if len(pb):
                    best = P.mean(axis=1)
                    best = best[best == best.max()].index.min()
                    if best in pb.index:
                        trows.append({"ds": "pbmc14k", "method": meth, "seed": int(s), "setting": float(best),
                                      "n_clusters": int(pb.loc[best, "n_clusters"]), "ari": float(pb.loc[best, "ari"])})
    if trows:
        T = pd.DataFrame(trows)
        T.to_csv(RESULTS / "a3_transfer.csv", index=False)
        tm = T.groupby(["ds", "method"]).agg(ari=("ari", "median"), ari_min=("ari", "min"), ari_max=("ari", "max"),
                                             k=("n_clusters", "median"), k_min=("n_clusters", "min"),
                                             k_max=("n_clusters", "max")).reset_index()
        tb = []
        for x in [d_ for d_ in rep.ORDER if d_ in set(tm.ds)]:
            cells = []
            for meth in ("mica", "leiden_pca"):
                sel = tm[(tm.ds == x) & (tm.method == meth)]
                if not len(sel):
                    cells.append("n/a")
                    continue
                v = sel.iloc[0]
                k = f"{int(v.k)}" if v.k_min == v.k_max else f"{int(v.k_min)}\\textrm{{-}}{int(v.k_max)}"
                cells.append(f"{v.ari:.3f}\\,({k})")
            tb.append(f"{rep.NAMES[x]} ({truth_k[x]}) & " + " & ".join(cells) + " \\\\")
        tmean = tm[tm.ds.isin(small)].groupby("method").ari.mean()
        tb.append("\\midrule\n\\textit{Mean, nine small} & " +
                  " & ".join(f"\\textbf{{{tmean[meth]:.3f}}}" for meth in ("mica", "leiden_pca")) + " \\\\")
        (TABLES / "a3_transfer.tex").write_text(
            "\\begin{table}[t]\n\\centering\n\\caption{The transfer comparison, which is not label-free. Each small "
            "dataset takes the setting with the best mean ARI over the other eight at the same seed (MICA: $k$, "
            "Leiden: a resolution), and PBMC14K takes the baseline's resolution best on the nine. MICA's graph mode "
            "has no setting that carries over from $k$. ARI with the resulting cluster count in parentheses, medians "
            "over three seeds.}\n\\label{tab:transfer}\n\\footnotesize\n\\setlength{\\tabcolsep}{4pt}\n"
            "\\begin{tabular}{lcc}\n\\toprule\nDataset (types) & MICA & Leiden \\\\\n\\midrule\n"
            + "\n".join(tb) + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n")
        mk = T[(T.method == "mica")]
        nums["transfer"] = {
            "means_nine_small": {meth: round(float(v), 4) for meth, v in tmean.items()},
            "per_seed_means": {meth: {int(s): round(float(v), 4) for s, v in
                                      T[(T.method == meth) & T.ds.isin(small)].groupby("seed").ari.mean().items()}
                               for meth in ("mica", "leiden_pca")},
            "mica_k_chosen": sorted(set(mk.setting.astype(int))),
            "mica_k_equals_default_4": f"{int((mk.setting == 4).sum())} of {len(mk)} held-out choices",
            "leiden_resolutions_chosen": sorted(set(T[T.method == "leiden_pca"].setting.round(4))),
            "pbmc14k_leiden": {f"seed{int(r.seed)}": {"k": int(r.n_clusters), "ari": round(r.ari, 4)}
                               for _, r in T[T.ds == "pbmc14k"].iterrows()}}
        print("transfer: " + str(nums["transfer"]["means_nine_small"]) + ", MICA k " + str(nums["transfer"]["mica_k_chosen"]))

    # ---------- the summary table: every rule, both methods, no win count ----------
    o = rep.load_stability()
    rows = [("Default", o, "default"), ("Silhouette, shared space", o, "silhouette_shared"),
            ("Stability F, original grids", o, "stability_finest"), ("Stability M, original grids", o, "stability_most"),
            ("Stability F, matched", m, "stability_finest"), ("Stability M, matched", m, "stability_most")]
    lines, summ = [], {}

    def rng(d, r, meth):
        per = seed_means(d, r, meth, small)
        return f" [{per.min():.2f}, {per.max():.2f}]" if len(per) > 1 and per.max() - per.min() >= 0.005 else ""

    def secs(d, r, meth):
        s = d[(d.k_setting == r) & (d.method == meth) & d.ds.isin(small)]
        return f"{s.seconds.median():.0f}" if "seconds" in s and s.seconds.notna().any() else ("none" if r == "default" else "n/a")

    def pbmc(d, r, meth):
        s = d[(d.k_setting == r) & (d.method == meth) & (d.ds == "pbmc14k")]
        return f"{s.ari.median():.3f}" if len(s) else "n/a"

    for label, d, r in rows:
        cells, rec = [], {}
        for meth in ("mica", "leiden_pca"):
            a = rep.agg(d, r)
            sel = a[a.ds.isin(small) & (a.method == meth)] if len(a) else a
            mu = float(sel.ari.mean()) if len(sel) == len(small) else None
            cells.append("n/a" if mu is None else f"{mu:.3f}{rng(d, r, meth)}")
            rec[meth] = None if mu is None else round(mu, 4)
        cells += [pbmc(d, r, "mica"), pbmc(d, r, "leiden_pca"), secs(d, r, "mica"), secs(d, r, "leiden_pca")]
        lines.append(f"{label} & " + " & ".join(cells) + " \\\\")
        summ[label] = rec
    if trows:
        T_ = pd.DataFrame(trows)
        tp = T_[T_.ds == "pbmc14k"]
        lines.append("\\midrule\n\\multicolumn{7}{l}{\\textit{Uses the other datasets' labels}} \\\\")
        tr = {meth: T_[(T_.method == meth) & T_.ds.isin(small)].groupby("ds").ari.median().mean() for meth in ("mica", "leiden_pca")}
        lines.append(f"Transfer & {tr['mica']:.3f} & {tr['leiden_pca']:.3f} & n/a & "
                     f"{tp.ari.median():.3f} & none & none \\\\" if len(tp) else
                     f"Transfer & {tr['mica']:.3f} & {tr['leiden_pca']:.3f} & n/a & n/a & none & none \\\\")
        summ["Transfer"] = {k: round(float(v), 4) for k, v in tr.items()}
    (TABLES / "a3_summary.tex").write_text(
        "\\begin{table*}[t]\n\\centering\n\\caption{Every rule in one place. Mean ARI over the nine small datasets of "
        "the per-dataset medians, with the range of that mean over the three seeds where a rule was run at three. "
        "PBMC14K is reported apart, as medians. Selection time is the median per small dataset in seconds, for the "
        "rules that select by rerunning. The default rule and the transfer select nothing at run time, and the "
        "silhouette rule was not timed. "
        "There is no win count: the point is how far each method's mean moves with the rule.}\n"
        "\\label{tab:summary}\n\\footnotesize\n\\setlength{\\tabcolsep}{4pt}\n"
        "\\begin{tabular}{lcccccc}\n\\toprule\n & \\multicolumn{2}{c}{Nine small datasets} & "
        "\\multicolumn{2}{c}{PBMC14K} & \\multicolumn{2}{c}{Selection (s)} \\\\\n"
        "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\n"
        "Rule & MICA & Leiden & MICA & Leiden & MICA & Leiden \\\\\n\\midrule\n"
        + "\n".join(lines) + "\n\\bottomrule\n\\end{tabular}\n\\end{table*}\n")
    nums["summary"] = summ
    (RESULTS / "a3_matched_numbers.json").write_text(json.dumps(nums, indent=1, default=str))
    print(f"wrote {RESULTS / 'a3_matched_numbers.json'} and the tables")
    return 0


if __name__ == "__main__":
    sys.exit(main())
