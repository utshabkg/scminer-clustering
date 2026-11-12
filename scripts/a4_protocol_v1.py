"""Protocol v1's criteria and exploratory comparisons, from the saved picks and curves (U-44).

Usage: python scripts/a4_protocol_v1.py [--seeds 3-9] [--pbmc-seeds 0-4]
Reads results/v1/{picks,curves,picks_chooser,chooser_r}/ and results/v1/regression.csv, and
writes results/a4_protocol_v1_numbers.json. The criteria and their thresholds are fixed in
docs/protocol_matched_v1.md. Summary of a design: the per-dataset median ARI over the seeds,
then the mean over the nine small datasets, as in Assignment 3 (scripts/a3_report.py, agg).

Scorers, as named in the protocol: S0 ours at B=25 (the original protocol), S1 ours at B=100,
S2 ours at B=100 on all cells (PBMC14K), S3-V, S3-G and S3-G2 chooseR's own code at B=100.
"""
import json
import sys
from glob import glob

import numpy as np
import pandas as pd

from scminer_clust import stability
from scminer_clust.paths import RESULTS

V1 = RESULTS / "v1"
SMALL = ["yan", "goolam", "buettner", "pollen", "chung", "usoskin", "kolodziejczyk", "klein", "zeisel"]
FIVE = ["yan", "goolam", "buettner", "pollen", "chung"]  # at most 317 cells
METHODS = ["mica", "leiden_pca", "kmeans_pca"]
F, M = "stability_finest", "stability_most"


def short(ds):
    return ds.replace("paper_", "")


def scorer_name(r):
    if r.scorer.startswith("chooseR"):
        return "S3-" + r.scorer.split("_")[1]
    if r.scorer == "ours_all":
        return "S2"
    return "S0" if int(r.b) == 25 else "S1"


def load(sub):
    files = sorted(glob(str(V1 / sub / "*.csv")))
    files = [f for f in files if "crosscheck" not in f]
    if not files:
        return pd.DataFrame()
    d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    d["ds"] = d.dataset.map(short)
    d["S"] = [scorer_name(r) for r in d.itertuples()]
    return d


def summary(p, S, design, reading, method, seeds, datasets=SMALL):
    s = p[(p.S == S) & (p.design == design) & (p.reading == reading) & (p.method == method)
          & p.seed.isin(seeds) & p.ds.isin(datasets)]
    if not len(s):
        return None
    per_ds = s.groupby("ds").ari.median()
    per_seed = s.groupby("seed").agg(mean=("ari", "mean"), n=("ari", "size"))
    per_seed = per_seed[per_seed.n == len(datasets)]["mean"]
    return {"mean": round(float(per_ds.mean()), 4), "n_datasets": int(per_ds.size),
            "n_runs": int(len(s)), "seeds": sorted(int(x) for x in s.seed.unique()),
            "per_seed": {int(k): round(float(v), 4) for k, v in per_seed.items()},
            "per_dataset": {k: round(float(v), 4) for k, v in per_ds.items()},
            "k_median": {k: float(v) for k, v in s.groupby("ds").n_clusters.median().items()}}


def lead(p, S, design, reading, a, b, seeds):
    x, y = summary(p, S, design, reading, a, seeds), summary(p, S, design, reading, b, seeds)
    if not x or not y:
        return None
    common = sorted(set(x["per_seed"]) & set(y["per_seed"]))
    leads = sum(x["per_seed"][s] > y["per_seed"][s] for s in common)
    return {"a": a, "b": b, "mean_a": x["mean"], "mean_b": y["mean"], "gap": round(x["mean"] - y["mean"], 4),
            "a_leads_seeds": f"{leads} of {len(common)}", "leads": leads, "n_seeds": len(common)}


def verdict_c1(r, need):
    if r is None or r["n_seeds"] < need:
        return "not yet computable"
    if r["gap"] > 0 and r["leads"] >= 5:
        return "holds"
    return "weak" if r["gap"] > 0 else "fails"


def reading_change(p, S, design, method, seeds):
    s = p[(p.S == S) & (p.design == design) & (p.method == method) & p.seed.isin(seeds) & p.ds.isin(SMALL)]
    med = s.pivot_table(index="ds", columns="reading", values="ari", aggfunc="median")
    if F not in med or M not in med:
        return None
    return round(float((med[F] - med[M]).abs().mean()), 4)


def flat_top(c, seeds, design="mid", b=25):
    """Share of settings within 0.01 of the best, pooling each dataset's curves over the seeds,
    as in scripts/a3_matched.py (mechanism)."""
    s = c[(c.scorer == "ours") & (c.b == b) & c[f"in_{design}"] & c.seed.isin(seeds) & c.ds.isin(SMALL)]
    rows = []
    for (ds, m), g in s.groupby(["ds", "method"]):
        t = g[g.stability >= g.stability.max() - 0.01]
        rows.append({"ds": ds, "method": m, "share_tied": len(t) / len(g),
                     "span": int(t.n_clusters.max() - t.n_clusters.min()),
                     "share_at_one": float((g.stability >= 0.9999).mean())})
    K = pd.DataFrame(rows)
    if not len(K):
        return {}
    K["size"] = np.where(K.ds.isin(FIVE), "five_smallest", "four_larger")
    out = {f"{sz}/{m}": {"share_tied": round(float(g.share_tied.mean()), 3), "span": round(float(g.span.mean()), 2),
                         "share_at_one": round(float(g.share_at_one.mean()), 3)} for (sz, m), g in K.groupby(["size", "method"])}
    out["all"] = {m: {"share_tied": round(float(g.share_tied.mean()), 3), "span": round(float(g.span.mean()), 2)}
                  for m, g in K.groupby("method")}
    return out


def repick(curve_rows, settings):
    """stability.pick on a subset of one run's curve, in the given order."""
    sc = {}
    for s in settings:
        r = curve_rows.loc[s]
        sc[s] = {"stability": r.stability, "ci_low": r.ci_low, "ci_high": r.ci_high, "n_clusters": int(r.n_clusters)}
    return stability.pick(sc)


def decomposition(c, seeds, b=25):
    """E3: Leiden's mean under four candidate sets that change one thing at a time."""
    s = c[(c.scorer == "ours") & (c.b == b) & (c.method == "leiden_pca") & c.seed.isin(seeds) & c.ds.isin(SMALL)]
    out = []
    for (ds, seed), g in s.groupby(["ds", "seed"]):
        g = g.drop_duplicates("setting").set_index("setting").sort_index()
        orig = g[g.in_original]
        old_counts = sorted(orig.n_clusters.unique())
        by = {}
        for st, r in orig.iterrows():
            by.setdefault(int(r.n_clusters), []).append(st)
        orig_one = [v[(len(v) - 1) // 2] for k, v in sorted(by.items())]
        mid = g[g.in_mid].sort_values("n_clusters")
        mid_old = mid[mid.n_clusters.isin(old_counts)]
        sets = {"original": list(orig.index), "original_one_per_count": orig_one,
                "mid_old_counts": list(mid_old.index), "mid": list(mid.index)}
        for name, settings in sets.items():
            if not settings:
                continue
            for reading, info in repick(g, settings).items():
                out.append({"ds": ds, "seed": seed, "set": name, "reading": reading,
                            "ari": float(g.loc[info["chosen_setting"], "ari"]), "k": info["n_clusters"]})
    D = pd.DataFrame(out)
    if not len(D):
        return {}
    res = {}
    for (name, reading), g in D.groupby(["set", "reading"]):
        res[f"{name}/{reading}"] = round(float(g.groupby("ds").ari.median().mean()), 4)
    return res


def chooser_agreement(p, seeds, design="mid"):
    """E4: S1 against S3-G, finest reading: same chosen count, and the ARI change."""
    a = p[(p.S == "S1") & (p.reading == F) & (p.design == design) & p.seed.isin(seeds)]
    b = p[(p.S == "S3-G") & (p.design == design) & p.seed.isin(seeds)]
    j = a.merge(b, on=["ds", "seed", "method"], suffixes=("_s1", "_s3"))
    if not len(j):
        return {}
    return {m: {"runs": int(len(g)), "same_count": int((g.n_clusters_s1 == g.n_clusters_s3).sum()),
                "mean_abs_ari_change": round(float((g.ari_s1 - g.ari_s3).abs().mean()), 4)}
            for m, g in j.groupby("method")}


def seeds_of(arg, default):
    if not arg:
        return default
    a, b = arg.split("-")
    return list(range(int(a), int(b) + 1))


def main(argv):
    opts = dict(zip(argv[::2], argv[1::2]))
    conf = seeds_of(opts.get("--seeds"), list(range(3, 10)))
    pseeds = seeds_of(opts.get("--pbmc-seeds"), list(range(0, 5)))
    allseeds = list(range(0, 10))
    p = pd.concat([load("picks"), load("picks_chooser")], ignore_index=True)
    c = load("curves")
    out = {"seeds_confirmatory": conf, "pbmc_seeds": pseeds,
           "runs_present": {f"{ds}": sorted(int(s) for s in g.seed.unique()) for ds, g in p.groupby("ds")}}

    # R0
    reg = V1 / "regression.csv"
    if reg.exists():
        r = pd.read_csv(reg)
        out["R0"] = r.fillna("").to_dict(orient="records")
    # C1 and C5
    out["C1"] = {}
    for reading in (F, M):
        lr = lead(p, "S0", "mid", reading, "mica", "leiden_pca", conf)
        out["C1"][reading] = {**(lr or {}), "verdict": verdict_c1(lr, len(conf))}
    lr = lead(p, "S0", "original", M, "leiden_pca", "mica", conf)
    out["C5"] = {**(lr or {}), "verdict": "not yet computable" if lr is None else ("holds" if lr["gap"] > 0 else "fails")}
    # C2
    rc = {m: reading_change(p, "S0", "mid", m, conf) for m in ("mica", "leiden_pca", "kmeans_pca")}
    out["C2"] = {"mean_abs_reading_change": rc,
                 "verdict": "not yet computable" if None in (rc["mica"], rc["leiden_pca"])
                 else ("holds" if rc["mica"] > rc["leiden_pca"] else "fails")}
    # C3
    ft = flat_top(c, conf)
    if ft:
        a, b_, cl = (ft.get("five_smallest/mica", {}).get("share_tied"), ft.get("four_larger/mica", {}).get("share_tied"),
                     ft.get("five_smallest/leiden_pca", {}).get("share_tied"))
        ok = None if None in (a, b_, cl) else (a > b_ and a > cl)
        out["C3"] = {**ft, "verdict": "not yet computable" if ok is None else ("holds" if ok else "fails")}
    # C4
    pb = p[(p.ds == "pbmc14k") & (p.S == "S0") & (p.design == "mid") & (p.reading == M) & p.seed.isin(pseeds)]
    if len(pb):
        w = pb.pivot_table(index="seed", columns="method", values="ari")
        if {"mica", "leiden_pca"} <= set(w.columns):
            w = w.dropna(subset=["mica", "leiden_pca"])
            ties = int(((w.mica - w.leiden_pca).abs() <= 0.02).sum())
            out["C4"] = {"per_seed": w.round(4).to_dict(orient="index"), "ties": f"{ties} of {len(w)}",
                         "verdict": "not yet computable" if len(w) < len(pseeds) else ("holds" if ties >= 3 else "fails")}
    # F1
    lr = lead(p, "S3-G", "mid", F, "mica", "leiden_pca", conf)
    out["F1"] = {**(lr or {}), "verdict": "not yet computable" if lr is None else ("holds" if lr["gap"] > 0 else "fails")}
    # FC
    cc = sorted(glob(str(V1 / "chooser_r" / "crosscheck_*.csv")))
    if cc:
        x = pd.concat([pd.read_csv(f) for f in cc])
        out["FC"] = {"matrices": int(len(x)), "identical": int(x.equal.sum()), "max_abs_diff": float(x.max_abs_diff.max()),
                     "verdict": "holds" if bool(x.equal.all()) else "fails"}
    # summaries for the tables: every scorer, design, reading and method, confirmatory and all seeds
    tab = {}
    for (S, design, reading, method), _ in p.groupby(["S", "design", "reading", "method"]):
        for label, seeds in (("conf", conf), ("all", allseeds)):
            s = summary(p, S, design, reading, method, seeds)
            if s:
                tab[f"{S}/{design}/{reading}/{method}/{label}"] = {k: s[k] for k in ("mean", "n_datasets", "n_runs", "per_seed")}
    out["summaries"] = tab
    out["E2_menu_vs_representation"] = {
        f"{d}/{r}": {"kmeans_minus_leiden": (lead(p, "S0", d, r, "kmeans_pca", "leiden_pca", conf) or {}).get("gap"),
                     "mica_minus_kmeans": (lead(p, "S0", d, r, "mica", "kmeans_pca", conf) or {}).get("gap")}
        for d in ("original", "mid") for r in (F, M)}
    out["E3_decomposition_leiden"] = decomposition(c, conf)
    out["E4_chooser_agreement"] = {"mid": chooser_agreement(p, allseeds, "mid"),
                                   "original": chooser_agreement(p, allseeds, "original")}
    st = p[p.scorer.str.startswith("chooseR")]
    if len(st):
        out["S3_failed_settings"] = {m: int(g.n_failed.max()) for m, g in st.groupby("method")}
        out["S3_no_choice_runs_V"] = {m: int(g.chosen_setting.isna().sum()) for m, g in st[st.S == "S3-V"].groupby("method")}
    (RESULTS / "a4_protocol_v1_numbers.json").write_text(json.dumps(out, indent=1, default=str))
    for k in ("C1", "C2", "C3", "C4", "C5", "F1", "FC"):
        if k in out:
            v = out[k]
            print(k, v.get("verdict") if isinstance(v, dict) and "verdict" in v else "",
                  {kk: vv for kk, vv in (v.items() if isinstance(v, dict) else []) if kk not in ("per_seed",)} if k != "C3" else "")
    print("wrote", RESULTS / "a4_protocol_v1_numbers.json")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
