"""R-7 independent recomputation of the paper's newer numbers.

Reads the committed run files directly and recomputes each value without importing the
generators (scripts/a4_protocol_v1.py, scripts/large_scale/make_tables.py). Each value is
compared with what the paper prints, taken from latex/numbers_*.tex.

Run from the repository root:
    python scripts/tools/r7_recompute.py            # every group
    python scripts/tools/r7_recompute.py allen      # one group

Prints one line per value: MATCH or DIFF, the paper's value and the recomputed one.
"""
import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
LS = ROOT / "results" / "large_scale"


def paper_values():
    vals = {}
    for name in ("numbers_a4.tex", "numbers_large.tex"):
        text = (ROOT / "latex" / name).read_text()
        for m in re.finditer(r"\\csname v@(.+?)\\endcsname\{(.*)\}\s*$", text, re.M):
            vals[m.group(1)] = m.group(2)
    return vals


PAPER = paper_values()
RESULTS = []


def check(key, value, note=""):
    printed = PAPER.get(key)
    ok = printed is not None and printed == value
    RESULTS.append((key, ok))
    status = "MATCH" if ok else "DIFF "
    print(f"  {status} {key:42s} paper {str(printed):>16s}   recomputed {value:>16s}  {note}")


def r3(x):
    return f"{x:.3f}"


def thousands(x):
    return f"{int(x):,}"


def allen():
    print("\nGroup A: Allen atlas (Section VI-H, Table large_scale)")
    dims = json.loads((LS / "dims.json").read_text())
    stages = [json.loads(line) for line in (LS / "stages.jsonl").read_text().splitlines() if line.strip()]
    manifest = next(s for s in stages if s["stage"] == "manifest")
    check("large.n_cells", thousands(dims["n_cells"]))
    check("large.n_genes", thousands(dims["n_genes_kept"]))
    check("large.n_genes_raw", thousands(dims["n_genes_raw"]))
    print(f"         (classes in the manifest: {manifest['n_classes']}, rows in the files: {manifest['n_rows_in_files']:,})")

    ps = pd.read_csv(LS / "partition_scores.csv")
    km = ps[ps.method == "PCA + k-means"].sort_values("seed")
    ld = ps[(ps.method == "PCA + Leiden") & (ps.setting == "resolution 0.05")].sort_values("seed")
    ld_def = ps[(ps.method == "PCA + Leiden") & (ps.setting == "resolution 1.00") & (ps.seed == 0)].iloc[0]
    check("large.kmeans.n_seeds", str(km.seed.nunique()))
    check("large.kmeans.ari", r3(km.ari.iloc[0]), "seed 0")
    check("large.kmeans.ami", r3(km.ami.iloc[0]), "seed 0")
    check("large.kmeans.ari_seeds", f"{r3(km.ari.min())} to {r3(km.ari.max())}")
    check("large.leiden.ari", r3(ld.ari.iloc[0]), "seed 0")
    check("large.leiden.ami", r3(ld.ami.iloc[0]), "seed 0")
    check("large.leiden.k", str(ld.n_clusters.iloc[0]), "seed 0")
    check("large.leiden.ari_seeds", f"{r3(ld.ari.min())} to {r3(ld.ari.max())}")
    check("large.leiden.k_seeds", f"{ld.n_clusters.min()} to {ld.n_clusters.max()}")
    check("large.leiden.default.k", str(ld_def.n_clusters))
    check("large.leiden.default.ari", r3(ld_def.ari))

    # The paper rule at seed 0: the resolution whose count is closest to 34, ARI breaking ties.
    s0 = ps[(ps.method == "PCA + Leiden") & (ps.seed == 0)].copy()
    s0["res"] = s0.setting.str.replace("resolution ", "").astype(float)
    s0 = s0[s0.res < 1.0]  # the default is reported apart
    s0["dist"] = (s0.n_clusters - 34).abs()
    tied = s0[s0.dist == s0.dist.min()].sort_values("ari", ascending=False)
    pick, runner = tied.iloc[0], tied.iloc[1]
    check("large.leiden.res", f"{pick.res:.2f}", f"{len(tied)} resolutions tie at distance {pick.dist}")
    check("large.leiden.runnerup.res", f"{runner.res:.2f}")
    check("large.leiden.runnerup.k", str(runner.n_clusters))
    check("large.leiden.runnerup.ari", r3(runner.ari))

    search = pd.read_csv(LS / "leiden_search.csv")
    check("large.leiden.n_runs", str(len(search)))
    check("large.leiden.n_res", str(search[search.seed == 0].resolution.nunique()))
    check("large.leiden.k_at_pilot_res", str(search[(search.seed == 0) & (search.resolution == 0.17)].n_clusters.iloc[0]))

    sa = pd.read_csv(LS / "seed_agreement.csv")
    check("large.kmeans.seed_ari_min", r3(sa[sa.method == "PCA + k-means"].ari.min()))
    check("large.leiden.seed_ari_min",
          r3(sa[(sa.method == "PCA + Leiden") & (sa.setting == "resolution 0.05")].ari.min()))

    met = pd.read_csv(LS / "metrics.csv")

    def m(method, data, setting=None):
        rows = met[(met.method == method) & (met.data == data)]
        if setting is not None:
            rows = rows[rows.setting == setting]
        return rows.iloc[0]

    check("large.sample.kmeans.ari", r3(m("PCA + k-means", "100k sample").ari))
    check("large.sample.leiden.ari", r3(m("PCA + Leiden", "100k sample").ari))
    check("large.sample.leiden.k", str(m("PCA + Leiden", "100k sample").n_clusters))
    check("large.sample.mica.ari", r3(m("MICA (graph mode)", "100k sample", "paper rule").ari))
    check("large.sample.mica.k", str(m("MICA (graph mode)", "100k sample", "paper rule").n_clusters))
    check("large.sample.mica_default.ari", r3(m("MICA (graph mode)", "100k sample", "resolution 1.8").ari))
    check("large.sample.mica_default.k", str(m("MICA (graph mode)", "100k sample", "resolution 1.8").n_clusters))
    check("large.assigned.kmeans.ari", r3(m("PCA + k-means", "sample + assigned").ari))
    check("large.assigned.leiden.ari", r3(m("PCA + Leiden", "sample + assigned").ari))
    check("large.assigned.mica.ari", r3(m("MICA (graph mode)", "sample + assigned", "paper rule").ari))
    check("large.kmeans.vs_sample", r3(m("PCA + k-means", "full data").ari_sample_cells))

    info = json.loads((LS / "assign_info.json").read_text())
    check("large.n_assigned", thousands(info["n_assigned"]))
    check("large.n_sample", thousands(info["n_reference"]))
    check("large.sample.mica.res", f"{info['mica']['resolution']}")
    check("large.unseen_pct", f"{100 * info['n_assigned'] / dims['n_cells']:.1f}")
    check("large.sample_pct", f"{100 * info['n_reference'] / dims['n_cells']:.1f}")
    check("large.assign.seconds", f"{info['knn_seconds'] + info['vote_seconds']['mica']:.0f}",
          "neighbour search + MICA vote")

    ho = pd.read_csv(LS / "assignment_holdout.csv")
    print(f"         (holdout: {ho.fold.nunique()} folds, mean accuracy over folds)")
    for key in ("kmeans", "leiden", "mica"):
        check(f"large.holdout.{key}.acc", f"{100 * ho[ho.partition == key].accuracy.mean():.1f}")

    # Time and memory. Elapsed time is the outer wall clock of each stage process.
    costs = pd.read_csv(LS / "stage_costs.csv").set_index("stage")
    through_leiden = ["manifest", "scan", "hvg", "normalize", "pca", "kmeans", "leiden"]
    cached_s = costs.loc[through_leiden, "wall_s_outer"].sum()
    validation_cold_read_s = 713  # first pass of the 100k validation, from PROGRESS.md Nov 19 (log not in git)
    disk_s = cached_s - costs.loc["scan", "wall_s"] + validation_cold_read_s
    check("large.full.minutes_cached", f"{cached_s / 60:.0f}", f"{cached_s:.0f} s, manifest to Leiden")
    check("large.full.minutes_disk", f"{disk_s / 60:.0f}", f"{disk_s:.0f} s, cached scan swapped for the 713 s cold read")
    check("large.read.minutes_cold", f"{validation_cold_read_s / 60:.0f}")
    check("large.max_process_gb", f"{costs.peak_process_gb.max():.1f}",
          f"stage {costs.peak_process_gb.idxmax()}")
    check("large.full.max_gb", f"{costs.peak_tree_gb.max():.0f}", f"stage {costs.peak_tree_gb.idxmax()}")
    prep = ["scan", "hvg", "normalize", "pca"]
    check("large.prep.minutes", f"{costs.loc[prep, 'wall_s'].sum() / 60:.0f}", "scan to PCA, cached")
    print(f"         (prep from disk, for the table caption's 17 min: "
          f"{(costs.loc[prep, 'wall_s'].sum() - costs.loc['scan', 'wall_s'] + validation_cold_read_s) / 60:.1f} min)")
    print(f"         (collect stage: wall_s {costs.loc['collect', 'wall_s']} against wall_s_outer "
          f"{costs.loc['collect', 'wall_s_outer']}, not used in the paper)")

    # Pilot values typed into make_tables.py, checked against the pilot's own metrics file.
    pilot = pd.read_csv(ROOT / "results" / "baselines_allen_100k_hvg2000_seurat_v3" / "metrics.csv").set_index("method")
    print(f"         (pilot k-means {pilot.loc['kmeans_pca', 'seconds']:.0f} s, Leiden "
          f"{pilot.loc['leiden_pca', 'seconds'] / 60:.0f} min, peak {pilot.loc['leiden_pca', 'peak_rss_gb']:.0f} GB: "
          f"table says 46 s, 17 min, 27 GB)")
    check("large.leiden.pilot_res", f"{pilot.loc['leiden_pca', 'resolution']:.2f}")


V1 = ROOT / "results" / "v1"
SMALL = ["buettner", "chung", "goolam", "klein", "kolodziejczyk", "pollen", "usoskin", "yan", "zeisel"]
CONFIRM = range(3, 10)


def load_v1(folder):
    df = pd.concat(pd.read_csv(f) for f in sorted((V1 / folder).glob("*.csv")))
    df["ds"] = df.dataset.str.replace("paper_", "", regex=False)
    return df


def readings(g):
    """Both readings, rebuilt from one run's stability curve with the definitions in Section IV.

    Finest stable: the threshold is the largest lower interval bound over the candidates, and the
    setting with the most clusters whose score clears it is kept. Most stable: the highest score,
    ties going to fewer clusters.
    """
    threshold = g.ci_low.max()
    passing = g[g.stability >= threshold]
    finest = passing.sort_values("n_clusters").iloc[-1]
    most = g.sort_values(["stability", "n_clusters"], ascending=[False, True]).iloc[0]
    return finest, most


def label_free():
    print("\nGroup B: label-free warning signs (Table checks, Discussion, Figure 1)")
    print("         Runs: nine small datasets, seeds 3 to 9, matched candidates, our scorer.")
    c = load_v1("curves")
    c = c[(c.scorer == "ours") & c.ds.isin(SMALL) & c.seed.isin(CONFIRM)]
    picks = load_v1("picks")
    picks = picks[(picks.scorer == "ours") & (picks.design == "mid")]

    rows, pick_mismatch = [], []
    for (ds, seed, method), g25 in c[c.in_mid & (c.b == 25)].groupby(["ds", "seed", "method"]):
        f0, m0 = readings(g25)
        f1, _ = readings(c[c.in_mid & (c.b == 100) & (c.ds == ds) & (c.seed == seed) & (c.method == method)])
        best = g25.stability.max()
        tied = g25[g25.stability >= best - 0.01]
        rows.append(dict(ds=ds, seed=seed, method=method, kF=f0.n_clusters, kM=m0.n_clusters, ariF=f0.ari,
                         ariM=m0.ari, kF100=f1.n_clusters, floor=g25.n_clusters.min(),
                         share=len(tied) / len(g25), tied_lo=tied.n_clusters.min(), tied_hi=tied.n_clusters.max()))
        # our picks against the ones the tagged runner wrote
        for b, reading, k in [(25, "stability_finest", f0.n_clusters), (25, "stability_most", m0.n_clusters),
                              (100, "stability_finest", f1.n_clusters)]:
            saved = picks[(picks.ds == ds) & (picks.seed == seed) & (picks.method == method) & (picks.b == b)
                          & (picks.reading == reading)].n_clusters
            if len(saved) != 1 or int(saved.iloc[0]) != int(k):
                pick_mismatch.append((ds, seed, method, b, reading, k, saved.tolist()))
    R = pd.DataFrame(rows)
    print(f"         Picks rebuilt from the curves differ from the runner's saved picks in "
          f"{len(pick_mismatch)} of {3 * len(R)} cases." + (f" First: {pick_mismatch[:3]}" if pick_mismatch else ""))

    check("L.runs", str(len(R)))
    gap = 5
    R["kgap"] = (R.kF - R.kM).abs()
    R["arigap"] = (R.ariF - R.ariM).abs()
    dis, agr = R[R.kgap >= gap], R[R.kgap < gap]
    check("L.gap", str(gap))
    check("L.disagree.n", str(len(dis)))
    check("L.disagree.mica", str((dis.method == "mica").sum()))
    check("L.disagree.arigap", r3(dis.arigap.mean()))
    check("L.agree.arigap", r3(agr.arigap.mean()))
    fl = R.kF == R.floor
    check("L.floor.F.n", str(fl.sum()))
    check("L.floor.F.ari", r3(R[fl].ariF.mean()))
    check("L.interior.F.ari", r3(R[~fl].ariF.mean()), f"{(~fl).sum()} interior runs")
    check("L.floor.M.n", str((R.kM == R.floor).sum()))

    def krange(s):
        lo, hi = int(s.min()), int(s.max())
        return str(lo) if lo == hi else f"{lo}\\textrm{{-}}{hi}"

    pm = R[(R.ds == "pollen") & (R.method == "mica")]
    check("L.pollen.mica.share", f"{100 * pm.share.median():.0f}", "median share within 0.01 of the best")
    check("L.pollen.mica.tiedlo", f"{pm.tied_lo.median():.0f}")
    check("L.pollen.mica.tiedhi", f"{pm.tied_hi.median():.0f}")
    check("L.pollen.mica.kF", krange(pm.kF))
    check("L.pollen.mica.kM", krange(pm.kM))
    check("L.pollen.mica.ariF", r3(pm.ariF.median()))
    check("L.pollen.mica.ariM", r3(pm.ariM.median()))
    facts = pd.read_csv(ROOT / "docs" / "dataset_facts.csv")
    print(f"         (Pollen types in docs/dataset_facts.csv: {facts.iloc[:, 5][facts.iloc[:, 0] == 'paper_pollen'].tolist()}, "
          f"paper says {PAPER.get('L.pollen.mica.truek')})")
    check("L.pollen.leiden.kM", krange(R[(R.ds == "pollen") & (R.method == "leiden_pca")].kM))

    bl = R[(R.ds == "buettner") & (R.method == "leiden_pca")]
    check("L.buettner.leiden.floorM", f"{(bl.kM == bl.floor).sum()} of {len(bl)}")
    check("L.buettner.leiden.ariF", r3(bl.ariF.median()))
    check("L.buettner.leiden.ariM", r3(bl.ariM.median()))

    # check 6: finest-stable pick at 25 against 100 subsamples
    for m, key in [("mica", "L.moved.mica"), ("leiden_pca", "L.moved.leiden")]:
        g = R[R.method == m]
        check(key, f"{(g.kF != g.kF100).sum()} of {len(g)}")
    # and chooseR's own code (scorer G) against our scorer at 100 subsamples, Leiden on Buettner
    ch = load_v1("picks_chooser")
    ch = ch[(ch.scorer == "chooseR_G") & (ch.design == "mid") & (ch.ds == "buettner") & (ch.method == "leiden_pca")
            & ch.seed.isin(CONFIRM)].set_index("seed").n_clusters
    b100 = bl.set_index("seed").kF100
    check("L.buettner.leiden.impl", f"{int((b100 != ch.reindex(b100.index)).sum())} of {len(b100)}")

    # check 1: on how many datasets Leiden's own grid never reaches 2 clusters
    og = c[c.in_original & (c.b == 25) & (c.method == "leiden_pca")]
    check("L.ownnever2.leiden", str(int((og.groupby("ds").n_clusters.min() > 2).sum())))


GROUPS = {"allen": allen, "labelfree": label_free}

if __name__ == "__main__":
    for name in sys.argv[1:] or GROUPS:
        GROUPS[name]()
    bad = [k for k, ok in RESULTS if not ok]
    print(f"\n{len(RESULTS) - len(bad)} of {len(RESULTS)} values match." + (f" Differ: {', '.join(bad)}" if bad else ""))
