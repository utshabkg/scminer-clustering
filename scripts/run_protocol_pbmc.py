"""Protocol v1 on PBMC14K, MICA's graph mode, one seed, in stages (U-42).

Usage: python scripts/run_protocol_pbmc.py <seed> <stage> [b] [--workers N] [--require-tag]
Stages, run in this order (mica and pca can run side by side once full is done):
  full   full-data sweeps for MICA (graph mode, every matched Louvain value), Leiden (original
         and matched grids) and PCA plus k-means, then the candidate sets. Writes meta.json.
  mica   MICA's graph pipeline on each of the b subsamples, one file per subsample, so a
         stopped run resumes where it left off. node2vec embeddings are kept.
  pca    Leiden and k-means on each subsample, sharing one PCA, also one file per subsample.
  score  every scorer on the saved partitions, the picks, and the old-layout B=25 output.

The design is the small-data runner's (scripts/run_protocol_small.py), with two differences.
MICA's node2vec and Louvain come from MICA's own code, which takes no seed, so MICA's
partitions are new draws even at seeds 0 to 2. And a dense Louvain set for MICA would cost
about 73 CPU-hours per seed, so MICA gets the middle, lowest and highest setting per count,
not all of them. Leiden gets all four variants.

Scorers written by the score stage (our scorer, stability.score_setting):
  S0  B=25 on the 3,000-cell reference subset, the original protocol
  S1  B=100 on the same reference subset
  S2  B=100 on all 13,605 cells (original, mid, low and high designs, not dense)
chooseR's own R code (S3) is applied to the same store by scripts/run_chooser_r.py.
"""
import os
import shutil
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.cluster import KMeans

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_rule_stability_graph as rrsg  # noqa: E402  (MICA's graph pipeline, reused unchanged)

from scminer_clust import cluster, grids, metrics, protocol, stability  # noqa: E402
from scminer_clust.paths import PROCESSED, RESULTS  # noqa: E402

NAME = rrsg.NAME
K_RANGE = list(range(2, 21))
N_INIT = 50
MICA_OLD = [float(r) for r in np.round(np.arange(0.1, grids.MICA_GRAPH_MAX_OLD + 1e-9, 0.1), 1)]
MICA_ALL = [float(r) for r in grids.MICA_GRAPH_MATCHED]
LEI_OLD = [float(r) for r in grids.LEIDEN_OLD]
LEI_M = [float(r) for r in grids.LEIDEN_MATCHED]
LEI_ALL = list(dict.fromkeys(LEI_OLD + LEI_M))
OUT = RESULTS / "v1"
FULL_DESIGNS = ("original", "mid", "low", "high")  # scored on all cells too (S2)


def store_dir(seed):
    return protocol.STORE / f"{NAME}_seed{seed}"


def load():
    a = ad.read_h5ad(PROCESSED / f"{NAME}.h5ad")
    assert not [k for k in a.obsm if k.startswith("X_pca")]
    X = np.asarray(a.X.todense() if sparse.issparse(a.X) else a.X, dtype=np.float32)
    return a, X, a.obs["label"].to_numpy(), a.n_obs


def pca_partitions(a, cells, seed, resolutions, ks):
    """rrsg.leiden_partitions plus k-means on the same principal components."""
    g, emb = cluster._graph(a[np.asarray(cells)].copy(), seed=seed)
    lei = {}
    for r in resolutions:
        sc.tl.leiden(g, resolution=float(r), random_state=seed, flavor="igraph", n_iterations=2, directed=False)
        lei[float(r)] = g.obs["leiden"].to_numpy().astype(str)
    km = {k: KMeans(n_clusters=k, n_init=N_INIT, random_state=seed).fit_predict(emb) for k in ks}
    return lei, km


def stage_full(seed, b, workers, fresh=False):
    st = store_dir(seed)
    if (st / "meta.json").exists():
        # MICA's full sweep is a new draw every time, so the subsample files of an earlier full
        # stage are never valid after a rerun. Refuse, or with --fresh clear them first.
        if not fresh:
            sys.exit(f"{st} already has a meta.json. Rerun with --fresh to clear it and start over.")
        (st / "meta.json").unlink()  # first, so a waiting mica or pca stage keeps waiting
        for d in ("mica_sub", "pca_sub", "emb"):
            shutil.rmtree(st / d, ignore_errors=True)
        for f in list(st.glob("*.npy")) + list(st.glob("*.csv")) + [st / "scored.json"]:
            f.unlink(missing_ok=True)
    a, X, truth, n = load()
    D = rrsg.full_distance(X)
    allc = np.arange(n)
    tim = {}
    t = time.time()
    f_mica = rrsg.mica_graph_partitions(D, allc, f"v1s{seed}_full", workers=workers, resolutions=MICA_ALL,
                                        keep=store_dir(seed) / "emb" / "full")
    f_mica = {r: f_mica[r] for r in MICA_ALL}
    tim["full_mica"] = time.time() - t
    t = time.time()
    f_lei, f_km = pca_partitions(a, allc, seed, LEI_ALL, K_RANGE)
    tim["full_pca"] = time.time() - t
    mica_old = {r: f_mica[r] for r in MICA_OLD}
    lei_old = {r: f_lei[r] for r in LEI_OLD}
    lei_m = {r: f_lei[r] for r in LEI_M}
    common = grids.common_counts(f_mica, lei_m)
    assert common == grids.common_counts(f_mica, lei_m, f_km), "k-means changed the matched counts"
    cand = {"mica": {"original": grids.original_settings(mica_old)},
            "leiden_pca": {"original": grids.original_settings(lei_old)},
            "kmeans_pca": {}}
    for which in ("mid", "low", "high"):
        cand["mica"][which] = list(grids.per_count(f_mica, common, which).values())
        cand["leiden_pca"][which] = list(grids.per_count(lei_m, common, which).values())
    cand["mica"]["dense"] = cand["mica"]["mid"]  # not affordable in graph mode, see the docstring
    cand["leiden_pca"]["dense"] = [s for v in grids.dense_per_count(lei_m, common).values() for s in v]
    one = list(grids.one_per_count(f_km, common).values())
    cand["kmeans_pca"] = {"original": grids.original_settings(f_km), "mid": one, "low": one, "high": one, "dense": one}
    need = {m: set().union(*map(set, cand[m].values())) for m in cand}
    sub_settings = {"mica": [r for r in MICA_ALL if r in need["mica"]],
                    "leiden_pca": [r for r in LEI_ALL if r in need["leiden_pca"]], "kmeans_pca": K_RANGE}
    full = {"mica": f_mica, "leiden_pca": f_lei, "kmeans_pca": f_km}
    st.mkdir(parents=True, exist_ok=True)
    for m in full:
        np.save(st / f"{m}_full.npy", np.stack([np.asarray(full[m][s]) for s in sub_settings[m]]))
        np.save(st / f"{m}_settings.npy", np.asarray(sub_settings[m], dtype=np.float64))
        # every setting's full-data count and ARI, for the old-layout full_sweep rows and the transfer
        pd.DataFrame([{"method": m, "setting": s, "n_clusters": len(np.unique(lab)), "ari": metrics.score(truth, lab)["ari"]}
                      for s, lab in full[m].items()]).to_csv(st / f"{m}_full_sweep.csv", index=False)
    subs = stability.subsamples(n, b=b, seed=seed)
    np.save(st / "cells.npy", np.stack(subs).astype(np.int32))
    meta = {"dataset": NAME, "seed": seed, "b": b, "n": n, "size": len(subs[0]), "common_counts": common,
            "candidates": cand, "sub_settings": sub_settings, "methods": list(full), "timings": tim}
    protocol.write_json(st / "meta.json", meta)
    print(f"  matched counts {common}. Settings per subsample: "
          + ", ".join(f"{m} {len(v)}" for m, v in sub_settings.items()), flush=True)


def read_meta(seed):
    import json
    p = store_dir(seed) / "meta.json"
    while not p.exists():  # the pca and mica stages may be started before full has finished
        print("  waiting for the full stage", flush=True)
        time.sleep(120)
    return json.loads(p.read_text())


def stage_subsamples(seed, kind, workers):
    meta = read_meta(seed)
    st = store_dir(seed)
    cells_all = np.load(st / "cells.npy")
    if kind == "mica":
        a, X, truth, n = load()
        D = rrsg.full_distance(X)
    else:
        a = ad.read_h5ad(PROCESSED / f"{NAME}.h5ad")
    out = st / f"{kind}_sub"
    out.mkdir(exist_ok=True)
    t0 = time.time()
    for i, cells in enumerate(cells_all):
        f = out / f"{i:03d}.npz"
        if f.exists():
            if sub_file_ok(f, meta, kind, cells):
                continue
            print(f"    {f.name} was made for other settings or cells, recomputing", flush=True)
        t = time.time()
        if kind == "mica":
            parts = rrsg.mica_graph_partitions(D, cells, f"v1s{seed}_sub{i + 1}", workers=workers,
                                               resolutions=meta["sub_settings"]["mica"], keep=st / "emb" / f"sub{i + 1:03d}")
            arrays = {"mica": np.stack([protocol.encode(parts[float(s)]) for s in meta["sub_settings"]["mica"]])}
        else:
            lei, km = pca_partitions(a, cells, seed, meta["sub_settings"]["leiden_pca"], K_RANGE)
            arrays = {"leiden_pca": np.stack([protocol.encode(lei[float(s)]) for s in meta["sub_settings"]["leiden_pca"]]),
                      "kmeans_pca": np.stack([protocol.encode(km[int(k)]) for k in K_RANGE])}
        tmp = out / f"{i:03d}.tmp.npz"  # written whole, then renamed, so a stopped run leaves no half file
        np.savez(tmp, cells=np.asarray(cells), seconds=time.time() - t, **arrays,
                 **{f"settings_{m}": np.asarray(meta["sub_settings"][m], dtype=np.float64) for m in arrays})
        os.replace(tmp, f)
        print(f"    {kind} subsample {i + 1} of {len(cells_all)} in {time.time() - t:.0f} s "
              f"({(time.time() - t0) / 3600:.2f} h)", flush=True)


def sub_file_ok(f, meta, kind, cells):
    """A saved subsample file is reused only if it was made for these settings and these cells."""
    try:
        z = np.load(f)
        methods = ["mica"] if kind == "mica" else ["leiden_pca", "kmeans_pca"]
        return np.array_equal(z["cells"], np.asarray(cells)) and all(
            np.array_equal(z[f"settings_{m}"], np.asarray(meta["sub_settings"][m], dtype=np.float64)) for m in methods)
    except Exception:
        return False


def codes_for(seed, method, b):
    kind = "mica" if method == "mica" else "pca"
    meta = read_meta(seed)
    cells_all = np.load(store_dir(seed) / "cells.npy")
    files = [store_dir(seed) / f"{kind}_sub" / f"{i:03d}.npz" for i in range(b)]
    missing = [f.name for f in files if not f.exists()]
    assert not missing, f"{method}: {len(missing)} subsamples missing, first {missing[:3]}"
    bad = [f.name for i, f in enumerate(files) if not sub_file_ok(f, meta, kind, cells_all[i])]
    assert not bad, f"{method}: {len(bad)} subsample files do not match meta.json, first {bad[:3]}"
    secs = [float(np.load(f)["seconds"]) for f in files]
    return np.stack([np.load(f)[method] for f in files]), secs


def old_rows(method, picks, full_d, scores, truth, seed, b, seconds):
    """The rows scripts/run_rule_stability_graph.py writes for one method, same order and layout."""
    rows = []
    for rule, info in picks.items():
        rec = {"dataset": NAME, "method": method, "k_setting": rule, "seed": seed, "genes": "all",
               "mode": "graph" if method == "mica" else "pca", "b_subsamples": b, "seconds": seconds, **info}
        rec.update(metrics.score(truth, full_d[info["chosen_setting"]]))
        rows.append(rec)
    for s, sc_ in scores.items():
        rows.append({"dataset": NAME, "method": method, "k_setting": "stability_curve", "setting": s,
                     "seed": seed, **sc_, "ari": metrics.score(truth, full_d[s])["ari"]})
    for s, lab in full_d.items():
        rows.append({"dataset": NAME, "method": method, "k_setting": "full_sweep", "setting": s, "seed": seed,
                     "n_clusters": len(np.unique(lab)), "ari": metrics.score(truth, lab)["ari"]})
    return rows


def stage_score(seed, b):
    t_start = time.time()
    meta = read_meta(seed)
    st = store_dir(seed)
    truth = ad.read_h5ad(PROCESSED / f"{NAME}.h5ad").obs["label"].to_numpy()
    n = meta["n"]
    cells = np.load(st / "cells.npy")[:b]
    cand = {m: {d: [float(s) for s in v] for d, v in meta["candidates"][m].items()} for m in meta["candidates"]}
    sub_settings = {m: [float(s) for s in v] for m, v in meta["sub_settings"].items()}
    full, codes, secs = {}, {}, {}
    for m in sub_settings:
        fl = np.load(st / f"{m}_full.npy")
        full[m] = {s: fl[j] for j, s in enumerate(sub_settings[m])}
        codes[m], secs[m] = codes_for(seed, m, b)
        np.save(st / f"{m}_codes.npy", codes[m])  # b x S x size, the layout scripts/run_chooser_r.py reads
    ref = stability.reference_cells(n, seed=seed)
    allc = np.arange(n)
    b_old = min(protocol.B_OLD, b)  # below 25 only in the pre-tag smoke test
    counters = {("ref", b_old): protocol.Counter(n, ref, cells, b_old),
                ("ref", b): protocol.Counter(n, ref, cells, b),
                ("all", b): protocol.Counter(n, allc, cells, b)}
    scores, curve_rows = {}, []
    for m in sub_settings:
        in_full = set().union(*(set(cand[m][d]) for d in FULL_DESIGNS))
        for j, s in enumerate(sub_settings[m]):
            sf = metrics.score(truth, full[m][s])
            for (where, bb), ctr in counters.items():
                if where == "all" and s not in in_full:
                    continue
                lab = np.asarray(full[m][s])[ctr.ref]
                tog = ctr.together(codes[m][:bb, j])
                sc_ = stability.score_setting(protocol.distance_ours(tog, ctr.seen), lab, seed=seed)
                scores[(m, where, bb, s)] = sc_
                curve_rows.append({"dataset": NAME, "seed": seed, "method": m, "setting": s, "b": bb,
                                   "scorer": f"ours_{where}", **sc_, "ari": sf["ari"], "ami": sf["ami"],
                                   **{f"in_{d}": s in cand[m][d] for d in protocol.DESIGNS}})
        print(f"  scored {m} ({(time.time() - t_start) / 60:.0f} min)", flush=True)
    pick_rows = []
    for m in sub_settings:
        for d in protocol.DESIGNS:
            for (where, bb) in counters:
                if where == "all" and d not in FULL_DESIGNS:
                    continue
                sc_d = {s: scores[(m, where, bb, s)] for s in cand[m][d]}
                for reading, info in stability.pick(sc_d).items():
                    sf = metrics.score(truth, full[m][info["chosen_setting"]])
                    pick_rows.append({"dataset": NAME, "seed": seed, "method": m, "design": d, "b": bb,
                                      "scorer": f"ours_{where}", "reading": reading,
                                      "chosen_setting": info["chosen_setting"], "n_clusters": info["n_clusters"],
                                      "threshold": info["threshold"], "n_candidates": info["n_candidates"],
                                      "n_passing": info["n_passing"], "stability": info["stability"],
                                      "ari": sf["ari"], "ami": sf["ami"]})
    protocol.write_json(st / "scored.json", {"b": b})  # the codes saved above hold b subsamples
    tag = f"{NAME}_seed{seed}"
    for sub, rows in (("curves", curve_rows), ("picks", pick_rows)):
        (OUT / sub).mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(OUT / sub / f"{tag}.csv", index=False)
    # the original protocol's output (S0), MICA and Leiden, in the graph runner's layout
    sweeps = {m: pd.read_csv(st / f"{m}_full_sweep.csv") for m in ("mica", "leiden_pca")}
    for design, folder, grid in (("original", "rule_stability", (MICA_OLD, LEI_OLD)),
                                 ("mid", "rule_matched_stability", (MICA_ALL, LEI_M))):
        rows = []
        for m, g in zip(("mica", "leiden_pca"), grid):
            settings = cand[m][design]
            sc_d = {s: scores[(m, "ref", b_old, s)] for s in settings}
            # the full-data labels are stored only for settings clustered on subsamples, so the
            # full_sweep rows come from the saved sweep table, the other rows from the labels
            full_d = {s: full[m][s] for s in settings}
            r = old_rows(m, stability.pick(sc_d), full_d, sc_d, truth, seed, b_old,
                         round(sum(secs[m][:b_old]) * len(settings) / len(sub_settings[m]), 1))
            r = [x for x in r if x["k_setting"] != "full_sweep"]
            sw = sweeps[m][sweeps[m].setting.round(6).isin(np.round(g, 6))]
            r += [{"dataset": NAME, "method": m, "k_setting": "full_sweep", "setting": float(x.setting), "seed": seed,
                   "n_clusters": int(x.n_clusters), "ari": float(x.ari)} for x in sw.itertuples()]
            rows += r
        metrics.write(f"v1/{folder}_{NAME}_graph_seed{seed}", rows,
                      {"dataset": NAME, "rule": "chooseR-style subsampling stability", "mica_mode": "graph",
                       "grid": "matched" if design == "mid" else "original", "b_subsamples": b_old,
                       "fraction": stability.FRACTION, "seed": seed,
                       "mica_resolutions": list(grid[0]), "leiden_resolutions": list(grid[1])})
    prov = protocol.provenance({"dataset": NAME, "seed": seed, "b": b, "timings": meta.get("timings"),
                                "seconds_per_subsample": {m: float(np.mean(v)) for m, v in secs.items()},
                                "score_minutes": (time.time() - t_start) / 60, "store": str(st),
                                "input_sha256": protocol.sha256(PROCESSED / f"{NAME}.h5ad")})
    protocol.write_json(OUT / "provenance" / f"{tag}.json", prov)
    print(f"  done: {len(pick_rows)} picks in {(time.time() - t_start) / 60:.0f} min", flush=True)


def main(seed, stage, b=protocol.B_FULL, workers=8, fresh=False):
    seed, b, workers = int(seed), int(b), int(workers)
    print(f"{NAME} seed {seed} stage {stage} B={b} workers {workers}", flush=True)
    if stage == "full":
        stage_full(seed, b, workers, fresh)
    elif stage in ("mica", "pca"):
        stage_subsamples(seed, stage, workers)
    elif stage == "score":
        stage_score(seed, b)
    else:
        sys.exit(f"unknown stage {stage}")
    return 0


if __name__ == "__main__":
    argv = sys.argv[1:]
    if "--require-tag" in argv:
        argv.remove("--require-tag")
        protocol.require_tag()
    fresh = "--fresh" in argv
    argv = [x for x in argv if x != "--fresh"]
    workers = 8
    if "--workers" in argv:
        i = argv.index("--workers")
        workers = argv[i + 1]
        del argv[i:i + 2]
    sys.exit(main(*argv, workers=workers, fresh=fresh))
