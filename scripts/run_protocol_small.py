"""Protocol v1 on one small dataset and seed: one run at B=100, scored many ways (U-41).

Usage: python scripts/run_protocol_small.py <dataset> <seed> [b] [--require-tag]

The design is fixed in docs/protocol_matched_v1.md. Three methods on the same cells:
MICA in small-data mode, PCA plus Leiden, and PCA plus k-means (Jian's textbook baseline,
added in protocol v1). MICA and Leiden reuse scripts/run_rule_stability.py unchanged, so
the first 25 subsamples reproduce the original B=25 runs exactly. Leiden and k-means share
one PCA per subsample, computed on a fresh copy of those cells.

Leiden is swept over the original grid and the matched grid together, so every candidate
set can be scored afterwards from the saved partitions:
  original   each method's own grid, capped at 2 to 20 clusters (as first specified)
  mid        the matched design: counts all methods reach, the middle setting per count
  low, high  the first or last setting per count, to test whether the choice matters
  dense      every matched-grid setting that gives a matched count
MICA's k and k-means' k have one setting per count, so their low, high and dense equal mid.

Writes:
  data/runs/protocol_v1/<dataset>_seed<seed>/  every subsample partition (the store)
  results/v1/curves/<dataset>_seed<seed>.csv  every setting's score at B=25 and B=100
  results/v1/picks/<dataset>_seed<seed>.csv   the choice per design, method, B and reading
  results/v1/rule_stability_<exp>/ and rule_matched_stability_<exp>/
      the original protocol's B=25 output for MICA and Leiden, in the old layout, so seeds
      0 to 2 can be compared with the saved results (scripts/check_v1_regression.py)
  results/v1/provenance/<dataset>_seed<seed>.json
The labels are read only to score partitions after the choices are made.
"""
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_rule_stability as rrs  # noqa: E402  (reused unchanged: MICA's pipeline and the Leiden sweep)

from scminer_clust import grids, metrics, protocol, stability  # noqa: E402
from scminer_clust.paths import PROCESSED, RESULTS  # noqa: E402

K_RANGE = list(rrs.K_RANGE)
N_INIT = rrs.N_INIT  # 50 restarts, as for MICA's k-means
LEI_OLD = [float(r) for r in grids.LEIDEN_OLD]
LEI_M = [float(r) for r in grids.LEIDEN_MATCHED]
LEI_ALL = list(dict.fromkeys(LEI_OLD + LEI_M))
OUT = RESULTS / "v1"


def kmeans_parts(emb, ks, seed):
    return {k: KMeans(n_clusters=k, n_init=N_INIT, random_state=seed).fit_predict(emb) for k in ks}


def old_rows(name, method, picks, full, scores, truth, seed, seconds, n_clusterings):
    """The rows scripts/run_rule_stability.py writes for one method, in the same order and layout."""
    rows = []
    for rule, info in picks.items():
        info = {**info, "seconds": seconds, "n_clusterings": n_clusterings}
        rec = {"dataset": name, "method": method, "k_setting": rule, "seed": seed,
               "genes": "all", "hvg_flavor": "",
               "b_subsamples": stability.B_SUBSAMPLES, "fraction": stability.FRACTION, **info}
        rec.update(metrics.score(truth, full[info["chosen_setting"]]))
        rows.append(rec)
    for s, lab in full.items():
        rows.append({"dataset": name, "method": method, "k_setting": "full_sweep", "setting": s,
                     "seed": seed, "genes": "all", "n_clusters": len(np.unique(lab)),
                     "ari": metrics.score(truth, lab)["ari"]})
    for s, sc_ in scores.items():
        rows.append({"dataset": name, "method": method, "k_setting": "stability_curve", "setting": s,
                     "seed": seed, "genes": "all", **sc_, "ari": metrics.score(truth, full[s])["ari"]})
    return rows


def main(name, seed, b=protocol.B_FULL):
    seed, b = int(seed), int(b)
    t_start = time.time()
    adata = ad.read_h5ad(PROCESSED / f"{name}.h5ad")
    assert not [k for k in adata.obsm if k.startswith("X_pca")], "a cached full-data PCA would leak into subsamples"
    truth = adata.obs["label"].to_numpy()
    n = adata.n_obs
    X = rrs.dense(adata)
    D = rrs.mica_distance(name, X)
    print(f"{name}: {n} cells, {X.shape[1]} genes, seed {seed}, B={b}", flush=True)
    allc = np.arange(n)

    # ---------- full-data sweeps, each timed ----------
    tim = {}
    t = time.time()
    _, f_mica = rrs.mica_partitions(D, allc, seed)
    tim["full_mica"] = time.time() - t
    t = time.time()
    emb_full, f_lei = rrs.leiden_partitions(adata, allc, seed, LEI_ALL)
    tim["full_leiden"] = time.time() - t
    t = time.time()
    f_km = kmeans_parts(emb_full, K_RANGE, seed)
    tim["full_kmeans"] = time.time() - t
    assert not [k for k in adata.obsm if k.startswith("X_pca")], "the full-data sweep cached a PCA on adata"
    lei_old = {r: f_lei[r] for r in LEI_OLD}
    lei_m = {r: f_lei[r] for r in LEI_M}
    full = {"mica": f_mica, "leiden_pca": f_lei, "kmeans_pca": f_km}

    # ---------- candidate sets ----------
    common = grids.common_counts(f_mica, lei_m)
    common3 = grids.common_counts(f_mica, lei_m, f_km)
    assert common == common3, f"k-means changed the matched counts: {common} against {common3}"
    cand = {}
    for m, f in (("mica", f_mica), ("kmeans_pca", f_km)):
        one = list(grids.one_per_count(f, common).values())
        cand[m] = {"original": grids.original_settings(f), "mid": one, "low": one, "high": one, "dense": one}
    mid = list(grids.per_count(lei_m, common, "mid").values())
    assert mid == list(grids.one_per_count(lei_m, common).values())
    cand["leiden_pca"] = {"original": grids.original_settings(lei_old), "mid": mid,
                          "low": list(grids.per_count(lei_m, common, "low").values()),
                          "high": list(grids.per_count(lei_m, common, "high").values()),
                          "dense": [s for v in grids.dense_per_count(lei_m, common).values() for s in v]}
    need = {m: set().union(*map(set, cand[m].values())) for m in cand}
    sub_settings = {"mica": K_RANGE, "kmeans_pca": K_RANGE,
                    "leiden_pca": [r for r in LEI_ALL if r in need["leiden_pca"]]}
    print(f"  matched counts: {common} ({len(common)} of 19). Leiden settings per subsample: "
          f"{len(sub_settings['leiden_pca'])} (dense {len(cand['leiden_pca']['dense'])})", flush=True)

    # ---------- subsamples ----------
    subs = stability.subsamples(n, b=b, seed=seed)
    size = len(subs[0])
    codes = {m: np.empty((b, len(sub_settings[m]), size), dtype=np.int16) for m in sub_settings}
    t_rep = {m: np.zeros(b) for m in sub_settings}
    for i, cells in enumerate(subs):
        t = time.time()
        _, pm = rrs.mica_partitions(D, cells, seed)
        t_rep["mica"][i] = time.time() - t
        t = time.time()
        emb, pl = rrs.leiden_partitions(adata, cells, seed, sub_settings["leiden_pca"])
        t_rep["leiden_pca"][i] = time.time() - t
        t = time.time()
        pk = kmeans_parts(emb, K_RANGE, seed)
        t_rep["kmeans_pca"][i] = time.time() - t
        for m, parts in (("mica", pm), ("leiden_pca", pl), ("kmeans_pca", pk)):
            for j, s in enumerate(sub_settings[m]):
                codes[m][i, j] = protocol.encode(parts[s])
        if (i + 1) % 10 == 0:
            print(f"    subsample {i + 1} of {b} ({(time.time() - t_start) / 60:.1f} min)", flush=True)
    assert not [k for k in adata.obsm if k.startswith("X_pca")]

    # ---------- the store ----------
    store = protocol.STORE / f"{name}_seed{seed}"
    protocol.save_store(store, np.stack(subs), {m: {"settings": sub_settings[m], "full": full[m], "codes": codes[m]}
                                                  for m in sub_settings})
    meta = {"dataset": name, "seed": seed, "b": b, "n": n, "size": size, "common_counts": common,
            "candidates": cand, "sub_settings": sub_settings, "methods": list(sub_settings)}
    protocol.write_json(store / "meta.json", meta)

    # ---------- scores: our scorer at B=25 and B=100 ----------
    ref = stability.reference_cells(n, seed=seed)
    t = time.time()
    counters = {bb: protocol.Counter(n, ref, subs, bb) for bb in (protocol.B_OLD, b)}
    scores, t_score = {}, {m: 0.0 for m in sub_settings}
    curve_rows = []
    for m in sub_settings:
        for j, s in enumerate(sub_settings[m]):
            lab_ref = np.asarray(full[m][s])[ref]
            sc_full = metrics.score(truth, full[m][s])
            for bb, ctr in counters.items():
                ts = time.time()
                tog = ctr.together(codes[m][:bb, j])
                sc_ = stability.score_setting(protocol.distance_ours(tog, ctr.seen), lab_ref, seed=seed)
                t_score[m] += time.time() - ts
                scores[(m, bb, s)] = sc_
                curve_rows.append({"dataset": name, "seed": seed, "method": m, "setting": s, "b": bb,
                                   "scorer": "ours", **sc_, "ari": sc_full["ari"], "ami": sc_full["ami"],
                                   **{f"in_{d}": s in cand[m][d] for d in protocol.DESIGNS}})
    tim["scoring"] = time.time() - t

    # ---------- picks per design, method, B and reading ----------
    pick_rows = []
    for m in sub_settings:
        for d in protocol.DESIGNS:
            for bb in counters:
                sc = {s: scores[(m, bb, s)] for s in cand[m][d]}
                for reading, info in stability.pick(sc).items():
                    s = info["chosen_setting"]
                    sf = metrics.score(truth, full[m][s])
                    pick_rows.append({"dataset": name, "seed": seed, "method": m, "design": d, "b": bb,
                                      "scorer": "ours", "reading": reading, "chosen_setting": s,
                                      "n_clusters": info["n_clusters"], "threshold": info["threshold"],
                                      "n_candidates": info["n_candidates"], "n_passing": info["n_passing"],
                                      "stability": info["stability"], "ari": sf["ari"], "ami": sf["ami"]})
    tag = f"{name}_seed{seed}"
    for sub, rows in (("curves", curve_rows), ("picks", pick_rows)):
        (OUT / sub).mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(OUT / sub / f"{tag}.csv", index=False)

    # ---------- the original protocol's output, B=25, MICA and Leiden, old layout ----------
    exp_tag = name + (f"_seed{seed}" if seed != rrs.SEED else "")
    per_set = {m: np.sum(t_rep[m][:protocol.B_OLD]) / len(sub_settings[m]) for m in sub_settings}
    for design, folder, lei_full, grid in (("original", "rule_stability", lei_old, grids.LEIDEN_OLD),
                                           ("mid", "rule_matched_stability", lei_m, grids.LEIDEN_MATCHED)):
        rows = []
        for m, fm in (("mica", f_mica), ("leiden_pca", lei_full)):
            settings = cand[m][design]
            sc = {s: scores[(m, protocol.B_OLD, s)] for s in settings}
            # Seconds are reconstructed from this run's timings: the full sweep's share for this
            # grid, the first 25 subsamples' share for these settings, and their scoring.
            share = len(fm) / len(full[m])
            secs = (tim[f"full_{'leiden' if m == 'leiden_pca' else 'mica'}"] * share + per_set[m] * len(settings)
                    + t_score[m] / (2 * len(sub_settings[m])) * len(settings))
            rows += old_rows(name, m, stability.pick(sc), fm, sc, truth, seed, round(secs, 1),
                             len(settings) * (protocol.B_OLD + 1))
        metrics.write(f"v1/{folder}_{exp_tag}", rows,
                      {"dataset": name, "rule": "chooseR-style subsampling stability", "n_hvg": None,
                       "grid": "matched" if design == "mid" else "original",
                       "b_subsamples": stability.B_SUBSAMPLES, "fraction": stability.FRACTION,
                       "k_range": [2, 20], "resolutions": grid.tolist(), "seed": seed})

    tim["total"] = time.time() - t_start
    prov = protocol.provenance({"dataset": name, "seed": seed, "b": b, "timings": tim,
                                "seconds_per_subsample": {m: float(np.mean(v)) for m, v in t_rep.items()},
                                "store": str(store), "input_sha256": protocol.sha256(PROCESSED / f"{name}.h5ad"),
                                "seconds_note": "old-layout seconds are reconstructed from this run's timings"})
    protocol.write_json(OUT / "provenance" / f"{tag}.json", prov)
    print(f"  done in {tim['total'] / 60:.1f} min. "
          + ", ".join(f"{r['method']} {r['design']} B={r['b']} {r['reading'][10:]}: k={r['n_clusters']} ARI {r['ari']:.3f}"
                      for r in pick_rows if r["design"] in ("original", "mid") and r["b"] == protocol.B_OLD), flush=True)
    return 0


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--require-tag"]
    if "--require-tag" in sys.argv:
        protocol.require_tag()
    sys.exit(main(*args))
