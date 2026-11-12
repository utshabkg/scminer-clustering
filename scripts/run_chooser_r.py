"""Score protocol v1's saved partitions with chooseR's own R code, scorer S3 (U-43).

Usage: python scripts/run_chooser_r.py <dataset> <seed> [--designs original,mid] [--workers 8]
                                       [--tally-max 704] [--crosscheck N]

For every candidate setting in the chosen designs, the co-clustering over all cells and all
100 subsamples goes to scripts/score_chooser_r.R, which applies chooseR's code unchanged:
silhouette in 1 minus the frequency, the per-cluster average, and boot_median (BCa, 25,000
replicates) over those averages. On datasets of up to --tally-max cells, R also builds the
frequency itself from the raw subsample labels with chooseR's find_matches and percent_match,
and the result is compared with our own float64 frequency (protocol.freq_chooser). With
--crosscheck N, the first N settings of larger datasets are tallied in R too.

chooseR's choice (protocol.pick_chooser) is then applied per design, three ways for settings
whose BCa interval cannot be computed (boot.ci stops with an error, or returns no interval
because every bootstrap median is equal), as declared in the protocol:
  V   as chooseR's code behaves: an error anywhere stops the run, so no choice is made
  G   (primary) a failed setting gets a point interval at its median
  G2  failed settings are dropped
Writes results/v1/chooser_r/<dataset>_seed<seed>.csv (one row per setting),
results/v1/picks_chooser/<dataset>_seed<seed>.csv and results/v1/chooser_r/crosscheck_<...>.csv.
"""
import json
import subprocess
import sys
import time
import zlib
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from scminer_clust import metrics, protocol
from scminer_clust.paths import DATA, PROCESSED, REPO, RESULTS

CHOOSER = DATA / "external" / "chooseR_410bdbfaebcbc57df45176491545a1165e3c6fc5" / "R" / "pipeline.R"
RSCRIPT = Path.home() / "anaconda3" / "envs" / "chooser-r" / "bin" / "Rscript"
SHM = Path("/dev/shm/scm_chooser")
OUT = RESULTS / "v1"


def args_of(argv):
    opts = {"designs": "original,mid", "workers": 8, "tally-max": 704, "crosscheck": 0}
    pos = []
    i = 0
    while i < len(argv):
        if argv[i].startswith("--"):
            opts[argv[i][2:]] = argv[i + 1]
            i += 2
        else:
            pos.append(argv[i])
            i += 1
    return pos, opts


def run_batch(jobs, workers, tag):
    """jobs: list of manifest rows whose files exist. Runs them in `workers` R processes."""
    chunks = [jobs[i::workers] for i in range(workers) if jobs[i::workers]]
    procs, outs = [], []
    for c, chunk in enumerate(chunks):
        man, out = SHM / f"{tag}_man{c}.csv", SHM / f"{tag}_out{c}.csv"
        pd.DataFrame(chunk).to_csv(man, index=False)
        procs.append(subprocess.Popen([str(RSCRIPT), str(REPO / "scripts" / "score_chooser_r.R"), str(CHOOSER),
                                       str(man), str(out)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True))
        outs.append(out)
    frames = []
    for p, out in zip(procs, outs):
        so, se = p.communicate()
        if p.returncode != 0:
            raise RuntimeError(f"R failed: {se[-1500:]}")
        frames.append(pd.read_csv(out))
    return pd.concat(frames, ignore_index=True)


def main(ds, seed, designs="original,mid", workers=8, tally_max=704, crosscheck=0):
    seed, workers, tally_max, crosscheck = int(seed), int(workers), int(tally_max), int(crosscheck)
    designs = designs.split(",")
    t0 = time.time()
    st = protocol.STORE / f"{ds}_seed{seed}"
    meta = json.loads((st / "meta.json").read_text())
    n = meta["n"]
    scored = st / "scored.json"  # PBMC14K: the score stage records how many subsamples it saved
    b = json.loads(scored.read_text())["b"] if scored.exists() else meta["b"]
    cells = np.load(st / "cells.npy")[:b]
    truth = ad.read_h5ad(PROCESSED / f"{ds}.h5ad").obs["label"].to_numpy()
    ctr = protocol.Counter(n, np.arange(n), cells, b)
    SHM.mkdir(parents=True, exist_ok=True)
    run_tag = f"{ds}_s{seed}"
    jobs, info = [], {}
    for m in meta["methods"]:
        _, codes, settings, full = protocol.load_store(st, m)
        assert codes.shape[0] == b, f"{m}: {codes.shape[0]} subsamples saved, expected {b}"
        want = set().union(*(set(map(float, meta["candidates"][m][d])) for d in designs))
        for j, s in enumerate(settings):
            if float(s) in want:
                jid = f"{run_tag}_{m}_{s:g}"
                info[jid] = (m, j, float(s), codes, full[j])
                jobs.append(jid)
    print(f"{ds} seed {seed}: {len(jobs)} settings to score with chooseR's code, {workers} R workers", flush=True)
    # chooseR's percent_match divides by a fixed 100, so R tallies only at b = 100. With fewer
    # subsamples (the declared PBMC14K fallback) matrix mode passes together / seen, which is
    # percent_match with n = b.
    tallied = (set(jobs) if n <= tally_max else set(jobs[:crosscheck])) if b == 100 else set()
    batch = max(workers, 400 if n <= 3000 else workers)
    res, checks = [], []
    for start in range(0, len(jobs), batch):
        rows, freqs = [], {}
        for jid in jobs[start:start + batch]:
            m, j, s, codes, lab = info[jid]
            base = SHM / jid
            np.asarray(lab).astype(np.int64).astype(np.int32).tofile(f"{base}_lab.bin")
            tog = ctr.together(codes[:, j])
            freq = protocol.freq_chooser(tog, ctr.seen)
            row = {"job_id": jid, "m": n, "b": b, "setting": s, "labels_file": f"{base}_lab.bin",
                   "r_seed": zlib.crc32(jid.encode()) % (2 ** 31 - 1), "dump_file": ""}
            if jid in tallied:
                tab = np.full((n, b), -1, dtype=np.int32)
                for t in range(b):
                    tab[cells[t], t] = codes[t, j]
                tab.flatten(order="F").tofile(f"{base}_tab.bin")
                row.update(mode="tally", data_file=f"{base}_tab.bin", dump_file=f"{base}_dump.bin")
                freqs[jid] = freq
            else:
                freq.astype(np.float64).tofile(f"{base}_freq.bin")
                row.update(mode="matrix", data_file=f"{base}_freq.bin")
            rows.append(row)
        out = run_batch(rows, workers, f"{run_tag}_{start}")
        for jid, freq in freqs.items():
            r = np.fromfile(SHM / f"{jid}_dump.bin", dtype=np.float64).reshape(n, n, order="F")
            checks.append({"job_id": jid, "max_abs_diff": float(np.max(np.abs(r - freq))),
                           "equal": bool(np.array_equal(r, freq))})
        res.append(out)
        for f in SHM.glob(f"{run_tag}_*"):
            f.unlink()
        print(f"  {min(start + batch, len(jobs))} of {len(jobs)} scored ({(time.time() - t0) / 60:.1f} min)", flush=True)
    res = pd.concat(res, ignore_index=True)
    res["method"] = [info[j][0] for j in res.job_id]
    res["setting"] = [info[j][2] for j in res.job_id]
    (OUT / "chooser_r").mkdir(parents=True, exist_ok=True)
    res.to_csv(OUT / "chooser_r" / f"{ds}_seed{seed}.csv", index=False)
    if checks:
        pd.DataFrame(checks).to_csv(OUT / "chooser_r" / f"crosscheck_{ds}_seed{seed}.csv", index=False)

    # chooseR's choice per design, three ways for failed bootstraps
    picks = []
    by = {(r.method, r.setting): r for r in res.itertuples()}
    for m in meta["methods"]:
        _, _, settings, full = protocol.load_store(st, m)
        col = {float(s): j for j, s in enumerate(settings)}
        for d in designs:
            cand = [float(s) for s in meta["candidates"][m][d]]
            for variant in ("V", "G", "G2"):
                rows = []
                for s in cand:
                    r = by[(m, s)]
                    failed = r.status in ("error", "no_interval")
                    if failed and variant == "G2":
                        continue
                    low = None if (failed and variant == "V") else (r.med_direct if failed else r.low_med)
                    rows.append({"setting": s, "n_clusters": int(r.n_clusters), "med": float(r.med_direct),
                                 "low_med": low})
                chosen, threshold = protocol.pick_chooser(rows)
                rec = {"dataset": ds, "seed": seed, "method": m, "design": d, "b": b, "scorer": f"chooseR_{variant}",
                       "reading": "stability_finest", "chosen_setting": chosen, "threshold": threshold,
                       "n_candidates": len(cand), "n_failed": int(sum(by[(m, s)].status in ("error", "no_interval") for s in cand))}
                if chosen is not None:
                    lab = full[col[chosen]]
                    sf = metrics.score(truth, lab)
                    rec.update(n_clusters=sf["n_clusters"], ari=sf["ari"], ami=sf["ami"])
                picks.append(rec)
    (OUT / "picks_chooser").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(picks).to_csv(OUT / "picks_chooser" / f"{ds}_seed{seed}.csv", index=False)
    rv = subprocess.run([str(RSCRIPT), "-e", 'cat(R.version.string, format(packageVersion("boot")), '
                         'format(packageVersion("cluster")), format(packageVersion("dplyr")), '
                         'format(packageVersion("tibble")), sep = " | ")'], capture_output=True, text=True).stdout
    protocol.write_json(OUT / "provenance" / f"chooser_{ds}_seed{seed}.json",
                        protocol.provenance({"dataset": ds, "seed": seed, "b": b, "designs": designs,
                                             "tally_max": tally_max, "crosscheck": crosscheck, "r_versions": rv,
                                             "chooser_pipeline_sha256": protocol.sha256(CHOOSER),
                                             "minutes": (time.time() - t0) / 60}))
    status = res.status.value_counts().to_dict()
    print(f"  done in {(time.time() - t0) / 60:.1f} min. Status: {status}. "
          + (f"Cross-check: {sum(c['equal'] for c in checks)} of {len(checks)} frequency matrices identical"
             if checks else ""), flush=True)
    return 0


if __name__ == "__main__":
    argv = sys.argv[1:]
    if "--require-tag" in argv:
        argv = [a for a in argv if a != "--require-tag"]
        protocol.require_tag()
    pos, o = args_of(argv)
    sys.exit(main(*pos, designs=o["designs"], workers=o["workers"], tally_max=o["tally-max"],
                  crosscheck=o["crosscheck"]))
