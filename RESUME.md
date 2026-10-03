# Resume note (for Claude, after a session limit or a restart)

Last updated: 2026-10-03 00:30 CDT. Read this first, then `PROGRESS.md` and `WORKPLAN.md`.

The user asked Claude to keep doing Utshab's tasks until something depends on Laraib or Rhea, to stop before any long GPU run and ask first, and to keep going after a session limit resets.

## Running right now (byobu sessions on the workstation, they survive a session limit)

| Session | Task | What it does | When it ends |
|---|---|---|---|
| `scm-allen-pilot` | U-16 | `scripts/run_allen_pilot.sh`: prepare 100,000 atlas cells (done, 7 min, 24 GB), baselines on all genes (done, 25 min, 55 GB), baselines on 2,000 seurat_v3 genes (running), then MICA graph mode on the 100,000 cells | Several hours. Timing and peak memory go to `logs/allen_pilot_*.log` (lines PEAK_RSS_GB and WALL_SECONDS) |

Check with `byobu list-sessions | grep scm-`. Never edit a bash script that a running job is executing.

## When the pilot finishes

```bash
source ~/anaconda3/etc/profile.d/conda.sh && conda activate scminer
export OMP_NUM_THREADS=8
python scripts/score_mica.py allen_100k results/reproduction/mica_allen_100k --exp pilot_mica_allen_100k
python scripts/mica_label_free.py allen_100k results/reproduction/mica_allen_100k --exp labelfree_mica_allen_100k
python scripts/update_tracker.py
```

Then write the timing and memory of every step into PROGRESS.md and extrapolate to 2.3 million cells (MICA's cost grows faster than linearly).

## Unfinished in-session work

1. **Verification of the Oct 3 claims.** A workflow (run `wf_217dd7e1-06d`) checks `docs/research/verification/claims_oct3.md` with three agents: numbers, skeptic, code review. If it was lost, rerun it with `docs/research/verification/verify_oct3_workflow.js`. Fix whatever it finds before writing anything up.
2. **Write-up, not done yet.**
   - Add an Oct 3 entry to PROGRESS.md.
   - Put a note on the Oct 2 entry: MICA small-data numbers changed on Oct 3 with 100 k-means restarts.
   - Update the Where we stand table.
   - Update WORKPLAN statuses. Done: U-12, U-13, U-14, U-15, U-17, U-27, U-28. In progress: U-16.
   - Make a simple results file for Rhea's R-1, because the tracker now has 1,410 rows with label-free and authors'-input rows mixed in.
   - Commit and push.
3. **Next CPU tasks** after the write-up:
   - U-29: the Kolodziejczyk side test.
   - A CPU prototype of the retrieve-then-rerank neighbor graph (U-19) on small data.
4. **Stop and ask Utshab** before U-18, the GPU kernel, and anything else that uses the GPUs. Also ask about U-30, contacting the authors.

## What we know so far (details in `results/summary_*.md`)

- **Reproduction.** Our MICA distance matches the authors' code to 1e-6. On the authors' inputs, MICA reproduces 9 of 10 published scores once its final k-means reaches the best fit (100 restarts). Kolodziejczyk does not reproduce (0.727 against 1.0).
- **Told the number of clusters.** MICA beats Leiden on the five small datasets in every gene setting. On the five mid-sized ones, the gene set decides the winner.
- **Choosing the number itself** (silhouette). MICA's average lead over Leiden falls from 0.13 to 0.02. It loses on small datasets and wins on mid-sized ones.
- **What MICA measures.** On Zeisel and Kolodziejczyk, MICA's similarity is mostly co-detection. On Usoskin and Klein it is not. MICA is robust to sequencing depth except on Usoskin.

## Environment

- Envs: `scminer` (our code) and `mica-legacy` (authors' MICA).
- Data: `data/` links to `/media/20TB/shared/sides_data/43.scMINER-Clustering`.
- Cap threads on every job (`OMP_NUM_THREADS` and similar, 4 to 16). The machine is shared, and uncapped jobs pushed the load to 392.
- MICA small-data scoring: `scripts/score_all_mica.sh` does seeded draws, scoring and label-free in one go.
