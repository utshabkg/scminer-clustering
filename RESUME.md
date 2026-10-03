# Resume note (for Claude, after a session limit or a restart)

Last updated: 2026-10-03 11:30 CDT. Read this first, then `PROGRESS.md` and `WORKPLAN.md`.

The user asked Claude to keep doing Utshab's tasks until something depends on Laraib or Rhea, to stop before any long GPU run and ask first, and to keep going after a session limit resets.

## Running right now

`scm-logbase-fix` runs `scripts/rerun_logbase_fix.sh` (log in `logs/rerun_logbase_fix.log`, ends with RERUN DONE, about an hour).

**Why.** The Oct 3 review found the authors' inputs are not all log2: Pollen and Buettner are log10, Chung and Usoskin are natural log (`PAPER_LOG_BASE` in `src/scminer_clust/data.py`). So their seurat and seurat_v3 gene sets were wrong.

**What it does.** Redoes, for those four datasets and those two methods:
- the gene exports
- MICA
- the baselines
- the label-free runs

Then it rescores all MICA runs and reruns the bootstrap and every summary. All-gene and cell_ranger runs were not affected.

Earlier work, already done: The Allen pilot (U-16) finished overnight and is scored:
- MICA graph mode on 100,000 cells: ARI 0.652 at class level, 7.1 h (5.3 h of it the MI neighbor graph), 140 GB peak.
- Leiden: 0.661 with 2,000 genes and 0.624 with all genes, about 17 min each, 27 and 55 GB peak.
- Prepare: 7 min, 24 GB.

U-29 (Kolodziejczyk side test) is done, in `results/kolod_side_test/`. Per-plate centering does not recover 1.0, it drops MICA to 0.16.

## Unfinished in-session work

1. **Verification of the Oct 3 claims is done.** Findings are saved in `docs/research/verification/oct3_findings.txt`. Bugs fixed: the log base, and graph-mode repeats leaking into the summaries.

   Wording to use in the write-up:
   - MICA leads Leiden in 19 of 20 small-dataset settings (Yan seurat is the exception).
   - 100 restarts give a near-best fit, not always the best. On Pollen the best fit scores 0.980.
   - Usoskin near-optimal fits span ARI 0.78 to 0.87.
   - All genes: MICA wins 7, ties Klein, loses Zeisel and PBMC14K. p = 0.027 does not survive correction for four gene settings.
   - Our Leiden is within 0.035 of the paper's Scanpy on 7 of 10.
   - The label-free comparison flips with the k range and the embedding used, so it is inconclusive.
   - Leiden on tiny data is sensitive to the number of PCs (Buettner goes from 0.29 to 0.81 with 10 PCs).
   - The bootstrap ignores rerun variance and cell dependence, so its intervals are too narrow.

   Recheck every quoted number after the rerun.
2. **Write-up, not done yet.**
   - Add an Oct 3 entry to PROGRESS.md.
   - Put a note on the Oct 2 entry: MICA small-data numbers changed on Oct 3 with 100 k-means restarts.
   - Update the Where we stand table.
   - Update WORKPLAN statuses. Done: U-12, U-13, U-14, U-15, U-16, U-17, U-27, U-28, U-29.
   - Make a simple results file for Rhea's R-1, because the tracker now has 1,410 rows with label-free and authors'-input rows mixed in.
   - Commit and push.
3. **Next CPU tasks** after the write-up:
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
