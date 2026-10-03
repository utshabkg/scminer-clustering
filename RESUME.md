# Resume note (for Claude, after a session limit or a restart)

Last updated: 2026-10-03 12:40 CDT. Read this first, then `PROGRESS.md` and `WORKPLAN.md`.

The user asked Claude to keep doing Utshab's tasks until something depends on Laraib or Rhea, to stop before any long GPU run and ask first, and to keep going after a session limit resets.

## Running right now (Oct 3, 12:35)

Three CPU jobs, all for U-19:
- **`scm-rq2`** runs `scripts/run_rq2_small.sh`. Usoskin, Kolodziejczyk and Klein are done, Zeisel is running. Log: `logs/rq2_small.log`, ends with RQ2 SMALL DONE.
- **`scm-rq2-v3`** runs `scripts/run_rq2_large.py paper_pbmc14k`. Log: `logs/rq2_large_v3.log`, ends with RQ2 V3 DONE. Output: `results/rq2_large_paper_pbmc14k/`. It overwrites the stopped v2's results. v2's numbers are kept in `logs/rq2_large_v2_partial.log`.
- **`scm-rq2-allen`** runs `scripts/run_rq2_allen.sh`, about 2 to 4 h. Log: `logs/rq2_allen.log`, ends with RQ2 ALLEN DONE. Output: `results/rq2_allen_100k/`.

**U-19 findings so far.** None of this is written into PROGRESS.md yet.
- **PCA candidates fail at scale.** On PBMC14K with 50 PCs, 400 candidates find 33% of the true 79 neighbors. Detected genes find 72%, and both merged find 76%. Neighbors-of-neighbors rounds level off near 86 to 88%.
- **Correlation on all genes works.** No PCA, rescored exactly. On PBMC14K it reaches 90.4% at 400 candidates, 97.4% at 800 and 99.6% at 1,600. MICA's HNSW gets 97.3% (taking each cell's 79 closest links in MICA's graph, which flatters MICA).
  - Timing for all cells: correlation itself takes 5 s, and rescoring 200 candidates per cell takes 23 s with 8 workers. HNSW took 903 s.
  - Brute force would take about 35 min for all 13,605 cells.
- **Correlation on the small sets:** Kolodziejczyk 99.5% at 200 candidates, Usoskin 97% at 400, Klein 92% at 800. Klein is the hardest.
- **Hub injection does not help.** Adding the graph's top hub cells to every list is worse than adding more correlation candidates. Hubs account for only 15% of the misses.
- **Hubs are real in the exact MI graph:**
  - PBMC14K: the top 1% of hub cells take 23% of all exact neighbor slots. Their links join the same label 42% of the time, against 76% for other links.
  - Klein: in-degree has a median of 7 and a maximum of 868.
  - The earlier claim that "entropy drives the ranking" was wrong. Across 50 queries the median Spearman is 0.09. Don't repeat it.
- **Leiden on the exact MI graph is a weak stand-in for MICA.** On Klein it scores 0.59, while the approximate graphs score 0.84 to 0.98. On PBMC14K it scores 0.64 on MICA's graph and 0.61 on the fast graph, against MICA's own 0.853. So graph recall may not be the right goal by itself. Hubness reduction is a possible new task: propose it to Utshab and don't start it unasked.
- **Exploration scripts** are in the scratchpad, not the repo. They show up in the logs: `rq2_proxies.log` and `rq2_hubs.log`.

**When the runs finish.**
1. Write a U-19 section into PROGRESS.md. Mainly cover correlation candidates, the 97% recall against HNSW, the time, the Allen numbers, the hubs, and the weak Leiden stand-in.
2. Commit and push.

## Done on Oct 3 (all pushed except as noted)

- The log-base fix rerun finished, and every quoted number was rechecked after it. The corrected wording is in the PROGRESS.md Oct 3 entry, so use that and not older notes. Key changes:
  - MICA leads Leiden in 18 of 20 small settings, not 19. Chung seurat is now a tie.
  - MICA beats all baselines on Chung in 3 of 4 settings, not 4.
  - The 2k Wilcoxon p ranges from 0.08 to 0.43.
  - The label-free lead falls from 0.12 to 0.03.
  - With fewer PCs, Buettner Leiden scores 0.51 to 0.78 (`results/pc_sensitivity/`).
- WORKPLAN statuses are updated: U-12 to U-17 and U-27 to U-29 are done, U-19 is in progress. Four decisions are added to the decision log.
- `docs/allen_pilot.md` holds the pilot timings and the 2.3M extrapolation. `docs/allen_files.csv` lists the raw files.
- `scripts/export_pivot.py` writes `docs/results_for_pivot.csv` for Rhea's R-1. Rerun it after `update_tracker.py`.

## Next tasks

1. U-19 on larger data. PBMC14K (13,605 cells) is the first place a saving should show. The exact all-pairs matrix there is 93M pairs, which numpy can't do in reasonable time. So measure recall on a random sample of about 500 query cells with exact search against all cells, using `MICodes.dist` row by row. Don't build the full matrix.
2. **Stop and ask Utshab** before U-18, the GPU kernel, or anything else that uses the GPUs. Also ask about U-30, contacting the authors.

## What we know (details in PROGRESS.md Oct 3)

- **Reproduction:** 9 of 10 published MICA scores are reproduced within 0.05 on the authors' inputs. Kolodziejczyk is the exception (0.727 against 1.0). Plate centering does not explain it.
- **Matched genes:** MICA wins clearly on the five tiny sets. On the five larger sets the gene set decides. Across datasets the result is not significant after correction.
- **Label-free:** inconclusive. The result flips with the k range and the embedding.
- **Similarity:** co-detection carries 57 to 62% of the MI on Zeisel, Kolodziejczyk and Usoskin, and 25% on Klein. MICA is depth-robust except on Usoskin.
- **Allen:** MICA as shipped can't do 2.3M cells on this machine. Memory is about 3.2 TB extrapolated.
- **U-19 so far:** merged PCA and detected-gene candidates reach 97.5 to 99.3% recall at k=80. The landmark proof doesn't work.

## Environment

- Envs: `scminer` (our code) and `mica-legacy` (authors' MICA).
- Data: `data/` links to `/media/20TB/shared/sides_data/43.scMINER-Clustering`.
- Cap threads on every job (`OMP_NUM_THREADS` and similar, 4 to 16). The machine is shared, and uncapped jobs pushed the load to 392.
- Never edit a bash script while a job is running it.
- `gh` is not installed. Push with plain git, which has the credentials.
- MICA small-data scoring: `scripts/score_all_mica.sh` does the seeded draws, scoring and label-free runs in one go.
