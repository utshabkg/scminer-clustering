# Resume note (for Claude, after a session limit or a restart)

Last updated: 2026-10-03 13:40 CDT. Read this first, then `PROGRESS.md` and `WORKPLAN.md`.

The user asked Claude to keep doing Utshab's tasks until something depends on Laraib or Rhea, to stop before any long GPU run and ask first, and to keep going after a session limit resets.

## Running right now (Oct 3, 13:40)

**`scm-rq2-allen`** runs `scripts/run_rq2_allen.sh`, which is v2 with the sparse exact path. It takes about 2.5 h. The log is `logs/rq2_allen.log` and ends with RQ2 ALLEN DONE. Output goes to `results/rq2_allen_100k/`, and the c=3200 graph is written in MICA format to `data/runs/rq2/allen_100k_corr3200/`. Brute force for the 500 queries is cached in `data/runs/rq2/allen_100k_exact_q500.npz`.

Everything else is finished:
- the small sets, Zeisel included
- PBMC14K v3
- the PBMC14K downstream runs

The U-19 section in the PROGRESS.md Oct 3 entry is written and pushed in 3df88c4. It covers the PBMC14K table, MICA downstream 0.853 against 0.852, the Allen query sweep, the failures, the hubs and the caveats.

**When Allen finishes:**
1. Replace "The full timed build at that size is running now." in PROGRESS.md with the measured times and the Leiden numbers.
2. Run MICA's own downstream on the Allen graph. In mica-legacy: `python scripts/mica_downstream.py data/runs/rq2/allen_100k_corr3200 data/runs/rq2/downstream_allen_100k_corr3200_rep1 8`. It takes about 80 min, about 40 for node2vec and 40 for Louvain. Then score it in scminer: `python scripts/score_downstream.py allen_100k <run_dir> corr3200 rq2_downstream_allen_100k_corr3200_rep1`. Compare with MICA's 0.652, which is in `results/pilot_mica_allen_100k`.
3. Commit and push. Then tell Utshab the results and ask about two things:
   - GPU time for U-18
   - whether to add a hubness-reduction task

**Measured on Allen so far:**
- Dense brute force took 1,832 s per 500 cells and sparse took 767 s.
- MICA HNSW recall is 0.903 in 5.2 h.
- Correlation recall: 200 gives 0.558, 400 gives 0.710, 800 gives 0.823, 1600 gives 0.898 and 3200 gives 0.946.
- Codes are 17.6% nonzero.

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
