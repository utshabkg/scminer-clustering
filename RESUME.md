# Resume note (for Claude, after a session limit or a restart)

Last updated: 2026-10-03 17:50 CDT. Read this first, then `PROGRESS.md` and `WORKPLAN.md`.

The user asked Claude to keep doing Utshab's tasks until something depends on Laraib or Rhea, to stop before any long GPU run and ask first, and to keep going after a session limit resets.

## Running right now (Oct 3, 17:50)

MICA's own node2vec and Louvain, run on the Allen 100K graphs. This is the last step of U-19 at this size.
- **`scm-downstream-allen`** runs `scripts/run_downstream_allen.sh`, log `logs/downstream_allen.log`. Done: our corr3200 graph scored ARI 0.602 and AMI 0.819. MICA's own pilot got 0.652 and 0.836. Running now: the same wrapper on MICA's own graph, about 70 min.
- **`scm-downstream-allen2`** waits for that to finish, then runs a second repeat on both graphs, about 2.3 h. Log `logs/downstream_allen_rep2.log`, ends with DOWNSTREAM ALLEN REP2 DONE.

**How to read the result.**
- If MICA's graph gives about 0.65 through our wrapper in both repeats while ours stays near 0.60, the graph costs ARI at this size, even though its recall is higher (94.6% against 90.3%). That would fit the hub and Klein findings.
- If both graphs vary by about 0.05, it is run-to-run noise.
- Write whichever holds into the Allen paragraph of PROGRESS.md. It currently ends with "MICA's own later steps on our graph are running now". Then commit, push, and report to Utshab.

The rest of U-19 is written in PROGRESS.md and pushed in 4a8811c. That covers the Allen table: HNSW 90.3% in 5.2 h, correlation 1,600 at 89.8% in 1.1 h, and correlation 3,200 at 94.6% in 2.2 h.

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
