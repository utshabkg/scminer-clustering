# Resume note (for Claude, after a session limit or a restart)

Last updated: 2026-10-14 04:10 CDT. Read this first, then `PROGRESS.md` and `WORKPLAN.md`.

The user asked Claude to keep doing Utshab's tasks until something depends on Laraib or Rhea, to stop before any long GPU run and ask first, and to keep going after a session limit resets.

## Current task (Oct 14): Assignment 2 report, due today

The report is written, verified and pushed. The user compiles it on Overleaf from `latex/`: `main.tex`, `sections/`, `tables/`, `figures/`, `references.bib` and `IEEEtran.cls`. A local compile gives 4 pages of body plus 1 page of references. The preview is `latex/a2_report_preview.pdf`, which is untracked.

- **How it was built.**
  - Tables and the figure come from `scripts/a2_report.py`, and every quoted number is in `results/a2_report_numbers.json`.
  - The hub numbers come from `scripts/rq2_hubs.py` and `results/rq2_hubs/`.
  - The notebook is `notebooks/assignment2_preliminary.ipynb`, built and run by `scripts/make_a2_notebook.py`.
  - The literature was verified by workflow wf_c52c8b84-309. The report itself was verified by wf_682bf22b-e8d (numbers, citations, requirements, writing, reviewer), and its confirmed findings were applied.
- **Left for the user before submitting:**
  1. Put the real Overleaf view link into `\overleaflink` in `latex/main.tex`.
  2. Make the GitHub repo accessible to the professor. It is private, so either make it public or add the professor as a collaborator.
  3. Merge PR #1 (Laraib's Zeisel card) so the card exists on main. Merging was blocked by the permission classifier, so the user must merge on GitHub or allow it.
  4. Merge PR #2 (Rhea's pivots).
  5. Optional: Rhea's Google Sheet is viewable by anyone with the link, while U-30 asks us not to discuss the Kolodziejczyk result outside the team yet.
- **After the merges** (follow-up fixes and PROGRESS/WORKPLAN, held back to avoid conflicts):
  - The PR1 card needs a correct 2,903-cell sentence, post-QC labels on the counts, log2(CPM + 1) and a title line.
  - U-4 goes to [x] and R-1 to [~].
  - Add a PROGRESS.md entry for Oct 14 covering the report, the label-free repeats, the hub script and the corrected hub numbers (44% and 75%, not 42% and 76%).

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
