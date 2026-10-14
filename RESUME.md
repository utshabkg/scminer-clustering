# Resume note (for Claude, after a session limit or a restart)

Last updated: 2026-10-14 05:30 CDT. Read this first, then `PROGRESS.md` and `WORKPLAN.md`.

The user asked Claude to keep doing Utshab's tasks until something depends on Laraib or Rhea, to stop before any long GPU run and ask first, and to keep going after a session limit resets.

## Where things stand (Oct 14, after the Assignment 2 submission)

Nothing is running. Everything is committed and pushed.

**Assignment 2 is done.** The report is in `latex/` (`main.tex` plus `sections/00_abstract` to `05_schedule`), 4 pages of body plus 1 page of references. Utshab compiles it on Overleaf. The local preview `latex/a2_report_preview.pdf` is untracked.
- Tables and the figure come from `scripts/a2_report.py`, and every quoted number is in `results/a2_report_numbers.json`.
- Hub numbers come from `scripts/rq2_hubs.py` and `results/rq2_hubs/`.
- The notebook is `notebooks/assignment2_preliminary.ipynb`, built by `scripts/make_a2_notebook.py`.
- Two verification workflows checked it: wf_c52c8b84-309 for the literature, wf_682bf22b-e8d for numbers, citations, requirements, writing and a reviewer read-through. Confirmed findings were applied.
- To rebuild the PDF locally: copy `latex/` to a scratch folder, then `conda activate texenv && tectonic -X compile main.tex`. PyMuPDF in that env counts pages and renders previews.

**Both teammate pull requests are merged.** PR #1 (Laraib's Zeisel card, L-1) and PR #2 (Rhea's pivots, R-1). The card's wrong 2,903-cell claim is fixed. R-1 stays `[~]` until the odd results are marked.

**One thing still open for Utshab:** the Overleaf link. `latex/main.tex` line 15 still has the placeholder `REPLACE-WITH-OUR-LINK`. Utshab pasted the real link into Overleaf directly, so if he ever copies the local `main.tex` over again, the placeholder would overwrite it. Ask whether to store the real link in the repo.

## Next tasks (no GPU, nothing blocked)

The plan was realigned to the instructor's targets: main experiments by Nov 4, full draft by Nov 18. Eight tasks were added, U-31 to U-38.

1. **U-32, the stability rule.** Pick the setting whose clusters best survive subsampling, for MICA and Leiden alike, in the spirit of MultiK and chooseR. It varies k in MICA's small-data mode. Start on PBMC14K and Zeisel. This is the rule that gives the final clustering we report.
2. **U-33, the significance rule.** Merge clusters a test cannot separate, in the style of sc-SHC.
3. **U-35 and U-36.** Rerun the silhouette rule in one shared embedding over a matched range of cluster counts, and run label-free Leiden on the Allen subset so it becomes a third anchor.
4. **U-34 and U-37.** The transfer rule (hold out one dataset, take the best resolution from the others) and the final comparison table across all eleven sets.

**Still needs Utshab's decision:** GPU time for U-18 (now gated behind a Nov 4 result), and U-30, whether to write to the scMINER authors about the Kolodziejczyk input.

**Rules fixed in the report, do not quietly change them:** every rule is applied the same way to both methods, its settings are committed before scoring, and a win for MICA means a higher median ARI on at least 7 of the 11 sets and a mean lead of at least 0.05.

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
