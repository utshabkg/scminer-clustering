# Resume note (for Claude, after a session limit or a restart)

Last updated: 2026-10-22 12:10 CDT. Read this first, then `PROGRESS.md` and `WORKPLAN.md`.

The user asked Claude to keep doing Utshab's tasks until something depends on Laraib or Rhea, to stop before any long GPU run and ask first, and to keep going after a session limit resets.

## Where things stand (Oct 22, Assignment 3 in progress)

**Assignment 3 is due Wednesday Oct 28.** 5 to 7 pages excluding references and appendix,
in the same Overleaf project. It needs new results, a revised introduction, related work and
methods, and a completion plan naming each person. The instructor's writing guide is
`docs/Project_Paper_Writing_Guide.docx`, and `scripts/lint_prose.py` enforces the team's
style rules across all prose.

**The manuscript is drafted** in `latex/`, restructured to the guide's order: abstract,
changes since Assignment 2, introduction, related work, problem and method, experimental
design, results, discussion, conclusion with the completion plan, then a technical appendix.
It compiles to 7 pages of body plus references and appendix. Assignment 2 is archived in
`latex/archive/assignment2`.

**The new experiment** is the chooseR-style subsampling stability rule applied identically to
MICA and to PCA plus Leiden, in `src/scminer_clust/stability.py` and
`scripts/run_rule_stability.py`. Seed 0 is complete on all nine small-data datasets.

**Headline result.** The ranking depends on how the rule breaks ties, not on the data.
- Read as chooseR specifies (finest setting above a threshold): MICA 0.658, Leiden 0.624, MICA wins 5 of 9.
- Read as the plain maximum: MICA 0.473, Leiden 0.613, MICA wins 3 of 9.
- MICA's mean moves 0.185 between the two readings, Leiden's 0.011.
- Mechanism, measured: among settings tied within 0.01 of the best, the span of cluster
  counts averages 6.7 for MICA and 2.3 for Leiden. On Yan the rule cannot separate 3
  clusters from 19.
- A hypothesis that did NOT hold: stability is not less informative for MICA. It correlates
  slightly better with ARI (0.68 against 0.53). The problem is the compressed range.
- Adding Zeisel weakened the spread-versus-overclustering correlation from -0.74 to -0.45,
  so that prediction is reported as not confirmed. Do not restore the stronger wording.

**Running now.**
- `scm-repeats`: `scripts/run_rules_repeats.sh`. MICA's shipped default (small-data mode
  defaults to k=4 whatever the dataset), the shared-space silhouette rule, then stability
  seeds 1 and 2. Log `logs/rules_repeats.log`, ends with RULES REPEATS DONE.
- A background smoke test of `scripts/run_rule_stability_graph.py 0 2`, which is PBMC14K in
  MICA's graph mode. It crosses into the mica-legacy environment for node2vec. If it works,
  run it with B=25 for the flagship dataset, about 3 hours.

**Before submitting, check:** the abstract claims three seeds, so confirm seeds 1 and 2
finished and regenerate with `scripts/a3_report.py`, which rewrites every table and
`results/a3_report_numbers.json`. Then run `scripts/lint_prose.py` and recompile.

**Teammates.** Laraib asked for preprocessing numbers weekly and for technical background:
`scripts/dataset_facts.py` writes `docs/dataset_facts.md` from the processed files, and
`docs/primer.md` explains the project for a non-CS reader. Rhea's PR #3 is merged, and the
six anomalies she raised were each checked and are genuine results, written up in the
manuscript appendix.

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
