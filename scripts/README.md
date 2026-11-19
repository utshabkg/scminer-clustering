# Scripts

Run every script from the repository root, for example `python scripts/a3_report.py`. Scripts that take
arguments say how to call them in their first lines.

The scripts the paper (Assignment 4) names, and the protocol v1 code, sit at the top of this folder.
Their paths stay fixed, because the paper cites them and the protocol v1 check compares some of them
with the tag `matched-protocol-v1`. Everything else is grouped below by job.

## Rebuild the paper's tables, figures and numbers (no experiments, about a minute)

| Script | What it does |
|---|---|
| `a4_protocol_v1.py` | Computes every protocol v1 verdict and comparison from `results/v1/`. With `--tables` it also writes the five `latex/tables/a4_*.tex`, Figure 2 and `latex/numbers_a4.tex` |
| `a3_report.py` | Writes most Assignment 3 tables, Figure 1 and `results/a3_report_numbers.json` from the saved runs |
| `a3_matched.py` | Writes the matched-run, transfer and summary tables and `results/a3_matched_numbers.json` |
| `assignment2/a2_report.py` | Writes the reproduction table and `results/a2_report_numbers.json` |
| `log_timings.py` | Copies the timings that exist only in run logs into `results/a3_log_timings.csv` (needs `logs/`) |
| `make_a3_notebook.py` | Builds and runs `notebooks/assignment3_results.ipynb` |

## Protocol v1: the confirmatory runs of Assignment 4

| Script | What it does |
|---|---|
| `run_protocol_v1.sh` | Runs one part per call: the nine small datasets, a PBMC14K stage, or chooseR's code |
| `run_protocol_small.py` | One small dataset and seed: MICA, PCA plus Leiden and PCA plus k-means on 100 subsamples |
| `run_protocol_pbmc.py` | PBMC14K in four restartable stages: full-data sweeps, MICA subsamples, PCA subsamples, scoring |
| `fetch_chooser.sh` | Downloads chooseR's R code at a pinned commit and checks its checksums |
| `run_chooser_r.py`, `score_chooser_r.R` | Score the saved partitions with chooseR's own R code (scorer S3) |
| `check_v1_regression.py` | Checks that seeds 0 to 2 reproduce the Assignment 3 results (criterion R0) |

## The stability rule: the exploratory runs of Assignment 3

| Script | What it does |
|---|---|
| `run_rule_stability.py` | The stability rule on one small dataset, MICA and Leiden |
| `run_rule_stability_graph.py` | The same on PBMC14K, with MICA in graph mode |
| `mica_downstream.py` | MICA's own steps after the neighbour graph (node2vec, then Louvain), used by graph mode |
| `run_rule_mode_size.py` | MICA's two modes and Leiden on one fixed pool of 2,500 PBMC14K cells |
| `run_rule_silhouette.py` | The silhouette rule, scored in one shared space |
| `run_rule_default.py` | MICA's shipped default of four clusters |
| `run_rules_small.sh`, `run_rules_repeats.sh`, `run_rules_pbmc_graph.sh`, `run_rules_silhouette_all.sh` | Batch drivers that ran the scripts above at seeds 0 to 2 |
| `run_rules_matched.sh` | Runs the matched-candidate version of each part, one part and seed per call |
| `probe_grid_coverage.py`, `probe_mica_graph_coverage.py` | Record which cluster counts each method's settings reach, which fixed the matched candidates |

## `large_scale/`: the Allen atlas on every cell

| Script | What it does |
|---|---|
| `run.sh` | Runs one stage, or `all`, with its time and peak memory recorded. Start it in a byobu session |
| `allen.py` | The stages: manifest, three passes over the raw files (gene filter, seurat_v3 genes, normalisation), PCA, k-means, Leiden, assignment of the cells outside MICA's sample, and scoring |
| `validate_100k.py` | Checks the pipeline against the in-memory pilot on the same 100,000 cells |
| `make_tables.py` | Writes the paper's Allen table and `latex/numbers_large.tex` |

## `setup/`: get the data and build the environments

| Script | What it does |
|---|---|
| `setup_envs.sh` | Builds the `scminer` environment (`main`) or the authors' `mica-legacy` environment (`mica`) |
| `prepare_paper_inputs.py` | Converts the scMINER authors' processed inputs for the ten benchmark datasets, values untouched |
| `download_pbmc14k.py`, `download_zeisel.py`, `download_allen.sh` | Download raw data into `data/raw/` |
| `prepare_dataset.py` | Quality control and normalisation of one raw dataset |
| `export_for_mica.py`, `export_hvg_subsets.py` | Write inputs in the format the authors' MICA can read |
| `make_mi_reference.py` | Reference distances from the authors' code, for `tests/test_mi_matches_mica.py` |
| `dataset_facts.py` | Cells, genes and labels per dataset, into `docs/dataset_facts.csv` and `.md` |

## `assignment2/`: reproduction, baselines and the first label-free runs (Oct 1 to 14)

| Script | What it does |
|---|---|
| `a2_report.py`, `make_a2_notebook.py` | Rebuild the Assignment 2 tables, figure and notebook |
| `run_mica_job.sh`, `score_mica.py`, `mica_seeded_kmeans.py`, `score_all_mica.sh` | Run the authors' MICA and score it as the paper did |
| `run_baselines.py` | PCA with k-means, Ward and Leiden, told the true number of clusters |
| `run_label_free.py`, `mica_label_free.py` | The same methods with the number of clusters chosen without labels |
| `bootstrap_gaps.py` | How sure we are of each gap between Leiden and MICA |
| `run_diagnostics.py`, `summarize_diagnostics.py` | What MICA's similarity measures, and read-depth thinning |
| `summarize_gene_sets.py`, `summarize_label_free.py` | Summary tables of those runs |
| `kolod_side_test.py`, `pc_sensitivity.py` | Two side checks, on Kolodziejczyk's batch effect and on the number of components |
| `rerun_logbase_fix.sh`, `run_allen_pilot.sh` | One-off batches: the Oct 3 log-base fix and the Allen 100,000-cell pilot |

## `graph/`: the faster neighbour graph (Appendix of the paper)

| Script | What it does |
|---|---|
| `measure_graph_cost.py`, `measure_mica_hnsw.py` | Time and memory to build our graph and MICA's |
| `run_rq2_graph.py`, `run_rq2_large.py`, `run_rq2_allen.py` (run by `run_rq2_small.sh`, `run_rq2_large.sh` and `run_rq2_allen.sh`) | How well our cheap graph recovers MICA's neighbours |
| `export_fast_graph.py`, `run_downstream_*.sh`, `score_downstream.py` | Run MICA's later steps on either graph and score them |
| `rq2_hubs.py` | Hub cells in MICA's neighbour graph |

## `tools/`: team tools

| Script | What it does |
|---|---|
| `lint_prose.py` | Checks prose for the writing rules we agreed on |
| `overleaf_package.py` | Packs the changed files of `latex/` for a manual upload to Overleaf |
| `make_paper_tables.py`, `export_pivot.py` | Rhea's R-1 and R-2 tables |
| `update_tracker.py` | Rebuilds `docs/results_tracker.csv` from every `results/*/metrics.csv` |
| `peak_mem.sh` | Runs a command and logs its peak memory |

**Do not rerun `update_tracker.py` casually.** `docs/results_tracker.csv` is frozen at Oct 3, and
`a3_report.py` reads it. A rebuild would pull in later runs and change the PBMC14K default value in two
of the paper's tables.

The subfolders were made on Nov 14. Before that, every script sat at the top of this folder, which is
the path older notes in `PROGRESS.md` use. `git log --follow <file>` shows a moved file's full history.
