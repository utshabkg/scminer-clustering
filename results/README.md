# Results

The saved output of every run, committed so that every table and number in the paper can be rebuilt
without rerunning an experiment (see the main `README.md`). Most folders hold one run: a `metrics.csv`
with one row per method and setting, and a `params.json` with the settings, seed and candidates used.

## How to read a folder name

- `paper_` means the scMINER authors' processed input. In the Assignment 2 folders a dataset name
  without it, such as `baselines_klein`, means our own preprocessing of the raw data. From Assignment 3
  on, every run on a benchmark dataset uses the authors' input, whatever the name says.
- `_seed1`, `_seed2` and so on give the seed when each seed has its own folder. In the
  `rule_stability_`, `rule_matched_stability_` and silhouette folders, a name with no seed is seed 0.
  Most other folders keep seeds 0 to 2 as rows of one `metrics.csv`, in its `seed` column.
- `_rep1` to `_rep5` are repeats of the same MICA job, which takes no seed.
- `_hvg2000_<flavor>` means only the 2,000 most variable genes, picked with that Scanpy flavor
  (`seurat`, `seurat_v3` or `cell_ranger`). A name ending in plain `_hvg2000` means `seurat`.
- `rule_matched_` means both methods were offered the same cluster counts. `rule_` alone means each
  method's own grid of settings.

## Assignment 4: protocol v1

| Folder | What it holds | Made by |
|---|---|---|
| `v1/` | Every protocol v1 run: seeds 0 to 9 on the nine small datasets (3 to 9 are the confirmatory seeds, 0 to 2 the exploratory ones rerun) and 0 to 4 on PBMC14K. `picks/` holds each method's choice per design, scorer and reading, `curves/` every setting's stability score and interval, `chooser_r/` and `picks_chooser/` the results of chooseR's own code, `provenance/` the code version and tag check of every run, and `regression.csv` criterion R0. Its `rule_*` folders are the 25-subsample scoring in the Assignment 3 layout, under the same names as the Assignment 3 folders one level up | `scripts/run_protocol_*.py`, `scripts/run_chooser_r.py`, `scripts/check_v1_regression.py` (`regression.csv`) |
| `v1_pretag/` | Smoke tests run before the protocol was tagged. No number in the paper comes from here | the same scripts |
| `a4_protocol_v1_numbers.json` | Every criterion's verdict and every value the paper prints through `latex/numbers_a4.tex`. The few numbers typed by hand in the text are not in it | `scripts/a4_protocol_v1.py` |

## Assignment 3: the stability rule at seeds 0 to 2

| Folders | What they hold | Made by |
|---|---|---|
| `rule_stability_paper_*` | The stability rule on each method's own grid, nine small datasets and PBMC14K | `scripts/run_rule_stability.py`, `scripts/run_rule_stability_graph.py` |
| `rule_matched_stability_paper_*` | The same on matched candidates | the same, with `--matched` |
| `rule_silhouette_paper_*`, `rule_matched_silhouette_paper_*` | The silhouette rule in one shared space | `scripts/run_rule_silhouette.py` |
| `rule_mode_size_*`, `rule_matched_mode_size_*` | MICA's two modes on one pool of 2,500 PBMC14K cells, three pools | `scripts/run_rule_mode_size.py` |
| `rule_default_mica_small_data` | MICA's shipped default of four clusters | `scripts/run_rule_default.py` |
| `grid_coverage_*` | Which cluster counts each method's settings reach, used to fix the matched candidates | `scripts/probe_*_coverage.py` |
| `a3_*.json`, `a3_*.csv` | The numbers and per-dataset tables the Assignment 3 report scripts collect | `scripts/a3_report.py`, `scripts/a3_matched.py`, `scripts/log_timings.py` |

## Assignment 2: reproduction, baselines and first label-free runs

| Folders | What they hold | Made by |
|---|---|---|
| `reproduction_*`, `repeat_paper_*` | The authors' MICA output, scored as the paper did (true number of clusters) | `scripts/assignment2/score_mica.py` |
| `baselines_*` | PCA with k-means, Ward and Leiden, told the true number of clusters | `scripts/assignment2/run_baselines.py` |
| `labelfree_*`, `labelfree_mica_*`, `repeat_labelfree_*` | The same methods choosing the number of clusters without labels | `scripts/assignment2/run_label_free.py`, `scripts/assignment2/mica_label_free.py` |
| `mica_*_hvg2000_*` | MICA on 2,000 variable genes | `scripts/assignment2/score_mica.py` |
| `rq1_diagnostics_*`, `kolod_side_test`, `pc_sensitivity`, `uncertainty`, `pilot_mica_allen_100k` | Side checks and the Allen pilot | `scripts/assignment2/` |
| `summary_*` | Summary tables of these runs | `scripts/assignment2/summarize_*.py` |
| `a2_report_numbers.json` | Every number the reproduction table and the Assignment 2 report quote | `scripts/assignment2/a2_report.py` |

## The Allen atlas at full scale

| Folder | What it holds | Made by |
|---|---|---|
| `large_scale/` | PCA with k-means and with Leiden on all 2,341,350 cells, MICA on the 100,000-cell sample and with the other cells assigned: scores (`metrics.csv`), time and memory per stage (`stage_costs.csv`), the Leiden resolution search, seed agreement and the holdout check of the assignment | `scripts/large_scale/allen.py`, `make_tables.py` |
| `large_scale_sample100k/` | The same pipeline on the pilot's 100,000 cells, and `validation.json`, its comparison with the in-memory pilot | `scripts/large_scale/validate_100k.py` |

## The neighbour graph comparison

| Folders | What they hold | Made by |
|---|---|---|
| `graph_cost_*` | Time and memory to build our graph and MICA's | `scripts/graph/measure_graph_cost.py`, `scripts/graph/measure_mica_hnsw.py` |
| `rq2_*` | How well our graph recovers MICA's neighbours, and MICA's later steps run on each graph | `scripts/graph/` |

Raw MICA output (21 GB) lives in `results/reproduction/` on the lab machine only, outside git. The
scored copies above are what the paper uses.
