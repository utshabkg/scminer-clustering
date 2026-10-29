# scMINER clustering: what a published advantage is worth without labels

Course project for Clustering Algorithms (Missouri S&T, Fall 2026). We study scMINER, a mutual
information based framework for clustering single-cell RNA data
([Pan et al., Nature Communications 2025](https://doi.org/10.1038/s41467-025-59620-6)), and its
clustering module MICA.

**The question.** Published benchmarks fix each method's setting with the reference labels, either
by supplying the true number of cell types or by keeping the run that reproduces it. A user has no
labels. So when the same label-free rule has to choose the setting for both methods, does MICA
still beat a standard pipeline of principal components and Leiden?

**What we found.** On the nine small benchmark datasets the answer depends on the rule, not the
method. Reading one published rule (chooseR) two defensible ways reverses the ranking: MICA leads
0.661 to 0.632 under the authors' reading and trails 0.473 to 0.635 under the other, on the same
scores. The cause is measurable. Among settings the rule treats as equivalent, the span of cluster
counts averages 7.3 for MICA against 2.7 for the baseline, and on two datasets the rule cannot
separate 2 clusters from 20.

That failure turns out to belong to one specific pipeline. Running both of MICA's modes on a single
fixed pool of 2,500 cells, holding the data, the size and the subsamples constant, the ambiguity is
present in small-data mode and absent in graph mode. It is not an effect of dataset size, and not a
property of mutual information, since both arms use the same distance on the same cells. It is
MICA's small-data mode, which is the mode nine of the ten benchmark datasets run in.

## Quick start

Rebuild every table and figure in the report from saved results, without rerunning an experiment:

```bash
conda env create -f envs/scminer.yml   # installs the project too, so scminer_clust imports
conda activate scminer
python scripts/a3_report.py            # rewrites latex/tables/*.tex, latex/figures/*, results/a3_*.json
python scripts/a3_matched.py           # the matched-grid, transfer and summary tables
pytest tests/ -q                       # 4 tests, about 6 seconds
```

Run both from the repository root. If `scminer_clust` fails to import, the project was not
installed into the environment: `pip install -e .` fixes it. `envs/*.lock.yml` record the exact
versions our runs used, for reference rather than for building: the torch pin is a CUDA 12.8
local wheel that no index serves.

`notebooks/assignment3_results.ipynb` runs the same thing and displays each table and figure next to
the rows behind it. Both read only from `results/`, which is committed.

## Layout

| Path | What is in it |
|---|---|
| `src/scminer_clust/` | project code: the MI distance (`mi.py`), the clustering pipelines (`cluster.py`), the selection rules (`selection.py`, `stability.py`), metrics |
| `scripts/` | one runner per experiment, each writing a folder under `results/` |
| `results/` | one folder per run, each with `metrics.csv` and a `params.json` recording its settings, seed and candidate grid |
| `notebooks/` | rebuild the report's displays from `results/` |
| `latex/` | the current report, with earlier submissions under `latex/archive/` |
| `docs/` | dataset cards, a primer for non-specialists, literature notes |
| `tests/` | the two properties the experiments rely on |
| `envs/` | conda environments |

## Reproducing the experiments

These rerun the clustering and take minutes to hours. The stability rule is the main one.

```bash
python scripts/run_rule_stability.py <dataset> [n_hvg] [flavor] [seed]   # the nine small datasets
python scripts/run_rule_stability_graph.py [seed] [b_subsamples]         # PBMC14K, MICA graph mode
python scripts/run_rule_mode_size.py [n_pool] [seed] [b_subsamples]      # mode against size, one fixed pool
python scripts/run_rule_silhouette.py <dataset> [n_hvg] [flavor]        # silhouette in one shared space
python scripts/run_rule_default.py                                       # each method's shipped default
python scripts/measure_graph_cost.py <dataset> <c> [workers] [dense]     # shortlist graph: time, memory, recall
python scripts/measure_mica_hnsw.py <dataset> [workers]                  # MICA's HNSW alone (mica-legacy env)
python scripts/probe_grid_coverage.py fine                               # cluster counts a Leiden grid reaches
python scripts/probe_mica_graph_coverage.py [max_resolution]              # the same for MICA's graph mode
```

The three stability runners take `--matched`, which offers both methods the same cluster
counts with one setting per count (`src/scminer_clust/grids.py`) and writes to
`results/rule_matched_*`. `scripts/run_rules_matched.sh small|pbmc|modesize <seed>` runs one
part at one seed.

`scripts/run_rules_small.sh` runs the nine small datasets at seed 0 and
`scripts/run_rules_repeats.sh` supplies seeds 1 and 2. MICA's graph mode needs the authors' own
environment, because the node2vec step runs there as a subprocess:

```bash
bash scripts/setup_envs.sh mica   # clones MICA 1.0.0 and builds the mica-legacy env
```

`envs/mica-legacy.lock.yml` records that environment rather than building it, since MICA is
installed from source.

Two things make the stability rule affordable. MICA's distance between two cells depends only on
those two cells, so the full matrix is computed once per dataset and subsampled, which is exact
rather than an approximation and is checked by `tests/test_subsample_invariance.py`. And our
implementation of the distance agrees with the authors' code to within `1e-6`
(`tests/test_mi_matches_mica.py`).

## Data

Datasets are not in this repository. They live on the lab workstation's shared disk at
`/media/20TB/shared/sides_data/43.scMINER-Clustering/`, which the gitignored `data` symlink points
to. The ten benchmark datasets are used in the processed form the scMINER authors released, so any
difference from their published numbers cannot come from our preprocessing. `scripts/` holds the
download and preprocessing steps, and `docs/dataset_facts.csv` records the cell and gene counts
after filtering for every dataset.

Long runs go in a detached byobu session with a log:

```bash
byobu new-session -d -s scm-<task> "bash -lc 'source ~/anaconda3/etc/profile.d/conda.sh && conda activate scminer && <command> 2>&1 | tee logs/<task>.log'"
```

## Where to look first

- [PROGRESS.md](PROGRESS.md): what has actually been done, newest first, with the numbers.
- [WORKPLAN.md](WORKPLAN.md): who does what, by when, and how a task counts as finished.
- [docs/primer.md](docs/primer.md): the data, the methods and the question, assuming no background.
- `notebooks/assignment3_results.ipynb`: every current result, rebuilt from saved files.

## Team

Utshab Kumar Ghosh (lead, all code and experiments), Laraib Khan, Rhea Eckman.
