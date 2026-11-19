# When Clustering Methods Must Choose for Themselves

**Label-free selection and single-cell clustering reliability.** A course project for Clustering
Algorithms at Missouri S&T, Fall 2026, by Utshab Kumar Ghosh, Laraib Khan and Rhea Eckman.

Benchmarks of single-cell clustering usually pick each method's setting with the reference labels. A
real user has no labels. We let one label-free rule pick the setting for every method and ask how that
changes the ranking, and when the ranking can be trusted. The case study is MICA, the mutual information
clustering in scMINER ([Pan et al., 2025](https://doi.org/10.1038/s41467-025-59620-6)), against principal
components followed by Leiden or by k-means, on the ten scMINER benchmark datasets. The rule is
subsampling stability in the style of chooseR.

## What we found

These are the Assignment 4 results: a protocol fixed in advance, seven new seeds on the nine small
datasets and five seeds on PBMC14K. The numbers are in `latex/numbers_a4.tex`.

- **The candidates decide the ranking.** On each method's own grid of settings, Leiden beats MICA on
  the nine small datasets at every new seed, 0.650 against 0.452 in mean ARI, when the rule takes the
  most stable setting. Offered the same cluster counts (matched candidates), the two are level, 0.454
  against 0.452.
- **So does the reading of the scores.** Taking the finest stable setting instead, as chooseR does,
  MICA leads on matched candidates at all seven seeds, 0.647 against 0.591. But the datasets split about
  evenly, and with 25 subsamples the lead rests on one dataset.
- **chooseR's own code often cannot choose.** When a coarse setting reproduces perfectly, its
  bootstrap interval fails and the code makes no choice. That happens in 106 of 190 MICA runs, counting
  every seed and dataset and both candidate designs.
- **PBMC14K**, the one larger dataset, gives a tie on matched candidates at four of five seeds, but the
  tie changes with which setting stands for each cluster count.
- PCA plus k-means scores lowest in every design on the nine small datasets.
- **At 2.3 million cells** (the Allen mouse brain atlas) both PCA pipelines run on every cell in about 20
  minutes. MICA would need about 3.2 TB of memory, so it runs on a 100,000-cell sample and the other
  cells are assigned to its clusters. The paper keeps that approximation apart from the full-data results.

The paper turns this into two kinds of checks. Six need no labels, so anyone can run them on new data:
they ask whether a rule's choice is determined by the data. Five need labels: they ask whether a
benchmark's ranking of methods holds up.

## Start here

| To | Look at |
|---|---|
| Read the paper | [`latex/`](latex/), starting from `latex/main.tex`. The PDF is built on Overleaf. `latex/README.md` says which file and script make each table |
| Trace a number in the paper | Most protocol v1 numbers come from `latex/numbers_a4.tex`, written by `scripts/a4_protocol_v1.py` from the run files in [`results/v1/`](results/v1/). The few still typed by hand are listed in `docs/r7_number_check_a4.md` |
| Read the protocol we fixed before the confirmatory runs | [`docs/protocol_matched_v1.md`](docs/protocol_matched_v1.md), tagged `matched-protocol-v1` |
| Learn the background, with no prior knowledge | [`docs/primer.md`](docs/primer.md) |
| See what was done, when and by whom | [`PROGRESS.md`](PROGRESS.md) (newest first) and [`WORKPLAN.md`](WORKPLAN.md) |

## Rebuild every table, figure and number

This takes about a minute and reruns no experiment. The scripts read only the committed files in
`results/` and `docs/`.

```bash
conda env create -f envs/scminer.yml
conda activate scminer
pip install -e .                          # this project, so scminer_clust imports
python scripts/a4_protocol_v1.py --tables # protocol v1: Tables I to V, Figure 2, latex/numbers_a4.tex
python scripts/a3_report.py               # Assignment 3 tables and Figure 1
python scripts/a3_matched.py              # the transfer table and the matched-run numbers
python scripts/assignment2/a2_report.py   # the reproduction table
pytest -q tests                           # 19 tests, about 10 seconds
git status                                # nothing should have changed
```

Run everything from the repository root. Tests that need the datasets skip themselves when the data are
missing. The figures come out byte for byte the same with the matplotlib version in
`envs/scminer.lock.yml`, and look the same with others.

## What is where

| Folder | Contents |
|---|---|
| [`latex/`](latex/) | The paper. The scripts write its tables, figures and most of its numbers. Earlier submissions are in `latex/archive/` and the git tag `a3-submitted` |
| [`scripts/`](scripts/) | The scripts the paper names and the protocol v1 code sit at the top. Data setup, Assignment 2, the graph comparison and team tools are in subfolders. See `scripts/README.md` |
| [`results/`](results/) | The saved output of every run, with protocol v1 in `results/v1/`. See `results/README.md` |
| [`src/scminer_clust/`](src/scminer_clust/) | Shared code: the mutual information distance, the clustering pipelines, stability scoring, the candidate grids and the protocol v1 helpers |
| [`tests/`](tests/) | Checks that our distance matches MICA's, that subsampling the distance matrix is exact, and that the protocol v1 pieces behave |
| [`docs/`](docs/) | The protocol, a primer, dataset cards and the team's trackers. See `docs/README.md` |
| [`notebooks/`](notebooks/) | Notebooks that rebuild and show the Assignment 2 and Assignment 3 tables |
| [`envs/`](envs/) | Conda environments: `scminer` (main), `chooser-r` (chooseR's R code) and a record of `mica-legacy` (the authors' MICA) |

## Rerunning the experiments

These take hours to days. Everything runs on CPU. Start long runs in a detached byobu session with a log.

```bash
bash scripts/setup/setup_envs.sh mica        # the authors' MICA 1.0.0 in its own environment
bash scripts/fetch_chooser.sh                # chooseR's R code at commit 410bdbf, checksummed
conda env create -f envs/chooser-r.yml
bash scripts/run_protocol_v1.sh small 3      # nine small datasets, one seed (about 2 h)
bash scripts/run_protocol_v1.sh pbmc 3 full  # PBMC14K, one seed, in four stages: full, then
bash scripts/run_protocol_v1.sh pbmc 3 mica  #   mica (8 to 10 h) and pca (about 2.5 h),
bash scripts/run_protocol_v1.sh pbmc 3 pca   #   which can run side by side,
bash scripts/run_protocol_v1.sh pbmc 3 score #   then score
bash scripts/run_protocol_v1.sh chooser 3    # chooseR's own code on the nine small datasets
bash scripts/run_protocol_v1.sh chooser-pbmc 3   # and on PBMC14K
```

Protocol v1 runs refuse to start unless the protocol code matches the tag `matched-protocol-v1`, and
each records that check in `results/v1/provenance/`. A small dataset that already has a finished run is
skipped, so in this checkout the `small` command reruns nothing.

The Assignment 3 runs used `scripts/run_rules_small.sh`, `scripts/run_rules_repeats.sh`,
`scripts/run_rules_pbmc_graph.sh` and `scripts/run_rules_matched.sh`. Scripts that take arguments say
how to call them in their first lines.

## Data

The datasets are not in this repository. They live on the lab workstation's shared disk at
`/media/20TB/shared/sides_data/43.scMINER-Clustering/`. Git ignores the `data` symlink, so a fresh
clone needs it made once, from the repository root:

```bash
ln -s /media/20TB/shared/sides_data/43.scMINER-Clustering data
```

The ten benchmark inputs are the scMINER authors' processed release
([Zenodo 15040179](https://doi.org/10.5281/zenodo.15040179)). No script downloads them. Their
`*.Processed_eset.rds` files go in `data/raw/paper_inputs/`, and `scripts/setup/prepare_paper_inputs.py`
converts them without changing a value. The protocol v1 subsample partitions (11 GB) are in
`data/runs/protocol_v1/` on the same disk and are shared on request.

## Versions

| Tag | What it marks |
|---|---|
| `a3-submitted` | The Assignment 3 submission |
| `matched-protocol-v1` | Protocol v1, fixed before the confirmatory runs |
| `a4-submitted` | The Assignment 4 submission |
