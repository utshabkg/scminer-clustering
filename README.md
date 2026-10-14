# scMINER clustering

Course project for Clustering Algorithms (Missouri S&T, Fall 2026). We study scMINER, a mutual information based framework for clustering single-cell RNA data ([Pan et al., Nature Communications 2025](https://doi.org/10.1038/s41467-025-59620-6)). Our main question: when each method must choose its own clustering without the reference labels, does scMINER's clustering (MICA) beat a standard pipeline of principal components and Leiden? The paper's benchmark lets the labels pick every method's setting, so we ask what is left without them. Scaling the mutual information graph and annotating clusters with a language model follow that question.

Team: Utshab Kumar Ghosh (lead), Laraib Khan and Rhea Eckman.

## Layout

- `envs/`: conda environment files
- `scripts/`: download and run scripts
- `src/`: project code
- `results/`: small tables and figures that everyone can browse on GitHub
- `notebooks/`: notebooks that rebuild the report's tables and figures from `results/`
- `latex/`: the current report, with the submitted proposal in `latex/archive/proposal/`
- `docs/dataset_cards/`: one page per dataset
- `docs/research/`: literature review and notes

Datasets live outside the repo, on the lab workstation's shared disk at `/media/20TB/shared/sides_data/43.scMINER-Clustering/`. The `data` link in the repo points there. They are never stored on GitHub.

## Running jobs

Every run happens on the lab workstation inside a conda environment. Long jobs go in a detached byobu session and write a log file:

```bash
byobu new-session -d -s scm-<task> "bash -lc 'source ~/anaconda3/etc/profile.d/conda.sh && conda activate scminer && <command> 2>&1 | tee logs/<task>.log'"
```

Use `byobu list-sessions` to see running jobs and `byobu attach -t scm-<task>` to look at one.

Setup steps for the environments will be added here once they are built.

## Where to look first

- [WORKPLAN.md](WORKPLAN.md): who does what, by when, and how a task counts as finished.
- [PROGRESS.md](PROGRESS.md): what has actually been done, newest first.
- `notebooks/assignment2_preliminary.ipynb`: rebuilds every table and figure in the Assignment 2 report from saved results.
