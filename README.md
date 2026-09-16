# scMINER clustering

Course project for Clustering Algorithms (Missouri S&T, Fall 2026). We study scMINER, a mutual information based framework for clustering single-cell RNA data ([Pan et al., Nature Communications 2025](https://doi.org/10.1038/s41467-025-59620-6)). We test its clustering part (MICA) on three dataset regimes and extend it to data with millions of cells.

Team: Utshab Kumar Ghosh (lead), Laraib Khan and Rhea.

## Layout

- `envs/`: conda environment files
- `scripts/`: download and run scripts
- `src/`: project code
- `results/`: small tables and figures that everyone can browse on GitHub
- `docs/research/`: literature review and notes

Datasets live outside the repo in `~/data/scminer-clustering/`.

## Running jobs

Every run happens on the lab workstation inside a conda environment. Long jobs go in a detached byobu session and write a log file:

```bash
byobu new-session -d -s scm-<task> "bash -lc 'source ~/anaconda3/etc/profile.d/conda.sh && conda activate scminer && <command> 2>&1 | tee logs/<task>.log'"
```

Use `byobu list-sessions` to see running jobs and `byobu attach -t scm-<task>` to look at one.

Setup steps for the environments will be added here once they are built.
