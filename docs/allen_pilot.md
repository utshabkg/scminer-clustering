# Allen atlas pilot: 100,000 cells (U-16)

Run on Oct 2 to 3, CPU only, on the shared workstation (MICA with 8 workers). The sample is 100,000 cells drawn at random (seed 0) from the 13 raw 10Xv3 region files, scored against the 34 cell classes. After filtering, 28,536 genes remain. The raw files are listed with their sizes in `docs/allen_files.csv` (97 GB in all).

## What we measured

| Run | Wall time | Peak memory | ARI (class level) |
|---|---|---|---|
| Prepare: read, sample, filter, normalize | 7 min | 24 GB | |
| Leiden, 2,000 genes (seurat_v3) | 23 min | 27 GB | 0.661 |
| Leiden, all genes | 25 min | 55 GB | 0.624 |
| MICA graph mode, all genes | 7.1 h | 140 GB | 0.652 |

Most of MICA's time goes to one step:

| MICA step | Time |
|---|---|
| Mutual information neighbor graph (HNSW, 80 neighbors) | 5.3 h |
| node2vec embedding | 40 min |
| Louvain over 40 settings | 39 min |
| UMAP plots | 72 min |

The Leiden times include the full sweep of 200 resolutions, which alone takes about 17 minutes.

## What 2.3 million cells would cost

The full atlas has 23 times as many cells. The numbers below are rough. They assume the neighbor search grows like n log n (as HNSW does), the other steps grow in proportion to the number of cells, and so does memory.

- **MICA time.** The neighbor graph would take about 155 hours (5.3 h × 23 × 1.27). node2vec and Louvain would add about 15 hours each. That comes to 7 to 8 days without the plots. (Added Oct 28: MICA's search alone took 3.8 h on the same cells when rerun for the graph comparison, which gives about 110 hours for the graph and about 6 days in all, so the paper now says 6 to 8 days.)
- **MICA memory.** About 3.2 TB if it grows in proportion. The workstation has 754 GB, shared with other users. Part of the reason is that MICA holds the whole matrix as dense numbers. Its 100,000 cell output file is 11.4 GB, which is exactly a dense 100,000 by 28,536 table, and one such copy at 2.3 million cells is about 260 GB.
- **Leiden.** About 9 to 10 hours for the whole run. Memory would be about 600 GB with 2,000 genes and 1.3 TB with all genes, because scaling the genes before PCA turns the sparse matrix dense. So the baseline also needs a sparse PCA at full size.

So MICA as shipped cannot run on the whole atlas on this machine. Time is a problem and memory is a hard stop. This is why RQ2 goes after the neighbor search, which is the largest time cost (U-18, U-19), and clusters a sketch of the atlas before assigning the remaining cells, which deals with memory (U-20).

Logs: `logs/allen_pilot_*.log` and `logs/mica_allen_100k.log`. Scores: `results/pilot_mica_allen_100k/`, `results/baselines_allen_100k*/`.
