# Protocol v1: matched candidates, fixed before the confirmatory runs

Status: written on 2026-11-12, before any confirmatory run. It is committed together with the code that runs it, and that commit carries the annotated tag `matched-protocol-v1`. Every confirmatory job (`--require-tag`) refuses to start unless the files it executes are exactly the tagged ones: `src/scminer_clust`, the protocol scripts, the R scorer, the R environment files and this document (`protocol.CODE`). Notes, the manuscript and the analysis script can change while the runs go on. Each run records the commit, whether its code matched the tag, and its package versions in a provenance file. If anything here has to change, the change goes in the deviations log at the end and gets a new tag (`matched-protocol-v1.1` and so on).

## Why this exists

Jian Liu's comments on Assignment 3 accepted the matched-candidate design as an exploratory improvement. He asked for the main conclusions to be checked on new independent runs under a protocol fixed in advance, with at least five seeds. He also asked for a check of our differences from chooseR, a PCA plus k-means baseline, and a test of whether the result depends on how candidate settings are chosen.

## What had already been seen

- The original design (each method's own grid, capped at 2 to 20 clusters): MICA and PCA plus Leiden, seeds 0 to 2, on all ten datasets. Its design and results arrived in one commit (`384d8f6`), so git alone does not show the order.
- The matched design: same two methods, seeds 0 to 2, on all ten datasets. It was introduced after the original results had been seen, and its grids and results also arrived in one commit (`faca18c`).
- Pre-tag tests of this code, which printed choices for all three methods: Yan seed 0, Kolodziejczyk seed 1, and PBMC14K seed 0 at 3 subsamples. chooseR's own code was also run on the Yan seed 0 partitions.
- Nothing at seeds 3 to 9 on any dataset, by any method or scorer.

## Data

The ten benchmark datasets as released by the scMINER authors, all genes, inputs in `data/processed/paper_<name>.h5ad`. Each run records the sha256 of its input file.

| Dataset | Cells | Genes | Types | Mode for MICA |
|---|---|---|---|---|
| Yan | 124 | 20,286 | 8 | small-data |
| Goolam | 124 | 23,571 | 5 | small-data |
| Buettner | 182 | 8,989 | 3 | small-data |
| Pollen | 249 | 6,982 | 11 | small-data |
| Chung | 317 | 10,897 | 4 | small-data |
| Usoskin | 622 | 16,383 | 4 | small-data |
| Kolodziejczyk | 704 | 13,698 | 3 | small-data |
| Klein | 2,717 | 24,175 | 4 | small-data |
| Zeisel | 2,903 | 14,890 | 7 | small-data |
| PBMC14K | 13,605 | 8,846 | 7 | graph |

## Methods

- **MICA, small-data mode**: the cached mutual information distance, classical MDS to 19 dimensions, k-means with 50 restarts (`scripts/run_rule_stability.py`, `mica_partitions`). Setting: k.
- **MICA, graph mode** (PBMC14K): exact 79 nearest neighbours from the cached distance, then MICA's own node2vec and Louvain (`scripts/run_rule_stability_graph.py`). Setting: Louvain resolution, where higher gives fewer clusters. MICA's node2vec and Louvain take no seed, so these partitions are new draws at every run, seeds 0 to 2 included.
- **PCA plus Leiden**: Scanpy, scaled data, 50 principal components (arpack), 15-neighbour graph, Leiden (igraph, 2 iterations), all seeded. Setting: resolution.
- **PCA plus k-means** (new in v1, the textbook baseline): k-means with 50 restarts on the same 50 principal components as Leiden, computed once per subsample. Setting: k.

## Candidate sets

| Design | MICA | Leiden | k-means |
|---|---|---|---|
| original | k from 2 to 20 (graph mode: MICA's own 0.1 to 4.0 grid), settings giving 2 to 20 clusters | the 0.05 to 2.0 grid, settings giving 2 to 20 clusters | k from 2 to 20 |
| mid (the matched design, primary) | the counts from 2 to 20 that MICA and Leiden both reach on the full data, one setting per count, the middle one in grid order (`grids.one_per_count`) | same | same |
| low, high | as mid, but the first or last setting per count (graph mode only) | first or last setting per count | as mid |
| dense | as mid (a dense Louvain set would cost about 73 CPU-hours per PBMC14K seed) | every matched-grid setting giving a matched count | as mid |

The matched grids are those of Assignment 3 (`grids.LEIDEN_MATCHED`, `grids.MICA_GRAPH_MATCHED`). k-means reaches every count, so it cannot change the matched counts, and each run asserts that it does not.

## The rule and the scorers

Each setting is clustered on the full data and on 100 subsamples of 80 percent of the cells, drawn by `stability.subsamples(n, 100, seed)`. The first 25 of these are exactly the 25 the original protocol drew. Every subsample partition is saved, and the same partitions are scored four ways:

- **S0**, the original protocol: our scorer (`stability.score_setting`, unchanged) on the first 25 subsamples, with co-clustering kept on a random 3,000-cell reference subset when a dataset has more than 3,000 cells (only PBMC14K). Percentile bootstrap over cells within clusters, 2,000 draws. **S0 is primary.**
- **S1**: our scorer on all 100 subsamples, same reference subset.
- **S2** (PBMC14K only): our scorer on all 100 subsamples and all 13,605 cells, original, mid, low and high designs.
- **S3**: chooseR's own R code (github.com/rbpatt2019/chooseR, commit `410bdbf`, fetched by `scripts/fetch_chooser.sh`) on all 100 subsamples and all cells, original and mid designs. Its tally, scoring, silhouette, grouping and `boot_median` (BCa, 25,000 replicates, over the per-cluster average silhouettes) run unchanged in R (`scripts/score_chooser_r.R`). Its choice is ported line for line in `protocol.pick_chooser`.

Two readings of each set of scores:

- **finest stable** (F, chooseR's): the threshold is the highest lower bound among the candidates, and the choice is the setting with the most clusters whose score clears it. Ties go to the higher score in our scorer, and to the higher setting value in chooseR's code.
- **most stable** (M): the highest score, ties to fewer clusters. It uses no interval, so it is the same under S1 and S3 up to float rounding, and S3 reports only F.

When BCa cannot give an interval for a setting (boot.ci stops with an error, or returns no interval because every bootstrap median is equal), S3 is reported three ways, decided now:

- **V**, as chooseR's code behaves: the threshold is undefined, so no choice is made.
- **G** (primary): that setting gets a point interval at its median.
- **G2**: that setting is dropped.

## Seeds

- Nine small datasets: seeds 0 to 9. Seeds 3 to 9 are the confirmatory set. Seeds 0 to 2 are rerun so S1, S3 and k-means exist there too, and their S0 output must reproduce the saved Assignment 3 results.
- PBMC14K: seeds 0 to 4, all new draws for MICA.
- Not rerun, and staying at three seeds or one, as secondary results: the mode-against-size test (three seeds), the silhouette rule (one), the shipped defaults (three) and the transfer comparison (from the full sweeps). The paper says which and why.

## Outcome and summary

ARI against the reference labels (AMI as a secondary measure), read only after the choices are made. A design's summary is the per-dataset median over the seeds in question, then the mean over the nine small datasets. Per-seed means are reported as support.

## Criteria (S0, mid design, nine small datasets, confirmatory seeds 3 to 9, unless stated)

The values from seeds 0 to 2 are in parentheses.

- **R0, reproduction.** Every column but `seconds` of the S0 output for MICA and Leiden at seeds 0 to 2 equals the saved result (`scripts/check_v1_regression.py`). On PBMC14K only original-grid Leiden can reproduce.
- **C1, the matched ranking.** MICA's mean ARI is above Leiden's under F (0.662 against 0.543) and under M (0.473 against 0.452). Holds if the summary mean is higher and MICA leads in at least 5 of the 7 per-seed means. Weak if only the summary mean is higher. Fails otherwise. Judged separately for F and M. We expect M to be fragile, since its gap was 0.021.
- **C2, the reading moves MICA more.** The mean over datasets of the absolute difference between the F and M medians is larger for MICA than for Leiden (0.197 against 0.118).
- **C3, the flat top sits in the smallest datasets.** Pooling each dataset's mid-design curves over the seven seeds, the share of MICA's settings within 0.01 of its best is higher, on average, in the five datasets of at most 317 cells than in the four larger ones (0.534 against 0.231). It is also higher than Leiden's share in those five (0.177).
- **C4, PBMC14K under M.** Seeds 0 to 4: MICA and Leiden are within 0.02 ARI of each other in at least 3 of the 5 seeds.
- **C5, the original-grid reversal** (secondary). On the original grids under M, Leiden's mean is above MICA's (0.635 against 0.473).
- **F1, fidelity to chooseR.** Under S3-G, MICA's F mean is above Leiden's on the mid design, as under S0, for seeds 3 to 9. Also reported, without a verdict: the same for seeds 0 to 9, for k-means, for the original design, and for PBMC14K.
- **FC, the frequency check.** On every dataset of at most 704 cells (every scored setting), on the first 5 settings of Klein and of Zeisel, and on the first 2 of PBMC14K, at every seed, the co-clustering frequency R tallies with chooseR's `find_matches` and `percent_match` equals ours exactly.

## Exploratory, reported with no verdict

- **E1, candidate choice.** S0 picks and means for the low, high and dense Leiden sets against mid.
- **E2, menu against representation.** PCA plus k-means against PCA plus Leiden (same representation, different menus), and MICA against PCA plus k-means (same menu of k, different representation), in both designs.
- **E3, why matching changes the ranking.** The change in Leiden's mean split into one setting per count, the finer grid, and the wider range of counts.
- **E4, implementation details.** S0 against S1 (25 against 100 subsamples), S1 against S2 (reference subset against all cells, PBMC14K), and S1 against S3 (our bootstrap against chooseR's).
- Seed-to-seed ranges, chosen counts, and AMI throughout.

## Fallbacks, decided now

- If PBMC14K falls behind: every PBMC14K seed is cut to the same number of subsamples that all seeds have reached by 2026-11-14 09:00 CST, at least 50. S0 is unaffected. Seeds 3 and 4 go first. S3 then uses the frequency together over seen for the subsamples kept, which is chooseR's `percent_match` with n equal to that number, and R does not tally (its `percent_match` assumes 100).
- If the R environment fails: S1 and S2 are reported, and S3 is reported as not run. A Python port of chooseR's bootstrap may be shown, labelled as such and never as chooseR's code.
- If time runs out, drop in this order: S3 on the low, high and dense sets (not planned anyway), then PBMC14K seeds 0 to 2.
- If R0 fails: up to two hours to find the cause in threads or versions. Otherwise v1 stands as a new protocol, and the paper reports the size of the mismatch.

## Deviations log

None yet.
