# Claims to verify (Oct 3)

Repo: /home/ugdf8/coursework/26Fall/Clustering_Algorithms/scminer-clustering. Raw results: results/*/metrics.csv and params.json; summaries in results/summary_*.md and .csv, results/uncertainty/; tracker docs/results_tracker.csv; published scores docs/published_scores.csv (from the authors' Benchmark_ARI_metrics files in external/scMINER_analysis). MICA raw outputs in results/reproduction/mica_*/ (clustering_*.txt, seeded/seed_*.txt, mi_reduced.pkl, node2vec embeddings). Baseline partitions in data/runs/. "Authors' input" = the authors' processed ExpressionSets from Zenodo 15040179 (data/processed/paper_<name>.h5ad). Scores use the paper's rule (true number of clusters; Leiden's best ARI among resolutions 0.01 to 2.0 that hit it) unless stated.

## Implementation
A1. Our MICA distance (src/scminer_clust/mi.py) matches the authors' functions to 1e-6 on 9,999 random Zeisel cell pairs in both modes (tests/test_mi_matches_mica.py, reference built by scripts/make_mi_reference.py).
A2. MICA's small-data mode uses D = 1 - I/sqrt(H_x H_y) (MICA/lib/mutual_info.py mi_norm), while graph mode and the paper's Eq. 3 use 1 - I/H(x,y) (MICA/lib/distance.py). Only the second is a proven metric.

## Reproduction
R1. On the authors' inputs, once MICA's final k-means reaches its best fit (100 restarts, scripts/mica_seeded_kmeans.py), MICA lands within 0.05 of the published ARI on 9 of 10 datasets, and matches to three decimals on Yan, Goolam, Buettner, Chung, Klein and PBMC14K. Kolodziejczyk is the exception (0.727 against 1.0).
R2. On Yan, Chung, Buettner and Goolam the published MICA score equals the ARI of the lowest-objective k-means fit on MICA's own embedding. MICA's default (10 restarts) can stop on worse fits: our earlier 10-restart draws gave Yan 0.77 to 0.94 and Chung 0.39 to 0.54.
R3. On Usoskin, k-means fits within 0.2 percent of the best objective span ARI 0.79 to 0.87.
R4. Graph-mode MICA is stable across runs: PBMC14K 0.850 to 0.857 over six runs, Zeisel (authors' input) 0.755 to 0.781 over six runs.

## MICA against Leiden on the authors' inputs, ten datasets (results/summary_authors_input_ari.csv, results/uncertainty/)
B1. With all genes, MICA is ahead on 8 of 10 (Klein a tie), mean Leiden minus MICA -0.17, paired Wilcoxon p = 0.027, and the 95 percent bootstrap interval excludes zero on 9 of 10.
B2. With 2,000 genes (any of three selection methods) MICA is ahead on 6 of 10, and the dataset-level test is not significant (p = 0.08 to 0.28).
B3. On the five small datasets (124 to 317 cells) MICA is ahead of Leiden in every gene setting, and the bootstrap interval excludes zero in 17 of 20. On Buettner and Chung MICA beats all three principal-component baselines (k-means, hierarchical, Leiden) in every gene setting. On Goolam, k-means and hierarchical beat MICA.
B4. Our Leiden on the authors' input with 2,000 seurat-method genes reproduces the paper's own Scanpy score closely on several datasets (Goolam 0.572 against 0.572, Kolodziejczyk 0.996 against 0.996, Pollen 0.872 against 0.874, Buettner 0.301 against 0.277, Zeisel 0.657 against 0.630).

## Choosing the number of clusters without labels (results/summary_label_free.*)
C1. Rule: each method keeps the partition with the best silhouette in its own embedding, over k = 2 to 20 or its resolution sweep. Mean ARI over 10 datasets and 4 gene settings: MICA 0.819 told k, 0.577 by silhouette; Leiden 0.692 and 0.556. MICA's lead shrinks from 0.13 to 0.02.
C2. On the small datasets MICA leads when told k (0.81 against 0.56) but trails by silhouette (0.49 against 0.61). On the mid-sized datasets the two tie when told k (0.83 against 0.82) and MICA leads by silhouette (0.66 against 0.51).
C3. 27.5 percent of MICA's silhouette choices land on the edge of the range (k = 2 or 20).

## What the similarity measures (results/summary_diagnostics.md)
D1. The co-detection share of MICA's mutual information (median over 20,000 random pairs) is 0.57 to 0.62 on Usoskin, Kolodziejczyk and Zeisel but 0.25 on Klein.
D2. MICA run on detection only nearly reproduces full MICA on Zeisel (0.698 against 0.706, agreement 0.87) and Kolodziejczyk (0.711 against 0.721, agreement 0.98), but not on Usoskin (0.53 against 0.885) or Klein (0.58 against 0.98).
D3. Thinning counts to 25 percent barely moves MICA on Zeisel, Kolodziejczyk and Klein. On Usoskin MICA falls from 0.885 to 0.57 at 50 percent and 0.27 at 25 percent. Leiden falls on Usoskin only at 25 percent.
D4. Our rebuilt MICA pipeline (our distance, MDS, k-means with 50 restarts) gives the same ARI as MICA's own runs on these datasets (e.g. Usoskin 0.885, Kolodziejczyk 0.721, Klein 0.98, Zeisel 0.706 against 0.705).
