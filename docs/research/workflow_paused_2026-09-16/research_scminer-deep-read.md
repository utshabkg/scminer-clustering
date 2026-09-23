# scMINER and MICA deep read (checked 2026-09-16)

**Primary sources used throughout**
- Pan Q et al., Nature Communications 16:4305 (2025), https://doi.org/10.1038/s41467-025-59620-6 (open access at PMC12062461).
- Supplementary Information: https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-59620-6/MediaObjects/41467_2025_59620_MOESM1_ESM.pdf
- Peer Review File (MOESM3) and Source Data zip (MOESM4), same path with the names swapped.
- Preprint: bioRxiv 2023, https://doi.org/10.1101/2023.01.26.523391
- Code: https://github.com/jyyulab/MICA (v1.0.0 and `main`), https://github.com/jyyulab/scMINER (R package v1.1.0), https://github.com/jyyulab/SJARACNe, https://github.com/jyyulab/NetBID, https://github.com/jyyulab/scMINER_analysis. Docs: https://jyyulab.github.io/scMINER/bookdown/
- All four repos were migrated on 2026-07-12. Each has a "Migration Log" issue, and the issue history was kept.

## 0. Bottom line

- MICA uses no new clustering algorithm. It computes a binned mutual information (MI) distance between cells over all genes. It then runs either classical MDS plus k-means (under 5,000 cells) or an HNSW MI-kNN graph, node2vec, a Euclidean kNN graph and Louvain (5,000 cells or more). Reviewer #2 wrote that the paper "does not propose fundamentally new methods" (Peer Review File).
- The accuracy numbers use an oracle. In MDS mode k was set to the true number of classes. For MICA-GE and the baselines, resolution was swept from 0.1 to 4 until the cluster count matched the truth, keeping the highest ARI.
- The largest run was Covid650K (647,366 cells): 13 h and 204.4 GB peak on 40 CPU cores. MICA has no GPU path.
- The `main` code would not run 2M cells as shipped. It densifies an n×n kNN adjacency, which is about 32 TB at n=2M. Even the sparse branches would likely take over 12 h and use 450 to 650 GB (estimate).
- I found no 2025 or 2026 paper that extends MICA's clustering itself. The citing papers apply the network and activity side, extend it (APNet, lncAPNet), or cite it in reviews and new clustering methods.

## 1. MICA clustering pipeline

### 1.1 Preprocessing (scMINER R side)
- **QC.** Metrics are nUMI, nFeature, pctMito (`MT-|mt-`), pctSpikeIn (`ERCC-`) and per-gene nCell. The "auto" cutoffs are Median ± 3×MAD (`drawSparseEsetQC()`, `filterSparseEset()`).
- **Normalization.** Counts go to CPM and then log2 (`normalizeSparseEset()`). TPM datasets were only log2-transformed.
- **No feature selection.** MICA uses all genes that pass QC, and the paper argues against HVG selection.
- **Export.** `generateMICAinput()` writes .txt or .h5ad. It first calls `as.matrix()`, which makes the matrix dense in R, and it has an optional `downSample_N`.

### 1.2 How MI between two cells is estimated
Source: Methods Eq. 1 to 3, plus `MICA/lib/distance.py`, `MICA/lib/mutual_info.py` and `mihnsw/hnswlib/space_mi.h`.
- **What the samples are.** Each cell is a vector over G genes. The MI between cells i and j treats the G genes as paired samples of two discrete variables.
- **Binning.** Bins are equal-width, computed separately for each cell between its own min and max expression. The number of bins is b = floor(G^(1/3)), so 20,000 genes give 27 bins. The flags `-bsz` (fixed b) and `-bpr` (power) override this.
- **Estimator.** A plug-in estimate from the 2D histogram, with no bias correction and no explicit handling of zeros.
- **Speed.** The GE path runs in a C++ HNSW fork ("mihnsw") with a custom `space='mi'`. The MDS path is a pure-Python double loop over `np.histogram2d` (`mi_calc`) and is effectively single-threaded.

### 1.3 MI to distance
- **Paper (Eq. 3).** D = 1 − I(Ci;Cj)/H(Ci,Cj). This matches the C++ HNSW code, which computes `2 − (hx+hy)/hxy`, and the numba distance.
- **Discrepancy.** In the MDS path of v1.0.0 and `main`, `mi_norm` computes D = 1 − I/sqrt(H_i·H_j) instead. That is a different normalization from the one in the paper. This comes from my reading of the code and was not tested.

### 1.4 MICA-MDS (small data)
- **Steps.** The full n×n MI distance matrix goes into classical MDS: double-centering, `np.linalg.eigh`, then keep the positive eigenvalues (`transform.py`). The first 19 dimensions are kept by default, and the paper says ARI plateaus from about 5 dimensions (Fig. 2f,g). k-means runs on that embedding.
- **Consensus.** The paper says k-means runs 10 times and the runs are combined by consensus clustering. In the v1.0.0 code the default is a single `sklearn.KMeans(n_clusters=k)`. With sklearn 0.24.1 that means n_init=10 and the best inertia is kept, which is not consensus.
  - A consensus over several `-dd` dimension values (co-association matrix plus average-linkage hierarchical clustering) was added on 2025-03-31 (commit dfe3ddc0).
  - The legacy 0.2.x code did SC3-style consensus over dimensions and bootstraps.
- **Choosing k.** The user must supply k (`-nck`, default 4). A multi-k option was added in 2024. The wrapper suggests `-nck 4 5 6 7 8 9 10 -dd 20`, and an optional silhouette score (`-sil`) can help pick among them.

### 1.5 MICA-GE (large data), `mica_ge.py` defaults
1. Approximate kNN graph under the MI distance with HNSW: k=80 (`-nnm`), ef=400, M=8.
2. node2vec through pecanpy (SparseOTF, weighted): 20 dimensions, walk length 60, 110 walks per node, window 10, p=2.8, q=0.4.
3. Exact kNN on the 20-D embedding (sklearn ball_tree, k=20, Euclidean or cosine).
4. Louvain with python-louvain (`community.best_partition`). Default resolution is 1.822, and a range can be swept in parallel (`-minr/-maxr/-ss`, `-nw` workers).
5. Optional consensus (`-cs CSPA`), off by default.
6. UMAP (min_dist 0.6) or t-SNE for display.

The mode switch is in `mica auto`: GE when cells ≥ `-cc` (5,000), MDS otherwise. The authors call the 5,000 threshold "an empirical guideline rather than a strict cutoff" (Peer Review File).

### 1.6 Two things I inferred from the code (not tested)
- **Distances used as weights.** In GE, MI distances are written as edge weights for weighted node2vec, and `kneighbors_graph(mode='distance')` feeds Euclidean distances to Louvain as weights. Farther neighbors therefore get larger weights. The effect on accuracy is not measured anywhere.
- **Dense adjacency.** Both `main` and tag 1.0.0 run `kneighbors_graph(...).toarray()` and then `nx.from_numpy_matrix`, which builds a dense n×n matrix. At 650K cells that is about 3.4 TB. It cannot be the code that produced the published 153 to 204 GB Covid650K run.
  - Branches `huge_dataset` and `huge_floating_k` (Dec 2024 to Jan 2025) switch to sparse matrices and float32. PR #29 was closed and PR #30 is still open.
  - So the exact code behind the large benchmark is UNVERIFIED.

## 2. SJARACNe networks and NetBID-style activity and hidden drivers

**Per-cluster network inputs** (`generateSJARACNeInput()`)
- One folder per group (a cluster or an annotated cell type). Each folder holds an expression file (genes × cells) plus hub lists.
- Hub lists are curated: about 1,500 TFs and more than 5,500 signaling genes (SIGs) for human or mouse (`inst/RData/tf_sigs_*.RData`). In the PBMC14K demo, the B-cell group had 835 TFs and 4,148 SIGs present among 8,572 genes.
- Default `downSample_N = 1000` cells per group. Optional SuperCell metacells: the paper used gamma=20 for benchmarks and superCell_gamma=50 for activity clustering.
- The code warns below 100 single cells or 50 metacells per group. The Discussion states a minimum of 500 cells overall.

**SJARACNe algorithm** (Methods, and Khatamian A et al., Bioinformatics 2019, https://doi.org/10.1093/bioinformatics/bty907)
1. Gene-gene MI by adaptive partitioning, as in ARACNe-AP (Lachmann A et al., Bioinformatics 2016, https://doi.org/10.1093/bioinformatics/btw216).
2. Edges filtered at `-pb` 1e-7, then Data Processing Inequality pruning.
3. `-n 100` bootstraps and a consensus network kept at `-pc`. The cutoff is 1e-2 for single cells, 1e-3 for metacells and 1e-5 for bulk data (docs).
4. Edges are annotated with Spearman and Pearson correlation, slope and p-value, and each target is signed.
- It runs as a C++ core behind a Python and CWL wrapper, locally with cwltool (needs node.js) or on LSF with cwlexec.

**Activity** (`getActivity_individual` and `getActivity_inBatch`)
- Expression is z-normalized per column. Activity is then the mean of each driver's targets.
- Methods: "mean" (default), "weightedmean" (MI × sign of Spearman, which is the NetBID2 default), "absmean" and "maxmean".
- Per-cluster activity matrices are concatenated, and missing values are filled with the minimum non-NA value.
- **Hidden drivers** are drivers that change in differential activity (`getDA()`, default t-test, with limma and Wilcoxon also available) but not in expression.
- scMINER does not use NetBID2's Bayesian BID model (Dong X et al., Nat Commun 14:2581, 2023, https://doi.org/10.1038/s41467-023-38335-6).

**Compute cost (not reported in the paper)**
- Work scales with hubs × genes × bootstraps, and each adaptive-partitioning MI costs about O(m log m) for m cells.
- The PBMC14K demo LSF script reserves 8 GB and 1 slot.
- Users report much more:
  - A 2025 comment on SJARACNe issue #50 says 1,361 cells on 32 cores needed 160 to 200 GB.
  - A ~300K-cell run failed with SIGTERM, and the maintainers' answer was to downsample.
  - A 2019 maintainer note (issue #3) puts 28K genes × 1,000 samples at about 8 GB.
- **Calibration problem (new, unmerged).** Draft PR #73 (Aug 2026, developer "adamdingliang") argues the legacy consensus p-value uses an uncorrected normal tail. On BRCA100 it understates the tail about 94-fold for TF and 452-fold for SIG networks. The PR proposes a fixed minimum-recurrence rule (K_edge=6), and PR #71 adds estimator-matched null models. This shows the default `-pc` threshold is not a calibrated p-value.

## 3. Complexity, runtime, memory, and 2M-cell feasibility

**Asymptotics (my analysis of the code)**
- **MDS mode.** O(n²·G) MI evaluations, O(n²) memory and an O(n³) eigendecomposition.
- **GE mode.** About O(n log n) HNSW distance calls, each O(G + b²). On top of that come the random-walk tokens (n × 110 × 60), a kNN graph on the 20-D embedding, and one Louvain run per resolution in the sweep.
- **Dense input.** Both modes load the whole cells × genes matrix dense (`adata.to_df()`), and HNSW stores another float32 copy.

**Published numbers (Source Data, Supp. Fig. 4, 40-core CPU node)**

| Dataset (cells) | MICA time / peak RAM | Scanpy | Seurat | SC3s | scVI | scDeepCluster |
|---|---|---|---|---|---|---|
| PBMC14K (13,605) | 336 s / 2.39 GB | 195 s / 1.26 | 137 s / 3.63 | 287 s / 1.79 | 1,444 s / 1.12 | 1,472 s / 6.59 |
| Klein (2,717, MDS) | 3,112 s / 0.85 GB | 43 s | 136 s | 79 s | 1,946 s | 7,854 s |
| HMC76K (76,533) | 0.7 h / 27.7 GB | 0.1 h / 27.3 | 0.3 h / 59.9 | 0.3 h / 18.2 | 0.4 h / 27.4 | 9 h / 111.4 |
| Covid97K (97,499) | 1.3 h / 26.6 GB | 0.1 h / 2.7 | 0.2 h / 13.6 | 0.4 h / 2.8 | 0.4 h / 2.7 | 2.5 h / 115.1 |
| Covid650K (647,366) | 13 h / 204.4 GB | 0.5 h / 15.1 | 1.1 h / 74.8 | 7.6 h / 16 | 4.7 h / 15 | 89.9 h / 869.2 |

Notes on these numbers:
- MICA used GPU for nothing. The baseline hardware is described only as CPU cores (scDeepCluster was run with a `_cpu` script).
- The repo copy `scMINER_analysis/.../scMINER_RunSummary.txt` gives different figures for Covid650K: 45,497 s and 153 GB. Those were one-shot runs at a single resolution.
- Answering a reviewer on 5,250 subsampled PBMCs, the authors reported MDS at ARI 0.833, 4,332 s and 1,235 MB, and GE at ARI 0.828, 734 s and 3,242 MB.
- The only downsampling in the paper is for SJARACNe (1,000 cells per type for GSE164378) and metacells. The large datasets were only checked for speed and visually (Sankey plots), with no ARI. An older issue on scMINER (#14, 2020) reports an 89K-cell run killed at 512 GB, with an earlier version.

**Would it run on 2,000,000 cells on the project server? Not as shipped. Probably not within 12 h even after fixes. All figures below are estimates.**
1. **`main` or v1.0.0.** The dense `.toarray()` on the final kNN graph needs about 2e6² × 8 bytes, or 32 TB. It fails at once. The R export `as.matrix()` for 2M × 20K genes needs about 320 GB and writes a huge text file.
2. **With the sparse `huge_*` branch and a Python-side h5ad.**
   - Memory: the dense float32 cell × gene matrix is about 160 GB at G=20K. HNSW adds a copy (about 160 GB) plus its links, before walks, networkx graphs and Louvain worker processes. Scaling Covid650K's 153 to 204 GB by 3.1× gives roughly 470 to 630 GB, near the 754 GB limit.
   - Time: n log n scaling from 12.6 to 13 h on 40 cores gives about 40 h on 40 cores. Even perfect scaling to 120 cores gives about 13 to 14 h, and python-louvain is single-threaded for each resolution. So a single run exceeds the 12 h window, before any resolution sweep.
3. **Workable options.**
   - Sketch or downsample to 100K to 300K cells, cluster, then propagate labels by kNN. scMINER has `downSample_N` but no propagation step.
   - Re-implement the MI kNN (see section 8).

## 4. Benchmarks

**Clustering datasets** (Supp. Table 1; preprocessed data are said to be on Zenodo at https://doi.org/10.5281/zenodo.15040179, whose contents are UNVERIFIED because Zenodo returned 403)

| Dataset | Protocol | Cells | Classes | Accession |
|---|---|---|---|---|
| Yan | Tang | 124 | 8 | GSE36552 |
| Goolam | Smart-Seq2 | 124 | 5 | E-MTAB-3321 |
| Buettner | C1 | 182 | 3 | E-MTAB-2805 |
| Pollen | SMARTer | 249 | 11 | SRP041736 |
| Chung | SMARTer | 317 | 5 | GSE75688 |
| Usoskin | STRT-seq | 622 | 4 | GSE59739 |
| Kolod | SMARTer | 704 | 3 | E-MTAB-2600 |
| Klein | inDrop | 2,717 | 4 | GSE65525 |
| Zeisel | STRT-seq | 2,903 | 7 | GSE60361 |
| PBMC14K (Zheng FACS-sorted) | 10x | 13,605 | 7 | SRP073767 |
| HMC76K (speed only) | 10x | 76,533 | 20 | Azimuth |
| Covid97K (speed only) | 10x | 97,499 | 28 | E-MTAB-10026 |
| Covid650K (speed only) | 10x | 647,366 | 51 | EGAD00001007982 (EGA, probably controlled access, UNVERIFIED) |

- Gene counts are not in the table. The PBMC14K source has 21,592 genes (docs).
- Labels come from the original studies. PBMC14K labels come from FACS sorting, with 3 of the 10 sorted populations dropped.
- The run scripts set k to 34 for Covid650K and 4 for Chung. Those disagree with Supp. Table 1, which lists 51 and 5.

**Baselines and settings.** Seurat 5.0.1 (Louvain), Scanpy 1.10.2 (Leiden), SC3s 0.1.1, scVI 1.0.3 (Leiden) and scDeepCluster (PyTorch, consensus k-means). Every baseline used 2,000 HVGs, and the oracle k or resolution rule from section 0.

**Metrics.** ARI, AMI, NMI, ASW, AvgBIO (mean of ARI, NMI, ASW), accuracy and purity.

**Headline clustering results** (mean over the 10 datasets, from Source Data)
- **ARI.** scMINER 0.844 (SD 0.151), Seurat 0.733, scVI 0.695, Scanpy 0.679, SC3s 0.679, scDeepCluster 0.654. The paired Wilcoxon test against Seurat gives p = 0.027 with n = 10.
- **Per-dataset ARI.** Yan 0.939, Goolam 0.544, Buettner 0.934, Pollen 0.942, Chung 0.655, Usoskin 0.847, Kolod 1.0, Klein 0.98, Zeisel 0.744, PBMC14K 0.853. scMINER lost on Usoskin (Seurat 0.94), Goolam (SC3s 0.678) and Pollen (scDeepCluster 0.943).
- **Other metrics against Seurat.**
  - AMI 0.862 vs 0.784 (p = 0.058) and NMI 0.866 vs 0.791 (p = 0.058). Neither is significant.
  - ASW is higher for Seurat (0.777 vs 0.722).
  - AvgBIO 0.810 vs 0.767 (p = 0.375).
- **PBMC14K T cells vs Seurat.** Treg purity 96% vs 78% and TCM purity 92% vs 61%. Treg accuracy 70% vs 52% and TCM accuracy 95% vs 88%.
- **Robustness.** Bin size from 5 to 100 had little effect in MDS mode. In GE mode, ARI depends mostly on Louvain resolution: on Zeisel it ranged from about 0.12 to 0.81 across the sweep.
- **Source Data error.** The Fig. 2g sheet reports two mean ARIs above 1 (1.657 and 1.61).

**Network and activity results**
- **Networks.** On a Perturb-seq knockout set (GSE221321), AUROC was 0.69, against GENIE3 0.65, GRNBoost2 0.57 and PIDC 0.52. Predicted drivers were correct for 17 of 18 exhausted T-cell TFs (94%), against PIDC 78%, GENIE3 39% and GRNBoost2 28%. The SIG network reached AUROC 0.75, and 207 of 294 SIGs showed significant GSEA enrichment.
- **Other ground truth.** ATAC-seq footprinting (GSE123236, GSE112731), CRISPRa/i (GSE218988) and CITE-seq (GSE213282, GSE164378).
- **Activity-based clustering ARI.** PBMC14K 0.87 vs SCENIC 0.69, and Treg dataset 0.998.

## 5. Software stack and install issues

**Components**
- **R package.** scMINER needs R 4.2.3 or newer and installs with `devtools::install_github("jyyulab/scMINER")`.
- **Python tools.** MICA and SJARACNe need Python 3.9.2 in conda, installed from source (`pip install .` for MICA, `python setup.py build/install` for SJARACNe).
- **Stale PyPI packages.** PyPI "MICA" is still 0.2.2 (Python ==3.7.6), so `pip install MICA` does not give v1.0.0. PyPI SJARACNe is 0.2.1.
- **Pinned dependencies.** MICA's `requirements.txt` pins anndata 0.8.0, pandas 1.2.3, scipy 1.7.3, scikit-learn 0.24.1, numba 0.53.1, networkx<3 (`from_numpy_matrix` was removed in networkx 3), pecanpy 1.0.1, python-louvain 0.15 and gensim 4.1.2. Python 3.9 reached end of life in October 2025.
- **Other requirements.** `setup.py` compiles the bundled `mihnsw`. SJARACNe needs node.js plus cwltool, or cwlexec for LSF.
- **Input formats.**
  - scMINER reads the 10x directory, 10x HDF5, h5ad (through the anndata R package) and text tables into a SparseEset (a sparse Biobase ExpressionSet).
  - MICA reads .txt, .csv or .h5ad with genes as rows by default (changed in June 2024).

**Known issues**
- The R and Python versions conflict in one conda env (scMINER #19, 2025). The maintainer advised separate environments.
- anndata 0.8.0 cannot read newer .h5ad files (MICA #33, open).
- Errors importing mihnsw ("Illegal instruction") depend on gcc and conda versions (MICA README).
- Apple M-series machines need a separate branch.
- SJARACNe hits "Did not find output file ... TF_run_001.adj" repeatedly (#44, #51, #52, #57), bootstraps get killed for lack of memory (#50), `-pb` 1e-7 was rounded to zero (#53, fixed 2025-11), and the DPI tolerance still cannot be set (#8, open since 2019).
- A MICAplot typo (#20) was fixed on 2026-07-23.
- **No GPU support anywhere.** A PyTorch Geometric node2vec was merged in 2023 but is commented out in `main`.

## 6. Limitations

**Stated in the Discussion**
1. It needs at least 500 cells, which limits rare-population studies. Suggested workarounds are signal amplification (computational or experimental) or mapping to references such as the Human Cell Atlas.
2. Activity is only estimated for the curated TF and SIG lists, so new non-TF, non-SIG drivers are missed. The proposed future work is to use both upstream and downstream first-order neighbors.
3. There is no ground truth for large datasets, so the accuracy benchmark covers only small datasets.
- The reviewers also raised the unclear rationale for the mode split and how to pick resolution without labels. The authors answered with resolution sweeps, silhouette scores and marker inspection. The docs pick PBMC14K resolution 2.05 "based on the knowledge about" the dataset.

**Inferred (from code and design, not tested)**
- **Dropout and zeros.** Equal-width per-cell bins put zeros and low values in bin 0, so the (0,0) cell of the joint histogram dominates. MI then largely reflects co-detection patterns and sequencing depth. There is no dropout model or bias correction, and MI between cells ignores gene identity beyond pairing.
- **Low-dimensional data breaks the estimator.** With p=10 or 20 features, b = floor(p^(1/3)) = 2 bins, and MI is estimated from 10 to 20 samples. That is close to meaningless. MICA's cell-to-cell MI only makes sense for high-dimensional data such as Dataset 3. Datasets 1 and 2 would need an MI over features instead.
- **Cost.** MDS mode is quadratic in time and memory with a cubic eigensolve. Both modes keep the whole matrix dense and use all genes.
- **Resolution and k.** There is no automatic model selection, and published accuracy depends on oracle k or resolution. python-louvain can produce badly connected communities, which Leiden fixes (Traag VA et al., Sci Rep 2019, https://doi.org/10.1038/s41598-019-41695-z).
- **Weights, seeds, batches.** Distances are used as weights (section 1.6). No random seeds are set for KMeans, Louvain, HNSW or walks. Harmony batch correction is commented out as "does not work well".
- **SJARACNe cutoff.** The consensus threshold is miscalibrated, per the developers' own PR #73.

## 7. Papers from 2025 to 2026 citing or extending scMINER or MICA
The lists come from Europe PMC "cited by" (13), OpenAlex (15) and Semantic Scholar (15 for the paper, 5 for the preprint). Google Scholar was not queried and my web-search budget ran out, so completeness is UNVERIFIED.

- **Extensions of the activity and network logic:**
  - APNet (Bioinformatics 2025): https://doi.org/10.1093/bioinformatics/btaf063
  - lncAPNet (Bioinformatics Advances 2026), which explicitly uses "NetBID2's and scMINER's activity logic" on lncRNA SJARACNe networks: https://doi.org/10.1093/bioadv/vbag204
  - SINBA drug-synergy inference (bioRxiv 2025): https://doi.org/10.1101/2025.11.20.689490
- **Protocol:** human B-cell gene activity landscape (STAR Protocols 2025): https://doi.org/10.1016/j.xpro.2025.103614
- **Applications, mostly St. Jude:**
  - Genes & Development 2025: https://doi.org/10.1101/gad.352632.125
  - Cancer Discovery 2025: https://doi.org/10.1158/2159-8290.cd-24-0841
  - Immunity 2026: https://doi.org/10.1016/j.immuni.2026.03.008
  - Cell Death & Disease 2026: https://doi.org/10.1038/s41419-026-08555-x
  - iScience 2026: https://doi.org/10.1016/j.isci.2026.117410
  - J Clin Invest 2026: https://doi.org/10.1172/jci196753
  - J Transl Med 2026: https://doi.org/10.1186/s12967-026-08098-x
  - bioRxiv 2026 on synthetic T-cell states: https://doi.org/10.64898/2026.01.07.698190
- **Related Yu-lab methods:** Spotiphy (Nature Methods 2025, https://doi.org/10.1038/s41592-025-02622-5) and transcriptome-size normalization (Nat Commun 2025, https://doi.org/10.1038/s41467-025-56623-1).
- **Clustering and MI methods or reviews that cite it:**
  - scHashFormer (Brief Bioinform 2026): https://doi.org/10.1093/bib/bbag402
  - Survey of single-cell clustering by encoder type (Arch Comput Methods Eng 2026): https://doi.org/10.1007/s11831-026-10709-1
  - MI-POG (Cancer Genomics Proteomics 2026): https://doi.org/10.21873/cgp.20606
  - Trends in Genetics 2026 review on fate TFs: https://doi.org/10.1016/j.tig.2026.05.005
  - SYMOL (Genome Research 2026): https://doi.org/10.1101/gr.281603.125
  - Whether scHashFormer benchmarks against MICA is UNVERIFIED.
- **Nothing found** that improves MICA's MI estimator, its scalability, or its choice of resolution.

## 8. What this means for the project
- **Framing.** "MI distance plus MDS or graph-embedding clustering" fits the course's "adapts a clustering method" rule. Novelty for a paper would have to come from the gaps above, not from running scMINER as-is.
- **Openings I think are defensible (not yet tested):**
  1. **Scale.** Pre-bin each cell once (uint8 codes) and compute the joint histogram only over genes that are nonzero in either cell, treating the rest as the zero bin (strictly, only when zeros are the only values in bin 0). Batch this on GPU to make an exact or HNSW MI-kNN that reaches 2M cells in under 12 h.
  2. **Dropout-aware estimators.** Try excluding genes that are zero in both cells, quantile bins, or bias-corrected or KSG-style MI, and measure the effect on ARI.
  3. **Label-free resolution selection** to replace the oracle, evaluated on atlases with labels.
  4. **Weight semantics.** Fix the distance-as-weight choice and ablate its effect.
- **RAG link.** Cluster-specific hidden-driver lists (drivers with differential activity but flat expression) are a natural query set for the planned literature-retrieval step.