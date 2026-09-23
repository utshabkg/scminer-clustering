# Literature sweep: MI-based similarity and million-cell clustering (as of 2026-09-16)

## 0. What scMINER/MICA actually computes (checked against the paper and the code)

- **Paper:** Pan et al., Nat Commun 2025, https://doi.org/10.1038/s41467-025-59620-6. The cell-cell distance is D = 1 − I(Ci;Cj)/H(Ci,Cj). Genes are the samples. Expression is put into equal-width bins, with b = n_genes^(1/3) bins. MICA-MDS (≤5,000 cells) builds the full MI matrix, runs MDS to 19 dims, then consensus k-means. MICA-GE (>5,000 cells) builds an MI kNN graph (k=80), runs node2vec to 20 dims, builds a Euclidean kNN graph (k=20) and runs Louvain.
- **Accuracy evidence is small-scale only.** The ARI benchmark uses 10 small ground-truth datasets (mean ARI 0.84). MI was compared with Pearson, Spearman and Euclidean on only 4 gold-standard sets (Yan, Pollen, Kolod, Buettner).
  - HMC76K, Covid97K and Covid650K were used for runtime and memory only, on a 40-core CPU. The numbers are only in Supplementary Fig. 4c.
  - The authors admit this: "absence of ground-truth datasets for large-scale datasets limits the accuracy benchmark". They also state that scMINER "requires a minimum of 500 cells".
- **The code differs from the paper** (github.com/jyyulab/MICA, master).
  - The Methods say HNSW (hnswlib). `MICA/mica_ge.py` actually calls `pynndescent.NNDescent(metric=numba_calc_mi_dis, n_neighbors=80)`.
  - `requirements.txt` pins pynndescent 0.5.2 and has no hnswlib. Dependencies are 2021-era (numpy 1.20).
  - The graph is built in networkx and clustered with python-louvain. At millions of cells these are likely bottlenecks (my inference, UNVERIFIED).
- **D is a true metric.** 1 − I/H(X,Y) is the normalized variation of information. Vinh, Epps & Bailey prove it is a metric (JMLR 2010, Theorem 3, https://www.jmlr.org/papers/volume11/vinh10a/vinh10a.pdf). Kraskov et al. state the same (EPL 2005, https://doi.org/10.1209/epl/i2004-10483-y).
  - MICA bins each cell into a partition of the gene set, so D is exactly a metric on those binned cells (a pseudometric on raw cells).
  - This matters for scaling. Triangle-inequality pruning and farthest-first sketching guarantees carry over to MI space.

## 1. Information-theoretic similarity for scRNA-seq

**Benchmarks of similarity metrics**
- **Skinnider, Squair & Foster**, Nat Methods 2019, https://doi.org/10.1038/s41592-019-0372-4.
  - Tested 17 measures. Proportionality (ρp, φs) was best across tasks.
  - MI (computed with WGCNA) "performed relatively poorly on one or more tasks".
  - Cell clustering was tested on only 561 cells from 7 cell lines.
  - This contradicts scMINER's claim that MI is better, and neither study is at scale.
- **Kim et al.**, Brief Bioinform 2019, https://doi.org/10.1093/bib/bby076. Correlation-based metrics beat distance-based ones. Whether MI was in their panel is UNVERIFIED.
- **Watson et al.**, Brief Bioinform 2022, https://doi.org/10.1093/bib/bbac387. Tested 17 proximity metrics. The best one depends on sparsity, rarity and whether structure is discrete or continuous. Euclidean did poorly. MI was not in their list.
- **Ji et al.**, bioRxiv 2023, https://doi.org/10.1101/2023.12.26.572833. Compares distances between cell populations, not cells. MSE and cosine were best.
- **scCluBench** (Xu et al., arXiv Dec 2025, https://arxiv.org/abs/2512.02471). Covers 14 methods and 36 datasets, the largest about 27.5K cells. No MI-based method is included.
- **Yin et al.**, Genome Biol 2025, https://doi.org/10.1186/s13059-025-03719-y. A clustering benchmark. Details are UNVERIFIED because PMC returned a captcha.

**Information-theoretic distances already run at multi-million scale**
- **minicore** (Baker, Dyjack, Braverman, Hicks, Langmead; bioRxiv 2021, https://doi.org/10.1101/2021.03.24.436859, and ACM BCB 2021).
  - Does k-means++, localsearch++ and mini-batch k-means with JSD, KL and Bhattacharyya on sparse counts.
  - Clusters 4.06M cells in minutes using under 10 GiB RAM.
  - This is the closest prior art for information-theoretic distances at multi-million scale. It treats each cell as a distribution and uses divergences, not MI between profiles.

**MI estimators relevant to sparse counts**
- **Plug-in histogram:** MICA.
- **Adaptive partitioning:** Darbellay & Vajda, IEEE TIT 1999, https://doi.org/10.1109/18.761290. Also ARACNe-AP (Lachmann et al., Bioinformatics 2016, https://doi.org/10.1093/bioinformatics/btw216), about 200x faster than fixed bandwidth. Also SJARACNe (Khatamian et al., Bioinformatics 2019, https://doi.org/10.1093/bioinformatics/bty907).
- **KSG kNN estimator:** Kraskov et al., Phys Rev E 2004, https://doi.org/10.1103/PhysRevE.69.066138.
  - Mixed discrete-continuous KSG handles point masses at zero: Gao et al., NeurIPS 2017, https://arxiv.org/abs/1709.06212.
  - Normi used mixed KSG specifically against dropouts (Zeng et al., Brief Bioinform 2023, https://doi.org/10.1093/bib/bbad326).
  - Shachaf et al. show kNN-MI beats fixed binning for GRN inference (BMC Bioinformatics 2023, https://doi.org/10.1186/s12859-022-05047-5).
- **MIMIQ** (O'Hanlon et al., bioRxiv 2026, https://doi.org/10.64898/2026.04.10.717770). Adaptive binning built for small integer counts, with a negative-binomial copula χ² test. It is gene-gene only.
- **Other gene-gene MI:** SINUM (BMC Bioinformatics 2024, https://doi.org/10.1186/s12859-024-05895-3) and PIDC (Chan et al., Cell Systems 2017, https://doi.org/10.1016/j.cels.2017.08.014).
- **Zero-inflation-aware MI for cell-cell distance:** none found.

**Fast, GPU and approximate MI**
- **CUDA B-spline MI**, up to 82x (Shi et al., BMC Res Notes 2011, https://doi.org/10.1186/1756-0500-4-189).
- **GPU-ARACNE** (He et al., BMC Syst Biol 2017, https://doi.org/10.1186/s12918-017-0458-5).
- **TINGe**, parallel whole-genome MI (Zola et al., IEEE TPDS 2010).
- **GPU joint-histogram MI** for image registration, about 25x (Shams et al., Comput Methods Programs Biomed 2010). This is the same kernel as cell-cell MI (paired vectors into a joint histogram).
- **EDGE** (Noshad, Zeng & Hero, ICASSP 2019, https://arxiv.org/abs/1801.09125). LSH-based MI estimator in O(N) time for one pair of high-dimensional variables.
- **Sliced MI** (Goldfeld & Greenewald, NeurIPS 2021, https://arxiv.org/abs/2110.05279). It has been criticized as saturating and deceptive (Semenenko et al., arXiv 2025, https://arxiv.org/abs/2506.04053).
- **Fast-mRMR** (Ramírez-Gallego et al., Int J Intell Syst 2017) covers MI feature selection on CPU, GPU and Spark.
- All of these speed up gene-gene MI or single-pair estimation. None builds a kNN graph over millions of objects under an MI distance.

## 2. Scalable clustering at atlas scale

**Sketching and subsampling**
- **geosketch** (Hie et al., Cell Systems 2019, https://doi.org/10.1016/j.cels.2019.05.003). Plaid box covering.
- **Sphetcher** (Do et al., iScience 2020, https://doi.org/10.1016/j.isci.2020.101126).
- **Hopper** (DeMeo & Berger, Bioinformatics 2020, https://doi.org/10.1093/bioinformatics/btaa408). Farthest-first traversal, a 2-approximation of the Hausdorff distance. It found 64 macrophages (0.004%) in a 5,000-cell sketch of 1.3M cells.
- **Kernel herding** (Baskaran et al., ACM BCB 2022, https://arxiv.org/abs/2207.00584).
- **scSampler**, a maximin design (Song et al., Bioinformatics 2022, https://doi.org/10.1093/bioinformatics/btac271).
- **Seurat v5 leverage-score sketch plus ProjectData** (Hao et al., Nat Biotechnol 2024, https://doi.org/10.1038/s41587-023-01767-y). The vignette sketches 50K of 1.3M cells.
- **scValue** (Huang et al., Brief Bioinform 2025, https://doi.org/10.1093/bib/bbaf279). Tested on 20 datasets up to 4M cells against 7 samplers.
- **SKIM** (PLoS Comput Biol 2026, https://doi.org/10.1371/journal.pcbi.1014638).
- **Spatial benchmark** (Gingerich et al., NAR 2026, https://doi.org/10.1093/nar/gkag434). Leverage scores beat geosketch. No theory.
- **Crawford et al.** (bioRxiv 2023, https://doi.org/10.1101/2023.09.15.557825). Sketch choice mattered more than graph construction for label propagation.
- **Metacells:** MetaQ (Nat Commun 2025, https://doi.org/10.1038/s41467-025-56424-6) and CAMP, a coreset method claiming "provable geometric guarantees" (bioRxiv Dec 2025, https://doi.org/10.64898/2025.12.11.693725). The CAMP guarantee details are UNVERIFIED.

**Landmark, Nyström and anchor graphs**
- **Nyström spectral grouping:** Fowlkes et al., TPAMI 2004.
- **LSC:** Chen & Cai, AAAI 2011.
- **U-SPEC:** Huang et al., IEEE TKDE 2020, https://arxiv.org/abs/1903.01057. Handles about 10M points on a 64 GB PC.
- **Specter:** Do et al., Genome Res 2021, https://doi.org/10.1101/gr.267906.120.
- **Secuer:** PLoS Comput Biol 2022, https://doi.org/10.1371/journal.pcbi.1010753.
- **SnapATAC Nyström landmarks:** Nat Commun 2021, https://doi.org/10.1038/s41467-021-21583-9.

**Approximate kNN**
- **HNSW:** Malkov & Yashunin, TPAMI 2020, https://doi.org/10.1109/TPAMI.2018.2889473.
- **NN-Descent** for generic similarities: Dong et al., WWW 2011. pynndescent accepts numba custom metrics, which is what MICA uses.
- **FAISS:** Douze et al., arXiv 2024, https://arxiv.org/abs/2401.08281.
- **CAGRA** on GPU: Ootomo et al., ICDE 2024, https://arxiv.org/abs/2308.15136.
- **Limitation:** FAISS and CAGRA support only L2 and inner-product style metrics, not arbitrary MI.

**GPU stacks**
- **rapids-singlecell** (Dicks et al., arXiv 2026, https://arxiv.org/abs/2603.02402).
  - 1.3M cells in about 25 s on an RTX PRO 6000, versus about 52 min on a 32-core CPU.
  - Leiden about 100x faster.
  - About 1M cells per GPU before hitting the int32 limit. Dask out-of-core processed Tahoe-100M in under 20 min.
- **ScaleSC** (Bioinform Adv 2025, https://doi.org/10.1093/bioadv/vbaf167). 13M simulated cells in about 51 min using about 60 GB GPU memory.
- **Billato et al.** (bioRxiv 2025, https://doi.org/10.1101/2025.10.28.681564). GPU about 15x faster than the best CPU pipeline.

**Graph clustering and k-means family**
- **Leiden:** Traag et al., Sci Rep 2019, https://doi.org/10.1038/s41598-019-41695-z.
- **GVE-Leiden** (ICPP 2024, https://arxiv.org/abs/2312.13936). 403M edges/s on a 3.8B-edge graph, 3x faster than cuGraph on an A100.
- **Modularity resolution limit:** Fortunato & Barthélemy, PNAS 2007, https://doi.org/10.1073/pnas.0605965104. Resolution-limit-free CPM: Traag et al., PRE 2011, https://doi.org/10.1103/PhysRevE.84.016114. This is directly relevant to whether rare clusters survive in huge graphs.
- **PARC:** 1.1M cells in 13 min (Bioinformatics 2020, https://doi.org/10.1093/bioinformatics/btaa042).
- **mbkmeans:** 1.3M cells, HDF5 on-disk (PLoS Comput Biol 2021, https://doi.org/10.1371/journal.pcbi.1008625).
- **SC3s:** up to 2M cells (BMC Bioinformatics 2022, https://doi.org/10.1186/s12859-022-05085-z).
- **SHARP:** 10M cells via random projection ensembles (Genome Res 2020, https://doi.org/10.1101/gr.254557.119).
- **SCALEX:** online projection into a fixed embedding (Nat Commun 2022, https://doi.org/10.1038/s41467-022-33758-z).

**How atlases actually cluster**
- **Allen whole mouse brain** (Yao et al., Nature 2023, https://doi.org/10.1038/s41586-023-06812-z). This is explicitly cluster-then-propagate.
  - About 4M cells were clustered with scrattch.bigcat iterative Jaccard-Leiden.
  - "a randomly subsampled set of cells ... to compute high variance genes and perform PCA, then projected to all the cells".
  - For integration they down-sampled by within-modality clusters "ensuring preservation of rare cell types". Leiden ran on the down-sample, and "cluster membership of other cells were imputed based on KNNs" (Annoy).
  - The result was 5,322 clusters. UMAP used at most 1,000 cells per cluster, with the rest imputed.
- **CELLxGENE** (NAR 2025, https://doi.org/10.1093/nar/gkae1142).
  - The Census holds 65M+ cells, and the platform "does not recluster". Labels come from submitters.
  - The Census docs list maintained scVI (trained on the 74.3M-cell 2024-07-01 release) and Geneformer embeddings.
- **Tabula Sapiens 2.0** (Cell 2026, https://doi.org/10.1016/j.cell.2026.08.010). Its clustering pipeline is UNVERIFIED.

## 3. Prior work scaling MI distances to millions of objects

I found nothing that computes a pairwise MI distance, or an MI kNN graph, over more than 1M objects in any field. The closest items:

- **MICA-GE at 650K cells.** It uses NN-Descent with an exact per-pair MI, reports runtime only, and has no accuracy check at that size.
- **minicore at 4M cells,** with divergences rather than MI.
- **LSH for f-divergences and "mutual information loss"** (Chen et al., NeurIPS 2019, https://arxiv.org/abs/1910.12414). It hashes distributions, not the MI between two profiles.
- **EDGE,** which is one MI estimate over many samples.
- **Learned embeddings where Euclidean distance mimics an expensive metric** (Courty et al., "Learning Wasserstein Embeddings", ICLR 2018, https://arxiv.org/abs/1710.07457). This is a template for distilling MI distance into a FAISS-searchable space. I found no MI version.
- **No random-projection, sketch or landmark (Nyström) approximation of a cell-cell MI matrix** was found.

## 4. Cluster-then-propagate and rare-population guarantees

**In practice:** Seurat v5 ProjectData (`refdata = list(cluster_full = "seurat_clusters")`), the ABC atlas kNN imputation, and SnapATAC Nyström extension. None of them reports label accuracy on the propagated cells for rare types. The Seurat vignette does not quantify it.

**Theory that exists, all generic:**
- **Uniform sampling.** P(miss a population with fraction p) ≤ exp(−ps). CURE (Guha et al., SIGMOD 1998, https://doi.org/10.1145/276304.276312) gives a Chernoff sample-size bound so that each cluster u gets at least f|u| points with probability 1−δ (Theorem 4.1). The formula is paraphrased from a poorly OCRed PDF, so check it against the original.
- **Farthest-first traversal** is a 2-approximation of the k-center/Hausdorff objective (as used in Hopper). It works in any metric space, including MI space.
- **Kernel herding** converges in MMD (Chen, Welling & Smola, UAI 2010, https://arxiv.org/abs/1203.3472). It preserves proportions, so it does not enrich rare types.
- **Sensitivity-sampling coresets** (e.g. https://arxiv.org/abs/2405.01339) bound k-means cost, not recovery of a specific small cluster.
- **Nyström spectral clustering** is said to lack a bounded approximation ratio (Entropy 2018, https://doi.org/10.3390/e20070519). I did not verify this against the original proof.
- **Leverage scores** favoring rare populations is a heuristic claim (Seurat vignette), not a theorem.

I found no bound connecting sketch size, rare-type fraction, metric and kNN propagation error for single-cell data.

## 5. Gap analysis

**G1. No accuracy evidence for MI clustering above about 10K labeled cells.** scMINER itself reports only efficiency at 650K. Skinnider 2019 found MI weak. So "MI beats Euclidean on PCA at atlas scale, especially for rare types" is untested and could well come out negative. Evaluation targets with deep hierarchical labels already exist: the ABC atlas (about 4M cells) and Cao 2M (used by minicore).

**G2. No sketching or landmarks in MI space.** Every sketcher works in PCA or Euclidean space. Because MICA's D is a proven metric, farthest-first sketching under D inherits Hopper's 2-approximation Hausdorff guarantee in MI geometry. Landmark lower bounds |D(x,l) − D(y,l)| ≤ D(x,y) are also valid for pruning kNN candidates. Neither use was found.

**G3. No GPU cell-cell MI.** GPU MI exists only for gene-gene networks (2011, 2017) and image registration. Two unused ideas:
- A batched joint histogram, using bincount on b·bin_i + bin_j, or a one-hot contingency table computed by matmul.
- Sparsity. The (0,0) cell of the table equals G minus the number of genes with bin > 0 in either cell, so a pair only needs its nonzero high-bin genes.

This is feasible on 2x RTX 6000 Ada, but it is engineering. On its own it is probably not enough for a paper.

**G4. The cell-cell MI estimator has never been examined.**
- MICA uses per-cell min-max equal-width binning with about 20 bins for ~8,000 genes (b = G^(1/3)). With 80-95% zeros, bin 0 dominates and MI mostly reflects co-detection, which likely depends on sequencing depth (UNVERIFIED hypothesis).
- The following have not been compared for cell-cell distances: NVI vs NID (1 − I/max H), adaptive partitioning, mixed KSG, and count-aware binning (MIMIQ-style).
- scMINER only varied bin count (5 to 100).

**G5. No theory for sketch-then-propagate on rare types under a non-Euclidean metric.** CURE, k-center and herding bounds have not been tied to kNN label-propagation error on rare types. A bound plus empirical check (population separation greater than 2× covering radius implies a representative is sampled) looks doable in a semester.

**G6. Reproducibility.** The paper says HNSW but the code uses NN-Descent. Dependencies are pinned to 2021, and there are no accuracy benchmarks at large scale.

**Cost (rough, UNVERIFIED).** NN-Descent with k=80 means roughly 1e10 to 1e11 MI evaluations for 2M cells, each over thousands of genes. Running MICA-GE unmodified at 2M cells is likely many CPU hours to days. Time it on 100K cells first.

**Candidate contribution, stated skeptically.** Farthest-first sketch in MI metric (with guarantee), exact GPU MI on the sketch, clustering, then propagation with landmark pruning, plus the rare-type bound (G2 + G3 + G5). Its value depends entirely on G1 showing MI helps. A careful negative result would still be publishable as a benchmark.

## Search log

**WebSearch queries** (the session's 200-search budget ran out mid-sweep, so the rest used the Europe PMC and arXiv APIs):
- MICA mutual information clustering bioRxiv
- Skinnider 2019 proportionality
- Kim 2019 similarity metrics
- benchmark distance metrics single-cell 2023 2024
- proportionality rho propr
- KSG GPU implementation
- PIDC
- GPU MI ARACNe 2023
- MI estimation dropout sparse 2023
- "mutual information" "cell-cell" 2024/2025
- adaptive partitioning ARACNe-AP
- Leiden GPU cugraph rapids-singlecell
- rapids-singlecell 2026
- geosketch
- Hopper
- kernel herding sketching
- scSampler
- sketchR benchmark rare cell
- Seurat v5 sketch leverage
- Specter anchor graph
- U-SPEC
- SHARP
- Allen whole mouse brain clustering
- CELLxGENE Census 2025
- Tabula Sapiens 2025
- Kraskov 2005 hierarchical MI
- LSH Jensen-Shannon f-divergence
- sliced MI
- MI metric index ANN
- Meila VI
- MIC RapidMic
- GVE-Leiden
- mbkmeans
- SC3s
- SCALEX
- SnapATAC Nyström
- coreset rare clusters
- subsample-then-classify rare cells
- CURE Chernoff
- sketch vs full clustering ARI
- Nyström guarantee
- LSC
- pynndescent custom metric
- information-theoretic cell distance 2025/2026
- scMINER follow-up GPU
- SINUM
- Vinh 2010 metric
- Noshad EDGE
- Wasserstein embeddings
- TINGe
- fast-mRMR
- GPU MI image registration
- MI pairwise million samples hashing
- scTenifold/MI graph million cells
- FAISS/CAGRA
- HNSW
- ANN benchmark single-cell
- metacells millions
- MI clustering atlas-scale
- Hopper k-center guarantee

**Europe PMC queries:**
- "mutual information" AND "zero-inflated" AND single-cell (7 hits, none cell-cell)
- TITLE "mutual information" AND single-cell, 2022 to 2026 (16 hits)
- MI estimator AND single-cell (11 hits)
- MI AND sketching/subsampling AND clustering AND single-cell (0 hits)
- MI AND landmark/Nystrom/anchor graph AND clustering (2 hits, irrelevant)
- MI AND "cell-cell distance" (0 hits)
- MI AND GPU AND single-cell (0 hits)
- sketch AND rare cell AND guarantee (0 hits)

**arXiv API queries:**
- MI AND locality sensitive hashing (4 hits)
- MI AND nearest neighbor search (0 hits)
- variation of information AND nearest neighbor (0 hits)
- single-cell AND MI AND clustering (3 hits, irrelevant)

**Code and full texts read directly:** MICA `lib/distance.py`, `lib/neighbor_graph.py`, `mica_ge.py`, `mica.py`, `requirements.txt`, and the full texts of scMINER, Skinnider 2019, Kraskov 2005, Vinh 2010, the ABC atlas and rapids-singlecell.

**Not covered because the budget ran out:** Scribe/RDI, how PIDC discretizes, and the full text of the Yin 2025 benchmark.

Local text copies of the key PDFs (scminer.txt, skinnider.txt, vinh.txt, rsc.txt) are in `/tmp/claude-1599031162/-home-ugdf8-coursework-26Fall-Clustering-Algorithms/03c4a410-85dd-436d-8b00-c479892460b0/scratchpad/`.