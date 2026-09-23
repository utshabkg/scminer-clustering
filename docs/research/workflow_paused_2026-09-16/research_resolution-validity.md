**Literature sweep: choosing cluster resolution, judging cluster validity and evaluation pitfalls in single-cell data (2017 to Sept 2026, weighted to 2023 to 2026)**

**Three findings that matter most for the project**
1. Several published, peer-reviewed methods now pick the number of clusters with a statistical test. None of them has been applied to mutual-information (MI) based clustering such as scMINER's MICA.
2. The LLM and knowledge-guided ways of choosing resolution are all recent (late 2025 to 2026), mostly preprints, and none controls for double dipping. The marker genes an LLM reads are computed from the same cells that defined the clusters.
3. No independent benchmark I found includes MICA or scMINER. scMINER's own evaluation has known weak points (details in section 6).

---

### 1. Choosing resolution with a statistical test
- **CHOIR** (Sant, Mucke, Corces; *Nature Genetics* 2025, 57:1309–1319; https://doi.org/10.1038/s41588-025-02148-8).
  - It first builds a deliberately over-split tree, raising resolution until silhouette is maximized.
  - It then works back up the tree. For each pair of sibling clusters, it trains balanced random forests to tell them apart and compares their accuracy with a permutation null (100 permutations by default). Pairs that cannot be told apart are merged.
  - It was benchmarked on 230 simulated and 5 real datasets, the largest a 3.29M-nucleus human brain atlas, which shows it scales to our Dataset 2 size.
  - Full-text extraction lists competitors including Seurat/Scanpy Louvain and Leiden, SCCAF, sc-SHC, Cytocipher, GiniClust3, CIDR, SAFE-clustering, scCAN, SHARP and Spectrum. The exact list of 15 is UNVERIFIED.
  - Reported results: Seurat over-split homogeneous snATAC populations and under-split a pool of 190 cell lines. sc-SHC avoided over-splitting but tended to under-split.
  - Stated limitations: it costs more than a single Louvain run, and it assumes batch is not confounded with biology.
- **sc-SHC** (Grabski, Street, Irizarry; *Nature Methods* 2023, 20:1196–1202; https://doi.org/10.1038/s41592-023-01933-9).
  - It extends Significance of Hierarchical Clustering (SHC) to counts, testing each split against a fitted single-population null with family-wise error control.
  - `testClusters()` can audit any existing labeling, including MICA's.
  - Defaults per the README: alpha 0.05 (0.25 suggested for discovery), 2,500 features, 30 PCs (https://github.com/igrabski/sc-SHC).
- **recall** (DenAdel et al.; *American Journal of Human Genetics* 2025; https://www.sciencedirect.com/science/article/pii/S0002929725000618; code at https://github.com/lcrawlab/recall).
  - It adds knockoff-style artificial null genes, clusters, and requires every cluster pair to show significant differential expression after calibration. If a pair fails, it re-clusters at lower resolution.
  - Discrepancy: the docs give the title "A knockoff calibration method to avoid over-clustering…" and DOI 10.1016/j.ajhg.2025.02.008. The journal title is "Artificial variables help to avoid over-clustering…" and Nature Index shows 10.1016/j.ajhg.2025.02.014. Exact DOI UNVERIFIED.
- **Cytocipher** (*Bioinformatics* 2023, 39(7):btad435; https://academic.oup.com/bioinformatics/article/39/7/btad435/7224247). Starts over-clustered and merges clusters whose marker-set enrichment scores do not differ significantly.
- **SCCAF** (Miao et al.; *Nature Methods* 2020, 17:621–628; https://doi.org/10.1038/s41592-020-0825-9). Merges clusters based on the confusion matrix of a classifier trained on cluster labels. CHOIR is a statistically calibrated descendant of this idea.

### 2. Stability, consensus and multi-resolution methods
- **SC3** (Kiselev et al.; *Nature Methods* 2017; https://doi.org/10.1038/nmeth.4236). Consensus k-means, with k estimated from eigenvalues that exceed the Tracy–Widom bound. **SC3s** extends it to millions of cells (*BMC Bioinformatics* 2022; https://doi.org/10.1186/s12859-022-05085-z).
- **chooseR** (Patterson-Cross, Levine, Menon; *BMC Bioinformatics* 2021, 22:39; https://doi.org/10.1186/s12859-021-03957-4). Subsamples 80% of cells, re-clusters, picks the resolution with near-maximal median co-clustering silhouette, and reports per-cluster robustness.
- **MultiK** (Liu, Thennavan, Garay, Marron, Perou; *Genome Biology* 2021; https://doi.org/10.1186/s13059-021-02445-5). Subsampling consensus across resolutions, returning several candidate K values.
- **scCCESS** (*BMC Bioinformatics* 2019; https://doi.org/10.1186/s12859-019-3179-5; https://github.com/PYangLab/scCCESS). Ensembles of autoencoders on random projections. K is the most stable choice across the ensemble.
- **IKAP** (Chen et al.; *GigaScience* 2019, 8(10):giz121; https://academic.oup.com/gigascience/article/8/10/giz121/5579995). Tunes PCs and resolution using gap statistics plus how well differentially expressed markers classify the clusters.
- **clustree** (Zappia and Oshlack; *GigaScience* 2018, 7(7):giy083; https://github.com/lazappi/clustree). A visualization of how clusters split across resolutions, not a selection rule.
- **MRtree** (Peng, …, Roeder; *Nucleic Acids Research* 2021, 49(16):e91; https://academic.oup.com/nar/article/49/16/e91/6300355). Reconciles clusterings at many resolutions into one hierarchy and cuts the tree where stability drops sharply.
- **scclusteval** (Tang et al.; *Bioinformatics* 2021, 37(15):2212; https://academic.oup.com/bioinformatics/article/37/15/2212/5962080). Per-cluster Jaccard stability under subsampling.
- **scICE** (Kim et al.; *Nature Communications* 2025; https://doi.org/10.1038/s41467-025-60702-8).
  - Uses an "inconsistency coefficient" computed across Leiden runs with different random seeds, so no consensus matrix is needed. A result counts as consistent when median IC < 1.005.
  - Claims up to 30x speedup over multiK and chooseR, tested on 48 datasets.
  - A press release (https://phys.org/news/2025-07-scientists-tackle-cell-reliability-crisis.html) says about two-thirds of existing analyses were unstable. That figure is the press release's claim, not verified in the paper.
- **Stability is not validity.** A stable split can still be driven by noise, batch or cell-cycle signal. CHOIR, sc-SHC and recall exist because stability alone lets in false splits.

### 3. Double dipping and inference after clustering
- **Count splitting** (Neufeld, Gao, Popp, Battle, Witten; *Biostatistics* 2024, 25(1):270–287; https://academic.oup.com/biostatistics/article-abstract/25/1/270/6893953).
  - Binomial thinning splits each count into independent training and test matrices (Poisson assumption, epsilon = 0.5 recommended).
  - You cluster on the training matrix and test on the test matrix.
- **Negative binomial count splitting** (Neufeld et al.; arXiv 2023; https://arxiv.org/abs/2307.12985). Needed because Poisson splitting breaks down under overdispersion. Journal status UNVERIFIED.
- **Data thinning** (Neufeld, Dharamshi, Gao, Witten; *JMLR* 2024, 25(57); https://arxiv.org/abs/2301.07276). Generalizes the idea to any convolution-closed distribution, and the `datathin` package includes a tutorial on choosing K.
- **Selective inference, exact p-values**
  - Hierarchical clustering: Gao, Bien, Witten, *JASA* 2024, 119(545):332–342 (https://doi.org/10.1080/01621459.2022.2116331).
  - k-means: Chen and Witten, *JMLR* 2023, 24(152) (https://www.jmlr.org/papers/v24/22-0371.html).
  - Single-feature tests (CADET): Chen and Gao, *Biostatistics* 2025, 26(1):kxae046 (https://arxiv.org/abs/2311.16375).
  - False discovery rate control via data splitting for testing after clustering: arXiv 2024 (https://arxiv.org/abs/2410.06451).
- **ClusterDE** (bioRxiv 2023, https://doi.org/10.1101/2023.07.21.550107; software paper in *Journal of Computational Biology* 2026, https://doi.org/10.1177/15578666251383562). A synthetic one-cluster null serves as a negative control for post-clustering differential expression. A search snippet says the full method appeared at RECOMB 2025 (UNVERIFIED).
- **Relevance to scMINER.** Its differential driver-activity tests between MICA clusters use the same cells that defined the clusters, which is exactly the double-dipping pattern above. This is my inference; the paper does not discuss it.

### 4. Silhouette and internal metrics
- **Rautenstrauch and Ohler**, "Shortcomings of silhouette in single-cell integration benchmarking" (*Nature Biotechnology* 2025; https://doi.org/10.1038/s41587-025-02743-4; preprint https://www.biorxiv.org/content/10.1101/2025.01.21.634098v1).
  - Silhouette assumes compact, convex clusters, which label-defined single-cell groups often are not.
  - It cannot reliably assess batch removal or biology conservation. They recommend CiLISI plus ARI-type metrics instead.
- **Yerbury et al.** (*ACM TKDD* 19(8), 2025; https://arxiv.org/abs/2404.10351). Based on 2.7M partitions, relative validity indices such as silhouette are unsuitable for choosing between similarity paradigms. That is exactly what "MI distance vs PCA-Euclidean judged by silhouette" does.
- **Tomašev and Radovanović**, "Clustering evaluation in high-dimensional data" (Springer, 2016; https://link.springer.com/chapter/10.1007/978-3-319-24211-8_4). Hub points sit near cluster centres and shrink between-cluster distances.
- **Counterpoint:** Powell (arXiv 2022; https://arxiv.org/abs/2201.05214) finds most Euclidean validity indices stay sensitive in high dimensions on synthetic data.
- **Sciaraffa et al.** (*Frontiers in Bioinformatics* 2025, 5:1562410; https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2025.1562410/full). Of 15 internal indices, only Banfield–Raftery and within-cluster dispersion predicted accuracy across datasets. Silhouette did worse.
- **Validity-index benchmark:** Simpson, Campello, Stojanovski (arXiv 2025; https://arxiv.org/abs/2511.05983) compare 26 internal indices on more than 16,000 datasets.

### 5. External metrics, rare cell types and hierarchy
- **ARI vs AMI.** Romano, Vinh, Bailey, Verspoor, *JMLR* 2016, 17 (https://arxiv.org/abs/1512.01286): prefer ARI when reference clusters are large and balanced, and AMI when the reference is unbalanced with small clusters.
- **Pair-counting decomposition.** Warrens and van der Hoef, *Journal of Classification* 2022, 39:487–509 (https://doi.org/10.1007/s00357-022-09413-z): ARI breaks down into per-cluster agreement dominated by large clusters, so rare types barely move it.
- **Choice of null model** also changes chance-adjusted scores: Gates and Ahn (arXiv 2017; https://arxiv.org/abs/1701.06508).
- **Rare cell types**
  - scCAD (*Nature Communications* 2024; https://doi.org/10.1038/s41467-024-51891-9) scores rare-type detection with F1 across 25 datasets. Best F1 was only about 0.42, which shows how hard this is.
  - CHOIR adds an entropy-of-accuracy measure to detect under-clustering.
- **Hierarchy-aware metrics**
  - wRI and wNMI weight disagreements by cell-type tree distance: Wu and Wu, *Genome Biology* 2020, 21:123 (https://doi.org/10.1186/s13059-020-02027-x).
  - LCAD (lowest-common-ancestor distance on the Cell Ontology) and scGraph-OntoRWR: Wu et al., *Genome Biology* 2025, 26:334 (https://doi.org/10.1186/s13059-025-03781-6).
  - Hop-F1: SOCAM preprint 2026 (https://www.biorxiv.org/content/10.64898/2026.01.13.699356v1.full).
  - Cell Ontology reference: Tan et al. (https://arxiv.org/abs/2506.10037). Search says it was published in *Scientific Data* 2026 (UNVERIFIED).

### 6. Knowledge-guided and LLM-guided resolution selection
- **GPTAnno** (Song, Tang, Liu, Wang, Qian, Zou, Hou; bioRxiv Dec 2025; https://doi.org/10.1101/2025.11.27.690951). This is the closest prior art to our RAG idea.
  - Runs Louvain at resolutions 0.1 to 0.7 and queries GPT-5 ten times per resolution with the top-10 Wilcoxon markers per cluster.
  - Maps the answers to Cell Ontology terms and scores each resolution by consistency, robustness and mean ontology distance. It picks the highest composite score.
  - Evaluated on 12 datasets with a weighted Cell Ontology agreement score. The chosen resolution was best in 8 of 12 and within 10% otherwise.
  - Gaps: no statistical null, no marker splitting, a coarse resolution grid, and a closed-source LLM.
- **HypoGeneAgent** (Yuan et al., Genentech; arXiv Sept 2025; https://arxiv.org/abs/2509.09740). The LLM proposes GO hypotheses for gene-program clusters. Resolution is scored by how similar hypotheses are within a cluster versus between clusters, beating silhouette and modularity on one K562 Perturb-seq dataset.
- **scMRMA** (*Nucleic Acids Research* 2022, 50(2):e7; https://academic.oup.com/nar/article/50/2/e7/6396893). Uses the hierarchical PanglaoDB marker database to drive iterative re-clustering. Authors not re-verified.
- **LLM annotators with split or merge ability, none of which selects resolution in a principled way**
  - GPTCelltype: Hou and Ji, *Nature Methods* 2024 (https://doi.org/10.1038/s41592-024-02235-4).
  - CASSIA (includes a Subclustering agent): Xie et al., *Nature Communications* 2025 (https://doi.org/10.1038/s41467-025-67084-x). Volume and year UNVERIFIED.
  - LICT (checks marker expression to judge credibility): *Communications Biology* 2025 (https://doi.org/10.1038/s42003-025-08745-x).
  - CellTypeAgent (verifies against CELLxGENE): arXiv 2025 (https://arxiv.org/abs/2505.08844).
  - mLLMCelltype (consensus across several LLMs): *Communications Biology* 2026 (https://doi.org/10.1038/s42003-026-10420-8).
  - scExtract (reads the source article to annotate; I found no evidence it sets resolution): *Genome Biology* 2025, 26:174 (https://doi.org/10.1186/s13059-025-03639-x).
- **scLLM-DSC** (arXiv June 2026; https://arxiv.org/abs/2606.13007) injects LLM and NCBI gene embeddings into deep clustering. It does not choose K.

### 7. Clustering benchmarks, and where MI-based methods stand
- **Yu, Cao, Yang, Yang** (*Genome Biology* 2022; https://doi.org/10.1186/s13059-022-02622-0). Tested 14 methods on estimating the number of cell types using Tabula Muris subsamples. Monocle3, scLCA, scCCESS-SIMLR, CIDR and Seurat had the best concordance.
- ***Briefings in Bioinformatics* 2024, 25(1):bbad497** (https://academic.oup.com/bib/article/25/1/bbad497/7505072). Seven methods: CosTaL, Seurat, DESC and SC3 beat SCCAF and scDeepCluster.
- ***Genome Biology* 2025** (https://doi.org/10.1186/s13059-025-03719-y). 28 methods on 10 paired transcriptomic and proteomic datasets. Top three for transcriptomics were scDCC, scAIDE and FlowSOM. Authors not verified.
- **scCluBench** (arXiv Dec 2025, https://arxiv.org/abs/2512.02471; listed in the AAAI-26 proceedings, https://ojs.aaai.org/index.php/AAAI/article/view/37110).
  - 15 methods on 36 datasets, run with the true K given.
  - Single-cell foundation models underperformed. Top methods were scCDCG, scMAE and scNAME.
  - scCDCG may come from the same group (UNVERIFIED).
- **arXiv July 2026** (https://arxiv.org/html/2607.25288v1). Deep autoencoders were not significantly better than PCA+Leiden or scVI.
- **MI-based methods:** none of the benchmarks above includes MICA or scMINER, and I found no independent benchmark that does. In scMINER's own evaluation:
  - The 2023 preprint (https://www.biorxiv.org/content/10.1101/2023.01.26.523391) compared against Seurat, Scanpy and SC3 on 11 datasets using ARI.
  - The published *Nature Communications* 2025, 16:4305 version reportedly used 10 datasets against Seurat, Scanpy, SC3s, scVI and scDeepCluster, adding AMI, NMI, ASW and AvgBIO. This comes from search snippets because the full text was blocked; only partly verified.
  - In the preprint, resolutions for MICA and Seurat were tuned so the number of clusters matched the known number of cell types.
  - Silhouette of 2D UMAP against reference labels was offered as evidence of purity.
  - K is chosen by ranking silhouette plus the user's prior knowledge. The authors call choosing K "a manual and problem-dependent task." I have not checked whether the published version changed this.

---

### GAP ANALYSIS

**Open gaps that look defensible for a publication**
1. **No significance-calibrated resolution selection for MI-distance clustering.** CHOIR, sc-SHC and recall all assume PCA/Euclidean pipelines. A feasible contribution: run CHOIR-style permutation tests or recall-style knockoff genes on MICA's graph mode, then compare against silhouette-based K. This runs on CPU and fits the 120-core machine.
2. **LLM and knowledge-guided resolution choice has a double-dipping problem nobody addresses.** GPTAnno, CASSIA and LICT all feed markers computed on the clusters being judged. Proposed approach:
   - Cluster on the count-split training fold.
   - Derive markers on the test fold.
   - Retrieve evidence from CellMarker, PanglaoDB, the Cell Ontology and literature (RAG).
   - Accept a split only when both the statistical test and the retrieved evidence support it.
   - Evaluate with LCAD or hierarchy-weighted metrics, rare-type F1 and AMI.
   
   This fits Utshab's retrieval and vLLM background, and a local open-weight LLM avoids API cost. I found nothing combining count splitting with LLM annotation.
3. **Evaluation protocol gap.** Benchmarks usually give methods the true K (scCluBench, scMINER). No study evaluates methods at their own automatically selected K using hierarchy-aware plus rare-type metrics. That would be a solid course deliverable and a possible workshop or short paper.
4. **Scale.** Pairwise MI costs O(n^2). Dataset 2 (2M cells) needs graph mode or sketching. CHOIR's 3.3M-cell run is a useful precedent.
5. **Risk.** GPTAnno already occupies "LLM picks resolution via ontology distance." Novelty has to come from statistical calibration, retrieval grounding and MI clustering, not from LLM resolution selection itself.

**Not found or unverified**
- No paper with the exact "clustering deserves more attention" title. The closest critiques are CHOIR, sc-SHC, and "From Cluster to Claim" (*Advanced Genetics* 2026, https://doi.org/10.1002/ggn2.70045). The last was blocked (403) and is known only from a search snippet.
- Full text of the published scMINER paper (PMC and Nature paywall/captcha, Europe PMC returned 503).
- Author lists for the *Genome Biology* 2025 benchmark and scMRMA.
- Publication status of NB count splitting.
- Exact volume numbers for CASSIA and *Nature Biotechnology* silhouette paper.
- The WebSearch budget (200 calls) ran out before I could check scCDCG authorship or atlas-scale benchmarks with rare-type metrics.

**Queries tried** (WebSearch unless noted):
- CHOIR Nature Genetics 2025
- Grabski scSHC Nature Methods 2023
- count splitting Neufeld Witten double dipping
- chooseR BMC Bioinformatics
- MultiK Genome Biology
- scCCESS autoencoder ensemble
- clustree GigaScience
- IKAP GigaScience
- "Benchmarking clustering algorithms on estimating the number of cell types"
- SC3 random matrix theory
- "clustering deserves more attention" critique (two variants)
- scICE inconsistency coefficient
- silhouette shortcomings Rautenstrauch
- ARI imbalanced / AMI Vinh
- Romano JMLR 2016
- Warrens van der Hoef
- Cell Ontology distance metric
- LCAD / scGraph-OntoRWR
- Wu wRI wNMI
- LLM merge/split clusters
- LLM agent resolution selection (found HypoGeneAgent)
- gene set enrichment for resolution selection
- ClusterDE
- Cytocipher
- CASSIA
- data thinning
- selective inference for hierarchical and k-means clustering
- NB count splitting
- CADET single-feature test
- count splitting to choose K
- recall / artificial variables
- MICA and scMINER benchmarks
- large-scale clustering benchmarks 2024–2026
- Genome Biology 2025 transcriptomic/proteomic benchmark
- Briefings 2024 critical assessment
- rare cell type detection F1
- hierarchy-aware ARI
- high-dimensional silhouette and hubness
- SCCAF
- GPTCelltype
- scMRMA / ScType resolution
- mLLMCelltype / LICT / CellTypeAgent
- scExtract
- WebFetch: GPTAnno PDF, scMINER 2023 preprint PDF, HypoGeneAgent, scLLM-DSC, scCluBench, arXiv 2607.25288, Yerbury, Powell, Simpson, sc-SHC GitHub, recall docs, CHOIR site, MRtree (NAR), silhouette bioRxiv, phys.org (scICE)

One departure from read-only: a failed curl wrote a 190-byte HTML error page to the session scratchpad (`/tmp/claude-1599031162/-home-ugdf8-coursework-26Fall-Clustering-Algorithms/03c4a410-85dd-436d-8b00-c479892460b0/scratchpad/scminer.xml`). No project files were touched.