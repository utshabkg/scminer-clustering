# Public datasets for the three course regimes (scMINER project), checked 2026-09-16

**Picks:** Dataset 1: Levine_32dim CyTOF (backup: PBMC14K). Dataset 2: Allen whole mouse brain 10Xv3, 2.3M cells (backup: Cao 2020 fetal atlas). Dataset 3: Zeisel 2015 (backup: Kolodziejczyk 2015 or Chu 2016). Everything below the heading "Verified" was checked against a primary source this session unless it is marked UNVERIFIED. I ran out of web searches near the end, so a few counts come from memory and are marked.

## 0. What scMINER itself used (verified)

The source is Supplementary Table 1 of Pan et al. 2025, *Nature Communications* 16:4305 ([supplement PDF](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-59620-6/MediaObjects/41467_2025_59620_MOESM1_ESM.pdf); [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12062461/)).

| Set | Cells | Classes | Accession |
|---|---|---|---|
| Yan / Goolam / Buettner / Pollen / Chung | 124 / 124 / 182 / 249 / 317 | 8 / 5 / 3 / 11 / 5 | GSE36552, E-MTAB-3321, E-MTAB-2805, SRP041736, GSE75688 |
| Usoskin | 622 | 4 | GSE59739 |
| Kolodziejczyk | 704 | 3 | E-MTAB-2600 |
| Klein | 2,717 | 4 | GSE65525 |
| Zeisel | 2,903 | 7 | GSE60361 |
| PBMC14K (Zheng 2017, FACS-sorted) | 13,605 | 7 | SRP073767 |
| HMC76K (human motor cortex) | 76,533 | 20 | Azimuth |
| Covid97K | 97,499 | 28 | E-MTAB-10026 |
| Covid650K | 647,366 | 51 | EGAD00001007982 |

- **Covid accessions:** EGAD… is controlled access at EGA. E-MTAB-10026 is Stephenson 2021, which reports 781,123 cells passing QC ([Nat Med 2021](https://www.nature.com/articles/s41591-021-01329-2)). The two accessions in the table may be swapped (UNVERIFIED).
- **Size limits:** The paper says scMINER "requires a minimum of 500 cells". Its clustering component (MICA) uses consensus k-means with MDS below 5,000 cells and a graph-embedding mode above that. The largest dataset tested was about 650K cells. Nothing shows that scMINER runs at 2M cells (UNVERIFIED).
- **PBMC14K:** The data ships inside the scMINER R package (17,986 genes, 14,000 cells before QC, 7 × 2,000 sorted populations) ([tutorial](https://jyyulab.github.io/scMINER/articles/quick_tutorial.html)). Both [scMINER](https://github.com/jyyulab/scMINER) and [MICA](https://github.com/jyyulab/MICA) are Apache-2.0 and accept h5ad input.
- **Implication:** Of the scMINER sets, only Usoskin, Kolodziejczyk, Klein, Zeisel and PBMC14K meet both the course thresholds and scMINER's 500-cell floor.

## 1. How trustworthy the labels are

From most to least trustworthy:
1. Labels from genotype demultiplexing (Tahoe-100M, CellBench) or from the experimental design, such as culture condition (Kolodziejczyk).
2. FACS sorting (Zheng sorted PBMC).
3. Manual gating (cytometry).
4. Curated labels made by clustering the same data (Allen atlas, Zeisel, Baron, Parse). These are circular: they favor methods that work like the one used to make them.
5. Label transfer (HLCA extended, MERFISH mapped to the 10x taxonomy).

The report should state which tier each dataset's labels come from.

## 2. Dataset 1 (at least 5,000 samples and 20 variables, "ordinary multivariate")

| Candidate | Size | Labels | Source / format | License | Notes |
|---|---|---|---|---|---|
| **Levine_32dim** (Levine 2015 CyTOF bone marrow) | 265,627 cells × 32 proteins | 14 manually gated populations, about 39% of cells labeled (about 104K) | [HDCytoData](https://bioconductor.org/packages/HDCytoData) (SummarizedExperiment/flowSet), [FlowRepository FR-FCM-ZZPH](https://flowrepository.org/id/FR-FCM-ZZPH) | UNVERIFIED | Classic benchmark ([Weber & Robinson 2016, Cytometry A](https://doi.org/10.1002/cyto.a.23030); [Weber & Soneson 2019, F1000Research](https://pubmed.ncbi.nlm.nih.gov/31857895/)). A 2025 Genome Biology benchmark of 28 clustering algorithms ranked FlowSOM among the top methods ([doi](https://doi.org/10.1186/s13059-025-03719-y)). |
| **Samusik_01** (mouse bone marrow CyTOF) | 86,864 × 39 | 24 gated populations, 53,173 labeled | HDCytoData | UNVERIFIED | Same caveats as Levine_32dim. |
| **PBMC14K / Zheng sorted** | 13,605 × 17,986 (full sorted collection about 94K cells, UNVERIFIED) | 7 FACS populations (10 in the full set) | scMINER package; 10x website per-population matrices ([example](https://www.10xgenomics.com/datasets/cd-19-plus-b-cells-1-standard-1-1-0)); 20K-cell subset on [Zenodo 3357167](https://doi.org/10.5281/zenodo.3357167) ([Abdelaal 2019](https://doi.org/10.1186/s13059-019-1795-z)) | Other 10x dataset pages say CC BY 4.0; not confirmed for these 2016 pages | Paper: [Zheng 2017, Nat Commun](https://doi.org/10.1038/ncomms14049). |
| **Hao 2021 CITE-seq PBMC** (protein panel only) | about 162K cells × 228 antibodies | l1/l2/l3 hierarchy | GSE164378; Seurat/Azimuth reference | UNVERIFIED | [Cell 2021](https://doi.org/10.1016/j.cell.2021.04.048). scMINER used this dataset for CITE-seq validation. Labels came from WNN clustering, so they are partly circular. |
| Baron human pancreas | 8,569 × about 17.5–20K genes | 14 types, from clustering | GSE84133; Bioconductor [scRNAseq](https://bioconductor.org/packages/scRNAseq) | GEO public | Four donors, so donor batch effects. |

**Pitfalls**
- **Cytometry scaling:** Apply arcsinh first (cofactor 5 for CyTOF). Score only the gated cells, or report unassigned cells separately.
- **MI on few features:** scMINER estimates mutual information between cells across features. With 32 to 39 features that estimate is noisy. This is an interesting adaptation question, but expect weak results.
- **PBMC14K overlap:** Each sorted population was a separate library, so library effects line up with the labels and can make clustering look easier than it is. The CD4 memory, naive and Treg sorts also overlap biologically, so ARI cannot reach 1. scMINER's own figure shows mixed CD4 memory and Treg clusters.

## 3. Dataset 2 (at least 2,000,000 samples and 10 variables)

| Candidate | Size | Labels (tier) | Download | License | Notes |
|---|---|---|---|---|---|
| **Allen whole mouse brain, 10Xv3 only** ([Yao 2023, Nature](https://doi.org/10.1038/s41586-023-06812-z)) | 2.3M cells × 32,285 genes (4.04M across all chemistries) | 34 classes, 338 subclasses, 1,201 supertypes, 5,322 clusters (tier 4, heavily curated, checked against MERFISH) | Public AWS S3, no login. I listed the bucket: 13 raw h5ad files split by brain region, 3.1 to 13.7 GB each, about 95 GB raw, and the same again for log2 versions ([docs](https://alleninstitute.github.io/abc_atlas_access/descriptions/WMB-10Xv3.html)) | CC BY-NC 4.0 | One chemistry, so no 10Xv2/v3 batch effect. Labels are in a separate cell-metadata CSV joined on `cell_label`. You concatenate 13 files. Neurons dominate the class balance. |
| **Cao 2020 human fetal atlas** ([Science](https://doi.org/10.1126/science.aba7721)) | 4,062,980 cells × 45,676 genes | 68 CELLxGENE cell types (77 main types in the paper) (tier 4) | **One h5ad file, 20.03 GB**, on [CELLxGENE](https://cellxgene.cziscience.com/collections/c114c20f-1ef4-49a5-9c2e-d965787fb90c); a 1M-cell subset is also there | CC BY 4.0 ([CELLxGENE policy](https://cellxgene.cziscience.com/docs/032__Contribute%20and%20Publish%20Data)) | sci-RNA-seq3 gives very shallow, sparse cells. Cell type and organ are confounded. |
| **Parse 10M PBMC cytokine** | 9,697,974 cells, 12 donors, 90 cytokines | 18 types (Harmony + Scanpy clustering + manual review, tier 4) plus Allen Institute (AIFI) Human Immune Health Atlas and CellTypist labels ([AIFI](https://apps.allenimmunology.org/aifi/resources/parse-10m-cytokines/)) | One 41 GB h5ad behind a web form ([Parse](https://www.parsebiosciences.com/datasets/10-million-human-pbmcs-in-a-single-experiment/)) | CC BY-NC 4.0 | Fits the PBMC theme best. Cytokine responses add structure that competes with cell type. The labels favor Leiden-style methods. |
| **Tahoe-100M** ([bioRxiv 2025](https://doi.org/10.1101/2025.02.20.639398)) | 100,648,790 cells × 62,710 genes, 50 cancer cell lines, 14 plates of 8.2 to 14.5M cells each | Cell line assigned by demuxlet on SNPs, reported accuracy above 98% (tier 1) | Parquet on [Hugging Face](https://huggingface.co/datasets/tahoebio/Tahoe-100M) (1.69 TB, streamable); h5ad per plate from the [Arc Virtual Cell Atlas](https://github.com/ArcInstitute/arc-virtual-cell-atlas) (Google Cloud, requester pays) | CC0 | The labels are not circular, but cell lines are not cell types and separate too easily. Drug effects confound clusters. You must take a 2M-cell slice. |
| HLCA full ([Nat Med 2023](https://doi.org/10.1038/s41591-023-02327-2)) | 2,282,447 cells | Only the 584,944-cell core has expert labels; the rest were transferred with scArches (tier 5) | One 21.86 GB h5ad, [CELLxGENE](https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293) | CC BY 4.0 | 49 studies with disease, so strong batch effects. |
| Allen MERFISH brain C57BL6J-638850 | 4,334,174 cells × 550 genes (500-gene panel) | Mapped to the 10x taxonomy (tier 5) | 7 GB log2 h5ad ([docs](https://alleninstitute.github.io/abc_atlas_access/descriptions/MERFISH-C57BL6J-638850.html)) | CC BY-NC 4.0 | Lower-dimensional and fairly dense, so it is a cheap large-scale option. |
| MOCA, Cao 2019 ([Nature](https://doi.org/10.1038/s41586-019-0969-x)) | 2,058,652 cells profiled | 38 main types, 655 subclusters (tier 4) | [Portal](https://oncoscape.v3.sttrcancer.org/atlas.gs.washington.edu.mouse.rna/downloads), text and RDS files, no h5ad (UNVERIFIED) | UNVERIFIED | Median 671 UMIs per cell. Only about 1.3 to 1.5M cells pass QC or carry labels ([PMC6434952](https://pmc.ncbi.nlm.nih.gov/articles/PMC6434952/)), so **it falls below 2M labeled cells**. Avoid. |

**Ruled out**
- **Tabula Sapiens:** about 1.1M cells across tissue-split files ([CELLxGENE](https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5)), below 2M.
- **Siletti 2023 human brain** ([Science](https://doi.org/10.1126/science.add7046)): split into many overlapping CELLxGENE files, no single download.
- **CELLxGENE Census** (2025-11-08 LTS, more than 125M unique cells, [release notes](https://chanzuckerberg.github.io/cellxgene-census/cellxgene_census_docsite_data_release_info.html)): usable for slices, but author labels vary between studies and you must filter on `is_primary_data`.

## 4. Dataset 3 (at least 500 samples and 700 variables)

| Candidate | Cells × genes | Labels (tier) | Accession / download | Used by |
|---|---|---|---|---|
| **Zeisel 2015** ([Science](https://doi.org/10.1126/science.aaa1934)) | 3,005 × 19,972 (scMINER used 2,903) | 9 level-1 classes (7 in scMINER), 47 subclasses (tier 4, from BackSPIN clustering) | GSE60361; Bioconductor scRNAseq; [Hemberg lab](https://hemberg-lab.github.io/scRNA.seq.datasets/) RDS | scMINER, SC3 and many others |
| **Kolodziejczyk 2015** ([PubMed](https://pubmed.ncbi.nlm.nih.gov/26431182/)) | 704 × about 38K (UNVERIFIED) | 3 culture conditions: serum, 2i, a2i (tier 1) | E-MTAB-2600 | scMINER gold standard |
| **Chu 2016** ([Genome Biol](https://doi.org/10.1186/s13059-016-1033-x)) | 1,018 × 19,097 | 7 cell lines or differentiated types, 69 to 212 cells each (tier 1 or 2) | GSE75748 | Many clustering papers, e.g. [scCluBench](https://arxiv.org/abs/2512.02471) |
| Usoskin 2015 ([Nat Neurosci](https://doi.org/10.1038/nn.3881)) | 622 × about 25K (UNVERIFIED) | 4 types, 11 subtypes (tier 4) | GSE59739 | scMINER |
| Klein 2015 | 2,717 × about 24K (UNVERIFIED) | 4 time points (tier 1, but biologically graded) | GSE65525 | scMINER |
| CellBench 10x ([Tian 2019, Nat Methods](https://doi.org/10.1038/s41592-019-0425-8)) | 3,803 × 11,778 (Abdelaal filtering) | 5 cell lines, genotype truth (tier 1) | GSE118767; [sc_mixology](https://github.com/LuyiTian/sc_mixology) | Abdelaal 2019, Szmigiel 2025 ([bioRxiv](https://doi.org/10.1101/2025.08.20.671270)) |
| Baron, Muraro, Segerstolpe, Xin | 8,569 / 2,122 / 2,133 / 1,449 | Clustering-based (tier 4) | GSE84133, GSE85241, E-MTAB-5061, GSE81608; [Zenodo 3357167](https://doi.org/10.5281/zenodo.3357167) | Abdelaal 2019 |

- **Too small:** Darmanis (466 cells), Pollen, Yan and Buettner all fall below 500 cells.
- **Pitfalls:**
  - Kolodziejczyk and Usoskin include ERCC spike-ins and use read counts from full-length protocols. Strip the ERCCs and do not use UMI-based normalization.
  - scMINER preprocessed Buettner with FPKM and Chung with TPM, following SC3.
  - Kolodziejczyk is easy: most methods reach an ARI near 1, so it does not separate methods well.
- **Ready-made h5ad collection:** scCluBench (arXiv Dec 2025; AAAI 2026, [proceedings](https://ojs.aaai.org/index.php/AAAI/article/view/37110)) has 36 h5ad datasets of 611 to 27,499 cells. It benchmarks scGPT, Geneformer and GeneCompass, which is a useful 2025–26 comparison point. Hosting link not verified.

## 5. Memory: 2M cells × 2,000 HVGs

- **Dense float32:** 2e6 × 2e3 × 4 B = **16.0 GB** (14.9 GiB); float64 is 32 GB.
  - Fits easily in 754 GB RAM.
  - Fits on one 48 GB GPU with about 30 GB left for work.
- **Sparse CSR** with float32 values and int32 indices costs 8 B per nonzero, plus 16 MB for row pointers:

  | Density | Size |
  |---|---|
  | 2% (shallow sci-RNA-seq) | 0.64 GB |
  | 10% | 3.2 GB |
  | 20% | 6.4 GB |
  | 50% | 16 GB, same as dense |

  int32 indices work up to 2^31 nonzeros, which is about 54% density here. The densities for HVG matrices are estimates, not measured.
- **Full matrix before HVG selection** (Allen 10Xv3, 2.3M × 32K genes):
  - Assuming about 3,000 detected genes per cell (UNVERIFIED), that is about 7e9 nonzeros. That exceeds the int32 limit, so indices go to int64 at 12 B per nonzero, about **84 GB**.
  - It fits in RAM but not on the GPU. Select HVGs in backed or chunked mode, or with rapids-singlecell using Dask or out-of-core ([arXiv 2026](https://arxiv.org/abs/2603.02402)).
- **Derived data:**
  - 50 PCs: 0.4 GB.
  - kNN graph with k=15: 30M edges, well under 1 GB.
  - k-means distance buffer at K=338 subclasses: 2.7 GB.
  - All fine on GPU.
- **Hard limit:** Any method that needs an n × n matrix is impossible at 2M cells: pairwise distances, pairwise MI, SC3-style consensus, dense spectral or agglomerative clustering. That matrix is 4e12 × 4 B = **16 TB**. Even 100K × 100K is 40 GB.
  - This is why MICA switches to graph mode above 5,000 cells.
  - Computing MI from each cell to m landmark cells costs n × m. At m = 2,000 that is 16 GB, which is feasible and could be a scalable-scMINER angle to explore.

## 6. Biology or non-biology?

**Dataset 1:** scRNA-seq is not "ordinary multivariate" in spirit: roughly 20K sparse count features, usually more genes than cells. Cytometry (CyTOF) data meets the wording and stays single-cell immune, so a non-biology dataset is unnecessary. Keep PBMC14K as the reproduction anchor either way.

**Dataset 2:** A 2M-cell atlas is a legitimate large-scale problem, and the MERFISH brain set (550 genes) is a cheaper way to satisfy a strict reading. Stay with biology unless the instructor objects.

**Fallbacks:**

| Regime | Dataset | Size | Labels | License / source |
|---|---|---|---|---|
| D1 | Covertype | 581,012 × 54 | 7 classes | CC BY 4.0 ([UCI](https://archive.ics.uci.edu/dataset/31/covertype)) |
| D1 | Optical digits | 5,620 × 64 | 10 classes | CC BY 4.0 ([UCI](https://archive.ics.uci.edu/dataset/80/optical+recognition+of+handwritten+digits)) |
| D2 | MNIST8M | 8.1M × 784 | 10 classes | [LIBSVM](https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/multiclass.html); license UNVERIFIED |
| D2 | KDD Cup 1999 (full) | about 4.9M × 41 (UCI lists about 4M) | 23 classes | CC BY 4.0 ([UCI](https://archive.ics.uci.edu/dataset/130/kdd+cup+1999+data)); many duplicate rows (UNVERIFIED share) |
| D2 | HIGGS | 11M × 28 | Binary only, weak for clustering | CC BY 4.0, 2.6 GB ([UCI](https://archive.ics.uci.edu/dataset/280/higgs)) |
| D3 | Fashion-MNIST | 70,000 × 784 | 10 classes | MIT ([GitHub](https://github.com/zalandoresearch/fashion-mnist)) |
| D3 | COIL-20 | 1,440 images × 16,384 pixels (128×128) | 20 objects | [Columbia](https://www.cs.columbia.edu/CAVE/software/softlib/coil-20.php); counts UNVERIFIED this session |

Smartphone activity recognition (10,299 × 561, 6 classes, CC BY 4.0, [UCI](https://archive.ics.uci.edu/dataset/240/human+activity+recognition+using+smartphones)) meets Dataset 1 but has fewer than 700 variables for Dataset 3.

## 7. Recommendations

**Dataset 1**
- **Primary: Levine_32dim.** Manually gated labels, truly ordinary multivariate, published baselines (FlowSOM, PhenoGraph), a small file, and immune cells like the rest of the project.
- **Backup: PBMC14K.** scMINER's flagship dataset with FACS labels, if the instructor accepts scRNA-seq here. Samusik_01 is a second protein-level option.

**Dataset 2**
- **Primary: Allen whole mouse brain 10Xv3 (2.3M cells).**
  - Best-curated hierarchical labels at scale, so you can score at the class, subclass and supertype levels.
  - One chemistry.
  - Public S3 download with no login.
  - Non-commercial license is fine for coursework and academic papers.
- **Backup: Cao 2020 fetal atlas.** One 20 GB h5ad, CC BY 4.0, 4.06M cells.
- **If label trust matters most** (for example, for a paper): take a 2M-cell slice of Tahoe-100M (CC0, genotype labels).
- **If the PBMC theme matters most:** take a 2M-cell subsample of Parse 10M, but it requires the form and has circular labels.
- **Avoid:** MOCA (under 2M labeled cells) and HLCA full (mostly transferred labels).

**Dataset 3**
- **Primary: Zeisel 2015.** It is in scMINER's benchmark, so you can compare directly with the reported ARI. It also has two label levels (7 to 9 and 47) for granularity analysis.
- **Backup: Kolodziejczyk 2015 or Chu 2016.** Both have design-based labels. Also run Usoskin and Klein as a cheap panel, since scMINER reported results on them.

**Not verified this session**
- Licenses for the Levine/Samusik cytometry sets and for MNIST8M.
- Exact gene counts for Kolodziejczyk, Usoskin and Klein.
- The Hao 2021 cell count of 161,764.
- Whether scMINER runs at 2M cells.
- HVG density values, which are estimates.