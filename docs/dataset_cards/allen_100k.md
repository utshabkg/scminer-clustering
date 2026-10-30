# Allen Brain Cell Atlas — WMB-10Xv3 Dataset Card

## Basics

- **Full name:** Allen Brain Cell Atlas – Whole Mouse Brain (ABC-WMB), specifically the WMB-10Xv3 dataset.
- **Short name used in our code:** `allen_wmb` (original dataset); `allen_100k` (processed experimental subset).
- **Paper it comes from:** Yao, Z., van Velthoven, C. T. J., Kunst, M., et al. (2023). *A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain.* *Nature*, 624, 317–332.
- **Link:** [https://doi.org/10.1038/s41586-023-06812-z](https://doi.org/10.1038/s41586-023-06812-z)
- **What was measured:** Single-cell RNA expression across the whole adult mouse brain, including neuronal and non-neuronal populations. Brain tissue was dissected, dissociated into individual cells, and analyzed using 10x Genomics v3 single-cell RNA sequencing (scRNA-seq). The broader atlas also incorporated MERFISH spatial transcriptomics for anatomical characterization.
- **Species:** Mouse (*Mus musculus*)
- **Tissue:** Whole adult mouse brain
- **Collection:** Brain tissue dissociation followed by scRNA-seq using 10x Genomics v3 chemistry.
- **Which of our three size regimes it covers:** Large scale. The original WMB-10Xv3 dataset contains 2,341,350 cells, while our project uses a randomly sampled subset of 100,000 cells for experimental evaluation.

## Size

- **Cells before our filtering:** 100,000 cells randomly sampled from the original 2,341,350-cell WMB-10Xv3 dataset using random seed 0.
- **Cells after our filtering:** 100,000 cells, as confirmed by `docs/dataset_facts.csv`. The processed experimental subset retains all 100,000 sampled cells.
- **Genes before and after:** The original WMB-10Xv3 dataset is described as containing 32,285 genes. The exact gene count immediately before our filtering is not recorded in `docs/dataset_facts.csv`. The processed 100,000-cell experimental subset contains 28,536 genes.
- **File size on disk, and where it lives:** The original WMB-10Xv3 expression-matrix archive is listed as 176.31 GB by the Allen Institute. Our project's 13 regional raw files total approximately 97 GB, as documented in `docs/allen_pilot.md`. The raw files are stored on the lab workstation at `/media/20TB/shared/sides_data/43.scMINER-Clustering/raw/allen_wmb/`. The exact disk size of the processed 100,000-cell subset is not stated.

## Labels

- **How many labels, and at how many levels:** The original atlas provides four hierarchical annotation levels: 34 classes, 338 subclasses, 1,201 supertypes, and 5,322 clusters (cell types). Our current `allen_100k` experiments evaluate clustering against the 34 broad cell classes. The finer levels are not the reported reference-label level in the current 100,000-cell pilot.
- **Where the labels came from:** Labels were produced through computational clustering of single-cell gene-expression profiles, followed by biological annotation using molecular markers, transcriptional similarities, and spatial information. They were not independently assigned through laboratory cell sorting.
- **Why that matters:** Because the reference labels originated from computational clustering, methods using similar gene-expression patterns and graph-based clustering strategies may obtain favorable agreement. Our project report acknowledges this potential bias and treats clustering comparisons at atlas scale as directional rather than absolute.
- **Smallest label group, with its cell count:** At the class level used in our experiments, the smallest population is "15 HY Gnrh1 Glut", containing 5 cells (0.005% of the 100,000-cell subset). The largest population is "31 OPC-Oligo", containing 20,283 cells (20.3%). These values are confirmed by `docs/dataset_facts.csv`.

## Source and terms

- **Download link we used:** Official WMB-10Xv3 dataset page: [https://alleninstitute.github.io/abc_atlas_access/descriptions/WMB-10Xv3.html](https://alleninstitute.github.io/abc_atlas_access/descriptions/WMB-10Xv3.html)
- **Original link, if different, and whether it still works:** Allen Brain Cell Atlas – Data Access: [https://alleninstitute.github.io/abc_atlas_access/intro.html](https://alleninstitute.github.io/abc_atlas_access/intro.html). The exact current availability of the link has not been independently rechecked.
- **Whole adult mouse brain taxonomy of cell types:** Not stated (separate link not supplied).
- **License and any restriction:** Creative Commons Attribution–NonCommercial 4.0 International (CC BY-NC 4.0). The dataset permits non-commercial research and redistribution subject to attribution and license conditions. Commercial use is restricted and may require separate permission from the rights holder. This is particularly relevant if the processed dataset or derived resources are redistributed.

## Preprocessing we applied

- **Quality control:** Our project applies a quality-control rule based on 4 median absolute deviations (MADs) of log-transformed total counts and detected genes per cell. No mitochondrial-content filter is applied. The final experimental dataset contains 100,000 cells.
- **Normalization:** Counts per million (CPM), followed by `log2(CPM + 1)`, as confirmed by `docs/dataset_facts.csv`.
- **Gene filtering:** The processed experimental subset contains 28,536 genes after preprocessing. The exact standalone gene-filtering threshold is not stated in `docs/dataset_facts.csv`.
- **Anything unusual about this dataset:** The original Allen WMB-10Xv3 dataset contains approximately 2.3 million cells, but our current experiments use a 100,000-cell subset. The processed dataset contains 28,536 genes and 34 broad cell classes, including rare classes with very few cells. The large gene-expression matrix creates substantial computational challenges for mutual-information-based clustering, particularly memory consumption and nearest-neighbor graph construction. Our Allen pilot documented these scalability limitations, motivating experiments with approximate nearest-neighbor methods and graph-based approaches.

## Notes

The Allen Brain Cell Atlas WMB-10Xv3 dataset is used to investigate the scalability and clustering performance of scMINER/MICA on large-scale single-cell RNA-seq data.

The original dataset contains 2,341,350 cells, while our current experimental subset contains 100,000 cells and 28,536 genes. Clustering performance is evaluated using the 34 broad cell-class annotations.

The reference labels originate from computational clustering and biological annotation rather than independent laboratory cell sorting. This introduces potential benchmarking bias when evaluating clustering algorithms against these labels.

Our completed Allen pilot compared MICA graph mode with Leiden-based clustering. The pilot revealed substantial computational costs associated with mutual-information nearest-neighbor graph construction, motivating subsequent scalability experiments.

The dataset is distributed under CC BY-NC 4.0, requiring attribution and restricting commercial use.
