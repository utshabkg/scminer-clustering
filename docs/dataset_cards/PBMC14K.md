# Peripheral Blood Mononuclear Cells 14K (PBMC14K) Dataset Card

## Basics

- **Full name:** Peripheral Blood Mononuclear Cells 14K.
- **Short name used in our code:** `paper_pbmc14k` (authors' processed benchmark); `pbmc14k` (separately processed version). PBMC14K is the dataset's common name.
- **Paper it comes from:** Zheng, G. X. Y., Terry, J. M., Belgrader, P., et al. (2017). *Massively parallel digital transcriptional profiling of single cells.* *Nature Communications*, 8, 14049.
- **Link:** [https://doi.org/10.1038/ncomms14049](https://doi.org/10.1038/ncomms14049)
- **What was measured:** Single-cell RNA expression (3′-end transcript counts/UMIs) in human peripheral blood mononuclear cells (PBMCs). Specific immune-cell populations were antibody-purified or sorted before droplet-based single-cell RNA sequencing (10x Genomics). These are blood immune cells, not solid tissue.
- **Species:** Human (*Homo sapiens*).
- **Tissue:** Peripheral blood.
- **Collection:** Antibody-based purification/sorting of immune-cell populations followed by droplet-based 3′ single-cell RNA sequencing.
- **Which of our three size regimes it covers:** Ordinary. The authors' released benchmark (`paper_pbmc14k`) contains 13,605 cells and 8,846 genes.

## Size

- **Cells before our filtering:** Our team applied no additional filtering to the authors' released `paper_pbmc14k` benchmark. Historically, the scMINER PBMC14K subset began with 14,000 cells (seven sorted populations of 2,000 each), selected from an earlier 20,000-cell reference matrix containing ten populations.
- **Cells after our filtering:** 13,605 cells in the authors' released processed benchmark (`paper_pbmc14k`), confirmed by `docs/dataset_facts.csv`. The authors' preprocessing removed 395 of the initial 14,000 cells (approximately 2.82%); our team did not further filter this benchmark.
- **Genes before and after:** The original scMINER preprocessing description reports 17,986 genes before filtering and 8,846 retained afterward. The earlier 20,000-cell source matrix reportedly contained 21,592 genes. Our team's authoritative `paper_pbmc14k` processed input has 8,846 genes; the original pre-filtering gene count is historical source information, not a count independently reproduced by our team's benchmark pipeline.
- **File size on disk, and where it lives:** The raw PBMC14K data are stored on the lab workstation's shared disk at `/media/20TB/shared/sides_data/43.scMINER-Clustering/raw/PBMC14K/` and are not stored in the GitHub repository. The exact on-disk size of our processed benchmark file is not stated. The public archive lists `10_PBMC14K.RawData.zip` (27.7 MB) and `10_PBMC14K.Processed_eset.rds` (22.2 MB); these download sizes should not be treated as the local processed file size.

## Labels

- **How many labels, and at how many levels:** Seven antibody-sorted immune-cell populations at one supplied reference-label level: CD14+ monocytes (Monocyte), CD19+ B cells (B), CD4+/CD25+ regulatory T cells (CD4Treg), CD4+/CD45RA+/CD25− naïve T cells (CD4TN), CD4+/CD45RO+ memory T cells (CD4TCM), CD56+ natural killer cells (NK), and CD8+/CD45RA+ naïve cytotoxic T cells (CD8TN). Full and abbreviated names do not constitute separate biological hierarchy levels. [scMINER documentation](https://jyyulab.github.io/scMINER/bookdown/get-started.html).
- **Where the labels came from:** Antibody-based laboratory sorting before sequencing, inherited from the Zheng et al. purified PBMC reference populations. The labels were not produced by clustering this PBMC14K expression matrix. scMINER retained seven of the original ten populations, omitting CD34+, broad CD4+ helper T, and broad CD8+ cytotoxic T populations because of sorting purity or phenotype overlap, as described in the source documentation.
- **Why that matters:** Clustering-derived reference labels can favor methods resembling the algorithms used to create them. PBMC14K's independently sorted populations provide a less circular benchmark for external clustering agreement (such as ARI and NMI). However, sorting errors and transcriptional similarities between T-cell subsets mean these labels are not infallible.
- **Smallest label group, with its cell count:** Monocyte (CD14+ monocytes), with 1,837 cells (13.5% of 13,605). The largest class is CD4TN, with 1,994 cells. These values correspond to the `paper_pbmc14k` row in `docs/dataset_facts.csv`. Before the original authors' quality control, each of the seven selected populations contained 2,000 cells.

## Source and terms

- **Download link we used:** [https://zenodo.org/records/15040179](https://zenodo.org/records/15040179)
- **Original link, if different, and whether it still works:** Original PBMC sequencing experiments: [SRA SRP073767](https://www.ncbi.nlm.nih.gov/sra/?term=SRP073767). Reformatted 20,000-cell sorted PBMC reference data: [Zenodo 3357167](https://zenodo.org/records/3357167). scMINER PBMC14K benchmarking archive: [Zenodo 15040179](https://zenodo.org/records/15040179). Current availability of each link has not been independently rechecked.
- **License and any restriction:** Not stated for the exact downloaded PBMC14K benchmark files; confirm the archive's licensing terms before redistribution.

## Preprocessing we applied

- **Quality control:** Our project uses the scMINER authors' released processed PBMC14K dataset (`paper_pbmc14k`) without additional QC filtering. The original scMINER procedure used `filterSparseEset(pbmc14k_raw.eset, filter_mode = "auto", filter_type = "both")` with reported thresholds of 458–3,694 UMIs, at least 221 detected genes, and mitochondrial UMI fraction at most 4.08%. The source description reports 395 removed cells and 13,605 retained. These are the *authors' historical QC settings*, not filters applied again by our team. The benchmark-processing stage is recorded as "as released by the scMINER authors" in `docs/dataset_facts.csv`.
- **Normalization:** Counts per million (CPM), followed by `log2(CPM + 1)`, as recorded in `docs/dataset_facts.csv`.
- **Gene filtering:** The original scMINER procedure reportedly retained genes detected in at least 70 cells. Our team applied no additional gene filtering to the authors' released `paper_pbmc14k` benchmark.
- **Anything unusual about this dataset:** The scMINER developers repaired gene symbols against GRCh37.82, selected seven populations from ten antibody-purified samples, and removed genes absent in the selected cells. "PBMC14K" refers to the original 14,000-cell pre-QC collection rather than the 13,605-cell released analysis matrix. Some selected T-cell groups are transcriptionally similar despite different sorting markers. The repository also contains a separately processed `pbmc14k` dataset (13,646 cells; 15,064 genes); those figures must not be mixed with the authors' processed `paper_pbmc14k` benchmark (13,605 cells; 8,846 genes).

## Notes

PBMC14K was used in scMINER/MICA clustering benchmarks and analyses of immune-cell regulatory-driver activity. The scMINER study highlights challenges in distinguishing related CD4 T-cell populations, including regulatory and central-memory T cells.

Unlike reference labels derived by clustering the same expression data, PBMC14K's labels originate from antibody-sorted sample identity. This reduces circularity in external clustering evaluation, while leaving limitations from sorting purity, transcriptional overlap, potential batch effects, and class imbalance after filtering.

For reporting experimental results, use the `paper_pbmc14k` counts and preprocessing details in `docs/dataset_facts.csv` and `docs/dataset_facts.md`. Do not substitute statistics from the separately processed `pbmc14k` version.
