## Basics

- **Full name and short name used in our code:**
    Full name: Zeisel 2015 mouse cortex and hippocampus
    Short name used in our code: zeisel
- **Paper it comes from:** Zeisel A. et al. (2015). Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. Science, 347(6226), 1138–1142. 
Link : https://doi.org/10.1126/science.aaa1934
- **What was measured:** single cells from the somatosensory cortex and hippocampus (CA1) of mice 20 to 31 days old, sequenced with STRT-seq, which tags each RNA molecule with a barcode (a UMI) so molecules can be counted
•	Species: Mus musculus
•	Tissue: Somatosensory cortex and hippocampal CA1
- **Which of our three size regimes it covers:** High dimensional. The downloaded matrix contains 20,006 genes and 3,005 cells, so the number of features is much larger than the number of observations.

## Size

- **Cells before our filtering:** 3,005 cells
- **Cells after our filtering:** 2,987 cells 
- **Genes before and after:**
    Before: 20,006 genes (19,972 regular genes plus 34 mitochondrial genes)
    After: 18,907 genes
- **File size on disk, and where it lives** On the lab workstation's shared disk, in /media/20TB/shared/sides_data/43.scMINER-Clustering/raw/zeisel/ (24 MB). It is not on GitHub

## Labels

- **How many labels, and at how many levels** 7 broad classes and 47 subclasses. The paper describes 9 major classes, but the released labels group astrocytes with ependymal cells and endothelial with mural cells, which gives 7. 
  184 cells have no subclass and are marked (none), 189 before filtering
  Regions: 1,675 cells from the somatosensory cortex and 1,312 from the hippocampus (CA1)
- **Where the labels came from:** The Zeisel et al. authors generated the cell-type groups by clustering the same gene-expression data using their BackSPIN bi-clustering method and then assigned biological cell-type names using marker genes and other biological evidence.
- **Why that matters:** These are clustering-derived reference labels rather than purely independent laboratory ground-truth labels. Consequently, an evaluated clustering method that produces partitions similar to BackSPIN may appear particularly favorable when compared against these annotations.
- **Smallest label group, with its cell count** microglia, 94 cells. Smallest subclass: ClauPyr, 5 cells

## Source and terms

- **Download link we used:** Bioconductor's scRNAseq collection, a public library of biology datasets. Dataset zeisel-brain-2015, version 2023-12-14: https://bioconductor.org/packages/scRNAseq/
- **Original link, if different, and whether it still works** http://linnarssonlab.org/cortex/ The page still loads, but its data files now give an access error (403). The raw counts are also on GEO as GSE60361: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60361
- **License and any restriction** The original Zeisel paper/Linnarsson data source does not state a specific dataset license. The Bioconductor scRNAseq package through which we obtained the processed dataset is distributed under CC0.

## Preprocessing we applied

- **Quality control:** We calculated each cell's total UMI count and number of detected genes on a log scale. Cells identified as >4 median absolute deviations (MAD) from the median outliers according to these measures were removed. This removed 18 cells, leaving 2,987 cells.
- **Normalization:** For each cell, raw counts were converted to counts per million (CPM), followed by a log2 transformation.
- **Gene filtering:** Genes were retained if they had nonzero expression in at least 3 retained cells. This reduced the feature count from 20,006 to 18,907.
- **Anything unusual about this dataset** Mitochondrial fraction was not used as an additional cell-level filtering criterion in our preprocessing.

## Notes

The original study analyzed 3,005 cells and described 9 major classes and 47 molecular subclasses. The broad labels distributed with the dataset combine some of the original classes, resulting in 7 broad labels in our downloaded metadata. Some cells lack a subclass annotation and are recorded as "(none)" (189 before our filtering and 184 afterward). The original paper's downstream analysis retained 2,903 cells, whereas our QC pipeline retains 2,987, because we apply our own common preprocessing procedure. The original Linnarsson laboratory page remains relevant as the historical source, while we use the versioned Bioconductor copy for reproducibility.
