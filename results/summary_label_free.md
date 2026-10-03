# Choosing the number of clusters without labels

Ten benchmark sets on the authors' own input. Told k is the paper's rule. Silhouette means each method tries k = 2 to 20 (or its resolution sweep) and keeps the best silhouette in its own embedding.

## Overall

| method       |   mean ARI told k |   mean ARI silhouette |   mean loss | cases chosen k equals true k   |
|:-------------|------------------:|----------------------:|------------:|:-------------------------------|
| Leiden       |             0.68  |                 0.542 |       0.138 | 0 of 40                        |
| MICA         |             0.802 |                 0.572 |       0.229 | 5 of 40                        |
| hierarchical |             0.605 |                 0.397 |       0.208 | 2 of 40                        |
| k-means      |             0.588 |                 0.417 |       0.172 | 2 of 40                        |

## Loss from not knowing k, by gene setting

|                                              |   mean |   median |   max |
|:---------------------------------------------|-------:|---------:|------:|
| ('Leiden', '2000 genes (cell_ranger)')       |  0.075 |    0.038 | 0.37  |
| ('Leiden', '2000 genes (seurat)')            |  0.169 |    0.194 | 0.509 |
| ('Leiden', '2000 genes (seurat_v3)')         |  0.148 |    0.046 | 0.542 |
| ('Leiden', 'all genes')                      |  0.16  |    0.132 | 0.481 |
| ('MICA', '2000 genes (cell_ranger)')         |  0.218 |    0.096 | 0.685 |
| ('MICA', '2000 genes (seurat)')              |  0.213 |    0.096 | 0.674 |
| ('MICA', '2000 genes (seurat_v3)')           |  0.182 |    0.095 | 0.479 |
| ('MICA', 'all genes')                        |  0.304 |    0.254 | 0.569 |
| ('hierarchical', '2000 genes (cell_ranger)') |  0.201 |    0.219 | 0.475 |
| ('hierarchical', '2000 genes (seurat)')      |  0.151 |    0.114 | 0.471 |
| ('hierarchical', '2000 genes (seurat_v3)')   |  0.21  |    0.188 | 0.638 |
| ('hierarchical', 'all genes')                |  0.269 |    0.198 | 0.846 |
| ('k-means', '2000 genes (cell_ranger)')      |  0.201 |    0.239 | 0.473 |
| ('k-means', '2000 genes (seurat)')           |  0.205 |    0.204 | 0.464 |
| ('k-means', '2000 genes (seurat_v3)')        |  0.168 |    0.057 | 0.561 |
| ('k-means', 'all genes')                     |  0.113 |    0.083 | 0.536 |

## ARI with silhouette-chosen k

|                                               |   Leiden |   MICA |   hierarchical |   k-means |
|:----------------------------------------------|---------:|-------:|---------------:|----------:|
| ('buettner', '2000 genes (cell_ranger)')      |    0.264 |  0.212 |          0.291 |     0.297 |
| ('buettner', '2000 genes (seurat)')           |    0.296 |  0.412 |          0.246 |     0.309 |
| ('buettner', '2000 genes (seurat_v3)')        |    0.331 |  0.435 |          0.27  |     0.333 |
| ('buettner', 'all genes')                     |    0.231 |  0.474 |          0.144 |     0.291 |
| ('chung', '2000 genes (cell_ranger)')         |    0.505 |  0.554 |          0.009 |     0.106 |
| ('chung', '2000 genes (seurat)')              |    0.241 |  0.059 |          0.001 |     0.001 |
| ('chung', '2000 genes (seurat_v3)')           |    0.165 |  0.545 |          0.009 |     0.022 |
| ('chung', 'all genes')                        |    0.046 |  0.557 |          0.109 |     0.158 |
| ('goolam', '2000 genes (cell_ranger)')        |    0.63  |  0.243 |          0.63  |     0.63  |
| ('goolam', '2000 genes (seurat)')             |    0.86  |  0.24  |          0.63  |     0.63  |
| ('goolam', '2000 genes (seurat_v3)')          |    0.86  |  0.208 |          0.272 |     0.272 |
| ('goolam', 'all genes')                       |    0.619 |  0.267 |          0.272 |     0.272 |
| ('klein', '2000 genes (cell_ranger)')         |    0.5   |  0.724 |          0.392 |     0.372 |
| ('klein', '2000 genes (seurat)')              |    0.5   |  0.784 |          0.428 |     0.323 |
| ('klein', '2000 genes (seurat_v3)')           |    0.5   |  0.786 |          0.384 |     0.336 |
| ('klein', 'all genes')                        |    0.501 |  0.501 |          0.415 |     0.413 |
| ('kolodziejczyk', '2000 genes (cell_ranger)') |    0.63  |  0.63  |          0.629 |     0.63  |
| ('kolodziejczyk', '2000 genes (seurat)')      |    0.629 |  0.63  |          0.629 |     0.63  |
| ('kolodziejczyk', '2000 genes (seurat_v3)')   |    0.63  |  0.63  |          0.629 |     0.63  |
| ('kolodziejczyk', 'all genes')                |    0.63  |  0.497 |          0.716 |     0.716 |
| ('pbmc14k', '2000 genes (cell_ranger)')       |    0.396 |  0.543 |          0.094 |     0.094 |
| ('pbmc14k', '2000 genes (seurat)')            |    0.396 |  0.59  |          0.094 |     0.094 |
| ('pbmc14k', '2000 genes (seurat_v3)')         |    0.397 |  0.574 |          0.094 |     0.094 |
| ('pbmc14k', 'all genes')                      |    0.397 |  0.623 |          0.225 |     0.214 |
| ('pollen', '2000 genes (cell_ranger)')        |    0.842 |  0.878 |          0.9   |     0.887 |
| ('pollen', '2000 genes (seurat)')             |    0.815 |  0.894 |          0.918 |     0.918 |
| ('pollen', '2000 genes (seurat_v3)')          |    0.815 |  0.908 |          0.912 |     0.916 |
| ('pollen', 'all genes')                       |    0.87  |  0.89  |          0.102 |     0.902 |
| ('usoskin', '2000 genes (cell_ranger)')       |    0.722 |  0.798 |          0.314 |     0.316 |
| ('usoskin', '2000 genes (seurat)')            |    0.354 |  0.903 |          0.249 |     0.126 |
| ('usoskin', '2000 genes (seurat_v3)')         |    0.36  |  0.898 |          0.36  |     0.1   |
| ('usoskin', 'all genes')                      |    0.379 |  0.266 |          0.125 |     0.089 |
| ('yan', '2000 genes (cell_ranger)')           |    0.765 |  0.498 |          0.796 |     0.798 |
| ('yan', '2000 genes (seurat)')                |    0.804 |  0.432 |          0.712 |     0.738 |
| ('yan', '2000 genes (seurat_v3)')             |    0.804 |  0.421 |          0.794 |     0.794 |
| ('yan', 'all genes')                          |    0.765 |  0.496 |          0.74  |     0.692 |
| ('zeisel', '2000 genes (cell_ranger)')        |    0.757 |  0.816 |          0.444 |     0.536 |
| ('zeisel', '2000 genes (seurat)')             |    0.416 |  0.763 |          0.099 |     0.102 |
| ('zeisel', '2000 genes (seurat_v3)')          |    0.689 |  0.797 |          0.45  |     0.53  |
| ('zeisel', 'all genes')                       |    0.367 |  0.518 |          0.35  |     0.364 |
