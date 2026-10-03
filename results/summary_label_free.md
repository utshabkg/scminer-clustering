# Choosing the number of clusters without labels

Ten benchmark sets on the authors' own input. Told k is the paper's rule. Silhouette means each method tries k = 2 to 20 (or its resolution sweep) and keeps the best silhouette in its own embedding.

## Overall

| method       |   mean ARI told k |   mean ARI silhouette |   mean loss | cases chosen k equals true k   |
|:-------------|------------------:|----------------------:|------------:|:-------------------------------|
| Leiden       |             0.692 |                 0.556 |       0.136 | 0 of 40                        |
| MICA         |             0.819 |                 0.577 |       0.242 | 5 of 40                        |
| hierarchical |             0.618 |                 0.399 |       0.219 | 2 of 40                        |
| k-means      |             0.601 |                 0.42  |       0.182 | 2 of 40                        |

## Loss from not knowing k, by gene setting

|                                              |   mean |   median |   max |
|:---------------------------------------------|-------:|---------:|------:|
| ('Leiden', '2000 genes (cell_ranger)')       |  0.075 |    0.038 | 0.37  |
| ('Leiden', '2000 genes (seurat)')            |  0.179 |    0.194 | 0.532 |
| ('Leiden', '2000 genes (seurat_v3)')         |  0.129 |    0.031 | 0.55  |
| ('Leiden', 'all genes')                      |  0.16  |    0.132 | 0.481 |
| ('MICA', '2000 genes (cell_ranger)')         |  0.218 |    0.096 | 0.685 |
| ('MICA', '2000 genes (seurat)')              |  0.227 |    0.1   | 0.674 |
| ('MICA', '2000 genes (seurat_v3)')           |  0.215 |    0.115 | 0.492 |
| ('MICA', 'all genes')                        |  0.308 |    0.254 | 0.569 |
| ('hierarchical', '2000 genes (cell_ranger)') |  0.201 |    0.219 | 0.475 |
| ('hierarchical', '2000 genes (seurat)')      |  0.172 |    0.116 | 0.471 |
| ('hierarchical', '2000 genes (seurat_v3)')   |  0.235 |    0.147 | 0.638 |
| ('hierarchical', 'all genes')                |  0.269 |    0.198 | 0.846 |
| ('k-means', '2000 genes (cell_ranger)')      |  0.201 |    0.239 | 0.473 |
| ('k-means', '2000 genes (seurat)')           |  0.216 |    0.196 | 0.464 |
| ('k-means', '2000 genes (seurat_v3)')        |  0.196 |    0.128 | 0.561 |
| ('k-means', 'all genes')                     |  0.113 |    0.083 | 0.536 |

## ARI with silhouette-chosen k

|                                               |   Leiden |   MICA |   hierarchical |   k-means |
|:----------------------------------------------|---------:|-------:|---------------:|----------:|
| ('buettner', '2000 genes (cell_ranger)')      |    0.264 |  0.212 |          0.291 |     0.297 |
| ('buettner', '2000 genes (seurat)')           |    0.281 |  0.315 |          0.293 |     0.318 |
| ('buettner', '2000 genes (seurat_v3)')        |    0.284 |  0.308 |          0.259 |     0.307 |
| ('buettner', 'all genes')                     |    0.231 |  0.474 |          0.144 |     0.291 |
| ('chung', '2000 genes (cell_ranger)')         |    0.505 |  0.554 |          0.009 |     0.106 |
| ('chung', '2000 genes (seurat)')              |    0.44  |  0.555 |          0.001 |     0.037 |
| ('chung', '2000 genes (seurat_v3)')           |    0.504 |  0.553 |          0.009 |     0.019 |
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
| ('pollen', '2000 genes (seurat)')             |    0.838 |  0.88  |          0.924 |     0.914 |
| ('pollen', '2000 genes (seurat_v3)')          |    0.879 |  0.874 |          0.937 |     0.91  |
| ('pollen', 'all genes')                       |    0.87  |  0.89  |          0.102 |     0.902 |
| ('usoskin', '2000 genes (cell_ranger)')       |    0.722 |  0.798 |          0.314 |     0.316 |
| ('usoskin', '2000 genes (seurat)')            |    0.36  |  0.898 |          0.254 |     0.184 |
| ('usoskin', '2000 genes (seurat_v3)')         |    0.362 |  0.86  |          0.356 |     0.144 |
| ('usoskin', 'all genes')                      |    0.379 |  0.266 |          0.125 |     0.089 |
| ('yan', '2000 genes (cell_ranger)')           |    0.765 |  0.498 |          0.796 |     0.798 |
| ('yan', '2000 genes (seurat)')                |    0.804 |  0.432 |          0.712 |     0.738 |
| ('yan', '2000 genes (seurat_v3)')             |    0.804 |  0.421 |          0.794 |     0.794 |
| ('yan', 'all genes')                          |    0.765 |  0.496 |          0.74  |     0.692 |
| ('zeisel', '2000 genes (cell_ranger)')        |    0.757 |  0.816 |          0.444 |     0.536 |
| ('zeisel', '2000 genes (seurat)')             |    0.416 |  0.763 |          0.099 |     0.102 |
| ('zeisel', '2000 genes (seurat_v3)')          |    0.689 |  0.797 |          0.45  |     0.53  |
| ('zeisel', 'all genes')                       |    0.367 |  0.518 |          0.35  |     0.364 |
