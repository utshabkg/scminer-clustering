# MICA against Leiden, by gene set

## ARI

| dataset       |   published MICA |   MICA, authors' input |   authors' gene count | MICA graph mode, ours / authors' input   |   MICA, all |   MICA, 2k (seurat) |   MICA, 2k (cell_ranger) |   MICA, 2k (seurat_v3) |   Leiden, all |   Leiden, 2k (seurat) |   Leiden, 2k (cell_ranger) |   Leiden, 2k (seurat_v3) |
|:--------------|-----------------:|-----------------------:|----------------------:|:-----------------------------------------|------------:|--------------------:|-------------------------:|-----------------------:|--------------:|----------------------:|---------------------------:|-------------------------:|
| zeisel        |            0.744 |                  0.695 |                 14890 | 0.766 / 0.781                            |       0.704 |               0.728 |                    0.806 |                  0.844 |         0.845 |                 0.615 |                      0.692 |                    0.795 |
| pbmc14k       |            0.853 |                  0.853 |                  8846 |                                          |       0.854 |               0.653 |                    0.531 |                  0.737 |         0.876 |                 0.656 |                      0.588 |                    0.755 |
| usoskin       |            0.847 |                  0.803 |                 16383 |                                          |       0.885 |               0.84  |                    0.859 |                  0.877 |         0.588 |                 0.885 |                      0.919 |                    0.922 |
| kolodziejczyk |            1     |                  0.727 |                 13698 |                                          |       0.721 |               0.988 |                    1     |                  0.992 |         0.562 |                 0.996 |                      1     |                    1     |
| klein         |            0.98  |                  0.98  |                 24175 |                                          |       0.98  |               0.793 |                    0.728 |                  0.839 |         0.975 |                 0.95  |                      0.766 |                    0.976 |

All cells are the mean ARI over seeds under the paper's selection rule (true number of clusters). A star means at least one seed could not produce the true number of clusters and was scored at the closest count. Published values are ARI.

## AMI

| dataset       |   published MICA |   MICA, authors' input |   authors' gene count | MICA graph mode, ours / authors' input   |   MICA, all |   MICA, 2k (seurat) |   MICA, 2k (cell_ranger) |   MICA, 2k (seurat_v3) |   Leiden, all |   Leiden, 2k (seurat) |   Leiden, 2k (cell_ranger) |   Leiden, 2k (seurat_v3) |
|:--------------|-----------------:|-----------------------:|----------------------:|:-----------------------------------------|------------:|--------------------:|-------------------------:|-----------------------:|--------------:|----------------------:|---------------------------:|-------------------------:|
| zeisel        |            0.744 |                  0.739 |                 14890 | 0.742 / 0.761                            |       0.754 |               0.701 |                    0.776 |                  0.797 |         0.798 |                 0.632 |                      0.683 |                    0.745 |
| pbmc14k       |            0.853 |                  0.871 |                  8846 |                                          |       0.871 |               0.737 |                    0.741 |                  0.798 |         0.9   |                 0.758 |                      0.747 |                    0.814 |
| usoskin       |            0.847 |                  0.772 |                 16383 |                                          |       0.837 |               0.811 |                    0.837 |                  0.853 |         0.665 |                 0.842 |                      0.882 |                    0.897 |
| kolodziejczyk |            1     |                  0.784 |                 13698 |                                          |       0.781 |               0.979 |                    1     |                  0.985 |         0.69  |                 0.992 |                      1     |                    1     |
| klein         |            0.98  |                  0.963 |                 24175 |                                          |       0.963 |               0.777 |                    0.744 |                  0.845 |         0.958 |                 0.926 |                      0.817 |                    0.96  |

All cells are the mean AMI over seeds under the paper's selection rule (true number of clusters). A star means at least one seed could not produce the true number of clusters and was scored at the closest count. Published values are ARI.
