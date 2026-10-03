# What MICA's similarity measures

## Co-detection

Share of the mutual information between two cells that comes only from which genes are detected in both (median over 20,000 random pairs), and how MICA does when it sees detection only.

| dataset       |   codetection_share_median |   share_same_type |   share_diff_type |   detected_values_in_bin0 |   ari_detection_only |   agreement_with_full_mica |   knn_overlap |   ari_full_mica |
|:--------------|---------------------------:|------------------:|------------------:|--------------------------:|---------------------:|---------------------------:|--------------:|----------------:|
| klein         |                      0.253 |             0.248 |             0.255 |                     0     |                0.582 |                      0.591 |         0.424 |           0.98  |
| kolodziejczyk |                      0.613 |             0.611 |             0.613 |                     0.096 |                0.711 |                      0.976 |         0.611 |           0.721 |
| usoskin       |                      0.574 |             0.58  |             0.573 |                     0.014 |                0.531 |                      0.6   |         0.612 |           0.885 |
| zeisel        |                      0.62  |             0.603 |             0.626 |                     0     |                0.698 |                      0.873 |         0.584 |           0.706 |

## Depth

ARI against the labels after thinning every cell's counts to a fraction of its reads.

|                                 |   100% |   75% |   50% |   25% |
|:--------------------------------|-------:|------:|------:|------:|
| ('klein', 'leiden_pca')         |  0.977 | 0.971 | 0.962 | 0.942 |
| ('klein', 'mica')               |  0.98  | 0.982 | 0.977 | 0.966 |
| ('kolodziejczyk', 'leiden_pca') |  0.562 | 0.562 | 0.562 | 0.562 |
| ('kolodziejczyk', 'mica')       |  0.721 | 0.721 | 0.721 | 0.721 |
| ('usoskin', 'leiden_pca')       |  0.588 | 0.597 | 0.623 | 0.407 |
| ('usoskin', 'mica')             |  0.885 | 0.866 | 0.572 | 0.268 |
| ('zeisel', 'leiden_pca')        |  0.846 | 0.849 | 0.862 | 0.881 |
| ('zeisel', 'mica')              |  0.706 | 0.704 | 0.704 | 0.694 |

Agreement (ARI) with the clusters found at full depth.

|                                 |   25% |   50% |   75% |
|:--------------------------------|------:|------:|------:|
| ('klein', 'leiden_pca')         | 0.956 | 0.973 | 0.983 |
| ('klein', 'mica')               | 0.982 | 0.993 | 0.994 |
| ('kolodziejczyk', 'leiden_pca') | 1     | 1     | 1     |
| ('kolodziejczyk', 'mica')       | 1     | 1     | 1     |
| ('usoskin', 'leiden_pca')       | 0.478 | 0.576 | 0.85  |
| ('usoskin', 'mica')             | 0.284 | 0.634 | 0.976 |
| ('zeisel', 'leiden_pca')        | 0.929 | 0.942 | 0.98  |
| ('zeisel', 'mica')              | 0.925 | 0.967 | 0.983 |
