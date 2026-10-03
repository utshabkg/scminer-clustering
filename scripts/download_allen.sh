#!/usr/bin/env bash
# Download the Allen Brain Cell Atlas whole mouse brain 10x v3 data (Yao et al. 2023).
# Raw counts only (13 region files, about 95 GB). The log2 copies are skipped to
# save disk, since we normalize ourselves. Public S3 bucket, no login. Resumable:
# rerun it and curl continues partial files. Files download 6 at a time, because a
# single connection was limited to about 2 MB/s. Usage: scripts/download_allen.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/data/raw/allen"; mkdir -p "$OUT"
B="https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com"
META=metadata/WMB-10X/20241115
for f in cell_metadata.csv gene.csv views/cell_metadata_with_cluster_annotation.csv; do
  curl -sS -C - --retry 5 -o "$OUT/$(basename $f)" "$B/$META/$f"
done
printf '%s\n' CB CTXsp HPF HY Isocortex-1 Isocortex-2 MB MY OLF P PAL STR TH \
  | xargs -P 6 -I{} sh -c 'curl -sS -C - --retry 10 -o "'"$OUT"'/WMB-10Xv3-{}-raw.h5ad" "'"$B"'/expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-{}-raw.h5ad" && echo "$(date +%H:%M:%S) done {}"'
du -sh "$OUT"; echo "all done"
