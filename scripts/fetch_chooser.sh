#!/usr/bin/env bash
# Fetch chooseR's code at the pinned commit into the data folder (not the repo, since it is
# GPL-3.0) and check the files against the hashes recorded when protocol v1 was fixed.
set -eu
SHA=410bdbfaebcbc57df45176491545a1165e3c6fc5
cd "$(dirname "$0")/.."
D=data/external/chooseR_$SHA
mkdir -p "$D/R" "$D/examples"
for f in R/pipeline.R examples/1_seurat_pipeline.R LICENSE; do
  [ -f "$D/$f" ] || curl -fsSL "https://raw.githubusercontent.com/rbpatt2019/chooseR/$SHA/$f" -o "$D/$f"
done
sha256sum -c - <<SUMS
ca0016318ffa431395ff20669861970294fda2c329ef40ec4b6275f1949201c3  $D/R/pipeline.R
c56f9e43ebe9799ed4c7d25e47390c14a6e8be00ed2a2496b97fca2045579bbb  $D/examples/1_seurat_pipeline.R
3972dc9744f6499f0f9b2dbf76696f2ae7ad8af9b23dde66d6af86c9dfb36986  $D/LICENSE
SUMS
