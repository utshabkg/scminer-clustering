# Protocol v1, pre-tag checks

Outputs of the checks run on 2026-11-12, before the tag `matched-protocol-v1`, with the same code. They are kept apart from `results/v1/`, which holds only the confirmatory runs.

- Yan seed 0 and Kolodziejczyk seed 1 at B=100. Their original-protocol output (the first 25 subsamples, MICA and Leiden) equals the saved Assignment 3 results on every column but seconds: `regression.csv`, made by `scripts/check_v1_regression.py`.
- chooseR's own R code on the Yan seed 0 partitions: `chooser_r/`, `picks_chooser/`. R's tallied co-clustering frequency equals ours on 93 of 93 settings (`chooser_r/crosscheck_paper_yan_seed0.csv`).

These two seeds are part of the exploratory set (seeds 0 to 2), and their results were seen before the tag. See `docs/protocol_matched_v1.md`.
