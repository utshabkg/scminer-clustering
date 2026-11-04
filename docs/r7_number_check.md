# R-7 number check, Assignment 3 (first pass)

Rhea, 2026-11-02. Flags only, nothing in `latex/` was changed.

**What was checked.** The `.tex` in the repo at commit 5888389 (last LaTeX change Oct 31, 13:50),
every section and every A3 table. Each aggregate number was traced to `results/a3_report_numbers.json`,
`results/a3_matched_numbers.json`, `results/a3_seed_variation.csv`,
`results/a3_matched_seed_variation.csv` or `results/a3_matched_mechanism.csv`. Numbers that are in
none of these were recomputed from the per-run `metrics.csv` files. Tables I and II were compared
cell by cell with the two seed-variation CSVs and match in every cell (four cells differ only
because the CSV stores 4 decimals; the unrounded values give the printed figure). Every A3
table (I to X) and Figure 1 is cited at least once in the text.

Line numbers refer to `latex/sections/*.tex`.

## Flags

- **The PDF of the paper is the old draft.** The PDF compiled from Overleaf on Nov 2 at 20:32
  (`Single_Cell_Genomic_Clustering_..._Interpretation.pdf`) is the version from before the matched
  runs. It has no matched table, summary table or transfer table, still has the "Plan through
  November 4" section, numbers the tables differently (IV is all four rules, V the common range,
  VI mode against size), and quotes the post hoc figures in the abstract ("0.12 against 0.07 per
  dataset"). The repo has the final text. Either Overleaf has not pulled the Oct 31 commits, or
  someone compiled an old version. This needs checking before anything is submitted on Nov 4,
  because if Overleaf is behind, every fix below lands in the wrong copy.

- **The subsample-draw paragraph is true only on the original grids** (`06_results.tex`
  l.198 to 203). "0.09 on average for both methods and by as much as 0.46" matches
  `a3_seed_variation.csv` (MICA 0.085, baseline 0.093, Buettner MICA plain maximum 0.46).
  On matched candidates, `a3_matched_seed_variation.csv` gives 0.12 overall, 0.089 for MICA
  and **0.157 for the baseline**, with a largest range of **0.50 for the baseline on Zeisel**
  (finest stable, 0.361 to 0.863). So on the runs the paper says it relies on, "for both
  methods" does not hold. The baseline moves almost twice as much as MICA, and the largest move
  is the baseline's. The same paragraph's "35 of 54 runs" is also the original-grid figure. On
  matched candidates it is 38 of 54. The subsection carries no grid label and sits after the
  matched subsections, so a reader will take it as matched.

- **One ARI range in the PBMC14K paragraph mixes the two grids** (`06_results.tex` l.140 to
  141, "On either grid, relaxing the threshold by 0.02 moves the baseline to 6 or 7 clusters at
  0.785 to 0.875"). On matched candidates, `a3_matched_numbers.json`
  (`threshold_sensitivity.pbmc14k`) gives 0.786 to 0.875. For the original grids,
  `a3_report_numbers.json` stores only the counts (6, 7, 7), not the ARIs. Recomputing them
  from the stability curves gives about 0.785 at the low end, and the earlier PDF draft said
  "0.785 to 0.873". So the low end comes from one grid and the high end from the other. The
  conclusion ("ahead of MICA") holds on both. Suggest giving each grid's range, and adding the
  original-grid ARIs to the JSON so the figure has a source.

- **Which grid "The baseline's best setting scores 0.876 to 0.878" means is unclear**
  (`06_results.tex` l.138). The sentence follows the matched description and ends "above any
  pick here", but 0.876 to 0.878 is the original grid (curve maxima 0.8777, 0.8763, 0.8762).
  On matched candidates the baseline's best is 0.871 to 0.875. Neither range is in either JSON.
  The sentence after it does label both grids for MICA's peak, so this one stands out.

- **Two different times are given for one MICA run on PBMC14K** (`06_results.tex` l.209 to
  211 against Table IV). The cost subsection says one run of MICA as shipped takes "about 11
  minutes". Table IV gives 13.2 min for the same pipeline (7.2 min build plus 6.0 min of
  node2vec and Louvain). The 11 minutes traces only to the Oct 31 entry in `PROGRESS.md` (HNSW
  429 s + node2vec 250 s + one Louvain 2.5 s), not to a results file. Both may be right, if
  Table IV's downstream time covers more than one Louvain resolution, but a reader will see two
  numbers for the same thing. Either reconcile them or say why they differ.

- **`a3_report_numbers.json` has wrong selection times** (not the paper). The JSON gives
  `stability_*.seconds_median` as 31.6 s for MICA and 25.1 s for the baseline. The per-run
  `metrics.csv` files give medians of 21.8 s and 18.7 s over the 27 small-dataset runs. Those
  are what Table III prints (22 and 19) and what its caption describes. The paper is right and
  the JSON is stale, but the JSON is meant to be the source, so it should be regenerated. In the
  same file, `default.seconds_median` is `NaN` for both methods.

## Minor

- **Rounding convention.** The matched MICA selection median is exactly 24.5 s and is printed
  as 24 (text l.218 and Table III). Elsewhere halves round up: 224.5 to 225 (l.214), 3.25 to
  3.3 and 1.25 to 1.3 (l.109, `07_discussion.tex` l.32). Pick one convention.
- **The Zeisel cost sentence has no grid label** (l.208 to 209). "103 s for MICA and 107 s for
  the baseline" are original-grid medians. The matched runs give 110 s and 113 s. The paragraph
  labels its matched figures later, so this is probably clear enough.
- **The JSON stores two figures for one mean.** `a3_report_numbers.json` has MICA's mean tied
  span on its own grid as 7.3 (`mechanism`) and as 6.89 (`common_range.means.mica_own`), and the
  baseline's as 2.7 and 2.33. The text uses 7.3 and 2.7, so nothing in the paper is wrong, but
  the second pair should be labelled or removed so nobody quotes it.

## Not in the five source files, but verified from per-run files

These are correct. They are listed so the next pass knows where they come from.

| Text | Value | Source |
|---|---|---|
| Intro l.5 | Zeisel baseline 3 or 4 clusters at 0.05 | `rule_stability_paper_zeisel*` curves |
| Results l.24 | 0.639 against 0.671 at the third seed | per-run metrics, original grids (no per-seed means in the JSON) |
| Table III | original-grid seed ranges [0.64, 0.70] etc. | per-run metrics, as above |
| Design l.82 to 84 | AMI above ARI by 0.17 / 0.12, too few in 35 of 54 | per-run metrics |
| Results l.142 | MICA peak 0.851 to 0.859, original grid | PBMC14K curves |
| Results l.138 to 140 | baseline default 0.853 to 0.867, MICA default 0.846 | `labelfree_paper_pbmc14k`, `labelfree_mica_paper_pbmc14k` |
| Results l.191, App. l.66 | baseline silhouette picks 2 on 7 of 9, 0.593 / 0.478 | `rule_matched_silhouette_*` |
| App. l.63, l.153 | space effect 0.13 / 0.16, Goolam 0.38, Chung 0.52 | `rule_silhouette_*` |
| Discussion l.39 | 0.267 on Buettner at k = 7 | Buettner curves |
| Results l.152 to 155 | floors 4, 3, 2; matched starts 4, 4, 3 | `grid_coverage_mica_graph/summary.csv`, `rule_matched_mode_size_*` |
| App. l.33 to 35 | Yan 19 to 7 at 0.925, Pollen 11 at 0.969, Klein +0.36 | `a3_common_range.csv` |
| Results l.248 to 257 | 2.4 min, 87 s, ratios 2.1 / 1.4 / 1.7 / 1.3 | `graph_cost_*`, `a2_report_numbers.json` |
| Results l.209 to 221 | 103 / 8 min, 93 / 7.5 min, 17 to 109 s, 14 to 225 s, 494, 468 to 1,040 | per-run metrics |
| Conclusion l.43 to 46 | 7.1 h, 140 GB, 3.2 TB, 1.8 TB | `docs/allen_pilot.md`, scaled by 23 |

## Could not verify from the repo

`logs/` is gitignored and is not on this machine, so these were not checked:

- MICA's distance matrix: "4 minutes" on Zeisel and "4 to 6 minutes" above 2,500 cells
  (`logs/rules_small.log`, per `PROGRESS.md`).
- The matched full-data sweep: "10 to 36 s" per small run and "about 16 minutes" per PBMC14K
  seed (`logs/matched_*.log`), and the "at most 36 s" in the Table III caption.
