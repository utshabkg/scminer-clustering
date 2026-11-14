# R-7 number check, Assignment 4 draft

Rhea, 2026-11-13. Flags only. Nothing in `latex/` or Overleaf was changed. Utshab applies the fixes.

**What was checked.** The `.tex` on GitHub `main` at 19188e3 (title "When Clustering Methods Must
Choose for Themselves"), every section and every table it inputs. Line numbers refer to
`latex/sections/*.tex`.

1. **Source chain.** `scripts/a4_protocol_v1.py --tables 1` was rerun on a clean clone from
   `results/v1`. It reproduces `results/a4_protocol_v1_numbers.json`, `latex/numbers_a4.tex` and
   all five A4 tables (`a4_criteria`, `a4_confirm`, `a4_decomp`, `a4_designs`, `a4_pbmc`) byte for
   byte. All 141 uses of `\V{}` (116 keys) are defined, and each value equals `latex_values` in
   the JSON.
2. **Independent recomputation.** The headline means, per-seed leads at all ten seeds, the
   leave-one-out gaps, the Klein picks and every PBMC14K pick were recomputed from
   `results/v1/picks*`, `curves`, `chooser_r` and `provenance` with a separate script, not the
   generator. All agree.
3. **Claims about each number** (higher, level, holds, every seed, without any one dataset). Each
   was checked against the values. Results below.
4. **Tables and figures.** All 13 tables and both figures are cited in the text.
5. **Provenance.** All 190 provenance files record commit 620f940, `code_matches_tag: true`, clean
   tree. The annotated tag `matched-protocol-v1` is on GitHub and dated Nov 12, 11:03.

## Flags to fix tonight

- **PBMC14K cost, `06_results.tex` l.207.** "a median of four and a half to six minutes" per MICA
  subsample. The provenance files give 293 to 363 s per seed, which is 4.9 to 6.0 minutes.
  Suggest "about five to six minutes". The rest of the sentence holds (74 to 90 s for the
  baselines, 8.1 to 10.1 h and 2.1 to 2.5 h per run of 100).
- **Two sentences describe the same scorer as if different, l.171 to 174.** "With all 100
  subsamples Leiden takes that partition at only two seeds" and "With the co-clustering on 3,000
  cells (S1) the two then tie at three or more seeds" are both S1. The counts are right (two
  seeds, then 3 ties under the lowest setting, 4 under every setting, 5 under the highest).
  Suggest joining them: "With all 100 subsamples on the same 3,000 cells (S1), Leiden takes that
  partition at only two seeds, the medians are level again and the two tie at three or more
  seeds under every variant."
- **"Much of MICA's lead comes from Buettner", l.41 to 42.** Buettner alone adds 0.081 to the
  mean gap, which is more than the whole lead of 0.056, and without it Leiden leads (-0.028).
  Suggest "all of MICA's lead", which matches the abstract ("rests on one dataset").
- **"Level" has no stated cut-off on the small datasets.** Matched plain maximum, gap 0.002, is
  "level". Original-grid finest stable, 0.647 against 0.662 (gap 0.015, MICA higher at 5 of 7
  seeds), is "level" (l.34, Abstract, Discussion l.38). S1 original-grid finest stable, gap 0.022
  (MICA higher at 7 of 7), is "a MICA lead" (l.107 to 110). 0.015 and 0.022 sit close together.
  Suggest one sentence in Section V saying a gap under 0.02 is called level, the same cut C4
  uses on PBMC14K.
- **"Barely moves", l.61.** MICA's mean is 0.452 on both designs, identical before rounding.
  Suggest "does not move".
- **Denominators differ from the stated default, Section VI-E.** Section VI opens by saying a
  figure is over seeds 3 to 9 on the nine small datasets unless stated. "106 of 190", "36 of
  190", "83 of 190", "282 of 8,923" and "84 of 95" / "66 of 95" count all ten seeds plus PBMC14K
  (and both designs for the 190). Suggest "of 190 runs, every seed and all ten datasets, both
  designs" at first use.

## Correct, but typed by hand

The request says every number comes from `numbers_a4.tex`. These protocol v1 figures do not.
Each was checked against the per-run files and is right. For Nov 30 they should move into
`write_numbers()`. For tonight nothing needs to change.

| Where | Text | Source checked |
|---|---|---|
| Abstract l.15, Results l.190, Conclusion l.9 | "four of five seeds" | `\V{C4.ties}` exists, use it |
| Results l.31 | gap at seeds 0 to 2 was 0.021 | `.s02` macros, 0.473 - 0.452 |
| Results l.15 to 18 | 0.635 to 0.473, 0.662 to 0.543, 0.473 to 0.452 | `\V{...s02}` macros exist and match |
| Results l.64 | Klein 2 clusters at 0.501, 4 at 0.973 at five seeds, 5 at 0.885 or 0.888 at two | `picks` |
| Results l.99 | MICA failures 2 or 3 clusters scoring 1.0, Leiden two thirds scoring 0 (30 of 45), $k$-means 13 to 20 at 0 (140 of 184) and 29 at 2 to 4 scoring 1.0, Klein and Kolodziejczyk every run, Buettner none | `chooser_r`, `picks_chooser` |
| Results l.167 to 178 | every seed count in the PBMC14K paragraph | `picks`, all scorers and variants |
| Results l.184 to 186 | medians 0.993 to 0.999, chooseR bound 0.48 to 0.66, ours within 0.004 | `chooser_r`, `curves` (lowest chooseR bound on the original grid is 0.477, 0.471 if matched settings are included) |
| Results l.187 to 189 | 8 clusters at every seed, $k$-means 2 in every run of Table VII but one | `picks` |
| Results l.204 to 210 | small-data cost, "under a minute" for chooseR scoring (max 0.99 min) | `provenance` |
| Intro l.5 to 6 | Zeisel 3 or 4 at 0.05, 20 or 21 at 1.0 | `curves` (20 at seeds 3, 5, 6, 7, above the cap at the rest) |
| Results l.14 | label-fixed 0.809 and 0.614 | mean of `a2_report_numbers.json` reproduction medians (0.8087, 0.6136) |
| Criteria table, Assignment 3 column | 0.197 vs 0.118, 0.534 vs 0.231, Leiden 0.177 | `a3_matched_numbers.json` (typed into the generator) |

## Minor

- **AMI on matched candidates (l.35 to 36).** The text gives the finest-stable AMI, where the
  order matches the ARI. Under the plain maximum the AMI puts Leiden ahead, 0.592 against 0.575,
  where the ARI has them level. The text says nothing false, but a reader may assume the AMI
  agrees under both readings. These AMI means appear in no table (the same flag as for
  Assignment 3). Suggest a clause, or a column in Table IV.
- **Appendix tables reached only from the appendix.** Tables `rulesall`, `common`, `transfer`
  and `ami` are cited, but only inside their own appendices, and the body never points to
  Appendices C, D, F, G or H. Section VI-G quotes transfer, default and silhouette figures. It
  could cite Appendices F and G as it already cites the mode-against-size appendix.
- **Stale pointer, `09_appendix.tex` l.24 to 25.** "so the spread across rules quoted in
  Section VI can be read off directly". Section VI no longer quotes a spread across rules.
- **Generator flag.** `main.tex` l.13 says the numbers come from `a4_protocol_v1.py --tables`.
  The argument parser pairs arguments, so `--tables` with nothing after it is dropped and the
  tables, `numbers_a4.tex` and `latex_values` are not written. `--tables 1` works. Worth fixing
  in the script, or the comment, before the Nov 30 rerun.
- **Figure 2 was not compared visually.** The regenerated `a4_seed_means.pdf` differs in bytes
  (21,511 against 20,064), likely the plotting library version on this machine. The numbers
  behind it match.
- **The PDF in Downloads is the Assignment 3 version** (`..._Interpretation-4.pdf`, compiled
  Nov 11). It is not the draft that was checked here.
