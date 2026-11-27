# R-7 number check, paper draft for Nov 30

Rhea, 2026-11-26. Flags only. Nothing in `latex/` or Overleaf was changed. Utshab applies the fixes.

**What was checked.** The `.tex` at `a8d2e64` (branch `rhea/r7-r9-r10-paper-draft`), every section and
every table it inputs, and the Overleaf PDF. Line numbers refer to `latex/sections/*.tex` and
`latex/tables/*.tex` at that commit.

1. **Source chain.** All five generators rerun on macOS with `envs/scminer.yml`: `a4_protocol_v1.py
   --tables`, `a3_report.py`, `a3_matched.py`, `assignment2/a2_report.py`, `large_scale/make_tables.py`.
   The JSONs, `numbers_a4.tex`, `numbers_large.tex` and every table reproduce byte for byte. Seven figure
   files differ only in rendering (newer matplotlib, no visible pixel change). `make_tables.py` took the
   Allen time and memory from the committed `stage_costs.csv`, since the server logs are not in git.
2. **Every value defined and matching.** All 288 uses of `\V{}` (200 keys) in the 29 input files are
   defined. All 1,125 values in `numbers_a4.tex` and all 67 in `numbers_large.tex` equal their JSON. No
   undefined `\ref` or `\cite`. The Overleaf PDF has no `??`, and its abstract gives the Nov 23 Allen time.
3. **Independent recomputation.** `scripts/tools/r7_recompute.py` recomputes 75 values from the run files
   without the generators: 50 Allen values from `results/large_scale/` and 25 label-free values from
   `results/v1/curves`, rebuilding both readings with the definitions in Section IV. All 75 match. The
   rebuilt picks equal the runner's saved picks in 567 of 567 cases. Taken as given: the 713 s cold read
   (PROGRESS, Nov 19) and the pilot's MICA cost (`docs/allen_pilot.md`).
4. **Claims about each number.** No claim is contradicted by the data. "Every seed" claims hold at all ten
   seeds. The 0.02 cut-off for "level" is now stated (`05_design.tex` l.62 to 63) and every use of "level"
   or "tie" on the benchmark datasets meets it. The flags below are about scope and phrasing.
5. **Nov 13 flags.** All six fixed (one needs a word, flag 7 below). See the follow-up section.
6. **Tables and figures.** Every table and figure the paper inputs is cited in the text. `latex/README.md`
   is out of date (last flag).

## Flags to fix

Each is "clarity" (true, but a reader could take it the wrong way) unless marked.

**Allen section and abstract**

1. **`06_results.tex` l.201 to 202, "The whole run ... took 19 minutes."** The 19 minutes runs from the
   manifest to the last Leiden run (1,130 s). The assignment (about 160 s) and collect stages are not in it.
   Suggest: "From reading the manifest to the last Leiden run, with 3 $k$-means seeds and 14 Leiden runs,
   took 19 minutes." At l.204 to 205, "takes about 29" is an estimate. Suggest "would take about 29".
2. **`00_abstract.tex` l.23 to 24, "no single process used more than 18.2 GB".** True, but it can read as the
   machine's whole need. Suggest adding "and the processes running side by side held up to 59 GB together",
   as Section VI-H does.
3. **`06_results.tex` l.224 and l.233, MICA and Leiden "level" on the sample (0.652 against 0.661).** These
   are single runs, while the cut-off in Section V speaks of mean or median scores. Suggest at
   `05_design.tex` l.62: "when their scores, or their mean or median scores, differ by less than 0.02".
4. **`06_results.tex` l.234, "Leiden run on every cell, at its own setting, reaches 0.748."** That setting is
   resolution 0.05, chosen with the labels, not the 0.17 used on the sample. Suggest: "Leiden run on every
   cell, at resolution 0.05, the count closest to 34 there, reaches 0.748."

**Denominators and scope**

5. **`07_discussion.tex` l.66 to 67 (and `00_abstract.tex` l.14), "chooseR's code as written makes no choice
   in most MICA runs".** The count behind "most" is 106 of 190, which pools both designs, every seed and
   PBMC14K. Here it is used to fail check 6 for the matched comparison only, and the matched-only count is
   not in the JSON. Suggest "in 106 of 190 MICA runs, both designs and every seed", or report the matched
   count.
6. **`06_results.tex` l.106, "84 of 95 MICA runs and 66 of 95 Leiden runs".** The 95 are matched runs only,
   but the subsection says just before that it counts both designs together. Suggest "on matched
   candidates, in 84 of 95 MICA runs ...". Same for `checks.tex` l.23 "(all seeds and datasets)": add
   "matched".
7. **`checks.tex` l.21, "49 of 189 finest-stable picks and 106 plain-maximum picks".** Add "of 189" after
   106. The bare 106 also collides with the 106 of 190 in the next row.
8. **`checks.tex` l.18, "matching moved Leiden's mean from 0.650 to 0.454".** Add "under the plain
   maximum" (under finest stable it moves from 0.662 to 0.591).
9. **`07_discussion.tex` l.36 to 38, "7 percent of Leiden's".** That is Leiden on the five smallest datasets
   only. Suggest "against 5 percent on the four larger datasets and 7 percent of Leiden's on the same five".

**Medians presented as single values**

10. **`07_discussion.tex` l.75 to 81 and `00_abstract.tex` l.21 to 22, MICA on Pollen.** "72 percent of the
    settings lie within 0.01 of the best", "tied settings running from 2 to 17" and the ARIs 0.906 and 0.087
    are all medians over the 7 confirmatory runs. Suggest "in the median run, 72 percent of the settings lie
    within 0.01 of the best, and the tied settings run from 2 to 17 clusters", and "the two readings' picks
    score a median ARI of 0.906 and 0.087".
11. **`06_results.tex` l.140, "picks differ by 0.292 per dataset".** This is the mean absolute change in ARI.
    Suggest "the two readings' picks differ in ARI by 0.292 on average per dataset".

**Checks and effects tables**

12. **`checks.tex` l.28, check 9 ("across implementations").** Its example compares 25 with 100 subsamples,
    which is check 6's subject. The Discussion (l.64 to 65) says the comparison "passes check 9 once
    chooseR's failed intervals are repaired", which refers to F1, not shown in the cell. Suggest: "Finest
    stable, matched: MICA ahead under our scorer, 0.647 against 0.591, and under chooseR's code with the
    point repair, 0.679 against 0.591."
13. **`effects.tex`.** All six labels match the criteria verdicts. But F1 and C4, both confirmed, have no row,
    and C4 is a headline in the Conclusion (l.11 to 13). Suggest adding: "MICA ahead under chooseR's code,
    matched | 0.679 against 0.591, at 7 of 7 new seeds | Confirmed (F1)" and "PBMC14K tie, plain maximum |
    within 0.02 at 4 of 5 seeds, three of them exploratory | Confirmed (C4), depends on the setting per
    count". In row 4, "at 7 of 7 seeds" should read "new seeds".

**Smaller points**

14. **`06_results.tex` l.78.** Under the plain maximum, MICA's lead over $k$-means (0.082) falls to 0.001
    without Usoskin. The Discussion (l.29) says so, but this sentence does not. Suggest adding "so under the
    plain maximum the gap rests on Usoskin".
15. **`06_results.tex` l.141, "(0.197 against 0.118 at seeds 0 to 2)".** These are the Assignment 3 values
    (`A3.C2.*`). The protocol rerun of the same seeds gives 0.198 (`C2.mica.s02`), while l.19 to 21 says the
    reruns reproduce exactly. Probably rounding or a small difference in definition. Uncertain. Suggest "as
    the exploratory report recorded".
16. **`06_results.tex` l.254, "75 to 90 seconds".** The Nov 13 check found 74 to 90 s from the provenance
    files. Recheck the lower bound.
17. **`latex/README.md` is out of date.** It lists Figure 1 as `a3_stability_curve.pdf` (now
    `a4_stability_curve.pdf`, from `a4_protocol_v1.py`), omits `checks.tex`, `effects.tex`, `datasets.tex`
    (typed by hand), `large_scale.tex` and `numbers_large.tex` (`make_tables.py`), numbers the tables from
    the old order, and still says "about 20" numbers are typed by hand.

## Follow-up on the Nov 13 check (`docs/r7_number_check_a4.md`)

**Flags.** All six fixed: the PBMC14K per-subsample time (l.253), the two S1 sentences (joined, l.175 to
177), "all of MICA's lead" (l.43), the 0.02 cut-off (`05_design.tex` l.62 to 63), "does not move" (l.63),
and the scopes of the Section VI-E counts (l.98 to 101), except the "of 95" in flag 6 above.

**Typed by hand.** Moved into `\V{}` since Nov 13: "four of five seeds" (`C4.ties`), the exploratory values
and the Assignment 3 column of the criteria table (`A3.*`, from `exploratory_values()`), and the small-data
cost medians. The old "gap 0.021" sentence is gone. Still typed, and correct where rechecked:

| Where | Text | Rechecked |
|---|---|---|
| `06_results.tex` l.66 | Klein picks (2 clusters at 0.501; 4 at 0.973; 5 at 0.885 to 0.888) | yes, from `results/v1/picks` |
| `06_results.tex` l.102 | chooseR's failed settings | consistent with the per-dataset counts, not fully rechecked |
| `06_results.tex` l.171 to 194 | PBMC14K seed counts, 0.993 to 0.999, 0.48 to 0.66, 0.004, "8 at every seed" | "every seed" at l.194 confirmed for S0, S1, S2 |
| `06_results.tex` l.250 to 256 | "under a minute", PBMC14K cost | see flag 16 |
| `07_discussion.tex` l.99 to 101 | "three ties", "one ties and one ... lead" | yes, from C4's per-seed values |
| `02_introduction.tex` l.4 to 6 | Zeisel example | no |
| `06_results.tex` l.14 | 0.809 and 0.614 | no |

## Not checked

"No comparison passes all eleven" checks (`07_discussion.tex` l.63). The stability curves themselves, which
would need the 11 GB of subsample partitions on the workstation.
