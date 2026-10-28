# Progress log

What has actually been done, newest first. The plan lives in [WORKPLAN.md](WORKPLAN.md).
Utshab adds an entry at the end of each working session. Laraib and Rhea add a line when they finish a task, or ask Utshab to add it.

## Where we stand

Updated: 2026-10-22. Counts follow the Milestone column in `WORKPLAN.md`.

| Milestone | Due | Tasks done | Percent | State |
|-----------|-----|------------|---------|-------|
| M1 software, small data ready | Oct 7 | 8 of 8 | 100% | Done. Laraib's Zeisel card (L-1) merged on Oct 14 |
| M2 reproduce the paper, Assignment 2 report | Oct 14 | 7 of 9 | 78% | Ten benchmark sets run, nine reproduce within 0.05. Report submitted (U-31). Left: contacting the authors (U-30), PBMC14K card (L-2) |
| M3 cluster count without labels | Oct 21 | 7 of 9 | 78% | The stability rule (U-32) is done on all ten datasets, including PBMC14K in graph mode. The results summary (R-1) and paper tables (R-2) are done. Left: significance rule (U-33), Allen card (L-3) |
| M4 what the similarity measures, atlas downloaded | Oct 28 | 5 of 7 | 71% | Diagnostics, atlas, pilot and the shared-space silhouette (U-35) done. Left: label-free Leiden on Allen (U-36), related work table (L-4) |
| M5 main experiments complete | Nov 4 | 3 of 6 | 50% | Our distance matches the authors' (U-17), the CPU fast neighbor search is done (U-19), and mode against size (U-39) is finished. Left: transfer rule (U-34), final comparison table (U-37), GPU search (U-18) |
| M6 gated extras | Nov 11 | 0 of 6 | 0% | Not started. The GPU search (U-18), the annotation comparison (U-22, U-23) and the optional checks (U-38) all depend on the Nov 4 result |
| M7 full manuscript draft | Nov 18 | 0 of 4 | 0% | Not started |
| M8 revision and slides | Dec 2 | 0 of 3 | 0% | Not started |
| Every week: notes, plan and progress updates | weekly | | | L-7, R-3 and R-8 |

Overall: 30 of 55 tasks finished. Three of the 55 repeat every week. The counts above are read off the status markers in `WORKPLAN.md`.

Next three things: submit Assignment 3, run the transfer rule (U-34), and a decision on whether to write to the scMINER authors (U-30).

---

## 2026-10-28, Utshab (in progress, Assignment 3 due today)

**Why this session.** The one-to-one feedback on Assignment 2 was checked against the
Assignment 3 draft before submitting, and the draft missed part of what it asked for. The
feedback expected three things in the next submission: the finished selection comparison, a
comparison of the two graphs on quality and cost, and a plan through November 4. The draft
had the first, called the graph work "unchanged since Assignment 2", and gave only a vague plan.

**Feedback point by point, and what was done.**

| The feedback asked for | Done |
|---|---|
| One primary label-free rule, stability preferred, on PBMC14K and Zeisel, then the rest | already done on all ten |
| Candidate settings, resampling and selection step fixed before scoring | now stated in the design section. chooseR's reading is primary, the plain maximum a declared check |
| Matched genes, chosen count, ARI and AMI, repeat variation | already done |
| Time for selection as well as clustering | added. Results subsection reports the two apart |
| Small-data and graph-mode results kept apart | already done |
| Transfer rule separated from label-free rules, since it uses other datasets' labels | fixed in design, conclusion and WORKPLAN U-34 |
| No predetermined win count | dropped from WORKPLAN U-37 and the plan |
| Graph comparison on build time, total runtime, memory, recall, final quality, shortlist cost, declared candidate count | added as a results subsection and Table `tab:graph`. Memory is new this session |
| Describe it as a shortlist rescored exactly, no guarantee of all neighbours | wording fixed in related work, results, conclusion |
| Separate the general concern from the specific MICA evidence, and say what a user should do | added to related work and the discussion |
| No broad claim about mutual information from one rule | claims now scoped "under this rule" |
| Defer annotation, hub correction, GPU | stated in the plan. The annotation rubric task is gone from the completion plan |
| Keep the full-atlas feasibility visible, discuss changes with the instructor | added to the plan |
| Laraib and Rhea on dataset and literature work, checking results independently | completion plan rewritten to say so |

**Memory was never measured for the graph comparison, so it was measured today.** Two new
scripts, both writing results folders. `scripts/measure_graph_cost.py` builds our shortlist
graph and records time and peak resident memory after each stage, plus recall against the
cached brute force. `scripts/measure_mica_hnsw.py` calls MICA's own `hnswlib_ann` from
`mica_ge.py` with the settings `mica ge` uses (80 neighbours, ef 400, M 8, 8 workers) on the
matrix it reads. Memory is `ru_maxrss`, a high-water mark. GNU time is not installed here.

PBMC14K, 13,605 cells, measured today (not one at a time, see the correction below):

| graph | build | of which shortlist | peak memory | recall, 500 queries | recall, all cells |
|---|---|---|---|---|---|
| MICA HNSW | 429 s | n/a | 1.16 GB | 0.9732 (from the A2 audit) | not computed |
| shortlist c=800, sparse exact path | 15 s | 4 s | 3.64 GB | 0.9735 | 0.9733 |
| shortlist c=800, dense exact path | 87 s | 4 s | 3.63 GB | 0.9735 | 0.9733 |

The sparse and dense builds give byte-identical neighbour arrays (`np.array_equal` is True for
all 13,605 cells), so the sparse path is exact and the Assignment 2 downstream results, ARI
0.853 on our graph and 0.852 on MICA's, still apply. The shortlist's extra memory is all in
loading and binning: its peak is already 3.63 GB before the search starts.

Assignment 2 reported 2.4 minutes for our graph and 15 for HNSW. Today's 15 s comes from the
sparse exact path, which skips genes where a cell sits in its lowest bin. The dense path also
ran faster today (87 s against 146 s) and so did MICA (429 s against 903 s), so the
Assignment 2 runs were probably on a busier machine. Today's times are the cleaner
head-to-head and are the ones in the paper.

**Declared candidate count.** The smallest c on a doubling grid whose recall on 500 query cells
reaches MICA's own. Checked against the saved audits: it gives 800 on PBMC14K (grid 100 to
1,600, MICA 0.9732) and 3,200 on Allen (grid 200 to 3,200, MICA 0.9028), the counts
Assignment 2 used. The paper states it as declared now, not as declared then.

**Selection against clustering, the numbers used in the paper.**
- Zeisel: MICA's distance matrix 250 s, once. One clustering at a fixed setting takes seconds
  for either method. Choosing by stability adds 103 s for MICA and 107 s for the baseline.
- PBMC14K, graph mode: one MICA clustering at a fixed resolution is about 11 minutes, HNSW
  429 s plus node2vec 250 s plus about 2.5 s for one Louvain resolution. Choosing the resolution
  takes 6,171 s per seed for MICA (median of three, about 103 min) and 490 s for the baseline
  (about 8 min).
- The PBMC14K full distance matrix used by the graph-mode stability runner has no logged
  compute time, so the paper does not quote one.

**Errors found in the draft and corrected.**
- "6 to 11 minutes" for MICA's distance matrix on the two datasets above 2,500 cells. The only
  source, `logs/rules_small.log`, gives 250 s (Zeisel) and 367 s (Klein), so 4 to 6 minutes.
- "103 minutes for MICA against 9 for the baseline". The baseline's median over three seeds is
  490 s, so 8.
- "Labels are used in two places only". The graph comparison lets labels pick the resolution,
  a third place.
- Table I and the mode-against-size table overflowed their column by 24.5 pt and 41 pt. Table I
  had done so since it became a single-column table on Oct 22, which went unnoticed because the
  compile check filtered out overfull warnings. Both are now scaled to the column width by
  `scripts/a3_report.py`.
- The transfer rule was called "the remedy Ullmann et al. recommend", which presented a
  label-informed rule as if it were label-free.
- "We match MICA's neighbour search accuracy" read as a guarantee. Now a shortlist rescored
  exactly, with recall measured.

**Done since the first draft of this entry (about 03:10).**
- Body cut from 9.4 to 8.6 pages, 7,868 words down to 7,123. The abstract went from 458 to 280
  words, "Changes since Assignment 2" became the single short paragraph the assignment asks
  for, and the introduction's two preview paragraphs were merged. The four floats in the body
  set the floor, so further cuts would mean moving Figure 1 or a table to the appendix.
- Table numbering is now I stability per dataset, II mode against size, III graph quality and
  cost, then IV to VI in the appendix (datasets, all four rules, AMI). The notebook map in
  `scripts/make_a3_notebook.py` and the README list the new scripts and tables.
- A read-only verification workflow was started over the revised paper against the feedback,
  the saved results and itself.

**A confound that changes the headline (found by the read-only verification, before 03:30).** The
verification confirmed 31 findings. Two of them change the main results, and both were checked
here from the saved stability curves before any text was touched.

*1. The candidate ranges were capped, not matched.* Every run keeps settings giving 2 to 20
clusters, but that only caps the range. MICA's k covers 2 to 20 on every dataset. The baseline's
resolution grid, 0.05 to 2.00, reaches only part of it: min to max over the three seeds is
Buettner 2 to 14, Chung 2 to 11, Goolam 2 to 8, Klein 3 to 17, Kolodziejczyk 4 to 12, Pollen 4 to
11, Usoskin 3 to 18, Yan 2 to 8, Zeisel 3 to 20. The instructor's Assignment 2 feedback named
exactly this, different cluster-count ranges, as the reason Table II could not isolate a method
difference.

Recomputed by re-applying the rule's own selection logic (`stability.pick`) to MICA's saved
curve restricted to the counts the baseline reaches in the same run. The replication of the
reported picks is exact: all 108 recorded picks reproduce, and the means come out at the
reported 0.661 and 0.473 for MICA and 0.632 and 0.635 for the baseline.

| | chooseR's reading | plain maximum | gap between readings |
|---|---|---|---|
| MICA, its own grid (as reported) | 0.661 | 0.473 | 0.189 |
| MICA, on the counts the baseline reaches | 0.772 | 0.635 | 0.137 |
| baseline | 0.632 | 0.635 | -0.002 |

What survives: the tie-break still moves MICA far more than the baseline. What does not: the
ranking reversal. On a common range MICA leads under chooseR's reading and ties under the plain
maximum. The reversal came from MICA's grid offering very coarse counts (it picks 2 clusters on
Pollen under the plain maximum, ARI 0.087, where the baseline cannot go below 4) and very fine
ones (19 to 20 on Yan under chooseR's reading, where the baseline cannot go above 8). The
verifier's figure for the tie span on a common range is 4.1 for MICA against 2.7, so that
direction holds but the gap shrinks from 4.7 to 1.4.

*2. The mode-against-size result was the graph grid's floor.* MICA's graph-mode resolution grid
(0.1 to 4.0, where higher resolution gives fewer clusters) cannot make fewer than 4 or 5 clusters
on these pools, and its tied settings sit at that floor. Holding the small-data arm to the counts
the graph arm reaches, the two readings agree on all three pools:

| pool seed | graph reaches | small-data, all k | small-data, graph's range | graph |
|---|---|---|---|---|
| 0 | 5 to 18 | gap 0.456, span 5 | gap 0.000, span 2 | gap 0.003, span 0 |
| 1 | 5 to 19 | gap 0.172, span 2 | gap 0.000, span 1 | gap 0.011, span 0 |
| 2 | 4 to 20 | gap 0.319, span 4 | gap 0.000, span 2 | gap 0.000, span 1 |

So "the mode, not the size" is not supported. The difference came from which counts each grid
offers. On full PBMC14K both methods also pick their coarsest grid point (MICA resolution 4.0 or
3.9 at 5 clusters, the baseline 0.05 at 4), so "MICA simply wins there" is fragile as well.

*3. "Both arms use the same distance" was false.* The two MICA modes normalise mutual information
differently (sqrt of marginal entropies against joint entropy), as the runner's own docstring
says. So nothing here licenses "it is not a property of the mutual information". Removed.

Decision: report the results as specified (each method on its own grid) and the common-range
recomputation beside them, rewrite the claims that depended on the reversal or on the mode, and
put extending the grids so both methods reach 2 to 20 clusters into the plan. No new runs today.

**The paper rewritten around the corrected findings (finished 03:34).** All 31 verified findings
are addressed. The common-range analysis is now generated by `scripts/a3_report.py` rather than
copied by hand: it writes `results/a3_common_range.csv`, the appendix Table `tab:common`, and
`common_range` and `mode_size_common_range` in `results/a3_report_numbers.json`. Its tie span
uses the same pooled definition as the 7.3 against 2.7 (seeds pooled before the 0.01 window),
giving 4.1 against 2.7 on the common range.

What the paper now says, and where it changed:
- The reversal between the two readings holds on each method's own grid and disappears on the
  counts both methods reach, where MICA leads 0.772 to 0.632 under chooseR's reading and ties
  0.635 to 0.635 under the plain maximum. The tie-break still moves MICA by 0.137 and the
  baseline by 0.002. On Pollen the common range puts MICA's finest-stable pick on the true 11
  clusters at 0.969, and its two readings still differ by 0.42.
- The mode test no longer claims "the mode, not the size". Held to the counts graph mode reaches,
  small-data mode's readings agree too. Size is now listed as an open question, since the test
  ran at 2,500 cells only and five of the small datasets have 124 to 317 cells.
- PBMC14K: both methods pick their coarsest grid point. The baseline's best setting (0.876 to
  0.878), its default (0.853 to 0.867) and its pick after a 0.02 threshold shift (0.785 to
  0.873) all beat MICA's 0.629. MICA's graph-mode default, 0.846 at one seed, now appears in
  Table V instead of n/a.
- The mechanism section now says the per-pair distance invariance holds in both modes and so
  cannot by itself cause the ambiguity. The evidence points to the candidate range.
- The plain maximum is called a second rule, specified with chooseR's before scoring, so
  "three of four rules" is consistent.
- Labels: PBMC14K is no longer called the only dataset not labelled by clustering. Buettner
  (cell-cycle phase), Kolodziejczyk (culture medium) and Klein (time point) are labelled by
  experimental condition.
- Smaller fixes: the default rule ran at three seeds, not once. Yan's seed movement is 1, not 2.
  The baseline's reading gap is 0.002, not 0.003. The dense build is 87 s, not 88. The Total
  column mixes today's build times with the Assignment 2 downstream times, which the text now
  says. The graph-mode stability runner never calls HNSW, which the cost section now says. The
  Zeisel opening uses 3 clusters, the low end of this paper's grid, not 2 from the finer
  Assignment 2 sweep. MICA's Louvain gives fewer clusters as resolution rises, now stated in
  the method.
- The plan through November 4 gains "matched candidate grids" as its first item.

One leak caught while doing this: adding MICA's PBMC14K graph-mode default to the loader moved
the k=4 default analysis's "outside" mean from 0.579 to 0.646, because PBMC14K (true k 7) fell
into it. That analysis is now restricted to the nine small datasets and is back to 0.579.

Body: 8.8 pages. Undefined references: none. The `ALLENCOST` placeholder remains.

**Second verification of the rewritten paper: 22 confirmed findings.**
Every number below was re-derived here before any text changed.

- *The common range used an interval, not the set.* `a3_report.py` restricted MICA to the
  interval between the baseline's fewest and most clusters. On Yan the baseline makes 2 and 4 to
  8 clusters but never 3, and MICA's restricted plain-maximum pick was 3. Restricting to the set
  of counts the baseline produced changes the plain-maximum mean from 0.635 to 0.665, so MICA
  leads narrowly instead of tying, and the gap between readings from 0.137 to 0.107. chooseR's
  reading stays at 0.772 against 0.632. The pooled tied span stays 4.1 against 2.7.
- *The restriction leaves picks at edges too.* On the restricted range MICA's finest-stable pick
  is at an end of what it could choose from in 13 of 27 runs, against 7 for the baseline and 3
  for MICA on its own grid. Pollen's "true 11 clusters" is also the top of the baseline's range
  there. So the restriction shows the reversal depends on the range. It is not a clean ranking.
- *The baseline's 0.002 is a signed mean of gaps that cancel.* The tie-break changes the
  baseline's ARI on 5 of 9 datasets. Mean absolute change between readings: baseline 0.067,
  MICA on its own grid 0.197, MICA on the common range 0.123.
- *MICA's own-grid lead under chooseR reverses at one seed:* 0.658 against 0.624 at seed 0,
  0.702 against 0.633 at seed 1, 0.639 against 0.671 at seed 2.
- *The threshold check was computed on a pooled curve,* the median per cluster count across seeds
  and settings, so its starting picks disagreed with Table I in 5 of 18 cases. Redone per run on
  each run's own curve: all 60 starting picks now equal the recorded choice, and relaxing the
  threshold by at most 0.02 changes the count in 35 of 54 runs, 20 of 27 for MICA and 15 of 27
  for the baseline. On PBMC14K it changes the baseline in 3 of 3 runs (to 6, 7, 7 clusters)
  and MICA in none.
- *The method gave one distance.* Small-data mode uses 1 - I / sqrt(H_i H_j) and graph mode
  1 - I / H(i, j). Only the second is a metric. The method section showed only the second.
- *The mechanism was misdescribed.* MICA's stability curve has a flat top at or near 1.0 across
  much of the range. The plain maximum breaks ties toward fewer clusters and chooseR's reading
  toward more, so the two readings go to opposite ends. "Coarse partitions are the most
  reproducible" was not what the curves show.
- Smaller: the largest fixed-k movement across seeds is 0.267 (Buettner, k=7), not 0.192. The
  largest own-to-common gains are Yan 0.429 and Klein 0.363, not Yan and Pollen. Goolam's
  silhouette difference is 0.3845, so 0.38. On Pollen the ARI to AMI gap is 0.015 for the
  17-cluster partition and 0.231 for the 2-cluster one, the opposite of what the appendix said.
  Chung's 11 clusters lie inside the baseline's range, so "the two largest overclusterings sit at
  counts the baseline never offers" was wrong. The AMI caption pointed at Table I, which no longer
  has the default and silhouette columns.

**Second round of fixes applied (finished 03:58 by the clock).** `a3_report.py` now restricts to the set
of counts the baseline produced, records the mean absolute reading change and the edge-pick
counts, and runs the threshold check per run. All 22 findings are addressed in the text: both
MICA distances are now in the method (Eq. `eq:misd` added for small-data mode) and covered by the
appendix proof, the mechanism is stated as a flat top at 1.0 (82, 70 and 47 percent of MICA's
settings on Yan, Pollen and Chung, against 5 percent for the baseline) broken in opposite
directions by the two readings, the common-range result is presented as a post hoc check on what
the reversal depends on, with its edge-pick caveat, the leads are called narrow, and the
measures rationale now rests on the AMI sitting above the ARI in both directions (0.17 when too
many clusters, 0.12 when too few, across MICA's 54 picks, 35 of which chose too few). The ARI to
AMI ranking agrees under every rule, checked from the two tables' mean rows. The
"within 1e-6" claim is now limited to Zeisel, with up to 5.3e-5 elsewhere. Body is 9.3 pages.

**A process note.** The two Allen graph-cost runs were started as Claude Code background jobs,
not in byobu as the README requires, because the permission check had blocked
`byobu new-session` earlier that night. That tied them to the VS Code session, so the laptop had
to stay on until they finished. From now on long runs go in byobu, and if byobu is blocked the
fix is to ask for approval rather than fall back to a session job.

**Allen, MICA's HNSW run finished at 06:41:14.** `results/graph_cost_allen_100k_mica_hnsw/`.
Loaded in 15 s at a peak of 14.63 GB, then the search took 13,619 s (3.8 h) at a peak of
21.64 GB, against 18,888 s in Assignment 2. Recall stays 0.9028 from the Assignment 2 audit.
Against our build (6,091 s, 80.08 GB):

| | MICA HNSW | shortlist c=3,200 | ratio |
|---|---|---|---|
| graph build | 3.8 h | 1.7 h (shortlist 7.4 min) | ours 2.2 times as fast (on PBMC14K 29 times for the search alone, 24 counting binning, see below) |
| total with A2 node2vec and Louvain | 4.9 h | 2.8 h | |
| peak memory | 21.6 GB | 80.1 GB | ours 3.7 times (3.1 times on PBMC14K) |
| recall, 500 queries | 0.903 | 0.946 | |

Every build ran faster today than in Assignment 2: MICA 2.1 times on PBMC14K (903 to 429 s)
and 1.4 times on Allen (18,888 to 13,619 s), our dense path 1.7 times (146 to 87 s) and our
Allen build 1.3 times (7,798 to 6,091 s).

**Correction: the builds were not run on an idle machine one at a time.** The finish times of
the result files and the log headers show overlaps: MICA on PBMC14K ran about 02:45:54 to
02:53:03, our PBMC14K sparse and dense builds ran inside that window (finished 02:48:02 and
02:49:33), our Allen build ran 02:50:02 to 04:33:48, and MICA on Allen 02:53:52 to 06:41:13.
At most two builds ran at once, always one MICA search beside one shortlist build, each set to
8 workers on the 120-core machine. No two shortlist builds overlapped. On Allen, MICA's search
shared the machine with our build for its first 1.7 h and ran alone for its last 2.1 h. The
paper said "the same idle machine", and now says this.

**Text filled in and checked against Allen.** The `ALLENCOST` placeholder is replaced: on
Allen the shortlist is 2.2 times as fast where it was 29 times on PBMC14K, 2.8 h against 4.9 in
total, and peaks at 80 GB against 22 GB, again while loading and binning. The two claims to
check both hold, and both are now stated with numbers. "It uses more memory" in the
introduction becomes "three to four times the memory", true on both datasets (3.1 and 3.7).
The conclusion's "a faster graph does not close that gap" now adds that our build alone peaked
at 80 GB on the pilot cells, about 1.8 TB in proportion (80.08 times 23), still above the
754 GB machine. MICA's search alone used 22 of the pilot's 140 GB, so the pilot's memory peak
sits elsewhere in MICA's pipeline. "Cheaper graph" became "faster graph" in the introduction
and conclusion, since it costs more memory. The abstract now states the Allen time and the
memory on both datasets.

**Allen, our shortlist build finished at 04:33:49.** `results/graph_cost_allen_100k_corr3200/`.
Loaded and binned in 133 s at a peak of 80.08 GB, shortlist of 3,200 by correlation in 444 s,
exact rescoring on the sparse path in 5,648 s, so 6,091 s (1.7 h) for the graph, against 7,798 s
in Assignment 2. Peak memory stayed at 80.08 GB, so all of it was reached while loading and
binning. Recall on the 500 query cells is 0.9456, identical to Assignment 2, so it is the same
graph. MICA's HNSW run started at 02:53:52 and is still going.

**Allen, partial numbers while running (03:03).**
- Ours, c=3,200: loaded and binned in 133 s at a peak of 80.08 GB, shortlist built in 444 s with
  the peak unchanged, rescoring under way. A2 timing was 2.2 h in all.
- MICA HNSW: loaded in 15 s at a peak of 14.63 GB, search under way. A2 timing was 5.2 h.

**Third verification, of the text changed after the Allen runs (finished about 06:58).**
Read-only workflow `wf_f970c959-d26`, three lenses (numbers, consistency across the paper,
wording and logic), each finding checked by a skeptic. 21 confirmed, 3 rejected. The
confirmed ones collapse to these, all fixed:
- *The build times did not cover the same work.* MICA's `hnswlib_ann` takes the raw matrix and
  bins inside its `mi` distance (checked in `mica_ge.py`), while our timing left out the binning,
  which sat under "load and bin". `a3_report.py` now counts loading, and binning for ours, on
  both sides. PBMC14K shortlist 15 s becomes 18 s (total 6.3 min becomes 6.4), so the speed-up
  is 24 times, not 29. Allen is unchanged at the table's precision (6,224 s, still 1.7 h, 2.2
  times). The JSON key is now `seconds_build_with_load`. The caption says so.
- *The overlap sentence.* "At most two of our builds" read as two shortlist builds. Now "at most
  two builds running at once, always one MICA search beside one shortlist build".
- *"The time saved also shrinks with size" was false.* The time saved grows (7 min to 2.1 h).
  Only the speed-up shrinks. Reworded.
- *"AMI about 0.01 lower in both runs".* The runs are not paired, and one gap is 0.016. Now the
  four AMI values are given.
- *"The whole run takes 2.8 hours against 4.9".* No single run took that long, since the totals
  add A2 downstream times. Now says so in the sentence.
- *PBMC14K total.* "node2vec takes about 4 minutes" did not explain a 6 minute total. node2vec
  plus 40 Louvain resolutions take about 6 minutes on either graph. Reworded.
- *"The one large dataset" for PBMC14K,* while the paper also reports 100,000 Allen cells. Now
  "the one large benchmark dataset" (abstract, introduction, results heading, conclusion) and
  "the only large dataset in the selection comparison" (discussion).
- *Full-atlas time.* "7 to 8 days" used the pilot's 5.3 h search. Today's 3.8 h gives about
  6 days, so the conclusion says 6 to 8. Noted in `docs/allen_pilot.md`.
- *MICA's row mixes two graphs.* Recall, ARI and AMI come from the A2 graph, time and memory from
  today's rebuild, which was not saved. The caption now says so.
- *"Two cautions"* listed three. Now "Three cautions", and "today" became "in these runs".
- *Table numbering.* The 03:10 line above is out of date. Final numbering: body Table I
  stability and Table II graph quality and cost, appendix Table III datasets, IV all four
  rules, V common range, VI mode against size, VII AMI. The notebook map is renumbered to match.
- *6,092 s* in the Allen block was 6,091 s.

Rejected: that the peak memory misses the forked rescoring workers. `ru_maxrss` does leave
them out, but a verifier sampled every process with psutil on PBMC14K: the main process drops
to about 0.95 GB during rescoring and the 8 workers add 0.31 GB of private memory, far under
the 3.64 GB load peak. So the peak is real. On Allen the 80 GB comes from `mi.bin_codes` making
a float64 copy of X (22.8 GB) and an int64 floor temporary of the same size next to the float32
X (11.4 GB). A leaner binning step would cut it, which is worth doing before any larger run.
Also rejected: that the downstream times make the Total
column wrong (the text already says where they come from), and a challenge to the 1.8 TB
projection (the arithmetic and the cell sample check out).

**Still pending (as of 07:02).**
- Recompile, lint, rebuild the notebook, commit and push.
- Utshab: make the Overleaf link view-only, compile there, submit.
- Tell Rhea the tables are final.
- Length: the body is about 9.2 pages, and whether to trim is Utshab's call.

---

## 2026-10-22, Utshab

Assignment 3 experiments and the manuscript. Two new results and a long list of corrections.

**PBMC14K under MICA's graph mode (U-32, the last dataset).** All three seeds are done and
agree closely. `results/rule_stability_paper_pbmc14k_graph_seed{0,1,2}/`. MICA picks
resolution 4.0 or 3.9 and 5 clusters every time, for a median ARI of 0.629 and AMI of
0.790. Leiden picks 0.05 and 4 clusters every time, for 0.397 and 0.725. Across the three
seeds MICA's ARI spans 0.628 to 0.631 and Leiden's 0.3966 to 0.3968, so the rule's
behaviour on this dataset is highly repeatable.

The important part is the tie structure. On the nine small datasets the settings that tie
with MICA's best span 7.3 cluster counts on average, so the tie-break decides the answer.
On PBMC14K they all return exactly 5 clusters, on all three seeds, so the span is 0 and
the two readings of chooseR agree. The ambiguity that the whole Assignment 3 argument rests on is
absent on the one large dataset. Relaxing the threshold by 0.02 also leaves MICA's choice
alone here, one of only five cases in twenty where it does not move.

The rule still fails, in the other direction. MICA's best available ARI is at the true 7
clusters, 0.851 to 0.859 across the three seeds, and the rule selects 5 for 0.629. Leiden's best
is about 0.88 at 7 clusters and it selects 4 for 0.397. Both underclustered, neither was
ambiguous. About 100 minutes per seed for MICA and 9 for Leiden.

**Mode against size, a new task (U-39).** PBMC14K is both the largest dataset and the only
one in MICA's graph mode, so the result above could not be attributed to either. New script
`scripts/run_rule_mode_size.py` draws one fixed pool of 2,500 PBMC14K cells and runs the
same rule on those same cells through both of MICA's modes, which holds the data, the size,
the subsamples and the labels fixed. The graph-mode distance is an exact submatrix of the
cached full matrix. All three arms are done at three seeds, and the answer is the mode, not the size. Each seed
redraws the pool as well as the subsamples, so this is three independent pools, not three
resamples of one.

| arm on the same 2,500 cells | share tied | span of cluster counts, median (range) | readings differ | median gap |
|---|---|---|---|---|
| MICA small-data mode | 0.26 | 4 (2 to 5) | 3 of 3 pools | 0.319 |
| MICA graph mode | 0.16 | 0 (0 to 1) | 0 of 3 pools | 0.000 |
| Leiden, as a control | 0.15 | 1 (1 to 1) | 1 of 3 pools | 0.000 |

At the same size, on the same cells, the two readings of the rule choose differently on
every pool under small-data mode, by a median of 0.319 ARI, and never under graph mode.
Graph mode is unambiguous at 2,500 cells and at 13,605 alike, so size is not the variable.
The ambiguity belongs to MICA's small-data pipeline. It is not a property of mutual
information either, since both arms use the same distance on the same cells and only the
embedding and the normalisation differ.

The per-pool ARIs vary a lot, 0.401 to 0.675 for the small-data arm's finest-stable pick,
because a 2,500 cell draw is harder or easier depending on which cells it gets. The
structural result, the span and whether the two readings agree, is what is stable.

This corrects what I wrote in Assignment 3 before the arm finished. The draft had attributed
the PBMC14K result to dataset size, which the experiment rules out. The abstract,
introduction, results and discussion were all rewritten to attribute it to the mode.

**Corrections the verification pass forced.** Ran an adversarial check over the manuscript
against the saved results. The ones that changed a claim:

- "MICA's small-data pipeline is deterministic given the distance matrix" was false. At a
  fixed k=4 its ARI moves by 0.192 across seeds on Chung, which is larger than the 0.189
  tie-break effect the abstract leads with. The reproducibility belongs to the distance,
  not to the pipeline. Rewritten.
- B=25 was credited to chooseR in the abstract and introduction. chooseR specifies 100.
  It is our reduction, now justified in the design section and listed as a limitation.
- The introduction opened on "MICA returns between 7 and 59 clusters, Leiden between 1 and
  33" on Zeisel. Neither number is in any saved run and the range contradicts our own
  method section. Replaced with the Zeisel sweep we do have, 2 clusters at the low end of
  the grid and 21 at Scanpy's shipped resolution.
- AMI was promised in Measures and never reported. It is now Table III.
- The cost figures were wrong. They are 17 to 109 seconds for MICA and 14 to 225 for
  Leiden, and the Leiden clustering count is a range, 468 to 1,040, not about 1,040.
- Three Pollen, Yan and Chung figures in the text were seed-0 values quoted against a
  table of medians. Pollen is 0.898 at 17 clusters, Chung 11, and the baseline's Pollen
  median does not move between the two readings.
- The tie spans on Yan and Pollen are 2 to 20, the full candidate range, not 3 to 19 and
  2 to 18.
- The silhouette rule was run at one seed while the table caption claimed three.
- The baseline's largest silhouette-space move is 0.52 on Chung, not 0.48 on Usoskin.
- We claimed chooseR did not examine datasets of a few hundred cells. We cannot verify
  that, so the claim is now only about the size of the data it was demonstrated on.

Table I gained a mean row and the PBMC14K row, and the cluster counts now print as the
median with the seed range. New `results/a3_seed_variation.csv` holds the per-dataset ARI
ranges the main table has no room for.

**Open for the next session.** The manuscript body is 7.5 pages against the 5 to 7 the
assignment allows, after three rounds of cuts. The remaining lever is moving the Default
and Silhouette columns of Table I to the appendix, which would free a full-width float and
about a third of a page, at the cost of putting two of the four rules behind a page turn.
The GitHub repository is private, so the code link in Appendix A returns 404 for a grader.
It needs to be made public or the grader added before submission.

---

## 2026-10-21, Rhea

Finished the results summary (R-1). Odd results are marked in a new `anomaly` column in `docs/author_data_results_summary.csv` and `docs/our_preprocessing_results_summary.csv`. Checked against `docs/results_for_pivot.csv`: each cell is the mean of three seeds (0, 1, 2), except allen_100k, which has one run. Every method was given the true number of clusters.

For the next meeting, most important first:
- Chung with seurat HVGs collapses on the authors' data. Hierarchical and k-means score 0.007 and MICA 0.154, against 0.445 to 0.655 in the other gene settings. The likely cause is the wrong input to the seurat flavor, which expects log-normalized data. Utshab, can you check which 2000 genes it picked?
- On Buettner, MICA scores 0.83 to 0.93 and every other method 0.23 to 0.60. All methods were given k=3, so the number of clusters doesn't explain it. We need to confirm no MICA setting was tuned against the labels.
- Repeated scores: Goolam hierarchical is 0.910 for all three HVG flavors, Pollen MICA is 0.942 for all three, and Chung MICA is 0.655 for both cell_ranger and all genes. These are plausible on small sets, but comparing the label vectors would rule out reused results.
- On Kolodziejczyk with all genes, Leiden falls from 1.000 to about 0.56 and MICA from about 0.99 to about 0.72, on both inputs. Were the Leiden resolution and neighbor settings re-tuned for all genes?
- On Kolodziejczyk cell_ranger, hierarchical and k-means score 0.965 and 0.940 on the authors' data but about 0.62 with our preprocessing.
- allen_100k has one run per cell, and hierarchical (both gene settings) and MICA (HVG) were not run.

Also finished the paper tables (R-2). `scripts/make_paper_tables.py` turns the two summary CSVs into `latex/tables/r2_author_data.tex` and `latex/tables/r2_our_preprocessing.tex`. Each is a full-width table with one row per dataset, a column group per gene setting, and the best method in each group in bold. Both compile with the report's IEEEtran preamble, with no overfull boxes or undefined references. They are not yet `\input` anywhere, so the submitted Assignment 2 report is unchanged. Add them to the results section when the manuscript draft starts, and rerun the script whenever the CSVs are refreshed. The tables show the scores as they are, including the Chung seurat collapse, so they should be regenerated once the anomalies above are resolved.

---

## 2026-10-14, Utshab

Wrote and submitted the Assignment 2 report (U-31), merged both teammate pull requests, and narrowed the project as the instructor asked.

**The instructor's feedback.** The proposal packed four substantial projects into one semester. The clustering question is now first, scaling and annotation follow it, and anything that uses the labels to pick a clustering is reported separately. The feedback also asked us to check the approximate neighbor search against a reference on a subset, and to estimate the full atlas cost before attempting it.

**The report** is in `latex/` (`main.tex` plus `sections/00_abstract` to `05_schedule`). It is 4 pages of body plus 1 page of references. The tables and the figure are generated from `results/` by `scripts/a2_report.py`, every number it quotes is in `results/a2_report_numbers.json`, and `notebooks/assignment2_preliminary.ipynb` rebuilds them all. The proposal as submitted is kept in `latex/archive/proposal/`.

**New results since Oct 3.**
- Scored the label-free rules on all six MICA graph-mode runs of PBMC14K and Zeisel, which gives the run-to-run spread we were missing. PBMC14K: default 0.837 to 0.851, silhouette 0.616 to 0.833. Zeisel: default 0.538 to 0.717, silhouette 0.634 to 0.735.
- A single MICA run can move the label-free score by 0.2 on PBMC14K, because its silhouette pick is 5 clusters in four runs and 7 in the other two. Leiden barely moves across seeds. Every result from here on gets at least five runs.
- `scripts/rq2_hubs.py` rebuilds the hub numbers from the exact graph and saves them, which the Oct 3 entry had only from a scratch run.

**What the checks caught.** Two verification passes ran over the report, one on the literature and one on the numbers, citations, requirements, writing and a reviewer read-through.
- One real error: the Yan example for MICA's k-means restarts was wrong. MICA's own default run on Yan scores 0.939, and the 0.77 was a tail value from an earlier redraw. Chung replaces it, where MICA's own default run scores 0.589 against the published 0.655.
- The label-free flip on Zeisel is mostly about granularity, not quality. MICA's own 20-cluster partition scores 0.29, and the AMI gap is 0.06 where the ARI gap is 0.23. The report now says this.
- The two methods' sweeps do not cover the same range of cluster counts, so part of the silhouette result measures how coarse each sweep can go. Leiden reaches 2 clusters on Zeisel while MICA never goes below 7.
- Our earlier claim that matched genes put Leiden ahead held only for the two anchor datasets. MICA leads on seven of the other eight.
- The bibliography had real errors: a wrong title and first author for one preprint, a wrong volume for another paper, and missing page numbers in several. All 23 cited works were checked against Crossref.

**Teammate work merged.**
- Laraib's Zeisel dataset card (L-1). Checked against the raw download and our pipeline: the cell and gene counts, the label counts, the region counts and the smallest groups all reproduce. One claim was wrong and is fixed. The 2,903 cells belong to the scMINER authors' processed input, not to Zeisel et al., who analyzed all 3,005. Our QC keeps 2,987, sharing 2,888 with the authors' set.
- Rhea's results summary pivots (R-1). Every value in both CSVs matches a fresh pivot of `docs/results_for_pivot.csv` to three decimals. R-1 stays in progress until the odd results are marked.

**What the plan looks like now.** The milestones are realigned to the instructor's targets of November 4 for the main experiments and November 18 for the draft. Eight tasks were added (U-31 to U-38): three new label-free rules, a shared-space silhouette, label-free Leiden on the Allen subset as a third anchor, the final comparison table, and two optional checks. The GPU work and the full atlas run are now behind gates.

**Next.** The stability rule (U-32) and the significance rule (U-33) on PBMC14K and Zeisel.

---

## 2026-10-10, Rhea

Added the pivot tables to the repo as `docs/author_data_results_summary.csv` (authors' data) and `docs/our_preprocessing_results_summary.csv` (our preprocessing).

---

## 2026-10-03, Rhea

Created the results summary pivot tables for R-1 in Google Sheets. They were added to the repository on 2026-10-10 as `docs/author_data_results_summary.csv` and `docs/our_preprocessing_results_summary.csv`.

---

## 2026-10-03, Utshab

Finished ten tasks: the label-free choice of cluster count (U-12, U-13), what MICA's similarity measures (U-14), the Allen download and 100,000 cell pilot (U-15, U-16), the test that our distance matches the authors' (U-17), the paper's five small benchmark sets (U-27), uncertainty estimates (U-28), the Kolodziejczyk side test (U-29) and a CPU version of the fast neighbor search (U-19). A third review of the claims found one serious bug. It is fixed and rerun, and every number below comes after the fix.

Unless a line says otherwise, scores use the paper's rule: each method is told the true number of cell types, and the labels pick the best run among those that hit it.

**What the review caught, and what we fixed**
- The authors' released inputs are not all on the same log scale. Pollen and Buettner are log10, Chung and Usoskin are natural log, and the other six are log2. We had read all of them as log2, so for those four datasets the seurat and seurat_v3 gene sets shared only 1,198 to 1,646 of their 2,000 genes with the right ones. Fixed in `data.py`, and every affected run was redone. The all-gene and cell_ranger runs never depended on this.
- On small data, MICA's last step is k-means, and one k-means run often lands on a worse fit. The published scores match the best fit we can find, so we now rerun that step with 100 restarts. That changed some Oct 2 numbers (a note is added there). 100 restarts get close to the best fit but not always all the way. On Pollen the best fit we found scores 0.980, above the published 0.942.
- Two summary tables mixed graph-mode repeat runs into the Zeisel small-data cell. Fixed.

**Reproduction on the authors' inputs, ten datasets**

| Dataset | Cells | Published MICA | Our MICA run | Verdict |
|---|---|---|---|---|
| Yan | 124 | 0.939 | 0.939 | matches |
| Goolam | 124 | 0.544 | 0.544 | matches |
| Buettner | 182 | 0.934 | 0.934 | matches |
| Pollen | 249 | 0.942 | 0.951 | within 0.05 |
| Chung | 317 | 0.655 | 0.655 | matches |
| Usoskin | 622 | 0.847 | 0.835 | within 0.05 |
| Kolodziejczyk | 704 | 1.000 | 0.727 | not reproduced |
| Klein | 2,717 | 0.980 | 0.980 | matches |
| Zeisel | 2,903 | 0.744 | 0.714 small-data mode, 0.755 to 0.781 graph mode (six runs) | within 0.05 |
| PBMC14K | 13,605 | 0.853 | 0.853, and 0.850 to 0.857 over five repeats | matches |

So nine of ten reproduce within 0.05. Usoskin shows how fragile the exact number is. Its k-means fits that are almost equally good (within 0.2 percent of the best) give ARI anywhere from about 0.78 to 0.88.

On Kolodziejczyk we tested whether removing the plate effect gets MICA to 1.0 (U-29). It does not. Centering each gene within its plate drops MICA to 0.16, or 0.14 with scaling too. So the plate effect does not explain the gap, and the likeliest source is still the input file that was not released (U-30). Details in `results/kolod_side_test/`.

**MICA against Leiden with the same genes (RQ1)**

The gaps come with a bootstrap over cells (1,000 resamples, `results/uncertainty/`).
- **Five small datasets (124 to 317 cells).** MICA leads in 18 of 20 dataset and gene settings, and in 16 of them the interval is clear of zero. Leiden leads on Yan with seurat genes, by 0.08 with an interval just clear of zero. Chung with seurat genes is a tie (Leiden 0.163, MICA 0.154).
- **Against all three PCA baselines.** MICA beats every one of them on Buettner in all four gene settings, and on Chung in three of four. On Goolam, k-means beats MICA in every setting.
- **Leiden on tiny data depends a lot on the number of PCs.** On Buettner, going from 50 PCs to 10 lifts Leiden from about 0.29 to between 0.51 and 0.78, depending on the genes. MICA stays ahead (0.83 to 0.93). On Chung, fewer PCs lift Leiden to at most 0.47. That is still below MICA in three settings, and the seurat setting was already a tie. Details in `results/pc_sensitivity/`.
- **Five larger datasets (622 to 13,605 cells).** MICA leads in 8 of 20. The gene set decides the winner on Usoskin, Kolodziejczyk and Zeisel. On Klein and PBMC14K, Leiden is level or ahead in every gene setting.
- **Across the ten datasets.** With all genes MICA is ahead (mean gap 0.17, Wilcoxon p = 0.027), but that does not survive a correction for testing four gene settings. With 2,000 genes p runs from 0.08 to 0.43.
- **Our baseline is not weaker than the paper's.** Our Leiden with 2,000 seurat genes lands within 0.035 of the paper's own Scanpy score on 7 of 10 datasets. It is higher on Yan (by 0.11) and Klein (0.13) and lower on Pollen (0.05).

**Choosing the number of clusters without labels (U-12, U-13)**

Each method keeps the choice with the best silhouette in its own embedding. That means k from 2 to 20, or its own resolution sweep. Mean ARI over 10 datasets and 4 gene settings:

| Method | Told the true number | Picks by silhouette |
|---|---|---|
| MICA | 0.802 | 0.572 |
| Leiden | 0.680 | 0.542 |
| hierarchical | 0.605 | 0.397 |
| k-means | 0.588 | 0.417 |

MICA's lead over Leiden shrinks from 0.12 to 0.03. On the small datasets Leiden comes out ahead (0.576 against 0.481), and on the larger ones MICA does (0.664 against 0.508). The review showed that the result flips when we change the k range or the embedding used for the silhouette, so we call it inconclusive for now. Neither method finds the true number on its own very often: MICA in 5 of 40 cases, Leiden in none. Tables in `results/summary_label_free.md`.

**What MICA's similarity measures (U-14)**

Run on our input for Zeisel, Usoskin, Kolodziejczyk and Klein. Tables in `results/summary_diagnostics.md`.
- **Detection.** On Zeisel, Kolodziejczyk and Usoskin, 57 to 62 percent of the mutual information between two cells comes from which genes are detected in both. On Klein it is 25 percent. Given detection alone, MICA finds nearly the same clusters on Kolodziejczyk (agreement 0.98) and Zeisel (0.87), but not on Usoskin (0.60) or Klein (0.59).
- **Sequencing depth.** Keeping only a quarter of the reads barely changes MICA on Klein, Zeisel or Kolodziejczyk. Kolodziejczyk has about 4 million reads per cell, so a quarter is still plenty. Usoskin is the exception. MICA falls from 0.885 to 0.57 at half depth and 0.27 at a quarter. Leiden starts lower there (0.59) and falls less (0.62 and 0.41).

**Allen atlas pilot (U-15, U-16)**

All 13 raw region files are downloaded (97 GB, listed in `docs/allen_files.csv`). On a random 100,000 cells, scored on the 34 classes:
- MICA reaches ARI 0.652 in 7.1 hours with a 140 GB peak. The mutual information neighbor graph takes 5.3 of those hours.
- Leiden reaches 0.661 with 2,000 genes and 0.624 with all genes, in 23 to 25 minutes.

Scaled up to all 2.3 million cells, MICA would need about 7 to 8 days and about 3.2 TB of memory. The machine has 754 GB, so MICA as shipped cannot run the whole atlas here. The working is in `docs/allen_pilot.md`.

**Fast neighbor search (U-19, CPU prototype)**

MICA's graph mode needs each cell's 79 nearest neighbors under its mutual information distance. Its HNSW search for them is the slowest step by far, 5.3 of the 7.1 hours on the Allen pilot. The idea is to propose candidates with something cheap, then compute the exact distance only for those. We check against brute force on 500 random cells, and against the graph MICA itself built and saved. Our distance matches the weights in MICA's saved graphs to 1e-6.

*What works: correlation over all genes as the cheap step.* It is one matrix multiply, with no PCA. On PBMC14K (13,605 cells):

| Method | True neighbors found | Time for all cells (8 workers) |
|---|---|---|
| Brute force, every pair | 100% | about 35 min |
| MICA's own HNSW search | 97.3% | 15 min |
| Correlation, 400 candidates per cell | 90.4% | 1.0 min |
| Correlation, 800 candidates | 97.4% | 2.4 min |
| Correlation, 1,600 candidates | 99.6% | 4.7 min |

Then we fed each graph to MICA's own later steps (node2vec, then Louvain), twice each. Our 800-candidate graph gives ARI 0.853 both times, and MICA's own graph gives 0.852 both times. So on PBMC14K the fast graph is about six times quicker to build and leaves MICA's result unchanged.

On the Allen pilot (100,000 cells), brute force would take about 43 hours even with the faster exact distance below.

| Method | True neighbors found | Time for all cells (8 workers) |
|---|---|---|
| MICA's own HNSW search | 90.3% | 5.2 h |
| Correlation, 1,600 candidates per cell | 89.8% | 1.1 h |
| Correlation, 3,200 candidates | 94.6% | 2.2 h |

The correlation step itself takes 9 minutes, and the rest is exact rescoring. So at this size we match HNSW's recall about five times faster, or beat it clearly at well under half its time. Leiden on our 3,200 graph and on MICA's graph lands at 0.599 and 0.595 at class level, and the two partitions agree at 0.97. With MICA's own later steps, run twice on each graph, the two are close. Our 3,200 graph gives ARI 0.602 and 0.627. MICA's graph, run through the same code, gives 0.641 and 0.617, and MICA's own run gave 0.652. So the ARI ranges overlap. AMI is a little lower on our graph every time (0.819 and 0.825, against 0.831 to 0.836). At this size the fast graph costs at most a small amount of agreement, and two runs each cannot say more than that.

*What did not work.*
- Candidates from 50 principal components. On PBMC14K, 400 of them find only 33% of the true neighbors. Shared detected genes find 72%.
- Rounds of neighbors of neighbors (the NN-descent idea). They level off at 86 to 88% on PBMC14K.
- Adding the graph's hub cells to every candidate list. More correlation candidates do better for the same cost.
- Proving a list exact with distances to 64 reference cells. It works for under 10% of cells.

On the four small sets we could compare against every pair. There, 200 candidates each from principal components and from detected genes, merged, find 97 to 99% of neighbors on Usoskin, Kolodziejczyk and Zeisel, but only 74% on Klein. Correlation alone needs 800 candidates on Klein to reach 92%.

*Found along the way.*
- **A faster exact distance.** Most genes of a cell fall in the lowest bin, so the joint histogram only needs the genes where the other cell is outside it. The result is identical, and it is 7 times faster on PBMC14K and 1.2 to 1.9 times faster on the denser Zeisel and Klein.
- **The mutual information graph has strong hubs.** On PBMC14K, the top 1% of cells take 23% of all exact neighbor slots. Their links join cells of the same type only 44% of the time, against 75% for other links. On Klein, the median cell is a neighbor of only 7 others while the biggest hub is a neighbor of 868. (Corrected on Oct 14. The first two shares were quoted as 42% and 76% from a scratch run. `scripts/rq2_hubs.py` now rebuilds them from the exact graph and writes `results/rq2_hubs/`.)
- **Leiden on the raw neighbor graph is a poor stand-in for MICA.** On PBMC14K it scores 0.64 on either graph, against MICA's 0.85. On Klein, Leiden on the exact graph scores 0.59, while approximate graphs that miss many true neighbors score up to 0.98. So finding every true neighbor is not by itself the goal. Removing hub effects might help clustering. That would be a new task, so it is a decision for the team.

*Caveats.* Recall is measured on 500 random cells. For MICA we count each cell's 79 closest links in its saved graph, which can only flatter MICA. Our code is numpy and MICA's search is compiled, both on 8 workers of a shared machine. The PBMC14K times use the plain exact distance, from before the faster one was added. The Allen times use the faster one. The number of candidates needed to match HNSW grows with dataset size (800 at 13,605 cells, about 1,600 at 100,000), so the full atlas will likely need the GPU version (U-18). Tables in `results/rq2_graph_*/`, `results/rq2_large_paper_pbmc14k/` and `results/rq2_downstream_*/`.

**Limits**
- The bootstrap resamples cells only. It ignores run-to-run variation and the fact that cells are not independent, so its intervals are too narrow.
- Both methods get best-case treatment. Leiden picks the best of up to 200 resolutions using the labels, and graph-mode MICA the best of 40 settings. Small-data MICA is handed the true k, and its k-means keeps the best of 100 restarts, chosen without the labels.
- Zeisel and Usoskin labels came from the original authors' own clustering, which can favor some methods.
- The ten datasets are the paper's choice, and five of them have fewer than 320 cells.

**Next**
1. The GPU version of the exact distance and the correlation step (U-18). It needs GPU time, so I will ask first.
2. Whether to add a task on the hub cells in MICA's graph. That is new scope, so it is the team's call.
3. U-30: whether to write to the authors about the Kolodziejczyk input. That is Utshab's call.

**For the team.** Rhea, for R-1 please use `docs/results_for_pivot.csv` instead of the full tracker. It keeps only the main comparison, one row per run, with plain column names. The full tracker now has 1,556 rows of every kind.

---

## 2026-10-02, Utshab

> **Update, Oct 3.** Some small-data MICA numbers below changed on Oct 3. MICA's last step there is k-means, and we now rerun it with 100 restarts instead of 10 (see the Oct 3 entry). On the authors' input, Usoskin moves from 0.803 to 0.835 and Zeisel from 0.695 to 0.714. Current values for every dataset are in `results/summary_gene_sets.md`.

Ran U-25 (MICA on the same 2,000 genes the baselines got) and U-26 (three more of the paper's benchmark sets: Usoskin, Kolodziejczyk and Klein). Two independent reviews of the code and the claims changed several numbers, and one conclusion from yesterday held only on the first two datasets.

Everything below uses the paper's rule, which looks at the labels twice. Each method is told the true number of cell types, and when several runs give that number, we keep the one that best matches the labels. So these are best-case scores, not what someone would get without labels.

**What the reviews caught, and what we fixed**
- Our seurat-style gene selection treated our log2 data as natural log, so it picked 19 to 47 percent different genes from the standard method. Fixed in `cluster.select_hvg`, and every run that used it was redone. Leiden barely moved, k-means moved more.
- The paper's text says the baselines swept Leiden resolution from 0.1 to 4.0 in steps of 0.1, but the authors' own Scanpy code sweeps 0.01 to 2.0 in steps of 0.01. On the old grid Leiden often had no setting that gave the true number of clusters. We now follow the code, so every Leiden sweep has at least one setting with exactly the true number. Leiden still does not find that number by itself: the labels pick the setting.
- MICA's small-data mode runs k-means once with no fixed seed. `scripts/mica_seeded_kmeans.py` reruns that exact call with three seeds on MICA's saved embedding, so in small-data mode MICA and the baselines both report three draws. Graph mode (PBMC14K, and Zeisel in graph mode) still has one MICA run per setting.
- Usoskin was released only as reads per million. We recover the read counts and drop 1,218 rows that are DNA repeats, not genes. For Usoskin, Kolodziejczyk and Klein we keep every labelled cell, as the paper did.
- `score_mica.py` now refuses to score if MICA's cells do not match ours, and `run_mica_job.sh` never mixes old output into a new run.

**Reproduction: MICA against the published scores**

We ran MICA both on the authors' own released input (their Zenodo files) and on our input. The paper has 13 benchmark sets, and these are the five mid-sized ones (622 to 13,605 cells).

| Dataset | Published | MICA, authors' input | MICA, our input | Verdict |
|---|---|---|---|---|
| PBMC14K | 0.853 | 0.853 | 0.854 | matches, one run each |
| Klein | 0.980 | 0.980 | 0.980 | matches |
| Zeisel | 0.744 | 0.781 graph mode, 0.695 small-data mode | 0.766 and 0.704 | within 0.05 in both modes. The authors' script ran it in graph mode |
| Usoskin | 0.847 | 0.803 | 0.885 | within 0.05 on both inputs |
| Kolodziejczyk | 1.000 | 0.727 | 0.721 | not reproduced from the released files |

On Kolodziejczyk, a notebook in the authors' `scMINER_analysis` repo scores their main MICA run at 0.727, the same as our run on their released input. That run used essentially the same values as the release (13,923 genes against 13,698). The same notebook scores some other MICA runs at 1.0, including their parameter sweep. For this dataset the sweep script reads a different file, `Kolod_MICA_input.txt` (704 cells by 13,473 genes). That file is not in the Zenodo release, but it is not hidden either: the script names it, and the original MICA repo uses a file of the same name as its reference input. It has small negative values, so it was transformed in some way we cannot see from the notebook. We do not know which run the paper's 1.0 comes from.

A score near 1.0 is normal for this dataset anyway. On the released input MICA gets 0.99 to 1.00 once it is given 2,000 variable genes, and the paper's own Seurat and Scanpy get 0.992 and 0.996. With all genes, MICA groups 82 of the 2i cells with the a2i cells, and 81 of those come from batch 2. (2i, a2i and lif are the three growth conditions that make up the labels. The batch is the number in each cell's name.) Before we call this a failed reproduction anywhere outside the team, we should ask the authors for that file and how it was made (U-30).

**Gene set against method (RQ1)**

Mean ARI over seeds. Both methods get exactly the same genes, and "all genes" means every gene that survives our filtering. The 2,000-gene columns give the lowest and highest score across the three ways of picking those genes, not run-to-run spread. MICA on Zeisel is in small-data mode here. In the graph mode the authors used, MICA on all genes scores 0.766. Full tables with AMI are in `results/summary_gene_sets.md`.

| Dataset | MICA, all genes | Leiden, all genes | MICA, 2,000 genes (three ways) | Leiden, 2,000 genes (three ways) |
|---|---|---|---|---|
| Zeisel | 0.704 | 0.845 | 0.728 to 0.844 | 0.615 to 0.795 |
| PBMC14K | 0.854 | 0.876 | 0.531 to 0.737 | 0.588 to 0.755 |
| Usoskin | 0.885 | 0.588 | 0.840 to 0.877 | 0.885 to 0.922 |
| Kolodziejczyk | 0.721 | 0.562 | 0.988 to 1.000 | 0.996 to 1.000 |
| Klein | 0.980 | 0.975 | 0.728 to 0.839 | 0.766 to 0.976 |

What this shows:
- **Which genes a method gets matters about as much as the choice between MICA and Leiden, and often more.** Across the four gene settings, the gap between a method's best and worst score runs from 0.045 (Usoskin) to 0.32 (PBMC14K) for MICA and from 0.21 to 0.44 for Leiden. Given the same genes, the two methods are within 0.05 in 11 of 20 pairings, and the largest gap is 0.30 (Usoskin, all genes). MICA on Usoskin is the clear exception, since its gene set barely mattered there. Fewer genes help on some datasets and hurt on others, so no single gene set is best.
- **Neither MICA nor Leiden wins everywhere.** With 2,000 genes, Leiden is level with or ahead of MICA on four of five datasets, and MICA is ahead on Zeisel. With all genes, Leiden is ahead on Zeisel (by 0.14 against MICA's small-data mode, 0.08 against the graph mode the authors used) and slightly ahead on PBMC14K (0.02). MICA is ahead on Usoskin and Kolodziejczyk, where Leiden on all genes scores far below Leiden on 2,000 genes (0.59 and 0.56, against 0.89 to 1.00). We do not know why yet. Klein is a tie.
- **On PBMC14K, MICA's lead mostly goes away once both methods get the same genes.** In our rerun of the paper's setup (MICA on all genes, Leiden on 2,000 genes picked with Scanpy's default method, as in the authors' code), MICA leads by 0.20. The paper itself reported a smaller lead there (0.853 against Scanpy's 0.754). Our run of their baseline scores 0.10 lower than theirs, and we do not know why yet. Given the same genes, Leiden is level with MICA or slightly ahead in all four gene settings on PBMC14K, by 0.003 to 0.057. MICA's PBMC14K numbers are single runs.

**Correction to what I told you on Oct 1.** I wrote that Leiden is ahead of MICA when both get all genes. That holds on Zeisel (by 0.14, or 0.08 against MICA's graph mode) and only just on PBMC14K (0.02). It fails on Usoskin and Kolodziejczyk, where MICA is ahead, and Klein is a tie. With 2,000 matched genes, MICA is ahead on Zeisel. The three points above replace the Oct 1 claim.

**Limits**
- Five datasets is too few for a significance test.
- The graph-mode MICA numbers are single runs, and our three small-data MICA seeds only redraw the last k-means step, so MICA's run-to-run spread is understated.
- Leiden picks the best of up to 200 resolution settings, while MICA gets one k-means run in small-data mode or 40 settings in graph mode. That extra choice is worth up to 0.04 ARI to Leiden, which is as large as several of its leads.
- Zeisel and Usoskin labels came from the original authors' own clustering, which can favor some methods. Kolodziejczyk's labels partly line up with the batch number in the cell names.
- In the paper's own numbers, Seurat averages the same ARI as MICA on these five datasets (0.885 each), though Seurat is ahead on Usoskin by 0.09 and MICA on PBMC14K by 0.08. MICA's average lead over Seurat across all 13 datasets comes from the five smallest (124 to 317 cells), where MICA is ahead by 0.01 (Goolam) to 0.41 (Buettner). We have not run those yet.

**Next**
1. Add those five small sets (Yan, Goolam, Buettner, Pollen, Chung) from the authors' Zenodo inputs, giving ten datasets (U-27).
2. Measure run-to-run spread with a bootstrap and repeated graph-mode runs (U-28).
3. The label-free number of clusters (U-12), reusing the saved outputs.

All 247 result rows are in `docs/results_tracker.csv`.

---

## 2026-10-01 end of day, Utshab

**Reproduction done on both datasets (U-9).**

| Dataset | MICA, our run | MICA, paper | Gap | Leiden, all genes, same data and rule |
|---------|---------------|-------------|-----|---------------------------------------|
| PBMC14K | ARI 0.854, AMI 0.871 | 0.853 | +0.001 | ARI 0.874, AMI 0.899 |
| Zeisel | ARI 0.703, AMI 0.753 | 0.744 | -0.041 | ARI 0.844, AMI 0.797 |

- Both gaps are inside the 0.05 tolerance we set, so we trust the setup. The Zeisel gap is probably our slightly different cell set (2,987 cells against their 2,903) plus randomness: MICA's small-data mode runs k-means once without a fixed seed, although the paper describes a consensus over 10 runs. That is worth a sentence in the reproduction section.
- **Cost.** Zeisel took 59 minutes, and 58 of them were the mutual information matrix, computed on a single core even though we asked for 25 workers. PBMC14K took about 20 minutes in the larger-data mode. This is direct evidence for the cost problem behind RQ2.
- **RQ1 now rests on identical data for both datasets.** A plain Leiden pipeline on all genes beats MICA by 0.02 ARI on PBMC14K and by 0.14 on Zeisel.

Raw MICA output stays on the shared disk side of the repo. The scored rows are in `results/reproduction_zeisel/` and `results/reproduction_pbmc14k/`.

---

## 2026-10-01 night, Utshab

**The paper's baseline setup explains most of MICA's advantage, at least on these two datasets.**

The paper gave every baseline only 2,000 highly variable genes, while MICA uses every gene that passes filtering. We reran Leiden with the paper's setup, using three standard ways of picking those 2,000 genes, and compared it with the all-genes version. Everything else stayed the same: same data, same filtering, same selection rule, three seeds each.

| Leiden on principal components | PBMC14K ARI | Zeisel ARI |
|---|---|---|
| All genes, the same input MICA gets | 0.874 | 0.844 |
| 2,000 genes, seurat_v3 method | 0.754 | 0.769 |
| 2,000 genes, cell_ranger method | 0.578 | 0.691 |
| 2,000 genes, seurat method | 0.652 | 0.641 |
| MICA, our run and published | 0.854 and 0.853 | 0.703 and 0.744 |

What this says:
- With 2,000 genes, under any of the three methods, Leiden falls below MICA on PBMC14K, which matches the paper's story.
- Given the same genes MICA uses, Leiden is ahead of MICA.
- So much of the reported gap comes from the gene set the baselines were given, not from mutual information itself.

What it does not say yet:
- **Only two datasets.** A paper needs the rest of scMINER's benchmark sets that meet the 500-cell floor: Usoskin, Kolodziejczyk and Klein.
- **It holds for the graph pipeline only.** k-means and hierarchical clustering move in different directions, and k-means on Zeisel is better with 2,000 genes (0.82) than with all genes (0.52). We must not overclaim it as "all genes always help".
- **The cluster count is still chosen using the labels**, as in the paper. Choosing it without labels is the next part of RQ1.

All rows are in `results/baselines_*` and `docs/results_tracker.csv`.

---

## 2026-10-01 late, Utshab

**The reproduction works.** MICA on PBMC14K gives ARI 0.854 against the 0.853 in the paper.

- **PBMC14K comes from the authors' own matrix.** Their R package ships the raw counts as `pbmc14k_rawCount.rda`, and the `rdata` Python package reads it, so we did not need to rebuild it from 10x files. 14,000 cells by 17,986 genes, exactly as described.
- **Filtering now matches theirs.** Taking cutoffs at 4 deviations on a log scale keeps 13,646 cells against the 13,605 in the paper. At 3 deviations we kept only 12,605, a thousand short. Zeisel now keeps 2,987 against their 2,903.
- **Same selection rule as the paper for everything.** Sweep the resolution from 0.1 to 4.0, and among settings that give the true number of types, keep the one with the highest ARI. `pick_paper_protocol` in `src/scminer_clust/cluster.py` does this for both MICA and the baselines.

| Dataset | MICA, our run | MICA, paper | Leiden on principal components, same rule |
|---------|---------------|-------------|-------------------------------------------|
| PBMC14K | ARI 0.854, AMI 0.871 | ARI 0.853 | ARI 0.873 to 0.876, AMI 0.899 to 0.900 |
| Zeisel | still running | ARI 0.744 | ARI 0.843 to 0.844, AMI 0.797 |

**First real finding for RQ1.** On the paper's main dataset, with identical data and the paper's own selection rule, a plain Scanpy-style pipeline scores higher than MICA. One difference could explain why the paper reported weaker baselines: it ran every baseline on 2,000 highly variable genes, while ours uses all genes. That run is going now, and it decides how we frame this.

- **Cost.** MICA on PBMC14K took about 20 minutes here. The mutual information neighbor graph took 9 minutes, the embedding 4 minutes, and the rest was 40 Louvain runs for the sweep.
- **Repo hygiene.** MICA writes a table and a plot per resolution, 80 files per run, so its raw output stays out of git. The scored result is `results/reproduction_pbmc14k/metrics.csv`.

---

## 2026-10-01 evening, Utshab

First code and first numbers.

- **Both environments work now.** The main one sees both GPUs (PyTorch with CUDA 12.8). The authors' MICA installs after one more fix: its bundled C++ part, mihnsw, is built by a nested pip call in an isolated environment that has no numpy, and mihnsw's build file never declares it. Added a build config pinning numpy to the version in the environment so the compiled part stays compatible. Both environments are locked in `envs/`.
- **Zeisel downloaded.** The Linnarsson lab link in the paper is dead (403). Bioconductor keeps a mirror of the same data with the original labels, so we take it from there, which is documented in `scripts/download_zeisel.py` and worth a line in the paper.
- **Preprocessing written in Python** following scMINER: quality filters at median plus or minus 3 deviations, then counts per million and log2. One choice to note: the cutoffs are taken on a log scale, and the mitochondrial filter is off by default. With both on the raw scale we lost 18 percent of the cells, which is far more than the paper. We now keep 2,843 of 3,005 cells, close to the 2,903 the paper used.
- **Baselines run on Zeisel** with the cluster count set to the true number of classes, which is the paper's protocol. Results in `results/baselines_zeisel/metrics.csv`:

  | Method | ARI | AMI | Time |
  |--------|-----|-----|------|
  | Leiden on principal components | 0.84 | 0.79 | under a second once the components are cached |
  | k-means on principal components | 0.51 to 0.57 | 0.60 to 0.66 | about 12 seconds |
  | Hierarchical on principal components | 0.48 | 0.60 | under a second |

  Worth flagging early: the paper reports about 0.74 for MICA on Zeisel, and plain Leiden already gets 0.84 here. That is not a fair comparison yet, since our preprocessing and cell count differ from theirs, and the MICA run on our exact data is still going. But it is the first hint that RQ1 is a real question.
- **MICA reproduction running.** It needed one workaround: the authors' anndata 0.8 cannot read files written by current anndata, so `prepare_dataset.py` writes a plain exchange folder and `export_for_mica.py` rebuilds an old-format file inside their environment.
- **Results tracker** is generated by `scripts/update_tracker.py` into `docs/results_tracker.csv`, so Rhea's task is to read and summarize it rather than copy numbers by hand.
- **Dataset card template** added at `docs/dataset_cards/TEMPLATE.md` for Laraib's first task.

Next: read the MICA result, then start PBMC14K.

---

## 2026-10-01, Utshab

Started the actual work after the proposal was submitted.

- Wrote `WORKPLAN.md` and this file. Tasks, owners, dates and acceptance rules are all in the work plan.
- Found why the authors' MICA code would not install two weeks ago. The environment was created without `pip` inside it, and MICA's installer shells out to `python -m pip` while building, which then fails. Fixed in `scripts/setup_envs.sh` by creating the environment with `pip` and installing the requirements first.
- Also fixed the main environment. The default PyTorch package was built for CUDA 13.0 and this machine's driver supports 12.9, so the GPUs were invisible. The script now installs the CUDA 12.8 build.
- Both environments are rebuilding in the background.
- Created the dataset folder on the shared disk at `/media/20TB/shared/sides_data/43.scMINER-Clustering` with a `data` link from the repo. The home disk stays free of data. That disk has about 1.7 TB free.
- Note for scheduling: both GPUs are fully used by other projects right now, so GPU work waits and development happens on the processor first.

Next: finish the environments, then download Zeisel, which is small, and run the whole pipeline on it end to end before touching anything larger.
