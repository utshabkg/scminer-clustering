# Work plan

Project: clustering single-cell data with scMINER, extending it, and writing it up as a paper.
Course: Clustering Algorithms, Missouri S&T, Fall 2026. Team: Utshab (lead), Laraib, Rhea.
Final deliverable: **early December 2026**. Weekly meeting: **Wednesday 2pm**.

## How to use this file

Each task has a status marker. Change the marker, commit, done.

- `[ ]` not started
- `[~]` in progress
- `[x]` finished
- `[!]` blocked, with a note saying why

We update this file in the Wednesday meeting. [PROGRESS.md](PROGRESS.md) is the other half of the picture: this file says what we plan to do, the progress file records what actually happened and how far along we are.

If a task is unclear, say so in the meeting instead of guessing. Every task below names the file it produces, so "finished" means that file exists and the person who asked for it has seen it.

## What we are trying to answer

We picked scMINER (Pan et al., Nature Communications 2025), which clusters cells using mutual information instead of ordinary distances. Our proposal raised three questions.

**RQ1. Does the method still win when nobody tells it how many clusters to find, and what is the similarity actually measuring?**
Answered by: a table of clustering scores with the number of clusters picked two ways, once using the true labels as the paper did and once without them, plus a measurement of how much of the similarity comes from which genes are detected rather than how much they are expressed.

**RQ2. Can it handle millions of cells on one machine, and what does the shortcut cost?**
Answered by: a run on the 2.3 million cell brain atlas with timing and memory, the fraction of true nearest neighbors our faster method finds, and how many small cell populations survive.

**RQ3. Does the regulatory information scMINER produces help explain the clusters?**
Answered by: cluster names proposed from retrieved literature using marker genes alone, drivers alone, and both, scored against the known labels.

## Who does what

Utshab writes all the code and runs every experiment, because only he has access to the workstation. Laraib and Rhea own work that needs no coding: documenting data, building comparison tables, checking biological claims, writing, and keeping the plan and results organized.

The first task for each of you (L-1, L-7, R-1) is set so we can get moving. Everything after that is a proposal: Rhea gets first pick at the Oct 7 meeting, and any task can be swapped between Laraib and Rhea. It should come to about 3 to 4 hours a week each.

## Milestones

There are 55 tasks: 39 for Utshab and 8 each for Laraib and Rhea. Three of them (L-7, R-3, R-8) repeat every week, and each of the other 52 belongs to the milestone in its row. U-25 and U-26 were added on Oct 2 to firm up the first finding, U-27 to U-30 came out of the Oct 2 review, and U-31 to U-38 come from the Assignment 2 feedback and our revised plan, and U-39 was added on Oct 22 because the PBMC14K result could not otherwise be attributed to size or to mode.

| ID | By | What is true when it is done |
|----|----|------------------------------|
| M1 | Oct 7 | Software works, the two small datasets are downloaded and cleaned, the repo structure is in place |
| M2 | Oct 14 | We reproduce the paper's numbers on those two datasets, the baselines and scoring code run, and the Assignment 2 report is submitted |
| M3 | Oct 21 | RQ1 part one: scores with and without using the labels to choose the number of clusters. The stability and significance rules run on the two anchor datasets |
| M4 | Oct 28 | RQ1 part two: what the similarity measures, plus the brain atlas downloaded and a timing pilot done. Every dataset run under the new rules, both gene settings, with repeats |
| M5 | Nov 4 | **Main experiments complete.** The new label-free rules have run on every dataset with repeats, and the final label-free comparison table is built |
| M6 | Nov 11 | Gated extras done or dropped: the GPU neighbor search, and the annotation comparison |
| M7 | Nov 18 | **Full manuscript draft**, with all figures and tables |
| M8 | Dec 2 | Revision from the feedback, slides, submitted |

## Tasks: Utshab

| ID | Task | Input | Output | Done when | Due | Milestone | Status |
|----|------|-------|--------|-----------|-----|-----------|--------|
| U-1 | Fix the main software environment so the GPUs are visible | `envs/scminer.yml` | `envs/scminer.lock.yml` | PyTorch reports the GPUs | Oct 7 | M1 | `[x]` |
| U-2 | Install the authors' MICA code | their GitHub repo | `envs/mica-legacy.lock.yml` | `mica -h` runs | Oct 7 | M1 | `[x]` |
| U-3 | Repo scaffolding: code package, run scripts, configs, results folders | plan | `src/`, `scripts/`, `configs/` | a test run writes a results folder | Oct 7 | M1 | `[x]` |
| U-4 | Add Laraib and Rhea as repo collaborators | their GitHub accounts | GitHub settings | both can open the repo | Oct 7 | M1 | `[x]` |
| U-5 | Download the Zeisel dataset and its labels | public source | `data/raw/zeisel/` | 3,005 cells load with labels attached | Oct 7 | M1 | `[x]` |
| U-6 | Assemble the PBMC14K dataset and labels | the authors' raw count file and label file | `data/raw/pbmc14k/` | about 13,605 cells with 7 labelled types | Oct 7 | M1 | `[x]` |
| U-7 | Preprocessing code: quality filters, normalization, log transform | raw data | `src/scminer_clust/preprocess.py`, processed files | both datasets processed, counts logged | Oct 7 | M1 | `[x]` |
| U-8 | Scoring code: ARI, AMI, NMI, silhouette, runtime, peak memory, repeated seeds | processed data | `src/scminer_clust/metrics.py` | a run writes `metrics.csv` with a seed column | Oct 14 | M2 | `[x]` |
| U-9 | Reproduce MICA on Zeisel and PBMC14K | processed data | `results/reproduction/metrics.csv` | our ARI is within about 0.05 of the published value, or the gap is explained | Oct 14 | M2 | `[x]` |
| U-10 | Baselines: k-means on principal components, Leiden, hierarchical | processed data | `results/baselines/metrics.csv` | all three run on both datasets | Oct 14 | M2 | `[x]` |
| U-11 | Export the first results table and figures for the team | results CSVs | `results/` committed | Rhea can build tables from them | Oct 14 | M2 | `[x]` |
| U-12 | Rule for choosing the number of clusters without labels | processed data | `src/scminer_clust/selection.py` | the rule picks a number on all methods, no labels used | Oct 21 | M3 | `[x]` |
| U-13 | Compare oracle setting against the label-free rule for every method | U-12 | `results/summary_label_free.md` | one table, both settings, both datasets | Oct 21 | M3 | `[x]` |
| U-14 | Measure what the similarity captures: detection versus expression, and the depth test | processed data | `results/rq1_diagnostics/` | shares reported, and cluster changes at 25, 50 and 75 percent of counts | Oct 28 | M4 | `[x]` |
| U-15 | Download the Allen brain atlas | public S3 bucket | `data/raw/allen/`, list in `docs/allen_files.csv` | file list and sizes recorded | Oct 28 | M4 | `[x]` |
| U-16 | Preprocess the atlas and run a 100,000 cell timing pilot | U-15 | `docs/allen_pilot.md` | measured time and memory, extrapolated to 2.3M | Oct 28 | M4 | `[x]` |
| U-17 | Our own mutual information distance, proven to match the authors' | their code | `tests/test_mi_matches_mica.py` | matches on 10,000 random cell pairs | Nov 4 | M5 | `[x]` |
| U-18 | GPU version of the distance | U-17 | `src/scminer_clust/mi.py` | same numbers as the CPU version, with a speedup recorded | Nov 11 | M6 | `[ ]` |
| U-19 | Fast neighbor graph: propose candidates cheaply, rescore exactly | U-17 (CPU first, U-18 for speed) | `results/rq2_graph_<dataset>/` | neighbor recall measured against brute force | Nov 4 | M5 | `[x]` CPU version done. Matches MICA's recall about 5 times faster on 100,000 cells |
| U-20 | Cluster a sketch of the atlas and assign the rest, full 2.3M run | U-19 | `results/rq2_atlas/` | the run finishes within a day and is scored | Nov 18 | M7 | `[ ]` |
| U-21 | Rare population analysis | U-20 | `results/rq2_atlas/rare_types.csv` | per-type recall, and which types were lost | Nov 18 | M7 | `[ ]` |
| U-22 | Driver activity in Python, including hidden drivers | PBMC14K clusters | `results/rq3_drivers/` | driver lists per cluster | Nov 11 | M6 | `[ ]` |
| U-23 | Retrieval and local language model annotation, three query types | U-22 | `results/rq3_rag/` | names and cited sources for every cluster | Nov 11 | M6 | `[ ]` |
| U-24 | Final figures, method and results writing, repo cleanup | everything | `latex/sections/`, `results/` | a fresh clone reproduces the main table | Dec 2 | M8 | `[ ]` |
| U-25 | Run MICA on the same 2,000 genes the baselines got, all three selection methods | processed data | `results/mica_*_hvg2000_*/` | a two by two table of method against gene set on every dataset | Oct 14 | M2 | `[x]` |
| U-26 | Add the rest of the paper's eligible benchmark sets: Usoskin, Kolodziejczyk, Klein | Bioconductor mirror and GEO | `data/raw/`, `results/reproduction_*/`, `results/baselines_*/` | MICA and baselines scored on all five datasets | Oct 14 | M2 | `[x]` |
| U-27 | Add the paper's five small benchmark sets (Yan, Goolam, Buettner, Pollen, Chung) from the authors' Zenodo inputs, MICA and Leiden at all four gene settings | Zenodo 15040179 | `results/` | ten datasets in the gene-set table | Oct 21 | M3 | `[x]` |
| U-28 | Estimate run-to-run spread: paired bootstrap of the ARI and AMI differences, five repeat graph-mode MICA runs on PBMC14K and Zeisel | saved partitions | `results/uncertainty/` | intervals next to every gap in the table | Oct 21 | M3 | `[x]` |
| U-29 | Kolodziejczyk side test: does removing the plate effect, or the authors' gene filter, recover ARI 1.0? Report any use of the labels as leakage | released input | `results/kolod_side_test/` | a short note on what recovers the published score | Oct 28 | M4 | `[x]` |
| U-30 | Decide whether to ask the scMINER authors for their unreleased Kolodziejczyk input before we describe the non-reproduction anywhere outside the team | R5 evidence in PROGRESS.md | email or a recorded decision | decision recorded in the decision log | Oct 14 | M2 | `[ ]` |
| U-31 | Assignment 2 report: revised question, closest work, preliminary results, plan and schedule | the feedback and our results | `latex/` and the submitted PDF | submitted by Oct 14 | Oct 14 | M2 | `[x]` |
| U-32 | Stability rule: pick the setting whose clusters best survive subsampling, for MICA and Leiden alike. Varies $k$ in MICA's small-data mode | processed data | `src/scminer_clust/selection.py`, `results/rule_stability/` | the rule picks a setting for both methods on all ten sets, no labels used | Oct 21 | M3 | `[ ]` |
| U-33 | Significance rule: merge clusters a test cannot separate, in the style of sc-SHC | U-32 | `results/rule_significance/` | the rule runs on both methods and reports the merged cluster count | Oct 21 | M3 | `[ ]` |
| U-34 | Transfer rule: give each method the resolution that scored best on the other datasets. Also score MICA's default on the Allen subset, outside the paper's benchmark | U-32 | `results/rule_transfer/` | one held-out score per dataset per method | Nov 4 | M5 | `[ ]` |
| U-35 | Rerun the silhouette rule in one shared embedding, over the same range of cluster counts for both methods | U-32 | `results/rule_silhouette_shared/` | both methods scored in the same space over the same range of k | Oct 28 | M4 | `[ ]` |
| U-36 | Label-free Leiden on the Allen 100,000 cell subset, so it becomes a third anchor | U-16 | `results/labelfree_allen_100k/` | Leiden scored under every rule on the Allen subset | Oct 28 | M4 | `[ ]` |
| U-37 | Final label-free comparison table across all eleven sets, with the win criterion from the report applied | U-32 to U-36 | `results/summary_label_free_rules.md` | one table, every rule, every dataset, with repeats | Nov 4 | M5 | `[ ]` |
| U-38 | Optional checks if time allows: Leiden on the correlation neighbor graph, and a hub correction on the mutual information graph | U-19 | `results/rq2_swap/`, `results/rq2_hubness/` | both run or formally dropped at the Nov 4 meeting | Nov 11 | M6 | `[ ]` |
| U-39 | Separate MICA's mode from the dataset size: one fixed pool of 2,500 PBMC14K cells run through both of MICA's modes under the same rule | U-32 | `scripts/run_rule_mode_size.py`, `results/rule_mode_size_pbmc14k_n2500_seed0/` | both arms report a tie span on the same cells, so mode and size are no longer confounded | Nov 4 | M5 | `[~]` |

## Tasks: Laraib

| ID | Task | What to do | Output | Done when | Due | Milestone | Status |
|----|------|-----------|--------|-----------|-----|-----------|--------|
| L-1 | Dataset card: Zeisel | Fill the template: where the data comes from, how many cells and genes, what the labels mean, who made them, the license, and what preprocessing we applied | `docs/dataset_cards/zeisel.md` | every field filled, with links | Oct 7 | M1 | `[x]` |
| L-2 | Dataset card: PBMC14K | Same template. Note that labels come from sorting cells with antibodies, not from clustering | `docs/dataset_cards/pbmc14k.md` | every field filled | Oct 14 | M2 | `[ ]` |
| L-3 | Dataset card: Allen brain atlas | Same template, plus the label levels (34 classes, 338 subclasses) and the license restriction | `docs/dataset_cards/allen_wmb.md` | every field filled | Oct 21 | M3 | `[ ]` |
| L-4 | Related work comparison table | Using the saved reviews in `docs/research/`, build a table of about 15 papers with columns: paper, year, venue, what it does, biggest data it handled, does it use mutual information, does it choose cluster count without labels, does it use a language model, how ours differs | `docs/related_work_table.md` | 15 rows, each with a link | Oct 28 | M4 | `[ ]` |
| L-5 | Rewrite background and related work in Overleaf | Use L-4 so the paper's story matches the table | `latex/sections/01_background.tex`, `05_related_work_contribution.tex` | sections updated and compiling | Nov 11 | M5 | `[ ]` |
| L-6 | Marker check | For each cluster name our system proposes, look the markers up in CellMarker 2.0 and PanglaoDB on their websites, and record whether they agree, with the link you used | `results/rq3_rag/marker_check.csv` | every proposed name checked | Nov 18 | M6 | `[ ]` |
| L-7 | Meeting minutes | Short notes after each Wednesday meeting: decisions, who owns what next | `docs/minutes/YYYY-MM-DD.md` | one file per meeting | weekly | every week | `[ ]` |
| L-8 | Proofread the paper | Check the writing rules we agreed on: no long dashes, no semicolons, plain words, consistent terms | the LaTeX sections | full read-through with fixes | Dec 2 | M8 | `[ ]` |

## Tasks: Rhea

| ID | Task | What to do | Output | Done when | Due | Milestone | Status |
|----|------|-----------|--------|-----------|-----|-----------|--------|
| R-1 | Results summary | Utshab refreshes `docs/results_for_pivot.csv` after every run. It holds only the main comparison, one row per run (the full `docs/results_tracker.csv` has every kind of row and is not meant for this). Download it from GitHub, open it in Excel or Google Sheets, and build a one-page pivot: rows are dataset and gene_setting, columns are method, values are the average of ari. Keep the two inputs (authors' data and our preprocessing) apart, for example with a filter on the input column. Mark anything that looks odd and bring it to the meeting | `docs/author_data_results_summary.csv` and `docs/our_preprocessing_results_summary.csv` | every run in `docs/results_for_pivot.csv` appears in the summary, with odd results marked | Oct 21 | M3 | `[x]` odd results marked in a new anomaly column of both CSVs |
| R-2 | Turn results into paper tables | Take the result CSVs and format them as LaTeX tables with captions | `latex/tables/r2_*.tex`, made by `scripts/make_paper_tables.py` | tables compile in the draft | Oct 21 | M3 | `[x]` both tables compile with the report preamble. Rerun the script when the summary CSVs change |
| R-3 | Keep this plan current | After each meeting, update statuses, the milestone table and the risk list | `WORKPLAN.md` | statuses match reality each week | weekly | every week | `[ ]` |
| R-4 | Rare population sheets | From the per-cluster counts, work out how many reference cell types survive clustering and which ones vanish, then summarize | `results/rq2_atlas/rare_summary.csv` | summary table with a short written finding | Nov 18 | M6 | `[ ]` |
| R-5 | Annotation audit | For about 50 sampled annotations, use the rubric to mark whether the quoted text really supports the proposed name, yes or no, with a reason | `results/rq3_rag/audit.csv` | 50 rows scored | Nov 25 | M7 | `[ ]` |
| R-6 | Slide outline | Draft 10 to 12 slides: problem, limitation, our questions, method, results, what we learned | `docs/slides_outline.md` | outline agreed in the meeting | Nov 25 | M7 | `[ ]` |
| R-7 | Number consistency check | Every number in the text matches the tables, and every table is referenced somewhere | the paper | one pass with corrections listed | Dec 2 | M8 | `[ ]` |
| R-8 | Update the progress file | After each meeting, refresh the milestone status table at the top of the progress log | `PROGRESS.md` | percentages match the task statuses | weekly | every week | `[ ]` |

## Weekly schedule

| Week | Meeting | Milestone | Utshab | Laraib | Rhea |
|------|---------|-----------|--------|--------|------|
| Oct 1 to 7 | Oct 7 | M1 | U-1 to U-7 | L-1 | R-1 (start) |
| Oct 8 to 14 | Oct 14 | M2 | U-8 to U-11, U-25, U-26 | L-2 | R-1 |
| Oct 15 to 21 | Oct 21 | M3 | U-32, U-33 | L-3 | R-1, R-2 |
| Oct 22 to 28 | Oct 28 | M4 | U-35, U-36 | L-4 | R-2 |
| Oct 29 to Nov 4 | Nov 4 | M5 | U-34, U-37 | L-4 | weekly tasks only |
| Nov 5 to 11 | Nov 11 | M6 | U-18, U-22, U-23, U-38 | L-5 | weekly tasks only |
| Nov 12 to 18 | Nov 18 | M7 | U-20, U-21, U-24 starts | L-6 | R-4 |
| Nov 19 to 25 | Nov 25 | M8 | U-24 | L-6 | R-5, R-6 |
| Nov 26 to Dec 2 | Dec 2 | M8 | U-24 | L-8 | R-7 |

Every week, starting Oct 7: Laraib takes the meeting notes (L-7), and after the meeting Rhea updates the statuses in this file (R-3) and the table at the top of `PROGRESS.md` (R-8). Utshab will show how to edit a file on GitHub at the Oct 7 meeting.

Standing meeting agenda: read the progress log, update statuses here, clear anything blocked, agree the next week's three headline tasks.

## How we hand work over

- Utshab pushes results as `results/<experiment>/metrics.csv` plus any figures, and adds a one-line note in `PROGRESS.md` saying what the run was.
- Laraib and Rhea read those files on GitHub. Nobody else needs to run code or log into the workstation.
- For spreadsheet work, download the CSV, work in Excel or Google Sheets, then send the file back and Utshab commits it, or edit directly on GitHub if you are comfortable with that.
- If a file you need is missing or confusing, say so in the meeting or message Utshab. Do not rebuild numbers by hand.
- Writing happens in Overleaf. The LaTeX files in `latex/sections/` are the same ones, one file per section, so two people can edit different sections without clashing.

## Where things live

```
scminer-clustering/
  WORKPLAN.md          this file
  PROGRESS.md          what has actually been done
  latex/               the paper, one file per section under latex/sections/
  src/scminer_clust/   our code
  scripts/             one runner per experiment
  configs/             settings for each run, so results can be repeated
  results/<id>/        metrics.csv, params.json, figures
  results/summary_*.md the main tables, rebuilt after every run
  docs/results_for_pivot.csv  the main comparison, one row per run (for R-1)
  docs/results_tracker.csv    every row from every experiment
  docs/allen_pilot.md  atlas pilot: time, memory, what 2.3 million cells would cost
  docs/dataset_cards/  one page per dataset
  docs/research/       the literature reviews we did before the proposal
  docs/minutes/        meeting notes
  data -> /media/20TB/shared/sides_data/43.scMINER-Clustering  (not in git)
```

Datasets never go into git. They live on the shared 20 TB disk, and the repo only holds code, settings, small result tables and figures.

## Risks and what we do about them

| Risk | Sign it is happening | What we do |
|------|----------------------|------------|
| Both GPUs busy with other projects | No free memory when a run is due | Develop and test on the processor, run GPU jobs overnight, and ask before starting a long one |
| The authors' MICA code will not install | Install errors keep coming back | Call their distance functions straight from the source checkout, which is all we need to prove our version matches |
| PBMC14K matrix is hard to rebuild | Cell barcodes do not line up with the label file | Use the authors' preprocessed copy or a published subset of the same data |
| The brain atlas is too big or too slow to download | Download runs for days, or the disk fills | Use one or two brain regions as a smaller version, or switch to the Cao 2020 atlas, which is one file |
| Network inference is too slow for every cluster | A single cluster takes hours | Restrict it to one branch of the label hierarchy, or group cells first, which is what the authors do |
| The language model will not fit while the GPUs are busy | Out of memory on load | Run it at night, or use a smaller model for the pilot |
| Someone publishes our retrieval idea first | A new preprint appears | Our main contribution stays the scaling and the evaluation, which are separate findings |
| Shared disk fills up | Less than 200 GB free | Delete intermediate files after each conversion and keep only what we need |
| A teammate gets busy | A task slips two weeks | Tasks are weekly sized and independent, so they can be swapped or dropped without blocking the runs |

## Decisions so far

| Date | Decision | Why |
|------|----------|-----|
| Sep 15 | Chose scMINER as the paper to extend | Mutual information clustering plus regulatory drivers suits both the course and our interest in retrieval |
| Sep 16 | Everything in Python, no R | Nobody on the team writes R, and the clustering and network tools are already Python |
| Sep 16 | Datasets: PBMC14K, Allen brain atlas, Zeisel | They match the three size regimes the course requires and all have labels |
| Sep 23 | Proposal reorganized around three research questions | The instructor asked for a problem to limitation to question to method line of reasoning |
| Oct 1 | Work order: reproduce, then RQ1, then RQ2, then RQ3 | The cheap processor-only results come first, so we have findings even if the big runs slip |
| Oct 1 | Data lives on the shared 20 TB disk | Keeps the home disk free |
| Oct 1 | Ask before long GPU runs | The GPUs are shared with other projects |
| Oct 1 | Filtering at 4 deviations on a log scale | Matches the paper's cell counts: 13,646 against 13,605 on PBMC14K |
| Oct 1 | Baselines and MICA both picked with the paper's rule | So the reproduction compares like with like. The label-free rule is added separately for RQ1 |
| Oct 2 | Test MICA on 2,000 genes too, and add three more benchmark sets | Leiden on all genes beat MICA on two datasets, and we need to know whether the genes or the method explain it before claiming anything |
| Oct 2 | Keep every labelled cell for Usoskin, Kolodziejczyk and Klein | The paper's cell counts for these sets equal every labelled cell, so the authors dropped none |
| Oct 2 | Leiden resolution swept in steps of 0.01, as in the authors' code | The paper's text says 0.1, but that grid often cannot produce the true number of clusters |
| Oct 2 | Report MICA as three seeded k-means draws on its own embedding | MICA's default is one unseeded draw, while the baselines report three seeds |
| Oct 2 | Treat Kolodziejczyk as not reproduced from the released files | MICA gets 0.727 on the authors' released input, as does their own notebook. Their 1.0 came from an input they did not release |
| Oct 3 | Read each of the authors' inputs at its own log base | Pollen and Buettner are log10, Chung and Usoskin natural log, the rest log2. Reading all of them as log2 picked the wrong seurat and seurat_v3 genes |
| Oct 3 | Rerun MICA's final k-means with 100 restarts | The published scores match the best k-means fit, and one run often misses it by a lot (Usoskin 0.66 against 0.89) |
| Oct 3 | Call the label-free comparison inconclusive for now | Which method wins flips with the k range and with the embedding used for the silhouette |
| Oct 3 | Drop the landmark proof of exact neighbor lists | Mutual information distances bunch together, so the bounds prove under 10 percent of lists exact |
| Oct 14 | Narrow the project to the label-free question, after the instructor's feedback | Reproduction, scaling, label-free choice and the annotation study were four projects at once. Scaling and annotation now follow the main question |
| Oct 14 | Judge the fast neighbor graph by MICA's own later steps, not by Leiden on the raw graph | Leiden on the raw graph scores 0.64 on PBMC14K against MICA's 0.85, so it is a poor stand-in |
| Oct 14 | Fix a win criterion before the runs | Without one, picking the rule after seeing the scores would repeat the selection bias we are testing for |
| Oct 14 | Move the full 2.3 million cell run to a stretch goal | The pilot projects 7 to 8 days and about 3.2 TB of memory, against our 754 GB |
| Oct 3 | Fast neighbor search proposes candidates by correlation over all genes, then scores them exactly | 50 principal components miss most of MICA's neighbors on PBMC14K (33% found), correlation finds 97% with 800 candidates |
| Oct 3 | Judge a fast graph by MICA's own later steps, not only by neighbor recall | Graphs that miss many exact neighbors can cluster as well or better, because MICA's graph has strong hub cells |
