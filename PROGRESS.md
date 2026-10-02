# Progress log

What has actually been done, newest first. The plan lives in [WORKPLAN.md](WORKPLAN.md).
Utshab adds an entry at the end of each working session. Laraib and Rhea add a line when they finish a task, or ask Utshab to add it.

## Where we stand

Updated: 2026-10-02. Counts follow the Milestone column in `WORKPLAN.md`.

| Milestone | Due | Tasks done | Percent | State |
|-----------|-----|------------|---------|-------|
| M1 software, small data ready | Oct 7 | 6 of 8 | 75% | Software works, data ready. Laraib has joined GitHub. Waiting on Rhea's invite (U-4) and the Zeisel card (L-1) |
| M2 reproduce the paper | Oct 14 | 6 of 9 | 67% | Five of the paper's 13 benchmark sets run, four reproduce within 0.05. Left: decide on contacting the authors (U-30), PBMC14K card (L-2), results summary (R-1) |
| M3 cluster count without labels | Oct 21 | 0 of 6 | 0% | Not started. Now also holds the five small benchmark sets (U-27) and uncertainty estimates (U-28) |
| M4 what the similarity measures, atlas downloaded | Oct 28 | 0 of 5 | 0% | Not started |
| M5 fast neighbor search | Nov 11 | 0 of 4 | 0% | Not started |
| M6 full 2.3 million cell run | Nov 18 | 0 of 4 | 0% | Not started |
| M7 retrieval and drivers | Nov 25 | 0 of 4 | 0% | Not started |
| M8 paper and slides | Dec 2 | 0 of 3 | 0% | Not started |
| Every week: notes, plan and progress updates | weekly | | | Start at the Oct 7 meeting (L-7, R-3, R-8) |

Overall: 12 of 46 tasks finished. Three of the 46 repeat every week.

Next three things: add the paper's five small benchmark sets so we have ten datasets (U-27), estimate how much the numbers move between runs (U-28), and the label-free rule for the number of clusters (U-12).

---

## 2026-10-02, Utshab

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
