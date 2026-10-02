# Progress log

What has actually been done, newest first. The plan lives in [WORKPLAN.md](WORKPLAN.md).
Utshab adds an entry at the end of each working session. Laraib and Rhea add a line when they finish a task, or ask Utshab to add it.

## Where we stand

Updated: 2026-10-01. Counts follow the Milestone column in `WORKPLAN.md`.

| Milestone | Due | Tasks done | Percent | State |
|-----------|-----|------------|---------|-------|
| M1 software, small data ready | Oct 7 | 6 of 8 | 75% | Software works, Zeisel and PBMC14K ready. Waiting on the GitHub invites (U-4) and the Zeisel card (L-1) |
| M2 reproduce the paper | Oct 14 | 4 of 6 | 67% | Both datasets reproduced within tolerance, baselines done. Left: PBMC14K card (L-2), results summary (R-1) |
| M3 cluster count without labels | Oct 21 | 0 of 4 | 0% | Not started |
| M4 what the similarity measures, atlas downloaded | Oct 28 | 0 of 4 | 0% | Not started |
| M5 fast neighbor search | Nov 11 | 0 of 4 | 0% | Not started |
| M6 full 2.3 million cell run | Nov 18 | 0 of 4 | 0% | Not started |
| M7 retrieval and drivers | Nov 25 | 0 of 4 | 0% | Not started |
| M8 paper and slides | Dec 2 | 0 of 3 | 0% | Not started |
| Every week: notes, plan and progress updates | weekly | | | Start at the Oct 7 meeting (L-7, R-3, R-8) |

Overall: 10 of 40 tasks finished. Three of the 40 repeat every week.

Next three things: run MICA on the same 2,000 genes the baselines got, so we know whether the gap comes from the genes or from the method, add Usoskin, Kolodziejczyk and Klein from the paper's benchmark, and the label-free rule for the number of clusters (U-12).

-----------|-----|------------|---------|-------|
| M1 software, small data ready | Oct 7 | 6 of 10 | 60% | Both environments work, Zeisel and PBMC14K downloaded and preprocessed, scaffolding in place |
| M2 reproduce the paper | Oct 14 | 4 of 6 | 67% | Both datasets reproduced within tolerance, baselines done, first results shared |
| M3 cluster count without labels | Oct 21 | 0 of 4 | 0% | Not started |
| M4 what the similarity measures, atlas downloaded | Oct 28 | 0 of 5 | 0% | Not started |
| M5 fast neighbor search | Nov 11 | 0 of 5 | 0% | Not started |
| M6 full 2.3 million cell run | Nov 18 | 0 of 4 | 0% | Not started |
| M7 retrieval and drivers | Nov 25 | 0 of 5 | 0% | Not started |
| M8 paper and slides | Dec 2 | 0 of 4 | 0% | Not started |

Overall: 10 of 43 tasks finished.

Next three things: the label-free rule for choosing the number of clusters (U-12), extend the gene-set comparison to Usoskin, Kolodziejczyk and Klein, get GitHub usernames from Laraib and Rhea (U-4).

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
