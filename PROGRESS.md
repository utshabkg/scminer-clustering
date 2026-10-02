# Progress log

What has actually been done, newest first. The plan lives in [WORKPLAN.md](WORKPLAN.md).
Utshab adds an entry at the end of each working session. Laraib and Rhea add a line when they finish a task, or ask Utshab to add it.

## Where we stand

Updated: 2026-10-01

| Milestone | Due | Tasks done | Percent | State |
|-----------|-----|------------|---------|-------|
| M1 software, small data ready | Oct 7 | 0 of 10 | 0% | Started. Environments being fixed, data folder created |
| M2 reproduce the paper | Oct 14 | 0 of 6 | 0% | Not started |
| M3 cluster count without labels | Oct 21 | 0 of 4 | 0% | Not started |
| M4 what the similarity measures, atlas downloaded | Oct 28 | 0 of 5 | 0% | Not started |
| M5 fast neighbor search | Nov 11 | 0 of 5 | 0% | Not started |
| M6 full 2.3 million cell run | Nov 18 | 0 of 4 | 0% | Not started |
| M7 retrieval and drivers | Nov 25 | 0 of 5 | 0% | Not started |
| M8 paper and slides | Dec 2 | 0 of 4 | 0% | Not started |

Overall: 0 of 43 tasks finished.

Next three things: finish the two software environments (U-1, U-2), download Zeisel and PBMC14K (U-5, U-6), get GitHub usernames from Laraib and Rhea so they can see this repo (U-4).

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
