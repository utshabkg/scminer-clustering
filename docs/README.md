# Docs

## Start with these

| File | What it is |
|---|---|
| `primer.md` | The data, the methods and the question in plain words, assuming no background |
| `protocol_matched_v1.md` | Protocol v1, written and tagged (`matched-protocol-v1`) before the confirmatory runs. The paper cites it |
| `dataset_facts.md`, `dataset_facts.csv` | Cells, genes and cell types per dataset after filtering. Written by `scripts/setup/dataset_facts.py` |
| `dataset_cards/` | Cards for Zeisel, PBMC14K and the Allen 100,000-cell sample: source, labels, licence, preprocessing. `TEMPLATE.md` is the blank form. The other datasets have their numbers in `dataset_facts.md` |

## Team trackers and checks

| File | Kept by | What it is |
|---|---|---|
| `provenance.csv` | Rhea (R-9) | For each result: its design, when the design was fixed, and whether git shows that order. The source of the paper's last appendix table |
| `seed_coverage.csv` | Rhea (R-10) | Seeds done and still needed for each comparison, with a reason for every gap |
| `r7_number_check.md`, `r7_number_check_a4.md` | Rhea (R-7) | Number checks of the Assignment 3 report and the Assignment 4 draft |
| `reading_list_a4.md` | Utshab, for Laraib | The checked references behind the Assignment 4 framing |
| `related_work_table.md` | Laraib (L-4) | A comparison table of the closest related papers |
| `minutes/` | Laraib (L-7) | Where meeting notes go, one file per meeting. None added yet |

## Inputs the scripts read

| File | What it is |
|---|---|
| `published_scores.csv` | The ARI and AMI the scMINER paper reports per dataset and method |
| `results_tracker.csv` | Every row of every `results/*/metrics.csv` as of Oct 3. `scripts/a3_report.py` reads it, so it is not rebuilt (see `scripts/README.md`) |
| `results_for_pivot.csv`, `author_data_results_summary.csv`, `our_preprocessing_results_summary.csv` | Rhea's R-1 pivot of the Assignment 2 runs and its two summaries |
| `allen_files.csv`, `allen_pilot.md` | The Allen atlas files, and the 100,000-cell pilot that sizes the 2.3 million cell run |

## History

`research/` holds the literature sweeps of September, before the project's question was fixed, and the
Oct 3 check of our first claims. It is kept as a record and is not updated.
