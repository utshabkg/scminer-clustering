# The paper

`main.tex` builds the full research paper, which grew out of Assignment 4. It reads
`sections/` (one file per section), `tables/`, `figures/`, `references.bib` and the two number files,
`numbers_a4.tex` and `numbers_large.tex`. The submitted PDF is compiled on Overleaf with pdfLaTeX,
so keep the source in plain ASCII and write accents as LaTeX commands, for example `Du{\`o}`.

## Where the tables and numbers come from

| In the paper | File | Written by |
|---|---|---|
| Table I (datasets) | `tables/datasets.tex` | by hand. The Allen atlas's sizes are `\V{}` values from `numbers_large.tex`, and the other datasets' sizes come from `docs/dataset_facts.csv` |
| Tables II to VI (criteria, confirmatory runs, why matching matters, candidate designs, PBMC14K) | `tables/a4_criteria.tex`, `a4_confirm.tex`, `a4_decomp.tex`, `a4_designs.tex`, `a4_pbmc.tex` | `scripts/a4_protocol_v1.py --tables` |
| Figure 1 (stability curves) and Figure 2 (per-seed means) | `figures/a4_stability_curve.pdf`, `a4_seed_means.pdf` | the same |
| Table VII (the Allen atlas) | `tables/large_scale.tex` | `scripts/large_scale/make_tables.py` |
| Table VIII (effects) and Table IX (checks) | `tables/effects.tex`, `checks.tex` | by hand. Every result in them is a `\V{}` value from `numbers_a4.tex` |
| Tables X, XII, XIII, XIV, XVI and XVII (appendix) | `tables/a3_datasets.tex`, `a3_rules_all.tex`, `a3_common_range.tex`, `a3_mode_size.tex`, `a3_ami.tex`, `a3_graph_cost.tex` | `scripts/a3_report.py` |
| Table XI (reproduction) | `tables/a2_reproduction.tex` | `scripts/assignment2/a2_report.py` |
| Table XV (transfer) | `tables/a3_transfer.tex` | `scripts/a3_matched.py` |
| Table XVIII (when each design was fixed) | typed in `sections/09_appendix.tex`, from `docs/provenance.csv` | by hand |

Tables are numbered in the order they appear in the source, so moving one renumbers the rest. Update
this list when that happens.

- `numbers_a4.tex` holds most protocol v1 numbers the text quotes. `scripts/a4_protocol_v1.py --tables`
  writes it. `numbers_large.tex` holds the Allen atlas numbers, and `scripts/large_scale/make_tables.py`
  writes it. The text prints a value with `\V{key}`, and an unknown key shows up in bold as `??key`.
  `main.tex` reads `numbers_large.tex` only if the file is there. Without it the paper still compiles,
  but the Allen numbers in the text and in Table I show as bold `??large.<key>`. Table VII keeps its
  numbers, because the script writes them straight into the table.
- Still typed by hand: a few protocol v1 numbers in the text, mostly values for one dataset or one seed
  quoted as examples (the Klein picks, the PBMC14K seed counts and the PBMC14K cost, for instance), and
  the Assignment 2 and 3 figures quoted in the text and the appendix. `docs/r7_number_check_a4.md` lists
  the typed numbers as of Nov 13, each with its source. Some have since moved into `\V{}`. The
  Assignment 3 column of Table II is no longer typed. `scripts/a4_protocol_v1.py` writes it from
  `results/a3_matched_numbers.json`, and the text quotes the same values through the `A3.*` keys.
- A few files in `tables/` and `figures/` are written by scripts but no longer used by the paper:
  `a2_selection.tex`, `a3_matched.tex`, `a3_summary.tex`, `r2_*.tex` (from
  `scripts/tools/make_paper_tables.py`), `a2_neighbor_recall.pdf`, `a3_stability_curve.pdf` (the old
  Figure 1) and the three `.png` figures. The notebooks show `a2_neighbor_recall.png` and
  `a3_stability_curve.png`. `a4_stability_curve.png` is written next to the PDF of Figure 1.

## Overleaf

Overleaf does not sync with this repository, so files move by hand. Download the team's Overleaf files
into `latex/updated_tex_<date>/` (ignored by git). The current paper is the folder
`3. Project Revision and Further Results`. Compare it with this folder, then run
`python scripts/tools/overleaf_package.py --overleaf "latex/updated_tex_<date>/3. Project Revision and Further Results"`.
It writes the changed files, a full zip and an upload list to `scratchpad/overleaf_upload/`. Upload them
into that same Overleaf folder.

## Earlier versions

`archive/proposal/` and `archive/assignment2/` are the proposal and Assignment 2 as submitted, and each
compiles on its own. Assignment 3 and Assignment 4 as submitted are the git tags `a3-submitted` and `a4-submitted`.
