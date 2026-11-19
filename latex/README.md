# The paper

`main.tex` builds the current paper, Assignment 4, from `sections/` (one file per section),
`tables/`, `figures/` and `references.bib`. The submitted PDF is compiled on Overleaf with pdfLaTeX,
so keep the source in plain ASCII and write accents as LaTeX commands, for example `Du{\`o}`.

## Where the tables and numbers come from

| In the paper | File | Written by |
|---|---|---|
| Tables I to V (criteria, confirmatory runs, why matching matters, candidate designs, PBMC14K) | `tables/a4_criteria.tex`, `a4_confirm.tex`, `a4_decomp.tex`, `a4_designs.tex`, `a4_pbmc.tex` | `scripts/a4_protocol_v1.py --tables` |
| Figure 2 (per-seed means) | `figures/a4_seed_means.pdf` | the same |
| Figure 1 (stability curves), Tables VI, VIII, IX, X, XII and XIII | `figures/a3_stability_curve.pdf`, `tables/a3_datasets.tex`, `a3_rules_all.tex`, `a3_common_range.tex`, `a3_mode_size.tex`, `a3_ami.tex`, `a3_graph_cost.tex` | `scripts/a3_report.py` |
| Table XI (transfer) | `tables/a3_transfer.tex` | `scripts/a3_matched.py` |
| Table VII (reproduction) | `tables/a2_reproduction.tex` | `scripts/assignment2/a2_report.py` |
| Table XIV (when each design was fixed) | typed in `sections/09_appendix.tex`, from `docs/provenance.csv` | by hand |

- `numbers_a4.tex` holds most protocol v1 numbers the text quotes. The text prints one with `\V{key}`,
  and an unknown key shows up in bold as `??key`. The same script writes it.
- Still typed by hand: about 20 protocol v1 numbers in the text (each listed with its source in
  `docs/r7_number_check_a4.md`), the Assignment 3 column of Table I, and the Assignment 2 and 3 figures
  quoted in the text.
- A few files in `tables/` and `figures/` are written by the scripts but no longer used by the
  paper: `a2_selection.tex`, `a3_matched.tex`, `a3_summary.tex`, `r2_*.tex`, `a2_neighbor_recall.pdf`
  and the two `.png` figures, which the notebooks show.

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
