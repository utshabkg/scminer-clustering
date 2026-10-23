# Research archive

Working notes from the planning and verification stages, kept so the choices in the report can
be traced back to what we read. Nothing here is a result. The results are in `results/`, the
report is in `latex/`, and what was actually done and when is in `PROGRESS.md`.

## Literature sweeps (September 2026)

Written while choosing the question. They are dated, and some of what they describe as open has
since been settled by our own runs.

| File | Covers |
|---|---|
| `workflow_paused_2026-09-16/research_scminer-deep-read.md` | scMINER and MICA in detail: the benchmark protocol, the two clustering modes, the peer review exchange, and where the shipped defaults come from |
| `workflow_paused_2026-09-16/research_resolution-validity.md` | choosing the number of clusters, judging whether a split is real, and the known pitfalls in single-cell evaluation |
| `workflow_paused_2026-09-16/research_scalable-mi-clustering.md` | mutual information similarity and clustering at million-cell scale |
| `workflow_paused_2026-09-16/research_datasets.md` | candidate public datasets for the three regimes we considered |
| `workflow_paused_2026-09-16/research_llm-rag-annotation.md` | retrieval and language model approaches to annotating clusters, for the component we have not yet run |

`ideate_*.json` in the same folder holds the three project directions we weighed before settling
on the one in the report. `journal.jsonl` and `workflow_script*.js` are the machine artifacts
that produced these notes, kept only for provenance.

## Verification (October 2026)

`verification/claims_oct3.md` lists the claims we wanted checked against the saved runs, and
`verification/oct3_findings.txt` records what the check found. It is worth reading: it is where
we caught two real bugs, a log-base error that silently changed which genes the `seurat`
selector picked, and a summary script that mixed MICA's small-data and graph-mode runs into one
Zeisel average. Both are fixed, and the second is written up in the report's anomalies appendix.
