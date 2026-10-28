"""Build and run notebooks/assignment3_results.ipynb, which rebuilds the Assignment 3 displays.

Usage: python scripts/make_a3_notebook.py
The notebook reads only saved result files, so it runs in under a minute. The
experiments themselves are run by the scripts it names.
"""
import sys

import nbformat
from nbclient import NotebookClient

from scminer_clust.paths import REPO

CELLS = [
    ("md", """# Assignment 3: label-free selection results

This notebook rebuilds every table and figure in the Assignment 3 manuscript from saved result files. It does not rerun any experiment. Each number traces back to a folder under `results/` holding a `metrics.csv` and a `params.json` with that run's settings, seed and candidate grid.

| Display | Produced by |
|---|---|
| Table I, ARI and chosen cluster count under the stability rule | `scripts/run_rule_stability.py`, `run_rule_stability_graph.py` |
| Table II, graph quality and cost | `scripts/measure_graph_cost.py`, `measure_mica_hnsw.py`, the `rq2_*` audits and downstream runs |
| Figure 1, stability against cluster count | the stability runs |
| Table III (appendix), the datasets | `scripts/dataset_facts.py` |
| Table IV (appendix), all four rules | the stability runs plus `run_rule_silhouette.py` and `run_rule_default.py` |
| Table V (appendix), the stability rule on a common range | `scripts/a3_report.py`, from the saved stability curves |
| Table VI (appendix), mode against size on one fixed pool of cells | `scripts/run_rule_mode_size.py` |
| Table VII (appendix), adjusted mutual information | the same runs as Table IV |
| The mechanism numbers | `scripts/a3_report.py` |

The environment is `envs/scminer.yml`. Data live outside the repository, as described in `README.md`."""),
    ("code", """import json, sys
import pandas as pd
from IPython.display import Image, display
from scminer_clust.paths import REPO, RESULTS
sys.path.insert(0, str(REPO / "scripts"))
import a3_report
a3_report.main()   # rewrites latex/tables/a3_*.tex, the figure, and results/a3_report_numbers.json
nums = json.loads((RESULTS / "a3_report_numbers.json").read_text())"""),
    ("md", "## Table IV (appendix). ARI of the partition each rule chose\n\nNo reference labels are used to make the choice. Median over the repeats, with the observed range. The means quoted in the paper cover the nine small-data datasets, so PBMC14K is listed but excluded from them."),
    ("code", """rows = []
for rule in ["default", "silhouette_shared", "stability_finest", "stability_most"]:
    if rule not in nums or "per_dataset" not in nums[rule]:
        continue
    per = nums[rule]["per_dataset"]
    for ds in per["mica"]:
        rows.append({"rule": rule, "dataset": ds, "MICA": per["mica"][ds],
                     "Leiden": per["leiden_pca"][ds],
                     "MICA minus Leiden": round(per["mica"][ds] - per["leiden_pca"][ds], 3)})
table2 = pd.DataFrame(rows).pivot(index="dataset", columns="rule",
                                  values=["MICA", "Leiden", "MICA minus Leiden"])
table2.round(3)"""),
    ("md", "## Summary across rules\n\nThe spread across rules is the result: it is larger for MICA than for the baseline, and larger than the gap between the methods under any one rule."),
    ("code", """summary = pd.DataFrame([
    {"rule": r, "MICA mean": nums[r]["mica_mean"], "Leiden mean": nums[r]["leiden_mean"],
     "MICA wins": f'{nums[r]["mica_wins"]} of {nums[r]["n_datasets"]}'}
    for r in ["default", "silhouette_shared", "stability_finest", "stability_most"]
    if r in nums and "mica_mean" in nums[r]])
spread = summary[["MICA mean", "Leiden mean"]].max() - summary[["MICA mean", "Leiden mean"]].min()
display(summary.round(3))
print("spread across rules:", spread.round(3).to_dict())"""),
    ("md", "## Table V (appendix). The stability rule on a common range\n\nThe 2 to 20 filter caps both grids without matching them. Here the rule is re-applied to MICA's saved scores restricted to the cluster counts the baseline reaches in the same run. Nothing is rerun and no labels are used. This is the analysis that removed the ranking reversal and the apparent mode effect."),
    ("code", """cr = nums["common_range"]
display(pd.DataFrame(cr["means"]).T.round(3))
print("baseline reaches:", cr["baseline_reach"])
print("pooled tied span:", cr["pooled_tied_span"])
print("mode test on a common range:", json.dumps(nums["mode_size_common_range"], indent=1))"""),
    ("md", "## Why the rules disagree\n\nAmong settings the stability rule scores within 0.01 of its best, how wide a span of cluster counts remains?"),
    ("code", """mech = pd.read_csv(RESULTS / "a3_mechanism.csv")
display(mech[["ds", "method", "share_tied", "k_tied_min", "k_tied_max", "chosen_k", "true_k"]])
print({k: v for k, v in nums["mechanism"].items() if k != "table"})"""),
    ("md", "## Figure 1. Stability and accuracy against the number of clusters"),
    ("code", """display(Image(filename=str(REPO / "latex" / "figures" / "a3_stability_curve.png"), width=900))"""),
    ("md", "## Table II. Graph quality and cost\n\nMICA's HNSW search against a correlation shortlist rescored exactly. The downstream ARIs let the labels pick the resolution, since this compares graphs rather than selection rules."),
    ("code", """print(json.dumps(nums.get("graph_cost", {}), indent=1))"""),
    ("md", "## MICA's shipped default against this benchmark's cluster counts"),
    ("code", """print(json.dumps(nums.get("mica_default", {}), indent=1))
print("true cluster counts:", sorted(nums["true_k"].values()))"""),
]


def main():
    nb = nbformat.v4.new_notebook()
    nb.cells = [nbformat.v4.new_markdown_cell(src) if kind == "md" else nbformat.v4.new_code_cell(src)
                for kind, src in CELLS]
    nb.metadata["kernelspec"] = {"name": "scminer", "display_name": "scminer", "language": "python"}
    out = REPO / "notebooks"
    out.mkdir(exist_ok=True)
    NotebookClient(nb, timeout=900, kernel_name="scminer", resources={"metadata": {"path": str(out)}}).execute()
    nbformat.write(nb, out / "assignment3_results.ipynb")
    print(f"wrote {out / 'assignment3_results.ipynb'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
