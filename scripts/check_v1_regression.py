"""Check that protocol v1 reproduces the saved original-protocol results (criterion R0).

Usage: python scripts/check_v1_regression.py [results_root]   (default results/v1)
For every old-layout folder protocol v1 wrote (rule_stability_* and rule_matched_stability_*),
compare it with the folder of the same name under results/, the Assignment 3 runs. Every
column but `seconds` must be equal, value for value, on the rows both files have: the two
picks, every stability_curve row and every full_sweep row. The original-grid files from Oct 22
predate the full_sweep rows, which is reported, not failed. Seconds are always reconstructed in
v1, so they are not compared.
Writes <results_root>/regression.csv and exits non-zero on any mismatch.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from scminer_clust.paths import RESULTS

KEY = ["method", "k_setting", "setting"]


def norm(df):
    df = df.copy()
    if "setting" not in df:
        df["setting"] = np.nan
    df["setting"] = pd.to_numeric(df["setting"], errors="coerce").round(6)
    df.loc[df.k_setting.isin(["stability_finest", "stability_most"]), "setting"] = -1.0
    return df


def compare(new_dir, old_dir):
    new, old = norm(pd.read_csv(new_dir / "metrics.csv")), norm(pd.read_csv(old_dir / "metrics.csv"))
    out = {"folder": new_dir.name, "rows_new": len(new), "rows_old": len(old)}
    if "pbmc14k" in new_dir.name:
        # MICA's node2vec and Louvain take no seed, so MICA cannot reproduce, and the matched
        # candidates depend on MICA's full-data sweep. Only original-grid Leiden is expected to.
        if new_dir.name.startswith("rule_matched_"):
            out["picks"] = "not expected to reproduce (MICA has no seed)"
            return out
        new, old = new[new.method == "leiden_pca"], old[old.method == "leiden_pca"]
        out["note"] = "Leiden only"
    for kind in ("picks", "stability_curve", "full_sweep"):
        sel = (lambda d: d.k_setting.isin(["stability_finest", "stability_most"])) if kind == "picks" \
            else (lambda d, k=kind: d.k_setting == k)
        a, b = new[sel(new)].set_index(KEY).sort_index(), old[sel(old)].set_index(KEY).sort_index()
        if not len(b):
            out[kind] = "absent in saved file"
            continue
        if len(a) != len(b) or not a.index.equals(b.index):
            out[kind] = f"different rows: {len(a)} new against {len(b)} saved"
            continue
        cols = [c for c in b.columns if c in a.columns and c != "seconds" and b[c].notna().any()]
        bad = []
        for c in cols:
            x, y = a[c], b[c]
            xn, yn = pd.to_numeric(x, errors="coerce"), pd.to_numeric(y, errors="coerce")
            if xn.notna().all() and yn.notna().all():
                same = np.array_equal(xn.to_numpy(float), yn.to_numpy(float))
            else:
                same = (x.astype(str).str.replace(r"\.0$", "", regex=True)
                        == y.astype(str).str.replace(r"\.0$", "", regex=True)).all()
            if not same:
                bad.append(c)
        missing = [c for c in b.columns if c not in a.columns and b[c].notna().any()]
        out[kind] = "equal" if not bad and not missing else f"differs in {bad} missing {missing}"
    return out


def main(root=None):
    root = Path(root) if root else RESULTS / "v1"
    rows = []
    for d in sorted(root.glob("rule_*stability_*")):
        old = RESULTS / d.name
        if not old.exists():
            rows.append({"folder": d.name, "picks": "no saved folder (new seed)"})
            continue
        rows.append(compare(d, old))
    df = pd.DataFrame(rows)
    df.to_csv(root / "regression.csv", index=False)
    print(df.to_string(index=False))
    checked = df[df.picks != "no saved folder (new seed)"]
    bad = checked[[c for c in ("picks", "stability_curve", "full_sweep") if c in checked]].map(
        lambda v: v not in ("equal", "absent in saved file", "not expected to reproduce (MICA has no seed)")
        and not pd.isna(v)).any(axis=1)
    print(f"{len(checked)} folders compared, {int(bad.sum())} with a mismatch")
    return int(bad.sum() > 0)


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
