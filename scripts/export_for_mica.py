"""Build an h5ad the authors' MICA can read. Run inside the mica-legacy env.

Current anndata writes files that anndata 0.8 cannot open (MICA issue #33), so
prepare_dataset.py writes a plain exchange folder and this script turns it into
an old-format h5ad.

Usage: conda run -n mica-legacy python scripts/export_for_mica.py zeisel
"""
import sys
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
from scipy import sparse

PROCESSED = Path(__file__).resolve().parents[1] / "data" / "processed"


def main(name, out=None):
    exch = PROCESSED / f"{name}_exchange"
    X = sparse.load_npz(exch / "X.npz")
    cells = (exch / "cells.txt").read_text().splitlines()
    genes = (exch / "genes.txt").read_text().splitlines()
    labels = pd.read_csv(exch / "labels.csv", index_col=0)

    adata = anndata.AnnData(
        X=np.asarray(X.todense(), dtype=np.float32),
        obs=labels.loc[cells],
        var=pd.DataFrame(index=genes),
    )
    out = Path(out) if out else PROCESSED / f"{name}_mica.h5ad"
    adata.write(out)
    print(f"wrote {out}: {adata.shape[0]} cells x {adata.shape[1]} genes, anndata {anndata.__version__}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None))
