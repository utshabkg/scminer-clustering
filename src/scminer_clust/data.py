"""Dataset loaders. Each returns an AnnData of cells x genes with labels in .obs."""
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from .paths import RAW


def _read_takane_frame(path):
    """Read a Bioconductor takane data frame (HDF5) into a pandas DataFrame."""
    with h5py.File(path) as h:
        g = h["data_frame"]
        names = [x.decode() for x in g["column_names"][:]]
        cols = {}
        for i, name in enumerate(names):
            v = g[f"data/{i}"][:]
            cols[name] = np.array([x.decode() for x in v]) if v.dtype.kind == "S" else v
        index = [x.decode() for x in g["row_names"][:]]
    return pd.DataFrame(cols, index=index)


def load_zeisel() -> ad.AnnData:
    """Zeisel 2015 mouse cortex and hippocampus, 3005 cells.

    Source: Bioconductor scRNAseq collection (zeisel-brain-2015), which mirrors
    the original expression_mRNA_17-Aug-2014.txt. The lab's own download link is
    dead as of October 2026. Labels: level1class (7 broad types) and level2class
    (47 subtypes), both from the original study.
    """
    d = RAW / "zeisel"
    with h5py.File(d / "matrix.h5") as h:
        g = h["compressed_sparse_matrix"]
        n_genes, n_cells = (int(x) for x in g["shape"][:])
        counts = sparse.csc_matrix(
            (g["data"][:].astype(np.float32), g["indices"][:].astype(np.int32), g["indptr"][:].astype(np.int64)),
            shape=(n_genes, n_cells),
        )
    obs = _read_takane_frame(d / "column_data.h5")
    var = _read_takane_frame(d / "row_data.h5")
    adata = ad.AnnData(X=counts.T.tocsr(), obs=obs, var=var)
    adata.var_names_make_unique()
    adata.obs = adata.obs.rename(columns={"level1class": "label", "level2class": "label_fine"})
    adata.uns["dataset"] = "zeisel"
    return adata


def load_pbmc14k() -> ad.AnnData:
    """PBMC14K, the main benchmark in the scMINER paper. 14,000 cells, 17,986 genes.

    This is the authors' own raw count matrix, shipped with their R package as
    pbmc14k_rawCount.rda, so our reproduction starts from exactly their numbers.
    The labels come from their PBMC14k_trueLabel.txt.gz. Cells were purified with
    antibody beads before sequencing (Zheng et al. 2017), so the labels do not
    come from clustering.
    """
    import gzip

    import rdata

    d = RAW / "pbmc14k"
    obj = rdata.conversion.convert(rdata.parser.parse_file(d / "pbmc14k_rawCount.rda"))["pbmc14k_rawCount"]
    n_genes, n_cells = (int(x) for x in obj.Dim)
    counts = sparse.csc_matrix(
        (obj.x.astype(np.float32), obj.i.astype(np.int32), obj.p.astype(np.int64)),
        shape=(n_genes, n_cells),
    )
    genes, cells = [np.asarray(x) for x in obj.Dimnames]

    with gzip.open(d / "PBMC14k_trueLabel.txt.gz", "rt") as fh:
        labels = pd.read_csv(fh, sep="\t", index_col=0)

    adata = ad.AnnData(X=counts.T.tocsr(), obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))
    adata.var_names_make_unique()
    adata.obs["label"] = labels["trueLabel"].reindex(adata.obs_names).to_numpy()
    adata.obs["label_fine"] = labels["trueLabel_full"].reindex(adata.obs_names).to_numpy()
    adata.uns["dataset"] = "pbmc14k"
    return adata


LOADERS = {"zeisel": load_zeisel, "pbmc14k": load_pbmc14k}


def load(name: str) -> ad.AnnData:
    return LOADERS[name]()
