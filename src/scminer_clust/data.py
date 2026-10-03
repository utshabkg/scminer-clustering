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
            if str(i) not in g["data"]:
                continue  # nested columns (a data frame inside a column) are stored elsewhere and skipped
            v = g[f"data/{i}"][:]
            cols[name] = np.array([x.decode() for x in v]) if v.dtype.kind == "S" else v
        index = [x.decode() for x in g["row_names"][:]]
    return pd.DataFrame(cols, index=index)


def _load_takane_sce(d):
    """Matrix plus cell and gene tables from a Bioconductor scRNAseq download.

    The matrix is stored genes x cells, so it is transposed to cells x genes.
    """
    with h5py.File(d / "matrix.h5") as h:
        g = h["compressed_sparse_matrix"]
        n_genes, n_cells = (int(x) for x in g["shape"][:])
        counts = sparse.csc_matrix(
            (g["data"][:].astype(np.float32), g["indices"][:].astype(np.int32), g["indptr"][:].astype(np.int64)),
            shape=(n_genes, n_cells),
        )
    obs = _read_takane_frame(d / "column_data.h5")
    var = _read_takane_frame(d / "row_data.h5")
    return ad.AnnData(X=counts.T.tocsr(), obs=obs, var=var)


def load_zeisel() -> ad.AnnData:
    """Zeisel 2015 mouse cortex and hippocampus, 3005 cells.

    Source: Bioconductor scRNAseq collection (zeisel-brain-2015), which mirrors
    the original expression_mRNA_17-Aug-2014.txt. The lab's own download link is
    dead as of October 2026. Labels: level1class (7 broad types) and level2class
    (47 subtypes), both from the original study.
    """
    adata = _load_takane_sce(RAW / "zeisel")
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


def load_usoskin() -> ad.AnnData:
    """Usoskin 2015 mouse sensory neurons, the 622 cells scMINER benchmarks on.

    Source: Bioconductor scRNAseq (usoskin-brain-2015), which holds the reads per
    million released on GEO (GSE59739). We keep the four neuron classes, which give
    exactly the 622 cells and 4 classes in the paper, and drop non-neuronal cells,
    outliers, empty wells and doublets. Known quirks, none of which affect
    clustering: 26 gene names were turned into dates by Excel upstream (March1 to
    March11, Sept1 to Sept15 appear as numbers like 40969), and the authors' own
    MICA input stops at 16,383 genes, probably Excel's column limit.
    """
    adata = _load_takane_sce(RAW / "usoskin")
    adata = adata[adata.obs["Level 1"].isin(["NP", "TH", "NF", "PEP"])].copy()
    # 1,218 rows starting with r_ are RepeatMasker repeats, not genes.
    adata = adata[:, ~adata.var_names.str.startswith("r_")].copy()
    # Recover read counts from reads per million. The authors divided by each
    # cell's total reads (the Reads column), so this is exact for shallow cells and
    # within half a read for deep ones. Count-based gene selection needs counts.
    reads = adata.obs["Reads"].to_numpy(dtype=float)
    adata.X = sparse.csr_matrix(np.round(sparse.diags(reads / 1e6) @ adata.X)).astype(np.float32)
    adata.X.eliminate_zeros()
    adata.obs["label"] = adata.obs["Level 1"].to_numpy()
    adata.obs["label_fine"] = adata.obs["Level 3"].to_numpy()
    adata.obs["label_level2"] = adata.obs["Level 2"].to_numpy()  # the authors' 8-group labels
    adata.uns["dataset"] = "usoskin"
    adata.uns["keep_all_cells"] = True  # the paper used every labelled cell
    return adata


def load_kolodziejczyk() -> ad.AnnData:
    """Kolodziejczyk 2015 mouse embryonic stem cells in three culture conditions.

    Source: Bioconductor scRNAseq (kolodziejczyk-esc-2015), raw read counts with
    spike-ins held separately. The labels are the culture condition each cell was
    grown in (2i, a2i, serum with LIF), so they come from the experiment design,
    not from clustering. 704 cells and 3 classes, as in the paper.
    """
    adata = _load_takane_sce(RAW / "kolodziejczyk")
    # Genes are Ensembl IDs, so flag the mitochondrial genes by ID instead of name.
    mt_ids = {f"ENSMUSG000000{n}" for n in range(64336, 64373)} | {"ENSMUSG00000065947"}
    adata.var["mito"] = adata.var_names.isin(mt_ids)
    adata.obs["label"] = adata.obs["culture"].to_numpy()
    adata.obs["label_fine"] = adata.obs["culture"].to_numpy()
    adata.uns["dataset"] = "kolodziejczyk"
    adata.uns["keep_all_cells"] = True  # the paper used every labelled cell
    return adata


def load_klein() -> ad.AnnData:
    """Klein 2015 mouse embryonic stem cells, days 0, 2, 4 and 7 after LIF removal.

    Source: GEO GSE65525, one inDrop UMI count table per day (genes as rows, no
    header). The label is the day the cell was collected, set by the experiment
    design. 933 + 303 + 683 + 798 = 2,717 cells and 4 classes, as in the paper.
    GEO's sample titles say 935, 301, 682 and 799, which also sum to 2,717. We
    follow the actual file columns, as the SC3 copy the authors used does.
    """
    d = RAW / "klein"
    days = {"d0": "GSM1599494_ES_d0_main", "d2": "GSM1599497_ES_d2_LIFminus",
            "d4": "GSM1599498_ES_d4_LIFminus", "d7": "GSM1599499_ES_d7_LIFminus"}
    parts = []
    for day, stem in days.items():
        df = pd.read_csv(d / f"{stem}.csv.bz2", header=None, index_col=0)
        df = df.loc[:, df.notna().any()]  # guard against a trailing empty column
        a = ad.AnnData(X=sparse.csr_matrix(df.to_numpy(dtype=np.float32).T),
                       obs=pd.DataFrame(index=[f"{day}_{i}" for i in range(df.shape[1])]),
                       var=pd.DataFrame(index=df.index.astype(str).rename(None)))
        a.obs["label"] = day
        parts.append(a)
    adata = ad.concat(parts)
    adata.obs["label_fine"] = adata.obs["label"].to_numpy()
    adata.var_names_make_unique()
    adata.uns["dataset"] = "klein"
    adata.uns["keep_all_cells"] = True  # the paper used every labelled cell
    return adata


LOADERS = {"zeisel": load_zeisel, "pbmc14k": load_pbmc14k, "usoskin": load_usoskin,
           "kolodziejczyk": load_kolodziejczyk, "klein": load_klein}


def load(name: str) -> ad.AnnData:
    return LOADERS[name]()


PAPER_INPUTS = {  # Zenodo 15040179 file stem for each of the paper's ten small and mid-sized benchmarks
    "yan": "01_Yan", "goolam": "02_Goolam", "buettner": "03_Buettner", "pollen": "04_Pollen",
    "chung": "05_Chung", "usoskin": "06_Usoskin", "kolodziejczyk": "07_Kolod", "klein": "08_Klein",
    "zeisel": "09_Zeisel", "pbmc14k": "10_PBMC14K",
}


PAPER_LOG_BASE = {"pollen": 10, "buettner": 10, "chung": np.e, "usoskin": np.e}  # all others log2


def load_paper_input(name) -> ad.AnnData:
    """The authors' own processed MICA input for a benchmark dataset.

    Their Zenodo release (doi 10.5281/zenodo.15040179) has one ExpressionSet per
    dataset, already filtered and log transformed, which is exactly what they fed
    to MICA. Running MICA on it separates "our preprocessing differs" from
    "MICA does not reproduce". The values are used as they are.
    """
    import rdata

    o = rdata.conversion.convert(rdata.parser.parse_file(RAW / "paper_inputs" / f"{PAPER_INPUTS[name]}.Processed_eset.rds"))
    ex = o.assayData["exprs"]
    if hasattr(ex, "x") and not hasattr(ex, "i"):  # dense dgeMatrix, values stored column by column
        n_genes, n_cells = (int(x) for x in ex.Dim)
        X = sparse.csc_matrix(np.asarray(ex.x, dtype=np.float32).reshape(n_cells, n_genes).T)
        genes, cells = [np.asarray(x).astype(str) for x in ex.Dimnames]
    elif hasattr(ex, "x"):  # stored as a sparse dgCMatrix
        n_genes, n_cells = (int(x) for x in ex.Dim)
        X = sparse.csc_matrix((ex.x.astype(np.float32), ex.i.astype(np.int32), ex.p.astype(np.int64)), shape=(n_genes, n_cells))
        genes, cells = [np.asarray(x).astype(str) for x in ex.Dimnames]
    else:  # stored as a dense matrix
        X = sparse.csc_matrix(np.asarray(ex.values, dtype=np.float32))
        genes, cells = (np.asarray(ex.coords[d].values).astype(str) for d in ex.dims)
    pheno = vars(o.phenoData)["data"]
    pheno.columns = [str(c) for c in pheno.columns]
    label_col = next(c for c in ("TrueLabel_4group", "TrueLabel", "trueLabel", "true_label") if c in pheno.columns)
    adata = ad.AnnData(X=X.T.tocsr(), obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))
    adata.obs["label"] = pheno[label_col].astype(str).reindex(adata.obs_names).to_numpy()
    adata.obs["label_fine"] = adata.obs["label"].to_numpy()
    adata.uns["dataset"] = f"paper_{name}"
    adata.uns["label_column"] = label_col
    # The authors' inputs are not all on the same log scale. Checked on Oct 3 by
    # which base makes each cell's values integer steps or CPM sums: Pollen and
    # Buettner are log10(x + 1), Chung and Usoskin are ln(x + 1), and the other six
    # are log2(CPM + 1). Gene selection must know the base, and the approximate
    # untransformed copy used by count-based selection (seurat_v3) is base**x - 1.
    base = PAPER_LOG_BASE.get(name, 2)
    adata.uns["log1p"] = {"base": base}
    counts = adata.X.copy()
    counts.data = np.power(float(base), counts.data) - 1
    adata.layers["counts"] = counts
    return adata


ALLEN_REGIONS = ["CB", "CTXsp", "HPF", "HY", "Isocortex-1", "Isocortex-2", "MB", "MY", "OLF", "P", "PAL", "STR", "TH"]


def load_allen(n_cells=None, seed=0) -> ad.AnnData:
    """Allen Brain Cell Atlas, whole mouse brain, 10x v3 (Yao et al., Nature 2023).

    About 2.3 million cells by 32,285 genes, raw counts in 13 region files.
    Labels come from the atlas taxonomy: class (34) is "label", subclass (338) is
    "label_fine", with supertype and cluster kept too. These labels were made by
    the atlas team's own clustering, so they are not independent ground truth.
    n_cells draws a uniform random subset (for pilots) before reading any counts.
    The atlas was already quality-filtered by its authors, so we keep every cell.
    """
    d = RAW / "allen"
    meta = pd.read_csv(d / "cell_metadata_with_cluster_annotation.csv",
                       usecols=["cell_label", "feature_matrix_label", "library_method", "class", "subclass", "supertype", "cluster"],
                       index_col="cell_label")
    meta = meta[meta["library_method"] == "10Xv3"]
    if n_cells:
        meta = meta.sample(n=n_cells, random_state=seed)
    parts = []
    for r in ALLEN_REGIONS:
        want = meta.index[meta["feature_matrix_label"] == f"WMB-10Xv3-{r}"]
        if len(want) == 0:
            continue
        f = ad.read_h5ad(d / f"WMB-10Xv3-{r}-raw.h5ad", backed="r")
        rows = np.flatnonzero(f.obs_names.isin(want))
        a = f[rows].to_memory()
        a.obs = a.obs[[]]
        parts.append(a)
        f.file.close()
    adata = ad.concat(parts)
    adata.var = parts[0].var[["gene_symbol"]].copy()
    m = meta.reindex(adata.obs_names)
    adata.obs["label"] = m["class"].to_numpy()
    adata.obs["label_fine"] = m["subclass"].to_numpy()
    adata.obs["supertype"] = m["supertype"].to_numpy()
    adata.obs["cluster"] = m["cluster"].to_numpy()
    adata.obs["region"] = m["feature_matrix_label"].str.replace("WMB-10Xv3-", "").to_numpy()
    adata.uns["dataset"] = "allen" if not n_cells else f"allen_{n_cells // 1000}k"
    adata.uns["keep_all_cells"] = True
    return adata


LOADERS["allen_100k"] = lambda: load_allen(100_000)
