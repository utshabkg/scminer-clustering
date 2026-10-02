"""Download the Zeisel 2015 dataset.

The Linnarsson lab download link (expression_mRNA_17-Aug-2014.txt) returns 403 as
of October 2026, so we take the copy Bioconductor keeps for its scRNAseq package.
It is the same data with the original level1class and level2class labels.
"""
import sys
import urllib.request

from scminer_clust.paths import RAW

BASE = "https://gypsum.artifactdb.com/file/scRNAseq%2Fzeisel-brain-2015%2F2023-12-14%2F"
FILES = {
    "matrix.h5": "assays%2F0%2Fmatrix.h5",
    "column_data.h5": "column_data%2Fbasic_columns.h5",
    "row_data.h5": "row_data%2Fbasic_columns.h5",
}


def main():
    out = RAW / "zeisel"
    out.mkdir(parents=True, exist_ok=True)
    for name, path in FILES.items():
        target = out / name
        if target.exists():
            print(f"have {name} ({target.stat().st_size / 1e6:.1f} MB)")
            continue
        print(f"downloading {name}")
        urllib.request.urlretrieve(BASE + path, target)
        print(f"  saved {target.stat().st_size / 1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
