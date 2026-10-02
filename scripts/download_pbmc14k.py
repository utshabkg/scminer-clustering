"""Download PBMC14K, the scMINER paper's main benchmark.

The authors ship their own raw count matrix inside their R package, so we take it
straight from their repository instead of rebuilding it from the 10x sorted
population files. That way the reproduction starts from the same matrix they used.
The file is an R data file, read in Python with the `rdata` package.
"""
import sys
import urllib.request

from scminer_clust.paths import RAW

FILES = {
    "pbmc14k_rawCount.rda": "https://raw.githubusercontent.com/jyyulab/scMINER/main/data/pbmc14k_rawCount.rda",
    "PBMC14k_trueLabel.txt.gz": "https://raw.githubusercontent.com/jyyulab/scMINER/main/inst/extdata/demo_pbmc14k/PBMC14k_trueLabel.txt.gz",
}


def main():
    out = RAW / "pbmc14k"
    out.mkdir(parents=True, exist_ok=True)
    for name, url in FILES.items():
        target = out / name
        if target.exists():
            print(f"have {name} ({target.stat().st_size / 1e6:.1f} MB)")
            continue
        print(f"downloading {name}")
        urllib.request.urlretrieve(url, target)
        print(f"  saved {target.stat().st_size / 1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
