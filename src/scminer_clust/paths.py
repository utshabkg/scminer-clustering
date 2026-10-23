"""Project paths. Data lives on the shared disk, results in the repo."""
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "data"          # symlink to /media/20TB/shared/sides_data/43.scMINER-Clustering
RAW = DATA / "raw"
PROCESSED = DATA / "processed"
RESULTS = REPO / "results"


def results_dir(exp_id: str) -> Path:
    d = RESULTS / exp_id
    d.mkdir(parents=True, exist_ok=True)
    return d
