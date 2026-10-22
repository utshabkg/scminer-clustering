"""Check prose against the team's writing rules before anything is handed in.

Usage: python scripts/lint_prose.py [path ...]
With no arguments it checks the manuscript sections, the markdown at the repo root
and everything in docs/. Exits non-zero if anything is flagged, so it can gate a
commit.

The rules: no em or en dashes, no semicolons, no comma placed inside quotation
marks, and none of a short list of stock phrases. In LaTeX a double hyphen renders
as an en dash and a triple as an em dash, so those are flagged too. Markdown table
separators and horizontal rules are skipped, as are BibTeX page ranges, where a
double hyphen is correct.
"""
import re
import sys
from pathlib import Path

from scminer_clust.paths import REPO

STOCK = ["moreover", "furthermore", "in conclusion", "it is worth noting", "delve",
         "leverage", "crucial", "comprehensive"]
PATTERN = re.compile(r"[—–]|;|,\"|,''|(?<![-%!])---?(?![->])|\b(" + "|".join(STOCK) + r")\b", re.I)
SEPARATOR = re.compile(r"[\s|:+-]*$")
DEFAULT = ["latex/sections", "latex/tables", "README.md", "PROGRESS.md", "WORKPLAN.md", "RESUME.md",
           "docs/primer.md", "docs/dataset_cards", "docs/allen_pilot.md", "docs/minutes"]
SKIP_SUFFIX = {".bib", ".csv", ".json", ".cls", ".png", ".pdf", ".docx", ".ipynb", ".xlsx", ".txt", ".js"}
# docs/research holds raw saved literature sweeps, not prose we write or hand in
SKIP_DIRS = {"archive", "research", "verification"}


MATH = re.compile(r"\$[^$]*\$|\\begin\{equation\}.*?\\end\{equation\}|\\\\[;,!:]", re.S)


def strip_math(text):
    """Blank out LaTeX math, where a semicolon is notation rather than punctuation."""
    return MATH.sub(lambda m: " " * len(m.group()), text)


def check(path):
    hits = []
    try:
        text = path.read_text()
    except (UnicodeDecodeError, OSError):
        return hits
    if path.suffix == ".tex":
        text = strip_math(text)
    for i, line in enumerate(text.splitlines(), 1):
        stripped = line.strip()
        if SEPARATOR.fullmatch(stripped) or stripped.startswith("%"):
            continue
        for m in PATTERN.finditer(line):
            hits.append((path, i, m.group(), stripped[:88]))
    return hits


def main(args):
    targets = [Path(a) for a in args] or [REPO / d for d in DEFAULT]
    files = []
    for t in targets:
        if t.is_dir():
            files += [p for p in sorted(t.rglob("*")) if p.is_file() and p.suffix not in SKIP_SUFFIX]
        elif t.is_file() and t.suffix not in SKIP_SUFFIX:
            files.append(t)
    files = [f for f in files if not SKIP_DIRS & set(f.parts)]
    hits = [h for f in files for h in check(f)]
    for path, line, found, context in hits:
        rel = path.relative_to(REPO) if REPO in path.parents else path
        print(f"{rel}:{line}: {found!r}  {context}")
    print(f"checked {len(files)} files: " + ("clean" if not hits else f"{len(hits)} to fix"))
    return 1 if hits else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
