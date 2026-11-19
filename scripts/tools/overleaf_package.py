"""Package latex/ for a manual upload to Overleaf.

Overleaf does not sync with this repository, so files move by hand. This script
compares latex/ with an export of the Overleaf project and writes:

  <out>/changed/    only the files that are new or differ, in their folders
  <out>/full.zip    the whole project, for Overleaf's "New project > Upload project"
  <out>/MANIFEST.txt  what to upload into which Overleaf folder

Usage:
  python scripts/tools/overleaf_package.py --overleaf "latex/updated_tex_Nov-13-2026/3. Project Revision and Further Results" --out scratchpad/overleaf_upload
"""

import argparse
import filecmp
import shutil
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LATEX = ROOT / "latex"
KEEP_TOP = ["main.tex", "IEEEtran.cls", "references.bib", "numbers_a4.tex", "numbers_large.tex"]
KEEP_DIRS = ["sections", "tables", "figures"]
SKIP_SUFFIX = {".aux", ".log", ".out", ".bbl", ".blg", ".synctex.gz", ".toc"}


def project_files():
    files = [LATEX / name for name in KEEP_TOP if (LATEX / name).exists()]
    for d in KEEP_DIRS:
        files += [p for p in sorted((LATEX / d).rglob("*")) if p.is_file() and p.suffix not in SKIP_SUFFIX]
    return [p.relative_to(LATEX) for p in files]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--overleaf", required=True, help="folder with the current Overleaf export")
    ap.add_argument("--out", default="scratchpad/overleaf_upload")
    args = ap.parse_args()
    overleaf = Path(args.overleaf)
    if not overleaf.is_absolute():
        overleaf = ROOT / overleaf
    out = Path(args.out)
    if not out.is_absolute():
        out = ROOT / out
    if out.exists():
        shutil.rmtree(out)
    (out / "changed").mkdir(parents=True)

    rel = project_files()
    changed, new = [], []
    for r in rel:
        theirs = overleaf / r
        if not theirs.exists():
            new.append(r)
        elif not filecmp.cmp(LATEX / r, theirs, shallow=False):
            changed.append(r)
        else:
            continue
        dest = out / "changed" / r
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(LATEX / r, dest)

    only_overleaf = sorted(
        p.relative_to(overleaf) for p in overleaf.rglob("*")
        if p.is_file() and p.relative_to(overleaf) not in set(rel)
    )

    with zipfile.ZipFile(out / "full.zip", "w", zipfile.ZIP_DEFLATED) as z:
        for r in rel:
            z.write(LATEX / r, r.as_posix())

    lines = [f"Upload into the existing Overleaf project ({len(changed)} changed, {len(new)} new).",
             "Select the folder in Overleaf's file tree, click Upload, pick the files, and confirm overwrite.", ""]
    groups = {}
    for r in sorted(changed + new):
        groups.setdefault(r.parent.as_posix(), []).append((r.name, "new" if r in new else "changed"))
    for folder, items in groups.items():
        lines.append(f"[{'project root' if folder == '.' else folder + '/'}]")
        lines += [f"  {name}  ({kind})" for name, kind in items]
        lines.append("")
    if only_overleaf:
        lines.append("Only in Overleaf, not in the repo (leave them or delete them by hand):")
        lines += [f"  {r.as_posix()}" for r in only_overleaf]
    (out / "MANIFEST.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
