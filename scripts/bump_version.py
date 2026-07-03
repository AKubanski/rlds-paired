#!/usr/bin/env python3
"""Synchronise the release version across all metadata files.

Single source of truth for the *number* is whatever you pass on the command
line (the release workflow passes the GitHub Release tag); this script then
writes it into every file that needs it, so they never drift:

    rlds_paired/_version.py   __version__           -> PyPI / runtime version
    .zenodo.json              "version"             -> Zenodo record version
    CITATION.cff              version, date-released -> GitHub citation box

Why this matters: Zenodo archives the *committed* tree of the release tag and
reads ``.zenodo.json`` from it. If the version there is stale, every Zenodo
record gets the same version label. The release workflow runs this BEFORE
committing and tagging, so the archived tree always carries the right version.

Usage:
    python scripts/bump_version.py 1.0.1
"""

from __future__ import annotations

import datetime
import json
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print("usage: bump_version.py X.Y.Z", file=sys.stderr)
        return 2

    version = argv[1].lstrip("v").strip()
    if not re.fullmatch(r"\d+\.\d+\.\d+([.\-+].+)?", version):
        print(
            f"error: '{version}' is not a valid version (expected X.Y.Z)",
            file=sys.stderr,
        )
        return 2

    today = datetime.date.today().isoformat()

    # 1) _version.py
    vf = ROOT / "rlds_paired" / "_version.py"
    txt = vf.read_text(encoding="utf-8")
    new_txt, n = re.subn(
        r'__version__\s*=\s*["\'][^"\']*["\']',
        f'__version__ = "{version}"',
        txt,
    )
    if n != 1:
        print("error: could not locate __version__ in _version.py", file=sys.stderr)
        return 1
    vf.write_text(new_txt, encoding="utf-8")

    # 2) .zenodo.json
    zf = ROOT / ".zenodo.json"
    data = json.loads(zf.read_text(encoding="utf-8"))
    data["version"] = version
    zf.write_text(
        json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    # 3) CITATION.cff (version + date-released)
    cf = ROOT / "CITATION.cff"
    c = cf.read_text(encoding="utf-8")
    c, nv = re.subn(r"(?m)^version:.*$", f'version: "{version}"', c)
    c, nd = re.subn(r"(?m)^date-released:.*$", f'date-released: "{today}"', c)
    if nv != 1 or nd != 1:
        print(
            "error: could not update version/date-released in CITATION.cff",
            file=sys.stderr,
        )
        return 1
    cf.write_text(c, encoding="utf-8")

    print(f"bumped version -> {version}  (date-released {today})")
    print("  updated: rlds_paired/_version.py, .zenodo.json, CITATION.cff")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
