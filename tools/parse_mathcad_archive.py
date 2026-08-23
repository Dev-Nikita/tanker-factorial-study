"""Extract scalar definitions and function bodies from the archived Mathcad worksheets.

The archive contains three XMCD (Mathcad 14 XML) worksheets. They are *different
computational cases*, not versions of one model, so their parameters must never be
mixed. This tool reads each worksheet programmatically and writes

    results/archive/archive_parameters.csv      one row per scalar definition
    results/archive/archive_functions.csv       one row per function definition
    results/archive/archive_manifest.csv        SHA-256 of every source file

so that every parameter used in the paper can be traced to a worksheet and a
region id rather than to a hand transcription.

Usage:  python3 tools/parse_mathcad_archive.py [/path/to/маткад]
"""
from __future__ import annotations

import hashlib
import re
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ARCHIVE = ROOT.parent / "маткад"
OUT = ROOT / "results" / "archive"

ML = "{http://schemas.mathsoft.com/math30}"
WS = "{http://schemas.mathsoft.com/worksheet30}"


# --------------------------------------------------------------------------- #
def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _text(el) -> str:
    return "".join(el.itertext()).strip()


def _to_infix(el) -> str:
    """Render an ml: expression tree as a compact infix string."""
    tag = el.tag.replace(ML, "")
    if tag == "id":
        return (el.text or "").strip()
    if tag == "real":
        return (el.text or "").strip()
    if tag == "str":
        return f'"{(el.text or "").strip()}"'
    if tag == "apply":
        kids = list(el)
        if not kids:
            return ""
        op = kids[0].tag.replace(ML, "")
        args = [_to_infix(k) for k in kids[1:]]
        binary = {"plus": " + ", "minus": " - ", "mult": "*", "div": "/",
                  "pow": "^", "lessThan": " < ", "greaterThan": " > "}
        if op in binary and len(args) == 2:
            return f"({args[0]}{binary[op]}{args[1]})"
        if op == "neg" and len(args) == 1:
            return f"(-{args[0]})"
        if op == "indexer" and len(args) == 2:
            return f"{args[0]}[{args[1]}]"
        if op == "parens" and len(args) == 1:
            return args[0]
        if op == "derivative":
            return f"d/d({','.join(args)})"
        if op == "id":
            return f"{_to_infix(kids[0])}({', '.join(args)})"
        return f"{op}({', '.join(args)})"
    if tag in ("parens", "sequence"):
        return ", ".join(_to_infix(k) for k in el)
    if tag == "function":
        kids = list(el)
        return f"{_to_infix(kids[0])}({', '.join(_to_infix(k) for k in kids[1:])})"
    kids = list(el)
    if len(kids) == 1:
        return _to_infix(kids[0])
    if kids:
        return f"{tag}({', '.join(_to_infix(k) for k in kids)})"
    return (el.text or "").strip()


def parse_worksheet(path: Path) -> tuple[list[dict], list[dict]]:
    tree = ET.parse(path)
    root = tree.getroot()
    scalars, functions = [], []

    for region in root.iter(f"{WS}region"):
        rid = region.get("region-id", "")
        for define in region.iter(f"{ML}define"):
            kids = list(define)
            if len(kids) < 2:
                continue
            lhs, rhs = kids[0], kids[1]

            # ---- function definition: f(args) := body ---------------------- #
            if lhs.tag == f"{ML}function":
                fk = list(lhs)
                name = _to_infix(fk[0])
                args = [_to_infix(k) for k in fk[1:]]
                functions.append({
                    "worksheet": path.name, "region_id": rid,
                    "name": name, "args": ", ".join(args),
                    "body": _to_infix(rhs)[:4000],
                })
                continue

            # ---- scalar definition: name := value --------------------------- #
            if lhs.tag == f"{ML}id":
                name = (lhs.text or "").strip()
                if rhs.tag == f"{ML}real":
                    try:
                        val = float((rhs.text or "").strip())
                    except ValueError:
                        continue
                    scalars.append({"worksheet": path.name, "region_id": rid,
                                    "name": name, "value": val,
                                    "expression": str(val)})
                else:
                    scalars.append({"worksheet": path.name, "region_id": rid,
                                    "name": name, "value": float("nan"),
                                    "expression": _to_infix(rhs)[:400]})
    return scalars, functions


def main():
    archive = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_ARCHIVE
    if not archive.exists():
        raise SystemExit(f"archive not found: {archive}")
    OUT.mkdir(parents=True, exist_ok=True)

    manifest, scalars, functions = [], [], []
    for p in sorted(archive.rglob("*")):
        if p.is_file() and p.suffix.lower() in (".xmcd", ".mcd", ".doc", ".docx", ".html"):
            manifest.append({"file": str(p.relative_to(archive)),
                             "bytes": p.stat().st_size, "sha256": sha256(p)})
    for p in sorted(archive.glob("*.xmcd")):
        s, f = parse_worksheet(p)
        scalars += s
        functions += f
        print(f"{p.name}: {len(s)} scalar definitions, {len(f)} functions")

    pd.DataFrame(manifest).to_csv(OUT / "archive_manifest.csv", index=False)
    ds = pd.DataFrame(scalars)
    ds.to_csv(OUT / "archive_parameters.csv", index=False)
    pd.DataFrame(functions).to_csv(OUT / "archive_functions.csv", index=False)
    print(f"\n-> {OUT}")
    return ds


if __name__ == "__main__":
    main()
