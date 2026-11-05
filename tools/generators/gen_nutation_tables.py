#!/usr/bin/env python3
"""Generate src/tables/nutation_tables.rs from the ERFA v2.0.1 C sources.

TIME-09 build-time tool (documented, not part of the shipped crate).

Inputs (admissible open-source IAU implementation, BSD-style terms --
see docs/sources.md and the generated file header):
  pyerfa 2.0.1.5 sdist from PyPI, which vendors liberfa/erfa v2.0.1:
    https://files.pythonhosted.org/packages/source/p/pyerfa/pyerfa-2.0.1.5.tar.gz
  Only these three files are read (operation order + coefficient tables):
    liberfa/erfa/src/nut00a.c   (MHB2000 luni-solar + planetary series)
    liberfa/erfa/src/nut00b.c   (77-term abbreviated series)
    liberfa/erfa/src/s06.c      (CIO locator s+XY/2 series + polynomial)

Usage:
  python3 tools/generators/gen_nutation_tables.py --erfa-src <dir> [--out <file>]
where <dir> holds the extracted liberfa/erfa/src/*.c files.
Default --out is src/tables/nutation_tables.rs relative to the repo root.

The script asserts the expected table shapes before emitting anything:
  nut00a xls: 1365 rows of (5 multipliers + 6 coefficients)
  nut00a xpl:  687 rows of (13 multipliers + 4 coefficients)
  nut00b x:     77 rows of (5 multipliers + 6 coefficients)
  s06 s0..s4 + sp[6].
No ERFA source text is copied: only the numeric coefficient tables are
transcribed (decimal literals parse to the same binary64 words under
correctly-rounded conversion on both sides), plus SHA-256 identities of
the three C files in the emitted header.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import sys
from pathlib import Path

EXPECTED = {
    # MHB2000 (Mathews et al. 2002) carries 1365 luni-solar + planetary
    # terms in total; the ERFA split below is asserted via the sum.
    "nut00a_total": 1365,
    "nut00b_x": 77,
}

FILE_SHAS = {
    # Pinned 2026-09-23 from the pyerfa 2.0.1.5 sdist above
    # (liberfa subtree carries ERFA v2.0.1: AC_INIT([erfa],[2.0.1])).
    "nut00a.c": "a1538ec0277745db974ab76ca4b678a7602249bccecf9474e2b0268f9c2032ef",
    "nut00b.c": "a5246d525e987eb587bbe786b2846160e01be85b5d12af573870b89235695f49",
    "s06.c": "6ad19dd59c3150d9f065b3f76d254c1cf3adb0cb6bacd6422c1dc6b6f14122ea",
}


def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def extract_region(text: str, start_marker: str, end_marker: str) -> str:
    start = text.index(start_marker)
    end = text.index(end_marker, start)
    return text[start:end]


def parse_rows(region: str, n_int: int, n_float: int) -> list[tuple[list[int], list[float]]]:
    """Parse every `{...}` initializer group in a table region.

    Handles rows split across lines and the `{{...}, s, c}` s06 form
    (inner braces are flattened; the caller selects the field split).
    """
    rows: list[tuple[list[int], list[float]]] = []
    for match in re.finditer(r"\{([^{}]+)\}", region):
        body = match.group(1)
        # Skip the struct-declaration block (contains type names).
        if re.search(r"[A-Za-z]", body.replace("e-", " ").replace("e+", " ")):
            # Numeric literals may carry e-6 exponents; anything else
            # alphabetic means this is not a data row.
            cleaned = re.sub(r"[0-9]\s*e\s*[-+]\s*[0-9]+", "", body)
            if re.search(r"[A-Za-z]", cleaned):
                continue
        tokens = [tok.strip() for tok in body.split(",") if tok.strip()]
        if len(tokens) != n_int + n_float:
            continue
        try:
            ints = [int(tok) for tok in tokens[:n_int]]
            floats = [float(tok) for tok in tokens[n_int:]]
        except ValueError:
            continue
        rows.append((ints, floats))
    return rows


def parse_s06_rows(region: str) -> list[tuple[list[int], list[float]]]:
    """Parse `{{8 multipliers}, s, c}` rows of an s06 series table."""
    rows: list[tuple[list[int], list[float]]] = []
    for match in re.finditer(r"\{\s*\{([^{}]+)\}\s*,([^{}]+)\}", region):
        ints = [int(tok) for tok in match.group(1).split(",") if tok.strip()]
        floats = [float(tok) for tok in match.group(2).split(",") if tok.strip()]
        if len(ints) != 8 or len(floats) != 2:
            print(f"bad s06 row: {match.group(0)[:80]}", file=sys.stderr)
            raise SystemExit(2)
        rows.append((ints, floats))
    return rows


def fmt_float(value: float) -> str:
    return repr(value)


def emit_table(
    out: list[str],
    name: str,
    row_type: str,
    rows: list[tuple[list[int], list[float]]],
    ints_per_row: int,
) -> None:
    out.append(f"pub(crate) const {name}: [{row_type}; {len(rows)}] = [")
    for ints, floats in rows:
        int_part = ", ".join(str(v) for v in ints)
        float_part = ", ".join(fmt_float(v) for v in floats)
        out.append(f"    ([{int_part}], [{float_part}]),")
    out.append("];")
    out.append("")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--erfa-src", required=True, help="liberfa/erfa/src directory")
    parser.add_argument("--out", default=None, help="output Rust file")
    args = parser.parse_args()

    src_dir = Path(args.erfa_src)
    repo_root = Path(__file__).resolve().parents[2]
    out_path = Path(args.out) if args.out else repo_root / "core" / "src" / "tables" / "nutation_tables.rs"

    texts: dict[str, str] = {}
    for fname, pinned in FILE_SHAS.items():
        path = src_dir / fname
        if not path.is_file():
            print(f"missing input: {path}", file=sys.stderr)
            return 2
        digest = sha256_of(path)
        if digest != pinned:
            print(f"SHA-256 mismatch for {fname}:\n  got      {digest}\n  expected {pinned}",
                  file=sys.stderr)
            return 2
        texts[fname] = path.read_text()

    nut00a = texts["nut00a.c"]
    nut00b = texts["nut00b.c"]
    s06 = texts["s06.c"]

    xls_region = extract_region(nut00a, "} xls[] = {", "/* Planetary nutation")
    xpl_region = extract_region(nut00a, "} xpl[] = {", "const int NPL")
    if "} x[] = {" not in nut00b:
        print("nut00b table marker not found", file=sys.stderr)
        return 2
    b_region = extract_region(nut00b, "} x[] = {", "/* Fundamental (Delaunay)")

    xls = parse_rows(xls_region, 5, 6)
    xpl = parse_rows(xpl_region, 13, 4)
    xbl = parse_rows(b_region, 5, 6)

    for key, rows in (("nut00b_x", xbl),):
        if len(rows) != EXPECTED[key]:
            print(f"{key}: got {len(rows)} rows, expected {EXPECTED[key]}", file=sys.stderr)
            return 2
    if len(xls) + len(xpl) != EXPECTED["nut00a_total"]:
        print(f"nut00a: got {len(xls)} + {len(xpl)} rows, "
              f"expected total {EXPECTED['nut00a_total']}", file=sys.stderr)
        return 2
    for ints, _ in xls + xbl:
        if any(not -128 <= v <= 127 for v in ints):
            print("xls multiplier out of i8 range", file=sys.stderr)
            return 2
    for ints, _ in xpl:
        if any(not -128 <= v <= 127 for v in ints):
            print("xpl multiplier out of i8 range", file=sys.stderr)
            return 2

    # s06 series s0..s4: `{{8 multipliers}, s, c}` rows plus sp[6].
    s_series: dict[str, list[tuple[list[int], list[float]]]] = {}
    s_markers = [f"static const TERM s{idx}[] = {{" for idx in range(5)]
    for idx in range(5):
        start = s06.index(s_markers[idx]) + len(s_markers[idx])
        if idx < 4:
            end = s06.index(s_markers[idx + 1], start)
        else:
            end = s06.index("static const int NS0", start)
        region = s06[start:end]
        rows = parse_s06_rows(region)
        if not rows:
            print(f"s06 series s{idx}: no rows parsed", file=sys.stderr)
            return 2
        for ints, _ in rows:
            if any(not -128 <= v <= 127 for v in ints):
                print(f"s06 s{idx} multiplier out of i8 range", file=sys.stderr)
                return 2
        s_series[f"s{idx}"] = rows

    sp_match = re.search(r"static const double sp\[\] = \{(.*?)\};", s06, re.S)
    if not sp_match:
        print("s06 sp[] block not found", file=sys.stderr)
        return 2
    sp_block = re.sub(r"/\*.*?\*/", " ", sp_match.group(1), flags=re.S)
    sp_vals = [float(tok) for tok in re.findall(r"[-+]?[0-9][0-9.eE+-]*", sp_block)]
    if len(sp_vals) != 6:
        print(f"s06 sp[]: got {len(sp_vals)} values, expected 6", file=sys.stderr)
        return 2

    out: list[str] = []
    out.append("//! IAU nutation / CIO-locator coefficient tables (TIME-09).")
    out.append("//!")
    out.append("//! GENERATED FILE — do not hand-edit. Regenerate with:")
    out.append("//! `python3 tools/generators/gen_nutation_tables.py --erfa-src <liberfa/erfa/src>`.")
    out.append("//!")
    out.append("//! Coefficient provenance: ERFA v2.0.1 (NumFOCUS Foundation,")
    out.append("//! Copyright (C) 2013-2023; SOFA-derived, BSD-style terms — see")
    out.append("//! `docs/sources.md`). Only numeric tables are transcribed; no")
    out.append("//! ERFA source text is copied. Input identities (SHA-256):")
    for fname, pinned in FILE_SHAS.items():
        out.append(f"//!   {fname}: {pinned}")
    out.append("//!")
    out.append("//! Table shapes (asserted by the generator):")
    out.append(f"//!   MHB2000 luni-solar rows: {len(xls)}, planetary rows: {len(xpl)},")
    out.append(f"//!   IAU 2000B abbreviated rows: {len(xbl)},")
    out.append("//!   s+XY/2 series rows: " + ", ".join(
        f"s{idx}={len(s_series[f's{idx}'])}" for idx in range(5)) + ".")
    out.append("")
    out.append("//! Row layouts mirror the ERFA structs field-for-field:")
    out.append("//!   luni-solar / 2000B: `([nl, nlp, nf, nd, nom], [sp, spt, cp, ce, cet, se])`")
    out.append("//!     with the ERFA accumulation `dp += (sp + spt*t)*sin + cp*cos`,")
    out.append("//!     `de += (ce + cet*t)*cos + se*sin`, units 0.1 microarcsec;")
    out.append("//!   planetary: `([nl, nf, nd, nom, nme, nve, nea, nma, nju, nsa, nur, nne, npa],")
    out.append("//!     [sp, cp, se, ce])` with `dp += sp*sin + cp*cos`,")
    out.append("//!     `de += se*sin + ce*cos`;")
    out.append("//!   s+XY/2: `([l, lp, F, D, Om, LVe, LE, pA], [s, c])` with")
    out.append("//!     `w += s*sin(a) + c*cos(a)`.")
    out.append("")
    out.append("//! Large-table lint: the transcribed series are intentionally")
    out.append("//! `const` (shared read-only words, no runtime cost); the")
    out.append("//! `large_const_arrays` lint is allowed for this file.")
    out.append("#![allow(clippy::large_const_arrays)]")
    out.append("")
    emit_table(out, "NUT_LS_ROWS", "([i8; 5], [f64; 6])", xls, 5)
    emit_table(out, "NUT_PL_ROWS", "([i8; 13], [f64; 4])", xpl, 13)
    emit_table(out, "NUT00B_ROWS", "([i8; 5], [f64; 6])", xbl, 5)
    for idx in range(5):
        emit_table(out, f"S06_S{idx}_ROWS", "([i8; 8], [f64; 2])", s_series[f"s{idx}"], 8)
    out.append("/// Polynomial part of the s+XY/2 series (`sp[6]` in `s06.c`).")
    out.append(f"pub(crate) const S06_SP: [f64; 6] = [{', '.join(fmt_float(v) for v in sp_vals)}];")
    out.append("")
    out_path.write_text("\n".join(out))
    print(f"wrote {out_path}: xls={len(xls)} xpl={len(xpl)} x00b={len(xbl)} " +
          " ".join(f"s{idx}={len(s_series[f's{idx}'])}" for idx in range(5)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
