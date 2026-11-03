#!/usr/bin/env python3
"""Generate src/tables/iau2000_tables.rs from Skyfield's nutation data.

Offline coefficient generator; not part of the runtime crate.

Input: the MIT-licensed numerical dataset in ``data/iau2000_nutation.npz``
beside this script. Its exact SHA-256 matches the dataset recorded in the
shipped Rust table. No v3 checkout, NumPy, Skyfield or network is required.
An explicit ``--npz`` must contain the same dataset; a changed realization
requires a separately reviewed source change, not a new header on old data.
Only the numeric coefficient arrays are transcribed; no Skyfield source
text is copied. Floats are emitted with ``repr`` (shortest round-trip, so
correctly-rounded parsing reproduces the exact binary64 words on both
sides); integer multipliers are emitted as exact decimal literals.

Usage:
  python3 tools/generators/gen_iau2000_tables.py [--npz <file>] [--out <file>]
Defaults: the bundled dataset, and src/tables/iau2000_tables.rs relative to
the repo root. Run rustfmt on the output to restore the Rust source layout.

The script verifies the source identity and expected shapes before writing:
  ke1 (14,), ke0_t (33, 14), se0_t_0/se0_t_1 (33,),
  nals_t (678, 5), lunisolar_longitude/obliquity_coefficients (678, 3),
  napl_t (687, 14), nutation_coefficients_longitude/obliquity (687, 2).
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import io
import os
import struct
import sys
import tempfile
import zipfile
from pathlib import Path

DEFAULT_INPUT = Path(__file__).with_name("data") / "iau2000_nutation.npz"
EXPECTED_SHA256 = "db231254da3bcbef24faa994eadb80f378d39460aa407026c156e6627c30f11a"

EXPECTED = {
    "ke1": ("int64", (14,)),
    "ke0_t": ("int64", (33, 14)),
    "se0_t_0": ("float64", (33,)),
    "se0_t_1": ("float64", (33,)),
    "nals_t": ("int64", (678, 5)),
    "lunisolar_longitude_coefficients": ("float64", (678, 3)),
    "lunisolar_obliquity_coefficients": ("float64", (678, 3)),
    "napl_t": ("int64", (687, 14)),
    "nutation_coefficients_longitude": ("float64", (687, 2)),
    "nutation_coefficients_obliquity": ("float64", (687, 2)),
}


def read_npy(data: bytes):
    """Parse a NumPy .npy v1.0 blob (C order only) into (dtype, shape, values)."""
    if len(data) < 10 or data[:6] != b"\x93NUMPY":
        raise ValueError("bad or truncated npy header")
    major, minor = data[6], data[7]
    if (major, minor) != (1, 0):
        raise ValueError(f"unsupported npy version {(major, minor)}")
    (hlen,) = struct.unpack("<H", data[8:10])
    header = ast.literal_eval(data[10 : 10 + hlen].decode("latin1"))
    descr, fortran, shape = header["descr"], header["fortran_order"], header["shape"]
    if fortran is not False:
        raise ValueError("Fortran order not supported")
    if not isinstance(shape, tuple) or any(type(dim) is not int or dim < 0 for dim in shape):
        raise ValueError("invalid npy shape")
    body = data[10 + hlen :]
    if descr in ("<f8", "=f8", "|f8"):
        fmt, size, conv = "<d", 8, float
        kind = "float64"
    elif descr in ("<i8", "=i8", "|i8"):
        fmt, size, conv = "<q", 8, int
        kind = "int64"
    else:
        raise ValueError(f"unsupported descr {descr}")
    count = 1
    for dim in shape:
        count *= dim
    if len(body) != count * size:
        raise ValueError("npy body length differs from shape")
    values = [conv(v[0]) for v in struct.iter_unpack(fmt, body[: count * size])]
    return kind, tuple(shape), values


def fmt_floats(values):
    return ", ".join(repr(v) for v in values)


def fmt_ints(values):
    return ", ".join(f"{v}_i32" if -2**31 <= v < 2**31 else str(v) for v in values)


def emit_array_f64(name: str, rows: int, cols: int | None, values) -> str:
    body = fmt_floats(values)
    if cols is None:
        return (
            f"pub const {name}: [f64; {rows}] = [\n    {body}\n];\n"
        )
    lines = []
    for r in range(rows):
        chunk = values[r * cols : (r + 1) * cols]
        lines.append("    [" + ", ".join(repr(v) for v in chunk) + "],")
    return f"pub const {name}: [[f64; {cols}]; {rows}] = [\n" + "\n".join(lines) + "\n];\n"


def emit_array_i32(name: str, rows: int, cols: int | None, values) -> str:
    if cols is None:
        return f"pub const {name}: [i32; {rows}] = [\n    {fmt_ints(values)}\n];\n"
    lines = []
    for r in range(rows):
        chunk = values[r * cols : (r + 1) * cols]
        lines.append("    [" + ", ".join(f"{v}_i32" for v in chunk) + "],")
    return f"pub const {name}: [[i32; {cols}]; {rows}] = [\n" + "\n".join(lines) + "\n];\n"


def generate(raw: bytes) -> str:
    """Serialize only the pinned numerical realization, before any output I/O."""
    digest = hashlib.sha256(raw).hexdigest()
    if digest != EXPECTED_SHA256:
        raise ValueError(f"input SHA-256 {digest} differs from {EXPECTED_SHA256}")
    arrays = {}
    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
        for item in zf.infolist():
            if item.filename.endswith(".npy"):
                arrays[item.filename[:-4]] = read_npy(zf.read(item.filename))
    for name, (kind, shape) in EXPECTED.items():
        if name not in arrays:
            raise ValueError(f"missing array {name}")
        got_kind, got_shape, _ = arrays[name]
        if got_kind != kind or got_shape != shape:
            raise ValueError(f"{name}: expected {kind} {shape}, got {got_kind} {got_shape}")

    return render_tables(arrays, digest)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    root = Path(__file__).resolve().parents[2]
    npz_path = (
        Path(args.npz)
        if args.npz
        else DEFAULT_INPUT
    )
    out_path = (
        Path(args.out) if args.out else root / "core" / "src" / "tables" / "iau2000_tables.rs"
    )
    try:
        if npz_path.resolve() == out_path.resolve() or (
            out_path.exists() and os.path.samefile(npz_path, out_path)
        ):
            raise ValueError("input/output paths alias")
        text = generate(npz_path.read_bytes())
        # Publish only a complete serialization; a failed write keeps the prior table.
        descriptor, temporary = tempfile.mkstemp(prefix=".iau2000-", dir=out_path.parent)
        try:
            with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as stream:
                stream.write(text)
            os.replace(temporary, out_path)
        finally:
            Path(temporary).unlink(missing_ok=True)
    except (OSError, ValueError, zipfile.BadZipFile) as error:
        print(f"IAU2000 generator: {error}", file=sys.stderr)
        return 2
    print(f"wrote {out_path} (sha256 {EXPECTED_SHA256})")
    return 0


def render_tables(arrays, digest):
    parts = []
    parts.append(
        "//! IAU 2000A nutation and complementary-term tables (generated, do not edit).\n"
        "//!\n"
        "//! Source: `nutation.npz` bundled with Skyfield 1.54\n"
        "//! (MIT licensed; see `docs/sources.md`). Numeric coefficient tables use\n"
        "//! shortest-round-trip decimal literals, integers exact decimals.\n"
        f"//!\n//! Input SHA-256: `{digest}`.\n"
        "//!\n"
        "//! Generated by `tools/generators/gen_iau2000_tables.py`\n"
        "//!\n"
        "//! `large_const_arrays` lint is allowed for this file.\n"
        "#![allow(clippy::large_const_arrays)]\n"
    )
    order = [
        ("ke1", None),
        ("ke0_t", 14),
        ("se0_t_0", None),
        ("se0_t_1", None),
        ("nals_t", 5),
        ("lunisolar_longitude_coefficients", 3),
        ("lunisolar_obliquity_coefficients", 3),
        ("napl_t", 14),
        ("nutation_coefficients_longitude", 2),
        ("nutation_coefficients_obliquity", 2),
    ]
    const_names = {
        "ke1": "IAU2000_KE1",
        "ke0_t": "IAU2000_KE0_T",
        "se0_t_0": "IAU2000_SE0_T0",
        "se0_t_1": "IAU2000_SE0_T1",
        "nals_t": "IAU2000_NALS_T",
        "lunisolar_longitude_coefficients": "IAU2000_LUNISOLAR_LON",
        "lunisolar_obliquity_coefficients": "IAU2000_LUNISOLAR_OBL",
        "napl_t": "IAU2000_NAPL_T",
        "nutation_coefficients_longitude": "IAU2000_NUT_LON",
        "nutation_coefficients_obliquity": "IAU2000_NUT_OBL",
    }
    for name, cols in order:
        kind, shape, values = arrays[name]
        const = const_names[name]
        if kind == "int64":
            for v in values:
                if not -2**31 <= v < 2**31:
                    raise ValueError(f"{name}: multiplier {v} outside i32")
            ivalues = [int(v) for v in values]
            parts.append(emit_array_i32(const, shape[0], cols, ivalues))
        else:
            parts.append(emit_array_f64(const, shape[0], cols, [float(v) for v in values]))
    return "\n".join(parts)


if __name__ == "__main__":
    sys.exit(main())
