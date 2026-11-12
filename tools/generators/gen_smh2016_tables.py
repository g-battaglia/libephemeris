#!/usr/bin/env python3
"""Generate `src/tables/smh2016_tables.rs`: the vendored SMH-2016 Delta-T realization.

Replicates the pinned release's default Delta-T construction
(`libephemeris/state.py::_build_enhanced_timescale` at peeled
`37ffb0555b7a9f5083b4f1b6dc1a0f10feeb85fa`, `refs/tags/v3.2.1`), which
merges Skyfield 1.54's bundled historical table (`historic_deltat.npy`,
semi-annual observed 1657-1984) with its bundled IERS daily table
(`iers.npz`, 1973-~2027) and feeds the merged pair to Skyfield's own
`timelib.build_delta_t` (linear `numpy.interp` inside the table span,
`Table-S15.2020.txt` splines + SMH-2016 long-term parabola + hand-tuned
patch splines outside, from bundled `delta_t.npz`).

Skyfield (MIT, https://github.com/skyfielders/python-skyfield) is an
owner-admissible reference for this table, not a Rust runtime
dependency: only its numeric table bytes and the deterministic,
already-pinned construction outputs are transcribed here. No Swiss
Ephemeris material is involved.

The construction runs through Skyfield's own functions on
hash-verified inputs, so the emitted constants are bitwise what the
release computes at runtime. The Rust consumer (`src/time/smh2016.rs`)
reimplements only the *evaluation* (`numpy.interp` linear segments,
cubic-Horner spline rows, the `ut1_jd` triple evaluation), which the
release Delta-T sweep checks bitwise.

Usage (from the repository root)::

    python3 tools/generators/gen_smh2016_tables.py

Regenerating against different Skyfield bytes changes the port's
behavior; the input SHAs below are the table identity. Changes require
scientific qualification in the separate comparison repository.
"""

from __future__ import annotations

import hashlib
import os
import re
import sys

SOURCE_COMMIT = "37ffb0555b7a9f5083b4f1b6dc1a0f10feeb85fa"
SOURCE_TAG = "refs/tags/v3.2.1"
SKYFIELD_VERSION = "1.54"
EXPECTED_SHAS = {
    "iers.npz": "c7d7536d898dfa9f8cd43e8044ff51e108cc8289675a13fee9822010a1c4935c",
    "historic_deltat.npy": "f5346b780b36a0325b1847dc6c0083d66edc7e88b7f648b4c98a67bbd02b5d3f",
    "delta_t.npz": "2d12bd3e789543b78a1f53c8b76ed7fecffdf7e5149cfb6a0aed21a8b3db5ff6",
}

OUT_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "src", "tables", "smh2016_tables.rs",
)


def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    here = os.path.dirname(os.path.realpath(__file__))
    repo_root = os.path.realpath(os.path.join(here, "..", ".."))
    if os.path.realpath(os.getcwd()) != repo_root:
        print(f"run from the repository root ({repo_root})", file=sys.stderr)
        return 2
    try:
        import numpy as np
        import skyfield
        import skyfield.data
        from skyfield.timelib import build_delta_t
    except ImportError as exc:
        print(f"requires numpy + skyfield {SKYFIELD_VERSION}: {exc}", file=sys.stderr)
        return 2
    if skyfield.__version__ != SKYFIELD_VERSION:
        print(
            f"skyfield {skyfield.__version__} != pinned {SKYFIELD_VERSION}; refusing",
            file=sys.stderr,
        )
        return 2
    data_dir = os.path.dirname(os.path.realpath(skyfield.data.__file__))
    for name, want in EXPECTED_SHAS.items():
        got = sha256_file(os.path.join(data_dir, name))
        if got != want:
            print(f"{name}: sha {got} != pinned {want}; refusing", file=sys.stderr)
            return 2

    # --- replicate _build_enhanced_timescale (same ops, same order) ---
    iers = dict(np.load(os.path.join(data_dir, "iers.npz")))
    n = len(iers["tt_jd_minus_arange"])
    iers_tt = iers["tt_jd_minus_arange"] + np.arange(n, dtype=np.float64)
    iers_dt = iers["delta_t_1e7"] / 1e7
    hist = np.load(os.path.join(data_dir, "historic_deltat.npy"))
    assert hist.shape == (2, 656), hist.shape
    mask = hist[0] < iers_tt[0]
    merged_tt = np.concatenate([hist[0][mask], iers_tt])
    merged_dt = np.concatenate([hist[1][mask], iers_dt])
    assert bool(np.all(np.diff(merged_tt) > 0)), "merged table must be strict"
    assert bool(np.all(np.isfinite(merged_tt))) and bool(
        np.all(np.isfinite(merged_dt))
    )

    # --- replicate build_delta_t construction (skyfield's own code) ---
    func = build_delta_t((merged_tt, merged_dt))
    curve_table = np.asarray(func.long_term_function.table, dtype=np.float64)
    assert curve_table.shape[0] == 6, curve_table.shape
    ncols = curve_table.shape[1]
    lower = curve_table[0]
    upper = curve_table[1]
    assert bool(np.all(upper - lower > 0)), "spline widths must be positive"

    # --- emit Rust ---
    def fmt_row(values) -> str:
        return ", ".join(repr(float(v)) for v in values)

    lines = []
    lines.append("//! Vendored SMH-2016 Delta-T realization tables (DT-1).")
    lines.append("//!")
    lines.append("//! GENERATED FILE — do not hand-edit. Regenerate with:")
    lines.append("//! `python3 tools/generators/gen_smh2016_tables.py`.")
    lines.append("//!")
    lines.append("//! Table provenance: the pinned release's default Delta-T")
    lines.append("//! (`libephemeris/state.py::_build_enhanced_timescale` @ peeled")
    lines.append(f"//! `{SOURCE_COMMIT}` (`{SOURCE_TAG}`)) merges Skyfield")
    lines.append(f"//! {SKYFIELD_VERSION}'s bundled historical table with its bundled")
    lines.append("//! IERS daily table and feeds the pair to Skyfield's own")
    lines.append("//! `timelib.build_delta_t`. Skyfield is MIT-licensed; only its")
    lines.append("//! numeric table bytes and the deterministic construction outputs")
    lines.append("//! are transcribed here (no Skyfield source text is copied).")
    lines.append("//! Input identities (SHA-256 of the Skyfield bundled files):")
    for name in ("historic_deltat.npy", "iers.npz", "delta_t.npz"):
        lines.append(f"//!   {name}: {EXPECTED_SHAS[name]}")
    lines.append("//!")
    lines.append("//! Table shapes (asserted by the generator):")
    lines.append(f"//!   merged daily table rows: {len(merged_tt)}")
    lines.append(f"//!     TT span: {merged_tt[0]!r} ..= {merged_tt[-1]!r} d")
    lines.append(f"//!   long-term curve spline rows: {ncols}")
    lines.append(f"//!     J-century span: {lower[0]!r} ..= {upper[-1]!r}")
    lines.append("//!")
    lines.append("//! Evaluation (in `src/smh2016.rs`, checked bitwise by the")
    lines.append("//! release Delta-T sweep): linear `numpy.interp` segments over")
    lines.append("//! the merged table inside its span, else the curve splines at")
    lines.append("//! `J = (tt - 1721045.0) / 365.25` with cubic-Horner rows")
    lines.append("//! `[x0, x1, a3, a2, a1, a0]`; the provider replicates the")
    lines.append("//! release `ut1_jd` triple evaluation at /86400.0 day scale.")
    lines.append("")
    lines.append("//! Large-table lint: the transcribed table is intentionally")
    lines.append("//! `const` (shared read-only words, no runtime cost); the")
    lines.append("//! `large_const_arrays` lint is allowed for this file.")
    lines.append("#![allow(clippy::large_const_arrays)]")
    lines.append("//! Some transcribed data words (e.g. `6.28` seconds near a")
    lines.append("//! zero-crossing) resemble well-known constants; they are")
    lines.append("//! measurements, not approximations, so `approx_constant`")
    lines.append("//! is allowed for this file (values must never be edited).")
    lines.append("#![allow(clippy::approx_constant)]")
    lines.append("")
    lines.append("/// Merged SMH-2016 table: TT Julian days, strictly increasing.")
    lines.append(f"pub(crate) const SMH2016_TABLE_TT: [f64; {len(merged_tt)}] = [")
    for v in merged_tt:
        lines.append(f"    {float(v)!r},")
    lines.append("];")
    lines.append("")
    lines.append("/// Merged SMH-2016 table: Delta T in **seconds**.")
    lines.append(f"pub(crate) const SMH2016_TABLE_DT: [f64; {len(merged_dt)}] = [")
    for v in merged_dt:
        lines.append(f"    {float(v)!r},")
    lines.append("];")
    lines.append("")
    lines.append("/// Long-term curve spline rows `[x0, x1, a3, a2, a1, a0]` in")
    lines.append("/// Julian centuries, evaluated with cubic Horner in")
    lines.append("/// `t = (x - x0) / (x1 - x0)` over the segment selected by")
    lines.append("/// clamped `numpy.interp` of `x` against the `x0` column.")
    lines.append(f"pub(crate) const SMH2016_CURVE_ROWS: [[f64; 6]; {ncols}] = [")
    # rustfmt keeps short array rows on one line and lays long ones out
    # vertically; the 13 pinned rows are 52-65 vs 95-105 chars, so a
    # 80-char cut reproduces its layout exactly and the committed file
    # is `cargo fmt --check` clean straight from the generator.
    for c in range(ncols):
        flat = f"    [{fmt_row(curve_table[:, c])}],"
        if len(flat) <= 80:
            lines.append(flat)
        else:
            lines.append("    [")
            for v in curve_table[:, c]:
                lines.append(f"        {float(v)!r},")
            lines.append("    ],")
    lines.append("];")
    lines.append("")
    text = "\n".join(lines)

    # --- verify the emitted text parses back bitwise ---
    got_tt = [float(v) for v in re.findall(r"^    (-?\d[\d.eE+-]*)[,]$", text, re.M)]
    assert len(got_tt) == len(merged_tt) + len(merged_dt), (len(got_tt), len(merged_tt))
    for want, got in zip(
        list(map(float, merged_tt)) + list(map(float, merged_dt)), got_tt
    ):
        assert want == got and np.signbit(want) == np.signbit(got), (want, got)

    with open(OUT_PATH, "w") as fh:
        fh.write(text)
    print(f"wrote {OUT_PATH} ({len(merged_tt)} table rows, {ncols} curve rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
