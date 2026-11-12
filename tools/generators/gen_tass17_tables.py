#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0
"""Regenerate the shipped TASS 1.7 Rust tables from the explicit v4 snapshot.

The default input is data/tass17.json beside this script; --input selects another
copy of that fixed realization. No Git checkout, v3 module, engine or network is
used. The historical numeric lineage and MIT attribution are retained; this is
not independent physical validation or permission to replace the theory data.
Every coefficient is emitted as a shortest round-trip binary64 literal. The
current Rust artifact, including its intentional documentation, is reproduced
byte-for-byte. --check compares without writing; --out selects the artifact.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile

DEFAULT_INPUT = Path(__file__).with_name("data") / "tass17.json"
OUTPUT = Path(__file__).resolve().parents[2] / "core/src/tables/tass17_tables.rs"
SOURCE_COMMIT = '37ffb0555b7a9f5083b4f1b6dc1a0f10feeb85fa'
BLOB_SHAS = {
    "libephemeris/moon_theories/tass17.py": "9703f929adaf5c2e8f9723da0f6a5501ce2bf401816dd6173a5baa3d9ef0db19",
    "libephemeris/moon_theories/tass17_data.py": "150ed1172433dc7e915d53042445b70ac9e2528875fd4d77f69bea615da3ab04"
}
CONCATENATED_SHA256 = '29d6a0372aa0a476ec1e1b03acea64f2229e795d56fedc809b199037f3a754f7'
EXPECTED_DATA_SHA256 = '40d14d3a507a78b8c4a6425d062d968f50581a3b3bcf37636547e52d05a4a0f3'
EXPECTED_BODIES = [
    "MIMAS",
    "ENCELADUS",
    "TETHYS",
    "DIONE",
    "RHEA",
    "TITAN",
    "IAPETUS",
    "HYPERION"
]
TERM_RE = re.compile(r"[A-Z]+_[0-3]_[0-9]+")

PREFIX = r"""// SPDX-License-Identifier: MIT
// Copyright (c) 2005 Johannes Gajdosik
// Third-party/adapted data — see THIRD_PARTY_NOTICES.md.
//! TASS 1.7 periodic-term tables (generated, do not edit).
//!
//! Stellarium TASS implementation data (Vienne-Duriez theory),
//! generated from the pinned numeric dataset by
//! `tools/generators/gen_tass17_tables.py` (`--check` regenerates byte-identical).
//! Numbers use shortest round-trip binary64 decimal literals.
//!
//! Theory version: the Vienne-Duriez (1995) paper documents the
//! TASS 1.6 solution; this table is the later IMCCE 1.7 revision
//! (no separate 1.7 publication), as transcribed by Stellarium.
//!
//! Source pin: `37ffb0555b7a9f5083b4f1b6dc1a0f10feeb85fa`
//! (`libephemeris/moon_theories/tass17_data.py`, 1535 periodic
//! terms in 299 term lists, 299 series rows).
//! Blob SHA-256: tass17_data.py `150ed1172433dc7e915d53042445b70ac9e2528875fd4d77f69bea615da3ab04`,
//! concatenated (tass17/tass17_data) `29d6a0372aa0a476ec1e1b03acea64f2229e795d56fedc809b199037f3a754f7`.
//!
//! Ordered layers: per-body per-series [`TassTerm`] leaf lists (amplitude,
//! phase, frequency triples, e.g. `DIONE_0_0`) grouped into the per-series
//! [`TassMultiTerm`] rows (8 index words weighting the seven mean
//! longitudes, e.g. `HYPERION_3`); then the [`tass_series`] body/series
//! accessor (catalog order, Mimas 0 through Hyperion 7);
//! then [`TASS17_BODIES`] (8 body-parameter rows: Saturn-relative `mu`,
//! average-motion scale, constant element seed); then
//! [`TASS17_TO_VSOP87`] (9 row-major rotation words), all consumed in file
//! order by the TASS Saturnian-theory kernels in `crate::planetary_moons`.
//!
//! Frames, units, epochs: none owned — these are table words, not computed
//! frames. Term amplitudes/phases are radians, frequencies are radians per
//! day, `mu` words are AU^3/day^2, seeds are element units, rotation words
//! are dimensionless; the day-count argument and frame rotation belong to
//! the consuming `crate::planetary_moons` kernels, never here.
//!
//! Refusals: none — constants never fail. Out-of-range body/series indices
//! return `None` from [`tass_series`]; unserved bodies and epochs are
//! refused fail-closed at the consuming `crate::planetary_moons` entries,
//! never here.
//!
//! Domains and limits: term, list, row and body counts are fixed by the
//! pinned release blob above; a revised theory needs a regenerated table,
//! never a hand edit here.

/// One periodic term: `(amplitude, phase, frequency)`.
///
/// `Term`: amplitude/phase in radians, frequency in
/// radians per day; the argument is `phase + freq * t + arg`.
pub type TassTerm = (f64, f64, f64);

/// One multiterm row: `(8 index values, term list)`.
///
/// `MultiTerm`: the index values weight the seven mean
/// longitudes (`arg = sum(indices[i] * lon[i])`; the eighth value
/// is carried unchanged and never read by the kernel).
pub type TassMultiTerm = ([i8; 8], &'static [TassTerm]);

/// TASS body parameters: `(name, mu, aam, s0[6])`.
///
/// `TASS17_BODIES` row: `mu` is the Saturn-relative
/// gravitational parameter in AU^3/day^2 units, `aam` the
/// average angular motion scale, `s0` the constant element
/// seed `[n, lambda, k, h, q, p]`. Index order is the model
/// body index (Mimas 0.. Hyperion 7).
pub struct TassBodyParams {
    /// Satellite name.
    /// Kept for row correspondence with the model table;
    /// the kernel indexes positionally, never by name.
    #[allow(dead_code)]
    pub name: &'static str,
    /// Gravitational parameter in AU^3/day^2 units.
    pub mu: f64,
    /// Average angular-motion scale (rad/day).
    pub aam: f64,
    /// Constant element seed `[n, lambda, k, h, q, p]`.
    pub s0: [f64; 6],
}

"""

def f64_lit(value):
    return repr(value)


def render_tables(catalog):
    term_lists = {row["name"]: row["terms"] for row in catalog["term_lists"]}
    series = {row["name"]: row["rows"] for row in catalog["series"]}
    bodies = [(row["name"], row["mu"], row["aam"], row["s0"], row["series"]) for row in catalog["bodies"]]
    rotation = catalog["rotation"]
    lines = PREFIX.splitlines()
    for list_name in sorted(term_lists):
        terms = term_lists[list_name]
        lines.append("#[rustfmt::skip]")
        if len(terms) == 1:
            amp, phase, freq = terms[0]
            flat = (
                f"pub static {list_name}: &[TassTerm] = "
                f"&[({f64_lit(amp)}, {f64_lit(phase)}, {f64_lit(freq)})];"
            )
            if len(flat) <= 100:
                lines.append(flat)
                lines.append("")
                continue
            lines.append(f"pub static {list_name}: &[TassTerm] =")
            lines.append(
                f"    &[({f64_lit(amp)}, {f64_lit(phase)}, {f64_lit(freq)})];"
            )
            lines.append("")
            continue
        lines.append(f"pub static {list_name}: &[TassTerm] = &[")
        for amp, phase, freq in terms:
            flat = f"    ({f64_lit(amp)}, {f64_lit(phase)}, {f64_lit(freq)}),"
            if len(flat) <= 100:
                lines.append(flat)
            else:
                lines.append("    (")
                lines.append(f"        {f64_lit(amp)},")
                lines.append(f"        {f64_lit(phase)},")
                lines.append(f"        {f64_lit(freq)},")
                lines.append("    ),")
        lines.append("];")
        lines.append("")

    for body, _mu, _aam, _s0, _names in bodies:
        for s in range(4):
            key = f"{body}_{s}"
            lines.append("#[rustfmt::skip]")
            lines.append(f"pub static {key}: &[TassMultiTerm] = &[")
            for idx, ref in series[key]:
                idx_text = ", ".join(str(v) for v in idx)
                lines.append(f"    ([{idx_text}], {ref}),")
            lines.append("];")
            lines.append("")

    lines.append("/// Body parameters in the model body-index order.")
    lines.append("#[rustfmt::skip]")
    lines.append("pub static TASS17_BODIES: &[TassBodyParams] = &[")
    for name, mu, aam, s0, _names in bodies:
        s0_text = ", ".join(f64_lit(v) for v in s0)
        lines.append("    TassBodyParams {")
        lines.append(f'        name: "{name}",')
        lines.append(f"        mu: {f64_lit(mu)},")
        lines.append(f"        aam: {f64_lit(aam)},")
        # Keep a seed list on one line only when it fits the 100-character width.
        if len(f"        s0: [{s0_text}],") <= 100:
            lines.append(f"        s0: [{s0_text}],")
        else:
            lines.append("        s0: [")
            for v in s0:
                lines.append(f"            {f64_lit(v)},")
            lines.append("        ],")
        lines.append("    },")
    lines.append("];")
    lines.append("")
    lines.append("/// TASS 17-to-VSOP 87 row-major rotation.")
    lines.append("#[rustfmt::skip]")
    lines.append("pub static TASS17_TO_VSOP87: &[f64; 9] = &[")
    # Same cut: the 9-word rotation always lays out vertically.
    rot_flat = "    " + ", ".join(f64_lit(v) for v in rotation) + ","
    if len(rot_flat) <= 100:
        lines.append(rot_flat)
    else:
        for v in rotation:
            lines.append(f"    {f64_lit(v)},")
    lines.append("];")
    lines.append("")
    lines.append("/// Fetch one body's series (`body` 0-7, `series` 0-3).")
    lines.append("///")
    lines.append("/// ; returns `None` for")
    lines.append("/// out-of-range indices instead of panicking.")
    lines.append("pub fn tass_series(body: usize, series: usize) -> Option<&'static [TassMultiTerm]> {")
    lines.append("    match (body, series) {")
    for index, (name, _mu, _aam, _s0, _names) in enumerate(bodies):
        for s in range(4):
            lines.append(f"        ({index}, {s}) => Some({name}_{s}),")
    lines.append("        _ => None,")
    lines.append("    }")
    lines.append("}")
    lines.append("")

    return "\n".join(lines)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def object_fields(value, expected, label):
    require(isinstance(value, dict) and set(value) == set(expected), f"invalid {label} fields")


def finite_word(value, label):
    require(type(value) is float and math.isfinite(value), f"{label} must be a finite binary64 word")


def reject_constant(value):
    raise ValueError(f"nonfinite JSON constant: {value}")


def unique_fields(pairs):
    result = {}
    for name, value in pairs:
        require(name not in result, f"duplicate JSON field: {name}")
        result[name] = value
    return result


def validate_catalog(catalog):
    object_fields(catalog, ("schema_version", "lineage", "term_lists", "series", "bodies", "rotation"), "catalog")
    require(type(catalog["schema_version"]) is int and catalog["schema_version"] == 1, "unsupported schema")
    lineage = catalog["lineage"]
    object_fields(lineage, ("seed_commit", "seed_blobs_sha256", "seed_concatenated_sha256",
                            "license", "copyright", "note"), "lineage")
    require(lineage["seed_commit"] == SOURCE_COMMIT and lineage["seed_blobs_sha256"] == BLOB_SHAS and
            lineage["seed_concatenated_sha256"] == CONCATENATED_SHA256, "incorrect pinned lineage")
    require(lineage["license"] == "MIT" and lineage["copyright"] == "Copyright (c) 2005 Johannes Gajdosik",
            "incorrect TASS data attribution")
    require(isinstance(lineage["note"], str) and lineage["note"].strip(), "missing lineage explanation")

    leaves = catalog["term_lists"]
    require(isinstance(leaves, list) and len(leaves) == 299, "expected 299 term lists")
    names, term_count = [], 0
    for leaf in leaves:
        object_fields(leaf, ("name", "terms"), "term list")
        name = leaf["name"]
        require(isinstance(name, str) and TERM_RE.fullmatch(name) and
                name.rsplit("_", 2)[0] in EXPECTED_BODIES, "invalid term-list name")
        names.append(name)
        terms = leaf["terms"]
        require(isinstance(terms, list) and terms, "empty or invalid term list")
        for term in terms:
            require(isinstance(term, list) and len(term) == 3, "expected amplitude/phase/frequency triple")
            for word in term:
                finite_word(word, "term coefficient")
        term_count += len(terms)
    require(names == sorted(names) and len(set(names)) == 299, "term-list order or uniqueness changed")
    require(term_count == 1535, "expected 1535 periodic terms")

    groups = catalog["series"]
    require(isinstance(groups, list) and len(groups) == 32, "expected 32 body series")
    expected_series = [f"{name}_{axis}" for name in EXPECTED_BODIES for axis in range(4)]
    references, series_names = [], []
    for group in groups:
        object_fields(group, ("name", "rows"), "series")
        series_names.append(group["name"])
        rows = group["rows"]
        require(isinstance(rows, list) and rows, "empty or invalid series")
        for row in rows:
            require(isinstance(row, list) and len(row) == 2, "invalid multiterm row")
            indices, reference = row
            require(isinstance(indices, list) and len(indices) == 8 and
                    all(type(word) is int and -128 <= word <= 127 for word in indices), "indices must be eight i8 words")
            require(isinstance(reference, str) and reference in names, "dangling term-list reference")
            require(reference.startswith(str(group["name"]) + "_"), "reference belongs to another series")
            references.append(reference)
    require(series_names == expected_series, "body-series order changed")
    require(len(references) == 299 and set(references) == set(names), "series references must retain all 299 leaf lists")

    bodies = catalog["bodies"]
    require(isinstance(bodies, list) and len(bodies) == 8, "expected eight bodies")
    for body, expected in zip(bodies, EXPECTED_BODIES, strict=True):
        object_fields(body, ("name", "mu", "aam", "s0", "series"), "body")
        require(body["name"] == expected, "body index order changed")
        for field in ("mu", "aam"):
            finite_word(body[field], field)
            require(body[field] > 0, f"{field} must be positive")
        require(isinstance(body["s0"], list) and len(body["s0"]) == 6, "expected six seed words")
        for word in body["s0"]:
            finite_word(word, "element seed")
        require(body["series"] == [f"{expected}_{axis}" for axis in range(4)], "body-series wiring changed")

    rotation = catalog["rotation"]
    require(isinstance(rotation, list) and len(rotation) == 9, "expected nine row-major rotation words")
    for word in rotation:
        finite_word(word, "rotation")
    # Schema 1 admits exactly this historical numeric realization. A different
    # theory cannot keep these source claims by editing otherwise valid words.
    payload = {name: catalog[name] for name in ("term_lists", "series", "bodies", "rotation")}
    encoded = json.dumps(payload, sort_keys=True, ensure_ascii=False, allow_nan=False, separators=(",", ":")).encode()
    require(hashlib.sha256(encoded).hexdigest() == EXPECTED_DATA_SHA256, "numeric dataset differs from pinned realization")


def load_catalog(path):
    data = Path(path).read_bytes()
    catalog = json.loads(data, parse_constant=reject_constant, object_pairs_hook=unique_fields)
    validate_catalog(catalog)
    return catalog, hashlib.sha256(data).hexdigest()


def write_output(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=".tass17-", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as stream:
            stream.write(text)
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--out", type=Path, default=OUTPUT)
    args = parser.parse_args(argv)
    try:
        require(args.input.resolve() != args.out.resolve(), "input/output paths alias")
        if args.input.exists() and args.out.exists():
            require(not os.path.samefile(args.input, args.out), "input/output paths alias")
        catalog, digest = load_catalog(args.input)
        text = render_tables(catalog)
        if args.check:
            if args.out.read_bytes() != text.encode("utf-8"):
                print(f"{args.out} differs from regenerated output", file=sys.stderr)
                return 1
            print(f"identity: 299 term lists, 1535 terms, 299 series rows; snapshot {digest}")
        else:
            write_output(args.out, text)
            print(f"wrote {args.out} (299 term lists, 1535 terms)")
        return 0
    except (ValueError, OSError, UnicodeError, OverflowError) as error:
        print(f"TASS17 generator: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
