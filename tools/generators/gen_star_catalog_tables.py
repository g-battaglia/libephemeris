#!/usr/bin/env python3
"""Generate fixed-star tables from the explicit, self-contained v4 snapshot.

The default JSON input stores 1,447 ICRS/J2000 rows, ordered aliases and lookup
maps. Its lineage records the one-time extraction of the already shipped data;
this serializer neither retrieves public catalogs nor independently revalidates
astronomy. It never reads another checkout, imports an engine or accesses Git.

Catalog order and alias collisions are preserved. Unicode 16.0 decimal-digit
ranges are frozen in the input, so regeneration does not depend on the host
Unicode database. Finite binary64 measurements use round-trip decimal literals.

Usage:
  python3 tools/generators/gen_star_catalog_tables.py --check
  python3 tools/generators/gen_star_catalog_tables.py [--input <json>] [--out <file>]

Exit status: 0 on success, 1 on regeneration drift, 2 on invalid input.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

EXPECTED_ROWS = 1447
EXPECTED_ALIASES = 820
EXPECTED_NAME_TO_HIP = 999
EXPECTED_HIP_FIX = 5
EXPECTED_GREEK = 24
EXPECTED_CONSTELLATION = 174
EXPECTED_ND_RUNS = 71
EXPECTED_UNICODE_VERSION = "16.0.0"
GENERATOR_REL = "tools/generators/gen_star_catalog_tables.py"
CANONICAL_REL = "tools/generators/data/star_catalog.json"
REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = REPO_ROOT / CANONICAL_REL


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON field: {key}")
        result[key] = value
    return result


def _reject_constant(value):
    raise ValueError(f"nonfinite JSON literal: {value}")


def _text(value):
    return (isinstance(value, str) and bool(value) and
            all(ord(ch) >= 32 and ord(ch) != 127 and not 0xD800 <= ord(ch) <= 0xDFFF
                for ch in value))


def _hex(value, length):
    return (isinstance(value, str) and len(value) == length and
            all(ch in "0123456789abcdef" for ch in value))


def _validate_lineage(lineage):
    """Check provenance fields, not the unavailable external catalog contents."""
    fields = {"seed_commit", "seed_blobs_sha256", "seed_concatenated_sha256",
              "unicode_ranges_source", "note", "public_catalog_lineage"}
    if not isinstance(lineage, dict) or set(lineage) != fields:
        raise ValueError("expected the complete canonical lineage schema")
    if not _hex(lineage["seed_commit"], 40):
        raise ValueError("invalid canonical seed commit")
    hashes = lineage["seed_blobs_sha256"]
    if (not isinstance(hashes, dict) or not hashes or
            any(not _text(path) or not _hex(digest, 64) for path, digest in hashes.items()) or
            not _hex(lineage["seed_concatenated_sha256"], 64)):
        raise ValueError("invalid canonical lineage SHA-256")
    sources = lineage["public_catalog_lineage"]
    if (not _text(lineage["note"]) or not _text(lineage["unicode_ranges_source"]) or
            not isinstance(sources, list) or not sources or any(not _text(item) for item in sources)):
        raise ValueError("invalid canonical catalog lineage descriptions")


def _pairs(catalog, key, count, value_type, *, unique=True):
    pairs = catalog[key]
    if not isinstance(pairs, list) or len(pairs) != count:
        raise ValueError(f"{key} requires exactly {count} pairs")
    keys = set()
    for pair in pairs:
        if not isinstance(pair, list) or len(pair) != 2 or not _text(pair[0]):
            raise ValueError(f"invalid {key} pair")
        if value_type is str:
            valid = _text(pair[1])
        else:
            # Retain the -1 sentinel for names without a HIP entry.
            valid = type(pair[1]) is int and -1 <= pair[1] <= 0x7FFFFFFF
        if not valid or (unique and pair[0] in keys):
            raise ValueError(f"invalid or duplicate {key} key/value: {pair[0]}")
        keys.add(pair[0])


def _string_map(value):
    return (isinstance(value, dict) and
            all(_text(key) and _text(item) for key, item in value.items()))


def _validate_public_build(build, row_count):
    """Admit explicit acquisition metadata; hashes do not authenticate an origin."""
    fields = {"builder", "seed_sha256", "sources", "propagation", "selection", "report"}
    if (not isinstance(build, dict) or set(build) != fields or
            # The first spelling belongs to immutable older review snapshots;
            # relocating the builder must not require changing their data hashes.
            build["builder"] not in ("tools/star_catalog_builder.py",
                                     "tools/generators/star_catalog_builder.py") or
            not _hex(build["seed_sha256"], 64)):
        raise ValueError("invalid public-catalog build identity")
    propagation = build["propagation"]
    if (not isinstance(propagation, dict) or set(propagation) !=
            {"function", "from_jd", "to_jd", "extension_sha256"} or
            propagation["function"] != "libephemeris.fixed_stars.propagate_proper_motion" or
            propagation["from_jd"] != 2448349.0625 or propagation["to_jd"] != 2451545.0 or
            not _hex(propagation["extension_sha256"], 64)):
        raise ValueError("invalid native catalog propagation identity")
    sources = build["sources"]
    if not isinstance(sources, list) or len(sources) != 5:
        raise ValueError("public build requires all five source identities")
    keys = set()
    for source in sources:
        if (not isinstance(source, dict) or not _text(source.get("key")) or
                not _text(source.get("input_path")) or not _hex(source.get("sha256"), 64)):
            raise ValueError("invalid public-catalog source identity")
        key = source["key"]
        if key in keys:
            raise ValueError("duplicate public-catalog source")
        keys.add(key)
        if key == "wgsn":
            if (set(source) != {"key", "url", "format", "input_path", "sha256"} or
                    not _text(source["url"]) or not _text(source["format"])):
                raise ValueError("invalid WGSN source metadata")
        elif (set(source) != {"key", "catalog", "columns", "constraints", "row_limit",
                             "units", "source_units", "input_path", "sha256"} or
                not _text(source["catalog"]) or source["row_limit"] != -1 or
                not isinstance(source["columns"], list) or not source["columns"] or
                any(not _text(column) for column in source["columns"]) or
                any(not _string_map(source[field]) for field in
                    ("constraints", "units", "source_units"))):
            raise ValueError("invalid VizieR source metadata")
    if keys != {"hip2", "hipmain", "cross", "xhip", "wgsn"}:
        raise ValueError("public build requires the declared catalog set")
    if build["selection"] != {"hp_limit": 5.2, "round_digits": [9, 9, 7, 7, 3, 2, 3]}:
        raise ValueError("unsupported public catalog selection/rounding policy")
    report = build["report"]
    if (not isinstance(report, dict) or set(report) !=
            {"total", "iau_named", "preserved", "with_rv"} or
            any(type(value) is not int or not 0 <= value <= row_count for value in report.values()) or
            report["total"] != row_count):
        raise ValueError("public build counts do not match its catalog")


def _validate_photometry(annotations, rows):
    """Validate sparse source declarations, not the origin of arbitrary rows."""
    if not isinstance(annotations, list):
        raise ValueError("photometry must be a list of individually qualified entries")
    by_id, seen = {row[0]: row for row in rows}, set()
    fields = {"star_id", "hip_number", "value", "band", "catalogue_field", "statistic",
              "epoch", "variability", "historical_unqualified_value", "source"}
    for item in annotations:
        if not isinstance(item, dict) or set(item) != fields:
            raise ValueError("invalid photometry annotation schema")
        ident = item["star_id"]
        if type(ident) is not int or ident not in by_id or ident in seen:
            raise ValueError("duplicate or absent photometry star ID")
        row = by_id[ident]
        if (type(item["hip_number"]) is not int or item["hip_number"] != row[3] or
                type(item["value"]) is not float or item["value"].hex() != row[10].hex()):
            raise ValueError("photometry identity/value does not match its catalog row")
        # Only the reviewed H5 realization is currently supported. H5 is not
        # the Hp median and establishes neither an epoch nor variability.
        if (item["band"] != "johnson_v" or item["catalogue_field"] != "hipparcos_main_vmag" or
                item["statistic"] != "unspecified" or item["epoch"] is not None or
                item["variability"] != "unknown"):
            raise ValueError("unsupported or inferred photometry qualification")
        historical = item["historical_unqualified_value"]
        if historical is not None and (type(historical) is not float or not math.isfinite(historical)):
            raise ValueError("historical unqualified magnitude must be finite or null")
        source = item["source"]
        if (not isinstance(source, dict) or set(source) !=
                {"url", "input_sha256", "input_value", "round_decimal_places"} or
                not _text(source["url"]) or not source["url"].startswith("https://") or
                not _hex(source["input_sha256"], 64) or
                type(source["input_value"]) is not float or not math.isfinite(source["input_value"]) or
                type(source["round_decimal_places"]) is not int or source["round_decimal_places"] != 3 or
                round(source["input_value"], 3).hex() != item["value"].hex()):
            raise ValueError("invalid photometry source/rounding declaration")
        seen.add(ident)


def load_catalog(path: Path):
    """Read a frozen v1 snapshot or explicit public-built v2 candidate as data."""
    raw = Path(path).read_bytes()
    catalog = json.loads(raw, object_pairs_hook=_unique_object, parse_constant=_reject_constant)
    fields = {"schema_version", "lineage", "unicode_version", "rows", "aliases",
              "name_to_hip", "hip_fix", "greek", "constellations", "decimal_ranges"}
    if not isinstance(catalog, dict):
        raise ValueError("expected the complete catalog schema")
    version = catalog.get("schema_version")
    if type(version) is not int or version not in (1, 2):
        raise ValueError("expected catalog schema version 1 or 2")
    if set(catalog) != fields | ({"build"} if version == 2 else set()) | (
            {"photometry"} if "photometry" in catalog else set()):
        raise ValueError("expected the complete catalog schema")
    _validate_lineage(catalog["lineage"])
    rows = catalog["rows"]
    if not isinstance(rows, list) or not rows or (version == 1 and len(rows) != EXPECTED_ROWS):
        raise ValueError(f"rows requires {EXPECTED_ROWS} entries for v1, a nonempty catalog for v2")
    if version == 2:
        _validate_public_build(catalog["build"], len(rows))
    ids = set()
    names = set()
    for row in rows:
        if (not isinstance(row, list) or len(row) != 11 or
                type(row[0]) is not int or not 1000000 < row[0] <= 0x7FFFFFFF or
                not _text(row[1]) or not _text(row[2]) or
                type(row[3]) is not int or not 0 <= row[3] <= 0x7FFFFFFF or
                any(type(value) is not float or not math.isfinite(value) for value in row[4:])):
            raise ValueError("invalid catalog row: expected IDs, names and seven finite binary64 values")
        if row[0] in ids or row[1].upper() in names:
            raise ValueError("duplicate catalog ID or uppercase traditional name")
        ids.add(row[0])
        names.add(row[1].upper())
    _validate_photometry(catalog.get("photometry", []), rows)
    _pairs(catalog, "aliases", EXPECTED_ALIASES, int, unique=False)
    if any(identifier not in ids for _, identifier in catalog["aliases"]):
        raise ValueError("alias refers to an absent catalog ID")
    _pairs(catalog, "name_to_hip", EXPECTED_NAME_TO_HIP, int)
    _pairs(catalog, "hip_fix", EXPECTED_HIP_FIX, int)
    _pairs(catalog, "greek", EXPECTED_GREEK, str)
    _pairs(catalog, "constellations", EXPECTED_CONSTELLATION, str)
    ranges = catalog["decimal_ranges"]
    if (catalog["unicode_version"] != EXPECTED_UNICODE_VERSION or
            not isinstance(ranges, list) or len(ranges) != EXPECTED_ND_RUNS):
        raise ValueError("expected the frozen Unicode 16.0 decimal grammar (71 ranges)")
    previous = -1
    for pair in ranges:
        if (not isinstance(pair, list) or len(pair) != 2 or
                any(type(value) is not int for value in pair) or
                not previous < pair[0] <= pair[1] <= 0x10FFFF or (pair[1] - pair[0] + 1) % 10):
            raise ValueError("invalid, overlapping or unordered decimal-digit range")
        previous = pair[1]
    return catalog, hashlib.sha256(raw).hexdigest()


def ref_star_key(name: str) -> str:
    """Space-free traditional name; the designation suffix keeps its case."""
    s = "".join(ch for ch in name if not ch.isspace())
    head, sep, tail = s.partition(",")
    return head.lower() + sep + tail


def nomen_key(nomenclature: str) -> str:
    """Space-free nomenclature key for the stable catalog index."""
    return "".join(nomenclature.split())


def rust_str(text: str) -> str:
    # Catalog keys contain no control characters. Escape Rust string syntax
    # while preserving the stored name and designation bytes.
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def rust_hex(value: int) -> str:
    return f"0x{value:X}"


def emit_wrapped_block(out: list[str], doc: list[str], decl: str, tokens: list[str]) -> None:
    # Greedy width packing: tokens joined with single spaces; start a line
    # when the next token would
    # push the line past 99 columns. Fixed-width id words always land
    # 10 per line; variable-width keys pack 7-11 per line.
    out.extend(doc)
    out.append(decl)
    current = "    "
    for token in tokens:
        piece = token if current == "    " else " " + token
        if current != "    " and len(current) + len(piece) > 99:
            out.append(current)
            current = "    " + token
        else:
            current += piece
    out.append(current)
    out.append("];")
    out.append("")


def emit_id_block(out: list[str], doc: list[str], decl: str, ids: list[int]) -> None:
    emit_wrapped_block(out, doc, decl, [str(v) + "," for v in ids])


def emit_str_block(out: list[str], doc: list[str], decl: str, keys: list[str]) -> None:
    emit_wrapped_block(out, doc, decl, [rust_str(k) + "," for k in keys])


def emit_pair_block(
    out: list[str], doc: list[str], decl: str, pairs: list[tuple[str, int]]
) -> None:
    out.extend(doc)
    out.append(decl)
    for key, value in pairs:
        out.append(f"    ({rust_str(key)}, {value}),")
    out.append("];")
    out.append("")


def emit_str_pair_block(
    out: list[str], doc: list[str], decl: str, pairs: list[tuple[str, str]]
) -> None:
    out.extend(doc)
    out.append(decl)
    for key, value in pairs:
        out.append(f"    ({rust_str(key)}, {rust_str(value)}),")
    out.append("];")
    out.append("")


def emit_range_block(
    out: list[str], doc: list[str], decl: str, runs: list[tuple[int, int]]
) -> None:
    out.extend(doc)
    out.append(decl)
    for start, end in runs:
        out.append(f"    ({rust_hex(start)}, {rust_hex(end)}),")
    out.append("];")
    out.append("")


def render_catalog(catalog, input_sha256: str) -> str:
    """Derive stable native lookup arrays without changing catalog measurements."""
    rows, aliases = catalog["rows"], catalog["aliases"]
    name_to_hip, hip_fix = dict(catalog["name_to_hip"]), dict(catalog["hip_fix"])
    greek, constellations = dict(catalog["greek"]), dict(catalog["constellations"])
    nd_runs = catalog["decimal_ranges"]
    ref_sorted_ids = [row[0] for row in sorted(rows, key=lambda r: ref_star_key(r[1]))]
    nomen_sorted_ids = [row[0] for row in sorted(rows, key=lambda r: nomen_key(r[2]))]

    name_to_id_upper = [(row[1].upper(), row[0]) for row in rows]
    aliases_upper: dict[str, int] = {}
    for key, value in aliases:
        aliases_upper[key.upper()] = value

    nomen_upper_keys = sorted({"".join(row[2].split()).upper() for row in rows})

    out: list[str] = [
        "//! Generated fixed-star catalog tables.",
        "#![allow(clippy::approx_constant)]",
        "//! Catalog measurements are data, including values that resemble named constants.",
        "//!",
        f"//! Generated by `{GENERATOR_REL}` from the v4 canonical snapshot." if
        catalog["schema_version"] == 1 else
        f"//! Generated by `{GENERATOR_REL}` from a public-catalog candidate; review required.",
        f"//! Default input: `{CANONICAL_REL}`." if catalog["schema_version"] == 1 else
        "//! Source hashes, queries, units and native propagation identity are recorded in the input JSON.",
        f"//! Canonical input SHA-256: `{input_sha256}`." if catalog["schema_version"] == 1 else
        f"//! Public-catalog input SHA-256: `{input_sha256}`.",
        f"//! Historical seed lineage: `{catalog['lineage']['seed_commit']}`.",
        "//! Sparse photometric source annotations do not qualify the remaining catalog measurements.",
        f"//! Rows: {len(rows)}; alias pairs: {len(aliases)}; name/HIP pairs: {len(name_to_hip)}.",
        "//! Astrometric/name sources include Hipparcos, IAU and WGSN; photometric metadata may be incomplete.",
        "//!",
        "//! Positions are ICRS/J2000 RA/Dec in degrees; proper motion is arcseconds/year;",
        "//! parallax is milliarcseconds; radial velocity is km/s; brightness is magnitudes.",
        "//! Epoch propagation and frame transformations belong to `crate::stars`.",
        "//! Name, alias and nomenclature maps provide the catalog's stable lookup order.",
        f"//! Unicode {EXPECTED_UNICODE_VERSION} decimal-digit ranges support catalog indices.",
        "//! Regenerate the tables rather than editing catalog measurements by hand.",
        "",
    ]
    if catalog.get("photometry"):
        out.append("use crate::stars::{")
        out.append("    PhotometricBand, PhotometricCatalogueField, PhotometricStatistic, PhotometricVariability,")
        out.append("    StarCatalogEntry, StarPhotometry,")
        out.append("};")
    else:
        out.append("use crate::stars::{StarCatalogEntry, StarPhotometry};")
    out.append("")
    out.append("/// Decoded catalog rows.")
    out.append("/// Built through the `const` row constructor so there is exactly one")
    out.append("/// row representation in the crate.")
    out.append("pub(crate) const STAR_CATALOG: &[StarCatalogEntry] = &[")
    for row in rows:
        out.append("    StarCatalogEntry::from_row(")
        out.append(f"        {row[0]},")
        out.append(f"        {rust_str(row[1])},")
        out.append(f"        {rust_str(row[2])},")
        out.append(f"        {row[3]},")
        for i in range(4, 11):
            out.append(f"        {repr(row[i])},")
        out.append("    ),")
    out.append("];")
    out.append("")

    out.append("/// Individually qualified magnitudes; absent entries remain explicitly unqualified.")
    out.append("pub(crate) const STAR_PHOTOMETRY: &[(i32, StarPhotometry)] = &[")
    for item in catalog.get("photometry", []):
        historical = item["historical_unqualified_value"]
        historical = "None" if historical is None else f"Some({historical!r})"
        out.extend([
            "    (",
            f"        {item['star_id']},",
            "        StarPhotometry {",
            f"            value: {item['value']!r},",
            "            band: PhotometricBand::JohnsonV,",
            "            catalogue_field: PhotometricCatalogueField::HipparcosMainVmag,",
            "            statistic: PhotometricStatistic::Unspecified,",
            "            epoch: None,",
            "            variability: PhotometricVariability::Unknown,",
            f"            historical_unqualified_value: {historical},",
            "        },",
            "    ),",
        ])
    out.append("];")
    out.append("")

    emit_id_block(
        out,
        [
            "/// Star ids in order (sorted by the",
            "/// traditional-name key; retained as the stable v 1 catalog index).",
        ],
        "pub(crate) const REF_SORTED_IDS: &[i32] = &[",
        ref_sorted_ids,
    )
    emit_id_block(
        out,
        [
            "/// Star ids in order (sorted by the",
            "/// space-stripped nomenclature; computed by the generator).",
            "/// The v 2 resolver uses this stable catalog index.",
        ],
        "pub(crate) const NOMEN_SORTED_IDS: &[i32] = &[",
        nomen_sorted_ids,
    )
    emit_pair_block(
        out,
        [
            "/// Upper traditional name -> id; names",
            "/// are unique under `.upper()`, asserted by the generator.",
        ],
        "pub(crate) const NAME_TO_ID_UPPER: &[(&str, i32)] = &[",
        name_to_id_upper,
    )
    emit_pair_block(
        out,
        [
            "/// Upper alias -> id; last wins on",
            "/// uppercase collisions.",
        ],
        "pub(crate) const ALIASES_UPPER: &[(&str, i32)] = &[",
        list(aliases_upper.items()),
    )
    emit_str_block(
        out,
        [
            "/// Sorted upper space-stripped nomenclature keys;",
            "/// binary-searched.",
            f"/// Shared designations are deduplicated ({len(nomen_upper_keys)} keys",
            f"/// from {len(rows)} rows; e.g. one designation shared by a bright pair), so the",
            "/// count difference against the catalog is expected.",
        ],
        "pub(crate) const NOMEN_UPPER_KEYS: &[&str] = &[",
        nomen_upper_keys,
    )
    emit_pair_block(
        out,
        [
            "/// Alias -> star-id pairs in canonical insertion order;",
            "/// duplicate keys retain their position.",
            "/// Linear scans preserve the canonical first-hit order.",
        ],
        "pub(crate) const STAR_ALIASES: &[(&str, i32)] = &[",
        aliases,
    )
    emit_pair_block(
        out,
        [
            "/// Upper name -> HIP pairs in canonical insertion order;",
            "/// keys are already uppercase.",
        ],
        "pub(crate) const STAR_NAME_TO_HIP: &[(&str, i32)] = &[",
        list(name_to_hip.items()),
    )
    emit_pair_block(
        out,
        ["/// Curated name -> HIP corrections.",
         "/// Five mappings are preserved from the canonical snapshot, separately",
         "/// from the catalog rows. Serialization does not revalidate their identities."],
        "pub(crate) const NAME_HIP_FIX: &[(&str, i32)] = &[",
        list(hip_fix.items()),
    )
    emit_str_pair_block(
        out,
        [
            "/// Greek letter name -> 2-letter abbreviation;",
            "/// keys are uppercase.",
            "/// Abbreviations compose Bayer nomenclatures (`\"Alpha Leonis\"` -> `\"al\"` +",
            "/// `\"Leo\"` -> `\"alLeo\"`). Omicron abbreviates to `\"omi\"`; only omega",
            "/// abbreviates to `\"om\"`.",
        ],
        "pub(crate) const GREEK_LETTER_ABBREV: &[(&str, &str)] = &[",
        list(greek.items()),
    )
    emit_str_pair_block(
        out,
        [
            "/// Constellation name -> 3-letter abbreviation;",
            "/// keys are uppercase.",
            "/// Genitive and nominative forms both map to the same abbreviation",
            "/// (supporting `\"Alpha Leonis\"` and `\"Alpha Leo\"` alike).",
        ],
        "pub(crate) const CONSTELLATION_ABBREV: &[(&str, &str)] = &[",
        list(constellations.items()),
    )
    emit_range_block(
        out,
        [
            "/// Unicode decimal-digit (Nd) runs `(start, end)` with values",
            "/// `0..=9` ascending from each run start, for native catalog indices.",
        ],
        "pub(crate) const UNICODE_ND_RANGES: &[(u32, u32)] = &[",
        nd_runs,
    )

    return "\n".join(out)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT,
                        help="explicit v1 snapshot or v2 public-catalog JSON (default: shipped snapshot)")
    parser.add_argument("--out", type=Path, default=REPO_ROOT / "core/src/tables/star_catalog_tables.rs")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    try:
        if args.input.resolve() == args.out.resolve() or (args.out.exists() and args.input.samefile(args.out)):
            raise ValueError("output must differ from the canonical input")
        catalog, digest = load_catalog(args.input)
        canonical_output = REPO_ROOT / "core/src/tables/star_catalog_tables.rs"
        protected = (DEFAULT_INPUT, canonical_output) if catalog["schema_version"] == 2 else (DEFAULT_INPUT,)
        if not args.check:
            if any(args.out.resolve() == path.resolve() or
                   (args.out.exists() and path.exists() and args.out.samefile(path)) for path in protected):
                raise ValueError("output must not replace shipped catalog data; use a separate --out file")
        text = render_catalog(catalog, digest)
        current = args.out.read_text(encoding="utf-8") if args.out.is_file() else None
        if args.check:
            if current != text:
                import difflib
                diff = difflib.unified_diff((current or "").splitlines(), text.splitlines(),
                                            str(args.out), "regenerated", lineterm="")
                print("regeneration drift:\n" + "\n".join(list(diff)[:40]), file=sys.stderr)
                return 1
            print(f"identical: {args.out} ({len(catalog['rows'])} rows, "
                  f"{len(catalog['aliases'])} aliases, {len(catalog['decimal_ranges'])} decimal digit runs)")
            return 0
        args.out.write_text(text, encoding="utf-8")
        print(f"wrote {args.out}: {len(catalog['rows'])} rows, {len(catalog['aliases'])} aliases")
        return 0
    except (OSError, UnicodeError, ValueError) as error:
        print(f"invalid catalog input/output: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
