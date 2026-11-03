#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0
"""Build a review-only fixed-star candidate from explicit public catalog inputs.

The cache contains hip2.json, hipmain.json, cross.json, xhip.json (normalized
VizieR query tables) and IAU-CSN.txt. Existing cache files are never downloaded
again. Missing files require --fetch and optional astroquery; WGSN acquisition
uses the standard library. No acquisition occurs merely by importing this tool.

The four table caches record catalog IDs, exact queries, declared source units
and normalized units. Local cache provenance is a declaration, not proof that a
file came from the named service. Every actual input is hashed in the candidate.
Rust performs J1991.25-to-J2000 proper-motion propagation. This offline tool
only joins records, converts declared units, applies naming policy and rounds
fields. The output is a schema-v2 review candidate; it never updates the shipped
catalog or Rust table. Run the separate serializer only after reviewing it.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import os
import re
import sys
import tempfile
from pathlib import Path

try:
    from . import gen_star_catalog_tables as serializer
except ImportError:  # Direct invocation from outside the checkout.
    import gen_star_catalog_tables as serializer


HP_LIMIT = 5.2
FROM_JD, TO_JD = 2448349.0625, 2451545.0
IAU_CSN_URL = "https://www.pas.rochester.edu/~emamajek/WGSN/IAU-CSN.txt"
WGSN_FORMAT = "IAU-CSN fixed-width name plus HIP/HD columns"
ROUND_DIGITS = [9, 9, 7, 7, 3, 2, 3]
QUERIES = {
    "hip2": {"catalog": "I/311/hip2",
             "columns": ["HIP", "RArad", "DErad", "Plx", "pmRA", "pmDE", "Hpmag"],
             "constraints": {"Hpmag": "<7.2"}, "row_limit": -1,
             "units": {"RArad": "deg", "DErad": "deg", "Plx": "mas", "pmRA": "mas/yr",
                       "pmDE": "mas/yr", "Hpmag": "mag"}},
    "hipmain": {"catalog": "I/239/hip_main", "columns": ["HIP", "HD", "Vmag"],
                "constraints": {"Vmag": "<7.7"}, "row_limit": -1, "units": {"Vmag": "mag"}},
    "cross": {"catalog": "IV/27A/catalog", "columns": ["HD", "Bayer", "Fl", "Cst"],
              "constraints": {}, "row_limit": -1, "units": {}},
    "xhip": {"catalog": "V/137D/XHIP", "columns": ["HIP", "RV"],
             "constraints": {"RV": "!= NULL"}, "row_limit": -1, "units": {"RV": "km/s"}},
}
GREEK = dict(zip(
    ["alf", "bet", "gam", "del", "eps", "zet", "eta", "the", "iot", "kap", "lam",
     "mu", "nu", "ksi", "xi", "omi", "pi", "rho", "sig", "tau", "ups", "phi", "chi", "psi", "ome"],
    ["al", "be", "ga", "de", "ep", "ze", "et", "th", "io", "ka", "la",
     "mu", "nu", "xi", "xi", "omi", "pi", "rh", "si", "ta", "up", "ph", "ch", "ps", "ome"],
    strict=True,
))


def _digest(path):
    with Path(path).open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def _read_json(path):
    raw = Path(path).read_bytes()
    table = json.loads(raw, object_pairs_hook=serializer._unique_object,
                       parse_constant=serializer._reject_constant)
    return table, hashlib.sha256(raw).hexdigest()


def _validate_table(table, key):
    """Refuse unlabelled/raw-unit tables instead of guessing their source units."""
    query = QUERIES[key]
    if not isinstance(table, dict) or set(table) != {"schema_version", *query, "source_units", "rows"}:
        raise ValueError(f"invalid normalized {key} cache schema")
    if type(table["schema_version"]) is not int or table["schema_version"] != 1:
        raise ValueError(f"unsupported {key} cache schema version")
    if any(table[field] != value for field, value in query.items()):
        raise ValueError(f"{key} cache does not match its public catalog/query/normalized units")
    source_units = table["source_units"]
    if (not isinstance(source_units, dict) or set(source_units) != set(query["units"]) or
            any(not isinstance(unit, str) or not unit for unit in source_units.values())):
        raise ValueError(f"missing declared {key} source units")
    rows = table["rows"]
    if not isinstance(rows, list) or not rows:
        raise ValueError(f"{key} cache must contain a nonempty query result")
    for row in rows:
        if not isinstance(row, dict) or set(row) != set(query["columns"]):
            raise ValueError(f"invalid {key} row columns")
        for field, value in row.items():
            if value is None:
                continue
            if field in ("Bayer", "Cst"):
                if not isinstance(value, str):
                    raise ValueError(f"{key}.{field} requires text or null")
            elif field in ("HIP", "HD", "Fl"):
                if type(value) is not int or not 0 <= value <= 0x7FFFFFFF:
                    raise ValueError(f"{key}.{field} requires a nonnegative integer or null")
            elif type(value) not in (int, float) or not math.isfinite(value):
                raise ValueError(f"{key}.{field} requires a finite number or null")
        if key == "hip2":
            ra, dec = row["RArad"], row["DErad"]
            if (ra is not None and not 0 <= ra < 360) or (dec is not None and not -90 <= dec <= 90):
                raise ValueError("normalized Hipparcos coordinates require degrees in their physical ranges")
    return rows


def bayer_key(bayer, flamsteed, constellation):
    """Use the catalog's Greek/component spelling, or its Flamsteed number."""
    constellation = (constellation or "").strip()
    if not constellation:
        return None
    bayer = (bayer or "").strip()
    if bayer:
        head = bayer.rstrip("0123456789")
        component = bayer[len(head):]
        abbreviation = GREEK.get(head.rstrip("."))
        if abbreviation is None:
            return None
        return abbreviation + (component.zfill(2) if component else "") + constellation
    return None if flamsteed is None else f"{flamsteed}{constellation}"


def _parse_wgsn(raw):
    """Read full 18-column proper names and explicit HIP identifications."""
    names = {}
    pattern = re.compile(r"\s(?P<hip>\d{1,6}|_)\s+(?P<hd>\d{1,6}|_)\s+"
                         r"[\d.]+\s+[-+]?[\d.]+\s")
    for line in raw.decode("utf-8").splitlines():
        if not line.strip() or line.startswith(("#", "$")) or len(line.split()) < 8:
            continue
        match, name = pattern.search(line), line[:18].strip()
        if match and name and match["hip"] != "_":
            names[int(match["hip"])] = name
    if not names:
        raise ValueError("WGSN input contains no recognized named HIP stars")
    return names


def read_wgsn(path):
    raw = Path(path).read_bytes()
    return _parse_wgsn(raw), hashlib.sha256(raw).hexdigest()


def _fetch_table(key):
    """Acquire one fixed public query; astronomy is not evaluated here."""
    try:
        from astroquery.vizier import Vizier
    except ImportError as error:
        raise ValueError("--fetch requires optional astroquery (and its Astropy dependency)") from error
    query = QUERIES[key]
    client = Vizier(columns=query["columns"], row_limit=query["row_limit"])
    client.TIMEOUT = 180
    result = client.query_constraints(catalog=query["catalog"], cache=False, **query["constraints"])
    if len(result) != 1:
        raise ValueError(f"expected exactly one public {key} query table")
    table, rows = result[0], []
    source_units = {}
    for field in query["units"]:
        unit = table[field].unit
        if unit is None:
            raise ValueError(f"public {key}.{field} did not declare its required unit")
        source_units[field] = str(unit)
    for row in table:
        converted = {}
        for field in query["columns"]:
            value = row[field]
            if value is None or bool(getattr(value, "mask", False)):
                converted[field] = None
            elif field in ("Bayer", "Cst"):
                converted[field] = str(value)
            elif field in ("HIP", "HD", "Fl"):
                converted[field] = int(value)
            else:
                # I/311 raw coordinates are radians. Trust the declared column
                # unit, never the misleading RArad/DErad field names alone.
                value = (float(value) * table[field].unit).to_value(query["units"][field])
                converted[field] = float(value) if math.isfinite(value) else None
        rows.append(converted)
    normalized = {"schema_version": 1, **copy.deepcopy(query), "source_units": source_units, "rows": rows}
    _validate_table(normalized, key)
    return normalized


def _publish(path, payload, *, overwrite=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(prefix=f".{path.name}-", dir=path.parent)
    temporary = Path(name)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(payload)
            output.flush()
            os.fsync(output.fileno())
        if overwrite:
            os.replace(temporary, path)
        else:
            os.link(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _json_bytes(value):
    return (json.dumps(value, ensure_ascii=False, allow_nan=False, indent=2) + "\n").encode()


def acquire_inputs(directory, *, fetch=False):
    """Reuse every existing input; explicit fetch only fills missing cache files."""
    directory = Path(directory)
    tables, sources = {}, []
    for key, query in QUERIES.items():
        path = directory / f"{key}.json"
        if not path.exists() and fetch:
            _publish(path, _json_bytes(_fetch_table(key)))
        table, digest = _read_json(path)
        tables[key] = _validate_table(table, key)
        sources.append({"key": key, **copy.deepcopy(query), "source_units": table["source_units"],
                        "input_path": str(path.resolve()), "sha256": digest})
    path = directory / "IAU-CSN.txt"
    if not path.exists() and fetch:
        import urllib.request
        request = urllib.request.Request(IAU_CSN_URL, headers={"User-Agent": "libephemeris"})
        with urllib.request.urlopen(request, timeout=120) as response:
            payload = response.read()
        # Validate these exact bytes before they become a reusable cache input.
        _parse_wgsn(payload)
        _publish(path, payload)
    names, digest = read_wgsn(path)
    sources.append({"key": "wgsn", "url": IAU_CSN_URL, "format": WGSN_FORMAT,
                    "input_path": str(path.resolve()), "sha256": digest})
    return tables, names, sources


def build_rows(tables, names, seed_rows):
    """Join published fields and apply native propagation plus stable-ID policy."""
    from libephemeris.fixed_stars import propagate_proper_motion

    curated = {row[3]: row for row in seed_rows if row[3]}
    if len(curated) != len(seed_rows):
        raise ValueError("seed rows require unique nonzero HIP IDs for lossless maintenance")
    preserved = {hip: row[0] for hip, row in curated.items()}
    hd_by_hip, magnitudes = {}, {}
    for row in tables["hipmain"]:
        if row["HIP"] is None:
            continue
        if row["HD"] is not None:
            hd_by_hip[row["HIP"]] = row["HD"]
        if row["Vmag"] is not None:
            magnitudes[row["HIP"]] = row["Vmag"]
    nomen_by_hd = {}
    for row in tables["cross"]:
        nomen = bayer_key(row["Bayer"], row["Fl"], row["Cst"])
        if row["HD"] is not None and nomen:
            nomen_by_hd[row["HD"]] = nomen
    velocities = {row["HIP"]: row["RV"] for row in tables["xhip"]
                  if row["HIP"] is not None and row["RV"] is not None}
    selected = {}
    # Keep the final physical row for duplicate public HIP records, as the
    # pinned maintainer did before selection; cache order is part of the input.
    hip2 = {row["HIP"]: row for row in tables["hip2"]
            if row["HIP"] is not None and row["RArad"] is not None and row["DErad"] is not None}
    for hip, row in hip2.items():
        nomen = nomen_by_hd.get(hd_by_hip.get(hip))
        hp = row["Hpmag"] if row["Hpmag"] is not None else 99.0
        if hip not in names and hip not in preserved and not (hp <= HP_LIMIT and nomen):
            continue
        pmra, pmdec = row["pmRA"] or 0.0, row["pmDE"] or 0.0
        ra, dec = propagate_proper_motion(row["RArad"], row["DErad"], pmra / 1000.0,
                                          pmdec / 1000.0, FROM_JD, TO_JD)
        name = names.get(hip) or nomen or f"HIP{hip}"
        magnitude = magnitudes.get(hip)
        selected[hip] = [0, name, nomen or name, hip, ra, dec, pmra / 1000.0, pmdec / 1000.0,
                         row["Plx"] or 0.0, velocities.get(hip, 0.0), hp if magnitude is None else magnitude]
    for hip, row in curated.items():
        if hip not in selected:
            selected[hip] = list(row)
        else:
            selected[hip][1:3] = row[1:3]
            # Preserve the reviewed seed magnitude, including any explicit
            # photometric correction; its sparse source annotation is copied
            # with the seed below, not inferred for other downloaded rows.
            selected[hip][10] = row[10]
            if selected[hip][9] == 0.0 and row[9] != 0.0:
                selected[hip][9] = row[9]
    by_name = {}
    for hip, row in selected.items():
        by_name.setdefault(row[1], []).append(hip)
    for name, group in by_name.items():
        if len(group) < 2:
            continue
        if len(group) > 16:
            raise ValueError(f"more than 16 same-name components need explicit review: {name}")
        for letter, hip in zip("ABCDEFGHIJKLMNOP", sorted(group)):
            if hip in preserved:
                continue
            row = selected[hip]
            if row[2] == row[1]:
                row[2] = name + letter
            row[1] = name + letter
    next_id = max(preserved.values(), default=1000000) + 1
    rows = []
    for hip, row in sorted(selected.items()):
        if hip in preserved:
            row[0] = preserved[hip]
        else:
            row[0], next_id = next_id, next_id + 1
        row[4:] = [round(float(value), digits) for value, digits in zip(row[4:], ROUND_DIGITS, strict=True)]
        if (not all(math.isfinite(value) for value in row[4:]) or not 0 <= row[4] < 360 or
                not -90 <= row[5] <= 90):
            raise ValueError(f"invalid propagated candidate row: HIP {hip}")
        rows.append(row)
    if len({row[1].upper() for row in rows}) != len(rows):
        raise ValueError("candidate traditional names collide under uppercase lookup")
    report = {"total": len(rows), "iau_named": sum(hip in names for hip in selected),
              "preserved": len(preserved), "with_rv": sum(row[9] != 0.0 for row in rows)}
    return rows, report


def _check_target(output, inputs, overwrite):
    output = Path(output)
    protected = [*inputs, serializer.DEFAULT_INPUT,
                 serializer.REPO_ROOT / "core/src/tables/star_catalog_tables.rs"]
    for path in protected:
        if output.resolve() == path.resolve() or (output.exists() and path.exists() and output.samefile(path)):
            raise ValueError("candidate output must differ from every input and shipped catalog artifact")
    if os.path.lexists(output) and not overwrite:
        raise FileExistsError(output)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", required=True, type=Path)
    parser.add_argument("--seed", type=Path, default=serializer.DEFAULT_INPUT)
    parser.add_argument("--out", required=True, type=Path, help="review-only candidate JSON")
    parser.add_argument("--fetch", action="store_true", help="fetch missing cache inputs; never replace existing files")
    parser.add_argument("--overwrite", action="store_true", help="replace an existing review candidate only")
    args = parser.parse_args(argv)
    try:
        inputs = [args.seed, serializer.DEFAULT_INPUT,
                  serializer.REPO_ROOT / "core/src/tables/star_catalog_tables.rs",
                  *(args.cache_dir / f"{key}.json" for key in QUERIES), args.cache_dir / "IAU-CSN.txt"]
        _check_target(args.out, inputs, args.overwrite)
        seed, seed_hash = serializer.load_catalog(args.seed)
        tables, names, sources = acquire_inputs(args.cache_dir, fetch=args.fetch)
        import libephemeris._core as core
        extension = Path(core.__file__)
        inputs.append(extension)
        _check_target(args.out, inputs, args.overwrite)
        extension_hash = _digest(extension)
        rows, report = build_rows(tables, names, seed["rows"])
        candidate = copy.deepcopy(seed)
        candidate["schema_version"], candidate["rows"] = 2, rows
        candidate["build"] = {
            "builder": "tools/generators/star_catalog_builder.py", "seed_sha256": seed_hash,
            "sources": sources,
            "propagation": {"function": "libephemeris.fixed_stars.propagate_proper_motion",
                            "from_jd": FROM_JD, "to_jd": TO_JD, "extension_sha256": extension_hash},
            "selection": {"hp_limit": HP_LIMIT, "round_digits": ROUND_DIGITS}, "report": report,
        }
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=f".{args.out.name}-review-", dir=args.out.parent) as name:
            staged = Path(name) / "candidate.json"
            _publish(staged, _json_bytes(candidate))
            serializer.load_catalog(staged)
            identities = [(args.seed, seed_hash), (extension, extension_hash)]
            identities.extend((Path(row["input_path"]), row["sha256"]) for row in sources)
            for path, expected in identities:
                if _digest(path) != expected:
                    raise ValueError(f"catalog input changed during build: {path}")
            _check_target(args.out, inputs, args.overwrite)
            if args.overwrite:
                os.replace(staged, args.out)
            else:
                os.link(staged, args.out)
        print(json.dumps(report, sort_keys=True))
        print(f"review candidate: {args.out}; shipped data unchanged")
        return 0
    except (OSError, UnicodeError, ValueError, ImportError) as error:
        print(f"star catalog build failed: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
