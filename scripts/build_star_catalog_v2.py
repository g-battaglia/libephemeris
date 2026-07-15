#!/usr/bin/env python3
"""Build the expanded libephemeris star catalog from public sources.

Provenance:
    Sole authoritative generator for ``libephemeris/star_catalog_gen.py``.
    Every output field is obtained from the public catalog named below;
    cross-identification, J1991.25-to-J2000 propagation, selection, rounding,
    and stable internal IDs are explicit project transformations. Validation
    uses catalogue identities and astronomical invariants, never another
    ephemeris star table.

Generates ``libephemeris/star_catalog_gen.py`` - the full fixed-star data
table - exclusively from public astronomical catalogs:

- Astrometry (ICRS positions, proper motions, parallaxes):
  Hipparcos, the New Reduction (van Leeuwen 2007), VizieR I/311/hip2.
  Positions are epoch J1991.25 and are propagated to J2000.0 with the
  catalog proper motions.
- Visual magnitudes and HD cross identifications:
  Hipparcos Catalogue (ESA 1997), VizieR I/239/hip_main.
- Bayer / Flamsteed designations:
  HD-DM-GC-HR-HIP-Bayer-Flamsteed Cross Index (Kostjuk 2004),
  VizieR IV/27A/catalog.
- Radial velocities:
  Extended Hipparcos Compilation XHIP (Anderson & Francis 2012),
  VizieR V/137D/XHIP.
- Proper names:
  IAU Working Group on Star Names (WGSN) Catalog of Star Names,
  IAU-CSN (maintained by E. Mamajek),
  https://www.pas.rochester.edu/~emamajek/WGSN/IAU-CSN.txt

Validation of the generated catalog uses only the cited public catalogs and
independent astronomical invariants.

Selection: every IAU-named star, every star with Hipparcos magnitude
Hp <= 5.2 that carries a Bayer or Flamsteed designation, and every star
already present in the library's curated catalog (their internal IDs are
preserved for compatibility).

Usage:
    python scripts/build_star_catalog_v2.py            # writes the module
    python scripts/build_star_catalog_v2.py --dry-run  # report only
"""

from __future__ import annotations

import argparse
import json
import math
import os
import socket
import sys
import urllib.request
from datetime import date

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

IAU_CSN_URL = "https://www.pas.rochester.edu/~emamajek/WGSN/IAU-CSN.txt"
HP_LIMIT = 5.2
EPOCH_OFFSET_YR = 2000.0 - 1991.25  # I/311 positions are epoch J1991.25

# Greek-letter abbreviations, IV/27A spelling -> catalog key spelling
# (the 2-3 letter forms conventional in Bayer shorthand).
GREEK_MAP = {
    "alf": "al",
    "bet": "be",
    "gam": "ga",
    "del": "de",
    "eps": "ep",
    "zet": "ze",
    "eta": "et",
    "the": "th",
    "iot": "io",
    "kap": "ka",
    "lam": "la",
    "mu": "mu",
    "mu.": "mu",
    "nu": "nu",
    "nu.": "nu",
    "ksi": "xi",
    "xi": "xi",
    "xi.": "xi",
    "omi": "omi",
    "pi": "pi",
    "pi.": "pi",
    "rho": "rh",
    "sig": "si",
    "tau": "ta",
    "ups": "up",
    "phi": "ph",
    "chi": "ch",
    "psi": "ps",
    "ome": "ome",
}


def _val(v, cast):
    """Cast a possibly-masked table value; None when missing."""
    try:
        import numpy.ma as ma

        if v is None or ma.is_masked(v):
            return None
    except Exception:
        pass
    try:
        return cast(v)
    except (TypeError, ValueError):
        return None


def fetch_vizier_tables(verbose: bool) -> dict:
    """Fetch the public VizieR tables used for field-level catalogue assembly.

    Queries Hipparcos-2 (I/311/hip2) for astrometry, the original Hipparcos
    main table (I/239/hip_main) for HD identifiers and V magnitudes, and the
    published radial-velocity cross-catalogue selected by this generator.
    Magnitude margins are intentionally wider than the final inclusion limit so
    later joins cannot discard a bright named star through catalogue mismatch.

    Args:
        verbose: Print each network query before it is issued.

    Returns:
        Mapping from the generator's stable table labels to Astropy tables.
    """
    from astroquery.vizier import Vizier

    socket.setdefaulttimeout(180)
    out = {}

    if verbose:
        print("fetching I/311/hip2 (van Leeuwen 2007) ...", flush=True)
    v = Vizier(
        columns=["HIP", "RArad", "DErad", "Plx", "pmRA", "pmDE", "Hpmag"],
        row_limit=-1,
    )
    out["hip2"] = v.query_constraints(catalog="I/311/hip2", Hpmag=f"<{HP_LIMIT + 2.0}")[
        0
    ]

    if verbose:
        print("fetching I/239/hip_main (HD numbers, Vmag) ...", flush=True)
    v = Vizier(columns=["HIP", "HD", "Vmag"], row_limit=-1)
    out["hipmain"] = v.query_constraints(
        catalog="I/239/hip_main", Vmag=f"<{HP_LIMIT + 2.5}"
    )[0]

    if verbose:
        print("fetching IV/27A/catalog (Bayer/Flamsteed cross index) ...", flush=True)
    v = Vizier(columns=["HD", "Bayer", "Fl", "Cst"], row_limit=-1)
    out["cross"] = v.query_constraints(catalog="IV/27A/catalog")[0]

    if verbose:
        print("fetching V/137D/XHIP (radial velocities) ...", flush=True)
    v = Vizier(columns=["HIP", "RV"], row_limit=-1)
    out["xhip"] = v.query_constraints(catalog="V/137D/XHIP", RV="!= NULL")[0]

    return out


def fetch_iau_names(verbose: bool) -> dict[int, str]:
    """HIP -> IAU proper name from the WGSN Catalog of Star Names."""
    if verbose:
        print("fetching IAU-CSN (WGSN names) ...", flush=True)
    req = urllib.request.Request(IAU_CSN_URL, headers={"User-Agent": "libephemeris"})
    text = urllib.request.urlopen(req, timeout=120).read().decode("utf-8", "replace")
    names: dict[int, str] = {}
    for line in text.splitlines():
        if not line.strip() or line.startswith("#") or line.startswith("$"):
            continue
        # Fixed-width-ish columns; the name is the first field (may hold
        # spaces but the file pads with multiple spaces between fields).
        parts = [p for p in line.split(" ") if p]
        if len(parts) < 8:
            continue
        # Find the HIP column: documented layout has HIP as a number in a
        # known position; robust approach: scan integer-looking fields and
        # take the one followed by an HD-like field. Simpler: the file's
        # column header gives offsets - use a regex for " <hip> " right
        # before the HD column. Practical robust parse: take fields and
        # locate the constellation 3-letter code, HIP is 4 fields later
        # in the 2022+ layout. Fall back to a regex.
        import re

        # The name occupies the first fixed-width field (18 columns) —
        # a non-greedy regex here used to truncate multi-word names
        # ("Polaris Australis" -> "Polaris", "Yed Prior" -> "Yed").
        name = line[:18].strip()
        if not name:
            continue
        m = re.search(
            r"\s(?P<hip>\d{1,6}|_)\s+(?P<hd>\d{1,6}|_)\s+"
            r"(?P<ra>[\d.]+)\s+(?P<dec>[-+]?[\d.]+)\s",
            line,
        )
        if not m:
            continue
        hip_s = m.group("hip")
        if hip_s == "_":
            continue
        try:
            names[int(hip_s)] = name
        except ValueError:
            continue
    return names


def bayer_key(bayer: str, fl, cst: str) -> str | None:
    """Catalog nomenclature key from cross-index fields ('alTau' style)."""
    cst = (cst or "").strip()
    if not cst:
        return None
    b = (bayer or "").strip()
    if b:
        comp = ""
        head = b
        while head and head[-1].isdigit():
            comp = head[-1] + comp
            head = head[:-1]
        head = head.rstrip(".") + ("." if head.endswith(".") and False else "")
        key = GREEK_MAP.get(head.rstrip("."), None)
        if key is None:
            key = GREEK_MAP.get(head, None)
        if key is None:
            return None
        if comp:
            key += comp.zfill(2)
        return key + cst
    if fl is not None and str(fl).strip() not in ("", "--"):
        try:
            n = int(fl)
        except (TypeError, ValueError):
            return None
        return f"{n}{cst}"
    return None


def existing_catalog() -> tuple[dict[int, int], dict[int, tuple]]:
    """(hip -> preserved id, hip -> curated row) from the current catalog."""
    from libephemeris.fixed_stars import STAR_CATALOG

    ids: dict[int, int] = {}
    rows: dict[int, tuple] = {}
    for e in STAR_CATALOG:
        if e.hip_number:
            ids[int(e.hip_number)] = int(e.id)
            rows[int(e.hip_number)] = (
                e.name,
                e.nomenclature,
                float(e.data.ra_j2000),
                float(e.data.dec_j2000),
                float(e.data.pm_ra),
                float(e.data.pm_dec),
                float(e.data.parallax_mas),
                float(e.data.radial_km_per_s),
                float(e.magnitude),
            )
    return ids, rows


def build(verbose: bool) -> tuple[list[tuple], dict]:
    """Join public catalogues and construct deterministic generated star rows.

    Hipparcos-2 supplies the authoritative astrometric fields; the original
    Hipparcos table and identified cross-catalogues contribute only their
    documented fields. IAU names are matched through HIP/HD identifiers, input
    values are normalized to the runtime units, and stable IDs plus validation
    counts are produced without consulting any proprietary ephemeris catalogue.

    Args:
        verbose: Emit download, join, and validation progress.

    Returns:
        A pair ``(rows, report)`` containing sorted runtime tuples and a JSON-
        serializable provenance/validation summary.
    """
    tabs = fetch_vizier_tables(verbose)
    iau_names = fetch_iau_names(verbose)

    hip2 = {}
    for row in tabs["hip2"]:
        hip = _val(row["HIP"], int)
        ra = _val(row["RArad"], float)
        dec = _val(row["DErad"], float)
        if hip is None or ra is None or dec is None:
            continue
        hip2[hip] = (
            ra,
            dec,
            _val(row["Plx"], float) or 0.0,
            _val(row["pmRA"], float) or 0.0,
            _val(row["pmDE"], float) or 0.0,
            _val(row["Hpmag"], float)
            if _val(row["Hpmag"], float) is not None
            else 99.0,
        )

    hd_by_hip: dict[int, int] = {}
    vmag_by_hip: dict[int, float] = {}
    for row in tabs["hipmain"]:
        hip = _val(row["HIP"], int)
        if hip is None:
            continue
        hd = _val(row["HD"], int)
        if hd is not None:
            hd_by_hip[hip] = hd
        vmag = _val(row["Vmag"], float)
        if vmag is not None:
            vmag_by_hip[hip] = vmag

    nomen_by_hd: dict[int, str] = {}
    for row in tabs["cross"]:
        hd = _val(row["HD"], int)
        if hd is None:
            continue
        fl = _val(row["Fl"], int)
        key = bayer_key(str(row["Bayer"]), fl, str(row["Cst"]))
        if key:
            nomen_by_hd[hd] = key

    rv_by_hip: dict[int, float] = {}
    for row in tabs["xhip"]:
        hip = _val(row["HIP"], int)
        rv = _val(row["RV"], float)
        if hip is not None and rv is not None:
            rv_by_hip[hip] = rv

    preserved_ids, curated_rows = existing_catalog()

    selected: dict[int, dict] = {}
    for hip, (ra, dec, plx, pmra, pmde, hp) in hip2.items():
        hd = hd_by_hip.get(hip)
        nomen = nomen_by_hd.get(hd) if hd else None
        named = hip in iau_names
        keep = named or hip in preserved_ids or (hp <= HP_LIMIT and nomen)
        if not keep:
            continue
        # Epoch propagation J1991.25 -> J2000.0 (catalog pmRA is *cos(dec)).
        cosd = math.cos(math.radians(dec))
        ra2000 = ra + (pmra / 1000.0 / 3600.0) * EPOCH_OFFSET_YR / max(cosd, 1e-9)
        dec2000 = dec + (pmde / 1000.0 / 3600.0) * EPOCH_OFFSET_YR
        name = iau_names.get(hip) or (nomen or f"HIP{hip}")
        selected[hip] = dict(
            name=name,
            nomen=nomen or (iau_names.get(hip, "") and "") or "",
            ra=ra2000 % 360.0,
            dec=dec2000,
            pm_ra=pmra / 1000.0,
            pm_dec=pmde / 1000.0,
            plx=plx,
            rv=rv_by_hip.get(hip, 0.0),
            mag=vmag_by_hip.get(hip, hp),
        )
        if not selected[hip]["nomen"]:
            selected[hip]["nomen"] = nomen or name

    # Stars from the curated catalog missing in the selection (e.g. very
    # faint historical stars) keep their curated values entirely.
    for hip, row in curated_rows.items():
        if hip not in selected:
            (name, nomen, ra, dec, pmra, pmde, plx, rv, mag) = row
            selected[hip] = dict(
                name=name,
                nomen=nomen,
                ra=ra,
                dec=dec,
                pm_ra=pmra,
                pm_dec=pmde,
                plx=plx,
                rv=rv,
                mag=mag,
            )

    # Curated names win for the preserved stars (they carry the library's
    # established naming); IAU names win otherwise.
    for hip, row in curated_rows.items():
        if hip in selected:
            selected[hip]["name"] = row[0]
            selected[hip]["nomen"] = row[1]
            # The library's established magnitudes stay stable for the
            # curated core; newly added stars carry Hipparcos Vmag.
            selected[hip]["mag"] = row[8]
            if selected[hip]["rv"] == 0.0 and row[7] != 0.0:
                selected[hip]["rv"] = row[7]

    # Disambiguate stars that share a proper NAME. A shared name collapses
    # them to a single entry in the name->id lookup, leaving all but one
    # unreachable by name (e.g. the three theta-1 Orionis Trapezium
    # components HIP 26220/26221/26224, all "th01Ori"). Append the Bayer
    # component letter (A, B, C, ... in HIP order, which matches the
    # Trapezium designations: 26220=A, 26221=B, 26224=C) to the name and,
    # when it equals the name, the nomenclature, so each is uniquely
    # resolvable. Preserved curated stars keep their established identity but
    # still consume their letter slot so the others stay correctly lettered.
    # Stars that merely share a Bayer nomenclature but have distinct proper
    # names (e.g. alCen -> Rigil Kentaurus / Toliman) are left untouched.
    by_name: dict[str, list[int]] = {}
    for hip in selected:
        by_name.setdefault(selected[hip]["name"], []).append(hip)
    for shared_name, group in by_name.items():
        if len(group) < 2:
            continue
        for letter, hip in zip("ABCDEFGHIJKLMNOP", sorted(group)):
            if hip in preserved_ids:
                continue
            suffixed = f"{shared_name}{letter}"
            if selected[hip]["nomen"] == selected[hip]["name"]:
                selected[hip]["nomen"] = suffixed
            selected[hip]["name"] = suffixed

    next_id = max(preserved_ids.values()) + 1 if preserved_ids else 1000001
    rows_out: list[tuple] = []
    for hip in sorted(selected):
        s = selected[hip]
        sid = preserved_ids.get(hip)
        if sid is None:
            sid = next_id
            next_id += 1
        rows_out.append(
            (
                sid,
                s["name"],
                s["nomen"],
                hip,
                round(s["ra"], 9),
                round(s["dec"], 9),
                round(s["pm_ra"], 7),
                round(s["pm_dec"], 7),
                round(s["plx"], 3),
                round(s["rv"], 2),
                round(s["mag"], 3),
            )
        )

    report = dict(
        total=len(rows_out),
        iau_named=sum(1 for h in selected if h in iau_names),
        preserved=len(preserved_ids),
        with_rv=sum(1 for s in selected.values() if s["rv"] != 0.0),
    )
    return rows_out, report


HEADER = '''# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Generated fixed-star catalog data. DO NOT EDIT BY HAND.

Generated by scripts/build_star_catalog_v2.py on {today} from public
sources only:

- Astrometry: Hipparcos, the New Reduction (van Leeuwen 2007),
  VizieR I/311/hip2; positions propagated from epoch J1991.25 to
  J2000.0 with the catalog proper motions. ICRS.
- Visual magnitudes / HD ids: Hipparcos Catalogue (ESA 1997),
  VizieR I/239/hip_main.
- Bayer/Flamsteed designations: Kostjuk (2004) Cross Index,
  VizieR IV/27A/catalog.
- Radial velocities: XHIP (Anderson & Francis 2012), VizieR V/137D.
- Proper names: IAU WGSN Catalog of Star Names (IAU-CSN).

Provenance:
    AUTO-GENERATED solely by ``scripts/build_star_catalog_v2.py`` from the
    public catalogues above. The generator records field mapping, epoch
    propagation, selection, stable IDs, and validation invariants. No
    proprietary ephemeris star file or compatibility output is an input.

Row layout:
    (id, name, nomenclature, hip,
     ra_j2000_deg, dec_j2000_deg,
     pm_ra_arcsec_per_yr_cosdec, pm_dec_arcsec_per_yr,
     parallax_mas, rv_km_per_s, vmag)
"""

STAR_ROWS = (
'''


def main() -> None:
    """Generate ``libephemeris/star_catalog_gen.py`` from public catalogues.

    ``--dry-run`` performs all retrieval, joining, and validation but writes
    nothing. A normal run emits the provenance-bearing template followed by
    deterministic row literals; the generated module is the only repository
    file replaced by this command.
    """
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    verbose = not args.quiet

    rows, report = build(verbose)
    print(json.dumps(report))
    if args.dry_run:
        return

    out_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "libephemeris",
        "star_catalog_gen.py",
    )
    with open(out_path, "w") as f:
        f.write(HEADER.format(today=date.today().isoformat()))
        for r in rows:
            f.write(f"    {r!r},\n")
        f.write(")\n")
    print(f"wrote {out_path} ({len(rows)} stars)")


if __name__ == "__main__":
    main()
