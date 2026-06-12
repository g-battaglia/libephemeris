"""Verification harness for the houses.py clean-room provenance work (WS1).

Compares libephemeris ``houses_armc`` / ``house_pos`` against pyswisseph over a
dense (system, ARMC, latitude) grid, and supports dump/check invariance gating
so that expression-level rewrites can be proven output-preserving.

Modes:
    compare            Side-by-side vs pyswisseph (informational summary).
    --dump FILE        Record every libephemeris output (incl. raises) as JSON.
    --check FILE       Re-run the grid and fail (exit 1) if any output drifted
                       from the dump by more than --tol degrees (default 1e-9).

Usage:
    .venv/bin/python scripts/verify_houses_clean_room.py compare
    .venv/bin/python scripts/verify_houses_clean_room.py --dump tasks/results/houses_baseline.json
    .venv/bin/python scripts/verify_houses_clean_room.py --check tasks/results/houses_baseline.json
    ... --systems J,U,Y,i,P   (restrict to specific house systems)

This script is standalone (project convention: verification via standalone
scripts, not the pytest suite). It needs no ephemeris files: houses_armc and
house_pos are pure functions of (armc, lat, eps[, ascmc9]).
"""

from __future__ import annotations

import argparse
import json
import math
import sys

# House systems accepted by pyswisseph 2.10 houses_armc/house_pos (probed).
ALL_SYSTEMS = list("ABCDEFGHIJKLMNOPQRSTUVWXY") + ["i"]
# 'Z' (unknown code) deliberately excluded from the default grid: SE treats it
# as Placidus and raises at polar latitudes while libephemeris falls back to
# Porphyry — a known, separately-tracked parity issue (WS10), not WS1 scope.

ARMC_GRID = [
    0.0, 30.5, 45.0, 90.0, 123.456, 135.0,
    180.0, 222.2, 270.0, 300.001, 315.0, 359.9999,
]
LAT_GRID = [
    0.0, 0.001, -0.001, 5.5, -5.5, 23.44, -23.44, 41.9, -41.9,
    51.5, -51.5, 60.0, -60.0, 66.0, -66.0, 66.55, -66.55, 67.0, -67.0,
    70.0, -70.0, 80.0, -80.0, 85.0, -85.0, 88.0, -88.0, 89.0, -89.0,
    89.99, -89.99,
]
EPS = 23.4392911
# Sunshine systems take the Sun's declination via ascmc9.
SUNSHINE_ASCMC9 = [0.0, 11.47, -23.44]

# house_pos sub-grid (heavier per-case cost, so smaller).
HP_ARMC = [0.0, 90.0, 123.456, 270.0]
HP_GEOLAT = [0.0, 23.44, -23.44, 51.5, -51.5, 66.55, -66.55, 80.0, -80.0, 89.0, -89.0]
HP_LON = [float(x) for x in range(0, 360, 30)] + [15.3, 123.456]
HP_LATB = [0.0, 1.2, -5.0, 17.0, -17.0]
# Systems whose house_pos is verified bit-perfect vs SE (review): gate tightly.
HP_GATED = set("PGKRCXMBOLQ")


def _koch_on_validity_edge(case) -> bool:
    """True when the Koch arc fraction sits exactly on 0 or 2."""
    import math as _m

    armc, geolat, blon, blat = case[2], case[3], case[4], case[5]
    er = _m.radians(EPS)
    br = _m.radians(blat)
    lr = _m.radians(blon)
    dec = _m.degrees(
        _m.asin(_m.sin(br) * _m.cos(er) + _m.cos(br) * _m.sin(lr) * _m.sin(er))
    )
    ra = _m.degrees(
        _m.atan2(_m.sin(lr) * _m.cos(er) - _m.tan(br) * _m.sin(er), _m.cos(lr))
    ) % 360.0
    md = (ra - armc) % 360.0
    if md >= 180.0:
        md -= 360.0
    prod = _m.tan(_m.radians(geolat)) * _m.tan(_m.radians(dec))
    ad = _m.degrees(_m.asin(max(-1.0, min(1.0, prod))))
    ad_mc = _m.degrees(
        _m.asin(
            max(
                -1.0,
                min(
                    1.0,
                    _m.tan(er)
                    * _m.tan(_m.radians(geolat))
                    * _m.sin(_m.radians(armc)),
                ),
            )
        )
    )
    sa = 90.0 + ad_mc
    if sa == 0:
        return True
    if md >= 0:
        frac = (md - ad + ad_mc) / sa
    else:
        frac = (md + 180.0 + ad + ad_mc) / sa
    return min(abs(frac - 0.0), abs(frac - 2.0)) < 1e-9


def _body_on_meridian(case) -> bool:
    """True when the case's body sits on the meridian (md ~ 0 or 180)."""
    import math as _m

    armc, blon = case[2], case[4]
    y = _m.sin(_m.radians(blon)) * _m.cos(_m.radians(EPS))
    x = _m.cos(_m.radians(blon))
    ra = _m.degrees(_m.atan2(y, x)) % 360.0
    md = abs((ra - armc + 540.0) % 360.0 - 180.0)
    return md < 1e-9 or abs(md - 180.0) < 1e-9


def wrap_delta(a: float, b: float) -> float:
    """Smallest absolute angular difference in degrees."""
    d = math.fmod(abs(a - b), 360.0)
    return min(d, 360.0 - d)


def house_pos_delta(a: float, b: float, hsys: str) -> float:
    """House-position difference, wrap-aware at the N->1 boundary (N=36 for G)."""
    n = 36.0 if hsys == "G" else 12.0
    d = math.fmod(abs(a - b), n)
    return min(d, n - d)


def case_key(kind: str, hsys: str, *params: float) -> str:
    return f"{kind}|{hsys}|" + "|".join(repr(p) for p in params)


def run_le_houses(hsys: str, armc: float, lat: float, ascmc9: float):
    """Run libephemeris houses_armc; return ('ok', cusps+ascmc) or ('raise', name)."""
    import libephemeris as le

    try:
        cusps, ascmc = le.houses_armc(armc, lat, EPS, ord(hsys), ascmc9)
        return ("ok", list(cusps) + list(ascmc))
    except Exception as exc:  # noqa: BLE001 - recorded, not hidden
        return ("raise", type(exc).__name__)


def run_se_houses(hsys: str, armc: float, lat: float, ascmc9: float):
    import swisseph as swe

    try:
        cusps, ascmc = swe.houses_armc(armc, lat, EPS, hsys.encode(), ascmc9)
        return ("ok", list(cusps) + list(ascmc))
    except Exception as exc:  # noqa: BLE001
        return ("raise", type(exc).__name__)


def run_le_house_pos(hsys: str, armc: float, geolat: float, lon: float, latb: float):
    import libephemeris as le

    try:
        return ("ok", le.house_pos(armc, geolat, EPS, (lon, latb), hsys))
    except Exception as exc:  # noqa: BLE001
        return ("raise", type(exc).__name__)


def run_se_house_pos(hsys: str, armc: float, geolat: float, lon: float, latb: float):
    import swisseph as swe

    try:
        return ("ok", swe.house_pos(armc, geolat, EPS, (lon, latb), hsys.encode()))
    except Exception as exc:  # noqa: BLE001
        return ("raise", type(exc).__name__)


def iter_cases(systems: list[str]):
    """Yield (key, runner_args) for every grid case."""
    for hsys in systems:
        ascmc9_values = SUNSHINE_ASCMC9 if hsys in ("I", "i") else [0.0]
        for armc in ARMC_GRID:
            for lat in LAT_GRID:
                for a9 in ascmc9_values:
                    yield ("houses", hsys, armc, lat, a9)
        for armc in HP_ARMC:
            for geolat in HP_GEOLAT:
                for lon in HP_LON:
                    for latb in HP_LATB:
                        yield ("house_pos", hsys, armc, geolat, lon, latb)


def evaluate_le(case):
    if case[0] == "houses":
        _, hsys, armc, lat, a9 = case
        return run_le_houses(hsys, armc, lat, a9)
    _, hsys, armc, geolat, lon, latb = case
    return run_le_house_pos(hsys, armc, geolat, lon, latb)


def evaluate_se(case):
    if case[0] == "houses":
        _, hsys, armc, lat, a9 = case
        return run_se_houses(hsys, armc, lat, a9)
    _, hsys, armc, geolat, lon, latb = case
    return run_se_house_pos(hsys, armc, geolat, lon, latb)


def do_dump(systems: list[str], path: str) -> int:
    data = {}
    for case in iter_cases(systems):
        kind, hsys = case[0], case[1]
        status, value = evaluate_le(case)
        data[case_key(kind, hsys, *case[2:])] = [status, value]
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh)
    print(f"dumped {len(data)} cases -> {path}")
    return 0


def do_check(systems: list[str], path: str, tol: float) -> int:
    with open(path, encoding="utf-8") as fh:
        baseline = json.load(fh)
    checked = drifted = missing = 0
    worst: tuple[float, str] = (0.0, "")
    for case in iter_cases(systems):
        key = case_key(case[0], case[1], *case[2:])
        if key not in baseline:
            missing += 1
            continue
        checked += 1
        b_status, b_value = baseline[key]
        status, value = evaluate_le(case)
        if status != b_status:
            drifted += 1
            print(f"DRIFT(status) {key}: {b_status}/{b_value!r} -> {status}/{value!r}")
            continue
        if status == "raise":
            if value != b_value:
                drifted += 1
                print(f"DRIFT(exc) {key}: {b_value} -> {value}")
            continue
        if case[0] == "houses":
            if len(value) != len(b_value):
                drifted += 1
                print(f"DRIFT(len) {key}: {len(b_value)} -> {len(value)}")
                continue
            delta = max(wrap_delta(a, b) for a, b in zip(value, b_value))
        else:
            delta = house_pos_delta(value, b_value, case[1])
        if delta > worst[0]:
            worst = (delta, key)
        if delta > tol:
            drifted += 1
            print(f"DRIFT(value) {key}: max |delta| = {delta:.3e}")
    print(
        f"check: {checked} cases, {drifted} drifted, {missing} not in baseline; "
        f"worst delta {worst[0]:.3e} @ {worst[1] or 'n/a'}"
    )
    return 1 if drifted else 0


def do_compare(systems: list[str], cusp_tol: float, hp_tol: float) -> int:
    per_system: dict[str, dict[str, int]] = {}
    worst: dict[str, tuple[float, str]] = {}
    for case in iter_cases(systems):
        kind, hsys = case[0], case[1]
        stats = per_system.setdefault(hsys, {"match": 0, "diff": 0, "raise_ok": 0, "raise_mismatch": 0, "info_diff": 0})
        le_status, le_value = evaluate_le(case)
        se_status, se_value = evaluate_se(case)
        key = case_key(kind, hsys, *case[2:])
        if le_status == "raise" or se_status == "raise":
            if le_status == se_status == "raise":
                stats["raise_ok"] += 1
            else:
                stats["raise_mismatch"] += 1
                print(f"RAISE-MISMATCH {key}: LE={le_status}:{le_value!r} SE={se_status}:{se_value!r}")
            continue
        if kind == "houses":
            n = min(len(le_value), len(se_value))
            delta = max(wrap_delta(a, b) for a, b in zip(le_value[:n], se_value[:n]))
            gated = not (hsys in ("I", "i") and case[4] != 0.0)
            tol = cusp_tol
            if hsys == "Q" and case[3] == 0.0 and (
                min(case[2] % 90.0, 90.0 - case[2] % 90.0) < 0.001
            ):
                # Pullen SR at a 90.000-degree degenerate quadrant: the
                # reference's ratio solver stops ~1.3e-5 deg from the
                # exact solution (same truncation family as Placidus
                # below); ours converges further.
                tol = max(tol, 5e-5)
            if hsys in ("P", "G") and abs(abs(case[3]) - 66.55) < 0.02:
                # Within ~0.01 deg of the polar circle the reference's
                # Placidus/Gauquelin iteration stops with a condition
                # residual of ~1e-6 deg while ours solves the cusp
                # equation to <1e-7 (verified: SE residual 1.06e-6,
                # LE 4e-8 at armc 30.5, lat -66.55).  Allow the
                # reference's own truncation there.
                tol = max(tol, 1e-5)
        else:
            delta = house_pos_delta(le_value, se_value, hsys)
            gated = hsys in HP_GATED
            tol = hp_tol
            if hsys == "K" and _koch_on_validity_edge(case):
                # The reference's Koch validity classification at a
                # fraction of exactly 0 or 2 (saturated ascensional
                # difference or exact meridian alignment) is decided by
                # its internal float noise: geometrically mirrored
                # cases get opposite answers (verified: SE returns 0.0
                # at armc 0, lat 23.44, body 180/-5 but 4.0 at armc 0,
                # lat -51.5, body 90/17 — both land fraction = 2).
                # Informational only.
                gated = False
            if hsys == "R" and _body_on_meridian(case):
                # The reference's Regiomontanus position is discontinuous
                # at exactly md = 0 (verified: 4.0006 -> 10.0000 -> 3.9994
                # across the meridian at armc 90, lat -51.5, blat 17).
                # We return the two-sided closed-form limit; treat the
                # reference's own jump point as informational.
                gated = False
        if delta <= tol:
            stats["match"] += 1
        elif gated:
            stats["diff"] += 1
            if delta > worst.get(hsys, (0.0, ""))[0]:
                worst[hsys] = (delta, key)
        else:
            stats["info_diff"] += 1
    total_diff = 0
    for hsys in sorted(per_system):
        s = per_system[hsys]
        total_diff += s["diff"] + s["raise_mismatch"]
        line = (
            f"  {hsys}: match={s['match']} diff={s['diff']} info_diff={s['info_diff']} "
            f"raise_ok={s['raise_ok']} raise_mismatch={s['raise_mismatch']}"
        )
        if hsys in worst:
            line += f"  worst={worst[hsys][0]:.3e} @ {worst[hsys][1]}"
        print(line)
    print(f"compare: {total_diff} gated mismatches vs pyswisseph")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", nargs="?", default="compare", choices=["compare"])
    ap.add_argument("--dump", metavar="FILE")
    ap.add_argument("--check", metavar="FILE")
    ap.add_argument("--systems", help="comma-separated house codes (default: all)")
    ap.add_argument("--tol", type=float, default=1e-9, help="check tolerance (deg)")
    ap.add_argument("--cusp-tol", type=float, default=1e-6)
    ap.add_argument("--hp-tol", type=float, default=1e-6)
    args = ap.parse_args()

    systems = args.systems.split(",") if args.systems else ALL_SYSTEMS
    for code in systems:
        if code not in ALL_SYSTEMS:
            ap.error(f"unknown house system {code!r}")

    if args.dump:
        return do_dump(systems, args.dump)
    if args.check:
        return do_check(systems, args.check, args.tol)
    return do_compare(systems, args.cusp_tol, args.hp_tol)


if __name__ == "__main__":
    sys.exit(main())
