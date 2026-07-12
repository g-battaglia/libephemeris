#!/usr/bin/env python3
"""Hypothetical-body orbital-element provenance gate (zero-mismatch).

Enforces, for the Hamburg-school Uranian planets and the trans-Plutonian
bodies in ``libephemeris/hypothetical.py``, that the elements actually used
at runtime trace to published sources rather than drifting silently:

  * The live Keplerian elements (``URANIAN_KEPLERIAN_ELEMENTS``) must equal
    ``data/fictitious_orbits.csv`` field-for-field per body — the CSV
    (Witte/Lefeldt 1928, refined by Neely 1988) is the single source of
    truth — and the mean motion must equal the Gaussian value n = k/a^1.5.
  * Proserpina's live elements must equal its CSV row (Abramov).
  * Transpluto's live elements must be reproducible from the published
    Strubell 1952 row by mean-anomaly propagation (epoch -> J2000) and
    general precession of the argument of perihelion (J1945 -> J2000).
  * The legacy display tables (``URANIAN_ELEMENTS`` and the per-body
    ``*_KEPLERIAN_ELEMENTS`` ``L0`` / Hades ``M0``) — historical
    pyswisseph-oracle fits, no longer consulted at runtime but pinned by
    tests for API stability — must match a frozen snapshot, so any future
    edit is a conscious reclassification rather than silent drift.

See docs/methodology/hypothetical-bodies.md for the per-body source table
and the manual CSV<->publication cross-check record.

Usage:
    python scripts/check_hypothetical_provenance.py            # gate
    python scripts/check_hypothetical_provenance.py --explain  # show derivations
"""

from __future__ import annotations

import argparse
import sys

from libephemeris import hypothetical as hyp

# IAU 2006 general precession in longitude at J2000 (Capitaine et al. 2003):
# 5028.796195 arcsec / Julian century.
_PRECESSION_ARCSEC_PER_CENTURY = 5028.796195
_J2000_JD = 2451545.0
_DAYS_PER_CENTURY = 36525.0

# Frozen snapshot of the legacy display tables (oracle-calibrated fits, not
# consulted at runtime). Drift here must be a conscious reclassification.
_LEGACY_URANIAN_L0 = {
    "Cupido": 241.2067,
    "Hades": 176.0581,
    "Zeus": 32.1893,
    "Kronos": 213.2096,
    "Apollon": 71.8925,
    "Admetos": 142.2269,
    "Vulkanus": 195.6753,
    "Poseidon": 274.4073,
}
_LEGACY_PERBODY = {
    # body: (table attr name on the module, field, frozen value)
    "CUPIDO_KEPLERIAN_ELEMENTS": ("L0", 105.301693),
    "HADES_KEPLERIAN_ELEMENTS": ("M0", 26.850162),
    "ZEUS_KEPLERIAN_ELEMENTS": ("L0", 104.289095),
    "KRONOS_KEPLERIAN_ELEMENTS": ("L0", 17.111353),
    "APOLLON_KEPLERIAN_ELEMENTS": ("L0", 138.565328),
    "ADMETOS_KEPLERIAN_ELEMENTS": ("L0", 350.613913),
    "VULKANUS_KEPLERIAN_ELEMENTS": ("L0", 55.397715),
    "POSEIDON_KEPLERIAN_ELEMENTS": ("L0", 166.140256),
}

# CSV body-name -> Uranian body id constant
_URANIAN_CSV_NAME = {
    "Cupido": hyp.CUPIDO,
    "Hades": hyp.HADES,
    "Zeus": hyp.ZEUS,
    "Kronos": hyp.KRONOS,
    "Apollon": hyp.APOLLON,
    "Admetos": hyp.ADMETOS,
    "Vulkanus": hyp.VULKANUS,
    "Poseidon": hyp.POSEIDON,
}


def _csv_by_name() -> dict[str, object]:
    return {e.name: e for e in hyp.load_bundled_fictitious_orbits()}


def check_uranian(rows, problems, explain):
    for csv_name, body_id in _URANIAN_CSV_NAME.items():
        row = rows[csv_name]
        live = hyp.URANIAN_KEPLERIAN_ELEMENTS[body_id]
        pairs = {
            "a": (live.a, row.semi_axis),
            "e": (live.e, row.eccentricity.evaluate(0.0)),
            "i": (live.i, row.inclination.evaluate(0.0)),
            "omega": (live.omega, row.arg_perihelion.evaluate(0.0)),
            "Omega": (live.Omega, row.asc_node.evaluate(0.0)),
            "M0": (live.M0, row.mean_anomaly.evaluate(0.0)),
            "epoch": (live.epoch, row.epoch_jd),
        }
        for field, (got, want) in pairs.items():
            if got != want:
                problems.append(
                    f"Uranian {csv_name}.{field}: live {got!r} != CSV {want!r}"
                )
        expected_n = hyp._K_GAUSS_DEG / live.a**1.5
        if live.n != expected_n:
            problems.append(
                f"Uranian {csv_name}.n: {live.n!r} != k/a^1.5 {expected_n!r}"
            )
        if explain:
            print(f"  {csv_name:9s} a={live.a} e={live.e} M0={live.M0} -> CSV match")


def check_proserpina(rows, problems, explain):
    row = rows["Proserpina"]
    live = hyp.HYPOTHETICAL_ELEMENTS[hyp.PROSERPINA]
    pairs = {
        "a": (live.a, row.semi_axis),
        "e": (live.e, row.eccentricity.evaluate(0.0)),
        "i": (live.i, row.inclination.evaluate(0.0)),
        "M0": (live.M0, row.mean_anomaly.evaluate(0.0)),
        "epoch": (live.epoch, row.epoch_jd),
    }
    for field, (got, want) in pairs.items():
        if got != want:
            problems.append(f"Proserpina.{field}: live {got!r} != CSV {want!r}")
    if explain:
        print(f"  Proserpina a={live.a} M0={live.M0} epoch={live.epoch} -> CSV match")


def check_transpluto(rows, problems, explain):
    # Strubell 1952 originals from the CSV row.
    row = rows["Isis-Transpluto"]
    strubell_epoch = row.epoch_jd  # 2368547.66
    strubell_equinox = row.equinox_jd  # 2431456.5 = J1945.0
    strubell_omega = row.arg_perihelion.evaluate(0.0)  # 0.7
    strubell_M0 = row.mean_anomaly.evaluate(0.0)  # 0.0
    live = hyp.HYPOTHETICAL_ELEMENTS[hyp.ISIS]

    # Mean-anomaly propagation epoch -> J2000 with the live Gaussian motion.
    derived_M0 = (strubell_M0 + live.n * (_J2000_JD - strubell_epoch)) % 360.0
    # Argument of perihelion precesses with the general precession in
    # longitude from the element equinox (J1945) to J2000.
    centuries = (_J2000_JD - strubell_equinox) / _DAYS_PER_CENTURY
    precession_deg = _PRECESSION_ARCSEC_PER_CENTURY * centuries / 3600.0
    derived_omega = strubell_omega + precession_deg

    if abs(derived_M0 - live.M0) > 1e-2:
        problems.append(
            f"Transpluto M0: derived {derived_M0:.6f} vs live {live.M0} "
            f"(|Δ| {abs(derived_M0 - live.M0):.2e} > 1e-2)"
        )
    if abs(derived_omega - live.omega) > 1e-4:
        problems.append(
            f"Transpluto omega: derived {derived_omega:.6f} vs live {live.omega} "
            f"(|Δ| {abs(derived_omega - live.omega):.2e} > 1e-4)"
        )
    if explain:
        print(
            f"  Transpluto: M0 {strubell_M0} +n·{_J2000_JD - strubell_epoch:.2f}d "
            f"= {derived_M0:.6f} (live {live.M0}); "
            f"omega {strubell_omega} +{precession_deg:.6f}° precession "
            f"= {derived_omega:.6f} (live {live.omega})"
        )


def check_legacy(problems, explain):
    for name, l0 in _LEGACY_URANIAN_L0.items():
        body_id = _URANIAN_CSV_NAME[name]
        got = hyp.URANIAN_ELEMENTS[body_id].L0
        if got != l0:
            problems.append(f"legacy URANIAN_ELEMENTS[{name}].L0: {got!r} != {l0!r}")
    for table_name, (field, frozen) in _LEGACY_PERBODY.items():
        got = getattr(getattr(hyp, table_name), field)
        if got != frozen:
            problems.append(f"legacy {table_name}.{field}: {got!r} != {frozen!r}")
    if explain:
        print("  legacy tables match frozen snapshot (oracle fits, unused at runtime)")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--explain", action="store_true", help="print derivations")
    args = ap.parse_args()

    rows = _csv_by_name()
    problems: list[str] = []
    if args.explain:
        print("Uranian (live vs data/fictitious_orbits.csv):")
    check_uranian(rows, problems, args.explain)
    check_proserpina(rows, problems, args.explain)
    check_transpluto(rows, problems, args.explain)
    check_legacy(problems, args.explain)

    for p in problems:
        print(f"MISMATCH: {p}")
    print(f"hypothetical provenance: {len(problems)} mismatch(es) (gate: 0)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
