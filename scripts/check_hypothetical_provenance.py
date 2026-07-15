#!/usr/bin/env python3
"""Verify the fail-closed, independently sourced hypothetical-data boundary.

This gate checks *semantic provenance invariants*, not merely that a file kept
the same bytes.  It pins the reviewed CSV, independently reconstructs every
enabled runtime row from the published quantities documented in
``libephemeris.hypothetical``, verifies the Neely source-rate consistency, and
reconstructs the published Selena convention before calling every public ID to
prove that unsupported bodies fail closed.

The reviewed sources themselves are not vendored.  Their identifiers and the
SHA-256 digests of the scans used during transcription are printed by
``--explain`` so a release reviewer can obtain and compare the public scans
without placing research copies in the repository.

Provenance:
    Project-authored integrity and derivation gate. Its expected values are
    independently recomputed from the page-level public quantities catalogued
    in ``docs/methodology/hypothetical-bodies.md``; scan hashes identify review
    evidence but do not license or embed it. The gate also proves unsupported
    records fail closed. It does not compare against or learn from another
    ephemeris implementation.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import sys
from dataclasses import fields

from libephemeris import hypothetical as hyp
from libephemeris.exceptions import UnknownBodyError


_CSV_SHA256 = "7b35add4e93b5b4dd5473d0eb4bdbcb6f914e7f22e99a2e8a0b8953dd6ed4308"

_REVIEWED_SOURCE_DIGESTS = {
    "Neely 1980, Matrix Magazine VII, complete 63-page scan": (
        "f3d2e834a4a684736383ec387c9c2786064ca4668256724865da2f3d8fa03020"
    ),
    "Harrington 1988, AJ 96, pp. 1476-1478": (
        "c0f74fda659cb5fbd7c4db388b74c3169dc269d481338a1f6c4f163f3ab0d49c"
    ),
    "Adams 1846, Explanation of the Observed Irregularities of Uranus": (
        "451b7668369911d6241128616a4dc6a108320e84a0b34bf9f4fbbfd6c6e873ec"
    ),
}

# The Velichko--Larin volume expressly restricts reproduction.  Review used
# browser page images only; no scan is vendored and therefore no repository
# hash is asserted for it.  The locator still makes the factual inputs
# independently reviewable: uniform seven-year motion (p. 17), exclusion of
# the alternative seven-years-minus-fifteen-days cycle (p. 18), January 1800
# (p. 20), March 1879 (p. 29), and January/December 2000 plus January 2007
# (p. 45).
_REVIEWED_SOURCE_LOCATORS = {
    "Velichko and Larin 2007, ISBN 5-900191-12-5, pp. 17, 18, 20, 29, 45": (
        "https://djvu.online/file/lBJKQiKf9n4Gd"
    ),
    "IAU 2015 Resolution B3, nominal terrestrial mass parameter": (
        "https://www.iau.org/static/resolutions/IAU2015_English.pdf"
    ),
    "IAU 2012 Resolution B2, exact astronomical unit": (
        "https://www.iau.org/static/resolutions/IAU2012_English.pdf"
    ),
}

_SUPPORTED = {
    *range(hyp.CUPIDO, hyp.POSEIDON + 1),
    hyp.HARRINGTON,
    hyp.NEPTUNE_LEVERRIER,
    hyp.NEPTUNE_ADAMS,
    hyp.PLUTO_LOWELL,
    hyp.WHITE_MOON,
}
_ALL_IDS = set(range(hyp.CUPIDO, hyp.WALDEMATH + 1))
_UNSUPPORTED = _ALL_IDS - _SUPPORTED


# Literal CSV expectations.  Each tuple is:
# (epoch, equinox, M, a, e, omega, node, i, geocentric).
_EXPECTED_CSV_ROWS = {
    "Cupido": (
        2415020.0,
        2415020.0,
        163.7409,
        40.99837,
        0.00460,
        171.4333,
        129.8325,
        1.0833,
        False,
    ),
    "Hades": (
        2415020.0,
        2415020.0,
        27.6496,
        50.66744,
        0.00245,
        148.1796,
        161.3339,
        1.0500,
        False,
    ),
    "Zeus": (
        2415020.0,
        2415020.0,
        165.1232,
        59.21436,
        0.00120,
        299.0440,
        0.0,
        0.0,
        False,
    ),
    "Kronos": (
        2415020.0,
        2415020.0,
        169.0193,
        64.81690,
        0.00305,
        208.8801,
        0.0,
        0.0,
        False,
    ),
    "Apollon": (2415020.0, 2415020.0, 0.0, 70.29949, 0.0, 138.0533, 0.0, 0.0, False),
    "Admetos": (2415020.0, 2415020.0, 0.0, 73.62765, 0.0, 351.3350, 0.0, 0.0, False),
    "Vulkanus": (2415020.0, 2415020.0, 0.0, 77.25568, 0.0, 55.8983, 0.0, 0.0, False),
    "Poseidon": (2415020.0, 2415020.0, 0.0, 83.66907, 0.0, 165.5163, 0.0, 0.0, False),
    "Harrington": (2374696.5, 2451545.0, 0.0, 101.2, 0.411, 208.5, 275.4, 32.4, False),
    "Leverrier-Neptune": (
        2395662.5,
        2395662.5,
        34.0 + 1.0 / 60.0 + 56.0 / 3600.0,
        36.1539,
        0.10761,
        284.75,
        0.0,
        0.0,
        False,
    ),
    "Adams-Neptune": (
        2395575.5,
        2395575.5,
        23.85,
        37.25,
        0.120615,
        299.0 + 11.0 / 60.0,
        0.0,
        0.0,
        False,
    ),
    "Lowell-Pluto": (2396757.5, 2396757.5, 178.3, 43.0, 0.202, 203.8, 0.0, 0.0, False),
}


def _check(problems: list[str], condition: bool, message: str) -> None:
    """Append ``message`` when a provenance invariant is false."""
    if not condition:
        problems.append(message)


def _close(actual: float, expected: float, *, atol: float = 1e-12) -> bool:
    """Compare transcription values without hiding material rounding drift."""
    return math.isclose(actual, expected, rel_tol=0.0, abs_tol=atol)


def _all_numeric_fields_are_nan(value: object) -> bool:
    """Return whether a public unsupported container is an explicit sentinel."""
    return all(
        item.name == "name"
        or (
            isinstance(getattr(value, item.name), float)
            and math.isnan(getattr(value, item.name))
        )
        for item in fields(value)
    )


def _row_tuple(row: hyp.OrbitalElements) -> tuple[float | bool, ...]:
    """Project a parsed CSV row onto the nine provenance-bearing fields."""
    assert row.equinox_jd is not None
    return (
        row.epoch_jd,
        row.equinox_jd,
        row.mean_anomaly.constant,
        row.semi_axis,
        row.eccentricity.constant,
        row.arg_perihelion.constant,
        row.asc_node.constant,
        row.inclination.constant,
        row.is_geocentric,
    )


def _check_csv(problems: list[str]) -> None:
    """Verify the packaged transcription byte-for-byte and field-for-field."""
    csv_path = hyp.get_bundled_fictitious_orbits_path()
    digest = hashlib.sha256(csv_path.read_bytes()).hexdigest()
    _check(problems, digest == _CSV_SHA256, f"unexpected CSV SHA-256: {digest}")

    rows = hyp.load_bundled_fictitious_orbits()
    _check(
        problems,
        len(rows) == len(_EXPECTED_CSV_ROWS),
        f"bundled row count is {len(rows)}, expected {len(_EXPECTED_CSV_ROWS)}",
    )
    actual_by_name = {row.name: row for row in rows}
    _check(
        problems,
        set(actual_by_name) == set(_EXPECTED_CSV_ROWS),
        "bundled row names differ from the reviewed set",
    )
    for name, expected in _EXPECTED_CSV_ROWS.items():
        row = actual_by_name.get(name)
        if row is None:
            continue
        actual = _row_tuple(row)
        _check(
            problems,
            all(
                a == e if isinstance(e, bool) else _close(float(a), float(e))
                for a, e in zip(actual, expected, strict=True)
            ),
            f"unexpected {name} row: {actual!r}",
        )


def _check_runtime_registries(problems: list[str]) -> None:
    """Reconstruct enabled runtime registries from their documented sources."""
    uranian_ids = set(range(hyp.CUPIDO, hyp.POSEIDON + 1))
    _check(
        problems,
        set(hyp.URANIAN_KEPLERIAN_ELEMENTS) == uranian_ids,
        "Uranian runtime registry differs from IDs 40-47",
    )
    _check(
        problems,
        set(hyp.URANIAN_ELEMENTS) == uranian_ids,
        "legacy Uranian source projection differs from IDs 40-47",
    )
    _check(
        problems,
        set(hyp.FICTITIOUS_ORBITAL_ELEMENTS)
        == {
            hyp.HARRINGTON,
            hyp.NEPTUNE_LEVERRIER,
            hyp.NEPTUNE_ADAMS,
            hyp.PLUTO_LOWELL,
        },
        "classical-prediction registry differs from IDs 50-53 excluding 54",
    )
    _check(
        problems,
        not hyp.HYPOTHETICAL_ELEMENTS,
        "unreviewed generic Keplerian registry is not empty",
    )

    # The primary Neely table prints rounded rates.  The runtime's standard
    # Gaussian propagation derived from a must agree within that print precision.
    for body_id, source in hyp._NEELY_SOURCE_ROWS.items():
        runtime = hyp.URANIAN_KEPLERIAN_ELEMENTS[body_id]
        _check(problems, runtime.name == source.name, f"body {body_id}: name drift")
        for field_name, actual, expected in (
            ("a", runtime.a, source.a),
            ("e", runtime.e, source.e),
            ("i", runtime.i, source.i),
            ("omega", runtime.omega, source.omega),
            ("node", runtime.Omega, source.node),
            ("M0", runtime.M0, source.mean_anomaly),
        ):
            _check(
                problems,
                _close(actual, expected),
                f"body {body_id}: {field_name} differs from Neely transcription",
            )
        rate_error = abs(runtime.n * 36525.0 - source.mean_anomaly_rate_century)
        _check(
            problems,
            rate_error < 0.003,
            f"body {body_id}: Gaussian/source rate error is {rate_error} deg/century",
        )

    # Independent arithmetic checks for the three transformed historical rows.
    leverrier = hyp.FICTITIOUS_ORBITAL_ELEMENTS[hyp.NEPTUNE_LEVERRIER]
    _check(
        problems,
        _close(leverrier.mean_anomaly.constant, 34.0 + 1.0 / 60.0 + 56.0 / 3600.0),
        "Le Verrier sexagesimal mean anomaly was not preserved",
    )
    adams = hyp.FICTITIOUS_ORBITAL_ELEMENTS[hyp.NEPTUNE_ADAMS]
    _check(
        problems,
        _close(
            adams.mean_anomaly.constant, (323.0 + 2.0 / 60.0) - (299.0 + 11.0 / 60.0)
        ),
        "Adams M != printed mean longitude minus printed perihelion longitude",
    )
    _check(
        problems,
        _close(adams.arg_perihelion.constant, 299.0 + 11.0 / 60.0),
        "Adams perihelion is not a correct sexagesimal conversion",
    )
    lowell = hyp.FICTITIOUS_ORBITAL_ELEMENTS[hyp.PLUTO_LOWELL]
    _check(
        problems,
        _close(lowell.mean_anomaly.constant, (22.1 - 203.8) % 360.0),
        "Lowell M != printed epsilon minus printed varpi",
    )

    # Selena is reconstructed only from the publication's stated uniform rule
    # and two long-baseline one-minute checkpoints.  The source says the point
    # has an approximately seven-year cycle (p. 17), while p. 18 explicitly
    # excludes the separate seven-years-minus-fifteen-days convention from its
    # ephemerides.  Pages 20 and 45 give January 1800 and January 2000 phases.
    # Their identical month labels make the 73,048-day separation independent
    # of any fixed within-month sampling convention.  The seven-year statement
    # uniquely selects 28 complete turns; 29 would imply only about 6.76 years.
    early_epoch = 2378496.5
    early_longitude = 30.0 + 6.0 + 8.0 / 60.0
    source_epoch = 2451544.5
    source_longitude = 240.0 + 2.0 + 18.0 / 60.0
    expected_elapsed_days = source_epoch - early_epoch
    expected_displacement = 28.0 * 360.0 + source_longitude - early_longitude
    expected_rate = expected_displacement / expected_elapsed_days
    expected_period_days = 360.0 / expected_rate
    _check(
        problems,
        _close(hyp._WHITE_MOON_EARLY_SOURCE_EPOCH_JD, early_epoch),
        "Selena early source epoch is not 1800-01-01 00:00 TT",
    )
    _check(
        problems,
        _close(hyp._WHITE_MOON_EARLY_SOURCE_LONGITUDE_DEG, early_longitude),
        "Selena early source phase is not p. 20's 6 deg 08 arcmin Taurus",
    )
    _check(
        problems,
        _close(hyp._WHITE_MOON_SOURCE_EPOCH_JD, source_epoch),
        "Selena source epoch is not 2000-01-01 00:00 TT",
    )
    _check(
        problems,
        _close(
            hyp._WHITE_MOON_SOURCE_LONGITUDE_DEG,
            source_longitude,
        ),
        "Selena source phase is not p. 45's 2 deg 18 arcmin Sagittarius",
    )
    _check(
        problems,
        hyp._WHITE_MOON_COMPLETE_TURNS == 28,
        "Selena source unwrap is not the unique 28-turn seven-year solution",
    )
    _check(
        problems,
        _close(hyp._WHITE_MOON_SOURCE_ELAPSED_DAYS, expected_elapsed_days),
        "Selena source checkpoint interval is not 73,048 calendar days",
    )
    _check(
        problems,
        _close(hyp._WHITE_MOON_SOURCE_DISPLACEMENT_DEG, expected_displacement),
        "Selena unwrapped source displacement is not 10,286 deg 10 arcmin",
    )
    _check(
        problems,
        _close(hyp._WHITE_MOON_RATE_DEG_PER_DAY, expected_rate),
        "Selena rate is not the documented 1800--2000 endpoint derivation",
    )
    _check(
        problems,
        _close(hyp._WHITE_MOON_PERIOD_DAYS, expected_period_days),
        "Selena period is not the reciprocal of the source-derived rate",
    )

    # The first two rows are the defining endpoints.  The last three are
    # independent checkpoints distributed across the source's date range.  A
    # correctly unwrapped uniform model must reproduce the printed values to
    # half an arcminute, the maximum rounding error of a nearest-minute table.
    for jd, published_longitude in (
        (early_epoch, early_longitude),
        (source_epoch, source_longitude),
        (2407409.5, 120.0 + 27.0 + 29.0 / 60.0),
        (2451879.5, 270.0 + 19.0 + 28.0 / 60.0),
        (2454101.5, 240.0 + 2.0 + 22.0 / 60.0),
    ):
        calculated = hyp.calc_white_moon_position(jd)[0]
        error = abs((calculated - published_longitude + 180.0) % 360.0 - 180.0)
        _check(
            problems,
            error <= 0.5 / 60.0,
            f"Selena published-checkpoint error is {error} deg at JD {jd}",
        )

    # The source supplies no physical radius.  Reconstruct the explicitly
    # conventional display distance from IAU nominal SI constants, rather than
    # accepting an unexplained compatibility-table distance.
    reconstructed_distance = (
        hyp._NOMINAL_EARTH_GM_M3_S2
        * (hyp._WHITE_MOON_PERIOD_DAYS * hyp._SECONDS_PER_DAY / (2.0 * math.pi)) ** 2
    ) ** (1.0 / 3.0) / hyp._ASTRONOMICAL_UNIT_M
    _check(
        problems,
        _close(hyp._WHITE_MOON_DISTANCE_AU, reconstructed_distance),
        "Selena display distance is not the documented IAU/Kepler derivation",
    )


def _check_fail_closed_behavior(problems: list[str]) -> None:
    """Exercise all recognised IDs and verify their declared support status."""
    _check(
        problems,
        set(hyp.HYPOTHETICAL_PROVENANCE) == _ALL_IDS,
        "provenance registry does not cover every compatibility ID",
    )
    for body_id in sorted(_SUPPORTED):
        status, reason = hyp.HYPOTHETICAL_PROVENANCE.get(body_id, ("", ""))
        expected_status = (
            "published-model" if body_id == hyp.WHITE_MOON else "primary-transcription"
        )
        _check(
            problems,
            status == expected_status and bool(reason.strip()),
            f"supported body {body_id} lacks its {expected_status} reason",
        )
        try:
            state = hyp.calc_hypothetical_position(body_id, 2451545.0)
        except Exception as exc:  # pragma: no cover - diagnostic gate
            problems.append(f"supported body {body_id} raised {type(exc).__name__}")
        else:
            _check(
                problems,
                len(state) == 6
                and all(math.isfinite(component) for component in state),
                f"supported body {body_id} did not return six finite values",
            )

    for body_id in sorted(_UNSUPPORTED):
        status, reason = hyp.HYPOTHETICAL_PROVENANCE.get(body_id, ("", ""))
        _check(
            problems,
            status == "unsupported" and bool(reason.strip()),
            f"body {body_id} is not explicitly unsupported with a reason",
        )
        try:
            hyp.calc_hypothetical_position(body_id, 2451545.0)
        except UnknownBodyError:
            pass
        except Exception as exc:  # pragma: no cover - diagnostic gate
            problems.append(
                f"body {body_id} raised {type(exc).__name__}, expected UnknownBodyError"
            )
        else:
            problems.append(f"unsupported body {body_id} unexpectedly calculated")

    public_sentinels = (
        hyp.TRANSPLUTO_KEPLERIAN_ELEMENTS,
        hyp.PICKERING_PLANET_X_ELEMENTS,
        hyp.VULCAN_ELEMENTS,
        hyp.WALDEMATH_ELEMENTS,
    )
    for value in public_sentinels:
        _check(
            problems,
            _all_numeric_fields_are_nan(value),
            f"unsupported public container {type(value).__name__} contains data",
        )


def main() -> int:
    """Run the provenance gate and return a process exit status."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--explain", action="store_true")
    args = parser.parse_args()
    problems: list[str] = []

    _check_csv(problems)
    _check_runtime_registries(problems)
    _check_fail_closed_behavior(problems)

    if args.explain:
        print(f"Pinned CSV: {_CSV_SHA256}")
        for source, digest in _REVIEWED_SOURCE_DIGESTS.items():
            print(f"Reviewed source: {source}; SHA-256 {digest}")
        for source, locator in _REVIEWED_SOURCE_LOCATORS.items():
            print(f"Reviewed browser source: {source}; {locator}")
        for body_id in sorted(_ALL_IDS):
            print(body_id, *hyp.HYPOTHETICAL_PROVENANCE[body_id])

    for problem in problems:
        print(f"MISMATCH: {problem}")
    print(f"hypothetical provenance: {len(problems)} mismatch(es) (gate: 0)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
