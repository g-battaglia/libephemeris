#!/usr/bin/env python3
"""Verify the exact hypothetical-body data and provenance classifications.

This is an integrity gate, not a blanket exception.  It pins the reviewed CSV,
checks every runtime element against an explicit expected value, requires a
per-ID provenance classification, and verifies that all compatibility IDs
40--58 calculate finite states.  ``built-in`` rows are the project's shipped
element sets; the gate must never silently relabel them ``primary-derived``.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import sys
from pathlib import Path

from libephemeris import hypothetical as hyp

_CSV_SHA256 = "e34bcc0508bda6e24763dddbb94d95028022f14763ea9235f5d705e607298c7b"

_ROW_NAMES = [
    "Cupido",
    "Hades",
    "Zeus",
    "Kronos",
    "Apollon",
    "Admetos",
    "Vulkanus",
    "Poseidon",
    "Isis-Transpluto",
    "Nibiru",
    "Harrington",
    "Leverrier-Neptune",
    "Adams-Neptune",
    "Lowell-Pluto",
    "Pickering-Pluto",
    "Vulcan",
    "Selena",
    "Proserpina",
    "Waldemath",
]


_URANIAN = {
    hyp.CUPIDO: (
        2415020.0,
        2415020.0,
        40.99837,
        0.00460,
        1.0833,
        171.4333,
        129.8325,
        163.7409,
    ),
    hyp.HADES: (
        2415020.0,
        2415020.0,
        50.66744,
        0.00245,
        1.0500,
        148.1796,
        161.3339,
        27.6496,
    ),
    hyp.ZEUS: (2415020.0, 2415020.0, 59.21436, 0.00120, 0.0, 299.0440, 0.0, 165.1232),
    hyp.KRONOS: (2415020.0, 2415020.0, 64.81690, 0.00305, 0.0, 208.8801, 0.0, 169.0193),
    hyp.APOLLON: (2415020.0, 2415020.0, 70.29949, 0.0, 0.0, 0.0, 0.0, 138.0533),
    hyp.ADMETOS: (2415020.0, 2415020.0, 73.62765, 0.0, 0.0, 0.0, 0.0, 351.3350),
    hyp.VULKANUS: (2415020.0, 2415020.0, 77.25568, 0.0, 0.0, 0.0, 0.0, 55.8983),
    hyp.POSEIDON: (2415020.0, 2415020.0, 83.66907, 0.0, 0.0, 0.0, 0.0, 165.5163),
}

_LEGACY_URANIAN = {
    hyp.CUPIDO: ("Cupido", 241.2067, 1.091437, 1.5, 21.6, 83.0),
    hyp.HADES: ("Hades", 176.0581, 0.736380, 1.5, 258.0, 56.0),
    hyp.ZEUS: ("Zeus", 32.1893, 0.532664, 1.0, 141.0, 41.0),
    hyp.KRONOS: ("Kronos", 213.2096, 0.420481, 1.0, 98.0, 32.0),
    hyp.APOLLON: ("Apollon", 71.8925, 0.341403, 1.0, 326.0, 26.0),
    hyp.ADMETOS: ("Admetos", 142.2269, 0.283756, 1.0, 250.0, 22.0),
    hyp.VULKANUS: ("Vulkanus", 195.6753, 0.240116, 1.0, 178.0, 18.0),
    hyp.POSEIDON: ("Poseidon", 274.4073, 0.207016, 1.0, 105.0, 16.0),
}

_LEGACY_KEPLERIAN_PUBLIC = {
    "CUPIDO_KEPLERIAN_ELEMENTS": (
        "CupidoKeplerianElements",
        {
            "name": "Cupido",
            "epoch": 2415020.0,
            "a": 40.99837,
            "e": 0.0,
            "i": 0.0,
            "omega": 0.0,
            "Omega": 0.0,
            "L0": 105.301693,
            "n": 0.0037945179,
        },
    ),
    "HADES_KEPLERIAN_ELEMENTS": (
        "HadesKeplerianElements",
        {
            "name": "Hades",
            "epoch": 2415020.0,
            "a": 50.66744,
            "e": 0.00245,
            "i": 1.05,
            "omega": 148.1796,
            "Omega": 161.3339,
            "M0": 26.850162,
            "n": 0.00278759,
        },
    ),
    "ZEUS_KEPLERIAN_ELEMENTS": (
        "ZeusKeplerianElements",
        {
            "name": "Zeus",
            "epoch": 2415020.0,
            "a": 59.21436,
            "e": 0.0,
            "i": 0.0,
            "omega": 0.0,
            "Omega": 0.0,
            "L0": 104.289095,
            "n": 0.002220375,
        },
    ),
    "KRONOS_KEPLERIAN_ELEMENTS": (
        "KronosKeplerianElements",
        {
            "name": "Kronos",
            "epoch": 2415020.0,
            "a": 64.8169,
            "e": 0.0,
            "i": 0.0,
            "omega": 0.0,
            "Omega": 0.0,
            "L0": 17.111353,
            "n": 0.0019351856,
        },
    ),
    "APOLLON_KEPLERIAN_ELEMENTS": (
        "ApollonKeplerianElements",
        {
            "name": "Apollon",
            "epoch": 2415020.0,
            "a": 70.36118,
            "e": 0.0,
            "i": 0.0,
            "omega": 0.0,
            "Omega": 0.0,
            "L0": 138.565328,
            "n": 0.0017177599,
        },
    ),
    "ADMETOS_KEPLERIAN_ELEMENTS": (
        "AdmetosKeplerianElements",
        {
            "name": "Admetos",
            "epoch": 2415020.0,
            "a": 73.736396,
            "e": 0.0,
            "i": 0.0,
            "omega": 0.0,
            "Omega": 0.0,
            "L0": 350.613913,
            "n": 0.0016016766,
        },
    ),
    "VULKANUS_KEPLERIAN_ELEMENTS": (
        "VulkanusKeplerianElements",
        {
            "name": "Vulkanus",
            "epoch": 2415020.0,
            "a": 77.445895,
            "e": 0.0,
            "i": 0.0,
            "omega": 0.0,
            "Omega": 0.0,
            "L0": 55.397715,
            "n": 0.0015069325,
        },
    ),
    "POSEIDON_KEPLERIAN_ELEMENTS": (
        "PoseidonKeplerianElements",
        {
            "name": "Poseidon",
            "epoch": 2415020.0,
            "a": 83.666307,
            "e": 0.0,
            "i": 0.0,
            "omega": 0.0,
            "Omega": 0.0,
            "L0": 166.140256,
            "n": 0.0013256078,
        },
    ),
}

_CLASSICAL = {
    hyp.NIBIRU: (
        1856113.380954,
        1856113.380954,
        0.0,
        234.8921,
        0.981092,
        103.966,
        -44.567,
        158.708,
    ),
    hyp.HARRINGTON: (2374696.5, 2451545.0, 0.0, 101.2, 0.411, 208.5, 275.4, 32.4),
    hyp.NEPTUNE_LEVERRIER: (
        2395662.5,
        2395662.5,
        34.05,
        36.15,
        0.10761,
        284.75,
        0.0,
        0.0,
    ),
    hyp.NEPTUNE_ADAMS: (2395662.5, 2395662.5, 24.28, 37.25, 0.12062, 299.11, 0.0, 0.0),
    hyp.PLUTO_LOWELL: (2425977.5, 2425977.5, 281.0, 43.0, 0.202, 204.9, 0.0, 0.0),
    hyp.PLUTO_PICKERING: (2425977.5, 2425977.5, 48.95, 55.1, 0.31, 280.1, 100.0, 15.0),
}


def _check_equal(problems: list[str], label: str, got: object, want: object) -> None:
    if got != want:
        problems.append(f"{label}: {got!r} != {want!r}")


def _csv_sources(path: Path) -> dict[str, str]:
    sources: dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(",", 10)
        if len(parts) != 11:
            raise ValueError(f"CSV row has {len(parts)} fields: {line!r}")
        sources[parts[0].strip()] = parts[10].strip()
    return sources


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--explain", action="store_true")
    args = parser.parse_args()
    problems: list[str] = []

    csv_path = hyp.get_bundled_fictitious_orbits_path()
    digest = hashlib.sha256(csv_path.read_bytes()).hexdigest()
    _check_equal(problems, "fictitious_orbits.csv SHA-256", digest, _CSV_SHA256)

    rows = hyp.load_bundled_fictitious_orbits()
    _check_equal(problems, "bundled row order", [row.name for row in rows], _ROW_NAMES)
    try:
        sources = _csv_sources(csv_path)
    except ValueError as exc:
        problems.append(str(exc))
        sources = {}
    for name, source in sources.items():
        if not source.startswith("built-in:") and not source.startswith(
            "primary-derived:"
        ):
            problems.append(
                f"CSV {name} must start with 'built-in:' or 'primary-derived:', got {source!r}"
            )

    all_ids = set(range(hyp.CUPIDO, hyp.WALDEMATH + 1))
    _check_equal(problems, "provenance IDs", set(hyp.HYPOTHETICAL_PROVENANCE), all_ids)
    for body_id in sorted(all_ids):
        status, citation = hyp.HYPOTHETICAL_PROVENANCE[body_id]
        if status not in ("built-in", "primary-derived"):
            problems.append(
                f"body {body_id} status must be 'built-in' or 'primary-derived', got {status!r}"
            )
        if not citation.strip():
            problems.append(f"body {body_id} has an empty citation/review note")

    _check_equal(
        problems,
        "legacy Uranian IDs",
        set(hyp.URANIAN_ELEMENTS),
        set(_LEGACY_URANIAN),
    )
    for body_id, expected_legacy_uranian in _LEGACY_URANIAN.items():
        legacy_element = hyp.URANIAN_ELEMENTS[body_id]
        _check_equal(
            problems,
            f"legacy Uranian body {body_id} type",
            type(legacy_element),
            hyp.UranianElements,
        )
        legacy_values = (
            legacy_element.name,
            legacy_element.L0,
            legacy_element.n,
            legacy_element.amplitude,
            legacy_element.phase,
            legacy_element.phase_rate,
        )
        _check_equal(
            problems,
            f"legacy Uranian body {body_id} values",
            legacy_values,
            expected_legacy_uranian,
        )

    for object_name, (
        class_name,
        expected_legacy_keplerian,
    ) in _LEGACY_KEPLERIAN_PUBLIC.items():
        legacy_class = getattr(hyp, class_name, None)
        legacy_object = getattr(hyp, object_name, None)
        if legacy_class is None:
            problems.append(f"missing public class {class_name}")
            continue
        if legacy_object is None:
            problems.append(f"missing public constant {object_name}")
            continue
        _check_equal(
            problems,
            f"{object_name} type",
            type(legacy_object),
            legacy_class,
        )
        _check_equal(
            problems,
            f"{object_name} values",
            vars(legacy_object),
            expected_legacy_keplerian,
        )

    _check_equal(
        problems, "Uranian IDs", set(hyp.URANIAN_KEPLERIAN_ELEMENTS), set(_URANIAN)
    )
    for body_id, expected_uranian in _URANIAN.items():
        uranian_element = hyp.URANIAN_KEPLERIAN_ELEMENTS[body_id]
        uranian_values = (
            uranian_element.epoch,
            uranian_element.equinox_jd,
            uranian_element.a,
            uranian_element.e,
            uranian_element.i,
            uranian_element.omega,
            uranian_element.Omega,
            uranian_element.M0,
        )
        _check_equal(
            problems, f"Uranian body {body_id}", uranian_values, expected_uranian
        )
        _check_equal(
            problems,
            f"Uranian body {body_id} n",
            uranian_element.n,
            0.9856076686 / uranian_element.a**1.5,
        )

    _check_equal(
        problems, "classical IDs", set(hyp.FICTITIOUS_ORBITAL_ELEMENTS), set(_CLASSICAL)
    )
    for body_id, expected_classical in _CLASSICAL.items():
        classical_element = hyp.FICTITIOUS_ORBITAL_ELEMENTS[body_id]
        classical_values = (
            classical_element.epoch_jd,
            classical_element.equinox_jd,
            classical_element.mean_anomaly.evaluate(0.0),
            classical_element.semi_axis,
            classical_element.eccentricity.evaluate(0.0),
            classical_element.arg_perihelion.evaluate(0.0),
            classical_element.asc_node.evaluate(0.0),
            classical_element.inclination.evaluate(0.0),
        )
        _check_equal(
            problems,
            f"classical body {body_id}",
            classical_values,
            expected_classical,
        )

    transpluto = hyp.TRANSPLUTO_KEPLERIAN_ELEMENTS
    _check_equal(
        problems,
        "Transpluto public type",
        type(transpluto),
        hyp.TransplutoKeplerianElements,
    )
    _check_equal(
        problems,
        "Transpluto",
        (
            transpluto.name,
            transpluto.epoch,
            transpluto.a,
            transpluto.e,
            transpluto.i,
            transpluto.omega,
            transpluto.Omega,
            transpluto.M0,
        ),
        ("Transpluto", 2451545.0, 77.775, 0.3, 0.0, 1.468282, 0.0, 119.265936),
    )
    _check_equal(
        problems,
        "generic IDs",
        set(hyp.HYPOTHETICAL_ELEMENTS),
        {hyp.ISIS, hyp.PROSERPINA},
    )
    _check_equal(
        problems,
        "Vulcan constant",
        vars(hyp.VULCAN_ELEMENTS),
        {
            "name": "Vulcan",
            "epoch": 2415020.0,
            "a": 0.13744,
            "e": 0.019,
            "i": 7.5,
            "M0": 252.8987988,
            "n_century": 707550.7341,
            "omega0": 322.212069,
            "omega_rate": 1670.056,
            "Omega0": 47.787931,
            "Omega_rate": -1670.056,
        },
    )
    _check_equal(
        problems,
        "Pickering public constant exists",
        isinstance(hyp.PICKERING_PLANET_X_ELEMENTS, hyp.PickeringPlanetXElements),
        True,
    )
    _check_equal(
        problems,
        "Waldemath public constant exists",
        isinstance(hyp.WALDEMATH_ELEMENTS, hyp.WaldemathElements),
        True,
    )

    jd = 2451545.0
    for body_id in sorted(all_ids):
        state = hyp.calc_hypothetical_position(body_id, jd)
        if len(state) != 6 or not all(math.isfinite(value) for value in state):
            problems.append(f"body {body_id} did not return a finite six-value state")
        elif not (0.0 <= state[0] < 360.0 and state[2] > 0.0):
            problems.append(
                f"body {body_id} returned invalid polar coordinates {state!r}"
            )

    for white_jd in (2415020.0, 2451545.0, 2488070.0):
        white = hyp.calc_white_moon_position(white_jd)
        centuries = (white_jd - 2451545.0) / 36525.0
        expected_white = (
            (242.2205555 + 5143.5418158 * centuries) % 360.0,
            0.0,
            0.05280098949,
            5143.5418158 / 36525.0,
            0.0,
            0.0,
        )
        _check_equal(problems, f"White Moon at {white_jd}", white, expected_white)

    if args.explain:
        print(f"Pinned CSV: {_CSV_SHA256}")
        for body_id in sorted(all_ids):
            print(body_id, *hyp.HYPOTHETICAL_PROVENANCE[body_id])
    for problem in problems:
        print(f"MISMATCH: {problem}")
    print(f"hypothetical provenance: {len(problems)} mismatch(es) (gate: 0)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
