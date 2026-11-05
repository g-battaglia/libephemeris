"""Crossing roots retain the same sidereal and observer coordinates (#15/#16)."""

from data_paths import ADMITTED_CORE_HASHES, CORE_LEB
from pathlib import Path
import hashlib
import json
import math
import subprocess

import libephemeris as eph
import pytest

FIXTURE = CORE_LEB
LIVE_MODES = (14, 17, 40, 36, 30, 31, 33, 32, 27, 35, 29, 28, 39)


@pytest.fixture(autouse=True)
def configured():
    assert hashlib.sha256(FIXTURE.read_bytes()).hexdigest() in ADMITTED_CORE_HASHES
    eph.set_calc_mode("leb")
    eph.set_precision_tier("base")
    eph.set_leb_file(str(FIXTURE))
    eph.set_delta_t_userdef(64.0 / 86400.0)
    eph.set_iers_delta_t_enabled(False)
    eph.set_sid_mode(eph.SIDM_LAHIRI)
    eph.set_topo(12.5, 41.9, 80.0)
    yield
    eph.set_sid_mode(0)


def residual(longitude, target):
    return abs(math.remainder(longitude - target, 360.0)) * 3600


@pytest.mark.parametrize(
    "search,calculate,body,extra",
    [
        ("solcross", "calc", 0, 0),
        ("solcross_ut", "calc_ut", 0, 0),
        ("mooncross", "calc", 1, 0),
        ("mooncross_ut", "calc_ut", 1, 0),
        ("helio_cross", "calc", 4, 8),
        ("helio_cross_ut", "calc_ut", 4, 8),
    ],
)
@pytest.mark.parametrize("frame", [65536, 32768, 65536 | 32768])
@pytest.mark.parametrize("backwards", [False, True])
def test_crossings_close_in_the_configured_coordinate_system(
    search, calculate, body, extra, frame, backwards
):
    flags = 2 | 256 | frame | extra
    calculate = getattr(eph, calculate)
    target = calculate(2448058.0, body, flags)[0][0]
    start = 2460310.5
    args = ([body] if search.startswith("helio") else []) + [
        target,
        start,
        flags,
        backwards,
    ]
    root = getattr(eph, search)(*args)
    assert root < start if backwards else root >= start
    assert residual(calculate(root, body, flags)[0][0], target) <= 0.001


@pytest.mark.parametrize("frame", [65536, 32768, 65536 | 32768])
def test_generic_crossing_uses_the_same_context(frame):
    flags = 2 | 256 | frame
    target = eph.calc_ut(2448058.0, eph.SUN, flags)[0][0]
    root = eph.cross_ut(eph.SUN, target, 2460310.5, flags)
    assert root >= 2460310.5
    assert residual(eph.calc_ut(root, eph.SUN, flags)[0][0], target) <= 0.001


@pytest.mark.parametrize(
    "search,calculate", [("mooncross_node", "calc"), ("mooncross_node_ut", "calc_ut")]
)
@pytest.mark.parametrize("backwards", [False, True])
def test_node_search_returns_its_actual_topocentric_sidereal_coordinates(
    search, calculate, backwards
):
    flags = 2 | 256 | 65536 | 32768
    time, longitude, latitude = getattr(eph, search)(2460310.5, flags, backwards)
    values = getattr(eph, calculate)(time, eph.MOON, flags)[0]
    assert abs(values[1]) * 3600 <= 0.001
    assert residual(values[0], longitude) <= 0.001
    assert abs(values[1] - latitude) * 3600 <= 0.001


@pytest.mark.parametrize("body,search", [(0, "solcross_ut"), (1, "mooncross_ut")])
def test_crossing_really_consumes_the_observers_altitude(body, search):
    flags = 2 | 256 | 32768
    target = eph.calc_ut(2448058.0, body, flags)[0][0]
    roots = []
    for altitude in (0.0, 10000.0):
        eph.set_topo(12.5, 41.9, altitude)
        root = getattr(eph, search)(target, 2460310.5, flags)
        assert residual(eph.calc_ut(root, body, flags)[0][0], target) <= 0.001
        roots.append(root)
    assert roots[0] != roots[1]


def test_crossing_module_exposes_actual_native_solver_constants():
    from libephemeris import crossing

    assert (
        crossing.NR_TOLERANCE
        == crossing.NR_TOLERANCE_SUN
        == crossing.NR_TOLERANCE_MOON
        == 0.001 / 3600
    )
    assert crossing.STATION_SPEED_THRESHOLD == 0.01
    assert crossing.STATION_VELOCITY_TOLERANCE == 1e-6
    assert [
        crossing.NR_MAX_ITER_SUN,
        crossing.NR_MAX_ITER_MOON,
        crossing.NR_MAX_ITER_PLANET,
        crossing.NR_MAX_ITER_HELIO,
        crossing.NR_MAX_ITER_VERY_SLOW,
        crossing.NR_MAX_ITER_STATION,
    ] == [50, 50, 60, 80, 80, 100]
    assert crossing.calc is eph.calc and crossing.calc_ut is eph.calc_ut
    assert (
        crossing.Error is eph.Error
        and crossing.CalculationError is eph.CalculationError
    )
    assert crossing.EphemerisRangeError is eph.EphemerisRangeError
    assert (
        crossing.FLG_SWIEPH,
        crossing.FLG_SPEED,
        crossing.FLG_HELCTR,
        crossing.FLG_RADIANS,
        crossing.SUN,
        crossing.MOON,
    ) == (2, 256, 8, 8192, 0, 1)
