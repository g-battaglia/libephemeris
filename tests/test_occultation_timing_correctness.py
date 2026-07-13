"""Correctness of lunar-occultation timing.

The occultation tests only check structure/length. This pins the time of greatest
occultation to an independent geometric criterion: it is the instant the
geocentric apparent Moon-body angular separation is minimum. We locate that
minimum independently (golden-section on the library's own apparent positions)
and require the reported time to coincide to a few seconds, with the separation
at a genuine local minimum.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as L
from libephemeris.constants import MOON, VENUS, MARS, FLG_SWIEPH, FLG_EQUATORIAL


def _sep(jd_ut, sid):
    m, _ = L.calc_ut(jd_ut, MOON, FLG_SWIEPH | FLG_EQUATORIAL)
    b, _ = L.calc_ut(jd_ut, sid, FLG_SWIEPH | FLG_EQUATORIAL)
    ra1, d1, ra2, d2 = map(math.radians, (m[0], m[1], b[0], b[1]))
    c = math.sin(d1) * math.sin(d2) + math.cos(d1) * math.cos(d2) * math.cos(ra1 - ra2)
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def _golden_min(jd0, sid, half=0.04):
    gr = (math.sqrt(5) - 1) / 2
    a, b = jd0 - half, jd0 + half
    c, d = b - gr * (b - a), a + gr * (b - a)
    fc, fd = _sep(c, sid), _sep(d, sid)
    for _ in range(40):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = _sep(c, sid)
        else:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = _sep(d, sid)
    jm = (a + b) / 2
    return jm, _sep(jm, sid)


@pytest.fixture(autouse=True)
def _skyfield():
    L.set_calc_mode("skyfield")
    yield
    L.set_calc_mode("auto")


@pytest.mark.slow
@pytest.mark.parametrize("body", [VENUS, MARS])
def test_lun_occult_time_is_separation_minimum(body):
    jd0 = L.julday(2020, 1, 1, 0.0, L.GREG_CAL)
    _, tret = L.lun_occult_when_glob(jd0, body, FLG_SWIEPH, 0, False)
    jdmax = tret[0]
    jm, minsep = _golden_min(jdmax, body)
    # Reported max == independent geocentric separation minimum (within ~30 s).
    assert abs(jdmax - jm) * 86400.0 < 30.0, f"body {body}: {(jdmax - jm) * 86400:.1f}s"
    # It is a genuine minimum: the separation grows on both sides.
    assert _sep(jdmax - 0.02, body) > minsep
    assert _sep(jdmax + 0.02, body) > minsep
    # A (parallax-enabled global) lunar occultation: geocentric approach within
    # ~1.5 deg (Moon radius + horizontal parallax).
    assert minsep < 1.5, f"body {body}: min separation {minsep:.3f} deg too large"
