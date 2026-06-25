"""Regression: solar-eclipse obscuration is a true fraction in [0, 1].

Obscuration is the fraction of the Sun's *area* covered by the Moon, so it is
physically bounded by 1.0 (100% = fully covered). Audit round v10 found that for
TOTAL eclipses several code paths returned the Moon/Sun disc *area ratio*
(R_moon/R_sun)^2 ~ 1.12 > 1, which is not a fraction (it is the magnitude's
domain) and contradicted the functions' own docstrings ("1.0 for total
eclipse"). They now clamp a total eclipse to 1.0. Annular eclipses keep the
area ratio (a ring of Sun remains, so < 1) and partial eclipses use the
two-disc lens overlap -- both already correct.

This deliberately diverges from the upstream reference for total eclipses (the
reference returns ~1.12); the divergence favours physical correctness and is
documented in docs/comparison/known-differences.md.

Checks (independent of any external implementation):
  * total: obscuration == 1.0 from every public entry point (cross-API),
  * obscuration <= 1.0 at every sampled location/time,
  * annular: 0 < obscuration < 1 and equals the squared diameter ratio,
  * partial: 0 < obscuration < 1.
"""

from __future__ import annotations

import pytest

from libephemeris import (
    FLG_SWIEPH,
    julday,
    sol_eclipse_how,
    sol_eclipse_how_details,
    sol_eclipse_obscuration_at_loc,
    sol_eclipse_when_glob,
    sol_eclipse_when_loc,
    sol_eclipse_where,
)

pytestmark = pytest.mark.slow


def _glob_max(year, month, day):
    jd0 = julday(year, month, day, 0.0)
    _rf, tret = sol_eclipse_when_glob(jd0, FLG_SWIEPH, 0, False)
    return tret[0]


def _local_max(jd_near, lon, lat):
    _rc, tloc, _attr = sol_eclipse_when_loc(
        jd_near - 0.2, (lon, lat, 0.0), FLG_SWIEPH, False
    )
    return tloc[0]


def test_total_eclipse_obscuration_is_one_all_apis():
    """A total eclipse fully covers the Sun: obscuration == 1.0 everywhere."""
    jm = _glob_max(2024, 4, 8)
    # Dallas is inside the path of totality.
    lon, lat = -96.80, 32.78
    tmax = _local_max(jm, lon, lat)
    geo = (lon, lat, 0.0)

    _rc, attr = sol_eclipse_how(tmax, geo, FLG_SWIEPH)
    assert attr[2] == pytest.approx(1.0, abs=1e-9)
    assert sol_eclipse_obscuration_at_loc(tmax, geo, FLG_SWIEPH) == pytest.approx(
        1.0, abs=1e-9
    )
    details = sol_eclipse_how_details(tmax, geo, FLG_SWIEPH)
    assert details["max_obscuration"] == pytest.approx(1.0, abs=1e-9)
    assert details["max_obscuration_percent"] == pytest.approx(100.0, abs=1e-6)
    # The eclipse magnitude (diameter coverage) still records the >1 excess.
    assert attr[0] > 1.0


def test_obscuration_never_exceeds_one():
    """Across the path and its surroundings, obscuration stays a fraction <= 1."""
    jm = _glob_max(2024, 4, 8)
    for lat in range(-10, 71, 8):
        _rc, attr = sol_eclipse_how(jm, (-100.0, float(lat), 0.0), FLG_SWIEPH)
        assert 0.0 <= attr[2] <= 1.0 + 1e-9, f"obscuration {attr[2]} at lat {lat}"


def test_annular_obscuration_is_area_ratio_below_one():
    """Annular: a ring of Sun remains, so obscuration = (ratio)^2 < 1."""
    jm = _glob_max(2023, 10, 14)
    _wrf, geopos, _a = sol_eclipse_where(jm, FLG_SWIEPH)
    geo = (geopos[0], geopos[1], 0.0)
    _rc, attr = sol_eclipse_how(jm, geo, FLG_SWIEPH)
    obsc = attr[2]
    assert 0.0 < obsc < 1.0
    # attr[1] is the lunar/solar diameter ratio; area ratio is its square.
    assert obsc == pytest.approx(attr[1] ** 2, rel=1e-6)


def test_partial_obscuration_is_fraction():
    """Partial: 0 < obscuration < 1 (sampled away from the central path)."""
    jm = _glob_max(2024, 4, 8)
    for lon, lat in [(-74.0, 40.7), (-87.6, 41.9)]:
        tmax = _local_max(jm, lon, lat)
        obsc = sol_eclipse_obscuration_at_loc(tmax, (lon, lat, 0.0), FLG_SWIEPH)
        assert 0.0 < obsc < 1.0
