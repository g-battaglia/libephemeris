"""Regression: solar-eclipse obscuration semantics per entry point.

For a TOTAL eclipse the reference API reports the Moon/Sun disc *area ratio*
``(R_moon/R_sun)^2 ~ 1.05-1.12 > 1`` in ``attr[2]``: it is the eclipse
magnitude's domain rather than a physically-bounded covered fraction. The
reference-compatible entry points (``sol_eclipse_how``,
``sol_eclipse_when_loc``, ``sol_eclipse_where``) mirror that behavior 1:1
(review round 2 removed an earlier clamp to 1.0 that diverged from the
reference; e.g. Dallas 2024-04-08 must report ~1.1167, verified against the
reference ephemeris 2.10.03).

The convenience extensions (``sol_eclipse_obscuration_at_loc``,
``sol_eclipse_how_details``) are not part of the reference API and keep
reporting the physically-bounded covered fraction (1.0 when the Sun is fully
covered).

Annular eclipses report the same area ratio in both conventions (< 1, a ring
of Sun remains) and partial eclipses use the two-disc lens overlap.

Checks (independent of any external implementation at runtime):
  * total: attr[2] equals the squared diameter ratio (> 1) from the
    reference-compatible entry points, 1.0 from the extensions,
  * attr[2] equals the squared diameter ratio wherever the eclipse is total,
    and stays a fraction <= 1.0 elsewhere,
  * annular: 0 < obscuration < 1 and equals the squared diameter ratio,
  * partial: 0 < obscuration < 1.
"""

from __future__ import annotations

import pytest

from libephemeris import (
    ECL_ANNULAR,
    ECL_TOTAL,
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

# Reference obscuration at the Dallas local maximum of the 2024-04-08 total
# eclipse (reference ephemeris 2.10.03: 1.1166797).
DALLAS_OBSCURATION_REF = 1.11668


def _glob_max(year, month, day):
    jd0 = julday(year, month, day, 0.0)
    _rf, tret = sol_eclipse_when_glob(jd0, FLG_SWIEPH, 0, False)
    return tret[0]


def _local_max(jd_near, lon, lat):
    _rc, tloc, _attr = sol_eclipse_when_loc(
        jd_near - 0.2, (lon, lat, 0.0), FLG_SWIEPH, False
    )
    return tloc[0]


def test_total_eclipse_obscuration_is_area_ratio_all_apis():
    """Total: attr[2] is the disc area ratio > 1 (reference behavior);
    the non-reference extensions keep the covered fraction 1.0."""
    jm = _glob_max(2024, 4, 8)
    # Dallas is inside the path of totality.
    lon, lat = -96.80, 32.78
    tmax = _local_max(jm, lon, lat)
    geo = (lon, lat, 0.0)

    rc, attr = sol_eclipse_how(tmax, geo, FLG_SWIEPH)
    assert rc & ECL_TOTAL
    # Reference-compatible attr[2]: the Moon/Sun disc area ratio, > 1.
    assert attr[2] == pytest.approx(attr[1] ** 2, rel=1e-9)
    assert attr[2] == pytest.approx(DALLAS_OBSCURATION_REF, abs=2e-3)
    # The eclipse magnitude (diameter coverage) also records the >1 excess.
    assert attr[0] > 1.0

    # when_loc reports the same attr layout at its own local maximum.
    _rc_loc, _tloc, attr_loc = sol_eclipse_when_loc(jm - 0.2, geo, FLG_SWIEPH, False)
    assert attr_loc[2] == pytest.approx(DALLAS_OBSCURATION_REF, abs=2e-3)

    # Extensions (not in the reference API): physically-bounded fraction.
    assert sol_eclipse_obscuration_at_loc(tmax, geo, FLG_SWIEPH) == pytest.approx(
        1.0, abs=1e-9
    )
    details = sol_eclipse_how_details(tmax, geo, FLG_SWIEPH)
    assert details["max_obscuration"] == pytest.approx(1.0, abs=1e-9)
    assert details["max_obscuration_percent"] == pytest.approx(100.0, abs=1e-6)


def test_obscuration_is_area_ratio_only_in_totality():
    """Across the path and its surroundings: area ratio (>1) exactly where
    the eclipse is total, a fraction <= 1 everywhere else."""
    jm = _glob_max(2024, 4, 8)
    for lat in range(-10, 71, 8):
        rc, attr = sol_eclipse_how(jm, (-100.0, float(lat), 0.0), FLG_SWIEPH)
        if rc & ECL_TOTAL:
            assert attr[2] == pytest.approx(attr[1] ** 2, rel=1e-9), (
                f"total-phase obscuration at lat {lat} must be the area ratio"
            )
            assert attr[2] > 1.0
        else:
            assert 0.0 <= attr[2] <= 1.0 + 1e-9, f"obscuration {attr[2]} at lat {lat}"


def test_annular_obscuration_is_area_ratio_below_one():
    """Annular: a ring of Sun remains, so obscuration = (ratio)^2 < 1."""
    jm = _glob_max(2023, 10, 14)
    _wrf, geopos, _a = sol_eclipse_where(jm, FLG_SWIEPH)
    geo = (geopos[0], geopos[1], 0.0)
    rc, attr = sol_eclipse_how(jm, geo, FLG_SWIEPH)
    assert rc & ECL_ANNULAR
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
