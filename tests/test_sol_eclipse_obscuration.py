"""Solar-eclipse obscuration invariants per entry point.

The reference-facing ``attr[2]`` follows the measured reference behavior:
while one disc lies inside the other it is the disc-area ratio
``(R_moon/R_sun)^2`` — greater than 1 during totality (the Moon overfills
the Sun) and less than 1 for an annular eclipse (a ring of Sun remains) —
and the two-disc lens-overlap fraction for a partial phase. The
library-specific extensions (``sol_eclipse_obscuration_at_loc`` and the
``max_obscuration`` field of ``sol_eclipse_how_details``) instead report
the physically bounded covered fraction in [0, 1].

Checks (independent of any external implementation at runtime):
  * total: attr[2] == (R_moon/R_sun)^2 > 1 from the compatibility entry
    points; the bounded extensions report 1.0,
  * attr[2] == the area ratio where the eclipse is total,
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


def _glob_max(year, month, day):
    jd0 = julday(year, month, day, 0.0)
    _rf, tret = sol_eclipse_when_glob(jd0, FLG_SWIEPH, 0, False)
    return tret[0]


def _local_max(jd_near, lon, lat):
    _rc, tloc, _attr = sol_eclipse_when_loc(
        jd_near - 0.2, (lon, lat, 0.0), FLG_SWIEPH, False
    )
    return tloc[0]


def test_total_eclipse_obscuration_is_bounded_all_apis():
    """Total: attr[2] is the disc-area ratio (> 1) from the reference-facing
    entry points; the bounded extensions report the covered fraction 1.0."""
    jm = _glob_max(2024, 4, 8)
    # Dallas is inside the path of totality.
    lon, lat = -96.80, 32.78
    tmax = _local_max(jm, lon, lat)
    geo = (lon, lat, 0.0)

    rc, attr = sol_eclipse_how(tmax, geo, FLG_SWIEPH)
    assert rc & ECL_TOTAL
    # attr[2]: the reference reports the disc-area ratio, which exceeds 1 in
    # totality and equals the squared diameter ratio.
    assert attr[2] == pytest.approx(attr[1] ** 2, rel=1e-9)
    assert attr[2] > 1.0
    # The eclipse magnitude (diameter coverage) also records the > 1 excess.
    assert attr[0] > 1.0

    # when_loc reports the same attr layout at its own local maximum.
    _rc_loc, _tloc, attr_loc = sol_eclipse_when_loc(jm - 0.2, geo, FLG_SWIEPH, False)
    assert attr_loc[2] == pytest.approx(attr_loc[1] ** 2, rel=1e-9)
    assert attr_loc[2] > 1.0

    # Library-specific extensions expose a physically bounded fraction.
    assert sol_eclipse_obscuration_at_loc(tmax, geo, FLG_SWIEPH) == pytest.approx(
        1.0, abs=1e-9
    )
    details = sol_eclipse_how_details(tmax, geo, FLG_SWIEPH)
    assert details["max_obscuration"] == pytest.approx(1.0, abs=1e-9)
    assert details["max_obscuration_percent"] == pytest.approx(100.0, abs=1e-6)


def test_obscuration_is_area_ratio_in_totality():
    """Across the path and its surroundings: the disc-area ratio where the
    eclipse is total (>= 1), and a fraction <= 1 everywhere else."""
    jm = _glob_max(2024, 4, 8)
    for lat in range(-10, 71, 8):
        rc, attr = sol_eclipse_how(jm, (-100.0, float(lat), 0.0), FLG_SWIEPH)
        if rc & ECL_TOTAL:
            assert attr[2] == pytest.approx(attr[1] ** 2, rel=1e-6), (
                f"total-phase obscuration at lat {lat} must be the area ratio"
            )
            assert attr[2] >= 1.0 - 1e-9, f"obscuration {attr[2]} at lat {lat}"
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
