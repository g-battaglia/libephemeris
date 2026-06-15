"""Characterization: solar-eclipse path geometry vs NASA's published canon.

Independent reference: the greatest-eclipse position and path widths published
in NASA's Five Millennium Canon of Solar Eclipses
(https://eclipse.gsfc.nasa.gov/solar.html).

  * Greatest-eclipse central position (sol_eclipse_where) matches the canon to
    ~0.05 deg -- the shadow-axis / ellipsoid intersection is accurate.
  * The northern/southern limits bracket the central line at every traced step.
  * calc_eclipse_path_width is a simplified geometric model (cone radius l2
    projected via 1/sin(altitude)); near greatest eclipse it runs ~5-8% wider
    than the canon's perpendicular-to-track path widths. This test pins that
    band so the agreement level is documented and gross regressions are caught.
"""

from __future__ import annotations

import pytest

from libephemeris import (
    FLG_SWIEPH,
    calc_eclipse_central_line,
    calc_eclipse_northern_limit,
    calc_eclipse_path_width,
    calc_eclipse_southern_limit,
    julday,
    sol_eclipse_when_glob,
    sol_eclipse_where,
)

pytestmark = pytest.mark.slow

# (year, month, day, NASA greatest-eclipse lat, lon (E<0 = W), path width km)
_CANON = [
    (2017, 8, 21, 36.967, -87.672, 114.7),
    (2024, 4, 8, 25.292, -104.138, 197.5),
]


def _greatest(year, month, day):
    jd0 = julday(year, month, day, 0.0)
    _rf, tret = sol_eclipse_when_glob(jd0, FLG_SWIEPH, 0, False)
    return tret[0]


@pytest.mark.parametrize("year,month,day,lat,lon,_w", _CANON)
def test_greatest_eclipse_position_matches_canon(year, month, day, lat, lon, _w):
    jmax = _greatest(year, month, day)
    _rf, geopos, _attr = sol_eclipse_where(jmax, FLG_SWIEPH)
    assert geopos[1] == pytest.approx(lat, abs=0.05), f"lat {geopos[1]} vs {lat}"
    assert geopos[0] == pytest.approx(lon, abs=0.05), f"lon {geopos[0]} vs {lon}"


@pytest.mark.parametrize("year,month,day,_lat,_lon,width", _CANON)
def test_path_width_within_documented_band(year, month, day, _lat, _lon, width):
    jmax = _greatest(year, month, day)
    w = calc_eclipse_path_width(jmax, flags=FLG_SWIEPH)
    ratio = w / width
    # Simplified model runs ~5-8% high near greatest eclipse (see docstring).
    assert 0.95 <= ratio <= 1.15, f"{year}: path width {w:.1f} km, canon {width} (ratio {ratio:.3f})"


@pytest.mark.parametrize("year,month,day,_lat,_lon,_w", _CANON)
def test_limits_bracket_central_line(year, month, day, _lat, _lon, _w):
    jmax = _greatest(year, month, day)
    a, b = jmax - 1.0 / 24.0, jmax + 1.0 / 24.0
    _tc, latc, lonc = calc_eclipse_central_line(a, b, 20.0, FLG_SWIEPH)
    _tn, latn, _lonn = calc_eclipse_northern_limit(a, b, 20.0, FLG_SWIEPH)
    _ts, lats, _lons = calc_eclipse_southern_limit(a, b, 20.0, FLG_SWIEPH)
    assert len(latc) >= 3
    assert len(latn) == len(latc) and len(lats) == len(latc)
    for s, c, n in zip(lats, latc, latn):
        # Northern limit is north of the central line, southern limit south.
        assert s < c < n, f"limits do not bracket central line: {s} < {c} < {n}"
        assert (n - s) < 5.0  # plausible band extent in latitude (deg)
