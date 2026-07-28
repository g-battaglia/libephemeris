"""ascmc angle-velocity stencil in houses_ex2 (finding AL-F2).

The angle rates (Asc, Vertex, co-ascendant, ...) carry a large third
derivative on the latitude-amplified slots, so a fixed centered difference at
1/4096 day leaves an f'''*h**2 truncation error that grows from a few
arcsec/day at mid-latitude to hundreds of arcsec/day near the polar circle.
houses_ex2 removes it with two-step Richardson extrapolation of the angle
positions at h = 1/4096 and h/2 = 1/8192 day (O(h**2) -> O(h**4)); the cusp
rates keep the plain centered difference.

These tests are self-contained: the "reference" derivative here is an
independent Richardson estimate at a different step pair (1/2048, 1/4096).
Two consistent extrapolations agreeing to a fraction of an arcsec/day is the
in-repo proxy for the external-reference agreement measured separately.
"""

from __future__ import annotations

import pytest

import libephemeris as swe

LON = 12.4964  # Rome
_H = 1.0 / 4096.0  # implementation's coarse half-step


def _wrap(d: float) -> float:
    return (d + 180.0) % 360.0 - 180.0


def _central_ascmc(jd: float, lat: float, hsys: int, slot: int, h: float) -> float:
    """Centered difference of an ascmc angle position, deg/day."""
    _, am = swe.houses_ex(jd - h, lat, LON, hsys, 0)
    _, ap = swe.houses_ex(jd + h, lat, LON, hsys, 0)
    return _wrap(ap[slot] - am[slot]) / (2.0 * h)


def _richardson_ascmc(
    jd: float, lat: float, hsys: int, slot: int, big_h: float
) -> float:
    """Two-step Richardson estimate combining central diffs at big_h, big_h/2."""
    d_big = _central_ascmc(jd, lat, hsys, slot, big_h)
    d_small = _central_ascmc(jd, lat, hsys, slot, big_h / 2.0)
    return (4.0 * d_small - d_big) / 3.0


class TestAscmcSpeedRichardson:
    @pytest.mark.unit
    def test_reported_equals_richardson_stencil(self):
        """Every ascmc rate is the documented (h, h/2) Richardson combination."""
        jd = swe.julday(2000, 6, 15, 9.0)
        lat = 41.9028
        _, _, _, ascmc_speed = swe.houses_ex2(jd, lat, LON, ord("P"), 0)
        for slot in range(8):
            expected = _richardson_ascmc(jd, lat, ord("P"), slot, _H)
            assert abs(ascmc_speed[slot] - expected) < 1e-9, (
                f"slot {slot}: reported {ascmc_speed[slot]} != Richardson {expected}"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("hour", [0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
    def test_midlat_angles_converged(self, hour):
        """At Rome every angle rate agrees with an independent Richardson pair
        to well under 0.25 arcsec/day across the day."""
        jd = swe.julday(1948, 3, 7, hour)
        lat = 41.9028
        _, _, _, ascmc_speed = swe.houses_ex2(jd, lat, LON, ord("P"), 0)
        for slot in range(8):
            indep = _richardson_ascmc(jd, lat, ord("P"), slot, 2.0 * _H)
            err_arcsec = abs(ascmc_speed[slot] - indep) * 3600.0
            assert err_arcsec < 0.25, (
                f'hour {hour} slot {slot}: cross-pair {err_arcsec:.4f}"/day >= 0.25'
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("lat", [64.0, 66.0])
    def test_high_latitude_ascendant_fix(self, lat):
        """Near the polar circle the Ascendant rate is (a) a converged value,
        and (b) far from the plain 1/4096 centered difference the pre-fix code
        reported -- confirming Richardson materially corrects the truncation."""
        jd = swe.julday(1948, 3, 7, 6.6)
        _, _, _, ascmc_speed = swe.houses_ex2(jd, lat, LON, ord("P"), 0)
        reported = ascmc_speed[0]  # Ascendant

        indep = _richardson_ascmc(jd, lat, ord("P"), 0, 2.0 * _H)
        converged_arcsec = abs(reported - indep) * 3600.0
        assert converged_arcsec < 2.0, (
            f'lat {lat}: Asc rate not converged, cross-pair {converged_arcsec:.3f}"/day'
        )

        naive = _central_ascmc(jd, lat, ord("P"), 0, _H)
        correction_arcsec = abs(reported - naive) * 3600.0
        assert correction_arcsec > 50.0, (
            f'lat {lat}: Richardson correction only {correction_arcsec:.1f}"/day; '
            "the plain 1/4096 centered difference would fail here"
        )

    @pytest.mark.unit
    def test_slow_angles_unperturbed(self):
        """MC/ARMC move quasi-uniformly and were already near the roundoff
        floor; Richardson must not perturb them beyond a small fraction of an
        arcsec/day, and ARMC keeps the sidereal rotation rate."""
        jd = swe.julday(1948, 3, 7, 6.6)
        lat = 41.9028
        _, _, _, ascmc_speed = swe.houses_ex2(jd, lat, LON, ord("P"), 0)
        # ARMC advances at ~360.986 deg/day.
        assert abs(ascmc_speed[2] - 360.9856) < 0.01
        for slot in (1, 2):  # MC, ARMC
            naive = _central_ascmc(jd, lat, ord("P"), slot, _H)
            drift_arcsec = abs(ascmc_speed[slot] - naive) * 3600.0
            assert drift_arcsec < 0.5, (
                f'slow slot {slot}: Richardson shifted the rate by {drift_arcsec:.4f}"/day'
            )
