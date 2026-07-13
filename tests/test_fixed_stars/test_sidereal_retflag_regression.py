"""Regression tests for fixed-star sidereal reduction and echoed flags.

* I2-F1 -- sidereal fixed-star positions are referred to the mean equinox
  of date: the ecliptic longitude uses the TRUE ayanamsha (mean + nutation in
  longitude) and the equatorial output is measured on the MEAN equator with the
  mean ayanamsha removed from right ascension.

* I2-F2 -- the UT entry points (fixstar_ut / fixstar2_ut) echo implied flags:
  FLG_NONUT for J2000/SIDEREAL and
  FLG_NOGDEFL | FLG_NOABERR for heliocentric / true-position output.

The tests check internal consistency against the library's independently
implemented nutation and ayanamsha models.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem


def _wrap_deg(x: float) -> float:
    """Wrap an angular difference (degrees) to (-180, 180]."""
    return ((x + 180.0) % 360.0) - 180.0


def _set_independent_sidereal_mode() -> None:
    """Select an arbitrary caller-defined zero point for frame invariants."""
    ephem.set_sid_mode(ephem.SIDM_USER, 2451545.0, 24.0)


# ---------------------------------------------------------------------------
# I2-F1: sidereal reduction referred to the mean equinox
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestSiderealStarMeanEquinox:
    STARS = ["Regulus", "Spica"]
    DATES = [
        (2000, 1, 1),
        (2050, 1, 1),
    ]

    def test_ecliptic_sidereal_subtracts_true_ayanamsha(self):
        """(trop - mean_ayanamsha - sid) must equal the nutation in longitude.

        Sidereal longitude = mean_lon - mean_ayanamsha = true_lon - (mean +
        dpsi) ayanamsha, so subtracting the *mean* ayanamsha from the true
        longitude leaves exactly dpsi. dpsi is read reference-free as the
        of-date minus mean-ecliptic (NONUT) longitude.
        """
        _set_independent_sidereal_mode()
        saw_nonzero = False
        for star in self.STARS:
            for y, m, d in self.DATES:
                jd = ephem.julday(y, m, d, 12.0)
                trop = ephem.fixstar2_ut(star, jd, 0)[0][0]
                nonut = ephem.fixstar2_ut(star, jd, ephem.FLG_NONUT)[0][0]
                sid = ephem.fixstar2_ut(star, jd, ephem.FLG_SIDEREAL)[0][0]
                mean_aya = ephem.get_ayanamsa_ut(jd)

                dpsi_arcsec = _wrap_deg(trop - nonut) * 3600.0
                resid_arcsec = _wrap_deg(trop - mean_aya - sid) * 3600.0

                if abs(dpsi_arcsec) > 1.0:
                    saw_nonzero = True
                assert abs(resid_arcsec - dpsi_arcsec) < 0.005, (
                    f'{star} {y}: residual {resid_arcsec:.4f}" should equal '
                    f'dpsi {dpsi_arcsec:.4f}" (true ayanamsha subtracted)'
                )
        assert saw_nonzero, "test dates should have non-trivial nutation"

    def test_equatorial_sidereal_uses_mean_equator(self):
        """SID|EQ right ascension = mean-equator RA - mean ayanamsha.

        Declination must match the mean equator (NONUT), unchanged by the
        ayanamsha subtraction. The pre-fix code used the true equator, leaving
        a nutation-in-RA error and a shifted declination.
        """
        _set_independent_sidereal_mode()
        for star in self.STARS:
            for y, m, d in [(2000, 1, 1), (1950, 1, 1)]:
                jd = ephem.julday(y, m, d, 12.0)
                mean_eq = ephem.fixstar2_ut(
                    star, jd, ephem.FLG_EQUATORIAL | ephem.FLG_NONUT
                )[0]
                sid_eq = ephem.fixstar2_ut(
                    star, jd, ephem.FLG_EQUATORIAL | ephem.FLG_SIDEREAL
                )[0]
                mean_aya = ephem.get_ayanamsa_ut(jd)

                resid_ra = _wrap_deg(mean_eq[0] - mean_aya - sid_eq[0]) * 3600.0
                resid_dec = (mean_eq[1] - sid_eq[1]) * 3600.0
                assert abs(resid_ra) < 0.01, (
                    f'{star} {y}: SID|EQ RA residual {resid_ra:.4f}"'
                )
                assert abs(resid_dec) < 0.01, (
                    f'{star} {y}: SID|EQ Dec residual {resid_dec:.4f}"'
                )

    def test_j2000_and_nonut_sidereal_use_mean_ayanamsha(self):
        """SID|J2000 and SID|NONUT subtract the mean ayanamsha (no dpsi)."""
        _set_independent_sidereal_mode()
        star = "Regulus"
        jd = ephem.julday(2000, 1, 1, 12.0)
        mean_aya = ephem.get_ayanamsa_ut(jd)
        for frame in (ephem.FLG_J2000, ephem.FLG_NONUT):
            trop = ephem.fixstar2_ut(star, jd, frame)[0][0]
            sid = ephem.fixstar2_ut(star, jd, frame | ephem.FLG_SIDEREAL)[0][0]
            resid = _wrap_deg(trop - mean_aya - sid) * 3600.0
            assert abs(resid) < 0.005, f'frame {frame}: residual {resid:.4f}"'

    def test_sidereal_longitude_speed_excludes_nutation_rate(self):
        """Sidereal SPEED must not carry the of-date nutation rate.

        Since sidereal longitude = mean_lon - mean_ayanamsha, its rate equals
        the mean-ecliptic (NONUT) longitude rate minus the mean-ayanamsha
        drift, with no d(dpsi)/dt term. Comparing against the NONUT speed and a
        finite-difference ayanamsha rate is independently computed; the pre-fix code
        differenced the true (of-date) longitude and so retained d(dpsi)/dt,
        ~0.16"/day at 1900/2200.
        """
        _set_independent_sidereal_mode()
        h = 1.0
        for star in self.STARS:
            for y in (1900, 2000, 2200):
                jd = ephem.julday(y, 1, 1, 12.0)
                sid_sp = ephem.fixstar2_ut(
                    star, jd, ephem.FLG_SPEED | ephem.FLG_SIDEREAL
                )[0][3]
                nonut_sp = ephem.fixstar2_ut(
                    star, jd, ephem.FLG_SPEED | ephem.FLG_NONUT
                )[0][3]
                aya_rate = _wrap_deg(
                    ephem.get_ayanamsa_ut(jd + h) - ephem.get_ayanamsa_ut(jd - h)
                ) / (2.0 * h)
                resid = (sid_sp - (nonut_sp - aya_rate)) * 3600.0
                assert abs(resid) < 0.01, (
                    f'{star} {y}: sid-speed residual {resid:.5f}"/day'
                )


# ---------------------------------------------------------------------------
# I2-F2: UT entry points echo the implied flag bits
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestFixstarUtRetflagImpliedBits:
    _JD = 2451545.0  # J2000.0

    # (input flags, expected echoed retflag):
    #   * default ephemeris bit FLG_SWIEPH added when none given
    #   * FLG_NONUT echoed for J2000 / SIDEREAL (mean-equinox output)
    #   * FLG_NOGDEFL | FLG_NOABERR echoed for HELCTR / TRUEPOS
    CASES = [
        (0, ephem.FLG_SWIEPH),
        (ephem.FLG_NONUT, ephem.FLG_NONUT | ephem.FLG_SWIEPH),
        (ephem.FLG_J2000, ephem.FLG_J2000 | ephem.FLG_SWIEPH | ephem.FLG_NONUT),
        (
            ephem.FLG_SIDEREAL,
            ephem.FLG_SIDEREAL | ephem.FLG_SWIEPH | ephem.FLG_NONUT,
        ),
        (
            ephem.FLG_SIDEREAL | ephem.FLG_EQUATORIAL,
            ephem.FLG_SIDEREAL
            | ephem.FLG_EQUATORIAL
            | ephem.FLG_SWIEPH
            | ephem.FLG_NONUT,
        ),
        (
            ephem.FLG_HELCTR,
            ephem.FLG_HELCTR | ephem.FLG_SWIEPH | ephem.FLG_NOGDEFL | ephem.FLG_NOABERR,
        ),
        (
            ephem.FLG_TRUEPOS,
            ephem.FLG_TRUEPOS
            | ephem.FLG_SWIEPH
            | ephem.FLG_NOGDEFL
            | ephem.FLG_NOABERR,
        ),
    ]

    @pytest.mark.parametrize("func_name", ["fixstar_ut", "fixstar2_ut"])
    @pytest.mark.parametrize("flags_in,expected", CASES)
    def test_ut_retflag_echoes_implied_bits(
        self, func_name: str, flags_in: int, expected: int
    ):
        _set_independent_sidereal_mode()
        func = getattr(ephem, func_name)
        retflag = func("Regulus", self._JD, flags_in)[2]
        assert retflag == expected, (
            f"{func_name}(flags={flags_in}) echoed {retflag}, expected {expected}"
        )
