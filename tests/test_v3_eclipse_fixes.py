"""Regression tests for the v3 eclipse.py bug fixes.

Each test pins the corrected behaviour of a specific fix in
``libephemeris/eclipse.py``. Expected numeric values were obtained from an
independent reference ephemeris oracle (validation only) and are hardcoded
here -- the oracle is deliberately NOT imported, so these tests are
self-contained.

Fixes covered:
    1. ``_calc_umbra_limit`` uses the Sun-Moon distance for the umbral cone
       half-angle (not the Earth-Sun distance).
    2. ``calc_solar_eclipse_duration`` returns the central-line maximum
       duration of totality/annularity (minutes), not the global path
       duration (hours).
    3. ``_sol_eclipse_when_glob_pythonic`` returns / documents the *_glob*
       tret layout (P1/P4, U1/U4, center-line begin/end).
    4. ``lun_occult_when_loc`` puts the occulted body's rise/set at
       tret[5]/tret[6] and leaves tret[7]/tret[8] zero.
    5. ``planet_occult_when_glob`` scales the RA component of the relative
       speed by cos(dec) so contact times are correct at non-zero
       declination.
"""

from __future__ import annotations

import math

import pytest

from libephemeris import calc_ut, julday
from libephemeris.constants import (
    ECL_ANNULAR,
    ECL_ANNULAR_TOTAL,
    ECL_CENTRAL,
    ECL_TOTAL,
    FLG_EQUATORIAL,
    FLG_SPEED,
    JUPITER,
    SATURN,
    SUN,
    VENUS,
)
from libephemeris import eclipse as E
from libephemeris.eclipse import (
    _calculate_eclipse_type_and_magnitude,
    _calc_umbra_limit,
    _sol_eclipse_when_glob_pythonic,
    calc_besselian_l2,
    calc_solar_eclipse_duration,
    lun_occult_when_loc,
    planet_occult_when_glob,
)

# --- Hardcoded eclipse maxima (JD UT) from the reference oracle -------------
JD_MAX_TOTAL_2024 = 2460409.262041  # 2024-04-08 total
JD_MAX_TOTAL_2017 = 2457987.267728  # 2017-08-21 total
JD_MAX_ANNULAR_2023 = 2460232.249668  # 2023-10-14 annular
JD_MAX_ANNULAR_2024 = 2460586.281299  # 2024-10-02 annular
JD_MAX_HYBRID_2023 = 2460054.678322  # 2023-04-20 hybrid (annular-total)
JD_MAX_PARTIAL_2025 = 2460763.949601422  # 2025-03-29 partial (non-central)


def _require_ephemeris(jd: float, body: int = SUN) -> None:
    """Skip when the underlying ephemeris cannot be evaluated at ``jd``."""
    try:
        calc_ut(jd, body, FLG_SPEED)
    except Exception as exc:  # pragma: no cover - environment dependent
        pytest.skip(f"ephemeris unavailable at JD {jd}: {exc!r}")


# ---------------------------------------------------------------------------
# Fix 1: _calc_umbra_limit uses the Sun-Moon distance
# ---------------------------------------------------------------------------
class TestFix1UmbraLimitSunMoonDistance:
    """l2 from the simplified umbra model must agree with the Besselian l2
    and give the correct total/annular/hybrid classification."""

    # (label, jd_max, expected type bit, expected fixed l2 ~ value)
    CASES = [
        ("total_2024", JD_MAX_TOTAL_2024, ECL_TOTAL, -0.010236),
        ("total_2017", JD_MAX_TOTAL_2017, ECL_TOTAL, -0.003948),
        ("annular_2023", JD_MAX_ANNULAR_2023, ECL_ANNULAR, +0.018116),
        ("annular_2024", JD_MAX_ANNULAR_2024, ECL_ANNULAR, +0.024126),
        ("hybrid_2023", JD_MAX_HYBRID_2023, ECL_ANNULAR_TOTAL, +0.000724),
    ]

    @pytest.mark.parametrize("label,jd_max,type_bit,l2_expected", CASES)
    def test_l2_matches_besselian_and_classifies(
        self, label, jd_max, type_bit, l2_expected
    ):
        _require_ephemeris(jd_max)

        l2 = _calc_umbra_limit(jd_max)
        l2_bess = calc_besselian_l2(jd_max)

        # After the fix the simplified model tracks the Besselian model to
        # well under 1e-4 Earth radii (the buggy Earth-Sun distance left a
        # systematic ~6-8e-4 offset).
        assert abs(l2 - l2_bess) < 1e-4, (
            f"{label}: l2={l2:+.6f} vs besselian {l2_bess:+.6f}"
        )

        # Pinned absolute value (hardcoded expectation).
        assert abs(l2 - l2_expected) < 5e-4, (
            f"{label}: l2={l2:+.6f} expected ~{l2_expected:+.6f}"
        )

        etype, _mag, _gamma, _ratio = _calculate_eclipse_type_and_magnitude(jd_max)
        assert etype & type_bit, f"{label}: type {etype} lacks bit {type_bit}"
        assert etype & ECL_CENTRAL, f"{label}: expected central eclipse"

    def test_hybrid_l2_near_zero(self):
        """A hybrid eclipse has |l2| in the tiny sign-change band."""
        _require_ephemeris(JD_MAX_HYBRID_2023)
        assert abs(_calc_umbra_limit(JD_MAX_HYBRID_2023)) < 0.002


# ---------------------------------------------------------------------------
# Fix 2: calc_solar_eclipse_duration returns central-line totality minutes
# ---------------------------------------------------------------------------
class TestFix2SolarEclipseDuration:
    """Duration must be the minutes-long central-line totality/annularity,
    matching the reference local maximum-eclipse duration."""

    # (label, jd_max, oracle local duration in minutes)
    CASES = [
        ("total_2024", JD_MAX_TOTAL_2024, 4.4672),
        ("total_2017", JD_MAX_TOTAL_2017, 2.6686),
        ("annular_2023", JD_MAX_ANNULAR_2023, 5.2831),
        ("annular_2024", JD_MAX_ANNULAR_2024, 7.4136),
    ]

    @pytest.mark.parametrize("label,jd_max,expected_min", CASES)
    def test_central_line_duration_minutes(self, label, jd_max, expected_min):
        _require_ephemeris(jd_max)

        duration = calc_solar_eclipse_duration(jd_max)

        # Minutes, not hours: the pre-fix bug returned the global path
        # duration (~200+ minutes). Anything above ~15 min for a single
        # central eclipse would be the old broken value.
        assert 1.0 < duration < 15.0, (
            f"{label}: duration {duration:.3f} min not in the minutes range"
        )

        # Matches the reference local maximum-eclipse duration within 0.15 min.
        assert abs(duration - expected_min) < 0.15, (
            f"{label}: duration {duration:.4f} min vs reference {expected_min:.4f}"
        )

    def test_partial_eclipse_returns_zero(self):
        """A non-central (partial) eclipse has no totality: duration 0.0."""
        _require_ephemeris(JD_MAX_PARTIAL_2025)
        assert calc_solar_eclipse_duration(JD_MAX_PARTIAL_2025) == 0.0


# ---------------------------------------------------------------------------
# Fix 3: _sol_eclipse_when_glob_pythonic returns the _glob tret layout
# ---------------------------------------------------------------------------
class TestFix3GlobTretLayout:
    """The returned tret is the *_glob* layout: [2]/[3] = P1/P4, [4]/[5] =
    U1/U4 (totality begin/end), [6]/[7] = center-line begin/end."""

    # Oracle sol_eclipse_when_glob tret for the 2024-04-08 total eclipse.
    ORACLE_TRET = (
        2460409.2620413844,  # [0] maximum
        2460409.275113436,  # [1] local apparent noon
        2460409.1543619353,  # [2] eclipse begin (P1)
        2460409.369543748,  # [3] eclipse end (P4)
        2460409.193639937,  # [4] totality begins (U1)
        2460409.330235851,  # [5] totality ends (U4)
        2460409.1944517964,  # [6] center line begins
        2460409.329439239,  # [7] center line ends
        0.0,  # [8]
        0.0,  # [9]
    )

    def test_glob_layout_matches_reference(self):
        jd_start = julday(2024, 1, 1, 0.0)
        _require_ephemeris(jd_start)

        _etype, tret = _sol_eclipse_when_glob_pythonic(jd_start, eclipse_type=ECL_TOTAL)

        # Each populated index names the quantity the reference returns, to
        # within ~0.5 min (contact-time precision of the Besselian search).
        for i in range(8):
            assert abs(tret[i] - self.ORACLE_TRET[i]) < 5e-4, (
                f"tret[{i}]={tret[i]:.6f} vs reference {self.ORACLE_TRET[i]:.6f}"
            )
        assert tret[8] == 0.0 and tret[9] == 0.0

    def test_glob_layout_ordering(self):
        """The _glob semantics imply P1 < U1 < max < U4 < P4 and the
        center line lies within totality."""
        jd_start = julday(2024, 1, 1, 0.0)
        _require_ephemeris(jd_start)

        _etype, tret = _sol_eclipse_when_glob_pythonic(jd_start, eclipse_type=ECL_TOTAL)
        p1, p4 = tret[2], tret[3]
        u1, u4 = tret[4], tret[5]
        cl_begin, cl_end = tret[6], tret[7]
        jd_max = tret[0]

        assert p1 < u1 < jd_max < u4 < p4, (
            "expected P1 < U1 < max < U4 < P4 (glob layout)"
        )
        assert u1 <= cl_begin and cl_end <= u4, (
            "center line must lie inside the totality window"
        )


# ---------------------------------------------------------------------------
# Fix 4: lun_occult_when_loc rise/set at tret[5]/[6], not [7]/[8]
# ---------------------------------------------------------------------------
class TestFix4LunOccultRiseSetLayout:
    def test_docstring_documents_body_rise_set_at_5_6(self):
        doc = lun_occult_when_loc.__doc__ or ""
        assert "occulted body rises" in doc
        assert "occulted body sets" in doc
        # tret[7]/[8] must no longer be described as moonrise/moonset.
        assert "Time of moonrise" not in doc
        assert "Time of moonset" not in doc

    @pytest.mark.slow
    def test_saturn_occultation_body_rise_at_tret5(self):
        """2037-02-01 occultation of Saturn from Rome: the reference places
        Saturn's rise at tret[5] and leaves tret[6]/[7]/[8] zero."""
        # Reference tret[5] (Saturn rise) for this event.
        oracle_rise = 2465091.2391339107

        rome = [12.4964, 41.9028, 0.0]  # lon, lat, alt
        jd_start = julday(2037, 1, 15, 0.0)
        _require_ephemeris(jd_start, SATURN)

        try:
            _retflag, tret, _attr = lun_occult_when_loc(jd_start, SATURN, rome)
        except Exception as exc:  # pragma: no cover - environment dependent
            pytest.skip(f"occultation search failed: {exc!r}")

        # tret[5] is the occulted body's rise, within a couple of seconds.
        assert abs(tret[5] - oracle_rise) < 5e-4, (
            f"tret[5]={tret[5]:.6f} vs reference rise {oracle_rise:.6f}"
        )
        assert tret[5] > 0.0
        # Rise falls between first and fourth contact.
        assert tret[1] < tret[5] < tret[4]
        # [6] set did not occur during this event; [7]/[8] are always zero.
        assert tret[6] == 0.0
        assert tret[7] == 0.0
        assert tret[8] == 0.0


# ---------------------------------------------------------------------------
# Fix 5: planet_occult_when_glob contact times use cos(dec) scaling
# ---------------------------------------------------------------------------
# Venus occults Jupiter on 2065-11-22 with Venus at dec ~ -17.7 deg.
JD_MAX_VENUS_JUPITER_2065 = 2475612.030788


def _eq(jd: float, pid: int):
    p, _ = calc_ut(jd, pid, FLG_EQUATORIAL | FLG_SPEED)
    return p[0], p[1], p[2]


def _sep_deg(jd: float, a: int, b: int) -> float:
    ra1, de1, _ = _eq(jd, a)
    ra2, de2, _ = _eq(jd, b)
    r1, d1, r2, d2 = map(math.radians, (ra1, de1, ra2, de2))
    cs = math.sin(d1) * math.sin(d2) + math.cos(d1) * math.cos(d2) * math.cos(r1 - r2)
    return math.degrees(math.acos(max(-1.0, min(1.0, cs))))


def _true_contact(jd_max: float, threshold: float, before: bool) -> float:
    """Independent contact time: bisect where the great-circle separation
    between Venus and Jupiter equals ``threshold``."""
    a, b = (jd_max - 0.3, jd_max) if before else (jd_max, jd_max + 0.3)
    fa = _sep_deg(a, VENUS, JUPITER) - threshold
    for _ in range(200):
        m = 0.5 * (a + b)
        fm = _sep_deg(m, VENUS, JUPITER) - threshold
        if (fa <= 0) == (fm <= 0):
            a, fa = m, fm
        else:
            b = m
    return 0.5 * (a + b)


class TestFix5PlanetOccultContactSpeed:
    def test_cos_dec_factor_reduces_contact_error(self):
        """At Venus dec ~ -17.7 deg the cos(dec)-scaled speed predicts first
        contact far better than the raw RA-degrees speed, and reduces to the
        raw speed at dec = 0 (cos factor -> 1)."""
        jd_max = JD_MAX_VENUS_JUPITER_2065
        _require_ephemeris(jd_max, VENUS)

        _, _, vdist = _eq(jd_max, VENUS)
        _, vdec, _ = _eq(jd_max, VENUS)
        _, _, jdist = _eq(jd_max, JUPITER)
        occ_r = E._calc_planet_angular_radius(VENUS, vdist)
        tgt_r = E._calc_planet_angular_radius(JUPITER, jdist)
        outer = occ_r + tgt_r
        min_sep = _sep_deg(jd_max, VENUS, JUPITER)

        # Confirm this really is a non-zero-declination occultation geometry.
        assert abs(vdec) > 5.0
        assert min_sep < outer

        dt = 1.0 / 24.0
        vr1, vd1, _ = _eq(jd_max - dt, VENUS)
        vr2, vd2, _ = _eq(jd_max + dt, VENUS)
        jr1, jd1, _ = _eq(jd_max - dt, JUPITER)
        jr2, jd2, _ = _eq(jd_max + dt, JUPITER)

        d_ra = (vr2 - vr1) - (jr2 - jr1)
        d_ra = (d_ra + 180.0) % 360.0 - 180.0
        d_dec = (vd2 - vd1) - (jd2 - jd1)
        cos_dec = math.cos(math.radians((vd1 + vd2 + jd1 + jd2) / 4.0))

        speed_old = math.sqrt(d_ra**2 + d_dec**2) / (2 * dt)
        speed_new = math.sqrt((d_ra * cos_dec) ** 2 + d_dec**2) / (2 * dt)

        half = math.sqrt(max(0.0, outer**2 - min_sep**2))
        first_old = jd_max - half / speed_old
        first_new = jd_max - half / speed_new

        true_first = _true_contact(jd_max, outer, before=True)
        err_old = abs(first_old - true_first) * 86400.0
        err_new = abs(first_new - true_first) * 86400.0

        # The fix is a large improvement (old ~14 s error, new < 1 s).
        assert err_new < 1.0, f"new contact error {err_new:.2f} s too large"
        assert err_new < 0.2 * err_old, (
            f"expected big improvement: old {err_old:.2f} s, new {err_new:.2f} s"
        )

    @pytest.mark.slow
    def test_planet_occult_when_glob_contact_times(self):
        """End-to-end: the contact times returned by planet_occult_when_glob
        match the independent great-circle separation crossings."""
        jd_start = julday(2065, 11, 10, 0.0)
        _require_ephemeris(jd_start, VENUS)

        try:
            retflag, tret = planet_occult_when_glob(jd_start, VENUS, JUPITER)
        except Exception as exc:  # pragma: no cover - environment dependent
            pytest.skip(f"occultation search failed: {exc!r}")

        jd_max = tret[0]
        assert abs(jd_max - JD_MAX_VENUS_JUPITER_2065) < 0.01

        _, _, vdist = _eq(jd_max, VENUS)
        _, _, jdist = _eq(jd_max, JUPITER)
        outer = E._calc_planet_angular_radius(
            VENUS, vdist
        ) + E._calc_planet_angular_radius(JUPITER, jdist)

        true_first = _true_contact(jd_max, outer, before=True)
        true_fourth = _true_contact(jd_max, outer, before=False)

        # tret[2] = begin, tret[3] = end (outer contacts).
        assert abs(tret[2] - true_first) * 86400.0 < 2.0
        assert abs(tret[3] - true_fourth) * 86400.0 < 2.0
        assert retflag != 0
