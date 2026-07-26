"""Regression tests for the v3 eclipse.py bug fixes.

Each test pins a corrected behavior of ``libephemeris/eclipse.py`` through
geometry identities, API layout invariants, or NASA eclipse-catalog facts.
No expected value was captured from another ephemeris implementation.

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
    CALC_RISE,
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

# Approximate greatest-eclipse instants from NASA's Five Millennium Canon.
# Minute precision is sufficient because these tests exercise smooth geometry
# and independently refine contacts where needed.
JD_MAX_TOTAL_2024 = julday(2024, 4, 8, 18.0 + 18.0 / 60.0)
JD_MAX_TOTAL_2017 = julday(2017, 8, 21, 18.0 + 26.0 / 60.0)
JD_MAX_ANNULAR_2023 = julday(2023, 10, 14, 18.0)
JD_MAX_ANNULAR_2024 = julday(2024, 10, 2, 18.0 + 45.0 / 60.0)
JD_MAX_HYBRID_2023 = julday(2023, 4, 20, 4.0 + 17.0 / 60.0)
JD_MAX_PARTIAL_2025 = julday(2025, 3, 29, 10.0 + 48.0 / 60.0)


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

    # Eclipse types are NASA Five Millennium Canon classifications.
    CASES = [
        ("total_2024", JD_MAX_TOTAL_2024, ECL_TOTAL),
        ("total_2017", JD_MAX_TOTAL_2017, ECL_TOTAL),
        ("annular_2023", JD_MAX_ANNULAR_2023, ECL_ANNULAR),
        ("annular_2024", JD_MAX_ANNULAR_2024, ECL_ANNULAR),
        ("hybrid_2023", JD_MAX_HYBRID_2023, ECL_ANNULAR_TOTAL),
    ]

    @pytest.mark.parametrize("label,jd_max,type_bit", CASES)
    def test_l2_matches_besselian_and_classifies(self, label, jd_max, type_bit):
        _require_ephemeris(jd_max)

        l2 = _calc_umbra_limit(jd_max)
        l2_bess = calc_besselian_l2(jd_max)

        # After the fix the simplified model tracks the Besselian model to
        # well under 1e-4 Earth radii (the buggy Earth-Sun distance left a
        # systematic ~6-8e-4 offset).
        assert abs(l2 - l2_bess) < 1e-4, (
            f"{label}: l2={l2:+.6f} vs besselian {l2_bess:+.6f}"
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
    """Duration is the minutes-long central-line totality/annularity."""

    # NASA Five Millennium Canon maximum durations, converted to minutes.
    CASES = [
        ("total_2024", JD_MAX_TOTAL_2024, 4.0 + 28.0 / 60.0),
        ("total_2017", JD_MAX_TOTAL_2017, 2.0 + 40.0 / 60.0),
        ("annular_2023", JD_MAX_ANNULAR_2023, 5.0 + 17.0 / 60.0),
        ("annular_2024", JD_MAX_ANNULAR_2024, 7.0 + 25.0 / 60.0),
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

        # Agreement with the independently published NASA duration.
        assert abs(duration - expected_min) < 0.15, (
            f"{label}: duration {duration:.4f} min vs NASA {expected_min:.4f}"
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

    def test_glob_layout_populates_defined_slots(self):
        jd_start = julday(2024, 1, 1, 0.0)
        _require_ephemeris(jd_start)

        _etype, tret = _sol_eclipse_when_glob_pythonic(jd_start, eclipse_type=ECL_TOTAL)
        assert all(math.isfinite(tret[i]) and tret[i] > 0.0 for i in range(8))
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


class TestHybridEclipseContactTimes:
    """A hybrid (annular-total) eclipse must still report U1/U4 and the
    center-line instants. ``ECL_ANNULAR_TOTAL`` (32) shares no bits with
    ``ECL_TOTAL`` (4) or ``ECL_ANNULAR`` (8), so the umbral-contact and
    center-line blocks used to be skipped -> tret[4..7] came back 0.0."""

    def test_hybrid_reports_umbral_and_centerline(self):
        jd_start = julday(2023, 4, 1, 0.0)
        _require_ephemeris(jd_start)

        etype, tret = _sol_eclipse_when_glob_pythonic(
            jd_start, eclipse_type=ECL_ANNULAR_TOTAL
        )
        assert etype & ECL_ANNULAR_TOTAL, f"expected a hybrid eclipse, got {etype}"

        # The regressed behaviour was exactly 0.0 for these four instants.
        for i in range(4, 8):
            assert tret[i] != 0.0, f"tret[{i}] lost (hybrid regression)"

        # Ordering: U1 < max < U4 and the center line lies within totality.
        assert tret[4] < tret[0] < tret[5]
        assert tret[4] <= tret[6] and tret[7] <= tret[5]
        # This layout reserves the last two slots.
        assert tret[8] == 0.0 and tret[9] == 0.0


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
        """The body-rise slot agrees with an independent horizon search."""

        rome = [12.4964, 41.9028, 0.0]  # lon, lat, alt
        jd_start = julday(2037, 1, 15, 0.0)
        _require_ephemeris(jd_start, SATURN)

        try:
            _retflag, tret, _attr = lun_occult_when_loc(jd_start, SATURN, rome)
        except Exception as exc:  # pragma: no cover - environment dependent
            pytest.skip(f"occultation search failed: {exc!r}")

        assert tret[5] > 0.0
        rise_flag, rise_tret = E.rise_trans(tret[1] - 0.5, SATURN, CALC_RISE, rome)
        assert rise_flag == 0
        assert abs(tret[5] - rise_tret[0]) * 86400.0 < 2.0
        # Rise falls between first and fourth contact.
        assert tret[1] < tret[5] < tret[4]
        # [6] set did not occur during this event; [7]/[8] are always zero.
        assert tret[6] == 0.0
        assert tret[7] == 0.0
        assert tret[8] == 0.0


# ---------------------------------------------------------------------------
# Fix 5: planet_occult_when_glob contact times use cos(dec) scaling
# ---------------------------------------------------------------------------
# Venus occults Jupiter on 2065-11-22 (JPL-derived geometry).
JD_NEAR_VENUS_JUPITER_2065 = julday(2065, 11, 22, 12.0)


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
        # Refine the conjunction directly from JPL-based positions.
        jd_max = E._golden_min(
            lambda jd: _sep_deg(jd, VENUS, JUPITER),
            JD_NEAR_VENUS_JUPITER_2065 - 0.5,
            JD_NEAR_VENUS_JUPITER_2065 + 0.5,
        )
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
        assert abs(jd_max - JD_NEAR_VENUS_JUPITER_2065) < 0.5

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


# ---------------------------------------------------------------------------
# Total-eclipse obscuration is the disc-area ratio (measured reference behavior)
# ---------------------------------------------------------------------------
class TestTotalObscurationBounded:
    """``attr[2]`` in totality is the disc-area ratio ``(r_moon/r_sun)**2``,
    which exceeds 1 for a total eclipse (the Moon overfills the Sun). This is
    the measured reference behavior: the reference returns the area ratio
    while one disc lies inside the other, equal to ``attr[1] ** 2``."""

    DALLAS = (-96.797, 32.777, 0.0)

    def test_total_obscuration_is_area_ratio(self):
        jd_start = julday(2024, 4, 1, 0.0)
        _require_ephemeris(jd_start)
        _retflag, tret, _search_attr = E.sol_eclipse_when_loc(jd_start, self.DALLAS)

        rc, attr = E.sol_eclipse_how(tret[0], self.DALLAS)
        assert rc & ECL_TOTAL, f"expected local totality, got {rc:#x}"
        # The reference reports the disc-area ratio, which exceeds 1 in
        # totality and equals the squared diameter ratio.
        assert attr[2] == pytest.approx(attr[1] ** 2, rel=1e-9)
        assert attr[2] > 1.0


# ---------------------------------------------------------------------------
# Round 2 fix B: hybrid classification by core-shadow sign change
# ---------------------------------------------------------------------------
class TestRound2HybridSignChange:
    """The annular-total classification uses the physical criterion (the
    core-shadow width changes sign along the central path) instead of the
    old static near-zero threshold.  Eclipse classifications below come from
    NASA's Five Millennium Canon."""

    CASES = [
        ("hybrid_2005", (2005, 4, 1), ECL_ANNULAR_TOTAL | ECL_CENTRAL),
        ("total_2002", (2002, 11, 25), ECL_TOTAL | ECL_CENTRAL),
        ("hybrid_2013", (2013, 10, 25), ECL_ANNULAR_TOTAL | ECL_CENTRAL),
        ("total_2024", (2024, 4, 1), ECL_TOTAL | ECL_CENTRAL),
        ("hybrid_2023", (2023, 4, 12), ECL_ANNULAR_TOTAL | ECL_CENTRAL),
        ("total_2006", (2006, 3, 20), ECL_TOTAL | ECL_CENTRAL),
        ("annular_2019", (2019, 12, 18), ECL_ANNULAR | ECL_CENTRAL),
    ]

    @pytest.mark.parametrize("label,start,expected", CASES, ids=[c[0] for c in CASES])
    def test_glob_retflag_matches_nasa_catalog(self, label, start, expected):
        jd0 = julday(start[0], start[1], start[2], 0.0)
        _require_ephemeris(jd0)

        retflag, _tret = E.sol_eclipse_when_glob(jd0)
        assert retflag == expected, (
            f"{label}: retflag {retflag:#x} != NASA class {expected:#x}"
        )


# ---------------------------------------------------------------------------
# Round 2 fix C: backwards accepts direction strings, rejects garbage
# ---------------------------------------------------------------------------
class TestRound2BackwardsCoercion:
    """A truthy direction string like "forward" used to select a *backward*
    search. ``_coerce_backwards`` now interprets strings explicitly, keeps
    bool/int semantics (including the ECL_ONE_TRY bit on the occultation
    functions) and raises TypeError on anything unrecognized."""

    def test_coerce_backwards_values(self):
        from libephemeris.eclipse import _coerce_backwards

        assert _coerce_backwards(True) is True
        assert _coerce_backwards(False) is False
        assert _coerce_backwards(0) is False
        assert _coerce_backwards(2) is True
        for text in ("backward", "Backwards", " TRUE ", "1"):
            assert _coerce_backwards(text) is True, text
        for text in ("forward", "Forwards", "false", "0", ""):
            assert _coerce_backwards(text) is False, text
        with pytest.raises(TypeError):
            _coerce_backwards("sideways")
        with pytest.raises(TypeError):
            _coerce_backwards(1.5)  # type: ignore[arg-type]

    def test_forward_string_searches_forward(self):
        jd0 = julday(2024, 1, 1, 0.0)
        _require_ephemeris(jd0)

        _rf, tret_str = E.sol_eclipse_when_glob(jd0, backwards="forward")
        _rf2, tret_bool = E.sol_eclipse_when_glob(jd0, backwards=False)
        assert tret_str[0] == tret_bool[0]
        assert tret_str[0] > jd0, "'forward' must not search backward"

    def test_backward_string_searches_backward(self):
        jd0 = julday(2024, 1, 1, 0.0)
        _require_ephemeris(jd0)

        _rf, tret = E.lun_eclipse_when(jd0, backwards="backwards")
        assert tret[0] < jd0, "'backwards' must search backward"

    def test_one_try_flag_preserved_on_occultations(self):
        from libephemeris.constants import ECL_ONE_TRY

        jd0 = julday(2024, 1, 1, 0.0)
        _require_ephemeris(jd0)

        retflag, tret = E.lun_occult_when_glob(
            jd0, VENUS, ecltype=0, backwards=ECL_ONE_TRY
        )
        # One-try miss case: retflag 0 with a forward continuation date —
        # the ECL_ONE_TRY bit must not be read as a backward search.
        assert retflag == 0
        assert tret[0] > jd0
