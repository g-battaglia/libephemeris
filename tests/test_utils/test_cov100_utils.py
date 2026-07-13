# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage completion tests for ``libephemeris.utils``.

These tests target the specific lines and branch arcs left uncovered by the
broader unit suite: the split-form of ``cotrans_sp``, its pole/degenerate
velocity branches, the ``@overload`` stub bodies, the refraction early-return
paths, the centisecond formatters, ``csroundsec`` boundary handling,
``split_deg``/nakshatra decomposition, ``d2l`` rounding, and ``calc_angles``.

Inputs are deterministic (fixed angles, J2000 Julian Day) so the assertions
remain stable and require no network access.
"""

from __future__ import annotations

import math
import typing

import pytest

from libephemeris.utils import (
    APP_TO_TRUE,
    EQU2HOR,
    SPLIT_DEG_KEEP_DEG,
    SPLIT_DEG_KEEP_SIGN,
    SPLIT_DEG_NAKSHATRA,
    SPLIT_DEG_ROUND_DEG,
    SPLIT_DEG_ROUND_MIN,
    SPLIT_DEG_ROUND_SEC,
    SPLIT_DEG_ZODIACAL,
    TRUE_TO_APP,
    azalt,
    calc_angles,
    cotrans_sp,
    cs2lonlatstr,
    cs2timestr,
    csnorm,
    csroundsec,
    d2l,
    degnorm,
    difcs2n,
    difcsn,
    difdegn,
    difrad2n,
    radnorm,
    refrac,
    split_deg,
)
from libephemeris.utils import _decompose_to_dms, _rounding_offset

TWO_PI = 2.0 * math.pi


class TestCotransSpOverloadStubs:
    """Cover the ``@overload`` stub bodies of ``cotrans_sp`` (lines 30, 36)."""

    @pytest.mark.unit
    def test_overload_stubs_callable(self):
        """Invoking the registered overload stubs executes their ``...`` body."""
        overloads = typing.get_overloads(cotrans_sp)
        assert len(overloads) == 2
        # First overload: 2-arg form (coord, eps). Body is `...` -> None.
        assert overloads[0]([1.0, 2.0, 3.0], 23.4) is None
        # Second overload: 3-arg form (coord, speed, eps). Body is `...` -> None.
        assert overloads[1]([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], 23.4) is None


class TestCotransSpSplitForm:
    """Cover the 3-arg split form and degenerate velocity branches."""

    @pytest.mark.unit
    def test_split_form_returns_pair(self):
        """3-arg form returns (coord_3, speed_3) and matches the flat form."""
        coord3 = (85.0, 5.0, 1.0)
        speed3 = (0.9, 0.1, 0.0)
        eps = 23.44

        coord_out, speed_out = cotrans_sp(coord3, speed3, eps)
        flat = cotrans_sp(
            (coord3[0], coord3[1], coord3[2], speed3[0], speed3[1], speed3[2]),
            eps,
        )

        assert len(coord_out) == 3
        assert len(speed_out) == 3
        for a, b in zip((*coord_out, *speed_out), flat):
            assert abs(a - b) < 1e-9

    @pytest.mark.unit
    def test_pole_branches_zero_latitude_speed(self):
        """Input latitude 90 hits cos_new_lat and cos_lat_sq pole branches.

        Lines 147 (latitude-speed undefined at pole) and 164 (longitude-speed
        derivative without the sec^2 term) are exercised here.
        """
        (lon, lat, dist), (lon_sp, lat_sp, dist_sp) = cotrans_sp(
            [10.0, 90.0, 1.0], [1.0, 1.0, 0.0], 0.0
        )
        assert abs(lat - 90.0) < 1e-9
        assert lat_sp == 0.0
        assert dist_sp == 0.0

    @pytest.mark.unit
    def test_denominator_degenerate_branch(self):
        """Choose inputs so x^2 + y^2 ~ 0, hitting the denom guard (line 170)."""
        eps = 23.0
        # cos(lon)=0 at lon=90; pick lat so the rotated y-component vanishes too.
        lat = math.degrees(math.atan(1.0 / math.tan(math.radians(-eps))))
        (_, _, _), (lon_sp, lat_sp, dist_sp) = cotrans_sp(
            [90.0, lat, 1.0], [1.0, 1.0, 0.0], eps
        )
        assert lon_sp == 0.0


class TestAzaltAndAzaltRevDegenerate:
    """Cover azalt azimuth wrap usage and azalt_rev pole branch (line 442)."""

    @pytest.mark.unit
    def test_azalt_basic_call(self):
        """A plain azalt call exercises the azimuth-normalisation path."""
        jd = 2451545.0
        geopos = (12.5, 41.9, 0.0)
        az, alt_true, alt_app = azalt(
            jd, EQU2HOR, geopos, 1013.25, 15.0, (90.0, 23.5, 1.0)
        )
        assert 0.0 <= az < 360.0
        assert -90.0 <= alt_true <= 90.0

    @pytest.mark.unit
    def test_azalt_rev_object_at_celestial_pole(self):
        """azalt_rev with object near the zenith of a pole observer.

        Observer at the geographic north pole looking straight up places the
        object at the celestial pole, where the hour angle is undefined and
        the code returns ha = 0.0 (line 442).
        """
        from libephemeris.utils import azalt_rev, HOR2EQU

        jd = 2451545.0
        geopos = (0.0, 90.0, 0.0)  # north pole observer
        ra, dec = azalt_rev(jd, HOR2EQU, geopos, 0.0, 90.0)
        assert abs(dec - 90.0) < 1e-6


class TestRefracEarlyReturns:
    """Cover plain refrac() compatibility clamps."""

    @pytest.mark.unit
    def test_true_to_app_below_horizon_returns_input(self):
        """Deep negative true altitude cannot be lifted above 0 -> returns input."""
        alt = -5.0
        out = refrac(alt, 1013.25, 15.0, TRUE_TO_APP)
        assert out == alt

    @pytest.mark.unit
    def test_app_to_true_below_horizon_threshold_returns_input(self):
        """A reverse candidate that is not positive leaves the input unchanged."""
        alt = -1.0
        out = refrac(alt, 1013.25, 15.0, APP_TO_TRUE)
        assert out == alt


class TestDegnormRadnormSnap:
    """Cover the [0, upper) snap-to-zero guards (lines 765 and 804)."""

    @pytest.mark.unit
    def test_degnorm_tiny_negative_snaps_to_zero(self):
        """(-1e-17) % 360.0 == 360.0 must snap to 0.0."""
        assert (-1e-17) % 360.0 == 360.0  # precondition for the branch
        assert degnorm(-1e-17) == 0.0

    @pytest.mark.unit
    def test_radnorm_tiny_negative_snaps_to_zero(self):
        """A tiny negative radian input whose modulo lands on 2*pi snaps to 0.0."""
        assert (-1e-16) % TWO_PI >= TWO_PI  # precondition for the branch
        assert radnorm(-1e-16) == 0.0


class TestDifferenceHelpers:
    """Cover difdegn, difrad2n, difcs2n, difcsn return paths."""

    @pytest.mark.unit
    def test_difdegn(self):
        """difdegn returns the positive [0, 360) difference (line 866)."""
        assert difdegn(10.0, 20.0) == 350.0
        assert difdegn(20.0, 10.0) == 10.0

    @pytest.mark.unit
    def test_difrad2n_both_branches(self):
        """difrad2n normalises to [-pi, pi] (lines 894-897)."""
        # diff >= pi triggers the subtraction branch (line 896).
        result = difrad2n(math.pi + 0.1, 0.0)
        assert abs(result - (math.pi + 0.1 - TWO_PI)) < 1e-9
        assert result < 0.0
        # small diff stays as-is.
        assert abs(difrad2n(0.2, 0.1) - 0.1) < 1e-9

    @pytest.mark.unit
    def test_difcs2n_both_branches(self):
        """difcs2n normalises to [-180, 180] in centiseconds (lines 941-944)."""
        # exactly 180 deg -> -180 deg per reference convention.
        assert difcs2n(64800000, 0) == -64800000
        # small positive difference unchanged.
        assert difcs2n(720000, 360000) == 360000

    @pytest.mark.unit
    def test_difcsn(self):
        """difcsn returns the positive centisecond difference (line 982)."""
        assert difcsn(360000, 720000) == 129240000
        assert difcsn(720000, 360000) == 360000


class TestCsnorm:
    """Cover csnorm return path (line 1022)."""

    @pytest.mark.unit
    def test_csnorm(self):
        assert csnorm(-360000) == 129240000
        assert csnorm(360000) == 360000
        assert csnorm(129600000) == 0


class TestCsroundsec:
    """Cover every branch of csroundsec (lines 1062-1088)."""

    @pytest.mark.unit
    def test_positive_simple_round(self):
        """Positive value, no boundary correction (lines 1063-1064)."""
        assert csroundsec(150) == 200
        assert csroundsec(149) == 100

    @pytest.mark.unit
    def test_positive_sign_boundary_rounddown(self):
        """Positive value just below a 30-deg boundary rounds down (line 1075)."""
        # cs in [10799950, 10799999] rounds up to 10800000, then -100.
        assert csroundsec(10799960) == 10799900

    @pytest.mark.unit
    def test_negative_simple_truncation(self):
        """Negative value uses truncation-toward-zero (lines 1066-1067)."""
        assert csroundsec(-150) == -100

    @pytest.mark.unit
    def test_negative_sign_boundary_rounddown(self):
        """Negative value at a 30-deg boundary rounds down further (line 1082)."""
        assert csroundsec(-10800120) == -10800100

    @pytest.mark.unit
    def test_negative_small_to_minus_100(self):
        """Negative value in (-150, -100] with zero result snaps to -100 (1086)."""
        assert csroundsec(-120) == -100

    @pytest.mark.unit
    def test_negative_small_to_zero(self):
        """Negative value in (-50, 0) rounds to 0 (line 1085 false path)."""
        assert csroundsec(-40) == 0


class TestCs2LonLatStr:
    """Cover bytes decode and the seconds carry-over branch."""

    @pytest.mark.unit
    def test_bytes_plus_and_minus_decoded(self):
        """Both ``plus`` and ``minus`` decode from bytes (lines 1187, 1189)."""
        assert cs2lonlatstr(360000, b"N", b"S") == "1N00"
        assert cs2lonlatstr(-360000, b"E", b"W") == "1W00"

    @pytest.mark.unit
    def test_seconds_carry_without_minute_overflow(self):
        """Seconds round to 60, carry to minutes, minutes stay < 60 (1213->1219)."""
        # remainder 5970 -> seconds=(6020)//100=60 -> carry, minutes 0+1=1 < 60.
        assert cs2lonlatstr(5970, b"N", b"S") == "0N01"

    @pytest.mark.unit
    def test_seconds_nonzero_output(self):
        """Non-zero seconds produce the full minute'second format."""
        assert cs2lonlatstr(12345678, b"E", b"W") == "34E17'37"


class TestCs2TimeStr:
    """Cover bytes separator decode and suppresszero branch."""

    @pytest.mark.unit
    def test_bytes_separator_decoded(self):
        """A bytes separator decodes to str (line 1262)."""
        assert cs2timestr(360000, b":") == "01:00:00"

    @pytest.mark.unit
    def test_suppresszero_drops_seconds(self):
        """suppresszero with zero seconds drops the seconds field (1299-1300)."""
        assert cs2timestr(360000, ":", True) == "01:00"

    @pytest.mark.unit
    def test_suppresszero_keeps_seconds_when_nonzero(self):
        """suppresszero with non-zero seconds keeps the seconds field."""
        assert cs2timestr(360150, ":", True) == "01:00:02"


class TestD2l:
    """Cover d2l rounding (lines 1449-1452)."""

    @pytest.mark.unit
    def test_positive_rounds_half_away(self):
        assert d2l(1.5) == 2
        assert d2l(1.4) == 1

    @pytest.mark.unit
    def test_negative_rounds_half_away(self):
        assert d2l(-1.5) == -2
        assert d2l(-1.4) == -1


class TestRoundingOffsetAndDecompose:
    """Cover _rounding_offset (1471-1477) and _decompose_to_dms (1490-1499)."""

    @pytest.mark.unit
    def test_rounding_offset_levels(self):
        """Each rounding flag selects its offset; no flag returns 0.0."""
        assert _rounding_offset(SPLIT_DEG_ROUND_DEG) == 0.5
        assert _rounding_offset(SPLIT_DEG_ROUND_MIN) == 0.5 / 60.0
        assert _rounding_offset(SPLIT_DEG_ROUND_SEC) == 0.5 / 3600.0
        assert _rounding_offset(0) == 0.0

    @pytest.mark.unit
    def test_decompose_without_rounding(self):
        """secfr keeps the true fractional part when has_rounding is False."""
        ideg, imin, isec, secfr = _decompose_to_dms(45.50, False)
        assert ideg == 45
        assert imin == 30
        assert isec == 0
        assert 0.0 <= secfr < 1.0

    @pytest.mark.unit
    def test_decompose_with_rounding(self):
        """secfr becomes the integer-seconds value when has_rounding is True."""
        ideg, imin, isec, secfr = _decompose_to_dms(10.0 + 34.56 / 3600.0, True)
        assert secfr == float(isec)


class TestSplitDegNakshatra:
    """Cover the nakshatra decomposition path (lines 1510-1545)."""

    @pytest.mark.unit
    def test_nakshatra_basic(self):
        """A plain nakshatra split returns the within-nakshatra position."""
        result = split_deg(45.5, SPLIT_DEG_NAKSHATRA)
        assert result[4] == 3  # nakshatra index

    @pytest.mark.unit
    def test_nakshatra_round_deg(self):
        """Rounding to degrees applies the offset and re-reduces."""
        result = split_deg(13.3, SPLIT_DEG_NAKSHATRA | SPLIT_DEG_ROUND_DEG)
        assert isinstance(result[0], int)

    @pytest.mark.unit
    def test_nakshatra_keep_deg_suppresses_offset(self):
        """KEEP_DEG suppresses an offset that would bump the integer degree."""
        # pos_in_nak 5.6, +0.5 would become 6.x -> suppressed, stays 5 deg.
        result = split_deg(
            5.6, SPLIT_DEG_NAKSHATRA | SPLIT_DEG_ROUND_DEG | SPLIT_DEG_KEEP_DEG
        )
        assert result[0] == 5

    @pytest.mark.unit
    def test_nakshatra_keep_sign_suppresses_offset(self):
        """KEEP_SIGN suppresses an offset that would cross into next nakshatra."""
        result = split_deg(
            13.0, SPLIT_DEG_NAKSHATRA | SPLIT_DEG_ROUND_DEG | SPLIT_DEG_KEEP_SIGN
        )
        assert result[4] == 0  # stays in nakshatra 0

    @pytest.mark.unit
    def test_nakshatra_keep_deg_offset_not_suppressed(self):
        """KEEP_DEG set but offset does not bump the integer (arc 1524->1535)."""
        # pos_in_nak ~5.0; +0.5 stays within the same integer degree.
        result = split_deg(
            5.0, SPLIT_DEG_NAKSHATRA | SPLIT_DEG_ROUND_DEG | SPLIT_DEG_KEEP_DEG
        )
        assert result[0] == 5

    @pytest.mark.unit
    def test_nakshatra_keep_sign_offset_not_suppressed(self):
        """KEEP_SIGN set but offset stays inside the nakshatra (arc 1528->1535)."""
        result = split_deg(
            5.0, SPLIT_DEG_NAKSHATRA | SPLIT_DEG_ROUND_DEG | SPLIT_DEG_KEEP_SIGN
        )
        assert result[4] == 0

    @pytest.mark.unit
    def test_nakshatra_index_wraps_at_27(self):
        """An offset pushing the index to 27 wraps it back to 0 (lines 1537-1538)."""
        result = split_deg(359.9, SPLIT_DEG_NAKSHATRA | SPLIT_DEG_ROUND_DEG)
        assert result[4] == 0

    @pytest.mark.unit
    def test_nakshatra_round_deg_uses_exact_segment_geometry(self):
        """Rounding follows the exact 360/27 segment geometry."""
        result = split_deg(0.09, SPLIT_DEG_NAKSHATRA | SPLIT_DEG_ROUND_DEG)
        assert result[4] == 0
        assert result[:2] == (0, 35)

    @pytest.mark.unit
    def test_nakshatra_decimal_boundary_normalizes(self):
        """A longitude beyond one turn is normalized before segmentation."""
        result = split_deg(373.3333333333333, SPLIT_DEG_NAKSHATRA)
        assert result[4] == 0

    @pytest.mark.unit
    def test_negative_nakshatra_flag_uses_ordinary_signed_split(self):
        """Negative input bypasses nakshatra segmentation by contract."""
        assert split_deg(-30.5, SPLIT_DEG_NAKSHATRA) == (30, 30, 0, 0.0, -1)

    @pytest.mark.unit
    def test_zodiac_indices_normalize_and_combined_negative_flag_precedence(self):
        """Zodiac indices normalize and negative combined flags stay defined."""
        assert split_deg(370.0, SPLIT_DEG_ZODIACAL)[4] == 0
        assert split_deg(390.0, SPLIT_DEG_ZODIACAL)[4] == 1
        assert split_deg(480.0, SPLIT_DEG_ZODIACAL)[4] == 4
        assert split_deg(-370.0, SPLIT_DEG_ZODIACAL)[4] == 0
        assert split_deg(-720.0, SPLIT_DEG_ZODIACAL)[4] == 0
        assert split_deg(-30.5, SPLIT_DEG_NAKSHATRA | SPLIT_DEG_ZODIACAL)[4] == 1

    @pytest.mark.unit
    def test_nakshatra_docstring_documents_normalized_indices(self):
        """The public docs describe the independent normalized semantics."""
        doc = split_deg.__doc__ or ""
        assert "27 exact equal segments" in doc
        assert "NAKSHATRA alone follows the ordinary signed split" in doc
        assert "non-negative NAKSHATRA split" in doc


class TestSplitDegMain:
    """Cover the main split_deg path (lines 1604-1645)."""

    @pytest.mark.unit
    def test_negative_input_sets_negative_sign(self):
        """Negative input -> sign_out = -1 (lines 1604-1606)."""
        result = split_deg(-30.5, 0)
        assert result[4] == -1
        assert result[0] == 30

    @pytest.mark.unit
    def test_positive_no_flags(self):
        """Positive input with no flags -> sign_out = +1 (lines 1610-1611)."""
        result = split_deg(45.5, 0)
        assert result == (45, 30, 0, 0.0, 1)

    @pytest.mark.unit
    def test_zodiacal(self):
        """ZODIACAL extracts a 0-11 sign (lines 1636-1640)."""
        result = split_deg(45.5, SPLIT_DEG_ZODIACAL)
        assert result == (15, 30, 0, 0.0, 1)

    @pytest.mark.unit
    def test_zodiacal_sign_twelve_wraps_to_zero(self):
        """A value rounding to 360 maps sign 12 back to 0 (lines 1638-1639)."""
        result = split_deg(359.9, SPLIT_DEG_ZODIACAL | SPLIT_DEG_ROUND_DEG)
        assert result[4] == 0

    @pytest.mark.unit
    def test_keep_deg_suppresses_offset(self):
        """KEEP_DEG suppresses an offset bumping the integer degree (1626-1627)."""
        result = split_deg(45.6, SPLIT_DEG_ROUND_DEG | SPLIT_DEG_KEEP_DEG)
        assert result[0] == 45

    @pytest.mark.unit
    def test_keep_sign_suppresses_offset(self):
        """KEEP_SIGN suppresses an offset crossing a 30-deg boundary (1630-1631)."""
        result = split_deg(
            29.9, SPLIT_DEG_ROUND_DEG | SPLIT_DEG_KEEP_SIGN | SPLIT_DEG_ZODIACAL
        )
        assert result[4] == 0  # still in sign 0

    @pytest.mark.unit
    def test_round_min_and_sec(self):
        """ROUND_MIN and ROUND_SEC each engage the offset path."""
        assert split_deg(45.999, SPLIT_DEG_ROUND_MIN)[1] == 0
        assert split_deg(45.99999, SPLIT_DEG_ROUND_SEC)[0] == 46


class TestCalcAngles:
    """Cover calc_angles (lines 1661-1686)."""

    @pytest.mark.unit
    def test_calc_angles_returns_dict_with_planets(self):
        """calc_angles caches and returns angles plus Sun/Moon/Mercury/Venus."""
        result = calc_angles(2451545.0, 41.9, 12.5)
        assert isinstance(result, dict)
        for key in ("Sun", "Moon", "Mercury", "Venus"):
            assert key in result
            assert isinstance(result[key], float)
