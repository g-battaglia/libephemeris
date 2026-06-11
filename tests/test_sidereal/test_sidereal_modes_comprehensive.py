"""
Comprehensive tests for sidereal ayanamsha modes.

Tests all 43 ayanamsha modes for consistency, value ranges, and
correct application to planetary positions and house cusps.
"""

from __future__ import annotations

import math
import random

import pytest

import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    MARS,
    JUPITER,
    SATURN,
    FLG_SIDEREAL,
    FLG_SPEED,
    SIDM_FAGAN_BRADLEY,
    SIDM_LAHIRI,
    SIDM_DELUCE,
    SIDM_RAMAN,
    SIDM_USHASHASHI,
    SIDM_KRISHNAMURTI,
    SIDM_DJWHAL_KHUL,
    SIDM_YUKTESHWAR,
    SIDM_JN_BHASIN,
    SIDM_BABYL_KUGLER1,
    SIDM_BABYL_KUGLER2,
    SIDM_BABYL_KUGLER3,
    SIDM_BABYL_HUBER,
    SIDM_BABYL_ETPSC,
    SIDM_ALDEBARAN_15TAU,
    SIDM_HIPPARCHOS,
    SIDM_SASSANIAN,
    SIDM_GALCENT_0SAG,
    SIDM_J2000,
    SIDM_J1900,
    SIDM_B1950,
    SIDM_SURYASIDDHANTA,
    SIDM_SURYASIDDHANTA_MSUN,
    SIDM_ARYABHATA,
    SIDM_ARYABHATA_MSUN,
    SIDM_SS_REVATI,
    SIDM_SS_CITRA,
    SIDM_TRUE_CITRA,
    SIDM_TRUE_REVATI,
    SIDM_TRUE_PUSHYA,
    SIDM_GALCENT_RGILBRAND,
    SIDM_GALEQU_IAU1958,
    SIDM_GALEQU_TRUE,
    SIDM_GALEQU_MULA,
    SIDM_GALALIGN_MARDYKS,
    SIDM_TRUE_MULA,
    SIDM_GALCENT_MULA_WILHELM,
    SIDM_ARYABHATA_522,
    SIDM_BABYL_BRITTON,
    SIDM_TRUE_SHEORAN,
    SIDM_GALCENT_COCHRANE,
    SIDM_GALEQU_FIORENZA,
    SIDM_VALENS_MOON,
)


ALL_MODES = [
    (SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
    (SIDM_LAHIRI, "Lahiri"),
    (SIDM_DELUCE, "De Luce"),
    (SIDM_RAMAN, "Raman"),
    (SIDM_USHASHASHI, "Ushashashi"),
    (SIDM_KRISHNAMURTI, "Krishnamurti"),
    (SIDM_DJWHAL_KHUL, "Djwhal Khul"),
    (SIDM_YUKTESHWAR, "Yukteshwar"),
    (SIDM_JN_BHASIN, "JN Bhasin"),
    (SIDM_BABYL_KUGLER1, "Babyl Kugler1"),
    (SIDM_BABYL_KUGLER2, "Babyl Kugler2"),
    (SIDM_BABYL_KUGLER3, "Babyl Kugler3"),
    (SIDM_BABYL_HUBER, "Babyl Huber"),
    (SIDM_BABYL_ETPSC, "Babyl ETPSC"),
    (SIDM_ALDEBARAN_15TAU, "Aldebaran 15Tau"),
    (SIDM_HIPPARCHOS, "Hipparchos"),
    (SIDM_SASSANIAN, "Sassanian"),
    (SIDM_GALCENT_0SAG, "Galcent 0Sag"),
    (SIDM_J2000, "J2000"),
    (SIDM_J1900, "J1900"),
    (SIDM_B1950, "B1950"),
    (SIDM_SURYASIDDHANTA, "Surya Siddhanta"),
    (SIDM_SURYASIDDHANTA_MSUN, "Surya MSun"),
    (SIDM_ARYABHATA, "Aryabhata"),
    (SIDM_ARYABHATA_MSUN, "Aryabhata MSun"),
    (SIDM_SS_REVATI, "SS Revati"),
    (SIDM_SS_CITRA, "SS Citra"),
    (SIDM_TRUE_CITRA, "True Citra"),
    (SIDM_TRUE_REVATI, "True Revati"),
    (SIDM_TRUE_PUSHYA, "True Pushya"),
    (SIDM_GALCENT_RGILBRAND, "Galcent Rgilbrand"),
    (SIDM_GALEQU_IAU1958, "Galequ IAU1958"),
    (SIDM_GALEQU_TRUE, "Galequ True"),
    (SIDM_GALEQU_MULA, "Galequ Mula"),
    (SIDM_GALALIGN_MARDYKS, "Galalign Mardyks"),
    (SIDM_TRUE_MULA, "True Mula"),
    (SIDM_GALCENT_MULA_WILHELM, "Galcent Mula Wilhelm"),
    (SIDM_ARYABHATA_522, "Aryabhata 522"),
    (SIDM_BABYL_BRITTON, "Babyl Britton"),
    (SIDM_TRUE_SHEORAN, "True Sheoran"),
    (SIDM_GALCENT_COCHRANE, "Galcent Cochrane"),
    (SIDM_GALEQU_FIORENZA, "Galequ Fiorenza"),
    (SIDM_VALENS_MOON, "Valens Moon"),
]

# Modes known to have unusual ayanamsha values (e.g., mode 40 ~357°)
UNUSUAL_MODES = {40}


class TestAyanamshaValues:
    """Tests for ayanamsha value computation across all modes."""

    @pytest.mark.unit
    @pytest.mark.parametrize("mode,name", ALL_MODES)
    def test_ayanamsha_returns_finite(self, mode: int, name: str):
        """Each ayanamsha mode returns a finite value."""
        swe.set_sid_mode(mode)
        jd = 2451545.0
        ayan = swe.get_ayanamsa_ut(jd)
        assert math.isfinite(ayan), f"{name} (mode {mode}): ayanamsha not finite"

    @pytest.mark.unit
    @pytest.mark.parametrize("mode,name", ALL_MODES)
    def test_ayanamsha_in_range(self, mode: int, name: str):
        """Ayanamsha at J2000 should be in a reasonable range."""
        swe.set_sid_mode(mode)
        jd = 2451545.0
        ayan = swe.get_ayanamsa_ut(jd)
        # Most ayanamshas are 20-30° at J2000; some are very different
        assert 0 <= ayan < 360, (
            f"{name} (mode {mode}): ayanamsha {ayan} out of [0, 360)"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "mode,name",
        [
            (SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
            (SIDM_LAHIRI, "Lahiri"),
            (SIDM_RAMAN, "Raman"),
            (SIDM_KRISHNAMURTI, "Krishnamurti"),
        ],
    )
    def test_major_ayanamsha_approximate_values(self, mode: int, name: str):
        """Major ayanamshas should be approximately 20-25° at J2000."""
        swe.set_sid_mode(mode)
        jd = 2451545.0
        ayan = swe.get_ayanamsa_ut(jd)
        assert 20 < ayan < 30, f"{name}: ayanamsha {ayan}° not in expected 20-30° range"

    @pytest.mark.unit
    @pytest.mark.parametrize("mode,name", ALL_MODES)
    def test_ayanamsha_increases_over_time(self, mode: int, name: str):
        """Ayanamsha should generally increase with time (precession)."""
        swe.set_sid_mode(mode)
        jd1 = 2451545.0  # J2000
        jd2 = jd1 + 365.25 * 100  # +100 years
        ayan1 = swe.get_ayanamsa_ut(jd1)
        ayan2 = swe.get_ayanamsa_ut(jd2)
        diff = (ayan2 - ayan1) % 360
        if diff > 180:
            diff -= 360
        # Precession is ~1.4°/century. Allow some modes to be zero-like
        assert diff > -5.0, f"{name}: ayanamsha decreased by {diff}° over 100 years"


class TestSiderealPositions:
    """Tests for sidereal planetary positions."""

    @pytest.mark.unit
    @pytest.mark.parametrize("mode,name", ALL_MODES)
    def test_sidereal_sun_valid(self, mode: int, name: str):
        """Sun sidereal position is valid for each mode."""
        swe.set_sid_mode(mode)
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, SUN, FLG_SIDEREAL)
        assert 0 <= result[0] < 360, (
            f"{name}: Sun sidereal lon {result[0]} out of range"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "mode,name",
        [
            (SIDM_LAHIRI, "Lahiri"),
            (SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
            (SIDM_RAMAN, "Raman"),
            (SIDM_TRUE_CITRA, "True Citra"),
            (SIDM_GALCENT_0SAG, "Galcent 0Sag"),
            (SIDM_J2000, "J2000"),
            (SIDM_BABYL_KUGLER1, "Babyl Kugler1"),
            (SIDM_HIPPARCHOS, "Hipparchos"),
            (SIDM_YUKTESHWAR, "Yukteshwar"),
            (SIDM_SURYASIDDHANTA, "Surya Siddhanta"),
        ],
    )
    def test_sidereal_offset_matches_ayanamsha(self, mode: int, name: str):
        """Tropical - sidereal should approximately equal ayanamsha."""
        swe.set_sid_mode(mode)
        jd = 2451545.0
        r_trop, _ = swe.calc_ut(jd, SUN, 0)
        r_sid, _ = swe.calc_ut(jd, SUN, FLG_SIDEREAL)
        ayan = swe.get_ayanamsa_ut(jd)

        diff = (r_trop[0] - r_sid[0]) % 360
        if diff > 350:
            diff -= 360
        ayan_norm = ayan % 360
        if ayan_norm > 350:
            ayan_norm -= 360

        assert abs(diff - ayan_norm) < 0.01, (
            f"{name}: tropical-sidereal diff {diff:.4f} vs ayanamsha {ayan_norm:.4f}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_sidereal_lahiri_multiple_bodies(self, body_id: int, body_name: str):
        """Lahiri sidereal positions valid for multiple bodies."""
        swe.set_sid_mode(SIDM_LAHIRI)
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, FLG_SIDEREAL | FLG_SPEED)
        assert 0 <= result[0] < 360
        assert math.isfinite(result[3])

    @pytest.mark.unit
    @pytest.mark.parametrize("mode,name", ALL_MODES)
    def test_sidereal_multiple_dates(self, mode: int, name: str):
        """Sidereal position valid across multiple dates for each mode."""
        swe.set_sid_mode(mode)
        jds = [2415020.0, 2440587.5, 2451545.0, 2460000.0]
        for jd in jds:
            result, _ = swe.calc_ut(jd, SUN, FLG_SIDEREAL)
            assert 0 <= result[0] < 360, f"{name}: lon {result[0]} at JD {jd}"
            assert math.isfinite(result[1])

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "mode,name",
        [
            (SIDM_LAHIRI, "Lahiri"),
            (SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
        ],
    )
    def test_sidereal_speed_close_to_tropical(self, mode: int, name: str):
        """Sidereal speed should be close to tropical speed (differ by ~precession)."""
        swe.set_sid_mode(mode)
        jd = 2451545.0
        r_trop, _ = swe.calc_ut(jd, SUN, FLG_SPEED)
        r_sid, _ = swe.calc_ut(jd, SUN, FLG_SIDEREAL | FLG_SPEED)
        # Speed difference should be tiny (~0.004°/day for precession)
        assert abs(r_trop[3] - r_sid[3]) < 0.01, (
            f"{name}: tropical speed {r_trop[3]} vs sidereal {r_sid[3]}"
        )


class TestSiderealHouses:
    """Tests for sidereal house cusps."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "mode,name",
        [
            (SIDM_LAHIRI, "Lahiri"),
            (SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
            (SIDM_RAMAN, "Raman"),
            (SIDM_KRISHNAMURTI, "Krishnamurti"),
            (SIDM_TRUE_CITRA, "True Citra"),
        ],
    )
    def test_sidereal_houses_valid(self, mode: int, name: str):
        """Sidereal houses return 12 valid cusps."""
        swe.set_sid_mode(mode)
        jd = 2451545.0
        cusps, ascmc = swe.houses_ex(jd, 41.9, 12.5, ord("P"), FLG_SIDEREAL)
        assert len(cusps) >= 12
        for i, cusp in enumerate(cusps[:12]):
            assert 0 <= cusp < 360, f"{name}: cusp {i + 1} = {cusp} out of range"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "mode,name",
        [
            (SIDM_LAHIRI, "Lahiri"),
            (SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
        ],
    )
    def test_sidereal_cusps_offset_from_tropical(self, mode: int, name: str):
        """Sidereal cusps should differ from tropical by ~ayanamsha."""
        swe.set_sid_mode(mode)
        jd = 2451545.0
        lat, lon = 41.9, 12.5
        cusps_trop, _ = swe.houses(jd, lat, lon, ord("P"))
        cusps_sid, _ = swe.houses_ex(jd, lat, lon, ord("P"), FLG_SIDEREAL)
        ayan = swe.get_ayanamsa_ut(jd)

        for i in range(12):
            diff = (cusps_trop[i] - cusps_sid[i]) % 360
            if diff > 350:
                diff -= 360
            ayan_norm = ayan % 360
            if ayan_norm > 350:
                ayan_norm -= 360
            assert abs(diff - ayan_norm) < 0.1, (
                f"{name}: cusp {i + 1} diff {diff:.3f} vs ayanamsha {ayan_norm:.3f}"
            )


class TestAyanamshaExUt:
    """Tests for the extended ayanamsha function."""

    @pytest.mark.unit
    @pytest.mark.parametrize("mode,name", ALL_MODES)
    def test_get_ayanamsa_ex_ut_returns_correct_retflag(self, mode: int, name: str):
        """get_ayanamsa_ex_ut should return FLG_SWIEPH (2) as retflag."""
        swe.set_sid_mode(mode)
        jd = 2451545.0
        retflag, ayan = swe.get_ayanamsa_ex_ut(jd, 0)
        # BUG-004 fix: retflag should be 2 (FLG_SWIEPH), not the input flags
        assert retflag == 2, f"{name}: retflag {retflag} != 2 (FLG_SWIEPH)"
        assert math.isfinite(ayan)

    @pytest.mark.unit
    @pytest.mark.parametrize("mode,name", ALL_MODES)
    def test_ayanamsa_ex_ut_matches_ayanamsa_ut(self, mode: int, name: str):
        """The ex variant returns the TRUE ayanamsha (mean + nutation).

        Reference API semantics (verified against pyswisseph): the plain
        function returns the mean ayanamsha; the _ex variant adds nutation
        in longitude unless FLG_NONUT is set.
        """
        swe.set_sid_mode(mode)
        jd = 2451545.0
        ayan_simple = swe.get_ayanamsa_ut(jd)
        _, ayan_ex = swe.get_ayanamsa_ex_ut(jd, 0)
        # The difference is the nutation in longitude (|dpsi| <= ~17.5");
        # wrap-aware because ayanamshas near 0 deg can wrap to ~360.
        delta = abs(ayan_simple - ayan_ex) % 360.0
        delta = min(delta, 360.0 - delta)
        assert delta < 0.006, f"{name}: simple {ayan_simple} vs ex {ayan_ex}"
        # With FLG_NONUT the ex variant returns the mean value exactly
        _, ayan_ex_nonut = swe.get_ayanamsa_ex_ut(jd, swe.FLG_NONUT)
        assert abs(ayan_simple - ayan_ex_nonut) < 1e-10, (
            f"{name}: simple {ayan_simple} vs ex(NONUT) {ayan_ex_nonut}"
        )
