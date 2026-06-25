"""
Flag combinations stress test.

Tests that libephemeris handles a large number of flag combinations
without crashing, producing NaN, or returning wrong-length tuples.
Covers single flags, flag pairs, and random multi-flag combinations
across multiple bodies.
"""

from __future__ import annotations

import math
import random

import pytest

import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    PLUTO,
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
    OSCU_APOG,
    CHIRON,
    INTP_APOG,
    FLG_SPEED,
    FLG_HELCTR,
    FLG_TRUEPOS,
    FLG_J2000,
    FLG_NONUT,
    FLG_NOGDEFL,
    FLG_NOABERR,
    FLG_EQUATORIAL,
    FLG_XYZ,
    FLG_RADIANS,
    FLG_SIDEREAL,
    SIDM_LAHIRI,
)

# All single flags that can be combined
COMBINABLE_FLAGS = [
    (FLG_SPEED, "SPEED"),
    (FLG_TRUEPOS, "TRUEPOS"),
    (FLG_J2000, "J2000"),
    (FLG_NONUT, "NONUT"),
    (FLG_NOGDEFL, "NOGDEFL"),
    (FLG_NOABERR, "NOABERR"),
    (FLG_EQUATORIAL, "EQUATORIAL"),
    (FLG_XYZ, "XYZ"),
    (FLG_RADIANS, "RADIANS"),
]

GEOCENTRIC_BODIES = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (PLUTO, "Pluto"),
    (MEAN_NODE, "MeanNode"),
    (TRUE_NODE, "TrueNode"),
    (MEAN_APOG, "MeanApog"),
    (OSCU_APOG, "OscuApog"),
    (CHIRON, "Chiron"),
    (INTP_APOG, "IntpApog"),
]


class TestSingleFlags:
    """Test each single flag individually for each body."""

    @pytest.mark.unit
    @pytest.mark.parametrize("flag,flag_name", COMBINABLE_FLAGS)
    @pytest.mark.parametrize("body_id,body_name", GEOCENTRIC_BODIES)
    def test_single_flag_no_crash(
        self, body_id: int, body_name: str, flag: int, flag_name: str
    ):
        """Single flag produces a valid 6-element result."""
        jd = 2451545.0
        result, retflag = swe.calc_ut(jd, body_id, flag)
        assert len(result) == 6, f"{body_name}+{flag_name}: got {len(result)} elements"
        for i, val in enumerate(result):
            assert math.isfinite(val), (
                f"{body_name}+{flag_name}: result[{i}] = {val} (not finite)"
            )


class TestFlagPairs:
    """Test pairs of compatible flags."""

    # Generate all unique pairs
    _pairs = []
    for i, (f1, n1) in enumerate(COMBINABLE_FLAGS):
        for f2, n2 in COMBINABLE_FLAGS[i + 1 :]:
            _pairs.append((f1 | f2, f"{n1}|{n2}"))

    @pytest.mark.unit
    @pytest.mark.parametrize("flags,desc", _pairs)
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
        ],
    )
    def test_flag_pair_no_crash(
        self, body_id: int, body_name: str, flags: int, desc: str
    ):
        """Flag pair produces valid result without crash."""
        jd = 2451545.0
        result, retflag = swe.calc_ut(jd, body_id, flags)
        assert len(result) == 6, f"{body_name}+{desc}: got {len(result)} elements"
        for i, val in enumerate(result):
            assert math.isfinite(val), f"{body_name}+{desc}: result[{i}] = {val}"


class TestHeliocentricFlags:
    """Test heliocentric flag with various combinations."""

    # Bodies valid for heliocentric
    HELIO_BODIES = [
        (MERCURY, "Mercury"),
        (VENUS, "Venus"),
        (MARS, "Mars"),
        (JUPITER, "Jupiter"),
        (SATURN, "Saturn"),
    ]

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "extra,desc",
        [
            (0, "helctr"),
            (FLG_SPEED, "helctr+speed"),
            (FLG_EQUATORIAL, "helctr+equatorial"),
            (FLG_J2000, "helctr+J2000"),
            (FLG_NOABERR, "helctr+noaberr"),
            (FLG_TRUEPOS, "helctr+truepos"),
            (FLG_SPEED | FLG_EQUATORIAL, "helctr+speed+equatorial"),
            (FLG_SPEED | FLG_J2000, "helctr+speed+J2000"),
            (FLG_XYZ, "helctr+XYZ"),
            (FLG_XYZ | FLG_SPEED, "helctr+XYZ+speed"),
        ],
    )
    @pytest.mark.parametrize("body_id,body_name", HELIO_BODIES)
    def test_heliocentric_with_extra_flags(
        self, body_id: int, body_name: str, extra: int, desc: str
    ):
        """Heliocentric with extra flags produces valid output."""
        jd = 2451545.0
        flags = FLG_HELCTR | extra
        result, retflag = swe.calc_ut(jd, body_id, flags)
        assert len(result) == 6
        for i, val in enumerate(result):
            assert math.isfinite(val), f"{body_name}+{desc}: result[{i}] = {val}"


class TestSiderealFlags:
    """Test sidereal flag combinations."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "extra,desc",
        [
            (0, "sidereal"),
            (FLG_SPEED, "sidereal+speed"),
            (FLG_EQUATORIAL, "sidereal+equatorial"),
            (FLG_NOABERR, "sidereal+noaberr"),
            (FLG_TRUEPOS, "sidereal+truepos"),
            (FLG_J2000, "sidereal+J2000"),
            (FLG_SPEED | FLG_EQUATORIAL, "sidereal+speed+equatorial"),
        ],
    )
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
        ],
    )
    def test_sidereal_with_extra_flags(
        self, body_id: int, body_name: str, extra: int, desc: str
    ):
        """Sidereal with extra flags produces valid output."""
        swe.set_sid_mode(SIDM_LAHIRI)
        jd = 2451545.0
        flags = FLG_SIDEREAL | extra
        result, retflag = swe.calc_ut(jd, body_id, flags)
        assert len(result) == 6
        for i, val in enumerate(result):
            assert math.isfinite(val), f"{body_name}+{desc}: result[{i}] = {val}"


class TestRandomFlagCombinations:
    """Stress test with random combinations of flags."""

    @staticmethod
    def _make_random_flag_combos(n: int, seed: int = 42) -> list[tuple[int, str]]:
        """Generate n random flag combinations."""
        rng = random.Random(seed)
        flag_list = [
            FLG_SPEED,
            FLG_TRUEPOS,
            FLG_J2000,
            FLG_NONUT,
            FLG_NOGDEFL,
            FLG_NOABERR,
            FLG_EQUATORIAL,
        ]
        combos = []
        for _ in range(n):
            # Pick 1-5 random flags
            k = rng.randint(1, min(5, len(flag_list)))
            chosen = rng.sample(flag_list, k)
            combined = 0
            names = []
            for f in chosen:
                combined |= f
                names.append(str(f))
            combos.append((combined, "|".join(names)))
        return combos

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "flags,desc",
        _make_random_flag_combos.__func__(200, seed=42),
    )
    def test_random_flag_combo_sun(self, flags: int, desc: str):
        """Random flag combo on Sun doesn't crash."""
        jd = 2451545.0
        result, retflag = swe.calc_ut(jd, SUN, flags)
        assert len(result) == 6
        for i, val in enumerate(result):
            assert math.isfinite(val), f"Sun flags={desc}: result[{i}] = {val}"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "flags,desc",
        _make_random_flag_combos.__func__(200, seed=100),
    )
    def test_random_flag_combo_moon(self, flags: int, desc: str):
        """Random flag combo on Moon doesn't crash."""
        jd = 2451545.0
        result, retflag = swe.calc_ut(jd, MOON, flags)
        assert len(result) == 6
        for i, val in enumerate(result):
            assert math.isfinite(val), f"Moon flags={desc}: result[{i}] = {val}"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "flags,desc",
        _make_random_flag_combos.__func__(200, seed=200),
    )
    def test_random_flag_combo_mars(self, flags: int, desc: str):
        """Random flag combo on Mars doesn't crash."""
        jd = 2451545.0
        result, retflag = swe.calc_ut(jd, MARS, flags)
        assert len(result) == 6
        for i, val in enumerate(result):
            assert math.isfinite(val), f"Mars flags={desc}: result[{i}] = {val}"


class TestXYZRadians:
    """Tests for XYZ and RADIANS output modes."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
        ],
    )
    def test_xyz_returns_cartesian(self, body_id: int, body_name: str):
        """XYZ flag returns Cartesian coordinates."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, FLG_XYZ)
        x, y, z = result[0], result[1], result[2]
        # Cartesian distance should match spherical distance
        r_xyz = math.sqrt(x * x + y * y + z * z)
        result_sph, _ = swe.calc_ut(jd, body_id, 0)
        r_sph = result_sph[2]
        assert abs(r_xyz - r_sph) < 1e-6, (
            f"{body_name}: XYZ distance {r_xyz} vs spherical {r_sph}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
        ],
    )
    def test_radians_consistent_with_degrees(self, body_id: int, body_name: str):
        """RADIANS output should be degrees * pi/180."""
        jd = 2451545.0
        r_deg, _ = swe.calc_ut(jd, body_id, 0)
        r_rad, _ = swe.calc_ut(jd, body_id, FLG_RADIANS)
        # Longitude — allow 1e-8 tolerance due to internal rounding
        expected_rad = math.radians(r_deg[0])
        assert abs(r_rad[0] - expected_rad) < 1e-8, (
            f"{body_name}: radians lon {r_rad[0]} vs expected {expected_rad}"
        )
        # Latitude
        expected_lat_rad = math.radians(r_deg[1])
        assert abs(r_rad[1] - expected_lat_rad) < 1e-8, (
            f"{body_name}: radians lat {r_rad[1]} vs expected {expected_lat_rad}"
        )
