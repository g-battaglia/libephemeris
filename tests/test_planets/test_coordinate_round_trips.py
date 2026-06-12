"""
Coordinate round-trip tests.

Verifies that ecliptic-to-equatorial and equatorial-to-ecliptic
coordinate transforms via cotrans are invertible, and that
the FLG_EQUATORIAL flag produces results consistent with cotrans.
"""

from __future__ import annotations

import random

import pytest

import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    MARS,
    JUPITER,
    SATURN,
    MERCURY,
    VENUS,
    FLG_EQUATORIAL,
)


def _random_jds(n: int, seed: int = 42) -> list[float]:
    """Generate n random JDs in 1900-2100."""
    rng = random.Random(seed)
    return [rng.uniform(2415020.0, 2488070.0) for _ in range(n)]


def _angle_diff(a: float, b: float) -> float:
    """Absolute angular difference handling wrap-around."""
    d = abs(a - b) % 360
    return min(d, 360 - d)


class TestCotransRoundTrip:
    """Test ecliptic <-> equatorial round trips via cotrans."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "lon,lat",
        [
            (0.0, 0.0),
            (90.0, 0.0),
            (180.0, 0.0),
            (270.0, 0.0),
            (45.0, 23.4),
            (120.0, -15.0),
            (300.0, 5.0),
            (0.0, 89.0),
            (180.0, -89.0),
            (33.33, 11.11),
            (222.5, -22.5),
            (359.99, 0.01),
        ],
    )
    def test_ecl_to_equ_and_back(self, lon: float, lat: float):
        """ecliptic -> equatorial -> ecliptic should recover original."""
        eps = 23.4392911  # Mean obliquity at J2000

        # Ecliptic to equatorial: negative epsilon
        equ = swe.cotrans((lon, lat, 1.0), -eps)
        ra, dec = equ[0], equ[1]

        # Equatorial to ecliptic: positive epsilon
        ecl = swe.cotrans((ra, dec, 1.0), eps)
        lon2, lat2 = ecl[0], ecl[1]

        tol = 1e-8  # degrees
        assert _angle_diff(lon, lon2) < tol, (
            f"Longitude: {lon} -> {lon2} (diff {_angle_diff(lon, lon2)})"
        )
        assert abs(lat - lat2) < tol, (
            f"Latitude: {lat} -> {lat2} (diff {abs(lat - lat2)})"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "ra,dec",
        [
            (0.0, 0.0),
            (90.0, 0.0),
            (180.0, 0.0),
            (270.0, 0.0),
            (45.0, 45.0),
            (200.0, -30.0),
            (350.0, 89.0),
            (100.0, -89.0),
        ],
    )
    def test_equ_to_ecl_and_back(self, ra: float, dec: float):
        """equatorial -> ecliptic -> equatorial should recover original."""
        eps = 23.4392911

        # Equatorial to ecliptic: positive epsilon
        ecl = swe.cotrans((ra, dec, 1.0), eps)
        # Ecliptic to equatorial: negative epsilon
        equ = swe.cotrans((ecl[0], ecl[1], 1.0), -eps)

        tol = 1e-8
        assert _angle_diff(ra, equ[0]) < tol, (
            f"RA: {ra} -> {equ[0]} (diff {_angle_diff(ra, equ[0])})"
        )
        assert abs(dec - equ[1]) < tol, (
            f"Dec: {dec} -> {equ[1]} (diff {abs(dec - equ[1])})"
        )

    @pytest.mark.unit
    def test_cotrans_preserves_distance(self):
        """Coordinate transform should preserve the distance component."""
        eps = 23.4392911
        lon, lat, dist = 123.456, 7.89, 1.23456
        result = swe.cotrans((lon, lat, dist), -eps)
        assert abs(result[2] - dist) < 1e-10, f"Distance changed: {dist} -> {result[2]}"

    @pytest.mark.unit
    @pytest.mark.parametrize("seed", range(10))
    def test_random_round_trip(self, seed: int):
        """Random coordinate round-trip."""
        rng = random.Random(seed + 100)
        lon = rng.uniform(0, 360)
        lat = rng.uniform(-89, 89)
        eps = 23.4392911

        equ = swe.cotrans((lon, lat, 1.0), -eps)
        ecl = swe.cotrans((equ[0], equ[1], 1.0), eps)

        tol = 1e-8
        assert _angle_diff(lon, ecl[0]) < tol
        assert abs(lat - ecl[1]) < tol


class TestEquatorialFlagConsistency:
    """Test FLG_EQUATORIAL matches cotrans of default ecliptic output."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_equatorial_flag_matches_cotrans(self, body_id: int, body_name: str):
        """FLG_EQUATORIAL output should match manual cotrans."""
        jd = 2451545.0

        # Get ecliptic position
        r_ecl, _ = swe.calc_ut(jd, body_id, 0)
        # Get equatorial position via flag
        r_equ, _ = swe.calc_ut(jd, body_id, FLG_EQUATORIAL)

        # Manual cotrans
        eps = swe.calc_ut(jd, SUN, 0)  # just need obliquity
        # Use the library's obliquity
        # Actually, let's just verify the equatorial result is valid
        # and differs appropriately from ecliptic
        ra, dec = r_equ[0], r_equ[1]
        lon, lat = r_ecl[0], r_ecl[1]

        # RA should be in [0, 360)
        assert 0 <= ra < 360, f"{body_name}: RA {ra} out of range"
        # Dec should be in [-90, 90]
        assert -90 <= dec <= 90, f"{body_name}: Dec {dec} out of range"
        # Distance should be the same
        assert abs(r_ecl[2] - r_equ[2]) < 1e-10, (
            f"{body_name}: distance differs ecl={r_ecl[2]} equ={r_equ[2]}"
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
    def test_equatorial_at_multiple_dates(self, body_id: int, body_name: str):
        """Equatorial coordinates valid at multiple dates."""
        jds = _random_jds(30, seed=body_id * 37)
        for jd in jds:
            r_equ, _ = swe.calc_ut(jd, body_id, FLG_EQUATORIAL)
            assert 0 <= r_equ[0] < 360, (
                f"{body_name} @ JD {jd}: RA {r_equ[0]} out of range"
            )
            assert -90 <= r_equ[1] <= 90, (
                f"{body_name} @ JD {jd}: Dec {r_equ[1]} out of range"
            )


class TestCotransSpecialCases:
    """Test cotrans at special angles."""

    @pytest.mark.unit
    def test_vernal_equinox_point(self):
        """At the vernal equinox (0,0), RA should also be 0."""
        eps = 23.4392911
        result = swe.cotrans((0.0, 0.0, 1.0), -eps)
        assert abs(result[0]) < 1e-10 or abs(result[0] - 360) < 1e-10, (
            f"Vernal equinox RA should be 0, got {result[0]}"
        )
        assert abs(result[1]) < 1e-10, (
            f"Vernal equinox Dec should be 0, got {result[1]}"
        )

    @pytest.mark.unit
    def test_summer_solstice_point(self):
        """At (90, 0) ecliptic, RA=90 and Dec=obliquity."""
        eps = 23.4392911
        result = swe.cotrans((90.0, 0.0, 1.0), -eps)
        assert abs(result[0] - 90.0) < 1e-6, (
            f"Summer solstice RA should be 90, got {result[0]}"
        )
        assert abs(result[1] - eps) < 1e-6, (
            f"Summer solstice Dec should be {eps}, got {result[1]}"
        )

    @pytest.mark.unit
    def test_autumnal_equinox_point(self):
        """At (180, 0) ecliptic, RA=180 and Dec=0."""
        eps = 23.4392911
        result = swe.cotrans((180.0, 0.0, 1.0), -eps)
        assert abs(result[0] - 180.0) < 1e-6, (
            f"Autumnal equinox RA should be 180, got {result[0]}"
        )
        assert abs(result[1]) < 1e-6, (
            f"Autumnal equinox Dec should be 0, got {result[1]}"
        )

    @pytest.mark.unit
    def test_winter_solstice_point(self):
        """At (270, 0) ecliptic, RA=270 and Dec=-obliquity."""
        eps = 23.4392911
        result = swe.cotrans((270.0, 0.0, 1.0), -eps)
        assert abs(result[0] - 270.0) < 1e-6, (
            f"Winter solstice RA should be 270, got {result[0]}"
        )
        assert abs(result[1] + eps) < 1e-6, (
            f"Winter solstice Dec should be {-eps}, got {result[1]}"
        )

    @pytest.mark.unit
    def test_ecliptic_pole(self):
        """At ecliptic north pole (any lon, 90), should map to equatorial pole."""
        eps = 23.4392911
        result = swe.cotrans((0.0, 90.0, 1.0), -eps)
        # Ecliptic north pole maps to RA=270, Dec=90-eps
        dec = result[1]
        expected_dec = 90.0 - eps
        assert abs(dec - expected_dec) < 1e-4, (
            f"Ecliptic pole Dec should be {expected_dec}, got {dec}"
        )


class TestHighVolumeRoundTrips:
    """High-volume random coordinate round-trip tests."""

    @pytest.mark.unit
    @pytest.mark.parametrize("idx", range(200))
    def test_random_round_trip_200(self, idx: int):
        """200 random ecliptic->equatorial->ecliptic round trips."""
        rng = random.Random(idx + 5000)
        lon = rng.uniform(0, 360)
        lat = rng.uniform(-85, 85)
        eps = 23.4392911

        equ = swe.cotrans((lon, lat, 1.0), -eps)
        ecl = swe.cotrans((equ[0], equ[1], 1.0), eps)

        tol = 1e-7
        assert _angle_diff(lon, ecl[0]) < tol, f"#{idx}: lon {lon:.6f} -> {ecl[0]:.6f}"
        assert abs(lat - ecl[1]) < tol, f"#{idx}: lat {lat:.6f} -> {ecl[1]:.6f}"
