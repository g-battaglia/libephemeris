"""
Fixed star proper motion tests.

Verifies that fixed stars show correct proper motion across centuries,
and that the fixstar_ut function handles various star names and dates.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.constants import FLG_SPEED


# Well-known bright stars with significant proper motion
BRIGHT_STARS = [
    "Sirius",
    "Arcturus",
    "Vega",
    "Capella",
    "Rigel",
    "Procyon",
    "Betelgeuse",
    "Altair",
    "Aldebaran",
    "Spica",
    "Antares",
    "Pollux",
    "Fomalhaut",
    "Deneb",
    "Regulus",
]

# Stars with large proper motion (arcseconds/year)
HIGH_PROPER_MOTION_STARS = [
    ("Sirius", 1.3),  # ~1.3"/yr total
    ("Arcturus", 2.3),  # ~2.3"/yr total
    ("Procyon", 1.3),  # ~1.3"/yr total
    ("Altair", 0.66),  # ~0.66"/yr total
    ("Pollux", 0.63),  # ~0.63"/yr total
]


class TestFixstarBasic:
    """Basic fixed star functionality tests."""

    @pytest.mark.unit
    @pytest.mark.parametrize("star_name", BRIGHT_STARS)
    def test_fixstar_returns_valid_position(self, star_name: str):
        """fixstar_ut returns valid position for known bright stars."""
        jd = 2451545.0
        result = swe.fixstar_ut(star_name, jd, 0)
        # Returns (name, (lon, lat, dist, slon, slat, sdist), retflag)
        assert len(result) >= 2
        pos = result[1] if isinstance(result[1], tuple) else result[0]
        # If result format is (pos_tuple, retflag)
        if isinstance(result[0], tuple):
            pos = result[0]
        elif isinstance(result[0], str):
            pos = result[1]

        assert len(pos) == 6, f"{star_name}: expected 6 elements, got {len(pos)}"
        lon = pos[0]
        assert 0 <= lon < 360, f"{star_name}: longitude {lon} out of range"

    @pytest.mark.unit
    @pytest.mark.parametrize("star_name", BRIGHT_STARS[:5])
    def test_fixstar_with_speed_flag(self, star_name: str):
        """fixstar_ut works with FLG_SPEED."""
        jd = 2451545.0
        result = swe.fixstar_ut(star_name, jd, FLG_SPEED)
        if isinstance(result[0], tuple):
            pos = result[0]
        elif isinstance(result[0], str):
            pos = result[1]
        else:
            pos = result[0]

        assert len(pos) == 6


class TestProperMotion:
    """Test that stars move due to proper motion over centuries."""

    @pytest.mark.unit
    @pytest.mark.parametrize("star_name,pm_arcsec_yr", HIGH_PROPER_MOTION_STARS)
    def test_star_position_changes_over_century(
        self, star_name: str, pm_arcsec_yr: float
    ):
        """Star position should change detectably over 100 years."""
        jd1 = 2451545.0  # J2000
        jd2 = jd1 + 365.25 * 100  # +100 years

        r1 = swe.fixstar_ut(star_name, jd1, 0)
        r2 = swe.fixstar_ut(star_name, jd2, 0)

        # Extract positions
        if isinstance(r1[0], str):
            pos1, pos2 = r1[1], r2[1]
        else:
            pos1, pos2 = r1[0], r2[0]

        lon1, lat1 = pos1[0], pos1[1]
        lon2, lat2 = pos2[0], pos2[1]

        # Total positional change over 100 years
        dlon = lon2 - lon1
        if dlon > 180:
            dlon -= 360
        elif dlon < -180:
            dlon += 360
        dlat = lat2 - lat1

        total_change_deg = math.sqrt(dlon**2 + dlat**2)
        total_change_arcsec = total_change_deg * 3600

        # Expected change from proper motion over 100 years
        # Plus precession (~1.4°/century = 5040"/century)
        # So total change should be at least precession
        assert total_change_arcsec > 1000, (
            f'{star_name}: total change {total_change_arcsec:.1f}" over '
            f"100 years (too small — expected precession + proper motion)"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("star_name", BRIGHT_STARS[:5])
    def test_star_position_10_year_change(self, star_name: str):
        """Star position should change over 10 years due to precession."""
        jd1 = 2451545.0
        jd2 = jd1 + 365.25 * 10

        r1 = swe.fixstar_ut(star_name, jd1, 0)
        r2 = swe.fixstar_ut(star_name, jd2, 0)

        if isinstance(r1[0], str):
            pos1, pos2 = r1[1], r2[1]
        else:
            pos1, pos2 = r1[0], r2[0]

        # Over 10 years, precession moves everything ~0.14°
        dlon = abs(pos2[0] - pos1[0])
        if dlon > 180:
            dlon = 360 - dlon

        # Should see at least 0.05° change over 10 years
        assert dlon > 0.05, f"{star_name}: only {dlon:.4f}° change over 10 years"


class TestFixstarMultipleDates:
    """Test fixed stars at multiple dates across centuries."""

    @pytest.mark.unit
    @pytest.mark.parametrize("star_name", BRIGHT_STARS)
    def test_fixstar_across_centuries(self, star_name: str):
        """Star position is valid from 1600 to 2500."""
        for year in [1600, 1700, 1800, 1900, 2000, 2100, 2200, 2400]:
            jd = swe.julday(year, 1, 1, 12.0)
            result = swe.fixstar_ut(star_name, jd, 0)
            if isinstance(result[0], str):
                pos = result[1]
            else:
                pos = result[0]

            assert 0 <= pos[0] < 360, f"{star_name} @ {year}: lon={pos[0]}"
            assert -90 <= pos[1] <= 90, f"{star_name} @ {year}: lat={pos[1]}"

    @pytest.mark.unit
    def test_star_longitude_increases_with_precession(self):
        """Star longitudes should generally increase due to precession."""
        star = "Spica"
        jd1 = swe.julday(1900, 1, 1, 12.0)
        jd2 = swe.julday(2100, 1, 1, 12.0)

        r1 = swe.fixstar_ut(star, jd1, 0)
        r2 = swe.fixstar_ut(star, jd2, 0)

        if isinstance(r1[0], str):
            pos1, pos2 = r1[1], r2[1]
        else:
            pos1, pos2 = r1[0], r2[0]

        # Precession moves longitudes forward by ~1.4°/century
        dlon = (pos2[0] - pos1[0]) % 360
        if dlon > 180:
            dlon -= 360

        # Over 200 years, should see ~2.8° increase
        assert 1.5 < dlon < 5.0, f"{star}: longitude change {dlon:.2f}° over 200 years"


class TestFixstarMagnitude:
    """Test fixed star magnitude function."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "star_name,expected_mag_range",
        [
            ("Sirius", (-2.0, -1.0)),
            ("Arcturus", (-0.5, 0.5)),
            ("Vega", (-0.5, 0.5)),
            ("Rigel", (-0.5, 0.5)),
            ("Betelgeuse", (0.0, 1.5)),
            ("Aldebaran", (0.5, 1.5)),
            ("Spica", (0.5, 1.5)),
            ("Regulus", (1.0, 2.0)),
        ],
    )
    def test_fixstar_magnitude_reasonable(
        self, star_name: str, expected_mag_range: tuple[float, float]
    ):
        """Star magnitudes should be within expected ranges."""
        result = swe.fixstar_mag(star_name)
        # fixstar_mag returns (magnitude, star_name) tuple
        mag = result[0] if isinstance(result, tuple) else result
        lo, hi = expected_mag_range
        assert lo <= mag <= hi, f"{star_name}: magnitude {mag} outside [{lo}, {hi}]"


class TestFixstarContinuity:
    """Test that star positions vary smoothly over time."""

    @pytest.mark.unit
    @pytest.mark.parametrize("star_name", ["Sirius", "Vega", "Arcturus"])
    def test_star_continuity_over_year(self, star_name: str):
        """Star position changes smoothly over 1 year (monthly samples)."""
        jd_start = 2451545.0
        prev_lon = None
        for i in range(12):
            jd = jd_start + i * 30.44
            result = swe.fixstar_ut(star_name, jd, 0)
            if isinstance(result[0], str):
                pos = result[1]
            else:
                pos = result[0]

            lon = pos[0]
            if prev_lon is not None:
                diff = abs(lon - prev_lon)
                if diff > 180:
                    diff = 360 - diff
                # Over 1 month, star should move < 0.2°
                assert diff < 0.2, f"{star_name}: jump {diff:.4f}° at month {i}"
            prev_lon = lon
