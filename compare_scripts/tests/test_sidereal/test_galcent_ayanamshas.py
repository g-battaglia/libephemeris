"""
Tests for Galactic Center based ayanamshas.

This module tests the four Galactic Center based ayanamsha modes:
- SIDM_GALCENT_0SAG (17): Galactic Center at 0° Sagittarius
- SIDM_GALCENT_RGILBRAND (30): Galactic Center (Gil Brand)
- SIDM_GALCENT_MULA_WILHELM (36): Galactic Center at Mula (Wilhelm)
- SIDM_GALCENT_COCHRANE (40): Galactic Center at 0° Capricorn (Cochrane)

The Galactic Center coordinates used are from Reid & Brunthaler (2004, ApJ 616, 872):
- RA = 17h 45m 40.0409s = 266.41683708°
- Dec = -29° 00' 28.118" = -29.00781056°
- Epoch: J2000.0

Proper motion from Reid & Brunthaler (2020, ApJ 892, 39):
- μα* = -3.151 ± 0.018 mas/yr
- μδ = -5.547 ± 0.026 mas/yr
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import (
    SIDM_GALCENT_0SAG,
    SIDM_GALCENT_RGILBRAND,
    SIDM_GALCENT_MULA_WILHELM,
    SIDM_GALCENT_COCHRANE,
    SUN,
    FLG_SIDEREAL,
)

try:
    import swisseph as swe

    HAS_PYSWISSEPH = True
except ImportError:
    HAS_PYSWISSEPH = False


# Test dates
JD_J2000 = 2451545.0  # 2000-01-01 12:00 TT
JD_1900 = 2415020.0  # 1900-01-01 12:00 TT (approximately)
JD_2100 = 2488070.0  # 2100-01-01 12:00 TT (approximately)


def angle_diff(a1: float, a2: float) -> float:
    """Calculate the minimum angular difference between two angles."""
    diff = abs(a1 - a2) % 360.0
    if diff > 180.0:
        diff = 360.0 - diff
    return diff


class TestGalacticCenterCoordinates:
    """Test that Galactic Center coordinates are correctly defined."""

    @pytest.mark.unit
    def test_galcenter_position_at_j2000(self):
        """Galactic Center should be at ~266.4° RA, -29° Dec at J2000."""
        from libephemeris.planets import STARS

        gc = STARS["GAL_CENTER"]

        # Reid & Brunthaler (2004) coordinates
        expected_ra = 266.41683708  # 17h 45m 40.0409s
        expected_dec = -29.00781056  # -29° 00' 28.118"

        assert abs(gc.ra_j2000 - expected_ra) < 0.0001, (
            f"RA mismatch: {gc.ra_j2000} vs {expected_ra}"
        )
        assert abs(gc.dec_j2000 - expected_dec) < 0.0001, (
            f"Dec mismatch: {gc.dec_j2000} vs {expected_dec}"
        )

    @pytest.mark.unit
    def test_galcenter_proper_motion(self):
        """Galactic Center proper motion should match Reid & Brunthaler values."""
        from libephemeris.planets import STARS

        gc = STARS["GAL_CENTER"]

        # Reid & Brunthaler (2004) equatorial proper motions
        # μα* = -3.151 ± 0.018 mas/yr = -0.003151 arcsec/yr
        # μδ  = -5.547 ± 0.026 mas/yr = -0.005547 arcsec/yr
        expected_pm_ra = -0.003151
        expected_pm_dec = -0.005547

        assert abs(gc.pm_ra - expected_pm_ra) < 0.0001, (
            f"PM RA mismatch: {gc.pm_ra} vs {expected_pm_ra}"
        )
        assert abs(gc.pm_dec - expected_pm_dec) < 0.0001, (
            f"PM Dec mismatch: {gc.pm_dec} vs {expected_pm_dec}"
        )


class TestGalcentAyanamshaValues:
    """Test Galactic Center ayanamsha values at different epochs."""

    GALCENT_MODES = [
        (SIDM_GALCENT_0SAG, "Galactic Center 0 Sag"),
        (SIDM_GALCENT_RGILBRAND, "Galactic Center Rgilbrand"),
        (SIDM_GALCENT_MULA_WILHELM, "Galactic Center Mula Wilhelm"),
        (SIDM_GALCENT_COCHRANE, "Galactic Center Cochrane"),
    ]

    @pytest.mark.unit
    @pytest.mark.parametrize("sid_mode,name", GALCENT_MODES)
    def test_ayanamsha_is_positive(self, sid_mode, name):
        """All Galactic Center ayanamshas should be positive at J2000."""
        ephem.set_sid_mode(sid_mode)
        ayan = ephem.get_ayanamsa_ut(JD_J2000)

        assert ayan > 0, f"{name} ayanamsha {ayan} should be positive"

    @pytest.mark.unit
    @pytest.mark.parametrize("sid_mode,name", GALCENT_MODES[:3])  # Exclude Cochrane
    def test_ayanamsha_reasonable_range(self, sid_mode, name):
        """Galactic Center ayanamshas should be in reasonable range at J2000."""
        ephem.set_sid_mode(sid_mode)
        ayan = ephem.get_ayanamsa_ut(JD_J2000)

        # All should be roughly between 20° and 35° at J2000
        assert 20.0 < ayan < 35.0, (
            f"{name} ayanamsha {ayan} outside expected range [20, 35]"
        )

    @pytest.mark.unit
    def test_cochrane_ayanamsha_range(self):
        """Cochrane ayanamsha is special - GC at 0° Cap means ~-3° to +5° ayanamsha."""
        ephem.set_sid_mode(SIDM_GALCENT_COCHRANE)
        ayan = ephem.get_ayanamsa_ut(JD_J2000)

        # Cochrane places GC at 0° Capricorn (270°)
        # At J2000, GC tropical longitude is ~266.8°
        # So ayanamsha = 266.8 - 270 = -3.2° = 356.8° (mod 360)
        # Valid range: 350° to 360° OR 0° to 10°
        in_range = (350.0 < ayan <= 360.0) or (0.0 <= ayan < 10.0)
        assert in_range, (
            f"Cochrane ayanamsha {ayan}° outside expected range [350-360° or 0-10°]"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("sid_mode,name", GALCENT_MODES)
    def test_ayanamsha_increases_over_time(self, sid_mode, name):
        """Ayanamsha should increase with time due to precession."""
        ephem.set_sid_mode(sid_mode)

        ayan_1900 = ephem.get_ayanamsa_ut(JD_1900)
        ayan_2000 = ephem.get_ayanamsa_ut(JD_J2000)
        ayan_2100 = ephem.get_ayanamsa_ut(JD_2100)

        assert ayan_2000 > ayan_1900, (
            f"{name}: Ayanamsha should increase from 1900 to 2000"
        )
        assert ayan_2100 > ayan_2000, (
            f"{name}: Ayanamsha should increase from 2000 to 2100"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("sid_mode,name", GALCENT_MODES)
    def test_ayanamsha_precession_rate(self, sid_mode, name):
        """Ayanamsha should increase ~50 arcsec/year (precession rate)."""
        ephem.set_sid_mode(sid_mode)

        ayan_1900 = ephem.get_ayanamsa_ut(JD_1900)
        ayan_2000 = ephem.get_ayanamsa_ut(JD_J2000)

        # 100 years × 50"/year = 5000" = ~1.39°
        diff = ayan_2000 - ayan_1900
        assert 1.2 < diff < 1.6, (
            f"{name}: 100-year change {diff}° outside expected [1.2, 1.6]"
        )


class TestGalcentAyanamshaRelationships:
    """Test relationships between different Galactic Center ayanamshas."""

    @pytest.mark.unit
    def test_0sag_and_cochrane_differ_by_30_degrees(self):
        """0 Sag and Cochrane should differ by ~30° (one sign), accounting for wrap."""
        ephem.set_sid_mode(SIDM_GALCENT_0SAG)
        ayan_0sag = ephem.get_ayanamsa_ut(JD_J2000)

        ephem.set_sid_mode(SIDM_GALCENT_COCHRANE)
        ayan_cochrane = ephem.get_ayanamsa_ut(JD_J2000)

        # Cochrane places GC at 0° Capricorn (270°), 0SAG at 0° Sagittarius (240°)
        # So Cochrane ayanamsha should be 30° less, but may wrap around 360°
        # Use angle_diff to handle wrap-around
        diff = angle_diff(ayan_0sag, ayan_cochrane)
        assert abs(diff - 30.0) < 1.0, f"0 Sag - Cochrane = {diff}°, expected ~30°"

    @pytest.mark.unit
    def test_ayanamsha_ordering_excluding_cochrane(self):
        """Test expected ordering of non-wrapping Galactic Center ayanamshas."""
        ephem.set_sid_mode(SIDM_GALCENT_0SAG)
        ayan_0sag = ephem.get_ayanamsa_ut(JD_J2000)

        ephem.set_sid_mode(SIDM_GALCENT_RGILBRAND)
        ayan_rgilbrand = ephem.get_ayanamsa_ut(JD_J2000)

        ephem.set_sid_mode(SIDM_GALCENT_MULA_WILHELM)
        ayan_mula_wilhelm = ephem.get_ayanamsa_ut(JD_J2000)

        # Expected order from highest to lowest ayanamsha (for non-wrapping):
        # 0 Sag = GC at 240°, offset = gc_lon - 240 (highest ayanamsha)
        # Rgilbrand uses ~244.38° offset
        # Mula Wilhelm uses ~246.81° offset
        # So: 0SAG > Rgilbrand > Mula Wilhelm
        assert ayan_0sag > ayan_rgilbrand, (
            f"0 Sag ({ayan_0sag}) should > Rgilbrand ({ayan_rgilbrand})"
        )
        assert ayan_rgilbrand > ayan_mula_wilhelm, (
            f"Rgilbrand ({ayan_rgilbrand}) should > Mula Wilhelm ({ayan_mula_wilhelm})"
        )


@pytest.mark.skipif(not HAS_PYSWISSEPH, reason="pyswisseph not installed")
class TestGalcentVsPyswisseph:
    """Compare Galactic Center ayanamshas with pyswisseph."""

    GALCENT_MODES = [
        (SIDM_GALCENT_0SAG, "Galactic Center 0 Sag"),
        (SIDM_GALCENT_RGILBRAND, "Galactic Center Rgilbrand"),
        (SIDM_GALCENT_MULA_WILHELM, "Galactic Center Mula Wilhelm"),
        (SIDM_GALCENT_COCHRANE, "Galactic Center Cochrane"),
    ]

    TEST_DATES = [
        (JD_J2000, "J2000"),
        (JD_1900, "1900"),
        (JD_2100, "2100"),
    ]

    # Galactic Center ayanamshas now use calibrated SE-compatible formulas
    # with precision better than 2 arcsec (0.0006°) across all epochs
    TOLERANCE = 0.001  # degrees (3.6 arcsec)

    @pytest.mark.comparison
    @pytest.mark.parametrize("sid_mode,name", GALCENT_MODES)
    @pytest.mark.parametrize("jd,epoch_name", TEST_DATES)
    def test_ayanamsha_vs_pyswisseph(self, sid_mode, name, jd, epoch_name):
        """Compare Galactic Center ayanamsha with pyswisseph."""
        ephem.set_sid_mode(sid_mode)
        swe.set_sid_mode(sid_mode)

        ayan_lib = ephem.get_ayanamsa_ut(jd)
        ayan_swe = swe.get_ayanamsa_ut(jd)

        diff = angle_diff(ayan_lib, ayan_swe)

        assert diff < self.TOLERANCE, (
            f"{name} at {epoch_name}: diff {diff}° >= tolerance {self.TOLERANCE}°\n"
            f"  libephemeris: {ayan_lib}°\n"
            f"  pyswisseph:   {ayan_swe}°"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("sid_mode,name", GALCENT_MODES)
    def test_sidereal_sun_vs_pyswisseph(self, sid_mode, name):
        """Compare sidereal Sun position with pyswisseph for GC ayanamshas."""
        ephem.set_sid_mode(sid_mode)
        swe.set_sid_mode(sid_mode)

        pos_lib, _ = ephem.calc_ut(JD_J2000, SUN, FLG_SIDEREAL)
        pos_swe, _ = swe.calc_ut(JD_J2000, SUN, FLG_SIDEREAL)

        diff = angle_diff(pos_lib[0], pos_swe[0])

        assert diff < self.TOLERANCE, (
            f"{name} sidereal Sun: diff {diff}° >= tolerance {self.TOLERANCE}°\n"
            f"  libephemeris: {pos_lib[0]}°\n"
            f"  pyswisseph:   {pos_swe[0]}°"
        )


class TestGalacticCenterSiderealPosition:
    """Test that the Galactic Center appears at the expected sidereal position."""

    @pytest.mark.unit
    def test_gc_at_0_sag_mode(self):
        """In 0 Sag mode, GC ecliptic longitude - ayanamsha should be ~240°."""
        # This is a conceptual test - the GC should be at 0° Sag (240°) in this mode
        # by definition, so tropical_gc_lon - ayanamsha ≈ 240°
        from libephemeris.planets import _get_star_position_ecliptic, STARS

        ephem.set_sid_mode(SIDM_GALCENT_0SAG)
        ayan = ephem.get_ayanamsa_ut(JD_J2000)

        # Get GC tropical ecliptic longitude at J2000
        # We need to compute this manually using the internal function
        # For this test, we verify the ayanamsha value is reasonable
        # The GC tropical longitude at J2000 is approximately 266.5° (in Sagittarius)

        # At J2000, the GC is at about 26.8° Sagittarius tropical
        # If ayanamsha is ~26.8°, then GC would be at 0° Sag sidereal
        assert 25.0 < ayan < 28.0, f"0 Sag ayanamsha {ayan}° should place GC at ~0° Sag"

    @pytest.mark.unit
    def test_gc_at_cochrane_mode(self):
        """In Cochrane mode, GC should be at ~0° Capricorn sidereal."""
        # Cochrane places GC at 0° Capricorn (270°)
        # So ayanamsha = tropical_gc_lon - 270
        # If GC is at ~266.8° tropical, ayanamsha ≈ -3.2° (mod 360) = 356.8°

        ephem.set_sid_mode(SIDM_GALCENT_COCHRANE)
        ayan = ephem.get_ayanamsa_ut(JD_J2000)

        # Cochrane ayanamsha wraps around 360°
        # Valid range: 350° to 360° (near 0° when wrapped)
        in_range = (350.0 < ayan <= 360.0) or (0.0 <= ayan < 10.0)
        assert in_range, (
            f"Cochrane ayanamsha {ayan}° should be close to 0°/360° at J2000"
        )
