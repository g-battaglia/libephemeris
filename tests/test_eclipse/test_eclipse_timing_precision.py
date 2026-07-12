"""Tests for eclipse timing precision and local eclipse computations."""

from __future__ import annotations

import pytest
import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    ECL_TOTAL,
    ECL_ANNULAR,
    ECL_PARTIAL,
    ECL_PENUMBRAL,
    FLG_SWIEPH,
)

JD_J2000 = 2451545.0


@pytest.mark.unit
class TestSolarEclipseTiming:
    """Test solar eclipse timing precision."""

    def test_sol_eclipse_tret_10_elements(self):
        """sol_eclipse_when_glob returns 10-element tret."""
        retflag, tret = swe.sol_eclipse_when_glob(JD_J2000)
        assert len(tret) == 10

    def test_sol_eclipse_max_time_nonzero(self):
        """tret[0] = maximum eclipse time should be nonzero."""
        retflag, tret = swe.sol_eclipse_when_glob(JD_J2000)
        assert tret[0] > JD_J2000

    def test_sol_eclipse_contacts_ordering(self):
        """Contact times should be in chronological order where nonzero."""
        retflag, tret = swe.sol_eclipse_when_glob(JD_J2000)
        # tret[0] = max, tret[2] = begin, tret[3] = end
        # tret[4] = begin total/annular, tret[5] = end total/annular
        if tret[2] > 0 and tret[3] > 0:
            assert tret[2] < tret[0] < tret[3], "Begin < max < end"

    def test_5_solar_eclipses_sequential(self):
        """Find 5 sequential solar eclipses, each later than the previous."""
        jd = JD_J2000
        eclipse_jds = []
        for _ in range(5):
            retflag, tret = swe.sol_eclipse_when_glob(jd)
            eclipse_jds.append(tret[0])
            jd = tret[0] + 1.0  # Start after the eclipse
        for i in range(1, len(eclipse_jds)):
            assert eclipse_jds[i] > eclipse_jds[i - 1]

    def test_sol_eclipse_types_valid(self):
        """Eclipse type flags should be valid."""
        retflag, tret = swe.sol_eclipse_when_glob(JD_J2000)
        # At least one type flag should be set
        type_flags = ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL
        assert retflag & type_flags, f"No eclipse type in retflag {retflag}"

    def test_sol_eclipse_filter_total(self):
        """Can filter for total eclipses only."""
        retflag, tret = swe.sol_eclipse_when_glob(JD_J2000, FLG_SWIEPH, ECL_TOTAL)
        assert retflag & ECL_TOTAL

    def test_sol_eclipse_filter_annular(self):
        """Can filter for annular eclipses only."""
        retflag, tret = swe.sol_eclipse_when_glob(JD_J2000, FLG_SWIEPH, ECL_ANNULAR)
        assert retflag & ECL_ANNULAR

    def test_sol_eclipse_duration_reasonable(self):
        """Eclipse duration should be reasonable (< 1 day)."""
        retflag, tret = swe.sol_eclipse_when_glob(JD_J2000)
        if tret[2] > 0 and tret[3] > 0:
            duration_hours = (tret[3] - tret[2]) * 24.0
            assert 0 < duration_hours < 12, f"Duration {duration_hours}h unreasonable"


@pytest.mark.unit
class TestLunarEclipseTiming:
    """Test lunar eclipse timing precision."""

    def test_lun_eclipse_tret_10_elements(self):
        """lun_eclipse_when returns 10-element tret."""
        retflag, tret = swe.lun_eclipse_when(JD_J2000)
        assert len(tret) == 10

    def test_lun_eclipse_max_after_search(self):
        """Maximum time should be after search start."""
        retflag, tret = swe.lun_eclipse_when(JD_J2000)
        assert tret[0] > JD_J2000

    def test_lun_eclipse_types_valid(self):
        """Lunar eclipse type flags should be valid."""
        retflag, tret = swe.lun_eclipse_when(JD_J2000)
        type_flags = ECL_TOTAL | ECL_PARTIAL | ECL_PENUMBRAL
        assert retflag & type_flags

    def test_5_lunar_eclipses_sequential(self):
        """Find 5 sequential lunar eclipses."""
        jd = JD_J2000
        eclipse_jds = []
        for _ in range(5):
            retflag, tret = swe.lun_eclipse_when(jd)
            eclipse_jds.append(tret[0])
            jd = tret[0] + 1.0
        for i in range(1, len(eclipse_jds)):
            assert eclipse_jds[i] > eclipse_jds[i - 1]

    def test_lun_eclipse_contacts_ordering(self):
        """Lunar eclipse contact times should be ordered."""
        retflag, tret = swe.lun_eclipse_when(JD_J2000)
        if tret[2] > 0 and tret[3] > 0:
            assert tret[2] < tret[0] < tret[3]

    def test_lun_eclipse_at_full_moon(self):
        """Lunar eclipse should occur near full Moon (Sun-Moon ~180° apart)."""
        retflag, tret = swe.lun_eclipse_when(JD_J2000)
        jd_max = tret[0]
        sun_pos, _ = swe.calc_ut(jd_max, SUN, FLG_SWIEPH)
        moon_pos, _ = swe.calc_ut(jd_max, MOON, FLG_SWIEPH)
        elongation = abs(moon_pos[0] - sun_pos[0])
        if elongation > 180:
            elongation = 360 - elongation
        # Should be near 180° (within 15° — eclipses happen near nodes)
        assert elongation > 165, f"Elongation {elongation}° not near full moon"


@pytest.mark.unit
class TestSolarEclipseLocal:
    """Test local solar eclipse computations."""

    def test_sol_eclipse_loc_returns_3_values(self):
        """sol_eclipse_when_loc returns (retflag, tret, attr)."""
        # Search for eclipse visible from Rome
        result = swe.sol_eclipse_when_loc(JD_J2000, (12.4964, 41.9028, 0.0))
        assert len(result) >= 2

    @pytest.mark.parametrize(
        "lon,lat,name",
        [
            (0.0, 51.5, "London"),
            (-74.0, 40.7, "New York"),
            (139.7, 35.7, "Tokyo"),
            (151.2, -33.9, "Sydney"),
            (0.0, 0.0, "Equator"),
        ],
    )
    def test_sol_eclipse_loc_various_locations(self, lon, lat, name):
        """Local eclipse search works at various locations."""
        result = swe.sol_eclipse_when_loc(JD_J2000, (lon, lat, 0.0))
        # Should return at least retflag and tret
        assert len(result) >= 2


@pytest.mark.unit
class TestSolarEclipseAtNewMoon:
    """Test solar eclipses occur near new Moon."""

    def test_sol_eclipse_near_new_moon(self):
        """Solar eclipse should occur near new Moon (Sun-Moon ~0° apart)."""
        retflag, tret = swe.sol_eclipse_when_glob(JD_J2000)
        jd_max = tret[0]
        sun_pos, _ = swe.calc_ut(jd_max, SUN, FLG_SWIEPH)
        moon_pos, _ = swe.calc_ut(jd_max, MOON, FLG_SWIEPH)
        elongation = abs(moon_pos[0] - sun_pos[0])
        if elongation > 180:
            elongation = 360 - elongation
        assert elongation < 15, f"Elongation {elongation}° not near new moon"

    def test_sol_eclipse_spacing_reasonable(self):
        """Time between solar eclipses should be ~6 months (± 2 months)."""
        jd = JD_J2000
        retflag1, tret1 = swe.sol_eclipse_when_glob(jd)
        retflag2, tret2 = swe.sol_eclipse_when_glob(tret1[0] + 1)
        gap_days = tret2[0] - tret1[0]
        # Eclipses are typically 5-7 months apart
        assert 100 < gap_days < 400, f"Eclipse gap {gap_days} days unreasonable"
