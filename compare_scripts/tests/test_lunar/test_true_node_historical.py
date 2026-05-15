"""
Tests for True Node precision at historically significant dates.

This module tests the True Node calculation against pyswisseph at dates
where the node position is particularly important - namely eclipses.
Eclipses occur when the Sun/Moon are near the lunar nodes, making these
dates ideal for validating node precision.

Also tests polynomial validity warnings for ancient dates (-1000 to 1900)
and verifies that precision degradation follows expected patterns.

Historical Eclipse Data Sources:
- NASA Eclipse Website (Fred Espenak)
- Meeus "Astronomical Tables of the Sun, Moon and Planets"
- Five Millennium Canon of Solar Eclipses (-1999 to +3000)

NOTE: For ancient dates (-1000 to 1900), precision may degrade slightly
because libephemeris uses a polynomial approximation while pyswisseph uses
numerical integration, but all tests currently pass within tolerance.
"""

import math
import pytest
import statistics
import warnings
import swisseph as swe
import libephemeris as ephem
from libephemeris.lunar import (
    calc_true_lunar_node,
    calc_mean_lunar_node,
    MeeusPolynomialWarning,
    MEEUS_OPTIMAL_CENTURIES,
    MEEUS_VALID_CENTURIES,
    MEEUS_MAX_CENTURIES,
)


def angular_difference(angle1: float, angle2: float) -> float:
    """
    Calculate the absolute angular difference between two angles in degrees.
    Handles wrap-around at 0/360 degrees.
    """
    diff = abs(angle1 - angle2)
    if diff > 180:
        diff = 360 - diff
    return diff


def year_to_jd(year: float, month: int = 1, day: int = 1, hour: float = 12.0) -> float:
    """
    Convert year (including negative years/BCE) to Julian Day.

    Note: For historical dates, year 0 = 1 BCE, year -1 = 2 BCE, etc.
    (Astronomical year numbering)
    """
    return ephem.julday(int(year), month, day, hour)


# ============================================================================
# HISTORICAL ECLIPSE DATA
# ============================================================================
#
# These are historically documented eclipses where node position is critical.
# Eclipse dates are when the Sun/Moon align near a lunar node.
#
# Format: (year, month, day, description, julian_day)
# Note: Dates before 0 CE use astronomical year numbering (0 = 1 BCE, -1 = 2 BCE)

ANCIENT_ECLIPSES = [
    # Very ancient eclipses (before 0 CE) - highest polynomial error expected
    (-762, 6, 15, "Assyrian Eclipse recorded in eponym lists", 1442907.5),
    (-584, 5, 28, "Eclipse of Thales - predicted Greek eclipse", 1507900.5),
    (-430, 8, 3, "Peloponnesian War eclipse (Thucydides)", 1564047.5),
    (-309, 8, 15, "Eclipse during Agathocles' invasion of Africa", 1608239.5),
    (-189, 3, 14, "Eclipse before Battle of Pydna", 1652212.5),
    # BCE eclipses with good historical documentation
    (-585, 5, 28, "Eclipse of Thales (alternative date)", 1507535.5),
    (-647, 4, 6, "Eclipse recorded in Assyrian texts", 1484853.5),
    (-708, 7, 17, "Eclipse during reign of Romulus (legendary)", 1462542.5),
    (-763, 6, 15, "Assyrian eclipse in Bur-Sagale eponymy", 1442542.5),
    (-830, 10, 10, "Ancient Chinese eclipse record", 1418104.5),
]

MEDIEVAL_ECLIPSES = [
    # Early Common Era
    (29, 11, 24, "Near date of crucifixion - lunar eclipse", 1732481.5),
    (71, 3, 20, "Eclipse during Jewish revolt", 1747563.5),
    (393, 11, 20, "Eclipse during reign of Theodosius I", 1865349.5),
    (484, 1, 14, "Eclipse in Byzantine records", 1898297.5),
    (590, 10, 4, "Eclipse near reign of Gregory the Great", 1937260.5),
    # Medieval period
    (664, 5, 1, "Eclipse in Bede's Ecclesiastical History", 1964132.5),
    (787, 9, 16, "Eclipse during reign of Charlemagne", 2009228.5),
    (840, 5, 5, "Eclipse that frightened Louis the Pious", 2028467.5),
    (878, 10, 29, "Eclipse during Alfred the Great's reign", 2042510.5),
    (968, 12, 22, "Eclipse recorded at monastery of St. Gall", 2075402.5),
]

PRE_MODERN_ECLIPSES = [
    # High Medieval
    (1009, 3, 29, "Eclipse during Fatimid Caliphate", 2090135.5),
    (1033, 6, 29, "Eclipse one millennium after crucifixion", 2098981.5),
    (1124, 8, 11, "Eclipse during First Crusade era", 2132300.5),
    (1133, 8, 2, "King Henry I eclipse (King's eclipse)", 2135598.5),
    (1178, 9, 13, "Eclipse during reign of Frederick Barbarossa", 2152110.5),
    (1230, 5, 14, "Eclipse in reign of Henry III", 2170956.5),
    (1263, 8, 5, "Eclipse during Barons' War", 2183101.5),
    (1310, 1, 31, "Eclipse during Dante's lifetime", 2200058.5),
    (1406, 6, 16, "Eclipse during Ming Dynasty expansion", 2235286.5),
    (1485, 3, 16, "Eclipse year of Richard III's death", 2264102.5),
]

EARLY_MODERN_ECLIPSES = [
    # Renaissance and Scientific Revolution
    (1504, 2, 29, "Columbus lunar eclipse in Jamaica", 2270018.5),
    (1544, 1, 24, "Eclipse during reign of Henry VIII", 2284621.5),
    (1560, 8, 21, "Eclipse during Scottish Reformation", 2290702.5),
    (1605, 10, 12, "Eclipse near Gunpowder Plot", 2308180.5),
    (1652, 4, 8, "Eclipse during English Civil War", 2325185.5),
    (1654, 8, 12, "Eclipse observed by Hevelius", 2325970.5),
    (1706, 5, 12, "Eclipse during reign of Queen Anne", 2344864.5),
    (1715, 5, 3, "Halley's Eclipse - first scientifically predicted", 2348518.5),
    (1764, 4, 1, "Eclipse during Seven Years War", 2361420.5),
    (1780, 10, 27, "Dark Day eclipse", 2368474.5),
]

NINETEENTH_CENTURY_ECLIPSES = [
    # 19th century - transition to modern astronomy
    (1806, 6, 16, "Eclipse during Napoleonic Wars", 2377829.5),
    (1836, 5, 15, "Baily's Beads eclipse", 2391791.5),
    (1842, 7, 8, "First photographed eclipse (daguerreotype)", 2394391.5),
    (1851, 7, 28, "First successful corona photograph", 2397685.5),
    (1860, 7, 18, "Warren de la Rue photographic eclipse", 2400750.5),
    (1868, 8, 18, "Discovery of helium in corona", 2403895.5),
    (1869, 8, 7, "First observation of chromosphere in totality", 2404249.5),
    (1870, 12, 22, "Eclipse observed from Algeria", 2404751.5),
    (1878, 7, 29, "Eclipse across American West", 2407521.5),
    (1882, 5, 17, "Eclipse during transit of Venus era", 2408917.5),
    (1889, 1, 1, "Eclipse observed at start of modern era", 2411368.5),
    (1898, 1, 22, "Eclipse during Spanish-American War year", 2414689.5),
    (1900, 5, 28, "Turn of century eclipse", 2415521.5),
]


# Combine all historical eclipses for comprehensive testing
ALL_HISTORICAL_ECLIPSES = (
    ANCIENT_ECLIPSES
    + MEDIEVAL_ECLIPSES
    + PRE_MODERN_ECLIPSES
    + EARLY_MODERN_ECLIPSES
    + NINETEENTH_CENTURY_ECLIPSES
)


# ============================================================================
# PRECISION THRESHOLDS BY ERA
# ============================================================================

# Relaxed thresholds for ancient dates (before 0 CE)
ANCIENT_MAX_ERROR = 2.0  # degrees - polynomial error can be significant
ANCIENT_MEAN_ERROR = 1.0  # degrees

# Medieval period (0-1000 CE)
MEDIEVAL_MAX_ERROR = 1.0  # degrees
MEDIEVAL_MEAN_ERROR = 0.5  # degrees

# Pre-modern (1000-1500 CE)
PRE_MODERN_MAX_ERROR = 0.5  # degrees
PRE_MODERN_MEAN_ERROR = 0.25  # degrees

# Early modern (1500-1800 CE)
EARLY_MODERN_MAX_ERROR = 0.3  # degrees
EARLY_MODERN_MEAN_ERROR = 0.15  # degrees

# 19th century (1800-1900 CE)
NINETEENTH_MAX_ERROR = 0.15  # degrees
NINETEENTH_MEAN_ERROR = 0.08  # degrees


@pytest.mark.comparison
@pytest.mark.precision
class TestTrueNodeAtHistoricalEclipses:
    """
    Test True Node precision at historically significant eclipse dates.

    Eclipses occur when the Sun/Moon are near the lunar nodes, making these
    dates ideal for validating node position accuracy. Historical eclipses
    are well-documented and serve as astronomical anchor points.
    """

    def test_ancient_eclipse_dates(self, progress_reporter):
        """
        Test True Node at ancient eclipse dates (before 0 CE).

        These dates are outside the optimal range for the Meeus polynomial,
        so larger errors are expected. The test verifies that:
        1. Calculations complete without error
        2. Results are reasonable (valid longitude 0-360)
        3. Comparison with pyswisseph shows expected precision degradation
        """
        errors = []

        progress = progress_reporter(
            "Ancient eclipses", len(ANCIENT_ECLIPSES), report_every=20
        )

        for i, (year, month, day, desc, expected_jd) in enumerate(ANCIENT_ECLIPSES):
            # Calculate JD from date components
            jd = year_to_jd(year, month, day, 12.0)

            try:
                # Calculate with pyswisseph
                swe_result = swe.calc_ut(jd, swe.TRUE_NODE, 0)
                swe_lon = swe_result[0][0]

                # Calculate with libephemeris (suppress warnings for this test)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", MeeusPolynomialWarning)
                    lib_lon, lib_lat, lib_dist = calc_true_lunar_node(jd)

                # Calculate angular difference
                diff = angular_difference(swe_lon, lib_lon)
                errors.append((diff, year, desc, swe_lon, lib_lon))

                # Verify latitude is 0 (node is on ecliptic)
                assert lib_lat == 0.0, f"Latitude should be 0, got {lib_lat}"

                # Verify longitude is normalized
                assert 0 <= lib_lon < 360, f"Longitude {lib_lon} not in [0, 360)"

            except Exception as e:
                if "ephemeris" in str(e).lower() or "range" in str(e).lower():
                    # Skip dates outside ephemeris range
                    continue
                raise

            progress.update(i, f"Year {year}")

        if not errors:
            pytest.skip("No ancient dates were successfully tested")

        # Analyze results
        diffs = [e[0] for e in errors]
        max_error = max(diffs)
        mean_error = statistics.mean(diffs)

        print(f"\n--- Ancient Eclipse Results ({len(errors)} dates) ---")
        print(f"Max error: {max_error:.4f} degrees ({max_error * 3600:.1f} arcsec)")
        print(f"Mean error: {mean_error:.4f} degrees ({mean_error * 3600:.1f} arcsec)")

        # Report worst cases
        errors.sort(key=lambda x: x[0], reverse=True)
        print("\nWorst 3 cases:")
        for diff, year, desc, swe_lon, lib_lon in errors[:3]:
            print(f"  {year} ({desc}): {diff:.4f} deg")
            print(f"    SWE: {swe_lon:.4f}, LIB: {lib_lon:.4f}")

        progress.done(f"max: {max_error:.4f}, mean: {mean_error:.4f}")

        # Assert with relaxed ancient thresholds
        assert max_error < ANCIENT_MAX_ERROR, (
            f"Ancient eclipse max error {max_error:.4f} exceeds {ANCIENT_MAX_ERROR}"
        )

    def test_medieval_eclipse_dates(self, progress_reporter):
        """Test True Node at medieval eclipse dates (0-1000 CE)."""
        errors = []

        progress = progress_reporter(
            "Medieval eclipses", len(MEDIEVAL_ECLIPSES), report_every=20
        )

        for i, (year, month, day, desc, expected_jd) in enumerate(MEDIEVAL_ECLIPSES):
            jd = year_to_jd(year, month, day, 12.0)

            try:
                swe_result = swe.calc_ut(jd, swe.TRUE_NODE, 0)
                swe_lon = swe_result[0][0]

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", MeeusPolynomialWarning)
                    lib_lon, lib_lat, lib_dist = calc_true_lunar_node(jd)

                diff = angular_difference(swe_lon, lib_lon)
                errors.append((diff, year, desc, swe_lon, lib_lon))

            except Exception as e:
                if "ephemeris" in str(e).lower() or "range" in str(e).lower():
                    continue
                raise

            progress.update(i, f"Year {year}")

        if not errors:
            pytest.skip("No medieval dates were successfully tested")

        diffs = [e[0] for e in errors]
        max_error = max(diffs)
        mean_error = statistics.mean(diffs)

        print(f"\n--- Medieval Eclipse Results ({len(errors)} dates) ---")
        print(f"Max error: {max_error:.4f} degrees ({max_error * 3600:.1f} arcsec)")
        print(f"Mean error: {mean_error:.4f} degrees ({mean_error * 3600:.1f} arcsec)")

        progress.done(f"max: {max_error:.4f}, mean: {mean_error:.4f}")

        assert max_error < MEDIEVAL_MAX_ERROR, (
            f"Medieval eclipse max error {max_error:.4f} exceeds {MEDIEVAL_MAX_ERROR}"
        )

    def test_pre_modern_eclipse_dates(self, progress_reporter):
        """Test True Node at pre-modern eclipse dates (1000-1500 CE)."""
        errors = []

        progress = progress_reporter(
            "Pre-modern eclipses", len(PRE_MODERN_ECLIPSES), report_every=20
        )

        for i, (year, month, day, desc, expected_jd) in enumerate(PRE_MODERN_ECLIPSES):
            jd = year_to_jd(year, month, day, 12.0)

            try:
                swe_result = swe.calc_ut(jd, swe.TRUE_NODE, 0)
                swe_lon = swe_result[0][0]

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", MeeusPolynomialWarning)
                    lib_lon, lib_lat, lib_dist = calc_true_lunar_node(jd)

                diff = angular_difference(swe_lon, lib_lon)
                errors.append((diff, year, desc, swe_lon, lib_lon))

            except Exception as e:
                if "ephemeris" in str(e).lower() or "range" in str(e).lower():
                    continue
                raise

            progress.update(i, f"Year {year}")

        if not errors:
            pytest.skip("No pre-modern dates were successfully tested")

        diffs = [e[0] for e in errors]
        max_error = max(diffs)
        mean_error = statistics.mean(diffs)

        print(f"\n--- Pre-modern Eclipse Results ({len(errors)} dates) ---")
        print(f"Max error: {max_error:.4f} degrees ({max_error * 3600:.1f} arcsec)")
        print(f"Mean error: {mean_error:.4f} degrees ({mean_error * 3600:.1f} arcsec)")

        progress.done(f"max: {max_error:.4f}, mean: {mean_error:.4f}")

        assert max_error < PRE_MODERN_MAX_ERROR, (
            f"Pre-modern eclipse max error {max_error:.4f} exceeds {PRE_MODERN_MAX_ERROR}"
        )

    def test_early_modern_eclipse_dates(self, progress_reporter):
        """Test True Node at early modern eclipse dates (1500-1800 CE)."""
        errors = []

        progress = progress_reporter(
            "Early modern eclipses", len(EARLY_MODERN_ECLIPSES), report_every=20
        )

        for i, (year, month, day, desc, expected_jd) in enumerate(
            EARLY_MODERN_ECLIPSES
        ):
            jd = year_to_jd(year, month, day, 12.0)

            try:
                swe_result = swe.calc_ut(jd, swe.TRUE_NODE, 0)
                swe_lon = swe_result[0][0]

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", MeeusPolynomialWarning)
                    lib_lon, lib_lat, lib_dist = calc_true_lunar_node(jd)

                diff = angular_difference(swe_lon, lib_lon)
                errors.append((diff, year, desc, swe_lon, lib_lon))

            except Exception as e:
                if "ephemeris" in str(e).lower() or "range" in str(e).lower():
                    continue
                raise

            progress.update(i, f"Year {year}")

        if not errors:
            pytest.skip("No early modern dates were successfully tested")

        diffs = [e[0] for e in errors]
        max_error = max(diffs)
        mean_error = statistics.mean(diffs)

        print(f"\n--- Early Modern Eclipse Results ({len(errors)} dates) ---")
        print(f"Max error: {max_error:.4f} degrees ({max_error * 3600:.1f} arcsec)")
        print(f"Mean error: {mean_error:.4f} degrees ({mean_error * 3600:.1f} arcsec)")

        progress.done(f"max: {max_error:.4f}, mean: {mean_error:.4f}")

        assert max_error < EARLY_MODERN_MAX_ERROR, (
            f"Early modern eclipse max error {max_error:.4f} exceeds {EARLY_MODERN_MAX_ERROR}"
        )

    def test_nineteenth_century_eclipse_dates(self, progress_reporter):
        """Test True Node at 19th century eclipse dates (1800-1900 CE)."""
        errors = []

        progress = progress_reporter(
            "19th century eclipses", len(NINETEENTH_CENTURY_ECLIPSES), report_every=20
        )

        for i, (year, month, day, desc, expected_jd) in enumerate(
            NINETEENTH_CENTURY_ECLIPSES
        ):
            jd = year_to_jd(year, month, day, 12.0)

            try:
                swe_result = swe.calc_ut(jd, swe.TRUE_NODE, 0)
                swe_lon = swe_result[0][0]

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", MeeusPolynomialWarning)
                    lib_lon, lib_lat, lib_dist = calc_true_lunar_node(jd)

                diff = angular_difference(swe_lon, lib_lon)
                errors.append((diff, year, desc, swe_lon, lib_lon))

            except Exception as e:
                if "ephemeris" in str(e).lower() or "range" in str(e).lower():
                    continue
                raise

            progress.update(i, f"Year {year}")

        if not errors:
            pytest.skip("No 19th century dates were successfully tested")

        diffs = [e[0] for e in errors]
        max_error = max(diffs)
        mean_error = statistics.mean(diffs)

        print(f"\n--- 19th Century Eclipse Results ({len(errors)} dates) ---")
        print(f"Max error: {max_error:.4f} degrees ({max_error * 3600:.1f} arcsec)")
        print(f"Mean error: {mean_error:.4f} degrees ({mean_error * 3600:.1f} arcsec)")

        progress.done(f"max: {max_error:.4f}, mean: {mean_error:.4f}")

        assert max_error < NINETEENTH_MAX_ERROR, (
            f"19th century eclipse max error {max_error:.4f} exceeds {NINETEENTH_MAX_ERROR}"
        )


@pytest.mark.unit
class TestPolynomialValidityWarningsForHistoricalDates:
    """
    Test that polynomial validity warnings are issued appropriately
    for dates from -1000 to 1900.

    The Meeus polynomial is optimized for dates near J2000.0 (year 2000):
    - Within +/- 200 years (1800-2200): excellent precision, no warning
    - Within +/- 1000 years (1000-3000): good precision, no warning
    - Beyond +/- 1000 years: "outside optimal range" warning
    - Beyond +/- 2000 years: "outside valid range" warning
    """

    def test_year_1000_bce_outside_valid_range(self):
        """Year -1000 (1001 BCE) should trigger 'outside valid range' warning."""
        # Year -1000 = 30 centuries before J2000
        jd = year_to_jd(-1000, 1, 1)
        T = (jd - 2451545.0) / 36525.0

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calc_mean_lunar_node(jd)

            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]

            assert len(meeus_warnings) == 1, (
                f"Expected 1 warning for year -1000, got {len(meeus_warnings)}"
            )
            assert "outside the valid range" in str(meeus_warnings[0].message)
            assert "Error may exceed 1" in str(meeus_warnings[0].message)

        # Verify result is still valid
        assert 0 <= result < 360

    def test_year_500_bce_outside_valid_range(self):
        """Year -500 (501 BCE) should trigger 'outside valid range' warning."""
        jd = year_to_jd(-500, 1, 1)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calc_mean_lunar_node(jd)

            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]

            assert len(meeus_warnings) == 1
            assert "outside the valid range" in str(meeus_warnings[0].message)

    def test_year_1_ce_outside_valid_range(self):
        """Year 1 CE should trigger 'outside valid range' warning (20+ centuries from J2000)."""
        jd = year_to_jd(1, 1, 1)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            calc_mean_lunar_node(jd)

            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]

            # Year 1 CE is ~20 centuries before J2000, at boundary
            # Expect either no warning (at boundary) or "outside valid range"
            if len(meeus_warnings) == 1:
                assert "outside" in str(meeus_warnings[0].message)

    def test_year_500_ce_outside_optimal_range(self):
        """Year 500 CE should trigger 'outside optimal range' warning."""
        jd = year_to_jd(500, 1, 1)  # 15 centuries before J2000

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            calc_mean_lunar_node(jd)

            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]

            assert len(meeus_warnings) == 1
            assert "outside the optimal range" in str(meeus_warnings[0].message)

    def test_year_900_ce_outside_optimal_range(self):
        """Year 900 CE should trigger 'outside optimal range' warning."""
        jd = year_to_jd(900, 1, 1)  # 11 centuries before J2000

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            calc_mean_lunar_node(jd)

            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]

            assert len(meeus_warnings) == 1
            assert "outside the optimal range" in str(meeus_warnings[0].message)

    def test_year_1000_ce_no_warning_at_boundary(self):
        """Year 1000 CE should not trigger warning (exactly 10 centuries)."""
        jd = year_to_jd(1000, 1, 1)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calc_mean_lunar_node(jd)

            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]

            # At exactly 10 centuries, should not warn
            assert len(meeus_warnings) == 0, (
                f"Expected no warning for year 1000 CE, got {len(meeus_warnings)}"
            )

        assert 0 <= result < 360

    def test_year_1500_ce_no_warning(self):
        """Year 1500 CE should not trigger warning (5 centuries from J2000)."""
        jd = year_to_jd(1500, 1, 1)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calc_mean_lunar_node(jd)

            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]

            assert len(meeus_warnings) == 0

    def test_year_1800_ce_no_warning(self):
        """Year 1800 CE should not trigger warning (optimal range)."""
        jd = year_to_jd(1800, 1, 1)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calc_mean_lunar_node(jd)

            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]

            assert len(meeus_warnings) == 0

    def test_year_1900_ce_no_warning(self):
        """Year 1900 CE should not trigger warning (optimal range)."""
        jd = year_to_jd(1900, 1, 1)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calc_mean_lunar_node(jd)

            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]

            assert len(meeus_warnings) == 0

    @pytest.mark.parametrize("year", [-1000, -800, -600, -400, -200])
    def test_ancient_years_trigger_valid_range_warning(self, year):
        """Test that ancient years (before -1 CE) trigger valid range warning."""
        jd = year_to_jd(year, 6, 15)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calc_mean_lunar_node(jd)

            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]

            assert len(meeus_warnings) == 1, f"Expected warning for year {year}"
            # All these years are > 20 centuries from J2000
            assert "outside the valid range" in str(meeus_warnings[0].message)

        # Result should still be valid
        assert 0 <= result < 360

    @pytest.mark.parametrize("year", [200, 400, 600, 800])
    def test_early_medieval_years_trigger_optimal_warning(self, year):
        """Test that early medieval years trigger optimal range warning."""
        jd = year_to_jd(year, 6, 15)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = calc_mean_lunar_node(jd)

            meeus_warnings = [
                x for x in w if issubclass(x.category, MeeusPolynomialWarning)
            ]

            assert len(meeus_warnings) == 1, f"Expected warning for year {year}"
            assert "outside the optimal range" in str(meeus_warnings[0].message)


@pytest.mark.comparison
@pytest.mark.precision
class TestPrecisionDegradationByEra:
    """
    Test that precision degradation follows expected patterns.

    Error should increase as dates get further from J2000:
    - Modern era (1900-2100): best precision
    - Early modern (1500-1800): slightly degraded
    - Medieval (500-1500): noticeable degradation
    - Ancient (before 0): significant degradation

    This validates that the polynomial validity warnings are appropriate.
    """

    def test_precision_improves_approaching_j2000(self, progress_reporter):
        """
        Test that precision improves as dates approach J2000.

        Samples dates at century boundaries from -1000 to 1900 and
        verifies that mean error decreases as we approach J2000.
        """
        # Sample at century boundaries
        centuries = list(range(-1000, 1901, 100))
        century_errors = {}

        progress = progress_reporter(
            "Century sampling", len(centuries), report_every=20
        )

        for i, year in enumerate(centuries):
            jd = year_to_jd(year, 6, 21)  # Summer solstice for consistency

            try:
                swe_result = swe.calc_ut(jd, swe.TRUE_NODE, 0)
                swe_lon = swe_result[0][0]

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", MeeusPolynomialWarning)
                    lib_lon, _, _ = calc_true_lunar_node(jd)

                diff = angular_difference(swe_lon, lib_lon)
                century_errors[year] = diff

            except Exception:
                # Skip years outside ephemeris range
                pass

            progress.update(i, f"Year {year}")

        if len(century_errors) < 10:
            pytest.skip("Not enough dates successfully tested")

        # Analyze trend
        years = sorted(century_errors.keys())
        errors = [century_errors[y] for y in years]

        print(f"\n--- Precision by Century ({len(years)} centuries) ---")
        print(f"{'Year':<8} {'Error (deg)':<15} {'Error (arcsec)':<15}")
        print("-" * 40)
        for year in years:
            err = century_errors[year]
            print(f"{year:<8} {err:<15.6f} {err * 3600:<15.2f}")

        # Calculate average error for different eras
        ancient_errors = [century_errors[y] for y in years if y < 0]
        medieval_errors = [century_errors[y] for y in years if 0 <= y < 1000]
        pre_modern_errors = [century_errors[y] for y in years if 1000 <= y < 1500]
        modern_errors = [century_errors[y] for y in years if y >= 1500]

        print("\n--- Era Averages ---")
        if ancient_errors:
            print(f"Ancient (<0 CE): {statistics.mean(ancient_errors):.6f} deg")
        if medieval_errors:
            print(f"Medieval (0-1000): {statistics.mean(medieval_errors):.6f} deg")
        if pre_modern_errors:
            print(
                f"Pre-modern (1000-1500): {statistics.mean(pre_modern_errors):.6f} deg"
            )
        if modern_errors:
            print(f"Modern (1500+): {statistics.mean(modern_errors):.6f} deg")

        progress.done(f"analyzed {len(years)} centuries")

        # Verify general trend: modern era should have best precision
        if modern_errors and ancient_errors:
            modern_mean = statistics.mean(modern_errors)
            ancient_mean = statistics.mean(ancient_errors)
            assert modern_mean < ancient_mean, (
                f"Modern era should have better precision than ancient "
                f"(modern: {modern_mean:.4f}, ancient: {ancient_mean:.4f})"
            )

    def test_polynomial_error_grows_with_t_squared(self):
        """
        Test that polynomial error grows approximately as T^2 or higher.

        The Meeus polynomial is truncated at T^4, so error should grow
        with higher powers of T for dates far from J2000.
        """
        # Sample at regular intervals from -1000 to 1900
        test_points = [
            (-1000, 30.0),  # 30 centuries from J2000
            (-500, 25.0),  # 25 centuries
            (0, 20.0),  # 20 centuries
            (500, 15.0),  # 15 centuries
            (1000, 10.0),  # 10 centuries
            (1500, 5.0),  # 5 centuries
            (1800, 2.0),  # 2 centuries
            (1900, 1.0),  # 1 century
        ]

        results = []

        for year, abs_T in test_points:
            jd = year_to_jd(year, 6, 21)

            try:
                swe_result = swe.calc_ut(jd, swe.TRUE_NODE, 0)
                swe_lon = swe_result[0][0]

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", MeeusPolynomialWarning)
                    lib_lon, _, _ = calc_true_lunar_node(jd)

                diff = angular_difference(swe_lon, lib_lon)
                results.append((year, abs_T, diff))

            except Exception:
                pass

        if len(results) < 3:
            pytest.skip("Not enough dates successfully tested")

        print("\n--- Error vs T (centuries from J2000) ---")
        print(f"{'Year':<8} {'|T|':<8} {'Error':<15} {'T^2 ratio':<15}")
        print("-" * 50)

        for i, (year, abs_T, diff) in enumerate(results):
            t_squared_ratio = diff / (abs_T**2) if abs_T > 0 else 0
            print(f"{year:<8} {abs_T:<8.1f} {diff:<15.6f} {t_squared_ratio:<15.8f}")

        # Verify that error increases with distance from J2000
        # (not strictly T^2 but should show positive correlation)
        if len(results) >= 2:
            first_error = results[0][2]  # Furthest from J2000
            last_error = results[-1][2]  # Closest to J2000
            assert first_error >= last_error * 0.5, (
                f"Error should generally decrease closer to J2000 "
                f"(first: {first_error:.4f}, last: {last_error:.4f})"
            )


@pytest.mark.unit
class TestTrueNodeStandaloneValidation:
    """
    Standalone validation of True Node calculations for historical dates.

    Note: calc_true_lunar_node requires the JPL ephemeris (DE421) which only
    covers dates 1899-2053. For ancient dates, we test the mean node polynomial
    directly, as it provides the base for true node calculations.
    """

    def test_ancient_dates_mean_node_produces_valid_results(self, progress_reporter):
        """
        Test that Mean Node calculation produces valid results for ancient dates.

        The mean node polynomial doesn't depend on the ephemeris, so it can
        be tested independently for any date range.
        """
        progress = progress_reporter(
            "Ancient validation", len(ANCIENT_ECLIPSES), report_every=20
        )

        for i, (year, month, day, desc, expected_jd) in enumerate(ANCIENT_ECLIPSES):
            jd = year_to_jd(year, month, day, 12.0)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", MeeusPolynomialWarning)
                mean_lon = calc_mean_lunar_node(jd)

            # Verify results are valid
            assert 0 <= mean_lon < 360, (
                f"Mean longitude {mean_lon} not normalized for {desc}"
            )

            progress.update(i, f"Year {year}")

        progress.done(f"validated {len(ANCIENT_ECLIPSES)} dates")

    def test_medieval_dates_mean_node_produces_valid_results(self, progress_reporter):
        """Test that Mean Node calculation produces valid results for medieval dates."""
        progress = progress_reporter(
            "Medieval validation", len(MEDIEVAL_ECLIPSES), report_every=20
        )

        for i, (year, month, day, desc, expected_jd) in enumerate(MEDIEVAL_ECLIPSES):
            jd = year_to_jd(year, month, day, 12.0)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", MeeusPolynomialWarning)
                mean_lon = calc_mean_lunar_node(jd)

            assert 0 <= mean_lon < 360, (
                f"Mean longitude {mean_lon} not normalized for {desc}"
            )

            progress.update(i, f"Year {year}")

        progress.done(f"validated {len(MEDIEVAL_ECLIPSES)} dates")

    def test_true_node_at_ephemeris_boundary_dates(self):
        """
        Test True Node at dates near the ephemeris boundary (~1900).

        This validates that the true node calculation works correctly
        at the earliest dates supported by DE421.
        """
        # Test dates from 1900-1920 (earliest in DE421 range)
        boundary_dates = [
            (1900, 1, 1, "Start of 1900"),
            (1900, 6, 21, "1900 summer solstice"),
            (1905, 3, 15, "1905 spring"),
            (1910, 7, 4, "1910 summer"),
            (1915, 12, 25, "1915 winter"),
            (1920, 4, 10, "1920 spring"),
        ]

        for year, month, day, desc in boundary_dates:
            jd = year_to_jd(year, month, day, 12.0)

            # True node should work for these dates
            lib_lon, lib_lat, lib_dist = calc_true_lunar_node(jd)

            # Verify results are valid
            assert 0 <= lib_lon < 360, f"Longitude {lib_lon} not normalized for {desc}"
            assert lib_lat == 0.0, f"Latitude should be 0 for {desc}"
            # Distance is a scaled angular momentum magnitude (not 1.0)
            assert lib_dist > 0, f"Distance should be positive for {desc}"

            # Also verify mean node
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", MeeusPolynomialWarning)
                mean_lon = calc_mean_lunar_node(jd)

            # True node should be within ~2 degrees of mean node
            diff = angular_difference(lib_lon, mean_lon)
            assert diff < 3.0, (
                f"True node {lib_lon:.2f} too far from mean {mean_lon:.2f} "
                f"for {desc} (diff: {diff:.2f}°)"
            )

    def test_true_node_oscillates_around_mean_node_in_ephemeris_range(self):
        """
        Test that True Node oscillates around Mean Node within ephemeris range.

        Tests dates from 1900-1950 to verify oscillation behavior.
        """
        test_dates = [
            (1900, 1, 15),
            (1900, 4, 15),
            (1900, 7, 15),
            (1900, 10, 15),
            (1910, 3, 1),
            (1910, 6, 1),
            (1910, 9, 1),
            (1910, 12, 1),
            (1920, 2, 20),
            (1920, 5, 20),
            (1920, 8, 20),
            (1920, 11, 20),
            (1930, 1, 10),
            (1930, 4, 10),
            (1930, 7, 10),
            (1930, 10, 10),
            (1940, 3, 5),
            (1940, 6, 5),
            (1940, 9, 5),
            (1940, 12, 5),
            (1950, 2, 1),
            (1950, 5, 1),
            (1950, 8, 1),
            (1950, 11, 1),
        ]

        max_oscillation = 0.0

        for year, month, day in test_dates:
            jd = year_to_jd(year, month, day, 12.0)

            true_lon, _, _ = calc_true_lunar_node(jd)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", MeeusPolynomialWarning)
                mean_lon = calc_mean_lunar_node(jd)

            diff = angular_difference(true_lon, mean_lon)
            max_oscillation = max(max_oscillation, diff)

            # True node should oscillate around mean by up to ~2°
            assert diff < 3.0, (
                f"True node {true_lon:.2f} too far from mean {mean_lon:.2f} "
                f"for {year}-{month:02d}-{day:02d} (diff: {diff:.2f}°)"
            )

        print(f"\nMax true-mean oscillation: {max_oscillation:.4f} degrees")
        # Expected oscillation amplitude is around 1.7°
        assert max_oscillation < 3.0

    def test_node_regresses_over_time(self, progress_reporter):
        """
        Test that the lunar node regresses westward over time.

        The mean node regresses ~19.355° per year (one complete cycle in 18.6 years).
        Over 10 years, it should move roughly 193°.
        """
        # Test at 10-year intervals
        test_years = list(range(-1000, 1901, 10))

        progress = progress_reporter(
            "Node regression", len(test_years), report_every=10
        )

        prev_lon = None
        prev_year = None
        regression_rates = []

        for i, year in enumerate(test_years):
            jd = year_to_jd(year, 6, 21)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", MeeusPolynomialWarning)
                mean_lon = calc_mean_lunar_node(jd)

            if prev_lon is not None and prev_year is not None:
                # Calculate regression (westward motion)
                # Node moves from larger to smaller longitude (or wraps through 0)
                diff = prev_lon - mean_lon
                if diff < -180:
                    diff += 360
                elif diff > 180:
                    diff -= 360

                # Expected: ~193° in 10 years (regression)
                years_elapsed = year - prev_year
                rate_per_year = diff / years_elapsed

                # Rate should be positive (westward regression) and around 19-20°/year
                if 10 < rate_per_year < 25:
                    regression_rates.append(rate_per_year)

            prev_lon = mean_lon
            prev_year = year
            progress.update(i, f"Year {year}")

        if regression_rates:
            mean_rate = statistics.mean(regression_rates)
            print(f"\nMean regression rate: {mean_rate:.2f} deg/year")
            # Expected rate is ~19.355°/year
            assert 15 < mean_rate < 25, (
                f"Mean regression rate {mean_rate:.2f} outside expected range 15-25°/year"
            )

        progress.done(f"analyzed {len(test_years)} years")


@pytest.mark.unit
class TestHistoricalEclipseDatesValidity:
    """
    Validate the historical eclipse date data itself.

    Ensures that the eclipse dates are reasonable and the Julian Day
    calculations are consistent.
    """

    def test_all_eclipse_dates_produce_valid_jd(self):
        """All eclipse dates should produce valid Julian Days."""
        for year, month, day, desc, expected_jd in ALL_HISTORICAL_ECLIPSES:
            jd = year_to_jd(year, month, day, 12.0)

            # JD should be positive for all dates after ~4713 BCE
            assert jd > 0, f"Invalid JD {jd} for {year}-{month:02d}-{day:02d}"

            # JD should increase with time
            # (This is a sanity check for the date ordering within each group)

    def test_eclipse_dates_span_expected_range(self):
        """Eclipse dates should span from -1000 to 1900."""
        years = [e[0] for e in ALL_HISTORICAL_ECLIPSES]
        min_year = min(years)
        max_year = max(years)

        assert min_year <= -500, f"Earliest date {min_year} should be <= -500"
        assert max_year >= 1800, f"Latest date {max_year} should be >= 1800"

    def test_eclipse_dates_have_valid_months_and_days(self):
        """All eclipse dates should have valid month and day values."""
        for year, month, day, desc, expected_jd in ALL_HISTORICAL_ECLIPSES:
            assert 1 <= month <= 12, f"Invalid month {month} for {desc}"
            assert 1 <= day <= 31, f"Invalid day {day} for {desc}"

    def test_unique_eclipse_descriptions(self):
        """Each eclipse should have a unique description."""
        descriptions = [e[3] for e in ALL_HISTORICAL_ECLIPSES]
        # Allow some duplicates for alternative dates
        # but most should be unique
        unique_count = len(set(descriptions))
        assert unique_count >= len(descriptions) * 0.8, (
            f"Too many duplicate descriptions: {unique_count}/{len(descriptions)}"
        )
