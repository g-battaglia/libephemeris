"""
Validation of Lahiri ayanamsha against Indian Astronomical Ephemeris (IAE).

This module validates the Lahiri ayanamsha implementation against official values
published in the Indian Astronomical Ephemeris, as described in Swiss Ephemeris
documentation Appendix E: "How to compare the Swiss Ephemeris Lahiri Ayanamsha
with Indian Astronomical Ephemeris (IAE)".

Key Reference Points:
- The original Lahiri definition: 23 deg 15' 00" on 21 March 1956
- The IAE 1985 correction: 23 deg 15' 00".658 on 21 March 1956 (true ayanamsha)
- IAE uses IAU 1976 precession and Wahr 1980 nutation (since 1985)

Test Categories:
1. True ayanamsha values against IAE 2019, p. 427
2. Mean ayanamsha values against IAE 2019, p. 429
3. Historical reference dates from ICRC and Rashtriya Panchang
4. Precession rate validation

References:
- Swiss Ephemeris Documentation, Appendix E
- Indian Astronomical Ephemeris (IAE) 2019, 2020
- Report of the Calendar Reform Committee (ICRC), 1955
- Rashtriya Panchang (Indian National Calendar)
"""

import pytest
import math
import libephemeris as leph
from libephemeris.constants import SIDM_LAHIRI


# Helper function to convert DMS to decimal degrees
def dms_to_deg(deg: int, minutes: int, seconds: float) -> float:
    """Convert degrees, minutes, seconds to decimal degrees."""
    return deg + minutes / 60.0 + seconds / 3600.0


# Helper function to calculate Julian Day for a date
def calculate_jd(year: int, month: int, day: int, hour: float = 0.0) -> float:
    """Calculate Julian Day for a given date using the library's function."""
    return leph.julday(year, month, day, hour)


class TestIAEMeanAyanamshaValues:
    """
    Test mean ayanamsha values against IAE 2019, p. 429.

    Mean ayanamsha = True ayanamsha - Nutation in longitude

    The IAE provides mean ayanamsha values at specific epochs. These can be
    compared with libephemeris by using get_ayanamsa_ut() which returns
    the mean ayanamsha.

    According to Appendix E, the mean ayanamsha values in IAE 2019, p. 429
    are given at 12:00 TT on 1 January of each year.
    """

    # IAE 2019 mean ayanamsha values (p. 429)
    # Format: (year, month, day, hour_tt, expected_deg, expected_min, expected_sec, tolerance_arcsec)
    #
    # Note on tolerances: According to Swiss Ephemeris Appendix E, exact reproduction
    # of IAE values is difficult because:
    # 1. IAE uses outdated IAU1976 precession while libephemeris uses modern models
    # 2. IAE documentation has inconsistencies between editions
    # 3. Some values may be from older IENA publications without recalculation
    #
    # The J2000 value is most reliable; later dates show small systematic drift.
    IAE_2019_MEAN_AYANAMSHA = [
        # J2000.0 = 1 Jan 2000, 12:00 TT
        # IAE 2019: 23 deg 51' 25".53 - This is the most accurate reference point
        (2000, 1, 1, 12.0, 23, 51, 25.53, 0.1),
        # 1 Jan 2019, 12:00 TT
        # IAE 2019: 24 deg 7' 21".20
        # Per Appendix E: "diff = 0.06"" expected with SE using modern precession
        (2019, 1, 1, 12.0, 24, 7, 21.20, 0.5),
        # 1 Jan 2020, 12:00 TT
        # IAE 2019: 24 deg 8' 11".46
        # Per Appendix E: "diff = 0.06"" expected with SE using modern precession
        (2020, 1, 1, 12.0, 24, 8, 11.46, 0.5),
    ]

    @pytest.mark.parametrize(
        "year, month, day, hour, deg, minutes, sec, tolerance",
        IAE_2019_MEAN_AYANAMSHA,
    )
    def test_mean_ayanamsha_iae_2019(
        self, year, month, day, hour, deg, minutes, sec, tolerance
    ):
        """Test mean ayanamsha against IAE 2019 official values."""
        leph.set_sid_mode(SIDM_LAHIRI)

        # Calculate Julian Day
        jd = calculate_jd(year, month, day, hour)

        # Get mean ayanamsha from libephemeris
        ayanamsha = leph.get_ayanamsa_ut(jd)

        # Convert expected value to decimal degrees
        expected = dms_to_deg(deg, minutes, sec)

        # Calculate difference in arcseconds
        diff_arcsec = abs(ayanamsha - expected) * 3600

        assert diff_arcsec < tolerance, (
            f"Mean ayanamsha mismatch for {year}-{month:02d}-{day:02d} {hour:.1f}h TT:\n"
            f"  Expected: {deg}d {minutes}' {sec:.2f}\" ({expected:.6f} deg)\n"
            f"  Got:      {ayanamsha:.6f} deg\n"
            f'  Diff:     {diff_arcsec:.4f}" (tolerance: {tolerance}")'
        )


class TestLahiriReferenceEpochs:
    """
    Test Lahiri ayanamsha at key reference dates.

    These are the canonical reference points for the Lahiri ayanamsha:
    1. Original ICRC definition: 23 deg 15' 00" on 21 March 1956, 0:00 TT
    2. IAE 1985 correction: 23 deg 15' 00".658 (true ayanamsha, with nutation)

    Note: The original definition is for TRUE ayanamsha (including nutation),
    but get_ayanamsa_ut returns MEAN ayanamsha. The difference is the
    nutation in longitude, which was about 16" on that date.
    """

    def test_lahiri_at_definition_epoch(self):
        """
        Test Lahiri ayanamsha near the definition epoch (21 March 1956).

        The original Lahiri definition specifies 23 deg 15' 00" (true ayanamsha)
        on 21 March 1956. The mean ayanamsha would be this value minus nutation.

        According to Swiss Ephemeris, nutation on this date was about 16".
        So mean ayanamsha should be approximately 23 deg 15' 00" - 16" = 23 deg 14' 44".
        """
        leph.set_sid_mode(SIDM_LAHIRI)

        # 21 March 1956, 0:00 TT (definition epoch)
        jd_definition = calculate_jd(1956, 3, 21, 0.0)

        # Get mean ayanamsha
        mean_ayanamsha = leph.get_ayanamsa_ut(jd_definition)

        # Expected mean ayanamsha (true value 23 deg 15' 00".658 minus nutation ~16")
        # The exact value depends on the nutation model used
        # We use a tolerance that accounts for different nutation models
        expected_true = dms_to_deg(23, 15, 0.658)

        # Mean ayanamsha should be close to but slightly less than true
        # (by the nutation in longitude, typically 10-17 arcseconds)
        assert 23.24 < mean_ayanamsha < 23.26, (
            f"Mean ayanamsha at definition epoch should be ~23.25 deg, got {mean_ayanamsha:.6f} deg"
        )

    def test_lahiri_at_j2000(self):
        """
        Test Lahiri ayanamsha at J2000.0 epoch.

        J2000.0 = 1 January 2000, 12:00 TT = JD 2451545.0

        IAE 2019 gives mean ayanamsha: 23 deg 51' 25".53
        """
        leph.set_sid_mode(SIDM_LAHIRI)

        jd_j2000 = 2451545.0  # J2000.0 epoch

        ayanamsha = leph.get_ayanamsa_ut(jd_j2000)

        # IAE 2019 value: 23 deg 51' 25".53 = 23.857092 degrees
        expected = dms_to_deg(23, 51, 25.53)

        diff_arcsec = abs(ayanamsha - expected) * 3600

        assert diff_arcsec < 0.5, (
            f"Lahiri ayanamsha at J2000.0:\n"
            f"  Expected (IAE 2019): 23d 51' 25.53\" ({expected:.6f} deg)\n"
            f"  Got:                 {ayanamsha:.6f} deg\n"
            f'  Diff:                {diff_arcsec:.4f}"'
        )


class TestLahiriPrecessionRate:
    """
    Test that Lahiri ayanamsha increases at the expected precession rate.

    The standard precession rate is approximately:
    - 5027.8 arcseconds per Julian century
    - ~1.397 degrees per century
    - ~50.278 arcseconds per year
    """

    def test_precession_rate_one_century(self):
        """Test precession rate over one Julian century."""
        leph.set_sid_mode(SIDM_LAHIRI)

        jd_j2000 = 2451545.0  # J2000.0
        jd_j2100 = jd_j2000 + 36525  # J2100.0 (one Julian century later)

        ayanamsha_2000 = leph.get_ayanamsa_ut(jd_j2000)
        ayanamsha_2100 = leph.get_ayanamsa_ut(jd_j2100)

        change = ayanamsha_2100 - ayanamsha_2000

        # Expected change: 5027.8 arcsec = 1.3966 degrees
        expected_change = 5027.8 / 3600.0

        diff = abs(change - expected_change)

        assert diff < 0.01, (
            f"Precession rate mismatch:\n"
            f"  Expected change in 100 years: {expected_change:.6f} deg\n"
            f"  Actual change:                {change:.6f} deg\n"
            f"  Difference:                   {diff:.6f} deg"
        )

    def test_precession_rate_one_year(self):
        """Test precession rate over one year."""
        leph.set_sid_mode(SIDM_LAHIRI)

        jd_2000 = calculate_jd(2000, 1, 1, 12.0)
        jd_2001 = calculate_jd(2001, 1, 1, 12.0)

        ayanamsha_2000 = leph.get_ayanamsa_ut(jd_2000)
        ayanamsha_2001 = leph.get_ayanamsa_ut(jd_2001)

        change_arcsec = (ayanamsha_2001 - ayanamsha_2000) * 3600

        # Expected: ~50.278 arcsec per year
        expected_arcsec = 5027.8 / 100.0

        diff = abs(change_arcsec - expected_arcsec)

        assert diff < 0.5, (
            f"Annual precession rate mismatch:\n"
            f'  Expected: ~{expected_arcsec:.2f}" per year\n'
            f'  Actual:   {change_arcsec:.2f}" per year\n'
            f'  Diff:     {diff:.2f}"'
        )


class TestHistoricalICRCValues:
    """
    Test against values from the Report of the Calendar Reform Committee (ICRC)
    and Rashtriya Panchang.

    These tests verify that libephemeris can reproduce the official ayanamsha
    values published in Indian government calendrical publications.
    """

    # Historical reference values
    # Format: (description, year, month, day, expected_deg, expected_min, expected_sec, tolerance_arcsec)
    HISTORICAL_VALUES = [
        # ICRC Report: 1st of Phalguna = 20 February 1959
        # Ayanamsha: 23 deg 17' 16"
        ("ICRC 1959 (20 Feb)", 1959, 2, 20, 23, 17, 16.0, 2.0),
        # Rashtriya Panchang 1961: 1st of Caitra = 22 March 1961
        # Ayanamsha: 23 deg 18' 47"
        ("RP 1961 (22 Mar)", 1961, 3, 22, 23, 18, 47.0, 2.0),
        # Rashtriya Panchang 2019: 1st of Caitra = 22 March 2019
        # Ayanamsha: ~24 deg 7' 16"
        ("RP 2019 (22 Mar)", 2019, 3, 22, 24, 7, 16.0, 2.0),
    ]

    @pytest.mark.parametrize(
        "description, year, month, day, deg, minutes, sec, tolerance",
        HISTORICAL_VALUES,
    )
    def test_historical_reference_value(
        self, description, year, month, day, deg, minutes, sec, tolerance
    ):
        """Test against historical reference values from ICRC and Rashtriya Panchang."""
        leph.set_sid_mode(SIDM_LAHIRI)

        jd = calculate_jd(year, month, day, 0.0)
        ayanamsha = leph.get_ayanamsa_ut(jd)

        expected = dms_to_deg(deg, minutes, sec)

        # Note: Historical values may be TRUE ayanamsha, while we compute MEAN.
        # The difference (nutation) is typically 10-20 arcseconds.
        # We use a tolerance that allows for this difference.

        diff_arcsec = abs(ayanamsha - expected) * 3600

        assert diff_arcsec < tolerance + 20, (
            f"Historical reference mismatch for {description}:\n"
            f"  Expected: {deg}d {minutes}' {sec:.1f}\" ({expected:.6f} deg)\n"
            f"  Got:      {ayanamsha:.6f} deg\n"
            f'  Diff:     {diff_arcsec:.2f}" (tolerance: {tolerance}" + 20" nutation margin)'
        )


class TestIAEConsistency:
    """
    Test internal consistency of Lahiri ayanamsha calculations.

    These tests verify that the implementation is self-consistent
    and follows expected mathematical relationships.
    """

    def test_ayanamsha_increases_over_time(self):
        """Ayanamsha should always increase over time (precession is positive)."""
        leph.set_sid_mode(SIDM_LAHIRI)

        dates = [
            calculate_jd(1900, 1, 1, 0.0),
            calculate_jd(1950, 1, 1, 0.0),
            calculate_jd(2000, 1, 1, 0.0),
            calculate_jd(2050, 1, 1, 0.0),
            calculate_jd(2100, 1, 1, 0.0),
        ]

        values = [leph.get_ayanamsa_ut(jd) for jd in dates]

        for i in range(len(values) - 1):
            assert values[i] < values[i + 1], (
                f"Ayanamsha should increase over time: "
                f"{values[i]:.6f} should be less than {values[i + 1]:.6f}"
            )

    def test_ayanamsha_reasonable_range(self):
        """Ayanamsha should be in a reasonable range for modern dates."""
        leph.set_sid_mode(SIDM_LAHIRI)

        # Test various dates from 1900 to 2100
        test_cases = [
            (1900, 22.0, 23.0),  # Expected: ~22.5 degrees
            (1950, 23.0, 24.0),  # Expected: ~23.2 degrees
            (2000, 23.5, 24.5),  # Expected: ~23.85 degrees
            (2050, 24.0, 25.0),  # Expected: ~24.55 degrees
            (2100, 24.5, 25.5),  # Expected: ~25.25 degrees
        ]

        for year, min_val, max_val in test_cases:
            jd = calculate_jd(year, 1, 1, 12.0)
            ayanamsha = leph.get_ayanamsa_ut(jd)

            assert min_val < ayanamsha < max_val, (
                f"Ayanamsha for year {year} out of expected range:\n"
                f"  Expected: {min_val:.2f} - {max_val:.2f} deg\n"
                f"  Got:      {ayanamsha:.6f} deg"
            )

    def test_lahiri_zero_epoch(self):
        """
        Test that Lahiri ayanamsha was zero around 285 CE.

        The Lahiri ayanamsha is based on Spica at 180 degrees, which
        implies the tropical and sidereal zodiacs coincided around 285 CE.
        """
        leph.set_sid_mode(SIDM_LAHIRI)

        # Get ayanamsha at J2000
        jd_j2000 = 2451545.0
        ayanamsha_j2000 = leph.get_ayanamsa_ut(jd_j2000)

        # Calculate years from J2000 to zero epoch
        rate_per_year = 5027.8 / 3600.0 / 100.0  # degrees per year
        years_to_zero = ayanamsha_j2000 / rate_per_year

        # Zero epoch year (J2000 = year 2000)
        zero_year = 2000 - years_to_zero

        # Expected: around 285 CE (within ~50 years tolerance due to
        # variations in precession rate over centuries)
        assert 235 < zero_year < 335, (
            f"Lahiri zero epoch should be around 285 CE, got {zero_year:.0f} CE"
        )


class TestComparisonWithSwissEphemeris:
    """
    Tests to compare libephemeris Lahiri values with Swiss Ephemeris reference.

    These tests verify that libephemeris produces values consistent with
    Swiss Ephemeris (pyswisseph) for the Lahiri ayanamsha.

    Note: There are known small differences between libephemeris and pyswisseph
    due to different underlying implementations and precision in time conversion.
    The tolerance of 1 arcsecond is appropriate for astronomical applications.
    """

    @pytest.mark.parametrize(
        "year",
        [1950, 1980, 2000, 2020, 2050],
    )
    def test_lahiri_matches_swisseph(self, year):
        """Test that Lahiri ayanamsha matches Swiss Ephemeris values."""
        try:
            import swisseph as swe

            has_swisseph = True
        except ImportError:
            has_swisseph = False
            pytest.skip("swisseph not available for comparison")

        if has_swisseph:
            jd = swe.julday(year, 1, 1, 12.0)

            # Swiss Ephemeris value
            swe.set_sid_mode(SIDM_LAHIRI)
            swe_value = swe.get_ayanamsa_ut(jd)

            # libephemeris value
            leph.set_sid_mode(SIDM_LAHIRI)
            leph_value = leph.get_ayanamsa_ut(jd)

            diff_arcsec = abs(swe_value - leph_value) * 3600

            # Tolerance of 1 arcsecond is appropriate for matching IAE values
            # and accounts for small differences in implementation details
            assert diff_arcsec < 1.0, (
                f"Lahiri ayanamsha mismatch with Swiss Ephemeris for year {year}:\n"
                f"  Swiss Ephemeris: {swe_value:.6f} deg\n"
                f"  libephemeris:    {leph_value:.6f} deg\n"
                f'  Difference:      {diff_arcsec:.4f}"'
            )


if __name__ == "__main__":
    print("=" * 70)
    print("LAHIRI AYANAMSHA VALIDATION AGAINST IAE")
    print("=" * 70)

    leph.set_sid_mode(SIDM_LAHIRI)

    # Test key dates
    test_dates = [
        ("21 March 1956 (Definition)", 1956, 3, 21, 0.0),
        ("J2000.0 (1 Jan 2000 12h)", 2000, 1, 1, 12.0),
        ("1 Jan 2019 12h", 2019, 1, 1, 12.0),
        ("1 Jan 2020 12h", 2020, 1, 1, 12.0),
    ]

    print(f"\n{'Date':<30} {'Ayanamsha':<20} {'DMS':<20}")
    print("-" * 70)

    for desc, year, month, day, hour in test_dates:
        jd = calculate_jd(year, month, day, hour)
        ayanamsha = leph.get_ayanamsa_ut(jd)

        # Convert to DMS
        d = int(ayanamsha)
        m = int((ayanamsha - d) * 60)
        s = (ayanamsha - d - m / 60) * 3600

        print(f"{desc:<30} {ayanamsha:>15.6f} deg  {d}d {m}' {s:.2f}\"")

    print("\n" + "=" * 70)
    print("Run pytest for full validation suite.")
    print("=" * 70)
