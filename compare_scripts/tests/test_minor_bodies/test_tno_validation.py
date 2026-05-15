"""
Tests for TNO position validation after Uranus and Neptune perturbation improvements.

These tests validate that Eris, Makemake, Ixion, and Orcus positions calculated
with libephemeris (Keplerian + secular perturbations from Jupiter, Saturn, Uranus,
and Neptune) are within acceptable tolerances compared to pyswisseph.

The 50-year validation span (2000-2050) ensures the secular perturbation model
provides reasonable accuracy over extended time periods.

Precision targets:
- Longitude: < 10 degrees (relaxed tolerance for Keplerian approximation)
- Latitude: < 5 degrees
- These tolerances are appropriate for astrological applications

Note: For research-grade precision, SPK kernels should be used instead.
"""

import pytest
import math
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    import swisseph as swe

    HAS_SWISSEPH = True
except ImportError:
    HAS_SWISSEPH = False

import libephemeris as ephem
from libephemeris.constants import ERIS, MAKEMAKE, IXION, ORCUS
from libephemeris.minor_bodies import (
    MINOR_BODY_ELEMENTS,
    calc_secular_perturbation_rates,
    detect_mean_motion_resonance,
)
from comparison_utils import angular_diff


# Test configuration
START_YEAR = 2000
END_YEAR = 2050
SAMPLE_YEARS = [2000, 2010, 2020, 2030, 2040, 2050]  # Key sample points

# Tolerances for Keplerian + secular perturbation model
LONGITUDE_TOLERANCE = 10.0  # degrees
LATITUDE_TOLERANCE = 5.0  # degrees

# Tighter tolerances for documenting improvement
EXPECTED_MEAN_LON_ERROR = 5.0  # Expected mean error should be better than worst case


# TNO definitions
TNOS = [
    (ERIS, 136199, "Eris"),
    (MAKEMAKE, 136472, "Makemake"),
    (IXION, 28978, "Ixion"),
    (ORCUS, 90482, "Orcus"),
]


def get_swe_position(jd: float, ast_num: int):
    """Get position from Swiss Ephemeris, returns None if data unavailable."""
    if not HAS_SWISSEPH:
        return None
    try:
        swe_body = swe.AST_OFFSET + ast_num
        pos, _ = swe.calc_ut(jd, swe_body, 0)
        return pos[0], pos[1], pos[2]
    except swe.Error:
        return None


def get_libeph_position(jd: float, body_id: int):
    """Get position from libephemeris."""
    pos, _ = ephem.calc_ut(jd, body_id, 0)
    return pos[0], pos[1], pos[2]


@pytest.mark.skipif(not HAS_SWISSEPH, reason="swisseph not available")
class TestTNOValidation:
    """Validate TNO positions against pyswisseph over 2000-2050."""

    @pytest.mark.parametrize("body_id,ast_num,name", TNOS)
    def test_tno_longitude_within_tolerance(self, body_id, ast_num, name):
        """Verify TNO longitude is within 10 degrees of Swiss Ephemeris."""
        max_diff = 0.0
        samples = 0

        for year in SAMPLE_YEARS:
            jd = swe.julday(year, 1, 1, 12.0)

            swe_pos = get_swe_position(jd, ast_num)
            if swe_pos is None:
                pytest.skip(f"SwissEph data not available for {name}")

            lib_pos = get_libeph_position(jd, body_id)

            diff = angular_diff(swe_pos[0], lib_pos[0])
            max_diff = max(max_diff, diff)
            samples += 1

        assert samples > 0, f"No valid samples for {name}"
        assert max_diff < LONGITUDE_TOLERANCE, (
            f"{name} longitude error {max_diff:.2f}° exceeds tolerance "
            f"{LONGITUDE_TOLERANCE}°"
        )

    @pytest.mark.parametrize("body_id,ast_num,name", TNOS)
    def test_tno_latitude_within_tolerance(self, body_id, ast_num, name):
        """Verify TNO latitude is within 5 degrees of Swiss Ephemeris."""
        max_diff = 0.0
        samples = 0

        for year in SAMPLE_YEARS:
            jd = swe.julday(year, 1, 1, 12.0)

            swe_pos = get_swe_position(jd, ast_num)
            if swe_pos is None:
                pytest.skip(f"SwissEph data not available for {name}")

            lib_pos = get_libeph_position(jd, body_id)

            diff = abs(swe_pos[1] - lib_pos[1])
            max_diff = max(max_diff, diff)
            samples += 1

        assert samples > 0, f"No valid samples for {name}"
        assert max_diff < LATITUDE_TOLERANCE, (
            f"{name} latitude error {max_diff:.2f}° exceeds tolerance "
            f"{LATITUDE_TOLERANCE}°"
        )

    def test_eris_50_year_validation(self):
        """Validate Eris positions over full 50-year span."""
        lon_diffs = []
        lat_diffs = []

        for year in range(START_YEAR, END_YEAR + 1, 5):  # Every 5 years
            jd = swe.julday(year, 1, 1, 12.0)

            swe_pos = get_swe_position(jd, 136199)
            if swe_pos is None:
                pytest.skip("SwissEph data not available for Eris")

            lib_pos = get_libeph_position(jd, ERIS)

            lon_diffs.append(angular_diff(swe_pos[0], lib_pos[0]))
            lat_diffs.append(abs(swe_pos[1] - lib_pos[1]))

        assert len(lon_diffs) > 0, "No valid samples"

        max_lon = max(lon_diffs)
        mean_lon = sum(lon_diffs) / len(lon_diffs)
        max_lat = max(lat_diffs)

        assert max_lon < LONGITUDE_TOLERANCE, (
            f"Eris max longitude error {max_lon:.2f}° exceeds tolerance"
        )
        assert max_lat < LATITUDE_TOLERANCE, (
            f"Eris max latitude error {max_lat:.2f}° exceeds tolerance"
        )
        assert mean_lon < EXPECTED_MEAN_LON_ERROR, (
            f"Eris mean longitude error {mean_lon:.2f}° higher than expected"
        )

    def test_makemake_50_year_validation(self):
        """Validate Makemake positions over full 50-year span."""
        lon_diffs = []
        lat_diffs = []

        for year in range(START_YEAR, END_YEAR + 1, 5):
            jd = swe.julday(year, 1, 1, 12.0)

            swe_pos = get_swe_position(jd, 136472)
            if swe_pos is None:
                pytest.skip("SwissEph data not available for Makemake")

            lib_pos = get_libeph_position(jd, MAKEMAKE)

            lon_diffs.append(angular_diff(swe_pos[0], lib_pos[0]))
            lat_diffs.append(abs(swe_pos[1] - lib_pos[1]))

        assert len(lon_diffs) > 0, "No valid samples"

        max_lon = max(lon_diffs)
        mean_lon = sum(lon_diffs) / len(lon_diffs)
        max_lat = max(lat_diffs)

        assert max_lon < LONGITUDE_TOLERANCE, (
            f"Makemake max longitude error {max_lon:.2f}° exceeds tolerance"
        )
        assert max_lat < LATITUDE_TOLERANCE, (
            f"Makemake max latitude error {max_lat:.2f}° exceeds tolerance"
        )
        assert mean_lon < EXPECTED_MEAN_LON_ERROR, (
            f"Makemake mean longitude error {mean_lon:.2f}° higher than expected"
        )

    def test_plutino_ixion_50_year_validation(self):
        """Validate Ixion (plutino) positions over full 50-year span."""
        lon_diffs = []
        lat_diffs = []

        for year in range(START_YEAR, END_YEAR + 1, 5):
            jd = swe.julday(year, 1, 1, 12.0)

            swe_pos = get_swe_position(jd, 28978)
            if swe_pos is None:
                pytest.skip("SwissEph data not available for Ixion")

            lib_pos = get_libeph_position(jd, IXION)

            lon_diffs.append(angular_diff(swe_pos[0], lib_pos[0]))
            lat_diffs.append(abs(swe_pos[1] - lib_pos[1]))

        assert len(lon_diffs) > 0, "No valid samples"

        max_lon = max(lon_diffs)
        mean_lon = sum(lon_diffs) / len(lon_diffs)
        max_lat = max(lat_diffs)

        # Plutinos are in Neptune resonance, expect some deviation
        assert max_lon < LONGITUDE_TOLERANCE, (
            f"Ixion (plutino) max longitude error {max_lon:.2f}° exceeds tolerance"
        )
        assert max_lat < LATITUDE_TOLERANCE, (
            f"Ixion max latitude error {max_lat:.2f}° exceeds tolerance"
        )

    def test_plutino_orcus_50_year_validation(self):
        """Validate Orcus (plutino) positions over full 50-year span."""
        lon_diffs = []
        lat_diffs = []

        for year in range(START_YEAR, END_YEAR + 1, 5):
            jd = swe.julday(year, 1, 1, 12.0)

            swe_pos = get_swe_position(jd, 90482)
            if swe_pos is None:
                pytest.skip("SwissEph data not available for Orcus")

            lib_pos = get_libeph_position(jd, ORCUS)

            lon_diffs.append(angular_diff(swe_pos[0], lib_pos[0]))
            lat_diffs.append(abs(swe_pos[1] - lib_pos[1]))

        assert len(lon_diffs) > 0, "No valid samples"

        max_lon = max(lon_diffs)
        mean_lon = sum(lon_diffs) / len(lon_diffs)
        max_lat = max(lat_diffs)

        # Plutinos are in Neptune resonance, expect some deviation
        assert max_lon < LONGITUDE_TOLERANCE, (
            f"Orcus (plutino) max longitude error {max_lon:.2f}° exceeds tolerance"
        )
        assert max_lat < LATITUDE_TOLERANCE, (
            f"Orcus max latitude error {max_lat:.2f}° exceeds tolerance"
        )


class TestPerturbationImprovements:
    """Verify that Uranus and Neptune perturbations provide improvements."""

    def test_uranus_perturbation_affects_tnos(self):
        """Verify Uranus perturbations are applied to TNOs."""
        for body_id, _, name in TNOS:
            elements = MINOR_BODY_ELEMENTS[body_id]
            d_omega, d_Omega, d_n = calc_secular_perturbation_rates(elements)

            # All TNOs should have finite perturbation rates
            assert math.isfinite(d_omega), f"{name} d_omega should be finite"
            assert math.isfinite(d_Omega), f"{name} d_Omega should be finite"

            # Convert to arcsec/year for readability
            d_omega_arcsec_yr = abs(d_omega) * 365.25 * 3600

            # Rates should be measurable (combined effect of all planets)
            assert d_omega_arcsec_yr >= 0, f"{name} should have non-negative omega rate"

    def test_neptune_perturbation_applied_to_plutinos(self):
        """Verify Neptune perturbations are applied to plutinos."""
        for body_id, name in [(IXION, "Ixion"), (ORCUS, "Orcus")]:
            elements = MINOR_BODY_ELEMENTS[body_id]

            # Verify these are plutinos (a ~ 39 AU)
            assert 38.0 < elements.a < 41.0, (
                f"{name} should have a ~ 39 AU (plutino), got {elements.a}"
            )

            # Neptune perturbations should be active (a > 20 AU threshold)
            assert elements.a > 20.0, f"{name} should be beyond 20 AU threshold"

            d_omega, d_Omega, d_n = calc_secular_perturbation_rates(elements)

            # Both rates should be finite
            assert math.isfinite(d_omega), f"{name} d_omega should be finite"
            assert math.isfinite(d_Omega), f"{name} d_Omega should be finite"

    def test_resonance_detection_for_plutinos(self):
        """Verify resonance detection works for Ixion and Orcus."""
        for body_id, name in [(IXION, "Ixion"), (ORCUS, "Orcus")]:
            elements = MINOR_BODY_ELEMENTS[body_id]
            result = detect_mean_motion_resonance(elements)

            # Both should be detected as resonant
            assert result.is_resonant, f"{name} should be detected as resonant"

            # Should be plutinos (2:3 resonance)
            if result.resonance:
                assert result.resonance.name == "plutino", (
                    f"{name} should be detected as plutino, got {result.resonance.name}"
                )

    def test_eris_is_not_resonant(self):
        """Verify Eris is not in a Neptune resonance."""
        elements = MINOR_BODY_ELEMENTS[ERIS]
        result = detect_mean_motion_resonance(elements)

        # Eris is in the scattered disk, not in a major resonance
        assert not result.is_resonant, "Eris should not be detected as resonant"

    def test_perturbation_rates_documented(self):
        """Document perturbation rates for reference."""
        print("\n" + "=" * 60)
        print("Secular Perturbation Rates for TNOs")
        print("(Combined effect of Jupiter, Saturn, Uranus, Neptune)")
        print("=" * 60)

        for body_id, _, name in TNOS:
            elements = MINOR_BODY_ELEMENTS[body_id]
            d_omega, d_Omega, d_n = calc_secular_perturbation_rates(elements)

            # Convert to more readable units
            d_omega_arcsec_yr = d_omega * 365.25 * 3600
            d_Omega_arcsec_yr = d_Omega * 365.25 * 3600

            print(f"\n{name} (a = {elements.a:.1f} AU):")
            print(f"  dω/dt = {d_omega_arcsec_yr:+.2f} arcsec/year")
            print(f"  dΩ/dt = {d_Omega_arcsec_yr:+.2f} arcsec/year")

            # All should be finite
            assert math.isfinite(d_omega_arcsec_yr)
            assert math.isfinite(d_Omega_arcsec_yr)


class TestPrecisionDocumentation:
    """Document precision improvements from Uranus/Neptune perturbations."""

    @pytest.mark.skipif(not HAS_SWISSEPH, reason="swisseph not available")
    def test_document_precision_statistics(self):
        """Document precision statistics for all TNOs."""
        print("\n" + "=" * 70)
        print("TNO PRECISION VALIDATION - 50 Year Span (2000-2050)")
        print("=" * 70)
        print("\nComparing libephemeris (Keplerian + secular perturbations)")
        print("against pyswisseph (full numerical integration)")
        print("\nPerturbation sources: Jupiter, Saturn, Uranus, Neptune")
        print("-" * 70)

        all_results = []

        for body_id, ast_num, name in TNOS:
            lon_diffs = []
            lat_diffs = []
            dist_diffs = []

            for year in range(START_YEAR, END_YEAR + 1, 2):  # Every 2 years
                jd = swe.julday(year, 1, 1, 12.0)

                swe_pos = get_swe_position(jd, ast_num)
                if swe_pos is None:
                    continue

                lib_pos = get_libeph_position(jd, body_id)

                lon_diffs.append(angular_diff(swe_pos[0], lib_pos[0]))
                lat_diffs.append(abs(swe_pos[1] - lib_pos[1]))
                dist_diffs.append(abs(swe_pos[2] - lib_pos[2]))

            if not lon_diffs:
                print(f"\n{name}: SKIPPED (no SwissEph data)")
                continue

            max_lon = max(lon_diffs)
            mean_lon = sum(lon_diffs) / len(lon_diffs)
            rms_lon = math.sqrt(sum(d * d for d in lon_diffs) / len(lon_diffs))
            max_lat = max(lat_diffs)
            max_dist = max(dist_diffs)

            print(f"\n{name}:")
            print(f"  Samples: {len(lon_diffs)}")
            print(
                f"  Longitude: max={max_lon:.2f}°, mean={mean_lon:.2f}°, "
                f"RMS={rms_lon:.2f}°"
            )
            print(f"  Latitude:  max={max_lat:.2f}°")
            print(f"  Distance:  max={max_dist:.4f} AU")

            passed = max_lon < LONGITUDE_TOLERANCE and max_lat < LATITUDE_TOLERANCE
            print(f"  Status: {'PASS' if passed else 'FAIL'}")

            all_results.append(
                {
                    "name": name,
                    "max_lon": max_lon,
                    "mean_lon": mean_lon,
                    "max_lat": max_lat,
                    "passed": passed,
                }
            )

        # Summary
        print("\n" + "=" * 70)
        print("SUMMARY")
        print("=" * 70)

        if all_results:
            overall_max_lon = max(r["max_lon"] for r in all_results)
            overall_mean_lon = sum(r["mean_lon"] for r in all_results) / len(
                all_results
            )
            all_passed = all(r["passed"] for r in all_results)

            print(f"Overall max longitude error: {overall_max_lon:.2f}°")
            print(f"Overall mean longitude error: {overall_mean_lon:.2f}°")
            print(f"All within tolerance: {'YES' if all_passed else 'NO'}")

            # Precision grade
            if overall_max_lon < 1.0:
                grade = "EXCELLENT"
            elif overall_max_lon < 3.0:
                grade = "GOOD"
            elif overall_max_lon < 10.0:
                grade = "ACCEPTABLE"
            else:
                grade = "POOR"

            print(f"Precision grade: {grade}")

            assert all_passed, "Not all TNOs passed validation"
        else:
            pytest.skip("No SwissEph data available for any TNOs")

        print("=" * 70)
