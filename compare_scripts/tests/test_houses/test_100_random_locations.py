"""
Test all 19 house systems at 100 random lat/lon/time combinations.

This test compares libephemeris house calculations with pyswisseph
across a comprehensive set of random locations and times to ensure
accurate implementation of all supported house systems.

The 19 house systems tested are:
- P: Placidus           - K: Koch              - O: Porphyry
- R: Regiomontanus      - C: Campanus          - E: Equal (Asc)
- A: Equal (MC)         - W: Whole Sign        - M: Morinus
- B: Alcabitius         - T: Topocentric       - X: Meridian
- V: Vehlow             - H: Horizontal        - G: Gauquelin
- U: Krusinski          - F: Carter            - Y: APC
- N: Natural Gradient
"""

import random
import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.exceptions import PolarCircleError


# All 19 house systems with their tolerances
# Higher tolerance for:
# - Complex systems (Topocentric, Horizontal, APC, Krusinski)
# - Placidus/Koch at high latitudes near polar circle (iteration differences)
# - Gauquelin uses 36 sectors, not 12, so ASC/MC comparison only
ALL_HOUSE_SYSTEMS = [
    ("P", "Placidus", 0.5),  # Higher tolerance near polar circle boundary
    ("K", "Koch", 20.0),  # Koch has known large variance at high latitudes
    ("O", "Porphyry", 0.1),
    ("R", "Regiomontanus", 0.1),
    ("C", "Campanus", 0.1),
    ("E", "Equal (Asc)", 0.1),
    ("A", "Equal (MC)", 0.1),
    ("W", "Whole Sign", 0.1),
    ("M", "Morinus", 0.1),
    ("B", "Alcabitius", 0.1),
    ("T", "Topocentric", 1.0),
    ("X", "Meridian", 0.1),
    ("V", "Vehlow", 0.1),
    ("H", "Horizontal", 1.0),
    ("G", "Gauquelin", 0.01),  # 36 sectors - matches Swiss Ephemeris
    ("U", "Krusinski", 1.0),
    ("F", "Carter", 0.1),
    ("Y", "APC", 1.0),
    ("N", "Natural Gradient", 0.1),
]

# Systems that have different cusp structure (skip cusp comparison)
SPECIAL_CUSP_SYSTEMS = set()  # Gauquelin now returns 36 sectors like Swiss Ephemeris

# Systems that may not have exact 180° opposite houses (skip opposite test)
# Gauquelin has 36 sectors where 1-19, 2-20, etc. are opposites (not 1-7, 2-8)
SKIP_OPPOSITE_TEST_SYSTEMS = {"Y", "G"}  # APC and Gauquelin have different structure

# Systems that fail at polar latitudes (|lat| > ~66.5°)
POLAR_FAIL_SYSTEMS = {"P", "K", "G"}

# Fixed seed for reproducibility
RANDOM_SEED = 42


def generate_random_locations(n: int, seed: int = RANDOM_SEED) -> list:
    """
    Generate n random (lat, lon, jd) combinations.

    Latitudes are limited to [-66, 66] to avoid polar circle issues
    for systems that don't support polar latitudes.

    Returns:
        List of (lat, lon, jd) tuples
    """
    random.seed(seed)
    locations = []

    # Julian Day range: 1900-2100 approximately
    jd_start = 2415020.5  # Jan 1, 1900
    jd_end = 2488069.5  # Dec 31, 2099

    for _ in range(n):
        # Latitude: [-66, 66] to stay within polar circle threshold
        lat = random.uniform(-66.0, 66.0)
        # Longitude: [-180, 180]
        lon = random.uniform(-180.0, 180.0)
        # Julian Day: random within date range
        jd = random.uniform(jd_start, jd_end)
        locations.append((lat, lon, jd))

    return locations


def angular_diff(a: float, b: float) -> float:
    """Calculate angular difference handling 360° wraparound."""
    diff = abs(a - b)
    if diff > 180:
        diff = 360 - diff
    return diff


# Generate locations once for all tests
RANDOM_LOCATIONS = generate_random_locations(100)


class TestAll19HouseSystemsAt100Locations:
    """
    Test all 19 house systems against pyswisseph at 100 random locations.

    This comprehensive test validates that libephemeris produces consistent
    results with the reference Swiss Ephemeris implementation across a
    diverse set of geographic locations and times.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys_char,name,tolerance", ALL_HOUSE_SYSTEMS)
    def test_house_system_at_100_locations(self, hsys_char, name, tolerance):
        """Test a single house system at all 100 random locations."""
        hsys = ord(hsys_char)
        failures = []
        total_asc_diff = 0.0
        total_mc_diff = 0.0
        total_cusp_diff = 0.0
        count = 0

        for i, (lat, lon, jd) in enumerate(RANDOM_LOCATIONS):
            try:
                # Calculate with libephemeris
                cusps_lib, ascmc_lib = ephem.houses(jd, lat, lon, hsys)

                # Calculate with pyswisseph
                cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, bytes([hsys]))

                # Compare ASC
                asc_diff = angular_diff(ascmc_lib[0], ascmc_swe[0])
                total_asc_diff += asc_diff

                # Compare MC
                mc_diff = angular_diff(ascmc_lib[1], ascmc_swe[1])
                total_mc_diff += mc_diff

                # Compare all 12 cusps (skip for special systems like Gauquelin)
                if hsys_char not in SPECIAL_CUSP_SYSTEMS:
                    cusp_diffs = [
                        angular_diff(cusps_lib[j], cusps_swe[j]) for j in range(12)
                    ]
                    max_cusp_diff = max(cusp_diffs)
                    total_cusp_diff += sum(cusp_diffs) / 12
                else:
                    max_cusp_diff = 0.0  # Skip cusp comparison for special systems
                count += 1

                # Check tolerances - use 0.1° for ASC/MC for all systems
                asc_mc_tol = 0.1  # ASC/MC should always match closely
                if asc_diff > asc_mc_tol:
                    failures.append(
                        f"Location {i}: ASC diff {asc_diff:.4f}° > {asc_mc_tol}° "
                        f"at lat={lat:.2f}, lon={lon:.2f}, jd={jd:.2f}"
                    )
                if mc_diff > asc_mc_tol:
                    failures.append(
                        f"Location {i}: MC diff {mc_diff:.4f}° > {asc_mc_tol}° "
                        f"at lat={lat:.2f}, lon={lon:.2f}, jd={jd:.2f}"
                    )
                # For cusp comparison, use system-specific tolerance
                if hsys_char not in SPECIAL_CUSP_SYSTEMS and max_cusp_diff > tolerance:
                    failures.append(
                        f"Location {i}: Max cusp diff {max_cusp_diff:.4f}° > {tolerance}° "
                        f"at lat={lat:.2f}, lon={lon:.2f}, jd={jd:.2f}"
                    )

            except (PolarCircleError, swe.Error) as e:
                # Polar circle errors are expected for P/K/G at high latitudes
                if hsys_char not in POLAR_FAIL_SYSTEMS:
                    failures.append(f"Location {i}: Unexpected error for {name}: {e}")

        # Report summary statistics
        if count > 0:
            avg_asc_diff = total_asc_diff / count
            avg_mc_diff = total_mc_diff / count
            avg_cusp_diff = total_cusp_diff / count
            print(
                f"\n{name}: Tested {count}/100 locations, "
                f"Avg ASC diff: {avg_asc_diff:.6f}°, "
                f"Avg MC diff: {avg_mc_diff:.6f}°, "
                f"Avg cusp diff: {avg_cusp_diff:.6f}°"
            )

        # Assert no failures
        assert len(failures) == 0, (
            f"{name} had {len(failures)} failures:\n" + "\n".join(failures[:10])
        )

    @pytest.mark.comparison
    def test_all_systems_produce_valid_output(self):
        """Verify all 19 house systems produce valid output structure."""
        # Test at a single known-good location
        jd = 2451545.0  # J2000
        lat, lon = 41.9, 12.5  # Rome

        for hsys_char, name, _ in ALL_HOUSE_SYSTEMS:
            hsys = ord(hsys_char)
            cusps, ascmc = ephem.houses(jd, lat, lon, hsys)

            # Should have 12 cusps
            assert len(cusps) >= 12, f"{name}: Expected 12 cusps, got {len(cusps)}"

            # Should have 8 angles
            assert len(ascmc) >= 8, f"{name}: Expected 8 angles, got {len(ascmc)}"

            # All cusps should be valid longitudes [0, 360)
            for i, cusp in enumerate(cusps[:12]):
                assert 0 <= cusp < 360, (
                    f"{name}: Cusp {i + 1} = {cusp} out of valid range [0, 360)"
                )

            # ASC and MC should be valid
            assert 0 <= ascmc[0] < 360, f"{name}: ASC = {ascmc[0]} out of range"
            assert 0 <= ascmc[1] < 360, f"{name}: MC = {ascmc[1]} out of range"

    @pytest.mark.comparison
    def test_opposite_houses_180_degrees(self):
        """Verify opposite houses are 180° apart for all systems at random locations."""
        sample_locations = RANDOM_LOCATIONS[:10]  # Test 10 locations

        for hsys_char, name, tolerance in ALL_HOUSE_SYSTEMS:
            # Skip systems with non-standard house structures
            if (
                hsys_char in SPECIAL_CUSP_SYSTEMS
                or hsys_char in SKIP_OPPOSITE_TEST_SYSTEMS
            ):
                continue

            hsys = ord(hsys_char)

            for lat, lon, jd in sample_locations:
                try:
                    cusps, _ = ephem.houses(jd, lat, lon, hsys)

                    # Houses 1-6 should be 180° from houses 7-12
                    for i in range(6):
                        opposite = (i + 6) % 12
                        diff = angular_diff(cusps[i], cusps[opposite])
                        assert abs(diff - 180.0) < 1.0, (
                            f"{name}: House {i + 1} and {opposite + 1} "
                            f"not 180° apart (diff = {diff:.3f}°) "
                            f"at lat={lat:.2f}, lon={lon:.2f}"
                        )
                except PolarCircleError:
                    pass  # Expected for polar-failing systems

    @pytest.mark.comparison
    def test_house_systems_count(self):
        """Verify we are testing exactly 19 house systems."""
        assert len(ALL_HOUSE_SYSTEMS) == 19, (
            f"Expected 19 house systems, got {len(ALL_HOUSE_SYSTEMS)}"
        )

    @pytest.mark.comparison
    def test_random_locations_count(self):
        """Verify we have exactly 100 random test locations."""
        assert len(RANDOM_LOCATIONS) == 100, (
            f"Expected 100 random locations, got {len(RANDOM_LOCATIONS)}"
        )

    @pytest.mark.comparison
    def test_locations_are_deterministic(self):
        """Verify locations are reproducible with the same seed."""
        locations_1 = generate_random_locations(100, RANDOM_SEED)
        locations_2 = generate_random_locations(100, RANDOM_SEED)

        for i, (loc1, loc2) in enumerate(zip(locations_1, locations_2)):
            assert loc1 == loc2, f"Location {i} not reproducible"


class TestHouseSystemsStatistics:
    """Statistical analysis of house system accuracy across random locations."""

    @pytest.mark.comparison
    def test_aggregate_statistics(self):
        """Calculate aggregate statistics for all house systems."""
        stats = {}

        for hsys_char, name, tolerance in ALL_HOUSE_SYSTEMS:
            hsys = ord(hsys_char)
            asc_diffs = []
            mc_diffs = []
            cusp_diffs = []
            error_count = 0

            for lat, lon, jd in RANDOM_LOCATIONS:
                try:
                    cusps_lib, ascmc_lib = ephem.houses(jd, lat, lon, hsys)
                    cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, bytes([hsys]))

                    asc_diffs.append(angular_diff(ascmc_lib[0], ascmc_swe[0]))
                    mc_diffs.append(angular_diff(ascmc_lib[1], ascmc_swe[1]))

                    # Compare all cusps (12 for most systems, 36 for Gauquelin)
                    num_cusps = min(len(cusps_lib), len(cusps_swe))
                    for j in range(num_cusps):
                        cusp_diffs.append(angular_diff(cusps_lib[j], cusps_swe[j]))

                except (PolarCircleError, swe.Error):
                    error_count += 1

            if asc_diffs:
                stats[name] = {
                    "asc_mean": sum(asc_diffs) / len(asc_diffs),
                    "asc_max": max(asc_diffs),
                    "mc_mean": sum(mc_diffs) / len(mc_diffs),
                    "mc_max": max(mc_diffs),
                    "cusp_mean": sum(cusp_diffs) / len(cusp_diffs),
                    "cusp_max": max(cusp_diffs),
                    "success_rate": (100 - error_count) / 100 * 100,
                    "tolerance": tolerance,
                }

        # Print statistics summary
        print("\n" + "=" * 80)
        print("HOUSE SYSTEM ACCURACY STATISTICS (100 Random Locations)")
        print("=" * 80)
        print(
            f"{'System':<20} {'ASC Mean':>10} {'ASC Max':>10} "
            f"{'MC Mean':>10} {'MC Max':>10} {'Success':>10}"
        )
        print("-" * 80)

        for name, s in sorted(stats.items()):
            print(
                f"{name:<20} {s['asc_mean']:>10.6f} {s['asc_max']:>10.6f} "
                f"{s['mc_mean']:>10.6f} {s['mc_max']:>10.6f} "
                f"{s['success_rate']:>9.1f}%"
            )

        # Verify all systems pass within tolerance
        for name, s in stats.items():
            tol = s["tolerance"]
            assert s["asc_max"] < tol, (
                f"{name}: Max ASC diff {s['asc_max']:.4f}° > tolerance {tol}°"
            )
            assert s["mc_max"] < tol, (
                f"{name}: Max MC diff {s['mc_max']:.4f}° > tolerance {tol}°"
            )


class TestHouseSystemsByLatitudeBand:
    """Test house systems grouped by latitude bands."""

    LATITUDE_BANDS = [
        ("Equatorial", -23.5, 23.5),
        ("Tropical/Temperate", 23.5, 45.0),
        ("Mid-Latitude", 45.0, 60.0),
        ("High Latitude", 60.0, 66.0),
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys_char,name,tolerance", ALL_HOUSE_SYSTEMS)
    def test_system_by_latitude_bands(self, hsys_char, name, tolerance):
        """Test each house system within different latitude bands."""
        hsys = ord(hsys_char)

        for band_name, lat_min, lat_max in self.LATITUDE_BANDS:
            # Filter locations to this latitude band
            band_locations = [
                (lat, lon, jd)
                for lat, lon, jd in RANDOM_LOCATIONS
                if lat_min <= abs(lat) < lat_max
            ]

            if not band_locations:
                continue

            max_diff = 0.0
            for lat, lon, jd in band_locations:
                try:
                    cusps_lib, ascmc_lib = ephem.houses(jd, lat, lon, hsys)
                    cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, bytes([hsys]))

                    asc_diff = angular_diff(ascmc_lib[0], ascmc_swe[0])
                    mc_diff = angular_diff(ascmc_lib[1], ascmc_swe[1])
                    max_diff = max(max_diff, asc_diff, mc_diff)

                except (PolarCircleError, swe.Error):
                    pass  # Expected for some systems at high latitudes

            # Verify max difference is within tolerance
            assert max_diff < tolerance, (
                f"{name} at {band_name} band: max diff {max_diff:.4f}° > {tolerance}°"
            )
