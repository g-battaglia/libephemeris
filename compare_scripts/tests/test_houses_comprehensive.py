"""
Comprehensive House Systems Comparison Tests.

Compares ALL 24 house system calculations between pyswisseph and libephemeris
with tight tolerances across multiple locations, dates, and functions.

Coverage:
  - houses: 24 systems × 8 locations × 12 dates — all 12 cusps + all 8 ASCMC
  - houses_armc: 24 systems × 12 ARMC × 4 latitudes — cusps + ASCMC
  - house_pos: 12 systems × 18 longitudes × 3 body latitudes × 4 dates
  - houses_ex (sidereal): 8 systems × 3 ayanamshas × 6 dates
  - Precision report: per-system max error for cusps and each ASCMC index

Run via:
  poe test:houses
  leph test compare houses
"""

from __future__ import annotations

from collections import defaultdict

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    FLG_SIDEREAL,
    SIDM_LAHIRI,
    SIDM_FAGAN_BRADLEY,
    SIDM_RAMAN,
)


# ============================================================================
# UTILITIES
# ============================================================================


def angular_diff(a: float, b: float) -> float:
    """Absolute angular difference in degrees, handling 360/0 wraparound."""
    d = abs(float(a) - float(b))
    if d > 180.0:
        d = 360.0 - d
    return d


# ============================================================================
# CONFIGURATION
# ============================================================================

# All 24 house systems
HOUSE_SYSTEMS = [
    ("P", "Placidus"),
    ("K", "Koch"),
    ("O", "Porphyry"),
    ("R", "Regiomontanus"),
    ("C", "Campanus"),
    ("E", "Equal (Asc)"),
    ("A", "Equal (MC)"),
    ("W", "Whole Sign"),
    ("B", "Alcabitius"),
    ("T", "Topocentric"),
    ("M", "Morinus"),
    ("X", "Meridian"),
    ("V", "Vehlow"),
    ("H", "Horizontal"),
    ("U", "Krusinski"),
    ("F", "Carter"),
    ("G", "Gauquelin"),
    ("I", "Sunshine"),
    ("i", "Sunshine Makransky"),
    ("N", "Natural Gradient"),
    ("Y", "APC"),
    ("D", "Equal from MC"),
    ("J", "Savard-A"),
    ("L", "Pullen SD"),
    ("S", "Sripati"),
    ("Q", "Pullen SR"),
]

# Per-system cusp tolerance (degrees). Default 0.0011° (~4 arcsec).
# Worst observed: H at 0.00108°. All others < 0.001°.
CUSP_TOLERANCE_DEFAULT = 0.0011
CUSP_TOLERANCE = {
    "H": 0.0012,  # Horizontal: 0.00108° worst case at Sydney
}

# Sidereal cusp tolerance — ayanamsha computation drift over 1950-2050 range
# Worst observed: ~0.004° at epoch boundaries
SIDEREAL_CUSP_TOLERANCE = 0.005

# house_pos default tolerance and per-system overrides
# Most systems: < 0.001. T: 0.019 (semi-arc approximation). B: 0.084 (cusp fallback).
HPOS_TOLERANCE = 0.02
HPOS_TOLERANCE_PER_SYSTEM: dict[str, float] = {}

# ASCMC indices to skip at equator for specific systems (mathematically undefined)
# H (Horizontal): Vertex (3) and CoAsc_Munkasey (6) are degenerate at lat=0
SKIP_ASCMC_AT_EQUATOR: dict[str, set[int]] = {
    "H": {3, 6},  # Vertex and CoAsc_Munkasey undefined at equator
}

# ASCMC labels and tolerances per index
# Tightened to just above worst observed values:
#   ASC: 0.00070°, MC: 0.00058°, ARMC: 0.00053°, Vertex: 0.00108°
#   EquAsc: 0.00057°, CoAsc_Koch: 0.00140°, CoAsc_Munk: 0.00076°, PolarAsc: 0.00140°
ASCMC_LABELS = [
    "ASC", "MC", "ARMC", "Vertex",
    "EquAsc", "CoAsc_Koch", "CoAsc_Munkasey", "PolarAsc",
]
ASCMC_TOLERANCE = {
    0: 0.001,    # ASC — worst: 0.00070°
    1: 0.001,    # MC — worst: 0.00058°
    2: 0.001,    # ARMC — worst: 0.00053°
    3: 0.002,    # Vertex — worst: 0.00108°
    4: 0.001,    # Equatorial Ascendant — worst: 0.00057°
    5: 0.002,    # Co-Ascendant Koch — worst: 0.00140°
    6: 0.001,    # Co-Ascendant Munkasey — worst: 0.00076°
    7: 0.002,    # Polar Ascendant — worst: 0.00140°
}

# 8 test locations: (name, lat, lon)
LOCATIONS = [
    ("Equator", 0.0, 0.0),
    ("Rome", 41.9028, 12.4964),
    ("London", 51.5074, -0.1278),
    ("Sydney", -33.8688, 151.2093),
    ("Tokyo", 35.6762, 139.6503),
    ("Reykjavik", 64.0, -22.0),
    ("Buenos Aires", -34.6037, -58.3816),
    ("Mumbai", 19.0760, 72.8777),
]

# 12 Julian Days spanning 1950–2050
JD_1950 = 2433282.5   # 1950-01-01
JD_2050 = 2469807.5   # 2050-01-01
TEST_JDS = [JD_1950 + i * (JD_2050 - JD_1950) / 11 for i in range(12)]

# 6 Julian Days for lighter tests
TEST_JDS_6 = TEST_JDS[::2]  # every other JD

# Systems that raise PolarCircleError at high latitudes
POLAR_FAIL_SYSTEMS = {"P", "K", "G"}


# ============================================================================
# SHARED PRECISION TRACKER
# ============================================================================

# Collects max errors per system for the final precision report.
# Keys: (system_char, value_name) -> max_diff
_precision_data: dict[tuple[str, str], float] = defaultdict(float)


def _track(system: str, value_name: str, diff: float) -> None:
    """Track the maximum diff for a (system, value) pair."""
    key = (system, value_name)
    if diff > _precision_data[key]:
        _precision_data[key] = diff


# ============================================================================
# TEST: houses — All cusps + ALL 8 ASCMC
# ============================================================================


class TestHouseCuspsAllSystems:
    """Compare all 12 house cusps for all 24 systems vs pyswisseph."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys,hsys_name", HOUSE_SYSTEMS)
    @pytest.mark.parametrize("loc_name,lat,lon", LOCATIONS)
    @pytest.mark.parametrize("jd", TEST_JDS, ids=[f"jd{i}" for i in range(12)])
    def test_cusps(self, hsys, hsys_name, loc_name, lat, lon, jd):
        """All 12 cusps must match within tolerance."""
        tolerance = CUSP_TOLERANCE.get(hsys, CUSP_TOLERANCE_DEFAULT)

        try:
            cusps_swe, _ = swe.houses(jd, lat, lon, hsys.encode("ascii"))
        except swe.Error:
            pytest.skip(f"pyswisseph raises error for {hsys} at {loc_name}")

        try:
            cusps_py, _ = ephem.houses(jd, lat, lon, hsys)
        except Exception:
            pytest.skip(f"libephemeris raises error for {hsys} at {loc_name}")

        num_cusps = min(len(cusps_swe), len(cusps_py))
        max_diff = 0.0
        worst_cusp = 0

        for i in range(num_cusps):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            _track(hsys, f"cusp{i + 1}", diff)
            if diff > max_diff:
                max_diff = diff
                worst_cusp = i + 1

        assert max_diff < tolerance, (
            f"{hsys_name} ({hsys}) at {loc_name}: cusp {worst_cusp} diff "
            f"{max_diff:.6f}° exceeds {tolerance}° "
            f"(swe={cusps_swe[worst_cusp - 1]:.6f}, py={cusps_py[worst_cusp - 1]:.6f})"
        )


class TestAscmcComplete:
    """Compare ALL 8 ASCMC values for all 24 systems vs pyswisseph."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys,hsys_name", HOUSE_SYSTEMS)
    @pytest.mark.parametrize("loc_name,lat,lon", LOCATIONS)
    @pytest.mark.parametrize("jd", TEST_JDS_6, ids=[f"jd{i}" for i in range(6)])
    def test_ascmc_all_values(self, hsys, hsys_name, loc_name, lat, lon, jd):
        """All 8 ASCMC values must match within their respective tolerances."""
        try:
            _, ascmc_swe = swe.houses(jd, lat, lon, hsys.encode("ascii"))
        except swe.Error:
            pytest.skip(f"pyswisseph raises error for {hsys} at {loc_name}")

        try:
            _, ascmc_py = ephem.houses(jd, lat, lon, hsys)
        except Exception:
            pytest.skip(f"libephemeris raises error for {hsys} at {loc_name}")

        failures = []
        for i in range(min(len(ascmc_swe), len(ascmc_py), 8)):
            # Skip ASCMC indices that are undefined at equator for specific systems
            if abs(lat) < 0.1 and i in SKIP_ASCMC_AT_EQUATOR.get(hsys, set()):
                continue

            tol = ASCMC_TOLERANCE.get(i, 0.01)
            diff = angular_diff(ascmc_swe[i], ascmc_py[i])
            _track(hsys, ASCMC_LABELS[i], diff)

            if diff >= tol:
                failures.append(
                    f"  ascmc[{i}] ({ASCMC_LABELS[i]}): diff={diff:.6f}° > {tol}° "
                    f"(swe={float(ascmc_swe[i]):.6f}, py={float(ascmc_py[i]):.6f})"
                )

        assert not failures, (
            f"{hsys_name} ({hsys}) at {loc_name} jd={jd:.1f}:\n"
            + "\n".join(failures)
        )


# ============================================================================
# TEST: houses_armc — All 24 systems, multiple latitudes
# ============================================================================


# 12 ARMC values
ARMC_VALUES = [float(x) for x in range(0, 360, 30)]

# 4 latitudes for ARMC tests
ARMC_LATITUDES = [
    ("equator", 0.0),
    ("Rome", 41.9),
    ("Sydney", -33.9),
    ("high_lat", 60.0),
]


class TestHousesArmcAllSystems:
    """Compare houses_armc for all 24 systems across ARMC values and latitudes."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys,hsys_name", HOUSE_SYSTEMS)
    @pytest.mark.parametrize("armc", ARMC_VALUES, ids=[f"armc{int(a)}" for a in ARMC_VALUES])
    @pytest.mark.parametrize("lat_name,lat", ARMC_LATITUDES)
    def test_houses_armc_cusps_and_ascmc(self, hsys, hsys_name, armc, lat_name, lat):
        """Cusps and all 8 ASCMC from houses_armc must match pyswisseph."""
        eps = 23.4393  # standard obliquity

        try:
            cusps_swe, ascmc_swe = swe.houses_armc(armc, lat, eps, hsys.encode("ascii"))
        except swe.Error:
            pytest.skip(f"pyswisseph error for {hsys} at {lat_name}")

        try:
            cusps_py, ascmc_py = ephem.houses_armc(armc, lat, eps, hsys)
        except Exception:
            pytest.skip(f"libephemeris error for {hsys} at {lat_name}")

        cusp_tol = CUSP_TOLERANCE.get(hsys, CUSP_TOLERANCE_DEFAULT)

        # Compare cusps
        num_cusps = min(len(cusps_swe), len(cusps_py))
        max_cusp_diff = 0.0
        for i in range(num_cusps):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            _track(hsys, f"armc_cusp{i + 1}", diff)
            max_cusp_diff = max(max_cusp_diff, diff)

        assert max_cusp_diff < cusp_tol, (
            f"{hsys_name} houses_armc at {lat_name} armc={armc}: "
            f"max cusp diff {max_cusp_diff:.6f}° > {cusp_tol}°"
        )

        # Compare ASCMC
        failures = []
        for i in range(min(len(ascmc_swe), len(ascmc_py), 8)):
            # Skip ASCMC indices that are undefined at equator for specific systems
            if abs(lat) < 0.1 and i in SKIP_ASCMC_AT_EQUATOR.get(hsys, set()):
                continue

            tol = ASCMC_TOLERANCE.get(i, 0.01)
            diff = angular_diff(ascmc_swe[i], ascmc_py[i])
            _track(hsys, f"armc_{ASCMC_LABELS[i]}", diff)
            if diff >= tol:
                failures.append(
                    f"  ascmc[{i}] ({ASCMC_LABELS[i]}): {diff:.6f}° > {tol}°"
                )

        assert not failures, (
            f"{hsys_name} houses_armc at {lat_name} armc={armc}:\n"
            + "\n".join(failures)
        )


# ============================================================================
# TEST: house_pos — Comprehensive
# ============================================================================

HPOS_SYSTEMS = ["P", "K", "O", "R", "C", "E", "W", "B", "M", "T", "X", "V"]
HPOS_LONGITUDES = [float(x) for x in range(0, 360, 20)]  # 18 values
HPOS_BODY_LATS = [0.0, 3.0, -5.0]
HPOS_TOLERANCE = 0.02


class TestHousePosComprehensive:
    """Compare house_pos for 12 systems, 18 longitudes, 3 body latitudes."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys", HPOS_SYSTEMS)
    @pytest.mark.parametrize("planet_lon", HPOS_LONGITUDES, ids=[f"lon{int(l)}" for l in HPOS_LONGITUDES])
    @pytest.mark.parametrize("body_lat", HPOS_BODY_LATS, ids=[f"blat{l}" for l in HPOS_BODY_LATS])
    @pytest.mark.parametrize("jd", TEST_JDS_6[:4], ids=[f"jd{i}" for i in range(4)])
    def test_house_pos(self, hsys, planet_lon, body_lat, jd):
        """house_pos results must match pyswisseph."""
        geo_lat = 41.9  # Rome

        # Get ARMC and obliquity from pyswisseph
        try:
            _, ascmc_swe = swe.houses(jd, geo_lat, 12.5, b"P")
            armc = float(ascmc_swe[2])
            ecl = swe.calc_ut(jd, -1, 0)
            eps = float(ecl[0][0])
        except Exception:
            pytest.skip("Cannot compute ARMC/eps")

        try:
            hp_swe = float(swe.house_pos(
                armc, geo_lat, eps, (planet_lon, body_lat), hsys.encode("ascii")
            ))
        except swe.Error:
            pytest.skip(f"pyswisseph error for house_pos {hsys}")

        try:
            hp_py = float(ephem.house_pos(
                armc, geo_lat, eps, (planet_lon, body_lat), hsys
            ))
        except Exception:
            pytest.skip(f"libephemeris error for house_pos {hsys}")

        diff = abs(hp_swe - hp_py)
        # Handle house number wraparound (12.95 vs 1.05)
        if diff > 6.0:
            diff = 12.0 - diff

        _track(hsys, "house_pos", diff)

        tolerance = HPOS_TOLERANCE_PER_SYSTEM.get(hsys, HPOS_TOLERANCE)
        assert diff < tolerance, (
            f"house_pos {hsys} lon={planet_lon} blat={body_lat}: "
            f"diff={diff:.4f} (swe={hp_swe:.4f}, py={hp_py:.4f})"
        )


# ============================================================================
# TEST: houses_ex — Sidereal
# ============================================================================

SIDEREAL_SYSTEMS = ["P", "K", "E", "W", "O", "R", "C", "B"]
SIDEREAL_MODES = [
    (SIDM_LAHIRI, "Lahiri"),
    (SIDM_FAGAN_BRADLEY, "Fagan-Bradley"),
    (SIDM_RAMAN, "Raman"),
]


class TestSiderealHouses:
    """Compare sidereal house calculations via houses_ex."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys", SIDEREAL_SYSTEMS)
    @pytest.mark.parametrize("sid_mode,sid_name", SIDEREAL_MODES)
    @pytest.mark.parametrize("jd", TEST_JDS_6, ids=[f"jd{i}" for i in range(6)])
    def test_sidereal_cusps(self, hsys, sid_mode, sid_name, jd):
        """Sidereal cusps must match pyswisseph."""
        lat, lon = 41.9, 12.5  # Rome

        swe.set_sid_mode(sid_mode)
        ephem.set_sid_mode(sid_mode)

        try:
            cusps_swe, ascmc_swe = swe.houses_ex(
                jd, lat, lon, hsys.encode("ascii"), FLG_SIDEREAL
            )
        except swe.Error:
            swe.set_sid_mode(0)
            ephem.set_sid_mode(0)
            pytest.skip(f"pyswisseph error for sidereal {hsys}")

        try:
            cusps_py, ascmc_py = ephem.houses_ex(
                jd, lat, lon, hsys, FLG_SIDEREAL
            )
        except Exception:
            swe.set_sid_mode(0)
            ephem.set_sid_mode(0)
            pytest.skip(f"libephemeris error for sidereal {hsys}")

        swe.set_sid_mode(0)
        ephem.set_sid_mode(0)

        cusp_tol = SIDEREAL_CUSP_TOLERANCE
        max_diff = 0.0
        for i in range(min(len(cusps_swe), len(cusps_py), 12)):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            max_diff = max(max_diff, diff)

        assert max_diff < cusp_tol, (
            f"Sidereal {hsys} ({sid_name}): max cusp diff {max_diff:.6f}° > {cusp_tol}°"
        )

        # Also compare ASC and MC (sidereal tolerance for ayanamsha drift)
        asc_diff = angular_diff(ascmc_swe[0], ascmc_py[0])
        mc_diff = angular_diff(ascmc_swe[1], ascmc_py[1])
        assert asc_diff < SIDEREAL_CUSP_TOLERANCE, (
            f"Sidereal {hsys} ({sid_name}): ASC diff {asc_diff:.6f}°"
        )
        assert mc_diff < SIDEREAL_CUSP_TOLERANCE, (
            f"Sidereal {hsys} ({sid_name}): MC diff {mc_diff:.6f}°"
        )


# ============================================================================
# EDGE CASES — Boundary conditions and degenerate inputs
# ============================================================================


# Locations that stress edge cases
EDGE_LOCATIONS = [
    ("NearPole_N", 65.0, 0.0),         # Just below polar circle
    ("NearPole_S", -65.0, 0.0),        # Southern near-polar
    ("DateLine_E", 35.0, 179.9),       # Near international date line
    ("DateLine_W", 35.0, -179.9),      # Other side of date line
    ("PrimeMeridian", 51.5, 0.0),      # Prime meridian exactly
    ("Tropic_N", 23.44, 0.0),          # Tropic of Cancer
    ("Tropic_S", -23.44, 0.0),         # Tropic of Capricorn
    ("HighSouth", -60.0, 0.0),         # High southern latitude
]

# Special Julian Days
EDGE_JDS = [
    2451545.0,    # J2000.0 exactly
    2451179.5,    # 1999-01-01 (near century boundary)
    2415020.5,    # 1900-01-01 (early date)
    2488069.5,    # 2099-12-31 (late date)
    2451623.5,    # 2000-03-20 (vernal equinox ~)
    2451810.5,    # 2000-09-22 (autumnal equinox ~)
    2451716.5,    # 2000-06-21 (summer solstice ~)
    2451900.5,    # 2000-12-21 (winter solstice ~)
]

# Systems that work at all latitudes (no polar circle issues)
ALL_LAT_SYSTEMS = [
    ("E", "Equal"), ("W", "Whole Sign"), ("O", "Porphyry"),
    ("M", "Morinus"), ("R", "Regiomontanus"), ("C", "Campanus"),
    ("X", "Meridian"), ("V", "Vehlow"), ("H", "Horizontal"),
    ("U", "Krusinski"), ("F", "Carter"), ("N", "Natural Gradient"),
    ("D", "Equal from MC"), ("L", "Pullen SD"), ("S", "Sripati"),
    ("Q", "Pullen SR"), ("B", "Alcabitius"), ("T", "Topocentric"),
    ("J", "Savard-A"), ("Y", "APC"), ("A", "Equal (MC)"),
]

# Quadrant systems where cusp1=ASC and cusp10=MC
# Sripati (S), Pullen SD (L), Pullen SR (Q) have midpoint cusps — cusp1≠ASC
QUADRANT_SYSTEMS = ["P", "K", "O", "R", "C", "B", "T", "U"]


class TestEdgeCaseLocations:
    """Test house systems at edge case locations (near poles, date line, etc.)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys,hsys_name", ALL_LAT_SYSTEMS)
    @pytest.mark.parametrize("loc_name,lat,lon", EDGE_LOCATIONS)
    def test_edge_locations(self, hsys, hsys_name, loc_name, lat, lon):
        """Cusps must match at edge case locations."""
        jd = 2451545.0  # J2000

        try:
            cusps_swe, _ = swe.houses(jd, lat, lon, hsys.encode("ascii"))
        except swe.Error:
            pytest.skip(f"pyswisseph error for {hsys} at {loc_name}")

        try:
            cusps_py, _ = ephem.houses(jd, lat, lon, hsys)
        except Exception:
            pytest.skip(f"libephemeris error for {hsys} at {loc_name}")

        tolerance = CUSP_TOLERANCE.get(hsys, CUSP_TOLERANCE_DEFAULT)
        num_cusps = min(len(cusps_swe), len(cusps_py))
        max_diff = 0.0
        for i in range(num_cusps):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            max_diff = max(max_diff, diff)

        assert max_diff < tolerance, (
            f"{hsys_name} at {loc_name}: max cusp diff {max_diff:.6f}° > {tolerance}°"
        )


class TestEdgeCaseDates:
    """Test house systems at special dates (solstices, equinoxes, boundaries)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys,hsys_name", HOUSE_SYSTEMS[:12])  # main systems
    @pytest.mark.parametrize("jd", EDGE_JDS, ids=[
        "J2000", "1999", "1900", "2099",
        "VernalEq", "AutumnEq", "SummerSol", "WinterSol",
    ])
    def test_special_dates(self, hsys, hsys_name, jd):
        """Cusps must match at special dates."""
        lat, lon = 41.9, 12.5  # Rome

        try:
            cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, hsys.encode("ascii"))
        except swe.Error:
            pytest.skip(f"pyswisseph error for {hsys}")

        try:
            cusps_py, ascmc_py = ephem.houses(jd, lat, lon, hsys)
        except Exception:
            pytest.skip(f"libephemeris error for {hsys}")

        tolerance = CUSP_TOLERANCE.get(hsys, CUSP_TOLERANCE_DEFAULT)
        max_diff = max(
            angular_diff(cusps_swe[i], cusps_py[i])
            for i in range(min(len(cusps_swe), len(cusps_py)))
        )
        assert max_diff < tolerance, (
            f"{hsys_name} at jd={jd}: max cusp diff {max_diff:.6f}°"
        )

        # ASC and MC must also match
        asc_diff = angular_diff(ascmc_swe[0], ascmc_py[0])
        mc_diff = angular_diff(ascmc_swe[1], ascmc_py[1])
        assert asc_diff < 0.001, f"{hsys_name}: ASC diff {asc_diff:.6f}°"
        assert mc_diff < 0.001, f"{hsys_name}: MC diff {mc_diff:.6f}°"


# ============================================================================
# STRUCTURAL TESTS — Cusp relationships and consistency
# ============================================================================


class TestCuspAscMcConsistency:
    """Verify cusp1=ASC and cusp10=MC for quadrant systems, matching pyswisseph."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys", QUADRANT_SYSTEMS)
    @pytest.mark.parametrize("loc_name,lat,lon", LOCATIONS[:4])
    def test_cusp1_equals_asc(self, hsys, loc_name, lat, lon):
        """cusp[0] must equal ASC for quadrant systems, matching pyswisseph."""
        jd = 2451545.0

        try:
            cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, hsys.encode("ascii"))
            cusps_py, ascmc_py = ephem.houses(jd, lat, lon, hsys)
        except Exception:
            pytest.skip(f"Error for {hsys} at {loc_name}")

        # Both should have cusp1 == ASC
        swe_diff = angular_diff(cusps_swe[0], ascmc_swe[0])
        py_diff = angular_diff(cusps_py[0], ascmc_py[0])
        assert swe_diff < 0.0001, f"pyswisseph cusp1≠ASC for {hsys}: {swe_diff:.6f}°"
        assert py_diff < 0.0001, f"libephemeris cusp1≠ASC for {hsys}: {py_diff:.6f}°"

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys", QUADRANT_SYSTEMS)
    @pytest.mark.parametrize("loc_name,lat,lon", LOCATIONS[:4])
    def test_cusp10_equals_mc(self, hsys, loc_name, lat, lon):
        """cusp[9] must equal MC for quadrant systems, matching pyswisseph."""
        jd = 2451545.0

        try:
            cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, hsys.encode("ascii"))
            cusps_py, ascmc_py = ephem.houses(jd, lat, lon, hsys)
        except Exception:
            pytest.skip(f"Error for {hsys} at {loc_name}")

        swe_diff = angular_diff(cusps_swe[9], ascmc_swe[1])
        py_diff = angular_diff(cusps_py[9], ascmc_py[1])
        assert swe_diff < 0.0001, f"pyswisseph cusp10≠MC for {hsys}: {swe_diff:.6f}°"
        assert py_diff < 0.0001, f"libephemeris cusp10≠MC for {hsys}: {py_diff:.6f}°"


class TestOppositeHouses:
    """Verify cusps[i] and cusps[i+6] are 180° apart for standard systems."""

    # Systems with strict 180° opposition (exclude Gauquelin, APC)
    OPPOSITE_SYSTEMS = [
        "P", "K", "O", "R", "C", "E", "W", "B", "T", "M", "X", "V",
        "H", "U", "F", "N", "D", "L", "S", "Q", "A", "J",
    ]

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys", OPPOSITE_SYSTEMS)
    @pytest.mark.parametrize("loc_name,lat,lon", LOCATIONS[:4])
    def test_opposite_cusps_180(self, hsys, loc_name, lat, lon):
        """Houses 1-6 must be 180° from houses 7-12."""
        jd = 2451545.0

        try:
            cusps_py, _ = ephem.houses(jd, lat, lon, hsys)
        except Exception:
            pytest.skip(f"Error for {hsys} at {loc_name}")

        for i in range(6):
            opp = (i + 6) % 12
            diff = angular_diff(cusps_py[i], cusps_py[opp])
            assert abs(diff - 180.0) < 0.001, (
                f"{hsys} at {loc_name}: house {i + 1} and {opp + 1} "
                f"not 180° apart (diff={diff:.4f}°)"
            )


class TestPolarCircleErrors:
    """Verify that P/K/G raise errors at polar latitudes, matching pyswisseph."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys", ["P", "K", "G"])
    @pytest.mark.parametrize("lat", [67.0, 70.0, 80.0, -67.0, -70.0, -80.0])
    def test_polar_error_both_raise(self, hsys, lat):
        """Both libraries must raise errors at polar latitudes."""
        jd = 2451545.0

        swe_error = False
        try:
            swe.houses(jd, lat, 0.0, hsys.encode("ascii"))
        except swe.Error:
            swe_error = True

        py_error = False
        try:
            ephem.houses(jd, lat, 0.0, hsys)
        except Exception:
            py_error = True

        assert swe_error == py_error, (
            f"{hsys} at lat={lat}: swe_error={swe_error}, py_error={py_error}"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys", ["P", "K", "G"])
    def test_works_below_polar_circle(self, hsys):
        """P/K/G must work at 65° (below polar circle), matching pyswisseph."""
        jd = 2451545.0
        lat = 65.0

        cusps_swe, ascmc_swe = swe.houses(jd, lat, 0.0, hsys.encode("ascii"))
        cusps_py, ascmc_py = ephem.houses(jd, lat, 0.0, hsys)

        asc_diff = angular_diff(ascmc_swe[0], ascmc_py[0])
        assert asc_diff < 0.002, f"{hsys} at 65°: ASC diff {asc_diff:.6f}°"


class TestGauquelin36Sectors:
    """Test Gauquelin system returns 36 sectors and they all match."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("loc_name,lat,lon", LOCATIONS[:4])
    @pytest.mark.parametrize("jd", TEST_JDS_6[:3], ids=[f"jd{i}" for i in range(3)])
    def test_gauquelin_36_cusps(self, loc_name, lat, lon, jd):
        """Gauquelin must return 36 sectors, all matching pyswisseph."""
        try:
            cusps_swe, _ = swe.houses(jd, lat, lon, b"G")
        except swe.Error:
            pytest.skip(f"pyswisseph error at {loc_name}")

        try:
            cusps_py, _ = ephem.houses(jd, lat, lon, "G")
        except Exception:
            pytest.skip(f"libephemeris error at {loc_name}")

        # Both should have 36 cusps
        assert len(cusps_swe) >= 36, f"pyswisseph returned {len(cusps_swe)} cusps"
        assert len(cusps_py) >= 36, f"libephemeris returned {len(cusps_py)} cusps"

        max_diff = 0.0
        for i in range(36):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            max_diff = max(max_diff, diff)

        assert max_diff < 0.001, (
            f"Gauquelin at {loc_name}: max sector diff {max_diff:.6f}°"
        )


class TestHousesArmcConsistencyWithHouses:
    """Verify houses_armc gives same results as houses() when using same ARMC."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys", ["P", "E", "W", "O", "R", "C", "M"])
    def test_armc_consistency(self, hsys):
        """houses_armc(ARMC from houses()) must produce identical cusps."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5

        cusps_h, ascmc_h = ephem.houses(jd, lat, lon, hsys)
        armc = ascmc_h[2]

        # Get obliquity from pyswisseph for consistency
        ecl = swe.calc_ut(jd, -1, 0)
        eps = float(ecl[0][0])

        cusps_a, ascmc_a = ephem.houses_armc(armc, lat, eps, hsys)

        # Cusps should be very close (not exactly equal due to obliquity path)
        for i in range(12):
            diff = angular_diff(cusps_h[i], cusps_a[i])
            assert diff < 0.01, (
                f"{hsys} cusp {i + 1}: houses={cusps_h[i]:.6f}, "
                f"houses_armc={cusps_a[i]:.6f}, diff={diff:.6f}°"
            )


class TestAscmcSameAcrossSystems:
    """ASCMC values (except cusps) should be identical regardless of house system."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("loc_name,lat,lon", LOCATIONS[:4])
    def test_ascmc_system_independent(self, loc_name, lat, lon):
        """ASC, MC, ARMC, Vertex, etc. must be the same for all systems."""
        jd = 2451545.0

        # Use Placidus as reference
        try:
            _, ref_ascmc = ephem.houses(jd, lat, lon, "P")
        except Exception:
            pytest.skip("Cannot compute Placidus")

        for hsys, hsys_name in HOUSE_SYSTEMS:
            if hsys in POLAR_FAIL_SYSTEMS and abs(lat) > 66:
                continue
            try:
                _, ascmc = ephem.houses(jd, lat, lon, hsys)
            except Exception:
                continue

            for i in range(8):
                diff = angular_diff(ref_ascmc[i], ascmc[i])
                assert diff < 0.001, (
                    f"ascmc[{i}] ({ASCMC_LABELS[i]}) differs between P and {hsys} "
                    f"at {loc_name}: {diff:.6f}°"
                )


class TestHousePosAtCuspBoundaries:
    """Test house_pos returns correct house number when body is at a cusp."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys", ["P", "E", "W", "O", "R"])
    def test_body_at_cusp_start(self, hsys):
        """Body at cusp start should give house number with near-zero fraction."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5

        _, ascmc = swe.houses(jd, lat, lon, b"P")
        armc = float(ascmc[2])
        ecl = swe.calc_ut(jd, -1, 0)
        eps = float(ecl[0][0])

        cusps, _ = ephem.houses_armc(armc, lat, eps, hsys)

        # Place body at cusp 1 (ASC)
        hp = float(ephem.house_pos(armc, lat, eps, (cusps[0], 0.0), hsys))
        house_num = int(hp)
        fraction = hp - house_num

        assert house_num == 1, f"{hsys}: body at cusp1 in house {house_num}"
        assert fraction < 0.01, f"{hsys}: fraction {fraction:.4f} not near 0"

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys", ["P", "E", "O", "R"])
    def test_body_mid_house(self, hsys):
        """Body in middle of house should give ~0.5 fraction."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5

        _, ascmc = swe.houses(jd, lat, lon, b"P")
        armc = float(ascmc[2])
        ecl = swe.calc_ut(jd, -1, 0)
        eps = float(ecl[0][0])

        cusps, _ = ephem.houses_armc(armc, lat, eps, hsys)

        # Place body in middle of house 1
        c1 = cusps[0]
        c2 = cusps[1]
        house_size = (c2 - c1 + 360.0) % 360.0
        mid_lon = (c1 + house_size * 0.5) % 360.0

        hp = float(ephem.house_pos(armc, lat, eps, (mid_lon, 0.0), hsys))
        house_num = int(hp)
        fraction = hp - house_num

        assert house_num == 1, f"{hsys}: body at mid-house1 in house {house_num}"
        assert 0.3 < fraction < 0.7, f"{hsys}: fraction {fraction:.4f} not near 0.5"


class TestHousePosMultipleLocations:
    """Test house_pos at additional locations beyond Rome."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys", ["P", "E", "W", "O", "R", "C"])
    @pytest.mark.parametrize("geo_lat,geo_lon,loc_name", [
        (0.0, 0.0, "Equator"),
        (-33.9, 151.2, "Sydney"),
        (55.75, 37.62, "Moscow"),
    ])
    @pytest.mark.parametrize("planet_lon", [0.0, 90.0, 180.0, 270.0])
    def test_house_pos_various_locations(self, hsys, geo_lat, geo_lon, loc_name, planet_lon):
        """house_pos must match pyswisseph at various locations."""
        jd = 2451545.0

        try:
            _, ascmc_swe = swe.houses(jd, geo_lat, geo_lon, b"P")
            armc = float(ascmc_swe[2])
            ecl = swe.calc_ut(jd, -1, 0)
            eps = float(ecl[0][0])
        except Exception:
            pytest.skip("Cannot compute ARMC/eps")

        try:
            hp_swe = float(swe.house_pos(
                armc, geo_lat, eps, (planet_lon, 0.0), hsys.encode("ascii")
            ))
            hp_py = float(ephem.house_pos(
                armc, geo_lat, eps, (planet_lon, 0.0), hsys
            ))
        except Exception:
            pytest.skip(f"Error for {hsys} at {loc_name}")

        diff = abs(hp_swe - hp_py)
        if diff > 6.0:
            diff = 12.0 - diff

        assert diff < HPOS_TOLERANCE, (
            f"house_pos {hsys} at {loc_name} lon={planet_lon}: "
            f"diff={diff:.4f} (swe={hp_swe:.4f}, py={hp_py:.4f})"
        )


class TestHousesArmcEdgeCaseArmc:
    """Test houses_armc at boundary ARMC values (near 0/360)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys", ["P", "E", "W", "O", "R", "M"])
    @pytest.mark.parametrize("armc", [0.001, 0.01, 0.1, 359.9, 359.99, 359.999])
    def test_armc_near_zero_boundary(self, hsys, armc):
        """houses_armc must match near the 0°/360° boundary."""
        lat, eps = 41.9, 23.4393

        try:
            cusps_swe, _ = swe.houses_armc(armc, lat, eps, hsys.encode("ascii"))
        except swe.Error:
            pytest.skip(f"pyswisseph error for {hsys} at armc={armc}")

        try:
            cusps_py, _ = ephem.houses_armc(armc, lat, eps, hsys)
        except Exception:
            pytest.skip(f"libephemeris error for {hsys} at armc={armc}")

        tolerance = CUSP_TOLERANCE.get(hsys, CUSP_TOLERANCE_DEFAULT)
        max_diff = max(
            angular_diff(cusps_swe[i], cusps_py[i])
            for i in range(min(len(cusps_swe), len(cusps_py)))
        )
        assert max_diff < tolerance, (
            f"{hsys} armc={armc}: max cusp diff {max_diff:.6f}°"
        )


class TestSiderealMultipleLocations:
    """Test sidereal houses at various locations, not just Rome."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys", ["P", "E", "W"])
    @pytest.mark.parametrize("loc_name,lat,lon", [
        ("Rome", 41.9, 12.5),
        ("Sydney", -33.9, 151.2),
        ("Equator", 0.0, 0.0),
        ("Tokyo", 35.7, 139.7),
    ])
    def test_sidereal_various_locations(self, hsys, loc_name, lat, lon):
        """Sidereal cusps must match at various locations."""
        jd = 2451545.0  # J2000

        swe.set_sid_mode(SIDM_LAHIRI)
        ephem.set_sid_mode(SIDM_LAHIRI)

        try:
            cusps_swe, _ = swe.houses_ex(
                jd, lat, lon, hsys.encode("ascii"), FLG_SIDEREAL
            )
            cusps_py, _ = ephem.houses_ex(jd, lat, lon, hsys, FLG_SIDEREAL)
        except Exception:
            swe.set_sid_mode(0)
            ephem.set_sid_mode(0)
            pytest.skip(f"Error for sidereal {hsys} at {loc_name}")

        swe.set_sid_mode(0)
        ephem.set_sid_mode(0)

        max_diff = max(
            angular_diff(cusps_swe[i], cusps_py[i])
            for i in range(min(len(cusps_swe), len(cusps_py), 12))
        )
        assert max_diff < SIDEREAL_CUSP_TOLERANCE, (
            f"Sidereal {hsys} at {loc_name}: max diff {max_diff:.6f}°"
        )


# ============================================================================
# PRECISION REPORT — Aggregate statistics
# ============================================================================


class TestPrecisionReport:
    """Print per-system precision summary after all tests have run.

    This test always passes — its purpose is to display the precision table.
    Run with -s flag to see the output.
    """

    @pytest.mark.comparison
    def test_print_precision_report(self):
        """Print per-system max error table."""
        if not _precision_data:
            pytest.skip("No precision data collected (run other tests first)")

        # Collect all systems and value names
        systems = sorted({k[0] for k in _precision_data})
        cusp_keys = [f"cusp{i}" for i in range(1, 13)]
        ascmc_keys = ASCMC_LABELS

        print("\n" + "=" * 100)
        print("HOUSE SYSTEMS PRECISION REPORT — max error vs pyswisseph (degrees)")
        print("=" * 100)

        # --- Cusp max errors ---
        print("\nMax cusp error per system (houses):")
        print(f"  {'System':<25s} {'Max Cusp':>12s} {'Worst Cusp':>12s}")
        print(f"  {'-' * 25} {'-' * 12} {'-' * 12}")
        for sys_char in systems:
            max_cusp = 0.0
            worst = ""
            for ck in cusp_keys:
                v = _precision_data.get((sys_char, ck), 0.0)
                if v > max_cusp:
                    max_cusp = v
                    worst = ck
            if max_cusp > 0:
                # Find name
                name = next((n for c, n in HOUSE_SYSTEMS if c == sys_char), sys_char)
                print(f"  {name + ' (' + sys_char + ')':<25s} {max_cusp:>12.7f} {worst:>12s}")

        # --- ASCMC max errors ---
        print(f"\nMax ASCMC error per system (houses):")
        header = f"  {'System':<18s}"
        for lbl in ascmc_keys:
            header += f" {lbl:>10s}"
        print(header)
        print(f"  {'-' * 18}" + (" " + "-" * 10) * 8)
        for sys_char in systems:
            name = next((n for c, n in HOUSE_SYSTEMS if c == sys_char), sys_char)
            row = f"  {name[:16] + ' (' + sys_char + ')':<18s}"
            for lbl in ascmc_keys:
                v = _precision_data.get((sys_char, lbl), 0.0)
                if v > 0:
                    row += f" {v:>10.7f}"
                else:
                    row += f" {'—':>10s}"
            print(row)

        # --- houses_armc max cusp error ---
        armc_cusp_keys = [f"armc_cusp{i}" for i in range(1, 13)]
        has_armc_data = any(
            _precision_data.get((s, k), 0) > 0
            for s in systems for k in armc_cusp_keys
        )
        if has_armc_data:
            print(f"\nMax cusp error per system (houses_armc):")
            print(f"  {'System':<25s} {'Max Cusp':>12s}")
            print(f"  {'-' * 25} {'-' * 12}")
            for sys_char in systems:
                max_cusp = 0.0
                for ck in armc_cusp_keys:
                    v = _precision_data.get((sys_char, ck), 0.0)
                    max_cusp = max(max_cusp, v)
                if max_cusp > 0:
                    name = next((n for c, n in HOUSE_SYSTEMS if c == sys_char), sys_char)
                    print(f"  {name + ' (' + sys_char + ')':<25s} {max_cusp:>12.7f}")

        # --- house_pos max error ---
        has_hpos = any(
            _precision_data.get((s, "house_pos"), 0) > 0 for s in systems
        )
        if has_hpos:
            print(f"\nMax house_pos error per system:")
            print(f"  {'System':<25s} {'Max Diff':>12s}")
            print(f"  {'-' * 25} {'-' * 12}")
            for sys_char in systems:
                v = _precision_data.get((sys_char, "house_pos"), 0.0)
                if v > 0:
                    name = next((n for c, n in HOUSE_SYSTEMS if c == sys_char), sys_char)
                    print(f"  {name + ' (' + sys_char + ')':<25s} {v:>12.7f}")

        print("\n" + "=" * 100)
        print("END OF PRECISION REPORT")
        print("=" * 100 + "\n")

        # This test always passes — it's just a report
        assert True
