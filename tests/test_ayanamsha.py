import pytest
import swisseph as swe
import libephemeris as pyephem
from libephemeris.constants import *

# Famous People for Ayanamsha Testing (born after 1899 for DE421 range)
FAMOUS_PEOPLE_AYANAMSHA = [
    ("Frida Kahlo", 1907, 7, 6, 8.5),
    ("Nelson Mandela", 1918, 7, 18, 14.0),
    ("Marilyn Monroe", 1926, 6, 1, 9.5),
    ("Queen Elizabeth II", 1926, 4, 21, 2.66),
    ("Barack Obama", 1961, 8, 4, 19.4),
]

# All 43 Ayanamsha modes
AYANAMSHA_MODES = [
    (SE_SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
    (SE_SIDM_LAHIRI, "Lahiri"),
    (SE_SIDM_DELUCE, "De Luce"),
    (SE_SIDM_RAMAN, "Raman"),
    (SE_SIDM_USHASHASHI, "Ushashashi"),
    (SE_SIDM_KRISHNAMURTI, "Krishnamurti"),
    (SE_SIDM_DJWHAL_KHUL, "Djwhal Khul"),
    (SE_SIDM_YUKTESHWAR, "Yukteshwar"),
    (SE_SIDM_JN_BHASIN, "JN Bhasin"),
    (SE_SIDM_BABYL_KUGLER1, "Babylonian (Kugler 1)"),
    (SE_SIDM_BABYL_KUGLER2, "Babylonian (Kugler 2)"),
    (SE_SIDM_BABYL_KUGLER3, "Babylonian (Kugler 3)"),
    (SE_SIDM_BABYL_HUBER, "Babylonian (Huber)"),
    (SE_SIDM_BABYL_ETPSC, "Babylonian (ETPSC)"),
    (SE_SIDM_ALDEBARAN_15TAU, "Aldebaran at 15 Tau"),
    (SE_SIDM_HIPPARCHOS, "Hipparchos"),
    (SE_SIDM_SASSANIAN, "Sassanian"),
    (SE_SIDM_GALCENT_0SAG, "Galactic Center at 0 Sag"),
    (SE_SIDM_J2000, "J2000"),
    (SE_SIDM_J1900, "J1900"),
    (SE_SIDM_B1950, "B1950"),
    (SE_SIDM_SURYASIDDHANTA, "Suryasiddhanta"),
    (SE_SIDM_SURYASIDDHANTA_MSUN, "Suryasiddhanta (mean Sun)"),
    (SE_SIDM_ARYABHATA, "Aryabhata"),
    (SE_SIDM_ARYABHATA_MSUN, "Aryabhata (mean Sun)"),
    (SE_SIDM_SS_REVATI, "SS Revati"),
    (SE_SIDM_SS_CITRA, "SS Citra"),
    (SE_SIDM_TRUE_CITRA, "True Citra"),
    (SE_SIDM_TRUE_REVATI, "True Revati"),
    (SE_SIDM_TRUE_PUSHYA, "True Pushya"),
    (SE_SIDM_GALCENT_RGILBRAND, "Galactic Center (Gil Brand)"),
    (SE_SIDM_GALEQU_IAU1958, "Galactic Equator (IAU 1958)"),
    (SE_SIDM_GALEQU_TRUE, "Galactic Equator (True)"),
    (SE_SIDM_GALEQU_MULA, "Galactic Equator at Mula"),
    (SE_SIDM_GALALIGN_MARDYKS, "Galactic Alignment (Mardyks)"),
    (SE_SIDM_TRUE_MULA, "True Mula"),
    (SE_SIDM_GALCENT_MULA_WILHELM, "Galactic Center at Mula (Wilhelm)"),
    (SE_SIDM_ARYABHATA_522, "Aryabhata 522"),
    (SE_SIDM_BABYL_BRITTON, "Babylonian (Britton)"),
    (SE_SIDM_TRUE_SHEORAN, "True Sheoran"),
    (SE_SIDM_GALCENT_COCHRANE, "Galactic Center (Cochrane)"),
    (SE_SIDM_GALEQU_FIORENZA, "Galactic Equator (Fiorenza)"),
    (SE_SIDM_VALENS_MOON, "Valens (Moon)"),
]


@pytest.mark.parametrize("sid_mode, mode_name", AYANAMSHA_MODES)
def test_ayanamsha_modes(sid_mode, mode_name):
    """Test all ayanamsha modes for accuracy"""
    jd = swe.julday(2000, 1, 1, 12.0)

    # Swiss Ephemeris
    swe.set_sid_mode(sid_mode)
    aya_swe = swe.get_ayanamsa_ut(jd)

    # Python Ephemeris
    pyephem.swe_set_sid_mode(sid_mode)
    aya_py = pyephem.swe_get_ayanamsa_ut(jd)

    diff = abs(aya_swe - aya_py)
    print(
        f"\n{mode_name} ({sid_mode}): SWE={aya_swe:.6f}°, PY={aya_py:.6f}°, Diff={diff:.6f}°"
    )

    # Tolerance: ayanamsha should match within 0.1 degrees
    # Star-based modes might have slightly larger differences due to proper motion
    tol = (
        1.0
        if sid_mode
        in [
            SE_SIDM_TRUE_CITRA,
            SE_SIDM_TRUE_REVATI,
            SE_SIDM_TRUE_PUSHYA,
            SE_SIDM_TRUE_MULA,
            SE_SIDM_GALCENT_0SAG,
            SE_SIDM_GALCENT_RGILBRAND,
        ]
        else 0.1
    )

    assert diff < tol, f"Ayanamsha {mode_name} mismatch: {diff}° > {tol}°"


@pytest.mark.parametrize("name, year, month, day, hour", FAMOUS_PEOPLE_AYANAMSHA)
@pytest.mark.parametrize(
    "sid_mode, mode_name", AYANAMSHA_MODES[:10]
)  # Test first 10 modes for speed
def test_sidereal_planets(name, year, month, day, hour, sid_mode, mode_name):
    """Test sidereal planetary positions with various ayanamsha"""
    jd = swe.julday(year, month, day, hour)

    # Swiss Ephemeris
    swe.set_sid_mode(sid_mode)
    pos_swe_tuple, _ = swe.calc_ut(jd, SE_SUN, SEFLG_SWIEPH | SEFLG_SIDEREAL)

    # Python Ephemeris
    pyephem.swe_set_sid_mode(sid_mode)
    pos_py, _ = pyephem.swe_calc_ut(jd, SE_SUN, SEFLG_SIDEREAL)

    # Extract longitudes (handle both tuple and array)
    lon_swe = float(pos_swe_tuple[0])
    lon_py = float(pos_py[0])

    diff_lon = abs(lon_swe - lon_py)
    if diff_lon > 180:
        diff_lon = 360 - diff_lon

    print(
        f"\n{name} ({mode_name}): SWE Sun={lon_swe:.4f}°, PY Sun={lon_py:.4f}°, Diff={diff_lon:.4f}°"
    )

    # Planetary positions should match within reasonable tolerance
    # Note: Small differences expected due to ayanamsha rate changes over time
    assert diff_lon < 0.5, f"Sidereal Sun position mismatch for {name} with {mode_name}"


def test_tropical_vs_sidereal():
    """Test that tropical and sidereal modes produce correct difference"""
    jd = swe.julday(2000, 1, 1, 12.0)

    # Tropical - swe.calc_ut returns ((lon, lat, dist, speed_lon, speed_lat, speed_dist), flags)
    pos_trop_swe_tuple, _ = swe.calc_ut(jd, SE_SUN, SEFLG_SWIEPH)
    pos_trop_py, _ = pyephem.swe_calc_ut(jd, SE_SUN, 0)

    # Sidereal (Lahiri)
    swe.set_sid_mode(SE_SIDM_LAHIRI)
    pyephem.swe_set_sid_mode(SE_SIDM_LAHIRI)
    pos_sid_swe_tuple, _ = swe.calc_ut(jd, SE_SUN, SEFLG_SWIEPH | SEFLG_SIDEREAL)
    pos_sid_py, _ = pyephem.swe_calc_ut(jd, SE_SUN, SEFLG_SIDEREAL)

    # Get ayanamsha
    aya_swe = swe.get_ayanamsa_ut(jd)
    aya_py = pyephem.swe_get_ayanamsa_ut(jd)

    # Check that tropical - sidereal = ayanamsha
    # Extract longitude from results (index 0 of the data tuple)
    diff_swe = (pos_trop_swe_tuple[0] - pos_sid_swe_tuple[0]) % 360.0
    diff_py = (pos_trop_py[0] - pos_sid_py[0]) % 360.0

    print(
        f"\nTropical-Sidereal: SWE={diff_swe:.4f}° (aya={aya_swe:.4f}°), PY={diff_py:.4f}° (aya={aya_py:.4f}°)"
    )

    assert abs(diff_swe - aya_swe) < 0.01, (
        "SWE: Tropical - Sidereal should equal Ayanamsha"
    )
    assert abs(diff_py - aya_py) < 0.01, (
        "PY: Tropical - Sidereal should equal Ayanamsha"
    )


if __name__ == "__main__":
    # Run detailed comparison tests
    print("=" * 80)
    print(" " * 20 + "AYANAMSHA COMPARISON TEST")
    print("=" * 80)

    # Test multiple ayanamsha modes
    test_modes = [
        (SE_SIDM_LAHIRI, "Lahiri"),
        (SE_SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
        (SE_SIDM_KRISHNAMURTI, "Krishnamurti"),
        (SE_SIDM_RAMAN, "Raman"),
        (SE_SIDM_TRUE_CITRA, "True Citra"),
    ]

    jd = swe.julday(2000, 1, 1, 12.0)

    print(f"\nDate: 2000-01-01 12:00 UT (JD={jd})")
    print(
        "\n{:<25} {:>15} {:>15} {:>15}".format(
            "Mode", "SwissEph", "PythonEph", "Difference"
        )
    )
    print("-" * 70)

    for mode, name in test_modes:
        swe.set_sid_mode(mode)
        aya_swe = swe.get_ayanamsa_ut(jd)

        pyephem.swe_set_sid_mode(mode)
        aya_py = pyephem.swe_get_ayanamsa_ut(jd)

        diff = abs(aya_swe - aya_py)
        status = "✓" if diff < 0.1 else "✗"

        print(f"{name:<25} {aya_swe:>14.6f}° {aya_py:>14.6f}° {diff:>14.6f}° {status}")

    print("\n" + "=" * 80)
    print("TROPICAL vs SIDEREAL COMPARISON")
    print("=" * 80)

    # Tropical
    pos_trop_swe_tuple, _ = swe.calc_ut(jd, SE_SUN, SEFLG_SWIEPH)
    pos_trop_py, _ = pyephem.swe_calc_ut(jd, SE_SUN, 0)

    # Sidereal (Lahiri)
    swe.set_sid_mode(SE_SIDM_LAHIRI)
    pyephem.swe_set_sid_mode(SE_SIDM_LAHIRI)
    pos_sid_swe_tuple, _ = swe.calc_ut(jd, SE_SUN, SEFLG_SWIEPH | SEFLG_SIDEREAL)
    pos_sid_py, _ = pyephem.swe_calc_ut(jd, SE_SUN, SEFLG_SIDEREAL)

    aya_swe = swe.get_ayanamsa_ut(jd)
    aya_py = pyephem.swe_get_ayanamsa_ut(jd)

    diff_swe = (pos_trop_swe_tuple[0] - pos_sid_swe_tuple[0]) % 360.0
    diff_py = (pos_trop_py[0] - pos_sid_py[0]) % 360.0

    print(f"\nSun Tropical (SWE):  {pos_trop_swe_tuple[0]:.6f}°")
    print(f"Sun Tropical (PY):   {pos_trop_py[0]:.6f}°")
    print(f"Sun Sidereal (SWE):  {pos_sid_swe_tuple[0]:.6f}°")
    print(f"Sun Sidereal (PY):   {pos_sid_py[0]:.6f}°")
    print(f"\nAyanamsha (SWE):     {aya_swe:.6f}°")
    print(f"Ayanamsha (PY):      {aya_py:.6f}°")
    print(f"\nTrop-Sid (SWE):      {diff_swe:.6f}°")
    print(f"\nTrop-Sid (PY):       {diff_py:.6f}°")
    print(
        f"\nDifference match: {'✓ PASS' if abs(diff_swe - aya_swe) < 0.01 and abs(diff_py - aya_py) < 0.01 else '✗ FAIL'}"
    )

    # Comprehensive planetary comparison for all subjects
    print("\n" + "=" * 80)
    print("COMPREHENSIVE PLANETARY POSITIONS - ALL SUBJECTS")
    print("=" * 80)

    planets = [
        (SE_SUN, "Sun"),
        (SE_MOON, "Moon"),
        (SE_MERCURY, "Mercury"),
        (SE_VENUS, "Venus"),
        (SE_MARS, "Mars"),
        (SE_JUPITER, "Jupiter"),
        (SE_SATURN, "Saturn"),
        (SE_URANUS, "Uranus"),
        (SE_NEPTUNE, "Neptune"),
        (SE_PLUTO, "Pluto"),
    ]

    for person_name, year, month, day, hour in FAMOUS_PEOPLE_AYANAMSHA:
        jd_person = swe.julday(year, month, day, hour)

        print(f"\n{'=' * 80}")
        print(f"SUBJECT: {person_name} ({year}-{month:02d}-{day:02d} {hour:.2f}h)")
        print(f"{'=' * 80}")

        # TROPICAL positions
        print(f"\n{'TROPICAL POSITIONS':-^80}")
        print(
            f"\n{'Planet':<12} {'SwissEph':>15} {'PythonEph':>15} {'Difference':>15} {'Status':>8}"
        )
        print("-" * 77)

        for planet_id, planet_name in planets:
            try:
                pos_swe_tup, _ = swe.calc_ut(jd_person, planet_id, SEFLG_SWIEPH)
                pos_py_trop, _ = pyephem.swe_calc_ut(jd_person, planet_id, 0)

                lon_swe = float(pos_swe_tup[0])
                lon_py = float(pos_py_trop[0])

                diff = abs(lon_swe - lon_py)
                if diff > 180:
                    diff = 360 - diff

                status = "✓" if diff < 0.01 else ("~" if diff < 0.1 else "✗")
                print(
                    f"{planet_name:<12} {lon_swe:>14.6f}° {lon_py:>14.6f}° {diff:>14.6f}° {status:>8}"
                )
            except Exception as e:
                print(f"{planet_name:<12} {'ERROR':>15} {str(e)[:30]:>15}")

        # SIDEREAL positions (Lahiri)
        print(f"\n{'SIDEREAL POSITIONS (Lahiri)':-^80}")
        print(
            f"\n{'Planet':<12} {'SwissEph':>15} {'PythonEph':>15} {'Difference':>15} {'Status':>8}"
        )
        print("-" * 77)

        swe.set_sid_mode(SE_SIDM_LAHIRI)
        pyephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        for planet_id, planet_name in planets:
            try:
                pos_swe_tup, _ = swe.calc_ut(
                    jd_person, planet_id, SEFLG_SWIEPH | SEFLG_SIDEREAL
                )
                pos_py_sid, _ = pyephem.swe_calc_ut(
                    jd_person, planet_id, SEFLG_SIDEREAL
                )

                lon_swe = float(pos_swe_tup[0])
                lon_py = float(pos_py_sid[0])

                diff = abs(lon_swe - lon_py)
                if diff > 180:
                    diff = 360 - diff

                status = "✓" if diff < 0.01 else ("~" if diff < 0.1 else "✗")
                print(
                    f"{planet_name:<12} {lon_swe:>14.6f}° {lon_py:>14.6f}° {diff:>14.6f}° {status:>8}"
                )
            except Exception as e:
                print(f"{planet_name:<12} {'ERROR':>15} {str(e)[:30]:>15}")

    print("\n" + "=" * 80)
    print("Legend: ✓ = diff < 0.01°  ~ = diff < 0.1°  ✗ = diff >= 0.1°")
    print("=" * 80 + "\n")
