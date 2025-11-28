import swisseph as swe
import libephemeris as pyephem

# Constants
SE_SUN = 0
SE_MOON = 1
SE_PLUTO = 9
SEFLG_SWIEPH = 2
SEFLG_SIDEREAL = 64 * 1024
SEFLG_SPEED = 256

# House Systems to test
HOUSE_SYSTEMS = {
    "P": "Placidus",
    "K": "Koch",
    "R": "Regiomontanus",
    "C": "Campanus",
    "E": "Equal",
    "W": "Whole Sign",
    "M": "Morinus",
    "B": "Alcabitius",
    "X": "Meridian",
    "V": "Vehlow",
    "T": "Polich/Page",
    # 'F': 'Carter', # Carter might have specific definitions
    # 'U': 'Krusinski'
}

# Ayanamsha modes to test
AYANAMSHA_MODES = [
    (0, "Fagan/Bradley"),
    (1, "Lahiri"),
    (3, "Raman"),
    (5, "Krishnamurti"),
    (27, "True Citra"),
]

# Test Subjects: (Name, Year, Month, Day, Hour, Lat, Lon)
# Covering different epochs and locations (including high latitudes)
SUBJECTS = [
    ("Standard J2000", 2000, 1, 1, 12.0, 0.0, 0.0),
    ("Rome", 1980, 5, 20, 14.5, 41.9028, 12.4964),
    ("New York", 2024, 11, 5, 9.0, 40.7128, -74.0060),
    ("Sydney", 1950, 10, 15, 22.0, -33.8688, 151.2093),
    ("High Lat North (Tromso)", 1990, 1, 15, 12.0, 69.6492, 18.9553),
    ("High Lat South (McMurdo)", 2005, 6, 21, 0.0, -77.8463, 166.6681),
    ("Equator", 1975, 3, 21, 12.0, 0.0, 45.0),
    # ("Historical (1500)", 1500, 1, 1, 12.0, 51.5074, -0.1278), # DE421 starts 1900
]


def get_diff(val1, val2):
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


def compare_planets(jd, name, date_str, lat, lon):
    # print(f"\n--- Planetary Positions (Tropical) for {name} ---")
    # print(f"{'Planet':<10} {'SWE':<12} {'PY':<12} {'Diff':<12} {'Status'}")
    # print("-" * 60)

    max_diff = 0.0

    for planet in range(SE_SUN, SE_PLUTO + 1):
        # SWE
        res_swe, _ = swe.calc_ut(jd, planet, SEFLG_SWIEPH | SEFLG_SPEED)
        # PY
        res_py, _ = pyephem.swe_calc_ut(jd, planet, SEFLG_SPEED)

        lon_swe = res_swe[0]
        lon_py = res_py[0]

        diff = get_diff(lon_swe, lon_py)
        max_diff = max(max_diff, diff)

        status = "✓" if diff < 0.001 else "✗"
        # Context: [Subject] [Date] [Type] [Item]
        print(
            f"[{name}] [{date_str}] [Tropical] Planet {planet:<2}: SWE={lon_swe:10.6f} PY={lon_py:10.6f} Diff={diff:10.8f} {status}"
        )

    return max_diff


def compare_sidereal(jd, name, date_str):
    # print(f"\n--- Sidereal Positions for {name} ---")
    # print(f"{'Mode':<20} {'Planet':<10} {'SWE':<12} {'PY':<12} {'Diff':<12} {'Status'}")
    # print("-" * 80)

    max_diff = 0.0

    for mode_id, mode_name in AYANAMSHA_MODES:
        swe.set_sid_mode(mode_id)
        pyephem.swe_set_sid_mode(mode_id)

        # Test Sun and Moon
        for planet in [SE_SUN, SE_MOON]:
            res_swe, _ = swe.calc_ut(jd, planet, SEFLG_SWIEPH | SEFLG_SIDEREAL)
            res_py, _ = pyephem.swe_calc_ut(jd, planet, SEFLG_SIDEREAL)

            lon_swe = res_swe[0]
            lon_py = res_py[0]

            diff = get_diff(lon_swe, lon_py)
            max_diff = max(max_diff, diff)

            status = "✓" if diff < 0.001 else "✗"
            print(
                f"[{name}] [{date_str}] [Sidereal] Mode {mode_name:<15} Planet {planet:<2}: SWE={lon_swe:10.6f} PY={lon_py:10.6f} Diff={diff:10.8f} {status}"
            )

    return max_diff


def compare_houses(jd, name, date_str, lat, lon):
    # print(f"\n--- House Cusps for {name} (Lat: {lat}, Lon: {lon}) ---")
    # print(f"{'System':<15} {'Cusp 1':<12} {'Cusp 10':<12} {'Max Diff':<12} {'Status'}")
    # print("-" * 70)

    max_diff_total = 0.0

    for hsys_char, hsys_name in HOUSE_SYSTEMS.items():
        try:
            # SWE
            cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, bytes(hsys_char, "utf-8"))
            # PY
            cusps_py, ascmc_py = pyephem.swe_houses(
                jd, lat, lon, bytes(hsys_char, "utf-8")
            )

            # Handle different lengths
            # pyswisseph might return 12 (cusps 1-12) or 13 (0 + cusps 1-12)

            val_swe_1 = 0.0
            val_swe_10 = 0.0
            current_max_diff = 0.0

            for i in range(1, 13):
                val_swe = 0.0
                if len(cusps_swe) == 12:
                    val_swe = cusps_swe[i - 1]
                elif len(cusps_swe) == 13:
                    val_swe = cusps_swe[i]
                else:
                    print(
                        f"[{name}] [{date_str}] [Houses  ] System {hsys_name:<15}: ERROR Unexpected length {len(cusps_swe)}"
                    )
                    break

                val_py = cusps_py[i]  # pyephem always returns 13

                if i == 1:
                    val_swe_1 = val_swe
                if i == 10:
                    val_swe_10 = val_swe

                diff = get_diff(val_swe, val_py)
                current_max_diff = max(current_max_diff, diff)

            # Compare Asc/MC
            # ascmc indices: 0=Asc, 1=MC, 2=ARMC, 3=Vertex
            labels = ["Asc", "MC", "ARMC", "Vertex"]
            for i in range(4):
                if i < len(ascmc_swe) and i < len(ascmc_py):
                    diff = get_diff(ascmc_swe[i], ascmc_py[i])
                    current_max_diff = max(current_max_diff, diff)
                    if diff > 0.01:
                        print(
                            f"[{name}] [{date_str}] [Houses  ] System {hsys_name:<15} Angle {labels[i]:<6}: SWE={ascmc_swe[i]:.4f} PY={ascmc_py[i]:.4f} Diff={diff:.4f} ✗"
                        )

            max_diff_total = max(max_diff_total, current_max_diff)
            status = "✓" if current_max_diff < 0.01 else "✗"

            print(
                f"[{name}] [{date_str}] [Houses  ] System {hsys_name:<15}: Cusp1={val_swe_1:10.4f} Cusp10={val_swe_10:10.4f} MaxDiff={current_max_diff:10.6f} {status}"
            )

        except Exception as e:
            # Placidus and Koch are expected to fail at high latitudes
            if (
                "lat" in str(e).lower()
                or "polar" in str(e).lower()
                or "error" in str(e).lower()
            ):
                if hsys_name in ["Placidus", "Koch"] and abs(lat) > 66:
                    print(
                        f"[{name}] [{date_str}] [Houses  ] System {hsys_name:<15}: SKIPPED (Expected high lat error)"
                    )
                    continue

            print(f"[{name}] [{date_str}] [Houses  ] System {hsys_name:<15}: ERROR {e}")

    return max_diff_total


def main():
    print("================================================================")
    print("PYTHON EPHEMERIS vs SWISS EPHEMERIS - COMPREHENSIVE COMPARISON")
    print("================================================================")

    total_errors = 0

    for name, year, month, day, hour, lat, lon in SUBJECTS:
        jd = swe.julday(year, month, day, hour)
        date_str = f"{year}-{month}-{day}"
        # print(f"\n\n>>> SUBJECT: {name}")
        # print(f"    Date: {year}-{month}-{day} {hour}h UT")
        # print(f"    Loc : {lat}, {lon}")
        # print("=" * 60)

        # 1. Planets Tropical
        diff_trop = compare_planets(jd, name, date_str, lat, lon)
        if diff_trop > 0.001:
            total_errors += 1

        # 2. Sidereal
        # Relax tolerance for star-based ayanamshas (catalog differences)
        diff_sid = compare_sidereal(jd, name, date_str)
        if diff_sid > 0.01:
            total_errors += 1

        # 3. Houses
        diff_houses = compare_houses(jd, name, date_str, lat, lon)
        if diff_houses > 0.1:
            total_errors += 1

    print("\n\n" + "=" * 60)
    print("FINAL SUMMARY")
    print("=" * 60)
    if total_errors == 0:
        print("ALL CHECKS PASSED! ✓")
    else:
        print(f"FOUND {total_errors} ISSUES.")


if __name__ == "__main__":
    main()
