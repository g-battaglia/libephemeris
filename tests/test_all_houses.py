import pytest
import swisseph as swe
import libephemeris as pyephem
from libephemeris.constants import *

# Famous People Data (Name, Year, Month, Day, Hour, Lat, Lon)
FAMOUS_PEOPLE = [
    ("Albert Einstein", 1879, 3, 14, 11.5, 48.4011, 10.0024),  # Ulm, Germany
    ("Marilyn Monroe", 1926, 6, 1, 9.5, 34.0522, -118.2437),  # Los Angeles
    ("Mahatma Gandhi", 1869, 10, 2, 7.75, 21.1702, 72.8311),  # Porbandar, India
    ("Pablo Picasso", 1881, 10, 25, 23.5, 36.7213, -4.4214),  # Málaga, Spain
    ("Nelson Mandela", 1918, 7, 18, 14.0, -31.9505, 28.7323),  # Mvezo, South Africa
    ("Marie Curie", 1867, 11, 7, 12.0, 52.2297, 21.0122),  # Warsaw, Poland
    ("Frida Kahlo", 1907, 7, 6, 8.5, 19.4326, -99.1332),  # Coyoacán, Mexico
    ("Leonardo da Vinci", 1452, 4, 15, 21.5, 43.7833, 10.9167),  # Vinci, Italy
    ("Bob Marley", 1945, 2, 6, 2.5, 18.0179, -76.8099),  # Nine Mile, Jamaica
    ("Equator Test", 2000, 1, 1, 12.0, 0.0, 0.0),  # Equator
]

HOUSE_SYSTEMS = [
    "P",
    "K",
    "R",
    "C",
    "E",
    "W",
    "O",  # Implemented initially
    "B",
    "T",
    "M",
    "X",
    "V",
    "F",  # Newly implemented
    # 'U', 'G' # Krusinski, Gauquelin (Placeholders/Complex)
]


@pytest.mark.parametrize("name, year, month, day, hour, lat, lon", FAMOUS_PEOPLE)
@pytest.mark.parametrize("hsys", HOUSE_SYSTEMS)
def test_houses_famous_people(name, year, month, day, hour, lat, lon, hsys):
    jd = swe.julday(year, month, day, hour)

    # Swiss Ephemeris
    try:
        cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, hsys.encode())
    except swe.Error as e:
        pytest.skip(f"SwissEph Error for {hsys} at {lat}: {e}")

    # Python Ephemeris
    try:
        cusps_py, ascmc_py = pyephem.swe_houses(jd, lat, lon, ord(hsys))
    except Exception as e:
        pytest.fail(f"PythonEphemeris Error for {hsys} at {lat}: {e}")

    # Compare Cusps
    print(f"\nTesting {name} ({lat}, {lon}) with System {hsys}")
    for i in range(1, 13):
        # pyswisseph: cusps[i-1] = house i (0-indexed)
        # libephemeris: cusps[i-1] = house i (0-indexed, now matching pyswisseph)
        diff = abs(cusps_swe[i - 1] - cusps_py[i - 1])
        if diff > 180:
            diff = 360 - diff

        # Tolerance varies by system complexity
        tol = 0.5
        if hsys in ["P", "K"]:
            tol = 0.1  # Iterative
        if hsys in ["T"]:
            tol = 1.0  # Topocentric might need tuning
        if hsys in ["F"]:
            tol = 1.0  # Carter might have small differences

        assert diff < tol, (
            f"House {i} mismatch: SWE={cusps_swe[i - 1]}, PY={cusps_py[i - 1]}, Diff={diff}"
        )

    # Compare Asc/MC
    # Note: Asc/MC should be same for most quadrant systems, but might differ for some (e.g. Vehlow, Morinus?)
    # swe_houses returns standard Asc/MC in ascmc array regardless of system?
    # Usually yes.

    diff_asc = abs(ascmc_swe[0] - ascmc_py[0])
    if diff_asc > 180:
        diff_asc = 360 - diff_asc
    assert diff_asc < 0.1, f"Asc mismatch: SWE={ascmc_swe[0]}, PY={ascmc_py[0]}"

    diff_mc = abs(ascmc_swe[1] - ascmc_py[1])
    if diff_mc > 180:
        diff_mc = 360 - diff_mc
    assert diff_mc < 0.1, f"MC mismatch: SWE={ascmc_swe[1]}, PY={ascmc_py[1]}"


if __name__ == "__main__":
    # Manual run
    for person in FAMOUS_PEOPLE:
        for hsys in HOUSE_SYSTEMS:
            try:
                test_houses_famous_people(*person, hsys)
                print(f"PASS: {person[0]} - {hsys}")
            except Exception as e:
                print(f"FAIL: {person[0]} - {hsys} -> {e}")
