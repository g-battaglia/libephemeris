"""Round 144: Fixed Star Catalog Completeness.

Tests that all stars in our catalog are found by fixstar2_ut and
that positions match pyswisseph for the full catalog at multiple epochs.
"""

from __future__ import annotations
import sys, os
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import FLG_SPEED, FLG_SWIEPH

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)
ephem.set_ephe_path(_REF_EPHE_PATH)

# Major stars to test
STARS = [
    "Aldebaran",
    "Algol",
    "Altair",
    "Antares",
    "Arcturus",
    "Betelgeuse",
    "Canopus",
    "Capella",
    "Castor",
    "Deneb",
    "Fomalhaut",
    "Pollux",
    "Procyon",
    "Regulus",
    "Rigel",
    "Sirius",
    "Spica",
    "Vega",
    "Achernar",
    "Acrux",
    "Bellatrix",
    "Denebola",
    "Elnath",
    "Hamal",
    "Markab",
    "Menkar",
    "Mirach",
    "Mira",
    "Nunki",
    "Rasalhague",
    "Sadalmelik",
    "Scheat",
    "Shaula",
    "Thuban",
    "Unukalhai",
    "Wezen",
    "Zubenelgenubi",
    "Zubeneschamali",
    "Alphard",
    "Alphecca",
    "Alpheratz",
    "Alnilam",
    "Alnitak",
    "Mintaka",
    "Saiph",
    "Rigil Kentaurus",
    "Toliman",
    "Acubens",
    "Vindemiatrix",
    "Dschubba",
    "Graffias",
    "Kaus Australis",
    "Kaus Borealis",
    "Kaus Media",
    "Nashira",
    "Sadalsuud",
    "Algenib",
    "Almach",
    "Ankaa",
    "Dabih",
    "Diphda",
    "Enif",
    "Izar",
    "Kochab",
    "Menkent",
    "Mimosa",
    "Peacock",
    "Phact",
    "Phecda",
    "Sabik",
    "Sargas",
    "Suhail",
    "Tejat",
    "Tureis",
    "Wasat",
    "Yed Prior",
    "Zaniah",
    "Zosma",
    "Acrab",
    "Agena",
    "Al Hecka",
    "Albireo",
    "Alcyone",
    "Alderamin",
    "Algieba",
    "Alhena",
    "Alioth",
    "Alkaid",
    "Alkes",
    "Alnair",
    "Alphard",
    "Asellus Australis",
    "Avior",
    "Baten Kaitos",
    "Cor Caroli",
    "Dubhe",
    "Gacrux",
    "Gienah",
    "Megrez",
    "Menkib",
    "Merak",
    "Miaplacidus",
    "Muphrid",
    "Naos",
    "Nihal",
    "Polaris",
    "Rastaban",
    "Rotanev",
    "Rukbat",
    "Schedar",
    "Syrma",
    "Unukalhai",
    "Vindemiatrix",
]

# Remove duplicates
STARS = list(dict.fromkeys(STARS))

TEST_JDS = [2451545.0, 2455197.5, 2458849.5, 2460676.5, 2466154.5]

flags = FLG_SWIEPH | FLG_SPEED

passed = 0
failed = 0
errors = 0
not_found_se = 0
not_found_le = 0
total = 0

for jd in TEST_JDS:
    for star_name in STARS:
        total += 1
        try:
            se_result = swe.fixstar2_ut(star_name, jd, flags)
            se_lon = se_result[0][0]
            se_lat = se_result[0][1]
        except Exception as e:
            not_found_se += 1
            continue

        try:
            le_result = ephem.fixstar2_ut(star_name, jd, flags)
            le_lon = le_result[0][0]
            le_lat = le_result[0][1]
        except Exception as e:
            not_found_le += 1
            print(f"NOT_FOUND_LE {star_name}: {e}")
            failed += 1
            continue

        lon_diff = abs(le_lon - se_lon) * 3600.0
        lat_diff = abs(le_lat - se_lat) * 3600.0

        # 10" tolerance for position
        if lon_diff < 10.0 and lat_diff < 10.0:
            passed += 1
        else:
            failed += 1
            print(
                f'FAIL {star_name:20s} JD={jd:.1f} lon_diff={lon_diff:.3f}" lat_diff={lat_diff:.3f}"'
            )

print(f"\n{'=' * 60}")
print(f"Round 144: Fixed Star Catalog Completeness")
print(f"{'=' * 60}")
print(f"Total:   {total}")
print(f"Passed:  {passed} ({100 * passed / max(total, 1):.1f}%)")
print(f"Failed:  {failed}")
print(f"Not found in SE: {not_found_se}")
print(f"Not found in LE: {not_found_le}")
print(f"Errors:  {errors}")
