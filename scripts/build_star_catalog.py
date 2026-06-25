#!/usr/bin/env python3
"""
Build Star Catalog from Hipparcos using astroquery.vizier.

SUPERSEDED — DO NOT USE TO REGENERATE THE SHIPPED CATALOG.
    The shipped catalog (libephemeris/star_catalog_gen.py) is produced by
    scripts/build_star_catalog_v2.py from van Leeuwen 2007 (CDS I/311).
    This v1 builder contains the known-wrong HIP name mappings that caused
    REVIEW-2026-06-10 item C9 (HIP 5742/7097/42911 mislabeled "Eta
    Piscium"/"Alrescha"/"Tarf"); it is retained only for the historical
    record and its unit tests, not for regeneration. See CLEAN.md.

This script fetches astrologically relevant stars from the Hipparcos I/239/hip_main
catalog via VizieR and generates StarCatalogEntry data compatible with
fixstar/fixstar2 functions in libephemeris.

The selection criteria for ~800-1000 stars:
1. All stars brighter than magnitude 4.0 (~900 stars) - naked-eye visible
2. All named stars with traditional/proper names
3. Stars from zodiacal constellations
4. Stars used in ayanamsha calculations

Usage:
    pip install astroquery  # Required dependency
    python scripts/build_star_catalog.py                    # Report mode (default)
    python scripts/build_star_catalog.py --output catalog   # Generate Python code
    python scripts/build_star_catalog.py --output csv       # Generate CSV file
    python scripts/build_star_catalog.py --output json      # Generate JSON file
    python scripts/build_star_catalog.py --mag-limit 4.5    # Custom magnitude limit
    python scripts/build_star_catalog.py --verbose          # Show progress details

Output Format (StarCatalogEntry):
    - id: Internal star ID (FIXSTAR_OFFSET + index)
    - name: Traditional star name or Bayer designation
    - nomenclature: Abbreviated Bayer designation (e.g., "alLeo")
    - hip_number: Hipparcos catalog number
    - data: StarData with RA/Dec J2000, proper motion (arcsec/year)
    - magnitude: Visual magnitude (Vmag from Hipparcos)

Data Sources:
    - Hipparcos Catalog: I/239/hip_main (ESA 1997)
    - Star names: Cross-referenced with IAU star name catalog
    - Constellations: Hipparcos constellation assignments

References:
    - Hipparcos Catalog (ESA SP-1200): https://www.cosmos.esa.int/web/hipparcos
    - VizieR Catalog Access: https://vizier.u-strasbg.fr/
    - IAU Star Names: https://www.iau.org/public/themes/naming_stars/
"""

import argparse
import json
import os
import sys
from dataclasses import dataclass, asdict
from datetime import datetime

# Ensure libephemeris can be imported from the project root
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Try to import astroquery
try:
    from astroquery.vizier import Vizier
    from astropy import units as u

    ASTROQUERY_AVAILABLE = True
except ImportError:
    ASTROQUERY_AVAILABLE = False
    Vizier = None  # type: ignore
    u = None  # type: ignore

# Fixed star ID offset (from constants.py)
FIXSTAR_OFFSET = 1000000

# IAU official star names with Hipparcos numbers
# Source: https://www.iau.org/public/themes/naming_stars/
IAU_STAR_NAMES: dict[int, str] = {
    677: "Alpheratz",
    746: "Caph",
    1067: "Algenib",
    1562: "Ankaa",
    2081: "Schedar",
    3179: "Diphda",
    3419: "Achernar",
    4427: "Mirach",
    5447: "Acamar",
    5742: "Eta Piscium",
    7097: "Alrescha",
    7588: "Sheratan",
    8102: "Hamal",
    8903: "Polaris",
    9640: "Mira",
    9884: "Mesarthim",
    11767: "Polaris",
    13847: "Rana",
    14135: "Azha",
    14576: "Algol",
    15863: "Mirfak",
    17358: "Menkar",
    17378: "Botein",
    17702: "Atik",
    18532: "Electra",
    18543: "Celaeno",
    18614: "Taygeta",
    18724: "Maia",
    18871: "Asterope",
    19587: "Merope",
    20042: "Alcyone",
    20205: "Atlas",
    20455: "Pleione",
    20711: "Ain",
    20889: "Prima Hyadum",
    20894: "Secunda Hyadum",
    21393: "Theta Tauri",
    21421: "Aldebaran",
    22449: "Cursa",
    23015: "Rigel",
    23875: "Capella",
    24186: "Bellatrix",
    24436: "Rigel",
    24608: "Capella",
    25336: "Phact",
    25428: "Mintaka",
    25930: "Arneb",
    26207: "Nihal",
    26311: "Alnilam",
    26727: "Alnitak",
    27072: "Saiph",
    27366: "Meissa",
    27628: "Betelgeuse",
    27989: "Betelgeuse",
    28360: "Canopus",
    28380: "Alhena",
    28614: "Furud",
    30122: "Mirzam",
    30324: "Canopus",
    30438: "Canopus",
    31681: "Sirius",
    32246: "Adhara",
    32349: "Sirius",
    33579: "Wezen",
    33856: "Aludra",
    34444: "Castor",
    35264: "Gomeisa",
    36188: "Avior",
    36850: "Procyon",
    37279: "Procyon",
    37826: "Pollux",
    39429: "Naos",
    39757: "Regor",
    39953: "Tureis",
    41037: "Alsuhail",
    42913: "Suhail",
    44816: "Miaplacidus",
    45238: "Alphard",
    45556: "Regulus",
    45941: "Merak",
    46390: "Dubhe",
    46733: "Algieba",
    49669: "Regulus",
    50099: "Adhafera",
    50335: "Tania Borealis",
    50372: "Alula Australis",
    50583: "Alula Borealis",
    51069: "Chertan",
    53910: "Denebola",
    54061: "Phecda",
    54872: "Alkes",
    55219: "Zosma",
    57632: "Megrez",
    57757: "Giausar",
    58001: "Alioth",
    59196: "Cor Caroli",
    59747: "Vindemiatrix",
    59774: "Mizar",
    60718: "Spica",
    61084: "Alcor",
    61281: "Algorab",
    61359: "Kraz",
    62434: "Gienah",
    63090: "Minelauva",
    63608: "Alkaid",
    65109: "Heze",
    65378: "Mizar",
    65474: "Spica",
    66657: "Syrma",
    67301: "Kang",
    67459: "Menkent",
    68002: "Arcturus",
    68702: "Zubenelgenubi",
    69673: "Arcturus",
    71683: "Rigil Kentaurus",
    71795: "Zubeneschamali",
    72105: "Unukalhai",
    72607: "Zubenelhakrabi",
    72622: "Yed Prior",
    73184: "Toliman",
    74785: "Yed Posterior",
    75097: "Nekkar",
    75141: "Ed Asich",
    76267: "Dschubba",
    76952: "Acrab",
    77070: "Graffias",
    78401: "Sabik",
    78820: "Rasalgethi",
    79593: "Sargas",
    80112: "Lesath",
    80331: "Shaula",
    80763: "Antares",
    81266: "Rastaban",
    81377: "Atria",
    82273: "Kornephoros",
    82396: "Marfik",
    83000: "Kaus Media",
    84012: "Kaus Australis",
    84143: "Kaus Borealis",
    84345: "Rasalhague",
    85670: "Cebalrai",
    85693: "Shaula",
    85927: "Nunki",
    86032: "Etamin",
    86228: "Albaldah",
    86742: "Alnasl",
    87073: "Ascella",
    87833: "Vega",
    88635: "Sheliak",
    89931: "Nunki",
    90185: "Sulafat",
    90496: "Arkab Posterior",
    91262: "Vega",
    92041: "Albireo",
    92855: "Algedi",
    93506: "Dabih",
    93864: "Kochab",
    95168: "Sadr",
    95947: "Altair",
    96295: "Nashira",
    97278: "Tarazed",
    97649: "Altair",
    98036: "Deneb Algedi",
    98337: "Giausar",
    100453: "Sadalsuud",
    100751: "Deneb",
    102098: "Deneb",
    102488: "Sadalmelik",
    105199: "Albali",
    107315: "Enif",
    107556: "Skat",
    109074: "Sadalmelik",
    109268: "Homam",
    109427: "Matar",
    112029: "Sadalbari",
    112158: "Biham",
    112440: "Ancha",
    112748: "Fomalhaut",
    113136: "Skat",
    113368: "Fomalhaut",
    113881: "Scheat",
    113963: "Markab",
}

# Additional traditional star names (not in IAU official list)
TRADITIONAL_STAR_NAMES: dict[int, str] = {
    # Ursa Major (Big Dipper)
    53910: "Denebola",
    # Orion's Belt
    25336: "Phact",
    # Additional Zodiacal stars
    20711: "Ain",
    20889: "Prima Hyadum",
    20894: "Secunda Hyadum",
    # Cancer constellation
    42806: "Acubens",
    42911: "Tarf",
    42516: "Asellus Borealis",
    43103: "Asellus Australis",
    # Sagittarius constellation
    87073: "Ascella",
    86742: "Alnasl",
    89931: "Nunki",
    # Capricornus constellation
    92855: "Algedi",
    93506: "Dabih",
    96295: "Nashira",
    # Aquarius constellation
    100453: "Sadalsuud",
    109074: "Sadalmelik",
    113136: "Skat",
    # Scorpius constellation
    79593: "Sargas",
    80112: "Lesath",
    80331: "Shaula",
    76267: "Dschubba",
    76952: "Acrab",
    77070: "Graffias",
    # Centaurus constellation
    68933: "Muhlifain",
    71352: "Epsilon Centauri",
    71865: "Eta Centauri",
    68282: "Zeta Centauri",
}

# Bayer designation Greek letters mapping
BAYER_GREEK_MAP = {
    "alf": "al",
    "bet": "be",
    "gam": "ga",
    "del": "de",
    "eps": "ep",
    "zet": "ze",
    "eta": "et",
    "the": "th",
    "iot": "io",
    "kap": "ka",
    "lam": "la",
    "mu.": "mu",
    "nu.": "nu",
    "ksi": "xi",
    "omi": "om",
    "pi.": "pi",
    "rho": "rh",
    "sig": "si",
    "tau": "ta",
    "ups": "up",
    "phi": "ph",
    "chi": "ch",
    "psi": "ps",
    "ome": "om",
}

# Constellation abbreviations (3-letter IAU codes)
CONSTELLATION_ABBREV = {
    "And": "And",
    "Ant": "Ant",
    "Aps": "Aps",
    "Aqr": "Aqr",
    "Aql": "Aql",
    "Ara": "Ara",
    "Ari": "Ari",
    "Aur": "Aur",
    "Boo": "Boo",
    "Cae": "Cae",
    "Cam": "Cam",
    "Cnc": "Cnc",
    "CVn": "CVn",
    "CMa": "CMa",
    "CMi": "CMi",
    "Cap": "Cap",
    "Car": "Car",
    "Cas": "Cas",
    "Cen": "Cen",
    "Cep": "Cep",
    "Cet": "Cet",
    "Cha": "Cha",
    "Cir": "Cir",
    "Col": "Col",
    "Com": "Com",
    "CrA": "CrA",
    "CrB": "CrB",
    "Crv": "Crv",
    "Crt": "Crt",
    "Cru": "Cru",
    "Cyg": "Cyg",
    "Del": "Del",
    "Dor": "Dor",
    "Dra": "Dra",
    "Equ": "Equ",
    "Eri": "Eri",
    "For": "For",
    "Gem": "Gem",
    "Gru": "Gru",
    "Her": "Her",
    "Hor": "Hor",
    "Hya": "Hya",
    "Hyi": "Hyi",
    "Ind": "Ind",
    "Lac": "Lac",
    "Leo": "Leo",
    "LMi": "LMi",
    "Lep": "Lep",
    "Lib": "Lib",
    "Lup": "Lup",
    "Lyn": "Lyn",
    "Lyr": "Lyr",
    "Men": "Men",
    "Mic": "Mic",
    "Mon": "Mon",
    "Mus": "Mus",
    "Nor": "Nor",
    "Oct": "Oct",
    "Oph": "Oph",
    "Ori": "Ori",
    "Pav": "Pav",
    "Peg": "Peg",
    "Per": "Per",
    "Phe": "Phe",
    "Pic": "Pic",
    "Psc": "Psc",
    "PsA": "PsA",
    "Pup": "Pup",
    "Pyx": "Pyx",
    "Ret": "Ret",
    "Sge": "Sge",
    "Sgr": "Sgr",
    "Sco": "Sco",
    "Scl": "Scl",
    "Sct": "Sct",
    "Ser": "Ser",
    "Sex": "Sex",
    "Tau": "Tau",
    "Tel": "Tel",
    "Tri": "Tri",
    "TrA": "TrA",
    "Tuc": "Tuc",
    "UMa": "UMa",
    "UMi": "UMi",
    "Vel": "Vel",
    "Vir": "Vir",
    "Vol": "Vol",
    "Vul": "Vul",
}


@dataclass
class StarCatalogData:
    """Data structure for a star catalog entry."""

    hip_number: int
    name: str
    nomenclature: str
    ra_j2000: float  # Right Ascension in degrees
    dec_j2000: float  # Declination in degrees
    pm_ra: float  # Proper motion RA (arcsec/year), includes cos(dec)
    pm_dec: float  # Proper motion Dec (arcsec/year)
    magnitude: float  # Visual magnitude


def check_astroquery() -> bool:
    """Check if astroquery library is available."""
    return ASTROQUERY_AVAILABLE


def fetch_hipparcos_catalog(
    mag_limit: float = 4.0, verbose: bool = False
) -> "list[StarCatalogData]":
    """
    Fetch stars from Hipparcos I/239/hip_main catalog via VizieR.

    Args:
        mag_limit: Maximum visual magnitude (default 4.0 for ~900 stars)
        verbose: Print progress information

    Returns:
        List of StarCatalogData objects
    """
    if not ASTROQUERY_AVAILABLE:
        print("Error: astroquery is not installed.", file=sys.stderr)
        print("Install with: pip install astroquery", file=sys.stderr)
        return []

    if verbose:
        print(f"Fetching Hipparcos catalog (Vmag <= {mag_limit})...")

    # Configure Vizier to get all columns we need
    v = Vizier(
        columns=[
            "HIP",  # Hipparcos number
            "RAJ2000",  # Right Ascension J2000 (degrees)
            "DEJ2000",  # Declination J2000 (degrees)
            "pmRA",  # Proper motion RA (mas/yr), includes cos(dec)
            "pmDE",  # Proper motion Dec (mas/yr)
            "Vmag",  # Visual magnitude
            "B-V",  # B-V color index
            "Plx",  # Parallax (mas)
            "SpType",  # Spectral type
        ],
        row_limit=-1,  # No limit
    )

    # Query Hipparcos catalog with magnitude constraint
    try:
        result = v.query_constraints(catalog="I/239/hip_main", Vmag=f"<={mag_limit}")
    except Exception as e:
        print(f"Error querying VizieR: {e}", file=sys.stderr)
        return []

    if not result or len(result) == 0:
        print("No results from VizieR query", file=sys.stderr)
        return []

    table = result[0]
    if verbose:
        print(f"  Retrieved {len(table)} stars from Hipparcos catalog")

    stars: list[StarCatalogData] = []

    for row in table:
        hip = int(row["HIP"])

        # Get position (required)
        ra = float(row["RAJ2000"])  # degrees
        dec = float(row["DEJ2000"])  # degrees

        # Get proper motion (convert mas/yr to arcsec/yr)
        # Handle masked/missing values
        try:
            pm_ra = float(row["pmRA"]) / 1000.0  # mas/yr -> arcsec/yr
        except (ValueError, TypeError):
            pm_ra = 0.0

        try:
            pm_dec = float(row["pmDE"]) / 1000.0  # mas/yr -> arcsec/yr
        except (ValueError, TypeError):
            pm_dec = 0.0

        # Get magnitude
        try:
            vmag = float(row["Vmag"])
        except (ValueError, TypeError):
            vmag = 99.0  # Unknown magnitude

        # Determine star name
        name = get_star_name(hip)

        # Generate nomenclature (Bayer designation abbreviation)
        nomenclature = generate_nomenclature(hip, name)

        stars.append(
            StarCatalogData(
                hip_number=hip,
                name=name,
                nomenclature=nomenclature,
                ra_j2000=ra,
                dec_j2000=dec,
                pm_ra=pm_ra,
                pm_dec=pm_dec,
                magnitude=vmag,
            )
        )

    if verbose:
        named_count = sum(1 for s in stars if not s.name.startswith("HIP "))
        print(f"  {named_count} stars have traditional names")

    return stars


def get_star_name(hip: int) -> str:
    """
    Get the traditional name for a star by Hipparcos number.

    Args:
        hip: Hipparcos catalog number

    Returns:
        Traditional name if known, otherwise "HIP {number}"
    """
    # Check IAU official names first
    if hip in IAU_STAR_NAMES:
        return IAU_STAR_NAMES[hip]

    # Check additional traditional names
    if hip in TRADITIONAL_STAR_NAMES:
        return TRADITIONAL_STAR_NAMES[hip]

    # Default to Hipparcos number
    return f"HIP {hip}"


def generate_nomenclature(hip: int, name: str) -> str:
    """
    Generate a short nomenclature string for a star.

    For named stars, this tries to generate a Bayer-style designation.
    For unnamed stars, uses the Hipparcos number.

    Args:
        hip: Hipparcos catalog number
        name: Star name

    Returns:
        Short nomenclature string (e.g., "alLeo", "bePer", "HIP32349")
    """
    # For now, use the HIP number as nomenclature for simplicity
    # A more complete implementation would cross-reference with Bayer/Flamsteed catalogs
    return f"HIP{hip}"


def output_python_catalog(stars: list[StarCatalogData], verbose: bool = False) -> str:
    """
    Generate Python code for StarCatalogEntry list.

    Args:
        stars: List of star data
        verbose: Print progress

    Returns:
        Python code string
    """
    lines = [
        '"""',
        "Auto-generated star catalog from Hipparcos I/239/hip_main.",
        f"Generated: {datetime.now().isoformat()}",
        f"Star count: {len(stars)}",
        '"""',
        "",
        "from libephemeris.fixed_stars import StarData, StarCatalogEntry",
        "from libephemeris.constants import FIXSTAR_OFFSET",
        "",
        "# Auto-generated star catalog entries",
        "HIPPARCOS_STAR_CATALOG: list[StarCatalogEntry] = [",
    ]

    for i, star in enumerate(sorted(stars, key=lambda s: s.hip_number)):
        star_id = FIXSTAR_OFFSET + 1000 + i  # Use offset to avoid conflicts
        lines.append("    StarCatalogEntry(")
        lines.append(f"        id={star_id},")
        lines.append(f'        name="{star.name}",')
        lines.append(f'        nomenclature="{star.nomenclature}",')
        lines.append(f"        hip_number={star.hip_number},")
        lines.append("        data=StarData(")
        lines.append(f"            ra_j2000={star.ra_j2000:.6f},")
        lines.append(f"            dec_j2000={star.dec_j2000:.6f},")
        lines.append(f"            pm_ra={star.pm_ra:.5f},")
        lines.append(f"            pm_dec={star.pm_dec:.5f},")
        lines.append("        ),")
        lines.append(f"        magnitude={star.magnitude:.2f},")
        lines.append("    ),")

    lines.append("]")
    lines.append("")

    return "\n".join(lines)


def output_csv(stars: list[StarCatalogData]) -> str:
    """
    Generate CSV output for star catalog.

    Args:
        stars: List of star data

    Returns:
        CSV string
    """
    lines = ["hip_number,name,nomenclature,ra_j2000,dec_j2000,pm_ra,pm_dec,magnitude"]

    for star in sorted(stars, key=lambda s: s.hip_number):
        # Escape commas in names
        name = star.name.replace('"', '""')
        if "," in name:
            name = f'"{name}"'

        lines.append(
            f"{star.hip_number},{name},{star.nomenclature},"
            f"{star.ra_j2000:.6f},{star.dec_j2000:.6f},"
            f"{star.pm_ra:.5f},{star.pm_dec:.5f},{star.magnitude:.2f}"
        )

    return "\n".join(lines)


def output_json(stars: list[StarCatalogData]) -> str:
    """
    Generate JSON output for star catalog.

    Args:
        stars: List of star data

    Returns:
        JSON string
    """
    data = {
        "catalog": "Hipparcos I/239/hip_main",
        "generated": datetime.now().isoformat(),
        "count": len(stars),
        "stars": [asdict(star) for star in sorted(stars, key=lambda s: s.hip_number)],
    }
    return json.dumps(data, indent=2)


def print_report(stars: list[StarCatalogData]) -> None:
    """
    Print a summary report of the catalog.

    Args:
        stars: List of star data
    """
    print("\n" + "=" * 60)
    print("HIPPARCOS STAR CATALOG SUMMARY")
    print("=" * 60)
    print(f"Total stars: {len(stars)}")

    named_stars = [s for s in stars if not s.name.startswith("HIP ")]
    print(f"Named stars: {len(named_stars)}")

    # Magnitude distribution
    mag_bins = {"<0": 0, "0-1": 0, "1-2": 0, "2-3": 0, "3-4": 0, ">4": 0}
    for star in stars:
        if star.magnitude < 0:
            mag_bins["<0"] += 1
        elif star.magnitude < 1:
            mag_bins["0-1"] += 1
        elif star.magnitude < 2:
            mag_bins["1-2"] += 1
        elif star.magnitude < 3:
            mag_bins["2-3"] += 1
        elif star.magnitude < 4:
            mag_bins["3-4"] += 1
        else:
            mag_bins[">4"] += 1

    print("\nMagnitude distribution:")
    for mag_range, count in mag_bins.items():
        print(f"  Vmag {mag_range}: {count} stars")

    # Show brightest stars
    print("\nBrightest stars (Vmag < 1.0):")
    brightest = sorted(
        [s for s in stars if s.magnitude < 1.0], key=lambda s: s.magnitude
    )
    for star in brightest[:10]:
        print(f"  {star.name:20s} HIP {star.hip_number:6d}  Vmag={star.magnitude:+.2f}")

    # Show some named stars
    print("\nSample of named stars:")
    named_sample = sorted(named_stars, key=lambda s: s.magnitude)[:15]
    for star in named_sample:
        print(
            f"  {star.name:20s} HIP {star.hip_number:6d}  "
            f"RA={star.ra_j2000:8.4f} Dec={star.dec_j2000:+8.4f}  Vmag={star.magnitude:.2f}"
        )

    print("\n" + "=" * 60)


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Build star catalog from Hipparcos via VizieR",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                        # Report mode (default)
  %(prog)s --output catalog       # Generate Python code
  %(prog)s --output csv           # Generate CSV file
  %(prog)s --output json          # Generate JSON file
  %(prog)s --mag-limit 4.5        # Include fainter stars
        """,
    )

    parser.add_argument(
        "--output",
        "-o",
        choices=["report", "catalog", "csv", "json"],
        default="report",
        help="Output format (default: report)",
    )

    parser.add_argument(
        "--mag-limit",
        "-m",
        type=float,
        default=4.0,
        help="Maximum visual magnitude (default: 4.0 for ~900 stars)",
    )

    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Print verbose output",
    )

    parser.add_argument(
        "--file",
        "-f",
        type=str,
        help="Output file path (default: stdout)",
    )

    args = parser.parse_args()

    # Check for astroquery
    if not check_astroquery():
        print("Error: astroquery is not installed.", file=sys.stderr)
        print("Install with: pip install astroquery", file=sys.stderr)
        return 1

    # Fetch catalog
    stars = fetch_hipparcos_catalog(mag_limit=args.mag_limit, verbose=args.verbose)

    if not stars:
        print("Error: No stars retrieved from catalog", file=sys.stderr)
        return 1

    # Generate output
    if args.output == "report":
        print_report(stars)
    elif args.output == "catalog":
        output = output_python_catalog(stars, verbose=args.verbose)
        if args.file:
            with open(args.file, "w") as f:
                f.write(output)
            print(f"Wrote Python catalog to {args.file}")
        else:
            print(output)
    elif args.output == "csv":
        output = output_csv(stars)
        if args.file:
            with open(args.file, "w") as f:
                f.write(output)
            print(f"Wrote CSV to {args.file}")
        else:
            print(output)
    elif args.output == "json":
        output = output_json(stars)
        if args.file:
            with open(args.file, "w") as f:
                f.write(output)
            print(f"Wrote JSON to {args.file}")
        else:
            print(output)

    return 0


if __name__ == "__main__":
    sys.exit(main())
