"""
Fixed star position calculations for libephemeris.

Computes ecliptic positions for bright fixed stars with:
- Proper motion correction (rigorous space motion approach)
- IAU 2006 precession from J2000 to date
- IAU 2000A nutation model (1365 terms) for sub-milliarcsecond precision
- Equatorial to ecliptic coordinate transformation

Supported stars:
- The four Royal Stars of Persia (watchers of the four directions):
  - Aldebaran (Alpha Tauri) - Watcher of the East (Tascheter)
  - Regulus (Alpha Leonis) - Watcher of the North (Venant)
  - Antares (Alpha Scorpii) - Watcher of the West (Satevis)
  - Fomalhaut (Alpha Piscis Austrini) - Watcher of the South (Hastorang)
- Spica (Alpha Virginis) - Used for ayanamsha calculations
- And many other bright stars...

Notes on precision:
    Proper motion is applied using the rigorous space motion approach from
    Hipparcos Vol. 1, Section 1.5.5 with second-order Taylor expansion.
    This method converts the position to a 3D unit vector, applies proper
    motion as angular velocity in the tangent plane with centripetal
    acceleration correction, and normalizes to account for spherical geometry.

    The second-order term (-0.5 * |V|² * P * t²) accounts for the curvature
    of the celestial sphere, significantly improving accuracy for high
    proper motion stars (e.g., Barnard's Star) over century-scale intervals.

    Limitations:
    - Ignores radial velocity (parallax causes small position shift)
    - Assumes constant proper motion (real stars accelerate slightly)
    - No annual parallax correction (distance effect negligible for distant stars)
    Typical error: <0.01 arcsec over ±100 years, <1 arcsec over ±500 years
    For research astronomy, use SIMBAD/Gaia catalogs.

Data Sources:
- Positions (RA/Dec at J2000.0): ESA Hipparcos Catalog (ESA SP-1200, 1997)
- Proper motions: van Leeuwen 2007 new Hipparcos reduction (A&A 474, 653-664)
  Independently verified via CDS/VizieR catalog I/311/hip2
- Star names: IAU Working Group on Star Names (WGSN, 2022)
- Visual magnitudes: Hipparcos photometry

References:
- Hipparcos Catalog Vol. 1, Section 1.5.5 (ESA SP-1200, 1997)
- van Leeuwen F., 2007, A&A 474, 653-664 (new Hipparcos reduction)
- IAU 2006 Precession: Capitaine et al. A&A 412, 567-586 (2003)
- IAU WGSN: https://www.iau.org/public/themes/naming_stars/
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple

from skyfield.api import Star
from skyfield.framelib import ecliptic_frame, ecliptic_J2000_frame

from .constants import (
    REGULUS,
    SPICA_STAR,
    ALGOL,
    SIRIUS,
    ALDEBARAN,
    ANTARES,
    VEGA,
    POLARIS,
    FOMALHAUT,
    BETELGEUSE,
    RIGEL,
    PROCYON,
    CAPELLA,
    ARCTURUS,
    DENEB,
    POLLUX,
    CASTOR,
    ALTAIR,
    ACHERNAR,
    CANOPUS,
    ACRUX,
    MIMOSA,
    GACRUX,
    HADAR,
    RIGIL_KENT,
    SHAULA,
    BELLATRIX,
    ELNATH,
    MIRA,
    ALNILAM,
    ALNITAK,
    MINTAKA,
    SAIPH,
    DIPHDA,
    ALPHARD,
    RASALHAGUE,
    ETAMIN,
    KOCHAB,
    ALKAID,
    DUBHE,
    MERAK,
    ALIOTH,
    MIZAR,
    ALCOR,
    VINDEMIATRIX,
    ZUBENELGENUBI,
    ZUBENESCHAMALI,
    UNUKALHAI,
    ALGIEBA,
    DENEBOLA,
    MARKAB,
    SCHEAT,
    ALCYONE,
    ALGORAB,
    ALPHECCA,
    DENEB_ALGEDI,
    ASTEROPE,
    CELAENO,
    ELECTRA,
    MAIA,
    MEROPE,
    TAYGETA,
    ATLAS,
    PLEIONE,
    PRIMA_HYADUM,
    SECUNDA_HYADUM,
    THETA_TAURI,
    AIN,
    MEISSA,
    PHECDA,
    MEGREZ,
    DELTA_CRUCIS,
    MENKENT,
    MUHLIFAIN,
    EPSILON_CENTAURI,
    ETA_CENTAURI,
    ZETA_CENTAURI,
    SARGAS,
    DSCHUBBA,
    GRAFFIAS,
    LESATH,
    ZOSMA,
    HAMAL,
    SHERATAN,
    MESARTHIM,
    ACUBENS,
    TARF,
    ASELLUS_BOREALIS,
    ASELLUS_AUSTRALIS,
    KAUS_AUSTRALIS,
    NUNKI,
    KAUS_MEDIA,
    KAUS_BOREALIS,
    ASCELLA,
    ALGEDI,
    DABIH,
    NASHIRA,
    SADALSUUD,
    SADALMELIK,
    SKAT,
    ETA_PISCIUM,
    ALRESCHA,
    ALPHERATZ,
    ALGENIB,
    PROPUS,
    ALHENA,
    TEJAT,
    WASAT,
    ADHARA,
    WEZEN,
    THUBAN,
    RASALGETHI,
    ALBIREO,
    MIRACH,
    ALMACH,
    MENKAR,
    FLG_SWIEPH,
    FLG_SPEED,
    FLG_SPEED3,
    FLG_NOABERR,
    FLG_NOGDEFL,
    FLG_EQUATORIAL,
    FLG_J2000,
    FLG_NONUT,
    FLG_SIDEREAL,
    FLG_MOSEPH,
    FLG_XYZ,
    FLG_RADIANS,
    FLG_TRUEPOS,
    FLG_TOPOCTR,
    J2000,
    DAYS_PER_JULIAN_YEAR,
)
from .utils import cotrans_sp
from .cache import get_true_obliquity, get_mean_obliquity
from .exceptions import Error


@dataclass
class StarData:
    """
    Fixed star astrometric data (ICRS J2000.0 epoch).

    Attributes:
        ra_j2000: Right Ascension at J2000.0 in degrees
        dec_j2000: Declination at J2000.0 in degrees
        pm_ra: Proper motion in RA (arcsec/year, includes cos(dec) factor)
        pm_dec: Proper motion in Dec (arcsec/year)
        parallax_mas: Trigonometric parallax in milliarcseconds (mas).
            Used to compute stellar distance: dist_AU = 206265000 / parallax_mas.
            A value of 0.0 means unknown parallax (distance defaults to 100000 AU).
        radial_km_per_s: Radial velocity in km/s (positive = receding).
            Affects distance variation over time and speed_dist output.

    Note:
        Proper motions are applied using rigorous space motion approach
        with second-order Taylor expansion (Hipparcos Vol. 1, Section 1.5.5).
        Accurate to <0.01 arcsec over ±100 years, <1 arcsec over ±500 years.
    """

    ra_j2000: float
    dec_j2000: float
    pm_ra: float
    pm_dec: float
    parallax_mas: float = 0.0
    radial_km_per_s: float = 0.0


@dataclass
class StarCatalogEntry:
    """
    Extended catalog entry for fixstar2 functions.

    Attributes:
        id: Internal star ID
        name: Traditional star name (e.g. "Regulus")
        nomenclature: Bayer/Flamsteed designation (e.g. "alLeo", "alVir")
        hip_number: Hipparcos catalog number (e.g. 49669)
        data: Astrometric data for position calculation
        magnitude: Visual magnitude (apparent brightness)
    """

    id: int
    name: str
    nomenclature: str
    hip_number: int
    data: StarData
    magnitude: float


def list_fixed_stars() -> Tuple[StarCatalogEntry, ...]:
    """Return the fixed-star catalog entries."""
    return tuple(STAR_CATALOG)


def propagate_proper_motion(
    ra_epoch: float,
    dec_epoch: float,
    pm_ra_cosdec: float,
    pm_dec: float,
    from_jd: float,
    to_jd: float,
) -> Tuple[float, float]:
    """
    Propagate star position using proper motion between two epochs.

    This function applies the standard astrometric proper motion formula to
    calculate a star's position at any target epoch, given its position and
    proper motion at a reference epoch. This is the fundamental calculation
    for converting Hipparcos positions (at J1991.25) to any other epoch.

    The formula includes the cos(dec) correction that accounts for the
    convergence of right ascension circles toward the celestial poles:

        RA(t) = RA(t0) + pm_ra_cosdec * dt / cos(dec)
        Dec(t) = Dec(t0) + pm_dec * dt

    where dt is the time difference in Julian years (365.25 days).

    Args:
        ra_epoch: Right Ascension at reference epoch (degrees, 0-360)
        dec_epoch: Declination at reference epoch (degrees, -90 to +90)
        pm_ra_cosdec: Proper motion in RA including cos(dec) factor
                      (arcseconds/year). This is mu_alpha* in Hipparcos notation.
        pm_dec: Proper motion in Declination (arcseconds/year).
                This is mu_delta in Hipparcos notation.
        from_jd: Reference epoch as Julian Day (e.g., J1991_25 = 2448349.0625)
        to_jd: Target epoch as Julian Day (e.g., J2000 = 2451545.0)

    Returns:
        Tuple[float, float]: (ra_target, dec_target) position at target epoch
            in degrees. RA is normalized to [0, 360), Dec is in [-90, +90].

    Example:
        # Propagate Sirius from J1991.25 to J2000.0
        >>> from libephemeris.fixed_stars import propagate_proper_motion
        >>> from libephemeris.constants import J1991_25, J2000
        >>> # Sirius at J1991.25 (from Hipparcos catalog)
        >>> ra_1991, dec_1991 = propagate_proper_motion(
        ...     ra_epoch=101.289,  # RA at J1991.25
        ...     dec_epoch=-16.713,  # Dec at J1991.25
        ...     pm_ra_cosdec=-0.54601,  # -546.01 mas/yr as arcsec/yr
        ...     pm_dec=-1.22307,  # -1223.07 mas/yr as arcsec/yr
        ...     from_jd=J1991_25,
        ...     to_jd=J2000,
        ... )

    Notes:
        - pm_ra_cosdec should already include the cos(dec) factor, as is
          standard in the Hipparcos catalog (column "pmRA" or mu_alpha*).
          This represents actual angular motion on the sky in the RA direction.

        - For stars near the poles (|dec| > 89.9 degrees), numerical precision
          may degrade due to the 1/cos(dec) factor. However, proper motion
          values are also smaller in the RA direction for polar stars.

        - This function uses the linear (first-order) proper motion model.
          For high proper motion stars over very long time spans (>500 years),
          the rigorous space motion approach (second-order Taylor expansion)
          is more accurate. See Hipparcos Vol. 1, Section 1.5.5.

        - The time difference is computed in Julian years (exactly 365.25 days).

    References:
        - Hipparcos and Tycho Catalogues, ESA SP-1200, Vol. 1, Section 1.2
        - The Hipparcos proper motion reference epoch is J1991.25 (JD 2448349.0625)
        - IAU SOFA Library: Proper motion transformations
    """
    import math

    # Time difference in Julian years
    dt_years = (to_jd - from_jd) / DAYS_PER_JULIAN_YEAR

    # Convert proper motions from arcsec/year to degrees/year
    pm_ra_deg_per_year = pm_ra_cosdec / 3600.0
    pm_dec_deg_per_year = pm_dec / 3600.0

    # Calculate declination at target epoch (simple linear motion)
    dec_target = dec_epoch + pm_dec_deg_per_year * dt_years

    # Clamp declination to valid range [-90, +90]
    if dec_target > 90.0:
        dec_target = 90.0
    elif dec_target < -90.0:
        dec_target = -90.0

    # Calculate RA at target epoch with cos(dec) correction
    # pm_ra_cosdec already includes cos(dec), so we divide by cos(dec)
    # to get the actual change in RA angle
    cos_dec = math.cos(math.radians(dec_epoch))

    # Avoid division by zero at poles (|dec| = 90 degrees)
    # At the poles, RA is undefined, so we keep it unchanged
    if abs(cos_dec) > 1e-10:
        ra_target = ra_epoch + (pm_ra_deg_per_year / cos_dec) * dt_years
    else:
        ra_target = ra_epoch

    # Normalize RA to [0, 360)
    ra_target = ra_target % 360.0
    if ra_target < 0:
        ra_target += 360.0

    return ra_target, dec_target


# Extended star catalog, built from the generated data table
# (libephemeris/star_catalog_gen.py). The table is produced by
# scripts/build_star_catalog_v2.py exclusively from public sources:
# Hipparcos New Reduction (van Leeuwen 2007, VizieR I/311) astrometry,
# ESA 1997 Hipparcos magnitudes (I/239), the Kostjuk (2004)
# Bayer/Flamsteed cross index (IV/27A), XHIP radial velocities
# (V/137D) and IAU WGSN proper names. Internal IDs of the library's
# original curated stars are preserved by the generator.
from .star_catalog_gen import STAR_ROWS as _STAR_ROWS  # noqa: E402

STAR_CATALOG: List[StarCatalogEntry] = [
    StarCatalogEntry(
        id=_r[0],
        name=_r[1],
        nomenclature=_r[2],
        hip_number=_r[3],
        data=StarData(
            ra_j2000=_r[4],
            dec_j2000=_r[5],
            pm_ra=_r[6],
            pm_dec=_r[7],
            parallax_mas=_r[8],
            radial_km_per_s=_r[9],
        ),
        magnitude=_r[10],
    )
    for _r in _STAR_ROWS
]

# Fixed star catalog (J2000.0 ICRS coordinates from Hipparcos)
# Legacy format for backward compatibility
FIXED_STARS = {entry.id: entry.data for entry in STAR_CATALOG}

# Build lookup from canonical name to star ID
_STAR_NAME_TO_ID = {entry.name.upper(): entry.id for entry in STAR_CATALOG}


# =============================================================================
# BAYER DESIGNATION PARSING
# =============================================================================
# Maps Greek letter names to their 2-letter Bayer abbreviations used in
# the nomenclature field of STAR_CATALOG entries.
#
# Format: "Alpha Leonis" -> "al" + "Leo" -> "alLeo"
# =============================================================================

# Greek letter names to 2-letter abbreviations
# These match the nomenclature format used in STAR_CATALOG
GREEK_LETTER_ABBREV: dict[str, str] = {
    "ALPHA": "al",
    "BETA": "be",
    "GAMMA": "ga",
    "DELTA": "de",
    "EPSILON": "ep",
    "ZETA": "ze",
    "ETA": "et",
    "THETA": "th",
    "IOTA": "io",
    "KAPPA": "ka",
    "LAMBDA": "la",
    "MU": "mu",
    "NU": "nu",
    "XI": "xi",
    "OMICRON": "omi",
    "PI": "pi",
    "RHO": "rh",
    "SIGMA": "si",
    "TAU": "ta",
    "UPSILON": "up",
    "PHI": "ph",
    "CHI": "ch",
    "PSI": "ps",
    "OMEGA": "om",  # omicron is "omi"; only omega abbreviates to "om"
}

# Constellation names (genitive and nominative forms) to 3-letter IAU abbreviations
# Both forms map to the same abbreviation to support "Alpha Leonis" and "Alpha Leo"
CONSTELLATION_ABBREV: dict[str, str] = {
    # Andromeda
    "ANDROMEDAE": "And",
    "ANDROMEDA": "And",
    # Antlia
    "ANTLIAE": "Ant",
    "ANTLIA": "Ant",
    # Apus
    "APODIS": "Aps",
    "APUS": "Aps",
    # Aquarius
    "AQUARII": "Aqr",
    "AQUARIUS": "Aqr",
    # Aquila
    "AQUILAE": "Aql",
    "AQUILA": "Aql",
    # Ara
    "ARAE": "Ara",
    "ARA": "Ara",
    # Aries
    "ARIETIS": "Ari",
    "ARIES": "Ari",
    # Auriga
    "AURIGAE": "Aur",
    "AURIGA": "Aur",
    # Bootes
    "BOOTIS": "Boo",
    "BOOTES": "Boo",
    # Caelum
    "CAELI": "Cae",
    "CAELUM": "Cae",
    # Camelopardalis
    "CAMELOPARDALIS": "Cam",
    # Cancer
    "CANCRI": "Cnc",
    "CANCER": "Cnc",
    # Canes Venatici
    "CANUM VENATICORUM": "CVn",
    "CANES VENATICI": "CVn",
    # Canis Major
    "CANIS MAJORIS": "CMa",
    "CANIS MAJOR": "CMa",
    # Canis Minor
    "CANIS MINORIS": "CMi",
    "CANIS MINOR": "CMi",
    # Capricornus
    "CAPRICORNI": "Cap",
    "CAPRICORNUS": "Cap",
    # Carina
    "CARINAE": "Car",
    "CARINA": "Car",
    # Cassiopeia
    "CASSIOPEIAE": "Cas",
    "CASSIOPEIA": "Cas",
    # Centaurus
    "CENTAURI": "Cen",
    "CENTAURUS": "Cen",
    # Cepheus
    "CEPHEI": "Cep",
    "CEPHEUS": "Cep",
    # Cetus
    "CETI": "Cet",
    "CETUS": "Cet",
    # Chamaeleon
    "CHAMAELEONTIS": "Cha",
    "CHAMAELEON": "Cha",
    # Circinus
    "CIRCINI": "Cir",
    "CIRCINUS": "Cir",
    # Columba
    "COLUMBAE": "Col",
    "COLUMBA": "Col",
    # Coma Berenices
    "COMAE BERENICES": "Com",
    "COMA BERENICES": "Com",
    # Corona Australis
    "CORONAE AUSTRALIS": "CrA",
    "CORONA AUSTRALIS": "CrA",
    # Corona Borealis
    "CORONAE BOREALIS": "CrB",
    "CORONA BOREALIS": "CrB",
    # Corvus
    "CORVI": "Crv",
    "CORVUS": "Crv",
    # Crater
    "CRATERIS": "Crt",
    "CRATER": "Crt",
    # Crux
    "CRUCIS": "Cru",
    "CRUX": "Cru",
    # Cygnus
    "CYGNI": "Cyg",
    "CYGNUS": "Cyg",
    # Delphinus
    "DELPHINI": "Del",
    "DELPHINUS": "Del",
    # Dorado
    "DORADUS": "Dor",
    "DORADO": "Dor",
    # Draco
    "DRACONIS": "Dra",
    "DRACO": "Dra",
    # Equuleus
    "EQUULEI": "Equ",
    "EQUULEUS": "Equ",
    # Eridanus
    "ERIDANI": "Eri",
    "ERIDANUS": "Eri",
    # Fornax
    "FORNACIS": "For",
    "FORNAX": "For",
    # Gemini
    "GEMINORUM": "Gem",
    "GEMINI": "Gem",
    # Grus
    "GRUIS": "Gru",
    "GRUS": "Gru",
    # Hercules
    "HERCULIS": "Her",
    "HERCULES": "Her",
    # Horologium
    "HOROLOGII": "Hor",
    "HOROLOGIUM": "Hor",
    # Hydra
    "HYDRAE": "Hya",
    "HYDRA": "Hya",
    # Hydrus
    "HYDRI": "Hyi",
    "HYDRUS": "Hyi",
    # Indus
    "INDI": "Ind",
    "INDUS": "Ind",
    # Lacerta
    "LACERTAE": "Lac",
    "LACERTA": "Lac",
    # Leo
    "LEONIS": "Leo",
    "LEO": "Leo",
    # Leo Minor
    "LEONIS MINORIS": "LMi",
    "LEO MINOR": "LMi",
    # Lepus
    "LEPORIS": "Lep",
    "LEPUS": "Lep",
    # Libra
    "LIBRAE": "Lib",
    "LIBRA": "Lib",
    # Lupus
    "LUPI": "Lup",
    "LUPUS": "Lup",
    # Lynx
    "LYNCIS": "Lyn",
    "LYNX": "Lyn",
    # Lyra
    "LYRAE": "Lyr",
    "LYRA": "Lyr",
    # Mensa
    "MENSAE": "Men",
    "MENSA": "Men",
    # Microscopium
    "MICROSCOPII": "Mic",
    "MICROSCOPIUM": "Mic",
    # Monoceros
    "MONOCEROTIS": "Mon",
    "MONOCEROS": "Mon",
    # Musca
    "MUSCAE": "Mus",
    "MUSCA": "Mus",
    # Norma
    "NORMAE": "Nor",
    "NORMA": "Nor",
    # Octans
    "OCTANTIS": "Oct",
    "OCTANS": "Oct",
    # Ophiuchus
    "OPHIUCHI": "Oph",
    "OPHIUCHUS": "Oph",
    # Orion
    "ORIONIS": "Ori",
    "ORION": "Ori",
    # Pavo
    "PAVONIS": "Pav",
    "PAVO": "Pav",
    # Pegasus
    "PEGASI": "Peg",
    "PEGASUS": "Peg",
    # Perseus
    "PERSEI": "Per",
    "PERSEUS": "Per",
    # Phoenix
    "PHOENICIS": "Phe",
    "PHOENIX": "Phe",
    # Pictor
    "PICTORIS": "Pic",
    "PICTOR": "Pic",
    # Pisces
    "PISCIUM": "Psc",
    "PISCES": "Psc",
    # Piscis Austrinus
    "PISCIS AUSTRINI": "PsA",
    "PISCIS AUSTRINUS": "PsA",
    # Puppis
    "PUPPIS": "Pup",
    # Pyxis
    "PYXIDIS": "Pyx",
    "PYXIS": "Pyx",
    # Reticulum
    "RETICULI": "Ret",
    "RETICULUM": "Ret",
    # Sagitta
    "SAGITTAE": "Sge",
    "SAGITTA": "Sge",
    # Sagittarius
    "SAGITTARII": "Sgr",
    "SAGITTARIUS": "Sgr",
    # Scorpius
    "SCORPII": "Sco",
    "SCORPIUS": "Sco",
    # Sculptor
    "SCULPTORIS": "Scl",
    "SCULPTOR": "Scl",
    # Scutum
    "SCUTI": "Sct",
    "SCUTUM": "Sct",
    # Serpens
    "SERPENTIS": "Ser",
    "SERPENS": "Ser",
    # Sextans
    "SEXTANTIS": "Sex",
    "SEXTANS": "Sex",
    # Taurus
    "TAURI": "Tau",
    "TAURUS": "Tau",
    # Telescopium
    "TELESCOPII": "Tel",
    "TELESCOPIUM": "Tel",
    # Triangulum
    "TRIANGULI": "Tri",
    "TRIANGULUM": "Tri",
    # Triangulum Australe
    "TRIANGULI AUSTRALIS": "TrA",
    "TRIANGULUM AUSTRALE": "TrA",
    # Tucana
    "TUCANAE": "Tuc",
    "TUCANA": "Tuc",
    # Ursa Major
    "URSAE MAJORIS": "UMa",
    "URSA MAJOR": "UMa",
    # Ursa Minor
    "URSAE MINORIS": "UMi",
    "URSA MINOR": "UMi",
    # Vela
    "VELORUM": "Vel",
    "VELA": "Vel",
    # Virgo
    "VIRGINIS": "Vir",
    "VIRGO": "Vir",
    # Volans
    "VOLANTIS": "Vol",
    "VOLANS": "Vol",
    # Vulpecula
    "VULPECULAE": "Vul",
    "VULPECULA": "Vul",
}


def _parse_bayer_designation(designation: str) -> str | None:
    """
    Parse a Bayer designation into nomenclature format.

    Converts designations like "Alpha Leonis", "Beta Persei", "Gamma Virginis"
    into the nomenclature format used in STAR_CATALOG (e.g., "alLeo", "bePer", "gaVir").

    Args:
        designation: Bayer designation string (e.g., "Alpha Leonis", "Beta Persei")

    Returns:
        Nomenclature string if parsed successfully (e.g., "alLeo"), None otherwise

    Examples:
        >>> _parse_bayer_designation("Alpha Leonis")
        'alLeo'
        >>> _parse_bayer_designation("Beta Persei")
        'bePer'
        >>> _parse_bayer_designation("Gamma Virginis")
        'gaVir'
        >>> _parse_bayer_designation("Invalid Name")
        None
    """
    if not designation:
        return None

    # Normalize: uppercase for matching
    parts = designation.upper().strip().split()

    if len(parts) < 2:
        return None

    # First part should be a Greek letter name
    greek_letter = parts[0]
    if greek_letter not in GREEK_LETTER_ABBREV:
        return None

    letter_abbrev = GREEK_LETTER_ABBREV[greek_letter]

    # Remaining parts form the constellation name
    constellation_name = " ".join(parts[1:])

    # Look up constellation abbreviation
    if constellation_name not in CONSTELLATION_ABBREV:
        return None

    const_abbrev = CONSTELLATION_ABBREV[constellation_name]

    # Build nomenclature: letter abbreviation (lowercase) + constellation abbreviation
    # e.g., "al" + "Leo" = "alLeo"
    return letter_abbrev + const_abbrev


def _parse_flamsteed_designation(designation: str) -> str | None:
    """
    Parse a Flamsteed designation into a normalized format for STAR_ALIASES lookup.

    Converts designations like "32 Leonis", "87 Virginis", "21 Tauri"
    into the format used in STAR_ALIASES (e.g., "32 LEO", "87 VIR", "21 TAU").

    Flamsteed designations consist of a number followed by the constellation
    name in genitive (Latin) form. This function normalizes them to the format
    "{number} {3-letter-abbrev}" used in STAR_ALIASES.

    Args:
        designation: Flamsteed designation string (e.g., "32 Leonis", "87 Virginis")

    Returns:
        Normalized alias string if parsed successfully (e.g., "32 LEO"), None otherwise

    Examples:
        >>> _parse_flamsteed_designation("32 Leonis")
        '32 LEO'
        >>> _parse_flamsteed_designation("87 Virginis")
        '87 VIR'
        >>> _parse_flamsteed_designation("21 Tauri")
        '21 TAU'
        >>> _parse_flamsteed_designation("Invalid Name")
        None
    """
    if not designation:
        return None

    # Normalize: strip and split
    parts = designation.strip().split()

    if len(parts) < 2:
        return None

    # First part should be a number (Flamsteed number)
    number_str = parts[0]
    if not number_str.isdigit():
        return None

    # Remaining parts form the constellation name
    constellation_name = " ".join(parts[1:]).upper()

    # Look up constellation abbreviation
    if constellation_name not in CONSTELLATION_ABBREV:
        return None

    const_abbrev = CONSTELLATION_ABBREV[constellation_name]

    # Build normalized alias: number + constellation abbreviation (uppercase)
    # e.g., "32" + "LEO" = "32 LEO"
    return f"{number_str} {const_abbrev.upper()}"


# STAR_ALIASES: Maps alternative star names to canonical SE_* constant IDs
# Includes: common names, Bayer designations (full and abbreviated),
# Flamsteed numbers, Arabic names, Latin names, Greek transliterations
STAR_ALIASES: dict[str, int] = {
    # ======== REGULUS (Alpha Leonis) - ROYAL STAR ========
    "COR LEONIS": REGULUS,
    "ALPHA LEONIS": REGULUS,
    "ALPHA LEO": REGULUS,
    "87 LEO": REGULUS,
    "α LEO": REGULUS,
    "QALB AL-ASAD": REGULUS,
    "BASILISKOS": REGULUS,
    "REX": REGULUS,
    "KALB AL-ASAD": REGULUS,
    "ALLEO": REGULUS,
    "WATCHER OF THE NORTH": REGULUS,
    "VENANT": REGULUS,
    # ======== SPICA (Alpha Virginis) ========
    "ALPHA VIRGINIS": SPICA_STAR,
    "ALPHA VIR": SPICA_STAR,
    "67 VIR": SPICA_STAR,
    "α VIR": SPICA_STAR,
    "AZIMECH": SPICA_STAR,
    "ALARAPH": SPICA_STAR,
    "ALVIR": SPICA_STAR,
    "SUNBULA": SPICA_STAR,
    "VIRGIN'S SPIKE": SPICA_STAR,
    "ARISTA": SPICA_STAR,
    # ======== ALGOL (Beta Persei) ========
    "DEMON STAR": ALGOL,
    "BETA PERSEI": ALGOL,
    "BETA PER": ALGOL,
    "26 PER": ALGOL,
    "β PER": ALGOL,
    "GORGONEA PRIMA": ALGOL,
    "RA'S AL-GHUL": ALGOL,
    "RAS AL-GHUL": ALGOL,
    "BEPER": ALGOL,
    "HEAD OF THE GHOUL": ALGOL,
    "GORGON'S HEAD": ALGOL,
    # ======== SIRIUS (Alpha Canis Majoris) ========
    "DOG STAR": SIRIUS,
    "ALPHA CANIS MAJORIS": SIRIUS,
    "ALPHA CMA": SIRIUS,
    "α CMA": SIRIUS,
    "9 CMA": SIRIUS,
    "CANICULA": SIRIUS,
    "ASCHERE": SIRIUS,
    "ALCMA": SIRIUS,
    "AL-SHIRA": SIRIUS,
    "SOTHIS": SIRIUS,
    # ======== ALDEBARAN (Alpha Tauri) - ROYAL STAR ========
    "EYE OF TAURUS": ALDEBARAN,
    "ALPHA TAURI": ALDEBARAN,
    "ALPHA TAU": ALDEBARAN,
    "87 TAU": ALDEBARAN,
    "α TAU": ALDEBARAN,
    "ALTAU": ALDEBARAN,
    "PARILICIUM": ALDEBARAN,
    "AL-DABARAN": ALDEBARAN,
    "FOLLOWER": ALDEBARAN,
    "ROHINI": ALDEBARAN,
    "WATCHER OF THE EAST": ALDEBARAN,
    "TASCHETER": ALDEBARAN,
    # ======== ANTARES (Alpha Scorpii) - ROYAL STAR ========
    "RIVAL OF MARS": ANTARES,
    "ALPHA SCORPII": ANTARES,
    "ALPHA SCO": ANTARES,
    "21 SCO": ANTARES,
    "α SCO": ANTARES,
    "ALSCO": ANTARES,
    "COR SCORPII": ANTARES,
    "CALB AL-AKRAB": ANTARES,
    "HEART OF SCORPION": ANTARES,
    "JYESHTHA": ANTARES,
    "WATCHER OF THE WEST": ANTARES,
    "SATEVIS": ANTARES,
    # ======== VEGA (Alpha Lyrae) ========
    "HARP STAR": VEGA,
    "ALPHA LYRAE": VEGA,
    "ALPHA LYR": VEGA,
    "3 LYR": VEGA,
    "α LYR": VEGA,
    "ALLYR": VEGA,
    "WEGA": VEGA,
    "AL-NASR AL-WAQI": VEGA,
    "FIDIS": VEGA,
    "ABHIJIT": VEGA,
    # ======== POLARIS (Alpha Ursae Minoris) ========
    "NORTH STAR": POLARIS,
    "POLE STAR": POLARIS,
    "ALPHA URSAE MINORIS": POLARIS,
    "ALPHA UMI": POLARIS,
    "1 UMI": POLARIS,
    "α UMI": POLARIS,
    "ALUMI": POLARIS,
    "CYNOSURA": POLARIS,
    "LODESTAR": POLARIS,
    "STELLA POLARIS": POLARIS,
    # ======== FOMALHAUT (Alpha Piscis Austrini) - ROYAL STAR ========
    "FISH'S MOUTH": FOMALHAUT,
    "ALPHA PISCIS AUSTRINI": FOMALHAUT,
    "ALPHA PSA": FOMALHAUT,
    "24 PSA": FOMALHAUT,
    "α PSA": FOMALHAUT,
    "ALPSA": FOMALHAUT,
    "FUM AL-HUT": FOMALHAUT,
    "OS PISCIS MERIDIANI": FOMALHAUT,
    "LONELY STAR": FOMALHAUT,
    "HASTORANG": FOMALHAUT,
    "WATCHER OF THE SOUTH": FOMALHAUT,
    # ======== BETELGEUSE (Alpha Orionis) ========
    "ARMPIT OF ORION": BETELGEUSE,
    "ALPHA ORIONIS": BETELGEUSE,
    "ALPHA ORI": BETELGEUSE,
    "58 ORI": BETELGEUSE,
    "α ORI": BETELGEUSE,
    "ALORI": BETELGEUSE,
    "BETELGEUZE": BETELGEUSE,
    "IBT AL-JAUZAH": BETELGEUSE,
    "ARDRA": BETELGEUSE,
    "YAD AL-JAWZA": BETELGEUSE,
    # ======== RIGEL (Beta Orionis) ========
    "FOOT OF ORION": RIGEL,
    "BETA ORIONIS": RIGEL,
    "BETA ORI": RIGEL,
    "19 ORI": RIGEL,
    "β ORI": RIGEL,
    "BEORI": RIGEL,
    "ALGEBAR": RIGEL,
    "RIJL JAUZAH AL-YUSRA": RIGEL,
    "ORION'S FOOT": RIGEL,
    # ======== PROCYON (Alpha Canis Minoris) ========
    "LITTLE DOG STAR": PROCYON,
    "ALPHA CANIS MINORIS": PROCYON,
    "ALPHA CMI": PROCYON,
    "10 CMI": PROCYON,
    "α CMI": PROCYON,
    "ALCMI": PROCYON,
    "ANTECANIS": PROCYON,
    "ELGOMAISA": PROCYON,
    "AL-GHUMAISA": PROCYON,
    # ======== CAPELLA (Alpha Aurigae) ========
    "SHE-GOAT": CAPELLA,
    "ALPHA AURIGAE": CAPELLA,
    "ALPHA AUR": CAPELLA,
    "13 AUR": CAPELLA,
    "α AUR": CAPELLA,
    "ALAUR": CAPELLA,
    "ALHAJOTH": CAPELLA,
    "AMALTHEA": CAPELLA,
    "GOAT STAR": CAPELLA,
    # ======== ARCTURUS (Alpha Bootis) ========
    "BEAR GUARD": ARCTURUS,
    "ALPHA BOOTIS": ARCTURUS,
    "ALPHA BOO": ARCTURUS,
    "16 BOO": ARCTURUS,
    "α BOO": ARCTURUS,
    "ALBOO": ARCTURUS,
    "AL-SIMAK AL-RAMIH": ARCTURUS,
    "HARIS AL-SAMA": ARCTURUS,
    "GUARDIAN OF BEAR": ARCTURUS,
    # ======== DENEB (Alpha Cygni) ========
    "TAIL OF HEN": DENEB,
    "ALPHA CYGNI": DENEB,
    "ALPHA CYG": DENEB,
    "50 CYG": DENEB,
    "α CYG": DENEB,
    "ALCYG": DENEB,
    "DHANAB AD-DAJAJAH": DENEB,
    "ARIDED": DENEB,
    "GALLINA": DENEB,
    # ======== POLLUX (Beta Geminorum) ========
    "TWIN STAR": POLLUX,
    "BETA GEMINORUM": POLLUX,
    "BETA GEM": POLLUX,
    "78 GEM": POLLUX,
    "β GEM": POLLUX,
    "BEGEM": POLLUX,
    "POLYDEUCES": POLLUX,
    "HEAD OF SECOND TWIN": POLLUX,
    # ======== CASTOR (Alpha Geminorum) ========
    "ALPHA GEMINORUM": CASTOR,
    "ALPHA GEM": CASTOR,
    "66 GEM": CASTOR,
    "α GEM": CASTOR,
    "ALGEM": CASTOR,
    "APOLLO": CASTOR,
    "HEAD OF FIRST TWIN": CASTOR,
    # ======== ALTAIR (Alpha Aquilae) ========
    "FLYING EAGLE": ALTAIR,
    "ALPHA AQUILAE": ALTAIR,
    "ALPHA AQL": ALTAIR,
    "53 AQL": ALTAIR,
    "α AQL": ALTAIR,
    "ALAQL": ALTAIR,
    "AL-NASR AL-TAIR": ALTAIR,
    "ATAIR": ALTAIR,
    "SRAVANA": ALTAIR,
    # ======== ACHERNAR (Alpha Eridani) ========
    "END OF RIVER": ACHERNAR,
    "ALPHA ERIDANI": ACHERNAR,
    "ALPHA ERI": ACHERNAR,
    "α ERI": ACHERNAR,
    "ALERI": ACHERNAR,
    "AKHIR AN-NAHR": ACHERNAR,
    "RIVER'S END": ACHERNAR,
    # ======== CANOPUS (Alpha Carinae) ========
    "SHIP'S PILOT": CANOPUS,
    "ALPHA CARINAE": CANOPUS,
    "ALPHA CAR": CANOPUS,
    "α CAR": CANOPUS,
    "ALCAR": CANOPUS,
    "SUHAIL": CANOPUS,
    "SUHAYL": CANOPUS,
    "AGASTYA": CANOPUS,
    # ======== ACRUX (Alpha Crucis) ========
    "ALPHA CRUCIS": ACRUX,
    "ALPHA CRU": ACRUX,
    "α CRU": ACRUX,
    "ALCRU": ACRUX,
    "CRUX ALPHA": ACRUX,
    "STAR OF BETHLEHEM": ACRUX,
    # ======== MIMOSA (Beta Crucis) ========
    "BETA CRUCIS": MIMOSA,
    "BETA CRU": MIMOSA,
    "β CRU": MIMOSA,
    "BECRU": MIMOSA,
    "BECRUX": MIMOSA,
    "CRUX BETA": MIMOSA,
    # ======== GACRUX (Gamma Crucis) ========
    "GAMMA CRUCIS": GACRUX,
    "GAMMA CRU": GACRUX,
    "γ CRU": GACRUX,
    "GACRU": GACRUX,
    "CRUX GAMMA": GACRUX,
    "RUBIDEA": GACRUX,
    # ======== DELTA CRUCIS (Delta Crucis) ========
    "DELTA CRUCIS": DELTA_CRUCIS,
    "DELTA CRU": DELTA_CRUCIS,
    "δ CRU": DELTA_CRUCIS,
    "DECRU": DELTA_CRUCIS,
    "CRUX DELTA": DELTA_CRUCIS,
    "DECRUX": DELTA_CRUCIS,
    # ======== HADAR (Beta Centauri) ========
    "BETA CENTAURI": HADAR,
    "BETA CEN": HADAR,
    "β CEN": HADAR,
    "BECEN": HADAR,
    "AGENA": HADAR,
    "KNEE OF CENTAUR": HADAR,
    # ======== RIGIL KENTAURUS (Alpha Centauri) ========
    "ALPHA CENTAURI": RIGIL_KENT,
    "ALPHA CEN": RIGIL_KENT,
    "α CEN": RIGIL_KENT,
    "ALCEN": RIGIL_KENT,
    # Toliman is the IAU name of alpha Cen B (catalog id 1000846); the
    # exact catalog name takes precedence, so this legacy alias only
    # documents the historical usage.
    "TOLIMAN": 1000846,
    "RIGIL KENT": RIGIL_KENT,
    "FOOT OF CENTAUR": RIGIL_KENT,
    "BUNGULA": RIGIL_KENT,
    # ======== MENKENT (Theta Centauri) ========
    "THETA CENTAURI": MENKENT,
    "THETA CEN": MENKENT,
    "θ CEN": MENKENT,
    "THCEN": MENKENT,
    "HARATAN": MENKENT,
    # ======== MUHLIFAIN (Gamma Centauri) ========
    "GAMMA CENTAURI": MUHLIFAIN,
    "GAMMA CEN": MUHLIFAIN,
    "γ CEN": MUHLIFAIN,
    "GACEN": MUHLIFAIN,
    # ======== EPSILON CENTAURI ========
    "EPSILON CENTAURI": EPSILON_CENTAURI,
    "EPSILON CEN": EPSILON_CENTAURI,
    "ε CEN": EPSILON_CENTAURI,
    "EPCEN": EPSILON_CENTAURI,
    # ======== ETA CENTAURI ========
    "ETA CENTAURI": ETA_CENTAURI,
    "ETA CEN": ETA_CENTAURI,
    "η CEN": ETA_CENTAURI,
    "ETCEN": ETA_CENTAURI,
    # ======== ZETA CENTAURI ========
    "ZETA CENTAURI": ZETA_CENTAURI,
    "ZETA CEN": ZETA_CENTAURI,
    "ζ CEN": ZETA_CENTAURI,
    "ZECEN": ZETA_CENTAURI,
    "ALNAIR": ZETA_CENTAURI,
    # ======== SHAULA (Lambda Scorpii) ========
    "LAMBDA SCORPII": SHAULA,
    "LAMBDA SCO": SHAULA,
    "λ SCO": SHAULA,
    "LASCO": SHAULA,
    "SCORPION'S STING": SHAULA,
    "35 SCO": SHAULA,
    # ======== BELLATRIX (Gamma Orionis) ========
    "AMAZON STAR": BELLATRIX,
    "GAMMA ORIONIS": BELLATRIX,
    "GAMMA ORI": BELLATRIX,
    "24 ORI": BELLATRIX,
    "γ ORI": BELLATRIX,
    "GAORI": BELLATRIX,
    "FEMALE WARRIOR": BELLATRIX,
    # ======== ELNATH (Beta Tauri) ========
    "BETA TAURI": ELNATH,
    "BETA TAU": ELNATH,
    "112 TAU": ELNATH,
    "β TAU": ELNATH,
    "BETAU": ELNATH,
    "NATH": ELNATH,
    "AL-NATH": ELNATH,
    "EL NATH": ELNATH,
    # ======== MIRA (Omicron Ceti) ========
    "WONDERFUL STAR": MIRA,
    "OMICRON CETI": MIRA,
    "OMICRON CET": MIRA,
    "68 CET": MIRA,
    "ο CET": MIRA,
    "OMCET": MIRA,
    "STELLA MIRA": MIRA,
    # ======== ALNILAM (Epsilon Orionis) ========
    "STRING OF PEARLS": ALNILAM,
    "EPSILON ORIONIS": ALNILAM,
    "EPSILON ORI": ALNILAM,
    "46 ORI": ALNILAM,
    "ε ORI": ALNILAM,
    "EPORI": ALNILAM,
    "AL-NIZAM": ALNILAM,
    # ======== ALNITAK (Zeta Orionis) ========
    "GIRDLE": ALNITAK,
    "ZETA ORIONIS": ALNITAK,
    "ZETA ORI": ALNITAK,
    "50 ORI": ALNITAK,
    "ζ ORI": ALNITAK,
    "ZEORI": ALNITAK,
    "AL-NITAK": ALNITAK,
    # ======== MINTAKA (Delta Orionis) ========
    "BELT STAR": MINTAKA,
    "DELTA ORIONIS": MINTAKA,
    "DELTA ORI": MINTAKA,
    "34 ORI": MINTAKA,
    "δ ORI": MINTAKA,
    "DEORI": MINTAKA,
    "MINTAKA": MINTAKA,
    # ======== SAIPH (Kappa Orionis) ========
    "SWORD OF GIANT": SAIPH,
    "KAPPA ORIONIS": SAIPH,
    "KAPPA ORI": SAIPH,
    "53 ORI": SAIPH,
    "κ ORI": SAIPH,
    "KAORI": SAIPH,
    "SAIF AL-JABBAR": SAIPH,
    # ======== MEISSA (Lambda Orionis) ========
    "LAMBDA ORIONIS": MEISSA,
    "LAMBDA ORI": MEISSA,
    "39 ORI": MEISSA,
    "λ ORI": MEISSA,
    "LAORI": MEISSA,
    "HEKA": MEISSA,
    "HEAD OF ORION": MEISSA,
    "AL-MAISAN": MEISSA,
    # ======== DIPHDA (Beta Ceti) ========
    "FROG": DIPHDA,
    "BETA CETI": DIPHDA,
    "BETA CET": DIPHDA,
    "16 CET": DIPHDA,
    "β CET": DIPHDA,
    "BECET": DIPHDA,
    "DENEB KAITOS": DIPHDA,
    "DIFDA AL-THANI": DIPHDA,
    # ======== ALPHARD (Alpha Hydrae) ========
    "SOLITARY ONE": ALPHARD,
    "ALPHA HYDRAE": ALPHARD,
    "ALPHA HYA": ALPHARD,
    "30 HYA": ALPHARD,
    "α HYA": ALPHARD,
    "ALHYA": ALPHARD,
    "COR HYDRAE": ALPHARD,
    "AL-FARD": ALPHARD,
    # ======== RASALHAGUE (Alpha Ophiuchi) ========
    "HEAD OF SERPENT HOLDER": RASALHAGUE,
    "ALPHA OPHIUCHI": RASALHAGUE,
    "ALPHA OPH": RASALHAGUE,
    "55 OPH": RASALHAGUE,
    "α OPH": RASALHAGUE,
    "ALOPH": RASALHAGUE,
    "RAS AL-HAWWA": RASALHAGUE,
    # ======== ETAMIN (Gamma Draconis) ========
    "DRAGON'S HEAD": ETAMIN,
    "GAMMA DRACONIS": ETAMIN,
    "GAMMA DRA": ETAMIN,
    "33 DRA": ETAMIN,
    "γ DRA": ETAMIN,
    "GADORA": ETAMIN,
    "ELTANIN": ETAMIN,
    "AL-TINNIN": ETAMIN,
    # ======== KOCHAB (Beta Ursae Minoris) ========
    "BETA URSAE MINORIS": KOCHAB,
    "BETA UMI": KOCHAB,
    "7 UMI": KOCHAB,
    "β UMI": KOCHAB,
    "BEUMI": KOCHAB,
    "KAUKAB": KOCHAB,
    # ======== ALKAID (Eta Ursae Majoris) ========
    "END OF TAIL": ALKAID,
    "ETA URSAE MAJORIS": ALKAID,
    "ETA UMA": ALKAID,
    "85 UMA": ALKAID,
    "η UMA": ALKAID,
    "ETUMA": ALKAID,
    "BENETNASH": ALKAID,
    "AL-QA'ID": ALKAID,
    # ======== DUBHE (Alpha Ursae Majoris) ========
    "BEAR'S BACK": DUBHE,
    "ALPHA URSAE MAJORIS": DUBHE,
    "ALPHA UMA": DUBHE,
    "50 UMA": DUBHE,
    "α UMA": DUBHE,
    "ALUMA": DUBHE,
    "THAHR AL-DUBB AL-AKBAR": DUBHE,
    # ======== MERAK (Beta Ursae Majoris) ========
    "BETA URSAE MAJORIS": MERAK,
    "BETA UMA": MERAK,
    "48 UMA": MERAK,
    "β UMA": MERAK,
    "BEUMA": MERAK,
    "AL-MARAKK": MERAK,
    # ======== ALIOTH (Epsilon Ursae Majoris) ========
    "EPSILON URSAE MAJORIS": ALIOTH,
    "EPSILON UMA": ALIOTH,
    "77 UMA": ALIOTH,
    "ε UMA": ALIOTH,
    "EPUMA": ALIOTH,
    "ALIATH": ALIOTH,
    # ======== MIZAR (Zeta Ursae Majoris) ========
    "ZETA URSAE MAJORIS": MIZAR,
    "ZETA UMA": MIZAR,
    "79 UMA": MIZAR,
    "ζ UMA": MIZAR,
    "ZEUMA": MIZAR,
    "HORSE AND RIDER": MIZAR,
    # ======== ALCOR (80 Ursae Majoris) ========
    "80 URSAE MAJORIS": ALCOR,
    "80 UMA": ALCOR,
    "G UMA": ALCOR,
    "SAIDAK": ALCOR,
    "SUHA": ALCOR,
    "AL-SAHJA": ALCOR,
    # ======== PHECDA (Gamma Ursae Majoris) ========
    "GAMMA URSAE MAJORIS": PHECDA,
    "GAMMA UMA": PHECDA,
    "64 UMA": PHECDA,
    "γ UMA": PHECDA,
    "GAUMA": PHECDA,
    "PHAD": PHECDA,
    "PHEKDA": PHECDA,
    "PHACD": PHECDA,
    # ======== MEGREZ (Delta Ursae Majoris) ========
    "DELTA URSAE MAJORIS": MEGREZ,
    "DELTA UMA": MEGREZ,
    "69 UMA": MEGREZ,
    "δ UMA": MEGREZ,
    "DEUMA": MEGREZ,
    "KAFFA": MEGREZ,
    # ======== VINDEMIATRIX (Epsilon Virginis) ========
    "GRAPE GATHERER": VINDEMIATRIX,
    "EPSILON VIRGINIS": VINDEMIATRIX,
    "EPSILON VIR": VINDEMIATRIX,
    "47 VIR": VINDEMIATRIX,
    "ε VIR": VINDEMIATRIX,
    "EPVIR": VINDEMIATRIX,
    "ALMUREDIN": VINDEMIATRIX,
    # ======== ZUBENELGENUBI (Alpha Librae) ========
    "SOUTHERN CLAW": ZUBENELGENUBI,
    "ALPHA LIBRAE": ZUBENELGENUBI,
    "ALPHA LIB": ZUBENELGENUBI,
    "9 LIB": ZUBENELGENUBI,
    "α LIB": ZUBENELGENUBI,
    "ALLIB": ZUBENELGENUBI,
    "KIFFA AUSTRALIS": ZUBENELGENUBI,
    # ======== ZUBENESCHAMALI (Beta Librae) ========
    "NORTHERN CLAW": ZUBENESCHAMALI,
    "BETA LIBRAE": ZUBENESCHAMALI,
    "BETA LIB": ZUBENESCHAMALI,
    "27 LIB": ZUBENESCHAMALI,
    "β LIB": ZUBENESCHAMALI,
    "BELIB": ZUBENESCHAMALI,
    "KIFFA BOREALIS": ZUBENESCHAMALI,
    # ======== UNUKALHAI (Alpha Serpentis) ========
    "SERPENT'S NECK": UNUKALHAI,
    "ALPHA SERPENTIS": UNUKALHAI,
    "ALPHA SER": UNUKALHAI,
    "24 SER": UNUKALHAI,
    "α SER": UNUKALHAI,
    "ALSER": UNUKALHAI,
    "COR SERPENTIS": UNUKALHAI,
    # ======== ALGIEBA (Gamma Leonis) ========
    "LION'S MANE": ALGIEBA,
    "GAMMA LEONIS": ALGIEBA,
    "GAMMA LEO": ALGIEBA,
    "41 LEO": ALGIEBA,
    "γ LEO": ALGIEBA,
    "GALEO": ALGIEBA,
    "AL-JABHAH": ALGIEBA,
    # ======== DENEBOLA (Beta Leonis) ========
    "LION'S TAIL": DENEBOLA,
    "BETA LEONIS": DENEBOLA,
    "BETA LEO": DENEBOLA,
    "94 LEO": DENEBOLA,
    "β LEO": DENEBOLA,
    "BELEO": DENEBOLA,
    "DHANAB AL-ASAD": DENEBOLA,
    # ======== MARKAB (Alpha Pegasi) ========
    "SADDLE": MARKAB,
    "ALPHA PEGASI": MARKAB,
    "ALPHA PEG": MARKAB,
    "54 PEG": MARKAB,
    "α PEG": MARKAB,
    "ALPEG": MARKAB,
    "MANKIB AL-FARAS": MARKAB,
    # ======== SCHEAT (Beta Pegasi) ========
    "LEG": SCHEAT,
    "BETA PEGASI": SCHEAT,
    "BETA PEG": SCHEAT,
    "53 PEG": SCHEAT,
    "β PEG": SCHEAT,
    "BEPEG": SCHEAT,
    "SAQ AL-FARAS": SCHEAT,
    # ======== ALCYONE (Eta Tauri - Pleiades) - BEHENIAN ========
    "PLEIADES": ALCYONE,
    "ETA TAURI": ALCYONE,
    "ETA TAU": ALCYONE,
    "25 TAU": ALCYONE,
    "η TAU": ALCYONE,
    "ETTAU": ALCYONE,
    "SEVEN SISTERS": ALCYONE,
    "KIMAH": ALCYONE,
    "AL-THURAYYA": ALCYONE,
    # ======== ALGORAB (Delta Corvi) - BEHENIAN ========
    "CROW'S WING": ALGORAB,
    "DELTA CORVI": ALGORAB,
    "DELTA CRV": ALGORAB,
    "7 CRV": ALGORAB,
    "δ CRV": ALGORAB,
    "DECRV": ALGORAB,
    "AL-GHIRAB": ALGORAB,
    "GIENAH CORVI": ALGORAB,
    # ======== ALPHECCA (Alpha Coronae Borealis) - BEHENIAN ========
    "GEMMA": ALPHECCA,
    "ALPHA CORONAE BOREALIS": ALPHECCA,
    "ALPHA CRB": ALPHECCA,
    "5 CRB": ALPHECCA,
    "α CRB": ALPHECCA,
    "ALCRB": ALPHECCA,
    "GNOSIA STELLA": ALPHECCA,
    "MUNIR AL-FAKKAH": ALPHECCA,
    "ASHTAROTH": ALPHECCA,
    # ======== DENEB ALGEDI (Delta Capricorni) - BEHENIAN ========
    "TAIL OF THE GOAT": DENEB_ALGEDI,
    "DELTA CAPRICORNI": DENEB_ALGEDI,
    "DELTA CAP": DENEB_ALGEDI,
    "49 CAP": DENEB_ALGEDI,
    "δ CAP": DENEB_ALGEDI,
    "DECAP": DENEB_ALGEDI,
    "SCHEDDI": DENEB_ALGEDI,
    "DHANAB AL-JADY": DENEB_ALGEDI,
    # ======== ASTEROPE (21 Tauri - Pleiades) ========
    "21 TAURI": ASTEROPE,
    "21 TAU": ASTEROPE,
    "STEROPE": ASTEROPE,
    "STEROPE I": ASTEROPE,
    # ======== CELAENO (16 Tauri - Pleiades) ========
    "16 TAURI": CELAENO,
    "16 TAU": CELAENO,
    "CELENO": CELAENO,
    # ======== ELECTRA (17 Tauri - Pleiades) ========
    "17 TAURI": ELECTRA,
    "17 TAU": ELECTRA,
    # ======== MAIA (20 Tauri - Pleiades) ========
    "20 TAURI": MAIA,
    "20 TAU": MAIA,
    # ======== MEROPE (23 Tauri - Pleiades) ========
    "23 TAURI": MEROPE,
    "23 TAU": MEROPE,
    # ======== TAYGETA (19 Tauri - Pleiades) ========
    "19 TAURI": TAYGETA,
    "19 TAU": TAYGETA,
    "TAYGETE": TAYGETA,
    # ======== ATLAS (27 Tauri - Pleiades) ========
    "27 TAURI": ATLAS,
    "27 TAU": ATLAS,
    # ======== PLEIONE (28 Tauri - Pleiades) ========
    "28 TAURI": PLEIONE,
    "28 TAU": PLEIONE,
    # ======== PRIMA HYADUM (Gamma Tauri - Hyades) ========
    "GAMMA TAURI": PRIMA_HYADUM,
    "GAMMA TAU": PRIMA_HYADUM,
    "54 TAU": PRIMA_HYADUM,
    "γ TAU": PRIMA_HYADUM,
    "GATAU": PRIMA_HYADUM,
    "HYADUM I": PRIMA_HYADUM,
    "FIRST HYAD": PRIMA_HYADUM,
    # ======== SECUNDA HYADUM (Delta^1 Tauri - Hyades) ========
    "DELTA TAURI": SECUNDA_HYADUM,
    "DELTA TAU": SECUNDA_HYADUM,
    "DELTA1 TAURI": SECUNDA_HYADUM,
    "DELTA1 TAU": SECUNDA_HYADUM,
    "61 TAU": SECUNDA_HYADUM,
    "δ TAU": SECUNDA_HYADUM,
    "DETAU": SECUNDA_HYADUM,
    "HYADUM II": SECUNDA_HYADUM,
    "SECOND HYAD": SECUNDA_HYADUM,
    # ======== THETA TAURI (Theta^2 Tauri - Hyades) ========
    "THETA2 TAURI": THETA_TAURI,
    "THETA2 TAU": THETA_TAURI,
    "78 TAU": THETA_TAURI,
    "θ TAU": THETA_TAURI,
    "THTAU": THETA_TAURI,
    "THETA^2 TAURI": THETA_TAURI,
    "THETA^2 TAU": THETA_TAURI,
    # ======== AIN (Epsilon Tauri - Hyades) ========
    "EPSILON TAURI": AIN,
    "EPSILON TAU": AIN,
    "74 TAU": AIN,
    "ε TAU": AIN,
    "EPTAU": AIN,
    "OCULUS BOREALIS": AIN,
    "BULL'S EYE": AIN,
    # ======== SARGAS (Theta Scorpii) ========
    "THETA SCORPII": SARGAS,
    "THETA SCO": SARGAS,
    "θ SCO": SARGAS,
    "THSCO": SARGAS,
    "GIRTAB": SARGAS,
    "SCORPION'S TAIL": SARGAS,
    # ======== DSCHUBBA (Delta Scorpii) ========
    "DELTA SCORPII": DSCHUBBA,
    "DELTA SCO": DSCHUBBA,
    "δ SCO": DSCHUBBA,
    "DESCO": DSCHUBBA,
    "ICLARCRAU": DSCHUBBA,
    "ICLARKRAV": DSCHUBBA,
    "SCORPION'S FOREHEAD": DSCHUBBA,
    # ======== GRAFFIAS (Beta Scorpii) ========
    "BETA SCORPII": GRAFFIAS,
    "BETA SCO": GRAFFIAS,
    "β SCO": GRAFFIAS,
    "BESCO": GRAFFIAS,
    "ACRAB": GRAFFIAS,
    "AKRAB": GRAFFIAS,
    "ELACRAB": GRAFFIAS,
    # ======== LESATH (Upsilon Scorpii) ========
    "UPSILON SCORPII": LESATH,
    "UPSILON SCO": LESATH,
    "υ SCO": LESATH,
    "UPSCO": LESATH,
    "STINGER": LESATH,
    "34 SCO": LESATH,
    # ======== ZOSMA (Delta Leonis) ========
    "DELTA LEONIS": ZOSMA,
    "DELTA LEO": ZOSMA,
    "68 LEO": ZOSMA,
    "δ LEO": ZOSMA,
    "DELEO": ZOSMA,
    "DHUR": ZOSMA,
    "DUHR": ZOSMA,
    "LION'S HIP": ZOSMA,
    "LION'S BACK": ZOSMA,
    # ======== ZODIACAL CONSTELLATION BRIGHT STARS ========
    # ======== HAMAL (Alpha Arietis) ========
    "ALPHA ARIETIS": HAMAL,
    "ALPHA ARI": HAMAL,
    "13 ARI": HAMAL,
    "α ARI": HAMAL,
    "ALARI": HAMAL,
    "RAM'S HEAD": HAMAL,
    # ======== SHERATAN (Beta Arietis) ========
    "BETA ARIETIS": SHERATAN,
    "BETA ARI": SHERATAN,
    "6 ARI": SHERATAN,
    "β ARI": SHERATAN,
    "BEARI": SHERATAN,
    "SHARATAN": SHERATAN,
    "AL-SHARATAIN": SHERATAN,
    # ======== MESARTHIM (Gamma Arietis) ========
    "GAMMA ARIETIS": MESARTHIM,
    "GAMMA ARI": MESARTHIM,
    "5 ARI": MESARTHIM,
    "γ ARI": MESARTHIM,
    "GAARI": MESARTHIM,
    "MESARTIM": MESARTHIM,
    "FIRST STAR OF ARIES": MESARTHIM,
    # ======== ACUBENS (Alpha Cancri) ========
    "ALPHA CANCRI": ACUBENS,
    "ALPHA CNC": ACUBENS,
    "65 CNC": ACUBENS,
    "α CNC": ACUBENS,
    "ALCNC": ACUBENS,
    "SERTAN": ACUBENS,
    "AL ZUBANAH": ACUBENS,
    # ======== TARF (Beta Cancri) ========
    "BETA CANCRI": TARF,
    "BETA CNC": TARF,
    "17 CNC": TARF,
    "β CNC": TARF,
    "BECNC": TARF,
    "AL TARF": TARF,
    # ======== ASELLUS BOREALIS (Gamma Cancri) ========
    "GAMMA CANCRI": ASELLUS_BOREALIS,
    "GAMMA CNC": ASELLUS_BOREALIS,
    "43 CNC": ASELLUS_BOREALIS,
    "γ CNC": ASELLUS_BOREALIS,
    "GACNC": ASELLUS_BOREALIS,
    "NORTHERN DONKEY": ASELLUS_BOREALIS,
    "NORTHERN ASS": ASELLUS_BOREALIS,
    # ======== ASELLUS AUSTRALIS (Delta Cancri) ========
    "DELTA CANCRI": ASELLUS_AUSTRALIS,
    "DELTA CNC": ASELLUS_AUSTRALIS,
    "47 CNC": ASELLUS_AUSTRALIS,
    "δ CNC": ASELLUS_AUSTRALIS,
    "DECNC": ASELLUS_AUSTRALIS,
    "SOUTHERN DONKEY": ASELLUS_AUSTRALIS,
    "SOUTHERN ASS": ASELLUS_AUSTRALIS,
    # ======== KAUS AUSTRALIS (Epsilon Sagittarii) ========
    "EPSILON SAGITTARII": KAUS_AUSTRALIS,
    "EPSILON SGR": KAUS_AUSTRALIS,
    "20 SGR": KAUS_AUSTRALIS,
    "ε SGR": KAUS_AUSTRALIS,
    "EPSGR": KAUS_AUSTRALIS,
    "SOUTHERN BOW": KAUS_AUSTRALIS,
    # ======== NUNKI (Sigma Sagittarii) ========
    "SIGMA SAGITTARII": NUNKI,
    "SIGMA SGR": NUNKI,
    "34 SGR": NUNKI,
    "σ SGR": NUNKI,
    "SISGR": NUNKI,
    "PELAGUS": NUNKI,
    # ======== KAUS MEDIA (Delta Sagittarii) ========
    "DELTA SAGITTARII": KAUS_MEDIA,
    "DELTA SGR": KAUS_MEDIA,
    "19 SGR": KAUS_MEDIA,
    "δ SGR": KAUS_MEDIA,
    "DESGR": KAUS_MEDIA,
    "MIDDLE BOW": KAUS_MEDIA,
    # ======== KAUS BOREALIS (Lambda Sagittarii) ========
    "LAMBDA SAGITTARII": KAUS_BOREALIS,
    "LAMBDA SGR": KAUS_BOREALIS,
    "22 SGR": KAUS_BOREALIS,
    "λ SGR": KAUS_BOREALIS,
    "LASGR": KAUS_BOREALIS,
    "NORTHERN BOW": KAUS_BOREALIS,
    # ======== ASCELLA (Zeta Sagittarii) ========
    "ZETA SAGITTARII": ASCELLA,
    "ZETA SGR": ASCELLA,
    "38 SGR": ASCELLA,
    "ζ SGR": ASCELLA,
    "ZESGR": ASCELLA,
    "ARMPIT": ASCELLA,
    # ======== ALGEDI (Alpha Capricorni) ========
    "ALPHA CAPRICORNI": ALGEDI,
    "ALPHA CAP": ALGEDI,
    "6 CAP": ALGEDI,
    "α CAP": ALGEDI,
    "ALCAP": ALGEDI,
    "GIEDI": ALGEDI,
    "PRIMA GIEDI": ALGEDI,
    "THE GOAT": ALGEDI,
    # ======== DABIH (Beta Capricorni) ========
    "BETA CAPRICORNI": DABIH,
    "BETA CAP": DABIH,
    "9 CAP": DABIH,
    "β CAP": DABIH,
    "BECAP": DABIH,
    "AL-DHABIH": DABIH,
    "LUCKY ONE OF SLAUGHTERER": DABIH,
    # ======== NASHIRA (Gamma Capricorni) ========
    "GAMMA CAPRICORNI": NASHIRA,
    "GAMMA CAP": NASHIRA,
    "40 CAP": NASHIRA,
    "γ CAP": NASHIRA,
    "GACAP": NASHIRA,
    "FORTUNATE ONE": NASHIRA,
    "SA'D NASHIRAH": NASHIRA,
    # ======== SADALSUUD (Beta Aquarii) ========
    "BETA AQUARII": SADALSUUD,
    "BETA AQR": SADALSUUD,
    "22 AQR": SADALSUUD,
    "β AQR": SADALSUUD,
    "BEAQR": SADALSUUD,
    "LUCKIEST OF LUCKY STARS": SADALSUUD,
    "SA'D AL-SU'UD": SADALSUUD,
    # ======== SADALMELIK (Alpha Aquarii) ========
    "ALPHA AQUARII": SADALMELIK,
    "ALPHA AQR": SADALMELIK,
    "34 AQR": SADALMELIK,
    "α AQR": SADALMELIK,
    "ALAQR": SADALMELIK,
    "LUCKY STAR OF KING": SADALMELIK,
    "SA'D AL-MALIK": SADALMELIK,
    # ======== SKAT (Delta Aquarii) ========
    "DELTA AQUARII": SKAT,
    "DELTA AQR": SKAT,
    "76 AQR": SKAT,
    "δ AQR": SKAT,
    "DEAQR": SKAT,
    "SCHEAT AQUARII": SKAT,
    "SHIN": SKAT,
    # ======== ETA PISCIUM ========
    "ETA PISCIUM": ETA_PISCIUM,
    "ETA PSC": ETA_PISCIUM,
    "99 PSC": ETA_PISCIUM,
    "η PSC": ETA_PISCIUM,
    "ETPSC": ETA_PISCIUM,
    "KULLAT NUNU": ETA_PISCIUM,
    # ======== ALRESCHA (Alpha Piscium) ========
    "ALPHA PISCIUM": ALRESCHA,
    "ALPHA PSC": ALRESCHA,
    "113 PSC": ALRESCHA,
    "α PSC": ALRESCHA,
    "ALPSC": ALRESCHA,
    "AL-RISHA": ALRESCHA,
    "THE KNOT": ALRESCHA,
    "THE CORD": ALRESCHA,
    # ======== ALTERNATE SPELLINGS / COMMON MISSPELLINGS ========
    # Betelgeuse variants (Arabic: yad al-jawza / ibt al-jawza)
    "BETELGEUX": BETELGEUSE,
    "BEETLEJUICE": BETELGEUSE,
    "BETELGUESE": BETELGEUSE,
    "BETELGUEUSE": BETELGEUSE,
    "BETELEGEUSE": BETELGEUSE,
    "BETEIGEUZE": BETELGEUSE,
    # Fomalhaut variants (Arabic: fum al-hut)
    "FORMALHAUT": FOMALHAUT,
    "FOMALHUT": FOMALHAUT,
    "FOMALAUT": FOMALHAUT,
    "FOMALHAULT": FOMALHAUT,
    "FUMALHAUT": FOMALHAUT,
    # Aldebaran variants (Arabic: al-dabaran)
    "ALDEBRAN": ALDEBARAN,
    "ALDEBERON": ALDEBARAN,
    "ALDEBERAN": ALDEBARAN,
    "ALDERBARAN": ALDEBARAN,
    "ALDEBARIN": ALDEBARAN,
    # Algol variants (Arabic: ra's al-ghul)
    "ALGOL": ALGOL,  # Already canonical but ensure alias exists
    "ALGHOL": ALGOL,
    "ALGOUL": ALGOL,
    "AL-GHUL": ALGOL,
    # Arcturus variants (Greek: arktos + ouros = bear guard)
    "ARCHTURUS": ARCTURUS,
    "ARTURUS": ARCTURUS,
    "ARKCTURUS": ARCTURUS,
    "ARCTUROS": ARCTURUS,
    # Antares variants (Greek: anti + Ares = rival of Mars)
    "ANTARIES": ANTARES,
    "ANTARAS": ANTARES,
    "ANTARRES": ANTARES,
    # Rigel variants (Arabic: rijl = foot)
    "RIEGEL": RIGEL,
    "RIJEL": RIGEL,
    "RIGIL": RIGEL,
    # Vega variants (Arabic: an-nasr al-waqi)
    "VEGHA": VEGA,
    # Polaris variants
    "POLARRIS": POLARIS,
    "POLARIUS": POLARIS,
    # Procyon variants (Greek: pro + kyon = before the dog)
    "PROCION": PROCYON,
    "PROCIAN": PROCYON,
    # Capella variants (Latin: she-goat)
    "CAPELA": CAPELLA,
    "CAPPELLA": CAPELLA,
    # Deneb variants (Arabic: dhanab = tail)
    "DANEB": DENEB,
    "DHENEB": DENEB,
    # Altair variants (Arabic: al-nasr al-ta'ir = flying eagle)
    "ALTIAR": ALTAIR,
    "ALTARE": ALTAIR,
    # Sirius variants (Greek: seirios = scorching)
    "SYRIUS": SIRIUS,
    "SIRUIS": SIRIUS,
    "SIRUS": SIRIUS,
    # Spica variants (Latin: ear of wheat)
    "SPIKA": SPICA_STAR,
    "SPICCA": SPICA_STAR,
    # Regulus variants (Latin: little king)
    "REGULAS": REGULUS,
    "REGULIS": REGULUS,
    # Canopus variants
    "CANOPIS": CANOPUS,
    "CANPOUS": CANOPUS,
    "CANOPOS": CANOPUS,
    # Achernar variants (Arabic: akhir an-nahr = end of river)
    "ACHENAR": ACHERNAR,
    "ARCHENAR": ACHERNAR,
    "ACHERNAN": ACHERNAR,
    # Castor/Pollux variants
    "CASTOR": CASTOR,  # Ensure canonical exists
    "KASTOR": CASTOR,
    "POLUX": POLLUX,
    "POLLUCKS": POLLUX,
    # Alpheratz variants (Alpha Andromedae)
    "ALPHERATZ": ALPHERATZ,
    "SIRRAH": ALPHERATZ,
    "SIRAH": ALPHERATZ,
    "ALPHERAT": ALPHERATZ,
    "ALPHA ANDROMEDAE": ALPHERATZ,
    # Algenib variants (Gamma Pegasi)
    "ALGENIB": ALGENIB,
    "GAMMA PEGASI": ALGENIB,
    # Propus variants (Eta Geminorum)
    "PROPUS": PROPUS,
    "ETA GEMINORUM": PROPUS,
    # Tejat variants (Mu Geminorum)
    "TEJAT": TEJAT,
    "TEJAT POSTERIOR": TEJAT,
    "MU GEMINORUM": TEJAT,
    # Alhena variants (Gamma Geminorum)
    "ALHENA": ALHENA,
    "ALMEISAN": ALHENA,
    "GAMMA GEMINORUM": ALHENA,
    # Wasat variants (Delta Geminorum)
    "WASAT": WASAT,
    "DELTA GEMINORUM": WASAT,
    # Adhara variants (Epsilon Canis Majoris)
    "ADHARA": ADHARA,
    "ADARA": ADHARA,
    "EPSILON CANIS MAJORIS": ADHARA,
    # Wezen variants (Delta Canis Majoris)
    "WEZEN": WEZEN,
    "WESEN": WEZEN,
    "DELTA CANIS MAJORIS": WEZEN,
    # Thuban variants (Alpha Draconis)
    "THUBAN": THUBAN,
    "ALPHA DRACONIS": THUBAN,
    # Rasalgethi variants (Alpha Herculis)
    "RASALGETHI": RASALGETHI,
    "RAS ALGETHI": RASALGETHI,
    "ALPHA HERCULIS": RASALGETHI,
    # Albireo variants (Beta Cygni)
    "ALBIREO": ALBIREO,
    "BETA CYGNI": ALBIREO,
    # Mirach variants (Beta Andromedae)
    "MIRACH": MIRACH,
    "BETA ANDROMEDAE": MIRACH,
    "BETA AND": MIRACH,
    "β AND": MIRACH,
    "BEAND": MIRACH,
    "MIRAK": MIRACH,
    "MIRAC": MIRACH,
    # Almach variants (Gamma Andromedae)
    "ALMACH": ALMACH,
    "ALAMAK": ALMACH,
    "ALMAAK": ALMACH,
    "GAMMA ANDROMEDAE": ALMACH,
    "GAMMA AND": ALMACH,
    "γ AND": ALMACH,
    "GA1AND": ALMACH,
    # Menkar variants (Alpha Ceti)
    "MENKAR": MENKAR,
    "MENKAB": MENKAR,
    "ALPHA CETI": MENKAR,
    "ALPHA CET": MENKAR,
    "α CET": MENKAR,
    "ALCET": MENKAR,
}


# Phonetic normalization for fuzzy matching of star common names
# Maps phonetically similar character sequences to canonical forms
_PHONETIC_REPLACEMENTS: list[tuple[str, str]] = [
    # Double consonants to single
    ("LL", "L"),
    ("TT", "T"),
    ("SS", "S"),
    ("RR", "R"),
    ("PP", "P"),
    ("CC", "C"),
    ("NN", "N"),
    ("MM", "M"),
    ("FF", "F"),
    # Common vowel confusions
    ("EU", "U"),  # Betelgeuse -> Betelguse
    ("EI", "I"),  # Beteigeuze -> Betigeuze
    ("AU", "A"),  # Fomalhaut -> Fomalhat
    ("AE", "E"),
    ("OE", "E"),
    ("OU", "O"),
    # Silent/ambiguous endings
    ("UE$", "U"),  # Betelgeuse -> Betelgeus
    ("E$", ""),  # trailing silent e
    # Common consonant confusions
    ("GH", "G"),  # Alghol -> Algol
    ("PH", "F"),
    ("CK", "K"),
    ("CH", "K"),
    ("SCH", "SH"),
    # X sounds
    ("X", "KS"),
    # Vowel simplification
    ("EE", "E"),
    ("AA", "A"),
    ("OO", "O"),
    ("II", "I"),
    ("UU", "U"),
    # Remove isolated vowels between consonants for approximate matching
    ("(?<=[BCDFGHJKLMNPQRSTVWXYZ])I(?=[BCDFGHJKLMNPQRSTVWXYZ])", ""),
]


def _normalize_phonetic(name: str) -> str:
    """
    Normalize a star name for phonetic fuzzy matching.

    Applies a series of transformations to reduce phonetically similar
    names to the same normalized form. This allows matching alternate
    spellings like:
    - Betelgeuse / Betelgeux / Beetlejuice
    - Fomalhaut / Formalhaut / Fomalaut
    - Aldebaran / Aldebran / Aldeberan

    Args:
        name: Star name to normalize

    Returns:
        Phonetically normalized name (uppercase, consonants preserved)

    Examples:
        >>> _normalize_phonetic("Betelgeuse")
        'BETLGUS'
        >>> _normalize_phonetic("Betelgeux")
        'BETLGUS'
        >>> _normalize_phonetic("Fomalhaut")
        'FOMALHAT'
    """
    import re

    result = name.upper().strip()

    # Remove non-alphabetic characters
    result = re.sub(r"[^A-Z]", "", result)

    # Apply phonetic replacements (order matters for some)
    for pattern, replacement in _PHONETIC_REPLACEMENTS:
        if pattern.endswith("$"):
            # End-of-string pattern
            actual_pattern = pattern[:-1]
            if result.endswith(actual_pattern):
                result = result[: -len(actual_pattern)] + replacement
        elif pattern.startswith("(?"):
            # Regex pattern
            result = re.sub(pattern, replacement, result)
        else:
            result = result.replace(pattern, replacement)

    # Remove vowels except leading vowel (keeps consonant skeleton)
    if len(result) > 1:
        first_char = result[0]
        rest = result[1:]
        rest = re.sub(r"[AEIOU]", "", rest)
        result = first_char + rest

    return result


def _fuzzy_match_star(name: str) -> int | None:
    """
    Find a star using fuzzy phonetic matching.

    This function normalizes the input name and compares it against
    normalized versions of all star aliases and canonical names.
    It's used as a fallback when exact matching fails.

    Args:
        name: Star name to search for (may be misspelled)

    Returns:
        Star ID if a unique match is found, None otherwise

    Examples:
        >>> _fuzzy_match_star("Betelgeux")  # Alternate spelling
        BETELGEUSE
        >>> _fuzzy_match_star("Formalhaut")  # Common misspelling
        FOMALHAUT
    """
    normalized_input = _normalize_phonetic(name)

    if len(normalized_input) < 3:
        # Too short for reliable fuzzy matching
        return None

    matches: list[int] = []
    matched_names: list[str] = []

    # Check against canonical names in STAR_CATALOG
    for entry in STAR_CATALOG:
        if _normalize_phonetic(entry.name) == normalized_input:
            if entry.id not in matches:
                matches.append(entry.id)
                matched_names.append(entry.name)

    # Check against aliases
    for alias, star_id in STAR_ALIASES.items():
        if _normalize_phonetic(alias) == normalized_input:
            if star_id not in matches:
                matches.append(star_id)
                matched_names.append(alias)

    if len(matches) == 1:
        return matches[0]
    elif len(matches) > 1:
        # Ambiguous match - return None to avoid false positives
        return None

    return None


def resolve_star_name(name: str) -> int | None:
    """
    Resolve a star name to its SE_* constant ID using reference API-compatible resolution.

    Implements the following resolution algorithm:
    1. Normalize input (uppercase, strip whitespace)
    2. If name starts with comma (reference API convention), do prefix search
    3. Try exact match in STAR_ALIASES dictionary
    4. Try exact match against canonical star names in STAR_CATALOG
    5. Try Bayer designation with Greek letter names (e.g., "Alpha Leonis")
    6. Try Flamsteed designation with number + constellation (e.g., "32 Leonis")
    7. Try fuzzy matching (alias contains search term for short inputs)
    8. Try phonetic fuzzy matching for common misspellings

    Args:
        name: Star name, alias, designation, or comma-prefixed partial name

    Returns:
        Star ID (SE_* constant) if found, None otherwise

    Examples:
        >>> resolve_star_name("Regulus")
        1000001
        >>> resolve_star_name(",alg")  # Prefix search for Algol
        1000003
        >>> resolve_star_name("Alpha Leo")
        1000001
        >>> resolve_star_name("Alpha Leonis")  # Full Bayer designation
        1000001
        >>> resolve_star_name("32 Leonis")  # Flamsteed designation
        1000001
        >>> resolve_star_name("SIRIUS")
        1000004
        >>> resolve_star_name("Betelgeux")  # Alternate spelling
        BETELGEUSE
        >>> resolve_star_name("Formalhaut")  # Common misspelling
        FOMALHAUT
    """
    if not name:
        return None

    # Normalize: uppercase and strip
    normalized = name.upper().strip()

    # Handle comma-prefix for partial match (reference API convention)
    if normalized.startswith(","):
        search_prefix = normalized[1:].strip()
        if not search_prefix:
            return None

        # Search for any alias or canonical name that starts with this prefix
        # First check canonical names
        for entry in STAR_CATALOG:
            if entry.name.upper().startswith(search_prefix):
                return entry.id

        # Then check aliases
        for alias, star_id in STAR_ALIASES.items():
            if alias.startswith(search_prefix):
                return star_id

        # Check nomenclature (e.g., alLeo, bePer)
        for entry in STAR_CATALOG:
            if entry.nomenclature.upper().startswith(search_prefix):
                return entry.id

        return None

    # Strip trailing comma and anything after (reference API format: "Regulus,alLeo")
    if "," in normalized:
        normalized = normalized.split(",")[0].strip()

    # 1. Try exact match against canonical star names — a star's own
    # catalog name outranks a legacy alias pointing elsewhere (e.g.
    # 'Suhail' is the IAU name of lambda Velorum but also a historical
    # alias of Canopus).
    if normalized in _STAR_NAME_TO_ID:
        return _STAR_NAME_TO_ID[normalized]

    # 2. Try exact match in STAR_ALIASES
    if normalized in STAR_ALIASES:
        return STAR_ALIASES[normalized]

    # 3. Try exact match against nomenclature (e.g., "ALLEO", "BEPER")
    for entry in STAR_CATALOG:
        if entry.nomenclature.upper() == normalized:
            return entry.id

    # 4. Try Bayer designation with Greek letter names (e.g., "Alpha Leonis")
    parsed_nomenclature = _parse_bayer_designation(name.strip())
    if parsed_nomenclature:
        for entry in STAR_CATALOG:
            if entry.nomenclature.upper() == parsed_nomenclature.upper():
                return entry.id

    # 5. Try Flamsteed designation (e.g., "32 Leonis", "87 Virginis")
    parsed_flamsteed = _parse_flamsteed_designation(name.strip())
    if parsed_flamsteed:
        # Look up in STAR_ALIASES using the normalized format (e.g., "32 LEO")
        if parsed_flamsteed in STAR_ALIASES:
            return STAR_ALIASES[parsed_flamsteed]

    # 6. Try prefix matching: check if any alias STARTS WITH the search term
    # Uses prefix matching instead of substring to avoid false positives
    # (e.g., "GENIB" matching "ZUBENELGENUBI" via substring)
    if len(normalized) >= 3:
        for alias, star_id in STAR_ALIASES.items():
            if alias.startswith(normalized):
                return star_id
        for entry in STAR_CATALOG:
            if entry.name.upper().startswith(normalized):
                return entry.id

    # 7. Try phonetic fuzzy matching for common misspellings
    fuzzy_result = _fuzzy_match_star(name)
    if fuzzy_result is not None:
        return fuzzy_result

    return None


def get_canonical_star_name(star_id: int) -> str | None:
    """
    Get the canonical star name for a star ID.

    Args:
        star_id: Star ID (SE_* constant)

    Returns:
        Canonical star name (e.g., "Regulus") or None if not found
    """
    for entry in STAR_CATALOG:
        if entry.id == star_id:
            return entry.name
    return None


def _calc_star_position_from_observer(
    star_id: int,
    earth_at_t,
    noaberr: bool = False,
    nogdefl: bool = False,
    j2000_frame: bool = False,
) -> Tuple[float, float, float]:
    """Calculate one fixed star from a precomputed observer position."""
    if star_id not in FIXED_STARS:
        raise ValueError(f"could not find star name {star_id}")

    star_data = FIXED_STARS[star_id]

    # Create Skyfield Star object
    # Note: pm_ra in StarData is already "pm_ra * cos(dec)" (proper motion in RA direction)
    # Skyfield expects ra_mas_per_year which is the same convention
    star = Star(
        ra_hours=star_data.ra_j2000 / 15.0,  # Convert degrees to hours
        dec_degrees=star_data.dec_j2000,
        ra_mas_per_year=star_data.pm_ra * 1000.0,  # arcsec/yr to mas/yr
        dec_mas_per_year=star_data.pm_dec * 1000.0,
        # Stars without a measured parallax get the reference default
        # (0.0001249 arcsec = 0.1249 mas).
        parallax_mas=star_data.parallax_mas if star_data.parallax_mas > 0.0 else 0.1249,
        radial_km_per_s=star_data.radial_km_per_s,  # Radial velocity for distance change
    )

    # Calculate position
    astrometric = earth_at_t.observe(star)

    if noaberr:
        pos = astrometric
    elif nogdefl:
        # Aberration without gravitational deflection:
        # Pass empty deflectors tuple to skip deflection while keeping aberration.
        pos = astrometric.apparent(deflectors=())
    else:
        pos = astrometric.apparent()

    # Transform to ecliptic coordinates
    frame = ecliptic_J2000_frame if j2000_frame else ecliptic_frame
    ecl = pos.frame_latlon(frame)

    # ecl returns (latitude, longitude, distance) as Skyfield Angle/Distance objects
    lat = ecl[0].degrees
    lon = ecl[1].degrees % 360.0
    # Distance from parallax (AU). Stars without measured parallax use
    # the reference's default of 0.0001249 arcsec (the Star object below
    # already received it), so the distance stays finite and consistent.
    dist = ecl[2].au

    return lon, lat, dist


def _calc_star_position_leb(
    star_id: int,
    jd_tt: float,
    noaberr: bool = False,
    nogdefl: bool = False,
    j2000_frame: bool = False,
) -> Tuple[float, float, float]:
    """Calculate fixed star position using LEB data (no Skyfield DE kernel).

    Same interface as _calc_star_position_skyfield but uses the LEB reader
    for Earth position/velocity instead of loading the full DE kernel.

    The pipeline mirrors _pipeline_icrs() in fast_calc.py adapted for stars:
    Earth from LEB → proper motion → ICRS Cartesian → geocentric vector →
    light-time → gravitational deflection → aberration → frame rotation.
    """
    import math

    from .fast_calc import (
        C_LIGHT_AU_DAY,
        _apply_aberration,
        _apply_gravitational_deflection,
        _cartesian_to_spherical,
        _fw2m,
        _get_leb_frame_data,
        _iau2006_precession_angles,
        _mat3_vec3,
        _rotate_equatorial_to_ecliptic,
        _rotate_icrs_to_ecliptic_j2000,
        _vec3_dist,
        _vec3_sub,
    )
    from .state import get_leb_reader

    reader = get_leb_reader()
    if reader is None:
        raise KeyError("No LEB reader available")

    if star_id not in FIXED_STARS:
        raise ValueError(f"could not find star name {star_id}")

    star_data = FIXED_STARS[star_id]
    from .constants import EARTH

    # 1. Earth position and velocity from LEB (ICRS barycentric)
    earth_pos, earth_vel = reader.eval_body(EARTH, jd_tt)

    # 2. Propagate proper motion from J2000 to observation date
    ra_date, dec_date = propagate_proper_motion(
        star_data.ra_j2000,
        star_data.dec_j2000,
        star_data.pm_ra,
        star_data.pm_dec,
        J2000,
        jd_tt,
    )

    # 3. RA/Dec → ICRS Cartesian
    if star_data.parallax_mas > 0.0:
        dist_internal = 206265000.0 / star_data.parallax_mas
    else:
        # Reference default parallax 0.0001249 arcsec for stars without
        # a measured value.
        dist_internal = 206265000.0 / 0.1249

    # Apply radial velocity correction to distance (matching Skyfield Star behavior)
    if star_data.radial_km_per_s != 0.0 and dist_internal < 1e11:
        dt_years = (jd_tt - J2000) / DAYS_PER_JULIAN_YEAR
        au_per_year = star_data.radial_km_per_s * 365.25 * 86400.0 / 149597870.7
        dist_internal += au_per_year * dt_years

    ra_rad = math.radians(ra_date)
    dec_rad = math.radians(dec_date)
    cos_dec = math.cos(dec_rad)
    star_icrs = (
        dist_internal * cos_dec * math.cos(ra_rad),
        dist_internal * cos_dec * math.sin(ra_rad),
        dist_internal * math.sin(dec_rad),
    )

    # 4. Geocentric vector
    geo = _vec3_sub(star_icrs, earth_pos)

    # 5. Light-time (zero for stars without parallax — infinite direction)
    # Note: for finite-distance stars, proper motion is not re-evaluated
    # at the retarded epoch (jd_tt - lt).  The error is < 0.001" for all
    # catalog stars because lt is at most ~few years and pm is < 10"/yr.
    if star_data.parallax_mas > 0.0:
        geo_dist = _vec3_dist(geo)
        lt = geo_dist / C_LIGHT_AU_DAY if geo_dist > 0 else 0.0
    else:
        lt = 0.0

    # 6. Gravitational deflection (skip if noaberr or nogdefl)
    if not noaberr and not nogdefl:
        geo = _apply_gravitational_deflection(geo, earth_pos, jd_tt, lt, reader)

    # 7. Aberration (skip if noaberr)
    if not noaberr:
        geo = _apply_aberration(geo, earth_vel, lt)

    # 8. Frame rotation
    if j2000_frame:
        ecl = _rotate_icrs_to_ecliptic_j2000(geo[0], geo[1], geo[2])
        lon, lat, dist = _cartesian_to_spherical(ecl[0], ecl[1], ecl[2])
    else:
        try:
            pn_mat, _dpsi, _deps, eps_true_rad = _get_leb_frame_data(reader, jd_tt)
        except (KeyError, ValueError, AttributeError):
            import erfa

            dpsi_rad, deps_rad = erfa.nut06a(2451545.0, jd_tt - 2451545.0)
            gamb, phib, psib, epsa = _iau2006_precession_angles(jd_tt)
            eps_true_rad = epsa + deps_rad
            pn_mat = _fw2m(gamb, phib, psib + dpsi_rad, eps_true_rad)
        geo_eq = _mat3_vec3(pn_mat, geo)
        ecl = _rotate_equatorial_to_ecliptic(
            geo_eq[0], geo_eq[1], geo_eq[2], eps_true_rad
        )
        lon, lat, dist = _cartesian_to_spherical(ecl[0], ecl[1], ecl[2])

    if star_data.parallax_mas == 0.0:
        dist = 100000.0

    return lon, lat, dist


def _calc_star_position_skyfield(
    star_id: int,
    jd_tt: float,
    noaberr: bool = False,
    nogdefl: bool = False,
    j2000_frame: bool = False,
    topo: "tuple | None" = None,
) -> Tuple[float, float, float]:
    """
    Calculate ecliptic position using Skyfield Star class with proper aberration.

    This function uses Skyfield's Star class which correctly handles:
    - Proper motion
    - Precession
    - Nutation
    - Annual aberration (unless noaberr=True)
    - Gravitational light deflection (unless nogdefl=True)

    Args:
        star_id: Star identifier (SE_* constant)
        jd_tt: Julian Day in Terrestrial Time (TT)
        noaberr: If True, skip aberration correction (astrometric position)
        nogdefl: If True, skip gravitational deflection but keep aberration
        j2000_frame: If True, return J2000 ecliptic coordinates instead of
            ecliptic of date.  This avoids the ~5" error from precessing
            Skyfield's ecliptic-of-date back to J2000 with a different
            precession model.

    Returns:
        Tuple of (longitude, latitude, distance) in degrees and AU

    Raises:
        ValueError: If star_id not in catalog
    """
    from .state import get_planets, get_timescale

    t = get_timescale().tt_jd(jd_tt)
    earth = get_planets()["earth"]
    if topo is not None:
        from skyfield.api import wgs84

        observer = earth + wgs84.latlon(topo[1], topo[0], topo[2])
        earth_at_t = observer.at(t)
    else:
        earth_at_t = earth.at(t)
    return _calc_star_position_from_observer(
        star_id, earth_at_t, noaberr, nogdefl, j2000_frame
    )


def calc_fixed_star_position(
    star_id: int,
    jd_tt: float,
    noaberr: bool = False,
    nogdefl: bool = False,
    j2000_frame: bool = False,
    topo: "tuple | None" = None,
) -> Tuple[float, float, float]:
    """
    Calculate ecliptic position of a fixed star at given date.

    Uses Skyfield Star class for accurate position calculation including:
    - Proper motion
    - Precession and nutation
    - Annual aberration (by default)
    - Gravitational light deflection (by default)

    Args:
        star_id: Star identifier (REGULUS, SPICA_STAR, etc.)
        jd_tt: Julian Day in Terrestrial Time (TT)
        noaberr: If True, skip aberration correction (astrometric position)
        nogdefl: If True, skip gravitational deflection but keep aberration
        j2000_frame: If True, return J2000 ecliptic coordinates

    Returns:
        Tuple[float, float, float]: (longitude, latitude, distance) where:
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees
            - distance: Arbitrary large value (AU) - stars are effectively infinite

    Raises:
        ValueError: If star_id not in catalog

    Notes:
        By default, includes annual aberration as is standard practice.
        Use noaberr=True for astrometric (geometric) position.

    References:
        IAU 2006 precession (Capitaine et al.)
        Skyfield library for apparent position calculation
    """
    from .state import get_leb_reader

    # The LEB star path is geocentric; topocentric requests go through
    # Skyfield (same convention as the planet pipeline).
    if topo is None and get_leb_reader() is not None:
        try:
            return _calc_star_position_leb(
                star_id, jd_tt, noaberr, nogdefl, j2000_frame
            )
        except KeyError:
            pass  # Body (EARTH) not in LEB file
        except ValueError as _leb_err:
            if "outside range" not in str(_leb_err).lower():
                raise  # Re-raise unexpected ValueError
            from .logging_config import get_logger

            get_logger().debug("LEB star fallback: %s", _leb_err)
    return _calc_star_position_skyfield(
        star_id, jd_tt, noaberr, nogdefl, j2000_frame, topo=topo
    )


def calc_fixed_star_velocity(
    star_id: int,
    jd_tt: float,
    noaberr: bool = False,
    nogdefl: bool = False,
    j2000_frame: bool = False,
    topo: "tuple | None" = None,
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate ecliptic position and velocity of a fixed star at given date.

    Uses central finite difference for velocity calculation.

    The velocity represents the rate of change of the star's ecliptic
    coordinates due to:
    1. Proper motion of the star itself through space
    2. Precession of the equinoxes (ecliptic coordinate grid rotation)
    3. Nutation (small periodic oscillations)
    4. Annual aberration (unless noaberr=True)

    The dominant contribution is precession at ~50.3 arcseconds/year = 0.0001378
    degrees/day in longitude.

    Args:
        star_id: Star identifier (REGULUS, SPICA_STAR, etc.)
        jd_tt: Julian Day in Terrestrial Time (TT)
        noaberr: If True, skip aberration correction (astrometric position)
        nogdefl: If True, skip gravitational deflection but keep aberration
        j2000_frame: If True, compute in J2000 ecliptic frame

    Returns:
        Tuple[float, float, float, float, float, float]:
            (longitude, latitude, distance, speed_lon, speed_lat, speed_dist)

    Raises:
        ValueError: If star_id not in catalog
    """
    if star_id not in FIXED_STARS:
        raise ValueError(f"could not find star name {star_id}")

    # Half-day step for central difference
    h = 0.5

    # Calculate position at current time (for return value)
    lon, lat, dist = calc_fixed_star_position(
        star_id, jd_tt, noaberr, nogdefl, j2000_frame, topo=topo
    )

    # Calculate positions at t-0.5 and t+0.5 for central difference
    lon_prev, lat_prev, dist_prev = calc_fixed_star_position(
        star_id, jd_tt - h, noaberr, nogdefl, j2000_frame, topo=topo
    )
    lon_next, lat_next, dist_next = calc_fixed_star_position(
        star_id, jd_tt + h, noaberr, nogdefl, j2000_frame, topo=topo
    )

    # Central difference: (f(t+h) - f(t-h)) / (2h) where 2h = 1.0 day
    speed_lon = lon_next - lon_prev

    # Handle wraparound at 360° (e.g., 359° -> 1° should give +2°, not -358°)
    if speed_lon > 180.0:
        speed_lon -= 360.0
    elif speed_lon < -180.0:
        speed_lon += 360.0

    # Latitude speed: pure finite difference (no wraparound needed for latitude)
    speed_lat = lat_next - lat_prev

    # Distance speed: radial motion (the star's radial velocity plus the
    # Earth's orbital radial component) - the reference reports it.
    speed_dist = dist_next - dist_prev

    return lon, lat, dist, speed_lon, speed_lat, speed_dist


def _se_star_key(name: str) -> str:
    """Reference search key for a traditional star name.

    The reference removes every whitespace character and lowercases the
    traditional-name part; the Bayer/Flamsteed part after a comma keeps
    its case.
    """
    s = "".join(ch for ch in name if not ch.isspace())
    head, sep, tail = s.partition(",")
    return head.lower() + sep + tail


_SE_SORTED_CATALOG: "list | None" = None


def _se_sorted_catalog():
    """Catalog entries sorted by their reference search key.

    Sequential star numbers ("1", "2", ...) index this order, mirroring
    the reference's sorted in-memory star list (its v2 behavior; the
    sequence is specific to the catalog shipped with this library).
    """
    global _SE_SORTED_CATALOG
    if _SE_SORTED_CATALOG is None:
        _SE_SORTED_CATALOG = sorted(STAR_CATALOG, key=lambda e: _se_star_key(e.name))
    return _SE_SORTED_CATALOG


def _resolve_star_se(star_name: str) -> tuple[int, str | None, str | None]:
    """Resolve a star with the reference's exact search semantics.

    Mirrors the reference resolution rules:
    - leading comma: the rest is a Bayer/Flamsteed nomenclature key,
      matched exactly (case preserved): ",alTau" -> Aldebaran;
    - "name,nomenclature": the traditional-name part is ignored and the
      nomenclature key decides;
    - leading digit: 1-based sequential number in the sorted catalog;
    - trailing '%' on a traditional name: prefix wildcard (a '%'
      anywhere else is an error);
    - otherwise: exact traditional-name match (whitespace removed,
      case-insensitive).

    No fuzzy matching, Bayer-word parsing or prefix guessing happens
    here - those remain available in the library's own search helpers
    (resolve_star_name, search helpers), not in the reference-named
    functions.

    Returns:
        (star_id, error_message, "Name,nomenclature"); on error the id
        is -1 and the canonical name is None.
    """
    sstar = _se_star_key(star_name)
    if not sstar:
        return -1, "star name empty", None

    def _found(entry) -> tuple[int, None, str]:
        return entry.id, None, f"{entry.name},{entry.nomenclature}"

    # ",nomenclature" or "name,nomenclature": nomenclature key search.
    if "," in sstar:
        key = sstar[sstar.index(",") + 1 :]
        if not key:
            return -1, f"could not find star name {sstar}", None
        for entry in STAR_CATALOG:
            if "".join(entry.nomenclature.split()) == key:
                return _found(entry)
        return -1, f"could not find star name ,{key}", None

    # Sequential star number.
    if sstar[0].isdigit():
        digits = ""
        for ch in sstar:
            if ch.isdigit():
                digits += ch
            else:
                break
        star_nr = int(digits)
        ordered = _se_sorted_catalog()
        if star_nr < 1 or star_nr > len(ordered):
            return (
                -1,
                f"sequential fixed star number {star_nr} is not available",
                None,
            )
        return _found(ordered[star_nr - 1])

    # Trailing '%' wildcard on the traditional name.
    if "%" in sstar:
        if not sstar.endswith("%") or sstar.count("%") != 1:
            return -1, f"invalid search string {sstar}", None
        prefix = sstar[:-1]
        for entry in _se_sorted_catalog():
            if _se_star_key(entry.name).startswith(prefix):
                return _found(entry)
        return -1, f"star search string {sstar} did not match", None

    # Exact traditional-name match. The alias table plays the role of
    # the reference catalog's additional name lines per star, so exact
    # alias keys resolve too (still no fuzzy or prefix matching).
    for entry in STAR_CATALOG:
        if _se_star_key(entry.name) == sstar:
            return _found(entry)
    for alias, star_id in STAR_ALIASES.items():
        if _se_star_key(alias).lower() == sstar:
            for entry in STAR_CATALOG:
                if entry.id == star_id:
                    return _found(entry)
    return -1, f"could not find star name {sstar}", None


def _resolve_star_id(star_name: str) -> tuple[int, str | None, str | None]:
    """
    Resolve a star name to its ID with reference API-compatible name resolution.

    Supports:
    - Exact star names: "Regulus", "Spica"
    - Bayer designations: "Alpha Leo", "Alpha Leonis"
    - Flamsteed numbers: "87 Leo"
    - Traditional/Arabic names: "Cor Leonis", "Dog Star"
    - Comma-prefix partial match: ",alg" finds Algol

    Args:
        star_name: Name of star (e.g. "Regulus", "Alpha Leo", ",alg")

    Returns:
        Tuple of (star_id, error_message, canonical_name).
        If error, star_id is -1 and canonical_name is None.
    """
    return _resolve_star_se(star_name)


def _preprocess_flags(iflag: int) -> int:
    """Preprocess calculation flags for fixed star functions.

    Strips ephemeris selection flags (MOSEPH) and converts SPEED3 to SPEED.
    FLG_TOPOCTR is kept: the observer position contributes diurnal
    aberration (and, for the nearest stars, a tiny diurnal parallax).

    Args:
        iflag: Raw input flags

    Returns:
        Cleaned flags
    """
    # Strip FLG_MOSEPH — accepted for compatibility, always uses Skyfield
    iflag = iflag & ~FLG_MOSEPH
    # FLG_SPEED3: treat as FLG_SPEED
    if iflag & FLG_SPEED3:
        iflag = (iflag & ~FLG_SPEED3) | FLG_SPEED
    return iflag


def _fixstar_ret_flags(flags_in: int) -> int:
    """Return flags echoed to the caller (reference convention).

    The input flags come back verbatim, with FLG_SWIEPH added when no
    ephemeris-selection bit was given (MOSEPH echoes as given).
    """
    from .constants import FLG_JPLEPH

    if not (flags_in & (FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH)):
        return flags_in | FLG_SWIEPH
    return flags_in


def _fixstar_topo() -> tuple:
    """(lon, lat, alt) for FLG_TOPOCTR star calculations.

    Raises Error when no geographic position has been set, like the
    planet path does.
    """
    from .state import get_topo

    topo = get_topo()
    if topo is None:
        raise Error(
            "topocentric position requested (FLG_TOPOCTR) but no "
            "geographic position set; call set_topo() first"
        )
    return (
        topo.longitude.degrees,
        topo.latitude.degrees,
        topo.elevation.m,
    )


def _apply_fixstar_flags(
    result: tuple, jd_tt: float, iflag: int, j2000_native: bool = False
) -> tuple:
    """Apply post-calculation flag transformations to fixed star results.

    Handles frame transformations (J2000, NONUT, ICRS), coordinate system
    changes (EQUATORIAL, SIDEREAL), and output format conversions (XYZ, RADIANS).

    Applied in order:
    1. J2000 / NONUT frame selection (ecliptic longitude adjustment)
    2. EQUATORIAL coordinate transformation
    3. SIDEREAL ayanamsha subtraction (from lon or RA)
    4. XYZ / RADIANS output format conversion

    Note: For fixed stars, SE applies sidereal correction AFTER equatorial
    conversion — subtracting ayanamsha from RA. This differs from planets
    where SE ignores sidereal for equatorial output entirely.

    Args:
        result: 6-tuple (lon, lat, dist, speed_lon, speed_lat, speed_dist)
        jd_tt: Julian Day in Terrestrial Time
        iflag: Calculation flags
        j2000_native: If True, positions are already in J2000 ecliptic frame
            (computed natively by Skyfield), so skip backward precession.

    Returns:
        Transformed 6-tuple
    """
    import math

    lon, lat, dist, speed_lon, speed_lat, speed_dist = result

    is_equatorial = bool(iflag & FLG_EQUATORIAL)

    # ---- 1. Frame selection ----
    # J2000: precess from ecliptic-of-date back to J2000 ecliptic
    # NONUT: remove nutation from ecliptic longitude (mean ecliptic of date)
    # These apply regardless of equatorial output (adjust ecliptic first).
    if iflag & FLG_J2000:
        if not j2000_native:
            from .astrometry import _precess_ecliptic

            lon, lat = _precess_ecliptic(lon, lat, jd_tt, J2000)
        # else: already in J2000 frame from Skyfield's ecliptic_J2000_frame
    elif iflag & FLG_NONUT:
        # Skyfield returns positions on the true ecliptic of date (with nutation).
        # NONUT means output on the mean ecliptic of date, so subtract dpsi
        # (nutation in longitude) from the ecliptic longitude.
        from .cache import get_cached_nutation

        dpsi_rad, _ = get_cached_nutation(jd_tt)
        lon = (lon - math.degrees(dpsi_rad)) % 360.0

    # ---- 2. Equatorial coordinate transformation ----
    if is_equatorial:
        if iflag & FLG_J2000:
            # J2000 obliquity for J2000 equatorial frame
            eps = 23.4392911  # IAU 2006 mean obliquity at J2000.0
        elif iflag & FLG_NONUT:
            # Mean equator of date: use mean obliquity (no nutation)
            eps = get_mean_obliquity(jd_tt)
        else:
            # True equator of date: use true obliquity
            eps = get_true_obliquity(jd_tt)

        result_sp = cotrans_sp((lon, lat, dist, speed_lon, speed_lat, speed_dist), -eps)
        lon, lat, dist = result_sp[0], result_sp[1], result_sp[2]
        speed_lon, speed_lat, speed_dist = result_sp[3], result_sp[4], result_sp[5]

    # ---- 3. Sidereal mode (ayanamsha subtraction) ----
    # For fixed stars, SE subtracts ayanamsha from the first coordinate
    # (ecliptic longitude or RA) AFTER equatorial conversion. This differs
    # from planets where SE ignores sidereal for equatorial output entirely.
    if iflag & FLG_SIDEREAL:
        from .state import get_timescale

        ts = get_timescale()
        t = ts.tt_jd(jd_tt)
        tjd_ut = t.ut1

        from .planets import get_ayanamsa_ut

        ayanamsa = get_ayanamsa_ut(tjd_ut)
        lon = (lon - ayanamsa) % 360.0

    # ---- 4. Output format conversion ----
    result = (lon, lat, dist, speed_lon, speed_lat, speed_dist)

    if iflag & FLG_XYZ:
        lon_rad = math.radians(lon)
        lat_rad = math.radians(lat)
        cos_lat = math.cos(lat_rad)
        sin_lat = math.sin(lat_rad)
        cos_lon = math.cos(lon_rad)
        sin_lon = math.sin(lon_rad)

        x = dist * cos_lat * cos_lon
        y = dist * cos_lat * sin_lon
        z = dist * sin_lat

        dlon_rad = math.radians(speed_lon)
        dlat_rad = math.radians(speed_lat)

        vx = (
            speed_dist * cos_lat * cos_lon
            - dist * sin_lat * cos_lon * dlat_rad
            - dist * cos_lat * sin_lon * dlon_rad
        )
        vy = (
            speed_dist * cos_lat * sin_lon
            - dist * sin_lat * sin_lon * dlat_rad
            + dist * cos_lat * cos_lon * dlon_rad
        )
        vz = speed_dist * sin_lat + dist * cos_lat * dlat_rad

        return (float(x), float(y), float(z), float(vx), float(vy), float(vz))

    if iflag & FLG_RADIANS:
        return (
            math.radians(lon),
            math.radians(lat),
            dist,
            math.radians(speed_lon),
            math.radians(speed_lat),
            speed_dist,
        )

    return result


def fixstar_ut(
    star: str, tjdut: float, flags: int = FLG_SWIEPH
) -> Tuple[Tuple[float, float, float, float, float, float], str, int]:
    """
    Calculate position of a fixed star for Universal Time.

    Reference API compatible function.

    Args:
        star: Name of star (e.g. "Regulus", "Spica")
        tjdut: Julian Day in Universal Time (UT1)
        flags: Calculation flags

    Returns:
        Tuple containing:
            - Position tuple: (lon, lat, dist, speed_lon, speed_lat, speed_dist)
            - star_name: Resolved star name (e.g. "Regulus")
            - retflags: Return flags (int)

    Raises:
        Error: If the star cannot be found.

    Note:
        UT (Universal Time) is converted to TT (Terrestrial Time) internally
        using Delta T before calculating the star position. For most modern dates,
        Delta T is about 69 seconds (as of 2020).

    Example:
        >>> pos, name, retflag = fixstar_ut("Regulus", 2451545.0, 0)
        >>> lon, lat, dist = pos[0], pos[1], pos[2]
    """
    ret_flags = _fixstar_ret_flags(flags)
    flags = _preprocess_flags(flags)

    star_id, error, canonical_name = _resolve_star_id(star)
    if error:
        raise Error(error)

    # Convert UT to TT using timescale (applies Delta T)
    from .state import get_timescale

    t = get_timescale().ut1_jd(tjdut)

    try:
        noaberr = bool(flags & FLG_NOABERR) or bool(flags & FLG_TRUEPOS)
        nogdefl = bool(flags & FLG_NOGDEFL)
        # Compute natively in J2000 ecliptic frame when requested.
        # This avoids the ~5" error from precessing Skyfield's ecliptic-of-date
        # back to J2000 with a different precession model.
        use_j2000 = bool(flags & FLG_J2000)
        topo = _fixstar_topo() if flags & FLG_TOPOCTR else None

        if flags & FLG_SPEED:
            lon, lat, dist, speed_lon, speed_lat, speed_dist = calc_fixed_star_velocity(
                star_id, t.tt, noaberr, nogdefl, j2000_frame=use_j2000, topo=topo
            )
            result = (lon, lat, dist, speed_lon, speed_lat, speed_dist)
        else:
            lon, lat, dist = calc_fixed_star_position(
                star_id, t.tt, noaberr, nogdefl, j2000_frame=use_j2000, topo=topo
            )
            result = (lon, lat, dist, 0.0, 0.0, 0.0)

        result = _apply_fixstar_flags(result, t.tt, flags, j2000_native=use_j2000)

        return (result, canonical_name or "", ret_flags)
    except Error:
        raise
    except (OSError, ValueError, KeyError) as e:
        raise Error(str(e)) from e


def batch_fixstars_ut(
    stars: Sequence[str],
    tjdut: float,
    flags: int = FLG_SWIEPH,
    *,
    skip_errors: bool = False,
) -> Tuple[
    Tuple[Tuple[float, float, float, float, float, float], str, int] | None, ...
]:
    """Calculate multiple fixed-star positions for Universal Time.

    The result order matches the input order. When ``skip_errors`` is True,
    unresolved stars keep their input slot as ``None``.
    """
    ret_flags = _fixstar_ret_flags(flags)
    flags = _preprocess_flags(flags)

    results: list[
        Tuple[Tuple[float, float, float, float, float, float], str, int] | None
    ] = [None] * len(stars)
    resolved: list[tuple[int, int, str]] = []
    for index, star_name in enumerate(stars):
        star_id, error, canonical_name = _resolve_star_id(star_name)
        if error:
            if skip_errors:
                continue
            raise Error(error)
        resolved.append((index, star_id, canonical_name or ""))

    if not resolved:
        return tuple(results)

    from .state import get_leb_reader, get_timescale

    noaberr = bool(flags & FLG_NOABERR) or bool(flags & FLG_TRUEPOS)
    nogdefl = bool(flags & FLG_NOGDEFL)
    use_j2000 = bool(flags & FLG_J2000)
    want_speed = bool(flags & FLG_SPEED)

    ts = get_timescale()
    t = ts.ut1_jd(tjdut)
    jd_tt = t.tt

    reader = get_leb_reader()
    _leb_ok = False
    if reader is not None:
        try:
            for index, star_id, canonical_name in resolved:
                lon, lat, dist = _calc_star_position_leb(
                    star_id, jd_tt, noaberr, nogdefl, use_j2000
                )
                if want_speed:
                    lon_prev, lat_prev, dist_prev = _calc_star_position_leb(
                        star_id, jd_tt - 0.5, noaberr, nogdefl, use_j2000
                    )
                    lon_next, lat_next, dist_next = _calc_star_position_leb(
                        star_id, jd_tt + 0.5, noaberr, nogdefl, use_j2000
                    )
                    speed_lon = lon_next - lon_prev
                    if speed_lon > 180.0:
                        speed_lon -= 360.0
                    elif speed_lon < -180.0:
                        speed_lon += 360.0
                    speed_lat = lat_next - lat_prev
                    speed_dist = dist_next - dist_prev
                else:
                    speed_lon = 0.0
                    speed_lat = 0.0
                    speed_dist = 0.0

                result = (lon, lat, dist, speed_lon, speed_lat, speed_dist)
                result = _apply_fixstar_flags(
                    result, jd_tt, flags, j2000_native=use_j2000
                )
                results[index] = (result, canonical_name, ret_flags)
            _leb_ok = True
        except KeyError:
            pass  # Body not in LEB file
        except ValueError as _leb_err:
            if "outside range" not in str(_leb_err).lower():
                raise

    if _leb_ok:
        return tuple(results)

    # Skyfield fallback path
    from .state import get_planets

    earth = get_planets()["earth"]
    earth_at_t = earth.at(t)

    if want_speed:
        t_prev = ts.tt_jd(jd_tt - 0.5)
        t_next = ts.tt_jd(jd_tt + 0.5)
        earth_at_prev = earth.at(t_prev)
        earth_at_next = earth.at(t_next)
    else:
        earth_at_prev = None
        earth_at_next = None

    for index, star_id, canonical_name in resolved:
        try:
            lon, lat, dist = _calc_star_position_from_observer(
                star_id, earth_at_t, noaberr, nogdefl, use_j2000
            )
            if want_speed:
                lon_prev, lat_prev, dist_prev = _calc_star_position_from_observer(
                    star_id, earth_at_prev, noaberr, nogdefl, use_j2000
                )
                lon_next, lat_next, dist_next = _calc_star_position_from_observer(
                    star_id, earth_at_next, noaberr, nogdefl, use_j2000
                )
                speed_lon = lon_next - lon_prev
                if speed_lon > 180.0:
                    speed_lon -= 360.0
                elif speed_lon < -180.0:
                    speed_lon += 360.0
                speed_lat = lat_next - lat_prev
                speed_dist = dist_next - dist_prev
            else:
                speed_lon = 0.0
                speed_lat = 0.0
                speed_dist = 0.0

            result = (lon, lat, dist, speed_lon, speed_lat, speed_dist)
            result = _apply_fixstar_flags(result, jd_tt, flags, j2000_native=use_j2000)
            results[index] = (result, canonical_name, ret_flags)
        except Error:
            if skip_errors:
                continue
            raise
        except (OSError, ValueError, KeyError) as e:
            if skip_errors:
                continue
            raise Error(str(e)) from e

    return tuple(results)


batch_fixstars_ut = batch_fixstars_ut


def fixstar(
    star: str, tjdet: float, flags: int = FLG_SWIEPH
) -> Tuple[Tuple[float, float, float, float, float, float], str, int]:
    """
    Calculate position of a fixed star for Terrestrial Time (TT).

    Reference API compatible function. Similar to fixstar_ut() but takes
    Terrestrial Time (TT, also known as Ephemeris Time) instead of Universal Time.

    Args:
        star: Name of star (e.g. "Regulus", "Spica")
        tjdet: Julian Day in Terrestrial Time (TT/ET)
        flags: Calculation flags

    Returns:
        Tuple containing:
            - Position tuple: (lon, lat, dist, speed_lon, speed_lat, speed_dist)
            - star_name: Resolved star name (e.g. "Regulus")
            - retflags: Return flags (int)

    Raises:
        Error: If the star cannot be found.

    Note:
        TT (Terrestrial Time) differs from UT (Universal Time) by Delta T,
        which varies from ~32 seconds (year 2000) to minutes (historical times).
        For most astrological applications, use fixstar_ut() instead.

    Example:
        >>> pos, name, retflag = fixstar("Regulus", 2451545.0, 0)
        >>> lon, lat, dist = pos[0], pos[1], pos[2]
    """
    ret_flags = _fixstar_ret_flags(flags)
    flags = _preprocess_flags(flags)

    star_id, error, canonical_name = _resolve_star_id(star)
    if error:
        raise Error(error)

    try:
        noaberr = bool(flags & FLG_NOABERR) or bool(flags & FLG_TRUEPOS)
        nogdefl = bool(flags & FLG_NOGDEFL)
        use_j2000 = bool(flags & FLG_J2000)
        topo = _fixstar_topo() if flags & FLG_TOPOCTR else None

        if flags & FLG_SPEED:
            lon, lat, dist, speed_lon, speed_lat, speed_dist = calc_fixed_star_velocity(
                star_id, tjdet, noaberr, nogdefl, j2000_frame=use_j2000, topo=topo
            )
            result = (lon, lat, dist, speed_lon, speed_lat, speed_dist)
        else:
            lon, lat, dist = calc_fixed_star_position(
                star_id, tjdet, noaberr, nogdefl, j2000_frame=use_j2000, topo=topo
            )
            result = (lon, lat, dist, 0.0, 0.0, 0.0)

        result = _apply_fixstar_flags(result, tjdet, flags, j2000_native=use_j2000)

        return (result, canonical_name or "", ret_flags)
    except Error:
        raise
    except (OSError, ValueError, KeyError) as e:
        raise Error(str(e)) from e


def _format_star_name(entry: StarCatalogEntry) -> str:
    """
    Format the full star name for return from fixstar2 functions.

    Returns format: "Name,Nomenclature" (e.g. "Regulus,alLeo")
    """
    return f"{entry.name},{entry.nomenclature}"


def _resolve_star2(star_name: str) -> Tuple[StarCatalogEntry | None, str | None]:
    """
    Resolve a star identifier with flexible lookup for fixstar2 functions.

    Supports multiple lookup methods:
    1. Exact star name (case-insensitive): "Regulus", "SPICA"
    2. Hipparcos catalog number (as string): "49669", ",49669"
    3. Hipparcos catalog number with HIP prefix: "HIP 49669", "HIP65474"
    4. Bayer designation with Greek letter names: "Alpha Leonis", "Beta Persei"
    5. Flamsteed designation with number + constellation: "32 Leonis", "87 Virginis"
    6. Partial name search (case-insensitive): "Reg", "pic"
    7. Bayer/Flamsteed nomenclature: "alLeo", "alVir"
    8. Format with comma: "Regulus,alLeo" (takes first part)
    9. Phonetic fuzzy matching for alternate spellings: "Betelgeux", "Formalhaut"

    Args:
        star_name: Star identifier - can be name, catalog number, or search string

    Returns:
        Tuple of (StarCatalogEntry, error_message). If error, entry is None.

    Examples:
        >>> entry, err = _resolve_star2("Regulus")         # Exact name
        >>> entry, err = _resolve_star2("49669")           # HIP number
        >>> entry, err = _resolve_star2(",49669")          # HIP with leading comma
        >>> entry, err = _resolve_star2("HIP 49669")       # HIP with prefix
        >>> entry, err = _resolve_star2("HIP65474")        # HIP with prefix (no space)
        >>> entry, err = _resolve_star2("Alpha Leonis")    # Bayer designation
        >>> entry, err = _resolve_star2("Beta Persei")     # Bayer designation
        >>> entry, err = _resolve_star2("32 Leonis")       # Flamsteed designation
        >>> entry, err = _resolve_star2("Reg")             # Partial match
        >>> entry, err = _resolve_star2("alLeo")           # Nomenclature
        >>> entry, err = _resolve_star2("Betelgeux")       # Alternate spelling
    """
    search = star_name.strip()

    if not search:
        return None, "Empty star name"

    # Check if it's a catalog number (numeric string, possibly with leading comma)
    number_search = search.lstrip(",").strip()
    if number_search.isdigit():
        hip_number = int(number_search)
        for entry in STAR_CATALOG:
            if entry.hip_number == hip_number:
                return entry, None
        return None, f"could not find star name {hip_number}"

    # Check for "HIP NNNNN" format (case-insensitive)
    search_upper_temp = search.upper()
    if search_upper_temp.startswith("HIP ") or search_upper_temp.startswith("HIP"):
        # Extract the numeric part after "HIP" (with or without space)
        hip_part = search[3:].strip()
        if hip_part.isdigit():
            hip_number = int(hip_part)
            for entry in STAR_CATALOG:
                if entry.hip_number == hip_number:
                    return entry, None
            return None, f"could not find star name HIP {hip_number}"

    # Handle comma-separated format (e.g., "Regulus,alLeo")
    if "," in search:
        search = search.split(",")[0].strip()

    search_upper = search.upper()

    # 1. Try exact name match (case-insensitive)
    for entry in STAR_CATALOG:
        if entry.name.upper() == search_upper:
            return entry, None

    # 2. Try exact nomenclature match (case-insensitive)
    for entry in STAR_CATALOG:
        if entry.nomenclature.upper() == search_upper:
            return entry, None

    # 3. Try Bayer designation with Greek letter names (e.g., "Alpha Leonis")
    parsed_nomenclature = _parse_bayer_designation(search)
    if parsed_nomenclature:
        for entry in STAR_CATALOG:
            if entry.nomenclature.upper() == parsed_nomenclature.upper():
                return entry, None

    # 3b. Try Flamsteed designation (e.g., "32 Leonis", "87 Virginis")
    parsed_flamsteed = _parse_flamsteed_designation(search)
    if parsed_flamsteed:
        # Look up in STAR_ALIASES using the normalized format (e.g., "32 LEO")
        if parsed_flamsteed in STAR_ALIASES:
            star_id = STAR_ALIASES[parsed_flamsteed]
            for entry in STAR_CATALOG:
                if entry.id == star_id:
                    return entry, None

    # 4. Try partial name match (prefix search, case-insensitive)
    matches: List[StarCatalogEntry] = []
    for entry in STAR_CATALOG:
        if entry.name.upper().startswith(search_upper):
            matches.append(entry)

    if len(matches) == 1:
        return matches[0], None
    elif len(matches) > 1:
        names = ", ".join(m.name for m in matches)
        return None, f"Ambiguous star name '{star_name}' matches: {names}"

    # 5. Try partial nomenclature match
    for entry in STAR_CATALOG:
        if entry.nomenclature.upper().startswith(search_upper):
            matches.append(entry)

    if len(matches) == 1:
        return matches[0], None
    elif len(matches) > 1:
        names = ", ".join(m.name for m in matches)
        return None, f"Ambiguous star name '{star_name}' matches: {names}"

    # 6. Try substring match in name (anywhere in the name)
    for entry in STAR_CATALOG:
        if search_upper in entry.name.upper():
            matches.append(entry)

    if len(matches) == 1:
        return matches[0], None
    elif len(matches) > 1:
        names = ", ".join(m.name for m in matches)
        return None, f"Ambiguous star name '{star_name}' matches: {names}"

    # 7. Try phonetic fuzzy matching for common misspellings
    fuzzy_result = _fuzzy_match_star(star_name)
    if fuzzy_result is not None:
        for entry in STAR_CATALOG:
            if entry.id == fuzzy_result:
                return entry, None

    return None, f"could not find star name {star_name.lower()}"


def fixstar2_ut(
    star: str, tjdut: float, flags: int = FLG_SWIEPH
) -> Tuple[Tuple[float, float, float, float, float, float], str, int]:
    """
    Calculate position of a fixed star for Universal Time with flexible lookup.

    Enhanced version of fixstar_ut() that supports flexible star lookup:
    - Star name (full or partial): "Regulus", "Reg"
    - Hipparcos catalog number: "49669", ",49669"
    - Hipparcos with HIP prefix: "HIP 49669", "HIP65474"
    - Bayer/Flamsteed designation: "alLeo", "alVir"

    Returns the position and full star name, allowing identification
    of which star was matched when using partial searches.

    Args:
        star: Star identifier (name, catalog number, or partial search)
        tjdut: Julian Day in Universal Time (UT1)
        flags: Calculation flags

    Returns:
        Tuple containing:
            - Position tuple: (lon, lat, dist, speed_lon, speed_lat, speed_dist)
            - star_name_out: Full star name "Name,Nomenclature"
              (e.g. "Regulus,alLeo")
            - retflags: Return flags (int)

    Raises:
        Error: If the star cannot be found.

    Note:
        UT (Universal Time) is converted to TT (Terrestrial Time) internally
        using Delta T before calculating the star position.

    Example:
        >>> pos, name, retflag = fixstar2_ut("Reg", 2451545.0, 0)
        >>> print(name)  # "Regulus,alLeo"
        >>> lon, lat, dist = pos[0], pos[1], pos[2]

        >>> pos, name, retflag = fixstar2_ut("49669", 2451545.0, 0)
        >>> print(name)  # "Regulus,alLeo" (looked up by HIP number)
    """
    flags = _preprocess_flags(flags)

    entry, error = _resolve_star2(star)
    if error or entry is None:
        raise Error(error or "could not find star name")

    from .state import get_timescale

    ts = get_timescale()
    t = ts.ut1_jd(tjdut)

    try:
        noaberr = bool(flags & FLG_NOABERR) or bool(flags & FLG_TRUEPOS)
        nogdefl = bool(flags & FLG_NOGDEFL)
        use_j2000 = bool(flags & FLG_J2000)

        if flags & FLG_SPEED:
            lon, lat, dist, speed_lon, speed_lat, speed_dist = calc_fixed_star_velocity(
                entry.id, t.tt, noaberr, nogdefl, j2000_frame=use_j2000
            )
            result = (lon, lat, dist, speed_lon, speed_lat, speed_dist)
        else:
            lon, lat, dist = calc_fixed_star_position(
                entry.id, t.tt, noaberr, nogdefl, j2000_frame=use_j2000
            )
            result = (lon, lat, dist, 0.0, 0.0, 0.0)

        star_name_out = _format_star_name(entry)
        result = _apply_fixstar_flags(result, t.tt, flags, j2000_native=use_j2000)

        return (result, star_name_out, flags)
    except Error:
        raise
    except (OSError, ValueError, KeyError) as e:
        raise Error(str(e)) from e


def fixstar2(
    star: str, tjdet: float, flags: int = FLG_SWIEPH
) -> Tuple[Tuple[float, float, float, float, float, float], str, int]:
    """
    Calculate position of a fixed star for Terrestrial Time with flexible lookup.

    Enhanced version of fixstar() that supports flexible star lookup:
    - Star name (full or partial): "Regulus", "Reg"
    - Hipparcos catalog number: "49669", ",49669"
    - Bayer/Flamsteed designation: "alLeo", "alVir"

    Returns the position and full star name, allowing identification
    of which star was matched when using partial searches.

    Args:
        star: Star identifier (name, catalog number, or partial search)
        tjdet: Julian Day in Terrestrial Time (TT/ET)
        flags: Calculation flags

    Returns:
        Tuple containing:
            - Position tuple: (lon, lat, dist, speed_lon, speed_lat, speed_dist)
            - star_name_out: Full star name "Name,Nomenclature"
              (e.g. "Spica,alVir")
            - retflags: Return flags (int)

    Raises:
        Error: If the star cannot be found.

    Note:
        TT (Terrestrial Time) differs from UT (Universal Time) by Delta T.
        For most astrological applications, use fixstar2_ut() instead.

    Example:
        >>> pos, name, retflag = fixstar2("Spica", 2451545.0, 0)
        >>> print(name)  # "Spica,alVir"
        >>> lon, lat, dist = pos[0], pos[1], pos[2]

        >>> pos, name, retflag = fixstar2("65474", 2451545.0, 0)
        >>> print(name)  # "Spica,alVir" (looked up by HIP number)
    """
    flags = _preprocess_flags(flags)

    entry, error = _resolve_star2(star)
    if error or entry is None:
        raise Error(error or "could not find star name")

    try:
        noaberr = bool(flags & FLG_NOABERR) or bool(flags & FLG_TRUEPOS)
        nogdefl = bool(flags & FLG_NOGDEFL)
        use_j2000 = bool(flags & FLG_J2000)

        if flags & FLG_SPEED:
            lon, lat, dist, speed_lon, speed_lat, speed_dist = calc_fixed_star_velocity(
                entry.id, tjdet, noaberr, nogdefl, j2000_frame=use_j2000
            )
            result = (lon, lat, dist, speed_lon, speed_lat, speed_dist)
        else:
            lon, lat, dist = calc_fixed_star_position(
                entry.id, tjdet, noaberr, nogdefl, j2000_frame=use_j2000
            )
            result = (lon, lat, dist, 0.0, 0.0, 0.0)

        star_name_out = _format_star_name(entry)
        result = _apply_fixstar_flags(result, tjdet, flags, j2000_native=use_j2000)

        return (result, star_name_out, flags)
    except Error:
        raise
    except (OSError, ValueError, KeyError) as e:
        raise Error(str(e)) from e


# Magnitude values for legacy _resolve_star_id lookup
# Built from STAR_CATALOG so every cataloged star has a magnitude available.
_STAR_MAGNITUDES = {entry.id: entry.magnitude for entry in STAR_CATALOG}


# =============================================================================
# STAR NAME TO HIP NUMBER MAPPING
# =============================================================================
# Mapping from common star names, Bayer designations, and Flamsteed numbers
# to Hipparcos (HIP) catalog numbers.
#
# Data sourced from the IAU Working Group on Star Names (WGSN) catalog:
# https://www.iau.org/public/themes/naming_stars/
# IAU-CSN (IAU Catalog of Star Names) last updated 2022-04-04
#
# This mapping provides direct lookup from star names to HIP numbers,
# independent of the internal SE_* star ID constants.
# =============================================================================

STAR_NAME_TO_HIP: dict[str, int] = {
    # =========================================================================
    # COMMON/PROPER STAR NAMES (IAU-approved names)
    # =========================================================================
    # Names are stored in uppercase for case-insensitive lookup
    # Source: IAU-CSN catalog
    # =========================================================================
    # A
    "ABSOLUTNO": -1,  # XO-5, no HIP (exoplanet host)
    "ACAMAR": 13847,  # Theta1 Eridani
    "ACHERNAR": 7588,  # Alpha Eridani
    "ACHIRD": 3821,  # Eta Cassiopeiae
    "ACRAB": 78820,  # Beta Scorpii (also Graffias)
    "ACRUX": 60718,  # Alpha Crucis
    "ACUBENS": 44066,  # Alpha Cancri
    "ADHAFERA": 50335,  # Zeta Leonis
    "ADHARA": 33579,  # Epsilon Canis Majoris
    "ADHIL": 6411,  # Xi Andromedae
    "AIN": 20889,  # Epsilon Tauri
    "AINALRAMI": 92761,  # Nu1 Sagittarii
    "ALADFAR": 94481,  # Eta Lyrae
    "ALASIA": 90004,  # HD 168746
    "ALBALDAH": 94141,  # Pi Sagittarii
    "ALBALI": 102618,  # Epsilon Aquarii
    "ALBIREO": 95947,  # Beta Cygni
    "ALCHIBA": 59199,  # Alpha Corvi
    "ALCOR": 65477,  # 80 Ursae Majoris
    "ALCYONE": 17702,  # Eta Tauri (Pleiades)
    "ALDEBARAN": 21421,  # Alpha Tauri
    "ALDERAMIN": 105199,  # Alpha Cephei
    "ALDHANAB": 108085,  # Gamma Gruis
    "ALDHIBAH": 83895,  # Zeta Draconis
    "ALDULFIN": 101421,  # Epsilon Delphini
    "ALFIRK": 106032,  # Beta Cephei
    "ALGEDI": 100064,  # Alpha2 Capricorni
    "ALGENIB": 1067,  # Gamma Pegasi
    "ALGIEBA": 50583,  # Gamma1 Leonis
    "ALGOL": 14576,  # Beta Persei
    "ALGORAB": 60965,  # Delta Corvi
    "ALHENA": 31681,  # Gamma Geminorum
    "ALIOTH": 62956,  # Epsilon Ursae Majoris
    "ALJANAH": 102488,  # Epsilon Cygni
    "ALKAID": 67301,  # Eta Ursae Majoris
    "ALKALUROPS": 75411,  # Mu1 Bootis
    "ALKAPHRAH": 44471,  # Kappa Ursae Majoris
    "ALKARAB": 115623,  # Upsilon Pegasi
    "ALKES": 53740,  # Alpha Crateris
    "ALMAAZ": 23416,  # Epsilon Aurigae
    "ALMACH": 9640,  # Gamma Andromedae
    "ALNAIR": 109268,  # Alpha Gruis
    "ALNASL": 88635,  # Gamma2 Sagittarii
    "ALNILAM": 26311,  # Epsilon Orionis
    "ALNITAK": 26727,  # Zeta Orionis
    "ALNIYAT": 80112,  # Sigma Scorpii
    "ALPHARD": 46390,  # Alpha Hydrae
    "ALPHECCA": 76267,  # Alpha Coronae Borealis
    "ALPHERATZ": 677,  # Alpha Andromedae
    "ALPHERG": 7097,  # Eta Piscium
    "ALRAKIS": 83608,  # Mu Draconis
    "ALRESCHA": 9487,  # Alpha Piscium (HIP 7097 is Alpherg/Eta Piscium)
    "ALRUBA": 86782,  # Draconis
    "ALSAFI": 96100,  # Sigma Draconis
    "ALSCIAUKAT": 41075,  # 31 Lyncis
    "ALSEPHINA": 42913,  # Delta Velorum
    "ALSHAIN": 98036,  # Beta Aquilae
    "ALSHAT": 100310,  # Nu Capricorni
    "ALTAIR": 97649,  # Alpha Aquilae
    "ALTAIS": 94376,  # Delta Draconis
    "ALTERF": 46750,  # Lambda Leonis
    "ALUDRA": 35904,  # Eta Canis Majoris
    "ALULA AUSTRALIS": 55203,  # Xi Ursae Majoris
    "ALULA BOREALIS": 55219,  # Nu Ursae Majoris
    "ALYA": 92946,  # Theta1 Serpentis
    "ALZIRR": 32362,  # Xi Geminorum
    "AMADIOHA": 29550,  # HD 43197
    "AMANSINAYA": -1,  # WASP-34, no HIP
    "ANADOLU": -1,  # WASP-52, no HIP
    "ANCHA": 110003,  # Theta Aquarii
    "ANGETENAR": 13288,  # Tau2 Eridani
    "ANIARA": 57820,  # HD 102956
    "ANKAA": 2081,  # Alpha Phoenicis
    "ANSER": 95771,  # Alpha Vulpeculae
    "ANTARES": 80763,  # Alpha Scorpii
    "ARCALIS": 72845,  # HD 131496
    "ARCTURUS": 69673,  # Alpha Bootis
    "ARKAB POSTERIOR": 95294,  # Beta2 Sagittarii
    "ARKAB PRIOR": 95241,  # Beta1 Sagittarii
    "ARNEB": 25985,  # Alpha Leporis
    "ASCELLA": 93506,  # Zeta Sagittarii
    "ASELLUS AUSTRALIS": 42911,  # Delta Cancri (consistent with STAR_CATALOG)
    "ASELLUS BOREALIS": 42806,  # Gamma Cancri (corrected from 43103=Iota Cnc)
    "ASHLESHA": 43109,  # Epsilon Hydrae
    "ASPIDISKE": 45556,  # Iota Carinae
    "ASTEROPE": 17579,  # 21 Tauri (Pleiades)
    "ATAKORAKA": -1,  # WASP-64, no HIP
    "ATHEBYNE": 80331,  # Eta Draconis
    "ATIK": 17448,  # Omicron Persei
    "ATLAS": 17847,  # 27 Tauri (Pleiades)
    "ATRIA": 82273,  # Alpha Trianguli Australis
    "AVIOR": 41037,  # Epsilon Carinae
    "AXOLOTL": 118319,  # HD 224693
    "AYEYARWADY": 13993,  # HD 18742
    "AZELFAFAGE": 107136,  # Pi1 Cygni
    "AZHA": 13701,  # Eta Eridani
    "AZMIDI": 38170,  # Xi Puppis
    # B
    "BAEKDU": 73136,  # 8 Ursae Minoris
    "BARNARD'S STAR": 87937,  # GJ 699
    "BATEN KAITOS": 8645,  # Zeta Ceti
    "BEEMIM": 20535,  # Upsilon3 Eridani
    "BEID": 19587,  # Omicron1 Eridani
    "BELEL": 95124,  # HD 181342
    "BELENOS": 6643,  # HD 8574
    "BELLATRIX": 25336,  # Gamma Orionis
    "BEREHYNIA": -1,  # HAT-P-15, no HIP
    "BETELGEUSE": 27989,  # Alpha Orionis
    "BHARANI": 13209,  # 41 Arietis
    "BIBHA": 48711,  # HD 86081
    "BIHAM": 109427,  # Theta Pegasi
    "BOSONA": 107251,  # HD 206610
    "BOTEIN": 14838,  # Delta Arietis
    "BRACHIUM": 73714,  # Sigma Librae
    "BUBUP": 26380,  # HD 38283
    "BUNA": 12191,  # HD 16175
    "BUNDA": 106786,  # Xi Aquarii
    # C
    "CANOPUS": 30438,  # Alpha Carinae
    "CAPELLA": 24608,  # Alpha Aurigae
    "CAPH": 746,  # Beta Cassiopeiae
    "CASTOR": 36850,  # Alpha Geminorum
    "CASTULA": 4422,  # Upsilon2 Cassiopeiae
    "CEBALRAI": 86742,  # Beta Ophiuchi
    "CEIBO": 37284,  # HD 63454
    "CELAENO": 17489,  # 16 Tauri (Pleiades)
    "CERVANTES": 86796,  # Mu Arae
    "CHALAWAN": 53721,  # 47 Ursae Majoris
    "CHAMUKUY": 20894,  # Theta2 Tauri
    "CHAOPHRAYA": -1,  # WASP-50, no HIP
    "CHARA": 61317,  # Beta Canum Venaticorum
    "CHASON": -1,  # HAT-P-5, no HIP
    "CHECHIA": 99894,  # HD 192699
    "CHERTAN": 54879,  # Theta Leonis
    "CITADELLE": 1547,  # HD 1502
    "CITALA": 33719,  # HD 52265
    "COCIBOLCA": 3479,  # HD 4208
    "COPERNICUS": 43587,  # 55 Cancri
    "COR CAROLI": 63125,  # Alpha2 Canum Venaticorum
    "CUJAM": 80463,  # Omega Herculis
    "CURSA": 23875,  # Beta Eridani
    # D
    "DABIH": 100345,  # Beta1 Capricorni
    "DALIM": 14879,  # Alpha Fornacis
    "DENEB": 102098,  # Alpha Cygni
    "DENEB ALGEDI": 107556,  # Delta Capricorni
    "DENEBOLA": 57632,  # Beta Leonis
    "DIADEM": 64241,  # Alpha Comae Berenices
    "DINGOLAY": 54158,  # HD 96063
    "DIPHDA": 3419,  # Beta Ceti
    "DIWO": -1,  # WASP-17, no HIP
    "DIYA": -1,  # WASP-72, no HIP
    "DOFIDA": 66047,  # HD 117618
    "DOMBAY": -1,  # HAT-P-3, no HIP
    "DSCHUBBA": 78401,  # Delta Scorpii
    "DUBHE": 54061,  # Alpha Ursae Majoris
    "DZIBAN": 86614,  # Psi1 Draconis
    # E
    "EBLA": 114322,  # HD 218566
    "EDASICH": 75458,  # Iota Draconis
    "ELECTRA": 17499,  # 17 Tauri (Pleiades)
    "ELGAFAR": 70755,  # Phi Virginis
    "ELKURUD": 29034,  # Theta Columbae
    "ELNATH": 25428,  # Beta Tauri
    "ELTANIN": 87833,  # Gamma Draconis
    "EMIW": 5529,  # HD 7199
    "ENIF": 107315,  # Epsilon Pegasi
    "ERRAI": 116727,  # Gamma Cephei
    # F
    "FAFNIR": 90344,  # 42 Draconis
    "FANG": 78265,  # Pi Scorpii
    "FAWARIS": 97165,  # Delta Cygni
    "FELIS": 48615,  # HR 3923
    "FELIXVARELA": 2247,  # BD-17 63
    "FLEGETONTE": 57370,  # HD 102195
    "FOMALHAUT": 113368,  # Alpha Piscis Austrini
    "FORMOSA": 56508,  # HD 100655
    "FRANZ": -1,  # HAT-P-14, no HIP
    "FULU": 2920,  # Zeta Cassiopeiae
    "FUMALSAMAKAH": 113889,  # Beta Piscium
    "FUNI": 61177,  # HD 109246
    "FURUD": 30122,  # Zeta Canis Majoris
    "FUYUE": 87261,  # Scorpii
    # G
    "GACRUX": 61084,  # Gamma Crucis
    "GAKYID": 42446,  # HD 73534
    "GEMINGA": -1,  # Pulsar, no HIP
    "GIAUSAR": 56211,  # Lambda Draconis
    "GIENAH": 59803,  # Gamma Corvi
    "GINAN": 60260,  # Epsilon Crucis
    "GLOAS": -1,  # WASP-13, no HIP
    "GOMEISA": 36188,  # Beta Canis Minoris
    "GRUMIUM": 87585,  # Xi Draconis
    "GUDJA": 77450,  # Kappa Serpentis
    "GUMALA": 94645,  # HD 179949
    "GUNIIBUU": 84405,  # 36 Ophiuchi
    # H
    "HADAR": 68702,  # Beta Centauri
    "HAEDUS": 23767,  # Eta Aurigae
    "HAMAL": 9884,  # Alpha Arietis
    "HASSALEH": 23015,  # Iota Aurigae
    "HATYSA": 26241,  # Iota Orionis
    "HELVETIOS": 113357,  # 51 Pegasi
    "HEZE": 66249,  # Zeta Virginis
    "HOGGAR": 21109,  # HD 28678
    "HOMAM": 112029,  # Zeta Pegasi
    "HORNA": -1,  # HAT-P-38, no HIP
    "HUNAHPU": 55174,  # HD 98219
    "HUNOR": 80076,  # HAT-P-2
    # I
    "IKLIL": 78104,  # Rho Scorpii
    "ILLYRIAN": 47087,  # HD 82886
    "IMAI": 59747,  # Delta Crucis
    "INQUILL": 84787,  # HD 156411
    "INTAN": 15578,  # HD 20868
    "INTERCRUS": 46471,  # HR 3743
    "IRENA": -1,  # WASP-38, no HIP
    "ITONDA": 108375,  # HD 208487
    "IZAR": 72105,  # Epsilon Bootis
    # J
    "JABBAH": 79374,  # Nu Scorpii
    "JISHUI": 37265,  # Omicron Geminorum
    # K
    "KAFFALJIDHMA": 12706,  # Gamma Ceti
    "KALAUSI": 47202,  # HD 83443
    "KAMUY": 79219,  # HD 145457
    "KANG": 69427,  # Kappa Virginis
    "KARAKA": 76351,  # HD 137388
    "KAUS AUSTRALIS": 90185,  # Epsilon Sagittarii
    "KAUS BOREALIS": 90496,  # Lambda Sagittarii
    "KAUS MEDIA": 89931,  # Delta Sagittarii
    "KAVEH": 92895,  # HD 175541
    "KEID": 19849,  # Omicron2 Eridani
    "KHAMBALIA": 69974,  # Lambda Virginis
    "KITALPHA": 104987,  # Alpha Equulei
    "KOCHAB": 72607,  # Beta Ursae Minoris
    "KOEIA": 12961,  # HIP 12961
    "KOIT": -1,  # XO-4, no HIP
    "KORNEPHOROS": 80816,  # Beta Herculis
    "KRAZ": 61359,  # Beta Corvi
    "KURHAH": 108917,  # Xi Cephei
    # L
    "LA SUPERBA": 62223,  # Y Canum Venaticorum
    "LARAWAG": 82396,  # Epsilon Scorpii
    "LERNA": -1,  # HAT-P-42, no HIP
    "LESATH": 85696,  # Upsilon Scorpii
    "LIBERTAS": 97938,  # Xi Aquilae
    "LICH": -1,  # PSR B1257+12, pulsar, no HIP
    "LIESMA": 66192,  # HD 118203
    "LILII BOREA": 13061,  # 39 Arietis
    "LIONROCK": 110813,  # HD 212771
    "LUCILINBURHUC": 30860,  # HD 45350
    "LUSITANIA": 30905,  # HD 45652
    # M
    "MAASYM": 85693,  # Lambda Herculis
    "MACONDO": 52521,  # HD 93083
    "MAGO": 24003,  # HD 32518
    "MAHASIM": 28380,  # Theta Aurigae
    "MAHSATI": 82651,  # HD 152581
    "MAIA": 17573,  # 20 Tauri (Pleiades)
    "MALMOK": -1,  # WASP-39, no HIP
    "MARFIK": 80883,  # Lambda Ophiuchi
    "MARKAB": 113963,  # Alpha Pegasi
    "MARKEB": 45941,  # Kappa Velorum
    "MAROHU": -1,  # WASP-6, no HIP
    "MARSIC": 79043,  # Kappa Herculis
    "MATAR": 112158,  # Eta Pegasi
    "MAZAALAI": -1,  # HAT-P-21, no HIP
    "MEBSUTA": 32246,  # Epsilon Geminorum
    "MEGREZ": 59774,  # Delta Ursae Majoris
    "MEISSA": 26207,  # Lambda Orionis
    "MEKBUDA": 34088,  # Zeta Geminorum
    "MELEPH": 42556,  # Epsilon Cancri
    "MENKALINAN": 28360,  # Beta Aurigae
    "MENKAR": 14135,  # Alpha Ceti
    "MENKENT": 68933,  # Theta Centauri
    "MENKIB": 18614,  # Xi Persei
    "MERAK": 53910,  # Beta Ursae Majoris
    "MERGA": 72487,  # 38 Bootis
    "MERIDIANA": 94114,  # Alpha Coronae Australis
    "MEROPE": 17608,  # 23 Tauri (Pleiades)
    "MESARTHIM": 8832,  # Gamma2 Arietis
    "MIAPLACIDUS": 45238,  # Beta Carinae
    "MIMOSA": 62434,  # Beta Crucis
    "MINCHIR": 42402,  # Sigma Hydrae
    "MINELAUVA": 63090,  # Delta Virginis
    "MINTAKA": 25930,  # Delta Orionis
    "MIRA": 10826,  # Omicron Ceti
    "MIRACH": 5447,  # Beta Andromedae
    "MIRAM": 13268,  # Eta Persei
    "MIRFAK": 15863,  # Alpha Persei
    "MIRZAM": 30324,  # Beta Canis Majoris
    "MISAM": 14668,  # Kappa Persei
    "MIZAR": 65378,  # Zeta Ursae Majoris
    "MOLDOVEANU": -1,  # XO-1, no HIP
    "MONCH": 72339,  # HD 130322
    "MONTUNO": -1,  # WASP-79, no HIP
    "MORAVA": -1,  # WASP-60, no HIP
    "MORIAH": -1,  # HAT-P-23, no HIP
    "MOTHALLAH": 8796,  # Alpha Trianguli
    "MOUHOUN": 22491,  # HD 30856
    "MPINGO": -1,  # WASP-71, no HIP
    "MULIPHEIN": 34045,  # Gamma Canis Majoris
    "MUPHRID": 67927,  # Eta Bootis
    "MUSCIDA": 41704,  # Omicron Ursae Majoris
    "MUSICA": 103527,  # 18 Delphini
    "MUSPELHEIM": -1,  # HAT-P-29, no HIP
    # N
    "NAHN": 44946,  # Xi Cancri
    "NALEDI": -1,  # WASP-62, no HIP
    "NAOS": 39429,  # Zeta Puppis
    "NASHIRA": 106985,  # Gamma Capricorni
    "NASTI": 40687,  # HD 68988
    "NATASHA": 48235,  # HD 85390
    "NEKKAR": 73555,  # Beta Bootis
    "NEMBUS": 7607,  # 51 Andromedae
    "NENQUE": 5054,  # HD 6434
    "NERVIA": 32916,  # HD 49674
    "NGANURGANITY": 33856,  # Sigma Canis Majoris
    "NIHAL": 25606,  # Beta Leporis
    "NIKAWIY": 74961,  # HD 136418
    "NOSAXA": 31895,  # HD 48265
    "NUNKI": 92855,  # Sigma Sagittarii
    "NUSAKAN": 75695,  # Beta Coronae Borealis
    "NUSHAGAK": 13192,  # HD 17156
    "NYAMIEN": -1,  # WASP-15, no HIP
    # O
    "OGMA": 80838,  # HD 149026
    "OKAB": 93747,  # Zeta Aquilae
    # P
    "PAIKAUHALE": 81266,  # Tau Scorpii
    "PARUMLEO": -1,  # WASP-32, no HIP
    "PEACOCK": 100751,  # Alpha Pavonis
    "PETRA": -1,  # WASP-80, no HIP
    "PHACT": 26634,  # Alpha Columbae
    "PHECDA": 58001,  # Gamma Ursae Majoris
    "PHERKAD": 75097,  # Gamma Ursae Minoris
    "PHOENICIA": 99711,  # HD 192263
    "PIAUTOS": 40881,  # Lambda Cancri
    "PINCOYA": 88414,  # HD 164604
    "PIPIRIMA": 82545,  # Mu2 Scorpii
    "PIPOLTR": -1,  # TrES-3, no HIP
    "PLEIONE": 17851,  # 28 Tauri (Pleiades)
    "POERAVA": 116084,  # HD 221287
    "POLARIS": 11767,  # Alpha Ursae Minoris
    "POLARIS AUSTRALIS": 104382,  # Sigma Octantis
    "POLIS": 89341,  # Mu Sagittarii
    "POLLUX": 37826,  # Beta Geminorum
    "PORRIMA": 61941,  # Gamma Virginis
    "PRAECIPUA": 53229,  # 46 Leonis Minoris
    "PRIMA HYADUM": 20205,  # Gamma Tauri
    "PROCYON": 37279,  # Alpha Canis Minoris
    "PROPUS": 29655,  # Eta Geminorum
    "PROXIMA CENTAURI": 70890,  # GJ 551, Alpha Centauri C
    # R
    "RAN": 16537,  # Epsilon Eridani
    "RANA": 17378,  # Delta Eridani
    "RAPETO": 83547,  # HD 153950
    "RASALAS": 48455,  # Mu Leonis
    "RASALGETHI": 84345,  # Alpha1 Herculis
    "RASALHAGUE": 86032,  # Alpha Ophiuchi
    "RASTABAN": 85670,  # Beta Draconis
    "REGULUS": 49669,  # Alpha Leonis
    "REVATI": 5737,  # Zeta Piscium
    "RIGEL": 24436,  # Beta Orionis
    "RIGIL KENTAURUS": 71683,  # Alpha Centauri A
    "ROSALIADECASTRO": 81022,  # HD 149143
    "ROTANEV": 101769,  # Beta Delphini
    "RUCHBAH": 6686,  # Delta Cassiopeiae
    "RUKBAT": 95347,  # Alpha Sagittarii
    # S
    "SABIK": 84012,  # Eta Ophiuchi
    "SACLATENI": 23453,  # Zeta Aurigae
    "SADACHBIA": 110395,  # Gamma Aquarii
    "SADALBARI": 112748,  # Mu Pegasi
    "SADALMELIK": 109074,  # Alpha Aquarii
    "SADALSUUD": 106278,  # Beta Aquarii
    "SADR": 100453,  # Gamma Cygni
    "SAIPH": 27366,  # Kappa Orionis
    "SALM": 115250,  # Tau Pegasi
    "SARGAS": 86228,  # Theta Scorpii
    "SARIN": 84379,  # Delta Herculis
    "SCHEAT": 113881,  # Beta Pegasi
    "SCHEDAR": 3179,  # Alpha Cassiopeiae
    "SECUNDA HYADUM": 20455,  # Delta1 Tauri
    "SEGIN": 8886,  # Epsilon Cassiopeiae
    "SEGINUS": 71075,  # Gamma Bootis
    "SHAULA": 85927,  # Lambda Scorpii
    "SHAMA": 69701,  # HD 99109
    "SHERATAN": 8903,  # Beta Arietis
    "SIKA": 50782,  # HD 99491
    "SIRIUS": 32349,  # Alpha Canis Majoris
    "SITULA": 111710,  # Kappa Aquarii
    "SKAT": 113136,  # Delta Aquarii
    "SPICA": 65474,  # Alpha Virginis
    "STRIBOR": 91085,  # HD 171028
    "SUBRA": 47508,  # Omicron Leonis
    "SUHAIL": 44816,  # Lambda Velorum
    "SULAFAT": 93194,  # Gamma Lyrae
    "SYRMA": 69701,  # Iota Virginis
    # T
    "TABIT": 22449,  # Pi3 Orionis
    "TAIYANGSHOU": 57399,  # Chi Ursae Majoris
    "TAIYI": 63076,  # 8 Draconis
    "TALITHA": 44127,  # Iota Ursae Majoris
    "TANIA AUSTRALIS": 50801,  # Mu Ursae Majoris
    "TANIA BOREALIS": 50372,  # Lambda Ursae Majoris
    "TARAZED": 97278,  # Gamma Aquilae
    "TAYGETA": 17531,  # 19 Tauri (Pleiades)
    "TEBERDA": 94256,  # HD 178813
    "TEGMINE": 40167,  # Zeta1 Cancri
    "TEJAT": 30343,  # Mu Geminorum (HIP 30343, corrected from erroneous 32362)
    "THUBAN": 68756,  # Alpha Draconis
    "TIAKI": 112122,  # Beta Gruis
    "TIANGUAN": 26451,  # Zeta Tauri
    "TIANYI": 62423,  # 7 Draconis
    "TITAWIN": 7513,  # Upsilon Andromedae (HIP 7513; 9683 was a typo)
    "TOLIMAN": 71681,  # Alpha Centauri B
    "TONATIUH": 58952,  # HD 104985
    "TORCULAR": 8198,  # Omicron Piscium
    "TUREIS": 39757,  # Rho Puppis
    "TYL": 91919,  # Epsilon Draconis
    # U
    "UKDAH": 47431,  # Iota Hydrae
    "UNUKALHAI": 77070,  # Alpha Serpentis
    "UNURGUNITE": 34444,  # Sigma Canis Majoris
    "URUK": 116076,  # HD 231701
    # V
    "VEGA": 91262,  # Alpha Lyrae
    "VERITATE": 116076,  # 14 Andromedae
    "VINDEMIATRIX": 63608,  # Epsilon Virginis
    "WASAT": 35550,  # Delta Geminorum
    "WAZN": 27628,  # Beta Columbae
    "WEZEN": 34444,  # Delta Canis Majoris
    # X
    "XAMIDIMURA": 82514,  # Mu1 Scorpii
    # Y
    "YED POSTERIOR": 79882,  # Epsilon Ophiuchi
    "YED PRIOR": 79593,  # Delta Ophiuchi
    "YILDUN": 85822,  # Delta Ursae Minoris
    # Z
    "ZANIAH": 60129,  # Eta Virginis
    "ZAURAK": 18543,  # Gamma Eridani
    "ZAVIJAVA": 57757,  # Beta Virginis
    "ZHANG": 48356,  # Upsilon1 Hydrae
    "ZIBAL": 15197,  # Zeta Eridani
    "ZOSMA": 54872,  # Delta Leonis
    "ZUBENELGENUBI": 72622,  # Alpha2 Librae
    "ZUBENELHAKRABI": 76333,  # Gamma Librae
    "ZUBENESCHAMALI": 74785,  # Beta Librae
    # =========================================================================
    # BAYER DESIGNATIONS (Greek letter + constellation)
    # =========================================================================
    # Format: "ALPHA CONSTELLATION" and abbreviated forms
    # =========================================================================
    # Alpha stars
    "ALPHA ANDROMEDAE": 677,
    "ALPHA AND": 677,
    "ALPHA AQUARII": 109074,
    "ALPHA AQR": 109074,
    "ALPHA AQUILAE": 97649,
    "ALPHA AQL": 97649,
    "ALPHA ARIETIS": 9884,
    "ALPHA ARI": 9884,
    "ALPHA AURIGAE": 24608,
    "ALPHA AUR": 24608,
    "ALPHA BOOTIS": 69673,
    "ALPHA BOO": 69673,
    "ALPHA CANIS MAJORIS": 32349,
    "ALPHA CMA": 32349,
    "ALPHA CANIS MINORIS": 37279,
    "ALPHA CMI": 37279,
    "ALPHA CAPRICORNI": 100064,
    "ALPHA CAP": 100064,
    "ALPHA CARINAE": 30438,
    "ALPHA CAR": 30438,
    "ALPHA CASSIOPEIAE": 3179,
    "ALPHA CAS": 3179,
    "ALPHA CENTAURI": 71683,
    "ALPHA CEN": 71683,
    "ALPHA CEPHEI": 105199,
    "ALPHA CEP": 105199,
    "ALPHA CETI": 14135,
    "ALPHA CET": 14135,
    "ALPHA CORONAE BOREALIS": 76267,
    "ALPHA CRB": 76267,
    "ALPHA CRUCIS": 60718,
    "ALPHA CRU": 60718,
    "ALPHA CYGNI": 102098,
    "ALPHA CYG": 102098,
    "ALPHA DRACONIS": 68756,
    "ALPHA DRA": 68756,
    "ALPHA ERIDANI": 7588,
    "ALPHA ERI": 7588,
    "ALPHA GEMINORUM": 36850,
    "ALPHA GEM": 36850,
    "ALPHA GRUIS": 109268,
    "ALPHA GRU": 109268,
    "ALPHA HERCULIS": 84345,
    "ALPHA HER": 84345,
    "ALPHA HYDRAE": 46390,
    "ALPHA HYA": 46390,
    "ALPHA LEONIS": 49669,
    "ALPHA LEO": 49669,
    "ALPHA LIBRAE": 72622,
    "ALPHA LIB": 72622,
    "ALPHA LYRAE": 91262,
    "ALPHA LYR": 91262,
    "ALPHA OPHIUCHI": 86032,
    "ALPHA OPH": 86032,
    "ALPHA ORIONIS": 27989,
    "ALPHA ORI": 27989,
    "ALPHA PAVONIS": 100751,
    "ALPHA PAV": 100751,
    "ALPHA PEGASI": 113963,
    "ALPHA PEG": 113963,
    "ALPHA PERSEI": 15863,
    "ALPHA PER": 15863,
    "ALPHA PHOENICIS": 2081,
    "ALPHA PHE": 2081,
    "ALPHA PISCIS AUSTRINI": 113368,
    "ALPHA PSA": 113368,
    "ALPHA PISCIUM": 9487,
    "ALPHA PSC": 9487,
    "ALPHA SAGITTARII": 95347,
    "ALPHA SGR": 95347,
    "ALPHA SCORPII": 80763,
    "ALPHA SCO": 80763,
    "ALPHA SERPENTIS": 77070,
    "ALPHA SER": 77070,
    "ALPHA TAURI": 21421,
    "ALPHA TAU": 21421,
    "ALPHA TRIANGULI": 8796,
    "ALPHA TRI": 8796,
    "ALPHA URSAE MAJORIS": 54061,
    "ALPHA UMA": 54061,
    "ALPHA URSAE MINORIS": 11767,
    "ALPHA UMI": 11767,
    "ALPHA VIRGINIS": 65474,
    "ALPHA VIR": 65474,
    # Beta stars
    "BETA ANDROMEDAE": 5447,
    "BETA AND": 5447,
    "BETA AQUARII": 106278,
    "BETA AQR": 106278,
    "BETA AQUILAE": 98036,
    "BETA AQL": 98036,
    "BETA ARIETIS": 8903,
    "BETA ARI": 8903,
    "BETA AURIGAE": 28360,
    "BETA AUR": 28360,
    "BETA BOOTIS": 73555,
    "BETA BOO": 73555,
    "BETA CANIS MAJORIS": 30324,
    "BETA CMA": 30324,
    "BETA CANIS MINORIS": 36188,
    "BETA CMI": 36188,
    "BETA CAPRICORNI": 100345,
    "BETA CAP": 100345,
    "BETA CARINAE": 45238,
    "BETA CAR": 45238,
    "BETA CASSIOPEIAE": 746,
    "BETA CAS": 746,
    "BETA CENTAURI": 68702,
    "BETA CEN": 68702,
    "BETA CEPHEI": 106032,
    "BETA CEP": 106032,
    "BETA CETI": 3419,
    "BETA CET": 3419,
    "BETA CORONAE BOREALIS": 75695,
    "BETA CRB": 75695,
    "BETA CRUCIS": 62434,
    "BETA CRU": 62434,
    "BETA CYGNI": 95947,
    "BETA CYG": 95947,
    "BETA DELPHINI": 101769,
    "BETA DEL": 101769,
    "BETA DRACONIS": 85670,
    "BETA DRA": 85670,
    "BETA ERIDANI": 23875,
    "BETA ERI": 23875,
    "BETA GEMINORUM": 37826,
    "BETA GEM": 37826,
    "BETA GRUIS": 112122,
    "BETA GRU": 112122,
    "BETA HERCULIS": 80816,
    "BETA HER": 80816,
    "BETA LEONIS": 57632,
    "BETA LEO": 57632,
    "BETA LEPORIS": 25606,
    "BETA LEP": 25606,
    "BETA LIBRAE": 74785,
    "BETA LIB": 74785,
    "BETA LYRAE": 92420,
    "BETA LYR": 92420,
    "BETA OPHIUCHI": 86742,
    "BETA OPH": 86742,
    "BETA ORIONIS": 24436,
    "BETA ORI": 24436,
    "BETA PEGASI": 113881,
    "BETA PEG": 113881,
    "BETA PERSEI": 14576,
    "BETA PER": 14576,
    "BETA PISCIUM": 113889,
    "BETA PSC": 113889,
    "BETA SAGITTARII": 95241,
    "BETA SGR": 95241,
    "BETA SCORPII": 78820,
    "BETA SCO": 78820,
    "BETA TAURI": 25428,
    "BETA TAU": 25428,
    "BETA URSAE MAJORIS": 53910,
    "BETA UMA": 53910,
    "BETA URSAE MINORIS": 72607,
    "BETA UMI": 72607,
    "BETA VIRGINIS": 57757,
    "BETA VIR": 57757,
    # Gamma stars
    "GAMMA ANDROMEDAE": 9640,
    "GAMMA AND": 9640,
    "GAMMA AQUARII": 110395,
    "GAMMA AQR": 110395,
    "GAMMA AQUILAE": 95501,
    "GAMMA AQL": 95501,
    "GAMMA ARIETIS": 8832,
    "GAMMA ARI": 8832,
    "GAMMA BOOTIS": 71075,
    "GAMMA BOO": 71075,
    "GAMMA CANCRI": 42806,
    "GAMMA CNC": 42806,
    "GAMMA CAPRICORNI": 106985,
    "GAMMA CAP": 106985,
    "GAMMA CASSIOPEIAE": 4427,
    "GAMMA CAS": 4427,
    "GAMMA CEPHEI": 116727,
    "GAMMA CEP": 116727,
    "GAMMA CORVI": 59803,
    "GAMMA CRV": 59803,
    "GAMMA CRUCIS": 61084,
    "GAMMA CRU": 61084,
    "GAMMA CYGNI": 100453,
    "GAMMA CYG": 100453,
    "GAMMA DRACONIS": 87833,
    "GAMMA DRA": 87833,
    "GAMMA GEMINORUM": 31681,
    "GAMMA GEM": 31681,
    "GAMMA GRUIS": 108085,
    "GAMMA GRU": 108085,
    "GAMMA LEONIS": 50583,
    "GAMMA LEO": 50583,
    "GAMMA LYRAE": 93194,
    "GAMMA LYR": 93194,
    "GAMMA ORIONIS": 25336,
    "GAMMA ORI": 25336,
    "GAMMA PEGASI": 1067,
    "GAMMA PEG": 1067,
    "GAMMA SAGITTARII": 88635,
    "GAMMA SGR": 88635,
    "GAMMA URSAE MAJORIS": 58001,
    "GAMMA UMA": 58001,
    "GAMMA URSAE MINORIS": 75097,
    "GAMMA UMI": 75097,
    "GAMMA VIRGINIS": 61941,
    "GAMMA VIR": 61941,
    # Delta stars
    "DELTA AQUARII": 113136,
    "DELTA AQR": 113136,
    "DELTA BOOTIS": 72659,
    "DELTA BOO": 72659,
    "DELTA CANCRI": 42911,
    "DELTA CNC": 42911,
    "DELTA CAPRICORNI": 107556,
    "DELTA CAP": 107556,
    "DELTA CASSIOPEIAE": 6686,
    "DELTA CAS": 6686,
    "DELTA CEPHEI": 110991,
    "DELTA CEP": 110991,
    "DELTA CORVI": 60965,
    "DELTA CRV": 60965,
    "DELTA CRUCIS": 59747,
    "DELTA CRU": 59747,
    "DELTA CYGNI": 97165,
    "DELTA CYG": 97165,
    "DELTA DRACONIS": 94376,
    "DELTA DRA": 94376,
    "DELTA ERIDANI": 17378,
    "DELTA ERI": 17378,
    "DELTA GEMINORUM": 35550,
    "DELTA GEM": 35550,
    "DELTA HERCULIS": 79992,
    "DELTA HER": 79992,
    "DELTA LEONIS": 54872,
    "DELTA LEO": 54872,
    "DELTA ORIONIS": 25930,
    "DELTA ORI": 25930,
    "DELTA SAGITTARII": 89931,
    "DELTA SGR": 89931,
    "DELTA SCORPII": 78401,
    "DELTA SCO": 78401,
    "DELTA TAURI": 20455,
    "DELTA TAU": 20455,
    "DELTA URSAE MAJORIS": 59774,
    "DELTA UMA": 59774,
    "DELTA URSAE MINORIS": 85822,
    "DELTA UMI": 85822,
    "DELTA VIRGINIS": 63090,
    "DELTA VIR": 63090,
    # Epsilon stars
    "EPSILON AQUARII": 102618,
    "EPSILON AQR": 102618,
    "EPSILON AURIGAE": 23416,
    "EPSILON AUR": 23416,
    "EPSILON BOOTIS": 72105,
    "EPSILON BOO": 72105,
    "EPSILON CANIS MAJORIS": 33579,
    "EPSILON CMA": 33579,
    "EPSILON CANCRI": 42556,
    "EPSILON CNC": 42556,
    "EPSILON CARINAE": 41037,
    "EPSILON CAR": 41037,
    "EPSILON CENTAURI": 66657,
    "EPSILON CEN": 66657,
    "EPSILON CYGNI": 102488,
    "EPSILON CYG": 102488,
    "EPSILON DRACONIS": 91919,
    "EPSILON DRA": 91919,
    "EPSILON ERIDANI": 16537,
    "EPSILON ERI": 16537,
    "EPSILON GEMINORUM": 32246,
    "EPSILON GEM": 32246,
    "EPSILON HYDRAE": 43109,
    "EPSILON HYA": 43109,
    "EPSILON LEONIS": 47908,
    "EPSILON LEO": 47908,
    "EPSILON OPHIUCHI": 86284,
    "EPSILON OPH": 86284,
    "EPSILON ORIONIS": 26311,
    "EPSILON ORI": 26311,
    "EPSILON PEGASI": 107315,
    "EPSILON PEG": 107315,
    "EPSILON SAGITTARII": 90185,
    "EPSILON SGR": 90185,
    "EPSILON SCORPII": 82396,
    "EPSILON SCO": 82396,
    "EPSILON TAURI": 20889,
    "EPSILON TAU": 20889,
    "EPSILON URSAE MAJORIS": 62956,
    "EPSILON UMA": 62956,
    "EPSILON VIRGINIS": 63608,
    "EPSILON VIR": 63608,
    # Zeta stars
    "ZETA AQUARII": 110960,
    "ZETA AQR": 110960,
    "ZETA AQUILAE": 93747,
    "ZETA AQL": 93747,
    "ZETA AURIGAE": 23453,
    "ZETA AUR": 23453,
    "ZETA CANIS MAJORIS": 30122,
    "ZETA CMA": 30122,
    "ZETA CASSIOPEIAE": 2920,
    "ZETA CAS": 2920,
    "ZETA CENTAURI": 68002,
    "ZETA CEN": 68002,
    "ZETA CEPHEI": 109492,
    "ZETA CEP": 109492,
    "ZETA DRACONIS": 83895,
    "ZETA DRA": 83895,
    "ZETA ERIDANI": 20535,
    "ZETA ERI": 20535,
    "ZETA GEMINORUM": 34088,
    "ZETA GEM": 34088,
    "ZETA HERCULIS": 81693,
    "ZETA HER": 81693,
    "ZETA LEONIS": 50335,
    "ZETA LEO": 50335,
    "ZETA OPHIUCHI": 81377,
    "ZETA OPH": 81377,
    "ZETA ORIONIS": 26727,
    "ZETA ORI": 26727,
    "ZETA PEGASI": 112029,
    "ZETA PEG": 112029,
    "ZETA PISCIUM": 5737,
    "ZETA PSC": 5737,
    "ZETA PUPPIS": 39429,
    "ZETA PUP": 39429,
    "ZETA SAGITTARII": 93506,
    "ZETA SGR": 93506,
    "ZETA SCORPII": 82671,
    "ZETA SCO": 82671,
    "ZETA TAURI": 26451,
    "ZETA TAU": 26451,
    "ZETA URSAE MAJORIS": 65378,
    "ZETA UMA": 65378,
    "ZETA VIRGINIS": 66249,
    "ZETA VIR": 66249,
    # Eta stars
    "ETA AQUARII": 102618,
    "ETA AQR": 102618,
    "ETA AURIGAE": 23767,
    "ETA AUR": 23767,
    "ETA BOOTIS": 67927,
    "ETA BOO": 67927,
    "ETA CANIS MAJORIS": 35904,
    "ETA CMA": 35904,
    "ETA CASSIOPEIAE": 3821,
    "ETA CAS": 3821,
    "ETA CENTAURI": 71352,
    "ETA CEN": 71352,
    "ETA CEPHEI": 102422,
    "ETA CEP": 102422,
    "ETA DRACONIS": 80331,
    "ETA DRA": 80331,
    "ETA GEMINORUM": 29655,
    "ETA GEM": 29655,
    "ETA HERCULIS": 81833,
    "ETA HER": 81833,
    "ETA LEONIS": 47908,
    "ETA LEO": 47908,
    "ETA LYRAE": 94481,
    "ETA LYR": 94481,
    "ETA OPHIUCHI": 84012,
    "ETA OPH": 84012,
    "ETA ORIONIS": 25281,
    "ETA ORI": 25281,
    "ETA PEGASI": 112158,
    "ETA PEG": 112158,
    "ETA PERSEI": 13268,
    "ETA PER": 13268,
    "ETA PISCIUM": 7097,  # (consistent with STAR_CATALOG)
    "ETA PSC": 7097,
    "ETA SAGITTARII": 89642,
    "ETA SGR": 89642,
    "ETA SCORPII": 84143,
    "ETA SCO": 84143,
    "ETA TAURI": 17702,
    "ETA TAU": 17702,
    "ETA URSAE MAJORIS": 67301,
    "ETA UMA": 67301,
    "ETA VIRGINIS": 60129,
    "ETA VIR": 60129,
    # Theta stars
    "THETA AQUARII": 110003,
    "THETA AQR": 110003,
    "THETA AURIGAE": 28380,
    "THETA AUR": 28380,
    "THETA BOOTIS": 70497,
    "THETA BOO": 70497,
    "THETA CENTAURI": 68933,
    "THETA CEN": 68933,
    "THETA DRACONIS": 78527,
    "THETA DRA": 78527,
    "THETA ERIDANI": 13847,
    "THETA ERI": 13847,
    "THETA LEONIS": 54879,
    "THETA LEO": 54879,
    "THETA OPHIUCHI": 83000,
    "THETA OPH": 83000,
    "THETA PEGASI": 109427,
    "THETA PEG": 109427,
    "THETA SCORPII": 86228,
    "THETA SCO": 86228,
    "THETA TAURI": 20894,
    "THETA TAU": 20894,
    "THETA URSAE MAJORIS": 46853,
    "THETA UMA": 46853,
    "THETA VIRGINIS": 64238,
    "THETA VIR": 64238,
    # Iota stars
    "IOTA AURIGAE": 23015,
    "IOTA AUR": 23015,
    "IOTA CARINAE": 45556,
    "IOTA CAR": 45556,
    "IOTA DRACONIS": 75458,
    "IOTA DRA": 75458,
    "IOTA ORIONIS": 26241,
    "IOTA ORI": 26241,
    "IOTA URSAE MAJORIS": 44127,
    "IOTA UMA": 44127,
    "IOTA VIRGINIS": 71957,
    "IOTA VIR": 71957,
    # Kappa stars
    "KAPPA AQUILAE": 96483,
    "KAPPA AQL": 96483,
    "KAPPA BOOTIS": 69481,
    "KAPPA BOO": 69481,
    "KAPPA HERCULIS": 79043,
    "KAPPA HER": 79043,
    "KAPPA OPHIUCHI": 86032,
    "KAPPA OPH": 86032,
    "KAPPA ORIONIS": 27366,
    "KAPPA ORI": 27366,
    "KAPPA SCORPII": 86670,
    "KAPPA SCO": 86670,
    "KAPPA URSAE MAJORIS": 44471,
    "KAPPA UMA": 44471,
    "KAPPA VELORUM": 45941,
    "KAPPA VEL": 45941,
    # Lambda stars
    "LAMBDA DRACONIS": 56211,
    "LAMBDA DRA": 56211,
    "LAMBDA GEMINORUM": 32362,
    "LAMBDA GEM": 32362,
    "LAMBDA HERCULIS": 85693,
    "LAMBDA HER": 85693,
    "LAMBDA LEONIS": 46750,
    "LAMBDA LEO": 46750,
    "LAMBDA OPHIUCHI": 80883,
    "LAMBDA OPH": 80883,
    "LAMBDA ORIONIS": 26207,
    "LAMBDA ORI": 26207,
    "LAMBDA SAGITTARII": 90496,
    "LAMBDA SGR": 90496,
    "LAMBDA SCORPII": 85927,
    "LAMBDA SCO": 85927,
    "LAMBDA URSAE MAJORIS": 50801,
    "LAMBDA UMA": 50801,
    "LAMBDA VELORUM": 44816,
    "LAMBDA VEL": 44816,
    "LAMBDA VIRGINIS": 69974,
    "LAMBDA VIR": 69974,
    # Mu stars
    "MU ARAE": 86796,
    "MU ARA": 86796,
    "MU BOOTIS": 75411,
    "MU BOO": 75411,
    "MU CEPHEI": 107259,
    "MU CEP": 107259,
    "MU DRACONIS": 83608,
    "MU DRA": 83608,
    "MU GEMINORUM": 30343,
    "MU GEM": 30343,
    "MU HERCULIS": 86974,
    "MU HER": 86974,
    "MU LEONIS": 48455,
    "MU LEO": 48455,
    "MU PEGASI": 112748,
    "MU PEG": 112748,
    "MU SAGITTARII": 89341,
    "MU SGR": 89341,
    "MU SCORPII": 82514,
    "MU SCO": 82514,
    "MU URSAE MAJORIS": 51250,
    "MU UMA": 51250,
    "MU VIRGINIS": 71957,
    "MU VIR": 71957,
    # Nu stars
    "NU ANDROMEDAE": 4436,
    "NU AND": 4436,
    "NU OPHIUCHI": 88048,
    "NU OPH": 88048,
    "NU SCORPII": 79374,
    "NU SCO": 79374,
    "NU URSAE MAJORIS": 55219,
    "NU UMA": 55219,
    # Xi stars
    "XI AQUILAE": 97938,
    "XI AQL": 97938,
    "XI BOOTIS": 72659,
    "XI BOO": 72659,
    "XI DRACONIS": 87585,
    "XI DRA": 87585,
    "XI GEMINORUM": 32362,
    "XI GEM": 32362,
    "XI PERSEI": 18614,
    "XI PER": 18614,
    "XI PUPPIS": 38170,
    "XI PUP": 38170,
    "XI URSAE MAJORIS": 55203,
    "XI UMA": 55203,
    # Omicron stars
    "OMICRON ANDROMEDAE": 3092,
    "OMICRON AND": 3092,
    "OMICRON CETI": 10826,
    "OMICRON CET": 10826,
    "OMICRON LEONIS": 47508,
    "OMICRON LEO": 47508,
    "OMICRON PERSEI": 17448,
    "OMICRON PER": 17448,
    "OMICRON URSAE MAJORIS": 41704,
    "OMICRON UMA": 41704,
    # Pi stars
    "PI SAGITTARII": 94141,
    "PI SGR": 94141,
    "PI SCORPII": 78265,
    "PI SCO": 78265,
    # Rho stars
    "RHO PUPPIS": 42913,
    "RHO PUP": 42913,
    "RHO SCORPII": 78104,
    "RHO SCO": 78104,
    # Sigma stars
    "SIGMA DRACONIS": 96100,
    "SIGMA DRA": 96100,
    "SIGMA LIBRAE": 73714,
    "SIGMA LIB": 73714,
    "SIGMA OCTANTIS": 104382,
    "SIGMA OCT": 104382,
    "SIGMA SAGITTARII": 92855,
    "SIGMA SGR": 92855,
    "SIGMA SCORPII": 80112,
    "SIGMA SCO": 80112,
    # Tau stars
    "TAU PEGASI": 98066,
    "TAU PEG": 98066,
    "TAU SCORPII": 81266,
    "TAU SCO": 81266,
    # Upsilon stars
    "UPSILON ANDROMEDAE": 9683,
    "UPSILON AND": 9683,
    "UPSILON PEGASI": 115623,
    "UPSILON PEG": 115623,
    "UPSILON SCORPII": 85696,
    "UPSILON SCO": 85696,
    # Phi/Chi/Psi/Omega stars
    "PHI VIRGINIS": 70755,
    "PHI VIR": 70755,
    "CHI URSAE MAJORIS": 54539,
    "CHI UMA": 54539,
    "PSI DRACONIS": 86614,
    "PSI DRA": 86614,
    "OMEGA HERCULIS": 80463,
    "OMEGA HER": 80463,
    # =========================================================================
    # FLAMSTEED DESIGNATIONS (Number + constellation)
    # =========================================================================
    # Selected bright stars with Flamsteed numbers
    # =========================================================================
    "16 TAURI": 17489,  # Celaeno
    "17 TAURI": 17499,  # Electra
    "19 TAURI": 17531,  # Taygeta
    "20 TAURI": 17573,  # Maia
    "21 TAURI": 17579,  # Asterope
    "23 TAURI": 17608,  # Merope
    "25 TAURI": 17702,  # Alcyone (Eta Tau)
    "27 TAURI": 17847,  # Atlas
    "28 TAURI": 17851,  # Pleione
    "38 BOOTIS": 72487,  # Merga
    "46 LEONIS MINORIS": 53229,  # Praecipua
    "47 URSAE MAJORIS": 53721,  # Chalawan
    "51 ANDROMEDAE": 7607,  # Nembus
    "51 PEGASI": 113357,  # Helvetios
    "55 CANCRI": 43587,  # Copernicus
    "80 URSAE MAJORIS": 65477,  # Alcor
    "85 URSAE MAJORIS": 67301,  # Alkaid (same as Eta UMa)
    # =========================================================================
    # ALTERNATIVE SPELLINGS AND HISTORICAL NAMES
    # =========================================================================
    "AGENA": 68702,  # Alternative for Hadar (Beta Centauri)
    "BECRUX": 62434,  # Alternative for Mimosa (Beta Crucis)
    "BENETNASH": 67301,  # Alternative for Alkaid
    "BUNGULA": 71683,  # Historical name for Rigil Kentaurus
    "COR LEONIS": 49669,  # Heart of the Lion = Regulus
    "DOG STAR": 32349,  # Common name for Sirius
    "NORTH STAR": 11767,  # Common name for Polaris
    "POLE STAR": 11767,  # Common name for Polaris
    "WEGA": 91262,  # Historical spelling of Vega
}


def get_hip_from_star_name(name: str) -> int | None:
    """
    Look up the Hipparcos (HIP) catalog number for a star name.

    Supports common/proper star names, Bayer designations, Flamsteed numbers,
    and alternative spellings. The lookup is case-insensitive.

    Args:
        name: Star name, Bayer designation, or Flamsteed number
              Examples: "Regulus", "Alpha Leonis", "Alpha Leo", "51 Pegasi"

    Returns:
        HIP catalog number if found, None if star not in mapping.
        Returns -1 for stars without HIP numbers (e.g., exoplanet hosts
        discovered by transit surveys).

    Examples:
        >>> get_hip_from_star_name("Regulus")
        49669
        >>> get_hip_from_star_name("Alpha Leo")
        49669
        >>> get_hip_from_star_name("alpha leonis")
        49669
        >>> get_hip_from_star_name("51 Pegasi")
        113357
        >>> get_hip_from_star_name("Unknown Star")
        None

    Note:
        Data sourced from IAU Working Group on Star Names (WGSN) catalog.
        See: https://www.iau.org/public/themes/naming_stars/
    """
    if not name:
        return None

    # Normalize: uppercase and strip whitespace
    normalized = name.upper().strip()

    # Direct lookup
    if normalized in STAR_NAME_TO_HIP:
        return STAR_NAME_TO_HIP[normalized]

    return None


def fixstar_mag(star: str) -> Tuple[float, str]:
    """
    Get the visual magnitude of a fixed star without calculating position.

    Lightweight function that returns only the magnitude, useful for
    visibility calculations where position is not needed.

    Compatible with pyswisseph: returns (magnitude, star_name) on success,
    raises Error if the star is not found.

    Args:
        star: Name of star (e.g. "Regulus", "Spica")

    Returns:
        Tuple containing:
            - magnitude: Visual magnitude (float, e.g. 1.40 for Regulus)
            - star_name_out: Full star name "Name,Nomenclature"
              (e.g. "Regulus,alLeo")

    Raises:
        Error: If the star cannot be found or magnitude is unavailable.

    Example:
        >>> mag, name = fixstar_mag("Regulus")
        >>> print(f"{name}: {mag}")  # "Regulus,alLeo: 1.40"
    """
    star_id, error, canonical_name = _resolve_star_id(star)
    if error:
        raise Error(error)

    if star_id not in _STAR_MAGNITUDES:
        raise Error(f"Magnitude not available for star ID {star_id}")

    # Build "Name,Nomenclature" format matching pyswisseph
    for entry in STAR_CATALOG:
        if entry.id == star_id:
            star_name_out = _format_star_name(entry)
            return (_STAR_MAGNITUDES[star_id], star_name_out)

    # Fallback: use canonical name if catalog entry not found
    name_out = canonical_name or star
    return (_STAR_MAGNITUDES[star_id], name_out)


def fixstar2_mag(star: str) -> Tuple[float, str]:
    """
    Get the visual magnitude of a fixed star with flexible lookup.

    Enhanced version that supports flexible star lookup like fixstar2:
    - Star name (full or partial): "Regulus", "Reg"
    - Hipparcos catalog number: "49669", ",49669"
    - Bayer/Flamsteed designation: "alLeo", "alVir"

    Returns the magnitude and the full star name, useful for
    visibility calculations where position is not needed.

    Compatible with pyswisseph: returns (magnitude, star_name) on success,
    raises Error if the star is not found.

    Args:
        star: Star identifier (name, catalog number, or partial search)

    Returns:
        Tuple containing:
            - magnitude: Visual magnitude (float)
            - star_name_out: Full star name "Name,Nomenclature"
              (e.g. "Regulus,alLeo")

    Raises:
        Error: If the star cannot be found.

    Example:
        >>> mag, name = fixstar2_mag("Reg")
        >>> print(f"{name}: {mag}")  # "Regulus,alLeo: 1.40"

        >>> mag, name = fixstar2_mag("49669")
        >>> print(f"{name}: {mag}")  # "Regulus,alLeo: 1.40"
    """
    entry, error = _resolve_star2(star)
    if error or entry is None:
        raise Error(error or "could not find star name")

    star_name_out = _format_star_name(entry)
    return (entry.magnitude, star_name_out)
