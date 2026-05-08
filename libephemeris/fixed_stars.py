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
    SE_REGULUS,
    SE_SPICA_STAR,
    SE_ALGOL,
    SE_SIRIUS,
    SE_ALDEBARAN,
    SE_ANTARES,
    SE_VEGA,
    SE_POLARIS,
    SE_FOMALHAUT,
    SE_BETELGEUSE,
    SE_RIGEL,
    SE_PROCYON,
    SE_CAPELLA,
    SE_ARCTURUS,
    SE_DENEB,
    SE_POLLUX,
    SE_CASTOR,
    SE_ALTAIR,
    SE_ACHERNAR,
    SE_CANOPUS,
    SE_ACRUX,
    SE_MIMOSA,
    SE_GACRUX,
    SE_HADAR,
    SE_RIGIL_KENT,
    SE_SHAULA,
    SE_BELLATRIX,
    SE_ELNATH,
    SE_MIRA,
    SE_ALNILAM,
    SE_ALNITAK,
    SE_MINTAKA,
    SE_SAIPH,
    SE_DIPHDA,
    SE_ALPHARD,
    SE_RASALHAGUE,
    SE_ETAMIN,
    SE_KOCHAB,
    SE_ALKAID,
    SE_DUBHE,
    SE_MERAK,
    SE_ALIOTH,
    SE_MIZAR,
    SE_ALCOR,
    SE_VINDEMIATRIX,
    SE_ZUBENELGENUBI,
    SE_ZUBENESCHAMALI,
    SE_UNUKALHAI,
    SE_ALGIEBA,
    SE_DENEBOLA,
    SE_MARKAB,
    SE_SCHEAT,
    SE_ALCYONE,
    SE_ALGORAB,
    SE_ALPHECCA,
    SE_DENEB_ALGEDI,
    SE_ASTEROPE,
    SE_CELAENO,
    SE_ELECTRA,
    SE_MAIA,
    SE_MEROPE,
    SE_TAYGETA,
    SE_ATLAS,
    SE_PLEIONE,
    SE_PRIMA_HYADUM,
    SE_SECUNDA_HYADUM,
    SE_THETA_TAURI,
    SE_AIN,
    SE_MEISSA,
    SE_PHECDA,
    SE_MEGREZ,
    SE_DELTA_CRUCIS,
    SE_MENKENT,
    SE_MUHLIFAIN,
    SE_EPSILON_CENTAURI,
    SE_ETA_CENTAURI,
    SE_ZETA_CENTAURI,
    SE_SARGAS,
    SE_DSCHUBBA,
    SE_GRAFFIAS,
    SE_LESATH,
    SE_ZOSMA,
    SE_HAMAL,
    SE_SHERATAN,
    SE_MESARTHIM,
    SE_ACUBENS,
    SE_TARF,
    SE_ASELLUS_BOREALIS,
    SE_ASELLUS_AUSTRALIS,
    SE_KAUS_AUSTRALIS,
    SE_NUNKI,
    SE_KAUS_MEDIA,
    SE_KAUS_BOREALIS,
    SE_ASCELLA,
    SE_ALGEDI,
    SE_DABIH,
    SE_NASHIRA,
    SE_SADALSUUD,
    SE_SADALMELIK,
    SE_SKAT,
    SE_ETA_PISCIUM,
    SE_ALRESCHA,
    SE_ALPHERATZ,
    SE_ALGENIB,
    SE_PROPUS,
    SE_ALHENA,
    SE_TEJAT,
    SE_WASAT,
    SE_ADHARA,
    SE_WEZEN,
    SE_THUBAN,
    SE_RASALGETHI,
    SE_ALBIREO,
    SE_MIRACH,
    SE_ALMACH,
    SE_MENKAR,
    SEFLG_SWIEPH,
    SEFLG_SPEED,
    SEFLG_SPEED3,
    SEFLG_NOABERR,
    SEFLG_NOGDEFL,
    SEFLG_EQUATORIAL,
    SEFLG_J2000,
    SEFLG_NONUT,
    SEFLG_SIDEREAL,
    SEFLG_MOSEPH,
    SEFLG_XYZ,
    SEFLG_RADIANS,
    SEFLG_ICRS,
    SEFLG_TRUEPOS,
    SEFLG_TOPOCTR,
    J2000,
    J1991_25,
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


# Extended star catalog with names and catalog numbers
STAR_CATALOG: List[StarCatalogEntry] = [
    StarCatalogEntry(
        id=SE_REGULUS,
        name="Regulus",
        nomenclature="alLeo",
        hip_number=49669,
        data=StarData(
            ra_j2000=152.092958,  # 10h 08m 22.3s (Alpha Leonis)
            dec_j2000=11.967208,  # +11° 58' 02"
            pm_ra=-0.24873,  # -248.73 mas/yr (westward)
            pm_dec=0.00559,  # +5.59 mas/yr (northward)
            parallax_mas=41.13,
            radial_km_per_s=-18.06,
        ),
        magnitude=1.40,  # Visual magnitude
    ),
    StarCatalogEntry(
        id=SE_SPICA_STAR,
        name="Spica",
        nomenclature="alVir",
        hip_number=65474,
        data=StarData(
            ra_j2000=201.298247,  # 13h 25m 11.6s (Alpha Virginis)
            dec_j2000=-11.161319,  # -11° 09' 41"
            pm_ra=-0.04235,  # -42.35 mas/yr
            pm_dec=-0.03067,  # -30.67 mas/yr
            parallax_mas=13.06,
            radial_km_per_s=-22.85,
        ),
        magnitude=1.04,  # Visual magnitude
    ),
    StarCatalogEntry(
        id=SE_ALGOL,
        name="Algol",
        nomenclature="bePer",
        hip_number=14576,
        data=StarData(
            ra_j2000=47.042219,  # 03h 08m 10.1s
            dec_j2000=40.955647,  # +40° 57' 20"
            pm_ra=0.00299,  # 2.99 mas/yr
            pm_dec=-0.00166,  # -1.66 mas/yr
            parallax_mas=36.27,
            radial_km_per_s=21.84,
        ),
        magnitude=2.12,
    ),
    StarCatalogEntry(
        id=SE_SIRIUS,
        name="Sirius",
        nomenclature="alCMa",
        hip_number=32349,
        data=StarData(
            ra_j2000=101.287155,  # 06h 45m 08.9s
            dec_j2000=-16.716116,  # -16° 42' 58"
            pm_ra=-0.54601,  # -546.01 mas/yr
            pm_dec=-1.22307,  # -1223.07 mas/yr
            parallax_mas=379.21,
            radial_km_per_s=-7.14,
        ),
        magnitude=-1.46,
    ),
    StarCatalogEntry(
        id=SE_ALDEBARAN,
        name="Aldebaran",
        nomenclature="alTau",
        hip_number=21421,
        data=StarData(
            ra_j2000=68.980163,  # 04h 35m 55.2s
            dec_j2000=16.509302,  # +16° 30' 33"
            pm_ra=0.06345,  # 63.45 mas/yr
            pm_dec=-0.18894,  # -188.94 mas/yr
            parallax_mas=48.94,
            radial_km_per_s=68.1,
        ),
        magnitude=0.85,
    ),
    StarCatalogEntry(
        id=SE_ANTARES,
        name="Antares",
        nomenclature="alSco",
        hip_number=80763,
        data=StarData(
            ra_j2000=247.351915,  # 16h 29m 24.5s
            dec_j2000=-26.432003,  # -26° 25' 55"
            pm_ra=-0.01211,  # -12.11 mas/yr
            pm_dec=-0.02330,  # -23.30 mas/yr
            parallax_mas=5.89,
            radial_km_per_s=-13.31,
        ),
        magnitude=1.06,
    ),
    StarCatalogEntry(
        id=SE_VEGA,
        name="Vega",
        nomenclature="alLyr",
        hip_number=91262,
        data=StarData(
            ra_j2000=279.234735,  # 18h 36m 56.3s
            dec_j2000=38.783689,  # +38° 47' 01"
            pm_ra=0.20094,  # 200.94 mas/yr
            pm_dec=0.28623,  # 286.23 mas/yr
            parallax_mas=130.23,
            radial_km_per_s=-19.21,
        ),
        magnitude=0.03,
    ),
    StarCatalogEntry(
        id=SE_POLARIS,
        name="Polaris",
        nomenclature="alUMi",
        hip_number=11767,
        data=StarData(
            ra_j2000=37.954561,  # 02h 31m 49.1s
            dec_j2000=89.264109,  # +89° 15' 51"
            pm_ra=0.04448,  # 44.48 mas/yr
            pm_dec=-0.01185,  # -11.85 mas/yr
            parallax_mas=7.54,
            radial_km_per_s=-14.06,
        ),
        magnitude=1.98,
    ),
    StarCatalogEntry(
        id=SE_FOMALHAUT,
        name="Fomalhaut",
        nomenclature="alPsA",
        hip_number=113368,
        data=StarData(
            ra_j2000=344.412693,  # 22h 57m 39.0s
            dec_j2000=-29.622237,  # -29° 37' 20"
            pm_ra=0.32895,  # 328.95 mas/yr
            pm_dec=-0.16467,  # -164.67 mas/yr
            parallax_mas=129.81,
            radial_km_per_s=28.9,
        ),
        magnitude=1.16,
    ),
    StarCatalogEntry(
        id=SE_BETELGEUSE,
        name="Betelgeuse",
        nomenclature="alOri",
        hip_number=27989,
        data=StarData(
            ra_j2000=88.792939,  # 05h 55m 10.3s
            dec_j2000=7.407064,  # +07° 24' 25"
            pm_ra=0.02754,  # 27.54 mas/yr
            pm_dec=0.01130,  # 11.30 mas/yr
            parallax_mas=6.55,
            radial_km_per_s=26.74,
        ),
        magnitude=0.42,
    ),
    StarCatalogEntry(
        id=SE_RIGEL,
        name="Rigel",
        nomenclature="beOri",
        hip_number=24436,
        data=StarData(
            ra_j2000=78.634467,  # 05h 14m 32.3s
            dec_j2000=-8.201638,  # -08° 12' 06"
            pm_ra=0.00131,  # 1.31 mas/yr
            pm_dec=0.00050,  # 0.50 mas/yr
            parallax_mas=3.78,
            radial_km_per_s=25.0,
        ),
        magnitude=0.12,
    ),
    StarCatalogEntry(
        id=SE_PROCYON,
        name="Procyon",
        nomenclature="alCMi",
        hip_number=37279,
        data=StarData(
            ra_j2000=114.825493,  # 07h 39m 18.1s
            dec_j2000=5.224993,  # +05° 13' 30"
            pm_ra=-0.71459,  # -714.59 mas/yr
            pm_dec=-1.03680,  # -1036.80 mas/yr
            parallax_mas=284.56,
            radial_km_per_s=-11.18,
        ),
        magnitude=0.34,
    ),
    StarCatalogEntry(
        id=SE_CAPELLA,
        name="Capella",
        nomenclature="alAur",
        hip_number=24608,
        data=StarData(
            ra_j2000=79.172328,  # 05h 16m 41.4s
            dec_j2000=45.997991,  # +45° 59' 53"
            pm_ra=0.07525,  # 75.25 mas/yr
            pm_dec=-0.42689,  # -426.89 mas/yr
            parallax_mas=76.2,
            radial_km_per_s=38.45,
        ),
        magnitude=0.08,
    ),
    StarCatalogEntry(
        id=SE_ARCTURUS,
        name="Arcturus",
        nomenclature="alBoo",
        hip_number=69673,
        data=StarData(
            ra_j2000=213.915300,  # 14h 15m 39.7s
            dec_j2000=19.182409,  # +19° 10' 57"
            pm_ra=-1.09339,  # -1093.39 mas/yr
            pm_dec=-2.00006,  # -2000.06 mas/yr
            parallax_mas=88.83,
            radial_km_per_s=-31.79,
        ),
        magnitude=-0.04,
    ),
    StarCatalogEntry(
        id=SE_DENEB,
        name="Deneb",
        nomenclature="alCyg",
        hip_number=102098,
        data=StarData(
            ra_j2000=310.357980,  # 20h 41m 25.9s
            dec_j2000=45.280339,  # +45° 16' 49"
            pm_ra=0.00201,  # 2.01 mas/yr
            pm_dec=0.00185,  # 1.85 mas/yr
            parallax_mas=2.31,
            radial_km_per_s=2.55,
        ),
        magnitude=1.25,
    ),
    StarCatalogEntry(
        id=SE_POLLUX,
        name="Pollux",
        nomenclature="beGem",
        hip_number=37826,
        data=StarData(
            ra_j2000=116.328958,  # 07h 45m 18.9s
            dec_j2000=28.026199,  # +28° 01' 34"
            pm_ra=-0.62655,  # -626.55 mas/yr
            pm_dec=-0.04580,  # -45.80 mas/yr
            parallax_mas=96.54,
            radial_km_per_s=-4.11,
        ),
        magnitude=1.14,
    ),
    StarCatalogEntry(
        id=SE_CASTOR,
        name="Castor",
        nomenclature="alGem",
        hip_number=36850,
        data=StarData(
            ra_j2000=113.649428,  # 07h 34m 35.9s
            dec_j2000=31.888276,  # +31° 53' 18"
            pm_ra=-0.19145,  # -191.45 mas/yr
            pm_dec=-0.14519,  # -145.19 mas/yr
            parallax_mas=64.12,
            radial_km_per_s=0.6,
        ),
        magnitude=1.58,
    ),
    StarCatalogEntry(
        id=SE_ALTAIR,
        name="Altair",
        nomenclature="alAql",
        hip_number=97649,
        data=StarData(
            ra_j2000=297.695827,  # 19h 50m 47.0s
            dec_j2000=8.868321,  # +08° 52' 06"
            pm_ra=0.53623,  # 536.23 mas/yr
            pm_dec=0.38529,  # 385.29 mas/yr
            parallax_mas=194.95,
            radial_km_per_s=-17.0,
        ),
        magnitude=0.77,
    ),
    StarCatalogEntry(
        id=SE_ACHERNAR,
        name="Achernar",
        nomenclature="alEri",
        hip_number=7588,
        data=StarData(
            ra_j2000=24.428523,  # 01h 37m 42.8s
            dec_j2000=-57.236753,  # -57° 14' 12"
            pm_ra=0.08700,  # 87.00 mas/yr
            pm_dec=-0.03824,  # -38.24 mas/yr
            parallax_mas=23.39,
            radial_km_per_s=28.46,
        ),
        magnitude=0.46,
    ),
    StarCatalogEntry(
        id=SE_CANOPUS,
        name="Canopus",
        nomenclature="alCar",
        hip_number=30438,
        data=StarData(
            ra_j2000=95.987958,  # 06h 23m 57.1s
            dec_j2000=-52.695661,  # -52° 41' 44"
            pm_ra=0.01993,  # 19.93 mas/yr
            pm_dec=0.02324,  # 23.24 mas/yr
            parallax_mas=10.55,
            radial_km_per_s=21.67,
        ),
        magnitude=-0.72,
    ),
    StarCatalogEntry(
        id=SE_ACRUX,
        name="Acrux",
        nomenclature="alCru",
        hip_number=60718,
        data=StarData(
            ra_j2000=186.649563,  # 12h 26m 35.9s
            dec_j2000=-63.099093,  # -63° 05' 57"
            pm_ra=-0.03583,  # -35.83 mas/yr
            pm_dec=-0.01486,  # -14.86 mas/yr
            parallax_mas=10.13,
            radial_km_per_s=-0.64,
        ),
        magnitude=0.76,
    ),
    StarCatalogEntry(
        id=SE_MIMOSA,
        name="Mimosa",
        nomenclature="beCru",
        hip_number=62434,
        data=StarData(
            ra_j2000=191.930263,  # 12h 47m 43.3s
            dec_j2000=-59.688764,  # -59° 41' 19"
            pm_ra=-0.04297,  # -42.97 mas/yr
            pm_dec=-0.01618,  # -16.18 mas/yr
            parallax_mas=11.71,
            radial_km_per_s=-7.57,
        ),
        magnitude=1.25,
    ),
    StarCatalogEntry(
        id=SE_GACRUX,
        name="Gacrux",
        nomenclature="gaCru",
        hip_number=61084,
        data=StarData(
            ra_j2000=187.791498,  # 12h 31m 09.9s
            dec_j2000=-57.113213,  # -57° 06' 48"
            pm_ra=0.02823,  # 28.23 mas/yr
            pm_dec=-0.26508,  # -265.08 mas/yr
            parallax_mas=36.83,
            radial_km_per_s=4.97,
        ),
        magnitude=1.64,
    ),
    StarCatalogEntry(
        id=SE_HADAR,
        name="Hadar",
        nomenclature="beCen",
        hip_number=68702,
        data=StarData(
            ra_j2000=210.955856,  # 14h 03m 49.4s (ICRS, Hipparcos HIP 68702)
            dec_j2000=-60.373035,  # -60° 22' 22.9"
            pm_ra=-0.03327,  # -33.27 mas/yr (Hipparcos, van Leeuwen 2007)
            pm_dec=-0.02316,  # -23.16 mas/yr (Hipparcos, van Leeuwen 2007)
            parallax_mas=8.32,
            radial_km_per_s=-3.74,
        ),
        magnitude=0.61,
    ),
    StarCatalogEntry(
        id=SE_RIGIL_KENT,
        name="Rigil Kentaurus",
        nomenclature="alCen",
        hip_number=71683,
        data=StarData(
            ra_j2000=219.902066,  # 14h 39m 36.5s
            dec_j2000=-60.833976,  # -60° 50' 02"
            pm_ra=-3.67925,  # -3679.25 mas/yr (Hipparcos, van Leeuwen 2007)
            pm_dec=0.47367,  # 473.67 mas/yr (Hipparcos, van Leeuwen 2007)
            parallax_mas=742.0,
            radial_km_per_s=-37.06,
        ),
        magnitude=-0.27,
    ),
    StarCatalogEntry(
        id=SE_SHAULA,
        name="Shaula",
        nomenclature="laSco",
        hip_number=85927,
        data=StarData(
            ra_j2000=263.402167,  # 17h 33m 36.5s
            dec_j2000=-37.103824,  # -37° 06' 14"
            pm_ra=-0.00853,  # -8.53 mas/yr
            pm_dec=-0.03080,  # -30.80 mas/yr
            parallax_mas=5.71,
            radial_km_per_s=-3.34,
        ),
        magnitude=1.62,
    ),
    StarCatalogEntry(
        id=SE_BELLATRIX,
        name="Bellatrix",
        nomenclature="gaOri",
        hip_number=25336,
        data=StarData(
            ra_j2000=81.282764,  # 05h 25m 07.9s
            dec_j2000=6.349703,  # +06° 20' 59"
            pm_ra=-0.00811,  # -8.11 mas/yr
            pm_dec=-0.01288,  # -12.88 mas/yr
            parallax_mas=12.92,
            radial_km_per_s=26.8,
        ),
        magnitude=1.64,
    ),
    StarCatalogEntry(
        id=SE_ELNATH,
        name="Elnath",
        nomenclature="beTau",
        hip_number=25428,
        data=StarData(
            ra_j2000=81.572971,  # 05h 26m 17.5s
            dec_j2000=28.607452,  # +28° 36' 27"
            pm_ra=0.02276,  # 22.76 mas/yr
            pm_dec=-0.17358,  # -173.58 mas/yr
            parallax_mas=24.36,
            radial_km_per_s=21.7,
        ),
        magnitude=1.65,
    ),
    StarCatalogEntry(
        id=SE_MIRA,
        name="Mira",
        nomenclature="omCet",
        hip_number=10826,
        data=StarData(
            ra_j2000=34.836617,  # 02h 19m 20.8s
            dec_j2000=-2.977640,  # -02° 58' 40"
            pm_ra=0.00933,  # 9.33 mas/yr
            pm_dec=-0.23736,  # -237.36 mas/yr
            parallax_mas=10.91,
            radial_km_per_s=88.96,
        ),
        magnitude=3.04,
    ),
    StarCatalogEntry(
        id=SE_ALNILAM,
        name="Alnilam",
        nomenclature="epOri",
        hip_number=26311,
        data=StarData(
            ra_j2000=84.053389,  # 05h 36m 12.8s
            dec_j2000=-1.201919,  # -01° 12' 07"
            pm_ra=0.00144,  # 1.44 mas/yr
            pm_dec=-0.00078,  # -0.78 mas/yr
            parallax_mas=1.65,
            radial_km_per_s=32.26,
        ),
        magnitude=1.69,
    ),
    StarCatalogEntry(
        id=SE_ALNITAK,
        name="Alnitak",
        nomenclature="zeOri",
        hip_number=26727,
        data=StarData(
            ra_j2000=85.189694,  # 05h 40m 45.5s
            dec_j2000=-1.942574,  # -01° 56' 33"
            pm_ra=0.00319,  # 3.19 mas/yr
            pm_dec=0.00203,  # 2.03 mas/yr
            parallax_mas=4.43,
            radial_km_per_s=-0.8,
        ),
        magnitude=1.74,
    ),
    StarCatalogEntry(
        id=SE_MINTAKA,
        name="Mintaka",
        nomenclature="deOri",
        hip_number=25930,
        data=StarData(
            ra_j2000=83.001667,  # 05h 32m 00.4s
            dec_j2000=-0.299095,  # -00° 17' 57"
            pm_ra=0.00064,  # 0.64 mas/yr
            pm_dec=-0.00069,  # -0.69 mas/yr
            parallax_mas=4.71,
            radial_km_per_s=24.07,
        ),
        magnitude=2.23,
    ),
    StarCatalogEntry(
        id=SE_SAIPH,
        name="Saiph",
        nomenclature="kaOri",
        hip_number=27366,
        data=StarData(
            ra_j2000=86.939119,  # 05h 47m 45.4s
            dec_j2000=-9.669605,  # -09° 40' 11"
            pm_ra=0.00146,  # 1.46 mas/yr
            pm_dec=-0.00128,  # -1.28 mas/yr
            parallax_mas=5.04,
            radial_km_per_s=27.05,
        ),
        magnitude=2.06,
    ),
    StarCatalogEntry(
        id=SE_MEISSA,
        name="Meissa",
        nomenclature="laOri",
        hip_number=26207,
        data=StarData(
            ra_j2000=83.784486,  # 05h 35m 08.3s
            dec_j2000=9.934156,  # +09° 56' 03"
            pm_ra=-0.00034,  # -0.34 mas/yr
            pm_dec=-0.00294,  # -2.94 mas/yr
            parallax_mas=2.97,
            radial_km_per_s=38.31,
        ),
        magnitude=3.33,
    ),
    StarCatalogEntry(
        id=SE_DIPHDA,
        name="Diphda",
        nomenclature="beCet",
        hip_number=3419,
        data=StarData(
            ra_j2000=10.897379,  # 00h 43m 35.4s
            dec_j2000=-17.986606,  # -17° 59' 12"
            pm_ra=0.23255,  # 232.55 mas/yr
            pm_dec=0.03199,  # 31.99 mas/yr
            parallax_mas=33.86,
            radial_km_per_s=41.98,
        ),
        magnitude=2.04,
    ),
    StarCatalogEntry(
        id=SE_ALPHARD,
        name="Alphard",
        nomenclature="alHya",
        hip_number=46390,
        data=StarData(
            ra_j2000=141.896847,  # 09h 27m 35.2s
            dec_j2000=-8.658602,  # -08° 39' 31"
            pm_ra=-0.01523,  # -15.23 mas/yr
            pm_dec=0.03437,  # 34.37 mas/yr
            parallax_mas=18.09,
            radial_km_per_s=-24.73,
        ),
        magnitude=1.98,
    ),
    StarCatalogEntry(
        id=SE_RASALHAGUE,
        name="Rasalhague",
        nomenclature="alOph",
        hip_number=86032,
        data=StarData(
            ra_j2000=263.733627,  # 17h 34m 56.1s
            dec_j2000=12.560035,  # +12° 33' 36"
            pm_ra=0.10807,  # 108.07 mas/yr
            pm_dec=-0.22157,  # -221.57 mas/yr
            parallax_mas=67.13,
            radial_km_per_s=4.76,
        ),
        magnitude=2.07,
    ),
    StarCatalogEntry(
        id=SE_ETAMIN,
        name="Etamin",
        nomenclature="gaDra",
        hip_number=87833,
        data=StarData(
            ra_j2000=269.151541,  # 17h 56m 36.4s
            dec_j2000=51.488896,  # +51° 29' 20"
            pm_ra=-0.00848,  # -8.48 mas/yr
            pm_dec=-0.02279,  # -22.79 mas/yr
            parallax_mas=21.14,
            radial_km_per_s=-24.85,
        ),
        magnitude=2.23,
    ),
    StarCatalogEntry(
        id=SE_KOCHAB,
        name="Kochab",
        nomenclature="beUMi",
        hip_number=72607,
        data=StarData(
            ra_j2000=222.676357,  # 14h 50m 42.3s
            dec_j2000=74.155504,  # +74° 09' 20"
            pm_ra=-0.03261,  # -32.61 mas/yr
            pm_dec=0.01142,  # 11.42 mas/yr
            parallax_mas=24.91,
            radial_km_per_s=16.17,
        ),
        magnitude=2.08,
    ),
    StarCatalogEntry(
        id=SE_ALKAID,
        name="Alkaid",
        nomenclature="etUMa",
        hip_number=67301,
        data=StarData(
            ra_j2000=206.885157,  # 13h 47m 32.4s
            dec_j2000=49.313267,  # +49° 18' 48"
            pm_ra=-0.12117,  # -121.17 mas/yr
            pm_dec=-0.01491,  # -14.91 mas/yr
            parallax_mas=31.38,
            radial_km_per_s=-29.2,
        ),
        magnitude=1.86,
    ),
    StarCatalogEntry(
        id=SE_DUBHE,
        name="Dubhe",
        nomenclature="alUMa",
        hip_number=54061,
        data=StarData(
            ra_j2000=165.931965,  # 11h 03m 43.7s
            dec_j2000=61.751035,  # +61° 45' 03"
            pm_ra=-0.13411,  # -134.11 mas/yr
            pm_dec=-0.03470,  # -34.70 mas/yr
            parallax_mas=26.54,
            radial_km_per_s=-22.81,
        ),
        magnitude=1.79,
    ),
    StarCatalogEntry(
        id=SE_MERAK,
        name="Merak",
        nomenclature="beUMa",
        hip_number=53910,
        data=StarData(
            ra_j2000=165.460319,  # 11h 01m 50.5s
            dec_j2000=56.382426,  # +56° 22' 57"
            pm_ra=0.08143,  # 81.43 mas/yr
            pm_dec=0.03349,  # 33.49 mas/yr
            parallax_mas=40.9,
            radial_km_per_s=-25.75,
        ),
        magnitude=2.37,
    ),
    StarCatalogEntry(
        id=SE_ALIOTH,
        name="Alioth",
        nomenclature="epUMa",
        hip_number=62956,
        data=StarData(
            ra_j2000=193.507290,  # 12h 54m 01.7s
            dec_j2000=55.959823,  # +55° 57' 35"
            pm_ra=0.11191,  # 111.91 mas/yr
            pm_dec=-0.00824,  # -8.24 mas/yr
            parallax_mas=39.51,
            radial_km_per_s=-25.76,
        ),
        magnitude=1.77,
    ),
    StarCatalogEntry(
        id=SE_MIZAR,
        name="Mizar",
        nomenclature="zeUMa",
        hip_number=65378,
        data=StarData(
            ra_j2000=200.981429,  # 13h 23m 55.5s
            dec_j2000=54.925362,  # +54° 55' 31"
            pm_ra=0.11901,  # 119.01 mas/yr
            pm_dec=-0.02597,  # -25.97 mas/yr
            parallax_mas=38.01,
            radial_km_per_s=-23.08,
        ),
        magnitude=2.23,
    ),
    StarCatalogEntry(
        id=SE_ALCOR,
        name="Alcor",
        nomenclature="80UMa",
        hip_number=65477,
        data=StarData(
            ra_j2000=201.306403,  # 13h 25m 13.5s
            dec_j2000=54.987954,  # +54° 59' 17"
            pm_ra=0.12021,  # 120.21 mas/yr
            pm_dec=-0.01604,  # -16.04 mas/yr
            parallax_mas=39.91,
            radial_km_per_s=-26.57,
        ),
        magnitude=3.99,
    ),
    StarCatalogEntry(
        id=SE_PHECDA,
        name="Phecda",
        nomenclature="gaUMa",
        hip_number=58001,
        data=StarData(
            ra_j2000=178.457679,  # 11h 53m 49.8s (Gamma Ursae Majoris)
            dec_j2000=53.694758,  # +53° 41' 41"
            pm_ra=0.10768,  # 107.68 mas/yr
            pm_dec=0.01101,  # 11.01 mas/yr
            parallax_mas=39.21,
            radial_km_per_s=-28.76,
        ),
        magnitude=2.44,
    ),
    StarCatalogEntry(
        id=SE_MEGREZ,
        name="Megrez",
        nomenclature="deUMa",
        hip_number=59774,
        data=StarData(
            ra_j2000=183.856503,  # 12h 15m 25.6s (Delta Ursae Majoris)
            dec_j2000=57.032615,  # +57° 01' 57"
            pm_ra=0.10411,  # 104.11 mas/yr
            pm_dec=0.00730,  # 7.30 mas/yr
            parallax_mas=40.51,
            radial_km_per_s=-32.48,
        ),
        magnitude=3.31,
    ),
    StarCatalogEntry(
        id=SE_VINDEMIATRIX,
        name="Vindemiatrix",
        nomenclature="epVir",
        hip_number=63608,
        data=StarData(
            ra_j2000=195.544157,  # 13h 02m 10.6s
            dec_j2000=10.959149,  # +10° 57' 33"
            pm_ra=-0.27380,  # -273.80 mas/yr
            pm_dec=0.01996,  # 19.96 mas/yr
            parallax_mas=29.76,
            radial_km_per_s=-43.32,
        ),
        magnitude=2.83,
    ),
    StarCatalogEntry(
        id=SE_ZUBENELGENUBI,
        name="Zubenelgenubi",
        nomenclature="alLib",
        hip_number=72622,
        data=StarData(
            ra_j2000=222.719638,  # 14h 50m 52.7s
            dec_j2000=-16.041778,  # -16° 02' 30"
            pm_ra=-0.10568,  # -105.68 mas/yr
            pm_dec=-0.06840,  # -68.40 mas/yr
            parallax_mas=43.03,
            radial_km_per_s=-36.82,
        ),
        magnitude=2.75,
    ),
    StarCatalogEntry(
        id=SE_ZUBENESCHAMALI,
        name="Zubeneschamali",
        nomenclature="beLib",
        hip_number=74785,
        data=StarData(
            ra_j2000=229.251724,  # 15h 17m 00.4s
            dec_j2000=-9.382914,  # -09° 22' 58"
            pm_ra=-0.09810,  # -98.10 mas/yr
            pm_dec=-0.01965,  # -19.65 mas/yr
            parallax_mas=17.62,
            radial_km_per_s=-60.47,
        ),
        magnitude=2.61,
    ),
    StarCatalogEntry(
        id=SE_UNUKALHAI,
        name="Unukalhai",
        nomenclature="alSer",
        hip_number=77070,
        data=StarData(
            ra_j2000=236.067089,  # 15h 44m 16.1s
            dec_j2000=6.425628,  # +06° 25' 32"
            pm_ra=0.13384,  # 133.84 mas/yr
            pm_dec=0.04481,  # 44.81 mas/yr
            parallax_mas=44.1,
            radial_km_per_s=-18.12,
        ),
        magnitude=2.65,
    ),
    StarCatalogEntry(
        id=SE_ALGIEBA,
        name="Algieba",
        nomenclature="gaLeo",
        hip_number=50583,
        data=StarData(
            ra_j2000=154.993144,  # 10h 19m 58.4s
            dec_j2000=19.841489,  # +19° 50' 30"
            pm_ra=0.30430,  # 304.30 mas/yr (Hipparcos, van Leeuwen 2007)
            pm_dec=-0.15428,  # -154.28 mas/yr (Hipparcos, van Leeuwen 2007)
            parallax_mas=25.96,
            radial_km_per_s=-57.44,
        ),
        magnitude=2.08,
    ),
    StarCatalogEntry(
        id=SE_DENEBOLA,
        name="Denebola",
        nomenclature="beLeo",
        hip_number=57632,
        data=StarData(
            ra_j2000=177.264910,  # 11h 49m 03.6s
            dec_j2000=14.572058,  # +14° 34' 19"
            pm_ra=-0.49768,  # -497.68 mas/yr
            pm_dec=-0.11467,  # -114.67 mas/yr
            parallax_mas=90.91,
            radial_km_per_s=-28.35,
        ),
        magnitude=2.14,
    ),
    StarCatalogEntry(
        id=SE_MARKAB,
        name="Markab",
        nomenclature="alPeg",
        hip_number=113963,
        data=StarData(
            ra_j2000=346.190223,  # 23h 04m 45.7s
            dec_j2000=15.205267,  # +15° 12' 19"
            pm_ra=0.06040,  # 60.40 mas/yr
            pm_dec=-0.04130,  # -41.30 mas/yr
            parallax_mas=24.46,
            radial_km_per_s=24.32,
        ),
        magnitude=2.49,
    ),
    StarCatalogEntry(
        id=SE_SCHEAT,
        name="Scheat",
        nomenclature="bePeg",
        hip_number=113881,
        data=StarData(
            ra_j2000=345.943514,  # 23h 03m 46.5s
            dec_j2000=28.082785,  # +28° 04' 58"
            pm_ra=0.18765,  # 187.65 mas/yr
            pm_dec=0.13693,  # 136.93 mas/yr
            parallax_mas=16.64,
            radial_km_per_s=35.66,
        ),
        magnitude=2.42,
    ),
    # ======== BEHENIAN FIXED STARS (additional) ========
    StarCatalogEntry(
        id=SE_ALCYONE,
        name="Alcyone",
        nomenclature="etTau",
        hip_number=17702,
        data=StarData(
            ra_j2000=56.871152,  # 03h 47m 29.1s (Eta Tauri, brightest Pleiad)
            dec_j2000=24.105136,  # +24° 06' 18"
            pm_ra=0.01934,  # 19.34 mas/yr
            pm_dec=-0.04367,  # -43.67 mas/yr
            parallax_mas=8.09,
            radial_km_per_s=34.09,
        ),
        magnitude=2.87,
    ),
    StarCatalogEntry(
        id=SE_ALGORAB,
        name="Algorab",
        nomenclature="deCrv",
        hip_number=60965,
        data=StarData(
            ra_j2000=187.466063,  # 12h 29m 51.9s (Delta Corvi)
            dec_j2000=-16.515431,  # -16° 30' 56"
            pm_ra=-0.21049,  # -210.49 mas/yr
            pm_dec=-0.13874,  # -138.74 mas/yr
            parallax_mas=37.55,
            radial_km_per_s=-18.15,
        ),
        magnitude=2.95,
    ),
    StarCatalogEntry(
        id=SE_ALPHECCA,
        name="Alphecca",
        nomenclature="alCrB",
        hip_number=76267,
        data=StarData(
            ra_j2000=233.671953,  # 15h 34m 41.3s (Alpha Coronae Borealis)
            dec_j2000=26.714693,  # +26° 42' 53"
            pm_ra=0.12027,  # 120.27 mas/yr
            pm_dec=-0.08958,  # -89.58 mas/yr
            parallax_mas=43.46,
            radial_km_per_s=-17.74,
        ),
        magnitude=2.23,
    ),
    StarCatalogEntry(
        id=SE_DENEB_ALGEDI,
        name="Deneb Algedi",
        nomenclature="deCap",
        hip_number=107556,
        data=StarData(
            ra_j2000=326.760184,  # 21h 47m 02.4s (Delta Capricorni)
            dec_j2000=-16.127287,  # -16° 07' 38"
            pm_ra=0.26170,  # 261.70 mas/yr (Hipparcos, van Leeuwen 2007)
            pm_dec=-0.29670,  # -296.70 mas/yr (Hipparcos, van Leeuwen 2007)
            parallax_mas=84.27,
            radial_km_per_s=16.22,
        ),
        magnitude=2.81,
    ),
    # ======== PLEIADES CLUSTER STARS ========
    # The Pleiades (M45) is an open star cluster in Taurus
    # These are the 9 brightest named stars visible to the naked eye
    StarCatalogEntry(
        id=SE_ASTEROPE,
        name="Asterope",
        nomenclature="21Tau",
        hip_number=17579,
        data=StarData(
            ra_j2000=56.476958,  # 03h 45m 54.5s (21 Tauri)
            dec_j2000=24.554722,  # +24° 33' 17"
            pm_ra=0.02018,  # 20.18 mas/yr
            pm_dec=-0.04487,  # -44.87 mas/yr
            parallax_mas=8.77,
            radial_km_per_s=21.45,
        ),
        magnitude=5.76,
    ),
    StarCatalogEntry(
        id=SE_CELAENO,
        name="Celaeno",
        nomenclature="16Tau",
        hip_number=17489,
        data=StarData(
            ra_j2000=56.200830,  # 03h 44m 48.2s (16 Tauri)
            dec_j2000=24.289389,  # +24° 17' 22"
            pm_ra=0.02038,  # 20.38 mas/yr
            pm_dec=-0.04481,  # -44.81 mas/yr
            parallax_mas=8.65,
            radial_km_per_s=24.78,
        ),
        magnitude=5.45,
    ),
    StarCatalogEntry(
        id=SE_ELECTRA,
        name="Electra",
        nomenclature="17Tau",
        hip_number=17499,
        data=StarData(
            ra_j2000=56.218908,  # 03h 44m 52.5s (17 Tauri)
            dec_j2000=24.113336,  # +24° 06' 48"
            pm_ra=0.02084,  # 20.84 mas/yr
            pm_dec=-0.04606,  # -46.06 mas/yr
            parallax_mas=8.06,
            radial_km_per_s=19.85,
        ),
        magnitude=3.70,
    ),
    StarCatalogEntry(
        id=SE_MAIA,
        name="Maia",
        nomenclature="20Tau",
        hip_number=17573,
        data=StarData(
            ra_j2000=56.456819,  # 03h 45m 49.6s (20 Tauri)
            dec_j2000=24.367750,  # +24° 22' 04"
            pm_ra=0.02095,  # 20.95 mas/yr
            pm_dec=-0.04598,  # -45.98 mas/yr
            parallax_mas=8.51,
            radial_km_per_s=33.1,
        ),
        magnitude=3.87,
    ),
    StarCatalogEntry(
        id=SE_MEROPE,
        name="Merope",
        nomenclature="23Tau",
        hip_number=17608,
        data=StarData(
            ra_j2000=56.581502,  # 03h 46m 19.6s (23 Tauri)
            dec_j2000=23.948353,  # +23° 56' 54"
            pm_ra=0.02113,  # 21.13 mas/yr
            pm_dec=-0.04365,  # -43.65 mas/yr
            parallax_mas=8.58,
            radial_km_per_s=23.12,
        ),
        magnitude=4.14,
    ),
    StarCatalogEntry(
        id=SE_TAYGETA,
        name="Taygeta",
        nomenclature="19Tau",
        hip_number=17531,
        data=StarData(
            ra_j2000=56.302063,  # 03h 45m 12.5s (19 Tauri)
            dec_j2000=24.467278,  # +24° 28' 02"
            pm_ra=0.02124,  # 21.24 mas/yr
            pm_dec=-0.04056,  # -40.56 mas/yr
            parallax_mas=7.97,
            radial_km_per_s=37.47,
        ),
        magnitude=4.30,
    ),
    StarCatalogEntry(
        id=SE_ATLAS,
        name="Atlas",
        nomenclature="27Tau",
        hip_number=17847,
        data=StarData(
            ra_j2000=57.290596,  # 03h 49m 09.7s (27 Tauri)
            dec_j2000=24.053417,  # +24° 03' 12"
            pm_ra=0.01770,  # 17.70 mas/yr
            pm_dec=-0.04418,  # -44.18 mas/yr
            parallax_mas=8.53,
            radial_km_per_s=25.28,
        ),
        magnitude=3.62,
    ),
    StarCatalogEntry(
        id=SE_PLEIONE,
        name="Pleione",
        nomenclature="28Tau",
        hip_number=17851,
        data=StarData(
            ra_j2000=57.296738,  # 03h 49m 11.2s (28 Tauri)
            dec_j2000=24.136750,  # +24° 08' 12"
            pm_ra=0.01807,  # 18.07 mas/yr
            pm_dec=-0.04720,  # -47.20 mas/yr
            parallax_mas=8.54,
            radial_km_per_s=34.14,
        ),
        magnitude=5.09,
    ),
    # ======== HYADES CLUSTER STARS ========
    # The Hyades is an open star cluster in Taurus, one of the nearest to Earth
    # These are the brightest named members visible to the naked eye
    StarCatalogEntry(
        id=SE_PRIMA_HYADUM,
        name="Prima Hyadum",
        nomenclature="gaTau",
        hip_number=20205,
        data=StarData(
            ra_j2000=64.948349,  # 04h 19m 47.6s (Gamma Tauri)
            dec_j2000=15.627643,  # +15° 37' 40"
            pm_ra=0.11546,  # 115.46 mas/yr
            pm_dec=-0.02342,  # -23.42 mas/yr
            parallax_mas=20.19,
            radial_km_per_s=54.6,
        ),
        magnitude=3.65,
    ),
    StarCatalogEntry(
        id=SE_SECUNDA_HYADUM,
        name="Secunda Hyadum",
        nomenclature="de1Tau",
        hip_number=20455,
        data=StarData(
            ra_j2000=65.733719,  # 04h 22m 56.1s (Delta^1 Tauri)
            dec_j2000=17.542514,  # +17° 32' 33"
            pm_ra=0.10656,  # 106.56 mas/yr
            pm_dec=-0.02918,  # -29.18 mas/yr
            parallax_mas=20.96,
            radial_km_per_s=58.87,
        ),
        magnitude=3.77,
    ),
    StarCatalogEntry(
        id=SE_THETA_TAURI,
        name="Theta Tauri",
        nomenclature="th2Tau",
        hip_number=20894,
        data=StarData(
            ra_j2000=67.165586,  # 04h 28m 39.7s (Theta^2 Tauri)
            dec_j2000=15.870882,  # +15° 52' 15"
            pm_ra=0.10842,  # 108.42 mas/yr
            pm_dec=-0.02674,  # -26.74 mas/yr
        ),
        magnitude=3.40,
    ),
    StarCatalogEntry(
        id=SE_AIN,
        name="Ain",
        nomenclature="epTau",
        hip_number=20889,
        data=StarData(
            ra_j2000=67.154163,  # 04h 28m 37.0s (Epsilon Tauri)
            dec_j2000=19.180560,  # +19° 10' 50"
            pm_ra=0.10619,  # 106.19 mas/yr
            pm_dec=-0.03784,  # -37.84 mas/yr
            parallax_mas=22.24,
            radial_km_per_s=53.5,
        ),
        magnitude=3.53,
    ),
    # ======== SOUTHERN CROSS (CRUX) CONSTELLATION ========
    # Completing the Crux constellation - Acrux, Mimosa, Gacrux already defined above
    StarCatalogEntry(
        id=SE_DELTA_CRUCIS,
        name="Delta Crucis",
        nomenclature="deCru",
        hip_number=59747,
        data=StarData(
            ra_j2000=183.786301,  # 12h 15m 08.7s (Delta Crucis)
            dec_j2000=-58.748927,  # -58° 44' 56"
            pm_ra=-0.03581,  # -35.81 mas/yr
            pm_dec=-0.01036,  # -10.36 mas/yr
        ),
        magnitude=2.80,
    ),
    # ======== CENTAURUS CONSTELLATION ========
    # Completing the bright stars of Centaurus (Alpha and Beta already defined above)
    StarCatalogEntry(
        id=SE_MENKENT,
        name="Menkent",
        nomenclature="thCen",
        hip_number=68933,
        data=StarData(
            ra_j2000=211.670528,  # 14h 06m 40.9s (Theta Centauri)
            dec_j2000=-36.369958,  # -36° 22' 12"
            pm_ra=-0.52053,  # -520.53 mas/yr
            pm_dec=-0.51806,  # -518.06 mas/yr
            parallax_mas=55.45,
            radial_km_per_s=-24.68,
        ),
        magnitude=2.06,
    ),
    StarCatalogEntry(
        id=SE_MUHLIFAIN,
        name="Muhlifain",
        nomenclature="gaCen",
        hip_number=61932,
        data=StarData(
            ra_j2000=190.379200,  # 12h 41m 31.0s (Gamma Centauri)
            dec_j2000=-48.959889,  # -48° 57' 36"
            pm_ra=-0.18572,  # -185.72 mas/yr
            pm_dec=0.00579,  # 5.79 mas/yr
            parallax_mas=25.06,
            radial_km_per_s=-28.46,
        ),
        magnitude=2.17,
    ),
    StarCatalogEntry(
        id=SE_EPSILON_CENTAURI,
        name="Epsilon Centauri",
        nomenclature="epCen",
        hip_number=66657,
        data=StarData(
            ra_j2000=204.971958,  # 13h 39m 53.3s (Epsilon Centauri)
            dec_j2000=-53.466389,  # -53° 27' 59"
            pm_ra=-0.01530,  # -15.30 mas/yr
            pm_dec=-0.01172,  # -11.72 mas/yr
        ),
        magnitude=2.30,
    ),
    StarCatalogEntry(
        id=SE_ETA_CENTAURI,
        name="Eta Centauri",
        nomenclature="etCen",
        hip_number=71352,
        data=StarData(
            ra_j2000=218.876841,  # 14h 35m 30.4s (Eta Centauri)
            dec_j2000=-42.157811,  # -42° 09' 28"
            pm_ra=-0.03473,  # -34.73 mas/yr
            pm_dec=-0.03272,  # -32.72 mas/yr
        ),
        magnitude=2.31,
    ),
    StarCatalogEntry(
        id=SE_ZETA_CENTAURI,
        name="Zeta Centauri",
        nomenclature="zeCen",
        hip_number=68002,
        data=StarData(
            ra_j2000=208.885225,  # 13h 55m 32.5s (Zeta Centauri)
            dec_j2000=-47.288375,  # -47° 17' 18"
            pm_ra=-0.05737,  # -57.37 mas/yr
            pm_dec=-0.04455,  # -44.55 mas/yr
        ),
        magnitude=2.55,
    ),
    # Scorpius constellation stars
    StarCatalogEntry(
        id=SE_SARGAS,
        name="Sargas",
        nomenclature="thSco",
        hip_number=86228,
        data=StarData(
            ra_j2000=264.329711,  # 17h 37m 19.1s
            dec_j2000=-42.997824,  # -42° 59' 52"
            pm_ra=0.00554,  # 5.54 mas/yr
            pm_dec=-0.00312,  # -3.12 mas/yr
            parallax_mas=10.86,
            radial_km_per_s=1.55,
        ),
        magnitude=1.87,
    ),
    StarCatalogEntry(
        id=SE_DSCHUBBA,
        name="Dschubba",
        nomenclature="deSco",
        hip_number=78401,
        data=StarData(
            ra_j2000=240.083359,  # 16h 00m 20.0s
            dec_j2000=-22.621710,  # -22° 37' 18"
            pm_ra=-0.01021,  # -10.21 mas/yr
            pm_dec=-0.03541,  # -35.41 mas/yr
            parallax_mas=6.64,
            radial_km_per_s=-24.62,
        ),
        magnitude=2.32,
    ),
    StarCatalogEntry(
        id=SE_GRAFFIAS,
        name="Graffias",
        nomenclature="beSco",
        hip_number=78820,
        data=StarData(
            ra_j2000=241.359296,  # 16h 05m 26.2s
            dec_j2000=-19.805453,  # -19° 48' 20"
            pm_ra=-0.00520,  # -5.20 mas/yr
            pm_dec=-0.02404,  # -24.04 mas/yr
            parallax_mas=8.07,
            radial_km_per_s=-22.73,
        ),
        magnitude=2.56,
    ),
    StarCatalogEntry(
        id=SE_LESATH,
        name="Lesath",
        nomenclature="upSco",
        hip_number=85696,
        data=StarData(
            ra_j2000=262.690901,  # 17h 30m 45.8s
            dec_j2000=-37.295811,  # -37° 17' 45"
            pm_ra=-0.00237,  # -2.37 mas/yr
            pm_dec=-0.03009,  # -30.09 mas/yr
            parallax_mas=5.66,
            radial_km_per_s=-18.89,
        ),
        magnitude=2.70,
    ),
    # Leo constellation stars
    StarCatalogEntry(
        id=SE_ZOSMA,
        name="Zosma",
        nomenclature="deLeo",
        hip_number=54872,
        data=StarData(
            ra_j2000=168.527089,  # 11h 14m 06.5s (Delta Leonis)
            dec_j2000=20.523611,  # +20° 31' 25"
            pm_ra=0.14342,  # 143.42 mas/yr
            pm_dec=-0.12988,  # -129.88 mas/yr
            parallax_mas=55.82,
            radial_km_per_s=-48.24,
        ),
        magnitude=2.56,
    ),
    # ======== ZODIACAL CONSTELLATION BRIGHT STARS ========
    # Stars from zodiacal constellations used in astrological interpretation
    # ======== ARIES CONSTELLATION ========
    # The Ram - first sign of the zodiac
    StarCatalogEntry(
        id=SE_HAMAL,
        name="Hamal",
        nomenclature="alAri",
        hip_number=9884,
        data=StarData(
            ra_j2000=31.793357,  # 02h 07m 10.4s (Alpha Arietis)
            dec_j2000=23.462418,  # +23° 27' 45"
            pm_ra=0.18855,  # 188.55 mas/yr
            pm_dec=-0.14808,  # -148.08 mas/yr
            parallax_mas=49.56,
            radial_km_per_s=11.2,
        ),
        magnitude=2.00,
    ),
    StarCatalogEntry(
        id=SE_SHERATAN,
        name="Sheratan",
        nomenclature="beAri",
        hip_number=8903,
        data=StarData(
            ra_j2000=28.660046,  # 01h 54m 38.4s (Beta Arietis)
            dec_j2000=20.808031,  # +20° 48' 29"
            pm_ra=0.09874,  # 98.74 mas/yr
            pm_dec=-0.11041,  # -110.41 mas/yr
            parallax_mas=55.6,
            radial_km_per_s=23.46,
        ),
        magnitude=2.64,
    ),
    StarCatalogEntry(
        id=SE_MESARTHIM,
        name="Mesarthim",
        nomenclature="gaAri",
        hip_number=8832,
        data=StarData(
            ra_j2000=28.382551,  # 01h 53m 31.8s (Gamma Arietis)
            dec_j2000=19.293852,  # +19° 17' 38"
            pm_ra=0.07920,  # 79.20 mas/yr
            pm_dec=-0.09763,  # -97.63 mas/yr
            parallax_mas=19.88,
            radial_km_per_s=23.69,
        ),
        magnitude=3.88,
    ),
    # ======== CANCER CONSTELLATION ========
    # The Crab - features the Beehive Cluster (M44)
    StarCatalogEntry(
        id=SE_ACUBENS,
        name="Acubens",
        nomenclature="alCnc",
        hip_number=44066,
        data=StarData(
            ra_j2000=134.621761,  # 08h 58m 29.2s (Alpha Cancri)
            dec_j2000=11.857700,  # +11° 51' 28"
            pm_ra=0.04323,  # 43.23 mas/yr (Hipparcos, van Leeuwen 2007)
            pm_dec=-0.02963,  # -29.63 mas/yr (Hipparcos, van Leeuwen 2007)
            parallax_mas=17.32,
            radial_km_per_s=-28.86,
        ),
        magnitude=4.25,
    ),
    StarCatalogEntry(
        id=SE_TARF,
        name="Tarf",
        nomenclature="beCnc",
        hip_number=42911,
        data=StarData(
            ra_j2000=130.821442,  # 08h 43m 17.1s (Beta Cancri)
            dec_j2000=9.185544,  # +09° 11' 08"
            pm_ra=-0.01767,  # -17.67 mas/yr
            pm_dec=-0.22926,  # -229.26 mas/yr
        ),
        magnitude=3.52,
    ),
    StarCatalogEntry(
        id=SE_ASELLUS_BOREALIS,
        name="Asellus Borealis",
        nomenclature="gaCnc",
        hip_number=42806,
        data=StarData(
            ra_j2000=130.821446,  # 08h 43m 17.1s (ICRS, Hipparcos HIP 42806)
            dec_j2000=21.468499,  # +21° 28' 06.6"
            pm_ra=-0.10351,  # -103.51 mas/yr (Hipparcos, van Leeuwen 2007)
            pm_dec=-0.03948,  # -39.48 mas/yr (Hipparcos, van Leeuwen 2007)
            parallax_mas=18.0,
            radial_km_per_s=16.56,
        ),
        magnitude=4.66,
    ),
    StarCatalogEntry(
        id=SE_ASELLUS_AUSTRALIS,
        name="Asellus Australis",
        nomenclature="deCnc",
        hip_number=42911,  # Corrected HIP number for Delta Cancri
        data=StarData(
            ra_j2000=131.171247,  # 08h 44m 41.1s (ICRS, Hipparcos HIP 42911)
            dec_j2000=18.154306,  # +18° 09' 15.5"
            pm_ra=-0.01767,  # -17.67 mas/yr (Hipparcos, van Leeuwen 2007)
            pm_dec=-0.22926,  # -229.26 mas/yr (Hipparcos, van Leeuwen 2007)
            parallax_mas=24.98,
            radial_km_per_s=1.72,
        ),
        magnitude=3.94,
    ),
    # ======== SAGITTARIUS CONSTELLATION ========
    # The Archer - prominent in the summer sky, contains galactic center
    StarCatalogEntry(
        id=SE_KAUS_AUSTRALIS,
        name="Kaus Australis",
        nomenclature="epSgr",
        hip_number=90185,
        data=StarData(
            ra_j2000=276.042993,  # 18h 24m 10.3s (Epsilon Sagittarii)
            dec_j2000=-34.384616,  # -34° 23' 05"
            pm_ra=-0.03942,  # -39.42 mas/yr
            pm_dec=-0.12420,  # -124.20 mas/yr
            parallax_mas=22.76,
            radial_km_per_s=-19.63,
        ),
        magnitude=1.85,
    ),
    StarCatalogEntry(
        id=SE_NUNKI,
        name="Nunki",
        nomenclature="siSgr",
        hip_number=92855,
        data=StarData(
            ra_j2000=283.816360,  # 18h 55m 15.9s (Sigma Sagittarii)
            dec_j2000=-26.296724,  # -26° 17' 48"
            pm_ra=0.01514,  # 15.14 mas/yr
            pm_dec=-0.05343,  # -53.43 mas/yr
            parallax_mas=14.32,
            radial_km_per_s=-10.42,
        ),
        magnitude=2.02,
    ),
    StarCatalogEntry(
        id=SE_KAUS_MEDIA,
        name="Kaus Media",
        nomenclature="deSgr",
        hip_number=89931,
        data=StarData(
            ra_j2000=275.248508,  # 18h 20m 59.6s (Delta Sagittarii)
            dec_j2000=-29.828104,  # -29° 49' 41"
            pm_ra=0.03254,  # 32.54 mas/yr
            pm_dec=-0.02557,  # -25.57 mas/yr
            parallax_mas=9.38,
            radial_km_per_s=-34.58,
        ),
        magnitude=2.70,
    ),
    StarCatalogEntry(
        id=SE_KAUS_BOREALIS,
        name="Kaus Borealis",
        nomenclature="laSgr",
        hip_number=90496,
        data=StarData(
            ra_j2000=276.992681,  # 18h 27m 58.2s (Lambda Sagittarii)
            dec_j2000=-25.421701,  # -25° 25' 18"
            pm_ra=-0.04476,  # -44.76 mas/yr
            pm_dec=-0.18566,  # -185.66 mas/yr
            parallax_mas=41.72,
            radial_km_per_s=-46.4,
        ),
        magnitude=2.81,
    ),
    StarCatalogEntry(
        id=SE_ASCELLA,
        name="Ascella",
        nomenclature="zeSgr",
        hip_number=93506,
        data=StarData(
            ra_j2000=285.653043,  # 19h 02m 36.7s (Zeta Sagittarii)
            dec_j2000=-29.880063,  # -29° 52' 48"
            pm_ra=0.01079,  # 10.79 mas/yr
            pm_dec=0.02111,  # 21.11 mas/yr
            parallax_mas=36.98,
            radial_km_per_s=26.64,
        ),
        magnitude=2.59,
    ),
    # ======== CAPRICORNUS CONSTELLATION ========
    # The Sea Goat - complementing Deneb Algedi already defined above
    StarCatalogEntry(
        id=SE_ALGEDI,
        name="Algedi",
        nomenclature="alCap",
        hip_number=100064,
        data=StarData(
            ra_j2000=304.513565,  # 20h 18m 03.3s (ICRS, Hipparcos HIP 100064)
            dec_j2000=-12.544852,  # -12° 32' 41.5"
            pm_ra=0.06263,  # 62.63 mas/yr (Hipparcos, van Leeuwen 2007)
            pm_dec=0.00266,  # 2.66 mas/yr (Hipparcos, van Leeuwen 2007)
            parallax_mas=5.73,
            radial_km_per_s=-21.67,
        ),
        magnitude=3.57,
    ),
    StarCatalogEntry(
        id=SE_DABIH,
        name="Dabih",
        nomenclature="beCap",
        hip_number=100345,
        data=StarData(
            ra_j2000=305.252803,  # 20h 21m 00.7s (Beta Capricorni)
            dec_j2000=-14.781405,  # -14° 46' 53"
            pm_ra=0.04492,  # 44.92 mas/yr
            pm_dec=0.00738,  # 7.38 mas/yr
            parallax_mas=9.98,
            radial_km_per_s=-5.74,
        ),
        magnitude=3.08,
    ),
    StarCatalogEntry(
        id=SE_NASHIRA,
        name="Nashira",
        nomenclature="gaCap",
        hip_number=106985,
        data=StarData(
            ra_j2000=325.022735,  # 21h 40m 05.5s (Gamma Capricorni)
            dec_j2000=-16.662308,  # -16° 39' 44"
            pm_ra=0.18756,  # 187.56 mas/yr
            pm_dec=-0.02245,  # -22.45 mas/yr
            parallax_mas=20.77,
            radial_km_per_s=-10.0,
        ),
        magnitude=3.68,
    ),
    # ======== AQUARIUS CONSTELLATION ========
    # The Water Bearer
    StarCatalogEntry(
        id=SE_SADALSUUD,
        name="Sadalsuud",
        nomenclature="beAqr",
        hip_number=106278,
        data=StarData(
            ra_j2000=322.889715,  # 21h 31m 33.5s (Beta Aquarii)
            dec_j2000=-5.571172,  # -05° 34' 16"
            pm_ra=0.01877,  # 18.77 mas/yr
            pm_dec=-0.00821,  # -8.21 mas/yr
            parallax_mas=6.07,
            radial_km_per_s=29.66,
        ),
        magnitude=2.87,
    ),
    StarCatalogEntry(
        id=SE_SADALMELIK,
        name="Sadalmelik",
        nomenclature="alAqr",
        hip_number=109074,
        data=StarData(
            ra_j2000=331.445983,  # 22h 05m 47.0s (Alpha Aquarii)
            dec_j2000=-0.319849,  # -00° 19' 11"
            pm_ra=0.01825,  # 18.25 mas/yr
            pm_dec=-0.00939,  # -9.39 mas/yr
            parallax_mas=6.23,
            radial_km_per_s=37.65,
        ),
        magnitude=2.96,
    ),
    StarCatalogEntry(
        id=SE_SKAT,
        name="Skat",
        nomenclature="deAqr",
        hip_number=113136,
        data=StarData(
            ra_j2000=343.662556,  # 22h 54m 39.0s (Delta Aquarii)
            dec_j2000=-15.820827,  # -15° 49' 15"
            pm_ra=-0.04260,  # -42.60 mas/yr (Hipparcos, van Leeuwen 2007)
            pm_dec=-0.02789,  # -27.89 mas/yr (Hipparcos, van Leeuwen 2007)
            parallax_mas=20.31,
            radial_km_per_s=44.84,
        ),
        magnitude=3.27,
    ),
    # ======== PISCES CONSTELLATION ========
    # The Fishes - where the vernal equinox currently resides
    StarCatalogEntry(
        id=SE_ETA_PISCIUM,
        name="Eta Piscium",
        nomenclature="etPsc",
        hip_number=5742,
        data=StarData(
            ra_j2000=18.437089,  # 01h 13m 44.9s (Eta Piscium)
            dec_j2000=15.345823,  # +15° 20' 45"
            pm_ra=0.01750,  # 17.50 mas/yr
            pm_dec=-0.02204,  # -22.04 mas/yr
        ),
        magnitude=3.62,
    ),
    StarCatalogEntry(
        id=SE_ALRESCHA,
        name="Alrescha",
        nomenclature="alPsc",
        hip_number=7097,
        data=StarData(
            ra_j2000=30.511749,  # 02h 02m 02.8s (ICRS, Hipparcos HIP 7097)
            dec_j2000=2.763761,  # +02° 45' 49.5"
            pm_ra=0.02714,  # 27.14 mas/yr (Hipparcos, van Leeuwen 2007)
            pm_dec=-0.00264,  # -2.64 mas/yr (Hipparcos, van Leeuwen 2007)
            parallax_mas=21.66,
            radial_km_per_s=31.77,
        ),
        magnitude=3.82,
    ),
    # ======== ANDROMEDA CONSTELLATION ========
    # Data: ESA Hipparcos catalogue (ESA SP-1200, 1997) via CDS/SIMBAD
    StarCatalogEntry(
        id=SE_ALPHERATZ,
        name="Alpheratz",
        nomenclature="alAnd",
        hip_number=677,
        data=StarData(
            ra_j2000=2.096916,  # 00h 08m 23.3s (Alpha Andromedae, HIP 677)
            dec_j2000=29.090431,  # +29° 05' 25.6"
            pm_ra=0.13746,  # 137.46 mas/yr (mu_alpha*, includes cos(dec))
            pm_dec=-0.16344,  # -163.44 mas/yr
            parallax_mas=33.62,
            radial_km_per_s=15.84,
        ),
        magnitude=2.06,
    ),
    # ======== PEGASUS CONSTELLATION - additional stars ========
    StarCatalogEntry(
        id=SE_ALGENIB,
        name="Algenib",
        nomenclature="gaPeg",
        hip_number=1067,
        data=StarData(
            ra_j2000=3.308968,  # 00h 13m 14.2s (Gamma Pegasi, HIP 1067)
            dec_j2000=15.183598,  # +15° 11' 01.0"
            pm_ra=0.00198,  # 1.98 mas/yr
            pm_dec=-0.00928,  # -9.28 mas/yr
            parallax_mas=8.33,
            radial_km_per_s=23.62,
        ),
        magnitude=2.84,
    ),
    # ======== GEMINI CONSTELLATION - additional stars ========
    # Data: ESA Hipparcos catalogue via CDS/SIMBAD
    StarCatalogEntry(
        id=SE_PROPUS,
        name="Propus",
        nomenclature="etGem",
        hip_number=29655,
        data=StarData(
            ra_j2000=93.719355,  # 06h 14m 52.6s (Eta Geminorum, HIP 29655)
            dec_j2000=22.506787,  # +22° 30' 24.4"
            pm_ra=-0.06246,  # -62.46 mas/yr
            pm_dec=-0.01212,  # -12.12 mas/yr
            parallax_mas=8.48,
            radial_km_per_s=26.52,
        ),
        magnitude=3.28,
    ),
    StarCatalogEntry(
        id=SE_TEJAT,
        name="Tejat",
        nomenclature="muGem",
        hip_number=30343,
        data=StarData(
            ra_j2000=95.740112,  # 06h 22m 57.6s (Mu Geminorum, HIP 30343)
            dec_j2000=22.513583,  # +22° 30' 48.9"
            pm_ra=0.05639,  # 56.39 mas/yr
            pm_dec=-0.11003,  # -110.03 mas/yr
            parallax_mas=14.08,
            radial_km_per_s=57.27,
        ),
        magnitude=2.87,
    ),
    StarCatalogEntry(
        id=SE_ALHENA,
        name="Alhena",
        nomenclature="gaGem",
        hip_number=31681,
        data=StarData(
            ra_j2000=99.427960,  # 06h 37m 42.7s (Gamma Geminorum, HIP 31681)
            dec_j2000=16.399280,  # +16° 23' 57.4"
            pm_ra=0.01381,  # 13.81 mas/yr
            pm_dec=-0.05496,  # -54.96 mas/yr
            parallax_mas=29.84,
            radial_km_per_s=-11.68,
        ),
        magnitude=1.92,
    ),
    StarCatalogEntry(
        id=SE_WASAT,
        name="Wasat",
        nomenclature="deGem",
        hip_number=35550,
        data=StarData(
            ra_j2000=110.030727,  # 07h 20m 07.4s (Delta Geminorum, HIP 35550)
            dec_j2000=21.982304,  # +21° 58' 56.3"
            pm_ra=-0.01513,  # -15.13 mas/yr
            pm_dec=-0.00979,  # -9.79 mas/yr
            parallax_mas=53.94,
            radial_km_per_s=-1.6,
        ),
        magnitude=3.53,
    ),
    # ======== CANIS MAJOR CONSTELLATION - additional stars ========
    StarCatalogEntry(
        id=SE_ADHARA,
        name="Adhara",
        nomenclature="epCMa",
        hip_number=33579,
        data=StarData(
            ra_j2000=104.656453,  # 06h 58m 37.5s (Epsilon CMa, HIP 33579)
            dec_j2000=-28.972086,  # -28° 58' 19.5"
            pm_ra=0.00324,  # 3.24 mas/yr
            pm_dec=0.00133,  # 1.33 mas/yr
            parallax_mas=8.05,
            radial_km_per_s=4.02,
        ),
        magnitude=1.50,
    ),
    StarCatalogEntry(
        id=SE_WEZEN,
        name="Wezen",
        nomenclature="deCMa",
        hip_number=34444,
        data=StarData(
            ra_j2000=107.097850,  # 07h 08m 23.5s (Delta CMa, HIP 34444)
            dec_j2000=-26.393200,  # -26° 23' 35.5"
            pm_ra=-0.00312,  # -3.12 mas/yr
            pm_dec=0.00331,  # 3.31 mas/yr
            parallax_mas=2.03,
            radial_km_per_s=53.82,
        ),
        magnitude=1.84,
    ),
    # ======== DRACO CONSTELLATION ========
    StarCatalogEntry(
        id=SE_THUBAN,
        name="Thuban",
        nomenclature="alDra",
        hip_number=68756,
        data=StarData(
            ra_j2000=211.097323,  # 14h 04m 23.4s (Alpha Draconis, HIP 68756)
            dec_j2000=64.375870,  # +64° 22' 33.1"
            pm_ra=-0.05634,  # -56.34 mas/yr
            pm_dec=0.01721,  # 17.21 mas/yr
            parallax_mas=10.76,
            radial_km_per_s=-21.48,
        ),
        magnitude=3.68,
    ),
    # ======== HERCULES CONSTELLATION ========
    StarCatalogEntry(
        id=SE_RASALGETHI,
        name="Rasalgethi",
        nomenclature="alHer",
        hip_number=84345,
        data=StarData(
            ra_j2000=258.661909,  # 17h 14m 38.9s (Alpha1 Herculis, HIP 84345)
            dec_j2000=14.390341,  # +14° 23' 25.2"
            pm_ra=-0.00732,  # -7.32 mas/yr
            pm_dec=0.03607,  # 36.07 mas/yr
            parallax_mas=9.07,
            radial_km_per_s=-35.22,
        ),
        magnitude=3.06,  # Variable (SRC), mean V from Bright Star Catalogue
    ),
    # ======== CYGNUS CONSTELLATION ========
    StarCatalogEntry(
        id=SE_ALBIREO,
        name="Albireo",
        nomenclature="beCyg",
        hip_number=95947,
        data=StarData(
            ra_j2000=292.680315,  # 19h 30m 43.3s (Beta1 Cygni, HIP 95947)
            dec_j2000=27.959674,  # +27° 57' 34.8"
            pm_ra=-0.00717,  # -7.17 mas/yr
            pm_dec=-0.00615,  # -6.15 mas/yr
            parallax_mas=7.51,
            radial_km_per_s=1.29,
        ),
        magnitude=3.08,
    ),
    # ======== ANDROMEDA CONSTELLATION (additional) ========
    StarCatalogEntry(
        id=SE_MIRACH,
        name="Mirach",
        nomenclature="beAnd",
        hip_number=5447,
        data=StarData(
            ra_j2000=17.433013,  # 01h 09m 43.9s (Beta Andromedae, HIP 5447)
            dec_j2000=35.620557,  # +35° 37' 14.0"
            pm_ra=0.17590,  # 175.90 mas/yr
            pm_dec=-0.11220,  # -112.20 mas/yr
            parallax_mas=16.52,
            radial_km_per_s=28.07,
        ),
        magnitude=2.05,
    ),
    StarCatalogEntry(
        id=SE_ALMACH,
        name="Almach",
        nomenclature="ga1And",
        hip_number=9640,
        data=StarData(
            ra_j2000=30.974803,  # 02h 03m 53.9s (Gamma1 Andromedae, HIP 9640)
            dec_j2000=42.329725,  # +42° 19' 47.0"
            pm_ra=0.04232,  # 42.32 mas/yr
            pm_dec=-0.04930,  # -49.30 mas/yr
            parallax_mas=9.19,
            radial_km_per_s=-1.9,
        ),
        magnitude=2.10,
    ),
    # ======== CETUS CONSTELLATION ========
    StarCatalogEntry(
        id=SE_MENKAR,
        name="Menkar",
        nomenclature="alCet",
        hip_number=14135,
        data=StarData(
            ra_j2000=45.569885,  # 03h 02m 16.8s (Alpha Ceti, HIP 14135)
            dec_j2000=4.089737,  # +04° 05' 23.1"
            pm_ra=-0.01041,  # -10.41 mas/yr
            pm_dec=-0.07685,  # -76.85 mas/yr
            parallax_mas=5.66,
            radial_km_per_s=35.32,
        ),
        magnitude=2.54,
    ),
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
    "OMICRON": "om",
    "PI": "pi",
    "RHO": "rh",
    "SIGMA": "si",
    "TAU": "ta",
    "UPSILON": "up",
    "PHI": "ph",
    "CHI": "ch",
    "PSI": "ps",
    "OMEGA": "om",  # Note: same abbrev as omicron in some catalogs
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
    "COR LEONIS": SE_REGULUS,
    "ALPHA LEONIS": SE_REGULUS,
    "ALPHA LEO": SE_REGULUS,
    "87 LEO": SE_REGULUS,
    "α LEO": SE_REGULUS,
    "QALB AL-ASAD": SE_REGULUS,
    "BASILISKOS": SE_REGULUS,
    "REX": SE_REGULUS,
    "KALB AL-ASAD": SE_REGULUS,
    "ALLEO": SE_REGULUS,
    "WATCHER OF THE NORTH": SE_REGULUS,
    "VENANT": SE_REGULUS,
    # ======== SPICA (Alpha Virginis) ========
    "ALPHA VIRGINIS": SE_SPICA_STAR,
    "ALPHA VIR": SE_SPICA_STAR,
    "67 VIR": SE_SPICA_STAR,
    "α VIR": SE_SPICA_STAR,
    "AZIMECH": SE_SPICA_STAR,
    "ALARAPH": SE_SPICA_STAR,
    "ALVIR": SE_SPICA_STAR,
    "SUNBULA": SE_SPICA_STAR,
    "VIRGIN'S SPIKE": SE_SPICA_STAR,
    "ARISTA": SE_SPICA_STAR,
    # ======== ALGOL (Beta Persei) ========
    "DEMON STAR": SE_ALGOL,
    "BETA PERSEI": SE_ALGOL,
    "BETA PER": SE_ALGOL,
    "26 PER": SE_ALGOL,
    "β PER": SE_ALGOL,
    "GORGONEA PRIMA": SE_ALGOL,
    "RA'S AL-GHUL": SE_ALGOL,
    "RAS AL-GHUL": SE_ALGOL,
    "BEPER": SE_ALGOL,
    "HEAD OF THE GHOUL": SE_ALGOL,
    "GORGON'S HEAD": SE_ALGOL,
    # ======== SIRIUS (Alpha Canis Majoris) ========
    "DOG STAR": SE_SIRIUS,
    "ALPHA CANIS MAJORIS": SE_SIRIUS,
    "ALPHA CMA": SE_SIRIUS,
    "α CMA": SE_SIRIUS,
    "9 CMA": SE_SIRIUS,
    "CANICULA": SE_SIRIUS,
    "ASCHERE": SE_SIRIUS,
    "ALCMA": SE_SIRIUS,
    "AL-SHIRA": SE_SIRIUS,
    "SOTHIS": SE_SIRIUS,
    # ======== ALDEBARAN (Alpha Tauri) - ROYAL STAR ========
    "EYE OF TAURUS": SE_ALDEBARAN,
    "ALPHA TAURI": SE_ALDEBARAN,
    "ALPHA TAU": SE_ALDEBARAN,
    "87 TAU": SE_ALDEBARAN,
    "α TAU": SE_ALDEBARAN,
    "ALTAU": SE_ALDEBARAN,
    "PARILICIUM": SE_ALDEBARAN,
    "AL-DABARAN": SE_ALDEBARAN,
    "FOLLOWER": SE_ALDEBARAN,
    "ROHINI": SE_ALDEBARAN,
    "WATCHER OF THE EAST": SE_ALDEBARAN,
    "TASCHETER": SE_ALDEBARAN,
    # ======== ANTARES (Alpha Scorpii) - ROYAL STAR ========
    "RIVAL OF MARS": SE_ANTARES,
    "ALPHA SCORPII": SE_ANTARES,
    "ALPHA SCO": SE_ANTARES,
    "21 SCO": SE_ANTARES,
    "α SCO": SE_ANTARES,
    "ALSCO": SE_ANTARES,
    "COR SCORPII": SE_ANTARES,
    "CALB AL-AKRAB": SE_ANTARES,
    "HEART OF SCORPION": SE_ANTARES,
    "JYESHTHA": SE_ANTARES,
    "WATCHER OF THE WEST": SE_ANTARES,
    "SATEVIS": SE_ANTARES,
    # ======== VEGA (Alpha Lyrae) ========
    "HARP STAR": SE_VEGA,
    "ALPHA LYRAE": SE_VEGA,
    "ALPHA LYR": SE_VEGA,
    "3 LYR": SE_VEGA,
    "α LYR": SE_VEGA,
    "ALLYR": SE_VEGA,
    "WEGA": SE_VEGA,
    "AL-NASR AL-WAQI": SE_VEGA,
    "FIDIS": SE_VEGA,
    "ABHIJIT": SE_VEGA,
    # ======== POLARIS (Alpha Ursae Minoris) ========
    "NORTH STAR": SE_POLARIS,
    "POLE STAR": SE_POLARIS,
    "ALPHA URSAE MINORIS": SE_POLARIS,
    "ALPHA UMI": SE_POLARIS,
    "1 UMI": SE_POLARIS,
    "α UMI": SE_POLARIS,
    "ALUMI": SE_POLARIS,
    "CYNOSURA": SE_POLARIS,
    "LODESTAR": SE_POLARIS,
    "STELLA POLARIS": SE_POLARIS,
    # ======== FOMALHAUT (Alpha Piscis Austrini) - ROYAL STAR ========
    "FISH'S MOUTH": SE_FOMALHAUT,
    "ALPHA PISCIS AUSTRINI": SE_FOMALHAUT,
    "ALPHA PSA": SE_FOMALHAUT,
    "24 PSA": SE_FOMALHAUT,
    "α PSA": SE_FOMALHAUT,
    "ALPSA": SE_FOMALHAUT,
    "FUM AL-HUT": SE_FOMALHAUT,
    "OS PISCIS MERIDIANI": SE_FOMALHAUT,
    "LONELY STAR": SE_FOMALHAUT,
    "HASTORANG": SE_FOMALHAUT,
    "WATCHER OF THE SOUTH": SE_FOMALHAUT,
    # ======== BETELGEUSE (Alpha Orionis) ========
    "ARMPIT OF ORION": SE_BETELGEUSE,
    "ALPHA ORIONIS": SE_BETELGEUSE,
    "ALPHA ORI": SE_BETELGEUSE,
    "58 ORI": SE_BETELGEUSE,
    "α ORI": SE_BETELGEUSE,
    "ALORI": SE_BETELGEUSE,
    "BETELGEUZE": SE_BETELGEUSE,
    "IBT AL-JAUZAH": SE_BETELGEUSE,
    "ARDRA": SE_BETELGEUSE,
    "YAD AL-JAWZA": SE_BETELGEUSE,
    # ======== RIGEL (Beta Orionis) ========
    "FOOT OF ORION": SE_RIGEL,
    "BETA ORIONIS": SE_RIGEL,
    "BETA ORI": SE_RIGEL,
    "19 ORI": SE_RIGEL,
    "β ORI": SE_RIGEL,
    "BEORI": SE_RIGEL,
    "ALGEBAR": SE_RIGEL,
    "RIJL JAUZAH AL-YUSRA": SE_RIGEL,
    "ORION'S FOOT": SE_RIGEL,
    # ======== PROCYON (Alpha Canis Minoris) ========
    "LITTLE DOG STAR": SE_PROCYON,
    "ALPHA CANIS MINORIS": SE_PROCYON,
    "ALPHA CMI": SE_PROCYON,
    "10 CMI": SE_PROCYON,
    "α CMI": SE_PROCYON,
    "ALCMI": SE_PROCYON,
    "ANTECANIS": SE_PROCYON,
    "ELGOMAISA": SE_PROCYON,
    "AL-GHUMAISA": SE_PROCYON,
    # ======== CAPELLA (Alpha Aurigae) ========
    "SHE-GOAT": SE_CAPELLA,
    "ALPHA AURIGAE": SE_CAPELLA,
    "ALPHA AUR": SE_CAPELLA,
    "13 AUR": SE_CAPELLA,
    "α AUR": SE_CAPELLA,
    "ALAUR": SE_CAPELLA,
    "ALHAJOTH": SE_CAPELLA,
    "AMALTHEA": SE_CAPELLA,
    "GOAT STAR": SE_CAPELLA,
    # ======== ARCTURUS (Alpha Bootis) ========
    "BEAR GUARD": SE_ARCTURUS,
    "ALPHA BOOTIS": SE_ARCTURUS,
    "ALPHA BOO": SE_ARCTURUS,
    "16 BOO": SE_ARCTURUS,
    "α BOO": SE_ARCTURUS,
    "ALBOO": SE_ARCTURUS,
    "AL-SIMAK AL-RAMIH": SE_ARCTURUS,
    "HARIS AL-SAMA": SE_ARCTURUS,
    "GUARDIAN OF BEAR": SE_ARCTURUS,
    # ======== DENEB (Alpha Cygni) ========
    "TAIL OF HEN": SE_DENEB,
    "ALPHA CYGNI": SE_DENEB,
    "ALPHA CYG": SE_DENEB,
    "50 CYG": SE_DENEB,
    "α CYG": SE_DENEB,
    "ALCYG": SE_DENEB,
    "DHANAB AD-DAJAJAH": SE_DENEB,
    "ARIDED": SE_DENEB,
    "GALLINA": SE_DENEB,
    # ======== POLLUX (Beta Geminorum) ========
    "TWIN STAR": SE_POLLUX,
    "BETA GEMINORUM": SE_POLLUX,
    "BETA GEM": SE_POLLUX,
    "78 GEM": SE_POLLUX,
    "β GEM": SE_POLLUX,
    "BEGEM": SE_POLLUX,
    "POLYDEUCES": SE_POLLUX,
    "HEAD OF SECOND TWIN": SE_POLLUX,
    # ======== CASTOR (Alpha Geminorum) ========
    "ALPHA GEMINORUM": SE_CASTOR,
    "ALPHA GEM": SE_CASTOR,
    "66 GEM": SE_CASTOR,
    "α GEM": SE_CASTOR,
    "ALGEM": SE_CASTOR,
    "APOLLO": SE_CASTOR,
    "HEAD OF FIRST TWIN": SE_CASTOR,
    # ======== ALTAIR (Alpha Aquilae) ========
    "FLYING EAGLE": SE_ALTAIR,
    "ALPHA AQUILAE": SE_ALTAIR,
    "ALPHA AQL": SE_ALTAIR,
    "53 AQL": SE_ALTAIR,
    "α AQL": SE_ALTAIR,
    "ALAQL": SE_ALTAIR,
    "AL-NASR AL-TAIR": SE_ALTAIR,
    "ATAIR": SE_ALTAIR,
    "SRAVANA": SE_ALTAIR,
    # ======== ACHERNAR (Alpha Eridani) ========
    "END OF RIVER": SE_ACHERNAR,
    "ALPHA ERIDANI": SE_ACHERNAR,
    "ALPHA ERI": SE_ACHERNAR,
    "α ERI": SE_ACHERNAR,
    "ALERI": SE_ACHERNAR,
    "AKHIR AN-NAHR": SE_ACHERNAR,
    "RIVER'S END": SE_ACHERNAR,
    # ======== CANOPUS (Alpha Carinae) ========
    "SHIP'S PILOT": SE_CANOPUS,
    "ALPHA CARINAE": SE_CANOPUS,
    "ALPHA CAR": SE_CANOPUS,
    "α CAR": SE_CANOPUS,
    "ALCAR": SE_CANOPUS,
    "SUHAIL": SE_CANOPUS,
    "SUHAYL": SE_CANOPUS,
    "AGASTYA": SE_CANOPUS,
    # ======== ACRUX (Alpha Crucis) ========
    "ALPHA CRUCIS": SE_ACRUX,
    "ALPHA CRU": SE_ACRUX,
    "α CRU": SE_ACRUX,
    "ALCRU": SE_ACRUX,
    "CRUX ALPHA": SE_ACRUX,
    "STAR OF BETHLEHEM": SE_ACRUX,
    # ======== MIMOSA (Beta Crucis) ========
    "BETA CRUCIS": SE_MIMOSA,
    "BETA CRU": SE_MIMOSA,
    "β CRU": SE_MIMOSA,
    "BECRU": SE_MIMOSA,
    "BECRUX": SE_MIMOSA,
    "CRUX BETA": SE_MIMOSA,
    # ======== GACRUX (Gamma Crucis) ========
    "GAMMA CRUCIS": SE_GACRUX,
    "GAMMA CRU": SE_GACRUX,
    "γ CRU": SE_GACRUX,
    "GACRU": SE_GACRUX,
    "CRUX GAMMA": SE_GACRUX,
    "RUBIDEA": SE_GACRUX,
    # ======== DELTA CRUCIS (Delta Crucis) ========
    "DELTA CRUCIS": SE_DELTA_CRUCIS,
    "DELTA CRU": SE_DELTA_CRUCIS,
    "δ CRU": SE_DELTA_CRUCIS,
    "DECRU": SE_DELTA_CRUCIS,
    "CRUX DELTA": SE_DELTA_CRUCIS,
    "DECRUX": SE_DELTA_CRUCIS,
    # ======== HADAR (Beta Centauri) ========
    "BETA CENTAURI": SE_HADAR,
    "BETA CEN": SE_HADAR,
    "β CEN": SE_HADAR,
    "BECEN": SE_HADAR,
    "AGENA": SE_HADAR,
    "KNEE OF CENTAUR": SE_HADAR,
    # ======== RIGIL KENTAURUS (Alpha Centauri) ========
    "ALPHA CENTAURI": SE_RIGIL_KENT,
    "ALPHA CEN": SE_RIGIL_KENT,
    "α CEN": SE_RIGIL_KENT,
    "ALCEN": SE_RIGIL_KENT,
    "TOLIMAN": SE_RIGIL_KENT,
    "RIGIL KENT": SE_RIGIL_KENT,
    "FOOT OF CENTAUR": SE_RIGIL_KENT,
    "BUNGULA": SE_RIGIL_KENT,
    # ======== MENKENT (Theta Centauri) ========
    "THETA CENTAURI": SE_MENKENT,
    "THETA CEN": SE_MENKENT,
    "θ CEN": SE_MENKENT,
    "THCEN": SE_MENKENT,
    "HARATAN": SE_MENKENT,
    # ======== MUHLIFAIN (Gamma Centauri) ========
    "GAMMA CENTAURI": SE_MUHLIFAIN,
    "GAMMA CEN": SE_MUHLIFAIN,
    "γ CEN": SE_MUHLIFAIN,
    "GACEN": SE_MUHLIFAIN,
    # ======== EPSILON CENTAURI ========
    "EPSILON CENTAURI": SE_EPSILON_CENTAURI,
    "EPSILON CEN": SE_EPSILON_CENTAURI,
    "ε CEN": SE_EPSILON_CENTAURI,
    "EPCEN": SE_EPSILON_CENTAURI,
    # ======== ETA CENTAURI ========
    "ETA CENTAURI": SE_ETA_CENTAURI,
    "ETA CEN": SE_ETA_CENTAURI,
    "η CEN": SE_ETA_CENTAURI,
    "ETCEN": SE_ETA_CENTAURI,
    # ======== ZETA CENTAURI ========
    "ZETA CENTAURI": SE_ZETA_CENTAURI,
    "ZETA CEN": SE_ZETA_CENTAURI,
    "ζ CEN": SE_ZETA_CENTAURI,
    "ZECEN": SE_ZETA_CENTAURI,
    "ALNAIR": SE_ZETA_CENTAURI,
    # ======== SHAULA (Lambda Scorpii) ========
    "LAMBDA SCORPII": SE_SHAULA,
    "LAMBDA SCO": SE_SHAULA,
    "λ SCO": SE_SHAULA,
    "LASCO": SE_SHAULA,
    "SCORPION'S STING": SE_SHAULA,
    "35 SCO": SE_SHAULA,
    # ======== BELLATRIX (Gamma Orionis) ========
    "AMAZON STAR": SE_BELLATRIX,
    "GAMMA ORIONIS": SE_BELLATRIX,
    "GAMMA ORI": SE_BELLATRIX,
    "24 ORI": SE_BELLATRIX,
    "γ ORI": SE_BELLATRIX,
    "GAORI": SE_BELLATRIX,
    "FEMALE WARRIOR": SE_BELLATRIX,
    # ======== ELNATH (Beta Tauri) ========
    "BETA TAURI": SE_ELNATH,
    "BETA TAU": SE_ELNATH,
    "112 TAU": SE_ELNATH,
    "β TAU": SE_ELNATH,
    "BETAU": SE_ELNATH,
    "NATH": SE_ELNATH,
    "AL-NATH": SE_ELNATH,
    "EL NATH": SE_ELNATH,
    # ======== MIRA (Omicron Ceti) ========
    "WONDERFUL STAR": SE_MIRA,
    "OMICRON CETI": SE_MIRA,
    "OMICRON CET": SE_MIRA,
    "68 CET": SE_MIRA,
    "ο CET": SE_MIRA,
    "OMCET": SE_MIRA,
    "STELLA MIRA": SE_MIRA,
    # ======== ALNILAM (Epsilon Orionis) ========
    "STRING OF PEARLS": SE_ALNILAM,
    "EPSILON ORIONIS": SE_ALNILAM,
    "EPSILON ORI": SE_ALNILAM,
    "46 ORI": SE_ALNILAM,
    "ε ORI": SE_ALNILAM,
    "EPORI": SE_ALNILAM,
    "AL-NIZAM": SE_ALNILAM,
    # ======== ALNITAK (Zeta Orionis) ========
    "GIRDLE": SE_ALNITAK,
    "ZETA ORIONIS": SE_ALNITAK,
    "ZETA ORI": SE_ALNITAK,
    "50 ORI": SE_ALNITAK,
    "ζ ORI": SE_ALNITAK,
    "ZEORI": SE_ALNITAK,
    "AL-NITAK": SE_ALNITAK,
    # ======== MINTAKA (Delta Orionis) ========
    "BELT STAR": SE_MINTAKA,
    "DELTA ORIONIS": SE_MINTAKA,
    "DELTA ORI": SE_MINTAKA,
    "34 ORI": SE_MINTAKA,
    "δ ORI": SE_MINTAKA,
    "DEORI": SE_MINTAKA,
    "MINTAKA": SE_MINTAKA,
    # ======== SAIPH (Kappa Orionis) ========
    "SWORD OF GIANT": SE_SAIPH,
    "KAPPA ORIONIS": SE_SAIPH,
    "KAPPA ORI": SE_SAIPH,
    "53 ORI": SE_SAIPH,
    "κ ORI": SE_SAIPH,
    "KAORI": SE_SAIPH,
    "SAIF AL-JABBAR": SE_SAIPH,
    # ======== MEISSA (Lambda Orionis) ========
    "LAMBDA ORIONIS": SE_MEISSA,
    "LAMBDA ORI": SE_MEISSA,
    "39 ORI": SE_MEISSA,
    "λ ORI": SE_MEISSA,
    "LAORI": SE_MEISSA,
    "HEKA": SE_MEISSA,
    "HEAD OF ORION": SE_MEISSA,
    "AL-MAISAN": SE_MEISSA,
    # ======== DIPHDA (Beta Ceti) ========
    "FROG": SE_DIPHDA,
    "BETA CETI": SE_DIPHDA,
    "BETA CET": SE_DIPHDA,
    "16 CET": SE_DIPHDA,
    "β CET": SE_DIPHDA,
    "BECET": SE_DIPHDA,
    "DENEB KAITOS": SE_DIPHDA,
    "DIFDA AL-THANI": SE_DIPHDA,
    # ======== ALPHARD (Alpha Hydrae) ========
    "SOLITARY ONE": SE_ALPHARD,
    "ALPHA HYDRAE": SE_ALPHARD,
    "ALPHA HYA": SE_ALPHARD,
    "30 HYA": SE_ALPHARD,
    "α HYA": SE_ALPHARD,
    "ALHYA": SE_ALPHARD,
    "COR HYDRAE": SE_ALPHARD,
    "AL-FARD": SE_ALPHARD,
    # ======== RASALHAGUE (Alpha Ophiuchi) ========
    "HEAD OF SERPENT HOLDER": SE_RASALHAGUE,
    "ALPHA OPHIUCHI": SE_RASALHAGUE,
    "ALPHA OPH": SE_RASALHAGUE,
    "55 OPH": SE_RASALHAGUE,
    "α OPH": SE_RASALHAGUE,
    "ALOPH": SE_RASALHAGUE,
    "RAS AL-HAWWA": SE_RASALHAGUE,
    # ======== ETAMIN (Gamma Draconis) ========
    "DRAGON'S HEAD": SE_ETAMIN,
    "GAMMA DRACONIS": SE_ETAMIN,
    "GAMMA DRA": SE_ETAMIN,
    "33 DRA": SE_ETAMIN,
    "γ DRA": SE_ETAMIN,
    "GADORA": SE_ETAMIN,
    "ELTANIN": SE_ETAMIN,
    "AL-TINNIN": SE_ETAMIN,
    # ======== KOCHAB (Beta Ursae Minoris) ========
    "BETA URSAE MINORIS": SE_KOCHAB,
    "BETA UMI": SE_KOCHAB,
    "7 UMI": SE_KOCHAB,
    "β UMI": SE_KOCHAB,
    "BEUMI": SE_KOCHAB,
    "KAUKAB": SE_KOCHAB,
    # ======== ALKAID (Eta Ursae Majoris) ========
    "END OF TAIL": SE_ALKAID,
    "ETA URSAE MAJORIS": SE_ALKAID,
    "ETA UMA": SE_ALKAID,
    "85 UMA": SE_ALKAID,
    "η UMA": SE_ALKAID,
    "ETUMA": SE_ALKAID,
    "BENETNASH": SE_ALKAID,
    "AL-QA'ID": SE_ALKAID,
    # ======== DUBHE (Alpha Ursae Majoris) ========
    "BEAR'S BACK": SE_DUBHE,
    "ALPHA URSAE MAJORIS": SE_DUBHE,
    "ALPHA UMA": SE_DUBHE,
    "50 UMA": SE_DUBHE,
    "α UMA": SE_DUBHE,
    "ALUMA": SE_DUBHE,
    "THAHR AL-DUBB AL-AKBAR": SE_DUBHE,
    # ======== MERAK (Beta Ursae Majoris) ========
    "BETA URSAE MAJORIS": SE_MERAK,
    "BETA UMA": SE_MERAK,
    "48 UMA": SE_MERAK,
    "β UMA": SE_MERAK,
    "BEUMA": SE_MERAK,
    "AL-MARAKK": SE_MERAK,
    # ======== ALIOTH (Epsilon Ursae Majoris) ========
    "EPSILON URSAE MAJORIS": SE_ALIOTH,
    "EPSILON UMA": SE_ALIOTH,
    "77 UMA": SE_ALIOTH,
    "ε UMA": SE_ALIOTH,
    "EPUMA": SE_ALIOTH,
    "ALIATH": SE_ALIOTH,
    # ======== MIZAR (Zeta Ursae Majoris) ========
    "ZETA URSAE MAJORIS": SE_MIZAR,
    "ZETA UMA": SE_MIZAR,
    "79 UMA": SE_MIZAR,
    "ζ UMA": SE_MIZAR,
    "ZEUMA": SE_MIZAR,
    "HORSE AND RIDER": SE_MIZAR,
    # ======== ALCOR (80 Ursae Majoris) ========
    "80 URSAE MAJORIS": SE_ALCOR,
    "80 UMA": SE_ALCOR,
    "G UMA": SE_ALCOR,
    "SAIDAK": SE_ALCOR,
    "SUHA": SE_ALCOR,
    "AL-SAHJA": SE_ALCOR,
    # ======== PHECDA (Gamma Ursae Majoris) ========
    "GAMMA URSAE MAJORIS": SE_PHECDA,
    "GAMMA UMA": SE_PHECDA,
    "64 UMA": SE_PHECDA,
    "γ UMA": SE_PHECDA,
    "GAUMA": SE_PHECDA,
    "PHAD": SE_PHECDA,
    "PHEKDA": SE_PHECDA,
    "PHACD": SE_PHECDA,
    # ======== MEGREZ (Delta Ursae Majoris) ========
    "DELTA URSAE MAJORIS": SE_MEGREZ,
    "DELTA UMA": SE_MEGREZ,
    "69 UMA": SE_MEGREZ,
    "δ UMA": SE_MEGREZ,
    "DEUMA": SE_MEGREZ,
    "KAFFA": SE_MEGREZ,
    # ======== VINDEMIATRIX (Epsilon Virginis) ========
    "GRAPE GATHERER": SE_VINDEMIATRIX,
    "EPSILON VIRGINIS": SE_VINDEMIATRIX,
    "EPSILON VIR": SE_VINDEMIATRIX,
    "47 VIR": SE_VINDEMIATRIX,
    "ε VIR": SE_VINDEMIATRIX,
    "EPVIR": SE_VINDEMIATRIX,
    "ALMUREDIN": SE_VINDEMIATRIX,
    # ======== ZUBENELGENUBI (Alpha Librae) ========
    "SOUTHERN CLAW": SE_ZUBENELGENUBI,
    "ALPHA LIBRAE": SE_ZUBENELGENUBI,
    "ALPHA LIB": SE_ZUBENELGENUBI,
    "9 LIB": SE_ZUBENELGENUBI,
    "α LIB": SE_ZUBENELGENUBI,
    "ALLIB": SE_ZUBENELGENUBI,
    "KIFFA AUSTRALIS": SE_ZUBENELGENUBI,
    # ======== ZUBENESCHAMALI (Beta Librae) ========
    "NORTHERN CLAW": SE_ZUBENESCHAMALI,
    "BETA LIBRAE": SE_ZUBENESCHAMALI,
    "BETA LIB": SE_ZUBENESCHAMALI,
    "27 LIB": SE_ZUBENESCHAMALI,
    "β LIB": SE_ZUBENESCHAMALI,
    "BELIB": SE_ZUBENESCHAMALI,
    "KIFFA BOREALIS": SE_ZUBENESCHAMALI,
    # ======== UNUKALHAI (Alpha Serpentis) ========
    "SERPENT'S NECK": SE_UNUKALHAI,
    "ALPHA SERPENTIS": SE_UNUKALHAI,
    "ALPHA SER": SE_UNUKALHAI,
    "24 SER": SE_UNUKALHAI,
    "α SER": SE_UNUKALHAI,
    "ALSER": SE_UNUKALHAI,
    "COR SERPENTIS": SE_UNUKALHAI,
    # ======== ALGIEBA (Gamma Leonis) ========
    "LION'S MANE": SE_ALGIEBA,
    "GAMMA LEONIS": SE_ALGIEBA,
    "GAMMA LEO": SE_ALGIEBA,
    "41 LEO": SE_ALGIEBA,
    "γ LEO": SE_ALGIEBA,
    "GALEO": SE_ALGIEBA,
    "AL-JABHAH": SE_ALGIEBA,
    # ======== DENEBOLA (Beta Leonis) ========
    "LION'S TAIL": SE_DENEBOLA,
    "BETA LEONIS": SE_DENEBOLA,
    "BETA LEO": SE_DENEBOLA,
    "94 LEO": SE_DENEBOLA,
    "β LEO": SE_DENEBOLA,
    "BELEO": SE_DENEBOLA,
    "DHANAB AL-ASAD": SE_DENEBOLA,
    # ======== MARKAB (Alpha Pegasi) ========
    "SADDLE": SE_MARKAB,
    "ALPHA PEGASI": SE_MARKAB,
    "ALPHA PEG": SE_MARKAB,
    "54 PEG": SE_MARKAB,
    "α PEG": SE_MARKAB,
    "ALPEG": SE_MARKAB,
    "MANKIB AL-FARAS": SE_MARKAB,
    # ======== SCHEAT (Beta Pegasi) ========
    "LEG": SE_SCHEAT,
    "BETA PEGASI": SE_SCHEAT,
    "BETA PEG": SE_SCHEAT,
    "53 PEG": SE_SCHEAT,
    "β PEG": SE_SCHEAT,
    "BEPEG": SE_SCHEAT,
    "SAQ AL-FARAS": SE_SCHEAT,
    # ======== ALCYONE (Eta Tauri - Pleiades) - BEHENIAN ========
    "PLEIADES": SE_ALCYONE,
    "ETA TAURI": SE_ALCYONE,
    "ETA TAU": SE_ALCYONE,
    "25 TAU": SE_ALCYONE,
    "η TAU": SE_ALCYONE,
    "ETTAU": SE_ALCYONE,
    "SEVEN SISTERS": SE_ALCYONE,
    "KIMAH": SE_ALCYONE,
    "AL-THURAYYA": SE_ALCYONE,
    # ======== ALGORAB (Delta Corvi) - BEHENIAN ========
    "CROW'S WING": SE_ALGORAB,
    "DELTA CORVI": SE_ALGORAB,
    "DELTA CRV": SE_ALGORAB,
    "7 CRV": SE_ALGORAB,
    "δ CRV": SE_ALGORAB,
    "DECRV": SE_ALGORAB,
    "AL-GHIRAB": SE_ALGORAB,
    "GIENAH CORVI": SE_ALGORAB,
    # ======== ALPHECCA (Alpha Coronae Borealis) - BEHENIAN ========
    "GEMMA": SE_ALPHECCA,
    "ALPHA CORONAE BOREALIS": SE_ALPHECCA,
    "ALPHA CRB": SE_ALPHECCA,
    "5 CRB": SE_ALPHECCA,
    "α CRB": SE_ALPHECCA,
    "ALCRB": SE_ALPHECCA,
    "GNOSIA STELLA": SE_ALPHECCA,
    "MUNIR AL-FAKKAH": SE_ALPHECCA,
    "ASHTAROTH": SE_ALPHECCA,
    # ======== DENEB ALGEDI (Delta Capricorni) - BEHENIAN ========
    "TAIL OF THE GOAT": SE_DENEB_ALGEDI,
    "DELTA CAPRICORNI": SE_DENEB_ALGEDI,
    "DELTA CAP": SE_DENEB_ALGEDI,
    "49 CAP": SE_DENEB_ALGEDI,
    "δ CAP": SE_DENEB_ALGEDI,
    "DECAP": SE_DENEB_ALGEDI,
    "SCHEDDI": SE_DENEB_ALGEDI,
    "DHANAB AL-JADY": SE_DENEB_ALGEDI,
    # ======== ASTEROPE (21 Tauri - Pleiades) ========
    "21 TAURI": SE_ASTEROPE,
    "21 TAU": SE_ASTEROPE,
    "STEROPE": SE_ASTEROPE,
    "STEROPE I": SE_ASTEROPE,
    # ======== CELAENO (16 Tauri - Pleiades) ========
    "16 TAURI": SE_CELAENO,
    "16 TAU": SE_CELAENO,
    "CELENO": SE_CELAENO,
    # ======== ELECTRA (17 Tauri - Pleiades) ========
    "17 TAURI": SE_ELECTRA,
    "17 TAU": SE_ELECTRA,
    # ======== MAIA (20 Tauri - Pleiades) ========
    "20 TAURI": SE_MAIA,
    "20 TAU": SE_MAIA,
    # ======== MEROPE (23 Tauri - Pleiades) ========
    "23 TAURI": SE_MEROPE,
    "23 TAU": SE_MEROPE,
    # ======== TAYGETA (19 Tauri - Pleiades) ========
    "19 TAURI": SE_TAYGETA,
    "19 TAU": SE_TAYGETA,
    "TAYGETE": SE_TAYGETA,
    # ======== ATLAS (27 Tauri - Pleiades) ========
    "27 TAURI": SE_ATLAS,
    "27 TAU": SE_ATLAS,
    # ======== PLEIONE (28 Tauri - Pleiades) ========
    "28 TAURI": SE_PLEIONE,
    "28 TAU": SE_PLEIONE,
    # ======== PRIMA HYADUM (Gamma Tauri - Hyades) ========
    "GAMMA TAURI": SE_PRIMA_HYADUM,
    "GAMMA TAU": SE_PRIMA_HYADUM,
    "54 TAU": SE_PRIMA_HYADUM,
    "γ TAU": SE_PRIMA_HYADUM,
    "GATAU": SE_PRIMA_HYADUM,
    "HYADUM I": SE_PRIMA_HYADUM,
    "FIRST HYAD": SE_PRIMA_HYADUM,
    # ======== SECUNDA HYADUM (Delta^1 Tauri - Hyades) ========
    "DELTA TAURI": SE_SECUNDA_HYADUM,
    "DELTA TAU": SE_SECUNDA_HYADUM,
    "DELTA1 TAURI": SE_SECUNDA_HYADUM,
    "DELTA1 TAU": SE_SECUNDA_HYADUM,
    "61 TAU": SE_SECUNDA_HYADUM,
    "δ TAU": SE_SECUNDA_HYADUM,
    "DETAU": SE_SECUNDA_HYADUM,
    "HYADUM II": SE_SECUNDA_HYADUM,
    "SECOND HYAD": SE_SECUNDA_HYADUM,
    # ======== THETA TAURI (Theta^2 Tauri - Hyades) ========
    "THETA2 TAURI": SE_THETA_TAURI,
    "THETA2 TAU": SE_THETA_TAURI,
    "78 TAU": SE_THETA_TAURI,
    "θ TAU": SE_THETA_TAURI,
    "THTAU": SE_THETA_TAURI,
    "THETA^2 TAURI": SE_THETA_TAURI,
    "THETA^2 TAU": SE_THETA_TAURI,
    # ======== AIN (Epsilon Tauri - Hyades) ========
    "EPSILON TAURI": SE_AIN,
    "EPSILON TAU": SE_AIN,
    "74 TAU": SE_AIN,
    "ε TAU": SE_AIN,
    "EPTAU": SE_AIN,
    "OCULUS BOREALIS": SE_AIN,
    "BULL'S EYE": SE_AIN,
    # ======== SARGAS (Theta Scorpii) ========
    "THETA SCORPII": SE_SARGAS,
    "THETA SCO": SE_SARGAS,
    "θ SCO": SE_SARGAS,
    "THSCO": SE_SARGAS,
    "GIRTAB": SE_SARGAS,
    "SCORPION'S TAIL": SE_SARGAS,
    # ======== DSCHUBBA (Delta Scorpii) ========
    "DELTA SCORPII": SE_DSCHUBBA,
    "DELTA SCO": SE_DSCHUBBA,
    "δ SCO": SE_DSCHUBBA,
    "DESCO": SE_DSCHUBBA,
    "ICLARCRAU": SE_DSCHUBBA,
    "ICLARKRAV": SE_DSCHUBBA,
    "SCORPION'S FOREHEAD": SE_DSCHUBBA,
    # ======== GRAFFIAS (Beta Scorpii) ========
    "BETA SCORPII": SE_GRAFFIAS,
    "BETA SCO": SE_GRAFFIAS,
    "β SCO": SE_GRAFFIAS,
    "BESCO": SE_GRAFFIAS,
    "ACRAB": SE_GRAFFIAS,
    "AKRAB": SE_GRAFFIAS,
    "ELACRAB": SE_GRAFFIAS,
    # ======== LESATH (Upsilon Scorpii) ========
    "UPSILON SCORPII": SE_LESATH,
    "UPSILON SCO": SE_LESATH,
    "υ SCO": SE_LESATH,
    "UPSCO": SE_LESATH,
    "STINGER": SE_LESATH,
    "34 SCO": SE_LESATH,
    # ======== ZOSMA (Delta Leonis) ========
    "DELTA LEONIS": SE_ZOSMA,
    "DELTA LEO": SE_ZOSMA,
    "68 LEO": SE_ZOSMA,
    "δ LEO": SE_ZOSMA,
    "DELEO": SE_ZOSMA,
    "DHUR": SE_ZOSMA,
    "DUHR": SE_ZOSMA,
    "LION'S HIP": SE_ZOSMA,
    "LION'S BACK": SE_ZOSMA,
    # ======== ZODIACAL CONSTELLATION BRIGHT STARS ========
    # ======== HAMAL (Alpha Arietis) ========
    "ALPHA ARIETIS": SE_HAMAL,
    "ALPHA ARI": SE_HAMAL,
    "13 ARI": SE_HAMAL,
    "α ARI": SE_HAMAL,
    "ALARI": SE_HAMAL,
    "RAM'S HEAD": SE_HAMAL,
    # ======== SHERATAN (Beta Arietis) ========
    "BETA ARIETIS": SE_SHERATAN,
    "BETA ARI": SE_SHERATAN,
    "6 ARI": SE_SHERATAN,
    "β ARI": SE_SHERATAN,
    "BEARI": SE_SHERATAN,
    "SHARATAN": SE_SHERATAN,
    "AL-SHARATAIN": SE_SHERATAN,
    # ======== MESARTHIM (Gamma Arietis) ========
    "GAMMA ARIETIS": SE_MESARTHIM,
    "GAMMA ARI": SE_MESARTHIM,
    "5 ARI": SE_MESARTHIM,
    "γ ARI": SE_MESARTHIM,
    "GAARI": SE_MESARTHIM,
    "MESARTIM": SE_MESARTHIM,
    "FIRST STAR OF ARIES": SE_MESARTHIM,
    # ======== ACUBENS (Alpha Cancri) ========
    "ALPHA CANCRI": SE_ACUBENS,
    "ALPHA CNC": SE_ACUBENS,
    "65 CNC": SE_ACUBENS,
    "α CNC": SE_ACUBENS,
    "ALCNC": SE_ACUBENS,
    "SERTAN": SE_ACUBENS,
    "AL ZUBANAH": SE_ACUBENS,
    # ======== TARF (Beta Cancri) ========
    "BETA CANCRI": SE_TARF,
    "BETA CNC": SE_TARF,
    "17 CNC": SE_TARF,
    "β CNC": SE_TARF,
    "BECNC": SE_TARF,
    "AL TARF": SE_TARF,
    # ======== ASELLUS BOREALIS (Gamma Cancri) ========
    "GAMMA CANCRI": SE_ASELLUS_BOREALIS,
    "GAMMA CNC": SE_ASELLUS_BOREALIS,
    "43 CNC": SE_ASELLUS_BOREALIS,
    "γ CNC": SE_ASELLUS_BOREALIS,
    "GACNC": SE_ASELLUS_BOREALIS,
    "NORTHERN DONKEY": SE_ASELLUS_BOREALIS,
    "NORTHERN ASS": SE_ASELLUS_BOREALIS,
    # ======== ASELLUS AUSTRALIS (Delta Cancri) ========
    "DELTA CANCRI": SE_ASELLUS_AUSTRALIS,
    "DELTA CNC": SE_ASELLUS_AUSTRALIS,
    "47 CNC": SE_ASELLUS_AUSTRALIS,
    "δ CNC": SE_ASELLUS_AUSTRALIS,
    "DECNC": SE_ASELLUS_AUSTRALIS,
    "SOUTHERN DONKEY": SE_ASELLUS_AUSTRALIS,
    "SOUTHERN ASS": SE_ASELLUS_AUSTRALIS,
    # ======== KAUS AUSTRALIS (Epsilon Sagittarii) ========
    "EPSILON SAGITTARII": SE_KAUS_AUSTRALIS,
    "EPSILON SGR": SE_KAUS_AUSTRALIS,
    "20 SGR": SE_KAUS_AUSTRALIS,
    "ε SGR": SE_KAUS_AUSTRALIS,
    "EPSGR": SE_KAUS_AUSTRALIS,
    "SOUTHERN BOW": SE_KAUS_AUSTRALIS,
    # ======== NUNKI (Sigma Sagittarii) ========
    "SIGMA SAGITTARII": SE_NUNKI,
    "SIGMA SGR": SE_NUNKI,
    "34 SGR": SE_NUNKI,
    "σ SGR": SE_NUNKI,
    "SISGR": SE_NUNKI,
    "PELAGUS": SE_NUNKI,
    # ======== KAUS MEDIA (Delta Sagittarii) ========
    "DELTA SAGITTARII": SE_KAUS_MEDIA,
    "DELTA SGR": SE_KAUS_MEDIA,
    "19 SGR": SE_KAUS_MEDIA,
    "δ SGR": SE_KAUS_MEDIA,
    "DESGR": SE_KAUS_MEDIA,
    "MIDDLE BOW": SE_KAUS_MEDIA,
    # ======== KAUS BOREALIS (Lambda Sagittarii) ========
    "LAMBDA SAGITTARII": SE_KAUS_BOREALIS,
    "LAMBDA SGR": SE_KAUS_BOREALIS,
    "22 SGR": SE_KAUS_BOREALIS,
    "λ SGR": SE_KAUS_BOREALIS,
    "LASGR": SE_KAUS_BOREALIS,
    "NORTHERN BOW": SE_KAUS_BOREALIS,
    # ======== ASCELLA (Zeta Sagittarii) ========
    "ZETA SAGITTARII": SE_ASCELLA,
    "ZETA SGR": SE_ASCELLA,
    "38 SGR": SE_ASCELLA,
    "ζ SGR": SE_ASCELLA,
    "ZESGR": SE_ASCELLA,
    "ARMPIT": SE_ASCELLA,
    # ======== ALGEDI (Alpha Capricorni) ========
    "ALPHA CAPRICORNI": SE_ALGEDI,
    "ALPHA CAP": SE_ALGEDI,
    "6 CAP": SE_ALGEDI,
    "α CAP": SE_ALGEDI,
    "ALCAP": SE_ALGEDI,
    "GIEDI": SE_ALGEDI,
    "PRIMA GIEDI": SE_ALGEDI,
    "THE GOAT": SE_ALGEDI,
    # ======== DABIH (Beta Capricorni) ========
    "BETA CAPRICORNI": SE_DABIH,
    "BETA CAP": SE_DABIH,
    "9 CAP": SE_DABIH,
    "β CAP": SE_DABIH,
    "BECAP": SE_DABIH,
    "AL-DHABIH": SE_DABIH,
    "LUCKY ONE OF SLAUGHTERER": SE_DABIH,
    # ======== NASHIRA (Gamma Capricorni) ========
    "GAMMA CAPRICORNI": SE_NASHIRA,
    "GAMMA CAP": SE_NASHIRA,
    "40 CAP": SE_NASHIRA,
    "γ CAP": SE_NASHIRA,
    "GACAP": SE_NASHIRA,
    "FORTUNATE ONE": SE_NASHIRA,
    "SA'D NASHIRAH": SE_NASHIRA,
    # ======== SADALSUUD (Beta Aquarii) ========
    "BETA AQUARII": SE_SADALSUUD,
    "BETA AQR": SE_SADALSUUD,
    "22 AQR": SE_SADALSUUD,
    "β AQR": SE_SADALSUUD,
    "BEAQR": SE_SADALSUUD,
    "LUCKIEST OF LUCKY STARS": SE_SADALSUUD,
    "SA'D AL-SU'UD": SE_SADALSUUD,
    # ======== SADALMELIK (Alpha Aquarii) ========
    "ALPHA AQUARII": SE_SADALMELIK,
    "ALPHA AQR": SE_SADALMELIK,
    "34 AQR": SE_SADALMELIK,
    "α AQR": SE_SADALMELIK,
    "ALAQR": SE_SADALMELIK,
    "LUCKY STAR OF KING": SE_SADALMELIK,
    "SA'D AL-MALIK": SE_SADALMELIK,
    # ======== SKAT (Delta Aquarii) ========
    "DELTA AQUARII": SE_SKAT,
    "DELTA AQR": SE_SKAT,
    "76 AQR": SE_SKAT,
    "δ AQR": SE_SKAT,
    "DEAQR": SE_SKAT,
    "SCHEAT AQUARII": SE_SKAT,
    "SHIN": SE_SKAT,
    # ======== ETA PISCIUM ========
    "ETA PISCIUM": SE_ETA_PISCIUM,
    "ETA PSC": SE_ETA_PISCIUM,
    "99 PSC": SE_ETA_PISCIUM,
    "η PSC": SE_ETA_PISCIUM,
    "ETPSC": SE_ETA_PISCIUM,
    "KULLAT NUNU": SE_ETA_PISCIUM,
    # ======== ALRESCHA (Alpha Piscium) ========
    "ALPHA PISCIUM": SE_ALRESCHA,
    "ALPHA PSC": SE_ALRESCHA,
    "113 PSC": SE_ALRESCHA,
    "α PSC": SE_ALRESCHA,
    "ALPSC": SE_ALRESCHA,
    "AL-RISHA": SE_ALRESCHA,
    "THE KNOT": SE_ALRESCHA,
    "THE CORD": SE_ALRESCHA,
    # ======== ALTERNATE SPELLINGS / COMMON MISSPELLINGS ========
    # Betelgeuse variants (Arabic: yad al-jawza / ibt al-jawza)
    "BETELGEUX": SE_BETELGEUSE,
    "BEETLEJUICE": SE_BETELGEUSE,
    "BETELGUESE": SE_BETELGEUSE,
    "BETELGUEUSE": SE_BETELGEUSE,
    "BETELEGEUSE": SE_BETELGEUSE,
    "BETEIGEUZE": SE_BETELGEUSE,
    # Fomalhaut variants (Arabic: fum al-hut)
    "FORMALHAUT": SE_FOMALHAUT,
    "FOMALHUT": SE_FOMALHAUT,
    "FOMALAUT": SE_FOMALHAUT,
    "FOMALHAULT": SE_FOMALHAUT,
    "FUMALHAUT": SE_FOMALHAUT,
    # Aldebaran variants (Arabic: al-dabaran)
    "ALDEBRAN": SE_ALDEBARAN,
    "ALDEBERON": SE_ALDEBARAN,
    "ALDEBERAN": SE_ALDEBARAN,
    "ALDERBARAN": SE_ALDEBARAN,
    "ALDEBARIN": SE_ALDEBARAN,
    # Algol variants (Arabic: ra's al-ghul)
    "ALGOL": SE_ALGOL,  # Already canonical but ensure alias exists
    "ALGHOL": SE_ALGOL,
    "ALGOUL": SE_ALGOL,
    "AL-GHUL": SE_ALGOL,
    # Arcturus variants (Greek: arktos + ouros = bear guard)
    "ARCHTURUS": SE_ARCTURUS,
    "ARTURUS": SE_ARCTURUS,
    "ARKCTURUS": SE_ARCTURUS,
    "ARCTUROS": SE_ARCTURUS,
    # Antares variants (Greek: anti + Ares = rival of Mars)
    "ANTARIES": SE_ANTARES,
    "ANTARAS": SE_ANTARES,
    "ANTARRES": SE_ANTARES,
    # Rigel variants (Arabic: rijl = foot)
    "RIEGEL": SE_RIGEL,
    "RIJEL": SE_RIGEL,
    "RIGIL": SE_RIGEL,
    # Vega variants (Arabic: an-nasr al-waqi)
    "VEGHA": SE_VEGA,
    # Polaris variants
    "POLARRIS": SE_POLARIS,
    "POLARIUS": SE_POLARIS,
    # Procyon variants (Greek: pro + kyon = before the dog)
    "PROCION": SE_PROCYON,
    "PROCIAN": SE_PROCYON,
    # Capella variants (Latin: she-goat)
    "CAPELA": SE_CAPELLA,
    "CAPPELLA": SE_CAPELLA,
    # Deneb variants (Arabic: dhanab = tail)
    "DANEB": SE_DENEB,
    "DHENEB": SE_DENEB,
    # Altair variants (Arabic: al-nasr al-ta'ir = flying eagle)
    "ALTIAR": SE_ALTAIR,
    "ALTARE": SE_ALTAIR,
    # Sirius variants (Greek: seirios = scorching)
    "SYRIUS": SE_SIRIUS,
    "SIRUIS": SE_SIRIUS,
    "SIRUS": SE_SIRIUS,
    # Spica variants (Latin: ear of wheat)
    "SPIKA": SE_SPICA_STAR,
    "SPICCA": SE_SPICA_STAR,
    # Regulus variants (Latin: little king)
    "REGULAS": SE_REGULUS,
    "REGULIS": SE_REGULUS,
    # Canopus variants
    "CANOPIS": SE_CANOPUS,
    "CANPOUS": SE_CANOPUS,
    "CANOPOS": SE_CANOPUS,
    # Achernar variants (Arabic: akhir an-nahr = end of river)
    "ACHENAR": SE_ACHERNAR,
    "ARCHENAR": SE_ACHERNAR,
    "ACHERNAN": SE_ACHERNAR,
    # Castor/Pollux variants
    "CASTOR": SE_CASTOR,  # Ensure canonical exists
    "KASTOR": SE_CASTOR,
    "POLUX": SE_POLLUX,
    "POLLUCKS": SE_POLLUX,
    # Alpheratz variants (Alpha Andromedae)
    "ALPHERATZ": SE_ALPHERATZ,
    "SIRRAH": SE_ALPHERATZ,
    "SIRAH": SE_ALPHERATZ,
    "ALPHERAT": SE_ALPHERATZ,
    "ALPHA ANDROMEDAE": SE_ALPHERATZ,
    # Algenib variants (Gamma Pegasi)
    "ALGENIB": SE_ALGENIB,
    "GAMMA PEGASI": SE_ALGENIB,
    # Propus variants (Eta Geminorum)
    "PROPUS": SE_PROPUS,
    "ETA GEMINORUM": SE_PROPUS,
    # Tejat variants (Mu Geminorum)
    "TEJAT": SE_TEJAT,
    "TEJAT POSTERIOR": SE_TEJAT,
    "MU GEMINORUM": SE_TEJAT,
    # Alhena variants (Gamma Geminorum)
    "ALHENA": SE_ALHENA,
    "ALMEISAN": SE_ALHENA,
    "GAMMA GEMINORUM": SE_ALHENA,
    # Wasat variants (Delta Geminorum)
    "WASAT": SE_WASAT,
    "DELTA GEMINORUM": SE_WASAT,
    # Adhara variants (Epsilon Canis Majoris)
    "ADHARA": SE_ADHARA,
    "ADARA": SE_ADHARA,
    "EPSILON CANIS MAJORIS": SE_ADHARA,
    # Wezen variants (Delta Canis Majoris)
    "WEZEN": SE_WEZEN,
    "WESEN": SE_WEZEN,
    "DELTA CANIS MAJORIS": SE_WEZEN,
    # Thuban variants (Alpha Draconis)
    "THUBAN": SE_THUBAN,
    "ALPHA DRACONIS": SE_THUBAN,
    # Rasalgethi variants (Alpha Herculis)
    "RASALGETHI": SE_RASALGETHI,
    "RAS ALGETHI": SE_RASALGETHI,
    "ALPHA HERCULIS": SE_RASALGETHI,
    # Albireo variants (Beta Cygni)
    "ALBIREO": SE_ALBIREO,
    "BETA CYGNI": SE_ALBIREO,
    # Mirach variants (Beta Andromedae)
    "MIRACH": SE_MIRACH,
    "BETA ANDROMEDAE": SE_MIRACH,
    "BETA AND": SE_MIRACH,
    "β AND": SE_MIRACH,
    "BEAND": SE_MIRACH,
    "MIRAK": SE_MIRACH,
    "MIRAC": SE_MIRACH,
    # Almach variants (Gamma Andromedae)
    "ALMACH": SE_ALMACH,
    "ALAMAK": SE_ALMACH,
    "ALMAAK": SE_ALMACH,
    "GAMMA ANDROMEDAE": SE_ALMACH,
    "GAMMA AND": SE_ALMACH,
    "γ AND": SE_ALMACH,
    "GA1AND": SE_ALMACH,
    # Menkar variants (Alpha Ceti)
    "MENKAR": SE_MENKAR,
    "MENKAB": SE_MENKAR,
    "ALPHA CETI": SE_MENKAR,
    "ALPHA CET": SE_MENKAR,
    "α CET": SE_MENKAR,
    "ALCET": SE_MENKAR,
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
        SE_BETELGEUSE
        >>> _fuzzy_match_star("Formalhaut")  # Common misspelling
        SE_FOMALHAUT
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
        SE_BETELGEUSE
        >>> resolve_star_name("Formalhaut")  # Common misspelling
        SE_FOMALHAUT
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

    # 1. Try exact match in STAR_ALIASES
    if normalized in STAR_ALIASES:
        return STAR_ALIASES[normalized]

    # 2. Try exact match against canonical star names
    if normalized in _STAR_NAME_TO_ID:
        return _STAR_NAME_TO_ID[normalized]

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
        parallax_mas=star_data.parallax_mas,  # Trigonometric parallax for distance
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
    # Distance from parallax (AU). Stars without parallax data (parallax_mas=0)
    # will have a very large distance from Skyfield; cap at 100000 AU for those.
    dist = ecl[2].au
    if star_data.parallax_mas == 0.0:
        dist = 100000.0

    return lon, lat, dist


def _calc_star_position_skyfield(
    star_id: int,
    jd_tt: float,
    noaberr: bool = False,
    nogdefl: bool = False,
    j2000_frame: bool = False,
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
) -> Tuple[float, float, float]:
    """
    Calculate ecliptic position of a fixed star at given date.

    Uses Skyfield Star class for accurate position calculation including:
    - Proper motion
    - Precession and nutation
    - Annual aberration (by default)
    - Gravitational light deflection (by default)

    Args:
        star_id: Star identifier (SE_REGULUS, SE_SPICA_STAR, etc.)
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
    return _calc_star_position_skyfield(star_id, jd_tt, noaberr, nogdefl, j2000_frame)


def calc_fixed_star_velocity(
    star_id: int,
    jd_tt: float,
    noaberr: bool = False,
    nogdefl: bool = False,
    j2000_frame: bool = False,
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
        star_id: Star identifier (SE_REGULUS, SE_SPICA_STAR, etc.)
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
        star_id, jd_tt, noaberr, nogdefl, j2000_frame
    )

    # Calculate positions at t-0.5 and t+0.5 for central difference
    lon_prev, lat_prev, dist_prev = calc_fixed_star_position(
        star_id, jd_tt - h, noaberr, nogdefl, j2000_frame
    )
    lon_next, lat_next, dist_next = calc_fixed_star_position(
        star_id, jd_tt + h, noaberr, nogdefl, j2000_frame
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

    # pyswisseph returns 0.0 for fixed star distance speed
    speed_dist = 0.0

    return lon, lat, dist, speed_lon, speed_lat, speed_dist


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
    star_id = resolve_star_name(star_name)
    if star_id is not None:
        canonical_name = get_canonical_star_name(star_id)
        return star_id, None, canonical_name

    return -1, f"could not find star name {star_name.lower()}", None


def _preprocess_flags(iflag: int) -> int:
    """Preprocess calculation flags for fixed star functions.

    Strips ephemeris selection flags (MOSEPH) and converts SPEED3 to SPEED.
    SEFLG_TOPOCTR is accepted silently (stars at infinite distance).

    Args:
        iflag: Raw input flags

    Returns:
        Cleaned flags
    """
    # Strip SEFLG_MOSEPH — accepted for compatibility, always uses Skyfield
    iflag = iflag & ~SEFLG_MOSEPH
    # SEFLG_SPEED3: treat as SEFLG_SPEED
    if iflag & SEFLG_SPEED3:
        iflag = (iflag & ~SEFLG_SPEED3) | SEFLG_SPEED
    # SEFLG_TOPOCTR: silently accept (stars at infinite distance, no parallax)
    iflag = iflag & ~SEFLG_TOPOCTR
    return iflag


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

    is_equatorial = bool(iflag & SEFLG_EQUATORIAL)

    # ---- 1. Frame selection ----
    # J2000: precess from ecliptic-of-date back to J2000 ecliptic
    # NONUT: remove nutation from ecliptic longitude (mean ecliptic of date)
    # These apply regardless of equatorial output (adjust ecliptic first).
    if iflag & SEFLG_J2000:
        if not j2000_native:
            from .astrometry import _precess_ecliptic

            lon, lat = _precess_ecliptic(lon, lat, jd_tt, J2000)
        # else: already in J2000 frame from Skyfield's ecliptic_J2000_frame
    elif iflag & SEFLG_NONUT:
        # Skyfield returns positions on the true ecliptic of date (with nutation).
        # NONUT means output on the mean ecliptic of date, so subtract dpsi
        # (nutation in longitude) from the ecliptic longitude.
        from .cache import get_cached_nutation

        dpsi_rad, _ = get_cached_nutation(jd_tt)
        lon = (lon - math.degrees(dpsi_rad)) % 360.0

    # ---- 2. Equatorial coordinate transformation ----
    if is_equatorial:
        if iflag & SEFLG_J2000:
            # J2000 obliquity for J2000 equatorial frame
            eps = 23.4392911  # IAU 2006 mean obliquity at J2000.0
        elif iflag & SEFLG_NONUT:
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
    if iflag & SEFLG_SIDEREAL:
        from .state import get_timescale

        ts = get_timescale()
        t = ts.tt_jd(jd_tt)
        tjd_ut = t.ut1

        from .planets import swe_get_ayanamsa_ut

        ayanamsa = swe_get_ayanamsa_ut(tjd_ut)
        lon = (lon - ayanamsa) % 360.0

    # ---- 4. Output format conversion ----
    result = (lon, lat, dist, speed_lon, speed_lat, speed_dist)

    if iflag & SEFLG_XYZ:
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

    if iflag & SEFLG_RADIANS:
        return (
            math.radians(lon),
            math.radians(lat),
            dist,
            math.radians(speed_lon),
            math.radians(speed_lat),
            speed_dist,
        )

    return result


def swe_fixstar_ut(
    star: str, tjdut: float, flags: int = SEFLG_SWIEPH
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
        >>> pos, name, retflag = swe_fixstar_ut("Regulus", 2451545.0, 0)
        >>> lon, lat, dist = pos[0], pos[1], pos[2]
    """
    flags = _preprocess_flags(flags)

    star_id, error, canonical_name = _resolve_star_id(star)
    if error:
        raise Error(error)

    # Convert UT to TT using timescale (applies Delta T)
    from .state import get_timescale

    t = get_timescale().ut1_jd(tjdut)

    try:
        noaberr = bool(flags & SEFLG_NOABERR) or bool(flags & SEFLG_TRUEPOS)
        nogdefl = bool(flags & SEFLG_NOGDEFL)
        # Compute natively in J2000 ecliptic frame when requested.
        # This avoids the ~5" error from precessing Skyfield's ecliptic-of-date
        # back to J2000 with a different precession model.
        use_j2000 = bool(flags & SEFLG_J2000)

        if flags & SEFLG_SPEED:
            lon, lat, dist, speed_lon, speed_lat, speed_dist = calc_fixed_star_velocity(
                star_id, t.tt, noaberr, nogdefl, j2000_frame=use_j2000
            )
            result = (lon, lat, dist, speed_lon, speed_lat, speed_dist)
        else:
            lon, lat, dist = calc_fixed_star_position(
                star_id, t.tt, noaberr, nogdefl, j2000_frame=use_j2000
            )
            result = (lon, lat, dist, 0.0, 0.0, 0.0)

        result = _apply_fixstar_flags(result, t.tt, flags, j2000_native=use_j2000)

        return (result, canonical_name or "", flags)
    except Error:
        raise
    except (OSError, ValueError, KeyError) as e:
        raise Error(str(e)) from e


def swe_batch_fixstars_ut(
    stars: Sequence[str],
    tjdut: float,
    flags: int = SEFLG_SWIEPH,
    *,
    skip_errors: bool = False,
) -> Tuple[
    Tuple[Tuple[float, float, float, float, float, float], str, int] | None, ...
]:
    """Calculate multiple fixed-star positions for Universal Time.

    The result order matches the input order. When ``skip_errors`` is True,
    unresolved stars keep their input slot as ``None``.
    """
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

    from .state import get_planets, get_timescale

    noaberr = bool(flags & SEFLG_NOABERR) or bool(flags & SEFLG_TRUEPOS)
    nogdefl = bool(flags & SEFLG_NOGDEFL)
    use_j2000 = bool(flags & SEFLG_J2000)
    want_speed = bool(flags & SEFLG_SPEED)

    ts = get_timescale()
    t = ts.ut1_jd(tjdut)
    earth = get_planets()["earth"]
    earth_at_t = earth.at(t)

    if want_speed:
        t_prev = ts.tt_jd(t.tt - 0.5)
        t_next = ts.tt_jd(t.tt + 0.5)
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
                lon_prev, lat_prev, _dist_prev = _calc_star_position_from_observer(
                    star_id, earth_at_prev, noaberr, nogdefl, use_j2000
                )
                lon_next, lat_next, _dist_next = _calc_star_position_from_observer(
                    star_id, earth_at_next, noaberr, nogdefl, use_j2000
                )
                speed_lon = lon_next - lon_prev
                if speed_lon > 180.0:
                    speed_lon -= 360.0
                elif speed_lon < -180.0:
                    speed_lon += 360.0
                speed_lat = lat_next - lat_prev
                speed_dist = 0.0
            else:
                speed_lon = 0.0
                speed_lat = 0.0
                speed_dist = 0.0

            result = (lon, lat, dist, speed_lon, speed_lat, speed_dist)
            result = _apply_fixstar_flags(result, t.tt, flags, j2000_native=use_j2000)
            results[index] = (result, canonical_name, flags)
        except Error:
            if skip_errors:
                continue
            raise
        except (OSError, ValueError, KeyError) as e:
            if skip_errors:
                continue
            raise Error(str(e)) from e

    return tuple(results)


batch_fixstars_ut = swe_batch_fixstars_ut


def swe_fixstar(
    star: str, tjdet: float, flags: int = SEFLG_SWIEPH
) -> Tuple[Tuple[float, float, float, float, float, float], str, int]:
    """
    Calculate position of a fixed star for Terrestrial Time (TT).

    Reference API compatible function. Similar to swe_fixstar_ut() but takes
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
        For most astrological applications, use swe_fixstar_ut() instead.

    Example:
        >>> pos, name, retflag = swe_fixstar("Regulus", 2451545.0, 0)
        >>> lon, lat, dist = pos[0], pos[1], pos[2]
    """
    flags = _preprocess_flags(flags)

    star_id, error, canonical_name = _resolve_star_id(star)
    if error:
        raise Error(error)

    try:
        noaberr = bool(flags & SEFLG_NOABERR) or bool(flags & SEFLG_TRUEPOS)
        nogdefl = bool(flags & SEFLG_NOGDEFL)
        use_j2000 = bool(flags & SEFLG_J2000)

        if flags & SEFLG_SPEED:
            lon, lat, dist, speed_lon, speed_lat, speed_dist = calc_fixed_star_velocity(
                star_id, tjdet, noaberr, nogdefl, j2000_frame=use_j2000
            )
            result = (lon, lat, dist, speed_lon, speed_lat, speed_dist)
        else:
            lon, lat, dist = calc_fixed_star_position(
                star_id, tjdet, noaberr, nogdefl, j2000_frame=use_j2000
            )
            result = (lon, lat, dist, 0.0, 0.0, 0.0)

        result = _apply_fixstar_flags(result, tjdet, flags, j2000_native=use_j2000)

        return (result, canonical_name or "", flags)
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


def swe_fixstar2_ut(
    star: str, tjdut: float, flags: int = SEFLG_SWIEPH
) -> Tuple[Tuple[float, float, float, float, float, float], str, int]:
    """
    Calculate position of a fixed star for Universal Time with flexible lookup.

    Enhanced version of swe_fixstar_ut() that supports flexible star lookup:
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
        >>> pos, name, retflag = swe_fixstar2_ut("Reg", 2451545.0, 0)
        >>> print(name)  # "Regulus,alLeo"
        >>> lon, lat, dist = pos[0], pos[1], pos[2]

        >>> pos, name, retflag = swe_fixstar2_ut("49669", 2451545.0, 0)
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
        noaberr = bool(flags & SEFLG_NOABERR) or bool(flags & SEFLG_TRUEPOS)
        nogdefl = bool(flags & SEFLG_NOGDEFL)
        use_j2000 = bool(flags & SEFLG_J2000)

        if flags & SEFLG_SPEED:
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


def swe_fixstar2(
    star: str, tjdet: float, flags: int = SEFLG_SWIEPH
) -> Tuple[Tuple[float, float, float, float, float, float], str, int]:
    """
    Calculate position of a fixed star for Terrestrial Time with flexible lookup.

    Enhanced version of swe_fixstar() that supports flexible star lookup:
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
        For most astrological applications, use swe_fixstar2_ut() instead.

    Example:
        >>> pos, name, retflag = swe_fixstar2("Spica", 2451545.0, 0)
        >>> print(name)  # "Spica,alVir"
        >>> lon, lat, dist = pos[0], pos[1], pos[2]

        >>> pos, name, retflag = swe_fixstar2("65474", 2451545.0, 0)
        >>> print(name)  # "Spica,alVir" (looked up by HIP number)
    """
    flags = _preprocess_flags(flags)

    entry, error = _resolve_star2(star)
    if error or entry is None:
        raise Error(error or "could not find star name")

    try:
        noaberr = bool(flags & SEFLG_NOABERR) or bool(flags & SEFLG_TRUEPOS)
        nogdefl = bool(flags & SEFLG_NOGDEFL)
        use_j2000 = bool(flags & SEFLG_J2000)

        if flags & SEFLG_SPEED:
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
    "ALRESCHA": 7097,  # Alpha Piscium (IAU: HIP 9487 for A component, HIP 7097 is the catalog entry)
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
    "SALM": 98066,  # Tau Pegasi
    "SARGAS": 86228,  # Theta Scorpii
    "SARIN": 79992,  # Delta Herculis
    "SCHEAT": 113881,  # Beta Pegasi
    "SCHEDAR": 3179,  # Alpha Cassiopeiae
    "SECUNDA HYADUM": 20455,  # Delta1 Tauri
    "SEGIN": 4427,  # Epsilon Cassiopeiae
    "SEGINUS": 71075,  # Gamma Bootis
    "SHAULA": 85927,  # Lambda Scorpii
    "SHAMA": 69701,  # HD 99109
    "SHERATAN": 8903,  # Beta Arietis
    "SIKA": 50782,  # HD 99491
    "SIRIUS": 32349,  # Alpha Canis Majoris
    "SITULA": 110672,  # Kappa Aquarii
    "SKAT": 113136,  # Delta Aquarii
    "SPICA": 65474,  # Alpha Virginis
    "STRIBOR": 91085,  # HD 171028
    "SUBRA": 47508,  # Omicron Leonis
    "SUHAIL": 44816,  # Lambda Velorum
    "SULAFAT": 93194,  # Gamma Lyrae
    "SYRMA": 71957,  # Iota Virginis
    # T
    "TABIT": 22449,  # Pi3 Orionis
    "TAIYANGSHOU": 54539,  # Chi Ursae Majoris
    "TAIYI": 53759,  # 8 Draconis
    "TALITHA": 44127,  # Iota Ursae Majoris
    "TANIA AUSTRALIS": 51250,  # Mu Ursae Majoris
    "TANIA BOREALIS": 50801,  # Lambda Ursae Majoris
    "TARAZED": 95501,  # Gamma Aquilae
    "TAYGETA": 17531,  # 19 Tauri (Pleiades)
    "TEBERDA": 94256,  # HD 178813
    "TEGMINE": 43103,  # Zeta1 Cancri
    "TEJAT": 30343,  # Mu Geminorum (HIP 30343, corrected from erroneous 32362)
    "THUBAN": 68756,  # Alpha Draconis
    "TIAKI": 23015,  # Beta Gruis
    "TIANGUAN": 25930,  # Zeta Tauri
    "TIANYI": 52403,  # 7 Draconis
    "TITAWIN": 9683,  # Upsilon Andromedae
    "TOLIMAN": 71681,  # Alpha Centauri B
    "TONATIUH": 43177,  # HD 104985
    "TORCULAR": 6193,  # Omicron Piscium
    "TUREIS": 42913,  # Rho Puppis
    "TYL": 91919,  # Epsilon Draconis
    # U
    "UKDAH": 52863,  # Iota Hydrae
    "UNUKALHAI": 77070,  # Alpha Serpentis
    "UNURGUNITE": 34444,  # Sigma Canis Majoris
    "URUK": 116076,  # HD 231701
    # V
    "VEGA": 91262,  # Alpha Lyrae
    "VERITATE": 74793,  # 14 Andromedae
    "VINDEMIATRIX": 63608,  # Epsilon Virginis
    "WASAT": 35550,  # Delta Geminorum
    "WAZN": 27628,  # Beta Columbae
    "WEZEN": 34444,  # Delta Canis Majoris
    # X
    "XAMIDIMURA": 82514,  # Mu1 Scorpii
    # Y
    "YED POSTERIOR": 86284,  # Epsilon Ophiuchi
    "YED PRIOR": 83000,  # Delta Ophiuchi
    "YILDUN": 85822,  # Delta Ursae Minoris
    # Z
    "ZANIAH": 60129,  # Eta Virginis
    "ZAURAK": 18543,  # Gamma Eridani
    "ZAVIJAVA": 57757,  # Beta Virginis
    "ZHANG": 49641,  # Upsilon1 Hydrae
    "ZIBAL": 20535,  # Zeta Eridani
    "ZOSMA": 54872,  # Delta Leonis
    "ZUBENELGENUBI": 72622,  # Alpha2 Librae
    "ZUBENELHAKRABI": 76470,  # Gamma Librae
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
    "ETA PISCIUM": 5742,  # (consistent with STAR_CATALOG)
    "ETA PSC": 5742,
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


def swe_fixstar_mag(star: str) -> Tuple[float, str]:
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
        >>> mag, name = swe_fixstar_mag("Regulus")
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


def swe_fixstar2_mag(star: str) -> Tuple[float, str]:
    """
    Get the visual magnitude of a fixed star with flexible lookup.

    Enhanced version that supports flexible star lookup like swe_fixstar2:
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
        >>> mag, name = swe_fixstar2_mag("Reg")
        >>> print(f"{name}: {mag}")  # "Regulus,alLeo: 1.40"

        >>> mag, name = swe_fixstar2_mag("49669")
        >>> print(f"{name}: {mag}")  # "Regulus,alLeo: 1.40"
    """
    entry, error = _resolve_star2(star)
    if error or entry is None:
        raise Error(error or "could not find star name")

    star_name_out = _format_star_name(entry)
    return (entry.magnitude, star_name_out)
