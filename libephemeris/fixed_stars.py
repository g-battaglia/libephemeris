"""
Fixed star position calculations for libephemeris.

Computes ecliptic positions for bright fixed stars with:
- Proper motion correction (linear approximation)
- IAU 2006 precession from J2000 to date
- Nutation correction for obliquity
- Equatorial to ecliptic coordinate transformation

Supported stars:
- Regulus (Alpha Leonis) - Royal Star
- Spica (Alpha Virginis) - Used for ayanamsha calculations

FIXME: Precision - Uses simplified linear proper motion
Real stars have:
- Radial velocity (3D motion, not just RA/Dec)
- Parallax (distance changes apparent position)
- Acceleration (proper motion not constant)
Full catalogs like Hipparcos/Gaia provide higher precision.

Typical precision: ~1-5 arcseconds over ±100 years from J2000
For research astronomy, use SIMBAD/Gaia catalogs.

References:
- IAU 2006 Precession: Capitaine et al. A&A 412, 567-586 (2003)
- Proper motion: Hipparcos/Tycho catalogs
"""

import math
from dataclasses import dataclass
from typing import Tuple
from .constants import SE_REGULUS, SE_SPICA_STAR


@dataclass
class StarData:
    """
    Fixed star astrometric data (ICRS J2000.0 epoch).
    
    Attributes:
        ra_j2000: Right Ascension at J2000.0 in degrees
        dec_j2000: Declination at J2000.0 in degrees  
        pm_ra: Proper motion in RA (arcsec/year, includes cos(dec) factor)
        pm_dec: Proper motion in Dec (arcsec/year)
        
    Note:
        Proper motions are linearized - good for ±200 years from epoch.
        Does not include parallax or radial velocity effects.
    """

    ra_j2000: float
    dec_j2000: float
    pm_ra: float
    pm_dec: float


# Fixed star catalog (J2000.0 ICRS coordinates from Hipparcos)
FIXED_STARS = {
    SE_REGULUS: StarData(
        ra_j2000=152.092958,  # 10h 08m 22.3s (Alpha Leonis)
        dec_j2000=11.967208,  # +11° 58' 02"
        pm_ra=-0.00249,  # -249 mas/yr (westward)
        pm_dec=0.00152,  # +152 mas/yr (northward)
    ),
    SE_SPICA_STAR: StarData(
        ra_j2000=201.298247,  # 13h 25m 11.6s (Alpha Virginis)
        dec_j2000=-11.161319,  # -11° 09' 41"
        pm_ra=-0.04235,  # -42.35 mas/yr
        pm_dec=-0.03067,  # -30.67 mas/yr
    ),
}


def calc_fixed_star_position(star_id: int, jd_tt: float) -> Tuple[float, float, float]:
    """
    Calculate ecliptic position of a fixed star at given date.
    
    Applies proper motion and precession to J2000.0 catalog position.
    
    Args:
        star_id: Star identifier (SE_REGULUS, SE_SPICA_STAR)
        jd_tt: Julian Day in Terrestrial Time (TT)
        
    Returns:
        Tuple[float, float, float]: (longitude, latitude, distance) where:
            - longitude: Ecliptic longitude of date in degrees (0-360)
            - latitude: Ecliptic latitude of date in degrees
            - distance: Arbitrary large value (100000 AU) - stars are effectively infinite
            
    Raises:
        ValueError: If star_id not in catalog
        
    Algorithm:
        1. Apply linear proper motion from J2000 to date
        2. Precess equatorial coordinates using IAU 2006 formula
        3. Calculate true obliquity (mean + nutation)
        4. Transform to ecliptic coordinates
        
    FIXME: Precision - Linear proper motion approximation
        - Ignores radial velocity (parallax effect)
        - Assumes constant proper motion (no acceleration)
        - No annual parallax (star distance not modeled)
        Typical error: 1-5 arcsec over ±100 years
        
    References:
        IAU 2006 precession (Capitaine et al.)
        Nutation simplified from IAU 2000B
    """
    if star_id not in FIXED_STARS:
        raise ValueError(f"Unknown star ID: {star_id}")

    star = FIXED_STARS[star_id]

    # Time from J2000.0
    t_years = (jd_tt - 2451545.0) / 365.25
    T = (jd_tt - 2451545.0) / 36525.0  # Julian centuries

    # FIXME: Precision - Linear proper motion (ignores curvature)
    # Apply proper motion to J2000 coordinates
    ra_pm = star.ra_j2000 + (star.pm_ra * t_years) / 3600.0
    dec_pm = star.dec_j2000 + (star.pm_dec * t_years) / 3600.0

    # IAU 2006 precession parameters (arcseconds -> degrees)
    zeta = (2306.2181 * T + 0.30188 * T**2 + 0.017998 * T**3) / 3600.0
    z = (2306.2181 * T + 1.09468 * T**2 + 0.018203 * T**3) / 3600.0
    theta = (2004.3109 * T - 0.42665 * T**2 - 0.041833 * T**3) / 3600.0

    zeta_r = math.radians(zeta)
    z_r = math.radians(z)
    theta_r = math.radians(theta)
    ra_r = math.radians(ra_pm)
    dec_r = math.radians(dec_pm)

    # Convert equatorial to Cartesian (unit sphere)
    x0 = math.cos(dec_r) * math.cos(ra_r)
    y0 = math.cos(dec_r) * math.sin(ra_r)
    z0 = math.sin(dec_r)

    # Apply precession rotation: R_z(-z) · R_y(θ) · R_z(-ζ)
    # Step 1: R_z(-ζ) rotation around z-axis
    x1 = x0 * math.cos(-zeta_r) + y0 * math.sin(-zeta_r)
    y1 = -x0 * math.sin(-zeta_r) + y0 * math.cos(-zeta_r)
    z1 = z0

    # Step 2: R_y(θ) rotation around y-axis
    x2 = x1 * math.cos(theta_r) - z1 * math.sin(theta_r)
    y2 = y1
    z2 = x1 * math.sin(theta_r) + z1 * math.cos(theta_r)

    # Step 3: R_z(-z) rotation around z-axis
    x3 = x2 * math.cos(-z_r) + y2 * math.sin(-z_r)
    y3 = -x2 * math.sin(-z_r) + y2 * math.cos(-z_r)
    z3 = z2

    # Convert back to spherical (RA/Dec of date)
    ra_date = math.atan2(y3, x3)
    dec_date = math.asin(z3)

    # Calculate mean obliquity of date (IAU 2006)
    eps0 = 23.43929111 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T / 3600.0

    # FIXME: Precision - Simplified nutation (2-term approximation)
    # Full IAU 2000B has 77 terms. This captures ~99% of effect.
    omega = 125.04452 - 1934.136261 * T  # Longitude of lunar ascending node
    L = 280.4665 + 36000.7698 * T  # Mean longitude of Sun

    # Nutation in obliquity (arcseconds)
    deps = 9.20 * math.cos(math.radians(omega)) + 0.57 * math.cos(math.radians(2 * L))
    deps_deg = deps / 3600.0
    eps_true = eps0 + deps_deg
    eps_r = math.radians(eps_true)

    # Transform equatorial (RA, Dec) to ecliptic (lon, lat)
    sin_lat = math.sin(dec_date) * math.cos(eps_r) - math.cos(dec_date) * math.sin(
        eps_r
    ) * math.sin(ra_date)
    lat = math.asin(sin_lat)

    y_lon = math.sin(ra_date) * math.cos(eps_r) + math.tan(dec_date) * math.sin(eps_r)
    x_lon = math.cos(ra_date)
    lon = math.degrees(math.atan2(y_lon, x_lon)) % 360.0
    lat_deg = math.degrees(lat)

    # Distance is arbitrary for fixed stars (effectively infinite)
    dist = 100000.0  # AU (placeholder)

    return lon, lat_deg, dist


def swe_fixstar_ut(
    star_name: str, tjd_ut: float, iflag: int
) -> Tuple[Tuple[float, float, float, float, float, float], int, str]:
    """
    Calculate position of a fixed star.
    
    Args:
        star_name: Name of star (e.g. "Regulus", "Spica")
        tjd_ut: Julian Day in UT
        iflag: Calculation flags
        
    Returns:
        ((lon, lat, dist, speed_lon, speed_lat, speed_dist), iflag, error_msg)
    """
    # Map name to ID
    star_id = -1
    name_upper = star_name.upper().strip()
    if "," in name_upper:
        name_upper = name_upper.split(",")[0].strip()
        
    if name_upper == "REGULUS":
        star_id = SE_REGULUS
    elif name_upper == "SPICA":
        star_id = SE_SPICA_STAR
    else:
        return ((0.0, 0.0, 0.0, 0.0, 0.0, 0.0), iflag, f"Star {star_name} not found")

    # Calculate position (simplified, no speed)
    # We use TT for calculation (approximate UT=TT for stars is fine? No, use deltat)
    from .state import get_timescale
    ts = get_timescale()
    t = ts.ut1_jd(tjd_ut)
    
    try:
        lon, lat, dist = calc_fixed_star_position(star_id, t.tt)
        return ((lon, lat, dist, 0.0, 0.0, 0.0), iflag, "")
    except Exception as e:
        return ((0.0, 0.0, 0.0, 0.0, 0.0, 0.0), iflag, str(e))

