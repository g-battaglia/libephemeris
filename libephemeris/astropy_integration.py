# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Astropy integration evaluation for libephemeris.

This module evaluates how astropy.coordinates and astropy.time could supplement
calculations in libephemeris. It provides comparison utilities to validate
libephemeris calculations against astropy's implementations.

Key Areas of Integration:
-------------------------

1. astropy.time:
   - Extended time scale support (TAI, TDB, TCG, TCB, UT1, UTC, GPS, etc.)
   - High-precision leap second handling
   - Time format parsing (ISO, JD, MJD, Unix, etc.)
   - Delta T calculations using IERS data
   - Comparison with Skyfield's time handling

2. astropy.coordinates:
   - Coordinate frame transformations (ICRS, GCRS, ITRS, etc.)
   - Geodetic/geocentric coordinate conversions
   - Proper motion and parallax corrections
   - Precession and nutation models
   - Barycentric corrections

Usage:
------
    >>> from libephemeris.astropy_integration import (
    ...     check_astropy_available,
    ...     compare_time_conversions,
    ...     compare_coordinate_transforms,
    ... )
    >>> if check_astropy_available():
    ...     results = compare_time_conversions(jd=2451545.0)
    ...     print(results)

Note:
-----
Astropy is an optional dependency. The functions in this module will raise
ImportError if astropy is not installed. Install with:
    pip install astropy

This module is designed for evaluation purposes to determine how astropy
could supplement (not replace) the existing Skyfield-based calculations.
"""

from __future__ import annotations

from typing import Any


def check_astropy_available() -> bool:
    """
    Check if astropy is available for use.

    Returns:
        bool: True if astropy is installed and can be imported, False otherwise.

    Example:
        >>> from libephemeris.astropy_integration import check_astropy_available
        >>> if check_astropy_available():
        ...     print("Astropy is available!")
    """
    try:
        import astropy.time
        import astropy.coordinates  # noqa: F401 (availability probe)

        return True
    except ImportError:
        return False


def _require_astropy() -> None:
    """Raise ImportError if astropy is not available."""
    if not check_astropy_available():
        raise ImportError(
            "Astropy is required for this function. Install with: pip install astropy"
        )


# =============================================================================
# TIME SCALE CONVERSIONS
# =============================================================================


def compare_time_conversions(jd_utc: float) -> dict[str, Any]:
    """
    Compare time scale conversions between libephemeris/Skyfield and astropy.

    This function demonstrates how astropy.time could supplement time handling
    by comparing Delta T calculations and time scale conversions.

    Args:
        jd_utc: Julian Day in UTC.

    Returns:
        dict: Comparison results with keys:
            - 'jd_utc': Input JD (UTC)
            - 'skyfield_delta_t_seconds': Delta T from Skyfield (seconds)
            - 'astropy_delta_t_seconds': Delta T from Astropy (seconds)
            - 'delta_t_difference_ms': Difference in milliseconds
            - 'skyfield_tt': JD in TT from Skyfield
            - 'astropy_tt': JD in TT from Astropy
            - 'tt_difference_seconds': Difference in TT (seconds)
            - 'astropy_tai': JD in TAI from Astropy
            - 'astropy_tdb': JD in TDB from Astropy
            - 'astropy_tcg': JD in TCG from Astropy
            - 'astropy_tcb': JD in TCB from Astropy

    Raises:
        ImportError: If astropy is not installed.

    Note:
        astropy provides additional time scales (TDB, TCG, TCB) that are not
        directly available in Skyfield. These are useful for:
        - TDB (Barycentric Dynamical Time): Solar system barycenter calculations
        - TCG (Geocentric Coordinate Time): Relativistic corrections at Earth
        - TCB (Barycentric Coordinate Time): Relativistic corrections at SSB

    Example:
        >>> from libephemeris.astropy_integration import compare_time_conversions
        >>> results = compare_time_conversions(2451545.0)  # J2000.0
        >>> print(f"Delta T difference: {results['delta_t_difference_ms']:.2f} ms")
    """
    _require_astropy()

    from astropy.time import Time
    from .state import get_timescale

    # Scale caveat: the input jd_utc is handed to Skyfield as UT1 and to
    # astropy as both UTC and UT1.  |UT1 - UTC| is kept below 0.9 s by
    # leap seconds, so the comparison values below carry an inherent
    # uncertainty of up to ~0.9 s (delta_t_difference_ms reflects the
    # two libraries' delta-T tables *plus* this scale conflation).  For
    # exact work convert the input explicitly before comparing.
    ts = get_timescale()
    t_sf = ts.ut1_jd(jd_utc)
    skyfield_delta_t = float(t_sf.delta_t)  # seconds
    skyfield_tt = float(t_sf.tt)

    # Astropy calculations
    t_ap = Time(jd_utc, format="jd", scale="utc")

    # Delta T from astropy (TT - UT1), treating the input as UT1
    t_ut1 = Time(jd_utc, format="jd", scale="ut1")
    astropy_delta_t = (t_ap.tt.jd - t_ut1.jd) * 86400.0  # seconds

    # Time scale conversions
    astropy_tt = t_ap.tt.jd
    astropy_tai = t_ap.tai.jd
    astropy_tdb = t_ap.tdb.jd
    astropy_tcg = t_ap.tcg.jd
    astropy_tcb = t_ap.tcb.jd

    return {
        "jd_utc": jd_utc,
        "skyfield_delta_t_seconds": skyfield_delta_t,
        "astropy_delta_t_seconds": astropy_delta_t,
        "delta_t_difference_ms": (skyfield_delta_t - astropy_delta_t) * 1000,
        "skyfield_tt": skyfield_tt,
        "astropy_tt": astropy_tt,
        "tt_difference_seconds": (skyfield_tt - astropy_tt) * 86400,
        "astropy_tai": astropy_tai,
        "astropy_tdb": astropy_tdb,
        "astropy_tcg": astropy_tcg,
        "astropy_tcb": astropy_tcb,
    }


def get_extended_time_scales(jd_utc: float) -> dict[str, float]:
    """
    Get time values in all scales available via astropy.

    This demonstrates time scales that astropy provides beyond what Skyfield offers.

    Args:
        jd_utc: Julian Day in UTC.

    Returns:
        dict: Julian Day values in various time scales:
            - 'utc': UTC (input)
            - 'ut1': UT1 (Earth rotation time)
            - 'tai': TAI (International Atomic Time)
            - 'tt': TT (Terrestrial Time)
            - 'tdb': TDB (Barycentric Dynamical Time)
            - 'tcg': TCG (Geocentric Coordinate Time)
            - 'tcb': TCB (Barycentric Coordinate Time)
            - 'gps': GPS time

    Raises:
        ImportError: If astropy is not installed.

    Example:
        >>> from libephemeris.astropy_integration import get_extended_time_scales
        >>> scales = get_extended_time_scales(2451545.0)
        >>> for name, jd in scales.items():
        ...     print(f"{name}: {jd:.10f}")
    """
    _require_astropy()

    from astropy.time import Time

    t = Time(jd_utc, format="jd", scale="utc")

    # Note: t.gps is a scalar (seconds since GPS epoch), not a Time object
    # We need to create a GPS Time separately to get JD
    gps_time = Time(t.gps, format="gps")

    return {
        "utc": float(t.utc.jd),
        "ut1": float(t.ut1.jd),
        "tai": float(t.tai.jd),
        "tt": float(t.tt.jd),
        "tdb": float(t.tdb.jd),
        "tcg": float(t.tcg.jd),
        "tcb": float(t.tcb.jd),
        "gps": float(gps_time.jd),
    }


def parse_time_string(time_string: str, format: str | None = None) -> float:
    """
    Parse a time string using astropy's flexible time parsing.

    This demonstrates astropy's time format parsing capabilities that could
    supplement libephemeris's time utilities.

    Args:
        time_string: Time string to parse (e.g., "2000-01-01T12:00:00").
        format: Optional format hint (e.g., "isot", "iso", "fits").
               If None, astropy will attempt to auto-detect.

    Returns:
        float: Julian Day in UTC.

    Raises:
        ImportError: If astropy is not installed.
        ValueError: If the time string cannot be parsed.

    Supported formats:
        - ISO: "2000-01-01 12:00:00" or "2000-01-01T12:00:00"
        - FITS: "2000-01-01T12:00:00.000"
        - Unix timestamp: seconds since 1970-01-01
        - MJD: Modified Julian Day
        - GPS: GPS seconds

    Note:
        Auto-detection (format=None) only recognizes the string-based date
        formats (ISO, FITS). The numeric formats (unix, mjd, gps) are NOT
        auto-detected and require an explicit ``format=`` argument.

    Example:
        >>> from libephemeris.astropy_integration import parse_time_string
        >>> jd = parse_time_string("2000-01-01T12:00:00")
        >>> print(f"JD: {jd:.6f}")  # ~2451545.0
    """
    _require_astropy()

    from astropy.time import Time

    if format:
        t = Time(time_string, format=format, scale="utc")
    else:
        t = Time(time_string, scale="utc")

    return float(t.jd)


# =============================================================================
# COORDINATE TRANSFORMATIONS
# =============================================================================


def compare_coordinate_transforms(
    ra_deg: float, dec_deg: float, jd_utc: float, lon_deg: float, lat_deg: float
) -> dict[str, Any]:
    """
    Compare coordinate transformations between libephemeris/Skyfield and astropy.

    This demonstrates how astropy.coordinates could supplement coordinate
    transformations in libephemeris.

    Args:
        ra_deg: Right Ascension in degrees (ICRS J2000).
        dec_deg: Declination in degrees (ICRS J2000).
        jd_utc: Julian Day in UTC for the observation.
        lon_deg: Observer longitude in degrees (East positive).
        lat_deg: Observer latitude in degrees (North positive).

    Returns:
        dict: Comparison results including:
            - 'icrs_ra', 'icrs_dec': Input ICRS coordinates
            - 'skyfield_alt', 'skyfield_az': Altitude/Azimuth from Skyfield
            - 'astropy_alt', 'astropy_az': Altitude/Azimuth from Astropy
            - 'alt_difference_arcsec': Difference in altitude (arcsec)
            - 'az_difference_arcsec': Difference in azimuth (arcsec)
            - 'astropy_galactic_l', 'astropy_galactic_b': Galactic coords
            - 'astropy_ecliptic_lon', 'astropy_ecliptic_lat': Ecliptic coords

    Raises:
        ImportError: If astropy is not installed.

    Note:
        astropy provides coordinate frames not directly available in Skyfield:
        - Galactic coordinates (l, b)
        - Supergalactic coordinates
        - Various ecliptic frames (BarycentricMeanEcliptic, etc.)
        - ITRS (geocentric terrestrial reference)

    Example:
        >>> from libephemeris.astropy_integration import compare_coordinate_transforms
        >>> # Compare for Sirius at Rome
        >>> results = compare_coordinate_transforms(
        ...     ra_deg=101.2875, dec_deg=-16.7161,
        ...     jd_utc=2451545.0, lon_deg=12.5, lat_deg=41.9
        ... )
        >>> print(f"Alt diff: {results['alt_difference_arcsec']:.2f} arcsec")
    """
    _require_astropy()

    from astropy.time import Time
    from astropy.coordinates import (
        SkyCoord,
        EarthLocation,
        AltAz,
        BarycentricMeanEcliptic,
    )
    import astropy.units as u
    from skyfield.api import Star, wgs84
    from .state import get_planets, get_timescale

    # Skyfield calculation
    ts = get_timescale()
    t_sf = ts.ut1_jd(jd_utc)
    planets = get_planets()
    earth = planets["earth"]

    # Create star and observer
    star = Star(ra_hours=ra_deg / 15.0, dec_degrees=dec_deg)
    observer = earth + wgs84.latlon(lat_deg, lon_deg)

    # Calculate alt/az with Skyfield
    apparent = observer.at(t_sf).observe(star).apparent()
    alt_sf, az_sf, _ = apparent.altaz()
    skyfield_alt = alt_sf.degrees
    skyfield_az = az_sf.degrees

    # Astropy calculation
    t_ap = Time(jd_utc, format="jd", scale="utc")
    location = EarthLocation(lon=lon_deg * u.deg, lat=lat_deg * u.deg, height=0 * u.m)
    coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")

    # Transform to AltAz
    altaz_frame = AltAz(obstime=t_ap, location=location)
    altaz = coord.transform_to(altaz_frame)
    astropy_alt = altaz.alt.deg
    astropy_az = altaz.az.deg

    # Additional coordinate frames from astropy
    galactic = coord.galactic
    ecliptic = coord.transform_to(BarycentricMeanEcliptic(equinox=t_ap))

    return {
        "icrs_ra": ra_deg,
        "icrs_dec": dec_deg,
        "jd_utc": jd_utc,
        "observer_lon": lon_deg,
        "observer_lat": lat_deg,
        "skyfield_alt": skyfield_alt,
        "skyfield_az": skyfield_az,
        "astropy_alt": astropy_alt,
        "astropy_az": astropy_az,
        "alt_difference_arcsec": (skyfield_alt - astropy_alt) * 3600,
        # Wrap the azimuth delta to (-180, 180] so a star near due north
        # (azimuths straddling 360/0) does not report a ~360 deg difference.
        "az_difference_arcsec": ((skyfield_az - astropy_az + 180.0) % 360.0 - 180.0)
        * 3600,
        "astropy_galactic_l": galactic.l.deg,
        "astropy_galactic_b": galactic.b.deg,
        "astropy_ecliptic_lon": ecliptic.lon.deg,
        "astropy_ecliptic_lat": ecliptic.lat.deg,
    }


def icrs_to_galactic(ra_deg: float, dec_deg: float) -> tuple[float, float]:
    """
    Convert ICRS (J2000) coordinates to Galactic coordinates using astropy.

    This is a coordinate transformation not directly available in Skyfield.

    Args:
        ra_deg: Right Ascension in degrees (ICRS).
        dec_deg: Declination in degrees (ICRS).

    Returns:
        tuple: (galactic_l, galactic_b) in degrees.
            - galactic_l: Galactic longitude (0-360)
            - galactic_b: Galactic latitude (-90 to +90)

    Raises:
        ImportError: If astropy is not installed.

    Example:
        >>> from libephemeris.astropy_integration import icrs_to_galactic
        >>> l, b = icrs_to_galactic(266.417, -29.008)  # Near Galactic center
        >>> print(f"l={l:.2f}, b={b:.2f}")  # ~(0, 0)
    """
    _require_astropy()

    from astropy.coordinates import SkyCoord
    import astropy.units as u

    coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
    galactic = coord.galactic

    return float(galactic.l.deg), float(galactic.b.deg)


def galactic_to_icrs(l_deg: float, b_deg: float) -> tuple[float, float]:
    """
    Convert Galactic coordinates to ICRS (J2000) using astropy.

    Args:
        l_deg: Galactic longitude in degrees.
        b_deg: Galactic latitude in degrees.

    Returns:
        tuple: (ra, dec) in degrees (ICRS).

    Raises:
        ImportError: If astropy is not installed.

    Example:
        >>> from libephemeris.astropy_integration import galactic_to_icrs
        >>> ra, dec = galactic_to_icrs(0, 0)  # Galactic center
        >>> print(f"RA={ra:.2f}, Dec={dec:.2f}")
    """
    _require_astropy()

    from astropy.coordinates import SkyCoord
    import astropy.units as u

    coord = SkyCoord(l=l_deg * u.deg, b=b_deg * u.deg, frame="galactic")
    icrs = coord.icrs

    return float(icrs.ra.deg), float(icrs.dec.deg)


def get_barycentric_correction(
    ra_deg: float, dec_deg: float, jd_utc: float, lon_deg: float, lat_deg: float
) -> dict[str, float]:
    """
    Calculate barycentric corrections using astropy.

    Barycentric corrections are important for high-precision radial velocity
    measurements and pulsar timing. This demonstrates astropy's capability
    for these specialized calculations.

    Args:
        ra_deg: Right Ascension in degrees (ICRS).
        dec_deg: Declination in degrees (ICRS).
        jd_utc: Julian Day in UTC.
        lon_deg: Observer longitude in degrees.
        lat_deg: Observer latitude in degrees.

    Returns:
        dict: Barycentric correction data:
            - 'bjd_tdb': Barycentric Julian Date in TDB
            - 'light_time_seconds': Light travel time correction (seconds)
            - 'rv_correction_km_s': Radial velocity correction (km/s)

    Raises:
        ImportError: If astropy is not installed.

    Note:
        These corrections account for:
        - Earth's orbital motion around the Sun
        - Earth's rotation
        - Relativistic effects (time dilation)
        - Light travel time to the Solar System barycenter

    Example:
        >>> from libephemeris.astropy_integration import get_barycentric_correction
        >>> result = get_barycentric_correction(
        ...     ra_deg=180.0, dec_deg=45.0,
        ...     jd_utc=2451545.0, lon_deg=0.0, lat_deg=50.0
        ... )
        >>> print(f"Light time: {result['light_time_seconds']:.2f} s")
    """
    _require_astropy()

    from astropy.time import Time
    from astropy.coordinates import SkyCoord, EarthLocation
    import astropy.units as u

    t = Time(jd_utc, format="jd", scale="utc")
    location = EarthLocation(lon=lon_deg * u.deg, lat=lat_deg * u.deg)
    coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")

    # Calculate barycentric correction
    # Light travel time to barycenter
    ltt = t.light_travel_time(coord, kind="barycentric", location=location)

    # BJD in TDB
    bjd_tdb = (t.tdb + ltt).jd

    # Radial velocity correction
    rv_corr = coord.radial_velocity_correction(
        obstime=t, location=location, kind="barycentric"
    )

    return {
        "bjd_tdb": float(bjd_tdb),
        "light_time_seconds": float(ltt.to(u.s).value),
        "rv_correction_km_s": float(rv_corr.to(u.km / u.s).value),
    }


# =============================================================================
# GEODETIC UTILITIES
# =============================================================================


def geodetic_to_geocentric(
    lon_deg: float, lat_deg: float, height_m: float
) -> tuple[float, float, float]:
    """
    Convert geodetic coordinates to geocentric Cartesian using astropy.

    Geodetic coordinates are based on the WGS84 ellipsoid (latitude is the
    angle from the equatorial plane along the surface normal). Geocentric
    coordinates are Cartesian with origin at Earth's center.

    Args:
        lon_deg: Geodetic longitude in degrees.
        lat_deg: Geodetic latitude in degrees.
        height_m: Height above the WGS84 ellipsoid in meters.

    Returns:
        tuple: (x, y, z) in meters (ITRS geocentric Cartesian).

    Raises:
        ImportError: If astropy is not installed.

    Example:
        >>> from libephemeris.astropy_integration import geodetic_to_geocentric
        >>> x, y, z = geodetic_to_geocentric(0, 45, 0)  # On ellipsoid at 45N
        >>> print(f"x={x/1e6:.2f} Mm, z={z/1e6:.2f} Mm")
    """
    _require_astropy()

    from astropy.coordinates import EarthLocation
    import astropy.units as u

    loc = EarthLocation(lon=lon_deg * u.deg, lat=lat_deg * u.deg, height=height_m * u.m)

    return (
        float(loc.x.to(u.m).value),
        float(loc.y.to(u.m).value),
        float(loc.z.to(u.m).value),
    )


def geocentric_to_geodetic(
    x_m: float, y_m: float, z_m: float
) -> tuple[float, float, float]:
    """
    Convert geocentric Cartesian to geodetic coordinates using astropy.

    Args:
        x_m: X coordinate in meters (ITRS).
        y_m: Y coordinate in meters (ITRS).
        z_m: Z coordinate in meters (ITRS).

    Returns:
        tuple: (longitude, latitude, height) where:
            - longitude in degrees
            - latitude in degrees (geodetic)
            - height in meters above WGS84 ellipsoid

    Raises:
        ImportError: If astropy is not installed.

    Example:
        >>> from libephemeris.astropy_integration import geocentric_to_geodetic
        >>> lon, lat, h = geocentric_to_geodetic(4e6, 0, 5e6)
        >>> print(f"Lon={lon:.2f}, Lat={lat:.2f}, h={h:.0f} m")
    """
    _require_astropy()

    from astropy.coordinates import EarthLocation
    import astropy.units as u

    loc = EarthLocation.from_geocentric(x_m * u.m, y_m * u.m, z_m * u.m)

    return (
        float(loc.lon.deg),
        float(loc.lat.deg),
        float(loc.height.to(u.m).value),
    )


# =============================================================================
# EVALUATION SUMMARY
# =============================================================================


def evaluate_integration_potential() -> dict[str, Any]:
    """
    Evaluate the potential benefits of astropy integration.

    This function provides a comprehensive evaluation of what astropy.coordinates
    and astropy.time could add to libephemeris.

    Returns:
        dict: Evaluation summary with keys:
            - 'astropy_available': Whether astropy is installed
            - 'time_features': List of time-related features astropy provides
            - 'coordinate_features': List of coordinate features astropy provides
            - 'recommended_integrations': Specific recommendations
            - 'compatibility_notes': Notes on Skyfield/astropy compatibility

    Example:
        >>> from libephemeris.astropy_integration import evaluate_integration_potential
        >>> evaluation = evaluate_integration_potential()
        >>> for feature in evaluation['time_features']:
        ...     print(f"- {feature}")
    """
    result: dict[str, Any] = {
        "astropy_available": check_astropy_available(),
        "time_features": [
            "Extended time scales: TDB, TCG, TCB, GPS (beyond Skyfield's TT/UT1/UTC)",
            "Flexible time string parsing (ISO, FITS, Unix, GPS seconds)",
            "High-precision leap second handling via IERS data",
            "Light travel time calculations for barycentric corrections",
            "BJD (Barycentric Julian Date) computation",
        ],
        "coordinate_features": [
            "Galactic coordinate system (l, b)",
            "Supergalactic coordinate system",
            "Multiple ecliptic frames (Mean, True, Barycentric)",
            "ITRS (Earth-fixed geocentric frame)",
            "Proper motion propagation with radial velocity",
            "Barycentric/heliocentric velocity corrections",
            "WGS84 ellipsoid geodetic calculations",
        ],
        "recommended_integrations": [
            "Optional Galactic coordinate conversions (not in Skyfield)",
            "TDB time scale for SSB calculations (more accurate than TT)",
            "Barycentric corrections for high-precision work",
            "Time string parsing as user convenience",
            "Geodetic-geocentric conversions using proper ellipsoid model",
        ],
        "compatibility_notes": [
            "Skyfield and astropy use similar underlying data (IERS, JPL)",
            "Both use ICRS as fundamental reference frame",
            "Delta T values may differ slightly due to different models",
            "Nutation models may have small differences (~0.1 arcsec)",
            "For most astrological uses, differences are negligible",
        ],
        "implementation_approach": [
            "Make astropy optional (install via 'pip install libephemeris[astropy]')",
            "Provide astropy-based alternatives in this module",
            "Use Skyfield for core calculations (performance, existing codebase)",
            "Use astropy for specialized features (Galactic coords, barycentric)",
            "Validate consistency between implementations via tests",
        ],
    }

    return result
