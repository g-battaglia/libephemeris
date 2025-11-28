"""
Utility functions for libephemeris.

Provides helper functions compatible with pyswisseph API including
angular calculations and other mathematical utilities.
"""


def difdeg2n(p1: float, p2: float) -> float:
    """
    Calculate distance in degrees p1 - p2 normalized to [-180;180].
    
    Compatible with pyswisseph's swe.difdeg2n() function.
    Computes the signed angular difference, handling 360° wrapping.
    
    Args:
        p1: First angle in degrees
        p2: Second angle in degrees
        
    Returns:
        Normalized difference in range [-180, 180]
        
    Examples:
        >>> difdeg2n(10, 20)
        -10.0
        >>> difdeg2n(350, 10)
        -20.0
        >>> difdeg2n(10, 350)
        20.0
        >>> difdeg2n(180, 0)
        180.0
    """
    diff = (p1 - p2) % 360.0
    if diff > 180.0:
        diff -= 360.0
    return diff


def swe_calc_angles(jd_ut: float, lat: float, lon: float):
    """
    Pre-calculate and cache astrological angles and planet positions
    for use with Arabic parts.

    Args:
        jd_ut: Julian Day (UT)
        lat: Latitude (degrees)
        lon: Longitude (degrees)

    Returns:
        Dictionary with calculated positions
    """
    from .state import set_angles_cache, set_topo
    from .angles import calc_angles
    from .planets import swe_calc_ut
    from .constants import SE_SUN, SE_MOON, SE_MERCURY, SE_VENUS

    # Set observer location
    set_topo(lon, lat, 0)

    # Calculate angles
    angles_dict = calc_angles(jd_ut, lat, lon)
    
    # Calculate and add planet positions for Arabic parts
    sun_pos, _ = swe_calc_ut(jd_ut, SE_SUN, 0)
    moon_pos, _ = swe_calc_ut(jd_ut, SE_MOON, 0)
    mercury_pos, _ = swe_calc_ut(jd_ut, SE_MERCURY, 0)
    venus_pos, _ = swe_calc_ut(jd_ut, SE_VENUS, 0)
    
    angles_dict['Sun'] = sun_pos[0]
    angles_dict['Moon'] = moon_pos[0]
    angles_dict['Mercury'] = mercury_pos[0]
    angles_dict['Venus'] = venus_pos[0]
    
    # Cache for Arabic parts
    set_angles_cache(angles_dict)
    
    return angles_dict
