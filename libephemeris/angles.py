"""
Astrological angles and chart points for libephemeris.

This module calculates the primary astrological angles:
- Ascendant (ASC): Ecliptic degree rising on eastern horizon  
- Midheaven (MC): Ecliptic degree culminating on meridian
- Descendant (DSC): Western horizon point (ASC + 180°)
- Imum Coeli (IC): Lower meridian (MC + 180°)
- Vertex: Intersection of prime vertical with western ecliptic
- Equatorial Ascendant: Equator crossing eastern horizon

All angles require observer geographic location (swe_set_topo).
Calculations are based on spherical astronomy and house system cusps.

Note:
    Angles are independent of house system choice (though computed via houses).
    They represent fundamental horizon/meridian intersections with the ecliptic.
"""

from typing import Dict
from .constants import (
    SE_ASCENDANT,
    SE_MC,
    SE_DESCENDANT,
    SE_IC,
    SE_VERTEX,
    SE_ANTIVERTEX,
)
from .houses import swe_houses


def calc_angles(jd_ut: float, lat: float, lon: float) -> Dict[str, float]:
    """
    Calculate all astrological angles for a given time and location.
    
    Args:
        jd_ut: Julian Day in Universal Time (UT1)
        lat: Geographic latitude in degrees (North positive)
        lon: Geographic longitude in degrees (East positive)
        
    Returns:
        Dict[str, float]: Dictionary mapping angle names to ecliptic longitudes:
            - "Ascendant": Rising degree on eastern horizon
            - "MC": Midheaven (culminating degree)
            - "ARMC": Right Ascension of MC in degrees
            - "Vertex": Prime vertical intersection (west)
            - "Equatorial_Ascendant": Equator on eastern horizon
            - "Descendant": Setting degree (Asc + 180°)
            - "IC": Lower meridian (MC + 180°)
            - "AntiVertex": Prime vertical intersection (east)
            
    Note:
        Uses Placidus house system internally, but angle values are
        system-independent. They are purely geometric intersections.
        
        Angles are sensitive to observer location. Polar regions may
        have undefined values for some angles.
    """
    _, ascmc = swe_houses(jd_ut, lat, lon, ord("P"))

    return {
        "Ascendant": ascmc[0],
        "MC": ascmc[1],
        "ARMC": ascmc[2],
        "Vertex": ascmc[3],
        "Equatorial_Ascendant": ascmc[4] if len(ascmc) > 4 else 0.0,
        "Descendant": (ascmc[0] + 180.0) % 360.0,
        "IC": (ascmc[1] + 180.0) % 360.0,
        "AntiVertex": (ascmc[3] + 180.0) % 360.0,
    }


def get_angle_value(angle_id: int, jd_ut: float, lat: float, lon: float) -> float:
    """
    Get ecliptic longitude for a specific angle.
    
    Args:
        angle_id: Angle identifier constant (SE_ASCENDANT, SE_MC, etc.)
        jd_ut: Julian Day in Universal Time (UT1)
        lat: Geographic latitude in degrees (North positive)
        lon: Geographic longitude in degrees (East positive)
        
    Returns:
        float: Ecliptic longitude of the angle in degrees (0-360)
               Returns 0.0 if angle_id is not recognized
               
    Example:
        >>> asc = get_angle_value(SE_ASCENDANT, jd, 41.9, 12.5)
        >>> mc = get_angle_value(SE_MC, jd, 41.9, 12.5)
    """
    angles = calc_angles(jd_ut, lat, lon)

    angle_map = {
        SE_ASCENDANT: "Ascendant",
        SE_MC: "MC",
        SE_DESCENDANT: "Descendant",
        SE_IC: "IC",
        SE_VERTEX: "Vertex",
        SE_ANTIVERTEX: "AntiVertex",
    }

    if angle_id in angle_map:
        return angles[angle_map[angle_id]]

    return 0.0

