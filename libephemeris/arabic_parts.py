"""
Arabic Parts (Lots) calculations for libephemeris.

Arabic Parts are classical astrological points calculated from combinations
of planetary positions and angles. They represent specific life areas and
are used in Hellenistic, Medieval, and modern astrology.

Formulas follow traditional methods from:
- Dorotheus of Sidon (1st century CE)
- Vettius Valens (2nd century CE)  
- Al-Biruni "Book of Instruction" (11th century)
- Robert Schmidt translations (Project Hindsight)

Standard formula structure: ASC + [Point A] - [Point B]
Many parts have day/night variations for sect consideration.

Supported Parts:
- Part of Fortune (Pars Fortunae): Body, health, material success
- Part of Spirit (Pars Spiritus): Soul, intellect, spiritual matters
- Part of Love (Pars Amoris/Eros): Romance, desire, relationships
- Part of Faith (Pars Fidei): Belief, religion, trust
"""

from typing import Dict
from .constants import (
    SE_PARS_FORTUNAE,
    SE_PARS_SPIRITUS,
    SE_PARS_AMORIS,
    SE_PARS_FIDEI,
)


def calc_arabic_part_of_fortune(
    asc: float, sun: float, moon: float, is_diurnal: bool = True
) -> float:
    """
    Calculate Part of Fortune (Pars Fortunae / Lot of Fortune).
    
    The most important Arabic Part, representing body, health, and material fortune.
    
    Formula (Vettius Valens, Al-Biruni):
        - Day chart (Sun above horizon): ASC + Moon - Sun
        - Night chart (Sun below horizon): ASC + Sun - Moon
        
    Args:
        asc: Ascendant (ecliptic) longitude in degrees
        sun: Sun ecliptic longitude in degrees
        moon: Moon ecliptic longitude in degrees
        is_diurnal: True if day chart (Sun above horizon), False if night
        
    Returns:
        float: Part of Fortune ecliptic longitude in degrees (0-360)
        
    Note:
        The formula reflects sect: in day charts, emphasize the Moon (nocturnal);
        in night charts, emphasize the Sun (diurnal). This balances opposites.
        
        Some modern astrologers use only the day formula regardless of sect.
        This implementation follows classical tradition.
    """
    if is_diurnal:
        lot = (asc + moon - sun) % 360.0
    else:
        lot = (asc + sun - moon) % 360.0

    return lot


def calc_arabic_part_of_spirit(
    asc: float, sun: float, moon: float, is_diurnal: bool = True
) -> float:
    """
    Calculate Part of Spirit (Pars Spiritus / Lot of Spirit / Daimon).
    
    Represents the soul, intellect, character, and spiritual development.
    
    Formula (opposite of Part of Fortune by sect):
        - Day chart: ASC + Sun - Moon
        - Night chart: ASC + Moon - Sun
        
    Args:
        asc: Ascendant longitude in degrees
        sun: Sun longitude in degrees
        moon: Moon longitude in degrees
        is_diurnal: True if day chart
        
    Returns:
        float: Part of Spirit longitude in degrees (0-360)
        
    Note:
        In Hellenistic astrology (Valens), Part of Spirit was called "Daimon"
        and considered complementary to Fortune. Together they represent
        body (Fortune) and soul (Spirit).
    """
    if is_diurnal:
        lot = (asc + sun - moon) % 360.0
    else:
        lot = (asc + moon - sun) % 360.0

    return lot


def calc_arabic_part_of_love(asc: float, venus: float, sun: float) -> float:
    """
    Calculate Part of Love (Pars Amoris / Lot of Eros).
    
    Represents romantic love, desire, sexual attraction, and relationships.
    
    Formula: ASC + Venus - Sun
    
    Args:
        asc: Ascendant longitude in degrees
        venus: Venus longitude in degrees
        sun: Sun longitude in degrees
        
    Returns:
        float: Part of Love longitude in degrees (0-360)
        
    Note:
        This formula is not sect-dependent. Venus naturally signifies
        love and attraction in all charts. Some medieval sources use
        alternate formulas with Moon or Mars for male/female charts.
    """
    return (asc + venus - sun) % 360.0


def calc_arabic_part_of_faith(asc: float, mercury: float, moon: float) -> float:
    """
    Calculate Part of Faith (Pars Fidei / Lot of Faith).
    
    Represents religious belief, trust, philosophical convictions.
    
    Formula: ASC + Mercury - Moon
    
    Args:
        asc: Ascendant longitude in degrees
        mercury: Mercury longitude in degrees
        moon: Moon longitude in degrees
        
    Returns:
        float: Part of Faith longitude in degrees (0-360)
        
    Note:
        Mercury represents rational thought and communication of beliefs.
        The Moon represents unconscious reception and emotional faith.
        This part shows how one's beliefs are communicated/manifested.
    """
    return (asc + mercury - moon) % 360.0


def is_day_chart(sun_lon: float, asc: float) -> bool:
    """
    Determine if chart is diurnal (day) or nocturnal (night) based on sect.
    
    A chart is diurnal if the Sun is above the horizon (in houses 7-12).
    
    Args:
        sun_lon: Sun ecliptic longitude in degrees (0-360)
        asc: Ascendant ecliptic longitude in degrees (0-360)
        
    Returns:
        bool: True if day chart (Sun above horizon), False if night chart
        
    Algorithm:
        The horizon runs from Ascendant (east) to Descendant (west).
        Points between ASC and DSC (going counter-clockwise through MC)
        are above the horizon. This is the zodiacal arc from ASC to ASC+180°.
        
    Note:
        This is a simplified 2D calculation using ecliptic longitude only.
        It assumes the horizon plane intersects the ecliptic at ASC/DSC.
        
        For extreme latitudes or precise calculations, use 3D horizon
        coordinates (altitude/azimuth). This method is traditional and
        sufficient for most astrological purposes.
    """
    desc = (asc + 180.0) % 360.0

    # Check if Sun is in upper hemisphere (ASC to DSC counter-clockwise)
    if asc < desc:
        # Normal case: e.g., ASC=0°, DSC=180°
        return asc <= sun_lon <= desc
    else:
        # Wrapped case: e.g., ASC=350°, DSC=170°
        # Sun is above if >= ASC (e.g., 350-360°) OR <= DSC (e.g., 0-170°)
        return sun_lon >= asc or sun_lon <= desc


def calc_all_arabic_parts(positions: Dict[str, float]) -> Dict[str, float]:
    """
    Calculate all standard Arabic parts from a position dictionary.
    
    Args:
        positions: Dictionary of celestial positions in ecliptic longitude degrees.
                  Required keys: 'Asc', 'Sun', 'Moon', 'Mercury', 'Venus'
                  All values should be in range 0-360 degrees.
                  
    Returns:
        Dict[str, float]: Dictionary mapping part names to longitudes:
            - "Pars_Fortunae": Part of Fortune
            - "Pars_Spiritus": Part of Spirit
            - "Pars_Amoris": Part of Love
            - "Pars_Fidei": Part of Faith
            
    Example:
        >>> positions = {
        ...     'Asc': 15.5, 'Sun': 120.0, 'Moon': 240.0,
        ...     'Mercury': 130.0, 'Venus': 110.0
        ... }
        >>> parts = calc_all_arabic_parts(positions)
        >>> print(parts['Pars_Fortunae'])
        135.5
        
    Note:
        Missing keys will default to 0.0. Ensure all required positions
        are present for accurate results.
    """
    asc = positions.get("Asc", 0.0)
    sun = positions.get("Sun", 0.0)
    moon = positions.get("Moon", 0.0)
    mercury = positions.get("Mercury", 0.0)
    venus = positions.get("Venus", 0.0)

    is_diurnal = is_day_chart(sun, asc)

    return {
        "Pars_Fortunae": calc_arabic_part_of_fortune(asc, sun, moon, is_diurnal),
        "Pars_Spiritus": calc_arabic_part_of_spirit(asc, sun, moon, is_diurnal),
        "Pars_Amoris": calc_arabic_part_of_love(asc, venus, sun),
        "Pars_Fidei": calc_arabic_part_of_faith(asc, mercury, moon),
    }

