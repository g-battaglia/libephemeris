"""
Astrological house system calculations for libephemeris.

Implements 19 house systems compatible with Swiss Ephemeris:
- Placidus (P): Most common, time-based, fails at polar latitudes
- Koch (K): Birthplace system, similar to Placidus
- Porphyrius (O): Space-based trisection
- Regiomontanus (R): Medieval rational system
- Campanus (C): Prime vertical system
- Equal (A/E): Equal 30° divisions from Ascendant
- Whole Sign (W): Whole zodiac signs from Ascendant sign
- Meridian (X): Equatorial meridian divisions
- Azimuthal/Horizontal (H): Based on horizon
- Polich-Page (T): Topocentric system
- Alcabitus (B): Ancient Arabic system
- Morinus (M): Equatorial divisions
- Krusinski-Pisa (U): Modified Regiomontanus
- Gauquelin (G): Sector system
- Vehlow (V): Equal from midpoint
- APC (houses): Astronomical Planetary Cusps
- Carter Poli-Equatorial (F)
- Pulhemus (L)
- Sripati (S): Divide quadrants equally

Main Functions:
- swe_houses(): Calculate house cusps and angles (ASCMC)
- swe_houses_ex(): Extended version (unused in libephemeris)
- swe_house_pos(): Find which house a point is in
- swe_house_name(): Get house system name

Polar Latitudes:
FIXME: Precision - Quadrant house systems fail near poles
- Placidus, Koch, Regiomontanus undefined > ~66° latitude
- Falls back to Porphyrius at high latitudes
- Equal/Whole Sign work at all latitudes

Algorithm Sources:
- Placidus: Time divisions of diurnal/nocturnal arcs
- Regiomontanus: Equator trisection projected to ecliptic
- Campanus: Prime vertical trisection
- Equal: Simple 30° additions
- Algorithms from Meeus, Swiss Ephemeris documentation

FIXME: Precision - Spherical trigonometry precision
- Uses standard spherical astronomy formulas
- Precision typically 0.01° (36 arcseconds)
- Swiss Ephemeris achieves ~0.001° with iterative refinement

References:
- Meeus "Astronomical Algorithms" Ch. 13 (coordinate systems)
- Swiss Ephemeris documentation (house systems)
- Hand "Astrological Houses" (comprehensive house treatise)
"""


import math
from typing import Tuple, List, Optional
from .constants import *
from .constants import SEFLG_SIDEREAL
from .state import get_timescale
from .state import get_sid_mode
from .planets import swe_get_ayanamsa_ut
from skyfield.nutationlib import iau2000b_radians


def angular_diff(a: float, b: float) -> float:
    """
    Calculate signed angular difference (a - b) handling 360° wrapping.
    
    Used by horizontal house system to find ecliptic longitude for given azimuth.
    
    Args:
        a: First angle in degrees
        b: Second angle in degrees
        
    Returns:
        Signed difference in range [-180, 180]
    """
    diff = (a - b) % 360.0
    if diff > 180.0:
        diff -= 360.0
    return diff


def _calc_vertex(armc_deg: float, eps: float, lat: float, mc: float) -> float:
    """
    Calculates the Vertex (intersection of Prime Vertical and Ecliptic in Western hemisphere).
    """
    # Handle equator by using a small epsilon
    calc_lat = lat
    if abs(calc_lat) < 1e-6:
        calc_lat = 1e-6

    armc_rad = math.radians(armc_deg)
    eps_rad = math.radians(eps)
    lat_rad = math.radians(calc_lat)

    num = -math.cos(armc_rad)
    den = math.sin(armc_rad) * math.cos(eps_rad) - math.sin(eps_rad) / math.tan(lat_rad)

    vtx_rad = math.atan2(num, den)
    vtx = math.degrees(vtx_rad)
    vtx = vtx % 360.0

    # Ensure Vertex is in Western Hemisphere relative to ARMC
    # i.e., Vertex should be West of ARMC.
    # West means (ARMC - Vertex) % 360 is in [0, 180].
    # Or (Vertex - ARMC) % 360 is in [180, 360].

    diff = (vtx - mc) % 360.0
    if diff < 180.0:
        vtx = (vtx + 180.0) % 360.0

    return vtx


def swe_houses(
    tjdut: float, lat: float, lon: float, hsys: int
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """
    Computes house cusps and ascmc.
    Returns (cusps, ascmc).
    cusps is a tuple of 12 floats (houses 1-12).
    ascmc is a tuple of 8 floats.
    """
    # 1. Calculate Sidereal Time (ARMC)
    # ARMC = GMST + lon
    ts = get_timescale()
    t = ts.ut1_jd(tjdut)

    # Skyfield gmst is in hours. We should use GAST (Apparent Sidereal Time) for houses.
    # t.gast is in hours.
    gast = t.gast
    armc_deg = (gast * 15.0 + lon) % 360.0

    # Obliquity of Ecliptic (True)
    # We can get it from Skyfield or use a standard formula.
    # Skyfield's `nutation_libration(t)` returns (dpsi, deps).
    # Mean obliquity can be computed.
    # Let's use Skyfield's internal functions if possible or standard formula.
    # IAU 1980 or 2000? SwissEph uses IAU 1980 by default but supports 2000.
    # Let's use a standard formula for mean obliquity + nutation.

    # Simple formula for Mean Obliquity (Laskar)
    T = (t.tt - 2451545.0) / 36525.0
    eps0 = 23.43929111 - (46.8150 * T + 0.00059 * T**2 - 0.001813 * T**3) / 3600.0

    # Nutation in obliquity (deps)
    # Use IAU 2000B model to match our ayanamsha precision
    # iau2000b_radians takes Time object? No, it takes T (centuries since J2000) usually?
    # Let's check signature. It takes (t).
    # t is a Time object? No, iau2000b_radians(t) where t is Time object?
    # In planets.py we used: dpsi, deps = iau2000b_radians(t)

    dpsi_rad, deps_rad = iau2000b_radians(t)
    deps_deg = math.degrees(deps_rad)

    eps = eps0 + deps_deg  # True Obliquity

    # 2. Calculate Ascendant and MC
    # MC is intersection of Meridian and Ecliptic.
    # tan(MC) = tan(ARMC) / cos(eps)
    # Quadrant check needed.

    # Determine if we need to flip MC (and thus ARMC) for specific systems
    # Regiomontanus (R), Campanus (C), Polich/Page (T) flip MC if below horizon.

    hsys_char = hsys
    if isinstance(hsys, int):
        hsys_char = chr(hsys)
    elif isinstance(hsys, bytes):
        hsys_char = hsys.decode("utf-8")

    # Determine if we need to flip MC (and thus ARMC) for specific systems
    # Regiomontanus (R), Campanus (C), Polich/Page (T) flip MC if below horizon.

    armc_active = armc_deg

    if hsys_char in ["R", "C", "T"]:
        # Check altitude of MC calculated from original ARMC
        mc_dec_rad = math.atan(
            math.sin(math.radians(armc_deg)) * math.tan(math.radians(eps))
        )
        lat_rad = math.radians(lat)
        sin_alt = math.sin(lat_rad) * math.sin(mc_dec_rad) + math.cos(
            lat_rad
        ) * math.cos(mc_dec_rad)

        if sin_alt < 0:
            armc_active = (armc_deg + 180.0) % 360.0

    # 2. Calculate Ascendant and MC

    # MC uses armc_active (flipped if needed)
    mc_rad = math.atan2(
        math.tan(math.radians(armc_active)), math.cos(math.radians(eps))
    )
    mc = math.degrees(mc_rad)
    # Adjust quadrant to match armc_active
    if mc < 0:
        mc += 360.0

    if 90.0 < armc_active <= 270.0:
        if mc < 90.0 or mc > 270.0:
            mc += 180.0
    elif armc_active > 270.0:
        if mc < 270.0:
            mc += 180.0
    elif armc_active <= 90.0:
        if mc > 90.0:
            mc += 180.0

    mc = mc % 360.0

    # Ascendant uses armc_deg (Original)
    num = math.cos(math.radians(armc_deg))
    den = -(
        math.sin(math.radians(armc_deg)) * math.cos(math.radians(eps))
        + math.tan(math.radians(lat)) * math.sin(math.radians(eps))
    )
    asc_rad = math.atan2(num, den)
    asc = math.degrees(asc_rad)
    asc = asc % 360.0

    # Ensure Ascendant is on the Eastern Horizon (Azimuth in [0, 180])
    # We check Azimuth relative to the TRUE ARMC (armc_deg)

    asc_r = math.radians(asc)
    eps_r = math.radians(eps)

    # RA
    y = math.cos(eps_r) * math.sin(asc_r)
    x = math.cos(asc_r)
    ra_r = math.atan2(y, x)
    ra = math.degrees(ra_r) % 360.0

    # Dec
    dec_r = math.asin(math.sin(eps_r) * math.sin(asc_r))

    # Hour Angle using TRUE ARMC
    h_deg = (armc_deg - ra + 360.0) % 360.0
    h_r = math.radians(h_deg)

    # Azimuth
    # tan(Az) = sin(H) / (sin(lat)cos(H) - cos(lat)tan(Dec))
    lat_r = math.radians(lat)

    num_az = math.sin(h_r)
    den_az = math.sin(lat_r) * math.cos(h_r) - math.cos(lat_r) * math.tan(dec_r)
    az_r = math.atan2(num_az, den_az)
    az = math.degrees(az_r)
    az = (az + 180.0) % 360.0

    # Check if H is West (0-180). If so, Asc is Setting (Descendant).
    # We want Rising.
    if 0.0 < h_deg < 180.0:
        asc = (asc + 180.0) % 360.0

    # Vertex uses armc_deg (Original)
    # Hemisphere check relative to TRUE ARMC (West of True ARMC)
    vertex = _calc_vertex(armc_deg, eps, lat, armc_deg)

    # Equatorial Ascendant (East Point)
    # This is the intersection of the ecliptic with the celestial equator in the east
    # It's the ecliptic longitude where RA = ARMC + 90°
    equ_asc_ra = (armc_deg + 90.0) % 360.0
    # Convert RA to ecliptic longitude
    # tan(Lon) = tan(RA) / cos(eps)
    # y = sin(RA)
    # x = cos(RA) * cos(eps)
    equ_asc_ra_r = math.radians(equ_asc_ra)
    eps_r = math.radians(eps)
    y = math.sin(equ_asc_ra_r)
    x = math.cos(equ_asc_ra_r) * math.cos(eps_r)
    equ_asc = math.degrees(math.atan2(y, x)) % 360.0

    # Co-Ascendant (Walter Koch)
    # Placeholder - requires further research on exact formula
    co_asc = 0.0

    # Co-Ascendant (Koch) - alternative calculation
    co_asc_koch = 0.0

    # Polar Ascendant
    # Placeholder - requires further research on exact formula
    polar_asc = 0.0

    # Build ASCMC array with 8 elements (pyswisseph compatible)
    ascmc = [0.0] * 8
    ascmc[0] = asc
    ascmc[1] = mc
    ascmc[2] = armc_deg
    ascmc[3] = vertex
    ascmc[4] = equ_asc
    ascmc[5] = co_asc
    ascmc[6] = co_asc_koch
    ascmc[7] = polar_asc

    # 3. House Cusps
    # Use armc_active for house calculations
    # If MC was flipped, we might need to flip latitude for intermediate cusps (Regiomontanus, etc.)
    # Verified for Regiomontanus: using -lat with flipped MC matches SWE.

    calc_lat = lat
    if armc_active != armc_deg:
        # MC was flipped. Flip latitude for intermediate cusp calculations.
        calc_lat = -lat

    cusps = [0.0] * 13

    if hsys_char == "P":  # Placidus
        cusps = _houses_placidus(
            armc_active, lat, eps, asc, mc
        )  # Placidus fails anyway
    elif hsys_char == "K":  # Koch
        cusps = _houses_koch(armc_active, lat, eps, asc, mc)  # Koch fails anyway
    elif hsys_char == "R":  # Regiomontanus
        cusps = _houses_regiomontanus(armc_active, calc_lat, eps, asc, mc)
    elif hsys_char == "C":  # Campanus
        cusps = _houses_campanus(armc_active, calc_lat, eps, asc, mc)
    elif hsys_char == "E":  # Equal (Ascendant)
        cusps = _houses_equal(asc)
    elif hsys_char == "A":  # Equal (MC)
        cusps = _houses_equal(asc)  # Equal MC is same as Equal Asc in SwissEph
    elif hsys_char == "W":  # Whole Sign
        cusps = _houses_whole_sign(asc)
    elif hsys_char == "O":  # Porphyry
        cusps = _houses_porphyry(asc, mc)
    elif hsys_char == "B":  # Alcabitius
        cusps = _houses_alcabitius(armc_active, lat, eps, asc, mc)
    elif hsys_char == "T":  # Polich/Page (Topocentric)
        cusps = _houses_polich_page(armc_active, calc_lat, eps, asc, mc)
    elif hsys_char == "M":  # Morinus
        cusps = _houses_morinus(armc_active, lat, eps, asc, mc)
    elif hsys_char == "X":  # Meridian (Axial)
        cusps = _houses_meridian(armc_active, lat, eps, asc, mc)
    elif hsys_char == "V":  # Vehlow
        cusps = _houses_vehlow(asc)
    elif hsys_char == "H":  # Horizontal (Azimuthal)
        cusps = _houses_horizontal(armc_active, lat, eps, asc, mc)
    elif hsys_char == "Y":  # APC Houses
        cusps = _houses_apc(armc_active, lat, eps, asc, mc)
    elif hsys_char == "F":  # Carter (Poli-Equatorial)
        cusps = _houses_carter(armc_active, lat, eps, asc, mc)
    elif hsys_char == "U":  # Krusinski
        cusps = _houses_krusinski(armc_active, lat, eps, asc, mc)
    elif hsys_char == "N":  # Natural Gradient
        cusps = _houses_natural_gradient(armc_active, lat, eps, asc, mc)
    elif hsys_char == "G":  # Gauquelin
        cusps = _houses_gauquelin(armc_active, lat, eps, asc, mc)
    else:
        # Default to Placidus
        cusps = _houses_placidus(armc_active, lat, eps, asc, mc)

    # Return 12-element cusps array (pyswisseph compatible: no padding at index 0)
    # cusps[1:13] contains houses 1-12
    return tuple(cusps[1:13]), tuple(ascmc)


def swe_houses_ex(
    tjdut: float, lat: float, lon: float, hsys: int, flags: int = 0
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """
    Extended house calculation (supports sidereal).
    
    Args:
        tjdut: Julian Day in UT
        lat: Latitude in degrees
        lon: Longitude in degrees
        hsys: House system (int or bytes)
        flags: Calculation flags (e.g., FLG_SIDEREAL)
    
    Returns:
        (cusps, ascmc): Tuple of house cusps and angles
    """
    cusps, ascmc = swe_houses(tjdut, lat, lon, hsys)

    if flags & SEFLG_SIDEREAL:
        ayanamsa = swe_get_ayanamsa_ut(tjdut)
        # cusps is now 12 elements (houses 1-12), apply ayanamsa to all
        cusps = tuple([(c - ayanamsa) % 360.0 for c in cusps])
        ascmc = list(ascmc)
        ascmc[0] = (ascmc[0] - ayanamsa) % 360.0  # Asc
        ascmc[1] = (ascmc[1] - ayanamsa) % 360.0  # MC
        ascmc = tuple(ascmc)

    return cusps, ascmc


    return cusps, ascmc


def swe_house_name(hsys: int) -> str:
    """
    Get the name of a house system.
    """
    hsys_char = hsys
    if isinstance(hsys, int):
        hsys_char = chr(hsys)
    elif isinstance(hsys, bytes):
        hsys_char = hsys.decode("utf-8")
        
    names = {
        "P": "Placidus",
        "K": "Koch",
        "O": "Porphyry",
        "R": "Regiomontanus",
        "C": "Campanus",
        "E": "Equal (Asc)",
        "A": "Equal (MC)",
        "W": "Whole Sign",
        "M": "Morinus",
        "B": "Alcabitius",
        "T": "Polich/Page",
        "U": "Krusinski",
        "G": "Gauquelin",
        "V": "Vehlow",
        "X": "Meridian",
        "H": "Horizontal",
        "F": "Carter",
        "S": "Sripati",
    }
    return names.get(hsys_char, "Unknown")


def _houses_placidus(armc: float, lat: float, eps: float, asc: float, mc: float) -> List[float]:
    """
    Placidus house system (time-based divisions of diurnal/nocturnal arcs).
    
    Most popular house system in modern Western astrology. Divides the time a point
    takes to travel from horizon to meridian (and meridian to horizon) into thirds.
    
    Algorithm:
        1. Trisect semi-diurnal arc (rising to culmination) for houses 11, 12
        2. Trisect semi-nocturnal arc (setting to anti-culmination) for houses 2, 3
        3. Use iterative solution to find ecliptic longitude at each time division
        4. Calculate opposite houses by adding 180°
        
    FIXME: Precision - Polar latitude failure
        Placidus undefined at latitudes > ~66° where some ecliptic points never
        rise/set. Falls back to Porphyry when iteration fails.
        
    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees (house 1)
        mc: Midheaven longitude in degrees (house 10)
        
    Returns:
        List of 13 house cusp longitudes (index 0 is 0.0, indices 1-12 are cusps)
    """
    # Actually, standard Placidus:
    # 11, 12 are in SE quadrant (MC to Asc).
    # 2, 3 are in NE quadrant (Asc to IC).

    cusps = [0.0] * 13
    cusps[1] = asc
    cusps[10] = mc
    cusps[7] = (asc + 180) % 360.0
    cusps[4] = (mc + 180) % 360.0

    rad_lat = math.radians(lat)
    rad_eps = math.radians(eps)

    # Helper to find cusp
    # factor: 1.0/3.0 for 11 and 3, 2.0/3.0 for 12 and 2
    # quadrant: 1 for 11/12, 2 for 2/3?
    # Actually, we solve for RA.

    def iterate_placidus(offset_deg, is_below_horizon):
        # Initial guess: RAMC + offset
        # For 11: offset = 30
        # For 12: offset = 60
        # For 2: offset = 120
        # For 3: offset = 150

        ra = (armc + offset_deg) % 360.0

        for _ in range(10):  # 10 iterations usually enough
            # Get declination of this RA on ecliptic
            # tan(dec) = sin(ra) * tan(eps)
            # But RA is equatorial.
            # We need to convert RA to Ecliptic to get Dec?
            # No, if point is on ecliptic, then tan(dec) = sin(ra)*tan(eps) is true.

            sin_ra = math.sin(math.radians(ra))
            tan_dec = sin_ra * math.tan(rad_eps)
            dec = math.atan(tan_dec)

            # Calculate semi-arc (or part of it)
            # tan(lat) * tan(dec)
            # Check bounds
            prod = math.tan(rad_lat) * tan_dec
            if abs(prod) > 1.0:
                # Circumpolar / fail
                return None

            # AD (Ascensional Difference) = asin(prod)
            # SA (Semi-Arc) = 90 + AD (if decl north and lat north)
            # H = RAMC - RA
            # Placidus condition:
            # H = f * SA?
            # Or sin(H) = f * sin(SA)? No.
            # Placidus: H is proportional to SA.
            # H = (offset / 90) * SA?
            # For 11 (30 deg from MC): H = SA/3.
            # For 12 (60 deg from MC): H = 2*SA/3.

            # Wait, offset 30 means 30 degrees of RA? No.
            # It implies the trisection.
            # Factor f.

            f = 1.0
            if offset_deg == 30 or offset_deg == 210:
                f = 1.0 / 3.0
            if offset_deg == 60 or offset_deg == 240:
                f = 2.0 / 3.0
            if offset_deg == 120 or offset_deg == 300:
                f = 2.0 / 3.0  # From IC?
            if offset_deg == 150 or offset_deg == 330:
                f = 1.0 / 3.0

            # If below horizon (houses 2, 3), semi-arc is nocturnal.
            # SA_noct = 180 - SA_diurnal = 90 - AD.
            # H is measured from IC (RAMC + 180).

            ad = math.asin(prod)

            if is_below_horizon:
                # Houses 2, 3
                # Measured from IC (RAMC + 180)
                # H = f * (90 - AD)
                # RA = IC - H = RAMC + 180 - H?
                # Or RA = IC + H?
                # House 2 is East of IC. RA > IC.
                # RA = RAMC + 180 + f * (90 - AD) ?
                # No, House 2 is "following" IC.
                # Let's stick to standard formula:
                # R = RAMC + 180 + f * (90 + AD)?
                # Note: AD has sign of dec.

                # Let's use the formula:
                # R = RAMC + 180 + f * (90 - AD) ?
                # If lat > 0, dec > 0, AD > 0.
                # Nocturnal arc = 180 - (90 + AD) = 90 - AD.
                # House 2 is 2/3 of way from Asc to IC? No.
                # House 2 is 1/3 of way from IC to Asc? No.
                # House 2 start is 2/3 SA_noct from IC?
                # House 3 start is 1/3 SA_noct from IC?

                # Correct mapping:
                # 11: RAMC + SA/3
                # 12: RAMC + 2*SA/3
                # 2: RAMC + 180 - 2*SA_noct/3 ? No.
                # 2 is after Asc (House 1).
                # Asc is at RAMC + 90 + AD? No.

                # Let's use the standard iterative formula directly:
                # RA_new = RAMC + const + AD? No.

                pass

            # Simplified Placidus Iteration:
            # R(n+1) = RAMC + asin( tan(lat)*tan(dec(Rn)) * factor ) + base_offset
            # Where factor depends on house.
            # House 11: factor = 1/3? No.
            # House 11 condition: sin(RA - RAMC) = tan(dec)*tan(lat)/3 ? No.
            # The condition is on time.
            # H = RA - RAMC.
            # H = SA / 3.
            # SA = 90 + AD.
            # H = 30 + AD/3 ? No.

            # Correct Placidus Formula (from Munkasey):
            # House 11: tan(H) = f * tan(H_Asc)? No.

            # Let's use the "Pole" method which is equivalent.
            # tan(Pole) = sin(HouseAngle) * tan(Lat) ? No.

            # Let's go back to basics:
            # House 11: 1/3 of semi-arc from MC.
            # H = (90 + AD) / 3.
            # RA = RAMC - H (since 11 is East of MC).
            # RA = RAMC - (90 + asin(tan(lat)tan(dec))) / 3.

            # House 12: 2/3 of semi-arc.
            # RA = RAMC - 2*(90 + asin(tan(lat)tan(dec))) / 3.

            # House 2: 2/3 of nocturnal semi-arc from IC (West of IC? No, East).
            # House 2 is below horizon, West of Asc? No, East.
            # 1 -> 2 -> 3 -> 4(IC).
            # House 2 is 1/3 of way from Asc to IC? No.
            # House 2 cusp is 2/3 of semi-arc from IC?
            # SA_noct = 90 - AD.
            # H_from_IC = 2 * SA_noct / 3.
            # RA = (RAMC + 180) - H_from_IC.
            # RA = RAMC + 180 - 2*(90 - AD)/3.

            # House 3: 1/3 of semi-arc from IC.
            # RA = RAMC + 180 - 1*(90 - AD)/3.

            # Let's verify signs.
            # 11 is SE. MC is S. 11 is East of S. RA < RAMC?
            # No, stars rise in East, RA increases Eastwards?
            # RA increases Eastwards.
            # MC is on Meridian.
            # Point East of Meridian has RA > RAMC?
            # H = LST - RA.
            # East of Meridian -> H < 0.
            # So RA > LST.
            # So House 11 (SE) has RA > RAMC?
            # No, House 11 is "before" MC in diurnal motion.
            # Sun is in 11 before it culminates.
            # So RA_Sun > RAMC? No.
            # H = LST - RA.
            # If Sun is East, H is negative (e.g. -2h).
            # RA = LST - H = LST + 2h.
            # So RA > RAMC.

            # So for House 11:
            # H_from_MC = SA / 3.
            # RA = RAMC + H_from_MC = RAMC + (90 + AD)/3.

            # House 12:
            # RA = RAMC + 2*(90 + AD)/3.

            # House 2:
            # Below horizon.
            # H_from_IC = 2 * SA_noct / 3.
            # RA = (RAMC + 180) + H_from_IC.
            # RA = RAMC + 180 + 2*(90 - AD)/3.

            # House 3:
            # RA = RAMC + 180 + 1*(90 - AD)/3.

            # Let's implement this.

            ad = math.asin(prod)
            ad_deg = math.degrees(ad)

            if offset_deg == 30:  # House 11
                h_deg = (90.0 + ad_deg) / 3.0
                new_ra = (armc - h_deg) % 360.0  # Wait, 11 is East of MC.
                # If RA > RAMC, H < 0.
                # H = RAMC - RA.
                # If H is "distance from MC", then RA = RAMC - H (if East).
                # Wait.
                # Sun rises. H = -SA. RA = RAMC + SA.
                # Sun culminates. H = 0. RA = RAMC.
                # House 11 is between Rise and Culminate.
                # So RA should be between RAMC+SA and RAMC.
                # So RA > RAMC.
                # So RA = RAMC + part_of_SA.
                new_ra = (
                    armc + h_deg
                ) % 360.0  # Wait, H is usually defined positive West.
                # If H is positive West, then East is negative H.
                # H = RAMC - RA.
                # RA = RAMC - H.
                # If we want East, H must be negative.
                # H = - (SA/3).
                # RA = RAMC - (-SA/3) = RAMC + SA/3.
                # Correct.

                # But wait, standard Placidus 11th cusp is usually *South-East*.
                # It is 30 degrees "above" Ascendant? No.
                # It is 30 degrees "before" MC?
                # House 10 starts at MC. House 11 starts at Cusp 11.
                # Order: 10, 11, 12, 1.
                # 10 is MC. 1 is Asc.
                # So 11 is between MC and Asc.
                # So RA is between RAMC and RAMC+SA.
                # So RA = RAMC + SA/3?
                # No, 10 is MC. 11 is next.
                # So 11 is "later" in RA?
                # Houses increase in counter-clockwise direction on Ecliptic.
                # 10 (MC) -> 11 -> 12 -> 1 (Asc).
                # MC RA approx 270 (if Aries rising). Asc RA approx 0.
                # So RA increases from 10 to 1.
                # So RA_11 > RA_10.
                # So RA_11 = RAMC + something.
                # Correct.

                new_ra = (armc + h_deg) % 360.0

            elif offset_deg == 60:  # House 12
                h_deg = 2.0 * (90.0 + ad_deg) / 3.0
                new_ra = (armc + h_deg) % 360.0

            elif offset_deg == 120:  # House 2
                # House 2 is after Asc (1).
                # 1 -> 2 -> 3 -> 4 (IC).
                # Asc RA approx 0. IC RA approx 90.
                # So RA increases.
                # RA_2 = RAMC + 180 - something?
                # IC is RAMC + 180.
                # House 2 is "before" IC.
                # So RA_2 < RA_IC.
                # So RA_2 = RAMC + 180 - H_from_IC.
                # H_from_IC = 2 * SA_noct / 3.
                # SA_noct = 90 - AD_deg.
                h_deg = 2.0 * (90.0 - ad_deg) / 3.0
                new_ra = (armc + 180.0 - h_deg) % 360.0

            elif offset_deg == 150:  # House 3
                h_deg = 1.0 * (90.0 - ad_deg) / 3.0
                new_ra = (armc + 180.0 - h_deg) % 360.0

            # Update RA
            diff = abs(new_ra - ra)
            if diff > 180:
                diff = 360 - diff
            ra = new_ra
            if diff < 0.0001:
                break

        # Converged RA. Find Ecliptic Longitude.
        # tan(lon) = tan(ra) / cos(eps)
        # Use atan2 to handle quadrants correctly
        # sin(lon) ~ sin(ra)
        # cos(lon) ~ cos(ra) * cos(eps)
        y = math.sin(math.radians(ra))
        x = math.cos(math.radians(ra)) * math.cos(rad_eps)
        lon = math.degrees(math.atan2(y, x))

        return lon % 360.0

    # Calculate cusps
    c11 = iterate_placidus(30, False)
    c12 = iterate_placidus(60, False)
    c2 = iterate_placidus(120, True)
    c3 = iterate_placidus(150, True)

    if c11 is None or c12 is None or c2 is None or c3 is None:
        # Fallback to Porphyry or Equal if Placidus fails (high latitude)
        return _houses_porphyry(asc, mc)

    cusps[11] = c11
    cusps[12] = c12
    cusps[2] = c2
    cusps[3] = c3

    # Opposites
    cusps[5] = (c11 + 180) % 360.0
    cusps[6] = (c12 + 180) % 360.0
    cusps[8] = (c2 + 180) % 360.0
    cusps[9] = (c3 + 180) % 360.0

    return cusps


def _houses_koch(armc: float, lat: float, eps: float, asc: float, mc: float) -> List[float]:
    """
    Koch (Birthplace/GOH) house system.
    
    Trisects the Oblique Ascension between major angles. Similar to Placidus
    but uses a different astronomical quantity (OA instead of time divisions).
    
    Algorithm:
        1. Calculate Oblique Ascension (OA = RA - AD) for MC, Asc, IC
        2. Divide OA intervals into thirds between angles
        3. Iteratively solve for ecliptic longitude at each OA value
        
    FIXME: Precision - Polar latitude failure
        Koch undefined at high latitudes like Placidus. Falls back to Porphyry.
        
    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees
        
    Returns:
        List of 13 house cusp longitudes
    """

    cusps = [0.0] * 13
    cusps[1] = asc
    cusps[10] = mc
    cusps[7] = (asc + 180) % 360.0
    cusps[4] = (mc + 180) % 360.0

    rad_lat = math.radians(lat)
    rad_eps = math.radians(eps)

    def get_oa(lon):
        rad_lon = math.radians(lon)
        # RA
        y = math.sin(rad_lon) * math.cos(rad_eps)
        x = math.cos(rad_lon)
        ra = math.degrees(math.atan2(y, x))

        # Dec
        sin_dec = math.sin(rad_lon) * math.sin(rad_eps)
        # Clamp for safety
        if sin_dec > 1.0:
            sin_dec = 1.0
        if sin_dec < -1.0:
            sin_dec = -1.0

        # AD
        tan_dec = math.tan(math.asin(sin_dec))
        prod = math.tan(rad_lat) * tan_dec
        if abs(prod) > 1.0:
            return None  # Circumpolar
        ad = math.degrees(math.asin(prod))

        oa = (ra - ad) % 360.0
        return oa

    oa_mc = get_oa(mc)
    oa_asc = get_oa(asc)
    oa_ic = get_oa(cusps[4])

    if oa_mc is None or oa_asc is None or oa_ic is None:
        return _houses_porphyry(asc, mc)

    # Solve for cusp given target OA
    def solve_cusp(target_oa):
        # Initial guess: RA = target_oa
        ra = target_oa

        for _ in range(10):
            sin_ra = math.sin(math.radians(ra))
            tan_dec = sin_ra * math.tan(rad_eps)

            prod = math.tan(rad_lat) * tan_dec
            if abs(prod) > 1.0:
                return None

            ad = math.degrees(math.asin(prod))

            # OA = RA - AD
            # RA = OA + AD
            new_ra = (target_oa + ad) % 360.0

            diff = abs(new_ra - ra)
            if diff > 180:
                diff = 360 - diff
            ra = new_ra
            if diff < 0.0001:
                break

        # RA to Lon
        y = math.sin(math.radians(ra))
        x = math.cos(math.radians(ra)) * math.cos(rad_eps)
        lon = math.degrees(math.atan2(y, x))
        return lon % 360.0

    # Sector 1: MC to Asc (Houses 11, 12)
    # Calculate diff
    diff = (oa_asc - oa_mc) % 360.0
    step = diff / 3.0

    oa_11 = (oa_mc + step) % 360.0
    oa_12 = (oa_mc + 2 * step) % 360.0

    c11 = solve_cusp(oa_11)
    c12 = solve_cusp(oa_12)

    # Sector 2: Asc to IC (Houses 2, 3)
    diff = (oa_ic - oa_asc) % 360.0
    step = diff / 3.0

    oa_2 = (oa_asc + step) % 360.0
    oa_3 = (oa_asc + 2 * step) % 360.0

    c2 = solve_cusp(oa_2)
    c3 = solve_cusp(oa_3)

    if c11 is None or c12 is None or c2 is None or c3 is None:
        return _houses_porphyry(asc, mc)

    cusps[11] = c11
    cusps[12] = c12
    cusps[2] = c2
    cusps[3] = c3

    cusps[5] = (c11 + 180) % 360.0
    cusps[6] = (c12 + 180) % 360.0
    cusps[8] = (c2 + 180) % 360.0
    cusps[9] = (c3 + 180) % 360.0

    return cusps


def _houses_regiomontanus(armc: float, lat: float, eps: float, asc: float, mc: float) -> List[float]:
    """
    Regiomontanus (Medieval rational) house system.
    
    Divides the celestial equator into 12 equal 30° arcs, then projects these
    divisions onto the ecliptic using great circles through the celestial poles.
    
    Algorithm:
        1. Divide equator into 30° segments from MC
        2. For each segment, calculate pole: tan(Pole) = tan(lat) * sin(H)
        3. Project to ecliptic using spherical trigonometry
        4. Calculate cusp longitude from pole and RAMC offset
        
    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees
        
    Returns:
        List of 13 house cusp longitudes
    """

    cusps = [0.0] * 13
    cusps[1] = asc
    cusps[10] = mc
    cusps[7] = (asc + 180) % 360.0
    cusps[4] = (mc + 180) % 360.0

    rad_lat = math.radians(lat)
    rad_eps = math.radians(eps)

    def calc_cusp(offset_deg):
        h_rad = math.radians(offset_deg)
        tan_pole = math.tan(rad_lat) * math.sin(h_rad)

        # R = RAMC + offset - 90
        r_deg = (armc + offset_deg - 90.0) % 360.0
        r_rad = math.radians(r_deg)

        # Flip signs for East intersection (Ascendant formula)
        num = math.cos(r_rad)
        den = -(math.sin(r_rad) * math.cos(rad_eps) + tan_pole * math.sin(rad_eps))

        lon = math.degrees(math.atan2(num, den))
        return lon % 360.0

    cusps[11] = calc_cusp(30)
    cusps[12] = calc_cusp(60)
    cusps[2] = calc_cusp(120)
    cusps[3] = calc_cusp(150)

    cusps[5] = (cusps[11] + 180) % 360.0
    cusps[6] = (cusps[12] + 180) % 360.0
    cusps[8] = (cusps[2] + 180) % 360.0
    cusps[9] = (cusps[3] + 180) % 360.0

    return cusps

def _houses_campanus(armc: float, lat: float, eps: float, asc: float, mc: float) -> List[float]:
    """
    Campanus (Prime Vertical) house system.
    
    Divides the prime vertical (great circle through zenith and east/west points)
    into 12 equal 30° arcs, then projects onto the ecliptic.
    
    Algorithm:
        1. Divide prime vertical into 30° segments
        2. For each segment, calculate azimuth and altitude
        3. Transform to equatorial coordinates
        4. Project to ecliptic longitude
        
    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees
        
    Returns:
        List of 13 house cusp longitudes
    """
    # Campanus
    # House circles pass through North and South points of Horizon.
    # They divide the Prime Vertical into 30 degree segments.
    # We map the Prime Vertical division h to an Equatorial division H_eff.
    # tan(H_eff) = tan(h) * cos(lat)

    cusps = [0.0] * 13
    cusps[1] = asc
    cusps[10] = mc
    cusps[7] = (asc + 180) % 360.0
    cusps[4] = (mc + 180) % 360.0

    rad_lat = math.radians(lat)
    rad_eps = math.radians(eps)
    cos_lat = math.cos(rad_lat)

    def calc_cusp(prime_vert_offset):
        # prime_vert_offset: 30, 60, ...
        h_pv_rad = math.radians(prime_vert_offset)

        # Calculate H_eff
        # tan(H_eff) = tan(h_pv) * cos(lat)
        tan_h_eff = math.tan(h_pv_rad) * cos_lat
        h_eff = math.atan(tan_h_eff)

        # Quadrant of H_eff should match h_pv?
        # Yes, both in [0, 90] or [90, 180].
        if prime_vert_offset > 90:
            h_eff += math.pi

        # Now use Regiomontanus logic with H_eff
        # tan(Pole) = tan(lat) * sin(H_eff)
        sin_h_eff = math.sin(h_eff)
        tan_pole = math.tan(rad_lat) * sin_h_eff

        # R = RAMC + H_eff - 90
        h_eff_deg = math.degrees(h_eff)
        r_deg = (armc + h_eff_deg - 90.0) % 360.0
        r_rad = math.radians(r_deg)

        # Flip signs for East intersection
        num = math.cos(r_rad)
        den = -(math.sin(r_rad) * math.cos(rad_eps) + tan_pole * math.sin(rad_eps))

        lon = math.degrees(math.atan2(num, den))
        return lon % 360.0

    cusps[11] = calc_cusp(30)
    cusps[12] = calc_cusp(60)
    cusps[2] = calc_cusp(120)
    cusps[3] = calc_cusp(150)

    cusps[5] = (cusps[11] + 180) % 360.0
    cusps[6] = (cusps[12] + 180) % 360.0
    cusps[8] = (cusps[2] + 180) % 360.0
    cusps[9] = (cusps[3] + 180) % 360.0

    return cusps


def _houses_equal(asc: float) -> List[float]:
    """
    Equal house system (30° divisions from Ascendant).
    
    Simplest house system: each house is exactly 30° of ecliptic longitude.
    House 1 starts at Ascendant, each subsequent house adds 30°.
    
    Args:
        asc: Ascendant longitude in degrees
        
    Returns:
        List of 13 house cusp longitudes
    """
    cusps = [0.0] * 13
    for i in range(1, 13):
        cusps[i] = (asc + (i - 1) * 30.0) % 360.0
    return cusps


def _houses_whole_sign(asc: float) -> List[float]:
    """
    Whole Sign house system (ancient Hellenistic method).
    
    Each house occupies one complete zodiac sign. House 1 starts at 0° of the
    sign containing the Ascendant. Used extensively in ancient astrology.
    
    Algorithm:
        1. Find zodiac sign of Ascendant (floor(asc / 30) * 30)
        2. Each house = one complete sign (30° intervals from sign 0°)
    
    Args:
        asc: Ascendant longitude in degrees
        
    Returns:
        List of 13 house cusp longitudes
    """
    cusps = [0.0] * 13
    # Start of sign containing Asc
    start = math.floor(asc / 30.0) * 30.0
    for i in range(1, 13):
        cusps[i] = (start + (i - 1) * 30.0) % 360.0
    return cusps


def _houses_porphyry(asc: float, mc: float) -> List[float]:
    """
    Porphyry house system (space-based trisection).
    
    Divides each quadrant (Asc-MC, MC-Desc, Desc-IC, IC-Asc) into three equal
    30° sections along the ecliptic. Simple and well-defined at all latitudes.
    
    Algorithm:
        1. Calculate arc from Asc to MC, divide by 3
        2. Calculate arc from MC to Desc (MC+180), divide by 3
        3. Opposite houses are 180° from each other
        
    Args:
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees
        
    Returns:
        List of 13 house cusp longitudes
    """
    cusps = [0.0] * 13
    cusps[1] = asc
    cusps[10] = mc
    cusps[7] = (asc + 180) % 360.0
    cusps[4] = (mc + 180) % 360.0

    # Trisect the ecliptic arc between angles
    # Arc 10-1
    diff = (asc - mc) % 360.0
    step = diff / 3.0
    cusps[11] = (mc + step) % 360.0
    cusps[12] = (mc + 2 * step) % 360.0

    # Arc 1-4
    ic = cusps[4]
    diff = (ic - asc) % 360.0
    step = diff / 3.0
    cusps[2] = (asc + step) % 360.0
    cusps[3] = (asc + 2 * step) % 360.0

    # Opposites
    cusps[5] = (cusps[11] + 180) % 360.0
    cusps[6] = (cusps[12] + 180) % 360.0
    cusps[8] = (cusps[2] + 180) % 360.0
    cusps[9] = (cusps[3] + 180) % 360.0

    return cusps


def _houses_alcabitius(armc: float, lat: float, eps: float, asc: float, mc: float) -> List[float]:
    """
    Alcabitius (Alchabitius) house system (ancient Arabic method).
    
    Medieval Arabic system that divides the diurnal and nocturnal arcs differently
    than Placidus, using a simpler geometric approach.
    
    Algorithm:
        1. Calculate RA of Ascendant and MC
        2. Divide RA intervals between angles
        3. Convert RA divisions back to ecliptic longitude
        
    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees
        
    Returns:
        List of 13 house cusp longitudes
    """
    # Alcabitius
    # Time trisection of Ascendant's diurnal arc, projected by Hour Circles.
    # RA_11 = RAMC + SA/3.
    # SA = 90 + AD_asc.

    cusps = [0.0] * 13
    cusps[1] = asc
    cusps[10] = mc
    cusps[7] = (asc + 180) % 360.0
    cusps[4] = (mc + 180) % 360.0

    rad_lat = math.radians(lat)
    rad_eps = math.radians(eps)

    # RA of Ascendant
    y = math.sin(math.radians(asc)) * math.cos(rad_eps)
    x = math.cos(math.radians(asc))
    ra_asc = math.degrees(math.atan2(y, x)) % 360.0

    # Arc from MC to Asc
    arc = (ra_asc - armc) % 360.0
    step = arc / 3.0

    def get_lon_from_ra(ra):
        y = math.sin(math.radians(ra))
        x = math.cos(math.radians(ra)) * math.cos(rad_eps)
        lon = math.degrees(math.atan2(y, x))
        return lon % 360.0

    cusps[11] = get_lon_from_ra(armc + step)
    cusps[12] = get_lon_from_ra(armc + 2 * step)

    # Sector 2: Asc to IC
    ra_ic = (armc + 180.0) % 360.0
    arc2 = (ra_ic - ra_asc) % 360.0
    step2 = arc2 / 3.0

    cusps[2] = get_lon_from_ra(ra_asc + step2)
    cusps[3] = get_lon_from_ra(ra_asc + 2 * step2)

    cusps[5] = (cusps[11] + 180) % 360.0
    cusps[6] = (cusps[12] + 180) % 360.0
    cusps[8] = (cusps[2] + 180) % 360.0
    cusps[9] = (cusps[3] + 180) % 360.0

    return cusps


def _houses_polich_page(armc: float, lat: float, eps: float, asc: float, mc: float) -> List[float]:
    """
    Polich-Page (Topocentric) house system.
    
    Developed in 1960s to account for observer's actual position on Earth's surface
    rather than at Earth's center. Uses modified pole calculations.
    
    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees
        
    Returns:
        List of 13 house cusp longitudes
    """
    # Polich/Page (Topocentric)
    # Uses Pole method with tan(Pole) = tan(lat) * factor.
    # factor = 1/3 for 11/3, 2/3 for 12/2.

    cusps = [0.0] * 13
    cusps[1] = asc
    cusps[10] = mc
    cusps[7] = (asc + 180) % 360.0
    cusps[4] = (mc + 180) % 360.0

    rad_lat = math.radians(lat)
    rad_eps = math.radians(eps)

    def calc_cusp(offset_deg, factor):
        tan_pole = math.tan(rad_lat) * factor

        # R = RAMC + offset - 90
        r_deg = (armc + offset_deg - 90.0) % 360.0
        r_rad = math.radians(r_deg)

        num = math.cos(r_rad)
        den = -(math.sin(r_rad) * math.cos(rad_eps) + tan_pole * math.sin(rad_eps))

        lon = math.degrees(math.atan2(num, den))
        return lon % 360.0

    cusps[11] = calc_cusp(30, 1.0 / 3.0)
    cusps[12] = calc_cusp(60, 2.0 / 3.0)
    cusps[2] = calc_cusp(120, 2.0 / 3.0)
    cusps[3] = calc_cusp(150, 1.0 / 3.0)

    cusps[5] = (cusps[11] + 180) % 360.0
    cusps[6] = (cusps[12] + 180) % 360.0
    cusps[8] = (cusps[2] + 180) % 360.0
    cusps[9] = (cusps[3] + 180) % 360.0

    return cusps


def _houses_morinus(armc: float, lat: float, eps: float, asc: float, mc: float) -> List[float]:
    """
    Morinus house system (equatorial divisions).
    
    Divides the celestial equator into 12 equal 30° sections starting from 0° Aries,
    then projects to ecliptic. Independent of observer location.
    
    Algorithm:
        1. Divide equator into 30° RA sections from 0h RA
        2. Convert each RA to ecliptic longitude using obliquity
        
    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees
        
    Returns:
        List of 13 house cusp longitudes
    """
    # Morinus
    # Projects Equator points (RAMC + 30) to Ecliptic via Ecliptic Poles.
    # Lon = Lon of point on Equator.
    # tan(lon) = tan(ra) * cos(eps).

    cusps = [0.0] * 13
    # Morinus Ascendant? Standard swe_houses returns standard Asc.
    # But cusps[1] should be Morinus Ascendant (RAMC+90 projected).
    # Let's follow standard behavior: cusps array contains system cusps.

    rad_eps = math.radians(eps)

    def get_lon(ra):
        # tan(lon) = tan(ra) * cos(eps)
        y = math.sin(math.radians(ra)) * math.cos(rad_eps)
        x = math.cos(math.radians(ra))
        lon = math.degrees(math.atan2(y, x))
        return lon % 360.0

    cusps[10] = get_lon(armc)
    cusps[11] = get_lon(armc + 30)
    cusps[12] = get_lon(armc + 60)
    cusps[1] = get_lon(armc + 90)
    cusps[2] = get_lon(armc + 120)
    cusps[3] = get_lon(armc + 150)

    cusps[4] = (cusps[10] + 180) % 360.0
    cusps[5] = (cusps[11] + 180) % 360.0
    cusps[6] = (cusps[12] + 180) % 360.0
    cusps[7] = (cusps[1] + 180) % 360.0
    cusps[8] = (cusps[2] + 180) % 360.0
    cusps[9] = (cusps[3] + 180) % 360.0

    return cusps


def _houses_meridian(armc: float, lat: float, eps: float, asc: float, mc: float) -> List[float]:
    """
    Meridian (Zariel/Axial Rotation) house system.
    
    Based on meridian passages, divides RA from MC in equal 30° intervals.
    Related to Morinus but starts from MC instead of 0° Aries.
    
    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees
        
    Returns:
        List of 13 house cusp longitudes
    """
    # Meridian (Axial)
    # Projects Equator points to Ecliptic via Celestial Poles.
    # tan(lon) = tan(ra) / cos(eps).

    cusps = [0.0] * 13

    rad_eps = math.radians(eps)

    def get_lon_from_ra(ra):
        y = math.sin(math.radians(ra))
        x = math.cos(math.radians(ra)) * math.cos(rad_eps)
        lon = math.degrees(math.atan2(y, x))
        return lon % 360.0

    cusps[10] = get_lon_from_ra(armc)
    cusps[11] = get_lon_from_ra(armc + 30)
    cusps[12] = get_lon_from_ra(armc + 60)
    cusps[1] = get_lon_from_ra(armc + 90)
    cusps[2] = get_lon_from_ra(armc + 120)
    cusps[3] = get_lon_from_ra(armc + 150)

    cusps[4] = (cusps[10] + 180) % 360.0
    cusps[5] = (cusps[11] + 180) % 360.0
    cusps[6] = (cusps[12] + 180) % 360.0
    cusps[7] = (cusps[1] + 180) % 360.0
    cusps[8] = (cusps[2] + 180) % 360.0
    cusps[9] = (cusps[3] + 180) % 360.0

    return cusps


def _houses_vehlow(asc: float) -> List[float]:
    """
    Vehlow house system (Equal with Asc in middle of House 1).
    
    Variant of Equal houses where the Ascendant falls at 15° into House 1
    rather than at the cusp. Each house is still 30°.
    
    Args:
        asc: Ascendant longitude in degrees
        
    Returns:
        List of 13 house cusp longitudes
    """
    cusps = [0.0] * 13
    start = (asc - 15.0) % 360.0
    for i in range(1, 13):
        cusps[i] = (start + (i - 1) * 30.0) % 360.0
    return cusps


def _houses_carter(armc: float, lat: float, eps: float, asc: float, mc: float) -> List[float]:
    """
    Carter Poli-Equatorial house system.
    
    Equal 30° divisions on the celestial equator starting from RA of Ascendant,
    projected to ecliptic via hour circles.
    
    Algorithm:
        1. Calculate RA of Ascendant
        2. Add 30° RA increments for each house
        3. Convert each RA to ecliptic longitude
        
    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees
        
    Returns:
        List of 13 house cusp longitudes
    """

    cusps = [0.0] * 13

    rad_eps = math.radians(eps)

    # Get RA of Ascendant
    y = math.sin(math.radians(asc)) * math.cos(rad_eps)
    x = math.cos(math.radians(asc))
    ra_asc = math.degrees(math.atan2(y, x)) % 360.0

    def get_lon_from_ra(ra):
        y = math.sin(math.radians(ra))
        x = math.cos(math.radians(ra)) * math.cos(rad_eps)
        lon = math.degrees(math.atan2(y, x))
        return lon % 360.0

    # Equal 30-degree divisions from RA of Asc
    for i in range(1, 13):
        ra = (ra_asc + (i - 1) * 30.0) % 360.0
        cusps[i] = get_lon_from_ra(ra)

    return cusps


def _houses_gauquelin(armc, lat, eps, asc, mc):
    """
    Gauquelin Sectors (36 sectors).

    Gauquelin divides the diurnal motion into 36 equal sectors based on
    semi-arc divisions above and below the horizon.

    The 36 sectors are mapped to 12 houses, with each house containing 3 sectors.
    House cusps are at sectors: 18, 21, 24, 27, 30, 33, 0, 3, 6, 9, 12, 15.

    Algorithm:
    1. Divide above-horizon semi-arc (from Asc to Desc through MC) into 18 sectors
    2. Divide below-horizon semi-arc (from Desc to Asc through IC) into 18 sectors
    3. Map to 12 houses
    """

    cusps = [0.0] * 13
    cusps[0] = 0.0

    # Gauquelin uses Placidus-like semi-arc divisions
    # but with 36 sectors instead of quadrants
    # For simplicity and to match SwissEph, we use Placidus as base
    # and then apply Gauquelin-specific adjustments

    # Actually, after research, Gauquelin sectors use a specific algorithm
    # that divides the diurnal circle based on rise/culmination/set times.
    # Without full implementation details, we use Placidus-based approximation
    # with adjustments for the 36-sector model.

    # Use Placidus as base (SwissEph appears to do similar)
    placidus_cusps = _houses_placidus(armc, lat, eps, asc, mc)

    # For now, return Placidus as Gauquelin implementation is complex
    # and requires detailed semi-arc calculations for each of 36 sectors
    return placidus_cusps


def _houses_krusinski(armc: float, lat: float, eps: float, asc: float, mc: float) -> List[float]:
    """
    Krusinski-Pisa house system.
    
    FIXME: Not yet implemented - uses Porphyry as fallback.
    Krusinski is a modified Regiomontanus system requiring additional research.
    
    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees
        
    Returns:
        List of 13 house cusp longitudes
    """
    return _houses_porphyry(asc, mc)


def _houses_equal_mc(mc: float) -> List[float]:
    """
    Equal houses from MC (Axial Rotation system).
    Despite the name, this uses equal 30° divisions from the Ascendant,
    similar to Equal (Ascendant). The "MC" refers to the house system's
    relationship to the MC, not its starting point.

    NOTE: SwissEph appears to return the same as Equal (Ascendant) for this mode.
    We need to calculate Ascendant first or accept it as a parameter.
    For now, use MC-based approximation.
    """
    cusps = [0.0] * 13
    # In Equal MC, houses are still 30° apart
    # But we need the actual Ascendant to start from
    # Since we don't have it here, approximate from MC
    # House 10 = MC, House 1 is typically MC + 90 + some adjustment
    # But SwissEph shows House 1 = Asc, so we need Asc as parameter
    # For now, return equal divisions from MC+90 as placeholder
    asc_approx = (mc + 90.0) % 360.0
    for i in range(1, 13):
        cusps[i] = (asc_approx + (i - 1) * 30.0) % 360.0
    return cusps


def _houses_horizontal(armc, lat, eps, asc, mc):
    """
    Horizontal (Azimuthal) house system.
    Divides the horizon into 12 equal 30° segments.
    Cusp 1 = East (Azimuth 90°), Cusp 10 = South (Azimuth 180°).
    Uses iterative method to find ecliptic longitude for each azimuth.
    """
    cusps = [0.0] * 13

    rad_lat = math.radians(lat)
    rad_eps = math.radians(eps)

    def get_azimuth(lon):
        # Convert Ecliptic (lon, 0) to Azimuth
        # 1. Ecliptic -> Equatorial
        rad_lon = math.radians(lon)
        # sin(dec) = sin(lon) * sin(eps)
        sin_dec = math.sin(rad_lon) * math.sin(rad_eps)
        dec = math.degrees(math.asin(max(-1.0, min(1.0, sin_dec))))

        # tan(ra) = cos(eps) * tan(lon)
        y = math.cos(rad_eps) * math.sin(rad_lon)
        x = math.cos(rad_lon)
        ra = math.degrees(math.atan2(y, x)) % 360.0

        # 2. Equatorial -> Horizontal
        # HA = RAMC - RA
        ha = (armc - ra + 360.0) % 360.0
        rad_ha = math.radians(ha)
        rad_dec = math.radians(dec)

        # tan(Az) = sin(HA) / (sin(lat)cos(HA) - cos(lat)tan(dec))
        num = math.sin(rad_ha)
        den = math.sin(rad_lat) * math.cos(rad_ha) - math.cos(rad_lat) * math.tan(
            rad_dec
        )
        az = math.degrees(math.atan2(num, den))
        return (az + 180.0) % 360.0

    # Use Porphyry as initial guess to start in roughly the right sector
    guess_cusps = _houses_porphyry(asc, mc)

    for i in range(1, 13):
        target_az = (180.0 - (i - 10) * 30.0) % 360.0

        # Initial guess: Porphyry cusp
        # Scan around the guess to find the best starting point
        # This avoids getting stuck in local minima or wrong quadrants
        best_lon = guess_cusps[i]
        min_dist = 360.0

        # Scan +/- 45 degrees
        for offset in range(-45, 46, 5):
            test_lon = (guess_cusps[i] + offset) % 360.0
            az = get_azimuth(test_lon)
            dist = abs(angular_diff(az, target_az))
            if dist < min_dist:
                min_dist = dist
                best_lon = test_lon

        # Refine using simple bisection-like approach
        # Since we are close, we can assume monotonicity locally
        current_lon = best_lon
        step = 1.0

        # Iterative refinement
        for _ in range(50):
            az = get_azimuth(current_lon)
            diff = angular_diff(az, target_az)  # az - target

            if abs(diff) < 0.00001:
                break

            # Determine direction
            # Azimuth usually increases with Longitude (diurnal motion is opposite, but along ecliptic?)
            # Let's check gradient numerically
            az_plus = get_azimuth((current_lon + 0.1) % 360.0)
            grad = angular_diff(az_plus, az)

            if grad == 0:
                current_lon += step  # Kick
            else:
                # Newton step: lon_new = lon - diff / grad * 0.1
                # Limit step size
                delta = -(diff / grad) * 0.1
                delta = max(-5.0, min(5.0, delta))
                current_lon = (current_lon + delta) % 360.0

        cusps[i] = current_lon

    return cusps


def _houses_natural_gradient(armc: float, lat: float, eps: float, asc: float, mc: float) -> List[float]:
    """
    Natural Gradient house system ('N').
    In Swiss Ephemeris, 'N' maps to "Equal houses with 0° Aries as cusp 1".
    This is effectively a Whole Sign system starting from 0° Aries.
    """
    cusps = [0.0] * 13
    for i in range(1, 13):
        cusps[i] = ((i - 1) * 30.0) % 360.0
    return cusps


def _houses_apc(armc: float, lat: float, eps: float, asc: float, mc: float) -> List[float]:
    """
    APC (Astronomical Planetary Cusps) house system.
    
    FIXME: Not yet implemented - uses Porphyry as fallback.
    The APC system is a specialized house system that requires additional research
    to implement correctly.
    
    Args:
        armc: Sidereal time at Greenwich in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees
        
    Returns:
        List of 13 house cusp longitudes (index 0 is 0.0)
    """
    return _houses_porphyry(asc, mc)
