"""
Crossing event calculations for libephemeris.

Finds exact times when the Sun or Moon cross specific ecliptic longitudes.
Uses Newton-Raphson iteration for sub-arcsecond precision.

Functions:
- swe_solcross_ut: Sun crossing events (e.g., ingresses, equinoxes)
- swe_mooncross_ut: Moon crossing events (for lunar mansion calculations)
- swe_cross_ut: Generic planet crossing

FIXME: Precision - Newton-Raphson convergence tolerance
Current tolerance: 1 arcsecond (1/3600°)
Iterations: 20-30 max
Typical convergence: 3-6 iterations for Sun, 5-10 for Moon

For sub-arcsecond precision, consider iterative refinement or
use Swiss Ephemeris swe_solcross_ut2 with higher tolerances.

Algorithm: Initial linear estimate + Newton-Raphson refinement
References: Meeus "Astronomical Algorithms" Ch. 5 (interpolation)
"""

from typing import Tuple
from .constants import SEFLG_SWIEPH, SEFLG_SPEED, SE_SUN, SE_MOON
from .planets import swe_calc_ut


def swe_solcross_ut(
    x2cross: float, jd_ut: float, flag: int = SEFLG_SWIEPH
) -> float:
    """
    Find when the Sun crosses a specific ecliptic longitude.
    
    Searches FORWARD in time for the next crossing after jd_ut.
    
    Args:
        x2cross: Target ecliptic longitude in degrees (0-360)
        jd_ut: Julian Day (UT) to start search from
        flag: Calculation flags (SEFLG_SWIEPH, etc.)
        
    Returns:
        float: Julian Day of crossing (UT)
        
    Raises:
        RuntimeError: If convergence fails or calculation error occurs
        
    Algorithm:
        1. Get current Sun position and velocity
        2. Linear estimate: dt = (target - current) / velocity
        3. Refine with Newton-Raphson: jd_new = jd + (target - actual) / velocity
        4. Converge to < 1 arcsecond (~2.8 seconds of time)
        
    Precision:
        Typically < 1 arcsecond (< 5 seconds of time for Sun)
        
    FIXME: Precision - Convergence tolerance 1 arcsec
        For higher precision, reduce tolerance in line 74
        Swiss Ephemeris achieves ~0.001 arcsec with tighter iteration
        
    Example:
        >>> # Find next Aries ingress (0°)
        >>> jd_ingress = swe_solcross_ut(0.0, jd_now)
        >>> # Find summer solstice (90°)
        >>> jd_solstice = swe_solcross_ut(90.0, jd_now)
    """
    x2cross = x2cross % 360.0

    try:
        pos, _ = swe_calc_ut(jd_ut, SE_SUN, flag | SEFLG_SPEED)
        lon_start = pos[0]
        speed = pos[3]  # degrees/day
    except Exception as e:
        raise RuntimeError(f"Failed to calculate Sun position: {e}")

    # Calculate angular distance to target (always forward)
    diff = (x2cross - lon_start) % 360.0

    # Handle retrograde motion
    if speed < 0 and diff > 0:
        diff -= 360.0

    # If already very close, look for next complete crossing
    if abs(diff) < 1e-5:
        if speed > 0:
            diff += 360.0
        else:
            diff -= 360.0

    # Initial time estimate (linear approximation)
    if speed == 0:
        speed = 0.9856  # Average Sun motion ~1°/day

    dt_guess = diff / speed
    jd_guess = jd_ut + dt_guess

    # FIXME: Precision - Newton-Raphson max iterations = 20
    # May fail for very slow bodies or near stationary points
    jd = jd_guess
    for iteration in range(20):
        try:
            pos, _ = swe_calc_ut(jd, SE_SUN, flag | SEFLG_SPEED)
            lon = pos[0]
            speed = pos[3]
        except Exception as e:
            raise RuntimeError(f"Failed to calculate Sun position during iteration: {e}")

        # Angular difference to target
        diff = (x2cross - lon) % 360.0
        if diff > 180:
            diff -= 360

        # FIXME: Precision - Convergence tolerance 1 arcsecond
        # Check convergence (< 1 arcsecond = 1/3600 degree)
        if abs(diff) < 1.0 / 3600.0:
            return jd

        # Newton-Raphson step
        if abs(speed) < 0.01:
            speed = 0.9856

        jd += diff / speed

        # Safety: prevent divergence
        if abs(jd - jd_guess) > 366:
            raise RuntimeError("Solar crossing search diverged")

    raise RuntimeError("Maximum iterations reached in solar crossing calculation")


def swe_mooncross_ut(
    x2cross: float, jd_ut: float, flag: int = SEFLG_SWIEPH
) -> float:
    """
    Find when the Moon crosses a specific ecliptic longitude.
    
    Searches FORWARD in time for the next crossing after jd_ut.
    
    Args:
        x2cross: Target ecliptic longitude in degrees (0-360)
        jd_ut: Julian Day (UT) to start search from
        flag: Calculation flags (SEFLG_SWIEPH, etc.)
        
    Returns:
        float: Julian Day of crossing (UT)
        
    Raises:
        RuntimeError: If convergence fails or calculation error occurs
        
    Note:
        Moon moves ~13° per day (27.3 day cycle).
        More variable speed than Sun due to orbit eccentricity.
        Uses 30 iterations vs Sun's 20 for added robustness.
        
    FIXME: Precision - Same 1 arcsec tolerance as Sun
        Moon's higher speed means ~0.1 seconds time precision
        
    Example:
        >>> # Find next new moon (Sun-Moon conjunction at same longitude)
        >>> sun_pos, _ = swe_calc_ut(jd_now, SE_SUN, SEFLG_SWIEPH)
        >>> jd_new_moon = swe_mooncross_ut(sun_pos[0], jd_now)
    """
    x2cross = x2cross % 360.0

    try:
        pos, _ = swe_calc_ut(jd_ut, SE_MOON, flag | SEFLG_SPEED)
        lon_start = pos[0]
        speed = pos[3]  # degrees/day
    except Exception as e:
        raise RuntimeError(f"Failed to calculate Moon position: {e}")

    # Calculate initial guess for NEXT crossing
    diff = (x2cross - lon_start) % 360.0

    if speed < 0 and diff > 0:
        diff -= 360.0

    if abs(diff) < 1e-5:
        if speed > 0:
            diff += 360.0
        else:
            diff -= 360.0

    # Initial time estimate
    if speed == 0:
        speed = 13.176  # Average Moon motion ~13.18°/day

    dt_guess = diff / speed
    jd_guess = jd_ut + dt_guess

    # FIXME: Precision - Moon uses 30 iterations (vs Sun's 20)
    # Higher iteration count for Moon's variable speed
    jd = jd_guess
    for iteration in range(30):
        try:
            pos, _ = swe_calc_ut(jd, SE_MOON, flag | SEFLG_SPEED)
            lon = pos[0]
            speed = pos[3]
        except Exception as e:
            raise RuntimeError(f"Failed to calculate Moon position during iteration: {e}")

        # Difference to target
        diff = (x2cross - lon) % 360.0
        if diff > 180:
            diff -= 360

        # Check convergence (< 1 arcsecond)
        if abs(diff) < 1.0 / 3600.0:
            return jd

        # Newton-Raphson step
        if abs(speed) < 0.1:
            speed = 13.176

        jd += diff / speed

        # Safety check
        if abs(jd - jd_guess) > 31:  # More than a month
            raise RuntimeError("Moon crossing search diverged")

    raise RuntimeError("Maximum iterations reached in moon crossing calculation")


def swe_cross_ut(
    planet_id: int, x2cross: float, jd_ut: float, flag: int = SEFLG_SWIEPH
) -> float:
    """
    Find when any planet crosses a specific ecliptic longitude.
    
    Generic version for all planets (Mercury, Venus, Mars, etc.).
    
    Args:
        planet_id: Planet ID (SE_MERCURY, SE_VENUS, etc.)
        x2cross: Target ecliptic longitude in degrees (0-360)
        jd_ut: Julian Day (UT) to start search from
        flag: Calculation flags
        
    Returns:
        float: Julian Day of crossing (UT)
        
    Raises:
        RuntimeError: If convergence fails or calculation error occurs
        
    Note:
        Uses adaptive iteration count based on typical planet speed.
        Slower planets (Jupiter, Saturn) may need more iterations.
        
    FIXME: Precision - Not optimized for retrograde stations
        Near stationary points (speed ~ 0), convergence is slower
        May need adaptive tolerance or damped Newton-Raphson
        
    Example:
        >>> # Mars ingress into Aries
        >>> jd_mars_aries = swe_cross_ut(SE_MARS, 0.0, jd_now)
    """
    x2cross = x2cross % 360.0

    try:
        pos, _ = swe_calc_ut(jd_ut, planet_id, flag | SEFLG_SPEED)
        lon_start = pos[0]
        speed = pos[3]
    except Exception as e:
        raise RuntimeError(f"Failed to calculate planet position: {e}")

    # Estimate typical speed if near zero
    typical_speeds = {
        2: 1.4,  # Mercury
        3: 1.2,  # Venus
        4: 0.5,  # Mars
        5: 0.08,  # Jupiter
        6: 0.03,  # Saturn
        7: 0.01,  # Uranus
        8: 0.006,  # Neptune
    }
    speed_default = typical_speeds.get(planet_id, 0.5)

    # Calculate initial guess
    diff = (x2cross - lon_start) % 360.0

    if speed < 0 and diff > 0:
        diff -= 360.0

    if abs(diff) < 1e-5:
        if speed > 0:
            diff += 360.0
        else:
            diff -= 360.0

    if abs(speed) < 0.001:
        speed = speed_default

    dt_guess = diff / speed
    jd_guess = jd_ut + dt_guess

    # Adaptive iteration count
    max_iter = 40 if abs(speed) < 0.1 else 20

    # FIXME: Precision - Adaptive iterations needed for slow planets
    # Outer planets may need 40+ iterations near stationary points
    jd = jd_guess
    for iteration in range(max_iter):
        try:
            pos, _ = swe_calc_ut(jd, planet_id, flag | SEFLG_SPEED)
            lon = pos[0]
            speed = pos[3]
        except Exception as e:
            raise RuntimeError(f"Failed to calculate planet position during iteration: {e}")

        diff = (x2cross - lon) % 360.0
        if diff > 180:
            diff -= 360

        if abs(diff) < 1.0 / 3600.0:
            return jd

        if abs(speed) < 0.001:
            speed = speed_default

        jd += diff / speed

        # Safety: longer range for slower planets
        max_range = 400 if abs(speed_default) < 0.1 else 366
        if abs(jd - jd_guess) > max_range:
            raise RuntimeError("Planet crossing search diverged")

    raise RuntimeError("Maximum iterations reached in planet crossing calculation")
