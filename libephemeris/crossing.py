# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Crossing event calculations for libephemeris.

Finds exact times when the Sun or Moon cross specific ecliptic longitudes.
Uses Newton-Raphson iteration for sub-milliarcsecond precision, with Brent's method
as a fallback near retrograde stations.

Functions:
- solcross_ut: Sun crossing events (e.g., ingresses, equinoxes)
- mooncross_ut: Moon crossing events (for lunar mansion calculations)
- cross_ut: Generic planet crossing

Precision: Newton-Raphson convergence tolerance
Tolerance: 0.001 arcsecond (sub-milliarcsecond) for all bodies (Sun, Moon, planets)
Iterations: Adaptive based on planet speed:
    - 50 for Sun/fast planets (speed >= 0.1°/day)
    - 30 for Moon (rapid convergence due to high speed ~13°/day)
    - 60 for slow planets like Saturn (0.01 <= speed < 0.1°/day)
    - 80 for very slow planets like Pluto (0.001 <= speed < 0.01°/day)
    - 100 near retrograde stations (speed < 0.001°/day)
Typical convergence: 5-8 iterations for Sun, 7-12 for Moon (up to 20 near nodes)

Station Handling:
    When a planet is near a retrograde station (velocity ~0°/day), Newton-Raphson
    can have convergence problems due to dividing by near-zero derivatives. In these
    cases, the algorithm automatically switches to Brent's method, which is more
    robust as it only requires bracketing the root, not computing derivatives.

Algorithm: Initial linear estimate + Newton-Raphson refinement + Brent's fallback
References:
    - Meeus "Astronomical Algorithms" Ch. 5 (interpolation)
    - Brent, R.P. (1973) "Algorithms for Minimization without Derivatives"
"""

from __future__ import annotations

import math
from typing import Callable, Tuple

from .exceptions import Error
from .constants import FLG_SWIEPH, FLG_SPEED, FLG_HELCTR, FLG_RADIANS, SUN, MOON
from .eclipse import _coerce_backwards
from .exceptions import EphemerisRangeError, CalculationError
from .planets import calc_ut, calc


def _normalize_cross_target(x2cross: float, flags: int) -> Tuple[float, int]:
    """Resolve ``FLG_RADIANS`` on a crossing target longitude.

    The reference interprets ``x2cross`` in radians when ``FLG_RADIANS`` is
    set (measured: ``solcross_ut(radians(90), ..., FLG_RADIANS)`` returns the
    same instant as the plain 90-degree call). The crossing math here is
    degree-based, so the target is converted up front and the flag stripped
    from the flags handed to calc()/calc_ut() — otherwise the returned
    longitude would come back in radians and corrupt the Newton/wrap steps.
    """
    if flags & FLG_RADIANS:
        x2cross = math.degrees(x2cross)
        flags &= ~FLG_RADIANS
    return x2cross % 360.0, flags


# Station detection threshold: speed below this indicates proximity to retrograde station
# At stations, Newton-Raphson can fail due to near-zero derivative (speed)
# Typical station speeds: Mercury ~0.05°/day slowing to 0, outer planets much slower
# Jupiter near station can be ~0.003°/day, Saturn ~0.001°/day — both need Brent's method
STATION_SPEED_THRESHOLD = 0.01  # degrees/day

# Newton-Raphson convergence constants
# Sub-milliarcsecond tolerance for all bodies — ensures timing precision
# limited only by the underlying ephemeris model, not by NR convergence.
# At Moon's ~13°/day, 0.001" ≈ 0.006ms timing; at Saturn's ~0.03°/day,
# 0.001" ≈ 0.3ms timing. This eliminates NR-induced timing jitter.
NR_TOLERANCE = 0.001 / 3600.0  # 0.001 arcsecond in degrees (all planets)
NR_TOLERANCE_SUN = 0.001 / 3600.0  # 0.001 arcsecond in degrees (Sun)
NR_TOLERANCE_MOON = 0.001 / 3600.0  # 0.001 arcsecond in degrees (Moon)
NR_MAX_ITER_SUN = 50  # Max iterations for Sun
# Moon iterations: 50 provides sufficient margin for the tighter 0.001"
# tolerance. Longitude crossings: Moon moves ~13°/day, typical 7-12 iters.
# Node crossings: latitude speed ~1°/day at nodes, up to 20-25 iterations.
NR_MAX_ITER_MOON = 50  # Max iterations for Moon
NR_MAX_ITER_PLANET = 60  # Max iterations for generic planets
NR_MAX_ITER_HELIO = 80  # Max iterations for heliocentric (slow planets)
NR_MAX_ITER_VERY_SLOW = 80  # Max iterations for very slow planets (Pluto, TNOs)
NR_MAX_ITER_STATION = 100  # Max iterations near retrograde stations


def _get_adaptive_max_iterations(speed: float) -> int:
    """
    Calculate adaptive iteration limit based on planet speed.

    For slow planets like Pluto (~0.004°/day) or near retrograde stations
    (speed approaching 0), more iterations are needed for Newton-Raphson
    convergence.

    Args:
        speed: Planet speed in degrees/day

    Returns:
        int: Maximum iterations to use

    Speed thresholds:
        |speed| >= 0.1: 50 iterations (normal planets)
        0.01 <= |speed| < 0.1: 60 iterations (slow planets like outer planets)
        0.001 <= |speed| < 0.01: 80 iterations (very slow: Pluto, TNOs)
        |speed| < 0.001: 100 iterations (near station/retrograde)
    """
    abs_speed = abs(speed)

    if abs_speed >= 0.1:
        return NR_MAX_ITER_PLANET
    elif abs_speed >= 0.01:
        return NR_MAX_ITER_HELIO
    elif abs_speed >= 0.001:
        return NR_MAX_ITER_VERY_SLOW
    else:
        return NR_MAX_ITER_STATION


def _is_near_station(speed: float) -> bool:
    """
    Check if a planet is near a retrograde station.

    Near stations, the planet's speed approaches zero, causing Newton-Raphson
    to have convergence issues (division by near-zero). This function detects
    when we should switch to a more robust method like Brent's.

    Args:
        speed: Planet speed in degrees/day

    Returns:
        bool: True if near station (speed below threshold)
    """
    return abs(speed) < STATION_SPEED_THRESHOLD


# Minimal bracket width (days) below which Brent switches to bisection.
# This is the step-size guard of Brent's method and is a TIME quantity
# (~1 ms); the function tolerances (degrees or deg/day) must not be
# compared against day intervals.
_BRENT_TIME_EPS_DAYS = 1e-8


def _brent_find_crossing(
    get_position_func: Callable[[float], Tuple[float, float]],
    x2cross: float,
    jd_a: float,
    jd_b: float,
    tolerance: float,
    max_iter: int = 100,
) -> float:
    """
    Find exact crossing time using Brent's method.

    Brent's method combines bisection, secant, and inverse quadratic interpolation.
    It's more robust than Newton-Raphson when derivatives are unreliable (near stations).

    Args:
        get_position_func: Function(jd) -> (longitude, speed) in degrees
        x2cross: Target ecliptic longitude in degrees
        jd_a: Start of bracket (Julian Day)
        jd_b: End of bracket (Julian Day)
        tolerance: Convergence tolerance in degrees
        max_iter: Maximum iterations

    Returns:
        float: Julian Day of crossing

    Raises:
        Error: If no root is bracketed or convergence fails

    Algorithm:
        Brent's method is guaranteed to converge if the function changes sign
        in [a, b]. It uses:
        - Bisection for safety
        - Secant method for faster convergence when applicable
        - Inverse quadratic interpolation when even faster

    References:
        Brent, R. P. (1973). Algorithms for Minimization without Derivatives.
    """

    def f(jd: float) -> float:
        """Return angular difference from target (signed, -180 to 180)."""
        lon, _ = get_position_func(jd)
        diff = (lon - x2cross) % 360.0
        if diff > 180:
            diff -= 360
        return diff

    fa = f(jd_a)
    fb = f(jd_b)

    # Check if root is bracketed
    if fa * fb > 0:
        # Root not bracketed - try to expand the bracket
        # This can happen if we guessed the bracket wrong
        raise Error(f"Brent's method: root not bracketed. f(a)={fa:.6f}, f(b)={fb:.6f}")

    # Ensure |f(a)| >= |f(b)| (b is the better guess)
    if abs(fa) < abs(fb):
        jd_a, jd_b = jd_b, jd_a
        fa, fb = fb, fa

    c = jd_a
    fc = fa
    mflag = True
    d = 0.0  # Will be set before use

    for _ in range(max_iter):
        if abs(fb) < tolerance:
            return jd_b

        if fa != fc and fb != fc:
            # Inverse quadratic interpolation
            s = (
                jd_a * fb * fc / ((fa - fb) * (fa - fc))
                + jd_b * fa * fc / ((fb - fa) * (fb - fc))
                + c * fa * fb / ((fc - fa) * (fc - fb))
            )
        else:
            # Secant method
            s = jd_b - fb * (jd_b - jd_a) / (fb - fa)

        # Conditions for using bisection instead
        cond1 = not (
            (3 * jd_a + jd_b) / 4 < s < jd_b or jd_b < s < (3 * jd_a + jd_b) / 4
        )
        cond2 = mflag and abs(s - jd_b) >= abs(jd_b - c) / 2
        cond3 = not mflag and abs(s - jd_b) >= abs(c - d) / 2
        cond4 = mflag and abs(jd_b - c) < _BRENT_TIME_EPS_DAYS
        cond5 = not mflag and abs(c - d) < _BRENT_TIME_EPS_DAYS

        if cond1 or cond2 or cond3 or cond4 or cond5:
            # Bisection
            s = (jd_a + jd_b) / 2
            mflag = True
        else:
            mflag = False

        fs = f(s)
        d = c
        c = jd_b
        fc = fb

        if fa * fs < 0:
            jd_b = s
            fb = fs
        else:
            jd_a = s
            fa = fs

        # Ensure |f(a)| >= |f(b)|
        if abs(fa) < abs(fb):
            jd_a, jd_b = jd_b, jd_a
            fa, fb = fb, fa

    # Return best estimate if max iterations reached
    return jd_b


def _find_bracket_for_crossing(
    get_position_func: Callable[[float], Tuple[float, float]],
    x2cross: float,
    jd_start: float,
    jd_end: float,
    num_samples: int = 20,
    from_end: bool = False,
) -> Tuple[float, float]:
    """
    Find a bracket [jd_a, jd_b] where the crossing occurs.

    Samples the position at regular intervals to find where the difference
    from target changes sign (indicating a crossing).

    Args:
        get_position_func: Function(jd) -> (longitude, speed) in degrees
        x2cross: Target ecliptic longitude in degrees
        jd_start: Start of search interval
        jd_end: End of search interval
        num_samples: Number of samples to take
        from_end: If False (default), scan forward from ``jd_start`` and return
            the FIRST (oldest) crossing in the interval. If True, scan backward
            from ``jd_end`` and return the LAST (most recent) crossing, i.e. the
            one nearest ``jd_end``. A backward search wants the nearest PAST
            crossing, which is the one adjacent to ``jd_end`` — not the oldest.

    Returns:
        Tuple[float, float]: (jd_a, jd_b) bracket containing the crossing,
        always ordered so jd_a < jd_b.

    Raises:
        Error: If no crossing is found in the interval
    """
    step = (jd_end - jd_start) / num_samples

    def get_diff(jd: float) -> float:
        lon, _ = get_position_func(jd)
        diff = (lon - x2cross) % 360.0
        if diff > 180:
            diff -= 360
        return diff

    # The wrapping filter below rejects sign changes at the antipodal point
    # (target+180°): a genuine crossing passes smoothly through 0 (small values
    # on each side), while antipodal wrapping shows diffs near ±180° that jump
    # by ~360°. Reject sign changes where the jump exceeds 180°.
    if from_end:
        # Walk backward from jd_end so the first sign change we hit is the
        # crossing closest to jd_end (the nearest past crossing).
        next_jd = jd_end
        next_diff = get_diff(jd_end)
        for i in range(num_samples - 1, -1, -1):
            curr_jd = jd_start + i * step
            curr_diff = get_diff(curr_jd)

            if curr_diff * next_diff <= 0 and abs(curr_diff - next_diff) < 180:
                return (curr_jd, next_jd)

            next_jd = curr_jd
            next_diff = curr_diff

        raise Error(
            f"No crossing found in interval [{jd_start}, {jd_end}] for target {x2cross}°"
        )

    prev_jd = jd_start
    prev_diff = get_diff(jd_start)

    for i in range(1, num_samples + 1):
        curr_jd = jd_start + i * step
        curr_diff = get_diff(curr_jd)

        # Check for sign change (crossing).
        if prev_diff * curr_diff <= 0 and abs(prev_diff - curr_diff) < 180:
            return (prev_jd, curr_jd)

        prev_jd = curr_jd
        prev_diff = curr_diff

    raise Error(
        f"No crossing found in interval [{jd_start}, {jd_end}] for target {x2cross}°"
    )


def _forward_first_crossing(
    get_position_func: Callable[[float], Tuple[float, float]],
    x2cross: float,
    tjdut: float,
    max_range: float,
    tolerance: float,
) -> float:
    """First crossing at or after ``tjdut`` via a forward bracket scan + Brent.

    Forward-only by construction (scans only [tjdut, tjdut+max_range] and returns
    the first sign change). The sample step scales with ``max_range`` so the scan
    stays cheap even for slow outer planets whose first crossing of a target lying
    behind them can be a full orbit (years–centuries) ahead. Raises ``Error`` if
    no crossing lies within ``max_range``.
    """
    num_samples = min(4000, max(int(max_range / 0.5), 200))
    jd_a, jd_b = _find_bracket_for_crossing(
        get_position_func, x2cross, tjdut, tjdut + max_range, num_samples=num_samples
    )
    return _brent_find_crossing(get_position_func, x2cross, jd_a, jd_b, tolerance, 100)


def _nearest_past_helio_crossing(
    get_position_func: Callable[[float], Tuple[float, float]],
    x2cross: float,
    tjd: float,
    jd_newton: float,
    period_days: float,
    tolerance: float,
) -> float:
    """Nearest strictly-past heliocentric crossing, validating Newton's result.

    Newton-Raphson sized from a single linear guess can converge into the basin
    of a crossing a whole revolution too far back: near aphelion the instantaneous
    speed runs well below the mean, so ``diff / speed`` overshoots into the wrong
    period. Validate ``jd_newton`` by scanning the *interior* of ``(jd_newton,
    tjd)`` — both endpoints guarded by a small ``margin`` so that a body sitting
    *at* the target at ``tjd`` (which must resolve a full cycle back) and a
    just-passed crossing (seconds ago) are both preserved — for a genuine crossing
    nearer to ``tjd``. If one exists, return it; otherwise keep ``jd_newton``.
    """
    margin = max(1.0, 0.01 * period_days)
    lo = min(jd_newton, tjd) + margin
    hi = tjd - margin
    if hi - lo >= margin:
        try:
            num_samples = max(40, int((hi - lo) / 15.0))
            jd_a, jd_b = _find_bracket_for_crossing(
                get_position_func,
                x2cross,
                lo,
                hi,
                num_samples=num_samples,
                # Nearest past crossing: scan from the tjd end and take the first
                # (most recent) sign change.
                from_end=True,
            )
            return _brent_find_crossing(
                get_position_func, x2cross, jd_a, jd_b, tolerance, 100
            )
        except (Error, RuntimeError):
            pass
    return jd_newton


# Retrograde-loop arc (deg) each planet sweeps: a target within this arc of the
# planet is reachable via an imminent retrograde dip rather than a full orbit
# ahead. Also reused as the "fine scan band" of the forward first-crossing guard.
_RETRO_ARC = {2: 16.0, 3: 18.0, 4: 22.0, 5: 12.0, 6: 9.0, 7: 6.0, 8: 5.0, 9: 4.0}
# Synodic period (days) — the spacing between successive retrograde loops.
_SYNODIC = {
    2: 116.0,
    3: 584.0,
    4: 780.0,
    5: 399.0,
    6: 378.0,
    7: 370.0,
    8: 367.0,
    9: 367.0,
}

# Planets whose Newton seed (a single linear estimate diff/speed) can converge
# into a crossing basin PAST the first true crossing, so the converged result is
# re-validated against an earlier crossing (see ``_guarded_first_crossing``).
# Slow outer planets (Jupiter..Pluto) skip whole retrograde loops. Mercury, Venus
# and Mars were previously excluded on the assumption that the forward-scan
# window already spanned their first crossing, but a momentarily slow (pre-
# station) inner planet still lets Newton overshoot diff/speed into the NEXT
# synodic cycle's crossing (measured skips of tens-to-hundreds of days). The
# guard's forward re-validation is cheap (it early-exits on the first crossing),
# so it is applied to them too.
_SLOW_GUARD_PLANETS = (2, 3, 4, 5, 6, 7, 8, 9)  # Mercury..Pluto
# Safe UPPER bound on each body's peak geocentric daily motion (deg/day; measured
# peaks Mer 2.202, Ven 1.259, Mars 0.776, Jup 0.242, Sat 0.130, Ura 0.061,
# Nep/Plu 0.039, padded ~10-25%). Used to size a forward-scan step that provably
# cannot step over a crossing: while the body is more than ``_RETRO_ARC`` degrees
# from the target, a step of (gap - arc) / max_rate advances it by at most
# (gap - arc) degrees, so it cannot reach the target mid-step.
_SLOW_MAX_RATE = {2: 2.5, 3: 1.35, 4: 0.85, 5: 0.30, 6: 0.16, 7: 0.08, 8: 0.05, 9: 0.05}


def _first_forward_crossing_adaptive(
    get_position_func: Callable[[float], Tuple[float, float]],
    x2cross: float,
    jd_lo: float,
    jd_hi: float,
    max_rate: float,
    band: float,
    fine_step: float = 2.0,
    coarse_cap: float = 300.0,
) -> float | None:
    """Earliest longitude crossing of ``x2cross`` in ``[jd_lo, jd_hi]``, or None.

    Adaptive step, cheap yet provably skip-free:

    * While the planet is more than ``band`` degrees from the target, the step
      covers at most the time to close that gap at ``max_rate`` (a safe upper
      bound on the body's peak daily motion), so the target cannot be reached
      mid-step and no crossing is jumped.
    * Within ``band`` — the neighbourhood of every crossing, including a
      retrograde triple-crossing clustered around a station — the step drops to
      ``fine_step``. Since d(delta)/dt equals the planet's speed away from the
      antipodal wrap (guaranteed in-band, because band << 180), delta is
      monotone between fine samples unless the speed changes sign there: a
      plain sign change thus brackets exactly one crossing. The only crossings
      a fixed grid could still miss are *grazing pairs* around a station where
      delta pokes past zero and back between two samples (arbitrarily narrow —
      e.g. a 0.1-arcsec, half-day graze at a station landing on the target).
      Those are caught exactly: when the speed flips sign between fine samples,
      the station (delta's extremum) is refined by bisection on the speed and
      the extremum's sign decides whether the pair exists.

    Returns the FIRST crossing, refined by bisection.
    """

    def _sample(jd_time: float) -> Tuple[float, float]:
        lon, spd = get_position_func(jd_time)
        return (lon - x2cross + 180.0) % 360.0 - 180.0, spd

    def _delta(jd_time: float) -> float:
        return _sample(jd_time)[0]

    def _bisect(lo: float, hi: float, f_lo: float) -> float:
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            f_mid = _delta(mid)
            if f_lo * f_mid <= 0:
                hi = mid
            else:
                lo, f_lo = mid, f_mid
        return 0.5 * (lo + hi)

    d_prev, v_prev = _sample(jd_lo)
    t_prev = jd_lo
    t = jd_lo
    while t < jd_hi:
        gap = abs(d_prev)
        in_band = gap <= band
        if in_band:
            step = fine_step
        else:
            step = min(coarse_cap, max(fine_step, (gap - band) / max_rate))
        t = min(t + step, jd_hi)
        d_cur, v_cur = _sample(t)
        if d_prev * d_cur < 0 and abs(d_cur - d_prev) < 180.0:
            return _bisect(t_prev, t, d_prev)
        if in_band and v_prev * v_cur < 0:
            # A station lies inside this fine interval, so delta has an
            # extremum here that may graze past the target and back between
            # the two samples (a hidden double crossing). Refine the station
            # by bisection on the speed, then let the extremum's sign decide.
            s_lo, s_hi, sv_lo = t_prev, t, v_prev
            for _ in range(50):
                s_mid = 0.5 * (s_lo + s_hi)
                v_mid = _sample(s_mid)[1]
                if sv_lo * v_mid <= 0:
                    s_hi = s_mid
                else:
                    s_lo, sv_lo = s_mid, v_mid
            t_station = 0.5 * (s_lo + s_hi)
            d_station = _delta(t_station)
            if d_station == 0.0:
                return t_station
            if d_prev * d_station < 0:
                # The extremum overshoots the target: the first of the
                # grazing pair lies in [t_prev, t_station].
                return _bisect(t_prev, t_station, d_prev)
        t_prev, d_prev, v_prev = t, d_cur, v_cur
    return None


def _guarded_first_crossing(
    get_position_func: Callable[[float], Tuple[float, float]],
    planet: int,
    x2cross: float,
    tjdut: float,
    jd_candidate: float,
) -> float:
    """Return the earliest crossing at or after ``tjdut``, re-validating a slow
    planet's far Newton result against an earlier crossing it may have skipped.

    The contract of ``cross_ut`` is the FIRST crossing for ``t >= tjdut``. Newton
    seeded from a single linear estimate ``diff / speed`` can land in the wrong
    retrograde loop for a slow outer planet with a distant target — the net
    progress is retarded by the retrograde loops, so the seed under-/over-shoots,
    and Newton then converges on a later crossing (or the 2nd/3rd of a retrograde
    triple), silently violating "first". The post-Newton guards only checked
    ">= tjdut" and "within max_range", never "earliest".

    This guard engages only for the slow on-table planets (Mercury..Mars move
    fast enough that the existing forward-scan window already spans their first
    crossing). A near-start retrograde triple-crossing can place the skipped
    first crossing only a few hundred days before the candidate, so — unlike a
    fixed synodic gate, which let those slip through — the guard runs whenever the
    candidate is more than a few days ahead, scanning ``[tjdut, jd_candidate)``
    for a genuinely earlier crossing. The scan is cheap (a coarse skip-free
    approach plus a short fine pass near the target that early-exits on the first
    crossing), so re-validating even the common single-crossing case is
    inexpensive; a more precise Newton value is kept when no earlier crossing
    exists.
    """
    if planet not in _SLOW_GUARD_PLANETS:
        return jd_candidate
    # Degenerate/trivial window: Newton is unambiguous for an imminent crossing,
    # and a zero-width scan would be meaningless.
    if jd_candidate - tjdut <= 5.0:
        return jd_candidate
    earlier = _first_forward_crossing_adaptive(
        get_position_func,
        x2cross,
        tjdut,
        jd_candidate - 0.5,  # exclude the candidate crossing itself
        _SLOW_MAX_RATE.get(planet, 0.16),
        _RETRO_ARC.get(planet, 9.0),
    )
    if earlier is not None and tjdut - 1e-6 <= earlier < jd_candidate - 1.0:
        return earlier
    return jd_candidate


def _cross_max_range(planet: int, speed_default: float) -> float:
    """Forward search horizon (days) for a planet's longitude crossing, sized to
    cover at least one orbital period so a target lying behind a slow planet is
    still found ahead."""
    if abs(speed_default) < 0.01:
        return 100000.0  # very slow (Pluto ~248 yr period)
    if abs(speed_default) < 0.05:
        return 40000.0  # Saturn/Uranus/Neptune
    if abs(speed_default) < 0.1:
        return 5000.0  # Jupiter (~12 yr)
    if planet in (2, 3):
        return 500.0  # Mercury, Venus
    return 800.0  # Mars and others


def solcross_ut(
    x2cross: float,
    tjdut: float,
    flags: int = FLG_SWIEPH,
    backwards: bool | int | str = False,
) -> float:
    """
    Find when the Sun crosses a specific ecliptic longitude.

    Searches forward (or backward if ``backwards=True``) in time from ``tjdut``.

    Args:
        x2cross: Target ecliptic longitude in degrees (0-360)
        tjdut: Julian Day (UT) to start search from
        flags: Calculation flags (FLG_SWIEPH, etc.)
        backwards: If True, search backward in time instead of forward.
            Also accepts ints and the direction strings understood by
            :func:`libephemeris.eclipse._coerce_backwards` (e.g.
            ``"forward"``, ``"backward"``).

    Returns:
        float: Julian Day of crossing (UT)

    Raises:
        Error: If convergence fails or calculation error occurs

    Algorithm:
        1. Get current Sun position and velocity
        2. Linear estimate: dt = (target - current) / velocity
        3. Refine with Newton-Raphson: jd_new = jd + (target - actual) / velocity
        4. Converge to < 0.001 arcsecond (sub-milliarcsecond precision)

    Precision:
        Typically < 0.001 arcsecond (< 0.03 seconds of time for Sun)

    Example:
        >>> # Find next Aries ingress (0°)
        >>> jd_ingress = solcross_ut(0.0, jd_now)
        >>> # Find summer solstice (90°)
        >>> jd_solstice = solcross_ut(90.0, jd_now)
        >>> # Find PREVIOUS Aries ingress
        >>> jd_prev = solcross_ut(0.0, jd_now, backwards=True)
    """
    backwards = _coerce_backwards(backwards)
    x2cross, flags = _normalize_cross_target(x2cross, flags)

    try:
        pos, _ = calc_ut(tjdut, SUN, flags | FLG_SPEED)
        lon_start = pos[0]
        speed = pos[3]  # degrees/day
    except (EphemerisRangeError, CalculationError, ValueError) as e:
        raise Error(f"Failed to calculate Sun position: {e}") from e

    # Initial time estimate (linear approximation)
    if speed == 0:
        speed = 0.9856  # Average Sun motion ~1°/day

    # Calculate angular distance to target — direction-aware
    diff = (x2cross - lon_start) % 360.0  # in [0, 360)

    if backwards:
        # Map diff from [0, 360) into (-360, 0].
        # When the Newton-Raphson forward result is passed back in, lon_start
        # lands within NR_TOLERANCE (~2.78e-7°) of the target — and can sit on
        # EITHER side, which means `diff` wraps close to 0 *or* close to 360.
        # Normalize to signed space first so the "at-crossing" guard catches
        # both sides. Threshold 1e-6° keeps a safety margin above
        # NR_TOLERANCE yet stays well below 1 second of apparent motion
        # (~1.14e-5°/s for Sun, ~1.52e-4°/s for Moon), so callers querying
        # the previous crossing moments after an event still get the nearby
        # event rather than a full cycle back.
        signed = diff - 360.0 if diff > 180.0 else diff
        if abs(signed) < 1e-6:
            diff = -360.0
        else:
            diff -= 360.0
    else:
        # Handle retrograde motion (Sun is never retrograde geocentrically,
        # but keep the original guard for robustness)
        if speed < 0 and diff > 0:
            diff -= 360.0

        # No dead-band skip here: (x2cross - lon) % 360 already maps a
        # just-passed target to ~360 (next cycle), while a body exactly
        # at or seconds before the target legitimately crosses now —
        # the reference returns the immediate crossing (verified:
        # swe.solcross_ut from the exact crossing instant returns that
        # instant, not one cycle later).

    dt_guess = diff / speed
    jd_guess = tjdut + dt_guess

    # Newton-Raphson iteration
    jd = jd_guess
    for iteration in range(NR_MAX_ITER_SUN):
        try:
            pos, _ = calc_ut(jd, SUN, flags | FLG_SPEED)
            lon = pos[0]
            speed = pos[3]
        except (EphemerisRangeError, CalculationError, ValueError) as e:
            raise Error(
                f"Failed to calculate Sun position during iteration: {e}"
            ) from e

        # Angular difference to target
        diff = (x2cross - lon) % 360.0
        if diff > 180:
            diff -= 360

        # Check convergence (< 0.001 arcsecond = 0.001/3600 degree for Sun)
        if abs(diff) < NR_TOLERANCE_SUN:
            return jd

        # Newton-Raphson step
        if abs(speed) < 0.01:
            speed = 0.9856

        jd += diff / speed

        # Safety: prevent divergence
        if abs(jd - jd_guess) > 366:
            raise Error("Solar crossing search diverged")

    raise Error("Maximum iterations reached in solar crossing calculation")


def solcross(
    x2cross: float,
    tjdet: float,
    flags: int = FLG_SWIEPH,
    backwards: bool | int | str = False,
) -> float:
    """
    Find when the Sun crosses a specific ecliptic longitude (TT version).

    This is the Terrestrial Time version of solcross_ut(). Takes Julian Day
    in TT (Terrestrial Time, also known as Ephemeris Time) instead of UT.

    Searches forward (or backward if ``backwards=True``) in time from ``tjdet``.

    Args:
        x2cross: Target ecliptic longitude in degrees (0-360)
        tjdet: Julian Day in Terrestrial Time (TT/ET) to start search from
        flags: Calculation flags (FLG_SWIEPH, etc.)
        backwards: If True, search backward in time instead of forward.
            Also accepts ints and the direction strings understood by
            :func:`libephemeris.eclipse._coerce_backwards` (e.g.
            ``"forward"``, ``"backward"``).

    Returns:
        float: Julian Day of crossing (TT)

    Raises:
        Error: If convergence fails or calculation error occurs

    Note:
        TT (Terrestrial Time) differs from UT (Universal Time) by Delta T,
        which varies from ~32 seconds (year 2000) to minutes (historical times).
        For most astrological applications, use solcross_ut() instead.

    Algorithm:
        1. Get current Sun position and velocity
        2. Linear estimate: dt = (target - current) / velocity
        3. Refine with Newton-Raphson: jd_new = jd + (target - actual) / velocity
        4. Converge to < 0.001 arcsecond (sub-milliarcsecond precision)

    Precision:
        Typically < 0.001 arcsecond (< 0.03 seconds of time for Sun)

    Example:
        >>> # Find next Aries ingress (0°) using TT
        >>> jd_ingress_tt = solcross(0.0, jd_tt_now)
        >>> # Find summer solstice (90°)
        >>> jd_solstice_tt = solcross(90.0, jd_tt_now)
        >>> # Find previous Aries ingress using TT
        >>> jd_prev_tt = solcross(0.0, jd_tt_now, backwards=True)
    """
    backwards = _coerce_backwards(backwards)
    x2cross, flags = _normalize_cross_target(x2cross, flags)

    try:
        pos, _ = calc(tjdet, SUN, flags | FLG_SPEED)
        lon_start = pos[0]
        speed = pos[3]  # degrees/day
    except (EphemerisRangeError, CalculationError, ValueError) as e:
        raise Error(f"Failed to calculate Sun position: {e}") from e

    # Initial time estimate (linear approximation)
    if speed == 0:
        speed = 0.9856  # Average Sun motion ~1°/day

    # Calculate angular distance to target — direction-aware
    diff = (x2cross - lon_start) % 360.0

    if backwards:
        # See solcross_ut for the signed-normalization rationale.
        signed = diff - 360.0 if diff > 180.0 else diff
        if abs(signed) < 1e-6:
            diff = -360.0
        else:
            diff -= 360.0
    else:
        # Handle retrograde motion (Sun/Moon prograde, guard for safety)
        if speed < 0 and diff > 0:
            diff -= 360.0

        # No dead-band skip here: (x2cross - lon) % 360 already maps a
        # just-passed target to ~360 (next cycle), while a body exactly
        # at or seconds before the target legitimately crosses now —
        # the reference returns the immediate crossing (verified:
        # swe.solcross_ut from the exact crossing instant returns that
        # instant, not one cycle later).

    dt_guess = diff / speed
    jd_guess = tjdet + dt_guess

    jd = jd_guess
    for iteration in range(NR_MAX_ITER_SUN):
        try:
            pos, _ = calc(jd, SUN, flags | FLG_SPEED)
            lon = pos[0]
            speed = pos[3]
        except (EphemerisRangeError, CalculationError, ValueError) as e:
            raise Error(
                f"Failed to calculate Sun position during iteration: {e}"
            ) from e

        # Angular difference to target
        diff = (x2cross - lon) % 360.0
        if diff > 180:
            diff -= 360

        # Check convergence (< 0.001 arcsecond for Sun)
        if abs(diff) < NR_TOLERANCE_SUN:
            return jd

        # Newton-Raphson step
        if abs(speed) < 0.01:
            speed = 0.9856

        jd += diff / speed

        # Safety: prevent divergence
        if abs(jd - jd_guess) > 366:
            raise Error("Solar crossing search diverged")

    raise Error("Maximum iterations reached in solar crossing calculation")


def mooncross_ut(
    x2cross: float,
    tjdut: float,
    flags: int = FLG_SWIEPH,
    backwards: bool | int | str = False,
) -> float:
    """
    Find when the Moon crosses a specific ecliptic longitude.

    Searches forward (or backward if ``backwards=True``) in time from ``tjdut``.

    Args:
        x2cross: Target ecliptic longitude in degrees (0-360)
        tjdut: Julian Day (UT) to start search from
        flags: Calculation flags (FLG_SWIEPH, etc.)
        backwards: If True, search backward in time instead of forward.
            Also accepts ints and the direction strings understood by
            :func:`libephemeris.eclipse._coerce_backwards` (e.g.
            ``"forward"``, ``"backward"``).

    Returns:
        float: Julian Day of crossing (UT)

    Raises:
        Error: If convergence fails or calculation error occurs

    Note:
        Moon moves ~13° per day (27.3 day cycle).
        More variable speed than Sun due to orbit eccentricity.

    Precision:
        Typically < 0.001 arcsecond (sub-milliarcsecond)

    Example:
        >>> # Find next new moon (Sun-Moon conjunction at same longitude)
        >>> sun_pos, _ = calc_ut(jd_now, SUN, FLG_SWIEPH)
        >>> jd_new_moon = mooncross_ut(sun_pos[0], jd_now)
        >>> # Find PREVIOUS Moon crossing of 0°
        >>> jd_prev = mooncross_ut(0.0, jd_now, backwards=True)
    """
    backwards = _coerce_backwards(backwards)
    x2cross, flags = _normalize_cross_target(x2cross, flags)

    try:
        pos, _ = calc_ut(tjdut, MOON, flags | FLG_SPEED)
        lon_start = pos[0]
        speed = pos[3]  # degrees/day
    except (EphemerisRangeError, CalculationError, ValueError) as e:
        raise Error(f"Failed to calculate Moon position: {e}") from e

    # Initial time estimate
    if speed == 0:
        speed = 13.176  # Average Moon motion ~13.18°/day

    # Calculate initial guess — direction-aware
    diff = (x2cross - lon_start) % 360.0

    if backwards:
        # See solcross_ut for the signed-normalization rationale.
        signed = diff - 360.0 if diff > 180.0 else diff
        if abs(signed) < 1e-6:
            diff = -360.0
        else:
            diff -= 360.0
    else:
        if speed < 0 and diff > 0:
            diff -= 360.0

        # No dead-band skip here: (x2cross - lon) % 360 already maps a
        # just-passed target to ~360 (next cycle), while a body exactly
        # at or seconds before the target legitimately crosses now —
        # the reference returns the immediate crossing (verified:
        # swe.solcross_ut from the exact crossing instant returns that
        # instant, not one cycle later).

    dt_guess = diff / speed
    jd_guess = tjdut + dt_guess

    # Newton-Raphson iteration
    jd = jd_guess
    for iteration in range(NR_MAX_ITER_MOON):
        try:
            pos, _ = calc_ut(jd, MOON, flags | FLG_SPEED)
            lon = pos[0]
            speed = pos[3]
        except (EphemerisRangeError, CalculationError, ValueError) as e:
            raise Error(
                f"Failed to calculate Moon position during iteration: {e}"
            ) from e

        # Difference to target
        diff = (x2cross - lon) % 360.0
        if diff > 180:
            diff -= 360

        # Check convergence (< 0.001 arcsecond for Moon)
        if abs(diff) < NR_TOLERANCE_MOON:
            return jd

        # Newton-Raphson step
        if abs(speed) < 0.1:
            speed = 13.176

        jd += diff / speed

        # Safety check
        if abs(jd - jd_guess) > 31:  # More than a month
            raise Error("Moon crossing search diverged")

    raise Error("Maximum iterations reached in moon crossing calculation")


def mooncross(
    x2cross: float,
    tjdet: float,
    flags: int = FLG_SWIEPH,
    backwards: bool | int | str = False,
) -> float:
    """
    Find when the Moon crosses a specific ecliptic longitude (TT version).

    This is the Terrestrial Time version of mooncross_ut(). Takes Julian Day
    in TT (Terrestrial Time, also known as Ephemeris Time) instead of UT.

    Searches forward (or backward if ``backwards=True``) in time from ``tjdet``.

    Args:
        x2cross: Target ecliptic longitude in degrees (0-360)
        tjdet: Julian Day in Terrestrial Time (TT/ET) to start search from
        flags: Calculation flags (FLG_SWIEPH, etc.)
        backwards: If True, search backward in time instead of forward.
            Also accepts ints and the direction strings understood by
            :func:`libephemeris.eclipse._coerce_backwards` (e.g.
            ``"forward"``, ``"backward"``).

    Returns:
        float: Julian Day of crossing (TT)

    Raises:
        Error: If convergence fails or calculation error occurs

    Note:
        TT (Terrestrial Time) differs from UT (Universal Time) by Delta T,
        which varies from ~32 seconds (year 2000) to minutes (historical times).
        For most astrological applications, use mooncross_ut() instead.

        Moon moves ~13° per day (27.3 day cycle).
        More variable speed than Sun due to orbit eccentricity.

    Algorithm:
        1. Get current Moon position and velocity
        2. Linear estimate: dt = (target - current) / velocity
        3. Refine with Newton-Raphson: jd_new = jd + (target - actual) / velocity
        4. Converge to < 0.001 arcsecond (sub-milliarcsecond precision)

    Precision:
        Typically < 0.001 arcsecond (sub-milliarcsecond)

    Example:
        >>> # Find next new moon (Sun-Moon conjunction at same longitude) using TT
        >>> sun_pos, _ = calc(jd_tt_now, SUN, FLG_SWIEPH)
        >>> jd_new_moon_tt = mooncross(sun_pos[0], jd_tt_now)
        >>> # Find previous Moon crossing of 0° using TT
        >>> jd_prev_tt = mooncross(0.0, jd_tt_now, backwards=True)
    """
    backwards = _coerce_backwards(backwards)
    x2cross, flags = _normalize_cross_target(x2cross, flags)

    try:
        pos, _ = calc(tjdet, MOON, flags | FLG_SPEED)
        lon_start = pos[0]
        speed = pos[3]  # degrees/day
    except (EphemerisRangeError, CalculationError, ValueError) as e:
        raise Error(f"Failed to calculate Moon position: {e}") from e

    # Initial time estimate
    if speed == 0:
        speed = 13.176  # Average Moon motion ~13.18°/day

    # Calculate initial guess — direction-aware
    diff = (x2cross - lon_start) % 360.0

    if backwards:
        # See solcross_ut for the signed-normalization rationale.
        signed = diff - 360.0 if diff > 180.0 else diff
        if abs(signed) < 1e-6:
            diff = -360.0
        else:
            diff -= 360.0
    else:
        if speed < 0 and diff > 0:
            diff -= 360.0

        # No dead-band skip here: (x2cross - lon) % 360 already maps a
        # just-passed target to ~360 (next cycle), while a body exactly
        # at or seconds before the target legitimately crosses now —
        # the reference returns the immediate crossing (verified:
        # swe.solcross_ut from the exact crossing instant returns that
        # instant, not one cycle later).

    dt_guess = diff / speed
    jd_guess = tjdet + dt_guess

    jd = jd_guess
    for iteration in range(NR_MAX_ITER_MOON):
        try:
            pos, _ = calc(jd, MOON, flags | FLG_SPEED)
            lon = pos[0]
            speed = pos[3]
        except (EphemerisRangeError, CalculationError, ValueError) as e:
            raise Error(
                f"Failed to calculate Moon position during iteration: {e}"
            ) from e

        # Difference to target
        diff = (x2cross - lon) % 360.0
        if diff > 180:
            diff -= 360

        # Check convergence (< 0.001 arcsecond for Moon)
        if abs(diff) < NR_TOLERANCE_MOON:
            return jd

        # Newton-Raphson step
        if abs(speed) < 0.1:
            speed = 13.176

        jd += diff / speed

        # Safety check
        if abs(jd - jd_guess) > 31:  # More than a month
            raise Error("Moon crossing search diverged")

    raise Error("Maximum iterations reached in moon crossing calculation")


def mooncross_node_ut(
    tjdut: float,
    flags: int = FLG_SWIEPH,
    backwards: bool | int | str = False,
) -> Tuple[float, float, float]:
    """
    Find when the Moon crosses its own orbital node (ascending or descending).

    The Moon crosses a node when its ecliptic latitude becomes zero - i.e., when
    it crosses the ecliptic plane. This is important for eclipse calculations,
    as eclipses can only occur when the Sun and Moon are near the lunar nodes.

    Searches forward (or backward if ``backwards=True``) in time from ``tjdut``.

    Args:
        tjdut: Julian Day (UT) to start search from
        flags: Calculation flags (FLG_SWIEPH, etc.)
        backwards: If True, search backward in time instead of forward.
            Also accepts ints and the direction strings understood by
            :func:`libephemeris.eclipse._coerce_backwards` (e.g.
            ``"forward"``, ``"backward"``).

    Returns:
        tuple: (jd_cross, xlon, xlat) - Julian Day of crossing, ecliptic longitude,
               and ecliptic latitude (should be ~0) at the crossing moment

    Raises:
        Error: If convergence fails or calculation error occurs

    Algorithm:
        1. Convert UT to TT for internal calculation
        2. Delegate to mooncross_node (TT version) for Newton-Raphson search
        3. Convert result TT back to UT for return value

    Note:
        The function finds the NEXT (or PREVIOUS when ``backwards=True``) node
        crossing regardless of whether it's ascending (latitude going from
        negative to positive) or descending (positive to negative).

        Moon crosses each node approximately every 13.6 days (half the nodal
        month of ~27.2 days).

        Both engines solve the SAME event — the Moon's ecliptic-latitude-zero
        crossing (measured black-box: identical instants). The only
        divergence is the time frame of the *_ut return value: the reference
        reports the TT/ET instant unconverted, while libephemeris reports
        the true UT instant, so the gap is exactly Delta-T (~64 s at J2000,
        growing at extreme epochs). See known-differences.md §11.

    Example:
        >>> # Find next lunar node crossing
        >>> jd_node, xlon, xlat = mooncross_node_ut(jd_now)
        >>> # Check if ascending or descending by examining latitude velocity
        >>> pos, _ = calc_ut(jd_node, MOON, FLG_SPEED)
        >>> is_ascending = pos[4] > 0  # positive lat velocity = ascending
        >>> # Find PREVIOUS lunar node crossing
        >>> jd_prev, _, _ = mooncross_node_ut(jd_now, backwards=True)
    """
    backwards = _coerce_backwards(backwards)

    from .state import get_timescale

    # Convert UT to TT for the TT-based search
    ts = get_timescale()
    t_start = ts.ut1_jd(tjdut)
    jd_tt_start = t_start.tt

    # Delegate to the TT version
    jd_tt_cross, xlon, xlat = mooncross_node(jd_tt_start, flags, backwards=backwards)

    # Convert result from TT back to UT
    t_cross = ts.tt_jd(jd_tt_cross)
    jd_ut_cross = t_cross.ut1

    return (float(jd_ut_cross), xlon, xlat)


def mooncross_node(
    tjdet: float,
    flags: int = FLG_SWIEPH,
    backwards: bool | int | str = False,
) -> Tuple[float, float, float]:
    """
    Find when the Moon crosses its own orbital node (TT version).

    This is the Terrestrial Time version of mooncross_node_ut(). Takes Julian Day
    in TT (Terrestrial Time, also known as Ephemeris Time) instead of UT.

    The Moon crosses a node when its ecliptic latitude becomes zero - i.e., when
    it crosses the ecliptic plane. This is important for eclipse calculations,
    as eclipses can only occur when the Sun and Moon are near the lunar nodes.

    Searches forward (or backward if ``backwards=True``) in time from ``tjdet``.

    Args:
        tjdet: Julian Day in Terrestrial Time (TT/ET) to start search from
        flags: Calculation flags (FLG_SWIEPH, etc.)
        backwards: If True, search backward in time instead of forward.
            Also accepts ints and the direction strings understood by
            :func:`libephemeris.eclipse._coerce_backwards` (e.g.
            ``"forward"``, ``"backward"``).

    Returns:
        tuple: (jd_cross, xlon, xlat) - Julian Day of crossing (TT), ecliptic longitude,
               and ecliptic latitude (should be ~0) at the crossing moment

    Raises:
        Error: If convergence fails or calculation error occurs

    Note:
        TT (Terrestrial Time) differs from UT (Universal Time) by Delta T,
        which varies from ~32 seconds (year 2000) to minutes (historical times).
        For most astrological applications, use mooncross_node_ut() instead.

        Both engines solve the SAME event — the Moon's ecliptic-latitude-zero
        crossing (measured black-box: identical instants). The only
        divergence is the time frame of the *_ut return value: the reference
        reports the TT/ET instant unconverted, while libephemeris reports
        the true UT instant, so the gap is exactly Delta-T (~64 s at J2000,
        growing at extreme epochs). See known-differences.md §11.

    Example:
        >>> # Find next lunar node crossing using TT
        >>> jd_node_tt, xlon, xlat = mooncross_node(jd_tt_now)
        >>> # Find previous lunar node crossing using TT
        >>> jd_prev_tt, _, _ = mooncross_node(jd_tt_now, backwards=True)
    """
    backwards = _coerce_backwards(backwards)

    # Half nodal month - time between successive node crossings
    HALF_NODAL_MONTH = 13.6
    direction = -1.0 if backwards else 1.0

    try:
        pos, _ = calc(tjdet, MOON, flags | FLG_SPEED)
        lat = pos[1]  # ecliptic latitude
    except (EphemerisRangeError, CalculationError, ValueError) as e:
        raise Error(f"Failed to calculate Moon position: {e}") from e

    # If starting essentially at a node crossing, the linear estimate is
    # dominated by FP noise and NR would converge on the start point
    # (failing the direction check and bouncing forever). Take a full
    # half-nodal-month step in the requested direction instead.
    if abs(lat) < 10 * NR_TOLERANCE_MOON:
        jd_guess = tjdet + HALF_NODAL_MONTH * direction
    else:
        # Locate the FIRST latitude sign change in the requested direction with
        # a bracket scan whose step provably cannot skip a root: successive Moon
        # node crossings are >= ~6.8 days apart, so a <= 3-day step samples
        # every interval and no root is jumped. A linear seed (-lat/lat_speed)
        # plus a coarse "needs_scan" gate could instead land Newton in the
        # *second* root's basin — skipping the first crossing by ~13 days — or
        # diverge outright; scanning from tjdet itself avoids both. Once
        # bracketed, the Newton/bisection loop below refines to full precision
        # and applies the strict-direction epsilon.
        SCAN_STEP = 2.0  # days; below the min root spacing, so no root is skipped
        SCAN_SPAN = HALF_NODAL_MONTH + 4.0  # first crossing lies within a half-cycle
        step_days = SCAN_STEP * direction
        lat_prev = lat
        t_prev = tjdet
        # Default seed if the scan finds no sign change within a half-cycle
        # (should not happen for the Moon): fall back to a half-nodal-month step.
        jd_guess = tjdet + HALF_NODAL_MONTH * direction
        n_steps = int(SCAN_SPAN / SCAN_STEP) + 1

        def _moon_lat(t: float) -> float:
            try:
                pos_t, _ = calc(t, MOON, flags | FLG_SPEED)
            except (EphemerisRangeError, CalculationError, ValueError) as exc:
                raise Error(f"Failed to calculate Moon position: {exc}") from exc
            return pos_t[1]

        for i in range(1, n_steps + 1):
            t_cur = tjdet + i * step_days
            lat_cur = _moon_lat(t_cur)
            if lat_prev * lat_cur < 0.0:
                # First sign change brackets the first crossing in
                # [t_prev, t_cur] (ordered by the search direction). Bisect the
                # bracket to the root before seeding Newton: a bare bracket
                # midpoint can land Newton on the far side of a crossing that
                # sits at the bracket edge, whose step then walks past tjdet and
                # the direction clamp bounces it into an adjacent latitude
                # extremum (division by a near-zero rate) and it never
                # converges. A root-refined seed keeps the Newton loop below on
                # this (first) crossing.
                lo, hi, f_lo = t_prev, t_cur, lat_prev
                for _ in range(60):
                    mid = 0.5 * (lo + hi)
                    f_mid = _moon_lat(mid)
                    if f_lo * f_mid <= 0.0:
                        hi = mid
                    else:
                        lo, f_lo = mid, f_mid
                jd_guess = 0.5 * (lo + hi)
                break
            t_prev, lat_prev = t_cur, lat_cur

    jd = jd_guess
    for iteration in range(NR_MAX_ITER_MOON):
        try:
            pos, _ = calc(jd, MOON, flags | FLG_SPEED)
            lat = pos[1]
            lat_speed = pos[4]
        except (EphemerisRangeError, CalculationError, ValueError) as e:
            raise Error(
                f"Failed to calculate Moon position during iteration: {e}"
            ) from e

        # Check convergence (< 0.001 arcsecond for Moon)
        if abs(lat) < NR_TOLERANCE_MOON:
            # Must be strictly past tjdet in the requested direction. The
            # epsilon (~0.09 s) is tight enough to accept legitimate results
            # seconds before/after tjdet yet wide enough to reject the
            # numerically identical "same event" case when the caller starts
            # at the exact node crossing itself.
            if (not backwards and jd > tjdet + 1e-6) or (
                backwards and jd < tjdet - 1e-6
            ):
                return (jd, pos[0], lat)
            # Converged on the wrong side of tjdet — nudge to the next
            # (or previous) crossing and keep iterating.
            jd = jd + (HALF_NODAL_MONTH / 2) * direction
            continue

        # Newton-Raphson step
        if abs(lat_speed) < 0.1:
            lat_speed = 1.0 if lat >= 0 else -1.0

        jd -= lat / lat_speed

        # Safety check — clamp to the requested direction
        if not backwards and jd < tjdet:
            jd = tjdet + HALF_NODAL_MONTH / 2
        elif backwards and jd > tjdet:
            jd = tjdet - HALF_NODAL_MONTH / 2
        elif abs(jd - tjdet) > 30:
            raise Error("Moon node crossing search diverged")

    raise Error("Maximum iterations reached in moon node crossing calculation")


def cross_ut(
    planet: int, x2cross: float, tjdut: float, flags: int = FLG_SWIEPH
) -> float:
    """
    Find when any planet crosses a specific ecliptic longitude.

    Generic version for all planets (Mercury, Venus, Mars, etc.).

    Args:
        planet: Planet ID (MERCURY, VENUS, etc.)
        x2cross: Target ecliptic longitude in degrees (0-360)
        tjdut: Julian Day (UT) to start search from
        flags: Calculation flags

    Returns:
        float: Julian Day of crossing (UT)

    Raises:
        Error: If convergence fails or calculation error occurs

    Note:
        Uses adaptive iteration count based on typical planet speed.
        Slower planets (Jupiter, Saturn) may need more iterations.

    Precision:
        Typically < 0.001 arcsecond

    Example:
        >>> # Mars ingress into Aries
        >>> jd_mars_aries = cross_ut(MARS, 0.0, jd_now)
    """
    x2cross, flags = _normalize_cross_target(x2cross, flags)

    try:
        pos, _ = calc_ut(tjdut, planet, flags | FLG_SPEED)
        lon_start = pos[0]
        speed = pos[3]
    except (EphemerisRangeError, CalculationError, ValueError) as e:
        raise Error(f"Failed to calculate planet position: {e}") from e

    # Estimate typical speed if near zero
    # Geocentric average speeds (°/day) - slower planets need more iterations
    # Note: geocentric speeds are affected by retrograde motion, values are
    # approximate averages during direct motion
    typical_speeds = {
        2: 1.4,  # Mercury
        3: 1.2,  # Venus
        4: 0.5,  # Mars
        5: 0.08,  # Jupiter
        6: 0.03,  # Saturn
        7: 0.01,  # Uranus
        8: 0.006,  # Neptune
        9: 0.004,  # Pluto (very slow, ~0.004°/day)
    }
    speed_default = typical_speeds.get(planet, 0.5)

    # Bodies outside the typical-speed table (e.g. Mean Node ~0.053°/day
    # retrograde with an ~18.6 yr period, Mean Apogee ~0.111°/day with an
    # ~8.85 yr period) would otherwise inherit the generic 0.5°/day fallback,
    # whose implied ~800 day search horizon is far shorter than their real
    # period — genuine crossings thousands of days out then look like
    # divergence. For a genuinely slow off-table body, derive a representative
    # speed and a search horizon (~1.2 periods) from the actual instantaneous
    # rate at tjdut. Table bodies are left exactly as before, and faster
    # off-table bodies (whose real period already fits the generic horizon) are
    # untouched. Sign is handled downstream: the forward scan below tracks the
    # wrapped longitude regardless of direction, so perpetually-retrograde
    # bodies work without special-casing.
    off_table_horizon: float | None = None
    if planet not in typical_speeds and 1e-6 < abs(speed) < 0.45:
        speed_default = abs(speed)
        off_table_horizon = min(1.2 * 360.0 / speed_default, 40000.0)

    # Calculate initial guess
    # For forward-only search, always calculate forward angular distance
    diff = (x2cross - lon_start) % 360.0

    # When retrograde, the planet will eventually turn direct and reach the target
    # So we estimate using typical (prograde) speed, not current retrograde speed
    if speed < 0:
        # Use average speed for time estimate (planet will turn direct)
        effective_speed = speed_default
    else:
        effective_speed = speed if abs(speed) > 0.001 else speed_default

    # No dead-band skip: (x2cross - lon) % 360 maps a just-passed
    # target to ~360 already, and a body exactly at the target crosses
    # NOW (same convention as solcross/mooncross).

    # Newton-Raphson assumes near-monotonic prograde motion and mishandles
    # crossings near a retrograde loop, in several ways that all return the wrong
    # crossing:
    #   * target just BEHIND now, reached via an imminent retrograde dip (even
    #     while still prograde, approaching a station) -> NR projects a full orbit
    #     ahead and returns a crossing months-to-decades too late;
    #   * currently RETROGRADE with a target it last passed in the recent past ->
    #     NR can step backward and converge to a PAST crossing (before tjdut),
    #     violating the forward-only contract;
    #   * near a station, the step diff/speed explodes and the search diverges.
    # In all of these the chronologically first crossing lies within ~one synodic
    # period, so scan forward from tjdut for the first wrapped sign change and
    # refine. The scan is forward-only by construction and early-exits on the
    # first crossing, so it is cheap whenever one exists soon. The fast common
    # case (prograde, target ahead, away from a station) skips it and uses NR.
    # (``_RETRO_ARC`` / ``_SYNODIC`` are module-level constants.)
    diff_back = (lon_start - x2cross) % 360.0
    if speed < 0 or _is_near_station(speed) or diff_back < _RETRO_ARC.get(planet, 22.0):
        scan_window = _SYNODIC.get(planet, 800.0) * 1.2
        if off_table_horizon is not None:
            # A slow off-table body (e.g. the perpetually-retrograde Mean Node)
            # can reach the target only after a large fraction of its long
            # period; widen the forward scan to cover it. The scan early-exits
            # on the first crossing, so this stays cheap when one lies near.
            scan_window = max(scan_window, off_table_horizon)

        def _delta(jd_time: float) -> float:
            lon_t = calc_ut(jd_time, planet, flags)[0][0]
            return (lon_t - x2cross + 180.0) % 360.0 - 180.0

        d_prev = _delta(tjdut)
        step = 0.5
        t_prev = tjdut
        t_scan = tjdut + step
        while t_scan <= tjdut + scan_window:
            d_cur = _delta(t_scan)
            if d_prev == 0.0:
                return float(t_prev)
            if d_prev * d_cur < 0 and abs(d_cur - d_prev) < 180.0:
                lo, hi = t_prev, t_scan
                f_lo = _delta(lo)
                for _ in range(60):
                    mid = 0.5 * (lo + hi)
                    f_mid = _delta(mid)
                    if f_lo * f_mid <= 0:
                        hi = mid
                    else:
                        lo = mid
                        f_lo = f_mid
                return float(0.5 * (lo + hi))
            t_prev, d_prev = t_scan, d_cur
            t_scan += step

    dt_guess = diff / effective_speed
    jd_guess = tjdut + dt_guess

    # Adaptive iteration count based on planet speed
    max_iter = _get_adaptive_max_iterations(speed)

    # Helper function for Brent's method fallback
    def get_position(jd_time: float) -> Tuple[float, float]:
        pos_result, _ = calc_ut(jd_time, planet, flags | FLG_SPEED)
        return pos_result[0], pos_result[3]

    # "First crossing" guard: any computed crossing is re-validated (for the slow
    # on-table planets, and only when it lies far ahead) against an earlier
    # crossing that Newton or a too-narrow bracket scan may have skipped. A near
    # result is returned untouched, so the fast common path is unaffected.
    def _first(jd_result: float) -> float:
        return _guarded_first_crossing(get_position, planet, x2cross, tjdut, jd_result)

    # Check if we're near a retrograde station - use Brent's method for robustness
    if _is_near_station(speed):
        # Near station: Newton-Raphson may fail due to division by near-zero speed
        # Use Brent's method which only requires bracketing, not derivatives
        try:
            # Estimate search window based on typical speeds.
            # Near stations or retrograde, the planet may need to complete a
            # full retrograde cycle before reaching the target. Use a generous
            # window: dt_guess * 2.0 covers the expected travel plus margin
            # for retrograde detours that can add ~50% more time.
            if speed_default > 0:
                min_window = 500.0 if speed_default < 0.1 else 60.0
                search_window = max(
                    min_window, abs(dt_guess) * 2.0, abs(diff / speed_default) * 3.0
                )
            else:
                search_window = max(500.0, abs(dt_guess) * 2.0)
            jd_bracket_start = tjdut
            jd_bracket_end = tjdut + search_window

            # Scale samples so we check at least every ~30 days
            bracket_samples = max(40, int(search_window / 30))

            # Find bracket containing the crossing
            jd_a, jd_b = _find_bracket_for_crossing(
                get_position,
                x2cross,
                jd_bracket_start,
                jd_bracket_end,
                num_samples=bracket_samples,
            )

            # Use Brent's method to find the exact crossing
            return _first(
                _brent_find_crossing(
                    get_position, x2cross, jd_a, jd_b, NR_TOLERANCE, max_iter
                )
            )
        except (Error, RuntimeError):
            # If Brent's method fails, fall through to Newton-Raphson as last resort
            pass

    # Newton-Raphson iteration
    jd = jd_guess
    station_fallback_triggered = False
    max_range = _cross_max_range(planet, speed_default)
    if off_table_horizon is not None:
        # Give the divergence guard and the forward-scan catch-all enough room to
        # reach a slow off-table body's genuine (far) crossing.
        max_range = max(max_range, off_table_horizon)

    for iteration in range(max_iter):
        try:
            pos, _ = calc_ut(jd, planet, flags | FLG_SPEED)
            lon = pos[0]
            speed = pos[3]
        except EphemerisRangeError as range_err:
            # Newton wandered outside the ephemeris coverage before any in-range
            # guard could fire. For a very slow outer planet ``max_range`` can
            # exceed the ephemeris span, so the ``|jd - tjdut| > max_range`` guard
            # never trips and Newton can step off the covered interval (e.g.
            # backward past the lower edge) mid-iteration. A genuine first
            # crossing may still lie in range ahead of tjdut, so fall back to the
            # forward bracket scan, which samples only [tjdut, tjdut+max_range]
            # and returns the earliest crossing. If no in-range crossing exists
            # (the scan itself runs off coverage), re-raise the ephemeris-range
            # error rather than masking it as a generic divergence.
            try:
                return _first(
                    _forward_first_crossing(
                        get_position, x2cross, tjdut, max_range, NR_TOLERANCE
                    )
                )
            except EphemerisRangeError:
                raise
            except (Error, RuntimeError):
                raise range_err
        except (CalculationError, ValueError) as e:
            raise Error(
                f"Failed to calculate planet position during iteration: {e}"
            ) from e

        diff = (x2cross - lon) % 360.0
        if diff > 180:
            diff -= 360

        # Check convergence (< 0.001 arcsecond)
        if abs(diff) < NR_TOLERANCE:
            if jd >= tjdut - 1e-6:
                return _first(jd)
            # Newton converged to a PAST crossing — the wrapped step diff/speed
            # walked backward (e.g. a target far behind a prograde planet, whose
            # first forward crossing is a whole orbit ahead). Forward-only is a
            # hard contract: return the first crossing at or after tjdut instead.
            try:
                return _first(
                    _forward_first_crossing(
                        get_position, x2cross, tjdut, max_range, NR_TOLERANCE
                    )
                )
            except EphemerisRangeError:
                # The forward scan ran off the ephemeris: no in-range crossing
                # exists, so surface the ephemeris-range condition rather than
                # masking it as a generic divergence.
                raise
            except (Error, RuntimeError):
                raise Error("Planet crossing search diverged")

        # Detect if we've encountered a station or retrograde during iteration.
        # Retrograde (speed < 0) causes NR to step backward in time, which
        # can oscillate indefinitely for forward-only searches.
        if (_is_near_station(speed) or speed < 0) and not station_fallback_triggered:
            station_fallback_triggered = True
            # Switch to Brent's method for robustness.
            # Search from the original start time forward — the planet may
            # need to complete a full retrograde cycle before reaching target.
            try:
                if speed_default > 0:
                    min_window = 500.0 if speed_default < 0.1 else 60.0
                    bracket_window = max(
                        min_window,
                        abs(dt_guess) * 2.0,
                        abs(diff / speed_default) * 3.0,
                    )
                else:
                    bracket_window = max(500.0, abs(dt_guess) * 2.0)
                bracket_samples = max(40, int(bracket_window / 30))
                jd_a, jd_b = _find_bracket_for_crossing(
                    get_position,
                    x2cross,
                    tjdut,
                    tjdut + bracket_window,
                    num_samples=bracket_samples,
                )
                return _first(
                    _brent_find_crossing(
                        get_position,
                        x2cross,
                        jd_a,
                        jd_b,
                        NR_TOLERANCE,
                        max_iter - iteration,
                    )
                )
            except (Error, RuntimeError):
                # If Brent fails, continue with Newton-Raphson
                pass

        # Update max_iter based on current speed (may change near stations)
        max_iter = max(max_iter, _get_adaptive_max_iterations(speed))

        if abs(speed) < 0.001:
            speed = speed_default

        jd += diff / speed

        # Safety: if Newton wandered beyond one orbital period from tjdut
        # (typically a near-station overshoot where diff/speed explodes as
        # speed -> 0), fall back to a forward bracket scan + Brent — a crossing
        # exists at or after tjdut within max_range. Forward-only by construction.
        if abs(jd - tjdut) > max_range:
            try:
                return _first(
                    _forward_first_crossing(
                        get_position, x2cross, tjdut, max_range, NR_TOLERANCE
                    )
                )
            except EphemerisRangeError:
                raise
            except (Error, RuntimeError):
                raise Error("Planet crossing search diverged")

    # Newton-Raphson exhausted its iterations without converging — e.g. a far
    # crossing of a very slow planet that NR keeps oscillating around within
    # max_range (so neither the divergence nor the backward guard fired). The
    # forward bracket scan is the robust catch-all.
    try:
        return _first(
            _forward_first_crossing(
                get_position, x2cross, tjdut, max_range, NR_TOLERANCE
            )
        )
    except EphemerisRangeError:
        raise
    except (Error, RuntimeError):
        raise Error("Maximum iterations reached in planet crossing calculation")


def helio_cross_ut(
    planet: int,
    x2cross: float,
    tjdut: float,
    flags: int = FLG_SWIEPH,
    backwards: bool | int | str = False,
) -> float:
    """
    Find when a planet crosses a specific heliocentric ecliptic longitude.

    Calculates the time when a planet, as seen from the Sun (heliocentric
    coordinates), crosses a specific ecliptic longitude. Useful for heliocentric
    astrology calculations.

    Searches forward (or backward if backwards=True) in time from tjdut.

    Args:
        planet: Planet ID (MERCURY, VENUS, EARTH, MARS, etc.)
                   Note: EARTH can be used to find when Earth crosses a longitude
                   as seen from the Sun.
        x2cross: Target heliocentric ecliptic longitude in degrees (0-360)
        tjdut: Julian Day (UT) to start search from
        flags: Calculation flags (FLG_SWIEPH, etc.). FLG_HELCTR is automatically added.
        backwards: If True, search backward in time instead of forward.
            Also accepts ints and the direction strings understood by
            :func:`libephemeris.eclipse._coerce_backwards` (e.g.
            ``"forward"``, ``"backward"``).

    Returns:
        float: Julian Day of crossing (UT)

    Raises:
        Error: If convergence fails or calculation error occurs

    Note:
        Heliocentric positions show where planets are relative to the Sun,
        not Earth. This is useful for:
        - Heliocentric astrology systems
        - Planetary synodic cycles
        - Finding planetary positions in solar-centered coordinates

    Algorithm:
        1. Get current heliocentric position and velocity
        2. Linear estimate: dt = (target - current) / velocity
        3. Refine with Newton-Raphson: jd_new = jd + (target - actual) / velocity
        4. Converge to < 0.001 arcsecond

    Example:
        >>> # Find when Mars crosses 0° heliocentric longitude
        >>> jd_cross = helio_cross_ut(MARS, 0.0, jd_now)
        >>> # Find when Earth crosses 90° as seen from the Sun
        >>> jd_earth_cross = helio_cross_ut(EARTH, 90.0, jd_now)
    """
    backwards = _coerce_backwards(backwards)
    x2cross, flags = _normalize_cross_target(x2cross, flags)

    # Always add FLG_HELCTR for heliocentric calculations
    helio_flag = flags | FLG_HELCTR | FLG_SPEED

    try:
        pos, _ = calc_ut(tjdut, planet, helio_flag)
        lon_start = pos[0]
        speed = pos[3]
    except (EphemerisRangeError, CalculationError, ValueError) as e:
        raise Error(f"Failed to calculate heliocentric planet position: {e}") from e

    # Estimate typical heliocentric speed if near zero
    # Heliocentric speeds are different from geocentric due to no retrograde
    typical_speeds = {
        2: 4.09,  # Mercury (heliocentric, no retro)
        3: 1.60,  # Venus
        14: 0.986,  # Earth
        4: 0.524,  # Mars
        5: 0.083,  # Jupiter
        6: 0.034,  # Saturn
        7: 0.012,  # Uranus
        8: 0.006,  # Neptune
        9: 0.004,  # Pluto
    }
    speed_default = typical_speeds.get(planet, 0.5)

    # Calculate initial guess
    diff = (x2cross - lon_start) % 360.0

    # Heliocentric planets don't go retrograde (except for very minor perturbations)
    # so we search forward by default, or backward if requested
    if backwards:
        # For backward search, map diff from [0, 360) into (-360, 0] so the
        # guess points into the past. A target essentially at the current
        # longitude resolves to a full cycle back (the nearest strictly-past
        # crossing).
        if diff > 1e-5:
            diff -= 360.0
        elif diff > -1e-5:
            diff -= 360.0
        # Size the backward guess with the MEAN heliocentric rate, not the
        # instantaneous speed. Near aphelion the instantaneous speed can run
        # ~30% below the mean (Mercury ~2.75 vs ~4.09 °/day), which places the
        # guess a half-revolution too far back and lets Newton converge into the
        # wrong basin — the crossing a whole revolution earlier. A post-
        # convergence check (_nearest_past_helio_crossing) still guards the
        # result, but the mean-rate guess keeps Newton in the right basin.
        guess_speed = speed_default if speed_default > 1e-9 else 1.0
        dt_guess = diff / guess_speed  # diff <= 0, guess_speed > 0 -> negative dt
    else:
        # No forward dead-band: a body a hair before (or exactly at) the target
        # legitimately crosses now, and the reference returns that imminent
        # crossing — not one orbital period later. This mirrors solcross/
        # mooncross/cross_ut, which have no dead-band either. (x2cross - lon)
        # % 360 already maps a just-passed target to ~360, i.e. the next cycle.
        if abs(speed) < 0.0001:
            speed = speed_default
        dt_guess = diff / speed
    jd_guess = tjdut + dt_guess

    # Adaptive iteration count based on planet speed
    max_iter = _get_adaptive_max_iterations(speed)

    # Helper function for Brent's method fallback (heliocentric)
    def get_helio_position(jd_time: float) -> Tuple[float, float]:
        pos_result, _ = calc_ut(jd_time, planet, helio_flag)
        return pos_result[0], pos_result[3]

    # Check if we're dealing with a very slow planet - use Brent's method for robustness
    # Heliocentric planets don't have true stations but very slow planets can still
    # benefit from the more robust method
    if _is_near_station(speed):
        try:
            search_window = (
                max(60.0, abs(dt_guess) * 2.0, abs(diff / speed_default) * 1.5)
                if speed_default > 0
                else max(120.0, abs(dt_guess) * 2.0)
            )
            if backwards:
                jd_bracket_start = tjdut - search_window
                jd_bracket_end = tjdut
            else:
                jd_bracket_start = tjdut
                jd_bracket_end = tjdut + search_window
            bracket_samples = max(40, int(search_window / 30))

            jd_a, jd_b = _find_bracket_for_crossing(
                get_helio_position,
                x2cross,
                jd_bracket_start,
                jd_bracket_end,
                num_samples=bracket_samples,
                # Backward search must return the NEAREST past crossing (the
                # one adjacent to tjdut), so scan the bracket from the tjdut end.
                from_end=backwards,
            )

            return _brent_find_crossing(
                get_helio_position, x2cross, jd_a, jd_b, NR_TOLERANCE, max_iter
            )
        except (Error, RuntimeError):
            pass  # Fall through to Newton-Raphson

    jd = jd_guess
    for iteration in range(max_iter):
        try:
            pos, _ = calc_ut(jd, planet, helio_flag)
            lon = pos[0]
            speed = pos[3]
        except (EphemerisRangeError, CalculationError, ValueError) as e:
            raise Error(
                f"Failed to calculate heliocentric position during iteration: {e}"
            ) from e

        diff = (x2cross - lon) % 360.0
        if diff > 180:
            diff -= 360

        # Check convergence (< 0.001 arcsecond)
        if abs(diff) < NR_TOLERANCE:
            if backwards:
                # Guarantee the nearest strictly-past crossing (Newton may have
                # landed a whole revolution too far back near aphelion).
                return _nearest_past_helio_crossing(
                    get_helio_position,
                    x2cross,
                    tjdut,
                    jd,
                    360.0 / speed_default if speed_default > 1e-9 else 720.0,
                    NR_TOLERANCE,
                )
            return jd

        # Update max_iter based on current speed (may change during iteration)
        max_iter = max(max_iter, _get_adaptive_max_iterations(speed))

        if abs(speed) < 0.0001:
            speed = speed_default

        jd += diff / speed

        # Safety: scale max_range with dt_guess so slow planets have enough room.
        # Heliocentric orbits are smooth (no retrograde), so NR converges well —
        # but the initial guess may be off by a significant fraction of dt_guess
        # due to elliptical orbit speed variation.
        base_range = 500.0 if abs(speed_default) < 0.05 else 400.0
        max_range = max(base_range, abs(dt_guess) * 2.0)
        if abs(jd - tjdut) > max_range:
            raise Error("Heliocentric crossing search diverged")

    raise Error("Maximum iterations reached in heliocentric crossing calculation")


def helio_cross(
    planet: int,
    x2cross: float,
    tjdet: float,
    flags: int = FLG_SWIEPH,
    backwards: bool | int | str = False,
) -> float:
    """
    Find when a planet crosses a specific heliocentric ecliptic longitude (TT version).

    This is the Terrestrial Time version of helio_cross_ut(). Takes Julian Day
    in TT (Terrestrial Time, also known as Ephemeris Time) instead of UT.

    Calculates the time when a planet, as seen from the Sun (heliocentric
    coordinates), crosses a specific ecliptic longitude.

    Searches forward (or backward if backwards=True) in time from tjdet.

    Args:
        planet: Planet ID (MERCURY, VENUS, EARTH, MARS, etc.)
        x2cross: Target heliocentric ecliptic longitude in degrees (0-360)
        tjdet: Julian Day in Terrestrial Time (TT/ET) to start search from
        flags: Calculation flags (FLG_SWIEPH, etc.). FLG_HELCTR is automatically added.
        backwards: If True, search backward in time instead of forward.
            Also accepts ints and the direction strings understood by
            :func:`libephemeris.eclipse._coerce_backwards` (e.g.
            ``"forward"``, ``"backward"``).

    Returns:
        float: Julian Day of crossing (TT)

    Raises:
        Error: If convergence fails or calculation error occurs

    Note:
        TT (Terrestrial Time) differs from UT (Universal Time) by Delta T,
        which varies from ~32 seconds (year 2000) to minutes (historical times).
        For most astrological applications, use helio_cross_ut() instead.

    Algorithm:
        1. Get current heliocentric position and velocity
        2. Linear estimate: dt = (target - current) / velocity
        3. Refine with Newton-Raphson: jd_new = jd + (target - actual) / velocity
        4. Converge to < 0.001 arcsecond

    Example:
        >>> # Find when Mars crosses 0° heliocentric longitude using TT
        >>> jd_cross_tt = helio_cross(MARS, 0.0, jd_tt_now)
    """
    backwards = _coerce_backwards(backwards)
    x2cross, flags = _normalize_cross_target(x2cross, flags)

    # Always add FLG_HELCTR for heliocentric calculations
    helio_flag = flags | FLG_HELCTR | FLG_SPEED

    try:
        pos, _ = calc(tjdet, planet, helio_flag)
        lon_start = pos[0]
        speed = pos[3]
    except (EphemerisRangeError, CalculationError, ValueError) as e:
        raise Error(f"Failed to calculate heliocentric planet position: {e}") from e

    # Estimate typical heliocentric speed if near zero
    typical_speeds = {
        2: 4.09,  # Mercury
        3: 1.60,  # Venus
        14: 0.986,  # Earth
        4: 0.524,  # Mars
        5: 0.083,  # Jupiter
        6: 0.034,  # Saturn
        7: 0.012,  # Uranus
        8: 0.006,  # Neptune
        9: 0.004,  # Pluto
    }
    speed_default = typical_speeds.get(planet, 0.5)

    # Calculate initial guess
    diff = (x2cross - lon_start) % 360.0

    if backwards:
        # For backward search, map diff from [0, 360) into (-360, 0] so the
        # guess points into the past. A target essentially at the current
        # longitude resolves to a full cycle back (the nearest strictly-past
        # crossing).
        if diff > 1e-5:
            diff -= 360.0
        elif diff > -1e-5:
            diff -= 360.0
        # Size the backward guess with the MEAN heliocentric rate, not the
        # instantaneous speed (which near aphelion runs well below the mean and
        # would place the guess a half-revolution too far back). See the UT twin
        # helio_cross_ut for the full rationale.
        guess_speed = speed_default if speed_default > 1e-9 else 1.0
        dt_guess = diff / guess_speed  # diff <= 0, guess_speed > 0 -> negative dt
    else:
        # No forward dead-band: a body a hair before (or exactly at) the target
        # legitimately crosses now, and the reference returns that imminent
        # crossing — not one orbital period later. This mirrors solcross/
        # mooncross/cross_ut, which have no dead-band either. (x2cross - lon)
        # % 360 already maps a just-passed target to ~360, i.e. the next cycle.
        if abs(speed) < 0.0001:
            speed = speed_default
        dt_guess = diff / speed

    jd_guess = tjdet + dt_guess

    # Adaptive iteration count based on planet speed
    max_iter = _get_adaptive_max_iterations(speed)

    # Helper function for Brent's method fallback (heliocentric TT)
    def get_helio_position_tt(jd_time: float) -> Tuple[float, float]:
        pos_result, _ = calc(jd_time, planet, helio_flag)
        return pos_result[0], pos_result[3]

    # Check if we're dealing with a very slow planet - use Brent's method for robustness
    if _is_near_station(speed):
        try:
            search_window = (
                max(60.0, abs(dt_guess) * 2.0, abs(diff / speed_default) * 1.5)
                if speed_default > 0
                else max(120.0, abs(dt_guess) * 2.0)
            )
            if backwards:
                jd_bracket_start = tjdet - search_window
                jd_bracket_end = tjdet
            else:
                jd_bracket_start = tjdet
                jd_bracket_end = tjdet + search_window
            bracket_samples = max(40, int(search_window / 30))

            jd_a, jd_b = _find_bracket_for_crossing(
                get_helio_position_tt,
                x2cross,
                jd_bracket_start,
                jd_bracket_end,
                num_samples=bracket_samples,
                # Backward search must return the NEAREST past crossing (the
                # one adjacent to tjdet), so scan the bracket from the tjdet end.
                from_end=backwards,
            )

            return _brent_find_crossing(
                get_helio_position_tt, x2cross, jd_a, jd_b, NR_TOLERANCE, max_iter
            )
        except (Error, RuntimeError):
            pass  # Fall through to Newton-Raphson

    jd = jd_guess
    for iteration in range(max_iter):
        try:
            pos, _ = calc(jd, planet, helio_flag)
            lon = pos[0]
            speed = pos[3]
        except (EphemerisRangeError, CalculationError, ValueError) as e:
            raise Error(
                f"Failed to calculate heliocentric position during iteration: {e}"
            ) from e

        diff = (x2cross - lon) % 360.0
        if diff > 180:
            diff -= 360

        # Check convergence (< 0.001 arcsecond)
        if abs(diff) < NR_TOLERANCE:
            if backwards:
                # Guarantee the nearest strictly-past crossing (Newton may have
                # landed a whole revolution too far back near aphelion).
                return _nearest_past_helio_crossing(
                    get_helio_position_tt,
                    x2cross,
                    tjdet,
                    jd,
                    360.0 / speed_default if speed_default > 1e-9 else 720.0,
                    NR_TOLERANCE,
                )
            return jd

        # Update max_iter based on current speed (may change during iteration)
        max_iter = max(max_iter, _get_adaptive_max_iterations(speed))

        if abs(speed) < 0.0001:
            speed = speed_default

        jd += diff / speed

        # Safety: scale max_range with dt_guess so slow planets have enough room
        base_range = 500.0 if abs(speed_default) < 0.05 else 400.0
        max_range = max(base_range, abs(dt_guess) * 2.0)
        if abs(jd - tjdet) > max_range:
            raise Error("Heliocentric crossing search diverged")

    raise Error("Maximum iterations reached in heliocentric crossing calculation")


# =============================================================================
# RETROGRADE STATION HANDLING
# =============================================================================
#
# Retrograde stations occur when a planet appears to stop and reverse direction
# as seen from Earth. At station points, the planet's velocity approaches zero,
# which can cause numerical instabilities in calculations.
#
# Two types of stations:
#   - Retrograde Station (SR): Planet slows, stops, and begins retrograde motion
#   - Direct Station (SD): Planet stops retrograde motion and resumes direct motion
#
# The functions below provide robust handling for calculations near these points.
# =============================================================================

# Tolerance for station finding (velocity threshold)
# A planet is considered "stationary" when |velocity| < this value
STATION_VELOCITY_TOLERANCE = 1e-6  # degrees/day


def is_retrograde(planet_id: int, jd_ut: float, flag: int = FLG_SWIEPH) -> bool:
    """
    Check if a planet is currently in retrograde motion.

    A planet is retrograde when its geocentric longitude is decreasing
    (negative velocity in longitude).

    Args:
        planet_id: Planet ID (MERCURY, VENUS, MARS, etc.)
        jd_ut: Julian Day (UT) to check
        flag: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        bool: True if planet is retrograde (velocity < 0), False otherwise

    Note:
        Sun and Moon never go retrograde from geocentric perspective.
        Only planets beyond Earth (Mars-Pluto) and Mercury/Venus can appear retrograde.

    Example:
        >>> is_retrograde(MERCURY, jd_now)
        True  # Mercury is currently retrograde
    """
    if planet_id in (SUN, MOON):
        # Sun and Moon are never retrograde from geocentric view
        return False

    try:
        pos, _ = calc_ut(jd_ut, planet_id, flag | FLG_SPEED)
        return pos[3] < 0
    except (EphemerisRangeError, CalculationError, ValueError):
        return False


def get_station_type(planet_id: int, jd_ut: float, flag: int = FLG_SWIEPH) -> str:
    """
    Determine if a planet is near a station and what type.

    Args:
        planet_id: Planet ID (MERCURY, VENUS, MARS, etc.)
        jd_ut: Julian Day (UT) to check
        flag: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        str: One of:
            - "direct": Planet is in direct (prograde) motion
            - "retrograde": Planet is in retrograde motion
            - "stationary_retrograde": Near station, about to go retrograde
            - "stationary_direct": Near station, about to resume direct motion

    Example:
        >>> get_station_type(MARS, jd_station)
        'stationary_retrograde'
    """
    if planet_id in (SUN, MOON):
        return "direct"

    try:
        pos, _ = calc_ut(jd_ut, planet_id, flag | FLG_SPEED)
        speed = pos[3]

        if _is_near_station(speed):
            # Near station - determine type by checking velocity trend
            # Look at velocity slightly before and after
            dt = 1.0  # 1 day
            pos_before, _ = calc_ut(jd_ut - dt, planet_id, flag | FLG_SPEED)
            pos_after, _ = calc_ut(jd_ut + dt, planet_id, flag | FLG_SPEED)

            speed_before = pos_before[3]
            speed_after = pos_after[3]

            # If speed is decreasing (becoming more negative), heading to retrograde
            if speed_before > speed_after:
                return "stationary_retrograde"
            else:
                return "stationary_direct"
        elif speed < 0:
            return "retrograde"
        else:
            return "direct"

    except (EphemerisRangeError, CalculationError, ValueError):
        return "direct"


def _brent_find_station(
    get_speed_func,
    jd_a: float,
    jd_b: float,
    tolerance: float = STATION_VELOCITY_TOLERANCE,
    max_iter: int = 100,
) -> float:
    """
    Find exact station time using Brent's method for root finding.

    Finds when velocity = 0 (station point) within the bracket [jd_a, jd_b].

    Args:
        get_speed_func: Function(jd) -> velocity in degrees/day
        jd_a: Start of bracket (Julian Day)
        jd_b: End of bracket (Julian Day)
        tolerance: Convergence tolerance for velocity
        max_iter: Maximum iterations

    Returns:
        float: Julian Day when velocity = 0

    Raises:
        Error: If root not bracketed or convergence fails
    """
    fa = get_speed_func(jd_a)
    fb = get_speed_func(jd_b)

    # Check if root is bracketed (velocity changes sign)
    if fa * fb > 0:
        raise Error(f"Station not bracketed: speed(a)={fa:.6f}, speed(b)={fb:.6f}")

    # Ensure |f(a)| >= |f(b)|
    if abs(fa) < abs(fb):
        jd_a, jd_b = jd_b, jd_a
        fa, fb = fb, fa

    c = jd_a
    fc = fa
    mflag = True
    d = 0.0

    for _ in range(max_iter):
        if abs(fb) < tolerance:
            return jd_b

        if fa != fc and fb != fc:
            # Inverse quadratic interpolation
            s = (
                jd_a * fb * fc / ((fa - fb) * (fa - fc))
                + jd_b * fa * fc / ((fb - fa) * (fb - fc))
                + c * fa * fb / ((fc - fa) * (fc - fb))
            )
        else:
            # Secant method
            s = jd_b - fb * (jd_b - jd_a) / (fb - fa)

        # Conditions for bisection fallback
        cond1 = not (
            (3 * jd_a + jd_b) / 4 < s < jd_b or jd_b < s < (3 * jd_a + jd_b) / 4
        )
        cond2 = mflag and abs(s - jd_b) >= abs(jd_b - c) / 2
        cond3 = not mflag and abs(s - jd_b) >= abs(c - d) / 2
        cond4 = mflag and abs(jd_b - c) < _BRENT_TIME_EPS_DAYS
        cond5 = not mflag and abs(c - d) < _BRENT_TIME_EPS_DAYS

        if cond1 or cond2 or cond3 or cond4 or cond5:
            s = (jd_a + jd_b) / 2
            mflag = True
        else:
            mflag = False

        fs = get_speed_func(s)
        d = c
        c = jd_b
        fc = fb

        if fa * fs < 0:
            jd_b = s
            fb = fs
        else:
            jd_a = s
            fa = fs

        if abs(fa) < abs(fb):
            jd_a, jd_b = jd_b, jd_a
            fa, fb = fb, fa

    return jd_b


def _find_station_bracket(
    planet_id: int,
    jd_start: float,
    jd_end: float,
    flag: int,
    num_samples: int = 50,
) -> Tuple[float, float]:
    """
    Find a bracket [jd_a, jd_b] containing a station (velocity sign change).

    Args:
        planet_id: Planet ID
        jd_start: Start of search interval
        jd_end: End of search interval
        flag: Calculation flags
        num_samples: Number of samples to take

    Returns:
        Tuple[float, float]: (jd_a, jd_b) bracket containing a station

    Raises:
        Error: If no station found in interval
    """
    step = (jd_end - jd_start) / num_samples

    pos_prev, _ = calc_ut(jd_start, planet_id, flag | FLG_SPEED)
    speed_prev = pos_prev[3]
    jd_prev = jd_start

    for i in range(1, num_samples + 1):
        jd_curr = jd_start + i * step
        pos_curr, _ = calc_ut(jd_curr, planet_id, flag | FLG_SPEED)
        speed_curr = pos_curr[3]

        # Check for sign change (station)
        if speed_prev * speed_curr <= 0:
            return (jd_prev, jd_curr)

        jd_prev = jd_curr
        speed_prev = speed_curr

    raise Error(f"No station found for planet {planet_id} in [{jd_start}, {jd_end}]")


def find_station_ut(
    planet_id: int,
    jd_ut: float,
    station_type: str = "any",
    flag: int = FLG_SWIEPH,
) -> Tuple[float, str]:
    """
    Find the next planetary station (stationary point) after a given date.

    A station occurs when a planet's apparent velocity reaches zero, marking
    the transition between direct and retrograde motion (or vice versa).

    Args:
        planet_id: Planet ID (MERCURY through PLUTO)
        jd_ut: Julian Day (UT) to start search from
        station_type: Type of station to find:
            - "any": Find next station regardless of type (default)
            - "retrograde" or "SR": Find next station where planet goes retrograde
            - "direct" or "SD": Find next station where planet resumes direct motion
        flag: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple[float, str]: (Julian Day of station, station type "SR" or "SD")

    Raises:
        ValueError: If planet_id is Sun or Moon (never station)
        Error: If convergence fails or no station found

    Note:
        Sun and Moon never go retrograde from Earth's perspective.

    Station Periods (approximate time between stations):
        - Mercury: ~116 days (synodic period)
        - Venus: ~584 days
        - Mars: ~780 days
        - Jupiter: ~399 days
        - Saturn: ~378 days
        - Uranus: ~370 days
        - Neptune: ~367 days
        - Pluto: ~367 days

    Precision:
        Typically finds station to within ~1 minute of arc in velocity

    Example:
        >>> jd_station, stype = find_station_ut(MERCURY, jd_now)
        >>> print(f"Mercury stations {stype} on JD {jd_station}")
    """
    if planet_id in (SUN, MOON):
        raise ValueError("Sun and Moon do not have retrograde stations")

    # Normalize station_type
    if station_type.upper() in ("SR", "RETROGRADE"):
        station_type = "retrograde"
    elif station_type.upper() in ("SD", "DIRECT"):
        station_type = "direct"
    else:
        station_type = "any"

    # Typical synodic periods (days) - approximate time between same-type stations
    synodic_periods = {
        2: 116,  # Mercury
        3: 584,  # Venus
        4: 780,  # Mars
        5: 399,  # Jupiter
        6: 378,  # Saturn
        7: 370,  # Uranus
        8: 367,  # Neptune
        9: 367,  # Pluto
    }

    # Search window is half synodic period (time between SR and SD)
    search_window = synodic_periods.get(planet_id, 400) / 2

    def get_speed(jd: float) -> float:
        pos, _ = calc_ut(jd, planet_id, flag | FLG_SPEED)
        return pos[3]

    # Get current motion direction
    get_speed(jd_ut)

    jd_search_start = jd_ut
    max_attempts = 4  # Allow searching up to 2 full synodic periods

    for attempt in range(max_attempts):
        try:
            # Find the next station bracket
            jd_search_end = jd_search_start + search_window
            jd_a, jd_b = _find_station_bracket(
                planet_id, jd_search_start, jd_search_end, flag
            )

            # Find exact station time
            jd_station = _brent_find_station(
                get_speed, jd_a, jd_b, STATION_VELOCITY_TOLERANCE
            )

            # Determine station type by looking at motion before and after
            speed_before = get_speed(jd_station - 1.0)
            speed_after = get_speed(jd_station + 1.0)

            if speed_before > 0 and speed_after < 0:
                found_type = "SR"  # Stationary Retrograde
            else:
                found_type = "SD"  # Stationary Direct

            # Check if this matches requested type
            if station_type == "any":
                return (jd_station, found_type)
            elif station_type == "retrograde" and found_type == "SR":
                return (jd_station, found_type)
            elif station_type == "direct" and found_type == "SD":
                return (jd_station, found_type)

            # Not the right type, search for next station
            jd_search_start = jd_station + 1.0

        except (Error, RuntimeError):
            # No station found in this window, extend search
            jd_search_start = jd_search_start + search_window

    raise Error(f"Could not find {station_type} station for planet {planet_id}")


def next_retrograde_ut(
    planet_id: int, jd_ut: float, flag: int = FLG_SWIEPH
) -> Tuple[float, float]:
    """
    Find the next retrograde period for a planet.

    Returns the start (SR) and end (SD) Julian Days for the next retrograde period.

    Args:
        planet_id: Planet ID (MERCURY through PLUTO)
        jd_ut: Julian Day (UT) to start search from
        flag: Calculation flags

    Returns:
        Tuple[float, float]: (JD of retrograde start, JD of retrograde end)

    Raises:
        ValueError: If planet is Sun or Moon
        Error: If stations cannot be found

    Example:
        >>> jd_sr, jd_sd = next_retrograde_ut(MERCURY, jd_now)
        >>> print(f"Mercury Rx: {jd_sr} to {jd_sd}")
    """
    # Check if currently retrograde
    if is_retrograde(planet_id, jd_ut, flag):
        # Find the end of current retrograde (SD), then find next SR
        jd_sd, _ = find_station_ut(planet_id, jd_ut, "direct", flag)
        jd_sr, _ = find_station_ut(planet_id, jd_sd + 1.0, "retrograde", flag)
        jd_sd_next, _ = find_station_ut(planet_id, jd_sr + 1.0, "direct", flag)
        return (jd_sr, jd_sd_next)
    else:
        # Find next SR and then the following SD
        jd_sr, _ = find_station_ut(planet_id, jd_ut, "retrograde", flag)
        jd_sd, _ = find_station_ut(planet_id, jd_sr + 1.0, "direct", flag)
        return (jd_sr, jd_sd)


def calc_velocity_at_station(
    planet_id: int,
    jd_station: float,
    flag: int = FLG_SWIEPH,
    dt: float = 1.0,
) -> Tuple[float, float, float]:
    """
    Calculate stable velocity components near a station point.

    Near stations, the standard 1-second numerical differentiation can have
    numerical noise issues. This function uses a wider timestep for more
    stable velocity estimation.

    Args:
        planet_id: Planet ID
        jd_station: Julian Day near the station
        flag: Calculation flags
        dt: Timestep in days for differentiation (default 1.0 day)

    Returns:
        Tuple[float, float, float]: (velocity_lon, velocity_lat, velocity_dist)
            in degrees/day, degrees/day, AU/day respectively

    Note:
        The velocity_lon will be very close to 0 at the actual station point.
        This function is useful for getting stable readings when the standard
        calc_ut velocity might show numerical noise.

    Example:
        >>> vlon, vlat, vdist = calc_velocity_at_station(MARS, jd_station)
        >>> print(f"Station velocity: {vlon:.6f} deg/day")  # Should be ~0
    """
    # Get positions at t-dt and t+dt
    pos_before, _ = calc_ut(jd_station - dt, planet_id, flag)
    pos_after, _ = calc_ut(jd_station + dt, planet_id, flag)

    # Calculate velocity using central difference
    lon_before, lat_before, dist_before = pos_before[0], pos_before[1], pos_before[2]
    lon_after, lat_after, dist_after = pos_after[0], pos_after[1], pos_after[2]

    # Handle longitude wraparound
    dlon = lon_after - lon_before
    if dlon > 180:
        dlon -= 360
    elif dlon < -180:
        dlon += 360

    v_lon = dlon / (2 * dt)
    v_lat = (lat_after - lat_before) / (2 * dt)
    v_dist = (dist_after - dist_before) / (2 * dt)

    return (v_lon, v_lat, v_dist)


def _find_previous_station_ut(
    planet_id: int,
    jd_ut: float,
    flag: int = FLG_SWIEPH,
) -> Tuple[float, str]:
    """
    Find the most recent planetary station strictly before ``jd_ut``.

    The backward counterpart of :func:`find_station_ut`: it scans backward from
    ``jd_ut`` and returns the nearest PAST station (velocity sign change),
    together with its type. Used by :func:`get_station_info` so a just-passed
    station can be reported with a negative ``days_to_station``.

    Args:
        planet_id: Planet ID (MERCURY through PLUTO)
        jd_ut: Julian Day (UT) to search backward from
        flag: Calculation flags

    Returns:
        Tuple[float, str]: (Julian Day of station, station type "SR" or "SD")

    Raises:
        ValueError: If planet_id is Sun or Moon (never station)
        Error: If no station is found within the search window
    """
    if planet_id in (SUN, MOON):
        raise ValueError("Sun and Moon do not have retrograde stations")

    # Same synodic periods as find_station_ut. The gap to the previous station
    # is at most a bit under one synodic period, so scan back a full period.
    synodic_periods = {
        2: 116,  # Mercury
        3: 584,  # Venus
        4: 780,  # Mars
        5: 399,  # Jupiter
        6: 378,  # Saturn
        7: 370,  # Uranus
        8: 367,  # Neptune
        9: 367,  # Pluto
    }
    search_window = float(synodic_periods.get(planet_id, 400))

    def get_speed(jd: float) -> float:
        pos, _ = calc_ut(jd, planet_id, flag | FLG_SPEED)
        return pos[3]

    # ~2-day resolution keeps us well below the shortest retrograde arc
    # (Mercury's, ~21 days) so no station pair is stepped over.
    num_samples = max(60, int(search_window / 2))
    step = search_window / num_samples

    jd_next = jd_ut
    speed_next = get_speed(jd_ut)

    for i in range(1, num_samples + 1):
        jd_curr = jd_ut - i * step
        speed_curr = get_speed(jd_curr)

        # First velocity sign change scanning backward = nearest past station.
        if speed_next * speed_curr <= 0:
            jd_station = _brent_find_station(
                get_speed, jd_curr, jd_next, STATION_VELOCITY_TOLERANCE
            )

            speed_before = get_speed(jd_station - 1.0)
            speed_after = get_speed(jd_station + 1.0)
            if speed_before > 0 and speed_after < 0:
                found_type = "SR"  # Stationary Retrograde
            else:
                found_type = "SD"  # Stationary Direct

            return (jd_station, found_type)

        jd_next = jd_curr
        speed_next = speed_curr

    raise Error(f"No previous station found for planet {planet_id}")


def get_station_info(planet_id: int, jd_ut: float, flag: int = FLG_SWIEPH) -> dict:
    """
    Get comprehensive information about the nearest station.

    Finds the nearest station to the given date and returns detailed
    information about it.

    Args:
        planet_id: Planet ID (MERCURY through PLUTO)
        jd_ut: Julian Day (UT)
        flag: Calculation flags

    Returns:
        dict: Station information containing:
            - "jd_station": Julian Day of nearest station
            - "station_type": "SR" (stationary retrograde) or "SD" (stationary direct)
            - "days_to_station": Days until station (negative if past)
            - "longitude_at_station": Ecliptic longitude at station
            - "is_currently_retrograde": Current retrograde status
            - "velocity": Current velocity in deg/day

    Example:
        >>> info = get_station_info(MERCURY, jd_now)
        >>> print(f"Next station: {info['station_type']} in {info['days_to_station']:.1f} days")
    """
    if planet_id in (SUN, MOON):
        raise ValueError("Sun and Moon do not have stations")

    # Get current state
    pos, _ = calc_ut(jd_ut, planet_id, flag | FLG_SPEED)
    current_velocity = pos[3]
    current_retrograde = current_velocity < 0

    # Find the NEAREST station in either direction. find_station_ut searches
    # strictly forward, so on its own it never reports a just-passed station;
    # pair it with the backward search and keep whichever is closer in time.
    candidates = []  # (jd_station, station_type)
    forward_error = None
    try:
        candidates.append(find_station_ut(planet_id, jd_ut, "any", flag))
    except (Error, RuntimeError) as e:
        forward_error = e
    try:
        candidates.append(_find_previous_station_ut(planet_id, jd_ut, flag))
    except (Error, RuntimeError) as e:
        forward_error = forward_error or e

    if not candidates:
        # Return partial info if both searches fail
        return {
            "jd_station": None,
            "station_type": None,
            "days_to_station": None,
            "longitude_at_station": None,
            "is_currently_retrograde": current_retrograde,
            "velocity": current_velocity,
            "error": str(forward_error),
        }

    # Nearest by absolute time distance. days_to_station is signed: negative
    # when the nearest station is in the past.
    jd_station, stype = min(candidates, key=lambda c: abs(c[0] - jd_ut))
    pos_station, _ = calc_ut(jd_station, planet_id, flag)

    return {
        "jd_station": jd_station,
        "station_type": stype,
        "days_to_station": jd_station - jd_ut,
        "longitude_at_station": pos_station[0],
        "is_currently_retrograde": current_retrograde,
        "velocity": current_velocity,
    }
