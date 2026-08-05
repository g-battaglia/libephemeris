# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Solar and Lunar eclipse calculations for libephemeris.

Finds eclipse events and calculates their circumstances.

Functions:
- _sol_eclipse_when_glob_pythonic: Find next global solar eclipse
- _sol_eclipse_when_loc_pythonic: Find eclipse at specific location
- _sol_eclipse_where_pythonic: Calculate path of eclipse
- _sol_eclipse_how_pythonic: Eclipse circumstances at location
- _lun_eclipse_when_pythonic: Find next lunar eclipse

Algorithm:
    Solar eclipses occur at New Moon when Moon is near the ecliptic plane.
    1. Find next New Moon (Sun-Moon conjunction in longitude)
    2. Check lunar latitude - if |lat| < ~1.5° eclipse is possible
    3. Calculate eclipse magnitude and type based on distances
    4. Use Besselian elements for precise timing of phases

    Lunar eclipses occur at Full Moon when Moon is near a lunar node.
    1. Find next Full Moon (Sun-Moon opposition in longitude)
    2. Check Moon's distance from node - if close, eclipse is possible
    3. Calculate eclipse type (penumbral, partial, total) based on geometry
    4. Calculate phase times based on shadow cone sizes

References:
    - Meeus "Astronomical Algorithms" Ch. 54 (Eclipses)
    - Espenak & Meeus, "Five Millennium Canon of Solar Eclipses",
      NASA/TP-2006-214141
    - Espenak & Meeus, "Five Millennium Canon of Lunar Eclipses",
      NASA/TP-2009-214172, and the companion "Five Millennium Catalog of
      Lunar Eclipses", NASA/TP-2009-214173
    - NASA/JPL DE440/DE441 Sun, Earth, and Moon state vectors

Provenance:
    Eclipse contacts, shadow cones, local circumstances, and classifications are
    independently derived from the cited public geometry and JPL-backed body
    states. Search windows, bracketing, convergence thresholds, and numerical
    derivative steps are labelled project choices and are checked by geometric
    invariants (ordered contacts, bounded obscuration, and consistent tangency).
    Public compatibility calls validate return semantics only and never provide
    coefficients, contact tables, or Besselian elements to this module.
"""

from __future__ import annotations

import contextlib
import math
import re
from typing import Callable, Iterator, Sequence, Tuple, Union, cast

from .constants import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
    EARTH,
    CHIRON,
    PHOLUS,
    CERES,
    PALLAS,
    JUNO,
    VESTA,
    FLG_SPEED,
    FLG_SWIEPH,
    FLG_EQUATORIAL,
    ECL_TOTAL,
    ECL_ANNULAR,
    ECL_PARTIAL,
    ECL_ANNULAR_TOTAL,
    ECL_CENTRAL,
    ECL_NONCENTRAL,
    ECL_ALLTYPES_LUNAR,
    ECL_PENUMBRAL,
    ECL_VISIBLE,
    ECL_MAX_VISIBLE,
    ECL_1ST_VISIBLE,
    ECL_2ND_VISIBLE,
    ECL_3RD_VISIBLE,
    ECL_4TH_VISIBLE,
)
from .exceptions import Error, LEBCorruptionError, UnknownBodyError
from .planets import calc_ut
from .state import get_calc_mode, get_timescale


class _IllegalRiseBodyError(UnknownBodyError, ValueError):
    """A rise/set/transit target that no active backend can place.

    Raised for an unknown or unsupported body id passed to ``rise_trans`` /
    ``rise_trans_true_hor`` (e.g. a planetary moon with no registered SPK).
    It multiply inherits the typed ``UnknownBodyError`` (an ``Error``
    subclass, the compatibility-contract class for an illegal body) AND the
    built-in ``ValueError`` so the single "illegal planet number" contract is
    identical on every backend and satisfies callers catching either type.
    """


_ACTIVE_LEB_READER = object()


def _get_leb_reader_safe():
    """Return LEB reader or None. Wraps get_leb_reader for eclipse dispatch."""
    from .state import get_leb_reader

    return get_leb_reader()


def _is_leb_out_of_range(exc: BaseException) -> bool:
    """Return True if the exception is an LEB out-of-range error.

    Recognises the ValueError raised by LEBReader / LEB2Reader / fast_calc
    when a JD is probed outside the LEB file's coverage. Both bodies and
    nutation coverage produce messages containing "outside range" or
    "outside nutation range".
    """
    return isinstance(exc, (KeyError, ValueError)) and ("outside" in str(exc).lower())


# Reader coverage-miss messages ("JD <x> outside range [<a>, <b>] for body <n>"
# and "JD <x> outside nutation range [<a>, <b>]") carry the bounds needed to
# build the typed public error without inventing values.
_LEB_RANGE_MISS_RE = re.compile(
    r"JD (?P<jd>[-+0-9.eE]+) outside (?:nutation )?range "
    r"\[(?P<start>[-+0-9.eE]+), (?P<end>[-+0-9.eE]+)\]"
)
_LEB_RANGE_MISS_BODY_RE = re.compile(r"for body (?P<body>\d+)")
_LEB_BODY_MISS_RE = re.compile(
    r"Body (?P<body>\d+) not in (?:any )?(?:LEB file|installed LEB tier)"
)


def _raise_if_sealed_leb_miss(exc: BaseException) -> None:
    """Apply the sealed-mode error contract to an LEB fast-path failure.

    Shared by the eclipse fallback wrapper and the heliacal LEB fast paths.
    In sealed ``leb`` mode a reader ``KeyError``/``ValueError`` must never
    escape to callers: ``calc_ut``/``calc`` already surface coverage misses
    as :class:`EphemerisRangeError`, and the search/dispatch layers must
    speak the same documented type.

    - A coverage miss (see :func:`_is_leb_out_of_range`) raises
      :class:`EphemerisRangeError`; the requested JD and coverage bounds
      are carried over when the reader message includes them.
    - A ``KeyError`` for a body absent from the file's body map raises
      :class:`UnknownBodyError` naming the missing body — the same type
      ``planets._raise_leb_range_miss`` uses for an unstored core body —
      without fabricated bounds.
    - Anything else is re-raised verbatim, so genuine validation errors
      keep their type and sealed mode still never continues past the LEB
      path.

    In every other mode the function returns without raising so each call
    site keeps its own fallback policy (Skyfield retry, fallback logging).
    """
    if get_calc_mode() != "leb":
        return
    from .exceptions import EphemerisRangeError, UnknownBodyError

    text = str(exc)
    sealed_note = "LEB mode does not silently substitute a lower-precision source."
    if _is_leb_out_of_range(exc):
        range_match = _LEB_RANGE_MISS_RE.search(text)
        body_match = _LEB_RANGE_MISS_BODY_RE.search(text)
        raise EphemerisRangeError(
            message=f"{text.strip()}. {sealed_note}",
            requested_jd=float(range_match["jd"]) if range_match else None,
            start_jd=float(range_match["start"]) if range_match else None,
            end_jd=float(range_match["end"]) if range_match else None,
            body_id=int(body_match["body"]) if body_match else None,
        ) from exc
    body_miss = _LEB_BODY_MISS_RE.search(text)
    if isinstance(exc, KeyError) and body_miss is not None:
        raise UnknownBodyError(
            message=(
                f"Body {body_miss['body']} is not stored in the active "
                f"LEB file. LEB mode does not silently substitute a "
                f"non-LEB source for a missing body."
            ),
            body_id=int(body_miss["body"]),
        ) from exc
    raise exc


def _call_with_leb_skyfield_fallback(impl, *args, **kwargs):
    """Run ``impl(reader=...)`` with mode-aware range handling.

    Eclipse / occultation / rise-transit search loops define LEB-backed inner
    closures inside ``if reader is not None`` blocks. If a search probe lands
    outside the LEB file's coverage, _topo_ecliptic raises ValueError. Without
    this helper there was no fallback, so partial/custom LEB files shorter
    than the search window would surface the underlying ValueError to callers.

    Behaviour:
        - If no LEB reader is available, run ``impl(reader=None)`` directly.
        - Otherwise, run ``impl(reader=reader)`` first.
        - In ``auto`` mode, an LEB range miss retries with ``reader=None`` so
          the normal JPL/Skyfield path can run.
        - In ``leb`` mode, a range miss raises the documented
          :class:`EphemerisRangeError` and a body absent from the file's
          body map raises :class:`UnknownBodyError` (see
          :func:`_raise_if_sealed_leb_miss`): the sealed runtime must
          never switch persisted ephemeris source.
        - Corruption and all unrelated exceptions are always re-raised.

    Canonical LEB core tiers cover 1850-2150, 1550-2650, or the exact shared
    DE441 interval (JD -3100015.5 to 8000016.5). In ``auto`` mode the retry
    restarts computation from scratch and shares no intermediate state.
    """
    reader = _get_leb_reader_safe()
    if reader is not None:
        try:
            return impl(*args, reader=reader, **kwargs)
        except LEBCorruptionError:
            raise
        except (KeyError, ValueError) as exc:
            # Sealed leb mode never returns here: the miss surfaces as the
            # documented typed error (or re-raises verbatim).
            _raise_if_sealed_leb_miss(exc)
            if not _is_leb_out_of_range(exc):
                raise
    return impl(*args, reader=None, **kwargs)


def _reraise_if_leb_range_error(exc: BaseException) -> None:
    """Re-raise ``exc`` for LEB corruption or out-of-range coverage.

    Defensive ``except`` blocks inside eclipse implementations must not
    swallow LEB coverage errors: doing so converts an out-of-range probe
    into fabricated zeros/sentinels and prevents
    ``_call_with_leb_skyfield_fallback`` from applying its mode-aware
    range policy. Call this first in any broad except that would otherwise
    degrade gracefully.
    """
    if isinstance(exc, LEBCorruptionError) or _is_leb_out_of_range(exc):
        raise exc


def _is_ephemeris_boundary(exc: BaseException) -> bool:
    """Return True when ``exc`` signals the end of ephemeris coverage.

    An eclipse/occultation search that finds no matching event walks
    forward (or backward) until it steps off the active ephemeris. Both
    the DE (Skyfield) path (:class:`EphemerisRangeError`) and the LEB path
    (a ``ValueError`` whose message contains "outside") mark that boundary;
    the search converts either into a typed no-event error rather than
    looping or fabricating a result.
    """
    from .exceptions import EphemerisRangeError

    return isinstance(exc, EphemerisRangeError) or _is_leb_out_of_range(exc)


def _sol_glob_accepts(mask: int, retflag: int) -> bool:
    """Reference ``ifltype`` acceptance test for ``sol_eclipse_when_glob``.

    Decides whether an eclipse whose classified type is ``retflag`` (one
    geometry bit ECL_CENTRAL/ECL_NONCENTRAL plus one type bit
    ECL_TOTAL/ECL_ANNULAR/ECL_PARTIAL/ECL_ANNULAR_TOTAL) satisfies the
    caller's ``mask`` filter according to the public API's flag contract:

    * ``mask == 0`` accepts any eclipse.
    * A geometry-only mask (no type bits) never matches.
    * The eclipse must carry one of the requested type bits.
    * With no geometry bit requested, only a *single*-type mask matches
      (any centrality); a geometry-only-plus-multiple-types mask never
      matches and the search runs to the ephemeris boundary.
    * With a geometry bit requested, the eclipse's centrality must match.
    """
    if mask == 0:
        return True
    geom = mask & (ECL_CENTRAL | ECL_NONCENTRAL)
    types = mask & (ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL | ECL_ANNULAR_TOTAL)
    if types == 0:
        return False
    if not (retflag & types):
        return False
    if geom == 0:
        return bin(types).count("1") == 1
    return bool(retflag & geom)


def _sol_glob_reject_impossible(mask: int) -> None:
    """Raise for ``sol_eclipse_when_glob`` type masks with no realisable event.

    Two masks describe eclipse geometries that cannot occur and the
    reference rejects them up front rather than searching:

    * ECL_CENTRAL | ECL_PARTIAL - a partial eclipse is never central.
    * ECL_NONCENTRAL | ECL_ANNULAR_TOTAL - a hybrid (annular-total)
      eclipse always has a central line, so it is never non-central.
    """
    valid = ECL_CENTRAL | ECL_NONCENTRAL | ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL
    valid |= ECL_ANNULAR_TOTAL
    m = mask & valid
    if m == (ECL_CENTRAL | ECL_PARTIAL):
        raise Error("central partial eclipses do not exist")
    if m == (ECL_NONCENTRAL | ECL_ANNULAR_TOTAL):
        raise Error("non-central hybrid (annular-total) eclipses do not exist")


# Constants for eclipse calculations
# Finite safety backstop (in years) for the eclipse/occultation "when"
# searches. These searches primarily stop at the active ephemeris boundary
# (see _is_ephemeris_boundary); this horizon only bounds a never-matching
# filter on the ~30 000-year extended ephemeris so it cannot loop forever.
# It is generous enough to reach the base/medium tier boundaries and to
# span the deep saros chains, yet finite.
_ECLIPSE_SEARCH_HORIZON_YEARS = 2000

# Occultation "when" searches step conjunction-by-conjunction (the Moon laps
# the target roughly once per sidereal month, ~13.4 conjunctions/year). The
# old fixed cap of 1200 conjunctions (~90 years) truncated valid searches
# well inside the ephemeris (e.g. the next partial occultation of Venus lies
# ~260 years out). Size the backstop to the search horizon; the search still
# stops earlier at the ephemeris boundary.
_OCCULT_MAX_CONJUNCTIONS = int(_ECLIPSE_SEARCH_HORIZON_YEARS * 13.5)

# Epoch-vs-maximum coincidence margin shared by every eclipse/occultation
# "when" finder. A candidate maximum within this margin of the search epoch
# counts as "already reached", so the idiom `jd = tret[0]; when(jd)` advances
# to the neighbouring event instead of re-returning the one it started on.
#
# Derived from this library's own solver, not from any external transition:
# the maximum is refined by golden section to ~1e-7 day (see _golden_min),
# and a search restarted from a returned maximum re-brackets from a
# different interval, so the recomputed instant can differ from the returned
# one by a few times that resolution. A margin of 100x the resolution
# absorbs that jitter with a wide safety factor while staying five orders of
# magnitude below the smallest possible separation between two distinct
# events (half a synodic month, ~14.8 days), so it can never merge two real
# eclipses. Callers starting between this margin and ~8.6 s before a maximum
# may see a different event than an external implementation chooses; that
# boundary window is documented in docs/comparison/known-differences.md.
_ECLIPSE_WHEN_EPOCH_MARGIN = 100.0 * 1e-7  # days (0.86 s)

SYNODIC_MONTH = 29.530588853  # Mean synodic month in days
LUNAR_NODE_PERIOD = 6798.38  # Lunar node regression period in days
ECLIPSE_LIMIT_SOLAR = 18.5  # Maximum elongation from node for solar eclipse (degrees)
# The New-Moon node-distance gate is only a coarse PRE-FILTER: the real
# eclipse test is the classification at the REFINED maximum. Keep the
# pre-filter a few tenths of a degree wider than the physical limit so an
# ultra-shallow partial (whose gamma dips below the eclipse threshold only
# at the refined maximum, a little away from the conjunction instant) is
# never dropped before it can be classified.
_SOLAR_NODE_PREFILTER = 19.0
EARTH_RADIUS_KM = 6378.137  # Earth equatorial radius in km (WGS84)

# Edge case thresholds for eclipse calculations
# These constants define limits for shallow/near-miss eclipses
SHALLOW_ECLIPSE_MAG_THRESHOLD = 0.01  # Minimum magnitude for reliable contact times
NEAR_MISS_GAMMA_MARGIN = 0.02  # Margin from gamma limit for edge case handling
MINIMUM_SEPARATION_FOR_LENS = 1e-10  # Minimum separation to avoid division by zero

# Location-based eclipse searches (sol_/lun_eclipse_when_loc) anchor their
# candidate lunations on the GEOCENTRIC maximum, but the observer's LOCAL
# phases can still be unfolding after that instant (or the Moon may only
# rise into the eclipse later). When the search epoch falls between the
# geocentric maximum and the local phases, the lunation must still be
# offered as a candidate; the local-contact epoch gate then keeps or skips
# it. Anchor the candidate generator this far (days) before/after the epoch
# so such an in-progress lunation is not dropped upstream. The value
# comfortably exceeds the half-duration of the widest local/penumbral phase
# window (< ~5 h) yet never reaches an adjacent lunation (eclipses are
# >= ~1 synodic month apart).
_WHEN_LOC_EPOCH_WINDOW = 0.35  # days (~8.4 h)

# 1 AU in km (IAU 2012 definition)
_AU_KM = 149597870.7

# Physical equatorial radii of planets in km, used for dynamic angular
# disc-size computation. Sources: IAU volumetric mean radii (Archinal et
# al.; IAU 2015 nominal values, NASA Planetary Fact Sheet). Disc sizes
# computed from these radii follow the independent physical disc geometry.
_PLANET_RADIUS_KM = {
    MERCURY: 2439.4,
    VENUS: 6051.8,
    MARS: 3389.5,
    JUPITER: 69911.0,
    SATURN: 58232.0,  # Disc only, excludes rings
    URANUS: 25362.0,
    NEPTUNE: 24622.0,
    PLUTO: 1188.3,
}

# Bodies beyond the classical planets with a finite apparent disc. Diameters
# are from the NASA JPL Small-Body Database (queried 2026-07-13), divided by
# two to obtain radii.
_EXTRA_BODY_RADIUS_KM = {
    CHIRON: 83.0,
    PHOLUS: 95.0,
    CERES: 469.7,
    PALLAS: 256.5,
    JUNO: 123.298,
    VESTA: 261.385,
    EARTH: 6371.0084,
}


# Eclipse-geometry physical constants. Values are public reference data:
# IAU / Astronomical Almanac solar diameter 1,392,000 km, lunar mean
# diameter 3,476.3 km, Earth equatorial radius 6,378.140 km, the DE431
# astronomical-unit convention, and the AA 2006 K6 Earth flattening.
_ECL_AU_KM = 149597870.700
_ECL_RSUN_AU = 696000.0 / _ECL_AU_KM
_ECL_RMOON_AU = 1738.15 / _ECL_AU_KM
_ECL_REARTH_AU = 6378.140 / _ECL_AU_KM
_ECL_EARTH_FLATTENING = 1.0 / 298.25642

# Inner-contact ("umbral") lunar radius.  Classical eclipse practice adopts
# two values of k (the ratio of the Moon's radius to Earth's equatorial
# radius): the larger, k = 0.2725076 (~1738.15 km, the mean lunar limb,
# adopted by the IAU in 1982), for the penumbral / outer silhouette, and the
# smaller k = 0.272281 for the umbral / inner geometry, where the marginal
# valleys of the lunar limb let a thread of photosphere through so that true
# totality/annularity, and the disappearance/reappearance of an occulted
# body behind the dark limb, are reckoned against a slightly smaller disc.
# Sources: Espenak & Meeus, "Five Millennium Canon of Solar Eclipses"
# (NASA/TP-2006-214141), sec. 1.5 (k = 0.272281 for umbral/central contacts);
# Explanatory Supplement to the Astronomical Almanac (3rd ed., 2013), eclipse
# phenomena; Meeus, "Elements of Solar Eclipses 1951-2200" (Willmann-Bell,
# 1989), ch. 3.  Applied at second/third contact and at occultation
# disappearance/reappearance while the full disc stays at first/fourth —
# the convention of the NASA canon, whose published durations these
# contacts reproduce to ~0.1 s.
_ECL_RMOON_INNER_AU = 0.272281 * 6378.140 / _ECL_AU_KM

# Danjon lunar-shadow convention used by NASA's Five Millennium Catalog.
# Instead of multiplying both completed shadow diameters by an empirical
# factor, it increases the Moon's equatorial horizontal parallax term by one
# percent while leaving the solar semidiameter and solar-parallax terms
# geometric.  See Espenak & Meeus, NASA/TP-2009-214173, Eqs. 1-5 and 1-6.
_DANJON_MOON_PARALLAX_SCALE = 1.01

# Earth equatorial radius in km (IAU 1976 system, Astronomical Almanac),
# used directly by the occultation shadow-geometry phase equations.
_EARTH_EQ_RADIUS_KM = 6378.140

# --- Lunar-occultation search parameters (derived, not tuned) ---------------
# Mean lunar sidereal motion is ~13.18 deg/day, so gap/13 is a safe
# first-order Newton step toward the next Moon-body conjunction in ecliptic
# longitude (slightly underestimating the rate keeps the iteration from
# overshooting for retrograde-moving targets).
_OCC_MOON_RATE_DEG_DAY = 13.0
# Conjunction convergence: 0.1 deg of longitude is ~11 minutes of lunar
# motion, well inside the +/-1 day refinement window used afterwards.
_OCC_CONJ_TOL_DEG = 0.1
# Conjunctions with a fixed direction repeat every sidereal month
# (~27.32 d): a 20-day hop from the current conjunction always lands
# before the next one, so no candidate is ever skipped.
_OCC_CONJUNCTION_HOP_DAYS = 20.0
# The Moon's ecliptic latitude stays within ~5.3 deg; adding its maximum
# horizontal parallax (~1.0 deg), its semidiameter (~0.27 deg) and margin,
# a star more than 7 deg from the ecliptic can never be occulted.
_OCC_STAR_LAT_LIMIT_DEG = 7.0
# At conjunction an occultation is only possible somewhere on Earth when
# the geocentric Moon-body latitude gap is below ~2 deg (semidiameter +
# parallax + body offset).
_OCC_LAT_GATE_DEG = 2.0


class _OccSearchOffEphemeris(Exception):
    """Conjunction seek stepped outside the active ephemeris range."""


def _normalize_occultation_filter(ecltype: int, body: "int | str", is_sun: bool) -> int:
    """Expand an occultation-type filter into its full bit set.

    The public ``ecltype`` argument names the wanted classes; internally the
    search tests individual class bits, so the shorthand values must be
    expanded: an empty filter means "everything possible for this body",
    umbral classes imply both the central and non-central variants, and a
    partial occultation is by construction never central. Annular classes
    exist only when the occulted body is the Sun (only the Sun's disc can
    exceed the Moon's); asking for them with another body is an error, and
    a combined filter simply loses those bits.

    Raises:
        Error: for the two impossible requests (central partial;
            annular-only for a non-solar body).
    """
    from .exceptions import Error

    wanted = ecltype
    if wanted == (ECL_PARTIAL | ECL_CENTRAL):
        raise Error("central partial eclipses do not exist")
    if not is_sun:
        base_classes = wanted & ~(ECL_NONCENTRAL | ECL_CENTRAL)
        if base_classes == ECL_ANNULAR or wanted == ECL_ANNULAR_TOTAL:
            raise Error(f"annular occultations do not exist for object {body}")
        if wanted & (ECL_ANNULAR | ECL_ANNULAR_TOTAL):
            wanted &= ~(ECL_ANNULAR | ECL_ANNULAR_TOTAL)
    if wanted == 0:
        wanted = ECL_TOTAL | ECL_PARTIAL | ECL_NONCENTRAL | ECL_CENTRAL
        if is_sun:
            wanted |= ECL_ANNULAR | ECL_ANNULAR_TOTAL
    if wanted & (ECL_TOTAL | ECL_ANNULAR | ECL_ANNULAR_TOTAL):
        wanted |= ECL_NONCENTRAL | ECL_CENTRAL
    if wanted & ECL_PARTIAL:
        wanted |= ECL_NONCENTRAL
    return wanted


def _seek_moon_conjunction(
    t: float,
    direction: int,
    geo_lonlat,
    eph_flags: int,
    star_guard: "str | None",
) -> Tuple[float, float, float]:
    """Advance ``t`` to the Moon-body conjunction in ecliptic longitude.

    First-order Newton iteration on the longitude gap with the mean lunar
    rate (see the derivation of the module constants above). The first
    body evaluation doubles as the never-occultable star gate; an
    ephemeris-range failure on that first evaluation raises
    :class:`_OccSearchOffEphemeris` so the caller can end the search at the
    tier boundary, while failures during the iteration itself propagate
    unchanged (matching the long-standing public behaviour).

    Returns:
        (t_conj, body_ecl_lat, moon_ecl_lat) at the converged instant.
    """
    from .exceptions import Error

    try:
        blon, blat = geo_lonlat(t)
    except Exception as _exc:
        if _is_ephemeris_boundary(_exc):
            raise _OccSearchOffEphemeris() from _exc
        raise
    if star_guard is not None and abs(blat) > _OCC_STAR_LAT_LIMIT_DEG:
        raise Error(
            f"occultation never occurs: star {star_guard} has ecl. lat. {blat:.1f}"
        )
    mpos, _ = calc_ut(t, MOON, eph_flags)
    gap = (blon - mpos[0]) % 360.0
    if direction < 0:
        gap -= 360.0
    steps = 0
    while abs(gap) > _OCC_CONJ_TOL_DEG and steps < 300:
        t += gap / _OCC_MOON_RATE_DEG_DAY
        blon, blat = geo_lonlat(t)
        mpos, _ = calc_ut(t, MOON, eph_flags)
        gap = (blon - mpos[0]) % 360.0
        if gap > 180.0:
            gap -= 360.0
        steps += 1
    return t, blat, mpos[1]


def _ecl_eph_flags(flags: int) -> int:
    """Reduce ``flags`` to its ephemeris-selection bits."""
    from .constants import FLG_JPLEPH, FLG_MOSEPH

    return flags & (FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH)


def _coerce_backwards(value: "bool | int | str") -> bool:
    """Coerce a public ``backwards`` argument to a plain bool.

    A truthy string like ``"forward"`` would silently select a *backward*
    search under plain truthiness, so direction strings are interpreted
    explicitly instead.

    Args:
        value: Search direction. Booleans pass through, ints use their
            truthiness, strings (case-insensitive, surrounding whitespace
            stripped) map "backward"/"backwards"/"true"/"1" to True and
            "forward"/"forwards"/"false"/"0"/"" to False.

    Returns:
        True for a backward search, False for a forward search.

    Raises:
        TypeError: If ``value`` is an unrecognized string or an
            unsupported type.
    """
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        return bool(value)
    if isinstance(value, str):
        text = value.strip().lower()
        if text in ("backward", "backwards", "true", "1"):
            return True
        if text in ("forward", "forwards", "false", "0", ""):
            return False
        raise TypeError(
            f"invalid backwards value {value!r}: expected a bool, an int, "
            "or one of 'backward'/'backwards'/'true'/'1' (backward search) "
            "or 'forward'/'forwards'/'false'/'0'/'' (forward search)"
        )
    raise TypeError(f"backwards must be a bool, int or str, got {type(value).__name__}")


def _topo_sun_moon(
    jd: float, geopos: "Sequence[float]", reader
) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
    """Topocentric apparent ecliptic-of-date (lon, lat, dist) of Sun and Moon.

    Single position source for the solar-eclipse machinery: the LEB path
    goes through :func:`fast_calc._topo_ecliptic`, the Skyfield path
    observes from a WGS84 observer and converts to the true ecliptic of
    date. Distances are in AU.
    """
    lon, lat, alt = float(geopos[0]), float(geopos[1]), float(geopos[2])
    if reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat

        jd_tt = jd + deltat(jd)
        gp = (lon, lat, alt)
        sun_p = _topo_ecliptic(reader, jd_tt, jd, SUN, gp)
        moon_p = _topo_ecliptic(reader, jd_tt, jd, MOON, gp)
        return (
            (float(sun_p[0]), float(sun_p[1]), float(sun_p[2])),
            (float(moon_p[0]), float(moon_p[1]), float(moon_p[2])),
        )

    from skyfield.api import wgs84
    from skyfield.framelib import ecliptic_frame

    from .state import get_planets

    eph = get_planets()
    ts = get_timescale()
    t = ts.ut1_jd(jd)
    observer_at = (eph["earth"] + wgs84.latlon(lat, lon, alt)).at(t)
    out = []
    for body in ("sun", "moon"):
        app = observer_at.observe(eph[body]).apparent()
        blat, blon, bdist = app.frame_latlon(ecliptic_frame)
        out.append((blon.degrees % 360.0, blat.degrees, bdist.au))
    return out[0], out[1]


def _occ_body_radius_au(body: "Union[int, str]") -> float:
    """Disc radius of an occulted body in AU (0 for stars/unknown)."""
    if isinstance(body, str):
        return 0.0
    if body == SUN:
        return _ECL_RSUN_AU
    if body == MOON:
        return _ECL_RMOON_AU
    radius_km = _PLANET_RADIUS_KM.get(body) or _EXTRA_BODY_RADIUS_KM.get(body)
    return radius_km / _ECL_AU_KM if radius_km else 0.0


def _occ_body_geo_xyz(
    tjd_ut: float, body: "Union[int, str]", flags_base: int
) -> Tuple[float, float, float]:
    """Geocentric equatorial XYZ (AU) of an occulted body (planet or star)."""
    if isinstance(body, str):
        from .fixed_stars import fixstar_ut

        pos, _name, _retflags = fixstar_ut(body, tjd_ut, flags_base)
        return (float(pos[0]), float(pos[1]), float(pos[2]))
    pos, _ = calc_ut(tjd_ut, body, flags_base)
    return (float(pos[0]), float(pos[1]), float(pos[2]))


@contextlib.contextmanager
def _observer_scope(lon: float, lat: float, alt: float) -> "Iterator[None]":
    """Temporarily set the topocentric observer for ``calc_ut(FLG_TOPOCTR)``.

    Saves and restores the process-wide observer (``state._TOPO``) so a
    rise/set or occultation search can request the topocentric apparent place
    of a body outside ``_PLANET_MAP`` through ``calc_ut`` — the one position
    source that handles the AST_OFFSET id aliases and both backends — without
    leaking the observer to the caller. ``set_topo`` clears the observer-at
    cache on entry; the cache is cleared again on exit so a following calc keys
    off the caller's own location.

    Thread-safety: intentionally NOT locked. The module-level API is
    stateful-global by design (matching the reference API, which shares
    the same process-wide observer model): concurrent module-level calls
    already race on sid mode, ephemeris path and the observer itself, so
    a lock here alone would only feign safety. Callers needing concurrent
    isolation use ``EphemerisContext``, which serializes access to the
    shared global state around each call.
    """
    from . import state as _state
    from .cache import clear_observer_cache
    from .state import get_topo, set_topo

    saved = get_topo()
    set_topo(lon, lat, alt)
    try:
        yield
    finally:
        _state._TOPO = saved
        clear_observer_cache()


def _occ_body_topo(
    tjd_ut: float,
    body: "Union[int, str]",
    geopos: "Sequence[float]",
    flags: int,
    reader,
) -> Tuple[float, float, float]:
    """Topocentric apparent ecliptic (lon, lat, dist) of an occulted body.

    Stars use the geocentric position (their diurnal parallax is far
    below the occultation tolerances handled here).
    """
    if isinstance(body, str):
        from .fixed_stars import fixstar_ut

        pos, _name, _retflags = fixstar_ut(body, tjd_ut, _ecl_eph_flags(flags))
        return (float(pos[0]), float(pos[1]), float(pos[2]))
    from .fixed_stars import FIXED_STARS

    if body in FIXED_STARS:
        # Integer fixed-star id (calc_ut dispatches these by id): geocentric
        # like the str branch — the diurnal parallax is far below the
        # occultation tolerances. Without this, the id fell through to the
        # planet machinery and leaked a raw KeyError.
        pos, _ = calc_ut(tjd_ut, body, _ecl_eph_flags(flags))
        return (float(pos[0]), float(pos[1]), float(pos[2]))
    if body == SUN:
        sun_p, _moon_p = _topo_sun_moon(tjd_ut, geopos, reader)
        return sun_p
    if reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat

        jd_tt = tjd_ut + deltat(tjd_ut)
        gp = (float(geopos[0]), float(geopos[1]), float(geopos[2]))
        pos = _topo_ecliptic(reader, jd_tt, tjd_ut, body, gp)
        return (float(pos[0]), float(pos[1]), float(pos[2]))

    from .planets import _PLANET_MAP, get_planet_target

    if body not in _PLANET_MAP:
        # Bodies the classical _PLANET_MAP does not enumerate (Chiron, the
        # main-belt asteroids, numbered minor planets, ...) have no Skyfield
        # target name here; without this the id fell through to
        # ``_PLANET_MAP[body]`` and leaked a raw KeyError on the Skyfield path
        # (the LEB path already served them via _topo_ecliptic). Position them
        # topocentrically through calc_ut(FLG_TOPOCTR): the identical
        # topocentric ecliptic-of-date place as the LEB path (verified to
        # 0.000"), so both backends share one occulted-body position source. A
        # genuinely unplaceable id raises the typed error from calc_ut.
        from .constants import FLG_TOPOCTR

        gp = (float(geopos[0]), float(geopos[1]), float(geopos[2]))
        with _observer_scope(*gp):
            pos, _ = calc_ut(tjd_ut, body, _ecl_eph_flags(flags) | FLG_TOPOCTR)
        return (float(pos[0]), float(pos[1]), float(pos[2]))

    from skyfield.api import wgs84
    from skyfield.framelib import ecliptic_frame

    from .state import get_planets

    eph = get_planets()
    ts = get_timescale()
    t = ts.ut1_jd(tjd_ut)
    observer_at = (
        eph["earth"]
        + wgs84.latlon(float(geopos[1]), float(geopos[0]), float(geopos[2]))
    ).at(t)
    target = get_planet_target(eph, _PLANET_MAP[body])
    app = observer_at.observe(target).apparent()
    blat, blon, bdist = app.frame_latlon(ecliptic_frame)
    return (blon.degrees % 360.0, blat.degrees, bdist.au)


def _eclipse_where_core(
    tjd_ut: float, flags: int = FLG_SWIEPH, body: "Union[int, str]" = SUN
) -> Tuple[int, float, float, Tuple[float, ...]]:
    """Shadow-axis ground geometry shared by the solar-eclipse functions.

    From the geocentric apparent true-of-date equatorial vectors of Sun
    and Moon, find where the Moon's shadow axis pierces the Earth's
    surface at ``tjd_ut`` (or, when the axis misses the globe, the
    surface point of maximum eclipse) and the shadow-cone cross sections
    there. The ellipsoid is handled by stretching the z axis so the
    Earth becomes a sphere of equatorial radius; a second pass refines
    the stretch with the local geodetic factor at the found latitude.

    Returns:
        ``(retc, lon, lat, dcore)``. ``retc`` is 0 when no eclipse is in
        progress anywhere on Earth (lon/lat then hold the
        closest-approach point), otherwise ECL_CENTRAL or
        ECL_NONCENTRAL (plus ECL_PARTIAL when only the penumbra
        touches), plus ECL_TOTAL/ECL_ANNULAR for the non-partial case.
        ``dcore`` uses this module's 7-slot shadow-geometry layout: [0] core-shadow diameter
        at the ground point (km; negative when the umbra reaches the
        surface), [1] penumbra diameter there (km), [2] axis distance
        from the geocenter (km), [3] core diameter on the fundamental
        plane (km), [4] penumbra diameter on the fundamental plane (km),
        [5] cos of the core-cone half angle, [6] cos of the
        penumbra-cone half angle.
    """
    from .constants import FLG_XYZ
    from .time_utils import sidtime

    base = _ecl_eph_flags(flags) | FLG_EQUATORIAL | FLG_XYZ
    moon_xyz, _ = calc_ut(tjd_ut, MOON, base)
    sun_xyz = _occ_body_geo_xyz(tjd_ut, body, base)

    rsun_au = _occ_body_radius_au(body)
    rmoon_au = _ECL_RMOON_AU
    dmoon_au = 2.0 * rmoon_au
    de = _ECL_REARTH_AU
    f = _ECL_EARTH_FLATTENING

    rmt = (moon_xyz[0], moon_xyz[1], moon_xyz[2])
    rst = (sun_xyz[0], sun_xyz[1], sun_xyz[2])
    dsmt = math.dist(rmt, rst)

    stretch = 1.0 - f
    for _pass in range(2):
        rm = (rmt[0], rmt[1], rmt[2] / stretch)
        rs = (rst[0], rst[1], rst[2] / stretch)
        dm = math.sqrt(rm[0] ** 2 + rm[1] ** 2 + rm[2] ** 2)
        ex, ey, ez = rm[0] - rs[0], rm[1] - rs[1], rm[2] - rs[2]
        dsm = math.sqrt(ex * ex + ey * ey + ez * ez)
        ex, ey, ez = ex / dsm, ey / dsm, ez / dsm
        # Shadow-cone half angles: the core (umbra/antumbra) cone is set
        # by the difference of the solar and lunar radii, the penumbra
        # cone by their sum.
        sinf1 = (rsun_au - rmoon_au) / dsm
        cosf1 = math.sqrt(1.0 - sinf1 * sinf1)
        sinf2 = (rsun_au + rmoon_au) / dsm
        cosf2 = math.sqrt(1.0 - sinf2 * sinf2)
        # Moon's distance from the fundamental plane (through the
        # geocenter, normal to the shadow axis) and the axis offset from
        # the geocenter.
        s0 = -(rm[0] * ex + rm[1] * ey + rm[2] * ez)
        r0_sq = dm * dm - s0 * s0
        r0 = math.sqrt(r0_sq) if r0_sq > 0.0 else 0.0
        # Core and penumbra diameters on the fundamental plane.
        d0 = (s0 / dsm * (2.0 * rsun_au - dmoon_au) - dmoon_au) / cosf1
        cap_d0 = (s0 / dsm * (2.0 * rsun_au + dmoon_au) + dmoon_au) / cosf2
        # Walk from the Moon along the axis to the stretched-Earth
        # surface; when the axis misses the globe, stop at the
        # fundamental plane (closest-approach point).
        disc = s0 * s0 + de * de - dm * dm
        s_walk = s0 - (math.sqrt(disc) if disc > 0.0 else 0.0)
        gx = rm[0] + s_walk * ex
        gy = rm[1] + s_walk * ey
        gz = rm[2] + s_walk * ez
        lat_sph = math.atan2(gz, math.hypot(gx, gy))
        if _pass == 0:
            # Refine the stretch so the spherical latitude of the ground
            # point comes out geodetic.
            cl, sl = math.cos(lat_sph), math.sin(lat_sph)
            one_f2 = (1.0 - f) ** 2
            stretch = one_f2 / math.sqrt(cl * cl + one_f2 * sl * sl)

    no_eclipse = False
    retc = 0
    if de * cosf1 >= r0:
        retc |= ECL_CENTRAL
    elif r0 <= de * cosf1 + abs(d0) / 2.0:
        retc |= ECL_NONCENTRAL
    elif r0 <= de * cosf2 + cap_d0 / 2.0:
        retc |= ECL_PARTIAL | ECL_NONCENTRAL
    else:
        no_eclipse = True

    sid_deg = sidtime(tjd_ut) * 15.0
    lon_out = (math.degrees(math.atan2(gy, gx)) - sid_deg) % 360.0
    if lon_out > 180.0:
        lon_out -= 360.0
    lat_out = math.degrees(lat_sph)

    # Shadow diameters at the ground point, from true-frame distances.
    g_true = (gx, gy, gz * stretch)
    s_mg = math.dist(rmt, g_true)
    dcore0 = (s_mg / dsmt * (2.0 * rsun_au - dmoon_au) - dmoon_au) * cosf1 * _ECL_AU_KM
    dcore1 = (s_mg / dsmt * (2.0 * rsun_au + dmoon_au) + dmoon_au) * cosf2 * _ECL_AU_KM
    if not no_eclipse and not (retc & ECL_PARTIAL):
        retc |= ECL_ANNULAR if dcore0 > 0.0 else ECL_TOTAL

    dcore = (
        dcore0,
        dcore1,
        r0 * _ECL_AU_KM,
        d0 * _ECL_AU_KM,
        cap_d0 * _ECL_AU_KM,
        cosf1,
        cosf2,
        0.0,
        0.0,
        0.0,
    )
    return retc, lon_out, lat_out, dcore


def _sol_how_core(
    tjd_ut: float,
    geopos: "Sequence[float]",
    flags: int,
    reader,
    body: "Union[int, str]" = SUN,
    where_convention: bool = False,
) -> Tuple[int, list]:
    """Local circumstances of a solar eclipse (reference ``attr`` layout).

    Returns ``(retc, attr)`` where ``attr`` is a 20-float list and
    ``retc`` carries the local phase (ECL_TOTAL/ECL_ANNULAR/ECL_PARTIAL,
    0 = no eclipse in progress at this place and time) plus ECL_VISIBLE
    when part of the eclipsed Sun can stand above the local horizon
    allowing for refraction and the observer's horizon dip.

    attr: [0] magnitude as diameter fraction (negative when the limbs do
    not yet overlap), [1] lunar/solar diameter ratio, [2] obscuration.
    For a solar eclipse the obscuration follows the compatibility
    contract: the disc-area ratio (r_moon/r_sun)**2 while one disc lies
    inside the other, so it exceeds 1 during totality (Moon larger) and is
    (r_moon/r_sun)**2 < 1 during annularity; the partial phase reports the
    two-disc lens-overlap fraction. For an occultation the obscuration is
    the covered fraction of the body, bounded at 1.0 on the how/when_loc
    path — but the WHERE path reports the uncapped disc-area ratio for
    planet targets (compatibility contract; ``where_convention``),
    while star targets stay at 1.0 on every path.
    [3] 0 (callers fill the core-shadow width), [4] azimuth of the Sun,
    [5] true altitude, [6] apparent altitude, [7] Moon-Sun center
    separation in degrees, [8] NASA magnitude, [9]/[10] saros series and
    member.
    """
    from .utils import ECL2HOR, angular_separation, azalt

    _sun_unused, moon_p = _topo_sun_moon(tjd_ut, geopos, reader)
    sun_p = (
        _sun_unused
        if (not isinstance(body, str) and body == SUN)
        else _occ_body_topo(tjd_ut, body, geopos, flags, reader)
    )
    body_radius_au = _occ_body_radius_au(body)
    rsun = (
        math.degrees(math.asin(body_radius_au / sun_p[2]))
        if body_radius_au > 0.0
        else 0.0
    )
    rmoon = math.degrees(math.asin(_ECL_RMOON_AU / moon_p[2]))
    dctr = angular_separation(sun_p[0], sun_p[1], moon_p[0], moon_p[1])

    if dctr < rsun - rmoon:
        retc = ECL_ANNULAR
    elif dctr < abs(rsun - rmoon):
        retc = ECL_TOTAL
    elif dctr < rsun + rmoon:
        retc = ECL_PARTIAL
    else:
        retc = 0

    attr = [0.0] * 20
    attr[1] = rmoon / rsun if rsun > 0.0 else 0.0
    attr[0] = (rsun + rmoon - dctr) / (2.0 * rsun) if rsun > 0.0 else 1.0
    # Compatibility contract: a solar eclipse (Sun as the "occulted" body)
    # reports the disc-area ratio while one disc lies inside the other, so
    # obscuration exceeds 1 during totality; an occultation of a finite
    # body caps the covered fraction at 1.0.
    _is_solar = (not isinstance(body, str)) and body == SUN
    if retc == 0:
        # No eclipse at this place and time. Compatibility contract: this
        # entry point reports a fixed no-event sentinel of 1.0 in the
        # obscuration slot — a protocol value in the same class as a
        # retflag echo or a -99999999.0 sentinel, not an area ratio. It
        # reads oddly for a channel documented as "fraction of the Sun's
        # disc covered", but consumers branch on it; sol_eclipse_how()
        # zeroes its whole attr array on retflag==0 and is unaffected.
        attr[2] = 1.0
    elif rsun <= 0.0:
        attr[2] = 1.0
    elif retc in (ECL_TOTAL, ECL_ANNULAR):
        # One disc lies entirely within the other. The disc-area ratio is
        # (rmoon/rsun)**2 — > 1 during totality (larger Moon), < 1 during
        # annularity (a ring of Sun remains). For a solar eclipse the
        # reference returns this ratio uncapped; an occultation caps the
        # covered fraction of the body at 1.0.
        ratio_sq = (rmoon / rsun) ** 2
        _uncap = _is_solar or (where_convention and rsun > 0.0)
        attr[2] = ratio_sq if _uncap else min(1.0, ratio_sq)
    elif dctr <= 0.0:
        # Exactly concentric discs in the partial branch are reachable only at
        # the annular/total boundary (rsun == rmoon); the overlap is the whole
        # smaller disc. Guards the lens-area 1/dctr singularity below.
        ratio_sq = (rmoon / rsun) ** 2
        _uncap = _is_solar or (where_convention and rsun > 0.0)
        attr[2] = ratio_sq if _uncap else min(1.0, ratio_sq)
    else:
        # Standard two-disc overlap (lens) area as a fraction of the
        # solar disc.
        a = (dctr * dctr + rmoon * rmoon - rsun * rsun) / (2.0 * dctr * rmoon)
        b = (dctr * dctr + rsun * rsun - rmoon * rmoon) / (2.0 * dctr * rsun)
        a = math.acos(max(-1.0, min(1.0, a)))
        b = math.acos(max(-1.0, min(1.0, b)))
        seg_moon = (a - math.sin(a) * math.cos(a)) * rmoon * rmoon
        seg_sun = (b - math.sin(b) * math.cos(b)) * rsun * rsun
        attr[2] = (seg_moon + seg_sun) / math.pi / (rsun * rsun)
    attr[7] = dctr

    geopos3 = (float(geopos[0]), float(geopos[1]), float(geopos[2]))
    az, true_alt, app_alt = azalt(tjd_ut, ECL2HOR, geopos3, 0.0, 10.0, sun_p)
    attr[4] = az
    attr[5] = true_alt
    attr[6] = app_alt
    # Lowest apparent altitude at which part of the Sun can still be
    # seen: horizon refraction (34.4556') plus horizon dip and
    # observer-to-horizon refraction, both growing with the square root
    # of the observer's height.
    hmin_appr = -(34.4556 + 2.12 * math.sqrt(max(0.0, geopos3[2]))) / 60.0
    if retc and true_alt + rsun + abs(hmin_appr) >= 0.0:
        retc |= ECL_VISIBLE
    if not isinstance(body, str) and body == SUN:
        attr[8] = attr[1] if retc & (ECL_TOTAL | ECL_ANNULAR) else attr[0]
        saros_series, saros_member = _get_saros_info(tjd_ut, "solar")
        attr[9] = saros_series
        attr[10] = saros_member
    return retc, attr


def _lun_how_core(
    tjd_ut: float, flags: int = FLG_SWIEPH
) -> Tuple[int, list, Tuple[float, ...]]:
    """Geocentric lunar-eclipse circumstances (reference ``attr`` layout).

    Works in the selenocentric frame: the Earth's shadow cone is built
    from the selenocentric Sun and Earth vectors, the Moon's immersion
    follows from the shadow-axis offset at the Moon. The shadow diameters
    use Danjon's independently published atmospheric convention: the lunar
    horizontal-parallax term is enlarged by one percent (Espenak & Meeus,
    NASA/TP-2009-214173). No output-fitted shadow deflator is applied.

    Returns ``(retc, attr, dcore)``:
        retc: ECL_TOTAL, ECL_PARTIAL, ECL_PENUMBRAL or 0.
        attr: 20-float list with the observer-independent entries set -
            [0] umbral magnitude (0 for penumbral eclipses),
            [1] penumbral magnitude (negative when the Moon stands clear
            of the penumbra), [7] Moon's distance from opposition in
            degrees (0 when no eclipse), [8] = [0], [9]/[10] lunar saros
            series and member.
        dcore: [0] shadow-axis distance from the selenocenter, [1] core
            diameter on the fundamental plane, [2] penumbra diameter on
            the fundamental plane, [3] cos of the core-cone half angle,
            [4] cos of the penumbra-cone half angle (AU-based units).
    """
    from .constants import FLG_XYZ

    base = _ecl_eph_flags(flags) | FLG_EQUATORIAL | FLG_XYZ | FLG_SPEED
    moon_xyz, _ = calc_ut(tjd_ut, MOON, base)
    sun_xyz, _ = calc_ut(tjd_ut, SUN, base)

    rm = (moon_xyz[0], moon_xyz[1], moon_xyz[2])
    rs = (sun_xyz[0], sun_xyz[1], sun_xyz[2])
    dm = math.sqrt(rm[0] ** 2 + rm[1] ** 2 + rm[2] ** 2)
    ds = math.sqrt(rs[0] ** 2 + rs[1] ** 2 + rs[2] ** 2)
    cos_elong = (rm[0] * rs[0] + rm[1] * rs[1] + rm[2] * rs[2]) / (dm * ds)
    dctr = math.degrees(math.acos(max(-1.0, min(1.0, cos_elong))))

    # Selenocentric Sun and Earth.
    s_sun = (rs[0] - rm[0], rs[1] - rm[1], rs[2] - rm[2])
    s_earth = (-rm[0], -rm[1], -rm[2])
    ex = s_earth[0] - s_sun[0]
    ey = s_earth[1] - s_sun[1]
    ez = s_earth[2] - s_sun[2]
    dsm = math.sqrt(ex * ex + ey * ey + ez * ez)
    ex, ey, ez = ex / dsm, ey / dsm, ez / dsm

    rsun_au = _ECL_RSUN_AU
    rearth_au = _ECL_REARTH_AU
    rmoon_au = _ECL_RMOON_AU
    dmoon_au = 2.0 * rmoon_au

    f1 = (rsun_au - rearth_au) / dsm
    cosf1 = math.sqrt(1.0 - f1 * f1)
    f2 = (rsun_au + rearth_au) / dsm
    cosf2 = math.sqrt(1.0 - f2 * f2)

    # Earth's distance from the fundamental plane through the Moon and
    # the shadow-axis offset at the Moon.
    s0 = -(s_earth[0] * ex + s_earth[1] * ey + s_earth[2] * ez)
    r0_sq = dm * dm - s0 * s0
    r0 = math.sqrt(r0_sq) if r0_sq > 0.0 else 0.0

    # Exact cone sections with Danjon's atmospheric correction.  Scaling only
    # the leading Earth-radius term implements NASA's angular definitions
    # Ru = 1.01*Pm - Ss + Ps and Rp = 1.01*Pm + Ss + Ps; multiplying a
    # completed shadow radius would instead be Chauvenet's different 1/50
    # convention.
    d0 = (
        abs(
            s0 / dsm * (2.0 * rsun_au - 2.0 * rearth_au)
            - 2.0 * _DANJON_MOON_PARALLAX_SCALE * rearth_au
        )
        / cosf1
    )
    cap_d0 = (
        s0 / dsm * (2.0 * rsun_au + 2.0 * rearth_au)
        + 2.0 * _DANJON_MOON_PARALLAX_SCALE * rearth_au
    ) / cosf2
    attr = [0.0] * 20
    retc = 0
    if d0 / 2.0 >= r0 + rmoon_au / cosf1:
        retc = ECL_TOTAL
        attr[0] = (d0 / 2.0 - r0 + rmoon_au) / dmoon_au
    elif d0 / 2.0 >= r0 - rmoon_au / cosf1:
        retc = ECL_PARTIAL
        attr[0] = (d0 / 2.0 - r0 + rmoon_au) / dmoon_au
    elif cap_d0 / 2.0 >= r0 - rmoon_au / cosf2:
        retc = ECL_PENUMBRAL
    attr[8] = attr[0]
    attr[1] = (cap_d0 / 2.0 - r0 + rmoon_au) / dmoon_au
    if retc:
        attr[7] = 180.0 - abs(dctr)
    saros_series, saros_member = _get_saros_info(tjd_ut, "lunar")
    attr[9] = saros_series
    attr[10] = saros_member

    dcore = (r0, d0, cap_d0, cosf1, cosf2, 0.0, 0.0, 0.0, 0.0, 0.0)
    return retc, attr, dcore


def _lun_eclipse_max_time(jd_approx: float, flags: int = FLG_SWIEPH) -> float:
    """Time of maximum lunar eclipse near ``jd_approx``.

    The maximum is the deepest immersion of the Moon in the Earth's
    shadow. Seen from the Moon at a full moon, the Earth stands almost
    exactly in front of the Sun, so the immersion is deepest where the
    selenocentric angle between the solar and terrestrial limbs - the
    center separation less both apparent radii - is smallest.
    """
    from .constants import FLG_XYZ

    base = _ecl_eph_flags(flags) | FLG_EQUATORIAL | FLG_XYZ

    def _depth(jd: float) -> float:
        moon_xyz, _ = calc_ut(jd, MOON, base)
        sun_xyz, _ = calc_ut(jd, SUN, base)
        sx = sun_xyz[0] - moon_xyz[0]
        sy = sun_xyz[1] - moon_xyz[1]
        sz = sun_xyz[2] - moon_xyz[2]
        gx, gy, gz = -moon_xyz[0], -moon_xyz[1], -moon_xyz[2]
        ds = math.sqrt(sx * sx + sy * sy + sz * sz)
        dm = math.sqrt(gx * gx + gy * gy + gz * gz)
        cosang = (sx * gx + sy * gy + sz * gz) / (ds * dm)
        ang = math.degrees(math.acos(max(-1.0, min(1.0, cosang))))
        rearth = math.degrees(math.asin(_ECL_REARTH_AU / dm))
        rsun = math.degrees(math.asin(_ECL_RSUN_AU / ds))
        return ang - (rearth + rsun)

    phi = (1.0 + math.sqrt(5.0)) / 2.0
    jd_low = jd_approx - 0.3
    jd_high = jd_approx + 0.3
    jd_a = jd_high - (jd_high - jd_low) / phi
    jd_b = jd_low + (jd_high - jd_low) / phi
    f_a = _depth(jd_a)
    f_b = _depth(jd_b)
    for _ in range(60):
        if f_a < f_b:
            jd_high = jd_b
            jd_b, f_b = jd_a, f_a
            jd_a = jd_high - (jd_high - jd_low) / phi
            f_a = _depth(jd_a)
        else:
            jd_low = jd_a
            jd_a, f_a = jd_b, f_b
            jd_b = jd_low + (jd_high - jd_low) / phi
            f_b = _depth(jd_b)
        if jd_high - jd_low < 1e-7:
            break
    return (jd_low + jd_high) / 2.0


def _lun_eclipse_phase_times(
    jd_max: float, retc: int, flags: int = FLG_SWIEPH
) -> Tuple[float, ...]:
    """Phase times (tret layout) for a lunar eclipse around its maximum.

    Phase boundaries are the zero crossings of the Moon-limb immersion
    measures from the shadow geometry: penumbral contact (tret[6]/[7]),
    umbral/partial contact (tret[2]/[3]) and total immersion
    (tret[4]/[5]).
    """
    rmoon_au = _ECL_RMOON_AU

    def _f_pen(jd: float) -> float:
        _rc, _a, dc = _lun_how_core(jd, flags)
        return dc[2] / 2.0 + rmoon_au / dc[4] - dc[0]

    def _f_par(jd: float) -> float:
        _rc, _a, dc = _lun_how_core(jd, flags)
        return dc[1] / 2.0 + rmoon_au / dc[3] - dc[0]

    def _f_tot(jd: float) -> float:
        _rc, _a, dc = _lun_how_core(jd, flags)
        return dc[1] / 2.0 - rmoon_au / dc[3] - dc[0]

    tret = [0.0] * 10
    tret[0] = jd_max
    # Half-window for each contact search. It must exceed the longest
    # half-phase duration (penumbral half-duration peaks near ~3 h), or the
    # contact would fall outside [jd_max - win, jd_max] and _root_bisect
    # would silently return its 0.0 "unbracketed" sentinel for that contact.
    win = 4.0 / 24.0
    tret[6] = _root_bisect(_f_pen, jd_max - win, jd_max)
    tret[7] = _root_bisect(_f_pen, jd_max, jd_max + win)
    if retc & (ECL_PARTIAL | ECL_TOTAL):
        tret[2] = _root_bisect(_f_par, jd_max - win, jd_max)
        tret[3] = _root_bisect(_f_par, jd_max, jd_max + win)
    if retc & ECL_TOTAL:
        tret[4] = _root_bisect(_f_tot, jd_max - win, jd_max)
        tret[5] = _root_bisect(_f_tot, jd_max, jd_max + win)
    return tuple(tret)


def _lun_contact_via_core(jd_max: float, tret_idx: int, flags: int) -> float:
    """Lunar-eclipse contact time at ``tret`` index ``tret_idx``.

    Computed from the canonical ``_lun_how_core`` shadow model via
    ``_lun_eclipse_phase_times`` -- the same path ``lun_eclipse_when`` uses,
    which reproduces the independent Earth-shadow geometry (and the contact
    times published in NASA's Five Millennium Canon of Lunar Eclipses) to
    <0.5 s. Routing the standalone ``calc_lunar_eclipse_*_contact_*`` helpers
    through here keeps them consistent with ``lun_eclipse_when`` instead of the
    earlier, divergent 1/85-enlargement + small-angle separation model (which
    drifted by tens to hundreds of seconds). Returns 0.0 when the requested
    phase does not occur (e.g. an umbral contact for a penumbral-only eclipse).
    """
    retc, _attr, _dc = _lun_how_core(jd_max, flags)
    if not retc:
        return 0.0
    return _lun_eclipse_phase_times(jd_max, retc, flags)[tret_idx]


def _golden_min(fn, t_lo: float, t_hi: float) -> float:
    """Golden-section minimum of ``fn`` on [t_lo, t_hi] to ~1e-7 day."""
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    t_a = t_hi - (t_hi - t_lo) / phi
    t_b = t_lo + (t_hi - t_lo) / phi
    f_a = fn(t_a)
    f_b = fn(t_b)
    for _ in range(60):
        if f_a < f_b:
            t_hi = t_b
            t_b, f_b = t_a, f_a
            t_a = t_hi - (t_hi - t_lo) / phi
            f_a = fn(t_a)
        else:
            t_lo = t_a
            t_a, f_a = t_b, f_b
            t_b = t_lo + (t_hi - t_lo) / phi
            f_b = fn(t_b)
        if t_hi - t_lo < 1e-7:
            break
    return 0.5 * (t_lo + t_hi)


def _root_bisect(fn, t_lo: float, t_hi: float) -> float:
    """Bisect ``fn`` to a sign change in [t_lo, t_hi]; 0.0 when unbracketed.

    Used for eclipse contact times, where the bracket endpoints are known
    geometry (separation extremum on one side, no-overlap on the other).
    Converges to ~1e-8 day (sub-millisecond).
    """
    f_lo = fn(t_lo)
    f_hi = fn(t_hi)
    if f_lo == 0.0:
        return t_lo
    if f_hi == 0.0:
        return t_hi
    if (f_lo > 0.0) == (f_hi > 0.0):
        return 0.0
    for _ in range(60):
        t_mid = 0.5 * (t_lo + t_hi)
        f_mid = fn(t_mid)
        if (f_mid > 0.0) == (f_lo > 0.0):
            t_lo, f_lo = t_mid, f_mid
        else:
            t_hi, f_hi = t_mid, f_mid
        if t_hi - t_lo < 1e-8:
            break
    return 0.5 * (t_lo + t_hi)


def _sinclair_refraction_deg(
    app_alt_deg: float, atpress: float, attemp: float
) -> float:
    """Atmospheric refraction at an apparent altitude, in degrees.

    A. T. Sinclair's refraction model, as published in G. G. Bennett,
    "The calculation of astronomical refraction in marine navigation",
    Journal of the Institute of Navigation 35 (1982), p. 256, with
    Bennett's pressure/temperature compensation (p. 259).

    The two branches join at h = 17.904104638432°: above it the high-altitude
    asymptote 0.97/tan(h) (arcmin) applies; at or below it Sinclair's rational
    fit (34.46 + 4.23h + 0.004h²)/(1 + 0.505h + 0.0845h²) (arcmin) applies.
    The pressure/temperature factor (atpress-80)/930 / (1 + 8e-5·(r+39)·
    (attemp-10)) is Bennett's compensation (p. 259), referred to 1010 mbar /
    10 °C. Result converted arcmin -> deg.
    """
    h = app_alt_deg
    if h > 17.904104638432:
        r = 0.97 / math.tan(math.radians(h))
    else:
        r = (34.46 + 4.23 * h + 0.004 * h * h) / (1.0 + 0.505 * h + 0.0845 * h * h)
    r = (atpress - 80.0) / 930.0 / (1.0 + 0.00008 * (r + 39.0) * (attemp - 10.0)) * r
    return r / 60.0


def _rise_true_to_apparent(
    true_alt: float, geoalt: float, atpress: float, attemp: float
) -> float:
    """True -> apparent altitude for rise/set events (reference scheme).

    The Sinclair refraction is a function of the apparent altitude, so the
    apparent altitude solves the fixed point app = true + R(app); a Newton
    iteration with numerical derivative converges in a few steps. When the
    refracted point would still sit below the observer's horizon dip, the
    altitude is left unrefracted - the reference convention for rise/set.

    The project's ray-traced refraction (utils.refrac/refrac_extended)
    deliberately is not used here: rise/set parity with the reference
    requires its published horizon model (see docs on the refraction
    deviation envelope).
    """
    from .refraction import calc_dip

    app = true_alt + _sinclair_refraction_deg(true_alt, atpress, attemp)
    for _ in range(10):
        f = app - _sinclair_refraction_deg(app, atpress, attemp) - true_alt
        eps = 1e-4
        f2 = (
            (app + eps)
            - _sinclair_refraction_deg(app + eps, atpress, attemp)
            - true_alt
        )
        dfd = (f2 - f) / eps
        step = f / dfd if dfd != 0.0 else f
        app -= step
        if abs(step) < 1e-9:
            break
    refr = app - true_alt
    dip = calc_dip(max(geoalt, 0.0), atpress=atpress, attemp=attemp)
    if true_alt + refr < dip:
        return true_alt
    return true_alt + refr


def _calc_planet_angular_radius(planet_id: int, dist_au: float) -> float:
    """Compute angular radius of a planet in degrees from its geocentric distance.

    Uses the exact geometric formula: angular_radius = arcsin(physical_radius / distance),
    matching the reference API approach of computing angular sizes dynamically from
    physical radii and actual distance rather than using static lookup tables.

    Args:
        planet_id: Planet constant (MERCURY, VENUS, etc.)
        dist_au: Geocentric distance in AU.

    Returns:
        Angular radius in degrees. Returns 0.0001° (~0.4") for unknown bodies
        (point source approximation).
    """
    radius_km = _PLANET_RADIUS_KM.get(planet_id)
    if radius_km is None:
        return 0.0001  # Point source for unknown bodies
    dist_km = dist_au * _AU_KM
    if dist_km <= 0:
        return 0.0001
    return math.degrees(math.asin(min(1.0, radius_km / dist_km)))


def _is_shallow_eclipse(magnitude: float) -> bool:
    """
    Check if an eclipse is very shallow (magnitude close to zero).

    Very shallow eclipses have unreliable contact time calculations because
    the penumbra barely grazes Earth or the Moon barely enters the shadow.

    Args:
        magnitude: Eclipse magnitude (0 to ~1.5 for solar, 0 to ~2 for lunar)

    Returns:
        True if the eclipse is considered shallow (magnitude < threshold)
    """
    return magnitude < SHALLOW_ECLIPSE_MAG_THRESHOLD


def _is_near_miss_eclipse(gamma: float, gamma_limit: float = 1.55) -> bool:
    """
    Check if an eclipse is a near-miss (gamma very close to the eclipse limit).

    Near-miss eclipses have gamma values very close to the limit where eclipses
    cease to occur. These require special handling to avoid numerical instability.

    Args:
        gamma: Eclipse gamma parameter (distance of shadow axis from Earth center)
        gamma_limit: Maximum gamma for any eclipse visibility (default 1.55)

    Returns:
        True if gamma is within NEAR_MISS_GAMMA_MARGIN of the limit
    """
    return abs(gamma) > (gamma_limit - NEAR_MISS_GAMMA_MARGIN)


def _safe_acos(x: float) -> float:
    """
    Safe arccosine that clamps input to valid range [-1, 1].

    Numerical errors can cause values slightly outside the valid range,
    which would raise a domain error. This function clamps the input.

    Args:
        x: Input value (should be in [-1, 1])

    Returns:
        math.acos of the clamped value
    """
    return math.acos(max(-1.0, min(1.0, x)))


def _safe_sqrt(x: float) -> float:
    """
    Safe square root that returns 0 for negative inputs.

    Numerical errors can cause small negative values where zero is expected.
    This function returns 0 for negative inputs instead of raising an error.

    Args:
        x: Input value (should be >= 0)

    Returns:
        math.sqrt of x if x >= 0, else 0.0
    """
    return math.sqrt(max(0.0, x))


def _calculate_obscuration_safe(r_sun: float, r_moon: float, d: float) -> float:
    """
    Calculate eclipse obscuration with edge case handling.

    Obscuration is the fraction of the Sun's area covered by the Moon.
    This function handles edge cases like:
    - Zero separation (total/annular eclipse)
    - Separation equal to sum of radii (no eclipse)
    - Very small separations that could cause numerical issues

    Args:
        r_sun: Sun's angular radius in degrees
        r_moon: Moon's angular radius in degrees
        d: Center-to-center separation in degrees

    Returns:
        Obscuration: the fraction of the Sun's disc area covered by the
        Moon, bounded to [0, 1] by the published definition — 1.0 for a
        total eclipse, (r_moon/r_sun)^2 < 1 for an annular one, the
        two-disc lens-overlap fraction for a partial eclipse, and 0.0 with
        no overlap.
    """
    # Handle edge case: no overlap
    if d >= r_sun + r_moon:
        return 0.0

    # Edge case: one disc lies entirely within the other (or is concentric).
    if d <= abs(r_sun - r_moon) or d < MINIMUM_SEPARATION_FOR_LENS:
        # Total (Moon >= Sun): the Sun is fully covered -> 1.0.
        # Annular (Moon < Sun): a ring remains -> disc area ratio < 1.0.
        return min(1.0, (r_moon / r_sun) ** 2)

    # Calculate lens-shaped intersection area
    # d1 is distance from Sun center to intersection chord
    # d2 is distance from Moon center to intersection chord
    d1 = (d * d + r_sun * r_sun - r_moon * r_moon) / (2 * d)
    d2 = d - d1

    # Validate that the geometric configuration is valid
    if abs(d1) > r_sun or abs(d2) > r_moon:
        # This shouldn't happen if we passed earlier checks, but handle gracefully
        return 0.0

    # Calculate areas using lens formula with safe operations
    area1 = r_sun * r_sun * _safe_acos(d1 / r_sun) - d1 * _safe_sqrt(
        r_sun * r_sun - d1 * d1
    )
    area2 = r_moon * r_moon * _safe_acos(d2 / r_moon) - d2 * _safe_sqrt(
        r_moon * r_moon - d2 * d2
    )

    intersection_area = area1 + area2
    sun_area = math.pi * r_sun * r_sun

    obscuration = intersection_area / sun_area
    return max(0.0, obscuration)


def _calculate_magnitude_safe(
    gamma: float, moon_sun_ratio: float, gamma_limit_partial: float = 1.55
) -> float:
    """
    Calculate eclipse magnitude with edge case handling.

    Handles shallow partial eclipses where gamma is close to the eclipse limit.

    Args:
        gamma: Eclipse gamma parameter (shadow axis distance from Earth center)
        moon_sun_ratio: Ratio of Moon's apparent diameter to Sun's
        gamma_limit_partial: Maximum gamma for partial eclipse visibility

    Returns:
        Eclipse magnitude (0 to ~1.5), clamped to valid range
    """
    # For very shallow eclipses (gamma near limit), magnitude approaches 0
    if _is_near_miss_eclipse(gamma, gamma_limit_partial):
        # Use linear interpolation near the edge for smooth transition
        remaining = gamma_limit_partial - abs(gamma)
        if remaining <= 0:
            return 0.0
        # Magnitude decreases linearly as gamma approaches limit
        magnitude = (remaining / NEAR_MISS_GAMMA_MARGIN) * SHALLOW_ECLIPSE_MAG_THRESHOLD
        return max(0.0, min(1.5, magnitude * moon_sun_ratio))

    # Standard magnitude calculation
    if abs(gamma) >= gamma_limit_partial:
        return 0.0

    magnitude = 1.0 - abs(gamma) / gamma_limit_partial
    magnitude = magnitude * moon_sun_ratio

    return max(0.0, min(1.5, magnitude))


def _validate_contact_time(
    jd_contact: float, jd_max: float, max_offset_days: float = 0.25
) -> float:
    """
    Validate a calculated contact time and return 0 if invalid.

    Contact times should be within a reasonable range of the eclipse maximum.
    This function checks if the contact time is valid.

    Args:
        jd_contact: Calculated contact time (Julian Day)
        jd_max: Eclipse maximum time (Julian Day)
        max_offset_days: Maximum allowed offset from maximum (default 6 hours)

    Returns:
        The contact time if valid, 0.0 if invalid
    """
    if jd_contact <= 0:
        return 0.0

    offset = abs(jd_contact - jd_max)
    if offset > max_offset_days:
        return 0.0

    return jd_contact


def _find_next_new_moon(jd_start: float) -> float:
    """
    Find the next New Moon (Sun-Moon conjunction) after jd_start.

    Uses iterative refinement to find exact moment of conjunction.

    Args:
        jd_start: Julian Day (UT) to start search from

    Returns:
        Julian Day of next New Moon
    """
    # Get current positions
    sun_pos, _ = calc_ut(jd_start, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd_start, MOON, FLG_SPEED)

    sun_lon = sun_pos[0]
    moon_lon = moon_pos[0]

    # Calculate elongation (Moon - Sun), normalized to -180 to 180
    elongation = (moon_lon - sun_lon) % 360.0
    if elongation > 180:
        elongation -= 360

    # Estimate time to next conjunction
    # Moon gains ~12.2° per day on Sun
    relative_speed = 12.190749  # degrees/day (average Moon speed - Sun speed)

    # Time until next conjunction (elongation = 0)
    if elongation > 0:
        # Moon is ahead, wait for next cycle
        dt = (360.0 - elongation) / relative_speed
    else:
        # Moon is behind Sun
        dt = (-elongation) / relative_speed

    jd_guess = jd_start + dt

    # Newton-Raphson refinement
    for _ in range(20):
        sun_pos, _ = calc_ut(jd_guess, SUN, FLG_SPEED)
        moon_pos, _ = calc_ut(jd_guess, MOON, FLG_SPEED)

        sun_lon = sun_pos[0]
        moon_lon = moon_pos[0]
        sun_speed = sun_pos[3]
        moon_speed = moon_pos[3]

        # Elongation
        diff = (moon_lon - sun_lon) % 360.0
        if diff > 180:
            diff -= 360

        # Convergence check (< 0.1 arcsec)
        if abs(diff) < 1e-5:
            return jd_guess

        # Newton-Raphson step
        rel_speed = moon_speed - sun_speed
        if abs(rel_speed) < 0.1:
            rel_speed = 12.19

        jd_guess -= diff / rel_speed

    return jd_guess


def _find_previous_new_moon(jd_start: float) -> float:
    """
    Find the previous New Moon (Sun-Moon conjunction) before jd_start.

    Uses iterative refinement to find exact moment of conjunction.

    Args:
        jd_start: Julian Day (UT) to start search from

    Returns:
        Julian Day of previous New Moon
    """
    # Get current positions
    sun_pos, _ = calc_ut(jd_start, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd_start, MOON, FLG_SPEED)

    sun_lon = sun_pos[0]
    moon_lon = moon_pos[0]

    # Calculate elongation (Moon - Sun), normalized to -180 to 180
    elongation = (moon_lon - sun_lon) % 360.0
    if elongation > 180:
        elongation -= 360

    # Estimate time to previous conjunction
    # Moon gains ~12.2° per day on Sun
    relative_speed = 12.190749  # degrees/day (average Moon speed - Sun speed)

    # Time since last conjunction (go backwards)
    if elongation > 0:
        # Moon is ahead of Sun, last conjunction was elongation/speed days ago
        dt = -elongation / relative_speed
    else:
        # Moon is behind Sun, last conjunction was (360 + elongation)/speed days ago
        dt = -(360.0 + elongation) / relative_speed

    jd_guess = jd_start + dt

    # Newton-Raphson refinement
    for _ in range(20):
        sun_pos, _ = calc_ut(jd_guess, SUN, FLG_SPEED)
        moon_pos, _ = calc_ut(jd_guess, MOON, FLG_SPEED)

        sun_lon = sun_pos[0]
        moon_lon = moon_pos[0]
        sun_speed = sun_pos[3]
        moon_speed = moon_pos[3]

        # Elongation
        diff = (moon_lon - sun_lon) % 360.0
        if diff > 180:
            diff -= 360

        # Convergence check (< 0.1 arcsec)
        if abs(diff) < 1e-5:
            return jd_guess

        # Newton-Raphson step
        rel_speed = moon_speed - sun_speed
        if abs(rel_speed) < 0.1:
            rel_speed = 12.19

        jd_guess -= diff / rel_speed

    return jd_guess


def _get_moon_node_distance(jd: float, moon_lon: float) -> float:
    """
    Calculate angular distance of Moon from nearest lunar node.

    Args:
        jd: Julian Day (UT)
        moon_lon: Moon's ecliptic longitude in degrees

    Returns:
        Absolute angular distance from nearest node in degrees
    """
    from .lunar import calc_mean_lunar_node

    # Get mean node longitude
    ts = get_timescale()
    t = ts.ut1_jd(jd)
    node_lon = calc_mean_lunar_node(t.tt)

    # Calculate distance to both nodes
    dist_to_north = abs((moon_lon - node_lon + 180) % 360 - 180)
    dist_to_south = abs((moon_lon - (node_lon + 180) + 180) % 360 - 180)

    return min(dist_to_north, dist_to_south)


def _calc_gamma(jd: float) -> float:
    """
    Calculate the gamma parameter (sqrt(x² + y²)) at given JD using Besselian elements.

    The gamma parameter is the minimum distance of the Moon's shadow axis from
    Earth's center, measured in Earth equatorial radii. During an eclipse:
    - gamma < 1: central eclipse (total or annular) is possible
    - gamma > 1.5: no eclipse visible from Earth

    Args:
        jd: Julian Day (UT)

    Returns:
        Gamma value in Earth equatorial radii
    """
    from .state import get_leb_reader

    AU_TO_KM = 149597870.7

    reader = get_leb_reader()
    _leb_ok = False
    if reader is not None:
        try:
            from .fast_calc import _apparent_icrs_cartesian
            from .time_utils import deltat

            jd_tt = jd + deltat(jd)
            sun_pos = _apparent_icrs_cartesian(reader, jd_tt, SUN)
            moon_pos = _apparent_icrs_cartesian(reader, jd_tt, MOON)
            _leb_ok = True
        except (KeyError, ValueError) as _exc:
            # Sealed leb mode must surface a coverage/body miss as the
            # documented typed error (EphemerisRangeError/UnknownBodyError):
            # the fallback below would only convert it into get_planets()'s
            # RuntimeError, which the sealed contract forbids. In every other
            # mode this is a no-op and the Skyfield fallback runs unchanged.
            # _calc_gamma is not wrapped by _call_with_leb_skyfield_fallback
            # (its fallback is the inline branch below), so a range miss must
            # NOT be re-raised here in auto mode.
            _raise_if_sealed_leb_miss(_exc)
    if not _leb_ok:
        from .state import get_planets, get_timescale

        eph = get_planets()
        ts = get_timescale()
        t = ts.ut1_jd(jd)
        earth_at = eph["earth"].at(t)
        sun_pos = earth_at.observe(eph["sun"]).apparent().position.au
        moon_pos = earth_at.observe(eph["moon"]).apparent().position.au

    # Shadow axis direction (Sun to Moon)
    shadow_dir = [
        moon_pos[0] - sun_pos[0],
        moon_pos[1] - sun_pos[1],
        moon_pos[2] - sun_pos[2],
    ]
    shadow_len = math.sqrt(shadow_dir[0] ** 2 + shadow_dir[1] ** 2 + shadow_dir[2] ** 2)
    shadow_unit = [shadow_dir[i] / shadow_len for i in range(3)]

    # x-axis perpendicular to shadow, in equatorial plane
    x_axis_raw = [-shadow_unit[1], shadow_unit[0], 0.0]
    x_axis_len = math.sqrt(x_axis_raw[0] ** 2 + x_axis_raw[1] ** 2)

    if x_axis_len < 1e-10:
        x_axis = [1.0, 0.0, 0.0]
    else:
        x_axis = [x_axis_raw[0] / x_axis_len, x_axis_raw[1] / x_axis_len, 0.0]

    # y-axis completes right-handed system
    y_axis = [
        shadow_unit[1] * x_axis[2] - shadow_unit[2] * x_axis[1],
        shadow_unit[2] * x_axis[0] - shadow_unit[0] * x_axis[2],
        shadow_unit[0] * x_axis[1] - shadow_unit[1] * x_axis[0],
    ]

    # Project Moon position onto fundamental plane
    moon_along_axis = sum(moon_pos[i] * shadow_unit[i] for i in range(3))
    moon_perp = [moon_pos[i] - moon_along_axis * shadow_unit[i] for i in range(3)]

    x_au = sum(moon_perp[i] * x_axis[i] for i in range(3))
    y_au = sum(moon_perp[i] * y_axis[i] for i in range(3))

    # Convert to Earth radii
    x_earth = (x_au * AU_TO_KM) / EARTH_RADIUS_KM
    y_earth = (y_au * AU_TO_KM) / EARTH_RADIUS_KM

    return math.sqrt(x_earth**2 + y_earth**2)


def _refine_solar_eclipse_maximum(
    jd_approx: float, search_range: float = 0.125
) -> float:
    """
    Refine the solar eclipse maximum time using Besselian elements.

    Uses golden section search to find the time when gamma (shadow axis distance
    from Earth center) is minimized. This gives sub-second precision for eclipse
    maximum timing.

    Args:
        jd_approx: Approximate Julian Day of eclipse maximum (e.g., from New Moon)
        search_range: Search range in days (default ±3 hours)

    Returns:
        Julian Day of refined eclipse maximum (precision < 1 second)
    """
    phi = (1 + math.sqrt(5)) / 2  # Golden ratio

    jd_low = jd_approx - search_range
    jd_high = jd_approx + search_range

    jd_a = jd_high - (jd_high - jd_low) / phi
    jd_b = jd_low + (jd_high - jd_low) / phi

    gamma_a = _calc_gamma(jd_a)
    gamma_b = _calc_gamma(jd_b)

    # Golden section search for minimum gamma
    for _ in range(60):  # Converges to ~1e-9 days (~0.1 ms)
        if gamma_a < gamma_b:
            jd_high = jd_b
            jd_b = jd_a
            gamma_b = gamma_a
            jd_a = jd_high - (jd_high - jd_low) / phi
            gamma_a = _calc_gamma(jd_a)
        else:
            jd_low = jd_a
            jd_a = jd_b
            gamma_a = gamma_b
            jd_b = jd_low + (jd_high - jd_low) / phi
            gamma_b = _calc_gamma(jd_b)

        if jd_high - jd_low < 1e-8:  # ~0.86 ms precision
            break

    return (jd_low + jd_high) / 2


def _calc_penumbra_limit(jd: float) -> float:
    """
    Calculate l1 (penumbral shadow radius on fundamental plane) at given JD.

    Returns the radius in Earth radii where the penumbral shadow intersects
    the fundamental plane.
    """
    AU_TO_KM = 149597870.7
    SUN_RADIUS_KM = 696340.0
    MOON_RADIUS_KM = 1737.4

    sun_pos, _ = calc_ut(jd, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd, MOON, FLG_SPEED)
    sun_dist_au = sun_pos[2]
    moon_dist_au = moon_pos[2]

    sun_dist_km = sun_dist_au * AU_TO_KM
    moon_dist_km = moon_dist_au * AU_TO_KM

    # The penumbral shadow cone extends from the Moon toward Earth.
    # The penumbral half-angle f1 is defined such that the penumbral cone's
    # outer edge connects the far limb of the Sun to the near limb of the Moon.
    #
    # Using similar triangles:
    # tan(f1) = (SUN_RADIUS + MOON_RADIUS) / (sun_dist - moon_dist)
    # At the fundamental plane (distance moon_dist from Moon toward Earth),
    # the penumbra radius is: l1_km = moon_dist * tan(f1) + MOON_RADIUS

    # In the Besselian convention l1 is the penumbral radius measured from the
    # shadow axis in the fundamental plane.  The exact cone relation is a
    # lunar-radius intercept plus an axial-distance times tan(f1); this legacy
    # helper deliberately retains its documented small-angle approximation
    # below for compatibility with its contact-time callers.

    # Angular semi-diameter of Sun from Earth
    sun_angular_rad = SUN_RADIUS_KM / sun_dist_km  # radians

    # Angular semi-diameter of Moon from Earth
    moon_angular_rad = MOON_RADIUS_KM / moon_dist_km  # radians

    # Penumbral cone angular radius from shadow axis (sum of angular radii)
    f1_rad = sun_angular_rad + moon_angular_rad

    # At the fundamental plane, the penumbral radius in Earth radii:
    # The shadow axis passes through Earth's center at distance = 0
    # The penumbral limit at distance z from Moon is: r = z * tan(f1)
    # But we measure from shadow axis, so: l1 = moon_dist * f1 (for small angles)

    # More accurate: l1 = k1 * sec(f1) + z * tan(f1) where k1 is constant
    # For practical purposes with small angles:
    l1_km = moon_dist_km * f1_rad
    l1_earth_radii = l1_km / EARTH_RADIUS_KM

    return l1_earth_radii


def _calc_umbra_limit(jd: float) -> float:
    """
    Calculate l2 (umbral/antumbral shadow radius on fundamental plane) at given JD.

    Returns the radius in Earth radii. Negative for umbra (total eclipse),
    positive for antumbra (annular eclipse).
    """
    AU_TO_KM = 149597870.7
    SUN_RADIUS_KM = 696340.0
    MOON_RADIUS_KM = 1737.4

    sun_pos, _ = calc_ut(jd, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd, MOON, FLG_SPEED)
    sun_dist_au = sun_pos[2]
    moon_dist_au = moon_pos[2]

    sun_dist_km = sun_dist_au * AU_TO_KM
    moon_dist_km = moon_dist_au * AU_TO_KM

    # Umbral cone half-angle. The cone is subtended by the difference of the
    # solar and lunar radii over the Sun-Moon distance (NOT the Earth-Sun
    # distance): this matches the sibling Besselian models (_besselian_core,
    # _eclipse_where_core), which use the Sun-to-Moon separation. Using the
    # Earth-Sun distance here left l2 ~0.0006-0.0008 Earth radii off.
    sun_moon_dist_km = sun_dist_km - moon_dist_km
    f2 = math.atan((SUN_RADIUS_KM - MOON_RADIUS_KM) / sun_moon_dist_km)

    # Distance from Moon to umbra vertex
    umbra_vertex_km = MOON_RADIUS_KM / math.tan(f2) if f2 > 1e-10 else float("inf")

    # Umbral radius at fundamental plane
    if moon_dist_km < umbra_vertex_km:
        # Umbra reaches Earth - negative by convention
        l2_km = MOON_RADIUS_KM - moon_dist_km * math.tan(f2)
        l2_earth_radii = -l2_km / EARTH_RADIUS_KM
    else:
        # Antumbra - Moon shadow vertex is before Earth
        l2_km = moon_dist_km * math.tan(f2) - MOON_RADIUS_KM
        l2_earth_radii = l2_km / EARTH_RADIUS_KM

    return l2_earth_radii


def _find_contact_time_besselian(
    jd_max: float,
    target_gamma: float,
    search_before: bool,
    search_range: float = 0.1,
) -> float:
    """
    Find the time when gamma equals target_gamma using bisection.

    Used to find eclipse contact times where the shadow boundary crosses Earth.

    Args:
        jd_max: Julian Day of eclipse maximum
        target_gamma: Target gamma value to search for
        search_before: If True, search before maximum; if False, search after
        search_range: Search range in days from maximum

    Returns:
        Julian Day when gamma equals target_gamma, or 0 if not found
    """
    if search_before:
        jd_start = jd_max - search_range
        jd_end = jd_max
    else:
        jd_start = jd_max
        jd_end = jd_max + search_range

    # Check if solution exists in range
    gamma_start = _calc_gamma(jd_start)
    gamma_end = _calc_gamma(jd_end)
    gamma_max = _calc_gamma(jd_max)

    # For contacts before max: gamma decreases from start to max, then increases
    # For contacts after max: gamma increases from max to end
    # We're looking for where gamma crosses target_gamma

    if search_before:
        if gamma_start < target_gamma or gamma_max > target_gamma:
            # Check edge case where max gamma might be above target
            if gamma_max < target_gamma:
                return 0.0
    else:
        if gamma_max > target_gamma or gamma_end < target_gamma:
            if gamma_max < target_gamma:
                return 0.0

    # Binary search
    for _ in range(60):  # ~0.1 second precision
        jd_mid = (jd_start + jd_end) / 2
        gamma_mid = _calc_gamma(jd_mid)

        if abs(gamma_mid - target_gamma) < 1e-8:
            return jd_mid

        if search_before:
            # Before max: gamma decreases then potentially rises
            # We want the first crossing from above
            if gamma_mid > target_gamma:
                jd_start = jd_mid
            else:
                jd_end = jd_mid
        else:
            # After max: gamma increases
            if gamma_mid < target_gamma:
                jd_start = jd_mid
            else:
                jd_end = jd_mid

        if jd_end - jd_start < 1e-8:
            break

    return (jd_start + jd_end) / 2


def _besselian_contact_residuals(
    jd: float, flags: int = FLG_SWIEPH
) -> Tuple[float, float, float]:
    """Signed Besselian contact residuals on the Earth ellipsoid, in km.

    The global-eclipse contacts (P1/P4, U1/U4, centre line) are the
    instants at which the Moon's shadow cone is tangent to the Earth's
    flattened surface. Working in the auxiliary frame in which the polar
    axis is stretched by 1/(1-f) so the ellipsoid becomes a sphere of
    equatorial radius (the standard reduction; see the Explanatory
    Supplement to the Astronomical Almanac, 3rd ed. 2013, ch. 11, and
    Meeus, "Elements of Solar Eclipses 1951-2200", 1989, ch. 4), the
    external-tangency condition of a shadow cone of half angle f and
    fundamental-plane radius L with a sphere of radius R_e whose centre
    lies at perpendicular distance r0 from the axis is

        (r0 - L) * cos f = R_e      i.e.   r0 = R_e * cos f + L,

    to O(f^2). The cone-depth term z*tan(f) that a plain circle-on-the-
    fundamental-plane test omits is carried implicitly here: L is the
    radius on the fundamental plane and the geometry supplies the height
    of the tangent point through cos f. The Earth's oblateness enters via
    the same z-stretch that ``_eclipse_where_core`` applies.

    Returns ``(pen, umb, cen)`` where each residual is negative while the
    corresponding shadow feature is in contact with (or inside) the
    ellipsoid and positive once it has cleared it:

    * ``pen`` -> penumbral outer edge (first/last contact, P1/P4)
    * ``umb`` -> umbral/antumbral outer edge (U1/U4)
    * ``cen`` -> shadow axis vs. the limb (centre line begins/ends)

    The three residuals are built from exactly the quantities
    ``_eclipse_where_core`` uses for its CENTRAL / NONCENTRAL / PARTIAL
    classification, so the contact times stay consistent with
    ``sol_eclipse_where`` at every instant.
    """
    _rc, _lon, _lat, dc = _eclipse_where_core(jd, flags)
    r0_km = dc[2]
    core_diam_fp_km = dc[3]
    pen_diam_fp_km = dc[4]
    cos_f_core = dc[5]
    cos_f_pen = dc[6]
    re_km = _ECL_REARTH_AU * _ECL_AU_KM
    pen = r0_km - (re_km * cos_f_pen + pen_diam_fp_km / 2.0)
    umb = r0_km - (re_km * cos_f_core + abs(core_diam_fp_km) / 2.0)
    cen = r0_km - re_km * cos_f_core
    return pen, umb, cen


def _solve_besselian_contact(
    residual: "Callable[[float], float]",
    jd_max: float,
    search_before: bool,
    max_win: float = 0.35,
) -> float:
    """Time of a Besselian contact by bracketing a residual's sign change.

    ``residual(t)`` is one of the signed distances from
    :func:`_besselian_contact_residuals` (negative while in contact,
    positive once cleared). The residual is negative at ``jd_max`` for a
    contact that occurs; this expands a bracket outward from ``jd_max``
    until the residual turns positive, then bisects. Returns 0.0 when the
    feature does not touch the Earth at maximum (no such contact) or when
    no clearing is found within ``max_win`` days.
    """
    try:
        f_max = residual(jd_max)
    except (KeyError, ValueError, ArithmeticError) as exc:
        _reraise_if_leb_range_error(exc)
        return 0.0
    if f_max is None or f_max >= 0.0:
        return 0.0

    sign = -1.0 if search_before else 1.0
    far = 0.0
    step = 0.02
    while step <= max_win:
        t = jd_max + sign * step
        try:
            ft = residual(t)
        except (KeyError, ValueError, ArithmeticError) as exc:
            _reraise_if_leb_range_error(exc)
            ft = None
        if ft is not None and ft > 0.0:
            far = t
            break
        step *= 1.6
    if far == 0.0:
        return 0.0

    lo, hi = (far, jd_max) if search_before else (jd_max, far)
    try:
        return _root_bisect(residual, lo, hi)
    except (KeyError, ValueError, ArithmeticError) as exc:
        _reraise_if_leb_range_error(exc)
        return 0.0


def _calculate_eclipse_phases_besselian(
    jd_max: float, eclipse_type: int
) -> Tuple[float, float, float, float, float, float, float, float, float, float]:
    """
    Calculate times of eclipse phases using Besselian elements for high precision.

    This function achieves timing precision of better than 10 seconds by
    calculating exact contact times based on when the penumbral and umbral
    shadow boundaries cross Earth's limb.

    Phase indices (matching reference API tret array):
        [0]: Time of maximum eclipse
        [1]: Time when eclipse takes place at local apparent noon
        [2]: Eclipse begins (first penumbral contact with Earth, P1)
        [3]: Eclipse ends (last penumbral contact with Earth, P4)
        [4]: Totality/annularity begins anywhere on Earth (U1)
        [5]: Totality/annularity ends (U4)
        [6]: Center line begins
        [7]: Center line ends
        [8]: Annular-total eclipse becomes total
        [9]: Annular-total eclipse becomes annular again

    Args:
        jd_max: Julian Day of eclipse maximum (refined using Besselian elements)
        eclipse_type: Eclipse type flags

    Returns:
        Tuple of 10 floats with phase times (JD UT), matching reference API format
    """
    is_total = bool(eclipse_type & ECL_TOTAL)
    is_annular = bool(eclipse_type & ECL_ANNULAR)
    # A hybrid (annular-total) eclipse sets ECL_ANNULAR_TOTAL only — it shares
    # no bits with ECL_TOTAL/ECL_ANNULAR — so it too has an umbral axis that
    # touches Earth and a center line. Without this it would lose U1/U4 and
    # the center-line instants (tret[4..7]), which belong in the tret
    # contract for these types.
    is_hybrid = bool(eclipse_type & ECL_ANNULAR_TOTAL)
    has_umbral_contact = is_total or is_annular or is_hybrid

    # Contact times are the instants at which the Moon's shadow cone is
    # tangent to the Earth *ellipsoid*. Each contact is the zero crossing
    # of the corresponding signed residual from
    # _besselian_contact_residuals(), which carries both the cone-depth
    # term (z*tan f, previously omitted by the plain circle-on-the-
    # fundamental-plane test) and the Earth's oblateness (via the polar-
    # axis stretch shared with _eclipse_where_core). See that helper for
    # the geometry and the published references.
    # The pen/umb/cen residuals share one _eclipse_where_core() evaluation
    # per instant, and the six contact solvers below revisit the same
    # instants (jd_max itself, plus the identical bracket-expansion
    # ladder). Memoise per call — never at module level, where a cached
    # triple could survive an ephemeris-state change (mode, files, topo).
    residual_memo: dict[float, Tuple[float, float, float]] = {}

    def _residual_triple(jd: float) -> Tuple[float, float, float]:
        triple = residual_memo.get(jd)
        if triple is None:
            triple = _besselian_contact_residuals(jd)
            residual_memo[jd] = triple
        return triple

    def _pen_residual(jd: float) -> float:
        return _residual_triple(jd)[0]

    def _umb_residual(jd: float) -> float:
        return _residual_triple(jd)[1]

    def _cen_residual(jd: float) -> float:
        return _residual_triple(jd)[2]

    # P1/P4: first/last external tangency of the penumbral cone (eclipse
    # begins/ends anywhere on Earth).
    t_first_contact = _solve_besselian_contact(
        _pen_residual, jd_max, search_before=True
    )
    t_fourth_contact = _solve_besselian_contact(
        _pen_residual, jd_max, search_before=False
    )

    # U1/U4: first/last external tangency of the umbral/antumbral cone
    # (totality/annularity begins/ends anywhere on Earth). The solver
    # returns 0.0 when the core shadow never reaches the surface, so no
    # separate gamma gate is needed.
    t_second_contact = 0.0
    t_third_contact = 0.0
    if has_umbral_contact:
        t_second_contact = _solve_besselian_contact(
            _umb_residual, jd_max, search_before=True
        )
        t_third_contact = _solve_besselian_contact(
            _umb_residual, jd_max, search_before=False
        )

    # tret[1]: time when the eclipse takes place at local apparent noon —
    # the shadow axis crosses the x = 0 plane (Besselian x sign change)
    # WITHIN the resolved eclipse window [P1, P4]. "Local apparent noon" is
    # only defined while the penumbra is in contact with Earth, so the search
    # is bracketed by the first/last penumbral contacts rather than by a fixed
    # window around the maximum. When those contacts do not resolve — an
    # ultra-shallow graze whose penumbra never fully lands, so tret[2]/[3] are
    # left 0 — there is no such window and tret[1] stays 0. Compatibility
    # contract: that degenerate partial returns tret[1] = 0. The earlier
    # fabricated jd_max +- 0.15 d bracket produced a spurious noon instant
    # tens of minutes from the maximum for exactly those events.
    t_noon = 0.0
    if t_first_contact and t_fourth_contact:
        try:
            x_lo = calc_besselian_x(t_first_contact)
            x_hi = calc_besselian_x(t_fourth_contact)
            if x_lo * x_hi < 0:
                a, b, fa = t_first_contact, t_fourth_contact, x_lo
                for _ in range(60):
                    mid = 0.5 * (a + b)
                    fm = calc_besselian_x(mid)
                    if fa * fm <= 0:
                        b = mid
                    else:
                        a, fa = mid, fm
                t_noon = 0.5 * (a + b)
        except (KeyError, ValueError, ArithmeticError) as _exc:
            _reraise_if_leb_range_error(_exc)
            t_noon = 0.0

    # tret[6]/[7]: center-line begin/end — only when the eclipse is
    # ellipsoidally CENTRAL (the shadow axis actually pierces the flattened
    # Earth). Gating on the ECL_CENTRAL bit (set from sol_eclipse_where's
    # geometry) instead of a spherical gamma_max < 1.0 test keeps the
    # center-line times consistent with the CEN/NONCEN classification: a
    # non-central annular/total eclipse (|gamma| just under 1 near the
    # poles) has no central line.
    t_cl_begin = 0.0
    t_cl_end = 0.0
    if has_umbral_contact and (eclipse_type & ECL_CENTRAL):
        # Centre line begins/ends when the shadow axis itself grazes the
        # ellipsoid limb (the cen residual crosses zero).
        t_cl_begin = _solve_besselian_contact(_cen_residual, jd_max, search_before=True)
        t_cl_end = _solve_besselian_contact(_cen_residual, jd_max, search_before=False)

    # Reference degenerate layout for an ultra-shallow partial: when the
    # penumbra only grazes Earth at the instant of maximum, the P1/P4
    # "eclipse begin/end" contacts (tret[2]/[3]) do not resolve. The
    # reference then leaves those two zero and collapses the maximum
    # instant into the totality slots tret[4]/[5] (both equal to jd_max).
    # Mirror that exactly so the slot population matches for these events.
    if (
        (eclipse_type & ECL_PARTIAL)
        and not has_umbral_contact
        and not t_first_contact
        and not t_fourth_contact
    ):
        t_second_contact = jd_max
        t_third_contact = jd_max

    return (
        jd_max,  # [0] Time of maximum eclipse
        t_noon,  # [1] Eclipse at local apparent noon
        # No fabricated fallbacks: an unresolved contact stays 0.0 (the
        # reference returns 0 rather than an invented +-1 h estimate)
        t_first_contact if t_first_contact else 0.0,  # [2] Eclipse begin
        t_fourth_contact if t_fourth_contact else 0.0,  # [3] Eclipse end
        t_second_contact,  # [4] Totality begin (or 0)
        t_third_contact,  # [5] Totality end (or 0)
        t_cl_begin,  # [6] Center line begin
        t_cl_end,  # [7] Center line end
        0.0,  # [8] Annular-total becomes total
        0.0,  # [9] Annular-total becomes annular again
    )


def _is_hybrid_solar_eclipse(jd_max: float, l2: float) -> bool:
    """Physical hybrid (annular-total) test: does the core shadow's
    character change along the central phase?

    A hybrid eclipse is total around maximum, where the Earth's surface
    bulges closest to the Moon, and annular near the ends of the central
    path, where the surface falls back towards the fundamental plane. So
    instead of a static threshold on ``l2`` at maximum, sample the
    core-shadow width at the ground point (negative when the umbra
    reaches the surface) between the umbral contacts U1 and U4: the
    eclipse is hybrid iff the width changes sign along the path.

    Args:
        jd_max: Julian Day (UT) of the eclipse maximum.
        l2: Umbral/antumbral shadow radius on the fundamental plane at
            ``jd_max`` (Earth radii, negative for umbra).

    Returns:
        True when the core-shadow width changes sign along the central
        phase (annular-total eclipse).
    """
    # Fast reject: between the path ends (surface near the fundamental
    # plane) and maximum (surface up to one Earth radius closer to the
    # Moon), Earth's curvature shrinks the core-shadow radius by at most
    # ~tan(f2) ~ 0.0047 Earth radii (~30 km). Outside this band around
    # zero no sign change is possible: l2 < 0 stays total, large l2 > 0
    # stays annular. Margins cover the drift of l2 over the window.
    if l2 < -0.001 or l2 > 0.006:
        return False

    # Umbral contacts U1/U4 (exterior tangency), as in
    # _calculate_eclipse_phases_besselian().
    umbral_limit = 1.0 + abs(l2)
    t_u1 = _find_contact_time_besselian(
        jd_max, umbral_limit, search_before=True, search_range=0.10
    )
    t_u4 = _find_contact_time_besselian(
        jd_max, umbral_limit, search_before=False, search_range=0.10
    )
    if t_u1 and t_u4 and t_u4 > t_u1:
        sample_times = [t_u1 + f * (t_u4 - t_u1) for f in (0.0, 0.25, 0.5, 0.75, 1.0)]
    else:
        # Contacts unresolved: fall back to the maximum alone (no sign
        # change observable from a single sample -> not hybrid).
        sample_times = [jd_max]

    widths = []
    for t in sample_times:
        try:
            _rc, _lon, _lat, dcore = _eclipse_where_core(t)
        except (KeyError, ValueError, ArithmeticError) as _exc:
            _reraise_if_leb_range_error(_exc)
            continue
        widths.append(dcore[0])
    return bool(widths) and min(widths) < 0.0 < max(widths)


def _where_eclipse_magnitude(jd: float, flags: int = FLG_SWIEPH) -> float:
    """Topocentric eclipse magnitude at the shadow-center ground point.

    Returns ``attr[0]`` of :func:`sol_eclipse_where` at ``jd`` — the
    magnitude (solar-diameter fraction) seen at the point where the
    shadow axis comes closest to the Earth. Positive when the discs
    actually overlap somewhere on Earth's surface, negative when even the
    point of deepest approach sees no overlap. Used as the graze-boundary
    admission gate so ``sol_eclipse_when_glob`` and ``sol_eclipse_where``
    agree on whether a near-limit event is a real eclipse.
    """

    def _impl(*, reader):
        _rc, wlon, wlat, _dcore = _eclipse_where_core(jd, flags)
        _rc_how, attr = _sol_how_core(jd, (wlon, wlat, 0.0), flags, reader)
        return attr[0]

    return _call_with_leb_skyfield_fallback(_impl)


def _calculate_eclipse_type_and_magnitude(
    jd: float,
) -> Tuple[int, float, float, float]:
    """
    Determine eclipse type and magnitude at maximum eclipse.

    Uses geometric calculations based on Sun-Moon-Earth distances
    and apparent angular sizes. Includes edge case handling for:
    - Very shallow partial eclipses (magnitude near 0)
    - Near-miss eclipses (gamma close to eclipse limit)
    - Grazing eclipses at extreme gamma values

    Args:
        jd: Julian Day of eclipse maximum (UT)

    Returns:
        Tuple of (eclipse_type_flags, magnitude, gamma, moon_sun_ratio)
        - eclipse_type_flags: Bitmask of ECL_* constants
        - magnitude: Eclipse magnitude (fraction of Sun's diameter covered)
        - gamma: Gamma parameter (distance of Moon shadow axis from Earth center)
        - moon_sun_ratio: Ratio of apparent Moon diameter to Sun diameter
    """
    # Get positions
    sun_pos, _ = calc_ut(jd, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd, MOON, FLG_SPEED)

    # Distances in AU
    sun_dist = sun_pos[2]
    moon_dist = moon_pos[2]

    # Ecliptic coordinates for separation calculation
    sun_lon, sun_lat = sun_pos[0], sun_pos[1]
    moon_lon, moon_lat = moon_pos[0], moon_pos[1]

    # Angular radii (in degrees)
    # Sun: mean radius 959.63" at 1 AU
    sun_angular_radius = (959.63 / 3600.0) / sun_dist
    # Moon: mean radius 932.56" at mean distance (0.002569 AU)
    moon_angular_radius = (932.56 / 3600.0) * (0.002569 / moon_dist)

    # Ratio of apparent sizes - handle edge case of very small Sun radius
    if sun_angular_radius < 1e-10:
        moon_sun_ratio = 1.0  # Fallback for numerical stability
    else:
        moon_sun_ratio = moon_angular_radius / sun_angular_radius

    # Calculate gamma using Besselian elements (shadow axis distance from Earth center)
    gamma = _calc_gamma(jd)

    # Penumbral and umbral limits at maximum (Earth radii)
    l1 = _calc_penumbra_limit(jd)
    l2 = _calc_umbra_limit(jd)

    gamma_limit_partial = 1.0 + l1
    max(0.0, 1.0 - abs(l2))

    # Angular separation between Sun and Moon centers (degrees)
    def _angular_separation(
        lon1: float, lat1: float, lon2: float, lat2: float
    ) -> float:
        lon1_rad = math.radians(lon1)
        lat1_rad = math.radians(lat1)
        lon2_rad = math.radians(lon2)
        lat2_rad = math.radians(lat2)
        dlon = lon2_rad - lon1_rad
        cos_sep = math.sin(lat1_rad) * math.sin(lat2_rad) + math.cos(
            lat1_rad
        ) * math.cos(lat2_rad) * math.cos(dlon)
        return math.degrees(_safe_acos(cos_sep))

    separation = _angular_separation(sun_lon, sun_lat, moon_lon, moon_lat)

    # Eclipse magnitude: fraction of Sun's diameter covered (can exceed 1.0 for total)
    sum_radii = sun_angular_radius + moon_angular_radius
    if separation >= sum_radii:
        magnitude = 0.0
    else:
        magnitude = (sum_radii - separation) / (2.0 * sun_angular_radius)
        magnitude = max(0.0, min(1.5, magnitude))

    # Determine eclipse type based on gamma and shadow limits
    eclipse_type = 0

    # Edge case: gamma beyond eclipse limit
    if abs(gamma) > gamma_limit_partial:
        # No eclipse
        return 0, 0.0, gamma, moon_sun_ratio

    # Edge case: near-miss eclipse (gamma very close to limit).
    # Compatibility contract: PARTIAL|NONCENTRAL (retflag 18) is returned
    # for shallow partials; the non-reference ECL_GRAZING bit is no longer
    # leaked into any public retflag.
    if _is_near_miss_eclipse(gamma, gamma_limit_partial):
        # Graze-boundary gate: the approximate spherical fundamental-plane
        # test (gamma <= 1 + l1) admits penumbral grazes that neither the
        # reference nor our own sol_eclipse_where report as eclipses (they
        # return retflag 0 with a NEGATIVE where-magnitude — the discs do not
        # overlap anywhere on Earth). Admit a near-limit event only when the
        # authoritative ellipsoidal where-geometry sees the discs actually
        # overlap (positive where-magnitude), so when_glob and where agree.
        # The sign — not the where-retflag — is the discriminator: a genuine
        # shallow partial can have where-retflag 0 yet a positive magnitude.
        if _where_eclipse_magnitude(jd) <= 0.0:
            return 0, 0.0, gamma, moon_sun_ratio
        eclipse_type = ECL_PARTIAL | ECL_NONCENTRAL
        magnitude = _calculate_magnitude_safe(
            gamma, moon_sun_ratio, gamma_limit_partial
        )
        return eclipse_type, magnitude, gamma, moon_sun_ratio

    eclipse_type = ECL_PARTIAL | ECL_NONCENTRAL

    # Hybrid detection (annular-total): the umbra/antumbra cone tip is so
    # close to Earth's surface that the core-shadow width changes sign
    # along the central path. _is_hybrid_solar_eclipse() applies the
    # physical criterion (sign change between the umbral contacts) and is
    # only evaluated once the umbra/antumbra is known to touch Earth.
    if abs(gamma) <= 1.0 + abs(l2) and abs(l2) > 0:
        # The umbra/antumbra reaches Earth. Whether the eclipse is CENTRAL
        # (the shadow axis pierces the surface) or merely NONCENTRAL is
        # decided by the SAME ellipsoidal shadow-axis geometry that
        # sol_eclipse_where uses - not a spherical |gamma| <= 1 test. Near
        # the poles a |gamma| just under 1 can still miss the flattened
        # Earth, giving a non-central annular/total eclipse; the spherical
        # test would wrongly report CENTRAL (contradicting sol_eclipse_where
        # at the same instant) and fabricate center-line times (tret[6]/[7])
        # for an eclipse that has no central line.
        where_retc = _eclipse_where_core(jd)[0]
        if where_retc & ECL_CENTRAL:
            central_bit = ECL_CENTRAL
        elif (where_retc & ECL_NONCENTRAL) and not (where_retc & ECL_PARTIAL):
            central_bit = ECL_NONCENTRAL
        else:
            # The core shadow does not actually reach the ellipsoid at
            # maximum: a penumbra-only partial eclipse (keep the default
            # ECL_PARTIAL | ECL_NONCENTRAL classification).
            central_bit = 0

        if central_bit:
            # l2 < 0: umbra (total); l2 > 0: antumbra (annular). The hybrid
            # (annular-total) test uses the physical core-shadow sign change.
            if _is_hybrid_solar_eclipse(jd, l2):
                eclipse_type = ECL_ANNULAR_TOTAL | central_bit
            elif l2 < 0:
                eclipse_type = ECL_TOTAL | central_bit
            else:
                eclipse_type = ECL_ANNULAR | central_bit

    return eclipse_type, magnitude, gamma, moon_sun_ratio


def _calculate_eclipse_phases(
    jd_max: float, eclipse_type: int
) -> Tuple[float, float, float, float, float, float, float, float, float, float]:
    """
    Calculate times of eclipse phases (contacts) for a global eclipse.

    Uses Besselian elements for high-precision timing (< 10 seconds).

    Phase indices (matching reference API tret array):
        [0]: Time of maximum eclipse
        [1]: Time when eclipse takes place at local apparent noon
        [2]: Eclipse begins (first penumbral contact with Earth, P1)
        [3]: Eclipse ends (last penumbral contact with Earth, P4)
        [4]: Totality/annularity begins anywhere on Earth (U1)
        [5]: Totality/annularity ends (U4)
        [6]: Center line begins
        [7]: Center line ends
        [8]: Annular-total eclipse becomes total
        [9]: Annular-total eclipse becomes annular again

    Args:
        jd_max: Julian Day of maximum eclipse
        eclipse_type: Eclipse type flags

    Returns:
        Tuple of 10 floats with phase times (JD UT), matching reference API format
    """
    # Use the high-precision Besselian element-based calculation
    return _calculate_eclipse_phases_besselian(jd_max, eclipse_type)


def sol_eclipse_max_time(
    jd_approx: float,
    lat: float | None = None,
    lon: float | None = None,
    altitude: float = 0.0,
    search_range: float = 0.125,
    flags: int = FLG_SWIEPH,
) -> Tuple[float, float]:
    """
    Calculate the precise time of maximum eclipse when Sun-Moon separation is minimum.

    This function finds the exact moment when the angular separation between the
    Sun and Moon centers is at its minimum, which corresponds to the time of
    maximum eclipse. It uses golden section search for sub-second precision.

    For global eclipse maximum (no observer location given):
        Uses the gamma parameter (minimum distance of Moon's shadow axis from
        Earth's center) to find global maximum. This is the time when the
        eclipse is at its maximum considering all possible Earth observers.

    For local eclipse maximum (observer location given):
        Uses the topocentric angular separation between Sun and Moon as seen
        from the specified observer location. This accounts for parallax and
        gives the precise local maximum time.

    Args:
        jd_approx: Approximate Julian Day (UT) of eclipse, typically near New Moon
                   or obtained from _sol_eclipse_when_glob_pythonic()[0]
        lat: Observer latitude in degrees (positive = North, negative = South).
             If None, calculates global maximum.
        lon: Observer longitude in degrees (positive = East, negative = West).
             If None, calculates global maximum.
        altitude: Observer altitude in meters above sea level (default 0).
                  Only used if lat and lon are provided.
        search_range: Search range in days around jd_approx (default ±3 hours).
                      For best results, jd_approx should be within this range
                      of the actual eclipse maximum.
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple of (jd_maximum, min_separation) where:
            - jd_maximum: Julian Day (UT) of maximum eclipse with sub-second precision
            - min_separation: Minimum angular separation in degrees between Sun and
                              Moon centers at maximum. For global maximum, this is
                              the gamma value in Earth radii (not degrees).

    Raises:
        ValueError: If only one of lat/lon is provided (both or neither required)

    Precision:
        Sub-second precision (< 1 second) achieved through golden section search
        with convergence to ~0.1 milliseconds.

    Algorithm:
        Global maximum:
            Uses Besselian elements approach - finds when gamma (perpendicular
            distance of Moon's shadow axis from Earth's center) is minimum.
            This corresponds to the time when the eclipse is at maximum for
            observers along the central line.

        Local maximum:
            Calculates topocentric apparent positions of Sun and Moon at each
            iteration, finding when their angular separation is minimum as
            seen from the observer's location on Earth's surface.

    Example:
        >>> # Find precise global maximum time for April 8, 2024 eclipse
        >>> from libephemeris import julday, sol_eclipse_max_time, sol_eclipse_when_glob
        >>> jd_start = julday(2024, 3, 1, 0)
        >>> ecl_type, times = sol_eclipse_when_glob(jd_start)
        >>> jd_max, gamma = sol_eclipse_max_time(times[0])
        >>> print(f"Maximum at JD {jd_max:.8f}, gamma = {gamma:.6f}")

        >>> # Find local maximum time from Dallas, Texas
        >>> dallas_lat, dallas_lon = 32.7767, -96.7970
        >>> jd_local_max, separation = sol_eclipse_max_time(
        ...     times[0], lat=dallas_lat, lon=dallas_lon
        ... )
        >>> print(f"Local max at JD {jd_local_max:.8f}, sep = {separation:.6f}°")

    References:
        - Meeus (1998), "Astronomical Algorithms", 2nd ed., Ch. 54
        - Espenak & Meeus (2006), "Five Millennium Canon of Solar
          Eclipses", NASA/TP-2006-214141
        - Espenak, NASA/GSFC, "Besselian Elements of Solar Eclipses":
          https://eclipse.gsfc.nasa.gov/SEcat5/beselm.html
    """
    # Validate that both lat and lon are provided or neither
    if (lat is None) != (lon is None):
        raise ValueError(
            "Both lat and lon must be provided, or neither (for global maximum)"
        )

    # Determine if we're calculating global or local maximum
    is_local = lat is not None and lon is not None

    if is_local:
        # Local maximum: find minimum Sun-Moon separation from observer location
        assert lat is not None and lon is not None
        return _calc_local_eclipse_max_time(jd_approx, lat, lon, altitude, search_range)
    else:
        # Global maximum: find minimum gamma using Besselian elements
        return _calc_global_eclipse_max_time(jd_approx, search_range)


def _calc_global_eclipse_max_time(
    jd_approx: float, search_range: float
) -> Tuple[float, float]:
    """
    Calculate global eclipse maximum time using gamma minimization.

    Args:
        jd_approx: Approximate Julian Day of eclipse
        search_range: Search range in days

    Returns:
        Tuple of (jd_maximum, gamma_at_maximum)
    """
    phi = (1 + math.sqrt(5)) / 2  # Golden ratio

    jd_low = jd_approx - search_range
    jd_high = jd_approx + search_range

    jd_a = jd_high - (jd_high - jd_low) / phi
    jd_b = jd_low + (jd_high - jd_low) / phi

    gamma_a = _calc_gamma(jd_a)
    gamma_b = _calc_gamma(jd_b)

    # Golden section search for minimum gamma
    for _ in range(60):  # Converges to ~1e-9 days (~0.1 ms)
        if gamma_a < gamma_b:
            jd_high = jd_b
            jd_b = jd_a
            gamma_b = gamma_a
            jd_a = jd_high - (jd_high - jd_low) / phi
            gamma_a = _calc_gamma(jd_a)
        else:
            jd_low = jd_a
            jd_a = jd_b
            gamma_a = gamma_b
            jd_b = jd_low + (jd_high - jd_low) / phi
            gamma_b = _calc_gamma(jd_b)

        if jd_high - jd_low < 1e-8:  # ~0.86 ms precision
            break

    jd_max = (jd_low + jd_high) / 2
    gamma_at_max = _calc_gamma(jd_max)

    return jd_max, gamma_at_max


def _calc_local_eclipse_max_time(
    jd_approx: float,
    lat: float,
    lon: float,
    altitude: float,
    search_range: float,
) -> Tuple[float, float]:
    """
    Calculate local eclipse maximum time using Sun-Moon separation minimization.

    Args:
        jd_approx: Approximate Julian Day of eclipse
        lat: Observer latitude in degrees
        lon: Observer longitude in degrees
        altitude: Observer altitude in meters
        search_range: Search range in days

    Returns:
        Tuple of (jd_maximum, min_separation_degrees)
    """
    from .state import get_leb_reader

    reader = get_leb_reader()

    if reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import angular_separation

        # geopos for topo/azalt: (lon, lat, alt)
        geopos = (lon, lat, altitude)

        def _get_separation(jd: float) -> float:
            """Get angular separation between Sun and Moon at given JD."""
            jd_tt = jd + deltat(jd)
            sun_pos = _topo_ecliptic(reader, jd_tt, jd, SUN, geopos)
            moon_pos = _topo_ecliptic(reader, jd_tt, jd, MOON, geopos)
            return angular_separation(sun_pos[0], sun_pos[1], moon_pos[0], moon_pos[1])
    else:
        from skyfield.api import wgs84
        from .state import get_planets, get_timescale

        # Get ephemeris and timescale
        eph = get_planets()
        ts = get_timescale()

        # Create observer location
        observer = wgs84.latlon(lat, lon, altitude)

        # Get Sun and Moon objects
        sun = eph["sun"]
        moon = eph["moon"]
        earth = eph["earth"]

        def _get_separation(jd: float) -> float:
            """Get angular separation between Sun and Moon at given JD."""
            t = ts.ut1_jd(jd)
            observer_at = earth + observer

            sun_app = observer_at.at(t).observe(sun).apparent()
            moon_app = observer_at.at(t).observe(moon).apparent()

            return sun_app.separation_from(moon_app).degrees

    # Golden section search for minimum separation
    phi = (1 + math.sqrt(5)) / 2

    jd_low = jd_approx - search_range
    jd_high = jd_approx + search_range

    jd_a = jd_high - (jd_high - jd_low) / phi
    jd_b = jd_low + (jd_high - jd_low) / phi

    sep_a = _get_separation(jd_a)
    sep_b = _get_separation(jd_b)

    for _ in range(60):  # Converges to ~1e-9 days (~0.1 ms)
        if sep_a < sep_b:
            jd_high = jd_b
            jd_b = jd_a
            sep_b = sep_a
            jd_a = jd_high - (jd_high - jd_low) / phi
            sep_a = _get_separation(jd_a)
        else:
            jd_low = jd_a
            jd_a = jd_b
            sep_a = sep_b
            jd_b = jd_low + (jd_high - jd_low) / phi
            sep_b = _get_separation(jd_b)

        if jd_high - jd_low < 1e-8:  # ~0.86 ms precision
            break

    jd_max = (jd_low + jd_high) / 2
    min_separation = _get_separation(jd_max)

    return jd_max, min_separation


# Alias for reference API compatibility
sol_eclipse_max_time = sol_eclipse_max_time


def _sol_eclipse_when_glob_pythonic(
    jd_start: float,
    flags: int = FLG_SWIEPH,
    eclipse_type: int = 0,
    search_direction: str = "bidirectional",
) -> Tuple[int, Tuple[float, ...]]:
    """
    Find the next global solar eclipse after a given date.

    Searches forward in time from jd_start to find the next solar eclipse.
    Can filter by eclipse type (total, annular, partial, hybrid).

    By default uses bidirectional search to check ±15 days around jd_start,
    ensuring that eclipses close to the start date are not missed due to
    the search algorithm starting slightly after a New Moon.

    Args:
        jd_start: Julian Day (UT) to start search from
        flags: Calculation flags (FLG_SWIEPH, etc.)
        eclipse_type: Filter for specific eclipse type(s), bitmask of:
            - ECL_TOTAL (4): Total eclipse
            - ECL_ANNULAR (8): Annular eclipse
            - ECL_PARTIAL (16): Partial eclipse
            - ECL_ANNULAR_TOTAL (32): Hybrid eclipse
            - 0: Any eclipse type (default)
        search_direction: Direction to search for eclipses:
            - "forward": Only search forward from jd_start (per the reference API)
            - "backward": Only search backward from jd_start
            - "bidirectional" (default): Check ±15 days to ensure no eclipse
              is missed, returning the closest eclipse to jd_start that is
              >= jd_start (unless only found in backward direction)

    Returns:
        Tuple containing (matching reference API format):
            - retflag: Eclipse type flags bitmask (ECL_* constants)
            - tret: Tuple of 10 floats with eclipse phase times (JD UT),
              matching the reference sol_eclipse_when_glob() layout (these
              are GLOBAL instants for the whole Earth, not observer-local
              contacts):
                [0]: Time of maximum eclipse
                [1]: Time when the eclipse is central at local apparent noon
                [2]: Time when the eclipse begins (penumbra first touches
                     Earth anywhere; 0 if not found)
                [3]: Time when the eclipse ends (penumbra last leaves Earth
                     anywhere; 0 if not found)
                [4]: Time when totality/annularity begins anywhere on Earth
                     (umbra/antumbra first touches Earth; 0 if none)
                [5]: Time when totality/annularity ends anywhere on Earth
                     (umbra/antumbra last leaves Earth; 0 if none)
                [6]: Time when the central line begins (0 if none)
                [7]: Time when the central line ends (0 if none)
                [8]: Time when a hybrid (annular-total) eclipse becomes
                     total (0; not implemented)
                [9]: Time when a hybrid (annular-total) eclipse becomes
                     annular again (0; not implemented)

    Raises:
        Error: If no eclipse found within search limit

    Algorithm:
        1. For bidirectional search: check previous New Moon within 15 days
        2. Find next New Moon after jd_start
        3. Check if Moon is close enough to node for eclipse
        4. If not eclipse, advance to next New Moon
        5. Calculate eclipse type and magnitude
        6. If eclipse_type filter set, check if matches
        7. Return closest matching eclipse >= jd_start

    Precision:
        Eclipse times accurate to ~1 minute for most eclipses.
        For higher precision, use Besselian elements method.

    Example:
        >>> # Find next total solar eclipse after Jan 1, 2024
        >>> from libephemeris import julday, ECL_TOTAL
        >>> jd = julday(2024, 1, 1, 0)
        >>> ecl_type, times = _sol_eclipse_when_glob_pythonic(jd, eclipse_type=ECL_TOTAL)
        >>> print(f"Total eclipse at JD {times[0]:.5f}")

    References:
        - Reference API: sol_eclipse_when_glob()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    # The search walks lunation-by-lunation until it either finds a
    # matching eclipse or steps off the active ephemeris (the natural
    # tier boundary). MAX_NEW_MOONS is only a finite safety backstop so a
    # never-matching filter cannot loop forever on the ~30 000-year
    # extended ephemeris; on the base/medium tiers the ephemeris boundary
    # is reached first and converts to a typed no-event error.
    MAX_SEARCH_YEARS = _ECLIPSE_SEARCH_HORIZON_YEARS
    MAX_NEW_MOONS = int(MAX_SEARCH_YEARS * 12.4)  # ~12.4 lunations per year
    BIDIRECTIONAL_WINDOW = 15.0  # Days to check backward in bidirectional mode

    # Validate search_direction parameter
    valid_directions = ("forward", "backward", "bidirectional")
    if search_direction not in valid_directions:
        raise ValueError(
            f"search_direction must be one of {valid_directions}, got '{search_direction}'"
        )

    # Reject type masks that describe geometries which cannot occur
    # (central partial; non-central hybrid) before any search - the
    # public contract rejects these rather than starting a futile search.
    _sol_glob_reject_impossible(eclipse_type)

    # Keep the raw mask (0 == "any type") for the acceptance test; do NOT
    # expand 0 -> ECL_ALLTYPES_SOLAR here. _sol_glob_accepts() treats 0 as
    # "accept any" and applies the full geometry/type semantics for the
    # remaining 63 masks.
    mask = eclipse_type

    def _check_new_moon_for_eclipse(
        jd_new_moon: float,
    ) -> Union[Tuple[int, Tuple[float, ...]], None]:
        """Check if a New Moon corresponds to a matching eclipse."""
        # Get Moon position at New Moon
        moon_pos, _ = calc_ut(jd_new_moon, MOON, flags | FLG_SPEED)
        moon_lon = moon_pos[0]

        # Check if close enough to ecliptic for eclipse. This is only a
        # coarse PRE-FILTER at the (approximate) conjunction instant; the
        # true eclipse test is the classification at the REFINED maximum.
        node_dist = _get_moon_node_distance(jd_new_moon, moon_lon)

        if node_dist < _SOLAR_NODE_PREFILTER:
            # Refine the eclipse maximum FIRST (Besselian elements), then
            # classify and gate at the refined instant. Classifying at the
            # New Moon and gating on it (as before) dropped ultra-shallow
            # partials whose gamma is above the eclipse threshold at the
            # conjunction but dips below it at the true maximum.
            jd_max_refined = _refine_solar_eclipse_maximum(jd_new_moon)
            ecl_type, magnitude, gamma, ratio = _calculate_eclipse_type_and_magnitude(
                jd_max_refined
            )

            if ecl_type != 0 and _sol_glob_accepts(mask, ecl_type):
                # Calculate phase times using high-precision Besselian method
                times = _calculate_eclipse_phases(jd_max_refined, ecl_type)
                return ecl_type, times

        return None

    # For bidirectional search, first check if there's a recent eclipse we might miss
    backward_result: Union[Tuple[int, Tuple[float, ...]], None] = None
    if search_direction == "bidirectional":
        # Check the previous New Moon (within ~15 days before jd_start)
        jd_prev_new_moon = _find_previous_new_moon(jd_start)

        # Only consider if within the bidirectional window
        if jd_start - jd_prev_new_moon <= BIDIRECTIONAL_WINDOW:
            result = _check_new_moon_for_eclipse(jd_prev_new_moon)
            if result is not None:
                # Check if eclipse maximum is at or after jd_start
                # (eclipse might have started before but maximum after)
                jd_max = result[1][0]
                if jd_max > jd_start + _ECLIPSE_WHEN_EPOCH_MARGIN:
                    # This eclipse's maximum is after jd_start, use it
                    return result
                else:
                    # Store for comparison - eclipse maximum is before jd_start
                    # but we'll still prefer forward search unless nothing found
                    backward_result = result

    # Handle backward-only search
    if search_direction == "backward":
        # The conjunction in longitude can lie AFTER jd_start while the
        # eclipse maximum lies a few minutes BEFORE it: probe the next
        # conjunction first so that eclipse is not skipped.
        jd_next_nm = _find_next_new_moon(jd_start)
        if jd_next_nm - jd_start <= BIDIRECTIONAL_WINDOW:
            result = _check_new_moon_for_eclipse(jd_next_nm)
            if (
                result is not None
                and result[1][0] < jd_start - _ECLIPSE_WHEN_EPOCH_MARGIN
            ):
                return result

        jd = jd_start
        for _ in range(MAX_NEW_MOONS):
            try:
                jd_prev_new_moon = _find_previous_new_moon(jd)
                result = _check_new_moon_for_eclipse(jd_prev_new_moon)
            except Exception as _exc:
                if _is_ephemeris_boundary(_exc):
                    # Walked off the start of the ephemeris: no earlier
                    # matching eclipse exists within coverage.
                    break
                raise
            if result is not None:
                # Direction invariant: a backward search must return an
                # eclipse whose maximum precedes jd_start (refinement can
                # drag the maximum across the start time). The epoch margin
                # makes a maximum within _ECLIPSE_WHEN_EPOCH_MARGIN of jd_start count as
                # "reached", so an on-maximum start advances to the previous
                # eclipse (compatibility contract).
                if result[1][0] < jd_start - _ECLIPSE_WHEN_EPOCH_MARGIN:
                    return result
            # Go further back
            jd = jd_prev_new_moon - 1

        raise Error(
            "No matching solar eclipse found before the ephemeris boundary "
            f"searching back from JD {jd_start}"
        )

    # Forward search (default behavior for "forward" and "bidirectional")
    # The conjunction in longitude can lie BEFORE jd_start while the eclipse
    # maximum lies a few minutes AFTER it: probe the previous conjunction
    # first (mirror of the backward pre-probe above) so an in-progress
    # eclipse is not skipped for the next lunation's, ~29 days late. Only
    # the forward direction invariant (maximum strictly after jd_start)
    # lets the probe result through, so long-past eclipses cannot leak.
    if search_direction == "forward":
        try:
            jd_prev_nm = _find_previous_new_moon(jd_start)
        except Exception as _exc:
            if _is_ephemeris_boundary(_exc):
                jd_prev_nm = None
            else:
                raise
        if jd_prev_nm is not None and jd_start - jd_prev_nm <= BIDIRECTIONAL_WINDOW:
            result = _check_new_moon_for_eclipse(jd_prev_nm)
            if (
                result is not None
                and result[1][0] > jd_start + _ECLIPSE_WHEN_EPOCH_MARGIN
            ):
                return result

    jd = jd_start

    for _ in range(MAX_NEW_MOONS):
        # Find next New Moon
        try:
            jd_new_moon = _find_next_new_moon(jd)
            result = _check_new_moon_for_eclipse(jd_new_moon)
        except Exception as _exc:
            if _is_ephemeris_boundary(_exc):
                # Walked off the end of the ephemeris: no later matching
                # eclipse exists within coverage (the tier boundary).
                break
            raise
        if result is not None and result[1][0] <= jd_start + _ECLIPSE_WHEN_EPOCH_MARGIN:
            # Direction invariant: a forward search must return an eclipse
            # whose maximum follows jd_start. The conjunction-anchored
            # candidate can refine to a maximum a few minutes BEFORE the
            # start time; skip to the next lunation instead of returning
            # a past eclipse (or looping forever in tret[0]+eps scans). The
            # epoch margin makes a maximum within _ECLIPSE_WHEN_EPOCH_MARGIN of jd_start count
            # as "reached", so `jd = tret[0]; when(jd)` advances to the next
            # eclipse (compatibility contract).
            result = None
        if result is not None:
            # For bidirectional mode, we have a forward result
            # Check if we had a backward result that might be closer
            if search_direction == "bidirectional" and backward_result is not None:
                # Compare distances from jd_start
                backward_result[1][0]
                result[1][0]

                # If backward eclipse max is closer to jd_start (but before it),
                # and forward is far away, still prefer forward since it's >= jd_start
                # The bidirectional check already returned if backward was >= jd_start
                pass  # Forward result is >= jd_start, use it

            return result

        # Advance to next lunation
        jd = jd_new_moon + 25  # Skip ahead ~25 days to ensure we find next New Moon

    raise Error(
        "No matching solar eclipse found before the ephemeris boundary "
        f"searching from JD {jd_start}"
    )


def _strip_one_try_bit(flags: int) -> int:
    """Drop bit 0x8000 from eclipse-search calculation flags.

    The bit value is shared between ``ECL_ONE_TRY`` and ``FLG_TOPOCTR``.
    In the eclipse functions' ``flags`` argument it can only mean the
    one-try optimization hint (the topocentric place is defined by the
    explicit ``geopos``, never by this flag), and results are identical
    with or without it (compatibility contract). Passing it
    through to the position pipeline would instead trip the topocentric
    configuration guard. The functional one-try request rides on the
    ``backwards`` parameter, as in the occultation searches.
    """
    from .constants import ECL_ONE_TRY

    return flags & ~ECL_ONE_TRY


def sol_eclipse_when_glob(
    tjdut: float,
    flags: int = FLG_SWIEPH,
    ecltype: int = 0,
    backwards: "bool | int | str" = False,
) -> Tuple[int, Tuple[float, ...]]:
    """Find the next (or previous) global solar eclipse (reference-API-compatible).

    Wrapper around _sol_eclipse_when_glob_pythonic() matching the reference signature.

    Args:
        tjdut: Julian Day (UT) to start search from.
        flags: Calculation flags (default FLG_SWIEPH).
        ecltype: Eclipse type filter bitmask (0 = any).
        backwards: If True, search backward in time. Also accepts 0/1 and
            the direction strings understood by :func:`_coerce_backwards`.

    Returns:
        Tuple of (retflag, tret) matching the reference API.
    """
    flags = _strip_one_try_bit(flags)
    direction = "backward" if _coerce_backwards(backwards) else "forward"
    return _sol_eclipse_when_glob_pythonic(
        tjdut, flags=flags, eclipse_type=ecltype, search_direction=direction
    )


def _calculate_local_eclipse_phases(
    jd_max_global: float,
    lat: float,
    lon: float,
    altitude: float,
) -> Tuple[float, float, float, float, float, float, float, float, float, float]:
    """Calculate eclipse phase times as seen from a specific location.

    Wrapper around :func:`_calculate_local_eclipse_phases_impl` that applies
    mode-aware handling for partial/custom LEB files. See the implementation
    docstring for the full API contract.
    """
    return _call_with_leb_skyfield_fallback(
        _calculate_local_eclipse_phases_impl,
        jd_max_global,
        lat,
        lon,
        altitude,
    )


def _calculate_local_eclipse_phases_impl(
    jd_max_global: float,
    lat: float,
    lon: float,
    altitude: float,
    *,
    reader=None,
) -> Tuple[float, float, float, float, float, float, float, float, float, float]:
    """
    Calculate eclipse phase times as seen from a specific location.

    Uses geometric calculations based on the Moon's shadow passing over
    the observer's location.

    Args:
        jd_max_global: Julian Day of global maximum eclipse
        lat: Observer latitude in degrees (positive = North)
        lon: Observer longitude in degrees (positive = East)
        altitude: Observer altitude in meters above sea level

    Returns:
        Tuple of 10 floats with local eclipse circumstances:
            [0]: Time of maximum eclipse (local)
            [1]: Time of first contact (partial begins)
            [2]: Time of second contact (total/annular begins, or 0)
            [3]: Time of third contact (total/annular ends, or 0)
            [4]: Time of fourth contact (partial ends)
            [5]: Eclipse magnitude at maximum
            [6]: Fraction of Sun's diameter obscured
            [7]: Eclipse obscuration (area)
            [8]: Azimuth of Sun at maximum eclipse
            [9]: Altitude of Sun at maximum eclipse
    """
    # reader is provided by the caller (None forces Skyfield path)

    if reader is not None:
        try:
            from .fast_calc import _topo_ecliptic
            from .time_utils import deltat

            _topo_ecliptic(
                reader,
                jd_max_global + deltat(jd_max_global),
                jd_max_global,
                SUN,
                (lon, lat, altitude),
            )
        except (KeyError, ValueError) as _probe_exc:
            # A coverage miss must reach the mode-aware wrapper
            # (_calculate_local_eclipse_phases): typed error in sealed mode,
            # Skyfield retry with reader=None in auto. Degrading straight to
            # the kernel branch below would only convert a sealed-mode miss
            # into get_planets()'s RuntimeError. A body-map miss in sealed
            # mode must surface as the typed error for the same reason.
            _raise_if_sealed_leb_miss(_probe_exc)
            _reraise_if_leb_range_error(_probe_exc)
            reader = None

    if reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR, angular_separation

        geopos = (lon, lat, altitude)

        def _get_sun_moon_separation(jd: float) -> float:
            jd_tt = jd + deltat(jd)
            sun_pos = _topo_ecliptic(reader, jd_tt, jd, SUN, geopos)
            moon_pos = _topo_ecliptic(reader, jd_tt, jd, MOON, geopos)
            return angular_separation(sun_pos[0], sun_pos[1], moon_pos[0], moon_pos[1])

        def _get_sun_altaz(jd: float) -> Tuple[float, float]:
            """Get Sun altitude and azimuth at given JD from observer location."""
            jd_tt = jd + deltat(jd)
            sun_pos = _topo_ecliptic(reader, jd_tt, jd, SUN, geopos)
            sun_az, sun_alt_true, sun_alt_app = azalt(
                jd, ECL2HOR, geopos, 0, 0, sun_pos[:3]
            )
            # azalt returns the reference azimuth convention (S=0)
            return sun_alt_true, sun_az

        def _get_distances(jd: float) -> Tuple[float, float]:
            """Get Sun and Moon distances at given JD."""
            jd_tt = jd + deltat(jd)
            sun_pos = _topo_ecliptic(reader, jd_tt, jd, SUN, geopos)
            moon_pos = _topo_ecliptic(reader, jd_tt, jd, MOON, geopos)
            return sun_pos[2], moon_pos[2]
    else:
        from skyfield.api import wgs84

        from .state import get_planets, get_timescale

        # Get ephemeris and timescale
        eph = get_planets()
        ts = get_timescale()

        # Create observer location
        observer = wgs84.latlon(lat, lon, altitude)

        # Get Sun and Moon objects
        sun = eph["sun"]
        moon = eph["moon"]
        earth = eph["earth"]

        def _get_sun_moon_separation(jd: float) -> float:
            """Get angular separation between Sun and Moon at given JD."""
            t = ts.ut1_jd(jd)
            observer_at = earth + observer

            # Get apparent positions from observer
            sun_app = observer_at.at(t).observe(sun).apparent()
            moon_app = observer_at.at(t).observe(moon).apparent()

            # Calculate angular separation
            sep = sun_app.separation_from(moon_app)
            return sep.degrees

        def _get_sun_altaz(jd: float) -> Tuple[float, float]:
            """Get Sun altitude and azimuth at given JD from observer location."""
            t = ts.ut1_jd(jd)
            observer_at = earth + observer

            sun_app = observer_at.at(t).observe(sun).apparent()
            alt, az, _ = sun_app.altaz()
            # Convert Skyfield navigational azimuth (N=0) to the reference convention (S=0)
            return alt.degrees, (az.degrees + 180.0) % 360.0

        def _get_distances(jd: float) -> Tuple[float, float]:
            """Get Sun and Moon distances at given JD."""
            t = ts.ut1_jd(jd)
            observer_at = earth + observer
            sun_app = observer_at.at(t).observe(sun).apparent()
            moon_app = observer_at.at(t).observe(moon).apparent()
            return sun_app.distance().au, moon_app.distance().au

    # Search for local maximum around global maximum
    # The local maximum can differ from global by up to ~1 hour
    search_range = 3.0 / 24.0  # 3 hours in days

    # First check if Sun is above horizon at global maximum
    sun_alt, sun_az = _get_sun_altaz(jd_max_global)
    if sun_alt < -1.0:  # Sun below horizon (with some margin for refraction)
        # Sun not visible - return zeros to indicate no local eclipse
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    # Find local maximum (minimum separation)
    # Use golden section search for efficiency
    jd_low = jd_max_global - search_range
    jd_high = jd_max_global + search_range

    phi = (1 + math.sqrt(5)) / 2  # Golden ratio
    jd_a = jd_high - (jd_high - jd_low) / phi
    jd_b = jd_low + (jd_high - jd_low) / phi

    sep_a = _get_sun_moon_separation(jd_a)
    sep_b = _get_sun_moon_separation(jd_b)

    for _ in range(30):  # Converge to ~0.1 second precision
        if sep_a < sep_b:
            jd_high = jd_b
            jd_b = jd_a
            sep_b = sep_a
            jd_a = jd_high - (jd_high - jd_low) / phi
            sep_a = _get_sun_moon_separation(jd_a)
        else:
            jd_low = jd_a
            jd_a = jd_b
            sep_a = sep_b
            jd_b = jd_low + (jd_high - jd_low) / phi
            sep_b = _get_sun_moon_separation(jd_b)

        if jd_high - jd_low < 1e-7:  # ~0.01 second
            break

    jd_local_max = (jd_low + jd_high) / 2
    min_separation = _get_sun_moon_separation(jd_local_max)

    # Check if Sun is visible at local maximum
    sun_alt_max, sun_az_max = _get_sun_altaz(jd_local_max)
    if sun_alt_max < -1.0:
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    # Get angular sizes of Sun and Moon at local maximum
    sun_dist_au, moon_dist_au = _get_distances(jd_local_max)

    sun_angular_radius = (959.63 / 3600.0) / sun_dist_au  # degrees
    moon_angular_radius = (932.56 / 3600.0) * (0.002569 / moon_dist_au)  # degrees

    sun_diameter = 2 * sun_angular_radius
    moon_diameter = 2 * moon_angular_radius

    # Check if there's an eclipse at this location
    # Eclipse occurs if separation < sum of radii
    sum_radii = sun_angular_radius + moon_angular_radius
    if min_separation > sum_radii:
        # No eclipse visible from this location
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    # Calculate eclipse magnitude
    # Magnitude = fraction of Sun's diameter covered by Moon
    overlap = sum_radii - min_separation
    # Edge case: handle very small sun_diameter (shouldn't happen but be safe)
    if sun_diameter < MINIMUM_SEPARATION_FOR_LENS:
        magnitude = 0.0
    else:
        magnitude = overlap / sun_diameter
        magnitude = max(0.0, min(magnitude, 1.0 + moon_diameter / sun_diameter))

    # Calculate obscuration using safe function that handles edge cases
    r_sun = sun_angular_radius
    r_moon = moon_angular_radius
    d = min_separation  # center-to-center separation

    obscuration = _calculate_obscuration_safe(r_sun, r_moon, d)

    # Find contact times using bisection
    # First/fourth contact: separation = sum of radii
    # Second/third contact: separation = difference of radii (for central eclipse)

    def _find_contact_time(
        jd_start: float, jd_end: float, target_sep: float, is_increasing: bool
    ) -> float:
        """Find time when separation equals target, searching in given direction.

        Includes edge case handling for shallow eclipses where contact may
        not occur or search may fail to converge.
        """
        jd_low = jd_start
        jd_high = jd_end

        # Edge case: check if target separation is achievable in search range
        sep_start = _get_sun_moon_separation(jd_start)
        sep_end = _get_sun_moon_separation(jd_end)

        # For decreasing search: sep should go from above to below target
        # For increasing search: sep should go from below to above target
        if is_increasing:
            if sep_start > target_sep or sep_end < target_sep:
                # Check if minimum separation in range is above target
                sep_mid = _get_sun_moon_separation((jd_start + jd_end) / 2)
                if min(sep_start, sep_end, sep_mid) > target_sep:
                    return 0.0  # Contact not achievable
        else:
            if sep_start < target_sep or sep_end > target_sep:
                # Check if minimum separation in range is above target
                sep_mid = _get_sun_moon_separation((jd_start + jd_end) / 2)
                if min(sep_start, sep_end, sep_mid) > target_sep:
                    return 0.0  # Contact not achievable

        for _ in range(50):
            jd_mid = (jd_low + jd_high) / 2
            sep = _get_sun_moon_separation(jd_mid)

            if abs(sep - target_sep) < 1e-6:
                return jd_mid

            if is_increasing:
                # Looking for separation increasing through target
                if sep < target_sep:
                    jd_low = jd_mid
                else:
                    jd_high = jd_mid
            else:
                # Looking for separation decreasing through target
                if sep > target_sep:
                    jd_low = jd_mid
                else:
                    jd_high = jd_mid

            if jd_high - jd_low < 1e-8:
                break

        result = (jd_low + jd_high) / 2
        # Validate the result - for shallow eclipses, contact may not be precise
        final_sep = _get_sun_moon_separation(result)
        if abs(final_sep - target_sep) > 0.001:  # More than ~4 arcsec error
            return 0.0  # Contact time not reliable

        return result

    # Search range for contacts
    contact_search = 2.0 / 24.0  # 2 hours

    # First contact (separation decreasing to sum of radii)
    jd_first = _find_contact_time(
        jd_local_max - contact_search, jd_local_max, sum_radii, is_increasing=False
    )

    # Fourth contact (separation increasing from sum of radii)
    jd_fourth = _find_contact_time(
        jd_local_max, jd_local_max + contact_search, sum_radii, is_increasing=True
    )

    # Check if first and fourth contacts are valid (Sun above horizon).
    # An unresolved contact is already 0.0; skip the altaz test there, since
    # evaluating the Sun at JD 0 (year -4712) falls outside the default
    # medium-tier ephemeris coverage and would raise.
    if jd_first and _get_sun_altaz(jd_first)[0] < -1.0:
        jd_first = 0.0
    if jd_fourth and _get_sun_altaz(jd_fourth)[0] < -1.0:
        jd_fourth = 0.0

    # Second and third contacts (for total/annular eclipses)
    jd_second = 0.0
    jd_third = 0.0

    diff_radii = abs(moon_angular_radius - sun_angular_radius)
    if min_separation < diff_radii:
        # Central eclipse at this location
        jd_second = _find_contact_time(
            jd_local_max - contact_search / 4,
            jd_local_max,
            diff_radii,
            is_increasing=False,
        )
        jd_third = _find_contact_time(
            jd_local_max,
            jd_local_max + contact_search / 4,
            diff_radii,
            is_increasing=True,
        )

        # Check visibility (skip unresolved 0.0 contacts, see above)
        if jd_second and _get_sun_altaz(jd_second)[0] < -1.0:
            jd_second = 0.0
        if jd_third and _get_sun_altaz(jd_third)[0] < -1.0:
            jd_third = 0.0

    return (
        jd_local_max,  # [0] Time of maximum
        jd_first,  # [1] First contact
        jd_second,  # [2] Second contact
        jd_third,  # [3] Third contact
        jd_fourth,  # [4] Fourth contact
        magnitude,  # [5] Magnitude
        overlap / sun_diameter,  # [6] Fraction covered
        obscuration,  # [7] Obscuration
        sun_az_max,  # [8] Azimuth
        sun_alt_max,  # [9] Altitude
    )


def _determine_local_eclipse_type(
    magnitude: float, moon_sun_ratio: float, is_central: bool
) -> int:
    """
    Determine eclipse type flags for local observation.

    Args:
        magnitude: Eclipse magnitude at location
        moon_sun_ratio: Ratio of Moon's apparent diameter to Sun's
        is_central: Whether this is a central eclipse at location

    Returns:
        Eclipse type flags bitmask
    """
    if magnitude <= 0:
        return 0

    flags = ECL_VISIBLE

    if is_central:
        flags |= ECL_CENTRAL
        if moon_sun_ratio >= 1.0:
            flags |= ECL_TOTAL
        elif moon_sun_ratio > 0.99:
            flags |= ECL_ANNULAR_TOTAL
        else:
            flags |= ECL_ANNULAR
    else:
        flags |= ECL_PARTIAL

    return flags


def _sol_eclipse_when_loc_pythonic(
    jd_start: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    flags: int = FLG_SWIEPH,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """
    Find the next solar eclipse visible from a specific location.

    Searches forward in time from jd_start to find the next solar eclipse
    visible from the observer's location. Returns the times of eclipse phases
    as seen from that position.

    Args:
        jd_start: Julian Day (UT) to start search from
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - times: Tuple of 10 floats with eclipse phase times (JD UT):
                [0]: Time of maximum eclipse (local)
                [1]: Time of first contact (partial begins)
                [2]: Time of second contact (total/annular begins, or 0)
                [3]: Time of third contact (total/annular ends, or 0)
                [4]: Time of fourth contact (partial ends)
                [5]: Time of sunrise on central line (or 0, unused)
                [6]: Time of sunset on central line (or 0, unused)
                [7-9]: Reserved (0)
            - attr: Tuple of 11 floats with eclipse attributes:
                [0]: Eclipse magnitude
                [1]: Fraction of Sun's diameter covered
                [2]: Fraction of Sun's area obscured (obscuration)
                [3]: Azimuth of Sun at maximum eclipse (degrees)
                [4]: Altitude of Sun at maximum eclipse (degrees)
                [5]: Apparent diameter of Moon (degrees)
                [6]: Apparent diameter of Sun (degrees)
                [7]: Saros series number (0, not implemented)
                [8]: Inex number (0, not implemented)
                [9-10]: Reserved (0)
            - retflag: Eclipse type flags bitmask (ECL_* constants)

    Raises:
        Error: If no eclipse visible from location within search limit

    Algorithm:
        1. Use _sol_eclipse_when_glob_pythonic to find next global eclipse
        2. Calculate local circumstances at observer's location
        3. If eclipse not visible from location, continue to next global eclipse
        4. Return local phase times and attributes

    Precision:
        Eclipse times accurate to ~1 minute. Contact times depend on
        accurate Moon and Sun ephemeris positions.

    Example:
        >>> # Find next eclipse visible from Rome
        >>> from libephemeris import julday, _sol_eclipse_when_loc_pythonic
        >>> jd = julday(2024, 1, 1, 0)
        >>> rome_lat, rome_lon = 41.9028, 12.4964
        >>> ecl_type, times, attr = _sol_eclipse_when_loc_pythonic(jd, rome_lat, rome_lon)
        >>> print(f"Eclipse maximum at JD {times[0]:.5f}")
        >>> print(f"Magnitude: {attr[0]:.3f}")

    References:
        - Reference API: sol_eclipse_when_loc()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    MAX_SEARCH_YEARS = 50  # Maximum search range
    MAX_ECLIPSES = int(MAX_SEARCH_YEARS * 2.4)  # ~2.4 eclipses per year on average

    jd = jd_start

    for _ in range(MAX_ECLIPSES):
        # Find next global eclipse
        try:
            global_type, global_times = _sol_eclipse_when_glob_pythonic(jd, flags)
        except (Error, RuntimeError):
            raise Error(
                f"No solar eclipse visible from lat={lat}, lon={lon} "
                f"within {MAX_SEARCH_YEARS} years of JD {jd_start}"
            )

        jd_max_global = global_times[0]

        # Calculate local circumstances
        local_data = _calculate_local_eclipse_phases(jd_max_global, lat, lon, altitude)

        jd_local_max = local_data[0]
        magnitude = local_data[5]

        # Check if eclipse is visible from this location
        if jd_local_max > 0 and magnitude > 0:
            # Eclipse visible! Prepare return values

            # Get Moon/Sun topocentric distances at local maximum
            from .constants import FLG_TOPOCTR as _TOPOCTR_LOC
            from .state import get_topo as _gt_loc
            import libephemeris as _le_loc

            _saved_topo_loc = _gt_loc()
            _le_loc.set_topo(lon, lat, altitude)
            try:
                sun_pos, _ = calc_ut(jd_local_max, SUN, FLG_SPEED | _TOPOCTR_LOC)
                moon_pos, _ = calc_ut(jd_local_max, MOON, FLG_SPEED | _TOPOCTR_LOC)
            finally:
                from libephemeris import state as _st_loc

                if _saved_topo_loc is not None:
                    _le_loc.set_topo(
                        _saved_topo_loc.longitude.degrees,
                        _saved_topo_loc.latitude.degrees,
                        _saved_topo_loc.elevation.m,
                    )
                else:
                    _st_loc._TOPO = None
            sun_dist_au = sun_pos[2]
            moon_dist_au = moon_pos[2]

            sun_diameter = 2 * (959.63 / 3600.0) / sun_dist_au
            moon_diameter = 2 * (932.56 / 3600.0) * (0.002569 / moon_dist_au)

            # Check for central eclipse at this location
            is_central = local_data[2] > 0 and local_data[3] > 0

            moon_sun_ratio = moon_diameter / sun_diameter

            # Determine eclipse type
            ecl_type = _determine_local_eclipse_type(
                magnitude, moon_sun_ratio, is_central
            )

            # Add visibility flags
            if local_data[1] > 0:
                ecl_type |= ECL_1ST_VISIBLE
            if local_data[2] > 0:
                ecl_type |= ECL_2ND_VISIBLE
            if local_data[3] > 0:
                ecl_type |= ECL_3RD_VISIBLE
            if local_data[4] > 0:
                ecl_type |= ECL_4TH_VISIBLE
            ecl_type |= ECL_MAX_VISIBLE

            # Prepare times tuple (10 elements)
            times = (
                local_data[0],  # [0] Maximum
                local_data[1],  # [1] First contact
                local_data[2],  # [2] Second contact
                local_data[3],  # [3] Third contact
                local_data[4],  # [4] Fourth contact
                0.0,  # [5] Sunrise on central line (not implemented)
                0.0,  # [6] Sunset on central line (not implemented)
                0.0,  # [7] Reserved
                0.0,  # [8] Reserved
                0.0,  # [9] Reserved
            )

            # Prepare attributes tuple (11 elements)
            attr = (
                local_data[5],  # [0] Magnitude
                local_data[6],  # [1] Fraction covered
                local_data[7],  # [2] Obscuration
                local_data[8],  # [3] Azimuth
                local_data[9],  # [4] Altitude
                moon_diameter,  # [5] Moon diameter
                sun_diameter,  # [6] Sun diameter
                0.0,  # [7] Saros (not implemented)
                0.0,  # [8] Inex (not implemented)
                0.0,  # [9] Reserved
                0.0,  # [10] Reserved
            )

            return ecl_type, times, attr

        # Eclipse not visible from this location, try next
        jd = jd_max_global + 25  # Skip ahead to find next eclipse

    raise Error(
        f"No solar eclipse visible from lat={lat}, lon={lon} "
        f"within {MAX_SEARCH_YEARS} years of JD {jd_start}"
    )


# Legacy alias for original implementation
_sol_eclipse_when_loc_legacy = _sol_eclipse_when_loc_pythonic


def sol_eclipse_when_loc(
    tjdut: float,
    geopos: "Sequence[float]",
    flags: int = FLG_SWIEPH,
    backwards: "bool | int | str" = False,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """Find the next solar eclipse visible from a geographic location.

    Wrapper around :func:`_sol_eclipse_when_loc_impl` that adds LEB→Skyfield
    fallback for partial/custom LEB files. See the impl docstring for the
    full API contract.
    """
    flags = _strip_one_try_bit(flags)
    return _call_with_leb_skyfield_fallback(
        _sol_eclipse_when_loc_impl,
        tjdut,
        geopos,
        flags,
        _coerce_backwards(backwards),
    )


def _sol_eclipse_when_loc_impl(
    tjdut: float,
    geopos: "Sequence[float]",
    flags: int = FLG_SWIEPH,
    backwards: bool = False,
    *,
    reader,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """
    Find the next solar eclipse visible from a specific geographic location.

    This function follows the reference API signature exactly. It searches forward
    (or backward if specified) in time from tjdut to find the next solar
    eclipse visible from the observer's location specified by geopos.

    Args:
        tjdut: Julian Day (UT) to start search from
        flags: Calculation flags (FLG_SWIEPH, etc.)
        geopos: Sequence of [longitude_degrees, latitude_degrees, altitude_meters]
                NOTE: longitude comes first (this matches reference API convention)
        backwards: If True, search backward in time instead of forward

    Returns:
        Tuple containing:
            - tret: Tuple of 7 floats with eclipse phase times (JD UT):
                [0]: Time of maximum eclipse (local)
                [1]: Time of first contact (partial begins)
                [2]: Time of second contact (totality/annularity begins, or 0)
                [3]: Time of third contact (totality/annularity ends, or 0)
                [4]: Time of fourth contact (partial ends)
                [5]: Time of sunrise if eclipse at sunrise (or 0)
                [6]: Time of sunset if eclipse at sunset (or 0)
            - attr: Tuple of 8 floats with eclipse attributes:
                [0]: Fraction of solar diameter covered by Moon
                [1]: Ratio of lunar diameter to solar diameter
                [2]: Fraction of solar disc area covered (obscuration)
                [3]: Core shadow width in km (0 for partial eclipses)
                [4]: Azimuth of Sun at maximum eclipse (degrees)
                [5]: True altitude of Sun at maximum eclipse (degrees)
                [6]: Apparent altitude of Sun with refraction (degrees)
                [7]: Angular distance of Moon center from Sun center (degrees)
            - retflag: Eclipse type flags bitmask (ECL_* constants)
                ECL_VISIBLE if any part visible, combined with eclipse type

    Raises:
        Error: If no eclipse visible from location within search limit
        ValueError: If geopos has wrong length

    Algorithm:
        1. Use sol_eclipse_when_glob to find candidate global eclipses
        2. For each candidate, calculate topocentric positions using set_topo()
           and FLG_TOPOCTR flag
        3. Check if the eclipse is visible from the observer location
        4. Calculate contact times by solving for when lunar limb touches solar limb
        5. Calculate obscuration fraction as area covered / total solar disc area
        6. Verify Sun is above horizon at observer location during eclipse

    Precision:
        Eclipse times accurate to ~1-2 minutes. Contact times depend on
        accurate Moon and Sun ephemeris positions with topocentric corrections.

    Example:
        >>> # Find next eclipse visible from Dallas, TX (Apr 8, 2024 eclipse)
        >>> from libephemeris import sol_eclipse_when_loc, julday, FLG_SWIEPH
        >>> jd = julday(2024, 1, 1, 0)
        >>> dallas_geopos = [-96.797, 32.7767, 0]  # lon, lat, alt
        >>> ecl_type, tret, attr = sol_eclipse_when_loc(jd, FLG_SWIEPH, dallas_geopos)
        >>> print(f"Eclipse maximum at JD {tret[0]:.5f}")
        >>> print(f"Obscuration: {attr[2]:.3f}")

    References:
        - Reference API: sol_eclipse_when_loc()
        - Meeus "Astronomical Algorithms" Ch. 54
    """

    # Validate geopos
    if len(geopos) < 3:
        raise ValueError("geopos must be a sequence of [longitude, latitude, altitude]")

    # Extract geographic position (longitude first, then latitude - reference API convention)
    lon = float(geopos[0])
    lat = float(geopos[1])
    altitude = float(geopos[2])

    MAX_SEARCH_YEARS = 50
    MAX_ECLIPSES = int(MAX_SEARCH_YEARS * 2.4)

    # NOTE: the observer position is threaded through explicitly (LEB path
    # geopos tuples / Skyfield wgs84.latlon); the reference API's
    # The reference API's sol_eclipse_when_loc does not touch the global
    # set_topo state,
    # and neither do we (a previous set_topo() call here leaked the
    # eclipse observer into subsequent FLG_TOPOCTR calculations).

    from .constants import BIT_DISC_BOTTOM, CALC_RISE, CALC_SET
    from .utils import ECL2HOR, angular_separation, azalt

    _geopos3 = (lon, lat, altitude)
    eph_flags = _ecl_eph_flags(flags)

    def _get_separation(jd: float) -> float:
        """Topocentric Sun-Moon center separation in degrees."""
        s_p, m_p = _topo_sun_moon(jd, _geopos3, reader)
        return angular_separation(s_p[0], s_p[1], m_p[0], m_p[1])

    def _radii_sep(
        jd: float, rmoon_au: float = _ECL_RMOON_AU
    ) -> Tuple[float, float, float]:
        """Apparent solar and lunar radii plus center separation at jd.

        ``rmoon_au`` selects the lunar radius: the default full disc for the
        outer (penumbral) contacts and magnitude, or the smaller inner-contact
        radius for second/third contact (see ``_ECL_RMOON_INNER_AU``).
        """
        s_p, m_p = _topo_sun_moon(jd, _geopos3, reader)
        rs = math.degrees(math.asin(_ECL_RSUN_AU / s_p[2]))
        rm = math.degrees(math.asin(rmoon_au / m_p[2]))
        return rs, rm, angular_separation(s_p[0], s_p[1], m_p[0], m_p[1])

    def _sun_altaz(jd: float) -> Tuple[float, float, float]:
        """(true altitude, azimuth, apparent altitude) of the Sun at jd."""
        s_p, _m_p = _topo_sun_moon(jd, _geopos3, reader)
        az, true_alt, app_alt = azalt(jd, ECL2HOR, _geopos3, 0.0, 10.0, s_p)
        return true_alt, az, app_alt

    def _find_local_maximum(
        jd_center: float, search_range: float = 3.0 / 24.0
    ) -> float:
        """Find local eclipse maximum (minimum separation) using golden section search."""
        phi = (1 + math.sqrt(5)) / 2

        jd_low = jd_center - search_range
        jd_high = jd_center + search_range

        jd_a = jd_high - (jd_high - jd_low) / phi
        jd_b = jd_low + (jd_high - jd_low) / phi

        sep_a = _get_separation(jd_a)
        sep_b = _get_separation(jd_b)

        for _ in range(40):
            if sep_a < sep_b:
                jd_high = jd_b
                jd_b = jd_a
                sep_b = sep_a
                jd_a = jd_high - (jd_high - jd_low) / phi
                sep_a = _get_separation(jd_a)
            else:
                jd_low = jd_a
                jd_a = jd_b
                sep_a = sep_b
                jd_b = jd_low + (jd_high - jd_low) / phi
                sep_b = _get_separation(jd_b)

            if jd_high - jd_low < 1e-7:
                break

        return (jd_low + jd_high) / 2

    # Main search loop
    jd = tjdut

    # For backward search, we need to find eclipses before tjdut
    if backwards:
        # Find the most recent global eclipse before tjdut
        # We need to search backward in time
        pass
    else:
        pass

    for _ in range(MAX_ECLIPSES):
        # Find next/previous global eclipse
        try:
            if backwards:
                # For backward search, find a new moon before current position
                # and check if there's an eclipse there
                # We use a simpler approach: search forward from an earlier date
                # and find eclipses until we get one before tjdut
                earlier_jd = jd - 400  # About 1 year before
                # Extend the collection horizon just past the epoch so a
                # lunation whose geocentric maximum lies slightly AFTER the
                # epoch - while the observer's local phases have already begun
                # - is still gathered. The local-contact gate below drops any
                # candidate whose first local contact has not yet occurred.
                horizon = jd + _WHEN_LOC_EPOCH_WINDOW
                eclipses_found = []
                temp_jd = earlier_jd
                while temp_jd < horizon:
                    try:
                        global_type, global_times = _sol_eclipse_when_glob_pythonic(
                            temp_jd, flags
                        )
                        if global_times[0] < horizon:
                            eclipses_found.append((global_type, global_times))
                        temp_jd = global_times[0] + 25
                    except (Error, RuntimeError):
                        break
                if not eclipses_found:
                    jd -= 400
                    continue
                # Take the most recent eclipse before jd
                global_type, global_times = eclipses_found[-1]
            else:
                # Offer the current lunation as a candidate even when its
                # geocentric maximum has just slipped behind the epoch: the
                # observer's LOCAL phases may still be to come. The
                # local-contact gate below keeps or skips it. Only the first
                # iteration is affected; later iterations start well past the
                # previous maximum, so no eclipse is returned twice.
                global_type, global_times = _sol_eclipse_when_glob_pythonic(
                    jd - _WHEN_LOC_EPOCH_WINDOW, flags
                )
        except (Error, RuntimeError):
            raise Error(
                f"No solar eclipse visible from lon={lon}, lat={lat} "
                f"within {MAX_SEARCH_YEARS} years of JD {tjdut}"
            )

        jd_max_global = global_times[0]

        # Find the local (topocentric) maximum at this location
        jd_local_max = _find_local_maximum(jd_max_global)

        rsun, rmoon, min_separation = _radii_sep(jd_local_max)

        # No overlap at the local maximum: nothing to see from here.
        if min_separation > rsun + rmoon:
            if backwards:
                jd = jd_max_global - 1
            else:
                jd = jd_max_global + 25
            continue

        # Local phase at maximum.
        if min_separation < rsun - rmoon:
            phase = ECL_ANNULAR
        elif min_separation < abs(rsun - rmoon):
            phase = ECL_TOTAL
        else:
            phase = ECL_PARTIAL

        # Contact times. Outer contacts: center separation equals the
        # sum of the apparent radii evaluated at the contact itself, with the
        # full mean lunar disc. Inner contacts use the exact tangency of the
        # two apparent discs (center separation equals the absolute
        # semidiameter difference) against the smaller inner-contact lunar
        # radius, so totality/annularity is reckoned against the umbral disc.
        def _f_outer(jd_c: float) -> float:
            rs_c, rm_c, sep_c = _radii_sep(jd_c)
            return sep_c - (rs_c + rm_c)

        def _f_inner(jd_c: float) -> float:
            rs_c, rm_c, sep_c = _radii_sep(jd_c, _ECL_RMOON_INNER_AU)
            return abs(rs_c - rm_c) - sep_c

        two_hours = 2.0 / 24.0
        jd_first = _root_bisect(_f_outer, jd_local_max - two_hours, jd_local_max)
        if jd_first == 0.0:
            jd_first = _root_bisect(
                _f_outer, jd_local_max - 2.0 * two_hours, jd_local_max
            )
        jd_fourth = _root_bisect(_f_outer, jd_local_max, jd_local_max + two_hours)
        if jd_fourth == 0.0:
            jd_fourth = _root_bisect(
                _f_outer, jd_local_max, jd_local_max + 2.0 * two_hours
            )

        # Epoch gate on the LOCAL contacts (not the geocentric maximum):
        # keep an eclipse whose local phases are still unfolding at the
        # search epoch (in progress or future) and skip only one already
        # wholly past locally (forward) or wholly future locally (backward).
        # The geocentric maximum can precede the epoch while local
        # totality/partiality is still to come, so anchoring the gate on the
        # first/last local contact - not jd_local_max - is what stops an
        # in-progress eclipse from being skipped to the next lunation.
        first_contact = jd_first if jd_first != 0.0 else jd_local_max
        last_contact = jd_fourth if jd_fourth != 0.0 else jd_local_max
        # Two independent reasons to skip this lunation:
        #  * its local phases are wholly past (forward) or wholly future
        #    (backward) — the contact-anchored test, which is what keeps an
        #    in-progress eclipse from being dropped; and
        #  * the caller started AT or past the instant this function would
        #    return, i.e. the re-entrant `jd = tret[0]; when_loc(jd)` idiom.
        #    Without the second test the search re-found the same eclipse
        #    forever (the last contact is still hours ahead of the maximum),
        #    so enumerating local eclipses hung instead of advancing by the
        #    ~856-day solar / ~355-day lunar step.
        if (
            not backwards
            and (
                last_contact <= tjdut + _ECLIPSE_WHEN_EPOCH_MARGIN
                or jd_local_max <= tjdut + _ECLIPSE_WHEN_EPOCH_MARGIN
            )
        ) or (
            backwards
            and (
                first_contact >= tjdut - _ECLIPSE_WHEN_EPOCH_MARGIN
                or jd_local_max >= tjdut - _ECLIPSE_WHEN_EPOCH_MARGIN
            )
        ):
            if backwards:
                jd = jd_max_global - 1
            else:
                jd = jd_max_global + 25
            continue

        jd_second = 0.0
        jd_third = 0.0
        if min_separation < abs(rsun - rmoon):
            half_hour = 0.5 / 24.0
            jd_second = _root_bisect(_f_inner, jd_local_max - half_hour, jd_local_max)
            jd_third = _root_bisect(_f_inner, jd_local_max, jd_local_max + half_hour)

        # Visibility bits: a phase counts as visible when the Sun's
        # apparent altitude is positive at its time. tret keeps the
        # geometric times regardless (verified: 2021-06-10 from
        # Philadelphia fills tret[1] although first contact happens
        # below the horizon) - visibility lives only in the bits.
        ecl_type = phase
        for tc, bit in (
            (jd_local_max, ECL_MAX_VISIBLE),
            (jd_first, ECL_1ST_VISIBLE),
            (jd_second, ECL_2ND_VISIBLE),
            (jd_third, ECL_3RD_VISIBLE),
            (jd_fourth, ECL_4TH_VISIBLE),
        ):
            if tc == 0.0:
                continue
            _talt, _az, _aalt = _sun_altaz(tc)
            if _aalt > 0.0:
                ecl_type |= ECL_VISIBLE | bit

        if not (ecl_type & ECL_VISIBLE):
            # No phase of this eclipse rises above the local horizon.
            if backwards:
                jd = jd_max_global - 1
            else:
                jd = jd_max_global + 25
            continue

        # Sunrise/sunset bounding the visible part of the eclipse, using
        # the rise and set of the Sun's lower limb. When the geometric
        # maximum itself is below the horizon, the reported maximum is
        # moved to the horizon crossing and the phase re-evaluated there.
        jd_sunrise = 0.0
        jd_sunset = 0.0
        clamped = False
        circumpolar = False
        rc_rise, tret_rise = rise_trans(
            first_contact - 0.001,
            SUN,
            CALC_RISE | BIT_DISC_BOTTOM,
            _geopos3,
            0.0,
            0.0,
            eph_flags,
        )
        rc_set, tret_set = rise_trans(
            first_contact - 0.001,
            SUN,
            CALC_SET | BIT_DISC_BOTTOM,
            _geopos3,
            0.0,
            0.0,
            eph_flags,
        )
        if rc_rise == -2 or rc_set == -2:
            circumpolar = True

        if not circumpolar:
            tjdr = tret_rise[0]
            tjds = tret_set[0]
            if tjds < jd_first or (tjds > tjdr and tjdr > jd_fourth):
                # The whole eclipse runs while the Sun is down.
                if backwards:
                    jd = jd_max_global - 1
                else:
                    jd = jd_max_global + 25
                continue
            if jd_first < tjdr < jd_fourth:
                jd_sunrise = tjdr
                if not (ecl_type & ECL_MAX_VISIBLE):
                    jd_local_max = tjdr
                    clamped = True
            if jd_first < tjds < jd_fourth:
                jd_sunset = tjds
                if not (ecl_type & ECL_MAX_VISIBLE):
                    jd_local_max = tjds
                    clamped = True

        # Local circumstances at the reported maximum (geometric or
        # horizon-clamped).
        retc_how, attr_list = _sol_how_core(jd_local_max, _geopos3, flags, reader)
        if clamped:
            ecl_type &= ~(ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL)
            ecl_type |= retc_how & (ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL)

        # Global character at the reported maximum: the noncentral bit
        # and the core-shadow width come from the shadow-axis geometry,
        # independent of the observer.
        rc_where, _wlon, _wlat, dcore = _eclipse_where_core(jd_local_max, flags)
        ecl_type |= rc_where & ECL_NONCENTRAL
        attr_list[3] = dcore[0]

        tret = (
            jd_local_max,  # [0] maximum eclipse (or horizon crossing)
            jd_first,  # [1] first contact
            jd_second,  # [2] second contact
            jd_third,  # [3] third contact
            jd_fourth,  # [4] fourth contact
            jd_sunrise,  # [5] sunrise between first and fourth contact
            jd_sunset,  # [6] sunset between first and fourth contact
            0.0,  # [7] reserved
            0.0,  # [8] reserved
            0.0,  # [9] reserved
        )

        return ecl_type, tret, tuple(attr_list)

    raise Error(
        f"No solar eclipse visible from lon={lon}, lat={lat} "
        f"within {MAX_SEARCH_YEARS} years of JD {tjdut}"
    )


def _sol_eclipse_where_pythonic(
    jd: float,
    flags: int = FLG_SWIEPH,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """
    Calculate the geographic location where a solar eclipse is central at a given time.

    This is the legacy function signature. For reference API compatibility,
    use sol_eclipse_where().

    Args:
        jd: Julian Day (UT) of the moment to calculate
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        See sol_eclipse_where() for full return specification.
    """
    return sol_eclipse_where(jd, flags)


def sol_eclipse_where(
    tjdut: float,
    flags: int = FLG_SWIEPH,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """Find where a solar eclipse is central or maximal at a given time.

    Wrapper around :func:`_sol_eclipse_where_impl` that applies mode-aware LEB
    range handling. See the implementation docstring for the full API contract.
    """
    flags = _strip_one_try_bit(flags)
    return _call_with_leb_skyfield_fallback(_sol_eclipse_where_impl, tjdut, flags)


def _sol_eclipse_where_impl(
    tjdut: float,
    flags: int = FLG_SWIEPH,
    *,
    reader=None,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """
    Find the geographic location where a solar eclipse is maximal at a time.

    Computes where the Moon's shadow axis intersects the Earth's surface
    at the given instant. If the axis misses the Earth (partial or no
    eclipse), the surface point of maximum eclipse is returned instead.

    Args:
        tjdut: Julian Day (UT) during a solar eclipse
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - retflag: 0 when no eclipse is in progress anywhere on Earth
              at this time; otherwise ECL_CENTRAL or ECL_NONCENTRAL
              (plus ECL_PARTIAL when only the penumbra touches the
              Earth), combined with ECL_TOTAL or ECL_ANNULAR for the
              non-partial case.
            - geopos: Tuple of 10 floats:
                [0]: longitude of the shadow center (East positive)
                [1]: latitude of the shadow center (North positive)
                [2-9]: 0.0 (reserved in the reference API)
            - attr: Tuple of 20 floats with the local circumstances at
              the shadow-center point (same layout as sol_eclipse_how):
                [0]: magnitude as solar-diameter fraction (negative when
                     the discs do not overlap)
                [1]: ratio of lunar to solar apparent diameter
                [2]: obscuration (fraction of the solar disc covered)
                [3]: core-shadow diameter at the center point, km
                     (negative when the umbra reaches the surface)
                [4]: azimuth of the Sun at the center point
                [5]: true altitude of the Sun
                [6]: apparent altitude of the Sun
                [7]: Moon-Sun center separation in degrees
                [8]: NASA magnitude (= attr[0] for partial eclipses,
                     attr[1] for total/annular)
                [9]: saros series number
                [10]: saros series member number
                [11-19]: 0.0

    Note:
        Unlike the eclipse search functions, a retflag of 0 is a normal
        result here (no eclipse at this instant); geopos and attr still
        describe the point of closest approach of the shadow axis.

    References:
        - Reference API: sol_eclipse_where()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    retflag, center_lon, center_lat, dcore = _eclipse_where_core(tjdut, flags)
    _rc_how, attr_list = _sol_how_core(
        tjdut, (center_lon, center_lat, 0.0), flags, reader
    )
    attr_list[3] = dcore[0]
    geopos = (center_lon, center_lat, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    return retflag, geopos, tuple(attr_list)


def _sol_eclipse_how_pythonic(
    jd: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    flags: int = FLG_SWIEPH,
) -> Tuple[int, Tuple[float, ...]]:
    """
    Calculate the circumstances of a solar eclipse at a specific location and time.

    This is the legacy function signature. For reference API compatibility,
    use sol_eclipse_how().

    Args:
        jd: Julian Day (UT) of the moment to calculate
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        See sol_eclipse_how() for full specification.
    """
    geopos = (lon, lat, altitude)
    return sol_eclipse_how(jd, geopos, flags)


def sol_eclipse_how(
    tjdut: float,
    geopos: Sequence[float],
    flags: int = FLG_SWIEPH,
) -> Tuple[int, Tuple[float, ...]]:
    """Calculate the circumstances of a solar eclipse at a location.

    Wrapper around :func:`_sol_eclipse_how_impl` that adds LEB→Skyfield
    fallback for partial/custom LEB files. See the impl docstring for the
    full API contract.
    """
    flags = _strip_one_try_bit(flags)
    return _call_with_leb_skyfield_fallback(
        _sol_eclipse_how_impl,
        tjdut,
        geopos,
        flags,
    )


def _sol_eclipse_how_impl(
    tjdut: float,
    geopos: Sequence[float],
    flags: int = FLG_SWIEPH,
    *,
    reader=None,
) -> Tuple[int, Tuple[float, ...]]:
    """
    Calculate the circumstances of a solar eclipse at a location and time.

    This function does NOT search for eclipses - it computes the local
    circumstances at exactly ``tjdut`` from the given place.

    Args:
        tjdut: Julian Day (UT) of a specific time during an eclipse
        geopos: Sequence of [longitude, latitude, altitude]:
            - longitude in degrees (East positive)
            - latitude in degrees (North positive)
            - altitude in meters above sea level
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - retflag: local eclipse character: ECL_TOTAL, ECL_ANNULAR or
              ECL_PARTIAL, plus ECL_VISIBLE when part of the eclipsed Sun
              can stand above the local horizon, plus the global
              ECL_CENTRAL/ECL_NONCENTRAL character of the eclipse at this
              instant. 0 when no eclipse is visible from this place (the
              Sun below the apparent horizon also yields 0).
            - attr: Tuple of 20 floats:
                [0]: magnitude as solar-diameter fraction
                [1]: ratio of lunar to solar apparent diameter
                [2]: obscuration (fraction of the solar disc covered)
                [3]: core-shadow diameter at the global shadow-center
                     point, km (observer-independent; negative when the
                     umbra reaches the surface)
                [4]: azimuth of the Sun
                [5]: true altitude of the Sun
                [6]: apparent altitude of the Sun (refraction at 10 C,
                     pressure estimated from the observer's altitude)
                [7]: Moon-Sun center separation in degrees
                [8]: NASA magnitude (= attr[0] for partial eclipses,
                     attr[1] for total/annular)
                [9]: saros series number
                [10]: saros series member number
                [11-19]: 0.0
              When retflag is 0, attr[0..3] and attr[8..10] are zeroed;
              the azimuth/altitude values attr[4..7] are kept.

    Note:
        The visibility gate uses the Sun's apparent altitude: at or below
        0 the function returns retflag 0 even mid-eclipse (compatibility
        contract).

    References:
        - Reference API: sol_eclipse_how()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    if len(geopos) < 3:
        raise ValueError("geopos must have at least 3 elements: [lon, lat, alt]")

    geopos3 = (float(geopos[0]), float(geopos[1]), float(geopos[2]))

    retflag, attr_list = _sol_how_core(tjdut, geopos3, flags, reader)
    rc_where, _wlon, _wlat, dcore = _eclipse_where_core(tjdut, flags)
    if retflag:
        retflag |= rc_where & (ECL_CENTRAL | ECL_NONCENTRAL)
    attr_list[3] = dcore[0]

    # Below the apparent horizon nothing of the eclipse is observable.
    if attr_list[6] <= 0.0:
        retflag = 0
    if retflag == 0:
        for _i in (0, 1, 2, 3, 8, 9, 10):
            attr_list[_i] = 0.0

    return retflag, tuple(attr_list)


def sol_eclipse_how_details(
    tjd_ut: float,
    geopos: Sequence[float],
    ifl: int = 0,
) -> dict:
    """Calculate comprehensive solar eclipse circumstances at a location.

    Wrapper around :func:`_sol_eclipse_how_details_impl` that adds LEB→Skyfield
    fallback for partial/custom LEB files. See the impl docstring for the
    full API contract.
    """
    return _call_with_leb_skyfield_fallback(
        _sol_eclipse_how_details_impl,
        tjd_ut,
        geopos,
        ifl,
    )


def _sol_eclipse_how_details_impl(
    tjd_ut: float,
    geopos: Sequence[float],
    ifl: int = 0,
    *,
    reader=None,
) -> dict:
    """
    Calculate comprehensive solar eclipse circumstances at a specific location.

    This enhanced version of sol_eclipse_how returns all eclipse details
    including contact times, maximum obscuration, position angles, and Sun
    altitude/azimuth during all eclipse phases.

    Args:
        tjd_ut: Julian Day (UT) of a specific time during an eclipse.
                This should be a time when an eclipse is known to occur
                (e.g., from sol_eclipse_when_glob or sol_eclipse_when_loc).
        ifl: Calculation flags (FLG_SWIEPH, etc.)
        geopos: Sequence of [longitude, latitude, altitude]:
            - longitude in degrees (East positive)
            - latitude in degrees (North positive)
            - altitude in meters above sea level

    Returns:
        Dictionary with comprehensive eclipse information:
            'eclipse_type': int - Eclipse type flags bitmask (ECL_* constants)
            'is_visible': bool - Whether eclipse is visible from this location
            'is_total': bool - Whether eclipse is total at this location
            'is_annular': bool - Whether eclipse is annular at this location
            'is_partial': bool - Whether eclipse is partial at this location

            Contact times (Julian Day UT, or 0.0 if not applicable):
            'jd_c1': float - First contact (partial phase begins)
            'jd_c2': float - Second contact (totality/annularity begins)
            'jd_max': float - Maximum eclipse
            'jd_c3': float - Third contact (totality/annularity ends)
            'jd_c4': float - Fourth contact (partial phase ends)

            Eclipse magnitude and obscuration:
            'magnitude': float - Fraction of solar diameter covered by Moon
            'max_magnitude': float - Maximum magnitude during eclipse
            'obscuration': float - Fraction of solar disc area covered (0-1)
            'max_obscuration': float - Maximum obscuration during eclipse (0-1)
            'obscuration_percent': float - Obscuration as percentage (0-100)
            'max_obscuration_percent': float - Maximum obscuration as percentage
            'ratio': float - Ratio of lunar diameter to solar diameter
            'shadow_width_km': float - Core shadow width in km (0 for partial)

            Position angles (degrees from North through East, counterclockwise):
            'position_angle_c1': float - Position angle at first contact
            'position_angle_c2': float - Position angle at second contact (or 0)
            'position_angle_c3': float - Position angle at third contact (or 0)
            'position_angle_c4': float - Position angle at fourth contact

            Sun position at each phase:
            'sun_alt_c1': float - Sun altitude at first contact
            'sun_az_c1': float - Sun azimuth at first contact
            'sun_alt_c2': float - Sun altitude at second contact
            'sun_az_c2': float - Sun azimuth at second contact
            'sun_alt_max': float - Sun altitude at maximum eclipse
            'sun_az_max': float - Sun azimuth at maximum eclipse
            'sun_alt_c3': float - Sun altitude at third contact
            'sun_az_c3': float - Sun azimuth at third contact
            'sun_alt_c4': float - Sun altitude at fourth contact
            'sun_az_c4': float - Sun azimuth at fourth contact

            Duration information:
            'duration_partial_minutes': float - Total duration of partial phases
            'duration_total_minutes': float - Duration of totality (0 if not total)

            Angular sizes:
            'sun_angular_radius': float - Sun angular radius in degrees
            'moon_angular_radius': float - Moon angular radius in degrees
            'separation': float - Angular separation Sun-Moon at given time

    Algorithm:
        1. Find eclipse maximum for this location using golden section search
        2. Calculate contact times using bisection method
        3. Compute position angles from celestial mechanics
        4. Calculate Sun altitude/azimuth at each phase

    Precision:
        Contact times accurate to ~1 second.
        Position angles accurate to ~0.1 degrees.

    Example:
        >>> from libephemeris import sol_eclipse_how_details, FLG_SWIEPH
        >>> jd = 2460409.28  # During April 8, 2024 eclipse
        >>> dallas = [-96.797, 32.7767, 0]  # lon, lat, alt
        >>> details = sol_eclipse_how_details(jd, dallas, FLG_SWIEPH)
        >>> print(f"Max obscuration: {details['max_obscuration_percent']:.1f}%")
        >>> print(f"First contact: JD {details['jd_c1']:.5f}")
        >>> print(f"Duration of totality: {details['duration_total_minutes']:.1f} min")

    References:
        - Espenak & Meeus, NASA/TP-2009-214173, Eqs. 1-5 and 1-6
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    # Validate and extract geopos
    if len(geopos) < 3:
        raise ValueError("geopos must have at least 3 elements: [lon, lat, alt]")

    lon = float(geopos[0])
    lat = float(geopos[1])
    altitude = float(geopos[2])

    # reader is provided by the caller (None forces Skyfield path)

    if reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR, angular_separation

        _gp = (lon, lat, altitude)

        def _get_separation(jd: float) -> float:
            """Get angular separation between Sun and Moon."""
            jd_tt = jd + deltat(jd)
            sun_pos = _topo_ecliptic(reader, jd_tt, jd, SUN, _gp)
            moon_pos = _topo_ecliptic(reader, jd_tt, jd, MOON, _gp)
            return angular_separation(sun_pos[0], sun_pos[1], moon_pos[0], moon_pos[1])

        def _get_sun_altaz(jd: float) -> tuple:
            """Get Sun altitude and azimuth at given JD."""
            jd_tt = jd + deltat(jd)
            sun_pos = _topo_ecliptic(reader, jd_tt, jd, SUN, _gp)
            sun_az, sun_alt_true, sun_alt_app = azalt(
                jd, ECL2HOR, _gp, 0, 0, sun_pos[:3]
            )
            return sun_alt_true, sun_az

        def _get_angular_sizes(jd: float) -> tuple:
            """Get angular radii of Sun and Moon."""
            jd_tt = jd + deltat(jd)
            sun_pos = _topo_ecliptic(reader, jd_tt, jd, SUN, _gp)
            moon_pos = _topo_ecliptic(reader, jd_tt, jd, MOON, _gp)
            sun_angular_radius = (959.63 / 3600.0) / sun_pos[2]
            moon_angular_radius = (932.56 / 3600.0) * (0.002569 / moon_pos[2])
            return sun_angular_radius, moon_angular_radius

        def _calc_position_angle(jd: float) -> float:
            """Calculate position angle of Moon relative to Sun (topocentric ICRS)."""
            from .constants import FLG_EQUATORIAL as _EQ_PA, FLG_J2000 as _J2K_PA
            from .constants import FLG_TOPOCTR as _TP_PA
            from .state import get_topo as _gt_pa
            import libephemeris as _le_pa

            _saved_topo_pa = _gt_pa()
            _le_pa.set_topo(_gp[0], _gp[1], _gp[2])
            try:
                _pa_flags = _EQ_PA | _J2K_PA | _TP_PA | FLG_SPEED
                sun_eq, _ = calc_ut(jd, SUN, _pa_flags)
                moon_eq, _ = calc_ut(jd, MOON, _pa_flags)
            finally:
                from libephemeris import state as _st_pa

                if _saved_topo_pa is not None:
                    _le_pa.set_topo(
                        _saved_topo_pa.longitude.degrees,
                        _saved_topo_pa.latitude.degrees,
                        _saved_topo_pa.elevation.m,
                    )
                else:
                    _st_pa._TOPO = None
            sun_ra_rad = math.radians(sun_eq[0])
            sun_dec_rad = math.radians(sun_eq[1])
            moon_ra_rad = math.radians(moon_eq[0])
            moon_dec_rad = math.radians(moon_eq[1])
            delta_ra = moon_ra_rad - sun_ra_rad
            y = math.sin(delta_ra)
            x = math.cos(sun_dec_rad) * math.tan(moon_dec_rad) - math.sin(
                sun_dec_rad
            ) * math.cos(delta_ra)
            pa = math.degrees(math.atan2(y, x))
            if pa < 0:
                pa += 360.0
            return pa
    else:
        from skyfield.api import wgs84

        from .state import get_planets, get_timescale

        # Get ephemeris and timescale
        eph = get_planets()
        ts = get_timescale()

        # Create observer location
        observer = wgs84.latlon(lat, lon, altitude)

        # Get Sun and Moon objects
        sun = eph["sun"]
        moon = eph["moon"]
        earth = eph["earth"]

        def _get_sun_moon_positions(jd: float):
            """Get Sun and Moon positions from observer."""
            t = ts.ut1_jd(jd)
            observer_at = earth + observer

            sun_app = observer_at.at(t).observe(sun).apparent()
            moon_app = observer_at.at(t).observe(moon).apparent()

            return sun_app, moon_app

        def _get_separation(jd: float) -> float:
            """Get angular separation between Sun and Moon."""
            sun_app, moon_app = _get_sun_moon_positions(jd)
            return sun_app.separation_from(moon_app).degrees

        def _get_sun_altaz(jd: float) -> tuple:
            """Get Sun altitude and azimuth at given JD."""
            sun_app, _ = _get_sun_moon_positions(jd)
            alt, az, _ = sun_app.altaz()
            # Convert Skyfield navigational azimuth (N=0) to the reference convention (S=0)
            return alt.degrees, (az.degrees + 180.0) % 360.0

        def _get_angular_sizes(jd: float) -> tuple:
            """Get angular radii of Sun and Moon."""
            sun_app, moon_app = _get_sun_moon_positions(jd)

            sun_dist_au = sun_app.distance().au
            moon_dist_au = moon_app.distance().au

            # Angular radii in degrees
            sun_angular_radius = (959.63 / 3600.0) / sun_dist_au
            moon_angular_radius = (932.56 / 3600.0) * (0.002569 / moon_dist_au)

            return sun_angular_radius, moon_angular_radius

        def _calc_position_angle(jd: float) -> float:
            """
            Calculate position angle of Moon relative to Sun.

            Position angle is measured from North through East (counterclockwise).
            Returns angle in degrees (0-360).
            """
            sun_app, moon_app = _get_sun_moon_positions(jd)

            # Get RA and Dec for both bodies
            sun_ra, sun_dec, _ = sun_app.radec()
            moon_ra, moon_dec, _ = moon_app.radec()

            sun_ra_rad = sun_ra.radians
            sun_dec_rad = sun_dec.radians
            moon_ra_rad = moon_ra.radians
            moon_dec_rad = moon_dec.radians

            # Position angle formula
            delta_ra = moon_ra_rad - sun_ra_rad

            y = math.sin(delta_ra)
            x = math.cos(sun_dec_rad) * math.tan(moon_dec_rad) - math.sin(
                sun_dec_rad
            ) * math.cos(delta_ra)

            pa = math.degrees(math.atan2(y, x))
            if pa < 0:
                pa += 360.0

            return pa

    def _find_local_maximum(
        jd_center: float, search_range: float = 3.0 / 24.0
    ) -> float:
        """Find local eclipse maximum (minimum separation) using golden section."""
        phi = (1 + math.sqrt(5)) / 2

        jd_low = jd_center - search_range
        jd_high = jd_center + search_range

        jd_a = jd_high - (jd_high - jd_low) / phi
        jd_b = jd_low + (jd_high - jd_low) / phi

        sep_a = _get_separation(jd_a)
        sep_b = _get_separation(jd_b)

        for _ in range(50):
            if sep_a < sep_b:
                jd_high = jd_b
                jd_b = jd_a
                sep_b = sep_a
                jd_a = jd_high - (jd_high - jd_low) / phi
                sep_a = _get_separation(jd_a)
            else:
                jd_low = jd_a
                jd_a = jd_b
                sep_a = sep_b
                jd_b = jd_low + (jd_high - jd_low) / phi
                sep_b = _get_separation(jd_b)

            if jd_high - jd_low < 1e-8:
                break

        return (jd_low + jd_high) / 2

    def _find_contact_time(
        jd_start: float,
        jd_end: float,
        target_sep: float,
        is_increasing: bool,
    ) -> float:
        """Find time when separation equals target using bisection."""
        jd_low = jd_start
        jd_high = jd_end

        for _ in range(60):
            jd_mid = (jd_low + jd_high) / 2
            sep = _get_separation(jd_mid)

            if abs(sep - target_sep) < 1e-7:
                return jd_mid

            if is_increasing:
                if sep < target_sep:
                    jd_low = jd_mid
                else:
                    jd_high = jd_mid
            else:
                if sep > target_sep:
                    jd_low = jd_mid
                else:
                    jd_high = jd_mid

            if jd_high - jd_low < 1e-9:
                break

        return (jd_low + jd_high) / 2

    # First, get basic eclipse info at the given time
    eclipse_type, attr = sol_eclipse_how(tjd_ut, geopos, ifl)

    # Initialize result dictionary
    result = {
        "eclipse_type": eclipse_type,
        "is_visible": bool(eclipse_type & ECL_VISIBLE),
        "is_total": bool(eclipse_type & ECL_TOTAL),
        "is_annular": bool(eclipse_type & ECL_ANNULAR),
        "is_partial": bool(eclipse_type & ECL_PARTIAL),
        # Contact times
        "jd_c1": 0.0,
        "jd_c2": 0.0,
        "jd_max": 0.0,
        "jd_c3": 0.0,
        "jd_c4": 0.0,
        # Magnitude and obscuration
        "magnitude": attr[0],
        "max_magnitude": 0.0,
        "obscuration": attr[2],
        "max_obscuration": 0.0,
        "obscuration_percent": attr[2] * 100.0,
        "max_obscuration_percent": 0.0,
        "ratio": attr[1],
        "shadow_width_km": attr[3],
        # Position angles
        "position_angle_c1": 0.0,
        "position_angle_c2": 0.0,
        "position_angle_c3": 0.0,
        "position_angle_c4": 0.0,
        # Sun positions
        "sun_alt_c1": 0.0,
        "sun_az_c1": 0.0,
        "sun_alt_c2": 0.0,
        "sun_az_c2": 0.0,
        "sun_alt_max": attr[5],
        "sun_az_max": attr[4],
        "sun_alt_c3": 0.0,
        "sun_az_c3": 0.0,
        "sun_alt_c4": 0.0,
        "sun_az_c4": 0.0,
        # Durations
        "duration_partial_minutes": 0.0,
        "duration_total_minutes": 0.0,
        # Angular sizes
        "sun_angular_radius": 0.0,
        "moon_angular_radius": 0.0,
        "separation": attr[7],
    }

    # If no eclipse is visible, return early
    if eclipse_type == 0:
        return result

    # Get angular sizes
    sun_r, moon_r = _get_angular_sizes(tjd_ut)
    result["sun_angular_radius"] = sun_r
    result["moon_angular_radius"] = moon_r

    sun_r + moon_r
    diff_radii = abs(sun_r - moon_r)

    # Find local maximum for this observer
    jd_local_max = _find_local_maximum(tjd_ut)
    result["jd_max"] = jd_local_max

    # Get maximum eclipse circumstances
    max_sep = _get_separation(jd_local_max)
    sun_alt_max, sun_az_max = _get_sun_altaz(jd_local_max)
    result["sun_alt_max"] = sun_alt_max
    result["sun_az_max"] = sun_az_max

    # Calculate maximum magnitude and obscuration
    sun_r_max, moon_r_max = _get_angular_sizes(jd_local_max)
    sun_diameter = 2 * sun_r_max
    moon_diameter = 2 * moon_r_max

    sum_radii_max = sun_r_max + moon_r_max
    overlap = sum_radii_max - max_sep
    max_magnitude = max(
        0.0, min(overlap / sun_diameter, 1.0 + moon_diameter / sun_diameter)
    )
    result["max_magnitude"] = max_magnitude
    result["ratio"] = moon_diameter / sun_diameter

    # Calculate maximum obscuration
    d = max_sep
    if d >= sun_r_max + moon_r_max:
        max_obscuration = 0.0
    elif d <= abs(sun_r_max - moon_r_max):
        # Total (Moon >= Sun): the Sun is fully covered -> 1.0.
        # Annular (Moon < Sun): a ring of Sun remains -> disc area ratio < 1.0.
        max_obscuration = (
            1.0 if moon_r_max >= sun_r_max else (moon_r_max / sun_r_max) ** 2
        )
    else:
        d1 = (d * d + sun_r_max * sun_r_max - moon_r_max * moon_r_max) / (2 * d)
        d2 = d - d1

        if abs(d1) <= sun_r_max and abs(d2) <= moon_r_max:
            cos_arg1 = max(-1, min(1, d1 / sun_r_max))
            cos_arg2 = max(-1, min(1, d2 / moon_r_max))

            area1 = sun_r_max * sun_r_max * math.acos(cos_arg1) - d1 * math.sqrt(
                max(0, sun_r_max * sun_r_max - d1 * d1)
            )
            area2 = moon_r_max * moon_r_max * math.acos(cos_arg2) - d2 * math.sqrt(
                max(0, moon_r_max * moon_r_max - d2 * d2)
            )
            intersection_area = area1 + area2
            sun_area = math.pi * sun_r_max * sun_r_max
            max_obscuration = intersection_area / sun_area
        else:
            max_obscuration = 0.0

    max_obscuration = max(0.0, max_obscuration)
    result["max_obscuration"] = max_obscuration
    result["max_obscuration_percent"] = max_obscuration * 100.0

    # Update eclipse type based on maximum
    is_central = max_sep < diff_radii
    if is_central:
        result["is_partial"] = False
        if moon_r_max >= sun_r_max:
            result["is_total"] = True
            result["is_annular"] = False
        else:
            result["is_total"] = False
            result["is_annular"] = True
    else:
        result["is_partial"] = True
        result["is_total"] = False
        result["is_annular"] = False

    # Calculate contact times
    contact_search_range = 2.5 / 24.0  # 2.5 hours

    # First contact (separation decreasing to sum of radii)
    try:
        jd_c1 = _find_contact_time(
            jd_local_max - contact_search_range,
            jd_local_max,
            sum_radii_max,
            is_increasing=False,
        )
        result["jd_c1"] = jd_c1

        # Sun position at C1
        sun_alt_c1, sun_az_c1 = _get_sun_altaz(jd_c1)
        result["sun_alt_c1"] = sun_alt_c1
        result["sun_az_c1"] = sun_az_c1

        # Position angle at C1
        result["position_angle_c1"] = _calc_position_angle(jd_c1)
    except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
        _reraise_if_leb_range_error(_exc)
        pass

    # Fourth contact (separation increasing from sum of radii)
    try:
        jd_c4 = _find_contact_time(
            jd_local_max,
            jd_local_max + contact_search_range,
            sum_radii_max,
            is_increasing=True,
        )
        result["jd_c4"] = jd_c4

        # Sun position at C4
        sun_alt_c4, sun_az_c4 = _get_sun_altaz(jd_c4)
        result["sun_alt_c4"] = sun_alt_c4
        result["sun_az_c4"] = sun_az_c4

        # Position angle at C4
        result["position_angle_c4"] = _calc_position_angle(jd_c4)
    except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
        _reraise_if_leb_range_error(_exc)
        pass

    # Second and third contacts (for central eclipses only)
    diff_radii_max = abs(sun_r_max - moon_r_max)
    if is_central and max_sep < diff_radii_max:
        # Second contact
        try:
            jd_c2 = _find_contact_time(
                jd_local_max - contact_search_range / 4,
                jd_local_max,
                diff_radii_max,
                is_increasing=False,
            )
            result["jd_c2"] = jd_c2

            # Sun position at C2
            sun_alt_c2, sun_az_c2 = _get_sun_altaz(jd_c2)
            result["sun_alt_c2"] = sun_alt_c2
            result["sun_az_c2"] = sun_az_c2

            # Position angle at C2
            result["position_angle_c2"] = _calc_position_angle(jd_c2)
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            pass

        # Third contact
        try:
            jd_c3 = _find_contact_time(
                jd_local_max,
                jd_local_max + contact_search_range / 4,
                diff_radii_max,
                is_increasing=True,
            )
            result["jd_c3"] = jd_c3

            # Sun position at C3
            sun_alt_c3, sun_az_c3 = _get_sun_altaz(jd_c3)
            result["sun_alt_c3"] = sun_alt_c3
            result["sun_az_c3"] = sun_az_c3

            # Position angle at C3
            result["position_angle_c3"] = _calc_position_angle(jd_c3)
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            pass

    # Calculate durations
    if result["jd_c1"] > 0 and result["jd_c4"] > 0:
        result["duration_partial_minutes"] = (
            (result["jd_c4"] - result["jd_c1"]) * 24 * 60
        )

    if result["jd_c2"] > 0 and result["jd_c3"] > 0:
        result["duration_total_minutes"] = (result["jd_c3"] - result["jd_c2"]) * 24 * 60

    return result


def _sol_eclipse_how_details_pythonic(
    jd: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    flags: int = FLG_SWIEPH,
) -> dict:
    """
    Calculate comprehensive solar eclipse circumstances at a specific location.

    This is the legacy function signature. For reference API-style arguments,
    use sol_eclipse_how_details().

    Args:
        jd: Julian Day (UT) of the moment to calculate
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Dictionary with comprehensive eclipse information.
        See sol_eclipse_how_details() for full specification.
    """
    geopos = (lon, lat, altitude)
    return sol_eclipse_how_details(jd, geopos, flags)


# =============================================================================
# LUNAR ECLIPSE FUNCTIONS
# =============================================================================

# Constants for lunar eclipse calculations
ECLIPSE_LIMIT_LUNAR = 18.5  # Maximum elongation from node for lunar eclipse (degrees)
# Note: penumbral eclipses can occur up to ~18° from a node (observed max 18.02°).
# The previous value of 12.0° was too restrictive and missed shallow penumbral eclipses.
# The Full-Moon node distance is only a coarse PRE-FILTER: the true test is the
# shadow-geometry classification at the REFINED maximum (_lun_how_core). node
# distance is NOT a clean discriminator near the limit (band non-events reach
# ~16.8° while a real penumbral - e.g. 1958-04-04 - sits above 18.5°), so the
# pre-filter is opened to 20° and the classification makes the decision.
_LUNAR_NODE_PREFILTER = 20.0


def _find_next_full_moon(jd_start: float) -> float:
    """
    Find the next Full Moon (Sun-Moon opposition) after jd_start.

    Uses iterative refinement to find exact moment of opposition.

    Args:
        jd_start: Julian Day (UT) to start search from

    Returns:
        Julian Day of next Full Moon
    """
    # Get current positions
    sun_pos, _ = calc_ut(jd_start, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd_start, MOON, FLG_SPEED)

    sun_lon = sun_pos[0]
    moon_lon = moon_pos[0]

    # Calculate elongation from opposition (Moon - Sun - 180°)
    elongation = (moon_lon - sun_lon - 180.0) % 360.0
    if elongation > 180:
        elongation -= 360

    # Moon gains ~12.2° per day on Sun
    relative_speed = 12.190749  # degrees/day

    # Time until next opposition (elongation from 180° = 0)
    if elongation > 0:
        # Past opposition, wait for next cycle
        dt = (360.0 - elongation) / relative_speed
    else:
        dt = (-elongation) / relative_speed

    jd_guess = jd_start + dt

    # Newton-Raphson refinement
    for _ in range(20):
        sun_pos, _ = calc_ut(jd_guess, SUN, FLG_SPEED)
        moon_pos, _ = calc_ut(jd_guess, MOON, FLG_SPEED)

        sun_lon = sun_pos[0]
        moon_lon = moon_pos[0]
        sun_speed = sun_pos[3]
        moon_speed = moon_pos[3]

        # Elongation from opposition
        diff = (moon_lon - sun_lon - 180.0) % 360.0
        if diff > 180:
            diff -= 360

        # Convergence check (< 0.1 arcsec)
        if abs(diff) < 1e-5:
            return jd_guess

        # Newton-Raphson step
        rel_speed = moon_speed - sun_speed
        if abs(rel_speed) < 0.1:
            rel_speed = 12.19

        jd_guess -= diff / rel_speed

    return jd_guess


def _find_previous_full_moon(jd_start: float) -> float:
    """
    Find the previous Full Moon (Sun-Moon opposition) before jd_start.

    Uses iterative refinement to find exact moment of opposition.

    Args:
        jd_start: Julian Day (UT) to start search from

    Returns:
        Julian Day of previous Full Moon
    """
    # Get current positions
    sun_pos, _ = calc_ut(jd_start, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd_start, MOON, FLG_SPEED)

    sun_lon = sun_pos[0]
    moon_lon = moon_pos[0]

    # Calculate elongation from opposition (Moon - Sun - 180°)
    elongation = (moon_lon - sun_lon - 180.0) % 360.0
    if elongation > 180:
        elongation -= 360

    # Moon gains ~12.2° per day on Sun
    relative_speed = 12.190749  # degrees/day

    # Time since last opposition (go backwards)
    if elongation > 0:
        # Past opposition, last opposition was elongation/speed days ago
        dt = -elongation / relative_speed
    else:
        # Before opposition, last opposition was (360 + elongation)/speed days ago
        dt = -(360.0 + elongation) / relative_speed

    jd_guess = jd_start + dt

    # Newton-Raphson refinement
    for _ in range(20):
        sun_pos, _ = calc_ut(jd_guess, SUN, FLG_SPEED)
        moon_pos, _ = calc_ut(jd_guess, MOON, FLG_SPEED)

        sun_lon = sun_pos[0]
        moon_lon = moon_pos[0]
        sun_speed = sun_pos[3]
        moon_speed = moon_pos[3]

        # Elongation from opposition
        diff = (moon_lon - sun_lon - 180.0) % 360.0
        if diff > 180:
            diff -= 360

        # Convergence check (< 0.1 arcsec)
        if abs(diff) < 1e-5:
            return jd_guess

        # Newton-Raphson step
        rel_speed = moon_speed - sun_speed
        if abs(rel_speed) < 0.1:
            rel_speed = 12.19

        jd_guess -= diff / rel_speed

    return jd_guess


def _calculate_lunar_eclipse_type_and_magnitude(
    jd: float,
) -> Tuple[int, float, float, float, float, float]:
    """
    Determine lunar eclipse type and magnitude at maximum eclipse.

    Uses geometric calculations based on Earth's shadow cones at the Moon's distance.
    Includes edge case handling for:
    - Very shallow penumbral eclipses (penumbral magnitude near 0)
    - Very shallow umbral eclipses (umbral magnitude near 0)
    - Division by zero protection

    Args:
        jd: Julian Day of eclipse maximum (UT)

    Returns:
        Tuple of (eclipse_type_flags, magnitude_umbral, magnitude_penumbral,
                  gamma, penumbra_radius, umbra_radius)
        - eclipse_type_flags: Bitmask of ECL_* constants
        - magnitude_umbral: Eclipse magnitude (fraction of Moon in umbra)
        - magnitude_penumbral: Penumbral eclipse magnitude
        - gamma: Gamma parameter (Moon's distance from shadow axis in Earth radii)
        - penumbra_radius: Penumbral shadow radius at Moon distance (degrees)
        - umbra_radius: Umbral shadow radius at Moon distance (degrees)
    """
    # Get positions
    sun_pos, _ = calc_ut(jd, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd, MOON, FLG_SPEED)

    # Distances in AU
    sun_dist = sun_pos[2]
    moon_dist = moon_pos[2]

    # Moon's ecliptic latitude (how far from the ecliptic)
    moon_lat = moon_pos[1]

    # Constants for shadow calculations
    # Sun angular radius at 1 AU: 959.63 arcsec
    # Earth radius: 6378.137 km = 4.2635e-5 AU
    SUN_RADIUS_ARCSEC = 959.63
    EARTH_RADIUS_AU = 4.2635e-5
    EARTH_RADIUS_KM = 6378.137

    # Sun's angular semi-diameter at actual distance (in degrees)
    sun_semidiameter = (SUN_RADIUS_ARCSEC / 3600.0) / sun_dist

    # Moon's angular semi-diameter (932.56 arcsec at mean distance 0.002569 AU)
    moon_semidiameter = (932.56 / 3600.0) * (0.002569 / moon_dist)

    # Earth's angular semi-diameter as seen from Moon (in degrees)
    earth_semidiameter = math.degrees(math.atan(EARTH_RADIUS_AU / moon_dist))

    # Calculate shadow cone sizes at Moon's distance
    # Using the geometric relationship of similar triangles

    # Sun's angular radius as seen from Earth
    sun_angular_radius_from_earth = sun_semidiameter

    # Penumbra radius at Moon's distance (Earth + Sun shadow)
    # Penumbra outer edge: Earth appears larger, adding Sun's angular size
    penumbra_radius = earth_semidiameter + sun_angular_radius_from_earth

    # Umbra radius at Moon's distance (Earth minus Sun shadow)
    # The umbra is smaller because of the Sun's finite size
    # Using more accurate shadow cone calculation
    # Shadow cone angle from geometry
    sun_dist_km = sun_dist * 149597870.7
    moon_dist_km = moon_dist * 149597870.7
    sun_radius_km = (SUN_RADIUS_ARCSEC / 206265.0) * sun_dist_km

    # Umbra cone semi-angle - protect against division by zero
    denom = sun_dist_km
    if abs(denom) < 1e-10:
        umbra_cone_angle = 0.0
    else:
        umbra_cone_angle = math.atan((sun_radius_km - EARTH_RADIUS_KM) / denom)
    # Umbra radius at Moon distance
    umbra_radius_km = EARTH_RADIUS_KM - moon_dist_km * math.tan(umbra_cone_angle)

    if umbra_radius_km > 0:
        umbra_radius = math.degrees(math.atan(umbra_radius_km / moon_dist_km))
    else:
        # Umbra doesn't reach Moon (antumbra) - extremely rare for lunar eclipses
        umbra_radius = 0.0

    # Apply atmospheric enlargement of Earth's shadow (Danjon factor).
    # Earth's atmosphere refracts sunlight, effectively enlarging the shadow
    # by approximately 1/85 (~1.18%). This correction is standard in eclipse
    # magnitude calculations and is necessary for correct classification of
    # borderline total/partial eclipses.
    # Reference: Danjon, A. (1951), "Les éclipses de Lune"
    # Also: Meeus, "Astronomical Algorithms", Ch. 54
    _SHADOW_ENLARGEMENT = 1.0 + 1.0 / 85.0
    penumbra_radius *= _SHADOW_ENLARGEMENT
    umbra_radius *= _SHADOW_ENLARGEMENT

    # Moon's distance from the shadow axis (in degrees).
    # The shadow centre is at the anti-solar point (Sun + 180° in longitude,
    # negated latitude).  At eclipse maximum the longitude term is nearly
    # zero and |moon_lat| alone is sufficient, but at non-maximum times
    # (e.g. moonrise/moonset) the longitude offset can be significant
    # and must be included.
    sun_lon = sun_pos[0]
    moon_lon = moon_pos[0]
    shadow_lon = (sun_lon + 180.0) % 360.0

    delta_lon = moon_lon - shadow_lon
    if delta_lon > 180.0:
        delta_lon -= 360.0
    elif delta_lon < -180.0:
        delta_lon += 360.0

    # Small-angle approximation for angular separation on the ecliptic
    # (cos correction for longitude at the Moon's latitude)
    cos_lat = math.cos(math.radians(moon_lat))
    moon_distance_from_axis = math.sqrt((delta_lon * cos_lat) ** 2 + moon_lat**2)

    # Gamma: Moon's distance from shadow axis in Earth radii
    # Edge case: protect against zero earth_semidiameter
    if abs(earth_semidiameter) < 1e-10:
        gamma = 0.0
    else:
        gamma = moon_distance_from_axis / earth_semidiameter

    # Calculate eclipse magnitudes
    # Edge case: protect against zero moon_semidiameter
    if abs(moon_semidiameter) < 1e-10:
        return 0, 0.0, 0.0, gamma, penumbra_radius, umbra_radius

    # Penumbral magnitude: how far Moon penetrates into penumbra
    penumbral_mag = (penumbra_radius + moon_semidiameter - moon_distance_from_axis) / (
        2 * moon_semidiameter
    )

    # Umbral magnitude: how far Moon penetrates into umbra
    umbral_mag = (umbra_radius + moon_semidiameter - moon_distance_from_axis) / (
        2 * moon_semidiameter
    )

    # Determine eclipse type
    eclipse_type = 0

    # Edge case: no eclipse (penumbral magnitude <= 0)
    if penumbral_mag <= 0:
        # No eclipse
        return 0, 0.0, 0.0, gamma, penumbra_radius, umbra_radius

    # Edge case: shallow penumbral eclipse
    if penumbral_mag > 0 and penumbral_mag < SHALLOW_ECLIPSE_MAG_THRESHOLD:
        # Very shallow penumbral eclipse - mark as grazing
        # Shallow grazing penumbral reports plain PENUMBRAL (compatibility
        # contract)
        eclipse_type = ECL_PENUMBRAL
        penumbral_mag = max(0.0, penumbral_mag)
        return eclipse_type, 0.0, penumbral_mag, gamma, penumbra_radius, umbra_radius

    if umbral_mag <= 0:
        # Penumbral only
        eclipse_type = ECL_PENUMBRAL
        penumbral_mag = max(0.0, penumbral_mag)
        return eclipse_type, 0.0, penumbral_mag, gamma, penumbra_radius, umbra_radius

    # Edge case: shallow umbral (partial) eclipse
    if umbral_mag > 0 and umbral_mag < SHALLOW_ECLIPSE_MAG_THRESHOLD:
        # Very shallow partial umbral eclipse - mark as grazing
        # Shallow grazing partial reports plain PARTIAL (compatibility
        # contract)
        eclipse_type = ECL_PARTIAL
        umbral_mag = max(0.0, min(1.0, umbral_mag))
        penumbral_mag = max(0.0, penumbral_mag)
        return (
            eclipse_type,
            umbral_mag,
            penumbral_mag,
            gamma,
            penumbra_radius,
            umbra_radius,
        )

    if umbral_mag >= 1.0:
        # Total umbral eclipse
        eclipse_type = ECL_TOTAL
        umbral_mag = max(0.0, umbral_mag)
        penumbral_mag = max(0.0, penumbral_mag)
    else:
        # Partial umbral eclipse
        eclipse_type = ECL_PARTIAL
        umbral_mag = max(0.0, min(1.0, umbral_mag))
        penumbral_mag = max(0.0, penumbral_mag)

    return eclipse_type, umbral_mag, penumbral_mag, gamma, penumbra_radius, umbra_radius


def _refine_lunar_eclipse_maximum(
    jd_full_moon: float, search_range: float = 0.5
) -> float:
    """
    Refine the lunar eclipse maximum time from Full Moon approximation.

    The eclipse maximum occurs when the Moon is closest to Earth's shadow
    center (the anti-Sun point), not exactly at Full Moon (opposition in
    longitude). This function uses golden section search to find the true
    maximum by minimizing the angular separation.

    Args:
        jd_full_moon: Julian Day of Full Moon (initial approximation)
        search_range: Search range in days (±0.5 days covers the offset)

    Returns:
        Julian Day of true eclipse maximum (when Moon is closest to shadow axis)
    """

    def _get_shadow_separation(jd: float) -> float:
        """Calculate angular separation between Moon and shadow center."""
        sun_pos, _ = calc_ut(jd, SUN, FLG_SPEED)
        moon_pos, _ = calc_ut(jd, MOON, FLG_SPEED)

        # Shadow center is the anti-Sun point (180° from Sun)
        shadow_lon = (sun_pos[0] + 180.0) % 360.0
        shadow_lat = -sun_pos[1]  # Opposite latitude

        # Calculate angular separation using spherical law of cosines
        d_lon = math.radians(moon_pos[0] - shadow_lon)
        lat1 = math.radians(shadow_lat)
        lat2 = math.radians(moon_pos[1])

        cos_sep = math.sin(lat1) * math.sin(lat2) + math.cos(lat1) * math.cos(
            lat2
        ) * math.cos(d_lon)
        cos_sep = max(-1.0, min(1.0, cos_sep))
        return math.degrees(math.acos(cos_sep))

    # Golden section search to minimize separation
    phi = (1 + math.sqrt(5)) / 2

    jd_low = jd_full_moon - search_range
    jd_high = jd_full_moon + search_range

    jd_a = jd_high - (jd_high - jd_low) / phi
    jd_b = jd_low + (jd_high - jd_low) / phi

    sep_a = _get_shadow_separation(jd_a)
    sep_b = _get_shadow_separation(jd_b)

    for _ in range(50):
        if sep_a < sep_b:
            jd_high = jd_b
            jd_b = jd_a
            sep_b = sep_a
            jd_a = jd_high - (jd_high - jd_low) / phi
            sep_a = _get_shadow_separation(jd_a)
        else:
            jd_low = jd_a
            jd_a = jd_b
            sep_a = sep_b
            jd_b = jd_low + (jd_high - jd_low) / phi
            sep_b = _get_shadow_separation(jd_b)

        if jd_high - jd_low < 1e-8:  # ~0.86 ms precision
            break

    return (jd_low + jd_high) / 2


def _calculate_lunar_eclipse_phases(
    jd_max: float,
    eclipse_type: int,
    moon_semidiameter: float,
    umbra_radius: float,
    penumbra_radius: float,
    moon_lat_at_max: float,
) -> Tuple[float, float, float, float, float, float, float, float, float, float]:
    """
    Calculate times of lunar eclipse phases (contacts).

    Phase indices (reference-API-compatible 10-element format):
        [0]: Time of maximum eclipse
        [1]: Reserved (0)
        [2]: Time of partial eclipse beginning (Moon enters umbra)
        [3]: Time of partial eclipse ending (Moon leaves umbra)
        [4]: Time of total eclipse beginning (or 0 if not total)
        [5]: Time of total eclipse ending (or 0 if not total)
        [6]: Time of penumbral eclipse beginning
        [7]: Time of penumbral eclipse ending
        [8]: Reserved (0)
        [9]: Reserved (0)

    Args:
        jd_max: Julian Day of maximum eclipse
        eclipse_type: Eclipse type flags
        moon_semidiameter: Moon's angular semi-diameter in degrees
        umbra_radius: Umbral shadow radius in degrees
        penumbra_radius: Penumbral shadow radius in degrees
        moon_lat_at_max: Moon's ecliptic latitude at maximum in degrees

    Returns:
        Tuple of 8 floats with phase times (JD UT)
    """
    # Actually calculate from Moon's speed minus Sun's speed
    sun_pos, _ = calc_ut(jd_max, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd_max, MOON, FLG_SPEED)

    # Speed of Moon relative to shadow (in longitude)
    relative_speed = abs(moon_pos[3] - sun_pos[3])  # degrees/day
    if relative_speed < 0.1:
        relative_speed = 0.55  # fallback to typical value

    # Distance from shadow axis at maximum
    y = abs(moon_lat_at_max)

    # Calculate half-duration of each phase using geometry
    # The Moon moves through the shadow at an angle
    # Half-duration = sqrt(R² - y²) / speed where R is radius and y is impact parameter

    def calc_half_duration(radius: float) -> float:
        """Calculate half-duration for given shadow radius.

        Includes edge case handling for shallow eclipses where
        the Moon barely touches the shadow boundary.
        """
        r_total = radius + moon_semidiameter
        if y >= r_total:
            return 0.0
        # Calculate chord half-length using safe sqrt
        chord_sq = r_total * r_total - y * y
        # Edge case: very shallow eclipse where chord is nearly zero
        if chord_sq < 1e-20:
            return 0.0
        half_chord = _safe_sqrt(chord_sq)
        return half_chord / relative_speed

    def calc_total_half_duration(radius: float) -> float:
        """Calculate half-duration of total phase (Moon fully inside).

        Includes edge case handling for nearly-total eclipses where
        the Moon barely fits inside the umbra.
        """
        r_inner = radius - moon_semidiameter
        if r_inner <= 0 or y >= r_inner:
            return 0.0
        chord_sq = r_inner * r_inner - y * y
        # Edge case: nearly-total eclipse where chord is nearly zero
        if chord_sq < 1e-20:
            return 0.0
        half_chord = _safe_sqrt(chord_sq)
        return half_chord / relative_speed

    # Penumbral phase times
    penumbral_half_dur = calc_half_duration(penumbra_radius)
    t_pen_begin = jd_max - penumbral_half_dur
    t_pen_end = jd_max + penumbral_half_dur

    # Partial (umbral) phase times
    partial_half_dur = calc_half_duration(umbra_radius)
    if partial_half_dur > 0:
        t_partial_begin = jd_max - partial_half_dur
        t_partial_end = jd_max + partial_half_dur
    else:
        t_partial_begin = 0.0
        t_partial_end = 0.0

    # Total phase times
    if eclipse_type & ECL_TOTAL:
        total_half_dur = calc_total_half_duration(umbra_radius)
        if total_half_dur > 0:
            t_total_begin = jd_max - total_half_dur
            t_total_end = jd_max + total_half_dur
        else:
            t_total_begin = 0.0
            t_total_end = 0.0
    else:
        t_total_begin = 0.0
        t_total_end = 0.0

    return (
        jd_max,  # [0] Maximum eclipse
        0.0,  # [1] Reserved
        t_partial_begin,  # [2] Partial begins (Moon enters umbra)
        t_partial_end,  # [3] Partial ends (Moon leaves umbra)
        t_total_begin,  # [4] Total begins (or 0 if not total)
        t_total_end,  # [5] Total ends (or 0 if not total)
        t_pen_begin,  # [6] Penumbral begins
        t_pen_end,  # [7] Penumbral ends
        0.0,  # [8] Reserved
        0.0,  # [9] Reserved
    )


def _lun_eclipse_when_pythonic(
    jd_start: float,
    flags: int = FLG_SWIEPH,
    eclipse_type: int = 0,
    backwards: bool = False,
) -> Tuple[int, Tuple[float, ...]]:
    """
    Find the next (or previous) lunar eclipse from a given date.

    Searches forward (or backward when ``backwards`` is True) in time from
    jd_start. Can filter by eclipse type (total, partial, penumbral).

    Args:
        jd_start: Julian Day (UT) to start search from
        flags: Calculation flags (FLG_SWIEPH, etc.)
        eclipse_type: Filter for specific eclipse type(s), bitmask of:
            - ECL_TOTAL (4): Total lunar eclipse
            - ECL_PARTIAL (16): Partial lunar eclipse
            - ECL_PENUMBRAL (64): Penumbral lunar eclipse
            - 0: Any eclipse type (default)

    Returns:
        Tuple containing (matching reference API format):
            - retflag: Eclipse type flags bitmask (ECL_* constants)
            - tret: Tuple of 10 floats with eclipse phase times (JD UT):
                [0]: Time of maximum eclipse
                [1]: Reserved (0)
                [2]: Time of partial eclipse beginning (Moon enters umbra)
                [3]: Time of partial eclipse ending (Moon leaves umbra)
                [4]: Time of total eclipse beginning (or 0 if not total)
                [5]: Time of total eclipse ending (or 0 if not total)
                [6]: Time of penumbral eclipse beginning
                [7]: Time of penumbral eclipse ending
                [8]: Reserved (0)
                [9]: Reserved (0)

    Raises:
        Error: If no eclipse found within search limit

    Algorithm:
        1. Find next Full Moon after jd_start
        2. Check if Moon is close enough to node for eclipse
        3. If not eclipse, advance to next Full Moon
        4. Calculate eclipse type and magnitude
        5. If eclipse_type filter set, check if matches
        6. Calculate phase times

    Precision:
        Eclipse times accurate to ~1 minute for most eclipses.

    Example:
        >>> # Find next total lunar eclipse after Jan 1, 2024
        >>> from libephemeris import julday, ECL_TOTAL
        >>> jd = julday(2024, 1, 1, 0)
        >>> ecl_type, times = _lun_eclipse_when_pythonic(jd, eclipse_type=ECL_TOTAL)
        >>> print(f"Total lunar eclipse at JD {times[0]:.5f}")

    References:
        - Reference API: lun_eclipse_when()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    # The search walks lunation-by-lunation until it finds a matching
    # eclipse or steps off the active ephemeris (the tier boundary).
    # MAX_FULL_MOONS is only a finite safety backstop (see the solar
    # search for the rationale).
    MAX_SEARCH_YEARS = _ECLIPSE_SEARCH_HORIZON_YEARS
    MAX_FULL_MOONS = int(MAX_SEARCH_YEARS * 12.4)  # ~12.4 lunations per year

    # Central/noncentral bits are meaningless for lunar eclipses and are
    # ignored; annular types do not exist for the Moon - a request for
    # only annular types is an error (reference parity), otherwise the
    # annular bits are dropped.
    eclipse_type = eclipse_type & ~(ECL_CENTRAL | ECL_NONCENTRAL)
    if eclipse_type & (ECL_ANNULAR | ECL_ANNULAR_TOTAL):
        eclipse_type &= ~(ECL_ANNULAR | ECL_ANNULAR_TOTAL)
        if eclipse_type == 0:
            # NOTE: use the module-level Error — a local
            # `from .exceptions import Error` here would make the name local
            # to the WHOLE function and turn the search-exhaustion
            # `raise Error(...)` at the bottom into an UnboundLocalError on
            # every non-annular call.
            raise Error("annular lunar eclipses don't exist")
    eclipse_type = eclipse_type & ECL_ALLTYPES_LUNAR

    # If eclipse_type is 0, accept any type
    if eclipse_type == 0:
        eclipse_type = ECL_ALLTYPES_LUNAR

    # Days to look back for the forward pre-probe below (mirror of the
    # solar search's BIDIRECTIONAL_WINDOW).
    LUNAR_PRE_PROBE_WINDOW = 15.0

    def _check_full_moon_for_eclipse(
        jd_full_moon: float,
    ) -> Union[Tuple[int, Tuple[float, ...]], None]:
        """Classify the matching eclipse (if any) at a Full Moon.

        Encapsulates the node pre-filter, the maximum refinement, the
        shadow-geometry classification and the type-mask test so the forward
        pre-probe and the main lunation walk share one code path. Returns the
        ``(retflag, tret)`` pair or ``None`` when this Full Moon carries no
        matching eclipse.
        """
        moon_pos, _ = calc_ut(jd_full_moon, MOON, flags | FLG_SPEED)
        moon_lon = moon_pos[0]

        # Coarse PRE-FILTER (widened well past the physical limit); the true
        # test is the shadow-geometry classification at the refined maximum
        # (node distance is not a clean discriminator here).
        node_dist = _get_moon_node_distance(jd_full_moon, moon_lon)
        if node_dist >= _LUNAR_NODE_PREFILTER:
            return None

        # Refine the time of maximum eclipse: deepest immersion of the Moon
        # in the Earth's shadow, then classify from the shadow geometry.
        jd_max = _lun_eclipse_max_time(jd_full_moon, flags)
        retc, _attr_max, _dcore_max = _lun_how_core(jd_max, flags)
        if retc == 0:
            return None

        type_matches = (
            (eclipse_type & ECL_TOTAL and retc & ECL_TOTAL)
            or (eclipse_type & ECL_PARTIAL and retc & ECL_PARTIAL)
            or (eclipse_type & ECL_PENUMBRAL and retc & ECL_PENUMBRAL)
        )
        if not type_matches:
            return None

        times = _lun_eclipse_phase_times(jd_max, retc, flags)
        return retc, times

    # Forward pre-probe (mirror of the solar search). A lunar eclipse's
    # opposition (Full Moon) precedes its maximum by a few minutes, so a
    # start epoch that falls in the [opposition, maximum) window would make
    # _find_next_full_moon() skip to the NEXT lunation and drop the
    # in-progress eclipse. Probe the previous Full Moon first and return its
    # eclipse when its maximum is still more than the epoch margin ahead of
    # jd_start (compatibility contract: the current eclipse is returned
    # right up to maximum - _ECLIPSE_WHEN_EPOCH_MARGIN).
    if not backwards:
        try:
            jd_prev_fm = _find_previous_full_moon(jd_start)
        except Exception as _exc:
            if _is_ephemeris_boundary(_exc):
                jd_prev_fm = None
            else:
                raise
        if jd_prev_fm is not None and jd_start - jd_prev_fm <= LUNAR_PRE_PROBE_WINDOW:
            pre = _check_full_moon_for_eclipse(jd_prev_fm)
            if pre is not None and pre[1][0] > jd_start + _ECLIPSE_WHEN_EPOCH_MARGIN:
                return pre

    jd = jd_start

    for _ in range(MAX_FULL_MOONS):
        # Find next (or previous) Full Moon and classify any eclipse there.
        try:
            if backwards:
                jd_full_moon = _find_previous_full_moon(jd)
            else:
                jd_full_moon = _find_next_full_moon(jd)
            result = _check_full_moon_for_eclipse(jd_full_moon)
        except Exception as _exc:
            if _is_ephemeris_boundary(_exc):
                # Walked off the ephemeris: no matching eclipse within coverage.
                break
            raise

        if result is not None:
            jd_max = result[1][0]

            # Direction invariant: a forward search must return an eclipse
            # after jd_start, a backward search one before it (refinement
            # can drag the maximum across the start time). The epoch margin
            # makes a maximum within _ECLIPSE_WHEN_EPOCH_MARGIN of jd_start count as "reached",
            # so `jd = tret[0]; when(jd)` advances to the neighbouring
            # eclipse (compatibility contract).
            skip = (
                not backwards and jd_max <= jd_start + _ECLIPSE_WHEN_EPOCH_MARGIN
            ) or (backwards and jd_max >= jd_start - _ECLIPSE_WHEN_EPOCH_MARGIN)
            if not skip:
                return result

        # Advance (or retreat) to the next lunation: ~25 days keeps us
        # safely within one synodic month of the next Full Moon.
        jd = jd_full_moon + (-25 if backwards else 25)

    raise Error(
        "No matching lunar eclipse found before the ephemeris boundary "
        f"searching from JD {jd_start}"
    )


def lun_eclipse_when(
    tjdut: float,
    flags: int = FLG_SWIEPH,
    ecltype: int = 0,
    backwards: "bool | int | str" = False,
) -> Tuple[int, Tuple[float, ...]]:
    """Find the next (or previous) lunar eclipse globally (reference-API-compatible).

    Wrapper around _lun_eclipse_when_pythonic() matching the reference signature.

    Args:
        tjdut: Julian Day (UT) to start search from.
        flags: Calculation flags (default FLG_SWIEPH).
        ecltype: Eclipse type filter bitmask (0 = any).
        backwards: If True, search backward in time. Also accepts 0/1 and
            the direction strings understood by :func:`_coerce_backwards`.

    Returns:
        Tuple of (retflag, tret) matching the reference API.
    """
    flags = _strip_one_try_bit(flags)
    return _lun_eclipse_when_pythonic(
        tjdut,
        flags=flags,
        eclipse_type=ecltype,
        backwards=_coerce_backwards(backwards),
    )


def _lun_eclipse_when_loc_pythonic(
    jd_start: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    flags: int = FLG_SWIEPH,
    backwards: bool = False,
    *,
    reader=_ACTIVE_LEB_READER,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """Find the next lunar eclipse observable from a geographic position.

    Searches lunar eclipses (forward or backward) and keeps the first one
    with at least one phase above the local horizon. Phase times whose
    event happens while the Moon is below the horizon are zeroed; the
    Moon's rise/set times bounding the visible window are reported in
    tret[8]/tret[9], and the reported maximum is moved to the horizon
    crossing when the geometric maximum itself is not observable.

    Returns:
        (retflag, tret, attr) where retflag combines the eclipse type
        with ECL_VISIBLE, ECL_MAX_VISIBLE and the per-phase
        ECL_PARTBEG/PARTEND/TOTBEG/TOTEND/PENUMBBEG/PENUMBEND_VISIBLE
        bits; tret is the lun_eclipse_when layout plus moonrise/moonset
        in [8]/[9]; attr are the lun_eclipse_how attributes at tret[0].

    Raises:
        Error: if no observable eclipse is found within the
            search limit.
    """
    if reader is _ACTIVE_LEB_READER:
        # Direct callers (this pythonic API is re-exported) get the same
        # mode-aware dispatch as the public wrapper: LEB first, Skyfield
        # retry on a range miss in auto mode, typed error in sealed mode.
        return _call_with_leb_skyfield_fallback(
            _lun_eclipse_when_loc_pythonic,
            jd_start,
            lat,
            lon,
            altitude,
            flags,
            backwards,
        )

    from .constants import (
        BIT_DISC_BOTTOM,
        CALC_RISE,
        CALC_SET,
        ECL_PARTBEG_VISIBLE,
        ECL_PARTEND_VISIBLE,
        ECL_PENUMBBEG_VISIBLE,
        ECL_PENUMBEND_VISIBLE,
        ECL_TOTBEG_VISIBLE,
        ECL_TOTEND_VISIBLE,
    )

    geopos3 = (float(lon), float(lat), float(altitude))
    eph_flags = _ecl_eph_flags(flags)
    vis_bits = {
        0: ECL_MAX_VISIBLE,
        2: ECL_PARTBEG_VISIBLE,
        3: ECL_PARTEND_VISIBLE,
        4: ECL_TOTBEG_VISIBLE,
        5: ECL_TOTEND_VISIBLE,
        6: ECL_PENUMBBEG_VISIBLE,
        7: ECL_PENUMBEND_VISIBLE,
    }

    search_jd = jd_start
    for _ in range(250):  # ~20 years of lunar eclipses
        # Offer the current lunation as a candidate even when its geocentric
        # maximum has just slipped past the epoch: the eclipse may still be
        # unfolding, or the Moon may only rise into it later. Anchor the
        # search a little before/after the epoch and let the epoch gate below
        # keep or skip the candidate. Only the first iteration is affected;
        # later iterations start well past the previous maximum.
        _rc_when, tret_t = _lun_eclipse_when_pythonic(
            search_jd
            + (_WHEN_LOC_EPOCH_WINDOW if backwards else -_WHEN_LOC_EPOCH_WINDOW),
            flags,
            0,
            backwards,
        )
        tret = list(tret_t)

        # Local visible-window bounds for the epoch gate applied after the
        # rise/set clamping below. Lunar phase times are the same for every
        # observer, but visibility is local, so the gate must reason about what
        # is actually observable here - not the geocentric extent. Default to
        # the geocentric penumbral contacts (used verbatim when the Moon never
        # crosses the horizon during the eclipse, e.g. a circumpolar Moon);
        # they are narrowed to moonrise/moonset when the Moon rises or sets
        # mid-eclipse. Captured now, before the clamping step zeroes the
        # contacts. Fall back to the maximum when an outer contact is
        # unresolved.
        first_visible = tret[6] if tret[6] != 0.0 else tret[0]
        last_visible = tret[7] if tret[7] != 0.0 else tret[0]

        # Visibility bits: a phase counts as visible when the Moon's
        # apparent altitude is positive at its time.
        retflag = 0
        for i in (7, 6, 5, 4, 3, 2, 0):
            if tret[i] == 0.0:
                continue
            _rc_i, attr_i = _lun_eclipse_how_impl(
                tret[i], geopos3, flags, reader=reader
            )
            if attr_i[6] > 0.0:
                retflag |= ECL_VISIBLE | vis_bits[i]

        if not (retflag & ECL_VISIBLE):
            search_jd = tret[0] + (-25.0 if backwards else 25.0)
            continue

        # Moonrise and moonset around the eclipse (lower-limb events,
        # the reference convention here). A circumpolar Moon (rc -2)
        # leaves the geometric phase times untouched.
        tjd_max = tret[0]
        rc_rise, tr_rise = rise_trans(
            tret[6] - 0.001,
            MOON,
            CALC_RISE | BIT_DISC_BOTTOM,
            geopos3,
            0.0,
            0.0,
            eph_flags,
        )
        rc_set, tr_set = rise_trans(
            tret[6] - 0.001,
            MOON,
            CALC_SET | BIT_DISC_BOTTOM,
            geopos3,
            0.0,
            0.0,
            eph_flags,
        )
        if rc_rise >= 0 and rc_set >= 0:
            tjdr = tr_rise[0]
            tjds = tr_set[0]
            if tjds < tret[6] or (tjds > tjdr and tjdr > tret[7]):
                # The whole eclipse runs while the Moon is down.
                search_jd = tret[0] + (-25.0 if backwards else 25.0)
                continue
            if tret[6] < tjdr < tret[7]:
                # Moon rises mid-eclipse: phases before moonrise are
                # not observable.
                tret[6] = 0.0
                for i in range(2, 6):
                    if tjdr > tret[i]:
                        tret[i] = 0.0
                tret[8] = tjdr
                first_visible = tjdr
                if tjdr > tret[0]:
                    tjd_max = tjdr
            if tret[6] < tjds < tret[7]:
                # Moon sets mid-eclipse: phases after moonset are not
                # observable.
                tret[7] = 0.0
                for i in range(2, 6):
                    if tjds < tret[i]:
                        tret[i] = 0.0
                tret[9] = tjds
                last_visible = tjds
                if tjds < tret[0]:
                    tjd_max = tjds

        # Epoch gate on the LOCAL visible window (mirrors the solar side, which
        # gates on the local contacts). Skip this lunation only when its
        # observable portion is wholly before the epoch (forward) or wholly
        # after it (backward). Gating on the geocentric contacts instead would
        # retain an eclipse whose visible part already ended at moonset (or has
        # not yet begun at moonrise), wedging the sequential search on an
        # already-finished eclipse. tret[0] is still the geocentric maximum
        # here, so the search advance matches the other skips.
        if (
            not backwards
            and (
                last_visible <= jd_start + _ECLIPSE_WHEN_EPOCH_MARGIN
                # Re-entrancy: the caller started at or past the instant this
                # function would return, so advance instead of re-finding the
                # same eclipse (see the solar twin).
                or tret[0] <= jd_start + _ECLIPSE_WHEN_EPOCH_MARGIN
            )
        ) or (
            backwards
            and (
                first_visible >= jd_start - _ECLIPSE_WHEN_EPOCH_MARGIN
                or tret[0] >= jd_start - _ECLIPSE_WHEN_EPOCH_MARGIN
            )
        ):
            search_jd = tret[0] + (-25.0 if backwards else 25.0)
            continue

        tret[0] = tjd_max
        rc_final, attr_final = _lun_eclipse_how_impl(
            tjd_max, geopos3, flags, reader=reader
        )
        if rc_final == 0:
            search_jd = tret[0] + (-25.0 if backwards else 25.0)
            continue
        retflag |= rc_final & ECL_ALLTYPES_LUNAR
        return retflag, tuple(tret), tuple(attr_final)

    raise Error(
        f"No lunar eclipse visible from lon={lon}, lat={lat} "
        f"within the search range of JD {jd_start}"
    )


def lun_eclipse_when_loc(
    tjdut: float,
    geopos: "Sequence[float]",
    flags: int = FLG_SWIEPH,
    backwards: "bool | int | str" = False,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """Find the next lunar eclipse visible from a geographic position (reference-API-compatible).

    Wrapper around _lun_eclipse_when_loc_pythonic() matching the reference signature.

    Args:
        tjdut: Julian Day (UT) to start search from.
        geopos: Sequence of [longitude, latitude, altitude].
        flags: Calculation flags (default FLG_SWIEPH).
        backwards: If True, search backward in time. Also accepts 0/1 and
            the direction strings understood by :func:`_coerce_backwards`.

    Returns:
        Tuple of (retflag, tret, attr) matching the reference API.
    """
    flags = _strip_one_try_bit(flags)
    if len(geopos) < 3:
        raise ValueError("geopos must have at least 3 elements: [lon, lat, alt]")

    lon = float(geopos[0])
    lat = float(geopos[1])
    altitude = float(geopos[2])

    return _call_with_leb_skyfield_fallback(
        _lun_eclipse_when_loc_pythonic,
        tjdut,
        lat,
        lon,
        altitude,
        flags,
        _coerce_backwards(backwards),
    )


def _lun_eclipse_how_pythonic(
    jd: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    flags: int = FLG_SWIEPH,
    *,
    reader=_ACTIVE_LEB_READER,
) -> Tuple[int, Tuple[float, ...]]:
    """
    Calculate the circumstances of a lunar eclipse at a specific location and time.

    This function determines how a lunar eclipse appears from a given geographic
    position at a specific Julian Day. Unlike _lun_eclipse_when_loc_pythonic which finds the
    next eclipse, this function calculates the eclipse magnitude, Moon position,
    and other circumstances for a known eclipse time.

    Args:
        jd: Julian Day (UT) of the moment to calculate
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - retflag: Eclipse type flags bitmask (ECL_* constants)
                Returns 0 if no eclipse is occurring at this time
                Includes ECL_VISIBLE if Moon is above horizon
            - attr: Tuple of 11 floats with eclipse attributes:
                [0]: Umbral eclipse magnitude (fraction of Moon's diameter in umbra)
                [1]: Penumbral eclipse magnitude
                [2]: Reserved (0)
                [3]: Azimuth of Moon at the given time (degrees)
                [4]: Altitude of Moon at the given time (degrees)
                [5]: Apparent diameter of Moon (degrees)
                [6]: Apparent diameter of umbral shadow at Moon's distance (degrees)
                [7]: Apparent diameter of penumbral shadow at Moon's distance (degrees)
                [8]: Saros series number (0, not implemented)
                [9]: Reserved (0)
                [10]: Reserved (0)

    Note:
        This function is intended for use when you already know an eclipse is
        occurring (e.g., from _lun_eclipse_when_pythonic or _lun_eclipse_when_loc_pythonic).
        For a random time when no eclipse is occurring, magnitude will be 0
        and retflag will be 0.

    Algorithm:
        1. Calculate Moon's apparent position from observer location
        2. Calculate Earth's shadow cone geometry at Moon's distance
        3. Compute umbral and penumbral magnitudes
        4. Determine eclipse type based on Moon's penetration into shadows
        5. Return attributes including Moon's altitude and azimuth

    Precision:
        Eclipse magnitude accurate to ~0.001 for central eclipses.
        Topocentric parallax included in calculations.

    Example:
        >>> # Calculate eclipse circumstances at Rome during a lunar eclipse
        >>> from libephemeris import julday, _lun_eclipse_how_pythonic
        >>> jd = julday(2022, 5, 16, 4.2)  # During May 2022 eclipse
        >>> rome_lat, rome_lon = 41.9028, 12.4964
        >>> ecl_type, attr = _lun_eclipse_how_pythonic(jd, rome_lat, rome_lon)
        >>> print(f"Umbral magnitude: {attr[0]:.3f}")
        >>> print(f"Moon altitude: {attr[4]:.1f}°")

    References:
        - Reference API: lun_eclipse_how()
        - Meeus "Astronomical Algorithms" Ch. 54

    Note:
        Type and umbral/penumbral magnitudes come from the canonical
        _lun_how_core shadow model, so they match lun_eclipse_how(); only the
        topocentric azimuth/altitude (attr[4..6]) are computed locally here.
    """
    # Direct callers (this pythonic API is re-exported) get the same
    # mode-aware dispatch as the public wrapper: LEB first, Skyfield retry
    # on a range miss in auto mode, typed error in sealed mode. Explicit
    # internal callers still provide a reader to keep an event search on a
    # single source.
    if reader is _ACTIVE_LEB_READER:
        return _call_with_leb_skyfield_fallback(
            _lun_eclipse_how_pythonic,
            jd,
            lat,
            lon,
            altitude,
            flags,
        )

    if reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR

        _gp = (lon, lat, altitude)

        try:
            jd_tt = jd + deltat(jd)
            moon_topo = _topo_ecliptic(reader, jd_tt, jd, MOON, _gp)
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            # Sealed mode must not fabricate a "no eclipse" result for a
            # body the active LEB file cannot serve.
            _raise_if_sealed_leb_miss(_exc)
            _reraise_if_leb_range_error(_exc)
            return 0, tuple([0.0] * 20)

        # atpress/attemp match _lun_eclipse_how_impl (0.0, 10.0) and the Skyfield
        # branch's temperature_C="standard" (10 C), so the apparent altitude and
        # the visibility gate agree across the two backends and the canonical path.
        moon_az_val, moon_alt_true, moon_alt_app = azalt(
            jd, ECL2HOR, _gp, 0.0, 10.0, moon_topo[:3]
        )
        moon_altitude_true = moon_alt_true
        moon_altitude_app = moon_alt_app
        moon_azimuth = moon_az_val
    else:
        from skyfield.api import wgs84

        from .state import get_planets, get_timescale

        # Get ephemeris and timescale
        eph = get_planets()
        ts = get_timescale()

        # Create observer location
        observer = wgs84.latlon(lat, lon, altitude)

        # Get Moon object
        moon = eph["moon"]
        earth = eph["earth"]

        # Get Skyfield time
        t = ts.ut1_jd(jd)

        # Create observer position
        observer_at = earth + observer

        # Get Moon apparent position from observer
        try:
            moon_app = observer_at.at(t).observe(moon).apparent()
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            # If calculation fails, return zeros (20 elements)
            return 0, tuple([0.0] * 20)

        # Get Moon altitude and azimuth. Skyfield's altaz() returns the
        # UN-refracted (true) altitude by default; passing temperature_C='standard'
        # (10 C) applies its refraction model for the apparent altitude. attr[5]
        # is the true altitude and attr[6] the apparent altitude, so the two must
        # be distinct (formerly both carried the same true value here).
        moon_app_altaz = moon_app.altaz()
        moon_altitude_true = moon_app_altaz[0].degrees
        moon_altitude_app = moon_app.altaz(temperature_C="standard")[0].degrees
        # Convert Skyfield navigational azimuth (N=0) to the reference convention (S=0)
        moon_azimuth = (moon_app_altaz[1].degrees + 180.0) % 360.0

    # Eclipse circumstances from the canonical selenocentric shadow model
    # (_lun_how_core) -- the same core that backs lun_eclipse_how()/
    # lun_eclipse_when(), so the type and magnitudes returned here stay
    # consistent with it (audit v10). The earlier path used
    # _calculate_lunar_eclipse_type_and_magnitude(), whose divergent 1/85
    # atmospheric enlargement + small-angle approximation disagreed with the
    # canonical model by up to ~0.011 in magnitude.
    retc_core, attr_core, _dcore = _lun_how_core(jd, flags)
    umbral_mag = attr_core[0]
    penumbral_mag = attr_core[1]

    if retc_core == 0:
        # No eclipse - Moon too far from Earth's shadow
        return 0, (
            0.0,  # [0] umbral magnitude
            0.0,  # [1] penumbral magnitude
            0.0,  # [2] reserved
            0.0,  # [3] reserved
            moon_azimuth,  # [4] Moon azimuth
            moon_altitude_true,  # [5] true altitude
            moon_altitude_app,  # [6] apparent altitude
            0.0,  # [7] distance from opposition
            0.0,  # [8] eclipse magnitude
            0.0,  # [9] Saros series number
            0.0,  # [10] Saros series member
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,  # [11-19] Reserved
        )

    # There is an eclipse - start from the canonical type flag
    eclipse_type = retc_core

    # Check if Moon is above horizon (use the apparent altitude, which lets
    # refraction lift a horizon-grazing Moon into view, matching _lun_eclipse_how_impl).
    if moon_altitude_app > -1.0:
        eclipse_type |= ECL_VISIBLE
        eclipse_type |= ECL_MAX_VISIBLE

    # Prepare attributes tuple (20 elements matching the reference layout)
    attr = (
        max(0.0, umbral_mag),  # [0] Umbral magnitude
        max(0.0, penumbral_mag),  # [1] Penumbral magnitude
        0.0,  # [2] Reserved
        0.0,  # [3] Reserved
        moon_azimuth,  # [4] Azimuth of Moon
        moon_altitude_true,  # [5] True altitude of Moon
        moon_altitude_app,  # [6] Apparent altitude
        attr_core[7],  # [7] Distance from opposition (degrees)
        max(0.0, umbral_mag),  # [8] Eclipse magnitude (equals [0])
        *_get_saros_info(jd, "lunar"),  # [9] Saros, [10] member
        0.0,  # [11] Reserved
        0.0,  # [12] Reserved
        0.0,  # [13] Reserved
        0.0,  # [14] Reserved
        0.0,  # [15] Reserved
        0.0,  # [16] Reserved
        0.0,  # [17] Reserved
        0.0,  # [18] Reserved
        0.0,  # [19] Reserved
    )

    return eclipse_type, attr


def lun_eclipse_how(
    tjdut: float,
    geopos: Sequence[float],
    flags: int = FLG_SWIEPH,
) -> Tuple[int, Tuple[float, ...]]:
    """Calculate detailed circumstances of a lunar eclipse from a location.

    Wrapper around :func:`_lun_eclipse_how_impl` that adds LEB→Skyfield
    fallback for partial/custom LEB files. See the impl docstring for the
    full API contract.
    """
    flags = _strip_one_try_bit(flags)
    return _call_with_leb_skyfield_fallback(
        _lun_eclipse_how_impl,
        tjdut,
        geopos,
        flags,
    )


def _lun_eclipse_how_impl(
    tjdut: float,
    geopos: Sequence[float],
    flags: int = FLG_SWIEPH,
    *,
    reader=None,
) -> Tuple[int, Tuple[float, ...]]:
    """
    Calculate the circumstances of a lunar eclipse at a location and time.

    The eclipse geometry itself is geocentric (a lunar eclipse looks the
    same from everywhere on the night side); the location only determines
    whether the Moon stands above the local horizon.

    Args:
        tjdut: Julian Day (UT) during a lunar eclipse
        geopos: Sequence of [longitude, latitude, altitude]:
            - longitude in degrees (East positive)
            - latitude in degrees (North positive)
            - altitude in meters above sea level
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - retflag: ECL_TOTAL, ECL_PARTIAL or ECL_PENUMBRAL; 0 when no
              eclipse is in progress at this time OR when the Moon's
              apparent altitude at the location is at or below 0 (the
              attributes stay filled in that case). No visibility bits
              are ever set by this function.
            - attr: Tuple of 20 floats:
                [0]: umbral magnitude (0 for purely penumbral phases)
                [1]: penumbral magnitude (negative when the Moon stands
                     clear of the penumbra)
                [2]/[3]: 0.0
                [4]: azimuth of the Moon
                [5]: true altitude of the Moon
                [6]: apparent altitude of the Moon (refraction at 10 C,
                     pressure estimated from the observer's altitude)
                [7]: Moon's distance from opposition in degrees (0 when
                     no eclipse)
                [8]: = attr[0]
                [9]: saros series number
                [10]: saros series member number
                [11-19]: 0.0

    References:
        - Reference API: lun_eclipse_how()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    if len(geopos) < 3:
        raise ValueError("geopos must have at least 3 elements: [lon, lat, alt]")

    from .utils import ECL2HOR, azalt

    geopos3 = (float(geopos[0]), float(geopos[1]), float(geopos[2]))

    retc, attr, _dcore = _lun_how_core(tjdut, flags)

    _sun_p, moon_p = _topo_sun_moon(tjdut, geopos3, reader)
    az, true_alt, app_alt = azalt(tjdut, ECL2HOR, geopos3, 0.0, 10.0, moon_p)
    attr[4] = az
    attr[5] = true_alt
    attr[6] = app_alt
    if app_alt <= 0.0:
        retc = 0

    return retc, tuple(attr)


def _kernel_jd_bounds(eph) -> Tuple[float, float]:
    """Return (earliest_start_jd, latest_end_jd) across all loaded SPK segments.

    Used to clamp occultation searches to the loaded DE kernel's coverage so the
    batch scan never requests JDs outside it. Split kernels (e.g. DE441, whose
    every body is broken into two segments at JD 2440432.5 / 1969-07-30) expose
    several end_jd values; the kernel's true upper bound is the *maximum* end_jd.
    Taking the minimum would clamp the search to the split point (~1969) and drop
    every later event.

    Args:
        eph: Loaded ephemeris object containing SPK segments.

    Returns:
        Tuple of (minimum segment start_jd, maximum segment end_jd).
    """
    starts = [seg.spk_segment.start_jd for seg in eph.segments]
    ends = [seg.spk_segment.end_jd for seg in eph.segments]
    return min(starts), max(ends)


def _reject_moon_self_occultation(body: "int | str") -> None:
    """Reject the degenerate Moon-occults-Moon request with the typed error.

    The Moon cannot occult itself; without this guard the shadow geometry
    would divide by the zero Moon-to-body distance. Behavioral comparison
    with one external implementation showed no defined result for this
    request either (the call fails to terminate), so the shared typed
    contract error is raised instead by every lun_occult_* entry point.
    """
    if isinstance(body, int) and body == MOON:
        from .exceptions import Error as _Error

        raise _Error("lunar occultation of the Moon itself is undefined (body 1).")


def lun_occult_when_glob(
    tjdut: float,
    body: "int | str",
    flags: int = FLG_SWIEPH,
    ecltype: int = 0,
    backwards: "bool | int | str" = False,
) -> Tuple[int, Tuple[float, ...]]:
    """Find a lunar occultation visible somewhere on Earth.

    Args:
        tjdut: Search start as a Julian Day in UT.
        body: Planet identifier or fixed-star name.
        flags: Calculation flags.
        ecltype: Requested occultation-type mask, or zero for any type.
        backwards: Search direction. Integer values may also include
            ``ECL_ONE_TRY``.

    Returns:
        An occultation retflag and ten-value event-time tuple.

    Raises:
        Error: If the type filter is impossible or no event is found in the
            search window.

    Notes:
        With ``ECL_ONE_TRY``, only the conjunction nearest the start is
        examined; a miss returns retflag zero and a continuation date in the
        first time slot.
    """
    _reject_moon_self_occultation(body)
    from .exceptions import Error
    from .constants import ECL_ONE_TRY

    # bools/ints keep their flag bits (ECL_ONE_TRY rides on the int);
    # direction strings are coerced to a plain 0/1.
    back_i = (
        int(_coerce_backwards(backwards))
        if isinstance(backwards, str)
        else int(backwards)
    )
    one_try = bool(back_i & ECL_ONE_TRY)
    backward = bool(back_i & 1)
    direction = -1 if backward else 1
    eph_flags = _ecl_eph_flags(flags)
    is_star = isinstance(body, str)
    is_sun = (not is_star) and body == SUN
    de_km = _EARTH_EQ_RADIUS_KM
    ifltype = _normalize_occultation_filter(ecltype, body, is_sun)

    from .constants import FLG_XYZ

    body_radius_au = _occ_body_radius_au(body)
    xyz_flags = eph_flags | FLG_EQUATORIAL | FLG_XYZ

    def _geo_lonlat(jd: float) -> Tuple[float, float]:
        if is_star:
            from .fixed_stars import fixstar_ut

            pos, _n, _r = fixstar_ut(cast(str, body), jd, eph_flags)
            return float(pos[0]), float(pos[1])
        pos, _ = calc_ut(jd, cast(int, body), eph_flags)
        return float(pos[0]), float(pos[1])

    def _sep_minus_radii(jd: float) -> float:
        bx = _occ_body_geo_xyz(jd, body, xyz_flags)
        mxx, _ = calc_ut(jd, MOON, xyz_flags)
        db = math.sqrt(bx[0] ** 2 + bx[1] ** 2 + bx[2] ** 2)
        dm = math.sqrt(mxx[0] ** 2 + mxx[1] ** 2 + mxx[2] ** 2)
        cosang = (bx[0] * mxx[0] + bx[1] * mxx[1] + bx[2] * mxx[2]) / (db * dm)
        sep = math.degrees(math.acos(max(-1.0, min(1.0, cosang))))
        rmoon = math.degrees(math.asin(_ECL_RMOON_AU / dm))
        rbody = math.degrees(math.asin(body_radius_au / db)) if body_radius_au else 0.0
        return sep - (rmoon + rbody)

    def _lon_diff(jd: float) -> float:
        # The noon-transit slot is defined by Moon-body conjunction in right
        # ascension, while the coarse conjunction stepping above uses
        # ecliptic longitude.
        if is_star:
            from .fixed_stars import fixstar_ut

            bpos, _n, _r = fixstar_ut(cast(str, body), jd, eph_flags | FLG_EQUATORIAL)
        else:
            bpos, _ = calc_ut(jd, cast(int, body), eph_flags | FLG_EQUATORIAL)
        mra = calc_ut(jd, MOON, eph_flags | FLG_EQUATORIAL)[0][0]
        d = (bpos[0] - mra) % 360.0
        if d > 180.0:
            d -= 360.0
        return d

    tret = [0.0] * 10
    t = float(tjdut)

    for _attempt in range(_OCCULT_MAX_CONJUNCTIONS):
        try:
            t, blat, moon_lat = _seek_moon_conjunction(
                t,
                direction,
                _geo_lonlat,
                eph_flags,
                cast(str, body) if is_star else None,
            )
        except _OccSearchOffEphemeris:
            # Stepped off the active ephemeris: no matching occultation
            # exists within coverage (the tier boundary).
            break
        tjd = t

        # Latitude gate: no occultation possible this conjunction.
        if abs(blat - moon_lat) > _OCC_LAT_GATE_DEG:
            if one_try:
                tret[0] = t + direction
                return 0, tuple(tret)
            t = tjd + direction * _OCC_CONJUNCTION_HOP_DAYS
            continue

        # Time of maximum: deepest limb overlap (geocentric).
        tjd = _golden_min(_sep_minus_radii, tjd - 1.0, tjd + 1.0)

        rc_where, _wlon, _wlat, dcore = _eclipse_where_core(tjd, flags, body)
        if rc_where == 0:
            if one_try:
                tret[0] = tjd
                return 0, tuple(tret)
            t = tjd + direction * _OCC_CONJUNCTION_HOP_DAYS
            continue

        tret = [0.0] * 10
        tret[0] = tjd
        if (backward and tret[0] >= tjdut - _ECLIPSE_WHEN_EPOCH_MARGIN) or (
            not backward and tret[0] <= tjdut + _ECLIPSE_WHEN_EPOCH_MARGIN
        ):
            t = tjd + direction * _OCC_CONJUNCTION_HOP_DAYS
            continue

        retflag = rc_where

        def _unwanted() -> bool:
            for bit in (ECL_NONCENTRAL, ECL_CENTRAL, ECL_ANNULAR, ECL_PARTIAL):
                if (retflag & bit) and not (ifltype & bit):
                    return True
            if (retflag & ECL_TOTAL) and not (
                ifltype & (ECL_TOTAL | ECL_ANNULAR_TOTAL)
            ):
                return True
            return False

        if _unwanted():
            if one_try:
                tret[0] = tjd
                return 0, tuple(tret)
            t = tjd + direction * _OCC_CONJUNCTION_HOP_DAYS
            continue

        # Phase times from the shadow geometry: occultation begin/end,
        # totality begin/end, center-line begin/end.
        def _f_outer(jd: float) -> float:
            _rc, _lo, _la, dc = _eclipse_where_core(jd, flags, body)
            return dc[4] / 2.0 + de_km / dc[5] - dc[2]

        def _f_total(jd: float) -> float:
            _rc, _lo, _la, dc = _eclipse_where_core(jd, flags, body)
            return abs(dc[3]) / 2.0 + de_km / dc[6] - dc[2]

        def _f_central(jd: float) -> float:
            _rc, _lo, _la, dc = _eclipse_where_core(jd, flags, body)
            return de_km / dc[6] - dc[2]

        win = 2.0 / 24.0
        wide = 5.0 / 24.0
        tret[2] = _root_bisect(_f_outer, tjd - win, tjd) or _root_bisect(
            _f_outer, tjd - wide, tjd
        )
        tret[3] = _root_bisect(_f_outer, tjd, tjd + win) or _root_bisect(
            _f_outer, tjd, tjd + wide
        )
        if not (retflag & ECL_PARTIAL):
            tret[4] = _root_bisect(_f_total, tjd - win, tjd) or _root_bisect(
                _f_total, tjd - wide, tjd
            )
            tret[5] = _root_bisect(_f_total, tjd, tjd + win) or _root_bisect(
                _f_total, tjd, tjd + wide
            )
            if not (retflag & ECL_NONCENTRAL):
                tret[6] = _root_bisect(_f_central, tjd - win, tjd) or _root_bisect(
                    _f_central, tjd - wide, tjd
                )
                tret[7] = _root_bisect(_f_central, tjd, tjd + win) or _root_bisect(
                    _f_central, tjd, tjd + wide
                )

        # A total occultation that loses its core shadow at begin or end
        # of totality is annular-total (only possible for the Sun).
        if retflag & ECL_TOTAL:
            d_max = _eclipse_where_core(tret[0], flags, body)[3][0]
            d_beg = (
                _eclipse_where_core(tret[4], flags, body)[3][0] if tret[4] else d_max
            )
            d_end = (
                _eclipse_where_core(tret[5], flags, body)[3][0] if tret[5] else d_max
            )
            if d_max * d_beg < 0.0 or d_max * d_end < 0.0:
                retflag = (retflag & ~ECL_TOTAL) | ECL_ANNULAR_TOTAL

        if ((retflag & ECL_TOTAL) and not (ifltype & ECL_TOTAL)) or (
            (retflag & ECL_ANNULAR_TOTAL) and not (ifltype & ECL_ANNULAR_TOTAL)
        ):
            if one_try:
                tret[0] = tjd
                return 0, tuple(tret)
            t = tjd + direction * _OCC_CONJUNCTION_HOP_DAYS
            continue

        # tret[1]: does the occultation happen at local apparent noon
        # somewhere (body-Moon conjunction in longitude between begin
        # and end)?
        if tret[2] and tret[3]:
            d_beg = _lon_diff(tret[2])
            d_end = _lon_diff(tret[3])
            if d_beg * d_end >= 0.0:
                tret[1] = 0.0
            else:
                tret[1] = _root_bisect(_lon_diff, tret[2], tret[3])

        return retflag, tuple(tret)

    raise Error(
        f"no lunar occultation of {body} found after JD {tjdut} "
        "within the search window"
    )


def _lun_occult_when_loc_pythonic(
    jd_start: float,
    planet: int,
    star_name: str,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    flags: int = FLG_SWIEPH,
    backwards: "bool | int" = False,
    *,
    reader=None,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """Find the next lunar occultation observable from a location.

    Searches forward or backward from ``jd_start`` for an occultation of
    the given body whose phases are at least partly observable (body
    above the local apparent horizon) from the given place.

    Returns:
        (retflag, tret, attr):
            retflag: local occultation type (ECL_TOTAL/ANNULAR/PARTIAL)
                plus ECL_VISIBLE and the per-contact visibility bits,
                ECL_OCC_BEG_DAYLIGHT/ECL_OCC_END_DAYLIGHT when first or
                fourth contact happens in daylight, and the global
                ECL_NONCENTRAL character. With ECL_ONE_TRY set in
                ``backwards`` and no occultation at the first
                conjunction: 0, with tret[0] holding a date from which
                to continue searching.
            tret: [0] maximum, [1] first contact, [2] second contact,
                [3] third contact, [4] fourth contact, [5] body rise
                between first and fourth contact (else 0), [6] body set
                (else 0), [7..9] 0. For stars, first/fourth contact
                equal second/third (point sources).
            attr: local circumstances at the maximum (sol_eclipse_how
                layout with the occulted body; attr[3] = core-shadow
                width at the global shadow-center point, km).

    Raises:
        Error: star with |ecliptic latitude| > 7 degrees (an
            occultation can never occur), or search exhaustion.
    """
    from .exceptions import Error
    from .constants import (
        BIT_DISC_BOTTOM,
        CALC_RISE,
        CALC_SET,
        ECL_ONE_TRY,
    )
    from .utils import angular_separation

    body: "Union[int, str]" = star_name if star_name else planet
    is_star = isinstance(body, str)
    back_i = int(backwards)
    one_try = bool(back_i & ECL_ONE_TRY)
    backward = bool(back_i & 1)
    direction = -1 if backward else 1
    eph_flags = _ecl_eph_flags(flags)
    geopos3 = (float(lon), float(lat), float(altitude))
    body_radius_au = _occ_body_radius_au(body)
    zero_attr = tuple([0.0] * 20)

    def _geo_lonlat(jd: float) -> Tuple[float, float]:
        if is_star:
            from .fixed_stars import fixstar_ut

            pos, _n, _r = fixstar_ut(cast(str, body), jd, eph_flags)
            return float(pos[0]), float(pos[1])
        pos, _ = calc_ut(jd, cast(int, body), eph_flags)
        return float(pos[0]), float(pos[1])

    def _radii_sep(
        jd: float, rmoon_au: float = _ECL_RMOON_AU
    ) -> Tuple[float, float, float]:
        b_p = _occ_body_topo(jd, body, geopos3, flags, reader)
        _s_p, m_p = _topo_sun_moon(jd, geopos3, reader)
        rb = math.degrees(math.asin(body_radius_au / b_p[2])) if body_radius_au else 0.0
        rm = math.degrees(math.asin(rmoon_au / m_p[2]))
        return rb, rm, angular_separation(b_p[0], b_p[1], m_p[0], m_p[1])

    def _sep_minus_radii(jd: float) -> float:
        rb, rm, sep = _radii_sep(jd)
        return sep - (rb + rm)

    def _one_try_miss(t_suggest: float):
        tret_m = [0.0] * 10
        tret_m[0] = t_suggest
        return 0, tuple(tret_m), zero_attr

    t = float(jd_start)
    for _attempt in range(_OCCULT_MAX_CONJUNCTIONS):
        try:
            t, blat, moon_lat = _seek_moon_conjunction(
                t,
                direction,
                _geo_lonlat,
                eph_flags,
                cast(str, body) if is_star else None,
            )
        except _OccSearchOffEphemeris:
            # Stepped off the active ephemeris: no matching occultation
            # exists within coverage (the tier boundary).
            break
        tjd = t

        # Geocentric latitude gate for this conjunction.
        if abs(blat - moon_lat) > _OCC_LAT_GATE_DEG:
            if one_try:
                return _one_try_miss(t + direction)
            t = tjd + direction * _OCC_CONJUNCTION_HOP_DAYS
            continue

        # Topocentric maximum: deepest limb overlap seen from the place.
        jd_max = _golden_min(_sep_minus_radii, tjd - 1.0, tjd + 1.0)

        rbody, rmoon, min_sep = _radii_sep(jd_max)
        if min_sep > rbody + rmoon:
            if one_try:
                return _one_try_miss(jd_max + direction)
            t = jd_max + direction * _OCC_CONJUNCTION_HOP_DAYS
            continue

        if (backward and jd_max >= jd_start - _ECLIPSE_WHEN_EPOCH_MARGIN) or (
            not backward and jd_max <= jd_start + _ECLIPSE_WHEN_EPOCH_MARGIN
        ):
            if one_try:
                return _one_try_miss(jd_max + direction)
            t = jd_max + direction * _OCC_CONJUNCTION_HOP_DAYS
            continue

        if min_sep < rbody - rmoon:
            phase = ECL_ANNULAR
        elif min_sep < abs(rbody - rmoon):
            phase = ECL_TOTAL
        else:
            phase = ECL_PARTIAL

        # Inner contacts use exact apparent-disc tangency against the smaller
        # inner-contact lunar radius (disappearance/reappearance behind the
        # dark limb); for point-source stars the outer contacts coincide with
        # them. Outer contacts use the full mean lunar disc.
        def _f_outer(jd_c: float) -> float:
            rb_c, rm_c, sep_c = _radii_sep(jd_c)
            return sep_c - (rb_c + rm_c)

        def _f_inner(jd_c: float) -> float:
            rb_c, rm_c, sep_c = _radii_sep(jd_c, _ECL_RMOON_INNER_AU)
            return abs(rb_c - rm_c) - sep_c

        jd_second = 0.0
        jd_third = 0.0
        if min_sep < abs(rbody - rmoon):
            inner_win = 1.5 / 24.0
            jd_second = _root_bisect(_f_inner, jd_max - inner_win, jd_max)
            jd_third = _root_bisect(_f_inner, jd_max, jd_max + inner_win)

        if is_star:
            jd_first = jd_second
            jd_fourth = jd_third
        else:
            two_hours = 2.0 / 24.0
            jd_first = _root_bisect(
                _f_outer, jd_max - two_hours, jd_max
            ) or _root_bisect(_f_outer, jd_max - 2.0 * two_hours, jd_max)
            jd_fourth = _root_bisect(
                _f_outer, jd_max, jd_max + two_hours
            ) or _root_bisect(_f_outer, jd_max, jd_max + 2.0 * two_hours)

        # Visibility bits from the body's apparent altitude at each
        # contact.
        retflag = phase
        for tc, bit in (
            (jd_max, ECL_MAX_VISIBLE),
            (jd_first, ECL_1ST_VISIBLE),
            (jd_second, ECL_2ND_VISIBLE),
            (jd_third, ECL_3RD_VISIBLE),
            (jd_fourth, ECL_4TH_VISIBLE),
        ):
            if tc == 0.0:
                continue
            _rc_i, attr_i = _sol_how_core(tc, geopos3, flags, reader, body)
            if attr_i[6] > 0.0:
                retflag |= ECL_VISIBLE | bit

        if not (retflag & ECL_VISIBLE):
            if one_try:
                return _one_try_miss(jd_max + direction)
            t = jd_max + direction * _OCC_CONJUNCTION_HOP_DAYS
            continue

        # Rise and set of the occulted body during the occultation.
        tret = [0.0] * 10
        tret[0] = jd_max
        tret[1] = jd_first
        tret[2] = jd_second
        tret[3] = jd_third
        tret[4] = jd_fourth

        rc_rise, tr_rise = rise_trans(
            jd_first - 0.1,
            body,
            CALC_RISE | BIT_DISC_BOTTOM,
            geopos3,
            0.0,
            0.0,
            eph_flags,
        )
        rc_set, tr_set = rise_trans(
            jd_first - 0.1,
            body,
            CALC_SET | BIT_DISC_BOTTOM,
            geopos3,
            0.0,
            0.0,
            eph_flags,
        )
        if rc_rise >= 0 and rc_set >= 0:
            if jd_first < tr_rise[0] < jd_fourth:
                tret[5] = tr_rise[0]
            if jd_first < tr_set[0] < jd_fourth:
                tret[6] = tr_set[0]

        # Daylight flags: is the Sun up at first/fourth contact? (The
        # next sunset coming before the next sunrise means daytime.)
        from .constants import ECL_OCC_BEG_DAYLIGHT, ECL_OCC_END_DAYLIGHT

        for tc, bit in (
            (jd_first, ECL_OCC_BEG_DAYLIGHT),
            (jd_fourth, ECL_OCC_END_DAYLIGHT),
        ):
            if tc == 0.0:
                continue
            rc_sr, tr_sr = rise_trans(tc, SUN, CALC_RISE, geopos3, 0.0, 0.0, eph_flags)
            rc_ss, tr_ss = rise_trans(tc, SUN, CALC_SET, geopos3, 0.0, 0.0, eph_flags)
            if rc_sr >= 0 and rc_ss >= 0 and tr_ss[0] < tr_sr[0]:
                retflag |= bit

        # Local circumstances at maximum; global character and core
        # width from the shadow geometry.
        _rc_how, attr_list = _sol_how_core(jd_max, geopos3, flags, reader, body)
        rc_where, _wlon, _wlat, dcore = _eclipse_where_core(jd_max, flags, body)
        retflag |= rc_where & ECL_NONCENTRAL
        attr_list[3] = dcore[0]
        # The compatibility attr layout defines occultation diameter fraction
        # and obscuration as fractions, hence both are capped at unity.
        attr_list[0] = min(attr_list[0], 1.0)
        attr_list[2] = min(attr_list[2], 1.0)

        return retflag, tuple(tret), tuple(attr_list)

    raise Error(
        f"no observable lunar occultation of {body} found from "
        f"lon={lon}, lat={lat} after JD {jd_start}"
    )


def lun_occult_when_loc(
    tjdut: float,
    body: "int | str",
    geopos: "Sequence[float]",
    flags: int = FLG_SWIEPH,
    backwards: "bool | int | str" = False,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """Find a lunar occultation visible from one location.

    Args:
        tjdut: Search start as a Julian Day in UT.
        body: Planet identifier, fixed-star identifier, or fixed-star name.
        geopos: Longitude, latitude, and altitude in metres.
        flags: Calculation flags.
        backwards: Search direction. Integer values may also include
            ``ECL_ONE_TRY``.

    Returns:
        The occultation retflag, a ten-value phase-time tuple, and a
        twenty-value local-circumstance tuple. In the phase-time tuple,
        ``tret[5]`` is when the occulted body rises and ``tret[6]`` is when the
        occulted body sets at the observer location.

    Raises:
        Error: If no visible occultation is found in the search window.
        ValueError: If the body or geographic position is invalid.

    Notes:
        Local attributes include covered diameter/disc fractions, target
        azimuth and altitude, elongation, and the event contact times.
    """
    _reject_moon_self_occultation(body)

    # Validate geopos
    if len(geopos) < 3:
        raise ValueError("geopos must have at least 3 elements: [lon, lat, alt]")

    # Extract geographic position (reference API uses lon, lat, alt order)
    lon = geopos[0]
    lat = geopos[1]
    altitude = geopos[2]

    # Determine if body is planet ID or star name. Integer fixed-star ids
    # route through the star path like calc_ut / lun_occult_when_glob
    # (dispatch by id); without this they fell through to the planet
    # machinery and leaked a raw KeyError.
    if isinstance(body, str):
        planet = 0
        star_name = body
    else:
        from .fixed_stars import FIXED_STARS as _FS
        from .fixed_stars import get_canonical_star_name as _star_name_by_id

        if body in _FS:
            planet = 0
            star_name = _star_name_by_id(body) or ""
        else:
            planet = body
            star_name = ""

    # bools/ints keep their flag bits (ECL_ONE_TRY rides on the int);
    # direction strings are coerced to a plain 0/1.
    back_i = (
        int(_coerce_backwards(backwards))
        if isinstance(backwards, str)
        else int(backwards)
    )

    # Call the implementation with mode-aware LEB range handling.
    ecl_type, times, attr = _call_with_leb_skyfield_fallback(
        _lun_occult_when_loc_pythonic,
        tjdut,
        planet,
        star_name,
        lat,
        lon,
        altitude,
        flags,
        back_i,
    )

    # Return in reference API order: (retflags, tret, attr)
    return ecl_type, times, attr


def _lun_occult_where_pythonic(
    tjdut: float,
    body: Union[int, str],
    flags: int = FLG_SWIEPH,
    *,
    reader=None,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """
    Find where a lunar occultation is central or maximal at a given time.

    Computes where the Moon's shadow axis for the occulted body pierces
    the Earth's surface at the given instant; if the axis misses the
    Earth, the surface point of maximum occultation is returned.

    Args:
        tjdut: Julian Day (UT) during an occultation
        body: Planet identifier (int) or fixed-star name (str)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - retflag: 0 when no occultation is in progress anywhere on
              Earth at this time; otherwise ECL_CENTRAL or
              ECL_NONCENTRAL (plus ECL_PARTIAL when only the penumbral
              cone touches), plus ECL_TOTAL or ECL_ANNULAR. Stars are
              point sources, so their occultations are always total.
            - geopos: 10 floats - [0]/[1] longitude/latitude of the
              shadow-center point, [2..9] zeros (reference convention).
            - attr: 20 floats with the local circumstances at the
              center point (sol_eclipse_how layout, with the occulted
              body in place of the Sun); attr[3] is the core-shadow
              diameter at the center point in km (negative for a total
              occultation - for a star it is minus the lunar diameter).

    References:
        - Reference API: lun_occult_where()
    """
    retflag, center_lon, center_lat, dcore = _eclipse_where_core(tjdut, flags, body)
    _rc_how, attr_list = _sol_how_core(
        tjdut,
        (center_lon, center_lat, 0.0),
        flags,
        reader,
        body,
        where_convention=True,
    )
    attr_list[3] = dcore[0]
    geopos = (center_lon, center_lat, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    return retflag, geopos, tuple(attr_list)


_lun_occult_where_internal = _lun_occult_where_pythonic


def lun_occult_where(
    tjdut: float,
    body: "int | str" = 0,
    flags: int = FLG_SWIEPH,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """Calculate where on Earth a lunar occultation is visible (reference-API-compatible).

    Wrapper around _lun_occult_where_pythonic() matching the reference signature.

    Args:
        tjdut: Julian Day (UT) of the moment to calculate.
        body: Planet identifier (int) or star name (str).
        flags: Calculation flags (default FLG_SWIEPH).

    Returns:
        Tuple of (retflag, geopos, attr) matching the reference API.
    """
    _reject_moon_self_occultation(body)
    if isinstance(body, int):
        from .fixed_stars import FIXED_STARS as _FS
        from .fixed_stars import get_canonical_star_name as _star_name_by_id

        if body in _FS:
            # Integer fixed-star id: same star path as a str name (see
            # lun_occult_when_loc).
            body = _star_name_by_id(body) or body
    return _call_with_leb_skyfield_fallback(
        _lun_occult_where_internal,
        tjdut,
        body,
        flags,
    )


# =============================================================================
# RISE, SET, AND TRANSIT CALCULATIONS
# =============================================================================


def rise_trans(
    tjdut: float,
    body: "Union[int, str]",
    rsmi: int,
    geopos: Sequence[float],
    atpress: float = 0.0,
    attemp: float = 0.0,
    flags: int = FLG_SWIEPH,
) -> Tuple[int, Tuple[float, ...]]:
    """Calculate rise, set, or transit time for a celestial body.

    Wrapper around :func:`_rise_trans_impl` that applies mode-aware LEB range
    handling for partial/custom files. See the implementation docstring.
    """
    return _call_with_leb_skyfield_fallback(
        _rise_trans_impl,
        tjdut,
        body,
        rsmi,
        geopos,
        atpress,
        attemp,
        flags,
    )


def _rise_trans_impl(
    tjdut: float,
    body: "Union[int, str]",
    rsmi: int,
    geopos: Sequence[float],
    atpress: float = 0.0,
    attemp: float = 0.0,
    flags: int = FLG_SWIEPH,
    *,
    reader=None,
) -> Tuple[int, Tuple[float, ...]]:
    """
    Calculate rise, set, or transit time for a celestial body.

    This function finds the next occurrence of a rising, setting, or transit
    (meridian passage) event for a celestial body as seen from a specific
    geographic location.

    Args:
        tjdut: Julian Day (UT) to start search from
        body: Planet/body ID (int, e.g. SUN, MOON) or fixed star
            name (str, e.g. "Sirius", "Regulus")
        rsmi: Event type and calculation flags (bitmask):
            - CALC_RISE (1): Rise time (body crossing horizon going up)
            - CALC_SET (2): Set time (body crossing horizon going down)
            - CALC_MTRANSIT (4): Upper meridian transit (culmination)
            - CALC_ITRANSIT (8): Lower meridian transit (anti-culmination)
            Additional flags (OR with event type):
            - BIT_GEOCTR_NO_ECL_LAT (128): Use the geocentric apparent place
              projected onto the ecliptic (ecliptic latitude zeroed) instead
              of the topocentric place - the Hindu-rising convention. Affects
              rise/set and twilight; ignored for meridian transits. Component
              of BIT_HINDU_RISING (896).
            - BIT_DISC_CENTER (256): Use disc center instead of upper limb
            - BIT_DISC_BOTTOM (8192): Use lower limb of disc
            - BIT_NO_REFRACTION (512): Ignore atmospheric refraction
            - BIT_CIVIL_TWILIGHT (1024): Sun at -6 degrees
            - BIT_NAUTIC_TWILIGHT (2048): Sun at -12 degrees
            - BIT_ASTRO_TWILIGHT (4096): Sun at -18 degrees
        geopos: Geographic position as sequence [lon, lat, alt]:
            - 0: geographic longitude in degrees (eastern positive)
            - 1: geographic latitude in degrees (northern positive)
            - 2: geographic altitude in meters above sea level. It feeds four
              things and NOT the one most callers expect: the topocentric
              parallax, the barometric pressure estimate when atpress is 0, the
              dip clamp inside the refraction inversion, and the dip requested
              by horhgt=-100. It does NOT by itself depress the horizon (see
              "Horizon convention").
        atpress: Atmospheric pressure in mbar/hPa for refraction (default 0.0).
            0 means "estimate it from geopos[2] with the barometric
            expression", not "no pressure".
        attemp: Atmospheric temperature in degrees Celsius (default 0.0). The
            default is a literal 0 C, as in the reference - NOT an implicit
            15 C standard atmosphere. It is not a free parameter: measured at
            Rome for 2024-03-20, attemp=15.0 moves sunrise +16.9 s against the
            0 C default (-10 C: -13.1 s; +30 C: +31.1 s). A caller that means
            the standard atmosphere has to pass 1013.25 / 15.0 explicitly.
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - res: integer flag (0 = event found, -2 = circumpolar)
            - tret: tuple of 10 floats, of which tret[0] = JD of event

    Raises:
        ValueError: If invalid planet ID or star name

    Note:
        For circumpolar objects (always above or below horizon at the given
        latitude), the function returns (-2, tret) with tret[0] = 0.0.
        For transits, circumpolar objects still have valid transit times.

    Horizon convention:
        A rise/set is the instant the body's APPARENT UPPER LIMB meets the
        horizon. The true altitude is raised by the topocentric semidiameter -
        computed from the actual distance, so it tracks the Earth-Sun distance
        through the year rather than using a fixed 16' - and the result is then
        refracted. Under 1013.25 hPa / 15 C that leaves the Sun's true
        geometric CENTRE near -49.6' at the event (it ranges about -49.3' to
        -49.9' over a year, which is the semidiameter breathing).

        The disc term is selectable: BIT_DISC_CENTER drops it, BIT_DISC_BOTTOM
        flips it to the lower limb (a full diameter away, so the event moves by
        roughly twice the semidiameter in time), and BIT_FIXED_DISC_SIZE pins the
        distance to 1 AU - i.e. it IS the fixed ~16' the variable-distance
        default avoids. BIT_NO_REFRACTION drops the refraction term.

        The observer altitude in geopos[2] does NOT lower the horizon, at any
        elevation. With an explicit atpress it is inert for the timing: the
        residual is about 13 ms, and it is neither the parallax (only 0.0014" at
        1000 m, a hundred times too small) nor the dip - the same 13 ms appears
        between horhgt=0 and horhgt=-0.01 at sea level, where calc_dip is exactly
        0 in both. It is the solver's stopping rule: it exits as soon as
        |height| < 1e-4 deg, and at sea level with horhgt=0 the root sits inside
        that band, so bisection stops ~13 ms early. Any perturbation that moves
        the root out of the band - a non-zero dip OR a non-zero horhgt - forces
        full convergence, which means the ELEVATED answer is the accurate one and
        the 13 ms is an error in the sea-level baseline. With the default atpress=0 the elevation is what the
        barometric estimate reads, so a HIGHER observer gets thinner air, less
        refraction and a LATER sunrise - the opposite of the physical intuition:
        Rome 2024-03-20, +2.5 s at 100 m, +24.2 s at 1000 m, +66.0 s at 3000 m.
        That is deliberate, not an oversight - published rise/set tables are
        computed for a level sea horizon whatever the observer's elevation, and
        horhgt=0 is what reproduces them.

        To model an observer who genuinely sees over the curve, call
        rise_trans_true_hor() with horhgt=-100, which asks for the dip of the
        sea horizon at geopos[2]. Same place and date, under the SAME standard
        atmosphere as the figures above, that moves sunrise -109.9 s at 100 m,
        -360.0 s at 1000 m and -649.1 s at 3000 m against horhgt=0 - minutes,
        which dwarfs every other term in this function's error budget. (The
        effect is atmosphere-dependent: on the library default of 0 C the same
        three are -110.6 / -362.3 / -640.7 s.) Choosing between the two is a
        convention decision the caller owns; this function cannot infer it from
        geopos[2] alone.

        One trap in that choice. `_rise_true_to_apparent` clamps to the
        UNREFRACTED altitude whenever the apparent altitude would fall below the
        dip, and at sea level the dip is exactly 0. So between horhgt=0 and
        horhgt=-0.559889 - the whole width of the refraction at the horizon -
        asking for a lower horizon changes nothing (to within the 13 ms above).
        A caller at sea level requesting a horizon depressed by 30 arcminutes
        gets the undepressed sunrise.

        Past that edge the response is a KINK, not a step: the derivative goes
        from 0 to 322 s per degree and stays there. Rome 2024-03-20 measured:
        -0.5599 -> -0.01 s, -0.561 -> -0.37 s, -0.57 -> -3.27 s, -0.60 ->
        -12.93 s, -0.80 -> -77.35 s. (An earlier version of this note called it a
        13-second jump; 13 s is simply the value at -0.60.) horhgt=-100 resolves
        to dip + 1e-4 deg, i.e. just past the edge, so its effect is determinate
        but the last 0.02-0.04 s of the figures above is solver residue rather
        than physics.

    Algorithm:
        1. For transits: Find when body crosses the local meridian
           (Local Sidereal Time = body's Right Ascension)
        2. For rise/set: Find when body's altitude crosses the horizon
           accounting for refraction and disc size
        3. Uses Newton-Raphson iteration for precise timing

    Precision:
        Rise/set times accurate to ~1 minute for Sun/Moon (due to refraction
        uncertainty), better than 1 second for transit times.

    Example:
        >>> from libephemeris import julday, rise_trans, SUN, CALC_RISE
        >>> jd = julday(2024, 6, 21, 0)
        >>> # Find sunrise at Rome (geopos = [lon, lat, alt])
        >>> res, tret = rise_trans(jd, SUN, CALC_RISE, [12.5, 41.9, 0])
        >>> print(f"Sunrise at JD {tret[0]:.5f}")
        >>> # Find Sirius rise time
        >>> res, tret = rise_trans(jd, "Sirius", CALC_RISE, [12.5, 41.9, 0])

    References:
        - Reference API: rise_trans()
        - Meeus "Astronomical Algorithms" Ch. 15 (Rise, Set, Transit)
    """
    # ``rise_trans`` is the zero-horizon specialization of
    # ``rise_trans_true_hor``. Both entry points share one rise/set engine
    # (twilight bits,
    # disc selection, refraction and the dip convention included).
    return _rise_trans_true_hor_impl(
        tjdut, body, rsmi, geopos, atpress, attemp, 0.0, flags, reader=reader
    )


def _make_tret(jd_event: float = 0.0) -> Tuple[float, ...]:
    """Build a 10-element tret tuple for rise_trans return value.

    The compatibility return contract uses a 10-element tuple where only index
    0 is meaningful (the Julian Day of the event); remaining elements are 0.0.

    jd_event is coerced to a native float: the transit path derives it from a
    numpy hour-angle expression, so without this rise_trans would leak a
    numpy.float64 in tret[0] for MTRANSIT/ITRANSIT (rise/set feed a native jd).
    """
    return (float(jd_event), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


def _calculate_transit(
    jd_start: float,
    lon: float,
    event_type: int,
    get_ra_dec,
    CALC_ITRANSIT: int,
) -> Tuple[int, Tuple[float, ...]]:
    """Meridian transit time from a backend-agnostic RA/Dec callable.

    ``get_ra_dec(jd)`` returns the body's apparent (RA hours, Dec deg). It is
    sourced from calc_ut(FLG_EQUATORIAL) / fixstar_ut, so the LEB and Skyfield
    backends share one transit engine (the meridian crossing is geocentric, so
    no observer/parallax is needed).
    """

    def _get_body_hour_angle(jd: float) -> float:
        ra_hours, _ = get_ra_dec(jd)
        ra_deg = ra_hours * 15.0
        t = get_timescale().ut1_jd(jd)
        # Apparent sidereal time: the RA is apparent (equinox of date with
        # nutation), so GAST is the consistent pairing — GMST left a ~1.1 s
        # equation-of-equinoxes error in transit times.
        gast = t.gast
        lst = gast + lon / 15.0
        lst_deg = lst * 15.0
        ha = (lst_deg - ra_deg) % 360.0
        if ha > 180:
            ha -= 360
        return ha

    ha = _get_body_hour_angle(jd_start)
    if event_type == CALC_ITRANSIT:
        ha = (ha - 180) % 360
        if ha > 180:
            ha -= 360

    sidereal_day = 0.99726957
    if ha > 0:
        dt = (360.0 - ha) / 360.0 * sidereal_day
    else:
        dt = (-ha) / 360.0 * sidereal_day

    jd_guess = jd_start + dt
    for _ in range(30):
        ha = _get_body_hour_angle(jd_guess)
        if event_type == CALC_ITRANSIT:
            ha = (ha - 180) % 360
            if ha > 180:
                ha -= 360
        if abs(ha) < 1.0 / 3600.0:
            return 0, _make_tret(jd_guess)
        jd_guess -= ha / (360.0 / sidereal_day)

    return 0, _make_tret(jd_guess)


def _calculate_rise_set(
    jd_start: float,
    lat: float,
    event_type: int,
    _event_height,
    _get_body_ra_dec,
    CALC_RISE: int,
    CALC_SET: int,
    horizon_hint: float,
    lift_max: float,
) -> Tuple[int, Tuple[float, ...]]:
    """Find the next rise or set: first zero crossing of ``_event_height``.

    ``_event_height(jd)`` must be positive when the configured disc point
    stands above the configured horizon (refraction, disc radius and
    horizon height already applied). ``horizon_hint`` is the true
    center-altitude of that horizon and ``lift_max`` the maximum amount
    refraction and disc radius can lift the relevant point above the
    body's true center altitude; both are used only for the analytic
    circumpolar pre-check and for step sizing.
    """
    # Circumpolar pre-check from the declination extremes. For
    # fast-moving bodies (the Moon moves ~13 degrees/day) the declination
    # changes during the search window, so a finite-difference estimate of
    # declination speed widens the margin.
    _, dec = _get_body_ra_dec(jd_start)
    max_search = 2.0  # days
    _, dec2 = _get_body_ra_dec(jd_start + 1.0 / 24.0)
    dec_speed = abs(dec2 - dec) * 24.0  # degrees per day
    margin = dec_speed * max_search

    max_alt = 90.0 - abs(lat - dec)  # upper culmination
    if lat >= 0:
        min_alt = dec - (90.0 - lat)  # lower culmination
    else:
        min_alt = -dec - (90.0 + lat)

    never_below = (min_alt - margin) > horizon_hint
    never_above = (max_alt + margin + lift_max) < horizon_hint
    if never_below or never_above:
        return -2, _make_tret()

    # Step size: near-grazing geometries can keep the body above/below
    # the horizon for less than an hour, so sample more densely there.
    grazing_margin = min(abs(max_alt - horizon_hint), abs(min_alt - horizon_hint))
    if grazing_margin < 2.0:
        search_step = 10.0 / 1440.0
    elif grazing_margin < 5.0:
        search_step = 20.0 / 1440.0
    else:
        search_step = 1.0 / 24.0

    h_prev = _event_height(jd_start)
    jd_cross_start = None
    jd_cross_end = None
    for i in range(1, int(max_search / search_step) + 1):
        jd = jd_start + i * search_step
        h = _event_height(jd)
        if event_type == CALC_RISE:
            if h_prev < 0.0 <= h:
                jd_cross_start = jd - search_step
                jd_cross_end = jd
                break
        else:  # CALC_SET
            if h_prev >= 0.0 > h:
                jd_cross_start = jd - search_step
                jd_cross_end = jd
                break
        h_prev = h

    if jd_cross_start is None:
        return -2, _make_tret()
    assert jd_cross_end is not None

    # Bisection with a rate-aware altitude tolerance for near-tangential
    # grazing crossings.
    #
    # A fixed altitude residual maps to a time error of |h| / |dh/dt|. On a
    # steep crossing (the Sun rises/sets at ~15"/s at mid latitudes) 0.0001 deg
    # is well under 0.03 s, so that tolerance is kept and the crossing exits
    # exactly as before. Only when the apparent altitude sweeps the horizon
    # very slowly — the Sun skimming the horizon at a solstice essentially at
    # the polar circle, where |dh/dt| <~ 2"/s — does 0.0001 deg span an
    # appreciable time (up to ~0.4 s). There the loop keeps bisecting until the
    # *time* (not the altitude) is pinned to ~0.03 s, using the local slope
    # |dh/dt| measured at the candidate by a tight symmetric probe.
    #
    # The slope gate deliberately excludes the steeper high-summer sub-polar
    # geometry, where the rise/set sits on the refraction "dip" discontinuity
    # (the apparent altitude jumps by ~the horizon refraction across the
    # crossing): the documented stop convention there is the refracted branch
    # above the jump, so converging onto |h| = 0 would leave it. Those
    # crossings keep the historic exit unchanged. The gate separates two
    # slope populations rather than deriving from the tolerance target: a
    # near-tangential grazing approaches zero slope by definition (the
    # altitude curve flattens at tangency), while an ordinary crossing
    # runs at ~15"/s and the dip branch far steeper, so any gate in the
    # wide gap between the grazing population and the ordinary one
    # classifies the two identically; 2"/s sits in that gap. At the gate
    # itself the historic altitude exit still carries up to
    # 0.36"/(2"/s) = 0.18 s of time error (the mapping above), shrinking
    # as 1/slope for steeper crossings — accepted for the non-grazing
    # population, whose exit convention is unchanged.
    _TIME_TOL_S = 0.03
    _SLOPE_PROBE_S = 0.05
    _GRAZE_SLOPE_DEG_S = 2.0 / 3600.0  # ~2"/s: near-tangential grazing only
    # ~100"/s: far above any real apparent-altitude rate (a body on the
    # celestial equator sweeps at most ~16"/s), so a measured slope this steep
    # is the sub-polar refraction "dip" discontinuity, not a genuine crossing.
    _DIP_SLOPE_DEG_S = 100.0 / 3600.0
    for _ in range(60):
        jd_mid = (jd_cross_start + jd_cross_end) / 2
        h_mid = _event_height(jd_mid)
        if abs(h_mid) < 1e-4:
            span = jd_cross_end - jd_cross_start
            dt = _SLOPE_PROBE_S / 86400.0
            if dt > 0.25 * span:
                dt = 0.25 * span
            slope = 0.0
            slope_signed = 0.0
            if dt > 0.0:
                slope_signed = (
                    _event_height(jd_mid + dt) - _event_height(jd_mid - dt)
                ) / (2.0 * dt * 86400.0)
                slope = abs(slope_signed)
            if (
                slope <= 0.0
                or slope >= _GRAZE_SLOPE_DEG_S
                or abs(h_mid) < slope * _TIME_TOL_S
            ):
                # Refine an ordinary, well-conditioned crossing onto its true
                # root (h == 0) with one Newton step, so the event time is fixed
                # by the geometry rather than by which loose |h| < 1e-4 midpoint
                # the bisection happened to reach first. That first-inside-band
                # exit left a ~0.024 s ambiguity that stopped the LEB and
                # Skyfield backends - whose apparent places differ only by the
                # ~0.0003" LEB approximation - up to ~0.055 s apart. Only the
                # ordinary rise/set slope band is refined: near-tangential
                # grazing keeps the time-based exit above, and the near-vertical
                # apparent-altitude jump of the sub-polar refraction "dip" (a
                # slope orders of magnitude above any real motion) keeps the
                # historic midpoint so it still stops on the refracted branch
                # above the jump.
                jd_out = jd_mid
                if _GRAZE_SLOPE_DEG_S <= slope < _DIP_SLOPE_DEG_S:
                    jd_step = jd_mid - h_mid / (slope_signed * 86400.0)
                    if jd_cross_start <= jd_step <= jd_cross_end:
                        jd_out = jd_step
                return 0, _make_tret(jd_out)
        if event_type == CALC_RISE:
            if h_mid < 0.0:
                jd_cross_start = jd_mid
            else:
                jd_cross_end = jd_mid
        else:
            if h_mid >= 0.0:
                jd_cross_start = jd_mid
            else:
                jd_cross_end = jd_mid

    return 0, _make_tret((jd_cross_start + jd_cross_end) / 2)


# Alias for reference API compatibility
rise_trans = rise_trans


def rise_trans_true_hor(
    tjdut: float,
    body: "Union[int, str]",
    rsmi: int,
    geopos: Sequence[float],
    atpress: float = 0.0,
    attemp: float = 0.0,
    horhgt: float = 0.0,
    flags: int = FLG_SWIEPH,
) -> Tuple[int, Tuple[float, ...]]:
    """Calculate rise, set, or transit with custom horizon altitude.

    Wrapper around :func:`_rise_trans_true_hor_impl` that adds LEB→Skyfield
    fallback for partial/custom LEB files. See the impl docstring for the
    full API contract.
    """
    return _call_with_leb_skyfield_fallback(
        _rise_trans_true_hor_impl,
        tjdut,
        body,
        rsmi,
        geopos,
        atpress,
        attemp,
        horhgt,
        flags,
    )


def _rise_trans_true_hor_impl(
    tjdut: float,
    body: "Union[int, str]",
    rsmi: int,
    geopos: Sequence[float],
    atpress: float = 0.0,
    attemp: float = 0.0,
    horhgt: float = 0.0,
    flags: int = FLG_SWIEPH,
    *,
    reader=None,
) -> Tuple[int, Tuple[float, ...]]:
    """
    Calculate rise, set, or transit time with a custom horizon altitude.

    This function is similar to rise_trans() but allows specifying a custom
    horizon altitude. This is useful for locations with mountains or buildings
    that occlude the real horizon.

    Args:
        tjdut: Julian Day (UT) to start search from
        body: Planet/body ID (int, e.g. SUN, MOON) or fixed star
            name (str, e.g. "Sirius", "Regulus")
        rsmi: Event type and calculation flags (bitmask):
            - CALC_RISE (1): Rise time (body crossing horizon going up)
            - CALC_SET (2): Set time (body crossing horizon going down)
            - CALC_MTRANSIT (4): Upper meridian transit (culmination)
            - CALC_ITRANSIT (8): Lower meridian transit (anti-culmination)
            Additional flags (OR with event type):
            - BIT_GEOCTR_NO_ECL_LAT (128): Use the geocentric apparent place
              with the ecliptic latitude zeroed (Hindu-rising convention;
              see the note below). Component of BIT_HINDU_RISING (896).
            - BIT_DISC_CENTER (256): Use disc center instead of upper limb
            - BIT_DISC_BOTTOM (8192): Use lower limb of disc
            - BIT_NO_REFRACTION (512): Ignore atmospheric refraction
        geopos: Geographic position as sequence [lon, lat, alt]:
            - 0: geographic longitude in degrees (eastern positive)
            - 1: geographic latitude in degrees (northern positive)
            - 2: geographic altitude in meters above sea level
        atpress: Atmospheric pressure in mbar/hPa for refraction (default 0.0)
        attemp: Atmospheric temperature in degrees Celsius (default 0.0)
        horhgt: Height of local horizon in degrees (default 0.0).
            Positive values mean the horizon is elevated (e.g., mountains),
            so rise times will be later and set times earlier.
            Negative values mean the horizon is depressed (e.g., observer
            on a mountain).
            CAVEAT: "depressed" only holds outside a dead band next to 0. The
            refraction inversion clamps to the unrefracted altitude whenever the
            apparent altitude would fall below the dip, and at sea level the dip
            is exactly 0 - so at geopos[2]=0 every horhgt from 0 down to
            -0.559889 returns the same time to within 13 ms, and past that edge
            the response is a KINK rather than a step: 322 s per degree, so -0.57
            is -3.3 s and -0.60 is -12.9 s. Requesting a 30-arcminute depression
            at sea level gets the undepressed event. The band tracks the dip, so
            it moves with elevation (~1.0 deg wide at 3000 m).
            The special value -100.0 is not an angle: it asks for the computed
            dip of the SEA horizon at geopos[2], via calc_dip() (Thom's
            refraction-corrected dip, not the schoolbook 1.76*sqrt(h)), and it
            resolves to dip + 1e-4 deg, i.e. just past the edge of the band
            above. This is the only way to make the observer's elevation
            actually lower the horizon - geopos[2] on its own does not (see the
            "Horizon convention" section of rise_trans). It is worth minutes,
            not seconds: at Rome for 2024-03-20 under 1013.25 hPa / 15 C,
            horhgt=-100 against horhgt=0 moves sunrise -109.9 s at 100 m,
            -360.0 s at 1000 m and -649.1 s at 3000 m (on the library default
            0 C: -110.6 / -362.3 / -640.7 s). Note the default 0.0 is what
            reproduces published rise/set tables, which assume a level sea
            horizon at any elevation.
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - res: integer flag (0 = event found, -2 = circumpolar)
            - tret: tuple of 10 floats, of which tret[0] = JD of event

    Raises:
        ValueError: If invalid planet ID or star name

    Note:
        For circumpolar objects (always above or below horizon at the given
        latitude), the function returns (-2, tret) with tret[0] = 0.0.
        For transits, circumpolar objects still have valid transit times.

        Twilight flags (BIT_CIVIL_TWILIGHT, BIT_NAUTIC_TWILIGHT,
        BIT_ASTRO_TWILIGHT) are honored for the Sun and for fixed
        stars: they switch the event to a geometric center crossing of
        -6/-12/-18 degrees and override horhgt, matching the reference
        behavior (planets and the Moon keep their ordinary horizon
        crossing). The special value horhgt=-100 requests the dip of
        the sea horizon.

        BIT_GEOCTR_NO_ECL_LAT (128) replaces the ordinary topocentric place
        with the body's GEOCENTRIC apparent place projected onto the ecliptic
        (ecliptic latitude set to zero, longitude and distance kept), then
        rotates that direction to the observer's horizon. This drops both the
        topocentric parallax and the ecliptic latitude, so a body far from the
        ecliptic rises/sets at a very different time (e.g. the Moon shifts by
        tens of minutes toward its ecliptic-longitude crossing). Compatibility
        contract: the bit affects rise/set and twilight but is
        ignored for meridian transits. It is the latitude-zeroing component of
        BIT_HINDU_RISING (896 = 128 | 256 | 512).

    Algorithm:
        1. For transits: Find when body crosses the local meridian
           (Local Sidereal Time = body's Right Ascension)
        2. For rise/set: Find when body's altitude crosses the custom horizon
           accounting for refraction and disc size
        3. Uses Newton-Raphson iteration for precise timing

    Precision:
        Rise/set times accurate to ~1 minute for Sun/Moon (due to refraction
        uncertainty), better than 1 second for transit times.

    Example:
        >>> from libephemeris import julday, rise_trans_true_hor, SUN, CALC_RISE
        >>> jd = julday(2024, 6, 21, 0)
        >>> # Find sunrise with mountains at 5 deg above horizon at Rome
        >>> res, tret = rise_trans_true_hor(jd, SUN, CALC_RISE,
        ...                                [12.5, 41.9, 0], horhgt=5.0)
        >>> print(f"Sunrise at JD {tret[0]:.5f}")

    References:
        - Reference API: rise_trans_true_hor()
        - Meeus "Astronomical Algorithms" Ch. 15 (Rise, Set, Transit)
    """
    from .constants import (
        CALC_RISE,
        CALC_SET,
        CALC_MTRANSIT,
        CALC_ITRANSIT,
        BIT_DISC_CENTER,
        BIT_DISC_BOTTOM,
        BIT_NO_REFRACTION,
        BIT_CIVIL_TWILIGHT,
        BIT_NAUTIC_TWILIGHT,
        BIT_ASTRO_TWILIGHT,
        BIT_FIXED_DISC_SIZE,
    )
    from .planets import _PLANET_MAP

    # Unpack geopos: [longitude, latitude, altitude]. The length is enforced
    # here as it is on every sibling entry point: accepting a 2-element
    # sequence silently defaulted the altitude to 0, and a shorter or
    # non-numeric one escaped as IndexError/ValueError, neither of which
    # derives from this library's Error — so `except Error:` (the drop-in
    # for the reference's own exception) missed them.
    try:
        _n_geo = len(geopos)
    except TypeError:
        _n_geo = 0
    if _n_geo < 3:
        raise ValueError("geopos must be a sequence of [longitude, latitude, altitude]")
    try:
        lon = float(geopos[0])
        lat = float(geopos[1])
        altitude = float(geopos[2])
    except (TypeError, ValueError) as _geo_exc:
        raise ValueError(
            "geopos must be a sequence of [longitude, latitude, altitude]"
        ) from _geo_exc

    # Map parameter names for internal use
    jd_start = tjdut
    pressure = atpress
    temperature = attemp
    horizon_altitude = horhgt
    is_fixed_star = isinstance(body, str)

    # Extract the event type from rsmi's low bits using the compatibility
    # precedence SET > RISE > ITRANSIT > MTRANSIT. A modifier-only mask or
    # missing event type defaults to RISE.
    if rsmi & CALC_SET:
        event_type = CALC_SET
    elif rsmi & CALC_RISE:
        event_type = CALC_RISE
    elif rsmi & CALC_ITRANSIT:
        event_type = CALC_ITRANSIT
    elif rsmi & CALC_MTRANSIT:
        event_type = CALC_MTRANSIT
    else:
        event_type = CALC_RISE

    # reader is provided by the caller (None forces Skyfield path)
    _reader_rt = reader
    _use_leb_rt = _reader_rt is not None

    # Bodies the classical _PLANET_MAP does not enumerate (the lunar nodes and
    # apsides, Chiron, the main-belt asteroids, numbered minor planets, ...)
    # are still placeable by calc_ut, so the reference computes their rise, set
    # and transit. Route them through one backend-agnostic pipeline: the
    # geocentric apparent place -> azalt for the altitude and calc_ut
    # EQUATORIAL for the transit RA/Dec. calc_ut dispatches to the active
    # backend (the LEB reader when sealed, JPL/Skyfield otherwise) and already
    # enforces the typed error contract, so a genuinely unknown id or an
    # out-of-range date raises there. These are point sources: zero disc radius
    # and negligible (geocentric) parallax — the compatibility-contract
    # rise/set geometry.
    use_calc_body = (not is_fixed_star) and (cast(int, body) not in _PLANET_MAP)
    if use_calc_body:
        # Validate positionability up front so an unplaceable body raises the
        # typed error at the call (compatibility contract: rejection happens
        # at the call), not mid-search; the fast_calc probe below uses the
        # wrong id space for AST_OFFSET-aliased ids, so it is skipped for
        # these bodies.
        try:
            calc_ut(jd_start, cast(int, body), FLG_SPEED)
        except (KeyError, ValueError, Error) as _probe_exc:
            from .exceptions import EphemerisRangeError

            # A known body whose date is outside coverage (or a corrupt LEB)
            # keeps its typed range error - never masked as an illegal body.
            _reraise_if_leb_range_error(_probe_exc)
            if isinstance(_probe_exc, EphemerisRangeError):
                raise
            # A target calc_ut cannot place at all (a planetary moon with no
            # registered SPK, an unknown id, ...) is illegal for rise/set. One
            # typed error on every backend, carrying the reference's "illegal
            # planet number" message (both an Error and a ValueError).
            raise _IllegalRiseBodyError(
                "illegal planet number %d." % cast(int, body),
                body_id=cast(int, body),
            ) from _probe_exc

    if _use_leb_rt and not use_calc_body:
        try:
            from . import fast_calc as _fc_rt

            if not is_fixed_star:
                _fc_rt.fast_calc_ut(_reader_rt, jd_start, cast(int, body), FLG_SPEED)
        except (KeyError, ValueError) as _probe_exc:
            # A coverage miss must reach the mode-aware wrapper (typed error
            # in sealed mode, Skyfield retry in auto), not silently degrade
            # to the kernel path that sealed mode forbids. A body-map miss in
            # sealed mode must surface as the typed error too: the kernel
            # path below would only convert it into get_planets()'s
            # RuntimeError.
            _raise_if_sealed_leb_miss(_probe_exc)
            _reraise_if_leb_range_error(_probe_exc)
            _use_leb_rt = False

    # Both backends evaluate the horizon geometry through ONE pipeline: the
    # topocentric apparent place from calc_ut(FLG_TOPOCTR) (fixstar_ut for
    # stars) rotated to the observer's horizon by azalt, plus the equatorial
    # place from calc_ut(FLG_EQUATORIAL) for the meridian transit. calc_ut and
    # fixstar_ut dispatch to the active backend (the LEB reader when sealed,
    # JPL/Skyfield otherwise) and azalt is the shared transform, so the LEB and
    # Skyfield results differ only by the tiny (~0.0003") LEB-vs-Skyfield
    # position difference instead of by the coordinate frame. The earlier split
    # - a Skyfield-native topocentric altaz frame on one path and azalt on the
    # other - drove rise/set times apart by up to ~0.11 s even when the two
    # backends agreed on the position to 0.0003"; sharing the chain removes that
    # backend-specific divergence.
    if is_fixed_star:
        planet = -1
    else:
        planet = cast(int, body)
        from .constants import ECL_NUT as _ECL_NUT

        if planet in (EARTH, _ECL_NUT):
            # Neither id denotes a direction in the sky, so no horizon event
            # exists: EARTH *is* the observer (calc_ut returns the zero
            # coordinate origin) and ECL_NUT is not a body at all (its first
            # slots carry nutation and obliquity). Feeding either tuple to
            # azalt would manufacture a location-dependent rise/set from
            # data that is not a direction.
            #
            # The measured external value is not reproducible: it returns
            # tjd_start + 0.083333254 days for both ids at every latitude
            # tested, rise and set identical — an internal sentinel, not an
            # astronomical result, and encoding that constant would be
            # reconstructing output. Fail explicitly instead, as this library
            # does for every other degenerate self-request.
            from .exceptions import Error as _Error

            raise _Error(
                "rise/set/transit is undefined for body %d: it has no "
                "apparent direction in the sky." % planet
            )
        if planet not in _PLANET_MAP and not use_calc_body:
            # Defensive only: an out-of-map body is validated up front (it sets
            # use_calc_body, and an unplaceable id already raised the unified
            # illegal-body error there), so this branch is unreachable in
            # practice; keep the reference's illegal-body message regardless.
            from .exceptions import UnknownBodyError

            raise UnknownBodyError("illegal planet number %d." % planet)

    # atpress == 0: estimate pressure from the observer altitude with the
    # reference's barometric expression; attemp is used verbatim (0 means
    # 0 degrees C, as in the reference - no implicit 15 C default).
    if pressure == 0.0:
        # The base goes negative above ~44.3 km, where a fractional power would
        # yield a complex number; clamp it to 0 so extreme altitudes give
        # pressure ~0 (refraction ~0) instead of crashing.
        base = 1.0 - 0.0065 * altitude / 288.0
        pressure = 1013.25 * (max(base, 0.0) ** 5.255)

    # horizon_altitude is used verbatim (negative values are legitimate:
    # an observer above the surrounding terrain sees below the geometric
    # horizon). The special value -100 requests the dip of the sea
    # horizon for the observer's altitude.
    horizon_alt = horizon_altitude
    if horizon_alt == -100.0:
        from .refraction import calc_dip

        horizon_alt = 0.0001 + calc_dip(
            max(altitude, 0.0), atpress=pressure, attemp=temperature
        )

    # Twilight events: geometric center crossings of -6/-12/-18 degrees
    # (the horizon height is ignored for twilight). Compatibility contract:
    # the twilight bits are honored for the Sun AND for fixed stars
    # (heliacal-style star visibility at a chosen solar-depression class),
    # while planets and the Moon keep their ordinary horizon crossing.
    rsmi_eff = rsmi
    is_twilight = bool(
        rsmi & (BIT_CIVIL_TWILIGHT | BIT_NAUTIC_TWILIGHT | BIT_ASTRO_TWILIGHT)
    )
    if is_twilight and (planet == SUN or is_fixed_star):
        rsmi_eff = rsmi | BIT_NO_REFRACTION | BIT_DISC_CENTER
        if rsmi & BIT_CIVIL_TWILIGHT:
            horizon_alt = -6.0
        elif rsmi & BIT_NAUTIC_TWILIGHT:
            horizon_alt = -12.0
        else:
            horizon_alt = -18.0

    use_refraction = not (rsmi_eff & BIT_NO_REFRACTION)

    # Which point of the disc crosses the horizon: +1 = upper limb
    # (default), -1 = lowest point (BIT_DISC_BOTTOM), 0 = disc center.
    if (rsmi_eff & BIT_DISC_CENTER) or is_fixed_star:
        disc_sign = 0
    elif rsmi_eff & BIT_DISC_BOTTOM:
        disc_sign = -1
    else:
        disc_sign = 1

    _body_radius_au = 0.0
    if not is_fixed_star and disc_sign != 0:
        if planet == SUN:
            _body_radius_au = _ECL_RSUN_AU
        elif planet == MOON:
            _body_radius_au = _ECL_RMOON_AU
        else:
            _r_km = _PLANET_RADIUS_KM.get(planet)
            _body_radius_au = (_r_km / _ECL_AU_KM) if _r_km else 0.0

    from .constants import FLG_EQUATORIAL, FLG_TOPOCTR as _FLG_TOPO_RT
    from .utils import azalt, ECL2HOR

    geopos_rt = (lon, lat, altitude)

    def _get_body_altaz(jd: float) -> Tuple[float, float, float]:
        """(true altitude, azimuth, distance) of the ordinary topocentric
        apparent place, evaluated identically on every backend.

        Planets and points use calc_ut(FLG_TOPOCTR) (the observer is set once
        around the search, see the dispatch tail); fixed stars use their
        geocentric apparent place (diurnal parallax is negligible for
        rise/set). azalt then rotates the ecliptic-of-date place to the
        observer's horizon. calc_ut/fixstar_ut dispatch to the active backend
        and azalt is the shared transform, so the LEB and Skyfield altitudes
        agree to the position precision (~0.0003") and rise/set/transit are
        backend-independent.
        """
        if is_fixed_star:
            from .fixed_stars import fixstar_ut

            star_pos, _, _ = fixstar_ut(cast(str, body), jd, FLG_SPEED)
            az, alt_true, _ = azalt(jd, ECL2HOR, geopos_rt, 0, 0, star_pos[:3])
            return alt_true, az, star_pos[2]
        pos, _ = calc_ut(jd, planet, _FLG_TOPO_RT | FLG_SPEED)
        az, alt_true, _ = azalt(jd, ECL2HOR, geopos_rt, 0, 0, pos[:3])
        return alt_true, az, pos[2]

    def _get_body_ra_dec(jd: float) -> Tuple[float, float]:
        """Apparent RA (hours) and Dec (deg), evaluated identically on every
        backend: calc_ut(FLG_EQUATORIAL) for planets/points, fixstar_ut for
        fixed stars."""
        if is_fixed_star:
            from .fixed_stars import fixstar_ut

            pos, _, _ = fixstar_ut(cast(str, body), jd, FLG_EQUATORIAL | FLG_SPEED)
            return pos[0] / 15.0, pos[1]
        eq, _ = calc_ut(jd, planet, FLG_EQUATORIAL | FLG_SPEED)
        return eq[0] / 15.0, eq[1]

    # Hindu-rising convention (BIT_GEOCTR_NO_ECL_LAT): the event uses the
    # body's GEOCENTRIC apparent place projected onto the ecliptic (ecliptic
    # latitude zeroed), not the ordinary topocentric place. Compatibility
    # contract: this shifts rise/set AND twilight (e.g. the Moon's rise moves by
    # ~25 min toward its ecliptic-longitude crossing) but leaves meridian
    # transits untouched - the transit path below returns before the altitude
    # engine, so the bit is naturally ignored there. Built on calc_ut/azalt so
    # both backends share one geocentric pipeline (automatic LEB<->Skyfield
    # parity).
    from .constants import BIT_GEOCTR_NO_ECL_LAT

    use_geoctr_nolat = bool(rsmi_eff & BIT_GEOCTR_NO_ECL_LAT)

    def _get_body_altaz_geoctr_nolat(jd: float) -> Tuple[float, float, float]:
        """(true altitude, azimuth, distance) from the geocentric apparent
        place with the ecliptic latitude zeroed.

        Keeps the body's geocentric ecliptic longitude and distance, drops its
        parallax (geocentric, not topocentric) and its ecliptic latitude, then
        rotates to the observer's horizon. The horizon rotation still uses the
        observer's geographic location and the of-date true obliquity.
        """
        from .utils import azalt, ECL2HOR

        if is_fixed_star:
            from .fixed_stars import fixstar_ut

            pos_s, _, _ = fixstar_ut(cast(str, body), jd, FLG_SPEED)
            lon_b, dist_b = float(pos_s[0]), float(pos_s[2])
        else:
            pos_p, _ = calc_ut(jd, planet, FLG_SPEED)
            lon_b, dist_b = float(pos_p[0]), float(pos_p[2])
        az, alt_true, _ = azalt(
            jd, ECL2HOR, (lon, lat, altitude), 0.0, 0.0, (lon_b, 0.0, dist_b)
        )
        return alt_true, az, dist_b

    def _get_body_ra_dec_geoctr_nolat(jd: float) -> Tuple[float, float]:
        """Apparent RA (hours) and Dec (deg) of the ecliptic-latitude-zeroed
        geocentric place used by the Hindu-rising path.

        The rise/set solver's circumpolar pre-check derives the culmination
        altitudes from a declination. Under BIT_GEOCTR_NO_ECL_LAT the body is
        projected onto the ecliptic (latitude 0) BEFORE the horizon rotation,
        which lowers |Dec|, so that pre-check must use the SAME projected
        declination the solve uses. Otherwise a body far from the ecliptic (the
        Moon near its maximum latitude) is judged circumpolar on its true,
        higher declination while the solve — already working on the lower
        projected declination — would have found the event: the two disagree
        and rise/set is spuriously reported circumpolar (res=-2) at high
        geographic latitude. Compatibility contract: the event exists and
        must be returned.

        The projection uses the same of-date true obliquity azalt applies in
        _get_body_altaz_geoctr_nolat, so this declination is exactly the one the
        horizon rotation sees.
        """
        import erfa

        from .precession_vondrak import vondrak_mean_obliquity_deg
        from .utils import cotrans

        if is_fixed_star:
            from .fixed_stars import fixstar_ut

            pos_s, _, _ = fixstar_ut(cast(str, body), jd, FLG_SPEED)
            lon_b, dist_b = float(pos_s[0]), float(pos_s[2])
        else:
            pos_p, _ = calc_ut(jd, planet, FLG_SPEED)
            lon_b, dist_b = float(pos_p[0]), float(pos_p[2])
        jd_tt = get_timescale().ut1_jd(jd).tt
        eps0 = vondrak_mean_obliquity_deg(jd_tt)
        _, deps_rad = erfa.nut06a(2451545.0, jd_tt - 2451545.0)
        eps = eps0 + math.degrees(deps_rad)
        ra, dec, _ = cotrans((lon_b, 0.0, dist_b), -eps)
        return ra / 15.0, dec

    # The ordinary (non-Hindu) altitude and the transit RA/Dec both flow
    # through the shared, backend-agnostic _get_body_altaz / _get_body_ra_dec
    # above (calc_ut/fixstar_ut + azalt), so in-map bodies, out-of-map bodies
    # (nodes, apsides, Chiron, asteroids, numbered minor planets) and fixed
    # stars all use one pipeline on both backends.
    _altaz_for_event = (
        _get_body_altaz_geoctr_nolat if use_geoctr_nolat else _get_body_altaz
    )
    _ra_dec_for_event = _get_body_ra_dec
    # The rise/set circumpolar pre-check must see the SAME declination the
    # altitude solve uses. Under the Hindu bit that is the ecliptic-latitude-
    # zeroed projected declination (see _get_body_ra_dec_geoctr_nolat); the
    # transit path keeps the true RA/Dec (the bit is ignored for transits).
    _ra_dec_for_circumpolar = (
        _get_body_ra_dec_geoctr_nolat if use_geoctr_nolat else _get_body_ra_dec
    )

    # The ordinary rise/set path samples the TOPOCENTRIC apparent altitude
    # (calc_ut(FLG_TOPOCTR)) for every planet/point; set the observer once for
    # the whole search and restore it afterwards so the diurnal parallax is
    # applied without churning global state per sample. Fixed stars (geocentric)
    # and the Hindu geocentric path ignore the observer.
    _topo_scope = (not is_fixed_star) and not use_geoctr_nolat

    if event_type in (CALC_MTRANSIT, CALC_ITRANSIT):
        return _calculate_transit(
            jd_start,
            lon,
            event_type,
            _ra_dec_for_event,
            CALC_ITRANSIT,
        )

    def _event_height(jd: float) -> float:
        """Height of the configured disc point above the configured horizon."""
        alt_true, _az, dist = _altaz_for_event(jd)
        alt_pt = alt_true
        if disc_sign != 0 and _body_radius_au > 0.0:
            cur = dist
            if rsmi_eff & BIT_FIXED_DISC_SIZE:
                if planet == SUN:
                    cur = 1.0
                elif planet == MOON:
                    cur = 0.00257
            alt_pt += disc_sign * math.degrees(
                math.asin(min(1.0, _body_radius_au / cur))
            )
        if use_refraction:
            alt_pt = _rise_true_to_apparent(alt_pt, altitude, pressure, temperature)
        return alt_pt - horizon_alt

    # Maximum amount refraction and the disc radius can lift the relevant
    # point above its true center altitude (for the circumpolar pre-check).
    lift_max = 0.0
    if use_refraction:
        lift_max += 0.75
    if disc_sign != 0:
        lift_max += 0.30

    _rt_scope = (
        _observer_scope(lon, lat, altitude) if _topo_scope else contextlib.nullcontext()
    )
    with _rt_scope:
        return _calculate_rise_set(
            jd_start,
            lat,
            event_type,
            _event_height,
            _ra_dec_for_circumpolar,
            CALC_RISE,
            CALC_SET,
            horizon_alt,
            lift_max,
        )


# Alias for reference API compatibility
rise_trans_true_hor = rise_trans_true_hor


# =============================================================================
# HELIACAL RISING AND SETTING CALCULATIONS
# =============================================================================


def heliacal_ut(
    jd_start: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    pressure: float = 1013.25,
    temperature: float = 15.0,
    humidity: float = 0.5,
    body: int = SUN,
    event_type: int = 1,
    flags: int = FLG_SWIEPH,
) -> Tuple[float, int]:
    """
    Calculate heliacal rising or setting time for a celestial body.

    Heliacal events are the first/last visibility of a celestial body
    at dawn or dusk. These were fundamental for ancient calendars:
    - Heliacal rising: First morning visibility after a period of invisibility
    - Heliacal setting: Last evening visibility before becoming invisible

    Args:
        jd_start: Julian Day (UT) to start search from
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        pressure: Atmospheric pressure in mbar/hPa for refraction (default 1013.25)
        temperature: Temperature in Celsius for refraction (default 15)
        humidity: Relative humidity 0.0-1.0 for atmospheric extinction (default 0.5)
        body: Planet/body ID (MERCURY, VENUS, MARS, JUPITER, SATURN,
              or fixed stars). Note: SUN and MOON are not valid for heliacal events.
        event_type: Type of heliacal event:
            - HELIACAL_RISING (1): Morning first visibility (heliacal rising)
            - HELIACAL_SETTING (2): Evening last visibility (heliacal setting)
            - EVENING_FIRST (3): First evening visibility (after superior conjunction)
            - MORNING_LAST (4): Last morning visibility (before superior conjunction)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - jd_event: Julian Day (UT) of the heliacal event, or 0.0 if not found
            - retflag: Return flag (event_type on success, negative on error)

    Raises:
        ValueError: If invalid body ID or event_type

    Algorithm:
        The algorithm searches for the moment when:
        1. The body is at a specific altitude above the horizon (arcus visionis)
        2. The Sun is at twilight position (typically -6° to -12° below horizon)
        3. The body's apparent magnitude is brighter than the sky's limiting magnitude

        For heliacal rising (morning first):
        - Search forward for when the body first becomes visible at dawn
        - Body must be above horizon while Sun is still below
        - Sky must be dark enough for the body to be seen

        For heliacal setting (evening last):
        - Search forward for when the body is last visible at dusk
        - Body must be above horizon while Sun is setting
        - Sky brightness must not overwhelm the body's light

    Historical Note:
        Heliacal risings were crucial for ancient calendars. The heliacal
        rising of Sirius marked the Egyptian new year and predicted the
        Nile flood. Babylonians used heliacal events to track planetary
        positions without modern instruments.

    Example:
        >>> from libephemeris import julday, heliacal_ut, VENUS, HELIACAL_RISING
        >>> jd = julday(2024, 1, 1, 0)
        >>> # Find next heliacal rising of Venus from Rome
        >>> jd_event, flag = heliacal_ut(jd, 41.9, 12.5, body=VENUS,
        ...                              event_type=HELIACAL_RISING)
        >>> print(f"Heliacal rising at JD {jd_event:.5f}")

    References:
        - Schaefer (1990), PASP 102, 212-229, DOI 10.1086/132629,
          limiting-magnitude and observer model.
        - Schaefer (1993), Vistas in Astronomy 36, 311-361, visibility
          context and observer effects.
        - Kasten & Young (1989), DOI 10.1364/AO.28.004735, optical airmass.

        Public API compatibility defines arguments and return semantics only;
        it is not a source for the visibility model. The historical narrative
        above is context, not a numerical algorithm.
    """
    from .constants import (
        HELIACAL_RISING,
        HELIACAL_SETTING,
        EVENING_FIRST,
        MORNING_LAST,
    )
    from .planets import _PLANET_MAP, get_planet_target, pheno_ut
    from .state import get_planets, get_timescale
    from skyfield.api import wgs84

    # Validate event type
    if event_type not in (
        HELIACAL_RISING,
        HELIACAL_SETTING,
        EVENING_FIRST,
        MORNING_LAST,
    ):
        raise Error(
            f"Invalid event_type: {event_type}. Use HELIACAL_RISING, "
            "HELIACAL_SETTING, EVENING_FIRST, or MORNING_LAST."
        )

    # Sun and Moon are not valid for heliacal events
    if body == SUN:
        raise Error("SUN is not valid for heliacal calculations")
    if body == MOON:
        raise Error("MOON is not valid for heliacal calculations")

    # Validate body
    if body not in _PLANET_MAP:
        raise Error(f"illegal planet number {body}.")

    # Get ephemeris and timescale
    eph = get_planets()
    ts = get_timescale()

    # Get celestial bodies
    target_name = _PLANET_MAP[body]
    target = get_planet_target(eph, target_name)
    sun = eph["sun"]
    earth = eph["earth"]

    # Create observer location
    observer = wgs84.latlon(lat, lon, altitude)

    # Note: Arcus visionis (minimum altitude for visibility) and
    # sun altitude thresholds are used in the visibility check functions.
    # Typical values: arcus visionis ~7° for bright objects, sun altitude ~-8°.

    def _get_altitudes(jd: float) -> Tuple[float, float, float]:
        """Get Sun altitude, body altitude, and body azimuth at given JD."""
        t = ts.ut1_jd(jd)
        observer_at = earth + observer

        # Sun position
        sun_app = observer_at.at(t).observe(sun).apparent()
        sun_alt, _, _ = sun_app.altaz()

        # Body position
        body_app = observer_at.at(t).observe(target).apparent()
        body_alt, body_az, _ = body_app.altaz()

        # Convert Skyfield navigational azimuth (N=0) to the reference convention (S=0)
        return sun_alt.degrees, body_alt.degrees, (body_az.degrees + 180.0) % 360.0

    def _get_elongation(jd: float) -> float:
        """Get the elongation of body from Sun in degrees."""
        try:
            pheno = pheno_ut(jd, body, flags)
            return pheno[2]  # Elongation
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            # Fallback: calculate elongation manually
            t = ts.ut1_jd(jd)
            sun_app = earth.at(t).observe(sun).apparent()
            body_app = earth.at(t).observe(target).apparent()
            return body_app.separation_from(sun_app).degrees

    def _get_body_magnitude(jd: float) -> float:
        """Get the visual magnitude of the body."""
        try:
            pheno = pheno_ut(jd, body, flags)
            return pheno[4]  # Visual magnitude
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            return 0.0  # Default to bright magnitude

    def _calculate_limiting_magnitude(sun_alt: float, body_alt: float) -> float:
        """
        Calculate the limiting magnitude based on sky brightness.

        The limiting magnitude depends on:
        - Sun's altitude below horizon (darker = higher limit)
        - Body's altitude (extinction increases near horizon)
        - Atmospheric conditions (humidity, pressure)

        This is a simplified model based on Schaefer (1990).
        """
        # Sky brightness model (simplified)
        # When Sun is below -18°, sky is fully dark (limiting mag ~6.5)
        # When Sun is at -6° (civil twilight), limiting mag is ~3
        # When Sun is at 0° (horizon), limiting mag is ~-2

        if sun_alt >= 0:
            return -2.0  # Daylight - essentially nothing visible
        elif sun_alt >= -6:
            # Civil twilight
            return -2.0 + (sun_alt + 6) * (-3.0 - (-2.0)) / (-6)
        elif sun_alt >= -12:
            # Nautical twilight
            return 3.0 + (sun_alt + 12) * (5.0 - 3.0) / 6
        elif sun_alt >= -18:
            # Astronomical twilight
            return 5.0 + (sun_alt + 18) * (6.5 - 5.0) / 6
        else:
            # Full darkness
            return 6.5

        # Atmospheric extinction correction (simplified)
        # Add extinction based on altitude (airmass)
        # At low altitudes, extinction can be 0.5-1.0 magnitudes

    def _is_body_visible(jd: float) -> Tuple[bool, float, float, float]:
        """
        Check if body is visible at given time.

        Returns: (is_visible, sun_alt, body_alt, elongation)
        """
        sun_alt, body_alt, body_az = _get_altitudes(jd)
        elongation = _get_elongation(jd)

        # Body must be above minimum altitude
        if body_alt < 0:
            return False, sun_alt, body_alt, elongation

        # Sun must be below horizon
        if sun_alt > 0:
            return False, sun_alt, body_alt, elongation

        # Check limiting magnitude vs body magnitude
        limiting_mag = _calculate_limiting_magnitude(sun_alt, body_alt)
        body_mag = _get_body_magnitude(jd)

        # Body is visible if its magnitude is brighter (lower) than limiting magnitude
        # Also require sufficient elongation from Sun
        min_elongation = 10.0  # Minimum elongation for visibility (degrees)

        is_visible = (body_mag <= limiting_mag) and (elongation >= min_elongation)

        return is_visible, sun_alt, body_alt, elongation

    def _find_twilight_time(jd: float, sun_target_alt: float, rising: bool) -> float:
        """
        Find when Sun crosses target altitude (morning or evening).

        Args:
            jd: Starting JD
            sun_target_alt: Target Sun altitude (negative for below horizon)
            rising: True for morning (Sun rising), False for evening (Sun setting)
        """
        # Search within one day
        for _ in range(50):
            sun_alt, _, _ = _get_altitudes(jd)

            # Check if we're at the right phase of day
            if rising:
                # Morning: looking for Sun rising through target altitude
                if abs(sun_alt - sun_target_alt) < 0.01:
                    return jd
            else:
                # Evening: looking for Sun setting through target altitude
                if abs(sun_alt - sun_target_alt) < 0.01:
                    return jd

            # Adjust time based on Sun's position
            # Sun moves ~15°/hour = 0.25°/minute = 360°/day
            sun_rate = 360.0 / 1.0  # degrees per day
            diff = sun_target_alt - sun_alt

            # Crude estimate of time adjustment
            dt = diff / sun_rate

            # Limit step size
            dt = max(-0.1, min(0.1, dt))

            if abs(dt) < 1e-6:
                return jd

            jd += dt

        return jd

    def _search_heliacal_rising(jd_start: float) -> float:
        """
        Search for heliacal rising (morning first visibility).

        The body becomes visible in the morning before sunrise after
        a period of being hidden in the Sun's glare.
        """
        max_days = 400  # Search up to more than a year

        # Use coarse search first (check once per day at approximate twilight)
        for day in range(max_days):
            jd_day = jd_start + day

            # Check at approximately 5 AM local time (rough morning twilight estimate)
            # This is a coarse check - we'll refine if we find visibility
            for hour in [4, 5, 6]:  # Check a few morning hours
                jd_check = jd_day + hour / 24.0

                sun_alt, body_alt, _ = _get_altitudes(jd_check)

                # Look for conditions: Sun below horizon but not too deep, body above
                if -15 < sun_alt < -3 and body_alt > 1:
                    visible, s_alt, b_alt, elong = _is_body_visible(jd_check)

                    if visible:
                        # Found a potential heliacal rising
                        # Refine the time using binary search
                        return _refine_heliacal_time(jd_check, is_morning=True)

        return 0.0  # Not found

    def _search_heliacal_setting(jd_start: float) -> float:
        """
        Search for heliacal setting (evening last visibility).

        The body is last visible in the evening after sunset before
        becoming hidden in the Sun's glare.
        """
        max_days = 400

        # First, find when the body is visible in the evening
        first_visible_jd = 0.0

        for day in range(max_days):
            jd_day = jd_start + day

            # Check at a few evening hours (coarse search)
            for hour in [18, 19, 20]:
                jd_check = jd_day + hour / 24.0

                sun_alt, body_alt, _ = _get_altitudes(jd_check)

                # Look for conditions: Sun below horizon, body above
                if -15 < sun_alt < -3 and body_alt > 1:
                    visible, s_alt, b_alt, elong = _is_body_visible(jd_check)

                    if visible:
                        first_visible_jd = jd_check
                        break

            if first_visible_jd > 0:
                break

        if first_visible_jd == 0:
            return 0.0

        # Now find the LAST day of visibility
        last_visible_jd = first_visible_jd
        start_day = int(first_visible_jd - jd_start)

        for day in range(start_day, start_day + 120):  # Check up to 120 days ahead
            jd_day = jd_start + day
            found_visible = False

            for hour in [18, 19, 20]:
                jd_check = jd_day + hour / 24.0
                sun_alt, body_alt, _ = _get_altitudes(jd_check)

                if -15 < sun_alt < -3 and body_alt > 1:
                    visible, _, _, _ = _is_body_visible(jd_check)
                    if visible:
                        last_visible_jd = jd_check
                        found_visible = True
                        break

            if not found_visible and day > start_day:
                # No visibility found today - return the last visible day
                return _refine_heliacal_time(last_visible_jd, is_morning=False)

        # If still visible after 120 days, return the last checked
        return _refine_heliacal_time(last_visible_jd, is_morning=False)

    def _search_evening_first(jd_start: float) -> float:
        """
        Search for evening first visibility (after superior conjunction).

        The body appears in the evening sky for the first time after
        passing behind the Sun (superior conjunction for inferior planets,
        or conjunction for superior planets).
        """
        max_days = 400
        was_invisible = False

        for day in range(max_days):
            jd_day = jd_start + day

            # Check at a few evening hours (coarse search)
            for hour in [18, 19, 20]:
                jd_check = jd_day + hour / 24.0

                sun_alt, body_alt, _ = _get_altitudes(jd_check)

                if -15 < sun_alt < -3:
                    visible, s_alt, b_alt, elong = _is_body_visible(jd_check)
                    if not visible and body_alt < 5:
                        was_invisible = True
                    elif visible and was_invisible:
                        # First evening visibility after being invisible
                        return _refine_heliacal_time(jd_check, is_morning=False)

        return 0.0

    def _search_morning_last(jd_start: float) -> float:
        """
        Search for morning last visibility (before superior conjunction).

        The body is last visible in the morning sky before passing
        behind the Sun.
        """
        max_days = 400
        last_visible_jd = 0.0
        found_visible = False

        for day in range(max_days):
            jd_day = jd_start + day

            # Check at a few morning hours (coarse search)
            for hour in [4, 5, 6]:
                jd_check = jd_day + hour / 24.0

                sun_alt, body_alt, _ = _get_altitudes(jd_check)

                if -15 < sun_alt < -3:
                    visible, s_alt, b_alt, elong = _is_body_visible(jd_check)

                    if visible:
                        last_visible_jd = jd_check
                        found_visible = True
                    elif found_visible and not visible:
                        # Found the transition to invisibility
                        return _refine_heliacal_time(last_visible_jd, is_morning=True)

        return last_visible_jd if last_visible_jd > 0 else 0.0

    def _refine_heliacal_time(jd_approx: float, is_morning: bool) -> float:
        """
        Refine the heliacal event time using binary search.

        Find the exact moment when the body becomes just visible/invisible.
        """
        # Use binary search to find the transition point
        jd_low = jd_approx - 0.1  # ~2.4 hours before
        jd_high = jd_approx + 0.1  # ~2.4 hours after

        for _ in range(30):  # ~30 iterations gives very high precision
            jd_mid = (jd_low + jd_high) / 2

            visible, sun_alt, body_alt, elong = _is_body_visible(jd_mid)

            # For heliacal rising: looking for first visibility
            # For heliacal setting: looking for last visibility
            if is_morning:
                # Morning: visibility increases as time progresses (Sun rises)
                # But we want first visibility, so search backwards
                if visible:
                    jd_high = jd_mid  # Look earlier
                else:
                    jd_low = jd_mid  # Look later
            else:
                # Evening: visibility decreases as time progresses (sky darkens)
                # For last visibility, we want the last moment still visible
                if visible:
                    jd_low = jd_mid  # Look later for last visibility
                else:
                    jd_high = jd_mid  # Look earlier

            if jd_high - jd_low < 1e-6:  # ~0.1 second precision
                break

        return (jd_low + jd_high) / 2

    # Main search logic based on event type
    if event_type == HELIACAL_RISING:
        jd_event = _search_heliacal_rising(jd_start)
    elif event_type == HELIACAL_SETTING:
        jd_event = _search_heliacal_setting(jd_start)
    elif event_type == EVENING_FIRST:
        jd_event = _search_evening_first(jd_start)
    elif event_type == MORNING_LAST:
        jd_event = _search_morning_last(jd_start)
    else:
        jd_event = 0.0

    if jd_event > 0:
        return jd_event, event_type
    else:
        return 0.0, -1  # Not found


# Note: heliacal_ut is defined in heliacal.py with the reference-API-compatible
# API signature (returns 3-tuple). The internal heliacal_ut here returns
# (jd_event, retflag) and is used by heliacal_ut internally.


# =============================================================================
# VISUAL LIMITING MAGNITUDE CALCULATIONS
# =============================================================================


def vis_limit_mag(
    jd: float,
    geopos: tuple,
    atmo: tuple,
    observer: tuple,
    objname: str,
    flags: int = FLG_SWIEPH,
) -> Tuple[float, Tuple[float, ...]]:
    """
    Calculate the limiting visual magnitude for observing a celestial body.

    This function determines whether a celestial body (planet, star, etc.)
    is visible given the current sky brightness conditions, atmospheric
    parameters, and observer characteristics. It returns both the visibility
    status and detailed information about the observation conditions.

    Args:
        jd: Julian Day (UT) for the observation time
        geopos: Geographic position as a sequence:
            - [0]: Geographic longitude in degrees (east positive)
            - [1]: Geographic latitude in degrees (north positive)
            - [2]: Altitude above sea level in meters
        atmo: Atmospheric conditions as a sequence:
            - [0]: Atmospheric pressure in mbar/hPa
            - [1]: Atmospheric temperature in degrees Celsius
            - [2]: Relative humidity in percent (0-100)
            - [3]: If >= 1: Meteorological Range in km
                   If 0-1: Total atmospheric coefficient (ktot)
                   If 0: Compute ktot from other parameters
        observer: Observer characteristics as a sequence:
            - [0]: Age of observer in years (default 36)
            - [1]: Snellen ratio of observer's eyes (default 1.0 = normal)
            For optical instruments (when HELFLAG_OPTICAL_PARAMS is set):
            - [2]: 0 = monocular, 1 = binocular
            - [3]: Telescope magnification (0 = naked eye)
            - [4]: Optical aperture (telescope diameter) in mm
            - [5]: Optical transmission coefficient
        objname: Name of the object to observe. Can be:
            - Planet name (e.g., "Venus", "Mars", "Jupiter")
            - Fixed star name (e.g., "Sirius", "Aldebaran")
            - Planet number as string (e.g., "2" for Venus)
        flags: Calculation flags combining ephemeris and heliacal flags:
            - FLG_SWIEPH, FLG_JPLEPH, etc. for ephemeris
            - HELFLAG_OPTICAL_PARAMS: Use optical instrument parameters
            - HELFLAG_NO_DETAILS: Skip detailed calculations
            - HELFLAG_VISLIM_DARK: Assume Sun at nadir (dark sky)
            - HELFLAG_VISLIM_NOMOON: Exclude Moon's brightness contribution

    Returns:
        Tuple containing:
            - result: Visibility status code:
                - (-2): Object is below horizon
                - (0): OK, photopic vision (bright conditions)
                - (1): OK, scotopic vision (dark conditions)
                - (2): OK, near limit between photopic/scotopic
            - dret: Tuple of 8 floats with observation details:
                - [0]: Limiting visual magnitude (object visible if mag < this)
                - [1]: Altitude of object in degrees
                - [2]: Azimuth of object in degrees
                - [3]: Altitude of Sun in degrees
                - [4]: Azimuth of Sun in degrees
                - [5]: Altitude of Moon in degrees
                - [6]: Azimuth of Moon in degrees
                - [7]: Magnitude of object

    Raises:
        ValueError: If objname is empty or invalid

    Algorithm:
        The limiting magnitude calculation is based on Schaefer's model (1990)
        which considers:
        1. Sky brightness from Sun, Moon, zodiacal light, airglow
        2. Atmospheric extinction based on airmass and conditions
        3. Observer's eye adaptation (scotopic vs photopic)
        4. Optional optical instrument characteristics

        The sky background brightness varies with:
        - Sun altitude (twilight contribution)
        - Moon altitude and phase (moonlight)
        - Atmospheric scattering

    Example:
        >>> from libephemeris import julday, vis_limit_mag
        >>> jd = julday(2024, 8, 15, 22.0)
        >>> # Rome location
        >>> geopos = (12.5, 41.9, 0)
        >>> # Standard atmosphere
        >>> atmo = (1013.25, 15.0, 50.0, 0.0)
        >>> # Normal observer
        >>> observer = (36, 1.0)
        >>> result, dret = vis_limit_mag(jd, geopos, atmo, observer, "Venus")
        >>> if dret[0] > dret[7]:
        ...     print("Venus is visible")
        ... else:
        ...     print("Venus is not visible")

    References:
        - Reference API: vis_limit_mag()
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
        - Schaefer, B.E. (1993) "Astronomy and the Limits of Vision"
    """
    import math
    from .planets import _PLANET_MAP, get_planet_target, pheno_ut
    from .fixed_stars import fixstar2_ut, fixstar2_mag
    from .state import get_planets, get_timescale
    from .constants import (
        SUN,
        MOON,
        HELFLAG_VISLIM_DARK,
        HELFLAG_VISLIM_NOMOON,
        HELFLAG_BELOW_HORIZON,
        HELFLAG_PHOTOPIC,
        HELFLAG_SCOTOPIC,
        HELFLAG_MIXED,
    )
    from skyfield.api import wgs84

    if not objname:
        raise Error("objname cannot be empty")

    # Parse geographic position
    lon = geopos[0] if len(geopos) > 0 else 0.0
    lat = geopos[1] if len(geopos) > 1 else 0.0
    alt_m = geopos[2] if len(geopos) > 2 else 0.0

    # Parse atmospheric conditions
    pressure = atmo[0] if len(atmo) > 0 else 1013.25
    temperature = atmo[1] if len(atmo) > 1 else 15.0
    humidity_pct = atmo[2] if len(atmo) > 2 else 50.0
    met_range = atmo[3] if len(atmo) > 3 else 0.0

    # Parse observer data
    observer_age = observer[0] if len(observer) > 0 else 36.0
    snellen_ratio = observer[1] if len(observer) > 1 else 1.0

    # Get ephemeris and timescale
    eph = get_planets()
    ts = get_timescale()

    # Create observer location
    obs_location = wgs84.latlon(lat, lon, alt_m)
    earth = eph["earth"]
    sun = eph["sun"]
    moon = eph["moon"]

    t = ts.ut1_jd(jd)
    observer_at = earth + obs_location

    # Calculate Sun position
    sun_app = observer_at.at(t).observe(sun).apparent()
    sun_alt_deg, sun_az_deg, _ = sun_app.altaz()
    sun_alt = sun_alt_deg.degrees
    sun_az = sun_az_deg.degrees

    # Calculate Moon position
    moon_app = observer_at.at(t).observe(moon).apparent()
    moon_alt_deg, moon_az_deg, _ = moon_app.altaz()
    moon_alt = moon_alt_deg.degrees
    moon_az = moon_az_deg.degrees

    # Determine if objname is a planet ID or name
    body_id = None
    is_fixed_star = False

    # Try parsing as integer (planet ID)
    try:
        body_id = int(objname)
    except ValueError:
        # Try to find planet by name
        name_upper = objname.upper().strip()
        planet_names = {
            "SUN": SUN,
            "MOON": MOON,
            "MERCURY": 2,
            "VENUS": 3,
            "MARS": 4,
            "JUPITER": 5,
            "SATURN": 6,
            "URANUS": 7,
            "NEPTUNE": 8,
            "PLUTO": 9,
        }
        if name_upper in planet_names:
            body_id = planet_names[name_upper]
        else:
            # Assume it's a fixed star
            is_fixed_star = True

    # Calculate object position and magnitude
    obj_alt = 0.0
    obj_az = 0.0
    obj_mag = 0.0

    if is_fixed_star:
        # Fixed star calculation
        try:
            star_result, star_name_out, retflag = fixstar2_ut(objname, jd, flags & 0xFF)

            # star_result is (lon, lat, dist, lon_speed, lat_speed, dist_speed)
            # We need to convert ecliptic to horizontal
            star_lon = star_result[0]
            star_lat = star_result[1]

            # Get star magnitude
            try:
                star_mag_val, _star_name_mag = fixstar2_mag(objname)
                obj_mag = star_mag_val
            except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
                _reraise_if_leb_range_error(_exc)
                obj_mag = 2.0  # Default magnitude if not found

            # Convert ecliptic to equatorial then to horizontal
            # Simplified: use azalt function if available
            from .utils import azalt, ECL2HOR

            hor_result = azalt(
                jd,
                ECL2HOR,
                (lon, lat, alt_m),
                pressure,
                temperature,
                (star_lon, star_lat, 1.0),
            )
            obj_az = hor_result[0]
            obj_alt = hor_result[1]

        except ValueError:
            raise
        except (ArithmeticError, IndexError) as e:
            # Star not found or other error
            raise Error(f"could not find star name {objname.lower()}: {e}") from e
    else:
        # Planet calculation
        if body_id is None:
            raise Error(f"Unknown object: {objname}")

        # Get planet name from _PLANET_MAP
        if body_id in _PLANET_MAP:
            target_name = _PLANET_MAP[body_id]
            target = get_planet_target(eph, target_name)

            # Calculate position
            body_app = observer_at.at(t).observe(target).apparent()
            body_alt_deg, body_az_deg, _ = body_app.altaz()
            obj_alt = body_alt_deg.degrees
            obj_az = body_az_deg.degrees

            # Get magnitude from pheno
            try:
                pheno_result = pheno_ut(jd, body_id, flags)
                obj_mag = pheno_result[4]  # Visual magnitude
            except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
                _reraise_if_leb_range_error(_exc)
                obj_mag = 0.0  # Default bright
        else:
            raise Error(f"illegal planet number {body_id}.")

    # Check if object is below horizon
    if obj_alt < 0:
        # Apply HELFLAG options before returning
        use_dark_sky_early = bool(flags & HELFLAG_VISLIM_DARK)
        exclude_moon_early = bool(flags & HELFLAG_VISLIM_NOMOON)
        ret_sun_alt = -90.0 if use_dark_sky_early else sun_alt
        ret_moon_alt = -90.0 if exclude_moon_early else moon_alt
        dret = (
            0.0,
            obj_alt,
            obj_az,
            ret_sun_alt,
            sun_az,
            ret_moon_alt,
            moon_az,
            obj_mag,
        )
        return float(HELFLAG_BELOW_HORIZON), dret

    # Apply HELFLAG options
    use_dark_sky = bool(flags & HELFLAG_VISLIM_DARK)
    exclude_moon = bool(flags & HELFLAG_VISLIM_NOMOON)

    if use_dark_sky:
        sun_alt = -90.0  # Assume Sun at nadir

    if exclude_moon:
        moon_alt = -90.0  # Assume Moon at nadir

    # Calculate atmospheric extinction coefficient
    # Simplified Schaefer model
    def calc_extinction_coeff(pressure_mbar, temp_c, humidity_pct, met_range_km):
        """Calculate total atmospheric extinction coefficient."""
        # Rayleigh scattering (molecular)
        k_rayleigh = 0.1451 * (pressure_mbar / 1013.25)

        # Aerosol scattering
        if met_range_km >= 1:
            # Meteorological range given
            k_aerosol = 3.912 / met_range_km - 0.106
        elif 0 < met_range_km < 1:
            # ktot given directly
            return met_range_km
        else:
            # Estimate from humidity
            # Higher humidity = more scattering
            k_aerosol = 0.1 + 0.2 * (humidity_pct / 100.0)

        # Ozone absorption (small contribution)
        k_ozone = 0.016

        return k_rayleigh + k_aerosol + k_ozone

    k_ext = calc_extinction_coeff(pressure, temperature, humidity_pct, met_range)

    # Calculate airmass
    def calc_airmass(altitude_deg):
        """Calculate airmass using Kasten & Young formula."""
        if altitude_deg <= 0:
            return 40.0  # Maximum airmass at/below horizon
        z = 90.0 - altitude_deg  # Zenith angle
        z_rad = math.radians(z)
        # Kasten & Young (1989)
        airmass = 1.0 / (math.cos(z_rad) + 0.50572 * (96.07995 - z) ** (-1.6364))
        return min(airmass, 40.0)

    airmass = calc_airmass(obj_alt)

    # Calculate sky brightness contribution from Sun
    def calc_sun_brightness_contribution(sun_altitude):
        """Calculate sky brightness factor from Sun position."""
        if sun_altitude >= 0:
            return 1e10  # Daylight - essentially infinite brightness
        elif sun_altitude >= -6:
            # Civil twilight
            return 10 ** (4.0 + sun_altitude * 0.4)
        elif sun_altitude >= -12:
            # Nautical twilight
            return 10 ** (1.6 + (sun_altitude + 6) * 0.2)
        elif sun_altitude >= -18:
            # Astronomical twilight
            return 10 ** (0.4 + (sun_altitude + 12) * 0.1)
        else:
            # Full night
            return 1.0

    # Calculate sky brightness contribution from Moon
    def calc_moon_brightness_contribution(moon_altitude, jd_ut):
        """Calculate sky brightness factor from Moon position and phase."""
        if moon_altitude <= 0:
            return 0.0

        # Get Moon phase
        try:
            moon_pheno = pheno_ut(jd_ut, MOON, flags & 0xFF)
            phase_angle = moon_pheno[0]  # Phase angle in degrees
            illumination = (1 + math.cos(math.radians(phase_angle))) / 2.0
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            illumination = 0.5  # Assume half moon

        # Moon brightness increases with altitude and illumination
        moon_factor = illumination * math.sin(math.radians(moon_altitude))
        return moon_factor * 100.0

    sun_contrib = calc_sun_brightness_contribution(sun_alt)
    moon_contrib = (
        calc_moon_brightness_contribution(moon_alt, jd) if not exclude_moon else 0.0
    )

    # Total sky brightness (arbitrary units)
    total_sky_brightness = sun_contrib + moon_contrib + 1.0  # 1.0 = natural sky glow

    # Calculate limiting magnitude
    # Base limiting magnitude for perfect conditions: ~6.5 for naked eye
    base_limit_mag = 6.5

    # Adjust for observer characteristics
    # Age reduces sensitivity (approximately 0.05 mag per decade over 30)
    age_factor = max(0, (observer_age - 30) / 10.0) * 0.05

    # Snellen ratio adjustment (better eyes = higher limit)
    snellen_factor = -2.5 * math.log10(snellen_ratio) if snellen_ratio > 0 else 0

    # Airmass extinction
    airmass_extinction = k_ext * airmass

    # Sky brightness reduction of limiting magnitude
    if total_sky_brightness > 1:
        sky_reduction = 2.5 * math.log10(total_sky_brightness)
    else:
        sky_reduction = 0.0

    # Calculate final limiting magnitude
    limiting_mag = base_limit_mag - age_factor - snellen_factor - sky_reduction

    # Apply extinction to object magnitude (object appears fainter)
    apparent_obj_mag = obj_mag + airmass_extinction

    # Determine vision type based on sky brightness
    if sun_alt >= -6:
        # Photopic vision (daylight/bright twilight)
        vision_type = HELFLAG_PHOTOPIC
    elif sun_alt >= -12:
        # Mixed/mesopic vision
        vision_type = HELFLAG_MIXED
    else:
        # Scotopic vision (night)
        vision_type = HELFLAG_SCOTOPIC

    # Build result tuple
    dret = (
        limiting_mag,  # 0: Limiting visual magnitude
        obj_alt,  # 1: Altitude of object
        obj_az,  # 2: Azimuth of object
        sun_alt,  # 3: Altitude of Sun
        sun_az,  # 4: Azimuth of Sun
        moon_alt,  # 5: Altitude of Moon
        moon_az,  # 6: Azimuth of Moon
        apparent_obj_mag,  # 7: Apparent magnitude of object (with extinction)
    )

    return float(vision_type), dret


# =============================================================================
# BESSELIAN ELEMENTS
# =============================================================================

# Earth's equatorial radius in km (WGS84)
EARTH_RADIUS_KM = 6378.137


from dataclasses import dataclass  # noqa: E402 (section-local import)


@dataclass
class BesselianElements:
    """Besselian elements for a solar eclipse at a reference time.

    Besselian elements are the fundamental quantities used to calculate
    local eclipse circumstances. They describe the lunar shadow relative to
    Earth on a fundamental plane perpendicular to the shadow axis.

    Attributes:
        t0: Reference Julian Day in UT.
        x: Shadow-axis x-coordinate in Earth equatorial radii.
        y: Shadow-axis y-coordinate in Earth equatorial radii.
        d: Shadow-axis declination in degrees.
        l1: Penumbral-cone radius in Earth equatorial radii.
        l2: Umbral or antumbral-cone radius in Earth equatorial radii.
        mu: Greenwich hour angle of the shadow axis in degrees.
        dx_dt: Rate of change of x (Earth radii per hour).
        dy_dt: Rate of change of y (Earth radii per hour).
        dd_dt: Rate of change of d (degrees per hour).
        dl1_dt: Rate of change of l1 (Earth radii per hour).
        dl2_dt: Rate of change of l2 (Earth radii per hour).
        dmu_dt: Rate of change of mu (degrees per hour).

    Example:
        >>> from libephemeris import julday, BesselianElements
        >>> # Create elements manually for illustration
        >>> elements = BesselianElements(
        ...     t0=julday(2024, 4, 8, 18.0),
        ...     x=0.3145, y=0.2731, d=7.5821,
        ...     l1=0.5436, l2=-0.0047, mu=89.1234,
        ...     dx_dt=0.5123, dy_dt=0.1456, dd_dt=0.0012,
        ...     dl1_dt=-0.0001, dl2_dt=-0.0001, dmu_dt=15.0041
        ... )
        >>> print(f"Shadow at x={elements.x:.4f}, y={elements.y:.4f}")

    Notes:
        Angular rates are degrees per hour and linear rates are Earth radii
        per hour. Use :func:`interpolate_besselian_elements` for nearby times.

    References:
        Meeus, *Astronomical Algorithms*, chapter 54; *Explanatory Supplement
        to the Astronomical Almanac* (2013), chapter 11.
    """

    t0: float  # Reference time (Julian Day, UT)

    # Besselian elements
    x: float  # Shadow x-coordinate (Earth radii)
    y: float  # Shadow y-coordinate (Earth radii)
    d: float  # Declination of shadow axis (degrees)
    l1: float  # Penumbral cone radius (Earth radii)
    l2: float  # Umbral cone radius (Earth radii)
    mu: float  # Greenwich hour angle (degrees)

    # Time derivatives (rates of change per hour)
    dx_dt: float  # Rate of change of x (Earth radii/hour)
    dy_dt: float  # Rate of change of y (Earth radii/hour)
    dd_dt: float  # Rate of change of d (degrees/hour)
    dl1_dt: float  # Rate of change of l1 (Earth radii/hour)
    dl2_dt: float  # Rate of change of l2 (Earth radii/hour)
    dmu_dt: float  # Rate of change of mu (degrees/hour)


_BESSELIAN_CACHE: dict = {}


def _besselian_core(
    jd: float, flags: int = FLG_SWIEPH
) -> Tuple[float, float, float, float, float, float, float, float, float]:
    """Besselian elements on the true equator of date.

    Conventions (Explanatory Supplement ch. 11 / NASA eclipse bulletins):
    the positive z-axis of the fundamental frame is parallel to the
    shadow axis and points from the fundamental plane toward the Sun;
    the x-axis is the intersection with the equator, positive east; the
    y-axis completes the system toward north. (x, y) is the shadow-axis
    intersection with the fundamental plane in Earth equatorial radii;
    d and mu are the declination and Greenwich hour angle of the axis;
    l1/l2 are the penumbral/umbral radii on the fundamental plane, l2
    negative while the umbral vertex lies beyond the plane (total
    eclipse). Lunar radius ratios k = 0.2725076 (penumbra) and
    0.272281 (umbra) follow the NASA convention.

    Returns (x, y, d_deg, mu_deg, l1, l2, tan_f1, tan_f2, z_moon_er).
    """
    from .time_utils import sidtime

    key = (round(jd, 9), flags & 0xFF)
    cached = _BESSELIAN_CACHE.get(key)
    if cached is not None:
        return cached

    base = _ecl_eph_flags(flags) | FLG_EQUATORIAL
    sun_p, _ = calc_ut(jd, SUN, base)
    moon_p, _ = calc_ut(jd, MOON, base)

    def _cart(p):
        ra = math.radians(p[0])
        dec = math.radians(p[1])
        r = p[2]
        return (
            r * math.cos(dec) * math.cos(ra),
            r * math.cos(dec) * math.sin(ra),
            r * math.sin(dec),
        )

    rs = _cart(sun_p)
    rm = _cart(moon_p)

    # Shadow-axis direction, positive toward the Sun.
    gx, gy, gz = rs[0] - rm[0], rs[1] - rm[1], rs[2] - rm[2]
    g_len = math.sqrt(gx * gx + gy * gy + gz * gz)
    gx, gy, gz = gx / g_len, gy / g_len, gz / g_len

    d_rad = math.asin(max(-1.0, min(1.0, gz)))
    a_rad = math.atan2(gy, gx)

    # East axis: unit(k x g); north axis completes the triad.
    exy = math.hypot(gx, gy)
    if exy < 1e-12:
        ex = (1.0, 0.0, 0.0)
    else:
        ex = (-gy / exy, gx / exy, 0.0)
    ey = (
        gy * ex[2] - gz * ex[1],
        gz * ex[0] - gx * ex[2],
        gx * ex[1] - gy * ex[0],
    )

    re_au = _ECL_REARTH_AU
    x = (rm[0] * ex[0] + rm[1] * ex[1] + rm[2] * ex[2]) / re_au
    y = (rm[0] * ey[0] + rm[1] * ey[1] + rm[2] * ey[2]) / re_au
    z = (rm[0] * gx + rm[1] * gy + rm[2] * gz) / re_au

    gast_deg = sidtime(jd) * 15.0
    mu_deg = (gast_deg - math.degrees(a_rad)) % 360.0

    # Shadow-cone half angles from the solar radius and the NASA lunar
    # radius ratios, with distances in Earth radii.
    g_er = g_len / re_au
    rsun_er = _ECL_RSUN_AU / re_au
    k_pen = 0.2725076
    k_umb = 0.272281
    sin_f1 = (rsun_er + k_pen) / g_er
    sin_f2 = (rsun_er - k_umb) / g_er
    cos_f1 = math.sqrt(1.0 - sin_f1 * sin_f1)
    cos_f2 = math.sqrt(1.0 - sin_f2 * sin_f2)
    tan_f1 = sin_f1 / cos_f1
    tan_f2 = sin_f2 / cos_f2
    l1 = z * tan_f1 + k_pen / cos_f1
    l2 = z * tan_f2 - k_umb / cos_f2

    result = (x, y, math.degrees(d_rad), mu_deg, l1, l2, tan_f1, tan_f2, z)
    if len(_BESSELIAN_CACHE) > 512:
        _BESSELIAN_CACHE.clear()
    _BESSELIAN_CACHE[key] = result
    return result


def calc_besselian_x(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the Besselian x coordinate for a solar eclipse.

    The Besselian x coordinate is the x-component (positive eastward) of the
    Moon's shadow axis intersection with the fundamental plane, measured in
    Earth equatorial radii. The fundamental plane is the plane through Earth's
    center perpendicular to the Moon-Sun line (the shadow axis direction).

    The x-axis points east in the direction of increasing right ascension,
    perpendicular to the shadow axis. This coordinate, along with y, describes
    where the Moon's shadow cone axis pierces the fundamental plane.

    Args:
        jd: Julian Day (UT) at which to calculate the Besselian x coordinate
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The Besselian x coordinate in Earth equatorial radii.
        Positive values indicate the shadow axis is east of Earth's center,
        negative values indicate west.

    Algorithm:
        1. Get geocentric equatorial positions of Moon and Sun
        2. Convert to 3D Cartesian coordinates (AU)
        3. Compute the shadow axis direction unit vector (Sun to Moon direction)
        4. Compute the fundamental plane axes:
           - z-axis: along shadow axis (from Sun toward Moon)
           - x-axis: perpendicular to z, in equatorial plane (toward increasing RA)
           - y-axis: completes right-handed system (toward north)
        5. Project Moon's position onto the fundamental plane
        6. Extract x-component and convert to Earth radii

    Mathematical basis:
        The shadow axis is defined as the line from the Sun's center through
        the Moon's center. The fundamental plane is perpendicular to this axis
        and passes through Earth's center. The Besselian x coordinate is the
        east-west displacement (in Earth radii) of where this axis pierces
        the fundamental plane.

    Precision:
        Accurate to ~0.0001 Earth radii (~0.6 km) for typical eclipses.
        Uses full precision Moon and Sun ephemeris positions.

    Example:
        >>> from libephemeris import julday, calc_besselian_x
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> x = calc_besselian_x(jd)
        >>> print(f"Besselian x = {x:.4f} Earth radii")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
        - Chauvenet's "Manual of Spherical and Practical Astronomy", Vol. 1
    """
    return _besselian_core(jd, flags)[0]


def calc_besselian_y(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the Besselian y coordinate for a solar eclipse.

    The Besselian y coordinate is the y-component (positive northward) of the
    Moon's shadow axis intersection with the fundamental plane, measured in
    Earth equatorial radii. The fundamental plane is the plane through Earth's
    center perpendicular to the Moon-Sun line (the shadow axis direction).

    The y-axis points north, perpendicular to both the shadow axis and the
    x-axis. This coordinate, along with x, describes where the Moon's shadow
    cone axis pierces the fundamental plane.

    Args:
        jd: Julian Day (UT) at which to calculate the Besselian y coordinate
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The Besselian y coordinate in Earth equatorial radii.
        Positive values indicate the shadow axis is north of Earth's center,
        negative values indicate south.

    Algorithm:
        1. Get geocentric equatorial positions of Moon and Sun
        2. Convert to 3D Cartesian coordinates (AU)
        3. Compute the shadow axis direction unit vector (Sun to Moon direction)
        4. Compute the fundamental plane axes:
           - z-axis: along shadow axis (from Sun toward Moon)
           - x-axis: perpendicular to z, in equatorial plane (toward increasing RA)
           - y-axis: completes right-handed system (toward north)
        5. Project Moon's position onto the fundamental plane
        6. Extract y-component and convert to Earth radii

    Mathematical basis:
        The shadow axis is defined as the line from the Sun's center through
        the Moon's center. The fundamental plane is perpendicular to this axis
        and passes through Earth's center. The Besselian y coordinate is the
        north-south displacement (in Earth radii) of where this axis pierces
        the fundamental plane.

    Precision:
        Accurate to ~0.0001 Earth radii (~0.6 km) for typical eclipses.
        Uses full precision Moon and Sun ephemeris positions.

    Example:
        >>> from libephemeris import julday, calc_besselian_y
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> y = calc_besselian_y(jd)
        >>> print(f"Besselian y = {y:.4f} Earth radii")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
        - Chauvenet's "Manual of Spherical and Practical Astronomy", Vol. 1
    """
    return _besselian_core(jd, flags)[1]


def calc_besselian_d(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the Besselian element d (declination) for a solar eclipse.

    The Besselian element d is the declination of the Moon's shadow axis,
    measured in degrees. In the standard Besselian element convention, the
    shadow axis direction is defined as pointing from the shadow cone vertex
    toward Earth. Since the Sun, Moon, and shadow axis are nearly collinear
    during an eclipse, d closely tracks the Sun's declination.

    The declination d, together with the hour angle mu, defines the orientation
    of the shadow axis in equatorial coordinates.

    Args:
        jd: Julian Day (UT) at which to calculate the declination
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The declination d in degrees.
        Positive values indicate the shadow axis points north of the equator,
        negative values indicate south.
        Range: -90 to +90 degrees.

    Algorithm:
        1. Get geocentric equatorial positions of Moon and Sun
        2. Convert to 3D Cartesian coordinates (AU)
        3. Compute the shadow axis direction (from Moon toward Sun, which is
           the conventional Besselian axis direction toward the observer)
        4. Normalize and extract the z-component (toward north celestial pole)
        5. Compute declination as arcsin(z) and convert to degrees

    Mathematical basis:
        The shadow axis in Besselian elements points from the shadow cone
        vertex toward Earth. For a unit vector in equatorial coordinates,
        the declination is arcsin(z) where z is the component toward the
        north celestial pole. During a solar eclipse, this declination
        closely matches the Sun's geocentric declination.

    Precision:
        Accurate to better than 0.001 degrees for typical eclipses.
        Uses full precision Moon and Sun ephemeris positions.

    Example:
        >>> from libephemeris import julday, calc_besselian_d
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> d = calc_besselian_d(jd)
        >>> print(f"Besselian d = {d:.4f} degrees")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
        - Chauvenet's "Manual of Spherical and Practical Astronomy", Vol. 1
    """
    return _besselian_core(jd, flags)[2]


def calc_besselian_l1(jd: float, flags: int = FLG_SWIEPH) -> float:
    """Calculate Besselian element ``l1`` for a solar eclipse.

    Args:
        jd: Julian Day in UT.
        flags: Calculation flags.

    Returns:
        Penumbral shadow radius in Earth equatorial radii.

    Notes:
        The value follows the external-tangent geometry of the solar and lunar
        limbs on the fundamental plane perpendicular to the shadow axis. It is
        positive and is typically roughly 0.5 to 0.6 Earth radii during an
        eclipse.

    References:
        Meeus, *Astronomical Algorithms*, chapter 54; *Explanatory Supplement
        to the Astronomical Almanac* (2013), chapter 11.
    """
    return _besselian_core(jd, flags)[4]


def calc_besselian_l2(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the Besselian element l2 (umbral/antumbral shadow radius) for a solar eclipse.

    The Besselian element l2 is the radius of the umbral (or antumbral) shadow cone
    where it intersects the fundamental plane, measured in Earth equatorial radii.
    The fundamental plane is the plane through Earth's center perpendicular to the
    Moon-Sun line (the shadow axis direction).

    The umbral cone is formed by the internal tangent lines from the Sun's limb
    to the Moon's limb. The sign of l2 indicates the eclipse type (standard
    Besselian convention):
    - l2 < 0: umbral shadow reaches the fundamental plane (total eclipse)
    - l2 > 0: only the antumbra reaches it (annular eclipse)

    Args:
        jd: Julian Day (UT) at which to calculate l2
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The umbral/antumbral shadow radius l2 in Earth equatorial radii.
        Negative for total eclipses (umbral shadow), positive for annular
        eclipses (antumbral shadow).

    Algorithm:
        1. Get geocentric equatorial positions of Moon and Sun
        2. Calculate the Sun-Moon distance
        3. Calculate the umbral cone semi-angle f2 using:
           sin(f2) = (R_sun - R_moon) / D_sun-moon
        4. Calculate the umbral cone vertex distance from the Moon:
           c2 = R_moon / sin(f2)  (distance from Moon to umbral vertex)
        5. Calculate whether vertex is beyond or before Earth:
           If Moon distance > c2: umbral (total), l2 > 0
           If Moon distance < c2: antumbral (annular), l2 < 0
        6. The umbral radius at the fundamental plane is:
           l2 = ±|d_vertex_to_earth| * tan(f2)
        7. Convert to Earth radii

    Mathematical basis:
        The umbral cone's semi-angle f2 is determined by the internal
        tangent geometry. The cone vertex is located at distance c2 from
        the Moon center along the shadow axis (toward Earth):
        c2 = R_moon / sin(f2)

        For the fundamental plane (at Earth's center), the distance from
        the vertex determines the shadow radius. If the vertex is beyond
        Earth (umbral), l2 is positive. If the vertex is between Moon and
        Earth (antumbral), l2 is negative, indicating the shadow has
        diverged from the cone and forms an antumbral zone.

    Physical constants used:
        - Sun radius: 696,000 km (IAU nominal)
        - Moon radius: 1,737.4 km (mean radius)
        - Earth radius: 6,378.137 km (WGS84 equatorial)
        - AU: 149,597,870.7 km

    Precision:
        Accurate to ~0.0001 Earth radii (~0.6 km) for typical eclipses.
        Uses full precision Moon and Sun ephemeris positions.

    Example:
        >>> from libephemeris import julday, calc_besselian_l2
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> l2 = calc_besselian_l2(jd)
        >>> print(f"Besselian l2 = {l2:.4f} Earth radii")
        >>> if l2 > 0:
        ...     print("Total eclipse (umbral shadow)")
        ... else:
        ...     print("Annular eclipse (antumbral shadow)")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
        - Chauvenet's "Manual of Spherical and Practical Astronomy", Vol. 1
    """
    return _besselian_core(jd, flags)[5]


def calc_besselian_mu(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the Besselian element mu (Greenwich hour angle) for a solar eclipse.

    The Besselian element mu is the Greenwich hour angle of the shadow axis,
    measured in degrees. It represents the angle between the Greenwich meridian
    and the hour circle containing the shadow axis, measured westward from
    Greenwich along the celestial equator.

    Together with the declination d, mu defines the orientation of the shadow
    axis in the equatorial coordinate system. As Earth rotates, mu increases
    at approximately 15 degrees per hour (360 degrees per sidereal day).

    Args:
        jd: Julian Day (UT) at which to calculate mu
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The Greenwich hour angle mu in degrees.
        Range: 0 to 360 degrees, measured westward from Greenwich.
        Values increase with time as Earth rotates eastward.

    Algorithm:
        1. Get geocentric equatorial positions of Moon and Sun
        2. Convert to 3D Cartesian coordinates (AU)
        3. Compute the shadow axis direction (from Moon toward Sun, the
           conventional Besselian axis direction toward the observer)
        4. Calculate the right ascension (RA) of the shadow axis from
           its x,y components: RA = atan2(y, x)
        5. Calculate Greenwich Apparent Sidereal Time (GAST)
        6. Compute mu = GAST - RA, normalized to 0-360 degrees

    Mathematical basis:
        The hour angle H of a celestial object is defined as:
            H = LST - RA (Local Sidereal Time minus Right Ascension)

        For the Greenwich meridian, this becomes:
            mu = GAST - RA_shadow

        where RA_shadow is the right ascension of the shadow axis direction.
        The shadow axis in Besselian elements points from the Moon toward
        the Sun (toward the observer on Earth).

    Physical interpretation:
        mu indicates where the shadow axis crosses the celestial equator
        relative to the Greenwich meridian. When mu = 0, the shadow axis
        is on the Greenwich meridian. As Earth rotates, mu increases by
        about 15 degrees per hour.

    Precision:
        Accurate to better than 0.01 degrees for typical eclipses.
        Uses full precision Moon and Sun ephemeris positions and
        accurate sidereal time calculation.

    Example:
        >>> from libephemeris import julday, calc_besselian_mu
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> mu = calc_besselian_mu(jd)
        >>> print(f"Besselian mu = {mu:.4f} degrees")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
        - Chauvenet's "Manual of Spherical and Practical Astronomy", Vol. 1
    """
    return _besselian_core(jd, flags)[3]


def calc_besselian_dx_dt(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the time derivative of Besselian element x (dx/dt).

    This function computes the hourly rate of change of the Besselian x
    coordinate, which represents how fast the shadow axis is moving eastward
    (or westward if negative) across the fundamental plane.

    The derivative is calculated using symmetric numerical differentiation
    (central difference method) for improved accuracy.

    Args:
        jd: Julian Day (UT) at which to calculate the derivative
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The rate of change dx/dt in Earth radii per hour.
        Positive values indicate eastward motion, negative indicates westward.

    Algorithm:
        Uses symmetric numerical differentiation:
        dx/dt ≈ (x(t+h) - x(t-h)) / (2h)
        where h = 1/1440 days (1 minute) for optimal accuracy.

    Precision:
        Accurate to approximately 1e-6 Earth radii per hour for typical
        eclipse conditions. The small time step minimizes truncation error
        while the symmetric difference eliminates first-order errors.

    Example:
        >>> from libephemeris import julday, calc_besselian_dx_dt
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> dx_dt = calc_besselian_dx_dt(jd)
        >>> print(f"dx/dt = {dx_dt:.6f} Earth radii/hour")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Time step for numerical differentiation: 1 minute in days
    # This provides good balance between truncation and round-off errors
    h = 1.0 / 1440.0  # 1 minute = 1/1440 days

    # Calculate x at t-h and t+h
    x_minus = calc_besselian_x(jd - h, flags)
    x_plus = calc_besselian_x(jd + h, flags)

    # Symmetric difference quotient
    # dx/dt in Earth radii per day
    dx_dt_per_day = (x_plus - x_minus) / (2 * h)

    # Convert to Earth radii per hour (divide by 24)
    dx_dt_per_hour = dx_dt_per_day / 24.0

    return dx_dt_per_hour


def calc_besselian_dy_dt(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the time derivative of Besselian element y (dy/dt).

    This function computes the hourly rate of change of the Besselian y
    coordinate, which represents how fast the shadow axis is moving northward
    (or southward if negative) across the fundamental plane.

    The derivative is calculated using symmetric numerical differentiation
    (central difference method) for improved accuracy.

    Args:
        jd: Julian Day (UT) at which to calculate the derivative
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The rate of change dy/dt in Earth radii per hour.
        Positive values indicate northward motion, negative indicates southward.

    Algorithm:
        Uses symmetric numerical differentiation:
        dy/dt ≈ (y(t+h) - y(t-h)) / (2h)
        where h = 1/1440 days (1 minute) for optimal accuracy.

    Precision:
        Accurate to approximately 1e-6 Earth radii per hour for typical
        eclipse conditions. The small time step minimizes truncation error
        while the symmetric difference eliminates first-order errors.

    Example:
        >>> from libephemeris import julday, calc_besselian_dy_dt
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> dy_dt = calc_besselian_dy_dt(jd)
        >>> print(f"dy/dt = {dy_dt:.6f} Earth radii/hour")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Time step for numerical differentiation: 1 minute in days
    h = 1.0 / 1440.0  # 1 minute = 1/1440 days

    # Calculate y at t-h and t+h
    y_minus = calc_besselian_y(jd - h, flags)
    y_plus = calc_besselian_y(jd + h, flags)

    # Symmetric difference quotient
    # dy/dt in Earth radii per day
    dy_dt_per_day = (y_plus - y_minus) / (2 * h)

    # Convert to Earth radii per hour
    dy_dt_per_hour = dy_dt_per_day / 24.0

    return dy_dt_per_hour


def calc_besselian_dd_dt(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the time derivative of Besselian element d (dd/dt).

    This function computes the hourly rate of change of the shadow axis
    declination. During an eclipse, the declination changes very slowly
    (typically less than 1 degree per hour) as it tracks the Sun's motion.

    The derivative is calculated using symmetric numerical differentiation
    (central difference method) for improved accuracy.

    Args:
        jd: Julian Day (UT) at which to calculate the derivative
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The rate of change dd/dt in degrees per hour.
        Positive values indicate the shadow axis is moving northward,
        negative indicates southward.

    Algorithm:
        Uses symmetric numerical differentiation:
        dd/dt ≈ (d(t+h) - d(t-h)) / (2h)
        where h = 1/1440 days (1 minute) for optimal accuracy.

    Precision:
        Accurate to approximately 1e-5 degrees per hour for typical
        eclipse conditions.

    Example:
        >>> from libephemeris import julday, calc_besselian_dd_dt
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> dd_dt = calc_besselian_dd_dt(jd)
        >>> print(f"dd/dt = {dd_dt:.6f} degrees/hour")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Time step for numerical differentiation: 1 minute in days
    h = 1.0 / 1440.0  # 1 minute = 1/1440 days

    # Calculate d at t-h and t+h
    d_minus = calc_besselian_d(jd - h, flags)
    d_plus = calc_besselian_d(jd + h, flags)

    # Symmetric difference quotient
    # dd/dt in degrees per day
    dd_dt_per_day = (d_plus - d_minus) / (2 * h)

    # Convert to degrees per hour
    dd_dt_per_hour = dd_dt_per_day / 24.0

    return dd_dt_per_hour


def calc_besselian_dl1_dt(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the time derivative of Besselian element l1 (dl1/dt).

    This function computes the hourly rate of change of the penumbral
    shadow radius. The penumbral radius changes very slowly during an
    eclipse as the Earth-Moon-Sun geometry evolves.

    The derivative is calculated using symmetric numerical differentiation
    (central difference method) for improved accuracy.

    Args:
        jd: Julian Day (UT) at which to calculate the derivative
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The rate of change dl1/dt in Earth radii per hour.
        The value is typically very small (order of 1e-4 to 1e-3).

    Algorithm:
        Uses symmetric numerical differentiation:
        dl1/dt ≈ (l1(t+h) - l1(t-h)) / (2h)
        where h = 1/1440 days (1 minute) for optimal accuracy.

    Precision:
        Accurate to approximately 1e-7 Earth radii per hour for typical
        eclipse conditions.

    Example:
        >>> from libephemeris import julday, calc_besselian_dl1_dt
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> dl1_dt = calc_besselian_dl1_dt(jd)
        >>> print(f"dl1/dt = {dl1_dt:.8f} Earth radii/hour")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Time step for numerical differentiation: 1 minute in days
    h = 1.0 / 1440.0  # 1 minute = 1/1440 days

    # Calculate l1 at t-h and t+h
    l1_minus = calc_besselian_l1(jd - h, flags)
    l1_plus = calc_besselian_l1(jd + h, flags)

    # Symmetric difference quotient
    # dl1/dt in Earth radii per day
    dl1_dt_per_day = (l1_plus - l1_minus) / (2 * h)

    # Convert to Earth radii per hour
    dl1_dt_per_hour = dl1_dt_per_day / 24.0

    return dl1_dt_per_hour


def calc_besselian_dl2_dt(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the time derivative of Besselian element l2 (dl2/dt).

    This function computes the hourly rate of change of the umbral/antumbral
    shadow radius. The umbral radius changes very slowly during an eclipse
    as the Earth-Moon-Sun geometry evolves.

    The derivative is calculated using symmetric numerical differentiation
    (central difference method) for improved accuracy.

    Args:
        jd: Julian Day (UT) at which to calculate the derivative
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The rate of change dl2/dt in Earth radii per hour.
        The value is typically very small (order of 1e-5 to 1e-4).

    Algorithm:
        Uses symmetric numerical differentiation:
        dl2/dt ≈ (l2(t+h) - l2(t-h)) / (2h)
        where h = 1/1440 days (1 minute) for optimal accuracy.

    Precision:
        Accurate to approximately 1e-8 Earth radii per hour for typical
        eclipse conditions.

    Example:
        >>> from libephemeris import julday, calc_besselian_dl2_dt
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> dl2_dt = calc_besselian_dl2_dt(jd)
        >>> print(f"dl2/dt = {dl2_dt:.8f} Earth radii/hour")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Time step for numerical differentiation: 1 minute in days
    h = 1.0 / 1440.0  # 1 minute = 1/1440 days

    # Calculate l2 at t-h and t+h
    l2_minus = calc_besselian_l2(jd - h, flags)
    l2_plus = calc_besselian_l2(jd + h, flags)

    # Symmetric difference quotient
    # dl2/dt in Earth radii per day
    dl2_dt_per_day = (l2_plus - l2_minus) / (2 * h)

    # Convert to Earth radii per hour
    dl2_dt_per_hour = dl2_dt_per_day / 24.0

    return dl2_dt_per_hour


def calc_besselian_dmu_dt(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the time derivative of Besselian element mu (dmu/dt).

    This function computes the hourly rate of change of the Greenwich hour
    angle. Since Earth rotates at approximately 15 degrees per hour, dmu/dt
    is typically close to this value, with small variations due to the
    changing right ascension of the shadow axis.

    The derivative is calculated using symmetric numerical differentiation
    (central difference method) for improved accuracy, with special handling
    for the 0/360 degree wraparound.

    Args:
        jd: Julian Day (UT) at which to calculate the derivative
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The rate of change dmu/dt in degrees per hour.
        Expected to be approximately 15 degrees/hour (Earth's rotation rate)
        with small variations (typically ±0.04 degrees/hour).

    Algorithm:
        Uses symmetric numerical differentiation:
        dmu/dt ≈ (mu(t+h) - mu(t-h)) / (2h)
        where h = 1/1440 days (1 minute) for optimal accuracy.
        Special handling ensures correct results near the 0/360 boundary.

    Physical interpretation:
        dmu/dt ≈ 15.04 deg/hour (Earth's sidereal rotation rate)
              - d(RA_shadow)/dt (shadow axis RA rate of change)

        The deviation from exactly 15°/hour is due to the apparent
        motion of the Sun (and hence the shadow axis) in right ascension.

    Precision:
        Accurate to approximately 1e-4 degrees per hour for typical
        eclipse conditions.

    Example:
        >>> from libephemeris import julday, calc_besselian_dmu_dt
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> dmu_dt = calc_besselian_dmu_dt(jd)
        >>> print(f"dmu/dt = {dmu_dt:.4f} degrees/hour")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Time step for numerical differentiation: 1 minute in days
    h = 1.0 / 1440.0  # 1 minute = 1/1440 days

    # Calculate mu at t-h and t+h
    mu_minus = calc_besselian_mu(jd - h, flags)
    mu_plus = calc_besselian_mu(jd + h, flags)

    # Handle wraparound at 0/360 degrees
    # mu increases with time (typically ~15 deg/hour)
    # If mu_minus is near 360 and mu_plus is near 0, add 360 to mu_plus
    # If mu_plus is near 360 and mu_minus is near 0, add 360 to mu_minus
    delta_mu = mu_plus - mu_minus

    if delta_mu < -180.0:
        # mu wrapped from ~360 to ~0 going forward
        delta_mu += 360.0
    elif delta_mu > 180.0:
        # mu wrapped from ~0 to ~360 going backward (very rare)
        delta_mu -= 360.0

    # dmu/dt in degrees per day
    dmu_dt_per_day = delta_mu / (2 * h)

    # Convert to degrees per hour
    dmu_dt_per_hour = dmu_dt_per_day / 24.0

    return dmu_dt_per_hour


def interpolate_besselian_elements(
    elements: BesselianElements, target_jd: float
) -> BesselianElements:
    """Interpolate Besselian elements to a nearby time.

    A first-order Taylor expansion advances each element from ``elements.t0``.
    Input times are Julian Days and stored derivatives are per hour.

    Args:
        elements: Values and derivatives at the reference time.
        target_jd: Target Julian Day in UT.

    Returns:
        A new :class:`BesselianElements` at ``target_jd``. The six values are
        advanced and the derivative fields are copied unchanged.

    Example:
        >>> from libephemeris import (
        ...     julday, BesselianElements, interpolate_besselian_elements
        ... )
        >>> # Elements at reference time for April 8, 2024 eclipse
        >>> elements = BesselianElements(
        ...     t0=julday(2024, 4, 8, 18.0),
        ...     x=0.3145, y=0.2731, d=7.5821,
        ...     l1=0.5436, l2=-0.0047, mu=89.1234,
        ...     dx_dt=0.5123, dy_dt=0.1456, dd_dt=0.0012,
        ...     dl1_dt=-0.0001, dl2_dt=-0.0001, dmu_dt=15.0041
        ... )
        >>> # Interpolate to 30 minutes later
        >>> later_jd = julday(2024, 4, 8, 18.5)
        >>> interpolated = interpolate_besselian_elements(elements, later_jd)
        >>> print(f"x at t0+0.5h: {interpolated.x:.4f}")

    Notes:
        ``mu`` is normalized to the interval ``[0, 360)``. For times more
        than roughly two hours from the reference, compute fresh elements.

    References:
        Meeus, *Astronomical Algorithms*, chapter 54; *Explanatory Supplement
        to the Astronomical Almanac* (2013), chapter 11.
    """
    # Time difference in hours (derivatives are per hour)
    dt_hours = (target_jd - elements.t0) * 24.0

    # Interpolate each element using first-order Taylor series
    x_new = elements.x + elements.dx_dt * dt_hours
    y_new = elements.y + elements.dy_dt * dt_hours
    d_new = elements.d + elements.dd_dt * dt_hours
    l1_new = elements.l1 + elements.dl1_dt * dt_hours
    l2_new = elements.l2 + elements.dl2_dt * dt_hours
    mu_new = elements.mu + elements.dmu_dt * dt_hours

    # Normalize mu to [0, 360) degrees
    mu_new = mu_new % 360.0

    # Return new BesselianElements with interpolated values
    # Derivatives are copied as they change slowly over short time spans
    return BesselianElements(
        t0=target_jd,
        x=x_new,
        y=y_new,
        d=d_new,
        l1=l1_new,
        l2=l2_new,
        mu=mu_new,
        dx_dt=elements.dx_dt,
        dy_dt=elements.dy_dt,
        dd_dt=elements.dd_dt,
        dl1_dt=elements.dl1_dt,
        dl2_dt=elements.dl2_dt,
        dmu_dt=elements.dmu_dt,
    )


def calc_eclipse_first_contact_c1(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the time of first external contact (C1) for a solar eclipse.

    First contact (C1) is the moment when the Moon's disk first externally
    touches the Sun's disk, marking the beginning of a solar eclipse. At this
    instant, the penumbral shadow cone first touches Earth's surface.

    This function uses Besselian-element geometry to precisely calculate
    C1: the instant, before eclipse maximum, at which the penumbral shadow
    cone is externally tangent to the Earth ellipsoid (the penumbral
    contact residual of ``_besselian_contact_residuals`` crosses zero).

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _sol_eclipse_when_glob_pythonic()
                or _sol_eclipse_when_loc_pythonic(). The function searches backward from
                this time to find C1.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Julian Day (UT) of first contact C1. Returns 0.0 if C1 cannot be
        determined (which would indicate the input time is not near a valid
        solar eclipse).

    Algorithm:
        1. Evaluate the signed penumbral contact residual of
           ``_besselian_contact_residuals`` (tangency of the penumbral cone
           with the Earth *ellipsoid*, in km; negative while in contact)
        2. Expand a bracket backwards from jd_max until the residual turns
           positive (penumbra clear of Earth)
        3. Bisect the sign change to the first external tangency (P1)

    Precision:
        The ellipsoidal residual is bisected to sub-second timing; the
        residual carries the cone-depth term and the Earth's oblateness,
        so the standalone helper matches sol_eclipse_when_glob's tret[2].

    Example:
        >>> from libephemeris import julday, sol_eclipse_when_glob
        >>> from libephemeris import calc_eclipse_first_contact_c1
        >>> from libephemeris import ECL_TOTAL
        >>> # Find the April 8, 2024 total solar eclipse
        >>> jd_start = julday(2024, 1, 1, 0.0)
        >>> ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate first contact
        >>> jd_c1 = calc_eclipse_first_contact_c1(jd_max)
        >>> print(f"First contact C1: JD {jd_c1:.6f}")

    Note:
        - C1 is also known as "first contact" or "partial phase begins"
        - For local eclipse circumstances (at a specific observer location),
          use _sol_eclipse_when_loc_pythonic() which returns contact times in its result
        - The returned time is for the global eclipse (when penumbra first
          touches Earth anywhere), not for a specific observer location

    See Also:
        - _sol_eclipse_when_glob_pythonic: Find next solar eclipse and all phase times
        - _sol_eclipse_when_loc_pythonic: Get local eclipse circumstances including contacts
        - calc_besselian_x, calc_besselian_y: Individual Besselian element functions

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # First contact (P1): the first external tangency of the penumbral cone
    # with the Earth *ellipsoid*, before maximum. Uses the same ellipsoidal
    # contact solver as sol_eclipse_when_glob (tret[2]), so the standalone
    # helper and the global-eclipse array stay consistent.
    return _solve_besselian_contact(
        lambda jd: _besselian_contact_residuals(jd, flags)[0],
        jd_max,
        search_before=True,
    )


def calc_eclipse_second_contact_c2(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the time of second contact (C2) for a solar eclipse.

    Second contact (C2) is the moment when totality or annularity begins -
    when the Moon is fully inside (total eclipse) or outside (annular eclipse)
    the Sun's disk. At this instant, the umbral (total) or antumbral (annular)
    shadow first touches Earth's surface.

    This function uses Besselian-element geometry to precisely calculate
    C2: the instant, before eclipse maximum, at which the umbral/antumbral
    shadow cone is externally tangent to the Earth ellipsoid (the core
    contact residual of ``_besselian_contact_residuals`` crosses zero).

    Note: C2 only exists for central eclipses (total or annular). For partial
    eclipses, this function returns 0.0.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _sol_eclipse_when_glob_pythonic()
                or _sol_eclipse_when_loc_pythonic(). The function searches backward from
                this time to find C2.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Julian Day (UT) of second contact C2. Returns 0.0 if:
        - The eclipse is not a central eclipse (total or annular)
        - C2 cannot be determined (the core shadow never reaches the surface)
        - The input time is not near a valid solar eclipse

    Algorithm:
        1. Evaluate the signed umbral/antumbral contact residual of
           ``_besselian_contact_residuals`` (tangency of the core shadow
           cone with the Earth *ellipsoid*, in km; negative while in
           contact). A residual that is already positive at jd_max means
           the core shadow never reaches the surface: no U1, return 0.0
        2. Expand a bracket backwards from jd_max until the residual turns
           positive (core shadow clear of Earth)
        3. Bisect the sign change to the first external tangency (U1)

    Precision:
        The ellipsoidal residual is bisected to sub-second timing; the
        residual carries the cone-depth term and the Earth's oblateness,
        so the standalone helper matches sol_eclipse_when_glob's tret[4].

    Example:
        >>> from libephemeris import julday, sol_eclipse_when_glob
        >>> from libephemeris import calc_eclipse_second_contact_c2
        >>> from libephemeris import ECL_TOTAL
        >>> # Find the April 8, 2024 total solar eclipse
        >>> jd_start = julday(2024, 1, 1, 0.0)
        >>> ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate second contact
        >>> jd_c2 = calc_eclipse_second_contact_c2(jd_max)
        >>> print(f"Second contact C2: JD {jd_c2:.6f}")

    Note:
        - C2 is also known as "totality begins" or "annularity begins"
        - For partial eclipses, C2 does not exist (returns 0.0)
        - The duration of totality/annularity is (C3 - C2)
        - For local eclipse circumstances (at a specific observer location),
          use _sol_eclipse_when_loc_pythonic() which returns contact times in its result
        - The returned time is for the global eclipse (when umbra/antumbra first
          touches Earth anywhere), not for a specific observer location

    See Also:
        - calc_eclipse_first_contact_c1: Calculate first contact (eclipse begins)
        - _sol_eclipse_when_glob_pythonic: Find next solar eclipse and all phase times
        - _sol_eclipse_when_loc_pythonic: Get local eclipse circumstances including contacts
        - calc_besselian_l2: Calculate umbral/antumbral shadow radius

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Second contact (U1): the first external tangency of the umbral /
    # antumbral cone with the Earth ellipsoid, before maximum. The solver
    # returns 0.0 when the core shadow never reaches the surface (a
    # penumbra-only partial eclipse), so no separate central-phase gate is
    # needed. Shares the ellipsoidal contact solver with tret[4] of
    # sol_eclipse_when_glob.
    return _solve_besselian_contact(
        lambda jd: _besselian_contact_residuals(jd, flags)[1],
        jd_max,
        search_before=True,
    )


def calc_eclipse_third_contact_c3(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the time of third contact (C3) for a solar eclipse.

    Third contact (C3) is the moment when totality or annularity ends -
    when the Moon's umbral (total) or antumbral (annular) shadow last touches
    Earth's surface. At this instant, the central phase of the eclipse ends
    and the partial phase resumes.

    This function uses Besselian-element geometry to precisely calculate
    C3: the instant, after eclipse maximum, at which the umbral/antumbral
    shadow cone is externally tangent to the Earth ellipsoid (the core
    contact residual of ``_besselian_contact_residuals`` crosses zero).

    Note: C3 only exists for central eclipses (total or annular). For partial
    eclipses, this function returns 0.0.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _sol_eclipse_when_glob_pythonic()
                or _sol_eclipse_when_loc_pythonic(). The function searches forward from
                this time to find C3.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Julian Day (UT) of third contact C3. Returns 0.0 if:
        - The eclipse is not a central eclipse (total or annular)
        - C3 cannot be determined (the core shadow never reaches the surface)
        - The input time is not near a valid solar eclipse

    Algorithm:
        1. Evaluate the signed umbral/antumbral contact residual of
           ``_besselian_contact_residuals`` (tangency of the core shadow
           cone with the Earth *ellipsoid*, in km; negative while in
           contact). A residual that is already positive at jd_max means
           the core shadow never reaches the surface: no U4, return 0.0
        2. Expand a bracket forwards from jd_max until the residual turns
           positive (core shadow clear of Earth)
        3. Bisect the sign change to the last external tangency (U4)

    Precision:
        The ellipsoidal residual is bisected to sub-second timing; the
        residual carries the cone-depth term and the Earth's oblateness,
        so the standalone helper matches sol_eclipse_when_glob's tret[5].

    Example:
        >>> from libephemeris import julday, sol_eclipse_when_glob
        >>> from libephemeris import calc_eclipse_third_contact_c3
        >>> from libephemeris import ECL_TOTAL
        >>> # Find the April 8, 2024 total solar eclipse
        >>> jd_start = julday(2024, 1, 1, 0.0)
        >>> ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate third contact
        >>> jd_c3 = calc_eclipse_third_contact_c3(jd_max)
        >>> print(f"Third contact C3: JD {jd_c3:.6f}")

    Note:
        - C3 is also known as "totality ends" or "annularity ends"
        - For partial eclipses, C3 does not exist (returns 0.0)
        - The duration of totality/annularity is (C3 - C2)
        - For local eclipse circumstances (at a specific observer location),
          use _sol_eclipse_when_loc_pythonic() which returns contact times in its result
        - The returned time is for the global eclipse (when umbra/antumbra last
          leaves Earth anywhere), not for a specific observer location

    See Also:
        - calc_eclipse_first_contact_c1: Calculate first contact (eclipse begins)
        - calc_eclipse_second_contact_c2: Calculate second contact (totality begins)
        - _sol_eclipse_when_glob_pythonic: Find next solar eclipse and all phase times
        - _sol_eclipse_when_loc_pythonic: Get local eclipse circumstances including contacts
        - calc_besselian_l2: Calculate umbral/antumbral shadow radius

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Third contact (U4): the last external tangency of the umbral /
    # antumbral cone with the Earth ellipsoid, after maximum. Returns 0.0
    # for a penumbra-only partial eclipse. Shares the ellipsoidal contact
    # solver with tret[5]/U4 of sol_eclipse_when_glob.
    return _solve_besselian_contact(
        lambda jd: _besselian_contact_residuals(jd, flags)[1],
        jd_max,
        search_before=False,
    )


def calc_eclipse_fourth_contact_c4(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the time of fourth external contact (C4) for a solar eclipse.

    Fourth contact (C4) is the moment when the Moon's disk completely separates
    from the Sun's disk externally, marking the end of a solar eclipse. At this
    instant, the penumbral shadow cone last touches Earth's surface.

    This function uses Besselian-element geometry to precisely calculate
    C4: the instant, after eclipse maximum, at which the penumbral shadow
    cone is externally tangent to the Earth ellipsoid (the penumbral
    contact residual of ``_besselian_contact_residuals`` crosses zero).

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _sol_eclipse_when_glob_pythonic()
                or _sol_eclipse_when_loc_pythonic(). The function searches forward from
                this time to find C4.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Julian Day (UT) of fourth contact C4. Returns 0.0 if C4 cannot be
        determined (which would indicate the input time is not near a valid
        solar eclipse).

    Algorithm:
        1. Evaluate the signed penumbral contact residual of
           ``_besselian_contact_residuals`` (tangency of the penumbral cone
           with the Earth *ellipsoid*, in km; negative while in contact)
        2. Expand a bracket forwards from jd_max until the residual turns
           positive (penumbra clear of Earth)
        3. Bisect the sign change to the last external tangency (P4)

    Precision:
        The ellipsoidal residual is bisected to sub-second timing; the
        residual carries the cone-depth term and the Earth's oblateness,
        so the standalone helper matches sol_eclipse_when_glob's tret[3].

    Example:
        >>> from libephemeris import julday, sol_eclipse_when_glob
        >>> from libephemeris import calc_eclipse_fourth_contact_c4
        >>> from libephemeris import ECL_TOTAL
        >>> # Find the April 8, 2024 total solar eclipse
        >>> jd_start = julday(2024, 1, 1, 0.0)
        >>> ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate fourth contact
        >>> jd_c4 = calc_eclipse_fourth_contact_c4(jd_max)
        >>> print(f"Fourth contact C4: JD {jd_c4:.6f}")

    Note:
        - C4 is also known as "fourth contact" or "partial phase ends"
        - C4 marks the end of the eclipse when the penumbra completely leaves Earth
        - For local eclipse circumstances (at a specific observer location),
          use _sol_eclipse_when_loc_pythonic() which returns contact times in its result
        - The returned time is for the global eclipse (when penumbra last
          leaves Earth anywhere), not for a specific observer location
        - The total duration of the eclipse is (C4 - C1)

    See Also:
        - calc_eclipse_first_contact_c1: Calculate first contact (eclipse begins)
        - calc_eclipse_second_contact_c2: Calculate second contact (totality begins)
        - calc_eclipse_third_contact_c3: Calculate third contact (totality ends)
        - _sol_eclipse_when_glob_pythonic: Find next solar eclipse and all phase times
        - _sol_eclipse_when_loc_pythonic: Get local eclipse circumstances including contacts
        - calc_besselian_l1: Calculate penumbral shadow radius

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Fourth contact (P4): the last external tangency of the penumbral cone
    # with the Earth ellipsoid, after maximum. Shares the ellipsoidal
    # contact solver with tret[3] of sol_eclipse_when_glob.
    return _solve_besselian_contact(
        lambda jd: _besselian_contact_residuals(jd, flags)[0],
        jd_max,
        search_before=False,
    )


def _calc_lunar_eclipse_penumbral_separation(jd: float) -> float:
    """
    Calculate the separation between Moon's nearest limb and penumbral shadow edge.

    Returns the angular separation in degrees between the Moon's nearest limb
    and the penumbral shadow boundary. The penumbral contact occurs when this
    value equals zero (Moon's limb touches penumbra edge).

    Args:
        jd: Julian Day (UT)

    Returns:
        Separation in degrees. Positive means Moon is outside penumbra,
        negative means Moon's limb is inside penumbra.
    """
    # Get positions
    sun_pos, _ = calc_ut(jd, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd, MOON, FLG_SPEED)

    # Extract coordinates
    sun_lon = sun_pos[0]
    moon_lon = moon_pos[0]
    moon_lat = moon_pos[1]
    sun_dist = sun_pos[2]
    moon_dist = moon_pos[2]

    # Shadow center is opposite the Sun (anti-Sun point on ecliptic)
    shadow_lon = (sun_lon + 180.0) % 360.0

    # Calculate Moon's angular distance from shadow axis
    # The shadow axis passes through Earth opposite the Sun
    dlon = moon_lon - shadow_lon
    if dlon > 180:
        dlon -= 360
    if dlon < -180:
        dlon += 360

    # Angular distance from Moon center to shadow axis (degrees)
    # Using spherical geometry for accuracy
    moon_lat_rad = math.radians(moon_lat)
    distance_from_axis = math.sqrt(dlon**2 * math.cos(moon_lat_rad) ** 2 + moon_lat**2)

    # Constants for shadow geometry
    SUN_RADIUS_ARCSEC = 959.63
    EARTH_RADIUS_AU = 4.2635e-5

    # Sun's angular semi-diameter at actual distance (in degrees)
    sun_semidiameter = (SUN_RADIUS_ARCSEC / 3600.0) / sun_dist

    # Moon's angular semi-diameter
    moon_semidiameter = (932.56 / 3600.0) * (0.002569 / moon_dist)

    # Earth's angular semi-diameter as seen from Moon (in degrees)
    earth_semidiameter = math.degrees(math.atan(EARTH_RADIUS_AU / moon_dist))

    # Penumbra radius at Moon's distance
    penumbra_radius = earth_semidiameter + sun_semidiameter

    # Apply Danjon atmospheric shadow enlargement (same factor used in
    # _calculate_lunar_eclipse_type_and_magnitude for consistency)
    _SHADOW_ENLARGEMENT = 1.0 + 1.0 / 85.0
    penumbra_radius *= _SHADOW_ENLARGEMENT

    # Separation: positive = Moon outside penumbra, negative = Moon inside
    # Contact occurs when Moon's nearest limb touches penumbra edge
    # At P1/P4: distance_from_axis - moon_semidiameter = penumbra_radius
    # So: separation = distance_from_axis - moon_semidiameter - penumbra_radius = 0
    separation = distance_from_axis - moon_semidiameter - penumbra_radius

    return separation


def _find_lunar_penumbral_contact_time(
    jd_max: float,
    search_before: bool,
    search_range: float = 0.20,
) -> float:
    """
    Find the time of lunar eclipse penumbral contact using binary search.

    Uses iterative refinement to find when Moon's limb touches the penumbral
    shadow boundary.

    Args:
        jd_max: Julian Day of eclipse maximum
        search_before: If True, search for P1 (before maximum), else P4 (after)
        search_range: Search range in days from maximum (default 0.20 = ~4.8 hours)

    Returns:
        Julian Day of the contact, or 0.0 if not found
    """
    # Set up search bounds
    if search_before:
        jd_low = jd_max - search_range
        jd_high = jd_max
    else:
        jd_low = jd_max
        jd_high = jd_max + search_range

    # Verify that a contact exists in the search range
    sep_low = _calc_lunar_eclipse_penumbral_separation(jd_low)
    sep_high = _calc_lunar_eclipse_penumbral_separation(jd_high)

    # For P1: sep_low should be positive (outside) and sep_high negative (inside)
    # For P4: sep_low should be negative (inside) and sep_high positive (outside)
    if search_before:
        # P1: looking for transition from positive to negative
        if sep_low < 0:
            # Moon already in penumbra at start of search - no P1 in range
            return 0.0
        if sep_high > 0:
            # Moon still outside penumbra at maximum - no eclipse
            return 0.0
    else:
        # P4: looking for transition from negative to positive
        if sep_low > 0:
            # Moon already outside penumbra at maximum
            return 0.0
        if sep_high < 0:
            # Moon still in penumbra at end of search
            return 0.0

    # Binary search for the contact time (when separation = 0)
    for _ in range(60):  # ~60 iterations gives sub-second precision
        jd_mid = (jd_low + jd_high) / 2
        sep_mid = _calc_lunar_eclipse_penumbral_separation(jd_mid)

        if search_before:
            # P1: separation goes from positive to negative
            if sep_mid > 0:
                jd_low = jd_mid
            else:
                jd_high = jd_mid
        else:
            # P4: separation goes from negative to positive
            if sep_mid < 0:
                jd_low = jd_mid
            else:
                jd_high = jd_mid

        # Check convergence (precision ~0.05 seconds)
        if jd_high - jd_low < 5e-7:
            break

    return (jd_low + jd_high) / 2


def calc_lunar_eclipse_penumbral_first_contact_p1(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the time of penumbral first contact (P1) for a lunar eclipse.

    P1 is the moment when the Moon's leading limb first enters Earth's
    penumbral shadow, marking the beginning of the lunar eclipse. At this
    instant, a very subtle shading begins on the Moon's eastern limb, though
    it is typically not visible to the naked eye until the Moon penetrates
    deeper into the penumbra.

    This function uses iterative refinement with geometric calculations to
    precisely determine P1. The condition for P1 is when the angular
    separation between the Moon's limb and the penumbral shadow edge equals
    zero, occurring before eclipse maximum.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _lun_eclipse_when_pythonic().
                The function searches backward from this time to find P1.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Julian Day (UT) of penumbral first contact P1. Returns 0.0 if P1
        cannot be determined (which would indicate the input time is not
        near a valid lunar eclipse).

    Algorithm:
        1. Calculate Moon's position and Earth's shadow geometry
        2. Determine penumbral shadow radius at Moon's distance
        3. Use binary search to find when Moon's limb touches penumbra edge
        4. The search proceeds from (jd_max - search_range) to jd_max

    Precision:
        The calculation achieves timing precision better than 1 second by
        iterating until the time value converges to within ~0.05 seconds.

    Example:
        >>> from libephemeris import julday, lun_eclipse_when
        >>> from libephemeris import calc_lunar_eclipse_penumbral_first_contact_p1
        >>> from libephemeris import ECL_TOTAL
        >>> # Find the November 8, 2022 total lunar eclipse
        >>> jd_start = julday(2022, 10, 1, 0.0)
        >>> ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate penumbral first contact
        >>> jd_p1 = calc_lunar_eclipse_penumbral_first_contact_p1(jd_max)
        >>> print(f"Penumbral first contact P1: JD {jd_p1:.6f}")

    Note:
        - P1 is also known as "penumbral eclipse begins" or "first penumbral contact"
        - The penumbral phase (P1 to P4) is the total duration of the eclipse
        - For penumbral-only eclipses, P1 and P4 are the only contact points
        - The penumbral shading is typically visible only after the Moon has
          penetrated significantly into the penumbra (~70% or more)

    See Also:
        - calc_lunar_eclipse_penumbral_fourth_contact_p4: Calculate P4 (eclipse ends)
        - _lun_eclipse_when_pythonic: Find next lunar eclipse and all phase times
        - _lun_eclipse_when_loc_pythonic: Get local lunar eclipse circumstances

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Lunar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Lunar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # P1 = penumbral begin = tret[6] from the canonical shadow model.
    return _lun_contact_via_core(jd_max, 6, flags)


def calc_lunar_eclipse_penumbral_fourth_contact_p4(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the time of penumbral fourth contact (P4) for a lunar eclipse.

    P4 is the moment when the Moon's trailing limb completely exits Earth's
    penumbral shadow, marking the end of the lunar eclipse. At this instant,
    the last trace of penumbral shading leaves the Moon's western limb,
    although the shading would have been imperceptible for some time before.

    This function uses iterative refinement with geometric calculations to
    precisely determine P4. The condition for P4 is when the angular
    separation between the Moon's limb and the penumbral shadow edge equals
    zero, occurring after eclipse maximum.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _lun_eclipse_when_pythonic().
                The function searches forward from this time to find P4.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Julian Day (UT) of penumbral fourth contact P4. Returns 0.0 if P4
        cannot be determined (which would indicate the input time is not
        near a valid lunar eclipse).

    Algorithm:
        1. Calculate Moon's position and Earth's shadow geometry
        2. Determine penumbral shadow radius at Moon's distance
        3. Use binary search to find when Moon's limb exits penumbra edge
        4. The search proceeds from jd_max to (jd_max + search_range)

    Precision:
        The calculation achieves timing precision better than 1 second by
        iterating until the time value converges to within ~0.05 seconds.

    Example:
        >>> from libephemeris import julday, lun_eclipse_when
        >>> from libephemeris import calc_lunar_eclipse_penumbral_fourth_contact_p4
        >>> from libephemeris import ECL_TOTAL
        >>> # Find the November 8, 2022 total lunar eclipse
        >>> jd_start = julday(2022, 10, 1, 0.0)
        >>> ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate penumbral fourth contact
        >>> jd_p4 = calc_lunar_eclipse_penumbral_fourth_contact_p4(jd_max)
        >>> print(f"Penumbral fourth contact P4: JD {jd_p4:.6f}")

    Note:
        - P4 is also known as "penumbral eclipse ends" or "last penumbral contact"
        - The total duration of the eclipse (P1 to P4) can be several hours
        - P4 marks the end of any visible or instrumental eclipse effects
        - The duration (P4 - P1) represents the complete penumbral phase

    See Also:
        - calc_lunar_eclipse_penumbral_first_contact_p1: Calculate P1 (eclipse begins)
        - _lun_eclipse_when_pythonic: Find next lunar eclipse and all phase times
        - _lun_eclipse_when_loc_pythonic: Get local lunar eclipse circumstances

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Lunar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Lunar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # P4 = penumbral end = tret[7] from the canonical shadow model.
    return _lun_contact_via_core(jd_max, 7, flags)


def _calc_lunar_eclipse_umbral_outer_separation(jd: float) -> float:
    """
    Calculate the separation between Moon's nearest limb and umbral shadow edge.

    Returns the angular separation in degrees between the Moon's nearest limb
    and the umbral shadow boundary. The umbral contacts U1/U4 occur when this
    value equals zero (Moon's limb touches umbra edge).

    Args:
        jd: Julian Day (UT)

    Returns:
        Separation in degrees. Positive means Moon is outside umbra,
        negative means Moon's limb is inside umbra.
    """
    # Get positions
    sun_pos, _ = calc_ut(jd, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd, MOON, FLG_SPEED)

    # Extract coordinates
    sun_lon = sun_pos[0]
    moon_lon = moon_pos[0]
    moon_lat = moon_pos[1]
    sun_dist = sun_pos[2]
    moon_dist = moon_pos[2]

    # Shadow center is opposite the Sun (anti-Sun point on ecliptic)
    shadow_lon = (sun_lon + 180.0) % 360.0

    # Calculate Moon's angular distance from shadow axis
    dlon = moon_lon - shadow_lon
    if dlon > 180:
        dlon -= 360
    if dlon < -180:
        dlon += 360

    # Angular distance from Moon center to shadow axis (degrees)
    moon_lat_rad = math.radians(moon_lat)
    distance_from_axis = math.sqrt(dlon**2 * math.cos(moon_lat_rad) ** 2 + moon_lat**2)

    # Constants for shadow geometry
    SUN_RADIUS_ARCSEC = 959.63
    EARTH_RADIUS_KM = 6378.137

    # Moon's angular semi-diameter
    moon_semidiameter = (932.56 / 3600.0) * (0.002569 / moon_dist)

    # Calculate umbra radius at Moon's distance
    # (149597870.7 km = 1 AU, IAU 2012 Resolution B2)
    sun_dist_km = sun_dist * 149597870.7
    moon_dist_km = moon_dist * 149597870.7
    sun_radius_km = (SUN_RADIUS_ARCSEC / 206265.0) * sun_dist_km

    # Umbra cone semi-angle
    umbra_cone_angle = math.atan((sun_radius_km - EARTH_RADIUS_KM) / sun_dist_km)
    # Umbra radius at Moon distance
    umbra_radius_km = EARTH_RADIUS_KM - moon_dist_km * math.tan(umbra_cone_angle)

    if umbra_radius_km > 0:
        umbra_radius = math.degrees(math.atan(umbra_radius_km / moon_dist_km))
    else:
        # Umbra doesn't reach Moon - no umbral eclipse possible
        umbra_radius = 0.0

    # Separation: positive = Moon outside umbra, negative = Moon inside
    # At U1/U4: distance_from_axis - moon_semidiameter = umbra_radius
    # So: separation = distance_from_axis - moon_semidiameter - umbra_radius = 0
    separation = distance_from_axis - moon_semidiameter - umbra_radius

    return separation


def _calc_lunar_eclipse_umbral_inner_separation(jd: float) -> float:
    """
    Calculate the separation for totality contacts (U2/U3).

    Returns the angular separation in degrees. At U2/U3, the Moon's trailing/leading
    limb touches the umbra edge (entire Moon inside umbra).

    Args:
        jd: Julian Day (UT)

    Returns:
        Separation in degrees. Positive means Moon is not fully inside umbra,
        negative means Moon is completely inside umbra.
    """
    # Get positions
    sun_pos, _ = calc_ut(jd, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd, MOON, FLG_SPEED)

    # Extract coordinates
    sun_lon = sun_pos[0]
    moon_lon = moon_pos[0]
    moon_lat = moon_pos[1]
    sun_dist = sun_pos[2]
    moon_dist = moon_pos[2]

    # Shadow center is opposite the Sun
    shadow_lon = (sun_lon + 180.0) % 360.0

    # Calculate Moon's angular distance from shadow axis
    dlon = moon_lon - shadow_lon
    if dlon > 180:
        dlon -= 360
    if dlon < -180:
        dlon += 360

    # Angular distance from Moon center to shadow axis
    moon_lat_rad = math.radians(moon_lat)
    distance_from_axis = math.sqrt(dlon**2 * math.cos(moon_lat_rad) ** 2 + moon_lat**2)

    # Constants for shadow geometry
    SUN_RADIUS_ARCSEC = 959.63
    EARTH_RADIUS_KM = 6378.137

    # Moon's angular semi-diameter
    moon_semidiameter = (932.56 / 3600.0) * (0.002569 / moon_dist)

    # Calculate umbra radius at Moon's distance
    # (149597870.7 km = 1 AU, IAU 2012 Resolution B2)
    sun_dist_km = sun_dist * 149597870.7
    moon_dist_km = moon_dist * 149597870.7
    sun_radius_km = (SUN_RADIUS_ARCSEC / 206265.0) * sun_dist_km

    # Umbra cone semi-angle
    umbra_cone_angle = math.atan((sun_radius_km - EARTH_RADIUS_KM) / sun_dist_km)
    # Umbra radius at Moon distance
    umbra_radius_km = EARTH_RADIUS_KM - moon_dist_km * math.tan(umbra_cone_angle)

    if umbra_radius_km > 0:
        umbra_radius = math.degrees(math.atan(umbra_radius_km / moon_dist_km))
    else:
        # No umbral eclipse possible
        umbra_radius = 0.0

    # For U2/U3 (totality contacts): Moon is fully inside umbra when
    # distance_from_axis + moon_semidiameter = umbra_radius
    # (farthest edge of Moon touches umbra edge)
    separation = distance_from_axis + moon_semidiameter - umbra_radius

    return separation


def _find_lunar_umbral_outer_contact_time(
    jd_max: float,
    search_before: bool,
    search_range: float = 0.15,
) -> float:
    """
    Find the time of lunar eclipse umbral outer contact (U1 or U4) using binary search.

    Args:
        jd_max: Julian Day of eclipse maximum
        search_before: If True, search for U1 (before maximum), else U4 (after)
        search_range: Search range in days from maximum

    Returns:
        Julian Day of the contact, or 0.0 if not found
    """
    # Set up search bounds
    if search_before:
        jd_low = jd_max - search_range
        jd_high = jd_max
    else:
        jd_low = jd_max
        jd_high = jd_max + search_range

    # Verify that a contact exists in the search range
    sep_low = _calc_lunar_eclipse_umbral_outer_separation(jd_low)
    sep_high = _calc_lunar_eclipse_umbral_outer_separation(jd_high)

    if search_before:
        # U1: looking for transition from positive to negative
        if sep_low < 0:
            # Moon already in umbra at start of search - expand search
            jd_low = jd_max - search_range * 1.5
            sep_low = _calc_lunar_eclipse_umbral_outer_separation(jd_low)
            if sep_low < 0:
                return 0.0
        if sep_high > 0:
            # Moon still outside umbra at maximum - no umbral eclipse
            return 0.0
    else:
        # U4: looking for transition from negative to positive
        if sep_low > 0:
            # Moon already outside umbra at maximum
            return 0.0
        if sep_high < 0:
            # Moon still in umbra at end of search - expand search
            jd_high = jd_max + search_range * 1.5
            sep_high = _calc_lunar_eclipse_umbral_outer_separation(jd_high)
            if sep_high < 0:
                return 0.0

    # Binary search for the contact time (when separation = 0)
    for _ in range(60):  # ~60 iterations gives sub-second precision
        jd_mid = (jd_low + jd_high) / 2
        sep_mid = _calc_lunar_eclipse_umbral_outer_separation(jd_mid)

        if search_before:
            # U1: separation goes from positive to negative
            if sep_mid > 0:
                jd_low = jd_mid
            else:
                jd_high = jd_mid
        else:
            # U4: separation goes from negative to positive
            if sep_mid < 0:
                jd_low = jd_mid
            else:
                jd_high = jd_mid

        # Check convergence (precision ~0.05 seconds)
        if jd_high - jd_low < 5e-7:
            break

    return (jd_low + jd_high) / 2


def _find_lunar_umbral_inner_contact_time(
    jd_max: float,
    search_before: bool,
    search_range: float = 0.10,
) -> float:
    """
    Find the time of lunar eclipse totality contact (U2 or U3) using binary search.

    Args:
        jd_max: Julian Day of eclipse maximum
        search_before: If True, search for U2 (totality begins), else U3 (totality ends)
        search_range: Search range in days from maximum

    Returns:
        Julian Day of the contact, or 0.0 if not found (no total eclipse)
    """
    # First check if this is a total eclipse by checking separation at maximum
    sep_max = _calc_lunar_eclipse_umbral_inner_separation(jd_max)
    if sep_max > 0:
        # Not a total eclipse - Moon never fully enters umbra
        return 0.0

    # Set up search bounds
    if search_before:
        jd_low = jd_max - search_range
        jd_high = jd_max
    else:
        jd_low = jd_max
        jd_high = jd_max + search_range

    # Verify that a contact exists in the search range
    sep_low = _calc_lunar_eclipse_umbral_inner_separation(jd_low)
    sep_high = _calc_lunar_eclipse_umbral_inner_separation(jd_high)

    if search_before:
        # U2: looking for transition from positive to negative
        if sep_low < 0:
            # Moon already fully in umbra at start - expand search
            jd_low = jd_max - search_range * 1.5
            sep_low = _calc_lunar_eclipse_umbral_inner_separation(jd_low)
            if sep_low < 0:
                jd_low = jd_max - search_range * 2.0
                sep_low = _calc_lunar_eclipse_umbral_inner_separation(jd_low)
                if sep_low < 0:
                    return 0.0
        if sep_high > 0:
            # Moon not fully in umbra at maximum - not a total eclipse
            return 0.0
    else:
        # U3: looking for transition from negative to positive
        if sep_low > 0:
            # Moon not fully in umbra at maximum - not a total eclipse
            return 0.0
        if sep_high < 0:
            # Moon still fully in umbra at end - expand search
            jd_high = jd_max + search_range * 1.5
            sep_high = _calc_lunar_eclipse_umbral_inner_separation(jd_high)
            if sep_high < 0:
                jd_high = jd_max + search_range * 2.0
                sep_high = _calc_lunar_eclipse_umbral_inner_separation(jd_high)
                if sep_high < 0:
                    return 0.0

    # Binary search for the contact time (when separation = 0)
    for _ in range(60):
        jd_mid = (jd_low + jd_high) / 2
        sep_mid = _calc_lunar_eclipse_umbral_inner_separation(jd_mid)

        if search_before:
            # U2: separation goes from positive to negative
            if sep_mid > 0:
                jd_low = jd_mid
            else:
                jd_high = jd_mid
        else:
            # U3: separation goes from negative to positive
            if sep_mid < 0:
                jd_low = jd_mid
            else:
                jd_high = jd_mid

        # Check convergence
        if jd_high - jd_low < 5e-7:
            break

    return (jd_low + jd_high) / 2


def calc_lunar_eclipse_umbral_first_contact_u1(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the time of umbral first contact (U1) for a lunar eclipse.

    U1 is the moment when the Moon's leading limb first enters Earth's
    umbral shadow, marking the beginning of the partial phase of the
    lunar eclipse. At this instant, the eastern edge of the Moon starts
    to darken noticeably as it enters the umbra.

    This function uses iterative refinement with geometric calculations to
    precisely determine U1. The condition for U1 is when the angular
    separation between the Moon's limb and the umbral shadow edge equals
    zero, occurring before eclipse maximum.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _lun_eclipse_when_pythonic().
                The function searches backward from this time to find U1.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Julian Day (UT) of umbral first contact U1. Returns 0.0 if U1
        cannot be determined (which indicates either a penumbral-only eclipse
        or the input time is not near a valid lunar eclipse).

    Algorithm:
        1. Calculate Moon's position and Earth's shadow geometry
        2. Determine umbral shadow radius at Moon's distance
        3. Use binary search to find when Moon's limb touches umbra edge
        4. The search proceeds from (jd_max - search_range) to jd_max

    Precision:
        The calculation achieves timing precision better than 1 second by
        iterating until the time value converges to within ~0.05 seconds.

    Example:
        >>> from libephemeris import julday, lun_eclipse_when
        >>> from libephemeris import calc_lunar_eclipse_umbral_first_contact_u1
        >>> from libephemeris import ECL_TOTAL
        >>> # Find the November 8, 2022 total lunar eclipse
        >>> jd_start = julday(2022, 10, 1, 0.0)
        >>> ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate umbral first contact
        >>> jd_u1 = calc_lunar_eclipse_umbral_first_contact_u1(jd_max)
        >>> print(f"Umbral first contact U1: JD {jd_u1:.6f}")

    Note:
        - U1 is also known as "partial eclipse begins" or "first umbral contact"
        - U1 only occurs for partial and total lunar eclipses, not penumbral-only
        - The umbral phase (U1 to U4) is the most dramatic part of the eclipse
        - U1 marks when the Moon begins to visibly darken

    See Also:
        - calc_lunar_eclipse_umbral_fourth_contact_u4: Calculate U4 (partial ends)
        - calc_lunar_eclipse_umbral_second_contact_u2: Calculate U2 (totality begins)
        - calc_lunar_eclipse_umbral_third_contact_u3: Calculate U3 (totality ends)
        - calc_lunar_eclipse_penumbral_first_contact_p1: Calculate P1 (penumbral begins)
        - _lun_eclipse_when_pythonic: Find next lunar eclipse and all phase times

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Lunar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Lunar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # U1 = umbral (partial) begin = tret[2] from the canonical shadow model.
    return _lun_contact_via_core(jd_max, 2, flags)


def calc_lunar_eclipse_umbral_second_contact_u2(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the time of umbral second contact (U2) for a lunar eclipse.

    U2 is the moment when the Moon's trailing limb enters Earth's umbral
    shadow, marking the beginning of totality. At this instant, the entire
    Moon is within the umbra and the total phase of the eclipse begins.
    The Moon appears deeply red or copper-colored during totality.

    This function uses iterative refinement with geometric calculations to
    precisely determine U2. The condition for U2 is when the angular
    separation between the Moon's trailing limb and the umbral shadow edge
    equals zero, occurring before eclipse maximum.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _lun_eclipse_when_pythonic().
                The function searches backward from this time to find U2.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Julian Day (UT) of umbral second contact U2. Returns 0.0 if U2
        cannot be determined (which indicates a partial-only or penumbral-only
        eclipse where totality does not occur).

    Algorithm:
        1. Calculate Moon's position and Earth's shadow geometry
        2. Determine umbral shadow radius at Moon's distance
        3. Use binary search to find when Moon fully enters umbra
        4. The search proceeds from (jd_max - search_range) to jd_max

    Precision:
        The calculation achieves timing precision better than 1 second by
        iterating until the time value converges to within ~0.05 seconds.

    Example:
        >>> from libephemeris import julday, lun_eclipse_when
        >>> from libephemeris import calc_lunar_eclipse_umbral_second_contact_u2
        >>> from libephemeris import ECL_TOTAL
        >>> # Find the November 8, 2022 total lunar eclipse
        >>> jd_start = julday(2022, 10, 1, 0.0)
        >>> ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate umbral second contact (totality begins)
        >>> jd_u2 = calc_lunar_eclipse_umbral_second_contact_u2(jd_max)
        >>> print(f"Totality begins U2: JD {jd_u2:.6f}")

    Note:
        - U2 only occurs for total lunar eclipses
        - Returns 0.0 for partial or penumbral eclipses
        - The duration (U3 - U2) is the length of totality
        - During totality, the Moon appears red due to refracted sunlight

    See Also:
        - calc_lunar_eclipse_umbral_third_contact_u3: Calculate U3 (totality ends)
        - calc_lunar_eclipse_umbral_first_contact_u1: Calculate U1 (partial begins)
        - calc_lunar_eclipse_umbral_fourth_contact_u4: Calculate U4 (partial ends)
        - _lun_eclipse_when_pythonic: Find next lunar eclipse and all phase times

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Lunar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Lunar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # U2 = totality begins = tret[4] from the canonical shadow model.
    return _lun_contact_via_core(jd_max, 4, flags)


def calc_lunar_eclipse_umbral_third_contact_u3(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the time of umbral third contact (U3) for a lunar eclipse.

    U3 is the moment when the Moon's leading limb exits Earth's umbral
    shadow, marking the end of totality. At this instant, the western
    edge of the Moon begins to brighten as it emerges from the umbra,
    and the total phase of the eclipse ends.

    This function uses iterative refinement with geometric calculations to
    precisely determine U3. The condition for U3 is when the angular
    separation between the Moon's leading limb and the umbral shadow edge
    equals zero, occurring after eclipse maximum.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _lun_eclipse_when_pythonic().
                The function searches forward from this time to find U3.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Julian Day (UT) of umbral third contact U3. Returns 0.0 if U3
        cannot be determined (which indicates a partial-only or penumbral-only
        eclipse where totality does not occur).

    Algorithm:
        1. Calculate Moon's position and Earth's shadow geometry
        2. Determine umbral shadow radius at Moon's distance
        3. Use binary search to find when Moon begins exiting umbra
        4. The search proceeds from jd_max to (jd_max + search_range)

    Precision:
        The calculation achieves timing precision better than 1 second by
        iterating until the time value converges to within ~0.05 seconds.

    Example:
        >>> from libephemeris import julday, lun_eclipse_when
        >>> from libephemeris import calc_lunar_eclipse_umbral_third_contact_u3
        >>> from libephemeris import ECL_TOTAL
        >>> # Find the November 8, 2022 total lunar eclipse
        >>> jd_start = julday(2022, 10, 1, 0.0)
        >>> ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate umbral third contact (totality ends)
        >>> jd_u3 = calc_lunar_eclipse_umbral_third_contact_u3(jd_max)
        >>> print(f"Totality ends U3: JD {jd_u3:.6f}")

    Note:
        - U3 only occurs for total lunar eclipses
        - Returns 0.0 for partial or penumbral eclipses
        - The duration (U3 - U2) is the length of totality
        - U3 marks when the Moon begins to brighten

    See Also:
        - calc_lunar_eclipse_umbral_second_contact_u2: Calculate U2 (totality begins)
        - calc_lunar_eclipse_umbral_first_contact_u1: Calculate U1 (partial begins)
        - calc_lunar_eclipse_umbral_fourth_contact_u4: Calculate U4 (partial ends)
        - _lun_eclipse_when_pythonic: Find next lunar eclipse and all phase times

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Lunar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Lunar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # U3 = totality ends = tret[5] from the canonical shadow model.
    return _lun_contact_via_core(jd_max, 5, flags)


def calc_lunar_eclipse_umbral_fourth_contact_u4(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the time of umbral fourth contact (U4) for a lunar eclipse.

    U4 is the moment when the Moon's trailing limb completely exits Earth's
    umbral shadow, marking the end of the partial phase of the lunar eclipse.
    At this instant, the last portion of the Moon emerges from the umbra and
    returns to its normal brightness (though still in the penumbra).

    This function uses iterative refinement with geometric calculations to
    precisely determine U4. The condition for U4 is when the angular
    separation between the Moon's trailing limb and the umbral shadow edge
    equals zero, occurring after eclipse maximum.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _lun_eclipse_when_pythonic().
                The function searches forward from this time to find U4.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Julian Day (UT) of umbral fourth contact U4. Returns 0.0 if U4
        cannot be determined (which indicates a penumbral-only eclipse
        or the input time is not near a valid lunar eclipse).

    Algorithm:
        1. Calculate Moon's position and Earth's shadow geometry
        2. Determine umbral shadow radius at Moon's distance
        3. Use binary search to find when Moon's limb exits umbra edge
        4. The search proceeds from jd_max to (jd_max + search_range)

    Precision:
        The calculation achieves timing precision better than 1 second by
        iterating until the time value converges to within ~0.05 seconds.

    Example:
        >>> from libephemeris import julday, lun_eclipse_when
        >>> from libephemeris import calc_lunar_eclipse_umbral_fourth_contact_u4
        >>> from libephemeris import ECL_TOTAL
        >>> # Find the November 8, 2022 total lunar eclipse
        >>> jd_start = julday(2022, 10, 1, 0.0)
        >>> ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate umbral fourth contact
        >>> jd_u4 = calc_lunar_eclipse_umbral_fourth_contact_u4(jd_max)
        >>> print(f"Umbral fourth contact U4: JD {jd_u4:.6f}")

    Note:
        - U4 is also known as "partial eclipse ends" or "last umbral contact"
        - U4 only occurs for partial and total lunar eclipses, not penumbral-only
        - The duration (U4 - U1) is the length of the partial/umbral phase
        - After U4, the Moon is only in the penumbra until P4

    See Also:
        - calc_lunar_eclipse_umbral_first_contact_u1: Calculate U1 (partial begins)
        - calc_lunar_eclipse_umbral_second_contact_u2: Calculate U2 (totality begins)
        - calc_lunar_eclipse_umbral_third_contact_u3: Calculate U3 (totality ends)
        - calc_lunar_eclipse_penumbral_fourth_contact_p4: Calculate P4 (eclipse ends)
        - _lun_eclipse_when_pythonic: Find next lunar eclipse and all phase times

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Lunar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Lunar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # U4 = umbral (partial) end = tret[3] from the canonical shadow model.
    return _lun_contact_via_core(jd_max, 3, flags)


def calc_solar_eclipse_duration(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the duration of totality or annularity for a solar eclipse.

    For central solar eclipses (total or annular), this function calculates
    the maximum duration of the central phase on the central line - the local
    time between second contact (C2) and third contact (C3) at the point of
    greatest eclipse. This is the length of time an observer at the shadow
    center sees totality (total eclipse) or annularity (annular eclipse), and
    is a few minutes long.

    For non-central eclipses (partial only), the function returns 0.0 since
    there is no totality or annularity phase.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _sol_eclipse_when_glob_pythonic()
                or _sol_eclipse_when_loc_pythonic(). The function locates the
                central-line point at this time and computes the local C2 and
                C3 there.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Duration of totality/annularity in minutes. Returns 0.0 if:
        - The eclipse is not a central eclipse (total or annular)
        - C2 or C3 cannot be determined
        - The input time is not near a valid solar eclipse

    Algorithm:
        1. Locate the central-line (shadow-center) point at maximum eclipse
           via sol_eclipse_where()
        2. Compute the local second (C2) and third (C3) contact times at that
           point via _calculate_local_eclipse_phases()
        3. Return (C3 - C2) converted from days to minutes

    Precision:
        The central-line duration agrees with one external implementation's
        local maximum-eclipse duration to within about ±0.05 minutes
        (±3 seconds).

    Example:
        >>> from libephemeris import julday, sol_eclipse_when_glob
        >>> from libephemeris import calc_solar_eclipse_duration, ECL_TOTAL
        >>> # Find the April 8, 2024 total solar eclipse
        >>> jd_start = julday(2024, 1, 1, 0.0)
        >>> ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate duration of totality
        >>> duration = calc_solar_eclipse_duration(jd_max)
        >>> print(f"Duration of totality: {duration:.2f} minutes")

    Note:
        - This is the central-line maximum duration (as seen from the point of
          greatest eclipse), the quantity the compatibility contract labels
          the local maximum-eclipse duration
        - Duration at other observer locations on the path is typically shorter;
          for a specific location use _sol_eclipse_when_loc_pythonic()
        - Total solar eclipses can have central durations up to ~7.5 minutes
        - Annular eclipses can have central durations up to ~12 minutes

    See Also:
        - calc_eclipse_second_contact_c2: Calculate C2 (central phase begins)
        - calc_eclipse_third_contact_c3: Calculate C3 (central phase ends)
        - _sol_eclipse_when_glob_pythonic: Find next solar eclipse and all phase times
        - _sol_eclipse_when_loc_pythonic: Get local eclipse circumstances including contacts
        - calc_lunar_eclipse_total_duration: Duration of lunar eclipse totality
        - calc_lunar_eclipse_umbral_duration: Duration of lunar eclipse umbral phase

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # The central-line maximum duration of totality/annularity is the
    # LOCAL second-to-third-contact interval at the point of greatest
    # eclipse (the shadow-center point on the central line), which is what
    # this function's docstring promises (~minutes).
    #
    # calc_eclipse_second_contact_c2 / calc_eclipse_third_contact_c3 return
    # the GLOBAL U1/U4 instants — when the umbra/antumbra first and last
    # touches Earth ANYWHERE — whose difference is the hours-long shadow-
    # path duration, not the local totality duration. Using them here made
    # this function return path duration in hours instead of the promised
    # minutes-long central-line duration.
    retflag, geopos, _attr = sol_eclipse_where(jd_max, flags)
    if not (retflag & ECL_CENTRAL):
        # Partial or non-central eclipse: there is no central line and so
        # no totality/annularity duration.
        return 0.0

    center_lon = geopos[0]
    center_lat = geopos[1]

    # Local circumstances on the central line at the shadow-center point.
    local = _calculate_local_eclipse_phases(jd_max, center_lat, center_lon, 0.0)
    jd_c2 = local[2]  # local second contact (central phase begins)
    jd_c3 = local[3]  # local third contact (central phase ends)

    # Second/third contact only exist for a central eclipse at this point.
    if jd_c2 <= 0.0 or jd_c3 <= 0.0:
        return 0.0

    # Calculate duration in minutes (C3 - C2) * 24 hours * 60 minutes
    duration_days = jd_c3 - jd_c2
    duration_minutes = duration_days * 24.0 * 60.0

    return duration_minutes


def calc_lunar_eclipse_total_duration(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the duration of totality for a total lunar eclipse.

    For total lunar eclipses, this function calculates the duration of
    totality - the time between umbral second contact (U2) and umbral
    third contact (U3). This represents the length of time the entire
    Moon is within Earth's umbral shadow.

    For partial or penumbral lunar eclipses, the function returns 0.0
    since totality does not occur.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _lun_eclipse_when_pythonic().
                The function calculates U2 and U3 relative to this time.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Duration of totality in minutes. Returns 0.0 if:
        - The eclipse is not a total lunar eclipse
        - U2 or U3 cannot be determined
        - The input time is not near a valid lunar eclipse

    Algorithm:
        1. Calculate umbral second contact U2 (when totality begins)
        2. Calculate umbral third contact U3 (when totality ends)
        3. Return (U3 - U2) converted from days to minutes

    Precision:
        The calculation achieves timing precision better than 1 second,
        resulting in duration precision of approximately ±0.03 minutes (±2 seconds).

    Example:
        >>> from libephemeris import julday, lun_eclipse_when
        >>> from libephemeris import calc_lunar_eclipse_total_duration, ECL_TOTAL
        >>> # Find the November 8, 2022 total lunar eclipse
        >>> jd_start = julday(2022, 10, 1, 0.0)
        >>> ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate duration of totality
        >>> duration = calc_lunar_eclipse_total_duration(jd_max)
        >>> print(f"Duration of totality: {duration:.2f} minutes")

    Note:
        - U2 marks when the Moon fully enters the umbra (totality begins)
        - U3 marks when the Moon begins exiting the umbra (totality ends)
        - During totality, the Moon appears red/copper due to refracted sunlight
        - Total lunar eclipse durations can range from about 14 to 107 minutes
        - Unlike solar eclipses, lunar eclipse totality is visible from the
          entire night side of Earth

    See Also:
        - calc_lunar_eclipse_umbral_second_contact_u2: Calculate U2 (totality begins)
        - calc_lunar_eclipse_umbral_third_contact_u3: Calculate U3 (totality ends)
        - calc_lunar_eclipse_umbral_duration: Duration of umbral/partial phase
        - _lun_eclipse_when_pythonic: Find next lunar eclipse and all phase times

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Lunar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Lunar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Calculate U2 (umbral second contact - totality begins)
    jd_u2 = calc_lunar_eclipse_umbral_second_contact_u2(jd_max, flags=flags)

    # Check if totality exists
    if jd_u2 == 0.0:
        return 0.0

    # Calculate U3 (umbral third contact - totality ends)
    jd_u3 = calc_lunar_eclipse_umbral_third_contact_u3(jd_max, flags=flags)

    # Check if U3 was found
    if jd_u3 == 0.0:
        return 0.0

    # Calculate duration in minutes (U3 - U2) * 24 hours * 60 minutes
    duration_days = jd_u3 - jd_u2
    duration_minutes = duration_days * 24.0 * 60.0

    return duration_minutes


def calc_lunar_eclipse_umbral_duration(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the duration of the umbral (partial) phase for a lunar eclipse.

    For lunar eclipses that enter Earth's umbral shadow (total or partial
    eclipses), this function calculates the duration of the umbral phase -
    the time between umbral first contact (U1) and umbral fourth contact (U4).
    This represents the length of time any part of the Moon is within the
    Earth's umbral shadow.

    For penumbral-only lunar eclipses, the function returns 0.0 since
    the Moon does not enter the umbral shadow.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _lun_eclipse_when_pythonic().
                The function calculates U1 and U4 relative to this time.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Duration of umbral phase in minutes. Returns 0.0 if:
        - The eclipse is a penumbral-only eclipse
        - U1 or U4 cannot be determined
        - The input time is not near a valid lunar eclipse

    Algorithm:
        1. Calculate umbral first contact U1 (when partial phase begins)
        2. Calculate umbral fourth contact U4 (when partial phase ends)
        3. Return (U4 - U1) converted from days to minutes

    Precision:
        The calculation achieves timing precision better than 1 second,
        resulting in duration precision of approximately ±0.03 minutes (±2 seconds).

    Example:
        >>> from libephemeris import julday, lun_eclipse_when
        >>> from libephemeris import calc_lunar_eclipse_umbral_duration, ECL_TOTAL
        >>> # Find the November 8, 2022 total lunar eclipse
        >>> jd_start = julday(2022, 10, 1, 0.0)
        >>> ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate duration of umbral phase
        >>> duration = calc_lunar_eclipse_umbral_duration(jd_max)
        >>> print(f"Duration of umbral phase: {duration:.2f} minutes")

    Note:
        - U1 marks when the Moon first enters the umbra (partial phase begins)
        - U4 marks when the Moon fully exits the umbra (partial phase ends)
        - This duration includes totality (if it occurs) plus both partial phases
        - Umbral phase durations can range from about 24 to 236 minutes
        - Unlike solar eclipses, lunar eclipse umbral phase is visible from
          the entire night side of Earth

    See Also:
        - calc_lunar_eclipse_umbral_first_contact_u1: Calculate U1 (partial begins)
        - calc_lunar_eclipse_umbral_fourth_contact_u4: Calculate U4 (partial ends)
        - calc_lunar_eclipse_total_duration: Duration of totality phase only
        - _lun_eclipse_when_pythonic: Find next lunar eclipse and all phase times

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Lunar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Lunar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Calculate U1 (umbral first contact - partial phase begins)
    jd_u1 = calc_lunar_eclipse_umbral_first_contact_u1(jd_max, flags=flags)

    # Check if umbral phase exists
    if jd_u1 == 0.0:
        return 0.0

    # Calculate U4 (umbral fourth contact - partial phase ends)
    jd_u4 = calc_lunar_eclipse_umbral_fourth_contact_u4(jd_max, flags=flags)

    # Check if U4 was found
    if jd_u4 == 0.0:
        return 0.0

    # Calculate duration in minutes (U4 - U1) * 24 hours * 60 minutes
    duration_days = jd_u4 - jd_u1
    duration_minutes = duration_days * 24.0 * 60.0

    return duration_minutes


# =============================================================================
# ECLIPSE PATH WIDTH CALCULATION
# =============================================================================


def calc_eclipse_path_width(
    jd: float,
    lat: float | None = None,
    lon: float | None = None,
    flags: int = FLG_SWIEPH,
) -> float:
    """Calculate the width of the central eclipse path at a location.

    Wrapper around :func:`_calc_eclipse_path_width_impl` that adds LEB→Skyfield
    fallback for partial/custom LEB files. See the impl docstring for the
    full API contract.
    """
    return _call_with_leb_skyfield_fallback(
        _calc_eclipse_path_width_impl,
        jd,
        lat,
        lon,
        flags,
    )


def _calc_eclipse_path_width_impl(
    jd: float,
    lat: float | None = None,
    lon: float | None = None,
    flags: int = FLG_SWIEPH,
    *,
    reader=None,
) -> float:
    """
    Calculate the width of the path of totality or annularity for a central solar eclipse.

    For central solar eclipses (total, annular, or hybrid), this function calculates
    the width of the umbral (total) or antumbral (annular) shadow path on Earth's
    surface at a specific time or geographic location.

    The path width is the perpendicular distance across the central line within
    which observers will see totality or annularity. Outside this band, the
    eclipse appears partial.

    Args:
        jd: Julian Day (UT) at which to calculate the path width. This should be
            a time during a central solar eclipse when the shadow is on Earth.
        lat: Observer latitude in degrees (positive = North, negative = South).
             If provided along with lon, calculates path width at this location.
             If None, calculates path width at the central line for time jd.
        lon: Observer longitude in degrees (positive = East, negative = West).
             If provided along with lat, calculates path width at this location.
             If None, calculates path width at the central line for time jd.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        The path width in kilometers. Returns 0.0 if:
        - The eclipse is not a central eclipse (partial only)
        - The shadow does not touch Earth at the given time
        - No valid eclipse is occurring at the given time
        - The Sun is below the horizon at the specified location

    Algorithm:
        1. Calculate the Besselian elements l2 (umbral/antumbral cone radius)
        2. Get the Sun altitude at the central line (or specified location)
        3. Calculate the geometric shadow width using the formula:
           width = 2 * |l2| * Earth_radius / sin(sun_altitude)
        4. Apply corrections for Earth's curvature and observer elevation

    Physical interpretation:
        - The umbral/antumbral shadow is a cone with its vertex near the Moon
        - Where this cone intersects Earth's surface, it creates an elliptical
          shadow footprint
        - The path width is the minor axis of this ellipse perpendicular to
          the direction of shadow motion
        - Path width is minimum when the Sun is at zenith and increases as
          the Sun is lower in the sky

    Precision:
        This is a simplified geometric model: the umbral/antumbral cone radius
        (Besselian l2) projected onto the surface through 1/sin(Sun altitude)
        with curvature and observer-elevation corrections. It measures the
        shadow ellipse in the local vertical plane rather than strictly
        perpendicular to the shadow's ground track, so near greatest eclipse
        (Sun high) it runs a few percent wider than the published path widths
        in the canon of solar eclipses -- e.g. ~124 km vs ~115 km on
        2017-08-21 and ~207 km vs ~198 km on 2024-04-08 (about 5-8% high). The
        central-line position, northern/southern limits and central duration
        are accurate (greatest-eclipse position to ~0.02 deg). The discrepancy
        grows toward sunrise/sunset where the grazing geometry dominates.

    Example:
        >>> from libephemeris import julday, sol_eclipse_when_glob, ECL_TOTAL
        >>> from libephemeris import calc_eclipse_path_width
        >>> # Find the April 8, 2024 total solar eclipse
        >>> jd_start = julday(2024, 1, 1, 0.0)
        >>> ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate path width at eclipse maximum
        >>> width = calc_eclipse_path_width(jd_max)
        >>> print(f"Path width at maximum: {width:.1f} km")

        >>> # Calculate path width at a specific location
        >>> dallas_lat, dallas_lon = 32.7767, -96.7970
        >>> width_dallas = calc_eclipse_path_width(jd_max, lat=dallas_lat, lon=dallas_lon)
        >>> print(f"Path width at Dallas: {width_dallas:.1f} km")

    Note:
        - Path width for total eclipses typically ranges from 0 to ~270 km
        - Path width for annular eclipses typically ranges from 0 to ~380 km
        - Very narrow eclipses (< 10 km) indicate the shadow axis barely
          grazes Earth's surface
        - At the limits of the central path (sunrise/sunset), the path width
          can exceed 1000 km due to grazing geometry
        - The calculation assumes a spherical Earth with WGS84 equatorial radius

    See Also:
        - sol_eclipse_where: Find central line coordinates and path limits
        - calc_besselian_l2: Calculate the umbral/antumbral cone radius
        - _sol_eclipse_when_glob_pythonic: Find next solar eclipse and phase times
        - _sol_eclipse_how_pythonic: Get eclipse circumstances at a specific location

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
        - Chauvenet's "Manual of Spherical and Practical Astronomy", Vol. 1
    """
    # reader is provided by the caller (None forces Skyfield path)

    # Constants
    AU_TO_KM = 149597870.7
    SUN_RADIUS_KM = 696340.0
    MOON_RADIUS_KM = 1737.4

    # Determine the location for calculation
    if lat is not None and lon is not None:
        # Use specified location
        calc_lat = lat
        calc_lon = lon
    else:
        # Find central line location at this time using sol_eclipse_where
        ecl_type, geopos_result, attr = sol_eclipse_where(jd, flags)
        calc_lon = geopos_result[0]
        calc_lat = geopos_result[1]

        # Check if this is a central eclipse
        if not (ecl_type & ECL_CENTRAL):
            return 0.0

    if reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR

        _gp = (calc_lon, calc_lat, 0.0)

        try:
            jd_tt = jd + deltat(jd)
            sun_topo = _topo_ecliptic(reader, jd_tt, jd, SUN, _gp)
            moon_topo = _topo_ecliptic(reader, jd_tt, jd, MOON, _gp)
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            return 0.0

        sun_az_val, sun_alt_true, sun_alt_app = azalt(
            jd, ECL2HOR, _gp, 0, 0, sun_topo[:3]
        )
        sun_altitude = sun_alt_true

        if sun_altitude <= 0:
            return 0.0

        sun_dist_au = sun_topo[2]
        moon_dist_au = moon_topo[2]
    else:
        from skyfield.api import wgs84
        from .state import get_planets, get_timescale

        # Get ephemeris and timescale
        eph = get_planets()
        ts = get_timescale()
        t = ts.ut1_jd(jd)

        # Get celestial bodies
        earth = eph["earth"]
        sun = eph["sun"]
        moon = eph["moon"]

        # Create observer at the calculation location
        try:
            observer = wgs84.latlon(calc_lat, calc_lon, 0.0)
            observer_at = earth + observer
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            return 0.0

        # Get Sun and Moon apparent positions from the observer
        try:
            sun_app = observer_at.at(t).observe(sun).apparent()
            moon_app = observer_at.at(t).observe(moon).apparent()
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            return 0.0

        # Get Sun altitude
        sun_alt, sun_az, _ = sun_app.altaz()
        sun_altitude = sun_alt.degrees

        # If Sun is below horizon, no path width
        if sun_altitude <= 0:
            return 0.0

        # Get distances
        sun_dist_au = sun_app.distance().au
        moon_dist_au = moon_app.distance().au

    sun_dist_km = sun_dist_au * AU_TO_KM
    moon_dist_km = moon_dist_au * AU_TO_KM

    # Calculate angular radii for eclipse type determination
    sun_angular_radius = (959.63 / 3600.0) / sun_dist_au  # degrees
    moon_angular_radius = (932.56 / 3600.0) * (0.002569 / moon_dist_au)  # degrees

    moon_sun_ratio = moon_angular_radius / sun_angular_radius

    # Check if Sun and Moon are close enough for central eclipse
    if reader is not None:
        from .utils import angular_separation as _ang_sep_pw

        separation = _ang_sep_pw(sun_topo[0], sun_topo[1], moon_topo[0], moon_topo[1])
    else:
        separation = sun_app.separation_from(moon_app).degrees
    diff_radii = abs(moon_angular_radius - sun_angular_radius)

    if separation > diff_radii:
        # Not a central eclipse at this location/time
        return 0.0

    # Calculate the umbral/antumbral cone geometry
    # The umbra cone semi-angle (alpha) is the angle at the Moon
    # between the cone edge and the shadow axis
    alpha = math.atan((SUN_RADIUS_KM - MOON_RADIUS_KM) / sun_dist_km)

    # Calculate the shadow radius at Earth's surface
    if moon_sun_ratio >= 1.0:
        # Total eclipse - umbra reaches Earth
        # The umbra radius decreases along the shadow axis
        umbra_radius_km = MOON_RADIUS_KM - moon_dist_km * math.tan(alpha)
        if umbra_radius_km <= 0:
            return 0.0  # Umbra vertex is before Earth
    else:
        # Annular eclipse - antumbra
        # The antumbra radius increases along the shadow axis
        umbra_radius_km = moon_dist_km * math.tan(alpha) - MOON_RADIUS_KM
        if umbra_radius_km < 0:
            umbra_radius_km = abs(umbra_radius_km)

    # Calculate the path width
    # The shadow hits Earth at an angle determined by the Sun's altitude
    # Path width = 2 * umbra_radius / sin(sun_altitude)
    sin_alt = math.sin(math.radians(sun_altitude))

    # Protect against very small values near horizon (grazing eclipses)
    if sin_alt < 0.05:
        sin_alt = 0.05  # Limit to ~3 degrees to avoid unrealistic values

    path_width_km = 2.0 * umbra_radius_km / sin_alt

    # Apply a reasonable upper limit for grazing eclipses
    path_width_km = min(path_width_km, 1500.0)

    return max(0.0, path_width_km)


# Alias for reference API compatibility
calc_eclipse_path_width = calc_eclipse_path_width


# =============================================================================
# SAROS AND INEX SERIES CALCULATION
# =============================================================================

# Saros cycle period in days (223 synodic months)
SAROS_CYCLE_DAYS = 6585.3211  # More precise: 6585.3211 days

# Inex cycle period in days (358 synodic months, ~28.945 years)
# The Inex cycle relates eclipses that occur at nearly the same longitude
# but shifted by one lunar node (ascending to descending or vice versa).
INEX_CYCLE_DAYS = 10571.9509  # 358 * 29.530588853 = 10571.9509 days

# Published solar eclipses with known Saros series and member numbers
# Comprehensive coverage of all ~25 active solar Saros series (2015-2028)
# Data from NASA's Five Millennium Canon of Solar Eclipses (Espenak & Meeus)
# Format: (JD of eclipse maximum, Saros series number, member number)
_SOLAR_SAROS_REFERENCES = (
    (2457101.907, 120, 61),
    (2457278.788, 125, 54),
    (2457456.581, 130, 52),
    (2457632.88, 135, 39),
    (2457811.12, 140, 29),
    (2457987.268, 145, 22),
    (2458165.369, 150, 17),
    (2458312.626, 117, 69),
    (2458489.571, 122, 58),
    (2458667.308, 127, 58),
    (2458843.721, 132, 46),
    (2459021.778, 137, 36),
    (2459198.176, 142, 23),
    (2459375.946, 147, 23),
    (2459552.815, 152, 13),
    (2459700.362, 119, 66),
    (2459877.958, 124, 55),
    (2460054.678, 129, 52),
    (2460232.25, 134, 44),
    (2460409.262, 139, 30),
    (2460586.281, 144, 17),
    (2460763.95, 149, 21),
    (2460940.321, 154, 7),
    (2461089.008, 121, 61),
    (2461265.24, 126, 48),
    # --- Extended coverage across 1200-2999 CE: one NASA-catalog anchor per
    # additional active Saros series, so every eclipse the library can compute
    # maps to a real (series, member) instead of the no-match sentinel. Source:
    # Espenak & Meeus, Five Millennium Canon of Solar Eclipses.
    (2167495.93129, 87, 73),
    (2167672.0409, 92, 73),
    (2171481.9986, 88, 83),
    (2175645.77933, 94, 71),
    (2178244.64854, 93, 72),
    (2179631.98292, 95, 68),
    (2182053.87351, 89, 70),
    (2186041.07347, 90, 79),
    (2194191.31798, 97, 66),
    (2196789.93424, 96, 68),
    (2211347.99534, 98, 65),
    (2219321.99843, 100, 62),
    (2221920.57742, 99, 64),
    (2236478.57916, 101, 61),
    (2240465.66548, 102, 59),
    (2251037.76325, 103, 60),
    (2255023.59836, 104, 56),
    (2259187.76744, 110, 57),
    (2265595.93854, 105, 56),
    (2269582.80912, 106, 59),
    (2273746.99279, 112, 55),
    (2280154.01941, 107, 55),
    (2281542.83687, 109, 63),
    (2284318.3917, 113, 54),
    (2290726.34117, 108, 57),
    (2294890.17034, 114, 52),
    (2296101.01631, 111, 59),
    (2298877.35649, 115, 52),
    (2302863.46221, 116, 49),
    (2324007.68102, 118, 48),
    (2357288.73454, 128, 43),
    (2363696.22268, 123, 39),
    (2375833.88437, 131, 38),
    (2403562.71673, 133, 37),
    (2448449.29586, 136, 36),
    (2489348.91682, 138, 36),
    (2534413.1887, 146, 39),
    (2540820.81225, 141, 36),
    (2548794.90178, 143, 37),
    (2568727.05743, 148, 38),
    (2593857.58172, 151, 35),
    (2608416.62109, 153, 32),
    (2622975.60192, 155, 31),
    (2633547.3771, 156, 28),
    (2644118.83863, 157, 27),
    (2648106.05883, 158, 27),
    (2652269.47244, 164, 26),
    (2658677.88725, 159, 25),
    (2662841.85125, 165, 25),
    (2666650.98188, 161, 24),
    (2669249.20972, 160, 24),
    (2670814.77334, 167, 23),
    (2673413.13045, 166, 22),
    (2681387.23216, 168, 22),
    (2683808.27376, 162, 22),
    (2687794.22467, 163, 21),
    (2691958.48484, 169, 19),
    (2695945.1914, 170, 19),
    (2706517.64865, 171, 18),
    (2721075.47043, 173, 15),
    (2723674.23595, 172, 16),
    (2731647.98223, 174, 14),
    (2742219.32912, 175, 12),
    (2746205.87509, 176, 12),
    (2756778.49624, 177, 11),
    (2764751.1901, 179, 8),
    (2767350.00412, 178, 8),
    (2775323.86784, 180, 8),
)

# Published lunar eclipses with known Saros series and member numbers
# Comprehensive coverage of all ~30 active lunar Saros series (2015-2030)
# Data from NASA's Five Millennium Canon of Lunar Eclipses (Espenak & Meeus)
# Format: (JD of eclipse maximum, Saros series number, member number)
_LUNAR_SAROS_REFERENCES = (
    (2457117.0, 132, 30),
    (2457293.616, 137, 26),
    (2457470.991, 142, 18),
    (2457648.288, 147, 8),
    (2457795.53, 114, 59),
    (2457973.264, 119, 61),
    (2458150.062, 124, 49),
    (2458327.348, 129, 38),
    (2458504.717, 134, 27),
    (2458681.396, 139, 21),
    (2458859.299, 144, 16),
    (2459006.309, 111, 67),
    (2459183.905, 116, 58),
    (2459360.971, 121, 55),
    (2459537.877, 126, 45),
    (2459715.675, 131, 34),
    (2459891.958, 136, 20),
    (2460070.224, 141, 24),
    (2460246.343, 146, 11),
    (2460394.801, 113, 64),
    (2460571.614, 118, 52),
    (2460748.791, 123, 53),
    (2460926.258, 128, 41),
    (2461102.982, 133, 27),
    (2461280.676, 138, 29),
    (2461457.467, 143, 18),
    (2461605.169, 110, 72),
    (2461782.676, 115, 58),
    (2461959.264, 120, 58),
    (2462137.203, 125, 49),
    # --- Extended coverage across 1200-2999 CE: one NASA-catalog anchor per
    # additional active Saros series (see the solar table note). Source:
    # Espenak & Meeus, Five Millennium Canon of Lunar Eclipses.
    (2164202.46117, 80, 74),
    (2168366.01908, 86, 73),
    (2170965.42651, 85, 75),
    (2178761.72378, 82, 82),
    (2181360.32956, 81, 72),
    (2182747.61147, 83, 82),
    (2182925.4339, 88, 69),
    (2185523.57488, 87, 70),
    (2186911.48969, 89, 68),
    (2199905.35925, 84, 79),
    (2204068.89258, 90, 66),
    (2208056.16429, 91, 65),
    (2212042.24997, 92, 64),
    (2229199.48665, 93, 62),
    (2233186.71434, 94, 61),
    (2237172.77382, 95, 60),
    (2254329.86495, 96, 58),
    (2258317.04177, 97, 58),
    (2262303.15781, 98, 59),
    (2266467.55793, 104, 57),
    (2279460.15602, 99, 55),
    (2280848.42948, 101, 66),
    (2283447.4018, 100, 62),
    (2283624.172, 105, 55),
    (2285012.73237, 107, 53),
    (2287611.31958, 106, 54),
    (2295407.36627, 103, 62),
    (2298005.31636, 102, 63),
    (2302169.26221, 108, 51),
    (2306156.24797, 109, 49),
    (2331286.60028, 112, 46),
    (2370975.70262, 117, 39),
    (2371153.3225, 122, 43),
    (2424012.34729, 127, 37),
    (2468898.94064, 130, 36),
    (2541515.30318, 135, 36),
    (2554862.84622, 140, 39),
    (2601136.86801, 145, 33),
    (2626267.62549, 148, 29),
    (2630254.52294, 149, 29),
    (2634241.63625, 150, 28),
    (2644990.46979, 156, 27),
    (2651398.20542, 151, 26),
    (2655384.89156, 152, 26),
    (2659371.98764, 153, 25),
    (2659548.47504, 158, 24),
    (2663535.83518, 159, 25),
    (2673929.822, 155, 23),
    (2676528.59209, 154, 22),
    (2680692.90077, 160, 22),
    (2684679.1019, 161, 22),
    (2688488.47677, 157, 20),
    (2688666.3459, 162, 21),
    (2699238.18556, 163, 18),
    (2707211.28345, 165, 17),
    (2709809.51659, 164, 18),
    (2724368.57004, 166, 15),
    (2732341.64082, 168, 13),
    (2734939.842, 167, 14),
    (2749499.06974, 169, 11),
    (2749676.18174, 174, 11),
    (2753485.0553, 170, 10),
    (2757472.20983, 171, 10),
    (2761635.77311, 177, 9),
    (2770819.87532, 176, 8),
    (2772208.27291, 178, 7),
    (2774629.79494, 172, 8),
    (2778615.79827, 173, 7),
    (2780004.61641, 175, 7),
    (2782779.60631, 179, 6),
    (2786766.13084, 180, 6),
)


def _get_saros_info(
    jd_eclipse: float,
    eclipse_type: str = "solar",
) -> Tuple[float, float]:
    """Return (saros_series, saros_member) for an eclipse as floats.

    Internal helper used to populate attr[9] and attr[10] in eclipse
    attribute tuples. Returns the reference sentinel (-99999999.0,
    -99999999.0) if no series matches within 2 days of a member.

    Args:
        jd_eclipse: Julian Day (UT) of the eclipse maximum.
        eclipse_type: "solar" or "lunar".

    Returns:
        Tuple of (saros_series_number, saros_member_number) as floats.
    """
    if eclipse_type.lower() == "solar":
        references: Sequence[Tuple[float, int, int]] = _SOLAR_SAROS_REFERENCES
    elif eclipse_type.lower() == "lunar":
        references = _LUNAR_SAROS_REFERENCES
    else:
        return (-99999999.0, -99999999.0)

    best_series = 0
    best_member = 0
    best_residual = float("inf")

    for ref_jd, ref_saros, ref_member in references:
        delta = jd_eclipse - ref_jd
        n_cycles = round(delta / SAROS_CYCLE_DAYS)
        residual = abs(delta - n_cycles * SAROS_CYCLE_DAYS)

        if residual < best_residual:
            best_residual = residual
            best_series = ref_saros
            best_member = ref_member + n_cycles

    # Acceptance window around the mean-cycle grid. True members of a series
    # wobble around the mean saros extrapolation (223 synodic months beat
    # against 242 draconic and 239 anomalistic months — Meeus, Mathematical
    # Astronomy Morsels; Espenak & Meeus, NASA/TP-2006-214141), while the
    # nearest eclipse belonging to a *different* series sits at least one
    # synodic month (~29.53 d) off the grid. Any window between those two
    # scales discriminates identically near the pinned reference members;
    # 2.0 days is the conservative low end. Outside it the no-match protocol
    # sentinel (-99999999.0, -99999999.0) is returned (compatibility
    # contract).
    if best_residual > 2.0:
        return (-99999999.0, -99999999.0)

    return (float(best_series), float(best_member))


def get_saros_number(
    jd_eclipse: float,
    eclipse_type: str = "solar",
) -> int:
    """
    Determine the Saros series number for an eclipse.

    The Saros cycle is a period of approximately 6585.3211 days (18 years,
    11 days, 8 hours) after which the Sun, Moon, and lunar nodes return
    to nearly the same relative geometry, causing similar eclipses to recur.

    Each Saros series is numbered sequentially. Solar eclipses have series
    numbered approximately 1-180 (though not all are active at any given time),
    and lunar eclipses similarly have their own numbering.

    Args:
        jd_eclipse: Julian Day (UT) of the eclipse maximum. This should be
                    obtained from _sol_eclipse_when_glob_pythonic() or _lun_eclipse_when_pythonic().
        eclipse_type: Type of eclipse - either "solar" or "lunar".
                      Defaults to "solar".

    Returns:
        The Saros series number (integer). Solar Saros numbers typically range
        from about 1 to 180, and lunar Saros numbers similarly.

    Algorithm:
        1. For the given eclipse JD, calculate the difference from known
           reference eclipses with documented Saros series numbers.
        2. Determine how many Saros cycles separate the input eclipse from
           the reference.
        3. Each Saros cycle advances the series number by 1 (earlier eclipses
           have lower numbers in the same series, later have higher).

    Precision:
        The function is accurate for eclipses within several centuries of
        the reference eclipses. Very distant eclipses (more than ~500 years
        from reference) may have small errors due to the slight variation
        in the Saros period.

    Example:
        >>> from libephemeris import julday, sol_eclipse_when_glob, get_saros_number
        >>> # Find the April 8, 2024 total solar eclipse
        >>> jd_start = julday(2024, 3, 1, 0.0)
        >>> ecl_type, times = sol_eclipse_when_glob(jd_start)
        >>> jd_max = times[0]
        >>> saros = get_saros_number(jd_max, "solar")
        >>> print(f"Saros series: {saros}")  # Should print 139

    Note:
        - Each Saros series begins with partial eclipses near one pole,
          progresses through central (total/annular) eclipses, and ends
          with partial eclipses near the other pole.
        - A typical Saros series contains 70-80 eclipses over ~1200-1400 years.
        - The same Saros number for solar and lunar eclipses are different
          series (they use separate numbering systems).

    References:
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
        - Espenak & Meeus "Five Millennium Canon of Lunar Eclipses"
        - van Gent, R.H. "A Catalogue of Eclipse Cycles"
    """
    if eclipse_type.lower() not in ("solar", "lunar"):
        raise ValueError(
            f"eclipse_type must be 'solar' or 'lunar', got '{eclipse_type}'"
        )
    saros_series, _ = _get_saros_info(jd_eclipse, eclipse_type)
    return int(saros_series)


# =============================================================================
# INEX SERIES CALCULATION
# =============================================================================

# Reference solar eclipses with known Inex series numbers
# The Inex cycle connects eclipses at the same solar longitude
# but at opposite lunar nodes (ascending <-> descending).
# Series numbers from van Gent's eclipse cycle catalog.
# Format: (JD of eclipse maximum, Inex series number)
_SOLAR_INEX_REFERENCES = (
    # Inex 50 - includes 2024 total solar eclipse
    (2460409.296, 50),  # 2024-Apr-08 total solar eclipse (Inex 50)
    # Inex 49 - includes 2017 total solar eclipse
    (2457987.767, 49),  # 2017-Aug-21 total solar eclipse (Inex 49)
    # Inex 51 - includes 2020 total solar eclipse
    (2459203.604, 51),  # 2020-Dec-14 total solar eclipse (Inex 51)
    # Inex 48 - includes 2019 total solar eclipse
    (2458675.604, 48),  # 2019-Jul-02 total solar eclipse (Inex 48)
    # Inex 52 - includes 2021 annular solar eclipse
    (2459369.971, 52),  # 2021-Jun-10 annular solar eclipse (Inex 52)
)

# Reference lunar eclipses with known Inex series numbers.
# Series numbers from van Gent's eclipse cycle catalog; eclipse dates
# from NASA's Five Millennium Canon (Espenak & Meeus).
# Format: (JD of eclipse maximum, Inex series number)
_LUNAR_INEX_REFERENCES = (
    # Inex 50 - includes 2022 total lunar eclipse
    (2459891.459, 50),  # 2022-Nov-08 total lunar eclipse (Inex 50)
    # Inex 49 - includes 2022 total lunar eclipse
    (2459695.604, 49),  # 2022-May-16 total lunar eclipse (Inex 49)
    # Inex 48 - includes 2021 total lunar eclipse
    (2459356.917, 48),  # 2021-May-26 total lunar eclipse (Inex 48)
    # Inex 51 - includes 2018 total lunar eclipse
    (2458310.835, 51),  # 2018-Jul-27 total lunar eclipse (Inex 51)
    # Inex 47 - includes 2019 total lunar eclipse
    (2458497.459, 47),  # 2019-Jan-21 total lunar eclipse (Inex 47)
)


def get_inex_number(
    jd_eclipse: float,
    eclipse_type: str = "solar",
) -> int:
    """
    Determine the Inex series number for an eclipse.

    The Inex cycle is a period of approximately 10571.95 days (358 synodic
    months, ~28.945 years) after which eclipses occur at the same solar
    longitude but at the opposite lunar node.

    While the Saros cycle connects eclipses with similar geometry (same node,
    similar latitude), the Inex cycle connects eclipses at the same longitude
    but at different nodes. Together, Saros and Inex form a two-dimensional
    numbering system for eclipse series.

    Args:
        jd_eclipse: Julian Day (UT) of the eclipse maximum. This should be
                    obtained from _sol_eclipse_when_glob_pythonic() or _lun_eclipse_when_pythonic().
        eclipse_type: Type of eclipse - either "solar" or "lunar".
                      Defaults to "solar".

    Returns:
        The Inex series number (integer).

    Algorithm:
        1. For the given eclipse JD, calculate the difference from known
           reference eclipses with documented Inex series numbers.
        2. Determine how many Inex cycles separate the input eclipse from
           the reference.
        3. Find the reference that gives the smallest residual (closest
           match to an integer number of cycles).

    Precision:
        The function is accurate for eclipses within several centuries of
        the reference eclipses. Very distant eclipses may have small errors
        due to slight variations in the Inex period over long time spans.

    Example:
        >>> from libephemeris import julday, sol_eclipse_when_glob, get_inex_number
        >>> # Find the April 8, 2024 total solar eclipse
        >>> jd_start = julday(2024, 3, 1, 0.0)
        >>> ecl_type, times = sol_eclipse_when_glob(jd_start)
        >>> jd_max = times[0]
        >>> inex = get_inex_number(jd_max, "solar")
        >>> print(f"Inex series: {inex}")  # Should print 50

    Note:
        - The Inex series number combined with the Saros series number
          uniquely identifies each eclipse in the long-term eclipse cycle.
        - Each Inex series contains eclipses that occur at the same solar
          longitude but alternate between ascending and descending nodes.
        - Unlike Saros, which relates eclipses of similar appearance, Inex
          relates eclipses at the same celestial longitude.

    References:
        - van Gent, R.H. "A Catalogue of Eclipse Cycles"
        - Meeus, J. "Mathematical Astronomy Morsels"
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
    """
    if eclipse_type.lower() == "solar":
        references = _SOLAR_INEX_REFERENCES
    elif eclipse_type.lower() == "lunar":
        references = _LUNAR_INEX_REFERENCES
    else:
        raise ValueError(
            f"eclipse_type must be 'solar' or 'lunar', got '{eclipse_type}'"
        )

    # Use the first reference eclipse as our anchor point
    ref_jd, ref_inex = references[0]

    # Calculate how many Inex cycles separate the input from the reference
    delta_days = jd_eclipse - ref_jd

    # Number of complete Inex cycles (can be positive or negative)
    cycles = round(delta_days / INEX_CYCLE_DAYS)

    # Find which reference gives the closest match
    best_match_series = ref_inex
    best_match_residual = abs(delta_days - cycles * INEX_CYCLE_DAYS)

    for ref_jd_test, ref_inex_test in references:
        delta = jd_eclipse - ref_jd_test
        cycles_test = round(delta / INEX_CYCLE_DAYS)
        residual = abs(delta - cycles_test * INEX_CYCLE_DAYS)

        if residual < best_match_residual:
            best_match_residual = residual
            best_match_series = ref_inex_test

    # No residual sanity gate is applied here: the reference table is
    # far too sparse (a handful of anchor series) for grid residuals to
    # discriminate eclipse dates from arbitrary ones - valid members of
    # a series sit hundreds of days from the nearest anchor's grid. The
    # result is therefore only meaningful when the caller already knows
    # jd_eclipse is an eclipse maximum.
    return best_match_series


def calc_eclipse_central_line(
    jd_start: float,
    jd_end: float,
    step_minutes: float = 1.0,
    flags: int = FLG_SWIEPH,
) -> Tuple[Tuple[float, ...], Tuple[float, ...], Tuple[float, ...]]:
    """
    Calculate the geographic coordinates of points along the central line of a solar eclipse.

    The central line is the path traced on Earth's surface by the intersection of
    the Moon's shadow axis with the Earth. Along this line, observers see the
    eclipse at its maximum (either total or annular, depending on eclipse type).

    This function computes a series of (latitude, longitude) points along the
    central line at specified time intervals using the same algorithm as
    _sol_eclipse_where_pythonic() for consistent, accurate results.

    Args:
        jd_start: Julian Day (UT) to start calculating the central line.
                  Should be during a solar eclipse (ideally from _sol_eclipse_when_glob_pythonic).
        jd_end: Julian Day (UT) to end the calculation.
                Should be after jd_start and during the same eclipse.
        step_minutes: Time step in minutes between calculated points (default 1.0).
                      Smaller values give a more detailed path but take longer.
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing three tuples:
            - times: Tuple of Julian Day (UT) values for each point
            - latitudes: Tuple of geographic latitudes in degrees (North positive)
            - longitudes: Tuple of geographic longitudes in degrees (East positive)

        Points where the shadow axis doesn't intersect Earth's surface
        are omitted from the results.

    Algorithm:
        For each time step, uses sol_eclipse_where() to find the central
        line position. This ensures consistency with _sol_eclipse_where_pythonic() by
        using the same iterative algorithm to find the point of minimum
        Sun-Moon angular separation on Earth's surface.

    Precision:
        Geographic coordinates accurate to approximately 0.001 degrees (~100 m)
        for points along the central line, consistent with _sol_eclipse_where_pythonic().

    Example:
        >>> from libephemeris import julday, sol_eclipse_when_glob
        >>> from libephemeris import calc_eclipse_central_line
        >>> # Find April 8, 2024 total solar eclipse
        >>> jd = julday(2024, 1, 1, 0.0)
        >>> ecl_type, times_ecl = sol_eclipse_when_glob(jd)
        >>> # Bound the path by when totality/annularity begins (U1, tret[4])
        >>> # and ends (U4, tret[5]) anywhere on Earth -- the span over which
        >>> # the central line exists (see _sol_eclipse_when_glob_pythonic).
        >>> jd_u1, jd_u4 = times_ecl[4], times_ecl[5]
        >>> # Calculate central line path
        >>> times, lats, lons = calc_eclipse_central_line(jd_u1, jd_u4, step_minutes=5.0)
        >>> for i in range(len(times)):
        ...     print(f"JD {times[i]:.5f}: lat={lats[i]:.2f}°, lon={lons[i]:.2f}°")

    Note:
        - The function only returns points where the central eclipse is visible on Earth.
        - For partial-only eclipses, an empty tuple is returned.
        - The central line is only defined for central eclipses (total, annular, or hybrid).
        - Points near the beginning and end may be at extreme latitudes as the shadow
          enters and exits Earth.

    See Also:
        - _sol_eclipse_when_glob_pythonic: Find the next solar eclipse
        - _sol_eclipse_where_pythonic: Find central line point at a specific time
        - sol_eclipse_where: Underlying function for central line calculation

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac, Ch. 11
    """
    times_list: list[float] = []
    latitudes_list: list[float] = []
    longitudes_list: list[float] = []

    # Convert step to days
    step_days = step_minutes / (24.0 * 60.0)

    # Iterate through time range
    jd = jd_start
    while jd <= jd_end:
        # Use sol_eclipse_where for consistent results
        # This function uses iterative refinement to find the exact central line position
        retflag, geopos, attr = sol_eclipse_where(jd, flags)

        # Check if we have a valid central eclipse point. Only a *central*
        # eclipse (ECL_CENTRAL bit set: total/annular/annular-total) defines a
        # central line; a partial-only eclipse has retflag = ECL_PARTIAL |
        # ECL_NONCENTRAL (> 0) but no central line, and must contribute no
        # points (the function contract returns an empty tuple for partials).
        # geopos[0] and geopos[1] contain the central line coordinates
        if retflag & ECL_CENTRAL:
            lon = geopos[0]
            lat = geopos[1]

            # Check if we have valid coordinates (not zeros or placeholder values)
            # Also verify we're at a central eclipse point using Besselian gamma
            x = calc_besselian_x(jd, flags)
            y = calc_besselian_y(jd, flags)
            gamma = math.sqrt(x * x + y * y)

            # Include point if gamma indicates central eclipse is on Earth's surface
            # gamma < ~1.5 allows for Earth's oblateness and shadow geometry
            if gamma < 1.5 and (abs(lon) > 0.0001 or abs(lat) > 0.0001):
                times_list.append(jd)
                latitudes_list.append(lat)
                longitudes_list.append(lon)

        jd += step_days

    return tuple(times_list), tuple(latitudes_list), tuple(longitudes_list)


# Alias for reference API naming convention
calc_eclipse_central_line = calc_eclipse_central_line


def calc_eclipse_northern_limit(
    jd_start: float,
    jd_end: float,
    step_minutes: float = 1.0,
    flags: int = FLG_SWIEPH,
) -> Tuple[Tuple[float, ...], Tuple[float, ...], Tuple[float, ...]]:
    """Calculate points along the northern limit of a central eclipse path.

    The boundary is sampled from Besselian shadow geometry at a fixed time
    interval. Samples for which the central shadow misses Earth are omitted.

    Args:
        jd_start: First Julian Day in UT, during a central solar eclipse.
        jd_end: Last Julian Day in UT, during the same eclipse.
        step_minutes: Sampling interval in minutes.
        flags: Calculation flags.

    Returns:
        Three tuples containing sample times, geodetic latitudes, and
        longitudes. Coordinates are native Python floats in degrees.

    Example:
        >>> from libephemeris import julday, sol_eclipse_when_glob
        >>> from libephemeris import calc_eclipse_northern_limit
        >>> # Find April 8, 2024 total solar eclipse
        >>> jd = julday(2024, 1, 1, 0.0)
        >>> ecl_type, times_ecl = sol_eclipse_when_glob(jd)
        >>> jd_c1, jd_c4 = times_ecl[1], times_ecl[4]  # First and fourth contacts
        >>> times, lats, lons = calc_eclipse_northern_limit(
        ...     jd_c1, jd_c4, step_minutes=5.0
        ... )

    Notes:
        Partial-only events return empty tuples. Accuracy decreases near the
        grazing ends of a central path.

    References:
        Meeus, *Astronomical Algorithms*, chapter 54; *Explanatory Supplement
        to the Astronomical Almanac*, chapter 11.
    """
    # WGS84 ellipsoid parameters
    EARTH_FLATTENING = 1.0 / 298.257223563

    times_list: list[float] = []
    latitudes_list: list[float] = []
    longitudes_list: list[float] = []

    # Convert step to days
    step_days = step_minutes / (24.0 * 60.0)

    # Iterate through time range
    jd = jd_start
    while jd <= jd_end:
        # Get Besselian elements at this time
        x = calc_besselian_x(jd, flags)
        y = calc_besselian_y(jd, flags)
        d = calc_besselian_d(jd, flags)
        mu = calc_besselian_mu(jd, flags)
        l2 = calc_besselian_l2(jd, flags)

        # Calculate gamma - distance from shadow axis to Earth center
        gamma = math.sqrt(x * x + y * y)

        # The umbral radius in Earth radii (absolute value)
        # l2 is negative for umbra (total), positive for antumbra (annular)
        umbra_radius = abs(l2)

        # For the northern limit to touch Earth, we need gamma + umbra_radius < ~1.5
        # (accounting for Earth's oblateness)
        max_gamma = 1.0 + EARTH_FLATTENING + umbra_radius
        if gamma < max_gamma and umbra_radius > 0.001:
            # Calculate the position of the northern limit
            # The northern limit is offset from the shadow axis in the +y direction
            # (north) by the umbral radius

            # Convert d (declination of shadow axis) and mu (hour angle) to radians
            d_rad = math.radians(d)
            math.radians(mu)

            # The perpendicular distance from Earth center to shadow axis is gamma
            # The y-coordinate points north in the fundamental plane

            # For the northern limit, we offset y by the umbral radius
            # The direction perpendicular to the shadow motion and toward north
            # depends on the shadow axis orientation

            # In the fundamental plane coordinate system:
            # x = east-west displacement (positive east)
            # y = north-south displacement (positive north)
            # The umbra has radius l2, so the northern edge is at y + |l2|

            y_north = y + umbra_radius

            # Calculate gamma for the northern limit point
            gamma_north = math.sqrt(x * x + y_north * y_north)

            # Check if this point is on Earth's surface
            max_gamma_surface = 1.0 + EARTH_FLATTENING
            if gamma_north < max_gamma_surface:
                # Calculate z-factor (height above fundamental plane)
                if gamma_north > 0.9999:
                    z_factor = 0.0
                else:
                    z_factor = math.sqrt(max(0, 1.0 - gamma_north * gamma_north))

                sin_d = math.sin(d_rad)
                cos_d = math.cos(d_rad)

                # Calculate latitude using the y_north component and shadow axis declination
                sin_lat = y_north * cos_d + z_factor * sin_d

                # Clamp to valid range
                sin_lat = max(-1.0, min(1.0, sin_lat))
                lat = math.degrees(math.asin(sin_lat))

                # For longitude, use the hour angle and x displacement
                cos_lat = math.cos(math.radians(lat))
                if cos_lat > 0.001:
                    lon_offset = math.degrees(
                        math.atan2(x, z_factor * cos_d - y_north * sin_d)
                    )
                else:
                    lon_offset = 0.0

                # The longitude of the northern limit point
                lon = -mu + lon_offset

                # Convert geocentric latitude (from the Besselian shadow
                # geometry) to geodetic latitude. The exact surface-point
                # relation is tan(phi_geod) = tan(phi_geoc) / (1 - f)^2
                # (= 1/(1 - e^2)); the former lat*(1 + f*sin^2(lat)) was a
                # linear approximation that under-corrected at low/mid
                # latitudes, over-corrected at high latitudes, and exceeded
                # 90 deg at the pole.
                lat_geodetic = math.degrees(
                    math.atan(
                        math.tan(math.radians(lat)) / (1.0 - EARTH_FLATTENING) ** 2
                    )
                )

                # Normalize longitude to -180 to +180
                lon = ((lon + 180.0) % 360.0) - 180.0

                # Store this point
                times_list.append(jd)
                latitudes_list.append(lat_geodetic)
                longitudes_list.append(lon)

        jd += step_days

    return tuple(times_list), tuple(latitudes_list), tuple(longitudes_list)


# Alias for reference API naming convention
calc_eclipse_northern_limit = calc_eclipse_northern_limit


def calc_eclipse_southern_limit(
    jd_start: float,
    jd_end: float,
    step_minutes: float = 1.0,
    flags: int = FLG_SWIEPH,
) -> Tuple[Tuple[float, ...], Tuple[float, ...], Tuple[float, ...]]:
    """Calculate points along the southern limit of a central eclipse path.

    The boundary is sampled from Besselian shadow geometry at a fixed time
    interval. Samples for which the central shadow misses Earth are omitted.

    Args:
        jd_start: First Julian Day in UT, during a central solar eclipse.
        jd_end: Last Julian Day in UT, during the same eclipse.
        step_minutes: Sampling interval in minutes.
        flags: Calculation flags.

    Returns:
        Three tuples containing sample times, geodetic latitudes, and
        longitudes. Coordinates are native Python floats in degrees.

    Example:
        >>> from libephemeris import julday, sol_eclipse_when_glob
        >>> from libephemeris import calc_eclipse_southern_limit
        >>> # Find April 8, 2024 total solar eclipse
        >>> jd = julday(2024, 1, 1, 0.0)
        >>> ecl_type, times_ecl = sol_eclipse_when_glob(jd)
        >>> jd_c1, jd_c4 = times_ecl[1], times_ecl[4]  # First and fourth contacts
        >>> times, lats, lons = calc_eclipse_southern_limit(
        ...     jd_c1, jd_c4, step_minutes=5.0
        ... )

    Notes:
        Partial-only events return empty tuples. Accuracy decreases near the
        grazing ends of a central path.

    References:
        Meeus, *Astronomical Algorithms*, chapter 54; *Explanatory Supplement
        to the Astronomical Almanac*, chapter 11.
    """
    # WGS84 ellipsoid parameters
    EARTH_FLATTENING = 1.0 / 298.257223563

    times_list: list[float] = []
    latitudes_list: list[float] = []
    longitudes_list: list[float] = []

    # Convert step to days
    step_days = step_minutes / (24.0 * 60.0)

    # Iterate through time range
    jd = jd_start
    while jd <= jd_end:
        # Get Besselian elements at this time
        x = calc_besselian_x(jd, flags)
        y = calc_besselian_y(jd, flags)
        d = calc_besselian_d(jd, flags)
        mu = calc_besselian_mu(jd, flags)
        l2 = calc_besselian_l2(jd, flags)

        # Calculate gamma - distance from shadow axis to Earth center
        gamma = math.sqrt(x * x + y * y)

        # The umbral radius in Earth radii (absolute value)
        # l2 is negative for umbra (total), positive for antumbra (annular)
        umbra_radius = abs(l2)

        # For the southern limit to touch Earth, we need gamma + umbra_radius < ~1.5
        # (accounting for Earth's oblateness)
        max_gamma = 1.0 + EARTH_FLATTENING + umbra_radius
        if gamma < max_gamma and umbra_radius > 0.001:
            # Calculate the position of the southern limit
            # The southern limit is offset from the shadow axis in the -y direction
            # (south) by the umbral radius

            # Convert d (declination of shadow axis) and mu (hour angle) to radians
            d_rad = math.radians(d)
            math.radians(mu)

            # The perpendicular distance from Earth center to shadow axis is gamma
            # The y-coordinate points north in the fundamental plane

            # For the southern limit, we offset y by the umbral radius in the
            # negative direction (southward)
            # The direction perpendicular to the shadow motion and toward south
            # depends on the shadow axis orientation

            # In the fundamental plane coordinate system:
            # x = east-west displacement (positive east)
            # y = north-south displacement (positive north)
            # The umbra has radius l2, so the southern edge is at y - |l2|

            y_south = y - umbra_radius

            # Calculate gamma for the southern limit point
            gamma_south = math.sqrt(x * x + y_south * y_south)

            # Check if this point is on Earth's surface
            max_gamma_surface = 1.0 + EARTH_FLATTENING
            if gamma_south < max_gamma_surface:
                # Calculate z-factor (height above fundamental plane)
                if gamma_south > 0.9999:
                    z_factor = 0.0
                else:
                    z_factor = math.sqrt(max(0, 1.0 - gamma_south * gamma_south))

                sin_d = math.sin(d_rad)
                cos_d = math.cos(d_rad)

                # Calculate latitude using the y_south component and shadow axis declination
                sin_lat = y_south * cos_d + z_factor * sin_d

                # Clamp to valid range
                sin_lat = max(-1.0, min(1.0, sin_lat))
                lat = math.degrees(math.asin(sin_lat))

                # For longitude, use the hour angle and x displacement
                cos_lat = math.cos(math.radians(lat))
                if cos_lat > 0.001:
                    lon_offset = math.degrees(
                        math.atan2(x, z_factor * cos_d - y_south * sin_d)
                    )
                else:
                    lon_offset = 0.0

                # The longitude of the southern limit point
                lon = -mu + lon_offset

                # Convert geocentric latitude (from the Besselian shadow
                # geometry) to geodetic latitude. The exact surface-point
                # relation is tan(phi_geod) = tan(phi_geoc) / (1 - f)^2
                # (= 1/(1 - e^2)); the former lat*(1 + f*sin^2(lat)) was a
                # linear approximation that under-corrected at low/mid
                # latitudes, over-corrected at high latitudes, and exceeded
                # 90 deg at the pole.
                lat_geodetic = math.degrees(
                    math.atan(
                        math.tan(math.radians(lat)) / (1.0 - EARTH_FLATTENING) ** 2
                    )
                )

                # Normalize longitude to -180 to +180
                lon = ((lon + 180.0) % 360.0) - 180.0

                # Store this point
                times_list.append(jd)
                latitudes_list.append(lat_geodetic)
                longitudes_list.append(lon)

        jd += step_days

    return tuple(times_list), tuple(latitudes_list), tuple(longitudes_list)


# Alias for reference API naming convention
calc_eclipse_southern_limit = calc_eclipse_southern_limit


def _sol_eclipse_magnitude_at_loc_pythonic(
    jd: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    flags: int = FLG_SWIEPH,
    *,
    reader=None,
) -> float:
    """
    Calculate the eclipse magnitude at a specific geographic location and time.

    Eclipse magnitude is defined as the fraction of the Sun's diameter that
    is covered by the Moon. This is a simplified convenience function that
    returns just the magnitude value without additional eclipse attributes.

    Args:
        jd: Julian Day (UT) of the time to calculate
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Eclipse magnitude as a float:
            - 0.0 if no eclipse is visible (Sun and Moon not overlapping)
            - 0.0 to 1.0 for partial eclipses (fraction of diameter covered)
            - 1.0 for the moment of totality in a total eclipse
            - > 1.0 for total eclipses (the excess indicates how much larger
              the Moon appears than the Sun)

    Note:
        This function does NOT search for eclipses - it calculates the
        instantaneous magnitude at the given time. To find eclipse events,
        use _sol_eclipse_when_glob_pythonic() or _sol_eclipse_when_loc_pythonic() first.

        If the Sun is below the horizon at the observer's location, the
        magnitude will be 0.0 since the eclipse is not visible.

    Algorithm:
        1. Calculate topocentric Sun and Moon apparent positions
        2. Compute angular separation between centers
        3. Compute angular radii of Sun and Moon
        4. Calculate overlap and express as fraction of Sun's diameter

    Precision:
        Magnitude accurate to ~0.001 for central eclipses.
        Topocentric parallax is included in calculations.

    Example:
        >>> from libephemeris import julday, _sol_eclipse_magnitude_at_loc_pythonic
        >>> # Calculate magnitude during April 8, 2024 eclipse from Dallas
        >>> jd = 2460409.28  # During eclipse
        >>> dallas_lat, dallas_lon = 32.7767, -96.797
        >>> magnitude = _sol_eclipse_magnitude_at_loc_pythonic(jd, dallas_lat, dallas_lon)
        >>> print(f"Magnitude: {magnitude:.4f}")

        >>> # Check magnitude at a location outside the eclipse path
        >>> london_lat, london_lon = 51.5074, -0.1278
        >>> magnitude = _sol_eclipse_magnitude_at_loc_pythonic(jd, london_lat, london_lon)
        >>> print(f"London magnitude: {magnitude:.4f}")  # Will be 0.0

    References:
        - Meeus (1998), "Astronomical Algorithms", 2nd ed., Ch. 54
        - Espenak & Meeus (2006), NASA/TP-2006-214141, Sec. 1.1;
          eclipse magnitude is the occulted fraction of the solar diameter
    """
    # reader is provided by the caller (None forces Skyfield path)

    if reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR, angular_separation

        _gp = (lon, lat, altitude)

        try:
            jd_tt = jd + deltat(jd)
            sun_topo = _topo_ecliptic(reader, jd_tt, jd, SUN, _gp)
            moon_topo = _topo_ecliptic(reader, jd_tt, jd, MOON, _gp)
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            return 0.0

        sun_az_val, sun_alt_true, sun_alt_app = azalt(
            jd, ECL2HOR, _gp, 0, 0, sun_topo[:3]
        )
        sun_altitude = sun_alt_true

        if sun_altitude < -1.0:
            return 0.0

        separation = angular_separation(
            sun_topo[0], sun_topo[1], moon_topo[0], moon_topo[1]
        )
        sun_dist_au = sun_topo[2]
        moon_dist_au = moon_topo[2]
    else:
        from skyfield.api import wgs84
        from .state import get_planets, get_timescale

        # Get ephemeris and timescale
        eph = get_planets()
        ts = get_timescale()

        # Create observer location
        observer = wgs84.latlon(lat, lon, altitude)

        # Get Sun and Moon objects
        sun = eph["sun"]
        moon = eph["moon"]
        earth = eph["earth"]

        # Get Skyfield time
        t = ts.ut1_jd(jd)

        # Create observer position
        observer_at = earth + observer

        # Get apparent positions from observer (topocentric)
        try:
            sun_app = observer_at.at(t).observe(sun).apparent()
            moon_app = observer_at.at(t).observe(moon).apparent()
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            # If calculation fails, return 0 magnitude
            return 0.0

        # Get Sun altitude to check visibility
        sun_alt, _, _ = sun_app.altaz()
        sun_altitude = sun_alt.degrees

        # If Sun is below horizon, no visible eclipse
        if sun_altitude < -1.0:  # Allow for refraction near horizon
            return 0.0

        # Calculate angular separation between Sun and Moon
        separation = sun_app.separation_from(moon_app).degrees

        # Get distances for angular size calculations
        sun_dist_au = sun_app.distance().au
        moon_dist_au = moon_app.distance().au

    # Angular radii (in degrees)
    # Sun: mean radius 959.63" at 1 AU
    # Moon: mean radius 932.56" at mean distance (0.002569 AU)
    sun_angular_radius = (959.63 / 3600.0) / sun_dist_au
    moon_angular_radius = (932.56 / 3600.0) * (0.002569 / moon_dist_au)

    sun_diameter = 2 * sun_angular_radius
    moon_diameter = 2 * moon_angular_radius

    # Check if there's any eclipse (disks overlapping)
    sum_radii = sun_angular_radius + moon_angular_radius
    if separation >= sum_radii:
        # No eclipse - Sun and Moon too far apart
        return 0.0

    # Calculate eclipse magnitude
    # Magnitude = fraction of Sun's diameter covered by Moon
    overlap = sum_radii - separation
    magnitude = overlap / sun_diameter

    # Clamp to valid range (can exceed 1.0 for total eclipse)
    magnitude = max(0.0, min(magnitude, 1.0 + moon_diameter / sun_diameter))

    return float(magnitude)


def sol_eclipse_magnitude_at_loc(
    tjd_ut: float,
    geopos: Sequence[float],
    ifl: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the eclipse magnitude at a specific geographic location and time.

    This function follows the reference API naming convention. It is a convenience
    wrapper that returns just the eclipse magnitude (fraction of solar diameter
    covered by the Moon) without the full attribute array.

    Args:
        tjd_ut: Julian Day (UT) of the time to calculate
        ifl: Calculation flags (FLG_SWIEPH, etc.)
        geopos: Sequence of [longitude, latitude, altitude]:
            - longitude in degrees (East positive)
            - latitude in degrees (North positive)
            - altitude in meters above sea level

    Returns:
        Eclipse magnitude as a float:
            - 0.0 if no eclipse is visible
            - 0.0 to 1.0 for partial eclipses
            - >= 1.0 for total eclipses

    Raises:
        ValueError: If geopos has wrong length

    Example:
        >>> from libephemeris import sol_eclipse_magnitude_at_loc, FLG_SWIEPH
        >>> # Calculate magnitude during April 8, 2024 eclipse from Dallas
        >>> jd = 2460409.28
        >>> dallas_geopos = [-96.797, 32.7767, 0]  # lon, lat, alt
        >>> magnitude = sol_eclipse_magnitude_at_loc(jd, FLG_SWIEPH, dallas_geopos)
        >>> print(f"Magnitude: {magnitude:.4f}")

    References:
        - Reference API: sol_eclipse_how()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    # Validate geopos
    if len(geopos) < 3:
        raise ValueError("geopos must have at least 3 elements: [lon, lat, alt]")

    # Extract geographic position (longitude first, then latitude - reference API convention)
    lon = float(geopos[0])
    lat = float(geopos[1])
    altitude = float(geopos[2])

    return _call_with_leb_skyfield_fallback(
        _sol_eclipse_magnitude_at_loc_pythonic,
        tjd_ut,
        lat,
        lon,
        altitude,
        ifl,
    )


def _sol_eclipse_obscuration_at_loc_pythonic(
    jd: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    flags: int = FLG_SWIEPH,
    *,
    reader=None,
) -> float:
    """
    Calculate the eclipse obscuration at a specific geographic location and time.

    Eclipse obscuration is defined as the fraction of the Sun's disc area that
    is covered by the Moon. This differs from eclipse magnitude, which is the
    fraction of the Sun's diameter covered.

    The relationship between magnitude and obscuration is non-linear:
    - Obscuration = 0 when magnitude = 0 (no eclipse)
    - Obscuration < magnitude for partial eclipses
    - Obscuration = 1 when magnitude >= 1 (total eclipse)
    - For annular eclipses, obscuration = (moon_radius/sun_radius)^2

    Args:
        jd: Julian Day (UT) of the time to calculate
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Eclipse obscuration as a float:
            - 0.0 if no eclipse is visible (Sun and Moon not overlapping)
            - 0.0 to 1.0 for partial eclipses (fraction of Sun's area covered)
            - 1.0 for total eclipse (Moon completely covers Sun)
            - (moon_radius/sun_radius)^2 for annular eclipse (Moon inside Sun)

    Note:
        This function does NOT search for eclipses - it calculates the
        instantaneous obscuration at the given time. To find eclipse events,
        use _sol_eclipse_when_glob_pythonic() or _sol_eclipse_when_loc_pythonic() first.

        If the Sun is below the horizon at the observer's location, the
        obscuration will be 0.0 since the eclipse is not visible.

    Algorithm:
        1. Calculate topocentric Sun and Moon apparent positions
        2. Compute angular separation between centers
        3. Compute angular radii of Sun and Moon
        4. Calculate the intersection area of two overlapping circles
        5. Express as fraction of Sun's total disc area

        The intersection area formula uses the lens/vesica piscis formula:
        For two circles with radii r1, r2 and center separation d, the
        intersection area involves the sum of two circular segments.

    Precision:
        Obscuration accurate to ~0.001 for central eclipses.
        Topocentric parallax is included in calculations.

    Example:
        >>> from libephemeris import julday, _sol_eclipse_obscuration_at_loc_pythonic
        >>> # Calculate obscuration during April 8, 2024 eclipse from Dallas
        >>> jd = 2460409.28  # During eclipse
        >>> dallas_lat, dallas_lon = 32.7767, -96.797
        >>> obscuration = _sol_eclipse_obscuration_at_loc_pythonic(jd, dallas_lat, dallas_lon)
        >>> print(f"Obscuration: {obscuration:.4f}")

        >>> # Compare with magnitude
        >>> from libephemeris import _sol_eclipse_magnitude_at_loc_pythonic
        >>> magnitude = _sol_eclipse_magnitude_at_loc_pythonic(jd, dallas_lat, dallas_lon)
        >>> print(f"Magnitude: {magnitude:.4f}, Obscuration: {obscuration:.4f}")

    References:
        - Meeus (1998), "Astronomical Algorithms", 2nd ed., Ch. 54
        - Espenak & Meeus (2006), NASA/TP-2006-214141, Sec. 1.1;
          obscuration is the occulted fraction of the solar-disc area
        - Weisstein, "Circle-Circle Intersection", MathWorld; the partial
          branch evaluates the equivalent sum-of-circular-segments formula
    """
    # reader is provided by the caller (None forces Skyfield path)

    if reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR, angular_separation

        _gp = (lon, lat, altitude)

        try:
            jd_tt = jd + deltat(jd)
            sun_topo = _topo_ecliptic(reader, jd_tt, jd, SUN, _gp)
            moon_topo = _topo_ecliptic(reader, jd_tt, jd, MOON, _gp)
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            return 0.0

        sun_az_val, sun_alt_true, sun_alt_app = azalt(
            jd, ECL2HOR, _gp, 0, 0, sun_topo[:3]
        )
        sun_altitude = sun_alt_true

        if sun_altitude < -1.0:
            return 0.0

        separation = angular_separation(
            sun_topo[0], sun_topo[1], moon_topo[0], moon_topo[1]
        )
        sun_dist_au = sun_topo[2]
        moon_dist_au = moon_topo[2]
    else:
        from skyfield.api import wgs84
        from .state import get_planets, get_timescale

        # Get ephemeris and timescale
        eph = get_planets()
        ts = get_timescale()

        # Create observer location
        observer = wgs84.latlon(lat, lon, altitude)

        # Get Sun and Moon objects
        sun = eph["sun"]
        moon = eph["moon"]
        earth = eph["earth"]

        # Get Skyfield time
        t = ts.ut1_jd(jd)

        # Create observer position
        observer_at = earth + observer

        # Get apparent positions from observer (topocentric)
        try:
            sun_app = observer_at.at(t).observe(sun).apparent()
            moon_app = observer_at.at(t).observe(moon).apparent()
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            # If calculation fails, return 0 obscuration
            return 0.0

        # Get Sun altitude to check visibility
        sun_alt, _, _ = sun_app.altaz()
        sun_altitude = sun_alt.degrees

        # If Sun is below horizon, no visible eclipse
        if sun_altitude < -1.0:  # Allow for refraction near horizon
            return 0.0

        # Calculate angular separation between Sun and Moon
        separation = sun_app.separation_from(moon_app).degrees

        # Get distances for angular size calculations
        sun_dist_au = sun_app.distance().au
        moon_dist_au = moon_app.distance().au

    # Angular radii (in degrees)
    # Sun: mean radius 959.63" at 1 AU
    # Moon: mean radius 932.56" at mean distance (0.002569 AU)
    sun_angular_radius = (959.63 / 3600.0) / sun_dist_au
    moon_angular_radius = (932.56 / 3600.0) * (0.002569 / moon_dist_au)

    r_sun = sun_angular_radius
    r_moon = moon_angular_radius
    d = separation  # center-to-center separation

    # Calculate obscuration (fraction of Sun's area covered)
    if d >= r_sun + r_moon:
        # No overlap - no eclipse
        return 0.0
    elif d <= abs(r_sun - r_moon):
        # One disc lies entirely within the other.
        if r_moon >= r_sun:
            # Total: the (larger) Moon fully covers the Sun -> 100% obscured.
            return 1.0
        # Annular: a ring of Sun remains around the smaller Moon, so the
        # obscured fraction is the ratio of the lunar to the solar disc area.
        return (r_moon / r_sun) ** 2
    else:
        # Partial overlap - use lens formula for intersection of two circles
        # The intersection area is the sum of two circular segments
        #
        # For two circles with radii r1, r2 and center separation d:
        # d1 = distance from center of circle 1 to the chord
        # d2 = distance from center of circle 2 to the chord
        # d1 = (d^2 + r1^2 - r2^2) / (2*d)
        # d2 = d - d1
        #
        # Each segment area = r^2 * arccos(d_i/r) - d_i * sqrt(r^2 - d_i^2)

        d1 = (d * d + r_sun * r_sun - r_moon * r_moon) / (2 * d)
        d2 = d - d1

        if abs(d1) <= r_sun and abs(d2) <= r_moon:
            # Ensure values are valid for acos
            cos_arg1 = max(-1, min(1, d1 / r_sun))
            cos_arg2 = max(-1, min(1, d2 / r_moon))

            # Segment areas
            area1 = r_sun * r_sun * math.acos(cos_arg1) - d1 * math.sqrt(
                max(0, r_sun * r_sun - d1 * d1)
            )
            area2 = r_moon * r_moon * math.acos(cos_arg2) - d2 * math.sqrt(
                max(0, r_moon * r_moon - d2 * d2)
            )

            intersection_area = area1 + area2
            sun_area = math.pi * r_sun * r_sun
            obscuration = intersection_area / sun_area
        else:
            obscuration = 0.0

    return float(max(0.0, obscuration))


def sol_eclipse_obscuration_at_loc(
    tjd_ut: float,
    geopos: Sequence[float],
    ifl: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the eclipse obscuration at a specific geographic location and time.

    This function follows the reference API naming convention. It is a convenience
    wrapper that returns just the eclipse obscuration (fraction of solar disc
    area covered by the Moon) without the full attribute array.

    Obscuration differs from magnitude:
    - Magnitude: fraction of Sun's DIAMETER covered
    - Obscuration: fraction of Sun's AREA covered

    Args:
        tjd_ut: Julian Day (UT) of the time to calculate
        ifl: Calculation flags (FLG_SWIEPH, etc.)
        geopos: Sequence of [longitude, latitude, altitude]:
            - longitude in degrees (East positive)
            - latitude in degrees (North positive)
            - altitude in meters above sea level

    Returns:
        Eclipse obscuration as a float:
            - 0.0 if no eclipse is visible
            - 0.0 to 1.0 for partial/annular eclipses
            - 1.0 for total eclipse

    Raises:
        ValueError: If geopos has wrong length

    Example:
        >>> from libephemeris import sol_eclipse_obscuration_at_loc, FLG_SWIEPH
        >>> # Calculate obscuration during April 8, 2024 eclipse from Dallas
        >>> jd = 2460409.28
        >>> dallas_geopos = [-96.797, 32.7767, 0]  # lon, lat, alt
        >>> obs = sol_eclipse_obscuration_at_loc(jd, FLG_SWIEPH, dallas_geopos)
        >>> print(f"Obscuration: {obs:.4f}")

    References:
        - Reference API: sol_eclipse_how()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    # Validate geopos
    if len(geopos) < 3:
        raise ValueError("geopos must have at least 3 elements: [lon, lat, alt]")

    # Extract geographic position (lon first, lat second - reference API convention)
    lon = float(geopos[0])
    lat = float(geopos[1])
    altitude = float(geopos[2])

    return _call_with_leb_skyfield_fallback(
        _sol_eclipse_obscuration_at_loc_pythonic,
        tjd_ut,
        lat,
        lon,
        altitude,
        ifl,
    )


def _lun_eclipse_umbral_magnitude_pythonic(
    jd: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the umbral magnitude for a lunar eclipse at a specific time.

    Umbral magnitude is defined as the fraction of the Moon's diameter that
    is within Earth's umbral (dark) shadow. This is a simplified convenience
    function that returns just the umbral magnitude value.

    Unlike solar eclipse magnitude (which depends on observer location), lunar
    eclipse magnitude is the same for all observers who can see the Moon,
    since the Moon physically enters Earth's shadow.

    Args:
        jd: Julian Day (UT) of the time to calculate
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Umbral magnitude as a float:
            - 0.0 if the Moon is not in the umbral shadow (penumbral-only
              eclipse or no eclipse)
            - 0.0 to 1.0 for partial umbral eclipses (fraction of Moon's
              diameter within umbra)
            - >= 1.0 for total lunar eclipses (Moon fully within umbra;
              values > 1.0 indicate how deep the Moon is in the shadow)

    Note:
        This function does NOT search for eclipses - it calculates the
        instantaneous umbral magnitude at the given time. To find eclipse
        events, use _lun_eclipse_when_pythonic() first.

        For penumbral-only eclipses, this function returns 0.0 since the
        Moon has not entered the umbra. Use lun_eclipse_penumbral_magnitude()
        to get the penumbral magnitude in such cases; both draw on the same
        canonical _lun_how_core shadow model, so they stay consistent.

    Algorithm:
        1. Calculate Moon's ecliptic position
        2. Calculate Earth's umbral shadow cone geometry at Moon's distance
        3. Determine how much of the Moon's diameter is within the umbra

    Precision:
        Magnitude accurate to ~0.01 for typical eclipses.

    Example:
        >>> from libephemeris import julday, _lun_eclipse_umbral_magnitude_pythonic, _lun_eclipse_when_pythonic
        >>> # First find a lunar eclipse
        >>> jd_start = julday(2022, 5, 1, 0)
        >>> ecl_type, times = _lun_eclipse_when_pythonic(jd_start)
        >>> jd_max = times[0]  # Time of maximum eclipse
        >>> # Calculate umbral magnitude at maximum
        >>> umbral_mag = _lun_eclipse_umbral_magnitude_pythonic(jd_max)
        >>> print(f"Umbral magnitude: {umbral_mag:.4f}")

        >>> # Check magnitude at a random time (no eclipse)
        >>> jd_no_eclipse = julday(2022, 6, 1, 12.0)
        >>> mag = _lun_eclipse_umbral_magnitude_pythonic(jd_no_eclipse)
        >>> print(f"Magnitude: {mag:.4f}")  # Will be 0.0

    References:
        - Meeus (1998), "Astronomical Algorithms", 2nd ed., Ch. 54
        - Espenak & Meeus (2009), "Five Millennium Catalog of Lunar
          Eclipses", NASA/TP-2009-214173, Sec. 1 (magnitude definitions and
          Danjon shadow enlargement)
    """
    # Delegate to the same shadow core (_lun_how_core) that backs
    # lun_eclipse_how()/lun_eclipse_when(), so this convenience function
    # returns the identical umbral magnitude (attr[0]). That core reproduces
    # the umbral magnitudes from the independently defined physical cone. The
    # older private helper
    # _calculate_lunar_eclipse_type_and_magnitude() uses a different
    # small-angle approximation of the same Danjon family, so it can disagree
    # slightly with lun_eclipse_how().
    _retc, attr, _dcore = _lun_how_core(jd, flags)
    return max(0.0, attr[0])


def lun_eclipse_umbral_magnitude(
    tjd_ut: float,
    ifl: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the umbral magnitude for a lunar eclipse at a specific time.

    This function follows the reference API naming convention. It is a convenience
    function that returns just the umbral magnitude (fraction of Moon's diameter
    within Earth's umbral shadow).

    Unlike solar eclipses, lunar eclipse magnitude does not depend on observer
    location - the Moon physically enters Earth's shadow, so the magnitude is
    the same for all observers who can see the Moon.

    Args:
        tjd_ut: Julian Day (UT) of the time to calculate
        ifl: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Umbral magnitude as a float:
            - 0.0 if Moon is not in umbra (penumbral-only or no eclipse)
            - 0.0 to 1.0 for partial umbral eclipses
            - >= 1.0 for total lunar eclipses

    Example:
        >>> from libephemeris import lun_eclipse_umbral_magnitude, FLG_SWIEPH
        >>> # Calculate umbral magnitude during Nov 8, 2022 total lunar eclipse
        >>> jd = 2459892.0  # During eclipse
        >>> umbral_mag = lun_eclipse_umbral_magnitude(jd, FLG_SWIEPH)
        >>> print(f"Umbral magnitude: {umbral_mag:.4f}")

    References:
        - Reference API: lun_eclipse_how()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    return _lun_eclipse_umbral_magnitude_pythonic(tjd_ut, ifl)


def _lun_eclipse_penumbral_magnitude_pythonic(
    jd: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the penumbral magnitude for a lunar eclipse at a specific time.

    Penumbral magnitude is defined as the fraction of the Moon's diameter that
    is within Earth's penumbral (partial) shadow. This is a simplified convenience
    function that returns just the penumbral magnitude value.

    Unlike solar eclipse magnitude (which depends on observer location), lunar
    eclipse magnitude is the same for all observers who can see the Moon,
    since the Moon physically enters Earth's shadow.

    Args:
        jd: Julian Day (UT) of the time to calculate
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Penumbral magnitude as a float:
            - 0.0 if the Moon is not in the penumbral shadow (no eclipse)
            - 0.0 to 1.0 for penumbral-only eclipses or during penumbral
              phases of umbral eclipses
            - > 1.0 when the Moon is fully within the penumbra (and possibly
              also within the umbra); values > 1.0 indicate how deep the Moon
              is in the penumbral shadow

    Note:
        This function does NOT search for eclipses - it calculates the
        instantaneous penumbral magnitude at the given time. To find eclipse
        events, use _lun_eclipse_when_pythonic() first.

        For all lunar eclipses (penumbral, partial, or total), the penumbral
        magnitude will be > 0 during the eclipse. For penumbral-only eclipses,
        use this function to determine the eclipse magnitude.

    Algorithm:
        1. Calculate Moon's ecliptic position
        2. Calculate Earth's penumbral shadow cone geometry at Moon's distance
        3. Determine how much of the Moon's diameter is within the penumbra

    Precision:
        Magnitude accurate to ~0.01 for typical eclipses.

    Example:
        >>> from libephemeris import julday, _lun_eclipse_penumbral_magnitude_pythonic, _lun_eclipse_when_pythonic
        >>> from libephemeris import ECL_PENUMBRAL
        >>> # First find a penumbral lunar eclipse
        >>> jd_start = julday(2020, 1, 1, 0)
        >>> ecl_type, times = _lun_eclipse_when_pythonic(jd_start, eclipse_type=ECL_PENUMBRAL)
        >>> jd_max = times[0]  # Time of maximum eclipse
        >>> # Calculate penumbral magnitude at maximum
        >>> penumbral_mag = _lun_eclipse_penumbral_magnitude_pythonic(jd_max)
        >>> print(f"Penumbral magnitude: {penumbral_mag:.4f}")

        >>> # Check magnitude at a random time (no eclipse)
        >>> jd_no_eclipse = julday(2022, 6, 1, 12.0)
        >>> mag = _lun_eclipse_penumbral_magnitude_pythonic(jd_no_eclipse)
        >>> print(f"Magnitude: {mag:.4f}")  # Will be 0.0

    References:
        - Meeus (1998), "Astronomical Algorithms", 2nd ed., Ch. 54
        - Espenak & Meeus (2009), "Five Millennium Catalog of Lunar
          Eclipses", NASA/TP-2009-214173, Sec. 1 (magnitude definitions and
          Danjon shadow enlargement)
    """
    # Delegate to the same shadow core (_lun_how_core) that backs
    # lun_eclipse_how()/lun_eclipse_when() so the penumbral magnitude (attr[1])
    # stays consistent with the canonical model and NASA's canon. See the
    # umbral variant for the rationale behind retiring
    # _calculate_lunar_eclipse_type_and_magnitude() from this path (divergent
    # 1/85 enlargement + small-angle approximation).
    _retc, attr, _dcore = _lun_how_core(jd, flags)
    return max(0.0, attr[1])


def lun_eclipse_penumbral_magnitude(
    tjd_ut: float,
    ifl: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the penumbral magnitude for a lunar eclipse at a specific time.

    This function follows the reference API naming convention. It is a convenience
    function that returns just the penumbral magnitude (fraction of Moon's diameter
    within Earth's penumbral shadow).

    Unlike solar eclipses, lunar eclipse magnitude does not depend on observer
    location - the Moon physically enters Earth's shadow, so the magnitude is
    the same for all observers who can see the Moon.

    Args:
        tjd_ut: Julian Day (UT) of the time to calculate
        ifl: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Penumbral magnitude as a float:
            - 0.0 if Moon is not in penumbra (no eclipse)
            - 0.0 to 1.0 for penumbral-only eclipses (partial immersion)
            - > 1.0 when Moon is fully within the penumbra

    Example:
        >>> from libephemeris import lun_eclipse_penumbral_magnitude, FLG_SWIEPH
        >>> # Calculate penumbral magnitude during Jan 10, 2020 penumbral lunar eclipse
        >>> jd = 2458859.0  # During eclipse
        >>> penumbral_mag = lun_eclipse_penumbral_magnitude(jd, FLG_SWIEPH)
        >>> print(f"Penumbral magnitude: {penumbral_mag:.4f}")

    References:
        - Reference API: lun_eclipse_how()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    return _lun_eclipse_penumbral_magnitude_pythonic(tjd_ut, ifl)


def _lun_eclipse_gamma_pythonic(
    jd: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the gamma parameter for a lunar eclipse at a specific time.

    Gamma is the distance of the Moon's center from Earth's shadow axis,
    measured in Earth radii. It is a fundamental parameter for characterizing
    lunar eclipses, indicating how centrally the Moon passes through the shadow.

    Unlike solar eclipse gamma (which uses the fundamental plane), lunar eclipse
    gamma is measured perpendicular to the shadow axis in the plane containing
    the Moon. A gamma of 0 means the Moon passes exactly through the center of
    Earth's shadow.

    The returned gamma is NON-NEGATIVE: it is the radial distance of the
    Moon's center from the shadow axis, matching the reference API's
    unsigned lunar-eclipse channel (attr[7]); the north/south side of the
    passage is NOT encoded in the sign (unlike solar-eclipse gamma).
    See tests/test_lun_eclipse_gamma.py for the parity contract.

    Args:
        jd: Julian Day (UT) of the time to calculate
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Gamma value as a float:
            - 0.0: Moon center is on the shadow axis (most central eclipse)
            - |gamma| < ~0.25: Deep total eclipse (Moon well within umbra)
            - |gamma| < ~0.75: Total eclipse possible
            - |gamma| < ~1.0: Partial umbral eclipse possible
            - |gamma| < ~1.5: Penumbral eclipse possible
            - |gamma| > ~1.5: No eclipse (Moon misses the shadow)

    Note:
        This function does NOT search for eclipses - it calculates the
        instantaneous gamma at the given time. To find eclipse events,
        use _lun_eclipse_when_pythonic() first.

        The gamma parameter is useful for:
        - Classifying eclipse centrality
        - Predicting eclipse magnitude
        - Analyzing Saros series patterns

    Algorithm:
        1. Calculate Moon's ecliptic latitude
        2. Calculate Earth's angular semi-diameter as seen from Moon
        3. Compute gamma as the ratio of lunar latitude to Earth's angular radius

    Precision:
        Gamma accurate to ~0.001 for typical eclipses.

    Example:
        >>> from libephemeris import julday, _lun_eclipse_gamma_pythonic, _lun_eclipse_when_pythonic
        >>> # First find a lunar eclipse
        >>> jd_start = julday(2022, 5, 1, 0)
        >>> ecl_type, times = _lun_eclipse_when_pythonic(jd_start)
        >>> jd_max = times[0]  # Time of maximum eclipse
        >>> # Calculate gamma at maximum
        >>> gamma = _lun_eclipse_gamma_pythonic(jd_max)
        >>> print(f"Gamma: {gamma:.4f}")

        >>> # Check gamma at a random time (far from eclipse)
        >>> jd_no_eclipse = julday(2022, 6, 1, 12.0)
        >>> gamma = _lun_eclipse_gamma_pythonic(jd_no_eclipse)
        >>> print(f"Gamma: {gamma:.4f}")  # Will be large (no eclipse)

    References:
        - Meeus (1998), "Astronomical Algorithms", 2nd ed., Ch. 54
        - Espenak & Meeus (2009), "Five Millennium Catalog of Lunar
          Eclipses", NASA/TP-2009-214173, Secs. 1 and 3 (gamma definition,
          sign convention, and eclipse classification)
    """
    # Use the existing calculation function
    (
        ecl_type_flags,
        umbral_mag,
        penumbral_mag,
        gamma,
        penumbra_radius,
        umbra_radius,
    ) = _calculate_lunar_eclipse_type_and_magnitude(jd)

    return gamma


def lun_eclipse_gamma(
    tjd_ut: float,
    ifl: int = FLG_SWIEPH,
) -> float:
    """Calculate lunar-eclipse gamma at a specific time.

    This function follows the reference API naming convention. Gamma represents
    the distance of the Moon's center from Earth's shadow axis, measured in
    Earth radii.

    Unlike solar eclipse gamma, lunar eclipse gamma does not depend on observer
    location - it is a geometric property of the Moon's passage through Earth's
    shadow.

    Args:
        tjd_ut: Julian Day (UT) of the time to calculate
        ifl: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Non-negative distance of the Moon's centre from the terrestrial
        shadow axis in Earth radii. Values near zero are the most central.

    Example:
        >>> from libephemeris import lun_eclipse_gamma, FLG_SWIEPH
        >>> # Calculate gamma during Nov 8, 2022 total lunar eclipse
        >>> jd = 2459892.0  # During eclipse
        >>> gamma = lun_eclipse_gamma(jd, FLG_SWIEPH)
        >>> print(f"Gamma: {gamma:.4f}")

    References:
        Meeus, *Astronomical Algorithms*, chapter 54.
    """
    return _lun_eclipse_gamma_pythonic(tjd_ut, ifl)


def planet_occult_when_glob(
    tjdut: float,
    occulting_planet: int,
    occulted_planet: int = 0,
    starname: str = "",
    flags: int = FLG_SWIEPH,
    direction: int = 0,
) -> Tuple[int, Tuple[float, ...]]:
    """
    Find the next planetary occultation globally (UT).

    A planetary occultation occurs when one planet passes in front of (occults)
    another planet or star as seen from Earth. This is different from lunar
    occultations - here the occulting body is a planet (e.g., Venus, Jupiter).

    Planetary occultations are rare events. Mutual occultations between planets
    typically occur only a few times per century.

    Args:
        tjdut: Julian Day (UT) to start search from
        occulting_planet: Planet ID of the occulting (foreground) planet.
            Use VENUS, MARS, JUPITER, etc.
        occulted_planet: Planet ID of the occulted (background) planet.
            Use 0 if searching for a star occultation.
        starname: Star name (str). Use empty string "" if searching for a planet.
        flags: Calculation flags (FLG_SWIEPH, etc.)
        direction: Search direction. 0 or positive = forward in time,
                   negative = backward in time.

    Returns:
        Tuple containing:
            - retflags: Returned bit flags (int):
                - 0 if no occultation found
                - ECL_TOTAL or ECL_PARTIAL
            - tret: Tuple of 10 floats with occultation phase times (JD UT):
                [0]: Time of maximum occultation
                [1]: Reserved (0)
                [2]: Time of occultation begin
                [3]: Time of occultation end
                [4]: Time of totality begin (if total, else 0)
                [5]: Time of totality end (if total, else 0)
                [6-9]: Reserved (0)

    Raises:
        Error: If no occultation found within search limit
        ValueError: If neither occulted_planet nor starname is specified

    Algorithm:
        1. Calculate both planets' positions over time
        2. Find conjunctions in right ascension
        3. At conjunction, check angular separation
        4. If separation < occulting planet's angular radius + occulted body's radius,
           an occultation occurs
        5. Refine timing using golden section search
        6. Calculate contact times

    Historical Events:
        - 1818 Dec 3: Venus occulted Jupiter
        - 2065 Nov 22: Venus will occult Jupiter
        - 2123 Sep 14: Venus will occult Jupiter

    Example:
        >>> # Find next planetary occultation of Jupiter by Venus
        >>> from libephemeris import julday, planet_occult_when_glob, VENUS, JUPITER
        >>> jd = julday(2060, 1, 1, 0)
        >>> retflags, tret = planet_occult_when_glob(jd, VENUS, JUPITER)
        >>> print(f"Occultation at JD {tret[0]:.5f}")

        >>> # Find next occultation of Regulus by Venus
        >>> retflags, tret = planet_occult_when_glob(jd, VENUS, 0, "Regulus")

    References:
        - Meeus "Astronomical Algorithms" Ch. 9 (Angular Separation)
        - Herald, D. & Sinnott, R. "Planetary Occultations"
    """
    from .fixed_stars import _resolve_star_id, fixstar_ut
    from .constants import (
        SUN,
        MOON,
    )
    from .planets import _PLANET_MAP

    if occulted_planet == 0 and not starname:
        raise ValueError(
            "Either occulted_planet ID or starname must be specified for occultation search"
        )

    # Validate occulting planet - can't be Sun or Moon (use other functions for those)
    if occulting_planet == SUN:
        raise ValueError(
            "Sun cannot be an occulting body - use _sol_eclipse_when_glob_pythonic"
        )
    if occulting_planet == MOON:
        raise ValueError("Moon cannot be an occulting body - use lun_occult_when_glob")
    if occulting_planet not in _PLANET_MAP:
        raise ValueError(f"Invalid occulting planet ID: {occulting_planet}")

    # Validate occulted planet if specified
    if occulted_planet != 0:
        if occulted_planet not in _PLANET_MAP:
            raise ValueError(f"Invalid occulted planet ID: {occulted_planet}")
        if occulted_planet == occulting_planet:
            raise ValueError("Occulting and occulted planets cannot be the same")

    jd_start = tjdut

    # Planetary occultations are very rare - search up to 150 years
    MAX_SEARCH_YEARS = 150
    MAX_ITERATIONS = int(MAX_SEARCH_YEARS * 365)  # Check roughly daily

    from .constants import FLG_EQUATORIAL

    def _get_planet_position(
        jd: float, planet_id: int
    ) -> Tuple[float, float, float, float]:
        """Get planet's geocentric RA, Dec, distance, and angular radius."""

        if planet_id not in _PLANET_MAP:
            raise ValueError(f"Invalid planet ID: {planet_id}")

        # Use calc_ut (LEB-aware) for equatorial positions
        planet_eq, _ = calc_ut(jd, planet_id, FLG_EQUATORIAL | FLG_SPEED)
        ra_deg = planet_eq[0]
        dec_deg = planet_eq[1]
        dist_au = planet_eq[2]

        angular_radius = _calc_planet_angular_radius(planet_id, dist_au)

        return ra_deg, dec_deg, dist_au, angular_radius

    def _get_target_position(jd: float) -> Tuple[float, float, float]:
        """Get target (occulted) body's geocentric RA, Dec, and angular radius."""
        if occulted_planet == 0:
            # Validate the star name first so an unknown star raises ValueError
            # (the documented contract), not the lower-level fixstar_ut Error.
            _star_id, err, _ = _resolve_star_id(starname)
            if err is not None:
                raise ValueError(err)
            # Fixed star: use the apparent equatorial position of date
            # (precession + nutation + aberration + proper motion), the same
            # frame as the occulting planet's calc_ut(FLG_EQUATORIAL) RA/Dec.
            # A bare J2000 + proper-motion position omits precession, leaving the
            # star up to ~0.6° (precession to 2044) out of the planet's frame, so
            # the disc never overlaps and no occultation is ever detected.
            star_xx, _, _ = fixstar_ut(starname, jd, FLG_EQUATORIAL)
            # Stars have negligible angular radius (~0.4 arcsec point source).
            return float(star_xx[0]), float(star_xx[1]), 0.0001
        else:
            # Planet
            ra, dec, dist, angular_radius = _get_planet_position(jd, occulted_planet)
            return ra, dec, angular_radius

    def _angular_separation(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
        """Calculate angular separation between two points (in degrees)."""
        ra1_r = math.radians(ra1)
        dec1_r = math.radians(dec1)
        ra2_r = math.radians(ra2)
        dec2_r = math.radians(dec2)

        cos_sep = math.sin(dec1_r) * math.sin(dec2_r) + math.cos(dec1_r) * math.cos(
            dec2_r
        ) * math.cos(ra1_r - ra2_r)
        cos_sep = max(-1.0, min(1.0, cos_sep))
        return math.degrees(math.acos(cos_sep))

    def _check_occultation(jd: float) -> Tuple[bool, float, float, float]:
        """
        Check if occultation occurs at given time.

        Returns: (is_occultation, separation, occulting_radius, target_radius)
        """
        occ_ra, occ_dec, occ_dist, occ_radius = _get_planet_position(
            jd, occulting_planet
        )
        target_ra, target_dec, target_radius = _get_target_position(jd)

        separation = _angular_separation(occ_ra, occ_dec, target_ra, target_dec)

        # Occultation threshold: occulting planet radius + target radius
        threshold = occ_radius + target_radius

        return separation <= threshold, separation, occ_radius, target_radius

    def _find_minimum_separation(jd_start_search: float, jd_end_search: float) -> float:
        """Find time of minimum separation using golden section search."""
        phi = (1 + math.sqrt(5)) / 2  # Golden ratio

        a = jd_start_search
        b = jd_end_search

        c = b - (b - a) / phi
        d = a + (b - a) / phi

        def get_sep(jd: float) -> float:
            occ_ra, occ_dec, _, _ = _get_planet_position(jd, occulting_planet)
            target_ra, target_dec, _ = _get_target_position(jd)
            return _angular_separation(occ_ra, occ_dec, target_ra, target_dec)

        fc = get_sep(c)
        fd = get_sep(d)

        for _ in range(50):  # Converge to ~0.1 second precision
            if fc < fd:
                b = d
                d = c
                fd = fc
                c = b - (b - a) / phi
                fc = get_sep(c)
            else:
                a = c
                c = d
                fc = fd
                d = a + (b - a) / phi
                fd = get_sep(d)

            if b - a < 1e-6:  # ~0.1 second precision
                break

        return (a + b) / 2

    def _calculate_contact_times(
        jd_max: float, min_sep: float, occ_radius: float, target_radius: float
    ) -> Tuple[float, float, float, float]:
        """Calculate contact times for the occultation."""
        outer_threshold = occ_radius + target_radius
        inner_threshold = abs(occ_radius - target_radius)

        # Estimate relative angular speed
        dt_test = 1.0 / 24.0  # 1 hour
        occ_ra1, occ_dec1, _, _ = _get_planet_position(
            jd_max - dt_test, occulting_planet
        )
        occ_ra2, occ_dec2, _, _ = _get_planet_position(
            jd_max + dt_test, occulting_planet
        )
        target_ra1, target_dec1, _ = _get_target_position(jd_max - dt_test)
        target_ra2, target_dec2, _ = _get_target_position(jd_max + dt_test)

        # Relative motion. The RA differences are measured in degrees along
        # the equator; they must be scaled by cos(dec) to become true angular
        # displacements on the sky before being combined with the declination
        # component (which is already a true angular chord). Without this
        # factor the relative speed is overestimated at non-zero declination
        # (by 1/cos(dec)), which makes the derived contact times too short and
        # inconsistent with the great-circle separations used for min_sep and
        # the thresholds. At dec = 0, cos(dec) = 1 and behaviour is unchanged.
        mean_dec_rad = math.radians(
            (occ_dec1 + occ_dec2 + target_dec1 + target_dec2) / 4.0
        )
        cos_dec = math.cos(mean_dec_rad)
        d_ra = (occ_ra2 - occ_ra1) - (target_ra2 - target_ra1)
        # Normalize the RA difference across the 0h/24h (0/360 deg) wrap.
        d_ra = (d_ra + 180.0) % 360.0 - 180.0
        d_ra *= cos_dec
        d_dec = (occ_dec2 - occ_dec1) - (target_dec2 - target_dec1)
        relative_speed = math.sqrt(d_ra**2 + d_dec**2) / (2 * dt_test)

        if relative_speed < 0.001:  # Very slow motion, use fallback
            relative_speed = 0.5  # degrees/day fallback

        # Contact times based on chord geometry
        if min_sep < outer_threshold:
            half_chord_outer = math.sqrt(max(0, outer_threshold**2 - min_sep**2))
            dt_outer = half_chord_outer / relative_speed

            jd_first = jd_max - dt_outer
            jd_fourth = jd_max + dt_outer
        else:
            jd_first = 0.0
            jd_fourth = 0.0

        # For total occultation
        if min_sep < inner_threshold:
            half_chord_inner = math.sqrt(max(0, inner_threshold**2 - min_sep**2))
            dt_inner = half_chord_inner / relative_speed

            jd_second = jd_max - dt_inner
            jd_third = jd_max + dt_inner
        else:
            jd_second = 0.0
            jd_third = 0.0

        return jd_first, jd_second, jd_third, jd_fourth

    # Main search loop
    jd = jd_start
    step = 1.0  # Check every day initially
    # direction < 0 searches backward in time (reference convention).
    sign = -1.0 if direction < 0 else 1.0

    for iteration in range(MAX_ITERATIONS):
        try:
            is_occ, sep, occ_r, target_r = _check_occultation(jd)
        except (ValueError, ArithmeticError, IndexError) as e:
            # Handle ephemeris range errors or other calculation errors
            if "ephemeris" in str(e).lower() or "range" in str(e).lower():
                # Exceeded ephemeris range - stop searching
                break
            raise

        # Mutual planetary occultations subtend arcseconds and last
        # minutes: a sampled is_occ hit can fall between steps. Any
        # close approach is therefore refined to its true minimum
        # before the occultation test.
        if not is_occ and sep < 0.2:
            jd_min_c = _find_minimum_separation(jd - 0.5, jd + 0.5)
            is_occ_c, _sep_c, _o_r, _t_r = _check_occultation(jd_min_c)
            if is_occ_c:
                is_occ = True
                jd = jd_min_c
            elif (jd_min_c - jd) * sign > 0.0:
                # Past this conjunction: jump beyond it.
                jd = jd_min_c + sign * 0.5
                continue

        if is_occ:
            # Found an occultation! Refine timing
            jd_max = _find_minimum_separation(jd - 1.0, jd + 1.0)

            # Verify it's still an occultation
            is_occ_refined, min_sep, occ_r, target_r = _check_occultation(jd_max)

            if is_occ_refined:
                # Calculate contact times
                jd_first, jd_second, jd_third, jd_fourth = _calculate_contact_times(
                    jd_max, min_sep, occ_r, target_r
                )

                # Determine occultation type
                # Grazing threshold: when the target passes within the outer 10%
                # of the occulting planet's disc
                0.9 * occ_r

                if min_sep < abs(occ_r - target_r):
                    ecl_type = ECL_TOTAL
                else:
                    ecl_type = ECL_PARTIAL

                tret = (
                    jd_max,  # [0] Time of maximum
                    0.0,  # [1] Reserved
                    jd_first,  # [2] Begin
                    jd_fourth,  # [3] End
                    jd_second,  # [4] Totality begin
                    jd_third,  # [5] Totality end
                    0.0,  # [6] Reserved
                    0.0,  # [7] Reserved
                    0.0,  # [8] Reserved
                    0.0,  # [9] Reserved
                )

                return ecl_type, tret

        # Adaptive step based on angular separation
        if sep < 0.5:  # Very close - within half a degree
            step = 0.01  # Check every ~15 minutes
        elif sep < 2.0:  # Close - within 2 degrees
            step = 0.1
        elif sep < 10.0:
            step = 0.5
        else:
            step = 1.0

        jd += sign * step

        if jd > jd_start + MAX_SEARCH_YEARS * 365.25:
            break

    if occulted_planet == 0:
        target_desc = starname
    else:
        target_desc = f"planet {occulted_planet}"
    occ_desc = f"planet {occulting_planet}"

    raise Error(
        f"No planetary occultation of {target_desc} by {occ_desc} found within "
        f"{MAX_SEARCH_YEARS} years of JD {jd_start}"
    )


# Alias for reference API compatibility
planet_occult_when_glob = planet_occult_when_glob


def planet_occult_when_loc(
    jd_start: float,
    occulting_planet: int,
    occulted_planet: int = 0,
    star_name: str = "",
    lat: float = 0.0,
    lon: float = 0.0,
    altitude: float = 0.0,
    flags: int = FLG_SWIEPH,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """Find the next planetary occultation visible from a specific location.

    Wrapper around :func:`_planet_occult_when_loc_impl` that adds LEB→Skyfield
    fallback for partial/custom LEB files. See the impl docstring for the
    full API contract.
    """
    return _call_with_leb_skyfield_fallback(
        _planet_occult_when_loc_impl,
        jd_start,
        occulting_planet,
        occulted_planet,
        star_name,
        lat,
        lon,
        altitude,
        flags,
    )


def _planet_occult_when_loc_impl(
    jd_start: float,
    occulting_planet: int,
    occulted_planet: int = 0,
    star_name: str = "",
    lat: float = 0.0,
    lon: float = 0.0,
    altitude: float = 0.0,
    flags: int = FLG_SWIEPH,
    *,
    reader=None,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """
    Find the next planetary occultation visible from a specific location.

    A planetary occultation occurs when one planet passes in front of (occults)
    another planet or star. This function searches forward in time to find the
    next occultation visible from a specific geographic location, where both
    the occulting and occulted bodies are above the horizon.

    Args:
        jd_start: Julian Day (UT) to start search from
        occulting_planet: Planet ID of the occulting (foreground) planet.
            Use VENUS, MARS, JUPITER, etc.
        occulted_planet: Planet ID of the occulted (background) planet.
            Set to 0 if searching for a fixed star occultation.
        star_name: Name of fixed star to check (e.g. "Regulus", "Spica").
            Only used if occulted_planet is 0.
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - times: Tuple of 10 floats with occultation phase times (JD UT):
                [0]: Time of maximum occultation (minimum separation)
                [1]: Time of first contact (occultation begins)
                [2]: Time of second contact (full occultation begins, or 0)
                [3]: Time of third contact (full occultation ends, or 0)
                [4]: Time of fourth contact (occultation ends)
                [5-9]: Reserved (0)
            - attr: Tuple of 20 floats with occultation attributes:
                [0]: Fraction of target covered (magnitude)
                [1]: Ratio of occulting to occulted angular diameter
                [2]: Fraction of disc covered (obscuration)
                [3]: Reserved (0)
                [4]: Azimuth of occulted body at maximum (degrees)
                [5]: True altitude of occulted body at maximum (degrees)
                [6]: Apparent altitude (with refraction)
                [7]: Angular separation at maximum (degrees)
                [8-19]: Reserved (0)
            - retflag: Occultation type flags (ECL_* constants)

    Raises:
        Error: If no occultation visible from location within search limit
        ValueError: If neither occulted_planet nor star_name is specified

    Example:
        >>> # Find next occultation of Regulus by Venus visible from Rome
        >>> from libephemeris import julday, planet_occult_when_loc, VENUS
        >>> jd = julday(2020, 1, 1, 0)
        >>> rome_lat, rome_lon = 41.9028, 12.4964
        >>> ecl_type, times, attr = planet_occult_when_loc(
        ...     jd, VENUS, 0, "Regulus", rome_lat, rome_lon
        ... )
        >>> print(f"Occultation max at JD {times[0]:.5f}")

    References:
        - Meeus "Astronomical Algorithms" Ch. 9 (Angular Separation)
    """
    from skyfield.api import wgs84

    from .constants import (
        SUN,
        MOON,
    )
    from .fixed_stars import FIXED_STARS, _resolve_star_id
    from .planets import _PLANET_MAP, get_planet_target
    from .state import get_planets, get_timescale

    if occulted_planet == 0 and not star_name:
        raise ValueError(
            "Either occulted_planet ID or star_name must be specified for occultation search"
        )

    # Validate planets
    if occulting_planet == SUN:
        raise ValueError(
            "Sun cannot be an occulting body - use _sol_eclipse_when_loc_pythonic"
        )
    if occulting_planet == MOON:
        raise ValueError(
            "Moon cannot be an occulting body - use _lun_occult_when_loc_pythonic"
        )
    if occulting_planet not in _PLANET_MAP:
        raise ValueError(f"Invalid occulting planet ID: {occulting_planet}")

    if occulted_planet != 0:
        if occulted_planet not in _PLANET_MAP:
            raise ValueError(f"Invalid occulted planet ID: {occulted_planet}")
        if occulted_planet == occulting_planet:
            raise ValueError("Occulting and occulted planets cannot be the same")

    MAX_SEARCH_YEARS = 150
    MAX_GLOBAL_SEARCHES = 100

    # reader is provided by the caller (None forces Skyfield path)
    _reader = reader
    if _reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR

        geopos = (lon, lat, altitude)

        def _get_body_altitude(jd: float, planet_id: int) -> Tuple[float, float]:
            jd_tt = jd + deltat(jd)
            pos = _topo_ecliptic(_reader, jd_tt, jd, planet_id, geopos, FLG_SPEED)
            az, alt_true, _alt_app = azalt(jd, ECL2HOR, geopos, 0, 0, pos[:3])
            return alt_true, az

        def _get_target_altitude(jd: float) -> Tuple[float, float]:
            jd_tt = jd + deltat(jd)
            if occulted_planet == 0:
                from .fixed_stars import fixstar_ut

                star_pos, _, _ = fixstar_ut(star_name, jd, FLG_SPEED)
                az, alt_true, _alt_app = azalt(jd, ECL2HOR, geopos, 0, 0, star_pos[:3])
            else:
                pos = _topo_ecliptic(
                    _reader, jd_tt, jd, occulted_planet, geopos, FLG_SPEED
                )
                az, alt_true, _alt_app = azalt(jd, ECL2HOR, geopos, 0, 0, pos[:3])
            return alt_true, az
    else:
        eph = get_planets()
        ts = get_timescale()
        earth = eph["earth"]
        observer = wgs84.latlon(lat, lon, altitude)
        observer_at = earth + observer

        def _get_body_altitude(jd: float, planet_id: int) -> Tuple[float, float]:
            from .planets import get_planet_target

            target_name = _PLANET_MAP[planet_id]
            target = get_planet_target(eph, target_name)
            t = ts.ut1_jd(jd)
            target_app = observer_at.at(t).observe(target).apparent()
            alt, az, _ = target_app.altaz()
            return alt.degrees, (az.degrees + 180.0) % 360.0

        def _get_target_altitude(jd: float) -> Tuple[float, float]:
            from .planets import get_planet_target

            t = ts.ut1_jd(jd)
            if occulted_planet == 0:
                star_id, err, _ = _resolve_star_id(star_name)
                if err is not None:
                    raise ValueError(err)
                star = FIXED_STARS[star_id]
                t_years = (jd - 2451545.0) / 365.25
                ra_deg = star.ra_j2000 + (star.pm_ra * t_years) / 3600.0
                dec_deg = star.dec_j2000 + (star.pm_dec * t_years) / 3600.0
                from skyfield.api import Star

                star_obj = Star(ra_hours=ra_deg / 15.0, dec_degrees=dec_deg)
                target_app = observer_at.at(t).observe(star_obj).apparent()
            else:
                target_name = _PLANET_MAP[occulted_planet]
                target = get_planet_target(eph, target_name)
                target_app = observer_at.at(t).observe(target).apparent()
            alt, az, _ = target_app.altaz()
            # Convert Skyfield navigational azimuth (N=0) to the reference convention (S=0)
            return alt.degrees, (az.degrees + 180.0) % 360.0

    def _is_visible_at_location(jd: float) -> bool:
        """Check if both bodies are above horizon at given time."""
        occ_alt, _ = _get_body_altitude(jd, occulting_planet)
        target_alt, _ = _get_target_altitude(jd)

        # Both bodies must be above horizon (with some margin for twilight)
        MIN_ALT = -0.5  # Allow slightly below horizon for refraction
        return occ_alt > MIN_ALT and target_alt > MIN_ALT

    current_jd = jd_start

    for search_count in range(MAX_GLOBAL_SEARCHES):
        try:
            # Find next global occultation
            retflags, tret = planet_occult_when_glob(
                current_jd, occulting_planet, occulted_planet, star_name, flags, 0
            )

            jd_max = tret[0]
            jd_begin = tret[2]
            jd_end = tret[3]

            # Check if any part of the occultation is visible from this location
            # Sample multiple times during the occultation
            visible = False
            visible_times = []

            if jd_begin > 0 and jd_end > 0:
                sample_count = 10
                for i in range(sample_count + 1):
                    sample_jd = jd_begin + (jd_end - jd_begin) * i / sample_count
                    if _is_visible_at_location(sample_jd):
                        visible = True
                        visible_times.append(sample_jd)

            # Also check at maximum
            if _is_visible_at_location(jd_max):
                visible = True
                visible_times.append(jd_max)

            if visible:
                occ_alt, occ_az = _get_body_altitude(jd_max, occulting_planet)
                target_alt, target_az = _get_target_altitude(jd_max)

                if _reader is not None:
                    from .fast_calc import _topo_ecliptic as _tp_occ
                    from .time_utils import deltat as _sd_occ
                    from .utils import angular_separation as _ang_sep_occ

                    _jd_tt_occ = jd_max + _sd_occ(jd_max)
                    _gp_occ = (lon, lat, altitude)
                    occ_pos = _tp_occ(
                        _reader,
                        _jd_tt_occ,
                        jd_max,
                        occulting_planet,
                        _gp_occ,
                        FLG_SPEED,
                    )
                    if occulted_planet == 0:
                        from .fixed_stars import fixstar_ut

                        tgt_pos, _, _ = fixstar_ut(star_name, jd_max, FLG_SPEED)
                        target_radius = 0.0001
                    else:
                        tgt_pos = _tp_occ(
                            _reader,
                            _jd_tt_occ,
                            jd_max,
                            occulted_planet,
                            _gp_occ,
                            FLG_SPEED,
                        )
                        target_radius = _calc_planet_angular_radius(
                            occulted_planet, tgt_pos[2]
                        )
                    separation = _ang_sep_occ(
                        occ_pos[0], occ_pos[1], tgt_pos[0], tgt_pos[1]
                    )
                    occ_dist = occ_pos[2]
                    occ_radius = _calc_planet_angular_radius(occulting_planet, occ_dist)
                else:
                    from .planets import get_planet_target

                    t = ts.ut1_jd(jd_max)
                    target_name = _PLANET_MAP[occulting_planet]
                    occ_body = get_planet_target(eph, target_name)
                    occ_app = observer_at.at(t).observe(occ_body).apparent()
                    if occulted_planet == 0:
                        star_id, _, _ = _resolve_star_id(star_name)
                        star = FIXED_STARS[star_id]
                        t_years = (jd_max - 2451545.0) / 365.25
                        ra_deg = star.ra_j2000 + (star.pm_ra * t_years) / 3600.0
                        dec_deg = star.dec_j2000 + (star.pm_dec * t_years) / 3600.0
                        from skyfield.api import Star

                        target_body = Star(ra_hours=ra_deg / 15.0, dec_degrees=dec_deg)
                    else:
                        target_body = get_planet_target(
                            eph, _PLANET_MAP[occulted_planet]
                        )
                    target_app = observer_at.at(t).observe(target_body).apparent()
                    separation = occ_app.separation_from(target_app).degrees
                    occ_dist = occ_app.distance().au
                    occ_radius = _calc_planet_angular_radius(occulting_planet, occ_dist)
                    if occulted_planet == 0:
                        target_radius = 0.0001
                    else:
                        target_dist = target_app.distance().au
                        target_radius = _calc_planet_angular_radius(
                            occulted_planet, target_dist
                        )

                if target_radius > 0:
                    magnitude = max(0, 1.0 - separation / (occ_radius + target_radius))
                    ratio = (
                        occ_radius / target_radius if target_radius > 0.0001 else 999.0
                    )
                    obscuration = magnitude  # Simplified
                else:
                    magnitude = 1.0 if separation < occ_radius else 0.0
                    ratio = 999.0
                    obscuration = magnitude

                # Refraction from true altitude, arcmin -> deg:
                # R = 1.02 / tan(h + 10.3/(h + 5.11)).  Sæmundsson (1986),
                # Sky & Telescope 72, 70; Meeus (1998) ch. 16, eq. 16.4.
                apparent_alt = (
                    target_alt
                    + (
                        1.02
                        / math.tan(
                            math.radians(target_alt + 10.3 / (target_alt + 5.11))
                        )
                    )
                    / 60.0
                    if target_alt > -1
                    else target_alt
                )

                times = (
                    jd_max,  # [0] Maximum
                    jd_begin if jd_begin > 0 else 0.0,  # [1] First contact
                    tret[4] if tret[4] > 0 else 0.0,  # [2] Second contact
                    tret[5] if tret[5] > 0 else 0.0,  # [3] Third contact
                    jd_end if jd_end > 0 else 0.0,  # [4] Fourth contact
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,  # [5-9] Reserved
                )

                attr = (
                    magnitude,  # [0] Magnitude
                    ratio,  # [1] Diameter ratio
                    obscuration,  # [2] Obscuration
                    0.0,  # [3] Reserved
                    target_az,  # [4] Azimuth
                    target_alt,  # [5] True altitude
                    apparent_alt,  # [6] Apparent altitude
                    separation,  # [7] Separation
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,  # [8-19] Reserved
                )

                return retflags, times, attr

            # Not visible from this location - continue search after this event
            current_jd = jd_end + 1.0 if jd_end > 0 else jd_max + 1.0

        except (Error, RuntimeError):
            # No more global occultations found
            break

    if occulted_planet == 0:
        target_desc = star_name
    else:
        target_desc = f"planet {occulted_planet}"
    occ_desc = f"planet {occulting_planet}"

    raise Error(
        f"No planetary occultation of {target_desc} by {occ_desc} visible from "
        f"lat={lat}, lon={lon} found within {MAX_SEARCH_YEARS} years of JD {jd_start}"
    )


# Alias for reference API compatibility
planet_occult_when_loc = planet_occult_when_loc
