# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Planetary position calculations for libephemeris.

This is the core module providing planet position calculations
using NASA JPL DE440 ephemeris via Skyfield.

Supported Bodies:
- Classical planets: Sun, Moon, Mercury, Venus, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
- Earth position
- Lunar nodes (Mean/True) via lunar.py
- Lilith/Lunar apogee (Mean/True) via lunar.py
- Minor bodies (asteroids, TNOs) via minor_bodies.py
- Fixed stars via fixed_stars.py
- Astrological angles via angles.py
- Arabic parts via arabic_parts.py

Main Functions:
- calc_ut(): Calculate positions in Universal Time
- calc(): Calculate positions in Ephemeris Time
- set_sid_mode(): Set sidereal zodiac mode
- get_ayanamsa_ut(): Get ayanamsha value

Coordinate Systems:
- Geocentric tropical (default)
- Heliocentric (with FLG_HELCTR)
- Topocentric (requires set_topo)
- Sidereal (requires set_sid_mode)

Precision Notes:
- Nutation: IAU 2006/2000A model via pyerfa (~0.01-0.05 mas precision)
- Obliquity: Vondrak 2011 long-term mean ecliptic/equator geometry plus ERFA
  nutation, shared by all code paths
- Precession: Vondrak, Capitaine & Wallace (2011), with ERFA frame bias
- Ayanamsa: Properly converts ET to UT using Delta T
- Planet positions: JPL DE440; uncertainty is body-, observable-, and
  epoch-dependent as documented by Park et al. (2021), not one universal
  angular bound
- Planets use NAIF planet center IDs (599, 699, etc.) for accurate positions
- Ecliptic frame uses the registered true-ecliptic-of-date reduction

References:
- Park et al. (2021), JPL DE440/DE441, AJ 161:105
- IERS Conventions (2010), IAU 2006/2000A through ERFA
- Vondrak, Capitaine & Wallace (2011), A&A 534:A22

Provenance:
    Geometric states originate in public JPL DE kernels, JPL Horizons/SPK, or a
    separately registered analytical/catalog model. The apparent-place chain
    applies explicit light-time, gravitational deflection, aberration,
    precession/nutation, observer, and output-frame operations from public
    JPL/IAU/IERS conventions. Photometric and physical constants are attributed
    beside their tables. Backend dispatch, finite-difference steps, caches, and
    flag/error behavior are project choices. Public compatibility observations
    constrain API semantics only and never supply an algorithm or coefficient.
"""

from __future__ import annotations

from contextlib import contextmanager as _contextmanager
from contextvars import ContextVar
from functools import lru_cache
from collections.abc import Iterator

import math
import warnings
from typing import Optional, Tuple, TYPE_CHECKING

from .tracing import (
    _is_active,
    _record,
    _record_alias,
    _restore_record,
    _snapshot_record,
)
from jplephem.exceptions import OutOfRangeError
from skyfield.errors import EphemerisRangeError as SkyfieldRangeError
from skyfield.api import Star
from skyfield.framelib import ecliptic_frame
from skyfield.relativity import add_aberration
from dataclasses import dataclass
import erfa

from .precession_vondrak import (
    method_b_accumulated_precession,
    vondrak_mean_obliquity_deg,
    vondrak_mean_obliquity_rad,
    vondrak_pn_matrix,
    vondrak_precession_matrix,
)
from .ayanamsha_definitions import (
    ALDEBARAN_TARGET_LON,
    AYANAMSHA_DEFINING,
    DYNAMIC_AYANAMSHA_MODES,
    GALCENT_TARGET_LON,
    MARDYKS_DEFINING,
    MULA_WILHELM_TARGET_LON,
    SHEORAN_TARGET_LON,
    VALENS_MOON_AYAN_T0_DEG,
    VALENS_MOON_T0_UT,
)

from .exceptions import EphemerisRangeError, LEBCorruptionError, UnknownBodyError

if TYPE_CHECKING:
    pass

# Combined exception tuple for catching both jplephem and Skyfield range errors
_RANGE_ERRORS = (OutOfRangeError, SkyfieldRangeError)

# Central-difference half-steps chosen in whole SI seconds. The Moon uses a
# six-second interval to resolve its faster apparent motion; other bodies use
# one second. These are numerical-analysis choices, not ephemeris fit values.
_MOON_SPEED_HALF_STEP_DAYS = 6.0 / 86400.0
_BODY_SPEED_HALF_STEP_DAYS = 1.0 / 86400.0
# Geometric (FLG_TRUEPOS) place of an SPK type-21 minor body. The type-21
# difference-table evaluation carries ~4e-7" of per-sample jitter, which a
# one-second central difference amplifies into a spurious +0.007..+0.02"/day
# longitude-speed bias (the derivative stops matching the body's own geometric
# positions). The apparent place is unaffected — Skyfield's observe/apparent
# reduction smooths it — so this larger step is used ONLY for the geometric
# type-21 speed. A five-minute half-step drops the jitter contribution below
# 1e-8"/day while these slow bodies' angular acceleration keeps the O(h^2)
# truncation negligible, restoring reference-grade agreement.
_ASTEROID_TRUEPOS_SPEED_HALF_STEP_DAYS = 300.0 / 86400.0

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
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
    OSCU_APOG,
    INTP_APOG,
    INTP_PERG,
    CUPIDO,
    HADES,
    ZEUS,
    KRONOS,
    APOLLON,
    ADMETOS,
    VULKANUS,
    POSEIDON,
    ISIS,
    AST_OFFSET,
    NIBIRU,
    HARRINGTON,
    NEPTUNE_LEVERRIER,
    NEPTUNE_ADAMS,
    PLUTO_LOWELL,
    PLUTO_PICKERING,
    CHIRON,
    PHOLUS,
    CERES,
    PALLAS,
    JUNO,
    VESTA,
    FIXSTAR_OFFSET,
    WHITE_MOON,
    VULCAN,
    PROSERPINA,
    WALDEMATH,
    FLG_SPEED,
    FLG_SWIEPH,
    FLG_HELCTR,
    FLG_TOPOCTR,
    FLG_CENTER_BODY,
    FLG_ORBEL_AA,
    FLG_SIDEREAL,
    ANGLE_OFFSET,
    ARABIC_OFFSET,
    FLG_BARYCTR,
    FLG_TRUEPOS,
    FLG_NOABERR,
    FLG_EQUATORIAL,
    FLG_J2000,
    FLG_NONUT,
    FLG_XYZ,
    FLG_RADIANS,
    FLG_ICRS,
    FLG_SPEED3,
    FLG_NOGDEFL,
    PARS_FORTUNAE,
    PARS_SPIRITUS,
    PARS_AMORIS,
    PARS_FIDEI,
    SIDBIT_ECL_DATE,
    SIDBIT_ECL_T0,
    SIDBIT_NO_PREC_OFFSET,
    SIDBIT_PREC_ORIG,
    SIDBIT_SSY_PLANE,
    SIDBIT_USER_UT,
    SIDM_FAGAN_BRADLEY,
    SIDM_LAHIRI,
    SIDM_DELUCE,
    SIDM_RAMAN,
    SIDM_USHASHASHI,
    SIDM_KRISHNAMURTI,
    SIDM_DJWHAL_KHUL,
    SIDM_YUKTESHWAR,
    SIDM_JN_BHASIN,
    SIDM_BABYL_KUGLER1,
    SIDM_BABYL_KUGLER2,
    SIDM_BABYL_KUGLER3,
    SIDM_BABYL_HUBER,
    SIDM_BABYL_ETPSC,
    SIDM_BABYL_BRITTON,
    SIDM_ALDEBARAN_15TAU,
    SIDM_TRUE_CITRA,
    SIDM_TRUE_REVATI,
    SIDM_TRUE_PUSHYA,
    SIDM_TRUE_MULA,
    SIDM_TRUE_SHEORAN,
    SIDM_HIPPARCHOS,
    SIDM_SASSANIAN,
    SIDM_J2000,
    SIDM_J1900,
    SIDM_B1950,
    SIDM_SURYASIDDHANTA,
    SIDM_SURYASIDDHANTA_MSUN,
    SIDM_ARYABHATA,
    SIDM_ARYABHATA_MSUN,
    SIDM_ARYABHATA_522,
    SIDM_SS_REVATI,
    SIDM_SS_CITRA,
    SIDM_GALCENT_0SAG,
    SIDM_GALCENT_RGILBRAND,
    SIDM_GALCENT_MULA_WILHELM,
    SIDM_GALCENT_COCHRANE,
    SIDM_GALEQU_IAU1958,
    SIDM_GALEQU_TRUE,
    SIDM_GALEQU_MULA,
    SIDM_GALEQU_FIORENZA,
    SIDM_GALALIGN_MARDYKS,
    SIDM_VALENS_MOON,
    SIDM_LAHIRI_1940,
    SIDM_LAHIRI_VP285,
    SIDM_KRISHNAMURTI_VP291,
    SIDM_LAHIRI_ICRC,
    SIDM_USER,
    NODBIT_MEAN,
    NODBIT_OSCU,
    NODBIT_OSCU_BAR,
    NODBIT_FOPOINT,
)

# Persisted core identities are semantic, not a filename convention.  This is
# intentionally body-based so a monolithic ``ephemeris_<tier>.leb`` receives
# the same fail-closed range policy as ``<tier>_core.leb2``.
_LEB_CORE_BODY_IDS = frozenset(range(SUN, MEAN_APOG + 1)) | {EARTH}

# Core bodies whose only sealed source is their own stored LEB channel.  The
# node/apsis points are excluded: the vector pipeline derives them from Moon
# states, and the mean points have runtime analytical models, so a file that
# omits their legacy channels can still serve them without leaving LEB.
_LEB_STORED_STATE_BODY_IDS = frozenset(range(SUN, PLUTO + 1)) | {EARTH}

# Import all sidereal mode constants (SIDM_*)
from .constants import *  # noqa: F403, F401
from .state import (
    _get_computation_ephemeris,
    get_sid_mode,
    get_timescale,
    get_topo,
)

# Planet mapping: Primary names for planets
# For outer planets, uses planet center (NAIF x99) if available in ephemeris,
# otherwise falls back to system barycenter (NAIF x)
# DE421: has centers for Mercury/Venus/Earth/Mars (199/299/399/499), barycenters for Jupiter+
# DE430/440/441: has centers only for Mercury/Venus/Earth (199/299/399), barycenters for Mars+
_PLANET_MAP = {
    SUN: "sun",
    MOON: "moon",
    MERCURY: "mercury",  # 199
    VENUS: "venus",  # 299
    MARS: "mars",  # 499
    JUPITER: "jupiter",  # 599 if available, else barycenter 5
    SATURN: "saturn",  # 699 if available, else barycenter 6
    URANUS: "uranus",  # 799 if available, else barycenter 7
    NEPTUNE: "neptune",  # 899 if available, else barycenter 8
    PLUTO: "pluto",  # 999 if available, else barycenter 9
    EARTH: "earth",
}

_dispatch_source: ContextVar[str | None] = ContextVar(
    "libephemeris_dispatch_source", default=None
)
_source_log_capture: ContextVar[list[str] | None] = ContextVar(
    "libephemeris_source_log_capture", default=None
)
_spk_center_forced_barycenters: ContextVar[frozenset[str]] = ContextVar(
    "libephemeris_spk_center_forced_barycenters", default=frozenset()
)


def _mark_dispatch_source(source: str) -> None:
    """Capture a precise inner source for publication by the outer entry point."""
    _dispatch_source.set(source)


@_contextmanager
def _capture_source_logs() -> Iterator[list[str]]:
    """Defer nested success-source logs until the public caller succeeds."""
    captured: list[str] = []
    token = _source_log_capture.set(captured)
    try:
        yield captured
    finally:
        _source_log_capture.reset(token)


def _fallback_trace_source(ipl: int, source_hint: str | None = None) -> str | None:
    """Return the source tag for a successful ``_calc_body`` fallback.

    Minor bodies and registered planetary moons record their specific source
    inside their calculation branch. Fixed stars expose the source captured by
    their calculation helper. Standard planets and
    lunar points derived from JPL Moon states use the Skyfield pipeline; the
    remaining supported bodies are computed by local analytical formulae
    (some of which use JPL vectors for frame or observer reduction).
    """
    from . import fixed_stars, minor_bodies, planetary_moons

    if source_hint is not None:
        return source_hint

    if (
        ipl in minor_bodies.MINOR_BODY_ELEMENTS
        or AST_OFFSET < ipl < FIXSTAR_OFFSET
        or planetary_moons.is_planetary_moon(ipl)
    ):
        return None
    if ipl in fixed_stars.FIXED_STARS:
        return fixed_stars._get_fixed_star_source()
    if ipl in (TRUE_NODE, OSCU_APOG) or ipl in _PLANET_MAP:
        return "Skyfield"
    return "Analytical"


def _log_successful_source(
    logger, body_id: int, jd: float, source: str, *, context: bool = False
) -> None:
    """Log one source only after all public output processing has succeeded."""
    captured = _source_log_capture.get()
    if captured is not None:
        captured.append(source)
        return
    suffix = " (context)" if context else ""
    if source == "Keplerian":
        from . import minor_bodies

        entries = minor_bodies._merged_epoch_entries(body_id)
        if entries:
            nearest = min(entries, key=lambda item: abs(jd - item.epoch))
            offset_years = abs(jd - nearest.epoch) / 365.25
            detail = (
                "multi-epoch approximation, "
                f"nearest-anchor-offset={offset_years:.2f}y, "
                f"anchor-range=[{entries[0].epoch:.1f}, {entries[-1].epoch:.1f}]"
            )
        else:
            detail = "single-epoch approximation"
        logger.warning(
            "body=%d jd=%.1f source=Keplerian(%s)%s. "
            "Precision is body- and date-dependent; this source is not "
            "classified as ephemeris-grade.",
            body_id,
            jd,
            detail,
            suffix,
        )
    else:
        logger.debug("body=%d jd=%.1f source=%s%s", body_id, jd, source, suffix)


def _raise_leb_range_miss(
    body_id: int, jd: float, context_reader: object | None = None
) -> None:
    """Fail closed for core states the active LEB configuration cannot serve.

    Companion minor bodies can have scientifically meaningful stored windows
    narrower than the core tier.  Their miss is allowed to reach the declared
    local model; core states must never be synthesized beyond their LEB range.
    A core body absent from the active reader (custom/partial file) is a
    declared missing-body condition, not a date-range one: without this the
    sealed vector path fails later with an undeclared ``KeyError``.

    ``context_reader`` is the context-local reader that actually served the
    failing attempt, when it differs from the global one.  The global reader
    keeps deciding whether the miss may fall through (it can still serve the
    state), but when it cannot, the typed error must describe the file that
    failed rather than misreport the global configuration.
    """
    from .inventory import get_body_coverage, get_reader_body_coverage
    from .state import get_calc_mode

    if get_calc_mode() != "leb" or body_id not in _LEB_CORE_BODY_IDS:
        return
    body_coverage = get_body_coverage(body_id, jd)
    if body_coverage is not None and body_coverage.contains(jd):
        return
    if context_reader is not None:
        local_coverage = get_reader_body_coverage(context_reader, body_id, jd)
        if local_coverage is not None and local_coverage.contains(jd):
            return
        body_coverage = local_coverage
    if body_coverage is None:
        if body_id not in _LEB_STORED_STATE_BODY_IDS:
            return
        raise UnknownBodyError(
            message=(
                f"Body {body_id} ({_PLANET_NAMES.get(body_id, 'unnamed')}) is "
                "not stored in the active LEB file. LEB mode does not silently "
                "substitute a non-LEB source for a missing core body."
            ),
            body_id=body_id,
        )
    raise EphemerisRangeError(
        message=(
            f"Body {body_id} at JD {jd:.6f} is outside active LEB coverage range "
            f"[{body_coverage.jd_start:.6f}, {body_coverage.jd_end:.6f}]. "
            "LEB mode does not silently substitute a lower-precision source."
        ),
        requested_jd=jd,
        start_jd=body_coverage.jd_start,
        end_jd=body_coverage.jd_end,
        body_id=body_id,
        ephemeris_file=body_coverage.data_file,
    )


# Fallback mapping when planet center not available in ephemeris
# DE430/440/441 only contain barycenters for Mars and outer planets
_PLANET_FALLBACK = {
    "mars": "mars barycenter",
    "jupiter": "jupiter barycenter",
    "saturn": "saturn barycenter",
    "uranus": "uranus barycenter",
    "neptune": "neptune barycenter",
    "pluto": "pluto barycenter",
}

# Planet ID to human-readable name mapping for error messages and debugging
_PLANET_NAMES = {
    SUN: "Sun",
    MOON: "Moon",
    MERCURY: "Mercury",
    VENUS: "Venus",
    MARS: "Mars",
    JUPITER: "Jupiter",
    SATURN: "Saturn",
    URANUS: "Uranus",
    NEPTUNE: "Neptune",
    PLUTO: "Pluto",
    MEAN_NODE: "mean Node",
    TRUE_NODE: "true Node",
    MEAN_APOG: "mean Apogee",
    OSCU_APOG: "osc. Apogee",
    INTP_APOG: "intp. Apogee",
    INTP_PERG: "intp. Perigee",
    EARTH: "Earth",
    ISIS: "Isis-Transpluto",
    NIBIRU: "Nibiru",
    HARRINGTON: "Harrington",
    # Historical trans-Uranian/trans-Neptunian predictions carry a
    # parenthetical annotation naming the planet each prediction targeted,
    # matching the compatibility metadata contract (Le Verrier 1846 and Adams
    # 1845-46 for Neptune; Lowell 1915 and Pickering 1919 for Pluto).
    NEPTUNE_LEVERRIER: "Leverrier (Neptune)",
    NEPTUNE_ADAMS: "Adams (Neptune)",
    PLUTO_LOWELL: "Lowell (Pluto)",
    PLUTO_PICKERING: "Pickering (Pluto)",
    CUPIDO: "Cupido",
    HADES: "Hades",
    ZEUS: "Zeus",
    KRONOS: "Kronos",
    APOLLON: "Apollon",
    ADMETOS: "Admetos",
    # "Vulcanus" is the canonical Witte-Sieggrün (Hamburg School) spelling of
    # the seventh Uranian planet; the constant keeps its VULKANUS spelling for
    # API compatibility while the display name matches the metadata contract.
    VULKANUS: "Vulcanus",
    POSEIDON: "Poseidon",
    # Geocentric/heliocentric symbolic points the runtime models analytically
    # (ids 55-58). They resolve to their published astrological names because
    # calc_ut returns positions for them, so the metadata must agree.
    VULCAN: "Vulcan",
    WHITE_MOON: "Selena/White Moon",
    PROSERPINA: "Proserpina",
    WALDEMATH: "Waldemath",
    CHIRON: "Chiron",
    PHOLUS: "Pholus",
    CERES: "Ceres",
    PALLAS: "Pallas",
    JUNO: "Juno",
    VESTA: "Vesta",
}


def _cob_velocity_correction(barycenter_name: str, t):
    """Compute COB velocity correction via central difference.

    Args:
        barycenter_name: Name for lookup in moon_theories (e.g. 'jupiter barycenter')
        t: Skyfield Time object at which to evaluate

    Returns:
        numpy array of velocity offset in AU/day (3-element)
    """
    import numpy as np
    from .moon_theories import get_cob_offset
    from .state import get_timescale

    dt = 1.0 / 86400.0  # 1 second in days
    ts = get_timescale()
    t_prev = ts.tt_jd(t.tt - dt)
    t_next = ts.tt_jd(t.tt + dt)
    offset_prev = get_cob_offset(barycenter_name, t_prev)
    offset_next = get_cob_offset(barycenter_name, t_next)
    return (np.array(offset_next) - np.array(offset_prev)) / (2.0 * dt)


class _CobCorrectedTarget:
    """Wrapper that applies COB (Center of Body) correction to barycenter positions.

    DE440 returns system barycenter positions for outer planets (Jupiter, Saturn,
    Neptune, Pluto), but astrological calculations expect planet center positions. This
    wrapper applies analytical moon theory corrections to convert barycenter to
    center of body positions.

    The wrapper implements the Skyfield VectorFunction protocol so it can be used
    with observer.at(t).observe(target) seamlessly.
    """

    def __init__(self, barycenter, barycenter_name: str):
        """Initialize with a barycenter target and its name.

        Args:
            barycenter: Skyfield VectorFunction (e.g., planets['jupiter barycenter'])
            barycenter_name: Name like 'jupiter barycenter' for lookup in moon_theories
        """
        self._barycenter = barycenter
        self._barycenter_name = barycenter_name
        # Copy attributes needed by Skyfield
        self.center = getattr(barycenter, "center", 0)
        self.target = getattr(barycenter, "target", None)

    def at(self, t):
        """Return position at time t with COB correction applied.

        Args:
            t: Skyfield Time object

        Returns:
            Skyfield ICRF position with COB offset applied
        """
        import numpy as np
        from skyfield.positionlib import ICRF
        from .moon_theories import get_cob_offset

        bary_pos = self._barycenter.at(t)
        pos_au = bary_pos.position.au
        vel_au_per_d = bary_pos.velocity.au_per_d

        offset = get_cob_offset(self._barycenter_name, t)
        corrected_pos = pos_au + np.array(offset)
        corrected_vel = vel_au_per_d + _cob_velocity_correction(
            self._barycenter_name, t
        )

        return ICRF(corrected_pos, corrected_vel, t=t, center=self.center)

    def __repr__(self):
        return f"<CobCorrectedTarget {self._barycenter_name}>"

    def _observe_from_bcrs(self, observer):
        """Observe this target from an observer in BCRS coordinates.

        This method is called by Skyfield's observe() function to compute
        light-time corrected positions. We delegate to the barycenter's
        implementation and then apply the COB correction.

        Args:
            observer: Skyfield Barycentric position of the observer

        Returns:
            Tuple of (position_au, velocity_au_per_d, time, light_time_days)
        """
        import numpy as np
        from .moon_theories import get_cob_offset

        pos, vel, t, light_time = self._barycenter._observe_from_bcrs(observer)

        offset = get_cob_offset(self._barycenter_name, t)
        corrected_pos = pos + np.array(offset)
        corrected_vel = vel + _cob_velocity_correction(self._barycenter_name, t)

        return corrected_pos, corrected_vel, t, light_time


class _SpkCenterTarget:
    """Wrapper that combines barycenter position with SPK-based center offset.

    This class uses precise planet center segments from planet_centers.bsp to
    compute the offset from system barycenter to planet center. This provides
    high-precision planet-center positions while the segment covers the
    requested epoch.  Outside that coverage it deliberately returns the JPL
    system barycenter instead of synthesizing a center offset.

    The wrapper implements the Skyfield VectorFunction protocol so it can be used
    with observer.at(t).observe(target) seamlessly.
    """

    def __init__(
        self, barycenter, center_segment, planet_name: str, barycenter_name: str
    ):
        """Initialize with a barycenter target and center segment.

        Args:
            barycenter: Skyfield VectorFunction (e.g., planets['jupiter barycenter'])
            center_segment: Skyfield ChebyshevPosition segment from planet_centers.bsp
            planet_name: Planet name for debugging (e.g., 'jupiter')
        """
        self._barycenter = barycenter
        self._center_segment = center_segment
        self._planet_name = planet_name
        self._barycenter_name = barycenter_name
        # Copy attributes needed by Skyfield
        self.center = getattr(barycenter, "center", 0)
        self.target = getattr(barycenter, "target", None)
        self._last_center_epoch_jd: float | None = None
        self._last_used_center = False

    @staticmethod
    def _time_jd(t) -> float | None:
        """Return a scalar TDB JD while preserving unsupported-time fallback."""
        try:
            return float(t.whole) + float(t.tdb_fraction)
        except (AttributeError, TypeError, ValueError):
            return None

    def _center_is_forced_off(self) -> bool:
        """Return whether an enclosing speed stencil selected barycenter data."""
        return self._barycenter_name in _spk_center_forced_barycenters.get()

    def _record_center_choice(self, t, used_center: bool) -> None:
        """Remember the source used for the outer calculation's final state."""
        self._last_center_epoch_jd = self._time_jd(t)
        self._last_used_center = used_center

    def _coverage_contains(self, jd_tdb: float) -> bool:
        """Return whether any descriptor in the center target covers an epoch."""
        segments = getattr(self._center_segment, "segments", None)
        candidates = segments if segments is not None else (self._center_segment,)
        found_metadata = False
        for segment in candidates:
            coverage_segment = getattr(segment, "spk_segment", segment)
            start_jd = getattr(coverage_segment, "start_jd", None)
            end_jd = getattr(coverage_segment, "end_jd", None)
            if start_jd is None or end_jd is None:
                # Preserve the historical behavior for a custom VectorFunction
                # without exposed jplephem coverage metadata.
                return True
            found_metadata = True
            try:
                if float(start_jd) <= jd_tdb <= float(end_jd):
                    return True
            except (TypeError, ValueError):
                return True
        return not found_metadata

    def _speed_stencil_requires_fallback(self, half_step_days: float) -> bool:
        """Return whether both speed samples must use the system barycenter.

        The public speed is a central difference of two complete apparent-place
        evaluations. If one sample uses the planet-center SPK and the other is
        outside its coverage, the finite center/barycenter position offset is
        divided by a two-second span and appears as a huge sign-changing speed.
        Near a boundary, keep both samples on the same explicit barycenter
        fallback used outside coverage. The main position remains the physical
        center whenever its epoch is covered.
        """
        if not self._last_used_center or self._last_center_epoch_jd is None:
            return True
        # Retarded target time changes by essentially the observation-time
        # step; a 2*h guard covers the tiny light-time derivative as well as
        # floating-point endpoint rounding. Consecutive DAF descriptors are
        # treated as a union, so their shared internal boundary is safe.
        margin = 2.0 * half_step_days
        epoch = self._last_center_epoch_jd
        return not (
            self._coverage_contains(epoch - margin)
            and self._coverage_contains(epoch + margin)
        )

    def _segment_covers(self, t) -> bool:
        """Return whether the physical center segment covers a scalar time."""
        if self._center_is_forced_off():
            return False
        try:
            jd_tdb = float(t.whole) + float(t.tdb_fraction)
        except (AttributeError, TypeError, ValueError):
            return True
        return self._coverage_contains(jd_tdb)

    def at(self, t):
        """Return position at time t with SPK center offset applied.

        Args:
            t: Skyfield Time object

        Returns:
            Skyfield ICRF position with center offset applied
        """
        from skyfield.positionlib import ICRF

        # Get barycenter position (from SSB)
        bary_pos = self._barycenter.at(t)
        pos_au = bary_pos.position.au
        vel_au_per_d = bary_pos.velocity.au_per_d

        # A system barycenter is the explicit fallback when the JPL
        # planet-center segment does not cover this epoch.
        corrected_pos = pos_au
        corrected_vel = vel_au_per_d
        used_center = False
        if self._segment_covers(t):
            try:
                # Get center offset from SPK segment (relative to barycenter)
                center_offset = self._center_segment.at(t)
                offset_au = center_offset.position.au
                offset_vel = center_offset.velocity.au_per_d

                # Combine: barycenter + offset = center
                corrected_pos = pos_au + offset_au
                corrected_vel = vel_au_per_d + offset_vel
                used_center = True
            except _RANGE_ERRORS:
                pass
        self._record_center_choice(t, used_center)

        return ICRF(corrected_pos, corrected_vel, t=t, center=self.center)

    def __repr__(self):
        return f"<SpkCenterTarget {self._planet_name}>"

    def _observe_from_bcrs(self, observer):
        """Observe this target from an observer in BCRS coordinates.

        This method is called by Skyfield's observe() function to compute
        light-time corrected positions. We delegate to the barycenter's
        implementation and then apply the center offset.

        Args:
            observer: Skyfield Barycentric position of the observer

        Returns:
            Tuple of (position_au, velocity_au_per_d, time, light_time_days)
        """
        pos, vel, t, light_time = self._barycenter._observe_from_bcrs(observer)

        # Keep the retarded system-barycenter state returned above when no
        # physical JPL planet-center state is available.
        corrected_pos = pos
        corrected_vel = vel

        # Skyfield returns the observation time ``t`` alongside the barycenter
        # light time.  The center offset is a target state, so evaluate it when
        # the light left the target, not at the observer epoch.  Preserve
        # Skyfield's whole/fraction split to avoid quantizing the retarded epoch
        # through a single large JD float.
        retarded_t = t.ts.tdb_jd(t.whole, t.tdb_fraction - light_time)
        used_center = False
        if self._segment_covers(retarded_t):
            try:
                center_offset = self._center_segment.at(retarded_t)
                offset_au = center_offset.position.au
                offset_vel = center_offset.velocity.au_per_d

                corrected_pos = pos + offset_au
                corrected_vel = vel + offset_vel
                used_center = True
            except _RANGE_ERRORS:
                pass
        self._record_center_choice(retarded_t, used_center)

        return corrected_pos, corrected_vel, t, light_time


# Type alias for position result tuple
PositionResult = Tuple[float, float, float, float, float, float]


def _to_native_floats(values: tuple) -> PositionResult:
    """Convert numpy float types to native Python floats.

    Skyfield returns numpy.float64 values which can cause issues with:
    - JSON serialization
    - `is True` comparisons (speed < 0 returns np.bool_ instead of bool)
    - Type checking in downstream code

    This function ensures all values are native Python floats for
    compatibility with the reference API contract.

    Args:
        values: Tuple of 6 values (lon, lat, dist, speed_lon, speed_lat, speed_dist)

    Returns:
        Tuple of 6 native Python floats
    """
    return (
        float(values[0]),
        float(values[1]),
        float(values[2]),
        float(values[3]),
        float(values[4]),
        float(values[5]),
    )


def _apply_output_flags(result: PositionResult, iflag: int) -> PositionResult:
    """Apply output format flags (FLG_XYZ, FLG_RADIANS) to position result.

    Post-processes the standard (lon, lat, dist, dlon, dlat, ddist) output
    from _calc_body into the format requested by the caller.

    When FLG_XYZ is set, converts spherical coordinates to Cartesian:
        (lon°, lat°, dist_AU) → (x, y, z) in AU
        (dlon°/d, dlat°/d, ddist_AU/d) → (vx, vy, vz) in AU/day

    When FLG_RADIANS is set, converts angular values from degrees to radians.
    Distance values are unaffected.

    When both are set, FLG_XYZ takes priority (Cartesian output has no angles).

    Args:
        result: Standard 6-tuple (lon, lat, dist, dlon, dlat, ddist) in degrees/AU
        iflag: Calculation flags (may include FLG_XYZ, FLG_RADIANS)

    Returns:
        Transformed 6-tuple based on flags.
    """
    lon, lat, dist, dlon, dlat, ddist = result

    if iflag & FLG_XYZ:
        # Spherical(deg)→Cartesian+velocity via the shared helper, so this path
        # stays identical to the LEB/Skyfield, Horizons, and fixed-star XYZ
        # post-processing. Cast to native floats (inputs may be numpy.float64).
        from .fast_calc import _spherical_to_cartesian_with_velocity

        x, y, z, vx, vy, vz = _spherical_to_cartesian_with_velocity(
            lon, lat, dist, dlon, dlat, ddist
        )
        return (float(x), float(y), float(z), float(vx), float(vy), float(vz))

    if iflag & FLG_RADIANS:
        # Convert angular values from degrees to radians
        # Distance (index 2) and distance speed (index 5) are unchanged
        return (
            math.radians(lon),
            math.radians(lat),
            dist,
            math.radians(dlon),
            math.radians(dlat),
            ddist,
        )

    return result


def _degrade_jpl_echo_on_fallback(flags: int) -> int:
    """Echo FLG_SWIEPH for a JPLEPH request once the JPL file has fallen back.

    The reference, asked for FLG_JPLEPH with a JPL file it cannot open, computes
    from its default ephemeris and echoes FLG_SWIEPH. libephemeris mirrors that:
    when ``state.jpl_fallback_active()`` reports the configured JPL kernel is
    being silently replaced by the tier default, a requested JPLEPH bit is
    echoed (and normalized) as SWIEPH. The computation is unaffected -- it always
    uses the JPL DE kernel; this only keeps retflag 1:1 with the reference. The
    state check is skipped unless JPLEPH was actually requested, so the common
    path pays no filesystem cost.
    """
    from .constants import FLG_JPLEPH

    if not (flags & FLG_JPLEPH):
        return flags
    from . import state

    if state.jpl_fallback_active():
        flags = (flags & ~FLG_JPLEPH) | FLG_SWIEPH
    return flags


def consume_non_echoed_flags(flags: int, planet: int) -> int:
    """Drop request bits this library consumes rather than echoes back.

    Shared by every calc entry point — module-level ``calc_ut``/``calc`` and
    the ``EphemerisContext`` equivalents — because the strip has to happen
    BEFORE the raw request is captured for the retflag echo. Living only in
    the module functions, it made the two entry points disagree: a context
    call echoed bits the module call consumed.

    * ``FLG_CENTER_BODY`` for Sun..Mars: ipl 0-4 have no satellite-system
      barycenter to resolve (positions are body centers already) and the
      measured reference retflag does not echo the bit for these five
      bodies, while nodes, apogees and asteroids echo it back unchanged.
    * ``FLG_JPLHOR`` / ``FLG_JPLHOR_APPROX``: this library performs no
      JPL-Horizons dpsi/deps Earth-orientation reduction (the flags are
      accepted for API compatibility only, see docs/reference/flags.md), so
      echoing them would advertise a correction that was never applied.
      An external implementation consumes them the same way for
      SWIEPH/default requests but, on an explicit ``FLG_JPLEPH`` request,
      applies its approximate reduction (measured: -0.048" on Mars,
      -0.050" on the Moon at J2000) and echoes ``FLG_JPLHOR_APPROX`` to
      say so. Reproducing that echo without the reduction would make the
      retflag lie about what was computed, and the reduction itself is not
      reconstructible from published models — see
      docs/comparison/intentional-divergences.md.
    """
    from .constants import FLG_JPLHOR, FLG_JPLHOR_APPROX

    if (flags & FLG_CENTER_BODY) and 0 <= planet <= MARS:
        flags &= ~FLG_CENTER_BODY
    if flags & (FLG_JPLHOR | FLG_JPLHOR_APPROX):
        flags &= ~(FLG_JPLHOR | FLG_JPLHOR_APPROX)
    return flags


def _normalize_calc_flags(flags: int) -> int:
    """Normalize calculation flags according to the public API contract.

    Shared by the module-level calc()/calc_ut() and the EphemerisContext
    entry points so both expose the same compatibility behavior:

    - FLG_MOSEPH is accepted for compatibility but stripped (all
      calculations use the JPL DE440/DE441 path).
    - Exactly one ephemeris bit is echoed. Selection bits are mutually
      exclusive, with priority JPLEPH > SWIEPH; FLG_SWIEPH is the default. A
      requested JPLEPH degrades to SWIEPH when the configured JPL file has
      silently fallen back to the tier default (see
      _degrade_jpl_echo_on_fallback), matching the reference.
    - FLG_SPEED3 implies FLG_SPEED for the computation, while retaining the
      SPEED3 marker internally so eligible bodies can use the public centered
      three-position derivative. SPEED wins when callers pass both bits.
    """
    from .constants import FLG_JPLEPH, FLG_MOSEPH

    flags = flags & ~FLG_MOSEPH
    flags = _degrade_jpl_echo_on_fallback(flags)
    if flags & FLG_JPLEPH:
        # JPLEPH takes priority; ensure the ephemeris bits stay mutually
        # exclusive so retflag never echoes both (the reference never does).
        flags &= ~FLG_SWIEPH
    elif not (flags & FLG_SWIEPH):
        flags |= FLG_SWIEPH
    if flags & FLG_SPEED3:
        if flags & FLG_SPEED:
            flags &= ~FLG_SPEED3
        else:
            flags |= FLG_SPEED
    return _resolve_center_flags(flags)


def _resolve_center_flags(flags: int) -> int:
    """Resolve conflicting observation-center flags to compatibility priority.

    When more than one center is requested, keep exactly one with priority
    TOPOCTR > BARYCTR > HELCTR > geocentric and strip losing bits from both
    computation and echoed flags.
    """
    if flags & FLG_TOPOCTR:
        flags &= ~(FLG_HELCTR | FLG_BARYCTR)
    elif flags & FLG_BARYCTR:
        flags &= ~FLG_HELCTR
    return flags


def _validate_barycentric_moseph(flags: int, entry_point: str) -> None:
    """Reject the unsupported barycentric-MOSEPH request combination.

    This public validation rule applies to the raw request before body dispatch.
    Center priority still applies (TOPOCTR removes BARYCTR), as does ephemeris
    priority (JPLEPH or SWIEPH wins over MOSEPH when supplied).
    """
    from .constants import FLG_JPLEPH, FLG_MOSEPH

    resolved = _resolve_center_flags(flags)
    if (
        resolved & FLG_BARYCTR
        and resolved & FLG_MOSEPH
        and not (resolved & (FLG_JPLEPH | FLG_SWIEPH))
    ):
        from .exceptions import Error

        raise Error(
            f"{entry_point}: barycentric analytical positions are not supported."
        )


def _echo_request_bits(retflag: int, raw_flags: int) -> int:
    """Restore public retflag echoes for the raw request.

    _normalize_calc_flags adjusts bits for the computation in two ways that
    must not leak into the echoed retflag:

    - FLG_SPEED3 implies FLG_SPEED internally, but a SPEED3-only request echoes
      SPEED3; an explicit SPEED bit takes priority when both are present.
    - FLG_MOSEPH is stripped from the computation flags (every calculation
      uses the JPL DE440/DE441 path), but the echo retains one ephemeris bit
      with priority JPLEPH > SWIEPH > MOSEPH. This matches the ECL_NUT and
      calc_pctr paths, which retain MOSEPH via _exclusive_ephemeris_bit.

    Both corrections touch only the ECHO — computation flags stay as
    normalized.
    """
    from .constants import FLG_JPLEPH, FLG_MOSEPH

    if (raw_flags & FLG_SPEED3) and not (raw_flags & FLG_SPEED):
        retflag = (retflag & ~FLG_SPEED) | FLG_SPEED3
    if (raw_flags & FLG_MOSEPH) and not (raw_flags & (FLG_JPLEPH | FLG_SWIEPH)):
        retflag = (retflag & ~(FLG_SWIEPH | FLG_JPLEPH)) | FLG_MOSEPH
    return retflag


def _implied_retflag_bits(flags: int) -> int:
    """Retflag bits implied by public flag normalization.

    Beyond echoing the caller's bits, the API echoes implied flags:
    J2000 and SIDEREAL output is referred to a mean equinox, so FLG_NONUT is
    echoed; heliocentric, barycentric and true-position output skips light
    deflection and annual aberration, so FLG_NOGDEFL | FLG_NOABERR are echoed.

    Only the ECHOED retflag carries these bits — the computation flags must
    stay untouched (adding NONUT to a sidereal request would switch the
    ayanamsha realization used by the calculation).
    """
    extra = 0
    if flags & (FLG_J2000 | FLG_SIDEREAL):
        extra |= FLG_NONUT
    if flags & (FLG_HELCTR | FLG_BARYCTR | FLG_TRUEPOS):
        extra |= FLG_NOGDEFL | FLG_NOABERR
    return extra


def _calc_tt_epheflag_echo(retflag: int, raw_flags: int) -> int:
    """Apply the public calc() (TT) ephemeris-bit echo convention.

    ``calc()`` echoes only the ephemeris-selection bits the
    caller actually passed (calc(jd, body, 0) echoes retflag 0), while
    calc_ut() adds the default SWIEPH bit. Strip the default-injected
    ephemeris bit when the raw request carried none. Echo-only: the
    computation keeps the normalized (defaulted) flags.
    """
    from .constants import FLG_JPLEPH, FLG_MOSEPH

    if not (raw_flags & (FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH)):
        return retflag & ~(FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH)
    return retflag


def _exclusive_ephemeris_bit(flags: int) -> int:
    """Force exactly one ephemeris bit in the flag word.

    ``calc_pctr()`` uses mutually exclusive ephemeris bits, with FLG_SWIEPH as
    the default and priority JPLEPH > SWIEPH > MOSEPH. Unlike calc()/calc_ut(),
    FLG_MOSEPH is not stripped and FLG_SPEED3 is not remapped here.
    """
    from .constants import FLG_JPLEPH, FLG_MOSEPH

    flags = _degrade_jpl_echo_on_fallback(flags)
    ephmask = FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH
    if flags & FLG_JPLEPH:
        epheflag = FLG_JPLEPH
    elif flags & FLG_SWIEPH:
        epheflag = FLG_SWIEPH
    elif flags & FLG_MOSEPH:
        epheflag = FLG_MOSEPH
    else:
        epheflag = FLG_SWIEPH
    return (flags & ~ephmask) | epheflag


def _strip_output_flags(flags: int) -> int:
    """Remove output-format bits (FLG_XYZ, FLG_RADIANS) from calc flags.

    They are output-format flags, not calculation flags: entry points strip
    them before the calculation and re-apply them on the result via
    _finalize_output_flags(). Keeping the pair in shared helpers prevents an
    entry point from echoing the flags in retflag without honoring them.
    """
    return flags & ~FLG_XYZ & ~FLG_RADIANS


def _finalize_output_flags(
    pos: PositionResult, retflag: int, flags: int
) -> Tuple[PositionResult, int]:
    """Apply output-format flags to a result and echo them in retflag.

    Also echoes the flag bits the reference's flag normalization implies
    (see _implied_retflag_bits) so retflag stays 1:1 with the reference.
    """
    return (
        _apply_output_flags(pos, flags),
        retflag | (flags & (FLG_XYZ | FLG_RADIANS)) | _implied_retflag_bits(flags),
    )


def _rebuild_mean_point_speed3(
    current: PositionResult,
    t,
    planet: int,
    flags: int,
    calc_fn,
    *,
    timescale=None,
) -> PositionResult:
    """Differentiate a mean lunar point in its final requested frame.

    SPEED3 is a three-position numerical derivative. Sampling the complete
    no-speed reduction at TT ``t±h`` keeps nutation, sidereal, coordinate, and
    output-format transformations inside the same centered stencil.
    """
    if not (flags & FLG_SPEED3) or planet not in (MEAN_NODE, MEAN_APOG):
        return current

    from .fast_calc import _SPEED3_STEP_DAYS, _apply_central_speed3_stencil

    if timescale is None:
        timescale = get_timescale()
    sample_flags = flags & ~(FLG_SPEED | FLG_SPEED3)
    sample_calc_flags = _strip_output_flags(sample_flags)

    def _sample(offset: float) -> PositionResult:
        sample_t = timescale.tt_jd(float(t.tt), offset)
        sample_pos, sample_retflag = calc_fn(sample_t, planet, sample_calc_flags)
        return _finalize_output_flags(sample_pos, sample_retflag, sample_flags)[0]

    previous = _sample(-_SPEED3_STEP_DAYS)
    following = _sample(_SPEED3_STEP_DAYS)
    return _apply_central_speed3_stencil(current, previous, following, flags)


def _run_pctr_pipeline(calc_fn, flags: int) -> Tuple[PositionResult, int]:
    """normalize-bits -> strip -> compute -> finalize, shared by calc_pctr().

    Both the module-level calc_pctr() and EphemerisContext.calc_pctr()
    route through this so their flag semantics cannot drift apart;
    ``calc_fn`` receives the stripped calculation flags.
    """
    flags = _exclusive_ephemeris_bit(flags)
    pos, retflag = calc_fn(_strip_output_flags(flags))
    return _finalize_output_flags(pos, retflag, flags)


def _south_node_from_north(
    north_result: tuple[float, float, float, float, float, float], flags: int
) -> tuple[float, float, float, float, float, float]:
    """Derive the descending ("south") lunar node from the ascending node.

    The south node is the antipode of the north node on the celestial sphere,
    so it is obtained from the already-formatted north result.  The transform
    depends on the output representation:

    * ``FLG_XYZ`` (Cartesian): negate the position and velocity vectors.
    * ``FLG_RADIANS`` (spherical radians): longitude + pi (mod 2*pi), with the
      latitude and latitude-speed negated.
    * default (spherical degrees): longitude + 180 (mod 360), with the
      latitude and latitude-speed negated.

    The previous code always did ``(north[0] + 180.0) % 360.0``; under
    ``FLG_XYZ`` that added 180 to a Cartesian x-component, and under
    ``FLG_RADIANS`` it added 180 to a radian longitude, both producing garbage.
    """
    if flags & FLG_XYZ:
        # Cartesian: the antipode is the negated position (and velocity) vector.
        return (
            -north_result[0],
            -north_result[1],
            -north_result[2],
            -north_result[3],
            -north_result[4],
            -north_result[5],
        )
    half_turn = math.pi if (flags & FLG_RADIANS) else 180.0
    full_turn = 2.0 * math.pi if (flags & FLG_RADIANS) else 360.0
    return (
        (north_result[0] + half_turn) % full_turn,
        -north_result[1],
        north_result[2],
        north_result[3],
        -north_result[4],
        north_result[5],
    )


def _body_uses_jpl_ephemeris(ipl: int) -> bool:
    """Check if a body uses the JPL ephemeris for calculations.

    Bodies that use the JPL ephemeris are subject to the ephemeris date range
    limitations (e.g., DE440 covers 1550-2650). Other bodies like lunar nodes,
    Lilith, hypothetical planets, and minor bodies with Keplerian elements
    are calculated mathematically and don't depend on the JPL ephemeris range.

    Args:
        ipl: Planet/body ID

    Returns:
        True if the body uses JPL ephemeris, False otherwise
    """
    # Standard planets use JPL ephemeris
    return ipl in _PLANET_MAP


def _wrap_ephemeris_range_error(
    skyfield_error: Exception,
    jd: float,
    body_id: int | None = None,
) -> "EphemerisRangeError":
    """
    Wrap Skyfield's EphemerisRangeError with enhanced error details.

    Creates a libephemeris.EphemerisRangeError with detailed information about:
    - The requested Julian Day number
    - The supported date range in both JD and calendar format
    - The body being calculated (if known)
    - The ephemeris file in use

    Args:
        skyfield_error: The original Skyfield EphemerisRangeError
        jd: The Julian Day that was requested
        body_id: The body ID being calculated (optional)

    Returns:
        EphemerisRangeError with enhanced error message
    """
    from .exceptions import EphemerisRangeError
    from . import state

    # Get ephemeris info
    path, start_jd, end_jd, denum = state.get_current_file_data(0)

    # Extract date strings from the original error message if available
    original_msg = str(skyfield_error)
    start_date = None
    end_date = None

    # Try to parse dates from "ephemeris segment only covers dates YYYY-MM-DD through YYYY-MM-DD"
    import re

    date_match = re.search(
        r"covers dates (\d{4}-\d{2}-\d{2}) through (\d{4}-\d{2}-\d{2})", original_msg
    )
    if date_match:
        start_date = date_match.group(1)
        end_date = date_match.group(2)

    # Convert requested JD to calendar date for the message
    from .time_utils import revjul

    req_year, req_month, req_day, req_hour = revjul(jd, 1)  # Gregorian calendar
    req_date_str = f"{req_year}-{req_month:02d}-{req_day:02d}"

    # Get body name if available
    body_name = None
    if body_id is not None:
        body_name = get_planet_name(body_id)

    # Get ephemeris filename
    ephemeris_file = None
    if path:
        import os

        ephemeris_file = os.path.basename(path)
    elif denum:
        ephemeris_file = f"de{denum}.bsp"

    # Build enhanced error message
    parts = []

    if body_name and body_id is not None:
        parts.append(f"Cannot calculate {body_name} (ID {body_id})")
    else:
        parts.append("Calculation failed")

    parts.append(f"for JD {jd:.6f} ({req_date_str}):")
    parts.append("date is outside ephemeris range.")

    if start_jd and end_jd:
        parts.append(f"\n  Supported range: JD {start_jd:.1f} to {end_jd:.1f}")
        if start_date and end_date:
            parts.append(f" ({start_date} to {end_date})")
    elif start_date and end_date:
        parts.append(f"\n  Supported range: {start_date} to {end_date}")

    if ephemeris_file:
        parts.append(f"\n  Ephemeris file: {ephemeris_file}")

    message = " ".join(parts[:3]) + "".join(parts[3:])

    return EphemerisRangeError(
        message=message,
        requested_jd=jd,
        start_jd=start_jd if start_jd else None,
        end_jd=end_jd if end_jd else None,
        start_date=start_date,
        end_date=end_date,
        body_id=body_id,
        body_name=body_name,
        ephemeris_file=ephemeris_file,
    )


# NAIF IDs for planet centers
_PLANET_CENTER_NAIF_IDS = {
    "jupiter": 599,
    "saturn": 699,
    "uranus": 799,
    "neptune": 899,
    "pluto": 999,
}


def get_planet_target(planets, target_name: str):
    """
    Get planet target from ephemeris with fallback to barycenter.

    For outer planets (Jupiter, Saturn, Uranus, Neptune, Pluto), this function
    first checks if planet_centers.bsp is available. If so, it returns a
    _SpkCenterTarget that combines the barycenter position with the precise
    SPK-based center offset.

    If no planet-center segment is available, it explicitly falls back to the
    system barycenter from the main JPL ephemeris.

    Args:
        planets: Skyfield SpiceKernel ephemeris object
        target_name: Planet name from _PLANET_MAP (e.g., 'jupiter', 'saturn')

    Returns:
        Skyfield planet object (or wrapper with center offset)

    Raises:
        KeyError: If neither planet center nor barycenter found in ephemeris
    """
    from .state import get_planet_center_segment

    # For outer planets, try to use SPK-based planet centers
    if target_name in _PLANET_CENTER_NAIF_IDS:
        naif_id = _PLANET_CENTER_NAIF_IDS[target_name]
        center_segment = get_planet_center_segment(naif_id)

        if center_segment is not None:
            # We have the planet center SPK segment
            # Get the barycenter from the main ephemeris
            barycenter_name = _PLANET_FALLBACK.get(target_name)
            if barycenter_name:
                try:
                    barycenter = planets[barycenter_name]
                    return _SpkCenterTarget(
                        barycenter, center_segment, target_name, barycenter_name
                    )
                except KeyError:
                    pass

    # Standard path: try a center already present in the main ephemeris, then
    # retain the JPL system barycenter as the explicit fallback.
    try:
        return planets[target_name]
    except KeyError:
        if target_name in _PLANET_FALLBACK:
            barycenter_name = _PLANET_FALLBACK[target_name]
            return planets[barycenter_name]
        raise


def get_planet_name(planet: int) -> str:
    """
    Get the human-readable name of a planet given its ID.

    Useful for error messages and debugging output.

    Args:
        planet: Planet/body ID (SUN, MOON, etc.)

    Returns:
        Human-readable planet name as a string.
        Returns an empty string for unrecognized planet IDs, matching the
        compatibility API.

    Example:
        >>> get_planet_name(0)
        'Sun'
        >>> get_planet_name(1)
        'Moon'
        >>> get_planet_name(4)
        'Mars'
    """
    if planet in _PLANET_NAMES:
        return _PLANET_NAMES[planet]
    if planet > AST_OFFSET:
        # Numbered minor planets: the reference resolves the built-in
        # set by asteroid number and returns an empty string for the
        # rest (names normally come from the asteroid ephemeris files).
        _BUILTIN_AST_NAMES = {
            1: "Ceres",
            2: "Pallas",
            3: "Juno",
            4: "Vesta",
            2060: "Chiron",
            5145: "Pholus",
        }
        builtin = _BUILTIN_AST_NAMES.get(planet - AST_OFFSET, "")
        if builtin:
            return builtin
        # Registry-served exotic bodies resolve to their catalog names, the
        # same names the LEB/SPK pipeline serves them under.
        from .exotic_bodies import exotic_display_name

        return exotic_display_name(planet)
    from . import planetary_moons

    if planetary_moons.is_planetary_moon(planet):
        name = planetary_moons.get_moon_name(planet)
        if not name.startswith("Unknown"):
            return name
    return ""


# The reference ephemeris treats calc_ut(jd, AST_OFFSET + 1, flags) identically to
# calc_ut(jd, CERES, flags); same for the other built-in asteroids.
_AST_NUMBER_REMAP = {
    1: CERES,  # 17
    2: PALLAS,  # 18
    3: JUNO,  # 19
    4: VESTA,  # 20
    2060: CHIRON,  # 15
    5145: PHOLUS,  # 16
}


def _remap_ast_offset(ipl: int) -> int:
    """Map AST_OFFSET + N to the dedicated body id for built-in asteroids.

    Must run before the LEB/Horizons dispatch so e.g. AST_OFFSET + 1 is
    served from the same backend as CERES; other ids pass through.
    """
    if ipl > AST_OFFSET:
        return _AST_NUMBER_REMAP.get(ipl - AST_OFFSET, ipl)
    return ipl


def _nutation_rate_deg_per_day(jd_tt: float, dt: float = 0.5) -> float:
    """Central-difference rate of nutation in longitude (deg/day).

    The lunar node/apse bodies are output on the true ecliptic of date,
    so their reference-compatible speeds include d(dpsi)/dt (about
    2e-6 deg/day at J2000); under FLG_NONUT the rate is excluded.
    """
    from .cache import get_cached_nutation

    dpsi_prev, _ = get_cached_nutation(jd_tt - dt)
    dpsi_next, _ = get_cached_nutation(jd_tt + dt)
    return math.degrees(dpsi_next - dpsi_prev) / (2.0 * dt)


def _add_of_date_nutation(
    lon: float, dlon: float, jd_tt: float, iflag: int
) -> Tuple[float, float]:
    """Rotate a mean-ecliptic-of-date longitude onto the true ecliptic of date.

    Adds nutation in longitude (Δψ) to ``lon`` and, when FLG_SPEED is set, the
    nutation-in-longitude rate dΔψ/dt to ``dlon`` -- the same treatment the
    lunar-node/apse and predicted-planet paths use so every of-date body is
    reported on the true ecliptic, matching the reference API (which applies Δψ
    to all bodies uniformly).

    The rotation is skipped -- leaving the mean ecliptic -- for the
    nutation-free frames: FLG_NONUT, FLG_J2000 (the J2000 ecliptic carries no
    nutation), and SIDEREAL+EQUATORIAL (rotated to the equator with the mean
    obliquity, so the longitude must stay mean too). For plain SIDEREAL the Δψ
    is added here and the true ayanamsha (mean + Δψ) is subtracted downstream,
    leaving the intrinsic sidereal longitude.

    Args:
        lon: Ecliptic longitude in degrees (mean ecliptic of date).
        dlon: Longitude speed in degrees/day.
        jd_tt: Julian Day in Terrestrial Time.
        iflag: Calculation flags bitmask.

    Returns:
        ``(lon, dlon)`` on the true ecliptic of date, or unchanged in a
        nutation-free frame.
    """
    _sid_eq = bool(iflag & FLG_SIDEREAL) and bool(iflag & FLG_EQUATORIAL)
    if bool(iflag & (FLG_NONUT | FLG_J2000)) or _sid_eq:
        return lon, dlon

    from .cache import get_cached_nutation

    dpsi_rad, _ = get_cached_nutation(jd_tt)
    lon = (lon + math.degrees(dpsi_rad)) % 360.0
    if iflag & FLG_SPEED:
        dlon += _nutation_rate_deg_per_day(jd_tt)
    return lon, dlon


def _strict_source_better_than_keplerian(ipl: int) -> bool:
    """Return True if a better-than-Keplerian source is available for ``ipl``.

    Used by the strict-precision gate to decide whether to refuse the Keplerian
    fallback.  The only fallback consulted here is ASSIST (local N-body data):
    by the time the gate is reached the registered SPK kernel has already
    declined to serve the epoch.  Availability is a purely local filesystem /
    import probe -- it never touches the network and never depends on whether
    an (irrelevant) auto-download was attempted.

    Args:
        ipl: libephemeris body ID.

    Returns:
        True if ASSIST can integrate this body (curated elements present and
        the ASSIST package + data files are installed), False otherwise.
    """
    from . import minor_bodies

    if ipl not in minor_bodies.MINOR_BODY_ELEMENTS:
        return False
    try:
        from .rebound_integration import check_assist_data_available

        return bool(check_assist_data_available())
    except (ImportError, RuntimeError, ValueError, FileNotFoundError):
        return False


def _spk_out_of_coverage_error(ipl: int, spk_info, jd_tt: float):
    """Build the strict-mode error for a registered kernel that excludes ``jd_tt``.

    Unlike :meth:`SPKRequiredError.for_body` -- which instructs the caller to
    *register* an SPK -- this variant is raised when a kernel IS registered for
    the body but this epoch falls outside its coverage window and no other
    high-precision source (ASSIST) is available.  The message names the
    registered kernel and its coverage so the remedy is obvious: widen the
    coverage, not register a kernel from scratch.

    Args:
        ipl: libephemeris body ID.
        spk_info: ``(spk_file, naif_id)`` for the registered kernel.
        jd_tt: The requested epoch (JD TT) that fell outside coverage.

    Returns:
        SPKRequiredError with an accurate, coverage-citing message.
    """
    import os

    from . import spk
    from .exceptions import SPKRequiredError

    spk_file = spk_info[0]
    body_name = spk._get_body_name(ipl) or str(ipl)
    try:
        coverage = spk.get_spk_coverage(spk_file)
    except (OSError, ValueError, KeyError, RuntimeError):
        coverage = None

    lines = [
        f"SPK kernel for {body_name} does not cover JD {jd_tt:.1f} TT "
        "(strict precision mode is enabled).",
        "",
        f"A kernel is already registered ({os.path.basename(spk_file)}),",
    ]
    if coverage is not None:
        lines.append(
            f"but its coverage is [{coverage[0]:.1f}, {coverage[1]:.1f}] TT, "
            "which excludes this epoch."
        )
    else:
        lines.append("but its coverage does not include this epoch.")
    lines += [
        "",
        "The available Keplerian approximation is body- and date-dependent "
        "and is not classified as ephemeris-grade.",
        "For accurate calculations at this epoch, you must either:",
        "",
        "1. Register a kernel whose coverage includes this epoch,",
        "2. Outside sealed mode, enable automatic SPK download:",
        "   >>> eph.set_auto_spk_download(True)",
        "3. Install ASSIST N-body data for a high-precision fallback, or",
        "4. Disable strict precision mode (not recommended):",
        "   >>> eph.set_strict_precision(False)",
    ]
    return SPKRequiredError(message="\n".join(lines), body_id=ipl, body_name=body_name)


def _try_auto_spk_download(t, ipl: int, iflag: int):
    """
    Try to automatically download and use SPK for a minor body.

    This is called when auto SPK download is enabled and no SPK is registered
    for the given body. It downloads the SPK from JPL Horizons via direct HTTP
    (no external dependencies beyond urllib) and then calculates the position.

    Args:
        t: Skyfield Time object
        ipl: Planet/body ID
        iflag: Calculation flags

    Returns:
        Position tuple (lon, lat, dist, speed_lon, speed_lat, speed_dist)
        or None if download fails or body not supported.

    Note:
        This function uses spk.download_and_register_spk() which communicates
        directly with JPL Horizons API via HTTP. No astroquery dependency.
    """
    from . import minor_bodies, spk
    from .constants import SPK_AUTO_DOWNLOAD_BLOCKED, SPK_BODY_NAME_MAP
    from .logging_config import get_logger

    logger = get_logger()

    # Skip bodies where JPL blocks SPK generation
    if ipl in SPK_AUTO_DOWNLOAD_BLOCKED:
        return None

    if ipl in SPK_BODY_NAME_MAP:
        horizons_id, naif_id = SPK_BODY_NAME_MAP[ipl]
    elif AST_OFFSET < ipl < FIXSTAR_OFFSET:
        # Arbitrary numbered asteroid: Horizons small-body record syntax
        horizons_id = f"{ipl - AST_OFFSET};"
    else:
        return None

    try:
        # Use the date range from the current precision tier so the downloaded
        # SPK file matches the user's configured coverage level. The file
        # is cached on disk and will be reused on subsequent calls.
        from .state import get_spk_date_range_for_tier

        start_date, end_date = get_spk_date_range_for_tier()
        request_start = max(
            _pad_iso_date(start_date, -45),
            minor_bodies.HORIZONS_SPK_DATE_MIN,
        )
        request_end = min(
            _pad_iso_date(end_date, 45),
            minor_bodies.HORIZONS_SPK_DATE_MAX,
        )

        body_name = spk._get_body_name(ipl) or horizons_id
        logger.info("Auto-downloading SPK for %s from JPL Horizons...", body_name)

        # Download and register using direct HTTP (no astroquery needed).
        # Don't pass naif_id — let download_and_register_spk auto-detect from
        # the SPK file, since JPL Horizons uses the 20000000+N convention while
        # our constants may use the older 2000000+N convention.
        spk.download_and_register_spk(
            body=horizons_id,
            ipl=ipl,
            start=request_start,
            end=request_end,
        )

        # Now try to calculate using the newly registered SPK
        return spk.calc_spk_body_position(t, ipl, iflag)

    except (OSError, ValueError, KeyError, RuntimeError, TypeError) as e:
        logger.warning("Auto SPK download failed for body %d: %s", ipl, e)
        return None
    except EphemerisRangeError as e:
        # The freshly downloaded kernel does not cover the requested epoch;
        # fall back to the Keplerian path like any other SPK miss.
        logger.warning("Auto SPK coverage miss for body %d: %s", ipl, e)
        return None


def _pad_iso_date(date_s: str, days: int) -> str:
    """Return ``date_s`` shifted by ``days`` days (YYYY-MM-DD)."""
    from datetime import datetime, timedelta

    dt = datetime.strptime(date_s, "%Y-%m-%d").date()
    return (dt + timedelta(days=days)).isoformat()


def _is_transient_close_race(exc: BaseException) -> bool:
    """Return True for the transient errors a concurrent close()/reload raises.

    A concurrent ``close()`` (module state or an :class:`EphemerisContext`)
    empties the SPK kernel or LEB reader while a calculation is mid-read. The
    in-flight lookup then fails, before the resources reload lazily, with one
    of a small set of shapes:

    - a ``ValueError`` whose message contains ``"closed file"`` (e.g.
      ``'seek of closed file'``) — an SPK/mmap handle closed underneath the
      read;
    - any ``KeyError``. The known race shapes are a transiently emptied
      Skyfield kernel (``"kernel 'de440.bsp' is missing 'EARTH'"``) or LEB
      reader (``"Body 14 not in any installed LEB tier"``, ``"not in LEB
      file"``, a bare channel-id key such as ``KeyError(0)``), but the
      predicate deliberately matches every ``KeyError``: a race surfacing
      through a new lookup shape must stay retryable rather than become an
      intermittent concurrency failure. A genuine (non-race) ``KeyError``
      is retried and then propagates unchanged once the retry budget is
      exhausted.

    These are races, not genuine coverage/body misses: those are converted to
    typed :class:`EphemerisRangeError` / :class:`UnknownBodyError` upstream
    (``_raise_leb_range_miss`` / ``validate_jd_range``) before a raw error can
    reach the ``_calc_body`` call. One retry after the lazy reload resolves the
    race; a genuine residual error re-raises on the retry. Other ``ValueError``
    shapes (e.g. an out-of-range ``"Invalid Time"``) are deliberately excluded
    so they still reach the range-error converter.
    """
    if isinstance(exc, ValueError):
        return "closed file" in str(exc)
    return isinstance(exc, KeyError)


# A single reload usually resolves a close()/reload race, but under sustained
# concurrent close() pressure the freshly reloaded kernel/reader can be emptied
# again before the retry completes. Retry a small bounded number of times so a
# multi-thread burst of close() calls still converges instead of surfacing a
# transient error; a genuine (non-race) error still propagates immediately.
_MAX_CLOSE_RACE_RETRIES = 8


def _calc_body_race_safe(t, planet: int, calc_iflag: int):
    """Call ``_calc_body``, transparently retrying a concurrent close() race.

    See :func:`_is_transient_close_race`. Non-race errors propagate on their
    first occurrence; a transient close()/reload error is retried after the
    lazy reload, up to ``_MAX_CLOSE_RACE_RETRIES`` times.
    """
    for attempt in range(_MAX_CLOSE_RACE_RETRIES + 1):
        try:
            return _calc_body(t, planet, calc_iflag)
        except (ValueError, KeyError) as _race_err:
            if (
                not _is_transient_close_race(_race_err)
                or attempt == _MAX_CLOSE_RACE_RETRIES
            ):
                raise
    raise AssertionError("unreachable")  # pragma: no cover


# Sidereal modes for which SIDBIT_ECL_T0 / SIDBIT_SSY_PLANE are inert.
#
# These are the star- and galactic-frame "true"/live-catalog modes whose zero
# point is *defined* by the actual direction of a star or a galactic-frame axis
# projected onto the ecliptic (or equator) of date. That direction already
# fixes the frame in which the sidereal longitude is measured, so re-projecting
# it onto a different reference plane -- the mean ecliptic/equinox of the mode's
# t0 (SIDBIT_ECL_T0) or the solar-system invariable plane (SIDBIT_SSY_PLANE) --
# is not defined by the underlying model. Measured reference behavior confirms
# the projection is a no-op for exactly these modes (the sidereal longitude is
# identical to the un-projected baseline), while the epoch-anchored mean modes
# (Fagan/Bradley, Lahiri, ...) and the two live modes whose zero is not a bare
# ecliptic-of-date star/frame direction -- Aldebaran = 15 Tau (14, an absolute
# ecliptic longitude anchor) and Vettius Valens (42, a Method-B epoch mode) --
# do accept the projection. This is the same "calculated" star/galactic class
# that drops FLG_NONUT from get_ayanamsa_ex (see _AYANAMSA_EX_NONUT_DROP_MODES);
# Skydram/Mardyks (34) suppresses it too but reaches the same result through the
# fixed-epoch frame path, so it need not be listed here.
_SIDBIT_PROJECTION_SUPPRESS_MODES = frozenset(
    {
        SIDM_GALCENT_0SAG,
        SIDM_TRUE_CITRA,
        SIDM_TRUE_REVATI,
        SIDM_TRUE_PUSHYA,
        SIDM_GALCENT_RGILBRAND,
        SIDM_GALEQU_IAU1958,
        SIDM_GALEQU_TRUE,
        SIDM_GALEQU_MULA,
        SIDM_TRUE_MULA,
        SIDM_GALCENT_MULA_WILHELM,
        SIDM_TRUE_SHEORAN,
        SIDM_GALCENT_COCHRANE,
    }
)


def _sidbit_projection_calc(
    tjd: float,
    planet: int,
    flags: int,
    raw_flags: int,
    sid_mode: int,
    sid_bits: int,
    calc_fn,
    tt_echo: bool = False,
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """Compute a SIDBIT-projected *ecliptic* sidereal position.

    Rewrites the request to FLG_J2000|FLG_NONUT (ecliptic), runs it through the
    active backend via ``calc_fn`` (calc_ut or calc), then rotates the result
    onto the projection frame and removes the mode's sidereal zero point:

      * SIDBIT_SSY_PLANE -> the solar-system invariable plane (Souami &
        Souchay 2012); the zero point is the J2000 ayanamsha direction
        projected onto that plane.
      * SIDBIT_ECL_T0 -> the mean ecliptic and equinox of the mode's reference
        epoch t0; the zero point is the mode's ayanamsha at t0.

    Ecliptic output (spherical or FLG_XYZ) is projected onto the plane.
    Equatorial output under SIDBIT_ECL_T0 is the position reduced to the
    MEAN EQUATOR AND EQUINOX of the mode's t0 — compatibility contract:
    the J2000|NONUT sub-request is precessed to t0 with no ayanamsha
    subtraction, and the retflag keeps the plain sidereal echo (the same
    convention the fixed-star path applies for ECL_T0 — see
    fixed_stars._sidbit_star_call). Equatorial output under
    SIDBIT_SSY_PLANE alone keeps the base behavior (the invariable-plane
    projection leaves the equator unchanged); the callers' guard routes
    only ECL_T0 equatorial here.

    ``tt_echo`` applies the calc()-vs-calc_ut() ephemeris-bit echo rule to
    the returned retflag. The internal sub-request always carries the
    defaulted FLG_SWIEPH, so without it a TT projection echoed that bit even
    when the caller passed none (65602 instead of 65600), breaking the echo
    contract every other calc() path follows.
    """
    from .sidereal_epoch import (
        fixed_epoch_request_flags,
        fixed_epoch_retflag,
        sidbit_ecliptic_matrix,
        ssy_plane_zero_point_deg,
        transform_sidbit_result,
    )

    if flags & FLG_EQUATORIAL:
        # SIDBIT_ECL_T0 equatorial: the star of the request is reduced to the
        # MEAN EQUATOR AND EQUINOX of the mode's t0 — measured reference
        # behavior on both calc and fixstar paths: no ayanamsha subtraction,
        # position and velocity in the t0 mean frame, plain sidereal retflag
        # echo. Re-dispatch as the J2000|NONUT request and precess to t0.
        from .sidereal_epoch import transform_equatorial_epoch_result

        sub_xx, sub_rf = calc_fn(tjd, planet, fixed_epoch_request_flags(flags))
        xx_eq = transform_equatorial_epoch_result(
            sub_xx, flags, _ecl_t0_epoch_jd(sid_mode)
        )
        eq_retflag = _echo_request_bits(fixed_epoch_retflag(sub_rf, flags), raw_flags)
        if tt_echo:
            eq_retflag = _calc_tt_epheflag_echo(eq_retflag, raw_flags)
        return (_to_native_floats(xx_eq), eq_retflag)

    # Measured reference behavior: with BOTH projection bits set,
    # SIDBIT_ECL_T0 takes precedence over SIDBIT_SSY_PLANE.
    if sid_bits & SIDBIT_ECL_T0:
        t0_jd = _ecl_t0_epoch_jd(sid_mode)
        zero_point = _ecl_t0_zero_point_deg(sid_mode, t0_jd)
    else:  # SIDBIT_SSY_PLANE
        t0_jd = _J2000_JD
        zero_point = ssy_plane_zero_point_deg(_calc_ayanamsa(_J2000_JD, sid_mode))

    m_ecl = sidbit_ecliptic_matrix(sid_bits, t0_jd, zero_point)
    if m_ecl is None:
        # The caller's guard should keep us out of here, but an assert
        # disappears under python -O and would hand None to the transform.
        # Fall through to the unprojected sidereal path, as the nod_aps twin
        # already does.
        return calc_fn(tjd, planet, flags)
    sub_xx, sub_rf = calc_fn(tjd, planet, fixed_epoch_request_flags(flags))
    xx = transform_sidbit_result(sub_xx, flags, m_ecl)
    # Echo the caller's representation and SIDEREAL bit; drop the internal
    # J2000 rewrite bit but RETAIN FLG_NONUT. The SIDBIT projection is defined
    # on a MEAN ecliptic and equinox (of the mode's t0 for SIDBIT_ECL_T0, or the
    # invariable plane for SIDBIT_SSY_PLANE), so the projected longitude is
    # declared without nutation -- the same mean-frame convention the base
    # sidereal path echoes via _implied_retflag_bits / fixed_epoch_retflag.
    # sub_rf always carries FLG_NONUT because fixed_epoch_request_flags forces
    # it into the rewritten sub-request.
    # Drop the J2000 bit the internal rewrite forced on, then restore the
    # caller's own representation bits — including FLG_J2000 when it was
    # explicitly requested. Omitting it from this mask made an explicit
    # SIDEREAL|J2000 projection request echo a retflag 32 lower than the
    # measured one, on both calc_ut() and calc().
    retflag = _echo_request_bits(
        (sub_rf & ~FLG_J2000)
        | (flags & (FLG_SIDEREAL | FLG_RADIANS | FLG_XYZ | FLG_J2000)),
        raw_flags,
    )
    if tt_echo:
        retflag = _calc_tt_epheflag_echo(retflag, raw_flags)
    return (_to_native_floats(xx), retflag)


def calc_ut(
    tjdut: float, planet: int, flags: int = FLG_SWIEPH | FLG_SPEED
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """
    Calculate planetary position for Universal Time.

    Reference API compatible function.

    Args:
        tjdut: Julian Day in Universal Time (UT1)
        planet: Planet/body ID (SUN, MOON, etc.)
        flags: Calculation flags (FLG_SPEED, FLG_HELCTR, etc.)

    Returns:
        Tuple containing:
            - Position tuple: (longitude, latitude, distance, speed_lon, speed_lat, speed_dist)
            - Return flag: flags value on success

    Raises:
        EphemerisRangeError: If the date is outside the ephemeris coverage

    Coordinate Output:
        - longitude: Ecliptic longitude in degrees (0-360)
        - latitude: Ecliptic latitude in degrees
        - distance: Distance in AU
        - speed_*: Daily motion in respective coordinates

    Ephemeris Selection Flags:
        The library always uses NASA JPL DE440 (or DE441 via env var) via Skyfield.

        - FLG_SWIEPH / FLG_JPLEPH (default): Uses NASA JPL DE440 ephemeris
          via Skyfield. Valid range depends on loaded ephemeris (DE440: 1550-2650 CE,
          DE441: -13200 to +17191 CE). Highest precision (sub-arcsecond).
          Supports all bodies including asteroids via SPK kernels.

        - FLG_MOSEPH: Accepted for API compatibility but **ignored**.
          Calculations always use JPL ephemeris regardless of this flag.

    Other Flags:
        - FLG_SPEED: Include velocity (default, always calculated)
        - FLG_HELCTR: Heliocentric instead of geocentric
        - FLG_TOPOCTR: Topocentric (requires set_topo)
        - FLG_SIDEREAL: Sidereal zodiac (requires set_sid_mode)

    Example:
        >>> pos, retflag = calc_ut(2451545.0, MARS, FLG_SPEED)
        >>> lon, lat, dist = pos[0], pos[1], pos[2]

        # For extended date range, set LIBEPHEMERIS_EPHEMERIS=de441.bsp
        >>> pos, retflag = calc_ut(1000000.0, MARS, FLG_SPEED)
    """
    from skyfield.errors import EphemerisRangeError as SkyfieldRangeError
    from .exceptions import validate_jd_range
    from .constants import ECL_NUT

    _validate_barycentric_moseph(flags, "calc_ut")
    requested_planet = planet

    # Handle ECL_NUT (-1) - returns nutation and obliquity.
    # ECL_NUT normalizes to one ephemeris bit but retains FLG_MOSEPH and does
    # not remap FLG_SPEED3, so it must run on raw flags before normalization.
    if planet == ECL_NUT:
        # Resolve conflicting center bits before deriving implied bits, and
        # drop SPEED3 when an explicit SPEED bit is also present.
        res_flags = _resolve_center_flags(flags)
        if res_flags & FLG_SPEED:
            res_flags &= ~FLG_SPEED3
        pos_nut, rf_nut = _calc_nutation_obliquity(
            tjdut, _exclusive_ephemeris_bit(res_flags)
        )
        result = (pos_nut, rf_nut | _implied_retflag_bits(res_flags))
        from .logging_config import get_logger

        get_logger().debug("body=%d jd=%.1f source=ERFA", requested_planet, tjdut)
        _record(requested_planet, "ERFA")
        return result

    flags = consume_non_echoed_flags(flags, planet)

    raw_flags = flags
    flags = _normalize_calc_flags(flags)

    # --- Fixed-epoch sidereal modes (SIDM_J2000/J1900/B1950) ---
    # These modes are frame requests, not ayanamsha offsets: the reference
    # returns the FLG_J2000|FLG_NONUT position (precessed to t0 for
    # J1900/B1950) for every representation. Rewrite the request and map
    # the result back (see sidereal_epoch.py).
    if flags & FLG_SIDEREAL:
        from .sidereal_epoch import (
            fixed_epoch_request_flags,
            fixed_epoch_retflag,
            is_fixed_epoch_request,
            transform_fixed_epoch_result,
        )

        _sidm = get_sid_mode()
        if is_fixed_epoch_request(flags, _sidm):
            nested_trace = (
                _snapshot_record(requested_planet) if _is_active() else (False, None)
            )
            try:
                with _capture_source_logs() as nested_sources:
                    sub_xx, sub_rf = calc_ut(
                        tjdut, planet, fixed_epoch_request_flags(flags)
                    )
                    xx_t0 = transform_fixed_epoch_result(sub_xx, flags, _sidm)
                    result = (
                        _to_native_floats(xx_t0),
                        _echo_request_bits(
                            fixed_epoch_retflag(sub_rf, flags), raw_flags
                        ),
                    )
            except Exception:
                _restore_record(requested_planet, nested_trace)
                raise
            if nested_sources:
                from .logging_config import get_logger

                _log_successful_source(
                    get_logger(), requested_planet, tjdut, nested_sources[-1]
                )
            return result

        # SIDBIT_ECL_T0 / SIDBIT_SSY_PLANE frame projections. Ecliptic output
        # is projected onto the plane; equatorial output is routed only for
        # SIDBIT_ECL_T0 (measured: the J2000|NONUT mean frame), while
        # SSY_PLANE alone leaves the equator unchanged. The star/galactic
        # "true" modes define their zero on the ecliptic of date, so the
        # reference leaves the projection inert for them (see
        # _SIDBIT_PROJECTION_SUPPRESS_MODES): fall through to the base
        # sidereal path, which reproduces the un-projected baseline.
        _bits = _get_sidereal_bits()
        if (
            (_bits & (SIDBIT_ECL_T0 | SIDBIT_SSY_PLANE))
            and (not (flags & FLG_EQUATORIAL) or (_bits & SIDBIT_ECL_T0))
            and _sidm not in _SIDBIT_PROJECTION_SUPPRESS_MODES
        ):
            return _sidbit_projection_calc(
                tjdut, planet, flags, raw_flags, _sidm, _bits, calc_ut
            )

    # Built-in asteroids by AST_OFFSET number: remap before LEB/Horizons
    # dispatch so both id forms are served by the same backend.
    planet = _remap_ast_offset(planet)

    # --- South nodes: derive from north node via the same dispatch path ---
    # South node = north node + 180° longitude, negated latitude.
    # Must be handled here (before LEB/Horizons) so the north node sub-call
    # goes through whichever backend (LEB, Horizons, Skyfield) is active.
    # Otherwise, negative body IDs fall through LEB/Horizons (which don't
    # store them) to the Skyfield path, causing velocity mismatches when
    # LEB computes the north node with Chebyshev derivatives but Skyfield
    # computes the south node with numerical differentiation.
    from .constants import MEAN_NODE, TRUE_NODE

    if planet in (-MEAN_NODE, -TRUE_NODE):
        north_ipl = abs(planet)
        north_trace = _snapshot_record(north_ipl) if _is_active() else (False, None)
        try:
            with _capture_source_logs() as nested_sources:
                north_result, retflag = calc_ut(tjdut, north_ipl, flags)
                result = (
                    _south_node_from_north(north_result, flags),
                    _echo_request_bits(retflag, raw_flags),
                )
        except Exception:
            _restore_record(north_ipl, north_trace)
            raise
        _record_alias(requested_planet, north_ipl, north_trace)
        if nested_sources:
            from .logging_config import get_logger

            _log_successful_source(
                get_logger(), requested_planet, tjdut, nested_sources[-1]
            )
        return result

    # --- LEB fast path: use precomputed binary ephemeris if available ---
    from .state import get_leb_reader
    from .logging_config import get_logger

    reader = get_leb_reader()
    # Degenerate topocentric Earth: the reference returns the exact zero
    # vector (Earth is the coordinate origin; TOPOCTR is echoed but the
    # position stays zero). The Skyfield path already handles this; the LEB
    # fast path would instead return the geocentre-from-topocentre offset
    # (~1 Earth radius), so route the case past LEB.
    if planet == EARTH and (flags & FLG_TOPOCTR):
        reader = None
    if reader is not None:
        try:
            from . import fast_calc

            result = fast_calc.fast_calc_ut(reader, tjdut, planet, flags)
            pos_out = _to_native_floats(result[0])
            retflag_out = _echo_request_bits(
                result[1] | _implied_retflag_bits(flags), raw_flags
            )
            _log_successful_source(get_logger(), requested_planet, tjdut, "LEB")
            _record(requested_planet, "LEB")
            return pos_out, retflag_out
        except LEBCorruptionError:
            raise
        except (KeyError, ValueError) as _leb_err:
            _raise_leb_range_miss(planet, tjdut)
            from .leb_reader import log_leb_fallback

            log_leb_fallback(f"body={planet} jd={tjdut:.1f}", _leb_err)
    # --- END LEB fast path ---

    # --- Horizons path: use NASA JPL Horizons API when no local ephemeris ---
    from .state import get_horizons_client

    h_client = get_horizons_client()
    if h_client is not None:
        try:
            from . import horizons_backend

            result = horizons_backend.horizons_calc_ut(h_client, tjdut, planet, flags)
            pos_out = _to_native_floats(result[0])
            retflag_out = _echo_request_bits(
                result[1] | _implied_retflag_bits(flags), raw_flags
            )
            trace_source = horizons_backend._trace_source(planet, flags)
            _log_successful_source(get_logger(), requested_planet, tjdut, trace_source)
            _record(requested_planet, trace_source)
            return pos_out, retflag_out
        except KeyError as _hz_err:
            get_logger().debug(
                "body=%d jd=%.1f source=Horizons->fallback reason=%s",
                planet,
                tjdut,
                _hz_err,
            )
        except (ConnectionError, ValueError, OSError) as _hz_err:
            # Transient network/API failure, or a ValueError/OSError raised
            # from the HELCTR / velocity / analytical path — fall through to
            # the Skyfield path (same fallback the context API uses) rather
            # than let it escape calc_ut.
            get_logger().warning(
                "body=%d jd=%.1f source=Horizons->fallback (network): %s",
                planet,
                tjdut,
                _hz_err,
            )
    # --- END Horizons path ---

    # Validate JD range for bodies that use the JPL ephemeris
    if _body_uses_jpl_ephemeris(planet):
        validate_jd_range(tjdut, planet, "calc_ut")

    # Strip FLG_XYZ and FLG_RADIANS from the flags passed to _calc_body
    # since they are output format flags, not calculation flags.
    # We apply them after the calculation is complete.
    calc_iflag = _strip_output_flags(flags)
    trace_active = _is_active()
    requested_trace = (
        _snapshot_record(requested_planet) if trace_active else (False, None)
    )
    alias_trace = _snapshot_record(planet) if trace_active else (False, None)

    def _restore_failed_trace() -> None:
        _restore_record(planet, alias_trace)
        _restore_record(requested_planet, requested_trace)

    from .cache import get_cached_time_ut1

    t = get_cached_time_ut1(tjdut)
    source_token = _dispatch_source.set(None)
    try:
        # A concurrent close() (module state or an EphemerisContext) can empty
        # the SPK kernel or LEB reader mid-read: ValueError ("seek of closed
        # file") or a KeyError from a transiently emptied kernel/reader.
        # _calc_body_race_safe retries after the lazy reload; anything else
        # re-raises unchanged (see _is_transient_close_race).
        pos, retflag = _calc_body_race_safe(t, planet, calc_iflag)
        source_hint = _dispatch_source.get()
        pos_out, rf_out = _finalize_output_flags(pos, retflag, flags)
        pos_out = _rebuild_mean_point_speed3(
            pos_out, t, planet, flags, _calc_body_race_safe
        )
        retflag_out = _echo_request_bits(rf_out, raw_flags)
        logger = get_logger()
        if trace_active or logger.isEnabledFor(10) or source_hint == "Keplerian":
            fallback_source = _fallback_trace_source(planet, source_hint)
            if fallback_source is None:
                _record_alias(requested_planet, planet, alias_trace)
            else:
                _log_successful_source(logger, requested_planet, tjdut, fallback_source)
                _record(requested_planet, fallback_source)
        return pos_out, retflag_out
    except SkyfieldRangeError as e:
        _restore_failed_trace()
        raise _wrap_ephemeris_range_error(e, tjdut, planet) from e
    except ValueError as e:
        _restore_failed_trace()
        # SPK type 21 kernels (asteroids) raise ValueError when date is
        # outside the kernel's coverage. Convert to EphemerisRangeError.
        if "Invalid Time" in str(e) or "time" in str(e).lower():
            raise _wrap_ephemeris_range_error(e, tjdut, planet) from e
        raise
    except Exception:
        _restore_failed_trace()
        raise
    finally:
        _dispatch_source.reset(source_token)


def calc(
    tjdet: float, planet: int, flags: int = FLG_SWIEPH | FLG_SPEED
) -> tuple[tuple[float, float, float, float, float, float], int]:
    """
    Calculate planetary position for Ephemeris Time (ET/TT).

    Reference API compatible function. Similar to calc_ut() but takes
    Terrestrial Time (TT, also known as Ephemeris Time) instead of Universal Time.

    Args:
        tjdet: Julian Day in Terrestrial Time (TT/ET)
        planet: Planet/body ID (SUN, MOON, etc.)
        flags: Calculation flags (FLG_SPEED, FLG_HELCTR, etc.)

    Returns:
        Tuple containing:
            - Position tuple: (longitude, latitude, distance, speed_lon, speed_lat, speed_dist)
            - Return flag: flags value on success

    Raises:
        EphemerisRangeError: If the date is outside the ephemeris coverage

    Note:
        TT (Terrestrial Time) differs from UT (Universal Time) by Delta T,
        which varies from ~32 seconds (year 2000) to minutes (historical times).
        For most astrological applications, use calc_ut() instead.

    Example:
        >>> pos, retflag = calc(2451545.0, JUPITER, FLG_SPEED)
        >>> lon, lat, dist = pos[0], pos[1], pos[2]
    """
    from skyfield.errors import EphemerisRangeError as SkyfieldRangeError
    from .exceptions import validate_jd_range
    from .constants import ECL_NUT

    _validate_barycentric_moseph(flags, "calc")
    requested_planet = planet

    # Handle ECL_NUT (-1) — nutation and obliquity. The input is already
    # TT, so compute directly (calc_ut converts UT first; this mirror was
    # missing here and ECL_NUT fell through to UnknownBodyError). Apply
    # exclusive-ephemeris-bit handling to the raw flags, as calc_ut does.
    if planet == ECL_NUT:
        # Use the same normalization as calc_ut(): ECL_NUT is an ERFA-backed
        # pseudo-body, so TT versus UT changes only the time conversion, not
        # the public flag contract. In particular, a request with no
        # ephemeris-selection bit receives the normal FLG_SWIEPH default.
        res_flags = _resolve_center_flags(flags)
        if res_flags & FLG_SPEED:
            res_flags &= ~FLG_SPEED3
        nut_flags = _exclusive_ephemeris_bit(res_flags)
        pos_nut, rf_nut = _calc_nutation_obliquity_tt(tjdet, nut_flags)
        # calc() (TT) echoes only the ephemeris-selection bits the caller
        # actually passed: a request with no source bit returns retflag 0,
        # not the SWIEPH default that calc_ut() injects. Measured against the
        # reference: calc ECL_NUT 0->0, 64->64, 256->256, 320->320, while
        # calc_ut ECL_NUT 0->2, 64->66, ... (SWIEPH default). Requests that
        # carry an explicit source bit (2, 4, 6, 66, 258, ...) are unaffected.
        retflag = _calc_tt_epheflag_echo(
            rf_nut | _implied_retflag_bits(res_flags), flags
        )
        result = (pos_nut, retflag)
        from .logging_config import get_logger

        get_logger().debug("body=%d jd=%.1f source=ERFA", requested_planet, tjdet)
        _record(requested_planet, "ERFA")
        return result

    flags = consume_non_echoed_flags(flags, planet)

    raw_flags = flags
    flags = _normalize_calc_flags(flags)

    # --- Fixed-epoch sidereal modes: frame requests, see calc_ut ---
    if flags & FLG_SIDEREAL:
        from .sidereal_epoch import (
            fixed_epoch_request_flags,
            fixed_epoch_retflag,
            is_fixed_epoch_request,
            transform_fixed_epoch_result,
        )

        _sidm = get_sid_mode()
        if is_fixed_epoch_request(flags, _sidm):
            nested_trace = (
                _snapshot_record(requested_planet) if _is_active() else (False, None)
            )
            try:
                with _capture_source_logs() as nested_sources:
                    sub_xx, sub_rf = calc(
                        tjdet, planet, fixed_epoch_request_flags(flags)
                    )
                    xx_t0 = transform_fixed_epoch_result(sub_xx, flags, _sidm)
                    result = (
                        _to_native_floats(xx_t0),
                        _calc_tt_epheflag_echo(
                            _echo_request_bits(
                                fixed_epoch_retflag(sub_rf, flags), raw_flags
                            ),
                            raw_flags,
                        ),
                    )
            except Exception:
                _restore_record(requested_planet, nested_trace)
                raise
            if nested_sources:
                from .logging_config import get_logger

                _log_successful_source(
                    get_logger(), requested_planet, tjdet, nested_sources[-1]
                )
            return result

        # SIDBIT_ECL_T0 / SIDBIT_SSY_PLANE frame projections. Ecliptic output
        # is projected onto the plane; equatorial output is routed only for
        # SIDBIT_ECL_T0 (measured: the J2000|NONUT mean frame), while
        # SSY_PLANE alone leaves the equator unchanged. The star/galactic
        # "true" modes define their zero on the ecliptic of date, so the
        # reference leaves the projection inert for them (see
        # _SIDBIT_PROJECTION_SUPPRESS_MODES): fall through to the base
        # sidereal path, which reproduces the un-projected baseline.
        _bits = _get_sidereal_bits()
        if (
            (_bits & (SIDBIT_ECL_T0 | SIDBIT_SSY_PLANE))
            and (not (flags & FLG_EQUATORIAL) or (_bits & SIDBIT_ECL_T0))
            and _sidm not in _SIDBIT_PROJECTION_SUPPRESS_MODES
        ):
            return _sidbit_projection_calc(
                tjdet, planet, flags, raw_flags, _sidm, _bits, calc, tt_echo=True
            )

    # Built-in asteroids by AST_OFFSET number (see calc_ut)
    planet = _remap_ast_offset(planet)

    # --- South nodes: derive from the north node via the same dispatch
    # path (see calc_ut for the backend-consistency rationale) ---
    from .constants import MEAN_NODE, TRUE_NODE

    if planet in (-MEAN_NODE, -TRUE_NODE):
        north_ipl = abs(planet)
        north_trace = _snapshot_record(north_ipl) if _is_active() else (False, None)
        try:
            with _capture_source_logs() as nested_sources:
                north_result, retflag = calc(tjdet, north_ipl, flags)
                result = (
                    _south_node_from_north(north_result, flags),
                    _calc_tt_epheflag_echo(
                        _echo_request_bits(retflag, raw_flags), raw_flags
                    ),
                )
        except Exception:
            _restore_record(north_ipl, north_trace)
            raise
        _record_alias(requested_planet, north_ipl, north_trace)
        if nested_sources:
            from .logging_config import get_logger

            _log_successful_source(
                get_logger(), requested_planet, tjdet, nested_sources[-1]
            )
        return result

    # --- LEB fast path: use precomputed binary ephemeris if available ---
    from .state import get_leb_reader
    from .logging_config import get_logger

    reader = get_leb_reader()
    # Degenerate topocentric Earth: see calc_ut — the reference returns the
    # zero vector; keep this off the LEB fast path.
    if planet == EARTH and (flags & FLG_TOPOCTR):
        reader = None
    if reader is not None:
        try:
            from . import fast_calc

            result = fast_calc.fast_calc_tt(reader, tjdet, planet, flags)
            pos_out = _to_native_floats(result[0])
            retflag_out = _calc_tt_epheflag_echo(
                _echo_request_bits(result[1] | _implied_retflag_bits(flags), raw_flags),
                raw_flags,
            )
            _log_successful_source(get_logger(), requested_planet, tjdet, "LEB")
            _record(requested_planet, "LEB")
            return pos_out, retflag_out
        except LEBCorruptionError:
            raise
        except (KeyError, ValueError) as _leb_err:
            _raise_leb_range_miss(planet, tjdet)
            from .leb_reader import log_leb_fallback

            log_leb_fallback(f"body={planet} jd={tjdet:.1f}", _leb_err)
    # --- END LEB fast path ---

    # --- Horizons path ---
    from .state import get_horizons_client

    h_client = get_horizons_client()
    if h_client is not None:
        try:
            from . import horizons_backend
            from .time_utils import deltat

            # calc uses TT, convert to UT for horizons_calc_ut
            jd_ut_approx = tjdet - deltat(tjdet)
            result = horizons_backend.horizons_calc_ut(
                h_client, jd_ut_approx, planet, flags
            )
            pos_out = _to_native_floats(result[0])
            retflag_out = _calc_tt_epheflag_echo(
                _echo_request_bits(result[1] | _implied_retflag_bits(flags), raw_flags),
                raw_flags,
            )
            trace_source = horizons_backend._trace_source(planet, flags)
            _log_successful_source(get_logger(), requested_planet, tjdet, trace_source)
            _record(requested_planet, trace_source)
            return pos_out, retflag_out
        except KeyError as _hz_err:
            get_logger().debug(
                "body=%d jd=%.1f source=Horizons->fallback reason=%s",
                planet,
                tjdet,
                _hz_err,
            )
        except (ConnectionError, ValueError, OSError) as _hz_err:
            # Transient network/API failure, or a ValueError/OSError raised
            # from the HELCTR / velocity / analytical path — fall through to
            # the Skyfield path (same fallback the context API uses) rather
            # than let it escape calc.
            get_logger().warning(
                "body=%d jd=%.1f source=Horizons->fallback (network): %s",
                planet,
                tjdet,
                _hz_err,
            )
    # --- END Horizons path ---

    # Validate JD range for bodies that use the JPL ephemeris
    if _body_uses_jpl_ephemeris(planet):
        validate_jd_range(tjdet, planet, "calc")

    # Strip FLG_XYZ and FLG_RADIANS from the flags passed to _calc_body
    # since they are output format flags, not calculation flags.
    calc_iflag = _strip_output_flags(flags)
    trace_active = _is_active()
    requested_trace = (
        _snapshot_record(requested_planet) if trace_active else (False, None)
    )
    alias_trace = _snapshot_record(planet) if trace_active else (False, None)

    def _restore_failed_trace() -> None:
        _restore_record(planet, alias_trace)
        _restore_record(requested_planet, requested_trace)

    from .cache import get_cached_time_tt

    t = get_cached_time_tt(tjdet)
    source_token = _dispatch_source.set(None)
    try:
        # Concurrent close() emptied the SPK kernel or LEB reader mid-read; see
        # calc_ut and _is_transient_close_race. Retry after the lazy reload.
        pos, retflag = _calc_body_race_safe(t, planet, calc_iflag)
        source_hint = _dispatch_source.get()
        pos_out, rf_out = _finalize_output_flags(pos, retflag, flags)
        pos_out = _rebuild_mean_point_speed3(
            pos_out, t, planet, flags, _calc_body_race_safe
        )
        retflag_out = _calc_tt_epheflag_echo(
            _echo_request_bits(rf_out, raw_flags), raw_flags
        )
        logger = get_logger()
        if trace_active or logger.isEnabledFor(10) or source_hint == "Keplerian":
            fallback_source = _fallback_trace_source(planet, source_hint)
            if fallback_source is None:
                _record_alias(requested_planet, planet, alias_trace)
            else:
                _log_successful_source(logger, requested_planet, tjdet, fallback_source)
                _record(requested_planet, fallback_source)
        return pos_out, retflag_out
    except SkyfieldRangeError as e:
        _restore_failed_trace()
        raise _wrap_ephemeris_range_error(e, tjdet, planet) from e
    except ValueError as e:
        _restore_failed_trace()
        if "Invalid Time" in str(e) or "time" in str(e).lower():
            raise _wrap_ephemeris_range_error(e, tjdet, planet) from e
        raise
    except Exception:
        _restore_failed_trace()
        raise
    finally:
        _dispatch_source.reset(source_token)


def calc_pctr(
    tjdet: float, planet: int, center: int, flags: int = FLG_SWIEPH | FLG_SPEED
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """
    Calculate planetary position as seen from another planet (planet-centric).

    Reference API compatible function.

    This function calculates the position of a target body (planet) as observed
    from another body (center) rather than from Earth (geocentric) or Sun
    (heliocentric). Useful for calculating, e.g., the position of Moon as
    seen from Mars, or Venus as seen from Jupiter.

    Args:
        tjdet: Julian Day in Terrestrial Time (TT/ET)
        planet: Target planet/body ID (SUN, MOON, etc.)
        center: Observer/center planet ID (the body from which to observe)
        flags: Calculation flags (FLG_SPEED, etc.)

    Returns:
        Tuple containing:
            - Position tuple: (longitude, latitude, distance, speed_lon, speed_lat, speed_dist)
            - Return flag: flags value on success

    Raises:
        EphemerisRangeError: If the date is outside the ephemeris coverage

    Note:
        - FLG_HELCTR / FLG_BARYCTR do not change the observer (always the
          center body) but, like the reference API, they imply an
          astrometric place: stellar aberration and gravitational light
          deflection are skipped (retflag echoes NOABERR|NOGDEFL and the
          center bit itself is stripped, matching the reference)
        - FLG_TOPOCTR is ignored (no topocentric correction on other planets)
        - Distance is the distance from center to planet in AU

    Example:
        >>> # Position of Moon as seen from Mars
        >>> pos, retflag = calc_pctr(2451545.0, MOON, MARS, FLG_SPEED)
        >>> print(f"Moon longitude from Mars: {pos[0]:.2f}°")
    """
    from skyfield.errors import EphemerisRangeError as SkyfieldRangeError
    from .exceptions import validate_jd_range

    # HELCTR/BARYCTR on calc_pctr use the public astrometric-place convention:
    # place (no aberration, no deflection — light-time still applied) and
    # strips the center bit from the echoed retflag. Map the bits onto the
    # explicit NOABERR|NOGDEFL flags so the computation gates and the
    # retflag echo both follow the explicit NOABERR|NOGDEFL path.
    if flags & (FLG_HELCTR | FLG_BARYCTR):
        flags = (flags | FLG_NOABERR | FLG_NOGDEFL) & ~(FLG_HELCTR | FLG_BARYCTR)

    # Sidereal FRAME requests are whole-vector rotations, not scalar ayanamsha
    # subtractions, so they must be intercepted here exactly as calc()/calc_ut()
    # do. Without this the planet-centric surface silently returned the base
    # sidereal position: measured against the reference, Mars seen from Jupiter
    # at J2000 kept latitude 0.978545 deg where SIDBIT_SSY_PLANE gives -0.375689
    # and SIDM_J1900 gives 0.989033.
    if flags & FLG_SIDEREAL:
        from .sidereal_epoch import (
            fixed_epoch_request_flags,
            fixed_epoch_retflag,
            is_fixed_epoch_request,
            transform_fixed_epoch_result,
        )

        _pctr_raw_flags = flags
        _sidm_p = get_sid_mode()

        def _pctr_calc_fn(_tjd: float, _ipl: int, _flags: int):
            return calc_pctr(_tjd, _ipl, center, _flags)

        if is_fixed_epoch_request(flags, _sidm_p):
            sub_xx, sub_rf = _pctr_calc_fn(
                tjdet, planet, fixed_epoch_request_flags(flags)
            )
            xx_t0 = transform_fixed_epoch_result(sub_xx, flags, _sidm_p)
            return (
                _to_native_floats(xx_t0),
                _echo_request_bits(fixed_epoch_retflag(sub_rf, flags), _pctr_raw_flags),
            )

        _bits_p = _get_sidereal_bits()
        if (
            (_bits_p & (SIDBIT_ECL_T0 | SIDBIT_SSY_PLANE))
            and (not (flags & FLG_EQUATORIAL) or (_bits_p & SIDBIT_ECL_T0))
            and _sidm_p not in _SIDBIT_PROJECTION_SUPPRESS_MODES
        ):
            return _sidbit_projection_calc(
                tjdet,
                planet,
                flags,
                _pctr_raw_flags,
                _sidm_p,
                _bits_p,
                _pctr_calc_fn,
            )

    # Validate JD range for bodies that use the JPL ephemeris
    if _body_uses_jpl_ephemeris(planet) or _body_uses_jpl_ephemeris(center):
        validate_jd_range(tjdet, planet, "calc_pctr")

    from .cache import get_cached_time_tt

    t = get_cached_time_tt(tjdet)
    try:
        result = _run_pctr_pipeline(
            lambda calc_iflag: _calc_body_pctr(t, planet, center, calc_iflag),
            flags,
        )
        from .state import get_calc_mode
        from .logging_config import get_logger

        source = "LEB" if get_calc_mode() == "leb" else "Skyfield"
        _log_successful_source(get_logger(), planet, tjdet, source)
        _record(planet, source)
        return result
    except SkyfieldRangeError as e:
        raise _wrap_ephemeris_range_error(e, tjdet, planet) from e


def _calc_body_pctr(
    t, ipl: int, iplctr: int, iflag: int
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """
    Calculate planet-centric position (internal function).

    Calculates the position of target body (ipl) as seen from center body (iplctr).
    Uses vector subtraction: position = target_SSB - observer_SSB.

    Args:
        t: Skyfield Time object (UT1 or TT)
        ipl: Target planet/body ID
        iplctr: Observer/center planet ID
        iflag: Calculation flags (FLG_SPEED, etc.)

    Returns:
        Tuple of (position_tuple, flags)
    """
    from skyfield.positionlib import ICRF

    from .state import _get_computation_ephemeris

    planets = _get_computation_ephemeris()

    # Validate that both bodies are in _PLANET_MAP (standard planets only for now)
    if ipl not in _PLANET_MAP:
        # Target not supported - raise clear error
        from .exceptions import UnknownBodyError

        raise UnknownBodyError(
            message=(
                f"Unknown target body ID {ipl} for planet-centric calculation. "
                f"calc_pctr() only supports standard planets (Sun=0 through Earth=14). "
                f"See libephemeris.constants for body ID constants."
            ),
            body_id=ipl,
        )

    if iplctr not in _PLANET_MAP:
        # Observer not supported - raise clear error
        from .exceptions import UnknownBodyError

        raise UnknownBodyError(
            message=(
                f"Unknown observer/center body ID {iplctr} for planet-centric calculation. "
                f"calc_pctr() only supports standard planets (Sun=0 through Earth=14) "
                f"as the observer. See libephemeris.constants for body ID constants."
            ),
            body_id=iplctr,
        )

    target_name = _PLANET_MAP[ipl]
    observer_name = _PLANET_MAP[iplctr]

    # Use get_planet_target() to get planet center (with COB correction) for gas giants
    # This ensures we use planet center NAIF IDs (599, 699, 799, 899) rather than
    # barycenter IDs (5, 6, 7, 8), providing sub-arcsecond positional accuracy
    target = get_planet_target(planets, target_name)
    observer = get_planet_target(planets, observer_name)

    # Use a fresh Time object to avoid Skyfield reify descriptor corruption
    # when the same Time is shared across multiple callers. Preserve the
    # (whole, fraction) split: collapsing to a single float64 (tt_jd(float(t.tt)))
    # quantizes the epoch to the JD ULP (~4.6e-10 d), which destroys the exact
    # ±dt spacing of the FLG_SPEED central-difference samples and biased the
    # Moon's planet-centric speed by ~-0.4"/day.
    t_fresh = get_timescale().tt_jd(float(t.whole), float(t.tt_fraction))

    # Helper function to get position vector at time t_
    # NOTE: We do NOT use get_cached_observer_at here because the cache
    # can return positions computed with a corrupted Time object.
    def get_vector(t_):
        # Both positions relative to SSB (Solar System Barycenter)
        tgt = target.at(t_)
        tgt_pos = tgt.position.au
        tgt_vel = tgt.velocity.au_per_d
        obs = observer.at(t_)
        obs_pos = obs.position.au
        obs_vel = obs.velocity.au_per_d

        # Target position relative to observer
        p_ = tgt_pos - obs_pos
        v_ = tgt_vel - obs_vel
        return p_, v_

    # Get position vector with light-time correction
    # Iterative light-time: target position at t - light_travel_time
    import numpy as np

    # Speed of light in AU/day: NOVAS/Skyfield C_AUDAY realization
    # (c exact from SI, IAU 1976 au), shared with astrometry.py and fast_calc.py.
    C_AU_PER_DAY = 173.1446326847

    if iflag & FLG_TRUEPOS:
        # Geometric position (no light-time correction)
        p, v = get_vector(t_fresh)
    else:
        # Light-time corrected position
        p, v = get_vector(t_fresh)
        # Hoist observer.at(t) out of the loop - observer stays at current time
        obs_at_t = observer.at(t_fresh)
        obs_pos = obs_at_t.position.au
        obs_vel = obs_at_t.velocity.au_per_d
        for _ in range(3):
            dist_au = np.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2)
            light_time = dist_au / C_AU_PER_DAY
            ts_lt = get_timescale()
            # Retard only the target, keep observer at current time. Keep the
            # (whole, fraction) split here too: collapsing tdb - light_time to
            # one float64 re-quantizes the retarded epoch (ULP noise that the
            # speed central difference amplifies by 1/(2*dt)).
            tgt_ret = target.at(
                ts_lt.tdb_jd(
                    float(t_fresh.whole),
                    float(t_fresh.tdb_fraction) - light_time,
                )
            )
            p = tgt_ret.position.au - obs_pos
            v = tgt_ret.velocity.au_per_d - obs_vel

        # Gravitational light deflection by Sun/Jupiter/Saturn, applied before
        # aberration (Skyfield .apparent() order) so calc_pctr and calc share
        # the same IAU/Skyfield apparent-place reduction. Skipped under
        # FLG_NOGDEFL.
        if not (iflag & FLG_NOGDEFL):
            from .fast_calc import _apply_gravitational_deflection

            p = np.array(
                _apply_gravitational_deflection(
                    (float(p[0]), float(p[1]), float(p[2])),
                    (float(obs_pos[0]), float(obs_pos[1]), float(obs_pos[2])),
                    float(t_fresh.tt),
                    light_time,
                    _SkyfieldDeflectorSource(planets),
                )
            )

        # Stellar aberration about the observer's barycentric velocity, using
        # the same relativistic correction Skyfield's .apparent() applies in
        # calc(). FLG_NOABERR skips it; FLG_TRUEPOS already skips light time
        # and this entire block.
        if not (iflag & FLG_NOABERR):
            from .fast_calc import _apply_aberration

            p = np.array(
                _apply_aberration(
                    (float(p[0]), float(p[1]), float(p[2])),
                    (float(obs_vel[0]), float(obs_vel[1]), float(obs_vel[2])),
                    light_time,
                )
            )

    # Create position object for coordinate conversion
    pos = ICRF(p, v, t=t_fresh, center=399)

    # Extract coordinates based on flags
    is_equatorial = bool(iflag & FLG_EQUATORIAL)
    is_icrs = bool(iflag & FLG_ICRS)
    is_sidereal = bool(iflag & FLG_SIDEREAL)
    # Mean equator of date when nutation is suppressed (FLG_NONUT) or the
    # sidereal+equatorial frame is requested -- mirrors _calc_body.
    _use_mean_equator = bool(iflag & FLG_NONUT) or is_sidereal

    p1, p2, p3 = 0.0, 0.0, 0.0
    dp1, dp2, dp3 = 0.0, 0.0, 0.0

    if is_equatorial:
        # Equatorial coordinates (RA/Dec)
        if iflag & FLG_J2000:
            # Mean equator/equinox of J2000: ICRS rotated by IAU 2006 frame
            # bias.
            from .fast_calc import _rotate_icrs_to_j2000_mean_equator

            x, y, z = pos.position.au
            xe, ye, ze = _rotate_icrs_to_j2000_mean_equator(
                float(x), float(y), float(z)
            )
            dist_val = math.sqrt(xe * xe + ye * ye + ze * ze)
            p1 = math.degrees(math.atan2(ye, xe)) % 360.0
            p2 = (
                math.degrees(math.asin(max(-1.0, min(1.0, ze / dist_val))))
                if dist_val > 0
                else 0.0
            )
            p3 = dist_val
        elif is_icrs:
            # ICRS equator of date (no frame bias): Vondrák 2011 precession,
            # with Skyfield nutation unless the mean equator is requested.
            import numpy as np

            if _use_mean_equator:
                pn = vondrak_precession_matrix(t.tt, frame_bias=False)
            else:
                dpsi, deps = t._nutation_angles_radians
                pn, _eps_true = vondrak_pn_matrix(
                    t.tt, float(dpsi), float(deps), frame_bias=False
                )
            xyz_eq = np.array(pn) @ np.array(pos.position.au)
            xe, ye, ze = float(xyz_eq[0]), float(xyz_eq[1]), float(xyz_eq[2])
            dist_val = math.sqrt(xe * xe + ye * ye + ze * ze)
            p1 = math.degrees(math.atan2(ye, xe)) % 360.0
            p2 = (
                math.degrees(math.asin(max(-1.0, min(1.0, ze / dist_val))))
                if dist_val > 0
                else 0.0
            )
            p3 = dist_val
        else:
            # Equator of date (with frame bias): Vondrák 2011 precession, with
            # Skyfield nutation unless the mean equator is requested. Matches the
            # LEB fast path.
            import numpy as np

            if _use_mean_equator:
                pn = vondrak_precession_matrix(t.tt)
            else:
                dpsi, deps = t._nutation_angles_radians
                pn, _eps_true = vondrak_pn_matrix(t.tt, float(dpsi), float(deps))
            xyz_eq = np.array(pn) @ np.array(pos.position.au)
            xe, ye, ze = float(xyz_eq[0]), float(xyz_eq[1]), float(xyz_eq[2])
            dist_val = math.sqrt(xe * xe + ye * ye + ze * ze)
            p1 = math.degrees(math.atan2(ye, xe)) % 360.0
            p2 = (
                math.degrees(math.asin(max(-1.0, min(1.0, ze / dist_val))))
                if dist_val > 0
                else 0.0
            )
            p3 = dist_val
    else:
        # Ecliptic coordinates (default)
        if iflag & FLG_J2000:
            # Reference J2000 ecliptic frame (frame bias + IAU 2006 obliquity)
            from .fast_calc import _rotate_icrs_to_ecliptic_j2000

            x, y, z = pos.position.au
            xe, ye, ze = _rotate_icrs_to_ecliptic_j2000(float(x), float(y), float(z))

            dist = math.sqrt(xe * xe + ye * ye + ze * ze)
            lon = math.degrees(math.atan2(ye, xe)) % 360.0
            lat = (
                math.degrees(math.asin(max(-1.0, min(1.0, ze / dist))))
                if dist > 0
                else 0.0
            )

            p1, p2, p3 = lon, lat, dist
        elif is_icrs:
            # ICRS ecliptic of date (no frame bias): Vondrák 2011 precession,
            # with Skyfield nutation unless FLG_NONUT requests the mean ecliptic,
            # then rotate equator -> ecliptic by eps_true.
            import numpy as np

            if iflag & FLG_NONUT:
                pn = vondrak_precession_matrix(t.tt, frame_bias=False)
                eps_true = vondrak_mean_obliquity_rad(t.tt)
            else:
                dpsi, deps = t._nutation_angles_radians
                pn, eps_true = vondrak_pn_matrix(
                    t.tt, float(dpsi), float(deps), frame_bias=False
                )
            xq, yq, zq = np.array(pn) @ np.array(pos.position.au)
            ce, se = math.cos(eps_true), math.sin(eps_true)
            xe = float(xq)
            ye = float(yq) * ce + float(zq) * se
            ze = -float(yq) * se + float(zq) * ce
            dist = math.sqrt(xe * xe + ye * ye + ze * ze)
            p1 = math.degrees(math.atan2(ye, xe)) % 360.0
            p2 = (
                math.degrees(math.asin(max(-1.0, min(1.0, ze / dist))))
                if dist > 0
                else 0.0
            )
            p3 = dist
        else:
            # Ecliptic of date: Vondrák 2011 precession, with Skyfield nutation
            # unless FLG_NONUT requests the mean ecliptic, then rotate equator ->
            # ecliptic by eps_true. Matches the LEB fast path.
            import numpy as np

            if iflag & FLG_NONUT:
                pn = vondrak_precession_matrix(t.tt)
                eps_true = vondrak_mean_obliquity_rad(t.tt)
            else:
                dpsi, deps = t._nutation_angles_radians
                pn, eps_true = vondrak_pn_matrix(t.tt, float(dpsi), float(deps))
            xq, yq, zq = np.array(pn) @ np.array(pos.position.au)
            ce, se = math.cos(eps_true), math.sin(eps_true)
            xe = float(xq)
            ye = float(yq) * ce + float(zq) * se
            ze = -float(yq) * se + float(zq) * ce
            dist = math.sqrt(xe * xe + ye * ye + ze * ze)
            p1 = math.degrees(math.atan2(ye, xe)) % 360.0
            p2 = (
                math.degrees(math.asin(max(-1.0, min(1.0, ze / dist))))
                if dist > 0
                else 0.0
            )
            p3 = dist

    # Calculate speed using central difference numerical differentiation if requested
    # Central difference: f'(x) ≈ [f(x+h) - f(x-h)] / (2h) - error O(h²)
    #
    # Use the same SI-second half-steps as _calc_body.
    # The sample epochs use the two-argument tt_jd(whole, fraction) form:
    # collapsing t.tt ± dt into a single float64 quantizes the step to the JD
    # ULP (~4.6e-10 d at JD ~2.46e6), destroying the exact ±dt spacing that a
    # six-second step needs. Keeping dt in the fraction slot preserves it
    # through Skyfield's extended-precision time arithmetic.
    dt = _MOON_SPEED_HALF_STEP_DAYS if ipl == MOON else _BODY_SPEED_HALF_STEP_DAYS

    # calc_pctr's public flag convention gates velocity on FLG_SPEED only:
    # SPEED3 alone is echoed but leaves velocity slots zero. calc()/calc_ut()
    # use the separate SPEED3 three-position convention.
    if iflag & FLG_SPEED:
        ts_inner = get_timescale()
        t_prev = ts_inner.tt_jd(t.tt, -dt)
        t_next = ts_inner.tt_jd(t.tt, dt)

        # Sample flags: strip FLG_SPEED always; for ECLIPTIC sidereal strip
        # FLG_SIDEREAL so both samples are tropical (the ayanamsha drift is
        # applied to dlon below). For EQUATORIAL sidereal keep FLG_SIDEREAL so
        # the samples differentiate the reported (mean-equator) position — the
        # certified true-rate convention, matching _calc_body and the LEB path.
        sample_flags = iflag & ~(FLG_SPEED | FLG_SPEED3)
        if not is_equatorial:
            sample_flags &= ~FLG_SIDEREAL
        forced_barycenters = set(_spk_center_forced_barycenters.get())
        for center_target in (target, observer):
            if isinstance(center_target, _SpkCenterTarget) and (
                center_target._speed_stencil_requires_fallback(dt)
            ):
                forced_barycenters.add(center_target._barycenter_name)
        center_token = _spk_center_forced_barycenters.set(frozenset(forced_barycenters))
        try:
            result_prev, _ = _calc_body_pctr(t_prev, ipl, iplctr, sample_flags)
            result_next, _ = _calc_body_pctr(t_next, ipl, iplctr, sample_flags)
        finally:
            _spk_center_forced_barycenters.reset(center_token)
        p1_prev, p2_prev, p3_prev = result_prev[0], result_prev[1], result_prev[2]
        p1_next, p2_next, p3_next = result_next[0], result_next[1], result_next[2]

        # Calculate derivatives using central difference
        dp1 = (p1_next - p1_prev) / (2.0 * dt)
        dp2 = (p2_next - p2_prev) / (2.0 * dt)
        dp3 = (p3_next - p3_prev) / (2.0 * dt)

        # Handle longitude wrap-around (same threshold as _calc_body)
        wrap_threshold = 180.0 / (2.0 * dt)
        if dp1 > wrap_threshold:
            dp1 -= 360.0 / (2.0 * dt)
        elif dp1 < -wrap_threshold:
            dp1 += 360.0 / (2.0 * dt)

    # Apply sidereal offset if requested (ecliptic only).
    # Routed through _get_ayanamsa_for_flags so FLG_NONUT/FLG_J2000 select
    # the proper ayanamsha variant, consistent with _calc_body.
    if is_sidereal and not is_equatorial:
        ayanamsa = _get_ayanamsa_for_flags(t.ut1, iflag)
        p1 = (p1 - ayanamsa) % 360.0

        # Correct velocity for ayanamsha rate if speed was calculated.
        # Shortest-arc delta: the ayanamsha is mod 360 and the two samples
        # can straddle the 0/360 wrap (star-anchored modes cross 0 on
        # supported dates), which would spike the speed by ~360/(2*dt).
        if iflag & FLG_SPEED:
            ayanamsa_prev = _get_ayanamsa_for_flags(t.ut1 - dt, iflag)
            ayanamsa_next = _get_ayanamsa_for_flags(t.ut1 + dt, iflag)
            da = (ayanamsa_next - ayanamsa_prev + 180.0) % 360.0 - 180.0
            dp1 -= da / (2.0 * dt)

    return _to_native_floats((p1, p2, p3, dp1, dp2, dp3)), iflag


class NutationFallbackWarning(UserWarning):
    """Warning for degraded nutation precision (legacy, kept for API compatibility).

    .. deprecated::
        This warning class is retained for backward API compatibility but is
        no longer issued in normal operation. LibEphemeris now uses pyerfa
        (IAU 2006/2000A, ~0.01-0.05 mas) as a required dependency for all
        nutation calculations.

    See Also:
        get_nutation_model: Check which nutation model is currently active
    """

    pass


def get_nutation_model() -> dict:
    """Check which nutation model is currently active.

    LibEphemeris uses the IAU 2006/2000A nutation model via pyerfa
    (erfa.nut06a) for all nutation calculations. This provides the
    highest precision currently adopted by the IAU (~0.01-0.05 mas).

    Returns:
        dict: A dictionary containing:
            - ``model`` (str): "IAU2006_2000A" (pyerfa erfa.nut06a)
            - ``terms`` (int): Number of terms (1365+ lunisolar+planetary)
            - ``precision`` (str): Expected precision ("~0.01-0.05 mas")
            - ``source`` (str): "pyerfa" (the underlying C library)
            - ``skyfield_available`` (bool): Always True (kept for backward compatibility)

    Examples:
        >>> import libephemeris as eph
        >>> info = eph.get_nutation_model()
        >>> print(f"Model: {info['model']}, precision: {info['precision']}")
    """
    return {
        "model": "IAU2006_2000A",
        "terms": 1365,
        "precision": "~0.01-0.05 mas",
        "source": "pyerfa",
        # Backward compatibility: Skyfield is always available (required dep)
        "skyfield_available": True,
    }


def _calc_nutation_obliquity(
    jd: float, iflag: int
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """
    Calculate nutation and obliquity data for ECL_NUT (-1).

    Uses the IAU 2006/2000A nutation model via pyerfa (erfa.nut06a) and the
    Vondrák 2011 long-term mean obliquity from the directly evaluated
    epsilon_A series. This is the long-range companion to the Vondrák
    precession used throughout the reduction and reproduces the public
    ECL_NUT mean-obliquity convention at remote epochs without extrapolating
    the short-range IAU 2006 polynomial.

    Args:
        jd: Julian Day in UT
        iflag: Calculation flags (not used for nutation)

    Returns:
        Tuple containing:
            - Data tuple: (true_obliquity, mean_obliquity, nutation_longitude, nutation_obliquity, 0, 0)
            - Return flag
    """
    from .state import get_timescale

    ts = get_timescale()
    t = ts.ut1_jd(jd)
    return _calc_nutation_obliquity_tt(t.tt, iflag)


def _calc_nutation_obliquity_tt(
    jd_tt: float, iflag: int
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """ECL_NUT data for a TT Julian Day (see _calc_nutation_obliquity)."""
    import math

    # Mean obliquity of the ecliptic (direct Vondrák 2011 epsilon_A series)
    mean_obliquity = vondrak_mean_obliquity_deg(jd_tt)

    # Nutation IAU 2006/2000A via pyerfa (~0.01-0.05 mas precision)
    dpsi_rad, deps_rad = erfa.nut06a(2451545.0, jd_tt - 2451545.0)
    delta_psi = math.degrees(dpsi_rad)
    delta_eps = math.degrees(deps_rad)

    # True obliquity = mean obliquity + nutation in obliquity
    true_obliquity = mean_obliquity + delta_eps

    # Return format: (true_obliquity, mean_obliquity, nutation_lon, nutation_obl, 0, 0)
    return (true_obliquity, mean_obliquity, delta_psi, delta_eps, 0.0, 0.0), iflag


def _precess_ecliptic_state(
    lon: float,
    lat: float,
    dlon: float,
    dlat: float,
    jd_tt: float,
    to_j2000: bool,
) -> Tuple[float, float, float, float]:
    """Precess an ecliptic (lon, lat) position *and* its (dlon, dlat) rate.

    Transforms between the ecliptic of date at ``jd_tt`` and the J2000
    ecliptic. Unlike a bare position precession, the longitude/latitude
    rates are carried through consistently so the returned velocity is the
    exact time derivative of the returned position in the target frame.

    The rate is obtained by finite-differencing the full precession map
    with the *of-date epoch varying* between samples. That varying epoch is
    what makes the transform capture the ~50.29"/yr drift of the equinox
    (general precession) in addition to the fixed rotation of the axes: a
    body's longitude speed therefore differs between the of-date and J2000
    frames by the precession rate (~0.14"/day), matching the reference API.

    Args:
        lon: Ecliptic longitude in degrees (source frame).
        lat: Ecliptic latitude in degrees (source frame).
        dlon: Longitude rate in degrees/day (source frame).
        dlat: Latitude rate in degrees/day (source frame).
        jd_tt: Julian Day (TT) defining the of-date frame.
        to_j2000: True to map of-date -> J2000; False for J2000 -> of-date.

    Returns:
        Tuple ``(lon2, lat2, dlon2, dlat2)`` in the target frame.
    """
    from .astrometry import _precess_ecliptic

    J2000 = 2451545.0
    dt = 0.5
    if to_j2000:
        lon2, lat2 = _precess_ecliptic(lon, lat, jd_tt, J2000)
        lo_m, la_m = _precess_ecliptic(
            lon - dlon * dt, lat - dlat * dt, jd_tt - dt, J2000
        )
        lo_p, la_p = _precess_ecliptic(
            lon + dlon * dt, lat + dlat * dt, jd_tt + dt, J2000
        )
    else:
        lon2, lat2 = _precess_ecliptic(lon, lat, J2000, jd_tt)
        lo_m, la_m = _precess_ecliptic(
            lon - dlon * dt, lat - dlat * dt, J2000, jd_tt - dt
        )
        lo_p, la_p = _precess_ecliptic(
            lon + dlon * dt, lat + dlat * dt, J2000, jd_tt + dt
        )
    d = lo_p - lo_m
    if d > 180.0:
        d -= 360.0
    elif d < -180.0:
        d += 360.0
    return lon2, lat2, d / (2.0 * dt), (la_p - la_m) / (2.0 * dt)


def _apply_hypothetical_topocentric(
    result: tuple, jd_tt: float, iflag: int, result_is_j2000: bool, planets_dict
) -> tuple:
    """Apply the observer's geocentric parallax to a hypothetical body.

    The Uranian / Transpluto / Vulcan / Waldemath / predicted-planet paths
    build a geocentric ecliptic position by hand and previously ignored
    FLG_TOPOCTR entirely — silently returning the geocentric place labelled as
    topocentric. The parallax is real (~0.4° for the close Waldemath Moon,
    ~166" for the White Moon), so subtract the observer's geocentric offset
    (Earth centre -> topocentre), position and velocity, in the same ecliptic
    frame the ``result`` is currently expressed in.

    No-op for heliocentric / barycentric requests (parallax is undefined
    there). Raises when FLG_TOPOCTR is set but no observer location was
    configured, matching the reference-API behaviour on the planets.

    Args:
        result: (lon, lat, dist, dlon, dlat, ddist) geocentric ecliptic.
        jd_tt: Julian Day (TT).
        iflag: Calculation flags.
        result_is_j2000: True if ``result`` is J2000 ecliptic; False if it is
            the ecliptic of date (the frame governs how the observer offset is
            oriented before the subtraction).
        planets_dict: Skyfield planets dict (for the Earth vector).

    Returns:
        The topocentric (lon, lat, dist, dlon, dlat, ddist), or ``result``
        unchanged when FLG_TOPOCTR is not requested / not applicable.
    """
    if not (iflag & FLG_TOPOCTR):
        return result
    if iflag & (FLG_HELCTR | FLG_BARYCTR):
        return result

    observer_topo = get_topo()
    if observer_topo is None:
        from .exceptions import ConfigurationError

        raise ConfigurationError(
            "FLG_TOPOCTR requires a geographic position: call "
            "set_topo(lon, lat, alt) first",
            missing_config="observer_location",
            suggestion="Call set_topo(lon, lat, alt) first",
        )

    from .rebound_integration import _icrs_to_ecliptic_j2000
    from .astrometry import _precess_ecliptic

    # Observer geocentric offset (topocentre - geocentre) in equatorial ICRS,
    # position (AU) and velocity (AU/day, i.e. Earth's diurnal rotation).
    t = get_timescale().tt_jd(jd_tt)
    earth = planets_dict["earth"]
    obs = (earth + observer_topo).at(t)
    geo = earth.at(t)
    op = obs.position.au - geo.position.au
    ov = obs.velocity.au_per_d - geo.velocity.au_per_d
    ox, oy, oz = _icrs_to_ecliptic_j2000(float(op[0]), float(op[1]), float(op[2]))
    ovx, ovy, ovz = _icrs_to_ecliptic_j2000(float(ov[0]), float(ov[1]), float(ov[2]))

    def _rotate_to_of_date(vx: float, vy: float, vz: float) -> tuple:
        """Rotate an ecliptic-J2000 vector to the ecliptic of date, preserving
        its magnitude (exact for the direction; the offset is ~1 Earth radius,
        so any residual is negligible)."""
        mag = math.sqrt(vx * vx + vy * vy + vz * vz)
        if mag == 0.0:
            return (0.0, 0.0, 0.0)
        lon_v = math.degrees(math.atan2(vy, vx))
        lat_v = math.degrees(math.asin(max(-1.0, min(1.0, vz / mag))))
        lon_d, lat_d = _precess_ecliptic(lon_v, lat_v, 2451545.0, jd_tt)
        lr, br = math.radians(lon_d), math.radians(lat_d)
        cl = math.cos(br)
        return (mag * cl * math.cos(lr), mag * cl * math.sin(lr), mag * math.sin(br))

    if not result_is_j2000:
        ox, oy, oz = _rotate_to_of_date(ox, oy, oz)
        ovx, ovy, ovz = _rotate_to_of_date(ovx, ovy, ovz)

    lon, lat, dist, dlon, dlat, ddist = result
    lr, br = math.radians(lon), math.radians(lat)
    dlr, dbr = math.radians(dlon), math.radians(dlat)
    cl, sl = math.cos(lr), math.sin(lr)
    cb, sb = math.cos(br), math.sin(br)

    # Geocentric body state -> Cartesian (ecliptic of the result frame).
    bx = dist * cb * cl
    by = dist * cb * sl
    bz = dist * sb
    bvx = ddist * cb * cl - dist * sb * dbr * cl - dist * cb * sl * dlr
    bvy = ddist * cb * sl - dist * sb * dbr * sl + dist * cb * cl * dlr
    bvz = ddist * sb + dist * cb * dbr

    # Topocentric = geocentric - observer offset (position and velocity).
    tx, ty, tz = bx - ox, by - oy, bz - oz
    tvx, tvy, tvz = bvx - ovx, bvy - ovy, bvz - ovz

    # Cartesian -> spherical state.
    r_xy2 = tx * tx + ty * ty
    r = math.sqrt(r_xy2 + tz * tz)
    if r == 0.0:
        return result
    r_xy = math.sqrt(r_xy2)
    lon_t = math.degrees(math.atan2(ty, tx)) % 360.0
    lat_t = math.degrees(math.asin(max(-1.0, min(1.0, tz / r))))
    ddist_t = (tx * tvx + ty * tvy + tz * tvz) / r
    dlon_t = math.degrees((tx * tvy - ty * tvx) / r_xy2) if r_xy2 > 0 else 0.0
    dlat_t = math.degrees((tvz - (tz / r) * ddist_t) / r_xy) if r_xy > 0 else 0.0
    return (lon_t, lat_t, r, dlon_t, dlat_t, ddist_t)


def _maybe_equatorial_convert(
    result: tuple, jd_tt: float, iflag: int, *, already_j2000: bool = False
) -> tuple:
    """Convert ecliptic coordinates to equatorial if FLG_EQUATORIAL is set.

    For bodies computed in ecliptic coordinates (lunar nodes, Lilith, etc.),
    this applies the ecliptic-to-equatorial transformation using the true
    obliquity of the ecliptic at the given date.

    When FLG_J2000 is set, precesses from date to J2000.0 before any
    equatorial conversion.

    Args:
        result: 6-tuple of (lon, lat, dist, dlon, dlat, ddist) in ecliptic coords
        jd_tt: Julian Day in Terrestrial Time (for obliquity calculation)
        iflag: Calculation flags bitmask
        already_j2000: The input coordinates are already J2000 ecliptic
            (natural frame of the hypothetical-body orbits): skip the
            date→J2000 precession step but keep the J2000 obliquity for the
            equatorial rotation, so EQ+J2000 output is a consistent J2000
            frame rather than a J2000-ecliptic/of-date-equator mix.

    Returns:
        If FLG_EQUATORIAL is set: transformed (RA, Dec, dist, dRA, dDec, ddist)
        If FLG_J2000 is set: precessed coordinates in J2000 frame
        Otherwise: original result unchanged
    """
    lon, lat, dist, dlon, dlat, ddist = result
    if not (math.isfinite(lon) and math.isfinite(lat)):
        # Degenerate/undefined point (e.g. the interpolated lunar apsides in
        # nod_aps): no frame conversion is defined, so pass it through and
        # keep NaN as NaN instead of clamping to a pole declination.
        return result

    # Apply precession from date to J2000 if requested. The velocity is
    # precessed alongside the position (not left in the of-date frame) so the
    # J2000 speed is the true time derivative of the J2000 position — it
    # differs from the of-date speed by the general-precession rate. When the
    # equatorial conversion below runs, it then rotates a J2000-frame velocity
    # (not a frame-mixed one) through the J2000 obliquity.
    if (iflag & FLG_J2000) and not already_j2000:
        J2000 = 2451545.0
        if iflag & FLG_SPEED:
            lon, lat, dlon, dlat = _precess_ecliptic_state(
                lon, lat, dlon, dlat, jd_tt, to_j2000=True
            )
        else:
            from .astrometry import _precess_ecliptic

            lon, lat = _precess_ecliptic(lon, lat, jd_tt, J2000)

    if not (iflag & FLG_EQUATORIAL):
        return (lon, lat, dist, dlon, dlat, ddist)

    from .cache import get_true_obliquity
    from .utils import cotrans_sp

    # For J2000 frame, use J2000 obliquity; otherwise use obliquity of date.
    # Sidereal mode uses the mean obliquity and precession-only equatorial frame.
    if iflag & FLG_J2000:
        # IAU 2006 mean obliquity at J2000 (see
        # fast_calc.OBLIQUITY_J2000_REF_DEG).
        eps = 84381.406 / 3600.0
    elif (iflag & FLG_NONUT) or (iflag & FLG_SIDEREAL):
        # NONUT / SIDEREAL: use mean obliquity (no nutation correction)
        from .cache import get_mean_obliquity

        eps = get_mean_obliquity(jd_tt)
    else:
        eps = get_true_obliquity(jd_tt)

    # Negative obliquity rotates ecliptic coordinates to the equator.
    result = cotrans_sp((lon, lat, dist, dlon, dlat, ddist), -eps)
    return (
        result[0],
        result[1],
        result[2],
        result[3],
        result[4],
        result[5],
    )


def _keplerian_position_at(
    jd_tt: float, ipl: int, iflag: int, planets_dict
) -> Tuple[float, float, float]:
    """Compute ecliptic position of a minor body via Keplerian fallback.

    Returns (lon, lat, dist) in ecliptic coordinates, either heliocentric
    or geocentric depending on iflag.

    Args:
        jd_tt: Julian Day in Terrestrial Time.
        ipl: Minor body ID (CHIRON, CERES, etc.).
        iflag: Calculation flags (FLG_HELCTR checked).
        planets_dict: Vector ephemeris resolved for the active mode.

    Returns:
        Tuple of (longitude_deg, latitude_deg, distance_au).
    """
    from . import minor_bodies
    from .state import get_timescale

    if ipl in minor_bodies.MINOR_BODY_ELEMENTS:
        lon_hel, lat_hel, r_hel = minor_bodies.calc_minor_body_heliocentric(ipl, jd_tt)
    else:
        # Arbitrary numbered asteroid: Keplerian elements fetched from
        # JPL SBDB (cached after the first call).  No element source at
        # all means the body is genuinely unknown.
        from .exceptions import UnknownBodyError

        try:
            lon_hel, lat_hel, r_hel = minor_bodies.calc_asteroid_by_number(
                ipl - AST_OFFSET, jd_tt, use_spk=False
            )
        except (ValueError, ConnectionError, OSError) as e:
            raise UnknownBodyError(
                message=(
                    f"Asteroid {ipl - AST_OFFSET} (body ID {ipl}): no SPK "
                    f"registered and JPL SBDB lookup failed ({e}). Register "
                    f"an SPK file or enable auto-download, or check the "
                    f"asteroid number."
                ),
                body_id=ipl,
            ) from e

    # The minor-body helpers return the TRUE ecliptic of date. The
    # nutation-free frames (NONUT, J2000 — whose downstream precession is
    # nutation-blind — and SIDEREAL+EQUATORIAL) need Δψ stripped, exactly
    # like every other body path.
    from .cache import get_cached_nutation

    _sid_eq = bool(iflag & FLG_SIDEREAL) and bool(iflag & FLG_EQUATORIAL)
    _nut_free = bool(iflag & (FLG_NONUT | FLG_J2000)) or _sid_eq

    if iflag & FLG_HELCTR:
        if _nut_free:
            dpsi_rad, _ = get_cached_nutation(jd_tt)
            lon_hel = (lon_hel - math.degrees(dpsi_rad)) % 360.0
        return lon_hel, lat_hel, r_hel

    # Convert heliocentric ecliptic spherical to Cartesian
    lon_rad = math.radians(lon_hel)
    lat_rad = math.radians(lat_hel)
    x_hel_ecl = r_hel * math.cos(lat_rad) * math.cos(lon_rad)
    y_hel_ecl = r_hel * math.cos(lat_rad) * math.sin(lon_rad)
    z_hel_ecl = r_hel * math.sin(lat_rad)

    # Get Earth position in ecliptic frame
    ts = get_timescale()
    t = ts.tt_jd(jd_tt)
    earth = planets_dict["earth"]
    sun = planets_dict["sun"]
    earth_helio = sun.at(t).observe(earth)
    earth_xyz_ecl = earth_helio.frame_xyz(ecliptic_frame).au

    # Geocentric = minor body heliocentric - Earth heliocentric
    x_geo_ecl = x_hel_ecl - earth_xyz_ecl[0]
    y_geo_ecl = y_hel_ecl - earth_xyz_ecl[1]
    z_geo_ecl = z_hel_ecl - earth_xyz_ecl[2]

    # Convert geocentric Cartesian back to spherical
    r_geo = math.sqrt(x_geo_ecl**2 + y_geo_ecl**2 + z_geo_ecl**2)
    lon = math.degrees(math.atan2(y_geo_ecl, x_geo_ecl)) % 360.0
    lat = math.degrees(math.asin(z_geo_ecl / r_geo)) if r_geo > 0 else 0.0

    if _nut_free:
        dpsi_rad, _ = get_cached_nutation(jd_tt)
        lon = (lon - math.degrees(dpsi_rad)) % 360.0

    return lon, lat, r_geo


def _calc_keplerian_fallback(t, ipl: int, iflag: int, planets):
    """Return the existing curated code-model result with honest provenance."""
    jd_tt = t.tt
    lon, lat, dist = _keplerian_position_at(jd_tt, ipl, iflag, planets)

    speed_lon = 0.0
    speed_lat = 0.0
    speed_dist = 0.0
    if iflag & FLG_SPEED:
        dt = 1.0 / 86400.0
        lon_prev, lat_prev, dist_prev = _keplerian_position_at(
            jd_tt - dt, ipl, iflag, planets
        )
        lon_next, lat_next, dist_next = _keplerian_position_at(
            jd_tt + dt, ipl, iflag, planets
        )
        speed_lon = (lon_next - lon_prev) / (2.0 * dt)
        speed_lat = (lat_next - lat_prev) / (2.0 * dt)
        speed_dist = (dist_next - dist_prev) / (2.0 * dt)

        if speed_lon > 180.0 / (2.0 * dt):
            speed_lon -= 360.0 / (2.0 * dt)
        if speed_lon < -180.0 / (2.0 * dt):
            speed_lon += 360.0 / (2.0 * dt)

    if (iflag & FLG_SIDEREAL) and not (iflag & FLG_EQUATORIAL):
        lon, speed_lon = _apply_sidereal_correction(lon, speed_lon, t.ut1, iflag)

    result = _to_native_floats(
        _maybe_equatorial_convert(
            (lon, lat, dist, speed_lon, speed_lat, speed_dist), jd_tt, iflag
        )
    )
    _mark_dispatch_source("Keplerian")
    return result, iflag


def _assist_state_at(jd_tt: float, ipl: int):
    """Propagate a minor body with ASSIST and return its heliocentric state."""
    from . import minor_bodies
    from .rebound_integration import propagate_orbit_assist

    elements = minor_bodies.MINOR_BODY_ELEMENTS[ipl]
    return propagate_orbit_assist(elements, elements.epoch, jd_tt)


def _assist_shift_state(result, jd_tt: float):
    """Taylor-shift an ASSIST state over a short interval in days."""
    from .rebound_integration import PropagationResult

    dt = jd_tt - result.jd_tt
    return PropagationResult(
        x=result.x + result.vx * dt,
        y=result.y + result.vy * dt,
        z=result.z + result.vz * dt,
        vx=result.vx,
        vy=result.vy,
        vz=result.vz,
        jd_tt=jd_tt,
    )


def _assist_position_from_state(
    result, jd_tt: float, iflag: int, planets_dict
) -> Tuple[float, float, float]:
    """Convert an ASSIST heliocentric state to the requested ecliptic position."""
    from .state import get_timescale

    lon_hel = result.ecliptic_lon
    lat_hel = result.ecliptic_lat
    r_hel = result.distance

    if iflag & FLG_HELCTR:
        # result.ecliptic_lon/lat are heliocentric ecliptic J2000 (rebound_
        # integration rotates ICRS->ecliptic with the fixed J2000 obliquity).
        # Every other body branch hands the caller an of-date true-ecliptic
        # longitude (which _maybe_equatorial_convert then reframes to J2000 /
        # equatorial per the flags); the geocentric branch below does the same.
        # Precess to the ecliptic of date and add nutation in longitude under
        # the identical NONUT / J2000 / SIDEREAL+EQUATORIAL rules, so
        # heliocentric follows the same contract. Returning the raw J2000
        # longitude froze it at the J2000 equinox and drifted from every other
        # path by the accumulated general precession (~7 deg / 25000" by 2500).
        from .astrometry import precess_from_j2000, C_LIGHT_AU_DAY
        from .cache import get_cached_nutation

        # Light-time retardation (body -> Sun), matching the SPK heliocentric
        # path (get_vector retards the target by dist/c for HELCTR). Without it
        # the ASSIST heliocentric place is the *geometric* position and steps by
        # ~1-2" against the retarded SPK place at the SPK/ASSIST coverage
        # boundary. Skipped for FLG_TRUEPOS (geometric position requested).
        px, py, pz = result.x, result.y, result.z
        if not (iflag & FLG_TRUEPOS):
            vx, vy, vz = result.vx, result.vy, result.vz
            lt = 0.0
            for _ in range(2):
                rx, ry, rz = px - vx * lt, py - vy * lt, pz - vz * lt
                lt = math.sqrt(rx * rx + ry * ry + rz * rz) / C_LIGHT_AU_DAY
            px, py, pz = px - vx * lt, py - vy * lt, pz - vz * lt
            r_hel = math.sqrt(px * px + py * py + pz * pz)
            lon_hel = math.degrees(math.atan2(py, px)) % 360.0
            lat_hel = math.degrees(math.asin(pz / r_hel)) if r_hel > 0 else 0.0

        lon_od, lat_od = precess_from_j2000(lon_hel, lat_hel, jd_tt)
        _sid_eq = bool(iflag & FLG_SIDEREAL) and bool(iflag & FLG_EQUATORIAL)
        if not (iflag & FLG_NONUT) and not (iflag & FLG_J2000) and not _sid_eq:
            dpsi_rad, _ = get_cached_nutation(jd_tt)
            lon_od = (lon_od + math.degrees(dpsi_rad)) % 360.0
        return lon_od, lat_od, r_hel

    ts = get_timescale()
    t = ts.tt_jd(jd_tt)
    sun = planets_dict["sun"]
    earth = planets_dict["earth"]

    # The ASSIST result is heliocentric ecliptic J2000 (rebound_integration
    # rotates ICRS->ecliptic with the fixed J2000 obliquity), with position AND
    # velocity available. Build the geocentric apparent place exactly like the
    # SPK type-21 path: everything in ecliptic J2000, geometric heliocentric
    # Earth, light-time retardation, then stellar aberration, then precession to
    # of-date + nutation. The previous code differenced an of-date Earth
    # (Skyfield ecliptic_frame) from the J2000 asteroid (a ~12 arcmin frame-mix)
    # and applied neither light-time nor aberration (~25 arcsec more).
    from .rebound_integration import _icrs_to_ecliptic_j2000
    from .astrometry import (
        C_LIGHT_AU_DAY,
        apply_aberration_to_position,
        precess_from_j2000,
    )
    from .cache import get_cached_nutation

    ast_pos = [result.x, result.y, result.z]  # heliocentric ecliptic J2000, AU
    ast_vel = [result.vx, result.vy, result.vz]  # AU/day

    earth_helio_icrs = earth.at(t).position.au - sun.at(t).position.au
    earth_helio_ecl = _icrs_to_ecliptic_j2000(
        float(earth_helio_icrs[0]),
        float(earth_helio_icrs[1]),
        float(earth_helio_icrs[2]),
    )

    apply_light_time = not (iflag & FLG_TRUEPOS)
    apply_aberration = not (iflag & FLG_TRUEPOS) and not (iflag & FLG_NOABERR)

    # Geocentric vector with light-time retardation. The asteroid's heliocentric
    # velocity is effectively constant over the ~minutes-to-hours light time, so
    # a first-order retardation (pos - vel*lt) matches re-integrating to
    # jd_tt - lt without the extra N-body cost. Two iterations converge.
    geo = [ast_pos[i] - earth_helio_ecl[i] for i in range(3)]
    if apply_light_time:
        for _ in range(2):
            lt = math.sqrt(geo[0] ** 2 + geo[1] ** 2 + geo[2] ** 2) / C_LIGHT_AU_DAY
            geo = [
                (ast_pos[i] - ast_vel[i] * lt) - earth_helio_ecl[i] for i in range(3)
            ]

    if apply_aberration:
        # Aberration uses Earth's barycentric velocity in ecliptic J2000, per the
        # IAU apparent-place convention shared with the SPK/planets pipeline.
        earth_vel_bary_ecl = _icrs_to_ecliptic_j2000(
            *(float(v) for v in earth.at(t).velocity.au_per_d)
        )
        geo = list(
            apply_aberration_to_position(
                (geo[0], geo[1], geo[2]),
                (
                    earth_vel_bary_ecl[0],
                    earth_vel_bary_ecl[1],
                    earth_vel_bary_ecl[2],
                ),
            )
        )

    r_geo = math.sqrt(geo[0] ** 2 + geo[1] ** 2 + geo[2] ** 2)
    lon_j2000 = math.degrees(math.atan2(geo[1], geo[0])) % 360.0
    lat_j2000 = math.degrees(math.asin(geo[2] / r_geo)) if r_geo > 0 else 0.0

    # The geocentric vector is in the J2000 ecliptic; precess to the ecliptic of
    # date and add nutation in longitude (Δψ) so the return matches the of-date
    # true-ecliptic contract every other body branch follows (the caller then
    # applies the J2000 / equatorial framing via _maybe_equatorial_convert).
    # Δψ is omitted for NONUT (mean ecliptic) and for J2000 — there the caller
    # precesses back to J2000, whose ecliptic carries no nutation. It is also
    # omitted for SIDEREAL+EQUATORIAL: _maybe_equatorial_convert rotates that
    # case with the MEAN obliquity, so the longitude must be the mean ecliptic
    # too (mirrors the node/Lilith _sid_eq rule), else ~Δψ leaks into RA/Dec.
    lon, lat = precess_from_j2000(lon_j2000, lat_j2000, jd_tt)
    _sid_eq = bool(iflag & FLG_SIDEREAL) and bool(iflag & FLG_EQUATORIAL)
    if not (iflag & FLG_NONUT) and not (iflag & FLG_J2000) and not _sid_eq:
        dpsi_rad, _ = get_cached_nutation(jd_tt)
        lon = (lon + math.degrees(dpsi_rad)) % 360.0

    return lon, lat, r_geo


def _assist_position_at(
    jd_tt: float, ipl: int, iflag: int, planets_dict
) -> Tuple[float, float, float]:
    """Compute ecliptic position of a minor body via ASSIST N-body integration.

    Returns (lon, lat, dist) in ecliptic coordinates, either heliocentric
    or geocentric depending on iflag.
    """
    result = _assist_state_at(jd_tt, ipl)
    return _assist_position_from_state(result, jd_tt, iflag, planets_dict)


def _apply_sidereal_correction(
    lon: float, dlon: float, ut1: float, iflag: int, sid_mode: int | None = None
) -> Tuple[float, float]:
    """Apply sidereal ayanamsa correction to ecliptic longitude and its velocity.

    Subtracts the ayanamsa from the longitude and corrects the velocity for
    the ayanamsa precession rate using a central-difference derivative.

    Args:
        lon: Ecliptic longitude in degrees (tropical).
        dlon: Longitude velocity in degrees/day.
        ut1: Julian Day in UT1.
        iflag: Calculation flags (checked for FLG_SPEED).
        sid_mode: Sidereal mode override; reads the global state when None.
            An explicit value lets thread-safe callers (EphemerisContext,
            Horizons) avoid swapping the global sidereal mode.

    Returns:
        Tuple of (corrected_lon, corrected_dlon).
    """
    ayanamsa = _get_ayanamsa_for_flags(ut1, iflag, sid_mode)
    lon = (lon - ayanamsa) % 360.0
    if iflag & FLG_SPEED:
        dt_aya = 1.0 / 86400.0
        ayanamsa_prev = _get_ayanamsa_for_flags(ut1 - dt_aya, iflag, sid_mode)
        ayanamsa_next = _get_ayanamsa_for_flags(ut1 + dt_aya, iflag, sid_mode)
        # Unwrap the delta into (-180, 180]: the ayanamsha is mod 360, and the
        # two samples can straddle the 0/360 wrap, which would turn the raw
        # -360 difference into a ~1.6e7 deg/day speed spike.
        d_aya = (ayanamsa_next - ayanamsa_prev + 180.0) % 360.0 - 180.0
        dlon -= d_aya / (2.0 * dt_aya)
    return lon, dlon


def _degenerate_origin_result() -> Tuple[float, float, float, float, float, float]:
    """Position tuple for a body that sits at its own observation origin.

    Earth observed geocentrically and the Sun observed heliocentrically are
    both trivially the zero vector: latitude, distance and every speed are zero
    and the longitude direction is physically meaningless. Measured reference
    behavior: the reference API returns the exact zero vector here for
    every frame -- including a sidereal ecliptic request, where the ayanamsha
    is NOT subtracted from the (undefined) zero-length longitude. Mirror that;
    the LEB fast path skips the same subtraction on a zero-length state (see
    ``fast_calc._fast_calc_core``).

    Returns:
        The 6-tuple ``(0, 0, 0, 0, 0, 0)``.
    """
    return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


def _calc_body(
    t, ipl: int, iflag: int
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """
    Calculate position of any celestial body or point (internal dispatcher).

    This is the core calculation function that routes requests to appropriate
    sub-modules based on body type. Supports all reference API body types.

    Supported body types:
        - Classical planets (Sun, Moon, Mercury-Pluto) via JPL DE440 ephemeris
        - Lunar nodes (Mean/True NORTH) via lunar.py. The south nodes
          (negative ids) are NOT handled here: the public entry points
          (calc/calc_ut and EphemerisContext) intercept them and derive the
          antipode via planets._south_node_from_north, so _calc_body only ever
          receives the positive (north) node id.
        - Lilith/Lunar apogee (Mean/Osculating) via lunar.py
        - Minor bodies (asteroids, TNOs) via minor_bodies.py with rigorous geocentric conversion
        - Fixed stars (Regulus, Spica) via fixed_stars.py
        - Astrological angles (ASC, MC, Vertex, etc.) via angles.py
        - Arabic parts (Fortune, Spirit, etc.) via arabic_parts.py

    Args:
        t: Skyfield Time object (UT1 or TT)
        ipl: Planet/body ID (SUN, MOON, MARS, etc.)
        iflag: Calculation flags bitmask (FLG_SPEED, FLG_HELCTR, etc.)

    Returns:
        Tuple containing:
            - Position tuple: (lon, lat, dist, speed_lon, speed_lat, speed_dist)
            - Return flag: iflag value on success

    Coordinate Systems:
        - Ecliptic (default): Longitude/Latitude relative to ecliptic plane
        - Equatorial (FLG_EQUATORIAL): Right Ascension/Declination
        - Heliocentric (FLG_HELCTR): Sun-centered coordinates
        - Topocentric (FLG_TOPOCTR): Observer location on Earth surface
        - Sidereal (FLG_SIDEREAL): Fixed zodiac (requires set_sid_mode)

    Precision:
        - Minor body geocentric conversion uses Skyfield's frame transformations
        - Properly handles precession, nutation, and true obliquity of date
    """
    from . import (
        lunar,
        minor_bodies,
        fixed_stars,
        angles,
        arabic_parts,
        planetary_moons,
    )
    from .state import get_angles_cache

    from .state import _get_computation_ephemeris

    planets = _get_computation_ephemeris()

    # Sun heliocentric = Sun from Sun = trivially (0,0,0). A sidereal ecliptic
    # request keeps the zero vector: measured reference behavior does not
    # subtract the ayanamsha from the undefined zero-length longitude
    # (see _degenerate_origin_result).
    if ipl == SUN and (iflag & FLG_HELCTR):
        _mark_dispatch_source("Analytical")
        return _to_native_floats(_degenerate_origin_result()), iflag

    # Remap AST_OFFSET + N to dedicated body IDs for built-in asteroids
    # (idempotent safety net; calc/calc_ut already remap pre-dispatch).
    ipl = _remap_ast_offset(ipl)

    # Handle planetary moons (Galilean moons, Titan, etc.)
    if planetary_moons.is_planetary_moon(ipl):
        from .state import get_calc_mode

        if get_calc_mode() == "leb":
            from .exceptions import Error

            raise Error(
                "planetary-moon SPK access is disabled in calculation mode 'leb'; "
                "no declared LEB or analytical model serves this body"
            )
        result = planetary_moons.calc_moon_position(t, ipl, iflag)
        if result is not None:
            result = _maybe_equatorial_convert(result, t.tt, iflag)
            _mark_dispatch_source("SPK")
            return result, iflag
        # No registered satellite SPK covers this moon.
        from .exceptions import Error

        raise Error(
            f"no satellite SPK registered for "
            f"{planetary_moons.get_moon_name(ipl)} (body {ipl}); "
            f"call register_moon_spk() with e.g. jup365.bsp/sat441.bsp"
        )

    # Public convention for Moon-derived abstract points under HELCTR/BARYCTR:
    # return an all-zero position because nodes/apses have no heliocentric or
    # barycentric place, while echoing the normal retflag.
    if ipl in (MEAN_NODE, TRUE_NODE, MEAN_APOG, OSCU_APOG, INTP_APOG, INTP_PERG) and (
        iflag & (FLG_HELCTR | FLG_BARYCTR)
    ):
        _mark_dispatch_source("Analytical")
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0), iflag

    # Handle lunar nodes (Mean/True North/South)
    if ipl in [MEAN_NODE, TRUE_NODE]:
        jd_tt = t.tt
        is_sidereal = bool(iflag & FLG_SIDEREAL)
        if ipl == MEAN_NODE:
            lon, lat, dist, dlon, dlat, ddist = lunar.calc_mean_lunar_node_state(jd_tt)
            # The clean-room model returns the mean ecliptic of date.
            # The reference ephemeris outputs in the true ecliptic of date,
            # so add nutation in longitude (dpsi) unless NONUT is set.
            # When SIDEREAL+EQUATORIAL, the reference API outputs mean ecliptic
            # (no nutation) converted with mean obliquity, so skip dpsi.
            # J2000 output is likewise the mean ecliptic precessed to J2000:
            # the reference treats FLG_J2000 as implying no nutation here.
            _sid_eq = is_sidereal and bool(iflag & FLG_EQUATORIAL)
            _no_nut = bool(iflag & (FLG_NONUT | FLG_J2000)) or _sid_eq
            if not _no_nut:
                from .cache import get_cached_nutation

                dpsi_rad, _ = get_cached_nutation(jd_tt)
                lon = (lon + math.degrees(dpsi_rad)) % 360.0
            if iflag & FLG_SPEED:
                # True-ecliptic output includes the nutation rate just as the
                # longitude includes nutation itself.  The model supplies the
                # intrinsic mean-ecliptic speed directly.
                if not _no_nut:
                    dlon += _nutation_rate_deg_per_day(jd_tt, 0.001)
            else:
                dlon = dlat = ddist = 0.0
            result = (lon, lat, dist, dlon, dlat, ddist)
            # Mean points compose SIDEREAL+J2000 by first expressing the
            # tropical state on the J2000 ecliptic and only then subtracting
            # the mean ayanamsha. Public outputs show that this ordering is
            # specific to the mean node/apogee; the rotations do not commute.
            _sid_j2000_ecl = (
                is_sidereal
                and bool(iflag & FLG_J2000)
                and not bool(iflag & FLG_EQUATORIAL)
            )
            if _sid_j2000_ecl:
                result = _maybe_equatorial_convert(result, jd_tt, iflag)
                lon, lat, dist, dlon, dlat, ddist = result
                lon, dlon = _apply_sidereal_correction(lon, dlon, t.ut1, iflag)
                result = (lon, lat, dist, dlon, dlat, ddist)
            else:
                if is_sidereal and not (iflag & FLG_EQUATORIAL):
                    lon, dlon = _apply_sidereal_correction(lon, dlon, t.ut1, iflag)
                    result = (lon, lat, dist, dlon, dlat, ddist)
                result = _maybe_equatorial_convert(result, jd_tt, iflag)
            return _to_native_floats(result), iflag
        else:  # TRUE_NODE
            calc_iflag = iflag
            lon, lat, dist = lunar.calc_true_lunar_node(jd_tt)
            # TrueNode includes nutation effects in its perturbation terms.
            # When NONUT is set, subtract dpsi to get mean ecliptic position.
            # When SIDEREAL+EQUATORIAL, the reference API also outputs mean ecliptic
            # (no nutation) converted with mean obliquity, so strip dpsi too.
            # J2000 output is likewise nutation-free (mean ecliptic precessed
            # to J2000), matching the reference.
            _sid_eq = is_sidereal and bool(calc_iflag & FLG_EQUATORIAL)
            _no_nut = bool(calc_iflag & (FLG_NONUT | FLG_J2000)) or _sid_eq
            if _no_nut:
                from .cache import get_cached_nutation

                dpsi_rad, _ = get_cached_nutation(jd_tt)
                lon = (lon - math.degrees(dpsi_rad)) % 360.0
            # Calculate velocity via central difference numerical differentiation
            # The ordinary 0.05-day central-difference half-step resolves the
            # smooth osculating curve; SPEED3/topocentric requests use their
            # finer public sampling convention.
            dlon, dlat, ddist = 0.0, 0.0, 0.0
            if calc_iflag & FLG_SPEED:
                dt = 0.0001 if calc_iflag & (FLG_SPEED3 | FLG_TOPOCTR) else 0.05
                try:
                    lon_prev, lat_prev, dist_prev = lunar.calc_true_lunar_node(
                        jd_tt - dt
                    )
                    lon_next, lat_next, dist_next = lunar.calc_true_lunar_node(
                        jd_tt + dt
                    )
                    # Handle longitude wrap-around before computing velocity
                    lon_diff = lon_next - lon_prev
                    if lon_diff > 180:
                        lon_diff -= 360.0
                    elif lon_diff < -180:
                        lon_diff += 360.0
                    dlon = lon_diff / (2.0 * dt)
                    dlat = (lat_next - lat_prev) / (2.0 * dt)
                    ddist = (dist_next - dist_prev) / (2.0 * dt)
                    # The raw true-node curve contains nutation; under
                    # NONUT (mean-ecliptic output) remove its rate, the
                    # same way dpsi was removed from the position.
                    if _no_nut:
                        dlon -= _nutation_rate_deg_per_day(jd_tt, dt)
                except (
                    IndexError,
                    ValueError,
                    ArithmeticError,
                    SkyfieldRangeError,
                ):
                    # At ephemeris boundaries, speed calculation may fail
                    # (including Skyfield range errors from ±dt offsets)
                    # Return 0 for speed components
                    from .logging_config import get_logger

                    get_logger().debug(
                        "True Node velocity fallback to 0 at jd=%.1f", jd_tt
                    )
            # FLG_J2000 is honored uniformly for every lunar point: under
            # SIDEREAL|J2000 precess the mean-ecliptic state to J2000 first,
            # then subtract the mean ayanamsha in that frame — the same
            # pipeline MEAN_NODE/MEAN_APOG use. Ayanamsha and ecliptic-plane
            # precession are independent, composable rotations; dropping one
            # would let |TrueNode - MeanNode| grow unbounded away from J2000.
            _sid_j2000_ecl = (
                is_sidereal
                and bool(iflag & FLG_J2000)
                and not bool(iflag & FLG_EQUATORIAL)
            )
            result = (lon, lat, dist, dlon, dlat, ddist)
            if _sid_j2000_ecl:
                result = _maybe_equatorial_convert(result, jd_tt, iflag)
                lon, lat, dist, dlon, dlat, ddist = result
                lon, dlon = _apply_sidereal_correction(lon, dlon, t.ut1, iflag)
                result = (lon, lat, dist, dlon, dlat, ddist)
            else:
                if is_sidereal and not (iflag & FLG_EQUATORIAL):
                    lon, dlon = _apply_sidereal_correction(lon, dlon, t.ut1, iflag)
                    result = (lon, lat, dist, dlon, dlat, ddist)
                result = _maybe_equatorial_convert(result, jd_tt, iflag)
            return _to_native_floats(result), iflag

    # Handle Lilith (Mean/Osculating Apogee)
    if ipl in [MEAN_APOG, OSCU_APOG]:
        jd_tt = t.tt
        is_sidereal = bool(iflag & FLG_SIDEREAL)
        if ipl == MEAN_APOG:
            lon, lat, dist, dlon, dlat, ddist = lunar.calc_mean_lilith_state(jd_tt)
            # The clean-room model returns the mean ecliptic of date.
            # The reference ephemeris outputs in the true ecliptic of date,
            # so add nutation in longitude (dpsi) unless NONUT is set.
            # When SIDEREAL+EQUATORIAL, the reference API outputs mean ecliptic
            # (no nutation) converted with mean obliquity, so skip dpsi.
            # J2000 output is likewise the mean ecliptic precessed to J2000:
            # the reference treats FLG_J2000 as implying no nutation here.
            _sid_eq = is_sidereal and bool(iflag & FLG_EQUATORIAL)
            _no_nut = bool(iflag & (FLG_NONUT | FLG_J2000)) or _sid_eq
            if not _no_nut:
                from .cache import get_cached_nutation

                dpsi_rad, _ = get_cached_nutation(jd_tt)
                lon = (lon + math.degrees(dpsi_rad)) % 360.0
            if iflag & FLG_SPEED:
                if not _no_nut:
                    dlon += _nutation_rate_deg_per_day(jd_tt, 0.001)
            else:
                dlon = dlat = ddist = 0.0
            result = (lon, lat, dist, dlon, dlat, ddist)
            _sid_j2000_ecl = (
                is_sidereal
                and bool(iflag & FLG_J2000)
                and not bool(iflag & FLG_EQUATORIAL)
            )
            if _sid_j2000_ecl:
                result = _maybe_equatorial_convert(result, jd_tt, iflag)
                lon, lat, dist, dlon, dlat, ddist = result
                lon, dlon = _apply_sidereal_correction(lon, dlon, t.ut1, iflag)
                result = (lon, lat, dist, dlon, dlat, ddist)
            else:
                if is_sidereal and not (iflag & FLG_EQUATORIAL):
                    lon, dlon = _apply_sidereal_correction(lon, dlon, t.ut1, iflag)
                    result = (lon, lat, dist, dlon, dlat, ddist)
                result = _maybe_equatorial_convert(result, jd_tt, iflag)
            return _to_native_floats(result), iflag
        else:  # OSCU_APOG
            calc_iflag = iflag
            lon, lat, dist = lunar.calc_true_lilith(jd_tt)
            # OscuApog includes nutation effects in its orbital computation.
            # When NONUT is set, subtract dpsi to get mean ecliptic position.
            # When SIDEREAL+EQUATORIAL, the reference API also outputs mean ecliptic
            # (no nutation) converted with mean obliquity, so strip dpsi too.
            # J2000 output is likewise nutation-free (mean ecliptic precessed
            # to J2000), matching the reference.
            _sid_eq = is_sidereal and bool(calc_iflag & FLG_EQUATORIAL)
            _no_nut = bool(calc_iflag & (FLG_NONUT | FLG_J2000)) or _sid_eq
            if _no_nut:
                from .cache import get_cached_nutation

                dpsi_rad, _ = get_cached_nutation(jd_tt)
                lon = (lon - math.degrees(dpsi_rad)) % 360.0
            # Calculate velocity via central difference numerical differentiation.
            # The osculating (true) apogee carries fast short-period structure,
            # so a coarse 0.05-day half-step chords across the curvature and
            # biases the reported longitude speed by up to ~1"/day, leaving it
            # self-inconsistent with the reported position (and ~0.2-0.8"/day off
            # the reference). A 0.002-day half-step brings the central difference
            # self-consistent to <0.01"/day over 1896-2073 while staying well
            # clear of double-precision differencing noise (which only reappears
            # below ~5e-4 day). SPEED3/topocentric requests keep their still-finer
            # public sampling. Same two evaluations as before -- only dt changes.
            dlon, dlat, ddist = 0.0, 0.0, 0.0
            if calc_iflag & FLG_SPEED:
                dt = 0.0001 if calc_iflag & (FLG_SPEED3 | FLG_TOPOCTR) else 0.002
                try:
                    lon_prev, lat_prev, dist_prev = lunar.calc_true_lilith(jd_tt - dt)
                    lon_next, lat_next, dist_next = lunar.calc_true_lilith(jd_tt + dt)
                    # Handle longitude wrap-around before computing velocity
                    lon_diff = lon_next - lon_prev
                    if lon_diff > 180:
                        lon_diff -= 360.0
                    elif lon_diff < -180:
                        lon_diff += 360.0
                    dlon = lon_diff / (2.0 * dt)
                    dlat = (lat_next - lat_prev) / (2.0 * dt)
                    ddist = (dist_next - dist_prev) / (2.0 * dt)
                    # The raw osculating curve contains nutation; under
                    # NONUT (mean-ecliptic output) remove its rate.
                    if _no_nut:
                        dlon -= _nutation_rate_deg_per_day(jd_tt, dt)
                except (
                    IndexError,
                    ValueError,
                    ArithmeticError,
                    SkyfieldRangeError,
                ):
                    # At ephemeris boundaries, speed calculation may fail
                    # Return 0 for speed components
                    from .logging_config import get_logger

                    get_logger().debug(
                        "True Lilith velocity fallback to 0 at jd=%.1f", jd_tt
                    )
            # FLG_J2000 honored uniformly under SIDEREAL: precess to J2000,
            # then subtract the mean ayanamsha (MEAN_NODE/MEAN_APOG pipeline).
            _sid_j2000_ecl = (
                is_sidereal
                and bool(iflag & FLG_J2000)
                and not bool(iflag & FLG_EQUATORIAL)
            )
            result = (lon, lat, dist, dlon, dlat, ddist)
            if _sid_j2000_ecl:
                result = _maybe_equatorial_convert(result, jd_tt, iflag)
                lon, lat, dist, dlon, dlat, ddist = result
                lon, dlon = _apply_sidereal_correction(lon, dlon, t.ut1, iflag)
                result = (lon, lat, dist, dlon, dlat, ddist)
            else:
                if is_sidereal and not (iflag & FLG_EQUATORIAL):
                    lon, dlon = _apply_sidereal_correction(lon, dlon, t.ut1, iflag)
                    result = (lon, lat, dist, dlon, dlat, ddist)
                result = _maybe_equatorial_convert(result, jd_tt, iflag)
            return _to_native_floats(result), iflag

    # Handle Interpolated Apogee/Perigee
    if ipl in [INTP_APOG, INTP_PERG]:
        jd_tt = t.tt
        is_sidereal = bool(iflag & FLG_SIDEREAL)
        calc_iflag = iflag
        use_speed3_state = bool(calc_iflag & FLG_SPEED3)
        use_speed3_velocity = bool(calc_iflag & FLG_SPEED) and (
            use_speed3_state or bool(calc_iflag & FLG_TOPOCTR)
        )
        lon, lat, dist, dlon, dlat, ddist = lunar.calc_interpolated_apse_state(
            jd_tt, ipl, speed3=use_speed3_state
        )
        if use_speed3_velocity and not use_speed3_state:
            speed3_state = lunar.calc_interpolated_apse_state(jd_tt, ipl, speed3=True)
            dlon, dlat, ddist = speed3_state[3:]
        # The independent interpolated-apse model is defined on the true
        # ecliptic of date (nutation included). Under NONUT, J2000, or the
        # SIDEREAL+EQUATORIAL mean-frame convention, subtract dpsi before the
        # requested frame conversion.
        _sid_eq = is_sidereal and bool(calc_iflag & FLG_EQUATORIAL)
        _no_nut = bool(calc_iflag & (FLG_NONUT | FLG_J2000)) or _sid_eq
        if _no_nut:
            from .cache import get_cached_nutation

            dpsi_rad, _ = get_cached_nutation(jd_tt)
            lon = (lon - math.degrees(dpsi_rad)) % 360.0
        if calc_iflag & FLG_SPEED:
            # The interpolated curves are true-ecliptic; under NONUT
            # (mean-ecliptic output) remove the nutation rate.
            if _no_nut:
                dlon -= _nutation_rate_deg_per_day(jd_tt, 0.001)
        else:
            dlon = dlat = ddist = 0.0
        # FLG_J2000 honored uniformly under SIDEREAL: precess to J2000,
        # then subtract the mean ayanamsha (MEAN_NODE/MEAN_APOG pipeline).
        _sid_j2000_ecl = (
            is_sidereal and bool(iflag & FLG_J2000) and not bool(iflag & FLG_EQUATORIAL)
        )
        result = (lon, lat, dist, dlon, dlat, ddist)
        if _sid_j2000_ecl:
            result = _maybe_equatorial_convert(result, jd_tt, iflag)
            lon, lat, dist, dlon, dlat, ddist = result
            lon, dlon = _apply_sidereal_correction(lon, dlon, t.ut1, iflag)
            result = (lon, lat, dist, dlon, dlat, ddist)
        else:
            if is_sidereal and not (iflag & FLG_EQUATORIAL):
                lon, dlon = _apply_sidereal_correction(lon, dlon, t.ut1, iflag)
                result = (lon, lat, dist, dlon, dlat, ddist)
            result = _maybe_equatorial_convert(result, jd_tt, iflag)
        return _to_native_floats(result), iflag

    # Handle Uranian planets (Hamburg School hypothetical bodies, IDs 40-47)
    if CUPIDO <= ipl <= POSEIDON:
        from . import hypothetical
        from skyfield.framelib import ecliptic_J2000_frame
        from .astrometry import apply_aberration_to_position

        jd_tt = t.tt
        # Center priority TOPOCTR > BARYCTR > HELCTR (see _resolve_center_flags):
        # a barycentric request suppresses the heliocentric branch, and a
        # topocentric one suppresses both so the geocentric+parallax path runs.
        is_topo = bool(iflag & FLG_TOPOCTR)
        is_bary = bool(iflag & FLG_BARYCTR) and not is_topo
        is_helio = bool(iflag & FLG_HELCTR) and not is_bary and not is_topo
        is_j2000 = bool(iflag & FLG_J2000)
        is_sidereal = bool(iflag & FLG_SIDEREAL)
        # The geocentric apparent place carries stellar aberration unless the
        # caller asked for the astrometric (FLG_NOABERR) or true (FLG_TRUEPOS)
        # place. Barycentric and heliocentric places are astrometric by
        # convention (their retflag echoes NOABERR|NOGDEFL), so aberration is
        # applied only on the geocentric branch below.
        apply_aberration = not (iflag & FLG_TRUEPOS) and not (iflag & FLG_NOABERR)

        # Position sourcing. With an active LEB reader the vector-resolver
        # Earth evaluation — and, when the serving file byte-matches its
        # pinned manifest SHA-256, the body propagation itself — come from the
        # reader. Every frame transform and finite-difference realization
        # below is shared by both sources, so the public speed semantics do
        # not depend on where the raw positions come from; the residual source
        # delta is the LEB approximation class (~1e-12 deg for the body,
        # <0.001" for Earth). Out-of-range dates fall back per sample.
        from .state import get_leb_reader, leb_fictitious_source_trusted

        # In mode="leb" get_leb_reader() raises when no global LEB resolves.
        # A context calculation may carry its own LEB while the global one is
        # absent; treat that as "no reader here" and fall back to the runtime
        # model exactly as the pre-companion branch did, instead of surfacing
        # the RuntimeError to the public caller.
        try:
            _reader = get_leb_reader()
        except RuntimeError:
            _reader = None

        if _reader is not None and leb_fictitious_source_trusted(_reader, ipl):
            # The body position is served from the manifest-pinned companion,
            # so report LEB as the dispatch backend (mirrors the SPK branch).
            _mark_dispatch_source("LEB")

            def _uranian_pos(jd: float) -> Tuple[float, float, float]:
                try:
                    (u_lon, u_lat, u_dist), _ = _reader.eval_body(ipl, jd)
                except LEBCorruptionError:
                    raise
                except (KeyError, ValueError):
                    return hypothetical._calc_uranian_planet_raw(ipl, jd)
                return u_lon, u_lat, u_dist

        else:

            def _uranian_pos(jd: float) -> Tuple[float, float, float]:
                return hypothetical._calc_uranian_planet_raw(ipl, jd)

        def _uranian_helio_state(
            jd: float,
        ) -> Tuple[float, float, float, float, float, float]:
            # Same centered 1-day stencil and 0/360 unwrap as
            # hypothetical.calc_uranian_planet: with the pure-model source the
            # returned state is arithmetic-identical to that function.
            u_lon, u_lat, u_dist = _uranian_pos(jd)
            dt_step = 1.0
            pos_prev = _uranian_pos(jd - dt_step)
            pos_next = _uranian_pos(jd + dt_step)
            u_dlon = (pos_next[0] - pos_prev[0]) / (2.0 * dt_step)
            # 0/360 wrap: valid only because these bodies are slow
            # (<<90 deg/day at dt_step=1); see calc_uranian_planet.
            if u_dlon > 180.0 / (2.0 * dt_step):
                u_dlon -= 360.0 / (2.0 * dt_step)
            elif u_dlon < -180.0 / (2.0 * dt_step):
                u_dlon += 360.0 / (2.0 * dt_step)
            u_dlat = (pos_next[1] - pos_prev[1]) / (2.0 * dt_step)
            u_ddist = (pos_next[2] - pos_prev[2]) / (2.0 * dt_step)
            return (u_lon, u_lat, u_dist, u_dlon, u_dlat, u_ddist)

        def _earth_helio_ecl_j2000_vector(jd: float) -> Tuple[float, float, float]:
            ts_i = get_timescale()
            t_i = ts_i.tt_jd(jd)
            earth_h = planets["sun"].at(t_i).observe(planets["earth"])
            exyz = earth_h.frame_xyz(ecliptic_J2000_frame).au
            return float(exyz[0]), float(exyz[1]), float(exyz[2])

        # The exact ICRS -> J2000-ecliptic matrix of the Skyfield frame used by
        # the fallback path. It is constant for this inertial frame, so it is
        # evaluated once (at the request time t) and reused for every sample
        # date; every Sun/Earth sub-helper below rotates through this identical
        # matrix so the LEB and Skyfield sources agree to the last bit.
        _ecl_m = ecliptic_J2000_frame.rotation_at(t)

        def _rot_ecl_j2000(
            vx: float, vy: float, vz: float
        ) -> Tuple[float, float, float]:
            return (
                float(_ecl_m[0, 0] * vx + _ecl_m[0, 1] * vy + _ecl_m[0, 2] * vz),
                float(_ecl_m[1, 0] * vx + _ecl_m[1, 1] * vy + _ecl_m[1, 2] * vz),
                float(_ecl_m[2, 0] * vx + _ecl_m[2, 1] * vy + _ecl_m[2, 2] * vz),
            )

        def _sun_bary_ecl_j2000_vector(jd: float) -> Tuple[float, float, float]:
            # The Sun's position relative to the solar-system barycentre, in
            # ecliptic J2000 (the helio -> bary shift, ~0.008 AU).
            ts_i = get_timescale()
            t_i = ts_i.tt_jd(jd)
            sxyz = planets["sun"].at(t_i).frame_xyz(ecliptic_J2000_frame).au
            return float(sxyz[0]), float(sxyz[1]), float(sxyz[2])

        def _earth_bary_vel_ecl_j2000_vector(
            jd: float,
        ) -> Tuple[float, float, float]:
            # Earth's barycentric (SSB) velocity in ecliptic J2000 (AU/day),
            # the observer velocity for stellar aberration.
            ts_i = get_timescale()
            t_i = ts_i.tt_jd(jd)
            v = planets["earth"].at(t_i).velocity.au_per_d
            return _rot_ecl_j2000(float(v[0]), float(v[1]), float(v[2]))

        if _reader is not None and _reader.has_body(SUN) and _reader.has_body(EARTH):
            from .fast_calc import C_LIGHT_AU_DAY

            def _earth_helio_ecl_j2000(jd: float) -> Tuple[float, float, float]:
                try:
                    (ex, ey, ez), _ = _reader.eval_body(EARTH, jd)
                    (sx, sy, sz), _ = _reader.eval_body(SUN, jd)
                except LEBCorruptionError:
                    raise
                except (KeyError, ValueError):
                    return _earth_helio_ecl_j2000_vector(jd)
                # Astrometric Earth-from-Sun: iterate the light time on the
                # Earth sample (Sun frozen at jd), replicating the Skyfield
                # observe() semantics of the fallback path above.
                rx, ry, rz = ex - sx, ey - sy, ez - sz
                lt = math.sqrt(rx * rx + ry * ry + rz * rz) / C_LIGHT_AU_DAY
                for _ in range(10):
                    try:
                        (ex, ey, ez), _ = _reader.eval_body(EARTH, jd - lt)
                    except LEBCorruptionError:
                        raise
                    except (KeyError, ValueError):
                        # A composite reader closed mid-calculation (e.g. a
                        # concurrent set_leb_file/close) raises KeyError once
                        # its body map is cleared; use the mode-aware vector
                        # resolver rather than crossing the backend boundary.
                        return _earth_helio_ecl_j2000_vector(jd)
                    rx, ry, rz = ex - sx, ey - sy, ez - sz
                    lt_next = math.sqrt(rx * rx + ry * ry + rz * rz) / C_LIGHT_AU_DAY
                    converged = abs(lt_next - lt) < 1e-14
                    lt = lt_next
                    if converged:
                        break
                return _rot_ecl_j2000(rx, ry, rz)

            def _sun_bary_ecl_j2000(jd: float) -> Tuple[float, float, float]:
                # SSB-relative Sun position from the LEB (barycentric ICRS),
                # rotated to ecliptic J2000 — the helio -> bary shift.
                try:
                    (sx, sy, sz), _ = _reader.eval_body(SUN, jd)
                except LEBCorruptionError:
                    raise
                except (KeyError, ValueError):
                    return _sun_bary_ecl_j2000_vector(jd)
                return _rot_ecl_j2000(sx, sy, sz)

            def _earth_bary_vel_ecl_j2000(
                jd: float,
            ) -> Tuple[float, float, float]:
                # Earth's barycentric velocity from the LEB (analytical
                # Chebyshev derivative, AU/day), rotated to ecliptic J2000.
                try:
                    _, (evx, evy, evz) = _reader.eval_body(EARTH, jd)
                except LEBCorruptionError:
                    raise
                except (KeyError, ValueError):
                    return _earth_bary_vel_ecl_j2000_vector(jd)
                return _rot_ecl_j2000(evx, evy, evz)

        else:
            _earth_helio_ecl_j2000 = _earth_helio_ecl_j2000_vector
            _sun_bary_ecl_j2000 = _sun_bary_ecl_j2000_vector
            _earth_bary_vel_ecl_j2000 = _earth_bary_vel_ecl_j2000_vector

        def _uranian_bary_pos(jd: float) -> Tuple[float, float, float]:
            # Barycentric ecliptic-J2000 place: the body's heliocentric
            # position in ecliptic-J2000 Cartesian, shifted by the Sun's own
            # barycentric offset. bary = helio + Sun_bary. (The prior code left
            # BARYCTR to fall through to the geocentric branch, which returns
            # helio - Earth_from_Sun — a ~1 AU error for these distant bodies.)
            u_lon, u_lat, u_dist = _uranian_pos(jd)
            lon_r = math.radians(u_lon)
            lat_r = math.radians(u_lat)
            cl = math.cos(lat_r)
            xh = u_dist * cl * math.cos(lon_r)
            yh = u_dist * cl * math.sin(lon_r)
            zh = u_dist * math.sin(lat_r)
            sx, sy, sz = _sun_bary_ecl_j2000(jd)
            xb, yb, zb = xh + sx, yh + sy, zh + sz
            rb = math.sqrt(xb * xb + yb * yb + zb * zb)
            lon_b = math.degrees(math.atan2(yb, xb)) % 360.0
            sin_lat = max(-1.0, min(1.0, zb / rb)) if rb > 0 else 0.0
            lat_b = math.degrees(math.asin(sin_lat))
            return lon_b, lat_b, rb

        def _uranian_bary_state(
            jd: float,
        ) -> Tuple[float, float, float, float, float, float]:
            # Same centered 1-day stencil / 0-360 unwrap as the heliocentric
            # state; the Sun's slow barycentric drift enters the speed through
            # the finite difference.
            lon0, lat0, dist0 = _uranian_bary_pos(jd)
            dt_step = 1.0
            prev = _uranian_bary_pos(jd - dt_step)
            nxt = _uranian_bary_pos(jd + dt_step)
            dlon = (nxt[0] - prev[0]) / (2.0 * dt_step)
            if dlon > 180.0 / (2.0 * dt_step):
                dlon -= 360.0 / (2.0 * dt_step)
            elif dlon < -180.0 / (2.0 * dt_step):
                dlon += 360.0 / (2.0 * dt_step)
            dlat = (nxt[1] - prev[1]) / (2.0 * dt_step)
            ddist = (nxt[2] - prev[2]) / (2.0 * dt_step)
            return (lon0, lat0, dist0, dlon, dlat, ddist)

        if is_helio or is_bary:
            # Heliocentric, or (with the Sun's barycentric offset added)
            # barycentric, J2000 ecliptic state (position + centered velocity).
            pos = _uranian_bary_state(jd_tt) if is_bary else _uranian_helio_state(jd_tt)
            lon, lat, dist = pos[0], pos[1], pos[2]
            dlon, dlat, ddist = pos[3], pos[4], pos[5]
            # Zero the speed slots when FLG_SPEED is absent so callers never
            # receive unrequested velocity data (the reference API returns
            # zeros); the Δψ-rate and sidereal speed terms below are
            # themselves gated on FLG_SPEED, so the zeros survive.
            if not (iflag & FLG_SPEED):
                dlon = dlat = ddist = 0.0
            if not is_j2000:
                # Precess J2000 -> ecliptic of date. Carry the velocity
                # through the same precession (do not leave it in the J2000
                # frame) so the of-date speed is the time derivative of the
                # of-date position and picks up the general-precession rate.
                if iflag & FLG_SPEED:
                    lon, lat, dlon, dlat = _precess_ecliptic_state(
                        lon, lat, dlon, dlat, jd_tt, to_j2000=False
                    )
                else:
                    from .astrometry import _precess_ecliptic

                    lon, lat = _precess_ecliptic(lon, lat, 2451545.0, jd_tt)
            # Add nutation in longitude (Δψ) to reach the true ecliptic of date,
            # so the sidereal ayanamsha (mean + Δψ) subtracted below lands on the
            # intrinsic sidereal longitude. No-op for J2000, NONUT, and
            # SIDEREAL+EQUATORIAL frames.
            lon, dlon = _add_of_date_nutation(lon, dlon, jd_tt, iflag)
            # Apply sidereal correction if requested (not for equatorial output)
            if is_sidereal and not (iflag & FLG_EQUATORIAL):
                lon, dlon = _apply_sidereal_correction(lon, dlon, t.ut1, iflag)
            result = (lon, lat, dist, dlon, dlat, ddist)
            # Coordinates are already J2000 ecliptic when is_j2000 (precession
            # handled above): skip the precession inside but keep the J2000
            # obliquity for EQ+J2000 (consistent frame, like all other bodies).
            result = _maybe_equatorial_convert(result, jd_tt, iflag, already_j2000=True)
            return _to_native_floats(result), iflag

        # Geocentric: convert heliocentric Keplerian orbit to geocentric
        def _get_uranian_geo_j2000(jd):
            h = _uranian_pos(jd)
            lon_r = math.radians(h[0])
            lat_r = math.radians(h[1])
            cl = math.cos(lat_r)
            xh = h[2] * cl * math.cos(lon_r)
            yh = h[2] * cl * math.sin(lon_r)
            zh = h[2] * math.sin(lat_r)

            ex, ey, ez = _earth_helio_ecl_j2000(jd)
            xg = xh - ex
            yg = yh - ey
            zg = zh - ez
            # Stellar (annual) aberration for the geocentric apparent place,
            # using Earth's barycentric velocity in ecliptic J2000 — the same
            # IAU convention and helper as the real-planet and SPK pipelines.
            # Folding it into this sample (rather than only the central place)
            # keeps the FLG_SPEED stencil consistent — the speed picks up the
            # aberration rate. Skipped for the astrometric (FLG_NOABERR) or true
            # (FLG_TRUEPOS) place. Applied in the J2000 ecliptic frame before the
            # of-date precession, exactly like the SPK type-21 apparent path.
            if apply_aberration:
                evx, evy, evz = _earth_bary_vel_ecl_j2000(jd)
                xg, yg, zg = apply_aberration_to_position((xg, yg, zg), (evx, evy, evz))
            rg = math.sqrt(xg * xg + yg * yg + zg * zg)
            lon_g = math.degrees(math.atan2(yg, xg)) % 360.0
            sin_lat = max(-1.0, min(1.0, zg / rg)) if rg > 0 else 0.0
            lat_g = math.degrees(math.asin(sin_lat))
            return lon_g, lat_g, rg

        lon, lat, dist = _get_uranian_geo_j2000(jd_tt)

        # The central-difference stencil runs only when speeds are requested:
        # without FLG_SPEED the reference API returns zeroed speed slots, and
        # skipping the stencil saves two extra Earth-position evaluations.
        if iflag & FLG_SPEED:
            dt_v = 1.0
            prev = _get_uranian_geo_j2000(jd_tt - dt_v)
            nxt = _get_uranian_geo_j2000(jd_tt + dt_v)
            dlon = (nxt[0] - prev[0]) / (2.0 * dt_v)
            # Unwrap a 0/360 boundary crossing. The longitudes are in [0, 360),
            # so a crossing is a ~360 deg jump in the raw difference; after
            # dividing by 2*dt_v the threshold and correction carry the same
            # 1/(2*dt_v) factor. Valid only because these bodies are slow
            # (<<90 deg/day at dt_v=1): a real speed above 180/(2*dt_v) cannot
            # occur, so it can only be a wrap. Do NOT reuse this for fast
            # bodies -- it would misread real motion.
            if dlon > 180.0 / (2.0 * dt_v):
                dlon -= 360.0 / (2.0 * dt_v)
            elif dlon < -180.0 / (2.0 * dt_v):
                dlon += 360.0 / (2.0 * dt_v)
            dlat = (nxt[1] - prev[1]) / (2.0 * dt_v)
            ddist = (nxt[2] - prev[2]) / (2.0 * dt_v)
        else:
            dlon = dlat = ddist = 0.0

        if not is_j2000:
            # Precess J2000 -> ecliptic of date, carrying the velocity so the
            # of-date speed is the derivative of the of-date position (see the
            # heliocentric branch above).
            if iflag & FLG_SPEED:
                lon, lat, dlon, dlat = _precess_ecliptic_state(
                    lon, lat, dlon, dlat, jd_tt, to_j2000=False
                )
            else:
                from .astrometry import _precess_ecliptic

                lon, lat = _precess_ecliptic(lon, lat, 2451545.0, jd_tt)

        # Add nutation in longitude (Δψ) so the of-date place is on the true
        # ecliptic (matching the reference API's uniform Δψ), before topocentric
        # parallax and the sidereal ayanamsha rotation. No-op for J2000 / NONUT /
        # SIDEREAL+EQUATORIAL frames.
        lon, dlon = _add_of_date_nutation(lon, dlon, jd_tt, iflag)

        # Topocentric parallax first (tropical of-date / J2000), so the
        # sidereal ayanamsha rotation below acts on the topocentric longitude.
        result = (lon, lat, dist, dlon, dlat, ddist)
        result = _apply_hypothetical_topocentric(
            result, jd_tt, iflag, is_j2000, planets
        )
        lon, lat, dist, dlon, dlat, ddist = result

        if is_sidereal and not (iflag & FLG_EQUATORIAL):
            lon, dlon = _apply_sidereal_correction(lon, dlon, t.ut1, iflag)

        result = (lon, lat, dist, dlon, dlat, ddist)
        # Already J2000 ecliptic when is_j2000 (see the heliocentric branch):
        # keep the J2000 obliquity for EQ+J2000, skip the internal precession.
        result = _maybe_equatorial_convert(result, jd_tt, iflag, already_j2000=True)
        return _to_native_floats(result), iflag

    # Handle Transpluto (Isis) — ISIS = 48
    if ipl == ISIS:
        from . import hypothetical
        from skyfield.framelib import ecliptic_J2000_frame

        jd_tt = t.tt
        is_helio = bool(iflag & FLG_HELCTR)
        is_j2000 = bool(iflag & FLG_J2000)
        is_sidereal = bool(iflag & FLG_SIDEREAL)

        # Same apparent-place contract as the other predicted planets: the
        # body is taken at the retarded epoch and, geocentrically, the
        # direction is corrected for the observer's barycentric velocity.
        # FLG_TRUEPOS asks for the geometric place and drops both terms;
        # FLG_NOABERR drops the aberration alone. Speed of light in AU/day,
        # the shared NOVAS/Skyfield C_AUDAY realization.
        _C_AU_DAY = 173.1446326846693
        _tp_truepos = bool(iflag & FLG_TRUEPOS)
        _tp_noaberr = bool(iflag & (FLG_TRUEPOS | FLG_NOABERR))

        def _transpluto_light_time(jd: float) -> float:
            """Sun-to-body light time in days (0 for a geometric request)."""
            if _tp_truepos:
                return 0.0
            lt = 0.0
            for _ in range(3):
                lt = float(hypothetical.calc_transpluto(jd - lt)[2]) / _C_AU_DAY
            return lt

        def _transpluto_helio_xyz(jd: float) -> Tuple[float, float, float]:
            """Heliocentric J2000-ecliptic cartesian position (AU)."""
            h = hypothetical.calc_transpluto(jd)
            lon_r = math.radians(h[0])
            lat_r = math.radians(h[1])
            cl = math.cos(lat_r)
            return (
                h[2] * cl * math.cos(lon_r),
                h[2] * cl * math.sin(lon_r),
                h[2] * math.sin(lat_r),
            )

        if (iflag & FLG_BARYCTR) and not is_helio:
            # Barycentric place: the model is heliocentric, so add the Sun's
            # barycentric vector. Without this branch the request fell through
            # to the geocentric path and returned geocentric coordinates while
            # advertising BARYCTR.
            from .astrometry import _precess_ecliptic as _pe

            def _transpluto_bary_place(jd: float) -> Tuple[float, float, float]:
                """Barycentric (lon, lat, dist) in the requested output frame."""
                _x, _y, _z = _transpluto_helio_xyz(jd - _transpluto_light_time(jd))
                _sv = (
                    planets["sun"]
                    .at(get_timescale().tt_jd(jd))
                    .frame_xyz(ecliptic_J2000_frame)
                    .au
                )
                _x += float(_sv[0])
                _y += float(_sv[1])
                _z += float(_sv[2])
                _r = math.sqrt(_x * _x + _y * _y + _z * _z)
                _lo = math.degrees(math.atan2(_y, _x)) % 360.0
                _la = math.degrees(
                    math.asin(max(-1.0, min(1.0, _z / _r)) if _r > 0 else 0.0)
                )
                if not is_j2000:
                    _lo, _la = _pe(_lo, _la, 2451545.0, jd)
                return _lo, _la, _r

            lon, lat, dist = _transpluto_bary_place(jd_tt)
            dlon = dlat = ddist = 0.0
            if iflag & (FLG_SPEED | FLG_SPEED3):
                # Central difference of the SAME reported quantity, so the
                # speed slots integrate back to the position slots. Leaving
                # them hard-zeroed advertised speeds the request never
                # received while the position plainly moves with time.
                _h = 0.01
                _p0 = _transpluto_bary_place(jd_tt - _h)
                _p1 = _transpluto_bary_place(jd_tt + _h)
                dlon = ((_p1[0] - _p0[0] + 180.0) % 360.0 - 180.0) / (2.0 * _h)
                dlat = (_p1[1] - _p0[1]) / (2.0 * _h)
                ddist = (_p1[2] - _p0[2]) / (2.0 * _h)
            lon, dlon = _add_of_date_nutation(lon, dlon, jd_tt, iflag)
            if is_sidereal and not (iflag & FLG_EQUATORIAL):
                lon, dlon = _apply_sidereal_correction(lon, dlon, t.ut1, iflag)
            result = (lon, lat, dist, dlon, dlat, ddist)
            result = _maybe_equatorial_convert(result, jd_tt, iflag, already_j2000=True)
            return _to_native_floats(result), iflag

        if is_helio:
            pos = hypothetical.calc_transpluto(jd_tt - _transpluto_light_time(jd_tt))
            lon, lat, dist = pos[0], pos[1], pos[2]
            dlon, dlat, ddist = pos[3], pos[4], pos[5]
            # Zero the speed slots when FLG_SPEED is absent (reference API
            # contract — see the Uranian heliocentric branch above).
            if not (iflag & FLG_SPEED):
                dlon = dlat = ddist = 0.0
            if not is_j2000:
                # Precess J2000 -> ecliptic of date, carrying the velocity so
                # the of-date speed is the derivative of the of-date position
                # (see the Uranian heliocentric branch above).
                if iflag & FLG_SPEED:
                    lon, lat, dlon, dlat = _precess_ecliptic_state(
                        lon, lat, dlon, dlat, jd_tt, to_j2000=False
                    )
                else:
                    from .astrometry import _precess_ecliptic

                    lon, lat = _precess_ecliptic(lon, lat, 2451545.0, jd_tt)
            # Add nutation in longitude (Δψ) to reach the true ecliptic of date
            # (see the Uranian heliocentric branch). No-op for J2000 / NONUT /
            # SIDEREAL+EQUATORIAL frames.
            lon, dlon = _add_of_date_nutation(lon, dlon, jd_tt, iflag)
            # Apply sidereal correction if requested (not for equatorial output)
            if is_sidereal and not (iflag & FLG_EQUATORIAL):
                lon, dlon = _apply_sidereal_correction(lon, dlon, t.ut1, iflag)
            result = (lon, lat, dist, dlon, dlat, ddist)
            # Already J2000 ecliptic when is_j2000 (see the Uranian branch):
            # keep the J2000 obliquity for EQ+J2000, skip the internal precession.
            result = _maybe_equatorial_convert(result, jd_tt, iflag, already_j2000=True)
            return _to_native_floats(result), iflag

        # Geocentric conversion
        def _get_transpluto_geo_j2000(jd):
            """Geocentric J2000-ecliptic apparent position for Transpluto.

            Standard apparent-place reduction (Explanatory Supplement to the
            Astronomical Almanac, 3rd ed. 2013, ch. 7): retarded target plus
            annual aberration, both suppressed for a geometric request. The
            purely geometric reduction used before left the place short by
            the ~20" aberration and made FLG_TRUEPOS a no-op here.
            """
            ts_i = get_timescale()
            t_i = ts_i.tt_jd(jd)
            earth_h = planets["sun"].at(t_i).observe(planets["earth"])
            exyz = earth_h.frame_xyz(ecliptic_J2000_frame).au
            ex, ey, ez = float(exyz[0]), float(exyz[1]), float(exyz[2])

            lt = 0.0
            xg = yg = zg = 0.0
            for _ in range(3):
                xh, yh, zh = _transpluto_helio_xyz(jd - lt)
                xg, yg, zg = xh - ex, yh - ey, zh - ez
                if _tp_truepos:
                    break
                lt = math.sqrt(xg * xg + yg * yg + zg * zg) / _C_AU_DAY

            rg = math.sqrt(xg * xg + yg * yg + zg * zg)
            if not _tp_noaberr and rg > 0.0:
                evel = (
                    planets["earth"]
                    .at(t_i)
                    .frame_xyz_and_velocity(ecliptic_J2000_frame)[1]
                    .au_per_d
                )
                xg += float(evel[0]) * rg / _C_AU_DAY
                yg += float(evel[1]) * rg / _C_AU_DAY
                zg += float(evel[2]) * rg / _C_AU_DAY
                rg = math.sqrt(xg * xg + yg * yg + zg * zg)
            lon_g = math.degrees(math.atan2(yg, xg)) % 360.0
            sin_lat = max(-1.0, min(1.0, zg / rg)) if rg > 0 else 0.0
            lat_g = math.degrees(math.asin(sin_lat))
            return lon_g, lat_g, rg

        lon, lat, dist = _get_transpluto_geo_j2000(jd_tt)

        # Stencil only when speeds are requested (reference API returns zeroed
        # speed slots without FLG_SPEED — see the Uranian geocentric branch).
        if iflag & FLG_SPEED:
            dt_v = 1.0
            prev = _get_transpluto_geo_j2000(jd_tt - dt_v)
            nxt = _get_transpluto_geo_j2000(jd_tt + dt_v)
            dlon = (nxt[0] - prev[0]) / (2.0 * dt_v)
            # Unwrap a 0/360 boundary crossing (see the Uranian path above):
            # the threshold and correction must carry the 1/(2*dt_v) factor.
            # Valid only for slow bodies (<<90 deg/day at dt_v=1); do NOT
            # reuse for fast ones.
            if dlon > 180.0 / (2.0 * dt_v):
                dlon -= 360.0 / (2.0 * dt_v)
            elif dlon < -180.0 / (2.0 * dt_v):
                dlon += 360.0 / (2.0 * dt_v)
            dlat = (nxt[1] - prev[1]) / (2.0 * dt_v)
            ddist = (nxt[2] - prev[2]) / (2.0 * dt_v)
        else:
            dlon = dlat = ddist = 0.0

        if not is_j2000:
            # Precess J2000 -> ecliptic of date, carrying the velocity so the
            # of-date speed is the derivative of the of-date position (see the
            # Uranian branch above).
            if iflag & FLG_SPEED:
                lon, lat, dlon, dlat = _precess_ecliptic_state(
                    lon, lat, dlon, dlat, jd_tt, to_j2000=False
                )
            else:
                from .astrometry import _precess_ecliptic

                lon, lat = _precess_ecliptic(lon, lat, 2451545.0, jd_tt)

        # Add nutation in longitude (Δψ) so the of-date place is on the true
        # ecliptic (see the Uranian geocentric branch), before topocentric
        # parallax and the sidereal ayanamsha rotation. No-op for J2000 / NONUT /
        # SIDEREAL+EQUATORIAL frames.
        lon, dlon = _add_of_date_nutation(lon, dlon, jd_tt, iflag)

        # Topocentric parallax first (tropical of-date / J2000), so the
        # sidereal ayanamsha rotation below acts on the topocentric longitude.
        result = (lon, lat, dist, dlon, dlat, ddist)
        result = _apply_hypothetical_topocentric(
            result, jd_tt, iflag, bool(iflag & FLG_J2000), planets
        )
        lon, lat, dist, dlon, dlat, ddist = result

        # Apply sidereal correction if requested (not for equatorial output)
        if is_sidereal and not (iflag & FLG_EQUATORIAL):
            lon, dlon = _apply_sidereal_correction(lon, dlon, t.ut1, iflag)

        result = (lon, lat, dist, dlon, dlat, ddist)
        # Already J2000 ecliptic when is_j2000 (see the Uranian branch):
        # keep the J2000 obliquity for EQ+J2000, skip the internal precession.
        result = _maybe_equatorial_convert(result, jd_tt, iflag, already_j2000=True)
        return _to_native_floats(result), iflag

    # Route the non-Uranian historical compatibility IDs through the local
    # hypothetical-body layer.  Only the entries admitted by that layer's
    # reviewed provenance registry can return a state: Selena (56) is the
    # source-backed geocentric symbolic construction, and IDs 48, 50--55 and
    # 57 carry heliocentric predicted-planet models transcribed from their
    # published sources.  IDs 49 and 58 remain in this dispatcher so every
    # public constant reaches one deterministic error boundary; their
    # calculation functions raise UnknownBodyError before any unavailable
    # element could be used.
    #
    # The heliocentric/geocentric classification below is therefore also a
    # coordinate-routing declaration for future provenance-complete models,
    # not evidence that the currently unsupported names carry numerical data.
    _FICT_HELIO_IDS = (
        NIBIRU,
        HARRINGTON,
        NEPTUNE_LEVERRIER,
        NEPTUNE_ADAMS,
        PLUTO_LOWELL,
        PLUTO_PICKERING,
    )
    if ipl in (WHITE_MOON, VULCAN, PROSERPINA, WALDEMATH) or ipl in _FICT_HELIO_IDS:
        from . import hypothetical

        jd_tt = t.tt
        is_sidereal = bool(iflag & FLG_SIDEREAL)

        if ipl in (VULCAN, PROSERPINA) or ipl in _FICT_HELIO_IDS:
            from skyfield.framelib import ecliptic_J2000_frame
            from .astrometry import _precess_ecliptic

            # Speed of light in AU/day, NOVAS/Skyfield C_AUDAY realization
            # (c exact from SI with the IAU 1976 au), identical to the value
            # used by astrometry.py and fast_calc.py so every light-time
            # iteration in the library shares one constant.
            _C_AU_DAY = 173.1446326846693
            if ipl == VULCAN:
                _calc_helio = hypothetical.calc_vulcan
            elif ipl == PROSERPINA:
                _calc_helio = hypothetical.calc_proserpina
            else:

                def _calc_helio(jd, _ipl=ipl):
                    return hypothetical.calc_fictitious_position(_ipl, jd)

            # FLG_TRUEPOS asks for the TRUE (geometric) position: the body
            # where it actually is at the request epoch, not where it is seen
            # from. For these analytically propagated bodies the only
            # light-path term the model carries is the light-time
            # retardation, so TRUEPOS simply evaluates the orbit at jd
            # instead of at jd - light-time. Measured against the reference,
            # FLG_NOABERR and FLG_NOGDEFL leave this path unchanged (no
            # separate observer-motion or deflection term exists here), and
            # TRUEPOS changes the geocentric longitude by several arcseconds
            # (e.g. ~18" for the outer predicted planets, whose light-time is
            # about half a day) and the heliocentric one by the Sun-to-body
            # light-time.
            _truepos = bool(iflag & FLG_TRUEPOS)

            def _helio_light_time_days(jd: float) -> float:
                """Sun-to-body light time in days at ``jd`` (0 under TRUEPOS).

                Iterated like the geocentric case: the heliocentric distance
                is itself evaluated at the retarded epoch. Two passes are
                ample — these orbits move well under a degree per day.
                """
                if _truepos:
                    return 0.0
                lt = 0.0
                for _ in range(3):
                    lt = float(_calc_helio(jd - lt)[2]) / _C_AU_DAY
                return lt

            def _helio_xyz_j2000(jd: float) -> Tuple[float, float, float]:
                """Heliocentric J2000-ecliptic cartesian position (AU)."""
                h = _calc_helio(jd)
                if ipl in (VULCAN, PROSERPINA):
                    # These element sets use the equinox of date; rotate
                    # to J2000 so the Earth subtraction is
                    # frame-consistent.
                    lon_j, lat_j = _precess_ecliptic(h[0], h[1], jd, 2451545.0)
                else:
                    # The predicted-planet propagation already returns
                    # J2000 ecliptic coordinates.
                    lon_j, lat_j = h[0], h[1]
                lon_r = math.radians(lon_j)
                lat_r = math.radians(lat_j)
                cl = math.cos(lat_r)
                return (
                    h[2] * cl * math.cos(lon_r),
                    h[2] * cl * math.sin(lon_r),
                    h[2] * math.sin(lat_r),
                )

            # Annual aberration is suppressed by FLG_TRUEPOS (geometric
            # place) and by FLG_NOABERR (its documented purpose). The
            # library applies the same rule to every other body, so these
            # analytic models follow the shared contract rather than
            # carrying a second one.
            _noaberr = bool(iflag & (FLG_TRUEPOS | FLG_NOABERR))
            _nogdefl = bool(iflag & (FLG_TRUEPOS | FLG_NOGDEFL))

            def _geo_of_date(jd: float) -> Tuple[float, float, float]:
                """Geocentric mean-ecliptic-of-date apparent position.

                Standard apparent-place reduction (Explanatory Supplement to
                the Astronomical Almanac, 3rd ed. 2013, ch. 7): the body is
                taken at the retarded epoch (light time) and the direction is
                corrected for the observer's barycentric velocity (annual
                aberration, up to ~20.5"). Both terms are suppressed for a
                geometric request; without the aberration term the geocentric
                place was short by that ~20" and FLG_TRUEPOS had almost
                nothing to remove.
                """
                ts_i = get_timescale()
                t_i = ts_i.tt_jd(jd)
                earth_h = planets["sun"].at(t_i).observe(planets["earth"])
                exyz = earth_h.frame_xyz(ecliptic_J2000_frame).au
                ex, ey, ez = float(exyz[0]), float(exyz[1]), float(exyz[2])

                # Light-time iteration: target retarded, observer fixed.
                # Under FLG_TRUEPOS the retardation is suppressed and the
                # body is taken at the request epoch (geometric position).
                lt = 0.0
                xg = yg = zg = 0.0
                for _ in range(3):
                    bx, by, bz = _helio_xyz_j2000(jd - lt)
                    xg, yg, zg = bx - ex, by - ey, bz - ez
                    if _truepos:
                        break
                    lt = math.sqrt(xg * xg + yg * yg + zg * zg) / _C_AU_DAY

                rg = math.sqrt(xg * xg + yg * yg + zg * zg)
                if not _nogdefl and rg > 0.0:
                    # Solar gravitational light deflection, from the shared
                    # PPN routine the planetary paths use. Omitting it left
                    # FLG_NOGDEFL with nothing to switch off on this whole
                    # family, while the reference moves by up to 0.16" for a
                    # body near the Sun's line of sight.
                    from .fast_calc import _apply_gravitational_deflection

                    _src = _SkyfieldDeflectorSource(planets)
                    xg, yg, zg = _apply_gravitational_deflection(
                        (xg, yg, zg),
                        (ex, ey, ez),
                        jd,
                        rg / _C_AU_DAY,
                        _src,
                    )
                    rg = math.sqrt(xg * xg + yg * yg + zg * zg)
                if not _noaberr and rg > 0.0:
                    # First-order annual aberration: shift the geometric
                    # direction by the observer's barycentric velocity over c
                    # (v/c ~ 1e-4, so higher-order terms stay below a
                    # milliarcsecond and well under this path's model error).
                    evel = (
                        planets["earth"]
                        .at(t_i)
                        .frame_xyz_and_velocity(ecliptic_J2000_frame)[1]
                        .au_per_d
                    )
                    xg += float(evel[0]) * rg / _C_AU_DAY
                    yg += float(evel[1]) * rg / _C_AU_DAY
                    zg += float(evel[2]) * rg / _C_AU_DAY
                    rg = math.sqrt(xg * xg + yg * yg + zg * zg)
                lon_g = math.degrees(math.atan2(yg, xg)) % 360.0
                sin_lat = max(-1.0, min(1.0, zg / rg)) if rg > 0 else 0.0
                lat_g = math.degrees(math.asin(sin_lat))
                # Back to mean ecliptic of date
                lon_g, lat_g = _precess_ecliptic(lon_g, lat_g, 2451545.0, jd)
                return lon_g, lat_g, rg

            if (iflag & FLG_BARYCTR) and not (iflag & FLG_HELCTR):
                # Barycentric place: the model is heliocentric, so add the
                # Sun's barycentric vector (the same construction the Uranian
                # branch uses for its own analytic bodies). Without it the
                # request fell through to the geocentric Earth-subtraction
                # path and returned geocentric coordinates while advertising
                # BARYCTR.
                _lt_b = _helio_light_time_days(jd_tt)
                bx, by, bz = _helio_xyz_j2000(jd_tt - _lt_b)
                _sun_b = planets["sun"].at(get_timescale().tt_jd(jd_tt))
                _sxyz = _sun_b.frame_xyz(ecliptic_J2000_frame).au
                bx += float(_sxyz[0])
                by += float(_sxyz[1])
                bz += float(_sxyz[2])
                _rb = math.sqrt(bx * bx + by * by + bz * bz)
                _lon_b = math.degrees(math.atan2(by, bx)) % 360.0
                _lat_b = math.degrees(
                    math.asin(max(-1.0, min(1.0, bz / _rb)) if _rb > 0 else 0.0)
                )
                # The vectors above are J2000-ecliptic; rotate to the mean
                # ecliptic of date for the shared treatment below.
                _lon_b, _lat_b = _precess_ecliptic(_lon_b, _lat_b, 2451545.0, jd_tt)
                pos = (_lon_b, _lat_b, _rb, 0.0, 0.0, 0.0)
                if iflag & FLG_SPEED:
                    _h = 0.01
                    _prev_lt = _helio_light_time_days(jd_tt - _h)
                    _next_lt = _helio_light_time_days(jd_tt + _h)

                    def _bary_lon_lat_dist(_jd, _lt):
                        _x, _y, _z = _helio_xyz_j2000(_jd - _lt)
                        _sv = (
                            planets["sun"]
                            .at(get_timescale().tt_jd(_jd))
                            .frame_xyz(ecliptic_J2000_frame)
                            .au
                        )
                        _x += float(_sv[0])
                        _y += float(_sv[1])
                        _z += float(_sv[2])
                        _r = math.sqrt(_x * _x + _y * _y + _z * _z)
                        _lo = math.degrees(math.atan2(_y, _x)) % 360.0
                        _la = math.degrees(
                            math.asin(max(-1.0, min(1.0, _z / _r)) if _r > 0 else 0.0)
                        )
                        return _precess_ecliptic(_lo, _la, 2451545.0, _jd) + (_r,)

                    _p0 = _bary_lon_lat_dist(jd_tt - _h, _prev_lt)
                    _p1 = _bary_lon_lat_dist(jd_tt + _h, _next_lt)
                    _dlon = (_p1[0] - _p0[0] + 180.0) % 360.0 - 180.0
                    pos = (
                        _lon_b,
                        _lat_b,
                        _rb,
                        _dlon / (2.0 * _h),
                        (_p1[1] - _p0[1]) / (2.0 * _h),
                        (_p1[2] - _p0[2]) / (2.0 * _h),
                    )
            elif iflag & FLG_HELCTR:
                # Apparent heliocentric place: the body as seen from the Sun,
                # i.e. evaluated at jd - (Sun-to-body light time). TRUEPOS
                # zeroes that light time and yields the geometric place.
                pos = _calc_helio(jd_tt - _helio_light_time_days(jd_tt))
                if iflag & (FLG_SPEED | FLG_SPEED3):
                    # The model's own rate is dλ/dt at the RETARDED epoch; the
                    # reported quantity is λ(jd − lt(jd)), whose derivative
                    # carries the extra (1 − dlt/dt) factor. Keeping the raw
                    # model rate left the speed off the derivative of the
                    # reported position by 0.25"/day on the fastest body
                    # (Vulcan). Central-difference the reported quantity.
                    _hh = 0.01
                    _q0 = _calc_helio(
                        (jd_tt - _hh) - _helio_light_time_days(jd_tt - _hh)
                    )
                    _q1 = _calc_helio(
                        (jd_tt + _hh) - _helio_light_time_days(jd_tt + _hh)
                    )
                    _dl = ((_q1[0] - _q0[0] + 180.0) % 360.0 - 180.0) / (2.0 * _hh)
                    pos = (
                        pos[0],
                        pos[1],
                        pos[2],
                        _dl,
                        (_q1[1] - _q0[1]) / (2.0 * _hh),
                        (_q1[2] - _q0[2]) / (2.0 * _hh),
                    )
                if ipl in _FICT_HELIO_IDS:
                    # calc_fictitious_position returns J2000 ecliptic for the
                    # predicted planets; precess to the mean ecliptic of date so
                    # the shared "mean ecliptic of date (+ nutation)" treatment
                    # below is frame-consistent. Without this, a J2000 longitude
                    # would get Δψ added (and be reported) without ever being
                    # precessed to of-date. VULCAN/PROSERPINA already use the
                    # equinox of date and need no rotation.
                    if iflag & FLG_SPEED:
                        h_lon, h_lat, h_dlon, h_dlat = _precess_ecliptic_state(
                            pos[0],
                            pos[1],
                            pos[3],
                            pos[4],
                            jd_tt,
                            to_j2000=False,
                        )
                        pos = (h_lon, h_lat, pos[2], h_dlon, h_dlat, pos[5])
                    else:
                        h_lon, h_lat = _precess_ecliptic(
                            pos[0], pos[1], 2451545.0, jd_tt
                        )
                        pos = (h_lon, h_lat, pos[2], 0.0, 0.0, 0.0)
            else:
                g_lon, g_lat, g_dist = _geo_of_date(jd_tt)
                # Vulcan's ~18.6-day orbit makes it move ~3 deg/day
                # geocentrically, so a 0.1-day central difference over-smooths
                # its speed by ~1"/day; use a tighter step for it. The slow
                # predicted planets / Proserpina stay at 0.1 day.
                dt_v = 0.01 if ipl == VULCAN else 0.1
                prev = _geo_of_date(jd_tt - dt_v)
                nxt = _geo_of_date(jd_tt + dt_v)
                g_dlon = (nxt[0] - prev[0]) / (2.0 * dt_v)
                if g_dlon > 180.0 / (2.0 * dt_v):
                    g_dlon -= 360.0 / (2.0 * dt_v)
                elif g_dlon < -180.0 / (2.0 * dt_v):
                    g_dlon += 360.0 / (2.0 * dt_v)
                g_dlat = (nxt[1] - prev[1]) / (2.0 * dt_v)
                g_ddist = (nxt[2] - prev[2]) / (2.0 * dt_v)
                pos = (g_lon, g_lat, g_dist, g_dlon, g_dlat, g_ddist)
        else:
            # Selena (56) and Waldemath (58) are GEOCENTRIC symbolic models.
            # Their centre and light-path flags were inert while the retflag
            # echoed them, so a heliocentric request came back geocentric
            # (measured 140 deg off on 56, 72 deg on 58). Re-centre the
            # geocentric vector on the requested origin and honour the
            # aberration toggle, exactly as the Uranian lunar points do.
            pos = hypothetical.calc_hypothetical_position(ipl, jd_tt)
            if iflag & (FLG_HELCTR | FLG_BARYCTR) or not (
                iflag & (FLG_TRUEPOS | FLG_NOABERR)
            ):
                from skyfield.framelib import ecliptic_J2000_frame as _eclJ2

                _ts_h = get_timescale()
                _t_h = _ts_h.tt_jd(jd_tt)
                _lo = math.radians(pos[0])
                _la = math.radians(pos[1])
                _cb = math.cos(_la)
                _gx = pos[2] * _cb * math.cos(_lo)
                _gy = pos[2] * _cb * math.sin(_lo)
                _gz = pos[2] * math.sin(_la)
                if iflag & (FLG_HELCTR | FLG_BARYCTR):
                    # FLG_HELCTR adds Earth's heliocentric vector, FLG_BARYCTR
                    # its barycentric one; with both set the heliocentric view
                    # wins, mirroring the other analytic branches.
                    if iflag & FLG_HELCTR:
                        _ev = (
                            planets["sun"]
                            .at(_t_h)
                            .observe(planets["earth"])
                            .frame_xyz(_eclJ2)
                            .au
                        )
                    else:
                        _ev = planets["earth"].at(_t_h).frame_xyz(_eclJ2).au
                    _gx += float(_ev[0])
                    _gy += float(_ev[1])
                    _gz += float(_ev[2])
                elif not (iflag & (FLG_TRUEPOS | FLG_NOABERR)):
                    # Annual aberration on the geocentric direction (the
                    # observer's barycentric velocity over c).
                    _rg0 = math.sqrt(_gx * _gx + _gy * _gy + _gz * _gz)
                    _vel = (
                        planets["earth"]
                        .at(_t_h)
                        .frame_xyz_and_velocity(_eclJ2)[1]
                        .au_per_d
                    )
                    _c = 173.1446326846693
                    _gx += float(_vel[0]) * _rg0 / _c
                    _gy += float(_vel[1]) * _rg0 / _c
                    _gz += float(_vel[2]) * _rg0 / _c
                _rg = math.sqrt(_gx * _gx + _gy * _gy + _gz * _gz)
                _new_lon = math.degrees(math.atan2(_gy, _gx)) % 360.0
                _new_lat = math.degrees(
                    math.asin(max(-1.0, min(1.0, _gz / _rg)) if _rg > 0 else 0.0)
                )
                pos = (_new_lon, _new_lat, _rg, pos[3], pos[4], pos[5])

        lon, lat, dist = pos[0], pos[1], pos[2]
        # Zero the speed slots when FLG_SPEED is absent (the reference
        # returns zeros; the Uranian/Transpluto branches and the standard
        # planet path already gate this way).
        if iflag & FLG_SPEED:
            dlon, dlat, ddist = pos[3], pos[4], pos[5]
        else:
            dlon, dlat, ddist = 0.0, 0.0, 0.0

        # These functions return mean ecliptic of date. Add nutation unless
        # suppressed (NONUT, SIDEREAL+EQUATORIAL, or J2000 output, which is
        # the mean ecliptic precessed to J2000 like the reference).
        _sid_eq = is_sidereal and bool(iflag & FLG_EQUATORIAL)
        _no_nut = bool(iflag & (FLG_NONUT | FLG_J2000)) or _sid_eq
        if not _no_nut:
            from .cache import get_cached_nutation

            dpsi_rad, _ = get_cached_nutation(jd_tt)
            lon = (lon + math.degrees(dpsi_rad)) % 360.0
            # The true-ecliptic longitude carries Δψ, so its speed carries the
            # nutation-in-longitude rate dΔψ/dt (~0.05-0.17"/day). Add it, the
            # same way the MEAN_NODE / Uranian paths do, so the reported dlon is
            # the true derivative of the reported (nutated) longitude.
            if iflag & FLG_SPEED:
                dlon += _nutation_rate_deg_per_day(jd_tt)

        # Topocentric parallax is applied to the tropical of-date place first,
        # so any sidereal ayanamsha rotation below acts on the topocentric
        # longitude (keeping the observer offset in the correct frame).
        result = (lon, lat, dist, dlon, dlat, ddist)
        result = _apply_hypothetical_topocentric(result, jd_tt, iflag, False, planets)
        lon, lat, dist, dlon, dlat, ddist = result

        if is_sidereal and not (iflag & FLG_EQUATORIAL):
            lon, dlon = _apply_sidereal_correction(lon, dlon, t.ut1, iflag)

        result = (lon, lat, dist, dlon, dlat, ddist)
        result = _maybe_equatorial_convert(result, jd_tt, iflag)
        return _to_native_floats(result), iflag

    # Handle minor bodies (asteroids and TNOs)
    # Strategy: try to get a type21 VectorFunction target so we can route
    # through the Skyfield observe/apparent pipeline (same as planets).
    # This avoids the ~0.3" systematic error from the legacy ecliptic J2000
    # + manual precession/nutation approach in _calc_type21_position.
    # Arbitrary numbered asteroids (AST_OFFSET + N, below the fixed-star
    # range) take the same pipeline: registered SPK -> auto-SPK ->
    # Keplerian elements (curated table or JPL SBDB) -> UnknownBodyError.
    _spk_type21_target = None
    if ipl in minor_bodies.MINOR_BODY_ELEMENTS or (AST_OFFSET < ipl < FIXSTAR_OFFSET):
        from .state import get_calc_mode

        if get_calc_mode() == "leb":
            if ipl not in minor_bodies.MINOR_BODY_ELEMENTS:
                from .exceptions import UnknownBodyError

                raise UnknownBodyError(
                    message=(
                        f"Asteroid {ipl - AST_OFFSET} (body ID {ipl}) has no "
                        "curated local model and online/SPK sources are disabled "
                        "in calculation mode 'leb'."
                    ),
                    body_id=ipl,
                )
            return _calc_keplerian_fallback(t, ipl, iflag, planets)

        from . import spk
        from .state import get_auto_spk_download, get_strict_precision
        from .exceptions import SPKRequiredError
        from .constants import SPK_AUTO_DOWNLOAD_BLOCKED, SPK_BODY_NAME_MAP
        from .logging_config import get_logger

        # First check if already registered
        _spk_type21_target = spk.get_spk_type21_target(ipl, t.tt)
        auto_download_attempted = False

        if _spk_type21_target is None:
            # Not registered yet — try auto-download, then check again
            if get_auto_spk_download():
                auto_download_attempted = True
                try:
                    # _try_auto_spk_download registers the SPK as a side effect
                    _try_auto_spk_download(t, ipl, iflag)
                except (OSError, ValueError, KeyError, RuntimeError, TypeError):
                    pass
                # Re-check after download
                _spk_type21_target = spk.get_spk_type21_target(ipl, t.tt)

        if _spk_type21_target is not None:
            # Route through the planet pipeline below (observe/apparent)
            pass
        else:
            # Fallback: try the legacy calc_spk_body_position (non-type21 SPK)
            try:
                spk_result = spk.calc_spk_body_position(t, ipl, iflag)
            except EphemerisRangeError:
                # Outside the registered kernel's coverage — continue to
                # the Keplerian path (documented out-of-coverage behavior).
                spk_result = None
            if spk_result is not None:
                # _calc_type21_position performs the full frame reduction
                # internally (precession/nutation per FLG_J2000/FLG_NONUT),
                # so a type-21 result under FLG_J2000 is ALREADY J2000
                # ecliptic: pass already_j2000=True or the converter would
                # precess it date->J2000 a second time. Type-2/3 results are
                # of-date and must be precessed exactly once here.
                _spk_info = spk.get_spk_body_info(ipl)
                _already_j2000 = (
                    _spk_info is not None and spk._detect_spk_type(_spk_info[0]) == 21
                )
                spk_result = _maybe_equatorial_convert(
                    spk_result, t.tt, iflag, already_j2000=_already_j2000
                )
                result = _to_native_floats(spk_result)
                _mark_dispatch_source("SPK")
                return result, iflag

            # In strict precision mode, require an SPK-grade source when the
            # caller disabled auto-download.  ``auto``/``skyfield`` historically
            # attempt provisioning and then remain available through the local
            # model when the network is unavailable; preserve that public
            # contract outside sealed LEB mode.  LEB never reaches this branch.
            #
            # When no download was attempted, key the decision off whether a
            # source better than Keplerian exists for this epoch. If a kernel is
            # registered but out of coverage and ASSIST is unavailable, report
            # that exact condition instead of asking for a kernel already there.
            if (
                get_strict_precision()
                and not auto_download_attempted
                and ipl in SPK_BODY_NAME_MAP
            ):
                if ipl not in SPK_AUTO_DOWNLOAD_BLOCKED:
                    spk_info = spk.get_spk_body_info(ipl)
                    if _strict_source_better_than_keplerian(ipl):
                        # A high-precision N-body fallback (ASSIST) is available
                        # locally -- fall through to it below, no network needed.
                        pass
                    elif spk_info is not None:
                        # A kernel is registered but does not cover this epoch,
                        # and no ASSIST data is present: refuse, and say so
                        # accurately (do not tell the user to register an SPK
                        # that is already registered).
                        raise _spk_out_of_coverage_error(ipl, spk_info, t.tt)
                    else:
                        # No kernel registered at all and no ASSIST: the original
                        # strict-mode contract -- register an SPK (or enable
                        # auto-download).
                        horizons_id, _ = SPK_BODY_NAME_MAP[ipl]
                        body_name = spk._get_body_name(ipl) or str(ipl)
                        raise SPKRequiredError.for_body(ipl, body_name, horizons_id)

            jd_tt = t.tt

            # Try ASSIST N-body integration fallback if available
            # (only for bodies with curated elements; arbitrary SBDB
            # asteroids have no initial conditions for the integrator)
            try:
                from .rebound_integration import check_assist_data_available

                if ipl in minor_bodies.MINOR_BODY_ELEMENTS and (
                    check_assist_data_available()
                ):
                    assist_state = _assist_state_at(jd_tt, ipl)
                    lon, lat, dist = _assist_position_from_state(
                        assist_state, jd_tt, iflag, planets
                    )

                    speed_lon = 0.0
                    speed_lat = 0.0
                    speed_dist = 0.0
                    if iflag & FLG_SPEED:
                        dt = 1.0 / 86400.0
                        prev_state = _assist_shift_state(assist_state, jd_tt - dt)
                        next_state = _assist_shift_state(assist_state, jd_tt + dt)
                        lon_prev, lat_prev, dist_prev = _assist_position_from_state(
                            prev_state, jd_tt - dt, iflag, planets
                        )
                        lon_next, lat_next, dist_next = _assist_position_from_state(
                            next_state, jd_tt + dt, iflag, planets
                        )
                        speed_lon = (lon_next - lon_prev) / (2.0 * dt)
                        speed_lat = (lat_next - lat_prev) / (2.0 * dt)
                        speed_dist = (dist_next - dist_prev) / (2.0 * dt)

                        if speed_lon > 180.0 / (2.0 * dt):
                            speed_lon -= 360.0 / (2.0 * dt)
                        if speed_lon < -180.0 / (2.0 * dt):
                            speed_lon += 360.0 / (2.0 * dt)

                    # Sidereal offset (ecliptic only) — like every other body
                    # branch; without this the ASSIST fallback silently returned
                    # a tropical longitude for FLG_SIDEREAL requests.
                    if (iflag & FLG_SIDEREAL) and not (iflag & FLG_EQUATORIAL):
                        lon, speed_lon = _apply_sidereal_correction(
                            lon, speed_lon, t.ut1, iflag
                        )

                    result = _to_native_floats(
                        _maybe_equatorial_convert(
                            (lon, lat, dist, speed_lon, speed_lat, speed_dist),
                            jd_tt,
                            iflag,
                        )
                    )
                    _mark_dispatch_source("ASSIST")
                    return result, iflag
            except (ImportError, RuntimeError, ValueError, FileNotFoundError):
                pass

            # Keplerian as the final, explicitly traced approximation.
            return _calc_keplerian_fallback(t, ipl, iflag, planets)

    # Handle fixed stars — delegate to the fixstar_ut() computation so every
    # flag (FLG_SIDEREAL, FLG_EQUATORIAL, FLG_SPEED, ...) gets the exact same
    # handling as the dedicated entry point. Dispatch BY ID: resolving the
    # traditional name is ambiguous (Flamsteed names starting with digits,
    # e.g. "29Psc", parse as sequential catalog numbers).
    if ipl in fixed_stars.FIXED_STARS:
        star_name = fixed_stars.get_canonical_star_name(ipl)
        pos, _star, _retflag = fixed_stars._fixstar_ut_by_id(
            ipl, star_name, t.ut1, iflag, record_trace=False
        )
        result = _to_native_floats(pos)
        _mark_dispatch_source(fixed_stars._get_fixed_star_source())
        return result, iflag

    # Handle astrological angles (requires observer location)
    if ANGLE_OFFSET <= ipl < ARABIC_OFFSET:
        topo = get_topo()
        if topo is None:
            from .exceptions import ConfigurationError

            raise ConfigurationError(
                "Angles require observer location. Call set_topo() first.",
                missing_config="observer_location",
                suggestion="Call set_topo(lon, lat, alt) first",
            )

        # Extract lat/lon from topo
        lat = topo.latitude.degrees
        lon = topo.longitude.degrees
        jd_ut = t.ut1

        angle_val = angles.get_angle_value(ipl, jd_ut, lat, lon)
        return (angle_val, 0.0, 0.0, 0.0, 0.0, 0.0), iflag

    # Handle Arabic parts (requires cached planet positions)
    if ARABIC_OFFSET <= ipl < ARABIC_OFFSET + 100:
        cache = get_angles_cache()
        if not cache:
            from .exceptions import ConfigurationError

            raise ConfigurationError(
                "Arabic parts require pre-calculated positions. "
                "Call calc_angles() first.",
                missing_config="angles_cache",
                suggestion="Call calc_angles() first",
            )

        # Map part IDs to calculation functions
        if ipl == PARS_FORTUNAE:
            asc = cache.get("Asc", cache.get("Ascendant", 0))
            sun = cache.get("Sun", 0)
            moon = cache.get("Moon", 0)
            is_diurnal = arabic_parts.is_day_chart(sun, asc)
            lon = arabic_parts.calc_arabic_part_of_fortune(asc, sun, moon, is_diurnal)
        elif ipl == PARS_SPIRITUS:
            asc = cache.get("Asc", cache.get("Ascendant", 0))
            sun = cache.get("Sun", 0)
            moon = cache.get("Moon", 0)
            is_diurnal = arabic_parts.is_day_chart(sun, asc)
            lon = arabic_parts.calc_arabic_part_of_spirit(asc, sun, moon, is_diurnal)
        elif ipl == PARS_AMORIS:
            asc = cache.get("Asc", cache.get("Ascendant", 0))
            venus = cache.get("Venus", 0)
            sun = cache.get("Sun", 0)
            lon = arabic_parts.calc_arabic_part_of_love(asc, venus, sun)
        elif ipl == PARS_FIDEI:
            asc = cache.get("Asc", cache.get("Ascendant", 0))
            mercury = cache.get("Mercury", 0)
            moon = cache.get("Moon", 0)
            lon = arabic_parts.calc_arabic_part_of_faith(asc, mercury, moon)
        else:
            lon = 0.0

        return _to_native_floats((lon, 0.0, 0.0, 0.0, 0.0, 0.0)), iflag

    # Handle standard planets (and type21 asteroids routed through planet pipeline)
    if _spk_type21_target is not None:
        # Type21 asteroid: use the VectorFunction wrapper for Skyfield pipeline
        target = _spk_type21_target
    elif ipl in _PLANET_MAP:
        target_name = _PLANET_MAP[ipl]
        target = get_planet_target(planets, target_name)
    else:
        # Unknown body - raise clear error instead of returning zeros
        from .exceptions import UnknownBodyError

        raise UnknownBodyError(
            message=(
                f"Unknown body ID {ipl}. "
                f"Supported bodies include: standard planets (0-14), lunar nodes (10-11), "
                f"Lilith/apogee (12-13, 21-22), asteroids (15-20), "
                f"historical fictitious bodies (40-58), minor bodies (AST_OFFSET+number), "
                f"and fixed stars (FIXSTAR_OFFSET+number). "
                f"See libephemeris.constants for all body ID constants."
            ),
            body_id=ipl,
        )

    # 2. Identify Observer
    observer_topo = get_topo()

    is_barycentric = bool(iflag & FLG_BARYCTR)

    # Earth geocentric is trivially (0,0,0,0,0,0) regardless of frame flags.
    # Return early to avoid division-by-zero in J2000/ICRS coordinate transforms
    # where dist=0 would cause NaN from asin(ze/dist). A sidereal ecliptic
    # request keeps the zero vector: measured reference behavior does not
    # subtract the ayanamsha from the undefined zero-length longitude
    # (see _degenerate_origin_result).
    if ipl == EARTH and not (iflag & FLG_HELCTR) and not is_barycentric:
        _mark_dispatch_source("Analytical")
        return _to_native_floats(_degenerate_origin_result()), iflag

    if iflag & FLG_HELCTR:
        # Heliocentric
        observer = planets["sun"]
        icrf_center = 10  # Sun
    elif is_barycentric:
        # Barycentric - position relative to Solar System Barycenter (SSB)
        # In Skyfield, the SSB is the origin (center=0), so we don't need an observer
        observer = None
        icrf_center = 0  # SSB
    elif iflag & FLG_TOPOCTR:
        if observer_topo is None:
            from .exceptions import ConfigurationError

            # Matching the reference API: topocentric positions without a
            # prior set_topo() are an error, not a silent geocentric result.
            raise ConfigurationError(
                "FLG_TOPOCTR requires a geographic position: call "
                "set_topo(lon, lat, alt) first",
                missing_config="observer_location",
                suggestion="Call set_topo(lon, lat, alt) first",
            )
        earth = planets["earth"]
        observer = earth + observer_topo
        icrf_center = observer_topo
    else:
        # Geocentric
        observer = planets["earth"]
        icrf_center = 399  # Earth

    # 3. Compute Position
    from .cache import get_cached_observer_at

    # Public planet identifiers denote planet centres. ``get_planet_target``
    # supplies a JPL centre segment when available and otherwise reconstructs
    # the centre from the system barycentre and independently sourced satellite
    # theory.  That same target is used for every observer centre and flag set.
    out_target = target

    # Helper to get vector at time t
    def get_vector(t_):
        # Target position relative to SSB
        tgt = out_target.at(t_)
        tgt_pos = tgt.position.au
        tgt_vel = tgt.velocity.au_per_d

        if observer is None:
            # Barycentric: target position is already relative to SSB
            return tgt_pos, tgt_vel

        # Observer relative to SSB (cached per observer+JD)
        obs = get_cached_observer_at(observer, t_)
        obs_pos = obs.position.au
        obs_vel = obs.velocity.au_per_d

        p_ = tgt_pos - obs_pos
        v_ = tgt_vel - obs_vel
        return p_, v_

    if iflag & FLG_TRUEPOS:
        # Geometric position (instantaneous)
        p, v = get_vector(t)
        from skyfield.positionlib import ICRF

        pos = ICRF(p, v, t=t, center=icrf_center)
    else:
        # Apparent position
        if (iflag & FLG_HELCTR) or (iflag & FLG_BARYCTR):
            # For SSB or Heliocentric, we need to apply light-time correction
            # Light-time correction: position shows where object WAS
            # when light left it to reach the observer (Sun for heliocentric)
            import numpy as np

            # Speed of light in AU/day: NOVAS/Skyfield C_AUDAY realization
            # (c exact from SI, IAU 1976 au), shared with astrometry/fast_calc.
            C_AU_PER_DAY = 173.1446326846693

            # Get initial geometric position
            p, v = get_vector(t)

            # Iterate the physical observer-to-target light path. For HELCTR
            # the observer is the Sun and the baseline is target-Sun; for
            # BARYCTR the observer is the SSB and the baseline is target-SSB.
            lt_target = out_target

            # Light-time iteration: retard only the target; the observer
            # (Sun or SSB) stays at the observation time. Retarding both
            # shifts the result by the observer's barycentric motion over
            # the light time (up to ~10 mas for the outer planets) and
            # diverges from the LEB pipeline.
            if observer is not None:
                _obs_at = observer.at(t)
                _obs_pos_lt = _obs_at.position.au
                _obs_vel_lt = _obs_at.velocity.au_per_d
            else:
                _obs_pos_lt = np.zeros(3)
                _obs_vel_lt = np.zeros(3)

            for _ in range(3):
                dist = np.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2)
                light_time = dist / C_AU_PER_DAY
                ts_lt = get_timescale()
                # Retarded epoch as a two-part JD: subtract the light time from
                # the fractional part and keep t.whole intact. Collapsing
                # t.tdb - light_time into a single float64 quantizes the retarded
                # instant to the JD ULP (~5.5e-10 d at JD ~2.46e6). The velocity
                # central difference re-evaluates this block at t ± dt (dt = 1 s),
                # so that per-sample quantization is divided by 2·dt and shows up
                # as a spurious HELCTR/BARYCTR speed bias that grows with the
                # body's motion (up to ~0.28"/day for Mercury near perihelion).
                # The two-part form preserves the exact ±dt spacing so the
                # reported speed is the true derivative of the reported position,
                # matching the LEB backend. (Same ULP hazard the sample-time
                # construction above guards against for t ± dt.)
                _tgt_ret = lt_target.at(
                    ts_lt.tdb_jd(t.whole, t.tdb_fraction - light_time)
                )
                p = _tgt_ret.position.au - _obs_pos_lt
                v = _tgt_ret.velocity.au_per_d - _obs_vel_lt

            from skyfield.positionlib import ICRF

            pos = ICRF(p, v, t=t, center=icrf_center)
        else:
            # Cache observer.at(t) to avoid recomputing Earth's position
            # for every planet at the same JD
            obs_at_t = get_cached_observer_at(observer, t)
            if (iflag & FLG_NOABERR) and (iflag & FLG_NOGDEFL):
                # Astrometric (the reference FLG_ASTROMETRIC combination)
                pos = obs_at_t.observe(target)
            elif iflag & FLG_NOABERR:
                # Reference semantics: NOABERR disables only aberration;
                # gravitational light deflection still applies (that is
                # why FLG_ASTROMETRIC is defined as NOABERR|NOGDEFL).
                astrometric = obs_at_t.observe(target)
                pos = _apply_deflection_only(astrometric, t, icrf_center, planets)
            elif iflag & FLG_NOGDEFL:
                # Aberration without gravitational deflection:
                # Pass empty deflectors tuple to skip Sun/Jupiter/Saturn
                # deflection while still applying stellar aberration.
                pos = obs_at_t.observe(target).apparent(deflectors=())
            else:
                # Apparent = deflection + aberration. The deflection is taken
                # from the SHARED clamped routine rather than Skyfield's, so
                # both backends use one implementation: Skyfield's
                # add_deflection has no ERFA dlim limiter, so for a planet
                # geometrically inside the solar disc its term blew up and the
                # two backends split by 14.3" in latitude (measured on Mars at
                # JD 2422088.26225, 0.017 deg from the Sun's centre, where the
                # LEB side sits 0.62" from the reference and the Skyfield side
                # 14.9"). Same deflectors (Sun, Jupiter, Saturn).
                import numpy as _np
                from skyfield.constants import C_AUDAY as _C_AUDAY
                from skyfield.functions import length_of as _length_of
                from skyfield.positionlib import Apparent as _Apparent
                from skyfield.relativity import add_aberration as _add_aberration

                _astro = obs_at_t.observe(target)
                _defl = _apply_deflection_only(_astro, t, icrf_center, planets)
                _p = _np.array(_defl.position.au, dtype=float)
                _v_obs = _np.array(
                    _astro.center_barycentric.velocity.au_per_d, dtype=float
                )
                _add_aberration(_p, _v_obs, _length_of(_p) / _C_AUDAY)
                pos = _Apparent(_p, _defl.velocity.au_per_d, t=t, center=icrf_center)

    # 4. Coordinate System & Speeds
    is_equatorial = bool(iflag & FLG_EQUATORIAL)
    is_icrs = bool(iflag & FLG_ICRS)
    is_sidereal = bool(iflag & FLG_SIDEREAL)

    p1, p2, p3 = 0.0, 0.0, 0.0
    dp1, dp2, dp3 = 0.0, 0.0, 0.0

    # Get position and velocity vectors in AU and AU/day
    # We need them in the correct frame.
    # Skyfield's pos.position.au and pos.velocity.au_per_d are in ICRS (Equatorial J2000).
    # If we want Ecliptic, we need to rotate them.

    # Define rotation matrix or use Skyfield's frame transform
    # Skyfield doesn't easily rotate velocity vectors with frame_latlon.
    # We have to do it manually or use `frame_xyz(frame)`.

    if is_equatorial:
        # Equatorial coordinates (Right Ascension / Declination)
        # Frame options: ICRS (J2000) or True Equator of Date

        if iflag & FLG_J2000:
            # Mean equator/equinox of J2000: ICRS rotated by IAU 2006 frame
            # bias. FLG_ICRS combines with FLG_J2000 instead of being
            # overridden by it (measured: the bias is still dropped, 0.0065"
            # to 0.0102" on the bodies sampled), so under ICRS the vector is
            # already in the requested frame.
            from .fast_calc import _rotate_icrs_to_j2000_mean_equator

            _xi, _yi, _zi = pos.position.au
            if is_icrs:
                xe, ye, ze = float(_xi), float(_yi), float(_zi)
            else:
                xe, ye, ze = _rotate_icrs_to_j2000_mean_equator(
                    float(_xi), float(_yi), float(_zi)
                )
            dist = math.sqrt(xe * xe + ye * ye + ze * ze)
            p1 = math.degrees(math.atan2(ye, xe)) % 360.0
            p2 = (
                math.degrees(math.asin(max(-1.0, min(1.0, ze / dist))))
                if dist > 0
                else 0.0
            )
            p3 = dist
            # Speeds are computed once in the generic central-difference
            # block below (section 4); an in-branch duplicate previously
            # ran two extra full pipeline evaluations whose results were
            # unconditionally discarded.
        else:
            # Equator of Date (true or mean depending on NONUT/SIDEREAL flags)
            # Use the already-computed pos (with light-time correction) for the
            # main position. This is critical for HELCTR/BARYCTR where pos
            # includes iterative light-time correction applied in section 3.
            # SIDEREAL+EQUATORIAL uses the public mean-equator convention,
            # equivalent to NONUT frame selection.
            _use_mean_equator = bool(iflag & FLG_NONUT) or is_sidereal

            if is_icrs:
                # ICRS equator of date (no frame bias): Vondrák 2011 precession,
                # with nutation unless the mean equator is requested.
                import numpy as np

                xyz_icrs = np.array(pos.position.au)
                if _use_mean_equator:
                    pn = vondrak_precession_matrix(t.tt, frame_bias=False)
                else:
                    dpsi, deps = t._nutation_angles_radians
                    pn, _eps_true = vondrak_pn_matrix(
                        t.tt, float(dpsi), float(deps), frame_bias=False
                    )
                xyz_eq = np.array(pn) @ xyz_icrs
                xe, ye, ze = float(xyz_eq[0]), float(xyz_eq[1]), float(xyz_eq[2])
                dist = math.sqrt(xe * xe + ye * ye + ze * ze)
                p1 = math.degrees(math.atan2(ye, xe)) % 360.0
                p2 = math.degrees(math.asin(max(-1.0, min(1.0, ze / dist))))
                p3 = dist
            elif _use_mean_equator:
                # Mean equator of date (with frame bias): Vondrák precession.
                import numpy as np

                pn = vondrak_precession_matrix(t.tt)
                xyz_eq = np.array(pn) @ np.array(pos.position.au)
                xe, ye, ze = float(xyz_eq[0]), float(xyz_eq[1]), float(xyz_eq[2])
                dist = math.sqrt(xe * xe + ye * ye + ze * ze)
                p1 = math.degrees(math.atan2(ye, xe)) % 360.0
                p2 = math.degrees(math.asin(max(-1.0, min(1.0, ze / dist))))
                p3 = dist
            else:
                # True equator of date (with frame bias): Vondrák precession +
                # Skyfield nutation, matching the LEB fast path.
                import numpy as np

                dpsi, deps = t._nutation_angles_radians
                pn, _eps_true = vondrak_pn_matrix(t.tt, float(dpsi), float(deps))
                xyz_eq = np.array(pn) @ np.array(pos.position.au)
                xe, ye, ze = float(xyz_eq[0]), float(xyz_eq[1]), float(xyz_eq[2])
                dist = math.sqrt(xe * xe + ye * ye + ze * ze)
                p1 = math.degrees(math.atan2(ye, xe)) % 360.0
                p2 = math.degrees(math.asin(max(-1.0, min(1.0, ze / dist))))
                p3 = dist

            # Speeds are computed once in the generic central-difference
            # block below (section 4); the previous in-branch duplicate
            # ran two extra full pipeline evaluations whose results were
            # unconditionally discarded before that block.

    else:
        # Ecliptic (Long/Lat)
        if iflag & FLG_J2000:
            # IAU 2006 mean-ecliptic J2000 coordinates (frame bias plus mean
            # obliquity). FLG_ICRS drops the bias here too, leaving the
            # pole/equinox on the ICRS realization and applying the J2000
            # mean obliquity alone.
            from .fast_calc import (
                _rotate_equatorial_to_ecliptic,
                _rotate_icrs_to_ecliptic_j2000,
            )

            x, y, z = pos.position.au
            if is_icrs:
                xe, ye, ze = _rotate_equatorial_to_ecliptic(
                    float(x),
                    float(y),
                    float(z),
                    vondrak_mean_obliquity_rad(_J2000_JD),
                )
            else:
                xe, ye, ze = _rotate_icrs_to_ecliptic_j2000(
                    float(x), float(y), float(z)
                )

            # Convert to spherical
            dist = math.sqrt(xe * xe + ye * ye + ze * ze)
            lon = math.degrees(math.atan2(ye, xe)) % 360.0
            lat = math.degrees(math.asin(max(-1.0, min(1.0, ze / dist))))

            p1, p2, p3 = lon, lat, dist

        else:
            # Ecliptic of date: reduce ICRS -> equator of date with Vondrák 2011
            # precession (the same single source as the LEB fast path), then
            # rotate equator -> ecliptic about the equinox by the obliquity. The
            # mean variant (FLG_NONUT) drops nutation; the frame bias is dropped
            # for FLG_ICRS.
            import numpy as np

            if iflag & FLG_NONUT:
                # Mean ecliptic of date (precession only, no nutation).
                eps = vondrak_mean_obliquity_rad(t.tt)
                pn = vondrak_precession_matrix(t.tt, frame_bias=not is_icrs)
            elif is_icrs:
                # ICRS ecliptic of date (no frame bias).
                dpsi, deps = t._nutation_angles_radians
                pn, eps = vondrak_pn_matrix(
                    t.tt, float(dpsi), float(deps), frame_bias=False
                )
            else:
                # True ecliptic of date.
                dpsi, deps = t._nutation_angles_radians
                pn, eps = vondrak_pn_matrix(t.tt, float(dpsi), float(deps))

            xq, yq, zq = np.array(pn) @ np.array(pos.position.au)
            ce, se = math.cos(eps), math.sin(eps)
            xe = float(xq)
            ye = float(yq) * ce + float(zq) * se
            ze = -float(yq) * se + float(zq) * ce
            dist = math.sqrt(xe * xe + ye * ye + ze * ze)
            p1 = math.degrees(math.atan2(ye, xe)) % 360.0
            p2 = math.degrees(math.asin(max(-1.0, min(1.0, ze / dist))))
            p3 = dist

    # 4. Speed (Central Difference Numerical Differentiation if requested)
    # Using central differences: f'(x) ≈ [f(x+h) - f(x-h)] / (2h)
    # This has error O(h²) compared to O(h) for forward differences,
    # providing ~100x better precision for the same timestep.
    #
    # The Moon uses a six-second half-step to resolve its faster apparent
    # motion; other bodies use one second. The geometric (FLG_TRUEPOS) place of
    # an SPK type-21 minor body uses a five-minute half-step so the type-21
    # sample jitter does not leak into the reported speed (see the constant's
    # note); the apparent place is smooth and keeps the one-second step.
    if _spk_type21_target is not None and (iflag & FLG_TRUEPOS) and ipl != MOON:
        dt = _ASTEROID_TRUEPOS_SPEED_HALF_STEP_DAYS
    else:
        dt = _MOON_SPEED_HALF_STEP_DAYS if ipl == MOON else _BODY_SPEED_HALF_STEP_DAYS
    dp1, dp2, dp3 = 0.0, 0.0, 0.0

    if iflag & FLG_SPEED:
        # Get positions at t - dt and t + dt for central difference.
        # Two-argument tt_jd(whole, fraction): collapsing t.tt ± dt into a
        # single float64 quantizes the step to the JD ULP (~4.6e-10 d at
        # JD ~2.46e6), which biases the central difference by up to ~0.08"/day
        # for short steps. Keeping dt in the fraction slot preserves the exact
        # ±dt spacing through
        # Skyfield's extended-precision time arithmetic.
        ts_inner = get_timescale()
        t_prev = ts_inner.tt_jd(t.tt, -dt)
        t_next = ts_inner.tt_jd(t.tt, dt)

        # Sample flags for the central difference. Strip FLG_SPEED always.
        #
        # For ECLIPTIC sidereal output, strip FLG_SIDEREAL too so both samples
        # are tropical of-date; the ayanamsha drift is then applied to dlon
        # analytically below.
        #
        # For EQUATORIAL sidereal output, KEEP FLG_SIDEREAL so the samples are
        # computed in the very frame of the reported position — the *mean*
        # equator of date. The reported speed must be the exact derivative of
        # the reported position (the project's certified true-rate divergence);
        # the SID|EQ position uses the mean equator, so its RA/Dec speed is the
        # mean-equator rate, NOT the true-equator (plain-EQ) rate. Stripping
        # FLG_SIDEREAL here would have produced a true-equator speed that does
        # not differentiate the reported mean-equator position (a ~0.8"/day RA,
        # ~1.8"/day Dec inconsistency for the Moon), and disagreed with the LEB
        # backend, which differentiates the reported position directly.
        sample_flags = iflag & ~FLG_SPEED
        if not is_equatorial:
            sample_flags &= ~FLG_SIDEREAL
        forced_barycenters = set(_spk_center_forced_barycenters.get())
        if isinstance(out_target, _SpkCenterTarget) and (
            out_target._speed_stencil_requires_fallback(dt)
        ):
            forced_barycenters.add(out_target._barycenter_name)
        center_token = _spk_center_forced_barycenters.set(frozenset(forced_barycenters))
        try:
            try:
                result_prev, _ = _calc_body(t_prev, ipl, sample_flags)
                result_next, _ = _calc_body(t_next, ipl, sample_flags)
            except TypeError:
                # Skyfield Time reify descriptor corruption: clear cache and retry
                # with fresh Time objects (see Skyfield #xxx)
                from .cache import clear_observer_cache

                clear_observer_cache()
                t_prev = ts_inner.tt_jd(float(t.tt), -dt)
                t_next = ts_inner.tt_jd(float(t.tt), dt)
                result_prev, _ = _calc_body(t_prev, ipl, sample_flags)
                result_next, _ = _calc_body(t_next, ipl, sample_flags)
        finally:
            _spk_center_forced_barycenters.reset(center_token)
        p1_prev, p2_prev, p3_prev = result_prev[0], result_prev[1], result_prev[2]
        p1_next, p2_next, p3_next = result_next[0], result_next[1], result_next[2]

        # Calculate derivatives using central difference formula
        # dp/dt = (p(t+dt) - p(t-dt)) / (2*dt)
        dp1 = (p1_next - p1_prev) / (2.0 * dt)
        dp2 = (p2_next - p2_prev) / (2.0 * dt)
        dp3 = (p3_next - p3_prev) / (2.0 * dt)

        # Handle longitude wrap-around for dp1
        # The threshold depends on dt: 180° / (2*dt) in degrees/day
        wrap_threshold = 180.0 / (2.0 * dt)
        if dp1 > wrap_threshold:
            dp1 -= 360.0 / (2.0 * dt)
        elif dp1 < -wrap_threshold:
            dp1 += 360.0 / (2.0 * dt)

    # 5. Sidereal Mode
    # Sidereal correction is applied to ecliptic longitude only.
    # The reference ephemeris ignores sidereal flag when outputting equatorial coords.
    # Use NONUT-aware ayanamsha when FLG_NONUT is set.
    if is_sidereal and not is_equatorial:
        ayanamsa = _get_ayanamsa_for_flags(t.ut1, iflag)
        p1 = (p1 - ayanamsa) % 360.0

        # Correct velocity for ayanamsha rate if speed was calculated
        # Central difference: (f(t+h) - f(t-h)) / (2h) for O(h²) precision.
        # Shortest-arc delta: the ayanamsha is mod 360 and the two samples
        # can straddle the 0/360 wrap (star-anchored modes cross 0 on
        # supported dates), which would spike the speed by ~360/(2*dt).
        if iflag & FLG_SPEED:
            ayanamsa_prev = _get_ayanamsa_for_flags(t.ut1 - dt, iflag)
            ayanamsa_next = _get_ayanamsa_for_flags(t.ut1 + dt, iflag)
            da = (ayanamsa_next - ayanamsa_prev + 180.0) % 360.0 - 180.0
            dp1 -= da / (2.0 * dt)

    result = _to_native_floats((p1, p2, p3, dp1, dp2, dp3))
    if _spk_type21_target is not None:
        _mark_dispatch_source("SPK")
    return result, iflag


@_contextmanager
def _swapped_context_state(ctx):
    """Temporarily install an EphemerisContext's per-instance state globally.

    The core calculation code (_calc_body, horizons_calc_ut, ayanamsha
    helpers) reads observer/sidereal settings from module globals; context
    calls swap them in under _CONTEXT_SWAP_LOCK so the save-set-restore
    cycle is atomic across threads.
    """
    from . import state

    with state._CONTEXT_SWAP_LOCK:
        saved = (
            state._TOPO,
            state._SIDEREAL_MODE,
            state._SIDEREAL_T0,
            state._SIDEREAL_AYAN_T0,
            state._SIDEREAL_BITS,
            state._ANGLES_CACHE,
        )
        state._TOPO = ctx.topo
        state._SIDEREAL_MODE = ctx.sidereal_mode
        state._SIDEREAL_T0 = ctx.sidereal_t0
        state._SIDEREAL_AYAN_T0 = ctx.sidereal_ayan_t0
        # The SIDBIT projection flags qualify the mode above, so they must
        # travel with it: left on the module value, a module-level
        # set_sid_mode(mode | SIDBIT_ECL_DATE) changed this context's result
        # despite the context carrying its own sidereal configuration.
        state._SIDEREAL_BITS = ctx.sidereal_bits
        state._ANGLES_CACHE = ctx._angles_cache
        try:
            yield
        finally:
            (
                state._TOPO,
                state._SIDEREAL_MODE,
                state._SIDEREAL_T0,
                state._SIDEREAL_AYAN_T0,
                state._SIDEREAL_BITS,
                state._ANGLES_CACHE,
            ) = saved


def _calc_body_with_context(
    t, ipl: int, iflag: int, ctx
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """
    Calculate position using an explicit EphemerisContext (thread-safe).

    This is a context-aware wrapper around the core calculation logic.
    It temporarily sets global state from context, calls _calc_body, then
    restores global state. This allows context-based thread-safe usage while
    reusing the existing calculation code.

    Args:
        t: Skyfield Time object (UT1 or TT)
        ipl: Planet/body ID
        iflag: Calculation flags
        ctx: EphemerisContext instance containing state

    Returns:
        Same as _calc_body: ((lon, lat, dist, dlon, dlat, ddist), retflag)

    Thread Safety:
        This function acquires state._CONTEXT_SWAP_LOCK to ensure that the
        save-set-restore cycle is atomic across threads. Without this lock,
        concurrent calls could interleave and corrupt each other's state.
    """
    with _swapped_context_state(ctx):
        return _calc_body(t, ipl, iflag)


def _calc_body_pctr_with_context(
    t, ipl: int, iplctr: int, iflag: int, ctx
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """
    Calculate planet-centric position using an explicit EphemerisContext (thread-safe).

    This is a context-aware wrapper around _calc_body_pctr.

    Args:
        t: Skyfield Time object (UT1 or TT)
        ipl: Target planet/body ID
        iplctr: Observer/center planet ID
        iflag: Calculation flags
        ctx: EphemerisContext instance containing state

    Returns:
        Same as _calc_body_pctr: ((lon, lat, dist, dlon, dlat, ddist), retflag)

    Thread Safety:
        This function acquires state._CONTEXT_SWAP_LOCK to ensure that the
        save-set-restore cycle is atomic across threads.
    """
    with _swapped_context_state(ctx):
        return _calc_body_pctr(t, ipl, iplctr, iflag)


def get_ayanamsa_ut(tjdut: float) -> float:
    """
    Calculate ayanamsa (sidereal offset) for a given Universal Time date.

    Returns the ayanamsa in degrees for the currently set sidereal mode.
    The ayanamsa represents the longitudinal offset between tropical and
    sidereal zodiacs. Use set_sid_mode() to select the ayanamsa system.

    Args:
        tjdut: Julian Day in Universal Time (UT1)

    Returns:
        Ayanamsa value in degrees (tropical_longitude - sidereal_longitude)

    Example:
        >>> set_sid_mode(SIDM_J2000)
        >>> ayanamsa = get_ayanamsa_ut(2451545.0)  # J2000.0
        >>> print(f"J2000 ayanamsha: {ayanamsa:.6f}°")
    """
    sid_mode, sid_t0, sid_ayan_t0 = get_sid_mode(full=True)
    # get_sid_mode(full=True) returns (int, float, float)
    assert isinstance(sid_mode, int)

    # SIDM_USER + SIDBIT_USER_UT interprets the user t0 as a UT date; that
    # conversion lives in _calc_ayanamsa, so skip the LEB fast path (which takes
    # sid_t0 verbatim) and let the two agree.
    _user_ut = sid_mode == SIDM_USER and bool(_get_sidereal_bits() & SIDBIT_USER_UT)

    # --- LEB fast path: compute ayanamsa from precomputed data ---
    from .state import get_leb_reader

    reader = get_leb_reader()
    if reader is not None and not _user_ut:
        try:
            from .fast_calc import _calc_ayanamsa_from_leb

            delta_t = reader.delta_t(tjdut)
            jd_tt = tjdut + delta_t
            return float(
                _calc_ayanamsa_from_leb(
                    reader,
                    jd_tt,
                    sid_mode=sid_mode,
                    sid_t0=sid_t0,
                    sid_ayan_t0=sid_ayan_t0,
                )
            )
        except (KeyError, ValueError):
            pass  # star-based mode or out of range, fall through
    # --- END LEB fast path ---

    return float(_calc_ayanamsa(tjdut, sid_mode))


def get_ayanamsa_name(sidmode: int) -> str:
    """
    Get the name of a sidereal mode.

    Compatible with the reference get_ayanamsa_name(). Measured reference
    behavior: only the predefined modes 0--46 carry a name; every id without a
    predefined name -- the unassigned block above the last predefined mode and
    the user-defined mode SIDM_USER (255) -- returns the empty string rather
    than a placeholder.

    The SIDBIT_* projection flags (SIDBIT_ECL_T0=256, SIDBIT_SSY_PLANE=512,
    SIDBIT_USER_UT=1024, SIDBIT_ECL_DATE=2048, SIDBIT_NO_PREC_OFFSET=4096,
    SIDBIT_PREC_ORIG=8192) may be OR-ed onto the base mode id. They occupy
    bits >= 8, so the name lookup is driven by the low byte only: e.g.
    get_ayanamsa_name(SIDBIT_ECL_T0 | SIDM_LAHIRI) == "Lahiri". Masking with
    0xFF matches the measured reference behavior across every combination of
    projection flags.
    """
    names = {
        SIDM_FAGAN_BRADLEY: "Fagan/Bradley",
        SIDM_LAHIRI: "Lahiri",
        SIDM_DELUCE: "De Luce",
        SIDM_RAMAN: "Raman",
        SIDM_USHASHASHI: "Usha/Shashi",
        SIDM_KRISHNAMURTI: "Krishnamurti",
        SIDM_DJWHAL_KHUL: "Djwhal Khul",
        SIDM_YUKTESHWAR: "Yukteshwar",
        SIDM_JN_BHASIN: "J.N. Bhasin",
        SIDM_BABYL_KUGLER1: "Babylonian/Kugler 1",
        SIDM_BABYL_KUGLER2: "Babylonian/Kugler 2",
        SIDM_BABYL_KUGLER3: "Babylonian/Kugler 3",
        SIDM_BABYL_HUBER: "Babylonian/Huber",
        SIDM_BABYL_ETPSC: "Babylonian/Eta Piscium",
        SIDM_BABYL_BRITTON: "Babylonian/Britton",
        SIDM_ALDEBARAN_15TAU: "Babylonian/Aldebaran = 15 Tau",
        SIDM_TRUE_CITRA: "True Citra",
        SIDM_TRUE_REVATI: "True Revati",
        SIDM_TRUE_PUSHYA: "True Pushya (PVRN Rao)",
        SIDM_TRUE_MULA: "True Mula (Chandra Hari)",
        SIDM_TRUE_SHEORAN: '"Vedic"/Sheoran',
        SIDM_HIPPARCHOS: "Hipparchos",
        SIDM_SASSANIAN: "Sassanian",
        SIDM_J2000: "J2000",
        SIDM_J1900: "J1900",
        SIDM_B1950: "B1950",
        SIDM_SURYASIDDHANTA: "Suryasiddhanta",
        SIDM_SURYASIDDHANTA_MSUN: "Suryasiddhanta, mean Sun",
        SIDM_ARYABHATA: "Aryabhata",
        SIDM_ARYABHATA_MSUN: "Aryabhata, mean Sun",
        SIDM_ARYABHATA_522: "Aryabhata 522",
        SIDM_SS_REVATI: "SS Revati",
        SIDM_SS_CITRA: "SS Citra",
        SIDM_GALCENT_0SAG: "Galact. Center = 0 Sag",
        SIDM_GALCENT_RGILBRAND: "Galactic Center (Gil Brand)",
        SIDM_GALCENT_MULA_WILHELM: "Dhruva/Gal.Center/Mula (Wilhelm)",
        SIDM_GALCENT_COCHRANE: "Cochrane (Gal.Center = 0 Cap)",
        SIDM_GALEQU_IAU1958: "Galactic Equator (IAU1958)",
        SIDM_GALEQU_TRUE: "Galactic Equator",
        SIDM_GALEQU_MULA: "Galactic Equator mid-Mula",
        SIDM_GALEQU_FIORENZA: "Galactic Equator (Fiorenza)",
        SIDM_GALALIGN_MARDYKS: "Skydram (Mardyks)",
        SIDM_VALENS_MOON: "Vettius Valens",
        SIDM_LAHIRI_1940: "Lahiri 1940",
        SIDM_LAHIRI_VP285: "Lahiri VP285",
        SIDM_KRISHNAMURTI_VP291: "Krishnamurti-Senthilathiban",
        SIDM_LAHIRI_ICRC: "Lahiri ICRC",
    }
    # Strip the SIDBIT_* projection flags (all >= 256) so a mode combined with
    # ECL_T0/SSY_PLANE/ECL_DATE/... still resolves to its base name.
    return names.get(sidmode & 0xFF, "")


@dataclass
class StarData:
    """
    Fixed star astrometric data for ayanamsha calculations.

    All coordinates are ICRS J2000.0 epoch. Proper motion values are
    mu_alpha* (includes cos(dec) factor) and mu_delta, as in Hipparcos.

    Attributes:
        ra_j2000: Right Ascension at J2000.0 in degrees
        dec_j2000: Declination at J2000.0 in degrees
        pm_ra: Proper motion in RA (arcsec/year, includes cos(dec) factor)
        pm_dec: Proper motion in Dec (arcsec/year)
        parallax: Parallax in arcseconds (default 0.0)
        radial_velocity: Radial velocity in km/s (default 0.0)

    Note:
        For high-precision ayanamsha calculations, parallax and radial
        velocity are used in the rigorous space motion approach.
    """

    ra_j2000: float  # degrees
    dec_j2000: float  # degrees
    pm_ra: float  # arcsec/year (mu_alpha*, includes cos(dec) factor)
    pm_dec: float  # arcsec/year (mu_delta)
    parallax: float = 0.0  # arcseconds
    radial_velocity: float = 0.0  # km/s


# Star Coordinates (ICRS J2000)
# =============================================================================
# High-precision astrometric data from Hipparcos/Gaia catalogs for ayanamsha
# calculations. All values are at J2000.0 epoch in ICRS reference frame.
#
# References:
# - Hipparcos Catalogue (ESA SP-1200, 1997)
# - Gaia DR3 (where available)
# - SIMBAD Astronomical Database
#
# For Spica (Alpha Virginis, HIP 65474) - critical for True Citra ayanamsha:
#   Hipparcos values (transformed to J2000.0):
#   - RA: 13h 25m 11.5794s = 201.2982475°
#   - Dec: -11° 09' 40.759" = -11.1613219°
#   - pm_ra (mu_alpha*): -42.50 ± 0.86 mas/yr
#   - pm_dec (mu_delta): -31.73 ± 0.57 mas/yr
#   - Parallax: 12.44 ± 0.86 mas
#   - Radial velocity: +1.0 km/s
# =============================================================================
STARS = {
    # SPICA (Alpha Virginis, HIP 65474)
    # High-precision Hipparcos values for True Citra ayanamsha
    # Proper motion includes rigorous space motion corrections
    "SPICA": StarData(
        ra_j2000=201.2982475,  # 13h 25m 11.5794s (Hipparcos J2000.0)
        dec_j2000=-11.1613219,  # -11° 09' 40.759" (Hipparcos J2000.0)
        pm_ra=-0.04250,  # -42.50 mas/yr (Hipparcos mu_alpha*)
        pm_dec=-0.03173,  # -31.73 mas/yr (Hipparcos mu_delta)
        parallax=0.01244,  # 12.44 mas (Hipparcos)
        radial_velocity=1.0,  # +1.0 km/s (towards us is negative, away is positive)
    ),
    # REVATI (Zeta Piscium A, HIP 5737)
    # High-precision Gaia DR3 values for True Revati ayanamsha
    # Proper motion includes rigorous space motion corrections
    "REVATI": StarData(
        ra_j2000=18.4328583349,  # 01h 13m 43.8860s (Gaia DR3 J2000.0)
        dec_j2000=7.5753601597,  # +07° 34' 31.296" (Gaia DR3 J2000.0)
        pm_ra=0.142693,  # 142.693 mas/yr (Gaia DR3 mu_alpha*)
        pm_dec=-0.053051,  # -53.051 mas/yr (Gaia DR3 mu_delta)
        parallax=0.0244595,  # 24.4595 mas (Gaia DR3)
        radial_velocity=15.0,  # +15.0 km/s (SIMBAD - away from us)
    ),
    # ALDEBARAN (Alpha Tauri, HIP 21421) — anchor of SIDM_ALDEBARAN_15TAU
    "ALDEBARAN": StarData(
        ra_j2000=68.98016279,  # 04h 35m 55.239s (Hipparcos, van Leeuwen 2007)
        dec_j2000=16.509302351,  # +16° 30' 33.489" (Hipparcos, van Leeuwen 2007)
        pm_ra=0.06345,  # 63.45 mas/yr (Hipparcos mu_alpha*)
        pm_dec=-0.18894,  # -188.94 mas/yr (Hipparcos mu_delta)
        parallax=0.04894,  # 48.94 mas (Hipparcos)
        radial_velocity=54.26,  # +54.26 km/s (XHIP)
    ),
    # PUSHYA (Delta Cancri / Asellus Australis, HIP 42911)
    # High-precision Gaia DR3 values for True Pushya ayanamsha
    # Proper motion includes rigorous space motion corrections
    "PUSHYA": StarData(
        ra_j2000=131.1712460977,  # 08h 44m 41.0991810454s (Gaia DR3 J2000.0)
        dec_j2000=18.1543080691,  # +18° 09' 15.509048595" (Gaia DR3 J2000.0)
        pm_ra=-0.018435,  # -18.435 mas/yr (Gaia DR3 mu_alpha*)
        pm_dec=-0.227813,  # -227.813 mas/yr (Gaia DR3 mu_delta)
        parallax=0.0238271,  # 23.8271 mas (Gaia DR3)
        radial_velocity=17.14,  # +17.14 km/s (SIMBAD - away from us)
    ),
    # MULA (Lambda Scorpii / Shaula, HIP 85927)
    # High-precision Hipparcos 2007 values for True Mula ayanamsha
    # (Gaia DR3 unavailable - star too bright at V=1.63)
    # Proper motion includes rigorous space motion corrections
    "MULA": StarData(
        ra_j2000=263.40216717,  # 17h 33m 36.52012s (Hipparcos 2007 J2000.0)
        dec_j2000=-37.10382356,  # -37° 06' 13.7648" (Hipparcos 2007 J2000.0)
        pm_ra=-0.00853,  # -8.53 mas/yr (Hipparcos 2007 mu_alpha*)
        pm_dec=-0.03080,  # -30.80 mas/yr (Hipparcos 2007 mu_delta)
        parallax=0.00571,  # 5.71 mas (Hipparcos 2007)
        radial_velocity=-3.00,  # -3.00 km/s (SIMBAD - approaching us)
    ),
    # Galactic Center (Sgr A*) - the supermassive black hole at the center of the
    # Milky Way, detected as a compact radio source.
    #
    # ICRS J2000 position from Reid & Brunthaler (2004, ApJ 616, 872):
    #   RA  = 17h 45m 40.0409s = 266.41683708 deg
    #   Dec = -29° 00' 28.118" = -29.00781056 deg
    #   Epoch: J2000.0 (based on VLBA observations)
    #   Uncertainty: ~0.4 mas in each coordinate
    #
    # Apparent proper motion from Reid & Brunthaler (2020, ApJ 892, 39):
    #   μl = -6.411 ± 0.008 mas/yr (along Galactic plane)
    #   μb = -0.219 ± 0.007 mas/yr (toward North Galactic Pole)
    #
    # Conversion from Galactic to equatorial proper motion:
    #   At Sgr A* position (l = 359.944°, b = -0.046°):
    #   The position angle of Galactic North from celestial North is ~58.3°
    #   Using the transformation matrix:
    #     μα* = μl*sin(posang) + μb*cos(posang) = -6.411*sin(58.3°) + -0.219*cos(58.3°)
    #         = -6.411*0.8507 - 0.219*0.5257 = -5.456 - 0.115 = -5.571 mas/yr
    #   But note: the sign conventions differ; verified against:
    #     μα* ≈ -2.70 mas/yr (Wikipedia, from Reid & Brunthaler papers)
    #     μδ  ≈ -5.6 mas/yr
    #   Using the equatorial proper motions from direct VLBI observation:
    #     μα* = -3.151 ± 0.018 mas/yr (Reid & Brunthaler 2004)
    #     μδ  = -5.547 ± 0.026 mas/yr (Reid & Brunthaler 2004)
    #
    # Note: This apparent motion is dominated by the Sun's orbital motion around
    # the Galactic center (galactic rotation + solar peculiar velocity).
    # The intrinsic motion of Sgr A* is effectively zero (<1 km/s).
    #
    # References:
    # - Reid, M. J. & Brunthaler, A. 2004, ApJ, 616, 872 (position)
    # - Reid, M. J. & Brunthaler, A. 2020, ApJ, 892, 39 (proper motion update)
    # - Petrov et al. 2011, AJ, 142, 35 (SIMBAD reference, consistent position)
    "GAL_CENTER": StarData(
        ra_j2000=266.41683708,  # 17h 45m 40.0409s (Reid & Brunthaler 2004)
        dec_j2000=-29.00781056,  # -29° 00' 28.118" (Reid & Brunthaler 2004)
        # pm_ra must be the great-circle rate mu_alpha* (= d(alpha)/dt * cos dec),
        # because both StarData.pm_ra and Skyfield's ``ra_mas_per_year`` expect
        # the on-sky angular rate (Skyfield divides it back by cos dec to get the
        # coordinate rate internally). Reid & Brunthaler's tabulated -3.151 mas/yr
        # for Sgr A* is the RA-*coordinate* rate d(alpha)/dt (no cos dec), so the
        # great-circle rate is
        #   mu_alpha* = -3.151 * cos(-29.00781 deg) = -3.151 * 0.874554 = -2.7557 mas/yr.
        # Feeding the raw -3.151 as if it were mu_alpha* over-states the GC's
        # apparent RA motion and drifts the GC ayanamsha up to +3.4" at year 1200.
        # This cos(dec) conversion is needed for GAL_CENTER only: the other
        # stars store catalog mu_alpha* values that already include cos(dec).
        pm_ra=-0.003151 * math.cos(math.radians(-29.00781056)),
        pm_dec=-0.005547,  # -5.547 mas/yr -> arcsec/yr (Reid & Brunthaler 2004)
    ),
    # Galactic North Pole (J2000)
    "GAL_NORTH_POLE": StarData(
        ra_j2000=192.85948,
        dec_j2000=27.12825,
        pm_ra=0.0,
        pm_dec=0.0,
    ),
}


def _iau2006_general_precession_deg(jd_tt: float) -> float:
    """Return ERFA's IAU 2006 general precession p_A since J2000."""
    angles = erfa.p06e(2451545.0, jd_tt - 2451545.0)
    return math.degrees(float(angles[12]))


def _galcent_ayanamsha(
    tjd_tt: float,
    target_lon: float,
    noaberr: bool = False,
    nogdefl: bool = False,
) -> float:
    """Galactic-Center mode: Sgr A* held at its published sidereal longitude.

    The ayanamsha is the apparent mean-ecliptic longitude of Sgr A*
    (Reid & Brunthaler 2004 ICRS position and proper motion) minus the
    mode's published target longitude (e.g. 240° for "Galactic Center at
    0° Sagittarius").

    Under ``noaberr`` (FLG_TRUEPOS / FLG_NOABERR) the anchor's annual
    aberration is removed from the zero point, and under ``nogdefl``
    (FLG_TRUEPOS / FLG_NOGDEFL) its solar light deflection is removed;
    measured reference behavior.
    """
    longitude = _get_star_position_ecliptic(
        STARS["GAL_CENTER"],
        tjd_tt,
        0.0,
        nonut=True,
        noaberr=noaberr,
        nogdefl=nogdefl,
    )
    return (longitude - target_lon) % 360.0


def _galcent_ra_of_date(
    tjd_tt: float, noaberr: bool = False, nogdefl: bool = False
) -> float:
    """Return apparent Sgr A* RA on the Vondrak mean equator of date.

    ``noaberr`` removes only the observer-velocity aberration and ``nogdefl``
    removes only the Sun's gravitational light deflection, so the four flag
    combinations match the FLG_TRUEPOS / FLG_NOABERR / FLG_NOGDEFL anchor
    convention (FLG_TRUEPOS maps to both, i.e. the geometric direction).
    """
    import numpy as np

    center = STARS["GAL_CENTER"]
    star = Star(
        ra_hours=center.ra_j2000 / 15.0,
        dec_degrees=center.dec_j2000,
        ra_mas_per_year=center.pm_ra * 1000.0,
        dec_mas_per_year=center.pm_dec * 1000.0,
    )
    time = get_timescale().tt_jd(tjd_tt)
    try:
        astrometric = _get_computation_ephemeris()["earth"].at(time).observe(star)
        astro_au = np.array(astrometric.position.au, dtype=float)
        if noaberr and nogdefl:
            # Geometric direction (FLG_TRUEPOS): no deflection, no aberration.
            vector = astro_au
        elif nogdefl:
            # Aberration only: displace the astrometric direction by aberration.
            vector = astro_au.copy()
            add_aberration(
                vector,
                astrometric.center_barycentric.velocity.au_per_d,
                astrometric.light_time,
            )
        elif noaberr:
            # Deflection but no aberration: subtract apparent()'s final
            # aberration step (see _star_position_ecliptic_uncached).
            apparent_au = np.array(astrometric.apparent().position.au, dtype=float)
            aberrated_au = astro_au.copy()
            add_aberration(
                aberrated_au,
                astrometric.center_barycentric.velocity.au_per_d,
                astrometric.light_time,
            )
            vector = apparent_au - (aberrated_au - astro_au)
        else:
            vector = np.array(astrometric.apparent().position.au)
    except _RANGE_ERRORS as error:
        raise _wrap_ephemeris_range_error(error, tjd_tt) from error
    x_coord, y_coord, _ = np.array(vondrak_precession_matrix(tjd_tt)) @ vector
    return math.atan2(float(y_coord), float(x_coord))


def _mula_wilhelm_longitude(
    tjd_tt: float, noaberr: bool = False, nogdefl: bool = False
) -> float:
    """Project Sgr A*'s hour circle onto the mean ecliptic of date."""
    alpha = _galcent_ra_of_date(tjd_tt, noaberr=noaberr, nogdefl=nogdefl)
    epsilon = vondrak_mean_obliquity_rad(tjd_tt)
    return (
        math.degrees(math.atan2(math.sin(alpha), math.cos(alpha) * math.cos(epsilon)))
        % 360.0
    )


def _get_star_position_ecliptic(
    star: StarData,
    tjd_tt: float,
    eps_true: float,
    nonut: bool = False,
    geometric: bool = False,
    noaberr: bool = False,
    nogdefl: bool = False,
) -> float:
    """Interpolating wrapper: see _star_position_ecliptic_uncached.

    The anchor longitude varies slowly (precession ~0.14"/day, proper
    motion <0.01"/day, annual aberration <0.36"/day), so it is evaluated
    exactly only at memoized 8-day grid nodes and interpolated with a
    Lagrange quadratic in between. The fastest interpolated term (annual
    aberration, ~20.5" sinusoid) leaves a cubic residual below 0.004" over a
    16-day node span — negligible against the pipeline's own precision — while
    cutting the per-call cost of the star-based ayanamsha modes by the
    node-reuse factor (the full Skyfield pipeline runs only at nodes).

    The Sun's gravitational light deflection is the exception: near the annual
    conjunction it changes sign across the source and can reach ~1" over a
    single day for a near-ecliptic anchor, a feature the 8-day node grid
    cannot represent. When deflection is applied (neither ``geometric`` nor
    ``nogdefl``) the deflected apparent place is therefore evaluated *exactly*
    at ``tjd_tt`` rather than interpolated; the per-date LRU cache still
    collapses the many-bodies-at-one-date chart workload to a single pipeline
    run. Deflection-free longitudes (``geometric`` frame poles, or FLG_NOGDEFL
    / FLG_TRUEPOS) stay on the fast interpolation path, since without the
    deflection spike the longitude is a slow, smooth function of date.
    """
    # Under nonut the true obliquity argument is unused: normalize it so
    # the node cache key does not vary with the caller's request date.
    eps_key = 0.0 if nonut else eps_true

    if not geometric and not nogdefl:
        # Deflection is active and cannot be interpolated: evaluate exactly.
        try:
            return _star_position_ecliptic_cached(
                star.ra_j2000,
                star.dec_j2000,
                star.pm_ra,
                star.pm_dec,
                star.parallax,
                star.radial_velocity,
                tjd_tt,
                eps_key,
                nonut,
                geometric,
                noaberr,
                nogdefl,
            )
        except _RANGE_ERRORS as e:
            # Translate the raw skyfield range error to our typed
            # EphemerisRangeError, matching the interpolation path below so the
            # star/galactic ayanamsha modes surface the same exception as
            # calc_ut/fixstar for out-of-tier epochs.
            raise _wrap_ephemeris_range_error(e, tjd_tt) from e

    h = 8.0  # node spacing in days
    base = math.floor(tjd_tt / (2.0 * h)) * (2.0 * h)

    def _node(jd_node: float) -> float:
        return _star_position_ecliptic_cached(
            star.ra_j2000,
            star.dec_j2000,
            star.pm_ra,
            star.pm_dec,
            star.parallax,
            star.radial_velocity,
            jd_node,
            eps_key,
            nonut,
            geometric,
            noaberr,
            nogdefl,
        )

    try:
        y0, y1, y2 = _node(base), _node(base + h), _node(base + 2.0 * h)
    except _RANGE_ERRORS as e:
        # The Skyfield star pipeline (earth.at().observe()) raises the raw
        # skyfield range error when the epoch is outside the active tier.
        # Translate to our typed EphemerisRangeError so the star/galactic
        # ayanamsha modes surface the same exception as calc_ut/fixstar.
        raise _wrap_ephemeris_range_error(e, tjd_tt) from e
    # Unwrap around y1 so the quadratic never crosses the 0/360 seam
    y0 += (y1 - y0 + 180.0) % 360.0 - 180.0 - (y1 - y0)
    y2 += (y1 - y2 + 180.0) % 360.0 - 180.0 - (y1 - y2)
    x = (tjd_tt - base) / h  # in [0, 2]
    lon = (
        y0 * (x - 1.0) * (x - 2.0) / 2.0 - y1 * x * (x - 2.0) + y2 * x * (x - 1.0) / 2.0
    )
    return lon % 360.0


@lru_cache(maxsize=1024)
def _star_position_ecliptic_cached(
    ra_j2000: float,
    dec_j2000: float,
    pm_ra: float,
    pm_dec: float,
    parallax: float,
    radial_velocity: float,
    tjd_tt: float,
    eps_true: float,
    nonut: bool,
    geometric: bool,
    noaberr: bool,
    nogdefl: bool,
) -> float:
    """Cache the scalarized fixed-star longitude calculation by full input.

    Scalar fields make the key immutable; the uncached routine retains the
    catalogue-motion and frame provenance documented for fixed stars.
    ``noaberr`` and ``nogdefl`` are both part of the key so the aberration-free
    (FLG_NOABERR), deflection-free (FLG_NOGDEFL) and geometric (FLG_TRUEPOS)
    anchors never collide with each other or with the default apparent one.
    """
    star = StarData(
        ra_j2000=ra_j2000,
        dec_j2000=dec_j2000,
        pm_ra=pm_ra,
        pm_dec=pm_dec,
        parallax=parallax,
        radial_velocity=radial_velocity,
    )
    return _star_position_ecliptic_uncached(
        star,
        tjd_tt,
        eps_true,
        nonut=nonut,
        geometric=geometric,
        noaberr=noaberr,
        nogdefl=nogdefl,
    )


def _star_position_ecliptic_uncached(
    star: StarData,
    tjd_tt: float,
    eps_true: float,
    nonut: bool = False,
    geometric: bool = False,
    noaberr: bool = False,
    nogdefl: bool = False,
) -> float:
    """
    Calculate ecliptic longitude of a fixed star at given date.

    Uses the Skyfield astrometric pipeline for the apparent ICRS direction:
      - Proper motion propagation (rigorous space motion with radial velocity)
      - Annual aberration (~20.5" correction) + light deflection

    then rotates that apparent ICRS vector to the **Vondrák 2011 ecliptic of
    date** (the same frame `_calc_body` uses), rather than Skyfield's IAU 2006
    `ecliptic_frame`. This keeps the star-based ayanamsha modes (True Citra /
    Revati / Galactic) on the library's precession+obliquity model so they stay
    consistent with the planet longitudes at remote epochs. Nutation is IAU
    2006/2000A (via `erfa.nut06a`).

    Args:
        star: Star catalog data (J2000.0 ICRS coordinates, proper motion, parallax, radial velocity)
        tjd_tt: Julian Day in Terrestrial Time (TT)
        eps_true: True obliquity of date in degrees (Vondrák mean + nutation),
                  used for the equator-of-date -> ecliptic-of-date rotation.
        nonut: When True, return the longitude on the MEAN ecliptic of date
            (precession-only rotation, mean obliquity; ``eps_true`` is
            ignored). Aberration and light deflection stay applied. Star-based
            mean ayanamshas use this NONUT longitude; true ayanamshas add
            nutation in longitude exactly once.
        geometric: When True, skip aberration and light deflection (use the
            astrometric direction). Required for frame poles, which are
            geometric directions rather than light sources.
        noaberr: When True, remove only the annual aberration of the anchor
            while keeping light deflection (the astrometric-place convention of
            FLG_NOABERR). Ignored when ``geometric`` is set, since a geometric
            direction already carries no aberration. Used by the star/galactic
            ayanamsha modes under FLG_TRUEPOS / FLG_NOABERR (measured reference
            behavior removes the anchor's aberration from the zero point).
        nogdefl: When True, remove the Sun's gravitational light deflection of
            the anchor while keeping annual aberration (the FLG_NOGDEFL
            convention; general relativity, Einstein 1916 — the same reduction
            the fixed-star pipeline applies). ``noaberr`` and ``nogdefl`` are
            independent: FLG_NOABERR alone keeps deflection, FLG_NOGDEFL alone
            keeps aberration, and FLG_TRUEPOS (which maps to both being True)
            yields the geometric direction. Ignored when ``geometric`` is set.

    Returns:
        Ecliptic longitude of date in degrees (0-360)

    References:
        - Skyfield: Brandon Rhodes, skyfield.readthedocs.io
        - Vondrák, Capitaine & Wallace, A&A 534, A22 (2011) — long-term precession
        - IAU 2000A nutation: Mathews, Herring & Buffett, JGR 107 (2002)
    """
    # Convert StarData to Skyfield Star object
    # StarData uses arcsec/yr; Skyfield uses mas/yr
    # StarData uses arcsec for parallax; Skyfield uses mas
    star_obj = Star(
        ra_hours=star.ra_j2000 / 15.0,
        dec_degrees=star.dec_j2000,
        ra_mas_per_year=star.pm_ra * 1000.0,
        dec_mas_per_year=star.pm_dec * 1000.0,
        parallax_mas=star.parallax * 1000.0 if star.parallax > 0 else 0.0,
        radial_km_per_s=star.radial_velocity,
    )

    # Use Skyfield pipeline: observe -> apparent (includes aberration)
    # then transform to ecliptic of date (includes precession + nutation)
    from .state import _get_computation_ephemeris

    planets = _get_computation_ephemeris()
    ts = get_timescale()
    t = ts.tt_jd(tjd_tt)
    earth = planets["earth"]

    import numpy as np

    # earth.at(t).observe(star) applies proper motion + light-time; .apparent()
    # then applies gravitational light deflection followed by annual aberration
    # as its final step.
    astrometric = earth.at(t).observe(star_obj)
    astro_au = np.array(astrometric.position.au, dtype=float)
    if geometric or (noaberr and nogdefl):
        # Geometric (astrometric) direction: no deflection, no aberration.
        # FLG_TRUEPOS (and FLG_NOABERR | FLG_NOGDEFL together) request exactly
        # this place; measured reference behavior removes both reductions.
        position_au = astro_au
    elif nogdefl:
        # Aberration only (FLG_NOGDEFL keeps annual aberration but drops the
        # Sun's gravitational light deflection): displace the astrometric
        # direction by aberration alone (Bradley 1728; IAU aberration constant
        # kappa ~ 20.49552"), leaving the geometric direction otherwise intact.
        position_au = astro_au.copy()
        add_aberration(
            position_au,
            astrometric.center_barycentric.velocity.au_per_d,
            astrometric.light_time,
        )
    elif noaberr:
        # Deflection only (FLG_NOABERR keeps the Sun's gravitational light
        # deflection but drops the observer-velocity aberration). apparent()
        # applies deflection first and aberration last, so subtract that final
        # displacement. Evaluating the aberration displacement on the
        # astrometric direction rather than on the deflected one differs only at
        # second order in (deflection x beta), i.e. well below a micro-arcsecond.
        apparent_au = np.array(astrometric.apparent().position.au, dtype=float)
        aberrated_au = astro_au.copy()
        add_aberration(
            aberrated_au,
            astrometric.center_barycentric.velocity.au_per_d,
            astrometric.light_time,
        )
        position_au = apparent_au - (aberrated_au - astro_au)
    else:
        # Full apparent place: gravitational light deflection + annual
        # aberration (the default apparent anchor; general relativity plus
        # Bradley aberration, matching the fixed-star pipeline).
        position_au = np.array(astrometric.apparent().position.au, dtype=float)

    # Rotate the apparent ICRS direction to the Vondrák ecliptic of date (same
    # reduction as _calc_body), not Skyfield's IAU 2006 ecliptic_frame.
    if nonut:
        # Mean ecliptic of date: precession-only rotation, mean obliquity.
        pn = vondrak_precession_matrix(tjd_tt)
        eps_rad = vondrak_mean_obliquity_rad(tjd_tt)
    else:
        dpsi_rad, deps_rad = erfa.nut06a(2451545.0, tjd_tt - 2451545.0)
        pn, _eps = vondrak_pn_matrix(tjd_tt, float(dpsi_rad), float(deps_rad))
        eps_rad = math.radians(eps_true)
    xq, yq, zq = np.array(pn) @ position_au
    ce, se = math.cos(eps_rad), math.sin(eps_rad)
    xe = float(xq)
    ye = float(yq) * ce + float(zq) * se
    return math.degrees(math.atan2(ye, xe)) % 360.0


def _ecl_date_ayanamsha_delta(
    defining_value_deg: float, t0_tt: float, tjd_tt: float
) -> float:
    """Longitude shift (deg) for referring a mode to the ecliptic of date.

    This is the geometric realization of SIDBIT_ECL_DATE, applied by
    ``_calc_ayanamsa`` to every defining pair (epoch modes, Valens and
    SIDM_USER); the LEB fast path delegates here when the bit is set so both
    backends stay numerically identical. Live star/galactic modes are
    naturally inert: their value is already an of-date longitude.

    SIDBIT_ECL_DATE refers the sidereal longitude to the mean ecliptic and
    equinox **of date** instead of the mode's defining mean ecliptic of ``t0``
    (the default Method-B frame; see
    :func:`precession_vondrak.method_b_accumulated_precession`). Geometrically
    the fixed sidereal zero point ``Z`` is defined on the mean ecliptic of
    ``t0`` at ecliptic longitude ``defining_value_deg`` (latitude 0). This
    returns the difference between the tropical longitude of ``Z`` measured on
    the mean ecliptic of date and on the mean ecliptic of ``t0`` (the Method-B
    baseline), i.e. the amount to add to the Method-B ayanamsha.

    Both ecliptic frames are built from the Vondrák 2011 long-term
    ecliptic/equator poles (``erfa.ltpecl``/``ltpequ`` via
    :func:`precession_vondrak._ltp_ecliptic_frame`), so the construction is
    valid over the full Vondrák span and introduces no fitted term. The shift
    is zero at ``tjd_tt == t0_tt`` and grows with the ecliptic tilt between the
    two epochs (a few milli-arcseconds for modern epochs, up to ~0.5-0.7" for
    modes whose defining epoch is a millennium or more from the request).

    Args:
        defining_value_deg: The mode's ayanamsha at ``t0`` (its sidereal-zero
            longitude on the mean ecliptic of ``t0``), in degrees.
        t0_tt: The mode's defining epoch, Julian Date in TT.
        tjd_tt: The epoch of interest, Julian Date in TT.

    Returns:
        Signed longitude shift in degrees, unwrapped to ``(-180, 180]``.
    """
    from .precession_vondrak import _ltp_ecliptic_frame

    # Sidereal zero point Z on the mean ecliptic of t0 (longitude defining
    # value, latitude 0), expressed in the GCRS basis: (ecliptic-of-t0
    # basis rows)^T @ [cos, sin, 0].
    e0 = _ltp_ecliptic_frame(t0_tt)
    lon0 = math.radians(defining_value_deg)
    cx, sx = math.cos(lon0), math.sin(lon0)
    zx = e0[0][0] * cx + e0[1][0] * sx
    zy = e0[0][1] * cx + e0[1][1] * sx
    zz = e0[0][2] * cx + e0[1][2] * sx

    # Tropical longitude of Z on the mean ecliptic and equinox of date.
    ed = _ltp_ecliptic_frame(tjd_tt)
    wx = ed[0][0] * zx + ed[0][1] * zy + ed[0][2] * zz
    wy = ed[1][0] * zx + ed[1][1] * zy + ed[1][2] * zz
    lon_date = math.degrees(math.atan2(wy, wx))

    method_b = defining_value_deg + method_b_accumulated_precession(tjd_tt, t0_tt)
    return (lon_date - method_b + 180.0) % 360.0 - 180.0


def _ecl_t0_zero_point_deg(sid_mode: int, t0_jd: float) -> float:
    """Sidereal zero point on the mean ecliptic of ``t0_jd``.

    A single evaluation for every mode, so the result never depends on which
    precision tier happens to be installed. An earlier revision fell back to
    a catalog-only anchor when the plane epoch lay outside the loaded kernel;
    that made the same public call return two different longitudes 0.79 deg
    apart on the medium and extended tiers, which is a worse defect than the
    divergence it was meant to close.
    """
    return _calc_ayanamsa(t0_jd, sid_mode)


def _ecl_t0_epoch_jd(sid_mode: int) -> float:
    """Defining epoch (JD TT) of a sidereal mode for the ECL_T0 projection.

    Valens (SIDM_VALENS_MOON) keeps its defining pair outside
    AYANAMSHA_DEFINING because its epoch is UT-anchored; without this
    special case the lookup silently fell back to J2000 and the ECL_T0
    projection became a near no-op (up to ~4' of error at high ecliptic
    latitudes) instead of using the mode's ancient defining ecliptic.

    SIDM_USER takes the epoch from the user's ``set_sid_mode(SIDM_USER, t0,
    ayan_t0)`` call, NOT from the defining table (whose lookup would fall
    back to J2000 and freeze the projection plane). Measured reference
    behavior: the plane is the mean ecliptic of the stored t0 taken
    literally (t0 = 0.0 really means JD 0), the latitude shift vanishes
    when t0 equals the computation date and grows by ~0.36 arcsec per year
    of |date - t0|, and it is independent of ayan_t0 (rotating the zero
    point within the plane moves longitudes, never latitudes). With
    SIDBIT_USER_UT the stored t0 is a UT date and is converted to TT here,
    mirroring the _calc_ayanamsa anchoring convention.
    """
    if sid_mode == SIDM_USER:
        _, t0, _ = get_sid_mode(full=True)
        if _get_sidereal_bits() & SIDBIT_USER_UT:
            from .time_utils import deltat

            t0 = t0 + deltat(t0)
        return t0
    if sid_mode == SIDM_VALENS_MOON:
        from .time_utils import deltat

        return VALENS_MOON_T0_UT + deltat(VALENS_MOON_T0_UT)
    return _ECL_T0_CLASSICAL_EPOCHS.get(
        sid_mode, AYANAMSHA_DEFINING.get(sid_mode, (0.0, _J2000_JD))[1]
    )


# Classical defining epochs (TT JD) for the ECL_T0 projection plane, used
# where a mode's classical defining epoch — the mean ecliptic and equinox its
# authors referred the sidereal zero to — differs from the epoch of the mode's
# VALUE anchor in AYANAMSHA_DEFINING. The SIDBIT_ECL_T0 projection plane must
# follow the classical defining epoch, not the value anchor; resolving it from
# the value anchor mis-orients the mean equator/ecliptic (up to ~29° in RA and
# tens of arcsec in ecliptic latitude for these modes). Each entry cites the
# published defining statement of the mode and resolves to a round conventional
# epoch of that mode's tradition (the beginning of the Christian era = year
# 0.0, the standard modern epoch 1900.0, or the Indian Calendar Reform anchor).
#
# LAHIRI: the value realization deliberately follows the Indian Astronomical
# Ephemeris tabulation anchored at J2000 (see ayanamsha_definitions.py), but
# the mode's classical defining epoch — and therefore the ecliptic the
# ECL_T0 projection refers to — is the Calendar Reform Committee anchor:
# mean ayanamsha 23°15' at the vernal equinox of 1956 (Report of the
# Calendar Reform Committee, CSIR, New Delhi, 1955; restated in the Indian
# Astronomical Ephemeris introduction). JD 2435553.5 = 1956 March 21.0.
_ECL_T0_CLASSICAL_EPOCHS: dict[int, float] = {
    SIDM_LAHIRI: 2435553.5,
    # BABYL_KUGLER1/2/3: the three Kugler solutions (Sternkunde und
    # Sterndienst in Babel II) are alternative VALUE readings of one and the
    # same Babylonian sidereal zodiac; the classical norm epoch of that
    # zodiac is ~-100 (P. Huber, "Ueber den Nullpunkt der babylonischen
    # Ekliptik", Centaurus 5, 1958; cf. Fagan, Zodiacs Old and New, 1950),
    # so the whole family shares the mean ecliptic of year -100 as its
    # projection plane. JD 1684532.5 = -100 Jan 1.0 (Julian). Measured
    # reference behavior confirms one shared plane for the triplet.
    SIDM_BABYL_KUGLER1: 1684532.5,
    SIDM_BABYL_KUGLER2: 1684532.5,
    SIDM_BABYL_KUGLER3: 1684532.5,
    # HIPPARCHOS: the Hipparchan sidereal norm (R. Mercier, "Studies in the
    # Medieval Conception of Precession") is anchored to Hipparchus's own
    # observational era, ~-128 (the epoch of his catalog and equinox
    # observations quoted by Ptolemy, Almagest III); the projection plane is
    # the mean ecliptic of that era. JD 1674484.0 (~-128.5 Julian).
    SIDM_HIPPARCHOS: 1674484.0,
    # DELUCE: Robert De Luce, "Constellational Astrology According to the
    # Hindu System" (De Luce Publishing Co., Los Angeles, 1963) fixes the
    # coincidence of the constellational and sign zodiacs (ayanamsha zero) at
    # the beginning of the Christian era; the defining plane is therefore the
    # mean ecliptic of year 0.0. JD 1721057.5 = 0000 Jan 1.0 (Julian).
    SIDM_DELUCE: 1721057.5,
    # DJWHAL_KHUL: the Bailey/Djwhal Khul Aquarian-age doctrine (channelled
    # 1940) sets the start of the Age of Aquarius at 2117 (ayanamsha 30°); its
    # astronomical realization (Phillip Lindsay, "The Beginning of the Age of
    # Aquarius", 2006) specifies the value via the Synetic Vernal Point at the
    # standard epoch 1900.0, which is the defining plane epoch.
    # JD 2415020.0 = 1900 Jan 0.5 (1900.0).
    SIDM_DJWHAL_KHUL: 2415020.0,
    # BABYL_BRITTON: J.P. Britton, "Studies in Babylonian lunar theory: Part
    # III. The introduction of the uniform zodiac" (Archive for History of
    # Exact Sciences 64, 2010) expresses the Babylonian-tropical displacement
    # as Δλ* = C − 1.3828°·Y with C = 3.20° and Y in centuries from the
    # beginning of the common era; the constant term is defined at Y = 0, so
    # year 0.0 is the classical defining plane epoch (the value zero-crossing
    # implied by the formula, ~231 CE, is the separate value anchor).
    # JD 1721057.5 = 0000 Jan 1.0 (Julian).
    SIDM_BABYL_BRITTON: 1721057.5,
    # LAHIRI_1940: N.C. Lahiri, "Indian Ephemeris of Planets' Positions"
    # (1st ed. 1939/1940) — the early Lahiri tradition, before the 1948/1955
    # Calendar Reform Committee standardization on the 285 CE zero. Like the
    # contemporaneous Raman and Krishnamurti Indian ayanamshas it is tabulated
    # at the standard epoch 1900.0, which is the defining plane epoch.
    # JD 2415020.0 = 1900 Jan 0.5 (1900.0).
    SIDM_LAHIRI_1940: 2415020.0,
}


# Star / galactic-center ayanamsha modes whose sidereal zero point is anchored
# to a live light source, so its apparent place carries both annual aberration
# and the Sun's gravitational light deflection. The default anchor applies both
# reductions (Bradley 1728 aberration, IAU constant kappa ~ 20.49552"; Einstein
# 1916 deflection, up to ~1" near the annual conjunction and sign-changing
# across it — the same reduction the fixed-star pipeline applies). Under
# FLG_NOABERR the aberration is removed, under FLG_NOGDEFL the deflection is
# removed, and under FLG_TRUEPOS both are (the geometric anchor); measured
# reference behavior. This deliberately EXCLUDES the fixed-catalogue Aldebaran
# anchor (SIDM_ALDEBARAN_15TAU: the reference reports no aberration and no
# deflection for it — see the Aldebaran branch), the geometric galactic-pole
# GALEQU modes (already aberration- and deflection-free) and the
# precession-formula modes (Mardyks, Valens) — the reference reports a zero
# TRUEPOS/NOABERR/NOGDEFL shift for all of those.
_ABERRANT_ANCHOR_MODES = frozenset(
    {
        SIDM_TRUE_CITRA,
        SIDM_TRUE_REVATI,
        SIDM_TRUE_PUSHYA,
        SIDM_TRUE_MULA,
        SIDM_TRUE_SHEORAN,
        SIDM_GALCENT_0SAG,
        SIDM_GALCENT_RGILBRAND,
        SIDM_GALCENT_COCHRANE,
        SIDM_GALCENT_MULA_WILHELM,
    }
)


def _calc_ayanamsa(
    tjd_ut: float,
    sid_mode: int,
    noaberr: bool = False,
    nogdefl: bool = False,
    sid_bits: Optional[int] = None,
    sid_t0: Optional[float] = None,
    sid_ayan_t0: Optional[float] = None,
) -> float:
    """Calculate the selected predefined mean ayanamsha.

    Epoch modes share the defining epochs in :mod:`ayanamsha_definitions` and
    accrue Vondrak Method-B precession on their defining mean ecliptic.  True
    star and galactic modes are evaluated from the live independent catalog
    geometry.  All public predefined IDs 0--46 therefore remain operational;
    only an actually unknown ID warns and falls back to Lahiri.

    Args:
        tjd_ut: Julian Day in Universal Time (UT1).
        sid_mode: Sidereal mode constant, optionally combined with SIDBIT flags.
        noaberr: When True (request carried FLG_TRUEPOS or FLG_NOABERR), remove
            the anchor's annual aberration from the zero point for the modes in
            ``_ABERRANT_ANCHOR_MODES``. Inert for every other mode.
        nogdefl: When True (request carried FLG_TRUEPOS or FLG_NOGDEFL), remove
            the anchor's solar gravitational light deflection from the zero
            point for the modes in ``_ABERRANT_ANCHOR_MODES``. Inert for every
            other mode. FLG_TRUEPOS sets both toggles, yielding the geometric
            anchor.
        sid_bits: SIDBIT projection flags to apply. ``None`` reads the active
            ones (``_get_sidereal_bits()``), which is what the module-level
            entry points want; callers holding an explicit sidereal
            configuration (the LEB reducer, a delegated context call) pass
            their own so the projection cannot be taken from unrelated global
            state.
        sid_t0: SIDM_USER reference epoch (JD). ``None`` reads the active
            one. Same reasoning as ``sid_bits``: a caller with an explicit
            sidereal configuration must supply its own anchor, otherwise an
            unrelated module-level ``set_sid_mode()`` moves its result (a
            sealed-LEB context measured 19.03 degrees of drift).
        sid_ayan_t0: SIDM_USER ayanamsha value at ``sid_t0`` (degrees), with
            the same ``None`` semantics.

    Returns:
        Mean ayanamsha in degrees, normalized to [0, 360).
    """
    bits = _get_sidereal_bits() if sid_bits is None else sid_bits
    # Strip positive SIDBIT projection flags (>= 256) to recover the base mode,
    # but never mask a negative mode: -1 & 0xFF == 255 (SIDM_USER) would turn
    # an invalid ID into a valid one. Negatives fall through to the invalid-mode
    # fallback below.
    if sid_mode >= 0:
        sid_mode &= 0xFF
    # The aberration-free / deflection-free anchor only applies to the live
    # star / galactic-center modes whose reference realization removes them;
    # every other mode is unaffected regardless of the request flags.
    eff_noaberr = noaberr and sid_mode in _ABERRANT_ANCHOR_MODES
    eff_nogdefl = nogdefl and sid_mode in _ABERRANT_ANCHOR_MODES
    tjd_tt = float(get_timescale().ut1_jd(tjd_ut).tt)

    if sid_mode == SIDM_USER:
        if sid_t0 is None or sid_ayan_t0 is None:
            _, t0, ayan_t0 = get_sid_mode(full=True)
        else:
            t0, ayan_t0 = sid_t0, sid_ayan_t0
        if bits & SIDBIT_USER_UT:
            # SIDBIT_USER_UT: the user supplied t0 as a UT date. Convert it to
            # TT (t0 is used below as a TT epoch for the accumulated-precession
            # baseline), so the ayanamsha is anchored on a single timescale.
            # deltat() returns Delta T in days.
            from .time_utils import deltat

            t0 = t0 + deltat(t0)
        value = ayan_t0 + method_b_accumulated_precession(tjd_tt, t0)
        if bits & SIDBIT_ECL_DATE:
            # Refer the fixed zero point to the mean ecliptic of date instead
            # of the mean ecliptic of t0 (see _ecl_date_ayanamsha_delta).
            value += _ecl_date_ayanamsha_delta(ayan_t0, t0, tjd_tt)
        return float(value % 360.0)

    if sid_mode in DYNAMIC_AYANAMSHA_MODES:
        if sid_mode == SIDM_TRUE_CITRA:
            value = (
                _get_star_position_ecliptic(
                    STARS["SPICA"],
                    tjd_tt,
                    0.0,
                    nonut=True,
                    noaberr=eff_noaberr,
                    nogdefl=eff_nogdefl,
                )
                - 180.0
            )
        elif sid_mode == SIDM_TRUE_REVATI:
            # Usha & Shashi place Revati (zeta Piscium) at 29 deg 50 min
            # Pisces; the ayanamsha is therefore star longitude + 10 arcmin.
            value = (
                _get_star_position_ecliptic(
                    STARS["REVATI"],
                    tjd_tt,
                    0.0,
                    nonut=True,
                    noaberr=eff_noaberr,
                    nogdefl=eff_nogdefl,
                )
                + 10.0 / 60.0
            )
        elif sid_mode == SIDM_TRUE_PUSHYA:
            # P.V.R. Narasimha Rao's definition places delta Cancri at
            # 16 degrees Cancer (106 degrees absolute longitude).
            value = (
                _get_star_position_ecliptic(
                    STARS["PUSHYA"],
                    tjd_tt,
                    0.0,
                    nonut=True,
                    noaberr=eff_noaberr,
                    nogdefl=eff_nogdefl,
                )
                - 106.0
            )
        elif sid_mode == SIDM_TRUE_MULA:
            # Chandra Hari's Mula origin places lambda Scorpii at 0 Sagittarius.
            value = (
                _get_star_position_ecliptic(
                    STARS["MULA"],
                    tjd_tt,
                    0.0,
                    nonut=True,
                    noaberr=eff_noaberr,
                    nogdefl=eff_nogdefl,
                )
                - 240.0
            )
        elif sid_mode == SIDM_TRUE_SHEORAN:
            # Sheoran's rigorous definition (The Science of Time, 2017)
            # anchors the sidereal zero to the live star Asellus Australis
            # (delta Cancri) at the fixed sidereal longitude 103°29'32.9375",
            # explicitly so that "the actual rate of precession [is]
            # irrelevant"; his -60° at the 4174 BCE winter solstice and the
            # 1°/71.75 yr rule are derived checkpoints, not the anchor.
            value = (
                _get_star_position_ecliptic(
                    STARS["PUSHYA"],
                    tjd_tt,
                    0.0,
                    nonut=True,
                    noaberr=eff_noaberr,
                    nogdefl=eff_nogdefl,
                )
                - SHEORAN_TARGET_LON
            )
        elif sid_mode == SIDM_ALDEBARAN_15TAU:
            # Aldebaran held at 15 degrees Taurus (45 degrees absolute),
            # evaluated from the live Hipparcos catalog position. This fixed
            # catalogue anchor is absent from _ABERRANT_ANCHOR_MODES, so its
            # reduction never responds to FLG_TRUEPOS / FLG_NOABERR / FLG_NOGDEFL
            # (measured reference behavior reports a zero shift). The default
            # apparent place keeps its existing aberration+deflection semantics
            # (the toggles below stay at their False defaults), differing from
            # the previous interpolated value only by removing the deflection
            # interpolation artifact (well below a milli-arcsecond off
            # conjunction, ~0.015" at worst near it).
            value = (
                _get_star_position_ecliptic(STARS["ALDEBARAN"], tjd_tt, 0.0, nonut=True)
                - ALDEBARAN_TARGET_LON
            )
        elif sid_mode in GALCENT_TARGET_LON:
            value = _galcent_ayanamsha(
                tjd_tt,
                GALCENT_TARGET_LON[sid_mode],
                noaberr=eff_noaberr,
                nogdefl=eff_nogdefl,
            )
        elif sid_mode in (
            SIDM_GALEQU_IAU1958,
            SIDM_GALEQU_TRUE,
            SIDM_GALEQU_MULA,
        ):
            pole_longitude = _get_star_position_ecliptic(
                STARS["GAL_NORTH_POLE"],
                tjd_tt,
                0.0,
                nonut=True,
                geometric=True,
            )
            ascending_node = (pole_longitude + 90.0) % 360.0
            targets = {
                SIDM_GALEQU_IAU1958: 240.0,
                SIDM_GALEQU_TRUE: 239.94708,
                SIDM_GALEQU_MULA: 246.6137,
            }
            value = ascending_node - targets[sid_mode]
        elif sid_mode == SIDM_GALALIGN_MARDYKS:
            # Mardyks (Sacred Astronomy, 1991): ayanamsha exactly 30 degrees
            # at the September equinox 1998, propagated as a fixed-epoch
            # frame with the IAU 2006 general-precession angle.
            mardyks_value, mardyks_epoch = MARDYKS_DEFINING
            value = mardyks_value + (
                _iau2006_general_precession_deg(tjd_tt)
                - _iau2006_general_precession_deg(mardyks_epoch)
            )
        elif sid_mode == SIDM_GALCENT_MULA_WILHELM:
            # Wilhelm: Sgr A* hour-circle (dhruva) projection held at the
            # middle of Mula (246 degrees 40 minutes).
            value = (
                _mula_wilhelm_longitude(
                    tjd_tt, noaberr=eff_noaberr, nogdefl=eff_nogdefl
                )
                - MULA_WILHELM_TARGET_LON
            )
        elif sid_mode == SIDM_VALENS_MOON:
            from .time_utils import deltat

            defining_epoch_tt = VALENS_MOON_T0_UT + deltat(VALENS_MOON_T0_UT)
            value = VALENS_MOON_AYAN_T0_DEG + method_b_accumulated_precession(
                tjd_tt, defining_epoch_tt
            )
            if bits & SIDBIT_ECL_DATE:
                # Valens is an epoch-pair mode (it lives on this branch only
                # for its UT-anchored epoch), so the ecliptic-of-date
                # projection applies to it like to every defining pair.
                value += _ecl_date_ayanamsha_delta(
                    VALENS_MOON_AYAN_T0_DEG, defining_epoch_tt, tjd_tt
                )
        else:
            raise AssertionError(f"Unhandled dynamic sidereal mode {sid_mode}")
        return float(value % 360.0)

    if sid_mode not in AYANAMSHA_DEFINING:
        # Measured reference behavior: an unrecognized sidereal mode (including
        # negatives and unregistered positive IDs) falls back to the default
        # Fagan/Bradley ayanamsha (mode 0), not to Lahiri.
        warnings.warn(
            f"Unknown sidereal mode {sid_mode} is not recognized; "
            f"falling back to Fagan/Bradley (mode {SIDM_FAGAN_BRADLEY}).",
            UserWarning,
            stacklevel=2,
        )
        sid_mode = SIDM_FAGAN_BRADLEY

    # Registered defining pair: ayanamsha value at the defining epoch,
    # propagated with Method-B on that epoch's mean ecliptic — the same
    # numerical contract as SIDM_USER. Evidence is deliberately per-mode:
    # docs/reference/ayanamsha.md distinguishes primary and geometric anchors
    # from secondary attributions and explicit project conventions.
    defining_value, defining_epoch = AYANAMSHA_DEFINING[sid_mode]
    value = defining_value + method_b_accumulated_precession(tjd_tt, defining_epoch)
    if bits & SIDBIT_ECL_DATE:
        # SIDBIT_ECL_DATE: refer the mode's fixed zero point to the mean
        # ecliptic of date instead of its defining mean ecliptic (Vondrák
        # 2011 pole geometry; see _ecl_date_ayanamsha_delta). Live star and
        # galactic modes are naturally inert: their value is already an
        # of-date longitude, so no delta applies on the dynamic branch.
        value += _ecl_date_ayanamsha_delta(defining_value, defining_epoch, tjd_tt)
    return float(value % 360.0)


def _get_true_ayanamsa(
    tjd_ut: float,
    sid_mode: int | None = None,
    noaberr: bool = False,
    nogdefl: bool = False,
) -> float:
    """
    Get TRUE ayanamsha (mean + nutation) for sidereal planet position calculations.

    Sidereal planet positions require the true ayanamsha (including nutation),
    even though get_ayanamsa_ut() returns the mean ayanamsha. This is standard
    practice for accurate sidereal coordinate conversion.

    Args:
        tjd_ut: Julian Day in Universal Time (UT1)
        sid_mode: Sidereal mode override; reads the global state when None.
            An explicit value lets thread-safe callers (EphemerisContext,
            Horizons) avoid swapping the global sidereal mode.
        noaberr: Forwarded to :func:`_calc_ayanamsa` (FLG_TRUEPOS / FLG_NOABERR
            requests) so the anchor of the aberrant star/galactic modes is
            evaluated without annual aberration.
        nogdefl: Forwarded to :func:`_calc_ayanamsa` (FLG_TRUEPOS / FLG_NOGDEFL
            requests) so the anchor of the aberrant star/galactic modes is
            evaluated without solar gravitational light deflection.

    Returns:
        True ayanamsha in degrees: the mean ayanamsha normalized to
        [0, 360) plus the signed nutation in longitude. The sum is NOT
        re-normalized — near a zero crossing the value may be slightly
        negative or slightly above 360, matching the reference convention
        (e.g. -0.0039 at J2000 for SIDM_J2000 rather than 359.9961).
    """
    if sid_mode is None:
        sid_mode = get_sid_mode()
    assert isinstance(sid_mode, int)

    # Get mean ayanamsha
    mean_ayanamsa = _calc_ayanamsa(tjd_ut, sid_mode, noaberr=noaberr, nogdefl=nogdefl)

    # Add nutation in longitude (IAU 2006/2000A via pyerfa, ~0.01-0.05 mas)
    ts = get_timescale()
    t_obj = ts.ut1_jd(tjd_ut)
    dpsi_rad, _ = erfa.nut06a(2451545.0, t_obj.tt - 2451545.0)
    nutation_deg = math.degrees(dpsi_rad)

    return mean_ayanamsa + nutation_deg


def _get_ayanamsa_for_flags(
    tjd_ut: float, iflag: int, sid_mode: int | None = None
) -> float:
    """Get appropriate ayanamsha based on calculation flags.

    Returns mean ayanamsha (no nutation) when FLG_NONUT or FLG_J2000 is
    set.  J2000 ecliptic coordinates contain no nutation component, so the
    true ayanamsha (mean + Δψ) would introduce a spurious ~9-17″ offset.
    Otherwise returns true ayanamsha (mean + nutation in longitude).

    Args:
        tjd_ut: Julian Day in Universal Time (UT1)
        iflag: Calculation flags bitmask
        sid_mode: Sidereal mode override; reads the global state when None.
            An explicit value lets thread-safe callers (EphemerisContext,
            Horizons) avoid swapping the global sidereal mode.

    Returns:
        Ayanamsha in degrees
    """
    # FLG_TRUEPOS / FLG_NOABERR evaluate the anchor of the aberrant star/galactic
    # ayanamsha modes without annual aberration, and FLG_TRUEPOS / FLG_NOGDEFL
    # without solar light deflection, so the ayanamsha subtracted matches the
    # planet place the same request produces (measured reference behavior).
    # A heliocentric or barycentric request has no Earth to be aberrated by,
    # so the live anchor must be evaluated without annual aberration there
    # too: keeping the geocentric anchor left the subtracted zero point
    # 13.7"-21.3" away from the observer's own frame for the nine
    # _ABERRANT_ANCHOR_MODES (measured 21.09" on Mars, mode 36, FLG_HELCTR).
    # The same reasoning removes the solar deflection, which is likewise an
    # Earth-observer path effect.
    _centered = bool(iflag & (FLG_HELCTR | FLG_BARYCTR))
    noaberr = bool(iflag & (FLG_TRUEPOS | FLG_NOABERR)) or _centered
    nogdefl = bool(iflag & (FLG_TRUEPOS | FLG_NOGDEFL)) or _centered
    if (iflag & FLG_NONUT) or (iflag & FLG_J2000):
        if sid_mode is None:
            sid_mode = get_sid_mode()
        assert isinstance(sid_mode, int)
        return _calc_ayanamsa(tjd_ut, sid_mode, noaberr=noaberr, nogdefl=nogdefl)
    return _get_true_ayanamsa(tjd_ut, sid_mode, noaberr=noaberr, nogdefl=nogdefl)


# SIDBIT projection flags with a distinct realization on the sidereal path:
#   * SIDBIT_ECL_T0 / SIDBIT_SSY_PLANE -> frame projection (see calc_ut).
#   * SIDBIT_USER_UT -> interpret the SIDM_USER t0 as a UT date.
#   * SIDBIT_ECL_DATE -> refer the mode's fixed zero point to the mean
#     ecliptic of date instead of its defining ecliptic of t0 (see
#     _ecl_date_ayanamsha_delta): applied in _calc_ayanamsa for every
#     defining pair (epoch modes, Valens, SIDM_USER) and naturally inert
#     for the live star/galactic modes, whose value is already an of-date
#     longitude. The LEB fast path delegates to _calc_ayanamsa when this
#     bit is set so both backends stay numerically identical.
_APPLIED_SIDBITS = SIDBIT_ECL_T0 | SIDBIT_SSY_PLANE | SIDBIT_USER_UT | SIDBIT_ECL_DATE

# SIDBIT flags accepted for compatibility but WITHOUT a distinct realization in
# this version, so they currently reduce to the base ayanamsha value:
#   * SIDBIT_PREC_ORIG (use the mode's original precession model): the measured
#     effect is below the ~0.1" ayanamsha floor in the standard configuration,
#     so reducing to the default precession reproduces it within floor.
#   * SIDBIT_NO_PREC_OFFSET (drop the mode's precession-offset reconciliation):
#     this is a time-constant per-mode correction (non-zero for only six modes,
#     up to ~0.83") whose value derives from each mode's reference-internal
#     precession-model assignment; it is not reconstructible from published
#     models without fitting the reference output, so it is left as a documented
#     no-op pending an independent, citable derivation.
# They are retained (not warned) here so the public API accepts them silently,
# matching the reference, which emits no warning for any SIDBIT flag.
_ACCEPTED_SIDBITS = SIDBIT_NO_PREC_OFFSET | SIDBIT_PREC_ORIG

# Union retained by set_sid_mode(): bits NOT in this union are still forwarded to
# state.set_sid_mode(), which strips-and-warns genuinely unknown high bits.
_IMPLEMENTED_SIDBITS = _APPLIED_SIDBITS | _ACCEPTED_SIDBITS


def _get_sidereal_bits() -> int:
    """Return the implemented SIDBIT projection flags currently in effect.

    The value lives in ``state._SIDEREAL_BITS``, beside the sidereal mode and
    epoch it qualifies, so that ``reset_session()``/``close()`` clear it with
    them and ``_swapped_context_state()`` swaps it for an EphemerisContext.
    """
    from . import state

    return state._SIDEREAL_BITS


def set_sid_mode(mode: int, t0: float = 0.0, ayan_t0: float = 0.0):
    """
    Set the sidereal zodiac mode for calculations.

    Configures which ayanamsa system to use for sidereal calculations.
    Affects all subsequent position calculations with FLG_SIDEREAL flag.

    Args:
        mode: Sidereal mode constant (SIDM_LAHIRI, SIDM_FAGAN_BRADLEY, etc.),
            optionally combined with a SIDBIT projection flag.
        t0: Reference time (JD) for user-defined ayanamsa (only for SIDM_USER)
        ayan_t0: Ayanamsa value at reference time t0 in degrees (only for SIDM_USER)

    Notes:
        All predefined base IDs 0--46 are accepted and computed. The SIDBIT
        projection flags SIDBIT_ECL_T0, SIDBIT_SSY_PLANE and SIDBIT_USER_UT are
        applied (see calc_ut / _calc_ayanamsa). The remaining three,
        SIDBIT_ECL_DATE, SIDBIT_NO_PREC_OFFSET and SIDBIT_PREC_ORIG, are
        accepted silently for compatibility but currently reduce to the base
        ayanamsha value (see _ACCEPTED_SIDBITS for the per-flag rationale); no
        warning is emitted, matching the reference. Genuinely unknown high bits
        still warn.

    Example:
        >>> set_sid_mode(SIDM_J2000)
        >>> pos, _ = calc_ut(2451545.0, SUN, FLG_SIDEREAL)
        >>> print(f"Sidereal Sun: {pos[0]:.6f}°")

        >>> # Custom ayanamsa: 24° at J2000.0, precessing at standard rate
        >>> set_sid_mode(SIDM_USER, t0=2451545.0, ayan_t0=24.0)
    """
    from . import state
    from .state import set_sid_mode

    if mode >= 0:
        state._SIDEREAL_BITS = mode & _IMPLEMENTED_SIDBITS
        # Forward the base mode plus any still-unsupported projection bits to
        # the state setter, which strips-and-warns for those remaining bits.
        base_mode = mode & ~_IMPLEMENTED_SIDBITS
    else:
        # A negative mode is an invalid ID, not a SIDBIT composite. Extracting
        # projection bits from it (mode & _IMPLEMENTED_SIDBITS on the all-ones
        # two's-complement pattern) would spuriously enable ECL_T0/SSY_PLANE/
        # USER_UT and the `& ~bits` would map it onto a positive base. Pass it
        # through verbatim so the reducer applies the invalid-mode fallback.
        state._SIDEREAL_BITS = 0
        base_mode = mode
    set_sid_mode(base_mode, t0, ayan_t0)


def get_ayanamsa(tjdet: float) -> float:
    """
    Calculate ayanamsa for a given Ephemeris Time (ET/TT) date.

    Similar to get_ayanamsa_ut() but takes Terrestrial Time instead of UT.

    Args:
        tjdet: Julian Day in Ephemeris Time (TT/ET)

    Returns:
        Ayanamsa value in degrees

    Note:
        Properly converts TT to UT1 using Skyfield's timescale with Delta T correction.
        Delta T (TT - UT) varies from ~32s (year 2000) to minutes (historical times).
        While ayanamsa changes slowly (~50"/year), correct conversion ensures
        consistency with the reference API contract.
    """
    ts = get_timescale()
    t_tt = ts.tt_jd(tjdet)
    tjd_ut = t_tt.ut1  # Proper TT to UT1 conversion using Delta T
    return get_ayanamsa_ut(tjd_ut)


# The public "calculated" retflag class drops FLG_NONUT from the extended
# function's echo. Modes outside that class retain it. This classification is
# independent of whether a mode has a native numerical definition.
_AYANAMSA_EX_NONUT_DROP_MODES = frozenset(
    {
        SIDM_GALCENT_0SAG,
        SIDM_TRUE_CITRA,
        SIDM_TRUE_REVATI,
        SIDM_TRUE_PUSHYA,
        SIDM_GALCENT_RGILBRAND,
        SIDM_GALEQU_IAU1958,
        SIDM_GALEQU_TRUE,
        SIDM_GALEQU_MULA,
        SIDM_TRUE_MULA,
        SIDM_GALCENT_MULA_WILHELM,
        SIDM_TRUE_SHEORAN,
        SIDM_GALCENT_COCHRANE,
    }
)


def _ayanamsa_ex_retflag(flags: int, sid_mode: int) -> int:
    """Echoed return flag for get_ayanamsa_ex[_ut].

    Measured against the reference, the echo carries only the resolved
    ephemeris-selection bit plus (conditionally) FLG_NONUT — every other input
    bit (FLG_SPEED, FLG_SIDEREAL, spurious bits) is dropped:

    - The ephemeris bit is a single-winner resolution with priority
      JPLEPH > SWIEPH > MOSEPH; with no ephemeris bit set it defaults to
      SWIEPH. So 0 -> 2, 1 -> 1, 4 -> 4, 258(SWIEPH|SPEED) -> 2, 6(SWIEPH|
      MOSEPH) -> 2, 3(JPLEPH|SWIEPH) -> 1.
    - FLG_NONUT is echoed only when set AND the mode is not one of the
      star/galactic "calculated" modes that drop it (so 64 -> 66, 65 -> 65 in
      the normal modes; 64 -> 2, 65 -> 1 in a drop mode).
    """
    from .constants import FLG_JPLEPH, FLG_MOSEPH

    if flags & FLG_JPLEPH:
        retflag = FLG_JPLEPH
    elif flags & FLG_SWIEPH:
        retflag = FLG_SWIEPH
    elif flags & FLG_MOSEPH:
        retflag = FLG_MOSEPH
    else:
        retflag = FLG_SWIEPH
    if (flags & FLG_NONUT) and sid_mode not in _AYANAMSA_EX_NONUT_DROP_MODES:
        retflag |= FLG_NONUT
    return retflag


def get_ayanamsa_ex(tjdet: float, flags: int = 0) -> Tuple[int, float]:
    """
    Calculate ayanamsa with extended flags for Ephemeris Time.

    Uses the sidereal mode set via set_sid_mode(). Returns the ayanamsa
    value along with the return flags, matching the reference signature.

    Args:
        tjdet: Julian Day in Ephemeris Time (TT/ET)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple of (retflag, ayanamsa):
            - retflag: Return flags (int)
            - ayanamsa: Ayanamsa value in degrees (tropical_lon - sidereal_lon)

    Example:
        >>> from libephemeris import get_ayanamsa_ex, set_sid_mode, SIDM_J2000
        >>> set_sid_mode(SIDM_J2000)
        >>> flags, aya = get_ayanamsa_ex(2451545.0)
        >>> print(f"Ayanamsa: {aya:.6f}")

    Note:
        The sidereal mode must be set via set_sid_mode() before calling.
    """
    sid_mode = get_sid_mode()
    assert isinstance(sid_mode, int)
    ayanamsa = _calc_ayanamsa_ex_value(tjdet, sid_mode, flags)
    # retflag echoes the input flags with FLG_SWIEPH added (reference
    # behaviour: 0 -> 2, FLG_NONUT -> 66), except the star/galactic
    # "calculated" modes drop FLG_NONUT (0 -> 2, FLG_NONUT -> 2).
    return (_ayanamsa_ex_retflag(flags, sid_mode), float(ayanamsa))


def get_ayanamsa_ex_ut(tjdut: float, flags: int = 0) -> Tuple[int, float]:
    """
    Calculate ayanamsa with extended flags for Universal Time.

    This is the UT version of get_ayanamsa_ex(). It internally converts
    from UT to TT before calculating.

    Args:
        tjdut: Julian Day in Universal Time (UT1)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple of (retflag, ayanamsa):
            - retflag: Return flags (int)
            - ayanamsa: Ayanamsa value in degrees (tropical_lon - sidereal_lon)

    Example:
        >>> from libephemeris import get_ayanamsa_ex_ut, set_sid_mode, SIDM_J2000
        >>> set_sid_mode(SIDM_J2000)
        >>> flags, aya = get_ayanamsa_ex_ut(2451545.0)
        >>> print(f"Ayanamsa: {aya:.6f}")

    Note:
        Internally converts UT to TT using Delta T before calculation.
    """
    ts = get_timescale()
    t_ut = ts.ut1_jd(tjdut)
    tjd_tt = t_ut.tt  # Convert UT1 to TT
    sid_mode = get_sid_mode()
    assert isinstance(sid_mode, int)
    ayanamsa = _calc_ayanamsa_ex_value(tjd_tt, sid_mode, flags)
    # retflag echoes the input flags with FLG_SWIEPH added (reference
    # behaviour: 0 -> 2, FLG_NONUT -> 66), except the star/galactic
    # "calculated" modes drop FLG_NONUT (0 -> 2, FLG_NONUT -> 2).
    return (_ayanamsa_ex_retflag(flags, sid_mode), float(ayanamsa))


def _calc_ayanamsa_ex_value(tjd_tt: float, sid_mode: int, flags: int = 0) -> float:
    """
    Internal function to calculate the ayanamsa value for a given TT and sid_mode.

    The _ex variants of the reference API return the *true* ayanamsha
    (mean + nutation in longitude) unless FLG_NONUT is set — that is the
    distinction from the plain get_ayanamsa, which returns the mean value.

    Args:
        tjd_tt: Julian Day in Terrestrial Time (TT)
        sid_mode: Sidereal mode constant
        flags: Calculation flags; FLG_NONUT selects the mean ayanamsha.

    Returns:
        Ayanamsa value in degrees
    """
    ts = get_timescale()
    t_obj = ts.tt_jd(tjd_tt)
    tjd_ut = t_obj.ut1
    # FLG_TRUEPOS / FLG_NOABERR evaluate the aberrant star/galactic anchors
    # without annual aberration, and FLG_TRUEPOS / FLG_NOGDEFL without solar
    # light deflection (measured reference behavior).
    noaberr = bool(flags & (FLG_TRUEPOS | FLG_NOABERR))
    nogdefl = bool(flags & (FLG_TRUEPOS | FLG_NOGDEFL))
    if flags & FLG_NONUT:
        return _calc_ayanamsa(tjd_ut, sid_mode, noaberr=noaberr, nogdefl=nogdefl)
    return _get_true_ayanamsa(tjd_ut, sid_mode, noaberr=noaberr, nogdefl=nogdefl)


class _SkyfieldDeflectorSource:
    """Adapter exposing the Skyfield kernel via LEBReader's eval_body API.

    Lets the Skyfield path reuse fast_calc._apply_gravitational_deflection
    (the shared PPN formula) for the NOABERR-only flag combination, where
    deflection must be applied without aberration.
    """

    _NAMES = {0: "sun", 5: "jupiter barycenter", 6: "saturn barycenter"}

    def __init__(self, planets_kernel):
        self._planets = planets_kernel

    def eval_body(self, body_id: int, jd_tt: float):
        name = self._NAMES[body_id]  # KeyError -> deflector skipped upstream
        ts = get_timescale()
        t = ts.tt_jd(jd_tt)
        at = self._planets[name].at(t)
        p = at.position.au
        v = at.velocity.au_per_d
        return (
            (float(p[0]), float(p[1]), float(p[2])),
            (float(v[0]), float(v[1]), float(v[2])),
        )


def _apply_deflection_only(astrometric, t, icrf_center, planets_kernel):
    """Apply gravitational deflection to an astrometric position (no aberration).

    Skyfield has no built-in deflection-without-aberration mode, so the
    geocentric astrometric vector is deflected with the shared PPN routine
    and re-wrapped as an ICRF position for the downstream frame transforms.
    """
    from skyfield.positionlib import ICRF

    from .fast_calc import _apply_gravitational_deflection

    geo = astrometric.position.au
    geo_t = (float(geo[0]), float(geo[1]), float(geo[2]))
    obs = astrometric.center_barycentric
    obs_pos = obs.position.au
    obs_t = (float(obs_pos[0]), float(obs_pos[1]), float(obs_pos[2]))
    lt = float(astrometric.light_time)

    source = _SkyfieldDeflectorSource(planets_kernel)
    deflected = _apply_gravitational_deflection(
        geo_t,
        obs_t,
        t.tt,
        lt,
        source,
    )
    return ICRF(
        (deflected[0], deflected[1], deflected[2]),
        astrometric.velocity.au_per_d,
        t=t,
        center=icrf_center,
    )


# Position tuple type for nod_aps results
PosTuple = Tuple[float, float, float, float, float, float]


def _apply_nodaps_topocentric(
    result: PosTuple,
    jd_tt: float,
    flags: int,
) -> PosTuple:
    """Subtract the observer state in the exact ``nod_aps`` output frame.

    The generic hypothetical-body reducer starts from a fixed J2000 ecliptic
    frame.  Planetary nodes and apsides instead start in the Vondrák
    ecliptic-of-date frame and can subsequently be rotated to their J2000 or
    sidereal representation.  Transform the ICRS observer offset through that
    same chain before applying parallax.

    The observer velocity includes both the Skyfield/IERS terrestrial-rotation
    state and the time derivative of the selected celestial frame.
    """
    if not (flags & FLG_TOPOCTR):
        return result
    if flags & (FLG_HELCTR | FLG_BARYCTR):
        return result
    if result[2] == 0.0:
        # Preserve sentinel points such as the Sun's undefined node slots and
        # unsupported/Earth results; parallax must not manufacture a body.
        return result

    observer_topo = get_topo()
    if observer_topo is None:
        from .exceptions import ConfigurationError

        raise ConfigurationError(
            "FLG_TOPOCTR requires a geographic position: call "
            "set_topo(lon, lat, alt) first",
            missing_config="observer_location",
            suggestion="Call set_topo(lon, lat, alt) first",
        )

    ts = get_timescale()
    t = ts.tt_jd(jd_tt)
    earth = _get_computation_ephemeris()["earth"]
    topocenter = (earth + observer_topo).at(t)
    geocenter = earth.at(t)
    observer_position: tuple[float, float, float] = (
        float(topocenter.position.au[0] - geocenter.position.au[0]),
        float(topocenter.position.au[1] - geocenter.position.au[1]),
        float(topocenter.position.au[2] - geocenter.position.au[2]),
    )
    observer_velocity: tuple[float, float, float] = (
        float(topocenter.velocity.au_per_d[0] - geocenter.velocity.au_per_d[0]),
        float(topocenter.velocity.au_per_d[1] - geocenter.velocity.au_per_d[1]),
        float(topocenter.velocity.au_per_d[2] - geocenter.velocity.au_per_d[2]),
    )

    def _to_result_frame(
        vector: tuple[float, float, float], epoch_tt: float
    ) -> tuple[float, float, float]:
        if flags & FLG_NONUT:
            matrix = vondrak_precession_matrix(epoch_tt)
            obliquity = vondrak_mean_obliquity_rad(epoch_tt)
        else:
            dpsi, deps = erfa.nut06a(_J2000_JD, epoch_tt - _J2000_JD)
            matrix, obliquity = vondrak_pn_matrix(epoch_tt, float(dpsi), float(deps))
        equatorial = tuple(
            float(sum(matrix[row][column] * vector[column] for column in range(3)))
            for row in range(3)
        )
        ce = math.cos(obliquity)
        se = math.sin(obliquity)
        ecliptic = (
            equatorial[0],
            equatorial[1] * ce + equatorial[2] * se,
            -equatorial[1] * se + equatorial[2] * ce,
        )
        if flags & FLG_J2000:
            ecliptic = _nodaps_vector_to_j2000(ecliptic, epoch_tt)
        if flags & FLG_SIDEREAL:
            epoch = ts.tt_jd(epoch_tt)
            ayanamsha = math.radians(
                _get_ayanamsa_for_flags(float(epoch.ut1), flags & ~FLG_J2000)
            )
            ca = math.cos(ayanamsha)
            sa = math.sin(ayanamsha)
            ecliptic = (
                ecliptic[0] * ca + ecliptic[1] * sa,
                -ecliptic[0] * sa + ecliptic[1] * ca,
                ecliptic[2],
            )
        return (float(ecliptic[0]), float(ecliptic[1]), float(ecliptic[2]))

    observer_position_frame = _to_result_frame(observer_position, jd_tt)
    observer_velocity_physical = _to_result_frame(observer_velocity, jd_tt)
    frame_step = 0.001
    frame_before = _to_result_frame(observer_position, jd_tt - frame_step)
    frame_after = _to_result_frame(observer_position, jd_tt + frame_step)
    observer_velocity_frame = (
        observer_velocity_physical[0]
        + (frame_after[0] - frame_before[0]) / (2.0 * frame_step),
        observer_velocity_physical[1]
        + (frame_after[1] - frame_before[1]) / (2.0 * frame_step),
        observer_velocity_physical[2]
        + (frame_after[2] - frame_before[2]) / (2.0 * frame_step),
    )

    lon, lat, dist, dlon, dlat, ddist = result
    lon_rad = math.radians(lon)
    lat_rad = math.radians(lat)
    dlon_rad = math.radians(dlon)
    dlat_rad = math.radians(dlat)
    cl = math.cos(lon_rad)
    sl = math.sin(lon_rad)
    cb = math.cos(lat_rad)
    sb = math.sin(lat_rad)
    body_position = (
        dist * cb * cl,
        dist * cb * sl,
        dist * sb,
    )
    body_velocity = (
        ddist * cb * cl - dist * sb * dlat_rad * cl - dist * cb * sl * dlon_rad,
        ddist * cb * sl - dist * sb * dlat_rad * sl + dist * cb * cl * dlon_rad,
        ddist * sb + dist * cb * dlat_rad,
    )
    topocentric_position = (
        body_position[0] - observer_position_frame[0],
        body_position[1] - observer_position_frame[1],
        body_position[2] - observer_position_frame[2],
    )
    topocentric_velocity = (
        body_velocity[0] - observer_velocity_frame[0],
        body_velocity[1] - observer_velocity_frame[1],
        body_velocity[2] - observer_velocity_frame[2],
    )

    # The geocentric point already contains annual aberration.  At first
    # post-Newtonian order, changing from the geocentre to a rotating Earth
    # observer adds the aberration from the observer's geocentric velocity.
    # Differentiate that standard direction map locally so its contribution
    # is carried into all three public speed channels.
    if not (flags & (FLG_TRUEPOS | FLG_NOABERR)):
        from .astrometry import apply_aberration_to_position

        def _with_diurnal_aberration(offset: float) -> tuple[float, float, float]:
            position = (
                topocentric_position[0] + offset * topocentric_velocity[0],
                topocentric_position[1] + offset * topocentric_velocity[1],
                topocentric_position[2] + offset * topocentric_velocity[2],
            )
            epoch = ts.tt_jd(jd_tt + offset)
            epoch_topocenter = (earth + observer_topo).at(epoch)
            epoch_geocenter = earth.at(epoch)
            relative_velocity: tuple[float, float, float] = (
                float(
                    epoch_topocenter.velocity.au_per_d[0]
                    - epoch_geocenter.velocity.au_per_d[0]
                ),
                float(
                    epoch_topocenter.velocity.au_per_d[1]
                    - epoch_geocenter.velocity.au_per_d[1]
                ),
                float(
                    epoch_topocenter.velocity.au_per_d[2]
                    - epoch_geocenter.velocity.au_per_d[2]
                ),
            )
            relative_velocity = _to_result_frame(relative_velocity, jd_tt + offset)
            return apply_aberration_to_position(position, relative_velocity)

        topocentric_position = _with_diurnal_aberration(0.0)
        aberration_step = 0.0001
        aberration_before = _with_diurnal_aberration(-aberration_step)
        aberration_after = _with_diurnal_aberration(aberration_step)
        topocentric_velocity = (
            (aberration_after[0] - aberration_before[0]) / (2.0 * aberration_step),
            (aberration_after[1] - aberration_before[1]) / (2.0 * aberration_step),
            (aberration_after[2] - aberration_before[2]) / (2.0 * aberration_step),
        )

    tx, ty, tz = topocentric_position
    tvx, tvy, tvz = topocentric_velocity
    distance_xy_sq = tx * tx + ty * ty
    distance = math.sqrt(distance_xy_sq + tz * tz)
    if distance == 0.0:
        return result
    distance_xy = math.sqrt(distance_xy_sq)
    radial_speed = (tx * tvx + ty * tvy + tz * tvz) / distance
    longitude_speed = (
        math.degrees((tx * tvy - ty * tvx) / distance_xy_sq)
        if distance_xy_sq > 0.0
        else 0.0
    )
    latitude_speed = (
        math.degrees((tvz - tz * radial_speed / distance) / distance_xy)
        if distance_xy > 0.0
        else 0.0
    )
    if not (flags & FLG_SPEED):
        longitude_speed = latitude_speed = radial_speed = 0.0
    return (
        float(math.degrees(math.atan2(ty, tx)) % 360.0),
        float(math.degrees(math.asin(max(-1.0, min(1.0, tz / distance))))),
        float(distance),
        float(longitude_speed),
        float(latitude_speed),
        float(radial_speed),
    )


def _format_nodaps_output(
    points: tuple[PosTuple, PosTuple, PosTuple, PosTuple],
    flags: int,
    jd_tt: float,
) -> tuple[PosTuple, PosTuple, PosTuple, PosTuple]:
    """Apply observer, frame, and representation flags to all four points."""
    transformed: list[PosTuple] = []
    for point in points:
        point = _apply_nodaps_topocentric(point, jd_tt, flags)
        if flags & FLG_EQUATORIAL:
            # Drop FLG_SIDEREAL for the equatorial rotation: the sidereal
            # ayanamsha is already suppressed for EQUATORIAL requests in
            # _calc_nod_aps (SID+EQU == EQU alone), and _maybe_equatorial_convert
            # would otherwise switch to the MEAN obliquity for a sidereal
            # request, leaving a nutation-in-obliquity residual against the pure
            # equatorial reference. FLG_NONUT/FLG_J2000 still select their frames.
            point = _maybe_equatorial_convert(
                point,
                jd_tt,
                flags & ~FLG_SIDEREAL,
                already_j2000=bool(flags & FLG_J2000),
            )
        transformed.append(_to_native_floats(_apply_output_flags(point, flags)))
    formatted = tuple(transformed)
    return (formatted[0], formatted[1], formatted[2], formatted[3])


def _spherical_state_to_xyz(
    lon_deg: float, lat_deg: float, r: float, dlon: float, dlat: float, dr: float
) -> Tuple[float, float, float, float, float, float]:
    """Cartesian position AND velocity from a spherical state (deg, AU, /day).

    The six-slot contract keeps the velocity slots meaningful under FLG_XYZ,
    so the rates are transformed with the position instead of being zeroed.
    """
    lo = math.radians(lon_deg)
    la = math.radians(lat_deg)
    dlo = math.radians(dlon)
    dla = math.radians(dlat)
    cl, sl = math.cos(lo), math.sin(lo)
    cb, sb = math.cos(la), math.sin(la)
    return (
        r * cb * cl,
        r * cb * sl,
        r * sb,
        dr * cb * cl - r * sb * cl * dla - r * cb * sl * dlo,
        dr * cb * sl - r * sb * sl * dla + r * cb * cl * dlo,
        dr * sb + r * cb * dla,
    )


def _nodaps_sidereal_frame_projection(
    entry_fn, tjd: float, planet: int, method: int, flags: int, jd_tt: float
):
    """Frame-projected nod_aps result, or None when no projection applies.

    Sidereal FRAME requests (the fixed-epoch modes and the SIDBIT_ECL_T0 /
    SIDBIT_SSY_PLANE projections) rotate the whole node/apse vector; a scalar
    ayanamsha subtracted from longitude cannot move a latitude at all, which
    left every projected latitude at its unprojected value.

    Measured reference behavior: the ecliptic result is the OF-DATE node/apse
    set rotated rigidly into the target frame. That is what keeps the node
    pair exactly antipodal — the of-date nodes lie on the ecliptic of date
    (latitude ~1e-5 deg) and a rigid rotation preserves antipodality.
    Re-requesting the points in the J2000 frame instead returned that frame's
    own geocentric node solution, whose latitudes are neither antipodal nor a
    rotation of the of-date ones (measured 55.9" / 38.1" off on Jupiter at
    JD 2415020.5).

    Equatorial output follows the calc surfaces: SIDBIT_SSY_PLANE alone
    leaves the equator untouched, and SIDBIT_ECL_T0 reduces the J2000|NONUT
    request to the mean equator of t0 (applying the *ecliptic* matrix to an
    equatorial vector had put Jupiter's node 23.99 deg off in RA).
    """
    if not (flags & FLG_SIDEREAL):
        return None
    import numpy as np

    from .sidereal_epoch import (
        _ecliptic_of_t0_matrix,
        _epoch_matrices,
        _rotate_spherical,
        equatorial_epoch_matrix,
        is_fixed_epoch_request,
        sidbit_ecliptic_matrix,
        ssy_plane_zero_point_deg,
    )
    from .sidereal_epoch import FIXED_EPOCH_T0

    sid_mode = get_sid_mode()
    fixed = is_fixed_epoch_request(flags, sid_mode)
    bits = _get_sidereal_bits()
    sidbit = bool(bits & (SIDBIT_ECL_T0 | SIDBIT_SSY_PLANE)) and (
        sid_mode not in _SIDBIT_PROJECTION_SUPPRESS_MODES
    )
    if not fixed and not sidbit:
        return None

    if flags & FLG_EQUATORIAL:
        if sidbit and not fixed and not (bits & SIDBIT_ECL_T0):
            return None
        # Same rigid-rotation construction as the ecliptic branch, on the
        # equatorial side: the of-date mean-equator result rotated into the
        # mean equator of t0.
        t0_eq = FIXED_EPOCH_T0[sid_mode] if fixed else _ecl_t0_epoch_jd(sid_mode)
        # Same two constructions as the ecliptic branch below, selected by the
        # caller's FLG_J2000. Leaving the J2000 case out of this branch made
        # the flag inert here: the equatorial projection returned bit-identical
        # output with and without it (measured 72.3" of declination and 24.9"
        # of RA on Venus's node at JD 2488070).
        eq_base = (flags & ~(FLG_SIDEREAL | FLG_RADIANS | FLG_XYZ)) | FLG_NONUT
        if flags & FLG_J2000:
            sub = entry_fn(tjd, planet, method, eq_base)
            m_eq = np.asarray(equatorial_epoch_matrix(t0_eq))
        else:
            sub = entry_fn(tjd, planet, method, eq_base & ~FLG_J2000)
            m_eq = (
                equatorial_epoch_matrix(t0_eq)
                @ np.asarray(equatorial_epoch_matrix(jd_tt)).T
            )

        def _project_eq(pt):
            out = _rotate_spherical(tuple(pt), m_eq)
            if flags & FLG_XYZ:
                return _to_native_floats(_spherical_state_to_xyz(*out))
            if flags & FLG_RADIANS:
                out = (
                    math.radians(out[0]),
                    math.radians(out[1]),
                    out[2],
                    math.radians(out[3]),
                    math.radians(out[4]),
                    out[5],
                )
            return _to_native_floats(out)

        return tuple(_project_eq(pt) for pt in sub)

    if fixed:
        _, m_target = _epoch_matrices(sid_mode)
    else:
        if bits & SIDBIT_ECL_T0:
            t0_jd = _ecl_t0_epoch_jd(sid_mode)
            zero_point = _ecl_t0_zero_point_deg(sid_mode, t0_jd)
        else:
            t0_jd = _J2000_JD
            zero_point = ssy_plane_zero_point_deg(_calc_ayanamsa(_J2000_JD, sid_mode))
        m_sidbit = sidbit_ecliptic_matrix(bits, t0_jd, zero_point)
        if m_sidbit is None:
            return None
        m_target = m_sidbit

    # The sub-request is always SPHERICAL and of-date: the rigid rotation is
    # defined on that frame. FLG_XYZ and FLG_RADIANS are representations of
    # the ROTATED result, applied below — leaving FLG_XYZ in the sub-request
    # returned Cartesian components that _rotate_spherical then read as
    # (lon, lat, r). FLG_J2000 is likewise a property of the output frame:
    # the caller's J2000 reduction is composed into the rotation instead of
    # being taken from the sub-request, whose own J2000 solution is not a
    # rotation of the of-date one.
    # Two constructions, selected by the caller's FLG_J2000:
    #  * without it, the sub-request is the MEAN-OF-DATE set and the rotation
    #    carries it into the target frame — this is what keeps the node pair
    #    antipodal;
    #  * with it, the caller has asked for the J2000 reduction first (the
    #    reference honours its distance-dependent parallax term there), so
    #    the sub-request is taken in the J2000 frame and rotated J2000 -> t0.
    #    Reusing the of-date construction here dropped that reduction and
    #    left the latitude 74" off at 1900, 274" at 1600.
    base_flags = (flags & ~(FLG_SIDEREAL | FLG_RADIANS | FLG_XYZ)) | FLG_NONUT
    if flags & FLG_J2000:
        sub = entry_fn(tjd, planet, method, base_flags)
        m = m_target
    else:
        sub = entry_fn(tjd, planet, method, base_flags & ~FLG_J2000)
        m = m_target @ np.asarray(_ecliptic_of_t0_matrix(jd_tt)).T

    def _project(pt):
        out = _rotate_spherical(tuple(pt), m)
        if flags & FLG_XYZ:
            return _to_native_floats(_spherical_state_to_xyz(*out))
        if flags & FLG_RADIANS:
            out = (
                math.radians(out[0]),
                math.radians(out[1]),
                out[2],
                math.radians(out[3]),
                math.radians(out[4]),
                out[5],
            )
        return _to_native_floats(out)

    return tuple(_project(pt) for pt in sub)


def nod_aps_ut(
    tjdut: float,
    planet: int,
    method: int = NODBIT_MEAN,
    flags: int = FLG_SWIEPH | FLG_SPEED,
) -> Tuple[PosTuple, PosTuple, PosTuple, PosTuple]:
    """
    Calculate planetary nodes and apsides for Universal Time.

    Reference API compatible function.

    This function computes the orbital nodes (ascending/descending) and apsides
    (perihelion/aphelion) for any planet. The nodes are the points where the
    planet's orbital plane intersects the ecliptic plane. The apsides are the
    points of closest (perihelion) and farthest (aphelion) approach to the Sun.

    Args:
        tjdut: Julian Day in Universal Time (UT1)
        planet: Planet/body ID (SUN, MOON, etc.)
        method: Method for node/apse calculation:
            - NODBIT_MEAN (1): Mean orbital elements (averaged)
            - NODBIT_OSCU (2): Osculating elements (instantaneous)
            - NODBIT_OSCU_BAR (4): Barycentric osculating elements
            - NODBIT_FOPOINT (256): Include focal point
        flags: Calculation flags (FLG_SPEED, etc.)

    Returns:
        Tuple of 4 position tuples, each containing 6 floats:
            - xnasc: Ascending node (lon, lat, dist, speed_lon, speed_lat, speed_dist)
            - xndsc: Descending node (same format)
            - xperi: Perihelion (same format)
            - xaphe: Aphelion (same format)

    Example:
        >>> from libephemeris import nod_aps_ut, MARS, NODBIT_MEAN
        >>> nasc, ndsc, peri, aphe = nod_aps_ut(2451545.0, MARS, NODBIT_MEAN)
        >>> print(f"Mars ascending node: {nasc[0]:.4f}°")
        >>> print(f"Mars perihelion: {peri[0]:.4f}°")

    Note:
        This function uses mean orbital elements for reliable results.
        For planets, the mean elements provide smooth, predictable values.
        Osculating elements can show rapid variations due to perturbations.
    """
    ts = get_timescale()
    t = ts.ut1_jd(tjdut)
    projected = _nodaps_sidereal_frame_projection(
        nod_aps_ut, tjdut, planet, method, flags, float(t.tt)
    )
    if projected is not None:
        return projected
    return _format_nodaps_output(
        _calc_nod_aps(t, planet, flags, method), flags, float(t.tt)
    )


def nod_aps(
    tjdet: float,
    planet: int,
    method: int = NODBIT_MEAN,
    flags: int = FLG_SWIEPH | FLG_SPEED,
) -> Tuple[PosTuple, PosTuple, PosTuple, PosTuple]:
    """
    Calculate planetary nodes and apsides for Ephemeris Time (ET/TT).

    Reference API compatible function. Similar to nod_aps_ut() but takes
    Terrestrial Time (TT, also known as Ephemeris Time) instead of Universal Time.

    Args:
        tjdet: Julian Day in Terrestrial Time (TT/ET)
        planet: Planet/body ID (SUN, MOON, etc.)
        method: Method for node/apse calculation (NODBIT_MEAN, etc.)
        flags: Calculation flags (default: FLG_SWIEPH | FLG_SPEED)

    Returns:
        Same as nod_aps_ut: (xnasc, xndsc, xperi, xaphe)

    Example:
        >>> from libephemeris import nod_aps, JUPITER, NODBIT_OSCU
        >>> nasc, ndsc, peri, aphe = nod_aps(2451545.0, JUPITER, NODBIT_OSCU)
    """
    ts = get_timescale()
    t = ts.tt_jd(tjdet)
    projected = _nodaps_sidereal_frame_projection(
        nod_aps, tjdet, planet, method, flags, float(t.tt)
    )
    if projected is not None:
        return projected
    return _format_nodaps_output(
        _calc_nod_aps(t, planet, flags, method), flags, float(t.tt)
    )


class HeliocentricNodApsWarning(UserWarning):
    """Warning for heliocentric nod_aps methodology differences.

    .. deprecated::
        This warning is no longer emitted since libephemeris now uses geocentric
        osculating elements for nod_aps calculations, matching the reference API's
        approach. Kept for backward compatibility.
    """

    pass


_J2000_JD = 2451545.0

# GM_sun in AU^3/day^2 as k^2, where k = 0.01720209895 is the Gaussian
# gravitational constant of the IAU (1976) System of Astronomical Constants
# (Explanatory Supplement to the Astronomical Almanac, 3rd ed., Tables 8.1/15.1).
_GM_SUN = 0.01720209895**2

# Reciprocal system-mass ratios (Sun / planet system), from the DE440
# planetary masses (Park, Folkner, Williams & Boggs 2021, AJ 161, 105;
# cf. Astronomical Almanac Table K7). Osculating elements use the two-body
# GM = GM_sun * (1 + 1/ratio): the perihelion direction of a low-eccentricity
# orbit is sensitive to GM at the ~0.1 deg level for Jupiter, so the planet
# mass cannot be neglected.
_SUN_MASS_RATIOS = {
    MERCURY: 6023600.0,
    VENUS: 408523.71,
    EARTH: 328900.5596,  # Earth + Moon
    MARS: 3098708.0,
    JUPITER: 1047.3486,
    SATURN: 3497.898,
    URANUS: 22902.98,
    NEPTUNE: 19412.24,
    PLUTO: 135200000.0,
}

# Sun-to-Earth mass ratio (Earth alone, not the Earth+Moon system), IAU 2009
# nominal value (Luzum et al. 2011, "The IAU 2009 system of astronomical
# constants", Celestial Mechanics and Dynamical Astronomy 110, 293). Used only
# for the two geocentric symbolic points (White Moon, Waldemath), whose
# osculating ellipse is taken about the Earth rather than the Sun.
_SUN_EARTH_MASS_RATIO = 332946.0487

# --- Astronomical-Almanac-style osculating central mass (FLG_ORBEL_AA) ---
#
# Under FLG_ORBEL_AA the two-body osculating fit uses the hierarchical
# JACOBI-coordinate central mass of classical celestial mechanics: ordering
# the planetary systems outward from the Sun, the Kepler part of the Jacobi
# Hamiltonian for the k-th body carries the gravitational parameter
#     GM_eff(k) = G * (M_sun + m_1 + ... + m_k),
# i.e. the Sun plus every planetary system at or interior to the body's own
# orbit, the body's own system included — the hierarchical Jacobi
# decomposition of the N-body problem, a standard construction of
# celestial mechanics (see e.g. H.C. Plummer, "An Introductory Treatise
# on Dynamical Astronomy", 1918; A. Morbidelli, "Modern Celestial
# Mechanics", 2002; cf. the osculating-element conventions of the
# Explanatory Supplement to the Astronomical Almanac, 3rd ed.). For a minor body the same
# hierarchy gives the restricted-problem parameter Sun + interior planets
# with no self term, its own mass being negligible; minor bodies are not
# themselves perturbers in the hierarchy, whose members are the nine major
# planetary systems. The "mass of a planet" is the mass of its SYSTEM
# (planet plus satellites), the published convention of the DE ephemerides -
# Earth therefore enters as the Earth+Moon barycentre.
#
# Only the osculating GEOMETRY (a, e, q, Q and the perihelion-referred
# angles om/varpi/M/nu/E) responds to the changed parameter: a shrinks and
# the angles shift solidly about the perihelion while the mean longitude L
# stays invariant; the derived period family keeps the pure solar GM, so its
# slots move only through the changed a. The Moon's geocentric orbit is
# unaffected. Reciprocal system-mass ratios come from _SUN_MASS_RATIOS
# (DE440; Park et al. 2021, cf. Astronomical Almanac Table K7). Nominal
# semi-major axes only order the nine systems and decide which planets are
# interior to a minor body's orbit; the ordering itself is fixed.
_ORBEL_AA_ORDER: tuple[int, ...] = (
    MERCURY,
    VENUS,
    EARTH,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
)
_ORBEL_AA_NOMINAL_A_AU: dict[int, float] = {
    MERCURY: 0.387,
    VENUS: 0.723,
    EARTH: 1.000,
    MARS: 1.524,
    JUPITER: 5.204,
    SATURN: 9.582,
    URANUS: 19.19,
    NEPTUNE: 30.07,
    PLUTO: 39.48,
}


def _orbel_aa_cumulative_through_self() -> dict[int, float]:
    """Map each major planet id to the summed reciprocal system-mass ratios of
    the Sun's companions up to and including it (Sun + interior planets + self),
    the FLG_ORBEL_AA central-mass excess for that planet."""
    out: dict[int, float] = {}
    total = 0.0
    for pid in _ORBEL_AA_ORDER:
        total += 1.0 / _SUN_MASS_RATIOS[pid]
        out[pid] = total
    return out


_ORBEL_AA_CUM_THROUGH_SELF: dict[int, float] = _orbel_aa_cumulative_through_self()


def _orbel_aa_interior_mass_sum(a_au: float) -> float:
    """Summed reciprocal system-mass ratios of the major planets strictly
    interior to semi-major axis ``a_au`` -- the FLG_ORBEL_AA central-mass excess
    for a minor body: the restricted-problem Jacobi parameter (Sun + every
    major planet inside its orbit; the body's own negligible mass adds no
    term, and minor bodies are not members of the hierarchy)."""
    total = 0.0
    for pid in _ORBEL_AA_ORDER:
        if _ORBEL_AA_NOMINAL_A_AU[pid] < a_au:
            total += 1.0 / _SUN_MASS_RATIOS[pid]
    return total


# Fictitious/hypothetical bodies whose published construction is *geocentric*
# (they orbit the Earth, not the Sun): the White Moon Selena (56) and the
# Waldemath dark moon (58). Their osculating elements are reduced with GM_earth.
_FICT_GEOCENTRIC_IDS = frozenset({WHITE_MOON, WALDEMATH})

# Heliocentric fictitious body whose runtime model is exactly *circular* and is
# referred to the mean equinox of date: Proserpina (57). Its osculating elements
# are reduced in that native of-date frame rather than J2000, because the tiny
# Gaussian mean motion (~0.0014 deg/day) would otherwise be contaminated by the
# ~50.29"/yr equinox precession rate that a J2000 request injects into the
# tangential velocity (see _calc_orbital_elements_fictitious). The genuinely
# eccentric of-date model Vulcan (55) is not listed: its ~19 deg/day motion
# makes the precession-rate term negligible, so its J2000 reduction is unaffected.
_FICT_ECL_OF_DATE_CIRCULAR_IDS = frozenset({PROSERPINA})

# Astronomical points that carry no heliocentric orbit and for which the
# reference raises "object N not valid" from get_orbital_elements: the lunar
# nodes (10, 11), the lunar apogees (12, 13), and the interpolated lunar
# apsides (21, 22). The Sun (0) is added by get_orbital_elements but is a
# valid target for orbit_max_min_true_distance (its geocentric range is the
# Earth's orbit), so it is handled separately at the two entry points.
_ORBITAL_INVALID_OBJECT_IDS = frozenset(
    {MEAN_NODE, TRUE_NODE, MEAN_APOG, OSCU_APOG, INTP_APOG, INTP_PERG}
)


def _validate_orbital_object(caller: str, ipl: int, *, sun_valid: bool) -> None:
    """Reject ids that have no orbital elements, matching the reference.

    Measured reference behavior (behavioral comparison only, no source
    inspected): ``get_orbital_elements`` raises ``object N not valid`` for the
    Sun, the lunar nodes/apogees (10-13) and the interpolated apsides (21, 22),
    ``illegal planet number N`` for the undefined block 23-39, and
    ``object N not valid`` for negatives. ``orbit_max_min_true_distance`` shares
    the inner ``get_orbital_elements`` wording but treats the Sun as valid.

    Args:
        caller: Public function name to prefix the message with.
        ipl: Requested body id.
        sun_valid: When False the Sun (id 0) is rejected as ``not valid``.

    Raises:
        Error: With the reference-format message for an id without elements.
    """
    from .exceptions import Error

    if ipl < 0 or ipl in _ORBITAL_INVALID_OBJECT_IDS or (ipl == SUN and not sun_valid):
        raise Error(
            f"{caller}: error in get_orbital_elements(): object {ipl} not valid"
        )
    if 23 <= ipl <= 39:
        raise Error(f"{caller}: illegal planet number {ipl}.")
    if WALDEMATH < ipl < AST_OFFSET:
        # No element source exists between the last fictitious body (58) and
        # the numbered-asteroid offset: the measured reference raises for the
        # whole block (its fictitious-element lookup fails), and the position
        # pipeline rejects these ids too. Without this guard the element
        # calculator leaked a zero-initialized 50-tuple that read as a real
        # orbit at 0 AU.
        raise Error(
            f"{caller}: error in get_orbital_elements(): object {ipl} not valid"
        )


def _nodaps_to_j2000(lon_deg: float, lat_deg: float, jd_tt: float) -> tuple:
    """Precess a nod_aps of-date ecliptic coordinate to the J2000 frame.

    Rotate ecliptic→equatorial and back with the J2000 mean obliquity on
    both sides around the equatorial Vondrák precession. Nutation already
    present in the of-date coordinate is deliberately retained.
    """
    cl = math.cos(math.radians(lat_deg))
    v = _nodaps_vector_to_j2000(
        (
            cl * math.cos(math.radians(lon_deg)),
            cl * math.sin(math.radians(lon_deg)),
            math.sin(math.radians(lat_deg)),
        ),
        jd_tt,
    )
    lon_j = math.degrees(math.atan2(v[1], v[0])) % 360.0
    r = math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    lat_j = math.degrees(math.asin(max(-1.0, min(1.0, v[2] / r))))
    return lon_j, lat_j


def _nodaps_vector_to_j2000(
    vector: tuple[float, float, float], jd_tt: float
) -> tuple[float, float, float]:
    """Rotate an ecliptic-of-date vector into the nod_aps J2000 frame.

    The transformation is linear, so it can be applied identically to a
    position vector and to its velocity vector.  Keeping that property
    explicit is important for the public ``nod_aps(..., FLG_SPEED)`` state
    convention.
    """
    eps_j = vondrak_mean_obliquity_rad(_J2000_JD)
    p_date = vondrak_precession_matrix(jd_tt)
    p_j2k = vondrak_precession_matrix(_J2000_JD)
    v = [float(vector[0]), float(vector[1]), float(vector[2])]
    ce, se = math.cos(eps_j), math.sin(eps_j)
    # ecliptic → equatorial (J2000 obliquity)
    v = [v[0], v[1] * ce - v[2] * se, v[1] * se + v[2] * ce]
    # equator of date → ICRS (transpose) → equator J2000
    v = [
        p_date[0][0] * v[0] + p_date[1][0] * v[1] + p_date[2][0] * v[2],
        p_date[0][1] * v[0] + p_date[1][1] * v[1] + p_date[2][1] * v[2],
        p_date[0][2] * v[0] + p_date[1][2] * v[1] + p_date[2][2] * v[2],
    ]
    v = [
        p_j2k[0][0] * v[0] + p_j2k[0][1] * v[1] + p_j2k[0][2] * v[2],
        p_j2k[1][0] * v[0] + p_j2k[1][1] * v[1] + p_j2k[1][2] * v[2],
        p_j2k[2][0] * v[0] + p_j2k[2][1] * v[1] + p_j2k[2][2] * v[2],
    ]
    # equatorial → ecliptic (J2000 obliquity)
    v = [v[0], v[1] * ce + v[2] * se, -v[1] * se + v[2] * ce]
    return (float(v[0]), float(v[1]), float(v[2]))


# Minor bodies whose nodes/apsides the reference computes but this version
# does not yet: the curated asteroids/centaurs plus any numbered asteroid
# (AST_OFFSET..FIXSTAR_OFFSET). nod_aps raises for these rather than
# returning a plausible-looking zero. Planned for a future release; see NEXT.md.
_MINOR_BODY_NODAPS = frozenset({CHIRON, PHOLUS, CERES, PALLAS, JUNO, VESTA})

# Bodies whose NODBIT_MEAN branch is backed by the independent mean-element
# model.  Pluto has no entry in that model and intentionally remains on the
# osculating-state speed path even when callers pass NODBIT_MEAN.
# IAU 2015 nominal solar radius divided by the exact IAU astronomical unit.
_NOMINAL_SOLAR_RADIUS_AU = 6.957e8 / 149_597_870_700.0


def _nodaps_ray_crosses_solar_disk(
    observer_to_target: tuple[float, float, float],
    observer_barycentric: tuple[float, float, float],
    sun_barycentric: tuple[float, float, float],
) -> bool:
    """Return whether the finite target ray intersects the nominal solar disk.

    The point-mass PPN reduction is valid outside the gravitating body's
    photosphere. Without an independently sourced solar mass profile, applying
    its singularity to an interior ray would be physically invalid.
    """
    distance = math.sqrt(sum(component * component for component in observer_to_target))
    if distance == 0.0:
        return False
    direction = tuple(component / distance for component in observer_to_target)
    observer_to_sun = tuple(
        sun_barycentric[index] - observer_barycentric[index] for index in range(3)
    )
    along_ray = sum(observer_to_sun[index] * direction[index] for index in range(3))
    if along_ray <= 0.0 or along_ray >= distance:
        return False
    closest = tuple(
        observer_to_sun[index] - along_ray * direction[index] for index in range(3)
    )
    impact = math.sqrt(sum(component * component for component in closest))
    return impact < _NOMINAL_SOLAR_RADIUS_AU


def _calc_nod_aps(
    t, ipl: int, iflag: int, method: int
) -> Tuple[PosTuple, PosTuple, PosTuple, PosTuple]:
    """
    Calculate orbital nodes and apsides using geocentric osculating elements.

    Computes the heliocentric osculating orbital elements from JPL DE440 state
    vectors, determines the 3D heliocentric positions of the ascending node,
    descending node, perihelion, and aphelion on the instantaneous orbit, then
    converts these positions to geocentric ecliptic coordinates of date.

    The returned longitudes are the node/apse directions in the ecliptic as
    seen from the selected observer.

    Args:
        t: Skyfield Time object
        ipl: Planet ID (MERCURY, VENUS, MARS, etc.)
        iflag: Calculation flags
        method: Node/apse calculation method (NODBIT_MEAN, NODBIT_OSCU, etc.)
            Currently all methods use osculating elements from JPL ephemeris.

    Returns:
        Tuple of (ascending_node, descending_node, perihelion, aphelion)
        Each element is a PosTuple: (longitude, latitude, distance, dlon, dlat, ddist)
    """
    zero_pos: PosTuple = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    # Minor bodies (curated asteroids/centaurs 15-20 and any numbered asteroid
    # AST_OFFSET+n) carry a real heliocentric orbit, so their nodes/apsides are
    # computed from the osculating state served by the same position pipeline as
    # calc/get_orbital_elements (see the osculating branch below). Same
    # membership test as _calc_orbital_elements so the two stay in lockstep.
    is_minor = ipl in _MINOR_BODY_NODAPS or (AST_OFFSET < ipl < FIXSTAR_OFFSET)

    # Method-bit precedence (measured reference behavior): NODBIT_MEAN wins
    # whenever it is set, even alongside NODBIT_OSCU / NODBIT_OSCU_BAR, so
    # methods 3/5/7 track method 1 (mean), not 2/4/6 (osculating); method 0
    # (no bits) also defaults to mean. An osculating bit only takes effect when
    # the mean bit is absent. The barycentric center (NODBIT_OSCU_BAR) is an
    # independent choice that applies only on the osculating path (see bar_mode
    # below); a body without mean elements (Pluto, minor bodies) therefore still
    # honors OSCU_BAR under method 5/7 because it falls through to osculating.
    prefer_mean = bool(method & NODBIT_MEAN) or not (
        method & (NODBIT_OSCU | NODBIT_OSCU_BAR)
    )

    # Interpolated lunar apsides (INTP_APOG / INTP_PERG): the interpolated
    # point is not an orbital-element body, so it has no node/apse
    # decomposition. Measured reference behavior returns not-a-number in the
    # three position slots (longitude, latitude, distance) for every method;
    # the speed slots stay 0.0 on a plain request, while an explicit
    # FLG_SPEED request propagates NaN through the speed channel too.
    # Mirror both cases exactly so callers can detect the undefined place.
    if ipl in (INTP_APOG, INTP_PERG):
        # Measured reference behavior: the undefined point keeps ZERO speed
        # slots on a plain request, but an explicit FLG_SPEED request
        # propagates NaN through the speed channel too.
        _spd = math.nan if iflag & FLG_SPEED else 0.0
        nan_pos: PosTuple = (
            math.nan,
            math.nan,
            math.nan,
            _spd,
            _spd,
            _spd,
        )
        return (nan_pos, nan_pos, nan_pos, nan_pos)

    # Bodies for which nodes/apsides are undefined here (mean/true node and
    # apogee, Uranians/fictitious bodies, and any out-of-range id): measured
    # reference behavior raises an error rather than returning a silent zero
    # that would read as a node at 0° Aries. Earth stays in _PLANET_MAP and is
    # handled below with its nodes zeroed (its orbit defines the ecliptic).
    # Minor bodies are not in _PLANET_MAP but carry a real orbit (handled below).
    if ipl not in _PLANET_MAP and not is_minor:
        from .exceptions import Error

        raise Error(f"nod_aps: nodes/apsides for body id {ipl} are not implemented.")

    # Speeds are near-instantaneous central differences of the independently
    # calculated positions.  Every speed slot is therefore a physical rate of
    # the corresponding position or distance channel.
    if iflag & FLG_SPEED:
        ts = get_timescale()
        base_flags = iflag & ~FLG_SPEED
        now = _calc_nod_aps(t, ipl, base_flags, method)
        dt_spd = 0.001  # days (~86 s)
        prev = _calc_nod_aps(ts.tt_jd(t.tt - dt_spd), ipl, base_flags, method)
        nxt = _calc_nod_aps(ts.tt_jd(t.tt + dt_spd), ipl, base_flags, method)

        def _with_speed_c(cur: PosTuple, pre: PosTuple, nex: PosTuple) -> PosTuple:
            dlon = (nex[0] - pre[0] + 180.0) % 360.0 - 180.0
            return (
                cur[0],
                cur[1],
                cur[2],
                dlon / (2.0 * dt_spd),
                (nex[1] - pre[1]) / (2.0 * dt_spd),
                (nex[2] - pre[2]) / (2.0 * dt_spd),
            )

        return (
            _with_speed_c(now[0], prev[0], nxt[0]),
            _with_speed_c(now[1], prev[1], nxt[1]),
            _with_speed_c(now[2], prev[2], nxt[2]),
            _with_speed_c(now[3], prev[3], nxt[3]),
        )

    from .state import _get_computation_ephemeris

    planets = _get_computation_ephemeris()
    jd_tt = t.tt

    # --- ICRS → ecliptic of date rotation (Vondrák 2011, matching calc()) ---
    # Step 1: precession(-nutation) matrix (ICRS → equator of date) from the same
    # Vondrák long-term precession used by the main reduction. Under FLG_NONUT
    # the mean equator/obliquity are used (no nutation); otherwise nutation stays
    # IAU 2006/2000A. Step 2: the matching obliquity rotates equator-of-date →
    # ecliptic-of-date.
    # FLG_J2000 is handled at the end by precessing the finished of-date
    # coordinates to J2000: the reference does NOT strip nutation from
    # nod_aps J2000 output (unlike calc), so the of-date frame stays as-is.
    if iflag & FLG_NONUT:
        pnm = vondrak_precession_matrix(jd_tt)
        eps_rad = vondrak_mean_obliquity_rad(jd_tt)
    else:
        dpsi_n, deps_n = erfa.nut06a(_J2000_JD, jd_tt - _J2000_JD)
        pnm, eps_rad = vondrak_pn_matrix(jd_tt, float(dpsi_n), float(deps_n))
    cos_eps = math.cos(eps_rad)
    sin_eps = math.sin(eps_rad)

    def _icrs_to_ecliptic(vec):
        """Rotate an ICRS vector to ecliptic of date."""
        # Apply precession-nutation
        x_eq = pnm[0][0] * vec[0] + pnm[0][1] * vec[1] + pnm[0][2] * vec[2]
        y_eq = pnm[1][0] * vec[0] + pnm[1][1] * vec[1] + pnm[1][2] * vec[2]
        z_eq = pnm[2][0] * vec[0] + pnm[2][1] * vec[1] + pnm[2][2] * vec[2]
        # Rotate equatorial → ecliptic
        return (x_eq, y_eq * cos_eps + z_eq * sin_eps, -y_eq * sin_eps + z_eq * cos_eps)

    def _ecliptic_to_icrs(vec):
        """Rotate an ecliptic-of-date vector back to ICRS."""
        x_eq = vec[0]
        y_eq = vec[1] * cos_eps - vec[2] * sin_eps
        z_eq = vec[1] * sin_eps + vec[2] * cos_eps
        return (
            pnm[0][0] * x_eq + pnm[1][0] * y_eq + pnm[2][0] * z_eq,
            pnm[0][1] * x_eq + pnm[1][1] * y_eq + pnm[2][1] * z_eq,
            pnm[0][2] * x_eq + pnm[1][2] * y_eq + pnm[2][2] * z_eq,
        )

    # --- Get Sun and Earth positions in ICRS ---
    sun = planets["sun"]
    earth = planets["earth"]
    sun_pos = sun.at(t)
    earth_pos = earth.at(t)

    # Earth's heliocentric position in ecliptic of date
    r_earth_icrs = earth_pos.position.au - sun_pos.position.au
    r_earth_ecl = _icrs_to_ecliptic(r_earth_icrs)

    # --- Get target's state vector ---
    # Minor bodies are not Skyfield ephemeris segments; their heliocentric state
    # comes from the position pipeline in the osculating branch below, so skip
    # the _PLANET_MAP lookup here (target_pos is unused for them).
    if not is_minor:
        target_name = _PLANET_MAP[ipl]
        try:
            target = planets[target_name]
        except KeyError:
            if target_name in _PLANET_FALLBACK:
                target = planets[_PLANET_FALLBACK[target_name]]
            else:
                return (zero_pos, zero_pos, zero_pos, zero_pos)

        target_pos = target.at(t)

    # For Moon, use lunar theory bodies rather than computing osculating
    # elements from geocentric state vectors.
    # NODBIT_MEAN → MEAN_NODE + MEAN_APOG
    # NODBIT_OSCU → TRUE_NODE + OSCU_APOG
    if ipl == MOON:
        jd_ut = t.ut1
        # Strip output-format bits (FLG_RADIANS / FLG_XYZ) and FLG_SPEED:
        # we need calc_ut() to return spherical (lon_deg, lat_deg, dist)
        # below so we can re-assemble nod_aps output in those units.
        # FLG_J2000 is stripped too: the reference precesses the finished
        # true-of-date points to J2000 (nutation NOT removed), which is
        # applied at the end of this branch.
        # FLG_EQUATORIAL and FLG_SIDEREAL are stripped so calc_ut() returns
        # ecliptic-of-date longitudes: the equatorial rotation is applied
        # exactly once by _format_nodaps_output (leaving them in would
        # double-rotate the point by the obliquity), and the sidereal
        # subtraction is applied by _moon_point with the mean ayanamsha to
        # match the planet path (see below). FLG_NONUT is retained so the
        # of-date ecliptic longitude is mean-of-date when requested.
        calc_flags = (
            iflag
            & ~FLG_SPEED
            & ~FLG_RADIANS
            & ~FLG_XYZ
            & ~FLG_J2000
            & ~FLG_EQUATORIAL
            & ~FLG_SIDEREAL
            & ~FLG_HELCTR
            & ~FLG_BARYCTR
            & ~FLG_TOPOCTR
        )

        # Select node and apogee bodies based on method.
        # NODBIT_OSCU_BAR is treated as OSCU for the Moon, matching the
        # reference (a barycentric variant is meaningless for the
        # geocentric lunar orbit). NODBIT_MEAN takes precedence when set
        # (see prefer_mean above).
        if prefer_mean:
            node_body = MEAN_NODE
            apog_body = MEAN_APOG
        else:
            node_body = TRUE_NODE
            apog_body = OSCU_APOG

        # Node longitude from lunar theory
        node_pos, _ = calc_ut(jd_ut, node_body, calc_flags)
        node_lon = node_pos[0]
        node_lat = node_pos[1]

        # Apogee from lunar theory
        apog_pos, _ = calc_ut(jd_ut, apog_body, calc_flags)
        apog_lon = apog_pos[0]
        apog_lat = apog_pos[1]
        apog_dist = apog_pos[2]

        # Perigee: opposite point of the same orbit — longitude +180°,
        # latitude mirrored, with its physically distinct apsis distance.
        peri_lon = (apog_lon + 180.0) % 360.0
        peri_lat = -apog_lat
        if not prefer_mean:
            # Distance of the osculating perigee from the instantaneous
            # ellipse (lon/lat above already follow the apogee transform).
            from .lunar import calc_osculating_perigee

            peri_dist = calc_osculating_perigee(jd_tt)[2]
        else:
            # Mean ellipse identity: r_peri = 2a − r_apog with the textbook
            # mean lunar distance a = 0.002569555 AU (384,400 km).
            peri_dist = 2.0 * 0.002569555 - apog_dist

        # Node distances from the orbit ellipse: evaluate
        # r = a(1−e²)/(1+e·cos ν) at the node's true anomaly — the standard
        # conic-section radius. The node bodies' own distance channel is not
        # used.
        a_orb = 0.5 * (peri_dist + apog_dist)
        e_orb = (apog_dist - peri_dist) / (apog_dist + peri_dist)
        p_orb = a_orb * (1.0 - e_orb * e_orb)
        u_node = math.radians((node_lon - peri_lon) % 360.0)
        nasc_dist = p_orb / (1.0 + e_orb * math.cos(u_node))
        ndsc_dist = p_orb / (1.0 - e_orb * math.cos(u_node))

        def _moon_point(lon_d: float, lat_d: float, dist: float) -> PosTuple:
            """Assemble a lunar point, honoring HELCTR and J2000 output."""
            if iflag & (FLG_HELCTR | FLG_BARYCTR):
                # Measured reference behavior: FLG_HELCTR adds Earth's
                # heliocentric vector to the geocentric lunar point, and
                # FLG_BARYCTR its BARYCENTRIC vector (the point becomes
                # SSB-relative). With both bits set the heliocentric view
                # wins, mirroring the reference's nod_aps precedence.
                if iflag & FLG_HELCTR:
                    earth_center = r_earth_ecl
                else:
                    _e_b = earth_pos.position.au
                    earth_center = _icrs_to_ecliptic((_e_b[0], _e_b[1], _e_b[2]))
                cl = math.cos(math.radians(lat_d))
                gx = dist * cl * math.cos(math.radians(lon_d)) + earth_center[0]
                gy = dist * cl * math.sin(math.radians(lon_d)) + earth_center[1]
                gz = dist * math.sin(math.radians(lat_d)) + earth_center[2]
                dist = math.sqrt(gx**2 + gy**2 + gz**2)
                lon_d = math.degrees(math.atan2(gy, gx)) % 360.0
                lat_d = math.degrees(math.asin(max(-1.0, min(1.0, gz / dist))))
            if iflag & FLG_J2000:
                lon_d, lat_d = _nodaps_to_j2000(lon_d, lat_d, jd_tt)
            if (iflag & FLG_SIDEREAL) and not (iflag & FLG_EQUATORIAL):
                # Measured reference behavior: nod_aps sidereal subtracts the
                # MEAN ayanamsha (no nutation-in-longitude term), identical to
                # the planet branch below. Requesting NONUT semantics from
                # _get_ayanamsa_for_flags selects that mean value. FLG_EQUATORIAL
                # suppresses the ayanamsha entirely (SID+EQU == EQU alone): the
                # equatorial rotation runs afterward in _format_nodaps_output, so
                # subtracting here would double-shift the coordinate.
                aya = _get_ayanamsa_for_flags(jd_ut, (iflag & ~FLG_J2000) | FLG_NONUT)
                lon_d = float((lon_d - aya) % 360.0)
            return (lon_d, lat_d, dist, 0.0, 0.0, 0.0)

        # Build output: nodes and apsides from lunar theory
        moon_xnasc = _moon_point(node_lon, node_lat, nasc_dist)
        moon_xndsc = _moon_point((node_lon + 180.0) % 360.0, -node_lat, ndsc_dist)
        moon_xperi = _moon_point(peri_lon, peri_lat, peri_dist)
        if method & NODBIT_FOPOINT:
            # Second focus of the ellipse in the aphelion slot: same
            # direction as the apogee, at distance 2ae = r_apog − r_peri
            # from the Earth focus.
            moon_xaphe = _moon_point(apog_lon, apog_lat, apog_dist - peri_dist)
        else:
            moon_xaphe = _moon_point(apog_lon, apog_lat, apog_dist)

        return (moon_xnasc, moon_xndsc, moon_xperi, moon_xaphe)

    else:
        is_helio = bool(iflag & FLG_HELCTR)
        # Measured reference behavior for nod_aps: with BOTH center flags
        # set, the heliocentric view wins (unlike calc, where the
        # barycentric flag prevails — the reference itself differs between
        # the two surfaces and nod_aps mirrors its nod_aps contract).
        is_bary = bool(iflag & FLG_BARYCTR) and not is_helio
        is_centered = is_helio or is_bary
        # For the Sun the reference returns the apsides of the apparent
        # solar orbit: Earth's orbit mirrored through the origin.
        mirror_sun = ipl == SUN
        elem_ipl = EARTH if mirror_sun else ipl

        # NODBIT_MEAN (the default) uses the published mean-element
        # polynomials (Meeus Table 31.A / Simon et al. 1994); bodies
        # without mean elements (e.g. Pluto and the minor bodies) fall back to
        # osculating state vectors for every method, matching the reference.
        mean_el = None
        bar_mode = False
        if prefer_mean and not is_minor:
            from .planetary_mean_elements import mean_orbital_elements

            mean_el = mean_orbital_elements(elem_ipl, jd_tt)

    mean_frame = False
    if mean_el is not None:
        # --- Mean elements (mean ecliptic and equinox of date) ---
        a, e_mag, i_deg, om_deg, pi_deg = mean_el
        incl = math.radians(i_deg)
        Omega = math.radians(om_deg)
        omega = math.radians((pi_deg - om_deg) % 360.0)
        p = a * (1.0 - e_mag**2)
        mean_frame = True
    else:
        # --- Osculating elements from JPL state vectors ---
        # NODBIT_OSCU_BAR's public convention uses solar-system-barycentric
        # vectors for planets beyond Jupiter; other cases are heliocentric.
        bar_mode = bool(method & NODBIT_OSCU_BAR) and ipl in (
            SATURN,
            URANUS,
            NEPTUNE,
            PLUTO,
        )
        if mirror_sun or ipl == EARTH:
            # Earth's heliocentric orbit is the orbit of the Earth-Moon
            # barycentre, not of the geocentre (which additionally executes
            # the ~4700 km monthly wobble about that barycentre). Osculating
            # elements built from the wobbling geocentre state give the wrong
            # ellipse, so use the EMB state — the same source the Sun's
            # mirrored apparent orbit uses. Measured reference behavior: Earth
            # osculating apsides sit within ~0.02° of the mean ones (a smooth
            # EMB orbit), not the ~0.5° a geocentre osculation would produce.
            emb_pos = planets["earth barycenter"].at(t)
            r_icrs = emb_pos.position.au - sun_pos.position.au
            v_icrs = emb_pos.velocity.au_per_d - sun_pos.velocity.au_per_d
        elif is_minor:
            # Heliocentric geometric ICRS state from the body's own position
            # pipeline (SPK/LEB/Keplerian) — the same source calc() and
            # get_orbital_elements report for this body — reduced with GM_sun.
            r_helio, v_helio = _minor_body_helio_icrs_state(jd_tt, ipl)
            # A trans-jovian centaur/TNO re-references to the solar-system
            # barycentre under NODBIT_OSCU_BAR, exactly as the planets beyond
            # Jupiter do (measured, e.g. Chiron shifts ~0.07° while the belt
            # asteroids at a ~ 2.4-2.8 AU are unaffected). The barycentric state
            # is the heliocentric one plus the Sun's barycentric offset.
            if bool(method & NODBIT_OSCU_BAR) and (
                _osculating_semi_major_axis(
                    r_helio[0],
                    r_helio[1],
                    r_helio[2],
                    v_helio[0],
                    v_helio[1],
                    v_helio[2],
                    _GM_SUN,
                )
                > _TRANS_JOVIAN_A_AU
            ):
                bar_mode = True
                sp = sun_pos.position.au
                sv = sun_pos.velocity.au_per_d
                r_icrs = (
                    r_helio[0] + float(sp[0]),
                    r_helio[1] + float(sp[1]),
                    r_helio[2] + float(sp[2]),
                )
                v_icrs = (
                    v_helio[0] + float(sv[0]),
                    v_helio[1] + float(sv[1]),
                    v_helio[2] + float(sv[2]),
                )
            else:
                r_icrs = r_helio
                v_icrs = v_helio
        elif bar_mode:
            # SSB-relative vectors for the planets beyond Jupiter.
            r_icrs = target_pos.position.au
            v_icrs = target_pos.velocity.au_per_d
        else:
            r_icrs = target_pos.position.au - sun_pos.position.au
            v_icrs = target_pos.velocity.au_per_d - sun_pos.velocity.au_per_d
        mass_ratio = _SUN_MASS_RATIOS.get(elem_ipl)
        GM = _GM_SUN * (1.0 + 1.0 / mass_ratio) if mass_ratio else _GM_SUN

        # Convert to the output ecliptic frame (of date, or J2000)
        r_ecl = _icrs_to_ecliptic(r_icrs)
        v_ecl = _icrs_to_ecliptic(v_icrs)

        r_mag = math.sqrt(r_ecl[0] ** 2 + r_ecl[1] ** 2 + r_ecl[2] ** 2)
        v_mag = math.sqrt(v_ecl[0] ** 2 + v_ecl[1] ** 2 + v_ecl[2] ** 2)

        # Angular momentum vector h = r × v
        hx = r_ecl[1] * v_ecl[2] - r_ecl[2] * v_ecl[1]
        hy = r_ecl[2] * v_ecl[0] - r_ecl[0] * v_ecl[2]
        hz = r_ecl[0] * v_ecl[1] - r_ecl[1] * v_ecl[0]
        h_mag = math.sqrt(hx**2 + hy**2 + hz**2)

        # Inclination
        incl = math.acos(max(-1.0, min(1.0, hz / h_mag))) if h_mag > 0 else 0.0

        # Node vector n = k × h (k = ecliptic pole = [0, 0, 1])
        nx = -hy
        ny = hx
        n_mag = math.sqrt(nx**2 + ny**2)

        # Longitude of ascending node
        if n_mag > 1e-10:
            Omega = math.atan2(ny, nx)
            if Omega < 0:
                Omega += 2.0 * math.pi
        else:
            Omega = 0.0

        # Eccentricity vector e = (v × h) / μ - r̂
        r_dot_v = r_ecl[0] * v_ecl[0] + r_ecl[1] * v_ecl[1] + r_ecl[2] * v_ecl[2]
        coef1 = v_mag**2 / GM - 1.0 / r_mag
        coef2 = r_dot_v / GM

        ex = coef1 * r_ecl[0] - coef2 * v_ecl[0]
        ey = coef1 * r_ecl[1] - coef2 * v_ecl[1]
        ez = coef1 * r_ecl[2] - coef2 * v_ecl[2]
        e_mag = math.sqrt(ex**2 + ey**2 + ez**2)

        # Argument of perihelion (or perigee for Moon)
        if n_mag > 1e-10 and e_mag > 1e-10:
            n_dot_e = nx * ex + ny * ey
            cos_omega = max(-1.0, min(1.0, n_dot_e / (n_mag * e_mag)))
            omega = math.acos(cos_omega)
            if ez < 0:
                omega = 2.0 * math.pi - omega
        else:
            omega = 0.0

        # Semi-major axis
        a = 1.0 / (2.0 / r_mag - v_mag**2 / GM)

        # Semi-latus rectum
        p = a * (1.0 - e_mag**2) if e_mag < 1.0 else h_mag**2 / GM

    # --- Compute 3D positions of nodes and apsides in orbital frame ---
    cos_Omega = math.cos(Omega)
    sin_Omega = math.sin(Omega)
    cos_incl = math.cos(incl)
    sin_incl = math.sin(incl)
    cos_omega = math.cos(omega)
    sin_omega = math.sin(omega)

    def _orbit_pos_3d(nu: float):
        """Compute ecliptic 3D position for given true anomaly.

        For planets: heliocentric ecliptic position.
        For Moon: geocentric ecliptic position.
        """
        denom = 1.0 + e_mag * math.cos(nu)
        r_orb = p / denom if abs(denom) > 1e-10 else a
        x_orb = r_orb * math.cos(nu)
        y_orb = r_orb * math.sin(nu)
        # Perifocal → ecliptic rotation (3-1-3: Omega, incl, omega)
        xe = (cos_Omega * cos_omega - sin_Omega * sin_omega * cos_incl) * x_orb + (
            -cos_Omega * sin_omega - sin_Omega * cos_omega * cos_incl
        ) * y_orb
        ye = (sin_Omega * cos_omega + cos_Omega * sin_omega * cos_incl) * x_orb + (
            -sin_Omega * sin_omega + cos_Omega * cos_omega * cos_incl
        ) * y_orb
        ze = (sin_omega * sin_incl) * x_orb + (cos_omega * sin_incl) * y_orb
        return (xe, ye, ze, r_orb)

    def _focal_point_3d():
        """Compute the second focal point of the orbit ellipse.

        The second (empty) focus is located at distance 2ae from the
        primary focus (Sun or Earth) along the apse line, in the
        anti-perihelion direction. In the perifocal frame this is the
        point (-2ae, 0, 0).

        This matches the reference convention which returns the second
        focal point in the aphelion/apogee slot of nod_aps results.
        """
        # Distance from primary focus to second focus = 2 * a * e
        f_dist = 2.0 * a * e_mag
        # In perifocal frame, the second focus is at (-f_dist, 0, 0)
        # (negative x = anti-perihelion direction)
        fx_orb = -f_dist
        # Rotate to ecliptic
        fxe = (cos_Omega * cos_omega - sin_Omega * sin_omega * cos_incl) * fx_orb
        fye = (sin_Omega * cos_omega + cos_Omega * sin_omega * cos_incl) * fx_orb
        fze = (sin_omega * sin_incl) * fx_orb
        return (fxe, fye, fze, f_dist)

    # --- Frame alignment for mean-element points -------------------------
    # Mean elements refer to the mean equinox of date; the of-date output
    # frame is the true equinox of date, so rotate by the nutation in
    # longitude unless NONUT. (FLG_J2000 precesses the finished of-date
    # coordinates at the end, nutation included — reference behavior.)
    if mean_frame and not (iflag & FLG_NONUT):
        from .cache import get_cached_nutation

        _dpsi_rad, _ = get_cached_nutation(jd_tt)
        _cos_dpsi = math.cos(_dpsi_rad)
        _sin_dpsi = math.sin(_dpsi_rad)

        def _align_mean(vec):
            return (
                vec[0] * _cos_dpsi - vec[1] * _sin_dpsi,
                vec[0] * _sin_dpsi + vec[1] * _cos_dpsi,
                vec[2],
            )
    else:

        def _align_mean(vec):
            return vec

    # Observer barycentric velocity in the output frame, for the annual
    # aberration the reference applies to these points unless suppressed.
    apply_aberr = (
        not is_centered and not mirror_sun and not (iflag & (FLG_TRUEPOS | FLG_NOABERR))
    )
    if apply_aberr:
        from .astrometry import apply_aberration_to_position

        _v_obs = earth_pos.velocity.au_per_d
        v_obs_ecl = _icrs_to_ecliptic((_v_obs[0], _v_obs[1], _v_obs[2]))

    apply_deflection = (
        not is_centered and not mirror_sun and not (iflag & (FLG_TRUEPOS | FLG_NOGDEFL))
    )
    if apply_deflection:
        from .fast_calc import C_LIGHT_AU_DAY, _apply_gravitational_deflection

        earth_barycentric = tuple(float(value) for value in earth_pos.position.au)
        sun_barycentric = tuple(float(value) for value in sun_pos.position.au)
        deflector_source = _SkyfieldDeflectorSource(planets)

    # Observer position vector in the working frame. For OSCU_BAR the
    # points are SSB-relative, so the observer is the barycentric Earth
    # (or the barycentric Sun under FLG_HELCTR) for frame consistency;
    # otherwise the heliocentric Earth (None = no subtraction, points are
    # already observer-relative).
    if bar_mode and is_bary:
        # OSCU_BAR points are already SSB-relative; a barycentric request
        # needs no observer subtraction.
        r_obs_ecl = None
    elif bar_mode and is_centered:
        _s_bary = sun_pos.position.au
        r_obs_ecl = _icrs_to_ecliptic((_s_bary[0], _s_bary[1], _s_bary[2]))
    elif bar_mode:
        _e_bary = earth_pos.position.au
        r_obs_ecl = _icrs_to_ecliptic((_e_bary[0], _e_bary[1], _e_bary[2]))
    elif is_bary:
        # FLG_BARYCTR on heliocentric orbit points: measured reference
        # behavior re-references each node/apse point to the solar-system
        # barycentre (the point vector gains the Sun's barycentric offset),
        # so the observer sits at -r_sun in the Sun-centered working frame.
        _s_bary = sun_pos.position.au
        r_obs_ecl = _icrs_to_ecliptic((-_s_bary[0], -_s_bary[1], -_s_bary[2]))
    elif is_centered:
        r_obs_ecl = None
    else:
        r_obs_ecl = r_earth_ecl

    def _to_geo_lonlat(center_pos):
        """Convert an orbit-frame point to output lon/lat/dist.

        Planets: heliocentric point → geocentric (subtract Earth) unless
        FLG_HELCTR. Sun: the apparent solar orbit is Earth's mirrored
        through the origin, so the geocentric point is the negated
        heliocentric one. Aberration is applied per the flag gate above.
        """
        px, py, pz = _align_mean((center_pos[0], center_pos[1], center_pos[2]))
        if mirror_sun and not is_centered:
            # Geocentric solar-orbit point = heliocentric Earth-orbit
            # point mirrored through the origin. With HELCTR the
            # reference returns the Earth-orbit point itself.
            gx, gy, gz = -px, -py, -pz
        elif r_obs_ecl is None:
            gx, gy, gz = px, py, pz
        else:
            gx = px - r_obs_ecl[0]
            gy = py - r_obs_ecl[1]
            gz = pz - r_obs_ecl[2]
        if apply_deflection:
            geo_icrs = _ecliptic_to_icrs((gx, gy, gz))
            if not _nodaps_ray_crosses_solar_disk(
                geo_icrs, earth_barycentric, sun_barycentric
            ):
                distance = math.sqrt(gx * gx + gy * gy + gz * gz)
                geo_icrs = _apply_gravitational_deflection(
                    geo_icrs,
                    earth_barycentric,
                    jd_tt,
                    distance / C_LIGHT_AU_DAY,
                    deflector_source,
                )
                gx, gy, gz = _icrs_to_ecliptic(geo_icrs)
        if apply_aberr:
            gx, gy, gz = apply_aberration_to_position((gx, gy, gz), v_obs_ecl)
        r_geo = math.sqrt(gx**2 + gy**2 + gz**2)
        lon = math.degrees(math.atan2(gy, gx)) % 360.0
        lat = (
            math.degrees(math.asin(max(-1.0, min(1.0, gz / r_geo))))
            if r_geo > 0
            else 0.0
        )
        return lon, lat, r_geo

    # Ascending node: true anomaly = -omega (body crosses ecliptic northward)
    pos_asc = _orbit_pos_3d(-omega)
    # Descending node: true anomaly = pi - omega
    pos_dsc = _orbit_pos_3d(math.pi - omega)
    # Perihelion/perigee: true anomaly = 0
    pos_peri = _orbit_pos_3d(0.0)
    # Aphelion: true anomaly = π (farthest point from Sun)
    # With NODBIT_FOPOINT, return the second focal point instead
    if method & NODBIT_FOPOINT:
        pos_aphe = _focal_point_3d()
    else:
        pos_aphe = _orbit_pos_3d(math.pi)

    # Convert to output ecliptic coordinates
    if mirror_sun or ipl == EARTH:
        # The solar/terrestrial orbit lies in the ecliptic plane by
        # construction (Earth's orbit defines the ecliptic), so its nodes are
        # undefined. Measured reference behavior returns zeros for both node
        # slots of the Sun and of Earth while still reporting the apsides.
        geo_asc = (0.0, 0.0, 0.0)
        geo_dsc = (0.0, 0.0, 0.0)
    else:
        geo_asc = _to_geo_lonlat(pos_asc)
        geo_dsc = _to_geo_lonlat(pos_dsc)
    geo_peri = _to_geo_lonlat(pos_peri)
    geo_aphe = _to_geo_lonlat(pos_aphe)

    # J2000 output. Measured reference behavior: the FLG_J2000 node/apse are
    # NOT a rigid frame rotation of the finished of-date geocentric point. The
    # reference precesses only the ecliptic LONGITUDE of the heliocentric
    # node/apse (a rotation about the ecliptic pole by the accumulated general
    # precession in longitude), keeping the point's of-date ecliptic latitude,
    # then re-projects it geocentrically against the Earth expressed in the
    # J2000 ecliptic. Because a point referred to the of-date plane is projected
    # against an Earth that carries the of-date-to-J2000 ecliptic tilt, the
    # J2000 node/apse acquire a distance-dependent (parallactic) latitude --
    # near zero at J2000, up to ~130" of latitude shift at |t-J2000| ~ 1
    # century -- that a rotation of the geocentric direction cannot reproduce
    # (a rotation is distance-independent, but the measured shift scales with
    # 1/geocentric-distance). Algebraically the observed map is
    #     geo_j2000 = Rz(delta) . geo_of_date + C,
    #     C = Rz(delta) . r_obs_of_date - r_obs_j2000,
    # where r_obs is the observer subtracted in the geocentric projection above
    # and delta is the ecliptic longitude of the of-date equinox measured in the
    # J2000 ecliptic. C vanishes for a heliocentric/centred request (no observer
    # subtraction) and for the Sun's mirrored apparent orbit, which reduce to a
    # pure longitude precession.
    if iflag & FLG_J2000:
        # The longitude precession and the observer's J2000 offset are both
        # referred to the MEAN ecliptic of date and of J2000 (no nutation): the
        # measured reference longitude reduction is the accumulated general
        # precession only, so folding the of-date nutation into the rotation
        # would mis-precess the longitude by the nutation-in-longitude term
        # (~16" near 1900). The of-date node/apse coordinates themselves retain
        # their nutation; only the frame-to-frame rotation is mean.
        _pn_mean = vondrak_precession_matrix(jd_tt)
        _eps_mean = vondrak_mean_obliquity_rad(jd_tt)
        _cm, _sm = math.cos(_eps_mean), math.sin(_eps_mean)

        def _mean_ecl_to_icrs(vec):
            xq = vec[0]
            yq = vec[1] * _cm - vec[2] * _sm
            zq = vec[1] * _sm + vec[2] * _cm
            return (
                _pn_mean[0][0] * xq + _pn_mean[1][0] * yq + _pn_mean[2][0] * zq,
                _pn_mean[0][1] * xq + _pn_mean[1][1] * yq + _pn_mean[2][1] * zq,
                _pn_mean[0][2] * xq + _pn_mean[1][2] * yq + _pn_mean[2][2] * zq,
            )

        def _to_ecl_j2000(vec_icrs):
            pj = vondrak_precession_matrix(_J2000_JD)
            eps_j = vondrak_mean_obliquity_rad(_J2000_JD)
            cj, sj = math.cos(eps_j), math.sin(eps_j)
            xq = (
                pj[0][0] * vec_icrs[0] + pj[0][1] * vec_icrs[1] + pj[0][2] * vec_icrs[2]
            )
            yq = (
                pj[1][0] * vec_icrs[0] + pj[1][1] * vec_icrs[1] + pj[1][2] * vec_icrs[2]
            )
            zq = (
                pj[2][0] * vec_icrs[0] + pj[2][1] * vec_icrs[1] + pj[2][2] * vec_icrs[2]
            )
            return (xq, yq * cj + zq * sj, -yq * sj + zq * cj)

        _eqx = _to_ecl_j2000(_mean_ecl_to_icrs((1.0, 0.0, 0.0)))
        _delta = math.atan2(_eqx[1], _eqx[0])
        _cd, _sd = math.cos(_delta), math.sin(_delta)
        if r_obs_ecl is None or mirror_sun:
            _cvec = (0.0, 0.0, 0.0)
        else:
            _obs_j = _to_ecl_j2000(_mean_ecl_to_icrs(r_obs_ecl))
            _cvec = (
                _cd * r_obs_ecl[0] - _sd * r_obs_ecl[1] - _obs_j[0],
                _sd * r_obs_ecl[0] + _cd * r_obs_ecl[1] - _obs_j[1],
                r_obs_ecl[2] - _obs_j[2],
            )

        def _prec_j2000(g):
            if g[2] == 0.0:
                return g
            cl = math.cos(math.radians(g[1]))
            vx = g[2] * cl * math.cos(math.radians(g[0]))
            vy = g[2] * cl * math.sin(math.radians(g[0]))
            vz = g[2] * math.sin(math.radians(g[1]))
            rx = _cd * vx - _sd * vy + _cvec[0]
            ry = _sd * vx + _cd * vy + _cvec[1]
            rz = vz + _cvec[2]
            rr = math.sqrt(rx * rx + ry * ry + rz * rz)
            lon_j = math.degrees(math.atan2(ry, rx)) % 360.0
            lat_j = math.degrees(math.asin(max(-1.0, min(1.0, rz / rr))))
            return (lon_j, lat_j, rr)

        geo_asc = _prec_j2000(geo_asc)
        geo_dsc = _prec_j2000(geo_dsc)
        geo_peri = _prec_j2000(geo_peri)
        geo_aphe = _prec_j2000(geo_aphe)

    # Sidereal output: subtract the ayanamsha from the longitudes,
    # consistent with calc_ut (issue #29: the flag was silently ignored
    # for planetary nodes/apsides).
    #
    # Measured reference behavior: nod_aps sidereal longitudes are formed by
    # subtracting the MEAN ayanamsha (mean + no Δψ), even when the of-date
    # coordinate itself carries nutation. Subtracting the true ayanamsha
    # instead left a residual of the nutation-in-longitude term (~13.9" at
    # J2000). Requesting NONUT semantics selects the mean value; FLG_J2000 is
    # stripped so the ayanamsha epoch matches the of-date longitude.
    #
    # FLG_EQUATORIAL suppresses the sidereal correction entirely: measured
    # reference behavior returns pure equatorial coordinates for SID+EQU
    # (identical to EQU alone). _format_nodaps_output rotates the of-date
    # longitude to the equator afterward, so subtracting the ayanamsha here
    # would double-shift the coordinate (mirrors calc_ut, where EQUATORIAL
    # skips the sidereal reduction).
    if (iflag & FLG_SIDEREAL) and not (iflag & FLG_EQUATORIAL):
        aya = _get_ayanamsa_for_flags(t.ut1, (iflag & ~FLG_J2000) | FLG_NONUT)

        def _to_sidereal(g):
            # Skip the zero sentinel (e.g. the Sun's node slots, which the
            # reference returns as zeros): subtracting the ayanamsha would
            # fabricate a spurious node longitude — mirror the J2000 guard.
            if g[2] == 0.0:
                return g
            # float() guards against a numpy.float64 ayanamsha leaking into the
            # output tuple (the library contract is native Python floats).
            return (float((g[0] - aya) % 360.0), g[1], g[2])

        geo_asc = _to_sidereal(geo_asc)
        geo_dsc = _to_sidereal(geo_dsc)
        geo_peri = _to_sidereal(geo_peri)
        geo_aphe = _to_sidereal(geo_aphe)

    # Build output tuples (lon, lat, dist, speed_lon, speed_lat, speed_dist)
    xnasc: PosTuple = (geo_asc[0], geo_asc[1], geo_asc[2], 0.0, 0.0, 0.0)
    xndsc: PosTuple = (geo_dsc[0], geo_dsc[1], geo_dsc[2], 0.0, 0.0, 0.0)
    xperi: PosTuple = (geo_peri[0], geo_peri[1], geo_peri[2], 0.0, 0.0, 0.0)
    xaphe: PosTuple = (geo_aphe[0], geo_aphe[1], geo_aphe[2], 0.0, 0.0, 0.0)

    return (xnasc, xndsc, xperi, xaphe)


def get_orbital_elements(tjdet: float, planet: int, flags: int) -> Tuple[float, ...]:
    """
    Calculate Keplerian orbital elements for a celestial body.

    Reference API compatible function matching the reference signature.

    This function computes the osculating (instantaneous) orbital elements
    for a planet at a given time. The elements describe the elliptical orbit
    that the planet would follow if all perturbations ceased at that moment.

    Args:
        tjdet: Julian Day in Ephemeris Time (TT/ET)
        planet: Planet/body ID (SUN, MOON, etc.)
        flags: Calculation flags (FLG_HELCTR for heliocentric, etc.)

    Returns:
        Flat tuple of 50 floats with orbital elements:
            [0] a: Semi-major axis (AU)
            [1] e: Eccentricity (0=circle, <1=ellipse, 1=parabola, >1=hyperbola)
            [2] i: Inclination (degrees, relative to ecliptic)
            [3] Omega: Longitude of ascending node (degrees)
            [4] omega: Argument of perihelion (degrees)
            [5] varpi: Longitude of periapsis (degrees)
            [6] M: Mean anomaly at epoch (degrees)
            [7] nu: True anomaly at epoch (degrees)
            [8] E: Eccentric anomaly at epoch (degrees)
            [9] L: Mean longitude at epoch (degrees)
            [10] P_sid: Sidereal orbital period (tropical years)
            [11] n: Mean daily motion (degrees/day)
            [12] P_trop: Tropical period (years)
            [13] P_syn: Synodic period (days, negative for inner planets/Moon)
            [14] T: Time of perihelion passage (JD)
            [15] q: Perihelion distance (AU)
            [16] Q: Aphelion distance (AU)
            [17-49]: Reserved (0.0)

    Example:
        >>> from libephemeris import get_orbital_elements, MARS, FLG_HELCTR
        >>> elements = get_orbital_elements(2451545.0, MARS, FLG_HELCTR)
        >>> a, e, i = elements[0], elements[1], elements[2]
        >>> print(f"Mars: a={a:.4f} AU, e={e:.4f}, i={i:.4f}")

    Note:
        - For heliocentric calculations (default), elements are relative to the Sun
        - With FLG_BARYCTR the elements of the trans-jovian bodies (Saturn,
          Uranus, Neptune, Pluto, and centaurs/TNOs beyond Jupiter) are
          re-referenced to the solar-system barycentre; Jupiter and every
          interior body are unchanged (matching the reference)
        - Moon's elements are geocentric (relative to Earth)
        - Elements change constantly due to perturbations from other planets
        - With FLG_ORBEL_AA the osculating fit uses the Astronomical Almanac
          central mass (Sun plus every major planet at or interior to the
          body's orbit, Earth as the Earth+Moon barycentre): the semi-major
          axis shrinks and the perihelion-referred angles shift solidly while
          the mean longitude is left invariant. In this context bit 32768 is
          always read as FLG_ORBEL_AA (it shares its value with FLG_TOPOCTR,
          which has no meaning for orbital elements)

    Divergences (documented, intentional): the Sun, the lunar nodes/apogees
    (10-13) and the interpolated apsides (21, 22) raise ``object N not valid``;
    ids 23-39 raise ``illegal planet number``; negatives and the sourceless
    block between the last fictitious body (58) and the numbered-asteroid
    offset raise ``object N not valid`` (the reference rejects the same ids
    through its element-file lookup). The fictitious/hypothetical bodies (40-58) are
    reduced from the library's *own* runtime models rather than a foreign
    element set, so where the two models disagree by arbitration the elements
    differ (e.g. Isis/Transpluto a=77.755 here vs 77.775; the White Moon (56)
    and Waldemath (58) are geocentric circular constructions here, so their
    elements are Earth-relative). Nibiru (49) is intentionally not modeled and
    the position pipeline's typed error propagates.
    """
    _validate_orbital_object("get_orbital_elements", planet, sun_valid=False)
    ts = get_timescale()
    t = ts.tt_jd(tjdet)
    return _calc_orbital_elements(t, planet, flags)


def get_orbital_elements_ut(tjd_ut: float, ipl: int, iflag: int) -> Tuple[float, ...]:
    """
    Calculate Keplerian orbital elements for Universal Time.

    Reference API compatible function. Similar to get_orbital_elements()
    but takes Universal Time instead of Ephemeris Time.

    Args:
        tjd_ut: Julian Day in Universal Time (UT1)
        ipl: Planet/body ID (SUN, MOON, etc.)
        iflag: Calculation flags

    Returns:
        Flat tuple of 50 floats with orbital elements (same as get_orbital_elements).

    Example:
        >>> elements = get_orbital_elements_ut(2451545.0, JUPITER, 0)
        >>> print(f"Jupiter semi-major axis: {elements[0]:.4f} AU")
    """
    _validate_orbital_object("get_orbital_elements_ut", ipl, sun_valid=False)
    ts = get_timescale()
    t = ts.ut1_jd(tjd_ut)
    return _calc_orbital_elements(t, ipl, iflag)


# Semi-major-axis threshold for FLG_BARYCTR orbital elements. Measured
# reference behavior: under FLG_BARYCTR the two-body osculating fit is
# re-referenced from the Sun to the solar-system barycentre only for orbits
# beyond Jupiter; Jupiter (a ~ 5.20 AU) and every interior body (the terrestrial
# planets and the main-belt asteroids Ceres/Pallas/Juno/Vesta at ~2.4-2.8 AU)
# are unchanged. The 6.0 AU boundary cleanly separates Jupiter's stable
# semi-major axis from Saturn (~9.58 AU) and the trans-jovian centaurs/TNOs.
_TRANS_JOVIAN_A_AU = 6.0


def _osculating_semi_major_axis(
    rx: float, ry: float, rz: float, vx: float, vy: float, vz: float, GM: float
) -> float:
    """Two-body semi-major axis (AU) from a state vector and central GM."""
    r = math.sqrt(rx * rx + ry * ry + rz * rz)
    if not (r > 0.0):
        return 0.0
    v2 = vx * vx + vy * vy + vz * vz
    denom = 2.0 / r - v2 / GM
    return (1.0 / denom) if denom != 0.0 else 0.0


def _sun_barycentric_ecliptic_state(
    t,
) -> Tuple[float, float, float, float, float, float]:
    """Sun position/velocity relative to the SSB, in the mean J2000 ecliptic.

    Adding this to a heliocentric state (target - Sun) yields the barycentric
    state (target - SSB): (target - Sun) + (Sun - SSB) = target - SSB.
    """
    from .fast_calc import _rotate_icrs_to_ecliptic_j2000

    planets = _get_computation_ephemeris()
    sun = planets["sun"].at(t)
    sp = sun.position.au
    sv = sun.velocity.au_per_d
    sx, sy, sz = _rotate_icrs_to_ecliptic_j2000(
        float(sp[0]), float(sp[1]), float(sp[2])
    )
    svx, svy, svz = _rotate_icrs_to_ecliptic_j2000(
        float(sv[0]), float(sv[1]), float(sv[2])
    )
    return (sx, sy, sz, svx, svy, svz)


def _ecliptic_spherical_to_cartesian_state(
    pos: Tuple[float, float, float, float, float, float],
) -> Tuple[float, float, float, float, float, float]:
    """Ecliptic spherical state -> Cartesian state, both in the same frame.

    Converts ``(lon, lat, r, dlon, dlat, ddist)`` in degrees / degrees-per-day
    / AU (the calc() spherical convention) to a Cartesian ``(x, y, z, vx, vy,
    vz)`` state in AU and AU/day. The velocity is the analytic time-derivative
    of the position, so the pair is self-consistent for a two-body reduction.
    """
    lon, lat, r, dlon, dlat, ddist = pos
    lam = math.radians(lon)
    bet = math.radians(lat)
    dlam = math.radians(dlon)
    dbet = math.radians(dlat)
    cos_lam = math.cos(lam)
    sin_lam = math.sin(lam)
    cos_bet = math.cos(bet)
    sin_bet = math.sin(bet)

    x = r * cos_bet * cos_lam
    y = r * cos_bet * sin_lam
    z = r * sin_bet
    vx = (
        ddist * cos_bet * cos_lam
        - r * sin_bet * cos_lam * dbet
        - r * cos_bet * sin_lam * dlam
    )
    vy = (
        ddist * cos_bet * sin_lam
        - r * sin_bet * sin_lam * dbet
        + r * cos_bet * cos_lam * dlam
    )
    vz = ddist * sin_bet + r * cos_bet * dbet
    return (x, y, z, vx, vy, vz)


def _minor_body_helio_icrs_state(
    jd_tt: float, ipl: int
) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
    """Heliocentric geometric ICRS state (AU, AU/day) for a minor body.

    Sources the state from the same position pipeline as calc() and
    get_orbital_elements: a geometric heliocentric place in the mean J2000
    ecliptic (FLG_HELCTR | FLG_J2000 | FLG_TRUEPOS: no light-time, aberration,
    or deflection). The state is rotated from the mean J2000 ecliptic back to
    ICRS so it enters the shared node/apse reduction on exactly the same footing
    as the Skyfield planet states (which the reduction rotates ICRS -> ecliptic
    of date). FLG_SPEED is always requested because the osculating ellipse needs
    the velocity even when the caller did not ask for node/apse speeds.

    Raises:
        Error: If the pipeline returns a degenerate (zero/non-finite) distance.
            Range/coverage/unknown-body failures surface as the pipeline's own
            typed Error subclass (never a silent zero state).
    """
    from .exceptions import Error
    from .fast_calc import _ICRS_TO_ECL_J2000_REF

    pos, _ = calc(jd_tt, ipl, FLG_HELCTR | FLG_J2000 | FLG_TRUEPOS | FLG_SPEED)
    if not (pos[2] > 0.0) or not math.isfinite(pos[2]):
        raise Error(
            f"nod_aps: the position pipeline returned a degenerate state "
            f"(distance {pos[2]} AU) for body {ipl} at JD {jd_tt}."
        )
    x, y, z, vx, vy, vz = _ecliptic_spherical_to_cartesian_state(pos)

    # Mean J2000 ecliptic -> ICRS is the transpose of the ICRS -> ecliptic
    # frame matrix used to define the FLG_J2000 ecliptic output above.
    m = _ICRS_TO_ECL_J2000_REF

    def _to_icrs(ax: float, ay: float, az: float) -> Tuple[float, float, float]:
        return (
            m[0][0] * ax + m[1][0] * ay + m[2][0] * az,
            m[0][1] * ax + m[1][1] * ay + m[2][1] * az,
            m[0][2] * ax + m[1][2] * ay + m[2][2] * az,
        )

    return _to_icrs(x, y, z), _to_icrs(vx, vy, vz)


def _calc_orbital_elements_minor(t, ipl: int, iflag: int = 0) -> Tuple[float, ...]:
    """Osculating orbital elements for a minor body (asteroid/centaur/TNO).

    The heliocentric state is taken from the same position pipeline that
    serves calc/calc_ut for this body (SPK, LEB, or the documented Keplerian
    fallback), requested as a *geometric* place (FLG_TRUEPOS: no light-time,
    aberration, or deflection) in the J2000 ecliptic (FLG_HELCTR | FLG_J2000).
    The resulting elements are therefore consistent with the positions the
    library reports for the body, whatever backend is active. The state is
    then reduced by the same state->elements math used for the planets.

    Raises:
        Error: If no usable heliocentric state can be obtained for the body
            (e.g. outside SPK coverage with no fallback, or an unknown
            asteroid number). Never returns a silent zero-element tuple.
    """
    GM = 0.01720209895**2  # Gaussian gravitational constant k^2 (IAU)
    # Geometric heliocentric J2000-ecliptic state (light-time/aberration/
    # deflection removed) reduced with GM_sun. Under FLG_BARYCTR a trans-jovian
    # minor body (e.g. a centaur/TNO beyond Jupiter) is re-referenced to the SSB.
    flags = FLG_HELCTR | FLG_J2000 | FLG_TRUEPOS | FLG_SPEED
    return _osculating_from_calc_state(
        t,
        ipl,
        flags,
        GM,
        want_barycentric=bool(iflag & FLG_BARYCTR),
        orbel_aa=bool(iflag & FLG_ORBEL_AA),
    )


def _osculating_from_calc_state(
    t,
    ipl: int,
    flags: int,
    GM: float,
    want_barycentric: bool = False,
    orbel_aa: bool = False,
) -> Tuple[float, ...]:
    """Osculating elements from the geometric state the pipeline serves.

    Requests ``calc(..., flags)`` (expected to be a geometric J2000-ecliptic
    place, heliocentric or geocentric depending on ``flags``), converts the
    spherical state to a Cartesian state vector, and reduces it with the shared
    state->elements math using the supplied central-body ``GM``. The position
    pipeline raises a typed Error subclass (UnknownBodyError,
    EphemerisRangeError, SPKRequiredError, ...) when it cannot serve the body;
    that clear message is allowed to propagate rather than being masked.

    Args:
        t: Skyfield Time object at the osculating epoch.
        ipl: Body id.
        flags: calc() flags selecting the frame/centre of the state.
        GM: Gravitational parameter of the central body (AU^3/day^2).

    Returns:
        The shared 50-slot osculating-element tuple.

    Raises:
        Error: If the pipeline returns a degenerate (zero/non-finite) state.
    """
    from .exceptions import Error

    tt_jd = float(t.tt)
    pos, _ = calc(tt_jd, ipl, flags)
    r = pos[2]

    if not (r > 0.0) or not math.isfinite(r):
        raise Error(
            f"orbital elements unavailable for body {ipl}: the position "
            f"pipeline returned a degenerate state "
            f"(distance {r} AU) at JD {tt_jd}."
        )

    # Ecliptic spherical (degrees, degrees/day) -> Cartesian state (AU, AU/day).
    x, y, z, vx, vy, vz = _ecliptic_spherical_to_cartesian_state(pos)

    # Barycentric re-reference for trans-jovian orbits: shift the heliocentric
    # state (target - Sun) to the barycentric state (target - SSB) by adding the
    # Sun's barycentric position/velocity. Only orbits beyond Jupiter respond
    # (measured), so the belt asteroids and any interior body stay heliocentric.
    if (
        want_barycentric
        and _osculating_semi_major_axis(x, y, z, vx, vy, vz, GM) > _TRANS_JOVIAN_A_AU
    ):
        sx, sy, sz, svx, svy, svz = _sun_barycentric_ecliptic_state(t)
        x, y, z = x + sx, y + sy, z + sz
        vx, vy, vz = vx + svx, vy + svy, vz + svz

    GM_slots: float | None = None
    if orbel_aa:
        # Astronomical Almanac convention (FLG_ORBEL_AA): reduce the osculating
        # geometry with the Sun plus every major planet interior to this body's
        # orbit as the central mass. The interior set is chosen from the default
        # (solar-GM) osculating a, which is ample to place the body between
        # planet orbits. The derived period family keeps the pure solar GM
        # (GM_slots), so it follows only through the changed a. In the orbital-
        # elements context bit 32768 always means ORBEL_AA, never TOPOCTR.
        a_default = _osculating_semi_major_axis(x, y, z, vx, vy, vz, GM)
        GM_slots = GM
        GM = _GM_SUN * (1.0 + _orbel_aa_interior_mass_sum(a_default))

    return _orbital_elements_from_ecliptic_state(
        x, y, z, vx, vy, vz, GM, ipl, tt_jd, GM_slots=GM_slots
    )


def _calc_orbital_elements_fictitious(
    t, ipl: int, orbel_aa: bool = False
) -> Tuple[float, ...]:
    """Osculating elements for a fictitious/hypothetical body (40-58).

    The elements are reduced from the library's *own* runtime model for the
    body, using the same state->elements machinery as the minor bodies. Most
    fictitious constructions are heliocentric (Uranians 40-47, Transpluto 48,
    the predicted planets 50-54, Vulcan 55, Proserpina 57) and use GM_sun. The
    two geocentric symbolic points -- White Moon Selena (56) and the Waldemath
    dark moon (58) -- orbit the Earth, so their osculating ellipse is taken
    from the geometric *geocentric* J2000-ecliptic state with GM_earth. Bodies
    the library does not model (e.g. Nibiru 49) raise UnknownBodyError from the
    position pipeline, which propagates.

    ``orbel_aa`` carries FLG_ORBEL_AA through to the shared reducer, exactly
    as the minor-body path does. Dropping it here made the flag a no-op for
    every fictitious id, while the heliocentric constructions sit beyond
    several planetary orbits and therefore have a genuinely different
    interior-mass sum. The geocentric symbolic points keep GM_earth: the
    Astronomical Almanac interior-mass convention is defined for
    heliocentric orbits, and no planet lies interior to an Earth satellite.
    """
    if ipl in _FICT_GEOCENTRIC_IDS:
        GM_earth = _GM_SUN / _SUN_EARTH_MASS_RATIO
        flags = FLG_J2000 | FLG_TRUEPOS | FLG_SPEED  # geocentric geometric place
        return _osculating_from_calc_state(t, ipl, flags, GM_earth)
    if ipl in _FICT_ECL_OF_DATE_CIRCULAR_IDS:
        # Reduce this exactly-circular of-date model in its native mean ecliptic
        # of date. A J2000 request rotates the position onto the of-date-circular
        # path but injects the equinox precession rate into the tangential
        # velocity, a ~2.7% deficit against the Gaussian circular speed
        # (n = k/a^1.5) that the reducer would read as a spurious e ~ 0.054
        # aphelion. In the mean-of-date frame (FLG_NONUT, no FLG_J2000) position
        # and velocity are self-consistent, so the circular orbit yields e ~ 0
        # with a = q = Q. The angular elements are degenerate for e = i = 0, so
        # the of-date reference changes no reported quantity.
        flags = FLG_HELCTR | FLG_NONUT | FLG_TRUEPOS | FLG_SPEED
        return _osculating_from_calc_state(t, ipl, flags, _GM_SUN, orbel_aa=orbel_aa)
    flags = FLG_HELCTR | FLG_J2000 | FLG_TRUEPOS | FLG_SPEED
    return _osculating_from_calc_state(t, ipl, flags, _GM_SUN, orbel_aa=orbel_aa)


def _calc_orbital_elements(t, ipl: int, iflag: int) -> Tuple[float, ...]:
    """
    Internal function to calculate orbital elements.

    Computes osculating Keplerian elements from the body's current position
    and velocity vectors. Uses state vectors from JPL ephemeris.

    Args:
        t: Skyfield Time object
        ipl: Planet ID
        iflag: Calculation flags

    Returns:
        Flat tuple of 50 floats with orbital elements (padded with 0.0)
    """
    planets = _get_computation_ephemeris()
    zero_elements = tuple([0.0] * 50)

    # Sun and Earth don't have heliocentric orbital elements
    if ipl == SUN:
        return zero_elements

    # Minor bodies (curated asteroids/centaurs 15-20 and numbered asteroids
    # AST_OFFSET+n): compute real osculating elements from the geometric
    # heliocentric state served by the position pipeline. Same membership
    # test as _calc_nod_aps so the two paths stay in lockstep, and never a
    # silent zero-element tuple, which would read as a real orbit at 0 AU.
    if ipl in _MINOR_BODY_NODAPS or (AST_OFFSET <= ipl < FIXSTAR_OFFSET):
        # AST_OFFSET itself ("asteroid 0") is included so the position
        # pipeline's typed UnknownBodyError propagates (the reference raises
        # for it too) instead of falling through to the zero-tuple below.
        return _calc_orbital_elements_minor(t, ipl, iflag)

    # Fictitious/hypothetical bodies (40-58): osculating elements from the
    # library's own runtime models (heliocentric for most, geocentric for the
    # White Moon and Waldemath), never a silent zero-element tuple. Unmodeled
    # ids (e.g. Nibiru 49) let the pipeline's typed error propagate.
    if CUPIDO <= ipl <= WALDEMATH:
        return _calc_orbital_elements_fictitious(
            t, ipl, orbel_aa=bool(iflag & FLG_ORBEL_AA)
        )

    # Get target and center bodies. Every public caller has already passed
    # _validate_orbital_object, so an id missing from the planet map here is
    # an internal-routing bug: raise the same typed error instead of leaking
    # a zero-initialized tuple that would read as a real orbit at 0 AU.
    if ipl not in _PLANET_MAP:
        from .exceptions import Error

        raise Error(f"error in get_orbital_elements(): object {ipl} not valid")

    target_name = _PLANET_MAP[ipl]
    # Represent Earth's heliocentric orbit by the Earth-Moon system barycentre,
    # separating the system's solar orbit from the geocenter's monthly motion.
    if ipl == EARTH:
        target = planets["earth barycenter"]
    else:
        # Try planet center first, fall back to barycenter if not available
        try:
            target = planets[target_name]
        except KeyError:
            if target_name in _PLANET_FALLBACK:
                target = planets[_PLANET_FALLBACK[target_name]]
            else:
                raise

    # For Moon, use geocentric orbit (around Earth)
    if ipl == MOON:
        center = planets["earth"]
        # Two-body mu for the Moon's geocentric osculating orbit must use
        # GM(Earth+Moon): Sun/(Earth+Moon) mass ratio = 328900.5596
        # (Earth alone, 332946, biases the semi-major axis ~+1.2%).
        GM = 0.01720209895**2 / 328900.5596
    else:
        center = planets["sun"]
        if (iflag & FLG_ORBEL_AA) and ipl in _ORBEL_AA_CUM_THROUGH_SELF:
            # Astronomical Almanac convention (FLG_ORBEL_AA): the central mass is
            # the Sun plus every planetary system at or interior to this planet's
            # orbit (Earth as the Earth+Moon barycentre). This shrinks the
            # osculating a and shifts the perihelion-referred angles solidly
            # while leaving the mean longitude invariant; the period family below
            # keeps GM_slots = GM_sun, so it follows only through the changed a.
            # In this context bit 32768 always means ORBEL_AA, never TOPOCTR.
            GM = _GM_SUN * (1.0 + _ORBEL_AA_CUM_THROUGH_SELF[ipl])
        else:
            # Default two-body mu = G(M_sun + M_planet_system). The planet mass
            # cannot be neglected: it shifts a/Q by up to ~1e-5 AU for the giants
            # (mirrors _calc_nod_aps, same _SUN_MASS_RATIOS table).
            mass_ratio = _SUN_MASS_RATIOS.get(ipl)
            GM = _GM_SUN * (1.0 + 1.0 / mass_ratio) if mass_ratio else _GM_SUN

    # Get heliocentric (or geocentric for Moon) position and velocity. In
    # sealed leb mode an out-of-range request evaluates the heliocentric center
    # (the Sun, body 0) first, so the raw EphemerisRangeError would name
    # "Body 0" rather than the planet the caller asked for. Re-raise naming the
    # requested body (mirrors the minor-body escalation, which already does so).
    from .exceptions import EphemerisRangeError

    try:
        center_pos = center.at(t)
        target_pos = target.at(t)
    except EphemerisRangeError as exc:
        if getattr(exc, "body_id", None) == ipl:
            raise
        _start = getattr(exc, "start_jd", None)
        _end = getattr(exc, "end_jd", None)
        _jd = getattr(exc, "requested_jd", float(t.tt))
        if _start is None or _end is None:
            raise
        raise EphemerisRangeError(
            message=(
                f"Body {ipl} at JD {_jd:.6f} is outside active LEB coverage "
                f"range [{_start:.6f}, {_end:.6f}]."
            ),
            requested_jd=_jd,
            start_jd=_start,
            end_jd=_end,
            body_id=ipl,
        ) from exc

    # Heliocentric state (target - Sun); Moon uses the geocentric state above.
    r_helio = target_pos.position.au - center_pos.position.au
    v_helio = target_pos.velocity.au_per_d - center_pos.velocity.au_per_d

    # Barycentric orbital elements (FLG_BARYCTR): the reference re-references the
    # two-body fit from the Sun to the SSB for the trans-jovian planets only
    # (Saturn dM0 ~ -2.89 deg, Uranus -1.54, Neptune -0.40, Pluto -0.31; Mercury
    # through Jupiter show no shift). The Moon is unaffected. Same GM as helio.
    _bary = (
        bool(iflag & FLG_BARYCTR)
        and ipl != MOON
        and _osculating_semi_major_axis(
            float(r_helio[0]),
            float(r_helio[1]),
            float(r_helio[2]),
            float(v_helio[0]),
            float(v_helio[1]),
            float(v_helio[2]),
            GM,
        )
        > _TRANS_JOVIAN_A_AU
    )
    if _bary:
        # Target relative to the SSB (Skyfield positions are BCRS/SSB-relative).
        r_icrs = target_pos.position.au
        v_icrs = target_pos.velocity.au_per_d
    else:
        r_icrs = r_helio
        v_icrs = v_helio

    # Convert from ICRS to the mean J2000 ecliptic using the IAU 2006 frame
    # bias and mean obliquity.
    from .fast_calc import _rotate_icrs_to_ecliptic_j2000

    x, y, z = _rotate_icrs_to_ecliptic_j2000(
        float(r_icrs[0]), float(r_icrs[1]), float(r_icrs[2])
    )
    vx, vy, vz = _rotate_icrs_to_ecliptic_j2000(
        float(v_icrs[0]), float(v_icrs[1]), float(v_icrs[2])
    )

    # Convert the ecliptic heliocentric (Moon: geocentric) state to the
    # 50-slot osculating-element tuple via the shared reducer. The derived
    # period slots use the central mass alone (measured reference
    # convention): the pure Gaussian solar GM for heliocentric orbits, and
    # for the Moon the Earth-only GM from the IAU Sun/Earth mass ratio
    # 332946.0487 (Astronomical Almanac / Explanatory Supplement, 3rd ed.).
    # The external implementation's geocentric GM realization for the Moon
    # is not published; the Earth-only value leaves a ~8e-4 relative offset
    # on the Moon's period family (down from 5.3e-3 with GM(Earth+Moon)).
    if ipl == MOON:
        GM_slots = 0.01720209895**2 / 332946.0487
    else:
        GM_slots = _GM_SUN
    return _orbital_elements_from_ecliptic_state(
        x, y, z, vx, vy, vz, GM, ipl, float(t.tt), GM_slots=GM_slots
    )


def _orbital_elements_from_ecliptic_state(
    x: float,
    y: float,
    z: float,
    vx: float,
    vy: float,
    vz: float,
    GM: float,
    ipl: int,
    tt_jd: float,
    GM_slots: float | None = None,
) -> Tuple[float, ...]:
    """Reduce an ecliptic state vector to the 50-slot orbital-element tuple.

    Shared by the planet/Moon path (heliocentric or, for the Moon, geocentric
    ICRS state rotated to the J2000 ecliptic) and the minor-body path
    (geometric heliocentric J2000 ecliptic state from the position pipeline),
    so both produce identical state->elements math.

    Args:
        x, y, z: Ecliptic position components (AU).
        vx, vy, vz: Ecliptic velocity components (AU/day).
        GM: Gravitational parameter of the central body (AU^3/day^2).
        ipl: Body ID (used only for the Earth special-case in the synodic slot).
        tt_jd: Epoch (JD TT) for the time-of-perihelion-passage slot.
        GM_slots: Optional central-mass-only gravitational parameter for the
            derived period/rate slots [10]-[14]. The osculating GEOMETRY
            (a, e, angles, q, Q) always uses the physical two-body ``GM``;
            the measured reference convention derives the period family from
            the central mass alone, so callers pass the pure Gaussian solar
            GM here (and the Earth-only GM for the Moon). ``None`` keeps
            ``GM`` for every slot.

    Returns:
        Flat tuple of 50 floats with orbital elements (padded with 0.0).
    """
    # Calculate orbital elements from state vectors
    r_mag = math.sqrt(x**2 + y**2 + z**2)
    v_mag = math.sqrt(vx**2 + vy**2 + vz**2)

    # Radial velocity
    r_dot_v = x * vx + y * vy + z * vz
    r_dot_v / r_mag

    # Specific angular momentum vector h = r x v
    hx = y * vz - z * vy
    hy = z * vx - x * vz
    hz = x * vy - y * vx
    h_mag = math.sqrt(hx**2 + hy**2 + hz**2)

    # Node vector n = k x h (k is unit vector in z-direction)
    nx = -hy
    ny = hx
    n_mag = math.sqrt(nx**2 + ny**2)

    # Eccentricity vector e = (v x h)/GM - r/|r|
    e_coef = v_mag**2 / GM - 1.0 / r_mag
    rv_coef = r_dot_v / GM
    ex = e_coef * x - rv_coef * vx
    ey = e_coef * y - rv_coef * vy
    ez = e_coef * z - rv_coef * vz
    e_mag = math.sqrt(ex**2 + ey**2 + ez**2)

    # Specific orbital energy
    epsilon = v_mag**2 / 2 - GM / r_mag

    # Semi-major axis
    if abs(epsilon) > 1e-10:
        a = -GM / (2 * epsilon)
    else:
        a = float("inf")  # Parabolic orbit

    # Handle near-circular and near-equatorial orbits
    e = e_mag
    if e < 1e-10:
        e = 0.0

    # Inclination
    if h_mag > 1e-10:
        i = math.acos(max(-1.0, min(1.0, hz / h_mag)))
    else:
        i = 0.0

    # Longitude of ascending node (Omega)
    if n_mag > 1e-10:
        Omega = math.atan2(ny, nx)
        if Omega < 0:
            Omega += 2 * math.pi
    else:
        Omega = 0.0

    # Argument of perihelion (omega)
    if n_mag > 1e-10 and e_mag > 1e-10:
        n_dot_e = nx * ex + ny * ey
        cos_omega = n_dot_e / (n_mag * e_mag)
        cos_omega = max(-1.0, min(1.0, cos_omega))
        omega = math.acos(cos_omega)
        if ez < 0:
            omega = 2 * math.pi - omega
    elif e_mag > 1e-10:
        # Near-equatorial orbit: measure from x-axis
        omega = math.atan2(ey, ex)
        if omega < 0:
            omega += 2 * math.pi
    else:
        omega = 0.0

    # True anomaly (nu)
    if e_mag > 1e-10:
        e_dot_r = ex * x + ey * y + ez * z
        cos_nu = e_dot_r / (e_mag * r_mag)
        cos_nu = max(-1.0, min(1.0, cos_nu))
        nu = math.acos(cos_nu)
        if r_dot_v < 0:
            nu = 2 * math.pi - nu
    else:
        # Circular orbit: measure from ascending node or x-axis
        if n_mag > 1e-10:
            n_dot_r = nx * x + ny * y
            cos_nu = n_dot_r / (n_mag * r_mag)
            cos_nu = max(-1.0, min(1.0, cos_nu))
            nu = math.acos(cos_nu)
            if z < 0:
                nu = 2 * math.pi - nu
        else:
            nu = math.atan2(y, x)
            if nu < 0:
                nu += 2 * math.pi

    # Eccentric anomaly (E) from true anomaly
    if e < 1.0 and e > 1e-10:
        tan_half_nu = math.tan(nu / 2)
        sqrt_term = math.sqrt((1 - e) / (1 + e))
        E = 2 * math.atan(sqrt_term * tan_half_nu)
        if E < 0:
            E += 2 * math.pi
    else:
        E = nu  # For circular orbits, E = nu

    # Mean anomaly (M) from eccentric anomaly (Kepler's equation)
    if e < 1.0:
        M = E - e * math.sin(E)
        if M < 0:
            M += 2 * math.pi
    else:
        M = E  # For circular or hyperbolic orbits

    # Mean longitude (L)
    L = (Omega + omega + M) % (2 * math.pi)

    # Longitude of perihelion (varpi)
    varpi = (Omega + omega) % (2 * math.pi)

    # Mean daily motion (n) in radians/day. The period/rate slots use the
    # central-mass-only GM when the caller supplies it (measured reference
    # convention); the geometry above already used the physical two-body GM.
    GM_for_slots = GM_slots if GM_slots is not None else GM
    if a > 0 and a < float("inf"):
        n = math.sqrt(GM_for_slots / (a**3))  # radians/day
    else:
        n = 0.0

    # Perihelion and aphelion distances
    if a > 0 and a < float("inf"):
        q = a * (1 - e)  # Perihelion
        Q = a * (1 + e)  # Aphelion
    else:
        q = r_mag  # For parabolic/hyperbolic
        Q = float("inf")

    # Sidereal orbital period (Keplerian period relative to stars)
    # Expressed in tropical years (365.24219 days/year)
    if n > 0:
        P_days = 2 * math.pi / n  # n is in radians/day
        P_sid_years = P_days / 365.24219
    else:
        P_days = float("inf")
        P_sid_years = float("inf")

    # Mean daily motion in degrees/day
    n_deg = math.degrees(n)

    # Tropical orbital period (return to same ecliptic longitude)
    # Accounts for general precession in longitude (~50.29"/year)
    _PRECESSION_DEG_PER_DAY = 50.2882 / 3600.0 / 365.25
    # Tropical period: the time to return to the same equinox-referred
    # longitude, i.e. 360 deg divided by the mean motion augmented by the
    # IAU general precession in longitude (~50.29"/yr) — the standard
    # equinox-referred period (Explanatory Supplement to the Astronomical
    # Almanac, 3rd ed., glossary, "tropical") — expressed in tropical
    # years. External implementations carry a constant extra
    # sidereal/tropical-year factor (~3.9e-5 relative) on this slot whose
    # unit convention has no published basis; that factor is intentionally
    # NOT reproduced (see docs/comparison/known-differences.md).
    if n_deg > 0:
        P_trop_days = 360.0 / (n_deg + _PRECESSION_DEG_PER_DAY)
        P_trop_years = P_trop_days / 365.24219
    else:
        P_trop_years = float("inf")

    # Time of perihelion passage (T)
    # T = t - M/n where M is in radians and n is radians/day
    if n > 0:
        T_jd = tt_jd - M / n
    else:
        T_jd = 0.0

    # Synodic period (relative to Earth's orbital motion)
    # P_syn = P_earth * P_planet / (P_planet - P_earth)
    # Negative for inner planets and Moon (P_planet < P_earth), positive for outer
    # Synodic period from the standard relation 1/P_syn = 1/P_E - 1/P
    # between SIDEREAL periods (textbook definition; e.g. Murray &
    # Dermott, "Solar System Dynamics", 1999): Earth's period in that
    # relation is therefore the sidereal year (365.256363 days, IAU/DE
    # convention).
    P_earth_days = 365.256363
    if P_days > 0 and P_days < float("inf") and ipl != EARTH:
        denom = P_days - P_earth_days
        if abs(denom) > 1e-10:
            P_syn = P_earth_days * P_days / denom
        else:
            P_syn = float("inf")
    else:
        P_syn = 0.0

    # Convert angles to degrees
    i_deg = math.degrees(i)
    Omega_deg = math.degrees(Omega)
    omega_deg = math.degrees(omega)
    nu_deg = math.degrees(nu)
    E_deg = math.degrees(E)
    M_deg = math.degrees(M)
    L_deg = math.degrees(L)
    varpi_deg = math.degrees(varpi)

    # Build the 17-element tuple matching the reference index layout
    elements = (
        a,  # [0] Semi-major axis (AU)
        e,  # [1] Eccentricity
        i_deg,  # [2] Inclination (degrees)
        Omega_deg,  # [3] Longitude of ascending node (degrees)
        omega_deg,  # [4] Argument of perihelion (degrees)
        varpi_deg,  # [5] Longitude of periapsis (degrees)
        M_deg,  # [6] Mean anomaly at epoch (degrees)
        nu_deg,  # [7] True anomaly at epoch (degrees)
        E_deg,  # [8] Eccentric anomaly at epoch (degrees)
        L_deg,  # [9] Mean longitude at epoch (degrees)
        P_sid_years,  # [10] Sidereal orbital period (tropical years)
        n_deg,  # [11] Mean daily motion (degrees/day)
        P_trop_years,  # [12] Tropical period (years)
        P_syn,  # [13] Synodic period (days, negative for inner/Moon)
        T_jd,  # [14] Time of perihelion passage (JD)
        q,  # [15] Perihelion distance (AU)
        Q,  # [16] Aphelion distance (AU)
    )

    # Pad to 50 elements for reference-API compatibility
    elements_50 = elements + tuple([0.0] * (50 - len(elements)))
    return elements_50


def orbit_max_min_true_distance(
    tjdet: float, planet: int, flags: int
) -> Tuple[float, float, float]:
    """
    Calculate the maximum, minimum, and current true geocentric distances.

    Reference API compatible function matching the reference
    ``orbit_max_min_true_distance``.

    This function computes the maximum and minimum true distances from Earth
    that a planet can reach during its orbital motion, plus the current true
    distance at the given time.

    For outer planets (Mars-Pluto), the minimum distance occurs around opposition
    and the maximum around conjunction with the Sun.

    For inner planets (Mercury, Venus), the minimum distance occurs near inferior
    conjunction and the maximum near superior conjunction.

    Args:
        tjdet: Julian Day in Terrestrial Time (TT/ET) - used to determine current
                orbital elements and current true distance.
        planet: Planet/body ID (SUN, MOON, etc.)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple of (max_distance, min_distance, true_distance) in AU.
        Order matches the reference API: (max, min, true).

    Example:
        >>> from libephemeris import orbit_max_min_true_distance, MARS
        >>> max_dist, min_dist, true_dist = orbit_max_min_true_distance(
        ...     2451545.0, MARS, 0
        ... )
        >>> print(f"Mars: {min_dist:.4f} - {max_dist:.4f} AU (now {true_dist:.4f})")

    Note:
        - For geocentric calculations, these represent the range of Earth-planet
          distances possible during the synodic cycle.
        - The Moon's distance is calculated as geocentric (Earth-Moon distance).
        - Sun returns the Earth-orbit range (aphelion/perihelion, ~1.017/0.983
          AU): the reference ignores FLG_HELCTR for the Sun and Moon (hel == geo)
          and reports the terrestrial orbit for the Sun.
        - The lunar nodes/apogees (10-13) and interpolated apsides (21, 22)
          raise ``object N not valid``; ids 23-39 raise ``illegal planet
          number``; negatives raise ``object N not valid``.
        - The White Moon (56) and Waldemath (58) are geocentric circular
          constructions here, so max == min == true == the model radius; this
          intentionally diverges from the reference's eccentric element set.
    """
    _validate_orbital_object("orbit_max_min_true_distance", planet, sun_valid=True)
    ts = get_timescale()
    t = ts.tt_jd(tjdet)
    return _calc_orbit_max_min_true_distance(t, planet, flags, tjdet)


def _emb_osculating_ecliptic_elements(t, jd_tt: float) -> Tuple[float, ...]:
    """Osculating heliocentric ecliptic elements of the Earth-Moon barycenter.

    For a geocentric extreme-distance query the observer is the Earth, whose
    osculating orbit about the Sun is (to the precision this function targets)
    the EMB ellipse. So the EMB provides the second ellipse for every
    geocentric two-body distance extremum, and its own perihelion/aphelion for
    the Sun and Earth special cases.

    Args:
        t: Skyfield Time object at the osculating epoch.
        jd_tt: Epoch (JD TT); only feeds the (unused) perihelion-passage slot.

    Returns:
        The shared 50-slot osculating-element tuple for the EMB.
    """
    planets = _get_computation_ephemeris()
    sun = planets["sun"]
    emb = planets["earth barycenter"]

    sun_pos = sun.at(t)
    emb_pos = emb.at(t)
    r_icrs = emb_pos.position.au - sun_pos.position.au
    v_icrs = emb_pos.velocity.au_per_d - sun_pos.velocity.au_per_d

    from .fast_calc import _rotate_icrs_to_ecliptic_j2000

    x, y, z = _rotate_icrs_to_ecliptic_j2000(
        float(r_icrs[0]), float(r_icrs[1]), float(r_icrs[2])
    )
    vx, vy, vz = _rotate_icrs_to_ecliptic_j2000(
        float(v_icrs[0]), float(v_icrs[1]), float(v_icrs[2])
    )

    GM = 0.01720209895**2
    return _orbital_elements_from_ecliptic_state(x, y, z, vx, vy, vz, GM, EARTH, jd_tt)


def _ellipse_positions(elements: Tuple[float, float, float, float, float], anom):
    """Ecliptic Cartesian positions (AU) sampled along an osculating ellipse.

    Args:
        elements: ``(a, e, inc_rad, node_rad, argperi_rad)``.
        anom: numpy array of eccentric anomalies (radians).

    Returns:
        ``(N, 3)`` numpy array of ecliptic positions, one per eccentric anomaly.
    """
    import numpy as np

    a, e, inc, node, argp = elements
    cos_e = np.cos(anom)
    sin_e = np.sin(anom)
    # Perifocal plane (focus at the Sun).
    x_pf = a * (cos_e - e)
    y_pf = a * math.sqrt(max(0.0, 1.0 - e * e)) * sin_e
    # Rotate by argument of perihelion about the orbit normal.
    cos_w = math.cos(argp)
    sin_w = math.sin(argp)
    x1 = x_pf * cos_w - y_pf * sin_w
    y1 = x_pf * sin_w + y_pf * cos_w
    # Incline about the node line, then rotate by the ascending node.
    cos_i = math.cos(inc)
    sin_i = math.sin(inc)
    y2 = y1 * cos_i
    z2 = y1 * sin_i
    cos_o = math.cos(node)
    sin_o = math.sin(node)
    x = x1 * cos_o - y2 * sin_o
    y = x1 * sin_o + y2 * cos_o
    return np.stack([x, y, z2], axis=-1)


def _two_ellipse_min_max(
    el1: Tuple[float, float, float, float, float],
    el2: Tuple[float, float, float, float, float],
    n_coarse: int = 1440,
    refine_iters: int = 25,
) -> Tuple[float, float]:
    """True 3D minimum and maximum distance between two osculating ellipses.

    Independent numerical search of ``|P1(E1) - P2(E2)|`` over both eccentric
    anomalies: a coarse anomaly grid locates the global extremum cells, then a
    local shrinking-window grid refines each to well under 1e-8 AU. The result
    is a geometric property of the two ellipses, independent of the search
    method.

    Args:
        el1: ``(a, e, inc_rad, node_rad, argperi_rad)`` of the first ellipse.
        el2: Elements of the second ellipse (same layout).
        n_coarse: Samples per anomaly on the coarse grid.
        refine_iters: Local-refinement zoom iterations.

    Returns:
        ``(max_distance, min_distance)`` in AU.
    """
    import numpy as np

    two_pi = 2.0 * math.pi
    grid = np.linspace(0.0, two_pi, n_coarse, endpoint=False)
    pa = _ellipse_positions(el1, grid)
    pb = _ellipse_positions(el2, grid)
    a_sq = np.sum(pa * pa, axis=1)
    b_sq = np.sum(pb * pb, axis=1)
    d_sq = a_sq[:, None] + b_sq[None, :] - 2.0 * (pa @ pb.T)
    np.maximum(d_sq, 0.0, out=d_sq)
    i_min = np.unravel_index(int(np.argmin(d_sq)), d_sq.shape)
    i_max = np.unravel_index(int(np.argmax(d_sq)), d_sq.shape)

    def _refine(c1: float, c2: float, want_max: bool) -> float:
        half = two_pi / n_coarse
        best = 0.0
        for _ in range(refine_iters):
            g1 = np.linspace(c1 - half, c1 + half, 9)
            g2 = np.linspace(c2 - half, c2 + half, 9)
            qa = _ellipse_positions(el1, g1)
            qb = _ellipse_positions(el2, g2)
            qa_sq = np.sum(qa * qa, axis=1)
            qb_sq = np.sum(qb * qb, axis=1)
            dd = qa_sq[:, None] + qb_sq[None, :] - 2.0 * (qa @ qb.T)
            np.maximum(dd, 0.0, out=dd)
            j = np.unravel_index(
                int(np.argmax(dd) if want_max else np.argmin(dd)), dd.shape
            )
            c1 = float(g1[j[0]])
            c2 = float(g2[j[1]])
            best = float(dd[j])
            half *= 0.25
        return math.sqrt(best)

    d_max = _refine(float(grid[i_max[0]]), float(grid[i_max[1]]), True)
    d_min = _refine(float(grid[i_min[0]]), float(grid[i_min[1]]), False)
    return d_max, d_min


def _calc_orbit_max_min_true_distance(
    t, ipl: int, iflag: int, tjd_et: float = 0.0
) -> Tuple[float, float, float]:
    """
    Internal function to calculate max/min/true geocentric distances.

    The extremes are the true 3D minimum and maximum distance between two
    osculating heliocentric ellipses at the epoch ``t``:

        * geocentric planet/asteroid: the body's ellipse vs the Earth-Moon
          barycenter's ellipse, solved numerically (see
          :func:`_two_ellipse_min_max`);
        * geocentric Sun: distance from the Sun (a point) to the EMB ellipse,
          i.e. the EMB perihelion/aphelion;
        * geocentric Earth: the observer against its own orbit, so minimum 0
          (coincident) and maximum the orbit diameter (2a);
        * geocentric Moon: the geocentric lunar osculating perigee/apogee;
        * heliocentric (``FLG_HELCTR``): the body's own perihelion/aphelion.

    Args:
        t: Skyfield Time object at the osculating epoch.
        ipl: Body ID.
        iflag: Calculation flags (``FLG_HELCTR`` selects the heliocentric case).
        tjd_et: Julian Day TT/ET for the current true distance.

    Returns:
        Tuple of (max_distance, min_distance, true_distance) in AU.
    """
    # The reference ignores FLG_HELCTR for the Sun and the Moon (hel == geo):
    # the Sun reports the Earth-orbit range and the Moon stays geocentric.
    if ipl in (SUN, MOON):
        iflag &= ~FLG_HELCTR

    # FLG_BARYCTR re-references the orbit to the SSB only for bodies beyond
    # Jupiter; for every interior body the reference returns the heliocentric
    # result for all three channels (max, min and true distance). Normalize the
    # interior case to FLG_HELCTR so the true-distance calc() below and the
    # element branch both stay heliocentric and agree with the reference.
    if (iflag & FLG_BARYCTR) and ipl not in (SUN, MOON):
        try:
            _a_helio = _calc_orbital_elements(t, ipl, FLG_HELCTR)[0]
        except (EphemerisRangeError, UnknownBodyError):
            # A typed source error is the caller's answer, not a reason to
            # silently classify the body as interior to Jupiter and downgrade
            # an explicit FLG_BARYCTR request: the calculation below fails on
            # the same cause anyway, with the diagnosis lost.
            raise
        except (ValueError, ArithmeticError, KeyError):
            # Genuinely unresolvable geometry: fall back to the heliocentric
            # convention, the documented default for a non-trans-jovian body.
            _a_helio = 0.0
        if not (_a_helio > _TRANS_JOVIAN_A_AU):
            iflag = (iflag & ~FLG_BARYCTR) | FLG_HELCTR

    # Current true distance at the requested ephemeris time (ET/TT input).
    # "True" here is the GEOMETRIC distance (no light-time correction):
    # measured reference behavior returns the FLG_TRUEPOS distance in this
    # slot, and the function's own name promises the true distance. Without
    # the flag the slot silently carried the apparent, light-time-corrected
    # value (up to ~1.6e-4 AU short for Mercury).
    true_dist = 0.0
    if tjd_et > 0:
        try:
            pos, _ = calc(tjd_et, ipl, iflag | FLG_TRUEPOS)
            true_dist = float(pos[2])
        except (IndexError, TypeError, ValueError):
            pass

    def _ellipse_from_elements(
        el: Tuple[float, ...],
    ) -> Tuple[float, float, float, float, float]:
        return (
            el[0],
            el[1],
            math.radians(el[2]),
            math.radians(el[3]),
            math.radians(el[4]),
        )

    # Heliocentric / barycentric: extremes are simply the body's own
    # perihelion/aphelion. The reference honors FLG_BARYCTR here too, so a
    # trans-jovian body reports its barycentric aphelion/perihelion (measured).
    if iflag & (FLG_HELCTR | FLG_BARYCTR):
        el = _calc_orbital_elements(t, ipl, iflag)
        if el[0] <= 0.0:
            return (0.0, 0.0, true_dist)
        return (el[16], el[15], true_dist)

    # Geocentric Sun: distance from the Sun to the EMB ellipse.
    if ipl == SUN:
        emb = _emb_osculating_ecliptic_elements(t, tjd_et)
        return (emb[16], emb[15], true_dist)

    # Earth is the observer: coincident minimum (0) and orbit-diameter maximum.
    if ipl == EARTH:
        emb = _emb_osculating_ecliptic_elements(t, tjd_et)
        return (2.0 * emb[0], 0.0, true_dist)

    # Geocentric Moon: the geocentric lunar osculating perigee/apogee.
    if ipl == MOON:
        el = _calc_orbital_elements(t, ipl, iflag)
        if el[0] <= 0.0:
            return (0.0, 0.0, true_dist)
        return (el[16], el[15], true_dist)

    # Geocentric symbolic points (White Moon 56, Waldemath 58): these orbit the
    # Earth, so the extremes are their own geocentric apogee/perigee (Earth
    # apogee/perigee), not a two-ellipse solve against the EMB. For the
    # near-circular published constructions this yields max == min == radius.
    if ipl in _FICT_GEOCENTRIC_IDS:
        el = _calc_orbital_elements(t, ipl, iflag)
        if el[0] <= 0.0:
            return (0.0, 0.0, true_dist)
        return (el[16], el[15], true_dist)

    # Geocentric planet/asteroid: true 3D extremum between the body's ellipse
    # and the EMB's.
    el = _calc_orbital_elements(t, ipl, iflag)
    if el[0] <= 0.0 or el[1] >= 1.0:  # invalid or non-elliptical osculating orbit
        return (0.0, 0.0, true_dist)
    emb = _emb_osculating_ecliptic_elements(t, tjd_et)
    max_dist, min_dist = _two_ellipse_min_max(
        _ellipse_from_elements(el), _ellipse_from_elements(emb)
    )
    return (max_dist, min_dist, true_dist)


def _calc_nod_aps_osculating(
    t, ipl: int, iflag: int
) -> Tuple[PosTuple, PosTuple, PosTuple, PosTuple]:
    """
    Calculate orbital nodes and apsides using osculating (instantaneous) elements.

    Uses the planet's current position and velocity to compute orbital elements.
    These are the "true" instantaneous orbital parameters that include short-term
    perturbations from other planets.
    """
    planets = _get_computation_ephemeris()
    zero_pos: PosTuple = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    target_name = _PLANET_MAP.get(ipl)
    if not target_name:
        return (zero_pos, zero_pos, zero_pos, zero_pos)

    # Try planet center first, fall back to barycenter if not available
    try:
        target = planets[target_name]
    except KeyError:
        if target_name in _PLANET_FALLBACK:
            target = planets[_PLANET_FALLBACK[target_name]]
        else:
            raise
    sun = planets["sun"]

    # Get heliocentric position and velocity
    sun_pos = sun.at(t)
    target_pos = target.at(t)

    # Heliocentric vectors in AU and AU/day (ICRS frame)
    r_icrs = target_pos.position.au - sun_pos.position.au
    v_icrs = target_pos.velocity.au_per_d - sun_pos.velocity.au_per_d

    # Convert from ICRS to the mean J2000 ecliptic using the IAU 2006 frame
    # bias and mean obliquity.
    from .fast_calc import _rotate_icrs_to_ecliptic_j2000

    x_ecl, y_ecl, z_ecl = _rotate_icrs_to_ecliptic_j2000(
        float(r_icrs[0]), float(r_icrs[1]), float(r_icrs[2])
    )
    vx_ecl, vy_ecl, vz_ecl = _rotate_icrs_to_ecliptic_j2000(
        float(v_icrs[0]), float(v_icrs[1]), float(v_icrs[2])
    )

    # GM_sun in AU^3/day^2
    GM_sun = 0.01720209895**2

    r_mag = math.sqrt(x_ecl**2 + y_ecl**2 + z_ecl**2)
    v_mag = math.sqrt(vx_ecl**2 + vy_ecl**2 + vz_ecl**2)

    # Angular momentum
    hx = y_ecl * vz_ecl - z_ecl * vy_ecl
    hy = z_ecl * vx_ecl - x_ecl * vz_ecl
    hz = x_ecl * vy_ecl - y_ecl * vx_ecl
    h_mag = math.sqrt(hx**2 + hy**2 + hz**2)

    # Node vector
    nx = -hy
    ny = hx
    n_mag = math.sqrt(nx**2 + ny**2)

    # Eccentricity vector
    r_dot_v = x_ecl * vx_ecl + y_ecl * vy_ecl + z_ecl * vz_ecl
    coef1 = v_mag**2 / GM_sun - 1.0 / r_mag
    coef2 = r_dot_v / GM_sun

    ex = coef1 * x_ecl - coef2 * vx_ecl
    ey = coef1 * y_ecl - coef2 * vy_ecl
    ez = coef1 * z_ecl - coef2 * vz_ecl
    e_mag = math.sqrt(ex**2 + ey**2 + ez**2)

    # Semi-major axis
    a = 1.0 / (2.0 / r_mag - v_mag**2 / GM_sun)

    # Inclination
    incl = math.acos(hz / h_mag) if h_mag > 0 else 0.0

    # Longitude of ascending node
    if n_mag > 1e-10:
        Omega = math.atan2(ny, nx)
        if Omega < 0:
            Omega += 2 * math.pi
    else:
        Omega = 0.0

    # Argument of perihelion
    if n_mag > 1e-10 and e_mag > 1e-10:
        n_dot_e = nx * ex + ny * ey
        cos_omega = n_dot_e / (n_mag * e_mag)
        cos_omega = max(-1.0, min(1.0, cos_omega))
        omega = math.acos(cos_omega)
        if ez < 0:
            omega = 2 * math.pi - omega
    else:
        omega = 0.0

    # Longitude of perihelion
    varpi = Omega + omega

    # Calculate positions
    Omega_deg = math.degrees(Omega) % 360.0
    varpi_deg = math.degrees(varpi) % 360.0

    # Distances
    if e_mag < 1.0:
        p = a * (1 - e_mag**2)
        r_asc = (
            p / (1 + e_mag * math.cos(-omega))
            if abs(1 + e_mag * math.cos(-omega)) > 1e-10
            else a
        )
        r_dsc = (
            p / (1 + e_mag * math.cos(math.pi - omega))
            if abs(1 + e_mag * math.cos(math.pi - omega)) > 1e-10
            else a
        )
    else:
        r_asc = r_dsc = a

    r_peri = a * (1 - e_mag)
    r_aphe = a * (1 + e_mag)

    # Latitudes at apsides
    lat_peri = (
        math.degrees(math.asin(math.sin(incl) * math.sin(omega))) if incl > 0 else 0.0
    )
    lat_aphe = -lat_peri

    xnasc: PosTuple = (Omega_deg, 0.0, r_asc, 0.0, 0.0, 0.0)
    xndsc: PosTuple = ((Omega_deg + 180.0) % 360.0, 0.0, r_dsc, 0.0, 0.0, 0.0)
    xperi: PosTuple = (varpi_deg, lat_peri, r_peri, 0.0, 0.0, 0.0)
    xaphe: PosTuple = ((varpi_deg + 180.0) % 360.0, lat_aphe, r_aphe, 0.0, 0.0, 0.0)

    return (xnasc, xndsc, xperi, xaphe)


# =============================================================================
# PLANETARY PHENOMENA: Phase, Elongation, Magnitude
# =============================================================================

# Physical radii of celestial bodies in kilometers
# Sources: NASA Planetary Fact Sheet volumetric mean radii
# For gas/ice giants (Jupiter, Saturn, Uranus, Neptune), volumetric mean radii
# are used rather than equatorial radii, as these better represent the
# sphere-equivalent size for apparent diameter calculations. This matches
# the convention used by standard ephemeris implementations.
# Reference: https://nssdc.gsfc.nasa.gov/planetary/factsheet/
_BODY_RADIUS_KM = {
    SUN: 696000.0,  # Solar radius (NASA fact sheet)
    MOON: 1737.5,  # Lunar mean radius (NASA fact sheet)
    MERCURY: 2439.4,  # Mercury volumetric mean radius
    VENUS: 6051.8,  # Venus volumetric mean radius
    MARS: 3389.5,  # Mars volumetric mean radius
    JUPITER: 69911.0,  # Jupiter volumetric mean radius
    SATURN: 58232.0,  # Saturn volumetric mean radius (disk only, excludes rings)
    URANUS: 25362.0,  # Uranus volumetric mean radius
    NEPTUNE: 24622.0,  # Neptune volumetric mean radius
    PLUTO: 1188.3,  # Pluto mean radius
}

# Asteroid photometric data for pheno: IAU H-G system (H absolute magnitude,
# G slope parameter) and diameters from NASA/JPL's Small-Body Database,
# retrieved 2026-07-13:
# https://ssd-api.jpl.nasa.gov/doc/sbdb.html
#
# Chiron and Pholus have no fitted G value in that data, so the conventional
# H-G default G=0.15 is used.  The phase law is from Bowell et al. (1989),
# "Application of photometric models to asteroids", Asteroids II.
_ASTEROID_HG = {
    CHIRON: (5.54, 0.15),
    PHOLUS: (7.07, 0.15),
    CERES: (3.34, 0.12),
    PALLAS: (4.12, 0.11),
    JUNO: (5.19, 0.32),
    VESTA: (3.25, 0.32),
}
_ASTEROID_DIAMETER_KM = {
    CHIRON: 166.0,
    PHOLUS: 190.0,
    CERES: 939.4,
    PALLAS: 513.0,
    JUNO: 246.596,
    VESTA: 522.77,
}

# Numbered minor bodies the library can position (the SPK registry in
# exotic_bodies.py) that also carry published photometry. Each row is
# (H, G, diameter_km): H the absolute magnitude, G the IAU H-G phase slope,
# and the effective diameter. Values are transcribed from NASA/JPL's
# Small-Body Database (https://ssd-api.jpl.nasa.gov/doc/sbdb.html), retrieved
# 2026-07-27. Where SBDB lists no fitted G, the conventional H-G default
# G = 0.15 is used (same convention as the curated Chiron/Pholus rows above).
# Keyed by the runtime ``AST_OFFSET + minor-planet-number`` id — the same id
# the SPK/registry path serves the position under.
#
# Positionable bodies deliberately absent here — e.g. 7066 Nessus and the
# distant TNOs (Eris, Sedna, Haumea, Makemake, Ixion, Orcus, Quaoar,
# Gonggong), for which SBDB publishes no effective diameter — keep the
# existing point-body behavior: pheno's diameter/magnitude slots ([3]/[4])
# report 0.0. Only bodies with BOTH a magnitude and a diameter are listed, so
# _ASTEROID_HG and _ASTEROID_DIAMETER_KM stay in lockstep.
_NUMBERED_ASTEROID_PHOTOMETRY: dict[int, tuple[float, float, float]] = {
    # --- main-belt ---
    AST_OFFSET + 10: (5.65, 0.15, 407.12),  # 10 Hygiea (SBDB: no G)
    AST_OFFSET + 16: (6.20, 0.20, 222.0),  # 16 Psyche
    AST_OFFSET + 52: (6.65, 0.18, 303.918),  # 52 Europa
    AST_OFFSET + 55: (7.82, 0.15, 84.794),  # 55 Pandora (SBDB: no G)
    AST_OFFSET + 80: (8.06, 0.15, 68.563),  # 80 Sappho (SBDB: no G)
    AST_OFFSET + 87: (6.93, 0.15, 253.051),  # 87 Sylvia (SBDB: no G)
    AST_OFFSET + 511: (6.42, 0.16, 270.327),  # 511 Davida
    AST_OFFSET + 704: (6.34, -0.02, 306.313),  # 704 Interamnia
    AST_OFFSET + 1181: (11.36, 0.15, 20.492),  # 1181 Lilith (SBDB: no G)
    # --- centaurs ---
    AST_OFFSET + 944: (10.55, 0.15, 38.0),  # 944 Hidalgo (SBDB: no G)
    AST_OFFSET + 8405: (9.20, 0.15, 66.0),  # 8405 Asbolus (SBDB: no G)
    AST_OFFSET + 10199: (6.55, 0.15, 302.0),  # 10199 Chariklo (SBDB: no G)
    # --- near-Earth asteroids ---
    AST_OFFSET + 433: (10.40, 0.46, 16.84),  # 433 Eros
    AST_OFFSET + 1221: (17.37, 0.15, 1.0),  # 1221 Amor (SBDB: no G)
    AST_OFFSET + 1566: (16.53, 0.15, 1.0),  # 1566 Icarus (SBDB: no G)
    AST_OFFSET + 1685: (14.26, 0.15, 3.4),  # 1685 Toro (SBDB: no G)
    AST_OFFSET + 4179: (15.29, 0.10, 5.4),  # 4179 Toutatis
    AST_OFFSET + 25143: (19.26, 0.15, 0.33),  # 25143 Itokawa (SBDB: no G)
    AST_OFFSET + 99942: (19.09, 0.24, 0.34),  # 99942 Apophis
    AST_OFFSET + 162173: (19.55, 0.15, 0.896),  # 162173 Ryugu (SBDB: no G)
    # --- trans-Neptunian (only those with a published diameter) ---
    AST_OFFSET + 20000: (3.79, 0.15, 900.0),  # 20000 Varuna (SBDB: no G)
}
for _bid, (_h_abs, _g_slope, _diam_km) in _NUMBERED_ASTEROID_PHOTOMETRY.items():
    _ASTEROID_HG[_bid] = (_h_abs, _g_slope)
    _ASTEROID_DIAMETER_KM[_bid] = _diam_km


def _asteroid_hg_magnitude(
    ipl: int, phase_angle_deg: float, helio_dist: float, geo_dist: float
) -> float:
    """Visual magnitude from the IAU H-G phase system.

    mag = H + 5*log10(r*Delta) - 2.5*log10((1-G)*phi1 + G*phi2) with the
    standard phi1/phi2 phase functions (Bowell et al. 1989).
    """
    h_mag, g_slope = _ASTEROID_HG[ipl]
    if helio_dist <= 0 or geo_dist <= 0:
        return h_mag
    tan_half = math.tan(math.radians(phase_angle_deg) / 2.0)
    phi1 = math.exp(-3.33 * tan_half**0.63)
    phi2 = math.exp(-1.87 * tan_half**1.22)
    blend = (1.0 - g_slope) * phi1 + g_slope * phi2
    mag = h_mag + 5.0 * math.log10(helio_dist * geo_dist)
    if blend > 0:
        mag -= 2.5 * math.log10(blend)
    return mag


# 1 AU in kilometers (IAU 2012 definition)
_AU_KM = 149597870.7

# Earth equatorial radius (km), IAU/IERS 2003 (a = 6378136.6 m; Astronomical
# Almanac 2006). Used for the Moon's equatorial horizontal parallax reported in
# pheno slot [5]: sin(pi) = R_eq / d, with d the geometric geocentric distance
# (Explanatory Supplement to the Astronomical Almanac; Meeus, Astronomical
# Algorithms, ch. 40 "Parallax"). This is the same constant already used by the
# lunar distance modulus in _calc_moon_magnitude (6378136.6 m).
_EARTH_EQ_RADIUS_KM = 6378.1366

# Conversion factor: radians to arcseconds
_RAD_TO_ARCSEC = 206264.80624709636  # (180/pi) * 3600


def _calc_apparent_diameter(radius_km: float, distance_au: float) -> float:
    """
    Calculate apparent angular diameter in degrees.

    Uses the small-angle approximation which is accurate for all solar system bodies
    as seen from Earth (maximum angular size is ~0.5 degrees for Sun/Moon).

    The formula is:
        diameter_deg = 2 * radius_km / distance_km * RAD_TO_DEG
                     = 2 * radius_km / (distance_au * AU_KM) * (180 / pi)

    Returns degrees for reference-API compatibility (pheno_ut attr[3]).

    Args:
        radius_km: Physical radius of the body in kilometers
        distance_au: Distance from observer to body in AU

    Returns:
        Apparent angular diameter in degrees
    """
    if distance_au <= 0:
        return 0.0
    distance_km = distance_au * _AU_KM
    # Angular diameter = 2 * arctan(radius/distance) ≈ 2 * radius/distance (small angle)
    # Convert from radians to degrees (not arcseconds) for reference-API compatibility
    return 2.0 * radius_km / distance_km * _RAD_TO_ARCSEC / 3600.0


# Visual magnitude parameters for outer planets (generic fallback:
# V = V(1,0) + 5*log10(r*d) + B1*a + B2*a^2 + B3*a^3).
# Every classical planet now has a dedicated Mallama & Hilton (2018) branch in
# _calc_planet_magnitude (Mercury, Venus, Mars, Jupiter, Saturn, Uranus,
# Neptune, Pluto), so this table is empty and only kept as an extension hook.
_PLANET_MAG_PARAMS: dict[int, tuple[float, float, float, float]] = {}

# --- Uranus photometry: Mallama & Hilton (2018) rotational-pole geometry ------
# Uranus reflects light more strongly when a pole faces the observer (its polar
# regions are depleted in light-absorbing methane), so its apparent magnitude
# depends on the planetographic sub-Earth and sub-solar latitudes. The IAU
# WGCCRE (2009/2015) north-pole direction of Uranus is constant in time:
#   right ascension alpha0 = 257.311 deg, declination delta0 = -15.175 deg
# (equatorial J2000). Flattening f = 0.0022927 (Mallama & Hilton 2018 Eq. 13,
# after Schmude et al. 2015).
_URANUS_POLE_RA = math.radians(257.311)
_URANUS_POLE_DEC = math.radians(-15.175)
_URANUS_FLATTENING = 0.0022927
_URANUS_ONE_MINUS_F_SQ = (1.0 - _URANUS_FLATTENING) ** 2

# J2000 epoch for Saturn ring calculations
_J2000 = 2451545.0


def _calc_moon_magnitude(
    phase_angle: float, geo_dist_au: float, helio_dist_au: float
) -> float:
    """
    Calculate Moon's visual magnitude using a piecewise photometric model.

    For moderate phase angles (α ≤ 147.14°), uses the Allen (1976) formula
    from Astrophysical Quantities with a linear phase coefficient and a
    quartic brightening/dimming term.

    For large phase angles (α > 147.14°), switches to the Samaha et al. (1969)
    cube-phase model which correctly captures the rapid dimming of the thin
    crescent Moon approaching new Moon. The stitch angle (147.14°) is chosen
    so the two formulas produce equal magnitudes at the transition.

    Distance correction uses both geocentric and heliocentric (Sun-Moon)
    distances, converting the geocentric distance to Earth radii for the
    standard 5·log10 distance modulus.

    Reference values:
    - Full Moon at mean distance: V ≈ -12.73 mag
    - Quarter Moon (α≈90°): V ≈ -10.0 to -10.5
    - New Moon (α≈180°): V → +∞ (vanishingly dim)

    Args:
        phase_angle: Sun-Moon-Earth angle in degrees (0° = full, 180° = new)
        geo_dist_au: Earth-Moon distance in AU
        helio_dist_au: Sun-Moon distance in AU

    Returns:
        Visual magnitude (more negative = brighter)

    References:
        - Allen, C.W., 1976, Astrophysical Quantities
        - Samaha, A.E., Asaad, A.S., Mikhail, J.S. (1969), "Visibility of
          the New Moon", Bulletin of Observatory Helwan, 84
    """
    # Constants for distance correction
    # AUNIT and EARTH_RADIUS from IAU/AA standards
    aunit_m = 1.49597870700e11  # 1 AU in meters (DE431)
    earth_radius_m = 6378136.6  # Earth equatorial radius in meters (AA 2006)

    alpha = abs(phase_angle)

    # Distance correction: 5 * log10(d_geo_earthradii * d_helio_au)
    # Converts geocentric distance from AU to Earth radii, then applies
    # the standard distance modulus with heliocentric distance in AU
    if geo_dist_au > 0 and helio_dist_au > 0:
        dist_correction = 5.0 * math.log10(
            geo_dist_au * helio_dist_au * aunit_m / earth_radius_m
        )
    else:
        dist_correction = 0.0

    # Stitch angle where Allen and Samaha formulas produce equal magnitudes
    _stitch_angle = 147.1385465

    if alpha <= _stitch_angle:
        # Allen (1976) formula from Astrophysical Quantities
        # V = -21.62 + 0.026·|α| + 4e-9·α⁴ + dist_correction
        base = -21.62 + 0.026 * alpha + 4.0e-9 * alpha**4
    else:
        # Samaha et al. (1969) cube-phase model for thin crescent
        # V = -4.5444 - 2.5·log10((180-α)³) + dist_correction
        # This properly diverges to +∞ as α→180° (new Moon)
        remainder = 180.0 - alpha
        if remainder > 0:
            base = -4.5444 - 2.5 * math.log10(remainder**3)
        else:
            # At exactly α=180° (geometrically impossible in practice)
            base = 50.0
    return base + dist_correction


def _moon_horizontal_parallax_deg(geo_dist_geometric_au: float) -> float:
    """Equatorial horizontal parallax of the Moon for pheno slot [5].

    The reference reports the Moon's geocentric equatorial horizontal parallax
    in pheno attribute [5] (0.0 for every other body). The convention, verified
    to 0.000000" against the reference over perigee-to-apogee dates, is

        sin(pi) = R_eq / d

    with R_eq the Earth equatorial radius (6378.1366 km, IAU/IERS 2003) and d
    the *geometric* geocentric distance (light-time removed), independent of the
    request's FLG_TRUEPOS bit. See the Explanatory Supplement to the
    Astronomical Almanac and Meeus, Astronomical Algorithms, ch. 40.

    Args:
        geo_dist_geometric_au: Geometric geocentric Moon distance in AU.

    Returns:
        Equatorial horizontal parallax in degrees (0.0 if the distance is
        non-positive).
    """
    d_km = geo_dist_geometric_au * _AU_KM
    if d_km <= 0.0:
        return 0.0
    return math.degrees(math.asin(_EARTH_EQ_RADIUS_KM / d_km))


def _calc_pheno_leb(tjd_ut: float, ipl: int, iflag: int) -> Tuple[float, ...]:
    """Compute planetary phenomena using the LEB fast path.

    Calls ``fast_calc.fast_calc_ut`` directly to obtain positions without
    loading the Skyfield DE kernel.  The results are numerically equivalent
    to ``_calc_pheno`` for all bodies
    supported by the LEB reader.

    Args:
        tjd_ut: Julian Day in Universal Time (UT1).
        ipl: Planet/body ID (SUN, MOON, MERCURY, etc.).
        iflag: Calculation flags.

    Returns:
        Tuple of 20 floats (matching the reference API).

    Raises:
        KeyError: If the body is not available in the LEB reader.
        ValueError: If the date is outside the LEB coverage.
    """
    from . import fast_calc
    from .state import get_leb_reader
    from .utils import angular_separation

    # Built-in asteroids requested via the AST_OFFSET + number alias resolve to
    # their dedicated id (AST_OFFSET + 4 -> Vesta, AST_OFFSET + 5145 -> Pholus),
    # exactly as calc_ut serves the position, so their pheno matches the
    # built-in id. Registry-served numbered asteroids keep their own id.
    ipl = _remap_ast_offset(ipl)

    reader = get_leb_reader()
    if reader is None:
        raise KeyError("No LEB reader available for _calc_pheno_leb")

    def _leb_calc(jd, body, flags):
        """Call fast_calc directly — no Skyfield fallback."""
        return fast_calc.fast_calc_ut(reader, jd, body, flags)

    # Preserve FLG_TRUEPOS to match _calc_pheno() which uses geometric
    # positions when TRUEPOS is set, and FLG_TOPOCTR so the observer's
    # topocentric place drives every distance-derived quantity (the reference
    # honors FLG_TOPOCTR in pheno). NOABERR/NOGDEFL are NOT propagated because
    # their semantics differ between fast_calc and Skyfield.
    base_flags = FLG_SPEED | (iflag & (FLG_TRUEPOS | FLG_TOPOCTR))

    # Unsupported bodies (nodes, apogees, etc.) — match _calc_pheno() behavior
    _PHENO_SUPPORTED = {
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
    }
    if ipl not in _PHENO_SUPPORTED and ipl not in _ASTEROID_HG:
        # Point and model bodies without a physical disc: measured reference
        # behavior fills the geometric phase triplet from their positions.
        return _pheno_positional(tjd_ut, ipl, iflag)

    # ------------------------------------------------------------------
    # Special case: Sun
    # ------------------------------------------------------------------
    if ipl == SUN:
        sun_pos, _ = _leb_calc(tjd_ut, SUN, base_flags)
        sun_dist_au = float(sun_pos[2])

        sun_radius_km = _BODY_RADIUS_KM.get(SUN, 695700.0)
        diameter = _calc_apparent_diameter(sun_radius_km, sun_dist_au)
        magnitude = (
            -26.86 + 5.0 * math.log10(sun_dist_au) if sun_dist_au > 0 else -26.86
        )

        # Phase quantities are inapplicable to the self-luminous Sun: the
        # phase triplet (angle, fraction, elongation) reports 0.0.
        return (0.0, 0.0, 0.0, diameter, magnitude) + (0.0,) * 15

    # ------------------------------------------------------------------
    # Geocentric ecliptic positions of target and Sun
    # ------------------------------------------------------------------
    target_pos, _ = _leb_calc(tjd_ut, ipl, base_flags)
    sun_pos, _ = _leb_calc(tjd_ut, SUN, base_flags)

    target_lon = float(target_pos[0])
    target_lat = float(target_pos[1])
    geo_dist = float(target_pos[2])

    sun_lon = float(sun_pos[0])
    sun_lat = float(sun_pos[1])
    sun_dist = float(sun_pos[2])  # Earth-Sun distance in AU

    # ------------------------------------------------------------------
    # Elongation (angular separation on the ecliptic sphere)
    # ------------------------------------------------------------------
    elongation = angular_separation(target_lon, target_lat, sun_lon, sun_lat)

    # ------------------------------------------------------------------
    # Heliocentric distance of the target
    # ------------------------------------------------------------------
    helio_pos, _ = _leb_calc(tjd_ut, ipl, (base_flags & ~FLG_TOPOCTR) | FLG_HELCTR)
    helio_dist = float(helio_pos[2])

    # ------------------------------------------------------------------
    # Phase angle
    # ------------------------------------------------------------------
    # Phase angle from the apparent geocentric triangle: both the body→Earth
    # and body→Sun legs come from the apparent geocentric Sun and body vectors
    # (geometric under FLG_TRUEPOS). This mirrors the Skyfield _calc_pheno()
    # composition exactly, so the two backends agree to <0.01"; see the note
    # there and docs/comparison/intentional-divergences.md for why the
    # reference's own (unpublished) composition is not tracked.
    try:
        # body→Earth: antipode of the body's apparent ecliptic direction
        _lat_r = math.radians(target_lat)
        _lon_r = math.radians(target_lon)
        be = (
            -math.cos(_lat_r) * math.cos(_lon_r),
            -math.cos(_lat_r) * math.sin(_lon_r),
            -math.sin(_lat_r),
        )
        # Sun−body from the spherical outputs (apparent, or geometric under
        # FLG_TRUEPOS, since target_pos/sun_pos honour the flag above)
        _slat = math.radians(sun_lat)
        _slon = math.radians(sun_lon)
        s_vec = (
            sun_dist * math.cos(_slat) * math.cos(_slon),
            sun_dist * math.cos(_slat) * math.sin(_slon),
            sun_dist * math.sin(_slat),
        )
        b_vec = (
            -geo_dist * be[0],
            -geo_dist * be[1],
            -geo_dist * be[2],
        )
        bs = (
            s_vec[0] - b_vec[0],
            s_vec[1] - b_vec[1],
            s_vec[2] - b_vec[2],
        )
        dot = bs[0] * be[0] + bs[1] * be[1] + bs[2] * be[2]
        mag_bs = math.sqrt(bs[0] ** 2 + bs[1] ** 2 + bs[2] ** 2)
        mag_be = math.sqrt(be[0] ** 2 + be[1] ** 2 + be[2] ** 2)
        if mag_bs > 0 and mag_be > 0:
            cos_pa = max(-1.0, min(1.0, dot / (mag_bs * mag_be)))
            phase_angle = math.degrees(math.acos(cos_pa))
        else:
            phase_angle = 0.0
    except (KeyError, ValueError):
        # Fallback: scalar law-of-cosines (sufficient for planets)
        if geo_dist > 0 and helio_dist > 0:
            cos_pa = (geo_dist**2 + helio_dist**2 - sun_dist**2) / (
                2.0 * geo_dist * helio_dist
            )
            cos_pa = max(-1.0, min(1.0, cos_pa))
            phase_angle = math.degrees(math.acos(cos_pa))
        else:
            phase_angle = 0.0

    # ------------------------------------------------------------------
    # Phase (illuminated fraction)
    # ------------------------------------------------------------------
    phase = (1.0 + math.cos(math.radians(phase_angle))) / 2.0

    # ------------------------------------------------------------------
    # Apparent diameter
    # ------------------------------------------------------------------
    if ipl in _ASTEROID_DIAMETER_KM:
        body_radius_km = _ASTEROID_DIAMETER_KM[ipl] / 2.0
    else:
        body_radius_km = _BODY_RADIUS_KM.get(ipl, 1000.0)
    diameter = _calc_apparent_diameter(body_radius_km, geo_dist)

    # ------------------------------------------------------------------
    # Visual magnitude
    # ------------------------------------------------------------------
    if ipl == MOON:
        magnitude = _calc_moon_magnitude(phase_angle, geo_dist, helio_dist)
    elif ipl in _ASTEROID_HG:
        magnitude = _asteroid_hg_magnitude(ipl, phase_angle, helio_dist, geo_dist)
    else:
        # For Saturn we need ecliptic coordinates for ring tilt
        geo_lon = float(target_pos[0])
        geo_lat = float(target_pos[1])
        helio_lon = float(helio_pos[0])
        helio_lat = float(helio_pos[1])

        # TT approximation for Saturn ring / Neptune secular brightness
        from .time_utils import deltat

        tjd = tjd_ut + deltat(tjd_ut)

        magnitude = _calc_planet_magnitude(
            ipl,
            helio_dist,
            geo_dist,
            phase_angle,
            geo_lon,
            geo_lat,
            helio_lon,
            helio_lat,
            tjd,
        )

    # Slot [5]: the Moon's horizontal parallax (0.0 for every other body). The
    # reference reports the geocentric equatorial horizontal parallax when no
    # observer is set, and the actual geocentric->topocentric parallactic
    # displacement under FLG_TOPOCTR (both verified to 0.000000" against the
    # reference); neither depends on the request's FLG_TRUEPOS bit.
    slot5 = 0.0
    if ipl == MOON:
        if iflag & FLG_TOPOCTR:
            geo_app, _ = _leb_calc(
                tjd_ut, MOON, base_flags & ~(FLG_TOPOCTR | FLG_TRUEPOS)
            )
            topo_app, _ = _leb_calc(
                tjd_ut, MOON, (base_flags | FLG_TOPOCTR) & ~FLG_TRUEPOS
            )
            slot5 = angular_separation(
                float(geo_app[0]),
                float(geo_app[1]),
                float(topo_app[0]),
                float(topo_app[1]),
            )
        else:
            geom_pos, _ = _leb_calc(
                tjd_ut, MOON, (base_flags & ~FLG_TOPOCTR) | FLG_TRUEPOS
            )
            slot5 = _moon_horizontal_parallax_deg(float(geom_pos[2]))

    return (phase_angle, phase, elongation, diameter, magnitude, slot5) + (0.0,) * 14


def _calc_pheno_asteroid(t, ipl: int, iflag: int) -> Tuple[float, ...]:
    """Phenomena for the curated asteroids (Chiron..Vesta, Pholus).

    Positions come from calc_ut (whatever backend serves the body); the
    phase angle uses the law of cosines on the geocentric, heliocentric
    and Earth-Sun distances, and the magnitude the IAU H-G system with
    NASA/JPL SBDB photometric parameters.
    """
    from .utils import angular_separation

    jd_ut = t.ut1
    base_flags = iflag & (FLG_TRUEPOS | FLG_NONUT | FLG_TOPOCTR)

    target_pos, _ = calc_ut(jd_ut, ipl, base_flags)
    sun_pos, _ = calc_ut(jd_ut, SUN, base_flags)
    helio_pos, _ = calc_ut(jd_ut, ipl, (base_flags & ~FLG_TOPOCTR) | FLG_HELCTR)

    geo_dist = float(target_pos[2])
    helio_dist = float(helio_pos[2])
    sun_dist = float(sun_pos[2])

    elongation = angular_separation(
        float(target_pos[0]), float(target_pos[1]), float(sun_pos[0]), float(sun_pos[1])
    )

    if geo_dist > 0 and helio_dist > 0:
        # Vector apparent-triangle phase angle (Sun-body-Earth angle at the
        # body), from the apparent geocentric and heliocentric directions -
        # the same convention as the planet path and the LEB backend. The
        # scalar law of cosines on the three apparent distances is NOT
        # consistent here (light time and aberration make them a
        # non-Euclidean triangle) and split the backends by up to ~10".
        def _vec(lon_deg: float, lat_deg: float, r: float) -> tuple:
            lo = math.radians(lon_deg)
            la = math.radians(lat_deg)
            return (
                r * math.cos(la) * math.cos(lo),
                r * math.cos(la) * math.sin(lo),
                r * math.sin(la),
            )

        b_vec = _vec(float(target_pos[0]), float(target_pos[1]), geo_dist)
        s_vec = _vec(float(sun_pos[0]), float(sun_pos[1]), sun_dist)
        be = (-b_vec[0], -b_vec[1], -b_vec[2])
        bs = (s_vec[0] - b_vec[0], s_vec[1] - b_vec[1], s_vec[2] - b_vec[2])
        cx = be[1] * bs[2] - be[2] * bs[1]
        cy = be[2] * bs[0] - be[0] * bs[2]
        cz = be[0] * bs[1] - be[1] * bs[0]
        phase_angle = math.degrees(
            math.atan2(
                math.sqrt(cx * cx + cy * cy + cz * cz),
                be[0] * bs[0] + be[1] * bs[1] + be[2] * bs[2],
            )
        )
    else:
        phase_angle = 0.0
    phase = (1.0 + math.cos(math.radians(phase_angle))) / 2.0

    diameter = _calc_apparent_diameter(_ASTEROID_DIAMETER_KM[ipl] / 2.0, geo_dist)
    magnitude = _asteroid_hg_magnitude(ipl, phase_angle, helio_dist, geo_dist)

    return (phase_angle, phase, elongation, diameter, magnitude) + (0.0,) * 15


# Pseudo- and degenerate-body phenomena tuples. Neither id has a physical
# disc seen from the Earth's centre, so no phase triangle exists for them:
# EARTH *is* the observer (its geocentric vector is the zero origin), and
# ECL_NUT is not a body at all — calc() returns nutation and obliquity in its
# first slots. Feeding either tuple to the geometric pipeline produced
# invented values (a zero-length vector blew the apparent diameter up to
# ~1.9e13 degrees under FLG_TRUEPOS, and the nutation angles were read as a
# longitude/latitude/distance triple).
#
# Measured reference behavior, invariant under every flag and date tested:
# EARTH returns an all-zero tuple with attr[3] = 180.0, and ECL_NUT returns
# NaN for the phase triplet. The reference's own attr[3] for ECL_NUT
# alternates between 180.0 and ~1e-10 with no dependence on the nutation
# values it was computed from, so it is not a documented quantity and is
# reported as 0.0 here (see docs/comparison/intentional-divergences.md).
_PHENO_EARTH = (0.0, 0.0, 0.0, 180.0) + (0.0,) * 16
_PHENO_ECL_NUT = (float("nan"),) * 3 + (0.0,) * 17


def _degenerate_pheno(planet: int) -> "Tuple[float, ...] | None":
    """Fixed phenomena tuple for a body with no phase geometry, else None."""
    from .constants import ECL_NUT

    if planet == EARTH:
        return _PHENO_EARTH
    if planet == ECL_NUT:
        return _PHENO_ECL_NUT
    return None


def pheno_ut(tjdut: float, planet: int, flags: int = FLG_SWIEPH) -> Tuple[float, ...]:
    """
    Compute planetary phenomena for Universal Time.

    Calculates phase angle, illuminated fraction, elongation, apparent diameter,
    and visual magnitude for planets and the Moon.

    Args:
        tjdut: Julian Day in Universal Time (UT1)
        planet: Planet/body ID (SUN, MOON, MERCURY, etc.)
        flags: Calculation flags (FLG_TRUEPOS, FLG_HELCTR, etc.)

    Returns:
        Tuple of 20 floats using the public compatibility layout:
            - [0]: Phase angle (Earth-planet-Sun) in degrees
            - [1]: Phase (illuminated fraction of disc, 0.0 to 1.0)
            - [2]: Elongation of planet from Sun in degrees
            - [3]: Apparent diameter of disc in degrees
            - [4]: Apparent visual magnitude
            - [5-19]: Reserved (0.0)

    Note:
        - For the Sun: phase angle = 0, phase = 1.0, elongation = 0
        - For the Moon: Hapke photometric model with opposition surge correction
        - Phase = 0.0 means new (fully dark), Phase = 1.0 means full (fully illuminated)
        - Elongation is measured from the Sun (0° = conjunction, 180° = opposition)

    Example:
        >>> attr = pheno_ut(2451545.0, MARS, 0)
        >>> print(f"Phase angle: {attr[0]:.2f}°")
        >>> print(f"Illumination: {attr[1]*100:.1f}%")
        >>> print(f"Elongation: {attr[2]:.2f}°")
        >>> print(f"Diameter: {attr[3]:.4f} deg")
        >>> print(f"Magnitude: {attr[4]:.2f}")
    """
    _degenerate = _degenerate_pheno(planet)
    if _degenerate is not None:
        return _degenerate

    # --- LEB fast path ---
    from .state import get_leb_reader

    _unsupported_pheno_flags = FLG_NOABERR | FLG_NOGDEFL
    reader = get_leb_reader()
    if reader is not None and not (flags & _unsupported_pheno_flags):
        try:
            return _calc_pheno_leb(tjdut, planet, flags)
        except (KeyError, ValueError) as _leb_err:
            # Fall back for missing bodies, out-of-range AND corrupted LEB
            # data, aligned with the calc_ut()/fixed-star convention.
            from .leb_reader import log_leb_fallback

            log_leb_fallback("pheno", _leb_err)
    # --- END LEB fast path ---

    ts = get_timescale()
    t = ts.ut1_jd(tjdut)
    return _calc_pheno(t, planet, flags)


def pheno(tjdet: float, planet: int, flags: int = FLG_SWIEPH) -> Tuple[float, ...]:
    """
    Compute planetary phenomena for Ephemeris Time (TT/ET).

    Same as pheno_ut() but takes Terrestrial Time instead of Universal Time.

    Args:
        tjdet: Julian Day in Ephemeris Time (TT/ET)
        planet: Planet/body ID (SUN, MOON, MERCURY, etc.)
        flags: Calculation flags

    Returns:
        Tuple of 20 floats (matching the reference API):
            - [0]: Phase angle, [1]: Phase, [2]: Elongation,
            - [3]: Diameter, [4]: Magnitude, [5-19]: Reserved (0.0)

    See Also:
        pheno_ut: Same function for Universal Time input
    """
    _degenerate = _degenerate_pheno(planet)
    if _degenerate is not None:
        return _degenerate

    # --- LEB fast path ---
    from .state import get_leb_reader

    _unsupported_pheno_flags_tt = FLG_NOABERR | FLG_NOGDEFL
    reader = get_leb_reader()
    if reader is not None and not (flags & _unsupported_pheno_flags_tt):
        try:
            from .time_utils import deltat

            tjd_ut = tjdet - deltat(tjdet)
            return _calc_pheno_leb(tjd_ut, planet, flags)
        except (KeyError, ValueError) as _leb_err:
            # Fall back for missing bodies, out-of-range AND corrupted LEB
            # data (see pheno_ut)
            from .leb_reader import log_leb_fallback

            log_leb_fallback("pheno", _leb_err)
    # --- END LEB fast path ---

    ts = get_timescale()
    t = ts.tt_jd(tjdet)
    return _calc_pheno(t, planet, flags)


_PHENO_POINT_BODIES = frozenset(
    (MEAN_NODE, TRUE_NODE, MEAN_APOG, OSCU_APOG, INTP_APOG, INTP_PERG)
)


def _pheno_positional(jd_ut: float, ipl: int, iflag: int) -> Tuple[float, ...]:
    """Phase triplet for point and model bodies without a physical disc.

    Measured reference behavior: the lunar node/apsis points report only
    their elongation from the Sun (phase angle and illuminated fraction stay
    0.0), while fictitious and registry-served bodies report the full
    geometric phase triplet from the apparent geocentric triangle. Apparent
    diameter and magnitude stay 0.0 for both families. Unknown ids propagate
    the typed error from the position pipeline.
    """
    pos_flags = iflag & (FLG_TOPOCTR | FLG_TRUEPOS)
    body = calc_ut(jd_ut, ipl, pos_flags)[0]
    sun = calc_ut(jd_ut, SUN, pos_flags)[0]

    def _cart(lon_deg: float, lat_deg: float, dist: float) -> Tuple[float, ...]:
        lon = math.radians(lon_deg)
        lat = math.radians(lat_deg)
        r = dist if dist > 0.0 else 1.0
        return (
            r * math.cos(lat) * math.cos(lon),
            r * math.cos(lat) * math.sin(lon),
            r * math.sin(lat),
        )

    P = _cart(body[0], body[1], body[2])
    S = _cart(sun[0], sun[1], sun[2])
    mag_p = math.sqrt(sum(c * c for c in P))
    mag_s = math.sqrt(sum(c * c for c in S))
    cos_el = sum(a * b for a, b in zip(P, S)) / (mag_p * mag_s)
    elongation = math.degrees(math.acos(max(-1.0, min(1.0, cos_el))))

    if ipl in _PHENO_POINT_BODIES:
        if ipl in (21, 22):
            # Interpolated apsides: the reference marks the inapplicable
            # phase slots NaN (same convention as its nod_aps output).
            return (float("nan"), float("nan"), elongation) + (0.0,) * 17
        return (0.0, 0.0, elongation) + (0.0,) * 17

    # Apparent geocentric triangle, like the disc bodies.
    ps = tuple(s - p for s, p in zip(S, P))
    pe = tuple(-p for p in P)
    mag_ps = math.sqrt(sum(c * c for c in ps))
    mag_pe = math.sqrt(sum(c * c for c in pe))
    if mag_ps > 0.0 and mag_pe > 0.0:
        cos_pa = sum(a * b for a, b in zip(ps, pe)) / (mag_ps * mag_pe)
        phase_angle = math.degrees(math.acos(max(-1.0, min(1.0, cos_pa))))
    else:
        phase_angle = 0.0
    phase = (1.0 + math.cos(math.radians(phase_angle))) / 2.0
    return (phase_angle, phase, elongation) + (0.0,) * 17


def _calc_pheno(t, ipl: int, iflag: int) -> Tuple[float, ...]:
    """
    Internal function to compute planetary phenomena.

    Calculates:
    1. Phase angle: angle Sun-Planet-Earth (for inner planets) or Earth-Planet-Sun (for outer)
    2. Phase: illuminated fraction using (1 + cos(phase_angle)) / 2
    3. Elongation: angular separation between planet and Sun as seen from Earth
    4. Apparent diameter: angular size of planet's disc
    5. Visual magnitude: brightness using empirical formulas

    Args:
        t: Skyfield Time object
        ipl: Planet/body ID
        iflag: Calculation flags

    Returns:
        Tuple of 20 floats (bare tuple, matching the reference API).
    """

    from .cache import get_cached_observer_at

    planets = _get_computation_ephemeris()

    # Built-in asteroids requested via the AST_OFFSET + number alias resolve to
    # their dedicated id (AST_OFFSET + 4 -> Vesta, AST_OFFSET + 5145 -> Pholus),
    # matching calc_ut's position aliasing and the built-in id's pheno.
    # Registry-served numbered asteroids keep their own id.
    ipl = _remap_ast_offset(ipl)

    # The reference ignores FLG_HELCTR in pheno: every quantity is
    # Earth-based regardless of the observer flag (measured behavior).
    iflag &= ~FLG_HELCTR

    earth = planets["earth"]
    sun = planets["sun"]

    # The reference DOES honor FLG_TOPOCTR in pheno: every distance-derived
    # quantity (apparent diameter, phase angle, elongation, magnitude) is taken
    # from the observer's topocentric place rather than the geocentre (measured
    # behavior; the effect is ~arcseconds of diameter and up to ~0.4 deg of
    # phase angle for the Moon, sub-milli for the far planets). Build the
    # topocentric observer once and use it for both the target and Sun legs.
    if iflag & FLG_TOPOCTR:
        observer_topo = get_topo()
        if observer_topo is None:
            from .exceptions import ConfigurationError

            raise ConfigurationError(
                "FLG_TOPOCTR requires a geographic position: call "
                "set_topo(lon, lat, alt) first",
                missing_config="observer_location",
                suggestion="Call set_topo(lon, lat, alt) first",
            )
        observer = earth + observer_topo
    else:
        observer = earth

    # Initialize return values
    phase_angle = 0.0
    phase = 1.0
    elongation = 0.0
    diameter = 0.0
    magnitude = 0.0

    # Special case: Sun
    if ipl == SUN:
        # Phase quantities are inapplicable to the self-luminous Sun: every
        # slot of the phase triplet reports 0.0 (phase angle, illuminated
        # fraction, elongation).
        # Get Sun distance from the (possibly topocentric) observer.
        sun_pos = get_cached_observer_at(observer, t).observe(sun).apparent()
        _, _, sun_dist = sun_pos.radec()

        phase_angle = 0.0
        phase = 0.0
        elongation = 0.0

        # Apparent diameter of Sun based on physical radius
        sun_dist_au = float(sun_dist.au)
        sun_radius_km = _BODY_RADIUS_KM.get(SUN, 695700.0)
        diameter = _calc_apparent_diameter(sun_radius_km, sun_dist_au)

        # Sun magnitude (V(1,0) = -26.86 at 1 AU, Mallama & Hilton 2018)
        magnitude = (
            -26.86 + 5.0 * math.log10(sun_dist_au) if sun_dist_au > 0 else -26.86
        )

        attr = (phase_angle, phase, elongation, diameter, magnitude) + (0.0,) * 15
        return attr

    # Get geocentric (or topocentric under FLG_TOPOCTR) position of target
    if ipl == MOON:
        target = planets["moon"]
    elif ipl in _PLANET_MAP:
        target_name = _PLANET_MAP[ipl]
        # Try planet center first, fall back to barycenter if not available
        try:
            target = planets[target_name]
        except KeyError:
            if target_name in _PLANET_FALLBACK:
                target = planets[_PLANET_FALLBACK[target_name]]
            else:
                raise
    elif ipl in _ASTEROID_HG:
        # Curated asteroids: positions via calc_ut (SPK/Keplerian backend),
        # photometry via the IAU H-G system (reference behavior).
        return _calc_pheno_asteroid(t, ipl, iflag)
    else:
        # Point and model bodies without a physical disc: measured reference
        # behavior fills the geometric phase triplet from their positions.
        return _pheno_positional(float(t.ut1), ipl, iflag)

    # Get observed positions (observer is geocentric or, under FLG_TOPOCTR,
    # the topocentre; set once near the top of this function).
    obs_at_t = get_cached_observer_at(observer, t)
    if iflag & FLG_TRUEPOS:
        # True geometric positions at t on BOTH legs: observe() would still
        # retard the target by light time, which the reference removes.
        target_pos_geo = (target - observer).at(t)
        sun_pos_geo = (sun - observer).at(t)
    else:
        # Apparent positions
        target_pos_geo = obs_at_t.observe(target).apparent()
        sun_pos_geo = obs_at_t.observe(sun).apparent()

    # Get heliocentric position of target for phase calculations
    target_helio = get_cached_observer_at(sun, t).observe(target)
    target_helio_dist = math.sqrt(sum(x**2 for x in target_helio.position.au))

    # Get geocentric distance
    target_geo_dist = math.sqrt(sum(x**2 for x in target_pos_geo.position.au))

    # Special handling for Moon
    if ipl == MOON:
        # Sun position from Earth, same flavor (geometric/apparent) as the
        # Moon leg above.
        sun_from_earth = sun_pos_geo

        # Get RA/Dec of Moon and Sun
        moon_ra, moon_dec, moon_dist = target_pos_geo.radec()
        sun_ra, sun_dec, sun_dist = sun_from_earth.radec()

        # Elongation: angular distance between Moon and Sun
        # Using spherical trigonometry
        moon_ra_rad = moon_ra.radians
        moon_dec_rad = moon_dec.radians
        sun_ra_rad = sun_ra.radians
        sun_dec_rad = sun_dec.radians

        cos_elong = math.sin(moon_dec_rad) * math.sin(sun_dec_rad) + math.cos(
            moon_dec_rad
        ) * math.cos(sun_dec_rad) * math.cos(moon_ra_rad - sun_ra_rad)
        cos_elong = max(-1.0, min(1.0, cos_elong))
        elongation = math.degrees(math.acos(cos_elong))

        # Phase angle for Moon using 3D vector approach
        # The phase angle is the angle at the Moon vertex between the
        # Sun-Moon and Earth-Moon directions. Using position vectors
        # is more numerically stable than law-of-cosines for the Moon's
        # extremely elongated triangle (Sun~1AU, Moon~0.003AU from Earth).
        r_moon = float(moon_dist.au)  # Earth-Moon distance (native float, not numpy)

        # Vectors in geocentric frame (Earth at origin)
        M = target_pos_geo.position.au  # Moon position
        S = sun_from_earth.position.au  # Sun position

        # Sun-Moon distance (for the magnitude model below)
        vec_moon_to_sun = S - M
        vec_moon_to_earth = -M
        mag_ms = math.sqrt(sum(x**2 for x in vec_moon_to_sun))
        mag_me = math.sqrt(sum(x**2 for x in vec_moon_to_earth))

        # Phase angle from the apparent geocentric triangle: both legs are
        # taken from the apparent geocentric Sun and Moon vectors (geometric
        # under FLG_TRUEPOS). See the note in the planet branch below for why
        # this composition is used in place of the reference's.
        bs = vec_moon_to_sun

        dot_prod = sum(a * b for a, b in zip(bs, vec_moon_to_earth))
        mag_bs = math.sqrt(sum(x**2 for x in bs))

        if mag_bs > 0 and mag_me > 0:
            cos_phase = dot_prod / (mag_bs * mag_me)
            cos_phase = max(-1.0, min(1.0, cos_phase))
            phase_angle = math.degrees(math.acos(cos_phase))
        else:
            phase_angle = 180.0 - elongation

        # Phase (illuminated fraction)
        phase = (1.0 + math.cos(math.radians(phase_angle))) / 2.0

        # Moon's apparent diameter based on physical radius
        moon_radius_km = _BODY_RADIUS_KM.get(MOON, 1737.4)
        diameter = _calc_apparent_diameter(moon_radius_km, r_moon)

        # Moon's magnitude using piecewise Allen/Samaha photometric model
        # mag_ms is the Sun-Moon distance (heliocentric distance of Moon) in AU
        magnitude = _calc_moon_magnitude(phase_angle, r_moon, mag_ms)

        # Slot [5]: the Moon's horizontal parallax. Without an observer the
        # reference reports the geocentric equatorial horizontal parallax from
        # the geometric geocentric distance; under FLG_TOPOCTR it reports the
        # actual geocentric->topocentric parallactic displacement (the on-sky
        # angle between the geocentric-apparent and topocentric-apparent Moon
        # directions). Both were verified to 0.000000" against the reference and
        # are independent of the request's FLG_TRUEPOS bit.
        if iflag & FLG_TOPOCTR:
            # Both legs are the APPARENT place (even under FLG_TRUEPOS, where
            # target_pos_geo would be geometric): the reference's topocentric
            # parallax is the apparent geo->topo displacement.
            geo_app = get_cached_observer_at(earth, t).observe(target).apparent()
            topo_app = obs_at_t.observe(target).apparent()
            v_geo = geo_app.position.au
            v_topo = topo_app.position.au
            dot_p = sum(a * b for a, b in zip(v_geo, v_topo))
            n_geo = math.sqrt(sum(a * a for a in v_geo))
            n_topo = math.sqrt(sum(a * a for a in v_topo))
            if n_geo > 0 and n_topo > 0:
                cos_par = max(-1.0, min(1.0, dot_p / (n_geo * n_topo)))
                horizontal_parallax = math.degrees(math.acos(cos_par))
            else:
                horizontal_parallax = 0.0
        else:
            moon_geom = (target - earth).at(t)
            dist_geom_au = math.sqrt(sum(x**2 for x in moon_geom.position.au))
            horizontal_parallax = _moon_horizontal_parallax_deg(dist_geom_au)

        attr = (
            phase_angle,
            phase,
            elongation,
            diameter,
            magnitude,
            horizontal_parallax,
        ) + (0.0,) * 14
        return attr

    # For planets: calculate elongation, phase angle, etc.
    # Get RA/Dec of planet and Sun
    planet_ra, planet_dec, planet_dist = target_pos_geo.radec()
    sun_ra, sun_dec, sun_dist = sun_pos_geo.radec()

    # Elongation: angular distance between planet and Sun
    planet_ra_rad = planet_ra.radians
    planet_dec_rad = planet_dec.radians
    sun_ra_rad = sun_ra.radians
    sun_dec_rad = sun_dec.radians

    cos_elong = math.sin(planet_dec_rad) * math.sin(sun_dec_rad) + math.cos(
        planet_dec_rad
    ) * math.cos(sun_dec_rad) * math.cos(planet_ra_rad - sun_ra_rad)
    cos_elong = max(-1.0, min(1.0, cos_elong))
    elongation = math.degrees(math.acos(cos_elong))

    # Phase angle: angle at the planet vertex between the planet→Sun and
    # planet→Earth directions, both taken from the apparent geocentric triangle
    # (geometric under FLG_TRUEPOS): planet→Earth = −P and planet→Sun = S − P,
    # from the apparent geocentric Sun and planet vectors already computed.
    #
    # The reference ephemeris composes this angle from a proprietary blend of
    # apparent/astrometric legs that no published single- or double-light-time
    # formula reproduces to <0.1" for every body: the per-body best fits
    # contradict one another and the angle even responds to aberration
    # (toggling NOABERR shifts it ~13-19"). Rather than track an unpublished
    # composition, this uses the fully documented apparent geocentric triangle.
    # The residual against the reference is arc-second level for the inner
    # planets (Mercury worst) and model-dependent; FLG_TRUEPOS stays exact.
    # See docs/comparison/intentional-divergences.md.
    P = target_pos_geo.position.au  # Planet position (geocentric)
    S = sun_pos_geo.position.au  # Sun position (geocentric)

    vec_planet_to_earth = -P
    mag_pe = math.sqrt(sum(x**2 for x in vec_planet_to_earth))
    vec_planet_to_sun = S - P

    dot_prod = sum(a * b for a, b in zip(vec_planet_to_sun, vec_planet_to_earth))
    mag_ps = math.sqrt(sum(x**2 for x in vec_planet_to_sun))

    if mag_ps > 0 and mag_pe > 0:
        cos_phase = dot_prod / (mag_ps * mag_pe)
        cos_phase = max(-1.0, min(1.0, cos_phase))
        phase_angle = math.degrees(math.acos(cos_phase))
    else:
        phase_angle = 0.0

    # Phase (illuminated fraction)
    phase = (1.0 + math.cos(math.radians(phase_angle))) / 2.0

    # Apparent diameter based on physical radius
    body_radius_km = _BODY_RADIUS_KM.get(
        ipl, 1000.0
    )  # Default 1000 km for unknown bodies
    diameter = _calc_apparent_diameter(body_radius_km, target_geo_dist)

    # Visual magnitude
    # For Saturn, we need ecliptic coordinates for ring tilt calculation
    # tjd is needed for Saturn (ring tilt) and Neptune (secular brightness)
    geo_lon = 0.0
    geo_lat = 0.0
    helio_lon = 0.0
    helio_lat = 0.0
    tjd = t.tt

    if ipl in (SATURN, URANUS):
        # Saturn needs the ecliptic coordinates for the ring-tilt term and
        # Uranus for the sub-Earth photometric latitude of Mallama & Hilton
        # (2018) Eq. 15 — without them the Uranus term degenerates to a
        # constant and the two backends split by tens of millimagnitudes.
        # Get geocentric ecliptic coordinates
        try:
            geo_ecl_lat, geo_ecl_lon, _ = target_pos_geo.frame_latlon(ecliptic_frame)
            geo_lon = geo_ecl_lon.degrees
            geo_lat = geo_ecl_lat.degrees
        except (AttributeError, ValueError, TypeError):
            geo_lon = 0.0
            geo_lat = 0.0

        # Get heliocentric ecliptic coordinates
        try:
            helio_ecl_lat, helio_ecl_lon, _ = target_helio.frame_latlon(ecliptic_frame)
            helio_lon = helio_ecl_lon.degrees
            helio_lat = helio_ecl_lat.degrees
        except (AttributeError, ValueError, TypeError):
            helio_lon = 0.0
            helio_lat = 0.0

    magnitude = _calc_planet_magnitude(
        ipl,
        target_helio_dist,
        target_geo_dist,
        phase_angle,
        geo_lon,
        geo_lat,
        helio_lon,
        helio_lat,
        tjd,
    )

    # Return tuple with at least 20 elements (reference-API compatibility)
    attr = (phase_angle, phase, elongation, diameter, magnitude) + (0.0,) * 15
    return attr


def _uranus_photometric_latitude(
    geo_lon: float, geo_lat: float, helio_lon: float, helio_lat: float
) -> float:
    """Photometric latitude phi' of Uranus for Mallama & Hilton (2018).

    phi' is the average of the absolute planetographic sub-Earth and sub-solar
    latitudes (Mallama & Hilton 2018, "Computing Apparent Planetary Magnitudes
    for The Astronomical Almanac"). Each planetocentric latitude is the angle
    between Uranus' IAU rotational pole and the Uranus->Earth / Uranus->Sun
    direction, then converted to planetographic latitude via Eq. 13
    (tan phi' = tan phi / (1 - f)^2).

    The IAU J2000 pole vector is dotted against apparent ecliptic-of-date
    directions; the resulting frame mixing shifts phi' by <=0.7 deg near the
    span extremes, i.e. <=0.6 mmag through the -8.4e-4 mag/deg coefficient,
    which is far below the model's structural agreement with any reference.

    Args:
        geo_lon: Geocentric apparent ecliptic longitude of Uranus (degrees).
        geo_lat: Geocentric apparent ecliptic latitude of Uranus (degrees).
        helio_lon: Heliocentric ecliptic longitude of Uranus (degrees).
        helio_lat: Heliocentric ecliptic latitude of Uranus (degrees).

    Returns:
        Photometric latitude phi' in degrees (>= 0).
    """
    # Uranus north-pole unit vector: equatorial J2000 -> J2000 ecliptic.
    eps = math.radians(23.43929111)  # J2000 mean obliquity (IAU 2006)
    cos_d = math.cos(_URANUS_POLE_DEC)
    pe = (
        cos_d * math.cos(_URANUS_POLE_RA),
        cos_d * math.sin(_URANUS_POLE_RA),
        math.sin(_URANUS_POLE_DEC),
    )
    pole = (
        pe[0],
        pe[1] * math.cos(eps) + pe[2] * math.sin(eps),
        -pe[1] * math.sin(eps) + pe[2] * math.cos(eps),
    )

    def _sub_latitude(lon_deg: float, lat_deg: float) -> float:
        # Direction Uranus->observer is the antipode of observer->Uranus.
        lon = math.radians(lon_deg)
        lat = math.radians(lat_deg)
        d = (
            -math.cos(lat) * math.cos(lon),
            -math.cos(lat) * math.sin(lon),
            -math.sin(lat),
        )
        sin_phi = max(-1.0, min(1.0, d[0] * pole[0] + d[1] * pole[1] + d[2] * pole[2]))
        phi_centric = math.asin(sin_phi)
        # Eq. 13: planetocentric -> planetographic latitude.
        phi_graphic = math.atan2(math.tan(phi_centric), _URANUS_ONE_MINUS_F_SQ)
        return abs(math.degrees(phi_graphic))

    phi_earth = _sub_latitude(geo_lon, geo_lat)
    phi_sun = _sub_latitude(helio_lon, helio_lat)
    return (phi_earth + phi_sun) / 2.0


def _calc_planet_magnitude(
    ipl: int,
    helio_dist: float,
    geo_dist: float,
    phase_angle: float,
    geo_lon: float = 0.0,
    geo_lat: float = 0.0,
    helio_lon: float = 0.0,
    helio_lat: float = 0.0,
    tjd: float = 0.0,
) -> float:
    """
    Calculate visual magnitude of a planet.

    Uses Mallama 2018 formulas for Mercury, Venus, Mars, Jupiter, Saturn
    for reference-API compatibility. These formulas are from:
    A. Mallama, J. Hilton, "Computing Apparent Planetary Magnitudes for
    The Astronomical Almanac" (2018).

    Args:
        ipl: Planet ID
        helio_dist: Heliocentric distance in AU
        geo_dist: Geocentric distance in AU
        phase_angle: Phase angle in degrees
        geo_lon: Geocentric ecliptic longitude in degrees (for Saturn)
        geo_lat: Geocentric ecliptic latitude in degrees (for Saturn)
        helio_lon: Heliocentric ecliptic longitude in degrees (for Saturn)
        helio_lat: Heliocentric ecliptic latitude in degrees (for Saturn)
        tjd: Julian day (for Saturn ring tilt calculation)

    Returns:
        Visual magnitude (smaller = brighter)
    """
    # Distance factor: 5 * log10(r * d)
    if helio_dist > 0 and geo_dist > 0:
        dist_factor = 5.0 * math.log10(helio_dist * geo_dist)
    else:
        dist_factor = 0.0

    a = phase_angle  # Phase angle in degrees

    # Mercury - Mallama 2018, 6th order polynomial
    if ipl == MERCURY:
        a2 = a * a
        a3 = a2 * a
        a4 = a3 * a
        a5 = a4 * a
        a6 = a5 * a
        magnitude = (
            -0.613
            + a * 6.3280e-02
            - a2 * 1.6336e-03
            + a3 * 3.3644e-05
            - a4 * 3.4265e-07
            + a5 * 1.6893e-09
            - a6 * 3.0334e-12
        )
        magnitude += dist_factor
        return magnitude

    # Venus - Mallama 2018, piecewise polynomial
    if ipl == VENUS:
        a2 = a * a
        a3 = a2 * a
        a4 = a3 * a
        if a <= 163.7:
            magnitude = (
                -4.384
                - a * 1.044e-03
                + a2 * 3.687e-04
                - a3 * 2.814e-06
                + a4 * 8.938e-09
            )
        else:
            magnitude = 236.05828 - a * 2.81914e00 + a2 * 8.39034e-03
        magnitude += dist_factor
        return magnitude

    # Mars - Mallama 2018, piecewise polynomial
    if ipl == MARS:
        a2 = a * a
        if a <= 50.0:
            magnitude = -1.601 + a * 0.02267 - a2 * 0.0001302
        else:
            magnitude = -0.367 - a * 0.02573 + a2 * 0.0003445
        magnitude += dist_factor
        return magnitude

    # Jupiter - Mallama 2018
    if ipl == JUPITER:
        a2 = a * a
        magnitude = -9.395 - a * 3.7e-04 + a2 * 6.16e-04
        magnitude += dist_factor
        return magnitude

    # Saturn - Mallama 2018 with ring tilt from Meeus
    if ipl == SATURN:
        # Ring tilt calculation from Meeus, p. 301ff
        # T is centuries from J2000
        T = (tjd - _J2000) / 36525.0

        # Ring plane inclination and ascending node
        incl = math.radians(28.075216 - 0.012998 * T + 0.000004 * T * T)
        omega = math.radians(169.508470 + 1.394681 * T + 0.000412 * T * T)

        # sinB is the sine of the ring tilt angle as seen from Earth/Sun
        # B is the "mean tilt of the ring plane to the Earth and Sun"
        # From Meeus formulae for ring visibility

        # Geocentric position contribution
        sin_B = math.sin(incl) * math.cos(math.radians(geo_lat)) * math.sin(
            math.radians(geo_lon) - omega
        ) - math.cos(incl) * math.sin(math.radians(geo_lat))

        # Heliocentric position contribution
        sin_B2 = math.sin(incl) * math.cos(math.radians(helio_lat)) * math.sin(
            math.radians(helio_lon) - omega
        ) - math.cos(incl) * math.sin(math.radians(helio_lat))

        # Mean of the two tilt angles
        sin_B = abs(math.sin((math.asin(sin_B) + math.asin(sin_B2)) / 2.0))

        # Mallama 2018 formula for Saturn with ring contribution
        magnitude = (
            -8.914 - 1.825 * sin_B + 0.026 * a - 0.378 * sin_B * math.exp(-2.25 * a)
        )
        magnitude += dist_factor
        return magnitude

    # Pluto - Mallama & Hilton 2018 formula
    # From "Computing Apparent Planetary Magnitudes for The Astronomical Almanac"
    # V(1,0) = -1.024 ± 0.003 mag (absolute magnitude at r=d=1 AU, α=0°)
    # Phase coefficient β = 0.0362 ± 0.0004 mag/degree
    # Formula: V = V(1,0) + 5*log10(r*d) + β*α
    # Note: Pluto also has a rotational light curve amplitude of ~±0.15 mag
    # with period 6.3872 days, but this requires sub-observer longitude data.
    if ipl == PLUTO:
        V0 = -1.024  # Absolute magnitude V(1,0)
        beta = 0.0362  # Phase coefficient in mag/degree
        phase_correction = beta * a  # Linear phase correction
        magnitude = V0 + dist_factor + phase_correction
        return magnitude

    # Uranus - Mallama & Hilton (2018) complete photometric law (Eq. 15)
    # From "Computing Apparent Planetary Magnitudes for The Astronomical Almanac"
    # (Astronomy and Computing 25, 2018, 10-24; arXiv:1808.01973), which folds
    # the sub-observer-latitude brightness variation and the phase-angle
    # dependence into:
    #     V = 5*log10(r*d) - 7.110 - 8.4e-4*phi' + 6.587e-3*a + 1.045e-4*a^2
    # where phi' is the photometric latitude (average of the absolute
    # planetographic sub-Earth and sub-solar latitudes; Eq. 14 gives the
    # phase-free geocentric form). Uranus brightens as a pole turns toward the
    # observer because its polar regions are depleted in light-absorbing
    # methane (Schmude et al. 2015). Valid to a = 154 deg; geocentric a <= 3.1
    # deg. The residual model divergence versus the reference API is documented
    # in docs/comparison/intentional-divergences.md ("Outer-planet visual
    # magnitude photometry").
    if ipl == URANUS:
        phi_prime = _uranus_photometric_latitude(geo_lon, geo_lat, helio_lon, helio_lat)
        magnitude = (
            dist_factor
            - 7.110
            - 8.4e-04 * phi_prime
            + 6.587e-03 * a
            + 1.045e-04 * a * a
        )
        return magnitude

    # Neptune - secular brightness variation
    # Neptune's albedo has been increasing since ~1980 due to seasonal changes
    # over its 165-year orbital period. The absolute magnitude V(1,0) transitions
    # linearly from -6.89 (pre-1980) to -7.00 (by J2000.0).
    # Reference: Lockwood & Thompson (1991), Sromovsky et al. (2003)
    if ipl == NEPTUNE:
        year = 2000.0 + (tjd - _J2000) / 365.25
        if year >= 2000.0:
            V0 = -7.00
        elif year <= 1980.0:
            V0 = -6.89
        else:
            # Linear interpolation: -6.89 at 1980 to -7.00 at 2000
            V0 = -6.89 + (year - 1980.0) * (-0.11 / 20.0)
        magnitude = V0 + dist_factor
        return magnitude

    # Outer planets using simplified formula
    if ipl in _PLANET_MAG_PARAMS:
        V0, B1, B2, B3 = _PLANET_MAG_PARAMS[ipl]
        phase_factor = B1 * a + B2 * a**2 + B3 * a**3
        magnitude = V0 + dist_factor + phase_factor
        return magnitude

    # Unknown planet - return approximate magnitude
    H = 10.0  # Assumed absolute magnitude
    if helio_dist > 0 and geo_dist > 0:
        return H + 5.0 * math.log10(helio_dist * geo_dist)
    return H


# Aliases for reference API compatibility
pheno_ut = pheno_ut
pheno = pheno


# =============================================================================
# Elongation Helper Functions
# =============================================================================


def get_elongation_from_sun(
    tjd_ut: float, ipl: int, iflag: int = 0
) -> Tuple[float, bool]:
    """
    Calculate the elongation of a planet from the Sun with morning/evening star distinction.

    This function returns the elongation (angular separation from the Sun) and
    determines whether the planet is a morning star (western elongation) or
    evening star (eastern elongation).

    Elongation convention:
        - Eastern elongation (positive): Planet is east of Sun, visible after sunset
          (evening star)
        - Western elongation (negative): Planet is west of Sun, visible before sunrise
          (morning star)

    Args:
        tjd_ut: Julian Day in Universal Time (UT1)
        ipl: Planet/body ID (MERCURY, VENUS, MARS, etc.)
        iflag: Calculation flags (default: 0)

    Returns:
        Tuple containing:
            - elongation: Signed elongation in degrees
                - Positive = eastern elongation (evening star)
                - Negative = western elongation (morning star)
            - is_evening_star: True if planet is an evening star (eastern elongation),
                               False if morning star (western elongation)

    Note:
        This function is most useful for Mercury and Venus (inferior planets),
        which alternate between morning and evening star status. Superior planets
        (Mars, Jupiter, Saturn, etc.) can also be classified this way, though
        they are visible for longer periods.

    Example:
        >>> from libephemeris import get_elongation_from_sun, VENUS
        >>> elongation, is_evening = get_elongation_from_sun(2451545.0, VENUS)
        >>> if is_evening:
        ...     print(f"Venus is an evening star at {abs(elongation):.1f}° eastern elongation")
        ... else:
        ...     print(f"Venus is a morning star at {abs(elongation):.1f}° western elongation")
    """
    # Get planet and Sun positions
    planet_pos, _ = calc_ut(tjd_ut, ipl, iflag)
    sun_pos, _ = calc_ut(tjd_ut, SUN, iflag)

    planet_lon = float(planet_pos[0])
    sun_lon = float(sun_pos[0])

    # Calculate the longitude difference (planet - Sun)
    # Normalize to -180 to +180 range
    lon_diff = planet_lon - sun_lon
    if lon_diff > 180.0:
        lon_diff -= 360.0
    elif lon_diff < -180.0:
        lon_diff += 360.0

    # Positive difference = planet is east of Sun = evening star
    # Negative difference = planet is west of Sun = morning star
    is_evening_star = bool(lon_diff > 0)

    return lon_diff, is_evening_star


def get_signed_elongation(tjd_ut: float, ipl: int, iflag: int = 0) -> float:
    """
    Calculate signed elongation of a planet from the Sun.

    Returns a signed value where:
        - Positive values indicate eastern elongation (evening star)
        - Negative values indicate western elongation (morning star)

    This is equivalent to calling get_elongation_from_sun() and returning
    only the first element.

    Args:
        tjd_ut: Julian Day in Universal Time (UT1)
        ipl: Planet/body ID (MERCURY, VENUS, MARS, etc.)
        iflag: Calculation flags (default: 0)

    Returns:
        Signed elongation in degrees:
            - Positive = eastern elongation (evening star)
            - Negative = western elongation (morning star)

    Example:
        >>> from libephemeris import get_signed_elongation, MERCURY
        >>> elong = get_signed_elongation(2451545.0, MERCURY)
        >>> print(f"Mercury elongation: {elong:+.1f}°")
    """
    elongation, _ = get_elongation_from_sun(tjd_ut, ipl, iflag)
    return elongation


def is_morning_star(tjd_ut: float, ipl: int, iflag: int = 0) -> bool:
    """
    Determine if a planet is a morning star (western elongation).

    A planet is a morning star when it is west of the Sun (negative elongation),
    meaning it rises before the Sun and is visible in the eastern sky before sunrise.

    Args:
        tjd_ut: Julian Day in Universal Time (UT1)
        ipl: Planet/body ID (MERCURY, VENUS, MARS, etc.)
        iflag: Calculation flags (default: 0)

    Returns:
        True if the planet is a morning star (western elongation),
        False if it is an evening star (eastern elongation)

    Note:
        For the Sun, Moon, or unsupported bodies, this returns False.

    Example:
        >>> from libephemeris import is_morning_star, VENUS
        >>> if is_morning_star(2451545.0, VENUS):
        ...     print("Venus is a morning star")
        ... else:
        ...     print("Venus is an evening star")
    """
    if ipl == SUN:
        return False  # Sun cannot be morning/evening star
    _, is_evening = get_elongation_from_sun(tjd_ut, ipl, iflag)
    return not is_evening


def is_evening_star(tjd_ut: float, ipl: int, iflag: int = 0) -> bool:
    """
    Determine if a planet is an evening star (eastern elongation).

    A planet is an evening star when it is east of the Sun (positive elongation),
    meaning it sets after the Sun and is visible in the western sky after sunset.

    Args:
        tjd_ut: Julian Day in Universal Time (UT1)
        ipl: Planet/body ID (MERCURY, VENUS, MARS, etc.)
        iflag: Calculation flags (default: 0)

    Returns:
        True if the planet is an evening star (eastern elongation),
        False if it is a morning star (western elongation)

    Note:
        For the Sun, Moon, or unsupported bodies, this returns False.

    Example:
        >>> from libephemeris import is_evening_star, VENUS
        >>> if is_evening_star(2451545.0, VENUS):
        ...     print("Venus is an evening star")
        ... else:
        ...     print("Venus is a morning star")
    """
    if ipl == SUN:
        return False  # Sun cannot be morning/evening star
    _, is_evening = get_elongation_from_sun(tjd_ut, ipl, iflag)
    return is_evening


def get_elongation_type(tjd_ut: float, ipl: int, iflag: int = 0) -> str:
    """
    Get the elongation type as a descriptive string.

    Args:
        tjd_ut: Julian Day in Universal Time (UT1)
        ipl: Planet/body ID (MERCURY, VENUS, MARS, etc.)
        iflag: Calculation flags (default: 0)

    Returns:
        One of: "eastern" (evening star), "western" (morning star), or "none" (Sun)

    Example:
        >>> from libephemeris import get_elongation_type, VENUS
        >>> elong_type = get_elongation_type(2451545.0, VENUS)
        >>> print(f"Venus has {elong_type} elongation")
    """
    if ipl == SUN:
        return "none"
    _, is_evening = get_elongation_from_sun(tjd_ut, ipl, iflag)
    return "eastern" if is_evening else "western"
