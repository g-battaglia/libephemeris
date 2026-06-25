# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""LibEphemeris -- high-precision astronomical ephemeris library.

Pure Python, API-compatible with the reference ephemeris, using NASA JPL
DE440/DE441 ephemerides via Skyfield and IAU 2006/2000A standards via pyerfa.

Calculation backend:
    The default calculation mode is ``"auto"``, which resolves data via
    LEB (bundled, auto-discovered, or auto-downloaded), Horizons API, or
    Skyfield -- no manual configuration required. For ~14x faster batch
    evaluation, configure a precomputed ``.leb`` binary ephemeris file with
    ``set_leb_file()`` or the ``LIBEPHEMERIS_LEB`` environment variable.

    The calculation mode (``set_calc_mode()`` / ``LIBEPHEMERIS_MODE``) controls
    backend selection:

    - ``"auto"`` (default): try LEB first, then Horizons API (if no local
      DE440), then Skyfield.
    - ``"skyfield"``: always use Skyfield/DE440.
    - ``"leb"``: require a valid LEB file (configured, auto-discovered, or
      auto-downloaded). Raises ``RuntimeError`` if none can be resolved.
      Bodies/flags not in the LEB file fall back to Skyfield.
    - ``"horizons"``: prefer NASA JPL Horizons API (requires internet).
      Bodies/flags not supported by Horizons fall back to Skyfield.
"""

from __future__ import annotations

# Auto-load .env file before any module reads environment variables.
# This must happen before importing logging_config (which reads
# LIBEPHEMERIS_LOG_LEVEL at module level).
from ._dotenv import load_dotenv as _load_dotenv

_load_dotenv()
del _load_dotenv

# Auto-load TOML config file (libephemeris-config.toml) after .env.
# TOML values have lower priority than env vars (which includes .env).
from ._config_toml import load_config as _load_config

_load_config()
del _load_config

from .constants import *
from .constants import FLG_SPEED, GREG_CAL, JUL_CAL
from .logging_config import (
    get_logger,
    set_log_level,
    disable_logging,
    enable_logging,
    format_file_size,
    LIBEPHEMERIS_LOG_LEVEL_ENV,
)
from .tracing import (
    start_tracing,
    get_trace_results,
)
from .exceptions import (
    # Base error class (reference API compatible)
    Error,
    # Category: Input validation errors
    InputValidationError,
    CoordinateError,
    InvalidBodyError,
    # Category: Data not found errors
    DataNotFoundError,
    UnknownBodyError,
    StarNotFoundError,
    SPKNotFoundError,
    SPKRequiredError,
    # Category: Calculation errors
    CalculationError,
    PolarCircleError,
    EphemerisRangeError,
    ConvergenceError,
    # Category: Configuration errors
    ConfigurationError,
    # Validation helpers
    validate_latitude,
    validate_longitude,
    validate_coordinates,
)
from .time_utils import (
    julday,
    revjul,
    deltat,
    deltat_ex,
    date_conversion as _date_conversion_calendars,
    day_of_week,
    utc_to_jd,
    jdet_to_utc,
    jdut1_to_utc,
    utc_time_zone,
    time_equ,
    lat_to_lmt,
    lmt_to_lat,
    sidtime,
    sidtime0,
    # TAI (International Atomic Time) functions
    TT_TAI_OFFSET_SECONDS,
    TT_TAI_OFFSET_DAYS,
    get_tai_utc_for_jd,
    utc_to_tai_jd,
    tai_jd_to_utc,
    tt_to_tai_jd,
    tai_to_tt_jd,
    tai_to_utc_jd,
    utc_jd_to_tai,
)
from .planets import (
    calc_ut,
    calc,
    calc_pctr,
    get_ayanamsa_ut,
    get_ayanamsa,
    get_ayanamsa_ex,
    get_ayanamsa_ex_ut,
    get_ayanamsa_name,
    set_sid_mode,
    nod_aps,
    nod_aps_ut,
    get_orbital_elements,
    get_orbital_elements_ut,
    orbit_max_min_true_distance,
    pheno,
    pheno_ut,
    get_planet_name,
    HeliocentricNodApsWarning,
    NutationFallbackWarning,
    get_nutation_model,
    # Elongation helper functions
    get_elongation_from_sun,
    get_signed_elongation,
    is_morning_star,
    is_evening_star,
    get_elongation_type,
)
from .houses import (
    houses,
    houses_armc,
    houses_armc_ex2,
    houses_ex,
    houses_ex2,
    house_name,
    house_pos,
    _house_pos_pythonic,
    _gauquelin_sector_pythonic,
    gauquelin_sector,
    # Polar latitude handling
    houses_with_fallback,
    houses_armc_with_fallback,
    get_polar_latitude_threshold,
    # Extreme latitude handling
    get_extreme_latitude_info,
    EXTREME_LATITUDE_THRESHOLD,
)
from .state import (
    set_topo,
    set_ephe_path,
    set_ephemeris_file,
    set_jpl_file,
    set_tid_acc,
    get_tid_acc,
    set_delta_t_userdef,
    get_delta_t_userdef,
    set_lapse_rate,
    get_lapse_rate,
    get_library_path,
    get_current_file_data,
    close,
    reset_session,
    # LEB binary ephemeris mode
    set_leb_file,
    get_leb_reader,
    set_calc_mode,
    get_calc_mode,
    set_auto_spk_download,
    get_auto_spk_download,
    set_spk_cache_dir,
    get_spk_cache_dir,
    set_spk_date_padding,
    get_spk_date_padding,
    # IERS Delta T configuration
    set_iers_delta_t_enabled,
    get_iers_delta_t_enabled,
    # Delta T model selection
    set_delta_t_model,
    get_delta_t_model,
    # Strict precision mode
    set_strict_precision,
    get_strict_precision,
    # Precision tier system
    PrecisionTier,
    TIERS,
    set_precision_tier,
    get_precision_tier,
    list_tiers,
    get_spk_date_range_for_tier,
    # Page cache management (containerised environments)
    release_data_cache,
)
from .iers_data import (
    # IERS data download functions
    download_iers_finals,
    download_leap_seconds,
    download_delta_t_data,
    load_iers_data,
    # IERS Delta T lookup functions
    get_observed_delta_t,
    get_observed_delta_t_data_range,
    is_observed_delta_t_available,
    get_delta_t_iers,
    get_ut1_utc,
    get_tai_utc,
    get_iers_data_range,
    is_iers_data_available,
    # IERS cache management
    clear_iers_cache,
    delete_iers_cache_files,
    get_iers_cache_info,
    # IERS configuration
    set_iers_cache_dir,
    get_iers_cache_dir,
    set_iers_auto_download,
    get_iers_auto_download,
)
from .crossing import (
    solcross_ut,
    solcross,
    mooncross_ut,
    mooncross,
    mooncross_node_ut,
    mooncross_node,
    cross_ut,
    helio_cross_ut,
    helio_cross,
    find_station_ut,
    next_retrograde_ut,
)
from .lunar import (
    # Meeus polynomial validity warnings/errors
    MeeusPolynomialWarning,
    MeeusRangeError,
    # Mean lunar functions
    calc_mean_lunar_node,
    calc_mean_lilith,
    calc_mean_lilith_with_latitude,
    # True lunar functions
    calc_true_lunar_node,
    calc_true_lilith,
    calc_interpolated_apogee,
    calc_interpolated_perigee,
    # Validity range constants
    MEEUS_OPTIMAL_CENTURIES,
    MEEUS_VALID_CENTURIES,
    MEEUS_MAX_CENTURIES,
)
from .eclipse import (
    BesselianElements,
    sol_eclipse_max_time,
    _sol_eclipse_when_glob_pythonic,
    sol_eclipse_when_glob,
    _sol_eclipse_when_loc_pythonic,
    sol_eclipse_when_loc,
    _sol_eclipse_where_pythonic,
    sol_eclipse_where,
    _sol_eclipse_how_pythonic,
    sol_eclipse_how,
    _sol_eclipse_how_details_pythonic,
    sol_eclipse_how_details,
    # Eclipse magnitude at location
    _sol_eclipse_magnitude_at_loc_pythonic,
    sol_eclipse_magnitude_at_loc,
    # Eclipse obscuration at location
    _sol_eclipse_obscuration_at_loc_pythonic,
    sol_eclipse_obscuration_at_loc,
    _lun_eclipse_when_pythonic,
    lun_eclipse_when,
    _lun_eclipse_when_loc_pythonic,
    lun_eclipse_when_loc,
    _lun_eclipse_how_pythonic,
    lun_eclipse_how,
    # Lunar eclipse umbral magnitude
    _lun_eclipse_umbral_magnitude_pythonic,
    lun_eclipse_umbral_magnitude,
    # Lunar eclipse penumbral magnitude
    _lun_eclipse_penumbral_magnitude_pythonic,
    lun_eclipse_penumbral_magnitude,
    # Lunar eclipse gamma parameter
    _lun_eclipse_gamma_pythonic,
    lun_eclipse_gamma,
    lun_occult_when_glob,
    _lun_occult_when_loc_pythonic,
    lun_occult_when_loc,
    _lun_occult_where_pythonic,
    lun_occult_where,
    # Planetary occultations
    planet_occult_when_glob,
    planet_occult_when_loc,
    rise_trans,
    rise_trans_true_hor,
    calc_besselian_x,
    calc_besselian_y,
    calc_besselian_d,
    calc_besselian_l1,
    calc_besselian_l2,
    calc_besselian_mu,
    calc_besselian_dx_dt,
    calc_besselian_dy_dt,
    calc_besselian_dd_dt,
    calc_besselian_dl1_dt,
    calc_besselian_dl2_dt,
    calc_besselian_dmu_dt,
    interpolate_besselian_elements,
    calc_eclipse_first_contact_c1,
    calc_eclipse_second_contact_c2,
    calc_eclipse_third_contact_c3,
    calc_eclipse_fourth_contact_c4,
    calc_lunar_eclipse_penumbral_first_contact_p1,
    calc_lunar_eclipse_penumbral_fourth_contact_p4,
    calc_lunar_eclipse_umbral_first_contact_u1,
    calc_lunar_eclipse_umbral_second_contact_u2,
    calc_lunar_eclipse_umbral_third_contact_u3,
    calc_lunar_eclipse_umbral_fourth_contact_u4,
    # Eclipse duration functions
    calc_solar_eclipse_duration,
    calc_lunar_eclipse_total_duration,
    calc_lunar_eclipse_umbral_duration,
    # Eclipse path width calculation
    calc_eclipse_path_width,
    # Eclipse central line coordinates
    calc_eclipse_central_line,
    # Eclipse northern limit coordinates
    calc_eclipse_northern_limit,
    # Eclipse southern limit coordinates
    calc_eclipse_southern_limit,
    # Saros series calculation
    get_saros_number,
    SAROS_CYCLE_DAYS,
    # Inex series calculation
    get_inex_number,
    INEX_CYCLE_DAYS,
)
from .heliacal import (
    _heliacal_ut_pythonic,
    heliacal_ut,
    _heliacal_pheno_ut_pythonic,
    heliacal_pheno_ut,
    vis_limit_mag,
    is_inner_planet,
    is_fixed_star,
    INNER_PLANETS,
)
from .extinction import (
    # Extinction calculation functions
    calc_airmass,
    calc_extinction_coefficient,
    calc_extinction_magnitude,
    calc_simple_extinction,
    apparent_magnitude_with_extinction,
    get_extinction_for_heliacal,
    # Individual component functions
    calc_rayleigh_coefficient,
    calc_aerosol_coefficient,
    calc_ozone_coefficient,
    calc_water_vapor_coefficient,
    # Data class for extinction components
    ExtinctionCoefficients,
    # Constants
    WAVELENGTH_U,
    WAVELENGTH_B,
    WAVELENGTH_V,
    WAVELENGTH_R,
    WAVELENGTH_I,
    # Twilight sky brightness model
    TwilightSkyBrightness,
    get_twilight_phase,
    calc_twilight_sky_brightness,
    calc_twilight_brightness_simple,
    calc_limiting_magnitude_twilight,
    # Twilight phase constants
    TWILIGHT_CIVIL_START,
    TWILIGHT_CIVIL_END,
    TWILIGHT_NAUTICAL_END,
    TWILIGHT_ASTRONOMICAL_END,
    DARK_SKY_BRIGHTNESS_V,
    # Schaefer visibility threshold model
    VisibilityResult,
    calc_eye_adaptation_state,
    calc_contrast_threshold,
    calc_visibility_threshold,
    is_object_visible,
    calc_limiting_magnitude_for_sky,
    OBSERVER_SKILL_INEXPERIENCED,
    OBSERVER_SKILL_AVERAGE,
    OBSERVER_SKILL_EXPERIENCED,
    OBSERVER_SKILL_EXPERT,
    EXPERIENCE_FACTORS,
)
from .utils import (
    degnorm,
    radnorm,
    difdeg2n,
    difdegn,
    difrad2n,
    difcs2n,
    difcsn,
    csnorm,
    csroundsec,
    cs2degstr,
    cs2lonlatstr,
    cs2timestr,
    deg_midp,
    rad_midp,
    d2l,
    split_deg,
    calc_angles,
    cotrans,
    cotrans_sp,
    azalt,
    azalt_rev,
    refrac,
    refrac_extended,
    ECL2HOR,
    EQU2HOR,
    HOR2ECL,
    HOR2EQU,
    TRUE_TO_APP,
    APP_TO_TRUE,
    SPLIT_DEG_ROUND_SEC,
    SPLIT_DEG_ROUND_MIN,
    SPLIT_DEG_ROUND_DEG,
    SPLIT_DEG_ZODIACAL,
    SPLIT_DEG_NAKSHATRA,
    SPLIT_DEG_KEEP_SIGN,
    SPLIT_DEG_KEEP_DEG,
)
from .fixed_stars import (
    fixstar_ut,
    fixstar,
    fixstar2,
    fixstar2_ut,
    batch_fixstars_ut,
    fixstar_mag,
    fixstar2_mag,
    list_fixed_stars,
    StarCatalogEntry,
    propagate_proper_motion,
)
from .context import EphemerisContext  # NEW: Thread-safe context API
from .spk import (  # SPK kernel support for high-precision minor body calculations
    download_spk,
    register_spk_body,
    unregister_spk_body,
    get_spk_body_info,
    list_spk_bodies,
    get_spk_coverage,
    download_and_register_spk,
)
from . import spk_auto  # Automatic SPK download and caching
from .spk_auto import (
    discover_local_spks,
    ensure_all_ephemerides,
)
from .download import download_for_tier, download_leb_for_tier
from .minor_bodies import (  # Generic asteroid lookup by number
    calc_asteroid_by_number,
    fetch_orbital_elements_from_sbdb,
    clear_asteroid_elements_cache,
    get_asteroid_number,
    clear_asteroid_name_cache,
    # Automatic SPK download for major asteroids
    auto_download_asteroid_spk,
    is_spk_available_for_body,
    ensure_major_asteroid_spk,
    is_major_asteroid,
    get_major_asteroid_info,
    list_major_asteroids,
    MAJOR_ASTEROID_SPK_INFO,
)
from .hypothetical import (  # Hamburg School Uranian planets
    calc_cupido,
    calc_hades,
    calc_zeus,
    calc_kronos,
    calc_apollon,
    calc_admetos,
    calc_vulkanus,
    calc_poseidon,
    calc_transpluto,  # Transpluto (Isis) - hypothetical trans-Plutonian planet
    calc_vulcan,  # Vulcan - hypothetical intramercurial planet
    calc_waldemath,  # Waldemath Moon (hypothetical secondary satellite)
    calc_proserpina,  # Proserpina - hypothetical trans-Plutonian planet
    calc_planet_x_pickering,  # Pickering's Planet O/X prediction (1919)
    calc_white_moon_position,  # White Moon (Selena) - opposite to Black Moon Lilith
    calc_uranian_planet,
    calc_hypothetical_position,
    URANIAN_KEPLERIAN_ELEMENTS,
    TRANSPLUTO_KEPLERIAN_ELEMENTS,
    VULCAN_ELEMENTS,  # Vulcan orbital elements
    WALDEMATH_ELEMENTS,  # Waldemath Moon orbital elements
    PICKERING_PLANET_X_ELEMENTS,  # Pickering's Planet O/X orbital elements
    # Orbital elements parser for custom hypothetical bodies (canonical names)
    parse_orbital_elements,
    OrbitalElements,
    get_orbital_body_by_name,
    calc_orbital_position,
    get_bundled_fictitious_orbits_path,
    load_bundled_fictitious_orbits,
    # Backward-compatible aliases (legacy SE-derived names)
    parse_seorbel,
    SeorbelElements,
    get_seorbel_body_by_name,
    calc_seorbel_position,
    get_bundled_seorbel_path,
    load_bundled_seorbel,
    TPolynomial,
    CUPIDO as CUPIDO_HYPO,  # Alias to avoid conflict with constants.py
    HADES as HADES_HYPO,  # Alias to avoid conflict with constants.py
    ZEUS as ZEUS_HYPO,  # Alias to avoid conflict with constants.py
    KRONOS as KRONOS_HYPO,  # Alias to avoid conflict with constants.py
    APOLLON as APOLLON_HYPO,  # Alias to avoid conflict with constants.py
    ADMETOS as ADMETOS_HYPO,  # Alias to avoid conflict with constants.py
    VULKANUS as VULKANUS_HYPO,  # Alias to avoid conflict with constants.py
    POSEIDON as POSEIDON_HYPO,  # Alias to avoid conflict with constants.py
    ISIS as ISIS_HYPO,  # Alias to avoid conflict with constants.py
    TRANSPLUTO as TRANSPLUTO_HYPO,  # Alias to avoid conflict with constants.py
    VULCAN as VULCAN_HYPO,  # Alias to avoid conflict with constants.py
    WALDEMATH as WALDEMATH_HYPO,  # Alias to avoid conflict with constants.py
    WHITE_MOON as WHITE_MOON_HYPO,  # Alias to avoid conflict with constants.py
    PROSERPINA as PROSERPINA_HYPO,  # Alias to avoid conflict with constants.py
    PLANET_X_PICKERING as PLANET_X_PICKERING_HYPO,  # Alias to avoid conflict with constants.py
)

# REBOUND/ASSIST n-body integration (optional dependency)
# Lazy import to avoid requiring rebound/assist at module load time
from .rebound_integration import (
    ReboundIntegrator,
    AssistEphemConfig,
    PropagationResult,
    check_rebound_available,
    check_assist_available,
    get_rebound_version,
    get_assist_version,
    propagate_orbit_rebound,
    propagate_orbit_assist,
    propagate_trajectory,
    compare_with_keplerian,
)

from .planetary_moons import (  # Planetary moons (Galilean moons, Titan, etc.)
    register_moon_spk,
    unregister_moon_spk,
    list_registered_moons,
    get_moon_name,
    is_planetary_moon,
    get_moon_coverage,
    calc_moon_position,
    close_moon_kernels,
    # NAIF IDs for planetary moons
    NAIF_IO,
    NAIF_EUROPA,
    NAIF_GANYMEDE,
    NAIF_CALLISTO,
    NAIF_MIMAS,
    NAIF_ENCELADUS,
    NAIF_TETHYS,
    NAIF_DIONE,
    NAIF_RHEA,
    NAIF_TITAN,
    NAIF_HYPERION,
    NAIF_IAPETUS,
    NAIF_MIRANDA,
    NAIF_ARIEL,
    NAIF_UMBRIEL,
    NAIF_TITANIA,
    NAIF_OBERON,
    NAIF_TRITON,
    NAIF_PHOBOS,
    NAIF_DEIMOS,
    NAIF_CHARON,
    # Moon ID to NAIF mapping
    MOON_NAIF_MAP,
    MOON_NAMES,
    MOON_PARENT_MAP,
)

# =============================================================================
# date_conversion wrapper (reference-API return shape)
# =============================================================================


def date_conversion(
    year: int,
    month: int,
    day: int,
    hour: float = 12.0,
    cal: "str | bytes" = b"g",
) -> "tuple[bool, float, tuple[int, int, int, float]]":
    """Convert and validate a calendar date, returning Julian Day number.

    Wrapper matching the reference ``date_conversion`` return convention:
    ``(valid, jd, (year, month, day, hour))``.

    Args:
        year: Calendar year.
        month: Month (1-12).
        day: Day of month (1-31).
        hour: Decimal hour (0.0-23.999...).
        cal: ``'g'`` / ``b'g'`` for Gregorian, ``'j'`` / ``b'j'`` for Julian.

    Returns:
        Tuple of ``(valid, jd, (year, month, day, hour))`` where *valid* is
        ``True`` when the date exists in the requested calendar, *jd* is the
        Julian Day number, and the inner tuple holds the (possibly normalised)
        date components.
    """
    if isinstance(cal, bytes):
        cal = cal.decode("ascii")
    cal_char = cal.lower()
    if cal_char not in ("j", "g"):
        raise ValueError(f"calendar must be 'j' or 'g', got: {cal!r}")

    cal_flag = JUL_CAL if cal_char == "j" else GREG_CAL
    jd = julday(year, month, day, hour, cal_flag)
    # Round-trip to check validity: convert JD back to calendar date
    y2, m2, d2, h2 = revjul(jd, cal_flag)
    valid = y2 == year and m2 == month and d2 == day
    return (valid, jd, (y2, m2, d2, h2))


# Helper for Arabic parts
from .arabic_parts import calc_all_arabic_parts

# .env file loader (public API for manual reloading)
from ._dotenv import load_dotenv

# Extended astrology helpers submodule
from . import contrib

__version__ = "3.0.0rc1"
version = __version__
__author__ = "Giacomo Battaglia"
__license__ = "AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial"

# Note: the original Pythonic variants of sol_eclipse_*, lun_eclipse_*,
# lun_occult_when_loc, heliacal_ut, heliacal_pheno_ut, and gauquelin_sector
# are imported above as _*_pythonic from their source modules. They are
# kept as private helpers, NOT re-exported as public bare names. The
# public bare names are the reference-API-compatible wrappers.

__all__ = [
    # .env file loader
    "load_dotenv",
    # Logging configuration
    "get_logger",
    "set_log_level",
    "disable_logging",
    "enable_logging",
    "format_file_size",
    # Computation tracing
    "start_tracing",
    "get_trace_results",
    # Exceptions - Base
    "Error",
    # Exceptions - Input Validation Category
    "InputValidationError",
    "CoordinateError",
    "InvalidBodyError",
    # Exceptions - Data Not Found Category
    "DataNotFoundError",
    "UnknownBodyError",
    "StarNotFoundError",
    "SPKNotFoundError",
    "SPKRequiredError",
    # Exceptions - Calculation Category
    "CalculationError",
    "PolarCircleError",
    "EphemerisRangeError",
    "ConvergenceError",
    # Exceptions - Configuration Category
    "ConfigurationError",
    # Lunar calculation warnings/errors
    "MeeusPolynomialWarning",
    "MeeusRangeError",
    # Lunar node and apogee functions
    "calc_mean_lunar_node",
    "calc_mean_lilith",
    "calc_mean_lilith_with_latitude",
    "calc_true_lunar_node",
    "calc_true_lilith",
    "calc_interpolated_apogee",
    "calc_interpolated_perigee",
    # Meeus validity range constants
    "MEEUS_OPTIMAL_CENTURIES",
    "MEEUS_VALID_CENTURIES",
    "MEEUS_MAX_CENTURIES",
    # Coordinate validation
    "validate_latitude",
    "validate_longitude",
    "validate_coordinates",
    # Thread-safe Context API
    "EphemerisContext",
    # Time functions (both swe_ and non-prefixed aliases)
    "julday",
    "revjul",
    "deltat",
    "deltat_ex",
    "date_conversion",
    "day_of_week",
    "utc_to_jd",
    "jdet_to_utc",
    "jdut1_to_utc",
    "utc_time_zone",
    "time_equ",
    "lat_to_lmt",
    "lmt_to_lat",
    "sidtime",
    "sidtime0",
    # TAI (International Atomic Time) functions
    "TT_TAI_OFFSET_SECONDS",
    "TT_TAI_OFFSET_DAYS",
    "get_tai_utc_for_jd",
    "utc_to_tai_jd",
    "tai_jd_to_utc",
    "tt_to_tai_jd",
    "tai_to_tt_jd",
    "tai_to_utc_jd",
    "utc_jd_to_tai",
    # Planet calculation
    "calc_ut",
    "calc",
    "calc_pctr",
    "nod_aps",
    "nod_aps_ut",
    "HeliocentricNodApsWarning",
    "get_orbital_elements",
    "get_orbital_elements_ut",
    "orbit_max_min_true_distance",
    # Planetary phenomena
    "pheno",
    "pheno_ut",
    # Elongation helper functions
    "get_elongation_from_sun",
    "get_signed_elongation",
    "is_morning_star",
    "is_evening_star",
    "get_elongation_type",
    # Houses
    "houses",
    "houses_armc",
    "houses_armc_ex2",
    "houses_ex",
    "houses_ex2",
    "house_name",
    "house_pos",
    "gauquelin_sector",
    # Ayanamsa (sidereal)
    "set_sid_mode",
    "get_ayanamsa_ut",
    "get_ayanamsa",
    "get_ayanamsa_ex",
    "get_ayanamsa_ex_ut",
    "get_ayanamsa_name",
    # Observer location
    "set_topo",
    "set_ephe_path",
    "set_ephemeris_file",
    "set_jpl_file",
    # Tidal acceleration
    "set_tid_acc",
    "get_tid_acc",
    # User-defined Delta T
    "set_delta_t_userdef",
    "get_delta_t_userdef",
    # Lapse rate for refraction
    "set_lapse_rate",
    "get_lapse_rate",
    # Library path
    "get_library_path",
    # Current file data
    "get_current_file_data",
    # Close and cleanup
    "close",
    # Crossings
    "solcross_ut",
    "solcross",
    "mooncross_ut",
    "mooncross",
    "mooncross_node_ut",
    "mooncross_node",
    "cross_ut",
    "helio_cross_ut",
    "helio_cross",
    "find_station_ut",
    "next_retrograde_ut",
    # Eclipses
    "sol_eclipse_max_time",
    "sol_eclipse_when_glob",
    "sol_eclipse_when_loc",
    "sol_eclipse_where",
    "sol_eclipse_how",
    "sol_eclipse_how_details",
    # Eclipse magnitude at location
    "sol_eclipse_magnitude_at_loc",
    # Eclipse obscuration at location
    "sol_eclipse_obscuration_at_loc",
    "lun_eclipse_when",
    "lun_eclipse_when_loc",
    "lun_eclipse_how",
    # Lunar eclipse umbral magnitude
    "lun_eclipse_umbral_magnitude",
    # Lunar eclipse penumbral magnitude
    "lun_eclipse_penumbral_magnitude",
    "lun_occult_when_glob",
    "lun_occult_when_loc",
    "lun_occult_where",
    # Besselian elements
    "calc_besselian_x",
    "calc_besselian_y",
    "calc_besselian_d",
    "calc_besselian_l1",
    "calc_besselian_l2",
    "calc_besselian_mu",
    # Besselian element time derivatives
    "calc_besselian_dx_dt",
    "calc_besselian_dy_dt",
    "calc_besselian_dd_dt",
    "calc_besselian_dl1_dt",
    "calc_besselian_dl2_dt",
    "calc_besselian_dmu_dt",
    # Solar eclipse contact points
    "calc_eclipse_first_contact_c1",
    "calc_eclipse_second_contact_c2",
    "calc_eclipse_third_contact_c3",
    "calc_eclipse_fourth_contact_c4",
    # Lunar eclipse penumbral contact points
    "calc_lunar_eclipse_penumbral_first_contact_p1",
    "calc_lunar_eclipse_penumbral_fourth_contact_p4",
    # Lunar eclipse umbral contact points
    "calc_lunar_eclipse_umbral_first_contact_u1",
    "calc_lunar_eclipse_umbral_second_contact_u2",
    "calc_lunar_eclipse_umbral_third_contact_u3",
    "calc_lunar_eclipse_umbral_fourth_contact_u4",
    # Eclipse duration functions
    "calc_solar_eclipse_duration",
    "calc_lunar_eclipse_total_duration",
    "calc_lunar_eclipse_umbral_duration",
    # Eclipse path width calculation
    "calc_eclipse_path_width",
    # Eclipse central line coordinates
    "calc_eclipse_central_line",
    # Saros series calculation
    "get_saros_number",
    "SAROS_CYCLE_DAYS",
    # Inex series calculation
    "get_inex_number",
    "INEX_CYCLE_DAYS",
    # Rise/Set/Transit
    "rise_trans",
    "rise_trans_true_hor",
    # Heliacal events
    "heliacal_ut",
    "heliacal_pheno_ut",
    "vis_limit_mag",
    "is_inner_planet",
    "is_fixed_star",
    "INNER_PLANETS",
    # Atmospheric extinction
    "calc_airmass",
    "calc_extinction_coefficient",
    "calc_extinction_magnitude",
    "calc_simple_extinction",
    "apparent_magnitude_with_extinction",
    "get_extinction_for_heliacal",
    "calc_rayleigh_coefficient",
    "calc_aerosol_coefficient",
    "calc_ozone_coefficient",
    "calc_water_vapor_coefficient",
    "ExtinctionCoefficients",
    "WAVELENGTH_U",
    "WAVELENGTH_B",
    "WAVELENGTH_V",
    "WAVELENGTH_R",
    "WAVELENGTH_I",
    # Twilight sky brightness model
    "TwilightSkyBrightness",
    "get_twilight_phase",
    "calc_twilight_sky_brightness",
    "calc_twilight_brightness_simple",
    "calc_limiting_magnitude_twilight",
    "TWILIGHT_CIVIL_START",
    "TWILIGHT_CIVIL_END",
    "TWILIGHT_NAUTICAL_END",
    "TWILIGHT_ASTRONOMICAL_END",
    "DARK_SKY_BRIGHTNESS_V",
    # Schaefer visibility threshold model
    "VisibilityResult",
    "calc_eye_adaptation_state",
    "calc_contrast_threshold",
    "calc_visibility_threshold",
    "is_object_visible",
    "calc_limiting_magnitude_for_sky",
    "OBSERVER_SKILL_INEXPERIENCED",
    "OBSERVER_SKILL_AVERAGE",
    "OBSERVER_SKILL_EXPERIENCED",
    "OBSERVER_SKILL_EXPERT",
    "EXPERIENCE_FACTORS",
    # Utilities
    "degnorm",
    "radnorm",
    "difdeg2n",
    "difdegn",
    "difrad2n",
    "difcs2n",
    "difcsn",
    "csnorm",
    "csroundsec",
    "cs2degstr",
    "cs2lonlatstr",
    "cs2timestr",
    "deg_midp",
    "rad_midp",
    "d2l",
    "split_deg",
    "calc_angles",
    "cotrans",
    "cotrans_sp",
    "azalt",
    "azalt_rev",
    "refrac",
    "refrac_extended",
    # swe_ prefixed utility aliases (reference-API compatibility)
    "get_planet_name",
    "version",
    "ECL2HOR",
    "EQU2HOR",
    "HOR2ECL",
    "HOR2EQU",
    "TRUE_TO_APP",
    "APP_TO_TRUE",
    "SPLIT_DEG_ROUND_SEC",
    "SPLIT_DEG_ROUND_MIN",
    "SPLIT_DEG_ROUND_DEG",
    "SPLIT_DEG_ZODIACAL",
    "SPLIT_DEG_NAKSHATRA",
    "SPLIT_DEG_KEEP_SIGN",
    "SPLIT_DEG_KEEP_DEG",
    # Helpers
    "calc_all_arabic_parts",
    # Fixed Stars
    "fixstar_ut",
    "fixstar",
    "fixstar2",
    "fixstar2_ut",
    "batch_fixstars_ut",
    "list_fixed_stars",
    "StarCatalogEntry",
    "fixstar_mag",
    "fixstar2_mag",
    # Proper motion propagation
    "propagate_proper_motion",
    # SPK kernel support
    "download_spk",
    "register_spk_body",
    "unregister_spk_body",
    "get_spk_body_info",
    "list_spk_bodies",
    "get_spk_coverage",
    "download_and_register_spk",
    # Automatic SPK download module
    "spk_auto",
    # Tier-aware download
    "download_for_tier",
    "download_leb_for_tier",
    # Auto SPK download configuration
    "set_auto_spk_download",
    "get_auto_spk_download",
    # SPK cache directory configuration
    "set_spk_cache_dir",
    "get_spk_cache_dir",
    # SPK date padding configuration
    "set_spk_date_padding",
    "get_spk_date_padding",
    # SPK body name mapping
    "SPK_BODY_NAME_MAP",
    "SPK_AUTO_DOWNLOAD_BLOCKED",
    "is_spk_auto_download_blocked",
    "get_horizons_id",
    "get_naif_id_from_ipl",
    "get_spk_body_info_from_map",
    # Required SPK bodies for high-precision calculations
    "REQUIRED_SPK_BODIES",
    # IERS Delta T configuration
    "set_iers_delta_t_enabled",
    "get_iers_delta_t_enabled",
    # Delta T model selection
    "set_delta_t_model",
    "get_delta_t_model",
    # Strict precision mode
    "set_strict_precision",
    "get_strict_precision",
    # Precision tier system
    "PrecisionTier",
    "TIERS",
    "set_precision_tier",
    "get_precision_tier",
    "list_tiers",
    "get_spk_date_range_for_tier",
    # Page cache management (containerised environments)
    "release_data_cache",
    # SPK discovery and ensure
    "discover_local_spks",
    "ensure_all_ephemerides",
    # IERS data download functions
    "download_iers_finals",
    "download_leap_seconds",
    "download_delta_t_data",
    "load_iers_data",
    # IERS Delta T lookup functions
    "get_observed_delta_t",
    "get_observed_delta_t_data_range",
    "is_observed_delta_t_available",
    "get_delta_t_iers",
    "get_ut1_utc",
    "get_tai_utc",
    "get_iers_data_range",
    "is_iers_data_available",
    # IERS cache management
    "clear_iers_cache",
    "delete_iers_cache_files",
    "get_iers_cache_info",
    # IERS cache directory configuration
    "set_iers_cache_dir",
    "get_iers_cache_dir",
    # IERS auto-download configuration
    "set_iers_auto_download",
    "get_iers_auto_download",
    # Asteroid name lookup
    "get_asteroid_number",
    "clear_asteroid_name_cache",
    "calc_asteroid_by_number",
    "fetch_orbital_elements_from_sbdb",
    "clear_asteroid_elements_cache",
    # Automatic SPK download for major asteroids
    "auto_download_asteroid_spk",
    "is_spk_available_for_body",
    "ensure_major_asteroid_spk",
    "is_major_asteroid",
    "get_major_asteroid_info",
    "list_major_asteroids",
    "MAJOR_ASTEROID_SPK_INFO",
    # Hypothetical bodies (Hamburg School Uranian planets)
    "calc_cupido",
    "calc_hades",
    "calc_zeus",
    "calc_kronos",
    "calc_apollon",
    "calc_admetos",
    "calc_vulkanus",
    "calc_poseidon",
    "calc_transpluto",  # Transpluto (Isis) - hypothetical trans-Plutonian planet
    "calc_vulcan",  # Vulcan - hypothetical intramercurial planet
    "calc_waldemath",  # Waldemath Moon - hypothetical second moon of Earth
    "calc_planet_x_pickering",  # Pickering's Planet O/X prediction (1919)
    "calc_white_moon_position",  # White Moon (Selena) - opposite to Black Moon Lilith
    "calc_uranian_planet",
    "calc_hypothetical_position",
    "URANIAN_KEPLERIAN_ELEMENTS",
    "TRANSPLUTO_KEPLERIAN_ELEMENTS",
    "VULCAN_ELEMENTS",
    "WALDEMATH_ELEMENTS",  # Waldemath Moon orbital elements
    "PICKERING_PLANET_X_ELEMENTS",  # Pickering's Planet O/X orbital elements
    # Orbital elements parser for custom hypothetical bodies (canonical names)
    "parse_orbital_elements",
    "OrbitalElements",
    "get_orbital_body_by_name",
    "calc_orbital_position",
    "get_bundled_fictitious_orbits_path",
    "load_bundled_fictitious_orbits",
    # Backward-compatible aliases (legacy SE-derived names)
    "parse_seorbel",
    "SeorbelElements",
    "get_seorbel_body_by_name",
    "calc_seorbel_position",
    "get_bundled_seorbel_path",
    "load_bundled_seorbel",
    "TPolynomial",
    # Planetary moons (Galilean moons, Titan, etc.)
    "register_moon_spk",
    "unregister_moon_spk",
    "list_registered_moons",
    "get_moon_name",
    "is_planetary_moon",
    "get_moon_coverage",
    "calc_moon_position",
    "close_moon_kernels",
    "MOON_NAIF_MAP",
    "MOON_NAMES",
    "MOON_PARENT_MAP",
    # REBOUND/ASSIST n-body integration
    "ReboundIntegrator",
    "AssistEphemConfig",
    "PropagationResult",
    "check_rebound_available",
    "check_assist_available",
    "get_rebound_version",
    "get_assist_version",
    "propagate_orbit_rebound",
    "propagate_orbit_assist",
    "propagate_trajectory",
    "compare_with_keplerian",
    # Ephemeris selection flags (from constants.py)
    "FLG_JPLEPH",
    "FLG_SWIEPH",
    "FLG_MOSEPH",
    # Bare-name constants (reference-API-compatible, without SE_/SEFLG_ prefix)
    "ACRONYCHAL_RISING",
    "ACRONYCHAL_SETTING",
    "ADMETOS",
    "APOLLON",
    "ARMC",
    "ASC",
    "ASTNAMFILE",
    "AST_OFFSET",
    "AUNIT_TO_KM",
    "AUNIT_TO_LIGHTYEAR",
    "AUNIT_TO_PARSEC",
    "BIT_ASTRO_TWILIGHT",
    "BIT_CIVIL_TWILIGHT",
    "BIT_DISC_BOTTOM",
    "BIT_DISC_CENTER",
    "BIT_FIXED_DISC_SIZE",
    "BIT_FORCE_SLOW_METHOD",
    "BIT_HINDU_RISING",
    "BIT_NAUTIC_TWILIGHT",
    "BIT_NO_REFRACTION",
    "CALC_ITRANSIT",
    "CALC_MTRANSIT",
    "CALC_RISE",
    "CALC_SET",
    "CERES",
    "CHIRON",
    "COASC1",
    "COASC2",
    "COMET_OFFSET",
    "COSMICAL_SETTING",
    "CUPIDO",
    "DE_NUMBER",
    "DELTAT_AUTOMATIC",
    "EARTH",
    "ECL_1ST_VISIBLE",
    "ECL_2ND_VISIBLE",
    "ECL_3RD_VISIBLE",
    "ECL_4TH_VISIBLE",
    "ECL_ALLTYPES_LUNAR",
    "ECL_ALLTYPES_SOLAR",
    "ECL_ANNULAR",
    "ECL_ANNULAR_TOTAL",
    "ECL_CENTRAL",
    "ECL_HYBRID",
    "ECL_MAX_VISIBLE",
    "ECL_NONCENTRAL",
    "ECL_NUT",
    "ECL_OCC_BEG_DAYLIGHT",
    "ECL_OCC_END_DAYLIGHT",
    "ECL_ONE_TRY",
    "ECL_PARTBEG_VISIBLE",
    "ECL_PARTEND_VISIBLE",
    "ECL_PARTIAL",
    "ECL_PENUMBBEG_VISIBLE",
    "ECL_PENUMBEND_VISIBLE",
    "ECL_PENUMBRAL",
    "ECL_TOTAL",
    "ECL_TOTBEG_VISIBLE",
    "ECL_TOTEND_VISIBLE",
    "ECL_VISIBLE",
    "EPHE_PATH",
    "EQUASC",
    "EVENING_FIRST",
    "EVENING_LAST",
    "FICTFILE",
    "FICT_MAX",
    "FICT_OFFSET",
    "FICT_OFFSET_1",
    "FIXSTAR",
    "FLG_ASTROMETRIC",
    "FLG_BARYCTR",
    "FLG_CENTER_BODY",
    "FLG_DEFAULTEPH",
    "FLG_DPSIDEPS_1980",
    "FLG_EQUATORIAL",
    "FLG_HELCTR",
    "FLG_ICRS",
    "FLG_J2000",
    "FLG_JPLHOR",
    "FLG_JPLHOR_APPROX",
    "FLG_NOABERR",
    "FLG_NOGDEFL",
    "FLG_NONUT",
    "FLG_ORBEL_AA",
    "FLG_RADIANS",
    "FLG_SIDEREAL",
    "FLG_SPEED",
    "FLG_SPEED3",
    "FLG_TEST_PLMOON",
    "FLG_TOPOCTR",
    "FLG_TROPICAL",
    "FLG_TRUEPOS",
    "FLG_XYZ",
    "FNAME_DE200",
    "FNAME_DE403",
    "FNAME_DE404",
    "FNAME_DE405",
    "FNAME_DE406",
    "FNAME_DFT",
    "FNAME_DFT2",
    "GREG_CAL",
    "HADES",
    "HARRINGTON",
    "HELFLAG_AV",
    "HELFLAG_AVKIND",
    "HELFLAG_AVKIND_MIN7",
    "HELFLAG_AVKIND_MIN9",
    "HELFLAG_AVKIND_PTO",
    "HELFLAG_AVKIND_VR",
    "HELFLAG_HIGH_PRECISION",
    "HELFLAG_LONG_SEARCH",
    "HELFLAG_NO_DETAILS",
    "HELFLAG_OPTICAL_PARAMS",
    "HELFLAG_SEARCH_1_PERIOD",
    "HELFLAG_VISLIM_DARK",
    "HELFLAG_VISLIM_NOMOON",
    "HELFLAG_VISLIM_PHOTOPIC",
    "HELFLAG_VISLIM_SCOTOPIC",
    "HELIACAL_RISING",
    "HELIACAL_SETTING",
    "INTP_APOG",
    "INTP_PERG",
    "ISIS",
    "JUL_CAL",
    "JUNO",
    "JUPITER",
    "KRONOS",
    "MARS",
    "MAX_STNAME",
    "MC",
    "MEAN_APOG",
    "MEAN_NODE",
    "MERCURY",
    "MIXEDOPIC_FLAG",
    "MODEL_BIAS",
    "MODEL_DELTAT",
    "MODEL_JPLHORA_MODE",
    "MODEL_JPLHOR_MODE",
    "MODEL_NUT",
    "MODEL_PREC_LONGTERM",
    "MODEL_PREC_SHORTTERM",
    "MODEL_SIDT",
    "MOD_BIAS_DEFAULT",
    "MOD_BIAS_IAU2000",
    "MOD_BIAS_IAU2006",
    "MOD_BIAS_NONE",
    "MOD_DELTAT_DEFAULT",
    "MOD_DELTAT_ESPENAK_MEEUS_2006",
    "MOD_DELTAT_STEPHENSON_1997",
    "MOD_DELTAT_STEPHENSON_ETC_2016",
    "MOD_DELTAT_STEPHENSON_MORRISON_1984",
    "MOD_DELTAT_STEPHENSON_MORRISON_2004",
    "MOD_JPLHORA_1",
    "MOD_JPLHORA_2",
    "MOD_JPLHORA_3",
    "MOD_JPLHORA_DEFAULT",
    "MOD_JPLHOR_DEFAULT",
    "MOD_JPLHOR_LONG_AGREEMENT",
    "MOD_NBIAS",
    "MOD_NDELTAT",
    "MOD_NJPLHOR",
    "MOD_NJPLHORA",
    "MOD_NNUT",
    "MOD_NPREC",
    "MOD_NUT_DEFAULT",
    "MOD_NUT_IAU_1980",
    "MOD_NUT_IAU_2000A",
    "MOD_NUT_IAU_2000B",
    "MOD_NUT_IAU_CORR_1987",
    "MOD_NUT_WOOLARD",
    "MOD_PREC_BRETAGNON_2003",
    "MOD_PREC_DEFAULT",
    "MOD_PREC_DEFAULT_SHORT",
    "MOD_PREC_IAU_1976",
    "MOD_PREC_IAU_2000",
    "MOD_PREC_IAU_2006",
    "MOD_PREC_LASKAR_1986",
    "MOD_PREC_NEWCOMB",
    "MOD_PREC_OWEN_1990",
    "MOD_PREC_SIMON_1994",
    "MOD_PREC_VONDRAK_2011",
    "MOD_PREC_WILLIAMS_1994",
    "MOD_PREC_WILL_EPS_LASK",
    "MOON",
    "MORNING_FIRST",
    "MORNING_LAST",
    "NALL_NAT_POINTS",
    "NASCMC",
    "NEPTUNE",
    "NEPTUNE_ADAMS",
    "NEPTUNE_LEVERRIER",
    "NFICT_ELEM",
    "NIBIRU",
    "NODBIT_FOPOINT",
    "NODBIT_MEAN",
    "NODBIT_OSCU",
    "NODBIT_OSCU_BAR",
    "NPLANETS",
    "NSE_MODELS",
    "NSIDM_PREDEF",
    "OSCU_APOG",
    "PALLAS",
    "PHOLUS",
    "PHOTOPIC_FLAG",
    "PLMOON_OFFSET",
    "PLUTO",
    "PLUTO_LOWELL",
    "PLUTO_PICKERING",
    "POLASC",
    "POSEIDON",
    "PROSERPINA",
    "SATURN",
    "SCOTOPIC_FLAG",
    "FNAME_DE431",
    "SIDBITS",
    "SIDBIT_ECL_DATE",
    "SIDBIT_ECL_T0",
    "SIDBIT_NO_PREC_OFFSET",
    "SIDBIT_PREC_ORIG",
    "SIDBIT_SSY_PLANE",
    "SIDBIT_USER_UT",
    "SIDM_ALDEBARAN_15TAU",
    "SIDM_ARYABHATA",
    "SIDM_ARYABHATA_522",
    "SIDM_ARYABHATA_MSUN",
    "SIDM_B1950",
    "SIDM_BABYL_BRITTON",
    "SIDM_BABYL_ETPSC",
    "SIDM_BABYL_HUBER",
    "SIDM_BABYL_KUGLER1",
    "SIDM_BABYL_KUGLER2",
    "SIDM_BABYL_KUGLER3",
    "SIDM_DELUCE",
    "SIDM_DJWHAL_KHUL",
    "SIDM_FAGAN_BRADLEY",
    "SIDM_GALALIGN_MARDYKS",
    "SIDM_GALCENT_0SAG",
    "SIDM_GALCENT_COCHRANE",
    "SIDM_GALCENT_MULA_WILHELM",
    "SIDM_GALCENT_RGILBRAND",
    "SIDM_GALEQU_FIORENZA",
    "SIDM_GALEQU_IAU1958",
    "SIDM_GALEQU_MULA",
    "SIDM_GALEQU_TRUE",
    "SIDM_HIPPARCHOS",
    "SIDM_J1900",
    "SIDM_J2000",
    "SIDM_JN_BHASIN",
    "SIDM_KRISHNAMURTI",
    "SIDM_KRISHNAMURTI_VP291",
    "SIDM_LAHIRI",
    "SIDM_LAHIRI_1940",
    "SIDM_LAHIRI_ICRC",
    "SIDM_LAHIRI_VP285",
    "SIDM_RAMAN",
    "SIDM_SASSANIAN",
    "SIDM_SS_CITRA",
    "SIDM_SS_REVATI",
    "SIDM_SURYASIDDHANTA",
    "SIDM_SURYASIDDHANTA_MSUN",
    "SIDM_TRUE_CITRA",
    "SIDM_TRUE_MULA",
    "SIDM_TRUE_PUSHYA",
    "SIDM_TRUE_REVATI",
    "SIDM_TRUE_SHEORAN",
    "SIDM_USER",
    "SIDM_USHASHASHI",
    "SIDM_VALENS_MOON",
    "SIDM_YUKTESHWAR",
    "SIMULATE_VICTORVB",
    "STARFILE",
    "STARFILE_OLD",
    "SUN",
    "TIDAL_26",
    "TIDAL_AUTOMATIC",
    "TIDAL_DE200",
    "TIDAL_DE403",
    "TIDAL_DE404",
    "TIDAL_DE405",
    "TIDAL_DE406",
    "TIDAL_DE421",
    "TIDAL_DE422",
    "TIDAL_DE430",
    "TIDAL_DE431",
    "TIDAL_DE441",
    "TIDAL_DEFAULT",
    "TIDAL_JPLEPH",
    "TIDAL_MOSEPH",
    "TIDAL_STEPHENSON_2016",
    "TIDAL_SWIEPH",
    "TJD_INVALID",
    "TRUE_NODE",
    "URANUS",
    "VARUNA",
    "VENUS",
    "VERTEX",
    "VESTA",
    "VULCAN",
    "VULKANUS",
    "WALDEMATH",
    "WHITE_MOON",
    "ZEUS",
]
