# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Exception classes for libephemeris.

This module provides a comprehensive exception hierarchy for ephemeris calculations.
The base Error class is reference API-compatible, and specialized exceptions provide
better categorization of error conditions.

Exception Hierarchy
-------------------
Error (base, reference API compatible)
├── InputValidationError (input data validation errors)
│   ├── CoordinateError (invalid latitude/longitude)
│   └── InvalidBodyError (body not valid for this operation)
├── DataNotFoundError (data lookup failures)
│   ├── UnknownBodyError (unknown celestial body ID)
│   ├── StarNotFoundError (fixed star not in catalog)
│   └── SPKNotFoundError (SPK kernel file not found)
├── CalculationError (calculation/algorithm failures)
│   ├── PolarCircleError (house calculation at polar latitudes)
│   ├── EphemerisRangeError (date outside ephemeris range)
│   └── ConvergenceError (numerical algorithm didn't converge)
└── ConfigurationError (missing configuration/state)

This hierarchy allows users to catch broad categories of errors or specific ones:

    # Catch any input validation error
    except InputValidationError:
        ...

    # Catch only coordinate errors
    except CoordinateError:
        ...

    # Catch all libephemeris errors
    except Error:
        ...

Provenance:
    Project-authored error taxonomy and input validation infrastructure. It
    contains no astronomical model, physical constant, or source-derived data.
    Compatibility affects only the public base exception and error behavior.
"""

from __future__ import annotations

import math


class Error(Exception):
    """Ephemeris calculation error.

    This exception is raised for ephemeris-related errors such as:
    - Ephemeris files not found
    - Unsupported planet or body ID
    - Dates out of range for available ephemeris data
    - Fixed star not found
    - Calculation failures

    This class is designed to be compatible with the reference API's Error type, allowing client code
    that catches it to work
    unchanged with libephemeris.Error.

    Example:
        >>> import libephemeris as ephem
        >>> try:
        ...     ephem.calc_ut(jd, invalid_planet_id)
        ... except ephem.Error as e:
        ...     print(f"Ephemeris error: {e}")
    """

    pass


# =============================================================================
# CATEGORY: INPUT VALIDATION ERRORS
# =============================================================================


class InputValidationError(Error):
    """Base class for input validation errors.

    This exception category covers all cases where user-provided input
    data is invalid, such as:
    - Coordinates outside valid ranges
    - Invalid body types for specific operations

    Catching this exception will catch all input validation errors.

    Example:
        >>> try:
        ...     ephem.set_topo(200.0, 91.0, 0.0)  # Invalid coordinates
        ... except ephem.InputValidationError as e:
        ...     print(f"Invalid input: {e}")
    """

    pass


class CoordinateError(InputValidationError):
    """Error raised when geographic coordinates are invalid.

    This exception is raised when latitude or longitude values are outside
    their valid ranges:
    - Latitude: must be in [-90, 90] degrees
    - Longitude: must be in [-180, 180] degrees

    Attributes:
        message: Human-readable error message
        coordinate_name: Name of the coordinate ("latitude" or "longitude")
        value: The invalid coordinate value
        min_value: Minimum allowed value for this coordinate
        max_value: Maximum allowed value for this coordinate

    Example:
        >>> import libephemeris as ephem
        >>> try:
        ...     ephem.set_topo(0.0, 91.0, 0.0)  # Invalid latitude
        ... except ephem.CoordinateError as e:
        ...     print(f"Invalid {e.coordinate_name}: {e.value}")
        ...     print(f"Valid range: [{e.min_value}, {e.max_value}]")
        Invalid latitude: 91.0
        Valid range: [-90, 90]

    See Also:
        set_topo: Sets observer location (validates coordinates)
        houses: House calculation (validates latitude)
    """

    def __init__(
        self,
        message: str,
        coordinate_name: str | None = None,
        value: float | None = None,
        min_value: float | None = None,
        max_value: float | None = None,
    ):
        super().__init__(message)
        self.message = message
        self.coordinate_name = coordinate_name
        self.value = value
        self.min_value = min_value
        self.max_value = max_value

    def __str__(self) -> str:
        return self.message

    def __repr__(self) -> str:
        return (
            f"CoordinateError({self.message!r}, "
            f"coordinate_name={self.coordinate_name!r}, value={self.value}, "
            f"min_value={self.min_value}, max_value={self.max_value})"
        )


class InvalidBodyError(InputValidationError):
    """Error raised when a celestial body is not valid for a specific operation.

    This exception is raised when a known, valid body ID is used in an
    operation that doesn't support that body type. For example:
    - Sun or Moon for heliacal rising calculations
    - Moon for retrograde station searches
    - Earth for geocentric calculations

    This is different from UnknownBodyError which is raised for completely
    unknown body IDs.

    Attributes:
        message: Human-readable error message
        body_id: The body ID that was rejected
        body_name: Human-readable name of the body (if known)
        operation: The operation that rejected the body

    Example:
        >>> import libephemeris as ephem
        >>> try:
        ...     ephem.heliacal_ut(jd, geo, atmo, obs, ephem.SUN, ...)
        ... except ephem.InvalidBodyError as e:
        ...     print(f"Cannot use {e.body_name} for {e.operation}")
        Cannot use Sun for heliacal rising

    See Also:
        UnknownBodyError: For completely unknown body IDs
    """

    def __init__(
        self,
        message: str,
        body_id: int | None = None,
        body_name: str | None = None,
        operation: str | None = None,
    ):
        super().__init__(message)
        self.message = message
        self.body_id = body_id
        self.body_name = body_name
        self.operation = operation

    def __str__(self) -> str:
        return self.message

    def __repr__(self) -> str:
        return (
            f"InvalidBodyError({self.message!r}, body_id={self.body_id}, "
            f"body_name={self.body_name!r}, operation={self.operation!r})"
        )


# =============================================================================
# CATEGORY: DATA NOT FOUND ERRORS
# =============================================================================


class DataNotFoundError(Error):
    """Base class for data lookup failures.

    This exception category covers all cases where required data could
    not be found, such as:
    - Unknown celestial body IDs
    - Fixed stars not in the catalog
    - SPK kernel files not found

    Catching this exception will catch all "not found" type errors.

    Example:
        >>> try:
        ...     ephem.calc_ut(jd, 99999, 0)  # Unknown body
        ... except ephem.DataNotFoundError as e:
        ...     print(f"Data not found: {e}")
    """

    pass


class UnknownBodyError(DataNotFoundError):
    """Error raised when an unknown or unsupported body ID is requested.

    This exception is raised when attempting to calculate the position of a
    celestial body using an unrecognized body ID. This helps users identify
    typos or unsupported body types early, rather than silently returning
    zero positions.

    Attributes:
        body_id: The unrecognized body ID that was requested
        message: Human-readable error message

    Example:
        >>> import libephemeris as ephem
        >>> try:
        ...     ephem.calc_ut(2451545.0, 99999, 0)  # Invalid body ID
        ... except ephem.UnknownBodyError as e:
        ...     print(f"Unknown body ID: {e.body_id}")
        Unknown body ID: 99999

    See Also:
        get_planet_name: Get name for a valid body ID
        SUN, MOON, etc.: Valid body ID constants
        InvalidBodyError: For valid bodies used in unsupported operations
    """

    def __init__(
        self,
        message: str,
        body_id: int | None = None,
    ):
        super().__init__(message)
        self.body_id = body_id
        self.message = message

    def __str__(self) -> str:
        return self.message

    def __repr__(self) -> str:
        return f"UnknownBodyError({self.message!r}, body_id={self.body_id})"


class StarNotFoundError(DataNotFoundError):
    """Error raised when a fixed star is not found in the catalog.

    This exception is raised when attempting to look up a fixed star by
    name, Hipparcos number, or other identifier that is not in the star
    catalog.

    Attributes:
        message: Human-readable error message
        star_id: The star identifier that was not found
        search_type: The type of search performed ("name", "hip", etc.)

    Example:
        >>> import libephemeris as ephem
        >>> try:
        ...     ephem.fixstar_ut("NonexistentStar", jd, 0)
        ... except ephem.StarNotFoundError as e:
        ...     print(f"Star not found: {e.star_id}")

    See Also:
        fixstar_ut: Fixed star position calculation
        fixstar2_ut: Fixed star position with HIP identifier
    """

    def __init__(
        self,
        message: str,
        star_id: str | None = None,
        search_type: str | None = None,
    ):
        super().__init__(message)
        self.message = message
        self.star_id = star_id
        self.search_type = search_type

    def __str__(self) -> str:
        return self.message

    def __repr__(self) -> str:
        return (
            f"StarNotFoundError({self.message!r}, "
            f"star_id={self.star_id!r}, search_type={self.search_type!r})"
        )


class SPKNotFoundError(DataNotFoundError):
    """Error raised when an SPK file is requested but not available.

    This exception is raised when attempting to register or load an SPK file
    that does not exist. It provides helpful instructions for obtaining the
    required SPK file.

    SPK (SPICE kernel) files contain high-precision ephemeris data for minor
    bodies (asteroids, TNOs, comets) and are downloaded from JPL Horizons.

    Attributes:
        filepath: Path to the SPK file that was not found
        body_name: Name of the celestial body (if known)
        body_id: JPL Horizons identifier (if known)
        message: Human-readable error message with instructions

    Example:
        >>> import libephemeris as ephem
        >>> try:
        ...     ephem.register_spk_body(ephem.CHIRON, "/path/to/chiron.bsp", 2002060)
        ... except ephem.SPKNotFoundError as e:
        ...     print(f"Missing SPK: {e.filepath}")
        ...     print("Instructions:", e.message)

    See Also:
        download_spk: Download SPK file from JPL Horizons
        download_and_register_spk: Download and register in one step
        enable_auto_spk: Enable automatic SPK downloading
    """

    def __init__(
        self,
        message: str,
        filepath: str | None = None,
        body_name: str | None = None,
        body_id: str | None = None,
    ):
        super().__init__(message)
        self.filepath = filepath
        self.body_name = body_name
        self.body_id = body_id
        self.message = message

    def __str__(self) -> str:
        return self.message

    def __repr__(self) -> str:
        return (
            f"SPKNotFoundError({self.message!r}, "
            f"filepath={self.filepath!r}, body_name={self.body_name!r}, "
            f"body_id={self.body_id!r})"
        )

    @classmethod
    def from_filepath(
        cls,
        filepath: str,
        body_name: str | None = None,
        body_id: str | None = None,
    ) -> "SPKNotFoundError":
        """Create an SPKNotFoundError with a helpful message.

        Args:
            filepath: Path to the missing SPK file
            body_name: Optional name of the celestial body
            body_id: Optional JPL Horizons identifier

        Returns:
            SPKNotFoundError with formatted instructions
        """
        lines = [f"SPK file not found: {filepath}"]
        lines.append("")
        lines.append("To obtain this SPK file, you can:")
        lines.append("")
        lines.append("1. Download manually using libephemeris.download_spk():")
        if body_id:
            lines.append("   >>> import libephemeris as eph")
            lines.append(f'   >>> eph.download_spk("{body_id}")')
        else:
            lines.append("   >>> import libephemeris as eph")
            lines.append('   >>> eph.download_spk("<body_id_or_name>")')
        lines.append("")
        lines.append("2. Download and register in one step:")
        if body_id and body_name:
            lines.append(
                f'   >>> eph.download_and_register_spk("{body_id}", eph.{body_name.upper()})'
            )
        else:
            lines.append(
                '   >>> eph.download_and_register_spk("<body_id>", <libephemeris_body_id>)'
            )
        lines.append("")
        lines.append("3. Enable automatic downloading (uses direct JPL Horizons HTTP):")
        lines.append("   >>> eph.set_auto_spk_download(True)")
        lines.append("")
        lines.append("4. Use the CLI download script:")
        lines.append("   $ python scripts/download_spk.py --bodies <body_name>")
        lines.append("")
        lines.append("For more information, see the libephemeris documentation.")

        message = "\n".join(lines)
        return cls(
            message=message,
            filepath=filepath,
            body_name=body_name,
            body_id=body_id,
        )


class SPKRequiredError(DataNotFoundError):
    """Error raised when SPK is required for precision but not available.

    This exception is raised in strict precision mode when calculating
    positions for major asteroids (Ceres, Pallas, Juno, Vesta, Chiron)
    without an SPK kernel registered. These bodies require SPK data for
    accurate calculations; the Keplerian fallback can have errors of
    several degrees.

    Attributes:
        body_id: The libephemeris body ID (CHIRON, CERES, etc.)
        body_name: Human-readable name of the body
        message: Human-readable error message with instructions

    Example:
        >>> import libephemeris as eph
        >>> eph.set_strict_precision(True)  # Enable strict mode
        >>> try:
        ...     eph.calc_ut(jd, eph.CHIRON, 0)  # No SPK registered
        ... except eph.SPKRequiredError as e:
        ...     print(f"SPK required for {e.body_name}")
        ...     # Either register SPK or disable strict mode

    See Also:
        set_strict_precision: Enable/disable strict precision mode
        download_and_register_spk: Download and register SPK for a body
        set_auto_spk_download: Enable automatic SPK downloading
    """

    def __init__(
        self,
        message: str,
        body_id: int | None = None,
        body_name: str | None = None,
    ):
        super().__init__(message)
        self.body_id = body_id
        self.body_name = body_name
        self.message = message

    def __str__(self) -> str:
        return self.message

    def __repr__(self) -> str:
        return (
            f"SPKRequiredError({self.message!r}, "
            f"body_id={self.body_id}, body_name={self.body_name!r})"
        )

    @classmethod
    def for_body(
        cls,
        body_id: int,
        body_name: str,
        horizons_id: str,
    ) -> "SPKRequiredError":
        """Create an SPKRequiredError with helpful instructions.

        Args:
            body_id: libephemeris body ID (CHIRON, etc.)
            body_name: Human-readable body name
            horizons_id: JPL Horizons identifier for downloading

        Returns:
            SPKRequiredError with formatted instructions
        """
        lines = [
            f"SPK kernel required for {body_name} (strict precision mode is enabled).",
            "",
            "The available Keplerian approximation is body- and date-dependent "
            "and is not classified as ephemeris-grade.",
            "For accurate calculations, you must either:",
            "",
            "1. Download and register SPK:",
            f'   >>> eph.download_and_register_spk("{horizons_id}", eph.{body_name.upper()}, '
            f'"1900-01-01", "2100-01-01")',
            "",
            "2. Outside sealed mode, enable automatic SPK download:",
            "   >>> eph.set_auto_spk_download(True)",
            "",
            "3. Disable strict precision mode (not recommended):",
            "   >>> eph.set_strict_precision(False)",
        ]
        message = "\n".join(lines)
        return cls(message=message, body_id=body_id, body_name=body_name)


# =============================================================================
# CATEGORY: CALCULATION ERRORS
# =============================================================================


class CalculationError(Error):
    """Base class for calculation and algorithm failures.

    This exception category covers all cases where a calculation or
    numerical algorithm failed, such as:
    - House calculations at polar latitudes
    - Dates outside ephemeris range
    - Iterative searches that didn't converge

    Catching this exception will catch all calculation-related errors.

    Example:
        >>> try:
        ...     ephem.houses(jd, 70.0, 0.0, ord('P'))  # Polar latitude
        ... except ephem.CalculationError as e:
        ...     print(f"Calculation failed: {e}")
    """

    pass


class PolarCircleError(CalculationError):
    """Error raised when house calculation fails at polar latitudes.

    Certain house systems (Placidus, Koch, Gauquelin) cannot be calculated
    at latitudes beyond the polar circle (~66.5°) because some ecliptic
    points never rise or set (circumpolar behavior).

    This exception provides detailed information about the polar circle
    condition to help users understand why the calculation failed and
    what alternatives are available.

    Attributes:
        latitude: The geographic latitude that triggered the error
        threshold: The polar circle threshold latitude for the given obliquity
        obliquity: The true obliquity of the ecliptic used in the calculation
        house_system: The requested house system that failed (e.g., 'P', 'K', 'G')
        message: Human-readable error message

    Example:
        >>> import libephemeris as ephem
        >>> try:
        ...     ephem.houses(jd, 70.0, 0.0, ord('P'))
        ... except ephem.PolarCircleError as e:
        ...     print(f"Polar error at {e.latitude}° (threshold: {e.threshold}°)")
        ...     # Use fallback house system
        ...     cusps, ascmc = ephem.houses(jd, 70.0, 0.0, ord('O'))

    See Also:
        houses_with_fallback: Automatically falls back to Porphyry
        get_polar_latitude_threshold: Returns the threshold for a given obliquity
    """

    def __init__(
        self,
        message: str,
        latitude: float | None = None,
        threshold: float | None = None,
        obliquity: float | None = None,
        house_system: str | None = None,
    ):
        super().__init__(message)
        self.latitude = latitude
        self.threshold = threshold
        self.obliquity = obliquity
        self.house_system = house_system
        self.message = message

    def __str__(self) -> str:
        return self.message

    def __repr__(self) -> str:
        return (
            f"PolarCircleError({self.message!r}, latitude={self.latitude}, "
            f"threshold={self.threshold}, obliquity={self.obliquity}, "
            f"house_system={self.house_system!r})"
        )


class EphemerisRangeError(CalculationError):
    """Error raised when a calculation date is outside the ephemeris range.

    This exception is raised when the requested Julian Day falls outside the
    date range covered by the loaded ephemeris file (e.g., DE440 covers
    1550-2650).

    This exception provides detailed information about the supported date range
    to help users understand why the calculation failed and what dates are valid.

    Attributes:
        requested_jd: The Julian Day that was requested
        start_jd: The start of the supported Julian Day range
        end_jd: The end of the supported Julian Day range
        start_date: The start date in calendar format (YYYY-MM-DD)
        end_date: The end date in calendar format (YYYY-MM-DD)
        body_id: The body ID that was being calculated (if available)
        body_name: Human-readable name of the body (if available)
        ephemeris_file: The ephemeris file in use (if available)
        message: Human-readable error message

    Example:
        >>> import libephemeris as ephem
        >>> try:
        ...     ephem.calc_ut(5000000.0, ephem.SUN, 0)  # Year ~8827 AD
        ... except ephem.EphemerisRangeError as e:
        ...     print(f"Date JD {e.requested_jd} is outside range")
        ...     print(f"Supported: JD {e.start_jd} to {e.end_jd}")
        ...     print(f"           ({e.start_date} to {e.end_date})")

    See Also:
        get_current_file_data: Returns information about the loaded ephemeris
        set_ephemeris_file: Load a different ephemeris with different date range
    """

    def __init__(
        self,
        message: str,
        requested_jd: float | None = None,
        start_jd: float | None = None,
        end_jd: float | None = None,
        start_date: str | None = None,
        end_date: str | None = None,
        body_id: int | None = None,
        body_name: str | None = None,
        ephemeris_file: str | None = None,
    ):
        super().__init__(message)
        self.requested_jd = requested_jd
        self.start_jd = start_jd
        self.end_jd = end_jd
        self.start_date = start_date
        self.end_date = end_date
        self.body_id = body_id
        self.body_name = body_name
        self.ephemeris_file = ephemeris_file
        self.message = message

    def __str__(self) -> str:
        return self.message

    def __repr__(self) -> str:
        return (
            f"EphemerisRangeError({self.message!r}, "
            f"requested_jd={self.requested_jd}, "
            f"start_jd={self.start_jd}, end_jd={self.end_jd}, "
            f"start_date={self.start_date!r}, end_date={self.end_date!r}, "
            f"body_id={self.body_id}, body_name={self.body_name!r}, "
            f"ephemeris_file={self.ephemeris_file!r})"
        )


class ConvergenceError(CalculationError):
    """Error raised when a numerical algorithm fails to converge.

    This exception is raised when an iterative search or calculation
    algorithm doesn't converge within the allowed number of iterations
    or tolerance limits. Common scenarios include:
    - Root-finding algorithms (e.g., Brent's method)
    - Eclipse/occultation searches
    - Crossing time calculations

    Attributes:
        message: Human-readable error message
        algorithm: Name of the algorithm that failed (e.g., "brent", "newton")
        iterations: Number of iterations attempted
        max_iterations: Maximum iterations allowed
        tolerance: Tolerance that was being sought
        last_value: Last computed value before failure (if available)

    Example:
        >>> import libephemeris as ephem
        >>> try:
        ...     ephem.solcross_ut(360.0, jd, 0)  # Search for crossing
        ... except ephem.ConvergenceError as e:
        ...     print(f"Algorithm '{e.algorithm}' failed after {e.iterations} iterations")

    See Also:
        solcross_ut: Sun crossing a longitude
        mooncross_ut: Moon crossing a longitude
    """

    def __init__(
        self,
        message: str,
        algorithm: str | None = None,
        iterations: int | None = None,
        max_iterations: int | None = None,
        tolerance: float | None = None,
        last_value: float | None = None,
    ):
        super().__init__(message)
        self.message = message
        self.algorithm = algorithm
        self.iterations = iterations
        self.max_iterations = max_iterations
        self.tolerance = tolerance
        self.last_value = last_value

    def __str__(self) -> str:
        return self.message

    def __repr__(self) -> str:
        return (
            f"ConvergenceError({self.message!r}, algorithm={self.algorithm!r}, "
            f"iterations={self.iterations}, max_iterations={self.max_iterations}, "
            f"tolerance={self.tolerance}, last_value={self.last_value})"
        )


# =============================================================================
# CATEGORY: CONFIGURATION ERRORS
# =============================================================================


class ConfigurationError(Error):
    """Error raised when required configuration or state is missing.

    This exception is raised when an operation requires specific
    configuration or state that hasn't been set up, such as:
    - Observer location not set for topocentric calculations
    - Ephemeris path not configured
    - Required dependencies not available

    Attributes:
        message: Human-readable error message
        missing_config: Name of the missing configuration item
        suggestion: Suggested action to resolve the error

    Example:
        >>> import libephemeris as ephem
        >>> try:
        ...     # Trying topocentric calc without setting location
        ...     ephem.calc_ut(jd, ephem.MOON, ephem.FLG_TOPOCTR)
        ... except ephem.ConfigurationError as e:
        ...     print(f"Missing: {e.missing_config}")
        ...     print(f"Suggestion: {e.suggestion}")
        Missing: observer_location
        Suggestion: Call set_topo(lon, lat, alt) first

    See Also:
        set_topo: Set observer location for topocentric calculations
        set_ephe_path: Set ephemeris file path
    """

    def __init__(
        self,
        message: str,
        missing_config: str | None = None,
        suggestion: str | None = None,
    ):
        super().__init__(message)
        self.message = message
        self.missing_config = missing_config
        self.suggestion = suggestion

    def __str__(self) -> str:
        return self.message

    def __repr__(self) -> str:
        return (
            f"ConfigurationError({self.message!r}, "
            f"missing_config={self.missing_config!r}, "
            f"suggestion={self.suggestion!r})"
        )


class NetworkSealedError(ConfigurationError, RuntimeError):
    """Raised before network I/O when the process policy is ``sealed``.

    ``RuntimeError`` is retained as a secondary base so existing optional
    downloader fallbacks that already catch runtime failures remain compatible.
    Callers that need to distinguish policy from transport failure should catch
    :class:`NetworkSealedError` or :class:`ConfigurationError`.
    """

    def __init__(self, message: str, purpose: str | None = None):
        super().__init__(
            message=message,
            missing_config="network_policy",
            suggestion="Provision data explicitly or select network_policy='allow'.",
        )
        self.purpose = purpose

    def __repr__(self) -> str:
        return f"NetworkSealedError({self.message!r}, purpose={self.purpose!r})"


class LEBCorruptionError(ValueError):
    """A LEB/LEB2 binary ephemeris file is corrupted or truncated.

    Deliberately subclasses ``ValueError`` for compatibility with existing
    reader callers. Source dispatchers nevertheless catch this distinct type
    first and treat it as a fatal provisioning error; corruption must never be
    confused with an ordinary range miss or trigger a source fallback.

    This is an internal robustness signal for the binary readers rather than
    a user-facing API exception; it is not re-exported from the package root.
    """


class FictitiousRuntimeDispatch(KeyError):
    """A fictitious body (40-58) was routed to its runtime analytical model.

    Deliberately subclasses ``KeyError`` so the LEB fast-path callers that
    catch ``(KeyError, ValueError)`` treat it as the ordinary miss signal and
    continue to the analytical dispatch. The distinct type lets the fallback
    logger recognize the miss as the by-design source selection for these
    bodies — their runtime model IS the declared local source — instead of a
    sealed-mode degradation worth a warning.

    This is an internal dispatch signal rather than a user-facing API
    exception; it is not re-exported from the package root.
    """


# =============================================================================
# VALIDATION HELPER FUNCTIONS
# =============================================================================


def validate_latitude(lat: float, func_name: str = "") -> None:
    """Validate that latitude is within valid range [-90, 90].

    Args:
        lat: Geographic latitude in degrees
        func_name: Optional function name for error message context

    Raises:
        CoordinateError: If latitude is outside [-90, 90] range

    Example:
        >>> validate_latitude(45.0)  # OK
        >>> validate_latitude(91.0)  # Raises CoordinateError
    """
    # NaN fails BOTH comparisons, so a bare range test lets it through and
    # the caller gets NaN-laced coordinates that look like a result. The
    # published contract for these surfaces promises a typed CoordinateError
    # for NaN and infinite input, so test finiteness explicitly.
    if not math.isfinite(lat) or lat < -90.0 or lat > 90.0:
        prefix = f"{func_name}: " if func_name else ""
        message = (
            f"{prefix}latitude {lat} is out of valid range. "
            f"Latitude must be between -90 and 90 degrees."
        )
        raise CoordinateError(
            message=message,
            coordinate_name="latitude",
            value=lat,
            min_value=-90.0,
            max_value=90.0,
        )


def validate_longitude(lon: float, func_name: str = "") -> None:
    """Validate that longitude is within valid range [-180, 360].

    Accepts both the signed convention ([-180, 180], west negative) and the
    east-positive convention ([0, 360)), since the reference API
    interprets geographic longitude modulo 360 degrees. The value is used
    as-is: the underlying sidereal-time and observer geometry are periodic in
    360 degrees, so 200 and -160 yield identical results.

    Args:
        lon: Geographic longitude in degrees
        func_name: Optional function name for error message context

    Raises:
        CoordinateError: If longitude is outside [-180, 360] range

    Example:
        >>> validate_longitude(12.5)   # OK (signed)
        >>> validate_longitude(200.0)  # OK (east-positive, = -160 deg)
        >>> validate_longitude(400.0)  # Raises CoordinateError
    """
    # Finiteness first, for the same reason as validate_latitude.
    if not math.isfinite(lon) or lon < -180.0 or lon > 360.0:
        prefix = f"{func_name}: " if func_name else ""
        message = (
            f"{prefix}longitude {lon} is out of valid range. "
            f"Longitude must be between -180 and 360 degrees."
        )
        raise CoordinateError(
            message=message,
            coordinate_name="longitude",
            value=lon,
            min_value=-180.0,
            max_value=360.0,
        )


def validate_coordinates(lat: float, lon: float, func_name: str = "") -> None:
    """Validate both latitude and longitude coordinates.

    Args:
        lat: Geographic latitude in degrees
        lon: Geographic longitude in degrees
        func_name: Optional function name for error message context

    Raises:
        CoordinateError: If latitude is outside [-90, 90] or
                        longitude is outside [-180, 360]

    Example:
        >>> validate_coordinates(41.9, 12.5)  # OK (Rome)
        >>> validate_coordinates(91.0, 0.0)  # Raises CoordinateError
    """
    validate_latitude(lat, func_name)
    validate_longitude(lon, func_name)


def validate_jd_range(
    jd: float, body_id: int | None = None, func_name: str = ""
) -> None:
    """Validate that Julian Day is within the supported ephemeris range.

    This function checks if the given Julian Day falls within the date range
    covered by the currently loaded ephemeris file. If the date is outside
    the range, an EphemerisRangeError is raised with detailed information.

    Args:
        jd: Julian Day number (UT1 or TT) to validate
        body_id: Optional body ID being calculated (for error message)
        func_name: Optional function name for error message context

    Raises:
        EphemerisRangeError: If the Julian Day is outside the ephemeris range

    Note:
        In LEB mode the stored per-body header is authoritative and no JPL
        kernel is opened. Other modes validate against the active DE kernel.

    Example:
        >>> from libephemeris.exceptions import validate_jd_range
        >>> validate_jd_range(2451545.0)  # J2000.0 - OK for DE440
        >>> validate_jd_range(3000000.0)  # Year ~3501 - raises EphemerisRangeError
    """
    import math

    from . import state
    from .planets import get_planet_name

    # Reject non-finite Julian Days up front. A NaN slips through every
    # ``jd < start`` / ``jd > end`` comparison below (all comparisons with NaN
    # are False), so without this guard it would sail past range validation and
    # surface as garbage output (e.g. calc_pctr) or an opaque "cannot convert
    # float NaN/infinity to integer" deep in the pipeline. Raise the same
    # EphemerisRangeError the out-of-range branch uses so callers see one
    # consistent, informative failure for every non-finite input.
    if not math.isfinite(jd):
        body_name = get_planet_name(body_id) if body_id is not None else None
        prefix = f"{func_name}: " if func_name else ""
        who = (
            f"{body_name} (ID {body_id})"
            if body_name and body_id is not None
            else "the requested body"
        )
        raise EphemerisRangeError(
            message=(
                f"{prefix}Cannot calculate {who} for JD {jd}: "
                f"Julian Day must be a finite number."
            ),
            requested_jd=jd,
            body_id=body_id,
            body_name=body_name,
        )

    if state.get_calc_mode() == "leb":
        if body_id is None:
            reader = state.get_leb_reader()
            if reader is None:
                return
            path = str(getattr(reader, "path", ""))
            start_jd, end_jd = reader.jd_range
            denum = 0
        else:
            from .inventory import get_body_coverage

            coverage = get_body_coverage(body_id, jd)
            if coverage is None:
                # The downstream LEB dispatcher distinguishes a legitimate
                # local model from a missing required state. Do not inspect a
                # DE kernel merely to validate a body it will not serve.
                return
            path = coverage.data_file or ""
            start_jd = coverage.jd_start
            end_jd = coverage.jd_end
            denum = 0
    else:
        # Ensure the selected JPL ephemeris is loaded before inspecting range.
        state.get_planets()
        path, start_jd, end_jd, denum = state.get_current_file_data(0)

    # If we couldn't get range information, skip validation
    # (this allows fallback behavior for special cases)
    if start_jd == 0.0 and end_jd == 0.0:
        return

    # Check if JD is within range
    if jd < start_jd or jd > end_jd:
        # Convert JD to calendar date for message
        from .time_utils import revjul

        req_year, req_month, req_day, req_hour = revjul(jd, 1)  # Gregorian
        req_date_str = f"{req_year}-{req_month:02d}-{req_day:02d}"

        start_year, start_month, start_day, _ = revjul(start_jd, 1)
        start_date = f"{start_year}-{start_month:02d}-{start_day:02d}"

        end_year, end_month, end_day, _ = revjul(end_jd, 1)
        end_date = f"{end_year}-{end_month:02d}-{end_day:02d}"

        # Get body name if available
        body_name = None
        if body_id is not None:
            body_name = get_planet_name(body_id)

        # Get ephemeris filename
        import os

        ephemeris_file = None
        if path:
            ephemeris_file = os.path.basename(path)
        elif denum:
            ephemeris_file = f"de{denum}.bsp"

        # Build error message
        parts = []
        prefix = f"{func_name}: " if func_name else ""

        if body_name and body_id is not None:
            parts.append(f"{prefix}Cannot calculate {body_name} (ID {body_id})")
        else:
            parts.append(f"{prefix}Calculation failed")

        parts.append(f"for JD {jd:.6f} ({req_date_str}):")
        parts.append("date is outside ephemeris range.")
        parts.append(f"\n  Supported range: JD {start_jd:.1f} to {end_jd:.1f}")
        parts.append(f" ({start_date} to {end_date})")

        if ephemeris_file:
            parts.append(f"\n  Ephemeris file: {ephemeris_file}")

        message = " ".join(parts[:3]) + "".join(parts[3:])

        raise EphemerisRangeError(
            message=message,
            requested_jd=jd,
            start_jd=start_jd,
            end_jd=end_jd,
            start_date=start_date,
            end_date=end_date,
            body_id=body_id,
            body_name=body_name,
            ephemeris_file=ephemeris_file,
        )
