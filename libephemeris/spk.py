# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
SPK kernel support for high-precision minor body calculations.

This module provides functionality to:
- Download SPK (SPICE kernel) files from JPL Horizons API
- Register mappings between libephemeris body IDs and SPK targets
- Calculate positions using SPK data instead of Keplerian approximations
- Support for both SPK type 2/3 (Skyfield) and type 21 (spktype21)

Using SPK kernels provides significantly higher precision for asteroids and TNOs:
- Keplerian model: ~10-30 arcseconds (asteroids), ~1-3 arcminutes (TNOs)
- SPK kernel: ~arcseconds to sub-arcsecond (within kernel coverage)

SPK Type Support:
- Type 2, 3: Chebyshev polynomials (supported by Skyfield)
- Type 21: Extended Modified Difference Arrays (supported by spktype21)

JPL Horizons returns type 21 for most asteroids and comets. This module
automatically detects the SPK type and uses the appropriate library.

Usage:
    >>> import libephemeris as eph
    >>> # Download and register in one step
    >>> eph.download_and_register_spk(
    ...     body="Chiron",
    ...     ipl=eph.CHIRON,
    ...     start="2000-01-01",
    ...     end="2100-01-01",
    ... )
    >>> # Now calc_ut uses SPK automatically
    >>> pos, _ = eph.calc_ut(2451545.0, eph.CHIRON, eph.FLG_SPEED)

References:
    - JPL Horizons API: https://ssd-api.jpl.nasa.gov/doc/horizons.html
    - NAIF SPICE: https://naif.jpl.nasa.gov/naif/
    - spktype21: https://pypi.org/project/spktype21/
"""

from __future__ import annotations

import json
import math
import os
import re
import ssl
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Optional, Union

import certifi
import numpy as np

from skyfield.framelib import ecliptic_frame

from .download import _is_valid_bsp
from .exceptions import SPKNotFoundError
from .logging_config import format_file_size, get_logger
from .state import get_timescale

# Vendored spktype21 for SPK type 21 support (upstream unmaintained since 2018).
# Includes numpy 2.x compatibility fix (.item() on map_array results).
_spktype21_module = None


def _get_spktype21():
    """Lazy load vendored spktype21 module."""
    global _spktype21_module
    if _spktype21_module is None:
        try:
            from .vendor.spktype21 import SPKType21

            _spktype21_module = SPKType21
        except ImportError:
            _spktype21_module = False  # Mark as unavailable
    return _spktype21_module if _spktype21_module is not False else None


# Threshold for showing progress bar (1 MB)
_PROGRESS_THRESHOLD_BYTES = 1 * 1024 * 1024

# Obliquity of J2000 for ICRS to ecliptic rotation
_OBLIQUITY_J2000_RAD = math.radians(23.4392911)


# =============================================================================
# MODULE-LEVEL STATE (managed by state.py, but defined here for type hints)
# =============================================================================
# Actual state is in state.py to maintain single source of truth

# NAIF ID convention for numbered asteroids:
# JPL Horizons uses: naif_id = asteroid_number + 20000000
# Some older references use: naif_id = asteroid_number + 2000000
# We support both by detecting the actual NAIF ID from SPK files
NAIF_ASTEROID_OFFSET = 2000000  # Legacy/default offset
NAIF_ASTEROID_OFFSET_HORIZONS = 20000000  # JPL Horizons SPK offset


# =============================================================================
# NAIF ID UTILITIES
# =============================================================================


def _extract_asteroid_number(body: str) -> Optional[int]:
    """
    Extract asteroid catalog number from body string.

    Args:
        body: Body identifier (e.g., "Chiron", "2060", "2060 Chiron", "(2060)")

    Returns:
        Asteroid number if found, None otherwise.

    Examples:
        >>> _extract_asteroid_number("2060")
        2060
        >>> _extract_asteroid_number("2060 Chiron")
        2060
        >>> _extract_asteroid_number("(136199) Eris")
        136199
        >>> _extract_asteroid_number("Chiron")
        None
    """
    # Try to find a number at the start or in parentheses
    # Pattern: optional parentheses, digits, optional name
    match = re.match(r"^\(?(\d+)\)?", body.strip())
    if match:
        return int(match.group(1))
    return None


def _deduce_naif_id(body: str, asteroid_number: Optional[int] = None) -> Optional[int]:
    """
    Deduce NAIF ID from body identifier.

    For numbered asteroids, NAIF ID = asteroid_number + 2000000.
    For unnumbered objects or named objects without numbers, returns None.

    Args:
        body: Body identifier string
        asteroid_number: Optional pre-extracted asteroid number

    Returns:
        NAIF ID if deducible, None otherwise.
    """
    if asteroid_number is None:
        asteroid_number = _extract_asteroid_number(body)

    if asteroid_number is not None:
        return asteroid_number + NAIF_ASTEROID_OFFSET

    return None


def _get_spk_targets(filepath: str) -> list[int]:
    """
    Get the list of target NAIF IDs in an SPK file.

    Args:
        filepath: Path to the SPK file

    Returns:
        List of target NAIF IDs in the file.
    """
    try:
        from jplephem.spk import SPK

        spk = SPK.open(filepath)
        try:
            targets = [int(seg.target) for seg in spk.segments]
        finally:
            spk.close()
        return targets
    except (OSError, ValueError, KeyError, IndexError):
        return []


def _find_naif_id_for_asteroid(filepath: str, asteroid_number: int) -> Optional[int]:
    """
    Find the NAIF ID for an asteroid number in an SPK file.

    JPL Horizons uses 20000000 + asteroid_number, while some older conventions
    use 2000000 + asteroid_number. This function checks the SPK file to find
    which convention is used.

    Args:
        filepath: Path to the SPK file
        asteroid_number: Asteroid catalog number (e.g., 2060 for Chiron)

    Returns:
        The actual NAIF ID found in the file, or None if not found.
    """
    targets = _get_spk_targets(filepath)

    # Try both conventions
    horizons_naif = asteroid_number + NAIF_ASTEROID_OFFSET_HORIZONS  # 20000000 + num
    legacy_naif = asteroid_number + NAIF_ASTEROID_OFFSET  # 2000000 + num

    if horizons_naif in targets:
        return horizons_naif
    if legacy_naif in targets:
        return legacy_naif

    # Return None if neither found
    return None


# =============================================================================
# BODY NAME UTILITIES
# =============================================================================

# Mapping from body IDs to human-readable names for helpful error messages
# This complements the SPK_BODY_NAME_MAP in constants.py which maps to Horizons IDs
_BODY_NAMES: dict[int, str] = {
    15: "Chiron",  # CHIRON
    16: "Pholus",  # PHOLUS
    17: "Ceres",  # CERES
    18: "Pallas",  # PALLAS
    19: "Juno",  # JUNO
    20: "Vesta",  # VESTA
    10000 + 136199: "Eris",  # ERIS
    10000 + 90377: "Sedna",  # SEDNA
    10000 + 136108: "Haumea",  # HAUMEA
    10000 + 136472: "Makemake",  # MAKEMAKE
    10000 + 28978: "Ixion",  # IXION
    10000 + 90482: "Orcus",  # ORCUS
    10000 + 50000: "Quaoar",  # QUAOAR
    10000 + 20000: "Varuna",  # VARUNA
    10000 + 7066: "Nessus",  # NESSUS
    10000 + 8405: "Asbolus",  # ASBOLUS
    10000 + 10199: "Chariklo",  # CHARIKLO
    10000 + 225088: "Gonggong",  # GONGGONG
    10000 + 99942: "Apophis",  # APOPHIS
    10000 + 10: "Hygiea",  # HYGIEA
    10000 + 433: "Eros",  # EROS
}


def _get_body_name(ipl: int) -> Optional[str]:
    """
    Get the human-readable name for a body ID.

    Args:
        ipl: libephemeris body ID

    Returns:
        Human-readable name if known, None otherwise.
    """
    return _BODY_NAMES.get(ipl)


def _get_horizons_id_for_body(ipl: int) -> Optional[str]:
    """
    Get the JPL Horizons target identifier for a body ID.

    Args:
        ipl: libephemeris body ID

    Returns:
        Horizons ID if available, None otherwise.
    """
    from .constants import SPK_BODY_NAME_MAP

    if ipl in SPK_BODY_NAME_MAP:
        return SPK_BODY_NAME_MAP[ipl][0]
    return None


# =============================================================================
# HORIZONS API DOWNLOAD
# =============================================================================

# JPL Horizons API endpoint
HORIZONS_API_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"


def _build_horizons_url(
    body: str,
    start: str,
    end: str,
    center: str = "500@0",
) -> str:
    """
    Build URL for JPL Horizons SPK file request.

    Args:
        body: Target body (name or number)
        start: Start date (YYYY-MM-DD)
        end: End date (YYYY-MM-DD)
        center: Reference center (default: 500@0 = Solar System Barycenter)

    Returns:
        URL string for Horizons API request.
    """
    # URL-encode parameters properly using urllib.parse.urlencode.
    # This is critical for body names containing special characters like
    # "Ceres;" where the semicolon must be encoded as %3B.
    params = {
        "format": "json",
        "COMMAND": f"'{body}'",
        "OBJ_DATA": "NO",
        "MAKE_EPHEM": "YES",
        "EPHEM_TYPE": "SPK",
        "CENTER": f"'{center}'",
        "START_TIME": f"'{start}'",
        "STOP_TIME": f"'{end}'",
    }

    query = urllib.parse.urlencode(params, safe="'@")
    return f"{HORIZONS_API_URL}?{query}"


def _sanitize_filename(body: str) -> str:
    """
    Create a safe filename from body name.

    Args:
        body: Body name/identifier

    Returns:
        Sanitized string safe for use in filenames.
    """
    # Remove/replace unsafe characters
    safe = re.sub(r"[^\w\-]", "_", body.lower())
    # Collapse multiple underscores
    safe = re.sub(r"_+", "_", safe)
    # Remove leading/trailing underscores
    safe = safe.strip("_")
    return safe or "body"


def download_spk(
    body: str,
    start: str,
    end: str,
    path: Optional[str] = None,
    center: str = "500@0",
    overwrite: bool = False,
    timeout: int = 120,
) -> str:
    """
    Download SPK ephemeris file from JPL Horizons.

    Requests an SPK (SPICE kernel) file for the specified body from JPL Horizons
    API. The file contains high-precision ephemeris data that can be used for
    accurate position calculations.

    Args:
        body: Target body identifier. Can be:
            - Asteroid number: "2060", "136199"
            - Name: "Chiron", "Eris"
            - Combined: "2060 Chiron", "(136199) Eris"
        start: Start date in YYYY-MM-DD format
        end: End date in YYYY-MM-DD format
        path: Directory to save the file. If None, uses ~/.libephemeris/spk/
        center: Reference center for ephemeris (default: "500@0" = SSB).
            Use "500@0" for compatibility with Skyfield/DE kernels.
        overwrite: If True, overwrite existing file. If False, skip if exists.
        timeout: Request timeout in seconds (default: 120)

    Returns:
        str: Full path to the downloaded SPK file

    Raises:
        ValueError: If body not found on Horizons or invalid parameters
        ConnectionError: If network request fails after retries
        IOError: If file cannot be written

    Example:
        >>> path = download_spk("Chiron", "2000-01-01", "2100-01-01")
        >>> print(path)
        /path/to/kernels/chiron_2000_2100.bsp
    """
    # Determine output directory
    if path is None:
        from .state import _get_data_dir

        path = os.path.join(_get_data_dir(), "spk")

    # Ensure directory exists
    os.makedirs(path, exist_ok=True)

    # Generate filename
    body_safe = _sanitize_filename(body)
    start_short = start.replace("-", "")[:6]  # YYYYMM
    end_short = end.replace("-", "")[:6]
    filename = f"{body_safe}_{start_short}_{end_short}.bsp"
    filepath = os.path.join(path, filename)

    logger = get_logger()

    # Check if already exists and valid
    if os.path.exists(filepath) and not overwrite:
        if _is_valid_bsp(filepath):
            return filepath
        logger.warning("Cached SPK file %s is corrupted, re-downloading", filepath)
        try:
            os.remove(filepath)
        except OSError:
            pass

    # Deduce NAIF ID for logging (if possible)
    naif_id = _deduce_naif_id(body)
    if naif_id is not None:
        logger.info(
            "Downloading SPK for %s (NAIF %d) from JPL Horizons...", body, naif_id
        )
    else:
        logger.info("Downloading SPK for %s from JPL Horizons...", body)

    # Build request URL
    url = _build_horizons_url(body, start, end, center)

    # Make request with retry
    max_retries = 2
    last_error: Optional[Exception] = None

    for attempt in range(max_retries + 1):
        try:
            # Create SSL context with certifi certificates
            ssl_context = ssl.create_default_context(cafile=certifi.where())

            # First request: get SPK file URL from Horizons
            req = urllib.request.Request(url)
            req.add_header("User-Agent", "libephemeris/0.1")

            with urllib.request.urlopen(
                req, timeout=timeout, context=ssl_context
            ) as response:
                data = json.loads(response.read().decode("utf-8"))

            # Check for errors in response
            if "error" in data:
                raise ValueError(f"Horizons API error: {data['error']}")

            # Get SPK file URL from response
            if "spk_file_id" not in data or "spk" not in data:
                # Check if it's a "result" response (body not found, etc.)
                if "result" in data:
                    raise ValueError(f"Horizons error: {data['result'][:500]}")
                raise ValueError(
                    f"Unexpected Horizons response format. Keys: {list(data.keys())}"
                )

            spk_data_b64 = data["spk"]

            import base64

            try:
                spk_data = base64.b64decode(spk_data_b64)
            except (OSError, ValueError, KeyError, IndexError) as e:
                raise ValueError(
                    f"Failed to decode SPK data from Horizons response: {e}"
                ) from e

            temp_fd, temp_path = tempfile.mkstemp(dir=path, suffix=".download")
            try:
                with os.fdopen(temp_fd, "wb") as f:
                    f.write(spk_data)

                if not _is_valid_bsp(temp_path):
                    os.unlink(temp_path)
                    raise ValueError(
                        "Downloaded SPK data failed validation - corrupt or incomplete"
                    )

                from .download import publish_temp_file

                publish_temp_file(temp_path, filepath)
            except (OSError, ValueError, KeyError, IndexError):
                if os.path.exists(temp_path):
                    os.unlink(temp_path)
                raise

            size_str = format_file_size(len(spk_data))
            logger.info("SPK saved: %s (%s)", filepath, size_str)
            return filepath

        except urllib.error.HTTPError as e:
            last_error = e
            if e.code == 400:
                # Bad request - likely invalid body name
                try:
                    error_data = json.loads(e.read().decode("utf-8"))
                    msg = error_data.get("message", str(e))
                except (OSError, ValueError, KeyError, IndexError):
                    msg = str(e)
                logger.warning("SPK download failed for %s: %s", body, msg)
                raise ValueError(f"Invalid body '{body}': {msg}") from e
            elif attempt < max_retries:
                logger.warning(
                    "SPK download attempt %d failed for %s: HTTP %d, retrying...",
                    attempt + 1,
                    body,
                    e.code,
                )
                time.sleep(1)  # Brief pause before retry
                continue
            else:
                logger.warning(
                    "SPK download failed for %s: HTTP %d after %d attempts",
                    body,
                    e.code,
                    max_retries + 1,
                )
                raise ConnectionError(
                    f"Failed to download SPK after {max_retries + 1} attempts: {e}"
                ) from e

        except urllib.error.URLError as e:
            last_error = e
            if attempt < max_retries:
                logger.warning(
                    "SPK download attempt %d failed for %s: %s, retrying...",
                    attempt + 1,
                    body,
                    e.reason,
                )
                time.sleep(1)
                continue
            else:
                logger.warning("SPK download failed for %s: %s", body, e.reason)
                raise ConnectionError(
                    f"Network error downloading SPK: {e.reason}"
                ) from e

    # Should not reach here, but just in case
    raise ConnectionError(
        f"Download failed: {last_error}"
    )  # pragma: no cover - retry loop always returns or raises


# =============================================================================
# SPK REGISTRATION
# =============================================================================


def register_spk_body(
    ipl: int,
    spk_file: str,
    naif_id: Union[int, str],
) -> None:
    """
    Register a mapping between a libephemeris body ID and an SPK kernel target.

    After registration, calc_ut() and calc() will automatically use the
    SPK kernel for this body instead of the Keplerian approximation.

    Args:
        ipl: libephemeris body ID (e.g., CHIRON, ERIS)
        spk_file: Path to the SPK file, or filename if in library path
        naif_id: NAIF ID of the body in the SPK kernel.
            For numbered asteroids: asteroid_number + 2000000
            (e.g., Chiron = 2060 -> naif_id = 2002060)

    Raises:
        SPKNotFoundError: If SPK file not found (with helpful instructions)
        ValueError: If naif_id not found in SPK kernel

    Example:
        >>> register_spk_body(CHIRON, "/path/to/chiron.bsp", 2002060)
        >>> # Now calc_ut(jd, CHIRON, ...) uses SPK data
    """
    from . import state

    # Convert naif_id to int if string
    if isinstance(naif_id, str):
        naif_id = int(naif_id)

    # Resolve file path
    if not os.path.isabs(spk_file):
        from .state import _get_data_dir

        data_dir = _get_data_dir()
        # Try data dir and spk subdirectory
        for search_dir in [os.path.join(data_dir, "spk"), data_dir]:
            full_path = os.path.join(search_dir, spk_file)
            if os.path.exists(full_path):
                spk_file = full_path
                break
        else:
            if not os.path.exists(spk_file):
                # Get helpful info for error message
                body_name = _get_body_name(ipl)
                body_id = _get_horizons_id_for_body(ipl)
                raise SPKNotFoundError.from_filepath(
                    filepath=spk_file,
                    body_name=body_name,
                    body_id=body_id,
                )

    if not os.path.exists(spk_file):
        # Get helpful info for error message
        body_name = _get_body_name(ipl)
        body_id = _get_horizons_id_for_body(ipl)
        raise SPKNotFoundError.from_filepath(
            filepath=spk_file,
            body_name=body_name,
            body_id=body_id,
        )

    # Detect SPK type
    spk_type = _detect_spk_type(spk_file)

    if spk_type == 21:
        # Type 21: Use spktype21 library
        kernel = _load_type21_kernel(spk_file)
        if kernel is None:
            raise ValueError(
                f"Failed to load type 21 SPK file {spk_file}. "
                "Install spktype21: pip install spktype21"
            )
        # Note: spktype21 doesn't provide a way to list segments easily,
        # so we skip validation for type 21 kernels. The computation will
        # fail at runtime if the NAIF ID is not in the kernel.
    else:
        # Type 2/3: Use Skyfield
        state._load_spk_kernel(spk_file)

        # Validate that naif_id exists in kernel
        kernel = state._SPK_KERNELS.get(spk_file)
        if kernel is not None:
            # Check if target exists in kernel
            try:
                _ = kernel[naif_id]
            except KeyError:
                # Try with string name
                try:
                    _ = kernel[str(naif_id)]
                except KeyError:
                    available = []
                    if hasattr(kernel, "names"):
                        available = list(kernel.names())[:10]
                    raise ValueError(
                        f"NAIF ID {naif_id} not found in SPK kernel {spk_file}. "
                        f"Available targets (first 10): {available}"
                    )

    # Register mapping (include SPK type for later use)
    state._SPK_BODY_MAP[ipl] = (spk_file, naif_id)


def unregister_spk_body(ipl: int) -> None:
    """
    Remove SPK registration for a body.

    After unregistration, the body will use Keplerian approximation again.

    Args:
        ipl: libephemeris body ID (e.g., CHIRON)
    """
    from . import state

    if ipl in state._SPK_BODY_MAP:
        del state._SPK_BODY_MAP[ipl]


def get_spk_body_info(ipl: int) -> Optional[tuple[str, int]]:
    """
    Get SPK registration info for a body.

    Args:
        ipl: libephemeris body ID

    Returns:
        Tuple of (spk_file, naif_id) if registered, None otherwise.
    """
    from . import state

    return state._SPK_BODY_MAP.get(ipl)


def list_spk_bodies() -> dict[int, tuple[str, int]]:
    """
    List all registered SPK body mappings.

    Returns:
        Dict mapping ipl -> (spk_file, naif_id) for all registered bodies.
    """
    from . import state

    return dict(state._SPK_BODY_MAP)


def get_spk_coverage(spk_file: str) -> Optional[tuple[float, float]]:
    """
    Get the time coverage of an SPK kernel.

    Supports both type 2/3 (Skyfield) and type 21 (spktype21) kernels.

    Args:
        spk_file: Path to SPK file

    Returns:
        Tuple of (start_jd, end_jd) Julian Day coverage, or None if error.
    """
    from . import state

    # Resolve path
    if not os.path.isabs(spk_file):
        from .state import _get_data_dir

        data_dir = _get_data_dir()
        for search_dir in [os.path.join(data_dir, "spk"), data_dir]:
            full_path = os.path.join(search_dir, spk_file)
            if os.path.exists(full_path):
                spk_file = full_path
                break

    if not os.path.exists(spk_file):
        return None

    # Detect SPK type
    spk_type = _detect_spk_type(spk_file)

    if spk_type == 21:
        # Use jplephem to get coverage for type 21
        try:
            from jplephem.spk import SPK

            spk = SPK.open(spk_file)
            try:
                start_jd = None
                end_jd = None
                for segment in spk.segments:
                    if hasattr(segment, "start_jd") and hasattr(segment, "end_jd"):
                        seg_start = float(segment.start_jd)
                        seg_end = float(segment.end_jd)
                        if start_jd is None or seg_start < start_jd:
                            start_jd = seg_start
                        if end_jd is None or seg_end > end_jd:
                            end_jd = seg_end
            finally:
                spk.close()
            if start_jd is not None and end_jd is not None:
                return (start_jd, end_jd)
        except (OSError, ValueError, KeyError, IndexError):
            pass
        return None

    # Type 2/3: Use Skyfield
    if spk_file not in state._SPK_KERNELS:
        try:
            state._load_spk_kernel(spk_file)
        except ValueError:
            return None

    kernel = state._SPK_KERNELS.get(spk_file)
    if kernel is None:
        return None

    try:
        if hasattr(kernel, "spk") and hasattr(kernel.spk, "segments"):
            segments = list(kernel.spk.segments)
            if segments:
                start_jd = min(
                    float(s.start_jd) for s in segments if hasattr(s, "start_jd")
                )
                end_jd = max(float(s.end_jd) for s in segments if hasattr(s, "end_jd"))
                return (start_jd, end_jd)
    except (OSError, ValueError, KeyError, IndexError):
        pass

    return None


# =============================================================================
# SPK TYPE DETECTION AND TYPE 21 SUPPORT
# =============================================================================


def _detect_spk_type(filepath: str) -> Optional[int]:
    """
    Detect the data type of an SPK file.

    SPK files can contain multiple segment types. This function checks
    for the presence of type 21 segments (Extended Modified Difference Arrays),
    which are commonly returned by JPL Horizons for asteroids and comets.

    Args:
        filepath: Path to the SPK file

    Returns:
        21 if file contains type 21 segments, 2 for type 2/3, None if error.
    """
    try:
        from jplephem.spk import SPK

        spk = SPK.open(filepath)
        try:
            for segment in spk.segments:
                if hasattr(segment, "data_type") and segment.data_type == 21:
                    return 21
            return 2  # Assume type 2/3 for Skyfield compatibility
        finally:
            spk.close()
    except (OSError, ValueError, KeyError, IndexError):
        return None


def _load_type21_kernel(filepath: str):
    """
    Load an SPK type 21 kernel using spktype21 library.

    Args:
        filepath: Path to the SPK file

    Returns:
        SPKType21 object or None if loading fails
    """
    from . import state

    # Check cache first
    if filepath in state._SPK_TYPE21_KERNELS:
        return state._SPK_TYPE21_KERNELS[filepath]

    SPKType21 = _get_spktype21()
    if SPKType21 is None:
        get_logger().warning(
            "spktype21 module not available. Install with: pip install spktype21"
        )
        return None

    try:
        kernel = SPKType21.open(filepath)
        state._SPK_TYPE21_KERNELS[filepath] = kernel
        return kernel
    except (OSError, ValueError, KeyError, IndexError) as e:
        get_logger().warning("Failed to load type 21 SPK %s: %s", filepath, e)
        return None


def _icrs_to_ecliptic_j2000(x: float, y: float, z: float) -> tuple:
    """
    Rotate coordinates from ICRS (equatorial J2000) to ecliptic J2000.

    SPK files from JPL Horizons contain coordinates in ICRS (equatorial).
    This function converts to ecliptic J2000 for consistency with
    libephemeris calculations.

    Args:
        x, y, z: Cartesian coordinates in ICRS

    Returns:
        Tuple of (x_ecl, y_ecl, z_ecl) in ecliptic J2000
    """
    ce = math.cos(_OBLIQUITY_J2000_RAD)
    se = math.sin(_OBLIQUITY_J2000_RAD)

    x_ecl = x
    y_ecl = y * ce + z * se
    z_ecl = -y * se + z * ce

    return (x_ecl, y_ecl, z_ecl)


# =============================================================================
# SKYFIELD VECTORFUNCTION WRAPPER FOR SPK TYPE 21
# =============================================================================


class _SpkType21Target:
    """Skyfield VectorFunction-compatible wrapper for SPK type 21 asteroids.

    Combines spktype21 heliocentric positions with the Sun's SSB position
    to produce SSB-centered ICRS positions, compatible with Skyfield's
    observe()/apparent() pipeline.

    This ensures asteroids use the same coordinate transformation path as
    planets (ICRS → Skyfield frame_latlon), avoiding the ~0.3" systematic
    error from the manual ecliptic J2000 + precession/nutation approach.
    """

    def __init__(self, kernel, naif_id: int, sun_target):
        """Initialize with spktype21 kernel and Sun target.

        Args:
            kernel: SPKType21 kernel object
            naif_id: NAIF ID of the asteroid in the kernel
            sun_target: Skyfield VectorFunction for the Sun (from SSB)
        """
        self._kernel = kernel
        self._naif_id = naif_id
        self._sun = sun_target
        self.center = 0  # SSB-centered (required for observe())
        self.target = naif_id

    def _at(self, t):
        """Compute SSB-centered ICRS position at time t.

        Matches the Skyfield VectorFunction._at() protocol:
        Returns (position_au, velocity_au_per_d, gcrs_position, message).

        Args:
            t: Skyfield Time object

        Returns:
            Tuple of (position, velocity, None, None)
        """
        AU_KM = 149597870.7

        # Get Sun SSB position (ICRS, AU)
        sun_pos = self._sun.at(t)
        sun_au = sun_pos.position.au
        sun_vel = sun_pos.velocity.au_per_d

        # Get heliocentric position from spktype21
        # Handle both scalar and array times
        jd_tdb = t.tdb
        if np.ndim(jd_tdb) == 0:
            # Scalar time
            pos_km, vel_km_s = self._kernel.compute_type21(
                10, self._naif_id, float(jd_tdb)
            )
            pos_au = np.array(pos_km) / AU_KM
            vel_au_d = np.array(vel_km_s) * 86400.0 / AU_KM
            # SSB = Sun + heliocentric
            position = (
                sun_au + pos_au.reshape(3, 1) if sun_au.ndim > 1 else sun_au + pos_au
            )
            velocity = (
                sun_vel + vel_au_d.reshape(3, 1)
                if sun_vel.ndim > 1
                else sun_vel + vel_au_d
            )
        else:
            # Array of times
            n = len(jd_tdb)
            position = np.empty((3, n))
            velocity = np.empty((3, n))
            for i in range(n):
                pos_km, vel_km_s = self._kernel.compute_type21(
                    10, self._naif_id, float(jd_tdb[i])
                )
                pos_au = np.array(pos_km) / AU_KM
                vel_au_d = np.array(vel_km_s) * 86400.0 / AU_KM
                position[:, i] = sun_au[:, i] + pos_au
                velocity[:, i] = sun_vel[:, i] + vel_au_d

        return position, velocity, None, None

    def at(self, t):
        """Return ICRF position at time t."""
        from skyfield.positionlib import ICRF

        position, velocity, _, _ = self._at(t)
        return ICRF(position, velocity, t=t, center=self.center)

    def _observe_from_bcrs(self, observer):
        """Observe this target from an observer in BCRS coordinates.

        Implements the Skyfield VectorFunction protocol for observe().
        Delegates to _correct_for_light_travel_time which iterates
        light-time using self._at().

        Args:
            observer: Skyfield Barycentric position of the observer

        Returns:
            Tuple of (position_au, velocity_au_per_d, time, light_time_days)
        """
        from skyfield.vectorlib import _correct_for_light_travel_time

        return _correct_for_light_travel_time(observer, self)

    def __repr__(self):
        return f"<SpkType21Target NAIF={self._naif_id}>"


def get_spk_type21_target(ipl: int, jd_tt: Optional[float] = None):
    """Get a Skyfield-compatible VectorFunction target for a type 21 asteroid.

    Returns a _SpkType21Target that can be used with Skyfield's observe/apparent
    pipeline, or None if the body is not registered or not type 21.

    Args:
        ipl: libephemeris body ID (e.g., CHIRON=15, CERES=17)
        jd_tt: Optional epoch (JD TT).  When given and outside the
            kernel's coverage, returns None so the caller falls through
            to the legacy path, which raises EphemerisRangeError and
            lands on the documented Keplerian out-of-coverage fallback —
            the same flow as type 2/3 kernels.  (A small margin covers
            the light-time retardation at the coverage edges.)

    Returns:
        _SpkType21Target or None
    """
    from . import state

    if ipl not in state._SPK_BODY_MAP:
        return None

    spk_file, naif_id = state._SPK_BODY_MAP[ipl]

    # Only handle type 21
    spk_type = _detect_spk_type(spk_file)
    if spk_type != 21:
        return None

    kernel = _load_type21_kernel(spk_file)
    if kernel is None:
        return None

    if jd_tt is not None:
        coverage = get_spk_coverage(spk_file)
        if coverage is not None:
            start_jd, end_jd = coverage
            # Outside the raw coverage: caller falls through to the direct /
            # Keplerian out-of-coverage path.
            if jd_tt <= start_jd or jd_tt >= end_jd:
                return None
            # Coverage-edge margin: Skyfield's observe() retards the target by
            # the one-way light-time, so jd_tt must sit at least that far inside
            # the kernel or observe() raises instead of falling through to the
            # documented Keplerian path. Size it from the body's actual
            # heliocentric distance — a fixed margin sized for Chiron (~0.11 d
            # near aphelion) is far too small for distant TNOs (Eris/Sedna reach
            # ~0.5 d). Floor at 0.2 d for the inner small bodies.
            _AU_KM = 149597870.7
            _C_AU_DAY = 173.1446326846693
            try:
                pos_km, _ = kernel.compute_type21(10, naif_id, float(jd_tt))
                helio_au = (
                    math.sqrt(pos_km[0] ** 2 + pos_km[1] ** 2 + pos_km[2] ** 2) / _AU_KM
                )
            except Exception:  # noqa: BLE001  # intentional best-effort probe
                # Probe failed near the edge: assume a distant body (~175 AU,
                # ~1 d light-time) so the margin is conservative, not too small.
                # Catch broadly — the probe is best-effort and the underlying
                # compute_type21/get_MDA_record/spke21 also raise RuntimeError
                # (e.g. "Invalid data" / segment errors); a probe failure must
                # never escape this margin sizing and break the documented
                # graceful fall-through to the Keplerian out-of-coverage path.
                helio_au = 175.0
            # Earth's orbit adds at most ~1 AU to the geocentric distance; a
            # 1.2 safety factor absorbs the light-time iteration and slop.
            margin = max(0.2, (helio_au + 1.0) / _C_AU_DAY * 1.2)
            if jd_tt < start_jd + margin or jd_tt > end_jd - margin:
                return None

    # Get the Sun target from the main ephemeris
    planets = state.get_planets()
    sun = planets["sun"]

    return _SpkType21Target(kernel, naif_id, sun)


def _calc_type21_position(
    kernel,
    naif_id: int,
    t,
    iflag: int,
) -> Optional[tuple]:
    """
    Calculate body position using SPK type 21 kernel.

    This function uses the spktype21 library to compute heliocentric
    coordinates, then converts to geocentric apparent positions with:
    - Light-time correction (unless FLG_TRUEPOS)
    - Aberration correction (unless FLG_NOABERR)
    - Precession to equinox of date (unless FLG_J2000)
    - Nutation (unless FLG_NONUT or FLG_J2000)

    Args:
        kernel: SPKType21 object
        naif_id: NAIF ID of the target body (e.g., 20002060 for Chiron)
        t: Skyfield Time object
        iflag: Calculation flags

    Returns:
        Position tuple (lon, lat, dist, speed_lon, speed_lat, speed_dist)
        or None if calculation fails.
    """
    from . import state
    from .constants import (
        FLG_HELCTR,
        FLG_J2000,
        FLG_NOABERR,
        FLG_NONUT,
        FLG_SIDEREAL,
        FLG_SPEED,
        FLG_TRUEPOS,
    )
    from .astrometry import (
        nutation_angles,
        precess_from_j2000,
        apply_aberration_to_position,
    )

    # Constants
    AU_KM = 149597870.7
    C_LIGHT_AU_DAY = 173.1446326846693  # Speed of light in AU/day

    # Check flags
    is_heliocentric = bool(iflag & FLG_HELCTR)
    # Light-time is applied for heliocentric output too (retarded by the
    # Sun->body light time, see the iteration below); the reference and the
    # rest of the library retard heliocentric positions. Aberration, however,
    # is a geocentric-observer effect and is excluded from heliocentric output.
    apply_light_time = not (iflag & FLG_TRUEPOS)
    apply_aberration = not (iflag & FLG_NOABERR) and not is_heliocentric
    apply_precession = not (iflag & FLG_J2000)
    apply_nutation = not (iflag & FLG_NONUT) and apply_precession

    # Earth/Sun sources for geocentric reduction and stellar aberration.
    planets = state.get_planets()
    sun = planets["sun"]
    earth = planets["earth"]

    def _of_date_lon_lat_dist(t_eval):
        """Tropical of-date (lon, lat, dist) for this type-21 body at t_eval.

        Runs the full apparent-place pipeline -- light-time iteration,
        geocentric reduction, annual aberration, precession J2000->of-date, and
        nutation of the equinox -- so the returned longitude/latitude live in
        the SAME of-date frame as this function's position output. Sampling it
        at t +/- dt lets FLG_SPEED central-difference the of-date longitude,
        which (unlike the previous J2000-frame analytic derivative) includes
        the frame-rotation rate -- precession plus nutation of the equinox,
        ~0.2"/day. This mirrors the type-2/3 path in calc_spk_body_position().
        Returns None on kernel-read failure.
        """
        jd_tdb_e = t_eval.tdb
        jd_tt_e = t_eval.tt

        # Earth heliocentric position (ecliptic J2000) for geocentric output.
        earth_ssb = np.array(earth.at(t_eval).position.au)
        sun_ssb = np.array(sun.at(t_eval).position.au)
        earth_helio_ecl = np.array(_icrs_to_ecliptic_j2000(*(earth_ssb - sun_ssb)))

        # Earth barycentric (SSB) velocity (ecliptic J2000) for stellar
        # aberration -- the IAU frame shared with the planets.py pipeline.
        if apply_aberration:
            earth_vel_bary_ecl = np.array(
                _icrs_to_ecliptic_j2000(*np.array(earth.at(t_eval).velocity.au_per_d))
            )
        else:
            earth_vel_bary_ecl = np.array([0.0, 0.0, 0.0])

        # Step 1: heliocentric position from SPK, with light-time iteration.
        # Light-time distance is Sun->body for heliocentric output and
        # Earth->body for geocentric output.
        jd_compute = jd_tdb_e
        if apply_light_time:
            for _ in range(3):
                try:
                    pos_km, _vel_km = kernel.compute_type21(10, naif_id, jd_compute)
                except (OSError, ValueError, KeyError, IndexError, RuntimeError) as e:
                    get_logger().debug("SPK type 21 computation failed: %s", e)
                    return None
                pos_ecl = np.array(_icrs_to_ecliptic_j2000(*(np.array(pos_km) / AU_KM)))
                if is_heliocentric:
                    dist_au = float(np.linalg.norm(pos_ecl))
                else:
                    dist_au = float(np.linalg.norm(pos_ecl - earth_helio_ecl))
                jd_compute = jd_tdb_e - dist_au / C_LIGHT_AU_DAY

        # Position at the (light-time corrected) emission epoch.
        try:
            pos_km, _vel_km = kernel.compute_type21(10, naif_id, jd_compute)
        except (OSError, ValueError, KeyError, IndexError, RuntimeError) as e:
            # RuntimeError included to match the documented graceful
            # fall-through (compute_type21/spke21 raise RuntimeError on
            # segment / "Invalid data" errors).
            get_logger().debug("SPK type 21 computation failed: %s", e)
            return None

        pos_ecl = np.array(_icrs_to_ecliptic_j2000(*(np.array(pos_km) / AU_KM)))

        # Step 2: geocentric reduction (unless heliocentric output requested).
        if is_heliocentric:
            pos_final = pos_ecl
        else:
            pos_final = pos_ecl - earth_helio_ecl

        # Step 3: annual aberration (geocentric apparent place only).
        if apply_aberration:
            pos_final = np.array(
                apply_aberration_to_position(
                    (float(pos_final[0]), float(pos_final[1]), float(pos_final[2])),
                    (
                        float(earth_vel_bary_ecl[0]),
                        float(earth_vel_bary_ecl[1]),
                        float(earth_vel_bary_ecl[2]),
                    ),
                )
            )

        # Step 4: spherical coordinates in the J2000 ecliptic.
        x, y, z = float(pos_final[0]), float(pos_final[1]), float(pos_final[2])
        r = math.sqrt(x**2 + y**2 + z**2)
        lon_j2000 = math.degrees(math.atan2(y, x)) % 360.0
        lat_j2000 = math.degrees(math.asin(z / r)) if r > 0 else 0.0

        # Step 5: precession J2000 -> equinox of date (if requested).
        if apply_precession:
            lon_d, lat_d = precess_from_j2000(lon_j2000, lat_j2000, jd_tt_e)
        else:
            lon_d, lat_d = lon_j2000, lat_j2000

        # Step 6: nutation in longitude (of-date equinox, if requested).
        if apply_nutation:
            delta_psi, _delta_eps = nutation_angles(jd_tt_e)
            # nutation_angles returns values in degrees
            lon_d = lon_d + delta_psi

        return (lon_d % 360.0, lat_d, r)

    base = _of_date_lon_lat_dist(t)
    if base is None:
        return None
    lon, lat, r = base

    # =========================================================================
    # Step 7: Calculate speeds if requested
    # =========================================================================
    # Central finite-difference of the OF-DATE longitude/latitude/distance, so
    # the of-date frame-rotation rate (precession + nutation of the equinox) is
    # included. The previous analytic derivative was evaluated in the fixed
    # J2000 ecliptic and omitted that ~0.2"/day term. Mirrors the type-2/3
    # path in calc_spk_body_position().
    if iflag & FLG_SPEED:
        dt = 1.0 / 86400.0  # 1 second, in days
        ts = get_timescale()
        prev = _of_date_lon_lat_dist(ts.tt_jd(t.tt - dt))
        nxt = _of_date_lon_lat_dist(ts.tt_jd(t.tt + dt))
        if prev is not None and nxt is not None:
            lon_prev, lat_prev, dist_prev = prev
            lon_next, lat_next, dist_next = nxt
            speed_lon = (lon_next - lon_prev) / (2.0 * dt)
            speed_lat = (lat_next - lat_prev) / (2.0 * dt)
            speed_dist = (dist_next - dist_prev) / (2.0 * dt)
            # Unwrap a 0/360 boundary crossing between the two samples.
            if speed_lon > 180.0 / (2.0 * dt):
                speed_lon -= 360.0 / (2.0 * dt)
            elif speed_lon < -180.0 / (2.0 * dt):
                speed_lon += 360.0 / (2.0 * dt)
            # At (numerically) the ecliptic pole the longitude -- atan2 of
            # two ~zero components -- is undefined, so its central difference
            # measures the arbitrary branch picked at each sample, not a
            # physical rate. Report 0.0, matching the xy_sq == 0 guard of the
            # analytic derivative this central difference replaced.
            if max(abs(lat), abs(lat_prev), abs(lat_next)) >= 90.0 - 1e-9:
                speed_lon = 0.0
        else:
            speed_lon, speed_lat, speed_dist = 0.0, 0.0, 0.0
    else:
        speed_lon, speed_lat, speed_dist = 0.0, 0.0, 0.0

    # =========================================================================
    # Step 8: Apply sidereal correction if requested
    # =========================================================================
    if iflag & FLG_SIDEREAL:
        # Flag-aware ayanamsa (true for of-date, mean for NONUT/J2000) and the
        # ayanamsha-rate speed correction, via the canonical helper. The
        # type-21 longitude already honors FLG_J2000/FLG_NONUT (precession /
        # nutation applied above), so the flag-aware ayanamsa matches it; the
        # helper also subtracts the ~50"/yr drift from speed_lon under
        # FLG_SPEED (previously left at the tropical value).
        from .planets import _apply_sidereal_correction

        lon, speed_lon = _apply_sidereal_correction(lon, speed_lon, t.ut1, iflag)

    return (lon, lat, r, speed_lon, speed_lat, speed_dist)


# =============================================================================
# CONVENIENCE FUNCTION
# =============================================================================


def download_and_register_spk(
    body: str,
    ipl: int,
    start: str,
    end: str,
    naif_id: Optional[int] = None,
    path: Optional[str] = None,
    center: str = "500@0",
    overwrite: bool = False,
    timeout: int = 120,
) -> str:
    """
    Download SPK from Horizons and register it for use in calculations.

    Convenience function that combines download_spk() and register_spk_body().

    Args:
        body: Target body identifier (see download_spk)
        ipl: libephemeris body ID (e.g., CHIRON)
        start: Start date (YYYY-MM-DD)
        end: End date (YYYY-MM-DD)
        naif_id: NAIF ID in the kernel. If None, auto-detected from SPK file.
            JPL Horizons uses: asteroid_number + 20000000
        path: Directory to save file (default: ~/.libephemeris/spk/)
        center: Reference center (default: "500@0" = SSB)
        overwrite: Overwrite existing file if True
        timeout: Request timeout in seconds

    Returns:
        str: Path to downloaded SPK file

    Raises:
        ValueError: If body not found, or naif_id cannot be deduced and not provided

    Example:
        >>> download_and_register_spk(
        ...     body="Chiron",
        ...     ipl=CHIRON,
        ...     start="2000-01-01",
        ...     end="2100-01-01",
        ... )
        >>> # Chiron now uses SPK for calculations
    """
    # Download
    spk_path = download_spk(
        body=body,
        start=start,
        end=end,
        path=path,
        center=center,
        overwrite=overwrite,
        timeout=timeout,
    )

    # Deduce NAIF ID if not provided
    if naif_id is None:
        # Try to get asteroid number
        asteroid_number = _extract_asteroid_number(body)
        if asteroid_number is not None:
            # Try to find the NAIF ID in the SPK file (handles both conventions)
            naif_id = _find_naif_id_for_asteroid(spk_path, asteroid_number)

        if naif_id is None:
            # Fall back to legacy convention
            naif_id = _deduce_naif_id(body)

        if naif_id is None:
            # Last resort: scan SPK file for target IDs.
            # If there's exactly one unique target (excluding common centers like
            # Sun=10, SSB=0), use it. This handles name-syntax bodies like "Ceres;"
            # where we can't extract an asteroid number from the name.
            targets = _get_spk_targets(spk_path)
            unique_targets = set(t for t in targets if t not in (0, 10))
            if len(unique_targets) == 1:
                naif_id = unique_targets.pop()

        if naif_id is None:
            raise ValueError(
                f"Cannot deduce NAIF ID for '{body}'. "
                f"Please provide naif_id explicitly. "
                f"For numbered asteroids: naif_id = asteroid_number + 20000000"
            )

    # Register
    register_spk_body(ipl, spk_path, naif_id)

    return spk_path


# =============================================================================
# SPK POSITION CALCULATION
# =============================================================================


def calc_spk_body_position(
    t,
    ipl: int,
    iflag: int,
) -> Optional[tuple[float, float, float, float, float, float]]:
    """
    Calculate body position using SPK kernel.

    Internal function called by planets._calc_body() when SPK is registered.
    Automatically detects SPK type and uses appropriate handler:
    - Type 21: Uses spktype21 library (JPL Horizons asteroids/comets)
    - Type 2/3: Uses Skyfield (Chebyshev polynomials)

    Args:
        t: Skyfield Time object
        ipl: libephemeris body ID
        iflag: Calculation flags (FLG_SPEED, FLG_HELCTR, etc.)

    Returns:
        Position tuple (lon, lat, dist, speed_lon, speed_lat, speed_dist)
        or None if body not registered or outside coverage.

    Raises:
        ValueError: If JD is outside SPK coverage
    """
    from . import state
    from .constants import (
        FLG_HELCTR,
        FLG_SPEED,
        FLG_SIDEREAL,
        FLG_NONUT,
        FLG_J2000,
    )

    # Check if body is registered
    if ipl not in state._SPK_BODY_MAP:
        return None

    spk_file, naif_id = state._SPK_BODY_MAP[ipl]

    # Detect SPK type
    spk_type = _detect_spk_type(spk_file)

    # Handle type 21 (JPL Horizons asteroids/comets)
    if spk_type == 21:
        # Out-of-coverage epochs raise EphemerisRangeError exactly like
        # the type 2/3 branch below; None stays reserved for genuine
        # computation failures (unreadable/corrupt kernel data).
        coverage = get_spk_coverage(spk_file)
        if coverage is not None:
            start_jd, end_jd = coverage
            if t.tt < start_jd or t.tt > end_jd:
                from .exceptions import EphemerisRangeError

                raise EphemerisRangeError(
                    f"JD {t.tt:.1f} outside SPK coverage "
                    f"[{start_jd:.1f}, {end_jd:.1f}] for body {ipl}"
                )
        kernel = _load_type21_kernel(spk_file)
        if kernel is not None:
            result = _calc_type21_position(kernel, naif_id, t, iflag)
            if result is not None:
                return result
        # Fall through to return None if type 21 calculation fails
        return None

    # Type 2/3: Use Skyfield
    kernel = state._SPK_KERNELS.get(spk_file)
    if kernel is None:
        return None

    # Check coverage
    coverage = get_spk_coverage(spk_file)
    if coverage is not None:
        start_jd, end_jd = coverage
        if t.tt < start_jd or t.tt > end_jd:
            from .exceptions import EphemerisRangeError

            raise EphemerisRangeError(
                f"JD {t.tt:.1f} outside SPK coverage [{start_jd:.1f}, {end_jd:.1f}] "
                f"for body {ipl}"
            )

    # Get target from kernel
    try:
        target = kernel[naif_id]
    except KeyError:
        return None

    # Get main ephemeris for Earth/Sun
    planets = state.get_planets()

    # Determine observer
    is_heliocentric = bool(iflag & FLG_HELCTR)

    if is_heliocentric:
        observer = planets["sun"]
    else:
        observer = planets["earth"]

    # Apparent pipeline, mirroring the type-21 path: light-time iteration
    # on the target (unless FLG_TRUEPOS) and, for geocentric output, annual
    # aberration (unless FLG_NOABERR). Heliocentric output is light-time
    # corrected too -- here the observer is the Sun, so the iteration below
    # retards by the Sun->body light time -- but carries no aberration,
    # matching the reference and the rest of the library's heliocentric paths.
    from .constants import FLG_NOABERR, FLG_TRUEPOS
    from .astrometry import apply_aberration_to_position

    C_AU_DAY = 173.144632674240
    apply_light_time = not (iflag & FLG_TRUEPOS)
    apply_aberr = not (iflag & FLG_NOABERR) and not is_heliocentric
    ts = get_timescale()

    def _rel_at(t_obs):
        """Light-time corrected target-minus-observer vector (AU)."""
        obs_at = observer.at(t_obs)
        obs_xyz = obs_at.position.au
        t_emit = t_obs
        rel = None
        for _ in range(3):
            tgt_xyz = target.at(t_emit).position.au
            rel = [tgt_xyz[i] - obs_xyz[i] for i in range(3)]
            if not apply_light_time:
                break
            dist_au = math.sqrt(rel[0] ** 2 + rel[1] ** 2 + rel[2] ** 2)
            t_emit = ts.tt_jd(t_obs.tt - dist_au / C_AU_DAY)
        if apply_aberr:
            obs_vel = obs_at.velocity.au_per_d
            rel = list(
                apply_aberration_to_position(
                    (rel[0], rel[1], rel[2]),
                    (obs_vel[0], obs_vel[1], obs_vel[2]),
                )
            )
        return rel

    from skyfield.positionlib import ICRF

    def _ecl_at(t_obs):
        rel = _rel_at(t_obs)
        rel_pos = ICRF(rel, t=t_obs, center=399)
        ecl = rel_pos.frame_latlon(ecliptic_frame)
        return (
            ecl[1].degrees % 360.0,
            ecl[0].degrees,
            ecl[2].au,
        )

    lon, lat, dist = _ecl_at(t)

    # Calculate speeds if requested
    speed_lon, speed_lat, speed_dist = 0.0, 0.0, 0.0

    if iflag & FLG_SPEED:
        dt = 1.0 / 86400.0  # 1 second in days
        lon_prev, lat_prev, dist_prev = _ecl_at(ts.tt_jd(t.tt - dt))
        lon_next, lat_next, dist_next = _ecl_at(ts.tt_jd(t.tt + dt))

        speed_lon = (lon_next - lon_prev) / (2.0 * dt)
        speed_lat = (lat_next - lat_prev) / (2.0 * dt)
        speed_dist = (dist_next - dist_prev) / (2.0 * dt)

        if speed_lon > 180.0 / (2.0 * dt):
            speed_lon -= 360.0 / (2.0 * dt)
        if speed_lon < -180.0 / (2.0 * dt):
            speed_lon += 360.0 / (2.0 * dt)

    # FLG_NONUT / FLG_J2000: frame_latlon(ecliptic_frame) above produces a TRUE
    # ecliptic-of-date longitude (with nutation in longitude Δψ). The mean
    # ecliptic (NONUT) and the J2000 frame carry no nutation, so strip Δψ here.
    # For J2000 the downstream _maybe_equatorial_convert precesses of-date→J2000
    # but does not remove nutation, so it must be removed first — matching the
    # type-21 path and the reference, whose J2000 ecliptic output is mean.
    if (iflag & FLG_NONUT) or (iflag & FLG_J2000):
        from .cache import get_cached_nutation

        dpsi_rad, _ = get_cached_nutation(t.tt)
        lon = (lon - math.degrees(dpsi_rad)) % 360.0
        if iflag & FLG_SPEED:
            from .planets import _nutation_rate_deg_per_day

            speed_lon -= _nutation_rate_deg_per_day(t.tt)

    # Apply sidereal correction if requested. Routed through the canonical
    # flag-aware helper so the ayanamsa variant matches the longitude frame set
    # above: TRUE ayanamsa (mean + Δψ) cancels Δψ on a true-of-date longitude,
    # while NONUT/J2000 (Δψ already stripped) get the MEAN ayanamsa. The helper
    # also removes the ~50"/yr ayanamsha drift from speed_lon under FLG_SPEED.
    if iflag & FLG_SIDEREAL:
        from .planets import _apply_sidereal_correction

        lon, speed_lon = _apply_sidereal_correction(lon, speed_lon, t.ut1, iflag)

    return (lon, lat, dist, speed_lon, speed_lat, speed_dist)
