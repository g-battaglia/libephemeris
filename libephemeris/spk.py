# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
SPK kernel support for high-precision minor body calculations.

Provenance:
    Kernel bytes and state vectors come from NASA/JPL Horizons and use the
    public NAIF SPK specification. Project-authored code handles HTTPS,
    validation, target registration, frame conversion, and type dispatch.
    Types 2/3 are read through Skyfield; type 21 uses the separately identified
    MIT ``spktype21`` vendor. Downloaded states are never fitted into a hidden
    runtime model.

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
from contextlib import contextmanager
from typing import Optional, Union

import certifi
import numpy as np

from skyfield.framelib import ecliptic_frame

from .download import _is_valid_bsp
from .exceptions import SPKNotFoundError
from .exotic_bodies import EXOTIC_IDS as _EXOTIC_IDS
from .exotic_bodies import exotic_display_name as _exotic_display_name
from .logging_config import format_file_size, get_logger
from .net import HTTPError, Request, URLError, open_url, urlencode
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

# Kilometres per astronomical unit (IAU 2012 Resolution B2) and the speed of
# light in AU/day (c = 299792.458 km/s, CODATA), shared by the type-21 state
# conversion and the light-time edge sizing.
_AU_KM = 149597870.7
_C_AU_DAY = 173.1446326846693

# Center assumed for a type-21 kernel whose reader exposes no segment
# summaries: JPL Horizons small-body kernels are heliocentric.
_TYPE21_DEFAULT_CENTER = 10

# Upper bounds (AU) on the Earth-to-center distance for the centers a type-21
# kernel may declare. Earth's aphelion is 1.017 AU and the Sun stays within
# 0.02 AU of the barycenter (NASA planetary fact sheet), so 1.1 AU bounds both
# the Earth-to-Sun and the Earth-to-barycenter distance. A planetary-system
# center (barycenter 1-9, or a body of that system such as 399 or 301) adds
# that system's aphelion, rounded up to the next 0.1 AU (same source).
_EARTH_REACH_FROM_SUN_OR_SSB_AU = 1.1
_PLANET_SYSTEM_APHELION_AU = {
    1: 0.5,
    2: 0.8,
    3: 1.1,
    4: 1.7,
    5: 5.5,
    6: 10.2,
    7: 20.2,
    8: 30.4,
    9: 49.4,
}

# Offset (days) of the light-time distance probe from the first segment
# start. The reader converts JD to seconds past J2000; probing a hair inside
# the segment keeps that round trip clear of the inclusive lower bound.
_TYPE21_PROBE_OFFSET_DAYS = 1.0e-6

# Half-step (days) of the FLG_SPEED central difference in the direct type-21
# path (_calc_type21_position). The planet pipeline's own half-steps live in
# planets.py; _type21_speed_stencil_days() takes the largest of them all.
_TYPE21_LEGACY_SPEED_HALF_STEP_DAYS = 1.0 / 86400.0

# Tolerances (days) when deciding whether a cached kernel covers a requested
# date range. The reported start is the usable start, one light-time plus
# the speed stencil inside the stored span: 1.5 days absorbs that band for
# any body within ~250 AU. The end tolerance covers the stencil and the
# sub-day rounding of Horizons coverage boundaries.
_SPK_COVERAGE_START_TOLERANCE_DAYS = 1.5
_SPK_COVERAGE_END_TOLERANCE_DAYS = 0.5


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

    Only a permanent catalog number is returned: a bare number, a
    parenthesized number, or a number followed by a proper name.  A
    provisional designation ("2003 UB313", "1998 KY26") begins with a
    discovery *year*, and a cometary designation ("1P/Halley", "C/1995 O1")
    carries an orbit-type letter and slash; neither leading token is a
    catalog number, so both return None.  Treating the year/comet prefix
    as a catalog number would fabricate a NAIF ID that is absent from the
    downloaded kernel and never resolves.

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
        >>> _extract_asteroid_number("2003 UB313") is None
        True
        >>> _extract_asteroid_number("1P/Halley") is None
        True
    """
    text = body.strip()

    # Cometary designations ("1P/Halley", "C/1995 O1", "2P/Encke") carry a
    # slash separating the orbit-type letter from the name/designation and
    # never denote a numbered-asteroid catalog number.
    if "/" in text:
        return None

    # Provisional designations ("2003 UB313", "1998 KY26", "2014 MU69"):
    # a four-digit discovery year followed by a packed survey code of one
    # or two uppercase letters plus an optional cycle count.  The leading
    # year is not a catalog number.  A genuine "number name" pair such as
    # "2003 Harding" is NOT matched here because the name continues with
    # lowercase letters rather than terminating after the uppercase code.
    if re.match(r"^\(?\d{4}\)?\s+[A-Z]{1,2}\d*\b", text):
        return None

    # Numbered asteroid: bare number, "(number)", or "number name".
    match = re.match(r"^\(?(\d+)\)?", text)
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

# Fill in every remaining precomputed exotic minor body from the single-source
# registry so SPK diagnostics name the full set (Interamnia, Davida, Europa,
# Sylvia, Psyche, Sappho, Pandora, Lilith, Amor, Itokawa, Ryugu, Toro,
# Toutatis, Icarus, Hidalgo) without hand-maintaining a second list. The
# hard-coded entries above win on any overlap.
for _exotic_id in _EXOTIC_IDS:
    _BODY_NAMES.setdefault(_exotic_id, _exotic_display_name(_exotic_id))
del _exotic_id


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
    # URL-encode parameters through the central network module.
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

    query = urlencode(params, safe="'@")
    return f"{HORIZONS_API_URL}?{query}"


def _sanitize_filename(body: Union[int, str]) -> str:
    """
    Create a safe filename stem from a body identifier.

    This is the canonical sanitizer shared by the SPK writer
    (:func:`download_spk`) and every cache-lookup path in
    :mod:`libephemeris.spk_auto` (``get_cache_path``, ``is_spk_cached``,
    ``_find_covering_spk``, ``_generate_spk_cache_filename``).  All of
    them must agree byte-for-byte so that a file written for ``"Ceres;"``
    (stem ``"ceres"``) is found again on lookup; a divergent sanitizer
    previously left trailing underscores (``"ceres_"``) and produced
    permanent cache misses plus redundant re-downloads.

    Args:
        body: Body name or numeric identifier (int or str).

    Returns:
        Sanitized lowercase string safe for use in filenames.
    """
    # Remove/replace unsafe characters
    safe = re.sub(r"[^\w\-]", "_", str(body).lower())
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
        path: Directory to save the file. If None, uses the configured SPK
            cache directory (set_spk_cache_dir() / LIBEPHEMERIS_SPK_DIR /
            TOML ``spk_dir``) when one is set, otherwise ``<data dir>/spk``.
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
    # Determine output directory. An explicit set_spk_cache_dir() /
    # LIBEPHEMERIS_SPK_DIR override wins, because that is the only directory
    # every spk_auto lookup consults; without one the file keeps landing in
    # <data dir>/spk, which in the default configuration is the same path the
    # lookups use. Honouring the override is what issue #55 asked for; moving
    # the cache anywhere else would hide kernels users already downloaded.
    if path is None:
        from .state import get_spk_cache_dir

        path = get_spk_cache_dir()
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

    # Check if already exists and valid. A corrupt cached kernel is NOT
    # removed here: the download below validates its temp file and only then
    # replaces this path, so a failed re-download leaves the existing bytes
    # untouched rather than deleting the user's only copy.
    if os.path.exists(filepath) and not overwrite:
        if _is_valid_bsp(filepath):
            return filepath
        logger.warning("Cached SPK file %s is corrupted, re-downloading", filepath)

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
            req = Request(url)
            req.add_header("User-Agent", "libephemeris/0.1")

            with open_url(
                req,
                purpose=f"Horizons SPK generation for {body}",
                timeout=timeout,
                context=ssl_context,
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

        except HTTPError as e:
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

        except URLError as e:
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


def _date_to_jd(date_s: str) -> float:
    """Convert a YYYY-MM-DD date to a Julian day at 00:00 UT."""
    from .time_utils import julday

    year_s, month_s, day_s = date_s.split("-")
    return julday(int(year_s), int(month_s), int(day_s), 0.0)


def _spk_covers_dates(spk_path: str, start: str, end: str) -> bool:
    """Return True when an SPK file covers the requested calendar interval."""
    coverage = get_spk_coverage(spk_path)
    if coverage is None:
        return True

    start_jd, end_jd = _date_to_jd(start), _date_to_jd(end)
    cov_start, cov_end = coverage
    # The coverage is the usable span: its start sits one light-time inside
    # the stored data, and Horizons reports boundaries at sub-day precision,
    # so allow the documented tolerances around the request boundaries.
    return (
        cov_start <= start_jd + _SPK_COVERAGE_START_TOLERANCE_DAYS
        and cov_end >= end_jd - _SPK_COVERAGE_END_TOLERANCE_DAYS
    )


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
        ValueError: If naif_id not found in SPK kernel, or (type 21) if the
            target's segments declare more than one center or a center the
            base ephemeris cannot resolve

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
        # Validate the target and its center against the reader's segment
        # summaries, as the type 2/3 branch validates against Skyfield. An
        # unknown NAIF ID or an unresolvable center would otherwise fail on
        # every evaluation; registration is where that belongs. A reader
        # that exposes no summaries cannot be checked and registers as is.
        if getattr(kernel, "segments", None):
            if not _type21_segments_for(kernel, naif_id):
                available = _type21_targets(kernel)
                raise ValueError(
                    f"NAIF ID {naif_id} has no type 21 segment in SPK kernel "
                    f"{spk_file}. Available type 21 targets (first 10): "
                    f"{available[:10]}"
                )
            _validate_type21_center(_type21_center(kernel, naif_id), naif_id, spk_file)
    else:
        # Type 2/3: Use Skyfield
        state._load_spk_kernel(spk_file)

        # Validate that naif_id exists in kernel
        kernel = state._SPK_KERNELS.get(spk_file)
        if kernel is not None:
            # Check if target exists in kernel. Skyfield raises KeyError for
            # unknown integer keys but ValueError ("unknown SPICE target") for
            # e.g. Sun-centered type-2 kernels, so catch both.
            try:
                _ = kernel[naif_id]
            except (KeyError, ValueError):
                # Try with string name
                try:
                    _ = kernel[str(naif_id)]
                except (KeyError, ValueError):
                    available = []
                    if hasattr(kernel, "names"):
                        available = list(kernel.names())[:10]
                    raise ValueError(
                        f"NAIF ID {naif_id} not found in SPK kernel {spk_file}. "
                        f"Available targets (first 10): {available}"
                    ) from None

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
    Get the usable time coverage of an SPK kernel.

    Supports both type 2/3 (Skyfield) and type 21 (spktype21) kernels. The
    bounds are Julian Days in TDB, the scale the kernels are evaluated in.

    For a type 21 kernel the span is the one every epoch of which can be
    served by the apparent-place pipeline: the stored span minus the
    light-time band at its start (the pipeline retards the evaluation epoch
    by the Earth-to-body light-time) and minus the FLG_SPEED stencil at both
    ends. An epoch inside the returned span is never answered from a
    lower-precision source.

    Args:
        spk_file: Path to SPK file

    Returns:
        Tuple of (start_jd, end_jd) Julian Day (TDB) coverage, or None if
        the kernel exposes no readable coverage.
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
        # Usable span across the file's type-21 targets, read through the
        # reader that evaluates them (registered kernels come from the cache;
        # a scanned file is opened for the probe and closed again).
        try:
            with _type21_kernel_for_probe(spk_file) as kernel:
                if kernel is None:
                    return None
                spans = [
                    span
                    for span in (
                        _type21_usable_span(kernel, target)
                        for target in _type21_targets(kernel)
                    )
                    if span is not None
                ]
        except (OSError, ValueError, KeyError, IndexError, RuntimeError):
            return None
        if not spans:
            return None
        return (min(span[0] for span in spans), max(span[1] for span in spans))

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
# SPK TYPE 21: SEGMENT METADATA, CENTER RESOLUTION, USABLE COVERAGE
# =============================================================================


def _type21_segments_for(kernel, naif_id: int) -> list:
    """Return the type-21 segments of ``kernel`` whose target is ``naif_id``.

    The reader exposes its segment summaries (``target``, ``center``,
    ``data_type``, ``start_jd``, ``end_jd``), so a kernel can be checked
    against the id it is registered for and its own center can be read
    instead of assumed. A reader without summaries yields an empty list.
    """
    segments = getattr(kernel, "segments", None) or []
    return [
        segment
        for segment in segments
        if getattr(segment, "target", None) == naif_id
        and getattr(segment, "data_type", None) == 21
    ]


def _type21_targets(kernel) -> list[int]:
    """Return the sorted NAIF ids that have type-21 segments in ``kernel``."""
    segments = getattr(kernel, "segments", None) or []
    return sorted(
        {
            int(segment.target)
            for segment in segments
            if getattr(segment, "data_type", None) == 21
            and getattr(segment, "target", None) is not None
        }
    )


def _type21_center(kernel, naif_id: int) -> int:
    """Return the center declared by ``naif_id``'s type-21 segments.

    Falls back to the heliocentric default when the reader exposes no
    summaries. The reader serves one (target, center) pair per evaluation,
    so segments declaring different centers for the same target are
    rejected here rather than half-served at runtime.

    Raises:
        ValueError: If the target's segments declare more than one center.
    """
    centers = {
        int(segment.center)
        for segment in _type21_segments_for(kernel, naif_id)
        if getattr(segment, "center", None) is not None
    }
    if not centers:
        return _TYPE21_DEFAULT_CENTER
    if len(centers) > 1:
        raise ValueError(
            f"NAIF ID {naif_id} has type 21 segments with different centers "
            f"{sorted(centers)}; a single center per target is required"
        )
    return centers.pop()


def _type21_earth_reach_au(center: int) -> Optional[float]:
    """Upper bound (AU) on the distance between Earth and ``center``.

    Returns None for a center that is neither the Sun, the barycenter nor a
    planetary-system id, i.e. one the base ephemeris cannot resolve.
    """
    if center in (0, 10):
        return _EARTH_REACH_FROM_SUN_OR_SSB_AU
    if center < 0:
        return None
    system = center if center < 10 else center // 100
    aphelion = _PLANET_SYSTEM_APHELION_AU.get(system)
    if aphelion is None:
        return None
    return aphelion + _EARTH_REACH_FROM_SUN_OR_SSB_AU


def _type21_center_target(center: int, planets):
    """Return the base-ephemeris VectorFunction of ``center`` (None for the SSB).

    Center 0 is the barycenter itself, so its state is identically zero. Any
    other center is looked up by NAIF id in the base ephemeris; the lookup
    errors propagate to the caller, which decides how to report them.
    """
    if center == 0:
        return None
    return planets[center]


def _validate_type21_center(center: int, naif_id: int, spk_file: str) -> None:
    """Reject, at registration, a type-21 center the runtime chain cannot use.

    The Sun and the barycenter need no lookup. Any other center must be a
    planetary-system id and must resolve in the base ephemeris now, so the
    failure is a registration error rather than a per-evaluation one.

    Raises:
        ValueError: If the center is unknown or absent from the base ephemeris.
    """
    if _type21_earth_reach_au(center) is None:
        raise ValueError(
            f"Center {center} of NAIF ID {naif_id} in SPK kernel {spk_file} is "
            "neither the Sun (10), the solar-system barycenter (0) nor a "
            "planetary-system id; it cannot be resolved against the base "
            "ephemeris"
        )
    if center in (0, 10):
        return
    from . import state

    try:
        _type21_center_target(center, state.get_planets())
    except (KeyError, ValueError, IndexError, TypeError, RuntimeError, OSError) as exc:
        raise ValueError(
            f"Center {center} of NAIF ID {naif_id} in SPK kernel {spk_file} is "
            f"not available in the base ephemeris: {exc}"
        ) from exc


def _type21_speed_stencil_days() -> float:
    """Largest half-step any FLG_SPEED central difference applies to a type-21 body."""
    from .planets import (
        _ASTEROID_TRUEPOS_SPEED_HALF_STEP_DAYS,
        _BODY_SPEED_HALF_STEP_DAYS,
    )

    return max(
        _ASTEROID_TRUEPOS_SPEED_HALF_STEP_DAYS,
        _BODY_SPEED_HALF_STEP_DAYS,
        _TYPE21_LEGACY_SPEED_HALF_STEP_DAYS,
    )


def _type21_request_stencil_days(iflag: int) -> float:
    """Half-step the planet pipeline samples around a FLG_SPEED request.

    Zero without FLG_SPEED (the stencil samples themselves arrive without
    the flag); otherwise the pipeline's own half-step for a type-21 body:
    five minutes for the geometric place, one second for the apparent one.
    """
    from .constants import FLG_SPEED, FLG_TRUEPOS

    if not iflag & FLG_SPEED:
        return 0.0
    from .planets import (
        _ASTEROID_TRUEPOS_SPEED_HALF_STEP_DAYS,
        _BODY_SPEED_HALF_STEP_DAYS,
    )

    if iflag & FLG_TRUEPOS:
        return _ASTEROID_TRUEPOS_SPEED_HALF_STEP_DAYS
    return _BODY_SPEED_HALF_STEP_DAYS


def _type21_stored_span(kernel, naif_id: int) -> Optional[tuple[float, float]]:
    """Return the TDB span ``[start, end)`` stored for ``naif_id``, or None.

    The reader evaluates a target over the union of its segments, from the
    earliest start (inclusive) to the latest end (exclusive).
    """
    segments = _type21_segments_for(kernel, naif_id)
    starts = [float(seg.start_jd) for seg in segments if hasattr(seg, "start_jd")]
    ends = [float(seg.end_jd) for seg in segments if hasattr(seg, "end_jd")]
    if not starts or not ends:
        return None
    return min(starts), max(ends)


def _type21_usable_span(
    kernel,
    naif_id: int,
    center: Optional[int] = None,
    *,
    stencil_days: Optional[float] = None,
) -> Optional[tuple[float, float]]:
    """Return the TDB span ``[usable_start, usable_end)`` servable for ``naif_id``.

    The apparent-place pipeline never evaluates the kernel at the requested
    epoch alone: light-time retardation moves the evaluation epoch backwards
    by up to one observer-to-body light-time, and FLG_SPEED samples a
    half-step on either side. The stored span therefore has an unusable band
    at its start (light-time plus stencil), while at its end only the stencil
    is lost, because retardation moves away from the end. The light-time
    bound is the body's distance from its center at the segment start plus
    the Earth-to-center bound for that center, divided by c.

    ``stencil_days`` selects the stencil to subtract. The default (None) is
    the largest half-step any path applies, which is what the reported
    coverage and the direct path (it samples its stencil internally) use.
    The planet pipeline reaches this module once per request AND once per
    stencil sample, so its gate passes the request's own half-step -- zero
    for a sample, which arrives without FLG_SPEED -- and never subtracts
    the stencil twice.

    Returns None when the kernel exposes no segment metadata for the target,
    when its center cannot be bounded, or when the distance probe fails.
    """
    stored = _type21_stored_span(kernel, naif_id)
    if stored is None:
        return None
    start_jd, end_jd = stored
    if center is None:
        center = _type21_center(kernel, naif_id)
    reach_au = _type21_earth_reach_au(center)
    if reach_au is None:
        return None
    try:
        pos_km, _vel_km_s = kernel.compute_type21(
            center, naif_id, start_jd + _TYPE21_PROBE_OFFSET_DAYS
        )
    except (OSError, ValueError, KeyError, IndexError, RuntimeError):
        return None
    distance_au = (
        math.sqrt(float(pos_km[0]) ** 2 + float(pos_km[1]) ** 2 + float(pos_km[2]) ** 2)
        / _AU_KM
    )
    stencil = _type21_speed_stencil_days() if stencil_days is None else stencil_days
    usable_start = start_jd + (distance_au + reach_au) / _C_AU_DAY + stencil
    usable_end = end_jd - stencil
    return usable_start, usable_end


@contextmanager
def _type21_kernel_for_probe(spk_file: str):
    """Yield a type-21 reader for ``spk_file`` without caching an extra handle.

    A registered kernel comes from the cache; an unregistered file (a cache
    directory scan) is opened for the duration of the block and closed
    again. Yields None when the reader module is unavailable.
    """
    from . import state

    kernel = state._SPK_TYPE21_KERNELS.get(spk_file)
    if kernel is not None:
        yield kernel
        return
    reader_cls = _get_spktype21()
    if reader_cls is None:
        yield None
        return
    kernel = reader_cls.open(spk_file)
    try:
        yield kernel
    finally:
        kernel.close()


# =============================================================================
# SKYFIELD VECTORFUNCTION WRAPPER FOR SPK TYPE 21
# =============================================================================


class _SpkType21Target:
    """Skyfield VectorFunction-compatible wrapper for an SPK type-21 body.

    The reader returns the body's state relative to the center its segments
    declare. This wrapper adds the barycentric state of that center -- zero
    for the barycenter itself, the base-ephemeris Sun for center 10, the
    base-ephemeris body for any other center -- so the result is an
    SSB-centered ICRS state compatible with Skyfield's observe()/apparent()
    pipeline. Reader failures surface as SPKEvaluationError: a kernel that
    cannot serve an epoch inside its usable coverage is never answered from
    a lower-precision source on the caller's behalf.

    Routing the body through the planet pipeline keeps the same ICRS ->
    frame_latlon reduction as the planets, avoiding the ~0.3" systematic
    error of the manual ecliptic J2000 + precession/nutation approach.
    """

    def __init__(
        self,
        kernel,
        naif_id: int,
        center_target=None,
        center: Optional[int] = None,
        *,
        spk_file: Optional[str] = None,
        body_id: Optional[int] = None,
    ):
        """Initialize with a reader and the base-ephemeris target of its center.

        Args:
            kernel: SPKType21 kernel object
            naif_id: NAIF ID of the body in the kernel
            center_target: Skyfield VectorFunction (SSB-relative) of the
                center the segments declare; None when that center is the
                barycenter itself.
            center: NAIF id of the declared center. Read from the kernel
                when omitted; heliocentric when the reader has no summaries.
            spk_file: Kernel path, for error messages.
            body_id: libephemeris body id, for error messages.

        Raises:
            ValueError: If a non-barycentric center comes without its target.
        """
        self._kernel = kernel
        self._naif_id = naif_id
        self._segment_center = (
            center if center is not None else _type21_center(kernel, naif_id)
        )
        if self._segment_center != 0 and center_target is None:
            raise ValueError(
                f"center_target is required for type 21 center {self._segment_center}"
            )
        self._center_target = None if self._segment_center == 0 else center_target
        self._spk_file = spk_file
        self._body_id = body_id
        self.center = 0  # SSB-centered (required for observe())
        self.target = naif_id

    @property
    def segment_center(self) -> int:
        """NAIF id of the center the kernel segments declare."""
        return self._segment_center

    def relative_state_au(
        self, jd_tdb: float, jd_fraction: float = 0.0
    ) -> tuple[np.ndarray, np.ndarray]:
        """Body state relative to the kernel center (ICRS, AU and AU/day).

        Args:
            jd_tdb: Evaluation epoch, Julian Day TDB (the reader's scale).
            jd_fraction: Optional second part of the epoch, added to
                ``jd_tdb`` by the reader. Passing a Skyfield time as
                ``(whole, tdb_fraction)`` keeps sub-ULP spacing, which a
                one-second speed stencil needs.

        Raises:
            SPKEvaluationError: If the reader cannot serve the epoch.
        """
        try:
            pos_km, vel_km_s = self._kernel.compute_type21(
                self._segment_center, self._naif_id, float(jd_tdb), float(jd_fraction)
            )
        except (OSError, ValueError, KeyError, IndexError, RuntimeError) as exc:
            raise self._evaluation_error(
                float(jd_tdb) + float(jd_fraction), exc
            ) from exc
        position = np.asarray(pos_km, dtype=float) / _AU_KM
        velocity = np.asarray(vel_km_s, dtype=float) * (86400.0 / _AU_KM)
        return position, velocity

    def state_at(self, t) -> tuple[np.ndarray, np.ndarray]:
        """SSB-centered ICRS state (AU, AU/day) at a scalar Skyfield time.

        The center chain: state relative to the declared center, plus the
        barycentric state of that center from the base ephemeris.
        """
        rel_pos, rel_vel = self.relative_state_au(float(t.whole), float(t.tdb_fraction))
        if self._center_target is None:
            return rel_pos, rel_vel
        center = self._center_target.at(t)
        return (
            np.asarray(center.position.au, dtype=float) + rel_pos,
            np.asarray(center.velocity.au_per_d, dtype=float) + rel_vel,
        )

    def _evaluation_error(self, jd_tdb: float, exc: BaseException):
        """Build the typed error for a reader failure at ``jd_tdb``."""
        from .exceptions import SPKEvaluationError

        stored = _type21_stored_span(self._kernel, self._naif_id)
        span_note = (
            f" (stored segments cover [{stored[0]:.6f}, {stored[1]:.6f}) TDB)"
            if stored is not None
            else ""
        )
        where = f" {self._spk_file}" if self._spk_file else ""
        return SPKEvaluationError(
            f"SPK type 21 kernel{where} could not evaluate NAIF ID "
            f"{self._naif_id} (center {self._segment_center}) at JD "
            f"{jd_tdb:.6f} TDB{span_note}: {exc}. The position is not degraded "
            "to the Keplerian approximation; register a kernel that serves "
            "this epoch or unregister the body.",
            body_id=self._body_id,
            naif_id=self._naif_id,
            spk_file=self._spk_file,
            requested_jd=jd_tdb,
        )

    def _at(self, t):
        """Compute SSB-centered ICRS state at time t.

        Matches the Skyfield VectorFunction._at() protocol:
        Returns (position_au, velocity_au_per_d, gcrs_position, message).

        Args:
            t: Skyfield Time object (scalar or array)

        Returns:
            Tuple of (position, velocity, None, None)
        """
        jd_tdb = t.tdb
        if np.ndim(jd_tdb) == 0:
            position, velocity = self.state_at(t)
            return position, velocity, None, None

        n = len(jd_tdb)
        whole = np.asarray(t.whole, dtype=float)
        fraction = np.asarray(t.tdb_fraction, dtype=float)
        position = np.empty((3, n))
        velocity = np.empty((3, n))
        for i in range(n):
            position[:, i], velocity[:, i] = self.relative_state_au(
                float(whole[i]), float(fraction[i])
            )
        if self._center_target is not None:
            center = self._center_target.at(t)
            position += np.asarray(center.position.au, dtype=float)
            velocity += np.asarray(center.velocity.au_per_d, dtype=float)
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
        return f"<SpkType21Target NAIF={self._naif_id} center={self._segment_center}>"


def _type21_declared_center(
    kernel, naif_id: int, spk_file: Optional[str], body_id: Optional[int]
) -> int:
    """Read the declared center of ``naif_id`` for the runtime paths.

    Raises:
        SPKEvaluationError: If the segments declare more than one center.
    """
    try:
        return _type21_center(kernel, naif_id)
    except ValueError as exc:
        from .exceptions import SPKEvaluationError

        raise SPKEvaluationError(
            str(exc), body_id=body_id, naif_id=naif_id, spk_file=spk_file
        ) from exc


def _type21_target_for(
    kernel,
    naif_id: int,
    spk_file: Optional[str] = None,
    body_id: Optional[int] = None,
    center: Optional[int] = None,
) -> _SpkType21Target:
    """Build the SSB-centered wrapper for ``naif_id`` with its center resolved.

    Raises:
        SPKEvaluationError: If the declared center is ambiguous or not
            available in the base ephemeris.
    """
    from . import state

    if center is None:
        center = _type21_declared_center(kernel, naif_id, spk_file, body_id)
    center_target = None
    if center != 0:
        try:
            center_target = _type21_center_target(center, state.get_planets())
        except (KeyError, ValueError, IndexError, TypeError) as exc:
            from .exceptions import SPKEvaluationError

            where = f" in {spk_file}" if spk_file else ""
            raise SPKEvaluationError(
                f"Center {center} declared by the SPK type 21 segments of NAIF "
                f"ID {naif_id}{where} is not available in the base ephemeris: "
                f"{exc}",
                body_id=body_id,
                naif_id=naif_id,
                spk_file=spk_file,
            ) from exc
    return _SpkType21Target(
        kernel, naif_id, center_target, center, spk_file=spk_file, body_id=body_id
    )


def get_spk_type21_target(ipl: int, jd_tt: Optional[float] = None, iflag: int = 0):
    """Get a Skyfield-compatible VectorFunction target for a type 21 asteroid.

    Returns a _SpkType21Target that can be used with Skyfield's observe/apparent
    pipeline, or None if the body is not registered or not type 21.

    Args:
        ipl: libephemeris body ID (e.g., CHIRON=15, CERES=17)
        jd_tt: Optional epoch (JD TT). It is converted to TDB, the reader's
            scale, and checked against the span the kernel can serve for this
            request: the stored span minus the light-time band at its start
            and minus the request's own FLG_SPEED half-step at both ends
            (see get_spk_coverage). Outside it, returns None so the caller
            falls through to calc_spk_body_position, which raises
            EphemerisRangeError and lands on the documented out-of-coverage
            handling -- the same flow as type 2/3 kernels.
        iflag: Calculation flags of the request; only FLG_SPEED and
            FLG_TRUEPOS matter (they size the stencil).

    Returns:
        _SpkType21Target or None

    Raises:
        SPKEvaluationError: If the kernel's declared center is ambiguous or
            cannot be resolved in the base ephemeris.
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

    center = _type21_declared_center(kernel, naif_id, spk_file, ipl)

    if jd_tt is not None:
        # Gate with the request's own stencil: the planet pipeline calls back
        # here for each FLG_SPEED sample (without the flag), so a request is
        # accepted exactly when its samples will be (see _type21_usable_span).
        span = _type21_usable_span(
            kernel, naif_id, center, stencil_days=_type21_request_stencil_days(iflag)
        )
        if span is not None:
            usable_start, usable_end = span
            jd_tdb = float(get_timescale().tt_jd(jd_tt).tdb)
            if jd_tdb < usable_start or jd_tdb >= usable_end:
                return None

    return _type21_target_for(kernel, naif_id, spk_file, ipl, center=center)


def _calc_type21_position(
    kernel,
    naif_id: int,
    t,
    iflag: int,
    target: Optional[_SpkType21Target] = None,
) -> tuple:
    """
    Calculate body position using SPK type 21 kernel.

    The reader's state (relative to the center its segments declare) is
    completed to a barycentric ICRS state by the center chain of
    :class:`_SpkType21Target`, then reduced to the requested observer with:
    - Light-time correction (unless FLG_TRUEPOS)
    - Aberration correction (unless FLG_NOABERR)
    - Precession to equinox of date (unless FLG_J2000)
    - Nutation (unless FLG_NONUT or FLG_J2000)

    Args:
        kernel: SPKType21 object
        naif_id: NAIF ID of the target body (e.g., 20002060 for Chiron)
        t: Skyfield Time object
        iflag: Calculation flags
        target: Prebuilt wrapper for ``naif_id``; built here when omitted.

    Returns:
        Position tuple (lon, lat, dist, speed_lon, speed_lat, speed_dist)

    Raises:
        SPKEvaluationError: If the kernel cannot serve an evaluation epoch or
            its declared center cannot be resolved in the base ephemeris.
    """
    from . import state
    from .constants import (
        FLG_EQUATORIAL,
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

    if target is None:
        target = _type21_target_for(kernel, naif_id)

    # Check flags
    is_heliocentric = bool(iflag & FLG_HELCTR)
    # Light-time is applied for heliocentric output too (retarded by the
    # Sun->body light time, see the iteration below); the rest of the library
    # retards heliocentric positions as well. Aberration, however, is a
    # geocentric-observer effect and is excluded from heliocentric output.
    apply_light_time = not (iflag & FLG_TRUEPOS)
    # FLG_TRUEPOS requests the geometric place: it suppresses aberration too,
    # exactly like the Keplerian path (planets.py) and planetary_moons.py.
    apply_aberration = (
        not (iflag & FLG_TRUEPOS) and not (iflag & FLG_NOABERR) and not is_heliocentric
    )
    apply_precession = not (iflag & FLG_J2000)
    # SIDEREAL|EQUATORIAL lands on the mean equator of date with no
    # ayanamsha (the shared _sid_eq convention), so nutation is skipped for
    # that combination exactly like FLG_NONUT.
    _sid_eq = bool(iflag & FLG_SIDEREAL) and bool(iflag & FLG_EQUATORIAL)
    apply_nutation = not (iflag & FLG_NONUT) and not _sid_eq and apply_precession

    # Observer for the geocentric / heliocentric reduction. Both are
    # barycentric states, like the target chain, so the subtraction is exact
    # whatever center the kernel declares.
    planets = state.get_planets()
    observer = planets["sun"] if is_heliocentric else planets["earth"]
    ts = get_timescale()

    def _of_date_lon_lat_dist(t_eval):
        """Tropical of-date (lon, lat, dist) for this type-21 body at t_eval.

        Runs the full apparent-place pipeline -- light-time iteration,
        observer reduction, annual aberration, precession J2000->of-date, and
        nutation of the equinox -- so the returned longitude/latitude live in
        the SAME of-date frame as this function's position output. Sampling it
        at t +/- dt lets FLG_SPEED central-difference the of-date longitude,
        which includes the frame-rotation rate -- precession plus nutation of
        the equinox, ~0.2"/day. This mirrors the type-2/3 path in
        calc_spk_body_position(). Reader failures raise SPKEvaluationError.
        """
        jd_tt_e = t_eval.tt

        # Observer barycentric state at the observation epoch (ICRS, AU).
        obs_at = observer.at(t_eval)
        obs_pos = np.asarray(obs_at.position.au, dtype=float)

        # Step 1: barycentric target state at the emission epoch. The
        # light-time iteration retards the epoch by the observer->body
        # distance (Sun->body for heliocentric output, the observer being
        # the Sun). The reader is evaluated in TDB, its own scale.
        # The retarded epoch keeps Skyfield's two-part form (whole, fraction)
        # so the speed stencil's one-second spacing survives the retardation.
        t_emit = t_eval
        if apply_light_time:
            for _ in range(3):
                tgt_pos, _tgt_vel = target.state_at(t_emit)
                dist_au = float(np.linalg.norm(tgt_pos - obs_pos))
                t_emit = ts.tdb_jd(
                    t_eval.whole, t_eval.tdb_fraction - dist_au / _C_AU_DAY
                )
        tgt_pos, _tgt_vel = target.state_at(t_emit)

        # Step 2: observer-relative vector in the J2000 ecliptic.
        pos_final = np.asarray(
            _icrs_to_ecliptic_j2000(*(tgt_pos - obs_pos)), dtype=float
        )

        # Step 3: annual aberration (geocentric apparent place only), from
        # the observer's barycentric velocity in the same ecliptic frame.
        if apply_aberration:
            obs_vel_ecl = _icrs_to_ecliptic_j2000(
                *np.asarray(obs_at.velocity.au_per_d, dtype=float)
            )
            pos_final = np.array(
                apply_aberration_to_position(
                    (float(pos_final[0]), float(pos_final[1]), float(pos_final[2])),
                    (
                        float(obs_vel_ecl[0]),
                        float(obs_vel_ecl[1]),
                        float(obs_vel_ecl[2]),
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

    lon, lat, r = _of_date_lon_lat_dist(t)

    # =========================================================================
    # Step 7: Calculate speeds if requested
    # =========================================================================
    # Central finite-difference of the OF-DATE longitude/latitude/distance, so
    # the of-date frame-rotation rate (precession + nutation of the equinox) is
    # included. Mirrors the type-2/3 path in calc_spk_body_position(). The
    # samples stay inside the kernel because the usable coverage subtracts
    # the stencil at both ends.
    if iflag & FLG_SPEED:
        dt = _TYPE21_LEGACY_SPEED_HALF_STEP_DAYS
        # Two-argument tt_jd(whole, fraction): collapsing t.tt +/- dt into one
        # float64 quantizes the step to the JD ULP (~4.6e-10 d), a ~4e-5
        # relative error on a one-second step; keeping dt in the fraction
        # slot preserves the exact spacing (as the planet pipeline does).
        lon_prev, lat_prev, dist_prev = _of_date_lon_lat_dist(ts.tt_jd(t.tt, -dt))
        lon_next, lat_next, dist_next = _of_date_lon_lat_dist(ts.tt_jd(t.tt, dt))
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

    # =========================================================================
    # Step 8: Apply sidereal correction if requested
    # =========================================================================
    if (iflag & FLG_SIDEREAL) and not (iflag & FLG_EQUATORIAL):
        # Ecliptic sidereal output only — SID|EQ keeps the plain mean-equator
        # RA (no ayanamsha; nutation already skipped above). Flag-aware
        # ayanamsa (true for of-date, mean for NONUT/J2000) and the
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
        path: Directory to save file (default: the configured SPK cache
            directory when one is set, otherwise ``<data dir>/spk``)
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
    if not overwrite and not _spk_covers_dates(spk_path, start, end):
        get_logger().warning(
            "Cached SPK file %s does not cover %s..%s, re-downloading",
            spk_path,
            start,
            end,
        )
        spk_path = download_spk(
            body=body,
            start=start,
            end=end,
            path=path,
            center=center,
            overwrite=True,
            timeout=timeout,
        )

    # Deduce NAIF ID if not provided.  The downloaded kernel's own target
    # IDs are authoritative, so resolve from the file first and use the
    # name only as a confirmation/fallback -- mirroring the already-correct
    # spk_auto._register_spk_after_download path.  A name-based deduction
    # can otherwise return a fabricated-but-non-None ID (e.g. 2000000+year
    # for a provisional designation like "2003 UB313" or a cometary
    # "1P/Halley"), which would bypass the single-target scan and register
    # a body whose NAIF ID is absent from the kernel, silently falling back
    # to the Keplerian approximation.
    if naif_id is None:
        # A permanent catalog number lets us match the exact convention the
        # kernel uses (20000000+N vs 2000000+N).
        asteroid_number = _extract_asteroid_number(body)
        if asteroid_number is not None:
            naif_id = _find_naif_id_for_asteroid(spk_path, asteroid_number)

        if naif_id is None:
            # Trust the file: a lone small-body target (NAIF > 1e6) is the
            # object we just downloaded.  This handles name-syntax bodies
            # like "Ceres;" and provisional/cometary designations where no
            # catalog number can be extracted from the name.
            targets = _get_spk_targets(spk_path)
            small_bodies = {t for t in targets if t > 1_000_000}
            if len(small_bodies) == 1:
                naif_id = small_bodies.pop()
            else:
                # Broader single-target scan excluding common centers
                # (Sun=10, SSB=0) for kernels using non-small-body numbering.
                unique_targets = {t for t in targets if t not in (0, 10)}
                if len(unique_targets) == 1:
                    naif_id = unique_targets.pop()

        if naif_id is None:
            # Last resort: legacy name-based deduction.
            naif_id = _deduce_naif_id(body)

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

    Two failure classes are kept apart, because every caller treats them
    differently:

    - ``EphemerisRangeError``: the epoch lies outside the kernel's usable
      coverage (see get_spk_coverage). Callers may continue with the
      documented lower-precision chain.
    - ``SPKEvaluationError`` (type 21): the epoch lies inside the usable
      coverage and the kernel still could not serve it. No caller catches
      it; the error reaches the user instead of a silently degraded value.

    Args:
        t: Skyfield Time object
        ipl: libephemeris body ID
        iflag: Calculation flags (FLG_SPEED, FLG_HELCTR, etc.)

    Returns:
        Position tuple (lon, lat, dist, speed_lon, speed_lat, speed_dist)
        or None if body not registered (type 2/3: kernel or target missing).

    Raises:
        EphemerisRangeError: If the epoch is outside the usable coverage.
        SPKEvaluationError: If a type 21 kernel cannot be opened, its center
            cannot be resolved, or it fails inside its usable coverage.
    """
    from . import state
    from .constants import (
        FLG_EQUATORIAL,
        FLG_HELCTR,
        FLG_SPEED,
        FLG_SIDEREAL,
        FLG_NONUT,
        FLG_J2000,
        FLG_TOPOCTR,
    )

    # Check if body is registered
    if ipl not in state._SPK_BODY_MAP:
        return None

    spk_file, naif_id = state._SPK_BODY_MAP[ipl]

    # Detect SPK type
    spk_type = _detect_spk_type(spk_file)

    # Handle type 21 (JPL Horizons asteroids/comets)
    if spk_type == 21:
        kernel = _load_type21_kernel(spk_file)
        if kernel is None:
            from .exceptions import SPKEvaluationError

            raise SPKEvaluationError(
                f"SPK type 21 kernel {spk_file} registered for body {ipl} "
                "could not be opened",
                body_id=ipl,
                naif_id=naif_id,
                spk_file=spk_file,
            )
        center = _type21_declared_center(kernel, naif_id, spk_file, ipl)

        # Gate on the usable coverage, in TDB (the scale the kernel is
        # evaluated in). Every epoch inside it can be served: the band the
        # light-time retardation and the speed stencil cannot reach is
        # already excluded, so nothing inside the reported span is refused.
        jd_tdb = float(t.tdb)
        span = _type21_usable_span(kernel, naif_id, center)
        if span is not None:
            usable_start, usable_end = span
            if jd_tdb < usable_start or jd_tdb >= usable_end:
                from .exceptions import EphemerisRangeError

                raise EphemerisRangeError(
                    f"JD {jd_tdb:.6f} (TDB) is outside the usable coverage "
                    f"[{usable_start:.6f}, {usable_end:.6f}) of SPK kernel "
                    f"{spk_file} for body {ipl}",
                    requested_jd=jd_tdb,
                    start_jd=usable_start,
                    end_jd=usable_end,
                    body_id=ipl,
                    ephemeris_file=spk_file,
                )

        target = _type21_target_for(kernel, naif_id, spk_file, ipl, center=center)
        return _calc_type21_position(kernel, naif_id, t, iflag, target=target)

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
    elif iflag & FLG_TOPOCTR:
        # Topocentric observer, like the type-21 and planetary-moon paths;
        # without this branch the diurnal parallax was silently dropped and
        # a topocentric request returned the geocentric place.
        topo = state.get_topo()
        if topo is None:
            from .exceptions import ConfigurationError

            raise ConfigurationError(
                "FLG_TOPOCTR requires a geographic position: "
                "call set_topo(lon, lat, alt) first",
                missing_config="observer_location",
                suggestion="Call set_topo(lon, lat, alt) first",
            )
        observer = planets["earth"] + topo
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
    # FLG_TRUEPOS requests the geometric place: it suppresses aberration too,
    # exactly like the Keplerian path (planets.py) and planetary_moons.py.
    apply_aberr = (
        not (iflag & FLG_TRUEPOS) and not (iflag & FLG_NOABERR) and not is_heliocentric
    )
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
    # SIDEREAL|EQUATORIAL also lands on the mean frame: the reported RA sits
    # on the mean equator of date with NO ayanamsha subtraction (the shared
    # convention of every other body class — see fast_calc's SID+EQ branch),
    # so the nutation term must be stripped here just like NONUT/J2000.
    _sid_eq = (iflag & FLG_SIDEREAL) and (iflag & FLG_EQUATORIAL)
    if (iflag & FLG_NONUT) or (iflag & FLG_J2000) or _sid_eq:
        from .cache import get_cached_nutation

        dpsi_rad, _ = get_cached_nutation(t.tt)
        lon = (lon - math.degrees(dpsi_rad)) % 360.0
        if iflag & FLG_SPEED:
            from .planets import _nutation_rate_deg_per_day

            speed_lon -= _nutation_rate_deg_per_day(t.tt)

    # Apply sidereal correction if requested — ecliptic output only (the
    # ayanamsha is a longitude offset; SID|EQ output keeps the plain mean-
    # equator RA). Routed through the canonical flag-aware helper so the
    # ayanamsa variant matches the longitude frame set above: TRUE ayanamsa
    # (mean + Δψ) cancels Δψ on a true-of-date longitude, while NONUT/J2000
    # (Δψ already stripped) get the MEAN ayanamsa. The helper also removes
    # the ~50"/yr ayanamsha drift from speed_lon under FLG_SPEED.
    if (iflag & FLG_SIDEREAL) and not (iflag & FLG_EQUATORIAL):
        from .planets import _apply_sidereal_correction

        lon, speed_lon = _apply_sidereal_correction(lon, speed_lon, t.ut1, iflag)

    return (lon, lat, dist, speed_lon, speed_lat, speed_dist)
