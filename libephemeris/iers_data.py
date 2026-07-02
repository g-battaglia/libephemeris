# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
IERS data download and parsing for observed Delta T values.

This module provides automatic download and parsing of Earth Orientation
Parameters (EOP) data from IERS (International Earth Rotation and Reference
Systems Service) for computing Delta T values for recent dates.

Delta T = TT - UT1, where:
    - TT = Terrestrial Time (atomic time + 32.184s)
    - UT1 = Universal Time (tied to Earth's rotation)

For recent dates (since ~1973), IERS provides observed UT1-UTC values which
can be used to compute high-precision Delta T values.

Data Sources:
    - IERS finals2000A.data: Earth Orientation Parameters (EOP)
      URL: https://datacenter.iers.org/data/latestVersion/finals2000A.data
      Alternative: https://maia.usno.navy.mil/ser7/finals2000A.data

    - IERS Leap Seconds: Required to compute Delta T from UT1-UTC
      URL: https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat

References:
    - IERS EOP Data: https://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html
    - IERS Bulletin A: https://www.iers.org/IERS/EN/Publications/Bulletins/Bulletins.html
"""

from __future__ import annotations

import math
import os
import threading
import time
import urllib.request
import urllib.error
from typing import Optional, Union
from dataclasses import dataclass


from .logging_config import get_logger

# =============================================================================
# MODULE CONFIGURATION
# =============================================================================

# Primary IERS data sources
IERS_FINALS_URL = "https://datacenter.iers.org/data/latestVersion/finals2000A.data"
IERS_FINALS_URL_BACKUP = "https://maia.usno.navy.mil/ser7/finals2000A.data"

# Leap seconds data
IERS_LEAP_SECONDS_URL = "https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat"
IERS_LEAP_SECONDS_URL_BACKUP = "https://www.ietf.org/timezones/data/leap-seconds.list"

# IERS observed Delta T data (monthly values from 1973-present)
IERS_DELTA_T_URL = "https://maia.usno.navy.mil/ser7/deltat.data"
IERS_DELTA_T_URL_BACKUP = "https://datacenter.iers.org/data/latestVersion/EOP_14_C04_IAU1980_one_file_1962-now.txt"

# Default cache settings
DEFAULT_CACHE_DIR = "iers_cache"
DEFAULT_MAX_AGE_DAYS = 30  # Re-download IERS data if older than this

# Lock for thread-safe operations
_IERS_LOCK = threading.RLock()

# =============================================================================
# DATA STRUCTURES
# =============================================================================


@dataclass
class IERSDataPoint:
    """A single IERS Earth Orientation Parameter data point.

    Attributes:
        mjd: Modified Julian Date
        year: Calendar year
        month: Calendar month (1-12)
        day: Calendar day (1-31)
        ut1_utc: UT1-UTC in seconds (observed value)
        ut1_utc_err: Error in UT1-UTC (seconds), if available
        is_prediction: True if this is a prediction, False if observed
    """

    mjd: float
    year: int
    month: int
    day: int
    ut1_utc: float
    ut1_utc_err: Optional[float] = None
    is_prediction: bool = False


@dataclass
class LeapSecondEntry:
    """A leap second entry.

    Attributes:
        mjd: Modified Julian Date when the leap second takes effect
        tai_utc: TAI-UTC in seconds after this date
        year: Calendar year
        month: Calendar month
        day: Calendar day
    """

    mjd: float
    tai_utc: float
    year: int
    month: int
    day: int


@dataclass
class DeltaTDataPoint:
    """An observed Delta T data point from IERS.

    Attributes:
        mjd: Modified Julian Date
        year: Calendar year
        month: Calendar month
        day: Calendar day
        delta_t: Delta T value in seconds (TT - UT1)
    """

    mjd: float
    year: int
    month: int
    day: int
    delta_t: float


# =============================================================================
# MODULE STATE
# =============================================================================

# Cached IERS data: {mjd: IERSDataPoint}
_IERS_DATA: dict[float, IERSDataPoint] = {}

# Cached leap seconds: list of LeapSecondEntry, sorted by MJD
_LEAP_SECONDS: list[LeapSecondEntry] = []

# Cached observed Delta T values: list of DeltaTDataPoint, sorted by MJD
_DELTA_T_DATA: list[DeltaTDataPoint] = []

# Last update timestamps
_IERS_DATA_TIMESTAMP: Optional[float] = None
_LEAP_SECONDS_TIMESTAMP: Optional[float] = None
_DELTA_T_DATA_TIMESTAMP: Optional[float] = None

# Cache directory override
_IERS_CACHE_DIR: Optional[str] = None

# Auto-download enabled flag
_IERS_AUTO_DOWNLOAD: Optional[bool] = None

# Environment variable for auto-download
_IERS_AUTO_DOWNLOAD_ENV_VAR = "LIBEPHEMERIS_IERS_AUTO_DOWNLOAD"


# =============================================================================
# CONFIGURATION FUNCTIONS
# =============================================================================


def set_iers_cache_dir(path: Optional[str]) -> None:
    """
    Set the directory for caching IERS data files.

    When set, IERS data files will be stored in this directory instead of the
    default cache location (within the library path).

    Args:
        path: Directory path for IERS cache, or None to use the default.
              The directory will be created if it doesn't exist.

    Example:
        >>> from libephemeris.iers_data import set_iers_cache_dir
        >>> set_iers_cache_dir("/data/ephemeris/iers_cache")
    """
    global _IERS_CACHE_DIR
    _IERS_CACHE_DIR = path


def get_iers_cache_dir() -> Optional[str]:
    """
    Get the current IERS cache directory override.

    Returns:
        The custom cache directory path if set, or None if using the default.
    """
    return _IERS_CACHE_DIR


def set_iers_auto_download(enabled: Optional[bool]) -> None:
    """
    Enable or disable automatic IERS data download.

    When enabled, the library will automatically download IERS data files
    when they are not available or are outdated.

    Args:
        enabled: True to enable automatic download, False to disable,
                 or None to use the environment variable.

    Environment Variable:
        LIBEPHEMERIS_IERS_AUTO_DOWNLOAD: Set to "1", "true", or "yes" to enable.

    Example:
        >>> from libephemeris.iers_data import set_iers_auto_download
        >>> set_iers_auto_download(True)
    """
    global _IERS_AUTO_DOWNLOAD
    _IERS_AUTO_DOWNLOAD = enabled


def get_iers_auto_download() -> bool:
    """
    Get the current IERS auto-download setting.

    Returns:
        True if automatic IERS download is enabled, False otherwise.
    """
    if _IERS_AUTO_DOWNLOAD is not None:
        return _IERS_AUTO_DOWNLOAD

    env_value = os.environ.get(_IERS_AUTO_DOWNLOAD_ENV_VAR, "").lower().strip()
    if env_value in ("1", "true", "yes", "on", "enabled"):
        return True
    if env_value in ("0", "false", "no", "off", "disabled"):
        return False

    # TOML config fallback
    from ._config_toml import get_bool as _toml_bool

    toml_value = _toml_bool("iers_auto_download")
    if toml_value is not None:
        return toml_value

    return False


# =============================================================================
# CACHE PATH HELPERS
# =============================================================================


def _get_cache_dir() -> str:
    """Get the IERS cache directory, creating it if necessary."""
    from .state import _get_data_dir

    if _IERS_CACHE_DIR is not None:
        cache_path = os.path.abspath(_IERS_CACHE_DIR)
    else:
        cache_path = os.path.join(_get_data_dir(), DEFAULT_CACHE_DIR)

    if not os.path.exists(cache_path):
        os.makedirs(cache_path, exist_ok=True)

    return cache_path


def _get_finals_cache_path() -> str:
    """Get the path for the cached finals2000A.data file."""
    return os.path.join(_get_cache_dir(), "finals2000A.data")


def _get_leap_seconds_cache_path() -> str:
    """Get the path for the cached leap seconds file."""
    return os.path.join(_get_cache_dir(), "leap_seconds.dat")


def _get_delta_t_cache_path() -> str:
    """Get the path for the cached IERS Delta T file."""
    return os.path.join(_get_cache_dir(), "deltat.data")


# =============================================================================
# DOWNLOAD FUNCTIONS
# =============================================================================


def _download_file(url: str, output_path: str, timeout: int = 30) -> bool:
    """
    Download a file from a URL to a local path.

    Thin returns-bool wrapper around download.download_file(), which
    provides the shared atomic pipeline (mkstemp -> chunked streaming ->
    publish_temp_file with world-readable permissions -> cleanup on
    error). Never raises.

    Args:
        url: URL to download from
        output_path: Local path to save the file
        timeout: Download timeout in seconds

    Returns:
        True if download succeeded, False otherwise
    """
    from pathlib import Path

    from .download import download_file

    try:
        return download_file(
            url,
            Path(output_path),
            description=os.path.basename(output_path),
            show_progress=False,
            timeout=timeout,
        )
    except (
        urllib.error.URLError,
        OSError,
        ValueError,
        KeyError,
        RuntimeError,
        TimeoutError,
    ):
        return False


def download_iers_finals(force: bool = False) -> str:
    """
    Download the IERS finals2000A.data file.

    This file contains observed and predicted Earth Orientation Parameters
    including UT1-UTC values needed for Delta T calculations.

    Args:
        force: If True, download even if a cached file exists

    Returns:
        Path to the downloaded/cached file

    Raises:
        ConnectionError: If download fails from all sources
    """
    logger = get_logger()
    cache_path = _get_finals_cache_path()

    # Check if we need to download
    if not force and os.path.exists(cache_path):
        # Check file age
        file_age_days = (time.time() - os.path.getmtime(cache_path)) / 86400
        if file_age_days < DEFAULT_MAX_AGE_DAYS:
            return cache_path

    logger.info("Downloading IERS finals2000A.data...")

    # Try primary source first
    with _IERS_LOCK:
        if _download_file(IERS_FINALS_URL, cache_path):
            logger.info("IERS finals2000A.data download complete")
            return cache_path

        # Try backup source
        logger.info("Primary source failed, trying backup...")
        if _download_file(IERS_FINALS_URL_BACKUP, cache_path):
            logger.info("IERS finals2000A.data download complete")
            return cache_path

    raise ConnectionError(
        "Failed to download IERS finals2000A.data from all sources. "
        f"Tried: {IERS_FINALS_URL}, {IERS_FINALS_URL_BACKUP}"
    )


def download_leap_seconds(force: bool = False) -> str:
    """
    Download the IERS leap seconds file.

    This file contains the history of leap seconds which is needed to
    convert between TAI and UTC.

    Args:
        force: If True, download even if a cached file exists

    Returns:
        Path to the downloaded/cached file

    Raises:
        ConnectionError: If download fails from all sources
    """
    logger = get_logger()
    cache_path = _get_leap_seconds_cache_path()

    # Check if we need to download
    if not force and os.path.exists(cache_path):
        # Check file age
        file_age_days = (time.time() - os.path.getmtime(cache_path)) / 86400
        if file_age_days < DEFAULT_MAX_AGE_DAYS:
            return cache_path

    logger.info("Downloading IERS leap seconds data...")

    # Try primary source first
    with _IERS_LOCK:
        if _download_file(IERS_LEAP_SECONDS_URL, cache_path):
            logger.info("Leap seconds data download complete")
            return cache_path

        # Try backup source
        logger.info("Primary source failed, trying backup...")
        if _download_file(IERS_LEAP_SECONDS_URL_BACKUP, cache_path):
            logger.info("Leap seconds data download complete")
            return cache_path

    raise ConnectionError(
        "Failed to download leap seconds data from all sources. "
        f"Tried: {IERS_LEAP_SECONDS_URL}, {IERS_LEAP_SECONDS_URL_BACKUP}"
    )


def download_delta_t_data(force: bool = False) -> str:
    """
    Download the IERS observed Delta T data file.

    This file contains observed Delta T (TT - UT1) values from IERS,
    providing monthly values from 1973 to present. These are the definitive
    observed values for recent dates, as opposed to computed/predicted values.

    Args:
        force: If True, download even if a cached file exists

    Returns:
        Path to the downloaded/cached file

    Raises:
        ConnectionError: If download fails from all sources

    Note:
        The data file contains observed Delta T values in seconds for the
        first of each month. Format: YEAR MONTH DAY DELTA_T

        This data is essential for high-precision calculations for recent
        dates (1973-present).
    """
    logger = get_logger()
    cache_path = _get_delta_t_cache_path()

    # Check if we need to download
    if not force and os.path.exists(cache_path):
        # Check file age
        file_age_days = (time.time() - os.path.getmtime(cache_path)) / 86400
        if file_age_days < DEFAULT_MAX_AGE_DAYS:
            return cache_path

    logger.info("Downloading IERS Delta T data...")

    # Try primary source first
    with _IERS_LOCK:
        if _download_file(IERS_DELTA_T_URL, cache_path):
            logger.info("IERS Delta T data download complete")
            return cache_path

        # Try backup source (note: different format, will need different parsing)
        logger.info("Primary source failed, trying backup...")
        if _download_file(IERS_DELTA_T_URL_BACKUP, cache_path):
            logger.info("IERS Delta T data download complete")
            return cache_path

    raise ConnectionError(
        "Failed to download IERS Delta T data from all sources. "
        f"Tried: {IERS_DELTA_T_URL}, {IERS_DELTA_T_URL_BACKUP}"
    )


# =============================================================================
# PARSING FUNCTIONS
# =============================================================================


def _parse_finals_data(filepath: str) -> dict[float, IERSDataPoint]:
    """
    Parse the IERS finals2000A.data file.

    The file uses fixed-width format with the following columns:
    - Columns 1-2: Year (last 2 digits)
    - Columns 3-4: Month
    - Columns 5-6: Day
    - Columns 8-15: MJD
    - Column 17: IERS (I) or Prediction (P) flag for polar motion
    - Columns 59-68: UT1-UTC (seconds)
    - Column 58: IERS (I) or Prediction (P) flag for UT1-UTC
    - Columns 69-78: Error in UT1-UTC (seconds)

    Args:
        filepath: Path to the finals2000A.data file

    Returns:
        Dictionary mapping MJD to IERSDataPoint
    """
    data: dict[float, IERSDataPoint] = {}

    with open(filepath, "r", encoding="latin-1") as f:
        for line in f:
            if len(line) < 68:
                continue

            try:
                # Parse date fields
                year_2digit = int(line[0:2].strip())
                month = int(line[2:4].strip())
                day = int(line[4:6].strip())

                # Convert 2-digit year to 4-digit
                # Assuming data from 1973-2099
                if year_2digit >= 73:
                    year = 1900 + year_2digit
                else:
                    year = 2000 + year_2digit

                # Parse MJD
                mjd_str = line[7:15].strip()
                if not mjd_str:
                    continue
                mjd = float(mjd_str)

                # Parse UT1-UTC value (columns 59-68, 1-indexed, so 58:68 in Python)
                ut1_utc_str = line[58:68].strip()
                if not ut1_utc_str:
                    continue
                ut1_utc = float(ut1_utc_str)

                # Check if this is observed (I) or predicted (P)
                # The flag is at column 58 (1-indexed), so 57 in Python
                is_prediction = False
                if len(line) > 57:  # pragma: no branch - earlier len<68 guard makes this always true
                    flag = line[57]
                    is_prediction = flag.upper() == "P"

                # Parse error if available (columns 69-78)
                ut1_utc_err = None
                if len(line) >= 78:
                    err_str = line[68:78].strip()
                    if err_str:
                        try:
                            ut1_utc_err = float(err_str)
                        except ValueError:
                            pass

                data[mjd] = IERSDataPoint(
                    mjd=mjd,
                    year=year,
                    month=month,
                    day=day,
                    ut1_utc=ut1_utc,
                    ut1_utc_err=ut1_utc_err,
                    is_prediction=is_prediction,
                )

            except (ValueError, IndexError):
                # Skip malformed lines
                continue

    return data


def _parse_leap_seconds(filepath: str) -> list[LeapSecondEntry]:
    """
    Parse leap seconds data file.

    Supports two formats:
    1. IERS Bulletin C format (Leap_Second.dat)
    2. IETF format (leap-seconds.list)

    Args:
        filepath: Path to the leap seconds file

    Returns:
        List of LeapSecondEntry sorted by MJD
    """
    entries: list[LeapSecondEntry] = []

    with open(filepath, "r", encoding="latin-1") as f:
        content = f.read()

    lines = content.strip().split("\n")

    # Detect format. The IETF leap-seconds.list is identified by its NTP-list
    # markers (#$ last-update, #@ expiry) or the literal "NTP"/filename. Do NOT
    # route on a bare leading '#': the primary IERS Leap_Second.dat is ALSO
    # '#'-commented but has numeric MJD data rows, and routing it to the IETF
    # parser silently discarded every downloaded IERS update (falling back to the
    # hardcoded list forever).
    is_ietf = (
        "NTP" in content
        or "#$" in content
        or "#@" in content
        or "leap-seconds.list" in content
    )
    if is_ietf:
        entries = _parse_leap_seconds_ietf(lines)
    else:
        entries = _parse_leap_seconds_iers(lines)

    # If the chosen parser found nothing, try the other before giving up — a
    # mis-detected file should still be parsed rather than silently dropped.
    if not entries:
        entries = (
            _parse_leap_seconds_iers(lines)
            if is_ietf
            else _parse_leap_seconds_ietf(lines)
        )

    # If parsing still failed, use hardcoded values as fallback
    if not entries:
        entries = _get_hardcoded_leap_seconds()

    return sorted(entries, key=lambda x: x.mjd)


def _parse_leap_seconds_iers(lines: list[str]) -> list[LeapSecondEntry]:
    """Parse the IERS leap-seconds file (Leap_Second.dat).

    Args:
        lines: Lines read from an IERS leap-second data file.

    Returns:
        Parsed leap-second entries (one per recognised data row).

    The primary IERS source uses all-numeric columns::

        #  MJD        Date        TAI-UTC (s)
        #           day month year
        41317.0    1  1 1972      10
        ...
        57754.0    1  1 2017      37

    i.e. ``MJD day month year TAI-UTC``. A textual ``year month-name day ...
    TAI-UTC[= Ns]`` variant is also tolerated so callers feeding the Bulletin-C
    prose form still parse.
    """
    months = {
        "jan": 1, "feb": 2, "mar": 3, "apr": 4, "may": 5, "jun": 6,
        "jul": 7, "aug": 8, "sep": 9, "oct": 10, "nov": 11, "dec": 12,
    }
    entries: list[LeapSecondEntry] = []

    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        try:
            parts = line.split()
            if len(parts) < 5:
                continue

            if parts[1].isdigit() and parts[2].isdigit() and parts[3].isdigit():
                # Numeric Leap_Second.dat row: MJD day month year TAI-UTC
                mjd = float(parts[0])
                day = int(parts[1])
                month = int(parts[2])
                year = int(parts[3])
                tai_utc = float(parts[4])
            else:
                # Textual variant: year month-name day ... TAI-UTC[= Ns]
                year = int(parts[0])
                month = months.get(parts[1].lower()[:3], 0)
                if month == 0:
                    continue
                day = int(parts[2])
                tai_utc = float(parts[4].rstrip("s"))
                mjd = _calendar_to_mjd(year, month, day)

            entries.append(
                LeapSecondEntry(
                    mjd=mjd,
                    tai_utc=tai_utc,
                    year=year,
                    month=month,
                    day=day,
                )
            )
        except (ValueError, IndexError):
            continue

    return entries


def _parse_leap_seconds_ietf(lines: list[str]) -> list[LeapSecondEntry]:
    """Parse IETF format leap seconds file (NTP timestamps)."""
    entries: list[LeapSecondEntry] = []

    # NTP epoch: January 1, 1900, 00:00:00 UTC
    # MJD of NTP epoch: 15020.0
    NTP_EPOCH_MJD = 15020.0

    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        try:
            parts = line.split()
            if len(parts) >= 2:
                ntp_seconds = int(parts[0])
                tai_utc = int(parts[1])

                # Convert NTP seconds to MJD
                mjd = NTP_EPOCH_MJD + ntp_seconds / 86400.0

                # Convert to calendar date
                year, month, day = _mjd_to_calendar(mjd)

                entries.append(
                    LeapSecondEntry(
                        mjd=mjd,
                        tai_utc=float(tai_utc),
                        year=year,
                        month=month,
                        day=day,
                    )
                )
        except (ValueError, IndexError):
            continue

    return entries


def _get_hardcoded_leap_seconds() -> list[LeapSecondEntry]:
    """
    Return hardcoded leap seconds as fallback.

    This is used when online data is unavailable. Last updated: January 2025.
    Note: Future leap seconds may be added by IERS; this list may become outdated.
    """
    # List of (year, month, day, TAI-UTC after)
    leap_second_data = [
        (1972, 1, 1, 10),
        (1972, 7, 1, 11),
        (1973, 1, 1, 12),
        (1974, 1, 1, 13),
        (1975, 1, 1, 14),
        (1976, 1, 1, 15),
        (1977, 1, 1, 16),
        (1978, 1, 1, 17),
        (1979, 1, 1, 18),
        (1980, 1, 1, 19),
        (1981, 7, 1, 20),
        (1982, 7, 1, 21),
        (1983, 7, 1, 22),
        (1985, 7, 1, 23),
        (1988, 1, 1, 24),
        (1990, 1, 1, 25),
        (1991, 1, 1, 26),
        (1992, 7, 1, 27),
        (1993, 7, 1, 28),
        (1994, 7, 1, 29),
        (1996, 1, 1, 30),
        (1997, 7, 1, 31),
        (1999, 1, 1, 32),
        (2006, 1, 1, 33),
        (2009, 1, 1, 34),
        (2012, 7, 1, 35),
        (2015, 7, 1, 36),
        (2017, 1, 1, 37),
    ]

    entries = []
    for year, month, day, tai_utc in leap_second_data:
        mjd = _calendar_to_mjd(year, month, day)
        entries.append(
            LeapSecondEntry(
                mjd=mjd,
                tai_utc=float(tai_utc),
                year=year,
                month=month,
                day=day,
            )
        )

    return entries


def _parse_delta_t_data(filepath: str) -> list[DeltaTDataPoint]:
    """
    Parse IERS Delta T data file.

    The file format is simple whitespace-separated columns:
        YEAR MONTH DAY DELTA_T

    Where DELTA_T is in seconds.

    Args:
        filepath: Path to the deltat.data file

    Returns:
        List of DeltaTDataPoint sorted by MJD
    """
    entries: list[DeltaTDataPoint] = []

    with open(filepath, "r", encoding="latin-1") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            try:
                parts = line.split()
                if len(parts) >= 4:
                    year = int(parts[0])
                    month = int(parts[1])
                    day = int(parts[2])
                    delta_t = float(parts[3])

                    # Validate ranges
                    if not (1900 <= year <= 2100):
                        continue
                    if not (1 <= month <= 12):
                        continue
                    if not (1 <= day <= 31):
                        continue

                    # Calculate MJD
                    mjd = _calendar_to_mjd(year, month, day)

                    entries.append(
                        DeltaTDataPoint(
                            mjd=mjd,
                            year=year,
                            month=month,
                            day=day,
                            delta_t=delta_t,
                        )
                    )
            except (ValueError, IndexError):
                # Skip malformed lines
                continue

    return sorted(entries, key=lambda x: x.mjd)


# =============================================================================
# DATE CONVERSION HELPERS
# =============================================================================


def _calendar_to_mjd(year: int, month: int, day: int) -> float:
    """Convert calendar date to Modified Julian Date."""
    # JD calculation (Meeus algorithm)
    if month <= 2:
        year -= 1
        month += 12

    a = int(year / 100)
    b = 2 - a + int(a / 4)

    jd = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + b - 1524.5

    # MJD = JD - 2400000.5
    return jd - 2400000.5


def _mjd_to_calendar(mjd: float) -> tuple[int, int, int]:
    """Convert Modified Julian Date to calendar date (year, month, day)."""
    jd = mjd + 2400000.5

    jd = jd + 0.5
    z = int(jd)
    f = jd - z

    if z < 2299161:
        a = z
    else:
        alpha = int((z - 1867216.25) / 36524.25)
        a = z + 1 + alpha - int(alpha / 4)

    b = a + 1524
    c = int((b - 122.1) / 365.25)
    d = int(365.25 * c)
    e = int((b - d) / 30.6001)

    day = int(b - d - int(30.6001 * e) + f)

    if e < 14:
        month = e - 1
    else:
        month = e - 13

    if month > 2:
        year = c - 4716
    else:
        year = c - 4715

    return year, month, day


def _jd_to_mjd(jd: float) -> float:
    """Convert Julian Date to Modified Julian Date."""
    return jd - 2400000.5


def _mjd_to_jd(mjd: float) -> float:
    """Convert Modified Julian Date to Julian Date."""
    return mjd + 2400000.5


# =============================================================================
# DATA LOADING AND ACCESS
# =============================================================================


def load_iers_data(force_download: bool = False) -> bool:
    """
    Load IERS data from cache, downloading if necessary.

    Args:
        force_download: If True, force re-download of data

    Returns:
        True if data was loaded successfully, False otherwise
    """
    global _IERS_DATA, _LEAP_SECONDS, _DELTA_T_DATA
    global _IERS_DATA_TIMESTAMP, _LEAP_SECONDS_TIMESTAMP, _DELTA_T_DATA_TIMESTAMP

    with _IERS_LOCK:
        finals_path = _get_finals_cache_path()
        leap_path = _get_leap_seconds_cache_path()
        delta_t_path = _get_delta_t_cache_path()

        # Check if we need to download/reload
        need_finals = force_download or not os.path.exists(finals_path)
        need_leap = force_download or not os.path.exists(leap_path)
        need_delta_t = force_download or not os.path.exists(delta_t_path)

        # Auto-download if enabled
        if get_iers_auto_download():
            try:
                if need_finals or force_download:
                    download_iers_finals(force=force_download)
                if need_leap or force_download:
                    download_leap_seconds(force=force_download)
                if need_delta_t or force_download:
                    download_delta_t_data(force=force_download)
            except ConnectionError:
                # If download fails, try to use existing cached data
                pass

        # Load finals data
        if os.path.exists(finals_path):
            try:
                _IERS_DATA = _parse_finals_data(finals_path)
                _IERS_DATA_TIMESTAMP = time.time()
            except (OSError, ValueError, ImportError):
                _IERS_DATA = {}

        # Load leap seconds data
        if os.path.exists(leap_path):
            try:
                _LEAP_SECONDS = _parse_leap_seconds(leap_path)
                _LEAP_SECONDS_TIMESTAMP = time.time()
            except (OSError, ValueError, ImportError):
                _LEAP_SECONDS = _get_hardcoded_leap_seconds()
        else:
            # Use hardcoded fallback
            _LEAP_SECONDS = _get_hardcoded_leap_seconds()
            _LEAP_SECONDS_TIMESTAMP = time.time()

        # Load observed Delta T data
        if os.path.exists(delta_t_path):
            try:
                _DELTA_T_DATA = _parse_delta_t_data(delta_t_path)
                _DELTA_T_DATA_TIMESTAMP = time.time()
            except (OSError, ValueError, ImportError):
                _DELTA_T_DATA = []
        else:
            _DELTA_T_DATA = []

    return bool(_IERS_DATA) or bool(_LEAP_SECONDS) or bool(_DELTA_T_DATA)


def _ensure_data_loaded() -> None:
    """Ensure IERS data is loaded, loading if necessary."""
    if not _IERS_DATA and not _LEAP_SECONDS and not _DELTA_T_DATA:
        load_iers_data()


def get_tai_utc(mjd: float) -> float:
    """
    Get TAI-UTC (leap seconds) for a given MJD.

    Args:
        mjd: Modified Julian Date

    Returns:
        TAI-UTC in seconds
    """
    _ensure_data_loaded()

    if not _LEAP_SECONDS:
        # Fallback: use approximate value for modern dates
        if mjd >= 57754.0:  # 2017-01-01
            return 37.0
        return 32.0

    # Find the applicable leap second entry
    tai_utc = 0.0
    for entry in _LEAP_SECONDS:
        if mjd >= entry.mjd:
            tai_utc = entry.tai_utc
        else:
            break

    return tai_utc


def get_ut1_utc(mjd: float) -> Optional[float]:
    """
    Get UT1-UTC for a given MJD from IERS data.

    Args:
        mjd: Modified Julian Date

    Returns:
        UT1-UTC in seconds, or None if data not available for this date
    """
    _ensure_data_loaded()

    # Snapshot: clear_iers_cache() rebinds the global under the lock, so a
    # local binding makes the reads race-free.
    data = _IERS_DATA
    if not data:
        return None

    # Interpolate between the bracketing daily entries. (The previous
    # nearest-midnight lookup short-circuited before interpolation could
    # ever run, turning DUT1 into a step function with ~ms-level error.)
    mjd_floor = math.floor(mjd)
    mjd_ceil = mjd_floor + 1
    frac = mjd - mjd_floor

    if frac == 0.0 and mjd_floor in data:
        return data[mjd_floor].ut1_utc

    if mjd_floor in data and mjd_ceil in data:
        ut1_utc_floor = data[mjd_floor].ut1_utc
        ut1_utc_ceil = data[mjd_ceil].ut1_utc
        return ut1_utc_floor + frac * (ut1_utc_ceil - ut1_utc_floor)

    # Edge of table: closest available date
    if mjd_floor in data:
        return data[mjd_floor].ut1_utc
    if mjd_ceil in data:
        return data[mjd_ceil].ut1_utc

    return None


def get_observed_delta_t(jd: float) -> Optional[float]:
    """
    Get observed Delta T value directly from IERS deltat.data file.

    This uses the observed Delta T values published by IERS, which are
    monthly values from 1973 to present. These are the definitive observed
    values and are more accurate than computing Delta T from UT1-UTC.

    For dates between data points, linear interpolation is used.

    Args:
        jd: Julian Date (UT1)

    Returns:
        Delta T in seconds, or None if data not available for this date

    Note:
        The IERS deltat.data file contains observed Delta T values for the
        first of each month. This function interpolates between these
        monthly values.
    """
    _ensure_data_loaded()

    if not _DELTA_T_DATA:
        return None

    mjd = _jd_to_mjd(jd)

    # Check if we're within the data range
    if mjd < _DELTA_T_DATA[0].mjd or mjd > _DELTA_T_DATA[-1].mjd:
        return None

    # Binary search for the bracket
    left = 0
    right = len(_DELTA_T_DATA) - 1

    while left < right - 1:
        mid = (left + right) // 2
        if _DELTA_T_DATA[mid].mjd <= mjd:
            left = mid
        else:
            right = mid

    # Linear interpolation between the two nearest points
    p1 = _DELTA_T_DATA[left]
    p2 = _DELTA_T_DATA[right]

    if p1.mjd == p2.mjd:
        return p1.delta_t

    # Fraction of the interval
    frac = (mjd - p1.mjd) / (p2.mjd - p1.mjd)

    # Interpolate
    delta_t = p1.delta_t + frac * (p2.delta_t - p1.delta_t)

    return delta_t


def get_observed_delta_t_data_range() -> Optional[tuple[float, float]]:
    """
    Get the date range covered by observed Delta T data.

    Returns:
        Tuple of (jd_start, jd_end), or None if no data loaded
    """
    _ensure_data_loaded()

    if not _DELTA_T_DATA:
        return None

    mjd_min = _DELTA_T_DATA[0].mjd
    mjd_max = _DELTA_T_DATA[-1].mjd

    return _mjd_to_jd(mjd_min), _mjd_to_jd(mjd_max)


def is_observed_delta_t_available(jd: float) -> bool:
    """
    Check if observed Delta T data is available for a given Julian Date.

    Args:
        jd: Julian Date

    Returns:
        True if observed Delta T data covers this date
    """
    _ensure_data_loaded()

    if not _DELTA_T_DATA:
        return False

    mjd = _jd_to_mjd(jd)
    return _DELTA_T_DATA[0].mjd <= mjd <= _DELTA_T_DATA[-1].mjd


def get_delta_t_iers(jd: float) -> Optional[float]:
    """
    Get Delta T from IERS data, using the best available source.

    This function tries multiple sources in order of preference:
    1. Observed Delta T values from IERS deltat.data (most accurate for recent dates)
    2. Computed Delta T from UT1-UTC in finals2000A.data

    Args:
        jd: Julian Date (UT1)

    Returns:
        Delta T in seconds, or None if data not available for this date

    Note:
        The observed Delta T values are preferred because they are the
        definitive published values. The computed values from UT1-UTC
        are used as a fallback for dates not covered by the observed data.
    """
    _ensure_data_loaded()

    # First try observed Delta T data (preferred source)
    observed_dt = get_observed_delta_t(jd)
    if observed_dt is not None:
        return observed_dt

    # Fall back to computing Delta T from UT1-UTC
    mjd = _jd_to_mjd(jd)

    # Get UT1-UTC from IERS data
    ut1_utc = get_ut1_utc(mjd)
    if ut1_utc is None:
        return None

    # Get TAI-UTC (leap seconds)
    tai_utc = get_tai_utc(mjd)

    # Calculate Delta T
    # TT = TAI + 32.184
    # Delta T = TT - UT1 = (UTC + TAI-UTC + 32.184) - (UTC + UT1-UTC)
    #         = TAI-UTC + 32.184 - UT1-UTC
    delta_t_seconds = tai_utc + 32.184 - ut1_utc

    return delta_t_seconds


def get_iers_data_range() -> Optional[tuple[float, float]]:
    """
    Get the date range covered by loaded IERS data.

    Returns:
        Tuple of (jd_start, jd_end), or None if no data loaded
    """
    _ensure_data_loaded()

    if not _IERS_DATA:
        return None

    mjd_values = list(_IERS_DATA.keys())
    mjd_min = min(mjd_values)
    mjd_max = max(mjd_values)

    return _mjd_to_jd(mjd_min), _mjd_to_jd(mjd_max)


def is_iers_data_available(jd: float) -> bool:
    """
    Check if IERS data is available for a given Julian Date.

    Args:
        jd: Julian Date

    Returns:
        True if IERS data is available for this date
    """
    # Load on demand like every other getter — otherwise this returns
    # False until some sibling function happens to load the table first.
    _ensure_data_loaded()
    mjd = _jd_to_mjd(jd)
    # Mirror get_ut1_utc's bracketing (floor / floor+1): availability must be
    # True exactly when get_ut1_utc would return a non-None value, otherwise the
    # two disagree for the last ~half day of the table (round() vs floor()).
    mjd_floor = math.floor(mjd)
    return mjd_floor in _IERS_DATA or (mjd_floor + 1) in _IERS_DATA


def clear_iers_cache() -> None:
    """Clear cached IERS data from memory."""
    global _IERS_DATA, _LEAP_SECONDS, _DELTA_T_DATA
    global _IERS_DATA_TIMESTAMP, _LEAP_SECONDS_TIMESTAMP, _DELTA_T_DATA_TIMESTAMP

    with _IERS_LOCK:
        _IERS_DATA = {}
        _LEAP_SECONDS = []
        _DELTA_T_DATA = []
        _IERS_DATA_TIMESTAMP = None
        _LEAP_SECONDS_TIMESTAMP = None
        _DELTA_T_DATA_TIMESTAMP = None


def delete_iers_cache_files() -> int:
    """
    Delete cached IERS data files from disk.

    Returns:
        Number of files deleted
    """
    deleted = 0
    cache_dir = _get_cache_dir()

    for filename in ["finals2000A.data", "leap_seconds.dat", "deltat.data"]:
        filepath = os.path.join(cache_dir, filename)
        if os.path.exists(filepath):
            try:
                os.remove(filepath)
                deleted += 1
            except OSError:
                pass

    # Also clear in-memory cache
    clear_iers_cache()

    return deleted


def get_iers_cache_info() -> dict[str, Union[str, int, float, bool, None]]:
    """
    Get information about the IERS data cache.

    Returns:
        Dictionary with cache information:
            - cache_dir: Path to cache directory
            - finals_exists: Whether finals data file exists
            - finals_age_days: Age of finals file in days, or None
            - leap_seconds_exists: Whether leap seconds file exists
            - leap_seconds_age_days: Age of leap seconds file in days, or None
            - delta_t_exists: Whether Delta T data file exists
            - delta_t_age_days: Age of Delta T file in days, or None
            - data_points: Number of IERS data points loaded
            - leap_second_entries: Number of leap second entries loaded
            - delta_t_entries: Number of Delta T entries loaded
            - data_range_start_jd: Start of data range (JD), or None
            - data_range_end_jd: End of data range (JD), or None
            - delta_t_range_start_jd: Start of Delta T data range (JD), or None
            - delta_t_range_end_jd: End of Delta T data range (JD), or None
    """
    cache_dir = _get_cache_dir()
    finals_path = _get_finals_cache_path()
    leap_path = _get_leap_seconds_cache_path()
    delta_t_path = _get_delta_t_cache_path()

    info: dict[str, Union[str, int, float, bool, None]] = {
        "cache_dir": cache_dir,
        "finals_exists": os.path.exists(finals_path),
        "finals_age_days": None,
        "leap_seconds_exists": os.path.exists(leap_path),
        "leap_seconds_age_days": None,
        "delta_t_exists": os.path.exists(delta_t_path),
        "delta_t_age_days": None,
        "data_points": len(_IERS_DATA),
        "leap_second_entries": len(_LEAP_SECONDS),
        "delta_t_entries": len(_DELTA_T_DATA),
        "data_range_start_jd": None,
        "data_range_end_jd": None,
        "delta_t_range_start_jd": None,
        "delta_t_range_end_jd": None,
    }

    if info["finals_exists"]:
        age = (time.time() - os.path.getmtime(finals_path)) / 86400
        info["finals_age_days"] = round(age, 1)

    if info["leap_seconds_exists"]:
        age = (time.time() - os.path.getmtime(leap_path)) / 86400
        info["leap_seconds_age_days"] = round(age, 1)

    if info["delta_t_exists"]:
        age = (time.time() - os.path.getmtime(delta_t_path)) / 86400
        info["delta_t_age_days"] = round(age, 1)

    data_range = get_iers_data_range()
    if data_range:
        info["data_range_start_jd"] = data_range[0]
        info["data_range_end_jd"] = data_range[1]

    delta_t_range = get_observed_delta_t_data_range()
    if delta_t_range:
        info["delta_t_range_start_jd"] = delta_t_range[0]
        info["delta_t_range_end_jd"] = delta_t_range[1]

    return info
