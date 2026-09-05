# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
REBOUND/ASSIST integration for ephemeris-quality asteroid orbit propagation.

This module provides high-precision n-body gravitational integration for asteroid
and comet orbit propagation using the REBOUND and ASSIST libraries.

REBOUND is a multi-purpose N-body integrator with several integrator choices:
- IAS15: High-accuracy adaptive integrator (15th order Gauss-Radau)
- WHFast: Fast symplectic Wisdom-Holman integrator
- MERCURIUS: Hybrid symplectic integrator for close encounters
- TRACE: Hybrid reversible integrator for arbitrary close encounters

ASSIST extends REBOUND for ephemeris-quality integrations of Solar System test
particles, matching JPL's small body integrator precision:
- Uses JPL DE440/DE441 ephemeris for planet positions
- Includes Sun, Moon, 8 planets, and 16 massive asteroids
- Gravitational harmonics (J2, J3, J4) for Earth and Sun
- General relativistic corrections
- Non-gravitational effects (Marsden model)
- Verified to machine precision against JPL Horizons

INSTALLATION:
    pip install rebound          # Base n-body integrator (~3 MB)
    pip install assist           # Ephemeris-quality extension (~3 MB)

    # ASSIST requires ephemeris data files (~1 GB total):
    mkdir -p data
    curl https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/linux_p1550p2650.440 -o data/linux_p1550p2650.440
    curl https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n16.bsp -o data/sb441-n16.bsp

USAGE:
    >>> from libephemeris.rebound_integration import (
    ...     propagate_orbit_assist,
    ...     propagate_orbit_rebound,
    ...     ReboundIntegrator,
    ... )
    >>> from libephemeris.minor_bodies import MINOR_BODY_ELEMENTS, CERES
    >>>
    >>> # Propagate Ceres orbit using ASSIST (ephemeris-quality)
    >>> elements = MINOR_BODY_ELEMENTS[CERES]
    >>> jd_start = 2460000.5  # Starting JD
    >>> jd_end = 2460365.5    # One year later
    >>> x, y, z, vx, vy, vz = propagate_orbit_assist(elements, jd_start, jd_end)

PRECISION:
    - ASSIST: Sub-arcsecond accuracy, matching JPL Horizons
    - REBOUND IAS15: ~0.001-1 arcsecond (depends on timestep and duration)
    - REBOUND WHFast: ~1-10 arcseconds (better energy conservation for long term)
    - Keplerian (minor_bodies.py): ~10-30 arcseconds for main belt asteroids

REFERENCES:
    - REBOUND: https://rebound.readthedocs.io/
    - ASSIST: https://assist.readthedocs.io/
    - Rein & Liu 2012, A&A 537, A128 (REBOUND)
    - Rein & Spiegel 2015, MNRAS 446, 1424 (IAS15)
    - Rein & Tamayo 2015, MNRAS 452, 376 (WHFast)
    - Holman et al. 2023 (ASSIST)

Provenance:
    Numerical integration is performed by the optional upstream REBOUND/ASSIST
    packages under their own licences; ASSIST reads public JPL DE and massive-
    asteroid data. This module is a project-authored adapter that converts the
    registered osculating elements, selects an explicitly requested integrator,
    validates data files, and converts units/frames. Default timestep, tolerance,
    cache, and fallback behavior are project choices. No upstream code or state
    vector table is vendored here.
"""

from __future__ import annotations

import math
import os
from collections import OrderedDict
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from threading import RLock
from typing import TYPE_CHECKING, Callable, Optional, Tuple, List

if TYPE_CHECKING:
    from .minor_bodies import OrbitalElements


class ReboundIntegrator(Enum):
    """Available integrators in REBOUND.

    Each integrator has different trade-offs between speed, accuracy, and
    use cases.
    """

    IAS15 = "ias15"
    """15th order Gauss-Radau integrator with adaptive timestep.
    
    Best for: High-accuracy integrations, close encounters, eccentric orbits.
    Accuracy: Machine precision for most problems.
    Speed: Moderate (adaptive, slows for close encounters).
    """

    WHFAST = "whfast"
    """Fast Wisdom-Holman symplectic integrator.
    
    Best for: Long-term integrations of planetary systems.
    Accuracy: O(dt^2) per step, but bounded errors (symplectic).
    Speed: Fast (fixed timestep).
    Note: Not suitable for close encounters.
    """

    MERCURIUS = "mercurius"
    """Hybrid symplectic integrator for close encounters.
    
    Best for: Systems with occasional close encounters.
    Accuracy: Symplectic most of time, IAS15 during encounters.
    Speed: Fast until close encounters occur.
    """

    TRACE = "trace"
    """Hybrid TRAnsits and Close Encounters integrator.
    
    Best for: Arbitrary close encounters, reversible integration.
    Accuracy: High accuracy even with many close encounters.
    Speed: Moderate.
    """

    LEAPFROG = "leapfrog"
    """Standard leapfrog/Verlet integrator.
    
    Best for: Simple problems, educational use.
    Accuracy: O(dt^2) per step.
    Speed: Fast.
    """


_ASSIST_PLANETS_FILENAMES = [
    "linux_p1550p2650.440",
    "linux_m13000p17000.441",
]

_ASSIST_ASTEROIDS_FILENAMES = [
    "sb441-n16.bsp",
]

# URLs for downloading ASSIST data files
_ASSIST_PLANETS_URL = (
    "https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/linux_p1550p2650.440"
)
# DE441 long planet ephemeris (-13000..+17000) — needed only for deep-time
# (extended-tier) N-body generation; ~2.6 GB, opt-in download.
_ASSIST_PLANETS_DE441_URL = (
    "https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de441/linux_m13000p17000.441"
)
_ASSIST_ASTEROIDS_URL = (
    "https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n16.bsp"
)

# Default directory for ASSIST data
_ASSIST_DEFAULT_DIR = Path(os.path.expanduser("~/.libephemeris/assist"))


def _search_assist_file(
    filenames: list[str],
    search_dirs: list[Path],
) -> Optional[str]:
    """Search for an ASSIST data file in multiple directories."""
    for d in search_dirs:
        if not d.exists():
            continue
        for fname in filenames:
            path = d / fname
            if path.exists():
                return str(path)
    return None


@dataclass
class AssistEphemConfig:
    """Configuration for ASSIST ephemeris data files.

    ASSIST requires two ephemeris data files:
    1. Planet ephemeris (DE440/DE441 in JPL Linux binary format): ~98 MB
    2. Asteroid perturber ephemeris (optional): ~616 MB

    Files are searched in this order:
    1. Explicit paths (planets_file, asteroids_file)
    2. data_dir parameter
    3. ASSIST_DIR environment variable
    4. ~/.libephemeris/assist/
    5. ./data/

    Download with:
        from libephemeris.rebound_integration import download_assist_data
        download_assist_data()  # Downloads to ~/.libephemeris/assist/

    Attributes:
        planets_file: Path to planet ephemeris file (e.g., linux_p1550p2650.440)
        asteroids_file: Path to asteroid perturber file (e.g., sb441-n16.bsp)
        data_dir: Directory containing ephemeris files (alternative to full paths)
    """

    planets_file: Optional[str] = None
    asteroids_file: Optional[str] = None
    data_dir: Optional[str] = None

    def __post_init__(self) -> None:
        """Resolve file paths from environment or defaults."""
        # Build search directories list
        search_dirs: list[Path] = []

        if self.data_dir:
            search_dirs.append(Path(self.data_dir))

        env_dir = os.environ.get("ASSIST_DIR")
        if env_dir:
            search_dirs.append(Path(env_dir))

        search_dirs.append(_ASSIST_DEFAULT_DIR)
        search_dirs.append(Path("data"))

        # Search for planet ephemeris
        if self.planets_file is None:
            self.planets_file = _search_assist_file(
                _ASSIST_PLANETS_FILENAMES, search_dirs
            )

        # Search for asteroid perturber ephemeris
        if self.asteroids_file is None:
            self.asteroids_file = _search_assist_file(
                _ASSIST_ASTEROIDS_FILENAMES, search_dirs
            )


# Expected sizes for ASSIST data files (bytes).
# Used for basic integrity verification after download.
_ASSIST_EXPECTED_SIZES = {
    "linux_p1550p2650.440": 102272352,  # ~98 MB
    "linux_m13000p17000.441": 2788676624,  # ~2.6 GB (DE441, deep-time)
    "sb441-n16.bsp": 645727232,  # ~616 MB
}

# Minimum acceptable file size (90% of expected) to catch truncated downloads.
_ASSIST_MIN_SIZE_RATIO = 0.90


def _verify_assist_file(path: Path, name: Optional[str] = None) -> bool:
    """Verify that an ASSIST data file looks valid.

    Checks that the file exists and its size is within the expected range.
    For BSP files, also attempts structural validation via jplephem.

    Args:
        path: Path to the file to verify.
        name: Logical filename to look the expectations up under. Defaults to
            ``path.name``; pass the final name explicitly when verifying a
            temporary download, whose own name carries no expectations. The
            size threshold is the only integrity check the planet files
            (``.440`` / ``.441``) have, so it is enforced for a fresh payload
            as well: a truncated download must never replace a good file.

    Returns:
        True if the file passes verification, False otherwise.
    """
    if not path.exists():
        return False

    logical_name = name if name is not None else path.name

    actual_size = path.stat().st_size
    if actual_size == 0:
        return False

    expected = _ASSIST_EXPECTED_SIZES.get(logical_name)
    if expected is not None:
        min_size = int(expected * _ASSIST_MIN_SIZE_RATIO)
        if actual_size < min_size:
            return False

    # For BSP files, try structural validation
    if logical_name.endswith(".bsp"):
        try:
            from jplephem.spk import SPK

            with SPK.open(str(path)) as kernel:
                for segment in kernel.segments:
                    _ = segment.center, segment.target
        except (ImportError, RuntimeError, ValueError):
            return False

    return True


def _create_ssl_context():
    """Create an SSL context using certifi certificates.

    Returns:
        ssl.SSLContext configured with certifi CA bundle.
    """
    import ssl

    try:
        import certifi

        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        return ssl.create_default_context()


def _assist_download_validator(name: str) -> Callable[[str], bool]:
    """Build the pre-publication check for one ASSIST download.

    The check runs on the temporary file, whose random name carries no size
    expectation, so the final ``name`` is passed explicitly. The size
    threshold is enforced there too: for the planet files it is the only
    integrity check, so a short payload is rejected before publication and
    the file already on disk survives. A legitimate upstream republication
    that shrinks a file is installed only after ``_ASSIST_EXPECTED_SIZES``
    is updated with it.

    Args:
        name: Final filename the download will be published under.

    Returns:
        A callable taking the temp path and returning True when the payload
        passes the structural checks.
    """

    def _validate(temp_path: str) -> bool:
        return _verify_assist_file(Path(temp_path), name=name)

    return _validate


def _download_single_file(
    url: str,
    dest: Path,
    description: str,
    show_progress: bool = True,
    quiet: bool = False,
    validator: Optional[Callable[[str], bool]] = None,
) -> None:
    """Download a single file with progress, atomic write, and SSL support.

    Args:
        url: URL to download from.
        dest: Destination file path.
        description: Human-readable description for progress output.
        show_progress: Whether to show a progress bar.
        quiet: Suppress all non-error output.
        validator: Optional check applied to the *temporary* file before it
            is published. Verifying before the atomic replace is what keeps a
            truncated multi-hundred-megabyte download from destroying the
            copy already on disk.

    Raises:
        ValueError: If ``validator`` rejects the downloaded payload.
    """
    import hashlib
    import tempfile

    ssl_context = _create_ssl_context()

    dest.parent.mkdir(parents=True, exist_ok=True)
    temp_fd, temp_path = tempfile.mkstemp(dir=dest.parent, suffix=".download")
    published = False

    try:
        from .net import Request

        req = Request(url, headers={"User-Agent": "libephemeris-download/1.0"})

        from .net import open_url

        with open_url(
            req,
            purpose=f"ASSIST data download: {description}",
            timeout=300,
            context=ssl_context,
        ) as response:
            total_size = int(response.headers.get("Content-Length", 0))

            if not quiet and total_size > 0:
                size_mb = total_size / (1024 * 1024)
                print(f"  {description}: {size_mb:.0f} MB")

            # Progress bar
            progress = None
            if show_progress and not quiet and total_size > 0:
                from .download import _get_progress_bar

                progress = _get_progress_bar(total_size, f"  {description}")

            sha256 = hashlib.sha256()
            chunk_size = 256 * 1024  # 256 KB chunks for large files

            with os.fdopen(temp_fd, "wb") as f:
                temp_fd = -1
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    sha256.update(chunk)
                    if progress:
                        progress.update(len(chunk))

            if progress:
                progress.close()

        if validator is not None and not validator(temp_path):
            raise ValueError(
                f"Downloaded file {dest.name} failed verification - "
                "corrupt or incomplete"
            )

        # Atomic move to final destination (also restores 0644 permissions)
        from .download import publish_temp_file

        publish_temp_file(temp_path, dest)
        published = True

        if not quiet:
            actual_mb = dest.stat().st_size / (1024 * 1024)
            print(f"  Done ({actual_mb:.1f} MB, sha256: {sha256.hexdigest()[:16]}...)")

    # The cleanup is unconditional. Narrow except clauses would miss anything
    # a caller-supplied validator raises (it is arbitrary code), orphaning a
    # fully written .download file next to a multi-gigabyte data set.
    finally:
        if not published:
            # The descriptor is only handed to os.fdopen once the response is
            # open; a failure before that would otherwise leak it.
            if temp_fd != -1:
                try:
                    os.close(temp_fd)
                except OSError:
                    pass
            try:
                os.unlink(temp_path)
            except OSError:
                pass


def download_assist_data(
    target_dir: Optional[str] = None,
    planets: bool = True,
    asteroids: bool = True,
    planets_de441: bool = False,
    force: bool = False,
    show_progress: bool = True,
    quiet: bool = False,
) -> Path:
    """Download ASSIST ephemeris data files from JPL.

    Downloads the JPL planet ephemeris and optionally the asteroid
    perturber file needed by ASSIST for ephemeris-quality integration.

    Files are downloaded atomically (temp file + rename) with SSL
    certificate verification via certifi. Existing files are verified
    for integrity before skipping.

    Args:
        target_dir: Directory to save files (default: ~/.libephemeris/assist/).
        planets: Download planet ephemeris (~98 MB).
        asteroids: Download asteroid perturbers (~616 MB).
        planets_de441: Download the DE441 long planet ephemeris (~2.6 GB),
            needed only for deep-time (extended-tier) N-body generation.
        force: Re-download even if files exist and are valid.
        show_progress: Show download progress bar.
        quiet: Suppress all non-error output.

    Returns:
        Path to the data directory.

    Raises:
        urllib.error.URLError: If download fails.
        OSError: If file system operations fail.

    Example:
        >>> from libephemeris.rebound_integration import download_assist_data
        >>> data_dir = download_assist_data()
        >>> # Now ASSIST will find files automatically
    """
    if target_dir is None:
        data_dir = _ASSIST_DEFAULT_DIR
    else:
        data_dir = Path(target_dir)

    data_dir.mkdir(parents=True, exist_ok=True)

    if not quiet:
        print(f"ASSIST data directory: {data_dir}")

    _ASSIST_FILES = [
        (planets, _ASSIST_PLANETS_URL, "linux_p1550p2650.440", "Planet ephemeris"),
        (asteroids, _ASSIST_ASTEROIDS_URL, "sb441-n16.bsp", "Asteroid perturbers"),
        (
            planets_de441,
            _ASSIST_PLANETS_DE441_URL,
            "linux_m13000p17000.441",
            "Planet ephemeris (DE441, deep-time, ~2.6 GB)",
        ),
    ]

    downloaded = 0
    skipped = 0

    for enabled, url, filename, description in _ASSIST_FILES:
        if not enabled:
            continue

        dest = data_dir / filename

        if not force and dest.exists():
            if _verify_assist_file(dest):
                if not quiet:
                    size_mb = dest.stat().st_size / (1024 * 1024)
                    print(f"  {description} already exists ({size_mb:.1f} MB): {dest}")
                skipped += 1
                continue
            elif not quiet:
                # The existing file stays until a verified replacement is
                # ready: these downloads are 98 MB to 2.6 GB, and unlinking
                # first turns a failed repair into an expensive loss.
                print(f"  {description} exists but failed verification, re-downloading")

        _download_single_file(
            url=url,
            dest=dest,
            description=description,
            show_progress=show_progress,
            quiet=quiet,
            validator=_assist_download_validator(filename),
        )
        downloaded += 1

    if not quiet:
        if downloaded == 0 and skipped > 0:
            print("All files already present and verified.")
        elif downloaded > 0:
            print(f"Downloaded {downloaded} file(s).")

    # Reset the cached availability check so it picks up the new files
    reset_assist_data_cache()

    return data_dir


# Mean obliquity at J2000.0 — the same 23.4392911 deg constant spk.py
# uses for its ICRS <-> ecliptic J2000 rotations, so states converted
# here live in the identical ecliptic frame as the SPK/Keplerian paths.
_EPS_J2000_RAD: float = math.radians(23.4392911)


def _ecliptic_j2000_to_icrs(x: float, y: float, z: float) -> Tuple[float, float, float]:
    """Rotate a vector from ecliptic J2000 to equatorial ICRS."""
    ce = math.cos(_EPS_J2000_RAD)
    se = math.sin(_EPS_J2000_RAD)
    return (x, y * ce - z * se, y * se + z * ce)


def _icrs_to_ecliptic_j2000(x: float, y: float, z: float) -> Tuple[float, float, float]:
    """Rotate a vector from equatorial ICRS to ecliptic J2000."""
    ce = math.cos(_EPS_J2000_RAD)
    se = math.sin(_EPS_J2000_RAD)
    return (x, y * ce + z * se, -y * se + z * ce)


# ASSIST integrates the test particle in the gravity of the Sun, the planets and
# the major-asteroid perturbers baked into the small-body ephemeris (sb441-n16:
# Ceres, Pallas, Juno, Vesta, Hygiea, Iris, ...). When the body being propagated
# *is* one of those perturbers, ASSIST adds that perturber's own (near-singular)
# acceleration onto the coincident test particle, which wrecks the integration
# (the orbit goes retrograde / diverges). ASSIST's Python API can only toggle the
# whole ASTEROIDS force, so when a self-coincidence is detected we drop it for
# that single integration. The omitted effect — the *other* major asteroids on
# the target — is < ~0.01" over a year, far below the catastrophic self term.
_ASSIST_SELF_AU: float = 1e-4  # ~15000 km: far below inter-asteroid separations,
#                                far above the ~0 self-coincidence of a perturber.


def _assist_disable_self_perturber(extras, ephem, x, y, z, t_rel) -> Optional[str]:
    """If the test particle at barycentric ICRF (x, y, z) AU coincides with an
    ASSIST asteroid perturber at time ``t_rel`` (= jd_tt - ephem.jd_ref), remove
    the ASTEROIDS force from ``extras`` to avoid a self-interaction singularity.

    Returns the matched perturber name, or None if there is no coincidence.
    """
    try:
        from assist.ephem import ASSIST_BODY_IDS
    except Exception:  # noqa: BLE001 (older/newer assist without the table)
        return None
    hit = None
    for bid, name in ASSIST_BODY_IDS.items():
        if bid < 11:  # 0-10 are Sun..Pluto (real bodies, never the "self" case)
            continue
        try:
            a = ephem.get_particle(name, t_rel)
        except Exception:  # noqa: BLE001
            continue
        if (x - a.x) ** 2 + (y - a.y) ** 2 + (z - a.z) ** 2 < _ASSIST_SELF_AU**2:
            hit = name
            break
    if hit is not None and "ASTEROIDS" in extras.forces:
        import warnings

        extras.forces = [f for f in extras.forces if f != "ASTEROIDS"]
        warnings.warn(
            f"propagate: the body being propagated coincides with the ASSIST "
            f"asteroid perturber '{hit}', so the ASTEROIDS force was disabled for "
            f"this integration to avoid a self-interaction singularity (the "
            f'perturbations from the other major asteroids, < ~0.01", are '
            f"omitted).",
            stacklevel=3,
        )
    return hit


@dataclass
class PropagationResult:
    """Result of orbit propagation.

    Contains position and velocity vectors in heliocentric ecliptic J2000
    coordinates at the target time.

    Attributes:
        x, y, z: Position components in AU
        vx, vy, vz: Velocity components in AU/day
        jd_tt: Julian Date (TT) of the result
        ecliptic_lon: Ecliptic longitude in degrees
        ecliptic_lat: Ecliptic latitude in degrees
        distance: Heliocentric distance in AU
    """

    x: float
    y: float
    z: float
    vx: float
    vy: float
    vz: float
    jd_tt: float

    @property
    def ecliptic_lon(self) -> float:
        """Ecliptic longitude in degrees (0-360)."""
        return math.degrees(math.atan2(self.y, self.x)) % 360.0

    @property
    def ecliptic_lat(self) -> float:
        """Ecliptic latitude in degrees (-90 to +90)."""
        dist = self.distance
        if dist == 0.0:
            return 0.0
        # Clamp: with subnormal components z/dist can round above 1.0
        return math.degrees(math.asin(max(-1.0, min(1.0, self.z / dist))))

    @property
    def distance(self) -> float:
        """Heliocentric distance in AU."""
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    def to_tuple(self) -> Tuple[float, float, float, float, float, float]:
        """Return (x, y, z, vx, vy, vz) tuple."""
        return (self.x, self.y, self.z, self.vx, self.vy, self.vz)


def check_rebound_available() -> bool:
    """Check if REBOUND is installed and available.

    Returns:
        bool: True if REBOUND can be imported, False otherwise.
    """
    try:
        import rebound  # noqa: F401 (availability probe)

        return True
    except ImportError:
        return False


def check_assist_available() -> bool:
    """Check if ASSIST is installed and available.

    Returns:
        bool: True if ASSIST can be imported, False otherwise.
    """
    try:
        import assist  # noqa: F401 (availability probe)

        return True
    except ImportError:
        return False


# Cached result of ASSIST data file availability check.
# None = not checked yet, True/False = cached result.
_assist_data_available: Optional[bool] = None
_ASSIST_PROPAGATION_CACHE_MAX = 256
_assist_propagation_cache: OrderedDict[tuple, "PropagationResult"] = OrderedDict()
_assist_propagation_cache_lock = RLock()
_ASSIST_ASTEROID_PERTURBER_NAMES = {
    "camilla",
    "ceres",
    "cybele",
    "davida",
    "eunomia",
    "euphrosyne",
    "europa",
    "hygiea",
    "interamnia",
    "iris",
    "juno",
    "pallas",
    "psyche",
    "sylvia",
    "thisbe",
    "vesta",
}


def check_assist_data_available() -> bool:
    """Check if ASSIST is installed AND its required data files are present.

    This is the recommended check before attempting ASSIST integration in
    hot paths. The result is cached after the first call to avoid repeated
    filesystem probes.

    The cache is cleared by ``reset_assist_data_cache()`` or
    ``libephemeris.close()``.

    Returns:
        bool: True if ASSIST is importable and the planet ephemeris file
              can be located, False otherwise.
    """
    global _assist_data_available
    if _assist_data_available is not None:
        return _assist_data_available

    if not check_assist_available():
        _assist_data_available = False
        return False

    # Check that at least the planet ephemeris file exists
    config = AssistEphemConfig()
    _assist_data_available = config.planets_file is not None
    return _assist_data_available


def reset_assist_data_cache() -> None:
    """Reset the cached ASSIST data availability check.

    Call this after downloading ASSIST data files or changing
    ``ASSIST_DIR`` so the next ``check_assist_data_available()``
    re-probes the filesystem.
    """
    global _assist_data_available
    _assist_data_available = None
    with _assist_propagation_cache_lock:
        _assist_propagation_cache.clear()


def get_rebound_version() -> Optional[str]:
    """Get the installed REBOUND version.

    Returns:
        str: Version string if installed, None otherwise.
    """
    try:
        import rebound  # noqa: F401 (availability probe)

        return rebound.__version__
    except ImportError:
        return None


def get_assist_version() -> Optional[str]:
    """Get the installed ASSIST version.

    Returns:
        str: Version string if installed, None otherwise.
    """
    try:
        import assist  # noqa: F401 (availability probe)

        return assist.__version__
    except ImportError:
        return None


def elements_to_rebound_orbit(
    elements: "OrbitalElements",
    jd_tt: float,
) -> dict:
    """Convert OrbitalElements to REBOUND orbital parameters.

    REBOUND uses different conventions for orbital elements:
    - Semi-major axis in simulation units (typically AU)
    - Angles in radians (we convert from degrees)
    - Mean anomaly or true anomaly

    Args:
        elements: OrbitalElements from minor_bodies module
        jd_tt: Julian Date (TT) at which to evaluate elements

    Returns:
        dict: Parameters suitable for rebound.Particle.add(...)
    """
    from .minor_bodies import apply_secular_perturbations

    # Get perturbed elements at target time
    omega, Omega, M, n, e_pert, i_pert = apply_secular_perturbations(
        elements, jd_tt, include_perturbations=True
    )

    return {
        "a": elements.a,  # Semi-major axis in AU
        "e": e_pert,  # Perturbed eccentricity
        "inc": math.radians(i_pert),  # Perturbed inclination in radians
        "omega": math.radians(omega),  # Argument of perihelion in radians
        "Omega": math.radians(Omega),  # Longitude of ascending node in radians
        "M": math.radians(M),  # Mean anomaly in radians
    }


def _elements_to_cartesian(
    elements: "OrbitalElements",
    jd_tt: float,
) -> dict:
    """Convert OrbitalElements to Cartesian state vector for ASSIST.

    ASSIST manages the Sun internally, so particles must be added with
    Cartesian coordinates (x, y, z, vx, vy, vz) rather than orbital
    elements.  This function creates a temporary REBOUND simulation with
    a unit-mass Sun, adds the particle via orbital elements, extracts
    the resulting Cartesian state, and returns it.

    The coordinates are heliocentric ecliptic J2000, in AU and AU/day.

    Args:
        elements: OrbitalElements from minor_bodies module.
        jd_tt: Julian Date (TT) at which to evaluate elements.

    Returns:
        dict with keys ``x, y, z, vx, vy, vz`` suitable for
        ``rebound.Simulation.add(m=0, ...)``.
    """
    import rebound  # noqa: F401 (availability probe)

    # Gaussian gravitational constant squared -> G*M_sun in AU^3/day^2
    k_squared = 0.00029591220828559

    tmp = rebound.Simulation()
    tmp.units = ("day", "AU", "Msun")
    tmp.G = k_squared
    tmp.add(m=1.0)  # Sun

    orb = elements_to_rebound_orbit(elements, jd_tt)
    tmp.add(m=0.0, primary=tmp.particles[0], **orb)
    tmp.move_to_hel()

    p = tmp.particles[1]
    return {
        "x": p.x,
        "y": p.y,
        "z": p.z,
        "vx": p.vx,
        "vy": p.vy,
        "vz": p.vz,
    }


def propagate_orbit_rebound(
    elements: "OrbitalElements",
    jd_start: float,
    jd_end: float,
    integrator: ReboundIntegrator = ReboundIntegrator.IAS15,
    dt: Optional[float] = None,
) -> PropagationResult:
    """Propagate an orbit using REBOUND N-body integration.

    Uses REBOUND to integrate the orbit in the gravitational field of the Sun
    only (2-body problem). For ephemeris-quality precision including planetary
    perturbations, use propagate_orbit_assist() instead.

    Args:
        elements: OrbitalElements at starting epoch
        jd_start: Starting Julian Date (TT)
        jd_end: Ending Julian Date (TT)
        integrator: Which REBOUND integrator to use (default: IAS15)
        dt: Timestep in days for fixed-step integrators (default: auto)

    Returns:
        PropagationResult: Position and velocity at jd_end

    Raises:
        ImportError: If REBOUND is not installed
        ValueError: If elements are invalid

    Example:
        >>> from libephemeris.minor_bodies import MINOR_BODY_ELEMENTS, CERES
        >>> elements = MINOR_BODY_ELEMENTS[CERES]
        >>> result = propagate_orbit_rebound(elements, 2460000.5, 2460365.5)
        >>> print(f"Ceres at {result.ecliptic_lon:.4f} deg, {result.distance:.4f} AU")
    """
    try:
        import rebound  # noqa: F401 (availability probe)
    except ImportError:
        raise ImportError(
            "REBOUND is required for n-body orbit propagation. "
            "Install with: pip install rebound"
        )

    # Gaussian gravitational constant squared (k^2)
    # This gives G*M_sun in AU^3/day^2
    k_squared = 0.00029591220828559  # AU^3/day^2

    # Create simulation
    sim = rebound.Simulation()
    sim.units = ("day", "AU", "Msun")
    sim.G = k_squared  # GM_sun = k^2 for Msun = 1

    # Set integrator
    sim.integrator = integrator.value

    # Set timestep for fixed-step integrators
    if dt is not None:
        sim.dt = dt
    elif integrator in (ReboundIntegrator.WHFAST, ReboundIntegrator.LEAPFROG):
        # Default to ~1/20 of orbital period
        period = 365.25 * elements.a**1.5  # Approximate period in days
        sim.dt = period / 20.0

    # Add the Sun at origin
    sim.add(m=1.0, hash="Sun")

    # Convert orbital elements to REBOUND format
    orb_params = elements_to_rebound_orbit(elements, jd_start)

    # Add test particle with orbital elements
    sim.add(
        m=0.0,  # Test particle (no mass)
        primary=sim.particles["Sun"],
        hash=elements.name,
        **orb_params,
    )

    # Move to heliocentric frame (Sun at origin)
    sim.move_to_hel()

    # Set initial time (relative, simulation starts at t=0)
    sim.t = 0.0

    # Integration time in days
    integration_time = jd_end - jd_start

    # Integrate
    sim.integrate(integration_time)

    # Get particle state
    p = sim.particles[1]  # The asteroid (index 1, after Sun)

    return PropagationResult(
        x=p.x,
        y=p.y,
        z=p.z,
        vx=p.vx,
        vy=p.vy,
        vz=p.vz,
        jd_tt=jd_end,
    )


def propagate_orbit_assist(
    elements: "OrbitalElements",
    jd_start: float,
    jd_end: float,
    ephem_config: Optional[AssistEphemConfig] = None,
    include_non_gravitational: bool = False,
    A1: float = 0.0,
    A2: float = 0.0,
    A3: float = 0.0,
) -> PropagationResult:
    """Propagate an orbit using ASSIST for ephemeris-quality integration.

    ASSIST provides the highest precision orbit propagation by including:
    - Sun, Moon, 8 planets from JPL DE440/DE441 ephemeris
    - 16 massive asteroids (Ceres, Vesta, Pallas, etc.)
    - Earth and Sun gravitational harmonics (J2, J3, J4)
    - General relativistic corrections
    - Optional non-gravitational forces (Marsden model)

    The precision matches JPL's small body integrator to machine precision.

    Args:
        elements: OrbitalElements at starting epoch
        jd_start: Starting Julian Date (TT)
        jd_end: Ending Julian Date (TT)
        ephem_config: Configuration for ephemeris data files
        include_non_gravitational: Include Marsden non-gravitational forces
        A1, A2, A3: Marsden parameters for non-gravitational acceleration
            (only used if include_non_gravitational=True)

    Returns:
        PropagationResult: Position and velocity at jd_end

    Raises:
        ImportError: If ASSIST or REBOUND is not installed
        FileNotFoundError: If ephemeris files are not found
        ValueError: If elements are invalid

    Example:
        >>> from libephemeris.minor_bodies import MINOR_BODY_ELEMENTS, CERES
        >>> from libephemeris.rebound_integration import propagate_orbit_assist
        >>> elements = MINOR_BODY_ELEMENTS[CERES]
        >>> result = propagate_orbit_assist(elements, 2460000.5, 2460365.5)
        >>> print(f"Ceres: lon={result.ecliptic_lon:.6f} deg")

    Note:
        Requires ~1 GB of ephemeris data files. See module docstring for
        download instructions.

        jd_start/jd_end are TT Julian Dates and are passed to ASSIST, which
        expects TDB. TT is used directly as an approximation of TDB; the
        TT-TDB difference is <= ~1.7 ms and is negligible for these
        integrations.
    """
    try:
        import rebound  # noqa: F401 (availability probe)
    except ImportError:
        raise ImportError(
            "REBOUND is required for ASSIST. Install with: pip install rebound"
        )

    try:
        import assist  # noqa: F401 (availability probe)
    except ImportError:
        raise ImportError(
            "ASSIST is required for ephemeris-quality integration. "
            "Install with: pip install assist\n"
            "See also: https://assist.readthedocs.io/"
        )

    cache_key = None
    if (
        ephem_config is None
        and not include_non_gravitational
        and A1 == 0.0
        and A2 == 0.0
        and A3 == 0.0
        and elements.name.lower() not in _ASSIST_ASTEROID_PERTURBER_NAMES
    ):
        # Key on the element VALUES, not just the name: two different orbits
        # sharing a name (and window) must not alias to the same slot.
        cache_key = (
            elements.name,
            float(elements.epoch),
            float(elements.a),
            float(elements.e),
            float(elements.i),
            float(elements.omega),
            float(elements.Omega),
            float(elements.M0),
            float(elements.n),
            float(jd_start),
            round(float(jd_end) * 86400.0),
        )
        with _assist_propagation_cache_lock:
            cached = _assist_propagation_cache.get(cache_key)
            if cached is not None:
                _assist_propagation_cache.move_to_end(cache_key)
                return cached

    # Configure ephemeris paths
    if ephem_config is None:
        ephem_config = AssistEphemConfig()

    # Verify ephemeris files exist
    if ephem_config.planets_file is None:
        raise FileNotFoundError(
            "Planet ephemeris file not found. Download from:\n"
            "curl https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/linux_p1550p2650.440 -o data/linux_p1550p2650.440\n"
            "Set ASSIST_DIR environment variable or pass ephem_config."
        )

    if not os.path.exists(ephem_config.planets_file):
        raise FileNotFoundError(
            f"Planet ephemeris file not found: {ephem_config.planets_file}"
        )

    if ephem_config.asteroids_file and not os.path.exists(ephem_config.asteroids_file):
        raise FileNotFoundError(
            f"Asteroid ephemeris file not found: {ephem_config.asteroids_file}"
        )

    # Create ASSIST ephemeris object
    if ephem_config.asteroids_file:
        ephem = assist.Ephem(ephem_config.planets_file, ephem_config.asteroids_file)
    else:
        ephem = assist.Ephem(ephem_config.planets_file)

    # Create REBOUND simulation
    sim = rebound.Simulation()

    # Attach ASSIST extras
    extras = assist.Extras(sim, ephem)

    # Convert orbital elements to heliocentric ecliptic-J2000 Cartesian
    # state.  ASSIST integrates in the *equatorial ICRF frame with the
    # solar-system barycenter at the origin* (time = TDB days relative
    # to ephem.jd_ref), so rotate ecliptic -> equatorial and shift the
    # origin by the Sun's barycentric state from ASSIST's own ephemeris.
    cart = _elements_to_cartesian(elements, jd_start)

    sun0 = ephem.get_particle("sun", jd_start - ephem.jd_ref)
    x0, y0, z0 = _ecliptic_j2000_to_icrs(cart["x"], cart["y"], cart["z"])
    vx0, vy0, vz0 = _ecliptic_j2000_to_icrs(cart["vx"], cart["vy"], cart["vz"])

    # Add test particle with barycentric equatorial state
    sim.add(
        m=0.0,
        x=x0 + sun0.x,
        y=y0 + sun0.y,
        z=z0 + sun0.z,
        vx=vx0 + sun0.vx,
        vy=vy0 + sun0.vy,
        vz=vz0 + sun0.vz,
    )

    # Avoid the self-interaction singularity when the propagated body is itself
    # one of ASSIST's asteroid perturbers (Ceres/Vesta/Pallas/Juno/...).
    _assist_disable_self_perturber(
        extras,
        ephem,
        x0 + sun0.x,
        y0 + sun0.y,
        z0 + sun0.z,
        jd_start - ephem.jd_ref,
    )

    # Configure non-gravitational forces if requested.  ASSIST's API
    # takes a flat array of 3 Marsden parameters (A1 radial, A2
    # transverse, A3 normal) per particle, after the particle exists.
    # The setter requires a numpy float64 array (it calls value.ctypes on
    # it); a plain Python list raises AttributeError.
    if include_non_gravitational:
        import numpy as np

        extras.particle_params = np.array([A1, A2, A3], dtype=np.float64)

    # Set initial time (JD - reference JD)
    sim.t = jd_start - ephem.jd_ref

    # Integration time
    t_end = jd_end - ephem.jd_ref

    # Integrate
    sim.integrate(t_end)

    # Get particle state (barycentric equatorial ICRF) and convert back
    # to the heliocentric ecliptic J2000 frame PropagationResult uses.
    p = sim.particles[0]
    sun1 = ephem.get_particle("sun", t_end)
    xh, yh, zh = _icrs_to_ecliptic_j2000(p.x - sun1.x, p.y - sun1.y, p.z - sun1.z)
    vxh, vyh, vzh = _icrs_to_ecliptic_j2000(
        p.vx - sun1.vx, p.vy - sun1.vy, p.vz - sun1.vz
    )

    result = PropagationResult(
        x=xh,
        y=yh,
        z=zh,
        vx=vxh,
        vy=vyh,
        vz=vzh,
        jd_tt=jd_end,
    )
    if cache_key is not None:
        with _assist_propagation_cache_lock:
            _assist_propagation_cache[cache_key] = result
            _assist_propagation_cache.move_to_end(cache_key)
            while len(_assist_propagation_cache) > _ASSIST_PROPAGATION_CACHE_MAX:
                _assist_propagation_cache.popitem(last=False)
    return result


def propagate_trajectory(
    elements: "OrbitalElements",
    jd_start: float,
    jd_end: float,
    num_points: int = 100,
    use_assist: bool = True,
    ephem_config: Optional[AssistEphemConfig] = None,
    integrator: ReboundIntegrator = ReboundIntegrator.IAS15,
) -> List[PropagationResult]:
    """Propagate an orbit and return positions at multiple times.

    Useful for generating ephemerides or plotting orbital trajectories.

    Args:
        elements: OrbitalElements at starting epoch
        jd_start: Starting Julian Date (TT)
        jd_end: Ending Julian Date (TT)
        num_points: Number of output points (evenly spaced in time)
        use_assist: Use ASSIST for ephemeris-quality (if available)
        ephem_config: ASSIST ephemeris configuration
        integrator: REBOUND integrator (only used if use_assist=False)

    Returns:
        List[PropagationResult]: Positions and velocities at each time

    Note:
        When use_assist=True, jd_start/jd_end (TT) are passed to ASSIST,
        which expects TDB. TT is used directly as an approximation of TDB;
        the TT-TDB difference is <= ~1.7 ms and is negligible here.

    Example:
        >>> trajectory = propagate_trajectory(elements, jd_start, jd_end, num_points=365)
        >>> for point in trajectory:
        ...     print(f"JD {point.jd_tt}: lon={point.ecliptic_lon:.4f}")
    """
    results = []

    # Time step between output points
    dt = (jd_end - jd_start) / (num_points - 1) if num_points > 1 else 0

    if use_assist and check_assist_available():
        try:
            import rebound  # noqa: F401 (availability probe)
            import assist  # noqa: F401 (availability probe)

            # Configure ephemeris
            if ephem_config is None:
                ephem_config = AssistEphemConfig()

            if ephem_config.planets_file and os.path.exists(ephem_config.planets_file):
                # Use ASSIST
                if ephem_config.asteroids_file:
                    ephem = assist.Ephem(
                        ephem_config.planets_file, ephem_config.asteroids_file
                    )
                else:
                    ephem = assist.Ephem(ephem_config.planets_file)

                sim = rebound.Simulation()
                extras = assist.Extras(sim, ephem)

                # ASSIST integrates in the equatorial ICRF frame with the
                # solar-system barycenter at the origin (same convention as
                # propagate_orbit_assist), so rotate ecliptic -> equatorial
                # and shift the origin by the Sun's barycentric state before
                # adding the particle, then undo both at each output epoch.
                cart = _elements_to_cartesian(elements, jd_start)
                sun0 = ephem.get_particle("sun", jd_start - ephem.jd_ref)
                x0, y0, z0 = _ecliptic_j2000_to_icrs(cart["x"], cart["y"], cart["z"])
                vx0, vy0, vz0 = _ecliptic_j2000_to_icrs(
                    cart["vx"], cart["vy"], cart["vz"]
                )
                sim.add(
                    m=0.0,
                    x=x0 + sun0.x,
                    y=y0 + sun0.y,
                    z=z0 + sun0.z,
                    vx=vx0 + sun0.vx,
                    vy=vy0 + sun0.vy,
                    vz=vz0 + sun0.vz,
                )

                # Avoid the self-interaction singularity when the body is itself
                # one of ASSIST's asteroid perturbers (Ceres/Vesta/Pallas/...).
                # Evaluated once before the loop on purpose: the match is by the
                # identity of the propagated body versus the fixed ASSIST
                # perturber set, which does not change along the trajectory, so a
                # per-step re-check would be redundant.
                _assist_disable_self_perturber(
                    extras,
                    ephem,
                    x0 + sun0.x,
                    y0 + sun0.y,
                    z0 + sun0.z,
                    jd_start - ephem.jd_ref,
                )

                sim.t = jd_start - ephem.jd_ref

                for i in range(num_points):
                    t = jd_start + i * dt
                    sim.integrate(t - ephem.jd_ref)
                    p = sim.particles[0]
                    sun_t = ephem.get_particle("sun", t - ephem.jd_ref)
                    xh, yh, zh = _icrs_to_ecliptic_j2000(
                        p.x - sun_t.x, p.y - sun_t.y, p.z - sun_t.z
                    )
                    vxh, vyh, vzh = _icrs_to_ecliptic_j2000(
                        p.vx - sun_t.vx, p.vy - sun_t.vy, p.vz - sun_t.vz
                    )
                    results.append(
                        PropagationResult(
                            x=xh, y=yh, z=zh, vx=vxh, vy=vyh, vz=vzh, jd_tt=t
                        )
                    )

                return results
        except (ImportError, FileNotFoundError, RuntimeError, ValueError):
            # Fall back to REBOUND. Besides a missing module/file, ASSIST and
            # rebound raise RuntimeError/ValueError when the integration epoch
            # falls outside the loaded planet-ephemeris date range; the 2-body
            # REBOUND path below needs no ephemeris file and can still succeed.
            pass

    # Use REBOUND without ASSIST
    if not check_rebound_available():
        raise ImportError(
            "REBOUND is required for orbit propagation. "
            "Install with: pip install rebound"
        )

    import rebound  # noqa: F401 (availability probe)

    k_squared = 0.00029591220828559
    sim = rebound.Simulation()
    sim.units = ("day", "AU", "Msun")
    sim.G = k_squared
    sim.integrator = integrator.value

    sim.add(m=1.0, hash="Sun")

    orb_params = elements_to_rebound_orbit(elements, jd_start)
    sim.add(m=0.0, primary=sim.particles["Sun"], **orb_params)
    sim.move_to_hel()
    sim.t = 0.0

    for i in range(num_points):
        t = i * dt
        sim.integrate(t)
        p = sim.particles[1]
        results.append(
            PropagationResult(
                x=p.x, y=p.y, z=p.z, vx=p.vx, vy=p.vy, vz=p.vz, jd_tt=jd_start + t
            )
        )

    return results


def compare_with_keplerian(
    elements: "OrbitalElements",
    jd_tt: float,
    use_assist: bool = False,
    ephem_config: Optional[AssistEphemConfig] = None,
) -> dict:
    """Compare REBOUND/ASSIST integration with Keplerian approximation.

    Useful for evaluating the improvement in precision from using numerical
    integration vs. the analytical Keplerian method with secular perturbations.

    Args:
        elements: OrbitalElements to propagate
        jd_tt: Target Julian Date (TT)
        use_assist: Use ASSIST if available (else REBOUND only)
        ephem_config: ASSIST ephemeris configuration

    Returns:
        dict: Comparison results including positions and differences

    Example:
        >>> comparison = compare_with_keplerian(elements, jd_tt)
        >>> print(f"Position difference: {comparison['angular_sep_arcsec']:.2f} arcsec")
    """
    from .minor_bodies import calc_minor_body_position

    # Keplerian result
    x_kep, y_kep, z_kep = calc_minor_body_position(elements, jd_tt)
    r_kep = math.sqrt(x_kep**2 + y_kep**2 + z_kep**2)
    lon_kep = math.degrees(math.atan2(y_kep, x_kep)) % 360.0
    # Guard the origin and clamp asin (r >= |z| but rounding can nudge past
    # +-1), mirroring PropagationResult.ecliptic_lat.
    if r_kep > 0.0:
        lat_kep = math.degrees(math.asin(max(-1.0, min(1.0, z_kep / r_kep))))
    else:
        lat_kep = 0.0

    # REBOUND/ASSIST result
    jd_start = elements.epoch

    if use_assist and check_assist_available():
        try:
            result = propagate_orbit_assist(elements, jd_start, jd_tt, ephem_config)
        except (ImportError, FileNotFoundError, RuntimeError, ValueError):
            # Same ASSIST->REBOUND fallback as propagate_trajectory: an epoch
            # outside the loaded planet-ephemeris range raises RuntimeError/
            # ValueError, which the 2-body REBOUND path can still handle.
            result = propagate_orbit_rebound(elements, jd_start, jd_tt)
    else:
        result = propagate_orbit_rebound(elements, jd_start, jd_tt)

    # Angular separation
    d_lon = result.ecliptic_lon - lon_kep
    if d_lon > 180:
        d_lon -= 360
    elif d_lon < -180:
        d_lon += 360

    d_lat = result.ecliptic_lat - lat_kep

    # Approximate angular separation in arcseconds
    cos_lat = math.cos(math.radians((lat_kep + result.ecliptic_lat) / 2))
    angular_sep = math.sqrt((d_lon * cos_lat) ** 2 + d_lat**2) * 3600  # arcsec

    return {
        "keplerian": {
            "x": x_kep,
            "y": y_kep,
            "z": z_kep,
            "lon": lon_kep,
            "lat": lat_kep,
            "dist": r_kep,
        },
        "nbody": {
            "x": result.x,
            "y": result.y,
            "z": result.z,
            "lon": result.ecliptic_lon,
            "lat": result.ecliptic_lat,
            "dist": result.distance,
        },
        "delta_lon_arcsec": d_lon * 3600,
        "delta_lat_arcsec": d_lat * 3600,
        "angular_sep_arcsec": angular_sep,
        "delta_dist_au": result.distance - r_kep,
        "method": "assist" if use_assist else "rebound",
        "propagation_days": jd_tt - jd_start,
    }
