# Precision Tuning Guide

This guide explains how to achieve maximum precision in LibEphemeris using optional dependencies and configuration options. LibEphemeris provides several mechanisms to enhance calculation precision beyond the default settings.

## Table of Contents

1. [Overview](#overview)
2. [SPK Kernels for Minor Bodies](#spk-kernels-for-minor-bodies)
3. [IERS Delta T Data](#iers-delta-t-data)
4. [Ephemeris File Selection](#ephemeris-file-selection)
5. [Tidal Acceleration](#tidal-acceleration)
6. [Automatic SPK Download](#automatic-spk-download)
7. [Configuration Summary](#configuration-summary)
8. [Best Practices](#best-practices)

---

## Overview

LibEphemeris offers different precision levels depending on your requirements:

| Feature | Default Precision | Enhanced Precision | How to Enable |
|---------|-------------------|-------------------|---------------|
| Planets | Sub-arcsecond | Sub-arcsecond | Built-in with DE440 |
| Minor Bodies | ~1-10 arcminutes | ~1-5 arcseconds | SPK kernels |
| Delta T | Skyfield model | IERS observed | IERS data download |
| Historical dates | DE440 (1550-2650) | Extended range | DE441/DE431 |

---

## SPK Kernels for Minor Bodies

By default, minor bodies (asteroids, centaurs, TNOs) are calculated using Keplerian orbital elements with secular perturbations. For significantly higher precision, you can use SPK kernels from JPL.

### Precision Comparison

| Method | Typical Accuracy | Use Case |
|--------|------------------|----------|
| Keplerian (default) | ~1-10 arcminutes | Quick estimates, historical charts |
| SPK kernel | ~1-5 arcseconds | Precise calculations, transit timing |

### Manual SPK Download and Registration

```python
import libephemeris as eph

# Download SPK for Chiron covering 2000-2100
spk_path = eph.download_spk(
    body="Chiron",
    start="2000-01-01",
    end="2100-01-01",
    directory="./spk_kernels"
)

# Register the SPK for Chiron calculations
eph.register_spk_body(
    ipl=eph.CHIRON,
    spk_path=spk_path,
    naif_id=eph.NAIF_CHIRON  # 2002060
)

# Now calc_ut automatically uses the SPK kernel
pos, _ = eph.calc_ut(2451545.0, eph.CHIRON, eph.FLG_SPEED)
print(f"Chiron (SPK): {pos[0]:.6f}°")
```

### One-Liner Convenience Function

```python
import libephemeris as eph

# Download and register in one call
eph.download_and_register_spk(
    body="Eris",
    ipl=eph.ERIS,
    start="2000-01-01",
    end="2100-01-01",
)

# Calculations now use SPK automatically
pos, _ = eph.calc_ut(2460000.0, eph.ERIS, 0)
```

### NAIF ID Constants

SPK kernels use NASA NAIF IDs. LibEphemeris provides constants for common bodies:

```python
# For numbered asteroids: NAIF_ID = asteroid_number + 2000000
NAIF_ASTEROID_OFFSET = 2000000
NAIF_CHIRON = 2002060    # Chiron (2060)
NAIF_CERES = 2000001     # Ceres (1)
NAIF_PALLAS = 2000002    # Pallas (2)
NAIF_JUNO = 2000003      # Juno (3)
NAIF_VESTA = 2000004     # Vesta (4)
NAIF_PHOLUS = 2005145    # Pholus (5145)
NAIF_NESSUS = 2007066    # Nessus (7066)
NAIF_ERIS = 2136199      # Eris (136199)
NAIF_SEDNA = 2090377     # Sedna (90377)
```

### SPK Management Functions

| Function | Description |
|----------|-------------|
| `download_spk(body, start, end, ...)` | Download SPK from JPL Horizons |
| `register_spk_body(ipl, spk_path, naif_id)` | Register SPK for a body |
| `unregister_spk_body(ipl)` | Remove SPK registration |
| `download_and_register_spk(...)` | Download and register in one call |
| `list_spk_bodies()` | List all registered SPK bodies |
| `get_spk_body_info(ipl)` | Get SPK info for a body |
| `get_spk_coverage(spk_path)` | Get date range covered by SPK |

### Thread-Safe SPK Usage

```python
from libephemeris import EphemerisContext, CHIRON

ctx = EphemerisContext()

# Register SPK in this context only
ctx.register_spk_body(CHIRON, "./chiron.bsp", 2002060)

pos, _ = ctx.calc_ut(2451545.0, CHIRON, 0)
```

---

## IERS Delta T Data

Delta T (TT - UT1) accounts for Earth's irregular rotation. For maximum precision on recent dates (1973-present), use observed values from IERS (International Earth Rotation and Reference Systems Service).

### Delta T Sources

| Source | Precision | Date Range | Notes |
|--------|-----------|------------|-------|
| Skyfield model | ~3.6 seconds | All dates | Default, based on historical/predicted values |
| IERS observed | ~0.1 seconds | 1973-present | Highest precision for recent dates |
| User-defined | Exact | All dates | For testing or special cases |

### Enabling IERS Delta T

```python
from libephemeris import set_iers_delta_t_enabled
from libephemeris.iers_data import download_delta_t_data

# Download IERS data (cached locally)
download_delta_t_data()

# Enable IERS Delta T
set_iers_delta_t_enabled(True)

# Now deltat() uses IERS observed values for recent dates
```

### Environment Variable Configuration

```bash
# Enable IERS Delta T via environment variable
export LIBEPHEMERIS_IERS_DELTA_T=1
```

### IERS Data Management

```python
from libephemeris.iers_data import (
    download_delta_t_data,
    download_iers_finals,
    download_leap_seconds,
    get_iers_cache_info,
    get_observed_delta_t,
    is_observed_delta_t_available,
    set_iers_cache_dir,
    set_iers_auto_download,
)

# Set custom cache directory
set_iers_cache_dir("/path/to/iers_cache")

# Enable automatic IERS data download
set_iers_auto_download(True)

# Check cache status
info = get_iers_cache_info()
print(f"Cache directory: {info['cache_dir']}")
print(f"Delta T entries: {info['delta_t_entries']}")

# Check if observed Delta T is available for a date
jd = 2451545.0  # J2000.0
if is_observed_delta_t_available(jd):
    delta_t = get_observed_delta_t(jd)
    print(f"Observed Delta T: {delta_t:.4f} seconds")
```

### User-Defined Delta T

For special cases (testing, very ancient dates, future predictions):

```python
from libephemeris import set_delta_t_userdef, get_delta_t_userdef, deltat

# Set a fixed Delta T of 65 seconds (in days)
set_delta_t_userdef(65.0 / 86400.0)

# All Delta T calculations now return this value
dt = deltat(2451545.0)
print(f"Delta T: {dt * 86400:.1f} seconds")  # 65.0 seconds

# Clear to resume computed values
set_delta_t_userdef(None)
```

---

## Ephemeris File Selection

Different JPL Development Ephemeris (DE) files provide different date ranges and precision levels.

### Available Ephemeris Files

| File | Date Range | Size | Notes |
|------|------------|------|-------|
| **de440s.bsp** | 1849-2150 | ~31 MB | Lightweight subset of DE440 |
| **de440.bsp** | 1550-2650 | ~114 MB | **Default**, ICRF 3.0 |
| **de441.bsp** | -13200-17191 | ~3.4 GB | Extended version of DE440 |

DE440 and DE441 have identical precision — DE441 is simply the extended-range version.
DE440s is a reduced-size subset of DE440 with no loss of precision within its range.

### Precision Tiers

LibEphemeris organizes the three current-generation files into precision tiers:

| Tier | File | Use Case |
|------|------|----------|
| `base` | de440s.bsp | Lightweight, modern-era usage (1849-2150) |
| `medium` | de440.bsp | General purpose **(DEFAULT)** (1550-2650) |
| `extended` | de441.bsp | Historical/far-future research (-13200 to +17191) |

```python
from libephemeris import set_precision_tier

set_precision_tier("extended")   # uses de441.bsp
```

```bash
export LIBEPHEMERIS_PRECISION=extended
```

### Selecting an Ephemeris File

You can also select a specific file directly, bypassing the tier system:

```python
from libephemeris import set_ephemeris_file, set_ephe_path

# Set custom ephemeris directory (optional)
set_ephe_path("/path/to/jpl-kernels")

# Select ephemeris file
set_ephemeris_file("de441.bsp")  # For extended date range
```

### Resolution Priority

1. `LIBEPHEMERIS_EPHEMERIS` environment variable
2. `set_ephemeris_file()` / `set_jpl_file()` programmatic call
3. Precision tier (`LIBEPHEMERIS_PRECISION` env var or `set_precision_tier()`)
4. Default: `de440.bsp` (medium tier)

### When to Use Different Files

- **de440s.bsp**: Lightweight option for modern dates (1849-2150), same precision as DE440
- **de440.bsp (default)**: Best for most modern applications (1550-2650)
- **de441.bsp**: For very ancient or far future dates

---

## Tidal Acceleration

The tidal acceleration affects Delta T calculations for dates far from the present. Different ephemeris files use different values.

### Tidal Acceleration Values

| Ephemeris | Tidal Acceleration |
|-----------|-------------------|
| DE421 | -25.85 arcsec/cy² |
| DE430 | -25.82 arcsec/cy² |
| DE431 | -25.80 arcsec/cy² |
| DE440/DE441 | -25.936 arcsec/cy² |

### Setting Tidal Acceleration

```python
from libephemeris import set_tid_acc, get_tid_acc
from libephemeris.constants import TIDAL_DE440

# Match tidal acceleration to ephemeris file
set_tid_acc(TIDAL_DE440)  # When using de440.bsp (default)
print(f"Tidal acceleration: {get_tid_acc()}")

# Use default (DE440-based)
set_tid_acc(0.0)  # TIDAL_AUTOMATIC
```

---

## Automatic SPK Download

For convenience, LibEphemeris can automatically download SPK kernels on demand using `astroquery`.

### Installation

```bash
pip install astroquery
```

### Enabling Automatic SPK Download

```python
import libephemeris as eph

# Enable automatic SPK download
eph.set_auto_spk_download(True)

# Set cache directory (optional)
eph.set_spk_cache_dir("./spk_cache")

# Set date padding (optional) - extends download range by N days
eph.set_spk_date_padding(365)  # Add 1 year buffer
```

### Environment Variable Configuration

```bash
export LIBEPHEMERIS_AUTO_SPK=1
```

### Using auto_get_spk() for On-Demand Downloads

```python
from libephemeris.spk_auto import auto_get_spk
from libephemeris.constants import CHIRON

# Automatically download and register SPK for a date range
jd_start = 2458849.5  # 2020-01-01
jd_end = 2462502.5    # 2030-01-01

spk_path = auto_get_spk(
    body_id="2060",  # Chiron
    jd_start=jd_start,
    jd_end=jd_end,
    ipl=CHIRON,  # Auto-register for calc_ut
)

# calc_ut now uses SPK data automatically
pos, _ = eph.calc_ut(2460000.0, CHIRON, 0)
```

### Enabling Common Bodies at Once

```python
from libephemeris.spk_auto import enable_common_bodies

# Enable auto-SPK for popular minor bodies
enable_common_bodies(
    start="2000-01-01",
    end="2100-01-01",
)
# Enables: Chiron, Pholus, Ceres, Pallas, Juno, Vesta, Eris, Sedna
```

### Cache Management

```python
from libephemeris.spk_auto import (
    list_cached_spk,
    get_cache_size,
    clear_spk_cache,
    prune_old_cache,
)

# List cached SPK files
for spk in list_cached_spk():
    print(f"{spk['filename']}: {spk['size_mb']:.2f} MB")
    if spk['date_start']:
        print(f"  Coverage: {spk['date_start']} to {spk['date_end']}")

# Check cache size
print(f"Total cache size: {get_cache_size():.2f} MB")

# Remove files not accessed in 30 days
pruned = prune_old_cache(max_age_days=30)
print(f"Removed {pruned} old cache files")

# Clear all cached SPK files
clear_spk_cache()
```

---

## Configuration Summary

### Global State Configuration

| Function | Purpose |
|----------|---------|
| `set_precision_tier(tier)` | Select precision tier (`"base"`, `"medium"`, `"extended"`) |
| `get_precision_tier()` | Get current precision tier name |
| `set_ephemeris_file(filename)` | Select JPL DE ephemeris file (overrides tier) |
| `set_ephe_path(path)` | Set ephemeris file directory |
| `set_tid_acc(value)` | Set tidal acceleration for Delta T |
| `set_delta_t_userdef(dt)` | Set user-defined Delta T value |
| `set_iers_delta_t_enabled(True)` | Enable IERS observed Delta T |
| `set_auto_spk_download(True)` | Enable automatic SPK downloads |
| `set_spk_cache_dir(path)` | Set SPK cache directory |
| `set_spk_date_padding(days)` | Set date padding for SPK downloads |

### Environment Variables

| Variable | Purpose | Values |
|----------|---------|--------|
| `LIBEPHEMERIS_EPHEMERIS` | Select ephemeris file (highest priority) | e.g., `de441.bsp` |
| `LIBEPHEMERIS_PRECISION` | Select precision tier | `base`, `medium`, `extended` |
| `LIBEPHEMERIS_IERS_DELTA_T` | Enable IERS Delta T | `1`, `true`, `yes` |
| `LIBEPHEMERIS_AUTO_SPK` | Enable auto SPK download | `1`, `true`, `yes` |
| `LIBEPHEMERIS_IERS_AUTO_DOWNLOAD` | Auto-download IERS data | `1`, `true`, `yes` |

---

## Best Practices

### For Maximum Precision (Modern Dates)

```python
import libephemeris as eph
from libephemeris.iers_data import download_delta_t_data, set_iers_auto_download

# 1. Enable IERS Delta T for recent dates
set_iers_auto_download(True)
download_delta_t_data()
eph.set_iers_delta_t_enabled(True)

# 2. Use SPK kernels for minor bodies
eph.set_auto_spk_download(True)
eph.set_spk_cache_dir("./spk_cache")
eph.set_spk_date_padding(365)

# 3. Use default DE440 ephemeris (already default)
# eph.set_ephemeris_file("de440.bsp")
```

### For Historical Calculations (Ancient Dates)

```python
import libephemeris as eph

# 1. Use extended ephemeris file
eph.set_ephemeris_file("de441.bsp")

# 2. Set matching tidal acceleration
eph.set_tid_acc(eph.TIDAL_DE440)  # DE441 uses same as DE440

# 3. Note: IERS Delta T not available before 1973
# Skyfield model is used automatically for historical dates
```

### For Offline/Reproducible Calculations

```python
import libephemeris as eph

# Disable automatic downloads
eph.set_auto_spk_download(False)
eph.set_iers_delta_t_enabled(False)

# Use local ephemeris files only
eph.set_ephe_path("/path/to/local/ephemeris")
eph.set_ephemeris_file("de440.bsp")

# Optionally set fixed Delta T for exact reproducibility
eph.set_delta_t_userdef(69.0 / 86400.0)  # Fixed 69 seconds
```

### Memory and Performance Considerations

```python
import libephemeris as eph

# Check current ephemeris info
path, start, end, denum = eph.get_current_file_data()
print(f"Using DE{denum}: JD {start:.1f} to {end:.1f}")

# Close files when done (optional, for long-running apps)
eph.close()

# SPK files are cached - consider pruning old files
from libephemeris.spk_auto import prune_old_cache
prune_old_cache(max_age_days=90)
```

---

## References

1. NASA JPL Horizons System: https://ssd.jpl.nasa.gov/horizons/
2. IERS Earth Orientation Data: https://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html
3. JPL Development Ephemeris: https://ssd.jpl.nasa.gov/planets/eph_export.html
4. NAIF SPICE: https://naif.jpl.nasa.gov/naif/
5. Skyfield Documentation: https://rhodesmill.org/skyfield/
