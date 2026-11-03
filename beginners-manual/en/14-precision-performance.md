# Chapter 14 — Precision, Performance, and LEB Mode

[Contents and setup](README.md) · [Previous](13-hypothetical-planets-arabic-parts.md) · [Next](15-cookbook.md) · [Italiano](../it/14-precisione-performance.md)

## What You Will Learn

In this chapter, you will discover how precise LibEphemeris is compared to reference astronomical sources, how precision behaves over different time scales (from antiquity to the distant future), what the LEB mode (precomputed binary ephemerides) is and how to optimize calculation speed, how to configure calculation modes and precision tiers, and how to manage data downloading and caching.

---

## 14.1 How precise is this library?

A legitimate question before entrusting your calculations to any software: how much can we trust the results?

The useful answer starts with **which body, date, data and observable?** An arcsecond is 1/3600 of a degree, but a numerical error below that value is not a universal physical-accuracy guarantee. Source uncertainty, Earth rotation, observer coordinates and atmospheric conditions can matter more. A zodiacal longitude, an occultation contact and naked-eye visibility have different error budgets.

Let's look at the details by category.

### Main Planets

LibEphemeris reads **JPL DE440/DE441** or local LEB coefficients. JPL ephemerides are numerical orbital solutions fitted to observations, not exact trajectories. LEB adds interpolation/storage choices; apparent positions also require light time and coordinate reductions. Consult the source publication and actual file coverage before assigning an error bound.

The example prints several decimals to illustrate output, not to certify every digit. J2000.0 is noon **TT**, so use `calc` rather than passing that same label to `calc_ut`:

```python
import libephemeris as ephem

# J2000.0 is January 1, 2000, noon TT, not noon UT1.
jd_tt = ephem.julday(2000, 1, 1, 12.0)
pos, _ = ephem.calc(jd_tt, ephem.SUN, ephem.FLG_SPEED)
print(f"Sun at J2000.0:")
print(f"  Longitude: {pos[0]:.9f}°")
print(f"  Latitude:  {pos[1]:.9f}°")
print(f"  Distance:    {pos[2]:.9f} AU")
print(f"  Speed:    {pos[3]:.6f}°/day")
```

```
Sun at J2000.0:
  Longitude: 280.368165285°
  Latitude:  0.000227407°
  Distance:    0.983327632 AU
  Speed:    1.019434°/day
```

### Moon

The Moon is the most complex celestial body to calculate. Its gravitational perturbations (due to the Sun, Jupiter, and the non-spherical shape of the Earth) require thousands of terms. The chosen lunar trajectory and apparent-place reductions both contribute to the result. Nearby lunar geometry also makes observer location and Earth orientation important: geocentric agreement alone does not validate a local eclipse contact.

### Fixed Stars

The compiled catalog combines **Hipparcos**, XHIP and identification sources. Hipparcos observations have epoch J1991.25; stored J2000 coordinates are a propagated representation, not new observations. Each star has its own uncertainty. Proper-motion uncertainty accumulates away from the reference epoch, and binary motion can make simple propagation inadequate. Do not assign a single bound to every star and century.

### Astrological Houses

Cusp calculations are purely geometric (depending on ARMC, latitude, and obliquity of the ecliptic). Reproducibility depends on matching the house definition, UT1, latitude and obliquity. Near polar singularities, small input changes can produce large shifts or make a construction unavailable. A difference between systems is not itself numerical error.

### Ayanamsha

An ayanamsha is a chosen sidereal convention, not a directly measured universal angle. Compare the same defining epoch, precession model and stellar anchor before interpreting a discrepancy. Agreement with another program measures consistency under those choices; it does not establish one physically correct zodiacal origin.

### Identify the physical model

Two useful examples of why source choice matters are:

- **Planetary magnitudes**: models based on Mallama & Hilton (2018), with body- and phase-dependent validity
- **Apparent diameters**: IAU 2015 values, updated from historical constants

### Delta T: Where Precision Drops

**Delta T** (the difference TT − UT) is the Achilles' heel of any historical ephemeris. For dates far from the present, Delta T is estimated using empirical formulas whose uncertainty grows rapidly:

```python
import libephemeris as ephem

# Delta T at different epochs
for anno in [1600, 1800, 1900, 1950, 2000, 2024]:
    jd = ephem.julday(anno, 1, 1, 12.0)
    dt = ephem.deltat(jd)
    print(f"Year {anno}: Delta T = {dt * 86400:7.2f} seconds")
```

```
Year 1600: Delta T =  109.11 seconds
Year 1800: Delta T =   12.60 seconds
Year 1900: Delta T =   -2.70 seconds
Year 1950: Delta T =   29.15 seconds
Year 2000: Delta T =   63.83 seconds
Year 2024: Delta T =   69.18 seconds
```

For the year 2000, Delta T is known with sub-millisecond precision (thanks to IERS data). For the year 1600, the uncertainty is several seconds. For the year -500 (ancient astronomy), the uncertainty can be **minutes** — which translates to significant errors in the Moon's position.

---

## 14.2 Time Scales of Error

Not all calculations have the same temporal validity. Every component of the library has a "precision horizon" beyond which the error grows.

### JPL Ephemerides

JPL ephemerides are generated by integrating the equations of motion of the Solar System. Each version covers a specific range:

- **DE440** (`medium` tier): **1550–2650** — full precision for over 1000 years
- **DE440s** (`base` tier): **1849–2150** — lightweight version, sufficient for modern use
- **DE441** (`extended` tier): **-13200 to +17191** — for historical and prehistoric research

Outside the chosen ephemeris' range, the library raises an error: it does not produce false results.

### Delta T

Before about 1600, the uncertainty of Delta T becomes the dominant error factor. For the Moon (which moves by ~0.5"/second), a 30-second error in Delta T produces an error of ~15" in position — still acceptable for astrology, but not for precise historical eclipses.

### Fixed Stars

The proper motion of stars is measured precisely by the Hipparcos catalog (epoch J2000.0). The propagation is linear and precise for **a few centuries** from the catalog. For very ancient or future dates, the actual proper motion might be non-linear (due to binary companions, gravitational perturbations).

### Minor Bodies (Keplerian)

For asteroids and comets calculated with Keplerian propagation (without SPK files), precision degrades **rapidly** — arcminutes in just a few months. For this reason, the library supports the automatic download of high-precision SPK files from JPL Horizons.

### Rule of Thumb

For a modern chart, start with base data, correct UT1/TT conversion and an explicit house/sidereal convention. For historical or scientific work, state the required observable and tolerated uncertainty. Compare aligned inputs and independent observations where available; a familiar date range does not remove the need to check the actual calculation.

---

## 14.3 Precision Tiers

LibEphemeris organizes ephemerides into three **precision tiers**, each with a different trade-off between temporal coverage and file size:

- **`base`** — uses `de440s.bsp` (~31 MB), covers **1849–2150**. Ideal for modern applications that do not require historical dates.

- **`medium`** (default) — uses `de440.bsp` (~114 MB), covers **1550–2650**. The best trade-off for most uses.

- **`extended`** — uses `de441.bsp` (~3.1 GB), covers **-13200 to +17191**. For historical research, astronomical archaeology, and calculations spanning millennia.

### Selecting a Tier

```python
import libephemeris as ephem

# View current tier
print(f"Active tier: {ephem.get_precision_tier()}")

# Change tier
ephem.set_precision_tier("base")
print(f"New tier: {ephem.get_precision_tier()}")

# Revert to default
ephem.set_precision_tier("medium")
print(f"Restored tier: {ephem.get_precision_tier()}")
```

```
Active tier: base
New tier: base
Restored tier: medium
```

The tier can also be set via the `LIBEPHEMERIS_PRECISION` environment variable:

```bash
export LIBEPHEMERIS_PRECISION=extended
```

### Information on Available Tiers

```python
import libephemeris as ephem
from libephemeris.state import list_tiers

for tier in list_tiers():
    print(f"{tier.name:10s}  {tier.ephemeris_file:14s}  {tier.description}")
```

```
base        de440s.bsp      Modern usage (1850-2150), ~31 MB
medium      de440.bsp       General purpose (1550-2650), ~114 MB
extended    de441.bsp       Extended range (-13200 to +17191), ~3.1 GB
```

---

## 14.4 LEB Mode: Precomputed Binary Ephemerides

### The Problem

The native Rust JPL/SPK route performs several steps for each `calc_ut()` request: reading segments, interpolating Chebyshev polynomials, rotating frames and applying apparent-place corrections. For a single calculation, this is imperceptible, but for thousands (a monthly ephemeris, a transit search over an entire year), the cost accumulates.

### The Solution: LEB

**LEB** (LibEphemeris Binary) is a precomputed ephemeris format that stores approximations with **Chebyshev polynomials** for every celestial body and every time interval. A `.leb` or `.leb2` file stores polynomial channels and auxiliary data for efficient local access. Raw channels still require appropriate apparent-place reductions; LEB does not skip all physics after one polynomial evaluation. Chunked LEB2 can reduce I/O, while indexed reading bounds decoded coefficient storage.

Compare LEB with the native JPL route using the same epoch, frame and settings. This example needs a local JPL kernel as well as LEB data. It measures one longitude difference, not the worst error over every body and date:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Native Rust JPL/SPK (historical mode name)
ephem.set_calc_mode("skyfield")
pos_sky, _ = ephem.calc_ut(jd, ephem.MARS, ephem.FLG_SPEED)

# Calculation with LEB
ephem.set_calc_mode("leb")
pos_leb, _ = ephem.calc_ut(jd, ephem.MARS, ephem.FLG_SPEED)

diff_arcsec = abs(pos_sky[0] - pos_leb[0]) * 3600
print(f"Mars (native JPL): {pos_sky[0]:.9f}°")
print(f"Mars (LEB):      {pos_leb[0]:.9f}°")
print(f"Difference:       {diff_arcsec:.9f}\"")

ephem.set_calc_mode("auto")
```

```
Mars (native JPL): 342.845795941°
Mars (LEB):      342.845795941°
Difference:       0.000001281"
```

### Activating LEB

There are three ways to activate LEB mode:

**1. Download the precomputed LEB file (recommended):**

```python
import libephemeris as ephem

# Download the LEB file for the medium tier (~175 MB)
ephem.download_leb_for_tier("medium")
# The file is saved in ~/.libephemeris/leb/ephemeris_medium.leb
# and automatically activated for this session
```

**2. Set the path manually:**

```python
import libephemeris as ephem

import os
from pathlib import Path

# Set LIBEPHEMERIS_DATA_DIR to a provisioned local root before starting Python.
path = Path(os.environ["LIBEPHEMERIS_DATA_DIR"]) / "leb" / "base_core.leb2"
assert path.is_file(), "Provision the base LEB2 group first"
ephem.set_leb_file(path)
position, _ = ephem.calc_ut(2460409.0, ephem.SUN, ephem.FLG_SPEED)
assert len(position) == 6
```

**3. Via environment variable:**

```bash
export LIBEPHEMERIS_LEB=/path/to/file/ephemeris_medium.leb
```

### Auto-discovery

If you have downloaded a LEB file with `download_leb_for_tier()`, the library automatically finds it in the standard location (`~/.libephemeris/leb/ephemeris_{tier}.leb`) without needing any configuration.

### Automatic Fallback

Native LEB supports topocentric, Cartesian, radian and non-nutated output where the body/model supplies the required inputs. These flags no longer imply a switch to Python Skyfield.

Forced `leb` stays on local coefficients and supported local models. Missing bodies, corrupt explicit files and uncovered epochs are errors, not permission to silently use JPL or Horizons. In `auto`, provider routing is a separate policy. A service may prefer an explicit failure over an unnoticed data change or network request. Inspect coverage and provision missing groups before computation.

### Three LEB Tiers

Like the JPL ephemerides, LEB files also exist in three versions:

- **`base`** — modern coverage, nominally 1850–2150 for the core.
- **`medium`** — a wider historical interval, nominally 1550–2650.
- **`extended`** — long-range products; inspect actual body and auxiliary bounds rather than assuming the JPL kernel's full span.

LEB1 and chunked LEB2 have different sizes, and optional asteroid/apogee/exotic
groups are separate products. Use `libephemeris leb info FILE --json` to see
which channels and epochs your actual image contains. A larger tier is not a
promise that every additional body is present or intrinsically more accurate.

---

## 14.5 Calculation Modes and Configuration

The library supports four calculation modes, which can be controlled with `set_calc_mode()`:

### `"auto"` (default)

Uses configured resources and native provider routing. Availability, body support and network policy govern the result; do not assume a universal LEB → Horizons → Skyfield order. For a reproducible offline exercise, explicit `leb` plus `sealed` is easier to reason about:

```python
import libephemeris as ephem

ephem.set_calc_mode("auto")
print(f"Mode: {ephem.get_calc_mode()}")
```

```
Mode: auto
```

### `"skyfield"`

Selects native Rust JPL/SPK, even when LEB is available. The string is a compatibility name, not a Python dependency. Prepare a suitable local JPL kernel first; keep the network sealed for local comparisons:

```python
import libephemeris as ephem

ephem.set_calc_mode("skyfield")
print(f"Mode: {ephem.get_calc_mode()}")

jd = ephem.julday(2024, 4, 8, 12.0)
pos, _ = ephem.calc_ut(jd, ephem.SUN, ephem.FLG_SPEED)
print(f"Sun (native JPL): {pos[0]:.6f}°")

ephem.set_calc_mode("auto")  # restore
```

```
Mode: skyfield
Sun (native JPL): 19.140437°
```

### `"leb"`

Selects local LEB coefficients and supported local models. Prepare the required image explicitly; the bundled base is sufficient for the first exercises. Missing/corrupt explicit files and uncovered requests raise errors. There is no implicit JPL/Horizons substitution, and sealed computation cannot download missing inputs:

```python
import libephemeris as ephem

ephem.set_calc_mode("leb")
print(f"Mode: {ephem.get_calc_mode()}")

ephem.set_calc_mode("auto")  # restore
```

```
Mode: leb
```

### `"horizons"`

Selects supported JPL Horizons network operations. This is an optional online example, not part of the offline exercises. Check unsupported combinations against the current function contract rather than assuming every body and flag has a remote equivalent. Set network policy deliberately and restore it afterwards:

```python
import libephemeris as ephem

ephem.set_network_policy("allow")  # optional online example
ephem.set_calc_mode("horizons")
print(f"Mode: {ephem.get_calc_mode()}")

jd = ephem.julday(2024, 4, 8, 12.0)
pos, _ = ephem.calc_ut(jd, ephem.SUN, ephem.FLG_SPEED)
print(f"Sun (Horizons): {pos[0]:.6f}°")

ephem.set_calc_mode("leb")
ephem.set_network_policy("sealed")  # restore the offline exercises
```

```
Mode: horizons
Sun (Horizons): 19.140437°
```

### Environment Variable

The mode can also be set via `LIBEPHEMERIS_MODE`:

```bash
export LIBEPHEMERIS_MODE=horizons
```

---

## 14.6 EphemerisContext: Thread-Safe Calculations

For multi-threaded applications (web servers, APIs, parallel calculations), the library offers `EphemerisContext` — a context that maintains its own isolated state (observer position, sidereal mode, angle cache) while sharing expensive resources (ephemeris files, timescale) in a thread-safe manner.

### Basic Use

```python
import libephemeris as ephem
from libephemeris import EphemerisContext
from libephemeris.constants import SUN, FLG_SPEED

ctx = EphemerisContext()
ctx.set_topo(12.5, 41.9, 0)  # Rome

jd = ephem.julday(2024, 4, 8, 12.0)
pos, flag = ctx.calc_ut(jd, SUN, FLG_SPEED)
print(f"Sun (from Rome, via context): {pos[0]:.4f}°")
```

```
Sun (from Rome, via context): 19.1404°
```

### Sidereal Calculation in Context

Each context can have its own sidereal mode:

```python
import libephemeris as ephem
from libephemeris import EphemerisContext
from libephemeris.constants import SUN, FLG_SPEED, FLG_SIDEREAL
from libephemeris.constants import SIDM_LAHIRI

ctx = EphemerisContext()
ctx.set_sid_mode(SIDM_LAHIRI)

jd = ephem.julday(2024, 4, 8, 12.0)
pos, _ = ctx.calc_ut(jd, SUN, FLG_SPEED | FLG_SIDEREAL)
print(f"Sidereal Sun (Lahiri): {pos[0]:.4f}°")
```

```
Sidereal Sun (Lahiri): 354.9458°
```

### Houses in Context

```python
import libephemeris as ephem
from libephemeris import EphemerisContext

ctx = EphemerisContext()
jd = ephem.julday(2024, 4, 8, 12.0)

cusps, ascmc = ctx.houses(jd, 41.9, 12.5, ord('P'))
print(f"Ascendant: {ascmc[0]:.4f}°")
print(f"Midheaven: {ascmc[1]:.4f}°")
```

```
Ascendant: 133.0806°
Midheaven: 31.9078°
```

### Parallel Calculations with Threads

```python
import libephemeris as ephem
from libephemeris import EphemerisContext
from libephemeris.constants import SUN, FLG_SPEED
import threading

results = {}

def calc_for_city(name, lon, lat):
    ctx = EphemerisContext()
    ctx.set_topo(lon, lat, 0)
    jd = ephem.julday(2024, 4, 8, 12.0)
    pos, _ = ctx.calc_ut(jd, SUN, FLG_SPEED)
    results[name] = pos[0]

cities = [
    ("Rome", 12.5, 41.9),
    ("New York", -74.0, 40.7),
    ("Tokyo", 139.7, 35.7),
]

threads = []
for name, lon, lat in cities:
    t = threading.Thread(target=calc_for_city, args=(name, lon, lat))
    threads.append(t)
    t.start()

for t in threads:
    t.join()

# Check every worker produced a result; completion order is not fixed.
assert len(results) == len(cities)
for name, lon in sorted(results.items()):
    print(f"{name}: Sun at {lon:.4f}°")
```

```
New York: Sun at 19.1404°
Rome: Sun at 19.1404°
Tokyo: Sun at 19.1404°
```

> **Note**: Geocentric positions are the same for all cities because the `FLG_TOPOCTR` flag was not used. The difference can be seen in the houses and topocentric positions.

### LEB in Context

Each context can have its own LEB file:

```python
from libephemeris import EphemerisContext

ctx = EphemerisContext()
import os
from pathlib import Path

path = Path(os.environ["LIBEPHEMERIS_DATA_DIR"]) / "leb" / "base_core.leb2"
assert path.is_file(), "Provision the base LEB2 group first"
ctx.set_leb_file(path)
# The mode was selected in the setup; choosing a path alone does not change it.
position, _ = ctx.calc_ut(2460409.0, ephem.SUN, ephem.FLG_SPEED)
assert len(position) == 6
```

### Closing Shared Resources

```python
from libephemeris import EphemerisContext

# Closes files and resources shared by all contexts
EphemerisContext.close()
```

---

## 14.7 Download and Data Management

Prepare required data explicitly before calculation. Acquisition depends on mode and network policy: the sealed exercises must not download missing files. The following provisioning examples are optional network operations, separate from ordinary computation.

### Downloading Data for a Tier

```python
import libephemeris as ephem

# Download everything needed for the "medium" tier:
# - de440.bsp (JPL ephemerides)
# - planet_centers_medium.bsp (planet centers)
# - SPK for all minor bodies
ephem.download_for_tier("medium")
```

### Downloading the LEB File

```python
import libephemeris as ephem

# Download the precomputed LEB ephemerides for the "medium" tier
# (~175 MB, automatically activated after download)
ephem.download_leb_for_tier("medium")
```

### Data Directory

By default, all files are saved in `~/.libephemeris/`. You can change this directory with the environment variable:

```bash
export LIBEPHEMERIS_DATA_DIR=/data/ephemerides
```

Or verify the current directory:

```python
import libephemeris as ephem

from pathlib import Path

data_directory = ephem.get_library_path()
print(f"Data directory is absolute: {Path(data_directory).is_absolute()}")
```

```
Data directory is absolute: True
```

### Complete Reset

`close()` releases default shared resources and resets transient observing state and some selections. It does not clear every override or invalidate retained handles. It is a lifecycle operation, not cleanup after each request. To compare before/after values, reapply the same configuration before calculating again:

```python
import libephemeris as ephem

# Do some calculations...
jd = ephem.julday(2024, 4, 8, 12.0)
pos, _ = ephem.calc_ut(jd, ephem.SUN, 0)
print(f"Before close: {pos[0]:.6f}°")

# Close everything
ephem.close()

# The next calculation automatically reloads the ephemerides
ephem.set_precision_tier("base")
ephem.set_iers_delta_t_enabled(False)
pos2, _ = ephem.calc_ut(jd, ephem.SUN, 0)
print(f"After close: {pos2[0]:.6f}°")
print(f"Identical result: {abs(pos[0] - pos2[0]) < 0.000001}")
```

```
Before close: 19.140437°
After close: 19.140437°
Identical result: True
```

### Environment Variables Summary

- `LIBEPHEMERIS_DATA_DIR` — base directory for all data (default: `~/.libephemeris`)
- `LIBEPHEMERIS_PRECISION` — precision tier: `base`, `medium`, `extended`
- `LIBEPHEMERIS_EPHEMERIS` — ephemeris file name (e.g., `de441.bsp`)
- `LIBEPHEMERIS_LEB` — path to the `.leb` file for binary mode
- `LIBEPHEMERIS_MODE` — calculation mode: `auto`, `skyfield`, `leb`, `horizons`
- `LIBEPHEMERIS_AUTO_SPK` — `1`/`0` to enable/disable automatic SPK download
- `LIBEPHEMERIS_SPK_DIR` — cache directory for minor bodies' SPK files

---

## Summary

- **Precision depends on the observable**, selected data, epoch, physical model and observer; printed decimals are not an accuracy guarantee.
- **Three precision tiers**: `base` (1849–2150, 31 MB), `medium` (1550–2650, 114 MB), `extended` (-13200 to +17191, 3.1 GB)
- **LEB**: local coefficients with their own coverage and interpolation budget; inspect the selected file rather than assuming identical results across routes.
- **Four calculation modes**: `auto` (configured routing), `skyfield` (native Rust JPL/SPK), `leb` (local coefficients/models), `horizons` (explicit network route).
- **EphemerisContext**: thread-safe context for parallel calculations, with isolated state for observer, sidereal mode, and cache
- **Explicit provisioning** with `download_for_tier()` and `download_leb_for_tier()`; sealed calculations cannot acquire missing data.
- **`close()`** releases default shared resources and resets selected state, not every configuration override.

## Functions Introduced

- `set_precision_tier(tier)` / `get_precision_tier()` — selects the precision tier (`"base"`, `"medium"`, `"extended"`)
- `set_calc_mode(mode)` / `get_calc_mode()` — sets the calculation mode (`"auto"`, `"skyfield"`, `"leb"`, `"horizons"`)
- `set_leb_file(filepath)` — selects a coefficient image; select the calculation mode separately
- `download_for_tier(tier)` — downloads all data for a tier
- `download_leb_for_tier(tier)` — downloads the precomputed LEB file
- `EphemerisContext()` — isolated context for thread-safe calculations
- `EphemerisContext.calc_ut(jd, body, flag)` — position calculation in context
- `EphemerisContext.houses(jd, lat, lon, hsys)` — houses in context
- `EphemerisContext.set_topo(lon, lat, alt)` — observer in context
- `EphemerisContext.set_sid_mode(mode)` — sidereal mode in context
- `EphemerisContext.set_leb_file(filepath)` — LEB for context
- `EphemerisContext.close()` — closes shared resources
- `close()` — releases shared resources and resets selected state; retained handles and explicit overrides can survive
- `get_library_path()` — path of the data directory