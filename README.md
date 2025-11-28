# LibEphemeris

**High-precision astronomical ephemeris library for Python (Swiss Ephemeris compatible, powered by Skyfield and JPL DE ephemerides).**

> [!WARNING]
> **Pre-Alpha Release**
>
> LibEphemeris is currently in an early pre-alpha stage. The public API is unstable and may change without notice. Do not rely on it in production yet.

LibEphemeris is an open-source alternative to Swiss Ephemeris that provides scientific-grade astronomical calculations using NASA JPL ephemerides via [Skyfield](https://rhodesmill.org/skyfield/). The goal is to mirror the Swiss Ephemeris API (as exposed by `pyswisseph`) while using modern, freely available data and a pure-Python implementation.

[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org)

---

## Features at a Glance

- **Planetary positions**: Sun, Moon, all major planets and Pluto.
- **High precision**: Based on NASA JPL DE421 by default (configurable to other DE files).
- **Multiple coordinate systems**: Ecliptic, equatorial, J2000 and of-date frames.
- **Observation modes**: Geocentric, topocentric, heliocentric, barycentric.
- **Velocities**: Full 6-component state vectors (position + velocity).
- **House systems (19)**: Placidus, Koch, Regiomontanus, Campanus, Equal, Whole Sign, Porphyry, Alcabitius, Polich/Page (Topocentric), Morinus, Meridian, Vehlow, Horizontal, Carter, Krusinski, Natural Gradient, and more.
- **Sidereal zodiac (43 ayanamshas)**: Fagan/Bradley, Lahiri, Raman, Krishnamurti, star-based and historical variants.
- **Extended points**: Lunar nodes, Lilith (mean and true), major asteroids (Chiron, Pholus, Ceres, Pallas, Juno, Vesta), TNOs (Orcus, Haumea, Quaoar, Makemake, Gonggong, Eris, Sedna), major fixed stars and Arabic parts.
- **Event finding**: Sun/Moon longitude crossings (e.g. ingress), with additional events planned (eclipses, etc.).
- **Swiss Ephemeris compatibility**: Same function names, flags and result structure as `pyswisseph` in most common use cases.

---

## Installation

Using `pip`:

```bash
pip install libephemeris
```

Using [`uv`](https://github.com/astral-sh/uv) (recommended for development):

```bash
uv pip install libephemeris
```

From source:

```bash
git clone https://github.com/g-battaglia/libephemeris.git
cd libephemeris
uv pip install -e .
```

### Requirements

- Python **3.10+**
- `skyfield>=1.53`
- `skyfield-data>=7.0.0`
- A JPL ephemeris file (DE421 by default, downloaded automatically on first use if not present locally)

---

## Quick Start

### Basic planetary positions

```python
import libephemeris as ephem
from libephemeris.constants import *

# Julian Day (J2000.0)
jd = ephem.swe_julday(2000, 1, 1, 12.0)

# Sun position (longitude, latitude, distance, and speeds)
sun, flags = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_SWIEPH | SEFLG_SPEED)
print(f"Sun longitude: {sun[0]:.6f}°")
print(f"Sun speed: {sun[3]:.6f}°/day")

# Moon position
moon, _ = ephem.swe_calc_ut(jd, SE_MOON, SEFLG_SWIEPH)
print(f"Moon longitude: {moon[0]:.6f}°")
```

### Houses and angles

```python
# Rome coordinates
lat, lon, alt = 41.9028, 12.4964, 0.0
jd = ephem.swe_julday(2024, 11, 5, 18.0)

# Placidus houses
cusps, ascmc = ephem.swe_houses(jd, lat, lon, b"P")

print(f"Ascendant: {ascmc[0]:.2f}°")
print(f"MC:        {ascmc[1]:.2f}°")
print(f"House 1:   {cusps[1]:.2f}°")
print(f"House 10:  {cusps[10]:.2f}°")
```

### Sidereal calculations

```python
# Lahiri ayanamsha
ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

ayanamsha = ephem.swe_get_ayanamsa_ut(jd)
print(f"Lahiri ayanamsha: {ayanamsha:.6f}°")

sun_sid, _ = ephem.swe_calc_ut(jd, SE_SUN, SEFLG_SWIEPH | SEFLG_SIDEREAL)
print(f"Sidereal Sun: {sun_sid[0]:.6f}°")
```

### Longitude crossings

```python
# Next time the Sun enters 0° Aries
next_cross = ephem.swe_solcross_ut(0.0, jd, SEFLG_SWIEPH)
print(f"Next Aries ingress (JD): {next_cross:.6f}")
```

---

## Configuring Ephemeris Files

By default, LibEphemeris uses the **JPL DE421** kernel (`de421.bsp`), which covers roughly **1900–2050**. The file is:

- loaded from a local path if already present, or
- automatically downloaded via Skyfield the first time it is needed.

You can control which ephemeris file is used and where it is loaded from.

### Choosing a different JPL kernel

Use `set_ephemeris_file()` to select a different `.bsp` file:

```python
from libephemeris import set_ephemeris_file

set_ephemeris_file("de431.bsp")  # very long time span, larger file
```

Supported JPL kernels include, for example:

- `de421.bsp`: 1900–2050 (default, ~16 MB)
- `de422.bsp`: −3000–3000 (~623 MB)
- `de430.bsp`: 1550–2650 (~128 MB)
- `de431.bsp`: −13200–17191 (~3.4 GB)

If the chosen file is not present locally, Skyfield will attempt to download it.

### Custom ephemeris directory

Use `set_ephe_path()` to point LibEphemeris to a directory containing JPL kernels:

```python
from libephemeris import set_ephe_path

set_ephe_path("/path/to/jpl-kernels")
```

Resolution order for the ephemeris file is:

1. The directory set via `set_ephe_path()`, if any.
2. The project/workspace root.
3. Download via Skyfield.

If you try to compute positions outside the date range covered by the selected kernel, LibEphemeris will raise an exception describing the supported range.

---

## Scientific Accuracy and Validation

### Ephemeris data

- **Source**: NASA JPL DE ephemerides (DE421 by default).
- **Time span**: 1900–2050 for DE421; extendable by selecting other kernels.
- **Precision**: Sub-arcsecond accuracy for major planets within the supported range.
- **Reference frame**: ICRS/J2000.0.

### Comparison with Swiss Ephemeris

LibEphemeris has been tested against Swiss Ephemeris using an automated test suite.

| Component               | Tests | Pass Rate | Max Difference |
| ----------------------- | ----- | --------- | -------------- |
| Planetary positions     | 229   | 100%      | < 0.001°       |
| House systems           | 113   | 100%      | < 0.001°       |
| Ayanamsha values        | 129   | 100%      | < 0.06°        |
| Lunar nodes / Lilith    | 40+   | 100%      | < 0.01°        |
| Velocities              | 100   | 100%      | < 0.01°/day    |

These comparisons are implemented in the `tests/` and `compare_scripts/` directories, and are run regularly during development.

---

## Swiss Ephemeris Compatibility

LibEphemeris aims to behave as a **drop-in replacement** for `pyswisseph` in many scenarios:

- Same function names (e.g. `swe_calc_ut`, `swe_houses`, `swe_julday`, `swe_revjul`, `swe_get_ayanamsa_ut`).
- Same integer constants and flags from `libephemeris.constants` (e.g. `SE_SUN`, `SEFLG_SWIEPH`, `SEFLG_SPEED`, `SE_SIDM_LAHIRI`).
- Similar return types and value ordering.

There are still differences and missing features compared to the full Swiss Ephemeris API, especially around:

- very long time ranges (beyond the chosen JPL kernel),
- eclipse and occultation functions,
- the full minor-planet and fixed-star catalogues.

Please open an issue if you hit a compatibility gap that is important for your use case.

---

## Development

### Project layout

```text
libephemeris/
├── libephemeris/
│   ├── __init__.py
│   ├── constants.py      # Constants and flags
│   ├── planets.py        # Planetary calculations
│   ├── houses.py         # House systems
│   ├── lunar.py          # Nodes and Lilith
│   ├── minor_bodies.py   # Asteroids and TNOs
│   ├── fixed_stars.py    # Fixed stars and points
│   ├── crossing.py       # Longitude crossing events
│   ├── angles.py         # Angle helpers (Asc, MC, etc.)
│   ├── arabic_parts.py   # Arabic parts calculations
│   ├── time_utils.py     # Time conversion helpers
│   └── state.py          # Global state (loader, ephemeris, sidereal mode)
├── tests/                # Comprehensive test suite
├── compare_scripts/      # Swiss Ephemeris comparison tools
└── README.md
```

### Development workflow

Install in editable mode with development dependencies:

```bash
uv pip install -e ".[dev]"
```

Run tests (via `poethepoet` tasks defined in `pyproject.toml`):

```bash
poe test       # run pytest
poe coverage   # run tests with coverage
poe lint       # run Ruff (lint)
poe format     # run Ruff formatter
```

---

## License

LibEphemeris is licensed under the **GNU Lesser General Public License v3.0 (LGPL-3.0)**.

See `LICENSE` for the full text.

---

Built for the astronomical and astrological communities.
