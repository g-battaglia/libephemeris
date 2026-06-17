# LibEphemeris

<div align="left">
    <img src="https://static.pepy.tech/badge/libephemeris/month" alt="PyPI Downloads">
    <img src="https://static.pepy.tech/badge/libephemeris/week" alt="PyPI Downloads">
    <img src="https://static.pepy.tech/personalized-badge/libephemeris?period=total&units=INTERNATIONAL_SYSTEM&left_color=GREY&right_color=BLUE&left_text=downloads/total" alt="PyPI Downloads">
    <img src="https://img.shields.io/pypi/v/libephemeris.svg" alt="PyPI Version">
    <img src="https://img.shields.io/pypi/pyversions/libephemeris.svg" alt="Python Versions">
    <img src="https://img.shields.io/github/license/g-battaglia/libephemeris.svg" alt="License">
</div>

A high-precision astronomical ephemeris library for Python, powered by NASA JPL DE440/DE441 ephemerides and IAU 2006/2000A standards.

**Drop-in replacement for PySwissEph** - readable Python algorithms, standard debugging, easy deployment on the scientific Python stack (NumPy, Skyfield, pyerfa).

---

## Features

- **NASA JPL DE440/DE441** - modern planetary ephemerides via Skyfield, with full-range DE441 support for deep-history and far-future work
- **IAU + Vondrák 2011 standards** - long-term precession and of-date mean obliquity (Vondrák 2011, valid ±200,000 years), nutation (IAU 2006/2000A) via the official ERFA library
- **Validated high precision** - planetary differences typically measured in fractions of an arcsecond, house cusps < 0.02", benchmarked across 4,400+ comparison rounds ([full report](https://github.com/g-battaglia/libephemeris/blob/main/docs/PRECISION.md))
- **Four backends, one API** - Skyfield, LEB (~14x speedup), Horizons API, and adaptive auto mode through the same `calc_ut()` interface
- **25 house systems, 43 ayanamsha modes** - independently verified against pyswisseph
- **Physical planet centers** - outer planets corrected from barycenters using JPL satellite ephemerides
- **Thread-safe contexts when you need them** - SwissEph-compatible globals for drop-in migration, `EphemerisContext` for concurrent workloads
- **15,000+ years of coverage** - `base`, `medium`, and `extended` precision tiers from modern use to -13200 / +17191 CE
- **Readable Python 3.12+** - the ephemeris algorithms are plain, inspectable Python; clean installs across CI, containers, and serverless from prebuilt scientific wheels

---

## Why LibEphemeris

Swiss Ephemeris is the industry standard for planetary calculations. But its Python binding (pyswisseph) wraps a large opaque C library - hard to build from source, hard to inspect or debug, tied to a single computation model.

LibEphemeris provides the **same API** with a modern foundation:

- **NASA JPL ephemerides** instead of semi-analytical theory - DE440/DE441 are the latest planetary ephemerides from the Jet Propulsion Laboratory, the same data used for spacecraft navigation.
- **IAU + Vondrák 2011 standards** - long-term precession and of-date mean obliquity (Vondrák, Capitaine & Wallace 2011, valid ±200,000 years instead of the IAU 2006 polynomial's few centuries), nutation (IAU 2006/2000A), all computed via the official ERFA library (the open-source implementation of IAU SOFA), not custom routines.
- **Physical planet centers** - Jupiter, Saturn, Uranus, Neptune corrected from system barycenters to actual body centers using JPL satellite ephemerides. Most libraries skip this.
- **Independently verified** - every function cross-validated against pyswisseph, JPL Horizons, and astropy/ERFA. [Precision report with full methodology](https://github.com/g-battaglia/libephemeris/blob/main/docs/PRECISION.md).
- **Readable Python algorithms** - plain, inspectable source and standard debugging instead of an opaque C library. Installs from prebuilt wheels (NumPy/Skyfield/pyerfa) across any platform, CI, or serverless environment.

**Switching from pyswisseph?** Your existing code works with minimal changes. [Migration guide](https://github.com/g-battaglia/libephemeris/blob/main/docs/guides/migration-guide.md).

### Where LibEphemeris is measurably *more* accurate

Matching the industry standard is the baseline; in two places the modern foundation does demonstrably better, verified against an independent physical criterion (not against any other engine):

- **House cusps over deep time.** House cusps are driven by sidereal time and obliquity, i.e. by precession. Engines that take sidereal time from an IAU 1976/2006 polynomial are only correct for a few centuries and drift by **degrees** at historical or far-future dates. LibEphemeris derives the cusps from the long-term Vondrák 2011 model (valid ±200,000 years) via a stable geometric construction, so cusps stay correct across the whole ±13,000-year ephemeris range — and, because houses and bodies share one obliquity and one ΔT, a chart is internally self-consistent at any epoch. Identical to the reference in the modern era; correct where the reference's own model diverges.

- **House cusp speeds (daily motion).** A cusp's speed is its true time derivative dλ/dt. LibEphemeris computes it as the genuine derivative of the full house solution, so the reported speed integrates to the cusp's actual motion. For the iteratively-solved systems (Placidus, Koch) a common analytic speed approximation mis-states the intermediate cusps by tens to hundreds of °/day near the polar circle; here LibEphemeris matches the real motion to **< 0.005 °/day**.

Full methodology: [Long-term sidereal time, precession & cusp speeds](https://github.com/g-battaglia/libephemeris/blob/main/docs/methodology/sidereal-time-longterm.md).

---

## Quick Start

```python
import libephemeris as swe
from libephemeris.constants import SUN, MOON, FLG_SPEED

jd = swe.julday(2000, 1, 1, 12.0)  # J2000.0

sun, _ = swe.calc_ut(jd, SUN, FLG_SPEED)
moon, _ = swe.calc_ut(jd, MOON, FLG_SPEED)

print(f"Sun:  {sun[0]:.4f} deg, speed {sun[3]:.4f} deg/day")
print(f"Moon: {moon[0]:.4f} deg, speed {moon[3]:.4f} deg/day")
```

```python
# House cusps (Placidus, Rome)
cusps, ascmc = swe.houses(swe.julday(2024, 11, 5, 18.0), 41.9028, 12.4964, b"P")
print(f"ASC: {ascmc[0]:.4f}, MC: {ascmc[1]:.4f}")
```

For concurrent or multi-threaded workloads, use `EphemerisContext` instead of the module-level global state.

More examples: [Getting Started](https://github.com/g-battaglia/libephemeris/blob/main/docs/guides/getting-started.md)

---

## Verified Precision

Every number independently measured against pyswisseph across 4,400+ comparison rounds. [Full report](https://github.com/g-battaglia/libephemeris/blob/main/docs/PRECISION.md).

| Category | Typical | Max | Scope |
|----------|---------|-----|-------|
| Sun-Pluto | 0.04-0.26" | 1.17" | Pluto max 0.75"; Neptune is the widest observed planetary delta |
| Moon | 0.70" | 3.32" | Different underlying lunar models |
| House cusps | < 0.01" | 0.02" | All 25 systems (modern); long-term model holds to ±13,000 yr |
| House cusp speeds | < 0.005°/day | — | True dλ/dt; more accurate than analytic approximations on Placidus/Koch |
| Fixed stars | < 0.1" | 0.51" | 116 Hipparcos catalog stars |
| Solar eclipses | - | < 6s | Timing accuracy |
| Lunar eclipses | - | < 8s | Timing accuracy |
| Ayanamsha | < 0.001 deg | 0.006 deg | All 43 sidereal modes |

---

## Four Backends, One API

Choose your trade-off between speed, locality, and setup. The same `calc_ut()` interface works across all four modes, from zero-install Horizons lookups to precomputed LEB throughput.

| Mode | Backend | Speed | Use case |
|------|---------|-------|----------|
| `"auto"` | LEB -> Horizons -> Skyfield | adaptive | **Default.** Best onboarding; resolves local or remote data transparently |
| `"skyfield"` | JPL DE440/DE441 via Skyfield | ~120 us | High-precision local JPL workflow |
| `"leb"` | Precomputed Chebyshev polynomials | ~5 us | Maximum throughput for repeated calculations |
| `"horizons"` | NASA JPL Horizons REST API | ~300 ms | No local ephemeris files required |

```python
from libephemeris import set_calc_mode
set_calc_mode("leb")  # or via env: LIBEPHEMERIS_MODE=leb
```

---

## Installation

```bash
pip install libephemeris
```

Out of the box, the wheel includes a bundled LEB2 base-tier core for the 14 core bodies (1850-2150). With the default `medium` tier, the library can auto-download the additional LEB2 data it needs on first use.

Recommended first-time setup:

```bash
libephemeris init                 # Optional but recommended interactive config
libephemeris download auto        # Download exactly what your config needs
libephemeris status               # Verify installed data and active setup
```

Prefer to install a tier directly? Use one of these:

```bash
libephemeris download base         # 1850-2150, lightweight
libephemeris download medium       # 1550-2650, ~200 MB (recommended)
libephemeris download extended     # -13200 to +17191 CE, full range
```

**Optional extras:** `pip install libephemeris[stars]` for star-catalog tooling, `[nbody]` for REBOUND/ASSIST n-body integration, `[all]` for everything. [Details](https://github.com/g-battaglia/libephemeris/blob/main/docs/guides/getting-started.md#optional-extras).

---

## Documentation

- [Getting Started](https://github.com/g-battaglia/libephemeris/blob/main/docs/guides/getting-started.md) - installation, ephemeris tiers, first calculations
- [Migration from PySwissEph](https://github.com/g-battaglia/libephemeris/blob/main/docs/guides/migration-guide.md) - API mapping, flag compatibility, known divergences
- [Precision Report](https://github.com/g-battaglia/libephemeris/blob/main/docs/PRECISION.md) - full methodology, comparison tables, verification process
- [Long-term sidereal time, precession & cusp speeds](https://github.com/g-battaglia/libephemeris/blob/main/docs/methodology/sidereal-time-longterm.md) - why houses and cusp speeds stay correct over ±13,000 years
- [Flag Reference](https://github.com/g-battaglia/libephemeris/blob/main/docs/reference/flags.md) - all supported flags with examples
- [House Systems](https://github.com/g-battaglia/libephemeris/blob/main/docs/reference/house-systems.md) - all 25 systems, verified against pyswisseph
- [Ayanamsha Modes](https://github.com/g-battaglia/libephemeris/blob/main/docs/reference/ayanamsha.md) - 43 sidereal modes
- [LEB Binary Ephemeris](https://github.com/g-battaglia/libephemeris/blob/main/docs/leb/guide.md) - format, generation, LEB2 compression
- [Horizons Backend](https://github.com/g-battaglia/libephemeris/blob/main/docs/architecture/horizons-backend.md) - HTTP client, pipeline, precision
- [Architecture](https://github.com/g-battaglia/libephemeris/blob/main/docs/development/architecture-overview.md) - internal design and data flow
- [Methodology](https://github.com/g-battaglia/libephemeris/blob/main/docs/methodology/overview.md) - planet centers, lunar apsides, pyerfa integration
- [CLI Reference](https://github.com/g-battaglia/libephemeris/blob/main/CLI.md) - full command reference
- [Changelog](https://github.com/g-battaglia/libephemeris/blob/main/CHANGELOG.md) - release history

---

## Contributing

```bash
git clone https://github.com/g-battaglia/libephemeris.git
cd libephemeris && uv pip install -e ".[dev]"
poe lint                           # Ruff lint + auto-fix
poe test:leb:fast                  # Recommended fast unit suite
poe test:compare:skyfield          # Cross-validate vs pyswisseph
```

---

## Part of the Kerykeion Ecosystem

LibEphemeris is the computation engine behind:

- **[Astrologer Studio](https://www.astrologerstudio.com)** — professional online astrology software (in production)
- **[Kerykeion](https://github.com/g-battaglia/kerykeion)** — Python astrology library (v6 alpha)
- **[Astrologer API](https://www.kerykeion.net/astrologer-api)** — hosted REST API for astrology data and SVG charts (upcoming)

Learn more at [kerykeion.net](https://kerykeion.net).

---

## License

Dual-licensed:

- **[AGPL-3.0-only](LICENSE)** — free for any use that complies with the
  AGPL (including network/SaaS source-disclosure obligations). PyPI
  releases carry this license.
- **Commercial license** — for closed-source products and services; see
  [LICENSING.md](LICENSING.md).

LibEphemeris is an independent, API-compatible implementation; it contains
no Swiss Ephemeris code ([NOTICE.md](NOTICE.md),
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)). "Swiss Ephemeris" is a
product of Astrodienst AG.
