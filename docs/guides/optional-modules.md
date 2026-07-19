# Optional Modules and Calculation Backends

LibEphemeris works out of the box for planets and luminaries using Skyfield
with NASA JPL DE440/441 ephemerides. For minor bodies (asteroids, TNOs,
centaurs), multiple calculation backends are available with different
precision levels. This guide explains what each optional module provides,
when to install it, and how the library selects the best available backend
at runtime.

## Table of Contents

- [Calculation Backends](#calculation-backends)
- [Minor Body Fallback Chain](#minor-body-fallback-chain)
- [Optional Extras](#optional-extras)
- [Data Requirements](#data-requirements)
- [Known Limitations](#known-limitations)
- [Choosing the Right Setup](#choosing-the-right-setup)

---

## Calculation Backends

Planets (Sun through Pluto), the Moon, and lunar nodes are served by the
active calculation mode (`auto`, `skyfield`, `leb`, or `horizons`). No
optional modules are needed for these bodies.

Minor bodies use a **fallback chain** of backends, each with different
precision and requirements:

| Backend | Precision | Requirements | Bodies |
|---------|-----------|--------------|--------|
| **SPK kernel** | Sub-arcsecond | Downloaded SPK file per body | All 37 mapped asteroids/TNOs |
| **ASSIST n-body** | Numerical model; validated accuracy is body/date dependent | `libephemeris[nbody]` + source files | Bodies with a supported orbit and source interval |
| **Keplerian** | ~1-10 arcminutes | None (built-in) | Any body with orbital elements |

For planets, the LEB (LibEphemeris Binary) and Horizons API backends are
also available. See [Getting Started](getting-started.md) for ephemeris
tier selection and [the LEB guide](../leb/guide.md) for binary ephemeris
details.

---

## Minor Body Fallback Chain

When you call `calc_ut()` for a minor body, the library tries each backend
in order and uses the first one that succeeds:

```
1. SPK kernel (local file)
     |
     v  not found?
2. Auto-download SPK from JPL Horizons (if enabled)
     |
     v  download failed or body blocked?
3. Strict precision check
     |  If strict_precision=True and body has a downloadable SPK:
     |  raise SPKRequiredError (prevents silent precision loss)
     |
     v  body not strictly required, or strict_precision=False?
4. ASSIST n-body integration (if installed + data available)
     |
     v  not installed?
5. Keplerian propagation (always available, lowest precision)
```

### Enabling auto-download

```python
import libephemeris as swe

swe.set_auto_spk_download(True)
pos, _ = swe.calc_ut(2460000.0, swe.CHIRON, 0)
```

SPK files are cached in `~/.libephemeris/spk/` and reused on subsequent
calls.

### Strict precision mode

Strict precision (enabled by default) prevents the library from silently
falling back to low-precision Keplerian calculations for bodies that have
downloadable SPK kernels. When an SPK is missing for such a body, the
library raises `SPKRequiredError` with instructions on how to obtain it.

Bodies for which JPL Horizons does not provide SPK generation (see
[Known Limitations](#known-limitations)) are exempt from this check and
fall through to ASSIST or Keplerian.

```python
swe.set_strict_precision(False)   # Allow Keplerian fallback for all bodies
swe.set_strict_precision(True)    # Require SPK where available (default)
```

---

## Optional Extras

### `nbody` -- REBOUND/ASSIST n-body integration

```bash
pip install libephemeris[nbody]
```

Installs [REBOUND](https://rebound.readthedocs.io/) (N-body integrator)
and [ASSIST](https://assist.readthedocs.io/) (JPL-quality small body
extension). Together they propagate asteroid orbits accounting for:

- Gravitational perturbations from Sun, Moon, and all 8 planets
- 16 massive asteroid perturbers (Ceres, Vesta, Pallas, etc.)
- Earth/Sun gravitational harmonics (J2, J3, J4)
- General relativistic corrections

**When you need it:** When a numerical perturbation model is preferable to a
Keplerian fallback and no direct body SPK is available. It is not equivalent to
a direct JPL minor-body trajectory, and accuracy cannot be generalized across
all bodies or arbitrarily long propagation spans. The complete source interval
is the intersection of the selected planet and asteroid-perturber files.

**Data files** (~714 MB total, downloaded separately):

```bash
libephemeris download assist
```

| File | Size | Content |
|------|------|---------|
| `linux_p1550p2650.440` | ~98 MB | Planet ephemeris for ASSIST |
| `sb441-n16.bsp` | ~616 MB | 16 asteroid perturbers |

Extended offline LEB generation also uses
`linux_m13000p17000.441` (~2.6 GB), but its complete N-body interval remains
bounded by `sb441-n16.bsp`. The data-v3 generator reads both headers and applies
an integration guard before publishing per-body coverage.

See [REBOUND Integration](../methodology/rebound-integration.md) for
technical details and precision benchmarks.

### `stars` -- Star catalog

```bash
pip install libephemeris[stars]
```

Installs [Astropy](https://www.astropy.org/) and Astroquery for star-catalog
building and related developer tooling. Runtime fixed-star calculations use the
catalog shipped by LibEphemeris and do not require this extra.

### `all` -- Permissive optional runtime features

```bash
pip install libephemeris[all]        # permissive-only; add ,nbody for the GPL extra (provenance-ok)
```

Installs the permissive optional runtime features (currently Astropy
integration). It does not include the Astroquery-based catalog tooling from
`stars`. The GPL-licensed `nbody` extra is a separate explicit opt-in. <!-- provenance-ok -->

### `dev` -- Development tools

```bash
pip install libephemeris[dev]
```

For contributors only. Includes testing tools (pytest and plugins), code
quality (ruff, mypy), SPK generation (spiceypy), and manual building
(ebooklib). Not needed for using the library. (Reference-comparison
tooling lives in the separate `validation/` repository, not in any extra
of this package.)

---

## Data Requirements

Different setups require different data downloads:

| Setup | What to download | Total size | Command |
|-------|-----------------|------------|---------|
| **Minimal** | Nothing (uses the bundled LEB2 base core where covered) | 0 up front | -- |
| **Recommended** | Medium tier DE440 + planet centers + SPKs | >300 MB plus minor-body SPKs | `libephemeris download medium` |
| **High precision asteroids** | Above + ASSIST data | ~900 MB | Above + `libephemeris download assist` |
| **Binary ephemeris (fast)** | Reviewed pinned LEB2 core | ~38 MB (medium) | `libephemeris download leb-medium` |
| **Offline everything** | DE + SPKs + IERS; generate optional LEB locally | depends on tiers | `libephemeris download all` |
| **Offline + n-body** | Above + ASSIST data | ~6-7 GB | Above + `libephemeris download assist` |

---

## Known Limitations

### Bodies without auto-downloadable SPK

Most of the 37 minor bodies in the SPK map can be automatically downloaded
from JPL Horizons. However, some bodies are blocked by JPL:

| Body | Reason | Alternatives |
|------|--------|-------------|
| **Bennu** (101955) | JPL Horizons blocks SPK generation | ASSIST n-body, Keplerian fallback, or manual SPK registration |

For Bennu, the library skips the auto-download attempt and the strict
precision check, allowing the fallback chain to proceed to ASSIST or
Keplerian. If you need sub-arcsecond Bennu positions, install the `nbody`
extra and download ASSIST data.

A mission-specific NAIF kernel (`bennu_refdrmc_v1.bsp` from OSIRIS-REx,
covering 2015-2023) exists and can be registered manually:

```python
import libephemeris as swe

swe.register_spk_body(
    ipl=swe.BENNU,
    spk_file="/path/to/bennu_refdrmc_v1.bsp",
    naif_id=2101955,
)
pos, _ = swe.calc_ut(2458849.5, swe.BENNU, 0)  # 2020-01-01
```

### Keplerian precision varies by body class

| Body class | Keplerian error | With ASSIST |
|------------|----------------|-------------|
| Main belt (Ceres, Vesta) | ~10-30 arcsec | <1 arcsec |
| Near-Earth (Bennu, Apophis) | ~30-100 arcsec | <1 arcsec |
| Trans-Neptunian (Eris, Sedna) | ~50-200 arcsec | <1 arcsec |
| Centaurs (Chiron, Pholus) | ~10-50 arcsec | <1 arcsec |

These errors grow with distance from the orbital elements epoch. For
applications requiring better than arcminute precision, use SPK kernels
or ASSIST.

---

## Choosing the Right Setup

**Astrology (natal charts, transits):**
Planets are always high-precision. For the major asteroids (Chiron, Ceres,
Pallas, Juno, Vesta), enable auto-download or download SPKs once:

```bash
libephemeris download medium
```

```python
import libephemeris as swe
swe.set_auto_spk_download(True)
```

**Research / high-precision asteroid work:**
Install the n-body extra for sub-arcsecond precision on all minor bodies:

```bash
pip install libephemeris[nbody]
libephemeris download medium
libephemeris download assist
```

**Offline / production deployment:**
Download everything upfront:

```bash
libephemeris download all
```

For fastest calculations beyond the bundled base core, generate a local LEB
from the selected JPL kernel:

```bash
leph leb generate medium groups
```

```python
import libephemeris as swe
from pathlib import Path

swe.set_leb_file(str(Path("~/.libephemeris/leb/medium_core.leb2").expanduser()))
```
