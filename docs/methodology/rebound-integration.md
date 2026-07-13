# REBOUND Integration

LibEphemeris integrates REBOUND (an N-body integrator for orbital dynamics)
with ASSIST (a small-body ephemeris extension) to provide high-precision orbit
propagation for asteroids and other minor bodies when an SPK does not cover the
requested epoch. It replaces the last-resort Keplerian approximation with full
N-body dynamics; strict precision and JPL SPKs remain the default path for
mapped bodies.

## Table of Contents

- [Background](#background)
- [Method](#method)
  - [Keplerian vs N-body Integration](#keplerian-vs-n-body-integration)
  - [Error Growth Over Time](#error-growth-over-time)
  - [Available Integrators](#available-integrators)
- [API Reference](#api-reference)
- [Precision and Validation](#precision-and-validation)
- [When to Use REBOUND/ASSIST](#when-to-use-reboundassist)
- [Installation](#installation)
- [References](#references)

## Background

Minor body positions (asteroids, comets, TNOs) are typically computed by propagating orbital elements forward or backward from a reference epoch. The simplest approach uses Keplerian (two-body) mechanics, which treats the minor body as orbiting the Sun in an unperturbed ellipse. This is fast but ignores gravitational perturbations from the major planets, leading to error that grows with time from epoch.

N-body integration solves the full equations of motion including perturbations from all major planets, the Moon, and optionally the largest asteroids. REBOUND is a widely used open-source N-body integrator, and ASSIST extends it with JPL ephemeris-based planet positions and non-gravitational force models. Together they enable LibEphemeris to achieve sub-arcsecond minor body positions over decades-long propagation spans.

## Method

### Keplerian vs N-body Integration

| Method | Precision | Use Case |
|--------|-----------|----------|
| Keplerian elements | ~10-30 arcsec (main belt) | Fast, good for near-epoch dates |
| N-body (REBOUND) | <1 arcsec | High precision, accounts for perturbations |

The Keplerian approximation provides the following typical precision by body class:
- **Main belt asteroids**: ~10-30 arcsec
- **Near-Earth asteroids**: ~30-100 arcsec
- **Trans-Neptunian objects**: ~50-200 arcsec

N-body integration with REBOUND/ASSIST reduces these errors by 10-100x depending on the body and time span.

### Error Growth Over Time

For a representative main belt asteroid (Ceres) at various intervals from epoch:

| Years from Epoch | Keplerian Error | N-body Error |
|------------------|-----------------|--------------|
| 0 | ~0.1 arcsec | ~0.001 arcsec |
| 1 | ~1 arcsec | ~0.01 arcsec |
| 10 | ~10 arcsec | ~0.1 arcsec |
| 100 | ~100 arcsec | ~1 arcsec |

### Available Integrators

REBOUND provides multiple integrators for different use cases:

| Integrator | Best For | Order |
|------------|----------|-------|
| IAS15 | Highest precision, adaptive timestep | 15th order |
| WHFast | Long-term stability, symplectic | 2nd order |
| MERCURIUS | Close encounters, hybrid | Adaptive |

#### IAS15 (Recommended for Precision)

- **Use for**: Highest precision requirements
- **Precision**: 15th order, adaptive timestep
- **Performance**: Slower but most accurate
- **Best for**: Research, close approaches, non-gravitational forces

#### WHFast (Recommended for Speed)

- **Use for**: Long-term stability studies
- **Precision**: 2nd order symplectic
- **Performance**: Very fast
- **Best for**: Century-scale propagation of stable orbits

#### MERCURIUS (Recommended for Close Encounters)

- **Use for**: Orbits with close planetary encounters
- **Precision**: Hybrid symplectic/Bulirsch-Stoer
- **Performance**: Moderate
- **Best for**: Near-Earth asteroids, potentially hazardous objects

## API Reference

LibEphemeris exposes the following REBOUND-related functions and classes:

| Function/Class | Description |
|----------------|-------------|
| `check_rebound_available()` | Returns True if REBOUND is installed |
| `check_assist_available()` | Returns True if ASSIST is installed |
| `get_rebound_version()` | Returns REBOUND version string |
| `get_assist_version()` | Returns ASSIST version string |
| `propagate_orbit_rebound()` | Propagate an orbit with REBOUND (2-body, Sun only) |
| `propagate_orbit_assist()` | Propagate an orbit with ASSIST (ephemeris-quality N-body) |
| `propagate_trajectory()` | Propagate an orbit to positions at multiple times |
| `ReboundIntegrator` | Enum selecting the REBOUND integrator (IAS15, WHFAST, MERCURIUS, TRACE, LEAPFROG) |
| `AssistEphemConfig` | Configuration for ASSIST ephemeris |
| `PropagationResult` | Result container for propagated orbits |
| `elements_to_rebound_orbit()` | Convert orbital elements to REBOUND orbit |

### Propagating an Orbit

```python
from libephemeris.minor_bodies import MINOR_BODY_ELEMENTS, CERES
from libephemeris.rebound_integration import (
    propagate_orbit_assist,
    propagate_orbit_rebound,
    ReboundIntegrator,
)

elements = MINOR_BODY_ELEMENTS[CERES]
jd_start = 2460000.5   # starting epoch (TT)
jd_end = 2460365.5     # one year later

# Ephemeris-quality N-body propagation (Sun, planets, Moon, major asteroids)
result = propagate_orbit_assist(elements, jd_start, jd_end)
print(f"Ceres: lon={result.ecliptic_lon:.6f} deg, dist={result.distance:.6f} AU")

# 2-body REBOUND propagation with an explicit integrator selection
result = propagate_orbit_rebound(
    elements, jd_start, jd_end, integrator=ReboundIntegrator.IAS15
)
x, y, z, vx, vy, vz = result.to_tuple()
```

## Precision and Validation

### Orbit Propagation Accuracy

The main precision benefits of N-body integration are:

1. **Perturbation handling**: N-body integration accounts for gravitational perturbations from major planets.
2. **Long-term stability**: Orbits stay bounded over centuries to millennia.
3. **Non-gravitational forces**: ASSIST can model solar radiation pressure, Yarkovsky effect, and outgassing.

### Precision at Epoch

For Ceres at its orbital epoch:
- Keplerian position error: ~0.1 arcsec
- N-body position error: ~0.001 arcsec (100x improvement)

## When to Use REBOUND/ASSIST

### Recommended Use Cases

1. **High-precision asteroid positions**: When arcsecond or better accuracy is required.

2. **Long-term orbit propagation**: Projecting orbits years to centuries from epoch.

3. **Close approach calculations**: Asteroids near Earth or other planets.

4. **Research applications**: Scientific studies requiring N-body dynamics.

5. **Non-gravitational effects**: Modeling solar radiation pressure, outgassing, etc.

### When Keplerian Is Sufficient

1. **Near-epoch dates**: Within a few years of the orbital element epoch.

2. **Main belt asteroids**: For dates close to epoch, Keplerian is often adequate.

3. **Performance-critical applications**: Keplerian is ~1000x faster.

4. **Approximate positions**: When arcsecond precision is not required.

## Installation

REBOUND and ASSIST are optional dependencies. Install them with:

```bash
pip install rebound assist
```

Or with LibEphemeris extras:

```bash
pip install libephemeris[nbody]
```

To verify availability:

```python
from libephemeris.rebound_integration import check_rebound_available, check_assist_available
if check_rebound_available() and check_assist_available():
    print("REBOUND/ASSIST available for high-precision orbit propagation")
```

When REBOUND/ASSIST is not installed, LibEphemeris falls back to Keplerian propagation.

## References

- Rein, H. & Liu, S.-F. "REBOUND: An open-source multi-purpose N-body code for collisional dynamics" (2012), Astronomy & Astrophysics, 537, A128
- Rein, H. & Spiegel, D.S. "IAS15: a fast, adaptive, high-order integrator for gravitational dynamics, accurate to machine precision over a billion orbits" (2015), MNRAS, 446, 1424-1437
- Holman, M.J. et al. "ASSIST: An Ephemeris-quality Test-particle Integrator" (2023), The Planetary Science Journal, 4, 69
- Rein, H. & Tamayo, D. "WHFast: a fast and unbiased implementation of a symplectic Wisdom-Holman integrator for long-term gravitational simulations" (2015), MNRAS, 452, 376-388
