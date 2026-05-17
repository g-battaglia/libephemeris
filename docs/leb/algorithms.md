# LEB Algorithms & Mathematical Foundations

> **Version:** 1.0 — March 2026
> **Audience:** Developers working on the LEB system, or anyone wanting to
> understand the mathematical and computational techniques behind the binary
> ephemeris format.

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Chebyshev Polynomials](#2-chebyshev-polynomials)
3. [The Clenshaw Algorithm](#3-the-clenshaw-algorithm)
4. [Chebyshev Fitting (Generation)](#4-chebyshev-fitting-generation)
5. [Analytical Velocity via Chebyshev Derivatives](#5-analytical-velocity-via-chebyshev-derivatives)
6. [Coordinate Systems and Storage Strategies](#6-coordinate-systems-and-storage-strategies)
7. [The ICRS Pipeline: From Barycentric to Apparent](#7-the-icrs-pipeline-from-barycentric-to-apparent)
8. [Light-Time Correction](#8-light-time-correction)
9. [Gravitational Deflection](#9-gravitational-deflection)
10. [Stellar Aberration](#10-stellar-aberration)
11. [Precession and Nutation](#11-precession-and-nutation)
12. [Center-of-Body (COB) Corrections](#12-center-of-body-cob-corrections)
13. [Longitude Unwrapping](#13-longitude-unwrapping)
14. [Delta-T (TT - UT1)](#14-delta-t-tt---ut1)
15. [Error Analysis and Precision Budget](#15-error-analysis-and-precision-budget)
16. [Historical Problems and Solutions](#16-historical-problems-and-solutions)

---

## 1. Introduction

The LEB (LibEphemeris Binary) system precomputes celestial body positions as
Chebyshev polynomial approximations and stores them in a compact binary file.
At runtime, evaluating a body's position reduces to:

1. An O(1) segment lookup (integer division)
2. A Clenshaw evaluation of the Chebyshev series (~0.5 us per component)
3. Coordinate transforms and physical corrections (light-time, aberration,
   deflection, precession-nutation)

This document explains each of these steps in mathematical detail.

### What is an Ephemeris?

An *ephemeris* (plural: *ephemerides*) is a table or function that gives the
positions of astronomical objects at specific times. The NASA JPL Development
Ephemeris (DE440/DE441) is the gold standard: it encodes planetary positions
as Chebyshev polynomial segments fitted to numerical integrations of the
solar system's equations of motion.

LEB applies the same technique one level higher: it fits Chebyshev polynomials
to the *output* of the Skyfield/JPL pipeline (which itself evaluates JPL's
Chebyshev polynomials, applies frame rotations, and computes apparent positions).
This creates a purpose-built cache of pre-transformed results.

### Why Chebyshev Polynomials?

Chebyshev polynomials are the optimal choice for function approximation because:

- They minimize the maximum approximation error (minimax property)
- They are numerically stable to evaluate (no Runge phenomenon)
- They can be evaluated efficiently via the Clenshaw algorithm
- Their derivatives have a simple closed-form recurrence
- They are the same technique used by JPL internally

The key insight is that celestial body positions are smooth, slowly-varying
functions of time (with some exceptions noted below). A degree-13 Chebyshev
polynomial over a 4-32 day interval can approximate these functions to
better than 0.001 arcsecond accuracy.

---

## 2. Chebyshev Polynomials

### Definition

The Chebyshev polynomial of the first kind, T_n(x), is defined on [-1, 1] by:

```
T_0(x) = 1
T_1(x) = x
T_n(x) = 2x * T_{n-1}(x) - T_{n-2}(x)    for n >= 2
```

Equivalently, T_n(cos(theta)) = cos(n * theta).

The first several polynomials are:

```
T_0(x) = 1
T_1(x) = x
T_2(x) = 2x^2 - 1
T_3(x) = 4x^3 - 3x
T_4(x) = 8x^4 - 8x^2 + 1
T_5(x) = 16x^5 - 20x^3 + 5x
```

### Chebyshev Series

A function f(x) on [-1, 1] can be approximated by a truncated Chebyshev series:

```
f(x) ~ c_0 * T_0(x) + c_1 * T_1(x) + ... + c_N * T_N(x)
```

where the coefficients c_k are chosen to minimize the approximation error.
This is analogous to a Fourier series but for polynomial approximation on
a finite interval.

### Domain Mapping

In LEB, each Chebyshev segment covers a time interval [jd_start, jd_end].
To evaluate the polynomial, the Julian Day must be mapped to the normalized
domain [-1, 1]:

```
tau = 2 * (jd - jd_mid) / interval_days
```

where:
- `jd_mid = jd_start + interval_days / 2` (segment midpoint)
- `interval_days = jd_end - jd_start` (segment width)

This mapping is crucial: the Clenshaw algorithm operates on tau in [-1, 1],
and the domain scaling factor appears in the velocity computation.

### Why Not Regular Polynomials?

Standard power-series polynomials (a_0 + a_1*x + a_2*x^2 + ...) suffer from:

1. **Runge's phenomenon**: High-degree polynomial interpolation on uniformly
   spaced points oscillates wildly near the interval boundaries.
2. **Numerical instability**: Evaluating x^13 in floating-point arithmetic
   accumulates significant rounding errors.
3. **Poor conditioning**: The Vandermonde matrix used for fitting becomes
   nearly singular at high degree.

Chebyshev polynomials avoid all three problems because:
- They use Chebyshev nodes (cosine-spaced), which cluster near the boundaries
- The Clenshaw recurrence uses only multiplication and addition (no powers)
- The Chebyshev basis is orthogonal, producing well-conditioned least-squares fits

---

## 3. The Clenshaw Algorithm

### Position Evaluation

The Clenshaw algorithm evaluates a Chebyshev series without explicitly
computing T_n(x). Given coefficients (c_0, c_1, ..., c_N) and evaluation
point tau in [-1, 1]:

```
Initialize:
    b_{N+1} = 0
    b_{N+2} = 0

Recurrence (k = N, N-1, ..., 1):
    b_k = c_k + 2*tau * b_{k+1} - b_{k+2}

Result:
    f(tau) = c_0 + tau * b_1 - b_2
```

This requires only 2 temporary variables and N multiply-add operations.
No arrays are allocated. For degree 13, this is 13 iterations -- about
0.5 microseconds in Python.

**Implementation** (`leb_reader.py:59-79`):

```python
def _clenshaw(coeffs: tuple[float, ...], tau: float) -> float:
    n = len(coeffs) - 1
    if n == 0:
        return coeffs[0]
    b_k1 = 0.0  # b_{k+1}
    b_k2 = 0.0  # b_{k+2}
    for k in range(n, 0, -1):
        b_k = coeffs[k] + 2.0 * tau * b_k1 - b_k2
        b_k2 = b_k1
        b_k1 = b_k
    return coeffs[0] + tau * b_k1 - b_k2
```

### Why Pure Python?

For single-point evaluation (the common case in `swe_calc_ut()`), numpy array
creation overhead (~5 us) would dominate the actual Clenshaw loop (~0.5 us for
degree 13). Pure Python `float` operations avoid this overhead entirely.

For batch evaluation during generation, numpy's vectorized `chebval()` is used
instead.

---

## 4. Chebyshev Fitting (Generation)

### Chebyshev Nodes

The generation process samples the reference function at **Chebyshev nodes**
(also called Chebyshev-Gauss points or Type I nodes):

```
x_k = cos(pi * (k + 0.5) / n)    for k = 0, 1, ..., n-1
```

where n = degree + 1 (the number of nodes equals the number of coefficients).

These nodes are not uniformly spaced -- they cluster near the boundaries of
[-1, 1]. This clustering is precisely what prevents Runge's phenomenon and
ensures near-optimal approximation.

**Mapped to Julian Days:**

```python
jd_nodes = 0.5 * (jd_end - jd_start) * chebyshev_nodes + 0.5 * (jd_start + jd_end)
```

### Fitting Procedure

For each segment, the fitting process is:

1. **Sample**: Evaluate the reference function (Skyfield, analytical formula,
   or SPK kernel) at the n = degree + 1 Chebyshev nodes.
2. **Fit**: Compute coefficients using `numpy.polynomial.chebyshev.chebfit()`.
   This solves the least-squares problem in the Chebyshev basis.
3. **Verify**: Evaluate the fitted polynomial at 10 uniformly-spaced test
   points (NOT on the Chebyshev nodes) and compare against the reference
   function. Track maximum error.

The verification step is critical: it catches problems like longitude
wrap-around, SPK boundary overshoot, or insufficient degree.

### Vectorized Evaluation

The major performance optimization in the generator is **batching all JDs
across all segments into a single Skyfield evaluation call**:

```
Step 1: Precompute ALL Julian Days needed:
  For each segment i:
    (degree+1) Chebyshev fit nodes
    10 uniform verification points
  Total: n_segments * (degree + 1 + 10) JDs

Step 2: Single vectorized Skyfield call:
  positions = target.at(ts.tt_jd(all_jds)).position.au.T  # (N, 3)

Step 3: Redistribute values and fit each segment independently
```

This eliminates Skyfield's per-call overhead (time conversion, SPK segment
lookup, Python function call dispatch) and achieves ~150x speedup for planets.

### Segment Width Selection

The segment width (interval_days) controls the trade-off between file size
and approximation accuracy:

| Shorter intervals | Longer intervals |
|-------------------|------------------|
| Better accuracy | Worse accuracy |
| More segments | Fewer segments |
| Larger file | Smaller file |

The optimal width depends on how rapidly the body's position changes:

- **Moon (4 days)**: Moves ~13 deg/day, orbital period 27.3 days.
  A 4-day segment spans ~52 degrees of lunar motion.
- **Sun/EMB (32 days)**: Moves ~1 deg/day, very smooth motion.
- **Uranus (64 days)**: Moves ~0.01 deg/day in barycentric coordinates.
  The apparent geocentric motion is faster due to Earth's parallax.
- **Mercury (16 days, degree 15)**: Most eccentric planet orbit (e=0.206).
  Needs higher degree to capture the non-uniform orbital velocity.

---

## 5. Analytical Velocity via Chebyshev Derivatives

### The Derivative Recurrence

The derivative of a Chebyshev series has a simple closed-form expression.
Given coefficients (c_0, c_1, ..., c_N), the derivative coefficients
(d_0, d_1, ..., d_{N-1}) are computed via the recurrence:

```
d_{N-1} = 2*N * c_N
d_k     = d_{k+2} + 2*(k+1) * c_{k+1}    for k = N-2, N-3, ..., 1
d_0     = d_2 / 2 + c_1
```

The derivative polynomial is then evaluated via the same Clenshaw algorithm.

**Implementation** (`leb_reader.py:82-112`):

```python
def _deriv_coeffs(coeffs: tuple[float, ...]) -> tuple[float, ...]:
    n = len(coeffs) - 1
    if n == 0:
        return (0.0,)
    d = [0.0] * n
    d[n - 1] = 2.0 * n * coeffs[n]
    for k in range(n - 2, 0, -1):
        d[k] = d[k + 2] + 2.0 * (k + 1) * coeffs[k + 1]
    d[0] = d[2] / 2.0 + coeffs[1] if n >= 2 else coeffs[1]
    return tuple(d)
```

### Domain Scaling

The derivative coefficients give d/d(tau), but we need d/d(jd). The chain
rule gives:

```
d/d(jd) = d/d(tau) * d(tau)/d(jd) = d/d(tau) * 2 / interval_days
```

So the velocity is:

```python
velocity = clenshaw(deriv_coeffs, tau) * 2.0 / body.interval_days
```

### Advantages Over Central Difference

The previous approach computed velocity via central difference:

```
v(t) = (f(t + dt) - f(t - dt)) / (2 * dt)
```

This required **two additional pipeline evaluations** per velocity component
(6 extra evaluations for 3D velocity with speed). The analytical derivative:

- Is **3x faster** (1 pipeline run instead of 3)
- Is **more precise** (no finite-difference truncation error)
- Does **not amplify Chebyshev fitting errors** (the derivative is computed
  from the polynomial itself, not from perturbed evaluations)

---

## 6. Coordinate Systems and Storage Strategies

### The Three Coordinate Types

LEB stores body positions in one of five coordinate frames, chosen to minimize
runtime computation and maximize Chebyshev fitting accuracy:

#### COORD_ICRS_BARY (type 0) — Planet Centers

**Used for:** Sun, Moon, Mercury, Venus, Mars, Earth, Chiron, Ceres-Vesta

Stores the **ICRS barycentric position** of the planet center in AU.
ICRS (International Celestial Reference System) is an inertial frame
centered at the Solar System Barycenter (SSB), with axes aligned to
distant quasars (effectively fixed in space).

**Why ICRS?** A single dataset in ICRS supports ALL output coordinate frames:
- Geocentric ecliptic of date (the default, most common)
- Geocentric equatorial of date
- J2000 ecliptic
- J2000 equatorial (ICRS itself)
- Heliocentric ecliptic
- Barycentric
- Sidereal (any ayanamsa)

If positions were stored in ecliptic-of-date coordinates, each output frame
would need its own precomputed dataset.

**Why planet centers, not barycenters?** For inner planets (Mercury-Mars) and
the Sun, the planet center IS the barycentric position (no moons or negligible
moon mass). For outer planets, see COORD_ICRS_BARY_SYSTEM below.

#### COORD_ICRS_BARY_SYSTEM (type 4) — System Barycenters with Runtime COB

**Used for:** Jupiter, Saturn, Uranus, Neptune, Pluto

Stores the **system barycenter** (the gravitational center of the planet plus
all its moons) in ICRS AU. The Center-of-Body (COB) correction — the offset
from system barycenter to planet center — is applied at runtime.

This was the critical innovation that enabled <0.001" precision for outer
planets. See [Section 12](#12-center-of-body-cob-corrections) for full details.

#### COORD_ECLIPTIC (type 1) — Ecliptic of Date

**Used for:** Mean Node, True Node, Mean Apogee (Lilith), Osculating Apogee,
Interpolated Apogee, Interpolated Perigee

Stores (longitude, latitude, distance) in degrees/degrees/AU in the ecliptic
coordinate system of the date. These bodies are computed by analytical
formulas that directly output ecliptic coordinates; storing them in ICRS
would require an unnecessary inverse rotation.

#### COORD_HELIO_ECL (type 2) — Heliocentric Ecliptic

**Used for:** Cupido, Hades, Zeus, Kronos, Apollon, Admetos, Vulkanus,
Poseidon, Transpluto (Uranian hypothetical bodies)

Stores (longitude, latitude, distance) in degrees/degrees/AU in heliocentric
ecliptic coordinates. These bodies are defined by Keplerian orbital elements
relative to the Sun.

#### COORD_GEO_ECLIPTIC (type 3) — Reserved

Defined in the format but **not used by any body**. Originally planned for
storing geocentric ecliptic positions, but abandoned because ecliptic
coordinates have discontinuities at retrograde stations that Chebyshev
polynomials cannot fit accurately.

### Why Not Store Everything in Ecliptic?

Early versions of LEB attempted to store planet positions in geocentric
ecliptic coordinates (`COORD_GEO_ECLIPTIC`). This failed because:

1. **Retrograde cusps**: When a planet stations (changes from direct to
   retrograde motion), its ecliptic longitude has a cusp — a sharp
   reversal. Chebyshev polynomials, being smooth, cannot fit cusps.
   The fitting error at stations was 3-5 arcseconds.

2. **Latitude sign changes**: Near opposition, a planet's ecliptic latitude
   can change sign rapidly, creating another fitting challenge.

The ICRS barycentric frame avoids both problems: planetary motion in ICRS
is always smooth (no retrograde, no cusps). The coordinate transforms that
produce the apparent ecliptic position are applied at runtime, after the
smooth ICRS position has been accurately recovered from the Chebyshev fit.

---

## 7. The ICRS Pipeline: From Barycentric to Apparent

The ICRS pipeline (Pipeline A/A') transforms raw barycentric positions into
apparent geocentric coordinates. This is the most complex part of LEB and
replicates what Skyfield does internally.

### Step-by-Step Pipeline

```
1. BODY POSITION
   Read (x, y, z) in AU from LEB Chebyshev data
   Read Earth position from LEB

2. GEOMETRIC VECTOR
   geo = body_position - earth_position

3. LIGHT-TIME CORRECTION (Section 8)
   Iterate: lt = |geo| / c, re-evaluate body at t - lt

4. GRAVITATIONAL DEFLECTION (Section 9)
   Apply PPN deflection by Sun, Jupiter, Saturn

5. STELLAR ABERRATION (Section 10)
   Apply classical aberration using Earth velocity

6. FRAME ROTATION (Section 11)
   ICRS -> equatorial of date (precession-nutation matrix)
   equatorial -> ecliptic (obliquity rotation)

7. SPHERICAL CONVERSION
   (x, y, z) -> (longitude, latitude, distance)

8. SIDEREAL CORRECTION (if requested)
   longitude -= ayanamsa
```

Each step is detailed in its own section below.

---

## 8. Light-Time Correction

### Physical Basis

Light travels at a finite speed (c = 173.14 AU/day = 299,792,458 m/s).
When we observe a planet at time t, we see it at the position it occupied
at time t - lt, where lt is the light travel time:

```
lt = |body(t - lt) - observer(t)| / c
```

This is an implicit equation: the light-time depends on the position, which
depends on the light-time. It is solved by fixed-point iteration.

### Implementation

```python
C_LIGHT_AU_DAY = 173.1446326846693  # AU/day

# Initial geometric vector (no light-time)
geo = body_pos - earth_pos

# 3 fixed-point iterations (converges to ~1e-15 AU)
for _ in range(3):
    dist = sqrt(geo[0]**2 + geo[1]**2 + geo[2]**2)
    lt = dist / C_LIGHT_AU_DAY
    retarded_pos, _ = reader.eval_body(ipl, jd_tt - lt)
    geo = retarded_pos - earth_pos
```

Three iterations are sufficient because:
- The speed of planets is ~1e-4 c (much less than light speed)
- Each iteration reduces the error by a factor of ~v/c ~ 1e-4
- After 3 iterations: error ~ (v/c)^3 ~ 1e-12, well below 0.001"

### Light-Time for System Barycenters

For bodies stored as `COORD_ICRS_BARY_SYSTEM` (Jupiter-Pluto), the light-time
iteration uses the raw system barycenter (smooth, easy to interpolate). The COB
correction is applied **after** the light-time iteration converges, at
**observer time** (not retarded time). This matches Skyfield's behavior in
`_SpkCenterTarget._observe_from_bcrs()`.

---

## 9. Gravitational Deflection

### Physical Basis

General relativity predicts that massive bodies deflect light rays passing
near them. The Sun deflects starlight by up to 1.75 arcseconds at the solar
limb. For planets observed from Earth, the maximum deflection is:

| Deflector | Maximum deflection |
|-----------|-------------------|
| Sun | ~0.004" (for bodies near the Sun) |
| Jupiter | ~0.017" (for bodies near Jupiter's limb) |
| Saturn | ~0.006" |

These are small but significant compared to the 0.001" precision target.

### PPN Formula

The Parameterized Post-Newtonian (PPN) deflection formula, as implemented
in Skyfield's `apparent()` method:

```
deflection = (1 + gamma) * G*M / (c^2 * d) * (e_hat x (e_hat x q_hat))
```

where:
- gamma = 1 (general relativity)
- G*M is the gravitational parameter of the deflector
- c is the speed of light
- d is the closest approach distance of the light ray to the deflector
- e_hat is the unit vector from deflector to body
- q_hat is the unit vector from deflector to observer
- x denotes the vector cross product

### Implementation in LEB

**Deflectors** (`fast_calc.py:262`):

```python
_DEFLECTORS = (
    (SE_SUN,  1.0),       # Sun, mass ratio = 1.0 (reference)
    (5,       1047.3486),  # Jupiter barycenter, Sun/Jupiter mass ratio
    (6,       3497.898),   # Saturn barycenter, Sun/Saturn mass ratio
)
```

The gravitational parameter of each deflector is computed as:
```
GM_deflector = GM_sun / mass_ratio
```

where GM_sun = 1.32712440017987 x 10^20 m^3/s^2.

**When applied:**
- Only for `COORD_ICRS_BARY` and `COORD_ICRS_BARY_SYSTEM` bodies
- Skipped for the Moon (too close; deflection formula breaks down)
- Skipped for heliocentric, barycentric, true position, and no-aberration modes

**Historical note:** The absence of gravitational deflection in early LEB
versions caused a 3.95" error for Saturn. Adding PPN deflection for Sun,
Jupiter, and Saturn reduced all planet errors to <0.001".

---

## 10. Stellar Aberration

### Physical Basis

Stellar aberration is a relativistic effect caused by the observer's velocity
relative to the incoming light. It shifts the apparent position of a body
in the direction of the observer's motion. The maximum annual aberration
(due to Earth's orbital velocity of ~30 km/s) is about 20.5 arcseconds.

### Classical First-Order Formula

LEB uses the classical (first-order in v/c) aberration formula:

```
u = geo / |geo|                   # unit vector to body
v = earth_vel / c                 # Earth velocity in natural units
u' = u + v - u * (u . v)         # aberrated direction
result = normalize(u') * |geo|    # scale back to original distance
```

This matches the pyswisseph implementation. The rigorous special-relativistic
formula differs by less than 1 milliarcsecond, which is negligible for our
0.001" precision target.

### When Applied

- Only for geocentric calculations (default mode)
- Skipped for heliocentric (`FLG_HELCTR`), barycentric (`FLG_BARYCTR`),
  true position (`FLG_TRUEPOS`), and no-aberration (`FLG_NOABERR`) modes

---

## 11. Precession and Nutation

### Physical Basis

The Earth's rotation axis is not fixed in space. It undergoes two motions:

1. **Precession**: A slow, smooth 25,772-year cycle caused by the Sun and
   Moon's gravitational torque on Earth's equatorial bulge. The rotation
   axis traces a cone in space.

2. **Nutation**: Short-period oscillations superimposed on precession, caused
   by the Moon's orbital plane precessing with an 18.6-year period.

Together, these determine the orientation of the "equator of date" and
"ecliptic of date" reference frames relative to the fixed ICRS frame.

### Precession-Nutation Matrix

The combined effect is encoded in a 3x3 rotation matrix that transforms
vectors from ICRS to the equatorial frame of a given date:

```python
# Primary method: PyERFA (IAU 2006 precession + IAU 2000A nutation)
import erfa
pn_matrix = erfa.pnm06a(2451545.0, jd_tt - 2451545.0)

# Fallback: libephemeris internal implementation
from astrometry import _precession_nutation_matrix as _pnm
pn_matrix = _pnm(jd_tt)
```

The matrix is applied to both position and velocity vectors.

### Nutation in LEB

Nutation angles (dpsi, deps) are stored as Chebyshev polynomial segments
in the LEB file (Section 2 of the binary format):

- **Interval**: 16 days
- **Degree**: 16
- **Components**: 2 (dpsi in radians, deps in radians)
- **Model**: IAU 2006/2000A (via `erfa.nut06a()` during generation)

These are used for:
1. Computing true obliquity: eps_true = eps_mean + deps
2. True ayanamsa: mean_ayanamsa + degrees(dpsi) (nutation in longitude)

### Mean Obliquity

The mean obliquity of the ecliptic (the angle between the equator and
ecliptic planes, ignoring nutation) is computed from the IAU 2006 polynomial:

```python
# Coefficients in arcseconds
_OBLIQUITY_COEFFS = (84381.406, -46.836769, -0.0001831,
                     0.00200340, -0.000000576, -0.0000000434)

def _mean_obliquity_iau2006(jd_tt: float) -> float:
    T = (jd_tt - 2451545.0) / 36525.0
    eps_arcsec = sum(c * T**i for i, c in enumerate(_OBLIQUITY_COEFFS))
    return eps_arcsec / 3600.0  # convert to degrees
```

### Ecliptic Rotation

The rotation from equatorial to ecliptic coordinates uses the true obliquity:

```
x_ecl =  x_eq
y_ecl =  y_eq * cos(eps) + z_eq * sin(eps)
z_ecl = -y_eq * sin(eps) + z_eq * cos(eps)
```

Then convert to spherical:
```
longitude = atan2(y_ecl, x_ecl)   # [0, 360) degrees
latitude  = asin(z_ecl / dist)    # [-90, 90] degrees
distance  = sqrt(x^2 + y^2 + z^2) # AU
```

---

## 12. Center-of-Body (COB) Corrections

### The Problem

JPL ephemerides provide positions for "planet barycenters" (NAIF IDs
5 = Jupiter barycenter, 6 = Saturn barycenter, etc.) — the gravitational
center of the planet plus all its moons. For astrological calculations,
we need the planet *center* — the physical center of the planet itself.

The offset between system barycenter and planet center is the COB correction.
For Jupiter, this offset can be up to ~0.01 AU (oscillating as Ganymede,
Callisto, Io, and Europa orbit). These oscillations have periods of 1.77-16.7
days and produce angular effects of 0.5-2 arcseconds.

### Why COB Breaks Chebyshev Fitting

The COB correction contains high-frequency oscillations from the inner moons.
Fitting planet_center = barycenter + COB with Chebyshev polynomials requires
very short intervals (< 1 day) to capture these oscillations. This would
produce enormous files and still leave residual fitting errors of 0.1-1".

### The Solution: Store Barycenters, Apply COB at Runtime

The `COORD_ICRS_BARY_SYSTEM` storage strategy separates the smooth and
oscillatory components:

1. **Generator**: Stores the pure system barycenter position (smooth, easy
   to fit with 32-64 day intervals and degree 13)
2. **Runtime**: Applies the COB correction using either:
   - `planet_centers.bsp` SPK segments (<0.001" precision), or
   - Analytical moon theory corrections (<0.01" precision) as fallback

### COB Evaluation Timing

**Critical**: The COB correction is evaluated at **observer time** (jd_tt),
not at **retarded time** (jd_tt - light_time). This matches Skyfield's
`_SpkCenterTarget._observe_from_bcrs()` behavior:

```python
# Skyfield's approach (simplified):
def _observe_from_bcrs(observer_pos, observer_t):
    # Light-time iteration on barycenter
    for _ in range(10):
        dist = |barycenter(t_retarded) - observer_pos|
        lt = dist / c
        t_retarded = observer_t - lt

    # COB correction at OBSERVER time (not retarded time)
    center_offset = planet_centers_segment.at(observer_t)
    target_pos = barycenter(t_retarded) + center_offset

    return target_pos - observer_pos
```

A previous bug in LEB evaluated COB at retarded time, producing ~0.002"
systematic errors for outer planets.

### SPK Planet Center Segments

High-precision COB data comes from `planet_centers_{tier}.bsp` files:

| NAIF ID | Planet | Typical offset magnitude |
|---------|--------|-------------------------|
| 599 | Jupiter center | ~0.001-0.01 AU |
| 699 | Saturn center | ~0.0001-0.001 AU |
| 799 | Uranus center | ~0.00001 AU |
| 899 | Neptune center | ~0.00001 AU |
| 999 | Pluto center | ~0.0001 AU |

These segments have limited date coverage per tier:

| Tier | File | Coverage varies by planet |
|------|------|---------------------------|
| Base | `planet_centers_base.bsp` | ~1849-2150 |
| Medium | `planet_centers_medium.bsp` | Varies (Jupiter: 1600-1997, Saturn: 1750-2014, etc.) |

When the SPK segment doesn't cover the requested date, the analytical
COB fallback is used automatically.

---

## 13. Longitude Unwrapping

### The Problem

Ecliptic longitude has a discontinuity at 0/360 degrees. When a body crosses
this boundary, the raw values jump from ~359 to ~1 (or vice versa). Fitting
a Chebyshev polynomial across such a jump would produce wildly incorrect
results — the polynomial would try to smoothly interpolate between 359 and 1,
passing through ~180 degrees.

### The Solution: Unwrap Before Fitting, Re-wrap After Evaluation

**During generation:**

```python
import numpy as np

# Raw longitude values with potential 360-degree jumps
raw_lon = [358.5, 359.2, 359.8, 0.3, 0.9, 1.5]

# Unwrap: remove 360-degree jumps
unwrapped = np.degrees(np.unwrap(np.radians(raw_lon)))
# Result: [358.5, 359.2, 359.8, 360.3, 360.9, 361.5]

# Now fit Chebyshev to the continuous unwrapped series
coeffs = chebfit(nodes, unwrapped, degree)
```

**During evaluation:**

```python
# Clenshaw gives the unwrapped value (can be > 360 or < 0)
raw_value = clenshaw(coeffs, tau)

# Re-wrap to [0, 360)
longitude = raw_value % 360.0
```

This approach is transparent: the user always sees longitude in [0, 360).

### Verification

The generator verifies that each segment correctly handles the wrap-around
by evaluating at 10 intermediate test points and comparing with the reference
function after re-wrapping. This catches cases where the unwrapping fails
(e.g., multiple 360-degree jumps within a single segment).

---

## 14. Delta-T (TT - UT1)

### What is Delta-T?

Delta-T is the difference between Terrestrial Time (TT, a uniform time scale
based on atomic clocks) and Universal Time (UT1, which tracks Earth's
irregular rotation). Since UT1 depends on Earth's actual rotation speed,
which is unpredictable, Delta-T cannot be computed analytically — it must
be measured and tabulated.

As of 2025, Delta-T is approximately 69.36 seconds. It was ~0 seconds around
1900 and is increasing irregularly.

### Storage in LEB

Delta-T values are stored as a sparse table of (JD, delta_t_days) pairs,
sampled every 30 days throughout the tier's date range:

```
Section 3: Delta-T Table
  Header: n_entries (uint32), reserved (uint32)
  Entries: [jd (float64), delta_t (float64)] x n_entries
```

### Interpolation

The reader uses **linear interpolation** between adjacent entries:

```python
idx = bisect_right(jds, jd) - 1
t = (jd - jds[idx]) / (jds[idx+1] - jds[idx])
delta_t = vals[idx] + t * (vals[idx+1] - vals[idx])
```

With 30-day spacing, linear interpolation introduces up to ~0.004 seconds
of error near epochs where Delta-T changes rapidly (e.g., around 1985).

### Delta-T Usage in fast_calc

**Critical**: The `fast_calc_ut()` entry point does **not** use
`reader.delta_t()` for the UT->TT conversion. Instead, it uses
`swe_deltat()` from `time_utils.py`, which provides the same high-precision
Delta-T model used by the Skyfield reference pipeline. This ensures exact
agreement with the reference:

```python
# fast_calc.py, fast_calc_ut():
from .time_utils import swe_deltat
delta_t = swe_deltat(tjd_ut)  # high-precision model
jd_tt = tjd_ut + delta_t
```

The reader's `delta_t()` method is still used by `fast_calc_tt()` for the
reverse (TT->UT) approximation needed by sidereal ayanamsa computation,
where the ~0.004s error is negligible.

---

## 15. Error Analysis and Precision Budget

### Error Sources

The total end-to-end error of an LEB position evaluation has several
independent contributions:

| Source | Typical magnitude | Affected bodies |
|--------|-------------------|-----------------|
| Chebyshev fitting residual | 1e-12 to 1e-8 AU | All |
| COB correction (SPK) | <0.001" | Jupiter-Pluto |
| COB correction (analytical fallback) | <0.01" | Jupiter-Pluto |
| Gravitational deflection omission | Up to 0.004" | **Fixed**: now included |
| Aberration formula (1st-order vs rigorous) | <0.001" | All ICRS bodies |
| Delta-T interpolation | 0.004s -> 0.002" (Moon only) | **Fixed**: use swe_deltat() |
| Precession-nutation model | <0.001" | All (IAU 2006/2000A) |
| Pipeline coordinate transform | <0.001" | All |

### Achieved Precision

After all fixes, the combined error for all 31 bodies across all three tiers:

| Category | Base (<0.001") | Medium (<0.001") | Extended |
|----------|---------------|------------------|----------|
| All planets (Sun-Pluto, Earth) | **<0.001"** | **<0.001"** | **<0.001"** |
| All asteroids (Chiron, Ceres-Vesta) | **<0.001"** | **<0.001"** | **<0.001"** |
| All ecliptic bodies (Nodes, Lilith) | **<0.001"** | **<0.001"** | <0.1" * |
| All hypothetical bodies (Uranians) | **~0.000"** | **~0.000"** | **~0.000"** |

\* Ecliptic body precision on the extended tier is limited by Meeus polynomial
degradation beyond ±20 centuries from J2000.0. Within ±1000 CE of J2000,
ecliptic body errors are <0.001".

### Error Amplification in Secondary Pipelines

When computing **heliocentric** or **equatorial** positions from the
geocentric ICRS data, errors can be amplified:

- **Moon heliocentric**: Geocentric error is amplified by the ratio
  heliocentric_distance / geocentric_distance ~ 390. A 0.0003" geocentric
  error becomes ~0.01" heliocentric.

- **Nearby asteroid latitude velocity**: The ICRS->ecliptic coordinate
  transform involves division by geocentric distance, amplifying errors
  for close-approach asteroids.

These are architectural limitations of the ICRS storage strategy and are
handled via relaxed tolerances in the secondary pipeline tests.

---

## 16. Historical Problems and Solutions

This section documents the major precision problems encountered during
LEB development and their solutions.

### Problem 1: COORD_GEO_ECLIPTIC — Retrograde Cusps (v1)

**Symptom**: 3-5 arcsecond errors at planetary stations.

**Cause**: Geocentric ecliptic coordinates have cusps at retrograde stations.
Chebyshev polynomials cannot fit discontinuities in the derivative.

**Solution**: Abandoned `COORD_GEO_ECLIPTIC` entirely. Store positions in
ICRS barycentric coordinates (always smooth) and apply coordinate transforms
at runtime.

### Problem 2: Missing Gravitational Deflection (v2)

**Symptom**: 3.95 arcsecond error for Saturn; systematic ~0.002" errors for
all planets.

**Cause**: The LEB pipeline did not include PPN gravitational deflection by
the Sun, Jupiter, and Saturn. Skyfield's `apparent()` method applies this
correction automatically.

**Solution**: Implemented `_apply_gravitational_deflection()` in `fast_calc.py`
with three deflectors (Sun, Jupiter barycenter, Saturn barycenter), using the
same PPN formula as Skyfield.

### Problem 3: COB Oscillations in Chebyshev Data (v2)

**Symptom**: 0.1-2 arcsecond errors for Jupiter-Pluto, persisting even with
short intervals and high degree.

**Cause**: Outer planet positions were stored as planet-center = barycenter +
COB correction. The COB correction contains high-frequency oscillations from
inner moons that Chebyshev polynomials cannot fit efficiently.

**Solution**: Created `COORD_ICRS_BARY_SYSTEM` (type 4). Store pure system
barycenters (smooth) in the LEB file; apply COB correction at runtime using
`planet_centers.bsp` or analytical moon theory.

### Problem 4: COB Evaluation Time Bug (v3)

**Symptom**: ~0.002" systematic error for outer planets, even with runtime COB.

**Cause**: `_apply_cob_correction()` was evaluating the COB correction at
**retarded time** (jd_tt - light_time), but Skyfield evaluates COB at
**observer time** (jd_tt).

**Solution**: One-line fix: changed `jd_tt - lt` to `jd_tt` in the COB
evaluation call.

### Problem 5: Asteroid Pipeline Mismatch (v3)

**Symptom**: 0.3-0.4 arcsecond systematic errors for all asteroids.

**Cause**: The reference asteroid pipeline used ecliptic J2000 + manual
precession + nutation, while the LEB pipeline used ICRS + Skyfield's
precession-nutation matrix. The two approaches differ at the 0.3" level.

**Solution**: Created `_SpkType21Target` VectorFunction wrapper in `spk.py`
that routes asteroids through Skyfield's observe/apparent pipeline (the same
path used by the reference code).

### Problem 6: Ecliptic Body Generation Bug (v3)

**Symptom**: Bodies 13, 21, 22 had larger-than-expected errors.

**Cause**: `generate_ecliptic_bodies_vectorized()` hardcoded `interval_days=8,
degree=13` for ALL ecliptic bodies, ignoring per-body BODY_PARAMS. Bodies
13/21/22 had been updated to `interval=4, degree=15` but the generator
still used the old values.

**Solution**: Split ecliptic bodies by their params in the generator.
Bodies with non-standard params are routed to the scalar generation path.

### Problem 7: Delta-T Discrepancy (v3)

**Symptom**: ~0.002 arcsecond Moon error, varying with epoch.

**Cause**: LEB reader's linearly-interpolated sparse Delta-T table introduced
up to ~0.004 seconds of error near 1985. The Moon moves ~13 deg/day =
0.00054 deg/sec, so a 0.004s error produces ~0.002" position error.

**Solution**: `fast_calc_ut()` now uses `swe_deltat()` (Skyfield's precise
model) instead of `reader.delta_t()` for the UT->TT conversion.
