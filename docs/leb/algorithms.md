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
17. [Sources and Project-Choice Boundary](#17-sources-and-project-choice-boundary)

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
positions of astronomical objects at specific times. NASA JPL DE440/DE441 are
the project's public high-precision state source. They encode planetary
positions as Chebyshev polynomial segments fitted to numerical integrations of
the Solar System's equations of motion.

LEB applies the same polynomial family one level higher: its generator samples
the registered JPL/IAU project pipeline and fits a project-native segmented
representation. The sampled values remain derived from the public JPL state
and named frame/correction standards; they are not observations of an API
compatibility target. This creates a purpose-built cache of reviewed state or
pre-transformed channels.

### Provenance contract

The scientific authority and the numerical representation are deliberately
separate:

- NASA/JPL DE440/DE441 and registered public analytical models supply the
  astronomical state being represented.
- IERS/IAU/ERFA and the cited physical models supply frame and apparent-place
  transformations.
- Clenshaw (1955) supplies the series-evaluation recurrence; standard
  Chebyshev identities supply the derivative recurrence.
- LEB segment widths, degrees, channel selection, binary layout, compression
  targets, chunk size, cache policy, verification grids, and failure behavior
  are project-authored choices disclosed here and in the generators.

Neither a segment parameter nor a retained coefficient is selected by fitting
the output of another ephemeris implementation. Compatibility comparisons are
regression checks only and are not generator inputs.

### Why Chebyshev Polynomials?

Chebyshev polynomials are a strong practical choice for smooth finite-interval
approximation because:

- Chebyshev-node interpolation often has near-minimax maximum-error behavior;
  it is **not** automatically the exact minimax polynomial (which would require
  a method such as Remez exchange).
- Endpoint-clustered nodes strongly suppress the Runge oscillation associated
  with high-degree interpolation on equally spaced nodes.
- They can be evaluated efficiently via the Clenshaw algorithm
- Their derivatives have a simple closed-form recurrence
- Piecewise Chebyshev representations are also documented by JPL's public SPK
  specification, although LEB's layout and fitted channels are project-native.

The key insight is that the registered barycentric state channels are smooth
over short intervals (with exceptions noted below). Degree and interval cannot
be assigned a universal angular accuracy: the reviewed artifact uses
body-specific values and accepts them only after native-component and
end-to-end angular validation.

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

where the coefficients `c_k` depend on the chosen approximation procedure. In
this project, `chebfit` solves the Chebyshev-basis system for values sampled at
`N + 1` Chebyshev nodes. That square-node fit should not be described as a
proof of the globally minimax polynomial. The series is analogous to a Fourier
series but uses a polynomial basis on a finite interval.

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

The selected Chebyshev construction mitigates those problems because:
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

This requires only two temporary variables and N multiply-add operations; no
arrays are allocated. Historical profiling on one development machine measured
about 0.5 microseconds for degree 13, but that number is workload, interpreter,
and hardware dependent rather than a format guarantee.

**Implementation** (`libephemeris.leb_reader._clenshaw`):

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

For single-point evaluation (the common case in `calc_ut()`), numpy array
creation overhead (~5 us) would dominate the actual Clenshaw loop (~0.5 us for
degree 13). Pure Python `float` operations avoid this overhead entirely.

For batch evaluation during generation, numpy's vectorized `chebval()` is used
instead.

---

## 4. Chebyshev Fitting (Generation)

### Chebyshev Nodes

The generation process samples the registered source function at **Chebyshev nodes**
(also called Chebyshev-Gauss points or Type I nodes):

```
x_k = cos(pi * (k + 0.5) / n)    for k = 0, 1, ..., n-1
```

where n = degree + 1 (the number of nodes equals the number of coefficients).

These nodes are not uniformly spaced -- they cluster near the boundaries of
[-1, 1]. This clustering suppresses the severe endpoint oscillation associated
with equally spaced high-degree interpolation. It supports good practical
maximum-error behavior for smooth channels, but the artifact's validation—not
the node name—decides whether a selected degree and interval are acceptable.

**Mapped to Julian Days:**

```python
jd_nodes = 0.5 * (jd_end - jd_start) * chebyshev_nodes + 0.5 * (jd_start + jd_end)
```

### Fitting Procedure

For each segment, the fitting process is:

1. **Sample**: Evaluate the registered source channel (direct JPL/Skyfield,
   independently sourced analytical model, or public SPK kernel) at the
   n = degree + 1 Chebyshev nodes.
2. **Fit**: Compute coefficients using `numpy.polynomial.chebyshev.chebfit()`.
   With `degree + 1` distinct nodes and the same number of coefficients, this
   is an exactly determined Chebyshev-basis fit evaluated by NumPy's
   least-squares machinery.
3. **Verify**: Evaluate the fitted polynomial at 10 uniformly-spaced test
   points (NOT on the Chebyshev nodes) and compare against the registered
   source channel. Track the maximum sampled error.

The ten interior points are a deterministic *segment screening grid*. They
catch many wrap, boundary, and insufficient-degree failures, but ten samples
cannot prove a continuous supremum-error bound. Release generation therefore
adds independent post-write sampling and end-to-end angular checks described
in [Section 15](#15-error-analysis-and-precision-budget).

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

This eliminates repeated per-call overhead (time conversion, SPK segment
lookup, and Python dispatch). Historical development measurements observed
roughly 150x faster generation for some planetary batches; that figure is a
hardware/workload measurement, not an algorithmic guarantee.

### Segment Width Selection

The segment width (interval_days) controls the trade-off between file size
and approximation accuracy:

| Shorter intervals | Longer intervals |
|-------------------|------------------|
| Better accuracy | Worse accuracy |
| More segments | Fewer segments |
| Larger file | Smaller file |

The selected width depends on how rapidly the stored channel changes. The
values below are current project parameters validated for the named channel,
not optima proved by a publication:

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

**Implementation** (`libephemeris.leb_reader._deriv_coeffs`):

```python
def _deriv_coeffs(coeffs: tuple[float, ...]) -> tuple[float, ...]:
    n = len(coeffs) - 1
    if n == 0:
        return (0.0,)
    if n == 1:
        return (coeffs[1],)
    d = [0.0] * n
    d[n - 1] = 2.0 * n * coeffs[n]
    for k in range(n - 2, 0, -1):
        d_k2 = d[k + 2] if k + 2 < n else 0.0
        d[k] = d_k2 + 2.0 * (k + 1) * coeffs[k + 1]
    d_2 = d[2] if n >= 3 else 0.0
    d[0] = d_2 / 2.0 + coeffs[1]
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

This required **two additional pipeline evaluations** for a vector velocity
request. The analytical derivative:

- reduces the required pipeline evaluations from three to one (the wall-clock
  speedup depends on body, cache state, degree, interpreter, and backend);
- avoids the selected finite-difference step and its truncation/cancellation
  trade-off; and
- differentiates the represented polynomial consistently.

Differentiation can amplify errors in high-order coefficients, so the
derivative is not declared more accurate merely because it is analytical.
Velocity precision requires its own source-state comparison and end-to-end
tests.

---

## 6. Coordinate Systems and Storage Strategies

### The Five Coordinate Types

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

This separation removes the high-frequency center offset from the stored
system-barycenter polynomial. Artifact-scoped measurements, rather than the
coordinate-type definition alone, establish the resulting precision. See
[Section 12](#12-center-of-body-cob-corrections) for full details.

#### COORD_ECLIPTIC (type 1) — Ecliptic of Date

**Used for:** Mean Node, True Node, Mean Apogee (Lilith), Osculating Apogee,
Interpolated Apogee, Interpolated Perigee

Stores (longitude, latitude, distance) in degrees/degrees/AU in the ecliptic
coordinate system of the date. These bodies are computed by analytical
formulas that directly output ecliptic coordinates; storing them in ICRS
would require an unnecessary inverse rotation.

#### Type 2 — Retired (was Heliocentric Ecliptic)

Value 2 carried the heliocentric-ecliptic channels of the pre-3.1.0
`uranians` companion. The companion and its evaluation pipeline are retired:
fictitious bodies are always served by their runtime analytical models, and
no body may be materialized with this coordinate type. The value stays
reserved so legacy files in the wild keep an unambiguous meaning.

#### COORD_GEO_ECLIPTIC (type 3) — Reserved

Defined in the format but **not used by any body**. It was originally planned
for pre-transformed geocentric ecliptic positions. The design was abandoned
because it couples the stored channel to observer, apparent-place, flag, and
longitude-unwrapping conventions and was less compact on the project's test
grid. A physical retrograde station is a smooth turning point, not a
discontinuity; the difficulty was the angular representation's curvature and
wrap handling, not an inability of polynomials to represent zero velocity.

### Why Not Store Everything in Ecliptic?

Early versions of LEB attempted to store planet positions in geocentric
ecliptic coordinates (`COORD_GEO_ECLIPTIC`). This was rejected because:

1. **Angular-channel conditioning**: Around a station, apparent longitude has
   a smooth turning point but can have substantially greater curvature than the
   underlying barycentric Cartesian state. A fixed degree/interval therefore
   used the coefficient budget less efficiently on the recorded test grid; the
   observed fitting error at stations was 3--5 arcseconds.

2. **Coordinate conventions**: Longitude requires explicit unwrapping, and a
   pre-transformed channel freezes observer/correction choices that callers may
   legitimately change with flags. Latitude sign changes are continuous and
   are not, by themselves, a mathematical discontinuity.

The ICRS barycentric frame keeps the stored state smooth and reusable. The
coordinate transforms that produce the apparent ecliptic position are applied
at runtime, after the Cartesian state has been recovered from the Chebyshev
fit.

---

## 7. The ICRS Pipeline: From Barycentric to Apparent

The ICRS pipeline (Pipeline A/A') transforms raw barycentric positions into
apparent geocentric coordinates. It implements the registered public
apparent-place chain and is validated independently against the direct
JPL/Skyfield path.

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

# 3 fixed-point iterations (project choice; validated over supported channels)
for _ in range(3):
    dist = sqrt(geo[0]**2 + geo[1]**2 + geo[2]**2)
    lt = dist / C_LIGHT_AU_DAY
    retarded_pos, _ = reader.eval_body(ipl, jd_tt - lt)
    geo = retarded_pos - earth_pos
```

Three iterations are the project's fixed runtime choice because:

- The speed of planets is ~1e-4 c (much less than light speed)
- Each iteration reduces the error by a factor of ~v/c ~ 1e-4
- The contraction estimate after three updates is of order `(v/c)^3` times
  the initial correction for ordinary Solar-System targets.

That estimate motivates the choice but is not treated as a universal error
certificate. Supported bodies/date ranges are checked against the registered
JPL-backed pipeline; a future faster or close-approach channel must revalidate
the iteration count.

### Light-Time for System Barycenters

For bodies stored as `COORD_ICRS_BARY_SYSTEM` (Jupiter-Pluto), runtime first
resolves the physical planet center from the matching JPL center segment when
that segment covers the requested epoch. Light time is then iterated on the
observer-to-body-center vector, including for `FLG_HELCTR` and `FLG_BARYCTR`.

If no center segment covers the epoch, the fallback is explicit: the stored
system barycenter is used as the target. No analytical satellite theory is
substituted for a missing JPL center segment.

---

## 9. Gravitational Deflection

### Physical Basis

General relativity predicts that massive bodies deflect light rays passing
near them. The Sun deflects a distant source by about 1.75 arcseconds at the
solar limb. Representative limiting scales are:

| Deflector | Maximum deflection |
|-----------|-------------------|
| Sun | ~1.75" at the limb for a distant source; finite-source geometry varies |
| Jupiter | ~0.017" (for bodies near Jupiter's limb) |
| Saturn | ~0.006" |

These effects matter in any sub-arcsecond error budget and therefore cannot be
silently omitted.

### Finite-source PPN reduction

The runtime follows the finite-source point-mass reduction used by the
MIT-licensed Skyfield relativity module and the cited IERS apparent-place
foundation. For each deflector it constructs unit vectors from the deflector
to observer and target, evaluates the deflector at the clamped closest-approach
epoch, and applies the explicit dot-product form implemented in
`fast_calc._apply_gravitational_deflection`. General relativity fixes the PPN
parameter to `gamma = 1`; public gravitational parameters and the exact SI
speed of light set the scale. The near-deflector guard (`0.01 AU`) is a
project safety policy for observer-at-deflector geometries, not a published
physical radius.

### Implementation in LEB

**Deflectors** (`libephemeris.fast_calc._DEFLECTORS`):

```python
_DEFLECTORS = (
    (SUN,  1.0),       # Sun, mass ratio = 1.0 (reference)
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
- Skipped for heliocentric, barycentric, true-position, and `FLG_NOGDEFL`
  modes. `FLG_NOABERR` suppresses aberration only and does not suppress
  deflection.

**Historical development measurement:** on the recorded Saturn stress case,
omitting gravitational deflection produced a 3.95" residual. Restoring the
documented Sun/Jupiter/Saturn reduction brought the reviewed finite validation
grid below its 0.001" acceptance threshold; this is not a universal bound for
arbitrary custom artifacts.

---

## 10. Stellar Aberration

### Physical Basis

Stellar aberration is a relativistic effect caused by the observer's velocity
relative to the incoming light. It shifts the apparent position of a body
in the direction of the observer's motion. The maximum annual aberration
(due to Earth's orbital velocity of ~30 km/s) is about 20.5 arcseconds.

### Relativistic formula and compatibility fallback

Normal runtime calls carry a positive light time and use the full
special-relativistic vector formula implemented by the MIT-licensed Skyfield
`add_aberration()` routine, including the inverse Lorentz factor. The local
implementation is written explicitly in `fast_calc._apply_aberration` and is
validated against the registered JPL/Skyfield channel.

Only an internal caller that omits light time takes the retained first-order
Bradley fallback:

```
u = geo / |geo|
beta = observer_velocity / c
u_first_order = normalize(u + beta - u * dot(u, beta))
```

That branch is documented for backward compatibility; it is not the primary
LEB apparent-place reduction and no compatibility-target implementation is its
scientific authority.

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
# Vondrák 2011 long-term precession (via pyerfa's ltp/ltpb) combined with
# the nutation angles — valid +/-200,000 years, matching the apparent-place
# reduction used everywhere else in the library.
from libephemeris.precession_vondrak import vondrak_pn_matrix

pn_matrix, eps_true = vondrak_pn_matrix(jd_tt, dpsi, deps)
```

On the LEB path the nutation angles (dpsi, deps) come from the LEB-stored
Chebyshev coefficients; on the Skyfield path they come from Skyfield's
IAU 2000A model. Both paths share the single Vondrák matrix source.

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

The active mean obliquity is the angle between the Vondrák (2011) mean
ecliptic and equator poles, evaluated by the shared long-term implementation
in `sidereal_longterm.mean_obliquity_rad`. The same realization is used for
body reductions and house geometry. The IAU 2006 polynomial remains only a
modern-era helper where explicitly named; it is not the LEB path's general
of-date obliquity over the extended tier.

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
It is computed as a public JPL planet-center state relative to the matching
system barycenter. Its amplitude and spectrum depend on the satellite-system
solution and epoch; exact bounds belong to the tier kernel and its verification
report, not to a hard-coded universal estimate. Although tiny beside a
planet's heliocentric distance, it is material at a milliarcsecond error budget.

### Why COB Breaks Chebyshev Fitting

The COB correction contains satellite-period oscillations much faster than the
smooth system-barycenter trajectory. Project sweeps found that fitting
`planet_center = barycenter + COB` at the long planetary segment widths either
required materially larger artifacts or failed the angular storage budget.
Separating the public center segment from the stored barycenter avoids making a
particular unpublished fit the physical source.

### The Solution: Store Barycenters, Apply COB at Runtime

The `COORD_ICRS_BARY_SYSTEM` storage strategy separates the smooth and
oscillatory components:

1. **Generator**: Stores the pure system barycenter position (smooth, easy
   to fit with 32-64 day intervals and degree 13)
2. **Runtime**: Uses the pinned tier-specific planet-center SPK segment when
   it covers the retarded epoch; otherwise it retains the stored system
   barycenter unchanged.

### COB Evaluation Timing

**Critical**: The planet-center offset is evaluated at the **retarded emission
time** (`jd_tt - light_time`) during each light-time iteration. This keeps the
observer-to-target vector tied to the physical center whenever a segment is
available:

```python
# Runtime approach (simplified):
def _observe_from_bcrs(observer_pos, observer_t):
    # Light-time iteration on the selected center/barycenter target
    for _ in range(3):  # project iteration count, independently validated
        target = center_if_covered_else_barycenter(t_retarded)
        dist = |target - observer_pos|
        lt = dist / c
        t_retarded = observer_t - lt

    return target - observer_pos
```

At coverage boundaries the position and speed stencil source-lock to the
system barycenter instead of mixing center and barycenter samples.

### SPK Planet Center Segments

High-precision COB data comes from `planet_centers_{tier}.bsp` files for NAIF
planet-center targets 599, 699, 799, 899, and 999 relative to barycenters
5, 6, 7, 8, and 9. Coverage is target- and tier-specific. The authoritative
intervals are the descriptors in the selected SPK, inspected at load time and
verified by `scripts/generate_planet_centers_spk.py`; prose dates are avoided
here because regenerated public kernels can legitimately change them.

When the SPK segment does not cover the requested date, runtime explicitly uses
the stored system barycenter. There is no analytical COB fallback.

The generator copies every matching DAF descriptor. JPL partitions several
center targets into consecutive segments; retaining only the first descriptor
caused the former 1997/2014/1989 cutoffs. Near any remaining physical coverage
boundary, both samples of an `FLG_SPEED` stencil use the same system-barycenter
fallback so the finite center offset cannot become a sign-changing velocity
spike.

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
by evaluating at 10 intermediate test points and comparing with the registered
source channel after re-wrapping. This catches cases where unwrapping fails
(e.g., multiple 360-degree jumps within a single segment).

---

## 14. Delta-T (TT - UT1)

### What is Delta-T?

Delta-T is the difference between Terrestrial Time (TT, a uniform time scale
based on atomic clocks) and Universal Time (UT1, which tracks Earth's
irregular rotation). Since UT1 depends on that irregular rotation, Delta-T is
assembled from observations for the measured era, historical reconstructions
where direct measurements do not exist, and an explicitly identified
extrapolation policy. It is not a timeless analytically predictable constant.

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

On the historical development probe that motivated the runtime change, 30-day
linear interpolation differed by about 0.004 seconds near 1985. That is a
sampled observation, not a universal interpolation bound.

### Delta-T Usage in fast_calc

**Critical**: The `fast_calc_ut()` entry point does **not** use
`reader.delta_t()` for the UT->TT conversion. Instead, it uses
`deltat()` from `time_utils.py`, which provides the same high-precision
Delta-T model used by the direct JPL/Skyfield source pipeline. This keeps both
project backends on one documented time realization without claiming that a
finite-precision conversion is mathematically exact:

```python
# fast_calc.py, fast_calc_ut():
from .time_utils import deltat
delta_t = deltat(tjd_ut)  # high-precision model
jd_tt = tjd_ut + delta_t
```

`fast_calc_tt()` also uses the canonical `deltat()` for the reverse
(TT->UT) conversion (one fixed-point inversion step); the reader's
`delta_t()` method survives only as a last-resort fallback when
`deltat()` is unavailable.

---

## 15. Error Analysis and Precision Budget

### Fitting, compression, and pipeline validation are different claims

LEB1 fitting and LEB2 compression have separate error mechanisms:

1. `generate_leb.py` samples `degree + 1` Chebyshev nodes and screens every
   segment at ten independent interior points. Its optional post-generation
   verifier defaults to 500 dates per body and compares the serialized reader
   against the registered state source.
2. `leb_compression.py` truncates each coefficient to a selected binary64
   precision. If every coefficient error is bounded by `e_k`, then on
   `x in [-1, 1]`, where `|T_k(x)| <= 1`, the native-component evaluation error
   is conservatively bounded by `sum(e_k)`. The LEB2 verifier uses the still
   more conservative `(degree + 1) * per_coefficient_target` when all orders
   share one target.
3. Reordering, byte shuffle, and Zstandard are reversible; they add no numeric
   error after mantissa truncation. Corruption is a separate integrity failure,
   not an accuracy tolerance.
4. Native-component bounds do not by themselves prove apparent angular
   accuracy: observer subtraction, small geocentric distance, differentiation,
   light time, and frame transformations can amplify them. The focused
   `scripts/test_leb2_precision.py` path therefore checks calculated angular
   results, and release attestations record the exact artifact and test scope.

The 10/500/200 sample counts (the last is the LEB2 verifier default), random
seed, coefficient targets, and acceptance thresholds are project choices.
They are evidence grids, not public astronomical constants and not a proof over
unsampled continuous time.

### Error Sources

The total end-to-end error of an LEB position evaluation has several
independent contributions:

| Source | Typical magnitude | Affected bodies |
|--------|-------------------|-----------------|
| Chebyshev fitting residual | 1e-12 to 1e-8 AU | All |
| COB correction (SPK) | <0.001" | Jupiter-Pluto |
| Missing center segment | Explicit system-barycenter fallback | Jupiter-Pluto |
| Gravitational deflection omission | Up to 0.004" | **Fixed**: now included |
| Aberration formula (1st-order vs rigorous) | <0.001" | All ICRS bodies |
| Delta-T interpolation | 0.004s -> 0.002" (Moon only) | **Fixed**: use deltat() |
| Precession-nutation model | <0.001" | All (IAU 2006/2000A) |
| Pipeline coordinate transform | <0.001" | All |

### Achieved Precision

The following table records development measurements associated with the
reviewed generators; it is not a timeless guarantee for an arbitrary custom
file. The authoritative claim for a distributed artifact is its hash-scoped
[build attestation](base-core-provenance.md), together with the exact command,
date range, channels, and tolerances used to verify it.

| Category | Base (<0.001") | Medium (<0.001") | Extended |
|----------|---------------|------------------|----------|
| All planets (Sun-Pluto, Earth) | **<0.001"** | **<0.001"** | **<0.001"** |
| All asteroids (Chiron, Ceres-Vesta) | **<0.001"** | **<0.001"** | **<0.001"** |
| All ecliptic bodies (Nodes, Lilith) | **<0.001"** | **<0.001"** | <0.1" * |

\* Ecliptic body precision on the extended tier is limited by Meeus polynomial
degradation beyond ±20 centuries from J2000.0. Within ±1000 CE of J2000,
ecliptic body errors are <0.001".

### Error Amplification in Secondary Pipelines

When computing **heliocentric** or **equatorial** positions from the
geocentric ICRS data, errors can be amplified:

- **Moon heliocentric**: A Cartesian error can be amplified by the ratio
  heliocentric distance / geocentric distance, roughly 390 near the ordinary
  Earth-Moon distance. A 0.0003" angular-equivalent scale can therefore become
  roughly 0.12", although the exact amplification depends on vector direction.

- **Nearby asteroid latitude velocity**: The ICRS->ecliptic coordinate
  transform involves division by geocentric distance, amplifying errors
  for close-approach asteroids.

These are architectural limitations of the ICRS storage strategy and are
handled via relaxed tolerances in the secondary pipeline tests.

---

## 16. Historical Problems and Solutions

This section documents the major precision problems encountered during
LEB development and their solutions.

### Problem 1: COORD_GEO_ECLIPTIC conditioning (v1)

**Symptom**: 3-5 arcsecond errors at planetary stations.

**Cause**: The stored angular channel combined longitude unwrapping and the
apparent observer geometry with relatively high curvature around stations. A
station is mathematically smooth; the failure was representation efficiency at
the selected fixed degree and interval, not a derivative discontinuity.

**Solution**: Abandoned `COORD_GEO_ECLIPTIC` entirely. Store positions in
ICRS barycentric Cartesian coordinates and apply coordinate transforms at
runtime.

### Problem 2: Missing Gravitational Deflection (v2)

**Symptom**: 3.95 arcsecond error for Saturn; systematic ~0.002" errors for
all planets.

**Cause**: The LEB pipeline did not include PPN gravitational deflection by
the Sun, Jupiter, and Saturn. Skyfield's `apparent()` method applies this
correction automatically.

**Solution**: Implemented `_apply_gravitational_deflection()` in `fast_calc.py`
with three deflectors (Sun, Jupiter barycenter, Saturn barycenter), from the
documented PPN/IAU apparent-place geometry. Skyfield comparisons validate the
independent implementation but do not supply its source expression.

### Problem 3: COB Oscillations in Chebyshev Data (v2)

**Symptom**: 0.1-2 arcsecond errors for Jupiter-Pluto, persisting even with
short intervals and high degree.

**Cause**: Outer planet positions were stored as planet-center = barycenter +
COB correction. The COB correction contains high-frequency oscillations from
inner moons that Chebyshev polynomials cannot fit efficiently.

**Solution**: Created `COORD_ICRS_BARY_SYSTEM` (type 4). Store pure system
barycenters (smooth) in the LEB file; use a pinned JPL planet-center segment at
runtime when it covers the requested epoch, otherwise retain the barycenter.

### Problem 4: Center-source timing and edge continuity

**Symptom**: center/barycenter mismatch in light time and velocity spikes at
coverage edges.

**Cause**: evaluating the offset at observer time does not describe the body at
the emission epoch, while sampling opposite sides of a coverage boundary from
different sources turns a finite center offset into an artificial speed.

**Solution**: evaluate the selected physical center at each retarded epoch and
source-lock both samples of the velocity stencil; use the system barycenter on
both sides when the stencil crosses an outer edge.

### Problem 5: Asteroid Pipeline Mismatch (v3)

**Symptom**: 0.3-0.4 arcsecond systematic errors for all asteroids.

**Cause**: an earlier project asteroid path used ecliptic J2000 plus separate
precession/nutation, while the LEB path used the registered ICRS/Skyfield
reduction. The two project paths differed at the 0.3" level on the recorded
development grid.

**Solution**: Created `_SpkType21Target` VectorFunction wrapper in `spk.py`
that routes asteroids through the project's registered Skyfield
observe/apparent pipeline, keeping both project backends on the same named
JPL/IAU reduction chain.

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

**Solution**: `fast_calc_ut()` now uses the project's canonical, separately
source-mapped `deltat()` model instead of `reader.delta_t()` for the UT->TT
conversion.

---

## 17. Sources and Project-Choice Boundary

- Park, Folkner, Williams, and Boggs (2021), “The JPL Planetary and Lunar
  Ephemerides DE440 and DE441,”
  [AJ 161:105](https://doi.org/10.3847/1538-3881/abd414): astronomical source
  states and DE440/DE441 context.
- NASA/JPL NAIF,
  [SPK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html):
  public segmented ephemeris types, frames, centers, and Chebyshev storage
  conventions. LEB is not an SPK clone and has its own format.
- Clenshaw (1955), “A note on the summation of Chebyshev series,”
  [MTAC 9, 118--120](https://doi.org/10.1090/S0025-5718-1955-0071856-0):
  stable series-evaluation recurrence.
- Petit and Luzum (eds.),
  [IERS Conventions (2010), TN 36](https://iers-conventions.obspm.fr/content/tn36.pdf),
  plus ERFA/IAU standards: time, frame, precession-nutation, and apparent-place
  foundations used by the represented project pipeline.
- IEEE 754 binary64 representation and the
  [Zstandard format](https://github.com/facebook/zstd/blob/dev/doc/zstd_compression_format.md):
  bit representation and final lossless byte compression.

Everything else specific to LEB -- file headers, body/group IDs, segment
parameters, coefficient layout, truncation policy, chunking, lazy mapping,
caches, verification grids, thresholds, and fallbacks -- is project-authored
and must remain documented and tested as such. The sources above explain the
scientific inputs and published numerical building blocks; they do not endorse
or prescribe this project's engineering parameters.
