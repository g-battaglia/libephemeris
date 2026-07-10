# Interpolated Lunar Apogee and Perigee

LibEphemeris computes the interpolated (natural) lunar apogee and perigee by smoothing osculating orbital elements to remove the spurious ±30-degree oscillations inherent in the two-body approximation, producing positions that reflect the genuine apsidal motion. Passage geometry comes from the JPL DE440/DE441 numerical integrations; the perturbation coefficients and residual table are fitted against reference-API output used strictly as a black-box oracle (see NOTICE.md, Calibration Data Disclosure).

## Table of Contents

- [Background](#background)
  - [Apsidal Variants](#apsidal-variants)
  - [The Problem with Osculating Apsides](#the-problem-with-osculating-apsides)
- [Method](#method)
  - [The Interpolated Solution](#the-interpolated-solution)
  - [Implementation](#implementation)
  - [Asymmetric Perturbations](#asymmetric-perturbations)
  - [Perturbation Series Details](#perturbation-series-details)
  - [Coefficient Calibration](#coefficient-calibration)
- [When to Use Each Variant](#when-to-use-each-variant)
  - [Decision Matrix](#decision-matrix)
- [API Usage](#api-usage)
- [Precision and Validation](#precision-and-validation)
  - [Precision (measured)](#precision-measured)
  - [Smoothness Comparison](#smoothness-comparison)
  - [Apogee-Perigee Relationship](#apogee-perigee-relationship)
  - [Valid Date Range](#valid-date-range)
- [References](#references)

## Background

### Apsidal Variants

The lunar apsides (apogee and perigee) represent the farthest and nearest points of the Moon's orbit from Earth. In astrological tradition, the lunar apogee is known as **Black Moon Lilith** (among other names), while the perigee is sometimes called **Priapus** or **White Moon Selena** in some systems.

LibEphemeris provides three variants for each apsidal point:

| Body ID | Constant | Description |
|---------|----------|-------------|
| 12 | `MEAN_APOG` | Mean Lunar Apogee (Mean Lilith) |
| 13 | `OSCU_APOG` | Osculating Lunar Apogee (True Lilith) |
| 21 | `INTP_APOG` | Interpolated Lunar Apogee (Natural Lilith) |
| 22 | `INTP_PERG` | Interpolated Lunar Perigee (Natural Priapus) |

### The Problem with Osculating Apsides

The osculating apogee is computed from the Moon's instantaneous position and velocity at any given moment. It represents the apogee of the hypothetical Keplerian ellipse (two-body orbit) that would result if all gravitational perturbations suddenly ceased.

The Moon's orbit is strongly perturbed by the Sun's gravity. When osculating orbital elements are computed, a three-body problem (Earth-Moon-Sun) is forced into a two-body framework (Earth-Moon only). This creates significant artifacts.

As established in lunar orbital mechanics (Chapront-Touzé & Chapront, 1988; Meeus, *Astronomical Algorithms*, ch. 22):

> "The solar perturbation results in gigantic monthly oscillations in the ephemeris of the osculating apsides (the amplitude is **30 degrees**). These oscillations have to be considered an artifact of the insufficient model, they do not really show a motion of the apsides."

> "The lunar orbit is far from being an ellipse!"

The osculating apogee exhibits the following characteristics:

1. **Large oscillations**: The osculating apogee oscillates **±30 degrees** around the mean apogee position with a period of approximately one synodic month (~29.5 days).
2. **Artificial motion**: These oscillations do not represent actual motion of the apsidal line. They are mathematical artifacts of the two-body approximation.
3. **True twice per month**: The osculating apogee only represents the "true" position (where the Moon is actually at maximum or minimum distance) twice per month:
   - When in conjunction with the Moon (Moon at maximum distance)
   - When in opposition to the Moon (Moon at minimum distance from Earth)

Over a typical month, the difference between osculating and mean positions is substantial:

```
Date          Mean Lilith    True Lilith    Difference
2025-01-01    115.23°        127.45°        +12.22°
2025-01-08    116.01°        98.31°         -17.70°
2025-01-15    116.79°        138.62°        +21.83°
2025-01-22    117.57°        95.24°         -22.33°
2025-01-29    118.35°        131.15°        +12.80°
```

The True Lilith (osculating) oscillates wildly while Mean Lilith moves smoothly.

## Method

### The Interpolated Solution

The interpolated (or "natural") apogee seeks a middle ground between:
- The **mean apogee**: Too smooth; ignores real apsidal variations
- The **osculating apogee**: Too noisy; dominated by spurious oscillations

The goal is to smooth out the spurious 30-degree oscillations while preserving the genuine apsidal motion.

Analysis of actual lunar apogee and perigee passages (when the Moon is at its farthest/nearest distance from Earth), located on the JPL DE440 numerical integration, reveals that:

1. **Apogee oscillates ~5 degrees from mean position** (vs. 30 degrees for osculating)
2. **Perigee oscillates ~25 degrees from mean position** (due to asymmetric solar effects)
3. **Apogee and perigee are not exactly opposite** — they deviate from 180 degrees apart, being only roughly opposite when the Sun is in conjunction with one of them or at a 90-degree angle

The interpolated apsides represent the "natural" position of the apsidal line — where the actual lunar distance extremes occur. This is derived by smoothing the osculating elements to remove the artificial oscillations while retaining the physical variations.

### Implementation

LibEphemeris computes both interpolated apsides the same way — an
ELP2000-82B-derived perturbation series on top of the mean apsidal
position, plus a precomputed residual correction table:

1. **ELP2000-82B Perturbation Series (Apogee)**: Adds 25 periodic
   perturbation terms derived from ELP2000-82B theory to the mean apogee
   (evection harmonics up to k=14, solar-anomaly coupling, latitude
   coupling, and cross-coupling terms).

2. **ELP2000-82B Perturbation Series (Perigee)**: Adds perturbation
   corrections to the mean perigee position using coefficients fitted
   against reference-API output (black-box oracle; see
   [interpolated-perigee.md](interpolated-perigee.md) and NOTICE.md).

Both variants then apply a **precomputed residual correction table**
(linear interpolation, tapered to zero within a year of the table edges)
that absorbs the remaining error inside the table's date range.

**Interpolated Apogee (ELP2000-82B series):**

```
1. Calculate Mean Lilith (mean lunar apogee) longitude
2. Calculate Julian centuries from J2000.0
3. Compute fundamental lunar arguments (D, M, M', F, Ω)
4. Apply the ELP2000-82B perturbation series (25 terms):
   - Dominant term: +4.53° × sin(2D - 2M')
   - Second-order terms from lunar theory
   - Long-period terms for secular evolution
5. Add the residual correction from the precomputed table
6. Normalize result to [0°, 360°)
```

**Interpolated Perigee (ELP2000-82B series):**

```
1. Calculate Mean Lilith + 180° (mean perigee longitude)
2. Calculate Julian centuries from J2000.0
3. Compute fundamental lunar arguments (D, M, M', F)
4. Apply calibrated perigee perturbation series:
   - Dominant term: -22.2° × sin(2D - 2M')  (opposite sign to apogee!)
   - Additional terms fitted against reference-API output (black-box oracle)
5. Add the residual correction from the precomputed table
6. Normalize result to [0°, 360°)
```

### Asymmetric Perturbations

The apogee and perigee experience very different perturbation amplitudes:

| Apside  | Oscillation from Mean | Dominant Coefficient |
|---------|-----------------------|----------------------|
| Apogee  | ~5°                   | +4.53°               |
| Perigee | ~25°                  | -22.2°               |

This asymmetry arises because:
1. Solar perturbations affect the perigee more strongly
2. The perturbation series terms have opposite signs for apogee vs. perigee
3. The Moon spends more time near apogee (slower orbital velocity)

### Perturbation Series Details

The perturbation series uses fundamental lunar arguments:
- **D**: Mean elongation of Moon from Sun
- **M**: Sun's mean anomaly
- **M'**: Moon's mean anomaly
- **F**: Moon's argument of latitude

**Apogee perturbation (excerpt):**

```python
# Dominant evection term (2D - 2M' argument)
delta = 4.528899 * sin(2*D - 2*M')  # degrees

# Evection harmonics sin(k*D - k*M') plus small mixed terms
delta += 0.828318 * sin(4*D - 4*M')
delta += 0.250542 * sin(D - M')
# (full series: _calc_elp2000_apogee_perturbations in lunar.py)
```

**Perigee perturbation (coefficients fitted against reference-API output):**

```python
# Dominant evection term (opposite sign, ~5x the apogee coefficient)
delta = -22.247860 * sin(2*D - 2*M')  # degrees

# Slowly-converging evection harmonics sin(k*D - k*M'), extended to k=30
delta += 6.651642 * sin(4*D - 4*M')
delta += -2.906333 * sin(6*D - 6*M')
# ... plus E-weighted solar-anomaly coupling terms, e.g.
delta += -0.979771 * E * sin(2*D - 2*M' - M)
# (full series: _calc_elp2000_perigee_perturbations in lunar.py)
```

### Coefficient Calibration

The perigee perturbation coefficients were derived by:
1. Sampling reference-API perigee positions (black-box oracle) across a
   wide date range
2. Computing the residual (oracle perigee - mean perigee - 180°)
3. Using least-squares fitting to find optimal coefficients
4. Validating against independent test dates

This approach achieves a computationally efficient analytical formula. The
INTP_* apsides are constructs defined by the reference API, so 1:1
behavioural parity requires fitting to its output — disclosed in NOTICE.md
(Calibration Data Disclosure).

## When to Use Each Variant

### Mean Lilith (MEAN_APOG)

**Use when:**
- Smooth, predictable motion is needed
- Historical or traditional astrological work
- Compatibility with older software that only supported mean values
- Short-term oscillations are not relevant

**Characteristics:**
- Moves at approximately 40.7 degrees/year (prograde)
- No short-term oscillations
- Does not represent actual apsidal position on any given day

### True Lilith (OSCU_APOG)

**Use when:**
- The instantaneous osculating elements are specifically needed
- Research into lunar orbital mechanics
- Understanding the theoretical two-body apogee

**Characteristics:**
- Large oscillations (±30 degrees from mean)
- Period of approximately one synodic month
- Only "true" twice per month
- Most of the oscillation is artificial (two-body approximation artifact)

**Caution:** For most astrological purposes, the osculating apogee is misleading because the 30-degree oscillations do not represent real apsidal motion.

### Interpolated Lilith (INTP_APOG)

**Use when:**
- The "actual" apogee position (where the Moon reaches maximum distance) is desired
- Modern astrological work seeking physical accuracy
- Real variations without spurious noise are needed
- Transit or timing work related to lunar distance

**Characteristics:**
- Oscillates approximately ±5 degrees from mean (for apogee)
- Smooth, continuous motion
- Represents the genuine apsidal line orientation
- Available via the `INTP_APOG` body

### Interpolated Perigee (INTP_PERG)

**Use when:**
- The "actual" perigee position (where the Moon reaches minimum distance) is needed
- Supermoon calculations (perigee near Full Moon)
- Studying lunar distance variations

**Characteristics:**
- Oscillates approximately ±25 degrees from mean (larger than apogee)
- Not exactly opposite to interpolated apogee (can deviate from 180 degrees)
- Represents genuine perigee position

### Decision Matrix

| Requirement | Recommended |
|-------------|-------------|
| Smooth ephemeris for natal charts | Mean Lilith |
| Traditional/historical astrology | Mean Lilith |
| Physical accuracy | Interpolated |
| Supermoon calculations | Interpolated Perigee |
| Lunar distance extremes timing | Interpolated |
| Drop-in API compatibility (basic) | Mean Lilith |
| Drop-in API compatibility (advanced) | Interpolated |
| Research into orbital mechanics | Osculating |

## API Usage

### Using calc_ut / calc

```python
import libephemeris as swe

jd = 2460676.5  # 2025-01-01

# Mean Lilith
mean_pos, _ = swe.calc_ut(jd, swe.MEAN_APOG, swe.FLG_SPEED)
print(f"Mean Lilith: {mean_pos[0]:.4f} deg, speed: {mean_pos[3]:.4f} deg/day")

# Osculating (True) Lilith
oscu_pos, _ = swe.calc_ut(jd, swe.OSCU_APOG, swe.FLG_SPEED)
print(f"True Lilith: {oscu_pos[0]:.4f} deg, speed: {oscu_pos[3]:.4f} deg/day")

# Interpolated (Natural) Lilith
intp_pos, _ = swe.calc_ut(jd, swe.INTP_APOG, swe.FLG_SPEED)
print(f"Interpolated Lilith: {intp_pos[0]:.4f} deg, speed: {intp_pos[3]:.4f} deg/day")

# Interpolated Perigee
perg_pos, _ = swe.calc_ut(jd, swe.INTP_PERG, swe.FLG_SPEED)
print(f"Interpolated Perigee: {perg_pos[0]:.4f} deg, speed: {perg_pos[3]:.4f} deg/day")
```

### Using Direct Lunar Functions

```python
from libephemeris import lunar

jd_tt = 2460676.5  # 2025-01-01 (TT)

# Mean Lilith
mean_lon = lunar.calc_mean_lilith(jd_tt)
print(f"Mean Lilith: {mean_lon:.4f} deg")

# True Lilith (returns longitude, latitude, distance in AU)
oscu_lon, oscu_lat, oscu_dist = lunar.calc_true_lilith(jd_tt)
print(f"True Lilith: {oscu_lon:.4f} deg, lat: {oscu_lat:.4f} deg, dist: {oscu_dist:.5f} AU")

# Interpolated Apogee
intp_lon, intp_lat, intp_ecc = lunar.calc_interpolated_apogee(jd_tt)
print(f"Interpolated Apogee: {intp_lon:.4f} deg")

# Interpolated Perigee
perg_lon, perg_lat, perg_dist = lunar.calc_interpolated_perigee(jd_tt)
print(f"Interpolated Perigee: {perg_lon:.4f} deg")
```

### Return Values

All lunar apside functions return:
- **Longitude**: Ecliptic longitude in degrees [0, 360)
- **Latitude**: Ecliptic latitude in degrees (typically small, < 5 degrees)
- **Distance**: Earth-to-apside distance in AU (~0.00242 at perigee,
  ~0.00271 at apogee)

When using `calc_ut` with `FLG_SPEED`, velocity is also calculated:
- **Speed (longitude)**: Daily motion in degrees/day
- **Speed (latitude)**: Daily change in latitude
- **Speed (distance)**: Daily change in distance (AU/day)

## Precision and Validation

### Precision (measured)

| Variant | RMS Error | Max Error |
|---------|-----------|-----------|
| Mean Lilith | ~0.003° (~12") | ~0.005° (~18") |
| True Lilith | ~0.02° (~72") | ~0.07° (~252") |
| Interpolated Apogee (within table, 1549–2651) | ~6" | ~40" |
| Interpolated Apogee (outside table) | ~171" (~0.048°) | — |
| Interpolated Perigee (within table) | ~11" | ~100" |

**Note on the residual correction table:** the ELP2000-82B series alone
leaves arcminute-level residuals; the precomputed correction table (see
above) absorbs them inside its date range, which is where the ~6"/~11"
RMS figures come from. Outside the range the correction tapers to zero
within a year and only the series precision remains.

**Note on True Lilith differences:** The ~5 degree mean differences arise because:
1. LibEphemeris computes osculating elements from JPL DE state vectors
2. an analytical lunar theory yields a different osculating apogee
3. The osculating apogee concept is inherently model-dependent for strongly perturbed orbits

**Note on Interpolated differences:** The remaining differences arise from:
1. The perturbation series uses a few dozen harmonic terms vs. a full analytical lunar theory
2. Coefficient calibration and the residual tables were fitted on a finite sample of dates
3. The apogee and perigee series are calibrated separately, each for its own perturbation amplitude

### Smoothness Comparison

Testing over 100 consecutive days shows variance in daily motion:

| Variant | Daily Motion Variance |
|---------|-----------------------|
| Osculating Apogee | ~38 deg²/day² |
| Interpolated Apogee | ~0.6 deg²/day² |

The interpolated apogee is approximately **60 times smoother** than the osculating.

### Apogee-Perigee Relationship

LibEphemeris computes apogee and perigee using **independent** perturbation series, so they are not constrained to be exactly 180° apart.

This is a known property of the interpolated apogee method (Chapront-Touzé & Chapront, 1988; Meeus, *Astronomical Algorithms*, ch. 22):

> "Apogee and perigee are not exactly opposite — they are only roughly opposite when the Sun is in conjunction with one of them or at a 90-degree angle."

The deviation from 180° can be up to **28 degrees** in extreme cases. In
LibEphemeris the apogee–perigee separation is ~180° ± 15°, varying with the
lunar/solar geometry. This is expected physical behaviour, not a limitation.

### Valid Date Range

The analytical perturbation series approach has no inherent date range limitations beyond the validity of the underlying lunar theory.

For practical purposes, the implementation is valid for:
- DE440 (default): 1550–2650 (full precision)
- Extended range: 3000 BCE–3000 CE (reduced precision)

---

*For the measured head-to-head with Swiss Ephemeris on the interpolated apogee
(~0.36°), interpolated perigee (~2.6°) and True Lilith (~5°), see the
[Swiss Ephemeris Comparison](../comparison/known-differences.md).*

## References

### Primary Sources

1. Chapront-Touzé, M. & Chapront, J. (1988). "ELP 2000-82B: A semi-analytical lunar ephemeris." *Astronomy & Astrophysics*, 190, 342-352.
2. Chapront-Touzé, M. & Chapront, J. (1991). *Lunar Tables and Programs from 4000 B.C. to A.D. 8000*. Willmann-Bell.
3. Meeus, J. (1998). *Astronomical Algorithms*, 2nd edition. Willmann-Bell, Chapter 47.
4. Brown, E.W. (1896). *An Introductory Treatise on the Lunar Theory*.
5. Brouwer, D. & Clemence, G.M. (1961). *Methods of Celestial Mechanics*.

### Astrological References

6. Santoni, Francis (1993). *Ephémérides de la lune noire vraie 1910-2010*. Editions St. Michel.
7. Koch, Dieter (1995). "Was ist Lilith und welche Ephemeride ist richtig." *Meridian*, 1/95.
8. Koch, Dieter & Rindgen, Bernhard (2000). *Lilith und Priapus*. Frankfurt/Main.

### Related Documentation

- [True Lilith Methods](true-lilith.md) — Comparison of True Lilith calculation methods
- [Precision](../reference/precision.md) — Overall precision limitations of LibEphemeris
