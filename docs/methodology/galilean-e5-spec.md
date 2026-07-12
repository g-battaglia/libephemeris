# Functional Specification — Galilean Satellite Module (Lieske E5)

Functional specification of `libephemeris/moon_theories/galilean.py`.
**Sources of all mathematics and coefficients below:**

- Lieske, J.H. (1998). "Galilean satellite ephemerides E5." *Astronomy &
  Astrophysics Supplement Series* 129, 205–217.
- Meeus, J. (1998). *Astronomical Algorithms*, 2nd ed., Willmann-Bell,
  Chapter 44 ("Positions of the satellites of Jupiter"), which prints the
  complete E5-based algorithm and the coefficient tables reproduced here.

This document restates the published theory plus the project-specific frame
adaptations as behavioral requirements. It deliberately contains **no
identifier names, comments, or code structure** from any prior
implementation. The implementer must work from this document alone (see
§9 "Process rules").

---

## 1. Purpose and I/O contract

Module path (fixed, required by the wheel audit):
`libephemeris/moon_theories/galilean.py`

Public API — exactly one function (this name is the project's own public
API and must be kept):

```
galilean_moon_positions(jd: float) -> tuple[
    tuple[float, float, float],
    tuple[float, float, float],
    tuple[float, float, float],
    tuple[float, float, float],
]
```

- **Input** `jd`: Julian Date, interpreted directly as JDE (TT/TDB). Do
  **not** apply any TT→TDB conversion. Do **not** apply light-time — the
  caller handles light-time upstream.
- **Output**: positions of Io, Europa, Ganymede, Callisto **in that
  order**, each as `(x, y, z)` in **kilometres**, **Jupiter-centric**, in
  the **J2000 ecliptic** frame (see §6 for the exact frame chain).
- The function must be a **pure, deterministic, smooth** function of `jd`
  (the caller central-differences it at ±1 s to obtain velocities).
- Return plain Python floats (stdlib `math` only; no numpy).
- The only name imported by the rest of the codebase is
  `galilean_moon_positions`; all other module contents are free choices.
- Inherent physical accuracy of the published theory in this form:
  ~100–200 km per moon (≈0.05″ at Jupiter's geocentric distance). The
  module is consumed for Jupiter's center-of-body offset, where this error
  is diluted by ~2×10⁻⁴.

## 2. Time variable

All time-dependent quantities below are functions of

```
t = jd − 2443000.5          (days from the E5 reference epoch)
```

## 3. Fundamental arguments (degrees)

Mean longitudes of the satellites:

| Symbol | Expression |
|---|---|
| ℓ₁ | 106.07719 + 203.488955790 · t |
| ℓ₂ | 175.73161 + 101.374724735 · t |
| ℓ₃ | 120.55883 + 50.317609207 · t |
| ℓ₄ | 84.44459 + 21.571071177 · t |

Longitudes of the proper perijoves:

| Symbol | Expression |
|---|---|
| π₁ | 97.0881 + 0.16138586 · t |
| π₂ | 154.8663 + 0.04726307 · t |
| π₃ | 188.1840 + 0.00712734 · t |
| π₄ | 335.2868 + 0.00184000 · t |

Longitudes of the proper nodes on Jupiter's equatorial plane:

| Symbol | Expression |
|---|---|
| ω₁ | 312.3346 − 0.13279386 · t |
| ω₂ | 100.4411 − 0.03263064 · t |
| ω₃ | 119.1942 − 0.00717703 · t |
| ω₄ | 322.6186 − 0.00175934 · t |

Principal inequality in Jupiter's longitude:

```
Γ = 0.33033 · sin(163.679 + 0.0010512 · t)
  + 0.03439 · sin(34.486 − 0.0161731 · t)
```

Phase of free libration:

```
Φλ = 199.6766 + 0.17379190 · t
```

Longitude of the node of Jupiter's equator on the ecliptic:

```
ψ = 316.5182 − 0.00000208 · t
```

Mean anomalies of Jupiter and Saturn:

```
G  = 30.23756 + 0.0830925701 · t + Γ
G′ = 31.97853 + 0.0334597339 · t
```

Longitude of Jupiter's perihelion (constant):

```
Π = 13.469942
```

All sine/cosine arguments in this document are in **degrees**.

## 4. Periodic series

### 4.1 Longitude corrections Σ₁…Σ₄ (degrees)

Each Σᵢ is the sum of the listed terms, `amplitude · sin(argument)`,
accumulated **in table order**.

**Σ₁ — Io (21 terms):**

| # | Amplitude | Argument |
|---|---|---|
| 1 | +0.47259 | 2(ℓ₁ − ℓ₂) |
| 2 | −0.03478 | π₃ − π₄ |
| 3 | +0.01081 | ℓ₂ − 2ℓ₃ + π₃ |
| 4 | +0.00738 | Φλ |
| 5 | +0.00713 | ℓ₂ − 2ℓ₃ + π₂ |
| 6 | −0.00674 | π₁ + π₃ − 2Π − 2G |
| 7 | +0.00666 | ℓ₂ − 2ℓ₃ + π₄ |
| 8 | +0.00445 | ℓ₁ − π₃ |
| 9 | −0.00354 | ℓ₁ − ℓ₂ |
| 10 | −0.00317 | 2ψ − 2Π |
| 11 | +0.00265 | ℓ₁ − π₄ |
| 12 | −0.00186 | G |
| 13 | +0.00162 | π₂ − π₃ |
| 14 | +0.00158 | 4(ℓ₁ − ℓ₂) |
| 15 | −0.00155 | ℓ₁ − ℓ₃ |
| 16 | −0.00138 | ψ + ω₃ − 2Π − 2G |
| 17 | −0.00115 | 2(ℓ₁ − 2ℓ₂ + ω₂) |
| 18 | +0.00089 | π₂ − π₄ |
| 19 | +0.00085 | ℓ₁ + π₃ − 2Π − 2G |
| 20 | +0.00083 | ω₂ − ω₃ |
| 21 | +0.00053 | ψ − ω₂ |

**Σ₂ — Europa (37 terms):**

| # | Amplitude | Argument |
|---|---|---|
| 1 | +1.06476 | 2(ℓ₂ − ℓ₃) |
| 2 | +0.04256 | ℓ₁ − 2ℓ₂ + π₃ |
| 3 | +0.03581 | ℓ₂ − π₃ |
| 4 | +0.02395 | ℓ₁ − 2ℓ₂ + π₄ |
| 5 | +0.01984 | ℓ₂ − π₄ |
| 6 | −0.01778 | Φλ |
| 7 | +0.01654 | ℓ₂ − π₂ |
| 8 | +0.01334 | ℓ₂ − 2ℓ₃ + π₂ |
| 9 | +0.01294 | π₃ − π₄ |
| 10 | −0.01142 | ℓ₂ − ℓ₃ |
| 11 | −0.01057 | G |
| 12 | −0.00775 | 2(ψ − Π) |
| 13 | +0.00524 | 2(ℓ₁ − ℓ₂) |
| 14 | −0.00460 | ℓ₁ − ℓ₃ |
| 15 | +0.00316 | ψ − 2G + ω₃ − 2Π |
| 16 | −0.00203 | π₁ + π₃ − 2Π − 2G |
| 17 | +0.00146 | ψ − ω₃ |
| 18 | −0.00145 | 2G |
| 19 | +0.00125 | ψ − ω₄ |
| 20 | −0.00115 | ℓ₁ − 2ℓ₃ + π₃ |
| 21 | −0.00094 | 2(ℓ₂ − ω₂) |
| 22 | +0.00086 | 2(ℓ₁ − 2ℓ₂ + ω₂) |
| 23 | −0.00086 | 5G′ − 2G + 52.225 |
| 24 | −0.00078 | ℓ₂ − ℓ₄ |
| 25 | −0.00064 | 3ℓ₃ − 7ℓ₄ + 4π₄ |
| 26 | +0.00064 | π₁ − π₄ |
| 27 | −0.00063 | ℓ₁ − 2ℓ₃ + π₄ |
| 28 | +0.00058 | ω₃ − ω₄ |
| 29 | +0.00056 | 2(ψ − Π − G) |
| 30 | +0.00056 | 2(ℓ₂ − ℓ₄) |
| 31 | +0.00055 | 2(ℓ₁ − ℓ₃) |
| 32 | +0.00052 | 3ℓ₃ − 7ℓ₄ + π₃ + 3π₄ |
| 33 | −0.00043 | ℓ₁ − π₃ |
| 34 | +0.00041 | 5(ℓ₂ − ℓ₃) |
| 35 | +0.00041 | π₄ − Π |
| 36 | +0.00032 | ω₂ − ω₃ |
| 37 | +0.00032 | 2(ℓ₃ − G − Π) |

**Σ₃ — Ganymede (43 terms):**

| # | Amplitude | Argument |
|---|---|---|
| 1 | +0.16490 | ℓ₃ − π₃ |
| 2 | +0.09081 | ℓ₃ − π₄ |
| 3 | −0.06907 | ℓ₂ − ℓ₃ |
| 4 | +0.03784 | π₃ − π₄ |
| 5 | +0.01846 | 2(ℓ₃ − ℓ₄) |
| 6 | −0.01340 | G |
| 7 | −0.01014 | 2(ψ − Π) |
| 8 | +0.00704 | ℓ₂ − 2ℓ₃ + π₃ |
| 9 | −0.00620 | ℓ₂ − 2ℓ₃ + π₂ |
| 10 | −0.00541 | ℓ₃ − ℓ₄ |
| 11 | +0.00381 | ℓ₂ − 2ℓ₃ + π₄ |
| 12 | +0.00235 | ψ − ω₃ |
| 13 | +0.00198 | ψ − ω₄ |
| 14 | +0.00176 | Φλ |
| 15 | +0.00130 | 3(ℓ₃ − ℓ₄) |
| 16 | +0.00125 | ℓ₁ − ℓ₃ |
| 17 | −0.00119 | 5G′ − 2G + 52.225 |
| 18 | +0.00109 | ℓ₁ − ℓ₂ |
| 19 | −0.00100 | 3ℓ₃ − 7ℓ₄ + 4π₄ |
| 20 | +0.00091 | ω₃ − ω₄ |
| 21 | +0.00080 | 3ℓ₃ − 7ℓ₄ + π₃ + 3π₄ |
| 22 | −0.00075 | 2ℓ₂ − 3ℓ₃ + π₃ |
| 23 | +0.00072 | π₁ + π₃ − 2Π − 2G |
| 24 | +0.00069 | π₄ − Π |
| 25 | −0.00058 | 2ℓ₃ − 3ℓ₄ + π₄ |
| 26 | −0.00057 | ℓ₃ − 2ℓ₄ + π₄ |
| 27 | +0.00056 | ℓ₃ + π₃ − 2Π − 2G |
| 28 | −0.00052 | ℓ₂ − 2ℓ₃ + π₁ |
| 29 | −0.00050 | π₂ − π₃ |
| 30 | +0.00048 | ℓ₃ − 2ℓ₄ + π₃ |
| 31 | −0.00045 | 2ℓ₂ − 3ℓ₃ + π₄ |
| 32 | −0.00041 | π₂ − π₄ |
| 33 | −0.00038 | 2G |
| 34 | −0.00037 | π₃ − π₄ + ω₃ − ω₄ |
| 35 | −0.00032 | 3ℓ₃ − 7ℓ₄ + 2π₃ + 2π₄ |
| 36 | +0.00030 | 4(ℓ₃ − ℓ₄) |
| 37 | +0.00029 | ℓ₃ + π₄ − 2Π − 2G |
| 38 | −0.00028 | ω₃ + ψ − 2Π − 2G |
| 39 | +0.00026 | ℓ₃ − Π − G |
| 40 | +0.00024 | ℓ₂ − 3ℓ₃ + 2ℓ₄ |
| 41 | +0.00021 | 2(ℓ₃ − Π − G) |
| 42 | −0.00021 | ℓ₃ − π₂ |
| 43 | +0.00017 | 2(ℓ₃ − π₃) |

**Σ₄ — Callisto (50 terms):**

| # | Amplitude | Argument |
|---|---|---|
| 1 | +0.84287 | ℓ₄ − π₄ |
| 2 | +0.03431 | π₄ − π₃ |
| 3 | −0.03305 | 2(ψ − Π) |
| 4 | −0.03211 | G |
| 5 | −0.01862 | ℓ₄ − π₃ |
| 6 | +0.01186 | ψ − ω₄ |
| 7 | +0.00623 | ℓ₄ + π₄ − 2G − 2Π |
| 8 | +0.00387 | 2(ℓ₄ − π₄) |
| 9 | −0.00284 | 5G′ − 2G + 52.225 |
| 10 | −0.00234 | 2(ψ − π₄) |
| 11 | −0.00223 | ℓ₃ − ℓ₄ |
| 12 | −0.00208 | ℓ₄ − Π |
| 13 | +0.00178 | ψ + ω₄ − 2π₄ |
| 14 | +0.00134 | π₄ − Π |
| 15 | +0.00125 | 2(ℓ₄ − G − Π) |
| 16 | −0.00117 | 2G |
| 17 | −0.00112 | 2(ℓ₃ − ℓ₄) |
| 18 | +0.00107 | 3ℓ₃ − 7ℓ₄ + 4π₄ |
| 19 | +0.00102 | ℓ₄ − G − Π |
| 20 | +0.00096 | 2ℓ₄ − ψ − ω₄ |
| 21 | +0.00087 | 2(ψ − ω₄) |
| 22 | −0.00085 | 3ℓ₃ − 7ℓ₄ + π₃ + 3π₄ |
| 23 | +0.00085 | ℓ₃ − 2ℓ₄ + π₄ |
| 24 | −0.00081 | 2(ℓ₄ − ψ) |
| 25 | +0.00071 | ℓ₄ + π₄ − 2Π − 3G |
| 26 | +0.00061 | ℓ₁ − ℓ₄ |
| 27 | −0.00056 | ψ − ω₃ |
| 28 | −0.00054 | ℓ₃ − 2ℓ₄ + π₃ |
| 29 | +0.00051 | ℓ₂ − ℓ₄ |
| 30 | +0.00042 | 2(ψ − G − Π) |
| 31 | +0.00039 | 2(π₄ − ω₄) |
| 32 | +0.00036 | ψ + Π − π₄ − ω₄ |
| 33 | +0.00035 | 2G′ − G + 188.37 |
| 34 | −0.00035 | ℓ₄ − π₄ + 2Π − 2ψ |
| 35 | −0.00032 | ℓ₄ + π₄ − 2Π − G |
| 36 | +0.00030 | 2G′ − 2G + 149.15 |
| 37 | +0.00029 | 3ℓ₃ − 7ℓ₄ + 2π₃ + 2π₄ |
| 38 | +0.00028 | ℓ₄ − π₄ + 2ψ − 2Π |
| 39 | −0.00028 | 2(ℓ₄ − ω₄) |
| 40 | −0.00027 | π₃ − π₄ + ω₃ − ω₄ |
| 41 | −0.00026 | 5G′ − 3G + 188.37 |
| 42 | +0.00025 | ω₄ − ω₃ |
| 43 | −0.00025 | ℓ₂ − 3ℓ₃ + 2ℓ₄ |
| 44 | −0.00023 | 3(ℓ₃ − ℓ₄) |
| 45 | +0.00021 | 2ℓ₄ − 2Π − 3G |
| 46 | −0.00021 | 2ℓ₃ − 3ℓ₄ + π₄ |
| 47 | +0.00019 | ℓ₄ − π₄ − G |
| 48 | −0.00019 | 2ℓ₄ − π₃ − π₄ |
| 49 | −0.00018 | ℓ₄ − π₄ + G |
| 50 | −0.00016 | ℓ₄ + π₃ − 2Π − 2G |

True longitudes (degrees):

```
Lᵢ = ℓᵢ + Σᵢ          (i = 1…4)
```

### 4.2 Latitude series (dimensionless; value of tan Bᵢ)

Each series sums `amplitude · sin(argument)`; then **Bᵢ = atan(series)**
(radians). Note carefully which arguments use the corrected longitudes Lᵢ,
which use the mean longitudes ℓᵢ, and the three composite arguments that
embed a **non-integer multiple of a longitude correction Σᵢ** — a uniform
integer-coefficient argument table cannot represent those three terms, so
they need explicit handling.

**tan B₁ — Io (7 terms):**

| # | Amplitude | Argument |
|---|---|---|
| 1 | +0.0006393 | L₁ − ω₁ |
| 2 | +0.0001825 | L₁ − ω₂ |
| 3 | +0.0000329 | L₁ − ω₃ |
| 4 | −0.0000311 | L₁ − ψ |
| 5 | +0.0000093 | L₁ − ω₄ |
| 6 | +0.0000075 | 3L₁ − 4ℓ₂ − 1.9927·Σ₁ + ω₂ |
| 7 | +0.0000046 | L₁ + ψ − 2Π − 2G |

**tan B₂ — Europa (9 terms):**

| # | Amplitude | Argument |
|---|---|---|
| 1 | +0.0081004 | L₂ − ω₂ |
| 2 | +0.0004512 | L₂ − ω₃ |
| 3 | −0.0003284 | L₂ − ψ |
| 4 | +0.0001160 | L₂ − ω₄ |
| 5 | +0.0000272 | ℓ₁ − 2ℓ₃ + 1.0146·Σ₂ + ω₂ |
| 6 | −0.0000144 | L₂ − ω₁ |
| 7 | +0.0000143 | L₂ + ψ − 2Π − 2G |
| 8 | +0.0000035 | L₂ − ψ + G |
| 9 | −0.0000028 | ℓ₁ − 2ℓ₃ + 1.0146·Σ₂ + ω₃ |

**tan B₃ — Ganymede (11 terms):**

| # | Amplitude | Argument |
|---|---|---|
| 1 | +0.0032402 | L₃ − ω₃ |
| 2 | −0.0016911 | L₃ − ψ |
| 3 | +0.0006847 | L₃ − ω₄ |
| 4 | −0.0002797 | L₃ − ω₂ |
| 5 | +0.0000321 | L₃ + ψ − 2Π − 2G |
| 6 | +0.0000051 | L₃ − ψ + G |
| 7 | −0.0000045 | L₃ − ψ − G |
| 8 | −0.0000045 | L₃ + ψ − 2Π |
| 9 | +0.0000037 | L₃ + ψ − 2Π − 3G |
| 10 | +0.0000030 | 2ℓ₂ − 3L₃ + 4.03·Σ₃ + ω₂ |
| 11 | −0.0000021 | 2ℓ₂ − 3L₃ + 4.03·Σ₃ + ω₃ |

**tan B₄ — Callisto (8 terms):**

| # | Amplitude | Argument |
|---|---|---|
| 1 | −0.0076579 | L₄ − ψ |
| 2 | +0.0044134 | L₄ − ω₄ |
| 3 | −0.0005112 | L₄ − ω₃ |
| 4 | +0.0000773 | L₄ + ψ − 2Π − 2G |
| 5 | +0.0000104 | L₄ − ψ + G |
| 6 | −0.0000102 | L₄ − ψ − G |
| 7 | +0.0000088 | L₄ + ψ − 2Π − 3G |
| 8 | −0.0000038 | L₄ + ψ − 2Π − G |

### 4.3 Radius series (dimensionless fractional corrections)

Each series sums `amplitude · cos(argument)`. All arguments use the
**mean** longitudes ℓᵢ (never Lᵢ).

**ΔR₁ — Io (7 terms):**

| # | Amplitude | Argument |
|---|---|---|
| 1 | −0.0041339 | 2(ℓ₁ − ℓ₂) |
| 2 | −0.0000387 | ℓ₁ − π₃ |
| 3 | −0.0000214 | ℓ₁ − π₄ |
| 4 | +0.0000170 | ℓ₁ − ℓ₂ |
| 5 | −0.0000131 | 4(ℓ₁ − ℓ₂) |
| 6 | +0.0000106 | ℓ₁ − ℓ₃ |
| 7 | −0.0000066 | ℓ₁ + π₃ − 2Π − 2G |

**ΔR₂ — Europa (11 terms):**

| # | Amplitude | Argument |
|---|---|---|
| 1 | +0.0093848 | ℓ₁ − ℓ₂ |
| 2 | −0.0003116 | ℓ₂ − π₃ |
| 3 | −0.0001744 | ℓ₂ − π₄ |
| 4 | −0.0001442 | ℓ₂ − π₂ |
| 5 | +0.0000553 | ℓ₂ − ℓ₃ |
| 6 | +0.0000523 | ℓ₁ − ℓ₃ |
| 7 | −0.0000290 | 2(ℓ₁ − ℓ₂) |
| 8 | +0.0000164 | 2(ℓ₂ − ω₂) |
| 9 | +0.0000107 | ℓ₁ − 2ℓ₃ + π₃ |
| 10 | −0.0000102 | ℓ₂ − π₁ |
| 11 | −0.0000091 | 2(ℓ₁ − ℓ₃) |

**ΔR₃ — Ganymede (10 terms):**

| # | Amplitude | Argument |
|---|---|---|
| 1 | −0.0014388 | ℓ₃ − π₃ |
| 2 | −0.0007919 | ℓ₃ − π₄ |
| 3 | +0.0006342 | ℓ₂ − ℓ₃ |
| 4 | −0.0001761 | 2(ℓ₃ − ℓ₄) |
| 5 | +0.0000294 | ℓ₃ − ℓ₄ |
| 6 | −0.0000156 | 3(ℓ₃ − ℓ₄) |
| 7 | +0.0000156 | ℓ₁ − ℓ₃ |
| 8 | −0.0000153 | ℓ₁ − ℓ₂ |
| 9 | +0.0000070 | 2ℓ₂ − 3ℓ₃ + π₃ |
| 10 | −0.0000051 | ℓ₃ + π₃ − 2Π − 2G |

**ΔR₄ — Callisto (16 terms):**

| # | Amplitude | Argument |
|---|---|---|
| 1 | −0.0073546 | ℓ₄ − π₄ |
| 2 | +0.0001621 | ℓ₄ − π₃ |
| 3 | +0.0000974 | ℓ₃ − ℓ₄ |
| 4 | −0.0000543 | ℓ₄ + π₄ − 2Π − 2G |
| 5 | −0.0000271 | 2(ℓ₄ − π₄) |
| 6 | +0.0000182 | ℓ₄ − Π |
| 7 | +0.0000177 | 2(ℓ₃ − ℓ₄) |
| 8 | −0.0000167 | 2ℓ₄ − ψ − ω₄ |
| 9 | +0.0000167 | ψ − ω₄ |
| 10 | −0.0000155 | 2(ℓ₄ − Π − G) |
| 11 | +0.0000142 | 2(ℓ₄ − ψ) |
| 12 | +0.0000105 | ℓ₁ − ℓ₄ |
| 13 | +0.0000092 | ℓ₂ − ℓ₄ |
| 14 | −0.0000089 | ℓ₄ − Π − G |
| 15 | −0.0000062 | ℓ₄ + π₄ − 2Π − 3G |
| 16 | +0.0000048 | 2(ℓ₄ − ω₄) |

Radius vectors, in units of Jupiter's equatorial radius:

```
R₁ = 5.90569  · (1 + ΔR₁)
R₂ = 9.39657  · (1 + ΔR₂)
R₃ = 14.98832 · (1 + ΔR₃)
R₄ = 26.36273 · (1 + ΔR₄)
```

Jupiter equatorial radius: **71492.0 km** (used only for the final
scaling in §6).

## 5. Precession adaptation (project deviation from the published recipe)

The E5 longitudes are referred to the mean equinox of B1950. The published
recipe precesses them to the equinox of date; **this module must not do
that**. Instead, apply the **fixed** precession arc from B1950.0
(JD 2433282.423) to J2000.0 (JD 2451545.0), so the output frame is
permanently J2000 (the consumer applies one fixed rotation to equatorial
ICRF afterwards):

```
T₀ = (2451545.0 − 2433282.423) / 36525
P  = 1.3966626 · T₀ + 0.0003088 · T₀²        (degrees)
```

Apply `+P` to the four true longitudes **and** to ψ:

```
Lᵢ′ = Lᵢ + P          (i = 1…4)
ψ′  = ψ + P
```

Use `Lᵢ′` and `ψ′` from here on. (Note the consequence: P cancels in the
difference Lᵢ′ − ψ′ used in §6.1, but ψ′ — not ψ — enters the node
rotation of §6.2. Both behaviors are required.)

## 6. Cartesian construction and frame chain

### 6.1 Jupiter-equatorial Cartesian coordinates

For each satellite i, with Bᵢ in radians from §4.2 and longitudes in
degrees:

```
Xᵢ = Rᵢ · cos(Lᵢ′ − ψ′) · cos(Bᵢ)
Yᵢ = Rᵢ · sin(Lᵢ′ − ψ′) · cos(Bᵢ)
Zᵢ = Rᵢ · sin(Bᵢ)
```

(Units: Jupiter radii. The X axis points toward the node ψ; Z toward
Jupiter's north pole.)

### 6.2 Rotation chain to the J2000 ecliptic

Constants, all held **fixed at their J2000 values** (project adaptation —
do not evaluate them at the date):

```
I  = 3.120262 + 0.0006 · (2451545.0 − 2415020.5) / 36525    (degrees)
       (inclination of Jupiter's equator on its orbital plane;
        the rate term is referred to epoch 1900.0 = JD 2415020.5
        and frozen at J2000)
iⱼ = 1.303267      (degrees; inclination of Jupiter's orbit on the ecliptic)
Ωⱼ = 100.464407    (degrees; ascending node of Jupiter's orbit)
```

Apply, in this exact order, to each satellite vector (x, y, z) starting
from (Xᵢ, Yᵢ, Zᵢ):

1. **Tilt equator → orbital plane** (rotation about the x-axis by I):
   ```
   x₁ = x
   y₁ = y·cos I − z·sin I
   z₁ = y·sin I + z·cos I
   ```
2. **Swing the node** (rotation about the z-axis by Φ = ψ′ − Ωⱼ):
   ```
   x₂ = x₁·cos Φ − y₁·sin Φ
   y₂ = x₁·sin Φ + y₁·cos Φ
   z₂ = z₁
   ```
3. **Tilt orbital plane → ecliptic** (rotation about the x-axis by iⱼ):
   ```
   x₃ = x₂
   y₃ = y₂·cos iⱼ − z₂·sin iⱼ
   z₃ = y₂·sin iⱼ + z₂·cos iⱼ
   ```
4. **Rotate to the equinox** (rotation about the z-axis by Ωⱼ):
   ```
   x₄ = x₃·cos Ωⱼ − y₃·sin Ωⱼ
   y₄ = x₃·sin Ωⱼ + y₃·cos Ωⱼ
   z₄ = z₃
   ```

### 6.3 Output scaling

Multiply each final vector by 71492.0 to convert Jupiter radii → km, and
return the four vectors as plain float 3-tuples in the order Io, Europa,
Ganymede, Callisto.

## 7. Floating-point guidance (required for baseline agreement)

The rewrite must reproduce the recorded numeric baseline to ≤1×10⁻³ km
per component (see §8). To stay at the centimetre level:

- Work in **degrees** end-to-end; convert to radians only at each trig
  call site, computing the conversion as `value * math.pi / 180.0`
  (multiply by π first, then divide by 180.0). Do **not** use
  `math.radians()` — its precomputed-constant rounding differs in the
  last ulp, which is visible at the ~0.1 m level through the large mean
  longitude arguments (|ℓ₁| reaches ~10⁷ degrees at the grid edges).
- Do **not** reduce any argument mod 360 anywhere.
- Build each trig argument with the multiplier structure exactly as
  printed in the tables — e.g. `2(ℓ₁ − ℓ₂)` as
  `2.0 * (l_io − l_europa)`-style "scale the difference", not
  `2·ℓ₁ − 2·ℓ₂` — and accumulate each series in table order.
- Use stdlib `math` floats throughout; no numpy, no `Fraction`.
- `Bᵢ` comes from `math.atan` in radians and feeds `sin`/`cos` directly.

## 8. Acceptance gates

Run from the repo root, in this order:

```
.venv/bin/python scripts/verify_galilean_clean_room.py --check tasks/results/galilean_baseline.json
```

must report **0 drifted** (worst moon delta ≤ 1e-3 km, worst COB delta
≤ 1e-12 AU over 9146 epochs spanning 1800–2200). Then:

```
pytest tests/test_planetary_moons.py tests/test_gas_giant_planet_center.py -q
poe lint
poe typecheck
```

must pass. If `--check` reports drift: a **periodic residual confined to
one satellite** means one mistranscribed coefficient in that satellite's
tables — recheck §4 against your implementation term by term; a
**broadband sub-centimetre residual** means floating-point re-association
— recheck §7.
