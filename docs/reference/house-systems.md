# House System Algorithms and Mathematical Formulas

This document provides mathematical documentation for the 25 house systems
(26 codes including the A/E equal-house alias) implemented in LibEphemeris.
Each section states the construction evaluated by this project. Historical
attribution and implementation provenance are separate questions: a system can
have a traditional name while its local Cartesian/spherical derivation is
project-authored from the public definition.

## Provenance contract and source status

No cusp is obtained from a copied cusp table or by fitting another program's
output. General coordinate transformations follow Meeus chapter 13 and the
IERS/ERFA frame conventions. ARMC, nutation, precession, and obliquity come from
the source-mapped IERS/ERFA/Vondrak pipeline. Those astronomy standards do not
define a house system; the division/projection rule must be accounted for
separately below.

The evidence labels are deliberately precise:

- **Author definition**: a public description by the system's author or named
  originator identifies the construction.
- **Published secondary definition**: an independent published treatment
  describes the construction, but this audit does not claim a primary-author
  transcription.
- **Project-defined variant**: the formula is completely disclosed here and in
  code, but no primary historical locator is asserted for the name/variant.
- **Mathematical identity**: the rule is only an equal division, midpoint, or
  other equation printed in full; no empirical coefficient exists.

| Code | Construction source/status | Local derivation boundary |
|---|---|---|
| `P` | **Published secondary definition:** Holden (1977); independent time-division description in Urania Trust, “The Astronomy of Houses.” | Project fixed-point solution of semi-arc trisection; convergence and polar failure policy are local choices. |
| `K` | **Published secondary definition:** Holden (1977) and Urania Trust's Koch description. | Project iterative oblique-ascension solution; tolerances/fallbacks are local. |
| `R` | **Published secondary definition:** Holden (1977); equator division projected through horizon north/south points. | Project spherical projection. |
| `C` | **Published secondary definition:** Holden (1977); prime-vertical division projected through horizon north/south points. | Project spherical projection. |
| `O` | **Published secondary definition / identity:** trisection of the four ecliptic angle arcs. | Direct modular-arc arithmetic. |
| `A`, `E` | **Mathematical identity:** twelve 30-degree sectors from the Ascendant. | `A` is an API alias of `E`. |
| `W` | **Mathematical identity:** sign containing the Ascendant becomes house 1. | Zero-based sign-boundary arithmetic. |
| `B` | **Published secondary definition:** Holden (1977) and independent Alcabitius histories. | Project semi-arc/oblique-ascension implementation and polar policy. |
| `T` | **Published definition:** Polich & Page (1961), *The Topocentric System of Houses*. | Project trigonometric evaluation of the published pole construction. |
| `M` | **Published secondary definition:** Holden (1977), equator division projected through the ecliptic poles. | Project spherical projection. |
| `X` | **Published secondary definition:** axial/meridian equatorial division. | Project RA-to-ecliptic mapping. |
| `H` | **Published secondary definition:** horizon/azimuth division projected by vertical circles. | Project vector/trigonometric projection and degeneracy policy. |
| `U` | **Author definition:** Bogdan Krusinski's public 1995 description. | Independent Cartesian construction in `house_constructions.py`. |
| `Y` | **Independent published description:** Ingmar de Boer, “APC Houses,” describing the Ascendant Parallel Circle attributed to L. Knegt. | Independent Cartesian construction in `house_constructions.py`; no claim that de Boer is the originator. |
| `I` | **Published modern description:** Sunshine system using the Sun's diurnal/nocturnal arc; the `I` branch is the documented Treindl solution. | Project numerical solution and fallback policy. |
| `i` | **Author definition:** Bob Makransky, “The Sunshine House System.” | Project spherical-trigonometry implementation of the upper/lower-meridian convention. |
| `J` | **Author definition:** John Savard, “Astrological House Systems,” Albategnus/Savard-A entry. | Independent Cartesian construction in `house_constructions.py`. |
| `L`, `Q` | **Author definition:** Walter D. Pullen's public “Sinusoidal House Systems” specification. | Project evaluation of the stated delta/ratio interpolation, including explicit narrow-quadrant policy. |
| `S` | **Project-defined API realization:** midpoint construction over Porphyry boundaries, conventionally named Sripati. | Formula is completely stated below; no primary-author locator is asserted here. |
| `V` | **Project-defined API realization:** equal sectors with the Ascendant at the middle of house 1. | Exact 15-degree offset is the complete rule; historical-priority claim is not needed for reproduction. |
| `F` | **Project-defined API realization:** 30-degree RA sectors from the Ascendant RA, projected along hour circles. | Exact projection is stated below; a page-level Carter/Cope locator remains unasserted. |
| `N` | **Project-defined API realization:** fixed natural-zodiac sectors beginning at 0 Aries. | Exact formula contains no observer-dependent astronomy. |
| `D` | **Project-defined API realization:** twelve equal sectors anchored to the MC. | Exact 30-degree formula is printed below. |
| `G` | **Published sector definition:** 36 Gauquelin diurnal/nocturnal sectors. | Project semi-arc solver, indexing, convergence, and polar failure behavior. |

“Project-defined” is not a hidden-source placeholder. It says exactly that the
equation shown in this document is the project's source of truth and that the
historical label has not been upgraded into an unsupported bibliographic claim.
See `house_constructions.py` for the fully commented vector frames, oriented
arc invariants, branch selection, and degeneracy handling of `J`, `U`, `Y`, and
`i`.

## Table of Contents

1. [Terminology and Definitions](#terminology-and-definitions)
2. [Placidus (P)](#placidus-p)
3. [Koch (K)](#koch-k)
4. [Regiomontanus (R)](#regiomontanus-r)
5. [Campanus (C)](#campanus-c)
6. [Porphyry (O)](#porphyry-o)
7. [Equal (A/E)](#equal-ae)
8. [Whole Sign (W)](#whole-sign-w)
9. [Alcabitius (B)](#alcabitius-b)
10. [Topocentric/Polich-Page (T)](#topocentricpolich-page-t)
11. [Morinus (M)](#morinus-m)
12. [Meridian/Axial Rotation (X)](#meridianaxial-rotation-x)
13. [Vehlow (V)](#vehlow-v)
14. [Carter Poli-Equatorial (F)](#carter-poli-equatorial-f)
15. [Gauquelin Sectors (G)](#gauquelin-sectors-g)
16. [Horizontal/Azimuthal (H)](#horizontalazimuthal-h)
17. [Krusinski-Pisa (U)](#krusinski-pisa-u)
18. [Natural Gradient (N)](#natural-gradient-n)
19. [APC (Y)](#apc-y)
20. [Sripati (S)](#sripati-s)
21. [Pullen SD / Neo-Porphyry (L)](#pullen-sd--neo-porphyry-l)
22. [Pullen SR (Q)](#pullen-sr-q)
23. [Equal from MC (D)](#equal-from-mc-d)
24. [Sunshine / Treindl (I)](#sunshine--treindl-i)
25. [Sunshine / Makransky (i)](#sunshine--makransky-i)
26. [Savard-A (J)](#savard-a-j)

---

## Terminology and Definitions

### Coordinate Systems

- **Ecliptic coordinates (λ, β)**: Longitude and latitude measured along the ecliptic plane
- **Equatorial coordinates (α, δ)**: Right Ascension (RA) and Declination
- **Horizontal coordinates (A, h)**: Azimuth and Altitude

### Key Variables

| Symbol | Definition | Units |
|--------|-----------|-------|
| λ | Ecliptic longitude | degrees (0-360) |
| ε | Obliquity of the ecliptic | degrees (~23.44°) |
| φ | Geographic latitude | degrees (-90 to +90) |
| θ | Local Sidereal Time (LST) | degrees |
| ARMC | Right Ascension of MC | degrees |
| α | Right Ascension | degrees |
| δ | Declination | degrees |
| AD | Ascensional Difference | degrees |
| SA | Semi-diurnal Arc | degrees |
| OA | Oblique Ascension | degrees |
| H | Hour Angle | degrees |

### Fundamental Formulas

**Ecliptic to Equatorial Conversion:**
```
tan(α) = sin(λ) · cos(ε) / cos(λ)
sin(δ) = sin(λ) · sin(ε)
```

**Equatorial to Ecliptic Conversion:**
```
tan(λ) = sin(α) / (cos(α) · cos(ε))
```

**Ascensional Difference:**
```
sin(AD) = tan(φ) · tan(δ)
```

**Semi-diurnal Arc:**
```
SA = 90° + AD  (diurnal, above horizon)
SA = 90° - AD  (nocturnal, below horizon)
```

---

## Placidus (P)

**Historical Origin:** Placidus de Titis (1603-1668), Italian mathematician and monk.

**Conceptual Basis:** Time-based division of the semi-diurnal and semi-nocturnal arcs into thirds.

### Algorithm

1. Calculate Semi-diurnal Arc (SA) for each point
2. Divide SA into three equal time segments
3. Iteratively find ecliptic longitude at each time division
4. Houses 11, 12 from MC to Asc (above horizon)
5. Houses 2, 3 from Asc to IC (below horizon)
6. Opposite houses (5, 6, 8, 9) are 180° from (11, 12, 2, 3)

### Mathematical Formulas

**House 11 (1/3 of semi-arc from MC):**
```
H₁₁ = SA/3 = (90° + AD)/3
RA₁₁ = ARMC + H₁₁
```

**House 12 (2/3 of semi-arc from MC):**
```
H₁₂ = 2·SA/3 = 2·(90° + AD)/3
RA₁₂ = ARMC + H₁₂
```

**House 2 (2/3 of nocturnal semi-arc from IC):**
```
H₂ = 2·(90° - AD)/3
RA₂ = ARMC + 180° - H₂
```

**House 3 (1/3 of nocturnal semi-arc from IC):**
```
H₃ = (90° - AD)/3
RA₃ = ARMC + 180° - H₃
```

**Iterative Solution:**
For each cusp, iterate until convergence (threshold: 1e-7°):
```
1. tan(δ) = sin(RA) · tan(ε)
2. AD = arcsin(tan(φ) · tan(δ))
3. RA_new = ARMC + H (computed from AD)
4. If |RA_new - RA| < 1e-7°, stop
```

**RA to Ecliptic Longitude:**
```
tan(λ) = sin(RA) / (cos(RA) · cos(ε))
λ = atan2(sin(RA), cos(RA) · cos(ε))
```

### Polar Limitation

Fails when |φ| + ε > 90° because tan(φ) · tan(δ) > 1, making AD undefined.

---

## Koch (K)

**Historical Origin:** Walter Koch (1895-1970), German astrologer. Also called GOH (Geburtsort-Häusersystem) or Birthplace system.

**Conceptual Basis:** Trisection of Oblique Ascension intervals between angles.

### Algorithm

1. Calculate Oblique Ascension for MC, Asc, and IC
2. Divide OA intervals into thirds
3. Iteratively solve for ecliptic longitude at each target OA
4. Calculate opposite houses by adding 180°

### Mathematical Formulas

**Oblique Ascension:**
```
OA = RA - AD
where AD = arcsin(tan(φ) · tan(δ))
```

**For houses 11, 12 (MC to Asc quadrant):**
```
OA₁₁ = OA_MC + (OA_Asc - OA_MC)/3
OA₁₂ = OA_MC + 2·(OA_Asc - OA_MC)/3
```

**For houses 2, 3 (Asc to IC quadrant):**
```
OA₂ = OA_Asc + (OA_IC - OA_Asc)/3
OA₃ = OA_Asc + 2·(OA_IC - OA_Asc)/3
```

**Iterative Solution for target OA:**
```
1. RA = target_OA (initial guess)
2. tan(δ) = sin(RA) · tan(ε)
3. AD = arcsin(tan(φ) · tan(δ))
4. RA_new = target_OA + AD
5. Repeat until |RA_new - RA| < 1e-7°
6. Convert final RA to ecliptic longitude
```

### Polar Limitation

Same as Placidus: fails when |φ| + ε > 90°.

---

## Regiomontanus (R)

**Historical Origin:** Johannes Müller von Königsberg (Regiomontanus, 1436-1476).

**Conceptual Basis:** Equal 30° divisions of the celestial equator, projected to ecliptic via great circles through the north and south points of the horizon.

### Algorithm

1. Divide equator into 30° segments from ARMC
2. Calculate pole of projection for each segment
3. Project equator point to ecliptic using spherical trigonometry

### Mathematical Formulas

**Pole of projection:**
```
tan(P) = tan(φ) · sin(H)
where H = 30°, 60°, 120°, 150° for houses 11, 12, 2, 3
```

**Right Ascension offset:**
```
R = ARMC + H - 90°
```

**Ecliptic longitude:**
```
λ = atan2(cos(R), -(sin(R)·cos(ε) + tan(P)·sin(ε)))
```

---

## Campanus (C)

**Historical Origin:** Giovanni di Campani (13th century).

**Conceptual Basis:** Equal 30° divisions of the prime vertical, projected onto the ecliptic.

### Algorithm

1. Divide prime vertical into 30° segments
2. Transform each segment to equatorial hour angle
3. Apply Regiomontanus-style projection

### Mathematical Formulas

**Prime vertical to hour angle transformation:**
```
tan(H_eff) = tan(h_pv) · cos(φ)
where h_pv = 30°, 60°, 120°, 150° (prime vertical offset)
```

**Pole calculation:**
```
tan(P) = tan(φ) · sin(H_eff)
```

**Ecliptic projection:**
```
R = ARMC + H_eff - 90°
λ = atan2(cos(R), -(sin(R)·cos(ε) + tan(P)·sin(ε)))
```

---

## Porphyry (O)

**Historical Origin:** Porphyry of Tyre (234-305 CE), Neoplatonic philosopher.

**Conceptual Basis:** Space-based trisection of ecliptic arc between angles.

### Algorithm

1. Divide ecliptic arc from MC to Asc into three equal parts
2. Divide ecliptic arc from Asc to IC into three equal parts
3. Opposite houses are 180° apart

### Mathematical Formulas

**Quadrant 10→1 (MC to Asc):**
```
step = (λ_Asc - λ_MC) mod 360° / 3
λ₁₁ = λ_MC + step
λ₁₂ = λ_MC + 2·step
```

**Quadrant 1→4 (Asc to IC):**
```
step = (λ_IC - λ_Asc) mod 360° / 3
λ₂ = λ_Asc + step
λ₃ = λ_Asc + 2·step
```

**Opposite houses:**
```
λ₅ = (λ₁₁ + 180°) mod 360°
λ₆ = (λ₁₂ + 180°) mod 360°
λ₈ = (λ₂ + 180°) mod 360°
λ₉ = (λ₃ + 180°) mod 360°
```

**Advantages:** Works at all latitudes, computationally simple.

---

## Equal (A/E)

**Historical Origin:** Ancient Hellenistic astrology.

**Conceptual Basis:** Each house is exactly 30° of ecliptic longitude.

### Mathematical Formula

```
λᵢ = (λ_Asc + (i-1)·30°) mod 360°
where i = 1, 2, ..., 12
```

**Properties:**
- Simplest calculation
- Works at all latitudes
- MC may not coincide with 10th house cusp

---

## Whole Sign (W)

**Historical Origin:** Oldest house system, used in Hellenistic, Indian, and Medieval astrology.

**Conceptual Basis:** Each house is one complete zodiac sign; House 1 is the sign containing the Ascendant.

### Mathematical Formula

```
start = floor(λ_Asc / 30°) × 30°
λᵢ = (start + (i-1)·30°) mod 360°
where i = 1, 2, ..., 12
```

**Example:** If Asc = 15° Taurus (45°), then:
- House 1 starts at 30° (0° Taurus)
- House 2 starts at 60° (0° Gemini)
- etc.

---

## Alcabitius (B)

**Historical Origin:** Abd al-Aziz al-Qabisi (Alcabitius, 10th century), Arab astrologer.

**Conceptual Basis:** Time trisection of Ascendant's diurnal arc, projected by hour circles.

### Algorithm

1. Calculate Right Ascension of Ascendant
2. Divide RA intervals between angles
3. Convert each RA division back to ecliptic longitude

### Mathematical Formulas

**Quadrant MC to Asc:**
```
arc = (RA_Asc - ARMC) mod 360°
step = arc / 3
RA₁₁ = ARMC + step
RA₁₂ = ARMC + 2·step
```

**Quadrant Asc to IC:**
```
arc = (RA_IC - RA_Asc) mod 360°
step = arc / 3
RA₂ = RA_Asc + step
RA₃ = RA_Asc + 2·step
```

**RA to ecliptic conversion:**
```
λ = atan2(sin(RA), cos(RA)·cos(ε))
```

---

## Topocentric/Polich-Page (T)

**Historical Origin:** Wendel Polich and A.P. Nelson Page (1961).

**Conceptual Basis:** Modified pole method to account for observer's actual position on Earth's surface.

### Algorithm

Similar to Regiomontanus but with modified pole factors.

### Mathematical Formulas

**Pole factors:**
```
For houses 11, 3: factor = 1/3
For houses 12, 2: factor = 2/3
```

**Pole calculation:**
```
tan(P) = tan(φ) · factor
```

**Ecliptic projection:**
```
R = ARMC + H - 90°
where H = 30°, 60°, 120°, 150° for houses 11, 12, 2, 3
λ = atan2(cos(R), -(sin(R)·cos(ε) + tan(P)·sin(ε)))
```

---

## Morinus (M)

**Historical Origin:** Jean-Baptiste Morin (1583-1656), French astrologer.

**Conceptual Basis:** Equal 30° divisions on celestial equator, projected to ecliptic via ecliptic poles.

### Algorithm

Project equator points (ARMC + n·30°) to ecliptic.

### Mathematical Formula

```
For each house cusp at RA = ARMC + (i-1)·30°:
tan(λ) = tan(RA) · cos(ε)
λ = atan2(sin(RA)·cos(ε), cos(RA))
```

**Properties:**
- Location-independent (ignores latitude)
- Morinus 10th cusp differs from standard MC

---

## Meridian/Axial Rotation (X)

**Historical Origin:** Also known as Zariel system.

**Conceptual Basis:** Equal 30° RA divisions from MC, projected to ecliptic via celestial poles.

### Mathematical Formula

```
For each house cusp at RA = ARMC + (i-10)·30°:
λ = atan2(sin(RA), cos(RA)·cos(ε))
```

**Difference from Morinus:**
```
Morinus: λ = atan2(sin(RA)·cos(ε), cos(RA))
Meridian: λ = atan2(sin(RA), cos(RA)·cos(ε))
```

---

## Vehlow (V)

**Source status:** Project-defined API realization carrying the conventional
Vehlow name. This audit does not claim a page-level Vehlow transcription for
the exact 15-degree offset below; the displayed identity is the complete local
definition.

**Conceptual Basis:** Equal houses with Ascendant at center (15°) of House 1 rather than at cusp.

### Mathematical Formula

```
start = (λ_Asc - 15°) mod 360°
λᵢ = (start + (i-1)·30°) mod 360°
```

**Result:** Ascendant falls at 15° into House 1.

---

## Carter Poli-Equatorial (F)

**Source status:** Project-defined API realization carrying the conventional
Carter/Cope name. A page-level primary locator for this exact projection is not
asserted; the equations below are the complete reproducible definition.

**Conceptual Basis:** Equal 30° divisions on equator starting from RA of Ascendant, projected to ecliptic.

### Mathematical Formula

```
RA_Asc = atan2(sin(λ_Asc)·cos(ε), cos(λ_Asc))

For each house i (1-12):
RA = RA_Asc + (i-1)·30°
λᵢ = atan2(sin(RA), cos(RA)·cos(ε))
```

---

## Gauquelin Sectors (G)

**Historical Origin:** Michel and Françoise Gauquelin (20th century), for statistical research.

**Conceptual Basis:** 36 clockwise sectors obtained by dividing each
diurnal/nocturnal quadrant into ninths. The boundaries are position circles,
not uniform right-ascension cuts on the equator.

### Independent numerical realization

The cardinal boundaries are fixed first:

```
sector 1  = Ascendant
sector 10 = Midheaven
sector 19 = Descendant
sector 28 = Imum Coeli
```

For each of the eight intermediate boundaries in a quadrant, let `k = 1..8`
be its ninth-arc index and compute the latitude/obliquity ascensional difference

```
a = asin(clamp(tan φ · tan ε, -1, +1)).
```

The trial meridian is `RAMC ± k·10°`, with the sign and base meridian selected
for the quadrant. The initial pole height of its position circle is

```
p₀ = atan(sin(k·a/9) / tan ε).
```

`Asc(r, ε, φ, p)` denotes the ecliptic intersection returned by the common
house-circle projection for trial meridian `r` and pole height `p`. Starting
with `λ₀ = Asc(r, ε, φ, p₀)`, the code refines

```
δⱼ     = asin(sin ε · sin λⱼ)
ADⱼ    = asin(clamp(tan φ · tan δⱼ, -1, +1))
pⱼ₊₁   = atan(sin(k·ADⱼ/9) / tan δⱼ)
λⱼ₊₁   = Asc(r, ε, φ, pⱼ₊₁)
```

until the wrapped longitude change is below `1e-8°`, with at most 100
iterations. Opposite boundaries are exact antipodes. The iteration limit,
convergence threshold, `1e-10` degeneracy guard, and `clamp` operations are
project numerical choices; the ninths-of-semi-arc construction is the sourced
Gauquelin definition.

### Polar policy

Within the polar circle (`|φ| >= 90° - ε`) ordinary rising/setting semi-arcs
are undefined. The API therefore returns an explicit project fallback of 36
clockwise equal 10-degree sectors from the Ascendant; it does not present that
fallback as the historical Gauquelin construction.

---

## Horizontal/Azimuthal (H)

**Historical Origin:** Traditional horizon-based system.

**Conceptual Basis:** House circles based on horizon plane divisions.

### Algorithm

The horizon is divided into twelve equal 30° arcs from the East Point.
Through each division point a vertical circle (great circle through zenith
and nadir) is drawn; its ecliptic intersection is the house cusp. Derived
from the geometric definition in R. W. Holden, *The Elements of House
Division* (1977), ch. 10, with the spherical triangle solved via Smart,
*Textbook on Spherical Astronomy*, ch. 3.

### Mathematical Formulas

**Co-latitude** (signed complement of geographic latitude):
```
φ' = 90° - φ   (for φ > 0)
φ' = -90° - φ  (for φ < 0)
```

**Vertical-circle parameters** at each horizon division angle θ (measured
from the meridian; θ = k·30°, k = 1..5), from the spherical triangle
(zenith, celestial pole, division point):
```
pole_height = arcsin(sin θ · sin φ')
ra_offset   = atan2(cos θ, sin θ · cos φ')
```

**House cusps** (base reference θ₀ = ARMC + 180°; Asc = rising-point
formula at the given RA and pole height):
```
cusp = Asc(θ₀ + 90° + ra_offset, ε, φ, pole_height)

θ = 30° → cusp 3    θ = 60° → cusp 2    θ = 90° → cusp 1 (East Point)
θ = 120° → cusp 12   θ = 150° → cusp 11
```
(At θ = 90°: pole_height = arcsin(sin φ') = φ' and ra_offset = 0, so
cusp 1 falls on the East Point. The remaining cusps follow by symmetry.)

---

## Krusinski-Pisa (U)

**Author source:** Bogdan Krusinski,
[“The Krusinski House System”](https://www.astrologia.pl/house-system.html)
(1995 public description): divide the Ascendant--zenith great circle into
30-degree arcs and project its division points to the ecliptic along meridian
(hour) circles.

**Conceptual Basis:** Great circle passing through Ascendant and Zenith, divided into 12 equal parts.

### Independent Cartesian construction

`house_constructions.py` evaluates that definition directly in the
right-handed equatorial frame of date. For Ascendant longitude `a`, obliquity
`ε`, ARMC `m`, and latitude `φ`, form the Ascendant and zenith unit vectors:

```
A = (cos a, sin a cos ε, sin a sin ε)
Z = (cos φ cos m, cos φ sin m, sin φ)
U = normalize(Z - (Z·A) A)
```

`U` is the zenith-facing tangent of the Ascendant--zenith great circle. Its
division point at signed arc `θ` is therefore

```
P(θ) = cos θ A + sin θ U
θ = 0°, +30°, +60°, +90°, -30°, -60°
      H1   H12    H11    H10    H2    H3
```

The hour circle through `P` has `RA = atan2(P_y, P_x)` and cuts the ecliptic at

```
λ = atan2(sin RA, cos RA cos ε).
```

The `atan2` expression preserves the quadrant without a copied case table.
Cusps 4--9 are exact antipodes of 10--3. The local Gram--Schmidt step above,
floating-point normalization, and antipodal construction are project-authored
numerical realization; the 30-degree division and hour-circle projection are
Krusinski's published definition.

---

## Natural Gradient (N)

**Source status:** Project-defined API realization. The fixed zodiac partition
printed below is the complete rule; no historical-origin claim is needed to
reproduce it.

**Conceptual Basis:** Equal houses starting from 0° Aries.

### Mathematical Formula

```
λᵢ = (i-1) × 30°
where i = 1, 2, ..., 12
```

House 1 = 0° Aries, House 2 = 0° Taurus, etc.

---

## APC (Y)

**Source:** Ingmar de Boer, [“APC Houses”](https://www.ingmardeboer.nl/index.php?title=APC_Houses), an independent construction description attributing the system to L. Knegt and the Dutch School of Ram.

**Conceptual Basis:** The small circle parallel to the celestial equator that
passes through the Ascendant is divided separately into six equal arcs above
and six below the horizon. House circles through the north and south points of
the horizon carry those division anchors to the ecliptic. A declination
parallel is a small circle, not a “great circle parallel to the horizon.”

### Algorithm

Uses the declination of the Ascendant to define the house circle.

The formula below is the spherical expression of that definition. The runtime
implementation in `house_constructions.py` constructs the same anchors and
planes with Cartesian vectors, selects the ecliptic intersection by oriented
arc containment, and explicitly validates cusp 1 = Ascendant and cusp 10 = MC.

### Mathematical Formulas

**Ascensional difference of Ascendant:**
```
kv = arctan(tan(φ)·tan(ε)·cos(ARMC) / (1 + tan(φ)·tan(ε)·sin(ARMC)))
```

**Declination of Ascendant:**
```
δ_Asc = arctan(sin(kv) / tan(φ))
```

**House cusp calculation (for house n):**
```
For n < 8 (below horizon):
  k = n - 1
  a = kv + ARMC + π/2 + k·(π/2 - kv)/3

For n >= 8 (above horizon):
  k = n - 13
  a = kv + ARMC + π/2 + k·(π/2 + kv)/3

λ = atan2(tan(δ_Asc)·tan(φ)·sin(ARMC) + sin(a),
          cos(ε)·(tan(δ_Asc)·tan(φ)·cos(ARMC) + cos(a))
          + sin(ε)·tan(φ)·sin(ARMC - a))
```

---

## Sripati (S)

**Source status:** Project-defined API realization conventionally named
Sripati. This document does not assert that the exact midpoint rule below is a
page-level transcription from Sripati.

**Conceptual Basis:** Midpoints of Porphyry cusps. Each cusp is the midpoint between consecutive Porphyry cusps, effectively shifting Porphyry houses by half a house.

### Mathematical Formula

```
1. Compute Porphyry cusps P₁..P₁₂
2. Sripati cusp Sᵢ = midpoint(Pᵢ, Pᵢ₊₁)
```

**Properties:** Works at all latitudes, same computational simplicity as Porphyry.
For public-API compatibility, `house_pos()` collapses the Sripati wrap interval
above house number 12 to exactly `1.0`; this includes house 12 and the first
half of house 1. The cusp-12 boundary itself remains `12.0`.

---

## Pullen SD / Neo-Porphyry (L)

**Author source:** Walter D. Pullen,
[“Sinusoidal House Systems”](https://www.astrolog.org/astrolog/astsine.htm),
“Sinusoidal Delta Houses” (first described in 1994).

**Conceptual Basis:** Sinusoidal delta interpolation between angles. Instead of linear trisection (Porphyry), uses a sinusoidal curve to distribute cusps more smoothly between the angles.

### Algorithm

For the smaller of a pair of opposite quadrant arcs, let its size be `q` and
write the three consecutive house widths as `x+n, x, x+n`. The opposite
quadrant widths are `x+3n, x+4n, x+3n`. Solving the two sum constraints gives:

```
n = (90 - q) / 4
x = (q - 30) / 2
```

The implementation evaluates the algebra symmetrically for every quadrant as:

```
d = quadrant_size - 90
side   = 30 + d/4
middle = 30 + d/2
```

When a quadrant is smaller than 30 degrees, the author specification says the
middle house is pinned to zero and the two side houses bisect the quadrant.
That narrow-quadrant rule is explicit code, not an inferred correction.

**Properties:** Works at all latitudes, smoother cusp distribution than Porphyry.

---

## Pullen SR (Q)

**Author source:** Walter D. Pullen,
[“Sinusoidal House Systems”](https://www.astrolog.org/astrolog/astsine.htm),
“Sinusoidal Ratio Houses” (first implemented in 2016).

**Conceptual Basis:** Sinusoidal ratio interpolation between angles. Similar to Pullen SD but uses a ratio-based sinusoidal formula.

### Algorithm

For quadrant size `q`, the three widths are `r*x, x, r*x`; the opposite
quadrant uses `r^3*x, r^4*x, r^3*x`. Therefore:

```
x * (2*r + 1) = q
x * r^3 * (2 + r) = 180 - q
r^3 * (2 + r) / (2*r + 1) = (180 - q) / q
```

LibEphemeris solves the positive `r` root with at most 20 Newton iterations,
an initial value of 1, a derivative guard of `1e-12`, a positive floor of
`0.01`, and a step convergence threshold of `1e-10`. Those are disclosed
project numerical choices; the width equations are Pullen's definition.

**Properties:** Works at all latitudes, alternative sinusoidal distribution to Pullen SD.

---

## Equal from MC (D)

**Source status:** Project-defined API realization. The exact MC-anchored
30-degree identity printed below is the complete rule.

**Conceptual Basis:** Equal 30-degree houses starting from the MC (10th house cusp) rather than the Ascendant.

### Mathematical Formula

```
λ₁₀ = MC
λᵢ = (MC + (i-10)·30°) mod 360°
where i = 1, 2, ..., 12
```

**Properties:** MC always coincides with the 10th house cusp. Ascendant may not coincide with the 1st house cusp. Works at all latitudes.

---

## Sunshine / Treindl (I)

**Source status:** Treindl numerical solution of Bob Makransky's published
Sunshine construction. This code identifier distinguishes the solution method,
not a different physical Sun path.

**Conceptual Basis:** Houses based on the Sun's diurnal arc, dividing the sky according to the Sun's actual path above and below the horizon. Uses the Sun's declination to define house boundaries.

### Algorithm

1. Compute the Sun's ascensional difference
   `AD = asin(tan(sun_declination) * tan(latitude))`.
2. Form the diurnal semi-arc `90 + AD` and nocturnal semi-arc `90 - AD`.
3. Trisect each semi-arc to obtain the eight intermediate anchors.
4. Project each anchor to the ecliptic through the corresponding house circle
   containing the horizon north/south points.
5. Apply the explicitly documented MC-below-horizon orientation branch.

If the requested ephemeris cannot supply solar declination, the caller uses the
Meeus chapter 25 low-precision apparent-Sun fallback documented in
`houses.py`; zero is never silently substituted.

**Properties:** Requires Sun's declination as input. Latitude-dependent.

---

## Sunshine / Makransky (i)

**Author source:** Bob Makransky, “The Sunshine House System,” public author
article reproduced in [Astrology House Systems](https://www.astro.nu/2017/03/03/astrology-house-systems/).

**Conceptual Basis:** Similar to Sunshine/Treindl but with an alternative division algorithm for the Sun's diurnal arc.

### Algorithm

The Sun's diurnal and nocturnal arcs are divided into six parts, then projected
through horizon-north/south house circles. `house_constructions.py` follows the
author's upper/lower-meridian spherical-trigonometry convention and documents:

- the equatorial frame, hour-angle sign, and horizon basis;
- which semi-arc supplies each cusp;
- the pole height and RA of each house-circle/equator intersection;
- ecliptic-intersection orientation; and
- polar/degenerate guards.

**Properties:** Requires Sun's declination as input. Latitude-dependent.

---

## Savard-A (J)

**Author source:** John Savard,
[“Astrological House Systems”](http://www.quadibloc.com/other/as01.htm),
Albategnus/Savard-A entry.

**Conceptual Basis:** Declination parallels cut the east-point-to-zenith arc of
the prime vertical at one third and two thirds of the observer latitude. Great
circles through each anchor and the horizon's north/south points intersect the
ecliptic to form the intermediate cusps.

### Algorithm

`house_constructions.py` derives every anchor in the equatorial frame of date,
forms the plane through the horizon north/south axis, intersects it with the
ecliptic, and selects the correct antipodal intersection with an oriented-arc
invariant. Opposite cusps are exact antipodes. No lookup or fitted branch table
is used.

**Properties:** Latitude-dependent.

---

## Summary Comparison

| Code | System | Calculation Type | Latitude Independent | Works at Poles |
|------|--------|-----------------|---------------------|----------------|
| P | Placidus | Iterative | No | No |
| K | Koch | Iterative | No | No |
| R | Regiomontanus | Direct | No | Yes* |
| C | Campanus | Direct | No | Yes* |
| O | Porphyry | Direct | No | Yes |
| E/A | Equal | Simple | Yes | Yes |
| W | Whole Sign | Simple | Yes | Yes |
| B | Alcabitius | Direct | No | Yes* |
| T | Topocentric | Direct | No | Yes* |
| M | Morinus | Direct | Yes | Yes |
| X | Meridian | Direct | Yes | Yes |
| V | Vehlow | Simple | Yes | Yes |
| F | Carter | Direct | No | Yes |
| G | Gauquelin | Time-based | No | No |
| H | Horizontal | Complex | No | Yes* |
| U | Krusinski | Complex | No | Yes* |
| N | Natural Gradient | Simple | Yes | Yes |
| Y | APC | Complex | No | Yes* |
| S | Sripati | Direct | No | Yes |
| L | Pullen SD | Direct | No | Yes |
| Q | Pullen SR | Direct | No | Yes |
| D | Equal from MC | Simple | Yes | Yes |
| I | Sunshine (Treindl) | Complex | No | No |
| i | Sunshine (Makransky) | Complex | No | No |
| J | Savard-A | Direct | No | Yes* |

25 distinct systems, 26 character codes (A is an alias for E; both compute equal houses from the Ascendant).

*May require special handling near poles.

---

## References

1. Meeus, Jean, *Astronomical Algorithms*, 2nd ed. (1998), chapter 13,
   [WorldCat record](https://search.worldcat.org/title/40521322). Used for
   general spherical-coordinate transformations, not as the definition of all
   house systems.
2. Petit, G. and Luzum, B. (eds.), *IERS Conventions (2010)*,
   [IERS Technical Note 36](https://iers-conventions.obspm.fr/content/tn36.pdf).
   Used for time/frame conventions surrounding the house geometry.
3. Vondrak, J., Capitaine, N., and Wallace, P. (2011), “New precession
   expressions, valid for long time intervals,”
   [A&A 534:A22](https://doi.org/10.1051/0004-6361/201117274), with the
   published corrigendum. Used for the common long-term frame, not a house
   division rule.
4. Holden, Ralph William, *The Elements of House Division* (1977),
   [WorldCat OCLC 4549972](https://search.worldcat.org/title/4549972).
5. Bates, Graham, [“The Astronomy of Houses”](https://www.uraniatrust.org/astrology/astronomy-of-houses),
   Urania Trust: independent construction descriptions and bibliography for
   the major ecliptic, space, and time systems.
6. Houlding, Deborah,
   [“House Division Principles”](https://www.skyscript.co.uk/glossary/house-division/),
   Skyscript: independent historical/construction overview.
7. Polich, Wendel, and A. P. Nelson Page, *The Topocentric System of Houses*
   (1961); the public Urania Trust treatment above identifies the radius/pole
   construction and the authors' 1964 *Spica* article.
8. Krusinski, Bogdan,
   [“The Krusinski House System”](https://www.astrologia.pl/house-system.html)
   (author's 1995 description).
9. de Boer, Ingmar,
   [“APC Houses”](https://www.ingmardeboer.nl/index.php?title=APC_Houses).
10. Makransky, Bob, “The Sunshine House System,” public author article
    [reproduced here](https://www.astro.nu/2017/03/03/astrology-house-systems/).
11. Savard, John,
    [“Astrological House Systems”](http://www.quadibloc.com/other/as01.htm),
    Albategnus/Savard-A entry.
12. Pullen, Walter D.,
    [“Sinusoidal House Systems”](https://www.astrolog.org/astrolog/astsine.htm)
    (author specification for SD and SR).

These sources define or independently describe the geometry; all prose and
code in this project are independently written. No source code, source
comments, documentation prose, algorithm transcription, or data file from the
API-compatibility target is used in this document or implementation.
