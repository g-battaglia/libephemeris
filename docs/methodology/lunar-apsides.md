# Lunar Nodes and Apsides

LibEphemeris exposes three related lunar-point families: mean points,
interpolated compatibility curves, and instantaneous osculating points. They
have different definitions and should not be interchanged.

## Mean node and mean apogee

`MEAN_NODE` uses the IERS 2003 Delaunay argument for the mean longitude of the
lunar ascending node. `MEAN_APOG` follows from the Delaunay identities
`F = L - Omega` and `l = L - varpi`, then is projected through a conventional
`5.145°` lunar-orbit inclination. ERFA supplies the standards implementation of
the IERS arguments.

The mean points start on the mean ecliptic of date. Their public state helpers
use centered half-day derivatives and native Python floats.

## Interpolated apogee and perigee

`INTP_APOG` and `INTP_PERG` use separate trigonometric series in the lunar
Delaunay arguments:

```text
D   mean elongation of the Moon from the Sun
M   mean anomaly of the Sun
M'  mean anomaly of the Moon
F   mean argument of latitude
```

The apogee series models the principal evection harmonics and solar/latitude
couplings. The perigee series has the larger asymmetric evection spectrum and
extends through harmonic 30. This is why the two interpolated directions are
not generally separated by exactly `180°`.

The analytical result is refined by `libephemeris/lunar_apse_model.py`,
generated from the DE440 apsis passages by the committed
`scripts/generate_lunar_apse_model.py` and pinned by SHA-256: the reported
curve is exact at every DE440 apsis passage. Latitude follows a two-harmonic
inclined-plane fit (within `0.25°` of the instantaneous passage latitude).
The distance channels carry the mean DE440 passage distances rather than the
Moon's instantaneous distance.

Interpolated points start on the true ecliptic of date. Speeds are centered
half-day derivatives of the complete reported curve, including longitude
unwrapping.

## Osculating points

`TRUE_NODE` and `OSCU_APOG` use the active NASA JPL lunar state. With
geocentric position `r`, velocity `v`, and Earth-Moon gravitational parameter
`μ`, LibEphemeris evaluates

```text
h = r × v
e = (v × h) / μ - r / |r|
p = |h|² / μ
```

The node follows from the orbital angular-momentum vector. Osculating
perigee/apogee directions follow from `±e`, with two-body radii
`p/(1±|e|)`. These points can move rapidly because they describe the
instantaneous orbit.

## Frames and flags

Mean points begin on the mean ecliptic of date. Interpolated and osculating
points begin on the true ecliptic of date. Shared reduction code applies
nutation, J2000, sidereal, equatorial, Cartesian, radians, topocentric, and
speed-output contracts.

## Sources

- [IERS Conventions (2010), Chapter 5, Eq. 5.43](https://iers-conventions.obspm.fr/content/chapter5/icc5.pdf), Delaunay fundamental arguments.
- [ERFA](https://github.com/liberfa/erfa), BSD-licensed implementation of the
  IAU SOFA standards routines.
- Chapront-Touzé & Chapront (1988), *ELP 2000-82B: a semi-analytical lunar
  ephemeris adequate for historical times*, A&A 190, 342–352.
- [Park et al. (2021), DE440 and DE441](https://ssd.jpl.nasa.gov/doc/de440_de441.html).

See [Scientific Precision](../reference/precision.md).
