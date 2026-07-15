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

The baseline apse and node come from the same ERFA/IERS model used by the mean
points. The analytical result is refined by `libephemeris/lunar_apse_model.py`,
generated from the DE440 apsis passages by the committed
`scripts/generate_lunar_apse_model.py` and pinned by SHA-256. The generator
uses `fad03`, `falp03`, `fal03`, and `faf03` directly; it has no compatibility
ephemeris target. On a two-year interval around J2000, passage-longitude errors
are below `12 arcsec` for apogee and `2 arcsec` for perigee. Latitude follows a
two-harmonic inclined-plane fit (within about `0.2°` of the instantaneous
passage latitude). The distance channels carry the mean DE440 passage
distances rather than the Moon's instantaneous distance.

Across all 14,581 passages per side, including both terminal intervals, the
fixed-grid residual has RMS/max error `2.66/18.39 arcsec` for apogee and
`0.80/6.35 arcsec` for perigee. The one-year terminal taper is part of the
generator and avoids a correction jump at the first and last DE440 passages.

Interpolated points start on the true ecliptic of date. Speeds are centered
half-day derivatives of the complete reported curve, including longitude
unwrapping.

The fitted correction tables cover the DE440 medium interval, approximately
1550--2650 CE. Outside the terminal taper, the standards-derived baseline plus
trigonometric series remains; accuracy then has the documented trig-only RMS
of about `0.045°` (apogee) and `0.112°` (perigee). This extrapolated curve is an
abstract coordinate, not a substitute for DE441 lunar distance extrema.

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
- Simon et al. (1994), *Numerical expressions for precession formulae and mean
  elements for the Moon and the planets*, A&A 282, 663--683, for the mean
  Earth-orbit eccentricity used to scale solar-anomaly harmonics.
- [Park et al. (2021), DE440 and DE441](https://ssd.jpl.nasa.gov/doc/de440_de441.html).

See [Scientific Precision](../reference/precision.md).
