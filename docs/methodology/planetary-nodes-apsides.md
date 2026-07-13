# Planetary Nodes and Apsides

The `nod_aps()` and `nod_aps_ut()` entry points return four points on an
orbit: ascending node, descending node, periapsis, and apapsis (or the second
focus when `NODBIT_FOPOINT` is selected). This page describes LibEphemeris's
independent calculation and the compatibility semantics of their six output
slots.

## Clean-room basis

The implementation uses NASA JPL states, published mean orbital elements,
and IAU frame standards. Public compatibility calls are used only as an
ephemeral behavioral check: no comparison-package source, prose, algorithm,
data file, recorded output table, or fitted comparison data enters the project.

For `NODBIT_MEAN`, Mercury through Neptune use the long-term mean-element
polynomials derived from Simon et al. (1994). The apparent solar orbit uses
the corresponding Earth orbit. Pluto has no entry in that mean-element set
and therefore uses the osculating-state path. For `NODBIT_OSCU`, a JPL
position/velocity state is reduced to the angular-momentum and eccentricity
vectors of its instantaneous two-body orbit. `NODBIT_OSCU_BAR` uses the
solar-system-barycentric state for the supported outer bodies.

Each orbit-frame point is converted to the requested observer center and then
to the ecliptic frame of date. Geocentric planetary points receive annual
aberration unless `FLG_NOABERR` or `FLG_TRUEPOS` is set. They receive the
finite-source PPN gravitational-deflection reduction unless `FLG_NOGDEFL` or
`FLG_TRUEPOS` is set. `FLG_NOABERR` does not disable deflection.

Topocentric output subtracts the Skyfield/IERS geocentric observer state after
rotating it through the same Vondrák precession-nutation and ecliptic frame as
the pseudo-point. The velocity includes the derivative of that frame, and the
apparent-place branch adds the standard first-order diurnal-aberration term
from the observer's terrestrial-rotation velocity.

## MEAN speed-slot compatibility

The compatibility API does not encode a conventional instantaneous derivative
for MEAN points. LibEphemeris preserves its six-slot convention as follows:

1. Evaluate each center-relative mean point at `t` and `t + 1 day`.
2. Put the wrapped longitude difference and distance difference in the
   corresponding speed slots.
3. For nodes, put the latitude difference in the latitude-speed slot.
4. For periapsis and apapsis, put the **latitude at `t + 1 day`** in that
   slot, without subtracting the latitude at `t`.
5. Treat the resulting spherical six-tuple as a position/velocity state,
   convert it to Cartesian coordinates, and add or subtract the instantaneous
   Earth Cartesian state for heliocentric or geocentric output.

Step 4 is deliberately reproduced even though its slot is not dimensionally a
rate. It also changes the geocentric longitude and distance speeds through the
spherical-to-Cartesian Jacobian; replacing only that value after center
conversion would not be compatible.

The Moon uses the same slot rule. Its starting points are geocentric, so
heliocentric output adds the Earth state. Bodies on the osculating path retain
the near-instantaneous central-difference speed calculation.

## Apparent-place limitations

Outside the solar photosphere, LibEphemeris applies the independently
implemented finite-source PPN point-mass reduction. `FLG_NOGDEFL` suppresses
that term, while `FLG_TRUEPOS` suppresses both gravitational deflection and
annual aberration.

A point-mass model is singular and physically invalid when the observer-to-
pseudo-point ray crosses the solar photosphere. A compatible extended-Sun
reduction would require a separately sourced and validated solar mass profile.
LibEphemeris currently detects such rays using the IAU 2015 nominal solar
radius and omits the invalid point-mass term instead of fitting an empirical
interior model. Results for these rare interior rays therefore remain a known
compatibility gap; `FLG_NOGDEFL` and `FLG_TRUEPOS` avoid it deterministically.

Annual aberration is computed from the independently implemented first-order
direction transformation. Its position channels are supported, but the legacy
six-slot convention has a distinct radial-speed behavior that is not derived
from an astronomical standard. The `ddist` slot for apparent geocentric MEAN
points therefore remains a documented compatibility gap. No externally derived
coefficient or correction table is used to conceal it.

## References

- Simon, J. L., Bretagnon, P., Chapront, J., Chapront-Touzé, M., Francou, G.,
  and Laskar, J. (1994), “Numerical expressions for precession formulae and
  mean elements for the Moon and the planets,” *Astronomy & Astrophysics*,
  **282**, 663–683.
- Vondrák, J., Capitaine, N., and Wallace, P. (2011), “New precession
  expressions, valid for long time intervals,” *Astronomy & Astrophysics*,
  **534**, A22.
- Petit, G. and Luzum, B., eds. (2010), *IERS Conventions (2010)*, IERS
  Technical Note 36, Chapter 11.
- IAU 2015 Resolution B3, nominal solar and planetary conversion constants.
