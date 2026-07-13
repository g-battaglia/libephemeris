# Osculating Lunar Apogee (True Lilith)

`OSCU_APOG` (body 13), often called True Lilith, is the apogee of the Moon's
instantaneous osculating two-body orbit. LibEphemeris derives it independently
from NASA JPL DE440/DE441 geocentric Moon state vectors.

## State-vector construction

For geocentric position `r`, velocity `v`, and the Earth-Moon gravitational
parameter `μ`:

1. Compute specific angular momentum `h = r × v`.
2. Compute the eccentricity vector `e = (v × h) / μ - r / |r|`.
3. The eccentricity vector points toward osculating perigee; the apogee
   direction is `-e`.
4. Convert `-e` to longitude and latitude in the true ecliptic frame of date.
5. With semi-latus rectum `p = |h|² / μ`, return apogee distance
   `p / (1 - |e|)`.

The same state gives the osculating perigee direction `+e` and distance
`p / (1 + |e|)`. The calculation uses the IAU 2015 Earth/Moon mass ratio and
nominal Earth gravitational parameter converted to AU³/day².

This is a geometric osculating construction. No fitted lunar perturbation
series is added: solar and planetary perturbations are already present in the
JPL-integrated state. The result can move rapidly because an osculating ellipse
is the instantaneous Keplerian orbit, not a smoothed physical track.

## Relation to other apogee points

- `MEAN_APOG` is the smooth mean curve built from ERFA/IERS Delaunay arguments
  and conventional inclined-orbit geometry.
- `INTP_APOG` is a compatibility-oriented Delaunay-argument perturbation series
  with a hash-pinned, versioned refinement.
- `OSCU_APOG` is the instantaneous eccentricity-vector result described here.

These points answer different questions and should not be interchanged. For an
actual lunar distance event, calculate the Moon's geocentric distance and
search for a local extremum.

## Frames and speeds

The raw state is evaluated in the true ecliptic of date. The shared output
pipeline removes nutation for `FLG_NONUT`/J2000 paths and applies sidereal,
equatorial, Cartesian, and radians formatting. Requested speeds are numerical
derivatives of the reported osculating curve with a short stencil chosen for
its rapid motion. Heliocentric and barycentric requests return the abstract
point zero vector for compatibility.

## Provenance

The implementation follows standard state-vector orbital mechanics and NASA
JPL/IAU data. It contains no fitted compatibility coefficients or sampled
third-party ephemeris output.

See [Lunar Nodes and Apsides](lunar-apsides.md) and
[Scientific Precision](../reference/precision.md).
