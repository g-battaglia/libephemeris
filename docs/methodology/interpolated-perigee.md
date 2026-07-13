# Interpolated Lunar Perigee

`INTP_PERG` (body 22) is the smooth interpolated lunar-perigee curve. It is
evaluated independently from `INTP_APOG`; the two directions are not forced
to be antipodal.

The model is anchored to the Moon's actual perigee passages in the JPL DE440
ephemeris: every geocentric-distance minimum over ~1549--2651 CE, refined to
about one second. The longitude starts from the mean perigee and applies a
Delaunay perturbation series fitted by least squares at those passages; the
larger evection spectrum extends through the 30th harmonic and includes
solar-anomaly and latitude couplings. A spline-interpolated residual table
makes the curve exact at every DE440 perigee passage. Tables and
coefficients live in `libephemeris/lunar_apse_model.py`, regenerated only by
the committed `scripts/generate_lunar_apse_model.py`; the integrity gate
accepts only the exact SHA-256-pinned module.

Latitude uses a two-harmonic inclined-plane model fitted to the same DE440
passages. The distance channel is the mean DE440 perigee-passage distance,
`0.0024236 AU`. Output starts on the true ecliptic of date, and the common
reduction pipeline applies frame, representation, sidereal, and nutation
flags. `FLG_SPEED` uses the centered half-day derivative.

```python
import libephemeris as eph

position, flags = eph.calc_ut(
    2451545.0,
    eph.INTP_PERG,
    eph.FLG_SPEED,
)
lon, lat, distance, dlon, dlat, ddistance = position
```

This coordinate is not a proxy for tides or a physical "supermoon" event. For
those questions, search the JPL Moon's geocentric distance for a local minimum.

See [Lunar Nodes and Apsides](lunar-apsides.md) for formulas, frames, and
sources.
