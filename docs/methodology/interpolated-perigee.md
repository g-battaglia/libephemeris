# Interpolated Lunar Perigee

`INTP_PERG` (body 22) is the smooth interpolated lunar-perigee curve. It is
evaluated independently from `INTP_APOG`; the two directions are not forced
to be antipodal.

The model is anchored to the Moon's actual perigee passages in the JPL DE440
ephemeris: every geocentric-distance minimum over ~1549--2651 CE, refined to
about one second. The longitude starts from the mean perigee and applies a
Delaunay perturbation series built from the ERFA/IERS 2003 arguments and fitted
by least squares at those passages; the larger evection spectrum extends
through the 30th harmonic and includes
solar-anomaly and latitude couplings. A spline-interpolated residual table
makes the curve accurate to `2 arcsec` at the sampled passages around J2000.
Across the full passage interval, its RMS/max residual is `0.80/6.35 arcsec`,
including the terminal taper. Tables and coefficients live in
`libephemeris/lunar_apse_model.py`, regenerated only by the committed
`scripts/generate_lunar_apse_model.py`; the integrity gate accepts only the
exact SHA-256-pinned module.

The generator reads only DE440 states and the BSD-licensed ERFA
implementations of the IAU SOFA fundamental arguments. Solar-anomaly
harmonics use the mean Earth-orbit eccentricity from Simon et al. (1994); no
comparison ephemeris supplies targets or coefficients. The table interval is
limited to DE440 medium coverage. A one-year terminal taper prevents jumps at
both ends; beyond it, trig-only accuracy is about `0.112°` RMS.

The residual figures above measure the model's fit to the DE440 apsis
passages, not agreement with any other implementation. Away from the
passages, interpolated apsides are definition-dependent: other published
realizations differ from this one by up to a few tenths of a degree in
longitude (and several arcminutes in latitude) at arbitrary dates, while
all agree at the physical DE440 passages themselves.

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
