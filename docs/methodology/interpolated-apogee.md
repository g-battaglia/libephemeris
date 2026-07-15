# Interpolated Lunar Apogee

`INTP_APOG` (body 21) is the smooth interpolated lunar-apogee curve. It is
distinct from both `MEAN_APOG` and the instantaneous, JPL-derived
`OSCU_APOG`.

The model is anchored to the Moon's actual apogee passages in the JPL DE440
ephemeris (Park et al. 2021, AJ 161:105): every geocentric-distance maximum
over ~1549--2651 CE, refined to about one second. The longitude starts from
the IERS 2003 mean lunar apse implemented by ERFA, adds a Delaunay-argument
series in `D`, `M`, `M'`, and `F` fitted by least squares at those passages
(dominant evection harmonics, solar-anomaly couplings, and a latitude
coupling), and finishes with a spline-interpolated residual table. Around J2000
the resulting longitude is within `12 arcsec` of the DE440 apogee-passage
longitude. Across the full passage interval, its RMS/max residual is
`2.66/18.39 arcsec`, including the terminal taper. Tables and coefficients live
in
`libephemeris/lunar_apse_model.py`, regenerated only by the committed
`scripts/generate_lunar_apse_model.py` and accepted by the integrity gate at
its exact SHA-256.

The generator reads only DE440 states and the BSD-licensed ERFA
implementations of the IAU SOFA fundamental arguments. Solar-anomaly
harmonics use the mean Earth-orbit eccentricity from Simon et al. (1994); no
comparison ephemeris supplies targets or coefficients. The table interval is
limited to DE440 medium coverage. A one-year terminal taper prevents jumps at
both ends; beyond it, trig-only accuracy is about `0.045°` RMS.

The latitude uses a two-harmonic inclined-plane model fitted to the same
DE440 passages, and the distance channel is the mean DE440 apogee-passage
distance, `0.0027100 AU`. The result is on the true ecliptic of date. Shared
reduction code handles `FLG_NONUT`, `FLG_J2000`, `FLG_EQUATORIAL`,
`FLG_SIDEREAL`, `FLG_XYZ`, and `FLG_RADIANS`. `FLG_SPEED` uses the centered
half-day derivative.

```python
import libephemeris as eph

position, flags = eph.calc_ut(
    2451545.0,
    eph.INTP_APOG,
    eph.FLG_SPEED,
)
lon, lat, distance, dlon, dlat, ddistance = position
```

This abstract curve is not a Moon-distance event finder. For physical apogee
timing, search the JPL Moon's geocentric distance for a local maximum.

See [Lunar Nodes and Apsides](lunar-apsides.md) for formulas, frames, and
sources.
