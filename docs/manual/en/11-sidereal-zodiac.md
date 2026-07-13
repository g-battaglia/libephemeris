# Chapter 11 — The sidereal zodiac and ayanamsha

## What you will learn

This chapter explains the tropical/sidereal distinction, the predefined modes
implemented by LibEphemeris, their source-audit status, and how to supply your
own sourced definition.

## 11.1 Two longitude origins

The tropical zodiac takes the moving vernal equinox as its longitude origin.
A sidereal zodiac instead uses an origin fixed to a stellar or other celestial
reference. Precession changes the angle between those origins over time.

The angle subtracted from tropical longitude is the **ayanamsha**:

```text
sidereal longitude = tropical longitude - ayanamsha
```

There is no single universal sidereal origin. Historical traditions choose
different stars, epochs, or geometric constructions. LibEphemeris therefore
requires a reproducible independent definition for every mode it computes.

## 11.2 Mode classification

All canonical `SIDM_*` constants are exported and every predefined base ID
0--46 computes without a fallback warning. `SIDM_USER` remains available for a
caller-defined reference epoch and offset.

The stellar modes propagate independently sourced catalogue astrometry through
the normal Skyfield/ERFA apparent-place pipeline. Galactic modes use the
published IAU galactic frame and their cited geometric construction. Fixed-epoch
modes are frame projections, not undocumented longitude tables.

The numerical-source audit is tracked separately from runtime support. The
authoritative [mode table](../../reference/ayanamsha.md) marks each mode's
defining epoch and primary-publication status.

## 11.3 Calculating a sidereal position

Choose a mode, then request `FLG_SIDEREAL`:

```python
import libephemeris as ephem

jd_ut = ephem.julday(2024, 4, 8, 12.0)
ephem.set_sid_mode(ephem.SIDM_TRUE_CITRA)

ayanamsha = ephem.get_ayanamsa_ut(jd_ut)
position, retflag = ephem.calc_ut(
    jd_ut,
    ephem.SUN,
    ephem.FLG_SIDEREAL | ephem.FLG_SPEED,
)
```

`FLG_SIDEREAL` can be combined with supported coordinate and speed flags. The
selected mode is process-global in the classic API; use `EphemerisContext` when
parallel calculations need different sidereal configurations.

```python
ctx = ephem.EphemerisContext()
ctx.set_sid_mode(ephem.SIDM_J2000)
position, retflag = ctx.calc_ut(
    jd_ut,
    ephem.MOON,
    ephem.FLG_SIDEREAL | ephem.FLG_SPEED,
)
```

## 11.4 User-defined ayanamsha

Use `SIDM_USER` when you have an independently sourced reference epoch and
offset:

```python
import libephemeris as ephem

ephem.set_sid_mode(
    ephem.SIDM_USER,
    t0=source_reference_jd,
    ayan_t0=source_offset_degrees,
)
```

LibEphemeris propagates that zero point with its independent long-term
precession implementation. You are responsible for the values' provenance and
time-scale convention. Pass the intended epoch explicitly rather than relying
on a sentinel.

## 11.5 Mean and true ayanamsha APIs

`get_ayanamsa_ut(jd_ut)` and `get_ayanamsa(jd_tt)` return the mean offset. The
extended variants also return a flag value and honor the requested nutation
convention:

```python
retflag, value = ephem.get_ayanamsa_ex_ut(jd_ut, ephem.FLG_SWIEPH)
```

Planetary, stellar, Skyfield, Horizons, and LEB paths use the same sidereal
dispatcher. Backend choice must not silently select a different anchor.

## 11.6 Fixed-epoch modes

`SIDM_J2000`, `SIDM_J1900`, and `SIDM_B1950` request coordinates in the named
mean reference frame. With ecliptic output, coordinates are projected onto that
epoch's mean ecliptic and equinox; equatorial output uses its mean equatorial
frame. These modes should not be interpreted as traditional astrological
ayanamshas.

## Summary

- Sidereal longitude subtracts an ayanamsha from tropical longitude.
- Every predefined base `SIDM_*` mode 0--46 computes without fallback.
- The source-audit status of each historical anchor is documented separately.
- Use `SIDM_USER` for a lawful, citable tradition-specific definition.
- External implementation output is used only for ephemeral comparison, never
  as an anchor source.

See [Sidereal modes](../../reference/ayanamsha.md) for the authoritative mode
list and provenance notes.
