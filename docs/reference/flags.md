# Flag Reference

Flags are bitmasks that control what is calculated and how results are returned.
`calc_ut()` / `calc()` return a 6-element tuple `(longitude, latitude, distance,
speed_lon, speed_lat, speed_dist)` and a return flag. Combine flags with `|`.

## Velocity

| Flag | Effect |
|------|--------|
| `FLG_SPEED` | Populate speed fields (pos[3]-pos[5]) with daily motion. Without this, speeds are zero. Almost every call should include it. |

## Observer

By default the observer is at Earth's center (geocentric).

| Flag | Effect |
|------|--------|
| `FLG_HELCTR` | Heliocentric: observer at the Sun. Distances become heliocentric AU. |
| `FLG_BARYCTR` | Barycentric: observer at the solar system barycenter. |
| `FLG_TOPOCTR` | Topocentric: observer on Earth's surface. Set position with `set_topo(lon, lat, alt)`. Matters most for the Moon (~1 deg parallax). |

## Coordinates

By default output is ecliptic longitude/latitude of date.

| Flag | Effect |
|------|--------|
| `FLG_EQUATORIAL` | Right Ascension (0-360 deg) and Declination (+/-90 deg) instead of ecliptic coordinates. |
| `FLG_XYZ` | Cartesian (x, y, z) in AU instead of spherical coordinates. |
| `FLG_RADIANS` | Angles in radians instead of degrees. |

## Reference Frame

By default positions are precessed to the equinox of date.

| Flag | Effect |
|------|--------|
| `FLG_J2000` | J2000.0 reference frame (no precession to date). |
| `FLG_NONUT` | Mean ecliptic/equator (no nutation applied). |
| `FLG_ICRS` | ICRS frame (no precession, no nutation, no frame bias). |

## Position Corrections

By default positions are apparent (light-time + aberration corrected).

| Flag | Effect |
|------|--------|
| `FLG_TRUEPOS` | Geometric position: no light-time correction. |
| `FLG_NOABERR` | Astrometric: light-time corrected, no aberration. |
| `FLG_NOGDEFL` | Skip gravitational deflection correction. |
| `FLG_ASTROMETRIC` | Shorthand for `FLG_NOABERR \| FLG_NOGDEFL`. |

## Sidereal Zodiac

| Flag | Effect |
|------|--------|
| `FLG_SIDEREAL` | Subtract ayanamsha from ecliptic longitude. Requires prior `set_sid_mode()` call to select ayanamsha (Lahiri, Fagan-Bradley, etc.). |
| `FLG_TROPICAL` | Tropical zodiac (value `0`); the default complement of `FLG_SIDEREAL`, listed for explicitness. |

## Compatibility

| Flag | Effect |
|------|--------|
| `FLG_JPLEPH` | Select the JPL ephemeris -- what libephemeris always uses, so effectively the default. |
| `FLG_SWIEPH` | Accepted for API compatibility, ignored -- all calculations use JPL DE440/DE441. |
| `FLG_MOSEPH` | Accepted for API compatibility, ignored. There is no Moshier fallback. |
| `FLG_DEFAULTEPH` | Default ephemeris selector (value `2`), equivalent to `FLG_SWIEPH`. |
| `FLG_SPEED3` | Converted to `FLG_SPEED` internally. |

## Advanced / specialized

These exist for pyswisseph API compatibility; most users never need them.

| Flag | Effect |
|------|--------|
| `FLG_CENTER_BODY` | Request the planet body center instead of the system barycenter. libephemeris already returns body centers by default (see [Planet Centers](../methodology/planet-centers-spk.md)). |
| `FLG_ORBEL_AA` | Use the Astronomical Almanac convention for osculating orbital elements (`get_orbital_elements`). |
| `FLG_JPLHOR` | JPL Horizons-consistent precession/nutation handling (legacy compatibility). |
| `FLG_JPLHOR_APPROX` | Approximate JPL Horizons mode (legacy compatibility). |
| `FLG_DPSIDEPS_1980` | IAU 1980 nutation-delta handling where applicable (legacy compatibility; shares its bit value with `FLG_JPLHOR`). |
| `FLG_TEST_PLMOON` | Internal planetary-moon test selector; not for general use. |

## Examples

```python
import libephemeris as swe
from libephemeris.constants import *

jd = swe.julday(2024, 3, 26, 12.0)

# Default: geocentric ecliptic apparent position with speed
pos, _ = swe.calc_ut(jd, MARS, FLG_SPEED)

# Heliocentric with speed
pos, _ = swe.calc_ut(jd, MARS, FLG_SPEED | FLG_HELCTR)

# Sidereal equatorial
swe.set_sid_mode(SIDM_LAHIRI)
pos, _ = swe.calc_ut(jd, SUN, FLG_SPEED | FLG_SIDEREAL | FLG_EQUATORIAL)

# J2000 ecliptic, no aberration
pos, _ = swe.calc_ut(jd, MOON, FLG_SPEED | FLG_J2000 | FLG_NOABERR)

# Topocentric (Rome)
swe.set_topo(12.4964, 41.9028, 0)
pos, _ = swe.calc_ut(jd, MOON, FLG_SPEED | FLG_TOPOCTR)
```
