# Python API

```python
import libephemeris as ephe
help(ephe.calc_ut)
help(ephe.EphemerisContext)
```

Installed docstrings document signatures, units, return layouts and errors.
The package includes `.pyi` stubs and `py.typed` for editors. The HTML
[API reference](api_reference.rst) is generated from the installed native wheel.
For a guided introduction, read the [manual](../beginners-manual/en/README.md).

## Time and coordinates

| Operation | Input/output convention |
| --- | --- |
| `julday` / `revjul` | Calendar ↔ Julian day; no change of time scale |
| `utc_to_jd` | Civil fields → `(jd_tt, jd_ut1)`, in that order |
| `calc` | TT Julian day |
| `calc_ut`, `houses` | UT1 Julian day |
| `deltat` | TT − UT1, returned in **days** |
| `calc_pctr` | TT at module level; **UT1** on `EphemerisContext` |
| `set_topo` | Longitude east, latitude north, height in metres |
| `houses` | Time, **latitude, longitude**, house selector |

The checked civil UTC conversion applies for Gregorian years 1972–2034.
Before 1972 and from 2035, `utc_to_jd` treats civil labels as UT1 and obtains TT
from the selected Delta-T model. It is not a future leap-second forecast.
Convert local wall time to UTC using its timezone and daylight-saving rules
before calling it. A two-part Julian date from a timescale handle is an
**additive pair**, not `(TT, UT1)`.

Ordinary `calc`/`calc_ut` return `(values, returned_flags)`; `values` contains:

| Index | Default spherical result | Unit |
| --- | --- | --- |
| 0, 1 | Ecliptic longitude, latitude | degrees |
| 2 | Distance | AU |
| 3, 4 | Longitude and latitude rates | degrees per TT day |
| 5 | Distance rate | AU per TT day |

The default position is apparent and geocentric. Request `FLG_SPEED` for rates.
`FLG_EQUATORIAL` replaces longitude/latitude with RA/declination; **RA is in
degrees**, not hours. `FLG_XYZ` selects Cartesian AU and AU/TT-day;
`FLG_RADIANS` changes spherical angles, not Cartesian components.
`FLG_J2000`, `FLG_NONUT`, `FLG_TRUEPOS`, `FLG_HELCTR`, `FLG_BARYCTR`,
`FLG_TOPOCTR` and `FLG_SIDEREAL` change frame, corrections or origin: consult
the function's flag rules before combining them. `ECL_NUT` is a special
four-angle result, not a position/velocity vector.

House cusps normally start at index zero (house one); selectors such as
Gauquelin have different layouts. Eclipse, visibility and phenomena arrays
have their own indexed contracts. Returned flags describe the served request,
not authenticated data provenance.

## Chart conventions

- **Tropical or sidereal:** ordinary positions are tropical. Select an ayanamsha
  with `set_sid_mode` and request `FLG_SIDEREAL` for sidereal positions; setting
  the mode alone does not change an unflagged calculation. Use `houses_ex` for
  house calculations with sidereal flags.
- **Retrograde motion:** request `FLG_SPEED` and inspect the longitude rate,
  `values[3]`. Its sign describes motion in the selected frame, not every frame.
- **Aspects:** compare longitudes in the same zodiac and frame. Wrap angular
  differences across 0°/360°; aspect definitions and orbs are application choices.
- **Houses and lunar points:** choose the house system and node/apogee definition
  explicitly. Mean, true, osculating and interpolated points are not synonyms;
  some house systems are undefined at polar latitudes.

The [cookbook](../beginners-manual/en/15-cookbook.md) develops these conventions
into natal-chart, transit, sidereal-chart and solar-return recipes.

## Independent observers

Configure shared defaults before creating workers; keep observer/sidereal
settings in separate contexts:

```python
import libephemeris as ephe

ephe.set_calc_mode("leb")
ephe.set_precision_tier("base")
ephe.set_network_policy("sealed")
ephe.set_iers_delta_t_enabled(False)
rome = ephe.EphemerisContext()
rome.set_topo(12.5, 41.9, 30.0)
_, jd_ut1 = ephe.utc_to_jd(2024, 4, 8, 12, 0, 0.0)
moon, flags = rome.calc_ut(jd_ut1, ephe.MOON, ephe.FLG_SPEED | ephe.FLG_TOPOCTR)
```

Native work releases the GIL but remains synchronous. In an async application,
offload complete calculations to bounded workers; cancellation of the waiting
request does not stop a native solver. Do not mutate global configuration
between concurrent requests. Global `close()` releases shared resources and
resets some defaults; it is **not per-request cleanup**. Retained handles may
outlive a resource selection. See the runnable
[worker example](../python/examples/async_worker.py).

## Find an operation

| Task | Entry points |
| --- | --- |
| Positions and stars | `calc_ut`, `calc_pctr`, `fixstar_ut` |
| Houses and sidereal coordinates | `houses`, `houses_ex`, `set_sid_mode`, `get_ayanamsa_ut` |
| Crossings and stations | `solcross_ut`, `mooncross_ut`, `find_station_ut` |
| Eclipses and occultations | `sol_eclipse_when_glob`, `lun_eclipse_when`, `lun_occult_when_glob` |
| Local sky and visibility | `azalt`, `rise_trans`, `heliacal_ut` |
| Resources | `EphemerisContext`, `get_leb_reader`, `download_for_tier` |

Some no-event results are status codes, `None` or zero contact times rather
than exceptions. Catch each operation's documented exceptions: not all derive
from `ephe.Error`. A zero/unbracketed contact does not necessarily mean the
physical event cannot occur.

For native diagnostics, inspect `help(ephe.tracing)` and its exported capture
helpers. For existing applications, read [migration](../MIGRATION.md).
Historical module names such as `rebound_integration` do not imply those
third-party engines are runtime dependencies.
