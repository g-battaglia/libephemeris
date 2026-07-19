# pyswisseph compatibility (drop-in) and gotchas

libephemeris targets the public pyswisseph interface: canonical function and
constant names, parameter order, tuple structure, flag encoding, state
transitions, and documented exceptions. The numerics are independently
implemented from JPL, IAU/ERFA, primary literature, and permissively licensed
catalogues — but the *contract* is the same. Most migrations are a one-line
import change.

```python
# before:  import swisseph as swe
import libephemeris as swe

pos, retflag = swe.calc_ut(jd_ut, swe.SUN, swe.FLG_SPEED)
```

## Naming rules

- **Bare canonical names only.** Functions: `calc_ut`, `houses_ex2`,
  `fixstar_ut`, `sol_eclipse_when_glob`, `nod_aps_ut`. Constants: `SUN`, `MOON`,
  `FLG_SPEED`, `FLG_SIDEREAL`, `SIDM_LAHIRI`, `NODBIT_MEAN`. There is **no**
  `swe_` prefix on functions and **no** `SE_` / `SEFLG_` prefix on constants.
- Only documented back-compat alias: `SE_FNAME_DE431` (== `FNAME_DE431`).
- If your legacy code used prefixed names via a shim, drop the prefixes.
- Every numeric return is a native Python `float`, not a numpy scalar — safe for
  JSON, `round()`, f-strings, and equality without `.item()`.

## Argument-order gotchas (these bite silently)

- **`houses(tjdut, lat, lon, hsys=ord('P'), iflag=0)`** — latitude **before**
  longitude, matching pyswisseph. Same order for `houses_ex`, `houses_ex2`,
  `houses_armc*`. Swapping lat/lon yields a wrong-but-plausible chart with no
  error.
- **`hsys` accepts an int code or a bytes char**: `ord('P')` and `b'P'` both
  select Placidus. (A bare `str` like `'P'` is not the documented form — use
  `ord('P')` or `b'P'`.)
- **`calc_ut` takes UT1; `calc` takes TT (ET).** They differ by ΔT (tens of
  seconds today, minutes in antiquity). Astrology almost always wants `calc_ut`.
- Time helpers mirror pyswisseph: `julday(y, m, d, hour, cal=GREG_CAL)`,
  `revjul`, `deltat`, `utc_to_jd`, `sidtime`. `date_conversion(y,m,d,hour,cal)`
  returns `(valid, jd, (y, m, d, hour))`.

## Return shapes (match pyswisseph)

| Call | Returns |
|------|---------|
| `calc_ut` / `calc` | `((lon, lat, dist, spd_lon, spd_lat, spd_dist), retflag)` |
| `houses` / `houses_ex` | `(cusps[12], ascmc[8])` |
| `houses_ex2` | `(cusps[12], ascmc[8], cusps_speed[12], ascmc_speed[8])` |
| `fixstar_ut` / `fixstar` | `(pos[6], resolved_name, retflag)` |
| `nod_aps_ut` / `nod_aps` | `(asc_node, desc_node, perihelion, aphelion)` — each `[6]` |
| `pheno_ut` / `pheno` | phenomena tuple (phase angle, phase, elongation, apparent diameter, magnitude, …) |

`ascmc` order: `[Asc, MC, ARMC, Vertex, EquAsc, CoAsc, CoAscKoch, PolarAsc]`.

## Flags

Common flags carry pyswisseph semantics: `FLG_SPEED`, `FLG_HELCTR`,
`FLG_TOPOCTR` (needs `set_topo`), `FLG_SIDEREAL` (needs `set_sid_mode`),
`FLG_EQUATORIAL`, `FLG_XYZ`, `FLG_RADIANS`, `FLG_J2000`, `FLG_NONUT`,
`FLG_TRUEPOS`, `FLG_NOABERR`, `FLG_NOGDEFL`.

**Ephemeris-selection flags are compatibility no-ops for source choice.**
`FLG_SWIEPH`, `FLG_JPLEPH`, and `FLG_MOSEPH` are accepted, but the actual source
is governed by the **calc mode** (see `calc-modes.md`), not by these flags —
`FLG_MOSEPH` in particular is ignored; calculations always use JPL data. The
`retflag` echoes the effective flags.

Sidereal work is unchanged: `set_sid_mode(SIDM_LAHIRI)` then pass
`FLG_SIDEREAL`; query the ayanamsha with `get_ayanamsa_ut(jd)`.

## Intentional divergences (documented, bounded)

libephemeris aims for indistinguishability from the reference for modern-era
positions, cusps, eclipses, and most sidereal modes (sub-arcsecond). A small set
of behaviors differ **on purpose**, each following a published definition or a
checkable physical invariant rather than fitting another implementation's
output. Do not "fix" these to match pyswisseph:

- `SIDEREAL | J2000` is honored uniformly on all lunar points (keeps
  `|TrueNode − MeanNode|` within its physical ±1.5° bound far from J2000).
- Of-date mean obliquity uses a self-consistent Vondrák frame (the Sun stays at
  ~0° ecliptic latitude at every epoch).
- House-cusp *speeds* are the true derivative of the reported cusp position
  (Whole Sign cusp speeds are `0.0`).
- Total-eclipse obscuration reports `attr[2] = 1.0` during totality (bounded
  covered-area fraction).
- Some historical hypothetical bodies (IDs 48, 49, 54, 55, 57, 58) keep their
  names but raise `UnknownBodyError` — no source-complete model was recovered,
  so they fail closed instead of returning an untraceable orbit.

The full list with rationale and measured bounds lives in the repo docs
(`docs/comparison/intentional-divergences.md`,
`docs/comparison/known-differences.md`).

## libephemeris-only extensions (not reference-compatible)

These have libephemeris-defined contracts — do not treat them as pyswisseph
entry points:

- Performance/state: `reset_session()` (fast state reset preserving readers and
  caches), idempotent `set_ephe_path()`, `release_data_cache()`.
- Offline/config: `set_calc_mode`/`get_calc_mode`, `set_leb_file`,
  `get_leb_reader`, `set_network_policy`, `set_precision_tier`.
- Introspection: `get_body_coverage`, `coverage`, `get_leb_inventory`,
  `get_runtime_data_requirements`, `inspect_leb_file`, `start_tracing` /
  `get_trace_results`.
- House fallbacks: `houses_with_fallback`, `houses_armc_with_fallback`.
- Thread-safe alternative to global state: `EphemerisContext` — use it for
  concurrent/multi-threaded workloads instead of the module-level globals.

## A note on sources

Do not consult, copy, or reference Swiss Ephemeris source or data files when
using libephemeris — it is a self-contained JPL/IAU/ERFA implementation. Install
the wheel, provision tiers with `libephemeris download`, and use the API above.
