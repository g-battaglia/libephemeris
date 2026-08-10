---
name: libephemeris
description: >-
  Use when writing Python that computes astronomical or astrological data with
  the `libephemeris` package (`import libephemeris`) — planet/point positions
  (calc_ut / calc), house cusps (houses / houses_ex2), fixed stars, eclipses,
  lunar nodes and apsides — or when porting `pyswisseph` / `swisseph` code to
  libephemeris as a drop-in. Covers the four calc modes (auto / skyfield / leb /
  horizons), the sealed offline LEB mode and its network policy, the cumulative
  base / medium / extended data tiers with best-by-date routing, per-body
  coverage (get_body_coverage), the typed error contract (EphemerisRangeError,
  UnknownBodyError, LEBCorruptionError, NetworkSealedError), and data
  provisioning (`libephemeris download`). Load this BEFORE calling any
  libephemeris function so field names, mode values, tier ranges, and error
  types are exact rather than recalled from pyswisseph memory.
---

# libephemeris

`libephemeris` is a pure-Python astronomical ephemeris library powered by NASA
JPL DE440 / DE441. It is a **drop-in replacement for and superset of
pyswisseph**: same function names, parameter order, tuple shapes, flag encoding,
and exception surface, plus its own performance/offline extensions. This skill
targets an agent writing Python **against** the library (not developing it).

Key facts to anchor on (verified against the installed package, `3.1.0`):

- The public API uses **canonical bare names** — `calc_ut`, `SUN`, `FLG_SPEED`,
  `SIDM_LAHIRI` — never the `swe_` / `SE_` / `SEFLG_` prefixes. (The lone
  documented back-compat alias is `SE_FNAME_DE431`.)
- Import convention: `import libephemeris as swe` — most pyswisseph code then
  runs unchanged.
- All numeric returns are **native Python `float`** (never numpy scalars).

## Quick start

```bash
pip install libephemeris         # 3.1.0 stable. The wheel already bundles
                                  # the base-tier core: the 14 core bodies,
                                  # JD range 1850–2150. No download needed
                                  # for that slice.
```

```python
import libephemeris as swe
from libephemeris import SUN, MOON, MARS, FLG_SPEED

jd = swe.julday(2000, 1, 1, 12.0)        # UT1 Julian Day (J2000.0)

# Positions. calc_ut takes UT1; calc takes TT (Terrestrial/Ephemeris Time).
pos, retflag = swe.calc_ut(jd, SUN, FLG_SPEED)
lon, lat, dist, slon, slat, sdist = pos  # deg, deg, AU, deg/day, deg/day, AU/day

# House cusps — NOTE argument order: (jd_ut, lat, lon, hsys). LAT BEFORE LON.
# hsys accepts an int code (ord('P')) OR a bytes char (b'P'); both work.
cusps, ascmc = swe.houses(jd, 41.9028, 12.4964, b"P")   # Placidus, Rome
asc, mc = ascmc[0], ascmc[1]

# Extended houses: cusps, angles, and their daily speeds (always computed).
cusps, ascmc, cusps_speed, ascmc_speed = swe.houses_ex2(jd, 41.9028, 12.4964, ord("P"))
```

Return shapes (memorize — they match pyswisseph):

- `calc_ut(tjdut, body, flags)` / `calc(tjdet, body, flags)` →
  `((lon, lat, dist, speed_lon, speed_lat, speed_dist), retflag)`.
- `houses(...)` → `(cusps[12], ascmc[8])`. `ascmc` =
  `[Asc, MC, ARMC, Vertex, EquAsc, CoAsc, CoAscKoch, PolarAsc]`.
- `houses_ex2(...)` → `(cusps[12], ascmc[8], cusps_speed[12], ascmc_speed[8])`.
- `fixstar_ut(name, tjdut, flags)` → `(pos[6], resolved_name, retflag)`.
- `nod_aps_ut(tjdut, body, method=NODBIT_MEAN, flags)` →
  `(asc_node, desc_node, perihelion, aphelion)`, each a 6-float tuple.

### Setting where data lives

- **`LIBEPHEMERIS_DATA_DIR`** (env var) relocates the downloaded data root
  (default `~/.libephemeris`). This is where tiers/SPKs/IERS land.
- **`set_ephe_path(dir)`** — pyswisseph-compatible: point at a directory of
  ephemeris files; `.bsp` SPKs found there are auto-registered (no network).
- **`set_leb_file(path)`** or **`LIBEPHEMERIS_LEB`** — pin one specific `.leb`
  binary ephemeris for the fast path.

## The four calculation modes

Select with `set_calc_mode(mode)` or the `LIBEPHEMERIS_MODE` env var
(case-insensitive). Read the current one with `get_calc_mode()`. The **only**
valid values are:

| Mode | Backend | Notes |
|------|---------|-------|
| `"auto"` (default) | best-by-date LEB → non-LEB fallback | Resolves LEB first; if no local DE440, tries Horizons; else Skyfield. Best onboarding. |
| `"skyfield"` | NASA JPL DE440/DE441 **local**, via Skyfield | High-precision local JPL. Never reads LEB implicitly. |
| `"leb"` | Sealed LEB + declared local models | Offline, deterministic, source-pure. **Never** opens JPL/BSP, Horizons, ASSIST, SPK, or a socket. |
| `"horizons"` | NASA JPL Horizons REST API | Requires internet. Bodies/flags Horizons can't serve fall back to Skyfield. |

There is **no `"jpl"` mode** — passing `"jpl"` raises `ValueError`. "Local JPL"
is the `"skyfield"` mode; "remote JPL" is the `"horizons"` mode.

`mode="skyfield"` and `mode="horizons"` **never select LEB implicitly**, and
`mode="leb"` **never opens a JPL/BSP source**. See
`references/calc-modes.md` for the sealed-mode contract, the network policy
(`auto` / `allow` / `sealed`), and computation tracing.

## Sealed LEB mode (offline, deterministic) — rc14 boundary

`mode="leb"` is a real source boundary. Persisted states come **only** from the
active LEB files. A missing or corrupt canonical group is a **provisioning
error**; a core request outside every available LEB interval raises
`EphemerisRangeError`. rc14's headline change: **no silent substitution to a
lower-precision or different source** on a miss — it fails loudly.

Network policy is orthogonal and enforced at the socket choke point:

- Policies: `"auto"` (default), `"allow"`, `"sealed"`. Set via
  `set_network_policy(...)`, `LIBEPHEMERIS_NETWORK_POLICY`, or TOML
  `network_policy`.
- Under `"auto"`: `leb` mode is **sealed**, every other mode is **allowed**.
- When sealed, any network attempt raises `NetworkSealedError` (re-exported at
  the package root; subclass of `ConfigurationError` and `RuntimeError`).

## Data tiers and per-body coverage

Three **cumulative** precision tiers (`set_precision_tier(...)` /
`LIBEPHEMERIS_PRECISION`):

| Tier | DE kernel | Range | Size |
|------|-----------|-------|------|
| `base` | DE440s | 1850–2150 | ~31 MB |
| `medium` (**default**) | DE440 | 1550–2650 | ~114 MB |
| `extended` | DE441 | ≈ −13200 to +17191 | ~3.1 GB |

Resolution is **best-by-date, per body and per date**, in fixed priority
`base → medium → extended`: a narrower, higher-priority artifact is used where
it covers the request; a broader tier is consulted only outside that stored
interval. Per-body ranges recorded in the companion files are authoritative.

Ask the library rather than assuming coverage:

```python
from libephemeris import get_body_coverage, CHIRON
cov = get_body_coverage(CHIRON, jd)     # BodyCoverage | None
# None => the active LEB reader does NOT contain this body. It does NOT imply
#         that any analytical or online fallback exists.
if cov:
    cov.jd_start, cov.jd_end, cov.data_file, cov.group
    cov.precision_class            # e.g. "ephemeris", "analytical", "numerical-model"
    cov.reviewed                   # manifest-verified?
    cov.contains(jd)               # bool
```

`get_leb_inventory()` returns active files, per-body ranges, mode, and effective
network policy — the fastest self-diagnostic. See
`references/tiers-and-coverage.md` for the `data-v3` artifact set, provisioning,
and the full typed error contract.

### Provisioning

```bash
libephemeris init          # optional interactive config -> libephemeris-config.toml
libephemeris download auto # what the active mode/tier needs (see the caveat below)
libephemeris download leb2-base|leb2-medium|leb2-extended
                           # the SHA-256-pinned LEB2 groups for one tier,
                           # cumulative: 5 files (base), 10 (medium), 15 (extended)
libephemeris status        # verify installed data + active config (add --json)
```

**`download auto` is mode-dependent.** Under `mode = "leb"` it fetches only the
5/10/15 LEB2 files for the configured tier. Under the default `mode = "auto"`
it *also* downloads the DE kernel, `planet_centers.bsp` and the minor-body SPKs
— multi-GB at `extended`. Set the mode before running it in a sealed
deployment.

**Do not use `libephemeris download base|medium|extended` to provision LEB2.**
Those subcommands download DE kernels + SPKs and install **zero** LEB2 groups,
so a sealed-`leb` runtime fetches hundreds of MB it will never open and still
fails provisioning. The LEB2 commands are the `leb2-` prefixed ones above.

In `leb` mode you do **not** need `planet_centers.bsp`: outer-planet channels
store system barycentres directly; only a non-sealed runtime may apply a JPL
centre offset on top.

## Traps to avoid

- **`houses` arg order is `(jd, lat, lon, hsys)`** — latitude before longitude.
  Swapping them silently produces a wrong-but-plausible chart.
- **`calc` vs `calc_ut`**: `calc_ut` takes UT1, `calc` takes TT. Astrology work
  almost always wants `calc_ut`. Mixing them shifts positions by ΔT.
- **A body absent from your custom `.leb` (or outside its stored range) raises a
  typed error** (`UnknownBodyError` for an unknown id; `EphemerisRangeError` for
  out-of-range) — it is not silently zero-filled or downgraded. Check first with
  `get_body_coverage`.
- **`FLG_MOSEPH` / `FLG_SWIEPH` / `FLG_JPLEPH`** are accepted for compatibility;
  the ephemeris source is governed by the **calc mode**, not these flags
  (`FLG_MOSEPH` is effectively ignored — calculations always use JPL data).
- **`LEBCorruptionError`** (a `ValueError` subclass) signals a truncated/corrupt
  LEB file; it is treated as fatal provisioning failure and **never** triggers a
  source fallback. It is *not* re-exported at the package root — import it from
  `libephemeris.exceptions` if you must catch it distinctly, else it surfaces as
  `ValueError`.
- **Never** seek Swiss Ephemeris files or data as a source; libephemeris is
  self-contained on JPL/IAU/ERFA. Just install and provision as above.

## References

- `references/calc-modes.md` — the four modes in depth, sealed-mode contract,
  network policy, and computation tracing.
- `references/tiers-and-coverage.md` — tiers, best-by-date routing, `data-v3`
  artifact set, `get_body_coverage` / `BodyCoverage`, `get_leb_inventory`,
  provisioning, and the full error contract.
- `references/reference-api-compat.md` — drop-in migration, naming, argument-order
  gotchas, return shapes, flags, and the documented intentional divergences.
