# Data tiers, coverage, provisioning, and the error contract

## The three cumulative tiers

`set_precision_tier(name)` / `LIBEPHEMERIS_PRECISION` / TOML `precision`. Read
with `get_precision_tier()`; enumerate with `list_tiers()` (each entry is a
`PrecisionTier` dataclass: `name`, `ephemeris_file`, `spk_date_range`,
`description`).

| Tier | DE kernel | Approx. range | Size |
|------|-----------|---------------|------|
| `base` | `de440s.bsp` | 1850–2150 | ~31 MB |
| `medium` (**default**) | `de440.bsp` | 1550–2650 | ~114 MB |
| `extended` | `de441.bsp` | ≈ −13200 to +17191 | ~3.1 GB |

The tiers are **cumulative**. Selecting `extended` layers `base` and `medium`
underneath it, and resolution is **best-by-date, per body and per date** in the
fixed priority order `base → medium → extended`:

1. `base` — DE440s, 1850–2150
2. `medium` — DE440, 1550–2650
3. `extended` — DE441, its exact source interval

A narrow, higher-priority artifact is used wherever it covers the request; a
broader tier is consulted only for dates outside that stored interval. **Per-body
coverage recorded in the companion files is authoritative** — do not assume a
tier's headline range applies to every body (minor-body and exotic ranges vary,
and extended exotics are limited to a guarded common interval).

Changing tier invalidates cached kernels and the LEB reader so the next call
reloads the correct data.

## The `data-v3` artifact set (rc14)

Each tier ships **five** LEB2 groups:
`("core", "asteroids", "exotics", "apogee", "uranians")`. Three tiers × five
groups = **15 immutable, SHA-256-pinned files**. `base_core.leb2` and
`base_uranians.leb2` also ship inside the wheel and are byte-identical to the
GitHub release assets.

- Shared nutation, ΔT, and fixed-star tables live once in each `core`; named
  companions carry only their own body channels (lightweight modular design).
- Outer-planet channels store the **system barycentres** supplied by the DE
  kernel. Optional planet-centre BSPs and analytical centre-of-body models never
  change artifact bytes. In `leb` mode you therefore **do not need
  `planet_centers.bsp`**; only a non-sealed runtime may apply an independently
  sourced JPL centre offset after reading the barycentric LEB state.
- Extended exotic trajectories generated via ASSIST are explicitly
  `precision_class="numerical-model"`, not direct minor-body JPL ephemerides,
  and their stored range is the guarded DE441 ∩ perturber interval. Outside it,
  sealed mode reports the deterministic local model it actually used instead of
  relabeling it as LEB coverage.

## Querying coverage at runtime

`get_body_coverage(body_id, jd=None)` (alias: `coverage(...)`) →
`BodyCoverage | None`.

```python
from libephemeris import get_body_coverage, MARS
cov = get_body_coverage(MARS, jd)
```

`None` means **the active LEB reader does not contain the body**. It does **not**
imply that an analytical or online fallback exists — treat `None` as "not served
here", then decide explicitly (switch tier/mode, or handle the gap).

`BodyCoverage` is a frozen dataclass:

| Field | Meaning |
|-------|---------|
| `body_id` | the requested body id (int) |
| `source` | `"LEB"` |
| `precision_class` | numerical origin — see below |
| `jd_start`, `jd_end` | closed JD interval that serves this body/date |
| `data_file` | the concrete selected file, or `None` when no candidate covers `jd` (then start/end are the union range) |
| `group` | canonical group inferred from the file (`core`, `asteroids`, …) |
| `reviewed` | whether the serving file is manifest-verified |

Methods: `cov.contains(jd)` → bool; `cov.to_dict()` → JSON-serializable.

`precision_class` values you will see: `"ephemeris"` (direct JPL-derived),
`"analytical"` (nodes/apsides/Lilith, hypothetical IDs 40–58), `"numerical-model"`
(ASSIST-integrated extended exotics), `"mixed"` (date-less union spanning classes),
`"unverified-local"` (not manifest-reviewed).

### Whole-runtime introspection

- `get_leb_inventory()` → dict with `mode`, `precision_tier`,
  `network_policy_configured`, `network_policy_effective`, `ready`,
  `reader_type`, per-file `files` (name/path/group/size/reviewed/body ranges),
  and `body_count`. The fastest single self-diagnostic.
- `get_runtime_data_requirements(tier=None)` → the manifest-derived file
  contract (`RuntimeDataRequirement` items: `name`, `kind`, `group`, `path`,
  `sha256`) that a sealed tier needs. A missing pin raises immediately rather
  than reporting partial readiness.
- `inspect_leb_file(path)` → body/range metadata for one file without activating
  it.

## Provisioning

The wheel already contains the base-tier core (14 core bodies) plus the Hamburg
(Uranian) companion for 1850–2150 — usable with zero downloads. For anything
else:

```bash
libephemeris init                 # optional interactive wizard -> libephemeris-config.toml
libephemeris download auto        # exactly what the configured mode/tier needs
libephemeris download base        # or a specific tier ...
libephemeris download medium
libephemeris download extended
libephemeris status               # verify installed data + active config
libephemeris status --json        # machine-readable
```

`download auto` installs every SHA-256-pinned LEB2 group cumulatively through
the configured tier: **5 files for base, 10 through medium, 15 through extended**.
The wheel supplies the two bundled base resources; the rest come from the
immutable `data-v3` release. Flags on every download command: `--force`,
`--no-progress`, `--quiet`. Relocate the data root with `LIBEPHEMERIS_DATA_DIR`
(default `~/.libephemeris`).

A generated `libephemeris-config.toml` (committable) captures the same settings:

```toml
[libephemeris]
precision = "medium"          # base | medium | extended
mode = "leb"                  # auto | skyfield | leb | horizons
network_policy = "sealed"     # auto | allow | sealed  (leb defaults to sealed)
auto_spk = true
strict_precision = true
```

## The typed error contract (rc14)

rc14's headline: **no silent substitution to a lower-precision or different
source on a miss.** Failures are typed and loud. All inherit from
`libephemeris.Error` (itself an `Exception`), so `except swe.Error` still catches
everything the reference API would.

| Exception | Raised when | Key attributes |
|-----------|-------------|----------------|
| `EphemerisRangeError` | requested JD is outside the ephemeris/LEB coverage (or non-finite) | `requested_jd`, `start_jd`, `end_jd`, `start_date`, `end_date`, `body_id`, `body_name`, `ephemeris_file` |
| `UnknownBodyError` | unknown/unsupported body id (rather than silently returning zeros) | `body_id` |
| `NetworkSealedError` | any network attempt while the effective policy is `sealed` | `purpose`; also subclasses `ConfigurationError` + `RuntimeError` |
| `SPKRequiredError` | strict precision on a major asteroid with no SPK registered | `body_id`, `body_name` |
| `PolarCircleError` | Placidus/Koch/Gauquelin house beyond the polar circle | `latitude`, `threshold`, `obliquity`, `house_system` |
| `LEBCorruptionError` | a LEB/LEB2 file is truncated/corrupt | — subclass of `ValueError`; **not** re-exported at the package root |

Notes:

- `LEBCorruptionError` is an internal robustness signal. Dispatchers catch it
  first and treat it as a **fatal provisioning error** — corruption is never
  confused with an ordinary range miss and never triggers a source fallback. To
  catch it distinctly, `from libephemeris.exceptions import LEBCorruptionError`;
  otherwise it surfaces as a `ValueError`.
- In `leb` mode, JD-range validation uses the **stored per-body LEB header** as
  authority — no JPL kernel is opened to validate a body it will not serve.
- A body that is absent from the active LEB reader (or a custom `.leb` you
  pinned) and has no traced local model raises a typed error — check
  `get_body_coverage(...)` first to branch cleanly.

```python
import libephemeris as swe

try:
    pos, _ = swe.calc_ut(jd, body, swe.FLG_SPEED)
except swe.EphemerisRangeError as e:
    # e.start_date / e.end_date tell you the served window; widen tier or clamp.
    ...
except swe.UnknownBodyError as e:
    # bad/unsupported body id: e.body_id
    ...
except swe.Error:
    # any other ephemeris-level failure
    ...
```
