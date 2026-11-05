# Data and configuration

## Choose a calculation route

| Mode | Uses |
| --- | --- |
| `leb` | Local polynomial coefficients and supported local models; no implicit JPL/Horizons substitution |
| `skyfield` | Native Rust JPL/SPK reader; the historical name does not load Python Skyfield |
| `auto` | Configured resources and provider routing |
| `horizons` | Explicit JPL Horizons network requests and supported routing |

Set Python defaults before creating contexts. A mode selects a route, not a
precision guarantee. Missing or corrupt explicit files are errors in forced
LEB mode; automatic routing does not make every body available.

## Included and optional data

The wheel includes the data-v4 base core LEB2 image for modern dates, with
stored core coverage JD 2396757.5–2506332.5. Additional body groups and tiers
use the same published data-v4 inventory; see the
[authenticated release manifest](../core/data/leb2/data-v4-manifest.json).
JPL kernels, IERS updates and N-body force files are separate. Planet-center
SPK downloads retain their independent data-v2rc1 release.

| Tier | Default JPL kernel | Nominal years |
| --- | --- | --- |
| `base` | `de440s.bsp` | 1850–2150 |
| `medium` | `de440.bsp` | 1550–2650 |
| `extended` | `de441.bsp` | −13200–17191 |

Check **actual per-file/body coverage**, especially for satellites, minor
bodies and event searches needing samples around the requested epoch.
Provision only the data you need:

```sh
libephemeris init
libephemeris download auto
libephemeris status --json
libephemeris leb info /path/to/base_core.leb2 --json
```

Downloads are explicit operations; `download all --yes` can be large. See
[CLI](cli.md) for targeted acquisition. JPL SPK segments also define their own
centres and frames; system-barycentre data are not automatically planet-centre
positions. A LEB image stores its own coordinates and model choices; see the
[binary format](data-format.md).

## Configuration

Precedence: explicit settings → environment → `.env` → TOML → defaults.
`reload_configuration()` reacquires captured settings. Existing Sessions and
retained handles keep their selected resources.

The conventional data root is `.libephemeris` under an explicit nonempty
`HOME`, or under the native user profile when `HOME` is absent. On Windows,
this uses `USERPROFILE` and the native profile API, not the working directory.
`LIBEPHEMERIS_DATA_DIR` remains the explicit data-directory override.

| Variable | Selects |
| --- | --- |
| `LIBEPHEMERIS_DATA_DIR` | Data root |
| `LIBEPHEMERIS_MODE` | `leb`, `skyfield`, `auto`, `horizons` |
| `LIBEPHEMERIS_PRECISION` | `base`, `medium`, `extended` |
| `LIBEPHEMERIS_LEB` | Explicit LEB path |
| `LIBEPHEMERIS_EPHEMERIS` | JPL kernel |
| `LIBEPHEMERIS_NETWORK_POLICY` | `auto`, `allow`, `sealed` |
| `LIBEPHEMERIS_IERS_DELTA_T` | IERS clock selection |
| `LIBEPHEMERIS_AUTO_SPK` | Automatic minor-body SPK acquisition |
| `LIBEPHEMERIS_STRICT_PRECISION` | Refuse approximate fallback |
| `LIBEPHEMERIS_LEB_READER` | `owned` or `indexed` storage |

For deterministic offline work, set the network policy to `sealed`, select
local resources and choose whether to use local IERS data. The
[installation example](installation.md#first-python-calculation) uses only the
bundled core and disables optional IERS acquisition. An explicit path is
validated when consumed, not necessarily when set.

`sealed` controls calculation-time acquisition; deliberately invoking a
provisioning command requests a download. `allow` permits acquisition but does
not itself start it. Avoid global configuration changes during concurrent work.

## Large LEB files

`owned` retains a complete immutable snapshot. `indexed` keeps a read-only file
descriptor, indexes and bounded decoded chunks:

```sh
export LIBEPHEMERIS_LEB_READER=indexed
```

Optional TOML limits (defaults shown):

```toml
[libephemeris]
leb_reader = "indexed"

[libephemeris.leb_limits]
metadata_bytes = 8388608
cache_bytes = 8388608
working_bytes = 67108864
compressed_bytes = 4194304
frame_bytes = 8388608
window_log_max = 23
```

The limits bound admitted metadata, retained coefficients, live decode work,
compressed buffers, decoded frames and zstd window size respectively. Window
logarithms must be 10–30. Zero cache bytes disables retention, not evaluation.
An insufficient limit raises a resource/configuration error; it does not
silently switch data or reduce precision. Legacy monolithic LEB2 v1 bodies can
need larger frame limits than chunked v2.

Keep an acquired inode immutable: replace a pathname atomically instead of
rewriting or truncating its open file. Use `owned` if in-place edits cannot be
excluded. Clones share ownership; the last owner releases the descriptor.
Indexed `bytes()`/borrowed contiguous reads are unavailable rather than
materializing the entire image. `cool()` releases retained coefficients, not
active borrowers; `warm()` does not load the whole indexed file.

Related runtime images share allocation-domain counters; do not sum them once
per image. Independent opens need not share a budget. These limits do **not**
bound process RSS, Python/SPK allocations or the OS cache. The static bundled
image is already zero-copy and is not converted into an indexed file.
