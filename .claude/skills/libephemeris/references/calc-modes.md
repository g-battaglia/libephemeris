# Calculation modes, sealed LEB, and network policy

Two independent switches govern *where numbers come from* and *whether the
process may touch the network*:

1. **Calc mode** — `set_calc_mode()` / `LIBEPHEMERIS_MODE` / TOML `mode`.
2. **Network policy** — `set_network_policy()` / `LIBEPHEMERIS_NETWORK_POLICY` /
   TOML `network_policy`.

They interact but are not the same knob. Understand both before you assume a
run is offline.

## The four calc modes

`set_calc_mode(mode)` accepts exactly `{"auto", "skyfield", "leb", "horizons"}`
(case-insensitive; leading/trailing space trimmed). Anything else raises
`ValueError`. `get_calc_mode()` resolves in this order: programmatic override →
`LIBEPHEMERIS_MODE` → TOML `mode` → default `"auto"`.

### `auto` (default)

Prefers the best-by-date LEB chain first, then keeps the configured non-LEB
fallback: if no local DE440 is available it tries the Horizons API, otherwise
Skyfield. Best for onboarding — it resolves local or remote data transparently.
Because the fallback can reach the network, `auto` is **not** guaranteed
offline unless you also seal the network policy or fully provision LEB data.

### `skyfield`

Always computes from local NASA JPL DE440/DE441 via Skyfield (~120 µs/call).
High-precision local JPL workflow. **Never reads LEB implicitly.** Fails if the
required DE kernel is not present locally.

### `leb` (sealed)

Sealed, offline, deterministic, source-pure (~5 µs/call). Persisted ephemeris
states come only from the active LEB files. In this mode the runtime **does not
open** DE/BSP kernels, planet-center subsets, registered or automatic SPKs,
Horizons, ASSIST, or a network socket.

- A missing or corrupt canonical group → **provisioning error** (see the error
  contract in `tiers-and-coverage.md`).
- A core request outside every available LEB interval → `EphemerisRangeError`.
- Curated bodies with **no meaningful LEB channel** (lunar nodes/apsides, Black
  Moon Lilith variants, historical hypothetical bodies) may use a **traced local
  analytical or Keplerian model**. Tracing labels that as its real source — it is
  never mislabeled as LEB.
- Calculations that internally need Skyfield's vector protocol (topocentric
  fixed stars, `calc_pctr`, `nod_aps`, related event paths) receive
  `LEBVectorEphemeris`, an in-memory adapter over the same LEB reader, so those
  paths stay sealed without duplicating astronomy.

### `horizons`

Prefers the NASA JPL Horizons REST API (~300 ms/call, requires internet). Bodies
or flags Horizons cannot serve fall back to Skyfield. Zero local ephemeris files
required.

> There is no `"jpl"` mode. If a spec says "JPL", map it to `skyfield` (local
> DE440/DE441) or `horizons` (remote Horizons) depending on whether you want a
> local kernel or an API call.

### rc14 compatibility guarantees

- `auto` prefers the same best-by-date LEB chain, then its configured non-LEB
  fallback.
- `skyfield` and `horizons` **do not select LEB implicitly**.
- `leb` **never opens a JPL/BSP source**.
- Public calculation signatures and the source/precision tracing vocabulary are
  unchanged.

## Network policy

The policy is the single choke point for outbound HTTPS (SPK generation, SBDB
lookups, IERS refreshes, Skyfield's implicit kernel download, Horizons).

- Values: `"auto"` (default), `"allow"`, `"sealed"`. Invalid spellings raise
  `ValueError` — a typo never silently widens the boundary.
- Resolution: `set_network_policy()` → `LIBEPHEMERIS_NETWORK_POLICY` → TOML
  `network_policy` → default `"auto"`.
- Under `"auto"`: **`leb` mode is sealed; every other mode is allowed.**
- `get_configured_network_policy()` returns the raw setting (may be `"auto"`);
  `get_network_policy()` returns the **effective** `"allow"` / `"sealed"` for the
  current mode.

When effective policy is `sealed`, the first network attempt raises
`NetworkSealedError` before any socket work:

```python
from libephemeris import NetworkSealedError
try:
    ...  # some path that would download
except NetworkSealedError as e:
    print(e.purpose)   # what was blocked
```

`NetworkSealedError` subclasses both `ConfigurationError` and `RuntimeError`, so
existing `except RuntimeError` downloader fallbacks still catch it; catch
`NetworkSealedError` or `ConfigurationError` to distinguish policy from transport
failure.

The `libephemeris download ...` CLI deliberately runs with `"allow"`: a download
is an explicit provisioning action, not calculation-time egress.

### Recipes

```python
import libephemeris as swe

# Fully offline & deterministic, fail loudly if data is missing:
swe.set_calc_mode("leb")            # -> effective network policy becomes 'sealed'

# Local high-precision JPL, but forbid any network anyway:
swe.set_calc_mode("skyfield")
swe.set_network_policy("sealed")

# Onboarding convenience, network allowed to fetch what's missing:
swe.set_calc_mode("auto")           # default; network_policy 'auto' -> allowed
```

## Computation tracing — which backend served each body

Use tracing to verify the source actually used (crucial when auditing a sealed
deployment or debugging a fallback):

```python
import libephemeris as swe
swe.start_tracing()
swe.calc_ut(jd, swe.SUN, swe.FLG_SPEED)
results = swe.get_trace_results()   # per-body source/precision records
```

This is how you confirm a body came from LEB vs. a traced analytical model vs.
Skyfield — the tracing vocabulary never labels an analytical/Keplerian model as
LEB.
