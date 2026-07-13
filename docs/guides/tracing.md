# Computation Tracing

LibEphemeris can record which sub-backend computed each celestial body during a
calculation session. This is useful for debugging, performance analysis, and
understanding the auto-mode fallback chain.

## Two Tracing Mechanisms

LibEphemeris offers two complementary ways to trace computation sources:

| Mechanism | API | Overhead | Output | Best for |
|-----------|-----|----------|--------|----------|
| **Programmatic tracing** | `start_tracing()` / `get_trace_results()` | A few context-local operations on fallback paths | `{body_id: "source"}` dict | Application code, automated checks |
| **DEBUG logging** | `LIBEPHEMERIS_LOG_LEVEL=DEBUG` | Log I/O per call | Structured log lines | Manual debugging, test runs |

## Programmatic Tracing (ContextVar-based)

The recommended approach for application code. Uses Python `ContextVar` for
thread-safe, low-overhead tracing.

The trace contract covers module-level and `EphemerisContext` `calc()` /
`calc_ut()` calls, the four `fixstar*()` entry points, and
`batch_fixstars_ut()`. Planetary moons are traced when calculated through
`calc()` or `calc_ut()`. It does not currently promise a direct entry for
`calc_pctr()` or for every utility that performs internal calculations.

### Basic Usage

```python
import libephemeris as swe
from libephemeris.constants import SUN, MOON, MARS, FLG_SPEED

# Start tracing -- returns a token for cleanup
token = swe.start_tracing()

# Compute positions as usual
swe.calc_ut(2451545.0, SUN, FLG_SPEED)
swe.calc_ut(2451545.0, MOON, FLG_SPEED)
swe.calc_ut(2451545.0, MARS, FLG_SPEED)

# Retrieve results: {body_id: "source_name"}
traces = swe.get_trace_results()
print(traces)
# Example: {0: "LEB", 1: "LEB", 4: "Skyfield"}

# Clean up
token.var.reset(token)
```

### Traced Backends

Each successfully returned body on the traced API surface records one of these
source tags:

| Source | Description |
|--------|-------------|
| `"LEB"` | Precomputed Chebyshev polynomials (fastest path) |
| `"Skyfield"` | Skyfield/JPL pipeline for major bodies, fixed stars, or JPL-derived lunar points |
| `"Horizons"` | NASA JPL Horizons REST API |
| `"SPK"` | Direct SPK kernel evaluation (minor bodies) |
| `"ASSIST"` | N-body integration via REBOUND/ASSIST |
| `"Keplerian"` | Analytical Keplerian orbit (last-resort fallback) |
| `"Analytical"` | Local formula or zero-origin convention for lunar points, hypothetical bodies, angles, Arabic parts, heliocentric Sun, or geocentric Earth |
| `"ERFA"` | IAU 2006/2000A nutation used for the `ECL_NUT` pseudo-body |
| `"Mixed"` | A fixed-star velocity stencil that crossed from LEB to Skyfield fallback |

### Thread Safety

`start_tracing()` uses Python's `contextvars.ContextVar`, so each thread
gets its own independent trace accumulator:

```python
import threading
import libephemeris as swe
from libephemeris.constants import SUN, MOON, FLG_SPEED

results = {}

def worker(name, body):
    token = swe.start_tracing()
    swe.calc_ut(2451545.0, body, FLG_SPEED)
    results[name] = swe.get_trace_results()
    token.var.reset(token)

t1 = threading.Thread(target=worker, args=("sun", SUN))
t2 = threading.Thread(target=worker, args=("moon", MOON))
t1.start(); t2.start()
t1.join(); t2.join()

# Each thread only sees its own body
print(results["sun"])   # {0: "LEB"}
print(results["moon"])  # {1: "LEB"}
```

### Nested Tracing

Calling `start_tracing()` inside an already-active session creates an
independent scope. Resetting the inner token restores the outer session:

```python
token_outer = swe.start_tracing()
swe.calc_ut(2451545.0, SUN, FLG_SPEED)

token_inner = swe.start_tracing()
swe.calc_ut(2451545.0, MOON, FLG_SPEED)
print(swe.get_trace_results())  # {1: "LEB"} -- inner only
token_inner.var.reset(token_inner)

print(swe.get_trace_results())  # {0: "LEB"} -- outer restored
token_outer.var.reset(token_outer)
```

### Overwrite Behavior

If the same body is computed multiple times, the last source wins:

```python
token = swe.start_tracing()
swe.calc_ut(jd1, SUN, FLG_SPEED)  # computed via LEB
swe.calc_ut(jd2, SUN, FLG_SPEED)  # computed via Skyfield
traces = swe.get_trace_results()
print(traces[SUN])  # "Skyfield" (last call wins)
token.var.reset(token)
```

### Return Value Safety

`get_trace_results()` always returns a **copy** of the internal dict.
Mutating the returned dict does not affect the tracing session.

### Performance

When tracing is **not active**, a normal `calc()` / `calc_ut()` fallback checks
the trace `ContextVar` once. It also opens a short-lived context-local source
hint scope (one set/get/reset); precise inner branches may update that hint with
SPK, ASSIST, Keplerian, planetary-moon, or fixed-star metadata. This defers
publication until public output processing has succeeded. Direct LEB/Horizons
branches instead check the nested-log capture context once and call the inactive
recorder once. Fixed-star calculations additionally keep lightweight
context-local metadata to distinguish LEB, Skyfield, and mixed velocity
stencils. When tracing is active, successful public dispatch adds a dictionary
assignment. No universal nanosecond figure is claimed because it depends on the
entry point, backend, Python version, and platform.

## DEBUG Log Tracing

For manual debugging and test runs, enable DEBUG-level logging to see
per-call source information in structured log lines:

```bash
LIBEPHEMERIS_LOG_LEVEL=DEBUG uv run pytest -s tests/test_tracing.py
```

Or programmatically:

```python
import logging
import libephemeris as swe

swe.set_log_level(logging.DEBUG)
```

Log output format:

```
[libephemeris] DEBUG: body=0 jd=2448045.9 source=LEB
[libephemeris] DEBUG: body=15 jd=2448045.9 source=SPK
[libephemeris] DEBUG: body=146199 jd=2448045.9 source=ASSIST
```

See [Testing -- Backend Isolation](../development/testing.md#backend-isolation)
for guidance on source-selection assertions.

## When to Use Which

- **Programmatic tracing** (`start_tracing()`): when you need structured data
  in application code -- e.g., to verify that LEB is being used, to build
  diagnostic endpoints, or to assert backend selection in tests.

- **DEBUG logging**: when you need human-readable output during manual debugging
  or CI test runs. More verbose (includes Julian Day per call), but requires
  log parsing to extract structured data.

For the traced API surface, both mechanisms publish the same final source under
the caller-facing body ID. Recursive alias and fixed-epoch work is held private
until the outer call's transformations and output conversion have succeeded.
