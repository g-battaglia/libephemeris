# Computation Tracing

LibEphemeris can record which sub-backend computed each celestial body during a
calculation session. This is useful for debugging, performance analysis, and
understanding the auto-mode fallback chain.

## Two Tracing Mechanisms

LibEphemeris offers two complementary ways to trace computation sources:

| Mechanism | API | Overhead | Output | Best for |
|-----------|-----|----------|--------|----------|
| **Programmatic tracing** | `start_tracing()` / `get_trace_results()` | ~50 ns when inactive | `{body_id: "source"}` dict | Application code, automated checks |
| **DEBUG logging** | `LIBEPHEMERIS_LOG_LEVEL=DEBUG` | Log I/O per call | Structured log lines | Manual debugging, test runs |

## Programmatic Tracing (ContextVar-based)

The recommended approach for application code. Uses Python `ContextVar` for
thread-safe, zero-overhead-when-inactive tracing.

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

Each successful computation records one of these source tags:

| Source | Description |
|--------|-------------|
| `"LEB"` | Precomputed Chebyshev polynomials (fastest path) |
| `"Skyfield"` | NASA JPL DE440/DE441 via Skyfield |
| `"Horizons"` | NASA JPL Horizons REST API |
| `"SPK"` | Direct SPK kernel evaluation (minor bodies) |
| `"ASSIST"` | N-body integration via REBOUND/ASSIST |
| `"Keplerian"` | Analytical Keplerian orbit (last-resort fallback) |

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

When tracing is **not active**, the overhead is a single
`ContextVar.get(None)` check per computation (~50 ns). When active, each
dispatch point adds one dict assignment (~100 ns). In both cases the cost
is negligible compared to ephemeris calculations (~5-120 us).

## DEBUG Log Tracing

For manual debugging and test runs, enable DEBUG-level logging to see
per-call source information in structured log lines:

```bash
LIBEPHEMERIS_LOG_LEVEL=DEBUG pytest -s tests/
```

Or programmatically:

```python
import logging
logging.getLogger("libephemeris").setLevel(logging.DEBUG)
```

Log output format:

```
[libephemeris] DEBUG: body=0 jd=2448045.9167 source=LEB
[libephemeris] DEBUG: body=15 jd=2448045.9167 source=SPK
[libephemeris] DEBUG: body=146199 jd=2448045.9167 source=ASSIST (n-body)
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

Both mechanisms track the same dispatch points (`_record()` call sites in
`planets.py` and `context.py`).
