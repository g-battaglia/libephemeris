# Testing Guide

This document describes the test infrastructure, verification methodology, and test results for LibEphemeris.

## Table of Contents

- [Test Results Overview](#test-results-overview)
- [Expected Failures (xfail)](#expected-failures-xfail)
- [Running All Tests (Without Skips)](#running-all-tests-without-skips)
- [Full Test Command](#full-test-command)
- [Categories of Skipped Tests](#categories-of-skipped-tests)
- [Considerations](#considerations)
- [Quick Reference](#quick-reference)
- [Testing with SPK Kernels (Downstream Projects)](#testing-with-spk-kernels-downstream-projects)

## Test Results Overview

When running `pytest`, the output may show results like:

```
6336 passed
   4 xpassed
   3 xfailed
  57 skipped
```

### What do these mean?

| Status | Meaning |
|--------|---------|
| **passed** | Test completed successfully |
| **xpassed** | Test was expected to fail, but passed (unexpected success) |
| **xfailed** | Test was expected to fail and did fail (expected behavior) |
| **skipped** | Test was intentionally skipped due to missing dependencies |

> **Note:** `xfailed` tests are NOT failures! They document known limitations or future goals.

---

## Expected Failures (xfail)

These tests are marked as "expected to fail" for documented reasons:

| Test | Reason |
|------|--------|
| `test_deltat_future_years_match_swe` | Delta T for future dates (>2030) is a prediction - different algorithms give divergent results |
| `test_interpolated_apogee_target_precision` | Target precision <0.1° requires implementing analytical Moshier method |
| `test_interpolated_perigee_target_precision` | Target precision <0.1° requires implementing analytical Moshier method |
| Benchmark tests | Performance comparisons between pure Python and C (pyswisseph) - Python is inherently slower |

---

## Running All Tests (Without Skips)

To run the complete test suite without skips, the following conditions must be satisfied.

### 1. Install pyswisseph

Most skipped tests require the C-based Swiss Ephemeris library for comparison:

```bash
pip install pyswisseph
```

### 2. Install Optional Dependencies

```bash
pip install requests  # Required for orbital elements update tests
```

### 3. Enable Network Tests

Some tests require network access to download data from JPL Horizons. Enable them with environment variables:

```bash
# Enable SPK auto-download tests
export LIBEPHEMERIS_TEST_SPK_AUTO_DOWNLOAD=1

# Enable SPK download tests
export LIBEPHEMERIS_TEST_SPK_DOWNLOAD=1
```

### 4. Ephemeris Data Files

Some tests require reference ephemeris data files. Point the
``REF_EPHE_PATH`` environment variable to a directory containing
the expected datasets:

- `fictitious_orbits.csv` - bundled orbital elements dataset (included in package)
- Asteroid/TNO data files
- `seorbel.txt` - orbital elements reference (optional, skipped if missing)

---

## Full Test Command

Run all tests with network tests enabled:

```bash
LIBEPHEMERIS_TEST_SPK_AUTO_DOWNLOAD=1 \
LIBEPHEMERIS_TEST_SPK_DOWNLOAD=1 \
pytest
```

Or set the variables permanently in the shell:

```bash
export LIBEPHEMERIS_TEST_SPK_AUTO_DOWNLOAD=1
export LIBEPHEMERIS_TEST_SPK_DOWNLOAD=1
pytest
```

---

## Categories of Skipped Tests

### pyswisseph-dependent (majority of skips)

Tests that compare LibEphemeris results against pyswisseph:

- `test_lunar/test_true_lilith_latitude.py`
- `test_lunar/test_true_lilith_precision.py`
- `test_minor_bodies/test_tno_validation.py`
- `test_sol_eclipse_where_how.py`
- And many more...

### Network-dependent

Tests that download data from external sources:

| Test File | Environment Variable |
|-----------|---------------------|
| `test_spk_auto.py` | `LIBEPHEMERIS_TEST_SPK_AUTO_DOWNLOAD=1` |
| `test_spk.py` | `LIBEPHEMERIS_TEST_SPK_DOWNLOAD=1` |
| `test_asteroid_by_number.py` | Always skipped (requires JPL SBDB) |

### Data-dependent (runtime skips)

These tests skip dynamically based on available data:

- Eclipse tests: skip if no eclipse found in search period
- Star tests: skip if star not in catalog
- TNO tests: skip if SwissEph data not available
- Historical tests: skip if date outside ephemeris range

---

## Considerations

1. **Network tests are slow** - they download data from JPL Horizons
2. **Some skips are unavoidable** - dynamic skips depend on runtime data availability
3. **pyswisseph versions vary** - some versions don't include eclipse functions
4. **Benchmark xfails are expected** - pure Python will always be slower than C

---

## Quick Reference

| Goal | Command |
|------|---------|
| Run all tests | `pytest` |
| Run with verbose output | `pytest -v` |
| Run with network tests | `LIBEPHEMERIS_TEST_SPK_AUTO_DOWNLOAD=1 LIBEPHEMERIS_TEST_SPK_DOWNLOAD=1 pytest` |
| Run specific test file | `pytest tests/test_file.py` |
| Run tests matching pattern | `pytest -k "pattern"` |
| Show skip reasons | `pytest -rs` |
| Show xfail reasons | `pytest -rx` |

---

## Testing with SPK Kernels (Downstream Projects)

Projects that depend on LibEphemeris (e.g., kerykeion) often need high-precision asteroid calculations using SPK kernels. This section provides ready-to-use pytest fixtures for configuring SPK auto-download in a test suite.

### The Problem

By default, LibEphemeris uses Keplerian orbital elements for minor bodies like Chiron, Ceres, Vesta, etc. This provides ~10-30 arcsecond accuracy. For sub-arcsecond precision, SPK kernels must be downloaded from JPL Horizons.

Without proper test configuration:
- First test run downloads SPK files (slow, network-dependent)
- Each test might trigger redundant download checks
- CI environments may fail if network is unavailable
- Test isolation becomes difficult

### Recommended Fixture Pattern

Add this fixture to the project's `conftest.py`:

```python
import pytest
import libephemeris as eph
from libephemeris.constants import CHIRON, CERES, PALLAS, JUNO, VESTA


@pytest.fixture(scope="session", autouse=True)
def setup_spk_for_tests():
    """
    Session-scoped fixture to pre-download SPK kernels for major asteroids.
    
    This ensures all SPK files are downloaded once at the start of the test
    session, rather than on-demand during individual tests. This improves
    test speed and reliability.
    """
    # Enable automatic SPK download
    eph.set_auto_spk_download(True)
    
    # Pre-download SPK kernels for commonly used major asteroids
    # These are the bodies with automatic SPK download support
    major_asteroids = [
        CHIRON,   # Centaur, commonly used in astrology
        CERES,    # Dwarf planet
        PALLAS,   # Main belt asteroid
        JUNO,     # Main belt asteroid
        VESTA,    # Main belt asteroid
    ]
    
    for body in major_asteroids:
        # ensure_major_asteroid_spk returns True if SPK is available
        # (either already cached or successfully downloaded)
        eph.ensure_major_asteroid_spk(body)
    
    yield
    
    # Optional: cleanup after all tests
    # eph.set_auto_spk_download(False)
```

### Minimal Fixture (Chiron Only)

For projects that only need Chiron (the most commonly used minor body):

```python
import pytest
import libephemeris as eph
from libephemeris.constants import CHIRON


@pytest.fixture(scope="session", autouse=True)
def setup_chiron_spk():
    """Pre-download Chiron SPK for high-precision calculations."""
    eph.set_auto_spk_download(True)
    eph.ensure_major_asteroid_spk(CHIRON)
    yield
```

### Environment Variables

LibEphemeris respects these environment variables for testing:

| Variable | Values | Description |
|----------|--------|-------------|
| `LIBEPHEMERIS_AUTO_SPK` | `1`, `true`, `yes` | Enable auto SPK download globally |
| `LIBEPHEMERIS_LOG_LEVEL` | `DEBUG`, `INFO`, `WARNING`, `ERROR` | Control logging verbosity |

### Source Tracing (DEBUG logs)

At `DEBUG` level, libephemeris logs which calculation backend was used for
every body at every dispatch point. The log format is:

```
body=<id> jd=<julian_day> source=<BACKEND>
```

Example output:

```
[libephemeris] DEBUG: body=0 jd=2448045.9167 source=LEB
[libephemeris] DEBUG: body=15 jd=2448045.9167 source=SPK
[libephemeris] DEBUG: body=146199 jd=2448045.9167 source=ASSIST (n-body)
```

Possible `source` values:

| Value | Description |
|-------|-------------|
| `LEB` | Precomputed `.leb` binary ephemeris (Chebyshev polynomial fast-path) |
| `LEB->fallback` | LEB lookup failed, retrying with Skyfield (the next log line shows the actual source) |
| `Skyfield` | NASA JPL DE440/DE441 via Skyfield |
| `Horizons` | NASA JPL Horizons online API |
| `SPK` | Direct SPK kernel evaluation (minor bodies) |
| `SPK (auto-downloaded)` | SPK kernel auto-downloaded for this body |
| `ASSIST (n-body)` | N-body integration via ASSIST |
| `Keplerian (fallback)` | Analytical Keplerian orbit (last resort, logged at WARNING) |

These log statements are emitted in `planets.py` (`calc_ut()`,
`calc()`, `_calc_body()`) and `context.py` (`calc_ut()`, `calc()`).

To enable source tracing during tests:

```bash
LIBEPHEMERIS_LOG_LEVEL=DEBUG pytest -s tests/
```

Or programmatically:

```python
import logging
logging.getLogger("libephemeris").setLevel(logging.DEBUG)
```

The Astrologer API uses this mechanism via an opt-in `X-Debug-Ephemeris: true`
HTTP header that temporarily lowers the logger to DEBUG, captures `source=`
entries, and injects them into the JSON response.

### Programmatic Tracing (ContextVar-based)

For structured, machine-readable source tracing in application code or tests,
use `start_tracing()` / `get_trace_results()` instead of DEBUG logging:

```python
import libephemeris as swe
from libephemeris.constants import SUN, FLG_SPEED

token = swe.start_tracing()
swe.calc_ut(2451545.0, SUN, FLG_SPEED)
traces = swe.get_trace_results()   # {0: "LEB"}
token.var.reset(token)
```

This returns a `{body_id: source_name}` dict and is thread-safe via
`ContextVar`. Overhead is ~50 ns when inactive. See the
[Tracing Guide](../guides/tracing.md) for full details, thread safety,
and nested tracing.

### Custom SPK Cache Directory

For custom storage locations:

```python
import pytest
import libephemeris as eph
from pathlib import Path


@pytest.fixture(scope="session", autouse=True)
def setup_spk_with_custom_cache(tmp_path_factory):
    """Setup SPK with custom cache directory."""
    # Use a persistent cache directory for CI
    cache_dir = Path.home() / ".cache" / "libephemeris-spk"
    cache_dir.mkdir(parents=True, exist_ok=True)
    
    eph.set_spk_cache_dir(str(cache_dir))
    eph.set_auto_spk_download(True)
    
    # Pre-download all major asteroids
    from libephemeris.constants import CHIRON, CERES, PALLAS, JUNO, VESTA
    for body in [CHIRON, CERES, PALLAS, JUNO, VESTA]:
        eph.ensure_major_asteroid_spk(body)
    
    yield
```

### Disabling Strict Precision Mode

LibEphemeris has a "strict precision" mode that raises errors when SPK kernels are not available for major asteroids. For tests that don't require SPK precision:

```python
import pytest
import libephemeris as eph


@pytest.fixture(autouse=True)
def disable_strict_precision():
    """Allow Keplerian fallback for minor body calculations."""
    original = eph.get_strict_precision()
    eph.set_strict_precision(False)
    yield
    eph.set_strict_precision(original if original else None)
```

### Complete Example for Downstream Projects

Here's a complete `conftest.py` for a project using libephemeris:

```python
"""
pytest configuration for a project using libephemeris.
"""
import pytest
import libephemeris as eph
from libephemeris.constants import (
    CHIRON, CERES, PALLAS, JUNO, VESTA,
    SIDM_FAGAN_BRADLEY,
)


@pytest.fixture(scope="session", autouse=True)
def setup_libephemeris():
    """
    Session-scoped setup for libephemeris.
    
    Downloads SPK kernels once at session start for consistent,
    high-precision asteroid calculations across all tests.
    """
    # Enable SPK auto-download
    eph.set_auto_spk_download(True)
    
    # Disable strict mode to allow Keplerian fallback if download fails
    eph.set_strict_precision(False)
    
    # Pre-download SPK for major asteroids
    for body in [CHIRON, CERES, PALLAS, JUNO, VESTA]:
        eph.ensure_major_asteroid_spk(body)
    
    yield


@pytest.fixture(autouse=True)
def reset_sidereal_mode():
    """Reset sidereal mode between tests."""
    eph.set_sid_mode(SIDM_FAGAN_BRADLEY)
    yield
    eph.set_sid_mode(SIDM_FAGAN_BRADLEY)
```

### Checking SPK Availability in Tests

To verify SPK kernels are properly configured:

```python
def test_spk_available_for_chiron():
    """Verify SPK is available for Chiron."""
    from libephemeris import is_spk_available_for_body
    from libephemeris.constants import CHIRON
    
    assert is_spk_available_for_body(CHIRON), (
        "Chiron SPK not available - ensure setup_spk fixture ran"
    )
```

### Supported Major Asteroids

The following bodies have automatic SPK download support:

| Constant | Body | Asteroid # | Description |
|----------|------|------------|-------------|
| `CERES` | Ceres | 1 | Dwarf planet in asteroid belt |
| `PALLAS` | Pallas | 2 | Second-largest asteroid |
| `JUNO` | Juno | 3 | Main belt asteroid |
| `VESTA` | Vesta | 4 | Second-most-massive asteroid |
| `CHIRON` | Chiron | 2060 | Centaur, important in astrology |

Other bodies (TNOs like Eris, Sedna, etc.) require manual SPK download using `download_and_register_spk()`.

### Troubleshooting

**SPK download fails in CI:**
- Ensure `astroquery` is installed: `pip install astroquery`
- Check network access to `ssd.jpl.nasa.gov`
- Use `LIBEPHEMERIS_LOG_LEVEL=DEBUG` to see download details

**Tests are slow on first run:**
- This is expected - SPK files are ~1-5 MB each
- Use CI caching for the `~/.libephemeris/spk/` directory
- Consider using `pytest-xdist` for parallel test execution

**Different results between runs:**
- Ensure `set_strict_precision(False)` when accepting Keplerian fallback
- Check that SPK coverage includes the test dates (default: 10 years around current date)
