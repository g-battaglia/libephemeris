# Testing Guide

LibEphemeris has backend-specific and topic-specific test commands. Test counts
change as coverage grows, so this guide documents stable command contracts
rather than snapshots of a particular run.

## Setup

Install the project and development dependencies in the active environment:

```bash
uv pip install -e ".[dev]"
```

Use `uv run ./leph ...` if the virtual environment is not activated. See the
[CLI reference](https://github.com/g-battaglia/libephemeris/blob/main/CLI.md)
for shell completion and the complete command
tree.

## Required test policy

**Do not run bare `pytest` or an all-project/full-suite command.** Select the
smallest backend subgroup or explicit test file that covers the change. This
keeps feedback fast and avoids unintentionally running network, slow, or
large-data jobs.

Recommended gates:

| Purpose | Command |
|---|---|
| Fast Skyfield sanity check | `uv run ./leph test skyfield essential` |
| Recommended broad unit gate | `uv run ./leph test leb-backend unit-fast` |
| One file | `pytest tests/test_file.py -v` |
| One test | `pytest tests/test_file.py::test_name -v` |
| Static checks | `uv run ./leph code lint`, `uv run ./leph code format`, `uv run ./leph code typecheck` |

The LEB backend gate requires the configured LEB data described by
`uv run ./leph test leb-backend --help`. If that data is unavailable, use the Skyfield
subgroup or a targeted file and report the limitation.

## Test groups

The developer CLI is the source of truth for available subcommands. Inspect it
with `uv run ./leph test --help` and `uv run ./leph test <group> --help`.

### Skyfield backend

```bash
uv run ./leph test skyfield essential    # quick, one representative test per module
uv run ./leph test skyfield smoke        # broader representative selection
uv run ./leph test skyfield unit-fast    # unit tests, parallel, excluding @slow
uv run ./leph test skyfield unit         # same scope, sequential and verbose
```

### LEB backend

```bash
uv run ./leph test leb-backend essential
uv run ./leph test leb-backend unit-fast  # recommended broad development gate
uv run ./leph test leb-backend unit       # sequential and verbose
```

### Targeted domains

```bash
uv run ./leph test lunar all
uv run ./leph test leb-format all
uv run ./leph test leb2-format all
uv run ./leph test horizons precision-quick  # requires network access
uv run ./leph test coverage run
```

Use the subgroup help before selecting a larger or data-dependent variant.

## Pytest markers

The registered markers are defined in `pytest.ini`:

| Marker | Scope |
|---|---|
| `slow` | Long-running tests |
| `network` | Live network access |
| `integration` | Integration behavior |
| `unit` | Unit behavior |
| `precision` | High-precision validation |
| `edge_case` | Boundary and error behavior |
| `benchmark` | Performance measurements |
| `leb_compare*` | LEB/Skyfield and cross-tier comparisons |

For a targeted direct run, exclusions can be made explicit:

```bash
pytest tests/test_file.py -v -m "not slow and not network"
```

An `xfail` is not a pass: it records a narrowly documented limitation. An
`xpass` should be investigated because the expectation may now be obsolete.
Always read the current marker reason in the test rather than relying on a
central list that can drift.

## Network-gated tests

Most tests are offline. Live SPK tests are opt-in because they contact NASA JPL
Horizons and write cache files:

| Area | Opt-in variable |
|---|---|
| Automatic SPK download | `LIBEPHEMERIS_TEST_SPK_AUTO_DOWNLOAD=1` |
| Direct SPK download | `LIBEPHEMERIS_TEST_SPK_DOWNLOAD=1` |
| Planetary-moon SPK data | `LIBEPHEMERIS_TEST_MOON_SPK=1` |

The runtime auto-SPK client sends direct HTTPS requests to JPL Horizons; it
does not require Astroquery. Keep network tests isolated and cache their data
outside the repository when practical.

Example targeted run:

```bash
LIBEPHEMERIS_TEST_SPK_AUTO_DOWNLOAD=1 \
pytest tests/test_spk_auto.py -v -m network
```

## Backend isolation

Tests that change global configuration must restore it. Prefer
`reset_session()` between independent calculations when ephemeris resources
can be retained; use `close()` only when a test specifically needs full
resource teardown. Context tests should use separate `EphemerisContext`
instances and must not assume that mutating module-global state changes an
existing context.

For source-selection assertions, use the structured tracing API:

```python
import libephemeris as swe

token = swe.start_tracing()
try:
    swe.calc_ut(2451545.0, swe.SUN, swe.FLG_SPEED)
    sources = swe.get_trace_results()
finally:
    token.var.reset(token)

assert swe.SUN in sources
```

See the [tracing guide](../guides/tracing.md) for source labels and nested or
threaded usage.

## Downstream SPK tests

Downstream projects should not make every ordinary unit test depend on a live
download. Either pre-populate a persistent CI cache, or isolate the network
setup in a session fixture and skip cleanly when JPL is unavailable:

```python
import pytest
import libephemeris as swe


@pytest.fixture(scope="session")
def chiron_spk():
    swe.set_auto_spk_download(True)
    try:
        if not swe.ensure_major_asteroid_spk(swe.CHIRON):
            pytest.skip("Chiron SPK is unavailable")
        yield
    finally:
        swe.set_auto_spk_download(False)
```

Tests that intentionally accept the lower-precision analytical fallback may
temporarily call `set_strict_precision(False)`, but must restore the prior
setting afterward. Precision tests should leave strict precision enabled so a
missing kernel fails explicitly instead of changing the numerical method.

## Reporting a verification run

Record the exact command, backend, selected scope, and result. Do not summarize
a targeted subgroup as the full suite, and do not copy test counts into durable
documentation: the command and exit status are the reproducible evidence.
