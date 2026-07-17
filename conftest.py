"""
Root pytest configuration for all LibEphemeris test suites.

Provides the --calc-mode CLI flag and an autouse fixture that enforces
the chosen calculation mode (auto / skyfield / leb) for every test,
regardless of which test directory is being run.
"""

import pytest
import libephemeris as ephem


@pytest.fixture
def allow_mocked_network():
    """Authorize unit tests whose transport is fully replaced by a fake."""
    ephem.set_network_policy("allow")
    yield
    ephem.set_network_policy(None)


# ============================================================================
# --calc-mode CLI flag
# ============================================================================


def pytest_addoption(parser):
    """Add --calc-mode option to pytest CLI."""
    parser.addoption(
        "--calc-mode",
        action="store",
        default=None,
        choices=["auto", "skyfield", "leb"],
        help="Force libephemeris calculation mode: auto (default), skyfield, or leb",
    )


# ============================================================================
# Session header – always visible at the top of every test run
# ============================================================================


def pytest_report_header(config):
    """Show active calc-mode in the pytest header."""
    mode = config.getoption("--calc-mode", default=None)
    if mode is not None:
        return [f"libephemeris calc-mode: {mode}"]
    # No explicit flag – report what will actually happen
    effective = ephem.get_calc_mode()
    reader = ephem.get_leb_reader()
    if reader is not None:
        return [f"libephemeris calc-mode: auto (LEB active: {reader})"]
    return [f"libephemeris calc-mode: {effective} (no LEB file loaded)"]


# ============================================================================
# Autouse fixture – enforces the mode for every single test
# ============================================================================


@pytest.fixture(autouse=True)
def enforce_calc_mode(request):
    """Enforce --calc-mode for every test, resetting after.

    When --calc-mode is given on the CLI, ``set_calc_mode()`` is called
    before each test and restored to its previous value afterwards.
    Without --calc-mode the fixture is a no-op (preserves current behaviour).
    """
    mode = request.config.getoption("--calc-mode", default=None)
    if mode is None:
        yield
        return

    original = ephem.get_calc_mode()
    ephem.set_calc_mode(mode)
    yield
    ephem.set_calc_mode(original)
