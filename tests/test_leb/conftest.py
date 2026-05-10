"""
Test fixtures for LEB (binary ephemeris) tests.

Provides:
- Small test .leb files generated on-the-fly
- LEBReader fixture for reader/pipeline tests
"""

from __future__ import annotations

import os
import sys
import pytest

# Ensure project root is on sys.path so 'from scripts.generate_leb import ...' works
_PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)


@pytest.fixture(scope="session")
def test_leb_file(tmp_path_factory):
    """Generate a small .leb file for testing (5-year range, subset of bodies)."""
    from scripts.generate_leb import assemble_leb
    from libephemeris.time_utils import swe_julday

    path = tmp_path_factory.mktemp("leb") / "test.leb"

    jd_start = swe_julday(2023, 1, 1, 0.0)
    jd_end = swe_julday(2028, 1, 1, 0.0)

    # Use a small subset of bodies for fast generation
    bodies = [
        0,  # SE_SUN
        1,  # SE_MOON
        3,  # SE_VENUS  (needed by heliacal no-kernel tests)
        4,  # SE_MARS
        14,  # SE_EARTH
        10,  # SE_MEAN_NODE
    ]

    assemble_leb(
        output=str(path),
        jd_start=jd_start,
        jd_end=jd_end,
        bodies=bodies,
        workers=1,
        verbose=False,
    )

    return str(path)


@pytest.fixture
def leb_reader(test_leb_file):
    """Open a LEBReader for the test .leb file."""
    from libephemeris.leb_reader import LEBReader

    reader = LEBReader(test_leb_file)
    yield reader
    reader.close()


@pytest.fixture(scope="session")
def test_leb_file_minimal(tmp_path_factory):
    """Generate a minimal .leb file with just Sun and Earth (1-year range)."""
    from scripts.generate_leb import assemble_leb
    from libephemeris.time_utils import swe_julday

    path = tmp_path_factory.mktemp("leb_min") / "minimal.leb"

    jd_start = swe_julday(2024, 1, 1, 0.0)
    jd_end = swe_julday(2025, 1, 1, 0.0)

    assemble_leb(
        output=str(path),
        jd_start=jd_start,
        jd_end=jd_end,
        bodies=[0, 14],  # Sun and Earth only
        workers=1,
        verbose=False,
    )

    return str(path)
