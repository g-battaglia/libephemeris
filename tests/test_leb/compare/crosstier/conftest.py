"""
Fixtures for cross-tier consistency tests.

Compares results from different LEB tier files at the same dates
to ensure consistency across base, medium, and extended tiers.
"""

from __future__ import annotations

from typing import Any, Callable, Generator

import pytest

import libephemeris as ephem

from tests.test_leb.compare.conftest import (
    CompareHelper,
    TierTolerances,
    generate_test_dates,
    leb_file_path,
    year_to_jd,
)


# =============================================================================
# OVERLAP RANGES
# =============================================================================

# Base-Medium overlap: 1860-2140 (full base range within medium)
_BASE_MEDIUM_START = year_to_jd(1860)
_BASE_MEDIUM_END = year_to_jd(2140)

# Medium-Extended overlap: 1560-2640 (full medium range within extended)
_MEDIUM_EXT_START = year_to_jd(1560)
_MEDIUM_EXT_END = year_to_jd(2640)

# All-tier overlap: 1860-2140 (base range, narrowest)
_ALL_OVERLAP_START = year_to_jd(1860)
_ALL_OVERLAP_END = year_to_jd(2140)

# Cross-tier tolerance: max of any individual tier tolerance
TOLS_CROSS = TierTolerances.for_tier(
    "medium",
    POSITION_ARCSEC=0.5,
    ECLIPTIC_ARCSEC=0.5,
    SPEED_LON_DEG_DAY=0.05,
    DISTANCE_AU=1e-4,
)


# =============================================================================
# HELPER
# =============================================================================


class CrossTierHelper:
    """Execute functions using two different LEB tier files."""

    def __init__(self, leb_path_a: str, leb_path_b: str):
        self.helper_a = CompareHelper(leb_path_a)
        self.helper_b = CompareHelper(leb_path_b)

    def setup(self) -> None:
        """Save current global state."""
        self.helper_a.setup()

    def teardown(self) -> None:
        """Restore saved global state."""
        self.helper_a.teardown()

    def tier_a(self, fn: Callable, *args: Any, **kwargs: Any) -> Any:
        """Call fn using tier A LEB file."""
        return self.helper_a.leb(fn, *args, **kwargs)

    def tier_b(self, fn: Callable, *args: Any, **kwargs: Any) -> Any:
        """Call fn using tier B LEB file."""
        return self.helper_b.leb(fn, *args, **kwargs)


# =============================================================================
# FIXTURES
# =============================================================================


@pytest.fixture(scope="session")
def base_medium_dates() -> list[float]:
    """50 dates in the base-medium overlap (1860-2140)."""
    return generate_test_dates(50, _BASE_MEDIUM_START, _BASE_MEDIUM_END)


@pytest.fixture(scope="session")
def medium_ext_dates() -> list[float]:
    """50 dates in the medium-extended overlap (1560-2640)."""
    return generate_test_dates(50, _MEDIUM_EXT_START, _MEDIUM_EXT_END)


@pytest.fixture(scope="session")
def all_overlap_dates() -> list[float]:
    """50 dates in the all-tier overlap (1860-2140)."""
    return generate_test_dates(50, _ALL_OVERLAP_START, _ALL_OVERLAP_END)


@pytest.fixture
def cross_base_medium(
    leb_file_base: str, leb_file: str
) -> Generator[CrossTierHelper, None, None]:
    """CrossTierHelper comparing base vs medium."""
    helper = CrossTierHelper(leb_file_base, leb_file)
    helper.setup()
    try:
        yield helper
    finally:
        helper.teardown()


@pytest.fixture
def cross_medium_extended(
    leb_file: str, leb_file_extended: str
) -> Generator[CrossTierHelper, None, None]:
    """CrossTierHelper comparing medium vs extended."""
    helper = CrossTierHelper(leb_file, leb_file_extended)
    helper.setup()
    try:
        yield helper
    finally:
        helper.teardown()


@pytest.fixture
def cross_base_extended(
    leb_file_base: str, leb_file_extended: str
) -> Generator[CrossTierHelper, None, None]:
    """CrossTierHelper comparing base vs extended."""
    helper = CrossTierHelper(leb_file_base, leb_file_extended)
    helper.setup()
    try:
        yield helper
    finally:
        helper.teardown()
