"""
Comprehensive tests for EphemerisContext thread-safe calculations.

Verifies that EphemerisContext provides independent state,
correct calculations, and thread safety.
"""

from __future__ import annotations

import math
import threading

import pytest

from libephemeris.context import EphemerisContext
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    MARS,
    JUPITER,
    FLG_SPEED,
    FLG_EQUATORIAL,
    FLG_SIDEREAL,
    FLG_HELCTR,
    FLG_TOPOCTR,
    SIDM_LAHIRI,
    SIDM_FAGAN_BRADLEY,
    SIDM_RAMAN,
)


class TestContextBasic:
    """Basic EphemerisContext functionality."""

    @pytest.mark.unit
    def test_context_creates_successfully(self):
        """EphemerisContext can be instantiated."""
        ctx = EphemerisContext()
        assert ctx is not None

    @pytest.mark.unit
    def test_context_calc_ut_returns_valid(self):
        """Context calc_ut returns valid 6-element tuple."""
        ctx = EphemerisContext()
        result, retflag = ctx.calc_ut(2451545.0, SUN, FLG_SPEED)
        assert len(result) == 6
        assert 0 <= result[0] < 360
        assert result[2] > 0

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,name",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MERCURY, "Mercury"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
        ],
    )
    def test_context_all_planets(self, body_id: int, name: str):
        """Context calc_ut works for multiple planets."""
        ctx = EphemerisContext()
        result, _ = ctx.calc_ut(2451545.0, body_id, FLG_SPEED)
        assert len(result) == 6
        assert 0 <= result[0] < 360, f"{name}: lon={result[0]}"

    @pytest.mark.unit
    def test_context_houses_returns_valid(self):
        """Context houses() returns valid cusps and ascmc."""
        ctx = EphemerisContext()
        cusps, ascmc = ctx.houses(2451545.0, 41.9, 12.5, ord("P"))
        assert len(cusps) >= 12
        assert len(ascmc) >= 8
        for i, c in enumerate(cusps[:12]):
            assert 0 <= c < 360, f"Cusp {i + 1}: {c}"


class TestContextSetTopo:
    """Test EphemerisContext set_topo for topocentric calculations."""

    @pytest.mark.unit
    def test_context_set_topo(self):
        """Context set_topo accepts valid coordinates."""
        ctx = EphemerisContext()
        ctx.set_topo(12.5, 41.9, 50.0)
        # Should not raise

    @pytest.mark.unit
    def test_context_topo_affects_result(self):
        """Setting topo in context should affect Moon position."""
        ctx = EphemerisContext()
        jd = 2451545.0

        geo_result, _ = ctx.calc_ut(jd, MOON, FLG_SPEED)

        ctx.set_topo(12.5, 41.9, 50.0)
        topo_result, _ = ctx.calc_ut(jd, MOON, FLG_TOPOCTR | FLG_SPEED)

        # Moon positions should differ between geo and topo
        lon_diff = abs(geo_result[0] - topo_result[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        dist_diff = abs(geo_result[2] - topo_result[2])

        total_diff = lon_diff + dist_diff
        assert total_diff > 1e-6, "Topo and geo Moon should differ"


class TestContextSidereal:
    """Test EphemerisContext sidereal mode."""

    @pytest.mark.unit
    def test_context_set_sid_mode(self):
        """Context set_sid_mode accepts valid mode IDs."""
        ctx = EphemerisContext()
        ctx.set_sid_mode(SIDM_LAHIRI)
        # Should not raise

    @pytest.mark.unit
    def test_context_sidereal_differs_from_tropical(self):
        """Sidereal longitude should differ from tropical."""
        ctx = EphemerisContext()
        jd = 2451545.0

        tropical, _ = ctx.calc_ut(jd, SUN, 0)

        ctx.set_sid_mode(SIDM_LAHIRI)
        sidereal, _ = ctx.calc_ut(jd, SUN, FLG_SIDEREAL)

        diff = abs(tropical[0] - sidereal[0])
        if diff > 180:
            diff = 360 - diff
        # Ayanamsha is ~23-24° for Lahiri at J2000
        assert 20 < diff < 30, f"Tropical-sidereal diff: {diff:.2f}°"

    @pytest.mark.unit
    def test_context_get_sid_mode(self):
        """Context get_sid_mode returns the set mode."""
        ctx = EphemerisContext()
        ctx.set_sid_mode(SIDM_FAGAN_BRADLEY)
        mode = ctx.get_sid_mode()
        assert mode == SIDM_FAGAN_BRADLEY


class TestContextFlagCombinations:
    """Test EphemerisContext with various flag combinations."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "flags,desc",
        [
            (0, "default"),
            (FLG_SPEED, "speed"),
            (FLG_EQUATORIAL, "equatorial"),
            (FLG_HELCTR, "heliocentric"),
            (FLG_SPEED | FLG_EQUATORIAL, "speed+equatorial"),
        ],
    )
    def test_context_flag_combos(self, flags: int, desc: str):
        """Context calc_ut works with various flag combinations."""
        ctx = EphemerisContext()
        result, _ = ctx.calc_ut(2451545.0, MARS, flags)
        assert len(result) == 6
        for i, val in enumerate(result):
            assert math.isfinite(val), f"{desc}: result[{i}]={val}"


class TestContextIndependence:
    """Test that multiple contexts maintain independent state."""

    @pytest.mark.unit
    def test_two_contexts_independent_sid_mode(self):
        """Two contexts should have independent sidereal modes."""
        ctx1 = EphemerisContext()
        ctx2 = EphemerisContext()

        ctx1.set_sid_mode(SIDM_LAHIRI)
        ctx2.set_sid_mode(SIDM_FAGAN_BRADLEY)

        assert ctx1.get_sid_mode() == SIDM_LAHIRI
        assert ctx2.get_sid_mode() == SIDM_FAGAN_BRADLEY

    @pytest.mark.unit
    def test_two_contexts_independent_results(self):
        """Two contexts with different sidereal modes give different results."""
        jd = 2451545.0

        ctx1 = EphemerisContext()
        ctx1.set_sid_mode(SIDM_LAHIRI)
        r1, _ = ctx1.calc_ut(jd, SUN, FLG_SIDEREAL)

        ctx2 = EphemerisContext()
        ctx2.set_sid_mode(SIDM_RAMAN)
        r2, _ = ctx2.calc_ut(jd, SUN, FLG_SIDEREAL)

        # Lahiri and Raman ayanamshas differ by ~1-2°
        diff = abs(r1[0] - r2[0])
        if diff > 180:
            diff = 360 - diff
        assert diff > 0.1, f"Different sid modes same result: diff={diff:.4f}°"


class TestContextThreadSafety:
    """Test that EphemerisContext is thread-safe."""

    @pytest.mark.unit
    def test_concurrent_calculations(self):
        """Multiple threads can use separate contexts concurrently."""
        results = {}
        errors = []

        def calc_in_thread(thread_id: int, body_id: int, sid_mode: int):
            try:
                ctx = EphemerisContext()
                ctx.set_sid_mode(sid_mode)
                jd = 2451545.0
                for _ in range(20):
                    result, _ = ctx.calc_ut(jd, body_id, FLG_SIDEREAL | FLG_SPEED)
                    assert 0 <= result[0] < 360
                    jd += 30.0
                results[thread_id] = True
            except Exception as e:
                errors.append((thread_id, str(e)))
                results[thread_id] = False

        threads = [
            threading.Thread(target=calc_in_thread, args=(0, SUN, SIDM_LAHIRI)),
            threading.Thread(
                target=calc_in_thread, args=(1, MOON, SIDM_FAGAN_BRADLEY)
            ),
            threading.Thread(target=calc_in_thread, args=(2, MARS, SIDM_RAMAN)),
            threading.Thread(
                target=calc_in_thread, args=(3, JUPITER, SIDM_LAHIRI)
            ),
        ]

        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=30)

        assert len(errors) == 0, f"Thread errors: {errors}"
        assert all(results.values()), f"Thread results: {results}"

    @pytest.mark.unit
    def test_concurrent_houses(self):
        """Multiple threads can calculate houses concurrently."""
        results = {}
        errors = []

        def calc_houses_in_thread(thread_id: int, lat: float, lon: float):
            try:
                ctx = EphemerisContext()
                jd = 2451545.0
                for _ in range(10):
                    cusps, ascmc = ctx.houses(jd, lat, lon, ord("P"))
                    assert len(cusps) >= 12
                    for c in cusps[:12]:
                        assert 0 <= c < 360
                    jd += 30.0
                results[thread_id] = True
            except Exception as e:
                errors.append((thread_id, str(e)))
                results[thread_id] = False

        threads = [
            threading.Thread(target=calc_houses_in_thread, args=(0, 41.9, 12.5)),
            threading.Thread(target=calc_houses_in_thread, args=(1, 40.7, -74.0)),
            threading.Thread(target=calc_houses_in_thread, args=(2, 35.7, 139.7)),
            threading.Thread(target=calc_houses_in_thread, args=(3, -33.9, 151.2)),
        ]

        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=30)

        assert len(errors) == 0, f"Thread errors: {errors}"
        assert all(results.values())


class TestContextDateRange:
    """Test EphemerisContext across date ranges."""

    @pytest.mark.unit
    @pytest.mark.parametrize("year", [1800, 1900, 2000, 2100, 2200, 2400])
    def test_context_across_centuries(self, year: int):
        """Context works across centuries."""
        ctx = EphemerisContext()
        jd = 2451545.0 + (year - 2000) * 365.25
        result, _ = ctx.calc_ut(jd, SUN, FLG_SPEED)
        assert 0 <= result[0] < 360
