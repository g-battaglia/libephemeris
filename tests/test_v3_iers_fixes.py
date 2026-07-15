# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Regression tests for correctness fixes in ``libephemeris.iers_data``.

These pin five bugs that were fixed:

1. ``get_ut1_utc`` interpolating across the 1 s leap-second step, producing a
   ~0.5 s Delta T dip at midday on leap-second days.
2. ``get_tai_utc`` returning 0.0 for every pre-1972 date instead of applying
   the rate-based UTC (piecewise-linear TAI-UTC) model.
3. ``_parse_delta_t_data`` reading the MJD column as Delta T for the C04 backup
   source (ingesting values of ~50000 s).
4. The finals2000A 2-digit year pivot misdating rows from 2073 onward.
5. Check-then-index races against ``clear_iers_cache`` in the range/availability
   helpers.

All tests inject synthetic data into the module globals, so none of them need
network access. The one test that reads real cached files is skipped when those
files are absent.
"""

from __future__ import annotations

import os

import pytest

import libephemeris.iers_data as iers


@pytest.fixture(autouse=True)
def _isolate_iers_state():
    """Snapshot and restore module globals so tests cannot leak into each other."""
    saved = (
        iers._IERS_DATA,
        iers._LEAP_SECONDS,
        iers._DELTA_T_DATA,
        iers._IERS_AUTO_DOWNLOAD,
    )
    # Never hit the network during these tests.
    iers.set_iers_auto_download(False)
    # Populate the leap-second table so no getter tries to (re)load from disk.
    iers._LEAP_SECONDS = iers._get_hardcoded_leap_seconds()
    try:
        yield
    finally:
        (
            iers._IERS_DATA,
            iers._LEAP_SECONDS,
            iers._DELTA_T_DATA,
            iers._IERS_AUTO_DOWNLOAD,
        ) = saved


# ---------------------------------------------------------------------------
# Fix 1: no Delta T dip at midday on a leap-second day.
# ---------------------------------------------------------------------------


def test_no_delta_t_dip_on_leap_second_day_midday():
    """Midday on a leap-second day must not dip ~0.5 s below the day boundary.

    2016-12-31 (MJD 57753) -> 2017-01-01 (MJD 57754) carried a +1 s leap.
    UT1-UTC steps by ~1 s across that boundary; interpolating UT1-UTC linearly
    used to smear the step across the day and drag Delta T ~0.5 s low at noon.
    """
    # Realistic bracketing UT1-UTC values around the 2017-01-01 leap second.
    iers._IERS_DATA = {
        57753.0: iers.IERSDataPoint(
            mjd=57753.0, year=2016, month=12, day=31, ut1_utc=-0.394
        ),
        57754.0: iers.IERSDataPoint(
            mjd=57754.0, year=2017, month=1, day=1, ut1_utc=0.606
        ),
    }

    boundary = iers.get_ut1_utc(57753.0)
    midday = iers.get_ut1_utc(57753.5)
    assert boundary is not None and midday is not None

    # Within the day the continuous quantity (UT1-TAI) barely drifts, so
    # UT1-UTC at noon stays close to the day-boundary value -- it must NOT jump
    # halfway toward the post-leap value (~+0.1, i.e. ~0.5 s too high).
    assert abs(midday - boundary) < 0.05

    # Express the same thing in Delta T terms: no ~0.5 s dip at midday.
    dt_boundary = iers.get_tai_utc(57753.0) + 32.184 - boundary
    dt_midday = iers.get_tai_utc(57753.5) + 32.184 - midday
    assert abs(dt_midday - dt_boundary) < 0.05
    # And specifically not the ~0.5 s dip the old linear scheme produced.
    assert dt_midday > dt_boundary - 0.1


def test_normal_day_interpolation_unchanged():
    """On a non-leap day the interpolation is still plain linear UT1-UTC."""
    iers._IERS_DATA = {
        58000.0: iers.IERSDataPoint(
            mjd=58000.0, year=2017, month=9, day=4, ut1_utc=0.2000
        ),
        58001.0: iers.IERSDataPoint(
            mjd=58001.0, year=2017, month=9, day=5, ut1_utc=0.2010
        ),
    }
    # No leap second between these days, so midday is the linear midpoint.
    assert iers.get_ut1_utc(58000.5) == pytest.approx(0.20050, abs=1e-9)


# ---------------------------------------------------------------------------
# Fix 2: pre-1972 TAI-UTC is non-zero (rate-based UTC model).
# ---------------------------------------------------------------------------


def test_tai_utc_pre_1972_is_nonzero():
    """Pre-1972 dates must not return 0.0; they follow the rate-based model."""
    # 1961-01-01 (MJD 37300): first tabulated offset ~1.42 s.
    assert iers.get_tai_utc(37300.0) == pytest.approx(1.4228180, abs=1e-6)
    # 1968-01-01 (MJD 39856): ~6.2 s.
    assert iers.get_tai_utc(39856.0) == pytest.approx(6.20533, abs=1e-3)
    # 1970-01-01 (MJD 40587): ~8.0 s.
    assert iers.get_tai_utc(40587.0) == pytest.approx(8.000082, abs=1e-3)


def test_tai_utc_pre_1972_generally_positive_and_bounded():
    """Across 1961-1971 TAI-UTC stays in the ~1.4-10 s band and never 0.

    Note: TAI-UTC is deliberately NOT monotonic here -- the pre-1972 UTC system
    used both rate changes and small offset "steps" (e.g. -0.05 s on
    1961-08-01), so the value can tick down slightly at a step boundary.
    """
    for mjd in range(37300, 41317, 30):
        val = iers.get_tai_utc(float(mjd))
        assert 1.0 <= val <= 10.0
        assert val != 0.0


def test_tai_utc_1972_boundary_is_integer_leap():
    """At 1972-01-01 the integer leap-second table takes over at exactly 10 s."""
    assert iers.get_tai_utc(41317.0) == pytest.approx(10.0, abs=1e-9)


# ---------------------------------------------------------------------------
# Fix 3: C04 backup parse yields a sane Delta T (not the MJD column).
# ---------------------------------------------------------------------------


def test_c04_backup_parse_yields_sane_delta_t(tmp_path):
    """C04 rows (year month day MJD x y UT1-UTC ...) must derive real Delta T."""
    c04 = (
        "# EOP 14 C04 IAU1980 one-file header line\n"
        "2016 1 1 57388 0.05 0.40 -0.1700000 0.0 0.0 0.0\n"
        "2016 1 2 57389 0.05 0.40 -0.1720000 0.0 0.0 0.0\n"
    )
    path = tmp_path / "c04.txt"
    path.write_text(c04)

    pts = iers._parse_delta_t_data(str(path))
    assert len(pts) == 2
    for p in pts:
        # Tens of seconds (~68 s for 2016), definitely not the ~57000 MJD.
        assert 30.0 < p.delta_t < 100.0
    # 36 (TAI-UTC) + 32.184 - (-0.17) = 68.354
    assert pts[0].delta_t == pytest.approx(68.354, abs=1e-3)


def test_primary_deltat_data_still_parses(tmp_path):
    """The primary YEAR MONTH DAY DELTA_T layout is unaffected by the C04 path."""
    primary = "1973  2  1   43.4724\n1973  3  1   43.5648\n"
    path = tmp_path / "deltat.data"
    path.write_text(primary)

    pts = iers._parse_delta_t_data(str(path))
    assert [round(p.delta_t, 4) for p in pts] == [43.4724, 43.5648]


def test_sanity_bound_rejects_mjd_sized_delta_t(tmp_path):
    """An MJD-sized value in the Delta T column is rejected, not ingested."""
    # A 4-column row whose "Delta T" is actually an MJD (~50000) must be dropped.
    data = "1973  2  1   43.4724\n1973  4  1   50000.0\n"
    path = tmp_path / "deltat.data"
    path.write_text(data)

    pts = iers._parse_delta_t_data(str(path))
    assert len(pts) == 1
    assert pts[0].delta_t == pytest.approx(43.4724, abs=1e-6)


# ---------------------------------------------------------------------------
# Fix 4: finals2000A 2-digit year pivot keys off MJD, not the 2-digit year.
# ---------------------------------------------------------------------------


def _make_finals_line(yy2, month, day, mjd, ut1_utc, flag="I"):
    """Build one fixed-width finals2000A.data row (>= 68 chars)."""
    line = list(" " * 80)
    line[0:2] = f"{yy2:02d}"
    line[2:4] = f"{month:02d}"
    line[4:6] = f"{day:02d}"
    line[7:15] = f"{mjd:8.2f}"
    line[57] = flag
    line[58:68] = f"{ut1_utc:10.7f}"
    return "".join(line)


def test_finals_year_pivot_uses_mjd(tmp_path):
    """2073 rows must decode as 2073, not 1973 (old 2-digit pivot bug)."""
    rows = [
        _make_finals_line(73, 1, 2, iers._calendar_to_mjd(1973, 1, 2), 0.8),
        _make_finals_line(0, 1, 1, iers._calendar_to_mjd(2000, 1, 1), 0.35),
        _make_finals_line(99, 12, 31, iers._calendar_to_mjd(1999, 12, 31), 0.7),
        _make_finals_line(73, 6, 15, iers._calendar_to_mjd(2073, 6, 15), 0.1, "P"),
    ]
    path = tmp_path / "finals2000A.data"
    path.write_text("\n".join(rows) + "\n")

    data = iers._parse_finals_data(str(path))
    years = {mjd: dp.year for mjd, dp in data.items()}
    assert years[iers._calendar_to_mjd(1973, 1, 2)] == 1973
    assert years[iers._calendar_to_mjd(2000, 1, 1)] == 2000
    assert years[iers._calendar_to_mjd(1999, 12, 31)] == 1999
    # The regression: MJD >= 51544 -> 2000+, so "73" here is 2073, not 1973.
    assert years[iers._calendar_to_mjd(2073, 6, 15)] == 2073


# ---------------------------------------------------------------------------
# Fix 5: range/availability helpers snapshot the global before indexing.
# ---------------------------------------------------------------------------


def test_range_helpers_do_not_raise_after_clear():
    """Helpers must return cleanly (not IndexError) when the cache is empty."""
    iers._IERS_DATA = {}
    iers._DELTA_T_DATA = []

    assert iers.get_iers_data_range() is None
    assert iers.get_observed_delta_t_data_range() is None
    assert iers.get_observed_delta_t(2451545.0) is None
    assert iers.is_observed_delta_t_available(2451545.0) is False
    assert iers.is_iers_data_available(2451545.0) is False


def test_observed_delta_t_snapshot_is_local():
    """get_observed_delta_t reads a local snapshot, immune to a mid-call clear."""
    iers._DELTA_T_DATA = [
        iers.DeltaTDataPoint(mjd=57388.0, year=2016, month=1, day=1, delta_t=68.10),
        iers.DeltaTDataPoint(mjd=57419.0, year=2016, month=2, day=1, delta_t=68.20),
    ]
    jd = iers._mjd_to_jd(57403.5)
    val = iers.get_observed_delta_t(jd)
    assert val is not None
    assert 68.10 <= val <= 68.20


@pytest.mark.skipif(
    not os.path.exists(
        os.path.join(iers._IERS_CACHE_DIR or "", "finals2000A.data")
        if iers._IERS_CACHE_DIR
        else "/nonexistent"
    ),
    reason="requires a pre-populated IERS finals cache (no network in tests)",
)
def test_real_cache_smoke():  # pragma: no cover - only runs with real data
    """Smoke test against a real cached finals file when one is present."""
    assert iers.load_iers_data() is True
    rng = iers.get_iers_data_range()
    assert rng is not None
