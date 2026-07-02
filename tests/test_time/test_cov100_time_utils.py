# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for :mod:`libephemeris.time_utils`.

These tests target specific uncovered lines and branch arcs in
``time_utils.py`` that are not exercised by the main suites: calendar
validation guards, the IERS Delta-T fall-through paths in ``deltat`` /
``deltat_ex``, the bytes-calendar decode in ``date_conversion``, the leap
second validation guards in ``utc_to_jd``, the IAU 1982 GMST fallback in
``_gmst06`` (erfa absent), the optional-argument branches in ``sidtime``,
and the sub-minute leap/second normalization in ``utc_time_zone``.

All tests run fully offline; ephemeris access uses the local DE440 kernel
already required by the surrounding suites.
"""

from __future__ import annotations

import sys

import pytest

from libephemeris import state
from libephemeris.constants import GREG_CAL
from libephemeris.exceptions import Error
from libephemeris import iers_data
from libephemeris import time_utils
from libephemeris.time_utils import (
    date_conversion as calendar_convert,
    deltat,
    deltat_ex,
    julday,
    revjul,
    sidtime,
    utc_time_zone,
    utc_to_jd,
)

# Deterministic epoch used throughout (J2000.0).
J2000 = 2451545.0


# ---------------------------------------------------------------------------
# julday / revjul / utc_to_jd / jdet_to_utc / jdut1_to_utc calendar guards
# ---------------------------------------------------------------------------


class TestCalendarGuards:
    """Invalid calendar flags must raise ValueError.

    The reference binding validates the calendar flag in all five
    functions (e.g. `swisseph.julday: invalid calendar (99)`); without
    the guard cal=99 would mean Julian in julday()/revjul() but Gregorian
    in the UTC conversion functions.
    """

    @pytest.mark.unit
    def test_julday_invalid_calendar(self):
        """julday rejects a calendar flag other than GREG_CAL/JUL_CAL."""
        with pytest.raises(ValueError, match="invalid calendar"):
            julday(2000, 1, 1, 12.0, cal=99)

    @pytest.mark.unit
    def test_revjul_invalid_calendar(self):
        """revjul rejects a calendar flag other than GREG_CAL/JUL_CAL."""
        with pytest.raises(ValueError, match="invalid calendar"):
            revjul(J2000, cal=99)

    @pytest.mark.unit
    def test_utc_to_jd_invalid_calendar(self):
        """utc_to_jd rejects a calendar flag other than GREG_CAL/JUL_CAL."""
        with pytest.raises(ValueError, match="invalid calendar"):
            time_utils.utc_to_jd(2000, 1, 1, 12, 0, 0.0, calendar=99)

    @pytest.mark.unit
    def test_jdet_to_utc_invalid_calendar(self):
        """jdet_to_utc rejects a calendar flag other than GREG_CAL/JUL_CAL."""
        with pytest.raises(ValueError, match="invalid calendar"):
            time_utils.jdet_to_utc(J2000, calendar=99)

    @pytest.mark.unit
    def test_jdut1_to_utc_invalid_calendar(self):
        """jdut1_to_utc rejects a calendar flag other than GREG_CAL/JUL_CAL."""
        with pytest.raises(ValueError, match="invalid calendar"):
            time_utils.jdut1_to_utc(J2000, calendar=99)


# ---------------------------------------------------------------------------
# deltat / deltat_ex IERS fall-through paths
# (209->216, 212-214, 294->301, 297-299)
# ---------------------------------------------------------------------------


class TestDeltatIersFallthrough:
    """Exercise the IERS-enabled branches of deltat and deltat_ex."""

    @pytest.mark.unit
    def test_deltat_iers_returns_none_falls_through(self, monkeypatch):
        """deltat: IERS enabled but no data -> fall through to Skyfield."""
        monkeypatch.setattr(state, "_IERS_DELTA_T_ENABLED", True)
        monkeypatch.setattr(state, "_DELTA_T_USERDEF", None)
        monkeypatch.setattr(iers_data, "get_delta_t_iers", lambda jd: None)

        # Falls through to the cached Skyfield value (~63.8 s at J2000).
        result = deltat(J2000)
        assert result == pytest.approx(63.83 / 86400.0, abs=1e-3)

    @pytest.mark.unit
    def test_deltat_iers_raises_falls_through(self, monkeypatch):
        """deltat: IERS raises ValueError -> except/pass -> Skyfield."""
        monkeypatch.setattr(state, "_IERS_DELTA_T_ENABLED", True)
        monkeypatch.setattr(state, "_DELTA_T_USERDEF", None)

        def _boom(jd):
            raise ValueError("simulated IERS failure")

        monkeypatch.setattr(iers_data, "get_delta_t_iers", _boom)

        result = deltat(J2000)
        assert result == pytest.approx(63.83 / 86400.0, abs=1e-3)

    @pytest.mark.unit
    def test_deltat_ex_iers_returns_none_falls_through(self, monkeypatch):
        """deltat_ex: IERS enabled but no data -> fall through to Skyfield."""
        monkeypatch.setattr(state, "_IERS_DELTA_T_ENABLED", True)
        monkeypatch.setattr(state, "_DELTA_T_USERDEF", None)
        monkeypatch.setattr(iers_data, "get_delta_t_iers", lambda jd: None)

        result = deltat_ex(J2000)
        assert result == pytest.approx(63.83 / 86400.0, abs=1e-3)

    @pytest.mark.unit
    def test_deltat_ex_iers_raises_falls_through(self, monkeypatch):
        """deltat_ex: IERS raises ArithmeticError -> except/pass -> Skyfield."""
        monkeypatch.setattr(state, "_IERS_DELTA_T_ENABLED", True)
        monkeypatch.setattr(state, "_DELTA_T_USERDEF", None)

        def _boom(jd):
            raise ArithmeticError("simulated IERS failure")

        monkeypatch.setattr(iers_data, "get_delta_t_iers", _boom)

        result = deltat_ex(J2000)
        assert result == pytest.approx(63.83 / 86400.0, abs=1e-3)


# ---------------------------------------------------------------------------
# date_conversion bytes-calendar decode (line 346)
# ---------------------------------------------------------------------------


class TestDateConversionBytes:
    """The calendar helper accepts a bytes calendar argument."""

    @pytest.mark.unit
    def test_bytes_calendar_decoded(self):
        """A bytes 'j' calendar is decoded and used like the str form."""
        y, m, d, h = calendar_convert(1582, 10, 15, 12.0, b"j")
        assert (y, m, d) == (1582, 10, 5)
        assert h == pytest.approx(12.0, abs=1e-10)


# ---------------------------------------------------------------------------
# utc_to_jd leap-second validation guards (lines 467, 472)
# ---------------------------------------------------------------------------


class TestUtcToJdLeapSecondGuards:
    """second >= 60 is only valid at 23:59:60 on real leap-second dates."""

    @pytest.mark.unit
    def test_second_60_wrong_clock_time_raises(self):
        """second=60 at a non-23:59 clock time raises Error (line 467)."""
        with pytest.raises(Error, match="no leap second"):
            utc_to_jd(2016, 12, 31, 12, 0, 60.0, GREG_CAL)

    @pytest.mark.unit
    def test_second_60_non_leap_date_raises(self):
        """second=60 at 23:59 but not a leap-second date raises (line 472)."""
        with pytest.raises(Error, match="no leap second"):
            utc_to_jd(2017, 6, 30, 23, 59, 60.0, GREG_CAL)

    @pytest.mark.unit
    def test_second_60_valid_leap_date_ok(self):
        """A genuine leap second (2016-12-31 23:59:60) is accepted."""
        jd_et, jd_ut = utc_to_jd(2016, 12, 31, 23, 59, 60.0, GREG_CAL)
        assert jd_et > jd_ut  # TT runs ahead of UT1


# ---------------------------------------------------------------------------
# _gmst06 IAU 1982 fallback (erfa absent) -> lines 946-959
# ---------------------------------------------------------------------------


class TestGmst06Fallback:
    """When erfa cannot be imported, the IAU 1982 polynomial is used."""

    @pytest.mark.unit
    def test_gmst06_falls_back_without_erfa(self, monkeypatch):
        """Force ImportError on 'import erfa' to drive the fallback path."""
        # sys.modules[name] = None makes `import name` raise ImportError.
        monkeypatch.setitem(sys.modules, "erfa", None)

        gmst_rad = time_utils._gmst06(J2000)
        assert isinstance(gmst_rad, float)
        # Sanity: GMST at J2000 noon is ~18.697 h -> ~4.894 rad (mod 2pi).
        import math

        gmst_hours = (math.degrees(gmst_rad) / 15.0) % 24.0
        assert gmst_hours == pytest.approx(18.697, abs=0.05)


# ---------------------------------------------------------------------------
# sidtime optional-argument branch arcs
# (1007->1021, 1016->1018, 1018->1021)
# ---------------------------------------------------------------------------


class TestSidtimeBranches:
    """Cover every combination of provided/omitted obliquity & nutation."""

    @pytest.mark.unit
    def test_both_provided_skips_ephemeris(self):
        """Both obliquity and nutation given -> skip the calc block (1007->1021)."""
        result = sidtime(J2000, 0.0, 23.44, 0.0)
        assert 0.0 <= result < 24.0

    @pytest.mark.unit
    def test_obliquity_given_nutation_none(self):
        """obliquity given, nutation None -> 1016->1018 then set nutation."""
        result = sidtime(J2000, 0.0, obliquity=23.44, nutation=None)
        assert 0.0 <= result < 24.0

    @pytest.mark.unit
    def test_nutation_given_obliquity_none(self):
        """nutation given, obliquity None -> set obliquity then 1018->1021."""
        result = sidtime(J2000, 0.0, obliquity=None, nutation=0.0)
        assert 0.0 <= result < 24.0


# ---------------------------------------------------------------------------
# utc_time_zone sub-minute second normalization
# (1396-1397, 1399-1400)
# ---------------------------------------------------------------------------


class TestUtcTimeZoneSubMinute:
    """Sub-minute timezone offsets push out_second over/under bounds."""

    @pytest.mark.unit
    def test_second_overflow_carries_minute(self):
        """A sub-minute west offset drives out_second >= 61 (1396-1397).

        offset = -0.5/60 h -> offset_min rounds to 0, residue -0.5 min,
        so offset_sec = -30 and out_second = second + 30.
        """
        offset = -0.5 / 60.0
        y, mo, d, hh, mm, ss = utc_time_zone(2024, 1, 15, 12, 0, 35.0, offset)
        # 35 + 30 = 65 -> -60 -> 5.0, minute carries +1.
        assert ss == pytest.approx(5.0, abs=1e-9)
        assert (y, mo, d, hh, mm) == (2024, 1, 15, 12, 1)

    @pytest.mark.unit
    def test_second_underflow_borrows_minute(self):
        """A sub-minute east offset drives out_second < 0 (1399-1400).

        offset = +0.5/60 h -> offset_min rounds to 0, residue +0.5 min,
        so offset_sec = +30 and out_second = second - 30.
        """
        offset = 0.5 / 60.0
        y, mo, d, hh, mm, ss = utc_time_zone(2024, 1, 15, 12, 5, 10.0, offset)
        # 10 - 30 = -20 -> +60 -> 40.0, minute borrows -1.
        assert ss == pytest.approx(40.0, abs=1e-9)
        assert (y, mo, d, hh, mm) == (2024, 1, 15, 12, 4)
