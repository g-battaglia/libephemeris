"""
Integration tests for LEB→Skyfield fallback in eclipse/occultation functions.

These tests verify the fallback path added by GitHub issue #19. The fallback
activates when the LEB file's coverage is smaller than the search window of
an eclipse/occultation/rise-transit function. For standard LEB files (base,
medium, extended) the fallback never fires; these tests simulate the
out-of-range condition by mocking _topo_ecliptic to raise ValueError.

Coverage:
- Unit-level: helper _call_with_leb_skyfield_fallback behaviour
- Integration: each refactored public function falls back correctly and
  produces a result close to the baseline (within search-precision tolerance).
"""

from __future__ import annotations

from unittest.mock import patch, MagicMock

import pytest

from libephemeris import julday
from libephemeris.eclipse import _call_with_leb_skyfield_fallback


# =============================================================================
# UNIT TESTS for _call_with_leb_skyfield_fallback
# =============================================================================


class TestCallWithLebSkyfieldFallback:
    """Unit tests for the helper that drives the fallback pattern."""

    def test_no_reader_calls_impl_with_none(self):
        """When no LEB is available, impl is called once with reader=None."""
        received = []

        def fake_impl(*args, reader, **kwargs):
            received.append(reader)
            return "ok"

        with patch(
            "libephemeris.eclipse._get_leb_reader_safe", return_value=None
        ):
            result = _call_with_leb_skyfield_fallback(fake_impl, 1, 2, 3)

        assert result == "ok"
        assert received == [None]

    def test_reader_available_no_error(self):
        """When LEB is available and impl succeeds, no retry happens."""
        received = []
        mock_reader = MagicMock(name="LEBReader")

        def fake_impl(*args, reader, **kwargs):
            received.append(reader)
            return "leb_ok"

        with patch(
            "libephemeris.eclipse._get_leb_reader_safe", return_value=mock_reader
        ):
            result = _call_with_leb_skyfield_fallback(fake_impl, 1, 2)

        assert result == "leb_ok"
        assert received == [mock_reader]

    def test_outside_range_triggers_skyfield_retry(self):
        """When impl raises 'outside range', retry with reader=None."""
        received = []
        mock_reader = MagicMock(name="LEBReader")

        def fake_impl(*args, reader, **kwargs):
            received.append(reader)
            if reader is not None:
                raise ValueError(
                    "JD 99999.0 outside range [2451545.0, 2469807.0] for body 10"
                )
            return "skyfield_ok"

        with patch(
            "libephemeris.eclipse._get_leb_reader_safe", return_value=mock_reader
        ):
            result = _call_with_leb_skyfield_fallback(fake_impl, 1)

        assert result == "skyfield_ok"
        assert received == [mock_reader, None]

    def test_outside_nutation_range_also_triggers_retry(self):
        """KeyError or ValueError mentioning nutation range also retries."""
        received = []
        mock_reader = MagicMock(name="LEBReader")

        def fake_impl(*args, reader, **kwargs):
            received.append(reader)
            if reader is not None:
                raise ValueError(
                    "JD 99999.0 outside nutation range [2451545.0, 2469807.0]"
                )
            return "skyfield_ok"

        with patch(
            "libephemeris.eclipse._get_leb_reader_safe", return_value=mock_reader
        ):
            result = _call_with_leb_skyfield_fallback(fake_impl, 1)

        assert result == "skyfield_ok"
        assert received == [mock_reader, None]

    def test_unrelated_error_is_reraised(self):
        """Non-range errors are propagated, not swallowed by fallback."""
        mock_reader = MagicMock(name="LEBReader")

        def fake_impl(*args, reader, **kwargs):
            if reader is not None:
                raise ValueError("Some unrelated error")
            return "should_not_reach"

        with patch(
            "libephemeris.eclipse._get_leb_reader_safe", return_value=mock_reader
        ):
            with pytest.raises(ValueError, match="Some unrelated error"):
                _call_with_leb_skyfield_fallback(fake_impl, 1)

    def test_keyerror_with_outside_message_triggers_retry(self):
        """KeyError mentioning 'outside' also triggers fallback."""
        received = []
        mock_reader = MagicMock(name="LEBReader")

        def fake_impl(*args, reader, **kwargs):
            received.append(reader)
            if reader is not None:
                raise KeyError("Body 10 outside range")
            return "ok"

        with patch(
            "libephemeris.eclipse._get_leb_reader_safe", return_value=mock_reader
        ):
            result = _call_with_leb_skyfield_fallback(fake_impl, 1)

        assert result == "ok"
        assert received == [mock_reader, None]


# =============================================================================
# INTEGRATION TESTS: simulate LEB out-of-range via _topo_ecliptic patching
# =============================================================================


def _raising_topo_ecliptic(reader, jd_tt, jd_ut1, ipl, geopos, iflag=0):
    """Stub for fast_calc._topo_ecliptic that always raises out-of-range."""
    raise ValueError(
        f"JD {jd_tt:.4f} outside range [0.0, 0.0] for body {ipl}"
    )


class TestEclipseFallback:
    """Integration tests: refactored public functions fall back to Skyfield."""

    def test_sol_eclipse_when_loc_falls_back_to_skyfield(self):
        """When LEB topocentric calculation raises out-of-range, Skyfield path runs."""
        from libephemeris import sol_eclipse_when_loc

        jd = julday(2024, 1, 1, 0)
        geopos = (-96.797, 32.7767, 0)  # Dallas TX

        # Baseline (uses LEB or Skyfield, whichever is configured)
        expected = sol_eclipse_when_loc(jd, geopos)

        # Now force LEB topocentric to raise out-of-range:
        # the fallback helper should catch this, retry with reader=None,
        # and produce a result via Skyfield.
        with patch(
            "libephemeris.fast_calc._topo_ecliptic",
            side_effect=_raising_topo_ecliptic,
        ):
            fallback = sol_eclipse_when_loc(jd, geopos)

        # Both produce a result (the fallback did not raise)
        assert fallback[0] == expected[0], "Eclipse types should match"
        # Maximum eclipse time matches within a few minutes
        # (LEB and Skyfield agree to <1 arcsec, so eclipse timing matches well)
        if expected[1][0] > 0 and fallback[1][0] > 0:
            assert abs(fallback[1][0] - expected[1][0]) < 1.0 / 24, (
                f"Max eclipse times differ by >1h: "
                f"{fallback[1][0]} vs {expected[1][0]}"
            )

    def test_rise_trans_falls_back_to_skyfield(self):
        """rise_trans falls back to Skyfield when LEB raises out-of-range."""
        from libephemeris import rise_trans
        from libephemeris.constants import SUN

        jd = julday(2024, 6, 21, 0)
        geopos = (12.5, 41.9, 0)  # Rome

        expected = rise_trans(jd, SUN, 1, geopos)  # 1 = rise

        with patch(
            "libephemeris.fast_calc._topo_ecliptic",
            side_effect=_raising_topo_ecliptic,
        ):
            fallback = rise_trans(jd, SUN, 1, geopos)

        # Both return (retflag, tret)
        assert fallback[0] == expected[0], "Return flags should match"
        # Rise time within 5 minutes
        if expected[1][0] > 0 and fallback[1][0] > 0:
            assert abs(fallback[1][0] - expected[1][0]) < 5.0 / (24 * 60), (
                f"Sun rise times differ by >5min: "
                f"{fallback[1][0]} vs {expected[1][0]}"
            )

    def test_lun_eclipse_when_loc_falls_back_to_skyfield(self):
        """lun_eclipse_when_loc falls back to Skyfield when LEB raises."""
        from libephemeris import lun_eclipse_when_loc

        jd = julday(2024, 1, 1, 0)
        geopos = (12.5, 41.9, 0)

        expected = lun_eclipse_when_loc(jd, geopos)

        with patch(
            "libephemeris.fast_calc._topo_ecliptic",
            side_effect=_raising_topo_ecliptic,
        ):
            fallback = lun_eclipse_when_loc(jd, geopos)

        assert fallback[0] == expected[0]
        # Max eclipse time agrees within hours
        if expected[1][0] > 0 and fallback[1][0] > 0:
            assert abs(fallback[1][0] - expected[1][0]) < 1.0
