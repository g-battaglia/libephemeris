"""Targeted regression tests for the v3.0.0 round-2 review fixes.

Each test pins the corrected behaviour of one review finding so a future
regression is caught. Grouped by subsystem; see release-notes for the
finding list.
"""

from __future__ import annotations

import os

import pytest

import libephemeris as le
from libephemeris import SUN
from libephemeris.extinction import calc_airmass
from libephemeris.fixed_stars import resolve_star_name
from libephemeris.state import _resolve_data_dir


class TestParityOneLiners:
    def test_pheno_sun_illuminated_fraction_is_one(self) -> None:
        """pheno_ut: the Sun is fully illuminated (attr[1] == 1.0)."""
        attr = le.pheno_ut(2451545.0, SUN, 0)
        assert attr[0] == 0.0  # phase angle
        assert attr[1] == 1.0  # illuminated fraction (was wrongly 0.0)

    def test_context_set_topo_alt_defaults_to_zero(self) -> None:
        """Context.set_topo parity with the module-level set_topo(alt=0.0)."""
        ctx = le.EphemerisContext()
        ctx.set_topo(12.5, 41.9)  # 2-arg form must not raise
        assert ctx.get_topo() is not None

    def test_regulus_flamsteed_alias_is_32_leo(self) -> None:
        """Regulus is 32 Leonis, not 87 Leonis."""
        assert resolve_star_name("32 Leonis") is not None
        assert resolve_star_name("87 Leonis") is None


class TestRobustness:
    def test_airmass_never_below_one_at_zenith(self) -> None:
        """Physical airmass floor is 1.0 at the zenith for every method."""
        assert calc_airmass(90.0, method="kasten_young") >= 1.0
        assert calc_airmass(90.0, method="rozenberg") >= 1.0
        assert calc_airmass(90.0, method="secant") >= 1.0

    def test_data_dir_expands_tilde_from_env(self) -> None:
        """A leading ~ in LIBEPHEMERIS_DATA_DIR is expanded, not literal."""
        prev = os.environ.get("LIBEPHEMERIS_DATA_DIR")
        os.environ["LIBEPHEMERIS_DATA_DIR"] = "~/customdata"
        try:
            resolved = _resolve_data_dir()
        finally:
            if prev is None:
                del os.environ["LIBEPHEMERIS_DATA_DIR"]
            else:
                os.environ["LIBEPHEMERIS_DATA_DIR"] = prev
        assert resolved == os.path.join(os.path.expanduser("~"), "customdata")
        assert "~" not in resolved


class TestStarResolution:
    def test_greek_single_letter_alias_resolves_in_all_resolvers(self) -> None:
        """'α Leo' must resolve to Regulus in resolve_star_name and _resolve_star2
        (the .upper() of a lowercase Greek letter differs from the key)."""
        from libephemeris.fixed_stars import _resolve_star2, resolve_star_name

        assert resolve_star_name("α Leo") == 1000001
        entry, err = _resolve_star2("α Leo")
        assert err is None
        assert entry is not None and entry.id == 1000001


class TestHypothetical:
    def test_undefined_fictitious_id_raises(self) -> None:
        """An id in the fictitious range but without an ephemeris must raise,
        not return a phantom (0,0,0) position."""
        from libephemeris.hypothetical import calc_hypothetical_position

        with pytest.raises(ValueError):
            calc_hypothetical_position(59, 2451545.0)


class TestTime:
    def test_date_conversion_reform_day_morning_not_shifted(self) -> None:
        """An instant in the 00:00-11:59 window of the Gregorian reform day
        (1582-10-15) must not be misclassified as Julian and shifted +10 days."""
        from libephemeris.time_utils import date_conversion

        year, month, day, hour = date_conversion(1582, 10, 15, 6.0, "g")
        assert (year, month, day) == (1582, 10, 15)

    def test_date_conversion_before_reform_is_julian(self) -> None:
        """A date before the reform is still treated as Julian input."""
        from libephemeris.time_utils import date_conversion

        year, month, day, hour = date_conversion(1582, 10, 4, 6.0, "g")
        # 1582-10-04 (Julian) -> 1582-10-14 (Gregorian)
        assert (year, month, day) == (1582, 10, 14)


class TestSplitDeg:
    def test_nakshatra_index_for_unnormalized_longitude(self) -> None:
        """split_deg with the nakshatra flag keeps the raw index for a
        longitude >= 373.33 deg, wrapping only the exact 360 rollover (idx 27),
        matching the reference API."""
        nak = le.SPLIT_DEG_NAKSHATRA
        assert le.split_deg(360.0, nak)[4] == 0  # exact rollover -> 0
        assert le.split_deg(365.0, nak)[4] == 0  # index 27 -> 0
        assert le.split_deg(373.4, nak)[4] == 28  # index 28 kept
        assert le.split_deg(476.58, nak)[4] == 35
        assert le.split_deg(700.0, nak)[4] == 52


class TestReaderResourceSafety:
    def test_zero_byte_leb_file_does_not_leak_fd(self, tmp_path) -> None:
        """A 0-byte stub raises without leaking the open file descriptor."""
        import gc
        import warnings

        from libephemeris.leb2_reader import LEB2Reader
        from libephemeris.leb_reader import LEBReader

        stub = tmp_path / "empty_base_core.leb"
        stub.write_bytes(b"")

        for reader_cls in (LEBReader, LEB2Reader):
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                with pytest.raises((ValueError, OSError)):
                    reader_cls(str(stub))
                gc.collect()
            assert not [w for w in caught if issubclass(w.category, ResourceWarning)]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
