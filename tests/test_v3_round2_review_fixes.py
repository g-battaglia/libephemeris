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


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
