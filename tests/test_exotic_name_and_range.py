"""Regression tests for exotic minor-body names and sealed-LEB range errors.

Covers two compatibility fixes:

* ``exotic_bodies.exotic_display_name`` resolves the public IAU minor-planet
  name for the 31 precomputed exotic bodies (feeds ``get_planet_name`` /
  SPK diagnostics), matching the reference, which reads the name from its
  asteroid ephemeris files.
* In sealed ``leb`` mode, an exotic requested outside every available LEB
  interval raises ``EphemerisRangeError`` naming the *requested* body, not the
  internal support Sun ("Body 0").
"""

from __future__ import annotations

import logging

import pytest

import libephemeris as le
from libephemeris.exceptions import EphemerisRangeError
from libephemeris.exotic_bodies import EXOTIC_IDS, exotic_display_name


class TestExoticDisplayName:
    """The public name resolver for AST_OFFSET registry bodies."""

    @pytest.mark.unit
    def test_all_registry_exotics_have_names(self):
        for body_id in EXOTIC_IDS:
            assert exotic_display_name(body_id), f"empty name for {body_id}"

    @pytest.mark.unit
    def test_known_names(self):
        assert exotic_display_name(146199) == "Eris"
        assert exotic_display_name(109942) == "Apophis"
        assert exotic_display_name(100377) == "Sedna"
        assert exotic_display_name(10010) == "Hygiea"

    @pytest.mark.unit
    def test_lilith_ast_suffix_normalized(self):
        """Asteroid 1181 (stored 'Lilith-ast') is reported as 'Lilith'."""
        assert exotic_display_name(11181) == "Lilith"

    @pytest.mark.unit
    def test_non_registry_returns_empty_string(self):
        # Empty string lets get_planet_name fall through to its default,
        # matching the reference (which returns '' for an unresolved asteroid).
        assert exotic_display_name(10012345) == ""
        assert exotic_display_name(0) == ""


class TestSealedLebExoticRangeError:
    """Sealed leb: an out-of-range exotic names the requested body."""

    def teardown_method(self):
        le.set_calc_mode(None)
        logging.disable(logging.NOTSET)

    @pytest.mark.unit
    def test_out_of_all_ranges_names_requested_body(self):
        le.set_calc_mode("leb")
        cov = le.get_body_coverage(146199, 2451545.0)
        if cov is None:
            pytest.skip("Eris not in the active LEB reader (core-only install)")
        logging.disable(logging.WARNING)
        # JD ~ year 8977: outside every stored interval, including the Sun's.
        with pytest.raises(EphemerisRangeError) as exc:
            le.calc_ut(5_000_000.0, 146199, le.FLG_SPEED)
        assert exc.value.body_id == 146199
        assert "146199" in str(exc.value)

    @pytest.mark.unit
    def test_calc_tt_path_also_names_requested_body(self):
        le.set_calc_mode("leb")
        cov = le.get_body_coverage(109942, 2451545.0)
        if cov is None:
            pytest.skip("Apophis not in the active LEB reader (core-only install)")
        logging.disable(logging.WARNING)
        with pytest.raises(EphemerisRangeError) as exc:
            le.calc(5_000_000.0, 109942, le.FLG_SPEED)
        assert exc.value.body_id == 109942

    @pytest.mark.unit
    def test_in_range_exotic_still_computed(self):
        le.set_calc_mode("leb")
        cov = le.get_body_coverage(146199, 2451545.0)
        if cov is None:
            pytest.skip("Eris not in the active LEB reader (core-only install)")
        # J2000 is inside every stored interval: must return a real position.
        pos, _ = le.calc_ut(2451545.0, 146199, le.FLG_SPEED)
        assert 0.0 <= pos[0] < 360.0
