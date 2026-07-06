"""Tests for coercion of the ``backwards`` parameter in crossing functions.

Under plain truthiness a string like ``"forward"`` would silently select a
BACKWARD search. All public crossing functions now interpret direction
strings explicitly via the shared coercion used by the eclipse search
functions (``bool`` passthrough, ``int`` truthiness, ``"forward"``/
``"backward"`` and friends by meaning, anything else raising TypeError).
"""

from __future__ import annotations

import pytest
import libephemeris as swe
from libephemeris.constants import EARTH, FLG_SWIEPH

JD_J2000 = 2451545.0


@pytest.mark.unit
class TestForwardStringDoesNotInvert:
    """A truthy "forward" string must NOT flip the search direction."""

    def test_solcross_ut_forward_string(self):
        jd_bool = swe.solcross_ut(55.0, JD_J2000, FLG_SWIEPH, backwards=False)
        jd_str = swe.solcross_ut(55.0, JD_J2000, FLG_SWIEPH, backwards="forward")
        assert jd_str == jd_bool
        assert jd_str > JD_J2000

    def test_solcross_forward_string(self):
        jd_bool = swe.solcross(55.0, JD_J2000, FLG_SWIEPH, backwards=False)
        jd_str = swe.solcross(55.0, JD_J2000, FLG_SWIEPH, backwards="forward")
        assert jd_str == jd_bool
        assert jd_str > JD_J2000

    def test_mooncross_ut_forward_string(self):
        jd_bool = swe.mooncross_ut(123.0, JD_J2000, FLG_SWIEPH, backwards=False)
        jd_str = swe.mooncross_ut(123.0, JD_J2000, FLG_SWIEPH, backwards="forward")
        assert jd_str == jd_bool
        assert jd_str > JD_J2000

    def test_mooncross_forward_string(self):
        jd_bool = swe.mooncross(123.0, JD_J2000, FLG_SWIEPH, backwards=False)
        jd_str = swe.mooncross(123.0, JD_J2000, FLG_SWIEPH, backwards="forward")
        assert jd_str == jd_bool
        assert jd_str > JD_J2000

    def test_mooncross_node_ut_forward_string(self):
        jd_bool, _, _ = swe.mooncross_node_ut(JD_J2000, FLG_SWIEPH, backwards=False)
        jd_str, _, _ = swe.mooncross_node_ut(JD_J2000, FLG_SWIEPH, backwards="forward")
        assert jd_str == jd_bool
        assert jd_str > JD_J2000

    def test_mooncross_node_forward_string(self):
        jd_bool, _, _ = swe.mooncross_node(JD_J2000, FLG_SWIEPH, backwards=False)
        jd_str, _, _ = swe.mooncross_node(JD_J2000, FLG_SWIEPH, backwards="forward")
        assert jd_str == jd_bool
        assert jd_str > JD_J2000

    def test_helio_cross_ut_forward_string(self):
        jd_bool = swe.helio_cross_ut(EARTH, 200.0, JD_J2000, FLG_SWIEPH, False)
        jd_str = swe.helio_cross_ut(
            EARTH, 200.0, JD_J2000, FLG_SWIEPH, backwards="forward"
        )
        assert jd_str == jd_bool
        assert jd_str > JD_J2000

    def test_helio_cross_forward_string(self):
        jd_bool = swe.helio_cross(EARTH, 200.0, JD_J2000, FLG_SWIEPH, False)
        jd_str = swe.helio_cross(
            EARTH, 200.0, JD_J2000, FLG_SWIEPH, backwards="forward"
        )
        assert jd_str == jd_bool
        assert jd_str > JD_J2000


@pytest.mark.unit
class TestBackwardStringsAndInts:
    """Backward direction strings and ints map like backwards=True."""

    @pytest.mark.parametrize("value", ["backward", "backwards", "true", "1", 1, True])
    def test_solcross_ut_backward_values(self, value):
        jd_bool = swe.solcross_ut(55.0, JD_J2000, FLG_SWIEPH, backwards=True)
        jd_val = swe.solcross_ut(55.0, JD_J2000, FLG_SWIEPH, backwards=value)
        assert jd_val == jd_bool
        assert jd_val < JD_J2000

    @pytest.mark.parametrize("value", ["forwards", "false", "0", "", 0, False])
    def test_solcross_ut_forward_values(self, value):
        jd_bool = swe.solcross_ut(55.0, JD_J2000, FLG_SWIEPH, backwards=False)
        jd_val = swe.solcross_ut(55.0, JD_J2000, FLG_SWIEPH, backwards=value)
        assert jd_val == jd_bool
        assert jd_val > JD_J2000

    def test_mooncross_ut_backward_string(self):
        jd_bool = swe.mooncross_ut(123.0, JD_J2000, FLG_SWIEPH, backwards=True)
        jd_str = swe.mooncross_ut(123.0, JD_J2000, FLG_SWIEPH, backwards="backward")
        assert jd_str == jd_bool
        assert jd_str < JD_J2000


@pytest.mark.unit
class TestInvalidBackwardsValues:
    """Unrecognized strings and unsupported types raise TypeError."""

    @pytest.mark.parametrize(
        "func_args",
        [
            lambda v: swe.solcross_ut(55.0, JD_J2000, FLG_SWIEPH, v),
            lambda v: swe.solcross(55.0, JD_J2000, FLG_SWIEPH, v),
            lambda v: swe.mooncross_ut(123.0, JD_J2000, FLG_SWIEPH, v),
            lambda v: swe.mooncross(123.0, JD_J2000, FLG_SWIEPH, v),
            lambda v: swe.mooncross_node_ut(JD_J2000, FLG_SWIEPH, v),
            lambda v: swe.mooncross_node(JD_J2000, FLG_SWIEPH, v),
            lambda v: swe.helio_cross_ut(EARTH, 200.0, JD_J2000, FLG_SWIEPH, v),
            lambda v: swe.helio_cross(EARTH, 200.0, JD_J2000, FLG_SWIEPH, v),
        ],
        ids=[
            "solcross_ut",
            "solcross",
            "mooncross_ut",
            "mooncross",
            "mooncross_node_ut",
            "mooncross_node",
            "helio_cross_ut",
            "helio_cross",
        ],
    )
    def test_invalid_string_raises_typeerror(self, func_args):
        with pytest.raises(TypeError):
            func_args("sideways")

    def test_unsupported_type_raises_typeerror(self):
        with pytest.raises(TypeError):
            swe.solcross_ut(55.0, JD_J2000, FLG_SWIEPH, backwards=1.5)
