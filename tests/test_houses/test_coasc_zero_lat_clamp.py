"""Regression: sign-aware CoAsc (Munkasey) clamp at |lat| < 1e-10.

The CoAsc formula degenerates as lat -> 0, and the clamp must preserve the
one-sided limits rather than collapsing both sides of zero onto one value:
non-H systems give 180 for lat -> 0+ and 0 for lat -> 0- (reversed for the
Horizontal system 'H'); exact lat == 0 lands on the positive-side limit.
ascmc[6] is coasc2 (M. Munkasey).
"""

from __future__ import annotations

import pytest

import libephemeris as ephem

ARMC = 100.0
EPS = 23.4392911


class TestCoAscZeroLatClamp:
    @pytest.mark.unit
    @pytest.mark.parametrize(
        "hsys,pos_expected,neg_expected",
        [
            (ord("P"), 180.0, 0.0),
            (ord("H"), 0.0, 180.0),
        ],
    )
    def test_houses_armc_clamp_is_sign_aware(
        self, hsys, pos_expected, neg_expected
    ):
        _, ascmc_pos = ephem.houses_armc(ARMC, 1e-12, EPS, hsys)
        _, ascmc_zero = ephem.houses_armc(ARMC, 0.0, EPS, hsys)
        _, ascmc_neg = ephem.houses_armc(ARMC, -1e-12, EPS, hsys)
        assert ascmc_pos[6] == pos_expected
        assert ascmc_zero[6] == pos_expected
        assert ascmc_neg[6] == neg_expected

    @pytest.mark.unit
    def test_houses_clamp_is_sign_aware(self):
        jd = 2460000.5
        _, ascmc_pos = ephem.houses(jd, 1e-12, 12.5, ord("P"))
        _, ascmc_neg = ephem.houses(jd, -1e-12, 12.5, ord("P"))
        assert ascmc_pos[6] == 180.0
        assert ascmc_neg[6] == 0.0
