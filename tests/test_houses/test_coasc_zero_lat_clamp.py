"""Analytic sign-aware CoAsc (Munkasey) limits at zero latitude."""

from __future__ import annotations

import pytest

import libephemeris as ephem

ARMC = 100.0
EPS = 23.4392911


class TestCoAscZeroLatClamp:
    @pytest.mark.unit
    @pytest.mark.parametrize("hsys", [ord("P"), ord("H")])
    def test_houses_armc_clamp_is_sign_aware(self, hsys):
        _, ascmc_pos = ephem.houses_armc(ARMC, 1e-12, EPS, hsys)
        _, ascmc_zero = ephem.houses_armc(ARMC, 0.0, EPS, hsys)
        _, ascmc_neg = ephem.houses_armc(ARMC, -1e-12, EPS, hsys)
        separation = abs((ascmc_pos[6] - ascmc_neg[6] + 180.0) % 360.0 - 180.0)
        assert separation > 179.0
        # The exact equator resolves to one of the two antipodal one-sided
        # limits. Measured reference behavior: the horizon system 'H' takes the
        # southern (negative-latitude) branch, every other system the northern
        # (positive-latitude) branch.
        if hsys == ord("H"):
            assert ascmc_zero[6] == ascmc_neg[6]
        else:
            assert ascmc_zero[6] == ascmc_pos[6]

    @pytest.mark.unit
    def test_houses_clamp_is_sign_aware(self):
        jd = 2460000.5
        _, ascmc_pos = ephem.houses(jd, 1e-12, 12.5, ord("P"))
        _, ascmc_neg = ephem.houses(jd, -1e-12, 12.5, ord("P"))
        separation = abs((ascmc_pos[6] - ascmc_neg[6] + 180.0) % 360.0 - 180.0)
        assert separation > 179.0
