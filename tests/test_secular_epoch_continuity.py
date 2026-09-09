from __future__ import annotations

import math

import pytest

from libephemeris import minor_bodies as mb


@pytest.fixture
def constant_rates(monkeypatch):
    elements = mb.OrbitalElements(
        name="synthetic",
        epoch=2451545.0,
        a=1.0,
        e=0.25,
        i=0.0,
        omega=40.0,
        Omega=20.0,
        M0=30.0,
        n=1.0,
    )
    monkeypatch.setattr(
        mb, "calc_secular_perturbation_rates", lambda _: (0.5, -0.25, 0.125)
    )
    monkeypatch.setattr(mb, "_calc_forced_elements", lambda *args, **kwargs: (0.0,) * 6)
    return elements


@pytest.mark.parametrize("offset", [-1.0, -0.5, -0.125, 0.0, 0.125, 0.5, 1.0])
def test_enabled_rates_apply_at_every_epoch_offset(constant_rates, offset):
    elements = constant_rates
    actual = mb.apply_secular_perturbations(elements, elements.epoch + offset)
    assert actual[:4] == (
        (elements.omega + 0.5 * offset) % 360.0,
        (elements.Omega - 0.25 * offset) % 360.0,
        (elements.M0 + 1.125 * offset) % 360.0,
        1.125,
    )


@pytest.mark.parametrize("side", [-1.0, 1.0])
def test_one_day_boundary_has_only_the_prescribed_motion(constant_rates, side):
    elements = constant_rates
    boundary = elements.epoch + side
    before = math.nextafter(boundary, -math.inf)
    after = math.nextafter(boundary, math.inf)
    first = mb.apply_secular_perturbations(elements, before)
    second = mb.apply_secular_perturbations(elements, after)
    span = after - before
    assert tuple(second[i] - first[i] for i in range(3)) == (
        0.5 * span,
        -0.25 * span,
        1.125 * span,
    )
    assert first[3] == second[3] == 1.125


@pytest.mark.parametrize("offset", [-1.0, -0.5, 0.0, 0.5, 1.0])
def test_explicitly_disabled_rates_preserve_unperturbed_motion(
    constant_rates, offset, monkeypatch
):
    def unexpected_rates(*args, **kwargs):
        pytest.fail("Disabled perturbations must not evaluate secular rates")

    monkeypatch.setattr(mb, "calc_secular_perturbation_rates", unexpected_rates)
    elements = constant_rates
    actual = mb.apply_secular_perturbations(
        elements, elements.epoch + offset, include_perturbations=False
    )
    assert actual == (
        elements.omega,
        elements.Omega,
        (elements.M0 + elements.n * offset) % 360.0,
        elements.n,
        elements.e,
        elements.i,
    )
