# SPDX-License-Identifier: AGPL-3.0
"""Small-amplitude frame harmonics require rate fidelity as well as positions."""

import math

import numpy as np
import pytest
from numpy.polynomial import chebyshev as cheb

from tools.leb_generation import format as wire
from tools.leb_generation.nutation import fit_nutation


def test_adaptive_nutation_preserves_analytic_harmonic_rates(tmp_path, monkeypatch):
    from libephemeris import cache

    start, end = 2451545.0, 2451577.0

    def harmonic(epoch):
        t = epoch - start
        return 1e-8 * math.sin(2 * math.pi * t / 3), 5e-9 * math.cos(2 * math.pi * t / 5)

    monkeypatch.setattr(cache, "get_cached_nutation", harmonic)
    section = fit_nutation(start, end, tmp_path)
    raw = section.path.read_bytes()
    header = wire.read_nutation_header(raw, 0)
    coefficients = np.frombuffer(raw, dtype="<f8", offset=wire.NUTATION_HEADER_SIZE)
    coefficients = coefficients.reshape(header.segment_count, 2, header.degree + 1)
    assert header.interval_days < 16.0
    for epoch in np.linspace(start, end, 129):
        index = min(int((epoch - start) / header.interval_days), header.segment_count - 1)
        node = 2 * (epoch - (start + index * header.interval_days)) / header.interval_days - 1
        values = cheb.chebval(node, coefficients[index].T)
        rate = cheb.chebval(node, cheb.chebder(coefficients[index].T, axis=0)) * (2 / header.interval_days)
        t = epoch - start
        expected = np.array([1e-8 * (2 * math.pi / 3) * math.cos(2 * math.pi * t / 3),
                             -5e-9 * (2 * math.pi / 5) * math.sin(2 * math.pi * t / 5)])
        assert np.max(np.abs(values - harmonic(epoch))) * (180 / math.pi) * 3600 <= .001
        assert np.linalg.norm(rate - expected) <= 1e-6 * np.linalg.norm(expected)


def test_undefined_nutation_samples_refuse_generation(tmp_path, monkeypatch):
    from libephemeris import cache

    monkeypatch.setattr(cache, "get_cached_nutation", lambda epoch: (float("nan"), 0.0))
    with pytest.raises(ValueError, match="invalid native nutation samples"):
        fit_nutation(2451545.0, 2451549.0, tmp_path)


@pytest.mark.parametrize("constant", [0.0, 1e-4])
def test_constant_frame_has_exactly_zero_fitted_rotation_rate(tmp_path, monkeypatch, constant):
    from libephemeris import cache

    monkeypatch.setattr(cache, "get_cached_nutation", lambda epoch: (constant, constant))
    section = fit_nutation(2451545.0, 2451577.0, tmp_path)
    raw = section.path.read_bytes()
    header = wire.read_nutation_header(raw, 0)
    coefficients = np.frombuffer(raw, dtype="<f8", offset=wire.NUTATION_HEADER_SIZE)
    coefficients = coefficients.reshape(header.segment_count, 2, header.degree + 1)
    assert np.all(coefficients[:, :, 0] == constant)
    assert np.all(coefficients[:, :, 1:] == 0.0)
