# SPDX-License-Identifier: AGPL-3.0
"""Fit native nutation while retaining angular-rate precision of date frames."""

import math

import numpy as np
from numpy.polynomial import chebyshev as cheb

from . import format as wire
from .generate import fit_grid
from .writer import SectionData


def fit_nutation(start, end, directory):
    """Refine a uniform TT grid until both value and derivative converge.

    Position accuracy alone does not control the rotation rate of a fitted
    frame. A degree-32 local interpolant supplies a denser independent sample
    grid for checking the published degree-16 derivative. The relative 1e-6
    vector-rate criterion matches the Cartesian source fitter; positions use
    PLAN's 0.001 arcsecond criterion against actual native off-grid samples.
    The higher-degree interpolant is a convergence estimate, not a rigorous
    source error bound. Complete public comparisons remain necessary.

    A failed segment halves the preferred interval for the entire uniform
    section. Eight refinements are allowed; undefined samples or failure to
    converge refuse generation. No source model or public tolerance changes.
    """
    from libephemeris.cache import get_cached_nutation

    degree = 16
    path = directory / "nutation.f64"

    def interpolate(left, interval, order):
        nodes = np.cos(np.pi * (np.arange(order + 1) + 0.5) / (order + 1))
        epochs = left + (nodes + 1.0) * interval / 2.0
        values = np.asarray([get_cached_nutation(float(epoch)) for epoch in epochs], dtype=np.float64)
        if values.shape != (order + 1, 2) or not np.isfinite(values).all():
            raise ValueError("invalid native nutation samples")
        actual_nodes = 2.0 * (epochs - left) / interval - 1.0
        coefficients = np.linalg.solve(cheb.chebvander(actual_nodes, order), values)
        if not np.isfinite(coefficients).all():
            raise ValueError("nonfinite native nutation fit")
        return coefficients

    for refinement in range(9):
        count, interval = fit_grid(start, end, 16.0 / 2**refinement)
        prefix = bytearray(wire.NUTATION_HEADER_SIZE)
        wire.write_nutation_header(prefix, 0, wire.NutationHeader(start, end, interval, degree, 2, count, 0))
        converged = True
        with path.open("wb") as destination:
            destination.write(prefix)
            for index in range(count):
                left = start + index * interval
                coefficients = interpolate(left, interval, degree)
                reference = interpolate(left, interval, 2 * degree)
                epochs = left + (np.linspace(-1.0, 1.0, 9) + 1.0) * interval / 2.0
                if index + 1 == count:
                    epochs[-1] = end
                nodes = 2.0 * (epochs - left) / interval - 1.0
                expected = np.asarray([get_cached_nutation(float(epoch)) for epoch in epochs], dtype=np.float64)
                if expected.shape != (9, 2) or not np.isfinite(expected).all():
                    raise ValueError("invalid off-grid native nutation samples")
                actual = cheb.chebval(nodes, coefficients).T
                rate = cheb.chebval(nodes, cheb.chebder(coefficients, axis=0)).T * (2.0 / interval)
                reference_rate = cheb.chebval(nodes, cheb.chebder(reference, axis=0)).T * (2.0 / interval)
                if not all(np.isfinite(values).all() for values in (actual, rate, reference_rate)):
                    raise ValueError("nonfinite fitted nutation position/rate")
                position_error = np.abs(actual - expected) * (180.0 / math.pi) * 3600.0
                difference = rate - reference_rate
                rate_error = np.hypot(difference[:, 0], difference[:, 1])
                reference_norm = np.hypot(reference_rate[:, 0], reference_rate[:, 1])
                if not (np.isfinite(rate_error).all() and np.isfinite(reference_norm).all()):
                    raise ValueError("nonfinite fitted nutation rate comparison")
                if (np.any(position_error > 0.001) or
                        np.any(rate_error > 1e-6 * reference_norm)):
                    converged = False
                    break
                destination.write(np.asarray(coefficients.T, dtype="<f8").tobytes())
        if converged:
            return SectionData(path, 0, path.stat().st_size)
    raise ValueError("native nutation fit did not converge in position and rate after eight refinements")
