# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for ``libephemeris.moon_theories.uranian``.

This module targets the loop-exhaustion branch in ``_solve_kepler`` where
Newton-Raphson iteration completes all 20 iterations without the
convergence ``break`` firing (arc 189->194). That path is exercised by
requesting a tolerance that can never be satisfied for a non-zero
eccentricity, forcing the loop to fall through to ``return E``.
"""

from __future__ import annotations

import math

from libephemeris.moon_theories import uranian


class TestSolveKeplerLoopExhaustion:
    """Cover the non-converging fall-through of ``_solve_kepler``."""

    def test_loop_exhausts_without_break(self):
        """With an unreachable tolerance the loop runs all 20 iterations.

        ``abs(delta) < tolerance`` is never true for tolerance 0.0 with a
        non-zero eccentricity, so the ``break`` (192->193) is skipped and
        the loop exits to ``return E`` (arc 189->194).
        """
        result = uranian._solve_kepler(1.0, 0.05, tolerance=0.0)

        assert isinstance(result, float)
        # Even without hitting the convergence break, Newton-Raphson has
        # effectively converged numerically; verify Kepler's equation holds.
        residual = result - 0.05 * math.sin(result) - 1.0
        assert abs(residual) < 1e-9

    def test_loop_exhausts_with_negative_tolerance(self):
        """A negative tolerance also guarantees the break never fires."""
        result = uranian._solve_kepler(2.5, 0.1, tolerance=-1.0)

        assert isinstance(result, float)
        residual = result - 0.1 * math.sin(result) - 2.5
        assert abs(residual) < 1e-9
