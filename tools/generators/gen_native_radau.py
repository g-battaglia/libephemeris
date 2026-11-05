#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0
"""Derive the native eight-node left-Radau integration tables from mathematics.

For n=8, the nodes on [-1,1] are the zeros of P_8 + P_7, including -1.
Transform them to [0,1], form the Newton interpolation basis, and integrate
that basis once and twice. No astronomical engine or implementation table is
an input. The construction follows Gauss-Radau quadrature and polynomial
integration, as used in Rein & Spiegel (2015), section 2.1.

Decimal arithmetic is confined to this offline generator. Independent builds
at 100 and 140 decimal digits must round every emitted word identically.
"""

from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "core/src/bodies/nbody/native/radau_tables.rs"
NODE_COUNT = 8


def legendre(order: int, argument: Decimal) -> Decimal:
    """Evaluate P_n by its defining three-term recurrence."""
    previous, current = Decimal(1), argument
    if order == 0:
        return previous
    for degree in range(1, order):
        previous, current = current, (
            (2 * degree + 1) * argument * current - degree * previous
        ) / (degree + 1)
    return current


def radau_polynomial(argument: Decimal) -> Decimal:
    """The left endpoint is a prescribed zero of P_8(x) + P_7(x)."""
    return legendre(NODE_COUNT, argument) + legendre(NODE_COUNT - 1, argument)


def nodes(precision: int) -> list[Decimal]:
    """Isolate each interior root, then bisect without approximate root seeds."""
    # The rational mesh only brackets sign changes. All refined roots are
    # required, independently of any published or upstream decimal literals.
    mesh = 512
    left = -Decimal(1) + Decimal(2) / mesh
    left_value = radau_polynomial(left)
    roots = [Decimal(0)]
    for index in range(2, mesh + 1):
        right = -Decimal(1) + Decimal(2 * index) / mesh
        right_value = radau_polynomial(right)
        if left_value * right_value < 0:
            low, high, low_value = left, right, left_value
            for _ in range(4 * precision):
                middle = (low + high) / 2
                if middle == low or middle == high:
                    break
                value = radau_polynomial(middle)
                if value == 0:
                    low = high = middle
                    break
                if low_value * value < 0:
                    high = middle
                else:
                    low, low_value = middle, value
            roots.append(((low + high) / 2 + 1) / 2)
        left, left_value = right, right_value
    if len(roots) != NODE_COUNT:
        raise ValueError(f"expected eight distinct Radau nodes, found {len(roots)}")
    return roots


def multiply_linear(coefficients: list[Decimal], root: Decimal) -> list[Decimal]:
    """Multiply an ascending-power polynomial by (s - root)."""
    result = [Decimal(0)] * (len(coefficients) + 1)
    for degree, value in enumerate(coefficients):
        result[degree] -= root * value
        result[degree + 1] += value
    return result


def integrated_basis(
    coefficients: list[Decimal], argument: Decimal, integrations: int
) -> Decimal:
    """Integrate from zero, with zero integration constants, once or twice."""
    total = Decimal(0)
    for degree, coefficient in enumerate(coefficients):
        denominator = Decimal(degree + 1)
        if integrations == 2:
            denominator *= degree + 2
        total += coefficient * argument ** (degree + integrations) / denominator
    return total


def derive(precision: int) -> tuple[list[Decimal], list[list[Decimal]], list[list[Decimal]]]:
    """Build nodes and integrals, checking all Radau polynomial moments."""
    with localcontext() as context:
        context.prec = precision
        points = nodes(precision)
        bases = [[Decimal(1)]]
        for root in points[:-1]:
            bases.append(multiply_linear(bases[-1], root))
        stages = points + [Decimal(1)]
        velocity = [[integrated_basis(basis, stage, 1) for basis in bases] for stage in stages]
        position = [[integrated_basis(basis, stage, 2) for basis in bases] for stage in stages]
        bound = Decimal(10) ** (-(precision - 20))
        # Use independently constructed Lagrange cardinal polynomials to
        # verify quadrature exactness through degree 2*n-2, not just n-1.
        weights = []
        for selected, node in enumerate(points):
            basis = [Decimal(1)]
            denominator = Decimal(1)
            for other, root in enumerate(points):
                if other != selected:
                    basis = multiply_linear(basis, root)
                    denominator *= node - root
            weights.append(integrated_basis(basis, Decimal(1), 1) / denominator)
        for degree in range(2 * NODE_COUNT - 1):
            moment = sum(weight * (node ** degree if degree else Decimal(1))
                         for weight, node in zip(weights, points, strict=True))
            if abs(moment - Decimal(1) / (degree + 1)) > bound:
                raise ValueError(f"Radau moment {degree} is not exact at {precision} digits")
        return points, velocity, position


def render(precision: int) -> str:
    """Emit correctly rounded binary64 literals with their mathematical origin."""
    points, velocity, position = derive(precision)
    lines = [
        "// SPDX-License-Identifier: AGPL-3.0",
        "//! Generated left-Radau nodes and integrated Newton basis polynomials.",
        "//! Source: zeros of `P_8 + P_7` on `[-1,1]`, mapped to `[0,1]`.",
        "//! No third-party implementation or numerical snapshot is an input.",
        "//! Regenerate: python3 tools/generators/gen_native_radau.py",
        "//! Both 100- and 140-digit derivations must emit identical binary64 words.",
        "//! Rows 0..8 are the eight nodes; row 8 is the requested step endpoint.",
        "//! Column k integrates `product_{j<k}(s - node[j])` once or twice.",
        "", "#[rustfmt::skip]",
        "pub(super) const NODES: [f64; 8] = [",
    ]
    lines.extend(f"    {float(point)!r}," for point in points)
    lines.append("];\n")
    for name, table in (("VELOCITY_INTEGRALS", velocity), ("POSITION_INTEGRALS", position)):
        lines.extend(["#[rustfmt::skip]", f"pub(super) const {name}: [[f64; 8]; 9] = ["])
        for row in table:
            lines.append("    [" + ", ".join(repr(float(value)) for value in row) + "],")
        lines.append("];\n")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="verify without writing")
    arguments = parser.parse_args()
    content = render(100)
    if render(140) != content:
        raise ValueError("independent precision runs disagree after rounding to binary64")
    if arguments.check:
        if not OUTPUT.is_file() or OUTPUT.read_text(encoding="utf-8") != content:
            raise ValueError(f"generated table is stale: {OUTPUT}")
        print("Native Radau tables match both precision runs and all 15 exact moments.")
    else:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(content, encoding="utf-8")
        print(f"Wrote {OUTPUT.relative_to(ROOT)} from both verified precision runs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
