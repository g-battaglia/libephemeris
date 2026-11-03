#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0
"""Generate analytic conic PV references independently at Decimal120/160.

Inputs are exact binary64 words, not astronomical measurements or another
engine's outputs. Bisection solves Kepler's equations independently of the
Rust safeguarded Newton iteration. Trigonometric functions use their defining
Taylor series with a Machin-formula pi; hyperbolic functions use Decimal.exp.
This is an offline numerical test generator, not runtime package mathematics.
"""

from __future__ import annotations

import argparse
from decimal import Decimal, ROUND_FLOOR, localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "core/src/bodies/nbody/native/conic_reference.rs"


def arctangent(value: Decimal) -> Decimal:
    total, term, index = value, value, 1
    while True:
        term *= -value * value
        new = total + term / (2 * index + 1)
        if new == total:
            return total
        total, index = new, index + 1


def pi() -> Decimal:
    return 16 * arctangent(Decimal(1) / 5) - 4 * arctangent(Decimal(1) / 239)


def sin_cos(value: Decimal, circle: Decimal) -> tuple[Decimal, Decimal]:
    value -= ((value + circle / 2) / circle).to_integral_value(rounding=ROUND_FLOOR) * circle
    sine, cosine, sine_term, cosine_term = value, Decimal(1), value, Decimal(1)
    for index in range(1, 1000):
        sine_term *= -value * value / ((2 * index) * (2 * index + 1))
        cosine_term *= -value * value / ((2 * index - 1) * (2 * index))
        next_sine, next_cosine = sine + sine_term, cosine + cosine_term
        if (next_sine, next_cosine) == (sine, cosine):
            return sine, cosine
        sine, cosine = next_sine, next_cosine
    raise ValueError("analytic trigonometric series did not converge")


def bisect(function, high: Decimal, precision: int) -> Decimal:
    low = Decimal(0)
    while function(high) < 0:
        high *= 2
    for _ in range(8 * precision):
        middle = (low + high) / 2
        if middle in (low, high):
            return middle
        value = function(middle)
        if value == 0:
            return middle
        if value < 0:
            low = middle
        else:
            high = middle
    raise ValueError("analytic root did not converge at requested precision")


def parameters() -> list[list[float]]:
    """a, e, M, GM, inclination, node, periapsis argument; radians and common units."""
    rows = []
    for eccentricity in (0.0, 0.2, 0.9, 0.99, 1.0 - 2.0**-52):
        for mean in (0.0, 0.1, -0.1, 1e-30, 1e12):
            rows.append([1.0, eccentricity, mean, 1.0, 0.0, 0.0, 0.0])
    for eccentricity in (1.0 + 2.0**-52, 1.01, 1.5, 3.0):
        for mean in (0.0, 0.1, -0.1, 1e-30, 100.0):
            rows.append([-1.0, eccentricity, mean, 1.0, 0.0, 0.0, 0.0])
    rows.extend([[2.7, 0.12, 0.5, 0.00029591220828559, 0.7, 1.2, -0.4],
                 [-2.0, 1.5, 0.2, 1.0, 2.4, -0.6, 0.9]])
    return rows


def state(words: list[float], precision: int) -> list[float]:
    """Reference parameterization uses direct high-precision physical formulas."""
    with localcontext() as context:
        context.prec = precision
        a, e, mean, gm, inclination, node, periapsis = map(Decimal.from_float, words)
        circle = 2 * pi()
        if a > 0:
            mean -= ((mean + circle / 2) / circle).to_integral_value(rounding=ROUND_FLOOR) * circle
            magnitude = abs(mean)
            anomaly = Decimal(0) if magnitude == 0 else bisect(
                lambda x: x - e * sin_cos(x, circle)[0] - magnitude,
                magnitude + e, precision,
            )
            if mean < 0:
                anomaly = -anomaly
            sine, cosine = sin_cos(anomaly, circle)
            shape = (1 - e*e).sqrt()
            radius = 1 - e*cosine
            x, y = a*(cosine-e), a*shape*sine
            speed = (gm/a).sqrt()
            vx, vy = -speed*sine/radius, speed*shape*cosine/radius
        else:
            magnitude = abs(mean)
            def equation(anomaly):
                exponential = anomaly.exp()
                return e*(exponential - 1/exponential)/2 - anomaly - magnitude
            anomaly = Decimal(0) if magnitude == 0 else bisect(equation, Decimal(1), precision)
            if mean < 0:
                anomaly = -anomaly
            exponential = anomaly.exp()
            sine, cosine = (exponential - 1/exponential)/2, (exponential + 1/exponential)/2
            shape = (e*e-1).sqrt()
            radius = e*cosine-1
            x, y = -a*(e-cosine), -a*shape*sine
            speed = (gm/-a).sqrt()
            vx, vy = -speed*sine/radius, speed*shape*cosine/radius
        si, ci = sin_cos(inclination, circle)
        sn, cn = sin_cos(node, circle)
        sp, cp = sin_cos(periapsis, circle)
        # Three elementary rotations are independent of the Rust P/Q basis.
        def rotate(first, second):
            along, across = cp*first-sp*second, sp*first+cp*second
            return [cn*along-sn*ci*across, sn*along+cn*ci*across, si*across]
        return list(map(float, rotate(x, y) + rotate(vx, vy)))


def render(precision: int) -> str:
    lines = [
        "// SPDX-License-Identifier: AGPL-3.0",
        "//! Generated analytic conic references, not observational ephemerides.",
        "//! Exact binary64 inputs; independent Decimal120/160 bisection agrees.",
        "//! Regenerate: python3 tools/generators/gen_native_conic_reference.py",
        "//! Parameters: [a, e, M, GM, inclination, ascending node, periapsis argument].",
        "//! Output: [x, y, z, vx, vy, vz], common length/time units and radians.",
        "", "// Preserve independently rounded literals, even when a row equals a familiar constant.",
        "#[allow(clippy::approx_constant)]", "#[rustfmt::skip]",
        "pub(super) const CASES: &[([f64; 7], [f64; 6])] = &[",
    ]
    for row in parameters():
        result = state(row, precision)
        lines.append("    ([" + ", ".join(map(repr, row)) + "], [" + ", ".join(map(repr, result)) + "]),")
    lines.append("];\n")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    text = render(120)
    if render(160) != text:
        raise ValueError("conic references disagree after binary64 rounding")
    if args.check:
        if not OUTPUT.is_file() or OUTPUT.read_text(encoding="utf-8") != text:
            raise ValueError("conic reference snapshot is stale")
        print(f"Conic references match Decimal120/160 for {len(parameters())} exact-input states.")
    else:
        OUTPUT.write_text(text, encoding="utf-8")
        print(f"Wrote {len(parameters())} independently checked conic references.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
