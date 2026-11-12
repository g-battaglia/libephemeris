# SPDX-License-Identifier: AGPL-3.0
"""Polynomial derivation and regeneration need no astronomical implementation."""

from decimal import Decimal, localcontext
from pathlib import Path
import subprocess
import sys

from tools.generators import gen_native_radau as generator


ROOT = Path(__file__).resolve().parents[2]


def test_distinct_precision_derivations_round_to_the_shipped_table():
    assert generator.render(100) == generator.render(140)
    assert generator.OUTPUT.read_text(encoding="utf-8") == generator.render(100)


def test_newton_basis_integrals_satisfy_independent_monomial_integrals():
    with localcontext() as context:
        context.prec = 80
        nodes, velocity, position = generator.derive(context.prec)
        # Recover Newton coefficients for each monomial by divided differences.
        # The generator's table construction itself never uses this recurrence.
        for degree in range(generator.NODE_COUNT):
            coefficients = [node ** degree if degree else Decimal(1) for node in nodes]
            for width in range(1, len(nodes)):
                for index in range(len(nodes) - 1, width - 1, -1):
                    coefficients[index] = (
                        coefficients[index] - coefficients[index - 1]
                    ) / (nodes[index] - nodes[index - width])
            for stage, argument in enumerate(nodes + [Decimal(1)]):
                first = sum(a * b for a, b in zip(coefficients, velocity[stage], strict=True))
                second = sum(a * b for a, b in zip(coefficients, position[stage], strict=True))
                assert abs(first - argument ** (degree + 1) / (degree + 1)) < Decimal("1e-60")
                assert abs(second - argument ** (degree + 2) / ((degree + 1) * (degree + 2))) < Decimal("1e-60")


def test_legendre_endpoint_identities_and_polynomial_product():
    with localcontext() as context:
        context.prec = 80
        for order in range(9):
            assert generator.legendre(order, Decimal(1)) == 1
            assert generator.legendre(order, Decimal(-1)) == (-1) ** order
        assert generator.multiply_linear([Decimal(1), Decimal(2)], Decimal(3)) == [
            Decimal(-3), Decimal(-5), Decimal(2)
        ]
        assert generator.radau_polynomial(Decimal(-1)) == 0


def test_check_is_read_only_and_runs_with_isolated_standard_library():
    before = generator.OUTPUT.read_bytes()
    result = subprocess.run(
        [sys.executable, "-I", "-B", str(ROOT / "tools/generators/gen_native_radau.py"), "--check"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "15 exact moments" in result.stdout
    assert generator.OUTPUT.read_bytes() == before
