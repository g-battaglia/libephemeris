# SPDX-License-Identifier: AGPL-3.0
"""Prepared independent conic-reference guards; no astronomical backend needed."""

from decimal import Decimal, localcontext
import math
from pathlib import Path
import subprocess
import sys

from tools.generators import gen_native_conic_reference as generator


ROOT = Path(__file__).resolve().parents[2]


def test_independent_precisions_round_to_the_exact_shipped_reference():
    first = generator.render(120)
    assert first == generator.render(160)
    assert generator.OUTPUT.read_text(encoding="utf-8") == first


def test_case_domains_retain_both_sides_of_the_parabolic_boundary():
    cases = generator.parameters()
    assert len(cases) == 47
    assert any(row[1] == 1.0 - 2.0**-52 for row in cases)
    assert any(row[1] == 1.0 + 2.0**-52 for row in cases)
    assert any(row[4] > math.pi / 2 for row in cases)
    assert any(abs(row[2]) > 1e10 for row in cases)
    for axis, eccentricity, _, gm, *_ in cases:
        assert gm > 0
        assert (axis > 0 and 0 <= eccentricity < 1) or (axis < 0 and eccentricity > 1)


def test_circular_and_hyperbolic_pericentre_parameterizations_are_independent():
    assert generator.state([1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0], 120) == [
        1.0, 0.0, 0.0, 0.0, 1.0, 0.0
    ]
    with localcontext() as context:
        context.prec = 120
        speed = float((Decimal(5) / 2).sqrt())
    assert generator.state([-2.0, 1.5, 0.0, 1.0, 0.0, 0.0, 0.0], 120) == [
        1.0, 0.0, 0.0, 0.0, speed, 0.0
    ]


def test_signed_anomaly_symmetry_does_not_need_a_second_engine():
    for axis, eccentricity in [(1.0, 0.99), (-1.0, 1.01)]:
        positive = generator.state([axis, eccentricity, 0.1, 1.0, 0.0, 0.0, 0.0], 120)
        negative = generator.state([axis, eccentricity, -0.1, 1.0, 0.0, 0.0, 0.0], 120)
        assert negative == [positive[0], -positive[1], positive[2],
                            -positive[3], positive[4], positive[5]]


def test_taylor_trigonometry_has_an_exact_angle_witness():
    with localcontext() as context:
        context.prec = 120
        circle = 2 * generator.pi()
        sine, cosine = generator.sin_cos(circle / 12, circle)
        assert abs(sine - Decimal("0.5")) < Decimal("1e-110")
        assert abs(cosine - Decimal(3).sqrt() / 2) < Decimal("1e-110")


def test_check_is_read_only_and_uses_isolated_standard_library():
    before = generator.OUTPUT.read_bytes()
    result = subprocess.run(
        [sys.executable, "-I", "-B", str(ROOT / "tools/generators/gen_native_conic_reference.py"), "--check"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "47 exact-input states" in result.stdout
    assert generator.OUTPUT.read_bytes() == before
