"""Exception contract for nonfinite angles and rounding flags."""

from __future__ import annotations

import pytest

from libephemeris.utils import split_deg


@pytest.mark.parametrize(
    "angle, flags", [(float("inf"), 17), (-float("inf"), 17), (-float("inf"), 1041)]
)
def test_keep_sign_rounding_infinite_domain_error(angle, flags):
    with pytest.raises(ValueError) as caught:
        split_deg(angle, flags)
    assert str(caught.value) == "math domain error"


@pytest.mark.parametrize(
    "angle, flags", [(float("inf"), 1041), (float("inf"), 33), (float("inf"), 16)]
)
def test_other_infinite_paths_keep_overflow_error(angle, flags):
    with pytest.raises(OverflowError) as caught:
        split_deg(angle, flags)
    assert str(caught.value) == "cannot convert float infinity to integer"


def test_keep_sign_nan_keeps_integer_conversion_error():
    with pytest.raises(ValueError) as caught:
        split_deg(float("nan"), 17)
    assert str(caught.value) == "cannot convert float NaN to integer"
