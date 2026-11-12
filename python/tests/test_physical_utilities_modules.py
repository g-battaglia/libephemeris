"""Public scalar modules call the native model and retain their v3 signatures."""

from pathlib import Path
import importlib
import inspect
import json
import math
import subprocess

import pytest


def test_physical_modules_keep_their_result_layout_and_public_error_alias():
    from libephemeris import house_constructions, schaefer, PolarCircleError

    assert house_constructions.PolarCircleError is PolarCircleError
    cusps = house_constructions.houses_savard_a(123, 42, 23.4, 20, 280)
    assert isinstance(cusps, list) and len(cusps) == 13 and cusps[0] == 0
    visibility = schaefer.is_object_visible(1, 30, -15)
    assert isinstance(visibility, tuple) and len(visibility) == 3
    conditions = schaefer.get_visibility_conditions(-15, 30, 1)
    assert isinstance(conditions, dict) and len(conditions) == 13
    with pytest.raises(PolarCircleError):
        house_constructions.houses_sunshine_makransky(123, 80, 23.4, 20, 280, 23)


@pytest.mark.parametrize("index", [-1, 0, 13])
def test_apc_invalid_cusp_indices_raise_key_error(index):
    from libephemeris.house_constructions import apc_cusp

    with pytest.raises(KeyError) as caught:
        apc_cusp(index, 0.5, 0.4, 1.0)
    assert caught.value.args == (index,)
