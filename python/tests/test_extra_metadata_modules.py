"""Native table records and orbital parameter conversion at old module paths."""

from data_paths import CORE_LEB

import importlib
import inspect
import json
import math
from pathlib import Path
import subprocess

import libephemeris as e
import pytest

FIXTURE = CORE_LEB


def test_julian_range_validation_keeps_stored_bounds_and_invalid_epoch_class():
    from libephemeris.exceptions import validate_jd_range

    e.set_calc_mode("leb")
    e.set_precision_tier("base")
    e.set_leb_file(str(FIXTURE))
    assert validate_jd_range(2451545.0) is None
    assert validate_jd_range(2451545.0, e.SUN, "unit-test") is None
    coverage = e.get_body_coverage(e.SUN)
    for bad in (coverage.jd_start - 1, coverage.jd_end + 1, float("nan"), float("inf")):
        with pytest.raises(e.EphemerisRangeError) as error:
            validate_jd_range(bad, e.SUN)
        assert error.value.body_id == e.SUN
        if math.isfinite(bad):
            assert error.value.start_jd == coverage.jd_start
            assert error.value.end_jd == coverage.jd_end
        else:
            assert error.value.start_jd is None and error.value.end_jd is None
    # A local analytical body has no stored bounds to reject here.
    assert validate_jd_range(1.0, e.CUPIDO) is None


def test_assist_module_helpers_use_native_discovery_and_empty_install(tmp_path):
    from libephemeris import rebound_integration as nbody

    assert isinstance(nbody.check_assist_data_available(), bool)
    assert nbody.reset_assist_data_cache() is None
    e.set_network_policy("sealed")
    try:
        installed = nbody.download_assist_data(
            str(tmp_path / "assist"),
            planets=False,
            asteroids=False,
            show_progress=False,
            quiet=True,
        )
        assert isinstance(installed, Path) and installed.is_dir()
        assert list(installed.iterdir()) == []
    finally:
        e.set_network_policy(None)
