"""Lunar/hypothetical native models retain selected main resources (#5)."""

from data_paths import CORE_LEB, DATA_ROOT
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import libephemeris as eph
import pytest

KERNEL = DATA_ROOT / "de440.bsp"
MODELS = (10, 11, 12, 13, 21, 22, 40, 48, 58)


@pytest.fixture(scope="module")
def kernel():
    if not KERNEL.is_file():
        pytest.skip("optional DE440 data unavailable")
    assert (
        hashlib.sha256(KERNEL.read_bytes()).hexdigest()
        == "a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5"
    )
    return KERNEL


@pytest.fixture(autouse=True)
def configured(kernel):
    eph.set_calc_mode("skyfield")
    eph.set_precision_tier("medium")
    eph.set_leb_file(None)
    eph.set_ephe_path(str(kernel.parent))
    eph.set_ephemeris_file(kernel.name)
    eph.set_network_policy("sealed")
    eph.set_auto_spk_download(False)
    eph.set_iers_delta_t_enabled(False)
    eph.set_delta_t_userdef(64.0 / 86400.0)
    yield
    eph.set_sid_mode(0)
    eph.set_calc_mode("leb")
    eph.set_precision_tier("base")
    eph.set_leb_file(str(CORE_LEB))


@pytest.mark.parametrize("body", MODELS)
@pytest.mark.parametrize("calculate", ["calc", "calc_ut"])
@pytest.mark.parametrize("flags", [258, 258 | 65536, 258 | 32])
@pytest.mark.parametrize("partial", [False, True])
def test_cold_and_warm_model_calls_share_the_available_selected_main(
    kernel, body, calculate, flags, partial
):
    context = eph.EphemerisContext(str(kernel.parent), kernel.name)
    context.set_sid_mode(eph.SIDM_LAHIRI)
    if partial:
        context.register_spk_body(9999999, str(kernel), 399)
    call = getattr(context, calculate)
    before = call(2451545.0, body, flags)
    assert len(before[0]) == 6 and all(math.isfinite(v) for v in before[0])
    call(2451545.0, eph.SUN, 258)
    call(2451545.0, eph.EARTH, 258)
    assert call(2451545.0, body, flags) == before
    if partial:
        assert context.get_spk_body_info(9999999) == (str(kernel), 399)


@pytest.mark.parametrize(
    "function",
    [
        "calc_true_lunar_node",
        "calc_true_lilith",
        "calc_mean_lunar_node",
        "calc_mean_lilith",
        "calc_mean_lilith_with_latitude",
        "calc_interpolated_apogee",
        "calc_interpolated_perigee",
    ],
)
def test_direct_scalar_lunar_helpers_work_before_and_after_a_planet_call(
    function, kernel
):
    # Each direct global helper gets a genuinely cold process, independent
    # of another test's retained default Runtime image.
    program = r"""
import json,sys
from pathlib import Path
import libephemeris as e
from libephemeris import lunar
kernel=Path(sys.argv[1])
e.set_calc_mode('skyfield');e.set_precision_tier('medium');e.set_leb_file(None)
e.set_ephe_path(str(kernel.parent));e.set_ephemeris_file(kernel.name)
e.set_network_policy('sealed');e.set_auto_spk_download(False);e.set_iers_delta_t_enabled(False)
call=getattr(lunar,sys.argv[2]);before=call(2451545.0)
e.calc_ut(2451545.0,e.SUN,258)
assert call(2451545.0)==before
print(json.dumps(before))
"""
    before = json.loads(
        subprocess.run(
            [sys.executable, "-c", program, str(kernel), function],
            check=True,
            text=True,
            capture_output=True,
            timeout=60,
        ).stdout
    )
    values = before if isinstance(before, list) else (before,)
    assert all(math.isfinite(value) for value in values)


def test_native_model_snapshots_remain_valid_across_concurrent_queries(kernel):
    context = eph.EphemerisContext(str(kernel.parent), kernel.name)
    requests = [(body, 2451545.0 + day) for body in MODELS for day in (0.0, 1.0)]
    expected = [context.calc(epoch, body, 258) for body, epoch in requests]
    with ThreadPoolExecutor(max_workers=8) as workers:
        actual = list(
            workers.map(lambda row: context.calc(row[1], row[0], 258), requests)
        )
    assert actual == expected


def test_missing_selected_main_keeps_sealed_policy_explicit(tmp_path):
    # A normal context may legitimately discover an available home kernel.
    # Make absence real in a new process, without disturbing that user cache.
    environment = {
        key: value
        for key, value in os.environ.items()
        if not key.startswith("LIBEPHEMERIS_")
    }
    environment.update(
        HOME=str(tmp_path),
        LIBEPHEMERIS_MODE="skyfield",
        LIBEPHEMERIS_NETWORK_POLICY="sealed",
        LIBEPHEMERIS_DATA_DIR=str(tmp_path),
        LIBEPHEMERIS_CONFIG=str(tmp_path / "absent.toml"),
        LIBEPHEMERIS_ENV_FILE=str(tmp_path / "absent.env"),
    )
    program = r"""
import sys
import libephemeris as e
e.set_calc_mode('skyfield');e.set_leb_file(None);e.set_network_policy('sealed');e.set_auto_spk_download(False)
context=e.EphemerisContext(sys.argv[1],'de440.bsp')
try:context.calc(2451545.,e.TRUE_NODE,258)
except e.NetworkSealedError:pass
else:raise AssertionError('a missing sealed main image must be refused')
"""
    subprocess.run(
        [sys.executable, "-c", program, str(tmp_path)],
        env=environment,
        check=True,
        text=True,
        capture_output=True,
        timeout=60,
    )


def test_available_selected_main_keeps_its_actual_model_coverage(kernel):
    context = eph.EphemerisContext(str(kernel.parent), kernel.name)
    with pytest.raises(eph.EphemerisRangeError):
        context.calc(1.0, eph.MEAN_NODE, 258)
