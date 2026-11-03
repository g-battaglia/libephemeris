"""The eclipse module retains its distinct scalar and eight-detail models."""

from data_paths import CORE_LEB, DATA_ROOT
from pathlib import Path
import hashlib
import inspect
import json
import os
import subprocess
import sys

import libephemeris as eph
import pytest

KERNEL = DATA_ROOT / "de440.bsp"
SHA256 = "a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5"


@pytest.fixture(scope="module")
def kernel(tmp_path_factory):
    if not KERNEL.is_file():
        pytest.skip("optional DE440 data unavailable")
    assert hashlib.sha256(KERNEL.read_bytes()).hexdigest() == SHA256
    # This comparison declares DE440 only, not neighboring center supplements.
    # Give both engines the same isolated primary-kernel directory.
    isolated = tmp_path_factory.mktemp("visibility-kernel") / KERNEL.name
    isolated.symlink_to(KERNEL)
    return isolated


@pytest.fixture(autouse=True)
def configured(kernel):
    eph.set_calc_mode("skyfield")
    eph.set_precision_tier("medium")
    eph.set_leb_file(None)
    eph.set_ephe_path(str(kernel.parent))
    eph.set_ephemeris_file(kernel.name)
    eph.set_network_policy("sealed")
    eph.set_auto_spk_download(False)
    eph.set_delta_t_userdef(64 / 86400)
    eph.set_iers_delta_t_enabled(False)
    yield
    eph.set_calc_mode("leb")
    eph.set_precision_tier("base")
    eph.set_leb_file(str(CORE_LEB))


@pytest.mark.parametrize(
    "function,args,bodies",
    [
        (
            "heliacal_ut",
            [2451545.0, [12.5, 41.9, 30.0], [1013.25, 15.0, 40.0, 0.0], [], "Venus", 1],
            (0, 3),
        ),
        (
            "vis_limit_mag",
            [2451545.0, [12.5, 41.9, 30.0], [1013.25, 15.0, 50.0, 0.0], [], "Venus"],
            (0, 3),
        ),
        ("sol_eclipse_magnitude_at_loc", [2460409.264, [-96.8, 32.8, 130.0]], (0, 1)),
        ("sol_eclipse_obscuration_at_loc", [2460409.264, [-96.8, 32.8, 130.0]], (0, 1)),
        ("sol_eclipse_max_time", [2460409.264], (0, 1)),
    ],
)
def test_direct_event_rays_trace_the_consumed_jpl_targets(function, args, bodies):
    token = eph.start_tracing()
    try:
        getattr(eph, function)(*args)
        providers = eph.get_trace_results()
        assert all(providers.get(body) == "Skyfield" for body in bodies)
    finally:
        token.reset()
