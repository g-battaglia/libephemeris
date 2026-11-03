"""Public station summaries share the native Session position route."""

from data_paths import ADMITTED_CORE_HASHES, CORE_LEB
from pathlib import Path
import hashlib
import inspect
import json
import math
import subprocess

import libephemeris as eph
import pytest

FIXTURE = CORE_LEB


@pytest.fixture(autouse=True)
def configured():
    assert hashlib.sha256(FIXTURE.read_bytes()).hexdigest() in ADMITTED_CORE_HASHES
    eph.set_calc_mode("leb")
    eph.set_precision_tier("base")
    eph.set_leb_file(str(FIXTURE))
    eph.set_delta_t_userdef(64 / 86400)
    eph.set_iers_delta_t_enabled(False)
    eph.set_sid_mode(0)


@pytest.mark.parametrize("body", [0, 1])
def test_station_summary_rejects_nonstationary_luminaries(body):
    from libephemeris import crossing

    with pytest.raises(ValueError):
        crossing.get_station_info(body, 2451545.0)


def test_station_velocity_rejects_a_zero_differencing_interval():
    from libephemeris import crossing

    with pytest.raises(ZeroDivisionError):
        crossing.calc_velocity_at_station(2, 2451545.0, dt=0.0)
