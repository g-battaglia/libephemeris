"""Native module-only chart lots and angular geometry match pinned public APIs."""

from data_paths import ADMITTED_CORE_HASHES, CORE_LEB
from pathlib import Path
import hashlib
import importlib
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


def test_sect_helpers_match_the_native_complete_chart_route():
    from libephemeris import arabic_parts

    positions = dict(Asc=20, Sun=90, Moon=240, Mercury=100, Venus=120, Sun_lat=0.01)
    kwargs = dict(jd=2460482.5, geo_lat=70, geo_lon=12.5)
    sect = arabic_parts.is_day_chart(90, 20, sun_lat=0.01, **kwargs)
    parts = arabic_parts.calc_all_arabic_parts(positions, **kwargs)
    assert parts["Pars_Fortunae"] == arabic_parts.calc_arabic_part_of_fortune(
        20, 90, 240, sect
    )
    assert parts["Pars_Spiritus"] == arabic_parts.calc_arabic_part_of_spirit(
        20, 90, 240, sect
    )
    assert parts["Pars_Amoris"] == arabic_parts.calc_arabic_part_of_love(20, 120, 90)
    assert parts["Pars_Fidei"] == arabic_parts.calc_arabic_part_of_faith(20, 100, 240)
