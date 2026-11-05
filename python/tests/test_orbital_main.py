"""Orbital reductions acquire their planetary image on a genuinely cold call."""

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest
from data_paths import DATA_ROOT


@pytest.fixture(scope="module")
def kernel():
    path = DATA_ROOT / "de440.bsp"
    if not path.is_file():
        pytest.skip("optional hash-verified DE440 is unavailable")
    with path.open("rb") as stream:
        assert hashlib.file_digest(stream, "sha256").hexdigest() == (
            "a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5"
        )
    return path


@pytest.mark.parametrize("function", [
    "get_orbital_elements", "get_orbital_elements_ut", "orbit_max_min_true_distance",
])
def test_orbital_calls_precede_any_position_or_main_query(function, kernel, tmp_path):
    program = r'''
import json, math, sys
from pathlib import Path
import libephemeris as e
kernel = Path(sys.argv[1])
e.set_calc_mode('skyfield'); e.set_precision_tier('medium'); e.set_leb_file(None)
e.set_ephe_path(str(kernel.parent)); e.set_ephemeris_file(kernel.name)
e.set_network_policy('sealed'); e.set_auto_spk_download(False)
e.set_iers_delta_t_enabled(False)
function = sys.argv[2]
call = getattr(e, function)
for invalid_jd in (float('nan'), float('inf'), -float('inf')):
    try:
        call(invalid_jd, e.MARS, 2)
    except ValueError:
        pass
    else:
        raise AssertionError('nonfinite epoch was accepted')
rows = []
for jd in (2451545.0, 2460310.5):
    for body in (e.MOON, e.MARS, e.EARTH):
        values = call(jd, body, 2)
        assert all(math.isfinite(x) for x in values)
        if function.startswith('get_orbital_elements'):
            assert len(values) == 50 and values[17:] == (0.0,) * 33
            a, eccentricity = values[:2]
            assert a > 0.0 and 0.0 <= eccentricity < 1.0
            assert abs(values[15] - a * (1.0 - eccentricity)) <= 1e-12
            assert abs(values[16] - a * (1.0 + eccentricity)) <= 1e-12
        else:
            assert len(values) == 3 and values[0] >= values[1] >= 0.0
            if body == e.EARTH:
                assert values[1:] == (0.0, 0.0)  # Same orbit and geocentric origin.
            else:
                assert values[1] > 0.0
        e.calc(jd, e.SUN, 2)
        assert call(jd, body, 2) == values
        rows.append(values)
print(json.dumps(rows))
'''
    result = subprocess.run(
        [sys.executable, "-I", "-B", "-c", program, str(kernel), function],
        text=True, capture_output=True, check=True, timeout=90,
        env={**os.environ, "LIBEPHEMERIS_DATA_DIR": str(tmp_path)},
    )
    assert len(json.loads(result.stdout)) == 6
