"""Cold native photometry owns its configured planetary image (#6 regression)."""

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


def test_sentinels_precede_missing_main_and_topocentric_observer(tmp_path):
    program = r"""
import json,math,sys
import libephemeris as e
e.set_calc_mode('skyfield');e.set_precision_tier('medium');e.set_leb_file(None)
e.set_ephe_path(sys.argv[1]);e.set_ephemeris_file('absent-owned-pheno.bsp');e.set_network_policy('sealed')
for function in (e.pheno,e.pheno_ut):
 earth=function(2451545.0,e.EARTH,e.FLG_TOPOCTR)
 assert earth[3]==180.0 and all(value==0.0 for slot,value in enumerate(earth) if slot!=3)
 nutation=function(2451545.0,e.ECL_NUT,e.FLG_TOPOCTR)
 assert all(math.isnan(value) for value in nutation[:3]) and all(value==0.0 for value in nutation[3:])
print('true')
"""
    result = subprocess.run(
        [sys.executable, "-c", program, str(tmp_path)],
        text=True,
        capture_output=True,
        check=True,
        timeout=30,
        env={**os.environ, "LIBEPHEMERIS_DATA_DIR": str(tmp_path)},
    )
    assert json.loads(result.stdout) is True
