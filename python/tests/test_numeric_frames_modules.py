"""Native numeric module paths, original signatures and physical comparisons."""

from data_paths import ADMITTED_CORE_HASHES, CORE_LEB

from concurrent.futures import ThreadPoolExecutor
import importlib
import hashlib
import inspect
import json
import math
from pathlib import Path
import subprocess

import libephemeris as e
import numpy as np
import pytest

MODULES = (
    "astrometry",
    "cache",
    "erfa_nutation",
    "mean_lunar_apse",
    "planetary_mean_elements",
    "precession_vondrak",
    "sidereal_longterm",
    "sidereal_epoch",
    "angles",
)


def test_native_cache_lifecycle_and_numpy_matrix_layouts():
    from libephemeris import cache, erfa_nutation, precession_vondrak, sidereal_epoch

    cache.clear_caches()
    expected = cache.get_cached_nutation(2451545.0)
    with ThreadPoolExecutor(max_workers=8) as pool:
        assert (
            list(pool.map(cache.get_cached_nutation, [2451545.0] * 16))
            == [expected] * 16
        )
    info = cache.get_cache_info()["nutation"]
    assert (
        info["hits"] >= 16
        and info["misses"] >= 1
        and info["currsize"] <= info["maxsize"]
    )
    erfa_nutation.clear_erfa_nutation_cache()
    with pytest.raises(ValueError):
        erfa_nutation.get_erfa_nutation_cached(2451545.0, "unknown")
    assert erfa_nutation.has_erfa()
    # The Vondrak epoch chain contains two first-order IERS frame-bias
    # factors B=I+[b]x. Transpose composition leaves (B B.T)^2, whose exact
    # Frobenius defect is sqrt(2)*(2|b|^2+|b|^4), rather than zero. These are
    # the published bias angles in mas, not fitted matrix-output constants.
    bias_squared = sum(
        math.radians(value / 3_600_000.0) ** 2 for value in (-14.6, -16.6170, -6.8192)
    )
    two_bias_defect = math.sqrt(2.0) * (
        2.0 * bias_squared + bias_squared * bias_squared
    )
    for value, model_defect in (
        (erfa_nutation.get_erfa_pnm06a_matrix(2460310.5), 0.0),
        (sidereal_epoch.equatorial_epoch_matrix(2460310.5), two_bias_defect),
        (
            sidereal_epoch.sidbit_ecliptic_matrix(e.SIDBIT_ECL_T0, 2460310.5, 23.0),
            two_bias_defect,
        ),
    ):
        assert (
            isinstance(value, np.ndarray)
            and value.dtype == np.float64
            and value.shape == (3, 3)
        )
        assert np.linalg.norm(value @ value.T - np.eye(3)) <= model_defect + 2e-14
    assert isinstance(precession_vondrak.vondrak_precession_matrix(2460310.5), tuple)
    assert all(
        isinstance(row, tuple)
        for row in precession_vondrak.vondrak_precession_matrix(2460310.5)
    )
    assert not hasattr(cache, "get_cached_observer_at")
    assert not hasattr(cache, "clear_observer_cache")


def test_generic_lunar_state_retains_physical_radial_velocity_and_domain_errors():
    from libephemeris.interpolated_lunar_apsides import state_to_polar, SUPPORTED_BODIES

    assert SUPPORTED_BODIES == (e.INTP_APOG, e.INTP_PERG)
    state = state_to_polar((2.0, 0.0, 0.0, 0.03, 0.01, 0.0))
    assert state[2] == 2.0 and state[5] == 0.03
    for vector in ((0, 0, 0, 0, 0, 0), (0, 0, 1, 0.01, 0, 0)):
        with pytest.raises(e.CalculationError):
            state_to_polar(vector)
