"""Tests that verify the DE kernel is NOT loaded when using LEB mode.

Phase 0 of the Skyfield→LEB porting plan. Tests are introduced alongside
each porting phase and start as xfail until the porting is complete.
"""

from __future__ import annotations

from unittest.mock import patch

import pytest


def _block_get_planets(*args, **kwargs):
    raise RuntimeError(
        "get_planets() was called in LEB mode — DE kernel bypass not working"
    )


@pytest.fixture
def leb_mode(test_leb_file):
    """Activate LEB mode and block get_planets() across all module aliases."""
    import libephemeris
    from libephemeris import state

    state._PLANETS = None
    libephemeris.set_leb_file(test_leb_file)

    targets = [
        "libephemeris.state.get_planets",
        "libephemeris.fixed_stars.get_planets",
        "libephemeris.planets.get_planets",
        "libephemeris.heliacal.get_planets",
        "libephemeris.eclipse.get_planets",
    ]

    patchers = []
    patched = []
    for target in targets:
        try:
            p = patch(target, side_effect=_block_get_planets)
            p.start()
            patchers.append(p)
            patched.append(target)
        except AttributeError:
            pass  # Module doesn't import get_planets directly
    assert "libephemeris.state.get_planets" in patched

    yield

    for p in patchers:
        p.stop()

    state._PLANETS = None
    libephemeris.set_leb_file("")


class TestFixedStarsNoKernel:
    """Phase 1: fixstar_ut and batch_fixstars_ut must not load DE kernel."""

    def test_swe_fixstar_ut_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.fixed_stars import fixstar_ut
        from libephemeris.time_utils import julday

        jd = julday(2024, 6, 15, 12.0)
        pos, _name, _flags = fixstar_ut("Regulus", jd, 0)
        assert state._PLANETS is None
        assert 140 < pos[0] < 160

    def test_swe_fixstar_ut_with_speed(self, leb_mode):
        from libephemeris import state
        from libephemeris.constants import FLG_SPEED
        from libephemeris.fixed_stars import fixstar_ut
        from libephemeris.time_utils import julday

        jd = julday(2024, 6, 15, 12.0)
        pos, _name, _flags = fixstar_ut("Regulus", jd, FLG_SPEED)
        assert state._PLANETS is None
        assert pos[3] != 0.0

    def test_swe_batch_fixstars_ut_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.fixed_stars import batch_fixstars_ut
        from libephemeris.time_utils import julday

        jd = julday(2024, 6, 15, 12.0)
        results = batch_fixstars_ut(["Regulus", "Spica", "Aldebaran"], jd, 0)
        assert state._PLANETS is None
        assert all(r is not None for r in results)

    def test_swe_fixstar_ut_j2000(self, leb_mode):
        from libephemeris import state
        from libephemeris.constants import FLG_J2000
        from libephemeris.fixed_stars import fixstar_ut
        from libephemeris.time_utils import julday

        jd = julday(2024, 6, 15, 12.0)
        _pos, _, _ = fixstar_ut("Regulus", jd, FLG_J2000)
        assert state._PLANETS is None

    def test_swe_fixstar_ut_noaberr(self, leb_mode):
        from libephemeris import state
        from libephemeris.constants import FLG_NOABERR
        from libephemeris.fixed_stars import fixstar_ut
        from libephemeris.time_utils import julday

        jd = julday(2024, 6, 15, 12.0)
        _pos, _, _ = fixstar_ut("Regulus", jd, FLG_NOABERR)
        assert state._PLANETS is None


class TestPhenoNoKernel:
    """Phase 3: pheno_ut must not load DE kernel."""

    def test_swe_pheno_ut_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.constants import MARS
        from libephemeris.planets import pheno_ut
        from libephemeris.time_utils import julday

        jd = julday(2024, 6, 15, 12.0)
        result = pheno_ut(jd, MARS, 0)
        assert state._PLANETS is None
        assert result[2] > 0  # elongation > 0

    def test_swe_pheno_ut_moon(self, leb_mode):
        from libephemeris import state
        from libephemeris.constants import MOON
        from libephemeris.planets import pheno_ut
        from libephemeris.time_utils import julday

        jd = julday(2024, 6, 15, 12.0)
        result = pheno_ut(jd, MOON, 0)
        assert state._PLANETS is None
        assert 0 < result[0] < 180  # phase angle

    def test_swe_pheno_ut_sun(self, leb_mode):
        from libephemeris import state
        from libephemeris.constants import SUN
        from libephemeris.planets import pheno_ut
        from libephemeris.time_utils import julday

        jd = julday(2024, 6, 15, 12.0)
        result = pheno_ut(jd, SUN, 0)
        assert state._PLANETS is None
        assert result[0] == 0.0  # phase angle for Sun = 0


class TestHeliacalNoKernel:
    """Phase 3: heliacal functions must not load DE kernel."""

    def test_heliacal_pheno_ut_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.constants import VENUS
        from libephemeris.heliacal import _heliacal_pheno_ut_pythonic
        from libephemeris.time_utils import julday

        jd = julday(2024, 6, 15, 12.0)
        _result, _retflag = _heliacal_pheno_ut_pythonic(
            jd, 41.9, 12.5, body=VENUS
        )
        assert state._PLANETS is None


class TestEclipseNoKernel:
    """Phase 4: eclipse core functions must not load DE kernel."""

    def test_calc_gamma_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.eclipse import _calc_gamma
        from libephemeris.time_utils import julday

        jd = julday(2024, 4, 8, 18.0)
        gamma = _calc_gamma(jd)
        assert state._PLANETS is None
        assert isinstance(gamma, float)

    def test_calc_penumbra_limit_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.eclipse import _calc_penumbra_limit
        from libephemeris.time_utils import julday

        jd = julday(2024, 4, 8, 18.0)
        l1 = _calc_penumbra_limit(jd)
        assert state._PLANETS is None
        assert l1 > 0

    def test_besselian_x_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.eclipse import calc_besselian_x
        from libephemeris.time_utils import julday

        jd = julday(2024, 4, 8, 18.0)
        x = calc_besselian_x(jd)
        assert state._PLANETS is None
        assert isinstance(x, float)

    def test_xyz_pipeline_a_no_double_convert(self, leb_mode):
        """FLG_XYZ for Pipeline A returns actual Cartesian, not spherical."""
        from libephemeris.constants import SUN, FLG_XYZ
        from libephemeris.planets import calc_ut
        from libephemeris.time_utils import julday

        jd = julday(2024, 6, 15, 12.0)
        pos_xyz, _ = calc_ut(jd, SUN, FLG_XYZ)
        pos_sph, _ = calc_ut(jd, SUN, 0)
        # XYZ distance should be ~1 AU (not degrees-sized values)
        import math
        r = math.sqrt(pos_xyz[0] ** 2 + pos_xyz[1] ** 2 + pos_xyz[2] ** 2)
        assert 0.98 < r < 1.02, f"r={r} not ~1 AU"
        # dist (spherical) should be ~1 AU
        assert 0.98 < pos_sph[2] < 1.02

    def test_xyz_pipeline_b_mean_node(self, leb_mode):
        """FLG_XYZ for Pipeline B (Mean Node) returns Cartesian."""
        from libephemeris.constants import MEAN_NODE, FLG_XYZ
        from libephemeris.planets import calc_ut
        from libephemeris.time_utils import julday

        jd = julday(2024, 6, 15, 12.0)
        pos_xyz, _ = calc_ut(jd, MEAN_NODE, FLG_XYZ)
        # Mean Node distance ~0.00257 AU, Cartesian components should be small
        import math
        r = math.sqrt(pos_xyz[0] ** 2 + pos_xyz[1] ** 2 + pos_xyz[2] ** 2)
        assert 0.001 < r < 0.01, f"r={r} unexpected for Mean Node"

    def test_xyz_radians_pipeline_a(self, leb_mode):
        """FLG_XYZ | FLG_RADIANS should NOT apply radians to Cartesian."""
        from libephemeris.constants import SUN, FLG_XYZ, FLG_RADIANS
        from libephemeris.planets import calc_ut
        from libephemeris.time_utils import julday

        jd = julday(2024, 6, 15, 12.0)
        pos_xyz, _ = calc_ut(jd, SUN, FLG_XYZ)
        pos_xyz_rad, _ = calc_ut(jd, SUN, FLG_XYZ | FLG_RADIANS)
        # Both should be identical — RADIANS is a no-op with XYZ
        assert abs(pos_xyz[0] - pos_xyz_rad[0]) < 1e-12
        assert abs(pos_xyz[1] - pos_xyz_rad[1]) < 1e-12

    def test_sol_eclipse_when_glob_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.eclipse import sol_eclipse_when_glob
        from libephemeris.time_utils import julday

        jd = julday(2024, 1, 1, 0.0)
        ecl_type, tret = sol_eclipse_when_glob(jd)
        assert state._PLANETS is None
        assert ecl_type != 0
        assert tret[0] > jd

    def test_lun_eclipse_when_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.eclipse import lun_eclipse_when
        from libephemeris.time_utils import julday

        jd = julday(2024, 1, 1, 0.0)
        ecl_type, _tret = lun_eclipse_when(jd)
        assert state._PLANETS is None
        assert ecl_type != 0

    def test_swe_sol_eclipse_how_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.eclipse import sol_eclipse_when_glob, sol_eclipse_how
        from libephemeris.time_utils import julday

        jd = julday(2024, 1, 1, 0.0)
        _, tret = sol_eclipse_when_glob(jd)
        jd_max = tret[0]
        _retflag, _attr = sol_eclipse_how(jd_max, [0.0, 0.0, 0.0])
        assert state._PLANETS is None

    def test_rise_trans_sun_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.constants import SUN, CALC_RISE
        from libephemeris.eclipse import rise_trans
        from libephemeris.time_utils import julday

        jd = julday(2024, 6, 15, 0.0)
        _retflag, tret = rise_trans(jd, SUN, CALC_RISE, [12.5, 41.9, 0])
        assert state._PLANETS is None
        assert tret[0] > jd

    def test_rise_trans_star_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.constants import CALC_RISE
        from libephemeris.eclipse import rise_trans
        from libephemeris.time_utils import julday

        jd = julday(2024, 6, 15, 0.0)
        _retflag, tret = rise_trans(jd, "Sirius", CALC_RISE, [12.5, 41.9, 0])
        assert state._PLANETS is None
        assert tret[0] > jd
