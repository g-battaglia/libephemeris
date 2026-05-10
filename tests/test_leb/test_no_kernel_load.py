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
    for target in targets:
        try:
            p = patch(target, side_effect=_block_get_planets)
            p.start()
            patchers.append(p)
        except AttributeError:
            pass

    yield

    for p in patchers:
        p.stop()

    state._PLANETS = None
    libephemeris.set_leb_file("")


class TestFixedStarsNoKernel:
    """Phase 1: swe_fixstar_ut and swe_batch_fixstars_ut must not load DE kernel."""

    def test_swe_fixstar_ut_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.fixed_stars import swe_fixstar_ut
        from libephemeris.time_utils import swe_julday

        jd = swe_julday(2024, 6, 15, 12.0)
        pos, _name, _flags = swe_fixstar_ut("Regulus", jd, 0)
        assert state._PLANETS is None
        assert 140 < pos[0] < 160

    def test_swe_fixstar_ut_with_speed(self, leb_mode):
        from libephemeris import state
        from libephemeris.constants import SEFLG_SPEED
        from libephemeris.fixed_stars import swe_fixstar_ut
        from libephemeris.time_utils import swe_julday

        jd = swe_julday(2024, 6, 15, 12.0)
        pos, _name, _flags = swe_fixstar_ut("Regulus", jd, SEFLG_SPEED)
        assert state._PLANETS is None
        assert pos[3] != 0.0

    def test_swe_batch_fixstars_ut_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.fixed_stars import swe_batch_fixstars_ut
        from libephemeris.time_utils import swe_julday

        jd = swe_julday(2024, 6, 15, 12.0)
        results = swe_batch_fixstars_ut(["Regulus", "Spica", "Aldebaran"], jd, 0)
        assert state._PLANETS is None
        assert all(r is not None for r in results)

    def test_swe_fixstar_ut_j2000(self, leb_mode):
        from libephemeris import state
        from libephemeris.constants import SEFLG_J2000
        from libephemeris.fixed_stars import swe_fixstar_ut
        from libephemeris.time_utils import swe_julday

        jd = swe_julday(2024, 6, 15, 12.0)
        _pos, _, _ = swe_fixstar_ut("Regulus", jd, SEFLG_J2000)
        assert state._PLANETS is None

    def test_swe_fixstar_ut_noaberr(self, leb_mode):
        from libephemeris import state
        from libephemeris.constants import SEFLG_NOABERR
        from libephemeris.fixed_stars import swe_fixstar_ut
        from libephemeris.time_utils import swe_julday

        jd = swe_julday(2024, 6, 15, 12.0)
        _pos, _, _ = swe_fixstar_ut("Regulus", jd, SEFLG_NOABERR)
        assert state._PLANETS is None


class TestPhenoNoKernel:
    """Phase 3: swe_pheno_ut must not load DE kernel."""

    def test_swe_pheno_ut_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.constants import SE_MARS
        from libephemeris.planets import swe_pheno_ut
        from libephemeris.time_utils import swe_julday

        jd = swe_julday(2024, 6, 15, 12.0)
        result = swe_pheno_ut(jd, SE_MARS, 0)
        assert state._PLANETS is None
        assert result[2] > 0  # elongation > 0

    def test_swe_pheno_ut_moon(self, leb_mode):
        from libephemeris import state
        from libephemeris.constants import SE_MOON
        from libephemeris.planets import swe_pheno_ut
        from libephemeris.time_utils import swe_julday

        jd = swe_julday(2024, 6, 15, 12.0)
        result = swe_pheno_ut(jd, SE_MOON, 0)
        assert state._PLANETS is None
        assert 0 < result[0] < 180  # phase angle

    def test_swe_pheno_ut_sun(self, leb_mode):
        from libephemeris import state
        from libephemeris.constants import SE_SUN
        from libephemeris.planets import swe_pheno_ut
        from libephemeris.time_utils import swe_julday

        jd = swe_julday(2024, 6, 15, 12.0)
        result = swe_pheno_ut(jd, SE_SUN, 0)
        assert state._PLANETS is None
        assert result[0] == 0.0  # phase angle for Sun = 0


class TestHeliacalNoKernel:
    """Phase 3: heliacal functions must not load DE kernel."""

    def test_heliacal_pheno_ut_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.constants import SE_VENUS
        from libephemeris.heliacal import heliacal_pheno_ut
        from libephemeris.time_utils import swe_julday

        jd = swe_julday(2024, 6, 15, 12.0)
        _result, _retflag = heliacal_pheno_ut(jd, 41.9, 12.5, body=SE_VENUS)
        assert state._PLANETS is None


class TestEclipseNoKernel:
    """Phase 4: eclipse core functions must not load DE kernel."""

    def test_calc_gamma_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.eclipse import _calc_gamma
        from libephemeris.time_utils import swe_julday

        jd = swe_julday(2024, 4, 8, 18.0)
        gamma = _calc_gamma(jd)
        assert state._PLANETS is None
        assert isinstance(gamma, float)

    def test_calc_penumbra_limit_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.eclipse import _calc_penumbra_limit
        from libephemeris.time_utils import swe_julday

        jd = swe_julday(2024, 4, 8, 18.0)
        l1 = _calc_penumbra_limit(jd)
        assert state._PLANETS is None
        assert l1 > 0

    def test_besselian_x_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.eclipse import calc_besselian_x
        from libephemeris.time_utils import swe_julday

        jd = swe_julday(2024, 4, 8, 18.0)
        x = calc_besselian_x(jd)
        assert state._PLANETS is None
        assert isinstance(x, float)

    def test_xyz_pipeline_a_no_double_convert(self, leb_mode):
        """SEFLG_XYZ for Pipeline A returns actual Cartesian, not spherical."""
        from libephemeris.constants import SE_SUN, SEFLG_XYZ
        from libephemeris.planets import swe_calc_ut
        from libephemeris.time_utils import swe_julday

        jd = swe_julday(2024, 6, 15, 12.0)
        pos_xyz, _ = swe_calc_ut(jd, SE_SUN, SEFLG_XYZ)
        pos_sph, _ = swe_calc_ut(jd, SE_SUN, 0)
        # XYZ distance should be ~1 AU (not degrees-sized values)
        import math
        r = math.sqrt(pos_xyz[0] ** 2 + pos_xyz[1] ** 2 + pos_xyz[2] ** 2)
        assert 0.98 < r < 1.02, f"r={r} not ~1 AU"
        # dist (spherical) should be ~1 AU
        assert 0.98 < pos_sph[2] < 1.02

    def test_xyz_pipeline_b_mean_node(self, leb_mode):
        """SEFLG_XYZ for Pipeline B (Mean Node) returns Cartesian."""
        from libephemeris.constants import SE_MEAN_NODE, SEFLG_XYZ
        from libephemeris.planets import swe_calc_ut
        from libephemeris.time_utils import swe_julday

        jd = swe_julday(2024, 6, 15, 12.0)
        pos_xyz, _ = swe_calc_ut(jd, SE_MEAN_NODE, SEFLG_XYZ)
        # Mean Node distance ~0.00257 AU, Cartesian components should be small
        import math
        r = math.sqrt(pos_xyz[0] ** 2 + pos_xyz[1] ** 2 + pos_xyz[2] ** 2)
        assert 0.001 < r < 0.01, f"r={r} unexpected for Mean Node"

    def test_xyz_radians_pipeline_a(self, leb_mode):
        """SEFLG_XYZ | SEFLG_RADIANS should NOT apply radians to Cartesian."""
        from libephemeris.constants import SE_SUN, SEFLG_XYZ, SEFLG_RADIANS
        from libephemeris.planets import swe_calc_ut
        from libephemeris.time_utils import swe_julday

        jd = swe_julday(2024, 6, 15, 12.0)
        pos_xyz, _ = swe_calc_ut(jd, SE_SUN, SEFLG_XYZ)
        pos_xyz_rad, _ = swe_calc_ut(jd, SE_SUN, SEFLG_XYZ | SEFLG_RADIANS)
        # Both should be identical — RADIANS is a no-op with XYZ
        assert abs(pos_xyz[0] - pos_xyz_rad[0]) < 1e-12
        assert abs(pos_xyz[1] - pos_xyz_rad[1]) < 1e-12

    def test_sol_eclipse_when_glob_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.eclipse import sol_eclipse_when_glob
        from libephemeris.time_utils import swe_julday

        jd = swe_julday(2024, 1, 1, 0.0)
        ecl_type, tret = sol_eclipse_when_glob(jd)
        assert state._PLANETS is None
        assert ecl_type != 0
        assert tret[0] > jd

    def test_lun_eclipse_when_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.eclipse import lun_eclipse_when
        from libephemeris.time_utils import swe_julday

        jd = swe_julday(2024, 1, 1, 0.0)
        ecl_type, tret = lun_eclipse_when(jd)
        assert state._PLANETS is None
        assert ecl_type != 0

    def test_swe_sol_eclipse_how_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.eclipse import sol_eclipse_when_glob, swe_sol_eclipse_how
        from libephemeris.time_utils import swe_julday

        jd = swe_julday(2024, 1, 1, 0.0)
        _, tret = sol_eclipse_when_glob(jd)
        jd_max = tret[0]
        retflag, attr = swe_sol_eclipse_how(jd_max, [0.0, 0.0, 0.0])
        assert state._PLANETS is None

    def test_rise_trans_sun_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.constants import SE_SUN, SE_CALC_RISE
        from libephemeris.eclipse import rise_trans
        from libephemeris.time_utils import swe_julday

        jd = swe_julday(2024, 6, 15, 0.0)
        retflag, tret = rise_trans(jd, SE_SUN, SE_CALC_RISE, [12.5, 41.9, 0])
        assert state._PLANETS is None
        assert tret[0] > jd

    def test_rise_trans_star_no_kernel(self, leb_mode):
        from libephemeris import state
        from libephemeris.constants import SE_CALC_RISE
        from libephemeris.eclipse import rise_trans
        from libephemeris.time_utils import swe_julday

        jd = swe_julday(2024, 6, 15, 0.0)
        retflag, tret = rise_trans(jd, "Sirius", SE_CALC_RISE, [12.5, 41.9, 0])
        assert state._PLANETS is None
        assert tret[0] > jd
