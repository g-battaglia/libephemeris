"""
Tests for tidal acceleration functions (set_tid_acc, get_tid_acc).

Tidal acceleration affects Delta T calculations for historical dates.
"""

import pytest
import libephemeris as ephem


class TestTidAccBasicFunctionality:
    """Test basic get/set functionality for tidal acceleration."""

    @pytest.fixture(autouse=True)
    def reset_tid_acc(self):
        """Reset tidal acceleration to default before and after each test."""
        # Reset to automatic/default before test
        ephem.set_tid_acc(ephem.TIDAL_AUTOMATIC)
        yield
        # Reset after test to avoid affecting other tests
        ephem.set_tid_acc(ephem.TIDAL_AUTOMATIC)

    @pytest.mark.unit
    def test_get_tid_acc_default(self):
        """get_tid_acc returns the named automatic zero point."""
        value = ephem.get_tid_acc()
        assert value == ephem.TIDAL_DEFAULT

    @pytest.mark.unit
    def test_set_tid_acc_de421(self):
        """set_tid_acc should accept DE421 tidal acceleration."""
        ephem.set_tid_acc(ephem.TIDAL_DE421)
        assert ephem.get_tid_acc() == ephem.TIDAL_DE421

    @pytest.mark.unit
    def test_set_tid_acc_de431(self):
        """set_tid_acc should accept DE431 tidal acceleration."""
        ephem.set_tid_acc(ephem.TIDAL_DE431)
        assert ephem.get_tid_acc() == ephem.TIDAL_DE431

    @pytest.mark.unit
    def test_set_tid_acc_de441(self):
        """set_tid_acc should accept DE441 tidal acceleration."""
        ephem.set_tid_acc(ephem.TIDAL_DE441)
        assert ephem.get_tid_acc() == ephem.TIDAL_DE441

    @pytest.mark.unit
    def test_set_tid_acc_custom_value(self):
        """set_tid_acc should accept custom numeric values."""
        custom_value = -26.5
        ephem.set_tid_acc(custom_value)
        assert ephem.get_tid_acc() == custom_value

    @pytest.mark.unit
    def test_set_tid_acc_automatic_resets_to_default(self):
        """TIDAL_AUTOMATIC resets to the public automatic default."""
        # First set a custom value
        ephem.set_tid_acc(ephem.TIDAL_DE421)
        assert ephem.get_tid_acc() == ephem.TIDAL_DE421

        # Reset with automatic
        ephem.set_tid_acc(ephem.TIDAL_AUTOMATIC)
        assert ephem.get_tid_acc() == ephem.TIDAL_DEFAULT

    @pytest.mark.unit
    def test_set_tid_acc_zero_stores_zero(self):
        """set_tid_acc with 0.0 should store 0.0 (reference-API compat)."""
        ephem.set_tid_acc(ephem.TIDAL_DE421)
        ephem.set_tid_acc(0.0)
        assert ephem.get_tid_acc() == 0.0


class TestTidAccConstants:
    """Test that all tidal acceleration constants are properly defined."""

    @pytest.mark.unit
    def test_se_tidal_constants_exist(self):
        """All TIDAL_* constants should exist."""
        assert hasattr(ephem, "TIDAL_DE200")
        assert hasattr(ephem, "TIDAL_DE403")
        assert hasattr(ephem, "TIDAL_DE404")
        assert hasattr(ephem, "TIDAL_DE405")
        assert hasattr(ephem, "TIDAL_DE406")
        assert hasattr(ephem, "TIDAL_DE421")
        assert hasattr(ephem, "TIDAL_DE422")
        assert hasattr(ephem, "TIDAL_DE430")
        assert hasattr(ephem, "TIDAL_DE431")
        assert hasattr(ephem, "TIDAL_DE440")
        assert hasattr(ephem, "TIDAL_DE441")
        assert hasattr(ephem, "TIDAL_DEFAULT")
        assert hasattr(ephem, "TIDAL_AUTOMATIC")

    @pytest.mark.unit
    def test_tidal_constants_aliases_exist(self):
        """TIDAL_* aliases (without SE_ prefix) should exist."""
        assert hasattr(ephem, "TIDAL_DE200")
        assert hasattr(ephem, "TIDAL_DE403")
        assert hasattr(ephem, "TIDAL_DE421")
        assert hasattr(ephem, "TIDAL_DE431")
        assert hasattr(ephem, "TIDAL_DE440")
        assert hasattr(ephem, "TIDAL_DE441")
        assert hasattr(ephem, "TIDAL_DEFAULT")
        assert hasattr(ephem, "TIDAL_AUTOMATIC")

    @pytest.mark.unit
    def test_tidal_constants_values_are_negative(self):
        """Tidal acceleration values should be negative."""
        assert ephem.TIDAL_DE200 < 0
        assert ephem.TIDAL_DE421 < 0
        assert ephem.TIDAL_DE431 < 0
        assert ephem.TIDAL_DE440 < 0
        assert ephem.TIDAL_DE441 < 0

    @pytest.mark.unit
    def test_tidal_automatic_is_nonphysical_sentinel(self):
        """TIDAL_AUTOMATIC is a nonphysical sentinel, not an acceleration."""
        assert isinstance(ephem.TIDAL_AUTOMATIC, int)
        assert ephem.TIDAL_AUTOMATIC > 0

    @pytest.mark.unit
    def test_tidal_default_is_de431(self):
        """The named historical default uses the published DE431 solution."""
        assert ephem.TIDAL_DEFAULT == ephem.TIDAL_DE431


class TestTidAccFunctionAliases:
    """Test that function aliases work correctly."""

    @pytest.fixture(autouse=True)
    def reset_tid_acc(self):
        """Reset tidal acceleration to default before and after each test."""
        ephem.set_tid_acc(ephem.TIDAL_AUTOMATIC)
        yield
        ephem.set_tid_acc(ephem.TIDAL_AUTOMATIC)

    @pytest.mark.unit
    def test_swe_prefixed_aliases_exist(self):
        """set_tid_acc and get_tid_acc should exist."""
        assert hasattr(ephem, "set_tid_acc")
        assert hasattr(ephem, "get_tid_acc")

    @pytest.mark.unit
    def test_swe_prefixed_aliases_work(self):
        """set_tid_acc and get_tid_acc should work identically."""
        ephem.set_tid_acc(ephem.TIDAL_DE421)
        assert ephem.get_tid_acc() == ephem.TIDAL_DE421

        # Cross-check with non-prefixed versions
        assert ephem.get_tid_acc() == ephem.TIDAL_DE421

    @pytest.mark.unit
    def test_non_prefixed_and_prefixed_are_same(self):
        """set_tid_acc and set_tid_acc should be the same function."""


class TestTidAccReturnTypes:
    """Test return types of tidal acceleration functions."""

    @pytest.fixture(autouse=True)
    def reset_tid_acc(self):
        """Reset tidal acceleration to default."""
        ephem.set_tid_acc(ephem.TIDAL_AUTOMATIC)
        yield
        ephem.set_tid_acc(ephem.TIDAL_AUTOMATIC)

    @pytest.mark.unit
    def test_get_tid_acc_returns_float(self):
        """get_tid_acc should return a float."""
        result = ephem.get_tid_acc()
        assert isinstance(result, float)

    @pytest.mark.unit
    def test_set_tid_acc_returns_none(self):
        """set_tid_acc should return None."""
        result = ephem.set_tid_acc(ephem.TIDAL_DE421)
        assert result is None


# ---------------------------------------------------------------------------
# Delta T adjustment for a non-default tidal acceleration
# ---------------------------------------------------------------------------
# The conversion is derived independently from the IERS Conventions (2010)
# mean-lunar-longitude rate.  A secular-acceleration difference accumulates
# 0.5*dacc*T**2 arcseconds; dividing by lunar mean motion converts that phase
# displacement into seconds of time.  Historical reductions use epoch 1955.5
# (Morrison & Stephenson 2004; LUNAR97).

# Sample dates before the 1955.5 reference epoch (years -3000 .. 1950)
_PRE_1955_YEARS = [-3000, -1000, 0, 1000, 1600, 1800, 1900, 1950]
# Sample dates after the 1955.5 epoch.
_POST_1955_YEARS = [1956, 1980, 2000, 2020, 2100, 3000]


def _jd(year):
    """Julian Day of Jan 1, 0h of the given year (Gregorian)."""
    return ephem.julday(year, 1, 1, 0.0)


class TestTidAccDeltaTAdjustment:
    """Physical Delta T rescaling for non-default tidal accelerations."""

    @pytest.fixture(autouse=True)
    def reset_state(self):
        """Reset tidal acceleration and userdef Delta T around each test."""
        ephem.set_tid_acc(ephem.TIDAL_AUTOMATIC)
        ephem.set_delta_t_userdef(None)
        yield
        ephem.set_tid_acc(ephem.TIDAL_AUTOMATIC)
        ephem.set_delta_t_userdef(None)

    @pytest.mark.unit
    def test_default_deltat_bit_exact_unchanged(self):
        """With the default tidal acceleration Delta T must be unchanged.

        The adjustment must be exactly 0.0 (not merely small) at the
        default, so the default Delta T path is bit-for-bit identical to
        the unadjusted model output.
        """
        from libephemeris.time_utils import _tid_acc_adjustment_seconds

        for year in _PRE_1955_YEARS + _POST_1955_YEARS:
            jd = _jd(year)
            assert _tid_acc_adjustment_seconds(jd, ephem.get_tid_acc()) == 0.0
            base = ephem.deltat(jd)
            ephem.set_tid_acc(ephem.TIDAL_DEFAULT)  # explicit default value
            assert ephem.deltat(jd) == base
            ephem.set_tid_acc(ephem.TIDAL_AUTOMATIC)
            assert ephem.deltat(jd) == base

    @pytest.mark.unit
    def test_adjustment_follows_iers_lunar_motion_conversion(self):
        """Delta T deltas follow phase acceleration divided by mean motion."""
        from libephemeris.time_utils import _TID_ACC_REFERENCE_JD

        # IERS TN 36, Eq. 5.43 linear terms for F=L-Omega and Omega,
        # arcseconds per Julian century.
        f_rate = 1739527262.8478
        omega_rate = -6962890.5431
        lunar_motion_arcsec_day = (f_rate + omega_rate) / 36525.0
        seconds_per_arcsec = 86400.0 / lunar_motion_arcsec_day

        for tid_acc in (-25.58, -26.5, -22.0):
            for year in _PRE_1955_YEARS:
                jd = _jd(year)
                ephem.set_tid_acc(ephem.TIDAL_AUTOMATIC)
                base = ephem.deltat(jd) * 86400.0
                ephem.set_tid_acc(tid_acc)
                delta = ephem.deltat(jd) * 86400.0 - base
                u = (jd - _TID_ACC_REFERENCE_JD) / 36525.0
                expected = (
                    -0.5 * seconds_per_arcsec * (tid_acc - ephem.TIDAL_DEFAULT) * u * u
                )
                assert delta == pytest.approx(expected, abs=1e-6), (
                    f"tid_acc={tid_acc}, year={year}"
                )

    @pytest.mark.unit
    def test_adjustment_zero_from_reference_epoch_onwards(self):
        """Modern Delta T does not depend on the historical adjustment."""
        for year in _POST_1955_YEARS:
            jd = _jd(year)
            ephem.set_tid_acc(ephem.TIDAL_AUTOMATIC)
            base = ephem.deltat(jd)
            ephem.set_tid_acc(-22.0)  # large offset from default
            assert ephem.deltat(jd) == base, f"year={year}"

    @pytest.mark.unit
    def test_userdef_delta_t_never_adjusted(self):
        """A user-defined Delta T is a hard override, never rescaled."""
        ephem.set_delta_t_userdef(0.3)  # days
        ephem.set_tid_acc(-22.0)
        jd = _jd(-1000)
        assert ephem.deltat(jd) == 0.3
        assert ephem.deltat_ex(jd, ephem.FLG_MOSEPH) == 0.3


class TestDeltatExEphemerisFlag:
    """deltat_ex() must honour the ephemeris flag's tidal acceleration."""

    @pytest.fixture(autouse=True)
    def reset_state(self):
        """Reset tidal acceleration and userdef Delta T around each test."""
        ephem.set_tid_acc(ephem.TIDAL_AUTOMATIC)
        ephem.set_delta_t_userdef(None)
        yield
        ephem.set_tid_acc(ephem.TIDAL_AUTOMATIC)
        ephem.set_delta_t_userdef(None)

    @pytest.mark.unit
    def test_moseph_uses_moshier_tidal_acceleration(self):
        """FLG_MOSEPH must equal an explicit set_tid_acc(TIDAL_MOSEPH)."""
        for year in (-3000, -1000, 0, 1600, 1900):
            jd = _jd(year)
            via_flag = ephem.deltat_ex(jd, ephem.FLG_MOSEPH)
            ephem.set_tid_acc(ephem.TIDAL_MOSEPH)  # -25.58
            via_setting = ephem.deltat(jd)
            ephem.set_tid_acc(ephem.TIDAL_AUTOMATIC)
            assert via_flag == via_setting, f"year={year}"
            # And it must differ from the default for pre-1955 dates
            assert via_flag != ephem.deltat(jd), f"year={year}"

    @pytest.mark.unit
    def test_moseph_modern_dates_unchanged(self):
        """FLG_MOSEPH matches deltat() for dates from 1955.0 on."""
        for year in (1956, 2000, 2020, 2100):
            jd = _jd(year)
            assert ephem.deltat_ex(jd, ephem.FLG_MOSEPH) == ephem.deltat(jd)

    @pytest.mark.unit
    def test_moseph_wins_over_swieph(self):
        """FLG_MOSEPH takes precedence when both bits are set."""
        jd = _jd(-1000)
        both = ephem.FLG_MOSEPH | ephem.FLG_SWIEPH
        assert ephem.deltat_ex(jd, both) == ephem.deltat_ex(jd, ephem.FLG_MOSEPH)

    @pytest.mark.unit
    def test_moseph_wins_over_jpleph(self):
        """FLG_MOSEPH takes precedence over FLG_JPLEPH."""
        jd = _jd(-1000)
        both = ephem.FLG_MOSEPH | ephem.FLG_JPLEPH
        assert ephem.deltat_ex(jd, both) == ephem.deltat_ex(jd, ephem.FLG_MOSEPH)

    @pytest.mark.unit
    def test_user_set_tid_acc_wins_over_moseph(self):
        """An explicit set_tid_acc() overrides the MOSEPH-implied value."""
        jd = _jd(-1000)
        ephem.set_tid_acc(-22.0)
        assert ephem.deltat_ex(jd, ephem.FLG_MOSEPH) == ephem.deltat(jd)

    @pytest.mark.unit
    def test_swieph_and_jpleph_match_deltat(self):
        """Non-MOSEPH ephemeris flags return the same Delta T as deltat()."""
        for year in (-1000, 1800, 2000):
            jd = _jd(year)
            base = ephem.deltat(jd)
            assert ephem.deltat_ex(jd, ephem.FLG_SWIEPH) == base
            assert ephem.deltat_ex(jd, ephem.FLG_JPLEPH) == base

    @pytest.mark.unit
    @pytest.mark.parametrize("year", [-1000, 1900, 2000])
    def test_moseph_query_does_not_mutate_tid_acc(self, year):
        """deltat_ex is pure: a MOSEPH query never leaks into get_tid_acc().

        The MOSEPH-implied acceleration is applied functionally to the
        returned Delta T only; the exposed state changes exclusively through
        set_tid_acc().
        """
        ephem.deltat_ex(_jd(year), ephem.FLG_MOSEPH)
        assert ephem.get_tid_acc() == ephem.TIDAL_DEFAULT

    @pytest.mark.unit
    @pytest.mark.parametrize("flag", [0, ephem.FLG_SWIEPH, ephem.FLG_JPLEPH])
    def test_no_deltat_query_mutates_tid_acc(self, flag):
        """Neither deltat_ex (any flag) nor deltat touches get_tid_acc()."""
        jd = _jd(1900)
        ephem.deltat_ex(jd, ephem.FLG_MOSEPH)
        ephem.deltat_ex(jd, flag)
        ephem.deltat(jd)
        assert ephem.get_tid_acc() == ephem.TIDAL_DEFAULT

    @pytest.mark.unit
    def test_moseph_adjustment_is_functional_only(self):
        """The MOSEPH acceleration adjusts the returned value, not the state."""
        jd = _jd(1900)
        via_flag = ephem.deltat_ex(jd, ephem.FLG_MOSEPH)
        assert ephem.get_tid_acc() == ephem.TIDAL_DEFAULT
        ephem.set_tid_acc(ephem.TIDAL_MOSEPH)
        assert via_flag == ephem.deltat(jd)

    @pytest.mark.unit
    def test_user_defined_delta_t_wins_and_state_is_untouched(self):
        """A hard Delta-T override wins; the exposed acceleration is stable."""
        jd = _jd(1900)
        ephem.deltat_ex(jd, ephem.FLG_MOSEPH)
        ephem.set_delta_t_userdef(0.123)

        assert ephem.deltat_ex(jd, ephem.FLG_SWIEPH) == 0.123
        assert ephem.deltat(jd) == 0.123
        assert ephem.get_tid_acc() == ephem.TIDAL_DEFAULT

    @pytest.mark.unit
    def test_explicit_tid_acc_is_not_reselected(self):
        """An explicit tidal acceleration wins over every ephemeris flag."""
        jd = _jd(1900)
        ephem.set_tid_acc(-22.0)

        for flag in (ephem.FLG_MOSEPH, ephem.FLG_SWIEPH, ephem.FLG_JPLEPH):
            ephem.deltat_ex(jd, flag)
            assert ephem.get_tid_acc() == -22.0
