"""
Tests for CompositeLEBReader.

Verifies from_directory, from_file_with_companions, body dispatch,
eval_body across multiple readers, jd_range merging, error handling,
and context manager protocol.
"""

from __future__ import annotations

import glob
import math
import os

import pytest

from libephemeris.leb_reader import open_leb
from libephemeris.leb_composite import CompositeLEBReader
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
    CHIRON,
    CERES,
    PALLAS,
    JUNO,
    VESTA,
    OSCU_APOG,
    INTP_APOG,
    INTP_PERG,
)

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

LEB1_BASE = os.path.join(PROJECT_ROOT, "data", "leb", "ephemeris_base.leb")
LEB2_DIR = os.path.join(PROJECT_ROOT, "data", "leb2")
LEB2_BASE_CORE = os.path.join(LEB2_DIR, "base_core.leb2")
LEB2_BASE_ASTEROIDS = os.path.join(LEB2_DIR, "base_asteroids.leb2")
LEB2_BASE_APOGEE = os.path.join(LEB2_DIR, "base_apogee.leb2")
LEB2_BASE_URANIANS = os.path.join(LEB2_DIR, "base_uranians.leb2")

SKIP_NO_LEB1 = pytest.mark.skipif(
    not os.path.exists(LEB1_BASE), reason="LEB1 base file not found"
)
SKIP_NO_LEB2_DIR = pytest.mark.skipif(
    not os.path.exists(LEB2_DIR), reason="LEB2 directory not found"
)
SKIP_NO_LEB2_CORE = pytest.mark.skipif(
    not os.path.exists(LEB2_BASE_CORE), reason="LEB2 base core not found"
)

JD_J2000 = 2451545.0  # 2000-01-01 12:00 TT


@pytest.fixture(scope="module")
def single_tier_dir(tmp_path_factory):
    """A directory holding exactly one tier's LEB2 group files.

    ``data/leb2`` mixes tiers (base/medium/extended), which ``from_directory``
    now rejects. The medium set is the only complete single tier present
    (medium_core/asteroids/apogee/uranians), so symlink those into a temp dir.
    """
    src = sorted(glob.glob(os.path.join(LEB2_DIR, "medium_*.leb2")))
    if not src:
        pytest.skip("medium_*.leb2 files not found")
    d = tmp_path_factory.mktemp("leb2_single_tier")
    for s in src:
        os.symlink(s, os.path.join(d, os.path.basename(s)))
    return str(d)


# ============================================================================
# from_directory
# ============================================================================


@SKIP_NO_LEB2_DIR
class TestFromDirectory:
    """Test CompositeLEBReader.from_directory."""

    @pytest.mark.unit
    def test_from_directory_opens_all_files(self, single_tier_dir):
        """from_directory discovers and opens all .leb files in a directory."""
        with CompositeLEBReader.from_directory(single_tier_dir) as reader:
            # Should have opened all the medium group files
            assert reader is not None
            # Core bodies must be present
            assert reader.has_body(SUN)
            assert reader.has_body(MOON)

    @pytest.mark.unit
    def test_from_directory_mixed_tiers_raises(self):
        """from_directory refuses to merge files from different tiers."""
        # data/leb2 holds base_*, medium_*, and extended_* files together;
        # merging them first-wins would silently corrupt positions.
        with pytest.raises(ValueError, match="mixed tiers"):
            CompositeLEBReader.from_directory(LEB2_DIR)

    @pytest.mark.unit
    def test_from_directory_has_core_bodies(self, single_tier_dir):
        """Core planets are accessible from composite reader."""
        with CompositeLEBReader.from_directory(single_tier_dir) as reader:
            for body in [
                SUN,
                MOON,
                MERCURY,
                VENUS,
                MARS,
                JUPITER,
                SATURN,
                URANUS,
                NEPTUNE,
                PLUTO,
            ]:
                assert reader.has_body(body), f"Missing core body {body}"

    @pytest.mark.unit
    def test_from_directory_has_asteroid_bodies(self, single_tier_dir):
        """Asteroid bodies are accessible when asteroids group is present."""
        with CompositeLEBReader.from_directory(single_tier_dir) as reader:
            for body in [CHIRON, CERES, PALLAS, JUNO, VESTA]:
                assert reader.has_body(body), f"Missing asteroid body {body}"

    @pytest.mark.unit
    def test_from_directory_has_apogee_bodies(self, single_tier_dir):
        """Apogee bodies are accessible when apogee group is present."""
        with CompositeLEBReader.from_directory(single_tier_dir) as reader:
            for body in [OSCU_APOG, INTP_APOG, INTP_PERG]:
                assert reader.has_body(body), f"Missing apogee body {body}"

    @pytest.mark.unit
    def test_from_directory_eval_body_returns_valid(self, single_tier_dir):
        """eval_body returns valid position/velocity tuples."""
        with CompositeLEBReader.from_directory(single_tier_dir) as reader:
            pos, vel = reader.eval_body(SUN, JD_J2000)
            assert len(pos) == 3
            assert len(vel) == 3
            # eval_body returns raw Chebyshev coords (may be cartesian AU or ecliptic)
            # Just verify finite values
            for v in pos + vel:
                assert math.isfinite(v)

    @pytest.mark.unit
    def test_from_directory_nonexistent_raises(self):
        """from_directory raises FileNotFoundError for missing directory."""
        with pytest.raises(FileNotFoundError):
            CompositeLEBReader.from_directory("/nonexistent/dir/path")

    @pytest.mark.unit
    def test_from_directory_empty_dir_raises(self, tmp_path):
        """from_directory raises FileNotFoundError if no .leb files found."""
        with pytest.raises(FileNotFoundError):
            CompositeLEBReader.from_directory(str(tmp_path))

    @pytest.mark.unit
    def test_from_directory_invalid_files_only_raises(self, tmp_path):
        """from_directory raises ValueError if all .leb files are invalid."""
        bad_file = tmp_path / "bad.leb"
        bad_file.write_bytes(b"not valid leb data here!!")
        with pytest.raises((ValueError, Exception)):
            CompositeLEBReader.from_directory(str(tmp_path))

    @pytest.mark.unit
    def test_from_directory_jd_range_is_widest(self, single_tier_dir):
        """jd_range is the widest range across all component readers."""
        with CompositeLEBReader.from_directory(single_tier_dir) as reader:
            jd_start, jd_end = reader.jd_range
            assert jd_start < JD_J2000 < jd_end
            # Medium tier covers ~1550-2650, so range is well over 100 years
            assert (jd_end - jd_start) > 365.25 * 100

    @pytest.mark.unit
    def test_from_directory_path_property(self, single_tier_dir):
        """path property returns first reader's path."""
        with CompositeLEBReader.from_directory(single_tier_dir) as reader:
            assert reader.path is not None
            assert reader.path.endswith(".leb2")


# ============================================================================
# from_file_with_companions
# ============================================================================


@SKIP_NO_LEB2_CORE
class TestFromFileWithCompanions:
    """Test CompositeLEBReader.from_file_with_companions."""

    @pytest.mark.unit
    def test_discovers_companion_files(self):
        """from_file_with_companions discovers other base_*.leb files."""
        with CompositeLEBReader.from_file_with_companions(LEB2_BASE_CORE) as reader:
            # Should have core bodies
            assert reader.has_body(SUN)
            assert reader.has_body(MOON)
            # If companion files exist, should also have those bodies
            if os.path.exists(LEB2_BASE_ASTEROIDS):
                assert reader.has_body(CHIRON)

    @pytest.mark.unit
    def test_core_only_when_no_companions(self, tmp_path):
        """Returns single-reader composite when no companions found."""
        import shutil

        # Copy just one file to a temp directory
        isolated = tmp_path / "base_core.leb2"
        shutil.copy2(LEB2_BASE_CORE, str(isolated))

        with CompositeLEBReader.from_file_with_companions(str(isolated)) as reader:
            assert reader.has_body(SUN)
            assert reader.has_body(MOON)
            # Should NOT have asteroid bodies since companion not present
            assert not reader.has_body(CHIRON)

    @pytest.mark.unit
    def test_eval_body_across_groups(self):
        """eval_body dispatches to correct reader for each body group."""
        if not os.path.exists(LEB2_BASE_ASTEROIDS):
            pytest.skip("base_asteroids.leb2 not found")

        with CompositeLEBReader.from_file_with_companions(LEB2_BASE_CORE) as reader:
            # Core body — eval_body returns raw coords (cartesian AU or ecliptic)
            pos_sun, vel_sun = reader.eval_body(SUN, JD_J2000)
            for v in pos_sun + vel_sun:
                assert math.isfinite(v)

            # Asteroid body (from different file)
            pos_chi, vel_chi = reader.eval_body(CHIRON, JD_J2000)
            for v in pos_chi + vel_chi:
                assert math.isfinite(v)

    @pytest.mark.unit
    def test_missing_body_raises_keyerror(self):
        """eval_body raises KeyError for body not in any reader."""
        with CompositeLEBReader.from_file_with_companions(LEB2_BASE_CORE) as reader:
            with pytest.raises(KeyError):
                reader.eval_body(99999, JD_J2000)

    @pytest.mark.unit
    def test_nonexistent_file_raises(self):
        """from_file_with_companions raises for nonexistent file."""
        with pytest.raises((FileNotFoundError, Exception)):
            CompositeLEBReader.from_file_with_companions("/nonexistent.leb")


# ============================================================================
# Body dispatch and eval_body correctness
# ============================================================================


@SKIP_NO_LEB2_DIR
class TestBodyDispatch:
    """Test that bodies are dispatched to correct readers."""

    CORE_BODIES = [
        SUN,
        MOON,
        MERCURY,
        VENUS,
        MARS,
        JUPITER,
        SATURN,
        URANUS,
        NEPTUNE,
        PLUTO,
        MEAN_NODE,
        TRUE_NODE,
        MEAN_APOG,
    ]
    ASTEROID_BODIES = [CHIRON, CERES, PALLAS, JUNO, VESTA]
    APOGEE_BODIES = [OSCU_APOG, INTP_APOG, INTP_PERG]

    TEST_DATES = [
        2451545.0,  # J2000
        2460000.0,  # 2023
        2458849.5,  # 2020-01-01
        2455197.5,  # 2010-01-01
    ]

    @pytest.mark.unit
    @pytest.mark.parametrize("body", CORE_BODIES)
    def test_core_body_eval(self, body, single_tier_dir):
        """Each core body can be evaluated at J2000."""
        with CompositeLEBReader.from_directory(single_tier_dir) as reader:
            if not reader.has_body(body):
                pytest.skip(f"Body {body} not in composite reader")
            pos, vel = reader.eval_body(body, JD_J2000)
            assert len(pos) == 3
            assert len(vel) == 3
            # Position components should be finite
            for v in pos + vel:
                assert math.isfinite(v), f"Non-finite value for body {body}"

    @pytest.mark.unit
    @pytest.mark.parametrize("body", ASTEROID_BODIES)
    def test_asteroid_body_eval(self, body, single_tier_dir):
        """Each asteroid body can be evaluated at J2000."""
        with CompositeLEBReader.from_directory(single_tier_dir) as reader:
            if not reader.has_body(body):
                pytest.skip(f"Body {body} not in composite reader")
            pos, vel = reader.eval_body(body, JD_J2000)
            assert len(pos) == 3
            assert len(vel) == 3
            for v in pos + vel:
                assert math.isfinite(v)

    @pytest.mark.unit
    @pytest.mark.parametrize("body", APOGEE_BODIES)
    def test_apogee_body_eval(self, body, single_tier_dir):
        """Each apogee body can be evaluated at J2000."""
        with CompositeLEBReader.from_directory(single_tier_dir) as reader:
            if not reader.has_body(body):
                pytest.skip(f"Body {body} not in composite reader")
            pos, vel = reader.eval_body(body, JD_J2000)
            assert len(pos) == 3
            assert len(vel) == 3

    @pytest.mark.unit
    @pytest.mark.parametrize("jd", TEST_DATES)
    def test_sun_across_dates(self, jd, single_tier_dir):
        """Sun position varies sensibly across dates."""
        with CompositeLEBReader.from_directory(single_tier_dir) as reader:
            pos, vel = reader.eval_body(SUN, jd)
            # eval_body returns raw cartesian AU coords, not ecliptic degrees
            for v in pos + vel:
                assert math.isfinite(v)
            # Sun geocentric distance should be small (< 0.1 AU for geocentric offset)
            dist = math.sqrt(pos[0] ** 2 + pos[1] ** 2 + pos[2] ** 2)
            assert dist < 0.02  # Sun geocentric Chebyshev is near 0 in geocentric frame

    @pytest.mark.unit
    def test_consistency_with_single_reader(self):
        """Composite reader produces same results as individual reader."""
        if not os.path.exists(LEB2_BASE_CORE):
            pytest.skip("LEB2 core not found")
        with open_leb(LEB2_BASE_CORE) as single:
            with CompositeLEBReader.from_file_with_companions(LEB2_BASE_CORE) as comp:
                for body in [SUN, MOON, MARS]:
                    pos_s, vel_s = single.eval_body(body, JD_J2000)
                    pos_c, vel_c = comp.eval_body(body, JD_J2000)
                    for i in range(3):
                        assert pos_s[i] == pytest.approx(pos_c[i], abs=1e-12)
                        assert vel_s[i] == pytest.approx(vel_c[i], abs=1e-12)


# ============================================================================
# Nutation, Delta-T, and auxiliary data
# ============================================================================


@SKIP_NO_LEB2_DIR
class TestAuxiliaryData:
    """Test nutation, delta_t, and other auxiliary data from composite."""

    @pytest.mark.unit
    def test_eval_nutation(self, single_tier_dir):
        """eval_nutation returns (dpsi, deps) in radians."""
        with CompositeLEBReader.from_directory(single_tier_dir) as reader:
            try:
                dpsi, deps = reader.eval_nutation(JD_J2000)
                # Nutation values are small angles in radians
                assert abs(dpsi) < 0.001  # < ~200 arcsec
                assert abs(deps) < 0.001
            except ValueError:
                pytest.skip("No nutation reader in composite")

    @pytest.mark.unit
    def test_delta_t(self, single_tier_dir):
        """delta_t returns reasonable values."""
        with CompositeLEBReader.from_directory(single_tier_dir) as reader:
            try:
                dt = reader.delta_t(JD_J2000)
                # At J2000, delta-T is about 63.8 seconds = 0.000738 days
                assert 0.0005 < dt < 0.001
            except ValueError:
                pytest.skip("No delta_t reader in composite")

    @pytest.mark.unit
    def test_has_body_false_for_unknown(self, single_tier_dir):
        """has_body returns False for nonexistent body IDs."""
        with CompositeLEBReader.from_directory(single_tier_dir) as reader:
            assert not reader.has_body(99999)
            assert not reader.has_body(-1)


# ============================================================================
# Context manager and close
# ============================================================================


@SKIP_NO_LEB2_DIR
class TestContextManager:
    """Test context manager and close behavior."""

    @pytest.mark.unit
    def test_context_manager_protocol(self, single_tier_dir):
        """Composite reader supports with-statement."""
        with CompositeLEBReader.from_directory(single_tier_dir) as reader:
            assert reader.has_body(SUN)
        # After exit, reader should be closed
        # Accessing may or may not raise -- just verify no crash on exit

    @pytest.mark.unit
    def test_explicit_close(self, single_tier_dir):
        """Explicit close() doesn't crash."""
        reader = CompositeLEBReader.from_directory(single_tier_dir)
        reader.eval_body(SUN, JD_J2000)
        reader.close()

    @SKIP_NO_LEB2_CORE
    @pytest.mark.unit
    def test_manual_construction(self):
        """Constructing CompositeLEBReader from list of readers."""
        r1 = open_leb(LEB2_BASE_CORE)
        readers = [r1]
        if os.path.exists(LEB2_BASE_ASTEROIDS):
            readers.append(open_leb(LEB2_BASE_ASTEROIDS))

        comp = CompositeLEBReader(readers)
        assert comp.has_body(SUN)
        pos, vel = comp.eval_body(SUN, JD_J2000)
        # eval_body returns raw cartesian AU — just verify finite values
        for v in pos + vel:
            assert math.isfinite(v)
        comp.close()
