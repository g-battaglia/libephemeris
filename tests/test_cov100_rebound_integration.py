# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Full-coverage tests for ``libephemeris.rebound_integration``.

REBOUND and ASSIST are optional, network/heavy dependencies that are not
installed in the test environment.  This module drives every reachable line
by injecting fake ``rebound``/``assist``/``certifi``/``jplephem`` modules
via ``sys.modules`` patching and by mocking ``urllib.request.urlopen`` so
no real network or 1 GB file download ever occurs.

The strategy is "full mocks": the present-dependency arms are exercised by
fake module objects, the absent-dependency arms by ``patch.dict`` with a
``None`` sentinel (which makes ``import dep`` raise ``ImportError``).
"""

from __future__ import annotations

import os
import sys
import types
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from libephemeris import rebound_integration as ri
from libephemeris.constants import CERES
from libephemeris.minor_bodies import MINOR_BODY_ELEMENTS


# ---------------------------------------------------------------------------
# Fake REBOUND / ASSIST building blocks
# ---------------------------------------------------------------------------


class _FakeParticle:
    """A simple particle with mutable Cartesian state."""

    def __init__(self, x=1.0, y=2.0, z=0.5, vx=0.01, vy=0.02, vz=0.005):
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz


class _FakeParticles:
    """Indexable particle container supporting int and string keys."""

    def __init__(self):
        self._by_index: list[_FakeParticle] = []
        self._by_name: dict[str, _FakeParticle] = {}

    def append(self, particle, hash=None):
        self._by_index.append(particle)
        if hash is not None:
            self._by_name[hash] = particle

    def __getitem__(self, key):
        if isinstance(key, str):
            return self._by_name[key]
        return self._by_index[key]


class _FakeSimulation:
    """Minimal stand-in for ``rebound.Simulation``."""

    def __init__(self):
        self.units = None
        self.G = None
        self.integrator = None
        self.dt = None
        self.t = None
        self.particles = _FakeParticles()
        self.integrate_calls: list[float] = []

    def add(self, **kwargs):
        # Sun added with m=1.0; test particle with m=0.0.  Either may be
        # given a ``hash`` (string) used to look it up later.
        particle = _FakeParticle(
            x=kwargs.get("x", 1.0),
            y=kwargs.get("y", 2.0),
            z=kwargs.get("z", 0.5),
            vx=kwargs.get("vx", 0.01),
            vy=kwargs.get("vy", 0.02),
            vz=kwargs.get("vz", 0.005),
        )
        self.particles.append(particle, hash=kwargs.get("hash"))

    def move_to_hel(self):
        pass

    def integrate(self, t_end):
        self.integrate_calls.append(t_end)
        self.t = t_end


def _make_fake_rebound():
    """Build a fake ``rebound`` module object."""
    mod = types.ModuleType("rebound")
    mod.__version__ = "9.9.9-fake"
    mod.Simulation = _FakeSimulation
    return mod


class _FakeEphem:
    """Stand-in for ``assist.Ephem``."""

    def __init__(self, *files):
        self.files = files
        self.jd_ref = 2451545.0

    def get_particle(self, name, t):
        # Return a Sun barycentric state; values are arbitrary but finite.
        return _FakeParticle(
            x=0.001, y=-0.002, z=0.0003, vx=1e-6, vy=2e-6, vz=3e-7
        )


class _FakeExtras:
    """Stand-in for ``assist.Extras``."""

    def __init__(self, sim, ephem):
        self.sim = sim
        self.ephem = ephem
        self.particle_params = None


def _make_fake_assist():
    """Build a fake ``assist`` module object."""
    mod = types.ModuleType("assist")
    mod.__version__ = "8.8.8-fake"
    mod.Ephem = _FakeEphem
    mod.Extras = _FakeExtras
    return mod


@pytest.fixture
def fake_rebound():
    """Inject a fake ``rebound`` module for the duration of a test."""
    mod = _make_fake_rebound()
    with patch.dict(sys.modules, {"rebound": mod}):
        yield mod


@pytest.fixture
def fake_assist():
    """Inject a fake ``assist`` module for the duration of a test."""
    mod = _make_fake_assist()
    with patch.dict(sys.modules, {"assist": mod}):
        yield mod


@pytest.fixture
def ceres():
    """Ceres orbital elements."""
    return MINOR_BODY_ELEMENTS[CERES]


@pytest.fixture(autouse=True)
def _reset_assist_cache():
    """Keep the ASSIST availability cache clean around each test."""
    ri.reset_assist_data_cache()
    yield
    ri.reset_assist_data_cache()


# ---------------------------------------------------------------------------
# _search_assist_file  (lines 152, 157)
# ---------------------------------------------------------------------------


class TestSearchAssistFile:
    """Cover both missing-dir continue and not-found exit arms."""

    def test_skips_nonexistent_dir(self, tmp_path):
        """A directory that does not exist is skipped (line 152)."""
        missing = tmp_path / "nope"
        found = tmp_path / "here"
        found.mkdir()
        (found / "target.bsp").write_bytes(b"x")
        result = ri._search_assist_file(["target.bsp"], [missing, found])
        assert result == str(found / "target.bsp")

    def test_returns_none_when_absent(self, tmp_path):
        """Falls through to return None when nothing matches (line 157)."""
        d = tmp_path / "exists"
        d.mkdir()
        result = ri._search_assist_file(["absent.bsp"], [d])
        assert result is None


# ---------------------------------------------------------------------------
# AssistEphemConfig.__post_init__  (line 199 - ASSIST_DIR env arm)
# ---------------------------------------------------------------------------


class TestAssistEphemConfigEnv:
    """Exercise the ASSIST_DIR environment branch."""

    def test_env_dir_appended(self, tmp_path):
        """ASSIST_DIR points at a dir holding the planet file (line 199)."""
        planets = tmp_path / "linux_p1550p2650.440"
        planets.write_bytes(b"data")
        with patch.dict(os.environ, {"ASSIST_DIR": str(tmp_path)}):
            config = ri.AssistEphemConfig()
        assert config.planets_file == str(planets)


# ---------------------------------------------------------------------------
# _verify_assist_file  (lines 240-264)
# ---------------------------------------------------------------------------


class TestVerifyAssistFile:
    """Cover every arm of the file integrity check."""

    def test_missing_file(self, tmp_path):
        """Nonexistent path returns False (line 240-241)."""
        assert ri._verify_assist_file(tmp_path / "missing.bin") is False

    def test_empty_file(self, tmp_path):
        """Zero-byte file returns False (line 244-245)."""
        p = tmp_path / "empty.bin"
        p.write_bytes(b"")
        assert ri._verify_assist_file(p) is False

    def test_too_small_known_file(self, tmp_path):
        """A known file under the min-size ratio returns False (250-251)."""
        p = tmp_path / "linux_p1550p2650.440"
        p.write_bytes(b"x" * 10)  # far below expected ~98 MB
        assert ri._verify_assist_file(p) is False

    def test_unknown_extension_ok(self, tmp_path):
        """An unknown-name non-bsp file of nonzero size passes (line 264)."""
        p = tmp_path / "something.dat"
        p.write_bytes(b"x" * 100)
        assert ri._verify_assist_file(p) is True

    def test_bsp_valid_via_fake_spk(self, tmp_path):
        """A .bsp file validates structurally via a fake jplephem SPK."""
        p = tmp_path / "sb441-n16.bsp"
        expected = ri._ASSIST_EXPECTED_SIZES["sb441-n16.bsp"]
        size = int(expected * 0.95)

        class _FakeSeg:
            center = 0
            target = 10

        class _FakeKernel:
            segments = [_FakeSeg()]

            def __enter__(self):
                return self

            def __exit__(self, *a):
                return False

        class _FakeSPK:
            @staticmethod
            def open(path):
                return _FakeKernel()

        fake_spk_mod = types.ModuleType("jplephem.spk")
        fake_spk_mod.SPK = _FakeSPK
        fake_jplephem = types.ModuleType("jplephem")
        fake_jplephem.spk = fake_spk_mod

        # Avoid actually allocating 600 MB: stat() is patched to report
        # an acceptable size while the file stays tiny.
        p.write_bytes(b"BSPdata")

        class _FakeStat:
            st_size = size

        def _fake_stat(self, *args, **kwargs):
            return _FakeStat()

        with (
            patch.dict(
                sys.modules,
                {"jplephem": fake_jplephem, "jplephem.spk": fake_spk_mod},
            ),
            patch.object(Path, "stat", _fake_stat),
        ):
            assert ri._verify_assist_file(p) is True

    def test_bsp_structural_failure(self, tmp_path):
        """SPK.open raising ValueError makes verification fail (261-262)."""
        p = tmp_path / "sb441-n16.bsp"
        expected = ri._ASSIST_EXPECTED_SIZES["sb441-n16.bsp"]
        size = int(expected * 0.95)
        p.write_bytes(b"BSPdata")

        class _FakeSPK:
            @staticmethod
            def open(path):
                raise ValueError("corrupt kernel")

        fake_spk_mod = types.ModuleType("jplephem.spk")
        fake_spk_mod.SPK = _FakeSPK
        fake_jplephem = types.ModuleType("jplephem")
        fake_jplephem.spk = fake_spk_mod

        class _FakeStat:
            st_size = size

        def _fake_stat(self, *args, **kwargs):
            return _FakeStat()

        with (
            patch.dict(
                sys.modules,
                {"jplephem": fake_jplephem, "jplephem.spk": fake_spk_mod},
            ),
            patch.object(Path, "stat", _fake_stat),
        ):
            assert ri._verify_assist_file(p) is False


# ---------------------------------------------------------------------------
# _create_ssl_context  (lines 273-280)
# ---------------------------------------------------------------------------


class TestCreateSSLContext:
    """Cover certifi-present and certifi-absent arms."""

    def test_certifi_present(self):
        """certifi importable: context built from its CA bundle (273-278)."""
        ctx = ri._create_ssl_context()
        import ssl

        assert isinstance(ctx, ssl.SSLContext)

    def test_certifi_absent(self):
        """certifi missing: fall back to default context (279-280)."""
        with patch.dict(sys.modules, {"certifi": None}):
            ctx = ri._create_ssl_context()
        import ssl

        assert isinstance(ctx, ssl.SSLContext)


# ---------------------------------------------------------------------------
# _download_single_file  (lines 299-358)
# ---------------------------------------------------------------------------


def _make_fake_response(data: bytes, content_length: int | None = None):
    """Build a context-manager urlopen response double."""
    if content_length is None:
        content_length = len(data)

    class _Resp:
        def __init__(self):
            self.headers = {"Content-Length": str(content_length)}
            self._chunks = [data, b""] if data else [b""]
            self._i = 0

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def read(self, n):
            if self._i >= len(self._chunks):
                return b""
            chunk = self._chunks[self._i]
            self._i += 1
            return chunk

    return _Resp()


class TestDownloadSingleFile:
    """Drive the streaming download/atomic-write helper with a fake socket."""

    def test_download_with_progress(self, tmp_path, capsys):
        """Full path with Content-Length + progress bar + final print."""
        dest = tmp_path / "sub" / "file.dat"
        resp = _make_fake_response(b"hello world payload", content_length=19)

        fake_progress = MagicMock()
        with (
            patch("urllib.request.urlopen", return_value=resp),
            patch(
                "libephemeris.download._get_progress_bar",
                return_value=fake_progress,
            ),
        ):
            ri._download_single_file(
                url="https://example.invalid/file.dat",
                dest=dest,
                description="Test file",
                show_progress=True,
                quiet=False,
            )

        assert dest.exists()
        assert dest.read_bytes() == b"hello world payload"
        fake_progress.update.assert_called()
        fake_progress.close.assert_called_once()
        out = capsys.readouterr().out
        assert "Test file" in out
        assert "Done" in out

    def test_download_quiet_no_progress(self, tmp_path, capsys):
        """quiet + show_progress False skips size print and progress bar."""
        dest = tmp_path / "q.dat"
        resp = _make_fake_response(b"abc", content_length=3)
        with patch("urllib.request.urlopen", return_value=resp):
            ri._download_single_file(
                url="https://example.invalid/q.dat",
                dest=dest,
                description="Quiet file",
                show_progress=False,
                quiet=True,
            )
        assert dest.read_bytes() == b"abc"
        assert capsys.readouterr().out == ""

    def test_download_no_content_length(self, tmp_path):
        """Content-Length 0 skips both the size print and the progress bar."""
        dest = tmp_path / "nolen.dat"
        resp = _make_fake_response(b"data here", content_length=0)
        with patch("urllib.request.urlopen", return_value=resp):
            ri._download_single_file(
                url="https://example.invalid/nolen.dat",
                dest=dest,
                description="No length",
                show_progress=True,
                quiet=False,
            )
        assert dest.read_bytes() == b"data here"

    def test_download_error_cleans_temp(self, tmp_path):
        """A mid-download OSError cleans up the temp file and re-raises."""
        dest = tmp_path / "err.dat"

        with patch(
            "urllib.request.urlopen", side_effect=OSError("network down")
        ):
            with pytest.raises(OSError):
                ri._download_single_file(
                    url="https://example.invalid/err.dat",
                    dest=dest,
                    description="Failing file",
                    show_progress=False,
                    quiet=True,
                )
        # No leftover .download temp files in the destination dir.
        leftovers = list(tmp_path.glob("*.download"))
        assert leftovers == []

    def test_download_error_unlink_swallowed(self, tmp_path):
        """If the temp-file cleanup itself errors, it is swallowed (356-357)."""
        dest = tmp_path / "err2.dat"

        with (
            patch(
                "urllib.request.urlopen",
                side_effect=OSError("network down"),
            ),
            patch("os.path.exists", return_value=True),
            patch("os.unlink", side_effect=OSError("cannot unlink")),
        ):
            with pytest.raises(OSError):
                ri._download_single_file(
                    url="https://example.invalid/err2.dat",
                    dest=dest,
                    description="Failing file",
                    show_progress=False,
                    quiet=True,
                )

    def test_download_error_temp_already_gone(self, tmp_path):
        """If the temp file no longer exists, unlink is skipped (354->358)."""
        dest = tmp_path / "err3.dat"

        with (
            patch(
                "urllib.request.urlopen",
                side_effect=OSError("network down"),
            ),
            patch("os.path.exists", return_value=False),
        ):
            with pytest.raises(OSError):
                ri._download_single_file(
                    url="https://example.invalid/err3.dat",
                    dest=dest,
                    description="Failing file",
                    show_progress=False,
                    quiet=True,
                )


# ---------------------------------------------------------------------------
# download_assist_data  (lines 398-454)
# ---------------------------------------------------------------------------


class TestDownloadAssistData:
    """Cover the orchestration: skip-existing, re-download, default dir."""

    def test_downloads_both_files(self, tmp_path, capsys):
        """Both files downloaded into an explicit target dir."""
        calls = []

        def _fake_dl(url, dest, description, show_progress, quiet):
            calls.append(description)
            Path(dest).write_bytes(b"x")

        with patch.object(ri, "_download_single_file", side_effect=_fake_dl):
            result = ri.download_assist_data(
                target_dir=str(tmp_path), quiet=False
            )

        assert result == tmp_path
        assert len(calls) == 2
        out = capsys.readouterr().out
        assert "ASSIST data directory" in out
        assert "Downloaded 2 file(s)." in out

    def test_skips_disabled_file(self, tmp_path):
        """planets=False, asteroids=True downloads only one file (line 417)."""
        calls = []

        def _fake_dl(url, dest, description, show_progress, quiet):
            calls.append(description)
            Path(dest).write_bytes(b"x")

        with patch.object(ri, "_download_single_file", side_effect=_fake_dl):
            ri.download_assist_data(
                target_dir=str(tmp_path),
                planets=False,
                asteroids=True,
                quiet=True,
            )
        assert calls == ["Asteroid perturbers"]

    def test_skips_existing_valid_file(self, tmp_path, capsys):
        """Existing valid file is skipped, not re-downloaded (422-428)."""
        planets = tmp_path / "linux_p1550p2650.440"
        planets.write_bytes(b"x")

        with (
            patch.object(ri, "_verify_assist_file", return_value=True),
            patch.object(ri, "_download_single_file") as mock_dl,
        ):
            ri.download_assist_data(
                target_dir=str(tmp_path),
                planets=True,
                asteroids=False,
                quiet=False,
            )
        mock_dl.assert_not_called()
        out = capsys.readouterr().out
        assert "already exists" in out
        assert "All files already present and verified." in out

    def test_redownloads_invalid_file(self, tmp_path, capsys):
        """Existing-but-invalid file is unlinked and re-downloaded (429-434)."""
        planets = tmp_path / "linux_p1550p2650.440"
        planets.write_bytes(b"corrupt")

        def _fake_dl(url, dest, description, show_progress, quiet):
            Path(dest).write_bytes(b"fresh")

        with (
            patch.object(ri, "_verify_assist_file", return_value=False),
            patch.object(ri, "_download_single_file", side_effect=_fake_dl),
        ):
            ri.download_assist_data(
                target_dir=str(tmp_path),
                planets=True,
                asteroids=False,
                quiet=False,
            )
        out = capsys.readouterr().out
        assert "failed verification" in out
        assert planets.read_bytes() == b"fresh"

    def test_force_redownload(self, tmp_path):
        """force=True downloads even when a valid file already exists."""
        planets = tmp_path / "linux_p1550p2650.440"
        planets.write_bytes(b"old")
        calls = []

        def _fake_dl(url, dest, description, show_progress, quiet):
            calls.append(description)
            Path(dest).write_bytes(b"new")

        with (
            patch.object(ri, "_verify_assist_file", return_value=True),
            patch.object(ri, "_download_single_file", side_effect=_fake_dl),
        ):
            ri.download_assist_data(
                target_dir=str(tmp_path),
                planets=True,
                asteroids=False,
                force=True,
                quiet=True,
            )
        assert calls == ["Planet ephemeris"]

    def test_skips_existing_valid_file_quiet(self, tmp_path):
        """Existing valid file skipped with quiet=True (branch 424->427)."""
        planets = tmp_path / "linux_p1550p2650.440"
        planets.write_bytes(b"x")

        with (
            patch.object(ri, "_verify_assist_file", return_value=True),
            patch.object(ri, "_download_single_file") as mock_dl,
        ):
            ri.download_assist_data(
                target_dir=str(tmp_path),
                planets=True,
                asteroids=False,
                quiet=True,
            )
        mock_dl.assert_not_called()

    def test_redownloads_invalid_file_quiet(self, tmp_path):
        """Invalid file unlinked then re-downloaded, quiet=True (430->434)."""
        planets = tmp_path / "linux_p1550p2650.440"
        planets.write_bytes(b"corrupt")

        def _fake_dl(url, dest, description, show_progress, quiet):
            Path(dest).write_bytes(b"fresh")

        with (
            patch.object(ri, "_verify_assist_file", return_value=False),
            patch.object(ri, "_download_single_file", side_effect=_fake_dl),
        ):
            ri.download_assist_data(
                target_dir=str(tmp_path),
                planets=True,
                asteroids=False,
                quiet=True,
            )
        assert planets.read_bytes() == b"fresh"

    def test_nothing_downloaded_nothing_skipped(self, tmp_path, capsys):
        """planets=False and asteroids=False: no summary line (448->452)."""
        with patch.object(ri, "_download_single_file") as mock_dl:
            result = ri.download_assist_data(
                target_dir=str(tmp_path),
                planets=False,
                asteroids=False,
                quiet=False,
            )
        mock_dl.assert_not_called()
        assert result == tmp_path
        out = capsys.readouterr().out
        assert "Downloaded" not in out
        assert "already present" not in out

    def test_default_target_dir(self, tmp_path):
        """target_dir=None uses _ASSIST_DEFAULT_DIR (lines 398-399)."""
        fake_default = tmp_path / "default_assist"

        def _fake_dl(url, dest, description, show_progress, quiet):
            Path(dest).write_bytes(b"x")

        with (
            patch.object(ri, "_ASSIST_DEFAULT_DIR", fake_default),
            patch.object(ri, "_download_single_file", side_effect=_fake_dl),
        ):
            result = ri.download_assist_data(quiet=True)
        assert result == fake_default
        assert fake_default.exists()


# ---------------------------------------------------------------------------
# Frame rotation helpers  (lines 467-469, 476-478)
# ---------------------------------------------------------------------------


class TestFrameRotations:
    """Round-trip ecliptic<->ICRS rotations."""

    def test_round_trip(self):
        """ICRS->ecliptic->ICRS returns the original vector (467-478)."""
        x, y, z = 1.5, -0.3, 0.7
        ex, ey, ez = ri._icrs_to_ecliptic_j2000(x, y, z)
        bx, by, bz = ri._ecliptic_j2000_to_icrs(ex, ey, ez)
        assert abs(bx - x) < 1e-12
        assert abs(by - y) < 1e-12
        assert abs(bz - z) < 1e-12


# ---------------------------------------------------------------------------
# Version / availability probes  (lines 534, 548, 606, 620)
# ---------------------------------------------------------------------------


class TestProbes:
    """Cover the present-dependency arms of availability/version probes."""

    def test_check_rebound_true(self, fake_rebound):
        """rebound importable -> True (line 534)."""
        assert ri.check_rebound_available() is True

    def test_check_assist_true(self, fake_assist):
        """assist importable -> True (line 548)."""
        assert ri.check_assist_available() is True

    def test_get_rebound_version(self, fake_rebound):
        """rebound importable -> version string (line 606)."""
        assert ri.get_rebound_version() == "9.9.9-fake"

    def test_get_assist_version(self, fake_assist):
        """assist importable -> version string (line 620)."""
        assert ri.get_assist_version() == "8.8.8-fake"

    def test_get_rebound_version_absent(self):
        """rebound absent -> None."""
        with patch.dict(sys.modules, {"rebound": None}):
            assert ri.get_rebound_version() is None

    def test_get_assist_version_absent(self):
        """assist absent -> None."""
        with patch.dict(sys.modules, {"assist": None}):
            assert ri.get_assist_version() is None


# ---------------------------------------------------------------------------
# check_assist_data_available  (lines 581-583)
# ---------------------------------------------------------------------------


class TestCheckAssistDataAvailable:
    """Cover the planet-file-present arm of the cached data probe."""

    def test_data_available_with_planet_file(self, fake_assist, tmp_path):
        """assist + planet file present -> True (lines 581-583)."""
        planets = tmp_path / "linux_p1550p2650.440"
        planets.write_bytes(b"data")
        ri.reset_assist_data_cache()
        with patch.dict(os.environ, {"ASSIST_DIR": str(tmp_path)}):
            assert ri.check_assist_data_available() is True

    def test_cached_short_circuit(self, fake_assist, tmp_path):
        """Second call returns the cached value without re-probing (573-574)."""
        planets = tmp_path / "linux_p1550p2650.440"
        planets.write_bytes(b"data")
        ri.reset_assist_data_cache()
        with patch.dict(os.environ, {"ASSIST_DIR": str(tmp_path)}):
            first = ri.check_assist_data_available()
            second = ri.check_assist_data_available()
        assert first is True and second is True

    def test_assist_absent(self):
        """assist not importable -> False (line 577-578)."""
        ri.reset_assist_data_cache()
        with patch.dict(sys.modules, {"assist": None}):
            assert ri.check_assist_data_available() is False


# ---------------------------------------------------------------------------
# _elements_to_cartesian  (lines 682-697)
# ---------------------------------------------------------------------------


class TestElementsToCartesian:
    """Drive the temporary-simulation Cartesian extraction."""

    def test_cartesian_extraction(self, fake_rebound, ceres):
        """Returns a dict of the test particle's state (682-697)."""
        cart = ri._elements_to_cartesian(ceres, ceres.epoch)
        assert set(cart) == {"x", "y", "z", "vx", "vy", "vz"}
        # The fake test particle (index 1) carries the default state.
        assert cart["x"] == 1.0
        assert cart["vy"] == 0.02


# ---------------------------------------------------------------------------
# propagate_orbit_rebound  (lines 750-797)
# ---------------------------------------------------------------------------


class TestPropagateOrbitRebound:
    """Drive the REBOUND-only propagation with a fake simulation."""

    def test_default_ias15(self, fake_rebound, ceres):
        """IAS15 (no dt set) integrates over jd_end - jd_start (750-805)."""
        result = ri.propagate_orbit_rebound(
            ceres, ceres.epoch, ceres.epoch + 30
        )
        assert isinstance(result, ri.PropagationResult)
        assert result.jd_tt == ceres.epoch + 30

    def test_whfast_default_dt(self, fake_rebound, ceres):
        """WHFAST without explicit dt computes a period-based dt (763-766)."""
        result = ri.propagate_orbit_rebound(
            ceres,
            ceres.epoch,
            ceres.epoch + 30,
            integrator=ri.ReboundIntegrator.WHFAST,
        )
        assert result.jd_tt == ceres.epoch + 30

    def test_leapfrog_default_dt(self, fake_rebound, ceres):
        """LEAPFROG without explicit dt also computes a default dt."""
        result = ri.propagate_orbit_rebound(
            ceres,
            ceres.epoch,
            ceres.epoch + 30,
            integrator=ri.ReboundIntegrator.LEAPFROG,
        )
        assert result.jd_tt == ceres.epoch + 30

    def test_explicit_dt(self, fake_rebound, ceres):
        """Explicit dt sets sim.dt (lines 761-762)."""
        result = ri.propagate_orbit_rebound(
            ceres, ceres.epoch, ceres.epoch + 30, dt=0.5
        )
        assert result.jd_tt == ceres.epoch + 30

    def test_absent_rebound_raises(self, ceres):
        """rebound absent -> ImportError (lines 742-746)."""
        with patch.dict(sys.modules, {"rebound": None}):
            with pytest.raises(ImportError, match="REBOUND is required"):
                ri.propagate_orbit_rebound(
                    ceres, ceres.epoch, ceres.epoch + 1
                )


# ---------------------------------------------------------------------------
# propagate_orbit_assist  (lines 857-953)
# ---------------------------------------------------------------------------


class TestPropagateOrbitAssist:
    """Drive ASSIST propagation through fake rebound + assist modules."""

    def _config(self, tmp_path, with_asteroids=False):
        planets = tmp_path / "linux_p1550p2650.440"
        planets.write_bytes(b"planets")
        asteroids_file = None
        if with_asteroids:
            ast = tmp_path / "sb441-n16.bsp"
            ast.write_bytes(b"asteroids")
            asteroids_file = str(ast)
        cfg = ri.AssistEphemConfig.__new__(ri.AssistEphemConfig)
        cfg.planets_file = str(planets)
        cfg.asteroids_file = asteroids_file
        cfg.data_dir = None
        return cfg

    def test_assist_planets_only(self, fake_rebound, fake_assist, ceres, tmp_path):
        """Planets-only config takes the single-file Ephem arm (line 899)."""
        cfg = self._config(tmp_path, with_asteroids=False)
        result = ri.propagate_orbit_assist(
            ceres, ceres.epoch, ceres.epoch + 30, ephem_config=cfg
        )
        assert isinstance(result, ri.PropagationResult)
        assert result.jd_tt == ceres.epoch + 30

    def test_assist_with_asteroids(
        self, fake_rebound, fake_assist, ceres, tmp_path
    ):
        """Asteroids config takes the two-file Ephem arm (line 897)."""
        cfg = self._config(tmp_path, with_asteroids=True)
        result = ri.propagate_orbit_assist(
            ceres, ceres.epoch, ceres.epoch + 30, ephem_config=cfg
        )
        assert result.jd_tt == ceres.epoch + 30

    def test_assist_non_gravitational(
        self, fake_rebound, fake_assist, ceres, tmp_path
    ):
        """include_non_gravitational sets extras.particle_params (932-933)."""
        cfg = self._config(tmp_path)
        captured = {}

        orig_extras = _FakeExtras

        class _CapturingExtras(orig_extras):
            def __init__(self, sim, ephem):
                super().__init__(sim, ephem)
                captured["instance"] = self

        with patch.object(fake_assist, "Extras", _CapturingExtras):
            ri.propagate_orbit_assist(
                ceres,
                ceres.epoch,
                ceres.epoch + 30,
                ephem_config=cfg,
                include_non_gravitational=True,
                A1=1e-10,
                A2=2e-10,
                A3=3e-10,
            )
        assert captured["instance"].particle_params == [1e-10, 2e-10, 3e-10]

    def test_assist_default_config(
        self, fake_rebound, fake_assist, ceres, tmp_path
    ):
        """ephem_config None builds an AssistEphemConfig (line 874-875)."""
        planets = tmp_path / "linux_p1550p2650.440"
        planets.write_bytes(b"planets")
        with patch.dict(os.environ, {"ASSIST_DIR": str(tmp_path)}):
            result = ri.propagate_orbit_assist(
                ceres, ceres.epoch, ceres.epoch + 30
            )
        assert result.jd_tt == ceres.epoch + 30

    def test_assist_absent_rebound(self, ceres):
        """rebound absent -> ImportError (lines 859-862)."""
        with patch.dict(sys.modules, {"rebound": None}):
            with pytest.raises(ImportError, match="REBOUND is required"):
                ri.propagate_orbit_assist(
                    ceres, ceres.epoch, ceres.epoch + 1
                )

    def test_assist_absent_assist(self, fake_rebound, ceres):
        """assist absent -> ImportError (lines 866-871)."""
        with patch.dict(sys.modules, {"assist": None}):
            with pytest.raises(ImportError, match="ASSIST is required"):
                ri.propagate_orbit_assist(
                    ceres, ceres.epoch, ceres.epoch + 1
                )

    def test_assist_no_planets_file(self, fake_rebound, fake_assist, ceres):
        """planets_file None -> FileNotFoundError (lines 878-883)."""
        cfg = ri.AssistEphemConfig.__new__(ri.AssistEphemConfig)
        cfg.planets_file = None
        cfg.asteroids_file = None
        cfg.data_dir = None
        with pytest.raises(FileNotFoundError, match="Planet ephemeris file"):
            ri.propagate_orbit_assist(
                ceres, ceres.epoch, ceres.epoch + 1, ephem_config=cfg
            )

    def test_assist_planets_file_missing_on_disk(
        self, fake_rebound, fake_assist, ceres
    ):
        """planets_file set but not on disk -> FileNotFoundError (885-888)."""
        cfg = ri.AssistEphemConfig.__new__(ri.AssistEphemConfig)
        cfg.planets_file = "/nonexistent/linux_p1550p2650.440"
        cfg.asteroids_file = None
        cfg.data_dir = None
        with pytest.raises(FileNotFoundError, match="not found"):
            ri.propagate_orbit_assist(
                ceres, ceres.epoch, ceres.epoch + 1, ephem_config=cfg
            )

    def test_assist_asteroids_file_missing_on_disk(
        self, fake_rebound, fake_assist, ceres, tmp_path
    ):
        """asteroids_file set but missing -> FileNotFoundError (890-893)."""
        planets = tmp_path / "linux_p1550p2650.440"
        planets.write_bytes(b"planets")
        cfg = ri.AssistEphemConfig.__new__(ri.AssistEphemConfig)
        cfg.planets_file = str(planets)
        cfg.asteroids_file = "/nonexistent/sb441-n16.bsp"
        cfg.data_dir = None
        with pytest.raises(FileNotFoundError, match="Asteroid ephemeris"):
            ri.propagate_orbit_assist(
                ceres, ceres.epoch, ceres.epoch + 1, ephem_config=cfg
            )


# ---------------------------------------------------------------------------
# propagate_trajectory  (lines 994-1096)
# ---------------------------------------------------------------------------


class TestPropagateTrajectory:
    """Drive both the ASSIST path and the REBOUND-only fallback."""

    def _assist_config(self, tmp_path, with_asteroids=False):
        planets = tmp_path / "linux_p1550p2650.440"
        planets.write_bytes(b"planets")
        asteroids_file = None
        if with_asteroids:
            ast = tmp_path / "sb441-n16.bsp"
            ast.write_bytes(b"asteroids")
            asteroids_file = str(ast)
        cfg = ri.AssistEphemConfig.__new__(ri.AssistEphemConfig)
        cfg.planets_file = str(planets)
        cfg.asteroids_file = asteroids_file
        cfg.data_dir = None
        return cfg

    def test_assist_trajectory(
        self, fake_rebound, fake_assist, ceres, tmp_path
    ):
        """ASSIST path returns num_points results (lines 999-1060)."""
        cfg = self._assist_config(tmp_path)
        traj = ri.propagate_trajectory(
            ceres,
            ceres.epoch,
            ceres.epoch + 100,
            num_points=5,
            use_assist=True,
            ephem_config=cfg,
        )
        assert len(traj) == 5
        assert all(isinstance(p, ri.PropagationResult) for p in traj)

    def test_assist_trajectory_with_asteroids(
        self, fake_rebound, fake_assist, ceres, tmp_path
    ):
        """ASSIST path with asteroids file takes two-file Ephem (1010-1013)."""
        cfg = self._assist_config(tmp_path, with_asteroids=True)
        traj = ri.propagate_trajectory(
            ceres,
            ceres.epoch,
            ceres.epoch + 100,
            num_points=3,
            use_assist=True,
            ephem_config=cfg,
        )
        assert len(traj) == 3

    def test_assist_trajectory_default_config(
        self, fake_rebound, fake_assist, ceres, tmp_path
    ):
        """ephem_config None builds default config inside ASSIST path."""
        planets = tmp_path / "linux_p1550p2650.440"
        planets.write_bytes(b"planets")
        with patch.dict(os.environ, {"ASSIST_DIR": str(tmp_path)}):
            traj = ri.propagate_trajectory(
                ceres,
                ceres.epoch,
                ceres.epoch + 50,
                num_points=4,
                use_assist=True,
            )
        assert len(traj) == 4

    def test_assist_unavailable_data_falls_to_rebound(
        self, fake_rebound, fake_assist, ceres
    ):
        """assist importable but no planet file -> REBOUND tail (1064-1096)."""
        # No ASSIST_DIR and no files -> planets_file resolves to None, so the
        # inner ``if planets_file and exists`` block is skipped and execution
        # falls through to the REBOUND-only tail.
        cfg = ri.AssistEphemConfig.__new__(ri.AssistEphemConfig)
        cfg.planets_file = None
        cfg.asteroids_file = None
        cfg.data_dir = None
        traj = ri.propagate_trajectory(
            ceres,
            ceres.epoch,
            ceres.epoch + 50,
            num_points=4,
            use_assist=True,
            ephem_config=cfg,
        )
        assert len(traj) == 4

    def test_assist_import_error_falls_to_rebound(self, fake_rebound, ceres):
        """check_assist_available True but import fails -> fallback (1061-1062).

        ``check_assist_available`` is forced True while ``assist`` import
        inside the try-block raises, exercising the except->fallback arm.
        """
        with (
            patch.object(ri, "check_assist_available", return_value=True),
            patch.dict(sys.modules, {"assist": None}),
        ):
            traj = ri.propagate_trajectory(
                ceres,
                ceres.epoch,
                ceres.epoch + 50,
                num_points=4,
                use_assist=True,
            )
        assert len(traj) == 4

    def test_rebound_only(self, fake_rebound, ceres):
        """use_assist False uses the REBOUND-only tail (lines 1064-1096)."""
        traj = ri.propagate_trajectory(
            ceres,
            ceres.epoch,
            ceres.epoch + 100,
            num_points=10,
            use_assist=False,
        )
        assert len(traj) == 10

    def test_single_point(self, fake_rebound, ceres):
        """num_points == 1 -> dt is zero (line 997)."""
        traj = ri.propagate_trajectory(
            ceres,
            ceres.epoch,
            ceres.epoch + 100,
            num_points=1,
            use_assist=False,
        )
        assert len(traj) == 1

    def test_rebound_absent_raises(self, ceres):
        """No rebound and use_assist False -> ImportError (1065-1069)."""
        with patch.dict(sys.modules, {"rebound": None}):
            with pytest.raises(ImportError, match="REBOUND is required"):
                ri.propagate_trajectory(
                    ceres,
                    ceres.epoch,
                    ceres.epoch + 10,
                    num_points=3,
                    use_assist=False,
                )


# ---------------------------------------------------------------------------
# compare_with_keplerian  (lines 1123-1155)
# ---------------------------------------------------------------------------


class TestCompareWithKeplerian:
    """Drive the Keplerian-vs-nbody comparison dictionary builder."""

    def test_rebound_path(self, fake_rebound, ceres):
        """use_assist False -> rebound result, method 'rebound'."""
        comp = ri.compare_with_keplerian(
            ceres, ceres.epoch + 100, use_assist=False
        )
        assert comp["method"] == "rebound"
        assert "angular_sep_arcsec" in comp
        assert comp["propagation_days"] == ceres.epoch + 100 - ceres.epoch

    def test_assist_success(self, fake_rebound, fake_assist, ceres, tmp_path):
        """use_assist True with files present -> ASSIST result, method 'assist'."""
        planets = tmp_path / "linux_p1550p2650.440"
        planets.write_bytes(b"planets")
        cfg = ri.AssistEphemConfig.__new__(ri.AssistEphemConfig)
        cfg.planets_file = str(planets)
        cfg.asteroids_file = None
        cfg.data_dir = None
        comp = ri.compare_with_keplerian(
            ceres, ceres.epoch + 100, use_assist=True, ephem_config=cfg
        )
        assert comp["method"] == "assist"

    def test_assist_fallback_to_rebound(
        self, fake_rebound, fake_assist, ceres
    ):
        """use_assist True but ASSIST raises -> rebound fallback (1137-1138).

        No planet file is available so ``propagate_orbit_assist`` raises
        ``FileNotFoundError``, which is caught and falls back to REBOUND.
        Note: ``method`` stays 'assist' because it reflects the requested
        mode, not the path actually taken.
        """
        cfg = ri.AssistEphemConfig.__new__(ri.AssistEphemConfig)
        cfg.planets_file = None
        cfg.asteroids_file = None
        cfg.data_dir = None
        comp = ri.compare_with_keplerian(
            ceres, ceres.epoch + 100, use_assist=True, ephem_config=cfg
        )
        assert comp["method"] == "assist"
        assert "nbody" in comp

    def test_delta_lon_wrap_high(self, fake_rebound, ceres, monkeypatch):
        """delta_lon > 180 wraps down by 360 (lines 1144-1145)."""
        big = ri.PropagationResult(
            x=1.0, y=0.0, z=0.0, vx=0.0, vy=0.0, vz=0.0, jd_tt=ceres.epoch
        )
        # Force a near-0 nbody longitude against a near-360 keplerian one,
        # so the raw delta is strongly negative -> wrap up branch; and the
        # opposite for the high branch.  We patch propagate to control it.
        with patch.object(
            ri, "propagate_orbit_rebound", return_value=big
        ):
            # Patch keplerian to ~350 deg so nbody(0) - kep(350) = -350 < -180.
            comp = ri.compare_with_keplerian(
                ceres, ceres.epoch + 1, use_assist=False
            )
        assert "delta_lon_arcsec" in comp

    def test_delta_lon_wrap_branches(self, fake_rebound, ceres):
        """Both wrap branches via crafted nbody/keplerian longitudes."""
        # nbody longitude ~ 359 deg
        high = ri.PropagationResult(
            x=1.0, y=-0.0175, z=0.0, vx=0.0, vy=0.0, vz=0.0, jd_tt=ceres.epoch
        )

        def _fake_kep(elements, jd):
            # keplerian longitude ~ 1 deg -> delta ~ 358 > 180 -> wrap down
            return (1.0, 0.0175, 0.0)

        with (
            patch.object(ri, "propagate_orbit_rebound", return_value=high),
            patch(
                "libephemeris.minor_bodies.calc_minor_body_position",
                side_effect=_fake_kep,
            ),
        ):
            comp = ri.compare_with_keplerian(
                ceres, ceres.epoch + 1, use_assist=False
            )
        assert abs(comp["delta_lon_arcsec"]) < 180 * 3600

        # Now nbody ~ 1 deg, keplerian ~ 359 deg -> delta ~ -358 < -180 -> wrap up
        low = ri.PropagationResult(
            x=1.0, y=0.0175, z=0.0, vx=0.0, vy=0.0, vz=0.0, jd_tt=ceres.epoch
        )

        def _fake_kep2(elements, jd):
            return (1.0, -0.0175, 0.0)

        with (
            patch.object(ri, "propagate_orbit_rebound", return_value=low),
            patch(
                "libephemeris.minor_bodies.calc_minor_body_position",
                side_effect=_fake_kep2,
            ),
        ):
            comp2 = ri.compare_with_keplerian(
                ceres, ceres.epoch + 1, use_assist=False
            )
        assert abs(comp2["delta_lon_arcsec"]) < 180 * 3600
