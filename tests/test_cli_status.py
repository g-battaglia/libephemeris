from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

import libephemeris.cli as cli_module
import libephemeris.dev_cli as dev_cli_module
import libephemeris.dev_cli.cmd_status as dev_status_module
import libephemeris.download as download_module
import libephemeris.iers_data as iers_data_module
import libephemeris.state as state_module


def _write_file(path: Path, size: int = 16) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"x" * size)


def _sample_dev_status_report() -> dict[str, object]:
    return {
        "version": "1.0.0-test",
        "project_root": "/repo",
        "data_directory": "/data",
        "download_prerequisites": {
            "spk": {
                "cache_dir": "/data/spk",
                "cache_files": [
                    {
                        "exists": True,
                        "name": "chiron_1_2.bsp",
                        "path": "/data/spk/chiron_1_2.bsp",
                        "size": 12,
                    }
                ],
                "tiers": {
                    "base": {
                        "status": "ok",
                        "de_kernel": {
                            "exists": True,
                            "name": "de440s.bsp",
                            "path": "/data/de440s.bsp",
                            "size": 10,
                        },
                        "range": {"start": "1900-01-01", "end": "2100-01-01"},
                        "cache_dir": "/data/spk",
                        "total_bodies": 2,
                        "required_bodies": 2,
                        "ready_bodies": 2,
                        "partial_bodies": 0,
                        "missing_bodies": 0,
                        "expected_fallback_bodies": 0,
                        "bodies": [
                            {
                                "name": "Chiron",
                                "status": "ready",
                                "path": "/data/spk/chiron_1_2.bsp",
                                "coverage": [1.0, 2.0],
                            },
                            {
                                "name": "Ceres",
                                "status": "ready",
                                "path": "/data/spk/ceres_1_2.bsp",
                                "coverage": [1.0, 2.0],
                            },
                        ],
                    },
                    "medium": {
                        "status": "partial",
                        "de_kernel": {
                            "exists": True,
                            "name": "de440.bsp",
                            "path": "/data/de440.bsp",
                            "size": 10,
                        },
                        "range": {"start": "1550-01-01", "end": "2650-01-01"},
                        "cache_dir": "/data/spk",
                        "total_bodies": 3,
                        "required_bodies": 2,
                        "ready_bodies": 1,
                        "partial_bodies": 0,
                        "missing_bodies": 1,
                        "expected_fallback_bodies": 1,
                        "bodies": [
                            {
                                "name": "Chiron",
                                "status": "ready",
                                "path": "/data/spk/chiron_1_2.bsp",
                                "coverage": [1.0, 2.0],
                            },
                            {
                                "name": "Ceres",
                                "status": "missing",
                                "path": None,
                                "coverage": None,
                            },
                            {
                                "name": "Bennu",
                                "status": "expected_fallback",
                                "path": None,
                                "coverage": None,
                                "reason": "JPL blocks SPK generation; Keplerian fallback is expected.",
                            },
                        ],
                    },
                    "extended": {
                        "status": "missing",
                        "de_kernel": {
                            "exists": False,
                            "name": "de441.bsp",
                            "path": "/data/de441.bsp",
                            "size": None,
                        },
                        "range": {"start": "1600-01-01", "end": "2500-01-01"},
                        "cache_dir": "/data/spk",
                        "total_bodies": 3,
                        "required_bodies": 2,
                        "ready_bodies": 0,
                        "partial_bodies": 0,
                        "missing_bodies": 2,
                        "expected_fallback_bodies": 1,
                        "bodies": [
                            {
                                "name": "Chiron",
                                "status": "missing",
                                "path": None,
                                "coverage": None,
                            },
                            {
                                "name": "Ceres",
                                "status": "missing",
                                "path": None,
                                "coverage": None,
                            },
                            {
                                "name": "Bennu",
                                "status": "expected_fallback",
                                "path": None,
                                "coverage": None,
                                "reason": "JPL blocks SPK generation; Keplerian fallback is expected.",
                            },
                        ],
                    },
                },
            },
            "iers": {
                "status": "partial",
                "cache_dir": "/data/iers",
                "present": 2,
                "total": 3,
                "files": {
                    "finals2000A.data": {
                        "exists": True,
                        "path": "/data/iers/finals2000A.data",
                        "age_days": 2.0,
                    },
                    "leap_seconds.dat": {
                        "exists": False,
                        "path": "/data/iers/leap_seconds.dat",
                        "age_days": None,
                    },
                    "deltat.data": {
                        "exists": True,
                        "path": "/data/iers/deltat.data",
                        "age_days": 4.0,
                    },
                },
            },
            "planet_center_sources": {
                "status": "partial",
                "cache_dir": "/data/tmp/planet_centers",
                "shared": {
                    "status": "ok",
                    "present": 1,
                    "total": 1,
                    "files": {
                        "naif0012.tls": {
                            "exists": True,
                            "name": "naif0012.tls",
                            "path": "/data/tmp/planet_centers/naif0012.tls",
                            "size": 8,
                        }
                    },
                },
                "tiers": {
                    "base": {
                        "status": "ok",
                        "present": 5,
                        "total": 5,
                        "files": {},
                    },
                    "medium": {
                        "status": "partial",
                        "present": 4,
                        "total": 5,
                        "files": {
                            "sat441.bsp": {
                                "exists": False,
                                "name": "sat441.bsp",
                                "path": "/data/tmp/planet_centers/sat441.bsp",
                                "size": None,
                            }
                        },
                    },
                    "extended": {
                        "status": "ok",
                        "present": 6,
                        "total": 6,
                        "files": {},
                    },
                },
            },
            "assist": {
                "status": "ok",
                "directory": "/data/assist",
                "present": 2,
                "total": 2,
                "files": {
                    "linux_p1550p2650.440": {
                        "exists": True,
                        "name": "linux_p1550p2650.440",
                        "path": "/data/assist/linux_p1550p2650.440",
                        "size": 10,
                    },
                    "sb441-n16.bsp": {
                        "exists": True,
                        "name": "sb441-n16.bsp",
                        "path": "/data/assist/sb441-n16.bsp",
                        "size": 10,
                    },
                },
            },
        },
        "generated_artifacts": {
            "leb1": {
                "base": {
                    "status": "ok",
                    "final": {
                        "exists": True,
                        "name": "ephemeris_base.leb",
                        "path": "/repo/data/leb/ephemeris_base.leb",
                        "size": 10,
                    },
                    "partial_present": 0,
                    "partial_total": 3,
                    "partials": {},
                },
                "medium": {
                    "status": "partial",
                    "final": {
                        "exists": False,
                        "name": "ephemeris_medium.leb",
                        "path": "/repo/data/leb/ephemeris_medium.leb",
                        "size": None,
                    },
                    "partial_present": 2,
                    "partial_total": 3,
                    "partials": {
                        "planets": {
                            "exists": True,
                            "name": "ephemeris_medium_planets.leb",
                            "path": "/repo/data/leb/ephemeris_medium_planets.leb",
                            "size": 10,
                        }
                    },
                },
                "extended": {
                    "status": "missing",
                    "final": {
                        "exists": False,
                        "name": "ephemeris_extended.leb",
                        "path": "/repo/data/leb/ephemeris_extended.leb",
                        "size": None,
                    },
                    "partial_present": 0,
                    "partial_total": 3,
                    "partials": {},
                },
            },
            "leb2": {
                "base": {
                    "status": "ok",
                    "present": 4,
                    "total": 4,
                    "files": {
                        "core": {
                            "exists": True,
                            "name": "base_core.leb2",
                            "path": "/repo/data/leb2/base_core.leb2",
                            "size": 10,
                        }
                    },
                },
                "medium": {
                    "status": "partial",
                    "present": 2,
                    "total": 4,
                    "files": {
                        "core": {
                            "exists": True,
                            "name": "medium_core.leb2",
                            "path": "/repo/data/leb2/medium_core.leb2",
                            "size": 10,
                        }
                    },
                },
                "extended": {
                    "status": "missing",
                    "present": 0,
                    "total": 4,
                    "files": {
                        "core": {
                            "exists": False,
                            "name": "extended_core.leb2",
                            "path": "/repo/data/leb2/extended_core.leb2",
                            "size": None,
                        }
                    },
                },
            },
            "planet_centers": {
                "base": {
                    "status": "ok",
                    "file": {
                        "exists": True,
                        "name": "planet_centers_base.bsp",
                        "path": "/data/planet_centers_base.bsp",
                        "size": 10,
                    },
                },
                "medium": {
                    "status": "missing",
                    "file": {
                        "exists": False,
                        "name": "planet_centers_medium.bsp",
                        "path": "/data/planet_centers_medium.bsp",
                        "size": None,
                    },
                },
                "extended": {
                    "status": "missing",
                    "file": {
                        "exists": False,
                        "name": "planet_centers_extended.bsp",
                        "path": "/data/planet_centers_extended.bsp",
                        "size": None,
                    },
                },
            },
        },
        "missing_commands": [
            "leph download spk-medium",
            "leph download iers",
            "leph leb generate medium groups",
            "leph generate planet-centers-medium",
        ],
    }


def test_download_all_includes_iers_downloads(monkeypatch):
    calls: list[tuple] = []

    def fake_download_leb2_for_tier(
        tier_name: str,
        force: bool,
        show_progress: bool,
        quiet: bool,
        activate: bool,
    ) -> None:
        calls.append(("leb2", tier_name, force, show_progress, quiet, activate))

    def fake_download_for_tier(
        tier_name: str,
        force: bool,
        show_progress: bool,
        quiet: bool,
    ) -> None:
        calls.append(("tier", tier_name, force, show_progress, quiet))

    monkeypatch.setattr(
        download_module,
        "download_leb2_for_tier",
        fake_download_leb2_for_tier,
    )
    monkeypatch.setattr(download_module, "download_for_tier", fake_download_for_tier)
    monkeypatch.setattr(
        iers_data_module,
        "download_iers_finals",
        lambda force=False: calls.append(("iers_finals", force)),
    )
    monkeypatch.setattr(
        iers_data_module,
        "download_leap_seconds",
        lambda force=False: calls.append(("leap_seconds", force)),
    )
    monkeypatch.setattr(
        iers_data_module,
        "download_delta_t_data",
        lambda force=False: calls.append(("delta_t", force)),
    )

    runner = CliRunner()
    result = runner.invoke(
        cli_module.cli,
        ["download", "all", "--force", "--no-progress", "--quiet", "--yes"],
    )

    assert result.exit_code == 0, result.output
    assert calls == [
        ("leb2", "base", True, False, True, False),
        ("tier", "base", True, False, True),
        ("leb2", "medium", True, False, True, False),
        ("tier", "medium", True, False, True),
        ("leb2", "extended", True, False, True, False),
        ("tier", "extended", True, False, True),
        ("iers_finals", True),
        ("leap_seconds", True),
        ("delta_t", True),
    ]


def test_download_auto_extended_installs_all_eligible_leb_tiers(monkeypatch):
    calls: list[tuple[object, ...]] = []

    monkeypatch.setattr(state_module, "get_precision_tier", lambda: "extended")
    monkeypatch.setattr(state_module, "get_calc_mode", lambda: "leb")
    monkeypatch.setattr(
        state_module,
        "set_leb_file",
        lambda path: calls.append(("activate", path)),
    )

    def fake_download_leb2_for_tier(**kwargs: object) -> None:
        calls.append(
            (
                "leb2",
                kwargs["tier_name"],
                kwargs["force"],
                kwargs["show_progress"],
                kwargs["quiet"],
                kwargs["activate"],
            )
        )

    monkeypatch.setattr(
        download_module,
        "download_leb2_for_tier",
        fake_download_leb2_for_tier,
    )

    result = CliRunner().invoke(
        cli_module.cli,
        ["download", "auto", "--force", "--no-progress", "--quiet"],
    )

    assert result.exit_code == 0, result.output
    assert calls == [
        ("leb2", "base", True, False, True, False),
        ("leb2", "medium", True, False, True, False),
        ("leb2", "extended", True, False, True, False),
        ("activate", None),
    ]


def test_status_verbose_levels(monkeypatch, tmp_path):
    home_dir = tmp_path / "home"
    data_dir = home_dir / ".libephemeris"
    leb_dir = data_dir / "leb"
    spk_dir = data_dir / "spk"
    assist_dir = home_dir / ".libephemeris" / "assist"
    iers_dir = data_dir / "iers_cache"

    _write_file(data_dir / "de440s.bsp", 32)
    _write_file(data_dir / "planet_centers_base.bsp", 48)
    _write_file(leb_dir / "ephemeris_medium.leb", 64)
    _write_file(leb_dir / "base_core.leb2", 12)
    _write_file(leb_dir / "base_asteroids.leb2", 12)
    _write_file(leb_dir / "base_apogee.leb2", 12)
    _write_file(leb_dir / "base_uranians.leb2", 12)
    _write_file(spk_dir / "alpha.bsp", 24)
    _write_file(spk_dir / "beta.bsp", 24)
    _write_file(assist_dir / "linux_p1550p2650.440", 20)
    _write_file(assist_dir / "sb441-n16.bsp", 20)
    _write_file(iers_dir / "finals2000A.data", 10)
    _write_file(iers_dir / "leap_seconds.dat", 10)
    _write_file(iers_dir / "deltat.data", 10)

    monkeypatch.setenv("HOME", str(home_dir))
    monkeypatch.setenv("LIBEPHEMERIS_LEB", str(leb_dir / "ephemeris_medium.leb"))
    monkeypatch.setattr(download_module, "get_data_dir", lambda: data_dir)
    monkeypatch.setattr(state_module, "get_calc_mode", lambda: "leb")
    monkeypatch.setattr(state_module, "get_precision_tier", lambda: "base")
    monkeypatch.setattr(state_module, "get_spk_cache_dir", lambda: str(spk_dir))
    monkeypatch.setattr(
        iers_data_module,
        "get_iers_cache_info",
        lambda: {
            "cache_dir": str(iers_dir),
            "finals_exists": True,
            "finals_age_days": 2.0,
            "leap_seconds_exists": True,
            "leap_seconds_age_days": 3.0,
            "delta_t_exists": True,
            "delta_t_age_days": 4.0,
        },
    )

    runner = CliRunner()

    base_result = runner.invoke(cli_module.cli, ["status"])
    assert base_result.exit_code == 0, base_result.output
    assert "Path:" not in base_result.output
    assert "leap_seconds.dat" in base_result.output

    verbose_result = runner.invoke(cli_module.cli, ["status", "-v"])
    assert verbose_result.exit_code == 0, verbose_result.output
    assert f"Path: {data_dir / 'de440s.bsp'}" in verbose_result.output
    assert f"Path: {leb_dir / 'base_core.leb2'}" in verbose_result.output
    assert "base_core.leb2" in verbose_result.output

    very_verbose_result = runner.invoke(cli_module.cli, ["status", "-vv"])
    assert very_verbose_result.exit_code == 0, very_verbose_result.output
    assert f"Path: {spk_dir / 'alpha.bsp'}" in very_verbose_result.output
    assert f"Path: {spk_dir / 'beta.bsp'}" in very_verbose_result.output

    json_result = runner.invoke(cli_module.cli, ["status", "--json"])
    assert json_result.exit_code == 0, json_result.output
    parsed = json.loads(json_result.output)
    assert parsed["data_directory"] == str(data_dir)
    assert parsed["calc_mode"] == "leb"
    assert parsed["precision_tier"] == "base"


def test_leph_status_reports_development_status(monkeypatch):
    monkeypatch.setattr(
        dev_status_module,
        "_build_dev_status_report",
        _sample_dev_status_report,
    )

    runner = CliRunner()
    result = runner.invoke(dev_cli_module.cli, ["status"])

    assert result.exit_code == 0, result.output
    assert "libephemeris dev 1.0.0-test" in result.output
    assert (
        "spk-medium: de440.bsp + 1/2 downloadable bodies, 1 expected fallback"
        in result.output
    )
    assert "LEB1 medium: ephemeris_medium.leb (partials 2/3)" in result.output
    assert "leph download iers" in result.output


def test_leph_status_verbose_modes_and_json(monkeypatch):
    monkeypatch.setattr(
        dev_status_module,
        "_build_dev_status_report",
        _sample_dev_status_report,
    )

    runner = CliRunner()

    verbose_result = runner.invoke(dev_cli_module.cli, ["status", "-v"])
    assert verbose_result.exit_code == 0, verbose_result.output
    assert "Kernel: /data/de440.bsp" in verbose_result.output
    assert "Missing bodies: Ceres" in verbose_result.output
    assert (
        "Expected fallback body: Bennu (JPL blocks SPK generation; Keplerian fallback is expected.)"
        in verbose_result.output
    )
    assert "Path: /repo/data/leb/ephemeris_medium.leb" in verbose_result.output

    very_verbose_result = runner.invoke(dev_cli_module.cli, ["status", "-vv"])
    assert very_verbose_result.exit_code == 0, very_verbose_result.output
    assert "Coverage: JD 1.0 -> 2.0" in very_verbose_result.output
    assert "Bennu (fallback)" in very_verbose_result.output
    assert "Path: /data/tmp/planet_centers/sat441.bsp" in very_verbose_result.output
    assert "Path: /repo/data/leb2/medium_core.leb2" in very_verbose_result.output

    json_result = runner.invoke(dev_cli_module.cli, ["status", "--json"])
    assert json_result.exit_code == 0, json_result.output
    parsed = json.loads(json_result.output)
    assert parsed["version"] == "1.0.0-test"
    assert parsed["project_root"] == "/repo"
    assert parsed["missing_commands"][0] == "leph download spk-medium"
