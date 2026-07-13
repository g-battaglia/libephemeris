"""Regression tests for backend selection in the developer test CLI."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner
import pytest

from libephemeris import dev_cli
from libephemeris.cli import cli as user_cli
from libephemeris.cli import init_wizard
from libephemeris.dev_cli import (
    cmd_completion,
    cmd_leb,
    cmd_leb2,
    cmd_release,
    cmd_test,
)
from libephemeris.exotic_bodies import EXOTIC_EXTENDED_IDS, EXOTIC_IDS
from libephemeris.leb_groups import LEB1_GROUPS, LEB2_GROUPS
from scripts import generate_leb2, release_leb, test_leb2_precision


def test_skyfield_runner_forces_skyfield_backend() -> None:
    """The command group must not inherit an auto-selected LEB backend."""
    with patch.object(cmd_test, "_pytest") as run_pytest:
        cmd_test._skyfield_pytest(["tests/example.py", "-q"])

    run_pytest.assert_called_once_with(
        ["tests/example.py", "-q", "--calc-mode", "skyfield"]
    )


def test_completion_scripts_use_repo_pinned_module_launcher() -> None:
    """Generated functions must not delegate to the removed wheel entry point."""
    runner = CliRunner()

    for shell in ("zsh", "bash", "fish"):
        result = runner.invoke(cmd_completion.completion_group, [shell])
        assert result.exit_code == 0
        assert "uv run leph" not in result.output
        assert "uv run --project" in result.output
        assert "python -m libephemeris.dev_cli" in result.output
        assert str(cmd_completion._PROJECT_ROOT) in result.output
        if shell == "bash":
            assert "IFS=',' read -r type value" in result.output
            assert 'COMPREPLY+=("$value")' in result.output
        elif shell == "fish":
            assert 'set -l metadata (string split "," $completion)' in result.output
            assert "echo $metadata[2]" in result.output


def test_generated_bash_completion_strips_click_metadata() -> None:
    """Bash TAB candidates contain values, never Click's ``plain,`` prefix."""
    script = cmd_completion._render_script("bash")
    probe = """
COMP_WORDS=(leph test '')
COMP_CWORD=2
_leph_completion
printf '%s\n' "${COMPREPLY[@]}"
"""

    result = subprocess.run(
        ["/bin/bash"],
        input=script + probe,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    candidates = result.stdout.splitlines()
    assert "skyfield" in candidates
    assert all(not candidate.startswith("plain,") for candidate in candidates)


def test_generated_zsh_completion_initializes_compinit_idempotently(tmp_path) -> None:
    """A pristine zsh can eval the script twice without touching the real home."""
    zsh = shutil.which("zsh")
    if zsh is None:
        pytest.skip("zsh is required for the generated zsh completion probe")
    script = cmd_completion._render_script("zsh")
    probe = """
eval "$LEPH_COMPLETION_SCRIPT"
first_rc=$?
eval "$LEPH_COMPLETION_SCRIPT"
second_rc=$?
print -r -- "first_rc=$first_rc"
print -r -- "second_rc=$second_rc"
whence -w compdef
whence -w _leph_completion
whence -w leph
"""
    env = {
        **os.environ,
        "HOME": str(tmp_path),
        "ZDOTDIR": str(tmp_path),
        "LEPH_COMPLETION_SCRIPT": script,
    }

    result = subprocess.run(
        [zsh, "-f"],
        input=probe,
        text=True,
        capture_output=True,
        check=False,
        env=env,
    )

    assert result.returncode == 0, result.stderr
    assert "first_rc=0" in result.stdout
    assert "second_rc=0" in result.stdout
    assert "compdef: function" in result.stdout
    assert "_leph_completion: function" in result.stdout
    assert "leph: function" in result.stdout


def test_module_launcher_uses_stable_click_program_name() -> None:
    """Click completion must consistently use the ``_LEPH_COMPLETE`` protocol."""
    with (
        patch.object(dev_cli.os, "chdir") as chdir,
        patch.object(dev_cli, "cli") as root_cli,
    ):
        dev_cli.main()

    chdir.assert_called_once_with(cmd_completion._PROJECT_ROOT)
    root_cli.assert_called_once_with(prog_name="leph")


def test_leb2_generator_group_order_matches_canonical_registry() -> None:
    """Conversion order cannot drift from runtime companion discovery order."""
    assert tuple(generate_leb2.LEB2_GROUPS) == LEB2_GROUPS


def test_leb1_generate_exposes_every_tier_group_combination() -> None:
    """Every canonical LEB1 group has a direct command for every tier."""
    for tier in ("base", "medium", "extended"):
        tier_group = cmd_leb.generate_group.commands[tier]
        assert set(LEB1_GROUPS) <= set(tier_group.commands)


def test_leb1_exotics_command_dispatches_registry_group() -> None:
    """The direct exotic command selects the matching tier and group."""
    runner = CliRunner()
    with patch.object(cmd_leb, "_gen") as run_generator:
        result = runner.invoke(cmd_leb.generate_group, ["medium", "exotics"])

    assert result.exit_code == 0
    run_generator.assert_called_once_with(["--tier", "medium", "--group", "exotics"])


def test_leb_help_reports_eligible_inventories_and_preserves_examples() -> None:
    """Generated help distinguishes full registries and keeps commands separate."""
    result = CliRunner().invoke(cmd_leb.leb_group, ["--help"])

    assert result.exit_code == 0
    assert "53 independently sourced bodies" in result.output
    assert "45 bodies (8 NEAs excluded)" in result.output
    assert (
        "leph leb generate medium groups   # all body groups + merge\n"
        "    leph leb verify medium"
    ) in result.output


def test_leb_generate_help_preserves_mode_lines() -> None:
    """Click must not collapse the generation-mode inventory into prose."""
    result = CliRunner().invoke(cmd_leb.generate_group, ["--help"])

    assert result.exit_code == 0
    assert "then merge\n    full" in result.output
    assert "macOS)\n    single" in result.output
    assert "slowest)\n    body <name>" in result.output


def test_direct_leb1_help_includes_exotics_in_group_and_merge_inventory() -> None:
    """The direct generator must not teach an incomplete three-group merge."""
    result = subprocess.run(
        [sys.executable, "scripts/generate_leb.py", "--help"],
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0
    assert "exotics    : Centaurs" in result.stdout
    assert "planets.leb asteroids.leb exotics.leb analytical.leb" in result.stdout


def test_leb2_convert_exposes_every_tier_group_combination() -> None:
    """The help promise of per-group control is backed by concrete commands."""
    expected = {
        f"{tier}-{group}"
        for tier in ("base", "medium", "extended")
        for group in LEB2_GROUPS
    }

    assert expected <= set(cmd_leb2.convert_group.commands)


def test_leb2_root_help_preserves_separate_examples() -> None:
    """Click must not reflow the convert and core-verification examples together."""
    result = CliRunner().invoke(cmd_leb2.leb2_group, ["--help"])

    assert result.exit_code == 0
    assert "all 4 body groups" in result.output
    assert ("exotics, apogee\n    leph leb2 verify base") in result.output


def test_leb2_dynamic_group_command_uses_matching_paths() -> None:
    """Dynamic group conversion passes the selected tier and group exactly."""
    runner = CliRunner()
    with patch.object(cmd_leb2, "_leb2") as run_leb2:
        result = runner.invoke(cmd_leb2.convert_group, ["medium-exotics"])

    assert result.exit_code == 0
    run_leb2.assert_called_once_with(
        [
            "convert",
            "data/leb/ephemeris_medium.leb",
            "-o",
            "data/leb2/medium_exotics.leb2",
            "--group",
            "exotics",
            "--tier",
            "medium",
        ]
    )


def test_leb2_uranian_conversion_is_retired_fail_closed() -> None:
    """Legacy hypothetical coefficients cannot be copied into a new artifact."""
    result = CliRunner().invoke(cmd_leb2.convert_group, ["base-uranians"])

    assert result.exit_code != 0
    assert "retired pending independently sourced elements" in result.output


def test_leb2_conversion_rejects_empty_named_group(monkeypatch, tmp_path) -> None:
    """A missing group fails instead of writing an auxiliary-only LEB2 file."""

    class EmptySource:
        bodies: dict[int, object] = {}
        closed = False

        def __init__(self, _path) -> None:
            pass

        def close(self) -> None:
            self.closed = True

    monkeypatch.setattr(generate_leb2, "LEB1Source", EmptySource)
    output = tmp_path / "empty.leb2"

    with pytest.raises(ValueError, match="no bodies.*exotics"):
        generate_leb2.convert_leb1_to_leb2(
            "legacy.leb", str(output), group="exotics", verbose=False
        )
    assert not output.exists()


def test_leb2_conversion_rejects_incomplete_nonempty_group(
    monkeypatch, tmp_path
) -> None:
    """A named output cannot silently contain a non-empty group subset."""

    class PartialSource:
        bodies = {0: object()}

        def __init__(self, _path) -> None:
            pass

        def close(self) -> None:
            pass

    monkeypatch.setattr(generate_leb2, "LEB1Source", PartialSource)
    output = tmp_path / "partial_core.leb2"

    with pytest.raises(ValueError, match="invalid 'core' body inventory"):
        generate_leb2.convert_leb1_to_leb2(
            "partial.leb", str(output), group="core", verbose=False
        )
    assert not output.exists()


def test_leb2_exotics_inventory_is_authenticated_by_tier() -> None:
    """Only an explicitly extended companion may omit the eight NEAs."""
    assert generate_leb2._group_inventory_error("exotics", EXOTIC_IDS) is None
    assert generate_leb2._group_inventory_error("exotics", EXOTIC_IDS, "base") is None
    assert generate_leb2._group_inventory_error("exotics", EXOTIC_IDS, "medium") is None
    assert generate_leb2._group_inventory_error("exotics", EXOTIC_EXTENDED_IDS)
    assert generate_leb2._group_inventory_error("exotics", EXOTIC_EXTENDED_IDS, "base")
    assert generate_leb2._group_inventory_error(
        "exotics", EXOTIC_EXTENDED_IDS, "medium"
    )
    assert (
        generate_leb2._group_inventory_error("exotics", EXOTIC_EXTENDED_IDS, "extended")
        is None
    )
    assert generate_leb2._group_inventory_error("exotics", EXOTIC_IDS, "extended")
    assert generate_leb2._group_inventory_error("exotics", EXOTIC_IDS[:-1])


def test_leb2_tier_verify_dispatches_core_file_verifier() -> None:
    """Tier verification explicitly checks the named core companion."""
    runner = CliRunner()
    with patch.object(cmd_leb2, "_leb2") as run_leb2:
        result = runner.invoke(cmd_leb2.verify_group, ["base"])

    assert result.exit_code == 0
    run_leb2.assert_called_once_with(
        [
            "verify",
            "data/leb2/base_core.leb2",
            "--reference",
            "data/leb/ephemeris_base.leb",
            "--samples",
            "500",
            "--group",
            "core",
            "--tier",
            "base",
        ]
    )


def test_leb2_verify_help_describes_core_component_check() -> None:
    """Verification help labels both its scope and the reported unit."""
    result = CliRunner().invoke(cmd_leb2.verify_group, ["--help"])

    assert result.exit_code == 0
    assert "core LEB2 companion" in result.output
    assert "native stored component units" in result.output
    assert "arcsecond" not in result.output


def test_leb2_verification_bound_uses_compression_target() -> None:
    """The verifier derives its native-unit bound from coefficient targets."""
    assert generate_leb2._verification_tolerance(1, 13) == 14e-12
    assert (
        generate_leb2._verification_tolerance(5, 13)
        == 14 * generate_leb2.DEFAULT_TARGET_AU
    )


def test_leb2_native_component_errors_preserve_units_and_wrap() -> None:
    """Ecliptic longitude wraps while latitude and distance keep native units."""
    errors, units = generate_leb2._native_component_errors(
        generate_leb2.COORD_ECLIPTIC,
        (359.99999999, 2.0, 1.25),
        (0.00000001, 2.5, 1.5),
    )
    assert errors == pytest.approx((2e-8, 0.5, 0.25))
    assert units == ("deg", "deg", "AU")


def test_leb2_verifier_rejects_empty_file(monkeypatch) -> None:
    """An empty but readable container cannot pass either verification mode."""
    from libephemeris import leb2_reader

    class EmptyReader:
        _bodies: dict[int, object] = {}
        closed = False

        def close(self) -> None:
            self.closed = True

    reader = EmptyReader()
    monkeypatch.setattr(leb2_reader, "LEB2Reader", lambda _path: reader)

    assert not generate_leb2.verify_leb2("empty.leb2", verbose=False)
    assert reader.closed


def test_leb2_verifier_rejects_incomplete_declared_core(monkeypatch) -> None:
    """The tier verifier requires the complete canonical 14-body core."""
    from libephemeris import leb2_reader

    class Body:
        pass

    class PartialReader:
        _bodies = {0: Body()}

        def close(self) -> None:
            pass

    monkeypatch.setattr(leb2_reader, "LEB2Reader", lambda _path: PartialReader())

    assert not generate_leb2.verify_leb2(
        "base_core.leb2", verbose=False, expected_group="core"
    )


def test_leb2_verifier_rejects_non_finite_positions(monkeypatch) -> None:
    """NaN position components cannot satisfy an AU error bound."""
    from libephemeris import leb2_reader, leb_reader

    class Body:
        jd_start = 2450000.0
        jd_end = 2460000.0
        degree = 0
        segment_count = 313
        interval_days = 32.0
        components = 3
        coord_type = generate_leb2.COORD_ICRS_BARY

    class LEB2Reader:
        _bodies = {0: Body()}
        jd_range = (2450000.0, 2460000.0)

        def eval_body(self, _bid, _jd):
            return (float("nan"), 0.0, 0.0), (0.0, 0.0, 0.0)

        def close(self) -> None:
            pass

    class LEB1Reader:
        _bodies = {0: Body()}
        jd_range = (2450000.0, 2460000.0)

        def has_body(self, _bid):
            return True

        def eval_body(self, _bid, _jd):
            return (1.0, 0.0, 0.0), (0.0, 0.0, 0.0)

        def close(self) -> None:
            pass

    monkeypatch.setattr(leb2_reader, "LEB2Reader", lambda _path: LEB2Reader())
    monkeypatch.setattr(leb_reader, "LEBReader", lambda _path: LEB1Reader())

    assert not generate_leb2.verify_leb2(
        "nan.leb2", "reference.leb", n_samples=1, verbose=False
    )


def test_leb2_verifier_rejects_wrong_reference_metadata(monkeypatch) -> None:
    """A forged coord_type cannot select the wrong units and still pass."""
    from libephemeris import leb2_reader, leb_reader

    class LEB2Body:
        jd_start = 2450000.0
        jd_end = 2460000.0
        degree = 13
        segment_count = 313
        interval_days = 32.0
        components = 3
        coord_type = generate_leb2.COORD_ECLIPTIC

    class LEB1Body:
        jd_start = 2450000.0
        jd_end = 2460000.0
        degree = 13
        segment_count = 313
        interval_days = 32.0
        components = 3
        coord_type = generate_leb2.COORD_ICRS_BARY

    class LEB2Reader:
        _bodies = {0: LEB2Body()}
        jd_range = (2450000.0, 2460000.0)
        closed = False

        def close(self) -> None:
            self.closed = True

    class LEB1Reader:
        _bodies = {0: LEB1Body()}
        jd_range = (2450000.0, 2460000.0)
        closed = False

        def has_body(self, _bid):
            return True

        def close(self) -> None:
            self.closed = True

    reader2 = LEB2Reader()
    reader1 = LEB1Reader()
    monkeypatch.setattr(leb2_reader, "LEB2Reader", lambda _path: reader2)
    monkeypatch.setattr(leb_reader, "LEBReader", lambda _path: reader1)

    assert not generate_leb2.verify_leb2(
        "wrong-frame.leb2", "reference.leb", n_samples=1, verbose=False
    )
    assert reader1.closed
    assert reader2.closed


def test_leb2_verifier_rejects_wrong_segment_count(monkeypatch) -> None:
    """Compressed per-body coverage cannot authenticate with fewer segments."""
    from libephemeris import leb2_reader, leb_reader

    class LEB2Body:
        jd_start = 2450000.0
        jd_end = 2460000.0
        degree = 13
        segment_count = 1
        interval_days = 32.0
        components = 3
        coord_type = generate_leb2.COORD_ICRS_BARY

    class LEB1Body(LEB2Body):
        segment_count = 313

    class LEB2Reader:
        _bodies = {0: LEB2Body()}
        jd_range = (2450000.0, 2460000.0)

        def close(self) -> None:
            pass

    class LEB1Reader:
        _bodies = {0: LEB1Body()}
        jd_range = (2450000.0, 2460000.0)

        def has_body(self, _bid):
            return True

        def close(self) -> None:
            pass

    monkeypatch.setattr(leb2_reader, "LEB2Reader", lambda _path: LEB2Reader())
    monkeypatch.setattr(leb_reader, "LEBReader", lambda _path: LEB1Reader())

    assert not generate_leb2.verify_leb2(
        "wrong-segments.leb2", "reference.leb", n_samples=1, verbose=False
    )


def test_leb2_verifier_rejects_wrong_file_coverage(monkeypatch) -> None:
    """The LEB2 file header must authenticate against the LEB1 header."""
    from libephemeris import leb2_reader, leb_reader

    class Body:
        jd_start = 2450000.0
        jd_end = 2460000.0
        degree = 13
        segment_count = 313
        interval_days = 32.0
        components = 3
        coord_type = generate_leb2.COORD_ICRS_BARY

    class LEB2Reader:
        _bodies = {0: Body()}
        jd_range = (2450001.0, 2460000.0)

        def eval_body(self, _bid, _jd):
            return (1.0, 0.0, 0.0), (0.0, 0.0, 0.0)

        def close(self) -> None:
            pass

    class LEB1Reader:
        _bodies = {0: Body()}
        jd_range = (2450000.0, 2460000.0)

        def has_body(self, _bid):
            return True

        def eval_body(self, _bid, _jd):
            return (1.0, 0.0, 0.0), (0.0, 0.0, 0.0)

        def close(self) -> None:
            pass

    monkeypatch.setattr(leb2_reader, "LEB2Reader", lambda _path: LEB2Reader())
    monkeypatch.setattr(leb_reader, "LEBReader", lambda _path: LEB1Reader())

    assert not generate_leb2.verify_leb2(
        "wrong-range.leb2", "reference.leb", n_samples=1, verbose=False
    )


def test_leb2_smoke_verifier_rejects_non_finite_values(monkeypatch) -> None:
    """The no-reference smoke path must also fail closed on NaN values."""
    from libephemeris import leb2_reader

    class Body:
        jd_start = 2450000.0
        jd_end = 2460000.0

    class LEB2Reader:
        _bodies = {0: Body()}
        jd_range = (2450000.0, 2460000.0)

        def eval_body(self, _bid, _jd):
            return (float("nan"), 0.0, 0.0), (0.0, 0.0, 0.0)

        def close(self) -> None:
            pass

    monkeypatch.setattr(leb2_reader, "LEB2Reader", lambda _path: LEB2Reader())

    assert not generate_leb2.verify_leb2("nan.leb2", verbose=False)


def test_leb2_precision_missing_inputs_fail_closed(monkeypatch) -> None:
    """A missing required ephemeris cannot be reported as a passing skip."""
    monkeypatch.setattr(test_leb2_precision.os.path, "isfile", lambda _path: False)

    assert not test_leb2_precision.run_test("base", n_dates=1)


def test_leb2_precision_uses_file_header_ranges() -> None:
    """Tier config contains paths only; date coverage comes from file headers."""
    for config in test_leb2_precision.TIER_CONFIG.values():
        assert set(config) == {"leb1", "leb2"}


def test_leb2_precision_wraps_arbitrary_finite_longitudes() -> None:
    """Out-of-range finite longitudes cannot turn into a negative error."""
    assert test_leb2_precision._angular_error_arcsec(1000.0, 0.0) == 80 * 3600
    assert test_leb2_precision._angular_error_arcsec(359.0, 1.0) == 2 * 3600


def test_leb2_precision_incomplete_core_fails_closed(monkeypatch) -> None:
    """A readable core with fewer than all 14 bodies cannot pass."""

    class IncompleteReader:
        _bodies = {0: object()}
        jd_range = (2450000.0, 2460000.0)

        def close(self) -> None:
            pass

    monkeypatch.setattr(test_leb2_precision.os.path, "isfile", lambda _path: True)
    monkeypatch.setattr(
        test_leb2_precision, "open_leb", lambda _path: IncompleteReader()
    )

    assert not test_leb2_precision.run_test("base", n_dates=1)


def test_leb2_precision_calculation_errors_fail_closed(monkeypatch) -> None:
    """Backend exceptions and zero completed comparisons are failures."""

    class CoreReader:
        _bodies = {bid: object() for bid in test_leb2_precision.CORE_BODY_IDS}
        jd_range = (2450000.0, 2460000.0)

        def close(self) -> None:
            pass

    monkeypatch.setattr(test_leb2_precision.os.path, "isfile", lambda _path: True)
    monkeypatch.setattr(test_leb2_precision, "open_leb", lambda _path: CoreReader())
    monkeypatch.setattr(test_leb2_precision.swe, "set_leb_file", lambda _path: None)
    monkeypatch.setattr(test_leb2_precision.swe, "set_calc_mode", lambda _mode: None)
    monkeypatch.setattr(test_leb2_precision.swe, "close", lambda: None)

    def fail_calc(*_args, **_kwargs):
        raise RuntimeError("synthetic calculation failure")

    monkeypatch.setattr(test_leb2_precision.swe, "calc_ut", fail_calc)

    assert not test_leb2_precision.run_test("base", n_dates=1)


def test_leb2_precision_non_finite_results_fail_closed(monkeypatch) -> None:
    """A complete body inventory still fails when calculations return NaN."""

    class CoreReader:
        _bodies = {bid: object() for bid in test_leb2_precision.CORE_BODY_IDS}
        jd_range = (2450000.0, 2460000.0)

        def close(self) -> None:
            pass

    monkeypatch.setattr(test_leb2_precision.os.path, "isfile", lambda _path: True)
    monkeypatch.setattr(test_leb2_precision, "open_leb", lambda _path: CoreReader())
    monkeypatch.setattr(test_leb2_precision.swe, "set_leb_file", lambda _path: None)
    monkeypatch.setattr(test_leb2_precision.swe, "set_calc_mode", lambda _mode: None)
    monkeypatch.setattr(test_leb2_precision.swe, "close", lambda: None)
    monkeypatch.setattr(
        test_leb2_precision.swe,
        "calc_ut",
        lambda *_args, **_kwargs: ((float("nan"), 0.0, 1.0, 0.0, 0.0, 0.0), 0),
    )

    assert not test_leb2_precision.run_test("base", n_dates=1)


def test_leb2_precision_short_position_shape_fails_closed(monkeypatch) -> None:
    """Two-component backend output cannot pass the three-component gate."""

    class CoreReader:
        _bodies = {bid: object() for bid in test_leb2_precision.CORE_BODY_IDS}
        jd_range = (2450000.0, 2460000.0)

        def close(self) -> None:
            pass

    monkeypatch.setattr(test_leb2_precision.os.path, "isfile", lambda _path: True)
    monkeypatch.setattr(test_leb2_precision, "open_leb", lambda _path: CoreReader())
    monkeypatch.setattr(test_leb2_precision.swe, "set_leb_file", lambda _path: None)
    monkeypatch.setattr(test_leb2_precision.swe, "set_calc_mode", lambda _mode: None)
    monkeypatch.setattr(test_leb2_precision.swe, "close", lambda: None)
    monkeypatch.setattr(
        test_leb2_precision.swe,
        "calc_ut",
        lambda *_args, **_kwargs: ((1.0, 2.0), 0),
    )

    assert not test_leb2_precision.run_test("base", n_dates=1)


def test_user_cli_installs_pinned_cores_and_limits_modular_assets(
    monkeypatch, tmp_path
) -> None:
    """Historical command names install current pinned cores."""
    import libephemeris.download as download_module

    calls = []

    def _download_core(tier_name, **kwargs):
        calls.append((tier_name, kwargs))
        return tmp_path / f"{tier_name}_core.leb2"

    monkeypatch.setattr(download_module, "download_leb_for_tier", _download_core)
    runner = CliRunner()
    config_result = runner.invoke(user_cli, ["config"])
    group_result = runner.invoke(user_cli, ["download", "--help"])
    leb1_help_result = runner.invoke(user_cli, ["download", "leb-base", "--help"])
    leb1_result = runner.invoke(user_cli, ["download", "leb-base"])
    download_result = runner.invoke(user_cli, ["download", "leb2-base", "--help"])

    assert config_result.exit_code == 0
    assert "All tiers have an exact SHA-256-pinned core asset" in config_result.output
    assert "Canonical LEB2 group names (not all are published)" in config_result.output
    assert "{tier}_exotics.leb2" in config_result.output
    assert group_result.exit_code == 0
    assert (
        "Historical leb-* commands install the pinned core in the current format"
        in group_result.output
    )
    assert leb1_help_result.exit_code == 0
    assert (
        "Install the reviewed, SHA-256-pinned core LEB asset" in leb1_help_result.output
    )
    assert "Reinstall even if already present" in leb1_help_result.output
    assert leb1_result.exit_code == 0
    assert calls and calls[0][0] == "base"
    assert download_result.exit_code == 0
    assert "Only SHA-256-pinned entries" in download_result.output
    normalized_download_help = " ".join(download_result.output.split())
    assert (
        "Other precomputed companions are skipped pending artifact-specific review."
        in normalized_download_help
    )


@pytest.mark.parametrize(
    "command",
    ["leb", "leb-base", "leb-medium", "leb-extended", "leb-dry-run"],
)
def test_dev_release_commands_are_retired(command: str) -> None:
    """Every legacy developer release command remains parseable but inert."""
    result = CliRunner().invoke(cmd_release.release_group, [command, "1.0.0"])

    assert result.exit_code == 1
    assert "retired pending clean-room regeneration" in result.output


def test_release_script_is_retired(capsys: pytest.CaptureFixture[str]) -> None:
    """The standalone release script rejects even dry-run invocations."""
    result = release_leb.main(["--version", "1.0.0", "--dry-run"])

    assert result == 1
    assert "retired pending clean-room regeneration" in capsys.readouterr().err


def test_skyfield_test_help_does_not_claim_default_backend() -> None:
    """The backend-specific test help must reflect default auto mode."""
    result = CliRunner().invoke(dev_cli.cli, ["test", "skyfield", "--help"])

    assert result.exit_code == 0
    assert "local reference backend" in result.output
    assert "default calculation engine" not in result.output


def test_user_cli_describes_every_toml_configuration_family() -> None:
    """The config guide covers Delta T selection and mmap preload keys."""
    result = CliRunner().invoke(user_cli, ["config"])

    assert result.exit_code == 0
    assert "LIBEPHEMERIS_DELTAT_MODEL" in result.output
    assert 'deltat_model = "smh2016"' in result.output
    assert "mmap_preload = false" in result.output
    assert "mmap_preload_start = 1800" in result.output
    assert "mmap_preload_end = 2200" in result.output


def test_init_wizard_generates_every_toml_configuration_family() -> None:
    """Generated TOML includes valid Delta T and mmap settings."""
    content = init_wizard._generate_toml({})

    assert 'deltat_model = "smh2016"' in content
    assert "mmap_preload = false" in content
    assert "mmap_preload_start = 1800" in content
    assert "mmap_preload_end = 2200" in content


def test_init_wizard_does_not_require_unpublished_exotics() -> None:
    """Wizard requires only the reviewed bundled base core."""
    files = init_wizard._get_required_files("base", "leb")
    leb2 = {
        item["filename"]: item["status"] for item in files if item["category"] == "LEB2"
    }

    assert leb2["base_exotics.leb2"] == "optional"
    assert sum(status == "required" for status in leb2.values()) == 1


def test_init_wizard_estimates_only_reviewed_artifacts() -> None:
    """Only manifest-pinned core sizes are presented as active estimates."""
    assert init_wizard._LEB2_SIZES == {
        "base": {"core": 10.8},
        "medium": {"core": 38.3},
        "extended": {"core": 334.9},
    }


def test_init_wizard_prefers_reviewed_leb2_in_selected_directory(
    tmp_path: Path,
) -> None:
    """Directory selection should prefer the current core over legacy LEB1."""
    reviewed = tmp_path / "medium_core.leb2"
    reviewed.touch()
    (tmp_path / "ephemeris_medium.leb").touch()

    assert init_wizard._resolve_leb_entry(str(tmp_path), "medium") == str(reviewed)
