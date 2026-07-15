"""Regression tests for the top-level LEB regeneration orchestrator."""

from __future__ import annotations

import os
from pathlib import Path
import re
import subprocess

from libephemeris.leb_groups import LEB2_GROUPS


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "regenerate-leb.sh"


def _run_sourced_script(command: str) -> subprocess.CompletedProcess[str]:
    """Load the Bash functions without invoking ``main`` and run a probe."""
    env = {**os.environ, "SCRIPT_UNDER_TEST": str(SCRIPT)}
    return subprocess.run(
        ["/bin/bash", "-c", 'source "$SCRIPT_UNDER_TEST"; ' + command],
        cwd=PROJECT_ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )


def _recorded_arguments(result: subprocess.CompletedProcess[str]) -> list[str]:
    """Return individual arguments emitted by the probe's mocked ``run``."""
    return [
        line.removeprefix("ARG=")
        for line in result.stdout.splitlines()
        if line.startswith("ARG=")
    ]


def test_shell_leb2_groups_match_canonical_python_registry() -> None:
    """The standalone Bash entrypoint cannot retain a retired companion."""
    source = SCRIPT.read_text(encoding="utf-8")
    match = re.search(r"^LEB2_GROUPS=\(([^)]*)\)$", source, flags=re.MULTILINE)

    assert match is not None
    assert tuple(match.group(1).split()) == LEB2_GROUPS


def test_all_selects_only_canonical_leb2_groups() -> None:
    """``all`` expands in the canonical conversion/distribution order."""
    result = _run_sourced_script("LEB2_GROUP_SPEC=all; selected_leb2_groups")

    assert result.returncode == 0, result.stderr
    assert tuple(result.stdout.splitlines()) == LEB2_GROUPS


def test_retired_uranians_group_is_rejected() -> None:
    """A stale caller must fail before attempting to generate a bogus file."""
    result = _run_sourced_script("parse_args --leb2-only --leb2-groups uranians")

    assert result.returncode == 2
    assert "invalid LEB2 group: uranians" in result.stderr


def test_analytical_leb1_affects_only_surviving_leb2_groups() -> None:
    """Nodes and apsides map to core/apogee after hypothetical-body removal."""
    result = _run_sourced_script(
        'LEB2_GROUP_SPEC=""; GROUP_SPEC=analytical; affected_leb2_groups'
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == ["core", "apogee"]


def test_conversion_propagates_group_and_extended_tier() -> None:
    """Tier-specific inventory authentication reaches the converter CLI."""
    result = _run_sourced_script(
        """
        PYTHON=python-under-test
        DRY_RUN=1
        QUIET=0
        mkdir() { :; }
        backup_path_if_exists() { :; }
        selected_leb2_groups() { printf '%s\n' exotics; }
        run() { printf 'ARG=%s\n' "$@"; }
        convert_leb2_tier extended
        """
    )

    assert result.returncode == 0, result.stderr
    assert _recorded_arguments(result) == [
        "python-under-test",
        "scripts/generate_leb2.py",
        "convert",
        "data/leb/ephemeris_extended.leb",
        "-o",
        "data/leb2/extended_exotics.leb2",
        "--group",
        "exotics",
        "--tier",
        "extended",
    ]


def test_verification_authenticates_group_and_extended_tier() -> None:
    """Verification rejects a readable file with the wrong companion inventory."""
    result = _run_sourced_script(
        """
        PYTHON=python-under-test
        VERIFY=1
        QUIET=0
        LEB2_VERIFY_SAMPLES=17
        selected_leb2_groups() { printf '%s\n' exotics; }
        run() { printf 'ARG=%s\n' "$@"; }
        verify_leb2_tier extended
        """
    )

    assert result.returncode == 0, result.stderr
    assert _recorded_arguments(result) == [
        "python-under-test",
        "scripts/generate_leb2.py",
        "verify",
        "data/leb2/extended_exotics.leb2",
        "--reference",
        "data/leb/ephemeris_extended.leb",
        "--samples",
        "17",
        "--group",
        "exotics",
        "--tier",
        "extended",
    ]
