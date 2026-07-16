"""Regression tests for the top-level LEB regeneration orchestrator."""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

from libephemeris.leb_groups import LEB2_GROUPS

# Groups the Bash orchestrator handles: exactly the ones converted from the
# merged main LEB1. Companion-only groups (uranians) are generated and
# converted from their standalone partial via the Python tooling instead.
MAIN_FILE_LEB2_GROUPS = tuple(g for g in LEB2_GROUPS if g != "uranians")

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


def test_script_parses_cleanly() -> None:
    """`bash -n` accepts the whole file (syntax gate for the heredocs too)."""
    result = subprocess.run(
        ["/bin/bash", "-n", str(SCRIPT)],
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr


def test_shell_leb2_groups_match_canonical_python_registry() -> None:
    """The standalone Bash entrypoint converts exactly the merged-main groups."""
    source = SCRIPT.read_text(encoding="utf-8")
    match = re.search(r"^LEB2_GROUPS=\(([^)]*)\)$", source, flags=re.MULTILINE)

    assert match is not None
    assert tuple(match.group(1).split()) == MAIN_FILE_LEB2_GROUPS


def test_all_selects_only_canonical_leb2_groups() -> None:
    """``all`` expands in the canonical conversion/distribution order."""
    result = _run_sourced_script(
        "LEB2_GROUP_SELECTION=all; resolve_leb2_group_selection; "
        'printf "%s\\n" "${SELECTED_LEB2_GROUPS[@]}"'
    )

    assert result.returncode == 0, result.stderr
    assert tuple(result.stdout.splitlines()) == MAIN_FILE_LEB2_GROUPS


def test_retired_uranians_group_is_rejected() -> None:
    """A stale caller must fail before attempting to generate a bogus file."""
    result = _run_sourced_script(
        "parse_arguments --leb2-only --leb2-groups uranians; resolve_configuration"
    )

    assert result.returncode == 2
    assert "invalid LEB2 group: uranians" in result.stderr


def test_analytical_leb1_affects_only_surviving_leb2_groups() -> None:
    """Nodes and apsides map to core/apogee after hypothetical-body removal."""
    result = _run_sourced_script(
        "parse_arguments --group analytical; resolve_configuration; "
        'printf "%s\\n" "${SELECTED_LEB2_GROUPS[@]}"'
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == ["core", "apogee"]


def test_no_merge_without_leb1_only_is_rejected_up_front() -> None:
    """LEB2 always converts from the merged main: refuse the doomed combo."""
    result = _run_sourced_script("parse_arguments --no-merge; resolve_configuration")

    assert result.returncode == 2
    assert "--leb1-only" in result.stderr


def test_no_merge_with_leb1_only_is_accepted() -> None:
    result = _run_sourced_script(
        "parse_arguments --no-merge --leb1-only; resolve_configuration; echo OK"
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().endswith("OK")


def test_leb2_only_does_not_require_the_nbody_stack() -> None:
    """Conversion from an existing LEB1 never imports rebound/assist."""
    result = _run_sourced_script(
        "parse_arguments extended --leb2-only; resolve_configuration; "
        "if nbody_required; then echo NBODY=yes; else echo NBODY=no; fi"
    )

    assert result.returncode == 0, result.stderr
    assert "NBODY=no" in result.stdout


def test_extended_exotics_generation_requires_the_nbody_stack() -> None:
    result = _run_sourced_script(
        "parse_arguments extended --group exotics; resolve_configuration; "
        "if nbody_required; then echo NBODY=yes; else echo NBODY=no; fi"
    )

    assert result.returncode == 0, result.stderr
    assert "NBODY=yes" in result.stdout
