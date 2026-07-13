"""Regression tests for the offline lunar-model provenance gate."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def test_lunar_model_provenance_gate_has_zero_mismatches() -> None:
    root = Path(__file__).resolve().parents[1]
    result = subprocess.run(
        [sys.executable, "scripts/check_lunar_model_provenance.py"],
        cwd=root,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "0 mismatch(es)" in result.stdout
