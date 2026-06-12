"""Lunar calibration commands: perigee perturbation coefficient fitting.

Replaces 2 poe tasks: calibrate-perigee, calibrate-perigee:quick.
"""

from __future__ import annotations

import subprocess
import sys

import click


def _python(args: list[str], use_dotenv: bool = True) -> None:
    """Run a python script, optionally loading .env."""
    import os

    run_env = dict(os.environ)
    if use_dotenv:
        # Load .env file if it exists (same behavior as poe's envfile = ".env")
        env_path = os.path.join(os.getcwd(), ".env")
        if os.path.exists(env_path):
            with open(env_path) as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#") and "=" in line:
                        key, _, value = line.partition("=")
                        run_env[key.strip()] = value.strip().strip("\"'")
    sys.exit(subprocess.call([sys.executable, *args], env=run_env))


@click.group(
    "calibrate",
    short_help="Fit lunar perigee perturbation coefficients against JPL DE441.",
    help="Calibrate lunar perigee perturbation coefficients against JPL DE441.\n\n"
    "Fits ELP2000-based harmonic perturbation coefficients to minimize the\n"
    "difference between analytical and geometric (JPL) perigee positions.\n\n"
    "Full workflow:\n\n"
    "  1. leph calibrate perigee              # Fit coefficients (~30 min)\n"
    "  2. Paste output into lunar.py           # _calc_elp2000_perigee_perturbations()\n"
    "  3. leph generate apse-corrections       # Regenerate residual tables\n"
    "  4. leph test lunar perigee              # Verify accuracy",
)
def calibrate_group() -> None:
    """Calibration commands."""


@calibrate_group.command(
    short_help="Full perigee calibration (1500-2500 CE, ~30 min).",
)
def perigee() -> None:
    """Calibrate perigee perturbation coefficients against JPL DE441 (~30 min).

    Full run: 1500-2500 CE range, passage-interpolated harmonic fit method.
    Output: /tmp/perigee_calibration.json

    Workflow:
      1. leph calibrate perigee
      2. Paste coefficients into _calc_elp2000_perigee_perturbations() in lunar.py
      3. leph generate apse-corrections
      4. leph test lunar perigee
    """
    _python(
        [
            "scripts/calibrate_perigee_perturbations.py",
            "--start-year",
            "1500",
            "--end-year",
            "2500",
            "--output",
            "/tmp/perigee_calibration.json",
        ]
    )


@calibrate_group.command(
    "perigee-quick",
    short_help="Quick perigee calibration (100-year range, ~2 min).",
)
def perigee_quick() -> None:
    """Quick perigee calibration (100-year range, ~2 min).

    For validation only -- use full 'perigee' for production coefficients.
    Output: /tmp/perigee_calibration_quick.json
    """
    _python(
        [
            "scripts/calibrate_perigee_perturbations.py",
            "--quick",
            "--output",
            "/tmp/perigee_calibration_quick.json",
        ]
    )
