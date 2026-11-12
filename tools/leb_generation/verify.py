# SPDX-License-Identifier: AGPL-3.0
"""Verify actual native LEB reads against the selected native source models."""

from pathlib import Path

import numpy as np

from .checks import check_positions
from .read import read_body_index
from .source import ANGULAR_IDS, check_input_hashes


def verify_source_file(path: Path, source, samples=64) -> dict[str, int]:
    """Check fixed PLAN position/rate criteria through the real native reader.

    Includes declared endpoints and deterministic off-grid interior epochs for
    every body. Cartesian source velocities are analytic. Angular channels
    receive position checks only; full apparent-frame/rate qualification is a
    separate public-pipeline comparison, not implied by these fit checks.
    """
    import libephemeris as engine
    if samples < 2:
        raise ValueError("verification needs at least two samples per body")
    _, entries = read_body_index(path)
    if not entries:
        raise ValueError("source verification requires a nonempty body inventory")
    previous = engine.get_calc_mode()
    try:
        engine.set_calc_mode("leb")
        context = engine.EphemerisContext(str(Path(path).resolve().parent))
        context.set_leb_file(str(Path(path).resolve()))
    finally:
        engine.set_calc_mode(previous)
    reader = context.get_leb_reader()
    counts = {"positions": 0, "cartesian_rates": 0}
    try:
        for entry in entries:
            source.require_domain(entry.body_id, entry.jd_start, entry.jd_end)
            epochs = np.linspace(entry.jd_start, entry.jd_end, samples)
            values = [reader.eval_body(entry.body_id, float(epoch)) for epoch in epochs]
            actual = np.asarray([value[0] for value in values])
            angular = entry.body_id in ANGULAR_IDS
            expected_states = None if angular else source.states(entry.body_id, epochs)
            expected = source.values(entry.body_id, epochs) if angular else expected_states[:, :3]
            check_positions(expected, actual, angular, entry.body_id)
            counts["positions"] += samples
            if not angular:
                rates = np.asarray([value[1] for value in values])
                expected_rates = expected_states[:, 3:]
                if not np.isfinite(rates).all() or not np.isfinite(expected_rates).all():
                    raise ValueError(f"body {entry.body_id}: nonfinite Cartesian rate")
                error = np.linalg.norm(rates - expected_rates, axis=1)
                if np.any(error > 1e-6 * np.linalg.norm(expected_rates, axis=1)):
                    raise ValueError(f"body {entry.body_id}: native-read Cartesian rate exceeds PLAN criterion")
                counts["cartesian_rates"] += samples
    finally:
        reader.close()
    # This validator also runs immediately before either writer publishes its
    # temporary file. Serial fitting must reject changed clock/image inputs just
    # as spawned fitting does, including newly appearing optional clock files.
    check_input_hashes(source.input_hashes().items())
    return counts
