# SPDX-License-Identifier: AGPL-3.0
"""Fixed PLAN position criteria for offline source fitting and native reads."""

import numpy as np


def check_positions(expected, actual, angular, body_id):
    """Check geometric vectors or wrapped angular channels in their native units.

    Cartesian vector error is relative to the source vector norm. Angular
    longitude/latitude use 0.001 arcsecond and range uses relative 1e-8. These
    raw-fit checks do not replace complete apparent-coordinate validation.
    """
    if not np.isfinite(expected).all() or not np.isfinite(actual).all():
        raise ValueError(f"body {body_id}: nonfinite verification position")
    difference = actual - expected
    if angular:
        difference[:, 0] = (difference[:, 0] + 180.0) % 360.0 - 180.0
        failed = (np.any(np.abs(difference[:, :2]) * 3600.0 > 0.001)
                  or np.any(np.abs(difference[:, 2]) > 1e-8 * np.abs(expected[:, 2])))
    else:
        failed = np.any(np.linalg.norm(difference, axis=1) > 1e-8 * np.linalg.norm(expected, axis=1))
    if failed:
        raise ValueError(f"body {body_id}: fitted position exceeds PLAN criterion")
