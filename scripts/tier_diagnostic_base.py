#!/usr/bin/env python3
"""
Tier diagnostic: base (de440s.bsp) -- SPK range 1850-2150.

Calculates all planets, lunar points, and minor bodies showing
ecliptic/equatorial coordinates, velocities, and data source
(DE440s / Analytical / SPK / Keplerian).

Usage:
    python scripts/tier_diagnostic_base.py
    python scripts/tier_diagnostic_base.py 1985-06-15
    python scripts/tier_diagnostic_base.py 2000-01-01 1950-03-15
    python scripts/tier_diagnostic_base.py --jd 2451545.0

Provenance:
    Project-authored entry point delegating to ``_tier_diagnostic.py`` and the
    registered base-tier DE440s/SPK/analytical channels. Dates are user-selected
    probes; printed coordinates and backend labels are diagnostics only and do
    not create scientific data or coefficients.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts._tier_diagnostic import run_diagnostic  # noqa: E402

if __name__ == "__main__":
    run_diagnostic("base")
