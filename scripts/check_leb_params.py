#!/usr/bin/env python3
"""Check LEB file body parameters for Saturn.

Provenance:
    Project-authored format diagnostic. It reads and reports project-native LEB
    metadata through the registered reader; it neither changes the artifact nor
    derives a coefficient. Reported interval and degree values are generation
    parameters, not an independent scientific source.
"""

from __future__ import annotations
import sys

sys.path.insert(0, ".")
from libephemeris.leb_reader import LEBReader

reader = LEBReader("data/leb/ephemeris_base.leb")
print(f"Version: {reader.version}")
print(f"Bodies: {len(reader.bodies)}")
for body_id in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]:
    b = reader.bodies[body_id]
    print(
        f"  Body {body_id:2d}: interval={b.interval_days:.1f}d  degree={b.degree}  "
        f"coord_type={b.coord_type}  n_components={b.n_components}  "
        f"segments={b.n_segments}  jd=[{b.jd_start:.1f}, {b.jd_end:.1f}]"
    )
reader.close()
