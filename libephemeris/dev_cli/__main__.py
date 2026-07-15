# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Allow ``python -m libephemeris.dev_cli`` to launch the leph CLI.

Provenance:
    Project-authored entry-point plumbing with no astronomical model, data, or
    coefficient. It delegates immediately to the documented development CLI.
"""

from . import main

main()
