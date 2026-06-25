#!/usr/bin/env python3
"""Verify the close() semantics fixes from the 2026-06 PR review.

Checks (standalone, no pytest):

1. LEB configuration survives close():
   set_calc_mode("leb") + set_leb_file(path) must keep working after
   close() — the next calculation reopens the SAME configured file
   instead of raising "no .leb file configured" or silently switching
   to an auto-discovered file.

2. EphemerisContext.close() resets the state-level timescale cache:
   the next get_timescale() must build a FRESH timescale object (the
   documented resource reset), and the context- and module-level
   timescales must remain the same singleton.

Usage:
    python scripts/verify_fix_close_state.py
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

FAILURES: list[str] = []


def check(cond: bool, label: str) -> None:
    status = "ok  " if cond else "FAIL"
    print(f"  [{status}] {label}")
    if not cond:
        FAILURES.append(label)


def main() -> int:
    import libephemeris as le
    from libephemeris import context as ctx_module
    from libephemeris import state
    from libephemeris.context import EphemerisContext

    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    leb_path = os.path.join(repo, "libephemeris", "data", "leb2", "base_core.leb2")

    print("1) LEB config survives close() in mode 'leb'")
    if not os.path.exists(leb_path):
        print(f"  [skip] bundled LEB2 file missing: {leb_path}")
    else:
        le.set_calc_mode("leb")
        le.set_leb_file(leb_path)
        pos1, _ = le.calc_ut(2451545.0, le.SUN, 0)
        le.close()
        check(state._CALC_MODE == "leb", "calc mode preserved across close()")
        check(state._LEB_FILE == leb_path, "configured .leb path preserved")
        try:
            pos2, _ = le.calc_ut(2451545.0, le.SUN, 0)
        except Exception as exc:  # noqa: BLE001 — report, don't crash
            check(False, f"calc after close() raised: {exc}")
        else:
            check(
                abs(pos1[0] - pos2[0]) < 1e-12,
                "same longitude before/after close()",
            )
            check(
                state._LEB_FILE == leb_path,
                "same file reopened (no auto-discovery switch)",
            )
        # Leave a clean slate for the next section.
        le.set_calc_mode("auto")
        le.set_leb_file(None)
        le.close()

    print("2) EphemerisContext.close() resets the state timescale cache")
    ts_state_1 = state.get_timescale()
    ctx = EphemerisContext()
    ts_ctx_1 = ctx.get_timescale()
    check(ts_ctx_1 is ts_state_1, "context shares the module-level singleton")
    EphemerisContext.close()
    check(ctx_module._SHARED_TS is None, "_SHARED_TS cleared by close()")
    ts_state_2 = state.get_timescale()
    check(
        ts_state_2 is not ts_state_1,
        "state timescale rebuilt after context close()",
    )
    ts_ctx_2 = ctx.get_timescale()
    check(ts_ctx_2 is ts_state_2, "context re-caches the fresh singleton")

    print()
    if FAILURES:
        print(f"FAILED: {len(FAILURES)} check(s):")
        for failure in FAILURES:
            print(f"  - {failure}")
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
