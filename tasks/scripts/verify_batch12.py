"""Fast standalone verification for Batch 12 work.

Covers:
- BUG-006 fix: Uranian bodies via Horizons analytical path
- 360° wrap-around correctness
- ECL_NUT special body
- azalt/azalt_rev/refrac round-trips
- TAI time function round-trips
- Horizons analytical bodies (Mean Node, Mean Apogee)
- calc_mode switching
"""

from __future__ import annotations

import math
import sys
import os

# Ensure we use the project
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    MARS,
    JUPITER,
    MEAN_NODE,
    MEAN_APOG,
    ECL_NUT,
    FLG_SWIEPH,
    FLG_SPEED,
    FLG_HELCTR,
    FLG_EQUATORIAL,
    FLG_SIDEREAL,
    FLG_TOPOCTR,
)
from libephemeris.utils import (
    ECL2HOR,
    EQU2HOR,
    HOR2EQU,
    TRUE_TO_APP,
    APP_TO_TRUE,
)
from libephemeris.horizons_backend import horizons_calc_ut

JD = 2451545.0  # J2000
passed = 0
failed = 0
errors = []


def check(name: str, condition: bool, detail: str = ""):
    global passed, failed
    if condition:
        passed += 1
    else:
        failed += 1
        msg = f"FAIL: {name}" + (f" — {detail}" if detail else "")
        errors.append(msg)
        print(msg)


# ── BUG-006: Uranian bodies via Horizons analytical path ──
print("=== BUG-006: Uranian bodies via Horizons analytical ===")
for body_id, name in [
    (40, "Cupido"),
    (41, "Hades"),
    (42, "Zeus"),
    (43, "Kronos"),
    (44, "Apollon"),
    (45, "Admetos"),
    (46, "Vulkanus"),
    (47, "Poseidon"),
]:
    try:
        result = horizons_calc_ut(None, JD, body_id, FLG_SWIEPH | FLG_HELCTR)
        data = result[0]
        check(f"{name} returns 6 values", len(data) == 6, f"got {len(data)}")
        check(f"{name} lon in [0,360)", 0.0 <= data[0] < 360.0, f"lon={data[0]}")
        check(f"{name} dist > 10 AU", data[2] > 10.0, f"dist={data[2]}")
    except Exception as e:
        check(f"{name} no error", False, str(e))

# Verify Uranian geocentric raises (expected)
try:
    horizons_calc_ut(None, JD, 40, FLG_SWIEPH)
    check("Uranian geocentric raises KeyError", False, "did not raise")
except KeyError:
    check("Uranian geocentric raises KeyError", True)

# TOPOCTR raises
try:
    horizons_calc_ut(None, JD, 0, FLG_TOPOCTR)
    check("TOPOCTR raises KeyError", False, "did not raise")
except KeyError:
    check("TOPOCTR raises KeyError", True)

# ── Horizons analytical: Mean Node, Mean Apogee ──
print("\n=== Horizons analytical bodies ===")
for body_id, name in [(10, "Mean Node"), (12, "Mean Apogee")]:
    result = horizons_calc_ut(None, JD, body_id, FLG_SWIEPH | FLG_SPEED)
    data = result[0]
    ref, _ = swe.calc_ut(JD, body_id, FLG_SWIEPH | FLG_SPEED)
    diff = abs(data[0] - ref[0])
    check(f"{name} horizons vs skyfield < 0.01°", diff < 0.01, f"diff={diff:.6f}°")
    check(f"{name} lon in [0,360)", 0.0 <= data[0] < 360.0)

# Mean Node speed should be negative (retrograde)
result = horizons_calc_ut(None, JD, 10, FLG_SWIEPH | FLG_SPEED)
check("Mean Node speed negative", result[0][3] < 0, f"speed={result[0][3]}")

# ── 360° wrap-around ──
print("\n=== 360° wrap-around ===")
for body in [SUN, MOON, MARS, JUPITER, MEAN_NODE]:
    for i in range(20):
        jd = JD + i * 18.25
        r, _ = swe.calc_ut(jd, body, FLG_SWIEPH)
        check(
            f"body {body} JD+{i * 18.25:.0f} in [0,360)",
            0.0 <= r[0] < 360.0,
            f"lon={r[0]}",
        )

# solcross near 0°
jd_cross = swe.solcross_ut(0.0, JD, 0)
r, _ = swe.calc_ut(jd_cross, SUN, FLG_SWIEPH)
check("solcross 0° Sun near 0°", r[0] < 0.01 or r[0] > 359.99, f"lon={r[0]}")

# Sidereal range
swe.set_sid_mode(1)
for body in [SUN, MOON, MARS]:
    r, _ = swe.calc_ut(JD, body, FLG_SWIEPH | FLG_SIDEREAL)
    check(f"sidereal body {body} in [0,360)", 0.0 <= r[0] < 360.0, f"lon={r[0]}")
swe.set_sid_mode(0)

# ── ECL_NUT ──
print("\n=== ECL_NUT special body ===")
r, _ = swe.calc_ut(JD, ECL_NUT, FLG_SWIEPH)
check("ECL_NUT 6 values", len(r) == 6)
check("ECL_NUT true obl ~23.4°", 23.0 < r[0] < 24.0, f"true_obl={r[0]}")
check("ECL_NUT mean obl ~23.4°", 23.0 < r[1] < 24.0, f"mean_obl={r[1]}")
check("ECL_NUT nut_lon small", abs(r[2]) < 0.01, f"nut_lon={r[2]}")
check("ECL_NUT nut_obl small", abs(r[3]) < 0.01, f"nut_obl={r[3]}")
check("ECL_NUT [4]=0", r[4] == 0.0)
check("ECL_NUT [5]=0", r[5] == 0.0)

# ── azalt / azalt_rev / refrac ──
print("\n=== azalt / azalt_rev / refrac ===")
sun, _ = swe.calc_ut(JD, SUN, FLG_SWIEPH)
az, true_alt, app_alt = swe.azalt(
    JD, ECL2HOR, (12.5, 41.9, 50.0), 1013.25, 15.0, (sun[0], sun[1], sun[2])
)
check("azalt az in [0,360)", 0.0 <= az < 360.0, f"az={az}")
check("azalt true_alt in [-90,90]", -90 <= true_alt <= 90, f"alt={true_alt}")
check("azalt all finite", all(math.isfinite(v) for v in [az, true_alt, app_alt]))

# azalt_rev round-trip
sun_eq, _ = swe.calc_ut(JD, SUN, FLG_SWIEPH | FLG_EQUATORIAL)
az2, ta2, _ = swe.azalt(
    JD, EQU2HOR, (12.5, 41.9, 50.0), 1013.25, 15.0, (sun_eq[0], sun_eq[1], sun_eq[2])
)
ra_out, dec_out = swe.azalt_rev(JD, 1, (12.5, 41.9, 50.0), az2, ta2)  # HOR2EQU=1
ra_diff = abs(ra_out - sun_eq[0])
if ra_diff > 180:
    ra_diff = 360 - ra_diff
check("azalt round-trip RA < 0.05°", ra_diff < 0.05, f"diff={ra_diff}")
check("azalt round-trip Dec < 0.05°", abs(dec_out - sun_eq[1]) < 0.05)

# refrac round-trip
for alt in [0.0, 5.0, 15.0, 30.0, 60.0, 85.0]:
    app = swe.refrac(alt, 1013.25, 15.0, TRUE_TO_APP)
    recovered = swe.refrac(app, 1013.25, 15.0, APP_TO_TRUE)
    check(
        f"refrac round-trip alt={alt}°",
        abs(recovered - alt) < 0.02,
        f"err={abs(recovered - alt)}",
    )

# Horizon refraction ~34'
ref_horizon = swe.refrac(0.0, 1013.25, 15.0, TRUE_TO_APP)
check("horizon refraction ~0.57°", 0.4 < ref_horizon < 0.7, f"ref={ref_horizon}")

# ── TAI time functions ──
print("\n=== TAI time functions ===")
# UTC -> TAI -> UTC round-trip
for y, m, d, h, mi, s in [
    (2000, 1, 1, 12, 0, 0.0),
    (2020, 6, 15, 6, 30, 45.5),
    (2010, 3, 1, 0, 0, 0.0),
]:
    jd_tai = swe.utc_to_tai_jd(y, m, d, h, mi, s)
    y2, m2, d2, h2, mi2, s2 = swe.tai_jd_to_utc(jd_tai)
    check(f"UTC-TAI-UTC {y}-{m}-{d} date match", y2 == y and m2 == m and d2 == d)
    check(f"UTC-TAI-UTC {y}-{m}-{d} sec match", abs(s2 - s) < 0.001, f"sec={s2}")

# TT -> TAI -> TT round-trip
for jd_tt in [JD, JD + 3652.5, JD + 7305.0]:
    jd_tai = swe.tt_to_tai_jd(jd_tt)
    jd_tt_out = swe.tai_to_tt_jd(jd_tai)
    check(f"TT-TAI-TT JD={jd_tt}", abs(jd_tt_out - jd_tt) < 1e-12)

# TT-TAI offset ~32.184s
jd_tai = swe.tt_to_tai_jd(JD)
diff_s = (JD - jd_tai) * 86400.0
check("TT-TAI = 32.184s", abs(diff_s - 32.184) < 0.001, f"diff={diff_s}")

# TAI-UTC at J2000 = 32s
leap = swe.get_tai_utc_for_jd(JD)
check("TAI-UTC at J2000 = 32s", abs(leap - 32.0) < 0.5, f"leap={leap}")

# UTC JD round-trip
for jd_utc in [JD, JD + 3652.5, JD + 7305.0]:
    jd_tai = swe.utc_jd_to_tai(jd_utc)
    jd_utc_out = swe.tai_to_utc_jd(jd_tai)
    err_us = abs(jd_utc_out - jd_utc) * 86400e6
    check(f"UTC-TAI-UTC JD {jd_utc}", err_us < 1.0, f"err={err_us:.3f}µs")

# ── calc_mode switching ──
print("\n=== calc_mode switching ===")
orig = swe.get_calc_mode()
swe.set_calc_mode("skyfield")
check("set_calc_mode skyfield", swe.get_calc_mode() == "skyfield")
swe.set_calc_mode("auto")
check("set_calc_mode auto", swe.get_calc_mode() == "auto")
try:
    swe.set_calc_mode("invalid_xyz")
    check("invalid mode raises", False)
except ValueError:
    check("invalid mode raises", True)
swe.set_calc_mode(orig)

# ── Summary ──
print(f"\n{'=' * 50}")
print(f"TOTAL: {passed} passed, {failed} failed out of {passed + failed}")
if errors:
    print("\nFailures:")
    for e in errors:
        print(f"  {e}")
    sys.exit(1)
else:
    print("ALL CHECKS PASSED ✓")
