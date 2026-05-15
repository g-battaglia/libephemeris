"""Round 1: Deep hypothetical/Uranian body precision audit."""

from __future__ import annotations
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

BODIES = [
    (40, swe.CUPIDO, "Cupido"),
    (41, swe.HADES, "Hades"),
    (42, swe.ZEUS, "Zeus"),
    (43, swe.KRONOS, "Kronos"),
    (44, swe.APOLLON, "Apollon"),
    (45, swe.ADMETOS, "Admetos"),
    (46, swe.VULKANUS, "Vulkanus"),
    (47, swe.POSEIDON, "Poseidon"),
    (48, swe.ISIS, "Transpluto"),
]

DATES = [
    (1900, 1, 1, 12.0),
    (1950, 1, 1, 12.0),
    (1980, 6, 15, 0.0),
    (2000, 1, 1, 12.0),
    (2010, 7, 1, 12.0),
    (2024, 1, 1, 0.0),
    (2024, 6, 21, 12.0),
    (2030, 1, 1, 0.0),
    (2050, 6, 1, 0.0),
    (2100, 1, 1, 12.0),
]


def angular_diff(a: float, b: float) -> float:
    d = abs(a - b)
    return min(d, 360 - d)


print("=" * 100)
print("ROUND 1: HYPOTHETICAL BODIES DEEP PRECISION AUDIT")
print("=" * 100)

max_diffs: dict[str, float] = {}
all_results = []

for body_le, body_se, name in BODIES:
    print(f"\n--- {name} (ID {body_le}) ---")
    for y, m, d, h in DATES:
        jd = swe.julday(y, m, d, h)
        try:
            pos_se, _ = swe.calc_ut(jd, body_se, swe.FLG_SPEED)
        except Exception as e:
            print(f"  SE skip {y}: {e}")
            continue
        try:
            pos_le, _ = ephem.calc_ut(jd, body_le, 256)  # FLG_SPEED
        except Exception as e:
            print(f"  LE skip {y}: {e}")
            continue

        dl = angular_diff(pos_se[0], pos_le[0])
        db = abs(pos_se[1] - pos_le[1])
        dd = abs(pos_se[2] - pos_le[2])
        dsl = abs(pos_se[3] - pos_le[3])

        all_results.append((name, y, dl, db, dd, dsl))
        if name not in max_diffs or dl > max_diffs[name]:
            max_diffs[name] = dl

        flag = "!!!" if dl > 0.5 else ("**" if dl > 0.1 else "")
        print(
            f'  {y}-{m:02d}-{d:02d}: dLon={dl:.6f}° ({dl * 3600:.1f}")  '
            f"dLat={db:.6f}°  dDist={dd:.6f} AU  dSpd={dsl:.6f}°/d {flag}"
        )

print("\n" + "=" * 100)
print("SUMMARY — MAX LONGITUDE DIFFERENCES:")
for name, mx in sorted(max_diffs.items(), key=lambda x: -x[1]):
    flag = (
        " *** NEEDS INVESTIGATION"
        if mx > 0.5
        else (" ** NOTABLE" if mx > 0.1 else " OK")
    )
    print(f'  {name:15s}: {mx:.6f}° ({mx * 3600:.1f}"){flag}')
