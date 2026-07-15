# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Generate the independent interpolated lunar-apse model from JPL DE440.

The interpolated (natural) apogee/perigee is defined here, independently and
reproducibly, as the smooth track through the Moon's actual apsis passages:

1. Find every geocentric-distance extremum of the DE440 Moon over the medium
   ephemeris range (~1549--2651 CE), refined to ~1 second.
2. Record the Moon's geometric true-ecliptic-of-date longitude, latitude and
   distance at each passage.
3. Subtract the mean apse derived from the IERS 2003 Delaunay arguments and
   fit a trigonometric series in the same IERS arguments at the passages.
4. Interpolate the remaining passage residual with a natural cubic spline and
   sample it onto the fixed runtime correction grid.
5. Refit the 2-harmonic latitude model and mean apsis distances.

Everything derives from the JPL DE440 kernel (Park et al. 2021, AJ 161:105),
the IAU SOFA/IERS 2003 fundamental arguments exposed by BSD-licensed ERFA,
the mean Earth-orbit eccentricity published by Simon et al. (1994), and this
script. No value is fitted to any external ephemeris implementation.

Outputs ``libephemeris/lunar_apse_model.py`` (tables + fitted coefficients).

Run:  uv run python scripts/generate_lunar_apse_model.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

J2000 = 2451545.0

# Correction-table grids (unchanged layout from the previous model).
APOGEE_CT_JD_START = 2286820.5
APOGEE_CT_STEP_DAYS = 10.0
APOGEE_CT_COUNT = 40250
PERIGEE_CT_JD_START = 2286820.5
PERIGEE_CT_STEP_DAYS = 2.0
PERIGEE_CT_COUNT = 201249
EDGE_TAPER_DAYS = 365.25

# --------------------------------------------------------------------------
# Delaunay trigonometric basis (project design):
# each term is (cD, cM, cMp, cF, e_power) -> E**e_power * sin(cD*D + cM*M +
# cMp*Mp + cF*F). The evection ladder k(D - Mp), its first and second solar
# couplings, and the 2(F - Mp) inclination term.
# --------------------------------------------------------------------------


def _basis(max_k: int, m_ks: list[int]) -> list[tuple[int, int, int, int, int]]:
    terms: list[tuple[int, int, int, int, int]] = []
    for k in range(1, max_k + 1):
        terms.append((k, 0, -k, 0, 0))
    terms.append((0, 1, 0, 0, 1))
    terms.append((1, 1, -1, 0, 1))
    terms.append((1, -1, -1, 0, 1))
    for k in m_ks:
        terms.append((k, -1, -k, 0, 1))
        terms.append((k, 1, -k, 0, 1))
    terms.append((0, 2, 0, 0, 2))
    terms.append((2, -2, -2, 0, 2))
    terms.append((2, 2, -2, 0, 2))
    terms.append((4, -2, -4, 0, 2))
    terms.append((4, 2, -4, 0, 2))
    terms.append((6, -2, -6, 0, 2))
    terms.append((0, 0, -2, 2, 0))
    return terms


APOGEE_TERMS = _basis(14, [2, 4, 6, 8, 10])
PERIGEE_TERMS = _basis(30, list(range(2, 15)))


def _fundamental_arguments(jd_tt: np.ndarray) -> tuple[np.ndarray, ...]:
    """Return the IERS 2003 Delaunay arguments in radians."""
    from libephemeris.mean_lunar_apse import lunar_delaunay_arguments

    D = np.empty_like(jd_tt)
    M = np.empty_like(jd_tt)
    Mp = np.empty_like(jd_tt)
    F = np.empty_like(jd_tt)
    for i, jd in enumerate(jd_tt):
        D[i], M[i], Mp[i], F[i] = lunar_delaunay_arguments(float(jd))
    return D, M, Mp, F


def _design_matrix(
    jd_tt: np.ndarray, terms: list[tuple[int, int, int, int, int]]
) -> np.ndarray:
    from libephemeris.mean_lunar_apse import (
        EARTH_ECCENTRICITY_J2000,
        EARTH_ECCENTRICITY_T,
        EARTH_ECCENTRICITY_T2,
    )

    D, M, Mp, F = _fundamental_arguments(jd_tt)
    T = (jd_tt - J2000) / 36525.0
    eccentricity = (
        EARTH_ECCENTRICITY_J2000
        + EARTH_ECCENTRICITY_T * T
        + EARTH_ECCENTRICITY_T2 * T**2
    )
    E = eccentricity / EARTH_ECCENTRICITY_J2000
    cols = []
    for cd, cm, cmp_, cf, epow in terms:
        arg = cd * D + cm * M + cmp_ * Mp + cf * F
        cols.append((E**epow) * np.sin(arg))
    return np.column_stack(cols)


# --------------------------------------------------------------------------
# DE440 apsis passages
# --------------------------------------------------------------------------


def _moon_states(jd_tt: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Geometric geocentric Moon in the true ecliptic of date (lon, lat, dist)."""
    from libephemeris.cache import get_cached_nutation
    from libephemeris.precession_vondrak import vondrak_pn_matrix
    from libephemeris.state import get_planets, get_timescale

    planets = get_planets()
    ts = get_timescale()
    t = ts.tt_jd(jd_tt)
    moon = planets["moon"]
    earth = planets["earth"]
    icrs = (moon - earth).at(t).position.au  # geometric vector, shape (3, n)

    lon = np.empty(jd_tt.shape)
    lat = np.empty(jd_tt.shape)
    dist = np.empty(jd_tt.shape)
    for i, jd in enumerate(jd_tt):
        dpsi, deps = get_cached_nutation(float(jd))
        pn, eps_true = vondrak_pn_matrix(float(jd), float(dpsi), float(deps))
        v = np.asarray(pn) @ icrs[:, i]
        ce, se = math.cos(eps_true), math.sin(eps_true)
        xe = v[0]
        ye = v[1] * ce + v[2] * se
        ze = -v[1] * se + v[2] * ce
        r = math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
        lon[i] = math.degrees(math.atan2(ye, xe)) % 360.0
        lat[i] = math.degrees(math.asin(ze / r))
        dist[i] = r
    return lon, lat, dist


def _moon_distance(jd_tt: np.ndarray) -> np.ndarray:
    from libephemeris.state import get_planets, get_timescale

    planets = get_planets()
    ts = get_timescale()
    t = ts.tt_jd(jd_tt)
    icrs = (planets["moon"] - planets["earth"]).at(t).position.au
    return np.sqrt((icrs**2).sum(axis=0))


def find_passages(jd_start: float, jd_end: float) -> tuple[np.ndarray, np.ndarray]:
    """Refined apogee and perigee passage times (TT JD) over [jd_start, jd_end]."""
    step = 0.5
    grid = np.arange(jd_start, jd_end + step, step)
    r = _moon_distance(grid)

    inner = np.arange(1, len(grid) - 1)
    is_max = (r[inner] >= r[inner - 1]) & (r[inner] > r[inner + 1])
    is_min = (r[inner] <= r[inner - 1]) & (r[inner] < r[inner + 1])

    def _refine(idx: np.ndarray) -> np.ndarray:
        t0 = grid[idx]
        h = step
        for _ in range(6):
            tm, tp = t0 - h, t0 + h
            rm = _moon_distance(tm)
            r0 = _moon_distance(t0)
            rp = _moon_distance(tp)
            denom = rm - 2.0 * r0 + rp
            shift = np.where(np.abs(denom) > 0.0, 0.5 * h * (rm - rp) / denom, 0.0)
            shift = np.clip(shift, -h, h)
            t0 = t0 + shift
            h = h / 4.0
        return t0

    apogees = _refine(inner[is_max])
    perigees = _refine(inner[is_min])
    return apogees, perigees


# --------------------------------------------------------------------------
# Natural cubic spline (vectorized Thomas solve)
# --------------------------------------------------------------------------


class NaturalCubicSpline:
    def __init__(self, x: np.ndarray, y: np.ndarray) -> None:
        n = len(x) - 1
        h = np.diff(x)
        alpha = np.zeros(n + 1)
        alpha[1:n] = 3.0 * ((y[2:] - y[1:-1]) / h[1:] - (y[1:-1] - y[:-2]) / h[:-1])
        line = np.ones(n + 1)
        mu = np.zeros(n + 1)
        z = np.zeros(n + 1)
        for i in range(1, n):
            line[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1]
            mu[i] = h[i] / line[i]
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / line[i]
        c = np.zeros(n + 1)
        b = np.zeros(n)
        d = np.zeros(n)
        for j in range(n - 1, -1, -1):
            c[j] = z[j] - mu[j] * c[j + 1]
            b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0
            d[j] = (c[j + 1] - c[j]) / (3.0 * h[j])
        self.x, self.y, self.b, self.c, self.d = x, y[:-1], b, c[:-1], d

    def __call__(self, xq: np.ndarray) -> np.ndarray:
        idx = np.clip(np.searchsorted(self.x, xq) - 1, 0, len(self.b) - 1)
        dx = xq - self.x[idx]
        return self.y[idx] + dx * (self.b[idx] + dx * (self.c[idx] + dx * self.d[idx]))


# --------------------------------------------------------------------------
# Model generation
# --------------------------------------------------------------------------


def build_side(
    passages: np.ndarray,
    mean_fn,
    terms: list[tuple[int, int, int, int, int]],
    ct_start: float,
    ct_step: float,
    ct_count: int,
    label: str,
):
    """Fit one apsis side; returns (coeffs, corrections, lat_fit, mean_dist, stats).

    The trigonometric series is fitted **directly at the DE440 passages** —
    the Delaunay harmonics carry the physically meaningful structure between
    passages, where an interpolant through 27.55-day knots could not. Only
    the small, slowly varying passage residuals are spline-interpolated onto
    the correction-table grid, so the model is exact (to the residual-spline
    error) at every actual apsis passage.
    """
    lon_p, lat_p, dist_p = _moon_states(passages)

    mean_lon = np.array([mean_fn(float(jd)) for jd in passages])
    delta = lon_p - mean_lon
    delta -= 360.0 * np.round(delta / 360.0)

    A = _design_matrix(passages, terms)
    coeffs, *_ = np.linalg.lstsq(A, delta, rcond=None)
    residual = delta - A @ coeffs

    # Interpolate the passage residuals onto the correction-table grid. The
    # grid has roughly one year of margin outside DE440 coverage; taper the
    # first/last observed residual through that margin so interpolation at the
    # terminal passages remains continuous and accurate.
    resid_spline = NaturalCubicSpline(passages, residual)
    grid = ct_start + ct_step * np.arange(ct_count)
    lo, hi = passages[0], passages[-1]
    mask = (grid >= lo) & (grid <= hi)
    corrections = np.zeros(ct_count)
    corrections[mask] = resid_spline(grid[mask]) * 3600.0
    before = grid < lo
    before_weight = np.clip(1.0 - (lo - grid[before]) / EDGE_TAPER_DAYS, 0.0, 1.0)
    corrections[before] = residual[0] * before_weight * 3600.0
    after = grid > hi
    after_weight = np.clip(1.0 - (grid[after] - hi) / EDGE_TAPER_DAYS, 0.0, 1.0)
    corrections[after] = residual[-1] * after_weight * 3600.0
    interpolated_residual = np.interp(passages, grid, corrections) / 3600.0
    table_error = residual - interpolated_residual

    # Latitude 2-harmonic refit at the passages: lat = a*sin(w) + b*sin(3w),
    # w = passage longitude - mean node longitude.
    from libephemeris.mean_lunar_apse import mean_lunar_node_position

    node = np.array([mean_lunar_node_position(float(jd))[0] for jd in passages])
    w = np.radians((lon_p % 360.0) - node)
    L = np.column_stack([np.sin(w), np.sin(3.0 * w)])
    lat_fit, *_ = np.linalg.lstsq(L, lat_p, rcond=None)
    lat_resid = lat_p - L @ lat_fit

    mean_dist = float(np.mean(dist_p))

    stats = {
        "passages": len(passages),
        "trig_rms_arcsec": float(np.sqrt(np.mean(residual**2))) * 3600.0,
        "trig_max_arcsec": float(np.max(np.abs(residual))) * 3600.0,
        "table_rms_arcsec": float(np.sqrt(np.mean(table_error**2))) * 3600.0,
        "table_max_arcsec": float(np.max(np.abs(table_error))) * 3600.0,
        "pre_trig_rms_deg": float(np.sqrt(np.mean(delta**2))),
        "lat_rms_arcsec": float(np.sqrt(np.mean(lat_resid**2))) * 3600.0,
        "mean_dist_au": mean_dist,
    }
    print(f"[{label}] {stats}")
    return coeffs, corrections, lat_fit, mean_dist, stats


def main() -> int:
    from libephemeris.mean_lunar_apse import (
        mean_lunar_apogee_position,
    )

    def mean_apogee(jd: float) -> float:
        return mean_lunar_apogee_position(jd)[0]

    def mean_perigee(jd: float) -> float:
        return (mean_lunar_apogee_position(jd)[0] + 180.0) % 360.0

    # Clamp to the DE440 medium kernel coverage (with a refinement margin).
    # The wider correction grids provide a one-year taper around the passage
    # interval without extrapolating the fitted spline.
    from libephemeris.state import get_current_file_data, get_planets

    spk = get_planets()["moon"].ephemeris.spk
    _, _, _, denum = get_current_file_data(0)
    if denum != 440:
        raise RuntimeError(
            f"lunar apse model generation requires DE440, active kernel is DE{denum}"
        )
    moon_segs = [s for s in spk.segments if s.target == 301]
    kernel_lo = max(s.start_jd for s in moon_segs)
    kernel_hi = min(s.end_jd for s in moon_segs)
    jd_lo = max(APOGEE_CT_JD_START + 5.0, kernel_lo + 3.0)
    jd_hi = min(
        APOGEE_CT_JD_START + APOGEE_CT_STEP_DAYS * (APOGEE_CT_COUNT - 1) - 5.0,
        kernel_hi - 3.0,
    )
    print(f"searching passages over JD [{jd_lo}, {jd_hi}] ...")
    apogees, perigees = find_passages(jd_lo, jd_hi)
    print(f"found {len(apogees)} apogee / {len(perigees)} perigee passages")

    apo = build_side(
        apogees,
        mean_apogee,
        APOGEE_TERMS,
        APOGEE_CT_JD_START,
        APOGEE_CT_STEP_DAYS,
        APOGEE_CT_COUNT,
        "apogee",
    )
    per = build_side(
        perigees,
        mean_perigee,
        PERIGEE_TERMS,
        PERIGEE_CT_JD_START,
        PERIGEE_CT_STEP_DAYS,
        PERIGEE_CT_COUNT,
        "perigee",
    )

    _write_module(apo, per)
    return 0


def _fmt_terms(terms, coeffs) -> str:
    rows = []
    for (cd, cm, cmp_, cf, epow), c in zip(terms, coeffs):
        rows.append(f"    ({float(c)!r}, {cd}, {cm}, {cmp_}, {cf}, {epow}),\n")
    return "".join(rows)


def _fmt_table(values: np.ndarray, per_line: int = 8) -> str:
    out = []
    line: list[str] = []
    for v in values:
        line.append(f"{v:.3f}")
        if len(line) == per_line:
            out.append("    " + ", ".join(line) + ",\n")
            line = []
    if line:
        out.append("    " + ", ".join(line) + ",\n")
    return "".join(out)


def _write_module(apo, per) -> None:
    apo_coeffs, apo_corr, apo_lat, apo_dist, apo_stats = apo
    per_coeffs, per_corr, per_lat, per_dist, per_stats = per
    header = f'''# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Independent interpolated lunar-apse model, derived from JPL DE440.

AUTO-GENERATED by scripts/generate_lunar_apse_model.py - DO NOT EDIT.
The output is deterministic for the reviewed DE440 kernel and generator.

Derivation (fully documented in the generator and in
docs/methodology/interpolated-perigee.md):
  1. Every apsis passage of the DE440 Moon (Park et al. 2021, AJ 161:105)
     over ~1549--2651 CE, refined to ~1 second.
  2. Geometric true-ecliptic-of-date longitude at each passage.
  3. Delaunay-argument trigonometric series fitted by least squares against
     track - the IERS 2003 mean apse.
  4. Natural-spline interpolation of the remaining passage residual, sampled
     onto fixed runtime grids in arcseconds.
  At runtime: longitude = mean + trig series + interpolated correction.

Sources: JPL DE440; IERS Conventions (2010), Eq. 5.43; ERFA's BSD-licensed
implementation of the IAU SOFA fundamental-argument routines; and the mean
Earth-orbit eccentricity from Simon et al. (1994), A&A 282, 663--683.
No value derives from any external ephemeris implementation.

Fit quality (vs the DE440 passage track):
  apogee : trig RMS {apo_stats["trig_rms_arcsec"]:.1f}\", table RMS {apo_stats["table_rms_arcsec"]:.2f}\",
           table max {apo_stats["table_max_arcsec"]:.2f}\",
           latitude RMS {apo_stats["lat_rms_arcsec"]:.1f}\", {apo_stats["passages"]} passages
  perigee: trig RMS {per_stats["trig_rms_arcsec"]:.1f}\", table RMS {per_stats["table_rms_arcsec"]:.2f}\",
           table max {per_stats["table_max_arcsec"]:.2f}\",
           latitude RMS {per_stats["lat_rms_arcsec"]:.1f}\", {per_stats["passages"]} passages
"""

from __future__ import annotations

# Trigonometric terms: (coefficient_deg, cD, cM, cMp, cF, e_power) ->
# coefficient * E**e_power * sin(cD*D + cM*M + cMp*Mp + cF*F).
APOGEE_TRIG_TERMS: tuple[tuple[float, int, int, int, int, int], ...] = (
{_fmt_terms(APOGEE_TERMS, apo_coeffs)})

PERIGEE_TRIG_TERMS: tuple[tuple[float, int, int, int, int, int], ...] = (
{_fmt_terms(PERIGEE_TERMS, per_coeffs)})

# Latitude 2-harmonic model: lat = a*sin(w) + b*sin(3w), w = lon - node.
APOGEE_LAT_COEFFS = ({float(apo_lat[0])!r}, {float(apo_lat[1])!r})
PERIGEE_LAT_COEFFS = ({float(per_lat[0])!r}, {float(per_lat[1])!r})

# Mean apsis distances from the DE440 passages (AU).
APOGEE_MEAN_DIST_AU = {apo_dist!r}
PERIGEE_MEAN_DIST_AU = {per_dist!r}

# ============================================================================
# APOGEE CORRECTION TABLE METADATA
# ============================================================================
APOGEE_CT_JD_START: float = {APOGEE_CT_JD_START!r}
APOGEE_CT_STEP_DAYS: float = {APOGEE_CT_STEP_DAYS!r}
APOGEE_CT_COUNT: int = {APOGEE_CT_COUNT}

# ============================================================================
# PERIGEE CORRECTION TABLE METADATA
# ============================================================================
PERIGEE_CT_JD_START: float = {PERIGEE_CT_JD_START!r}
PERIGEE_CT_STEP_DAYS: float = {PERIGEE_CT_STEP_DAYS!r}
PERIGEE_CT_COUNT: int = {PERIGEE_CT_COUNT}

# ============================================================================
# APOGEE CORRECTIONS (arcseconds)
# ============================================================================
APOGEE_CORRECTIONS: tuple = (
'''
    parts = [header, _fmt_table(apo_corr)]
    parts.append(
        ")\n\n"
        "# ============================================================================\n"
        "# PERIGEE CORRECTIONS (arcseconds)\n"
        "# ============================================================================\n"
        "PERIGEE_CORRECTIONS: tuple = (\n"
    )
    parts.append(_fmt_table(per_corr))
    parts.append(")\n")

    out = Path(__file__).resolve().parents[1] / "libephemeris" / "lunar_apse_model.py"
    out.write_text("".join(parts), encoding="utf-8")
    print(f"wrote {out} ({out.stat().st_size} bytes)")


if __name__ == "__main__":
    raise SystemExit(main())
