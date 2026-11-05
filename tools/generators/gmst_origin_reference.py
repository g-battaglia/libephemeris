#!/usr/bin/env python3
"""High-precision canonical VCW mean sidereal time, evaluated offline.

This independent diagnostic integrates the mean-pole connection from J2000:
    ds/dT = -(X*dY/dT - Y*dX/dT)/(1+Z), T in Julian centuries.
The equinox uses x=unit(q cross e), y=q cross x and
    EO = s - atan2(y dot a, x dot a),
    a=(1-X*X/(1+Z), -X*Y/(1+Z), -X).
The relative EO is zero at J2000. The canonical convention adds the published
IAU2006 datum and P03 Eq42 secular term:
    GMST = ERA - relative_EO + (0.014506 - 0.003872*T) arcseconds.
UT1 drives ERA alone; TT drives the mean celestial origin alone. No fitted
boundary, nutation, apparent sidereal time or Delta-T model enters this tool.

Decimal series: Vondrak, Capitaine & Wallace (2011), A&A 534 A22, with the
2012 corrigendum, as in ERFA ltpequ/ltpecl. Analytic derivatives avoid f64
finite-difference noise. Fixed Gauss-Legendre 16/32-node and 70/90-digit
configurations estimate numerical convergence, not physical model accuracy.
PyERFA poles are separately compared at requested endpoints; they do not
validate the transport integral. This is offline validation, not runtime code.

Sources:
https://www.aanda.org/articles/aa/pdf/2011/10/aa17274-11.pdf
https://github.com/liberfa/erfa/blob/v2.0.1/src/ltpequ.c
https://github.com/liberfa/erfa/blob/v2.0.1/src/ltpecl.c
The coefficients below are the same decimal physical data used by the native
Vondrak model; this is an independent high-precision evaluation/integration,
not a second independently observed precession model.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile

# Each row: period [centuries], X/P cosine, Y/Q cosine, X/P sine, Y/Q sine
# [arcseconds]. Decimal strings enter mpmath without a binary64 intermediate.
XY_TERMS = """
256.75 -819.940624 75004.344875 81491.287984 1558.515853
708.15 -8444.676815 624.033993 787.163481 7774.939698
274.20 2600.009459 1251.136893 1251.296102 -2219.534038
241.45 2755.175630 -1102.212834 -1257.950837 -2523.969396
2309.00 -167.659835 -2660.664980 -2966.799730 247.850422
492.20 871.855056 699.291817 639.744522 -846.485643
396.10 44.769698 153.167220 131.600209 -1393.124055
288.90 -512.313065 -950.865637 -445.040117 368.526116
231.10 -819.415595 499.754645 584.522874 749.045012
1610.00 -538.071099 -145.188210 -89.756563 444.704518
620.00 -189.793622 558.116553 524.429630 235.934465
157.87 -402.922932 -23.923029 -13.549067 374.049623
220.30 179.516345 -165.405086 -210.157124 -171.330180
1200.00 -9.814756 9.344131 -44.919798 -22.899655
"""
PQ_TERMS = """
708.15 -5486.751211 -684.661560 667.666730 -5523.863691
2309.00 -17.127623 2446.283880 -2354.886252 -549.747450
1620.00 -617.517403 399.671049 -428.152441 -310.998056
492.20 413.442940 -356.652376 376.202861 421.535876
1183.00 78.614193 -186.387003 184.778874 -36.776172
622.00 -180.732815 -316.800070 335.321713 -145.278396
882.00 -87.676083 198.296701 -185.138669 -34.744450
547.00 46.140315 101.135679 -120.972830 22.885731
"""
XY_POLY = (("5453.282155", "0.4252841", "-0.00037173", "-0.000000152"),
           ("-73750.930350", "-0.7675452", "-0.00018725", "0.000000231"))
PQ_POLY = (("5851.607687", "-0.1189000", "-0.00028913", "0.000000101"),
           ("-1600.886300", "1.1689818", "-0.00000020", "-0.000000437"))
DEFAULT_JDS = (2396757.5, 2396758.5, 2396759.5, 2451545.0, 2469806.5,
               2469807.5, 2469808.5, 2488147.5, 2561195.5, 2597720.5)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def request_jds(values):
    result = list(DEFAULT_JDS if not values else values)
    require(0 < len(result) <= 32, "require 1..32 explicit scalar TT requests")
    require(all(type(x) in (int, float) and math.isfinite(x) for x in result), "nonfinite TT request")
    require(len({float(x).hex() for x in result}) == len(result), "duplicate TT request")
    require(all(abs((x - 2451545.0) / 36525.0) <= 2000 for x in result),
            "outside published +/-200 kyr model interval; no extrapolation")
    return result


def parse_terms(mp, text):
    rows = [tuple(mp.mpf(word) for word in line.split()) for line in text.splitlines() if line.strip()]
    require(all(len(row) == 5 and row[0] > 0 for row in rows), "invalid decimal coefficient row")
    return rows


def series(mp, t, polynomial, terms):
    """Two pole coordinates and analytic derivatives per Julian century."""
    values, derivatives = [], []
    scale = mp.pi / mp.mpf(648000)
    for axis in range(2):
        poly = [mp.mpf(value) for value in polynomial[axis]]
        value = mp.fsum(coefficient * t**power for power, coefficient in enumerate(poly))
        derivative = mp.fsum(power * poly[power] * t**(power - 1) for power in range(1, len(poly)))
        for period, cos0, cos1, sin0, sin1 in terms:
            frequency = 2 * mp.pi / period
            angle = frequency * t
            cosine, sine = (cos0, sin0) if axis == 0 else (cos1, sin1)
            value += cosine * mp.cos(angle) + sine * mp.sin(angle)
            derivative += frequency * (-cosine * mp.sin(angle) + sine * mp.cos(angle))
        values.append(value * scale)
        derivatives.append(derivative * scale)
    return values, derivatives


def upper_pole(mp, x, y):
    radicand = 1 - x*x - y*y
    # ERFA clamps outside its spherical representation. This diagnostic refuses
    # that path: its smooth-pole derivative/connection would be unjustified.
    require(mp.isfinite(radicand) and radicand > 0, "pole reaches clamped/singular representation")
    return [x, y, mp.sqrt(radicand)]


def cross(a, b):
    return [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]


def dot(mp, a, b):
    return mp.fsum(x*y for x, y in zip(a, b))


def equinox_phase(mp, q, e):
    """EO gauge angle in the J2000 mean equatorial frame, radians."""
    require(all(mp.isfinite(value) for value in q + e), "nonfinite pole")
    for pole in (q, e):
        require(abs(dot(mp, pole, pole) - 1) < mp.mpf("1e-30"), "nonunit pole")
    require(1 + q[2] > mp.mpf("1e-30"), "south-pole connection singularity")
    x = cross(q, e)
    norm = mp.sqrt(dot(mp, x, x))
    require(norm > mp.mpf("1e-30"), "undefined equinox: coincident poles")
    x = [value / norm for value in x]
    y = cross(q, x)
    a = [1 - q[0]*q[0]/(1 + q[2]), -q[0]*q[1]/(1 + q[2]), -q[0]]
    p, r = dot(mp, x, a), dot(mp, y, a)
    require(p*p + r*r > mp.mpf("1e-30"), "undefined equinox gauge angle")
    return mp.atan2(r, p)


def poles_and_rate(mp):
    xy, pq = parse_terms(mp, XY_TERMS), parse_terms(mp, PQ_TERMS)
    epsilon = mp.mpf("84381.406") * mp.pi / 648000

    def equator(t):
        (x, y), (dx, dy) = series(mp, t, XY_POLY, xy)
        q = upper_pole(mp, x, y)
        return q, dx, dy

    def poles(t):
        q, _, _ = equator(t)
        (p, r), _ = series(mp, t, PQ_POLY, pq)
        z = upper_pole(mp, p, r)[2]
        e = [p, -r*mp.cos(epsilon) - z*mp.sin(epsilon),
             -r*mp.sin(epsilon) + z*mp.cos(epsilon)]
        return q, e

    def connection(t):
        q, dx, dy = equator(t)
        return -(q[0]*dy - q[1]*dx)/(1 + q[2])

    return poles, connection


def integrate(mp, function, end, order, panel_width):
    """Fixed composite GL rule, signed orientation, at most 800 panels."""
    if end == 0:
        return mp.mpf(0), 0
    panels = int(mp.ceil(abs(end) / panel_width))
    require(1 <= panels <= 800, "quadrature panel bound exceeded")
    nodes, weights = mp.gauss_quadrature(order, "legendre")
    half = end / (2 * panels)
    terms = []
    for index in range(panels):
        midpoint = (2*index + 1) * half
        terms.append(half * mp.fsum(weights[j] * function(midpoint + half*nodes[j]) for j in range(order)))
    return mp.fsum(terms), panels


def wrap(mp, angle):
    return (angle + mp.pi) % (2*mp.pi) - mp.pi


class MeanOriginReference:
    """Cached independent decimal-series reference in its own mpmath context.

    ``locator_centuries(T)`` returns s(T) in radians, with s(0)=0.
    ``relative_eo_centuries(T)`` returns relative EO modulo one turn.
    ``mean_sidereal_hours(ut1, tt)`` accepts represented scalar Julian days and
    returns an mpf in [0,24). Floats enter exactly, without decimal reformatting.
    The supported reference domain is |T|<=2000 Julian centuries. Ordinary
    fixed quadrature supplies a convergence reference, not interval bounds.
    """

    def __init__(self, dps=70, order=16, panel_width=10):
        import mpmath
        require(type(dps) is int and dps >= 50, "reference requires at least 50 digits")
        require(type(order) is int and order in (16, 32), "require GL order 16 or 32")
        self.mp = mpmath.mp.clone()
        self.mp.dps = dps
        self.dps, self.order = dps, order
        self.panel_width = self.mp.mpf(panel_width)
        require(self.mp.isfinite(self.panel_width) and self.panel_width > 0,
                "panel width must be finite and positive")
        require(self.mp.ceil(2000/self.panel_width) <= 800,
                "quadrature panel bound exceeded")
        nodes, weights = self.mp.gauss_quadrature(order, "legendre")
        self.gauss_nodes = tuple(nodes)
        self.gauss_weights = tuple(weights)
        self._poles, self._connection = poles_and_rate(self.mp)
        self._phase_zero = equinox_phase(self.mp, *self._poles(self.mp.mpf(0)))
        self._locator_nodes = {0: self.mp.mpf(0)}

    def _centuries(self, value):
        require(not isinstance(value, bool), "invalid century value")
        t = self.mp.mpf(value)
        require(self.mp.isfinite(t) and abs(t) <= 2000,
                "outside published +/-200 kyr model interval; no extrapolation")
        return t

    def centuries_from_jd(self, jd):
        """TT centuries from the exact represented scalar JD."""
        require(not isinstance(jd, bool), "invalid TT Julian day")
        return self._centuries((self.mp.mpf(jd)-2451545)/36525)

    def _panel(self, a, b):
        if a == b:
            return self.mp.mpf(0)
        half, centre = (b-a)/2, (a+b)/2
        return half*self.mp.fsum(w*self._connection(centre+half*x)
                                 for x, w in zip(self.gauss_nodes, self.gauss_weights))

    def locator_centuries(self, t):
        """Mean-pole connection integral from J2000, in radians.

        Whole lattice panels are cached. Only the correction from the nearest
        cached lattice node to the exact requested T is evaluated per request;
        halfway ties choose the lower node. No epoch is replaced by that node.
        """
        t = self._centuries(t)
        index = int(self.mp.ceil(t/self.panel_width-self.mp.mpf("0.5")))
        # At an edge an arbitrary panel width need not divide 2000. Keep the
        # anchor inside the reference domain and integrate the remaining tail.
        limit = int(self.mp.floor(2000/self.panel_width))
        index = max(-limit, min(limit, index))
        if index not in self._locator_nodes:
            direction = 1 if index > 0 else -1
            previous = max((k for k in self._locator_nodes if 0 <= k <= index), default=0) if direction > 0 else min(
                (k for k in self._locator_nodes if index <= k <= 0), default=0)
            for node in range(previous+direction, index+direction, direction):
                a, b = (node-direction)*self.panel_width, node*self.panel_width
                self._locator_nodes[node] = self._locator_nodes[node-direction]+self._panel(a, b)
        anchor = index*self.panel_width
        return self._locator_nodes[index]+self._panel(anchor, t)

    def poles_centuries(self, t):
        """Unit mean equator/ecliptic poles, using published decimal data."""
        return self._poles(self._centuries(t))

    def connection_centuries(self, t):
        """ds/dT in radians per Julian century, using analytic derivatives."""
        return self._connection(self._centuries(t))

    def relative_eo_centuries(self, t):
        """Relative mean equation of origins in radians, modulo a turn."""
        t = self._centuries(t)
        phase = equinox_phase(self.mp, *self._poles(t))
        return wrap(self.mp, self.locator_centuries(t)-(phase-self._phase_zero))

    def mean_origin_centuries(self, t):
        """GMST minus ERA, radians modulo a turn, including datum/P03 term."""
        t = self._centuries(t)
        conventional = (self.mp.mpf("0.014506")-self.mp.mpf("0.003872")*t)*self.mp.pi/648000
        return wrap(self.mp, -self.relative_eo_centuries(t)+conventional)

    def era_radians(self, ut1):
        """ERA for the exact represented UT1 JD, using published constants."""
        require(not isinstance(ut1, bool), "invalid UT1 Julian day")
        days = self.mp.mpf(ut1)-2451545
        require(self.mp.isfinite(days), "nonfinite UT1 Julian day")
        turns = self.mp.mpf("0.7790572732640")+self.mp.mpf("1.00273781191135448")*days
        return 2*self.mp.pi*(turns % 1)

    def mean_sidereal_hours(self, ut1, tt):
        """Canonical mean sidereal hours in [0,24), separate UT1 and TT."""
        t = self.centuries_from_jd(tt)
        return ((self.era_radians(ut1)+self.mean_origin_centuries(t))/(2*self.mp.pi)*24) % 24


def evaluate(mp, jd_tt, dps, order, panel_width):
    with mp.workdps(dps):
        # Preserve the exact scalar binary64 request, including its low bits.
        t = (mp.mpf(jd_tt) - mp.mpf("2451545.0")) / 36525
        poles, connection = poles_and_rate(mp)
        s, panels = integrate(mp, connection, t, order, mp.mpf(str(panel_width)))
        q, e = poles(t)
        phase = equinox_phase(mp, q, e)
        zero_phase = equinox_phase(mp, *poles(mp.mpf(0)))
        eo = wrap(mp, s - (phase - zero_phase))
        return {"relative_eo_rad": +eo, "connection_s_rad": +s,
                "equinox_phase_rad": +phase, "j2000_equinox_phase_rad": +zero_phase,
                "equator": [+v for v in q], "ecliptic": [+v for v in e], "panels": panels}


def reference_row(mp, erfa, jd_tt, numerical_tolerance, references=None):
    """Two-configuration convergence and optional same-epoch ERFA diagnostics."""
    coarse, fine = references or (MeanOriginReference(70, 16, 10), MeanOriginReference(90, 32, 5))
    t0, t1 = coarse.centuries_from_jd(jd_tt), fine.centuries_from_jd(jd_tt)
    s0, s1 = coarse.locator_centuries(t0), fine.locator_centuries(t1)
    eo0, eo1 = coarse.relative_eo_centuries(t0), fine.relative_eo_centuries(t1)
    m = fine.mp
    errors = {"locator": abs(s1-m.mpf(s0)), "relative_eo": abs(wrap(m, eo1-m.mpf(eo0)))}
    hours = fine.mean_sidereal_hours(jd_tt, jd_tt)
    erfa_angle = float(erfa.gmst06(2451545.0, jd_tt-2451545.0, 2451545.0, jd_tt-2451545.0))
    gmst_difference = wrap(m, hours*m.pi/12-m.mpf(erfa_angle))*648000/m.pi
    return {"jd_tt": jd_tt, "jd_tt_hex": float(jd_tt).hex(),
            "relative_eo_rad": m.nstr(eo1, 65), "connection_s_rad": m.nstr(s1, 65),
            "mean_sidereal_hours": m.nstr(hours, 65),
            "converged": all(value <= m.mpf(numerical_tolerance) for value in errors.values()),
            "estimated_numerical_errors_rad": {k: m.nstr(v, 12) for k, v in errors.items()},
            "pyerfa_gmst06_difference_arcsec": m.nstr(gmst_difference, 18)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jd-tt", type=float, action="append", help="repeat to replace the ten default TT JDs")
    parser.add_argument("--mpmath-site", type=Path, help="existing optional mpmath installation")
    parser.add_argument("--numerical-tolerance-rad", default="1e-26", help="convergence target, not physical accuracy")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    epochs = request_jds(args.jd_tt)
    require(not args.output.exists() and not args.output.is_symlink(), "output already exists")
    if args.mpmath_site:
        sys.path.insert(0, str(args.mpmath_site.resolve()))
    import mpmath as mp
    import erfa
    require(mp.isfinite(mp.mpf(args.numerical_tolerance_rad)) and mp.mpf(args.numerical_tolerance_rad) > 0,
            "invalid numerical convergence target")
    references = MeanOriginReference(70, 16, 10), MeanOriginReference(90, 32, 5)
    rows = [reference_row(mp, erfa, jd, args.numerical_tolerance_rad, references) for jd in epochs]
    result = {"schema": 2, "kind": "canonical_vcw_mean_sidereal_time", "diagnostic_only": True,
              "scientific_acceptance": False, "angle_branch": "relative EO modulo 2pi in [-pi,pi)",
              "initial_s": "zero; relative EO zero at J2000",
              "gmst_initial_arcsec_included": "0.014506", "p03_secular_arcsec_per_century": "-0.003872",
              "quadrature": {"orders": [16, 32], "panel_widths_centuries": [10, 5], "dps": [70, 90],
                             "numerical_tolerance_rad": args.numerical_tolerance_rad},
              "versions": {"mpmath": mp.__version__, "pyerfa": erfa.__version__},
              "tool_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "rows": rows,
              "all_converged": all(row["converged"] for row in rows)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".gmst-origin-", dir=args.output.parent)
    try:
        with os.fdopen(fd, "w") as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
            stream.write("\n")
        os.link(temporary, args.output)
    finally:
        Path(temporary).unlink()
    return 0 if result["all_converged"] else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError, ImportError) as error:
        print(f"GMST origin reference: {error}", file=sys.stderr)
        raise SystemExit(1) from error
