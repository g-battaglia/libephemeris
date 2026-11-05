"""Canonical mean origin: independent controls, convergence and table identity."""
import importlib.util
from pathlib import Path
import sys
import math
import re

import pytest

from tools.generators import gmst_origin_reference as reference
from tools.generators import gen_gmst_canonical_table as generator


@pytest.fixture
def mp():
    # Use the existing optional validation installation, never fetch a package.
    if importlib.util.find_spec("mpmath") is None:
        site = Path(__file__).resolve().parents[2] / "target/v4/mpmath-site"
        if not (site / "mpmath").is_dir():
            pytest.fail("existing mpmath validation installation is required")
        sys.path.insert(0, str(site))
    import mpmath
    with mpmath.workdps(50):
        yield mpmath


def test_requests_keep_clock_grid_and_only_explicit_domain_points():
    assert len(reference.request_jds(None)) == 10
    endpoints = [2396757.5, 2506332.5, 2287184.5, 2688953.5, -3100015.5, 8000016.5]
    planned = list(reference.DEFAULT_JDS) + [jd for jd in endpoints if jd not in reference.DEFAULT_JDS]
    assert len(reference.request_jds(planned)) == 15
    assert reference.request_jds([2451545.0]) == [2451545.0]


@pytest.mark.parametrize("values", [[float("nan")], [float("inf")], [True], ["2451545"],
                                    [2451545.0, 2451545.0], [2451545.0 + 36525*2001],
                                    list(range(33))])
def test_bad_domains_nonfinite_duplicate_and_unbounded_requests_refuse(values):
    with pytest.raises(ValueError):
        reference.request_jds(values)


def test_composite_gauss_rule_is_exact_for_independent_polynomial_and_orientation(mp):
    for endpoint in (mp.mpf(2), mp.mpf(-2)):
        integral, panels = reference.integrate(mp, lambda t: t**31, endpoint, 16, mp.mpf("0.75"))
        assert panels == 3
        assert abs(integral - endpoint**32/32) < mp.mpf("1e-40")
    result, panels = reference.integrate(mp, lambda t: pytest.fail("zero interval must not evaluate"), 0, 16, 1)
    assert result == 0 and panels == 0


def test_quadrature_panel_bound_refuses_before_evaluation(mp):
    with pytest.raises(ValueError, match="panel bound"):
        reference.integrate(mp, lambda t: pytest.fail("must refuse before evaluation"), 801, 16, 1)


def test_analytic_series_derivative_on_synthetic_polynomial_and_harmonic(mp):
    # A separate exact polynomial and quarter-cycle harmonic, not Vondrak data.
    polynomial = (("1", "2", "3"), ("4", "0", "-1"))
    terms = reference.parse_terms(mp, "4 2 3 5 7")
    values, derivatives = reference.series(mp, mp.mpf(1), polynomial, terms)
    scale = mp.pi/648000
    assert abs(values[0]/scale - 11) < mp.mpf("1e-45")
    assert abs(values[1]/scale - 10) < mp.mpf("1e-45")
    assert abs(derivatives[0]/scale - (8 - mp.pi)) < mp.mpf("1e-45")
    assert abs(derivatives[1]/scale - (-2 - 3*mp.pi/2)) < mp.mpf("1e-45")


def test_constant_latitude_connection_has_known_geometric_holonomy(mp):
    # q=(4/5 cos(t),4/5 sin(t),3/5) gives ds/dt=-(1-Z)=-2/5.
    x, y, z = mp.mpf(4)/5, mp.mpf(0), mp.mpf(3)/5
    dx, dy = mp.mpf(0), mp.mpf(4)/5
    rate = -(x*dy - y*dx)/(1+z)
    integral, _ = reference.integrate(mp, lambda t: rate, 2*mp.pi, 16, 1)
    assert abs(integral + 4*mp.pi/5) < mp.mpf("1e-45")


def test_equinox_gauge_has_independent_known_rotation(mp):
    q = [mp.mpf(0), mp.mpf(0), mp.mpf(1)]
    angle, epsilon = mp.mpf("0.3"), mp.mpf("0.4")
    e = [mp.sin(epsilon)*mp.sin(angle), -mp.sin(epsilon)*mp.cos(angle), mp.cos(epsilon)]
    assert abs(reference.equinox_phase(mp, q, e) + angle) < mp.mpf("1e-45")
    assert abs(reference.wrap(mp, 4*mp.pi + angle) - angle) < mp.mpf("1e-45")


def test_singular_clamped_nonunit_and_nonfinite_poles_refuse(mp):
    with pytest.raises(ValueError, match="clamped"):
        reference.upper_pole(mp, mp.mpf(1), mp.mpf(0))
    with pytest.raises(ValueError, match="clamped"):
        reference.upper_pole(mp, mp.mpf(2), mp.mpf(0))
    north = [mp.mpf(0), mp.mpf(0), mp.mpf(1)]
    with pytest.raises(ValueError, match="undefined equinox"):
        reference.equinox_phase(mp, north, north)
    with pytest.raises(ValueError, match="south-pole"):
        reference.equinox_phase(mp, [mp.mpf(0), mp.mpf(0), mp.mpf(-1)], north)
    with pytest.raises(ValueError, match="nonunit"):
        reference.equinox_phase(mp, [mp.mpf(0), mp.mpf(0), mp.mpf(2)], north)
    with pytest.raises(ValueError, match="nonfinite"):
        reference.equinox_phase(mp, [mp.nan, mp.mpf(0), mp.mpf(1)], north)


def test_published_decimal_table_shapes_are_explicit_without_evaluation(mp):
    assert len(reference.parse_terms(mp, reference.XY_TERMS)) == 14
    assert len(reference.parse_terms(mp, reference.PQ_TERMS)) == 8
    assert all(len(axis) == 4 for axis in reference.XY_POLY + reference.PQ_POLY)


def test_existing_output_refuses_before_optional_scientific_imports(tmp_path):
    output = tmp_path / "existing.json"
    output.write_text("preserve")
    with pytest.raises(ValueError, match="already exists"):
        reference.main(["--output", str(output)])
    assert output.read_text() == "preserve"


@pytest.fixture(scope="module")
def canonical():
    return reference.MeanOriginReference(70, 16, 10)


@pytest.fixture(scope="module")
def fine():
    return reference.MeanOriginReference(90, 32, 5)


def test_own_context_never_changes_global_mpmath_precision(canonical):
    import mpmath
    before = mpmath.mp.dps
    canonical.mean_sidereal_hours(2451545.0, 2451545.0)
    assert mpmath.mp.dps == before
    assert canonical.mp.dps == 70
    assert canonical.mp is not mpmath.mp


def test_exact_float_input_retains_low_bits_and_clock_separation(canonical):
    mp = canonical.mp
    tt = math.nextafter(2451545.0, math.inf)
    assert canonical.centuries_from_jd(tt) == mp.mpf(1)/(2**31*36525)
    fixed_tt = 2451545.125
    h0 = canonical.mean_sidereal_hours(2451545.0, fixed_tt)
    h1 = canonical.mean_sidereal_hours(2451545.25, fixed_tt)
    expected = mp.mpf("1.00273781191135448")*6
    assert abs(reference.wrap(mp, (h1-h0-expected)*mp.pi/12)) < mp.mpf("1e-65")
    # TT changes only the celestial origin; it does not advance ERA.
    origin0 = canonical.mean_origin_centuries(0)
    origin1 = canonical.mean_origin_centuries(mp.mpf(1)/100)
    actual = canonical.mean_sidereal_hours(2451545.0, 2451910.25)-canonical.mean_sidereal_hours(2451545.0, 2451545.0)
    assert abs(reference.wrap(mp, actual*mp.pi/12-(origin1-origin0))) < mp.mpf("1e-65")


def test_j2000_anchor_and_published_secular_term(canonical):
    mp = canonical.mp
    assert canonical.locator_centuries(0) == 0
    assert canonical.relative_eo_centuries(0) == 0
    assert abs(canonical.mean_origin_centuries(0)*648000/mp.pi-mp.mpf("0.014506")) < mp.mpf("1e-63")
    for t in (-1, mp.mpf("0.5"), 2):
        unwrapped = canonical.mean_origin_centuries(t)+canonical.relative_eo_centuries(t)
        assert abs(reference.wrap(mp, unwrapped-(mp.mpf("0.014506")-mp.mpf("0.003872")*t)*mp.pi/648000)) < mp.mpf("1e-63")


def test_cached_fixed_panels_have_analytic_orientation_and_reuse():
    r = reference.MeanOriginReference(60, 16, 10)
    calls = []
    r._connection = lambda t: calls.append(t) or t**4
    mp = r.mp
    for t in (mp.mpf(21), mp.mpf(22), mp.mpf(-15), mp.mpf(-21)):
        assert abs(r.locator_centuries(t)-t**5/5) < mp.mpf("1e-50")
    count = len(calls)
    r.locator_centuries(23)
    assert len(calls)-count == 16  # Cached 0->20 panels are not reacquired.
    assert {0, 1, 2, -1, -2} <= set(r._locator_nodes)


@pytest.mark.parametrize("kwargs", [{"dps": 49}, {"dps": True}, {"order": 8},
                                    {"panel_width": 0}, {"panel_width": -1},
                                    {"panel_width": "nan"}, {"panel_width": 1}])
def test_bad_reference_configurations_refuse(kwargs):
    with pytest.raises(ValueError):
        reference.MeanOriginReference(**kwargs)


@pytest.mark.parametrize("t", [2001, -2001, "nan", "inf", True])
def test_reference_domain_and_nonfinite_centuries_refuse(canonical, t):
    with pytest.raises(ValueError):
        canonical.locator_centuries(t)


def test_reference_nonfinite_clock_refuses(canonical):
    for ut1, tt in ((math.inf, 2451545.0), (2451545.0, math.nan), (True, 2451545.0)):
        with pytest.raises(ValueError):
            canonical.mean_sidereal_hours(ut1, tt)


def test_two_reference_configurations_at_far_modern_and_nonnode_epochs(canonical, fine):
    mp = fine.mp
    # 1e-6 microarcsecond agreement is the owner reference-convergence target.
    threshold = mp.mpf("1e-6")*mp.pi/(648000*10**6)
    for t in (-2000, -987.125, -1.5, 0, 0.5, 123.456789, 2000):
        s0, s1 = canonical.locator_centuries(t), fine.locator_centuries(t)
        eo0, eo1 = canonical.relative_eo_centuries(t), fine.relative_eo_centuries(t)
        assert abs(mp.mpf(s0)-s1) < threshold
        assert abs(reference.wrap(mp, mp.mpf(eo0)-eo1)) < threshold


def test_small_modern_gmst06_control_without_fitted_boundary(canonical):
    import erfa
    mp = canonical.mp
    for year in (1850, 1900, 2000, 2025, 2050):
        a, b = erfa.epj2jd(year)
        jd = float(a+b)
        actual = canonical.mean_sidereal_hours(jd, jd)*mp.pi/12
        expected = float(erfa.gmst06(2451545.0, jd-2451545.0, 2451545.0, jd-2451545.0))
        error = abs(reference.wrap(mp, actual-expected))*648000/mp.pi
        assert error <= mp.mpf("0.001"), (year, error)


def test_published_source_analytic_rate_matches_independent_high_precision_differentiation(canonical):
    mp = canonical.mp
    t = mp.mpf("-13.125")
    terms = reference.parse_terms(mp, reference.XY_TERMS)
    values, rates = reference.series(mp, t, reference.XY_POLY, terms)
    for axis in (0, 1):
        def coordinate(x):
            coefficients = [mp.mpf(c) for c in reference.XY_POLY[axis]]
            value = mp.mpf(0)
            for coefficient in reversed(coefficients):
                value = value*x+coefficient
            for row in terms:
                angle = 2*mp.pi*x/row[0]
                value += row[axis+1]*mp.cos(angle)+row[axis+3]*mp.sin(angle)
            return value*mp.pi/648000
        assert abs(values[axis]-coordinate(t)) < mp.mpf("1e-65")
        assert abs(rates[axis]-mp.diff(coordinate, t)) < mp.mpf("1e-65")


def test_generator_small_exact_grid_and_gauss_polynomial_moments(fine):
    text = generator.render_table(fine, start=-10, step=10, count=3)
    assert "CIO_S: [f64; 3]" in text
    assert re.search(r"\b0\.0,\s+// T=0", text)
    assert len(text.splitlines()) < 100
    for order in (16, 32):
        match = re.search(rf"GAUSS{order}:.*?= \[(.*?)\];", text, re.S)
        rows = [(float(x), float(w)) for x, w in re.findall(r"\(([^,]+), ([^)]+)\)", match[1])]
        assert len(rows) == order
        assert all(-1 < x < 1 and w > 0 for x, w in rows)
        for k in (0, 1, 2, 7, 2*order-2):
            target = 0 if k % 2 else 2/(k+1)
            assert abs(math.fsum(w*x**k for x, w in rows)-target) < 8e-16


def test_generated_full_table_shapes_nodes_and_binary64_rounding(fine):
    path = Path(__file__).resolve().parents[2]/"core/src/time/canonical_sidereal_table.rs"
    text = path.read_text()
    assert len(text.splitlines()) < 2000
    assert "CIO_S: [f64; 401]" in text
    assert "NODE_START_CENTURIES: f64 = -2000.0" in text
    assert "NODE_STEP_CENTURIES: f64 = 10.0" in text
    rows = re.findall(r"^    ([^,]+),\s+// T=(-?\d+)$", text, re.M)
    assert len(rows) == 401
    assert [int(t) for _, t in rows] == list(range(-2000, 2001, 10))
    for index in (0, 99, 199, 200, 201, 333, 400):
        value, t = rows[index]
        assert float(value) == float(fine.locator_centuries(int(t)))


def test_rust_float_roundtrips_tiny_negative_and_nonfinite_refusal():
    for value in (0, -0.0, -1.25, math.nextafter(0.0, 1), math.nextafter(1.0, 2)):
        assert float(generator.rust_float(value)).hex() == float(value).hex()
    with pytest.raises(ValueError):
        generator.rust_float(math.inf)
