"""Clean-room tests for hypothetical and user-supplied orbital elements."""

from __future__ import annotations

import ast
import math
from pathlib import Path

import pytest

from libephemeris.exceptions import UnknownBodyError
from libephemeris.hypothetical import (
    ADMETOS,
    ADMETOS_KEPLERIAN_ELEMENTS,
    APOLLON,
    APOLLON_KEPLERIAN_ELEMENTS,
    CUPIDO,
    CUPIDO_KEPLERIAN_ELEMENTS,
    FICTITIOUS_ORBITAL_ELEMENTS,
    HADES,
    HADES_KEPLERIAN_ELEMENTS,
    HARRINGTON,
    HYPOTHETICAL_ELEMENTS,
    HYPOTHETICAL_PROVENANCE,
    ISIS,
    KRONOS,
    KRONOS_KEPLERIAN_ELEMENTS,
    NIBIRU,
    PICKERING_PLANET_X_ELEMENTS,
    POSEIDON,
    POSEIDON_KEPLERIAN_ELEMENTS,
    PROSERPINA,
    TRANSPLUTO_KEPLERIAN_ELEMENTS,
    TransplutoKeplerianElements,
    URANIAN_KEPLERIAN_ELEMENTS,
    URANIAN_ELEMENTS,
    VULCAN,
    VULCAN_ELEMENTS,
    VULKANUS,
    VULKANUS_KEPLERIAN_ELEMENTS,
    WALDEMATH,
    WALDEMATH_ELEMENTS,
    WHITE_MOON,
    ZEUS,
    ZEUS_KEPLERIAN_ELEMENTS,
    AdmetosKeplerianElements,
    ApollonKeplerianElements,
    CupidoKeplerianElements,
    HadesKeplerianElements,
    HypotheticalElements,
    KronosKeplerianElements,
    OrbitalElements,
    PoseidonKeplerianElements,
    TPolynomial,
    UranianElements,
    VulkanusKeplerianElements,
    ZeusKeplerianElements,
    _calc_keplerian_position,
    _calc_keplerian_position_raw,
    _solve_kepler_equation,
    calc_admetos,
    calc_apollon,
    calc_cupido,
    calc_fictitious_position,
    calc_hades,
    calc_hypothetical_position,
    calc_kronos,
    calc_orbital_position,
    calc_planet_x_lowell,
    calc_planet_x_pickering,
    calc_poseidon,
    calc_proserpina,
    calc_transpluto,
    calc_transpluto_position,
    calc_uranian_longitude,
    calc_uranian_planet,
    calc_uranian_position,
    calc_vulcan,
    calc_vulkanus,
    calc_waldemath,
    calc_waldemath_position,
    calc_white_moon_position,
    calc_zeus,
    get_orbital_body_by_name,
    is_hypothetical_body,
    list_hypothetical_bodies,
    load_bundled_fictitious_orbits,
)


J2000 = 2451545.0

_RC7_LEGACY_PUBLIC_NAMES = {
    "ADMETOS_KEPLERIAN_ELEMENTS",
    "APOLLON_KEPLERIAN_ELEMENTS",
    "AdmetosKeplerianElements",
    "ApollonKeplerianElements",
    "CUPIDO_KEPLERIAN_ELEMENTS",
    "CupidoKeplerianElements",
    "HADES_KEPLERIAN_ELEMENTS",
    "HadesKeplerianElements",
    "KRONOS_KEPLERIAN_ELEMENTS",
    "KronosKeplerianElements",
    "POSEIDON_KEPLERIAN_ELEMENTS",
    "PoseidonKeplerianElements",
    "URANIAN_ELEMENTS",
    "UranianElements",
    "VULKANUS_KEPLERIAN_ELEMENTS",
    "VulkanusKeplerianElements",
    "ZEUS_KEPLERIAN_ELEMENTS",
    "ZeusKeplerianElements",
}


def test_compatibility_tables_cover_all_restored_fictitious_bodies() -> None:
    assert set(URANIAN_KEPLERIAN_ELEMENTS) == set(range(CUPIDO, POSEIDON + 1))
    assert set(HYPOTHETICAL_ELEMENTS) == {ISIS, PROSERPINA}
    assert TRANSPLUTO_KEPLERIAN_ELEMENTS is HYPOTHETICAL_ELEMENTS[ISIS]
    assert type(TRANSPLUTO_KEPLERIAN_ELEMENTS) is TransplutoKeplerianElements
    assert TRANSPLUTO_KEPLERIAN_ELEMENTS.name == "Transpluto"
    assert set(FICTITIOUS_ORBITAL_ELEMENTS) == set(range(NIBIRU, 55))
    assert VULCAN_ELEMENTS.name == "Vulcan"
    assert WALDEMATH_ELEMENTS.name == "Waldemath"
    assert set(HYPOTHETICAL_PROVENANCE) == set(range(CUPIDO, WALDEMATH + 1))

    rows = load_bundled_fictitious_orbits()
    assert len(rows) == 19
    row = get_orbital_body_by_name(rows, "Harrington")
    assert row is not None
    assert row.epoch_jd == 2374696.5
    assert row.mean_anomaly.evaluate(0.0) == 0.0
    assert row.semi_axis == 101.2
    assert row.eccentricity.evaluate(0.0) == 0.411
    assert row.arg_perihelion.evaluate(0.0) == 208.5
    assert row.asc_node.evaluate(0.0) == 275.4
    assert row.inclination.evaluate(0.0) == 32.4


def test_removed_element_constants_are_exported_at_package_import() -> None:
    import libephemeris as ephem

    assert ephem.VULCAN_ELEMENTS is VULCAN_ELEMENTS
    assert ephem.WALDEMATH_ELEMENTS is WALDEMATH_ELEMENTS
    assert ephem.PICKERING_PLANET_X_ELEMENTS is PICKERING_PLANET_X_ELEMENTS
    assert ephem.TransplutoKeplerianElements is TransplutoKeplerianElements
    assert URANIAN_KEPLERIAN_ELEMENTS
    assert TRANSPLUTO_KEPLERIAN_ELEMENTS is not None


def test_rc7_legacy_uranian_public_objects_keep_exact_types_and_values() -> None:
    expected_mean_elements = {
        CUPIDO: ("Cupido", 241.2067, 1.091437, 1.5, 21.6, 83.0),
        HADES: ("Hades", 176.0581, 0.736380, 1.5, 258.0, 56.0),
        ZEUS: ("Zeus", 32.1893, 0.532664, 1.0, 141.0, 41.0),
        KRONOS: ("Kronos", 213.2096, 0.420481, 1.0, 98.0, 32.0),
        APOLLON: ("Apollon", 71.8925, 0.341403, 1.0, 326.0, 26.0),
        ADMETOS: ("Admetos", 142.2269, 0.283756, 1.0, 250.0, 22.0),
        VULKANUS: ("Vulkanus", 195.6753, 0.240116, 1.0, 178.0, 18.0),
        POSEIDON: ("Poseidon", 274.4073, 0.207016, 1.0, 105.0, 16.0),
    }
    assert set(URANIAN_ELEMENTS) == set(expected_mean_elements)
    for body_id, expected in expected_mean_elements.items():
        elements = URANIAN_ELEMENTS[body_id]
        assert type(elements) is UranianElements
        assert (
            elements.name,
            elements.L0,
            elements.n,
            elements.amplitude,
            elements.phase,
            elements.phase_rate,
        ) == expected

    expected_keplerian = [
        (
            CUPIDO_KEPLERIAN_ELEMENTS,
            CupidoKeplerianElements,
            (
                "Cupido",
                2415020.0,
                40.99837,
                0.0,
                0.0,
                0.0,
                0.0,
                105.301693,
                0.0037945179,
            ),
        ),
        (
            HADES_KEPLERIAN_ELEMENTS,
            HadesKeplerianElements,
            (
                "Hades",
                2415020.0,
                50.66744,
                0.00245,
                1.05,
                148.1796,
                161.3339,
                26.850162,
                0.00278759,
            ),
        ),
        (
            ZEUS_KEPLERIAN_ELEMENTS,
            ZeusKeplerianElements,
            ("Zeus", 2415020.0, 59.21436, 0.0, 0.0, 0.0, 0.0, 104.289095, 0.002220375),
        ),
        (
            KRONOS_KEPLERIAN_ELEMENTS,
            KronosKeplerianElements,
            ("Kronos", 2415020.0, 64.8169, 0.0, 0.0, 0.0, 0.0, 17.111353, 0.0019351856),
        ),
        (
            APOLLON_KEPLERIAN_ELEMENTS,
            ApollonKeplerianElements,
            (
                "Apollon",
                2415020.0,
                70.36118,
                0.0,
                0.0,
                0.0,
                0.0,
                138.565328,
                0.0017177599,
            ),
        ),
        (
            ADMETOS_KEPLERIAN_ELEMENTS,
            AdmetosKeplerianElements,
            (
                "Admetos",
                2415020.0,
                73.736396,
                0.0,
                0.0,
                0.0,
                0.0,
                350.613913,
                0.0016016766,
            ),
        ),
        (
            VULKANUS_KEPLERIAN_ELEMENTS,
            VulkanusKeplerianElements,
            (
                "Vulkanus",
                2415020.0,
                77.445895,
                0.0,
                0.0,
                0.0,
                0.0,
                55.397715,
                0.0015069325,
            ),
        ),
        (
            POSEIDON_KEPLERIAN_ELEMENTS,
            PoseidonKeplerianElements,
            (
                "Poseidon",
                2415020.0,
                83.666307,
                0.0,
                0.0,
                0.0,
                0.0,
                166.140256,
                0.0013256078,
            ),
        ),
    ]
    for elements, expected_type, expected_values in expected_keplerian:
        assert type(elements) is expected_type
        assert tuple(vars(elements).values()) == expected_values

    # These are import-compatibility objects; the independently gated runtime
    # registry remains separate and unchanged.
    assert CUPIDO_KEPLERIAN_ELEMENTS is not URANIAN_KEPLERIAN_ELEMENTS[CUPIDO]


def test_rc7_legacy_public_names_are_real_module_level_definitions() -> None:
    import libephemeris.hypothetical as hypothetical

    module_path = Path(hypothetical.__file__)
    tree = ast.parse(module_path.read_text(encoding="utf-8"))
    defined_names: set[str] = set()
    for statement in tree.body:
        if isinstance(statement, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            defined_names.add(statement.name)
        elif isinstance(statement, ast.Assign):
            defined_names.update(
                target.id
                for target in statement.targets
                if isinstance(target, ast.Name)
            )
        elif isinstance(statement, ast.AnnAssign) and isinstance(
            statement.target, ast.Name
        ):
            defined_names.add(statement.target.id)

    assert _RC7_LEGACY_PUBLIC_NAMES <= defined_names


def test_harrington_position_is_finite_and_continuous() -> None:
    position = calc_fictitious_position(HARRINGTON, J2000)
    later = calc_hypothetical_position(HARRINGTON, J2000 + 1.0)

    assert len(position) == 6
    assert all(math.isfinite(value) for value in position)
    assert 0.0 <= position[0] < 360.0
    assert -90.0 <= position[1] <= 90.0
    assert position[2] > 0.0
    delta = ((later[0] - position[0] + 180.0) % 360.0) - 180.0
    assert abs(delta) < 1.0


@pytest.mark.parametrize("body_id", range(CUPIDO, WALDEMATH + 1))
def test_all_restored_ids_return_finite_states(
    body_id: int,
) -> None:
    assert is_hypothetical_body(body_id)
    state = calc_hypothetical_position(body_id, J2000)
    assert len(state) == 6
    assert all(math.isfinite(value) for value in state)
    assert 0.0 <= state[0] < 360.0
    assert state[2] > 0.0


@pytest.mark.parametrize(
    ("function", "body_id"),
    [
        (calc_cupido, CUPIDO),
        (calc_hades, HADES),
        (calc_zeus, ZEUS),
        (calc_kronos, KRONOS),
        (calc_apollon, APOLLON),
        (calc_admetos, ADMETOS),
        (calc_vulkanus, VULKANUS),
        (calc_poseidon, POSEIDON),
        (calc_transpluto, ISIS),
        (calc_transpluto_position, ISIS),
        (calc_vulcan, VULCAN),
        (calc_white_moon_position, WHITE_MOON),
        (calc_proserpina, PROSERPINA),
        (calc_waldemath, WALDEMATH),
        (calc_waldemath_position, WALDEMATH),
        (calc_planet_x_lowell, 53),
        (calc_planet_x_pickering, 54),
    ],
)
def test_named_entry_points_compute(function: object, body_id: int) -> None:
    assert is_hypothetical_body(body_id)
    result = function(J2000)  # type: ignore[operator]
    assert len(result) == 6
    assert all(math.isfinite(value) for value in result)


@pytest.mark.parametrize("function", [calc_uranian_planet, calc_uranian_position])
def test_generic_uranian_entry_points_compute(function: object) -> None:
    result = function(CUPIDO, J2000)  # type: ignore[operator]
    assert len(result) == 6


@pytest.mark.parametrize("function", [calc_uranian_planet, calc_uranian_position])
def test_generic_uranian_entry_points_preserve_body_id_on_error(
    function: object,
) -> None:
    with pytest.raises(UnknownBodyError):
        function(999, J2000)  # type: ignore[operator]


def test_uranian_longitude_matches_full_position() -> None:
    assert calc_uranian_longitude(CUPIDO, J2000) == calc_cupido(J2000)[0]


def test_missing_fictitious_elements_error_carries_body_id() -> None:
    with pytest.raises(UnknownBodyError) as raised:
        calc_fictitious_position(999, J2000)
    assert raised.value.body_id == 999


@pytest.mark.parametrize("jd", [2415020.0, J2000, 2488070.0])
def test_white_moon_pins_rc7_seven_year_convention(jd: float) -> None:
    centuries = (jd - J2000) / 36525.0
    expected = (
        (242.2205555 + 5143.5418158 * centuries) % 360.0,
        0.0,
        0.05280098949,
        5143.5418158 / 36525.0,
        0.0,
        0.0,
    )
    assert calc_white_moon_position(jd) == expected
    assert calc_white_moon_position(jd, use_true_lilith=True) == expected


def test_recognised_names_remain_available_for_api_compatibility() -> None:
    names = list_hypothetical_bodies()
    assert names[CUPIDO] == "Cupido"
    assert names[HARRINGTON] == "Harrington"
    assert names[NIBIRU] == "Nibiru"
    assert names[WALDEMATH] == "Waldemath"


def test_get_orbital_body_by_name_is_case_insensitive() -> None:
    rows = load_bundled_fictitious_orbits()
    body = get_orbital_body_by_name(rows, "hArRiNgToN")
    assert body is not None
    assert body.name == "Harrington"
    assert get_orbital_body_by_name(rows, "missing") is None


def test_synthetic_user_orbit_propagation_is_independent_of_builtins() -> None:
    elements = OrbitalElements(
        name="Synthetic circular test orbit",
        epoch_jd=J2000,
        equinox_jd=J2000,
        equinox_is_jdate=False,
        mean_anomaly=TPolynomial(0.0),
        semi_axis=2.0,
        eccentricity=TPolynomial(0.0),
        arg_perihelion=TPolynomial(0.0),
        asc_node=TPolynomial(0.0),
        inclination=TPolynomial(0.0),
    )
    at_epoch = calc_orbital_position(elements, J2000)
    later = calc_orbital_position(elements, J2000 + 10.0)

    assert at_epoch[:3] == pytest.approx((0.0, 0.0, 2.0), abs=1e-12)
    assert 0.0 < later[0] < 10.0
    assert later[2] == pytest.approx(2.0, abs=1e-12)


def test_reserved_keplerian_table_accepts_synthetic_elements(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    body_id = 999
    synthetic = HypotheticalElements(
        name="Synthetic",
        epoch=J2000,
        a=3.0,
        e=0.1,
        i=5.0,
        omega=20.0,
        Omega=30.0,
        M0=40.0,
        n=0.1,
    )
    monkeypatch.setitem(HYPOTHETICAL_ELEMENTS, body_id, synthetic)

    full = _calc_keplerian_position(body_id, J2000)
    raw = _calc_keplerian_position_raw(body_id, J2000)
    assert full[:3] == pytest.approx(raw)
    assert all(math.isfinite(value) for value in full)


@pytest.mark.parametrize("eccentricity", [0.0, 0.2, 0.9])
def test_kepler_solver_satisfies_equation(eccentricity: float) -> None:
    mean_anomaly = 1.2
    eccentric_anomaly = _solve_kepler_equation(mean_anomaly, eccentricity)
    assert eccentric_anomaly - eccentricity * math.sin(
        eccentric_anomaly
    ) == pytest.approx(mean_anomaly, abs=1e-10)


def test_non_hypothetical_id_is_rejected() -> None:
    with pytest.raises(ValueError):
        calc_hypothetical_position(0, J2000)
