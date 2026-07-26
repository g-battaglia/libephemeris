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
    LOWELL_PLANET_X_ELEMENTS,
    NEPTUNE_ADAMS,
    NEPTUNE_LEVERRIER,
    NIBIRU,
    PICKERING_PLANET_X_ELEMENTS,
    POSEIDON,
    POSEIDON_KEPLERIAN_ELEMENTS,
    PLUTO_LOWELL,
    PLUTO_PICKERING,
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
    calc_white_moon_position,
    calc_zeus,
    get_orbital_body_by_name,
    is_hypothetical_body,
    list_hypothetical_bodies,
    load_bundled_fictitious_orbits,
)


J2000 = 2451545.0

_SUPPORTED_IDS = {
    *range(CUPIDO, POSEIDON + 1),
    ISIS,
    HARRINGTON,
    NEPTUNE_LEVERRIER,
    NEPTUNE_ADAMS,
    PLUTO_LOWELL,
    PLUTO_PICKERING,
    VULCAN,
    WHITE_MOON,
    PROSERPINA,
    WALDEMATH,
}
_UNSUPPORTED_IDS = set(range(CUPIDO, WALDEMATH + 1)) - _SUPPORTED_IDS

# Reviewed literal transcriptions of one identified primary publication.
_PRIMARY_TRANSCRIPTION_IDS = {
    *range(CUPIDO, POSEIDON + 1),
    HARRINGTON,
    NEPTUNE_LEVERRIER,
    NEPTUNE_ADAMS,
    PLUTO_LOWELL,
    PLUTO_PICKERING,
}
# Published conventions realized from documented public sources.
_PUBLISHED_MODEL_IDS = {ISIS, VULCAN, WHITE_MOON, PROSERPINA, WALDEMATH}

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


def test_builtin_tables_contain_only_sourced_models() -> None:
    assert set(URANIAN_ELEMENTS) == set(range(CUPIDO, POSEIDON + 1))
    assert set(URANIAN_KEPLERIAN_ELEMENTS) == set(range(CUPIDO, POSEIDON + 1))
    assert set(HYPOTHETICAL_ELEMENTS) == {ISIS, PROSERPINA}
    assert type(TRANSPLUTO_KEPLERIAN_ELEMENTS) is TransplutoKeplerianElements
    assert TRANSPLUTO_KEPLERIAN_ELEMENTS.name == "Transpluto"
    assert HYPOTHETICAL_ELEMENTS[ISIS] is TRANSPLUTO_KEPLERIAN_ELEMENTS
    assert set(FICTITIOUS_ORBITAL_ELEMENTS) == {
        HARRINGTON,
        NEPTUNE_LEVERRIER,
        NEPTUNE_ADAMS,
        PLUTO_LOWELL,
        PLUTO_PICKERING,
    }
    assert VULCAN_ELEMENTS.name == "Vulcan"
    assert WALDEMATH_ELEMENTS.name == "Waldemath"
    assert set(HYPOTHETICAL_PROVENANCE) == set(range(CUPIDO, WALDEMATH + 1))
    assert _PRIMARY_TRANSCRIPTION_IDS | _PUBLISHED_MODEL_IDS == _SUPPORTED_IDS
    assert all(
        HYPOTHETICAL_PROVENANCE[body_id][0] == "primary-transcription"
        for body_id in _PRIMARY_TRANSCRIPTION_IDS
    )
    assert all(
        HYPOTHETICAL_PROVENANCE[body_id][0] == "published-model"
        for body_id in _PUBLISHED_MODEL_IDS
    )
    assert all(
        HYPOTHETICAL_PROVENANCE[body_id][0] == "unsupported"
        for body_id in _UNSUPPORTED_IDS
    )

    rows = load_bundled_fictitious_orbits()
    assert len(rows) == 12
    assert {row.name for row in rows} == {
        "Cupido",
        "Hades",
        "Zeus",
        "Kronos",
        "Apollon",
        "Admetos",
        "Vulkanus",
        "Poseidon",
        "Harrington",
        "Leverrier-Neptune",
        "Adams-Neptune",
        "Lowell-Pluto",
    }
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
    assert set(URANIAN_KEPLERIAN_ELEMENTS) == set(range(CUPIDO, POSEIDON + 1))
    assert set(HYPOTHETICAL_ELEMENTS) == {ISIS, PROSERPINA}


def test_legacy_public_element_objects_reflect_reviewed_or_unavailable_state() -> None:
    assert set(URANIAN_ELEMENTS) == set(range(CUPIDO, POSEIDON + 1))
    assert UranianElements.__name__ == "UranianElements"

    reviewed_keplerian = [
        (CUPIDO_KEPLERIAN_ELEMENTS, CupidoKeplerianElements, "Cupido"),
        (HADES_KEPLERIAN_ELEMENTS, HadesKeplerianElements, "Hades"),
        (ZEUS_KEPLERIAN_ELEMENTS, ZeusKeplerianElements, "Zeus"),
        (KRONOS_KEPLERIAN_ELEMENTS, KronosKeplerianElements, "Kronos"),
        (APOLLON_KEPLERIAN_ELEMENTS, ApollonKeplerianElements, "Apollon"),
        (ADMETOS_KEPLERIAN_ELEMENTS, AdmetosKeplerianElements, "Admetos"),
        (VULKANUS_KEPLERIAN_ELEMENTS, VulkanusKeplerianElements, "Vulkanus"),
        (POSEIDON_KEPLERIAN_ELEMENTS, PoseidonKeplerianElements, "Poseidon"),
    ]
    for elements, expected_type, expected_name in reviewed_keplerian:
        assert type(elements) is expected_type
        values = vars(elements).copy()
        assert values.pop("name") == expected_name
        assert values
        assert all(math.isfinite(value) for value in values.values())

    for finite_elements in (
        LOWELL_PLANET_X_ELEMENTS,
        TRANSPLUTO_KEPLERIAN_ELEMENTS,
        PICKERING_PLANET_X_ELEMENTS,
        VULCAN_ELEMENTS,
    ):
        assert all(
            math.isfinite(value)
            for name, value in vars(finite_elements).items()
            if name != "name"
        )
    # Waldemath carries the published uniform model (Sepharial 1918 rate and
    # anchor, Waltemath 1898 distance): every field is finite.
    values = vars(WALDEMATH_ELEMENTS).copy()
    assert isinstance(values.pop("name"), str)
    assert values
    assert all(math.isfinite(value) for value in values.values())


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


@pytest.mark.parametrize("body_id", sorted(_SUPPORTED_IDS))
def test_reviewed_position_is_finite_and_continuous(body_id: int) -> None:
    position = calc_hypothetical_position(body_id, J2000)
    later = calc_hypothetical_position(body_id, J2000 + 1.0)

    assert len(position) == 6
    assert all(math.isfinite(value) for value in position)
    assert 0.0 <= position[0] < 360.0
    assert -90.0 <= position[1] <= 90.0
    assert position[2] > 0.0
    delta = ((later[0] - position[0] + 180.0) % 360.0) - 180.0
    # Weston's Vulcan is the single fast body (~19.4 deg/day sidereal); all
    # other supported hypotheticals move well under a degree per day.
    # Vulcan races at ~19 deg/day, the Waldemath Dark Moon at its published
    # 3 deg/day; every other model moves below 1 deg/day.
    limit = 25.0 if body_id == VULCAN else 3.5 if body_id == WALDEMATH else 1.0
    assert abs(delta) < limit


@pytest.mark.parametrize(
    "body_id",
    sorted(_UNSUPPORTED_IDS),
)
def test_recognised_but_unverified_ids_raise_unknown_body(body_id: int) -> None:
    assert is_hypothetical_body(body_id)
    with pytest.raises(UnknownBodyError) as raised:
        calc_hypothetical_position(body_id, J2000)
    assert raised.value.body_id == body_id


def test_waldemath_matches_published_uniform_model() -> None:
    # Sepharial 1918: uniform 3 deg/day tropical longitude anchored at a
    # printed Lilith-Sun conjunction (here Waltemath's predicted 1898 Feb 2
    # transit, 00:00 GMT); Waltemath 1898: mean distance 1.03 million km.
    from libephemeris.hypothetical import (
        _WALDEMATH_ANCHOR_APPARENT_LON_DEG,
        _WALDEMATH_ANCHOR_JD_UT,
        _WALDEMATH_DISTANCE_AU,
    )
    from libephemeris.time_utils import deltat

    anchor_tt = _WALDEMATH_ANCHOR_JD_UT + deltat(_WALDEMATH_ANCHOR_JD_UT)
    state = calc_waldemath(anchor_tt)
    # Rate from the printed 177-day synodic period and the mean tropical
    # year; Sepharial's 126-year same-day return is exactly 260 synodic
    # revolutions, and the implied sidereal period is Waltemath's 119 days.
    rate = 360.0 / 177.0 + 360.0 / 365.2422
    assert 126.0 * 365.2422 / 177.0 == pytest.approx(260.0, abs=0.01)
    assert 360.0 / rate == pytest.approx(119.0, abs=0.3)
    later = calc_waldemath(anchor_tt + 1.0)
    advance = (later[0] - state[0]) % 360.0
    assert advance == pytest.approx(rate, abs=1e-12)
    assert state[1] == 0.0 and state[4] == 0.0 and state[5] == 0.0
    assert state[2] == pytest.approx(1.03e6 / 149_597_870.7, abs=1e-15)
    assert state[3] == pytest.approx(rate, abs=1e-12)
    assert _WALDEMATH_DISTANCE_AU == pytest.approx(0.00688512, abs=1e-8)
    assert 0.0 <= _WALDEMATH_ANCHOR_APPARENT_LON_DEG < 360.0


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
        (calc_planet_x_lowell, PLUTO_LOWELL),
        (calc_planet_x_pickering, PLUTO_PICKERING),
        (calc_vulcan, VULCAN),
        (calc_white_moon_position, WHITE_MOON),
        (calc_proserpina, PROSERPINA),
    ],
)
def test_reviewed_named_entry_points_return_finite_state(
    function: object, body_id: int
) -> None:
    assert is_hypothetical_body(body_id)
    state = function(J2000)  # type: ignore[operator]
    assert len(state) == 6
    assert all(math.isfinite(value) for value in state)


@pytest.mark.parametrize("function", [calc_uranian_planet, calc_uranian_position])
def test_generic_uranian_entry_points_calculate_reviewed_body(function: object) -> None:
    state = function(CUPIDO, J2000)  # type: ignore[operator]
    assert len(state) == 6
    assert all(math.isfinite(value) for value in state)


@pytest.mark.parametrize("function", [calc_uranian_planet, calc_uranian_position])
def test_generic_uranian_entry_points_preserve_body_id_on_error(
    function: object,
) -> None:
    with pytest.raises(UnknownBodyError):
        function(999, J2000)  # type: ignore[operator]


def test_uranian_longitude_calculates_reviewed_body() -> None:
    longitude = calc_uranian_longitude(CUPIDO, J2000)
    assert 0.0 <= longitude < 360.0


def test_missing_fictitious_elements_error_carries_body_id() -> None:
    with pytest.raises(UnknownBodyError) as raised:
        calc_fictitious_position(999, J2000)
    assert raised.value.body_id == 999


@pytest.mark.parametrize("use_true_lilith", [False, True])
def test_white_moon_matches_published_seven_year_model(
    use_true_lilith: bool,
) -> None:
    # Velichko--Larin (2007), p. 17 defines uniform zodiacal motion with a
    # seven-year cycle.  Pages 20 and 45 print January 1800 as 6 deg 08 arcmin
    # Taurus and January 2000 as 2 deg 18 arcmin Sagittarius.  LibEphemeris
    # maps monthly rows to the first civil day at 00:00 TT; the common January
    # label means any fixed alternative sample instant would cancel from the
    # rate's 200-year time difference.
    early_jd = 2378496.5
    early_longitude = 30.0 + 6.0 + 8.0 / 60.0
    source_jd = 2451544.5
    source_longitude = 240.0 + 2.0 + 18.0 / 60.0
    elapsed_days = source_jd - early_jd
    displacement = 28.0 * 360.0 + source_longitude - early_longitude
    daily_motion = displacement / elapsed_days
    period_days = 360.0 / daily_motion

    source_state = calc_white_moon_position(source_jd, use_true_lilith=use_true_lilith)
    assert source_state[0] == pytest.approx(source_longitude, abs=1e-12)
    early_state = calc_white_moon_position(early_jd, use_true_lilith=use_true_lilith)
    assert early_state[0] == pytest.approx(early_longitude, abs=1e-12)

    state = calc_white_moon_position(J2000, use_true_lilith=use_true_lilith)
    assert state[0] == pytest.approx(
        source_longitude + daily_motion * (J2000 - source_jd), abs=1e-12
    )
    assert state[1] == 0.0
    assert state[2] == pytest.approx(
        (3.986004e14 * (period_days * 86_400.0 / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0)
        / 149_597_870_700.0,
        abs=1e-15,
    )
    assert state[3] == pytest.approx(daily_motion, abs=1e-15)
    assert state[4:] == (0.0, 0.0)

    # These rows were not used to derive rate or phase.  Their maximum residual
    # is below half an arcminute, as required for nearest-minute printed data:
    # March 1879 is on p. 29; December 2000 and January 2007 are on p. 45.
    for jd, published in (
        (2407409.5, 120.0 + 27.0 + 29.0 / 60.0),
        (2451879.5, 270.0 + 19.0 + 28.0 / 60.0),
        (2454101.5, 240.0 + 2.0 + 22.0 / 60.0),
    ):
        calculated = calc_white_moon_position(jd)[0]
        difference = abs((calculated - published + 180.0) % 360.0 - 180.0)
        assert difference <= 0.5 / 60.0


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
