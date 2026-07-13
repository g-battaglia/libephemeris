"""Branch-focused clean-room tests for :mod:`libephemeris.hypothetical`."""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from libephemeris.exceptions import UnknownBodyError
from libephemeris.hypothetical import (
    CUPIDO,
    HARRINGTON,
    HYPOTHETICAL_ELEMENTS,
    ISIS,
    OrbitalElements,
    TPolynomial,
    UranianKeplerianElements,
    _calc_keplerian_position,
    _calc_keplerian_position_raw,
    _calc_orbital_position_raw,
    _calc_transpluto_raw,
    _keplerian_to_ecliptic_j2000,
    _parse_epoch_or_equinox,
    _parse_fictitious_orbits_csv,
    _parse_orbital_elements_line,
    _parse_t_polynomial,
    calc_fictitious_position,
    calc_hypothetical_position,
    get_bundled_fictitious_orbits_path,
    get_hypothetical_name,
    get_orbital_body_by_name,
    is_hypothetical_body,
    load_bundled_fictitious_orbits,
    parse_orbital_elements,
)


J2000 = 2451545.0


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("J2000", (2451545.0, False)),
        ("J1900", (2415020.0, False)),
        ("B1950", (2433282.42345905, False)),
        ("JDATE", (None, True)),
        ("2450000.5", (2450000.5, False)),
    ],
)
def test_parse_epoch_variants(text: str, expected: tuple[float | None, bool]) -> None:
    assert _parse_epoch_or_equinox(text) == expected


def test_parse_epoch_rejects_unknown_text() -> None:
    with pytest.raises(ValueError, match="Cannot parse"):
        _parse_epoch_or_equinox("not-an-epoch")


@pytest.mark.parametrize(
    ("expression", "constant", "linear"),
    [
        ("12.5", 12.5, 0.0),
        ("12.5 + 3.0 * T", 12.5, 3.0),
        ("12.5 - 3.0*T", 12.5, -3.0),
        ("3*T + 12", 12.0, 3.0),
        ("3*T - 12", -12.0, 3.0),
        ("T + 2", 2.0, 1.0),
    ],
)
def test_parse_t_polynomial_variants(
    expression: str, constant: float, linear: float
) -> None:
    value = _parse_t_polynomial(expression)
    assert value.constant == constant
    assert value.linear == linear


def test_parse_t_polynomial_rejects_non_numeric_text() -> None:
    with pytest.raises(ValueError):
        _parse_t_polynomial("invalid")


def test_parse_custom_file_and_inline_comment(tmp_path: Path) -> None:
    path = tmp_path / "orbits.txt"
    path.write_text(
        "# synthetic user data\n"
        "J2000,J2000,10 + 2*T,2.5,0.1,20,30,4,Synthetic # comment\n",
        encoding="utf-8",
    )
    rows = parse_orbital_elements(path)
    assert len(rows) == 1
    assert rows[0].name == "Synthetic"
    assert rows[0].mean_anomaly.linear == 2.0


def test_parse_custom_file_reports_line_number(tmp_path: Path) -> None:
    path = tmp_path / "bad.txt"
    path.write_text("J2000,too,few\n", encoding="utf-8")
    with pytest.raises(ValueError, match="line 1"):
        parse_orbital_elements(path)


def test_parse_custom_file_missing(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        parse_orbital_elements(tmp_path / "missing.txt")


def test_parse_legacy_line_geocentric_marker() -> None:
    row = _parse_orbital_elements_line("J2000,J2000,0,1,0,0,0,0,Synthetic, geo", 7)
    assert row is not None
    assert row.name == "Synthetic"
    assert row.is_geocentric is True
    assert row.line_number == 7


def test_parse_csv_validation(tmp_path: Path) -> None:
    path = tmp_path / "bad.csv"
    path.write_text("too,few,columns\n", encoding="utf-8")
    with pytest.raises(ValueError, match="at least 10"):
        _parse_fictitious_orbits_csv(path)

    path.write_text(
        "Synthetic,J2000,J2000,0,not-a-number,0,0,0,0,0\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="semi-major axis"):
        _parse_fictitious_orbits_csv(path)

    path.write_text(
        "Synthetic,J2000,J2000,0,1,0,0,0,0,not-a-flag\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="geocentric flag"):
        _parse_fictitious_orbits_csv(path)


def test_bundled_dataset_path_and_lookup() -> None:
    path = get_bundled_fictitious_orbits_path()
    assert path.name == "fictitious_orbits.csv"
    rows = load_bundled_fictitious_orbits()
    harrington = get_orbital_body_by_name(rows, "HARRINGTON")
    assert harrington is not None
    assert harrington.name == "Harrington"


def test_body_predicates_and_unknown_name() -> None:
    assert is_hypothetical_body(CUPIDO)
    assert not is_hypothetical_body(0)
    assert get_hypothetical_name(HARRINGTON) == "Harrington"
    assert get_hypothetical_name(999) == "Unknown (999)"


def test_fictitious_position_rejects_nonclassical_id_with_body_id() -> None:
    with pytest.raises(UnknownBodyError) as raised:
        calc_fictitious_position(ISIS, J2000)
    assert raised.value.body_id == ISIS


def test_transpluto_raw_computes_restored_elements() -> None:
    position = _calc_transpluto_raw(J2000)
    assert len(position) == 3
    assert all(math.isfinite(value) for value in position)


def test_dispatch_computes_restored_uranian_id() -> None:
    position = calc_hypothetical_position(CUPIDO, J2000)
    assert len(position) == 6
    assert all(math.isfinite(value) for value in position)


def test_synthetic_keplerian_rotation_and_precession() -> None:
    elements = UranianKeplerianElements(
        name="Synthetic",
        epoch=J2000,
        equinox_jd=2433282.42345905,
        a=4.0,
        e=0.2,
        i=10.0,
        omega=20.0,
        Omega=30.0,
        M0=40.0,
        n=0.05,
    )
    position = _keplerian_to_ecliptic_j2000(elements, J2000)
    assert all(math.isfinite(value) for value in position)
    assert 0.0 <= position[0] < 360.0
    assert position[2] > 0.0


def test_synthetic_reserved_table_velocity_and_raw(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    body_id = 999
    elements = __import__(
        "libephemeris.hypothetical", fromlist=["HypotheticalElements"]
    ).HypotheticalElements(
        name="Synthetic",
        epoch=J2000,
        a=1.5,
        e=0.0,
        i=0.0,
        omega=0.0,
        Omega=0.0,
        M0=359.999,
        n=1.0,
    )
    monkeypatch.setitem(HYPOTHETICAL_ELEMENTS, body_id, elements)
    full = _calc_keplerian_position(body_id, J2000)
    raw = _calc_keplerian_position_raw(body_id, J2000)
    assert full[:3] == pytest.approx(raw)
    assert abs(full[3]) < 2.0


def test_orbital_raw_handles_circular_and_elliptic_paths() -> None:
    circular = OrbitalElements(
        name="Circular",
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
    elliptic = OrbitalElements(
        name="Elliptic",
        epoch_jd=J2000,
        equinox_jd=J2000,
        equinox_is_jdate=False,
        mean_anomaly=TPolynomial(30.0),
        semi_axis=2.0,
        eccentricity=TPolynomial(0.2),
        arg_perihelion=TPolynomial(10.0),
        asc_node=TPolynomial(5.0),
        inclination=TPolynomial(3.0),
    )
    assert _calc_orbital_position_raw(circular, J2000) == pytest.approx((0.0, 0.0, 2.0))
    assert all(
        math.isfinite(value) for value in _calc_orbital_position_raw(elliptic, J2000)
    )
