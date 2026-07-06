"""Mean orbital elements of the major planets (mean equinox of date).

Polynomial expressions from Meeus, "Astronomical Algorithms" (2nd ed.,
1998), Table 31.A, which is derived from Simon et al. (1994),
"Numerical expressions for precession formulae and mean elements for
the Moon and planets", A&A 282, 663-683.

Each element is evaluated as ``c0 + c1*T + c2*T**2 + c3*T**3`` with T in
Julian centuries of 36525 days from J2000.0 TT, and refers to the mean
ecliptic and equinox **of date**.

These mean elements feed the NODBIT_MEAN branch of ``nod_aps`` /
``nod_aps_ut``: the reference API computes mean nodes and apsides of the
planets from this published element set (verified black-box to <0.1").
"""

from __future__ import annotations

from .constants import (
    EARTH,
    JUPITER,
    MARS,
    MERCURY,
    NEPTUNE,
    SATURN,
    URANUS,
    VENUS,
)

# Per planet: a (AU, polynomial), e, i (deg), Omega (deg), pi (deg).
# i / Omega / pi refer to the mean ecliptic and equinox of date.
# Earth has no node (its orbit defines the ecliptic): Omega is None.
_MEAN_ELEMENTS: dict[int, dict[str, tuple[float, ...] | None]] = {
    MERCURY: {
        "a": (0.387098310,),
        "e": (0.20563175, 0.000020407, -0.0000000283, -0.00000000018),
        "i": (7.004986, 0.0018215, -0.00001810, 0.000000056),
        "Omega": (48.330893, 1.1861890, 0.00017587, 0.000000211),
        "pi": (77.456119, 1.5564775, 0.00029589, 0.000000056),
    },
    VENUS: {
        "a": (0.723329820,),
        "e": (0.00677188, -0.000047766, 0.0000000975, 0.00000000044),
        "i": (3.394662, 0.0010037, -0.00000088, -0.000000007),
        "Omega": (76.679920, 0.9011190, 0.00040665, -0.000000080),
        "pi": (131.563707, 1.4022188, -0.00107337, -0.000005315),
    },
    EARTH: {
        "a": (1.000001018,),
        "e": (0.01670862, -0.000042037, -0.0000001236, 0.00000000004),
        "i": (0.0,),
        "Omega": None,
        "pi": (102.937348, 1.7195269, 0.00045962, 0.000000499),
    },
    MARS: {
        "a": (1.523679342,),
        "e": (0.09340062, 0.000090483, -0.0000000806, -0.00000000035),
        "i": (1.849726, -0.0006010, 0.00001276, -0.000000006),
        "Omega": (49.558093, 0.7720923, 0.00001605, 0.000002325),
        "pi": (336.060234, 1.8410331, 0.00013515, 0.000000318),
    },
    JUPITER: {
        "a": (5.202603191, 0.0000001913),
        "e": (0.04849485, 0.000163244, -0.0000004719, -0.00000000197),
        "i": (1.303270, -0.0054966, 0.00000465, -0.000000004),
        "Omega": (100.464441, 1.0209550, 0.00040117, 0.000000569),
        "pi": (14.331309, 1.6126668, 0.00103127, -0.000004569),
    },
    SATURN: {
        "a": (9.554909596, -0.0000021389),
        "e": (0.05550862, -0.000346818, -0.0000006456, 0.00000000338),
        "i": (2.488878, -0.0037363, -0.00001516, 0.000000089),
        "Omega": (113.665524, 0.8770979, -0.00012067, -0.000002380),
        "pi": (93.056787, 1.9637694, 0.00083757, 0.000004899),
    },
    URANUS: {
        "a": (19.218446062, -0.0000000372, 0.00000000098),
        "e": (0.04629590, -0.000027337, 0.0000000790, 0.00000000025),
        "i": (0.773196, 0.0007744, 0.00003749, -0.000000092),
        "Omega": (74.005947, 0.5211258, 0.00133982, 0.000018516),
        "pi": (173.005159, 1.4863784, 0.00021450, 0.000000433),
    },
    NEPTUNE: {
        "a": (30.110386869, -0.0000001663, 0.00000000069),
        "e": (0.00898809, 0.000006408, -0.0000000008, -0.00000000005),
        "i": (1.769952, -0.0093082, -0.00000708, 0.000000028),
        "Omega": (131.784057, 1.1022057, 0.00026006, -0.000000636),
        "pi": (48.123691, 1.4262677, 0.00037918, -0.000000003),
    },
}

_J2000_JD = 2451545.0


def _poly(coeffs: tuple[float, ...], t: float) -> float:
    """Evaluate a polynomial with coefficients ordered by ascending power."""
    result = 0.0
    for c in reversed(coeffs):
        result = result * t + c
    return result


def mean_orbital_elements(
    ipl: int, jd_tt: float
) -> tuple[float, float, float, float, float] | None:
    """Mean orbital elements of a planet at the mean equinox of date.

    Args:
        ipl: Planet ID (MERCURY..NEPTUNE or EARTH).
        jd_tt: Julian Day (TT).

    Returns:
        Tuple ``(a_au, e, i_deg, Omega_deg, pi_deg)`` where ``Omega_deg``
        is the longitude of the ascending node and ``pi_deg`` the
        longitude of perihelion, both referred to the mean ecliptic and
        equinox of date. For Earth (no node) ``Omega_deg`` is 0.0 and
        ``i_deg`` is 0.0. Returns None for bodies without mean elements
        (e.g. Pluto), for which callers fall back to osculating elements.
    """
    el = _MEAN_ELEMENTS.get(ipl)
    if el is None:
        return None
    t = (jd_tt - _J2000_JD) / 36525.0
    a = _poly(el["a"], t)  # type: ignore[arg-type]
    e = _poly(el["e"], t)  # type: ignore[arg-type]
    i_deg = _poly(el["i"], t)  # type: ignore[arg-type]
    omega_coeffs = el["Omega"]
    omega_deg = _poly(omega_coeffs, t) % 360.0 if omega_coeffs is not None else 0.0
    pi_deg = _poly(el["pi"], t) % 360.0  # type: ignore[arg-type]
    return (a, e, i_deg, omega_deg, pi_deg)
