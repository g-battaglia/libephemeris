"""Fixed-epoch sidereal modes (SIDM_J2000 / J1900 / B1950) as frame requests.

The reference API implements these three ayanamsha modes as *frame*
transformations (position expressed on the mean ecliptic/equinox of the
mode's epoch t0), not as a scalar longitude offset. Expected values below
were frozen from the black-box reference oracle (pyswisseph 2.10.3.2 on the
bundled reference ephemeris files).

Covers: planets (positions, speeds, retflag), both fixed-star families,
houses (plain, near-t0 quirk window, exact-t0, Sunshine collapse, Aries,
Gauquelin), and the mode-18 identity with FLG_J2000|FLG_NONUT.
"""

from __future__ import annotations

import pytest

import libephemeris as le
from libephemeris.constants import (
    FLG_EQUATORIAL,
    FLG_J2000,
    FLG_NONUT,
    FLG_SIDEREAL,
    FLG_SPEED,
    FLG_SWIEPH,
)


FROZEN = {
    "planets": {
        "18|2448000.5|1": [
            [
                307.084651906894,
                -0.5531095698839161,
                0.0025916489267754548,
                12.815424459164035,
                1.1261880086757028,
                -3.876382623306744e-05,
            ],
            65858,
        ],
        "18|2448000.5|5": [
            [
                95.17967779029192,
                0.10563447585538353,
                5.510347028590518,
                0.1465938644092365,
                0.0015150740285322913,
                0.01482691514777773,
            ],
            65858,
        ],
        "18|2415020.0|1": [
            [
                266.6888069313063,
                0.4484163699148131,
                0.0024765268649707996,
                14.151230053754094,
                1.2966981297747404,
                -2.9912792278011017e-05,
            ],
            65858,
        ],
        "18|2415020.0|5": [
            [
                242.42982687008384,
                0.8021065348883112,
                6.118255556219077,
                0.1960990420034192,
                0.00022370569102049434,
                -0.010330449015371088,
            ],
            65858,
        ],
        "19|2448000.5|1": [
            [
                305.68798743222914,
                -0.5433961764431808,
                0.0025916489267754556,
                12.815575259725001,
                1.1242338701598436,
                -3.876382623306738e-05,
            ],
            65858,
        ],
        "19|2448000.5|5": [
            [
                93.78309296412898,
                0.09277112908214999,
                5.510347028590518,
                0.14659374693745103,
                0.0015209154119477798,
                0.014826915147777727,
            ],
            65858,
        ],
        "19|2415020.0|1": [
            [
                265.2922289069144,
                0.4614758288622709,
                0.0024765268649707996,
                14.151263785525602,
                1.296609519023142,
                -2.9912792278011027e-05,
            ],
            65858,
        ],
        "19|2415020.0|5": [
            [
                241.03317492549246,
                0.8141602140737692,
                6.118255556219075,
                0.19609960348726474,
                0.0002409505306411119,
                -0.010330449015371095,
            ],
            65858,
        ],
        "20|2448000.5|1": [
            [
                306.3862397402475,
                -0.5482637599864082,
                0.002591648926775456,
                12.815500044745743,
                1.125208998359325,
                -3.87638262330674e-05,
            ],
            65858,
        ],
        "20|2448000.5|5": [
            [
                94.48130562941436,
                0.09920745090316069,
                5.510347028590519,
                0.14659380445154846,
                0.0015180283134106091,
                0.014826915147777725,
            ],
            65858,
        ],
        "20|2415020.0|1": [
            [
                265.99043835645705,
                0.45494344051298846,
                0.0024765268649708004,
                14.151247132003101,
                1.2966504476801437,
                -2.991279227801103e-05,
            ],
            65858,
        ],
        "20|2415020.0|5": [
            [
                241.73142156057858,
                0.8081365684380292,
                6.118255556219075,
                0.1960993220588384,
                0.00023228169813495708,
                -0.010330449015371095,
            ],
            65858,
        ],
    },
    "stars": {
        "18|Aldebaran|v1": [
            [
                69.78484920503931,
                -5.467151726372493,
                4214536.311398931,
                -6.411984493228058e-05,
                6.787012539446518e-06,
                0.0442609116998047,
            ],
            65858,
        ],
        "18|Aldebaran|v2": [
            [
                69.78484920503931,
                -5.467151726372493,
                4214536.311398931,
                -6.411984493228058e-05,
                6.787012539446518e-06,
                0.0442609116998047,
            ],
            65858,
        ],
        "18|Polaris|v1": [
            [
                88.56060048741556,
                66.10610288763652,
                27356108.917828392,
                -0.00020779629300933587,
                -4.4381422576845234e-05,
                -0.004206197905780325,
            ],
            65858,
        ],
        "18|Polaris|v2": [
            [
                88.56060048741556,
                66.10610288763652,
                27356108.917828392,
                -0.00020779629300933587,
                -4.4381422576845234e-05,
                -0.004206197905780325,
            ],
            65858,
        ],
        "19|Aldebaran|v1": [
            [
                68.38793736030283,
                -5.47975124027902,
                4214536.311398932,
                -6.412078301235929e-05,
                6.790877998652985e-06,
                0.04426091169980453,
            ],
            65858,
        ],
        "19|Aldebaran|v2": [
            [
                68.38793736030283,
                -5.47975124027902,
                4214536.311398932,
                -6.412078301235929e-05,
                6.790877998652985e-06,
                0.04426091169980453,
            ],
            65858,
        ],
        "19|Polaris|v1": [
            [
                87.16224833519443,
                66.09306210334928,
                27356108.917828396,
                -0.00020768589108259594,
                -4.4384268579140455e-05,
                -0.0042061979057802295,
            ],
            65858,
        ],
        "19|Polaris|v2": [
            [
                87.16224833519443,
                66.09306210334928,
                27356108.917828396,
                -0.00020768589108259594,
                -4.4384268579140455e-05,
                -0.0042061979057802295,
            ],
            65858,
        ],
        "20|Aldebaran|v1": [
            [
                69.08631504590612,
                -5.473452917217177,
                4214536.311398931,
                -6.412031542178558e-05,
                6.788929747310113e-06,
                0.04426091169980529,
            ],
            65858,
        ],
        "20|Aldebaran|v2": [
            [
                69.08631504590612,
                -5.473452917217177,
                4214536.311398931,
                -6.412031542178558e-05,
                6.788929747310113e-06,
                0.04426091169980529,
            ],
            65858,
        ],
        "20|Polaris|v1": [
            [
                87.86131386109012,
                66.09958560233854,
                27356108.917828392,
                -0.00020774103742281886,
                -4.4382895179701324e-05,
                -0.0042061979057789025,
            ],
            65858,
        ],
        "20|Polaris|v2": [
            [
                87.86131386109012,
                66.09958560233854,
                27356108.917828392,
                -0.00020774103742281886,
                -4.4382895179701324e-05,
                -0.0042061979057789025,
            ],
            65858,
        ],
    },
    "houses": {
        "18|P|1991": [
            [
                284.07013996159634,
                333.2023533527569,
                16.69727592181874,
                45.43435531528889,
                66.27733442691387,
                84.41220495920079,
                104.07013996159635,
                153.20235335275697,
                196.69727592181874,
                225.43435531528888,
                246.27733442691382,
                264.4122049592008,
            ],
            [
                284.07013996159634,
                225.43435531528888,
                222.83621939716525,
                143.38175722789265,
                310.51920408030566,
                325.5985223960679,
                290.04328367315895,
                145.59852239606792,
            ],
        ],
        "18|P|window": [
            [
                110.42381462590168,
                127.61015319580477,
                148.03546955489367,
                175.20015029891997,
                212.47016198875787,
                255.23381785828948,
                290.4238146259017,
                307.6101531958048,
                328.03546955489367,
                355.20015029891994,
                32.47016198875787,
                75.23381785828948,
            ],
            [
                110.42381462590168,
                355.20015029891994,
                355.5959578914171,
                246.73729336004996,
                85.95654567817378,
                62.807791099994894,
                106.09347441696548,
                242.8077910999949,
            ],
        ],
        "18|i|62": [
            [
                260.50645560957304,
                294.43112659134147,
                16.093286154030046,
                45.43435531528889,
                57.62804749588097,
                67.5179199483498,
                80.50645560957305,
                138.5811097611309,
                206.31662774650852,
                225.43435531528888,
                235.32782439743067,
                244.4378457350941,
            ],
            [
                260.50645560957304,
                225.43435531528888,
                222.83621939716525,
                138.85338420676098,
                310.51920408030566,
                332.0084298072968,
                299.47807928942706,
                152.00842980729684,
            ],
        ],
        "18|N|1991": [
            [
                0.0,
                30.0,
                60.0,
                90.0,
                120.0,
                150.0,
                180.0,
                210.0,
                240.0,
                270.0,
                300.0,
                330.0,
            ],
            [
                284.07013996159634,
                225.43435531528888,
                222.83621939716525,
                143.38175722789265,
                310.51920408030566,
                325.5985223960679,
                290.04328367315895,
                145.59852239606792,
            ],
        ],
        "19|P|t0": [
            [
                150.62354748212786,
                194.54403837270817,
                227.28863281763398,
                252.19004362003625,
                274.4306574038926,
                298.47721655889126,
                330.62354748212783,
                14.54403837270819,
                47.28863281763397,
                72.19004362003622,
                94.4306574038926,
                118.47721655889124,
            ],
            [
                150.62354748212786,
                72.19004362003622,
                70.69545161017224,
                347.02446108934345,
                159.10674854334604,
                163.8931739561716,
                131.9841445925957,
                343.89317395617155,
            ],
        ],
        "20|G|2035": [
            [
                105.68633880347359,
                94.12851228126013,
                79.53648728379619,
                60.72114720244815,
                39.43163328054025,
                19.563661360124986,
                2.912948268439087,
                349.43380921937364,
                338.5252866678583,
                329.57173441514425,
                322.0801626965543,
                315.6836811054753,
                310.1138933955037,
                305.1731245186592,
                300.71298732849084,
                296.6189969330337,
                292.7996094215836,
                289.17814492411594,
                285.6863388034736,
                274.12851228126016,
                259.5364872837962,
                240.72114720244815,
                219.43163328054024,
                199.563661360125,
                182.9129482684391,
                169.4338092193737,
                158.52528666785832,
                149.57173441514425,
                142.08016269655425,
                135.68368110547527,
                130.11389339550368,
                125.17312451865922,
                120.71298732849084,
                116.61899693303373,
                112.7996094215836,
                109.17814492411594,
            ],
            [
                105.68633880347359,
                329.57173441514425,
                332.8373509111466,
                232.69977638516076,
                63.58028715336404,
                37.55472711984811,
                76.78794969157775,
                217.55472711984814,
            ],
        ],
    },
}

POS_TOL = 0.15 / 3600.0  # deg; J2000-path floor vs the reference is ~0.05"
SPEED_TOL = 0.01  # deg/day channel tolerance (measured agreement ~1e-5)
HOUSE_TOL = 0.15 / 3600.0


def _dlon(a: float, b: float) -> float:
    return abs((a - b + 180.0) % 360.0 - 180.0)


@pytest.fixture(autouse=True)
def _restore_sid_mode():
    yield
    le.set_sid_mode(0, 0, 0)


@pytest.mark.parametrize("key", sorted(FROZEN["planets"]))
def test_planet_fixed_epoch_position_and_retflag(key):
    sm, jd, body = key.split("|")
    le.set_sid_mode(int(sm), 0, 0)
    expected, exp_rf = FROZEN["planets"][key]
    xx, rf = le.calc_ut(float(jd), int(body), FLG_SWIEPH | FLG_SIDEREAL | FLG_SPEED)
    assert rf == exp_rf
    # The Moon in 1900 carries a ~0.5" delta-T/ephemeris baseline vs the
    # reference (identical on the tropical path) — not a frame error.
    tol = 0.6 / 3600.0 if (body == "1" and jd == "2415020.0") else POS_TOL
    assert _dlon(xx[0], expected[0]) < tol
    assert abs(xx[1] - expected[1]) < POS_TOL
    assert abs(xx[2] - expected[2]) < 1e-6 * max(1.0, expected[2])
    assert abs(xx[3] - expected[3]) < SPEED_TOL
    assert abs(xx[4] - expected[4]) < SPEED_TOL


@pytest.mark.parametrize("key", sorted(FROZEN["stars"]))
def test_star_fixed_epoch_position_and_retflag(key):
    sm, star, ver = key.split("|")
    le.set_sid_mode(int(sm), 0, 0)
    fn = le.fixstar_ut if ver == "v1" else le.fixstar2_ut
    expected, exp_rf = FROZEN["stars"][key]
    res = fn(star, 2448000.5, FLG_SWIEPH | FLG_SIDEREAL | FLG_SPEED)
    xx, rf = res[0], res[-1]
    assert rf == exp_rf
    assert _dlon(xx[0], expected[0]) < POS_TOL
    assert abs(xx[1] - expected[1]) < POS_TOL
    assert abs(xx[3] - expected[3]) < SPEED_TOL
    assert abs(xx[4] - expected[4]) < SPEED_TOL


HOUSE_CASES = {
    "18|P|1991": (18, 2448000.5, 48.0, 16.0, "P"),
    "18|P|window": (18, 2451605.0, 48.0, 16.0, "P"),  # near-t0 sign quirk active
    "18|i|62": (18, 2448000.5, 62.0, 16.0, "i"),  # Sunshine 'i' collapses onto 'I'
    "18|N|1991": (18, 2448000.5, 48.0, 16.0, "N"),  # Aries stays unshifted
    "19|P|t0": (19, 2415020.0, -35.0, 151.0, "P"),  # exactly at the mode's t0
    "20|G|2035": (20, 2465000.5, 60.0, -70.0, "G"),  # Gauquelin 36 sectors
}


@pytest.mark.parametrize("key", sorted(HOUSE_CASES))
def test_houses_fixed_epoch(key):
    sm, jd, lat, lon, hs = HOUSE_CASES[key]
    le.set_sid_mode(sm, 0, 0)
    exp_c, exp_a = FROZEN["houses"][key]
    cusps, ascmc = le.houses_ex(jd, lat, lon, ord(hs), FLG_SIDEREAL)
    n = min(len(cusps), len(exp_c))
    for i in range(n):
        assert _dlon(cusps[i], exp_c[i]) < HOUSE_TOL, f"cusp {i + 1}"
    # ascmc[2] (ARMC) is a documented micro-divergence (<~15" per 5
    # centuries); every other slot must match.
    for i in (0, 1, 3, 4, 5, 6, 7):
        assert _dlon(ascmc[i], exp_a[i]) < HOUSE_TOL, f"ascmc {i}"
    assert _dlon(ascmc[2], exp_a[2]) < 20.0 / 3600.0


def test_mode18_is_j2000_nonut_identity():
    """SIDM_J2000 output is bit-identical to a FLG_J2000|FLG_NONUT request."""
    le.set_sid_mode(18, 0, 0)
    for body in (0, 1, 5, 11):
        sid, _ = le.calc_ut(2448000.5, body, FLG_SWIEPH | FLG_SIDEREAL | FLG_SPEED)
        j2k, _ = le.calc_ut(
            2448000.5, body, FLG_SWIEPH | FLG_J2000 | FLG_NONUT | FLG_SPEED
        )
        assert sid == j2k


def test_mode18_equatorial_matches_j2000_equatorial():
    le.set_sid_mode(18, 0, 0)
    sid, _ = le.calc_ut(
        2448000.5, 1, FLG_SWIEPH | FLG_SIDEREAL | FLG_EQUATORIAL | FLG_SPEED
    )
    j2k, _ = le.calc_ut(
        2448000.5,
        1,
        FLG_SWIEPH | FLG_J2000 | FLG_NONUT | FLG_EQUATORIAL | FLG_SPEED,
    )
    assert sid == j2k


def test_fixed_epoch_star_families_agree():
    """In a fixed frame the legacy and modern star speeds coincide."""
    le.set_sid_mode(19, 0, 0)
    r1 = le.fixstar_ut("Aldebaran", 2448000.5, FLG_SWIEPH | FLG_SIDEREAL | FLG_SPEED)
    r2 = le.fixstar2_ut("Aldebaran", 2448000.5, FLG_SWIEPH | FLG_SIDEREAL | FLG_SPEED)
    for a, b in zip(r1[0], r2[0]):
        assert abs(a - b) < 1e-9


def test_sunshine_pair_collapses_under_fixed_epoch():
    le.set_sid_mode(18, 0, 0)
    ci, _ = le.houses_ex(2448000.5, 62.0, 16.0, ord("i"), FLG_SIDEREAL)
    cI, _ = le.houses_ex(2448000.5, 62.0, 16.0, ord("I"), FLG_SIDEREAL)
    assert max(_dlon(a, b) for a, b in zip(ci, cI)) < 1e-9


def test_context_fixed_epoch_matches_module():
    """The context API mirrors the module-level fixed-epoch handling."""
    from libephemeris import EphemerisContext

    le.set_sid_mode(0, 0, 0)  # module state must not leak into the context
    ctx = EphemerisContext()
    ctx.set_sid_mode(19)
    xx_ctx, rf_ctx = ctx.calc_ut(2448000.5, 1, FLG_SWIEPH | FLG_SIDEREAL | FLG_SPEED)
    le.set_sid_mode(19, 0, 0)
    xx_mod, rf_mod = le.calc_ut(2448000.5, 1, FLG_SWIEPH | FLG_SIDEREAL | FLG_SPEED)
    assert rf_ctx == rf_mod
    for a, b in zip(xx_ctx, xx_mod):
        assert abs(a - b) < 1e-9
