"""Eclipse catalog validation against NASA eclipse data.

Cross-validates libephemeris eclipse predictions against the NASA Eclipse
Website (Espenak Canon) for historical and future eclipses.

Validation Plan v2, Section 5.

NASA Eclipse data sources:
    - https://eclipse.gsfc.nasa.gov/SEcat5/SE2001-2100.html (Five Millennium Canon)
    - https://eclipse.gsfc.nasa.gov/LEcat5/LE2001-2100.html (Lunar Eclipse Canon)

IMPORTANT: NASA catalog times are in TD (Terrestrial Dynamical Time), NOT UT.
Reference UT times below are computed as: UT = TD - DeltaT, using the DeltaT
values published in the NASA catalog for each eclipse. This gives ~0-2 second
agreement with libephemeris eclipse computations.

Note on excluded eclipses:
    - 2015-Apr-04 lunar: Excluded because it's the most borderline total lunar
      eclipse of the century (umbra magnitude 1.0008, total duration 4m43s).
      Different shadow enlargement models (NASA Danjon vs Swiss Ephemeris) may
      classify it differently.
    - 2027-Jul-18 lunar: Excluded because it's the smallest penumbral eclipse
      of the century (penumbral magnitude 0.0014), too marginal to reliably
      detect.
"""

from __future__ import annotations

import warnings

import pytest

import libephemeris as swe

warnings.filterwarnings("ignore")

pytestmark = pytest.mark.slow


def jd_to_datetime_str(jd: float) -> str:
    """Convert JD to human-readable string."""
    y, m, d, h = swe.revjul(jd)
    hh = int(h)
    mm = int((h - hh) * 60)
    ss = int(((h - hh) * 60 - mm) * 60)
    return f"{y}-{m:02d}-{d:02d} {hh:02d}:{mm:02d}:{ss:02d} UT"


def time_diff_seconds(jd1: float, jd2: float) -> float:
    """Difference between two JDs in seconds."""
    return abs(jd1 - jd2) * 86400


# ============================================================================
# §5.1 Solar Eclipses — NASA Espenak Canon (2001-2020)
# ============================================================================

# NASA Five Millennium Canon of Solar Eclipses, 2001-2100 catalog.
# Source: https://eclipse.gsfc.nasa.gov/SEcat5/SE2001-2100.html
#
# UT times computed from NASA TD times: UT = TD - DeltaT
# Each entry: (year, month, day, greatest_eclipse_ut_hours, type, description)
#
# Type codes: T=total, A=annular, H=hybrid (annular-total), P=partial
SOLAR_ECLIPSES_HISTORICAL = [
    # 09511: TD 12:04:46, DeltaT=64s -> UT 12:03:42
    (2001, 6, 21, 12.0617, "T", "Total - S. Atlantic, S. Africa, Madagascar"),
    # 09512: TD 20:53:01, DeltaT=64s -> UT 20:51:57
    (2001, 12, 14, 20.8658, "A", "Annular - Central America, Pacific"),
    # 09513: TD 23:45:22, DeltaT=64s -> UT 23:44:18
    (2002, 6, 10, 23.7383, "A", "Annular - Pacific, Mexico"),
    # 09514: TD 07:32:16, DeltaT=64s -> UT 07:31:12
    (2002, 12, 4, 7.5200, "T", "Total - S. Africa, S. Indian Ocean, Australia"),
    # 09515: TD 04:09:22, DeltaT=64s -> UT 04:08:18
    (2003, 5, 31, 4.1383, "A", "Annular - Iceland, Scotland"),
    # 09516: TD 22:50:22, DeltaT=64s -> UT 22:49:18
    (2003, 11, 23, 22.8217, "T", "Total - Antarctica"),
    # 09519: TD 20:36:51, DeltaT=65s -> UT 20:35:46
    (2005, 4, 8, 20.5961, "H", "Hybrid - Pacific, Panama, Colombia"),
    # 09521: TD 10:12:23, DeltaT=65s -> UT 10:11:18
    (2006, 3, 29, 10.1883, "T", "Total - Africa, Turkey, Russia"),
    # 09526: TD 10:22:12, DeltaT=66s -> UT 10:21:06
    (2008, 8, 1, 10.3517, "T", "Total - Arctic, Siberia, China"),
    # 09528: TD 02:36:25, DeltaT=66s -> UT 02:35:19
    (2009, 7, 22, 2.5886, "T", "Total - India, China, Pacific"),
    # 09529: TD 07:07:39, DeltaT=67s -> UT 07:06:32
    (2010, 1, 15, 7.1089, "A", "Annular - Africa, Indian Ocean, Asia"),
    # 09535: TD 23:53:54, DeltaT=68s -> UT 23:52:46
    (2012, 5, 20, 23.8794, "A", "Annular - China, Japan, Pacific, US"),
    # 09536: TD 22:12:55, DeltaT=68s -> UT 22:11:47
    (2012, 11, 13, 22.1964, "T", "Total - N. Australia, Pacific"),
    # 09541: TD 09:46:47, DeltaT=69s -> UT 09:45:38
    (2015, 3, 20, 9.7606, "T", "Total - N. Atlantic, Svalbard"),
    # 09543: TD 01:58:19, DeltaT=70s -> UT 01:57:09
    (2016, 3, 9, 1.9525, "T", "Total - Indonesia, Pacific"),
    # 09546: TD 18:26:40, DeltaT=70s -> UT 18:25:30
    (2017, 8, 21, 18.4250, "T", "Total - US coast to coast"),
    # 09551: TD 19:24:07, DeltaT=71s -> UT 19:22:56
    (2019, 7, 2, 19.3822, "T", "Total - S. Pacific, Chile, Argentina"),
    # 09552: TD 05:18:53, DeltaT=72s -> UT 05:17:41
    (2019, 12, 26, 5.2947, "A", "Annular - Saudi Arabia, India, SE Asia"),
    # 09553: TD 06:41:15, DeltaT=72s -> UT 06:40:03
    (2020, 6, 21, 6.6675, "A", "Annular - Africa, S. Asia, China"),
    # 09554: TD 16:14:39, DeltaT=72s -> UT 16:13:27
    (2020, 12, 14, 16.2242, "T", "Total - S. Pacific, Chile, Argentina"),
]

# Expected type mapping
TYPE_MAP = {
    "T": swe.ECL_TOTAL,
    "A": swe.ECL_ANNULAR,
    "H": swe.ECL_ANNULAR_TOTAL,
    "P": swe.ECL_PARTIAL,
}


class TestSolarEclipseCatalog:
    """§5.1 Validate solar eclipse predictions against NASA catalog."""

    @pytest.mark.parametrize(
        "year,month,day,expected_hour,ecl_type_str,desc",
        SOLAR_ECLIPSES_HISTORICAL,
        ids=[
            f"{y}-{m:02d}-{d:02d}-{t}" for y, m, d, _, t, _ in SOLAR_ECLIPSES_HISTORICAL
        ],
    )
    def test_solar_eclipse_timing_and_type(
        self,
        year: int,
        month: int,
        day: int,
        expected_hour: float,
        ecl_type_str: str,
        desc: str,
    ) -> None:
        """Verify solar eclipse timing within 60 seconds of NASA catalog."""
        # Search from 10 days before the expected date
        jd_search = swe.julday(year, month, day - 10, 0.0)
        ecl_type, times = swe.sol_eclipse_when_glob(jd_search)

        jd_max = times[0]
        y_found, m_found, d_found, h_found = swe.revjul(jd_max)

        # Verify correct date
        assert y_found == year, f"Wrong year: expected {year}, got {y_found} ({desc})"
        assert m_found == month, (
            f"Wrong month: expected {month}, got {m_found} ({desc})"
        )
        assert d_found == day, f"Wrong day: expected {day}, got {d_found} ({desc})"

        # Verify timing within 60 seconds of NASA UT (derived from TD - DeltaT)
        expected_jd = swe.julday(year, month, day, expected_hour)
        dt_sec = time_diff_seconds(jd_max, expected_jd)
        assert dt_sec < 60, (
            f"Timing error {dt_sec:.1f}s > 60s for {desc}\n"
            f"  Expected: {jd_to_datetime_str(expected_jd)}\n"
            f"  Got:      {jd_to_datetime_str(jd_max)}"
        )

        # Verify eclipse type classification
        expected_flag = TYPE_MAP[ecl_type_str]
        assert ecl_type & expected_flag, (
            f"Wrong type for {desc}: expected {ecl_type_str} "
            f"(flag {expected_flag}), got flags={ecl_type}"
        )

    def test_solar_eclipse_contacts_ordering(self) -> None:
        """Eclipse contact times should be in chronological order."""
        # 2024-04-08 total solar eclipse
        jd = swe.julday(2024, 4, 1, 0.0)
        ecl_type, times = swe.sol_eclipse_when_glob(jd, ecltype=swe.ECL_TOTAL)

        jd_c1 = times[2]  # First contact
        jd_c2 = times[4]  # Second contact (totality begins)
        jd_max = times[0]  # Maximum
        jd_c3 = times[5]  # Third contact (totality ends)
        jd_c4 = times[3]  # Fourth contact

        assert jd_c1 > 0, "C1 should be set for total eclipse"
        assert jd_c4 > 0, "C4 should be set for total eclipse"
        assert jd_c2 > 0, "C2 should be set for total eclipse"
        assert jd_c3 > 0, "C3 should be set for total eclipse"

        assert jd_c1 < jd_c2 < jd_max < jd_c3 < jd_c4, (
            f"Contact times not in order: C1={jd_c1} C2={jd_c2} "
            f"max={jd_max} C3={jd_c3} C4={jd_c4}"
        )


# ============================================================================
# §5.2 Lunar Eclipses — NASA catalog (2001-2022)
# ============================================================================

# NASA Five Millennium Catalog of Lunar Eclipses, 2001-2100 catalog.
# Source: https://eclipse.gsfc.nasa.gov/LEcat5/LE2001-2100.html
#
# UT times computed from NASA TD times: UT = TD - DeltaT
# Each entry: (year, month, day, greatest_eclipse_ut_hours, type, description)
#
# Note: 2015-Apr-04 excluded (borderline total, umbra mag 1.0008, 4m43s total)
# Note: Entries span 2001-2022 to cover 20 eclipses from the catalog
LUNAR_ECLIPSES_HISTORICAL = [
    # 09651: TD 20:21:40, DeltaT=64s -> UT 20:20:36
    (2001, 1, 9, 20.3433, "T", "Total lunar eclipse"),
    # 09657: TD 03:41:13, DeltaT=64s -> UT 03:40:09
    (2003, 5, 16, 3.6692, "T", "Total lunar eclipse"),
    # 09658: TD 01:19:38, DeltaT=64s -> UT 01:18:34
    (2003, 11, 9, 1.3094, "T", "Total lunar eclipse"),
    # 09659: TD 20:31:17, DeltaT=65s -> UT 20:30:12
    (2004, 5, 4, 20.5033, "T", "Total lunar eclipse"),
    # 09660: TD 03:05:11, DeltaT=65s -> UT 03:04:06
    (2004, 10, 28, 3.0683, "T", "Total lunar eclipse"),
    # 09665: TD 23:21:59, DeltaT=65s -> UT 23:20:54
    (2007, 3, 3, 23.3483, "T", "Total lunar eclipse"),
    # 09666: TD 10:38:27, DeltaT=66s -> UT 10:37:21
    (2007, 8, 28, 10.6225, "T", "Total lunar eclipse"),
    # 09667: TD 03:27:09, DeltaT=66s -> UT 03:26:03
    (2008, 2, 21, 3.4342, "T", "Total lunar eclipse"),
    # 09673: TD 11:39:34, DeltaT=67s -> UT 11:38:27
    (2010, 6, 26, 11.6408, "P", "Partial lunar eclipse"),
    # 09674: TD 08:18:04, DeltaT=67s -> UT 08:16:57
    (2010, 12, 21, 8.2825, "T", "Total lunar eclipse"),
    # 09675: TD 20:13:43, DeltaT=67s -> UT 20:12:36
    (2011, 6, 15, 20.2100, "T", "Total lunar eclipse"),
    # 09676: TD 14:32:56, DeltaT=68s -> UT 14:31:48
    (2011, 12, 10, 14.5300, "T", "Total lunar eclipse"),
    # 09682: TD 07:46:48, DeltaT=69s -> UT 07:45:39
    (2014, 4, 15, 7.7608, "T", "Total lunar eclipse"),
    # 09685: TD 02:48:17, DeltaT=69s -> UT 02:47:08
    (2015, 9, 28, 2.7856, "T", "Total lunar eclipse"),
    # 09690: TD 13:31:00, DeltaT=71s -> UT 13:29:49
    (2018, 1, 31, 13.4969, "T", "Total lunar eclipse"),
    # 09691: TD 20:22:54, DeltaT=71s -> UT 20:21:43
    (2018, 7, 27, 20.3619, "T", "Total lunar eclipse"),
    # 09692: TD 05:13:27, DeltaT=71s -> UT 05:12:16
    (2019, 1, 21, 5.2044, "T", "Total lunar eclipse"),
    # 09698: TD 11:19:53, DeltaT=72s -> UT 11:18:41
    (2021, 5, 26, 11.3114, "T", "Total lunar eclipse"),
    # 09700: TD 04:12:42, DeltaT=73s -> UT 04:11:29
    (2022, 5, 16, 4.1914, "T", "Total lunar eclipse"),
    # 09701: TD 11:00:22, DeltaT=73s -> UT 10:59:09
    (2022, 11, 8, 10.9858, "T", "Total lunar eclipse"),
]

LUNAR_TYPE_MAP = {
    "T": swe.ECL_TOTAL,
    "P": swe.ECL_PARTIAL,
    "Pen": swe.ECL_PENUMBRAL,
}


class TestLunarEclipseCatalog:
    """§5.2 Validate lunar eclipse predictions against NASA catalog."""

    @pytest.mark.parametrize(
        "year,month,day,expected_hour,ecl_type_str,desc",
        LUNAR_ECLIPSES_HISTORICAL,
        ids=[
            f"{y}-{m:02d}-{d:02d}-{t}" for y, m, d, _, t, _ in LUNAR_ECLIPSES_HISTORICAL
        ],
    )
    def test_lunar_eclipse_timing_and_type(
        self,
        year: int,
        month: int,
        day: int,
        expected_hour: float,
        ecl_type_str: str,
        desc: str,
    ) -> None:
        """Verify lunar eclipse timing within 60 seconds of NASA catalog."""
        jd_search = swe.julday(year, month, day - 10, 0.0)
        expected_flag = LUNAR_TYPE_MAP[ecl_type_str]
        ecl_type, times = swe.lun_eclipse_when(jd_search, ecltype=expected_flag)

        jd_max = times[0]
        y_found, m_found, d_found, h_found = swe.revjul(jd_max)

        # Verify correct date
        assert y_found == year, f"Wrong year: expected {year}, got {y_found} ({desc})"
        assert m_found == month, (
            f"Wrong month: expected {month}, got {m_found} ({desc})"
        )
        assert d_found == day, f"Wrong day: expected {day}, got {d_found} ({desc})"

        # Verify timing within 60 seconds of NASA UT (derived from TD - DeltaT)
        expected_jd = swe.julday(year, month, day, expected_hour)
        dt_sec = time_diff_seconds(jd_max, expected_jd)
        assert dt_sec < 60, (
            f"Timing error {dt_sec:.1f}s > 60s for {desc}\n"
            f"  Expected: {jd_to_datetime_str(expected_jd)}\n"
            f"  Got:      {jd_to_datetime_str(jd_max)}"
        )

        # Verify eclipse type
        assert ecl_type & expected_flag, (
            f"Wrong type for {desc}: expected {ecl_type_str}, got flags={ecl_type}"
        )

    def test_lunar_eclipse_contacts_ordering(self) -> None:
        """Lunar eclipse contact times in chronological order."""
        # 2022-05-16 total lunar eclipse
        jd = swe.julday(2022, 5, 10, 0.0)
        ecl_type, times = swe.lun_eclipse_when(jd, ecltype=swe.ECL_TOTAL)

        jd_max = times[0]
        jd_p1 = times[6]  # Penumbral begin
        jd_u1 = times[2]  # Partial begin (umbral)
        jd_u2 = times[4]  # Total begin
        jd_u3 = times[5]  # Total end
        jd_u4 = times[3]  # Partial end (umbral)
        jd_p4 = times[7]  # Penumbral end

        assert jd_p1 > 0, "P1 should be set for total lunar eclipse"
        assert jd_u1 > 0, "U1 should be set for total lunar eclipse"
        assert jd_u2 > 0, "U2 should be set for total lunar eclipse"
        assert jd_u3 > 0, "U3 should be set for total lunar eclipse"
        assert jd_u4 > 0, "U4 should be set for total lunar eclipse"
        assert jd_p4 > 0, "P4 should be set for total lunar eclipse"

        assert jd_p1 < jd_u1 < jd_u2 < jd_max, (
            f"Pre-max contacts not in order: P1={jd_p1} U1={jd_u1} "
            f"U2={jd_u2} max={jd_max}"
        )
        assert jd_max < jd_u3 < jd_u4 < jd_p4, (
            f"Post-max contacts not in order: max={jd_max} U3={jd_u3} "
            f"U4={jd_u4} P4={jd_p4}"
        )

    def test_lunar_eclipse_gamma_sign(self) -> None:
        """Verify gamma sign for a known lunar eclipse."""
        # 2022-11-08 total lunar eclipse
        jd = swe.julday(2022, 11, 1, 0.0)
        ecl_type, times = swe.lun_eclipse_when(jd, ecltype=swe.ECL_TOTAL)

        gamma = swe.lun_eclipse_gamma(times[0])
        # Gamma should be a finite number for a valid eclipse
        assert -2.0 < gamma < 2.0, f"Gamma out of range: {gamma}"


# ============================================================================
# §5.3 Future Eclipses (2025-2035)
# ============================================================================

# Future solar eclipses from NASA Five Millennium Canon
# Source: https://eclipse.gsfc.nasa.gov/SEcat5/SE2001-2100.html
FUTURE_SOLAR = [
    # (year, month, day, type)
    (2026, 2, 17, "A"),  # 09565: Annular - Antarctica
    (2026, 8, 12, "T"),  # 09566: Total - Arctic, Spain, Iceland
    (2027, 2, 6, "A"),  # 09567: Annular - S. America, Africa
    (2027, 8, 2, "T"),  # 09568: Total - Spain, N. Africa, Middle East
    (2028, 1, 26, "A"),  # 09569: Annular - S. America
    (2028, 7, 22, "T"),  # 09570: Total - Australia, New Zealand
    (2029, 1, 14, "P"),  # 09571: Partial
    (2029, 6, 12, "P"),  # 09572: Partial
    (2029, 7, 11, "P"),  # 09573: Partial
    (2029, 12, 5, "P"),  # 09574: Partial
]

# Future lunar eclipses from NASA Five Millennium Catalog
# Source: https://eclipse.gsfc.nasa.gov/LEcat5/LE2001-2100.html
#
# Note: 2027-Jul-18 excluded (penumbral magnitude 0.0014, smallest of the
# century, too marginal to reliably detect across ephemeris implementations)
FUTURE_LUNAR = [
    # (year, month, day, type)
    (2025, 3, 14, "T"),  # 09706: Total
    (2025, 9, 7, "T"),  # 09707: Total
    (2026, 3, 3, "T"),  # 09708: Total
    (2026, 8, 28, "P"),  # 09709: Partial
    (2027, 2, 20, "Pen"),  # 09710: Penumbral
    (2027, 8, 17, "Pen"),  # 09712: Penumbral
    (2028, 1, 12, "P"),  # 09713: Partial
    (2028, 7, 6, "P"),  # 09714: Partial
    (2028, 12, 31, "T"),  # 09715: Total
    (2029, 12, 20, "T"),  # 09717: Total
]


class TestFutureEclipses:
    """§5.3 Verify future eclipse predictions match NASA forecasts."""

    @pytest.mark.parametrize(
        "year,month,day,ecl_type_str",
        FUTURE_SOLAR,
        ids=[f"solar-{y}-{m:02d}-{d:02d}-{t}" for y, m, d, t in FUTURE_SOLAR],
    )
    def test_future_solar_eclipse_date_and_type(
        self, year: int, month: int, day: int, ecl_type_str: str
    ) -> None:
        """Future solar eclipse date and type match NASA prediction."""
        jd_search = swe.julday(year, month, day - 15, 0.0)

        if ecl_type_str == "P":
            # Partial eclipses: search for any type
            ecl_type, times = swe.sol_eclipse_when_glob(jd_search)
        else:
            expected_flag = TYPE_MAP[ecl_type_str]
            ecl_type, times = swe.sol_eclipse_when_glob(
                jd_search, ecltype=expected_flag
            )

        y_found, m_found, d_found, _ = swe.revjul(times[0])

        assert y_found == year and m_found == month and d_found == day, (
            f"Wrong date: expected {year}-{month:02d}-{day:02d}, "
            f"got {y_found}-{m_found:02d}-{d_found:02d}"
        )

        expected_flag = TYPE_MAP[ecl_type_str]
        assert ecl_type & expected_flag, (
            f"Wrong type: expected {ecl_type_str}, got flags={ecl_type}"
        )

    @pytest.mark.parametrize(
        "year,month,day,ecl_type_str",
        FUTURE_LUNAR,
        ids=[f"lunar-{y}-{m:02d}-{d:02d}-{t}" for y, m, d, t in FUTURE_LUNAR],
    )
    def test_future_lunar_eclipse_date_and_type(
        self, year: int, month: int, day: int, ecl_type_str: str
    ) -> None:
        """Future lunar eclipse date and type match NASA prediction."""
        jd_search = swe.julday(year, month, day - 15, 0.0)
        expected_flag = LUNAR_TYPE_MAP[ecl_type_str]

        ecl_type, times = swe.lun_eclipse_when(jd_search, ecltype=expected_flag)

        y_found, m_found, d_found, _ = swe.revjul(times[0])

        assert y_found == year and m_found == month and d_found == day, (
            f"Wrong date: expected {year}-{month:02d}-{day:02d}, "
            f"got {y_found}-{m_found:02d}-{d_found:02d}"
        )

        assert ecl_type & expected_flag, (
            f"Wrong type: expected {ecl_type_str}, got flags={ecl_type}"
        )
