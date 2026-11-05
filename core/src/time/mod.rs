//! Calendar, timescale and Earth-rotation kernels.
//!
//! [`julday`] and [`revjul`] use Meeus calendar arithmetic with astronomical
//! year numbering (year 0 is 1 BCE) and explicit proleptic calendars. Julian
//! days begin at noon. [`date_conversion`] instead selects its input calendar
//! at the 1582 reform, and the MJD calendar helpers have their own documented
//! conventions. A calendar label alone does not identify TT, UT1, TDB or UTC.
//!
//! The module provides UTC/TAI/TT/UT1 conversion, the exact 32.184-second TAI–TT
//! offset, observed and model Delta-T, periodic TT–TDB terms, mean/apparent
//! sidereal time, and equation-of-time/local-solar-time calculations. Delta-T
//! configuration and providers are explicit; unavailable selected models return
//! typed errors. Earth-orientation and leap-second inputs remain caller-owned.
//!
//! Time offsets are seconds unless an entry explicitly returns days; in
//! particular, configured Delta-T is returned in days. Dynamical entries use TT
//! and Earth rotation uses UT1. Angular outputs are degrees or named hours.
//! [`convert_date_checked`] reports impossible dates through `valid = false`,
//! not an error; UTC conversion separately rejects invalid civil labels.
//!
//! UTC calendar/JD helpers use ordinary 86400-second calendar-day arithmetic,
//! not a leap-second-preserving quasi-JD representation. Their civil-era and
//! leap-label limitations are documented at [`utc_to_jd_with_config`]. Two-part
//! epochs retain small offsets, but several conversions deliberately use a
//! collapsed scalar for lookup or calendar reconstruction.
//!
//! Checked entries describe their nonfinite/overflow errors individually;
//! unchecked formula helpers can propagate IEEE values. Borrowed callbacks run
//! synchronously and may perform I/O or panic; no general callback-purity or
//! panic-catching guarantee is made. Extrapolated models are estimates, not
//! guarantees of ephemeris accuracy.

use crate::constants::{FLG_MOSEPH, GREG_CAL, JUL_CAL, TIDAL_DEFAULT, TIDAL_MOSEPH};
use crate::errors::{Error, Result};
use crate::types::TwoPartEpoch;

mod iau2006_sidereal;
pub use iau2006_sidereal::sidereal_time_iau2006_hours;

#[deprecated(note = "Legacy boundary-calibrated convention; see issue #50 and MIGRATION.md")]
pub mod legacy_sidereal;

// Numerically qualified continuous VCW origin; model/clock roles stay separate.
mod canonical_sidereal;
mod canonical_sidereal_table;

/// Explicit proleptic calendar, without automatic historical reform switching.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Calendar {
    /// Julian calendar (`JUL_CAL == 0`).
    Julian,
    /// Gregorian calendar, proleptic for all dates (`GREG_CAL == 1`).
    Gregorian,
}

impl Calendar {
    /// Convert a calendar flag to [`Calendar`].
    ///
    /// Only [`JUL_CAL`] and [`GREG_CAL`] are accepted. Other flags return
    /// [`Error::InvalidCalendar`].
    pub fn from_flag(flag: i32) -> Result<Self> {
        match flag {
            f if f == JUL_CAL => Ok(Self::Julian),
            f if f == GREG_CAL => Ok(Self::Gregorian),
            _ => Err(Error::InvalidCalendar { value: flag }),
        }
    }

    /// Integer calendar flag for this system (`JUL_CAL` or `GREG_CAL`).
    pub fn flag(self) -> i32 {
        match self {
            Self::Julian => JUL_CAL,
            Self::Gregorian => GREG_CAL,
        }
    }
}

/// Calendar date with a fractional-hour component, as returned by
/// [`revjul`]: `(year, month, day, hour)` with `hour` in `[0, 24)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CalendarDate {
    /// Astronomical year numbering (year 0 == 1 BCE, negative == BCE).
    pub year: i32,
    /// Month `1..=12` for well-formed inputs.
    pub month: i32,
    /// Integer day of month for well-formed inputs.
    pub day: i32,
    /// Decimal hour `0.0..24.0`.
    pub hour: f64,
}

/// Convert a calendar date to a Julian Day number.
///
/// `year` uses astronomical numbering, and `hour` is decimal hours, not an
/// `hh.mm` encoding. The chosen calendar applies at every date; no reform gap
/// is inserted. The result has the same time scale as the supplied label.
///
/// Month/day/hour bounds are not validated: finite values flow through the
/// recurrence. Use [`convert_date_checked`] and inspect `valid` to check date
/// labels. A nonfinite hour returns [`Error::NonFiniteArgument`].
pub fn julday(year: i32, month: i32, day: i32, hour: f64, cal: Calendar) -> Result<f64> {
    if !hour.is_finite() {
        return Err(Error::NonFiniteArgument { value: hour });
    }
    // i32 -> f64 is exact (|v| < 2^53), so the float recurrence starts
    // from the integer inputs (narrowed
    // to binary64 by the first float multiply there).
    let mut y = year as f64;
    let mut m = month as f64;
    if m <= 2.0 {
        y -= 1.0;
        m += 12.0;
    }
    // Meeus INT is floor, not truncation: `f64::floor` here pins the
    // negative-year behavior (e.g. year -101 takes a = -2, not -1).
    let a = (y / 100.0).floor();
    let b = match cal {
        Calendar::Gregorian => 2.0 - a + (a / 4.0).floor(),
        Calendar::Julian => 0.0,
    };
    Ok((365.25 * (y + 4716.0)).floor()
        + (30.6001 * (m + 1.0)).floor()
        + day as f64
        + hour / 24.0
        + b
        - 1524.5)
}

/// Convert a Julian Day number to a calendar date.
///
/// The selected proleptic calendar applies even before 1582. Returned years
/// use astronomical numbering and hours are decimal hours. This relabels the
/// supplied epoch; it does not convert time scales or reconstruct leap seconds.
/// Scalar JD rounding limits the reconstructed time-of-day resolution.
/// Returns [`Error::NonFiniteArgument`] for nonfinite JD and
/// [`Error::IntegerOverflow`] when derived components do not fit `i32`.
pub fn revjul(jd: f64, cal: Calendar) -> Result<CalendarDate> {
    if !jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd });
    }
    let jd_shifted = jd + 0.5;
    let z = jd_shifted.floor();
    let f = jd_shifted - z;

    // Compatibility contract: always respect `cal`; no JD-threshold
    // auto-detection of the Gregorian reform.
    let a = match cal {
        Calendar::Gregorian => {
            let alpha = ((z - 1_867_216.25) / 36_524.25).floor();
            z + 1.0 + alpha - (alpha / 4.0).floor()
        }
        Calendar::Julian => z,
    };

    let b = a + 1524.0;
    let c = ((b - 122.1) / 365.25).floor();
    let d = (365.25 * c).floor();
    let e = ((b - d) / 30.6001).floor();

    let day_float = b - d - (30.6001 * e).floor() + f;
    let month_float = if e < 14.0 { e - 1.0 } else { e - 13.0 };
    let year_float = if month_float > 2.0 {
        c - 4716.0
    } else {
        c - 4715.0
    };

    // The integer day is obtained by truncating toward zero for the
    // non-negative day values this recurrence produces.
    let day_trunc = day_float.trunc();
    let hour = (day_float - day_trunc) * 24.0;
    let year = float_to_i32(year_float)?;
    let month = float_to_i32(month_float)?;
    let day = float_to_i32(day_trunc)?;
    Ok(CalendarDate {
        year,
        month,
        day,
        hour,
    })
}

/// Narrow an integral-valued `f64` to `i32`, erroring outside the range.
fn float_to_i32(value: f64) -> Result<i32> {
    if !value.is_finite() || value < i32::MIN as f64 || value > i32::MAX as f64 {
        return Err(Error::IntegerOverflow);
    }
    // `as` saturates, but the range check above already excludes the
    // saturation domain, so this conversion is exact.
    Ok(value as i32)
}

/// Julian Day of the Gregorian calendar reform: noon on 1582-10-15
/// (Gregorian), i.e. `julday(1582, 10, 15, 12.0, Gregorian)`.
pub const JD_GREGORIAN_REFORM: f64 = 2_299_161.0;

/// Gregorian dates on which a positive leap second was inserted at 23:59:60
/// UTC (end of June 30 or December 31).
///
/// IERS Bulletin C / BIPM lineage, sorted ascending through 2016-12-31.
/// This list controls civil leap-label admission independently of any supplied
/// [`LeapSecondEntry`] offset table; newer offset rows do not extend the list.
pub const LEAP_SECOND_DATES: [(i32, i32, i32); 27] = [
    (1972, 6, 30),
    (1972, 12, 31),
    (1973, 12, 31),
    (1974, 12, 31),
    (1975, 12, 31),
    (1976, 12, 31),
    (1977, 12, 31),
    (1978, 12, 31),
    (1979, 12, 31),
    (1981, 6, 30),
    (1982, 6, 30),
    (1983, 6, 30),
    (1985, 6, 30),
    (1987, 12, 31),
    (1989, 12, 31),
    (1990, 12, 31),
    (1992, 6, 30),
    (1993, 6, 30),
    (1994, 6, 30),
    (1995, 12, 31),
    (1997, 6, 30),
    (1998, 12, 31),
    (2005, 12, 31),
    (2008, 12, 31),
    (2012, 6, 30),
    (2015, 6, 30),
    (2016, 12, 31),
];

/// Return `true` when a positive leap second was inserted at the end of this
/// Gregorian date.
pub fn is_leap_second_date(year: i32, month: i32, day: i32) -> bool {
    LEAP_SECOND_DATES.contains(&(year, month, day))
}

/// Return the number of days in a calendar month (28--31).
///
/// February follows the selected calendar's leap-year rule with astronomical
/// year numbering. Any month other than February and the four 30-day months
/// returns 31, including out-of-range month numbers. This is not validation.
/// Use [`convert_date_checked`] and inspect its `valid` field for that purpose.
/// Divisibility applies equally to positive and negative years.
pub fn days_in_month(year: i32, month: i32, cal: Calendar) -> i32 {
    if month == 2 {
        let leap = if cal == Calendar::Gregorian {
            year % 4 == 0 && (year % 100 != 0 || year % 400 == 0)
        } else {
            year % 4 == 0
        };
        if leap { 29 } else { 28 }
    } else if month == 4 || month == 6 || month == 9 || month == 11 {
        30
    } else {
        31
    }
}

/// Target calendar for [`date_conversion`], using the same physical model `'j'` /
/// `'g'` calendar codes (case-insensitive).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CalendarTarget {
    /// Julian calendar (`'j'`).
    Julian,
    /// Gregorian calendar (`'g'`).
    Gregorian,
}

impl CalendarTarget {
    /// Parse a calendar code (`'j'`/`'J'` or `'g'`/`'G'`).
    ///
    /// Longer-string and non-string normalization (bytes decoding, type
    /// rejection) belongs to the language boundary; the core engine
    /// takes a single character, mirroring the validated code point.
    pub fn from_code(code: char) -> Result<Self> {
        match code.to_ascii_lowercase() {
            'j' => Ok(Self::Julian),
            'g' => Ok(Self::Gregorian),
            _ => Err(Error::InvalidCalendarCode { code }),
        }
    }

    /// the model calendar code for this target.
    pub fn code(self) -> char {
        match self {
            Self::Julian => 'j',
            Self::Gregorian => 'g',
        }
    }

    /// Corresponding [`Calendar`] flag value.
    pub fn calendar(self) -> Calendar {
        match self {
            Self::Julian => Calendar::Julian,
            Self::Gregorian => Calendar::Gregorian,
        }
    }
}

/// Convert a date between the Julian and Gregorian calendars.
///
/// Input dates before 1582-10-15 are interpreted as Julian, and dates on/after
/// it as Gregorian, using the date fields rather than an hour-inclusive JD.
/// Unlike [`julday`], this entry does not let the caller select the input
/// calendar. It does not reject historical reform-gap labels or invalid dates.
/// When input and target calendars coincide, finite-hour fields are returned
/// unchanged, without normalization. Years are astronomical; hours are decimal.
///
/// Worked reform examples:
/// `(1582, 10, 15, 12.0)` to Julian returns `(1582, 10, 5, 12.0)`, and
/// `(1582, 10, 5, 12.0)` to Gregorian returns `(1582, 10, 15, 12.0)`.
///
/// ```rust
/// use libephemeris_core::time::{CalendarDate, CalendarTarget, date_conversion};
///
/// let to_julian =
///     date_conversion(1582, 10, 15, 12.0, CalendarTarget::Julian).unwrap();
/// assert_eq!(
///     to_julian,
///     CalendarDate { year: 1582, month: 10, day: 5, hour: 12.0 }
/// );
/// let to_gregorian =
///     date_conversion(1582, 10, 5, 12.0, CalendarTarget::Gregorian).unwrap();
/// assert_eq!(
///     to_gregorian,
///     CalendarDate { year: 1582, month: 10, day: 15, hour: 12.0 }
/// );
/// ```
///
/// Non-finite `hour` is [`Error::NonFiniteArgument`].
pub fn date_conversion(
    year: i32,
    month: i32,
    day: i32,
    hour: f64,
    target: CalendarTarget,
) -> Result<CalendarDate> {
    if !hour.is_finite() {
        return Err(Error::NonFiniteArgument { value: hour });
    }
    let input = if (year, month, day) < (1582, 10, 15) {
        Calendar::Julian
    } else {
        Calendar::Gregorian
    };
    let want = target.calendar();
    if input == want {
        return Ok(CalendarDate {
            year,
            month,
            day,
            hour,
        });
    }
    revjul(julday(year, month, day, hour, input)?, want)
}

/// Date-validity flag, computed Julian day and reconstructed calendar fields.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CheckedDate {
    /// `true` when the input labels a real date in `cal`.
    pub valid: bool,
    /// Julian Day number of the (possibly normalized) input.
    pub jd: f64,
    /// Reconstructed fields. Date fields match the input when `valid`; the
    /// fractional hour can differ through scalar JD rounding.
    pub date: CalendarDate,
}

/// Convert and validate a calendar date, returning its Julian Day number.
///
/// `jd` comes from [`julday`]; `valid` compares reconstructed year/month/day
/// with the inputs, excluding the hour. An impossible date is therefore a
/// successful result with `valid = false`, not [`Error::InvalidDate`]. Hours
/// are decimal; their conventional `[0, 24)` range is not separately checked.
/// Nonfinite hours and unrepresentable reconstructed integer fields propagate
/// the errors of [`julday`] and [`revjul`]. No time-scale conversion occurs.
///
/// ```
/// use libephemeris_core::time::{Calendar, convert_date_checked};
///
/// let date = convert_date_checked(2024, 2, 30, 12.0, Calendar::Gregorian)?;
/// assert!(!date.valid); // Ok does not mean the original date existed
/// assert_eq!((date.date.month, date.date.day), (3, 1));
/// # Ok::<(), libephemeris_core::errors::Error>(())
/// ```
pub fn convert_date_checked(
    year: i32,
    month: i32,
    day: i32,
    hour: f64,
    cal: Calendar,
) -> Result<CheckedDate> {
    let jd = julday(year, month, day, hour, cal)?;
    let back = revjul(jd, cal)?;
    Ok(CheckedDate {
        valid: back.year == year && back.month == month && back.day == day,
        jd,
        date: back,
    })
}

/// Convert and validate a civil date selected by a calendar string.
///
/// `cal` must contain exactly one ASCII calendar code: `g`/`G` selects the
/// proleptic Gregorian calendar and `j`/`J` the proleptic Julian calendar.
/// Year uses astronomical numbering and hour is decimal hours. Invalid dates
/// are normalized by the calendar recurrence and reported with `valid=false`;
/// an invalid code or nonfinite hour is an error. This is a stateless operation
/// and reads no ephemeris or time-scale data.
pub fn convert_date_code(
    year: i32,
    month: i32,
    day: i32,
    hour: f64,
    cal: &str,
) -> Result<CheckedDate> {
    let mut codes = cal.chars();
    let code = codes.next().unwrap_or('\0');
    if codes.next().is_some() {
        return Err(Error::InvalidCalendarCode { code });
    }
    convert_date_checked(
        year,
        month,
        day,
        hour,
        CalendarTarget::from_code(code)?.calendar(),
    )
}

/// Day of the week for a Julian Day number: 0 = Monday through 6 = Sunday.
///
/// (`floor(jd + 0.5) % 7`, the reference-API
/// indexing). [`rem_euclid`](f64::rem_euclid) supplies the
/// non-negative `%` for pre-epoch (negative) Julian Days. J2000.0
/// (`2451545.0`, a Saturday) returns `5`. Non-finite `jd`
/// is [`Error::NonFiniteArgument`].
pub fn day_of_week(jd: f64) -> Result<i32> {
    if !jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd });
    }
    Ok((jd + 0.5).floor().rem_euclid(7.0) as i32)
}

/// Civil date/time fields returned by UTC and timezone conversion entries.
///
/// The producing function determines calendar and era conventions. Public
/// fields do not validate dates or prove that a leap-second label exists.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UtcFields {
    /// Calendar year.
    pub year: i32,
    /// Month `1..=12`.
    pub month: i32,
    /// Day of month.
    pub day: i32,
    /// Hour `0..=23`.
    pub hour: i32,
    /// Minute `0..=59`.
    pub minute: i32,
    /// Seconds within the minute, including a `60.x` label where the producing
    /// function admits one. Reconstruction and fractional timezone offsets
    /// can round or change this value.
    pub second: f64,
}

/// Convert local time to UTC by subtracting the timezone offset.
///
/// `timezone_offset` is hours east of UTC and is subtracted; this is a fixed
/// offset conversion, not a timezone-database or daylight-saving lookup.
/// Offset minutes round with ties to even, and any sub-minute residue changes
/// seconds. Day crossings use proleptic Gregorian noon-JD arithmetic.
///
/// Input civil bounds and actual leap dates are not checked. The normalization
/// rule preserves a `60.x` label where possible; it does not certify that the
/// label exists. Nonfinite seconds/offsets return [`Error::NonFiniteArgument`];
/// unrepresentable minute arithmetic returns [`Error::IntegerOverflow`].
///
/// Worked examples:
/// `2024-01-15 12:00:00` at offset `+1` (CET) returns `11:00:00 UTC`;
/// at offset `-5` (EST) it returns `17:00:00 UTC`.
///
/// ```rust
/// use libephemeris_core::time::utc_time_zone;
///
/// let cet = utc_time_zone(2024, 1, 15, 12, 0, 0.0, 1.0).unwrap();
/// assert_eq!((cet.year, cet.month, cet.day, cet.hour, cet.minute, cet.second),
///     (2024, 1, 15, 11, 0, 0.0));
/// let est = utc_time_zone(2024, 1, 15, 12, 0, 0.0, -5.0).unwrap();
/// assert_eq!((est.year, est.month, est.day, est.hour, est.minute, est.second),
///     (2024, 1, 15, 17, 0, 0.0));
/// ```
#[allow(clippy::too_many_arguments)]
pub fn utc_time_zone(
    year: i32,
    month: i32,
    day: i32,
    hour: i32,
    minute: i32,
    second: f64,
    timezone_offset: f64,
) -> Result<UtcFields> {
    if !timezone_offset.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: timezone_offset,
        });
    }
    if !second.is_finite() {
        return Err(Error::NonFiniteArgument { value: second });
    }
    // Component-wise arithmetic (no decimal-hour round trip): the seconds
    // value passes through bit-exact.
    let offset_min_f = timezone_offset * 60.0;
    if !offset_min_f.is_finite() {
        return Err(Error::IntegerOverflow);
    }
    if offset_min_f < i64::MIN as f64 || offset_min_f > i64::MAX as f64 {
        return Err(Error::IntegerOverflow);
    }
    // Use ties-to-even rounding for half-minute boundaries.
    let offset_min = offset_min_f.round_ties_even() as i64;
    // Sub-minute offset residue (rare, e.g. historical LMT zones).
    let offset_sec = (offset_min_f - offset_min as f64) * 60.0;

    let mut out_second = second - offset_sec;
    let mut minute_carry: i64 = 0;
    // Keep a leap second (60 <= s < 61) untouched by normalization.
    if out_second >= 61.0 || (out_second >= 60.0 && second < 60.0) {
        out_second -= 60.0;
        minute_carry += 1;
    } else if out_second < 0.0 {
        out_second += 60.0;
        minute_carry -= 1;
    }

    // Euclidean division keeps the remainder nonnegative.
    let total_min = (hour as i64)
        .checked_mul(60)
        .and_then(|h| h.checked_add(minute as i64))
        .and_then(|t| t.checked_sub(offset_min))
        .and_then(|t| t.checked_add(minute_carry))
        .ok_or(Error::IntegerOverflow)?;
    let day_shift = total_min.div_euclid(1440);
    let tod_min = total_min.rem_euclid(1440);

    let out_hour = (tod_min / 60) as i32;
    let out_minute = (tod_min % 60) as i32;

    let (out_year, out_month, out_day) = if day_shift == 0 {
        (year, month, day)
    } else {
        // Date arithmetic via JD at noon (immune to time-of-day rounding).
        // `day_shift` always fits in 2^53 (`|total_min| <= i64::MAX`), so
        // the `as f64` widening is exact.
        let jd_noon = julday(year, month, day, 12.0, Calendar::Gregorian)? + day_shift as f64;
        let shifted = revjul(jd_noon, Calendar::Gregorian)?;
        (shifted.year, shifted.month, shifted.day)
    };

    Ok(UtcFields {
        year: out_year,
        month: out_month,
        day: out_day,
        hour: out_hour,
        minute: out_minute,
        second: out_second,
    })
}

// Leap-second tables and UTC/TAI transitions.
//
// Parsers accept text contents and skip nonfinite rows. Callers supply the
// leap-second table explicitly; filesystem and refresh policy live in the
// IERS cache module. The packaged hardcoded table provides deterministic
// historical values. Nonfinite Julian days return NonFiniteArgument;
// calendar years outside i32 return IntegerOverflow.

/// TT - TAI offset in seconds: `TT = TAI + 32.184 s`, exactly.
pub const TT_TAI_OFFSET_SECONDS: f64 = 32.184;

/// Binary64 day representation of the defining 32.184-SI-second TT−TAI offset.
/// The seconds-to-days division rounds; the physical definition is exact.
pub const TT_TAI_OFFSET_DAYS: f64 = TT_TAI_OFFSET_SECONDS / 86400.0;

/// MJD of 1972-01-01, the first integer-leap-second epoch. Before this the
/// UTC system was rate-based and the integer table does not apply.
pub const FIRST_LEAP_SECOND_MJD: f64 = 41317.0;

/// MJD of 2017-01-01, used by the empty-table fallback in the source.
const EMPTY_TABLE_MODERN_MJD: f64 = 57754.0;

/// One leap-second table row: the MJD from which `tai_utc` applies, the
/// TAI-UTC offset in seconds after that instant, and the calendar label.
///
/// Fields are private so a consumer can
/// never observe a non-finite MJD/offset (see [`LeapSecondEntry::new`]);
/// the year/month/day labels are carried from the asset, exactly
/// as the source stores them without validation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LeapSecondEntry {
    mjd: f64,
    tai_utc: f64,
    year: i32,
    month: i32,
    day: i32,
}

impl LeapSecondEntry {
    /// Build a table row, rejecting non-finite MJD or offset.
    pub fn new(mjd: f64, tai_utc: f64, year: i32, month: i32, day: i32) -> Result<Self> {
        if !mjd.is_finite() {
            return Err(Error::NonFiniteArgument { value: mjd });
        }
        if !tai_utc.is_finite() {
            return Err(Error::NonFiniteArgument { value: tai_utc });
        }
        Ok(Self {
            mjd,
            tai_utc,
            year,
            month,
            day,
        })
    }

    /// MJD from which this offset applies.
    pub fn mjd(self) -> f64 {
        self.mjd
    }

    /// TAI-UTC in seconds from this row's MJD onward.
    pub fn tai_utc(self) -> f64 {
        self.tai_utc
    }

    /// Calendar year label of this row.
    pub fn year(self) -> i32 {
        self.year
    }

    /// Calendar month label of this row.
    pub fn month(self) -> i32 {
        self.month
    }

    /// Calendar day label of this row.
    pub fn day(self) -> i32 {
        self.day
    }
}

/// Packaged leap-second table: `(year, month, day, TAI-UTC after)`.
///
/// Returns an owned ascending table of effective-date rows, from 1972-01-01
/// through 2017-01-01 (TAI−UTC = 37 s). The MJD column is derived with
/// [`calendar_to_mjd`]. Future leap seconds require newer data; the final
/// offset is not a prediction that no further leap second will occur.
pub fn hardcoded_leap_seconds() -> Vec<LeapSecondEntry> {
    const ROWS: [(i32, i32, i32, f64); 28] = [
        (1972, 1, 1, 10.0),
        (1972, 7, 1, 11.0),
        (1973, 1, 1, 12.0),
        (1974, 1, 1, 13.0),
        (1975, 1, 1, 14.0),
        (1976, 1, 1, 15.0),
        (1977, 1, 1, 16.0),
        (1978, 1, 1, 17.0),
        (1979, 1, 1, 18.0),
        (1980, 1, 1, 19.0),
        (1981, 7, 1, 20.0),
        (1982, 7, 1, 21.0),
        (1983, 7, 1, 22.0),
        (1985, 7, 1, 23.0),
        (1988, 1, 1, 24.0),
        (1990, 1, 1, 25.0),
        (1991, 1, 1, 26.0),
        (1992, 7, 1, 27.0),
        (1993, 7, 1, 28.0),
        (1994, 7, 1, 29.0),
        (1996, 1, 1, 30.0),
        (1997, 7, 1, 31.0),
        (1999, 1, 1, 32.0),
        (2006, 1, 1, 33.0),
        (2009, 1, 1, 34.0),
        (2012, 7, 1, 35.0),
        (2015, 7, 1, 36.0),
        (2017, 1, 1, 37.0),
    ];
    ROWS.iter()
        .map(|&(year, month, day, tai_utc)| {
            // Inputs are small finite literals; construction is infallible.
            LeapSecondEntry::new(calendar_to_mjd(year, month, day), tai_utc, year, month, day)
                .expect("hardcoded leap-second row must be finite")
        })
        .collect()
}

/// Parse IERS Bulletin C leap-second content (`Leap_Second.dat`).
///
/// Accepts the all-numeric `MJD day month year TAI-UTC` rows and the
/// tolerated textual `year month-name day ... TAI-UTC[= Ns]` variant, in
/// that priority. Comment (`#`), blank, unparsable and nonfinite MJD/offset
/// rows are skipped. Rows retain input order; use [`parse_leap_seconds`] for
/// sorted output and the packaged empty-input fallback.
pub fn parse_leap_seconds_iers(text: &str) -> Vec<LeapSecondEntry> {
    const MONTHS: [&str; 12] = [
        "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec",
    ];
    let mut entries = Vec::new();
    for raw_line in text.trim().split('\n') {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 5 {
            continue;
        }
        // Numeric row iff columns 2-4 are all (ASCII) digits, as in the
        // source's `parts[1/2/3].isdigit` gate (asset files are ASCII).
        let numeric = parts[1].bytes().all(|b| b.is_ascii_digit())
            && !parts[1].is_empty()
            && parts[2].bytes().all(|b| b.is_ascii_digit())
            && !parts[2].is_empty()
            && parts[3].bytes().all(|b| b.is_ascii_digit())
            && !parts[3].is_empty();
        let row = if numeric {
            let parsed = (
                parts[0].parse::<f64>(),
                parts[1].parse::<i32>(),
                parts[2].parse::<i32>(),
                parts[3].parse::<i32>(),
                parts[4].parse::<f64>(),
            );
            match parsed {
                (Ok(mjd), Ok(day), Ok(month), Ok(year), Ok(tai_utc)) => {
                    LeapSecondEntry::new(mjd, tai_utc, year, month, day).ok()
                }
                _ => None,
            }
        } else {
            let year = parts[0].parse::<i32>().ok();
            let month_name = parts[1].to_ascii_lowercase();
            let month = MONTHS
                .iter()
                .position(|&m| month_name.get(..3) == Some(m))
                .map(|i| (i as i32) + 1)
                .unwrap_or(0);
            if month == 0 {
                continue;
            }
            let day = parts[2].parse::<i32>().ok();
            let tai_utc = parts[4].trim_end_matches('s').parse::<f64>().ok();
            match (year, day, tai_utc) {
                (Some(year), Some(day), Some(tai_utc)) => LeapSecondEntry::new(
                    calendar_to_mjd(year, month, day),
                    tai_utc,
                    year,
                    month,
                    day,
                )
                .ok(),
                _ => None,
            }
        };
        if let Some(entry) = row {
            entries.push(entry);
        }
    }
    entries
}

/// Parse IETF `leap-seconds.list` content (NTP timestamps).
///
/// Each data row is `ntp_seconds tai_utc`; MJD is
/// `15020.0 + ntp_seconds / 86400.0` with the NTP epoch of 1900-01-01.
/// Comment/blank/unparsable lines are skipped, as in the source.
pub fn parse_leap_seconds_ietf(text: &str) -> Vec<LeapSecondEntry> {
    /// MJD of the NTP epoch (1900-01-01).
    const NTP_EPOCH_MJD: f64 = 15020.0;
    let mut entries = Vec::new();
    for raw_line in text.trim().split('\n') {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 2 {
            continue;
        }
        let parsed = (parts[0].parse::<i64>(), parts[1].parse::<i64>());
        if let (Ok(ntp_seconds), Ok(tai_utc_int)) = parsed {
            // Integer-to-binary64 conversion rounds to nearest;
            // NTP/leap values in any real asset are far below 2^53.
            let tai_utc = tai_utc_int as f64;
            let mjd = NTP_EPOCH_MJD + ntp_seconds as f64 / 86400.0;
            if !mjd.is_finite() {
                continue;
            }
            // `mjd_to_calendar` fails only far outside `i32`; such a row
            // is skipped rather than admitted with a wrong label.
            if let Ok((year, month, day)) = mjd_to_calendar(mjd) {
                if let Ok(entry) = LeapSecondEntry::new(mjd, tai_utc, year, month, day) {
                    entries.push(entry);
                }
            }
        }
    }
    entries
}

/// Parse leap-second asset content in either supported format.
///
/// Selects IETF for text containing `NTP`, `#$`, `#@` or `leap-seconds.list`,
/// otherwise IERS. If the selected parser finds nothing, tries the other; if
/// both find nothing, returns [`hardcoded_leap_seconds`]. Malformed content is
/// therefore not an error here. Rows are stably sorted by ascending MJD.
/// Only supplied text is parsed; no filesystem or network access occurs.
pub fn parse_leap_seconds(text: &str) -> Vec<LeapSecondEntry> {
    let is_ietf = text.contains("NTP")
        || text.contains("#$")
        || text.contains("#@")
        || text.contains("leap-seconds.list");
    let (first, second) = if is_ietf {
        (parse_leap_seconds_ietf(text), parse_leap_seconds_iers(text))
    } else {
        (parse_leap_seconds_iers(text), parse_leap_seconds_ietf(text))
    };
    let mut entries = if first.is_empty() { second } else { first };
    if entries.is_empty() {
        entries = hardcoded_leap_seconds();
    }
    // All rows are finite by construction; `total_cmp` keeps the sort
    // stable for duplicate MJDs.
    entries.sort_by(|a, b| a.mjd.total_cmp(&b.mjd));
    entries
}

/// Convert a calendar date to Modified Julian Date.
///
/// Uses a Gregorian Meeus recurrence at midnight with truncation toward zero,
/// not the floor operations of [`julday`]. No reform switch or civil-field
/// validation is applied. The result is scale-free and uses
/// `MJD = JD - 2400000.5`. Before the reform, [`mjd_to_calendar`] follows a
/// different calendar convention, so these two helpers are not universal
/// inverses. For explicit proleptic-calendar conversion use [`julday`] and
/// [`jd_to_mjd`].
pub fn calendar_to_mjd(year: i32, month: i32, day: i32) -> f64 {
    let (mut y, mut m) = (year as f64, month as f64);
    if m <= 2.0 {
        y -= 1.0;
        m += 12.0;
    }
    let a = (y / 100.0).trunc();
    let b = 2.0 - a + (a / 4.0).trunc();
    (365.25 * (y + 4716.0)).trunc() + (30.6001 * (m + 1.0)).trunc() + day as f64 + b
        - 1524.5
        - 2_400_000.5
}

/// Convert a Modified Julian Date to a `(year, month, day)` calendar date.
///
/// Uses truncation toward zero and a historical reform switch: shifted
/// integer JD `z < 2299161` selects Julian, otherwise Gregorian. This differs
/// from the always-Gregorian [`calendar_to_mjd`] recurrence and is not its
/// inverse for all dates. The fractional day is discarded.
/// Nonfinite input returns [`Error::NonFiniteArgument`]; derived fields
/// outside `i32` return [`Error::IntegerOverflow`].
pub fn mjd_to_calendar(mjd: f64) -> Result<(i32, i32, i32)> {
    if !mjd.is_finite() {
        return Err(Error::NonFiniteArgument { value: mjd });
    }
    let jd = mjd + 2_400_000.5 + 0.5;
    let z = jd.trunc();
    let f = jd - z;
    let a = if z < 2_299_161.0 {
        z
    } else {
        let alpha = ((z - 1_867_216.25) / 36_524.25).trunc();
        z + 1.0 + alpha - (alpha / 4.0).trunc()
    };
    let b = a + 1524.0;
    let c = ((b - 122.1) / 365.25).trunc();
    let d = (365.25 * c).trunc();
    let e = ((b - d) / 30.6001).trunc();
    let day = (b - d - (30.6001 * e).trunc() + f).trunc();
    let month = if e < 14.0 { e - 1.0 } else { e - 13.0 };
    let year = if month > 2.0 { c - 4716.0 } else { c - 4715.0 };
    Ok((
        float_to_i32(year)?,
        float_to_i32(month)?,
        float_to_i32(day)?,
    ))
}

/// Convert Julian Day to Modified Julian Date (`MJD = JD - 2400000.5`).
///
/// , with the module's non-finite strictness.
pub fn jd_to_mjd(jd: f64) -> Result<f64> {
    if !jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd });
    }
    Ok(jd - 2_400_000.5)
}

/// Convert Modified Julian Date to Julian Day (`JD = MJD + 2400000.5`).
///
/// , with the module's non-finite strictness.
pub fn mjd_to_jd(mjd: f64) -> Result<f64> {
    if !mjd.is_finite() {
        return Err(Error::NonFiniteArgument { value: mjd });
    }
    Ok(mjd + 2_400_000.5)
}

/// IERS rate-based UTC coefficients for 1961–1971.
///
/// Each row is `(valid_from_mjd, offset_seconds, reference_mjd, seconds_per_day)`.
/// Select the last row whose start does not exceed the query MJD and evaluate
/// `offset + (mjd - reference_mjd) * rate`. Earlier dates clamp to the first
/// offset; UTC is physically undefined before its introduction. Modern UTC
/// requires the separately acquired leap-second table. The immutable rows also
/// are shared by all native pre-1972 UTC conversions.
pub const PRE_1972_TAI_UTC: [(f64, f64, f64, f64); 13] = [
    (37300.0, 1.4228180, 37300.0, 0.0012960),
    (37512.0, 1.3728180, 37300.0, 0.0012960),
    (37665.0, 1.8458580, 37665.0, 0.0011232),
    (38334.0, 1.9458580, 37665.0, 0.0011232),
    (38395.0, 3.2401300, 38761.0, 0.0012960),
    (38486.0, 3.3401300, 38761.0, 0.0012960),
    (38639.0, 3.4401300, 38761.0, 0.0012960),
    (38761.0, 3.5401300, 38761.0, 0.0012960),
    (38820.0, 3.6401300, 38761.0, 0.0012960),
    (38942.0, 3.7401300, 38761.0, 0.0012960),
    (39004.0, 3.8401300, 38761.0, 0.0012960),
    (39126.0, 4.3131700, 39126.0, 0.0025920),
    (39887.0, 4.2131700, 39126.0, 0.0025920),
];

/// TAI-UTC in seconds for a pre-1972 MJD via the rate-based model.
///
/// Dates before 1961-01-01 (before UTC existed) clamp to the earliest
/// tabulated offset; the value is an approximation there and TAI-UTC is
/// not physically defined. exactly,
/// including the clamp form `offset + (first_valid_from - ref) * rate`.
pub(crate) fn pre_1972_tai_utc(mjd: f64) -> f64 {
    let (first_from, first_off, first_ref, first_rate) = PRE_1972_TAI_UTC[0];
    if mjd < first_from {
        return first_off + (first_from - first_ref) * first_rate;
    }
    let mut selected = (first_off, first_ref, first_rate);
    for &(valid_from, off, rf, rate) in PRE_1972_TAI_UTC.iter() {
        if mjd >= valid_from {
            selected = (off, rf, rate);
        } else {
            break;
        }
    }
    selected.0 + (mjd - selected.1) * selected.2
}

/// TAI-UTC in seconds for a Modified Julian Date against `table`.
///
/// MJDs before [`FIRST_LEAP_SECOND_MJD`]
/// use the rate-based model independently of the table; an empty table
/// falls back to `37.0` at/after MJD 57754.0 (2017-01-01) and `32.0`
/// otherwise; otherwise the last row with `row.mjd <= mjd` wins (the
/// table is expected ascending, as [`parse_leap_seconds`] produces).
/// A nonempty table with no row at/before the query returns zero. After its
/// last row the final offset is retained, without a freshness check. This
/// lookup does not validate whether supplied offsets describe real leap events.
pub fn tai_utc_for_mjd(mjd: f64, table: &[LeapSecondEntry]) -> Result<f64> {
    if !mjd.is_finite() {
        return Err(Error::NonFiniteArgument { value: mjd });
    }
    if mjd < FIRST_LEAP_SECOND_MJD {
        return Ok(pre_1972_tai_utc(mjd));
    }
    if table.is_empty() {
        return Ok(if mjd >= EMPTY_TABLE_MODERN_MJD {
            37.0
        } else {
            32.0
        });
    }
    let mut tai_utc = 0.0;
    for row in table {
        if mjd >= row.mjd {
            tai_utc = row.tai_utc;
        } else {
            break;
        }
    }
    Ok(tai_utc)
}

/// TAI-UTC in seconds for a Julian Day against `table`.
///
/// Computes `jd - 2400000.5` then applies [`tai_utc_for_mjd`], including its
/// pre-1972 model, empty-table policy and final-offset retention. The supplied
/// day is the UTC calendar-day coordinate used for table lookup, not UT1.
pub fn tai_utc_for_jd(jd: f64, table: &[LeapSecondEntry]) -> Result<f64> {
    tai_utc_for_mjd(jd_to_mjd(jd)?, table)
}

/// Convert a UTC calendar label to Julian Day in TAI.
///
/// Forms decimal hours, obtains an ordinary calendar JD with [`julday`],
/// then adds the supplied table's TAI−UTC offset in days. Calendar and clock
/// bounds and leap-label validity are not checked here. A `60.x` label flows
/// into the calendar JD just like overflowing seconds, so it can share a
/// coordinate with the following day. This is not a leap-second-preserving
/// quasi-JD conversion and is not generally invertible at a leap boundary.
/// [`utc_to_jd_with_config`] adds civil-label admission, with its own stated
/// era and representation limitations. Nonfinite/overflow errors propagate
/// from the calendar and table helpers.
///
/// Worked shape: `2020-01-01 00:00:00`
/// UTC converts with TAI 37 s ahead (the modern packaged-table anchor).
///
/// ```rust
/// use libephemeris_core::time::{Calendar, hardcoded_leap_seconds, julday, utc_to_tai_jd};
///
/// let table = hardcoded_leap_seconds();
/// let tai = utc_to_tai_jd(2020, 1, 1, 0, 0, 0.0, Calendar::Gregorian, &table).unwrap();
/// let utc = julday(2020, 1, 1, 0.0, Calendar::Gregorian).unwrap();
/// assert_eq!(tai, utc + 37.0 / 86400.0);
/// ```
#[allow(clippy::too_many_arguments)] // positional arity mirrors the model
pub fn utc_to_tai_jd(
    year: i32,
    month: i32,
    day: i32,
    hour: i32,
    minute: i32,
    second: f64,
    cal: Calendar,
    table: &[LeapSecondEntry],
) -> Result<f64> {
    let decimal_hour = hour as f64 + minute as f64 / 60.0 + second / 3600.0;
    let jd_utc = julday(year, month, day, decimal_hour, cal)?;
    Ok(jd_utc + tai_utc_for_jd(jd_utc, table)? / 86400.0)
}

/// Convert a Julian Day in TAI to UTC calendar fields.
///
/// Looks up TAI−UTC at `jd_tai`, subtracts it, re-looks up at that estimate,
/// then subtracts again. [`revjul`] and truncation split the resulting ordinary
/// calendar JD in `cal`. This does not reconstruct an unambiguous `23:59:60`
/// label; scalar rounding also limits fractional seconds. Table/calendar
/// failures propagate, including nonfinite epochs and unrepresentable fields.
pub fn tai_jd_to_utc(jd_tai: f64, cal: Calendar, table: &[LeapSecondEntry]) -> Result<UtcFields> {
    let first = tai_utc_for_jd(jd_tai, table)?;
    let jd_approx = jd_tai - first / 86400.0;
    let refined = tai_utc_for_jd(jd_approx, table)?;
    let jd_utc = jd_tai - refined / 86400.0;
    let date = revjul(jd_utc, cal)?;
    // `revjul` guarantees `hour in [0, 24)`; truncation mirrors `int`.
    let hour = date.hour.trunc() as i32;
    let minute_frac = (date.hour - hour as f64) * 60.0;
    let minute = minute_frac.trunc() as i32;
    let second = (minute_frac - minute as f64) * 60.0;
    Ok(UtcFields {
        year: date.year,
        month: date.month,
        day: date.day,
        hour,
        minute,
        second,
    })
}

/// Convert Julian Day in TT to Julian Day in TAI (`TAI = TT - 32.184 s`).
///
/// The defining TT−TAI offset is exactly 32.184 SI seconds. Its binary64 days
/// representation and scalar subtraction round; this is not an exact
/// arbitrary-epoch round-trip guarantee. [`tt_to_tai_two_part`] retains a small
/// offset separately. Nonfinite input returns [`Error::NonFiniteArgument`].
///
/// ```rust
/// use libephemeris_core::time::{TT_TAI_OFFSET_DAYS, tt_to_tai_jd};
///
/// assert_eq!(tt_to_tai_jd(2451545.0).unwrap(), 2451545.0 - TT_TAI_OFFSET_DAYS);
/// assert_eq!(tt_to_tai_jd(2451545.0).unwrap(), 2451545.0 - 32.184 / 86400.0);
/// ```
pub fn tt_to_tai_jd(jd_tt: f64) -> Result<f64> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    Ok(jd_tt - TT_TAI_OFFSET_DAYS)
}

/// Convert Julian Day in TAI to Julian Day in TT (`TT = TAI + 32.184 s`).
///
/// Applies the defining 32.184-SI-second offset using rounded binary64 days.
/// Scalar subtraction/addition can discard information, so an arbitrary
/// [`tt_to_tai_jd`] round trip is not guaranteed exact. The J2000 example below
/// is exact at that input. Nonfinite input returns [`Error::NonFiniteArgument`].
///
/// ```rust
/// use libephemeris_core::time::{tai_to_tt_jd, tt_to_tai_jd};
///
/// let tai = tt_to_tai_jd(2451545.0).unwrap();
/// assert_eq!(tai_to_tt_jd(tai).unwrap(), 2451545.0);
/// ```
pub fn tai_to_tt_jd(jd_tai: f64) -> Result<f64> {
    if !jd_tai.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tai });
    }
    Ok(jd_tai + TT_TAI_OFFSET_DAYS)
}

/// Convert Julian Day in TAI to Julian Day in UTC.
///
/// Same two-step refinement as [`tai_jd_to_utc`], returning its ordinary
/// 86400-second calendar-day coordinate. The result is not an unambiguous
/// leap-label encoding; input/table errors propagate from that lookup chain.
pub fn tai_to_utc_jd(jd_tai: f64, table: &[LeapSecondEntry]) -> Result<f64> {
    let first = tai_utc_for_jd(jd_tai, table)?;
    let jd_approx = jd_tai - first / 86400.0;
    let refined = tai_utc_for_jd(jd_approx, table)?;
    Ok(jd_tai - refined / 86400.0)
}

/// Convert an ordinary UTC calendar JD to TAI (`TAI = UTC + TAI-UTC`).
///
/// Uses [`tai_utc_for_jd`] once and inherits its table/finiteness contract.
/// No leap-second label is represented separately from the scalar day, so this
/// is not a bijection through positive leap seconds.
pub fn utc_jd_to_tai(jd_utc: f64, table: &[LeapSecondEntry]) -> Result<f64> {
    Ok(jd_utc + tai_utc_for_jd(jd_utc, table)? / 86400.0)
}

// Two-part TAI/TT scale conversions.
//
// TwoPartEpoch retains small scale offsets without rounding them to a full
// Julian-day ULP. Leap-table selection and calendar output use the collapsed
// JD, as in the scalar functions. Retaining an offset does not repair the
// ordinary UTC JD representation's loss of leap-label identity.

/// TAI-UTC in seconds for a two-part Julian Day against `table`.
///
/// Row selection keys on the collapsed JD, matching [`tai_utc_for_jd`]
/// exactly for split inputs.
pub fn tai_utc_for_two_part(jd: TwoPartEpoch, table: &[LeapSecondEntry]) -> Result<f64> {
    tai_utc_for_jd(jd.value(), table)
}

/// Convert a two-part Julian Day in TT to TAI (`TAI = TT - 32.184 s`).
///
/// Same exact offset as [`tt_to_tai_jd`], applied via [`TwoPartEpoch::add_offset`].
pub fn tt_to_tai_two_part(jd_tt: TwoPartEpoch) -> Result<TwoPartEpoch> {
    jd_tt.add_offset(-TT_TAI_OFFSET_DAYS)
}

/// Convert a two-part Julian Day in TAI to TT (`TT = TAI + 32.184 s`).
///
/// Same exact offset as [`tai_to_tt_jd`], applied via [`TwoPartEpoch::add_offset`].
pub fn tai_to_tt_two_part(jd_tai: TwoPartEpoch) -> Result<TwoPartEpoch> {
    jd_tai.add_offset(TT_TAI_OFFSET_DAYS)
}

/// Convert a two-part Julian Day in UTC to TAI (`TAI = UTC + TAI-UTC`).
///
/// Same single lookup as [`utc_jd_to_tai`]; the offset is applied via
/// [`TwoPartEpoch::add_offset`].
pub fn utc_jd_to_tai_two_part(
    jd_utc: TwoPartEpoch,
    table: &[LeapSecondEntry],
) -> Result<TwoPartEpoch> {
    let offset_days = tai_utc_for_two_part(jd_utc, table)? / 86400.0;
    jd_utc.add_offset(offset_days)
}

/// Convert a two-part Julian Day in TAI to Julian Day in UTC.
///
/// Same two-step refinement as [`tai_to_utc_jd`], with the refinement keys
/// computed from the collapsed JD exactly as the scalar path computes them;
/// the final offset is applied via [`TwoPartEpoch::add_offset`].
pub fn tai_to_utc_jd_two_part(
    jd_tai: TwoPartEpoch,
    table: &[LeapSecondEntry],
) -> Result<TwoPartEpoch> {
    let collapsed = jd_tai.value();
    let first = tai_utc_for_jd(collapsed, table)?;
    let approx = collapsed - first / 86400.0;
    let refined = tai_utc_for_jd(approx, table)?;
    jd_tai.add_offset(-refined / 86400.0)
}

/// Convert a UTC calendar label to a two-part Julian Day in TAI.
///
/// Same `decimal_hour`/`julday`/lookup construction as [`utc_to_tai_jd`]
/// (identical rounded intermediate JD), then a split with the offset applied
/// via [`TwoPartEpoch::add_offset`]. The split cannot recover earlier calendar
/// rounding or leap-label information; the scalar helper's limitations apply.
#[allow(clippy::too_many_arguments)] // positional arity mirrors the model
pub fn utc_to_tai_two_part(
    year: i32,
    month: i32,
    day: i32,
    hour: i32,
    minute: i32,
    second: f64,
    cal: Calendar,
    table: &[LeapSecondEntry],
) -> Result<TwoPartEpoch> {
    let decimal_hour = hour as f64 + minute as f64 / 60.0 + second / 3600.0;
    let jd_utc = julday(year, month, day, decimal_hour, cal)?;
    let offset_days = tai_utc_for_jd(jd_utc, table)? / 86400.0;
    TwoPartEpoch::split(jd_utc)?.add_offset(offset_days)
}

/// Convert a two-part Julian Day in TAI to UTC calendar fields.
///
/// Collapses the input with [`TwoPartEpoch::value`] and delegates to
/// [`tai_jd_to_utc`]. Sub-JD-ulp corrections are therefore rounded before
/// reconstructing seconds; two-part input adds no calendar-output precision
/// here. The scalar helper's leap-label and error limitations apply.
pub fn tai_two_part_to_utc(
    jd_tai: TwoPartEpoch,
    cal: Calendar,
    table: &[LeapSecondEntry],
) -> Result<UtcFields> {
    tai_jd_to_utc(jd_tai.value(), cal, table)
}

// Delta-T models, user override and tidal adjustment.
//
// TT-minus-UT1 is returned in days; polynomial coefficients and tidal
// adjustments use seconds. The priority is an explicit user override, then
// enabled observed IERS data, then the selected SMH2016 or Espenak-Meeus model.
// A user override is never tidally adjusted. The tidal correction depends
// quadratically on centuries from 1955.5 and is zero after that reference epoch.
// Model selection is functional and does not mutate the caller configuration.

/// Delta-T model selector, mirroring `set_delta_t_model` names.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DeltaTModel {
    /// Stephenson-Morrison-Hohenkerk 2016, the model default, realized
    /// there through its ephemeris-backend integration (see [`DeltaTConfig`]).
    Smh2016,
    /// Classic NASA Espenak & Meeus (2006) polynomial set, natively
    /// reimplemented here from the published coefficients.
    EspenakMeeus,
}

impl DeltaTModel {
    /// Canonical model name (`"smh2016"` / `"espenak_meeus"`).
    pub fn as_name(self) -> &'static str {
        match self {
            Self::Smh2016 => "smh2016",
            Self::EspenakMeeus => "espenak_meeus",
        }
    }

    /// Parse a model name with the `set_delta_t_model`
    /// normalization (ASCII case/whitespace-insensitive). Returns `None`
    /// for unknown names; the language boundary maps that to
    /// `ValueError`.
    pub fn from_name(name: &str) -> Option<Self> {
        match name.trim().to_ascii_lowercase().as_str() {
            "smh2016" => Some(Self::Smh2016),
            "espenak_meeus" => Some(Self::EspenakMeeus),
            _ => None,
        }
    }
}

/// Julian Day of the 1955.5 reference epoch for the tidal-acceleration
/// adjustment, derived exactly as in the source (`2451545.0 - 44.5*365.25`,
/// i.e. 44.5 Julian years before J2000; Morrison & Stephenson 2004).
pub const TID_ACC_REFERENCE_JD: f64 = 2451545.0 - 44.5 * 365.25;

/// IERS Conventions (2010), TN 36, Eq. 5.43 linear rates in arcsec per
/// Julian century: `F = L - Omega` and `Omega`. Their sum is the mean
/// lunar longitude rate used to convert a lunar secular-acceleration
/// difference into clock seconds.
const IERS_F_RATE_ARCSEC_CY: f64 = 1739527262.8478;
const IERS_OMEGA_RATE_ARCSEC_CY: f64 = -6962890.5431;

/// Delta T in **seconds** via the Espenak & Meeus (2006) polynomial set
/// (NASA Five Millennium Canon, roughly -1999 to +3000 with a parabolic
/// tail outside).
///
/// `decimal_year` is a Gregorian decimal year. [`deltat_with_config`] uses the
/// monthly midpoint convention `year + (month - 0.5) / 12`. Exact branch
/// thresholds belong to the later branch. This unchecked formula returns
/// seconds; nonfinite inputs and arithmetic overflow propagate through IEEE
/// arithmetic. The configured entry checks its epoch unless a user override
/// returns first. A computed finite value is not a certified prediction.
pub fn deltat_espenak_meeus(decimal_year: f64) -> f64 {
    let y = decimal_year;
    if y < -500.0 {
        let u = (y - 1820.0) / 100.0;
        return -20.0 + 32.0 * u * u;
    }
    if y < 500.0 {
        let u = y / 100.0;
        return 10583.6 - 1014.41 * u + 33.78311 * u.powf(2.0)
            - 5.952053 * u.powf(3.0)
            - 0.1798452 * u.powf(4.0)
            + 0.022174192 * u.powf(5.0)
            + 0.0090316521 * u.powf(6.0);
    }
    if y < 1600.0 {
        let u = (y - 1000.0) / 100.0;
        return 1574.2 - 556.01 * u + 71.23472 * u.powf(2.0) + 0.319781 * u.powf(3.0)
            - 0.8503463 * u.powf(4.0)
            - 0.005050998 * u.powf(5.0)
            + 0.0083572073 * u.powf(6.0);
    }
    if y < 1700.0 {
        let t = y - 1600.0;
        return 120.0 - 0.9808 * t - 0.01532 * t.powf(2.0) + t.powf(3.0) / 7129.0;
    }
    if y < 1800.0 {
        let t = y - 1700.0;
        return 8.83 + 0.1603 * t - 0.0059285 * t.powf(2.0) + 0.00013336 * t.powf(3.0)
            - t.powf(4.0) / 1174000.0;
    }
    if y < 1860.0 {
        let t = y - 1800.0;
        return 13.72 - 0.332447 * t + 0.0068612 * t.powf(2.0) + 0.0041116 * t.powf(3.0)
            - 0.00037436 * t.powf(4.0)
            + 0.0000121272 * t.powf(5.0)
            - 0.0000001699 * t.powf(6.0)
            + 0.000000000875 * t.powf(7.0);
    }
    if y < 1900.0 {
        let t = y - 1860.0;
        return 7.62 + 0.5737 * t - 0.251754 * t.powf(2.0) + 0.01680668 * t.powf(3.0)
            - 0.0004473624 * t.powf(4.0)
            + t.powf(5.0) / 233174.0;
    }
    if y < 1920.0 {
        let t = y - 1900.0;
        return -2.79 + 1.494119 * t - 0.0598939 * t.powf(2.0) + 0.0061966 * t.powf(3.0)
            - 0.000197 * t.powf(4.0);
    }
    if y < 1941.0 {
        let t = y - 1920.0;
        return 21.20 + 0.84493 * t - 0.076100 * t.powf(2.0) + 0.0020936 * t.powf(3.0);
    }
    if y < 1961.0 {
        let t = y - 1950.0;
        return 29.07 + 0.407 * t - t.powf(2.0) / 233.0 + t.powf(3.0) / 2547.0;
    }
    if y < 1986.0 {
        let t = y - 1975.0;
        return 45.45 + 1.067 * t - t.powf(2.0) / 260.0 - t.powf(3.0) / 718.0;
    }
    if y < 2005.0 {
        let t = y - 2000.0;
        return 63.86 + 0.3345 * t - 0.060374 * t.powf(2.0)
            + 0.0017275 * t.powf(3.0)
            + 0.000651814 * t.powf(4.0)
            + 0.00002373599 * t.powf(5.0);
    }
    if y < 2050.0 {
        let t = y - 2000.0;
        return 62.92 + 0.32217 * t + 0.005589 * t.powf(2.0);
    }
    if y < 2150.0 {
        return -20.0 + 32.0 * ((y - 1820.0) / 100.0).powf(2.0) - 0.5628 * (2150.0 - y);
    }
    let u = (y - 1820.0) / 100.0;
    -20.0 + 32.0 * u * u
}

/// Delta-T adjustment (seconds) for a non-default tidal acceleration.
///
/// Derived from the IERS 2010 mean-lunar-longitude rate: linear in
/// `tid_acc - TIDAL_DEFAULT`, quadratic in elapsed Julian centuries from
/// [`TID_ACC_REFERENCE_JD`]. Exactly `0.0` at the default acceleration or
/// for epochs on/after the reference (modern observed Delta T is never
/// adjusted). `tjdut` is UT1 Julian days and `tid_acc` arcsec/Julian-century².
/// This unchecked helper uses its zero/reference guards before arithmetic;
/// nonfinite operands are not universally rejected. [`deltat_with_config`]
/// checks its non-overridden epoch, not every public configuration field.
pub fn tid_acc_adjustment_seconds(tjdut: f64, tid_acc: f64) -> f64 {
    let dacc = tid_acc - TIDAL_DEFAULT;
    if dacc == 0.0 || tjdut >= TID_ACC_REFERENCE_JD {
        return 0.0;
    }
    let u = (tjdut - TID_ACC_REFERENCE_JD) / 36525.0;
    tid_acc_adjustment_from_parts(dacc, u)
}

/// Clock-seconds coefficient for an acquired lunar acceleration difference.
/// The IERS longitude rates are added before conversion to arcseconds per day;
/// conversion to seconds and multiplication by minus one half retain that
/// association. The coefficient multiplies arcseconds per century squared
/// and squared elapsed Julian centuries to produce clock seconds.
pub fn tid_acc_adjustment_coefficient() -> f64 {
    let lunar_motion_arcsec_day = (IERS_F_RATE_ARCSEC_CY + IERS_OMEGA_RATE_ARCSEC_CY) / 36525.0;
    let seconds_per_arcsec = 86400.0 / lunar_motion_arcsec_day;
    -0.5 * seconds_per_arcsec
}

/// Finish a tidal adjustment from already consumed difference and centuries.
/// `dacc` is acceleration relative to the default, in arcsec/century squared;
/// `centuries` is elapsed Julian centuries from the reference epoch. The
/// caller owns subtraction, the zero/reference guards and century conversion.
/// Left-associated multiplication preserves the scalar clock arithmetic;
/// nonfinite operands propagate according to IEEE-754 without a new guard.
pub fn tid_acc_adjustment_from_parts(dacc: f64, centuries: f64) -> f64 {
    tid_acc_adjustment_coefficient() * dacc * centuries * centuries
}

/// Explicit-parameter Delta-T request (caller-supplied backends, no globals).
///
/// Carries resolved model, user override, IERS and tidal settings without
/// reading process configuration itself. Optional callbacks receive UT1 Julian
/// days and return TT−UT1 **seconds**, or `None` for unavailable coverage.
/// They are borrowed for `'a`, run synchronously when selected, and may acquire
/// their own data. Callback panics propagate; missing coverage should be `None`.
///
/// Public fields are not validated as a group. In particular, an override or
/// callback can supply nonfinite values, and this value type does not promise
/// pure computation, finite output or thread-safe callback implementations.
pub struct DeltaTConfig<'a> {
    /// Selected model (`smh2016` needs [`DeltaTConfig::smh2016_lookup_seconds`]).
    pub model: DeltaTModel,
    /// Hard user override in **days**, returned and never tidally
    /// adjusted.
    pub user_override_days: Option<f64>,
    /// Effective tidal acceleration in arcsec/century^2 for this
    /// computation (`TIDAL_DEFAULT` unless overridden/flag-selected).
    pub tidal_acceleration: f64,
    /// Whether observed IERS data outranks the selected model.
    pub iers_enabled: bool,
    /// Observed-IERS Delta-T provider (seconds), when `iers_enabled`.
    pub iers_lookup_seconds: Option<&'a dyn Fn(f64) -> Option<f64>>,
    /// SMH-2016 Delta-T provider (seconds) backing the default model.
    pub smh2016_lookup_seconds: Option<&'a dyn Fn(f64) -> Option<f64>>,
}

impl Default for DeltaTConfig<'_> {
    /// the model global defaults: `smh2016`, no override, default tidal
    /// acceleration, IERS disabled, no providers. A bare default carries
    /// no backends, so selecting `smh2016` on it reports
    /// [`Error::DeltaTModelUnavailable`]; contexts wire the vendored
    /// table provider through
    /// [`EphemerisContext::delta_t_config`](crate::context::EphemerisContext).
    fn default() -> Self {
        Self {
            model: DeltaTModel::Smh2016,
            user_override_days: None,
            tidal_acceleration: TIDAL_DEFAULT,
            iers_enabled: false,
            iers_lookup_seconds: None,
            smh2016_lookup_seconds: None,
        }
    }
}

/// Calculate Delta T (TT minus UT1) in **days** for a UT1 Julian Day.
///
/// A stored user override wins before epoch validation and is returned
/// unchanged, even if either value is nonfinite. Otherwise the epoch must be
/// finite ([`Error::NonFiniteArgument`]). Covering enabled IERS data wins over
/// the selected model; tidal adjustment is added to computed seconds before
/// division by 86400. An override is never tidally adjusted.
///
/// `None` from IERS falls through to the model. Selected SMH2016 with no usable
/// provider returns [`Error::DeltaTModelUnavailable`]; Espenak–Meeus uses
/// Gregorian calendar fields and can return calendar overflow errors. Callback
/// results, tidal settings and computed output are not independently checked
/// for finiteness. Callback I/O and panics retain their ordinary behavior.
pub fn deltat_with_config(tjdut: f64, cfg: &DeltaTConfig<'_>) -> Result<f64> {
    if let Some(user) = cfg.user_override_days {
        return Ok(user);
    }
    if !tjdut.is_finite() {
        return Err(Error::NonFiniteArgument { value: tjdut });
    }
    let adjust_seconds = tid_acc_adjustment_seconds(tjdut, cfg.tidal_acceleration);
    if cfg.iers_enabled {
        if let Some(lookup) = cfg.iers_lookup_seconds {
            if let Some(delta_t_seconds) = lookup(tjdut) {
                return Ok((delta_t_seconds + adjust_seconds) / 86400.0);
            }
        }
    }
    match cfg.model {
        DeltaTModel::EspenakMeeus => {
            let date = revjul(tjdut, Calendar::Gregorian)?;
            let decimal_year = date.year as f64 + (f64::from(date.month) - 0.5) / 12.0;
            Ok((deltat_espenak_meeus(decimal_year) + adjust_seconds) / 86400.0)
        }
        DeltaTModel::Smh2016 => match cfg.smh2016_lookup_seconds {
            Some(lookup) => match lookup(tjdut) {
                Some(delta_t_seconds) => Ok((delta_t_seconds + adjust_seconds) / 86400.0),
                None => Err(Error::DeltaTModelUnavailable {
                    model: cfg.model.as_name(),
                }),
            },
            None => Err(Error::DeltaTModelUnavailable {
                model: cfg.model.as_name(),
            }),
        },
    }
}

/// Resolve the effective tidal acceleration for one `deltat_ex` call,
/// using the same physical model flag selection: an explicitly set acceleration
/// (`Some`) is authoritative; otherwise the `FLG_MOSEPH` bit selects
/// `TIDAL_MOSEPH` and every other flag combination uses `TIDAL_DEFAULT`.
/// Only the MOSEPH bit is observed — other flag bits never affect Delta T.
pub fn resolve_tidal_acceleration(flag: i32, explicit: Option<f64>) -> f64 {
    if let Some(acc) = explicit {
        return acc;
    }
    if flag & FLG_MOSEPH != 0 {
        TIDAL_MOSEPH
    } else {
        TIDAL_DEFAULT
    }
}

/// Calculate Delta T in **days** with an explicit ephemeris-selection flag.
///
/// a user override wins before any flag
/// inspection (returned ); otherwise the flag-selected
/// acceleration ([`resolve_tidal_acceleration`]) applies to this
/// computation only and the shared [`deltat_with_config`] chain runs with
/// it. `explicit` carries the globally-set acceleration when one exists
/// (`None` while in automatic mode, letting the flag decide).
pub fn deltat_ex_with_config(
    tjdut: f64,
    flag: i32,
    cfg: &DeltaTConfig<'_>,
    tid_acc_explicit: Option<f64>,
) -> Result<f64> {
    if cfg.user_override_days.is_some() {
        return deltat_with_config(tjdut, cfg);
    }
    let resolved = DeltaTConfig {
        tidal_acceleration: resolve_tidal_acceleration(flag, tid_acc_explicit),
        ..*cfg
    };
    deltat_with_config(tjdut, &resolved)
}

// UTC calendar labels and UT1/TT conversions.
//
// Calendar and clock validation precede leap-label admission. Julian calendar
// labels are mapped to physical Gregorian dates for leap-table lookup. Before
// 1972 and from 2035 the civil label is treated as UT1 and the configured
// Delta-T supplies TT. Within the leap-second era, UTC passes through TAI and
// TT; the inverse reconstruction uses the same leap table and bounded rounding
// snap. Nonfinite epochs and invalid dates/times return typed errors.

/// Julian Day of 1972-01-01 00:00:00 UTC, the start of the leap-second era.
///
/// Instants within [`UTC_EPOCH_BAND_DAYS`]
/// of it take the leap-table UTC reconstruction path even when the
/// UT1 estimate still carries a 1971 label (UT1 leads/lags UTC by up to
/// 0.9 s, so the UT1 year is unreliable there).
pub const UTC_EPOCH_JD: f64 = 2441317.5;

/// Half-width in days of the UTC-epoch guard band (`2 s`, ).
pub const UTC_EPOCH_BAND_DAYS: f64 = 2.0 / 86400.0;

/// Seconds of round-off a calendar reconstruction from `jd` can carry.
///
/// two ulps of the epoch
/// (one JD ulp at modern magnitudes is ~40 us; the reconstruction performs
/// a handful of arithmetic steps on it), in seconds. Derived from the
/// representation itself, not fitted: a fixed larger threshold would flatten
/// valid sub-millisecond inputs to round seconds.
fn reconstruction_snap_tolerance(jd: f64) -> f64 {
    let magnitude = jd.abs();
    2.0 * (magnitude.next_up() - magnitude) * 86400.0
}

/// Split a Julian Day into UTC calendar fields with truncation-toward-zero
/// field splitting, then snap the floating-JD reconstruction
/// floor: a reading within [`reconstruction_snap_tolerance`] of a whole
/// second snaps to it, with the carry cascading through the calendar.
///
/// No leap seconds exist on these
/// paths (pre-1972/post-2035 labels and guard-band fallbacks), so the plain
/// carry applies.
fn jd_to_calendar_tuple(jd: f64, cal: Calendar) -> Result<UtcFields> {
    let date = revjul(jd, cal)?;
    // `revjul` guarantees `hour in [0, 24)`; truncation mirrors `int`.
    let hour = date.hour.trunc() as i32;
    let minute_frac = (date.hour - hour as f64) * 60.0;
    let minute = minute_frac.trunc() as i32;
    let second = (minute_frac - minute as f64) * 60.0;
    let sec_round = second.round_ties_even();
    if second != sec_round && (second - sec_round).abs() < reconstruction_snap_tolerance(jd) {
        let total = hour as i64 * 3600 + minute as i64 * 60 + sec_round as i64;
        if total < 86400 {
            return Ok(UtcFields {
                year: date.year,
                month: date.month,
                day: date.day,
                hour: (total / 3600) as i32,
                minute: ((total % 3600) / 60) as i32,
                second: (total % 60) as f64,
            });
        }
        let jd_next = julday(date.year, date.month, date.day, 12.0, cal)? + 1.0;
        let next = revjul(jd_next, cal)?;
        return Ok(UtcFields {
            year: next.year,
            month: next.month,
            day: next.day,
            hour: 0,
            minute: 0,
            second: 0.0,
        });
    }
    Ok(UtcFields {
        year: date.year,
        month: date.month,
        day: date.day,
        hour,
        minute,
        second,
    })
}

/// Snap a reconstructed Gregorian clock reading to its round instant.
///
/// a reading within
/// `tolerance` seconds of a whole second snaps to it with the carry
/// cascading through the calendar; a genuine leap-second minute keeps its
/// `60.0` reading (probed through [`utc_to_jd_with_config`] with the same
/// table/config), while `60.99998...`-style endings roll into the next day
/// like an ordinary midnight carry. Only the expected "not a real leap
/// second" rejections are swallowed; anything else (non-finite state, an
/// unavailable ΔT backend, integer overflow) propagates as a genuine error.
#[allow(clippy::too_many_arguments)] // positional arity mirrors the model
fn snap_reconstructed_second(
    year: i32,
    month: i32,
    day: i32,
    hour: i32,
    minute: i32,
    second: f64,
    tolerance: f64,
    table: &[LeapSecondEntry],
    dt_cfg: &DeltaTConfig<'_>,
) -> Result<UtcFields> {
    // Use ties-to-even rounding
    // for reconstructed seconds and timezone offsets.
    let sec_round = second.round_ties_even();
    if second == sec_round || (second - sec_round).abs() >= tolerance {
        return Ok(UtcFields {
            year,
            month,
            day,
            hour,
            minute,
            second,
        });
    }
    if sec_round < 60.0 {
        return Ok(UtcFields {
            year,
            month,
            day,
            hour,
            minute,
            second: sec_round,
        });
    }
    if sec_round == 60.0 {
        match utc_to_jd_with_config(
            year,
            month,
            day,
            hour,
            minute,
            60.0,
            Calendar::Gregorian,
            table,
            dt_cfg,
        ) {
            Ok(_) => {
                return Ok(UtcFields {
                    year,
                    month,
                    day,
                    hour,
                    minute,
                    second: 60.0,
                });
            }
            Err(
                Error::InvalidDate { .. } | Error::InvalidTime { .. } | Error::NoLeapSecond { .. },
            ) => {}
            Err(other) => return Err(other),
        }
    }
    let total = hour as i64 * 3600 + minute as i64 * 60 + sec_round as i64;
    if total < 86400 {
        return Ok(UtcFields {
            year,
            month,
            day,
            hour: (total / 3600) as i32,
            minute: ((total % 3600) / 60) as i32,
            second: (total % 60) as f64,
        });
    }
    let jd_next = julday(year, month, day, 12.0, Calendar::Gregorian)? + 1.0;
    let next = revjul(jd_next, Calendar::Gregorian)?;
    Ok(UtcFields {
        year: next.year,
        month: next.month,
        day: next.day,
        hour: 0,
        minute: 0,
        second: 0.0,
    })
}

/// Map a Julian-calendar civil label to its physical Gregorian date at noon.
///
/// The leap-second table and the UTC-era boundaries are Gregorian; the
/// calendar shift relabels the date, never the clock time, so converting at
/// noon yields the right date for any time within the day. Mirrors the
/// the model mapping in `utc_to_jd`.
fn julian_label_to_gregorian(year: i32, month: i32, day: i32) -> Result<(i32, i32, i32)> {
    let noon = julday(year, month, day, 12.0, Calendar::Julian)?;
    let greg = revjul(noon, Calendar::Gregorian)?;
    Ok((greg.year, greg.month, greg.day))
}

/// Convert a UTC date/time to Julian Day numbers `(jd_tt, jd_ut1)`.
///
/// Years use astronomical numbering in `cal`; hours/minutes are integer clock
/// fields and seconds may be fractional. Invalid dates are rejected before
/// invalid times, rather than normalized. Julian labels are mapped to their
/// physical Gregorian dates for era and leap-label decisions.
///
/// For Gregorian years `1972 <= year < 2035`, forms UTC → TAI → TT using
/// `table`, then makes two configured Delta-T substitutions to estimate UT1.
/// Outside that era, treats the civil label directly as UT1 and sets
/// `jd_tt = jd_ut1 + Delta-T`; those branches are not historical/future UTC
/// realizations. UTC before its introduction has no physical definition.
///
/// A `60.x` label is admitted only at 23:59 on a date in
/// [`LEAP_SECOND_DATES`]. The supplied offset table does not extend that list.
/// Admission does not make the ordinary calendar-JD TAI chain preserve a
/// distinct leap label: [`utc_to_tai_jd`] documents that limitation.
/// Scalar JD rounding limits returned time resolution.
///
/// `dt_cfg` supplies the model/override and optional observations; callbacks
/// may acquire their own data. Invalid dates/times, absent leap labels,
/// calendar overflow and selected-model failures retain their typed errors.
/// Use [`hardcoded_leap_seconds`] for a deterministic packaged offset table.
#[allow(clippy::too_many_arguments)] // positional arity mirrors the model
pub fn utc_to_jd_with_config(
    year: i32,
    month: i32,
    day: i32,
    hour: i32,
    minute: i32,
    second: f64,
    cal: Calendar,
    table: &[LeapSecondEntry],
    dt_cfg: &DeltaTConfig<'_>,
) -> Result<(f64, f64)> {
    // Date before time: a bad date wins even when the time is also invalid.
    if !(1..=12).contains(&month) || day < 1 || day > days_in_month(year, month, cal) {
        return Err(Error::InvalidDate { year, month, day });
    }
    if !(0..=23).contains(&hour) || !(0..=59).contains(&minute) || second < 0.0 {
        return Err(Error::InvalidTime {
            hour,
            minute,
            second,
        });
    }
    let (greg_year, greg_month, greg_day) = match cal {
        Calendar::Julian => julian_label_to_gregorian(year, month, day)?,
        Calendar::Gregorian => (year, month, day),
    };
    // A `60.x` label is provisionally admitted; this check decides whether
    // such a 23:59:60 label really exists on the mapped Gregorian date.
    if second >= 60.0 {
        if second >= 61.0 || hour != 23 || minute != 59 || greg_year < 1972 {
            return Err(Error::InvalidTime {
                hour,
                minute,
                second,
            });
        }
        if !is_leap_second_date(greg_year, greg_month, greg_day) {
            return Err(Error::NoLeapSecond {
                hour,
                minute,
                second,
            });
        }
    }
    // Outside the leap-second UTC era the civil label is treated as UT1
    // directly (both branches share the model form).
    if !(1972..2035).contains(&greg_year) {
        let decimal_hour = hour as f64 + minute as f64 / 60.0 + second / 3600.0;
        let jd_ut1 = julday(year, month, day, decimal_hour, cal)?;
        let jd_et = jd_ut1 + deltat_with_config(jd_ut1, dt_cfg)?;
        return Ok((jd_et, jd_ut1));
    }
    // Interim 1972-2035 native chain: UTC label -> TAI -> TT -> UT1.
    let (eff_year, eff_month, eff_day) = match cal {
        Calendar::Julian => (greg_year, greg_month, greg_day),
        Calendar::Gregorian => (year, month, day),
    };
    let jd_tai = utc_to_tai_jd(
        eff_year,
        eff_month,
        eff_day,
        hour,
        minute,
        second,
        Calendar::Gregorian,
        table,
    )?;
    let jd_tt = jd_tai + TT_TAI_OFFSET_DAYS;
    // Fixed-point UT1 refinement in `jdet_to_utc` form: the
    // first ΔT evaluation keys on TT (within ~90 s of the UT1 argument,
    // immaterial to every configured model), the second on the estimate.
    let est = jd_tt - deltat_with_config(jd_tt, dt_cfg)?;
    let jd_ut1 = jd_tt - deltat_with_config(est, dt_cfg)?;
    Ok((jd_tt, jd_ut1))
}

/// Convert a Julian Day in TT (Terrestrial Time, the modern successor to
/// Ephemeris Time) to UTC calendar fields.
///
/// Makes two configured Delta-T substitutions to estimate UT1. Before 1972
/// and from 2035 it returns that UT1 calendar label, not a UTC realization.
/// A two-second guard band at the 1972 boundary instead permits the TAI chain
/// to determine the era. Within the UTC era, subtracts TT−TAI and uses
/// [`tai_jd_to_utc`], then snaps near-integer seconds within two scalar-JD ulps
/// and relabels the date in `cal` if needed.
///
/// Calendar fields use astronomical years. Scalar rounding and the ordinary
/// UTC-JD representation prevent a universal exact or leap-label-preserving
/// round trip. At J2000 TT the ordinary reading is approximately
/// `2000-01-01 11:58:55.816`; microsecond-level reconstruction differences are
/// possible. Nonfinite epochs, calendar overflow and consumed Delta-T errors
/// propagate. See [`utc_to_jd_with_config`] for table and era limitations.
pub fn jdet_to_utc_with_config(
    jd_et: f64,
    cal: Calendar,
    table: &[LeapSecondEntry],
    dt_cfg: &DeltaTConfig<'_>,
) -> Result<UtcFields> {
    if !jd_et.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_et });
    }
    let est_first = jd_et - deltat_with_config(jd_et, dt_cfg)?;
    let jd_ut1_est = jd_et - deltat_with_config(est_first, dt_cfg)?;
    let near_epoch = (jd_ut1_est - UTC_EPOCH_JD).abs() < UTC_EPOCH_BAND_DAYS;
    if revjul(jd_ut1_est, Calendar::Gregorian)?.year < 1972 && !near_epoch {
        return jd_to_calendar_tuple(jd_ut1_est, cal);
    }
    if revjul(jd_ut1_est, Calendar::Gregorian)?.year >= 2035 {
        return jd_to_calendar_tuple(jd_ut1_est, cal);
    }
    // Interim 1972-2035 native chain: TT -> TAI -> UTC fields.
    let jd_tai = jd_et - TT_TAI_OFFSET_DAYS;
    let utc = tai_jd_to_utc(jd_tai, Calendar::Gregorian, table)?;
    // Instants the chain still resolves before 1972 have no UTC: return the
    // UT1 label directly like the model guard-band path.
    if utc.year < 1972 {
        return jd_to_calendar_tuple(jd_ut1_est, cal);
    }
    let snapped = snap_reconstructed_second(
        utc.year,
        utc.month,
        utc.day,
        utc.hour,
        utc.minute,
        utc.second,
        reconstruction_snap_tolerance(jd_ut1_est),
        table,
        dt_cfg,
    )?;
    if cal == Calendar::Julian {
        // Only the calendar date differs between the systems (a whole
        // number of days): remap the date at noon and carry the clock —
        // including a 60.x leap second — over.
        let jd_noon = julday(
            snapped.year,
            snapped.month,
            snapped.day,
            12.0,
            Calendar::Gregorian,
        )?;
        let julian = revjul(jd_noon, Calendar::Julian)?;
        return Ok(UtcFields {
            year: julian.year,
            month: julian.month,
            day: julian.day,
            hour: snapped.hour,
            minute: snapped.minute,
            second: snapped.second,
        });
    }
    Ok(snapped)
}

/// Convert a Julian Day in UT1 to UTC calendar fields.
///
/// Before 1972 and from 2035, returns the input's UT1 civil label directly.
/// Within the UTC era, forms TT with the caller's Delta-T, subtracts TT−TAI,
/// then uses [`tai_jd_to_utc`]. The 1972 guard band, near-integer-second snap
/// and calendar relabeling follow [`jdet_to_utc_with_config`].
///
/// The input is UT1, not UTC. Results depend on the selected Delta-T and offset
/// table, and carry the scalar-rounding and leap-label limitations documented
/// at [`utc_to_jd_with_config`]. No universal exact round trip is promised.
/// Nonfinite epochs, calendar overflow and consumed model errors propagate.
pub fn jdut1_to_utc_with_config(
    jd_ut1: f64,
    cal: Calendar,
    table: &[LeapSecondEntry],
    dt_cfg: &DeltaTConfig<'_>,
) -> Result<UtcFields> {
    if !jd_ut1.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut1 });
    }
    let near_epoch = (jd_ut1 - UTC_EPOCH_JD).abs() < UTC_EPOCH_BAND_DAYS;
    if revjul(jd_ut1, Calendar::Gregorian)?.year < 1972 && !near_epoch {
        return jd_to_calendar_tuple(jd_ut1, cal);
    }
    if revjul(jd_ut1, Calendar::Gregorian)?.year >= 2035 {
        return jd_to_calendar_tuple(jd_ut1, cal);
    }
    // Interim 1972-2035 native chain: UT1 -> TT -> TAI -> UTC fields.
    let jd_tt = jd_ut1 + deltat_with_config(jd_ut1, dt_cfg)?;
    let jd_tai = jd_tt - TT_TAI_OFFSET_DAYS;
    let utc = tai_jd_to_utc(jd_tai, Calendar::Gregorian, table)?;
    if utc.year < 1972 {
        return jd_to_calendar_tuple(jd_ut1, cal);
    }
    let snapped = snap_reconstructed_second(
        utc.year,
        utc.month,
        utc.day,
        utc.hour,
        utc.minute,
        utc.second,
        reconstruction_snap_tolerance(jd_ut1),
        table,
        dt_cfg,
    )?;
    if cal == Calendar::Julian {
        let jd_noon = julday(
            snapped.year,
            snapped.month,
            snapped.day,
            12.0,
            Calendar::Gregorian,
        )?;
        let julian = revjul(jd_noon, Calendar::Julian)?;
        return Ok(UtcFields {
            year: julian.year,
            month: julian.month,
            day: julian.day,
            hour: snapped.hour,
            minute: snapped.minute,
            second: snapped.second,
        });
    }
    Ok(snapped)
}

// TT/TDB conversion using USNO Circular 179 equation 2.6.
//
// The seven-term truncation of the Fairhead-Bretagnon 1990 series returns
// TDB-minus-TT in seconds. The argument is Julian centuries from J2000 TT;
// periodic terms use the published frequencies and phases. The rate helper
// analytically differentiates the same series. Scalar and two-part variants
// share the formula, with two-part subtraction retaining sub-JD-ULP offsets.
// Nonfinite inputs return NonFiniteArgument.

/// Series epoch for the TT/TDB conversion: J2000.0 = JD 2451545.0 TT.
///
/// The argument origin of USNO Circular 179 equation 2.6.
pub const TT_TDB_REFERENCE_JD: f64 = 2451545.0;

/// Julian centuries per series argument unit (the divisor of `T`).
const DAYS_PER_JULIAN_CENTURY: f64 = 36525.0;

/// TDB−TT in seconds at the additive TT Julian date `jd_whole + jd_fraction`.
///
/// Parts need not be integer-day/fractional-day normalized; both must be finite.
/// The J2000 subtraction precedes adding the second part. USNO Circular 179,
/// equation 2.6 supplies a seven-term Fairhead–Bretagnon (1990) truncation.
/// Near-present use of TDB rather than TT as argument changes this approximation
/// by less than a picosecond, not a universal deep-time bound.
/// Nonfinite parts return [`Error::NonFiniteArgument`]; overflowing derived
/// arithmetic is not separately checked and can produce an IEEE nonfinite result.
pub fn tdb_minus_tt_seconds(jd_whole: f64, jd_fraction: f64) -> Result<f64> {
    if !jd_whole.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_whole });
    }
    if !jd_fraction.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_fraction });
    }
    Ok(source_tdb_minus_tt_seconds(jd_whole, jd_fraction))
}

/// TDB-minus-TT seconds of an internally acquired source Time.
///
/// Separate TT day parts retain the seven-term Circular 179 order. Derived
/// nonfinite parts propagate through IEEE arithmetic; the owning state/frame
/// channel determines whether it can serve them. A public caller epoch must
/// already have passed its finite-input guard. This does not relax the direct
/// checked `tdb_minus_tt_seconds` API or certify extreme-epoch accuracy.
pub fn source_tdb_minus_tt_seconds(jd_whole: f64, jd_fraction: f64) -> f64 {
    let t = (jd_whole - TT_TDB_REFERENCE_JD + jd_fraction) / DAYS_PER_JULIAN_CENTURY;
    0.001657 * (628.3076 * t + 6.2401).sin()
        + 0.000022 * (575.3385 * t + 4.2970).sin()
        + 0.000014 * (1256.6152 * t + 6.1969).sin()
        + 0.000005 * (606.9777 * t + 4.0212).sin()
        + 0.000005 * (52.9691 * t + 0.4444).sin()
        + 0.000002 * (21.3299 * t + 5.5431).sin()
        + 0.000010 * t * (628.3076 * t + 4.2490).sin()
}

/// Convert additive TT date parts to a TDB second part, keeping `jd_whole` fixed.
///
/// Adds [`tdb_minus_tt_seconds`] / 86400 to `jd_fraction`; the returned part
/// is in days and is not normalized to `[0, 1)`. It may cross a day boundary.
/// Nonfinite input parts return [`Error::NonFiniteArgument`]; derived overflow
/// is not separately checked. The seven-term clock's model limits apply.
pub fn tt_to_tdb_fraction(jd_whole: f64, jd_fraction: f64) -> Result<f64> {
    Ok(jd_fraction + tdb_minus_tt_seconds(jd_whole, jd_fraction)? / 86400.0)
}

/// Convert a Julian Day in TT to TDB (`TDB = TT + (TDB-TT)`).
///
/// Evaluates the seven-term USNO Circular 179, equation 2.6 series once at TT,
/// without iteration. The near-present correction is about ±1.7 ms;
/// truncation accuracy is approximately 10 µs over 1600–2000 and degrades
/// outside that interval. Scalar JD addition rounds the correction; use
/// [`tt_to_tdb_fraction`] or [`tt_to_tdb_two_part`] to retain a small offset.
/// Nonfinite input returns [`Error::NonFiniteArgument`]; finite extreme epochs
/// can still produce nonfinite formula results. No exact inverse is promised.
pub fn tt_to_tdb_jd(jd_tt: f64) -> Result<f64> {
    Ok(jd_tt + tdb_minus_tt_seconds(jd_tt, 0.0)? / 86400.0)
}

/// Convert a Julian Day in TDB to TT (`TT = TDB - (TDB-TT)`).
///
/// Subtracts the same one-shot series as [`tt_to_tdb_jd`], evaluated at the
/// supplied TDB argument rather than solving an implicit inverse. Its model,
/// rounding and extreme-epoch limitations apply. Nonfinite input returns
/// [`Error::NonFiniteArgument`].
pub fn tdb_to_tt_jd(jd_tdb: f64) -> Result<f64> {
    Ok(jd_tdb - tdb_minus_tt_seconds(jd_tdb, 0.0)? / 86400.0)
}

/// TT-to-TDB conversion for a two-part epoch.
/// The series uses the collapsed argument; adding its correction through the
/// two-part type retains the sub-JD-ULP offset.
pub fn tt_to_tdb_two_part(jd_tt: TwoPartEpoch) -> Result<TwoPartEpoch> {
    let offset_days = tdb_minus_tt_seconds(jd_tt.value(), 0.0)? / 86400.0;
    jd_tt.add_offset(offset_days)
}

/// Convert a two-part Julian Day in TDB to TT.
///
/// Same series as [`tdb_to_tt_jd`], keyed on the collapsed JD, with the
/// offset applied via [`TwoPartEpoch::add_offset`].
pub fn tdb_to_tt_two_part(jd_tdb: TwoPartEpoch) -> Result<TwoPartEpoch> {
    let offset_days = tdb_minus_tt_seconds(jd_tdb.value(), 0.0)? / 86400.0;
    jd_tdb.add_offset(-offset_days)
}

/// Rate `d(TDB)/d(TT)` (dimensionless) at a TT Julian Day: `1` plus the
/// analytic derivative of [`tdb_minus_tt_seconds`] scaled to days per day.
///
/// Differentiates all seven USNO Circular 179, equation 2.6 terms, including
/// the term whose amplitude depends on TT centuries. Near the present the
/// correction to one is of order 1e-10; this is not an all-epoch bound.
/// Multiplying an AU-per-TDB-day velocity by this factor gives AU per TT day
/// for this clock model. It does not establish source or frame coverage.
///
/// Nonfinite TT input returns [`Error::NonFiniteArgument`]. Derived overflow
/// is not separately checked, and finite extrapolation does not extend the
/// seven-term model's physical accuracy interval. See [`tt_to_tdb_jd`].
pub fn tdb_rate(jd_tt: f64) -> Result<f64> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let t = (jd_tt - TT_TDB_REFERENCE_JD) / DAYS_PER_JULIAN_CENTURY;
    Ok(tdb_rate_at_century(t))
}

/// Rate of the acquired source clock, per TT day at the original TT parts.
///
/// The century operand matches `source_tdb_minus_tt_seconds`, including its
/// subtraction/addition order. Do not collapse the epoch before forming it.
/// Source-state acquisition owns finite-input and record-domain validation.
pub(crate) fn source_tdb_rate(jd_whole: f64, jd_fraction: f64) -> f64 {
    let t = (jd_whole - TT_TDB_REFERENCE_JD + jd_fraction) / DAYS_PER_JULIAN_CENTURY;
    tdb_rate_at_century(t)
}

/// Differentiate all seven clock terms, including the secular product.
fn tdb_rate_at_century(t: f64) -> f64 {
    let per_century = 0.001657 * 628.3076 * (628.3076 * t + 6.2401).cos()
        + 0.000022 * 575.3385 * (575.3385 * t + 4.2970).cos()
        + 0.000014 * 1256.6152 * (1256.6152 * t + 6.1969).cos()
        + 0.000005 * 606.9777 * (606.9777 * t + 4.0212).cos()
        + 0.000005 * 52.9691 * (52.9691 * t + 0.4444).cos()
        + 0.000002 * 21.3299 * (21.3299 * t + 5.5431).cos()
        + 0.000010 * ((628.3076 * t + 4.2490).sin() + t * 628.3076 * (628.3076 * t + 4.2490).cos());
    1.0 + per_century / (DAYS_PER_JULIAN_CENTURY * 86400.0)
}

/// Representable TT domain whose split TDB images lie inside a closed interval.
///
/// Bounds are Julian days. The same seven-term clock as [`tt_to_tdb_fraction`]
/// is inverted, then adjacent binary64 TT epochs select the inward endpoints.
/// No sub-ULP clock correction is discarded when checking a source boundary.
/// Nonfinite/inverted bounds, an unrepresentable interval or failure to settle
/// within eight neighboring epochs return a typed error. This transforms a
/// declared source domain; it does not certify the clock's deep-time accuracy.
pub fn tt_interval_from_tdb(low: f64, high: f64) -> Result<(f64, f64)> {
    if !low.is_finite() || !high.is_finite() || low > high {
        return Err(Error::InvalidSegment {
            begin: low,
            end: high,
        });
    }
    let endpoint = |bound: f64, lower: bool| -> Result<f64> {
        let difference = |tt: f64| -> Result<f64> {
            let whole = tt.floor();
            Ok((whole - bound) + tt_to_tdb_fraction(whole, tt - whole)?)
        };
        let mut epoch = tdb_to_tt_jd(bound)?;
        for _ in 0..8 {
            let residual = difference(epoch)?;
            let outside = if lower {
                residual < 0.0
            } else {
                residual > 0.0
            };
            if outside {
                epoch = if lower {
                    epoch.next_up()
                } else {
                    epoch.next_down()
                };
                continue;
            }
            let neighbor = if lower {
                epoch.next_down()
            } else {
                epoch.next_up()
            };
            let adjacent = difference(neighbor)?;
            if if lower {
                adjacent >= 0.0
            } else {
                adjacent <= 0.0
            } {
                epoch = neighbor;
                continue;
            }
            return Ok(epoch);
        }
        Err(Error::ConvergenceFailed {
            algorithm: "tdb_source_boundary_to_tt",
            iterations: 8,
            max_iterations: 8,
            tolerance: 0.0,
            last_value: epoch,
        })
    };
    let result = (endpoint(low, true)?, endpoint(high, false)?);
    if result.0 > result.1 {
        return Err(Error::InvalidSegment {
            begin: low,
            end: high,
        });
    }
    Ok(result)
}

/// Sample an explicit geometric TDB source on the LEB format's TT grid.
/// The callback receives split TDB days and returns ICRS AU/AU per TDB day.
/// The clock's analytic chain rule returns AU per TT day. Input order and
/// sub-JD-ULP clock offsets are retained; any source/clock error aborts the batch.
pub(crate) fn sample_geometric_tt(
    epochs: &[f64],
    mut source: impl FnMut(f64, f64) -> Result<([f64; 3], [f64; 3])>,
) -> Result<Vec<[f64; 6]>> {
    epochs
        .iter()
        .map(|&tt| {
            let whole = tt.floor();
            let fraction = tt_to_tdb_fraction(whole, tt - whole)?;
            let (p, v) = source(whole, fraction)?;
            let rate = tdb_rate(tt)?;
            Ok([p[0], p[1], p[2], v[0] * rate, v[1] * rate, v[2] * rate])
        })
        .collect()
}

// Sidereal time and the equation of the equinoxes.
//
// Inputs are UT1 Julian days and explicit TT-minus-UT1/nutation values. Mean
// sidereal time uses IAU 2006 inside 1850–2050: ERA plus the precession-in-RA
// polynomial of Capitaine, Wallace and Chapront 2003. Outside that interval,
// Simon 1994 Earth longitude and Vondrák 2011 precession provide the long-term
// geometric extension, joined through continuity offsets. Apparent time adds
// dpsi*cos(true obliquity). Outputs are degrees or hours as each function
// specifies; nonfinite epochs and required provider failures return errors.

use crate::positions::{
    ecliptic_pole_vondrak, equator_pole_vondrak, vondrak_mean_obliquity_deg,
    vondrak_mean_obliquity_rad,
};

/// Scalar-JD long-term mean-equatorial precession for sidereal time.
///
/// `jd_tt` is TT Julian Day; return row-major J2000-mean to of-date
/// mean-equatorial rotation acting on a Cartesian column vector. Evaluate
/// the independent scalar-JD equator/ecliptic poles, normalize their cross
/// product by component-wise division, then form the middle row as equator
/// cross equinox. This route deliberately differs from Julian-epoch ERFA
/// precession: acquisition and separate-product accumulation matter to values.
/// The Vondrák 2011/2012 model spans ±200,000 years about J2000 with
/// growing physical uncertainty; numerical identity is not an accuracy
/// certificate. Non-finite input, undefined equinox and arithmetic overflow
/// are refused.
pub(crate) fn sidereal_precession_matrix(jd_tt: f64) -> Result<[[f64; 3]; 3]> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let equator = equator_pole_vondrak(jd_tt);
    let ecliptic = ecliptic_pole_vondrak(jd_tt);
    let cross = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let v = cross(equator, ecliptic);
    let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if !norm.is_finite() {
        return Err(Error::MathDomain { value: norm });
    }
    if norm == 0.0 {
        return Err(Error::SingularMatrix);
    }
    let equinox = v.map(|value| value / norm);
    let matrix = [equinox, cross(equator, equinox), equator];
    for row in matrix {
        if let Some(&value) = row.iter().find(|value| !value.is_finite()) {
            return Err(Error::MathDomain { value });
        }
    }
    Ok(matrix)
}

/// TT Julian date corresponding to UT1 under the bundled default timescale.
///
/// Uses the merged observed Delta-T table and long-term spline already
/// supplied by the default SMH-2016 provider, retaining the timescale's
/// UT1-to-TT refinement and split-day operation sequence before producing a
/// scalar Julian date. No leap-second rounding, network or ephemeris asset
/// access is involved.
///
/// This acquisition is independent of caller-selected Delta-T models,
/// overrides and tidal acceleration. Those settings still apply to apparent
/// body calculations and sidereal-time model terms; automatically acquired
/// nutation/obliquity uses this timescale instead. The bundled observations
/// are a fixed dataset, not live DUT1/IERS updates or a guarantee of future
/// physical accuracy. Non-finite inputs or unrepresentable results return
/// typed errors rather than a fabricated epoch.
pub fn default_timescale_tt_from_ut1(jd_ut: f64) -> Result<f64> {
    crate::smh2016::timescale_tt_from_ut1(jd_ut)
}

/// Acquire a source-style scalar UT1 Time without collapsing its TT parts.
///
/// Returns `(whole, tt_fraction, original_ut1_fraction)` in Julian days. The
/// whole part is the entire supplied UT1 value, TT fraction comes from exactly
/// two bundled SMH-2016 evaluations, and UT1 fraction is `0.0 * jd_ut1`.
/// A fresh two-argument TT Time can therefore retain the exact acquired parts;
/// model Delta-T, user offsets, IERS configuration and a scalar TT inverse do
/// not enter this acquisition. No physical or convergence bound is implied.
///
/// IEEE nonfinite inputs/intermediates propagate, including the original
/// fraction's NaN and signed zero. This is an already-admitted Source primitive:
/// the enclosing public entry owns its actual input/error stage. Checked
/// [`default_timescale_tt_from_ut1`] retains its original refusal contract.
pub fn source_timescale_ut1_parts(jd_ut1: f64) -> (f64, f64, f64) {
    crate::smh2016::source_time_parts_from_ut1(jd_ut1)
}

/// Return bundled default TT and the directly evaluated TT-minus-UT1 value.
///
/// `jd_ut1` is a finite scalar UT1 Julian day. If `D(t)` is the bundled
/// historical/IERS and long-term TT-minus-UT1 curve in seconds, evaluate
/// `t1 = jd_ut1 + D(jd_ut1)/86400`, then
/// `tt = jd_ut1 + D(t1)/86400`. Return `(tt, D(tt))`: the second component
/// is a direct model observable in seconds, not the offset used to construct
/// TT and not a subtraction of rounded scalar Julian days. Exactly two
/// substitutions construct TT; there is no iterative convergence criterion.
///
/// Caller model, tidal acceleration and user overrides do not enter this
/// fixed timescale. The bundled table is fixed data, with extrapolation rather
/// than a guarantee of physical accuracy at remote epochs. Non-finite input
/// and nonrepresentable model arithmetic retain the existing typed time
/// errors. The two returned values are fully computed before caller delivery.
pub fn default_timescale_from_ut1(jd_ut1: f64) -> Result<(f64, f64)> {
    let tt = default_timescale_tt_from_ut1(jd_ut1)?;
    let delta = crate::smh2016::smh2016_lookup_seconds(jd_ut1)
        .ok_or(Error::NonFiniteArgument { value: jd_ut1 })?;
    if !delta.is_finite() {
        return Err(Error::MathDomain { value: jd_ut1 });
    }
    Ok((tt, delta))
}

/// Convert TT Julian Day to UT1 using the fixed bundled enhanced timescale.
///
/// Evaluate the bundled historical/IERS and long-term TT-minus-UT1 curve
/// once at TT, then subtract its days offset from the fractional day before
/// adding the whole day. This preserves the scalar inverse observable;
/// it is not a conversion through the caller's model/tidal/user overrides.
/// Non-finite inputs or arithmetic return typed errors. Custom tables,
/// live IERS and UTC leap labels are outside this fixed-data path; no
/// physical-accuracy or exact round-trip guarantee is asserted.
/// See [`default_timescale_tt_from_ut1`] for the forward timescale contract.
pub fn default_timescale_ut1_from_tt(jd_tt: f64) -> Result<f64> {
    crate::smh2016::timescale_ut1_from_tt(jd_tt)
}

/// J2000.0 epoch, in days.
pub const SIDEREAL_J2000_JD: f64 = 2451545.0;
/// Start of the modern IAU-2006 GMST window: 1 Jan 1850.
pub const SIDEREAL_MODERN_START_JD: f64 = 2396758.5;
/// End of the modern IAU-2006 GMST window: 1 Jan 2050.
pub const SIDEREAL_MODERN_END_JD: f64 = 2469807.5;
/// Sun-Earth light-time in days (`AU / c`, ).
pub const SUN_EARTH_LIGHT_TIME_DAYS: f64 = 1.495978707e11 / 299792458.0 / 86400.0;
/// Days per Julian millennium, the long-term-method time unit.
pub const DAYS_PER_JULIAN_MILLENNIUM: f64 = 365250.0;

/// Explicit TT-minus-UT1 and legacy sidereal boundary corrections.
///
/// Delta-T is in days; the two offsets are in degrees. Canonical mean time
/// consumes only Delta-T to form TT, independently of UT1 Earth rotation.
/// The offsets retain their actual values and effects in the deprecated
/// [`legacy_sidereal`] utilities; canonical routes never apply them.
#[derive(Debug, Clone, Copy)]
pub struct SiderealDeltaT {
    delta_t_days: f64,
    offset_lo_deg: f64,
    offset_hi_deg: f64,
}

impl SiderealDeltaT {
    /// Checked constructor: every value must be finite.
    pub fn try_new(delta_t_days: f64, offset_lo_deg: f64, offset_hi_deg: f64) -> Result<Self> {
        for v in [delta_t_days, offset_lo_deg, offset_hi_deg] {
            if !v.is_finite() {
                return Err(Error::NonFiniteArgument { value: v });
            }
        }
        Ok(Self {
            delta_t_days,
            offset_lo_deg,
            offset_hi_deg,
        })
    }

    /// TT-minus-UT1 at the epoch, in days.
    pub fn delta_t_days(&self) -> f64 {
        self.delta_t_days
    }

    /// Continuity offset at the 1850 boundary, in degrees.
    pub fn offset_lo_deg(&self) -> f64 {
        self.offset_lo_deg
    }

    /// Continuity offset at the 2050 boundary, in degrees.
    pub fn offset_hi_deg(&self) -> f64 {
        self.offset_hi_deg
    }
}

/// Greenwich Mean Sidereal Time in degrees via the IAU 2006 expression
///
/// `ERA(UT1)` uses the Capitaine, Guinot & McCarthy (2000) constants and
/// `p(T)` is the Capitaine, Wallace & Chapront (2003) precession-in-RA
/// polynomial. The result is not range-reduced here. Non-finite inputs are
/// [`Error::NonFiniteArgument`].
pub fn gmst_iau2006_deg(jd_ut1: f64, jd_tt: f64) -> Result<f64> {
    if !jd_ut1.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut1 });
    }
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let du = jd_ut1 - SIDEREAL_J2000_JD;
    // association (precedent: the sub-unity-rate ERA split
    // lives in `positions::earth_rotation_angle_rad`; this is the model
    // full-rate spelling used only inside the IAU 2006 polynomial branch).
    // the model spelling (`1.00273781191135448`): the trailing
    // digits are beyond f64 precision but the literal parses to the same
    // value as the truncated form in both languages (correctly rounded
    // decimal parsing; value verified `0x3ff00b36cdc9f32b`).
    #[allow(clippy::excessive_precision)]
    let era = 2.0 * std::f64::consts::PI * (0.7790572732640 + 1.00273781191135448 * du);
    let t = (jd_tt - SIDEREAL_J2000_JD) / DAYS_PER_JULIAN_CENTURY;
    let p = 0.014506
        + t * (4612.156534
            + t * (1.3915817 + t * (-0.00000044 + t * (-0.000029956 + t * -0.0000000368))));
    Ok(era.to_degrees() + p / 3600.0)
}

/// Rotate a 3-vector about the x-axis by `eps` radians, the model
/// `sidereal_longterm._rot_x` spelling (ecliptic-to-equator uses `-eps`).
///
/// Kept as a separate function (rather than reusing the position-pipeline
/// active rotations) because the component association is the model's own
/// and the sweep pins it exactly.
fn rot_x_release(v: [f64; 3], eps: f64) -> [f64; 3] {
    let c = eps.cos();
    let s = eps.sin();
    [v[0], v[1] * c + v[2] * s, -v[1] * s + v[2] * c]
}

/// Greenwich mean sidereal time in degrees via the long-term geometric
/// method.
///
/// `delta_t_days` is TT-minus-UT1 at `jd_ut1`. The result is range-reduced to `[0, 360)` with the floored
/// modulo. Non-finite inputs are [`Error::NonFiniteArgument`].
pub fn mean_sidereal_longterm_deg(jd_ut1: f64, delta_t_days: f64) -> Result<f64> {
    if !jd_ut1.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut1 });
    }
    if !delta_t_days.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: delta_t_days,
        });
    }
    let jd_tt = jd_ut1 + delta_t_days;
    let t = (jd_tt - SIDEREAL_J2000_JD) / DAYS_PER_JULIAN_MILLENNIUM;
    let t2 = t * t;
    let t3 = t * t2;
    // Simon et al. (1994) mean longitude of the Earth, association.
    let mut dlon = 100.46645683 + (1295977422.83429 * t - 2.04411 * t2 - 0.00523 * t3) / 3600.0;
    // Sun-Earth light-time correction, association.
    dlon = (dlon - SUN_EARTH_LIGHT_TIME_DAYS * 360.0 / 365.2425).rem_euclid(360.0);
    let rad = dlon.to_radians();
    let mut vec = [rad.cos(), rad.sin(), 0.0];
    // `_J2000 is already the J2000.0 TT epoch... adding Delta-T would
    // double-count the UT->TT offset`.
    let eps_j2000 = vondrak_mean_obliquity_rad(SIDEREAL_J2000_JD)?;
    vec = rot_x_release(vec, -eps_j2000);
    // `_precess_j2000_to_date`: row-major multiply, association.
    let r = sidereal_precession_matrix(jd_tt)?;
    vec = [
        r[0][0] * vec[0] + r[0][1] * vec[1] + r[0][2] * vec[2],
        r[1][0] * vec[0] + r[1][1] * vec[1] + r[1][2] * vec[2],
        r[2][0] * vec[0] + r[2][1] * vec[1] + r[2][2] * vec[2],
    ];
    let eps_date = vondrak_mean_obliquity_rad(jd_tt)?;
    vec = rot_x_release(vec, eps_date);
    let lon = vec[1].atan2(vec[0]).to_degrees();
    // `math.fmod` is the truncated remainder: native `%`.
    let dhour = (jd_ut1 - 0.5) % 1.0 * 360.0;
    Ok((lon + dhour).rem_euclid(360.0))
}

/// Continuity offset in degrees joining the long-term branch onto the
/// modern polynomial branch at `boundary_jd`.
///
/// `delta_t_days` is TT-minus-UT1 at the boundary. Computed from the live model rather than hardcoded, so it
/// can never drift out of sync when the Delta-T model changes. Non-finite inputs are [`Error::NonFiniteArgument`].
pub fn sidereal_continuity_offset_deg(boundary_jd: f64, delta_t_days: f64) -> Result<f64> {
    if !boundary_jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: boundary_jd });
    }
    if !delta_t_days.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: delta_t_days,
        });
    }
    let jd_tt = boundary_jd + delta_t_days;
    let diff = mean_sidereal_longterm_deg(boundary_jd, delta_t_days)?
        - gmst_iau2006_deg(boundary_jd, jd_tt)?.rem_euclid(360.0);
    Ok((diff + 180.0).rem_euclid(360.0) - 180.0)
}

/// Both window-boundary continuity offsets from the boundary Delta-T values
pub fn sidereal_continuity_offsets(
    delta_t_lo_days: f64,
    delta_t_hi_days: f64,
) -> Result<(f64, f64)> {
    Ok((
        sidereal_continuity_offset_deg(SIDEREAL_MODERN_START_JD, delta_t_lo_days)?,
        sidereal_continuity_offset_deg(SIDEREAL_MODERN_END_JD, delta_t_hi_days)?,
    ))
}

/// Canonical Greenwich mean sidereal time, in degrees `[0,360)`.
///
/// UT1 drives compensated ERA; TT=`jd_ut1+dt.delta_t_days()` drives the
/// continuous Vondrák mean-pole CIO transport. GMST is ERA minus relative
/// mean EO, plus the published IAU2006 datum0.014506arcsec and P03Eq42
/// secular term−0.003872arcsec per TT Julian century. A fixed-node locator
/// table and GL16 correction evaluate the transport without heap allocation.
/// No1850/2050branch or legacy boundary offset enters this calculation.
///
/// Numerical/model qualification covers±2000centuries around J2000; finite
/// dates beyond remain admitted but model accuracy is unqualified there.
/// Nonfinite UT1 returns [`Error::NonFiniteArgument`]. Undefined arithmetic
/// otherwise retains IEEE propagation. See Vondrák et al.(2011, A&A534A22)
/// and Capitaine et al.(2003, A&A412,567); this is the selected VCW convention,
/// while [`sidereal_time_iau2006_hours`] is the explicit standard utility.
pub fn mean_sidereal_time_deg(jd_ut1: f64, dt: &SiderealDeltaT) -> Result<f64> {
    canonical_sidereal::mean_degrees(jd_ut1, jd_ut1 + dt.delta_t_days())
}

/// Equation of the equinoxes in degrees: `dpsi * cos(eps_true)`.
///
/// This is the apparent-minus-mean correction the model applies inline in
/// both `apparent_sidereal_time_deg` and `_sidtime_internal`; naming it
/// lets the two paths share one spelling. Non-finite inputs are
/// [`Error::NonFiniteArgument`].
pub fn equation_of_equinoxes_deg(dpsi_deg: f64, eps_true_deg: f64) -> Result<f64> {
    if !dpsi_deg.is_finite() {
        return Err(Error::NonFiniteArgument { value: dpsi_deg });
    }
    if !eps_true_deg.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: eps_true_deg,
        });
    }
    Ok(dpsi_deg * eps_true_deg.to_radians().cos())
}

/// Local Apparent Sidereal Time at `jd_ut1` in degrees `[0, 360)`
///
/// This is the ARMC the house systems use when `longitude_deg` is the
/// observer's geographic longitude (East positive); longitude `0` gives
/// Greenwich apparent sidereal time. Nutation is caller-supplied so the
/// tuned nutation series is reused unchanged.
/// Non-finite inputs are [`Error::NonFiniteArgument`].
pub fn apparent_sidereal_time_deg(
    jd_ut1: f64,
    longitude_deg: f64,
    dt: &SiderealDeltaT,
    dpsi_deg: f64,
    eps_true_deg: f64,
) -> Result<f64> {
    if !jd_ut1.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut1 });
    }
    if !longitude_deg.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: longitude_deg,
        });
    }
    let gmst = mean_sidereal_time_deg(jd_ut1, dt)?;
    let eoe = equation_of_equinoxes_deg(dpsi_deg, eps_true_deg)?;
    Ok((gmst + eoe + longitude_deg).rem_euclid(360.0))
}

/// Local Apparent Sidereal Time in hours `[0, 24)`.
///
/// `jd_ut` is the Julian day in UT1, `longitude_deg` the geographic
/// longitude (East positive), `obliquity_deg` the (true) obliquity of the
/// ecliptic and `nutation_deg` the nutation in longitude. With explicit
/// values this is the full public computation; the `None` arms of the
/// `sidtime` (obliquity/nutation from `calc_ut ... ECL_NUT`) stay
/// with the serving layer. Non-finite inputs are
/// [`Error::NonFiniteArgument`].
///
/// Worked value (longitude `0`): at JD `2451545.0` this function
/// returns `18.697374828694365` h (approximately `18.6974` h) under the
/// canonical VCW convention. The legacy boundary offsets are not applied.
///
/// ```rust
/// use libephemeris_core::time::{SiderealDeltaT, greenwich_apparent_sidereal_time_hours};
///
/// // Pinned Delta-T state at J2000 (same values as the `sidereal_test_words`
/// // helper: Espenak-Meeus Delta-T with the 1850/2050 continuity offsets).
/// let dt = SiderealDeltaT::try_new(
///     0.0007392804723490653,
///     0.0003809980328242091,
///     0.0010661391453652413,
/// )
/// .unwrap();
/// let gast =
///     greenwich_apparent_sidereal_time_hours(2451545.0, 23.4392911, 0.0, &dt)
///         .unwrap();
/// assert_eq!(gast, 18.697374828694365);
/// ```
pub fn local_apparent_sidereal_time_hours(
    jd_ut: f64,
    longitude_deg: f64,
    obliquity_deg: f64,
    nutation_deg: f64,
    dt: &SiderealDeltaT,
) -> Result<f64> {
    if !jd_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut });
    }
    if !longitude_deg.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: longitude_deg,
        });
    }
    if !obliquity_deg.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: obliquity_deg,
        });
    }
    if !nutation_deg.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: nutation_deg,
        });
    }
    local_apparent_sidereal_time_hours_from_angles(
        jd_ut,
        longitude_deg,
        obliquity_deg,
        nutation_deg,
        dt,
    )
}

/// Combine the checked UT1/model clock with previously evaluated frame angles.
///
/// Epoch and east-positive longitude retain the public finite-input domains.
/// Obliquity and longitude nutation are degrees and retain their IEEE values:
/// NaN propagates through the equation of equinoxes and the returned local
/// apparent sidereal time in hours. This supports a caller which has already
/// consumed its own frame operations and must preserve their non-finite and
/// error chronology. The public checked function validates all four inputs.
pub(crate) fn local_apparent_sidereal_time_hours_from_angles(
    jd_ut: f64,
    longitude_deg: f64,
    obliquity_deg: f64,
    nutation_deg: f64,
    dt: &SiderealDeltaT,
) -> Result<f64> {
    if !jd_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut });
    }
    if !longitude_deg.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: longitude_deg,
        });
    }
    let gast =
        greenwich_apparent_sidereal_time_hours_from_angles(jd_ut, obliquity_deg, nutation_deg, dt)?;
    let last = gast + longitude_deg / 15.0;
    // `% 24` plus the (unreachable for finite inputs) negative
    // guard the model keeps with its `pragma: no cover` note.
    let mut last = last % 24.0;
    if last < 0.0 {
        last += 24.0;
    }
    Ok(last)
}

/// Unnormalized Greenwich rotation with already evaluated frame angles.
///
/// `jd_ut` is finite UT1 Julian days; the checked model `dt` supplies Delta-T
/// and continuity offsets. Obliquity and longitude nutation are degrees and
/// preserve their IEEE values. The mean sidereal angle retains its radians
/// round-trip before conversion to hours. The equation of equinoxes is then
/// added in hours, without longitude, modulo or further frame acquisition.
/// This closed-form transport reports the mean-time primitive's errors and
/// imposes no additional frame-angle domain checks.
pub(crate) fn greenwich_apparent_sidereal_time_hours_from_angles(
    jd_ut: f64,
    obliquity_deg: f64,
    nutation_deg: f64,
    dt: &SiderealDeltaT,
) -> Result<f64> {
    if !jd_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut });
    }
    let gmst_hours = mean_sidereal_time_deg(jd_ut, dt)?.to_radians().to_degrees() / 15.0;
    let equation_of_equinoxes = nutation_deg * obliquity_deg.to_radians().cos() / 15.0;
    Ok(gmst_hours + equation_of_equinoxes)
}

/// Greenwich Apparent Sidereal Time in hours `[0, 24)`.
///
/// Equivalent to [`local_apparent_sidereal_time_hours`] at longitude `0`.
/// Non-finite inputs are [`Error::NonFiniteArgument`].
pub fn greenwich_apparent_sidereal_time_hours(
    jd_ut: f64,
    obliquity_deg: f64,
    nutation_deg: f64,
    dt: &SiderealDeltaT,
) -> Result<f64> {
    local_apparent_sidereal_time_hours(jd_ut, 0.0, obliquity_deg, nutation_deg, dt)
}

/// Sidereal time at longitude 0 (= ARMC) and true obliquity for houses
///
/// `tjdut` is the Julian Day in UT1; `dpsi_deg`/`deps_deg` are the nutation
/// pair in degrees. Returns `(armc0_deg, eps_true_deg)`. The LEB-backed nutation
/// acquisition stays with the house serving layer. Non-finite
/// inputs are [`Error::NonFiniteArgument`].
pub fn house_armc_obliquity_deg(
    tjdut: f64,
    dt: &SiderealDeltaT,
    dpsi_deg: f64,
    deps_deg: f64,
) -> Result<(f64, f64)> {
    if !tjdut.is_finite() {
        return Err(Error::NonFiniteArgument { value: tjdut });
    }
    if !dpsi_deg.is_finite() {
        return Err(Error::NonFiniteArgument { value: dpsi_deg });
    }
    if !deps_deg.is_finite() {
        return Err(Error::NonFiniteArgument { value: deps_deg });
    }
    let jd_tt = tjdut + dt.delta_t_days();
    let eps_mean = vondrak_mean_obliquity_deg(jd_tt)?;
    let eps_true = eps_mean + deps_deg;
    let armc0 = apparent_sidereal_time_deg(tjdut, 0.0, dt, dpsi_deg, eps_true)?;
    Ok((armc0, eps_true))
}

pub(crate) mod public_sidereal;

#[cfg(test)]
mod tests {
    use super::*;

    const GREG: Calendar = Calendar::Gregorian;
    const JUL: Calendar = Calendar::Julian;

    #[test]
    fn calendar_flag_round_trip() {
        assert_eq!(Calendar::from_flag(0), Ok(Calendar::Julian));
        assert_eq!(Calendar::from_flag(1), Ok(Calendar::Gregorian));
        assert_eq!(Calendar::Julian.flag(), JUL_CAL);
        assert_eq!(Calendar::Gregorian.flag(), GREG_CAL);
        assert_eq!(
            Calendar::from_flag(99),
            Err(Error::InvalidCalendar { value: 99 })
        );
    }

    #[test]
    fn j2000_epoch_identity() {
        // JD 2451545.0 == 2000-01-01 12:00.
        assert_eq!(julday(2000, 1, 1, 12.0, GREG), Ok(2451545.0));
        assert_eq!(
            revjul(2451545.0, GREG),
            Ok(CalendarDate {
                year: 2000,
                month: 1,
                day: 1,
                hour: 12.0,
            })
        );
        assert_eq!(julday(2000, 1, 1, 0.0, GREG), Ok(2451544.5));
    }

    #[test]
    fn gregorian_reform_boundaries() {
        // Pinned against the model implementation (scratch oracle @ 37 ffb 05):
        // Oct 15 1582 Gregorian noon == JD 2299161 (JD_GREGORIAN_REFORM);
        // Oct 4 1582 Julian noon is one day earlier; the proleptic
        // Gregorian reading of Oct 4 trails by the 10-day gap.
        assert_eq!(julday(1582, 10, 15, 12.0, GREG), Ok(2299161.0));
        assert_eq!(julday(1582, 10, 4, 12.0, JUL), Ok(2299160.0));
        assert_eq!(julday(1582, 10, 4, 12.0, GREG), Ok(2299150.0));
        assert_eq!(julday(1582, 10, 15, 12.0, JUL), Ok(2299171.0));
        assert_eq!(
            revjul(2299160.5, GREG),
            Ok(CalendarDate {
                year: 1582,
                month: 10,
                day: 15,
                hour: 0.0,
            })
        );
        assert_eq!(
            revjul(2299160.5, JUL),
            Ok(CalendarDate {
                year: 1582,
                month: 10,
                day: 5,
                hour: 0.0,
            })
        );
    }

    #[test]
    fn julian_day_zero_and_year_zero() {
        // JD 0.0 == -4712-01-01 12:00 Julian (4713 BCE, astronomical year).
        assert_eq!(julday(-4712, 1, 1, 12.0, JUL), Ok(0.0));
        assert_eq!(
            revjul(0.0, JUL),
            Ok(CalendarDate {
                year: -4712,
                month: 1,
                day: 1,
                hour: 12.0,
            })
        );
        // Year 0 == 1 BCE; Gregorian/Julian proleptic readings differ by 2 d.
        assert_eq!(julday(0, 1, 1, 12.0, GREG), Ok(1721060.0));
        assert_eq!(julday(0, 1, 1, 12.0, JUL), Ok(1721058.0));
    }

    #[test]
    fn floor_not_truncation_for_negative_years() {
        // After the Jan shift, year -102 gives a = floor(-1.02) = -2;
        // truncation toward zero (-1) would move these JDs by a day.
        assert_eq!(julday(-101, 1, 15, 6.0, GREG), Ok(1684184.75));
        assert_eq!(julday(-101, 1, 15, 6.0, JUL), Ok(1684181.75));
        assert_eq!(julday(-44, 3, 15, 12.0, JUL), Ok(1705061.0));
        assert_eq!(julday(-100, 1, 1, 12.0, GREG), Ok(1684536.0));
    }

    #[test]
    fn leap_century_and_calendar_divergence() {
        // 1900 is not a Gregorian leap year (no Feb 29): consecutive days.
        assert_eq!(julday(1900, 2, 28, 0.0, GREG), Ok(2415078.5));
        assert_eq!(julday(1900, 3, 1, 0.0, GREG), Ok(2415079.5));
        // 2000 and 2400 are Gregorian leap years.
        assert_eq!(julday(2000, 2, 29, 12.0, GREG), Ok(2451604.0));
        assert_eq!(julday(2400, 2, 29, 12.0, GREG), Ok(2597701.0));
        // Same civil date drifts apart: 16-day gap by 2400.
        assert_eq!(julday(2400, 2, 29, 12.0, JUL), Ok(2597717.0));
        assert_eq!(julday(2100, 3, 1, 12.0, GREG), Ok(2488129.0));
        assert_eq!(julday(2100, 3, 1, 12.0, JUL), Ok(2488143.0));
        // Ordinary modern leap day.
        assert_eq!(julday(2024, 2, 29, 0.0, GREG), Ok(2460369.5));
        assert_eq!(julday(333, 1, 27, 12.0, GREG), Ok(1842712.0));
    }

    #[test]
    fn out_of_range_fields_flow_through_without_validation() {
        // owns date validation; julday computes through, exactly as
        // the source does (each pinned against the model oracle).
        assert_eq!(
            julday(2024, 13, 1, 12.0, GREG),
            julday(2025, 1, 1, 12.0, GREG)
        );
        assert_eq!(
            julday(2024, 0, 1, 12.0, GREG),
            julday(2023, 12, 1, 12.0, GREG)
        );
        assert_eq!(
            julday(2024, 3, 0, 12.0, GREG),
            julday(2024, 2, 29, 12.0, GREG)
        );
    }

    #[test]
    fn round_trip_is_exact_for_representable_hours() {
        // Hours 0/6/12/18 divide 24 exactly in binary64, so the day-fraction
        // round-trips bit-exactly across eras and both calendars.
        let dates = [
            (2000, 1, 1),
            (2024, 2, 29),
            (1900, 2, 28),
            (1582, 10, 15),
            (1582, 10, 4),
            (1, 1, 1),
            (0, 6, 30),
            (-100, 1, 1),
            (-101, 1, 15),
            (-4712, 1, 2),
            (5000, 7, 20),
        ];
        for &(year, month, day) in &dates {
            for &hour in &[0.0, 6.0, 12.0, 18.0] {
                for &cal in &[GREG, JUL] {
                    let jd = julday(year, month, day, hour, cal).unwrap();
                    assert!(jd.is_finite());
                    assert_eq!(
                        revjul(jd, cal),
                        Ok(CalendarDate {
                            year,
                            month,
                            day,
                            hour,
                        }),
                        "round trip failed for ({year}-{month}-{day} {hour} {cal:?})"
                    );
                }
            }
        }
    }

    #[test]
    fn non_finite_inputs_are_errors() {
        assert!(matches!(
            julday(2000, 1, 1, f64::NAN, GREG),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert_eq!(
            julday(2000, 1, 1, f64::INFINITY, GREG),
            Err(Error::NonFiniteArgument {
                value: f64::INFINITY,
            })
        );
        assert!(matches!(
            revjul(f64::NAN, GREG),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert_eq!(
            revjul(f64::INFINITY, JUL),
            Err(Error::NonFiniteArgument {
                value: f64::INFINITY,
            })
        );
    }

    #[test]
    fn days_in_month_follows_calendar_leap_rules() {
        // Gregorian century rules vs Julian every-4 th-year.
        assert_eq!(days_in_month(2000, 2, GREG), 29);
        assert_eq!(days_in_month(1900, 2, GREG), 28);
        assert_eq!(days_in_month(1900, 2, JUL), 29);
        assert_eq!(days_in_month(2000, 2, JUL), 29);
        assert_eq!(days_in_month(2024, 2, GREG), 29);
        assert_eq!(days_in_month(2023, 2, GREG), 28);
        assert_eq!(days_in_month(2024, 4, GREG), 30);
        assert_eq!(days_in_month(2024, 1, GREG), 31);
        // Astronomical year numbering: year 0 is divisible by 400.
        assert_eq!(days_in_month(0, 2, GREG), 29);
        assert_eq!(days_in_month(-100, 2, GREG), 28);
        assert_eq!(days_in_month(-100, 2, JUL), 29);
        // Out-of-range months flow through as 31 (validation is
        // `convert_date_checked`'s job, mirroring the source).
        assert_eq!(days_in_month(2024, 13, GREG), 31);
        assert_eq!(days_in_month(2024, 0, JUL), 31);
    }

    #[test]
    fn leap_second_table_matches_release_census() {
        assert_eq!(LEAP_SECOND_DATES.len(), 27);
        assert!(LEAP_SECOND_DATES.windows(2).all(|w| w[0] < w[1]));
        assert!(is_leap_second_date(1972, 6, 30));
        assert!(is_leap_second_date(2016, 12, 31));
        assert!(!is_leap_second_date(2017, 1, 1));
        assert!(!is_leap_second_date(1972, 1, 1));
        assert!(!is_leap_second_date(2024, 12, 31));
    }

    #[test]
    fn calendar_target_codes_are_case_insensitive() {
        assert_eq!(CalendarTarget::from_code('j'), Ok(CalendarTarget::Julian));
        assert_eq!(
            CalendarTarget::from_code('G'),
            Ok(CalendarTarget::Gregorian)
        );
        assert_eq!(CalendarTarget::Julian.code(), 'j');
        assert_eq!(CalendarTarget::Gregorian.calendar(), GREG);
        assert_eq!(
            CalendarTarget::from_code('x'),
            Err(Error::InvalidCalendarCode { code: 'x' })
        );
        // Observable message preserves the model ValueError wording.
        assert_eq!(
            format!("{}", Error::InvalidCalendarCode { code: 'x' }),
            "calendar must be 'j' or 'g', got: 'x'"
        );
    }

    #[test]
    fn date_conversion_crosses_the_reform_gap() {
        // Pinned against the model oracle @ 37 ffb 05 (ephemeral scratch).
        let jul = CalendarTarget::Julian;
        let greg = CalendarTarget::Gregorian;
        // First Gregorian date <-> Julian: the 10-day gap.
        assert_eq!(
            date_conversion(1582, 10, 15, 12.0, jul),
            Ok(CalendarDate {
                year: 1582,
                month: 10,
                day: 5,
                hour: 12.0,
            })
        );
        assert_eq!(
            date_conversion(1582, 10, 5, 12.0, greg),
            Ok(CalendarDate {
                year: 1582,
                month: 10,
                day: 15,
                hour: 12.0,
            })
        );
        // Julian Oct 4 lands on Gregorian Oct 14 (not Oct 15).
        assert_eq!(
            date_conversion(1582, 10, 4, 12.0, greg),
            Ok(CalendarDate {
                year: 1582,
                month: 10,
                day: 14,
                hour: 12.0,
            })
        );
        // Same-calendar input returns the original fields unchanged.
        assert_eq!(
            date_conversion(2024, 2, 29, 6.5, greg),
            Ok(CalendarDate {
                year: 2024,
                month: 2,
                day: 29,
                hour: 6.5,
            })
        );
        // Classification uses the calendar DATE, not the hour-inclusive JD:
        // a reform-day morning is Gregorian input, so the Julian reading is
        // Oct 5 (not a spurious ~10-day shift of a misclassified input).
        assert_eq!(
            date_conversion(1582, 10, 15, 0.0, jul),
            Ok(CalendarDate {
                year: 1582,
                month: 10,
                day: 5,
                hour: 0.0,
            })
        );
        assert!(matches!(
            date_conversion(2024, 1, 1, f64::NAN, greg),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn checked_date_conversion_reports_validity_without_rejecting() {
        // Valid leap day: echo matches, JD is noon-based.
        assert_eq!(
            convert_date_checked(2024, 2, 29, 12.0, GREG),
            Ok(CheckedDate {
                valid: true,
                jd: 2460370.0,
                date: CalendarDate {
                    year: 2024,
                    month: 2,
                    day: 29,
                    hour: 12.0,
                },
            })
        );
        // Impossible date: no error, but valid == false with the normalized
        // echo (Feb 29 2023 -> Mar 1), using the same physical model wrapper.
        let bad = convert_date_checked(2023, 2, 29, 12.0, GREG).unwrap();
        assert!(!bad.valid);
        assert_eq!(
            bad.date,
            CalendarDate {
                year: 2023,
                month: 3,
                day: 1,
                hour: 12.0,
            }
        );
        assert_eq!(bad.jd, julday(2023, 2, 29, 12.0, GREG).unwrap());
        // Month 13 flows through julday and fails the round trip.
        assert!(!convert_date_checked(2024, 13, 1, 12.0, GREG).unwrap().valid);
        // The hour is excluded from the validity comparison.
        assert!(convert_date_checked(2024, 1, 1, 6.5, GREG).unwrap().valid);
    }

    #[test]
    fn calendar_code_conversion_retains_reform_and_normalized_date_contracts() {
        let julian = convert_date_code(1582, 10, 5, 12.0, "J").unwrap();
        let gregorian = convert_date_code(1582, 10, 15, 12.0, "g").unwrap();
        assert!(julian.valid && gregorian.valid);
        assert_eq!(julian.jd, gregorian.jd);
        let invalid = convert_date_code(2024, 2, 30, 12.0, "G").unwrap();
        assert!(!invalid.valid);
        assert_eq!(
            (invalid.date.year, invalid.date.month, invalid.date.day),
            (2024, 3, 1)
        );
        assert!(convert_date_code(0, 2, 29, 0.0, "g").unwrap().valid);
        for code in ["", "gg", "x", "gregorian", "é"] {
            assert!(matches!(
                convert_date_code(2000, 1, 1, 12.0, code),
                Err(Error::InvalidCalendarCode { .. })
            ));
        }
        assert!(matches!(
            convert_date_code(2000, 1, 1, f64::NAN, "g"),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn day_of_week_matches_reference_indexing() {
        // Pinned: J2000 Saturday == 5; 2024-01-01 Monday == 0;
        // 1582-10-15 Gregorian Friday == 4; JD 0 Monday == 0.
        assert_eq!(day_of_week(2451545.0), Ok(5));
        assert_eq!(day_of_week(2460311.0), Ok(0));
        assert_eq!(day_of_week(2299160.5), Ok(4));
        assert_eq!(day_of_week(0.0), Ok(0));
        // Euclidean modulo stays nonnegative before the reference epoch.
        assert_eq!(day_of_week(-1.5), Ok(6));
        // Full-week sweep around J2000.
        for (offset, want) in [
            (-3.0, 2),
            (-2.0, 3),
            (-1.0, 4),
            (0.0, 5),
            (1.0, 6),
            (2.0, 0),
            (3.0, 1),
        ] {
            assert_eq!(day_of_week(2451545.0 + offset), Ok(want));
        }
        assert!(matches!(
            day_of_week(f64::NAN),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            day_of_week(f64::INFINITY),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn utc_time_zone_subtracts_offset_and_crosses_boundaries() {
        // the model doc examples: CET (+1) and EST (-5).
        assert_eq!(
            utc_time_zone(2024, 1, 15, 12, 0, 0.0, 1.0),
            Ok(UtcFields {
                year: 2024,
                month: 1,
                day: 15,
                hour: 11,
                minute: 0,
                second: 0.0,
            })
        );
        assert_eq!(
            utc_time_zone(2024, 1, 15, 12, 0, 0.0, -5.0),
            Ok(UtcFields {
                year: 2024,
                month: 1,
                day: 15,
                hour: 17,
                minute: 0,
                second: 0.0,
            })
        );
        // Day rollback / advance across month and year ends.
        assert_eq!(
            utc_time_zone(2024, 1, 1, 0, 30, 0.0, 1.0),
            Ok(UtcFields {
                year: 2023,
                month: 12,
                day: 31,
                hour: 23,
                minute: 30,
                second: 0.0,
            })
        );
        assert_eq!(
            utc_time_zone(2024, 12, 31, 23, 30, 0.0, -2.0),
            Ok(UtcFields {
                year: 2025,
                month: 1,
                day: 1,
                hour: 1,
                minute: 30,
                second: 0.0,
            })
        );
        // Fractional zone (IST +5:30).
        assert_eq!(
            utc_time_zone(2024, 1, 15, 12, 0, 0.0, 5.5),
            Ok(UtcFields {
                year: 2024,
                month: 1,
                day: 15,
                hour: 6,
                minute: 30,
                second: 0.0,
            })
        );
        // Sub-minute residue shifts only the seconds field...
        assert_eq!(
            utc_time_zone(2024, 1, 15, 12, 0, 30.0, 1.0 / 3600.0),
            Ok(UtcFields {
                year: 2024,
                month: 1,
                day: 15,
                hour: 12,
                minute: 0,
                second: 29.0,
            })
        );
        // ... borrowing a minute when the residue underflows...
        assert_eq!(
            utc_time_zone(2024, 1, 15, 12, 0, 0.5, 1.0 / 3600.0),
            Ok(UtcFields {
                year: 2024,
                month: 1,
                day: 15,
                hour: 11,
                minute: 59,
                second: 59.5,
            })
        );
        // ... and a real leap-second label passes through untouched.
        assert_eq!(
            utc_time_zone(2016, 12, 31, 23, 59, 60.0, 0.0),
            Ok(UtcFields {
                year: 2016,
                month: 12,
                day: 31,
                hour: 23,
                minute: 59,
                second: 60.0,
            })
        );
        assert!(matches!(
            utc_time_zone(2024, 1, 15, 12, 0, 0.0, f64::NAN),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            utc_time_zone(2024, 1, 15, 12, 0, f64::INFINITY, 0.0),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    // -- leap-second assets and UTC/TAI transitions ----------------

    fn leap_table() -> Vec<LeapSecondEntry> {
        hardcoded_leap_seconds()
    }

    #[test]
    fn hardcoded_table_matches_release_shape_and_bounds() {
        let table = leap_table();
        // Oracle pin: 28 rows, 1972/10.0 first, 2017/37.0 last.
        assert_eq!(table.len(), 28);
        assert_eq!(
            (
                table[0].year(),
                table[0].month(),
                table[0].day(),
                table[0].tai_utc()
            ),
            (1972, 1, 1, 10.0)
        );
        assert_eq!(
            (
                table[27].year(),
                table[27].month(),
                table[27].day(),
                table[27].tai_utc()
            ),
            (2017, 1, 1, 37.0)
        );
        let mut prev = f64::NEG_INFINITY;
        for row in &table {
            assert!(row.mjd() > prev, "table must be strictly ascending");
            prev = row.mjd();
            // MJD column is derived by the same rule the source uses.
            assert_eq!(
                row.mjd(),
                calendar_to_mjd(row.year(), row.month(), row.day())
            );
        }
        assert_eq!(table[0].mjd(), 41317.0);
    }

    #[test]
    fn tai_utc_lookup_matches_release_pins() {
        let table = leap_table();
        // Oracle pins (project-owned the model outputs, ephemeral scratch).
        assert_eq!(tai_utc_for_mjd(41317.0, &table), Ok(10.0));
        assert_eq!(tai_utc_for_mjd(57754.0, &table), Ok(37.0));
        assert_eq!(tai_utc_for_mjd(57753.0, &table), Ok(36.0));
        assert_eq!(tai_utc_for_jd(2451545.0, &table), Ok(32.0));
        assert_eq!(tai_utc_for_jd(2441317.5, &table), Ok(10.0));
        assert_eq!(tai_utc_for_jd(2457754.5, &table), Ok(37.0));
        // Exact float values from the oracle run.
        assert_eq!(
            tai_utc_for_mjd(41316.0, &table).unwrap().to_bits(),
            0x4023_c780_346d_c5d6
        );
        assert_eq!(tai_utc_for_mjd(30000.0, &table), Ok(1.4228180));
        // Pre-1961 clamp: UTC undefined, earliest offset held.
        assert_eq!(tai_utc_for_mjd(-100_000.0, &table), Ok(1.4228180));
        // Empty table: the source's 37/32 fallback.
        let empty: Vec<LeapSecondEntry> = Vec::new();
        assert_eq!(tai_utc_for_mjd(57754.0, &empty), Ok(37.0));
        assert_eq!(tai_utc_for_mjd(50000.0, &empty), Ok(32.0));
        // Documented strictness: the source would propagate these.
        assert!(matches!(
            tai_utc_for_jd(f64::NAN, &table),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            tai_utc_for_mjd(f64::INFINITY, &table),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn mjd_helpers_match_release_pins() {
        // Oracle pins.
        assert_eq!(calendar_to_mjd(1972, 1, 1), 41317.0);
        assert_eq!(calendar_to_mjd(2000, 1, 1), 51544.0);
        // Negative year: truncation (not floor) in the century term.
        assert_eq!(calendar_to_mjd(-100, 1, 1), -715465.0);
        assert_eq!(mjd_to_calendar(41317.0), Ok((1972, 1, 1)));
        assert_eq!(mjd_to_calendar(58000.75), Ok((2017, 9, 4)));
        assert_eq!(mjd_to_calendar(100000.0), Ok((2132, 9, 1)));
        // Pre-reform branch (z < 2299161).
        assert_eq!(mjd_to_calendar(-1000.0), Ok((1856, 2, 21)));
        // Round-trip holds where both directions use the Gregorian path.
        for &(y, m, d) in &[(1972, 1, 1), (2000, 2, 29), (2016, 12, 31), (1856, 2, 21)] {
            assert_eq!(mjd_to_calendar(calendar_to_mjd(y, m, d)), Ok((y, m, d)));
        }
        // Documented source asymmetry: `_calendar_to_mjd` is pure
        // proleptic-Gregorian while `_mjd_to_calendar` switches to the
        // Julian-style branch for JD < 2299161, so ancient round-trips
        // relabel the date (oracle pins — the Rust code reproduces them).
        assert_eq!(
            mjd_to_calendar(calendar_to_mjd(1582, 10, 4)),
            Ok((1582, 9, 24))
        );
        assert_eq!(
            mjd_to_calendar(calendar_to_mjd(-100, 1, 1)),
            Ok((-100, 1, 4))
        );
        assert_eq!(jd_to_mjd(2451545.0), Ok(51544.5));
        assert_eq!(mjd_to_jd(51544.5), Ok(2451545.0));
        assert!(matches!(
            mjd_to_calendar(f64::NAN),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            jd_to_mjd(f64::INFINITY),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn tt_tai_converters_are_exact_inverse_pairs() {
        // Oracle pin: J2000 TT -> TAI value.
        assert_eq!(
            tt_to_tai_jd(2451545.0).unwrap().to_bits(),
            0x4142_b42c_7ff3_cb3e
        );
        assert_eq!(TT_TAI_OFFSET_SECONDS, 32.184);
        assert_eq!(TT_TAI_OFFSET_DAYS, 32.184 / 86400.0);
        // Round-trip at several epochs (one rounding each way at most).
        for &jd in &[2400000.5, 2441317.5, 2451545.0, 2458849.5, 2460000.5] {
            let tai = tt_to_tai_jd(jd).unwrap();
            let back = tai_to_tt_jd(tai).unwrap();
            assert!(
                (back - jd).abs() <= f64::EPSILON * jd.abs(),
                "TT/TAI round-trip drifted at {jd}"
            );
        }
        assert!(matches!(
            tt_to_tai_jd(f64::NAN),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            tai_to_tt_jd(f64::INFINITY),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn utc_tai_converters_match_release_pins() {
        let table = leap_table();
        // Oracle pins (shortest-round-trip literals from the model run).
        assert_eq!(
            utc_to_tai_jd(2020, 1, 1, 0, 0, 0.0, GREG, &table)
                .unwrap()
                .to_bits(),
            0x4142_c270_c00e_0858
        );
        assert_eq!(
            utc_to_tai_jd(2000, 1, 1, 12, 0, 0.0, JUL, &table),
            Ok(2451558.00037037)
        );
        assert_eq!(
            tai_jd_to_utc(
                utc_to_tai_jd(2020, 1, 1, 0, 0, 0.0, GREG, &table).unwrap(),
                GREG,
                &table
            ),
            Ok(UtcFields {
                year: 2020,
                month: 1,
                day: 1,
                hour: 0,
                minute: 0,
                second: 0.0,
            })
        );
        // A 23:59:60 label on a leap date lands on the next midnight UTC.
        assert_eq!(
            tai_jd_to_utc(
                utc_to_tai_jd(2016, 12, 31, 23, 59, 60.0, GREG, &table).unwrap(),
                GREG,
                &table
            ),
            Ok(UtcFields {
                year: 2017,
                month: 1,
                day: 1,
                hour: 0,
                minute: 0,
                second: 0.0,
            })
        );
        assert_eq!(utc_jd_to_tai(2451545.0, &table), Ok(2451545.00037037));
        assert_eq!(
            tai_to_utc_jd(2451545.0 + 37.0 / 86400.0, &table),
            Ok(2451545.0000578705)
        );
        // JD-level and calendar-level paths agree.
        let jd_tai = utc_to_tai_jd(2023, 6, 15, 12, 30, 15.5, GREG, &table).unwrap();
        assert_eq!(
            tai_to_utc_jd(jd_tai, &table).unwrap(),
            jd_tai - 37.0 / 86400.0
        );
        assert_eq!(
            utc_jd_to_tai(tai_to_utc_jd(jd_tai, &table).unwrap(), &table),
            Ok(jd_tai)
        );
    }

    #[test]
    fn leap_parsers_match_release_behavior() {
        // Numeric IERS rows with comments and blanks.
        let rows = parse_leap_seconds_iers(
            "#  MJD        Date        TAI-UTC (s)\n\
             #           day month year\n\
             \n\
             41317.0    1  1 1972      10\n\
             41499.0    1  7 1972      11\n",
        );
        assert_eq!(rows.len(), 2);
        assert_eq!(
            (
                rows[0].mjd(),
                rows[0].tai_utc(),
                rows[0].year(),
                rows[0].month(),
                rows[0].day()
            ),
            (41317.0, 10.0, 1972, 1, 1)
        );
        assert_eq!(rows[1].mjd(), 41499.0);
        // Tolerated textual variant (oracle-verified sample shape).
        let text = parse_leap_seconds_iers("1972 JAN 1 foo 10s");
        assert_eq!(text.len(), 1);
        assert_eq!(
            (
                text[0].mjd(),
                text[0].tai_utc(),
                text[0].year(),
                text[0].month(),
                text[0].day()
            ),
            (41317.0, 10.0, 1972, 1, 1)
        );
        // IETF rows (oracle-verified MJDs and labels).
        let ietf =
            parse_leap_seconds_ietf("#$ 3814959996\n#@ 3876959996\n2272060800 10\n2287785600 11\n");
        assert_eq!(ietf.len(), 2);
        assert_eq!(
            (
                ietf[0].mjd(),
                ietf[0].tai_utc(),
                ietf[0].year(),
                ietf[0].month(),
                ietf[0].day()
            ),
            (41317.0, 10.0, 1972, 1, 1)
        );
        assert_eq!(ietf[1].mjd(), 41499.0);
        // Dispatcher: IETF markers route to the IETF parser...
        let via_ietf = parse_leap_seconds("#$ 3814959996\n2272060800 10\n");
        assert_eq!(via_ietf.len(), 1);
        assert_eq!(via_ietf[0].mjd(), 41317.0);
        // ... a misdetected IERS file still parses via the fallback...
        let misdetected =
            parse_leap_seconds("# NTP-looking comment\n41317.0    1  1 1972      10\n");
        assert_eq!(misdetected.len(), 1);
        assert_eq!(misdetected[0].tai_utc(), 10.0);
        // ... and unparsable content falls back to the packaged table.
        let fallback = parse_leap_seconds("garbage line one\ngarbage two\n");
        assert_eq!(fallback.len(), 28);
        assert_eq!(fallback[27].tai_utc(), 37.0);
        // Results are sorted ascending by MJD.
        let unsorted =
            parse_leap_seconds("41499.0    1  7 1972      11\n41317.0    1  1 1972      10\n");
        assert_eq!(unsorted[0].mjd(), 41317.0);
        assert_eq!(unsorted[1].mjd(), 41499.0);
        // Constructor strictness: non-finite rows are never admitted.
        assert!(matches!(
            LeapSecondEntry::new(f64::NAN, 10.0, 1972, 1, 1),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            LeapSecondEntry::new(41317.0, f64::INFINITY, 1972, 1, 1),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    // -- two-part TAI/TT scale conversions ---------------------------

    /// One collapsed ulp at a finite nonzero `x` (`next_up` needs Rust 1.86+,
    /// our MSRV). At JD magnitudes this is ~4.7 e-10 d (~40 us).
    fn ulp(x: f64) -> f64 {
        (x.next_up() - x).abs()
    }

    /// A priori bound on `|converted.diff(input) - offset|` for split inputs
    /// at epoch magnitudes: the double-double addition rounds once at the
    /// `lo`-part scale (`|lo| <= 2^-25 * |jd|` for split outputs, hence a
    /// `2^-78 * |jd|` rounding) and `diff` rounds its two subtractions at
    /// the same scale; `2^-75 * |jd|` covers all three with margin. This is
    /// an absolute-day bound (~3 e-17 d, ~3 ps, at JD 2.4 e 6), not a fitted
    /// tolerance.
    fn offset_identity_bound(jd: f64) -> f64 {
        2f64.powi(-75) * jd.abs()
    }

    #[test]
    fn two_part_tt_tai_offset_identity_and_round_trip() {
        let epochs = [
            2400000.5, 2441317.5, 2451545.0, 2458849.5, 2460000.5, 2500000.75,
        ];
        for &jd in &epochs {
            let split = TwoPartEpoch::split(jd).unwrap();
            // Exact-offset identity in both directions.
            let to_tt = tai_to_tt_two_part(split).unwrap();
            assert!(
                (to_tt.diff(split) - TT_TAI_OFFSET_DAYS).abs() <= offset_identity_bound(jd),
                "TAI->TT offset identity failed at {jd}"
            );
            let to_tai = tt_to_tai_two_part(split).unwrap();
            assert!(
                (to_tai.diff(split) + TT_TAI_OFFSET_DAYS).abs() <= offset_identity_bound(jd),
                "TT->TAI offset identity failed at {jd}"
            );
            // Collapsed agreement with the model-agreed scalar path.
            assert!(
                (to_tt.value() - tai_to_tt_jd(jd).unwrap()).abs() <= ulp(jd),
                "TAI->TT scalar agreement failed at {jd}"
            );
            assert!(
                (to_tai.value() - tt_to_tai_jd(jd).unwrap()).abs() <= ulp(jd),
                "TT->TAI scalar agreement failed at {jd}"
            );
            // Round-trips recover the epoch (two double-double legs plus one
            // collapse rounding stay within a single collapsed ulp).
            let round_tt = tai_to_tt_two_part(tt_to_tai_two_part(split).unwrap()).unwrap();
            assert!(
                (round_tt.value() - jd).abs() <= ulp(jd),
                "TT round-trip drifted at {jd}"
            );
            let round_tai = tt_to_tai_two_part(tai_to_tt_two_part(split).unwrap()).unwrap();
            assert!(
                (round_tai.value() - jd).abs() <= ulp(jd),
                "TAI round-trip drifted at {jd}"
            );
        }
        // J2000 pin: the two-part collapse lands on the scalar value's ulp.
        let j2000 = TwoPartEpoch::split(2451545.0).unwrap();
        assert_eq!(
            tt_to_tai_two_part(j2000).unwrap().value().to_bits(),
            tt_to_tai_jd(2451545.0).unwrap().to_bits()
        );
    }

    #[test]
    fn two_part_utc_tai_agrees_with_scalar_across_leap_boundaries() {
        let table = leap_table();
        let second = 1.0 / 86400.0;
        // Every packaged leap midnight: exact boundary, adjacent float
        // values, and +/-1 s / +/-1 d neighborhoods, both JD directions.
        for row in &table {
            let boundary = mjd_to_jd(row.mjd()).unwrap();
            let mut cases = vec![
                boundary,
                boundary - second,
                boundary + second,
                boundary - 1.0,
                boundary + 1.0,
            ];
            for &base in &[boundary, boundary - second, boundary + second] {
                cases.push(base.next_up());
                cases.push(base.next_down());
            }
            for jd in cases {
                let split = TwoPartEpoch::split(jd).unwrap();
                // Identical collapsed keys select identical rows by
                // construction (split is exact); pin that contract.
                assert_eq!(
                    tai_utc_for_two_part(split, &table),
                    tai_utc_for_jd(jd, &table),
                    "row selection diverged at {jd}"
                );
                assert!(
                    (tai_to_utc_jd_two_part(split, &table).unwrap().value()
                        - tai_to_utc_jd(jd, &table).unwrap())
                    .abs()
                        <= ulp(jd),
                    "TAI->UTC agreement failed at {jd}"
                );
                assert!(
                    (utc_jd_to_tai_two_part(split, &table).unwrap().value()
                        - utc_jd_to_tai(jd, &table).unwrap())
                    .abs()
                        <= ulp(jd),
                    "UTC->TAI agreement failed at {jd}"
                );
            }
        }
        // Monthly noon grid across the whole leap-second era, plus the
        // pre-1972 rate-model region and the empty-table fallback.
        let mut grid = Vec::new();
        for year in 1972..=2025 {
            for month in 1..=12 {
                grid.push(julday(year, month, 1, 12.0, GREG).unwrap());
            }
        }
        grid.push(julday(1960, 1, 1, 12.0, GREG).unwrap());
        grid.push(julday(1965, 6, 15, 12.0, GREG).unwrap());
        grid.push(julday(1971, 12, 31, 12.0, GREG).unwrap());
        for jd in grid {
            let split = TwoPartEpoch::split(jd).unwrap();
            assert!(
                (tai_to_utc_jd_two_part(split, &table).unwrap().value()
                    - tai_to_utc_jd(jd, &table).unwrap())
                .abs()
                    <= ulp(jd),
                "TAI->UTC grid agreement failed at {jd}"
            );
            assert!(
                (utc_jd_to_tai_two_part(split, &table).unwrap().value()
                    - utc_jd_to_tai(jd, &table).unwrap())
                .abs()
                    <= ulp(jd),
                "UTC->TAI grid agreement failed at {jd}"
            );
        }
        let empty: Vec<LeapSecondEntry> = Vec::new();
        for &jd in &[2451545.0, 2457754.5] {
            let split = TwoPartEpoch::split(jd).unwrap();
            assert!(
                (tai_to_utc_jd_two_part(split, &empty).unwrap().value()
                    - tai_to_utc_jd(jd, &empty).unwrap())
                .abs()
                    <= ulp(jd)
            );
        }
    }

    #[test]
    fn two_part_calendar_and_fields_paths() {
        let table = leap_table();
        // Calendar-label path agrees with the scalar path within one ulp.
        for &(y, mo, d, h, mi, s) in &[
            (2020, 1, 1, 0, 0, 0.0),
            (2023, 6, 15, 12, 30, 15.5),
            (2016, 12, 31, 23, 59, 60.0),
            (2000, 1, 1, 12, 0, 0.0),
        ] {
            let two = utc_to_tai_two_part(y, mo, d, h, mi, s, GREG, &table).unwrap();
            let one = utc_to_tai_jd(y, mo, d, h, mi, s, GREG, &table).unwrap();
            assert!(
                (two.value() - one).abs() <= ulp(one),
                "calendar path agreement failed at {y}-{mo}-{d}"
            );
        }
        // A 23:59:60 label on a leap date lands on the next midnight UTC.
        let leap_tai = utc_to_tai_two_part(2016, 12, 31, 23, 59, 60.0, GREG, &table).unwrap();
        assert_eq!(
            tai_two_part_to_utc(leap_tai, GREG, &table),
            Ok(UtcFields {
                year: 2017,
                month: 1,
                day: 1,
                hour: 0,
                minute: 0,
                second: 0.0,
            })
        );
        // Fields path is exactly-equal to the scalar path on the collapsed
        // JD (same operations on the same values), both calendars.
        for &jd in &[2451545.0, 2458849.5, 2441317.5] {
            for &cal in &[GREG, JUL] {
                assert_eq!(
                    tai_two_part_to_utc(TwoPartEpoch::split(jd).unwrap(), cal, &table),
                    tai_jd_to_utc(jd, cal, &table),
                    "fields path diverged at {jd} {cal:?}"
                );
            }
        }
    }

    #[test]
    fn two_part_scale_conversions_reject_unrepresentable() {
        let table = leap_table();
        // Hostile directly-constructed parts (fields are public): the
        // collapsed key is infinite, so lookup refuses it exactly as the
        // scalar path refuses infinite JD.
        let huge = TwoPartEpoch {
            hi: f64::MAX,
            lo: f64::MAX,
        };
        assert!(matches!(
            tai_utc_for_two_part(huge, &table),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            tai_to_utc_jd_two_part(huge, &table),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            utc_jd_to_tai_two_part(huge, &table),
            Err(Error::NonFiniteArgument { .. })
        ));
        // The fixed-offset legs overflow through `add_offset`'s guard
        // instead of collapsing to infinity.
        assert!(matches!(
            tt_to_tai_two_part(huge),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            tai_to_tt_two_part(huge),
            Err(Error::NonFiniteArgument { .. })
        ));
        let nan = TwoPartEpoch {
            hi: f64::NAN,
            lo: 0.0,
        };
        assert!(matches!(
            tai_utc_for_two_part(nan, &table),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            tai_to_tt_two_part(nan),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    // --- tests --------------------------------------------------
    //
    // exactly oracle pins below are project-owned outputs of the pinned
    // the model (`37ffb05...`) obtained ephemerally through the VAL-02
    // adapter (`deltat_espenak_meeus_poly`, `tid_acc_adjustment_seconds`,
    // `deltat_espenak`, `deltat_ex_espenak`, `deltat_espenak_tidacc`,
    // `deltat_ex_espenak_tidacc`, `deltat_userdef` probes); nothing is
    // stored outside these literals. They pin branch ownership (every
    // exact threshold year belongs to the later `<` branch), the C-`pow`
    // operation order, and the end-to-end priority chain.

    use crate::constants::{FLG_JPLEPH, FLG_SWIEPH, TIDAL_DE421, TIDAL_DE440, TIDAL_MOSEPH};

    fn espenak_cfg() -> DeltaTConfig<'static> {
        DeltaTConfig {
            model: DeltaTModel::EspenakMeeus,
            user_override_days: None,
            tidal_acceleration: TIDAL_DEFAULT,
            iers_enabled: false,
            iers_lookup_seconds: None,
            smh2016_lookup_seconds: None,
        }
    }

    #[test]
    fn delta_t_model_names() {
        assert_eq!(DeltaTModel::Smh2016.as_name(), "smh2016");
        assert_eq!(DeltaTModel::EspenakMeeus.as_name(), "espenak_meeus");
        assert_eq!(
            DeltaTModel::from_name("smh2016"),
            Some(DeltaTModel::Smh2016)
        );
        assert_eq!(
            DeltaTModel::from_name("espenak_meeus"),
            Some(DeltaTModel::EspenakMeeus)
        );
        // `set_delta_t_model` normalization (case/whitespace).
        assert_eq!(
            DeltaTModel::from_name("SMH2016"),
            Some(DeltaTModel::Smh2016)
        );
        assert_eq!(
            DeltaTModel::from_name("  Espenak_Meeus\n"),
            Some(DeltaTModel::EspenakMeeus)
        );
        assert_eq!(DeltaTModel::from_name("smh-2016"), None);
        assert_eq!(DeltaTModel::from_name(""), None);
    }

    #[test]
    fn espenak_meeus_oracle_pins() {
        // (decimal year, expected IEEE-754 bits of Delta T in seconds).
        // Every branch threshold appears exactly (later-branch ownership)
        // plus interiors on both sides of each seam.
        let pins: &[(f64, u64)] = &[
            (-2500.0, 0x40ed2675c28f5c2a),
            (-1000.0, 0x40d8d4eb851eb852),
            (-500.0, 0x40d0ccea01758e22),
            (-499.5, 0x40d0caa7718918b9),
            (0.0, 0x40c4abcccccccccd),
            (499.5, 0x40b65300fdcfcb2a),
            (500.0, 0x40b64e0b6f837b4a),
            (750.0, 0x40aa5b8654b020c6),
            (1599.5, 0x405e1e473f635d24),
            (1600.0, 0x405e000000000000),
            (1650.0, 0x404918d584182bea),
            (1699.5, 0x4021d5909b41d580),
            (1700.0, 0x4021a8f5c28f5c29),
            (1750.0, 0x402abd79d9cbd72f),
            (1799.5, 0x402bef7fa6f4cf08),
            (1800.0, 0x402b70a3d70a3d71),
            (1830.0, 0x401eb18a86d71f2c),
            (1859.5, 0x401da3b0d3d58c00),
            (1860.0, 0x401e7ae147ae147b),
            (1880.0, 0xc01408b0d02828c0),
            (1899.5, 0xc00b18f135f84b00),
            (1900.0, 0xc00651eb851eb852),
            (1910.0, 0x4024c6dc5d638865),
            (1919.5, 0x4035088f561408cd),
            (1920.0, 0x4035333333333333),
            (1930.0, 0x40382205bc01a36e),
            (1940.5, 0x403893a05d023997),
            (1941.0, 0x4038c5ec98d522d9),
            (1955.0, 0x403f0bf9da70d306),
            (1960.5, 0x4040a99401f4a22a),
            (1961.0, 0x4040ca398944c6d8),
            (1970.0, 0x404418b24940d118),
            (1985.5, 0x404b4eff70fdefac),
            (1986.0, 0x404b7059b4229b2e),
            (2000.0, 0x404fee147ae147ae),
            (2004.5, 0x4050271d8e7eeb51),
            (2005.0, 0x40502aeab367a0f9),
            (2010.0, 0x4050acd6a161e4f7),
            (2030.0, 0x4053675f6fd21ff3),
            (2049.5, 0x405723f58d1526d9),
            (2050.0, 0x40573ffffffffffe),
            (2100.0, 0x406957ae147ae147),
            (2149.5, 0x4074724b5dcc63f1),
            (2150.0, 0x407487ae147ae147),
            (2200.0, 0x407ba147ae147ae1),
            (3000.0, 0x40b153ae147ae148),
        ];
        assert_eq!(pins.len(), 46);
        for (year, bits) in pins {
            assert_eq!(
                deltat_espenak_meeus(*year).to_bits(),
                *bits,
                "decimal year {year}"
            );
        }
    }

    #[test]
    fn tid_acc_reference_and_adjustment_pins() {
        // Derivation identity (J2000 minus 44.5 Julian years), exact.
        assert_eq!(TID_ACC_REFERENCE_JD, 2435291.375);
        // Zero policy: default acceleration, or any epoch on/after the
        // reference, is exactly +0.0.
        assert_eq!(
            tid_acc_adjustment_seconds(1500000.0, TIDAL_DEFAULT).to_bits(),
            0x0000000000000000
        );
        assert_eq!(
            tid_acc_adjustment_seconds(2451545.0, -25.0).to_bits(),
            0x0000000000000000
        );
        assert_eq!(
            tid_acc_adjustment_seconds(TID_ACC_REFERENCE_JD, -25.0).to_bits(),
            0x0000000000000000
        );
        // Nonzero oracle pins (seconds): sign follows -dacc*u^2.
        let pins: &[(f64, f64, u64)] = &[
            (1500000.0, -25.0, 0xc07ddbc6be809db8),
            (1500000.0, TIDAL_DE421, 0x403ddbc6be809e28),
            (1500000.0, TIDAL_DE440, 0x40544dc4960589d2),
            (1500000.0, TIDAL_MOSEPH, 0xc0606c13b5938a20),
            (TID_ACC_REFERENCE_JD - 0.001, -25.0, 0xbcc3ad23c5c5397b),
        ];
        for (jd, acc, bits) in pins {
            assert_eq!(
                tid_acc_adjustment_seconds(*jd, *acc).to_bits(),
                *bits,
                "jd {jd} acc {acc}"
            );
        }
    }

    #[test]
    fn deltat_priority_chain_with_stub_providers() {
        use std::cell::Cell;
        let iers_called = Cell::new(false);
        let smh_called = Cell::new(false);
        let iers_hit = |_jd: f64| {
            iers_called.set(true);
            Some(70.0)
        };
        let smh_hit = |_jd: f64| {
            smh_called.set(true);
            Some(71.0)
        };
        let jd = 2451545.0;

        // User override wins; no provider is consulted.
        let user_cfg = DeltaTConfig {
            user_override_days: Some(65.0 / 86400.0),
            iers_enabled: true,
            iers_lookup_seconds: Some(&iers_hit),
            smh2016_lookup_seconds: Some(&smh_hit),
            ..espenak_cfg()
        };
        assert_eq!(deltat_with_config(jd, &user_cfg), Ok(65.0 / 86400.0));
        assert!(!iers_called.get() && !smh_called.get());

        // IERS hit outranks the model; the tidal adjustment still applies
        // (zero here: modern epoch).
        let iers_cfg = DeltaTConfig {
            iers_enabled: true,
            iers_lookup_seconds: Some(&iers_hit),
            smh2016_lookup_seconds: Some(&smh_hit),
            ..espenak_cfg()
        };
        assert_eq!(deltat_with_config(jd, &iers_cfg), Ok(70.0 / 86400.0));
        assert!(iers_called.get() && !smh_called.get());

        // IERS miss falls through to the model.
        iers_called.set(false);
        smh_called.set(false);
        let miss = |_jd: f64| {
            iers_called.set(true);
            None
        };
        let miss_cfg = DeltaTConfig {
            iers_enabled: true,
            iers_lookup_seconds: Some(&miss),
            ..espenak_cfg()
        };
        assert_eq!(
            deltat_with_config(jd, &miss_cfg),
            deltat_with_config(jd, &espenak_cfg())
        );
        assert!(iers_called.get());

        // Disabled IERS never consults its provider, even one that would hit.
        iers_called.set(false);
        let disabled_cfg = DeltaTConfig {
            iers_lookup_seconds: Some(&iers_hit),
            ..espenak_cfg()
        };
        assert_eq!(
            deltat_with_config(jd, &disabled_cfg),
            deltat_with_config(jd, &espenak_cfg())
        );
        assert!(!iers_called.get());

        // smh 2016 routes to its provider with the same adjustment rule.
        let smh_cfg = DeltaTConfig {
            smh2016_lookup_seconds: Some(&smh_hit),
            ..DeltaTConfig::default()
        };
        assert_eq!(deltat_with_config(jd, &smh_cfg), Ok(71.0 / 86400.0));
        assert!(smh_called.get());

        // smh 2016 with no usable backend is an honest error, never a
        // silent substitution of the other model.
        assert_eq!(
            deltat_with_config(jd, &DeltaTConfig::default()),
            Err(Error::DeltaTModelUnavailable { model: "smh2016" })
        );
        let smh_miss_cfg = DeltaTConfig {
            smh2016_lookup_seconds: Some(&miss),
            ..DeltaTConfig::default()
        };
        assert!(matches!(
            deltat_with_config(jd, &smh_miss_cfg),
            Err(Error::DeltaTModelUnavailable { .. })
        ));
        assert!(!format!("{}", Error::DeltaTModelUnavailable { model: "smh2016" }).is_empty());
    }

    #[test]
    fn deltat_end_to_end_oracle_pins() {
        // Delta T in days end to end (espenak model, default state).
        let days: &[(f64, u64)] = &[
            (2451545.0, 0x3f483988b9c48c47),
            (2378497.0, 0x3f24caf1d7dc58b3),
            (1500000.0, 0x3fcbda8ab8a30b5f),
            (2460000.5, 0x3f4bd2737ebc7eb1),
            (2435291.0, 0x3f37b62e8e1cf989),
        ];
        for (jd, bits) in days {
            assert_eq!(
                deltat_with_config(*jd, &espenak_cfg()).unwrap().to_bits(),
                *bits
            );
        }
        // deltat_ex: modern epochs are adjustment-free under every flag;
        // ancient epochs observe the MOSEPH acceleration.
        assert_eq!(
            deltat_ex_with_config(2451545.0, FLG_MOSEPH, &espenak_cfg(), None)
                .unwrap()
                .to_bits(),
            0x3f483988b9c48c47
        );
        assert_eq!(
            deltat_ex_with_config(2451545.0, FLG_SWIEPH, &espenak_cfg(), None)
                .unwrap()
                .to_bits(),
            0x3f483988b9c48c47
        );
        assert_eq!(
            deltat_ex_with_config(2378497.0, FLG_MOSEPH, &espenak_cfg(), None)
                .unwrap()
                .to_bits(),
            0x3f240eced640a8c1
        );
        assert_eq!(
            deltat_ex_with_config(2378497.0, FLG_SWIEPH, &espenak_cfg(), None)
                .unwrap()
                .to_bits(),
            0x3f24caf1d7dc58b3
        );
        assert_eq!(
            deltat_ex_with_config(1500000.0, FLG_MOSEPH, &espenak_cfg(), None)
                .unwrap()
                .to_bits(),
            0x3fcba8b73d804cdd
        );
        // Explicit acceleration threads through both entry points.
        let explicit_cfg = DeltaTConfig {
            tidal_acceleration: TIDAL_DE440,
            ..espenak_cfg()
        };
        assert_eq!(
            deltat_with_config(1500000.0, &explicit_cfg)
                .unwrap()
                .to_bits(),
            0x3fcbf957ed7c0377
        );
        assert_eq!(
            deltat_ex_with_config(1500000.0, FLG_MOSEPH, &espenak_cfg(), Some(-25.0))
                .unwrap()
                .to_bits(),
            0x3fcb255b276a852b
        );
        // The user override is returned (never adjusted).
        let user_cfg = DeltaTConfig {
            user_override_days: Some(65.0 / 86400.0),
            ..espenak_cfg()
        };
        assert_eq!(
            deltat_with_config(2451545.0, &user_cfg).unwrap().to_bits(),
            0x3f48a6dfc3518a6e
        );
        assert_eq!(
            deltat_ex_with_config(1500000.0, FLG_MOSEPH, &user_cfg, Some(-25.0))
                .unwrap()
                .to_bits(),
            0x3f48a6dfc3518a6e
        );
    }

    #[test]
    fn deltat_ex_flag_selection() {
        // Only the MOSEPH bit is observed; explicit settings win.
        assert_eq!(resolve_tidal_acceleration(FLG_MOSEPH, None), TIDAL_MOSEPH);
        assert_eq!(resolve_tidal_acceleration(0, None), TIDAL_DEFAULT);
        assert_eq!(resolve_tidal_acceleration(FLG_SWIEPH, None), TIDAL_DEFAULT);
        assert_eq!(resolve_tidal_acceleration(FLG_JPLEPH, None), TIDAL_DEFAULT);
        assert_eq!(
            resolve_tidal_acceleration(FLG_SWIEPH | FLG_MOSEPH, None),
            TIDAL_MOSEPH
        );
        assert_eq!(
            resolve_tidal_acceleration(FLG_JPLEPH | FLG_MOSEPH | FLG_SWIEPH, None),
            TIDAL_MOSEPH
        );
        assert_eq!(resolve_tidal_acceleration(FLG_MOSEPH, Some(-25.0)), -25.0);
        assert_eq!(
            resolve_tidal_acceleration(FLG_SWIEPH, Some(TIDAL_MOSEPH)),
            TIDAL_MOSEPH
        );
        // Functional selection agrees with the equivalent direct config.
        let jd = 1500000.0;
        let via_flag = deltat_ex_with_config(jd, FLG_MOSEPH, &espenak_cfg(), None).unwrap();
        let direct = DeltaTConfig {
            tidal_acceleration: TIDAL_MOSEPH,
            ..espenak_cfg()
        };
        assert_eq!(via_flag, deltat_with_config(jd, &direct).unwrap());
    }

    #[test]
    fn deltat_hostile_inputs() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(
                matches!(
                    deltat_with_config(bad, &espenak_cfg()),
                    Err(Error::NonFiniteArgument { .. })
                ),
                "epoch {bad}"
            );
            assert!(
                matches!(
                    deltat_ex_with_config(bad, FLG_MOSEPH, &espenak_cfg(), None),
                    Err(Error::NonFiniteArgument { .. })
                ),
                "ex epoch {bad}"
            );
        }
        // A stored override is returned even for a non-finite
        // epoch.
        let user_cfg = DeltaTConfig {
            user_override_days: Some(0.001),
            ..espenak_cfg()
        };
        assert_eq!(deltat_with_config(f64::NAN, &user_cfg), Ok(0.001));
        assert_eq!(
            deltat_ex_with_config(f64::NAN, FLG_MOSEPH, &user_cfg, Some(-25.0)),
            Ok(0.001)
        );
    }

    // -- UT/UT1 <-> TT pathways ------------------------------------

    fn packaged_table() -> Vec<LeapSecondEntry> {
        hardcoded_leap_seconds()
    }

    #[test]
    fn utc_epoch_constants_match_release() {
        // `_UTC_EPOCH_JD` is 1972-01-01 00:00 UTC and the guard band is an
        // exactly representable ~2 s window.
        assert_eq!(UTC_EPOCH_JD, 2441317.5);
        assert_eq!(UTC_EPOCH_BAND_DAYS, 2.0 / 86400.0);
        assert_eq!(revjul(UTC_EPOCH_JD, GREG).unwrap().year, 1972);
        assert_eq!(
            reconstruction_snap_tolerance(UTC_EPOCH_JD).to_bits(),
            (2.0 * ulp(UTC_EPOCH_JD) * 86400.0).to_bits()
        );
    }

    #[test]
    fn utc_to_jd_validation_order_and_messages() {
        let table = packaged_table();
        let cfg = espenak_cfg();
        // A bad date wins even when the time is also invalid.
        assert_eq!(
            utc_to_jd_with_config(2024, 13, 1, 25, 0, 0.0, GREG, &table, &cfg),
            Err(Error::InvalidDate {
                year: 2024,
                month: 13,
                day: 1
            })
        );
        assert_eq!(
            format!(
                "{}",
                Error::InvalidDate {
                    year: 2024,
                    month: 13,
                    day: 1
                }
            ),
            "invalid date: year = 2024, month = 13, day = 1"
        );
        // Out-of-range clock fields (hour 24 included).
        for (h, m, s) in [(24, 0, 0.0), (12, 60, 0.0), (12, 0, -1.0), (12, 0, 61.0)] {
            assert!(
                matches!(
                    utc_to_jd_with_config(2020, 1, 1, h, m, s, GREG, &table, &cfg),
                    Err(Error::InvalidTime { .. })
                ),
                "clock ({h}, {m}, {s})"
            );
        }
        assert_eq!(
            format!(
                "{}",
                Error::InvalidTime {
                    hour: 24,
                    minute: 0,
                    second: 0.0
                }
            ),
            "invalid time: 24:0:0.00"
        );
        // `60.x` labels: off-23:59 and pre-1972 are invalid times...
        assert!(matches!(
            utc_to_jd_with_config(2020, 1, 1, 12, 0, 60.0, GREG, &table, &cfg),
            Err(Error::InvalidTime { .. })
        ));
        assert!(matches!(
            utc_to_jd_with_config(1960, 6, 15, 23, 59, 60.0, GREG, &table, &cfg),
            Err(Error::InvalidTime { .. })
        ));
        // ... while a valid 23:59:60 label on a date with no insertion is
        // the distinct no-leap-second rejection.
        assert_eq!(
            utc_to_jd_with_config(2023, 5, 17, 23, 59, 60.0, GREG, &table, &cfg),
            Err(Error::NoLeapSecond {
                hour: 23,
                minute: 59,
                second: 60.0
            })
        );
        assert_eq!(
            format!(
                "{}",
                Error::NoLeapSecond {
                    hour: 23,
                    minute: 59,
                    second: 60.0
                }
            ),
            "invalid time (no leap second!): 23:59:60.00"
        );
        // Genuine leap labels are admitted, in both calendars: Julian
        // 2016-12-18 is the physical Gregorian 2016-12-31 insertion, while
        // Julian 2016-12-31 (Gregorian 2017-01-13) carries none.
        assert!(utc_to_jd_with_config(1972, 6, 30, 23, 59, 60.0, GREG, &table, &cfg).is_ok());
        assert!(utc_to_jd_with_config(2016, 12, 31, 23, 59, 60.5, GREG, &table, &cfg).is_ok());
        assert!(utc_to_jd_with_config(2016, 12, 18, 23, 59, 60.0, JUL, &table, &cfg).is_ok());
        assert!(matches!(
            utc_to_jd_with_config(2016, 12, 31, 23, 59, 60.0, JUL, &table, &cfg),
            Err(Error::NoLeapSecond { .. })
        ));
    }

    #[test]
    fn utc_to_jd_outer_branches_match_julday_plus_deltat() {
        let table = packaged_table();
        let cfg = espenak_cfg();
        // Pre-1972 and post-2035 share one structural form: the civil label
        // is UT1, TT adds ΔT. Pin the composition bit-exactly at ancient,
        // 1971, and far-future labels in both calendars.
        for (y, m, d, h, mi, s, cal) in [
            (1850, 3, 15, 12, 0, 0.0, GREG),
            (1971, 12, 31, 23, 59, 59.0, GREG),
            (1960, 6, 15, 0, 0, 0.0, JUL),
            (2040, 1, 1, 0, 0, 0.0, GREG),
            (3000, 7, 4, 18, 30, 15.25, JUL),
        ] {
            let (jd_tt, jd_ut1) =
                utc_to_jd_with_config(y, m, d, h, mi, s, cal, &table, &cfg).unwrap();
            let decimal_hour = h as f64 + mi as f64 / 60.0 + s / 3600.0;
            let expect_ut1 = julday(y, m, d, decimal_hour, cal).unwrap();
            assert_eq!(jd_ut1.to_bits(), expect_ut1.to_bits());
            assert_eq!(
                jd_tt.to_bits(),
                (expect_ut1 + deltat_with_config(expect_ut1, &cfg).unwrap()).to_bits()
            );
        }
        // The 2000-01-01 12:00 UTC label uses the leap-table UTC/TAI/TT chain:
        // TT/UT1 still differ by the configured ΔT within one rounding.
        let (jd_tt, jd_ut1) =
            utc_to_jd_with_config(2000, 1, 1, 12, 0, 0.0, GREG, &table, &cfg).unwrap();
        let dt = deltat_with_config(jd_ut1, &cfg).unwrap();
        assert!(
            (jd_tt - jd_ut1 - dt).abs() <= ulp(jd_tt),
            "middle-chain TT/UT1 split exceeds one ulp of ΔT"
        );
        // The configured backend threads through every branch: the default
        // `smh2016` with no provider is an honest error, never a silent
        // substitution (both civil-time branches call Delta-T).
        assert_eq!(
            utc_to_jd_with_config(
                1960,
                1,
                1,
                0,
                0,
                0.0,
                GREG,
                &table,
                &DeltaTConfig::default()
            ),
            Err(Error::DeltaTModelUnavailable { model: "smh2016" })
        );
        assert_eq!(
            utc_to_jd_with_config(
                2000,
                1,
                1,
                0,
                0,
                0.0,
                GREG,
                &table,
                &DeltaTConfig::default()
            ),
            Err(Error::DeltaTModelUnavailable { model: "smh2016" })
        );
    }

    #[test]
    fn utc_round_trips_recover_civil_labels() {
        let table = packaged_table();
        let cfg = espenak_cfg();
        // Forward then inverse recovers the label. The a priori chain bound
        // is 1 ms: each leg rounds a handful of times at the ~40 us JD-ulp
        // scale, so a dozen roundings stay an order of magnitude below it.
        // Labels stay clear of midnight by more than the bound so the date
        // fields are exact; leap-60 labels are excluded (the inverse spells
        // the same instant as next-midnight 00:00:00 — documented).
        let labels = [
            (1850, 3, 15, 12, 0, 0.0, GREG),
            (1971, 6, 30, 23, 59, 59.0, GREG),
            (1980, 1, 1, 0, 0, 1.0, GREG),
            (2000, 1, 1, 12, 0, 0.0, GREG),
            (2024, 2, 29, 12, 30, 15.25, GREG),
            (2030, 6, 15, 6, 45, 30.5, JUL),
            (2040, 1, 1, 0, 0, 1.0, GREG),
            (3000, 7, 4, 18, 30, 15.25, GREG),
        ];
        for (y, m, d, h, mi, s, cal) in labels {
            let (jd_tt, jd_ut1) =
                utc_to_jd_with_config(y, m, d, h, mi, s, cal, &table, &cfg).unwrap();
            for (name, back) in [
                (
                    "jdet",
                    jdet_to_utc_with_config(jd_tt, cal, &table, &cfg).unwrap(),
                ),
                (
                    "jdut1",
                    jdut1_to_utc_with_config(jd_ut1, cal, &table, &cfg).unwrap(),
                ),
            ] {
                assert_eq!(
                    (back.year, back.month, back.day, back.hour, back.minute),
                    (y, m, d, h, mi),
                    "{name} round-trip date drift at ({y}, {m}, {d})"
                );
                assert!(
                    (back.second - s).abs() <= 1e-3,
                    "{name} round-trip second drift {} at ({y}, {m}, {d})",
                    (back.second - s).abs()
                );
            }
        }
    }

    #[test]
    fn inverse_outer_branches_spell_the_ut1_estimate() {
        let table = packaged_table();
        let cfg = espenak_cfg();
        // Pre-1972 and post-2035 inverses return the calendar spelling of
        // the (refined, for TT) UT1 estimate: rebuilding a JD from the
        // fields recovers the estimate within the 1 ms chain bound.
        let bound_days = 1e-3 / 86400.0;
        for (jd_et, jd_ut1, cal) in [
            (2300000.0, 2300000.0, GREG),
            (2440000.0, 2440000.0, JUL),
            (2500000.5, 2500000.5, GREG),
        ] {
            for (name, back) in [
                (
                    "jdet",
                    jdet_to_utc_with_config(jd_et, cal, &table, &cfg).unwrap(),
                ),
                (
                    "jdut1",
                    jdut1_to_utc_with_config(jd_ut1, cal, &table, &cfg).unwrap(),
                ),
            ] {
                let rebuilt = julday(
                    back.year,
                    back.month,
                    back.day,
                    back.hour as f64 + back.minute as f64 / 60.0 + back.second / 3600.0,
                    cal,
                )
                .unwrap();
                let est = if name == "jdet" {
                    let first = jd_et - deltat_with_config(jd_et, &cfg).unwrap();
                    jd_et - deltat_with_config(first, &cfg).unwrap()
                } else {
                    jd_ut1
                };
                assert!(
                    (rebuilt - est).abs() <= bound_days,
                    "{name} spelling drift {} at {est}",
                    (rebuilt - est).abs()
                );
            }
        }
    }

    #[test]
    fn reconstruction_snap_fires_and_carries_at_midnight() {
        let table = packaged_table();
        let cfg = espenak_cfg();
        // One midnight JD in 1960 (pre-1972 inverse path, no ΔT iteration
        // on the input): a probe 43 us before it must snap up to the round
        // midnight, a probe 432 us before it must not.
        let midnight = julday(1960, 6, 16, 0.0, GREG).unwrap();
        assert_eq!(revjul(midnight, GREG).unwrap().hour, 0.0);
        let snap = jdut1_to_utc_with_config(midnight - 5e-10, GREG, &table, &cfg).unwrap();
        assert_eq!(
            (
                snap.year,
                snap.month,
                snap.day,
                snap.hour,
                snap.minute,
                snap.second
            ),
            (1960, 6, 16, 0, 0, 0.0)
        );
        let nosnap = jdut1_to_utc_with_config(midnight - 5e-9, GREG, &table, &cfg).unwrap();
        assert_eq!(
            (
                nosnap.year,
                nosnap.month,
                nosnap.day,
                nosnap.hour,
                nosnap.minute
            ),
            (1960, 6, 15, 23, 59)
        );
        assert!(
            (nosnap.second - 59.999568).abs() <= 1e-3,
            "unsnapped second {}",
            nosnap.second
        );
    }

    #[test]
    fn inverse_hostile_inputs_are_typed_errors() {
        let table = packaged_table();
        let cfg = espenak_cfg();
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                jdet_to_utc_with_config(bad, GREG, &table, &cfg),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                jdut1_to_utc_with_config(bad, GREG, &table, &cfg),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
        // A non-finite label second flows into `julday`, which reports it.
        assert!(matches!(
            utc_to_jd_with_config(2020, 1, 1, 0, 0, f64::NAN, GREG, &table, &cfg),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    // TT/TDB conversion using USNO Circular 179 equation 2.6.
    //
    // The seven-term truncation of the Fairhead-Bretagnon 1990 series returns
    // TDB-minus-TT in seconds. The argument is Julian centuries from J2000 TT;
    // periodic terms use the published frequencies and phases. The rate helper
    // analytically differentiates the same series. Scalar and two-part variants
    // share the formula, with two-part subtraction retaining sub-JD-ULP offsets.
    // Nonfinite inputs return NonFiniteArgument.

    #[test]
    fn tdb_series_epoch_and_magnitude_bound() {
        assert_eq!(TT_TDB_REFERENCE_JD, 2451545.0);
        // J2000 regression pin: the exact value agreed with the pinned
        // the model path in the sweep (TDB-TT ~= -95.7 us there;
        // the dominant annual term contributes ~= -71 us, the 4.2970-rad
        // phase term ~= -20 us). Guards transcription regressions.
        assert_eq!(
            tdb_minus_tt_seconds(2451545.0, 0.0).unwrap().to_bits(),
            0xbf19_1a2c_344d_6f44
        );
        assert_eq!(
            tt_to_tdb_fraction(2451545.0, 0.0).unwrap().to_bits(),
            0xbe13_0a5e_f843_d67c
        );
        // 1600-01-01 to 2200-12-31 monthly midpoints: the offset never
        // exceeds the coefficient-moduli bound (the ~2 ms the model claim).
        let mut jd = 2305507.5f64;
        while jd < 2524593.5 {
            let d = tdb_minus_tt_seconds(jd, 0.0).unwrap();
            assert!(d.abs() <= 0.00173, "TDB-TT out of bound at {jd}: {d} s");
            jd += 30.44;
        }
        // The series is odd-symmetric about no particular epoch, but it
        // must be smooth: adjacent-day steps stay far below the bound.
        let a = tdb_minus_tt_seconds(2451545.0, 0.0).unwrap();
        let b = tdb_minus_tt_seconds(2451546.0, 0.0).unwrap();
        assert!((a - b).abs() < 1e-4, "day step {a} -> {b}");
    }

    #[test]
    fn tdb_scalar_round_trip_within_one_ulp() {
        // Forward/backward across eras, seams and extremes; the fixed
        // point converges to far better than 1 ps, so one ulp of the
        // ~2.4 e 6 JD value is the tight honest bound.
        for &jd in &[
            0.5,
            100000.25,
            1500000.0,
            2305507.5,
            2400000.5,
            2433291.375,
            2451545.0,
            2451545.5,
            2460000.5,
            2500000.0,
            2800000.75,
            3000000.5,
        ] {
            // Per-point coefficient-moduli bound: periodic terms sum to at
            // most 0.001705 s and the mixed term to 0.000010*|T| s.
            let centuries = ((jd - TT_TDB_REFERENCE_JD) / DAYS_PER_JULIAN_CENTURY).abs();
            let bound_days = (0.001705 + 0.000010 * centuries) / 86400.0;
            let fwd = tt_to_tdb_jd(jd).unwrap();
            assert!(
                (fwd - jd).abs() <= bound_days,
                "forward offset out of bound at {jd}"
            );
            let back = tdb_to_tt_jd(fwd).unwrap();
            let ulp = (jd.next_up() - jd).max(jd - jd.next_down());
            assert!(
                (back - jd).abs() <= ulp,
                "round trip {jd} -> {fwd} -> {back} exceeds 1 ulp"
            );
            // Same from the TDB side.
            let there = tdb_to_tt_jd(jd).unwrap();
            let back_again = tt_to_tdb_jd(there).unwrap();
            assert!(
                (back_again - jd).abs() <= ulp,
                "reverse round trip {jd} exceeds 1 ulp"
            );
        }
    }

    #[test]
    fn tdb_fraction_matches_split_scalar_path() {
        // The adapter-compared kernel collapses onto the scalar path to 1
        // ulp of the full-JD value: `whole + fraction_result` vs
        // `tt_to_tdb_jd(whole + fraction)` differ only in the order of the
        // final additions. (Comparing fraction values directly would be
        // dishonest: the collapse itself quantizes the offset to JD
        // resolution.)
        for &(whole, frac) in &[
            (2451545.0, 0.0),
            (2451545.0, 0.5),
            (2451544.0, 0.999999),
            (2305507.0, 0.25),
            (2524593.0, 0.75),
            (100000.0, 0.1),
        ] {
            let via_frac = whole + tt_to_tdb_fraction(whole, frac).unwrap();
            let via_scalar = tt_to_tdb_jd(whole + frac).unwrap();
            let ulp = (via_scalar.next_up() - via_scalar).max(via_scalar - via_scalar.next_down());
            assert!(
                (via_frac - via_scalar).abs() <= ulp,
                "fraction path {whole}+{frac}: {via_frac} vs {via_scalar}"
            );
        }
    }

    #[test]
    fn tdb_two_part_collapse_within_one_ulp() {
        // Split epochs, including tailed ones whose offset lives in `lo`:
        // collapse stays within 1 ulp of the scalar value, and the applied
        // offset is recovered by `diff` within the 2^-75*|jd|
        // identity (same error-free machinery).
        for &jd in &[
            0.5,
            2451545.0,
            2451545.0 + 1e-10,
            2460000.5,
            2305507.5,
            2524593.5,
        ] {
            let epoch = TwoPartEpoch::split(jd).unwrap();
            let tdb = tt_to_tdb_two_part(epoch).unwrap();
            let scalar = tt_to_tdb_jd(jd).unwrap();
            let ulp = (scalar.next_up() - scalar).max(scalar - scalar.next_down());
            assert!(
                (tdb.value() - scalar).abs() <= ulp,
                "two-part collapse {jd}: {} vs {scalar}",
                tdb.value()
            );
            let recovered = tdb.diff(epoch);
            let applied = tdb_minus_tt_seconds(jd, 0.0).unwrap() / 86400.0;
            assert!(
                (recovered - applied).abs() <= 2f64.powi(-75) * jd.abs().max(1.0),
                "offset recovery at {jd}: {recovered} vs {applied}"
            );
            let tt = tdb_to_tt_two_part(epoch).unwrap();
            let scalar_back = tdb_to_tt_jd(jd).unwrap();
            let ulp_back =
                (scalar_back.next_up() - scalar_back).max(scalar_back - scalar_back.next_down());
            assert!(
                (tt.value() - scalar_back).abs() <= ulp_back,
                "two-part TT collapse {jd}"
            );
        }
    }

    #[test]
    fn tdb_rate_bound_and_finite_difference() {
        // A priori coefficient bound |rate - 1| <= 3.5 e-10 over |T| <= 4.1
        // (derivation in the comment above the magnitude test).
        let mut jd = 2305507.5f64;
        while jd < 2524593.5 {
            let rate = tdb_rate(jd).unwrap();
            assert!(
                (rate - 1.0).abs() <= 3.5e-10,
                "rate out of bound at {jd}: {rate}"
            );
            jd += 30.44;
        }
        // Central difference of the fraction kernel (h = 1 e-4 d): rounding
        // contributes ~1 e-12, truncation ~1 e-22, so 1 e-11 has an order of
        // magnitude of margin.
        for &jd in &[2451545.0, 2305507.5, 2524593.5, 2460000.0] {
            let h = 1e-4f64;
            let plus = tt_to_tdb_fraction(jd, 0.5 + h).unwrap();
            let minus = tt_to_tdb_fraction(jd, 0.5 - h).unwrap();
            let numeric = (plus - minus) / (2.0 * h);
            let analytic = tdb_rate(jd + 0.5).unwrap();
            assert!(
                (numeric - analytic).abs() <= 1e-11,
                "rate mismatch at {jd}: {numeric} vs {analytic}"
            );
        }
    }

    #[test]
    fn source_tt_windows_keep_exact_inward_representable_boundaries() {
        let delta = |tt: f64, bound: f64| {
            let whole = tt.floor();
            (whole - bound) + tt_to_tdb_fraction(whole, tt - whole).unwrap()
        };
        for (low, high) in [
            (2415020.5, 2488069.5),
            (2451545.0, 2451546.0),
            (-3100015.5, 8000016.5),
        ] {
            let (start, end) = tt_interval_from_tdb(low, high).unwrap();
            assert!(delta(start, low) >= 0.0);
            assert!(delta(start.next_down(), low) < 0.0);
            assert!(delta(end, high) <= 0.0);
            assert!(delta(end.next_up(), high) > 0.0);
        }
        for (low, high) in [(2.0, 1.0), (f64::NAN, 2.0), (1.0, f64::INFINITY)] {
            assert!(tt_interval_from_tdb(low, high).is_err());
        }
    }

    #[test]
    fn tt_geometric_source_uses_split_clock_and_physical_velocity_chain_rule() {
        let linear = |whole: f64, fraction: f64| {
            Ok(([(whole - 2451545.0) + fraction, 0.0, 0.0], [1.0, 0.0, 0.0]))
        };
        for tt in [2415108.6333333333, 2451545.125, 2460000.75] {
            let h = 0.0625;
            let rows = sample_geometric_tt(&[tt, tt - h, tt + h, tt], linear).unwrap();
            assert_eq!(rows[0], rows[3]);
            let finite_difference = (rows[2][0] - rows[1][0]) / (2.0 * h);
            assert!((finite_difference - rows[0][3]).abs() < 5e-11);
            let tdb_offset = rows[0][0] - (tt - 2451545.0);
            assert!((tdb_offset * 86400.0 - tdb_minus_tt_seconds(tt, 0.0).unwrap()).abs() < 5e-7);
        }
        assert!(sample_geometric_tt(&[], linear).unwrap().is_empty());
        assert!(matches!(
            sample_geometric_tt(&[f64::NAN], linear),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            sample_geometric_tt(&[2451545.0], |_, _| Err(Error::UnknownBody {
                body_id: 991
            })),
            Err(Error::UnknownBody { body_id: 991 })
        ));
    }

    // -- sidereal time / equation of the equinoxes ---------------

    /// Espenak-Meeus Delta-T values baked from the pinned the model (ephemeral
    /// probe, `set_delta_t_model("espenak_meeus")`, IERS off): J2000 and the
    /// two window boundaries, plus the model continuity offsets there.
    /// They pin the sidereal-time kernels against the model convention; the
    /// sweep in `time.json` covers the wider grid.
    fn sidereal_test_words() -> SiderealDeltaT {
        SiderealDeltaT::try_new(
            0.0007392804723490653,
            0.0003809980328242091,
            0.0010661391453652413,
        )
        .expect("finite test words")
    }

    #[test]
    fn sidereal_constants_match_release_identities() {
        assert_eq!(SIDEREAL_J2000_JD, 2451545.0);
        assert_eq!(SIDEREAL_MODERN_START_JD, 2396758.5);
        assert_eq!(SIDEREAL_MODERN_END_JD, 2469807.5);
        assert_eq!(
            SUN_EARTH_LIGHT_TIME_DAYS,
            1.495978707e11 / 299792458.0 / 86400.0
        );
        assert!((SUN_EARTH_LIGHT_TIME_DAYS - 0.0057755183).abs() < 1e-9);
        assert_eq!(DAYS_PER_JULIAN_MILLENNIUM, 365250.0);
    }

    #[test]
    fn gmst_iau2006_j2000_pin() {
        // Ephemeral the model probe (`_gmst_iau 2006_deg(2451545.0,
        // 2451545.0008)`); the IAU 2006 branch is pure closed form.
        let gmst = gmst_iau2006_deg(2451545.0, 2451545.0008).unwrap();
        assert_eq!(gmst, 280.46062243254534);
    }

    #[test]
    #[allow(deprecated)]
    fn legacy_mean_sidereal_time_j2000_pin() {
        // Ephemeral the model probe (`legacy_sidereal::mean_sidereal_time_deg(2451545.0)`
        // under the pinned Espenak-Meeus state).
        let dt = sidereal_test_words();
        assert_eq!(
            legacy_sidereal::mean_sidereal_time_deg(2451545.0, &dt).unwrap(),
            280.4606224304155
        );
    }

    #[test]
    #[allow(deprecated)]
    fn legacy_modern_branch_selected_strictly_inside_window() {
        let dt = sidereal_test_words();
        // Strictly inside: identical to the pure IAU 2006 branch.
        for jd in [
            SIDEREAL_MODERN_START_JD + 0.5,
            2451545.0,
            SIDEREAL_MODERN_END_JD - 0.5,
        ] {
            let jd_tt = jd + dt.delta_t_days();
            let expect = gmst_iau2006_deg(jd, jd_tt).unwrap().rem_euclid(360.0);
            assert_eq!(
                legacy_sidereal::mean_sidereal_time_deg(jd, &dt).unwrap(),
                expect
            );
        }
        // Exactly on the boundaries: the long-term branch with offsets.
        for (jd, off) in [
            (SIDEREAL_MODERN_START_JD, dt.offset_lo_deg()),
            (SIDEREAL_MODERN_END_JD, dt.offset_hi_deg()),
        ] {
            let expect = (mean_sidereal_longterm_deg(jd, dt.delta_t_days()).unwrap() - off)
                .rem_euclid(360.0);
            assert_eq!(
                legacy_sidereal::mean_sidereal_time_deg(jd, &dt).unwrap(),
                expect
            );
        }
    }

    #[test]
    #[allow(deprecated)]
    fn legacy_continuity_offsets_join_branches_smoothly() {
        // Boundary Delta-T values from the ephemeral the model probe (pinned
        // Espenak-Meeus state); the recomputed offsets must reproduce the
        // values, and each boundary value must join
        // the modern branch there (the offsets are ~1 e-3 deg, so the
        // (-180, 180] reduction is the identity and the join is exact up
        // to floating-point rounding).
        let dt = sidereal_test_words();
        let (off_lo, off_hi) =
            sidereal_continuity_offsets(8.231631738266498e-05, 0.0010773702417695439).unwrap();
        assert!((off_lo - dt.offset_lo_deg()).abs() <= 1e-12);
        assert!((off_hi - dt.offset_hi_deg()).abs() <= 1e-12);
        for (jd, dt_bdy) in [
            (SIDEREAL_MODERN_START_JD, 8.231631738266498e-05),
            (SIDEREAL_MODERN_END_JD, 0.0010773702417695439),
        ] {
            let jd_tt = jd + dt_bdy;
            let modern = gmst_iau2006_deg(jd, jd_tt).unwrap().rem_euclid(360.0);
            // `jd` itself takes the long-term branch (strict inequality),
            // evaluated with the boundary Delta-T value like the model.
            let dt_bdy_struct =
                SiderealDeltaT::try_new(dt_bdy, dt.offset_lo_deg(), dt.offset_hi_deg()).unwrap();
            let joined = legacy_sidereal::mean_sidereal_time_deg(jd, &dt_bdy_struct).unwrap();
            assert!(
                (joined - modern).abs() <= 1e-9,
                "branch join at {jd}: {joined} vs {modern}"
            );
        }
    }

    #[test]
    fn apparent_time_is_mean_plus_equation_of_equinoxes() {
        let dt = sidereal_test_words();
        for jd in [2396758.0, 2451545.0, 2469808.0, 2816787.5] {
            for (dpsi, eps) in [(0.0, 23.44), (0.0042, 23.4393), (-0.0031, 23.5)] {
                let mean = mean_sidereal_time_deg(jd, &dt).unwrap();
                let eoe = equation_of_equinoxes_deg(dpsi, eps).unwrap();
                let expect = (mean + eoe + 12.5).rem_euclid(360.0);
                assert_eq!(
                    apparent_sidereal_time_deg(jd, 12.5, &dt, dpsi, eps).unwrap(),
                    expect
                );
            }
        }
        // Zero nutation is the identity arm.
        assert_eq!(equation_of_equinoxes_deg(0.0, 23.44).unwrap(), 0.0);
        // Sign arm: positive dpsi advances apparent time.
        assert!(equation_of_equinoxes_deg(0.01, 23.44).unwrap() > 0.0);
    }

    #[test]
    fn sidtime_internal_j2000_pin_and_shape() {
        // Canonical VCW mean with the represented Espenak-Meeus TT offset.
        // Full-domain precision is checked independently against mpmath.
        let dt = sidereal_test_words();
        let gast = greenwich_apparent_sidereal_time_hours(2451545.0, 23.4392911, 0.0, &dt).unwrap();
        assert_eq!(gast, 18.697374828694365);
        // Greenwich is local time at longitude zero.
        assert_eq!(
            gast,
            local_apparent_sidereal_time_hours(2451545.0, 0.0, 23.4392911, 0.0, &dt).unwrap()
        );
        // Hours and degrees agree through the /15 scaling (up to the
        // radians round-trip, well inside 1 e-9 deg).
        let deg = apparent_sidereal_time_deg(2451545.0, 0.0, &dt, 0.0, 23.4392911).unwrap();
        assert!((gast * 15.0 - deg).abs() <= 1e-9);
        // Longitude wrap: +190 deg == -170 deg up to the last-ulp
        // association difference (`gast + 190/15` vs `gast - 170/15`
        // round the pre-modulo sum differently on both sides), and the
        // range is [0, 24).
        let a = local_apparent_sidereal_time_hours(2451545.0, 190.0, 23.44, 0.0, &dt).unwrap();
        let b = local_apparent_sidereal_time_hours(2451545.0, -170.0, 23.44, 0.0, &dt).unwrap();
        assert!((a - b).abs() <= 1e-14);
        assert!((0.0..24.0).contains(&a));
    }

    #[test]
    fn sibling_sidtime_composition_matches_release_words() {
        // The shared `_sidtime_internal` composition behind every sealed
        // GAST/LAST route (transit/star/lunar/kepler servings, solar-global
        // GAST, house ARMC): GMST from the caller Delta-T values plus the
        // equation of the equinoxes over caller nutation values. Pinned
        // against ephemeral the model `sidtime(jd)` values (v 3.2.1, ECL_NUT
        // true-obliquity/nutation-longitude inputs, Espenak-Meeus test
        // values for the model term). The 1 e-9 h ceiling (~0.05 mas) is the
        // sibling-composition agreement bound, not 1:1 acceptance: it guards
        // the untouched routes against speculative acquisition changes.
        let dt = sidereal_test_words();
        let cases: [(f64, f64, f64, f64); 3] = [
            (
                2451545.0,
                23.437676829667875,
                -0.0038699992874755788,
                18.69713811547336,
            ),
            (
                2460000.5,
                23.438418199581868,
                -0.0025762228160583125,
                10.306637826190714,
            ),
            (
                2460400.0,
                23.438685175576737,
                -0.0014548597103967106,
                0.5577813486623221,
            ),
        ];
        for (jd, eps_true, nut_lon, expect) in cases {
            let gast = greenwich_apparent_sidereal_time_hours(jd, eps_true, nut_lon, &dt).unwrap();
            assert!(
                (gast - expect).abs() <= 1e-9,
                "jd {jd}: gast {gast} vs release {expect}"
            );
            assert!((0.0..24.0).contains(&gast));
        }
    }

    #[test]
    fn house_armc_matches_apparent_at_zero_longitude() {
        let dt = sidereal_test_words();
        for jd in [2396758.0, 2451545.0, 2469808.0] {
            for (dpsi, deps) in [(0.0, 0.0), (0.0042, 0.0023), (-0.0031, -0.0041)] {
                let (armc0, eps_true) = house_armc_obliquity_deg(jd, &dt, dpsi, deps).unwrap();
                let jd_tt = jd + dt.delta_t_days();
                let expect_eps = vondrak_mean_obliquity_deg(jd_tt).unwrap() + deps;
                assert_eq!(eps_true, expect_eps);
                assert_eq!(
                    armc0,
                    apparent_sidereal_time_deg(jd, 0.0, &dt, dpsi, eps_true).unwrap()
                );
            }
        }
    }

    #[test]
    fn sidereal_hostile_inputs_are_typed_errors() {
        let dt = sidereal_test_words();
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                gmst_iau2006_deg(bad, 2451545.0),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                gmst_iau2006_deg(2451545.0, bad),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                mean_sidereal_longterm_deg(bad, 0.001),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                mean_sidereal_longterm_deg(2451545.0, bad),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                sidereal_continuity_offset_deg(bad, 0.001),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                mean_sidereal_time_deg(bad, &dt),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                equation_of_equinoxes_deg(bad, 23.44),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                equation_of_equinoxes_deg(0.001, bad),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                apparent_sidereal_time_deg(bad, 0.0, &dt, 0.0, 23.44),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                apparent_sidereal_time_deg(2451545.0, bad, &dt, 0.0, 23.44),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                local_apparent_sidereal_time_hours(bad, 0.0, 23.44, 0.0, &dt),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                local_apparent_sidereal_time_hours(2451545.0, 0.0, bad, 0.0, &dt),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                house_armc_obliquity_deg(bad, &dt, 0.0, 0.0),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                SiderealDeltaT::try_new(bad, 0.0, 0.0),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
    }

    #[test]
    fn tdb_hostile_inputs_are_typed_errors() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                tdb_minus_tt_seconds(bad, 0.0),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                tdb_minus_tt_seconds(2451545.0, bad),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                tt_to_tdb_fraction(bad, 0.5),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                tt_to_tdb_jd(bad),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                tdb_to_tt_jd(bad),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                tdb_rate(bad),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
        // Two-part epochs cannot hold non-finite values by construction.
        assert!(TwoPartEpoch::split(f64::NAN).is_err());
    }
}

// Equation of time and local apparent/mean solar time.
//
// Epochs are Julian days; geographic longitude is east-positive degrees.
// The equation of time is the normalized difference between apparent solar
// hour angle and uniform mean solar time, returned in days. Sun apparent RA
// and GAST come from supplied physical providers. LAT-to-LMT uses three
// fixed-point refinements; LMT-to-LAT evaluates once. Nonfinite input or
// provider values return NonFiniteArgument, and provider failures propagate.

/// Equation of Time in days.
///
/// `jd_ut` is UT1 Julian days. `sun_ra_deg_at` returns apparent Sun RA in
/// degrees and `gast_hours_at` Greenwich apparent sidereal hours; both receive
/// that same UT1 argument, Sun first. The Sun callback owns any required
/// UT1→TT conversion. These synchronous borrowed callbacks may perform I/O
/// or panic; callback errors propagate without replacement.
///
/// Same operation order as the model: `ra_hours = ra_deg / 15`,
/// `ut_hours = fmod(jd + 0.5, 1) * 24`,
/// `e_hours = fmod(gast - ra_hours + 12 - ut_hours, 24)` normalized into
/// `[-12, 12]` (exact `±12` is kept, not wrapped), returned as
/// `e_hours / 24`. Non-finite inputs or provider values are
/// [`Error::NonFiniteArgument`].
///
/// The sign convention is apparent minus mean solar time: LMT plus this
/// result gives LAT. Typical seasonal values are approximately −14 to +16
/// minutes; the actual value depends on the supplied solar/frame providers.
///
/// Wiring check with constant providers (J2000 noon, `ut_hours = 12`):
/// Sun RA `0` deg with GAST `18` h gives
/// `fmod(18 - 0 + 12 - 12, 24) = 18 > 12`, hence `-6` h:
///
/// ```rust
/// use libephemeris_core::time::time_equ_injected;
///
/// let eot = time_equ_injected(2451545.0, |_| Ok(0.0), |_| Ok(18.0)).unwrap();
/// assert_eq!(eot, -6.0 / 24.0);
/// ```
pub fn time_equ_injected(
    jd_ut: f64,
    sun_ra_deg_at: impl Fn(f64) -> Result<f64>,
    gast_hours_at: impl Fn(f64) -> Result<f64>,
) -> Result<f64> {
    if !jd_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut });
    }
    // the calculation order: Sun RA value first (`calc_ut(jd, 0, FLG_EQUATORIAL)`).
    let ra_hours = sun_ra_deg_at(jd_ut)? / 15.0;
    if !ra_hours.is_finite() {
        return Err(Error::NonFiniteArgument { value: ra_hours });
    }
    // Then the GAST value (`sidtime(jd)`).
    let gast = gast_hours_at(jd_ut)?;
    if !gast.is_finite() {
        return Err(Error::NonFiniteArgument { value: gast });
    }
    // Use the signed truncated remainder for the solar hour angle.
    let ut_hours = ((jd_ut + 0.5) % 1.0) * 24.0;
    let mut e_hours = (gast - ra_hours + 12.0 - ut_hours) % 24.0;
    if e_hours > 12.0 {
        e_hours -= 24.0;
    } else if e_hours < -12.0 {
        e_hours += 24.0;
    }
    Ok(e_hours / 24.0)
}

/// Local Apparent Time to Local Mean Time.
///
/// `jd_lat` is the Julian Day in LAT (sundial time), `longitude_deg` the
/// geographic longitude (East positive). `eot_days_at` is the injected
/// equation-of-time value. Same operation
/// order as the model: `offset = longitude / 360`, then three fixed-point
/// iterations of `ut_est = jd_lat - offset - eot(ut_est)` starting from
/// `ut_est = jd_lat - offset`, returning
/// `jd_lat - eot` using the final sampled EoT. Exactly three substitutions
/// are made; there is no convergence test or arbitrary-callback accuracy bound.
/// The EoT callback receives UT1 JD and returns days; errors and panics propagate.
/// Nonfinite inputs or provider values are [`Error::NonFiniteArgument`].
pub fn lat_to_lmt_injected(
    jd_lat: f64,
    longitude_deg: f64,
    eot_days_at: impl Fn(f64) -> Result<f64>,
) -> Result<f64> {
    if !jd_lat.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_lat });
    }
    if !longitude_deg.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: longitude_deg,
        });
    }
    let longitude_offset_days = longitude_deg / 360.0;
    let mut ut_est = jd_lat - longitude_offset_days;
    let mut eot = 0.0;
    for _ in 0..3 {
        eot = eot_days_at(ut_est)?;
        if !eot.is_finite() {
            return Err(Error::NonFiniteArgument { value: eot });
        }
        ut_est = jd_lat - longitude_offset_days - eot;
    }
    Ok(jd_lat - eot)
}

/// Local Mean Time to Local Apparent Time.
///
/// `jd_lmt` is the Julian Day in LMT, `longitude_deg` the geographic
/// longitude (East positive). `eot_days_at` is the injected
/// equation-of-time value, evaluated once at
/// `jd_lmt - longitude / 360`; returns `jd_lmt + eot`. Non-finite inputs
/// or provider values are [`Error::NonFiniteArgument`].
pub fn lmt_to_lat_injected(
    jd_lmt: f64,
    longitude_deg: f64,
    eot_days_at: impl Fn(f64) -> Result<f64>,
) -> Result<f64> {
    if !jd_lmt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_lmt });
    }
    if !longitude_deg.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: longitude_deg,
        });
    }
    let longitude_offset_days = longitude_deg / 360.0;
    let jd_ut_approx = jd_lmt - longitude_offset_days;
    let eot = eot_days_at(jd_ut_approx)?;
    if !eot.is_finite() {
        return Err(Error::NonFiniteArgument { value: eot });
    }
    Ok(jd_lmt + eot)
}

#[cfg(test)]
mod time11_tests {
    use super::*;

    /// Constant-value stub provider.
    fn const_word(c: f64) -> impl Fn(f64) -> Result<f64> {
        move |_| Ok(c)
    }

    #[test]
    fn time_equ_wiring_pin() {
        // J2000 noon: ut_hours = fmod(2451545.5, 1) * 24 = 12 exactly.
        // gast 18 h, RA 0 h -> fmod(18 - 0 + 12 - 12, 24) = 18 > 12 -> -6 h.
        let eot = time_equ_injected(2451545.0, const_word(0.0), const_word(18.0)).unwrap();
        assert_eq!(eot, -6.0 / 24.0);
    }

    #[test]
    fn time_equ_negative_arm() {
        // gast 0 h, RA 23 h (345 deg), UT12 h -> raw -23 h < -12 -> +24 = 1 h.
        let eot = time_equ_injected(2451545.0, const_word(345.0), const_word(0.0)).unwrap();
        assert!((eot - 1.0 / 24.0).abs() <= 1e-15);
    }

    #[test]
    fn time_equ_exact_half_day_boundaries_are_kept() {
        // Raw +12 h stays (guard is strict `>`): gast 12, RA 0, UT12.
        let plus = time_equ_injected(2451545.0, const_word(0.0), const_word(12.0)).unwrap();
        assert_eq!(plus, 0.5);
        // Raw -12 h stays (guard is strict `<`): gast 0, RA 12 h, UT12.
        let minus = time_equ_injected(2451545.0, const_word(180.0), const_word(0.0)).unwrap();
        assert_eq!(minus, -0.5);
    }

    #[test]
    fn time_equ_refusals() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                time_equ_injected(bad, const_word(0.0), const_word(18.0)),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                time_equ_injected(2451545.0, const_word(bad), const_word(18.0)),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                time_equ_injected(2451545.0, const_word(0.0), const_word(bad)),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
        // Provider errors propagate transparently.
        assert!(matches!(
            time_equ_injected(2451545.0, |_| Err(Error::IntegerOverflow), const_word(0.0)),
            Err(Error::IntegerOverflow)
        ));
    }

    #[test]
    fn lmt_to_lat_single_evaluation_wiring() {
        // The EoT value is read once at jd_lmt - longitude/360 (association), and LAT = LMT + EoT.
        let seen = std::cell::Cell::new(0.0f64);
        let probe = |ut: f64| {
            seen.set(ut);
            Ok(0.001)
        };
        let jd = 2451545.0;
        let lat = lmt_to_lat_injected(jd, 12.5, probe).unwrap();
        assert_eq!(seen.get(), jd - 12.5 / 360.0);
        assert!((lat - (jd + 0.001)).abs() <= 1e-15);
    }

    #[test]
    fn lat_lmt_zero_longitude_round_trip() {
        // Constant EoT stub: the three-iteration refinement is a fixed
        // point from the first pass, so the round trip closes to 1 ulp.
        let eot = const_word(0.005);
        for jd in [2451545.0, 2460311.0, 2378497.5, 2816787.5] {
            let lmt = lat_to_lmt_injected(jd, 0.0, &eot).unwrap();
            assert!((lmt - (jd - 0.005)).abs() <= 1e-12);
            let back = lmt_to_lat_injected(lmt, 0.0, &eot).unwrap();
            assert!((back - jd).abs() <= 1e-12);
        }
    }

    #[test]
    fn lat_to_lmt_iteration_matches_fixed_point() {
        // JD-periodic analytic stub exercises the refinement: the kernel
        // must agree with a 50-iteration fixed-point reference to 1 e-12.
        let stub = |ut: f64| Ok(0.002 * (ut.rem_euclid(1.0) * std::f64::consts::TAU).sin());
        let jd = 2460000.25;
        let lon = -73.99;
        let lmt = lat_to_lmt_injected(jd, lon, stub).unwrap();
        let off = lon / 360.0;
        // Exactly three refinement passes (the model count): the
        // kernel returns jd - eot_3, one evaluation past the twice-updated
        // UT, so an inline reference must agree to the last ulp.
        let mut ut2 = jd - off;
        for _ in 0..2 {
            ut2 = jd - off - stub(ut2).unwrap();
        }
        assert_eq!(lmt, jd - stub(ut2).unwrap());
        // Three passes already sit on the fixed point: a 50-iteration
        // reference agrees to 1 e-9 days.
        let mut ut = jd - off;
        for _ in 0..50 {
            ut = jd - off - stub(ut).unwrap();
        }
        assert!(((jd - lmt) - stub(ut).unwrap()).abs() <= 1e-9);
        // And lmt_to_lat inverts it to the iteration residual (~1 e-9 here).
        let back = lmt_to_lat_injected(lmt, lon, stub).unwrap();
        assert!((back - jd).abs() <= 1e-9);
    }

    #[test]
    fn longitude_seam_180() {
        // +/-180 deg (and 360 vs 0) shift the EoT lookup by whole days;
        // with a 1-day-periodic stub the conversions agree to 1 e-12 days.
        let stub = |ut: f64| Ok(0.003 * (ut.rem_euclid(1.0) * std::f64::consts::TAU).sin());
        for jd in [2451545.0, 2460311.7] {
            let a = lat_to_lmt_injected(jd, 180.0, stub).unwrap();
            let b = lat_to_lmt_injected(jd, -180.0, stub).unwrap();
            assert!((a - b).abs() <= 1e-12);
            let c = lmt_to_lat_injected(jd, 360.0, stub).unwrap();
            let d = lmt_to_lat_injected(jd, 0.0, stub).unwrap();
            assert!((c - d).abs() <= 1e-12);
        }
    }

    #[test]
    fn lat_lmt_nonfinite_refusals() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                lat_to_lmt_injected(bad, 0.0, const_word(0.0)),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                lat_to_lmt_injected(2451545.0, bad, const_word(0.0)),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                lat_to_lmt_injected(2451545.0, 0.0, const_word(bad)),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                lmt_to_lat_injected(bad, 0.0, const_word(0.0)),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                lmt_to_lat_injected(2451545.0, bad, const_word(0.0)),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                lmt_to_lat_injected(2451545.0, 0.0, const_word(bad)),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
    }
}

// Provider-bound equation-of-time and local solar-time conversion.
//
// The Delta-T provider returns TT-minus-UT1 days at UT1, the solar provider
// returns tropical apparent Sun RA in degrees at TT, and the GAST provider
// returns apparent sidereal hours at UT1. These entries compose the pure
// kernels above. They do not accept alternate bodies or coordinate flags.
// Nonfinite values and provider errors propagate without a substituted model.

/// Sealed Equation of Time in days.
///
/// `jd_ut` is UT1 Julian days. `delta_t_days_at` receives UT1 and returns
/// TT−UT1 days; `sun_ra_deg_at_tt` receives TT JD and returns apparent Sun RA
/// in degrees; `gast_hours_at_ut` receives the original UT1 and returns GAST
/// hours. The Delta-T callback runs first, then Sun, then GAST. Callbacks can
/// perform I/O or panic; their typed errors propagate. Arithmetic and sign
/// convention follow [`time_equ_injected`].
/// Non-finite epochs/values (including an overflowed `jd_ut + delta_t`) are
/// [`Error::NonFiniteArgument`]; provider errors propagate transparently.
pub fn time_equ_sealed(
    jd_ut: f64,
    delta_t_days_at: impl Fn(f64) -> Result<f64>,
    sun_ra_deg_at_tt: impl Fn(f64) -> Result<f64>,
    gast_hours_at_ut: impl Fn(f64) -> Result<f64>,
) -> Result<f64> {
    if !jd_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut });
    }
    let delta_t = delta_t_days_at(jd_ut)?;
    if !delta_t.is_finite() {
        return Err(Error::NonFiniteArgument { value: delta_t });
    }
    let jd_tt = jd_ut + delta_t;
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    // the calculation order: Sun RA value first, then the GAST value.
    let ra_deg = sun_ra_deg_at_tt(jd_tt)?;
    let gast = gast_hours_at_ut(jd_ut)?;
    time_equ_injected(jd_ut, |_| Ok(ra_deg), |_| Ok(gast))
}

/// Sealed Local Apparent Time to Local Mean Time.
///
/// `jd_lat` is the Julian Day in LAT, `longitude_deg` the geographic
/// longitude (East positive). The EoT value threads [`time_equ_sealed`]
/// (UT-to-TT plus Sun/GAST providers as above) through the three-iteration
/// refinement of [`lat_to_lmt_injected`] with the same operation order.
/// Non-finite inputs/values are [`Error::NonFiniteArgument`]; provider errors propagate.
pub fn lat_to_lmt_sealed(
    jd_lat: f64,
    longitude_deg: f64,
    delta_t_days_at: impl Fn(f64) -> Result<f64>,
    sun_ra_deg_at_tt: impl Fn(f64) -> Result<f64>,
    gast_hours_at_ut: impl Fn(f64) -> Result<f64>,
) -> Result<f64> {
    lat_to_lmt_injected(jd_lat, longitude_deg, |ut| {
        time_equ_sealed(ut, &delta_t_days_at, &sun_ra_deg_at_tt, &gast_hours_at_ut)
    })
}

/// Sealed Local Mean Time to Local Apparent Time.
///
/// `jd_lmt` is the Julian Day in LMT, `longitude_deg` the geographic
/// longitude (East positive). The EoT value threads [`time_equ_sealed`]
/// through the single evaluation of [`lmt_to_lat_injected`] with the same
/// operation order. Non-finite inputs/values are [`Error::NonFiniteArgument`]; provider
/// errors propagate.
pub fn lmt_to_lat_sealed(
    jd_lmt: f64,
    longitude_deg: f64,
    delta_t_days_at: impl Fn(f64) -> Result<f64>,
    sun_ra_deg_at_tt: impl Fn(f64) -> Result<f64>,
    gast_hours_at_ut: impl Fn(f64) -> Result<f64>,
) -> Result<f64> {
    lmt_to_lat_injected(jd_lmt, longitude_deg, |ut| {
        time_equ_sealed(ut, &delta_t_days_at, &sun_ra_deg_at_tt, &gast_hours_at_ut)
    })
}

#[cfg(test)]
mod time11_sealed_tests {
    use super::*;

    /// Constant-value stub provider.
    fn const_word(c: f64) -> impl Fn(f64) -> Result<f64> {
        move |_| Ok(c)
    }

    /// Physically plausible stub values: GAST advancing at the sidereal
    /// rate, Sun RA at the mean solar rate (rates chosen so the linear
    /// parts cancel in `GAST - RA - UT`, leaving only the millihour
    /// 1-day-periodic wobbles). A constant-value stub would give the EoT a
    /// slope of -1 in UT, under which the three-iteration fixed point
    /// oscillates instead of converging; these stubs keep
    /// `|dEoT/dUT| ~ 0.01`, the converged regime the model iterates in.
    const STUB_T0: f64 = 2451545.0;
    const STUB_GAST_RATE_HR_DAY: f64 = 24.06570982441908;
    const STUB_SUN_RATE_DEG_DAY: f64 = (24.06570982441908 - 24.0) * 15.0;
    fn stub_sun_ra_deg(tt: f64) -> Result<f64> {
        Ok(280.46
            + STUB_SUN_RATE_DEG_DAY * (tt - STUB_T0)
            + 0.5 * (tt.rem_euclid(1.0) * std::f64::consts::TAU).sin())
    }
    fn stub_gast_hours(ut: f64) -> Result<f64> {
        Ok(18.69
            + STUB_GAST_RATE_HR_DAY * (ut - STUB_T0)
            + 0.01 * (ut.rem_euclid(1.0) * std::f64::consts::TAU).sin())
    }

    #[test]
    fn sealed_time_equ_applies_ut_to_tt_shift() {
        // Delta-T shifts only the Sun-value epoch (TT); the GAST value stays
        // at the UT argument.
        let seen_sun = std::cell::Cell::new(0.0f64);
        let seen_gast = std::cell::Cell::new(0.0f64);
        let jd = 2451545.0;
        let eot = time_equ_sealed(
            jd,
            const_word(0.001),
            |tt| {
                seen_sun.set(tt);
                Ok(0.0)
            },
            |ut| {
                seen_gast.set(ut);
                Ok(18.0)
            },
        )
        .unwrap();
        assert_eq!(seen_sun.get(), jd + 0.001);
        assert_eq!(seen_gast.get(), jd);
        // Same values through the slice-1 kernel agree numerically.
        assert_eq!(
            eot,
            time_equ_injected(jd, const_word(0.0), const_word(18.0)).unwrap()
        );
    }

    #[test]
    fn sealed_time_equ_matches_injected_at_zero_delta_t() {
        let jd = 2460000.25;
        let sealed =
            time_equ_sealed(jd, const_word(0.0), const_word(37.5), const_word(9.25)).unwrap();
        let injected = time_equ_injected(jd, const_word(37.5), const_word(9.25)).unwrap();
        assert_eq!(sealed, injected);
    }

    #[test]
    fn sealed_time_equ_reads_sun_before_gast() {
        // the model acquisition order: Sun value first, then GAST.
        let order = std::cell::RefCell::new(Vec::new());
        let eot = time_equ_sealed(
            2451545.0,
            const_word(0.0),
            |_| {
                order.borrow_mut().push("sun");
                Ok(0.0)
            },
            |_| {
                order.borrow_mut().push("gast");
                Ok(18.0)
            },
        )
        .unwrap();
        assert_eq!(*order.borrow(), vec!["sun", "gast"]);
        assert_eq!(eot, -6.0 / 24.0);
    }

    #[test]
    fn sealed_delta_t_wiring_and_errors() {
        // A varying Delta-T threads per lookup UT (each EoT evaluation
        // re-reads it); errors propagate transparently.
        let eot = time_equ_sealed(
            2451545.0,
            |ut| Ok((ut - 2451545.0) * 0.0 + 0.002),
            const_word(0.0),
            const_word(18.0),
        )
        .unwrap();
        assert_eq!(eot, -6.0 / 24.0);
        assert!(matches!(
            time_equ_sealed(
                2451545.0,
                |_| Err(Error::IntegerOverflow),
                const_word(0.0),
                const_word(18.0)
            ),
            Err(Error::IntegerOverflow)
        ));
        for bad_dt in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                time_equ_sealed(
                    2451545.0,
                    const_word(bad_dt),
                    const_word(0.0),
                    const_word(18.0)
                ),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
    }

    #[test]
    fn sealed_nonfinite_refusals() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                time_equ_sealed(bad, const_word(0.0), const_word(0.0), const_word(0.0)),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                time_equ_sealed(2451545.0, const_word(0.0), const_word(bad), const_word(0.0)),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                time_equ_sealed(2451545.0, const_word(0.0), const_word(0.0), const_word(bad)),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                lat_to_lmt_sealed(bad, 0.0, const_word(0.0), const_word(0.0), const_word(0.0)),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                lat_to_lmt_sealed(
                    2451545.0,
                    bad,
                    const_word(0.0),
                    const_word(0.0),
                    const_word(0.0)
                ),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                lmt_to_lat_sealed(bad, 0.0, const_word(0.0), const_word(0.0), const_word(0.0)),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                lmt_to_lat_sealed(
                    2451545.0,
                    bad,
                    const_word(0.0),
                    const_word(0.0),
                    const_word(0.0)
                ),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
    }

    #[test]
    fn sealed_lat_lmt_round_trip_over_stubs() {
        // Physical-rate stubs (see helpers): the three-iteration
        // refinement sits on the fixed point, so the round trip closes.
        for jd in [2451545.0, 2460311.0, 2378497.5] {
            for lon in [0.0, 12.5, -73.99] {
                let lmt =
                    lat_to_lmt_sealed(jd, lon, const_word(0.001), stub_sun_ra_deg, stub_gast_hours)
                        .unwrap();
                let back = lmt_to_lat_sealed(
                    lmt,
                    lon,
                    const_word(0.001),
                    stub_sun_ra_deg,
                    stub_gast_hours,
                )
                .unwrap();
                assert!(
                    (back - jd).abs() <= 1e-9,
                    "jd={jd} lon={lon} residual={}",
                    (back - jd).abs()
                );
            }
        }
    }

    #[test]
    fn sealed_iterations_thread_through_delta_t_and_sun() {
        // lat_to_lmt performs three EoT evaluations (one Sun + one GAST
        // read each); lmt_to_lat performs one. Counts pin the orchestration
        // without duplicating the arithmetic.
        let sun_calls = std::cell::Cell::new(0);
        let gast_calls = std::cell::Cell::new(0);
        let sun = |_: f64| {
            sun_calls.set(sun_calls.get() + 1);
            Ok(10.0)
        };
        let gast = |_: f64| {
            gast_calls.set(gast_calls.get() + 1);
            Ok(20.0)
        };
        lat_to_lmt_sealed(2451545.0, 12.5, const_word(0.001), sun, gast).unwrap();
        assert_eq!(sun_calls.get(), 3);
        assert_eq!(gast_calls.get(), 3);
        sun_calls.set(0);
        gast_calls.set(0);
        lmt_to_lat_sealed(2451545.0, 12.5, const_word(0.001), sun, gast).unwrap();
        assert_eq!(sun_calls.get(), 1);
        assert_eq!(gast_calls.get(), 1);
    }

    #[test]
    fn sealed_seam_180() {
        // +/-180 deg shift the EoT lookup by whole days. The stub rates
        // cancel linearly and the wobbles are 1-day periodic, so both arms
        // read near-identical EoT values and agree to 1 e-9 days; the round
        // trip closes in the converged regime.
        for jd in [2451545.0, 2460311.7] {
            let a = lat_to_lmt_sealed(
                jd,
                180.0,
                const_word(0.0008),
                stub_sun_ra_deg,
                stub_gast_hours,
            )
            .unwrap();
            let b = lat_to_lmt_sealed(
                jd,
                -180.0,
                const_word(0.0008),
                stub_sun_ra_deg,
                stub_gast_hours,
            )
            .unwrap();
            assert!((a - b).abs() <= 1e-9, "seam {}", (a - b).abs());
            let back = lmt_to_lat_sealed(
                a,
                180.0,
                const_word(0.0008),
                stub_sun_ra_deg,
                stub_gast_hours,
            )
            .unwrap();
            assert!(
                (back - jd).abs() <= 1e-9,
                "round trip {}",
                (back - jd).abs()
            );
        }
    }
}

// Ephemeris-bound equation-of-time and local solar-time conversion.
//
// SealedTimeEquServing combines apparent solar RA from the bound position
// provider with sidereal time from the same configured Delta-T and nutation.
// All epochs are explicit UT1/TT Julian days; results use days and east-positive
// longitude degrees. Resources remain borrowed for the duration of the call.
// Data-range, frame and clock errors propagate through the shared kernels.

/// Concrete sealed Pipeline-A serving for the time-equation / LMT-LAT
/// family: serving-core Delta-T chain plus the sealed Sun-value and GAST
/// providers behind the sealed time-equation entries.
///
/// `serving` is the bound sealed LEB pair (geometric states +
/// TT-epoch nutation); `dt_cfg` is the serving-core Delta-T configuration
/// (the explicit Delta-T parameter contract; Espenak-Meeus is natively usable,
/// `smh2016` draws on the vendored table provider when the configuration
/// comes from a context, IERS still needs its caller-supplied provider).
/// Actual legacy offsets are resolved once for the published compatibility
/// getters; canonical GAST does not apply those values.
/// Fields are private with checked accessors; there is no global state
/// and no direct I/O by this adapter. Bound state, nutation and Delta-T
/// providers can perform I/O or panic; their own contracts still apply.
pub struct SealedTimeEquServing<'a, P> {
    serving: &'a crate::positions::SealedServing<'a, P>,
    dt_cfg: &'a DeltaTConfig<'a>,
    offset_lo_deg: f64,
    offset_hi_deg: f64,
}

impl<'a, P> SealedTimeEquServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Checked constructor: borrows the sealed serving pair and the
    /// Delta-T configuration, resolving the sidereal-continuity offsets
    /// from the chain's own modern-window boundary values. A chain with no
    /// usable backend reports [`Error::DeltaTModelUnavailable`] here
    /// (fail-closed at construction, never a silent model substitution).
    /// Resolving offsets invokes the bound Delta-T callbacks synchronously;
    /// they may perform I/O or panic. This adapter itself reads no global state.
    pub fn try_new(
        serving: &'a crate::positions::SealedServing<'a, P>,
        dt_cfg: &'a DeltaTConfig<'a>,
    ) -> Result<Self> {
        let dt_lo = deltat_with_config(SIDEREAL_MODERN_START_JD, dt_cfg)?;
        if !dt_lo.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_lo });
        }
        let dt_hi = deltat_with_config(SIDEREAL_MODERN_END_JD, dt_cfg)?;
        if !dt_hi.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_hi });
        }
        let (offset_lo_deg, offset_hi_deg) = sidereal_continuity_offsets(dt_lo, dt_hi)?;
        Ok(Self {
            serving,
            dt_cfg,
            offset_lo_deg,
            offset_hi_deg,
        })
    }

    /// Bound sealed serving pair.
    pub fn serving(&self) -> &crate::positions::SealedServing<'a, P> {
        self.serving
    }

    /// Bound serving-core Delta-T configuration.
    pub fn delta_t_config(&self) -> &DeltaTConfig<'a> {
        self.dt_cfg
    }

    /// Continuity offset at the 1850 boundary, in degrees.
    pub fn offset_lo_deg(&self) -> f64 {
        self.offset_lo_deg
    }

    /// Continuity offset at the 2050 boundary, in degrees.
    pub fn offset_hi_deg(&self) -> f64 {
        self.offset_hi_deg
    }

    /// TT−UT1 days at UT1 JD through the bound chain; nonfinite input/output
    /// returns [`Error::NonFiniteArgument`], and model errors propagate.
    pub fn delta_t_days_at_ut(&self, jd_ut: f64) -> Result<f64> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let dt = deltat_with_config(jd_ut, self.dt_cfg)?;
        if !dt.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt });
        }
        Ok(dt)
    }

    /// TT Julian date for UT1 JD through the bound chain.
    /// Nonfinite offsets or scalar sums are refused; no UTC conversion occurs.
    pub fn jd_tt_for_ut(&self, jd_ut: f64) -> Result<f64> {
        let dt = self.delta_t_days_at_ut(jd_ut)?;
        let jd_tt = jd_ut + dt;
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        Ok(jd_tt)
    }

    /// Sealed Sun apparent-RA value in degrees at a TT epoch. Sealed misses propagate
    /// fail-closed.
    pub fn sun_ra_deg_at_tt(&self, jd_tt: f64) -> Result<f64> {
        let (pos, _) = crate::positions::serving_compute_tt(
            self.serving,
            jd_tt,
            crate::constants::SUN,
            crate::constants::FLG_EQUATORIAL,
        )?;
        let ra_deg = pos[0];
        if !ra_deg.is_finite() {
            return Err(Error::NonFiniteArgument { value: ra_deg });
        }
        Ok(ra_deg)
    }

    /// Sealed GAST value in hours `[0, 24)` at a UT1 Julian date:
    /// the sidereal-time Greenwich composition over the sealed nutation pair at
    /// the TT epoch and the Vondrák mean obliquity the frame pipeline uses.
    /// The real construction-time legacy offsets remain available through
    /// their compatibility getters; they do not calibrate canonical GAST.
    pub fn gast_hours_at_ut(&self, jd_ut: f64) -> Result<f64> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let dt_days = self.delta_t_days_at_ut(jd_ut)?;
        let jd_tt = jd_ut + dt_days;
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let (dpsi_rad, deps_rad) = (self.serving.nutation_provider())(jd_tt)?;
        if !dpsi_rad.is_finite() {
            return Err(Error::NonFiniteArgument { value: dpsi_rad });
        }
        if !deps_rad.is_finite() {
            return Err(Error::NonFiniteArgument { value: deps_rad });
        }
        let dpsi_deg = dpsi_rad.to_degrees();
        let eps_mean_deg = crate::positions::vondrak_mean_obliquity_deg(jd_tt)?;
        let eps_true_deg = eps_mean_deg + deps_rad.to_degrees();
        let sdt = SiderealDeltaT::try_new(dt_days, self.offset_lo_deg, self.offset_hi_deg)?;
        greenwich_apparent_sidereal_time_hours(jd_ut, eps_true_deg, dpsi_deg, &sdt)
    }

    /// Apparent-minus-mean solar time in days at UT1 JD.
    /// Uses the bound Delta-T, Sun and GAST providers through [`time_equ_sealed`].
    pub fn time_equ(&self, jd_ut: f64) -> Result<f64> {
        time_equ_sealed(
            jd_ut,
            |ut| self.delta_t_days_at_ut(ut),
            |tt| self.sun_ra_deg_at_tt(tt),
            |ut| self.gast_hours_at_ut(ut),
        )
    }

    /// Convert local apparent JD to local mean JD at east-positive longitude
    /// in degrees. Uses three EoT substitutions, not a convergence test;
    /// borrowed-provider and numerical errors follow [`lat_to_lmt_sealed`].
    pub fn lat_to_lmt(&self, jd_lat: f64, longitude_deg: f64) -> Result<f64> {
        lat_to_lmt_sealed(
            jd_lat,
            longitude_deg,
            |ut| self.delta_t_days_at_ut(ut),
            |tt| self.sun_ra_deg_at_tt(tt),
            |ut| self.gast_hours_at_ut(ut),
        )
    }

    /// Convert local mean JD to local apparent JD at east-positive longitude
    /// in degrees. Uses one EoT evaluation; borrowed-provider and numerical
    /// errors follow [`lmt_to_lat_sealed`].
    pub fn lmt_to_lat(&self, jd_lmt: f64, longitude_deg: f64) -> Result<f64> {
        lmt_to_lat_sealed(
            jd_lmt,
            longitude_deg,
            |ut| self.delta_t_days_at_ut(ut),
            |tt| self.sun_ra_deg_at_tt(tt),
            |ut| self.gast_hours_at_ut(ut),
        )
    }
}

/// Equation of Time in days with the concrete serving-core Delta-T chain
///
/// `dt_cfg` is the serving-core Delta-T configuration (explicit Delta-T
/// parameter contract); `sun_ra_deg_at_tt` / `gast_hours_at_ut` are the
/// sealed Sun/GAST values. Thin orchestration over [`time_equ_sealed`]: each EoT
/// evaluation reads `deltat_with_config` at the lookup UT, the Sun value at
/// `+delta_t`, the GAST value at the UT argument, in the calculation order.
pub fn time_equ_pipeline_a(
    jd_ut: f64,
    dt_cfg: &DeltaTConfig<'_>,
    sun_ra_deg_at_tt: impl Fn(f64) -> Result<f64>,
    gast_hours_at_ut: impl Fn(f64) -> Result<f64>,
) -> Result<f64> {
    time_equ_sealed(
        jd_ut,
        |ut| deltat_with_config(ut, dt_cfg),
        sun_ra_deg_at_tt,
        gast_hours_at_ut,
    )
}

/// Local Apparent Time to Local Mean Time with the concrete serving-core
/// Delta-T chain.
///
/// Same composition as [`time_equ_pipeline_a`] threaded through the
/// three-iteration refinement of [`lat_to_lmt_sealed`] with the same
/// operation order.
pub fn lat_to_lmt_pipeline_a(
    jd_lat: f64,
    longitude_deg: f64,
    dt_cfg: &DeltaTConfig<'_>,
    sun_ra_deg_at_tt: impl Fn(f64) -> Result<f64>,
    gast_hours_at_ut: impl Fn(f64) -> Result<f64>,
) -> Result<f64> {
    lat_to_lmt_sealed(
        jd_lat,
        longitude_deg,
        |ut| deltat_with_config(ut, dt_cfg),
        sun_ra_deg_at_tt,
        gast_hours_at_ut,
    )
}

/// Local Mean Time to Local Apparent Time with the concrete serving-core
/// Delta-T chain.
///
/// Same composition as [`time_equ_pipeline_a`] threaded through the single
/// evaluation of [`lmt_to_lat_sealed`] with the same operation order.
pub fn lmt_to_lat_pipeline_a(
    jd_lmt: f64,
    longitude_deg: f64,
    dt_cfg: &DeltaTConfig<'_>,
    sun_ra_deg_at_tt: impl Fn(f64) -> Result<f64>,
    gast_hours_at_ut: impl Fn(f64) -> Result<f64>,
) -> Result<f64> {
    lmt_to_lat_sealed(
        jd_lmt,
        longitude_deg,
        |ut| deltat_with_config(ut, dt_cfg),
        sun_ra_deg_at_tt,
        gast_hours_at_ut,
    )
}

#[cfg(test)]
mod time11_pipeline_a_tests {
    use super::*;
    use crate::types::TwoPartEpoch;

    /// Provider that never serves anything: `has_body` is false so the
    /// sealed gate refuses before any state read. Enough for the
    /// binding-level refusal tests, which must fail without touching
    /// ephemeris data.
    struct Unserved;
    impl crate::positions::StateProvider for Unserved {
        fn provider_name(&self) -> &'static str {
            "time11-test/unserved"
        }
        fn has_body(&self, _body_id: i32) -> bool {
            false
        }
        fn body_coverage(&self, _body_id: i32) -> Option<(f64, f64)> {
            None
        }
        fn serving_coord_type(&self, _body_id: i32) -> Option<u32> {
            None
        }
        fn eval_raw(&self, body_id: i32, _jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
            Err(Error::LebBodyNotFound { body_id })
        }
        fn eval_raw_two_part(
            &self,
            body_id: i32,
            _epoch: TwoPartEpoch,
        ) -> Result<(Vec<f64>, Vec<f64>)> {
            Err(Error::LebBodyNotFound { body_id })
        }
        fn serving_entry_snapshot(&self, body_id: i32) -> Result<crate::positions::EntrySnapshot> {
            Err(Error::LebBodyNotFound { body_id })
        }
    }

    /// Natively usable Delta-T configuration (Espenak-Meeus, no
    /// providers): the serving-core chain in tests.
    fn espenak_cfg() -> DeltaTConfig<'static> {
        DeltaTConfig {
            model: DeltaTModel::EspenakMeeus,
            user_override_days: None,
            tidal_acceleration: TIDAL_DEFAULT,
            iers_enabled: false,
            iers_lookup_seconds: None,
            smh2016_lookup_seconds: None,
        }
    }

    /// Constant-value stub provider.
    fn const_word(c: f64) -> impl Fn(f64) -> Result<f64> {
        move |_| Ok(c)
    }

    /// Physically plausible stub values (slice-2 regime: linear rates
    /// cancel in `GAST - RA - UT`, leaving millihour 1-day-periodic
    /// wobbles so the three-iteration fixed point converges).
    const STUB_T0: f64 = 2451545.0;
    const STUB_GAST_RATE_HR_DAY: f64 = 24.06570982441908;
    const STUB_SUN_RATE_DEG_DAY: f64 = (24.06570982441908 - 24.0) * 15.0;
    fn stub_sun_ra_deg(tt: f64) -> Result<f64> {
        Ok(280.46
            + STUB_SUN_RATE_DEG_DAY * (tt - STUB_T0)
            + 0.5 * (tt.rem_euclid(1.0) * std::f64::consts::TAU).sin())
    }
    fn stub_gast_hours(ut: f64) -> Result<f64> {
        Ok(18.69
            + STUB_GAST_RATE_HR_DAY * (ut - STUB_T0)
            + 0.01 * (ut.rem_euclid(1.0) * std::f64::consts::TAU).sin())
    }

    #[test]
    fn pipeline_a_time_equ_success_matches_sealed() {
        // Concrete chain + stub sealed values: exactly identical to the
        // slice-2 entry over the same values and the same chain.
        let cfg = espenak_cfg();
        let jd = 2460000.25;
        let eot = time_equ_pipeline_a(jd, &cfg, const_word(37.5), const_word(9.25)).unwrap();
        let sealed = time_equ_sealed(
            jd,
            |ut| deltat_with_config(ut, &cfg),
            const_word(37.5),
            const_word(9.25),
        )
        .unwrap();
        assert_eq!(eot, sealed);
    }

    #[test]
    fn pipeline_a_epoch_rule_sun_at_tt_gast_at_ut() {
        // Distinctive user-override Delta-T pins the epoch rule: the Sun
        // value is read at jd_ut + delta_t, the GAST value at jd_ut.
        let cfg = DeltaTConfig {
            user_override_days: Some(0.007),
            ..espenak_cfg()
        };
        let jd = 2451545.0;
        let seen_sun = std::cell::Cell::new(0.0f64);
        let seen_gast = std::cell::Cell::new(0.0f64);
        let eot = time_equ_pipeline_a(
            jd,
            &cfg,
            |tt| {
                seen_sun.set(tt);
                Ok(0.0)
            },
            |ut| {
                seen_gast.set(ut);
                Ok(18.0)
            },
        )
        .unwrap();
        assert_eq!(seen_sun.get(), jd + 0.007);
        assert_eq!(seen_gast.get(), jd);
        assert_eq!(eot, -6.0 / 24.0);
        // lat_to_lmt threads the same rule through all three iterations:
        // each EoT evaluation reads Sun at UT + delta_t, GAST at UT.
        let pairs = std::cell::RefCell::new(Vec::new());
        lat_to_lmt_pipeline_a(
            jd,
            12.5,
            &cfg,
            |tt| {
                pairs.borrow_mut().push((tt, f64::NAN));
                Ok(0.0)
            },
            |ut| {
                let n = pairs.borrow().len();
                pairs.borrow_mut()[n - 1].1 = ut;
                Ok(18.0)
            },
        )
        .unwrap();
        let pairs = pairs.borrow();
        assert_eq!(pairs.len(), 3);
        for (tt, ut) in pairs.iter() {
            // Tolerance, not exactly: `(ut + 0.007) - ut` rounds at JD
            // magnitude (ulp ~ 5 e-10 there).
            assert!((tt - ut - 0.007).abs() <= 1e-9);
            // Lookups stay near the longitude-shifted epoch (the stub
            // EoT here is ~-0.25 d, so refinement legitimately walks).
            assert!((ut - (jd - 12.5 / 360.0)).abs() <= 0.5);
        }
    }

    #[test]
    fn pipeline_a_missing_delta_t_is_fail_closed() {
        // Bare default (smh 2016, no providers, IERS disabled) has no
        // usable backend: every entry reports DeltaTModelUnavailable,
        // never a substituted model.
        let cfg = DeltaTConfig::default();
        assert!(matches!(
            time_equ_pipeline_a(2451545.0, &cfg, const_word(0.0), const_word(18.0)),
            Err(Error::DeltaTModelUnavailable { .. })
        ));
        assert!(matches!(
            lat_to_lmt_pipeline_a(2451545.0, 0.0, &cfg, const_word(0.0), const_word(18.0)),
            Err(Error::DeltaTModelUnavailable { .. })
        ));
        assert!(matches!(
            lmt_to_lat_pipeline_a(2451545.0, 0.0, &cfg, const_word(0.0), const_word(18.0)),
            Err(Error::DeltaTModelUnavailable { .. })
        ));
        // Sealed-value errors propagate transparently with a usable chain.
        let good = espenak_cfg();
        assert!(matches!(
            time_equ_pipeline_a(
                2451545.0,
                &good,
                |_| Err(Error::IntegerOverflow),
                const_word(18.0),
            ),
            Err(Error::IntegerOverflow)
        ));
        assert!(matches!(
            time_equ_pipeline_a(2451545.0, &good, const_word(0.0), |_| Err(
                Error::IntegerOverflow
            ),),
            Err(Error::IntegerOverflow)
        ));
    }

    #[test]
    fn pipeline_a_lat_lmt_round_trip_over_physical_stubs() {
        // Physical-rate stubs + concrete Espenak-Meeus chain: the
        // three-iteration refinement sits on the fixed point, so the
        // round trip closes.
        let cfg = espenak_cfg();
        for jd in [2451545.0, 2460311.0, 2378497.5] {
            for lon in [0.0, 12.5, -73.99] {
                let lmt =
                    lat_to_lmt_pipeline_a(jd, lon, &cfg, stub_sun_ra_deg, stub_gast_hours).unwrap();
                let back = lmt_to_lat_pipeline_a(lmt, lon, &cfg, stub_sun_ra_deg, stub_gast_hours)
                    .unwrap();
                assert!(
                    (back - jd).abs() <= 1e-9,
                    "jd={jd} lon={lon} residual={}",
                    (back - jd).abs()
                );
            }
        }
    }

    #[test]
    fn pipeline_a_nonfinite_refusals() {
        let cfg = espenak_cfg();
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                time_equ_pipeline_a(bad, &cfg, const_word(0.0), const_word(0.0)),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                lat_to_lmt_pipeline_a(bad, 0.0, &cfg, const_word(0.0), const_word(0.0)),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                lat_to_lmt_pipeline_a(2451545.0, bad, &cfg, const_word(0.0), const_word(0.0)),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                lmt_to_lat_pipeline_a(bad, 0.0, &cfg, const_word(0.0), const_word(0.0)),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                lmt_to_lat_pipeline_a(2451545.0, bad, &cfg, const_word(0.0), const_word(0.0)),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
    }

    #[test]
    fn serving_try_new_resolves_offsets_and_refuses_missing_chain() {
        let provider = Unserved;
        let nut = |_: f64| Ok((0.0, 0.0));
        let serving = crate::positions::SealedServing::new(&provider, &nut);
        // Usable chain: offsets resolve from the chain's own boundary
        // values (user override pins them deterministically here).
        let cfg = DeltaTConfig {
            user_override_days: Some(0.001),
            ..espenak_cfg()
        };
        let bound = SealedTimeEquServing::try_new(&serving, &cfg).unwrap();
        let (lo, hi) = sidereal_continuity_offsets(0.001, 0.001).unwrap();
        assert_eq!(bound.offset_lo_deg(), lo);
        assert_eq!(bound.offset_hi_deg(), hi);
        assert!((bound.jd_tt_for_ut(2460000.0).unwrap() - 2460000.001).abs() <= 1e-12);
        assert!(matches!(
            bound.jd_tt_for_ut(f64::NAN),
            Err(Error::NonFiniteArgument { .. })
        ));
        // Unusable chain fails closed at construction.
        let missing = DeltaTConfig::default();
        assert!(matches!(
            SealedTimeEquServing::try_new(&serving, &missing),
            Err(Error::DeltaTModelUnavailable { .. })
        ));
    }

    #[test]
    fn serving_gast_reads_nutation_at_tt() {
        // User-override chain + recording nutation stub: the nutation
        // pair is read at jd_ut + delta_t, and the value matches the
        // composition over those exact values.
        let provider = Unserved;
        let seen_nut = std::cell::Cell::new(0.0f64);
        let nut = |tt: f64| {
            seen_nut.set(tt);
            Ok((0.0002, -0.0001))
        };
        let serving = crate::positions::SealedServing::new(&provider, &nut);
        let dt = 0.0011;
        let cfg = DeltaTConfig {
            user_override_days: Some(dt),
            ..espenak_cfg()
        };
        let bound = SealedTimeEquServing::try_new(&serving, &cfg).unwrap();
        let jd = 2460000.25;
        let gast = bound.gast_hours_at_ut(jd).unwrap();
        assert_eq!(seen_nut.get(), jd + dt);
        let sdt =
            SiderealDeltaT::try_new(dt, bound.offset_lo_deg(), bound.offset_hi_deg()).unwrap();
        let expected = greenwich_apparent_sidereal_time_hours(
            jd,
            vondrak_mean_obliquity_deg(jd + dt).unwrap() + (-0.0001f64).to_degrees(),
            0.0002f64.to_degrees(),
            &sdt,
        )
        .unwrap();
        assert_eq!(gast, expected);
        // A failing nutation value propagates (fail-closed, never mean).
        let bad_nut = |_: f64| -> Result<(f64, f64)> { Err(Error::IntegerOverflow) };
        let bad_serving = crate::positions::SealedServing::new(&provider, &bad_nut);
        let bad = SealedTimeEquServing::try_new(&bad_serving, &cfg).unwrap();
        assert!(matches!(
            bad.gast_hours_at_ut(jd),
            Err(Error::IntegerOverflow)
        ));
    }

    #[test]
    fn serving_sun_word_is_fail_closed_without_coverage() {
        // Unserved channel: the sealed gate refuses before any state
        // read, and the entries propagate the miss fail-closed.
        let provider = Unserved;
        let nut = |_: f64| Ok((0.0, 0.0));
        let serving = crate::positions::SealedServing::new(&provider, &nut);
        let cfg = DeltaTConfig {
            user_override_days: Some(0.001),
            ..espenak_cfg()
        };
        let bound = SealedTimeEquServing::try_new(&serving, &cfg).unwrap();
        assert!(matches!(
            bound.sun_ra_deg_at_tt(2460000.0),
            Err(Error::LebBodyNotFound { .. })
        ));
        assert!(matches!(
            bound.sun_ra_deg_at_tt(f64::NAN),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            bound.time_equ(2460000.0),
            Err(Error::LebBodyNotFound { .. })
        ));
        assert!(matches!(
            bound.lat_to_lmt(2460000.0, 0.0),
            Err(Error::LebBodyNotFound { .. })
        ));
        assert!(matches!(
            bound.lmt_to_lat(2460000.0, 0.0),
            Err(Error::LebBodyNotFound { .. })
        ));
    }
}

#[cfg(test)]
mod astropy_core_time_tests {
    use super::*;

    #[test]
    fn astropy_core_time_pair_preserves_refined_model_observable() {
        // The existing seconds provider independently performs the two
        // substitutions and reads the curve at the settled scalar TT.
        // These branch inputs include the UT/month boundary used by the
        // source-input declaration; no scientific expected value is invented.
        for ut in [2451544.4995, 2451545.0, 2460400.25] {
            let (tt, delta) = default_timescale_from_ut1(ut).unwrap();
            assert_eq!(
                tt.to_bits(),
                default_timescale_tt_from_ut1(ut).unwrap().to_bits()
            );
            assert_eq!(
                delta.to_bits(),
                crate::smh2016::smh2016_lookup_seconds(ut)
                    .unwrap()
                    .to_bits()
            );
        }
    }

    #[test]
    fn astropy_core_time_pair_rejects_nonfinite_public_input() {
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                default_timescale_from_ut1(value),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
    }
}
