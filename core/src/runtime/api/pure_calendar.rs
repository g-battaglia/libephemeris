//! Primary calendar records over the independent, stateless time producers.
//! Binding-specific scalar and tuple contracts continue to use `crate::time`.

use super::{Calendar, CalendarDate, Error};

/// Julian Day produced by calendar arithmetic, without a time-scale assignment.
/// No status, diagnostic or source flag is supplied by this producer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct JulianDay {
    /// Day number and fractional day, retaining the producer's binary64 value.
    pub value: f64,
}

/// Complete civil-date validation and normalization result.
#[derive(Debug, Clone, PartialEq)]
pub struct DateConversion {
    /// True when the original year, month and day survive the round trip.
    pub valid: bool,
    /// Normalized Julian Day in days; no UTC, UT1 or TT interpretation is added.
    pub julian_day: f64,
    /// Calendar fields from the producer's inverse conversion.
    pub date: CalendarDate,
    /// Producer integer status when available; this engine supplies none.
    pub status: Option<i32>,
    /// Producer text when available; this calculation has no diagnostic channel.
    pub diagnostic: Option<String>,
}

/// Calendar weekday, with no status or diagnostic channel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Weekday {
    /// Monday-based index in `0..=6`.
    pub index: i32,
}

/// Civil clock fields returned by a pure timezone shift.
/// This label conversion does not acquire leap-second or Earth-rotation data.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UtcDateTime {
    /// Astronomical year, including year zero.
    pub year: i32,
    /// Calendar month.
    pub month: i32,
    /// Calendar day of month.
    pub day: i32,
    /// Hour of day.
    pub hour: i32,
    /// Minute of hour.
    pub minute: i32,
    /// Second of minute, including a retained leap-second label when supplied.
    pub second: f64,
}

/// Convert a date in an explicit proleptic calendar to a typed Julian Day.
///
/// Astronomical years and decimal hours feed the existing calendar recurrence.
/// Out-of-range month/day labels are normalized, not validated. This is pure
/// calendar arithmetic: it neither loads data nor determines a time scale.
///
/// # Errors
/// Non-finite hours retain the producer's typed error. No status or diagnostic
/// is fabricated; [`date_conversion`] additionally reports civil-date validity.
pub fn julday(
    year: i32,
    month: i32,
    day: i32,
    hour: f64,
    calendar: Calendar,
) -> Result<JulianDay, Error> {
    crate::time::julday(year, month, day, hour, calendar).map(|value| JulianDay { value })
}

/// Validate and normalize a date in an explicit proleptic calendar.
///
/// Returns the existing checked producer's complete result, including invalid
/// dates as `valid=false`. Fractional hours and normalized date fields are moved
/// without rounding or reconstruction. No integer status or text is produced.
///
/// # Errors
/// Non-finite hours and unrepresentable inverse-calendar components propagate
/// the time engine's errors unchanged. No session, model or ephemeris is read.
pub fn date_conversion(
    year: i32,
    month: i32,
    day: i32,
    hour: f64,
    calendar: Calendar,
) -> Result<DateConversion, Error> {
    let result = crate::time::convert_date_checked(year, month, day, hour, calendar)?;
    Ok(DateConversion {
        valid: result.valid,
        julian_day: result.jd,
        date: result.date,
        status: None,
        diagnostic: None,
    })
}

/// Return the calendar weekday of a Julian Day: Monday zero through Sunday six.
///
/// This pure calculation does not assign a time scale or consult session data.
/// Non-finite input propagates the time producer's error; finite values retain
/// its actual domain, including dates outside ordinary ephemeris coverage.
pub fn day_of_week(julian_day: f64) -> Result<Weekday, Error> {
    crate::time::day_of_week(julian_day).map(|index| Weekday { index })
}

/// Shift a civil local-time label to UTC using an offset in hours east of UTC.
///
/// Uses the time engine's established Gregorian date rollover, minute rounding
/// and leap-label handling without introducing a new numerical path. This does
/// not transform UTC to UT1 or TT and consumes no configured leap-second data.
///
/// # Errors
/// Non-finite seconds/offsets and unrepresentable calendar or offset components
/// retain their original typed errors. The returned fields are moved directly
/// from the numerical producer; there is no status or diagnostic channel.
#[allow(clippy::too_many_arguments)]
pub fn utc_time_zone(
    year: i32,
    month: i32,
    day: i32,
    hour: i32,
    minute: i32,
    second: f64,
    offset: f64,
) -> Result<UtcDateTime, Error> {
    let fields = crate::time::utc_time_zone(year, month, day, hour, minute, second, offset)?;
    Ok(UtcDateTime {
        year: fields.year,
        month: fields.month,
        day: fields.day,
        hour: fields.hour,
        minute: fields.minute,
        second: fields.second,
    })
}
