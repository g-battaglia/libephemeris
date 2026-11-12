//! Pure primary records forward the existing independent producers unchanged.
use super::operations::native;
use super::*;

/// Convert calendar fields to a Julian day without assigning a time scale.
///
/// # Arguments
///
/// * `year` - Astronomical year, including year zero.
/// * `month` - Calendar month; normalization follows the calendar producer.
/// * `day` - Calendar day of month; normalization follows the calendar producer.
/// * `hour` - Civil hour, decimal hours.
/// * `calendar` - Explicit Julian or Gregorian calendar.
///
/// # Returns
///
/// Complete [`JulianDay`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// Returns the native input, configuration, data-coverage or calculation
/// failure. Its structured cause is available through `Error::source()`.
#[allow(clippy::too_many_arguments)]
pub fn julday(
    year: i32,
    month: i32,
    day: i32,
    hour: f64,
    calendar: Calendar,
) -> Result<JulianDay, Error> {
    let value = native(
        "julday",
        crate::engine_api::julday(year, month, day, hour, calendar.native()),
    )?;
    Ok(JulianDay { value: value.value })
}

/// Convert a Julian day to calendar fields.
///
/// # Arguments
///
/// * `jd` - Unassigned-scale Julian day.
/// * `calendar` - Explicit Julian or Gregorian calendar.
///
/// # Returns
///
/// Complete [`CalendarDate`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// Returns the native input, configuration, data-coverage or calculation
/// failure. Its structured cause is available through `Error::source()`.
#[allow(clippy::too_many_arguments)]
pub fn revjul(jd: f64, calendar: Calendar) -> Result<CalendarDate, Error> {
    let value = native("revjul", crate::time::revjul(jd, calendar.native()))?;
    Ok(CalendarDate {
        year: value.year,
        month: value.month,
        day: value.day,
        hour: value.hour,
    })
}

/// Validate and normalize civil calendar fields by round trip.
///
/// # Arguments
///
/// * `year` - Astronomical year, including year zero.
/// * `month` - Calendar month; normalization follows the calendar producer.
/// * `day` - Calendar day of month; normalization follows the calendar producer.
/// * `hour` - Civil hour, decimal hours.
/// * `calendar` - Explicit Julian or Gregorian calendar.
///
/// # Returns
///
/// Complete [`DateConversion`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// Returns the native input, configuration, data-coverage or calculation
/// failure. Its structured cause is available through `Error::source()`.
#[allow(clippy::too_many_arguments)]
pub fn date_conversion(
    year: i32,
    month: i32,
    day: i32,
    hour: f64,
    calendar: Calendar,
) -> Result<DateConversion, Error> {
    let value = native(
        "date_conversion",
        crate::engine_api::date_conversion(year, month, day, hour, calendar.native()),
    )?;
    Ok(DateConversion {
        valid: value.valid,
        julian_day: value.julian_day,
        date: CalendarDate {
            year: value.date.year,
            month: value.date.month,
            day: value.date.day,
            hour: value.date.hour,
        },
        status: value.status,
        diagnostic: value.diagnostic,
    })
}

/// Return the Monday-based weekday for a Julian day.
///
/// # Arguments
///
/// * `jd` - Unassigned-scale Julian day.
///
/// # Returns
///
/// Complete [`Weekday`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// Returns the native input, configuration, data-coverage or calculation
/// failure. Its structured cause is available through `Error::source()`.
#[allow(clippy::too_many_arguments)]
pub fn day_of_week(jd: f64) -> Result<Weekday, Error> {
    let value = native("day_of_week", crate::engine_api::day_of_week(jd))?;
    Ok(Weekday { index: value.index })
}

/// Shift civil clock fields by a time-zone offset.
///
/// # Arguments
///
/// * `year` - Astronomical year, including year zero.
/// * `month` - Calendar month; normalization follows the calendar producer.
/// * `day` - Calendar day of month; normalization follows the calendar producer.
/// * `hour` - UTC hour in 0..=23.
/// * `minute` - UTC minute in 0..=59.
/// * `second` - Input scalar in degrees.
/// * `offset` - Time-zone offset in hours, positive east.
///
/// # Returns
///
/// Complete [`UtcDateTime`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// Returns the native input, configuration, data-coverage or calculation
/// failure. Its structured cause is available through `Error::source()`.
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
    let value = native(
        "utc_time_zone",
        crate::engine_api::utc_time_zone(year, month, day, hour, minute, second, offset),
    )?;
    Ok(UtcDateTime {
        year: value.year,
        month: value.month,
        day: value.day,
        hour: value.hour,
        minute: value.minute,
        second: value.second,
    })
}

/// Normalize an angle to one positive turn in degrees.
///
/// # Arguments
///
/// * `value` - Input scalar in degrees.
///
/// # Returns
///
/// Complete [`Angle`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// This producer has no error channel; its native floating-point/domain
/// behavior is retained.
#[allow(clippy::too_many_arguments)]
pub fn degnorm(value: f64) -> Result<Angle, Error> {
    Ok(Angle {
        value: crate::utils::degnorm(value),
    })
}

/// Normalize an angle to one positive turn in radians.
///
/// # Arguments
///
/// * `value` - Input scalar in radians.
///
/// # Returns
///
/// Complete [`Angle`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// This producer has no error channel; its native floating-point/domain
/// behavior is retained.
#[allow(clippy::too_many_arguments)]
pub fn radnorm(value: f64) -> Result<Angle, Error> {
    Ok(Angle {
        value: crate::utils::radnorm(value),
    })
}

/// Return the positive angular difference in degrees.
///
/// # Arguments
///
/// * `first` - Input scalar in degrees.
/// * `second` - Input scalar in degrees.
///
/// # Returns
///
/// Complete [`Angle`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// This producer has no error channel; its native floating-point/domain
/// behavior is retained.
#[allow(clippy::too_many_arguments)]
pub fn difdegn(first: f64, second: f64) -> Result<Angle, Error> {
    Ok(Angle {
        value: crate::utils::difdegn(first, second),
    })
}

/// Return the signed shortest angular difference in degrees.
///
/// # Arguments
///
/// * `first` - Input scalar in degrees.
/// * `second` - Input scalar in degrees.
///
/// # Returns
///
/// Complete [`Angle`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// This producer has no error channel; its native floating-point/domain
/// behavior is retained.
#[allow(clippy::too_many_arguments)]
pub fn difdeg2n(first: f64, second: f64) -> Result<Angle, Error> {
    Ok(Angle {
        value: crate::utils::difdeg2n(first, second),
    })
}

/// Return the signed shortest angular difference in radians.
///
/// # Arguments
///
/// * `first` - Input scalar in radians.
/// * `second` - Input scalar in radians.
///
/// # Returns
///
/// Complete [`Angle`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// This producer has no error channel; its native floating-point/domain
/// behavior is retained.
#[allow(clippy::too_many_arguments)]
pub fn difrad2n(first: f64, second: f64) -> Result<Angle, Error> {
    Ok(Angle {
        value: crate::utils::difrad2n(first, second),
    })
}

/// Return the circular midpoint in degrees.
///
/// # Arguments
///
/// * `first` - Input scalar in degrees.
/// * `second` - Input scalar in degrees.
///
/// # Returns
///
/// Complete [`Angle`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// This producer has no error channel; its native floating-point/domain
/// behavior is retained.
#[allow(clippy::too_many_arguments)]
pub fn deg_midp(first: f64, second: f64) -> Result<Angle, Error> {
    Ok(Angle {
        value: crate::utils::deg_midp(first, second),
    })
}

/// Return the circular midpoint in radians.
///
/// # Arguments
///
/// * `first` - Input scalar in radians.
/// * `second` - Input scalar in radians.
///
/// # Returns
///
/// Complete [`Angle`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// This producer has no error channel; its native floating-point/domain
/// behavior is retained.
#[allow(clippy::too_many_arguments)]
pub fn rad_midp(first: f64, second: f64) -> Result<Angle, Error> {
    Ok(Angle {
        value: crate::utils::rad_midp(first, second),
    })
}

/// Round a binary64 scalar to the native signed integer domain.
///
/// # Arguments
///
/// * `value` - Input scalar in the native rounding domain.
///
/// # Returns
///
/// Complete [`Integer`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// Returns the native input, configuration, data-coverage or calculation
/// failure. Its structured cause is available through `Error::source()`.
#[allow(clippy::too_many_arguments)]
pub fn d2l(value: f64) -> Result<Integer, Error> {
    Ok(Integer {
        value: native("d2l", crate::utils::d2l_checked(value))?,
    })
}

/// Normalize centiseconds to one positive turn.
///
/// # Arguments
///
/// * `value` - Input scalar in centiseconds.
///
/// # Returns
///
/// Complete [`Integer`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// This producer has no error channel; its native floating-point/domain
/// behavior is retained.
#[allow(clippy::too_many_arguments)]
pub fn csnorm(value: i64) -> Result<Integer, Error> {
    Ok(Integer {
        value: crate::utils::csnorm(value),
    })
}

/// Round centiseconds to an integral angular second.
///
/// # Arguments
///
/// * `value` - Input scalar in centiseconds.
///
/// # Returns
///
/// Complete [`Integer`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// This producer has no error channel; its native floating-point/domain
/// behavior is retained.
#[allow(clippy::too_many_arguments)]
pub fn csroundsec(value: i64) -> Result<Integer, Error> {
    Ok(Integer {
        value: crate::utils::csroundsec(value),
    })
}

/// Return the positive angular difference in centiseconds.
///
/// # Arguments
///
/// * `first` - Input scalar in centiseconds.
/// * `second` - Input scalar in centiseconds.
///
/// # Returns
///
/// Complete [`Integer`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// This producer has no error channel; its native floating-point/domain
/// behavior is retained.
#[allow(clippy::too_many_arguments)]
pub fn difcsn(first: i64, second: i64) -> Result<Integer, Error> {
    Ok(Integer {
        value: crate::utils::difcsn(first, second),
    })
}

/// Return the signed shortest angular difference in centiseconds.
///
/// # Arguments
///
/// * `first` - Input scalar in centiseconds.
/// * `second` - Input scalar in centiseconds.
///
/// # Returns
///
/// Complete [`Integer`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// This producer has no error channel; its native floating-point/domain
/// behavior is retained.
#[allow(clippy::too_many_arguments)]
pub fn difcs2n(first: i64, second: i64) -> Result<Integer, Error> {
    Ok(Integer {
        value: crate::utils::difcs2n(first, second),
    })
}

/// Format an angle supplied in centiseconds.
///
/// # Arguments
///
/// * `value` - Input scalar in centiseconds.
///
/// # Returns
///
/// Complete [`Text`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// This producer has no error channel; its native floating-point/domain
/// behavior is retained.
#[allow(clippy::too_many_arguments)]
pub fn cs2degstr(value: i64) -> Result<Text, Error> {
    Ok(Text {
        value: crate::utils::cs2degstr(value),
    })
}

/// Format a geographic angle with caller-supplied direction markers.
///
/// # Arguments
///
/// * `value` - Input scalar in centiseconds.
/// * `plus` - Direction marker for a positive angle.
/// * `minus` - Direction marker for a negative angle.
///
/// # Returns
///
/// Complete [`Text`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// This producer has no error channel; its native floating-point/domain
/// behavior is retained.
#[allow(clippy::too_many_arguments)]
pub fn cs2lonlatstr(value: i64, plus: &str, minus: &str) -> Result<Text, Error> {
    Ok(Text {
        value: crate::utils::cs2lonlatstr(value, plus, minus),
    })
}

/// Format centiseconds as a clock string.
///
/// # Arguments
///
/// * `value` - Input scalar in centiseconds.
/// * `separator` - Clock-field separator.
/// * `suppress_zero` - Whether the native formatter suppresses zero seconds.
///
/// # Returns
///
/// Complete [`Text`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// This producer has no error channel; its native floating-point/domain
/// behavior is retained.
#[allow(clippy::too_many_arguments)]
pub fn cs2timestr(value: i64, separator: &str, suppress_zero: bool) -> Result<Text, Error> {
    Ok(Text {
        value: crate::utils::cs2timestr(value, separator, suppress_zero),
    })
}

/// Split an angle into integral and fractional degree components.
///
/// # Arguments
///
/// * `value` - Input scalar in degrees.
/// * `flags` - Calculation bitmask; frame, source and angular units follow the selected operation.
///
/// # Returns
///
/// Complete [`SplitDeg`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// Returns the native input, configuration, data-coverage or calculation
/// failure. Its structured cause is available through `Error::source()`.
#[allow(clippy::too_many_arguments)]
pub fn split_deg(value: f64, flags: i32) -> Result<SplitDeg, Error> {
    let (degrees, minutes, seconds, fraction, sign) =
        native("split_deg", crate::utils::split_deg_checked(value, flags))?;
    Ok(SplitDeg {
        degrees,
        minutes,
        seconds,
        fraction,
        sign,
    })
}

/// Rotate spherical coordinates by the supplied obliquity.
///
/// # Arguments
///
/// * `longitude` - East longitude, right ascension or crossing target in degrees; coordinate flags may select radians.
/// * `latitude` - North latitude or celestial latitude/declination in degrees.
/// * `distance` - Distance in caller-selected length units.
/// * `obliquity` - Obliquity in degrees.
///
/// # Returns
///
/// Complete [`Spherical`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// Returns the native input, configuration, data-coverage or calculation
/// failure. Its structured cause is available through `Error::source()`.
#[allow(clippy::too_many_arguments)]
pub fn cotrans(
    longitude: f64,
    latitude: f64,
    distance: f64,
    obliquity: f64,
) -> Result<Spherical, Error> {
    Ok(Spherical {
        values: native(
            "cotrans",
            crate::engine_api::cotrans(longitude, latitude, distance, obliquity),
        )?,
    })
}

/// Rotate spherical coordinates and their analytic daily rates.
///
/// # Arguments
///
/// * `longitude` - East longitude, right ascension or crossing target in degrees; coordinate flags may select radians.
/// * `latitude` - North latitude or celestial latitude/declination in degrees.
/// * `distance` - Distance in caller-selected length units.
/// * `longitude_speed` - Longitude or right-ascension rate in degrees/day.
/// * `latitude_speed` - Latitude or declination rate in degrees/day.
/// * `distance_speed` - Distance rate in the same length units/day.
/// * `obliquity` - Obliquity in degrees.
///
/// # Returns
///
/// Complete [`SphericalState`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// Returns the native input, configuration, data-coverage or calculation
/// failure. Its structured cause is available through `Error::source()`.
#[allow(clippy::too_many_arguments)]
pub fn cotrans_sp(
    longitude: f64,
    latitude: f64,
    distance: f64,
    longitude_speed: f64,
    latitude_speed: f64,
    distance_speed: f64,
    obliquity: f64,
) -> Result<SphericalState, Error> {
    let (position, velocity) = native(
        "cotrans_sp",
        crate::engine_api::cotrans_sp(
            longitude,
            latitude,
            distance,
            longitude_speed,
            latitude_speed,
            distance_speed,
            obliquity,
        ),
    )?;
    Ok(SphericalState { position, velocity })
}

/// Return the native name for a house-system selector.
///
/// # Arguments
///
/// * `system` - Single-byte house-system selector.
///
/// # Returns
///
/// Complete [`Text`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// This producer has no error channel; its native floating-point/domain
/// behavior is retained.
#[allow(clippy::too_many_arguments)]
pub fn house_name(system: u8) -> Result<Text, Error> {
    Ok(Text {
        value: crate::houses::house_name(char::from(system)).to_owned(),
    })
}

/// Return the name for a sidereal-mode selector.
///
/// # Arguments
///
/// * `mode` - Native sidereal-mode selector.
///
/// # Returns
///
/// Complete [`Text`] values. Optional metadata is populated only when the
/// producer supplies that channel. No array is truncated or recomputed.
///
/// # Errors
///
/// This producer has no error channel; its native floating-point/domain
/// behavior is retained.
#[allow(clippy::too_many_arguments)]
pub fn get_ayanamsa_name(mode: i32) -> Result<Text, Error> {
    Ok(Text {
        value: crate::sidereal::ayanamsha_name(mode).to_owned(),
    })
}
