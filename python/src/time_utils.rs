//! Civil calendars and time-scale functions over the native core clock.

use crate::errors::to_pyerr;
use crate::positions::runtime;
use libephemeris_core::errors::Result;
use libephemeris_core::python_support::Session;
use libephemeris_core::time::{self, Calendar, UtcFields};
use pyo3::prelude::*;

type Civil = (i32, i32, i32, i32, i32, f64);
type CheckedDateResult = (bool, f64, (i32, i32, i32, f64));

fn civil(value: UtcFields) -> Civil {
    (
        value.year,
        value.month,
        value.day,
        value.hour,
        value.minute,
        value.second,
    )
}

fn calendar(py: Python<'_>, flag: i32) -> PyResult<Calendar> {
    Calendar::from_flag(flag).map_err(|error| to_pyerr(py, error))
}

fn session_call<T: Send>(
    py: Python<'_>,
    call: impl FnOnce(&Session) -> Result<T> + Send,
) -> PyResult<T> {
    py.detach(|| call(runtime().shared_session()?.as_ref()))
        .map_err(|error| to_pyerr(py, error))
}

/// Convert a calendar date and decimal hour to a Julian day number.
///
/// year uses astronomical numbering (0 is 1 BCE); month/day are calendar
/// fields and hour is decimal hours, default 12.0 (noon). cal is GREG_CAL
/// (1, proleptic Gregorian) or JUL_CAL (0, Julian), without automatic reform.
/// This is calendar arithmetic, not a UTC/UT1/TT conversion. Out-of-range
/// date fields follow the calendar recurrence; use the root date_conversion
/// to check validity. Invalid cal or nonfinite hour raises ValueError.
///
/// Examples
/// --------
/// ```python
/// import libephemeris as eph
/// jd = eph.julday(2000, 1, 1, 12.0, eph.GREG_CAL)
/// date = eph.revjul(jd, eph.GREG_CAL)
/// ```
#[pyfunction(signature = (year, month, day, hour = 12.0, cal = 1))]
fn julday(py: Python<'_>, year: i32, month: i32, day: i32, hour: f64, cal: i32) -> PyResult<f64> {
    time::julday(year, month, day, hour, calendar(py, cal)?).map_err(|error| to_pyerr(py, error))
}

/// Check a calendar date and return its validity and normalized representation.
///
/// year uses astronomical numbering; month/day are civil fields, hour is
/// decimal hours. cal is 'g' or 'j' (str or ASCII bytes, case-insensitive).
/// Returns (valid, jd, (year, month, day, hour)); invalid calendar dates are
/// normalized and return valid=False, not a date-validation exception.
/// Invalid calendar codes/nonfinite hours raise ValueError; unrepresentable
/// output fields raise OverflowError. No physical time scale is converted.
///
/// Notes
/// -----
/// This is libephemeris.date_conversion. The separate function at
/// libephemeris.time_utils.date_conversion converts calendars and returns
/// only four civil fields; it does not return this checked-date triple.
///
/// Examples
/// --------
/// ```python
/// import libephemeris as eph
/// valid, jd, normalized_date = eph.date_conversion(2024, 2, 30, 0.0, b'g')
/// ```
#[pyfunction(signature = (year, month, day, hour = 12.0, cal = "g".to_owned()), text_signature = "(year, month, day, hour=12.0, cal=b'g')")]
fn date_conversion(
    py: Python<'_>,
    year: i32,
    month: i32,
    day: i32,
    hour: f64,
    #[pyo3(from_py_with = crate::utils::text)] cal: String,
) -> PyResult<CheckedDateResult> {
    py.detach(|| time::convert_date_code(year, month, day, hour, &cal))
        .map(|value| {
            (
                value.valid,
                value.jd,
                (
                    value.date.year,
                    value.date.month,
                    value.date.day,
                    value.date.hour,
                ),
            )
        })
        .map_err(|error| to_pyerr(py, error))
}

/// Convert a Julian day to (year, month, day, decimal_hour).
///
/// jd is a day number on the caller's time scale. cal is GREG_CAL (1,
/// proleptic Gregorian) or JUL_CAL (0); years use astronomical numbering.
/// This changes the calendar representation, not the physical time scale,
/// and does not reconstruct UTC leap-second labels. Invalid cal/nonfinite
/// jd raises ValueError; civil fields outside the supported integer range
/// raise OverflowError. See jdet_to_utc or jdut1_to_utc for UTC fields.
#[pyfunction(signature = (jd, cal = 1))]
fn revjul(py: Python<'_>, jd: f64, cal: i32) -> PyResult<(i32, i32, i32, f64)> {
    time::revjul(jd, calendar(py, cal)?)
        .map(|date| (date.year, date.month, date.day, date.hour))
        .map_err(|error| to_pyerr(py, error))
}

/// Return Delta-T (TT minus UT1) in days at UT1 Julian day tjdut.
///
/// A user override takes precedence; otherwise enabled observed IERS data
/// and the selected model/tidal acceleration determine the value. The unit
/// is days, not the seconds returned by IERS table helpers. Missing model/data
/// support can raise Error subclasses; checked nonfinite input raises
/// ValueError unless an override returns before computation. See deltat_ex
/// for per-call ephemeris-selection flags.
#[pyfunction]
fn deltat(py: Python<'_>, tjdut: f64) -> PyResult<f64> {
    py.detach(|| runtime().delta_t(tjdut))
        .map_err(|error| to_pyerr(py, error))
}

/// Return Delta-T in days at UT1 JD with an ephemeris-selection flag.
///
/// tjdut is UT1 Julian day; flag is the ephemeris bitmask (default
/// FLG_SWIEPH). A user Delta-T override wins first; otherwise an explicit
/// tidal acceleration wins over the flag-selected adjustment for this call.
/// No global setting is changed. Model/data errors propagate as for deltat.
#[pyfunction(signature = (tjdut, flag = 2))]
fn deltat_ex(py: Python<'_>, tjdut: f64, flag: i32) -> PyResult<f64> {
    session_call(py, |session| session.deltat_ex(tjdut, flag))
}

/// Return the weekday of Julian day jd: 0=Monday through 6=Sunday.
///
/// Civil days turn over at midnight, despite Julian days beginning at noon.
/// The input's physical time scale is not converted. Nonfinite jd raises
/// ValueError; no ephemeris or calendar setting is required.
#[pyfunction]
fn day_of_week(py: Python<'_>, jd: f64) -> PyResult<i32> {
    time::day_of_week(jd).map_err(|error| to_pyerr(py, error))
}

/// Convert checked UTC civil fields to (TT Julian day, UT1 Julian day).
///
/// Parameters
/// ----------
/// year, month, day : int
///     Calendar date; astronomical year numbering (0 is 1 BCE).
/// hour, minute : int
///     Hour 0-23 and minute 0-59.
/// second : float
///     Seconds, including fractions. Values in [60, 61) are accepted only
///     at 23:59 on an actual inserted leap-second date.
/// calendar : int, optional
///     GREG_CAL (1) or JUL_CAL (0); default proleptic Gregorian.
///
/// Returns
/// -------
/// tuple of float
///     (jd_tt, jd_ut1), in that order. For Gregorian years 1972 <= year <
///     2035 the conversion uses the bundled leap schedule and observed clock.
///     Before 1972 and from 2035, civil labels denote UT1 and configured
///     Delta-T supplies TT; this is not a forecast of future UTC leap seconds.
///
/// Raises
/// ------
/// Error
///     For invalid civil dates/times or model/data failures.
/// ValueError
///     For unsupported calendars or nonfinite seconds.
///
/// Examples
/// --------
/// ```python
/// import libephemeris as eph
/// jd_tt, jd_ut1 = eph.utc_to_jd(2000, 1, 1, 12, 0, 0.0)
/// utc_fields = eph.jdet_to_utc(jd_tt)
/// ```
#[pyfunction(signature = (year, month, day, hour, minute, second, calendar = 1))]
#[allow(clippy::too_many_arguments, reason = "Public civil-date signature")]
fn utc_to_jd(
    py: Python<'_>,
    year: i32,
    month: i32,
    day: i32,
    hour: i32,
    minute: i32,
    second: f64,
    calendar: i32,
) -> PyResult<(f64, f64)> {
    let calendar = self::calendar(py, calendar)?;
    session_call(py, |session| {
        session.utc_to_jd_numeric(year, month, day, hour, minute, second, calendar)
    })
}

/// Convert TT Julian day jd_et to UTC civil fields.
///
/// calendar is GREG_CAL (1) or JUL_CAL (0). Returns (year, month, day,
/// hour, minute, second), with fractional seconds and astronomical years.
/// The bundled UTC clock applies to reconstructed Gregorian years
/// 1972 <= year < 2035; outside that interval the configured Delta-T model
/// supplies the UT1 civil label. Leap-second labels are retained where
/// applicable. Nonfinite input/invalid calendar raises ValueError; output
/// overflow raises OverflowError and model/data failures propagate as Error.
#[pyfunction(signature = (jd_et, calendar = 1))]
fn jdet_to_utc(py: Python<'_>, jd_et: f64, calendar: i32) -> PyResult<Civil> {
    let calendar = self::calendar(py, calendar)?;
    session_call(py, |session| session.jdet_to_utc_numeric(jd_et, calendar)).map(civil)
}

/// Convert UT1 Julian day jd_ut1 to UTC civil fields.
///
/// calendar is GREG_CAL (1) or JUL_CAL (0). Returns (year, month, day,
/// hour, minute, second), with astronomical years and fractional seconds.
/// Modern conversion (Gregorian years 1972 <= year < 2035) uses the
/// independent bundled observed clock, not a user model Delta-T override.
/// Outside that interval the UT1 input is labelled directly. Nonfinite
/// input/invalid calendar raises ValueError; field overflow raises
/// OverflowError and clock/data failures propagate as Error.
#[pyfunction(signature = (jd_ut1, calendar = 1))]
fn jdut1_to_utc(py: Python<'_>, jd_ut1: f64, calendar: i32) -> PyResult<Civil> {
    let calendar = self::calendar(py, calendar)?;
    session_call(py, |session| session.jdut1_to_utc_numeric(jd_ut1, calendar)).map(civil)
}

/// Convert local Gregorian civil fields to UTC using a fixed hour offset.
///
/// year/month/day use astronomical numbering; hour/minute are integers and
/// second may be fractional. timezone_offset is hours east of UTC, so +1
/// subtracts one hour from the supplied local time. Returns (year, month,
/// day, hour, minute, second), carrying across date boundaries. Negative
/// offsets can perform the inverse conversion. No named-zone, DST or
/// historical timezone rules are consulted. Nonfinite operands raise
/// ValueError; unrepresentable output fields raise OverflowError.
#[pyfunction]
#[allow(clippy::too_many_arguments, reason = "Public civil-date signature")]
fn utc_time_zone(
    py: Python<'_>,
    year: i32,
    month: i32,
    day: i32,
    hour: i32,
    minute: i32,
    second: f64,
    timezone_offset: f64,
) -> PyResult<Civil> {
    time::utc_time_zone(year, month, day, hour, minute, second, timezone_offset)
        .map(civil)
        .map_err(|error| to_pyerr(py, error))
}

/// Return the equation of time: apparent minus mean solar time, in days.
///
/// jd is UT1 Julian day. A positive result means a sundial is ahead of
/// mean solar time. The Sun's apparent right ascension and sidereal clock
/// use the configured native calculation path; data, coverage and clock
/// errors propagate as Error subclasses. Checked nonfinite inputs raise
/// ValueError. See lat_to_lmt and lmt_to_lat for local solar-time conversion.
#[pyfunction]
fn time_equ(py: Python<'_>, jd: f64) -> PyResult<f64> {
    session_call(py, |session| session.time_equ(jd))
}

/// Convert local apparent solar time to local mean solar time.
///
/// jd_lat is a Julian day expressed in local apparent time, not UTC;
/// longitude is geographic degrees east. Returns the local mean-time
/// Julian day, retaining its local longitude offset. The equation-of-time
/// calculation can raise the same data/clock errors as time_equ; checked
/// nonfinite inputs raise ValueError. See lmt_to_lat for the inverse.
#[pyfunction]
fn lat_to_lmt(py: Python<'_>, jd_lat: f64, longitude: f64) -> PyResult<f64> {
    session_call(py, |session| session.lat_to_lmt(jd_lat, longitude))
}

/// Convert local mean solar time to local apparent solar time.
///
/// jd_lmt is a local mean-time Julian day; longitude is degrees east.
/// Returns the local apparent-time Julian day, not UT1 or UTC. Native
/// equation-of-time data/clock errors propagate as for time_equ; checked
/// nonfinite inputs raise ValueError. See lat_to_lmt for the inverse.
#[pyfunction]
fn lmt_to_lat(py: Python<'_>, jd_lmt: f64, longitude: f64) -> PyResult<f64> {
    session_call(py, |session| session.lmt_to_lat(jd_lmt, longitude))
}

/// Return local apparent sidereal time in hours at UT1 Julian day jd.
///
/// longitude is geographic degrees east (default Greenwich, 0). obliquity
/// is true obliquity and nutation is nutation in longitude, both degrees.
/// Each None component is acquired independently from the native ECL_NUT
/// result; a supplied component is retained. The result is reduced modulo
/// 24 hours. Frame/clock/data failures propagate as Error subclasses;
/// checked nonfinite numerical inputs raise ValueError. See sidtime0 for
/// explicitly supplied Greenwich obliquity/nutation.
#[pyfunction(signature = (jd, longitude = 0.0, obliquity = None, nutation = None))]
fn sidtime(
    py: Python<'_>,
    jd: f64,
    longitude: f64,
    obliquity: Option<f64>,
    nutation: Option<f64>,
) -> PyResult<f64> {
    session_call(py, |session| {
        session.sidtime_with(jd, longitude, obliquity, nutation)
    })
}

/// Return Greenwich apparent sidereal time in hours at UT1 Julian day jd.
///
/// obliquity is true obliquity and nutation is nutation in longitude, both
/// in degrees, not radians. Returns sidereal hours modulo 24, not degrees.
/// The configured model clock still applies. Clock/data errors propagate;
/// checked nonfinite numerical inputs raise ValueError. See sidtime for
/// optional native frame components and an east-positive site longitude.
#[pyfunction]
fn sidtime0(py: Python<'_>, jd: f64, obliquity: f64, nutation: f64) -> PyResult<f64> {
    session_call(py, |session| session.sidtime0(jd, obliquity, nutation))
}

/// Explicit Greenwich IAU2006 sidereal hours with independent UT1/TT Julian dates.
///
/// apparent=True uses the full IAU2006/2000A true frame and CIO locator;
/// False returns GMST. TT is limited to JD2396757.5..2469808.5 and UT1 to
/// 200000 Julian years either side of J2000. No session, Delta-T, data or
/// network is acquired. Existing sidtime/sidtime0 input domains are unchanged.
#[pyfunction(signature = (jd_ut1, jd_tt, *, apparent = true))]
fn sidtime_iau2006(py: Python<'_>, jd_ut1: f64, jd_tt: f64, apparent: bool) -> PyResult<f64> {
    py.detach(|| time::sidereal_time_iau2006_hours(jd_ut1, jd_tt, apparent))
        .map_err(|error| to_pyerr(py, error))
}

/// Return the leap-table TAI-minus-UTC offset in seconds at Julian day jd.
///
/// jd is the UTC calendar day coordinate used for table selection, not MJD.
/// The configured session supplies the leap table; this is a stepwise offset,
/// not Delta-T (TT minus UT1). Nonfinite jd raises ValueError; table loading
/// failures propagate as Error subclasses.
#[pyfunction]
fn get_tai_utc_for_jd(py: Python<'_>, jd: f64) -> PyResult<f64> {
    session_call(py, |session| {
        time::tai_utc_for_jd(jd, &session.leap_seconds()?)
    })
}

/// Convert UTC calendar fields to a scalar TAI Julian day using the leap table.
///
/// year/month/day use astronomical numbering; hour/minute are integer fields
/// and second may be fractional. calendar is GREG_CAL (1) or JUL_CAL (0).
/// Unlike utc_to_jd, this is direct calendar/leap-table arithmetic: it does
/// not perform that function's civil-field/leap-label validation or its
/// 1972-2035 observation-clock policy. A scalar UTC day cannot uniquely label
/// an inserted leap second; use utc_to_jd for checked civil UTC conversion.
/// Invalid calendar/nonfinite arithmetic raises ValueError; leap-table
/// acquisition failures propagate as Error subclasses.
#[pyfunction(signature = (year, month, day, hour, minute, second, calendar = 1))]
#[allow(clippy::too_many_arguments, reason = "Public civil-date signature")]
fn utc_to_tai_jd(
    py: Python<'_>,
    year: i32,
    month: i32,
    day: i32,
    hour: i32,
    minute: i32,
    second: f64,
    calendar: i32,
) -> PyResult<f64> {
    let calendar = self::calendar(py, calendar)?;
    session_call(py, |session| {
        time::utc_to_tai_jd(
            year,
            month,
            day,
            hour,
            minute,
            second,
            calendar,
            &session.leap_seconds()?,
        )
    })
}

/// Convert scalar TAI Julian day jd_tai to UTC calendar fields.
///
/// calendar is GREG_CAL (1) or JUL_CAL (0). Returns (year, month, day,
/// hour, minute, second) using the session leap table and stepwise offset
/// refinement, not utc_to_jd's modern observed-clock conversion. Scalar
/// calendar reconstruction does not preserve a distinct inserted leap-second
/// label. Invalid calendar/nonfinite input raises ValueError, field overflow
/// raises OverflowError and leap-table acquisition failures propagate as Error.
#[pyfunction(signature = (jd_tai, calendar = 1))]
fn tai_jd_to_utc(py: Python<'_>, jd_tai: f64, calendar: i32) -> PyResult<Civil> {
    let calendar = self::calendar(py, calendar)?;
    session_call(py, |session| {
        time::tai_jd_to_utc(jd_tai, calendar, &session.leap_seconds()?)
    })
    .map(civil)
}

/// Convert TT Julian day jd_tt to TAI JD by subtracting 32.184 seconds.
///
/// Returns a float in days on the TAI scale. No leap table or ephemeris is
/// needed; nonfinite input raises ValueError. See tai_to_tt_jd for the inverse.
#[pyfunction]
fn tt_to_tai_jd(py: Python<'_>, jd_tt: f64) -> PyResult<f64> {
    time::tt_to_tai_jd(jd_tt).map_err(|error| to_pyerr(py, error))
}

/// Convert TAI Julian day jd_tai to TT JD by adding 32.184 seconds.
///
/// Returns a float in days on the TT scale. No leap table or ephemeris is
/// needed; nonfinite input raises ValueError. See tt_to_tai_jd for the inverse.
#[pyfunction]
fn tai_to_tt_jd(py: Python<'_>, jd_tai: f64) -> PyResult<f64> {
    time::tai_to_tt_jd(jd_tai).map_err(|error| to_pyerr(py, error))
}

/// Convert TAI Julian day jd_tai to a scalar UTC Julian day.
///
/// Uses the session leap table and refined TAI-minus-UTC offset. The float
/// result is a calendar day coordinate, not a leap-second-preserving civil
/// label. Nonfinite input raises ValueError; leap-table acquisition errors
/// propagate. See jdet_to_utc for TT-to-UTC civil fields.
#[pyfunction]
fn tai_to_utc_jd(py: Python<'_>, jd_tai: f64) -> PyResult<f64> {
    session_call(py, |session| {
        time::tai_to_utc_jd(jd_tai, &session.leap_seconds()?)
    })
}

/// Convert scalar UTC Julian day jd_utc to TAI Julian day.
///
/// Adds the session leap-table TAI-minus-UTC offset in days. A scalar UTC
/// calendar day cannot distinguish an inserted leap-second label; use
/// utc_to_jd for checked civil UTC input. Nonfinite input raises ValueError;
/// leap-table acquisition errors propagate. See tai_to_utc_jd for the inverse.
#[pyfunction]
fn utc_jd_to_tai(py: Python<'_>, jd_utc: f64) -> PyResult<f64> {
    session_call(py, |session| {
        time::utc_jd_to_tai(jd_utc, &session.leap_seconds()?)
    })
}

/// Convert a historical civil date to the requested Julian/Gregorian calendar.
///
/// Exposed as libephemeris.time_utils.date_conversion. year/month/day/hour
/// describe the input civil label, with astronomical years and decimal hours;
/// its input calendar changes at 1582-10-15. calendar selects the output:
/// 'j' or 'g', str or ASCII bytes, case-insensitive. Returns only (year,
/// month, day, hour). No UTC/TT conversion is performed. Invalid calendar
/// codes or nonfinite hours raise ValueError; field overflow raises
/// OverflowError. Unlike libephemeris.date_conversion, this function does
/// not return (valid, jd, normalized_date).
#[pyfunction(name = "_time_date_conversion")]
fn module_date_conversion(
    py: Python<'_>,
    year: i32,
    month: i32,
    day: i32,
    hour: f64,
    #[pyo3(from_py_with=crate::utils::text)] calendar: String,
) -> PyResult<(i32, i32, i32, f64)> {
    py.detach(|| {
        libephemeris_core::python_support::module_positions::convert_calendar_date(
            year, month, day, hour, &calendar,
        )
    })
    .map(|date| (date.year, date.month, date.day, date.hour))
    .map_err(|error| to_pyerr(py, error))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(module_date_conversion, module)?)?;
    macro_rules! add {
        ($($name:ident),+ $(,)?) => { $(module.add_function(wrap_pyfunction!($name, module)?)?;)+ };
    }
    add!(
        julday,
        date_conversion,
        revjul,
        deltat,
        deltat_ex,
        day_of_week,
        utc_to_jd,
        jdet_to_utc,
        jdut1_to_utc,
        utc_time_zone,
        time_equ,
        lat_to_lmt,
        lmt_to_lat,
        sidtime,
        sidtime0,
        sidtime_iau2006,
        get_tai_utc_for_jd,
        utc_to_tai_jd,
        tai_jd_to_utc,
        tt_to_tai_jd,
        tai_to_tt_jd,
        tai_to_utc_jd,
        utc_jd_to_tai
    );
    Ok(())
}
