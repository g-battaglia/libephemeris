//! Calendar calculations over retained model configuration, without images.
use super::*;

impl Session {
    /// Convert a finite civil UTC label to `(TT JD, UT1 JD)` using ordinary
    /// bounded calendar fields. Years use astronomical numbering, angles are
    /// absent, and seconds may include an actual inserted leap second.
    ///
    /// Before 1972 and from 2035 the civil label denotes UT1 and this session's
    /// Delta-T supplies TT. Modern UTC uses the bundled leap schedule, the
    /// exact 32.184-second TT/TAI offset and the bundled observed Delta-T curve.
    /// Split day/fraction arithmetic retains subsecond precision. Invalid
    /// civil fields, nonfinite seconds and unsupported calendars are errors.
    #[allow(clippy::too_many_arguments)]
    pub fn utc_to_jd_numeric(
        &self,
        year: i32,
        month: i32,
        day: i32,
        hour: i32,
        minute: i32,
        second: f64,
        calendar: crate::time::Calendar,
    ) -> Result<(f64, f64), Error> {
        use crate::time::{Calendar, days_in_month, is_leap_second_date, julday, revjul};
        if !(1..=12).contains(&month) || day < 1 || day > days_in_month(year, month, calendar) {
            return Err(Error::InvalidDate { year, month, day });
        }
        if !second.is_finite() {
            return Err(Error::NonFiniteArgument { value: second });
        }
        if !(0..=23).contains(&hour) || !(0..=59).contains(&minute) || second < 0.0 {
            return Err(Error::InvalidTime {
                hour,
                minute,
                second,
            });
        }
        let gregorian = if calendar == Calendar::Julian {
            let date = revjul(
                julday(year, month, day, 12.0, calendar)?,
                Calendar::Gregorian,
            )?;
            (date.year, date.month, date.day)
        } else {
            (year, month, day)
        };
        if second >= 60.0
            && (second >= 61.0
                || hour != 23
                || minute != 59
                || gregorian.0 < 1972
                || !is_leap_second_date(gregorian.0, gregorian.1, gregorian.2))
        {
            return Err(Error::InvalidTime {
                hour,
                minute,
                second,
            });
        }
        if gregorian.0 < 1972 || gregorian.0 >= 2035 {
            let decimal_hour = hour as f64 + minute as f64 / 60.0 + second / 3600.0;
            let ut = julday(year, month, day, decimal_hour, calendar)?;
            return Ok((ut + self.deltat(ut)?, ut));
        }
        let midnight = julday(
            gregorian.0,
            gregorian.1,
            gregorian.2,
            0.0,
            Calendar::Gregorian,
        )?;
        let (mut seconds, fraction) = numeric_split(midnight * NUMERIC_DAY, 1.0);
        seconds += numeric_leap_offset(seconds, false);
        let clock = hour as f64 * 3600.0 + minute as f64 * 60.0 + second;
        let (clock_seconds, subsecond) = numeric_split(fraction + clock, 1.0);
        seconds += clock_seconds;
        let (whole, fraction) = numeric_split(seconds, NUMERIC_DAY);
        let tt_fraction = (fraction + subsecond) / NUMERIC_DAY + NUMERIC_TT_TAI;
        let tt = whole + tt_fraction;
        let delta = crate::smh2016::bundled_deltat_seconds(tt)?;
        Ok((tt, whole + (tt_fraction - delta / NUMERIC_DAY)))
    }

    /// Convert finite TT JD to UTC fields, with the bundled UTC clock during
    /// 1972–2035 and configured Delta-T outside that interval. Output seconds
    /// retain leap labels and a representation-derived whole-second snap.
    /// Calendar years must fit `i32`; model, calendar and finite-domain errors
    /// propagate without reconstructing language-specific exception state.
    pub fn jdet_to_utc_numeric(
        &self,
        jd_et: f64,
        calendar: crate::time::Calendar,
    ) -> Result<crate::time::UtcFields, Error> {
        use crate::time::{Calendar, UTC_EPOCH_BAND_DAYS, UTC_EPOCH_JD, revjul};
        if !jd_et.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_et });
        }
        let first = jd_et - self.deltat(jd_et)?;
        let estimate = jd_et - self.deltat(first)?;
        let year = revjul(estimate, Calendar::Gregorian)?.year;
        let near_epoch = (estimate - UTC_EPOCH_JD).abs() < UTC_EPOCH_BAND_DAYS;
        if (year < 1972 && !near_epoch) || year >= 2035 {
            return numeric_calendar_label(estimate, calendar);
        }
        let (whole, fraction) = if near_epoch {
            crate::smh2016::timescale_parts_from_ut1(estimate)?
        } else {
            numeric_split(jd_et, 1.0)
        };
        let fields = numeric_modern_fields(whole, fraction)?;
        if fields.year < 1972 {
            return numeric_calendar_label(estimate, calendar);
        }
        numeric_finish_fields(fields, estimate, calendar)
    }

    /// Convert finite UT1 JD to UTC calendar fields using the independent
    /// bundled observed clock. The session's model Delta-T does not change
    /// modern UTC. Before 1972 and from 2035 the input is labelled directly.
    /// Returns ordinary bounded fields; nonfinite or out-of-range dates fail.
    pub fn jdut1_to_utc_numeric(
        &self,
        jd_ut1: f64,
        calendar: crate::time::Calendar,
    ) -> Result<crate::time::UtcFields, Error> {
        use crate::time::{Calendar, UTC_EPOCH_BAND_DAYS, UTC_EPOCH_JD, revjul};
        let year = revjul(jd_ut1, Calendar::Gregorian)?.year;
        let near_epoch = (jd_ut1 - UTC_EPOCH_JD).abs() < UTC_EPOCH_BAND_DAYS;
        if (year < 1972 && !near_epoch) || year >= 2035 {
            return numeric_calendar_label(jd_ut1, calendar);
        }
        let (whole, fraction) = crate::smh2016::timescale_parts_from_ut1(jd_ut1)?;
        let fields = numeric_modern_fields(whole, fraction)?;
        if fields.year < 1972 {
            return numeric_calendar_label(jd_ut1, calendar);
        }
        numeric_finish_fields(fields, jd_ut1, calendar)
    }

    /// Delta-T days at UT1 JD using the explicit or flag-selected tidal term.
    /// A caller user override returns before flag/model/IERS consumption. An
    /// explicit tidal acceleration wins over ephemeris bits; otherwise the
    /// shared flag policy selects the adjustment for this call alone. No
    /// setting is mutated. The source-facing model chain preserves its actual
    /// primitive/error domains rather than adding a finite-epoch preflight.
    pub fn deltat_ex(&self, jd_ut: f64, flag: i32) -> Result<f64, Error> {
        let mut cfg = self
            .ctx
            .delta_t_config(self.delta_t_model(), self.ctx.effective_iers_delta_t());
        if let Some(value) = cfg.user_override_days {
            return Ok(value);
        }
        cfg.tidal_acceleration =
            crate::time::resolve_tidal_acceleration(flag, self.ctx.tidal_acceleration_override());
        crate::time::deltat_with_config(jd_ut, &cfg)
    }
}

/// SI seconds in a Julian day and the exact TT minus TAI day fraction.
const NUMERIC_DAY: f64 = 86400.0;
const NUMERIC_TT_TAI: f64 = 32.184 / NUMERIC_DAY;

/// Split a quantity into an integral quotient and a nonnegative remainder.
fn numeric_split(value: f64, divisor: f64) -> (f64, f64) {
    let quotient = (value / divisor).floor();
    (quotient, value - quotient * divisor)
}

/// Offset of the bundled UTC/TAI clock, with the inserted second represented
/// by its continuous one-second ramp. Absolute input seconds start at JD zero.
fn numeric_leap_offset(seconds: f64, tai: bool) -> f64 {
    let mut previous = None;
    for row in crate::time::hardcoded_leap_seconds() {
        let utc = (row.mjd() + 2400000.5) * NUMERIC_DAY;
        for minus in [1.0, 0.0] {
            let y = row.tai_utc() - minus;
            let x = utc - minus + if tai { y } else { 0.0 };
            if let Some((x0, y0)) = previous {
                if seconds < x {
                    return (y - y0) / (x - x0) * (seconds - x0) + y0;
                }
            } else if seconds <= x {
                return y;
            }
            if seconds == x {
                return y;
            }
            previous = Some((x, y));
        }
    }
    previous.expect("the bundled leap schedule is nonempty").1
}

/// UTC civil fields of a split TT day. Integral second arithmetic remains
/// exactly representable across the supported ephemeris date domain.
fn numeric_modern_fields(whole: f64, fraction: f64) -> Result<crate::time::UtcFields, Error> {
    use crate::time::{Calendar, UtcFields, revjul};
    let (seconds, remainder) = numeric_split(whole * NUMERIC_DAY, 1.0);
    let (more_seconds, subsecond) =
        numeric_split(remainder + (fraction - NUMERIC_TT_TAI) * NUMERIC_DAY, 1.0);
    let seconds = seconds + more_seconds;
    let (offset, leap_fraction) = numeric_split(numeric_leap_offset(seconds, true), 1.0);
    let leap = if leap_fraction > 0.0 { 1.0 } else { 0.0 };
    let raw_second = seconds - offset - leap;
    let (day, clock) = numeric_split(raw_second + 43200.0, NUMERIC_DAY);
    let date = revjul(day, Calendar::Gregorian)?;
    let hour = (clock / 3600.0).floor();
    let minute = ((clock - hour * 3600.0) / 60.0).floor();
    Ok(UtcFields {
        year: date.year,
        month: date.month,
        day: date.day,
        hour: hour as i32,
        minute: minute as i32,
        second: clock - hour * 3600.0 - minute * 60.0 + leap + subsecond,
    })
}

/// Split a scalar day into civil fields and remove whole-second roundoff.
fn numeric_calendar_label(
    jd: f64,
    calendar: crate::time::Calendar,
) -> Result<crate::time::UtcFields, Error> {
    let date = crate::time::revjul(jd, calendar)?;
    let hour = date.hour.trunc() as i32;
    let minutes = (date.hour - hour as f64) * 60.0;
    let minute = minutes.trunc() as i32;
    numeric_snap_fields(
        crate::time::UtcFields {
            year: date.year,
            month: date.month,
            day: date.day,
            hour,
            minute,
            second: (minutes - minute as f64) * 60.0,
        },
        jd,
        calendar,
        false,
    )
}

/// Snap within two input-day ulps, carrying an ordinary midnight into its
/// calendar date and retaining an actual 23:59:60 insertion when requested.
fn numeric_snap_fields(
    mut fields: crate::time::UtcFields,
    jd: f64,
    calendar: crate::time::Calendar,
    allow_leap: bool,
) -> Result<crate::time::UtcFields, Error> {
    let rounded = fields.second.round_ties_even();
    let magnitude = jd.abs();
    let tolerance = 2.0 * (magnitude.next_up() - magnitude) * NUMERIC_DAY;
    if fields.second != rounded && (fields.second - rounded).abs() < tolerance {
        let valid_leap = allow_leap
            && rounded == 60.0
            && fields.hour == 23
            && fields.minute == 59
            && crate::time::is_leap_second_date(fields.year, fields.month, fields.day);
        if rounded < 60.0 || valid_leap {
            fields.second = rounded;
        } else {
            let total = fields.hour * 3600 + fields.minute * 60 + rounded as i32;
            if total < 86400 {
                fields.hour = total / 3600;
                fields.minute = total % 3600 / 60;
                fields.second = (total % 60) as f64;
            } else {
                let next = crate::time::revjul(
                    crate::time::julday(fields.year, fields.month, fields.day, 12.0, calendar)?
                        + 1.0,
                    calendar,
                )?;
                fields.year = next.year;
                fields.month = next.month;
                fields.day = next.day;
                fields.hour = 0;
                fields.minute = 0;
                fields.second = 0.0;
            }
        }
    }
    Ok(fields)
}

/// Finish modern Gregorian UTC, then relabel the date in a requested Julian
/// calendar without changing its clock or physical instant.
fn numeric_finish_fields(
    fields: crate::time::UtcFields,
    estimate: f64,
    calendar: crate::time::Calendar,
) -> Result<crate::time::UtcFields, Error> {
    use crate::time::{Calendar, julday, revjul};
    let mut fields = numeric_snap_fields(fields, estimate, Calendar::Gregorian, true)?;
    if calendar == Calendar::Julian {
        let date = revjul(
            julday(
                fields.year,
                fields.month,
                fields.day,
                12.0,
                Calendar::Gregorian,
            )?,
            calendar,
        )?;
        fields.year = date.year;
        fields.month = date.month;
        fields.day = date.day;
    }
    Ok(fields)
}

#[cfg(test)]
mod numeric_calendar_tests {
    use super::*;
    use crate::context::EphemerisContext;
    use crate::time::{Calendar, julday};

    fn session(delta: f64) -> Session {
        Session::from_context(EphemerisContext::new().with_user_delta_t(Some(delta)))
    }

    #[test]
    fn modern_utc_uses_atomic_offset_and_an_independent_earth_rotation_clock() {
        let first = session(64.0 / NUMERIC_DAY);
        let second = session(1000.0 / NUMERIC_DAY);
        let output = first
            .utc_to_jd_numeric(2000, 1, 1, 12, 0, 0.0, Calendar::Gregorian)
            .unwrap();
        assert_eq!(output.0, 2451545.0 + (32.0 + 32.184) / NUMERIC_DAY);
        assert_eq!(
            output,
            second
                .utc_to_jd_numeric(2000, 1, 1, 12, 0, 0.0, Calendar::Gregorian)
                .unwrap()
        );
        let expected_delta = crate::smh2016::bundled_deltat_seconds(output.0).unwrap();
        assert!(((output.0 - output.1) * NUMERIC_DAY - expected_delta).abs() < 5e-5);
        assert_eq!(
            first
                .jdut1_to_utc_numeric(output.1, Calendar::Gregorian)
                .unwrap(),
            second
                .jdut1_to_utc_numeric(output.1, Calendar::Gregorian)
                .unwrap()
        );
    }

    #[test]
    fn inserted_second_is_preserved_through_tt_and_ut1_calendar_routes() {
        let session = session(64.0 / NUMERIC_DAY);
        for second in [59.0, 60.0, 60.25, 60.75] {
            let (tt, ut) = session
                .utc_to_jd_numeric(2016, 12, 31, 23, 59, second, Calendar::Gregorian)
                .unwrap();
            for fields in [
                session
                    .jdet_to_utc_numeric(tt, Calendar::Gregorian)
                    .unwrap(),
                session
                    .jdut1_to_utc_numeric(ut, Calendar::Gregorian)
                    .unwrap(),
            ] {
                assert_eq!(
                    (
                        fields.year,
                        fields.month,
                        fields.day,
                        fields.hour,
                        fields.minute
                    ),
                    (2016, 12, 31, 23, 59)
                );
                assert!(
                    (fields.second - second).abs() < 5e-5,
                    "second={second}, fields={fields:?}"
                );
            }
        }
        assert!(matches!(
            session.utc_to_jd_numeric(2017, 12, 31, 23, 59, 60.0, Calendar::Gregorian),
            Err(Error::InvalidTime { .. })
        ));
    }

    #[test]
    fn outer_era_uses_configured_delta_and_julian_labels_refer_to_same_instant() {
        let delta = 0.002;
        let session = session(delta);
        let (tt, ut) = session
            .utc_to_jd_numeric(1900, 1, 1, 0, 0, 0.0, Calendar::Gregorian)
            .unwrap();
        assert_eq!(ut, julday(1900, 1, 1, 0.0, Calendar::Gregorian).unwrap());
        assert_eq!(tt, ut + delta);
        assert_eq!(
            session
                .utc_to_jd_numeric(1999, 12, 19, 12, 0, 0.0, Calendar::Julian)
                .unwrap(),
            session
                .utc_to_jd_numeric(2000, 1, 1, 12, 0, 0.0, Calendar::Gregorian)
                .unwrap()
        );
    }

    #[test]
    fn numeric_calendar_rejects_nonfinite_and_invalid_civil_inputs() {
        let session = session(64.0 / NUMERIC_DAY);
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                session.utc_to_jd_numeric(2000, 1, 1, 0, 0, value, Calendar::Gregorian),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                session.jdet_to_utc_numeric(value, Calendar::Gregorian),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                session.jdut1_to_utc_numeric(value, Calendar::Gregorian),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
        assert!(matches!(
            session.utc_to_jd_numeric(2024, 2, 30, 12, 0, 0.0, Calendar::Gregorian),
            Err(Error::InvalidDate { .. })
        ));
        assert!(matches!(
            session.utc_to_jd_numeric(2024, 1, 1, 24, 0, 0.0, Calendar::Gregorian),
            Err(Error::InvalidTime { .. })
        ));
    }
}
