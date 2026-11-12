//! Native composition for the optional calculation convenience module.

use super::Session;
use crate::arabic_parts::{self, ArabicError, ArabicResult};
use crate::errors::{Error, Result};
use crate::time::{self, Calendar};

/// Public nakshatra identifiers in catalog order,0 through26.
pub const NAKSHATRA_IDS: &[(&str, i32)] = &[
    ("ASWINI", 0),
    ("BHARANI", 1),
    ("KRITIKHA", 2),
    ("ROHINI", 3),
    ("MRIGASIRA", 4),
    ("ARIDRA", 5),
    ("PUNARVASU", 6),
    ("PUSHYAMI", 7),
    ("ASLESHA", 8),
    ("MAKHA", 9),
    ("PUBBA", 10),
    ("UTTARA", 11),
    ("HASTA", 12),
    ("CHITTA", 13),
    ("SWATHI", 14),
    ("VISHAKA", 15),
    ("ANURADHA", 16),
    ("JYESTA", 17),
    ("MOOLA", 18),
    ("POORVASHADA", 19),
    ("UTTARASHADA", 20),
    ("SRAVANA", 21),
    ("DHANISHTA", 22),
    ("SATABISHA", 23),
    ("POORVABHADRA", 24),
    ("UTTARABHADRA", 25),
    ("REVATHI", 26),
];

/// Current UTC civil Julian day, measured by the operating-system clock.
/// The Unix epoch is JD2440587.5; no TT/UT1 correction belongs to this clock.
pub fn current_utc_jd() -> f64 {
    let seconds = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or_else(
            |error| -error.duration().as_secs_f64(),
            |duration| duration.as_secs_f64(),
        );
    2440587.5 + seconds / 86400.0
}

/// ISO Gregorian timestamp rounded to the nearest second, ties to even.
/// Civil years1..9999 match the public timestamp domain. Nonfinite epochs or
/// out-of-domain dates return a typed date error; no timezone is inferred.
pub fn iso_timestamp(jd: f64) -> ArabicResult<String> {
    let failure = || ArabicError::Parse {
        what: "Gregorian Julian day in years1..9999",
        input: jd.to_string(),
    };
    let date = time::revjul(jd, Calendar::Gregorian).map_err(|_| failure())?;
    if !(1..=9999).contains(&date.year) {
        return Err(failure());
    }
    let seconds = (date.hour * 3600.0).round_ties_even() as i64;
    let (year, month, day, seconds) = if seconds >= 86400 {
        let start = time::julday(date.year, date.month, date.day, 0.0, Calendar::Gregorian)
            .map_err(|_| failure())?;
        let next = time::revjul(start + 1.0, Calendar::Gregorian).map_err(|_| failure())?;
        (next.year, next.month, next.day, seconds - 86400)
    } else {
        (date.year, date.month, date.day, seconds)
    };
    if !(1..=9999).contains(&year) {
        return Err(failure());
    }
    Ok(format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}",
        seconds / 3600,
        seconds % 3600 / 60,
        seconds % 60
    ))
}

/// Human-readable house-system name for a complete one-character code.
/// Empty or multicharacter keys are unknown, rather than truncated.
pub fn house_name(code: &str) -> &'static str {
    let mut characters = code.chars();
    match (characters.next(), characters.next()) {
        (Some(code), None) => arabic_parts::house_system_name(code),
        _ => "Unknown",
    }
}

/// Fixed longitude or moving body for a bounded aspect search.
#[derive(Clone, Copy, Debug)]
pub enum AspectTarget {
    /// Ecliptic longitude in degrees.
    Longitude(f64),
    /// Public body identifier.
    Body(i32),
}

impl Session {
    /// Twelve Sripati house cusps in degrees at UT1 Julian day.
    /// Geographic longitude/latitude are east/north degrees. The native
    /// house model owns polar-domain failures and the selected time model.
    pub fn contrib_raman_houses(
        &self,
        jd: f64,
        longitude: f64,
        latitude: f64,
    ) -> Result<[f64; 12]> {
        self.houses(jd, latitude, longitude, b'S')
            .map(|value| value.cusps)
    }

    /// Apparent ecliptic longitudes of Aldebaran, Regulus, Antares, Fomalhaut.
    /// Each place uses the public fixed-star route at UT1 JD. An unavailable
    /// place has a NaN sentinel and a native debug log, preserving result shape.
    pub fn contrib_saturn_four_stars(&self, jd: f64) -> [f64; 4] {
        ["Aldebaran", "Regulus", "Antares", "Fomalhaut"].map(|name| {
            self.fixstar_ut(name, jd, 2)
                .map(|value| value.values[0])
                .unwrap_or_else(|error| {
                    log::debug!("fixed-star longitude {name} at JD{jd} unavailable: {error}");
                    f64::NAN
                })
        })
    }

    /// First sampled aspect within the specified orb over at most three years.
    ///
    /// UT1 days and ecliptic angular degrees use ordinary binary64 positions.
    /// The prescribed0.5-day sampling grid is a coarse convenience search,
    /// not an exact crossing solver. It returns(sampleJD,angular residual).
    /// A completed miss returns(NaN,+infinity), never a transport timeout.
    pub fn contrib_next_aspect(
        &self,
        body: i32,
        target: AspectTarget,
        aspect: f64,
        start: f64,
        orb: f64,
    ) -> Result<(f64, f64)> {
        for value in [start, aspect, orb] {
            if !value.is_finite() {
                return Err(Error::NonFiniteArgument { value });
            }
        }
        if let AspectTarget::Longitude(value) = target {
            if !value.is_finite() {
                return Err(Error::NonFiniteArgument { value });
            }
        }
        let mut jd = start;
        while jd - start < 365.25 * 3.0 {
            let first = self.calc_ut(jd, body, crate::constants::FLG_SPEED)?.values[0];
            let second = match target {
                AspectTarget::Longitude(value) => value,
                AspectTarget::Body(body) => {
                    self.calc_ut(jd, body, crate::constants::FLG_SPEED)?.values[0]
                }
            };
            let delta = arabic_parts::angular_separation(first, second)
                .map_err(|_| Error::NonFiniteArgument { value: first })?;
            let residual = (delta - aspect).abs();
            if residual <= orb {
                return Ok((jd, residual));
            }
            jd += 0.5;
        }
        Ok((f64::NAN, f64::INFINITY))
    }

    /// First0.5-day sample with negative longitude speed after a nonnegative
    /// preceding sample. Times are UT1 JD and the horizon is in days.
    /// A completed search without a change returns NaN; data errors propagate.
    pub fn contrib_next_retro(&self, body: i32, start: f64, max_days: f64) -> Result<f64> {
        for value in [start, max_days] {
            if !value.is_finite() {
                return Err(Error::NonFiniteArgument { value });
            }
        }
        let mut jd = start;
        let mut last = None;
        while jd - start < max_days {
            let speed = self.calc_ut(jd, body, crate::constants::FLG_SPEED)?.values[3];
            if last.is_some_and(|last| last >= 0.0) && speed < 0.0 {
                return Ok(jd);
            }
            last = Some(speed);
            jd += 0.5;
        }
        Ok(f64::NAN)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn timestamp_rounds_seconds_and_rolls_civil_midnight() {
        assert_eq!(iso_timestamp(2451545.0).unwrap(), "2000-01-01T12:00:00");
        assert_eq!(
            iso_timestamp(2451545.5 - 0.1 / 86400.0).unwrap(),
            "2000-01-02T00:00:00"
        );
        assert_eq!(house_name("SS"), "Unknown");
        assert_eq!(house_name("S"), "Sripati");
    }
}
