//! Native sequence inputs for heliacal visibility and event APIs.

use crate::engine_api::Session;
use crate::errors::Result;

fn site(words: &[f64]) -> Result<[f64; 3]> {
    let site = crate::heliacal::parse_pheno_geopos(words)?;
    Ok([site.lon_deg, site.lat_deg, site.altitude_m])
}

fn atmosphere(words: &[f64], search: bool) -> [f64; 4] {
    let mut atmosphere = [1013.25, 15.0, if search { 40.0 } else { 50.0 }, 0.0];
    for (target, value) in atmosphere.iter_mut().zip(words) {
        *target = *value;
    }
    if search && atmosphere[0].partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
        atmosphere[0] = 1013.25;
    }
    atmosphere
}

/// Next heliacal event using geographic/atmospheric numeric sequences.
///
/// Missing geographic words are zero (east/north degrees, metres). Missing
/// atmosphere words are 1013.25 mbar, 15 °C, 40% humidity and zero visibility
/// range in km. Nonpositive search pressure selects 1013.25 mbar. Observer
/// optics and event/flag domains belong to the native numerical solver.
#[allow(clippy::too_many_arguments)]
pub(crate) fn heliacal_ut_words(
    session: &Session,
    jd_ut: f64,
    geopos: &[f64],
    atmo: &[f64],
    observer: &[f64],
    name: &str,
    event: i32,
    flags: i32,
) -> Result<[f64; 3]> {
    session.native_heliacal_ut(
        jd_ut,
        site(geopos)?,
        atmosphere(atmo, true),
        observer,
        name,
        event,
        flags,
    )
}

/// Fifty heliacal phenomena slots from numeric input sequences.
///
/// Input units and missing-word defaults match [`heliacal_ut_words`].
/// Returned angular, magnitude and event-window slots retain the native
/// solver's units, physical domains and separate frame/model clocks.
#[allow(clippy::too_many_arguments)]
pub(crate) fn heliacal_pheno_ut_words(
    session: &Session,
    jd_ut: f64,
    geopos: &[f64],
    atmo: &[f64],
    observer: &[f64],
    name: &str,
    event: i32,
    flags: i32,
) -> Result<[f64; 50]> {
    session.native_heliacal_pheno_ut(
        jd_ut,
        site(geopos)?,
        atmosphere(atmo, true),
        observer,
        name,
        event,
        flags,
    )
}

/// Fixed-instant limiting visual magnitude and ten geometry/detail slots.
///
/// Geographic defaults are zero. Atmospheric defaults are 1013.25 mbar,
/// 15 °C, 50% humidity and zero meteorological range. An explicit pressure
/// remains explicit, including zero. Object, optics, data and solver errors
/// propagate; no timeout or failed calculation becomes a completed miss.
#[allow(clippy::too_many_arguments)]
pub(crate) fn vis_limit_mag_words(
    session: &Session,
    jd_ut: f64,
    geopos: &[f64],
    atmo: &[f64],
    observer: &[f64],
    name: &str,
    flags: i32,
) -> Result<(f64, [f64; 10])> {
    session.vis_limit_mag(
        jd_ut,
        site(geopos)?,
        atmosphere(atmo, false),
        observer,
        name,
        flags,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn private_sequence_entries_are_free_functions_not_session_methods() {
        let _ = heliacal_ut_words;
        let _ = heliacal_pheno_ut_words;
        let _ = vis_limit_mag_words;
    }

    #[test]
    fn numeric_sequences_keep_physical_defaults_and_explicit_words() {
        assert_eq!(site(&[]).unwrap(), [0.0, 0.0, 0.0]);
        assert_eq!(site(&[12.5, 41.9]).unwrap(), [12.5, 41.9, 0.0]);
        assert_eq!(atmosphere(&[], true), [1013.25, 15.0, 40.0, 0.0]);
        assert_eq!(atmosphere(&[], false), [1013.25, 15.0, 50.0, 0.0]);
        assert_eq!(
            atmosphere(&[0.0, 10.0, 75.0, 20.0], true),
            [1013.25, 10.0, 75.0, 20.0]
        );
        assert_eq!(atmosphere(&[0.0], false), [0.0, 15.0, 50.0, 0.0]);
    }
}
