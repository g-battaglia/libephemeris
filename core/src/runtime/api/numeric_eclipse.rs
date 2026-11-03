//! Numeric event targets and geographic-sequence composition for public events.

use super::{Session, SolarDurationError, StellarOccultationAcquisition};
use crate::errors::Result;
use crate::events;

/// Owned planetary identifier or fixed-star catalog key for an event.
#[derive(Clone, Debug)]
pub enum EventTarget {
    /// Numerical body identifier; registered and catalog identifiers retain
    /// the same position provider as an ordinary calculation.
    Body(i32),
    /// Fixed-star name or catalog query accepted by the native catalog.
    Star(String),
}

/// Geographic longitude east, latitude north (degrees) and height (metres).
///
/// At least three numeric words are required; extra words
/// are ignored. Coordinate domains are checked by the consuming calculation.
/// A short sequence returns `InvalidEclipseSearch` with entry `geopos`.
pub fn geographic_site(words: &[f64]) -> Result<[f64; 3]> {
    events::parse_geopos(words, "geopos", events::GEOPOS_SEQUENCE_MESSAGE)
        .map(|(lon, lat, height)| [lon, lat, height])
}

struct CatalogTarget<'a> {
    session: &'a Session,
    target: &'a EventTarget,
}

impl StellarOccultationAcquisition for CatalogTarget<'_> {
    fn position(&self, jd: f64, flags: i32) -> Result<[f64; 6]> {
        match self.target {
            EventTarget::Body(body) => self
                .session
                .calc_ut(jd, *body, flags)
                .map(|value| value.values),
            EventTarget::Star(key) => self
                .session
                .fixstar_ut(key, jd, flags)
                .map(|value| value.values),
        }
    }
}

#[doc(hidden)]
impl Session {
    /// Fixed-instant solar circumstances from a geographic numeric sequence.
    ///
    /// UT1 Julian day, geographic degrees/metres and the calculation bitmask
    /// retain [`Self::sol_eclipse_how`]'s numerical model and error domain.
    pub fn sol_eclipse_how_words(
        &self,
        jd: f64,
        site: &[f64],
        flags: i32,
    ) -> Result<(i32, [f64; 20])> {
        self.sol_eclipse_how(jd, geographic_site(site)?, flags)
    }

    /// Comprehensive fixed-instant solar circumstances at a numeric site.
    ///
    /// Units, contact windows and errors follow [`Self::sol_eclipse_how_details`].
    pub fn sol_eclipse_how_details_words(
        &self,
        jd: f64,
        site: &[f64],
        flags: i32,
    ) -> Result<events::EclipseHowDetails> {
        self.sol_eclipse_how_details(jd, geographic_site(site)?, flags)
    }

    /// Visible covered-diameter fraction at UT1 JD and a numeric site.
    ///
    /// Geographic words are degrees/metres; flags are retained as an API
    /// argument, while apparent local rays use their physical observer model.
    pub fn sol_eclipse_magnitude_words(&self, jd: f64, site: &[f64], _flags: i32) -> Result<f64> {
        self.sol_eclipse_magnitude_at_loc(jd, geographic_site(site)?)
    }

    /// Visible covered-area fraction at UT1 JD and a numeric site.
    ///
    /// Input units, visibility domain and errors follow the diameter query.
    pub fn sol_eclipse_obscuration_words(&self, jd: f64, site: &[f64], _flags: i32) -> Result<f64> {
        self.sol_eclipse_obscuration_at_loc(jd, geographic_site(site)?)
    }

    /// Solar maximum near UT1 JD, globally or at a fully specified site.
    ///
    /// Latitude/longitude are degrees and height is metres. A site is local
    /// only if both angular coordinates are supplied; supplying just one is an
    /// invalid geographic site. The bracket half-width
    /// is in days; returned maximum is UT1 JD and magnitude is dimensionless.
    /// Flags are accepted but do not alter the local apparent-ray model.
    #[allow(clippy::too_many_arguments)]
    pub fn sol_eclipse_maximum_words(
        &self,
        jd: f64,
        lat: Option<f64>,
        lon: Option<f64>,
        altitude: f64,
        search_range: f64,
        _flags: i32,
    ) -> Result<(f64, f64)> {
        if lat.is_some() != lon.is_some() {
            return Err(crate::errors::Error::InvalidEclipseSearch {
                entry: "geopos",
                detail: "Both lat and lon must be provided, or neither for a global maximum"
                    .to_owned(),
            });
        }
        let site = lat.zip(lon).map(|(lat, lon)| [lon, lat, altitude]);
        self.sol_eclipse_max_time(jd, site, search_range)
    }

    /// Next visible solar eclipse at a numeric geographic site.
    ///
    /// Inputs and output slots follow [`Self::sol_eclipse_when_loc`]; no
    /// observer setting is published to the caller's context.
    pub fn sol_eclipse_when_loc_words(
        &self,
        jd: f64,
        site: &[f64],
        flags: i32,
        backwards: bool,
    ) -> Result<(i32, [f64; 10], [f64; 20])> {
        self.sol_eclipse_when_loc(jd, geographic_site(site)?, flags, backwards)
    }

    /// Lunar circumstances at UT1 JD and a numeric geographic site.
    ///
    /// Units and geometry follow [`Self::lun_eclipse_how`].
    pub fn lun_eclipse_how_words(
        &self,
        jd: f64,
        site: &[f64],
        flags: i32,
    ) -> Result<(i32, [f64; 20])> {
        self.lun_eclipse_how(jd, geographic_site(site)?, flags)
    }

    /// Next visible lunar eclipse at a numeric geographic site.
    ///
    /// Units and phase/horizon slots follow [`Self::lun_eclipse_when_loc`].
    pub fn lun_eclipse_when_loc_words(
        &self,
        jd: f64,
        site: &[f64],
        flags: i32,
        backwards: bool,
    ) -> Result<(i32, [f64; 10], [f64; 20])> {
        self.lun_eclipse_when_loc(jd, geographic_site(site)?, flags, backwards)
    }

    /// Lunar occultation geography for a numeric body or fixed-star key.
    ///
    /// Numeric catalog identifiers use their normal public body position;
    /// names use the native fixed-star catalog. Both share the same shadow,
    /// horizon and observer geometry. Returned units and errors follow
    /// [`Self::lun_occult_where`] and [`Self::stellar_lun_occult_where`].
    pub fn lun_occult_where_target(
        &self,
        jd: f64,
        target: &EventTarget,
        flags: i32,
    ) -> Result<(i32, [f64; 10], [f64; 20])> {
        match target {
            EventTarget::Body(body) if crate::stars::canonical_star_name(*body).is_none() => {
                self.lun_occult_where(jd, *body, flags)
            }
            _ => self.stellar_lun_occult_where(
                jd,
                flags,
                &CatalogTarget {
                    session: self,
                    target,
                },
            ),
        }
    }

    /// Global lunar occultation search for a body or fixed-star key.
    ///
    /// UT1 days, calculation/type bitmasks and `backwards` bit0/ECL_ONE_TRY
    /// retain the native controllers' bounded search and contact domains.
    /// A completed miss and a failed calculation remain distinct outcomes.
    pub fn lun_occult_when_glob_target(
        &self,
        jd: f64,
        target: &EventTarget,
        flags: i32,
        mask: i32,
        backwards: i32,
    ) -> Result<(i32, [f64; 10])> {
        let (catalog_id, label) = match target {
            EventTarget::Body(body) if crate::stars::canonical_star_name(*body).is_none() => {
                return self.lun_occult_when_glob(jd, *body, flags, mask, backwards);
            }
            EventTarget::Body(body) => (Some(*body), body.to_string()),
            EventTarget::Star(key) => (None, key.clone()),
        };
        self.stellar_lun_occult_when_glob(
            jd,
            catalog_id,
            &label,
            flags,
            mask,
            backwards,
            &CatalogTarget {
                session: self,
                target,
            },
        )
    }

    /// Local lunar occultation search at east/north degrees and metres.
    ///
    /// Integer catalog identifiers resolve to their native canonical key for
    /// the local stellar route. Other bodies retain their numerical provider.
    /// Geographic length, data, position and solver errors propagate.
    pub fn lun_occult_when_loc_target(
        &self,
        jd: f64,
        target: &EventTarget,
        site: &[f64],
        flags: i32,
        backwards: i32,
    ) -> Result<(i32, [f64; 10], [f64; 20])> {
        let site = geographic_site(site)?;
        let key = match target {
            EventTarget::Body(body) => match crate::stars::canonical_star_name(*body) {
                Some(key) => key.to_owned(),
                None => return self.lun_occult_when_loc(jd, *body, site, flags, backwards),
            },
            EventTarget::Star(key) => key.clone(),
        };
        let target = EventTarget::Star(key.clone());
        self.stellar_lun_occult_when_loc(
            jd,
            &key,
            site,
            flags,
            backwards,
            &CatalogTarget {
                session: self,
                target: &target,
            },
        )
    }

    /// Global planetary occultation search with a zero sentinel for a star.
    ///
    /// Inputs use UT1 days and native body identifiers. Flags are retained in
    /// the public interface; the global apparent-place model determines its
    /// own position flags. Zero target with a nonempty key selects a star.
    #[allow(clippy::too_many_arguments)]
    pub fn planet_occult_when_glob_words(
        &self,
        jd: f64,
        occulting: i32,
        occulted: i32,
        star: &str,
        _flags: i32,
        direction: i32,
    ) -> Result<(i32, [f64; 10])> {
        self.planet_occult_when_glob(
            jd,
            occulting,
            (occulted != 0).then_some(occulted),
            star,
            direction,
        )
    }

    /// Local planetary occultation from explicit north/east degrees/metres.
    ///
    /// The zero target sentinel selects a star; phase and horizon algorithms,
    /// their bounded search and all failures are the native controller's.
    #[allow(clippy::too_many_arguments)]
    pub fn planet_occult_when_loc_words(
        &self,
        jd: f64,
        occulting: i32,
        occulted: i32,
        star: &str,
        lat: f64,
        lon: f64,
        altitude: f64,
        _flags: i32,
    ) -> Result<(i32, [f64; 10], [f64; 20])> {
        self.planet_occult_when_loc(
            jd,
            occulting,
            (occulted != 0).then_some(occulted),
            star,
            [lon, lat, altitude],
        )
    }

    /// Rise/set or meridian transit of a body or named fixed star.
    ///
    /// UT1 days, geographic degrees/metres, pressure mbar, temperature Celsius
    /// and horizon degrees use one native horizon solver. Numeric catalog ids
    /// keep the ordinary body route; only a name selects the stellar route.
    #[allow(clippy::too_many_arguments)]
    pub fn rise_trans_target(
        &self,
        jd: f64,
        target: &EventTarget,
        rsmi: i32,
        site: &[f64],
        pressure: f64,
        temperature: f64,
        horizon: f64,
        flags: i32,
    ) -> Result<(i32, [f64; 10])> {
        let site = geographic_site(site)?;
        match target {
            EventTarget::Body(body) => self.rise_trans_true_hor(
                jd,
                *body,
                rsmi,
                site,
                pressure,
                temperature,
                horizon,
                flags,
            ),
            EventTarget::Star(key) => self.star_rise_trans_true_hor(
                key,
                jd,
                rsmi,
                site,
                pressure,
                temperature,
                horizon,
                flags,
            ),
        }
    }

    /// Solar central-phase duration in minutes at a fixed UT1 maximum.
    ///
    /// This plain-error boundary preserves the actual physical failure from
    /// the observer/reader acquisition; no failure is converted into a miss.
    pub fn solar_eclipse_duration_minutes(&self, jd: f64, flags: i32) -> Result<f64> {
        self.solar_eclipse_duration_at_ut(jd, flags)
            .map_err(|error| match error {
                SolarDurationError::Where(error)
                | SolarDurationError::Native(error)
                | SolarDurationError::Reader { error, .. } => error,
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::errors::Error;

    #[test]
    fn geographic_sequence_requires_three_words_and_ignores_extras() {
        assert_eq!(
            geographic_site(&[12.5, 41.9, 30.0, 99.0]).unwrap(),
            [12.5, 41.9, 30.0]
        );
        assert!(matches!(
            geographic_site(&[12.5, 41.9]),
            Err(Error::InvalidEclipseSearch {
                entry: "geopos",
                ..
            })
        ));
    }
}
