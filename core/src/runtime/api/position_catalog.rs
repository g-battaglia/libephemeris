//! Catalog stars, terrestrial angles and cached lots share the public dispatcher.

use super::{InputEpoch, Position, PositionRequest, Session, completed_position};
use crate::constants::{
    ANGLE_OFFSET, ARABIC_OFFSET, PARS_AMORIS, PARS_FIDEI, PARS_FORTUNAE, PARS_SPIRITUS,
};
use crate::errors::Error;
use crate::positions;

impl Session {
    /// Dispatch nonplanetary identifiers before acquiring a planetary image.
    /// The requested epoch is TT or UT1; cached lots retain their cached inputs.
    pub(super) fn catalog_position(
        &self,
        request: PositionRequest,
    ) -> Result<Option<Position>, Error> {
        let PositionRequest {
            epoch,
            body,
            input_flags,
        } = request;
        let (jd, is_ut) = epoch.into_entry();
        // Catalog echo rules differ between TT and UT1. They do not replace
        // the raw input flags passed to the star calculation below.
        let returned_flags = match epoch {
            InputEpoch::Ut1(_) => {
                positions::normalize_calc_flags(input_flags, false)
                    | positions::implied_retflag_bits(input_flags)
            }
            InputEpoch::Tt(_) => input_flags | positions::implied_retflag_bits(input_flags),
        };
        let longitude = if crate::stars::star_entry_by_id(body).is_some() {
            let star = self.fixed_star_by_id(body, jd, input_flags, is_ut)?;
            return Ok(Some(completed_position(
                star.values,
                returned_flags,
                self.diagnostics.join("; "),
            )));
        } else if (ANGLE_OFFSET..ARABIC_OFFSET).contains(&body) {
            let observer = self.ctx.observer().ok_or(Error::ConfigurationMissing {
                missing_config: "observer_location",
                suggestion: "set an observer before calculating an angle identifier",
            })?;
            let ut = match epoch {
                InputEpoch::Ut1(jd) => jd,
                InputEpoch::Tt(jd) => crate::time::default_timescale_ut1_from_tt(jd)?,
            };
            self.numeric_angle_value(body, ut, observer.lat_deg(), observer.lon_deg())?
        } else if (ARABIC_OFFSET..ARABIC_OFFSET + 100).contains(&body) {
            self.cached_part_longitude(body)?
        } else {
            return Ok(None);
        };
        Ok(Some(completed_position(
            [longitude, 0.0, 0.0, 0.0, 0.0, 0.0],
            returned_flags,
            self.diagnostics.join("; "),
        )))
    }

    fn cached_part_longitude(&self, body: i32) -> Result<f64, Error> {
        use crate::arabic_parts;
        let cache = self.ctx.angles_cache();
        if cache.is_empty() {
            return Err(Error::ConfigurationMissing {
                missing_config: "angles_cache",
                suggestion: "calculate or set the context's angle cache first",
            });
        }
        let asc = cache
            .get("Asc")
            .or_else(|| cache.get("Ascendant"))
            .copied()
            .unwrap_or(0.0);
        let value = |name| cache.get(name).copied().unwrap_or(0.0);
        let result = match body {
            PARS_FORTUNAE | PARS_SPIRITUS => {
                let day = arabic_parts::is_day_chart_2d(value("Sun"), asc).map_err(part_error)?;
                if body == PARS_FORTUNAE {
                    arabic_parts::part_of_fortune(asc, value("Sun"), value("Moon"), day)
                } else {
                    arabic_parts::part_of_spirit(asc, value("Sun"), value("Moon"), day)
                }
            }
            PARS_AMORIS => arabic_parts::part_of_love(asc, value("Venus"), value("Sun")),
            PARS_FIDEI => arabic_parts::part_of_faith(asc, value("Mercury"), value("Moon")),
            _ => Ok(0.0),
        };
        result.map_err(part_error)
    }
}

fn part_error(error: crate::arabic_parts::ArabicError) -> Error {
    match error {
        crate::arabic_parts::ArabicError::NonFinite { value, .. } => {
            Error::NonFiniteArgument { value }
        }
        _ => Error::ConfigurationMissing {
            missing_config: "angles_cache",
            suggestion: "set valid finite ecliptic longitudes in the angle cache",
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{FLG_BARYCTR, FLG_MOSEPH, MARS};
    use crate::context::EphemerisContext;

    #[test]
    fn ordinary_body_declines_without_validating_epoch_or_flags() {
        let session = Session::from_context(EphemerisContext::new());
        for epoch in [InputEpoch::Tt(f64::NAN), InputEpoch::Ut1(f64::NAN)] {
            let result = session
                .catalog_position(PositionRequest {
                    epoch,
                    body: MARS,
                    input_flags: FLG_BARYCTR | FLG_MOSEPH,
                })
                .unwrap();
            assert!(result.is_none());
        }
    }

    #[test]
    fn angle_requires_observer_before_consuming_tt_or_ut1() {
        let session = Session::from_context(EphemerisContext::new());
        for epoch in [InputEpoch::Tt(f64::NAN), InputEpoch::Ut1(f64::NAN)] {
            let result = session.catalog_position(PositionRequest {
                epoch,
                body: ANGLE_OFFSET,
                input_flags: 0,
            });
            assert!(matches!(
                result,
                Err(Error::ConfigurationMissing {
                    missing_config: "observer_location",
                    ..
                })
            ));
        }
    }
}
