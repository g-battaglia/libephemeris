//! Supplied chart positions and geometric high-latitude sect.

use super::Session;
use crate::arabic_parts::{self, ArabicError, ArabicInputs};
use crate::errors::Error;
use std::collections::HashMap;
use std::fmt;

/// Invalid supplied lots map, or a native clock/horizon calculation failure.
#[derive(Debug)]
pub enum ArabicPartsError {
    /// Required longitude keys absent from the supplied map.
    MissingPositions(Vec<&'static str>),
    /// A non-finite position or invalid lots input.
    Input(ArabicError),
    /// Clock, observer or horizontal-coordinate calculation failed.
    Native(Error),
}

impl fmt::Display for ArabicPartsError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingPositions(keys) => {
                write!(formatter, "missing required positions: {}", keys.join(", "))
            }
            Self::Input(error) => error.fmt(formatter),
            Self::Native(error) => error.fmt(formatter),
        }
    }
}

impl std::error::Error for ArabicPartsError {}

impl From<ArabicError> for ArabicPartsError {
    fn from(error: ArabicError) -> Self {
        Self::Input(error)
    }
}

impl From<Error> for ArabicPartsError {
    fn from(error: Error) -> Self {
        Self::Native(error)
    }
}

impl Session {
    /// Determine sect from supplied solar and Ascendant ecliptic longitudes.
    ///
    /// Angles and geographic coordinates are degrees; `jd` is UT1. Complete
    /// date/location inputs with absolute latitude greater than 60 degrees use
    /// the actual unrefracted [`Self::azalt`] geometry. Otherwise the open
    /// Ascendant-to-Descendant arc defines night. Horizon boundaries count as
    /// day. Native horizon/clock failures and invalid consumed angles propagate.
    pub fn is_day_chart(
        &self,
        sun_lon: f64,
        asc: f64,
        jd: Option<f64>,
        geo_lat: Option<f64>,
        geo_lon: Option<f64>,
        sun_lat: f64,
    ) -> Result<bool, ArabicPartsError> {
        match (jd, geo_lat, geo_lon) {
            (Some(jd), Some(lat), Some(lon))
                if lat.abs() > arabic_parts::EXTREME_LATITUDE_THRESHOLD =>
            {
                Ok(self
                    .azalt(
                        jd,
                        crate::events::ECL2HOR,
                        lon,
                        lat,
                        0.0,
                        0.0,
                        0.0,
                        sun_lon,
                        sun_lat,
                    )?
                    .true_altitude
                    >= 0.0)
            }
            _ => Ok(arabic_parts::is_day_chart_2d(sun_lon, asc)?),
        }
    }

    /// Calculate four lots from supplied ecliptic positions in degrees.
    ///
    /// Requires Asc, Sun, Moon, Mercury and Venus; Sun_lat defaults to zero.
    /// At a complete observer/JD setting with |latitude|>60 degrees, sect is
    /// the sign of the Sun's geometric altitude using [`Self::azalt`]; otherwise
    /// the open Ascendant-to-Descendant ecliptic arc defines night. Horizon
    /// boundaries count as day. Output longitudes are wrapped to [0,360).
    /// Missing positions, non-finite angles and clock/projection errors fail.
    pub fn calc_all_arabic_parts(
        &self,
        positions: &HashMap<String, f64>,
        jd: Option<f64>,
        geo_lat: Option<f64>,
        geo_lon: Option<f64>,
    ) -> Result<HashMap<String, f64>, ArabicPartsError> {
        let missing: Vec<_> = ["Asc", "Sun", "Moon", "Mercury", "Venus"]
            .into_iter()
            .filter(|name| !positions.contains_key(*name))
            .collect();
        if !missing.is_empty() {
            return Err(ArabicPartsError::MissingPositions(missing));
        }
        let inputs = ArabicInputs::new(
            positions["Asc"],
            positions["Sun"],
            positions["Moon"],
            positions["Mercury"],
            positions["Venus"],
            positions.get("Sun_lat").copied().unwrap_or(0.0),
        )?;
        let diurnal =
            self.is_day_chart(inputs.sun, inputs.asc, jd, geo_lat, geo_lon, inputs.sun_lat)?;
        Ok(HashMap::from([
            (
                "Pars_Fortunae".into(),
                arabic_parts::part_of_fortune(inputs.asc, inputs.sun, inputs.moon, diurnal)?,
            ),
            (
                "Pars_Spiritus".into(),
                arabic_parts::part_of_spirit(inputs.asc, inputs.sun, inputs.moon, diurnal)?,
            ),
            (
                "Pars_Amoris".into(),
                arabic_parts::part_of_love(inputs.asc, inputs.venus, inputs.sun)?,
            ),
            (
                "Pars_Fidei".into(),
                arabic_parts::part_of_faith(inputs.asc, inputs.mercury, inputs.moon)?,
            ),
        ]))
    }
}

impl crate::config::Runtime {
    /// Sect from supplied angles using this runtime's configured UT1 clock.
    ///
    /// Arguments, units and errors follow [`Session::is_day_chart`]. This
    /// conversion does not open a positional ephemeris or acquire planets.
    pub fn is_day_chart(
        &self,
        sun_lon: f64,
        asc: f64,
        jd: Option<f64>,
        geo_lat: Option<f64>,
        geo_lon: Option<f64>,
        sun_lat: f64,
    ) -> Result<bool, ArabicPartsError> {
        let context = self.effective_context()?;
        Session::from_context(context).is_day_chart(sun_lon, asc, jd, geo_lat, geo_lon, sun_lat)
    }

    /// Calculate lots using this runtime's model clock without opening a
    /// positional ephemeris. Inputs and units follow [`Session::calc_all_arabic_parts`].
    /// Supplied positions define the chart; native horizon geometry supplies
    /// high-latitude sect. Configuration and calculation errors propagate.
    pub fn calc_all_arabic_parts(
        &self,
        positions: &HashMap<String, f64>,
        jd: Option<f64>,
        geo_lat: Option<f64>,
        geo_lon: Option<f64>,
    ) -> Result<HashMap<String, f64>, ArabicPartsError> {
        let context = self.effective_context()?;
        Session::from_context(context).calc_all_arabic_parts(positions, jd, geo_lat, geo_lon)
    }
}
