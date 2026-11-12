//! Horizon geometry with explicit session clocks and atmospheric settings.
use super::*;
pub use crate::events::{ECL2HOR, EQU2HOR, HOR2ECL, HOR2EQU};

/// Horizontal coordinates; azimuth is measured from South westward.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Horizontal {
    /// South-westward azimuth in degrees; North-eastward is `(azimuth+180)%360`.
    pub azimuth: f64,
    /// Geometric altitude in degrees.
    pub true_altitude: f64,
    /// Refracted altitude in degrees.
    pub apparent_altitude: f64,
}
/// Celestial coordinates from a reverse horizon projection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Celestial {
    /// Ecliptic longitude or equatorial RA in degrees, selected by the flag.
    pub longitude_or_ra: f64,
    /// Ecliptic latitude or equatorial declination in degrees.
    pub latitude_or_dec: f64,
}
impl Session {
    /// Horizon projection at UT. Geographic east longitude/latitude are degrees,
    /// altitude metres, pressure mbar (zero estimates it), temperature Celsius.
    /// Input celestial angles are degrees; distance is engine-ignored. ECL2HOR
    /// uses the independent bundled-timescale frame epoch, not model/body TT;
    /// sidereal terms use the explicitly configured model Delta-T. No engine or
    /// data fallback is introduced. Non-finite projection inputs/results refuse.
    #[allow(clippy::too_many_arguments)]
    pub fn azalt(
        &self,
        jd_ut: f64,
        flag: i32,
        lon: f64,
        lat: f64,
        alt: f64,
        pressure: f64,
        temp: f64,
        lon_or_ra: f64,
        lat_or_dec: f64,
    ) -> Result<Horizontal, Error> {
        let dt = self.sidereal_delta(jd_ut)?;
        let delta = |ut| self.deltat(ut);
        let lapse = self
            .ctx
            .lapse_rate_override()
            .unwrap_or(crate::events::RISE_SET_DIP_LAPSE_RATE_K_M);
        let (azimuth, true_altitude, apparent_altitude) = crate::events::azalt_deg(
            jd_ut,
            flag,
            lon,
            lat,
            alt,
            pressure,
            temp,
            (lon_or_ra, lat_or_dec, 1.0),
            &delta,
            dt.offset_lo_deg(),
            dt.offset_hi_deg(),
            lapse,
        )?;
        Ok(Horizontal {
            azimuth,
            true_altitude,
            apparent_altitude,
        })
    }

    /// Reverse projection of South-westward azimuth and geometric altitude
    /// (degrees) at UT; geographic height is accepted but ignored by the kernel.
    /// HOR2EQU selects equatorial RA/Dec, otherwise ecliptic longitude/latitude.
    #[allow(clippy::too_many_arguments)]
    pub fn azalt_rev(
        &self,
        jd_ut: f64,
        flag: i32,
        lon: f64,
        lat: f64,
        alt: f64,
        azimuth: f64,
        true_altitude: f64,
    ) -> Result<Celestial, Error> {
        let dt = self.sidereal_delta(jd_ut)?;
        let delta = |ut| self.deltat(ut);
        let (longitude_or_ra, latitude_or_dec) = crate::events::azalt_rev_deg(
            jd_ut,
            flag,
            lon,
            lat,
            alt,
            azimuth,
            true_altitude,
            &delta,
            dt.offset_lo_deg(),
            dt.offset_hi_deg(),
        )?;
        Ok(Celestial {
            longitude_or_ra,
            latitude_or_dec,
        })
    }
}
