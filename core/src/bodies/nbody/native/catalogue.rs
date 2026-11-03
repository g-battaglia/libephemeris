// SPDX-License-Identifier: AGPL-3.0
//! Physical parameters read from the selected JPL file, never a replacement
//! nominal catalogue or force coefficients fitted to another program's output.
//!
//! JPL classic-DE constants: AU, RE and ASUN are km, CLIGHT is km/s, GMs
//! are AU^3/day^2, EMRAT is the Earth/Moon mass ratio and Jn is unnormalized.
//! Actual values stay attached to the retained reader's identity/coverage.
//! Public configuration facts independently observed for the retained legacy
//! API select SSB coordinates, a Sun-only PPN outer sum and inverse-square
//! nongravitational distance scaling. Poles remain explicit physical inputs;
//! this module does not guess an unspecified orientation convention.

use super::classic_de::{ClassicDe, DeError};
use super::forces::{DistanceLaw, Eih, ForceModel, Nongravitational, RelativisticSources, Zonals};
use super::state::NumericalError;

/// Actual prescribed Newtonian/PPN major-body order: the EMB is replaced by
/// Earth and Moon, not added as a second copy of their combined mass.
pub(super) const MAJOR_TARGETS: [i32; 11] = [10, 1, 2, 399, 301, 4, 5, 6, 7, 8, 9];

/// The physical pole directions in the common ICRF frame at the force epoch.
/// The retained orientation provider, not a data reader, defines their model.
#[derive(Clone, Copy, Debug)]
pub(super) struct HarmonicPoles {
    pub(super) sun: [f64; 3],
    pub(super) earth: [f64; 3],
}

/// Parameters belonging to one selected classic-DE image.
#[derive(Clone, Debug)]
pub(super) struct Catalogue {
    pub(super) astronomical_unit_km: f64,
    pub(super) speed_of_light_au_day: f64,
    pub(super) major_masses: [f64; 11],
    sun_radius: f64,
    earth_radius: f64,
    sun_j2: f64,
    earth_j: [f64; 3],
}

impl Catalogue {
    pub(super) fn load(image: &ClassicDe) -> Result<Self, DeError> {
        Self::read(
            image.header.astronomical_unit_km,
            image.header.earth_moon_mass_ratio,
            |name| image.constant(name),
        )
    }

    fn read(
        astronomical_unit_km: f64,
        earth_moon_ratio: f64,
        read: impl Fn(&str) -> Result<f64, DeError>,
    ) -> Result<Self, DeError> {
        if ![astronomical_unit_km, earth_moon_ratio]
            .into_iter()
            .all(|value| value.is_finite() && value > 0.0)
        {
            return Err(DeError::Malformed("force units or Earth/Moon mass ratio"));
        }
        let positive = |name| -> Result<f64, DeError> {
            let value = read(name)?;
            if value.is_finite() && value > 0.0 {
                Ok(value)
            } else {
                Err(DeError::Malformed("nonpositive dimensional force constant"))
            }
        };
        let mass = |name| -> Result<f64, DeError> {
            let value = read(name)?;
            if value.is_finite() && value >= 0.0 {
                Ok(value)
            } else {
                Err(DeError::Malformed("negative or nonfinite GM"))
            }
        };
        let dimensionless = |name| -> Result<f64, DeError> {
            let value = read(name)?;
            if value.is_finite() {
                Ok(value)
            } else {
                Err(DeError::Malformed("nonfinite zonal coefficient"))
            }
        };
        let earth_moon_mass = mass("GMB")?;
        let moon_mass = earth_moon_mass / (1.0 + earth_moon_ratio);
        let earth_mass = earth_moon_mass - moon_mass;
        let major_masses = [
            mass("GMS")?,
            mass("GM1")?,
            mass("GM2")?,
            earth_mass,
            moon_mass,
            mass("GM4")?,
            mass("GM5")?,
            mass("GM6")?,
            mass("GM7")?,
            mass("GM8")?,
            mass("GM9")?,
        ];
        let speed_of_light_au_day = (positive("CLIGHT")? / astronomical_unit_km) * 86400.0;
        let sun_radius = positive("ASUN")? / astronomical_unit_km;
        let earth_radius = positive("RE")? / astronomical_unit_km;
        if ![speed_of_light_au_day, sun_radius, earth_radius]
            .into_iter()
            .all(|value| value.is_finite() && value > 0.0)
        {
            return Err(DeError::Malformed("unrepresentable force unit conversion"));
        }
        Ok(Self {
            astronomical_unit_km,
            speed_of_light_au_day,
            major_masses,
            sun_radius,
            earth_radius,
            sun_j2: dimensionless("J2SUN")?,
            earth_j: [
                dimensionless("J2E")?,
                dimensionless("J3E")?,
                dimensionless("J4E")?,
            ],
        })
    }

    /// Assemble all selected families without removing any source on an error.
    /// The API's alpha=1,r0=1,m=2,k=0 law reduces exactly to r^-2; its n is
    /// irrelevant. The PPN outer Sun selection does not truncate inner potentials.
    /// Pole vectors are mandatory here until an independently specified
    /// orientation provider supplies them at each actual evaluation epoch.
    pub(super) fn force_model(
        &self,
        poles: HarmonicPoles,
        nongravitational: Option<[f64; 3]>,
    ) -> Result<ForceModel, NumericalError> {
        if nongravitational.is_some_and(|parameters| !parameters.into_iter().all(f64::is_finite)) {
            return Err(NumericalError::NonFinite);
        }
        Ok(ForceModel {
            relativity: Some(Eih {
                speed_of_light: self.speed_of_light_au_day,
                outer: RelativisticSources::Sun,
            }),
            zonals: vec![
                Zonals::new(10, self.sun_radius, poles.sun, [self.sun_j2, 0.0, 0.0])?,
                Zonals::new(399, self.earth_radius, poles.earth, self.earth_j)?,
            ],
            nongravitational: nongravitational.map(|parameters| Nongravitational {
                parameters,
                law: DistanceLaw::InverseSquare,
            }),
        })
    }
}

#[cfg(test)]
#[path = "catalogue_tests.rs"]
mod tests;
