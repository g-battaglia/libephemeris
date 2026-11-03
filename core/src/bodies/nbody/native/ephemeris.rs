// SPDX-License-Identifier: AGPL-3.0
//! One retained selected-data owner for common-frame DE/SPK force snapshots.
//!
//! Major-body parameters/kinematics come from the same classic-DE image.
//! Small-body GMs use JPL physical-ID/header-label identities, including the
//! reserved SB441 KBO aliases, rather than file-order mass guesses.
//! SPK km/km-s vectors are translated from their declared major-body center
//! at the same two-part TDB epoch, then converted to the selected AU/day units.
//! No force family is removed when data, metadata or evaluation fails.

use super::catalogue::{Catalogue, MAJOR_TARGETS};
use super::classic_de::{ClassicDe, DeError, Kinematics};
use super::forces::Source;
use super::newtonian::PointMass;
use crate::errors::Error;
use crate::spk::{ChebState, ChebyshevFile, J2000_JD_TDB, SECONDS_PER_DAY};
use std::path::Path;

#[derive(Debug)]
pub(super) enum EphemerisError {
    De(DeError),
    Spk(Error),
    InvalidParticle(usize),
    InvalidSelection(&'static str),
}
impl From<DeError> for EphemerisError {
    fn from(error: DeError) -> Self {
        Self::De(error)
    }
}
impl From<Error> for EphemerisError {
    fn from(error: Error) -> Self {
        Self::Spk(error)
    }
}
impl std::fmt::Display for EphemerisError {
    fn fmt(&self, output: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::De(error) => write!(output, "{error}"),
            Self::Spk(error) => write!(output, "{error}"),
            Self::InvalidParticle(index) => {
                write!(output, "native force particle index {index} is unavailable")
            }
            Self::InvalidSelection(detail) => {
                write!(output, "invalid native force data selection: {detail}")
            }
        }
    }
}
impl std::error::Error for EphemerisError {}

#[derive(Clone, Copy)]
struct Asteroid {
    target: i32,
    mass: f64,
}

/// Immutable data/parameter owner; only bounded file cursors/coefficient caches
/// have internal mutexes. Trajectories keep their own state and frame scratch.
pub(super) struct Ephemeris {
    planets: ClassicDe,
    asteroids: Option<ChebyshevFile>,
    small_bodies: Vec<Asteroid>,
    pub(super) catalogue: Catalogue,
}

impl Ephemeris {
    pub(super) fn open(planets: &Path, asteroids: Option<&Path>) -> Result<Self, EphemerisError> {
        Self::from_planets(ClassicDe::open(planets)?, asteroids)
    }

    fn from_planets(planets: ClassicDe, asteroids: Option<&Path>) -> Result<Self, EphemerisError> {
        let catalogue = Catalogue::load(&planets)?;
        let asteroids = asteroids.map(ChebyshevFile::open).transpose()?;
        let mut small_bodies = Vec::new();
        if let Some(image) = &asteroids {
            for target in image.targets() {
                let mass = asteroid_mass(&planets, target)?;
                if mass < 0.0 || !mass.is_finite() {
                    return Err(DeError::Malformed("asteroid GM").into());
                }
                small_bodies.push(Asteroid { target, mass });
            }
            if small_bodies.is_empty() {
                return Err(EphemerisError::InvalidSelection(
                    "empty selected asteroid force file",
                ));
            }
        }
        Ok(Self {
            planets,
            asteroids,
            small_bodies,
            catalogue,
        })
    }

    pub(super) fn particle_count(&self) -> usize {
        MAJOR_TARGETS.len() + self.small_bodies.len()
    }

    /// Physical identity at a compatibility ordinal, without evaluating data.
    pub(super) fn physical_id(&self, index: usize) -> Option<i32> {
        MAJOR_TARGETS.get(index).copied().or_else(|| {
            index
                .checked_sub(MAJOR_TARGETS.len())
                .and_then(|offset| self.small_bodies.get(offset).map(|body| body.target))
        })
    }

    /// Legacy state-API order: Sun, Mercury, Venus, Earth, Moon, Mars through
    /// Pluto, then unique force-file targets in file occurrence order. The
    /// actual physical identifier is retained separately from this ordinal.
    pub(super) fn particle(&self, index: usize, epoch: [f64; 2]) -> Result<Source, EphemerisError> {
        if let Some(&target) = MAJOR_TARGETS.get(index) {
            let value = self.planets.state(target, epoch)?;
            return Ok(major(target, self.catalogue.major_masses[index], value));
        }
        let asteroid = self
            .small_bodies
            .get(index - MAJOR_TARGETS.len())
            .copied()
            .ok_or(EphemerisError::InvalidParticle(index))?;
        self.asteroid(asteroid, epoch, None)
    }

    /// Read a selected heliocentric asteroid seed from this same bounded owner.
    ///
    /// The relative state remains in SPK km/km-s, before any SSB translation.
    /// Curated initialization requires an actual Sun-centered descriptor;
    /// accepting a different center here would change the public seed model.
    /// A target absent from the selected image has no stored seed, not an error.
    pub(super) fn heliocentric_seed(
        &self,
        target: i32,
        epoch: [f64; 2],
    ) -> crate::errors::Result<Option<ChebState>> {
        if !self.small_bodies.iter().any(|body| body.target == target) {
            return Ok(None);
        }
        let file = self
            .asteroids
            .as_ref()
            .expect("selected target retains its file");
        let (center, relative) = file.state(target, epoch[0], epoch[1])?;
        if center != crate::spk::CENTER_SUN {
            return Err(Error::SpkAmbiguousCenter { target });
        }
        Ok(Some(relative))
    }

    /// Fill caller-owned reusable scratch with a consistent complete frame.
    /// No row is omitted or published as a successful partial result. The
    /// caller uses the vector only on Ok; its capacity may be reused on error.
    pub(super) fn frame(
        &self,
        epoch: [f64; 2],
        destination: &mut Vec<Source>,
    ) -> Result<(), EphemerisError> {
        let values = self.planets.states(MAJOR_TARGETS, epoch)?;
        destination.clear();
        for ((target, mass), value) in MAJOR_TARGETS
            .into_iter()
            .zip(self.catalogue.major_masses)
            .zip(values)
        {
            destination.push(major(target, mass, value));
        }
        for &asteroid in &self.small_bodies {
            destination.push(self.asteroid(asteroid, epoch, Some(&values))?);
        }
        Ok(())
    }

    fn asteroid(
        &self,
        asteroid: Asteroid,
        epoch: [f64; 2],
        majors: Option<&[Kinematics; 11]>,
    ) -> Result<Source, EphemerisError> {
        let file = self
            .asteroids
            .as_ref()
            .expect("asteroid identity belongs to retained file");
        let (center, relative) = file.state(asteroid.target, epoch[0], epoch[1])?;
        let center = match majors.and_then(|values| {
            MAJOR_TARGETS
                .iter()
                .position(|&target| target == center)
                .map(|index| values[index])
        }) {
            Some(value) => value,
            // SSB0 and the EMB3 also have exact definitions in this same
            // selected image. An unknown center is an explicit DE error,
            // never an assumed Sun or a cross-dataset replacement.
            None => self.planets.state(center, epoch)?,
        };
        let au = self.catalogue.astronomical_unit_km;
        let position = std::array::from_fn(|axis| {
            center.state.position[axis] + relative.position_km[axis] / au
        });
        let velocity = std::array::from_fn(|axis| {
            center.state.velocity[axis] + relative.velocity_km_s[axis] * (SECONDS_PER_DAY / au)
        });
        if !position.into_iter().chain(velocity).all(f64::is_finite) {
            return Err(DeError::Malformed("nonfinite composed asteroid state").into());
        }
        Ok(Source {
            point: PointMass {
                body: asteroid.target,
                position,
                gravitational_parameter: asteroid.mass,
            },
            velocity,
            post_newtonian_acceleration: None,
        })
    }

    /// Common declared TDB coverage of all selected physical perturbers, in
    /// Julian days. Disjoint windows remain disjoint. Coefficient corruption
    /// still fails at evaluation; a coverage label is not a numerical warranty.
    pub(super) fn coverage(&self) -> Vec<(f64, f64)> {
        let mut spans = vec![(self.planets.header.start, self.planets.header.end)];
        if let Some(file) = &self.asteroids {
            for asteroid in &self.small_bodies {
                let other: Vec<_> = file
                    .coverage(asteroid.target)
                    .into_iter()
                    .map(|(start, end)| {
                        (
                            J2000_JD_TDB + start / SECONDS_PER_DAY,
                            J2000_JD_TDB + end / SECONDS_PER_DAY,
                        )
                    })
                    .collect();
                spans = intersection(&spans, &other);
            }
        }
        spans
    }
}

/// Resolve a physical SPICE identity to its selected DE-header GM (AU³/day²).
///
/// JPL SB441 IOM392R-21-005 Table 1 gives the 343 numbered asteroids and 30
/// KBOs. Their published GMs uniquely identify MA8001..MA8030 in the supported
/// DE440/441 headers. Only these identity aliases are retained here; every
/// value is read from the selected header, with no substituted mass.
/// Source: <https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/SB441_IOM392R-21-005_perturbers.pdf>.
fn asteroid_mass(planets: &ClassicDe, target: i32) -> Result<f64, DeError> {
    let label = match target {
        2_136_199 => "MA8001", // Eris
        2_136_108 => "MA8002", // Haumea
        2_090_377 => "MA8003", // Sedna
        2_136_472 => "MA8004", // Makemake
        2_050_000 => "MA8005", // Quaoar
        2_084_522 => "MA8006", // 2002 TC302
        2_090_482 => "MA8007", // Orcus
        2_020_000 => "MA8008", // Varuna
        2_055_637 => "MA8009", // 2002 UX25
        2_028_978 => "MA8010", // Ixion
        2_307_261 => "MA8011", // 2002 MS4
        2_174_567 => "MA8012", // Varda
        3_361_580 => "MA8013", // 2006 QH181
        3_308_265 => "MA8014", // 2004 XR190
        2_055_565 => "MA8015", // 2002 AW197
        2_145_452 => "MA8016", // 2005 RN43
        2_090_568 => "MA8017", // 2004 GV9
        2_208_996 => "MA8018", // 2003 AZ84
        2_225_088 => "MA8019", // Gonggong
        2_019_521 => "MA8020", // Chaos
        2_120_347 => "MA8021", // Salacia
        2_278_361 => "MA8022", // 2007 JJ43
        3_525_142 => "MA8023", // 2010 KZ39
        2_230_965 => "MA8024", // 2004 XA192
        2_042_301 => "MA8025", // 2001 UR163
        2_455_502 => "MA8026", // 2003 UZ413
        3_545_742 => "MA8027", // 2010 RF43
        2_523_639 => "MA8028", // 2010 RE64
        2_528_381 => "MA8029", // 2008 ST291
        3_515_022 => "MA8030", // 2010 FX86
        _ => {
            let number = target
                .checked_sub(2_000_000)
                .filter(|number| *number > 0 && !(8001..=8030).contains(number))
                .ok_or(DeError::UnsupportedTarget(target))?;
            // These reserved labels denote the KBOs above, never the actual
            // numbered asteroids 8001..8030. Refuse that physical-ID collision.
            return planets.constant(&format!("MA{number:04}"));
        }
    };
    planets.constant(label)
}

fn major(target: i32, mass: f64, value: Kinematics) -> Source {
    Source {
        point: PointMass {
            body: target,
            position: value.state.position,
            gravitational_parameter: mass,
        },
        velocity: value.state.velocity,
        post_newtonian_acceleration: Some(value.acceleration),
    }
}
fn intersection(first: &[(f64, f64)], second: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut result = Vec::new();
    let (mut i, mut j) = (0, 0);
    while i < first.len() && j < second.len() {
        let start = first[i].0.max(second[j].0);
        let end = first[i].1.min(second[j].1);
        if start <= end {
            result.push((start, end));
        }
        if first[i].1 < second[j].1 {
            i += 1;
        } else {
            j += 1;
        }
    }
    result
}

#[cfg(test)]
#[path = "ephemeris_tests.rs"]
mod tests;
