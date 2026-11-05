// SPDX-License-Identifier: AGPL-3.0
//! Safe native adapter beneath the existing public N-body API.
//!
//! This is not a C ABI or a reconstruction of third-party simulation objects.
//! It owns only the two problems actually exposed: stationary Sun/test target,
//! or one raw-SSB test target in prescribed ephemerides. Compatibility ordinals
//! and method names are translated once; all computation uses original Rust.
//! The candidate explicitly selects ICRF-terrestrial/IAU-solar harmonic axes.
//! This is a model convention, not a claim of complete upstream equivalence;
//! all numerical/resource qualification remains explicit.

use super::super::{CartesianState, GAUSSIAN_K_SQUARED, Orbit, ReboundIntegrator};
use super::catalogue::HarmonicPoles;
use super::classic_de::DeError;
use super::conic::Elements;
use super::ephemeris::{Ephemeris as Dataset, EphemerisError};
use super::poles::PoleModel;
use super::radau::IntegrationError;
use super::solar::{Method, SolarTrajectory};
use super::state::{NumericalError, State};
use super::trajectory::{ForcedTrajectory, PropagationError};
use crate::errors::{Error, Result};
use std::path::Path;
use std::sync::Arc;

/// Native implementation/model identifiers, never upstream library versions.
pub(crate) fn versions() -> (&'static str, &'static str) {
    (
        "libephemeris-native-radau15/4.0.0-experimental",
        "libephemeris-native-ppn-icrf-terrestrial-iau-solar/4.0.0-experimental",
    )
}

pub(crate) struct Ephemeris {
    data: Arc<Dataset>,
}
impl Ephemeris {
    pub(crate) fn open(planets: &Path, asteroids: Option<&Path>) -> Result<Self> {
        for path in std::iter::once(planets).chain(asteroids) {
            if path.as_os_str().as_encoded_bytes().contains(&0) {
                return Err(Error::InvalidOrbitElements {
                    entry: "assist_ephemeris",
                    detail: "a path contains a NUL byte".into(),
                });
            }
        }
        let data =
            Dataset::open(planets, asteroids).map_err(|error| acquisition_error(error, planets))?;
        Ok(Self {
            data: Arc::new(data),
        })
    }
    /// Curated relative PV from the already retained bounded force image.
    /// No whole-file image, all-record coefficient array or epoch cache is added.
    pub(crate) fn curated_initial_state(
        &self,
        target: i32,
        epoch: f64,
    ) -> Result<Option<CartesianState>> {
        finite(epoch)?;
        self.data
            .heliocentric_seed(target, [epoch, 0.0])?
            .map(|state| {
                CartesianState {
                    position: super::super::equatorial_to_ecliptic(state.position_au()),
                    velocity: super::super::equatorial_to_ecliptic(state.velocity_au_per_day()),
                }
                .check()
                .map_err(output_error)
            })
            .transpose()
    }

    pub(crate) fn reference_jd(&self) -> f64 {
        2451545.0
    }
    /// Declared common windows of these exact retained force-file owners.
    pub(crate) fn coverage(&self) -> Vec<(f64, f64)> {
        self.data.coverage()
    }
    /// Selected asteroid compatibility ordinals and physical SPICE identities.
    /// Metadata iteration borrows the retained owner and reads no coefficients.
    #[cfg(feature = "nbody")]
    pub(crate) fn asteroid_targets(&self) -> impl Iterator<Item = (i32, i32)> + '_ {
        (11..self.data.particle_count()).map(|index| {
            (
                index as i32,
                self.data
                    .physical_id(index)
                    .expect("selected asteroid identity"),
            )
        })
    }
    fn particle(&self, body: i32, time: f64) -> Result<super::forces::Source> {
        finite(time)?;
        let index =
            usize::try_from(body).map_err(|_| native_error(3, "negative force body selector"))?;
        self.data
            .particle(index, [self.reference_jd(), time])
            .map_err(|error| native_error(3, error.to_string()))
    }
    pub(crate) fn state(&self, body: i32, time: f64) -> Result<CartesianState> {
        let source = self.particle(body, time)?;
        external(State {
            position: source.point.position,
            velocity: source.velocity,
        })
    }
    #[cfg(test)]
    pub(crate) fn position(&self, body: i32, time: f64) -> Result<[f64; 3]> {
        Ok(self.particle(body, time)?.point.position)
    }
    /// Analytic selected-DE Sun acceleration, SSB ICRF AU/day².
    pub(crate) fn sun_acceleration(&self, time: f64) -> Result<[f64; 3]> {
        self.particle(0, time)?
            .post_newtonian_acceleration
            .ok_or_else(|| native_error(3, "selected Sun has no kinematic acceleration"))
    }

    #[cfg(test)]
    pub(crate) fn gravitational_parameter(&self, body: i32, time: f64) -> Result<f64> {
        Ok(self.particle(body, time)?.point.gravitational_parameter)
    }
}

struct Prescribed<'a> {
    ephemeris: &'a Ephemeris,
    initial: Option<CartesianState>,
    start: f64,
    parameters: Option<[f64; 3]>,
    exclude: Option<i32>,
    trajectory: Option<ForcedTrajectory>,
}
impl Prescribed<'_> {
    /// Initialize the retained trajectory once, including for a zero-duration
    /// acceleration query. Merely asking for PV still performs no force work.
    fn trajectory_mut(&mut self) -> Result<&mut ForcedTrajectory> {
        if self.trajectory.is_none() {
            let initial = self
                .initial
                .ok_or_else(|| native_error(5, "prescribed target has not been initialized"))?
                .check()?;
            self.trajectory = Some(
                ForcedTrajectory::new(
                    Arc::clone(&self.ephemeris.data),
                    internal(initial),
                    [self.ephemeris.reference_jd(), self.start],
                    None,
                    self.parameters,
                    self.exclude,
                    candidate_poles,
                )
                .map_err(propagation_error)?,
            );
        }
        Ok(self
            .trajectory
            .as_mut()
            .expect("just initialized owned trajectory"))
    }
}

enum Route<'a> {
    Solar(SolarTrajectory),
    Prescribed(Prescribed<'a>),
}

/// Internal adapter; public callers never receive a general particle container.
pub(crate) struct Simulation<'a> {
    route: Route<'a>,
}

impl Simulation<'_> {
    pub(crate) fn two_body(
        orbit: Orbit,
        integrator: ReboundIntegrator,
        step: Option<f64>,
    ) -> Result<Self> {
        let elements = Elements {
            semimajor_axis: orbit.semi_major_axis,
            eccentricity: orbit.eccentricity,
            inclination: orbit.inclination,
            ascending_node: orbit.ascending_node,
            periapsis_argument: orbit.perihelion_argument,
            mean_anomaly: orbit.mean_anomaly,
        };
        let method = match integrator {
            ReboundIntegrator::Ias15 => Method::Adaptive,
            ReboundIntegrator::Whfast => Method::WisdomHolman,
            ReboundIntegrator::Leapfrog => Method::Leapfrog,
            ReboundIntegrator::Mercurius => Method::Mercurius,
            ReboundIntegrator::Trace => Method::Trace,
        };
        let trajectory = SolarTrajectory::new(elements, GAUSSIAN_K_SQUARED, method, step)
            .map_err(initial_orbit_error)?;
        Ok(Self {
            route: Route::Solar(trajectory),
        })
    }

    /// Stage the one raw target used by the existing safe public adapter.
    /// Multiple-target/general simulation mutation is deliberately not exposed.
    pub(crate) fn add(&mut self, state: CartesianState) {
        let Route::Prescribed(setup) = &mut self.route else {
            panic!("only a prescribed one-target adapter stages a target")
        };
        assert!(
            setup.initial.is_none() && setup.trajectory.is_none(),
            "prescribed target is initialized once"
        );
        setup.initial = Some(state);
    }
    pub(crate) fn set_time(&mut self, time: f64) {
        let Route::Prescribed(setup) = &mut self.route else {
            panic!("only a staged prescribed adapter changes its origin")
        };
        assert!(
            setup.trajectory.is_none(),
            "an active trajectory cannot relabel its origin"
        );
        setup.start = time;
    }
    pub(crate) fn state(&self, index: usize) -> Result<CartesianState> {
        match &self.route {
            Route::Solar(trajectory) => match index {
                0 => external(State {
                    position: [0.0; 3],
                    velocity: [0.0; 3],
                }),
                1 => external(trajectory.state()),
                _ => Err(native_error(5, "Solar target index is unavailable")),
            },
            Route::Prescribed(setup) => {
                if index != 0 {
                    return Err(native_error(5, "prescribed target index is unavailable"));
                }
                if let Some(trajectory) = &setup.trajectory {
                    external(trajectory.state())
                } else {
                    setup
                        .initial
                        .ok_or_else(|| {
                            native_error(5, "prescribed target has not been initialized")
                        })?
                        .check()
                }
            }
        }
    }
    pub(crate) fn integrate(&mut self, time: f64) -> Result<Vec<String>> {
        finite(time)?;
        match &mut self.route {
            Route::Solar(trajectory) => {
                if time != trajectory.time() {
                    trajectory
                        .advance(time)
                        .map_err(propagation_numerical_error)?;
                }
            }
            Route::Prescribed(setup) => {
                finite(setup.start)?;
                let target = time - setup.start;
                finite(target)?;
                setup
                    .trajectory_mut()?
                    .advance(target)
                    .map_err(propagation_error)?;
            }
        }
        // No fabricated upstream diagnostics: native arithmetic/domain failures
        // are typed errors; the applicable exact two-body maps need no warning
        // about absent interactions. Public warning qualification remains open.
        Ok(Vec::new())
    }
}

impl Simulation<'_> {
    /// Complete force acceleration of the staged/committed prescribed target,
    /// in SSB ICRF AU/day². This never advances or substitutes its trajectory.
    pub(crate) fn acceleration(&mut self) -> Result<[f64; 3]> {
        let Route::Prescribed(setup) = &mut self.route else {
            return Err(native_error(5, "acceleration requires a prescribed target"));
        };
        setup
            .trajectory_mut()?
            .acceleration()
            .map_err(propagation_error)
    }
}

impl<'a> Simulation<'a> {
    pub(crate) fn assist(ephemeris: &'a Ephemeris) -> Result<Self> {
        Ok(Self {
            route: Route::Prescribed(Prescribed {
                ephemeris,
                initial: None,
                start: 0.0,
                parameters: None,
                exclude: None,
                trajectory: None,
            }),
        })
    }
    pub(crate) fn exclude_asteroid(&mut self, body: i32) {
        let Route::Prescribed(setup) = &mut self.route else {
            return;
        };
        assert!(
            setup.trajectory.is_none(),
            "self identity is selected before integration"
        );
        if body >= 11 {
            setup.exclude = setup.ephemeris.data.physical_id(body as usize);
        }
    }
    pub(crate) fn non_gravitational_parameters(&mut self, parameters: [f64; 3]) {
        let Route::Prescribed(setup) = &mut self.route else {
            return;
        };
        assert!(
            setup.trajectory.is_none(),
            "nongravitational parameters precede integration"
        );
        setup.parameters = Some(parameters);
    }
    #[cfg(test)]
    pub(crate) fn asteroid_acceleration(&mut self) -> Result<[f64; 3]> {
        let state = self.state(0)?;
        let Route::Prescribed(setup) = &self.route else {
            return Err(native_error(5, "not a prescribed target"));
        };
        let elapsed = setup
            .trajectory
            .as_ref()
            .map_or(0.0, ForcedTrajectory::time);
        let mut frame = Vec::new();
        setup
            .ephemeris
            .data
            .frame(
                [setup.ephemeris.reference_jd(), setup.start + elapsed],
                &mut frame,
            )
            .map_err(|error| native_error(3, error.to_string()))?;
        super::newtonian::acceleration_from(
            state.position,
            frame
                .iter()
                .filter(|source| source.post_newtonian_acceleration.is_none())
                .map(|source| source.point),
            setup.exclude,
        )
        .map_err(propagation_numerical_error)
    }
}

fn candidate_poles(epoch: [f64; 2]) -> Result<HarmonicPoles> {
    // One explicit prescribed-force convention for every body/date. The
    // separate dated-pole model is not silently inserted into this API.
    // Numerical/resource acceptance still requires the full comparisons.
    PoleModel::IcrfTerrestrialIauSolar.at(epoch)
}
fn internal(state: CartesianState) -> State {
    State {
        position: state.position,
        velocity: state.velocity,
    }
}
fn external(state: State) -> Result<CartesianState> {
    CartesianState {
        position: state.position,
        velocity: state.velocity,
    }
    .check()
    .map_err(output_error)
}
fn finite(value: f64) -> Result<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(Error::NonFiniteArgument { value })
    }
}
fn native_error(stage: i32, diagnostic: impl Into<String>) -> Error {
    Error::NBodyNative {
        stage,
        code: 1,
        diagnostic: diagnostic.into(),
    }
}

/// A consumed native state that cannot form finite public output is a failed
/// calculation, not a new invalid caller argument eligible for another model.
fn output_error(error: Error) -> Error {
    native_error(5, error.to_string())
}

/// Preserve an actual missing-file condition only while selecting the dataset.
/// Header/record corruption, unsupported data and other I/O failures retain the
/// native initialization failure. The same missing-file error during a consumed
/// force evaluation goes through `propagation_error` and remains terminal.
fn acquisition_error(error: EphemerisError, planets: &Path) -> Error {
    match &error {
        EphemerisError::De(DeError::Io(io)) if io.kind() == std::io::ErrorKind::NotFound => {
            Error::EphemerisLoaderIo {
                path: planets.display().to_string(),
                kind: io.kind(),
                errno: io.raw_os_error(),
                message: error.to_string(),
            }
        }
        EphemerisError::Spk(Error::SpkFileNotFound { path }) => Error::EphemerisLoaderIo {
            path: path.clone(),
            kind: std::io::ErrorKind::NotFound,
            errno: None,
            message: error.to_string(),
        },
        _ => native_error(2, error.to_string()),
    }
}

/// Only initial orbit construction can classify a numerical domain as input.
/// Once a selected force or solver consumes the state, the same domain error
/// is a terminal calculation failure rather than permission to change models.
fn initial_orbit_error(error: NumericalError) -> Error {
    match error {
        NumericalError::InvalidOrbit(detail) => Error::InvalidOrbitElements {
            entry: "nbody",
            detail: detail.into(),
        },
        NumericalError::InvalidMass { body } => Error::InvalidOrbitElements {
            entry: "nbody",
            detail: format!("invalid GM for body {body}"),
        },
        _ => native_error(1, format!("native numerical failure: {error:?}")),
    }
}

/// Preserve the native propagation stage, status and numerical diagnostic.
/// Invalid mass/orbit domains here describe the consumed force/state, not the
/// caller's unvalidated orbital elements. They must remain terminal too.
fn propagation_numerical_error(error: NumericalError) -> Error {
    native_error(4, format!("native numerical failure: {error:?}"))
}

fn propagation_error(error: PropagationError) -> Error {
    match error {
        IntegrationError::Numerical(error) => propagation_numerical_error(error),
        IntegrationError::Force(error) => native_error(4, error.to_string()),
    }
}

#[cfg(test)]
#[path = "bridge_tests.rs"]
mod tests;
