// SPDX-License-Identifier: AGPL-3.0
//! Rust orbit propagation with explicit selected JPL force data.
//!
//! The engine integrates Solar and prescribed-ephemeris forces independently
//! in Rust, including relativistic, harmonic and optional nongravitational
//! accelerations. Historical REBOUND/ASSIST API names remain compatibility
//! names. Public states are heliocentric ecliptic J2000, AU and AU/day at scalar
//! TT Julian days; prescribed-force propagation uses TT as its TDB approximation.
//!
//! The optional `nbody` feature enables this engine; `native-nbody` is a
//! compatible alias. Without the feature availability is false and integration
//! returns a typed unavailable error. Inputs and derived states must be finite.
//! A zero-duration propagation returns the initial state without advancing the
//! solver. Nonfinite timesteps are rejected, as is zero for a nonzero interval.

use crate::errors::{Error, Result, Warning};
use crate::minor_bodies::{OrbitalElements, apply_secular_perturbations, calc_minor_body_position};
use std::path::PathBuf;
use std::sync::{Arc, Mutex, OnceLock};

#[cfg(feature = "nbody")]
use native::bridge as backend;

#[cfg(any(test, feature = "nbody"))]
#[path = "nbody/native/mod.rs"]
mod native;

#[path = "nbody_data.rs"]
mod data;
#[cfg(feature = "nbody")]
#[path = "nbody/source.rs"]
mod source;
pub use data::{
    ASSIST_ASTEROIDS_URL, ASSIST_EXPECTED_SIZES, ASSIST_MIN_SIZE_RATIO, ASSIST_PLANETS_DE441_URL,
    ASSIST_PLANETS_URL, verify_assist_file,
};
#[cfg(feature = "nbody")]
pub use source::NbodySpkSource;

/// Solar gravitational parameter in AU³/day², the square of the Gaussian constant.
pub const GAUSSIAN_K_SQUARED: f64 = 0.00029591220828559;
/// Planetary ephemeris names searched in precision-tier order.
pub const ASSIST_PLANET_FILENAMES: [&str; 2] = ["linux_p1550p2650.440", "linux_m13000p17000.441"];
/// Massive-asteroid ephemeris name.
pub const ASSIST_ASTEROID_FILENAME: &str = "sb441-n16.bsp";
/// Fixed mean J2000 obliquity used by the orbital ecliptic/ICRF convention, degrees.
pub const ECLIPTIC_J2000_OBLIQUITY_DEG: f64 = 23.4392911;

/// A compatibility-named method for the exposed Sun/test-target problem.
#[repr(i32)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ReboundIntegrator {
    /// Adaptive fifteenth-order Gauss–Radau integration.
    #[default]
    Ias15 = 0,
    /// Fixed-step Wisdom–Holman integration.
    Whfast = 1,
    /// Fixed-step leapfrog integration.
    Leapfrog = 4,
    /// Hybrid close-encounter integration.
    Mercurius = 9,
    /// Time-reversible hybrid close-encounter integration.
    Trace = 25,
}

impl ReboundIntegrator {
    /// Lowercase public name of this method.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Ias15 => "ias15",
            Self::Whfast => "whfast",
            Self::Leapfrog => "leapfrog",
            Self::Mercurius => "mercurius",
            Self::Trace => "trace",
        }
    }

    /// Parse a method name, returning an input error for an unknown name.
    pub fn parse(name: &str) -> Result<Self> {
        match name {
            "ias15" => Ok(Self::Ias15),
            "whfast" => Ok(Self::Whfast),
            "leapfrog" => Ok(Self::Leapfrog),
            "mercurius" => Ok(Self::Mercurius),
            "trace" => Ok(Self::Trace),
            _ => Err(Error::UnknownIntegrator { name: name.into() }),
        }
    }

    /// Whether the method uses the default orbital-period-based fixed step.
    pub fn is_fixed_step(self) -> bool {
        matches!(self, Self::Whfast | Self::Leapfrog)
    }
}

/// Heliocentric ecliptic-J2000 Cartesian position and velocity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CartesianState {
    /// Position components in AU.
    pub position: [f64; 3],
    /// Velocity components in AU/day.
    pub velocity: [f64; 3],
}

impl CartesianState {
    fn check(self) -> Result<Self> {
        for value in self.position.into_iter().chain(self.velocity) {
            if !value.is_finite() {
                return Err(Error::NonFiniteArgument { value });
            }
        }
        Ok(self)
    }
}

/// Osculating elements supplied to the native orbit constructor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Orbit {
    /// Semi-major axis in AU; negative for an allowed hyperbolic orbit.
    pub semi_major_axis: f64,
    /// Dimensionless eccentricity.
    pub eccentricity: f64,
    /// Inclination in radians.
    pub inclination: f64,
    /// Argument of perihelion in radians.
    pub perihelion_argument: f64,
    /// Longitude of ascending node in radians.
    pub ascending_node: f64,
    /// Mean anomaly at initialization in radians.
    pub mean_anomaly: f64,
}

/// Evaluate secularly perturbed elements at TT, converting only angular units.
///
/// Uses the native minor-body secular model. The input epoch is a finite Julian
/// day; the resulting orbit retains `a` in AU and angles in radians. Invalid
/// nonfinite elements are refused before native integration.
pub fn elements_to_rebound_orbit(elements: &OrbitalElements, jd_tt: f64) -> Result<Orbit> {
    finite(jd_tt)?;
    for value in [
        elements.epoch(),
        elements.a(),
        elements.e(),
        elements.i(),
        elements.omega(),
        elements.node(),
        elements.m0(),
        elements.n(),
    ] {
        finite(value)?;
    }
    let (perihelion, node, anomaly, _, eccentricity, inclination) =
        apply_secular_perturbations(elements, jd_tt, true);
    let orbit = Orbit {
        semi_major_axis: elements.a(),
        eccentricity,
        inclination: inclination.to_radians(),
        perihelion_argument: perihelion.to_radians(),
        ascending_node: node.to_radians(),
        mean_anomaly: anomaly.to_radians(),
    };
    for value in [
        orbit.semi_major_axis,
        orbit.eccentricity,
        orbit.inclination,
        orbit.perihelion_argument,
        orbit.ascending_node,
        orbit.mean_anomaly,
    ] {
        finite(value)?;
    }
    Ok(orbit)
}

/// A propagated state at a scalar TT Julian day.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PropagationResult {
    /// Heliocentric ecliptic-J2000 x, AU.
    pub x: f64,
    /// Heliocentric ecliptic-J2000 y, AU.
    pub y: f64,
    /// Heliocentric ecliptic-J2000 z, AU.
    pub z: f64,
    /// Heliocentric velocity x, AU/day.
    pub vx: f64,
    /// Heliocentric velocity y, AU/day.
    pub vy: f64,
    /// Heliocentric velocity z, AU/day.
    pub vz: f64,
    /// Output epoch, TT Julian day.
    pub jd_tt: f64,
}

impl PropagationResult {
    /// Construct a finite result; nonfinite state components or time are errors.
    pub fn from_state(state: CartesianState, jd_tt: f64) -> Result<Self> {
        let state = state.check()?;
        finite(jd_tt)?;
        let [x, y, z] = state.position;
        let [vx, vy, vz] = state.velocity;
        Ok(Self {
            x,
            y,
            z,
            vx,
            vy,
            vz,
            jd_tt,
        })
    }

    /// First-order state at a nearby TT epoch, with constant Cartesian velocity.
    ///
    /// This local approximation supplies apparent-place light-time and
    /// derivative samples. It does not replace orbit integration over a finite
    /// requested interval. Inputs and derived components must remain finite.
    pub fn shifted(self, jd_tt: f64) -> Result<Self> {
        finite(jd_tt)?;
        let interval = jd_tt - self.jd_tt;
        finite(interval)?;
        Self::from_state(
            CartesianState {
                position: [
                    self.x + self.vx * interval,
                    self.y + self.vy * interval,
                    self.z + self.vz * interval,
                ],
                velocity: [self.vx, self.vy, self.vz],
            },
            jd_tt,
        )
    }

    /// Euclidean heliocentric radius, AU; scaled hypot avoids intermediate overflow.
    pub fn distance(self) -> f64 {
        self.x.hypot(self.y).hypot(self.z)
    }
    /// Ecliptic longitude in degrees on [0,360).
    pub fn ecliptic_lon(self) -> f64 {
        self.y.atan2(self.x).to_degrees().rem_euclid(360.0)
    }
    /// Ecliptic latitude in degrees; the zero vector has latitude zero.
    pub fn ecliptic_lat(self) -> f64 {
        let radius = self.distance();
        if radius == 0.0 {
            0.0
        } else {
            (self.z / radius).clamp(-1.0, 1.0).asin().to_degrees()
        }
    }
    /// Position and velocity components in AU and AU/day.
    pub fn to_tuple(self) -> (f64, f64, f64, f64, f64, f64) {
        (self.x, self.y, self.z, self.vx, self.vy, self.vz)
    }
}

/// Owned values plus nonfatal native diagnostics, with no interpreter state.
#[derive(Clone, Debug)]
pub struct PropagationOutcome<T> {
    /// Successfully computed result.
    pub value: T,
    /// Actual upstream diagnostic messages.
    pub messages: Vec<String>,
    /// Native physical warnings such as self-interaction avoidance.
    pub warnings: Vec<Warning>,
}

/// Explicit ASSIST paths and the directory searched for missing selections.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AssistEphemConfig {
    /// Planet ephemeris path, or discover under the configured directories.
    pub planets_file: Option<PathBuf>,
    /// Optional massive-asteroid ephemeris path.
    pub asteroids_file: Option<PathBuf>,
    /// Highest-priority discovery directory after explicit file paths.
    pub data_dir: Option<PathBuf>,
}

impl AssistEphemConfig {
    /// Search explicit files, `data_dir`, ASSIST_DIR, the user data directory,
    /// then./data. This reads directory metadata only and never downloads.
    /// Explicit paths remain explicit even when missing; consumers report them.
    pub fn discover(&self) -> Self {
        let mut directories = Vec::new();
        if let Some(directory) = &self.data_dir {
            directories.push(directory.clone());
        }
        if let Some(directory) = std::env::var_os("ASSIST_DIR") {
            directories.push(directory.into());
        }
        if let Some(home) = crate::platform::user_home_directory() {
            directories.push(home.join(".libephemeris/assist"));
        }
        directories.push(PathBuf::from("data"));
        let find = |names: &[&str]| {
            directories
                .iter()
                .flat_map(|directory| names.iter().map(move |name| directory.join(name)))
                .find(|path| path.is_file())
        };
        Self {
            planets_file: self
                .planets_file
                .clone()
                .or_else(|| find(&ASSIST_PLANET_FILENAMES)),
            asteroids_file: self
                .asteroids_file
                .clone()
                .or_else(|| find(&[ASSIST_ASTEROID_FILENAME])),
            data_dir: self.data_dir.clone(),
        }
    }

    /// Whether the resolved planet path exists and an explicit asteroid path,
    /// when selected, also exists. This does not certify ephemeris coverage.
    pub fn data_available(&self) -> bool {
        let resolved = self.discover();
        resolved.planets_file.is_some_and(|path| path.is_file())
            && resolved.asteroids_file.is_none_or(|path| path.is_file())
    }
}

/// Non-gravitational radial, transverse and normal Marsden accelerations.
#[derive(Clone, Copy, Debug, Default)]
pub struct NonGravitationalParameters {
    /// Radial coefficient, AU/day² at the reference distance.
    pub radial: f64,
    /// Transverse coefficient, AU/day² at the reference distance.
    pub transverse: f64,
    /// Normal coefficient, AU/day² at the reference distance.
    pub normal: f64,
}

/// Whether independent Rust Solar integration is enabled by `nbody`.
pub fn check_rebound_available() -> bool {
    cfg!(feature = "nbody")
}
/// Whether independent Rust selected-ephemeris integration is enabled by `nbody`.
pub fn check_assist_available() -> bool {
    cfg!(feature = "nbody")
}
/// Compatibility identity of the Rust Solar engine, or `None` when disabled.
pub fn get_rebound_version() -> Option<&'static str> {
    #[cfg(feature = "nbody")]
    {
        Some(backend::versions().0)
    }
    #[cfg(not(feature = "nbody"))]
    {
        None
    }
}
/// Compatibility identity of the Rust selected-ephemeris engine, or `None` when disabled.
pub fn get_assist_version() -> Option<&'static str> {
    #[cfg(feature = "nbody")]
    {
        Some(backend::versions().1)
    }
    #[cfg(not(feature = "nbody"))]
    {
        None
    }
}
/// Whether ASSIST's planetary data are discoverable. No download or mapping occurs.
pub fn check_assist_data_available() -> bool {
    check_assist_available() && AssistEphemConfig::default().data_available()
}
/// Refresh live discovery and the model idle default propagation resources.
/// Active calculations retain their owned mappings until they finish.
pub fn reset_assist_data_cache() {
    release_default_assist_session();
}

fn finite(value: f64) -> Result<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(Error::NonFiniteArgument { value })
    }
}

fn invalid(detail: impl Into<String>) -> Error {
    Error::InvalidOrbitElements {
        entry: "nbody",
        detail: detail.into(),
    }
}

/// Approximate orbital period divided by twenty, in days, for fixed-step methods.
/// A finite positive semi-major axis is required; no replacement step is fitted.
pub fn default_timestep_days(semi_major_axis_au: f64) -> Result<f64> {
    finite(semi_major_axis_au)?;
    if semi_major_axis_au <= 0.0 {
        return Err(invalid(
            "a fixed-step default needs a positive semi-major axis",
        ));
    }
    let step = 365.25 * semi_major_axis_au.powf(1.5) / 20.0;
    finite(step)?;
    Ok(step)
}

fn timestep(
    elements: &OrbitalElements,
    integrator: ReboundIntegrator,
    explicit: Option<f64>,
) -> Result<Option<f64>> {
    if let Some(value) = explicit {
        finite(value)?;
        Ok(Some(value))
    } else if integrator.is_fixed_step() {
        default_timestep_days(elements.a()).map(Some)
    } else {
        Ok(None)
    }
}

/// Propagate a massless target about a unit-mass Sun with the selected backend.
///
/// Epochs are finite TT Julian days, `dt` is an optional finite nonzero step in
/// days (zero is allowed only for a zero interval). Fixed-step methods without `dt` use [`default_timestep_days`]. Native
/// orbit-construction and integration failures propagate; no analytical state
/// replaces an unsuccessful integration.
pub fn propagate_orbit_rebound(
    elements: &OrbitalElements,
    jd_start: f64,
    jd_end: f64,
    integrator: ReboundIntegrator,
    dt: Option<f64>,
) -> Result<PropagationResult> {
    Ok(propagate_orbit_rebound_with_warnings(elements, jd_start, jd_end, integrator, dt)?.value)
}

/// REBOUND propagation retaining actual native diagnostics with the result.
pub fn propagate_orbit_rebound_with_warnings(
    elements: &OrbitalElements,
    jd_start: f64,
    jd_end: f64,
    integrator: ReboundIntegrator,
    dt: Option<f64>,
) -> Result<PropagationOutcome<PropagationResult>> {
    finite(jd_start)?;
    finite(jd_end)?;
    let step = timestep(elements, integrator, dt)?;
    if jd_start != jd_end && step == Some(0.0) {
        return Err(invalid(
            "a nonzero integration interval needs a nonzero timestep",
        ));
    }
    let orbit = elements_to_rebound_orbit(elements, jd_start)?;
    #[cfg(feature = "nbody")]
    {
        let mut simulation = backend::Simulation::two_body(orbit, integrator, step)?;
        let interval = jd_end - jd_start;
        finite(interval)?;
        let messages = simulation.integrate(interval)?;
        Ok(PropagationOutcome {
            value: PropagationResult::from_state(simulation.state(1)?, jd_end)?,
            messages,
            warnings: Vec::new(),
        })
    }
    #[cfg(not(feature = "nbody"))]
    {
        let _ = (step, orbit);
        Err(Error::NBodyBackendUnavailable { backend: "rebound" })
    }
}

#[cfg(any(feature = "nbody", test))]
fn ecliptic_to_equatorial(vector: [f64; 3]) -> [f64; 3] {
    let (sine, cosine) = ECLIPTIC_J2000_OBLIQUITY_DEG.to_radians().sin_cos();
    [
        vector[0],
        vector[1] * cosine - vector[2] * sine,
        vector[1] * sine + vector[2] * cosine,
    ]
}

#[cfg(any(feature = "nbody", test))]
fn equatorial_to_ecliptic(vector: [f64; 3]) -> [f64; 3] {
    let (sine, cosine) = ECLIPTIC_J2000_OBLIQUITY_DEG.to_radians().sin_cos();
    [
        vector[0],
        vector[1] * cosine + vector[2] * sine,
        -vector[1] * sine + vector[2] * cosine,
    ]
}

#[cfg(feature = "nbody")]
fn shifted(state: CartesianState, center: CartesianState, add: bool) -> CartesianState {
    let direction = if add { 1.0 } else { -1.0 };
    CartesianState {
        position: std::array::from_fn(|index| {
            state.position[index] + direction * center.position[index]
        }),
        velocity: std::array::from_fn(|index| {
            state.velocity[index] + direction * center.velocity[index]
        }),
    }
}

#[cfg(feature = "nbody")]
fn rotated(state: CartesianState, to_equatorial: bool) -> CartesianState {
    let rotate = if to_equatorial {
        ecliptic_to_equatorial
    } else {
        equatorial_to_ecliptic
    };
    CartesianState {
        position: rotate(state.position),
        velocity: rotate(state.velocity),
    }
}

#[cfg(feature = "nbody")]
fn open_ephemeris(config: Option<&AssistEphemConfig>) -> Result<backend::Ephemeris> {
    let config = config.cloned().unwrap_or_default().discover();
    let planets = config
        .planets_file
        .as_deref()
        .ok_or_else(|| Error::EphemerisLoaderIo {
            path: "ASSIST planets".into(),
            kind: std::io::ErrorKind::NotFound,
            errno: None,
            message: "planet ephemeris was not found".into(),
        })?;
    for path in std::iter::once(planets).chain(config.asteroids_file.as_deref()) {
        std::fs::metadata(path).map_err(|error| Error::EphemerisLoaderIo {
            path: path.display().to_string(),
            kind: error.kind(),
            errno: error.raw_os_error(),
            message: error.to_string(),
        })?;
    }
    backend::Ephemeris::open(planets, config.asteroids_file.as_deref())
}

/// Reusable JPL force-file owners for the ASSIST-compatible interface.
///
/// A worker can retain this owner between propagations. It can move to another
/// worker, but simultaneous access requires an application mutex. Attached
/// simulations borrow it only during their synchronous calculations. Selected
/// files must remain immutable until this owner is dropped.
pub struct AssistEphemeris {
    #[cfg(feature = "nbody")]
    inner: backend::Ephemeris,
    _private: (),
}

impl AssistEphemeris {
    /// Resolve and open explicit or discovered JPL files; missing files are I/O
    /// errors, malformed ephemerides are native errors. No download occurs.
    pub fn open(config: Option<&AssistEphemConfig>) -> Result<Self> {
        #[cfg(feature = "nbody")]
        {
            Ok(Self {
                inner: open_ephemeris(config)?,
                _private: (),
            })
        }
        #[cfg(not(feature = "nbody"))]
        {
            let _ = config;
            Err(Error::NBodyBackendUnavailable { backend: "assist" })
        }
    }

    /// Propagate with these retained ephemerides, preserving the diagnostics of
    /// the native integration. Units and limits match [`propagate_orbit_assist`].
    pub fn propagate(
        &self,
        elements: &OrbitalElements,
        jd_start: f64,
        jd_end: f64,
        non_gravitational: Option<NonGravitationalParameters>,
    ) -> Result<PropagationOutcome<PropagationResult>> {
        finite(jd_start)?;
        finite(jd_end)?;
        #[cfg(feature = "nbody")]
        {
            let initial = assist_element_state(elements, jd_start)?;
            self.propagate_state_for_target(
                initial,
                jd_start,
                jd_end,
                non_gravitational,
                elements.physical_spice_id(),
            )
        }
        #[cfg(not(feature = "nbody"))]
        {
            let _ = (elements.physical_spice_id(), non_gravitational);
            Err(Error::NBodyBackendUnavailable { backend: "assist" })
        }
    }

    /// Integrate an explicit heliocentric ecliptic-J2000 initial state.
    ///
    /// Positions are AU, velocities AU/day, and epochs scalar TT JD under the
    /// documented TT-as-TDB approximation. The supplied state is used exactly;
    /// no catalogue elements or alternative source replaces it. Finite inputs,
    /// ephemeris coverage and native integration errors follow [`Self::propagate`].
    pub fn propagate_state(
        &self,
        initial: CartesianState,
        jd_start: f64,
        jd_end: f64,
        non_gravitational: Option<NonGravitationalParameters>,
    ) -> Result<PropagationOutcome<PropagationResult>> {
        self.propagate_state_for_target(initial, jd_start, jd_end, non_gravitational, None)
    }

    /// Only an actual body/source selector supplies the physical self identity.
    fn propagate_state_for_target(
        &self,
        initial: CartesianState,
        jd_start: f64,
        jd_end: f64,
        non_gravitational: Option<NonGravitationalParameters>,
        self_target: Option<i32>,
    ) -> Result<PropagationOutcome<PropagationResult>> {
        let initial = initial.check()?;
        finite(jd_start)?;
        finite(jd_end)?;
        #[cfg(feature = "nbody")]
        {
            let outcome = self.propagated_state_with(
                initial,
                jd_start,
                jd_end,
                non_gravitational,
                self_target,
                |_, _, _| Ok(()),
            )?;
            Ok(PropagationOutcome {
                value: outcome.value.0,
                messages: outcome.messages,
                warnings: outcome.warnings,
            })
        }
        #[cfg(not(feature = "nbody"))]
        {
            let _ = (initial, non_gravitational, self_target);
            Err(Error::NBodyBackendUnavailable { backend: "assist" })
        }
    }

    /// Shared integration route with an optional typed endpoint observation.
    /// Inputs have already passed the public finite/state checks. The callback
    /// borrows the same simulation before release; it never reintegrates it.
    #[cfg(feature = "nbody")]
    fn propagated_state_with<T>(
        &self,
        initial: CartesianState,
        jd_start: f64,
        jd_end: f64,
        non_gravitational: Option<NonGravitationalParameters>,
        self_target: Option<i32>,
        finish: impl FnOnce(&mut backend::Simulation<'_>, &backend::Ephemeris, f64) -> Result<T>,
    ) -> Result<PropagationOutcome<(PropagationResult, T)>> {
        let ephemeris = &self.inner;
        let (mut simulation, initial, warnings) =
            assist_initial_state(initial, jd_start, ephemeris, non_gravitational, self_target)?;
        let time = jd_end - ephemeris.reference_jd();
        let (state, messages) = if jd_start == jd_end {
            (initial, Vec::new())
        } else {
            let messages = simulation.integrate(time)?;
            let state = rotated(
                shifted(simulation.state(0)?, ephemeris.state(0, time)?, false),
                false,
            );
            (state, messages)
        };
        let value = PropagationResult::from_state(state, jd_end)?;
        let endpoint = finish(&mut simulation, ephemeris, time)?;
        Ok(PropagationOutcome {
            value: (value, endpoint),
            messages,
            warnings,
        })
    }

    /// Integrate a sequence of finite TT Julian days in one simulation while
    /// retaining the selected file owners. The initial sample requires no step.
    pub fn trajectory(
        &self,
        elements: &OrbitalElements,
        jd_start: f64,
        times: &[f64],
    ) -> Result<PropagationOutcome<Vec<PropagationResult>>> {
        finite(jd_start)?;
        for &time in times {
            finite(time)?;
        }
        #[cfg(feature = "nbody")]
        {
            let ephemeris = &self.inner;
            let (mut simulation, initial, warnings) =
                assist_initial(elements, jd_start, ephemeris, None, None)?;
            let mut values = Vec::with_capacity(times.len());
            let mut messages = Vec::new();
            for &time in times {
                let state = if time == jd_start {
                    initial
                } else {
                    let relative = time - ephemeris.reference_jd();
                    messages.extend(simulation.integrate(relative)?);
                    rotated(
                        shifted(simulation.state(0)?, ephemeris.state(0, relative)?, false),
                        false,
                    )
                };
                values.push(PropagationResult::from_state(state, time)?);
            }
            Ok(PropagationOutcome {
                value: values,
                messages,
                warnings,
            })
        }
        #[cfg(not(feature = "nbody"))]
        {
            let _ = elements;
            Err(Error::NBodyBackendUnavailable { backend: "assist" })
        }
    }
}

/// Thread-safe reusable ASSIST configuration and lazily acquired force files.
///
/// A session owns immutable discovered paths and opens them on its first
/// calculation. Concurrent users serialize access to the retained ephemeris owner.
/// Cloning an `Arc<AssistSession>` retains its files until active calls finish;
/// selecting another configuration therefore never invalidates a running call.
/// The selected files must remain unchanged while this session exists.
pub struct AssistSession {
    config: AssistEphemConfig,
    ephemeris: Mutex<Option<AssistEphemeris>>,
    #[cfg(feature = "nbody")]
    stored_trajectories: Mutex<Option<crate::spk::PositionalChebyshev>>,
}

/// Distinguish unavailable initial acquisition from a failed physical emission.
/// Once the central orbit is acquired, source failures must not select a
/// different trajectory through the ordinary optional-backend fallback.
#[derive(Debug)]
pub(crate) enum OwnedTrajectoryError {
    Initial(Error),
    Emission {
        error: Error,
        messages: Vec<String>,
        warnings: Vec<Warning>,
    },
}

/// A retained trajectory outcome with distinct acquisition/emission failures.
pub(crate) type OwnedTrajectoryResult<T> =
    std::result::Result<PropagationOutcome<T>, OwnedTrajectoryError>;

impl OwnedTrajectoryError {
    pub(crate) fn after_acquisition(error: Error) -> Self {
        Self::Emission {
            error,
            messages: Vec::new(),
            warnings: Vec::new(),
        }
    }
}

pub(crate) type BarycentricSampler<'a> = dyn FnMut(f64) -> Result<CartesianState> + 'a;
pub(crate) type SourceSunSampler<'a> = dyn Fn(f64) -> Result<CartesianState> + 'a;

impl AssistSession {
    /// Resolve a configuration without mapping files or downloading data.
    pub fn new(config: AssistEphemConfig) -> Self {
        Self {
            config: config.discover(),
            ephemeris: Mutex::new(None),
            #[cfg(feature = "nbody")]
            stored_trajectories: Mutex::new(None),
        }
    }

    /// Resolved read-only paths used by this session.
    pub fn config(&self) -> &AssistEphemConfig {
        &self.config
    }

    fn with_ephemeris<T>(
        &self,
        operation: impl FnOnce(&AssistEphemeris) -> Result<T>,
    ) -> Result<T> {
        let mut mapped = self
            .ephemeris
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        if mapped.is_none() {
            *mapped = Some(AssistEphemeris::open(Some(&self.config))?);
        }
        operation(mapped.as_ref().expect("ephemeris initialized above"))
    }

    /// Propagate using retained ephemerides; units and errors match
    /// [`propagate_orbit_assist_with_warnings`].
    pub fn propagate(
        &self,
        elements: &OrbitalElements,
        jd_start: f64,
        jd_end: f64,
        non_gravitational: Option<NonGravitationalParameters>,
    ) -> Result<PropagationOutcome<PropagationResult>> {
        self.with_ephemeris(|ephemeris| {
            ephemeris.propagate(elements, jd_start, jd_end, non_gravitational)
        })
    }

    /// Borrow the selected stored or integrated orbit and its force Sun.
    /// Covered Type2 targets in the selected asteroid image use their source
    /// records; other targets and segment types retain catalogue integration. Sample offsets are TT days from reception; the native force convention
    /// is TT-as-TDB. Keeping relative times avoids subtracting a delay from a
    /// large Julian date. The simulation, attached ephemeris and lock remain
    /// scoped together; no sampler or native pointer escapes the operation.
    pub(crate) fn with_barycentric_trajectory<T>(
        &self,
        body: i32,
        elements: &OrbitalElements,
        jd_tt: f64,
        operation: impl FnOnce(
            &PropagationResult,
            &mut BarycentricSampler<'_>,
            &SourceSunSampler<'_>,
        ) -> Result<T>,
    ) -> OwnedTrajectoryResult<T> {
        self.with_trajectory_observation(body, elements, jd_tt, false, |central, sample, sun, _| {
            operation(central, sample, sun)
        })
    }

    /// Acquire the same source for PV-only, curved and endpoint-force requests.
    /// A stored positional trajectory has no force-model acceleration claim.
    fn with_trajectory_observation<T>(
        &self,
        body: i32,
        elements: &OrbitalElements,
        jd_tt: f64,
        observe_acceleration: bool,
        operation: impl FnOnce(
            &PropagationResult,
            &mut BarycentricSampler<'_>,
            &SourceSunSampler<'_>,
            Option<[f64; 3]>,
        ) -> Result<T>,
    ) -> OwnedTrajectoryResult<T> {
        finite(jd_tt).map_err(OwnedTrajectoryError::Initial)?;
        #[cfg(feature = "nbody")]
        {
            self.with_ephemeris(|owner| {
                let ephemeris = &owner.inner;
                // A target already stored in the selected force image has its
                // own precise curved trajectory. Reintegrating catalogue
                // elements would discard that available source solution.
                if let (Some(path), Some((_, target))) = (
                    self.config.asteroids_file.as_deref(),
                    crate::constants::get_spk_body_info_from_map(body),
                ) {
                    let mut retained = self
                        .stored_trajectories
                        .lock()
                        .unwrap_or_else(|error| error.into_inner());
                    if retained.is_none() {
                        *retained = Some(crate::spk::PositionalChebyshev::open(path)?);
                    }
                    let reader = retained.as_mut().expect("stored source initialized above");
                    if reader.contains(target) && reader.trajectory_covers(target, jd_tt, 0.0)? {
                        let (state, descriptor) = reader.trajectory_state(target, jd_tt, 0.0)?;
                        if descriptor.center != crate::spk::CENTER_SUN {
                            return Err(Error::SpkAmbiguousCenter { target });
                        }
                        return Ok(stored_trajectory(
                            ephemeris, reader, target, jd_tt, state,
                            |central, sample, sun| operation(central, sample, sun, None),
                        ));
                    }
                    // Absent targets and genuine descriptor coverage misses
                    // retain the existing numerical catalogue trajectory.
                }
                let (mut simulation, _initial, warnings) =
                    assist_initial(
                        elements,
                        elements.epoch(),
                        ephemeris,
                        None,
                        crate::constants::get_spk_body_info_from_map(body).map(|(_, target)| target),
                    )?;
                #[cfg(test)]
                owned_input_observations::note(|| serde_json::json!({
                    "kind":"catalogue_initial", "name":elements.name(),
                    "elements":[elements.epoch(),elements.a(),elements.e(),elements.i(),elements.omega(),elements.node(),elements.m0(),elements.n()],
                    "initial_orbital_pv":owned_input_observations::pv(_initial),
                    "reception_tt":jd_tt,"reference_jd":ephemeris.reference_jd(),
                    "planets_file":self.config.planets_file,"asteroids_file":self.config.asteroids_file,
                    "source_clock":"TT-as-TDB"
                }));
                let reception_time = jd_tt - ephemeris.reference_jd();
                finite(reception_time)?;
                let mut messages = if jd_tt == elements.epoch() {
                    Vec::new()
                } else {
                    simulation.integrate(reception_time)?
                };
                let central_barycentric = simulation.state(0)?.check()?;
                let central_sun = ephemeris.state(0, reception_time)?;
                #[cfg(test)]
                owned_input_observations::note(|| serde_json::json!({
                    "kind":"reception_source", "reception_tt":jd_tt,
                    "reference_jd":ephemeris.reference_jd(),"relative_time":reception_time,
                    "target_ssb_icrf":owned_input_observations::pv(central_barycentric),
                    "sun_ssb_icrf":owned_input_observations::pv(central_sun)
                }));
                let central = PropagationResult::from_state(
                    rotated(shifted(central_barycentric, central_sun, false), false),
                    jd_tt,
                )?;
                let acceleration = if observe_acceleration {
                    match endpoint_acceleration(&mut simulation, ephemeris, reception_time) {
                        Ok(acceleration) => Some(acceleration),
                        Err(error) => {
                            return Ok(Err(OwnedTrajectoryError::Emission {
                                error,
                                messages,
                                warnings,
                            }));
                        }
                    }
                } else {
                    None
                };
                let mut sample = |offset: f64| {
                    finite(offset)?;
                    if offset == 0.0 {
                        #[cfg(test)]
                        owned_input_observations::note(|| serde_json::json!({
                            "kind":"target_sample", "offset":offset,"relative_time":reception_time,
                            "cached_reception":true,"ssb_icrf":owned_input_observations::pv(central_barycentric)
                        }));
                        return Ok(central_barycentric);
                    }
                    let time = reception_time + offset;
                    finite(time)?;
                    messages.extend(simulation.integrate(time)?);
                    let state = simulation.state(0)?.check()?;
                    #[cfg(test)]
                    owned_input_observations::note(|| serde_json::json!({
                        "kind":"target_sample", "offset":offset,"relative_time":time,
                        "cached_reception":false,"ssb_icrf":owned_input_observations::pv(state)
                    }));
                    Ok(state)
                };
                let sun = |offset: f64| {
                    finite(offset)?;
                    if offset == 0.0 {
                        #[cfg(test)]
                        owned_input_observations::note(|| serde_json::json!({
                            "kind":"force_sun_sample", "offset":offset,"relative_time":reception_time,
                            "cached_reception":true,"ssb_icrf":owned_input_observations::pv(central_sun)
                        }));
                        return Ok(central_sun);
                    }
                    let time = reception_time + offset;
                    finite(time)?;
                    let state = ephemeris.state(0, time)?.check()?;
                    #[cfg(test)]
                    owned_input_observations::note(|| serde_json::json!({
                        "kind":"force_sun_sample", "offset":offset,"relative_time":time,
                        "cached_reception":false,"ssb_icrf":owned_input_observations::pv(state)
                    }));
                    Ok(state)
                };
                let result = operation(&central, &mut sample, &sun, acceleration);
                Ok(match result {
                    Ok(value) => Ok(PropagationOutcome {
                        value,
                        messages,
                        warnings,
                    }),
                    Err(error) => Err(OwnedTrajectoryError::Emission {
                        error,
                        messages,
                        warnings,
                    }),
                })
            })
            .map_err(OwnedTrajectoryError::Initial)?
        }
        #[cfg(not(feature = "nbody"))]
        {
            let _ = (body, elements, observe_acceleration, operation);
            Err(OwnedTrajectoryError::Initial(
                Error::NBodyBackendUnavailable { backend: "assist" },
            ))
        }
    }

    /// Return the selected central PV and, for an integrated orbit, the actual
    /// endpoint acceleration in heliocentric ecliptic J2000 AU/day². Stored
    /// source PV retains its original opaque-supplier derivative contract.
    #[cfg(feature = "nbody")]
    pub(crate) fn state_with_acceleration(
        &self,
        body: i32,
        elements: &OrbitalElements,
        jd_tt: f64,
    ) -> OwnedTrajectoryResult<(PropagationResult, Option<[f64; 3]>)> {
        self.with_trajectory_observation(
            body,
            elements,
            jd_tt,
            true,
            |central, _, _, acceleration| Ok((*central, acceleration)),
        )
    }

    /// Integrate a curated minor body from its catalogue epoch to scalar TT JD.
    ///
    /// When this session's selected asteroid force file contains that target,
    /// its heliocentric J2000 position and analytic velocity initialize the
    /// integration at the catalogue epoch. The initial condition and asteroid
    /// forcing then use the same stored ephemeris solution. Bodies absent from
    /// that file retain their curated elements. Nonzero intervals are integrated;
    /// zero duration returns the selected initial state. No final-epoch SPK sample
    /// bypasses integration; no download or alternate file selection occurs.
    /// Analytic seed velocities use the same bounded asteroid image as the
    /// force model. Later calls share its resources, not cached epoch results.
    /// Runtime apparent positions may instead sample a covered stored curve.
    ///
    /// Output units/frame and TT-as-TDB convention follow [`PropagationResult`].
    /// Unknown bodies, malformed source frames/centers, uncovered initial epochs
    /// and native integration failures propagate as errors.
    pub fn propagate_curated(
        &self,
        body: i32,
        jd_tt: f64,
    ) -> Result<PropagationOutcome<PropagationResult>> {
        finite(jd_tt)?;
        let elements = crate::minor_bodies::minor_body_elements(body)
            .ok_or(Error::UnknownBody { body_id: body })?;
        self.with_ephemeris(|ephemeris| {
            #[cfg(feature = "nbody")]
            {
                let initial = match self.curated_initial_state(ephemeris, body, elements.epoch())? {
                    Some(initial) => initial,
                    None => assist_element_state(&elements, elements.epoch())?,
                };
                ephemeris.propagate_state_for_target(
                    initial,
                    elements.epoch(),
                    jd_tt,
                    None,
                    crate::constants::get_spk_body_info_from_map(body).map(|(_, target)| target),
                )
            }
            #[cfg(not(feature = "nbody"))]
            {
                let _ = (ephemeris, elements);
                Err(Error::NBodyBackendUnavailable { backend: "assist" })
            }
        })
    }

    #[cfg(feature = "nbody")]
    fn curated_initial_state(
        &self,
        ephemeris: &AssistEphemeris,
        body: i32,
        epoch: f64,
    ) -> Result<Option<CartesianState>> {
        let Some((_, target)) = crate::constants::get_spk_body_info_from_map(body) else {
            return Ok(None);
        };
        // Initialization and forcing share the native bounded file owner;
        // a separate registry would materialize the whole SPK for one seed.
        ephemeris.inner.curated_initial_state(target, epoch)
    }

    /// Sample a trajectory with retained ASSIST resources when requested.
    ///
    /// Epochs are TT Julian days. A nonpositive sample count gives no samples;
    /// one gives the initial epoch. Missing or unsuitable ASSIST data follow the
    /// documented Sun-only fallback; other native failures propagate.
    #[allow(clippy::too_many_arguments)]
    pub fn trajectory(
        &self,
        elements: &OrbitalElements,
        jd_start: f64,
        jd_end: f64,
        count: i64,
        use_assist: bool,
        integrator: ReboundIntegrator,
    ) -> Result<PropagationOutcome<Vec<PropagationResult>>> {
        let count = usize::try_from(count.max(0))
            .map_err(|_| invalid("trajectory sample count is too large"))?;
        trajectory_with(
            elements,
            jd_start,
            jd_end,
            count,
            use_assist,
            integrator,
            |times| {
                self.with_ephemeris(|ephemeris| ephemeris.trajectory(elements, jd_start, times))
            },
        )
    }

    /// Compare a propagated and secular Keplerian position using these retained
    /// files. Units, diagnostics and fallback rules match [`compare_with_keplerian`].
    pub fn compare(
        &self,
        elements: &OrbitalElements,
        jd_tt: f64,
        use_assist: bool,
    ) -> Result<PropagationOutcome<NbodyComparison>> {
        comparison_with(elements, jd_tt, use_assist, || {
            self.propagate(elements, elements.epoch(), jd_tt, None)
        })
    }
}

/// Apply one apparent-place operation to a covered, retained source trajectory.
/// The Sun and target use the same TT-as-TDB epoch convention. No finite
/// emission or speed sample is replaced with a held position/velocity.
#[cfg(feature = "nbody")]
fn stored_trajectory<T>(
    ephemeris: &backend::Ephemeris,
    reader: &mut crate::spk::PositionalChebyshev,
    target: i32,
    jd_tt: f64,
    reception: crate::spk::ChebState,
    operation: impl FnOnce(
        &PropagationResult,
        &mut BarycentricSampler<'_>,
        &SourceSunSampler<'_>,
    ) -> Result<T>,
) -> std::result::Result<PropagationOutcome<T>, OwnedTrajectoryError> {
    let time = jd_tt - ephemeris.reference_jd();
    let sun_at = |offset: f64| {
        finite(offset)?;
        let time = time + offset;
        finite(time)?;
        ephemeris.state(0, time)?.check()
    };
    let central_sun = sun_at(0.0).map_err(OwnedTrajectoryError::Initial)?;
    let heliocentric = CartesianState {
        position: reception.position_au(),
        velocity: reception.velocity_au_per_day(),
    }
    .check()
    .map_err(OwnedTrajectoryError::Initial)?;
    let central_barycentric = shifted(heliocentric, central_sun, true);
    let central = PropagationResult::from_state(rotated(heliocentric, false), jd_tt)
        .map_err(OwnedTrajectoryError::Initial)?;
    let mut sample = |offset: f64| {
        finite(offset)?;
        if offset == 0.0 {
            return Ok(central_barycentric);
        }
        let (state, descriptor) = reader.trajectory_state(target, jd_tt, offset)?;
        if descriptor.center != crate::spk::CENTER_SUN {
            return Err(Error::SpkAmbiguousCenter { target });
        }
        shifted(
            CartesianState {
                position: state.position_au(),
                velocity: state.velocity_au_per_day(),
            },
            sun_at(offset)?,
            true,
        )
        .check()
    };
    operation(&central, &mut sample, &sun_at)
        .map(|value| PropagationOutcome {
            value,
            messages: Vec::new(),
            warnings: Vec::new(),
        })
        .map_err(OwnedTrajectoryError::after_acquisition)
}

static DEFAULT_ASSIST_SESSION: OnceLock<Mutex<Option<Arc<AssistSession>>>> = OnceLock::new();

/// Reuse the current discovered default data selection across calls.
///
/// Directory discovery is refreshed on each call. A changed selection creates a
/// new session while existing callers retain their old mapping through `Arc`.
/// This never opens or downloads a file until an integration needs it.
pub fn default_assist_session() -> Arc<AssistSession> {
    let config = AssistEphemConfig::default().discover();
    let mut current = DEFAULT_ASSIST_SESSION
        .get_or_init(|| Mutex::new(None))
        .lock()
        .unwrap_or_else(|error| error.into_inner());
    if let Some(session) = current.as_ref() {
        if session.config == config {
            return Arc::clone(session);
        }
    }
    let session = Arc::new(AssistSession::new(config));
    *current = Some(Arc::clone(&session));
    session
}

/// the idle default ASSIST cache. Active `Arc` owners retain their
/// mappings until their synchronous calculations finish; the next new caller
/// discovers and lazily opens its configured files again.
pub fn release_default_assist_session() {
    if let Some(cache) = DEFAULT_ASSIST_SESSION.get() {
        cache
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .take();
    }
}

#[cfg(feature = "nbody")]
fn assist_initial<'a>(
    elements: &OrbitalElements,
    jd_start: f64,
    ephemeris: &'a backend::Ephemeris,
    non_gravitational: Option<NonGravitationalParameters>,
    self_target: Option<i32>,
) -> Result<(backend::Simulation<'a>, CartesianState, Vec<Warning>)> {
    assist_initial_state(
        assist_element_state(elements, jd_start)?,
        jd_start,
        ephemeris,
        non_gravitational,
        self_target.or(elements.physical_spice_id()),
    )
}

#[cfg(feature = "nbody")]
fn assist_element_state(elements: &OrbitalElements, jd_start: f64) -> Result<CartesianState> {
    backend::Simulation::two_body(
        elements_to_rebound_orbit(elements, jd_start)?,
        ReboundIntegrator::Ias15,
        None,
    )?
    .state(1)?
    .check()
}

#[cfg(feature = "nbody")]
fn assist_initial_state<'a>(
    initial: CartesianState,
    jd_start: f64,
    ephemeris: &'a backend::Ephemeris,
    non_gravitational: Option<NonGravitationalParameters>,
    self_target: Option<i32>,
) -> Result<(backend::Simulation<'a>, CartesianState, Vec<Warning>)> {
    let initial_time = jd_start - ephemeris.reference_jd();
    finite(initial_time)?;
    let barycentric = shifted(
        rotated(initial, true),
        ephemeris.state(0, initial_time)?,
        true,
    );
    let (simulation, warnings) = assist_barycentric_initial(
        barycentric,
        initial_time,
        ephemeris,
        non_gravitational,
        self_target,
    )?;
    Ok((simulation, initial, warnings))
}

#[cfg(feature = "nbody")]
fn assist_barycentric_initial<'a>(
    barycentric: CartesianState,
    initial_time: f64,
    ephemeris: &'a backend::Ephemeris,
    non_gravitational: Option<NonGravitationalParameters>,
    self_target: Option<i32>,
) -> Result<(backend::Simulation<'a>, Vec<Warning>)> {
    let mut simulation = backend::Simulation::assist(ephemeris)?;
    simulation.add(barycentric);
    #[cfg(test)]
    owned_input_observations::note(|| {
        serde_json::json!({
            "kind":"simulation_input", "relative_time":initial_time,
            "ssb_icrf":owned_input_observations::pv(barycentric)
        })
    });
    simulation.set_time(initial_time);
    let mut warnings = Vec::new();
    // Preserve standard N16 warning names, binding them to physical SPICE
    // identities. Other selected files may reorder or add asteroid records.
    const PERTURBERS: [(i32, &str); 16] = [
        (2_000_107, "Camilla"),
        (2_000_001, "Ceres"),
        (2_000_065, "Cybele"),
        (2_000_511, "Davida"),
        (2_000_015, "Eunomia"),
        (2_000_031, "Euphrosyne"),
        (2_000_052, "Europa"),
        (2_000_010, "Hygiea"),
        (2_000_704, "Interamnia"),
        (2_000_007, "Iris"),
        (2_000_003, "Juno"),
        (2_000_002, "Pallas"),
        (2_000_016, "Psyche"),
        (2_000_087, "Sylvia"),
        (2_000_088, "Thisbe"),
        (2_000_004, "Vesta"),
    ];
    // A small separation is a real encounter for an arbitrary test particle.
    // Exclude only a physical target selected by the caller, never a neighbour
    // inferred from its position or catalogue-element resemblance.
    if let Some((index, physical_id)) = self_perturber(self_target, ephemeris.asteroid_targets()) {
        simulation.exclude_asteroid(index);
        warnings.push(Warning::ReboundSelfInteraction {
            perturber: PERTURBERS
                .iter()
                .find(|(target, _)| *target == physical_id)
                .map_or_else(
                    || format!("SPICE {physical_id}"),
                    |(_, name)| (*name).into(),
                ),
        });
    }
    if let Some(parameters) = non_gravitational {
        let parameters = [parameters.radial, parameters.transverse, parameters.normal];
        for value in parameters {
            finite(value)?;
        }
        simulation.non_gravitational_parameters(parameters);
    }
    Ok((simulation, warnings))
}

/// Match physical identity without inspecting positions or encounter distance.
#[cfg(feature = "nbody")]
fn self_perturber(
    target: Option<i32>,
    selected: impl IntoIterator<Item = (i32, i32)>,
) -> Option<(i32, i32)> {
    let target = target?;
    selected
        .into_iter()
        .find(|&(_, physical_id)| physical_id == target)
}

/// Same retained simulation, self identity and force convention as its PV.
#[cfg(feature = "nbody")]
fn endpoint_acceleration(
    simulation: &mut backend::Simulation<'_>,
    ephemeris: &backend::Ephemeris,
    time: f64,
) -> Result<[f64; 3]> {
    let target = simulation.acceleration()?;
    let sun = ephemeris.sun_acceleration(time)?;
    let acceleration = equatorial_to_ecliptic(std::array::from_fn(|axis| target[axis] - sun[axis]));
    for component in acceleration {
        finite(component)?;
    }
    Ok(acceleration)
}

/// Propagate through ASSIST's JPL-driven force model.
///
/// Input/output frame and units follow [`PropagationResult`]. Optional Marsden
/// parameters are AU/day². Planetary data are required; asteroid perturbers are
/// optional. Existing immutable data files are opened read-only by the Rust engine and
/// must not be truncated or changed while a calculation is active. Arbitrary
/// element/state calls retain every configured perturber. Body-selected calls
/// exclude only that target's physical self-force and retain its advisory.
pub fn propagate_orbit_assist(
    elements: &OrbitalElements,
    jd_start: f64,
    jd_end: f64,
    config: Option<&AssistEphemConfig>,
    non_gravitational: Option<NonGravitationalParameters>,
) -> Result<PropagationResult> {
    Ok(
        propagate_orbit_assist_with_warnings(
            elements,
            jd_start,
            jd_end,
            config,
            non_gravitational,
        )?
        .value,
    )
}

/// ASSIST propagation retaining native diagnostics and self-interaction warnings.
pub fn propagate_orbit_assist_with_warnings(
    elements: &OrbitalElements,
    jd_start: f64,
    jd_end: f64,
    config: Option<&AssistEphemConfig>,
    non_gravitational: Option<NonGravitationalParameters>,
) -> Result<PropagationOutcome<PropagationResult>> {
    AssistEphemeris::open(config)?.propagate(elements, jd_start, jd_end, non_gravitational)
}

/// Equally spaced scalar TT Julian days including both endpoints.
/// Zero samples yields an empty vector; one sample yields the starting epoch.
/// Nonfinite epochs, an overflowing span or an unallocatable count are errors.
pub fn trajectory_epochs(jd_start: f64, jd_end: f64, count: usize) -> Result<Vec<f64>> {
    finite(jd_start)?;
    finite(jd_end)?;
    let mut times = Vec::new();
    times
        .try_reserve_exact(count)
        .map_err(|_| invalid("trajectory sample count cannot be allocated"))?;
    let span = jd_end - jd_start;
    finite(span)?;
    let step = if count > 1 {
        span / (count - 1) as f64
    } else {
        0.0
    };
    for index in 0..count {
        let epoch = if count > 1 && index == count - 1 {
            jd_end
        } else {
            jd_start + index as f64 * step
        };
        finite(epoch)?;
        times.push(epoch);
    }
    Ok(times)
}

/// Whether an ASSIST data/domain failure permits the documented Sun-only or
/// analytical fallback. Native selected-data, state and integration failures
/// are terminal, including code one; diagnostics/status remain unchanged.
pub fn assist_fallback(error: &Error) -> bool {
    matches!(
        error.primary_error(),
        Error::NBodyBackendUnavailable { .. }
            | Error::InvalidOrbitElements { .. }
            | Error::NonFiniteArgument { .. }
            | Error::EphemerisLoaderIo {
                kind: std::io::ErrorKind::NotFound,
                ..
            }
    )
}

/// Propagate a sampled trajectory while retaining one native simulation.
///
/// Uses ASSIST when requested and its data/model can serve the interval. Its
/// ordinary data/domain failures fall back to the documented Sun-only REBOUND
/// route; this adds no alternative numerical solver. Returned states keep their
/// actual sampled epochs. Fixed-step selectors use the Rust engine's default
/// timestep here, matching trajectory initialization rather than scalar setup.
pub fn propagate_trajectory(
    elements: &OrbitalElements,
    jd_start: f64,
    jd_end: f64,
    count: usize,
    use_assist: bool,
    config: Option<&AssistEphemConfig>,
    integrator: ReboundIntegrator,
) -> Result<Vec<PropagationResult>> {
    Ok(propagate_trajectory_with_warnings(
        elements, jd_start, jd_end, count, use_assist, config, integrator,
    )?
    .value)
}

/// Sampled propagation with all native messages and physical warnings retained.
/// Resources and integration stay owned until every output sample is complete.
pub fn propagate_trajectory_with_warnings(
    elements: &OrbitalElements,
    jd_start: f64,
    jd_end: f64,
    count: usize,
    use_assist: bool,
    config: Option<&AssistEphemConfig>,
    integrator: ReboundIntegrator,
) -> Result<PropagationOutcome<Vec<PropagationResult>>> {
    trajectory_with(
        elements,
        jd_start,
        jd_end,
        count,
        use_assist,
        integrator,
        |times| {
            AssistEphemeris::open(config)
                .and_then(|ephemeris| ephemeris.trajectory(elements, jd_start, times))
        },
    )
}

#[allow(clippy::too_many_arguments)]
fn trajectory_with(
    elements: &OrbitalElements,
    jd_start: f64,
    jd_end: f64,
    count: usize,
    use_assist: bool,
    integrator: ReboundIntegrator,
    assist: impl FnOnce(&[f64]) -> Result<PropagationOutcome<Vec<PropagationResult>>>,
) -> Result<PropagationOutcome<Vec<PropagationResult>>> {
    let times = trajectory_epochs(jd_start, jd_end, count)?;
    if times.is_empty() {
        return Ok(PropagationOutcome {
            value: Vec::new(),
            messages: Vec::new(),
            warnings: Vec::new(),
        });
    }
    #[cfg(feature = "nbody")]
    {
        if use_assist {
            match assist(&times) {
                Ok(result) => return Ok(result),
                Err(error) if assist_fallback(&error) => {}
                Err(error) => return Err(error),
            }
        }
        let mut simulation = backend::Simulation::two_body(
            elements_to_rebound_orbit(elements, jd_start)?,
            integrator,
            None,
        )?;
        let mut values = Vec::with_capacity(times.len());
        let mut messages = Vec::new();
        for time in times {
            messages.extend(simulation.integrate(time - jd_start)?);
            values.push(PropagationResult::from_state(simulation.state(1)?, time)?);
        }
        Ok(PropagationOutcome {
            value: values,
            messages,
            warnings: Vec::new(),
        })
    }
    #[cfg(not(feature = "nbody"))]
    {
        let _ = (elements, use_assist, assist, integrator);
        Err(Error::NBodyBackendUnavailable { backend: "rebound" })
    }
}

/// Cartesian and spherical state used by the Keplerian comparison.
#[derive(Clone, Copy, Debug)]
pub struct ComparisonPoint {
    /// Heliocentric position in AU.
    pub position: [f64; 3],
    /// Ecliptic longitude, degrees.
    pub longitude: f64,
    /// Ecliptic latitude, degrees.
    pub latitude: f64,
    /// Heliocentric distance, AU.
    pub distance: f64,
}

fn comparison_point(position: [f64; 3], jd_tt: f64) -> Result<ComparisonPoint> {
    let result = PropagationResult::from_state(
        CartesianState {
            position,
            velocity: [0.0; 3],
        },
        jd_tt,
    )?;
    Ok(ComparisonPoint {
        position,
        longitude: result.ecliptic_lon(),
        latitude: result.ecliptic_lat(),
        distance: result.distance(),
    })
}

/// Integrated and Keplerian positions with physical differences.
#[derive(Clone, Copy, Debug)]
pub struct NbodyComparison {
    /// Keplerian position.
    pub keplerian: ComparisonPoint,
    /// Native integrated position.
    pub nbody: ComparisonPoint,
    /// Longitude difference, arcseconds, wrapped to one half-turn.
    pub delta_lon_arcsec: f64,
    /// Latitude difference, arcseconds.
    pub delta_lat_arcsec: f64,
    /// Cosine-latitude tangent-plane separation, arcseconds.
    pub angular_sep_arcsec: f64,
    /// Heliocentric radius difference, AU.
    pub delta_dist_au: f64,
    /// Requested comparison mode: assist or rebound.
    pub method: &'static str,
    /// Integration interval, days.
    pub propagation_days: f64,
}

/// Compare a native propagation from the element epoch with its secular model.
///
/// Differences use ecliptic-J2000 angles and AU radii. The small-angle separation
/// applies the cosine of mean latitude. As in sampled propagation, unavailable
/// ASSIST data or a model-domain failure uses the documented Sun-only route.
pub fn compare_with_keplerian(
    elements: &OrbitalElements,
    jd_tt: f64,
    use_assist: bool,
    config: Option<&AssistEphemConfig>,
) -> Result<NbodyComparison> {
    Ok(compare_with_keplerian_with_warnings(elements, jd_tt, use_assist, config)?.value)
}

/// Keplerian comparison retaining integration and physical diagnostics.
pub fn compare_with_keplerian_with_warnings(
    elements: &OrbitalElements,
    jd_tt: f64,
    use_assist: bool,
    config: Option<&AssistEphemConfig>,
) -> Result<PropagationOutcome<NbodyComparison>> {
    comparison_with(elements, jd_tt, use_assist, || {
        propagate_orbit_assist_with_warnings(elements, elements.epoch(), jd_tt, config, None)
    })
}

fn comparison_with(
    elements: &OrbitalElements,
    jd_tt: f64,
    use_assist: bool,
    assist: impl FnOnce() -> Result<PropagationOutcome<PropagationResult>>,
) -> Result<PropagationOutcome<NbodyComparison>> {
    let keplerian = comparison_point(
        calc_minor_body_position(elements, jd_tt, true, None)?,
        jd_tt,
    )?;
    let rebound = || {
        propagate_orbit_rebound_with_warnings(
            elements,
            elements.epoch(),
            jd_tt,
            ReboundIntegrator::Ias15,
            None,
        )
    };
    let propagated = if use_assist && check_assist_available() {
        match assist() {
            Ok(result) => result,
            Err(error) if assist_fallback(&error) => rebound()?,
            Err(error) => return Err(error),
        }
    } else {
        rebound()?
    };
    let value = propagated.value;
    let nbody = comparison_point([value.x, value.y, value.z], jd_tt)?;
    let mut delta_lon = nbody.longitude - keplerian.longitude;
    if delta_lon > 180.0 {
        delta_lon -= 360.0;
    } else if delta_lon < -180.0 {
        delta_lon += 360.0;
    }
    let delta_lon = delta_lon * 3600.0;
    let delta_lat = (nbody.latitude - keplerian.latitude) * 3600.0;
    let middle_latitude = ((keplerian.latitude + nbody.latitude) / 2.0).to_radians();
    Ok(PropagationOutcome {
        value: NbodyComparison {
            keplerian,
            nbody,
            delta_lon_arcsec: delta_lon,
            delta_lat_arcsec: delta_lat,
            angular_sep_arcsec: (delta_lon * middle_latitude.cos()).hypot(delta_lat),
            delta_dist_au: nbody.distance - keplerian.distance,
            method: if use_assist { "assist" } else { "rebound" },
            propagation_days: jd_tt - elements.epoch(),
        },
        messages: propagated.messages,
        warnings: propagated.warnings,
    })
}

/// Download explicitly requested ASSIST files through the native download layer.
/// Existing files are validated and reused; nothing is downloaded implicitly.
pub use crate::download::download_assist_data;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn optional_assist_keeps_unavailable_acquisition_eligible() {
        for error in [
            Error::NBodyBackendUnavailable { backend: "assist" },
            Error::EphemerisLoaderIo {
                path: "/missing/selected-force-file".into(),
                kind: std::io::ErrorKind::NotFound,
                errno: None,
                message: "missing selected force file".into(),
            },
        ] {
            assert!(assist_fallback(&error), "{error:?}");
        }
    }

    #[test]
    fn optional_assist_never_downgrades_terminal_integration_or_force_failure() {
        for (stage, diagnostic) in [
            (2, "selected force header is corrupt"),
            (3, "selected force record is corrupt"),
            (4, "native solver did not converge"),
            (5, "native calculation state is unavailable"),
        ] {
            let error = Error::NBodyNative {
                stage,
                code: 1,
                diagnostic: diagnostic.into(),
            }
            .with_cause(Error::MathDomain { value: 17.0 });
            assert!(!assist_fallback(&error));
            assert!(matches!(
                error.primary_error(),
                Error::NBodyNative { stage: actual_stage, code: 1, diagnostic: actual }
                    if *actual_stage == stage && actual == diagnostic
            ));
            assert!(error.cause().is_some(), "native cause remains attached");
        }
    }

    #[cfg(feature = "nbody")]
    #[test]
    fn self_force_selection_requires_explicit_physical_identity() {
        // These ordinals deliberately differ from N16 catalogue order.
        let selected = [(11, 2_000_004), (12, 2_000_001), (13, 2_000_088)];
        assert_eq!(self_perturber(None, selected), None);
        assert_eq!(self_perturber(Some(2_000_999), selected), None);
        assert_eq!(
            self_perturber(Some(2_000_001), selected),
            Some((12, 2_000_001))
        );
        assert_eq!(
            self_perturber(Some(2_000_088), selected),
            Some((13, 2_000_088))
        );
    }

    #[test]
    fn frame_rotation_round_trip_and_units() {
        let original = [0.3, -0.5, 0.1];
        let restored = equatorial_to_ecliptic(ecliptic_to_equatorial(original));
        for (a, b) in original.into_iter().zip(restored) {
            assert!((a - b).abs() < 1e-15);
        }
        assert_eq!(
            trajectory_epochs(10.0, 20.0, 3).unwrap(),
            [10.0, 15.0, 20.0]
        );
        assert_eq!(trajectory_epochs(10.0, 20.0, 1).unwrap(), [10.0]);
        assert!(trajectory_epochs(10.0, 20.0, 0).unwrap().is_empty());
    }

    #[test]
    fn lazy_session_preserves_explicit_paths_and_empty_sampling_needs_no_data() {
        let config = AssistEphemConfig {
            planets_file: Some("/libephemeris-test-missing-planets.440".into()),
            asteroids_file: Some("/libephemeris-test-missing-asteroids.bsp".into()),
            data_dir: None,
        };
        let session = Arc::new(AssistSession::new(config.clone()));
        assert_eq!(session.config(), &config);
        let shared = Arc::clone(&session);
        let elements = crate::minor_bodies::minor_body_elements(crate::constants::CERES).unwrap();
        let empty = shared
            .trajectory(
                &elements,
                2451545.0,
                2451546.0,
                -1,
                true,
                ReboundIntegrator::Ias15,
            )
            .unwrap();
        assert!(empty.value.is_empty());
        assert!(session.ephemeris.lock().unwrap().is_none());
        assert_eq!(Arc::strong_count(&session), 2);
    }

    #[test]
    fn nonfinite_timestep_is_rejected_before_a_zero_interval() {
        let elements = crate::minor_bodies::minor_body_elements(crate::constants::CERES).unwrap();
        assert!(matches!(
            propagate_orbit_rebound(
                &elements,
                2451545.0,
                2451545.0,
                ReboundIntegrator::Ias15,
                Some(f64::NAN)
            ),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[cfg(feature = "nbody")]
    fn admitted_assist_test_config() -> AssistEphemConfig {
        use sha2::{Digest, Sha256};
        use std::io::Read;

        let directory = std::env::var_os("LIBEPHEMERIS_TEST_DATA_DIR")
            .map(std::path::PathBuf::from)
            .unwrap_or_else(|| std::path::PathBuf::from("/Volumes/ExFAT/libephemeris"))
            .join("assist");
        for (name, expected) in [
            (
                "linux_p1550p2650.440",
                "29915576d0a6555766b99485ac3056ee415e86df4fce282611c31afb329ad062",
            ),
            (
                "sb441-n16.bsp",
                "919d612ce3c72a78fc7158f9120156542d0f21e6b8b052e4c1339c759747fd90",
            ),
        ] {
            let path = directory.join(name);
            let mut file = std::fs::File::open(&path).expect("required local ASSIST fixture");
            let mut digest = Sha256::new();
            let mut buffer = [0u8; 65_536];
            loop {
                let length = file.read(&mut buffer).unwrap();
                if length == 0 {
                    break;
                }
                digest.update(&buffer[..length]);
            }
            let actual: String = digest
                .finalize()
                .iter()
                .map(|byte| format!("{byte:02x}"))
                .collect();
            assert_eq!(actual, expected, "{}", path.display());
        }
        AssistEphemConfig {
            planets_file: Some(directory.join("linux_p1550p2650.440")),
            asteroids_file: Some(directory.join("sb441-n16.bsp")),
            data_dir: None,
        }
    }

    #[cfg(feature = "nbody")]
    #[test]
    #[ignore = "explicit same-initial-state integration control; requires hashed ASSIST files and an output path"]
    fn eris_same_initial_state_native_replay_control() {
        let config = admitted_assist_test_config();
        let ephemeris = AssistEphemeris::open(Some(&config)).unwrap();
        let elements = crate::minor_bodies::minor_body_elements(crate::constants::ERIS).unwrap();
        let start = elements.epoch();
        let end = 2424473.5;
        let initial = assist_element_state(&elements, start).unwrap();
        let time = end - ephemeris.inner.reference_jd();
        let sun = ephemeris.inner.state(0, time).unwrap();
        let mut rows = Vec::new();
        let mut initial_words = None;
        for replay in 0..3 {
            let (mut simulation, retained, warnings) =
                assist_initial_state(initial, start, &ephemeris.inner, None, None).unwrap();
            assert_eq!(retained, initial);
            assert!(warnings.is_empty());
            let barycentric = simulation.state(0).unwrap();
            let words: Vec<_> = barycentric
                .position
                .into_iter()
                .chain(barycentric.velocity)
                .map(f64::to_bits)
                .collect();
            if let Some(expected) = &initial_words {
                assert_eq!(&words, expected);
            } else {
                initial_words = Some(words.clone());
            }
            let messages = simulation.integrate(time).unwrap();
            let barycentric_end = simulation.state(0).unwrap().check().unwrap();
            let final_words: Vec<_> = barycentric_end
                .position
                .into_iter()
                .chain(barycentric_end.velocity)
                .map(f64::to_bits)
                .collect();
            let heliocentric = rotated(shifted(barycentric_end, sun, false), false);
            rows.push(serde_json::json!({
                "replay": replay,
                "initial_barycentric_words": words,
                "final_barycentric_words": final_words,
                "final_barycentric_position": barycentric_end.position,
                "final_barycentric_velocity": barycentric_end.velocity,
                "final_ecliptic_position": heliocentric.position,
                "final_ecliptic_velocity": heliocentric.velocity,
                "messages": messages,
            }));
        }
        for pair in rows.windows(2) {
            assert_eq!(
                pair[0]["final_barycentric_words"],
                pair[1]["final_barycentric_words"]
            );
        }
        // The retired C backend's epsilon knob has no counterpart here. This
        // control measures native replay from identical initial words, not
        // sensitivity to a fictitious setting or a global integration bound.
        let output = std::env::var_os("LIBEPHEMERIS_ERIS_ACCURACY_OUTPUT")
            .expect("explicit diagnostic output path");
        let result = serde_json::json!({
            "complete": true,
            "public_qualified": false,
            "scope": "same-native-initial-state replay only; no absolute accuracy bound",
            "jd_start": start,
            "jd_end": end,
            "source_body": crate::constants::ERIS,
            "rows": rows,
        });
        let file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(output)
            .unwrap();
        serde_json::to_writer_pretty(file, &result).unwrap();
        let _ = admitted_assist_test_config();
    }

    #[cfg(feature = "nbody")]
    #[test]
    #[ignore = "requires the hash-verified local DE440 and sb441-n16 ASSIST files"]
    fn curated_ceres_excludes_only_its_own_asteroid_force() {
        let config = admitted_assist_test_config();
        let elements = crate::minor_bodies::minor_body_elements(crate::constants::CERES).unwrap();
        let epoch = elements.epoch();
        let ephemeris = AssistEphemeris::open(Some(&config)).unwrap();
        // The generic element-taking API retains this actual catalogue
        // object's identity without replacing its supplied initial elements.
        let result = ephemeris.propagate(&elements, epoch, epoch, None).unwrap();
        let initial = result.value;
        assert!(
            [
                initial.x, initial.y, initial.z, initial.vx, initial.vy, initial.vz
            ]
            .into_iter()
            .all(f64::is_finite),
            "initial={initial:?}"
        );
        assert_eq!(initial.jd_tt, epoch);
        let time = epoch - ephemeris.inner.reference_jd();
        let sun = ephemeris.inner.state(0, time).unwrap();
        let barycentric = shifted(
            rotated(
                CartesianState {
                    position: [initial.x, initial.y, initial.z],
                    velocity: [initial.vx, initial.vy, initial.vz],
                },
                true,
            ),
            sun,
            true,
        );
        let ceres = ephemeris.inner.position(12, time).unwrap_or_else(|error| {
            panic!("JPL Ceres state failed: {error}; initial={initial:?}; barycentric={barycentric:?}; warnings={:?}", result.warnings)
        });
        let displacement: [f64; 3] =
            std::array::from_fn(|index| barycentric.position[index] - ceres[index]);
        let separation = displacement[0]
            .hypot(displacement[1])
            .hypot(displacement[2]);
        assert!(
            separation < 1e-4,
            "separation={separation:.17e} AU; initial={initial:?}; barycentric={barycentric:?}; JPL_Ceres={ceres:?}; warnings={:?}",
            result.warnings
        );
        assert!(result.warnings.iter().any(|warning| matches!(warning, Warning::ReboundSelfInteraction { perturber } if perturber == "Ceres")),
            "coincident Ceres force remained active: separation={separation:.17e} AU; initial={initial:?}; warnings={:?}", result.warnings);

        // At the exact Ceres position its self-force is singular. Excluding
        // Ceres must leave a finite, nonzero acceleration from the other 15
        // asteroids, rather than turning the whole force family off.
        let mut simulation = backend::Simulation::assist(&ephemeris.inner).unwrap();
        simulation.add(CartesianState {
            position: ceres,
            velocity: [0.0; 3],
        });
        simulation.set_time(time);
        simulation.exclude_asteroid(12);
        let acceleration = simulation.asteroid_acceleration().unwrap();
        assert!(acceleration.into_iter().all(f64::is_finite));
        assert!(acceleration.into_iter().any(|value| value.abs() > 1e-15));
        let mut expected = [0.0; 3];
        let mut absolute_sum = [0.0; 3];
        for index in (11..27).filter(|index| *index != 12) {
            let position = ephemeris.inner.position(index, time).unwrap();
            let gm = ephemeris
                .inner
                .gravitational_parameter(index, time)
                .unwrap();
            let displacement: [f64; 3] = std::array::from_fn(|i| ceres[i] - position[i]);
            let radius = displacement[0]
                .hypot(displacement[1])
                .hypot(displacement[2]);
            for axis in 0..3 {
                let contribution = -gm * displacement[axis] / radius.powi(3);
                expected[axis] += contribution;
                absolute_sum[axis] += contribution.abs();
            }
        }
        for axis in 0..3 {
            // Independent inverse-cube sum; this operation-count bound covers
            // the two binary64 norm/summation evaluations, not ephemeris error.
            let bound = 128.0 * f64::EPSILON * absolute_sum[axis];
            assert!((acceleration[axis] - expected[axis]).abs() <= bound);
        }

        let session = AssistSession::new(config);
        let selected = session
            .curated_initial_state(&ephemeris, crate::constants::CERES, epoch)
            .unwrap()
            .expect("Ceres is stored in the selected N16 force image");
        let at_epoch = session
            .propagate_curated(crate::constants::CERES, epoch)
            .unwrap();
        assert_eq!(
            at_epoch.value,
            PropagationResult::from_state(selected, epoch).unwrap()
        );
        assert_ne!(
            at_epoch.value, initial,
            "the curated source must use the actual force-image state"
        );
        assert!(
            session
                .curated_initial_state(&ephemeris, crate::constants::CHIRON, epoch)
                .unwrap()
                .is_none()
        );
        assert!(matches!(
            session.propagate_curated(-1, epoch),
            Err(Error::UnknownBody { .. })
        ));
    }

    #[cfg(feature = "nbody")]
    #[test]
    #[ignore = "requires the hash-verified local DE440 and sb441-n16 ASSIST files"]
    fn covered_curated_targets_sample_the_retained_source_curve() {
        let config = admitted_assist_test_config();
        let independent = AssistEphemeris::open(Some(&config)).unwrap();
        let session = Arc::new(AssistSession::new(config));
        // These are all curated targets actually present in the selected N16
        // image. Source membership, rather than this test inventory, controls
        // the runtime branch.
        for (body, index) in [
            (17, 12),
            (18, 22),
            (19, 21),
            (20, 26),
            (10010, 18),
            (10016, 23),
            (10052, 17),
            (10087, 24),
            (10511, 14),
            (10704, 19),
        ] {
            let elements = crate::minor_bodies::minor_body_elements(body).unwrap();
            for epoch in [2451545.0, 2460310.5] {
                let outcome = session
                    .with_barycentric_trajectory(body, &elements, epoch, |_, sample, _| {
                        for offset in [-0.01, 0.0, 0.01] {
                            let actual = sample(offset)?;
                            let expected = independent.inner.position(
                                index,
                                epoch - independent.inner.reference_jd() + offset,
                            )?;
                            for (a, b) in actual.position.into_iter().zip(expected) {
                                assert!(
                                    (a - b).abs() <= 64.0 * f64::EPSILON * b.abs().max(1.0),
                                    "body={body} epoch={epoch} offset={offset} native={a} SDK={b}"
                                );
                            }
                            assert!(actual.velocity.into_iter().all(f64::is_finite));
                        }
                        Ok(())
                    })
                    .unwrap();
                assert!(outcome.messages.is_empty() && outcome.warnings.is_empty());
            }
        }
        assert!(session.stored_trajectories.lock().unwrap().is_some());
        let supplier = crate::engine_api::NativeMinorAssist::with_session(Arc::clone(&session));
        let geometric =
            crate::engine_api::MinorAssistAcquisition::state_at(&supplier, 18, 2451545.0).unwrap();
        let expected = session
            .with_barycentric_trajectory(
                18,
                &crate::minor_bodies::minor_body_elements(18).unwrap(),
                2451545.0,
                |central, _, _| Ok(*central),
            )
            .unwrap()
            .value;
        assert_eq!(
            geometric, expected,
            "TRUEPOS and apparent calls must select the same orbit"
        );
        let (richer, acceleration) =
            crate::engine_api::MinorAssistAcquisition::state_with_acceleration_at(
                &supplier, 18, 2451545.0,
            )
            .unwrap();
        assert_eq!(
            richer, geometric,
            "endpoint enrichment must retain selected PV"
        );
        assert!(
            acceleration.is_none(),
            "a stored curve is not the integrated force model"
        );

        // A target absent from N16 instead uses one integrated trajectory.
        // Endpoint enrichment must preserve that route's actual central PV and
        // diagnostics while retaining its real force derivative.
        let elements = crate::minor_bodies::minor_body_elements(crate::constants::ERIS).unwrap();
        let central = session
            .with_barycentric_trajectory(
                crate::constants::ERIS,
                &elements,
                elements.epoch(),
                |central, _, _| Ok(*central),
            )
            .unwrap();
        let enriched = session
            .state_with_acceleration(crate::constants::ERIS, &elements, elements.epoch())
            .unwrap();
        assert_eq!(enriched.value.0, central.value);
        assert_eq!(enriched.messages, central.messages);
        assert_eq!(enriched.warnings, central.warnings);
        let acceleration = enriched
            .value
            .1
            .expect("integrated source has its own force");
        assert!(acceleration.into_iter().all(f64::is_finite));
        assert!(acceleration.into_iter().any(|component| component != 0.0));

        let elements = crate::minor_bodies::minor_body_elements(18).unwrap();
        let failure = session.with_barycentric_trajectory(18, &elements, 2451545.0, |_, _, _| {
            Err::<(), _>(Error::NonFiniteArgument { value: f64::NAN })
        });
        assert!(matches!(
            failure,
            Err(OwnedTrajectoryError::Emission {
                error: Error::NonFiniteArgument { .. },
                ..
            })
        ));
    }

    #[cfg(feature = "nbody")]
    #[test]
    fn direct_c_integrators_conserve_a_circular_orbit() {
        let motion = GAUSSIAN_K_SQUARED.sqrt().to_degrees();
        let elements = OrbitalElements::new(
            "test".into(),
            2451545.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            motion,
        )
        .unwrap();
        let initial = propagate_orbit_rebound(
            &elements,
            elements.epoch(),
            elements.epoch(),
            ReboundIntegrator::Ias15,
            None,
        )
        .unwrap();
        assert!((initial.x - 1.0).abs() < 1e-14);
        assert!(initial.y.abs() < 1e-14);
        for integrator in [
            ReboundIntegrator::Ias15,
            ReboundIntegrator::Whfast,
            ReboundIntegrator::Leapfrog,
            ReboundIntegrator::Mercurius,
            ReboundIntegrator::Trace,
        ] {
            let result = propagate_orbit_rebound(
                &elements,
                elements.epoch(),
                elements.epoch() + 30.0,
                integrator,
                Some(0.05),
            )
            .unwrap();
            let radius_error = (result.distance() - 1.0).abs();
            assert!(
                radius_error < 1e-6,
                "{} radius error {radius_error}",
                integrator.as_str()
            );
        }
        #[cfg(feature = "nbody")]
        {
            assert_eq!(
                get_rebound_version(),
                Some("libephemeris-native-radau15/4.0.0-experimental")
            );
            assert_eq!(
                get_assist_version(),
                Some("libephemeris-native-ppn-icrf-terrestrial-iau-solar/4.0.0-experimental")
            );
        }
    }
}

/// Call-scoped test observations; no state or query is manufactured here.
#[cfg(test)]
pub(crate) mod owned_input_observations {
    use serde_json::Value;
    use std::cell::RefCell;

    std::thread_local! {
        static EVENTS: RefCell<Option<Vec<Value>>> = const { RefCell::new(None) };
    }

    pub(crate) struct Scope(Option<Vec<Value>>);
    impl Scope {
        pub(crate) fn enter() -> Self {
            Self(EVENTS.with(|slot| slot.borrow_mut().replace(Vec::new())))
        }

        pub(crate) fn finish(self) -> Vec<Value> {
            EVENTS.with(|slot| slot.borrow_mut().take().unwrap())
        }
    }
    impl Drop for Scope {
        fn drop(&mut self) {
            EVENTS.with(|slot| *slot.borrow_mut() = self.0.take());
        }
    }

    pub(crate) fn note(event: impl FnOnce() -> Value) {
        EVENTS.with(|slot| {
            if let Some(events) = slot.borrow_mut().as_mut() {
                events.push(event());
            }
        });
    }

    pub(crate) fn pv(value: super::CartesianState) -> Value {
        let values = [
            value.position[0],
            value.position[1],
            value.position[2],
            value.velocity[0],
            value.velocity[1],
            value.velocity[2],
        ];
        serde_json::json!({"values":values, "words":values.map(|v|format!("{:016x}",v.to_bits()))})
    }
}
