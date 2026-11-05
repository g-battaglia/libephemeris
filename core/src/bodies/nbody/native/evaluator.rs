// SPDX-License-Identifier: AGPL-3.0
//! Retained complete-force evaluation over one immutable selected dataset.
//!
//! Per-trajectory scratch and explicit physical self identity are separate
//! from shared ephemeris resources. Pole directions are supplied at each actual
//! epoch by the independently specified orientation model; this layer neither
//! guesses an orientation default nor learns one from a legacy force output.

use super::catalogue::HarmonicPoles;
use super::ephemeris::{Ephemeris, EphemerisError};
use super::forces::{ForceModel, Source};
use super::state::{NumericalError, State};
use std::sync::Arc;

/// Original data failures remain distinct from finite-state/force-domain errors.
/// This is the force callback's error type for Radau's typed-error route.
#[derive(Debug)]
pub(super) enum ForceError {
    Ephemeris(EphemerisError),
    /// An independently supplied orientation or split-epoch transformation failed.
    Frame(crate::errors::Error),
    Numerical(NumericalError),
}
impl From<EphemerisError> for ForceError {
    fn from(error: EphemerisError) -> Self {
        Self::Ephemeris(error)
    }
}
impl From<NumericalError> for ForceError {
    fn from(error: NumericalError) -> Self {
        Self::Numerical(error)
    }
}
impl std::fmt::Display for ForceError {
    fn fmt(&self, output: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Ephemeris(error) => write!(output, "{error}"),
            Self::Frame(error) => write!(output, "{error}"),
            Self::Numerical(error) => write!(output, "native force domain: {error:?}"),
        }
    }
}
impl std::error::Error for ForceError {}

/// One trajectory's force workspace, not a simulation-global mutable cache.
pub(super) struct Evaluator {
    ephemeris: Arc<Ephemeris>,
    model: ForceModel,
    sources: Vec<Source>,
    exclude: Option<i32>,
}

impl Evaluator {
    /// Local source crossing time in days, evaluated at the actual TDB epoch.
    /// This reuses the force workspace and selected ephemeris record cache.
    pub(super) fn maximum_step(
        &mut self,
        epoch: [f64; 2],
        state: State,
    ) -> Result<f64, ForceError> {
        state.checked()?;
        self.ephemeris.frame(epoch, &mut self.sources)?;
        Ok(super::encounter::maximum_step(
            state,
            &self.sources,
            self.exclude,
        )?)
    }

    /// Retain selected data and allocate reusable scratch once. All families,
    /// constants and optional A1/A2/A3 parameters are assembled before use.
    /// `exclude` is a physical source identifier, never an ordinal or distance.
    pub(super) fn new(
        ephemeris: Arc<Ephemeris>,
        initial_poles: HarmonicPoles,
        nongravitational: Option<[f64; 3]>,
        exclude: Option<i32>,
    ) -> Result<Self, ForceError> {
        let model = ephemeris
            .catalogue
            .force_model(initial_poles, nongravitational)?;
        let sources = Vec::with_capacity(ephemeris.particle_count());
        Ok(Self {
            ephemeris,
            model,
            sources,
            exclude,
        })
    }

    /// SSB ICRF acceleration in AU/day² at a two-part TDB Julian day.
    /// The supplied poles must describe that same epoch. Every nonself source
    /// and selected force family remains active. Data/domain failures stop this
    /// request; they never select another dataset, orientation or force model.
    pub(super) fn acceleration(
        &mut self,
        epoch: [f64; 2],
        state: State,
        poles: HarmonicPoles,
    ) -> Result<[f64; 3], ForceError> {
        state.checked()?;
        // Validate both new pole vectors before publishing either into the
        // retained model. The two slots are owned by Catalogue::force_model;
        // no external mutable model reference can change that invariant.
        let sun = self.model.zonals[0].with_pole(poles.sun)?;
        let earth = self.model.zonals[1].with_pole(poles.earth)?;
        self.ephemeris.frame(epoch, &mut self.sources)?;
        self.model.zonals[0] = sun;
        self.model.zonals[1] = earth;
        Ok(self
            .model
            .acceleration(state, &self.sources, self.exclude)?)
    }
}

#[cfg(test)]
#[path = "evaluator_tests.rs"]
mod tests;
