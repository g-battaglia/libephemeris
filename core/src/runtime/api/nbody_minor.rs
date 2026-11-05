// SPDX-License-Identifier: AGPL-3.0
//! Native minor-body precision fallback over retained numerical resources.

use super::MinorAssistAcquisition;
use crate::config::Runtime;
use crate::context::EphemerisContext;
use crate::errors::{Error, Result, Warning};
use crate::minor_bodies;
use crate::nbody::{self, AssistSession, PropagationResult};
use crate::spk::SpkRegistry;
use std::sync::{Arc, Mutex, OnceLock};

/// Actual diagnostics emitted by the numerical fallback, retained for callers.
#[derive(Debug, Default)]
pub struct MinorNbodyDiagnostics {
    /// Nonfatal upstream integration messages.
    pub messages: Vec<String>,
    /// Physical advisories, including asteroid self-interaction avoidance.
    pub warnings: Vec<Warning>,
}

/// Native supplier for all 37 curated minor-body element records.
///
/// Reuses a shared, lazily mapped ASSIST session. An optional runtime can publish
/// an automatically downloaded target registration; the caller's configuration
/// snapshot and every already active file owner remain unchanged. Construction
/// performs no discovery, mapping or download; data selection is lazy.
pub struct NativeMinorAssist<'a> {
    session: OnceLock<Arc<AssistSession>>,
    runtime: Option<&'a Runtime>,
    diagnostics: Mutex<MinorNbodyDiagnostics>,
}

impl NativeMinorAssist<'static> {
    /// Defer the shared default ASSIST selection without enabling auto-SPK work.
    pub fn discover() -> Self {
        Self {
            session: OnceLock::new(),
            runtime: None,
            diagnostics: Mutex::new(MinorNbodyDiagnostics::default()),
        }
    }

    /// Retain an explicitly configured ASSIST session without auto-SPK work.
    pub fn with_session(session: Arc<AssistSession>) -> Self {
        Self {
            session: OnceLock::from(session),
            runtime: None,
            diagnostics: Mutex::new(MinorNbodyDiagnostics::default()),
        }
    }
}

impl<'a> NativeMinorAssist<'a> {
    /// Attach native auto-SPK provisioning to this call's numerical fallback.
    /// Cached ephemerides are shared; no runtime mutex survives integration.
    pub fn for_runtime(runtime: &'a Runtime) -> Self {
        Self {
            session: OnceLock::new(),
            runtime: Some(runtime),
            diagnostics: Mutex::new(MinorNbodyDiagnostics::default()),
        }
    }

    fn session(&self) -> &AssistSession {
        self.session.get_or_init(nbody::default_assist_session)
    }

    /// Keep the concrete native orbit owner through its emission evaluations.
    /// This private entry is never inferred from an opaque public callback.
    pub(crate) fn with_native_trajectory<T>(
        &self,
        body: i32,
        tt: f64,
        operation: impl FnOnce(
            &PropagationResult,
            &mut nbody::BarycentricSampler<'_>,
            &nbody::SourceSunSampler<'_>,
        ) -> Result<T>,
    ) -> std::result::Result<T, nbody::OwnedTrajectoryError> {
        let elements = minor_bodies::minor_body_elements(body).ok_or(
            nbody::OwnedTrajectoryError::Initial(Error::UnknownBody { body_id: body }),
        )?;
        self.retain_trajectory_result(
            self.session()
                .with_barycentric_trajectory(body, &elements, tt, operation),
        )
    }

    fn retain_trajectory_result<T>(
        &self,
        result: std::result::Result<nbody::PropagationOutcome<T>, nbody::OwnedTrajectoryError>,
    ) -> std::result::Result<T, nbody::OwnedTrajectoryError> {
        let outcome = match result {
            Ok(outcome) => outcome,
            Err(nbody::OwnedTrajectoryError::Emission {
                error,
                messages,
                warnings,
            }) => {
                self.record_diagnostics(messages, warnings);
                return Err(nbody::OwnedTrajectoryError::after_acquisition(error));
            }
            Err(error) => return Err(error),
        };
        self.record_diagnostics(outcome.messages, outcome.warnings);
        Ok(outcome.value)
    }

    fn record_diagnostics(&self, messages: Vec<String>, warnings: Vec<Warning>) {
        let mut diagnostics = self
            .diagnostics
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        diagnostics.messages.extend(messages);
        diagnostics.warnings.extend(warnings);
    }

    /// Drain a locally owned supplier after its synchronous borrowers finish.
    ///
    /// Exclusive access avoids initializing a platform mutex solely to take
    /// diagnostics. Poisoned data remains available as with the shared drain.
    pub(crate) fn drain_diagnostics_exclusive(&mut self) -> MinorNbodyDiagnostics {
        std::mem::take(
            self.diagnostics
                .get_mut()
                .unwrap_or_else(|error| error.into_inner()),
        )
    }

    /// Drain this call's owned messages and physical advisories once.
    pub fn drain_diagnostics(&self) -> MinorNbodyDiagnostics {
        std::mem::take(
            &mut *self
                .diagnostics
                .lock()
                .unwrap_or_else(|error| error.into_inner()),
        )
    }
}

/// Selected tier's auto-SPK interval, padded by 45 days and clipped to the
/// Horizons generation domain. The returned ISO dates are Gregorian midnight.
pub fn automatic_spk_interval(context: &EphemerisContext) -> (String, String) {
    let (start, end) = crate::context::spk_date_range_for_tier(context.tier());
    let start = crate::spk::iso_to_jd(start).expect("native tier dates are canonical ISO");
    let end = crate::spk::iso_to_jd(end).expect("native tier dates are canonical ISO");
    let minimum = crate::spk::iso_to_jd(minor_bodies::HORIZONS_SPK_DATE_MIN)
        .expect("native Horizons boundary is a canonical ISO date");
    let maximum = crate::spk::iso_to_jd(minor_bodies::HORIZONS_SPK_DATE_MAX)
        .expect("native Horizons boundary is a canonical ISO date");
    (
        crate::spk::jd_to_iso_date((start - 45.0).max(minimum)),
        crate::spk::jd_to_iso_date((end + 45.0).min(maximum)),
    )
}

fn automatic_spk_target(body: i32) -> Option<String> {
    if crate::constants::is_spk_auto_download_blocked(body) {
        None
    } else if let Some((target, _)) = crate::constants::get_spk_body_info_from_map(body) {
        Some(target.into())
    } else if crate::constants::AST_OFFSET < body && body < crate::constants::FIXSTAR_OFFSET {
        Some(format!("{};", body - crate::constants::AST_OFFSET))
    } else {
        None
    }
}

impl MinorAssistAcquisition for NativeMinorAssist<'_> {
    fn refreshed_spk(
        &self,
        body: i32,
        _tt: f64,
        context: &EphemerisContext,
    ) -> Result<Option<SpkRegistry>> {
        let Some(runtime) = self.runtime else {
            return Ok(None);
        };
        if !context.effective_auto_spk_download() {
            return Ok(None);
        }
        let Some(target) = automatic_spk_target(body) else {
            return Ok(None);
        };
        let (start, end) = automatic_spk_interval(context);
        match runtime
            .download_and_register_spk(&target, body, &start, &end, None, None, "@0", false, 300.0)
        {
            Ok(_) => {
                let registered = runtime.effective_context()?;
                let mut retained = SpkRegistry::new();
                retained.retain_registration_from(registered.spk_registry(), body)?;
                Ok(Some(retained))
            }
            // Provisioning is optional here. Its transport/data failure does
            // not certify a numerical result; the next actual provider decides.
            // In-coverage coefficient evaluation happens after this hook and
            // therefore cannot be caught by this download-only boundary.
            Err(error) => {
                crate::logging::emit(
                    crate::logging::LogLevel::Warning,
                    &format!("Auto SPK download failed for body {body}: {error}"),
                );
                Ok(None)
            }
        }
    }

    fn contains(&self, body: i32) -> Result<bool> {
        Ok(minor_bodies::minor_body_elements(body).is_some())
    }
    fn available(&self) -> Result<bool> {
        Ok(nbody::check_assist_available() && self.session().config().data_available())
    }
    #[cfg(feature = "nbody")]
    fn state_with_acceleration_at(
        &self,
        body: i32,
        tt: f64,
    ) -> Result<(PropagationResult, Option<[f64; 3]>)> {
        let elements =
            minor_bodies::minor_body_elements(body).ok_or(Error::UnknownBody { body_id: body })?;
        if !tt.is_finite() {
            // Direct supplier calls retain the same resource-before-epoch
            // ordering as the original PV-only acquisition.
            return self.state_at(body, tt).map(|state| (state, None));
        }
        self.retain_trajectory_result(self.session().state_with_acceleration(body, &elements, tt))
            .map_err(|error| match error {
                nbody::OwnedTrajectoryError::Initial(error)
                | nbody::OwnedTrajectoryError::Emission { error, .. } => error,
            })
    }

    fn state_at(&self, body: i32, tt: f64) -> Result<PropagationResult> {
        let elements =
            minor_bodies::minor_body_elements(body).ok_or(Error::UnknownBody { body_id: body })?;
        if !tt.is_finite() {
            // Preserve the direct Rust supplier's resource-before-epoch error
            // ordering. Public Runtime calls reject invalid epochs earlier.
            return self
                .session()
                .propagate(&elements, elements.epoch(), tt, None)
                .map(|outcome| outcome.value);
        }
        self.with_native_trajectory(body, tt, |central, _, _| Ok(*central))
            .map_err(|error| match error {
                nbody::OwnedTrajectoryError::Initial(error)
                | nbody::OwnedTrajectoryError::Emission { error, .. } => error,
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::context::Tier;

    #[test]
    fn direct_supplier_preserves_unknown_body_and_invalid_epoch_error_order() {
        let session = Arc::new(AssistSession::new(nbody::AssistEphemConfig {
            planets_file: Some("/libephemeris-missing-test-planets.440".into()),
            asteroids_file: Some("/libephemeris-missing-test-asteroids.bsp".into()),
            data_dir: None,
        }));
        let supplier = NativeMinorAssist::with_session(Arc::clone(&session));
        assert!(matches!(
            supplier.state_at(-1, f64::NAN),
            Err(Error::UnknownBody { body_id: -1 })
        ));
        assert!(matches!(
            supplier.state_with_acceleration_at(-1, f64::NAN),
            Err(Error::UnknownBody { body_id: -1 })
        ));
        let elements = minor_bodies::minor_body_elements(18).unwrap();
        let original = session
            .propagate(&elements, elements.epoch(), f64::NAN, None)
            .unwrap_err();
        let current = supplier.state_at(18, f64::NAN).unwrap_err();
        assert_eq!(
            std::mem::discriminant(&current),
            std::mem::discriminant(&original)
        );
        assert!(
            !matches!(current, Error::NonFiniteArgument { .. }),
            "acquisition still precedes invalid TT for direct supplier calls"
        );
        let enriched = supplier
            .state_with_acceleration_at(18, f64::NAN)
            .unwrap_err();
        assert_eq!(
            std::mem::discriminant(&enriched),
            std::mem::discriminant(&original),
            "endpoint enrichment preserves direct acquisition error ordering"
        );
    }

    #[test]
    fn emission_error_retains_already_produced_central_and_emission_diagnostics() {
        let supplier = NativeMinorAssist::discover();
        let result: std::result::Result<(), _> =
            supplier.retain_trajectory_result(Err(nbody::OwnedTrajectoryError::Emission {
                error: Error::MathDomain { value: 777.0 },
                messages: vec!["central".into(), "emission".into()],
                warnings: vec![Warning::ReboundSelfInteraction {
                    perturber: "synthetic".into(),
                }],
            }));
        assert!(matches!(
            result,
            Err(nbody::OwnedTrajectoryError::Emission {
                error: Error::MathDomain { value: 777.0 },
                ..
            })
        ));
        let diagnostics = supplier.drain_diagnostics();
        assert_eq!(diagnostics.messages, ["central", "emission"]);
        assert_eq!(
            diagnostics.warnings,
            [Warning::ReboundSelfInteraction {
                perturber: "synthetic".into()
            }]
        );
        assert!(supplier.drain_diagnostics().messages.is_empty());
    }

    #[test]
    fn exclusive_diagnostic_drain_keeps_order_once_only_and_poison_recovery() {
        for poison in [false, true] {
            let shared = NativeMinorAssist::discover();
            let mut exclusive = NativeMinorAssist::discover();
            for supplier in [&shared, &exclusive] {
                let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    let mut diagnostics = supplier.diagnostics.lock().unwrap();
                    diagnostics
                        .messages
                        .extend(["first".into(), "second".into()]);
                    diagnostics.warnings.extend([
                        Warning::ReboundSelfInteraction {
                            perturber: "first".into(),
                        },
                        Warning::ReboundSelfInteraction {
                            perturber: "second".into(),
                        },
                    ]);
                    if poison {
                        panic!("controlled diagnostic writer unwind");
                    }
                }));
                assert_eq!(outcome.is_err(), poison);
            }
            let expected = shared.drain_diagnostics();
            let actual = exclusive.drain_diagnostics_exclusive();
            assert_eq!(actual.messages, ["first", "second"]);
            assert_eq!(actual.messages, expected.messages);
            assert_eq!(actual.warnings, expected.warnings);
            assert_eq!(actual.warnings.len(), 2);
            let again = exclusive.drain_diagnostics_exclusive();
            assert!(again.messages.is_empty());
            assert!(again.warnings.is_empty());
            let shared_again = shared.drain_diagnostics();
            assert!(shared_again.messages.is_empty());
            assert!(shared_again.warnings.is_empty());
        }
    }

    #[test]
    fn curated_membership_and_tier_download_domains_are_native() {
        let config = nbody::AssistEphemConfig {
            planets_file: Some("/libephemeris-missing-planets.440".into()),
            asteroids_file: Some("/libephemeris-missing-asteroids.bsp".into()),
            data_dir: None,
        };
        let supplier = NativeMinorAssist::with_session(Arc::new(AssistSession::new(config)));
        let elements = minor_bodies::ordered_minor_body_elements();
        assert_eq!(elements.len(), 37);
        for (body, _) in elements {
            assert!(supplier.contains(body).unwrap());
        }
        assert!(!supplier.contains(crate::constants::MARS).unwrap());
        assert!(!supplier.available().unwrap());
        for tier in [Tier::Base, Tier::Medium, Tier::Extended] {
            let context = EphemerisContext::new().with_tier(tier);
            let (start, end) = automatic_spk_interval(&context);
            assert!(start.as_str() >= minor_bodies::HORIZONS_SPK_DATE_MIN);
            assert!(end.as_str() <= minor_bodies::HORIZONS_SPK_DATE_MAX);
            assert!(start < end);
        }
        assert_eq!(
            automatic_spk_target(crate::constants::AST_OFFSET + 433),
            Some("433".into())
        );
        assert_eq!(
            automatic_spk_target(crate::constants::AST_OFFSET + 90_123),
            Some("90123;".into())
        );
    }

    #[test]
    fn central_state_linear_stencil_preserves_velocity_and_epoch() {
        let state = PropagationResult::from_state(
            nbody::CartesianState {
                position: [1.2, -0.5, 0.3],
                velocity: [0.01, -0.02, 0.005],
            },
            2451545.0,
        )
        .unwrap();
        let shifted = state.shifted(2451545.5).unwrap();
        assert_eq!([shifted.vx, shifted.vy, shifted.vz], [0.01, -0.02, 0.005]);
        for (actual, expected) in [shifted.x, shifted.y, shifted.z]
            .into_iter()
            .zip([1.205_f64, -0.51, 0.3025])
        {
            assert!((actual - expected).abs() <= 2.0 * f64::EPSILON * expected.abs());
        }
        assert!(state.shifted(f64::NAN).is_err());
    }
}
