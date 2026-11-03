//! Explicit SPK-seeded N-body source for offline coefficient generation.

use super::{AssistSession, CartesianState, PropagationOutcome, finite, invalid};
use crate::errors::Result;
use crate::spk::GeometricSpkSource;
use std::sync::Arc;

/// Guard for IAS15 force stages at the source boundaries, in TDB days.
/// This preserves the established maximum 32-day LEB segment guard.
const INTEGRATION_GUARD_DAYS: f64 = 32.0;

/// Geometric barycentric ICRS states integrated from one explicit SPK seed.
///
/// Epochs are TDB Julian days; positions and velocities are AU and AU/day.
/// The seed must be a direct heliocentric J2000 SPK target at the chosen anchor.
/// Its Sun offset comes from the same retained ASSIST force ephemeris used by
/// integration. No runtime clock, catalogue orbit, alternate image or download
/// participates. Each batch starts independent forward/backward simulations;
/// sorting never changes the returned request order. This is a source model,
/// not a claim that reduced-force integration reproduces a complete JPL orbit
/// across the entire advertised source interval.
pub struct NbodySpkSource {
    session: Arc<AssistSession>,
    anchor: f64,
    initial: CartesianState,
    target: i32,
    coverage: (f64, f64),
}

impl NbodySpkSource {
    /// Retain a selected ASSIST session and sample an explicit SPK at the anchor.
    ///
    /// Both planet and asteroid force files are required. Their actual common
    /// continuous coverage must contain the anchor with a 32-day stage guard.
    /// Files must remain immutable. Invalid frames, seed centers, coverage,
    /// nonfinite inputs, missing files and native failures propagate as errors.
    /// Construction opens only the selected force files and reads one seed;
    /// it performs no integration to a requested fitting epoch.
    /// The seed's explicit SPICE target identifies its self-force, if present;
    /// nearby perturbers are retained regardless of the initial separation.
    pub fn new(
        session: Arc<AssistSession>,
        seed: &GeometricSpkSource,
        anchor: f64,
    ) -> Result<Self> {
        finite(anchor)?;
        let planets =
            session.config().planets_file.as_deref().ok_or_else(|| {
                invalid("N-body generation requires an explicit planet force file")
            })?;
        let asteroids = session
            .config()
            .asteroids_file
            .as_deref()
            .ok_or_else(|| invalid("N-body generation requires an asteroid force file"))?;
        let windows = {
            // Both explicit paths are still mandatory above. Use the actual
            // retained owners, not a second path acquisition which could name
            // a replaced file after the selected session was opened.
            let _ = (planets, asteroids);
            session.with_ephemeris(|ephemeris| Ok(ephemeris.inner.coverage()))?
        };
        let coverage = windows
            .into_iter()
            .find_map(|(low, high)| {
                let guarded = (low + INTEGRATION_GUARD_DAYS, high - INTEGRATION_GUARD_DAYS);
                (guarded.0 <= anchor && anchor <= guarded.1).then_some(guarded)
            })
            .ok_or_else(|| invalid("N-body anchor lacks guarded common force-file coverage"))?;
        let (position, velocity) = seed.heliocentric_seed_tdb(anchor)?;
        let heliocentric = CartesianState { position, velocity }.check()?;
        let initial = session.with_ephemeris(|ephemeris| {
            let time = anchor - ephemeris.inner.reference_jd();
            super::shifted(heliocentric, ephemeris.inner.state(0, time)?, true).check()
        })?;
        Ok(Self {
            session,
            anchor,
            initial,
            target: seed.target(),
            coverage,
        })
    }

    /// Inclusive guarded continuous TDB domain shared by every force target.
    /// The seed file is required at the anchor only; its narrower target span
    /// does not limit a genuinely integrated trajectory.
    pub fn coverage_tdb(&self) -> (f64, f64) {
        self.coverage
    }

    /// TDB Julian day of the fixed initial SPK sample.
    pub fn anchor_tdb(&self) -> f64 {
        self.anchor
    }

    /// Integrate SSB ICRS rows on a TT grid, AU and AU per TT day.
    ///
    /// The native clock converts finite TT Julian days to this integrator's
    /// scalar TDB epochs; its analytic rate scales the resulting velocities.
    /// Seed, force coverage, diagnostics, request ordering and limitations are
    /// unchanged from [`Self::states_tdb`]. Empty input performs no integration;
    /// any clock, domain or native failure returns no partial result.
    pub fn states_tt(&self, epochs: &[f64]) -> Result<PropagationOutcome<Vec<[f64; 6]>>> {
        let tdb: Result<Vec<_>> = epochs
            .iter()
            .map(|&tt| crate::time::tt_to_tdb_jd(tt))
            .collect();
        let mut outcome = self.states_tdb(&tdb?)?;
        for (row, &tt) in outcome.value.iter_mut().zip(epochs) {
            let rate = crate::time::tdb_rate(tt)?;
            for value in &mut row[3..] {
                *value *= rate;
            }
        }
        Ok(outcome)
    }

    /// Integrate geometric `[x,y,z,vx,vy,vz]` rows at finite TDB epochs.
    ///
    /// Each nonzero direction starts from the identical barycentric seed and
    /// proceeds monotonically with IAS15. Duplicate epochs reuse the same state;
    /// arbitrary input ordering is restored. An empty batch opens no files.
    /// Out-of-domain/nonfinite requests and native failures return an error
    /// without a partial result. Diagnostics retain self-force exclusions.
    pub fn states_tdb(&self, epochs: &[f64]) -> Result<PropagationOutcome<Vec<[f64; 6]>>> {
        for &epoch in epochs {
            finite(epoch)?;
            if epoch < self.coverage.0 || epoch > self.coverage.1 {
                return Err(invalid(
                    "N-body sample is outside guarded common force-file coverage",
                ));
            }
        }
        if epochs.is_empty() {
            return Ok(PropagationOutcome {
                value: Vec::new(),
                messages: Vec::new(),
                warnings: Vec::new(),
            });
        }
        {
            let mut order: Vec<_> = (0..epochs.len()).collect();
            order.sort_by(|&a, &b| epochs[a].total_cmp(&epochs[b]));
            let split = order.partition_point(|&index| epochs[index] < self.anchor);
            self.session.with_ephemeris(|ephemeris| {
                let mut result = PropagationOutcome {
                    value: vec![[0.0; 6]; epochs.len()],
                    messages: Vec::new(),
                    warnings: Vec::new(),
                };
                for sweep in [
                    order[split..].to_vec(),
                    order[..split].iter().copied().rev().collect(),
                ] {
                    if sweep.is_empty() {
                        continue;
                    }
                    let start = self.anchor - ephemeris.inner.reference_jd();
                    let (mut simulation, warnings) = super::assist_barycentric_initial(
                        self.initial,
                        start,
                        &ephemeris.inner,
                        None,
                        Some(self.target),
                    )?;
                    result.warnings.extend(warnings);
                    let mut last_epoch = self.anchor;
                    let mut state = self.initial;
                    for index in sweep {
                        let epoch = epochs[index];
                        if epoch != last_epoch {
                            result.messages.extend(
                                simulation.integrate(epoch - ephemeris.inner.reference_jd())?,
                            );
                            state = simulation.state(0)?.check()?;
                            last_epoch = epoch;
                        }
                        let [x, y, z] = state.position;
                        let [vx, vy, vz] = state.velocity;
                        result.value[index] = [x, y, z, vx, vy, vz];
                    }
                }
                Ok(result)
            })
        }
    }
}
