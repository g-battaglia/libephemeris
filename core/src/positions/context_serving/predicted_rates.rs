//! Call-scoped differentiation of the existing Auto scalar support recipe.
//!
//! LEB records use TT ownership; actual center-offset reads use their converted
//! split TDB words. The scalar adapter retains sampled center velocities,
//! including the old side probes, and deliberately advertises no new capability.

use super::*;
use crate::positions::routing::{PredictedSpkReads, model_serving};

/// Only the scalar methods used by `icrs_state` are adapted. Other public or
/// split-epoch provider operations never receive this private borrowed view.
struct ScalarOffsets<'a, P> {
    provider: &'a P,
    center_at: &'a dyn Fn(i32, f64) -> Result<Option<[f64; 3]>>,
    velocity_center_at: &'a dyn Fn(i32, f64) -> Result<Option<[f64; 3]>>,
    sample_offset_velocity: bool,
    stored_reads: Option<&'a PredictedSpkReads>,
}

impl<P: StateProvider> ScalarOffsets<'_, P> {
    fn stored_query<T>(&self, body: i32, tt: f64, read: impl FnOnce() -> Result<T>) -> Result<T> {
        self.stored_query_two_part(body, TwoPartEpoch { hi: tt, lo: 0.0 }, read)
    }

    fn stored_query_two_part<T>(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
        read: impl FnOnce() -> Result<T>,
    ) -> Result<T> {
        let Some(reads) = self.stored_reads else {
            return read();
        };
        reads.query(
            body,
            epoch,
            &|id, hi, lo| {
                let epoch = hi + lo;
                if !matches!(self.provider.serving_coord_type(id), Some(0 | 4)) {
                    return Err(Error::MathDomain { value: epoch });
                }
                let Some((start, end)) = self.provider.smooth_interval_tt(id, epoch)? else {
                    return Err(Error::MathDomain { value: epoch });
                };
                // Inward bounds retain the split tail in the ownership test.
                // Each subtraction is rounded toward the record's interior.
                Ok((
                    ((start - hi).next_up() - lo).next_up(),
                    ((end - hi).next_down() - lo).next_down(),
                ))
            },
            read,
        )
    }

    fn pctr_state(&self, body: i32, epoch: TwoPartEpoch) -> Result<positions::GeometricState> {
        if self
            .stored_reads
            .is_some_and(PredictedSpkReads::is_refining)
        {
            // A small emission-delay drift can be below a scalar Julian-date
            // ULP. Preserve it while differentiating the source polynomial;
            // central and legacy acquisitions retain their original scalar path.
            self.geometric_state_two_part(body, epoch)
        } else {
            self.geometric_state(body, epoch.value())
        }
    }
}

impl<P: StateProvider> StateProvider for ScalarOffsets<'_, P> {
    fn provider_name(&self) -> &'static str {
        self.provider.provider_name()
    }
    fn has_body(&self, id: i32) -> bool {
        self.provider.has_body(id)
    }
    fn body_coverage(&self, id: i32) -> Option<(f64, f64)> {
        self.provider.body_coverage(id)
    }
    fn serving_coord_type(&self, id: i32) -> Option<u32> {
        self.provider.serving_coord_type(id)
    }
    fn eval_raw(&self, id: i32, tt: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        self.stored_query(id, tt, || self.provider.eval_raw(id, tt))
    }
    fn eval_raw_two_part(&self, id: i32, tt: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
        self.stored_query_two_part(id, tt, || self.provider.eval_raw_two_part(id, tt))
    }
    fn eval_raw_inline(&self, id: i32, tt: f64) -> Result<Option<crate::leb::CartesianState>> {
        self.stored_query(id, tt, || self.provider.eval_raw_inline(id, tt))
    }
    fn eval_raw_inline_two_part(
        &self,
        id: i32,
        epoch: TwoPartEpoch,
    ) -> Result<Option<crate::leb::CartesianState>> {
        self.stored_query_two_part(id, epoch, || {
            self.provider.eval_raw_inline_two_part(id, epoch)
        })
    }
    fn serving_entry_snapshot(&self, id: i32) -> Result<EntrySnapshot> {
        self.provider.serving_entry_snapshot(id)
    }
    fn selected_tier(&self, id: i32, tt: f64) -> Option<LebTier> {
        self.provider.selected_tier(id, tt)
    }
    fn smooth_interval_tt(&self, id: i32, tt: f64) -> Result<Option<(f64, f64)>> {
        self.provider.smooth_interval_tt(id, tt)
    }
    fn body_center_offset_at(&self, body: i32, tt: f64) -> Result<Option<[f64; 3]>> {
        (self.center_at)(body, tt)
    }
    fn body_center_offsets(&self, body: i32, tt: f64) -> Result<[Option<[f64; 3]>; 3]> {
        let central = (self.center_at)(body, tt)?;
        if !self.sample_offset_velocity || central.is_none() || central == Some([0.0; 3]) {
            return Ok([None, central, None]);
        }
        let step = crate::leb::CENTER_OFFSET_VELOCITY_STEP_DAYS;
        Ok([
            (self.velocity_center_at)(body, tt - step)?,
            central,
            (self.velocity_center_at)(body, tt + step)?,
        ])
    }
}

/// PCTR's complete scalar LEB/center recipe with actual ordered read ownership.
/// Physical center velocities retain all original side probes; no analytic
/// provider capability is inferred. Tiered stored domains include priority edges.
pub(in crate::positions) fn pctr_state_with_owned_support(
    provider: &BoundProvider<'_>,
    reads: &PredictedSpkReads,
    body: i32,
    epoch: TwoPartEpoch,
) -> Result<positions::GeometricState> {
    let read: &OffsetRead<'_> = &|registry, target, whole, fraction| {
        reads.relative_query(
            target,
            TwoPartEpoch {
                hi: whole,
                lo: fraction,
            },
            &|target, whole, fraction| registry.relative_record_offsets(target, whole, fraction),
            || registry.eval_body_relative_two_part_jd(target, whole, fraction),
        )
    };
    let center_at = |body, epoch| provider.center_offset_using(body, epoch, read);
    ScalarOffsets {
        provider,
        center_at: &center_at,
        velocity_center_at: &center_at,
        sample_offset_velocity: true,
        stored_reads: Some(reads),
    }
    .pctr_state(body, epoch)
}

/// Session-only composition: caller-supplied model callbacks and all other
/// body/flag/provider cases retain the original generic entry and acquisitions.
#[allow(clippy::too_many_arguments)]
pub(crate) fn model_position_with_owned_support(
    models: Option<&positions::AnalyticalModelCallback<'_>>,
    provider: &BoundProvider<'_>,
    request_ut1: f64,
    tt: f64,
    body: i32,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    delta_t: &dyn Fn(f64) -> Result<f64>,
    ayanamsha_at: &dyn Fn(f64, i32) -> Result<f64>,
    sidereal_mode: i32,
) -> Result<(PositionResult, i32)> {
    use crate::constants::{FLG_BARYCTR, FLG_HELCTR, FLG_SPEED};
    let eligible = models.is_none()
        && matches!(provider, BoundProvider::WithCenters(_, _))
        && model_serving::is_builtin_predicted_rate_body(body)
        && flags & FLG_SPEED != 0
        && flags & (FLG_HELCTR | FLG_BARYCTR) == 0;
    if !eligible {
        return positions::serving_model_position_with_models(
            models,
            provider,
            request_ut1,
            tt,
            body,
            flags,
            site,
            delta_t,
            ayanamsha_at,
            sidereal_mode,
        );
    }
    let offsets = PredictedSpkReads::new(true);
    let read: &OffsetRead<'_> = &|registry, target, whole, fraction| {
        offsets.relative_query(
            target,
            TwoPartEpoch {
                hi: whole,
                lo: fraction,
            },
            &|target, whole, fraction| registry.relative_record_offsets(target, whole, fraction),
            || registry.eval_body_relative_two_part_jd(target, whole, fraction),
        )
    };
    let center_at = |body, epoch| provider.center_offset_using(body, epoch, read);
    let velocity_read: &OffsetRead<'_> = &|registry, target, whole, fraction| {
        offsets.central_velocity_probe(
            TwoPartEpoch {
                hi: whole,
                lo: fraction,
            },
            || registry.eval_body_relative_two_part_jd(target, whole, fraction),
        )
    };
    let velocity_center_at = |body, epoch| provider.center_offset_using(body, epoch, velocity_read);
    let scalar = ScalarOffsets {
        provider,
        center_at: &center_at,
        velocity_center_at: &center_at,
        sample_offset_velocity: true,
        stored_reads: None,
    };
    let original_deflector = ScalarOffsets {
        provider,
        center_at: &center_at,
        velocity_center_at: &velocity_center_at,
        sample_offset_velocity: true,
        stored_reads: None,
    };
    // This view is consumed only by the private PPN-position callback. It is
    // never returned as a full velocity or used for Earth's aberration state.
    let position = ScalarOffsets {
        provider,
        center_at: &center_at,
        velocity_center_at: &velocity_center_at,
        sample_offset_velocity: false,
        stored_reads: None,
    };
    model_serving::owned_provider_predicted_position(
        &scalar,
        &original_deflector,
        &position,
        &offsets,
        request_ut1,
        tt,
        body,
        flags,
        site,
        ayanamsha_at,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;

    struct Barycenter {
        reads: RefCell<Vec<f64>>,
    }

    #[test]
    fn pctr_owned_scalar_reads_keep_center_velocity_probes_and_stored_domain() {
        let provider = Barycenter {
            reads: RefCell::new(Vec::new()),
        };
        let reads = PredictedSpkReads::new(true);
        let probes = RefCell::new(Vec::new());
        let center = |_: i32, tt: f64| {
            let whole = tt.floor();
            let fraction = crate::time::tt_to_tdb_fraction(whole, tt - whole)?;
            reads.relative_query(
                599,
                TwoPartEpoch {
                    hi: whole,
                    lo: fraction,
                },
                &|_, _, _| Ok((-1.0, 1.0)),
                || {
                    probes.borrow_mut().push(tt);
                    Ok(Some([1.0 + tt * tt, 0.0, 0.0]))
                },
            )
        };
        let scalar = ScalarOffsets {
            provider: &provider,
            center_at: &center,
            velocity_center_at: &center,
            sample_offset_velocity: true,
            stored_reads: Some(&reads),
        };
        assert!(!scalar.has_analytic_kinematics());
        let place = |tt| {
            positions::cartesian_to_spherical_deg(
                scalar.geometric_state(crate::constants::JUPITER, tt)?.pos(),
            )
        };
        let central = place(0.5).unwrap();
        reads.finish_central();
        let rate = reads
            .rates_with_owned_leb(&place, 0.5, 0.125, false)
            .unwrap()
            .unwrap();
        let expected_longitude = (-1.0f64 / (6.25 * 6.25 + 1.0)).to_degrees();
        let expected_radial = 6.25 / (6.25f64 * 6.25 + 1.0).sqrt();
        assert!((rate[0] - expected_longitude).abs() < expected_longitude.abs() * 1e-8);
        assert!((rate[2] - expected_radial).abs() < expected_radial.abs() * 1e-8);
        assert_eq!(central[2], (6.25f64 * 6.25 + 1.0).sqrt());
        // Every state preserves all three actual physical-velocity probes.
        assert_eq!(probes.borrow().len(), provider.reads.borrow().len() * 3);
    }

    struct SplitLinear;

    impl StateProvider for SplitLinear {
        fn provider_name(&self) -> &'static str {
            "closed-split-emission"
        }
        fn has_body(&self, _: i32) -> bool {
            true
        }
        fn body_coverage(&self, _: i32) -> Option<(f64, f64)> {
            Some((2451544.0, 2451546.0))
        }
        fn smooth_interval_tt(&self, _: i32, _: f64) -> Result<Option<(f64, f64)>> {
            Ok(self.body_coverage(2))
        }
        fn serving_coord_type(&self, _: i32) -> Option<u32> {
            Some(0)
        }
        fn serving_entry_snapshot(&self, _: i32) -> Result<EntrySnapshot> {
            Ok(EntrySnapshot { coord_type: 0 })
        }
        fn eval_raw(&self, _: i32, tt: f64) -> Result<(Vec<f64>, Vec<f64>)> {
            Ok((
                vec![1.0 + (tt - 2451545.0) / 64.0, 0.0, 0.0],
                vec![1.0 / 64.0, 0.0, 0.0],
            ))
        }
        fn eval_raw_two_part(&self, _: i32, tt: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
            Ok((
                vec![1.0 + ((tt.hi - 2451545.0) + tt.lo) / 64.0, 0.0, 0.0],
                vec![1.0 / 64.0, 0.0, 0.0],
            ))
        }
    }

    #[test]
    fn pctr_owned_derivative_keeps_sub_ulp_emission_delay_drift() {
        // A source moving linearly at v, observed at t - (delay + slope*t),
        // has exact derivative v*(1-slope). Collapsing the emitted Julian
        // date can hide the delay slope over several successive dyadic scales.
        // Here even the outermost initial probe moves the delay by only one
        // quarter of the Julian-date ULP, so scalar rounding removes its drift.
        let evaluate = |legacy: bool| {
            let source = SplitLinear;
            let reads = PredictedSpkReads::pctr(true);
            let absent = |_, _| Ok(None);
            let scalar = ScalarOffsets {
                provider: &source,
                center_at: &absent,
                velocity_center_at: &absent,
                sample_offset_velocity: true,
                stored_reads: Some(&reads),
            };
            let sample = |at: f64| {
                let offset = at - 2451545.0;
                let epoch = TwoPartEpoch {
                    hi: 2451545.0,
                    lo: offset - (1.0 / 65536.0 + offset / 268435456.0),
                };
                let state = if legacy {
                    scalar.geometric_state(2, epoch.value())?
                } else {
                    scalar.pctr_state(2, epoch)?
                };
                Ok([0.0, 0.0, state.pos()[0]])
            };
            let central = sample(2451545.0).unwrap();
            assert_eq!(central[2], 1.0 - 1.0 / (65536.0 * 64.0));
            reads.finish_central();
            (
                reads
                    .rates_with_owned_leb(&sample, 2451545.0, 0.125, false)
                    .unwrap()
                    .unwrap()[2],
                central,
            )
        };
        let expected = (1.0 - 1.0 / 268435456.0) / 64.0;
        let (resolved, central) = evaluate(false);
        let (collapsed, old_central) = evaluate(true);
        assert_eq!(central, old_central);
        // Floating evaluation/subtraction has a conservative operation bound
        // at the algorithm's initial 1/64-day step, independent of PLAN limits.
        let arithmetic_bound = 128.0 * f64::EPSILON / (1.0 / 64.0);
        assert!((resolved - expected).abs() <= arithmetic_bound);
        assert!((collapsed - expected).abs() > arithmetic_bound);
    }

    impl StateProvider for Barycenter {
        fn provider_name(&self) -> &'static str {
            "closed-system-barycenter"
        }
        fn has_body(&self, _: i32) -> bool {
            true
        }
        fn body_coverage(&self, _: i32) -> Option<(f64, f64)> {
            Some((-1e10, 1e10))
        }
        fn smooth_interval_tt(&self, _: i32, _: f64) -> Result<Option<(f64, f64)>> {
            Ok(Some((-1e10, 1e10)))
        }
        fn serving_coord_type(&self, _: i32) -> Option<u32> {
            Some(crate::leb::COORD_ICRS_BARY_SYSTEM)
        }
        fn serving_entry_snapshot(&self, _: i32) -> Result<EntrySnapshot> {
            Ok(EntrySnapshot {
                coord_type: crate::leb::COORD_ICRS_BARY_SYSTEM,
            })
        }
        fn eval_raw(&self, body: i32, tt: f64) -> Result<(Vec<f64>, Vec<f64>)> {
            self.reads.borrow_mut().push(tt);
            let position = match body {
                crate::constants::SUN => vec![0.0; 3],
                crate::constants::EARTH => vec![1.0, 0.0, 0.0],
                _ => vec![5.0, 1.0, 0.0],
            };
            Ok((position, vec![0.0; 3]))
        }
        fn eval_raw_two_part(&self, _: i32, _: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
            panic!("the adapter must retain the scalar route")
        }
    }

    #[test]
    fn scalar_adapter_preserves_actual_center_probes_and_their_sampled_velocity() {
        let provider = Barycenter {
            reads: RefCell::new(Vec::new()),
        };
        let offsets = RefCell::new(Vec::new());
        let center = |_: i32, tt: f64| {
            offsets.borrow_mut().push(tt);
            Ok(Some([1.0 + tt * tt, 0.0, 0.0]))
        };
        let scalar = ScalarOffsets {
            provider: &provider,
            center_at: &center,
            velocity_center_at: &center,
            sample_offset_velocity: true,
            stored_reads: None,
        };
        assert!(!scalar.has_analytic_kinematics());
        let state = scalar
            .geometric_state(crate::constants::JUPITER, 0.5)
            .unwrap();
        let step = crate::leb::CENTER_OFFSET_VELOCITY_STEP_DAYS;
        assert_eq!(*offsets.borrow(), [0.5, 0.5 - step, 0.5 + step]);
        assert_eq!(*provider.reads.borrow(), [0.5]);
        assert_eq!(state.pos(), [6.25, 1.0, 0.0]);
        assert!((state.vel()[0] - 1.0).abs() < 1e-10);
        assert!(state.is_planet_center());
    }

    #[test]
    fn retarded_scalar_queries_and_velocity_probes_keep_their_own_tdb_words() {
        let provider = Barycenter {
            reads: RefCell::new(Vec::new()),
        };
        let epochs = RefCell::new(Vec::new());
        let reads = PredictedSpkReads::new(true);
        let center = |_: i32, tt: f64| {
            let whole = tt.floor();
            let fraction = crate::time::tt_to_tdb_fraction(whole, tt - whole)?;
            reads.relative_query(
                599,
                TwoPartEpoch {
                    hi: whole,
                    lo: fraction,
                },
                &|_, _, _| Ok((-1.0, 1.0)),
                || {
                    epochs.borrow_mut().push((tt, whole, fraction));
                    Ok(Some([1.0, 0.0, 0.0]))
                },
            )
        };
        let scalar = ScalarOffsets {
            provider: &provider,
            center_at: &center,
            velocity_center_at: &center,
            sample_offset_velocity: true,
            stored_reads: None,
        };
        let reception = 2451545.25;
        let emission = reception - 0.0125;
        for tt in [reception, emission] {
            scalar
                .geometric_state(crate::constants::JUPITER, tt)
                .unwrap();
        }
        let step = crate::leb::CENTER_OFFSET_VELOCITY_STEP_DAYS;
        let expected = [
            reception,
            reception - step,
            reception + step,
            emission,
            emission - step,
            emission + step,
        ];
        assert_eq!(*provider.reads.borrow(), [reception, emission]);
        for ((tt, whole, fraction), expected_tt) in epochs.borrow().iter().zip(expected) {
            assert_eq!(tt.to_bits(), expected_tt.to_bits());
            assert_eq!(whole.to_bits(), expected_tt.floor().to_bits());
            let expected_fraction = crate::time::tt_to_tdb_fraction(
                expected_tt.floor(),
                expected_tt - expected_tt.floor(),
            )
            .unwrap();
            assert_eq!(fraction.to_bits(), expected_fraction.to_bits());
            assert_ne!(
                fraction.to_bits(),
                (expected_tt - expected_tt.floor()).to_bits()
            );
        }
        assert_eq!(epochs.borrow().len(), 6);
    }

    #[test]
    fn scalar_adapter_keeps_missing_zero_and_error_center_laziness() {
        for value in [
            Ok(None),
            Ok(Some([0.0; 3])),
            Err(Error::MathDomain { value: 123.0 }),
        ] {
            let provider = Barycenter {
                reads: RefCell::new(Vec::new()),
            };
            let epochs = RefCell::new(Vec::new());
            let center = |_: i32, tt| {
                epochs.borrow_mut().push(tt);
                value.clone()
            };
            let scalar = ScalarOffsets {
                provider: &provider,
                center_at: &center,
                velocity_center_at: &center,
                sample_offset_velocity: true,
                stored_reads: None,
            };
            let state = scalar.geometric_state(crate::constants::JUPITER, 0.5);
            assert_eq!(*epochs.borrow(), [0.5]);
            if value.is_err() {
                assert!(matches!(state, Err(Error::MathDomain { value: 123.0 })));
                assert!(provider.reads.borrow().is_empty());
            } else {
                let state = state.unwrap();
                assert_eq!(state.pos(), [5.0, 1.0, 0.0]);
                assert_eq!(state.vel(), [0.0; 3]);
                assert!(!state.is_planet_center());
                assert_eq!(*provider.reads.borrow(), [0.5]);
            }
        }
    }

    #[test]
    fn owned_position_rates_do_not_inherit_unused_offset_velocity_seam() {
        let provider = Barycenter {
            reads: RefCell::new(Vec::new()),
        };
        let reads = PredictedSpkReads::new(true);
        let epochs = RefCell::new(Vec::new());
        let physical = |tt: f64| {
            epochs.borrow_mut().push(tt);
            Ok(Some([1.0 + 0.01 * tt + 0.5 * tt * tt, 0.02 * tt, 0.5]))
        };
        let center = |_: i32, tt: f64| {
            reads.relative_query(
                599,
                TwoPartEpoch { hi: 0.0, lo: tt },
                &|_, _, epoch| {
                    Ok(if epoch < 0.0 {
                        (-1.0 - epoch, -epoch)
                    } else {
                        (-epoch, 1.0 - epoch)
                    })
                },
                || physical(tt),
            )
        };
        let velocity = |_: i32, tt| {
            reads.central_velocity_probe(TwoPartEpoch { hi: 0.0, lo: tt }, || physical(tt))
        };
        let scalar = ScalarOffsets {
            provider: &provider,
            center_at: &center,
            velocity_center_at: &velocity,
            sample_offset_velocity: true,
            stored_reads: None,
        };
        let position_source = ScalarOffsets {
            provider: &provider,
            center_at: &center,
            velocity_center_at: &velocity,
            sample_offset_velocity: false,
            stored_reads: None,
        };
        let epoch = -2.0_f64.powi(-30);
        scalar
            .geometric_state(crate::constants::JUPITER, epoch)
            .unwrap();
        reads.finish_central();
        let place = |tt| {
            let position = position_source
                .geometric_state(crate::constants::JUPITER, tt)?
                .pos();
            Ok(model_serving::model_spherical(position))
        };
        let rates = reads
            .rates_with_owned_leb(&place, epoch, 0.1, true)
            .unwrap()
            .expect("the consumed position has a wide left-owned polynomial");
        let p = [
            6.0 + 0.01 * epoch + 0.5 * epoch * epoch,
            1.0 + 0.02 * epoch,
            0.5,
        ];
        let v = [0.01 + epoch, 0.02, 0.0];
        let rxy2 = p[0] * p[0] + p[1] * p[1];
        let radius = (rxy2 + p[2] * p[2]).sqrt();
        let radial = (p[0] * v[0] + p[1] * v[1]) / radius;
        let expected = [
            ((p[0] * v[1] - p[1] * v[0]) / rxy2).to_degrees(),
            ((-p[2] * radial / radius) / rxy2.sqrt()).to_degrees(),
            radial,
        ];
        for (component, (actual, expected)) in rates.into_iter().zip(expected).enumerate() {
            // PLAN §4: relative speed resolution, with the absolute angular
            // criterion only near zero. Radial units never use that floor.
            let limit = if component < 2 && expected.abs() <= 1e-6 {
                1e-6
            } else {
                expected.abs() * 1e-6
            };
            assert!((actual - expected).abs() <= limit, "{actual} vs {expected}");
        }
        assert!(epochs.borrow()[3..].iter().all(|&tt| tt < 0.0));
    }

    #[test]
    fn owned_full_state_keeps_nonzero_earth_offset_velocity_roles() {
        let provider = Barycenter {
            reads: RefCell::new(Vec::new()),
        };
        let reads = PredictedSpkReads::new(true);
        let position_roles = RefCell::new(Vec::new());
        let physical = |tt: f64| {
            let dt = tt - 2451545.0;
            Ok(Some([0.001 + 0.0001 * dt * dt, 0.0001 * dt, 0.0]))
        };
        let center = |body, tt: f64| {
            if body != crate::constants::EARTH {
                return Ok(None);
            }
            let dt = tt - 2451545.0;
            reads.relative_query(
                399,
                TwoPartEpoch { hi: 0.0, lo: dt },
                &|_, _, epoch| Ok((-1.0 - epoch, 1.0 - epoch)),
                || physical(tt),
            )
        };
        let auxiliary = |body, tt: f64| {
            if body != crate::constants::EARTH {
                return Ok(None);
            }
            reads.central_velocity_probe(TwoPartEpoch { hi: tt, lo: 0.0 }, || physical(tt))
        };
        let position_at = |body, tt| {
            position_roles.borrow_mut().push(body);
            center(body, tt)
        };
        let full = ScalarOffsets {
            provider: &provider,
            center_at: &center,
            velocity_center_at: &center,
            sample_offset_velocity: true,
            stored_reads: None,
        };
        let original_deflector = ScalarOffsets {
            provider: &provider,
            center_at: &center,
            velocity_center_at: &auxiliary,
            sample_offset_velocity: true,
            stored_reads: None,
        };
        let position = ScalarOffsets {
            provider: &provider,
            center_at: &position_at,
            velocity_center_at: &auxiliary,
            sample_offset_velocity: false,
            stored_reads: None,
        };
        let state = model_serving::owned_provider_predicted_position(
            &full,
            &original_deflector,
            &position,
            &reads,
            2451545.0,
            2451545.0,
            crate::constants::NEPTUNE_ADAMS,
            crate::constants::FLG_SPEED | crate::constants::FLG_J2000,
            None,
            &|_, _| panic!("unused ayanamsha"),
        )
        .unwrap()
        .0;
        assert!(state.iter().all(|value| value.is_finite()));
        // Reaching the refined PPN callback proves that the full Earth
        // state, including both nonzero offset-velocity probes, was admitted.
        assert!(!position_roles.borrow().is_empty());
        assert!(
            position_roles
                .borrow()
                .iter()
                .all(|&body| body != crate::constants::EARTH)
        );
    }
}
