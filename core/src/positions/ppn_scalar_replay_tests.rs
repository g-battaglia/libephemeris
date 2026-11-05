//! Same-invocation scalar PPN reuse; the physical and Word formulas stay shared.

use super::*;
use crate::positions::correction_velocity::{ExpandedCenter, apparent_motion};
use crate::positions::derivative::AngularWord as Word;
use std::cell::{Cell, RefCell};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Read {
    Offsets(i32, u64),
    Entry(i32),
    Scalar(i32, u64),
    Split(i32, u64, u64),
    Word(i32, u64, u64),
    Tier(i32, u64),
    Domain(i32, u64),
    Relative(i32, u64, u64, i32, u64, u64),
    Acceleration(i32, u64),
}

struct Model {
    reads: RefCell<Vec<Read>>,
    all_deflectors: bool,
    fail_second_reception: bool,
    fail_word: Option<(i32, u64)>,
    fail_domains: Cell<bool>,
    fail_relative: bool,
    fail_center: bool,
}

impl Model {
    fn new() -> Self {
        Self {
            reads: RefCell::new(Vec::new()),
            all_deflectors: true,
            fail_second_reception: false,
            fail_word: None,
            fail_domains: Cell::new(false),
            fail_relative: false,
            fail_center: false,
        }
    }

    fn serves(&self, body: i32) -> bool {
        matches!(body, EARTH | MERCURY | SUN)
            || (self.all_deflectors && DEFLECTORS.iter().any(|entry| entry.0 == body))
    }

    fn values(&self, body: i32, tt: f64) -> Result<([f64; 3], [f64; 3])> {
        if !self.serves(body) {
            return Err(Error::LebBodyNotFound { body_id: body });
        }
        let (origin, velocity) = match body {
            EARTH => ([0.1, -0.2, 0.3], [0.002, -0.01, 0.003]),
            MERCURY => ([1.0, 0.2, 0.3], [0.02, -0.03, 0.01]),
            SUN => ([0.3, 1.0, 0.2], [0.001, 0.002, -0.003]),
            body if body == DEFLECTORS[1].0 => ([-2.0, 3.0, 0.5], [0.001, -0.002, 0.001]),
            _ => ([4.0, -1.0, 2.0], [-0.001, 0.001, 0.002]),
        };
        let offset = tt - FRAME_J2000;
        Ok((
            std::array::from_fn(|i| origin[i] + velocity[i] * offset),
            velocity,
        ))
    }

    fn scalars(&self) -> Vec<(i32, u64)> {
        self.reads
            .borrow()
            .iter()
            .filter_map(|read| match *read {
                Read::Scalar(body, tt) => Some((body, tt)),
                _ => None,
            })
            .collect()
    }

    fn words(&self) -> Vec<Read> {
        self.reads
            .borrow()
            .iter()
            .filter(|read| matches!(read, Read::Word(..)))
            .cloned()
            .collect()
    }
}

impl StateProvider for Model {
    fn provider_name(&self) -> &'static str {
        "test/affine-physical-ppn"
    }
    fn has_body(&self, body: i32) -> bool {
        self.serves(body)
    }
    fn body_coverage(&self, body: i32) -> Option<(f64, f64)> {
        self.serves(body)
            .then_some((FRAME_J2000 - 1.0, FRAME_J2000 + 1.0))
    }
    fn serving_coord_type(&self, body: i32) -> Option<u32> {
        self.serves(body).then_some(crate::leb::COORD_ICRS_BARY)
    }
    fn serving_entry_snapshot(&self, body: i32) -> Result<EntrySnapshot> {
        self.reads.borrow_mut().push(Read::Entry(body));
        if !self.serves(body) {
            return Err(Error::LebBodyNotFound { body_id: body });
        }
        Ok(EntrySnapshot {
            coord_type: crate::leb::COORD_ICRS_BARY,
        })
    }
    fn selected_tier(&self, body: i32, tt: f64) -> Option<LebTier> {
        self.reads.borrow_mut().push(Read::Tier(body, tt.to_bits()));
        None
    }
    fn body_center_offsets(&self, body: i32, tt: f64) -> Result<[Option<[f64; 3]>; 3]> {
        self.reads
            .borrow_mut()
            .push(Read::Offsets(body, tt.to_bits()));
        if self.fail_center && body == SUN {
            return Err(Error::MathDomain { value: 53.0 });
        }
        Ok([None; 3])
    }
    fn eval_raw(&self, body: i32, tt: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        self.reads
            .borrow_mut()
            .push(Read::Scalar(body, tt.to_bits()));
        if self.fail_second_reception
            && body == SUN
            && tt.to_bits() == FRAME_J2000.to_bits()
            && self
                .scalars()
                .iter()
                .filter(|key| **key == (body, tt.to_bits()))
                .count()
                == 2
        {
            return Err(Error::MathDomain { value: 73.0 });
        }
        let (position, velocity) = self.values(body, tt)?;
        Ok((position.to_vec(), velocity.to_vec()))
    }
    fn eval_raw_two_part(&self, body: i32, epoch: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
        self.reads
            .borrow_mut()
            .push(Read::Split(body, epoch.hi.to_bits(), epoch.lo.to_bits()));
        let (position, velocity) = self.values(body, epoch.value())?;
        Ok((position.to_vec(), velocity.to_vec()))
    }
    fn geometric_state_word_parts_two_part(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<CartesianWordParts> {
        self.reads
            .borrow_mut()
            .push(Read::Word(body, epoch.hi.to_bits(), epoch.lo.to_bits()));
        if self.fail_word == Some((body, epoch.value().to_bits())) {
            return Err(Error::MathDomain { value: 79.0 });
        }
        // This affine test source's public words are its complete supplied
        // model. This test does not replace native polynomial precision gates.
        let (position, velocity) = self.values(body, epoch.value())?;
        Ok((
            position.map(|value| [value, 0.0]),
            velocity.map(|value| [value, 0.0]),
        ))
    }
    fn has_analytic_kinematics(&self) -> bool {
        true
    }
    fn geometric_acceleration_tt(&self, body: i32, tt: f64) -> Result<Option<[f64; 3]>> {
        self.reads
            .borrow_mut()
            .push(Read::Acceleration(body, tt.to_bits()));
        Ok((body == EARTH).then_some([0.0; 3]))
    }
    fn nutation_rate_tt(&self, _: f64) -> Result<Option<(f64, f64)>> {
        Ok(Some((0.0, 0.0)))
    }
    fn eval_relative_velocity_two_part(
        &self,
        target: i32,
        target_epoch: TwoPartEpoch,
        center: i32,
        center_epoch: TwoPartEpoch,
    ) -> Result<Option<[f64; 3]>> {
        self.reads.borrow_mut().push(Read::Relative(
            target,
            target_epoch.hi.to_bits(),
            target_epoch.lo.to_bits(),
            center,
            center_epoch.hi.to_bits(),
            center_epoch.lo.to_bits(),
        ));
        if self.fail_relative {
            return Err(Error::MathDomain { value: 61.0 });
        }
        Ok(None)
    }
    fn smooth_interval_tt(&self, body: i32, tt: f64) -> Result<Option<(f64, f64)>> {
        self.reads
            .borrow_mut()
            .push(Read::Domain(body, tt.to_bits()));
        if self.fail_domains.get() {
            return Err(Error::MathDomain { value: 31.0 });
        }
        Ok(None)
    }
}

fn geometry() -> (StencilCenter, ExpandedCenter, f64, f64) {
    let ray = [0.9, 0.4, 0.0];
    let velocity = [0.018, -0.02, 0.007];
    let observer = [0.1, -0.2, 0.3];
    let observer_velocity = [0.002, -0.01, 0.003];
    let delay = vec3_dist(ray) / C_LIGHT_AU_DAY;
    let delay_rate =
        super::super::planet_range::RangeState::from_cartesian(ray, velocity).rate / C_LIGHT_AU_DAY;
    let center = StencilCenter::new(
        ray,
        velocity,
        observer,
        observer_velocity,
        [0.0; 3],
        MERCURY,
        true,
        true,
        None,
    )
    .unwrap();
    let expanded = ExpandedCenter {
        position: ray.map(Word::from),
        velocity: velocity.map(Word::from),
        observer: observer.map(Word::from),
        observer_velocity: observer_velocity.map(Word::from),
    };
    (center, expanded, delay, delay_rate)
}

fn motion(
    model: &Model,
    samples: Option<&PpnScalarSamples>,
    required: bool,
) -> Result<([Word; 3], [Word; 3])> {
    let (center, expanded, delay, rate) = geometry();
    apparent_motion(
        &super::state_reads::ProviderReads::new(model),
        &center,
        FRAME_J2000,
        delay,
        rate,
        [0.0; 3],
        required,
        expanded,
        samples,
    )
}

fn bits(value: ([Word; 3], [Word; 3])) -> ([[u64; 2]; 3], [[u64; 2]; 3]) {
    let convert = |words: [Word; 3]| words.map(|word| [word.hi.to_bits(), word.lo.to_bits()]);
    (convert(value.0), convert(value.1))
}

fn capture_scalar(
    model: &Model,
    samples: Option<&PpnScalarSamples>,
) -> Result<[Option<CorrectionRead>; DEFLECTORS.len() * 2]> {
    let reads = Cell::new([None; DEFLECTORS.len() * 2]);
    let (center, _, delay, _) = geometry();
    crate::positions::frames::apply_gravitational_deflection_with_correction(
        center.geo_geom,
        center.observer,
        FRAME_J2000,
        delay,
        |body, tt| capture_correction_read(model, body, tt, Some(&reads), samples),
    )?;
    Ok(reads.get())
}

#[test]
fn borrowed_provider_view_is_inert_and_forwards_original_split_words() {
    use super::state_reads::{CorrectionStateReads, ProviderReads};
    let model = Model::new();
    let source = ProviderReads::new(&model);
    assert!(model.reads.borrow().is_empty());
    let epoch = TwoPartEpoch {
        hi: FRAME_J2000,
        lo: -0.0,
    };
    let actual =
        CorrectionStateReads::geometric_state_word_parts_two_part(&source, MERCURY, epoch).unwrap();
    let (position, velocity) = model.values(MERCURY, epoch.value()).unwrap();
    let expected = (
        position.map(|value| [value, 0.0]),
        velocity.map(|value| [value, 0.0]),
    );
    let bits = |pair: CartesianWordParts| {
        (
            pair.0.map(|word| word.map(f64::to_bits)),
            pair.1.map(|word| word.map(f64::to_bits)),
        )
    };
    assert_eq!(bits(actual), bits(expected));
    assert_eq!(
        *model.reads.borrow(),
        vec![Read::Word(MERCURY, epoch.hi.to_bits(), epoch.lo.to_bits())]
    );
}

#[test]
fn successful_ppn_replay_removes_only_identical_scalar_queries_and_keeps_word_bits() {
    for required in [false, true] {
        let baseline = Model::new();
        let candidate = Model::new();
        let samples = PpnScalarSamples::default();
        let original_reads = capture_scalar(&baseline, None).unwrap();
        let captured_reads = capture_scalar(&candidate, Some(&samples)).unwrap();
        let keys = |reads: [Option<CorrectionRead>; DEFLECTORS.len() * 2]| {
            reads
                .into_iter()
                .flatten()
                .map(|read| (read.body, read.epoch.to_bits()))
                .collect::<Vec<_>>()
        };
        assert_eq!(keys(original_reads), keys(captured_reads));
        baseline.reads.borrow_mut().clear();
        candidate.reads.borrow_mut().clear();
        let expected = motion(&baseline, None, required).unwrap();
        let actual = motion(&candidate, Some(&samples), required).unwrap();
        assert_eq!(bits(actual), bits(expected));
        assert_eq!(candidate.words(), baseline.words());
        let original = baseline.scalars();
        let captured_keys = keys(captured_reads);
        let unmatched = original
            .iter()
            .copied()
            .filter(|key| !captured_keys.contains(key))
            .collect::<Vec<_>>();
        assert_eq!(candidate.scalars(), unmatched);
        assert!(original.len() - unmatched.len() >= 3);
    }
}

#[test]
fn actual_analytic_closest_epoch_must_be_acquired_when_only_its_neighbor_was_captured() {
    let baseline = Model::new();
    let expected = motion(&baseline, None, false).unwrap();
    let closest = baseline
        .scalars()
        .into_iter()
        .find(|key| key.0 == SUN && key.1 != FRAME_J2000.to_bits())
        .unwrap();
    let candidate = Model::new();
    let samples = PpnScalarSamples::default();
    let reads = Cell::new([None; DEFLECTORS.len() * 2]);
    capture_correction_read(&candidate, SUN, FRAME_J2000, Some(&reads), Some(&samples)).unwrap();
    let neighbor = f64::from_bits(closest.1).next_up();
    assert_ne!(neighbor.to_bits(), closest.1);
    capture_correction_read(&candidate, SUN, neighbor, Some(&reads), Some(&samples)).unwrap();
    candidate.reads.borrow_mut().clear();
    let actual = motion(&candidate, Some(&samples), false).unwrap();
    assert_eq!(bits(actual), bits(expected));
    assert!(candidate.scalars().contains(&closest));
    assert!(!candidate.scalars().contains(&(SUN, FRAME_J2000.to_bits())));
    assert_eq!(candidate.words(), baseline.words());
}

#[test]
fn scalar_sample_keys_keep_body_signed_zero_and_state_identity_exact() {
    let model = Model::new();
    let samples = PpnScalarSamples::default();
    let reads = Cell::new([None; DEFLECTORS.len() * 2]);
    capture_correction_read(&model, SUN, 0.0, Some(&reads), Some(&samples)).unwrap();
    assert!(samples.state_at(SUN, 0.0).is_some());
    assert!(samples.state_at(SUN, -0.0).is_none());
    assert!(samples.state_at(EARTH, 0.0).is_none());
    assert!(samples.state_at(SUN, 0.0_f64.next_up()).is_none());
    let wrong_state = model.geometric_state(EARTH, FRAME_J2000).unwrap();
    samples.record(1, SUN, FRAME_J2000, wrong_state);
    assert!(samples.state_at(SUN, FRAME_J2000).is_none());
}

#[test]
fn scalar_capture_does_not_advance_or_hide_domain_validation() {
    let model = Model::new();
    let samples = PpnScalarSamples::default();
    let reads = capture_scalar(&model, Some(&samples)).unwrap();
    assert!(
        !model
            .reads
            .borrow()
            .iter()
            .any(|read| matches!(read, Read::Domain(..)))
    );
    model.fail_domains.set(true);
    assert!(matches!(correction_read_domains(&model, reads),
        Err(Error::MathDomain { value }) if value.to_bits() == 31.0_f64.to_bits()));
    let first = reads.into_iter().flatten().next().unwrap();
    assert_eq!(
        model.reads.borrow().last(),
        Some(&Read::Domain(first.body, first.epoch.to_bits()))
    );
}

#[test]
fn optional_missing_states_are_not_recorded_or_reclassified() {
    let mut baseline = Model::new();
    baseline.all_deflectors = false;
    let mut candidate = Model::new();
    candidate.all_deflectors = false;
    let samples = PpnScalarSamples::default();
    capture_scalar(&candidate, Some(&samples)).unwrap();
    candidate.reads.borrow_mut().clear();
    assert_eq!(
        bits(motion(&candidate, Some(&samples), false).unwrap()),
        bits(motion(&baseline, None, false).unwrap())
    );
    for (body, _) in DEFLECTORS.into_iter().skip(1) {
        assert!(samples.state_at(body, FRAME_J2000).is_none());
        assert!(candidate.reads.borrow().contains(&Read::Entry(body)));
        assert!(
            !candidate
                .words()
                .iter()
                .any(|read| matches!(read, Read::Word(id, ..) if *id == body))
        );
    }
    assert!(matches!(motion(&candidate, Some(&samples), true),
        Err(Error::LebBodyNotFound { body_id }) if body_id == DEFLECTORS[1].0));
}

#[test]
fn full_native_serving_preserves_words_flags_trace_and_relative_acquisition() {
    for equatorial in [0, FLG_EQUATORIAL] {
        for xyz in [false, true] {
            let baseline_model = Model::new();
            let candidate_model = Model::new();
            let nutation = |_| Ok((0.0, 0.0));
            let mut baseline = SealedServing::new_stored(&baseline_model, &nutation);
            // The test source has no scalar/split reuse capability and takes
            // the analytic branch, so other uses of this existing flag do not
            // remove reads. Only the new scalar PPN replay is under comparison.
            baseline.reuse_immutable_central = false;
            let candidate = SealedServing::new_stored(&candidate_model, &nutation);
            let flags = FLG_SPEED | FLG_SWIEPH | FLG_J2000 | equatorial;
            let (expected, old_trace) = crate::trace::capture_events(|| {
                serving_position_core(&baseline, FRAME_J2000, MERCURY, flags, xyz, None).unwrap()
            });
            let (actual, new_trace) = crate::trace::capture_events(|| {
                serving_position_core(&candidate, FRAME_J2000, MERCURY, flags, xyz, None).unwrap()
            });
            assert_eq!(actual.0.map(f64::to_bits), expected.0.map(f64::to_bits));
            assert_eq!(actual.1, expected.1);
            assert_eq!(new_trace, old_trace);
            assert_eq!(candidate_model.words(), baseline_model.words());
            let relatives = |model: &Model| {
                model
                    .reads
                    .borrow()
                    .iter()
                    .filter(|read| matches!(read, Read::Relative(..)))
                    .cloned()
                    .collect::<Vec<_>>()
            };
            assert_eq!(relatives(&candidate_model), relatives(&baseline_model));
            assert_eq!(relatives(&candidate_model).len(), 1);
            let ppn_count = |model: &Model| {
                model
                    .scalars()
                    .iter()
                    .filter(|key| DEFLECTORS.iter().any(|entry| entry.0 == key.0))
                    .count()
            };
            assert!(ppn_count(&candidate_model) < ppn_count(&baseline_model));
        }
    }
}

#[test]
fn public_opaque_serving_keeps_second_reception_callback_failure() {
    let mut model = Model::new();
    model.fail_second_reception = true;
    let nutation = |_| Ok((0.0, 0.0));
    let serving = SealedServing::new(&model, &nutation);
    let error = serving_position_core(
        &serving,
        FRAME_J2000,
        MERCURY,
        FLG_SPEED | FLG_SWIEPH | FLG_J2000,
        false,
        None,
    )
    .unwrap_err();
    assert!(matches!(error, Error::MathDomain { value } if value.to_bits() == 73.0_f64.to_bits()));
    assert_eq!(
        model
            .scalars()
            .iter()
            .filter(|key| **key == (SUN, FRAME_J2000.to_bits()))
            .count(),
        2
    );
    assert!(
        !model
            .words()
            .iter()
            .any(|read| matches!(read, Read::Word(SUN, ..)))
    );
}

#[test]
fn reused_scalar_state_does_not_skip_later_word_failure() {
    let mut model = Model::new();
    model.fail_word = Some((SUN, FRAME_J2000.to_bits()));
    let nutation = |_| Ok((0.0, 0.0));
    let serving = SealedServing::new_stored(&model, &nutation);
    let error = serving_position_core(
        &serving,
        FRAME_J2000,
        MERCURY,
        FLG_SPEED | FLG_SWIEPH | FLG_J2000,
        false,
        None,
    )
    .unwrap_err();
    assert!(matches!(error, Error::MathDomain { value } if value.to_bits() == 79.0_f64.to_bits()));
    assert_eq!(
        model
            .scalars()
            .iter()
            .filter(|key| **key == (SUN, FRAME_J2000.to_bits()))
            .count(),
        1
    );
    assert!(
        model
            .words()
            .iter()
            .any(|read| matches!(read, Read::Word(SUN, ..)))
    );
}

#[test]
fn center_and_relative_errors_still_precede_analytic_ppn() {
    for fail_center in [false, true] {
        let mut model = Model::new();
        model.fail_center = fail_center;
        model.fail_relative = !fail_center;
        let nutation = |_| Ok((0.0, 0.0));
        let serving = SealedServing::new_stored(&model, &nutation);
        let (error, events) = crate::trace::capture_events(|| {
            serving_position_core(
                &serving,
                FRAME_J2000,
                MERCURY,
                FLG_SPEED | FLG_SWIEPH | FLG_J2000,
                false,
                None,
            )
            .unwrap_err()
        });
        let wanted = if fail_center { 53.0_f64 } else { 61.0_f64 };
        assert!(
            matches!(error, Error::MathDomain { value } if value.to_bits() == wanted.to_bits())
        );
        assert!(model.words().is_empty());
        assert!(events.is_empty());
    }
}

#[test]
fn force_sampled_serving_does_not_enter_word_or_replay_motion() {
    let model = Model::new();
    let nutation = |_| Ok((0.0, 0.0));
    let base = SealedServing::new_stored(&model, &nutation);
    let sampled = base.with_sampled_derivatives();
    serving_position_core(
        &sampled,
        FRAME_J2000,
        MERCURY,
        FLG_SPEED | FLG_SWIEPH | FLG_J2000,
        false,
        None,
    )
    .unwrap();
    assert!(model.words().is_empty());
    assert!(
        !model
            .reads
            .borrow()
            .iter()
            .any(|read| matches!(read, Read::Acceleration(..)))
    );
}
