//! Narrow-source tests: analytical motion requires no public provider policy.

use super::*;
use crate::positions::state_reads::{ApparentRateReads, CorrectionStateReads};
use crate::positions::{CartesianWordParts, GeometricState, StateFrame};
use std::cell::RefCell;

#[derive(Debug, PartialEq, Eq)]
enum Read {
    Scalar(i32, u64),
    Word(i32, u64, u64),
    Split(i32, u64, u64),
    Relative(i32, u64, u64, i32, u64, u64),
}

#[derive(Clone, Copy)]
enum PairedVelocity {
    Available,
    Unavailable,
    Failed,
}

/// This source intentionally does not implement StateProvider. In particular,
/// analytical eligibility and acceleration acquisition cannot be invoked here.
struct RateSource {
    reads: RefCell<Vec<Read>>,
    paired: PairedVelocity,
    fail_word_body: Option<i32>,
}

impl RateSource {
    fn new(paired: PairedVelocity) -> Self {
        Self {
            reads: RefCell::new(Vec::new()),
            paired,
            fail_word_body: None,
        }
    }

    fn state(body: i32, epoch: TwoPartEpoch) -> GeometricState {
        GeometricState {
            body_id: body,
            jd_tt: epoch.value(),
            pos: if body == 14 {
                [0.0; 3]
            } else {
                [2.0, 0.5, 0.25]
            },
            vel: [0.0; 3],
            frame: StateFrame::IcrsBarycentric,
            planet_center: false,
            tier: None,
        }
    }
}

impl CorrectionStateReads for RateSource {
    fn geometric_state(&self, body: i32, jd_tt: f64) -> Result<GeometricState> {
        self.reads
            .borrow_mut()
            .push(Read::Scalar(body, jd_tt.to_bits()));
        Ok(Self::state(body, TwoPartEpoch::split(jd_tt)?))
    }

    fn geometric_state_word_parts_two_part(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<CartesianWordParts> {
        self.reads
            .borrow_mut()
            .push(Read::Word(body, epoch.hi.to_bits(), epoch.lo.to_bits()));
        if self.fail_word_body == Some(body) {
            return Err(crate::errors::Error::MathDomain { value: 79.0 });
        }
        let state = Self::state(body, epoch);
        Ok((
            state.pos().map(|value| [value, 0.0]),
            state.vel().map(|value| [value, 0.0]),
        ))
    }
}

impl ApparentRateReads for RateSource {
    fn geometric_state_two_part(&self, body: i32, epoch: TwoPartEpoch) -> Result<GeometricState> {
        self.reads
            .borrow_mut()
            .push(Read::Split(body, epoch.hi.to_bits(), epoch.lo.to_bits()));
        Ok(Self::state(body, epoch))
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
        match self.paired {
            PairedVelocity::Available => Ok(Some([0.0; 3])),
            PairedVelocity::Unavailable => Ok(None),
            PairedVelocity::Failed => Err(crate::errors::Error::MathDomain { value: 61.0 }),
        }
    }
}

fn motion(
    source: &RateSource,
    epoch: TwoPartEpoch,
    observer: Option<i32>,
) -> Result<OwnedApparentMotion> {
    let position = [2.0, 0.5, 0.25];
    let light_time = LightTimeSolution {
        geo_pos: position,
        geo_vel: [0.0; 3],
        retarded_pos: position,
        retarded_vel: [0.0; 3],
        light_time_days: 0.125,
        light_time_rate: 0.0,
        retarded_jd_tt: epoch.value() - 0.125,
    };
    let stencil = StencilCenter::new(
        position, [0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3], 2, false, false, None,
    )?;
    analytic_apparent_motion(
        source,
        &RateObservation {
            jd_tt: epoch.value(),
            obs_split: Some((epoch.hi, epoch.lo)),
            observer_body: observer,
            light_time: &light_time,
            stencil: &stencil,
        },
        [0.0; 3],
        None,
    )
}

fn expected_reads(epoch: TwoPartEpoch, observer: Option<i32>) -> Vec<Read> {
    let emission = epoch.add_offset(-0.125).unwrap();
    let mut reads = Vec::new();
    if let Some(body) = observer {
        reads.push(Read::Relative(
            2,
            emission.hi.to_bits(),
            emission.lo.to_bits(),
            body,
            epoch.hi.to_bits(),
            epoch.lo.to_bits(),
        ));
    }
    reads.push(Read::Word(2, emission.hi.to_bits(), emission.lo.to_bits()));
    if let Some(body) = observer {
        reads.push(Read::Word(body, epoch.hi.to_bits(), epoch.lo.to_bits()));
        reads.push(Read::Split(body, epoch.hi.to_bits(), epoch.lo.to_bits()));
    }
    reads
}

#[test]
fn narrow_source_preserves_read_order_split_bits_and_absent_observer() {
    for epoch in [
        TwoPartEpoch {
            hi: 2451545.0,
            lo: 2.0_f64.powi(-40),
        },
        TwoPartEpoch { hi: -0.0, lo: -0.0 },
    ] {
        for observer in [Some(14), None] {
            for paired in [PairedVelocity::Available, PairedVelocity::Unavailable] {
                let source = RateSource::new(paired);
                let result = motion(&source, epoch, observer).unwrap();
                assert_eq!(*source.reads.borrow(), expected_reads(epoch, observer));
                assert_eq!(result.0.map(Word::rounded), [2.0, 0.5, 0.25]);
                assert_eq!(result.1.map(Word::rounded), [0.0; 3]);
            }
        }
    }
}

#[test]
fn relative_failure_precedes_every_word_read() {
    let source = RateSource::new(PairedVelocity::Failed);
    let epoch = TwoPartEpoch {
        hi: 2451545.0,
        lo: 2.0_f64.powi(-40),
    };
    assert!(
        matches!(motion(&source, epoch, Some(14)), Err(crate::errors::Error::MathDomain { value }) if value.to_bits() == 61.0_f64.to_bits())
    );
    assert_eq!(*source.reads.borrow(), expected_reads(epoch, Some(14))[..1]);
}

#[test]
fn word_failure_preserves_the_completed_read_prefix() {
    let epoch = TwoPartEpoch {
        hi: 2451545.0,
        lo: 2.0_f64.powi(-40),
    };
    for (body, completed) in [(2, 2), (14, 3)] {
        let mut source = RateSource::new(PairedVelocity::Unavailable);
        source.fail_word_body = Some(body);
        assert!(
            matches!(motion(&source, epoch, Some(14)), Err(crate::errors::Error::MathDomain { value }) if value.to_bits() == 79.0_f64.to_bits())
        );
        assert_eq!(
            *source.reads.borrow(),
            expected_reads(epoch, Some(14))[..completed]
        );
    }
}
