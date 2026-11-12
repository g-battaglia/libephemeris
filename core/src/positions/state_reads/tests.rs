//! The correction kernel can consume only its two scientific read operations.

use super::*;
use crate::positions::StencilCenter;
use crate::positions::correction_velocity::{ExpandedCenter, apparent_motion};
use crate::positions::derivative::AngularWord as Word;
use std::cell::Cell;

/// No StateProvider implementation or analytical-capability interface exists.
struct CorrectionOnly {
    scalar_reads: Cell<usize>,
    word_reads: Cell<usize>,
}

impl CorrectionStateReads for CorrectionOnly {
    fn geometric_state(&self, _: i32, _: f64) -> Result<GeometricState> {
        self.scalar_reads.set(self.scalar_reads.get() + 1);
        Err(crate::errors::Error::MathDomain { value: 13.0 })
    }

    fn geometric_state_word_parts_two_part(
        &self,
        _: i32,
        _: TwoPartEpoch,
    ) -> Result<CartesianWordParts> {
        self.word_reads.set(self.word_reads.get() + 1);
        Err(crate::errors::Error::MathDomain { value: 17.0 })
    }
}

#[test]
fn correction_kernel_accepts_two_reads_and_keeps_scalar_failure_before_words() {
    let source = CorrectionOnly {
        scalar_reads: Cell::new(0),
        word_reads: Cell::new(0),
    };
    let position = [1.0, 0.0, 0.0];
    let center = StencilCenter::new(
        position, [0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3], 2, true, false, None,
    )
    .unwrap();
    let result = apparent_motion(
        &source,
        &center,
        2451545.0,
        0.125,
        0.0,
        [0.0; 3],
        true,
        ExpandedCenter {
            position: position.map(Word::from),
            velocity: [Word::from(0.0); 3],
            observer: [Word::from(0.0); 3],
            observer_velocity: [Word::from(0.0); 3],
        },
        None,
    );
    assert!(
        matches!(result, Err(crate::errors::Error::MathDomain { value }) if value.to_bits() == 13.0_f64.to_bits())
    );
    assert_eq!(source.scalar_reads.get(), 1);
    assert_eq!(source.word_reads.get(), 0);
}
