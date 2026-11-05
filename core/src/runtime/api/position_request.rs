//! Original scalar position inputs, before any clock or flag preparation.

/// Time scale of the caller's Julian-day word.
///
/// Wrapping and unwrapping never validate, split or convert the value. Clock
/// selection remains with the calculation stage that actually consumes it.
#[derive(Clone, Copy, Debug)]
pub(super) enum InputEpoch {
    Tt(f64),
    Ut1(f64),
}

impl InputEpoch {
    /// Adapt an existing scalar/bool entry without changing its input bits.
    pub(super) fn from_entry(jd: f64, is_ut: bool) -> Self {
        if is_ut { Self::Ut1(jd) } else { Self::Tt(jd) }
    }

    /// Read the original Julian-day word without interpreting its scale.
    pub(super) fn scalar(self) -> f64 {
        match self {
            Self::Tt(jd) | Self::Ut1(jd) => jd,
        }
    }

    /// Adapt to numerical entry points that still accept a scalar/bool pair.
    pub(super) fn into_entry(self) -> (f64, bool) {
        (self.scalar(), matches!(self, Self::Ut1(_)))
    }
}

/// Unprepared operands for one position dispatch.
///
/// `input_flags` are the caller's raw flags, not normalized computation flags
/// or the flags returned by a completed calculation.
#[derive(Clone, Copy, Debug)]
pub(super) struct PositionRequest {
    pub(super) epoch: InputEpoch,
    pub(super) body: i32,
    pub(super) input_flags: i32,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn epoch_adapters_preserve_scale_and_every_scalar_bit() {
        let words = [
            0.0_f64.to_bits(),
            (-0.0_f64).to_bits(),
            2451545.125_f64.to_bits(),
            f64::INFINITY.to_bits(),
            f64::NEG_INFINITY.to_bits(),
            0x7ff8_0000_0000_1234, // Quiet NaN with a retained payload.
            0xfff0_0000_0000_0042, // Negative signaling NaN: no arithmetic.
        ];
        for word in words {
            for is_ut in [false, true] {
                let epoch = InputEpoch::from_entry(f64::from_bits(word), is_ut);
                assert_eq!(epoch.scalar().to_bits(), word);
                assert_eq!(matches!(epoch, InputEpoch::Ut1(_)), is_ut);
                let (jd, unwrapped_is_ut) = epoch.into_entry();
                assert_eq!(jd.to_bits(), word);
                assert_eq!(unwrapped_is_ut, is_ut);
            }
        }
    }

    #[test]
    fn request_retains_raw_flags_without_preparation() {
        let input_flags = i32::MIN | crate::constants::FLG_SPEED3 | crate::constants::FLG_SIDEREAL;
        let request = PositionRequest {
            epoch: InputEpoch::Ut1(-0.0),
            body: -crate::constants::TRUE_NODE,
            input_flags,
        };
        assert_eq!(request.input_flags, input_flags);
        assert_eq!(request.body, -crate::constants::TRUE_NODE);
        assert_eq!(request.epoch.scalar().to_bits(), (-0.0_f64).to_bits());
    }
}
