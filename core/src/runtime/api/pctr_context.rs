//! Planet-centric reduction over an acquired split UT1/TT epoch.
//!
//! The context supplies observer and sidereal settings; the calculation
//! retains the input time parts and immutable resource owners.

use super::*;

impl Session {
    /// Nonsidereal planet-centric position at an acquired split epoch.
    ///
    /// `whole` is the UT1 whole part, `tt_fraction` its TT offset and
    /// `original_ut1_fraction` the retained UT1 fraction. Epochs are Julian Days.
    /// `target` and `center` identify standard planets. The supplied flags govern
    /// apparent corrections and output representation.
    ///
    /// Results are longitude/latitude in degrees, range in AU and daily rates,
    /// or the requested Cartesian/radian form, plus flags and advisories.
    /// The caller owns file acquisition and request validation; this method
    /// retains the split epoch and propagates typed provider/range errors.
    #[allow(clippy::too_many_arguments)]
    pub fn calc_pctr_context_ut1_parts(
        &self,
        whole: f64,
        tt_fraction: f64,
        original_ut1_fraction: f64,
        target: i32,
        center: i32,
        calc_flags: i32,
        normalized_flags: i32,
    ) -> Result<Position, Error> {
        let (values, returned_flags) = crate::positions::calc_pctr_context_ut1_parts_in_context(
            &self.ctx,
            self.delta_t_model(),
            self.ctx.effective_iers_delta_t(),
            false,
            whole,
            tt_fraction,
            original_ut1_fraction,
            target,
            center,
            calc_flags,
            normalized_flags,
        )?;
        Ok(super::positions::pctr_position(
            self,
            values,
            returned_flags,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_parts_admit_both_aliases_before_nonfinite_frame_acquisition() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        let error = session
            .calc_pctr_context_ut1_parts(
                f64::NAN,
                f64::NAN,
                f64::NAN,
                MARS,
                EARTH,
                FLG_SWIEPH,
                FLG_SWIEPH,
            )
            .unwrap_err();
        assert!(matches!(
            error.primary_error(),
            Error::SpkUnregisteredBody { body: MARS }
        ));
        assert!(session.context().spk_registry().is_empty());
        assert!(matches!(
            session.calc_pctr(f64::NAN, MARS, EARTH, FLG_SWIEPH),
            Err(Error::NonFiniteArgument { .. })
        ));
    }
}
