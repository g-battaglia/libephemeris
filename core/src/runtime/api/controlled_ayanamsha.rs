//! Acquired TT mean offset for the reader-fast position helper.
//!
//! Declared as a child of `public_ayanamsha`; this shares the original selected
//! geometry, live-row factories and defining-epoch clocks without a reader-UT
//! table lookup or a second public wrapper conversion.

use super::*;

impl Session {
    /// Mean offset consumed by a reader-fast position at its acquired TT JD.
    ///
    /// Classical and ordinary USER anchors use this TT directly. The Source
    /// DYNAMIC15, ECL_DATE and USER_UT branches first obtain
    /// `ut = tt - configured_delta(tt)`, then independent observation TT. The
    /// consumed aberration/deflection Booleans qualify the actual live anchor;
    /// no caller flag operand is narrowed here. No reader delta_t table, main
    /// file or catalog field is consumed simply by entering the method.
    ///
    /// `configured_delta` is days and runs only for that delegated inversion.
    /// Actual live-row/clock/frame errors propagate; no mean-only/native/static
    /// replacement is installed after a refused callback.
    pub fn public_fast_mean_at_tt(
        &self,
        tt: f64,
        noaberr: bool,
        nogdefl: bool,
        rows: Option<&dyn AnchorRowAcquisition>,
        configured_delta: &dyn Fn(f64) -> Result<f64>,
    ) -> Result<(f64, Option<Warning>)> {
        let (selection, _) = crate::houses::sidereal_selection_for_context(&self.ctx);
        self.public_fast_mean_at_tt_with_selection(
            tt,
            noaberr,
            nogdefl,
            selection,
            rows,
            configured_delta,
        )
    }

    /// The same helper with the fast wrapper's actual retained mode/USER pair
    /// and current bits. Keeping this record explicit lets a Source module
    /// getter retain its outer fields while later flag callbacks replace live
    /// configuration. It does not replace the actual observation/provider owner.
    #[allow(clippy::too_many_arguments)]
    pub fn public_fast_mean_at_tt_with_selection(
        &self,
        tt: f64,
        noaberr: bool,
        nogdefl: bool,
        selection: SiderealSelection,
        rows: Option<&dyn AnchorRowAcquisition>,
        configured_delta: &dyn Fn(f64) -> Result<f64>,
    ) -> Result<(f64, Option<Warning>)> {
        let delegated = DYNAMIC_AYANAMSHA_MODES.contains(&selection.mode())
            || selection.bits() & SIDBIT_ECL_DATE != 0
            || selection.mode() == SIDM_USER && selection.bits() & SIDBIT_USER_UT != 0;
        let observation_tt = if delegated {
            let ut = tt - configured_delta(tt)?;
            crate::time::default_timescale_tt_from_ut1(ut)?
        } else {
            tt
        };
        let bits = FLG_NONUT
            | if noaberr { FLG_NOABERR } else { 0 }
            | if nogdefl { FLG_NOGDEFL } else { 0 };
        self.public_ayanamsha_selected_at_tt(observation_tt, bits, true, selection, rows)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn classical_reader_tt_does_not_consume_unused_inverse_or_image() {
        let session = Session::from_context(
            crate::context::EphemerisContext::new()
                .with_mode(crate::context::ProviderMode::NativeJpl),
        );
        let value = session.public_fast_mean_at_tt(2451545.0, true, true, None, &|_| {
            panic!("classical reader helper consumes no delegated inverse")
        });
        assert!(value.is_ok());
        assert!(session.context().spk_registry().is_empty());
    }
}
