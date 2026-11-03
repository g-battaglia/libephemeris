//! Ordered main-image acquisition for scalar-TT planetocentric calculations.
//!
//! Planetocentric calculations have no stored fast-attempt/logging handoff and
//! no Horizons-vector leg. Explicit sealed calculations keep their stored
//! vector provider. Auto/Horizons calculations select their actual main image,
//! then reuse the existing owned native vector calculation unchanged.

use super::*;
use crate::positions::{MainImageAcquisition, MainImageOwner};
use crate::types::TwoPartEpoch;

impl Session {
    /// Planetocentric position with a retained current main-image lifecycle.
    ///
    /// `jd_tt` is finite scalar TT, and `target`/`center` are standard planet
    /// IDs. Output is degrees/degrees/AU and daily rates; XYZ/RADIANS, SPEED,
    /// astrometric center-bit mapping and sidereal frames retain the existing
    /// planetocentric engine. SealedLeb uses its explicit stored-vector path
    /// and never invokes `main_image`. Other modes consume the actual current
    /// main provider directly; no LEB fast retry or Horizons vector request is
    /// introduced. The owner remains alive through all split samples.
    ///
    /// Finite-epoch refusal precedes acquisition. For a standard endpoint the
    /// selected-main envelope is checked before cached TT construction; for
    /// two unsupported endpoints, TT construction precedes main acquisition
    /// and their eventual unsupported-body refusal. A present/cached owner
    /// retains current configuration and unrelated registrations. Actual
    /// opener, range, alias, scalar/split-domain and numerical errors propagate.
    pub fn calc_pctr_with_main_image(
        &self,
        jd_tt: f64,
        target: i32,
        center: i32,
        flags: i32,
        main_image: &dyn MainImageAcquisition,
    ) -> Result<Position, Error> {
        self.calc_pctr_dispatch(jd_tt, target, center, flags, None, Some(main_image))
    }

    /// Planetocentric main acquisition with live rows only at the existing
    /// final sidereal tail. Reader/main acquisition and body SPEED geometry
    /// precede independent ayanamsha clocks; no row is prefetched here.
    #[allow(clippy::too_many_arguments)]
    pub fn calc_pctr_with_anchor_rows_and_main_image(
        &self,
        jd_tt: f64,
        target: i32,
        center: i32,
        flags: i32,
        rows: Option<&dyn crate::stars::AnchorRowAcquisition>,
        main_image: &dyn MainImageAcquisition,
    ) -> Result<Position, Error> {
        self.calc_pctr_dispatch(jd_tt, target, center, flags, rows, Some(main_image))
    }

    /// Shared Session dispatch. Explicit NativeJpl callers without a new main
    /// capability retain their already supplied registry contract; otherwise
    /// only a successful actual main acquisition identifies the native image.
    /// The strict low-level stored/native numerical factories are unchanged.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn calc_pctr_dispatch(
        &self,
        jd_tt: f64,
        target: i32,
        center: i32,
        flags: i32,
        rows: Option<&dyn crate::stars::AnchorRowAcquisition>,
        main_image: Option<&dyn MainImageAcquisition>,
    ) -> Result<Position, Error> {
        let anchors = rows.map(|rows| PositionAnchors {
            session: self,
            rows: Some(rows),
        });
        let anchors = anchors
            .as_ref()
            .map(|owner| owner as &dyn crate::sidereal::DynamicAnchorProvider);
        if self.ctx.mode() == ProviderMode::SealedLeb
            || (self.ctx.mode() == ProviderMode::NativeJpl
                && main_image.is_none()
                && !self.ctx.spk_registry().is_empty())
        {
            let (values, returned_flags) = positions::calc_pctr_in_context_with_anchors(
                &self.ctx,
                self.delta_t_model(),
                self.ctx.effective_iers_delta_t(),
                false,
                jd_tt,
                target,
                center,
                flags,
                anchors,
            )?;
            return Ok(pctr_position(self, values, returned_flags));
        }
        let owner = ordered_pctr_main(
            jd_tt,
            target,
            center,
            || {
                let acquired = match main_image {
                    Some(acquire) => acquire.acquire(&self.ctx)?,
                    None => positions::acquire_current_main_image(&self.ctx)?,
                };
                // Even a caller-supplied older receipt contributes only its
                // actual image/cache. All settings and minor mappings belong
                // to this current call, not the receipt's former Session.
                acquired.with_current_context(self.ctx.clone())
            },
            |owner| validate_main_envelope(owner, jd_tt, target),
            |tt| scalar_tt_constructor(tt).map(|_| ()),
        )?;
        let native = owner
            .native_provider_context()?
            .with_mode(ProviderMode::NativeJpl);
        let (values, returned_flags) = positions::calc_pctr_in_context_with_anchors(
            &native,
            self.delta_t_model(),
            self.ctx.effective_iers_delta_t(),
            false,
            jd_tt,
            target,
            center,
            flags,
            anchors,
        )?;
        Ok(pctr_position(self, values, returned_flags))
    }
}

/// Conditional public prelude without an eager frame/UT1/TDB conversion.
/// Standard endpoints cause main range validation before cached TT; the
/// otherwise unsupported pair reaches cached TT before its first main read.
/// Each acquisition runs once; any earlier failure prevents all later stages.
fn ordered_pctr_main<T>(
    jd_tt: f64,
    target: i32,
    center: i32,
    acquire: impl FnOnce() -> Result<T, Error>,
    validate: impl FnOnce(&T) -> Result<(), Error>,
    construct_tt: impl FnOnce(f64) -> Result<(), Error>,
) -> Result<T, Error> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    if positions::is_pctr_supported_body(target) || positions::is_pctr_supported_body(center) {
        let owner = acquire()?;
        validate(&owner)?;
        construct_tt(jd_tt)?;
        Ok(owner)
    } else {
        construct_tt(jd_tt)?;
        acquire()
    }
}

/// One-argument scalar TT constructor: integer day and retained fraction.
/// This performs no independent UT1, Delta-T, TDB, frame or data read. Fresh
/// vector/SPEED Time constructions remain with the existing numerical engine.
fn scalar_tt_constructor(tt: f64) -> Result<TwoPartEpoch, Error> {
    if !tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: tt });
    }
    let whole = tt.floor();
    Ok(TwoPartEpoch {
        hi: whole,
        lo: tt - whole,
    })
}

/// Range metadata belongs to the actual selected main, independently of a
/// target alias. Unavailable zero bounds skip the guard. Calendar constructor
/// failures precede the public typed range result in requested/start/end order.
fn validate_main_envelope(owner: &MainImageOwner, tt: f64, target: i32) -> Result<(), Error> {
    let session = Session::from_context(
        owner
            .native_provider_context()?
            .with_mode(ProviderMode::NativeJpl),
    );
    let metadata = session.current_file_metadata(0, Some(owner.registry_key()))?;
    let Some(crate::engine_api::CurrentFileMetadata::Main {
        start_jd, end_jd, ..
    }) = metadata
    else {
        return Ok(());
    };
    if start_jd == 0.0 && end_jd == 0.0 {
        return Ok(());
    }
    if tt < start_jd || tt > end_jd {
        return Err(Error::EphemerisRange {
            requested_jd: tt,
            start_jd,
            end_jd,
            body: Some(target),
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::{Cell, RefCell};

    #[test]
    fn standard_public_range_then_tt_then_computation_order() {
        let calls = RefCell::new(Vec::new());
        let owner = ordered_pctr_main(
            2451545.125,
            MARS,
            EARTH,
            || {
                calls.borrow_mut().push("main for range");
                Ok(17)
            },
            |owner| {
                assert_eq!(*owner, 17);
                calls.borrow_mut().push("main range");
                Ok(())
            },
            |tt| {
                assert_eq!(tt.to_bits(), 2451545.125f64.to_bits());
                calls.borrow_mut().push("cached TT constructor");
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(owner, 17);
        assert_eq!(
            *calls.borrow(),
            ["main for range", "main range", "cached TT constructor"]
        );
    }

    #[test]
    fn unsupported_pair_time_failure_precedes_first_main_acquisition() {
        let marker = Error::ConfigurationMissing {
            missing_config: "actual scalar Time factory",
            suggestion: "control failure",
        };
        let actual = ordered_pctr_main::<()>(
            2451545.0,
            90001,
            90002,
            || panic!("a failed earlier Time constructor cannot open main"),
            |_| panic!("unsupported pairs have no public main range prelude"),
            |_| Err(marker.clone()),
        )
        .unwrap_err();
        assert_eq!(actual, marker);
    }

    #[test]
    fn actual_range_refusal_prevents_tt_constructor() {
        let marker = Error::EphemerisRange {
            requested_jd: 100.0,
            start_jd: 200.0,
            end_jd: 300.0,
            body: Some(MARS),
        };
        let actual = ordered_pctr_main(
            100.0,
            MARS,
            90001,
            || Ok(()),
            |_| Err(marker.clone()),
            |_| panic!("public range refusal precedes cached TT"),
        )
        .unwrap_err();
        assert_eq!(actual, marker);
    }

    #[test]
    fn nonfinite_public_tt_never_acquires_main_or_constructs_time() {
        for jd in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                ordered_pctr_main::<()>(
                    jd,
                    MARS,
                    EARTH,
                    || panic!("nonfinite main"),
                    |_| panic!("nonfinite range"),
                    |_| panic!("nonfinite TT"),
                ),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
    }

    #[test]
    fn scalar_tt_constructor_does_not_evaluate_an_independent_clock() {
        let tt = scalar_tt_constructor(2451545.125).unwrap();
        assert_eq!(tt.hi.to_bits(), 2451545.0f64.to_bits());
        assert_eq!(tt.lo.to_bits(), 0.125f64.to_bits());
    }

    #[test]
    fn auto_and_horizons_pctr_use_main_owner_without_remote_or_leb_attempt() {
        for mode in [ProviderMode::Auto, ProviderMode::HorizonsRemote] {
            let session = Session::from_context(
                EphemerisContext::new()
                    .with_mode(mode)
                    .with_network_policy(Some(crate::context::NetworkPolicy::Sealed))
                    .with_delta_t_model(Some(crate::time::DeltaTModel::Smh2016)),
            );
            let count = Cell::new(0);
            let marker = Error::ConfigurationMissing {
                missing_config: "actual main acquisition",
                suggestion: "test owner",
            };
            let main = |current: &EphemerisContext| -> Result<MainImageOwner, Error> {
                count.set(count.get() + 1);
                assert_eq!(current.mode(), mode);
                Err(marker.clone())
            };
            let actual = session
                .calc_pctr_with_main_image(2451545.0, MARS, EARTH, FLG_SWIEPH, &main)
                .unwrap_err();
            assert_eq!(actual, marker);
            assert_eq!(count.get(), 1);
        }
    }

    #[test]
    fn sealed_pctr_never_consumes_an_attached_main_capability() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::SealedLeb));
        let main = |_: &EphemerisContext| -> Result<MainImageOwner, Error> {
            panic!("explicit sealed path cannot acquire JPL")
        };
        assert!(
            session
                .calc_pctr_with_main_image(2451545.0, MARS, EARTH, FLG_SWIEPH, &main)
                .is_err()
        );
    }

    #[test]
    fn selected_main_envelope_refuses_before_an_unsupported_target_alias() {
        // Project-authored descriptor fixture only: no astronomical expected
        // vector is created or consumed. The actual main metadata reader
        // supplies the envelope independently of target registrations.
        let mut loader = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        loader
            .spk_registry_mut()
            .register_image("selected actual main", crate::spk::moon_test_kernel_bytes())
            .unwrap();
        let current = EphemerisContext::new().with_mode(ProviderMode::Auto);
        let owner = MainImageOwner::from_loader_registry(
            current.clone(),
            loader.spk_registry(),
            "selected actual main",
            "requested.bsp".into(),
            "selected.bsp".into(),
            false,
        )
        .unwrap();
        let metadata = Session::from_context(
            owner
                .native_provider_context()
                .unwrap()
                .with_mode(ProviderMode::NativeJpl),
        )
        .current_file_metadata(0, Some(owner.registry_key()))
        .unwrap()
        .unwrap();
        let crate::engine_api::CurrentFileMetadata::Main { end_jd, .. } = metadata else {
            panic!("genuine owner must expose its main image metadata")
        };
        let calls = Cell::new(0);
        let main = |_: &EphemerisContext| {
            calls.set(calls.get() + 1);
            Ok(owner.clone())
        };
        let actual = Session::from_context(current)
            .calc_pctr_with_main_image(end_jd + 1.0, 90001, EARTH, FLG_SWIEPH, &main)
            .unwrap_err();
        assert_eq!(calls.get(), 1);
        assert!(matches!(
            actual.primary_error(),
            Error::EphemerisRange {
                requested_jd,
                body: Some(90001),
                ..
            } if requested_jd.to_bits() == (end_jd + 1.0).to_bits()
        ));
    }
}
