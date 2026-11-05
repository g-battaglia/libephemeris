//! Instantaneous, observer-independent lunar shadow scalars at UT1.
use super::*;
use crate::events;

/// Operation selection keeps canonical Cartesian shadow acquisition separate
/// from the small-angle gamma acquisition over public spherical positions.
#[derive(Clone, Copy)]
enum LunarScalar {
    Umbral,
    Penumbral,
    Gamma,
}

/// A configuration-only resolver borrowed synchronously at each consumed
/// public position call. Returned contexts own their registry snapshots.
type PositionContext<'a> = dyn Fn(&EphemerisContext, i32) -> Result<EphemerisContext, Error> + 'a;

impl Session {
    /// Fraction of the Moon's diameter immersed in Earth's umbral shadow at
    /// UT1 JD. Zero means no umbral immersion, 0–1 partial, >=1 total.
    /// The canonical selenocentric cone uses Moon then Sun public calc_ut at
    /// this original UT1, ephemeris bits plus EQUATORIAL|XYZ|SPEED and AU
    /// vectors. Danjon enlarges the lunar-parallax term by one percent;
    /// negative magnitudes are floored to positive zero. No site is acquired.
    /// Caller correction/frame flags are discarded; clocks and range/error
    /// acquisition follow each public calc_ut. Nonfinite input refuses before
    /// configuration. LEB/Auto retain sealed local serving and flag gates.
    /// This is an instantaneous scalar, without event search or accuracy bound.
    pub fn lun_eclipse_umbral_magnitude(&self, jd_ut: f64, flags: i32) -> Result<f64, Error> {
        self.lunar_scalar(jd_ut, flags, LunarScalar::Umbral, None)
    }

    /// Fraction of the lunar diameter immersed in the penumbral cone, floored
    /// to positive zero. Values >1 mean full immersion. Units, clocks, flags,
    /// acquisition, domain and errors follow [`Self::lun_eclipse_umbral_magnitude`];
    /// only the canonical penumbral channel differs.
    pub fn lun_eclipse_penumbral_magnitude(&self, jd_ut: f64, flags: i32) -> Result<f64, Error> {
        self.lunar_scalar(jd_ut, flags, LunarScalar::Penumbral, None)
    }

    /// Non-negative small-angle shadow-axis distance in Earth radii at UT1 JD.
    /// Every provider mode consumes Sun then Moon public calc_ut with fixed SPEED flags,
    /// irrespective of caller flags. Their ecliptic degrees and AU distances
    /// are used directly: sqrt((delta_lon*cos(latitude))²+latitude²) divided
    /// by atan(4.2635e-5/moon_distance) in degrees. Longitude difference is
    /// folded once at ±180 degrees; the near-zero angular-radius guard returns
    /// positive zero. The complete small-angle classification arithmetic runs
    /// before selecting gamma. No separate frame clock, obliquity or site is
    /// acquired. Each public calc_ut owns its provider clocks and ranges:
    /// stored bodies retain configured fast model TT; NativeJpl retains its
    /// independent bundled timescale.
    /// Nonfinite caller UT1 refuses before configuration; acquisition/arithmetic
    /// errors propagate in order. Explicit LEB retains its stored reader and refuses
    /// a missing channel without substituting JPL. This approximation supplies no
    /// certified physical-accuracy bound.
    pub fn lun_eclipse_gamma(&self, jd_ut: f64, flags: i32) -> Result<f64, Error> {
        self.lunar_scalar(jd_ut, flags, LunarScalar::Gamma, None)
    }

    /// Umbral magnitude with a lazily resolved, owned position context.
    /// The resolver supplies configuration only at each Moon then Sun call;
    /// argument, output, clocks, errors and flags follow the ordinary method.
    /// The callback is consumed at each actual public position acquisition.
    pub fn lun_eclipse_umbral_magnitude_with_position_context(
        &self,
        jd_ut: f64,
        flags: i32,
        resolve: &dyn Fn(&EphemerisContext, i32) -> Result<EphemerisContext, Error>,
    ) -> Result<f64, Error> {
        self.lunar_scalar(jd_ut, flags, LunarScalar::Umbral, Some(resolve))
    }

    /// Penumbral counterpart of the configuration-only umbral resolver route.
    pub fn lun_eclipse_penumbral_magnitude_with_position_context(
        &self,
        jd_ut: f64,
        flags: i32,
        resolve: &dyn Fn(&EphemerisContext, i32) -> Result<EphemerisContext, Error>,
    ) -> Result<f64, Error> {
        self.lunar_scalar(jd_ut, flags, LunarScalar::Penumbral, Some(resolve))
    }

    /// Gamma with lazy configuration at the Sun then Moon public calc_ut
    /// boundaries. Resolver lifetime and ownership match the magnitude routes;
    /// gamma retains its distinct fixed-SPEED spherical acquisition in every
    /// provider mode. The callback is consumed synchronously; each returned
    /// context owns its images and settings, and no callback borrow escapes.
    pub fn lun_eclipse_gamma_with_position_context(
        &self,
        jd_ut: f64,
        flags: i32,
        resolve: &dyn Fn(&EphemerisContext, i32) -> Result<EphemerisContext, Error>,
    ) -> Result<f64, Error> {
        self.lunar_scalar(jd_ut, flags, LunarScalar::Gamma, Some(resolve))
    }

    /// Check caller UT1 before dispatch. Gamma always acquires public spherical
    /// positions; native magnitudes acquire public Cartesian positions. Each
    /// call owns its configuration snapshot. Stored magnitudes consume the same public Moon then Sun Cartesian
    /// acquisitions, preserving the caller UT1 in sealed range wrappers.
    fn lunar_scalar(
        &self,
        jd_ut: f64,
        flags: i32,
        operation: LunarScalar,
        resolve: Option<&PositionContext<'_>>,
    ) -> Result<f64, Error> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let flags = events::strip_one_try_bit(flags);
        let position = |body, calc_flags| {
            let context = match resolve {
                Some(resolve) => resolve(&self.ctx, body)?,
                None => self.ctx.clone(),
            };
            Ok(Self::from_context(context)
                .calc_ut(jd_ut, body, calc_flags)?
                .values)
        };
        native_lunar_scalar_via(jd_ut, flags, operation, &position)
    }
}

/// Consume positions in operation-specific order, with no site or frame read.
/// Native magnitudes select only ephemeris bits; gamma ignores caller flags
/// entirely in every mode. The original UT1 is captured by the caller; clocks
/// and reader selection belong to each public position acquisition. No extra
/// obliquity, observer or vector conversion is performed for gamma.
/// The callback retains all six public position words without reconstruction.
fn native_lunar_scalar_via(
    jd_ut: f64,
    flags: i32,
    operation: LunarScalar,
    position: &dyn Fn(i32, i32) -> Result<[f64; 6], Error>,
) -> Result<f64, Error> {
    if matches!(operation, LunarScalar::Gamma) {
        let sun = position(SUN, FLG_SPEED)?;
        let moon = position(MOON, FLG_SPEED)?;
        return events::lun_gamma_ecliptic_via(sun, moon);
    }
    let base =
        (flags & (FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH)) | FLG_EQUATORIAL | FLG_XYZ | FLG_SPEED;
    let moon = position(MOON, base)?;
    let sun = position(SUN, base)?;
    let moon = (moon[0], moon[1], moon[2]);
    let sun = (sun[0], sun[1], sun[2]);
    let saros = events::saros_info(jd_ut, events::EclipseKind::Lunar);
    match operation {
        LunarScalar::Umbral => events::lun_umbral_magnitude_via(moon, sun, saros),
        LunarScalar::Penumbral => events::lun_penumbral_magnitude_via(moon, sun, saros),
        LunarScalar::Gamma => unreachable!(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;

    #[test]
    fn native_lunar_scalar_distinct_flags_order_and_no_extra_reads() {
        for operation in [LunarScalar::Umbral, LunarScalar::Penumbral] {
            let reads = RefCell::new(Vec::new());
            let position = |body, flags| {
                reads.borrow_mut().push(body);
                assert_eq!(flags, FLG_SWIEPH | FLG_EQUATORIAL | FLG_XYZ | FLG_SPEED);
                Ok(if body == MOON {
                    [-0.00257, 0.0, 0.0, 0.0, 0.0, 0.0]
                } else {
                    [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                })
            };
            assert!(
                native_lunar_scalar_via(
                    2451545.125,
                    FLG_SWIEPH | FLG_TOPOCTR | FLG_SIDEREAL | FLG_RADIANS,
                    operation,
                    &position
                )
                .unwrap()
                    > 1.0
            );
            assert_eq!(*reads.borrow(), [MOON, SUN]);
        }
        let reads = RefCell::new(Vec::new());
        let position = |body, flags| {
            assert_eq!(flags, FLG_SPEED);
            reads.borrow_mut().push(body);
            Ok(if body == SUN {
                [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
            } else {
                [180.0, 0.0, 0.00257, 0.0, 0.0, 0.0]
            })
        };
        assert_eq!(
            native_lunar_scalar_via(2451545.125, i32::MAX, LunarScalar::Gamma, &position)
                .unwrap()
                .to_bits(),
            0
        );
        assert_eq!(*reads.borrow(), [SUN, MOON]);
    }

    #[test]
    fn native_lunar_scalar_first_error_stops_second_acquisition() {
        for operation in [
            LunarScalar::Umbral,
            LunarScalar::Penumbral,
            LunarScalar::Gamma,
        ] {
            let reads = RefCell::new(Vec::new());
            let position = |body, _| {
                reads.borrow_mut().push(body);
                Err(Error::UnknownBody { body_id: body })
            };
            let first = if matches!(operation, LunarScalar::Gamma) {
                SUN
            } else {
                MOON
            };
            assert!(
                matches!(native_lunar_scalar_via(2451545.125, 2, operation, &position),
                Err(Error::UnknownBody { body_id }) if body_id == first)
            );
            assert_eq!(*reads.borrow(), [first]);
        }
    }

    #[test]
    fn native_lunar_scalar_nonfinite_precedes_configuration_and_gamma_starts_sun() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        let refused = |_: &EphemerisContext, body| Err(Error::UnknownBody { body_id: body });
        for operation in [
            LunarScalar::Umbral,
            LunarScalar::Penumbral,
            LunarScalar::Gamma,
        ] {
            assert!(matches!(
                session.lunar_scalar(f64::NAN, 2, operation, Some(&refused)),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
        assert!(matches!(
            session.lun_eclipse_gamma_with_position_context(2451545.125, 2, &refused),
            Err(Error::UnknownBody { body_id: SUN })
        ));
        assert!(matches!(
            session.lun_eclipse_umbral_magnitude_with_position_context(2451545.125, 2, &refused),
            Err(Error::UnknownBody { body_id: MOON })
        ));
    }
}

#[cfg(test)]
mod matched_leb_public_gamma_tests {
    use super::*;
    use std::cell::RefCell;

    #[test]
    fn gamma_all_modes_ignore_flags_before_public_position_resolution() {
        for mode in [
            ProviderMode::SealedLeb,
            ProviderMode::Auto,
            ProviderMode::NativeJpl,
            ProviderMode::HorizonsRemote,
        ] {
            let session = Session::from_context(EphemerisContext::new().with_mode(mode));
            for flags in [0, 2, 32770, 79442, -1] {
                let reads = RefCell::new(Vec::new());
                let resolve = |_: &EphemerisContext, body| {
                    reads.borrow_mut().push(body);
                    Err(Error::UnknownBody { body_id: body })
                };
                assert!(matches!(
                    session.lun_eclipse_gamma_with_position_context(2451545.125, flags, &resolve),
                    Err(Error::UnknownBody { body_id: SUN })
                ));
                assert_eq!(*reads.borrow(), [SUN]);
            }
        }
    }

    #[test]
    fn gamma_all_modes_nonfinite_ut_precedes_resolver() {
        for mode in [
            ProviderMode::SealedLeb,
            ProviderMode::Auto,
            ProviderMode::NativeJpl,
            ProviderMode::HorizonsRemote,
        ] {
            let session = Session::from_context(EphemerisContext::new().with_mode(mode));
            let resolve = |_: &EphemerisContext, _| -> Result<EphemerisContext, Error> {
                panic!("nonfinite caller UT1 must precede configuration")
            };
            for jd in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                assert!(matches!(
                    session.lun_eclipse_gamma_with_position_context(jd, 79442, &resolve),
                    Err(Error::NonFiniteArgument { .. })
                ));
            }
        }
    }

    #[test]
    fn gamma_sealed_position_error_stops_before_moon_and_has_no_jpl_substitute() {
        let context = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        let expected = Session::from_context(context.clone())
            .calc_ut(2451545.125, SUN, FLG_SPEED)
            .unwrap_err();
        let reads = RefCell::new(Vec::new());
        let resolve = |ctx: &EphemerisContext, body| {
            reads.borrow_mut().push(body);
            Ok(ctx.clone())
        };
        let actual = Session::from_context(context)
            .lun_eclipse_gamma_with_position_context(2451545.125, 79442, &resolve)
            .unwrap_err();
        assert_eq!(format!("{actual:?}"), format!("{expected:?}"));
        assert_eq!(*reads.borrow(), [SUN]);
    }

    #[test]
    fn gamma_original_before_month_source_positions_select_exact_source_word() {
        // Literal words are the independently recorded public source inputs,
        // not reconstructed orbital states or native output expectations.
        let sun = [
            0x40717dbd4114ed8a,
            0x3f2e1987f89cbf1d,
            0x3fef7774620b87d1,
            0x3ff04f641301a000,
            0xbec488cf6bc6de00,
            0xbee41a2627700000,
        ]
        .map(f64::from_bits);
        let moon = [
            0x406b29310a24e495,
            0x4014ece1cda96197,
            0x3f65f3fe8b53bcf7,
            0x402834dd015f8e00,
            0xbfb04231f9937600,
            0x3ef72044afb85800,
        ]
        .map(f64::from_bits);
        for flags in [0, 2, 32770, 79442, -1] {
            let reads = RefCell::new(Vec::new());
            let position = |body, actual_flags| {
                assert_eq!(actual_flags, FLG_SPEED);
                reads.borrow_mut().push(body);
                Ok(if body == SUN { sun } else { moon })
            };
            let gamma = native_lunar_scalar_via(
                f64::from_bits(0x4142b42c3fef9db2),
                flags,
                LunarScalar::Gamma,
                &position,
            )
            .unwrap();
            assert_eq!(gamma.to_bits(), 0x40600d9787d552cf);
            assert_eq!(*reads.borrow(), [SUN, MOON]);
        }
    }
}
