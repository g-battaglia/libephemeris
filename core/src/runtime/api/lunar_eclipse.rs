//! Context-owned lunar-eclipse circumstances at an original UT1/site query.
use super::*;
use crate::events;

/// Configuration-only callback borrowed for one synchronous circumstance call.
/// The input snapshot is shared; the returned context owns its images. No
/// position, clock value, input reference or callback escapes the operation.
type EclipsePositionResolver<'a> =
    dyn Fn(&EphemerisContext, i32) -> Result<EphemerisContext, Error> + 'a;

impl Session {
    /// Lunar-eclipse circumstances at UT1 JD and an east-longitude/latitude
    /// degree, height-metre site. Returns character bits and 20 attributes.
    /// Slots 0/1 are dimensionless umbral/penumbral magnitudes, 2/3 are zero,
    /// 4 is south-westward Moon azimuth, 5/6 true/refracted altitude in degrees,
    /// 7 opposition distance in degrees (zero outside an eclipse), 8 repeats slot 0,
    /// 9/10 are lunar Saros series/member and 11–19 are zero. The character
    /// bits are total, partial, penumbral or zero; apparent altitude <=0 clears
    /// only the bits and leaves every attribute filled. No visibility bits
    /// are added. Refraction uses estimated pressure and 10 Celsius.
    ///
    /// NativeJpl geometry consumes public [`Self::calc_ut`] Moon then Sun,
    /// with only ephemeris-selection bits plus EQUATORIAL|XYZ|SPEED. The
    /// Danjon shadow kernel then uses these geocentric AU vectors and Saros
    /// at the original UT1. Site acquisition independently constructs bundled
    /// global UT1 Time, WGS84 Earth observer and strict apparent Sun then Moon
    /// rays in the true ecliptic of date. It consumes no caller correction or
    /// frame flags and does not substitute public spherical TOPOCTR places.
    /// [`Self::azalt`] projects the Moon at the original UT1: configured Delta-T
    /// supplies sidereal terms; bundled global TT supplies the ecliptic frame.
    /// Scalar public TT windows and raw split-SPK directory ranges remain
    /// separate. The one-try bit is stripped before ephemeris-bit selection.
    ///
    /// Nonfinite caller epoch/site inputs refuse before acquisition. Native
    /// finite sites retain WGS84 arithmetic without a new latitude bound;
    /// body, support, main-image, range, frame, model and ten-step light-time
    /// errors propagate in consumed order. Explicit LEB/Auto retains the
    /// existing sealed local serving and flag/site gates. No other eclipse
    /// search or certified physical-accuracy guarantee is supplied here.
    pub fn lun_eclipse_how(
        &self,
        jd_ut: f64,
        geopos: [f64; 3],
        flags: i32,
    ) -> Result<(i32, [f64; 20]), Error> {
        self.lun_eclipse_how_inner(jd_ut, geopos, flags, None)
    }

    /// Circumstances with lazily resolved NativeJpl main-image configuration.
    /// Arguments, clocks, layout and errors follow [`Self::lun_eclipse_how`].
    /// `resolve` is consumed before each genuine public Moon/Sun calc_ut and
    /// before the direct WGS84 pair's main-image acquisition (Sun key). It
    /// returns configuration only, never a replacement position or epoch.
    /// Explicit supplied registries stay owned snapshots; a binding may fill
    /// process-selected paths only for an otherwise empty registry. No caller
    /// observer/context is changed. Other provider modes do not call it.
    pub fn lun_eclipse_how_with_position_context(
        &self,
        jd_ut: f64,
        geopos: [f64; 3],
        flags: i32,
        resolve: &dyn Fn(&EphemerisContext, i32) -> Result<EphemerisContext, Error>,
    ) -> Result<(i32, [f64; 20]), Error> {
        self.lun_eclipse_how_inner(jd_ut, geopos, flags, Some(resolve))
    }

    /// Dispatch checked caller inputs before borrowing a provider. LEB uses
    /// its unchanged local serving; Native keeps public scalar body acquisition
    /// and independent original-Time site acquisition as separate callbacks.
    fn lun_eclipse_how_inner(
        &self,
        jd_ut: f64,
        geopos: [f64; 3],
        flags: i32,
        resolve: Option<&EclipsePositionResolver<'_>>,
    ) -> Result<(i32, [f64; 20]), Error> {
        for value in [jd_ut, geopos[0], geopos[1], geopos[2]] {
            if !value.is_finite() {
                return Err(Error::NonFiniteArgument { value });
            }
        }
        let flags = events::strip_one_try_bit(flags);
        if self.ctx.mode() != ProviderMode::NativeJpl {
            return self.with_serving_stack(|stack| {
                let dt = |ut| crate::time::deltat_with_config(ut, stack.delta_cfg);
                let bound = events::SealedLunarLocalServing::try_new(
                    stack.serving,
                    &dt,
                    flags,
                    geopos[0],
                    geopos[1],
                    geopos[2],
                )?;
                events::LunarLocalEphemeris::how_attr(&bound, jd_ut)
            });
        }
        let position = |body, calc_flags| -> Result<[f64; 3], Error> {
            let context = resolve_eclipse_context(&self.ctx, body, resolve)?;
            let result = Self::from_context(context).calc_ut(jd_ut, body, calc_flags)?;
            Ok([result.values[0], result.values[1], result.values[2]])
        };
        let topo = || {
            let context = resolve_eclipse_context(&self.ctx, SUN, resolve)?;
            crate::positions::context_serving::native_lunar_eclipse_topocentric_pair(
                &context, jd_ut, geopos,
            )
        };
        let horizontal = |moon: [f64; 3]| {
            let alt = self.azalt(
                jd_ut,
                events::ECL2HOR,
                geopos[0],
                geopos[1],
                geopos[2],
                0.0,
                10.0,
                moon[0],
                moon[1],
            )?;
            Ok((alt.azimuth, alt.true_altitude, alt.apparent_altitude))
        };
        native_lunar_how_via(jd_ut, flags, &position, &topo, &horizontal)
    }
}

/// Own resolved configuration at a consumed position/main-image boundary.
/// The optional callback is synchronous and receives a shared snapshot;
/// absence clones that snapshot without environment or filesystem reads.
fn resolve_eclipse_context(
    context: &EphemerisContext,
    body: i32,
    resolve: Option<&EclipsePositionResolver<'_>>,
) -> Result<EphemerisContext, Error> {
    match resolve {
        Some(resolve) => resolve(context, body),
        None => Ok(context.clone()),
    }
}

/// Compose consumed geocentric and site callbacks without changing their
/// clocks. `position` returns equatorial-of-date Cartesian AU at original
/// UT1; `topo` returns apparent Sun/Moon (longitude_deg,latitude_deg,AU).
/// Shadow geometry/Saros precedes site reads even when its character is zero.
/// `horizontal` uses the Moon only; the Sun is nevertheless acquired first.
/// Callback errors abort before later reads; horizon gating preserves attrs.
fn native_lunar_how_via(
    jd_ut: f64,
    flags: i32,
    position: &dyn Fn(i32, i32) -> Result<[f64; 3], Error>,
    topo: &dyn Fn() -> Result<([f64; 3], [f64; 3]), Error>,
    horizontal: &dyn Fn([f64; 3]) -> Result<(f64, f64, f64), Error>,
) -> Result<(i32, [f64; 20]), Error> {
    let flags =
        (flags & (FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH)) | FLG_EQUATORIAL | FLG_XYZ | FLG_SPEED;
    let moon = position(MOON, flags)?;
    let sun = position(SUN, flags)?;
    let shadow = events::lun_how_core_via(
        (moon[0], moon[1], moon[2]),
        (sun[0], sun[1], sun[2]),
        events::saros_info(jd_ut, events::EclipseKind::Lunar),
    )?;
    let (_sun, moon) = topo()?;
    Ok(events::lun_eclipse_how_parts(&shadow, horizontal(moon)?))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;

    #[test]
    fn native_lunar_how_flags_and_acquisition_order() {
        let reads = RefCell::new(Vec::new());
        let position = |body, flags| {
            assert_eq!(flags, FLG_SWIEPH | FLG_EQUATORIAL | FLG_XYZ | FLG_SPEED);
            reads.borrow_mut().push(body);
            Ok(if body == MOON {
                [-0.00257, 0.0, 0.0]
            } else {
                [1.0, 0.0, 0.0]
            })
        };
        let topo = || {
            reads.borrow_mut().extend([10, 11]);
            Ok(([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]))
        };
        let horizontal = |moon| {
            assert_eq!(moon, [4.0, 5.0, 6.0]);
            reads.borrow_mut().push(12);
            Ok((100.0, -1.0, 0.0))
        };
        let (kind, attr) = native_lunar_how_via(
            2451545.125,
            FLG_SWIEPH | FLG_J2000 | FLG_SIDEREAL | FLG_NOABERR | FLG_NOGDEFL | FLG_TOPOCTR,
            &position,
            &topo,
            &horizontal,
        )
        .unwrap();
        assert_eq!(*reads.borrow(), [MOON, SUN, 10, 11, 12]);
        assert_eq!(kind, 0);
        assert!(attr[0] > 0.0);
        assert_eq!(attr[0].to_bits(), attr[8].to_bits());
        assert_eq!(&attr[4..7], &[100.0, -1.0, 0.0]);
        for i in [2, 3, 11, 12, 13, 14, 15, 16, 17, 18, 19] {
            assert_eq!(attr[i].to_bits(), 0);
        }
    }

    #[test]
    fn native_lunar_how_first_position_error_precedes_site() {
        let position = |body, _| Err(Error::UnknownBody { body_id: body });
        let topo = || -> Result<([f64; 3], [f64; 3]), Error> { panic!("unused site") };
        let horizontal = |_| -> Result<(f64, f64, f64), Error> { panic!("unused horizon") };
        assert!(matches!(
            native_lunar_how_via(2451545.125, 2, &position, &topo, &horizontal),
            Err(Error::UnknownBody { body_id: MOON })
        ));
    }

    #[test]
    fn native_lunar_how_no_eclipse_still_consumes_topocentric_pair() {
        let position = |body, _| {
            Ok(if body == MOON {
                [0.0, 0.00257, 0.0]
            } else {
                [1.0, 0.0, 0.0]
            })
        };
        let topo = || Err(Error::SpkUnregisteredBody { body: 599 });
        let horizontal = |_| -> Result<(f64, f64, f64), Error> { panic!("unused horizon") };
        assert!(matches!(
            native_lunar_how_via(2451545.125, 2, &position, &topo, &horizontal),
            Err(Error::SpkUnregisteredBody { body: 599 })
        ));
    }

    #[test]
    fn native_lunar_how_nonfinite_precedes_resolver() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        let resolve = |_: &EphemerisContext, _| -> Result<EphemerisContext, Error> {
            panic!("unused configuration")
        };
        assert!(matches!(
            session.lun_eclipse_how_with_position_context(f64::NAN, [0.0; 3], 2, &resolve),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            session.lun_eclipse_how_with_position_context(
                2451545.125,
                [0.0, 0.0, f64::INFINITY],
                2,
                &resolve
            ),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn native_lunar_how_resolver_consumes_moon_before_unused_site() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        let resolve = |context: &EphemerisContext, body| {
            assert!(context.spk_registry().is_empty());
            assert_eq!(body, MOON);
            Err(Error::ConfigurationMissing {
                missing_config: "test_main",
                suggestion: "supply image",
            })
        };
        assert!(matches!(
            session.lun_eclipse_how_with_position_context(
                2451545.125,
                [0.0, 91.0, 0.0],
                2,
                &resolve
            ),
            Err(Error::ConfigurationMissing {
                missing_config: "test_main",
                ..
            })
        ));
    }
}
