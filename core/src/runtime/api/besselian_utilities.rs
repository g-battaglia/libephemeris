//! UT1 Besselian shadow words over public scalar position acquisition.
use super::solar_circumstances::{PositionResolver, SelectedMain, solar_position};
use super::*;
use crate::events::{self, BesselianCore, BesselianRates};

/// Besselian scalar helpers consume caller-supplied anchor values.
/// No ephemeris acquisition or alternative arithmetic backend is involved.
pub struct EclipseCycleAnchors;

impl EclipseCycleAnchors {
    /// UT1 eclipse JD, Saros series and member, in strict first-minimum order.
    pub fn saros(kind: events::EclipseKind) -> &'static [(f64, i32, i32)] {
        match kind {
            events::EclipseKind::Solar => crate::eclipse_tables::SOLAR_SAROS_REFERENCES,
            events::EclipseKind::Lunar => crate::eclipse_tables::LUNAR_SAROS_REFERENCES,
        }
    }

    /// UT1 eclipse JD and Inex series, in strict first-minimum order.
    pub fn inex(kind: events::EclipseKind) -> &'static [(f64, i32)] {
        match kind {
            events::EclipseKind::Solar => crate::eclipse_tables::SOLAR_INEX_REFERENCES,
            events::EclipseKind::Lunar => crate::eclipse_tables::LUNAR_INEX_REFERENCES,
        }
    }
}

impl Session {
    /// Nine Besselian words at UT1 JD: shadow x/y in Earth radii, declination
    /// and Greenwich hour angle in degrees, penumbral/umbral radii, cone
    /// tangents and the Moon's axial coordinate in Earth radii.
    ///
    /// Sun then Moon public spherical calc_ut retain only caller ephemeris
    /// bits plus EQUATORIAL. Both places precede Cartesian conversion and
    /// axis/projection arithmetic; public sidtime is consumed afterward,
    /// before cone arithmetic. Native main selection is retained for the
    /// query. Stored calculations retain their scalar model clocks and reader
    /// owners. This checked geometry preserves its existing finite domain;
    /// arbitrary raw vector/operator protocols are not represented here.
    pub fn besselian_core_at_ut(&self, jd: f64, flags: i32) -> Result<BesselianCore, Error> {
        self.besselian_core_inner(jd, flags, None, &SelectedMain::new())
    }

    /// Acquire shadow words and retain the exact main-image owner selected by
    /// this sample. A following sample can use the returned Session without
    /// reopening the main image. Selection remains lazy: cache hits need not
    /// call this method, and a position failure returns before any later
    /// sample. Stored and already populated Native sessions retain their
    /// existing images and configuration unchanged.
    pub fn besselian_core_with_owner_at_ut(
        &self,
        jd: f64,
        flags: i32,
    ) -> Result<(BesselianCore, Self), Error> {
        let selected = SelectedMain::new();
        let value = self.besselian_core_inner(jd, flags, None, &selected)?;
        Ok((value, selected.into_inner().unwrap_or_else(|| self.clone())))
    }

    /// Besselian words with configuration resolved only when the Native main
    /// image is actually consumed. Returned contexts own all selected images;
    /// the callback is borrowed synchronously and never runs under a cache
    /// borrow. Stored modes use their retained Session directly.
    pub fn besselian_core_at_ut_with_position_context(
        &self,
        jd: f64,
        flags: i32,
        resolve: &PositionResolver<'_>,
    ) -> Result<BesselianCore, Error> {
        self.besselian_core_inner(jd, flags, Some(resolve), &SelectedMain::new())
    }

    /// Six per-hour rates from one-minute minus then plus samples at UT1 JD.
    /// x/y/l1/l2 rates are Earth radii/hour and d/mu rates are degrees/hour.
    /// One selected main owner survives both samples. The existing centered
    /// quotient divides by twice the step then by24, folding mu at±180°.
    /// Position, geometry and clock errors propagate at their consumed stage;
    /// no derivative or physical-accuracy bound is inferred from the stencil.
    pub fn besselian_rates_at_ut(&self, jd: f64, flags: i32) -> Result<BesselianRates, Error> {
        let selected = SelectedMain::new();
        let step = 1.0 / 1440.0;
        let minus = self.besselian_core_inner(jd - step, flags, None, &selected)?;
        let plus = self.besselian_core_inner(jd + step, flags, None, &selected)?;
        events::besselian_rates_via(&minus, &plus, step)
    }

    /// Dispatch only the acquisition channel; physical Besselian geometry is
    /// shared unchanged. No eager epoch or GAST preflight is introduced.
    fn besselian_core_inner(
        &self,
        jd: f64,
        flags: i32,
        resolve: Option<&PositionResolver<'_>>,
        selected: &SelectedMain,
    ) -> Result<BesselianCore, Error> {
        let position = |body, flags| {
            if self.ctx.mode() == ProviderMode::NativeJpl {
                solar_position(self, jd, body, flags, resolve, selected)
            } else {
                self.calc_ut(jd, body, flags).map(|value| value.values)
            }
        };
        besselian_from_positions(flags, &position, &mut || {
            self.sidtime_with(jd, 0.0, None, None).map(|h| h * 15.0)
        })
    }
}

/// Acquire both public places before reducing either one. Only the geometry
/// kernel consumes the lazy GAST callback, after its axis and projections.
fn besselian_from_positions(
    flags: i32,
    position: &dyn Fn(i32, i32) -> Result<[f64; 6], Error>,
    gast: &mut impl FnMut() -> Result<f64, Error>,
) -> Result<BesselianCore, Error> {
    let base = (flags & (FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH)) | FLG_EQUATORIAL;
    let sun = position(SUN, base)?;
    let moon = position(MOON, base)?;
    let cart = |p: [f64; 6]| {
        let ra = p[0].to_radians();
        let dec = p[1].to_radians();
        let r = p[2];
        (
            r * dec.cos() * ra.cos(),
            r * dec.cos() * ra.sin(),
            r * dec.sin(),
        )
    };
    events::besselian_core_lazy(cart(sun), cart(moon), gast)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;

    #[test]
    fn both_places_precede_the_consumed_gast_clock() {
        let trace = RefCell::new(Vec::new());
        let position = |body, flags| {
            trace.borrow_mut().push(body);
            assert_eq!(flags, FLG_SWIEPH | FLG_EQUATORIAL);
            Ok(if body == SUN {
                [0., 0., 1., 0., 0., 0.]
            } else {
                [0.1, 0.1, 0.00257, 0., 0., 0.]
            })
        };
        let mut gast = || {
            trace.borrow_mut().push(-1);
            Err(Error::UnknownBody { body_id: 918 })
        };
        let error =
            besselian_from_positions(FLG_SWIEPH | FLG_SPEED | FLG_J2000, &position, &mut gast)
                .unwrap_err();
        assert!(matches!(error, Error::UnknownBody { body_id: 918 }));
        assert_eq!(*trace.borrow(), [SUN, MOON, -1]);
    }

    #[test]
    fn a_later_moon_failure_precedes_geometry_of_the_sun() {
        let position = |body, _| {
            if body == SUN {
                Ok([f64::NAN, 0., 1., 0., 0., 0.])
            } else {
                Err(Error::UnknownBody { body_id: 919 })
            }
        };
        let error =
            besselian_from_positions(2, &position, &mut || panic!("unused GAST")).unwrap_err();
        assert!(matches!(error, Error::UnknownBody { body_id: 919 }));
    }
}
