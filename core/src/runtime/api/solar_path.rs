//! Public sampled solar shadow paths and fixed-epoch cone width. Besselian
//! scalar reads, public where and actual WGS84 rays retain distinct clocks.
use super::solar_circumstances::{PositionResolver, SelectedMain, solar_position, solar_selected};
use super::*;
use crate::events::{
    self as solar, BesselianCore, EclipseLimitEphemeris, EclipsePathEphemeris, SolarPathWidthSource,
};

/// Owned `(UT1 sample JDs, north latitudes, east longitudes)` component arrays.
/// Lengths agree; samples appear in walk order, with omitted points absent from
/// all three arrays. Coordinates are geodetic degrees, longitude in[-180,180).
pub type EclipsePath = (Vec<f64>, Vec<f64>, Vec<f64>);

/// Solar shadow-cone width in kilometres.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CapturedSolarPathWidth {
    /// Width, in kilometres, including the model's horizon projection and cap.
    pub value: f64,
}

impl Session {
    /// Complete a native cone-width continuation after apparent acquisition.
    ///
    /// `altitude` is unrefracted solar altitude in degrees, already above the
    /// horizon; Sun/Moon distances are AU. The late callback returns apparent
    /// separation in degrees, after angular-radius arithmetic. The result is
    /// kilometres. No file, Time or observer
    /// is acquired, and no observer-stage catch covers this arithmetic.
    pub fn eclipse_path_width_apparent_geometry(
        altitude: f64,
        sun_distance_au: f64,
        moon_distance_au: f64,
        separation: &dyn Fn() -> Result<f64, Error>,
    ) -> Result<CapturedSolarPathWidth, Error> {
        solar::solar_path_width_geometry(altitude, sun_distance_au, moon_distance_au, separation)
    }

    /// Calculate the solar shadow width in kilometres.
    pub fn calc_eclipse_path_width_captured(
        &self,
        jd: f64,
        lat: Option<f64>,
        lon: Option<f64>,
        flags: i32,
    ) -> Result<CapturedSolarPathWidth, Error> {
        self.solar_path_width_captured(jd, lat, lon, flags, None)
    }
    /// Captured counterpart with configuration resolved only at the actual
    /// Native main-image boundary. Ownership/error/clock stages are unchanged.
    pub fn calc_eclipse_path_width_captured_with_position_context(
        &self,
        jd: f64,
        lat: Option<f64>,
        lon: Option<f64>,
        flags: i32,
        resolve: &PositionResolver<'_>,
    ) -> Result<CapturedSolarPathWidth, Error> {
        self.solar_path_width_captured(jd, lat, lon, flags, Some(resolve))
    }
    /// Dispatch the shared width model over the selected provider.
    fn solar_path_width_captured(
        &self,
        jd: f64,
        lat: Option<f64>,
        lon: Option<f64>,
        flags: i32,
        resolve: Option<&PositionResolver<'_>>,
    ) -> Result<CapturedSolarPathWidth, Error> {
        if self.ctx.mode() != ProviderMode::NativeJpl {
            return Ok(CapturedSolarPathWidth {
                value: self.solar_path_width(jd, lat, lon, flags, resolve)?,
            });
        }
        solar::solar_path_width_captured_source_via(
            &NativeSolarPath::new(self, flags, resolve),
            jd,
            lat,
            lon,
        )
    }
    /// Sample the central solar shadow track on an inclusive UT1 JD interval.
    /// step_minutes must be finite, positive and greater than the binary64
    /// half-spacing limit throughout the window. The accumulator adds the day
    /// step repeatedly; no endpoint interpolation or sample-count repair occurs.
    /// Only CENTRAL samples with gamma<1.5 and a nonzero coordinate survive.
    /// Public where precedes Besselian x then y. Provider/clock failures retain
    /// their owned errors; invalid windows/steps refuse before native reads.
    pub fn calc_eclipse_central_line(
        &self,
        start: f64,
        end: f64,
        step_minutes: f64,
        flags: i32,
    ) -> Result<EclipsePath, Error> {
        self.solar_path_walk(start, end, step_minutes, flags, 0, None)
    }
    /// Sample the northern umbral/antumbral boundary, using x/y/d/mu/l2 in
    /// that order, WGS84 flattening and y+abs(l2). Grazing samples may be
    /// omitted; accuracy degrades at the ends. Interval/clock/error contracts
    /// and component units follow calc_eclipse_central_line.
    pub fn calc_eclipse_northern_limit(
        &self,
        start: f64,
        end: f64,
        step_minutes: f64,
        flags: i32,
    ) -> Result<EclipsePath, Error> {
        self.solar_path_walk(start, end, step_minutes, flags, 1, None)
    }
    /// Southern boundary, with y-abs(l2) and the same surface/geodetic gates.
    /// Returns the three owned component arrays, never a geographic polygon.
    /// Sampling, units, convergence and errors follow the northern boundary.
    pub fn calc_eclipse_southern_limit(
        &self,
        start: f64,
        end: f64,
        step_minutes: f64,
        flags: i32,
    ) -> Result<EclipsePath, Error> {
        self.solar_path_walk(start, end, step_minutes, flags, 2, None)
    }
    /// Width in kilometres of the total/annular cone at UT1 JD. An explicit
    /// site is used only when both latitude and longitude are Some; otherwise
    /// public where determines the central site and noncentral events return0.
    /// The true WGS84 Sun-then-Moon apparent pair and unrefracted horizon use
    /// bundled Time, independent of configured model/frame Delta-T. Caller
    /// flags affect only the public where branch. Below-horizon, off-track or
    /// exhausted-umbra samples return0. Width is capped at1500km; this local
    /// vertical-plane projection is approximate, especially near grazing ends.
    /// Main/Time/body errors propagate; selected observer/observation primitive
    /// errors return0, while altitude and subsequent arithmetic errors propagate.
    pub fn calc_eclipse_path_width(
        &self,
        jd: f64,
        lat: Option<f64>,
        lon: Option<f64>,
        flags: i32,
    ) -> Result<f64, Error> {
        self.solar_path_width(jd, lat, lon, flags, None)
    }
    /// Central sampler with a configuration-only resolver at the first actual
    /// Native main-image read. Supplied registries keep priority and ownership.
    pub fn calc_eclipse_central_line_with_position_context(
        &self,
        start: f64,
        end: f64,
        minutes: f64,
        flags: i32,
        resolve: &PositionResolver<'_>,
    ) -> Result<EclipsePath, Error> {
        self.solar_path_walk(start, end, minutes, flags, 0, Some(resolve))
    }
    /// Northern sampler with the same consumed-stage main resolver. No state,
    /// sample words, astronomical clock or observer is supplied by the callback.
    pub fn calc_eclipse_northern_limit_with_position_context(
        &self,
        start: f64,
        end: f64,
        minutes: f64,
        flags: i32,
        resolve: &PositionResolver<'_>,
    ) -> Result<EclipsePath, Error> {
        self.solar_path_walk(start, end, minutes, flags, 1, Some(resolve))
    }
    /// Southern sampler with call-scoped configuration acquisition; immutable
    /// images and cache results remain owned until the complete walk returns.
    pub fn calc_eclipse_southern_limit_with_position_context(
        &self,
        start: f64,
        end: f64,
        minutes: f64,
        flags: i32,
        resolve: &PositionResolver<'_>,
    ) -> Result<EclipsePath, Error> {
        self.solar_path_walk(start, end, minutes, flags, 2, Some(resolve))
    }
    /// Width with consumed-stage configuration acquisition; no eager finite
    /// epoch/site admission is added before the selected public branch.
    pub fn calc_eclipse_path_width_with_position_context(
        &self,
        jd: f64,
        lat: Option<f64>,
        lon: Option<f64>,
        flags: i32,
        resolve: &PositionResolver<'_>,
    ) -> Result<f64, Error> {
        self.solar_path_width(jd, lat, lon, flags, Some(resolve))
    }
    /// Keep sealed servings unchanged and use native acquisition.
    /// Validate native intervals before an otherwise unused main or clock.
    fn solar_path_walk(
        &self,
        start: f64,
        end: f64,
        minutes: f64,
        flags: i32,
        edge: u8,
        resolve: Option<&PositionResolver<'_>>,
    ) -> Result<EclipsePath, Error> {
        if self.ctx.mode() != ProviderMode::NativeJpl {
            return self.with_serving_stack(|stack| {
                let dt = |ut| crate::time::deltat_with_config(ut, stack.delta_cfg);
                let bound = solar::SealedSolarGlobalServing::try_new(stack.serving, &dt, flags)?;
                match edge {
                    0 => solar::sol_eclipse_central_line_sealed(&bound, start, end, minutes),
                    1 => solar::sol_eclipse_northern_limit_sealed(&bound, start, end, minutes),
                    _ => solar::sol_eclipse_southern_limit_sealed(&bound, start, end, minutes),
                }
            });
        }
        let native = NativeSolarPath::new(self, flags, resolve);
        match edge {
            0 => solar::eclipse_central_line_via(&native, start, end, minutes),
            _ => solar::solar_limit_source_via(&native, edge == 2, start, end, minutes),
        }
    }
    /// Stored width keeps its checked serving. Native width stages distinguish
    /// uncaught resources/Time/body lookup from caught observer/apparent rays.
    fn solar_path_width(
        &self,
        jd: f64,
        lat: Option<f64>,
        lon: Option<f64>,
        flags: i32,
        resolve: Option<&PositionResolver<'_>>,
    ) -> Result<f64, Error> {
        if self.ctx.mode() != ProviderMode::NativeJpl {
            return self.with_serving_stack(|stack| {
                let dt = |ut| crate::time::deltat_with_config(ut, stack.delta_cfg);
                let bound = solar::SealedSolarGlobalServing::try_new(stack.serving, &dt, flags)?;
                solar::sol_eclipse_path_width_sealed(&bound, jd, lat, lon)
            });
        }
        solar::solar_path_width_source_via(
            &NativeSolarPath::new(self, flags, resolve),
            jd,
            lat,
            lon,
        )
    }
}

/// Per-call owned main and Besselian cache. Decimal9/low8bits keys do not
/// round actual acquisition epochs. Callback execution occurs without a cache
/// borrow; only successful complete nine-word results enter the cache.
struct NativeSolarPath<'a> {
    session: &'a Session,
    flags: i32,
    resolve: Option<&'a PositionResolver<'a>>,
    selected: SelectedMain,
    cache: std::cell::RefCell<std::collections::HashMap<(u64, i32), BesselianCore>>,
}
impl<'a> NativeSolarPath<'a> {
    /// Construct empty call-local storage without native resource acquisition.
    fn new(session: &'a Session, flags: i32, resolve: Option<&'a PositionResolver<'a>>) -> Self {
        Self {
            session,
            flags,
            resolve,
            selected: SelectedMain::new(),
            cache: Default::default(),
        }
    }
    /// One complete Besselian read: Sun then Moon spherical equatorial scalar
    /// calc_ut using only ephemeris bits, Cartesian conversion, then public
    /// sidtime. Cache clear occurs only when more than512 entries existed.
    fn besselian(&self, jd: f64) -> Result<BesselianCore, Error> {
        let rounded = jd;
        let rounded = if rounded == 0.0 { 0.0 } else { rounded };
        let key = (rounded.to_bits(), self.flags & 0xff);
        if let Some(value) = self.cache.borrow().get(&key).copied() {
            return Ok(value);
        }
        let flags = (self.flags & (FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH)) | FLG_EQUATORIAL;
        let sun = solar_position(self.session, jd, SUN, flags, self.resolve, &self.selected)?;
        let moon = solar_position(self.session, jd, MOON, flags, self.resolve, &self.selected)?;
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
        let value = solar::besselian_core_lazy(cart(sun), cart(moon), &mut || {
            self.session.sidtime(jd).map(|h| h * 15.0)
        })?;
        let mut cache = self.cache.borrow_mut();
        if cache.len() > 512 {
            cache.clear();
        }
        cache.insert(key, value);
        Ok(value)
    }
    /// Public where including its center local stage; no retflag-only shortcut.
    fn where_point(&self, jd: f64) -> Result<(i32, f64, f64), Error> {
        let (flags, gp, _) = if let Some(resolve) = self.resolve {
            self.session
                .sol_eclipse_where_with_position_context(jd, self.flags, resolve)?
        } else {
            self.session.sol_eclipse_where(jd, self.flags)?
        };
        Ok((flags, gp[0], gp[1]))
    }
}
impl EclipsePathEphemeris for NativeSolarPath<'_> {
    fn path_where(&self, jd: f64) -> Result<(i32, f64, f64), Error> {
        self.where_point(jd)
    }
    fn path_besselian_xy(&self, jd: f64) -> Result<(f64, f64), Error> {
        let x = self.besselian(jd)?.x();
        let y = self.besselian(jd)?.y();
        Ok((x, y))
    }
}
impl EclipseLimitEphemeris for NativeSolarPath<'_> {
    fn limit_besselian(&self, jd: f64) -> Result<(f64, f64, f64, f64, f64), Error> {
        let x = self.besselian(jd)?.x();
        let y = self.besselian(jd)?.y();
        let d = self.besselian(jd)?.d_deg();
        let mu = self.besselian(jd)?.mu_deg();
        let l2 = self.besselian(jd)?.l2();
        Ok((x, y, d, mu, l2))
    }
}
impl SolarPathWidthSource for NativeSolarPath<'_> {
    type Sample = crate::positions::context_serving::NativeSolarPathRays;
    fn width_where(&self, jd: f64) -> Result<(i32, f64, f64), Error> {
        self.where_point(jd)
    }
    fn width_acquire(&self, jd: f64, lon: f64, lat: f64) -> Result<Option<Self::Sample>, Error> {
        let context = solar_selected(self.session, SUN, self.resolve, &self.selected)?;
        let time = crate::positions::context_serving::native_solar_path_time(&context.ctx, jd)?;
        match crate::positions::context_serving::native_solar_path_rays(
            &context.ctx,
            &time,
            lon,
            lat,
        ) {
            Ok(sample) => Ok(Some(sample)),
            Err(error) if width_observer_failure(&error) => Ok(None),
            Err(error) => Err(error),
        }
    }
    fn width_altitude(&self, sample: &Self::Sample) -> Result<f64, Error> {
        sample.altitude()
    }
    fn width_distances(&self, sample: &Self::Sample) -> (f64, f64) {
        sample.distances()
    }
    fn width_separation(&self, sample: &Self::Sample) -> f64 {
        sample.separation()
    }
}

/// Classify only the protected observer/apparent exception hierarchy.
/// LEB coverage/corruption remain terminal. An actual constructor precedes
/// coarse native families, so a library Error never becomes a numeric zero.
/// A raw SPK observer-range ValueError belongs to this catch; Main/Time/alias
/// failures do not pass through this predicate.
fn width_observer_failure(error: &Error) -> bool {
    matches!(
        error.primary_error(),
        Error::NonFiniteArgument { .. }
            | Error::InvalidObserver { .. }
            | Error::MathDomain { .. }
            | Error::NumericalOverflow { .. }
            | Error::IntegerOverflow
            | Error::RationalDivisionByZero
            | Error::LebBodyNotFound { .. }
            | Error::SpkOutOfRange { .. }
    )
}

#[cfg(test)]
#[path = "solar_path/tests.rs"]
mod tests;
