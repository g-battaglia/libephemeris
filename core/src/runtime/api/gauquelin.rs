//! Gauquelin sectors over real public scalar places and horizon events.
use super::*;
use crate::houses::GauquelinRiseSetQuery;

impl Session {
    /// Fixed-star Gauquelin sector at original UT1 Julian Day `jd_ut`.
    ///
    /// Methods 0/1 acquire the configured House ARMC/true-obliquity frame,
    /// then the public V1 scalar `fixstar_ut` place with `flags|SPEED`.
    /// TOPOCTR installs an owned per-call observer after the frame; the caller's
    /// observer and retained resources are unchanged. Only method1 sets body
    /// latitude to zero. Other methods outside2–5 follow method0. Output flags
    /// are passed literally to the star place; its first two public components
    /// feed the same accepted G geometry, even for XYZ/RADIANS spellings.
    ///
    /// Methods2–5 consume real `star_rise_trans` events, preserving its catalog,
    /// configured projection and independent bundled clocks. Disc center applies
    /// to2/3; no refraction to2/4. Pressure/temperature zero become1013.25mbar/
    /// 15C before event calls. Exact next/previous rise/set searches and 1.5/2.5
    /// day lookbacks supply the diurnal/nocturnal fraction in sectors[1,37).
    /// Stars remain geometric points; event acquisitions retain their fixed
    /// flags. No method0 fallback substitutes for missing rise/set events.
    ///
    /// Known builtin numeric/search-conversion failures are caught only within
    /// methods2–5, producing the owned base Error constructor and preserving
    /// their implicit context. Catalog misses, public ranges, coordinates and
    /// configuration failures propagate in consumed order. Missing constructor
    /// data stays an explicit graph gap. Main-image acquisition uses the
    /// Session's existing owner and configured or supplied immutable resources.
    /// Geography is east longitude/latitude degrees and metres. Keys and callback
    /// borrows are synchronous; no mutable context or provider borrow escapes.
    #[allow(clippy::too_many_arguments)]
    pub fn star_gauquelin_sector(
        &self,
        key: &str,
        jd_ut: f64,
        method: i32,
        geopos: [f64; 3],
        atpress: f64,
        attemp: f64,
        flags: i32,
    ) -> Result<f64, Error> {
        let frame = |ut| crate::time::public_sidereal::house_frame_at(ut, &|u| self.deltat(u));
        let place = |ut, actual_flags| {
            let star = if actual_flags & FLG_TOPOCTR != 0 {
                let observer = Observer::new(geopos[0], geopos[1], geopos[2])?;
                self.star_ut_at_observer(key, ut, actual_flags, observer)?
            } else {
                self.fixstar_ut_captured(key, ut, actual_flags)
                    .into_catalog_result()?
            };
            Ok((star.values[0], star.values[1]))
        };
        let events = |query: GauquelinRiseSetQuery| {
            self.star_rise_trans(
                key,
                query.jd_start_ut,
                query.rsmi,
                [query.lon_deg_east, query.lat_deg, query.height_m],
                query.atpress_mbar,
                query.attemp_c,
                query.flags,
            )
            .map(|(status, values)| (status, values[0]))
        };
        star_sector_via(
            key, jd_ut, method, geopos, atpress, attemp, flags, &frame, &place, &events,
        )
    }

    /// Integer-body entry over the same unchanged sector kernels. Public
    /// `calc_ut` supplies the scalar place; public rise_trans supplies events.
    /// The established integer C symbol/provider construction is not modified.
    /// This new SDK entry retains per-call topocentric ownership and forwards
    /// the original flags/UT1. The source does not normalize method6 into an
    /// event method. No string catalog dispatch or numeric expected word enters.
    #[allow(clippy::too_many_arguments)]
    pub fn gauquelin_sector(
        &self,
        jd_ut: f64,
        body: i32,
        method: i32,
        geopos: [f64; 3],
        atpress: f64,
        attemp: f64,
        flags: i32,
    ) -> Result<f64, Error> {
        let frame = |ut| crate::time::public_sidereal::house_frame_at(ut, &|u| self.deltat(u));
        let place = |ut, actual_flags| {
            let result = if actual_flags & FLG_TOPOCTR != 0 {
                let observer = Observer::new(geopos[0], geopos[1], geopos[2])?;
                Self::from_context(self.ctx.clone().with_observer(observer)).calc_ut(
                    ut,
                    body,
                    actual_flags,
                )?
            } else {
                self.calc_ut(ut, body, actual_flags)?
            };
            Ok((result.values[0], result.values[1]))
        };
        let events = |query: GauquelinRiseSetQuery| {
            self.rise_trans(
                query.jd_start_ut,
                body,
                query.rsmi,
                geopos,
                query.atpress_mbar,
                query.attemp_c,
                query.flags,
            )
            .map(|(status, values)| (status, values[0]))
        };
        star_sector_via(
            &body.to_string(),
            jd_ut,
            method,
            geopos,
            atpress,
            attemp,
            flags,
            &frame,
            &place,
            &events,
        )
    }
}

/// Positive-divisor remainder, retaining the positive zero and IEEE NaN.
fn modulo(value: f64, divisor: f64) -> f64 {
    let value = value % divisor;
    if value == 0.0 {
        0.0
    } else if value < 0.0 {
        value + divisor
    } else {
        value
    }
}
/// Own only the genuine sector constructor and its caught implicit context.
/// An absent caught constructor stays a graph gap; no class or field is inferred
/// from a coarse category, C status, message or expected/source result row.
fn not_found(_key: &str, method: i32, caught: Option<Error>) -> Error {
    let error = Error::GauquelinRiseSetNotFound { method };
    match caught {
        Some(cause) => error.with_cause(cause),
        None => error,
    }
}
/// Source event catch families represented by acquired builtin constructors or
/// dedicated producer roles. Library lookup/range/coordinate errors remain
/// outside this catch. Unknown native-domain variants are not blanket-caught.
fn caught_event_error(error: &Error) -> bool {
    matches!(
        error.primary_error(),
        Error::MathDomain { .. }
            | Error::NumericalOverflow { .. }
            | Error::InvalidStarSearch { .. }
    )
}
/// Synchronous selection and scalar composition; callbacks own real acquisition.
/// Methods2–5 consume no frame/place callback. The existing rise/set composition
/// owns all unchanged arithmetic/limits; -1 is only its unused query tag for a
/// star key, never an ephemeris ID or constructor body identity. Ordinary frame
/// acquisition precedes the site/place callback; only method1 zeroes latitude.
#[allow(clippy::too_many_arguments)]
fn star_sector_via(
    key: &str,
    jd: f64,
    method: i32,
    geo: [f64; 3],
    pressure: f64,
    temperature: f64,
    flags: i32,
    frame: &dyn Fn(f64) -> Result<(f64, f64), Error>,
    place: &dyn Fn(f64, i32) -> Result<(f64, f64), Error>,
    event: &dyn Fn(GauquelinRiseSetQuery) -> Result<(i32, f64), Error>,
) -> Result<f64, Error> {
    if matches!(method, 2..=5) {
        let selected = |query| {
            event(query).map_err(|error| {
                if caught_event_error(&error) {
                    not_found(key, method, Some(error))
                } else {
                    error
                }
            })
        };
        return crate::houses::gauquelin_sector_from_rise_set_via(
            jd,
            -1,
            method,
            geo[1],
            geo[0],
            geo[2],
            pressure,
            temperature,
            flags,
            &selected,
        )
        .map_err(|error| {
            if matches!(error, Error::GauquelinRiseSetNotFound { .. }) {
                not_found(key, method, None)
            } else {
                error
            }
        });
    }
    let (armc0, eps) = frame(jd)?;
    let armc = modulo(armc0 + geo[0], 360.0);
    let (longitude, latitude) = place(jd, flags | FLG_SPEED)?;
    let latitude = if method == 1 { 0.0 } else { latitude };
    crate::houses::house_pos(armc, geo[1], eps, 'G', longitude, latitude)
}

#[cfg(test)]
#[path = "gauquelin_tests.rs"]
mod tests;
