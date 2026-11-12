//! Numerical output reduction for scalar-clock Horizons states.
//!
//! The frame and ayanamsha callbacks are deliberately independent: frame
//! samples consume scalar TT, while ayanamsha samples consume UT and flags.
//! Callers must supply the configured selection and actual dynamic anchors;
//! this module neither substitutes a JPL clock nor fabricates anchor states.

use super::routing::model_serving::ModelRateOrigin;
use super::*;

/// IAU 2000A nutation and biased Vondrák precession at scalar TT.
/// Angles are radians; the TT Julian day must be finite.
pub(crate) fn scalar_frame(tt: f64) -> Result<HorizonsFrame> {
    let (psi, deps) = reference_frames::iau2000a_tenth_microarcsec(tt)?;
    let psi = psi * reference_frames::FRAME_TENTH_USEC_2_RAD;
    let deps = deps * reference_frames::FRAME_TENTH_USEC_2_RAD;
    let (pn, eps) = vondrak_pn_matrix(tt, psi, deps, true)?;
    Ok((pn, psi, eps, vondrak_precession_matrix(tt, true)?))
}

/// Precess spherical degrees/AU and per-day rates between date and J2000.
/// Angular rates use a centered half-day stencil with varying date epoch;
/// distance/radial rate are unchanged. Frame/domain errors propagate.
pub(crate) fn local_precess_state(
    mut s: PositionResult,
    tt: f64,
    to_j2000: bool,
) -> Result<PositionResult> {
    let transform = |lon, lat, epoch| {
        let (from, to) = if to_j2000 {
            (epoch, 2451545.0)
        } else {
            (2451545.0, epoch)
        };
        crate::hypothetical::precess_ecliptic(lon, lat, from, to)
    };
    let (lon, lat) = transform(s[0], s[1], tt)?;
    let (lm, bm) = transform(s[0] - s[3] * 0.5, s[1] - s[4] * 0.5, tt - 0.5)?;
    let (lp, bp) = transform(s[0] + s[3] * 0.5, s[1] + s[4] * 0.5, tt + 0.5)?;
    let mut dl = lp - lm;
    if dl > 180.0 {
        dl -= 360.0;
    } else if dl < -180.0 {
        dl += 360.0;
    }
    s[0] = lon;
    s[1] = lat;
    s[3] = dl;
    s[4] = bp - bm;
    Ok(s)
}

/// TT derivative of the actual Vondrák pole-angle mean obliquity, deg/day.
/// Each periodic/polynomial pole term is differentiated analytically; for
/// epsilon=acos(e·q), epsilon_dot=-(e_dot·q+e·q_dot)/sin(epsilon).
fn mean_obliquity_rate(tt: f64) -> Result<f64> {
    let t = (tt - crate::constants::J2000) / 36525.0;
    let pair_rate = |polynomial: &[(f64, f64); 4],
                     periods: &[f64],
                     cosines: [&[f64]; 2],
                     sines: [&[f64]; 2]| {
        let mut rates = [0.0; 2];
        for component in 0..2 {
            let coefficient = |i: usize| {
                if component == 0 {
                    polynomial[i].0
                } else {
                    polynomial[i].1
                }
            };
            rates[component] = (3.0 * coefficient(3))
                .mul_add(t, 2.0 * coefficient(2))
                .mul_add(t, coefficient(1));
            for (i, period) in periods.iter().enumerate() {
                let frequency = std::f64::consts::TAU / period;
                let (sine, cosine) = (frequency * t).sin_cos();
                rates[component] +=
                    frequency * (-sine * cosines[component][i] + cosine * sines[component][i]);
            }
            rates[component] *= ARCSEC_TO_RAD / 36525.0;
        }
        rates
    };
    let [pd, qd] = pair_rate(
        &VONDRAK_PQPOL,
        &VONDRAK_PQPER_PERIODS,
        [&VONDRAK_PQPER_P_COS, &VONDRAK_PQPER_Q_COS],
        [&VONDRAK_PQPER_P_SIN, &VONDRAK_PQPER_Q_SIN],
    );
    let [xd, yd] = pair_rate(
        &VONDRAK_XYPOL,
        &VONDRAK_XYPER_PERIODS,
        [&VONDRAK_XYPER_X_COS, &VONDRAK_XYPER_Y_COS],
        [&VONDRAK_XYPER_X_SIN, &VONDRAK_XYPER_Y_SIN],
    );
    let e = super::corrections::ecliptic_pole_vondrak(tt);
    let q = super::corrections::equator_pole_vondrak(tt);
    let (sine, cosine) = VONDRAK_EPS0_RAD.sin_cos();
    let ep = e[0];
    let eq = -e[1] * cosine - e[2] * sine;
    let ez = -e[1] * sine + e[2] * cosine;
    let ezd = -(ep * pd + eq * qd) / ez;
    let ed = [pd, -qd * cosine - ezd * sine, -qd * sine + ezd * cosine];
    let zd = -(q[0] * xd + q[1] * yd) / q[2];
    let qd = [xd, yd, zd];
    let dot_rate: f64 = (0..3).map(|i| ed[i] * q[i] + e[i] * qd[i]).sum();
    let epsilon = vondrak_mean_obliquity_rad(tt)?;
    let rate = -dot_rate / epsilon.sin();
    if !rate.is_finite() {
        return Err(Error::NonFiniteArgument { value: rate });
    }
    Ok(rate.to_degrees())
}

/// Apply date-to-J2000 precession unless already framed, then equatorial
/// conversion with fixed J2000 or Vondrák-2011 of-date obliquity. The latter
/// includes IAU nutation unless `NONUT` or `SIDEREAL` selects the mean.
/// Units are degrees,
/// AU and per-day rates. The equatorial velocity includes the analytic TT
/// derivative of mean/true obliquity; J2000 obliquity is fixed. Undefined angular
/// positions pass through unchanged.
pub(crate) fn local_convert(
    mut s: PositionResult,
    tt: f64,
    flags: i32,
    already_j2000: bool,
) -> Result<PositionResult> {
    if !s[0].is_finite() || !s[1].is_finite() {
        return Ok(s);
    }
    if flags & FLG_J2000 != 0 && !already_j2000 {
        if flags & FLG_SPEED != 0 {
            s = local_precess_state(s, tt, true)?;
        } else {
            (s[0], s[1]) = crate::hypothetical::precess_ecliptic(s[0], s[1], tt, 2451545.0)?;
        }
    }
    if flags & FLG_EQUATORIAL == 0 {
        return Ok(s);
    }
    let eps = if flags & FLG_J2000 != 0 {
        84381.406 / 3600.0
    } else {
        vondrak_mean_obliquity_deg(tt)?
            + if flags & (FLG_NONUT | FLG_SIDEREAL) == 0 {
                crate::erfa::nut06a(tt)?.1.to_degrees()
            } else {
                0.0
            }
    };
    let (p, mut v) = crate::utils::cotrans_sp(s[0], s[1], s[2], s[3], s[4], s[5], -eps)?;
    if flags & FLG_SPEED != 0 && flags & FLG_J2000 == 0 {
        let mut epsilon_rate = mean_obliquity_rate(tt)?;
        if flags & (FLG_NONUT | FLG_SIDEREAL) == 0 {
            epsilon_rate += crate::erfa::nut06a_obliquity_rate(tt)?.to_degrees();
        }
        // d(R_x(epsilon) p)/dt = R_x v + epsilon_dot * (e_x cross p).
        // The spherical Jacobian adds these RA/Dec rates without changing
        // the rotation-invariant radius or radial derivative.
        let (sin_ra, cos_ra) = p[0].to_radians().sin_cos();
        v[0] -= epsilon_rate * cos_ra * p[1].to_radians().tan();
        v[1] += epsilon_rate * sin_ra;
    }
    Ok([p[0], p[1], p[2], v[0], v[1], v[2]])
}

/// Local-model sidereal tail. All modes subtract the actual flag-aware UT
/// ayanamsha's wrapped one-second centered derivative, not polynomial drift.
pub(crate) fn local_sidereal(
    mut s: PositionResult,
    ut: f64,
    flags: i32,
    aya: &dyn Fn(f64, i32) -> Result<f64>,
) -> Result<PositionResult> {
    s[0] = (s[0] - aya(ut, flags)?).rem_euclid(360.0);
    if flags & FLG_SPEED != 0 {
        let dt = 1.0 / 86400.0;
        let previous = aya(ut - dt, flags)?;
        let next = aya(ut + dt, flags)?;
        s[3] -= ((next - previous + 180.0).rem_euclid(360.0) - 180.0) / (2.0 * dt);
    }
    Ok(s)
}

/// Fully bound reducers needing only caller-owned scale/anchor acquisition.
/// Missing dynamic anchors propagate; the callback must honor flags and the
/// configured snapshotted sidereal selection.
pub struct HorizonsNumerics<'a> {
    delta_t: &'a dyn Fn(f64) -> Result<f64>,
    aya: &'a dyn Fn(f64, i32) -> Result<f64>,
    star_based: bool,
    coverage: crate::lunar::LunarCoverage,
}

impl<'a> HorizonsNumerics<'a> {
    /// Bind Delta-T in days, flag-aware ayanamsha in degrees, base mode and
    /// explicit lunar coverage. Construction acquires no astronomical data.
    pub fn new(
        delta_t: &'a dyn Fn(f64) -> Result<f64>,
        aya: &'a dyn Fn(f64, i32) -> Result<f64>,
        base_mode: i32,
        coverage: crate::lunar::LunarCoverage,
    ) -> Self {
        Self {
            delta_t,
            aya,
            coverage,
            star_based: crate::sidereal::DYNAMIC_AYANAMSHA_MODES.contains(&base_mode),
        }
    }

    /// Reduce corrected ICRF AU/AU-day vectors using caller scalar TT/UT.
    pub fn output(
        &self,
        p: [f64; 3],
        v: [f64; 3],
        tt: f64,
        ut: f64,
        flags: i32,
    ) -> Result<(PositionResult, i32)> {
        HorizonsOutput {
            frame_at: &scalar_frame,
            ayanamsha_at: self.aya,
            star_based: self.star_based,
        }
        .output(p, v, tt, ut, flags)
    }

    /// Reduce a native-center local model, with representation shaping and
    /// the analytical-only final-representation SPEED3 stencil.
    pub fn local(&self, ut: f64, tt: f64, body: i32, flags: i32) -> Result<(PositionResult, i32)> {
        self.local_with_state(ut, tt, body, flags, None, None)
    }

    /// Reduce a supplied J2000-ecliptic hypothetical center state.
    /// Angles are degrees, distance AU and rates per day. The caller performs
    /// center composition and the model's stencil; frame/nutation/sidereal and
    /// representation reduction are identical to the native-center model.
    /// No observer parallax applies: the caller either serves a geocentric
    /// place or owns topocentric reduction separately.
    #[cfg(test)]
    pub(crate) fn hypothetical_state(
        &self,
        ut: f64,
        tt: f64,
        body: i32,
        flags: i32,
        state: PositionResult,
    ) -> Result<(PositionResult, i32)> {
        self.local_with_state(ut, tt, body, flags, Some(state), None)
    }

    /// Reduce a supplied hypothetical center state with geocentric parallax.
    ///
    /// `offset_j2000` is the observer geocentric offset (topocentre minus
    /// geocentre) and its rate in J2000-ecliptic AU and AU/day. The
    /// of-date-oriented offset subtracts
    /// on the tropical of-date place after nutation and before any sidereal
    /// rotation, and a requested J2000 frame precesses the topocentric place
    /// afterwards in the equatorial-conversion tail. Thus, for a topocentric request the early J2000 precession is
    /// skipped so the state stays of-date through the parallax, and the
    /// conversion runs with `already_j2000 = false` so it carries the late
    /// date-to-J2000 precession. Heliocentric/barycentric requests stay
    /// inert and non-topocentric requests are bitwise identical to
    /// [`HorizonsNumerics::hypothetical_state_using_ut`]. Angles are degrees,
    /// distance AU, rates per day; non-finite inputs refuse.
    #[cfg(test)]
    pub(crate) fn hypothetical_state_with_topo(
        &self,
        ut: f64,
        tt: f64,
        body: i32,
        flags: i32,
        state: PositionResult,
        offset_j2000: ([f64; 3], [f64; 3]),
    ) -> Result<(PositionResult, i32)> {
        self.local_with_state(ut, tt, body, flags, Some(state), Some(offset_j2000))
    }

    /// Reduce a supplied hypothetical state with UT1 acquired on demand.
    ///
    /// `ut_at` returns the entry's scalar UT1 JD and is invoked only at an
    /// ecliptic sidereal reduction, after model frame/nutation/parallax work.
    /// Tropical and equatorial outputs acquire no unused clock. State and
    /// optional J2000 observer offset use AU and per-day rates; angles are
    /// degrees. The ordinary hypothetical frame/format/error order is retained.
    pub(crate) fn hypothetical_state_using_ut(
        &self,
        ut_at: &dyn Fn() -> Result<f64>,
        tt: f64,
        body: i32,
        flags: i32,
        state: PositionResult,
        offset_j2000: Option<([f64; 3], [f64; 3])>,
    ) -> Result<(PositionResult, i32)> {
        self.local_with_state_using_ut(ut_at, tt, body, flags, Some(state), offset_j2000)
    }

    /// Reduce an internally acquired hypothetical state after epoch validation.
    /// Internal TT and cached-nutation arithmetic retain IEEE values, while
    /// the original finite epoch remains owned by the calculation entry. Units,
    /// lazy UT1 use, frame/parallax/format order and output flags are unchanged.
    /// Other local numerical entries retain their checked nutation provider.
    pub(crate) fn hypothetical_source_time_state_using_ut(
        &self,
        ut_at: &dyn Fn() -> Result<f64>,
        tt: f64,
        body: i32,
        flags: i32,
        state: PositionResult,
        offset_j2000: Option<([f64; 3], [f64; 3])>,
    ) -> Result<(PositionResult, i32)> {
        if tt.is_finite() {
            return self.hypothetical_state_using_ut(ut_at, tt, body, flags, state, offset_j2000);
        }
        self.hypothetical_source_time_state_with_rate_origin(
            ut_at,
            tt,
            body,
            flags,
            state,
            offset_j2000,
            ModelRateOrigin::LegacyStencil,
        )
    }

    /// Reduce an acquired state with its actual velocity method. Only a
    /// successful owned derivative admits instantaneous nutation; a declined
    /// refinement retains the original combined model/frame stencil.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn hypothetical_source_time_state_with_rate_origin(
        &self,
        ut_at: &dyn Fn() -> Result<f64>,
        tt: f64,
        body: i32,
        flags: i32,
        state: PositionResult,
        offset_j2000: Option<([f64; 3], [f64; 3])>,
        rate_origin: ModelRateOrigin,
    ) -> Result<(PositionResult, i32)> {
        self.local_with_time_state_using_ut(
            ut_at,
            tt,
            body,
            flags,
            Some(state),
            offset_j2000,
            !tt.is_finite(),
            rate_origin,
        )
    }

    /// Adapt an already acquired UT1 JD to the local state reducer.
    /// The constant callback preserves scalar remote-entry clock semantics.
    fn local_with_state(
        &self,
        ut: f64,
        tt: f64,
        body: i32,
        flags: i32,
        supplied: Option<PositionResult>,
        hypo_topo: Option<([f64; 3], [f64; 3])>,
    ) -> Result<(PositionResult, i32)> {
        self.local_with_state_using_ut(&|| Ok(ut), tt, body, flags, supplied, hypo_topo)
    }

    /// Bind local numerical providers without acquiring an unused UT1 clock.
    /// State/frame/nutation inputs and results retain degrees, AU and daily
    /// rates; provider errors occur only when their owning reduction uses them.
    fn local_with_state_using_ut(
        &self,
        ut_at: &dyn Fn() -> Result<f64>,
        tt: f64,
        body: i32,
        flags: i32,
        supplied: Option<PositionResult>,
        hypo_topo: Option<([f64; 3], [f64; 3])>,
    ) -> Result<(PositionResult, i32)> {
        self.local_with_time_state_using_ut(
            ut_at,
            tt,
            body,
            flags,
            supplied,
            hypo_topo,
            false,
            ModelRateOrigin::LegacyStencil,
        )
    }

    /// Bind the local reducer to checked or internally acquired TT nutation.
    /// All callbacks remain local to this synchronous reduction; no clock is
    /// acquired until its consuming stage. Internal nutation series propagate
    /// IEEE arithmetic after the calculation entry validates the original epoch.
    /// `tt` is scalar TT JD; nutation longitude `dpsi` is radians and its
    /// longitude rate is degrees/day. Arguments and shaping follow
    /// [`HorizonsLocal::local_with_state_using_ut`], returning six requested
    /// position/rate components plus the integer returned flags. Consumed
    /// synchronous clock, nutation, sidereal, frame and model errors propagate;
    /// no callback or borrowed state survives this call.
    #[allow(clippy::too_many_arguments)]
    fn local_with_time_state_using_ut(
        &self,
        ut_at: &dyn Fn() -> Result<f64>,
        tt: f64,
        body: i32,
        flags: i32,
        supplied: Option<PositionResult>,
        hypo_topo: Option<([f64; 3], [f64; 3])>,
        source_time: bool,
        rate_origin: ModelRateOrigin,
    ) -> Result<(PositionResult, i32)> {
        // Local model nutation uses the complete IAU 2006/2000A series.
        let dpsi = |t| -> Result<f64> {
            if source_time {
                Ok(crate::erfa::source_nut06a(t).0)
            } else {
                Ok(crate::erfa::nut06a(t)?.0)
            }
        };
        let rate = |t: f64, h: f64| -> Result<f64> {
            if rate_origin == ModelRateOrigin::LegacyStencil {
                Ok((dpsi(t + h)? - dpsi(t - h)?).to_degrees() / (2.0 * h))
            } else if source_time {
                Ok(crate::erfa::source_nut06a_longitude_rate(t).to_degrees())
            } else {
                crate::erfa::nut06a_longitude_rate(t).map(f64::to_degrees)
            }
        };
        let sid = |s, u, f| local_sidereal(s, u, f, self.aya);
        HorizonsLocal {
            coverage: self.coverage,
            delta_t_at: self.delta_t,
            dpsi_at: &dpsi,
            nutation_rate_at: &rate,
            sidereal: &sid,
            convert: &local_convert,
            precess_state: &local_precess_state,
        }
        .local_with_state_using_ut(ut_at, tt, body, flags, supplied, hypo_topo)
    }
}

/// Geocentric parallax for a tropical hypothetical place.
///
/// `state` is the nutated mean-ecliptic-of-date place in degrees/AU with
/// per-day rates; `tt` is scalar TT JD; `offset_j2000` is the observer
/// geocentric offset (topocentre minus geocentre) and its rate in
/// J2000-ecliptic AU and AU/day. Returns the topocentric place in the same
/// representation. Requests without `FLG_TOPOCTR` and heliocentric or
/// barycentric requests return `state` unchanged: parallax is undefined off
/// the geocentre. A zero topocentric radius returns the input unchanged.
/// Only finite inputs are accepted.
///
/// The J2000 offset is rotated to of-date ecliptic axes by precessing its
/// direction from J2000 to `tt` at preserved magnitude (exact for the
/// direction; the offset is about one Earth radius, so any residual is far
/// below the served model error). The body state converts through the
/// spherical/Cartesian relations with zero-distance and polar guards.
/// Through [`HorizonsNumerics::hypothetical_state_using_ut`] this helper
/// runs on the tropical of-date place, after nutation and before any
/// sidereal rotation or equatorial/J2000 conversion.
pub(crate) fn apply_hypothetical_topocentric_of_date(
    state: PositionResult,
    tt: f64,
    flags: i32,
    offset_j2000: ([f64; 3], [f64; 3]),
) -> Result<PositionResult> {
    if flags & FLG_TOPOCTR == 0 || flags & (FLG_HELCTR | FLG_BARYCTR) != 0 {
        return Ok(state);
    }
    for x in state.into_iter().chain([tt]) {
        if !x.is_finite() {
            return Err(Error::NonFiniteArgument { value: x });
        }
    }
    for v in offset_j2000.0.into_iter().chain(offset_j2000.1) {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    let (o, ov) = (
        rotate_offset_to_of_date(offset_j2000.0, tt)?,
        rotate_offset_to_of_date(offset_j2000.1, tt)?,
    );
    let (ox, oy, oz) = (o[0], o[1], o[2]);
    let (ovx, ovy, ovz) = (ov[0], ov[1], ov[2]);
    let (lon, lat, dist, dlon, dlat, ddist) =
        (state[0], state[1], state[2], state[3], state[4], state[5]);
    let (lr, br) = (lon.to_radians(), lat.to_radians());
    let (dlr, dbr) = (dlon.to_radians(), dlat.to_radians());
    let (cl, sl) = (lr.cos(), lr.sin());
    let (cb, sb) = (br.cos(), br.sin());
    // Geocentric body state in of-date ecliptic Cartesian axes.
    let bx = dist * cb * cl;
    let by = dist * cb * sl;
    let bz = dist * sb;
    let bvx = ddist * cb * cl - dist * sb * dbr * cl - dist * cb * sl * dlr;
    let bvy = ddist * cb * sl - dist * sb * dbr * sl + dist * cb * cl * dlr;
    let bvz = ddist * sb + dist * cb * dbr;
    // Topocentric state: geocentric minus observer offset, position and rate.
    let (tx, ty, tz) = (bx - ox, by - oy, bz - oz);
    let (tvx, tvy, tvz) = (bvx - ovx, bvy - ovy, bvz - ovz);
    let r_xy2 = tx * tx + ty * ty;
    let r = (r_xy2 + tz * tz).sqrt();
    if r == 0.0 {
        return Ok(state);
    }
    let r_xy = r_xy2.sqrt();
    let lon_t = ty.atan2(tx).to_degrees().rem_euclid(360.0);
    let lat_t = (tz / r).clamp(-1.0, 1.0).asin().to_degrees();
    let ddist_t = (tx * tvx + ty * tvy + tz * tvz) / r;
    let dlon_t = if r_xy2 > 0.0 {
        ((tx * tvy - ty * tvx) / r_xy2).to_degrees()
    } else {
        0.0
    };
    let dlat_t = if r_xy > 0.0 {
        ((tvz - (tz / r) * ddist_t) / r_xy).to_degrees()
    } else {
        0.0
    };
    Ok([lon_t, lat_t, r, dlon_t, dlat_t, ddist_t])
}

/// Rotate a J2000-ecliptic vector to the ecliptic of date at preserved
/// magnitude: the direction precesses from J2000 to `tt` while the length
/// is retained. Angles are degrees; `tt` is scalar TT JD. The zero vector
/// maps to itself. Non-finite components refuse.
fn rotate_offset_to_of_date(v: [f64; 3], tt: f64) -> Result<[f64; 3]> {
    for x in v.into_iter().chain([tt]) {
        if !x.is_finite() {
            return Err(Error::NonFiniteArgument { value: x });
        }
    }
    let mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if mag == 0.0 {
        return Ok([0.0; 3]);
    }
    let lon_v = v[1].atan2(v[0]).to_degrees();
    let lat_v = (v[2] / mag).clamp(-1.0, 1.0).asin().to_degrees();
    let (lon_d, lat_d) =
        crate::hypothetical::precess_ecliptic(lon_v, lat_v, crate::hypothetical::J2000_JD, tt)?;
    let (lr, br) = (lon_d.to_radians(), lat_d.to_radians());
    let cl = br.cos();
    Ok([mag * cl * lr.cos(), mag * cl * lr.sin(), mag * br.sin()])
}
/// true obliquity (radians), and precession-only matrix, respectively.
/// The matrices map ICRF into the true and mean equator of date.
pub(crate) type HorizonsFrame = ([[f64; 3]; 3], f64, f64, [[f64; 3]; 3]);

/// Bound numerical providers for Horizons output reduction.
///
/// `frame_at` supplies the true and mean equatorial date frames at scalar TT.
/// `ayanamsha_at` supplies the configured, flag-aware ayanamsha in degrees at UT.
/// Provider errors, including unavailable dynamic anchors, propagate unchanged.
pub(crate) struct HorizonsOutput<'a> {
    pub(crate) frame_at: &'a dyn Fn(f64) -> Result<HorizonsFrame>,
    pub(crate) ayanamsha_at: &'a dyn Fn(f64, i32) -> Result<f64>,
    /// Whether the snapshotted base mode uses an apparent stellar anchor.
    pub(crate) star_based: bool,
}

impl HorizonsOutput<'_> {
    /// Rotate an AU or AU/day vector at scalar TT. J2000 uses a fixed biased
    /// frame; SIDEREAL+EQUATORIAL uses precession only, without ayanamsha.
    fn rotate(&self, v: [f64; 3], tt: f64, flags: i32) -> Result<[f64; 3]> {
        if flags & FLG_J2000 != 0 {
            return if flags & FLG_EQUATORIAL != 0 {
                rotate_icrs_to_j2000_mean_equator(v)
            } else {
                rotate_icrs_to_ecliptic_j2000(v)
            };
        }
        let (pn, _, eps, p) = (self.frame_at)(tt)?;
        if flags & FLG_EQUATORIAL != 0 {
            mat3_vec3(if flags & FLG_SIDEREAL != 0 { p } else { pn }, v)
        } else {
            rotate_equatorial_to_ecliptic(mat3_vec3(pn, v)?, eps)
        }
    }

    /// Convert corrected ICRF position (AU) and velocity (AU/day) into the
    /// requested six-component representation, echoing flags unchanged.
    ///
    /// Moving-frame rates use linear state extrapolation through TT±0.0005 d.
    /// Stellar ayanamsha drift uses UT±1 second, with wrapped differences;
    /// other modes subtract polynomial precession and true-frame nutation
    /// drift. XYZ takes precedence over radians. Finite scalar epochs and
    /// vectors are required; frame/anchor acquisition failures propagate.
    /// Unsupported flags and center refusals belong to the dispatch caller.
    pub(crate) fn output(
        &self,
        position: [f64; 3],
        velocity: [f64; 3],
        tt: f64,
        ut: f64,
        flags: i32,
    ) -> Result<(PositionResult, i32)> {
        for x in position.into_iter().chain(velocity).chain([tt, ut]) {
            if !x.is_finite() {
                return Err(Error::NonFiniteArgument { value: x });
            }
        }
        let pos = self.rotate(position, tt, flags)?;
        let mut vel = self.rotate(velocity, tt, flags)?;
        let s = cartesian_to_spherical_deg(pos)?;
        let mut out = [s[0], s[1], s[2], 0.0, 0.0, 0.0];
        if flags & FLG_SPEED != 0 {
            if flags & FLG_J2000 != 0 {
                let rxy2 = pos[0] * pos[0] + pos[1] * pos[1];
                let rxy = rxy2.sqrt();
                let r2 = rxy2 + pos[2] * pos[2];
                let r = r2.sqrt();
                if r == 0.0 || rxy == 0.0 {
                    out[5] = if r == 0.0 {
                        (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]).sqrt()
                    } else {
                        0.0
                    };
                } else {
                    out[3] = ((pos[0] * vel[1] - pos[1] * vel[0]) / rxy2).to_degrees();
                    out[4] = ((vel[2] * rxy2 - pos[2] * (pos[0] * vel[0] + pos[1] * vel[1]))
                        / (r2 * rxy))
                        .to_degrees();
                    out[5] = (pos[0] * vel[0] + pos[1] * vel[1] + pos[2] * vel[2]) / r;
                }
            } else {
                let h = VEL_H_DAYS;
                let minus = self.rotate(
                    std::array::from_fn(|i| position[i] - velocity[i] * h),
                    tt - h,
                    flags,
                )?;
                let plus = self.rotate(
                    std::array::from_fn(|i| position[i] + velocity[i] * h),
                    tt + h,
                    flags,
                )?;
                let sm = cartesian_to_spherical_deg(minus)?;
                let sp = cartesian_to_spherical_deg(plus)?;
                let inv = 1.0 / (2.0 * h);
                let mut dl = (sp[0] - sm[0]).rem_euclid(360.0);
                if dl > 180.0 {
                    dl -= 360.0;
                }
                out[3] = dl * inv;
                out[4] = (sp[1] - sm[1]) * inv;
                out[5] = (sp[2] - sm[2]) * inv;
                if flags & FLG_XYZ != 0 {
                    vel = std::array::from_fn(|i| (plus[i] - minus[i]) * inv);
                }
            }
        } else {
            vel = [0.0; 3];
        }
        if flags & FLG_SIDEREAL != 0 && flags & FLG_EQUATORIAL == 0 {
            out[0] = (out[0] - (self.ayanamsha_at)(ut, flags)?).rem_euclid(360.0);
            if flags & FLG_SPEED != 0 {
                if self.star_based {
                    let dt = 1.0 / 86400.0;
                    let delta = ((self.ayanamsha_at)(ut + dt, flags)?
                        - (self.ayanamsha_at)(ut - dt, flags)?
                        + 180.0)
                        .rem_euclid(360.0)
                        - 180.0;
                    out[3] -= delta / (2.0 * dt);
                } else {
                    out[3] -= general_precession_rate_deg_day(tt)?;
                    if flags & FLG_J2000 == 0 {
                        out[3] -= ((self.frame_at)(tt + VEL_H_DAYS)?.1
                            - (self.frame_at)(tt - VEL_H_DAYS)?.1)
                            .to_degrees()
                            / (2.0 * VEL_H_DAYS);
                    }
                }
            }
            return Ok((apply_output_flags(out, flags), flags));
        }
        if flags & FLG_XYZ != 0 {
            return Ok(([pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]], flags));
        }
        Ok((apply_output_flags(out, flags), flags))
    }
}

/// Local-model tail providers. All states are spherical degrees/AU with
/// per-day rates. Conversion implements `_maybe_equatorial_convert`, with
/// `already_j2000` preventing a second hypothetical-body precession.
pub(crate) struct HorizonsLocal<'a> {
    pub(crate) coverage: crate::lunar::LunarCoverage,
    /// Configured `deltat(UT)` in days, used for the two UT inversion
    /// iterations when differentiating analytical output at exact TT samples.
    pub(crate) delta_t_at: &'a dyn Fn(f64) -> Result<f64>,
    pub(crate) dpsi_at: &'a dyn Fn(f64) -> Result<f64>,
    pub(crate) nutation_rate_at: &'a dyn Fn(f64, f64) -> Result<f64>,
    pub(crate) sidereal: &'a dyn Fn(PositionResult, f64, i32) -> Result<PositionResult>,
    pub(crate) convert: &'a dyn Fn(PositionResult, f64, i32, bool) -> Result<PositionResult>,
    /// `_precess_ecliptic_state`, including the moving-frame derivative.
    pub(crate) precess_state: &'a dyn Fn(PositionResult, f64, bool) -> Result<PositionResult>,
}

impl HorizonsLocal<'_> {
    /// Reduce a local mean-node/apogee or hypothetical native-center model.
    /// TT belongs to the entry: the analytical branch supplies configured
    /// Delta-T TT, while supplied hypothetical states retain their independently
    /// acquired model TT. This reducer does not replace that clock.
    /// Center/refusal gates belong to dispatch.
    /// Analytical SPEED3 differentiates the final representation at TT±h,
    /// solving UT with two Delta-T iterations. Hypothetical
    /// output does not acquire this analytical-only SPEED3 stencil.
    /// Missing conversion/sidereal providers and model-domain errors propagate.
    pub(crate) fn local(
        &self,
        ut: f64,
        tt: f64,
        body: i32,
        flags: i32,
    ) -> Result<(PositionResult, i32)> {
        self.local_with_state(ut, tt, body, flags, None, None)
    }

    /// Adapt an explicit scalar UT1 JD to the local model reduction.
    /// Lazy consumption preserves the existing supplied-entry semantics.
    fn local_with_state(
        &self,
        ut: f64,
        tt: f64,
        body: i32,
        flags: i32,
        supplied: Option<PositionResult>,
        hypo_topo: Option<([f64; 3], [f64; 3])>,
    ) -> Result<(PositionResult, i32)> {
        self.local_with_state_using_ut(&|| Ok(ut), tt, body, flags, supplied, hypo_topo)
    }

    /// Reduce a local state, acquiring scalar UT1 only for its owning tail.
    /// Ecliptic SID consumes the clock after frame/nutation/parallax; mean
    /// SPEED3 consumes it for its configured Delta-T inversion. Other routes
    /// preserve their TT-only domains and error order. Units follow `local`.
    fn local_with_state_using_ut(
        &self,
        ut_at: &dyn Fn() -> Result<f64>,
        tt: f64,
        body: i32,
        flags: i32,
        supplied: Option<PositionResult>,
        hypo_topo: Option<([f64; 3], [f64; 3])>,
    ) -> Result<(PositionResult, i32)> {
        let mean = body == MEAN_NODE || body == MEAN_APOG;
        if mean && flags & (FLG_HELCTR | FLG_BARYCTR) != 0 {
            return Ok(([0.0; 6], flags));
        }
        let mut state;
        if mean {
            let speed3 = flags & FLG_SPEED != 0 && flags & (FLG_SPEED3 | FLG_TOPOCTR) != 0;
            let s = if body == MEAN_NODE {
                crate::lunar::calc_mean_lunar_node_state(tt, self.coverage, speed3)?
            } else {
                crate::lunar::calc_mean_lilith_state(tt, self.coverage, speed3)?
            };
            let s = crate::lunar::shape_mean_point_state(
                s,
                tt,
                flags,
                self.dpsi_at,
                self.nutation_rate_at,
            )?;
            state = [s.0, s.1, s.2, s.3, s.4, s.5];
            if flags & FLG_SIDEREAL != 0 && flags & FLG_EQUATORIAL == 0 && flags & FLG_J2000 != 0 {
                state = (self.convert)(state, tt, flags, false)?;
                state = (self.sidereal)(state, ut_at()?, flags)?;
            } else {
                if flags & FLG_SIDEREAL != 0 && flags & FLG_EQUATORIAL == 0 {
                    state = (self.sidereal)(state, ut_at()?, flags)?;
                }
                state = (self.convert)(state, tt, flags, false)?;
            }
        } else {
            state = if let Some(s) = supplied {
                s
            } else {
                let s = crate::hypothetical::calc_hypothetical_position(body, tt)?;
                [s.lon_deg, s.lat_deg, s.dist_au, s.dlon, s.dlat, s.ddist]
            };
            if flags & FLG_SPEED == 0 {
                state[3..].fill(0.0);
            }
            let of_date = matches!(body, VULCAN | WHITE_MOON | PROSERPINA | WALDEMATH);
            // Apply parallax in the ecliptic-of-date frame. J2000 analytical models
            // first rotate to date; any requested J2000 output rotates back only after
            // the observer displacement has been subtracted.
            let topo_active = hypo_topo.is_some() && flags & FLG_TOPOCTR != 0;
            let to_of_date_first = (topo_active && !of_date)
                || (!topo_active && (of_date != (flags & FLG_J2000 == 0)));
            if to_of_date_first {
                if flags & FLG_SPEED != 0 {
                    state = (self.precess_state)(state, tt, of_date)?;
                } else {
                    let (from, to) = if of_date {
                        (tt, 2451545.0)
                    } else {
                        (2451545.0, tt)
                    };
                    (state[0], state[1]) =
                        crate::hypothetical::precess_ecliptic(state[0], state[1], from, to)?;
                }
            }
            let no_nut = flags & (FLG_NONUT | FLG_J2000) != 0
                || flags & (FLG_SIDEREAL | FLG_EQUATORIAL) == (FLG_SIDEREAL | FLG_EQUATORIAL);
            if !no_nut {
                state[0] = (state[0] + (self.dpsi_at)(tt)?.to_degrees()).rem_euclid(360.0);
                if flags & FLG_SPEED != 0 {
                    state[3] += (self.nutation_rate_at)(tt, 0.5)?;
                }
            }
            if let Some(offset) = hypo_topo {
                state = apply_hypothetical_topocentric_of_date(state, tt, flags, offset)?;
            }
            if flags & FLG_SIDEREAL != 0 && flags & FLG_EQUATORIAL == 0 {
                state = (self.sidereal)(state, ut_at()?, flags)?;
            }
            state = (self.convert)(state, tt, flags, hypo_topo.is_none())?;
        }
        let mut output = apply_output_flags(state, flags);
        if mean && flags & FLG_SPEED3 != 0 {
            let ut = ut_at()?;
            let sample_flags = flags & !(FLG_SPEED | FLG_SPEED3);
            let sample = |target_tt: f64| -> Result<PositionResult> {
                let mut sample_ut = target_tt - (self.delta_t_at)(ut)?;
                sample_ut = target_tt - (self.delta_t_at)(sample_ut)?;
                // The recursive calculation recomputes TT from the solved UT;
                // retaining target_tt here would erase its observable rounding.
                let sample_tt = sample_ut + (self.delta_t_at)(sample_ut)?;
                Ok(self.local(sample_ut, sample_tt, body, sample_flags)?.0)
            };
            output = apply_central_speed3_stencil(
                output,
                sample(tt - SPEED3_STEP_DAYS)?,
                sample(tt + SPEED3_STEP_DAYS)?,
                flags,
            );
        }
        Ok((output, flags))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn native_hypothetical_tail_acquires_ut_only_at_ecliptic_sidereal_reduction() {
        use std::cell::Cell;
        let calls = Cell::new(0);
        let ut_at = || {
            calls.set(calls.get() + 1);
            Err(Error::UnknownBody { body_id: 999 })
        };
        let delta = |_| panic!("hypothetical tail consumes no model Delta-T");
        let aya = |_, _| panic!("failed UT acquisition precedes ayanamsha");
        let numerics =
            HorizonsNumerics::new(&delta, &aya, 0, crate::lunar::LunarCoverage::Unenforced);
        let state = [120.0, 2.0, 1.0, 0.5, 0.25, 0.125];
        for flags in [0, FLG_EQUATORIAL, FLG_SIDEREAL | FLG_EQUATORIAL] {
            assert!(
                numerics
                    .hypothetical_state_using_ut(&ut_at, 2451545.0, VULCAN, flags, state, None)
                    .is_ok()
            );
            assert_eq!(calls.get(), 0);
        }
        let flags = FLG_SIDEREAL | FLG_TOPOCTR;
        let invalid_offset = Some(([f64::NAN, 0.0, 0.0], [0.0; 3]));
        assert!(matches!(
            numerics.hypothetical_state_using_ut(
                &ut_at,
                2451545.0,
                VULCAN,
                flags,
                state,
                invalid_offset
            ),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert_eq!(calls.get(), 0);
        assert_eq!(
            numerics.hypothetical_state_using_ut(
                &ut_at,
                2451545.0,
                VULCAN,
                FLG_SIDEREAL,
                state,
                None
            ),
            Err(Error::UnknownBody { body_id: 999 })
        );
        assert_eq!(calls.get(), 1);
    }

    #[test]
    fn fixed_frame_without_speed_and_xyz_priority() {
        let unavailable =
            |_: f64| -> Result<HorizonsFrame> { panic!("fixed frame must not acquire date data") };
        let anchor = |_: f64, _: i32| -> Result<f64> { panic!("tropical must not acquire anchor") };
        let reducer = HorizonsOutput {
            frame_at: &unavailable,
            ayanamsha_at: &anchor,
            star_based: false,
        };
        let flags = FLG_J2000 | FLG_XYZ | FLG_RADIANS;
        let (out, echoed) = reducer
            .output([1.0, 0.2, 0.3], [1.0; 3], 2451545.0, 2451545.0, flags)
            .unwrap();
        assert_eq!(echoed, flags);
        assert_eq!(&out[3..], &[0.0; 3]);
        assert_eq!(
            &out[..3],
            &rotate_icrs_to_ecliptic_j2000([1.0, 0.2, 0.3]).unwrap()
        );
    }

    #[test]
    fn rejects_nonfinite_before_provider_calls() {
        let frame = |_: f64| -> Result<HorizonsFrame> { panic!("invalid input") };
        let anchor = |_: f64, _: i32| -> Result<f64> { panic!("invalid input") };
        let reducer = HorizonsOutput {
            frame_at: &frame,
            ayanamsha_at: &anchor,
            star_based: false,
        };
        assert!(matches!(
            reducer.output([f64::NAN, 0.0, 0.0], [0.0; 3], 2451545.0, 2451545.0, 0),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn hypothetical_topocentric_is_inert_off_the_geocentre() {
        let state = [120.0, 5.0, 0.01, 1.0, 0.1, 0.001];
        let offset = ([4e-5, 0.0, 0.0], [0.0, 2e-4, 0.0]);
        // No topocentric flag, heliocentric and barycentric requests all
        // pass the state through bitwise unchanged.
        for flags in [0, FLG_HELCTR, FLG_BARYCTR, FLG_HELCTR | FLG_TOPOCTR] {
            assert_eq!(
                apply_hypothetical_topocentric_of_date(state, 2460000.5, flags, offset).unwrap(),
                state
            );
        }
        // A zero observer offset runs the spherical roundtrip, which is
        // near-identity (one-ulp trigonometric wobble, from ordinary trigonometric rounding),
        // not a bitwise copy.
        let roundtrip = apply_hypothetical_topocentric_of_date(
            state,
            2460000.5,
            FLG_TOPOCTR,
            ([0.0; 3], [0.0; 3]),
        )
        .unwrap();
        for i in 0..6 {
            assert!(
                (roundtrip[i] - state[i]).abs() < 1e-12,
                "slot {i}: {} vs {}",
                roundtrip[i],
                state[i]
            );
        }
    }

    #[test]
    fn hypothetical_topocentric_shifts_by_about_offset_over_distance() {
        // Body on the of-date +x axis; observer offset near +y (up to the
        // J2000-to-date direction precession). The topocentric longitude
        // shifts by about |offset|/dist radians negative, latitude and the
        // distance barely move, and zero input rates stay zero.
        let dist = 0.01;
        let off = 4e-5;
        let state = [0.0, 0.0, dist, 0.0, 0.0, 0.0];
        let out = apply_hypothetical_topocentric_of_date(
            state,
            2460000.5,
            FLG_TOPOCTR | FLG_SPEED,
            ([0.0, off, 0.0], [0.0; 3]),
        )
        .unwrap();
        let naive_lon = (-off / dist).atan().to_degrees().rem_euclid(360.0);
        let mut dlon = (out[0] - naive_lon).rem_euclid(360.0);
        if dlon > 180.0 {
            dlon -= 360.0;
        }
        assert!(dlon.abs() < 0.01, "lon {} vs naive {naive_lon}", out[0]);
        assert!(out[1].abs() < 0.01, "lat {}", out[1]);
        assert!((out[2] - (dist * dist + off * off).sqrt()).abs() < 1e-6);
        assert_eq!(&out[3..], &[0.0; 3]);
    }

    #[test]
    fn hypothetical_topocentric_zero_radius_returns_input() {
        // Zero body distance with a zero observer offset: the topocentric
        // radius is exactly zero, so the geocentric input returns rather
        // than dividing by zero.
        let state = [30.0, -10.0, 0.0, 0.0, 0.0, 0.0];
        let out = apply_hypothetical_topocentric_of_date(
            state,
            2460000.5,
            FLG_TOPOCTR,
            ([0.0; 3], [0.0; 3]),
        )
        .unwrap();
        assert_eq!(out, state);
    }

    #[test]
    fn hypothetical_topocentric_refuses_nonfinite() {
        let state = [0.0, 0.0, 0.01, 0.0, 0.0, 0.0];
        let offset = ([1e-5, 0.0, 0.0], [0.0; 3]);
        assert!(matches!(
            apply_hypothetical_topocentric_of_date(
                [f64::NAN, 0.0, 0.01, 0.0, 0.0, 0.0],
                2460000.5,
                FLG_TOPOCTR,
                offset
            ),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            apply_hypothetical_topocentric_of_date(state, f64::INFINITY, FLG_TOPOCTR, offset),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            apply_hypothetical_topocentric_of_date(
                state,
                2460000.5,
                FLG_TOPOCTR,
                ([f64::NAN, 0.0, 0.0], [0.0; 3])
            ),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    fn hypo_tail_for_test() -> (
        impl Fn(f64) -> Result<f64>,
        impl Fn(f64, i32) -> Result<f64>,
    ) {
        // Zero Delta-T and a constant 10-degree ayanamsha: the sidereal
        // tail subtracts exactly 10 degrees, so the parallax/sidereal
        // order is observable in the words.
        (|_| Ok(0.0), |_, _| Ok(10.0))
    }

    #[test]
    fn hypothetical_topo_tail_matches_legacy_post_tail_tropical_words() {
        // Tropical topocentric requests served through the in-tail stage
        // are bitwise identical to the previous post-tail application: the
        // conversion tail is the identity without J2000/equatorial bits.
        // HARRINGTON is natively J2000, so this also pins that its early
        // J2000-to-of-date precession still runs under parallax.
        let (delta_t, _) = hypo_tail_for_test();
        let aya = |_: f64, _: i32| -> Result<f64> { Ok(0.0) };
        let numerics =
            HorizonsNumerics::new(&delta_t, &aya, 0, crate::lunar::LunarCoverage::Unenforced);
        let tt = 2460000.5;
        let state = [120.0, 5.0, 0.5, 0.01, 0.002, 0.0001];
        let offset = ([4e-5, 1e-5, 2e-5], [1e-7, -2e-7, 3e-7]);
        for body in [VULCAN, HARRINGTON] {
            for flags in [FLG_TOPOCTR, FLG_TOPOCTR | FLG_SPEED] {
                let (legacy, legacy_flag) = numerics
                    .hypothetical_state(tt, tt, body, flags, state)
                    .unwrap();
                let legacy =
                    apply_hypothetical_topocentric_of_date(legacy, tt, flags, offset).unwrap();
                let (tail, tail_flag) = numerics
                    .hypothetical_state_with_topo(tt, tt, body, flags, state, offset)
                    .unwrap();
                assert_eq!(legacy_flag, tail_flag);
                assert_eq!(legacy, tail, "body {body} flags {flags:#x}");
            }
        }
    }

    #[test]
    fn hypothetical_topo_tail_applies_parallax_before_sidereal() {
        // With a constant 10-degree ayanamsha the sidereal longitude is the
        // topocentric longitude minus 10 degrees: parallax first, rotation
        // second, in the tropical frame. Applying the parallax after
        // the rotation would shift the words by the rotated offset.
        let (delta_t, aya) = hypo_tail_for_test();
        let numerics =
            HorizonsNumerics::new(&delta_t, &aya, 0, crate::lunar::LunarCoverage::Unenforced);
        let tt = 2460000.5;
        let state = [120.0, 5.0, 0.5, 0.0, 0.0, 0.0];
        let offset = ([4e-5, 1e-5, 2e-5], [0.0; 3]);
        let flags = FLG_TOPOCTR | FLG_SIDEREAL;
        let (out, _) = numerics
            .hypothetical_state_with_topo(tt, tt, VULCAN, flags, state, offset)
            .unwrap();
        // The tail nutates before the parallax, before the output-frame conversion:
        // rebuild the expectation with the same primitive order.
        let dpsi_deg = crate::erfa::nut06a(tt).unwrap().0.to_degrees();
        let nutated = [
            (state[0] + dpsi_deg).rem_euclid(360.0),
            state[1],
            state[2],
            state[3],
            state[4],
            state[5],
        ];
        let topo =
            apply_hypothetical_topocentric_of_date(nutated, tt, FLG_TOPOCTR, offset).unwrap();
        let expected_lon = (topo[0] - 10.0).rem_euclid(360.0);
        assert!(
            (out[0] - expected_lon).abs() < 1e-9,
            "lon {} vs parallax-first {expected_lon}",
            out[0]
        );
        // The parallax is real: the sidereal word differs from the
        // parallax-free sidereal word by arcseconds, not microarcseconds.
        let (geo, _) = numerics
            .hypothetical_state(tt, tt, VULCAN, flags, state)
            .unwrap();
        let mut shift = (out[0] - geo[0]).rem_euclid(360.0);
        if shift > 180.0 {
            shift -= 360.0;
        }
        assert!(
            shift.abs() > 1e-6 && shift.abs() < 1.0,
            "sidereal parallax shift {shift}"
        );
    }

    #[test]
    fn hypothetical_topo_tail_precesses_topocentric_place_to_j2000() {
        // J2000 output rotates the already-topocentric place after applying
        // parallax in the frame of date.
        let (delta_t, _) = hypo_tail_for_test();
        let aya = |_: f64, _: i32| -> Result<f64> { Ok(0.0) };
        let numerics =
            HorizonsNumerics::new(&delta_t, &aya, 0, crate::lunar::LunarCoverage::Unenforced);
        let tt = 2460000.5;
        let state = [120.0, 5.0, 0.5, 0.0, 0.0, 0.0];
        let offset = ([4e-5, 1e-5, 2e-5], [0.0; 3]);
        let flags = FLG_TOPOCTR | FLG_J2000;
        let (out, _) = numerics
            .hypothetical_state_with_topo(tt, tt, VULCAN, flags, state, offset)
            .unwrap();
        let topo = apply_hypothetical_topocentric_of_date(state, tt, FLG_TOPOCTR, offset).unwrap();
        let (lon_j2000, lat_j2000) = crate::hypothetical::precess_ecliptic(
            topo[0],
            topo[1],
            tt,
            crate::hypothetical::J2000_JD,
        )
        .unwrap();
        assert!(
            (out[0] - lon_j2000).abs() < 1e-9 && (out[1] - lat_j2000).abs() < 1e-9,
            "j2000 [{}, {}] vs late-precessed [{lon_j2000}, {lat_j2000}]",
            out[0],
            out[1]
        );
        assert!(
            (out[2] - topo[2]).abs() < 1e-12,
            "distance {} vs topocentric {}",
            out[2],
            topo[2]
        );
    }

    #[test]
    fn local_convert_uses_vondrak_of_date_obliquity() {
        // Recorded equatorial positions use the Vondrák long-term
        // pole-angle obliquity, distinct from the IAU short-term polynomial.
        let tt = 2440000.5;
        let state = [120.0, 5.0, 0.5, 0.01, 0.005, 0.0001];
        let expected_true = [
            123.35061438998025,
            25.038264569191384,
            0.5,
            0.011938635432020313,
            0.0026905860770085346,
            0.0001,
        ];
        let expected_mean = [
            123.35000805171084,
            25.036292191118033,
            0.5,
            0.011938369203255583,
            0.002690856555243848,
            0.0001,
        ];
        for (flags, expected) in [
            (FLG_SPEED | FLG_EQUATORIAL, expected_true),
            (FLG_EQUATORIAL, expected_true),
            (FLG_SPEED | FLG_EQUATORIAL | FLG_SIDEREAL, expected_mean),
            (FLG_EQUATORIAL | FLG_SIDEREAL, expected_mean),
            (FLG_SPEED | FLG_EQUATORIAL | FLG_NONUT, expected_mean),
        ] {
            let out = local_convert(state, tt, flags, false).unwrap();
            let count = if flags & FLG_SPEED != 0 { 3 } else { 6 };
            for (got, want) in out.iter().zip(expected.iter()).take(count) {
                assert!(
                    (got - want).abs() <= 1e-12,
                    "flags {flags}: {got} vs source {want}"
                );
            }
        }
    }

    #[test]
    fn equatorial_rates_differentiate_the_actual_moving_frame() {
        for tt in [2451545.0, 2460310.5, 2400000.5] {
            for flags in [
                FLG_EQUATORIAL,
                FLG_EQUATORIAL | FLG_NONUT,
                FLG_EQUATORIAL | FLG_SIDEREAL,
                FLG_EQUATORIAL | FLG_J2000,
            ] {
                let state = [120.0, 15.0, 2.0, 0.01, -0.003, 0.0001];
                let output = local_convert(state, tt, flags | FLG_SPEED, true).unwrap();
                let sample = |offset: f64| {
                    local_convert(
                        [
                            state[0] + state[3] * offset,
                            state[1] + state[4] * offset,
                            state[2] + state[5] * offset,
                            0.0,
                            0.0,
                            0.0,
                        ],
                        tt + offset,
                        flags,
                        true,
                    )
                    .unwrap()
                };
                let h = 1.0 / 32.0;
                let (m, p, m2, p2, m3, p3) = (
                    sample(-h),
                    sample(h),
                    sample(-2.0 * h),
                    sample(2.0 * h),
                    sample(-3.0 * h),
                    sample(3.0 * h),
                );
                for i in 0..2 {
                    let derivative =
                        (45.0 * (p[i] - m[i]) - 9.0 * (p2[i] - m2[i]) + p3[i] - m3[i]) / (60.0 * h);
                    assert!(
                        (output[i + 3] - derivative).abs() <= derivative.abs() * 1e-6,
                        "tt={tt}, flags={flags}, slot={i}: {} vs {derivative}",
                        output[i + 3]
                    );
                }
                assert_eq!(output[5], state[5]);
            }
        }
    }

    #[test]
    fn hypothetical_topo_tail_applies_parallax_before_equatorial() {
        // Equatorial topocentric requests convert the topocentric
        // ecliptic place: the output equals the equatorial conversion of
        // the parallax-corrected state, never a parallax of equatorial
        // coordinates (which would mix frames).
        let (delta_t, _) = hypo_tail_for_test();
        let aya = |_: f64, _: i32| -> Result<f64> { Ok(0.0) };
        let numerics =
            HorizonsNumerics::new(&delta_t, &aya, 0, crate::lunar::LunarCoverage::Unenforced);
        let tt = 2460000.5;
        let state = [120.0, 5.0, 0.5, 0.0, 0.0, 0.0];
        let offset = ([4e-5, 1e-5, 2e-5], [0.0; 3]);
        let flags = FLG_TOPOCTR | FLG_EQUATORIAL;
        let (out, _) = numerics
            .hypothetical_state_with_topo(tt, tt, VULCAN, flags, state, offset)
            .unwrap();
        // The tail nutates before the parallax, before the output-frame conversion:
        // rebuild the expectation with the same primitive order.
        let dpsi_deg = crate::erfa::nut06a(tt).unwrap().0.to_degrees();
        let nutated = [
            (state[0] + dpsi_deg).rem_euclid(360.0),
            state[1],
            state[2],
            state[3],
            state[4],
            state[5],
        ];
        let topo =
            apply_hypothetical_topocentric_of_date(nutated, tt, FLG_TOPOCTR, offset).unwrap();
        let converted = local_convert(topo, tt, flags, false).unwrap();
        assert_eq!(out, converted);
        // And the parallax survives the conversion: the equatorial word
        // differs from the parallax-free equatorial word.
        let (geo, _) = numerics
            .hypothetical_state(tt, tt, VULCAN, flags, state)
            .unwrap();
        assert!(
            (out[0] - geo[0]).abs() > 1e-9 || (out[1] - geo[1]).abs() > 1e-9,
            "equatorial parallax invisible: {out:?} vs {geo:?}"
        );
    }
}

#[cfg(test)]
mod source_time_domain_tests {
    use super::*;

    #[test]
    fn bound_tail_uses_instantaneous_nutation_only_for_owned_derivatives() {
        let numerics = HorizonsNumerics::new(
            &|_| panic!("unused Delta-T"),
            &|_, _| panic!("unused ayanamsha"),
            0,
            crate::lunar::LunarCoverage::Unenforced,
        );
        let tt = 2451545.0;
        let state = [120.0, 2.0, 1.0, 0.25, 0.125, 0.0625];
        let legacy = numerics
            .hypothetical_source_time_state_using_ut(
                &|| panic!("unused UT"),
                tt,
                VULCAN,
                FLG_SPEED,
                state,
                None,
            )
            .unwrap();
        for origin in [
            ModelRateOrigin::LegacyStencil,
            ModelRateOrigin::OwnedDerivative,
        ] {
            let actual = numerics
                .hypothetical_source_time_state_with_rate_origin(
                    &|| panic!("unused UT"),
                    tt,
                    VULCAN,
                    FLG_SPEED,
                    state,
                    None,
                    origin,
                )
                .unwrap();
            assert_eq!(actual.1, FLG_SPEED);
            assert_eq!(actual.0[..3], legacy.0[..3]);
            assert_eq!(actual.0[4..], state[4..]);
            let rate = match origin {
                ModelRateOrigin::LegacyStencil => (crate::erfa::nut06a(tt + 0.5).unwrap().0
                    - crate::erfa::nut06a(tt - 0.5).unwrap().0)
                    .to_degrees(),
                ModelRateOrigin::OwnedDerivative => {
                    crate::erfa::nut06a_longitude_rate(tt).unwrap().to_degrees()
                }
            };
            assert_eq!(actual.0[3].to_bits(), (state[3] + rate).to_bits());
            if origin == ModelRateOrigin::LegacyStencil {
                assert_eq!(actual.0.map(f64::to_bits), legacy.0.map(f64::to_bits));
            } else {
                assert_ne!(actual.0[3].to_bits(), legacy.0[3].to_bits());
            }
        }
    }

    #[test]
    fn rate_origin_preserves_lazy_ut_error_order_and_ieee_source_time() {
        use std::cell::Cell;
        let calls = Cell::new(0);
        let numerics = HorizonsNumerics::new(
            &|_| panic!("unused Delta-T"),
            &|_, _| panic!("UT error precedes ayanamsha"),
            0,
            crate::lunar::LunarCoverage::Unenforced,
        );
        for origin in [
            ModelRateOrigin::LegacyStencil,
            ModelRateOrigin::OwnedDerivative,
        ] {
            let ut = || {
                calls.set(calls.get() + 1);
                Err(Error::MathDomain { value: 456.0 })
            };
            let output = numerics.hypothetical_source_time_state_with_rate_origin(
                &ut,
                2451545.0,
                VULCAN,
                FLG_SPEED | FLG_SIDEREAL,
                [120.0, 2.0, 1.0, 0.25, 0.125, 0.0625],
                None,
                origin,
            );
            assert!(matches!(output, Err(Error::MathDomain { value: 456.0 })));
            assert_eq!(calls.replace(0), 1);
            let nan = numerics
                .hypothetical_source_time_state_with_rate_origin(
                    &ut,
                    f64::NAN,
                    VULCAN,
                    FLG_SPEED,
                    [f64::NAN, 0.0, f64::NAN, 0.25, 0.0, 0.0],
                    None,
                    origin,
                )
                .unwrap();
            assert!(nan.0[0].is_nan() && nan.0[3].is_nan());
            assert_eq!(calls.get(), 0);
        }
    }

    #[test]
    fn source_time_speed_preserves_ieee_without_consuming_a_clock_callback() {
        let unused_delta = |_: f64| -> Result<f64> { panic!("unused Delta-T") };
        let unused_aya = |_: f64, _: i32| -> Result<f64> { panic!("unused sidereal") };
        let numerics = HorizonsNumerics::new(
            &unused_delta,
            &unused_aya,
            0,
            crate::lunar::LunarCoverage::Unenforced,
        );
        let (state, flags) = numerics
            .hypothetical_source_time_state_using_ut(
                &|| panic!("unused UT inverse"),
                f64::NAN,
                VULCAN,
                FLG_SPEED,
                [f64::NAN, 0.0, f64::NAN, 0.25, 0.0, 0.0],
                None,
            )
            .unwrap();
        assert_eq!(flags, FLG_SPEED);
        assert!(state[0].is_nan() && state[3].is_nan());
        assert_eq!(state[4].to_bits(), 0.0_f64.to_bits());
        assert_eq!(state[5].to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn supplied_nutation_rate_callback_keeps_values_order_and_errors() {
        use std::cell::RefCell;
        for fail in [false, true] {
            let calls = RefCell::new(Vec::new());
            let dpsi = |_: f64| {
                calls.borrow_mut().push("angle");
                Ok(0.001)
            };
            let rate = |_: f64, h: f64| {
                calls.borrow_mut().push("rate");
                assert_eq!(h.to_bits(), 0.5_f64.to_bits());
                if fail {
                    Err(Error::MathDomain { value: 123.0 })
                } else {
                    Ok(0.125)
                }
            };
            let convert = |s, _, _, _| {
                calls.borrow_mut().push("convert");
                Ok(s)
            };
            let local = HorizonsLocal {
                coverage: crate::lunar::LunarCoverage::Unenforced,
                delta_t_at: &|_| panic!("unused Delta-T"),
                dpsi_at: &dpsi,
                nutation_rate_at: &rate,
                sidereal: &|_, _, _| panic!("unused sidereal"),
                convert: &convert,
                precess_state: &|_, _, _| panic!("already of date"),
            };
            let state = [120.0, 5.0, 0.5, 0.25, 0.0, 0.0];
            let output = local.local_with_state_using_ut(
                &|| panic!("unused UT"),
                2451545.0,
                VULCAN,
                FLG_SPEED,
                Some(state),
                None,
            );
            if fail {
                assert!(matches!(output, Err(Error::MathDomain { value: 123.0 })));
                assert_eq!(*calls.borrow(), ["angle", "rate"]);
            } else {
                let output = output.unwrap();
                assert_eq!(
                    output.0[0].to_bits(),
                    (120.0 + 0.001_f64.to_degrees()).to_bits()
                );
                assert_eq!(output.0[3].to_bits(), (0.25 + 0.125_f64).to_bits());
                assert_eq!(*calls.borrow(), ["angle", "rate", "convert"]);
            }
        }
    }

    #[test]
    fn symbolic_source_time_tail_keeps_nan_and_positive_zero_without_ut_inverse() {
        let unused_delta = |_: f64| -> Result<f64> { panic!("no defining anchor used") };
        let unused_aya = |_: f64, _: i32| -> Result<f64> { panic!("tropical output") };
        let numerics = HorizonsNumerics::new(
            &unused_delta,
            &unused_aya,
            0,
            crate::lunar::LunarCoverage::Unenforced,
        );
        let state = [f64::NAN, 0.0, f64::NAN, 0.0, 0.0, 0.0];
        let (result, flags) = numerics
            .hypothetical_source_time_state_using_ut(
                &|| panic!("unused UT inverse"),
                f64::NAN,
                WHITE_MOON,
                FLG_SWIEPH,
                state,
                None,
            )
            .unwrap();
        assert_eq!(flags, FLG_SWIEPH);
        assert_eq!(
            result.map(f64::to_bits),
            [0x7ff8_0000_0000_0000, 0, 0x7ff8_0000_0000_0000, 0, 0, 0]
        );
    }
}
