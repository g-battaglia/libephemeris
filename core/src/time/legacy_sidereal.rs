//! Explicit boundary-calibrated sidereal utilities.
//!
//! The historical convention joins IAU 2006 inside 1850–2050 to a geometric
//! Earth-longitude construction outside that window. Its supplied boundary
//! offsets remain real corrections, in degrees; they are never discarded.
//! This construction has the unresolved long-term UT1 clock defect described
//! in issue #50. It is retained for callers that explicitly need its former
//! convention, rather than as an independently qualified celestial origin.
//!
//! Canonical mean time uses separate UT1 ERA and TT celestial-origin clocks.
//! This explicitly deprecated namespace preserves the old convention, real
//! boundary-offset effects and its original finite-date admission.

use super::{
    SIDEREAL_MODERN_END_JD, SIDEREAL_MODERN_START_JD, SiderealDeltaT, gmst_iau2006_deg,
    mean_sidereal_longterm_deg,
};
use crate::errors::{Error, Result};

/// Compute the legacy branch-joining correction, in degrees `[-180,180)`.
///
/// `boundary_jd` is a scalar UT1 Julian day and `delta_t_days` is TT-minus-UT1
/// there, in days. The correction is geometric long-term GMST minus IAU 2006
/// GMST at that epoch. It depends on the supplied Delta-T; it is not an origin
/// constant. Nonfinite arguments return [`Error::NonFiniteArgument`]; errors
/// from the geometric frame propagate. No remote physical accuracy is implied.
pub fn sidereal_continuity_offset_deg(boundary_jd: f64, delta_t_days: f64) -> Result<f64> {
    if !boundary_jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: boundary_jd });
    }
    if !delta_t_days.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: delta_t_days,
        });
    }
    let jd_tt = boundary_jd + delta_t_days;
    let diff = mean_sidereal_longterm_deg(boundary_jd, delta_t_days)?
        - gmst_iau2006_deg(boundary_jd, jd_tt)?.rem_euclid(360.0);
    Ok((diff + 180.0).rem_euclid(360.0) - 180.0)
}

/// Return the legacy 1850 and 2050 corrections, in degrees.
///
/// The inputs are the separate TT-minus-UT1 values at those boundaries, in
/// days. Evaluation follows lower then upper boundary and propagates the
/// first error from [`sidereal_continuity_offset_deg`]. This utility acquires
/// no Delta-T model or configuration and guarantees no remote accuracy.
pub fn sidereal_continuity_offsets(
    delta_t_lo_days: f64,
    delta_t_hi_days: f64,
) -> Result<(f64, f64)> {
    Ok((
        sidereal_continuity_offset_deg(SIDEREAL_MODERN_START_JD, delta_t_lo_days)?,
        sidereal_continuity_offset_deg(SIDEREAL_MODERN_END_JD, delta_t_hi_days)?,
    ))
}

/// Return legacy Greenwich mean sidereal time, in degrees `[0,360)`.
///
/// IAU 2006 is used strictly inside the 1850–2050 window. At and below the
/// lower boundary the geometric branch subtracts `dt.offset_lo_deg()`;
/// at and above the upper boundary it subtracts `dt.offset_hi_deg()`.
/// Offsets may be any finite values accepted by [`SiderealDeltaT::try_new`].
/// `jd_ut1` is a scalar UT1 Julian day; TT is obtained from its supplied
/// Delta-T in days. Nonfinite UT1 returns [`Error::NonFiniteArgument`]; frame
/// and arithmetic errors retain the underlying legacy behavior. In particular
/// this utility does not repair issue #50 or promise correct independent UT1
/// rotation outside the modern window. No new finite date limit is imposed.
pub fn mean_sidereal_time_deg(jd_ut1: f64, dt: &SiderealDeltaT) -> Result<f64> {
    if !jd_ut1.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut1 });
    }
    if SIDEREAL_MODERN_START_JD < jd_ut1 && jd_ut1 < SIDEREAL_MODERN_END_JD {
        let jd_tt = jd_ut1 + dt.delta_t_days();
        return Ok(gmst_iau2006_deg(jd_ut1, jd_tt)?.rem_euclid(360.0));
    }
    let mut g = mean_sidereal_longterm_deg(jd_ut1, dt.delta_t_days())?;
    if jd_ut1 <= SIDEREAL_MODERN_START_JD {
        g -= dt.offset_lo_deg();
    } else {
        g -= dt.offset_hi_deg();
    }
    Ok(g.rem_euclid(360.0))
}
