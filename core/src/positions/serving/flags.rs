//! Request normalization and returned-label policy for position calculations.
//!
//! These helpers resolve provider, center, rate and representation flags without
//! acquiring data or calculating positions. Numerical flags and public labels
//! are separate: an accepted label does not select another ephemeris model.

use crate::constants::{
    FLG_BARYCTR, FLG_CENTER_BODY, FLG_HELCTR, FLG_J2000, FLG_JPLEPH, FLG_JPLHOR, FLG_JPLHOR_APPROX,
    FLG_MOSEPH, FLG_NOABERR, FLG_NOGDEFL, FLG_NONUT, FLG_RADIANS, FLG_SIDEREAL, FLG_SPEED,
    FLG_SPEED3, FLG_SWIEPH, FLG_TRUEPOS, FLG_XYZ, MARS,
};
use crate::errors::{Error, Result};
use crate::positions::resolve_center_flags;

/// Remove request bits that are consumed without being returned.
///
/// CENTER_BODY is inert for Sun through Mars, which have no separately
/// selected satellite-system center here. JPLHOR/JPLHOR_APPROX are accepted
/// but do not advertise an orientation correction that was not applied.
pub fn consume_non_echoed_flags(flags: i32, body_id: i32) -> i32 {
    let mut out = flags;
    if (out & FLG_CENTER_BODY) != 0 && (0..=MARS).contains(&body_id) {
        out &= !FLG_CENTER_BODY;
    }
    if (out & (FLG_JPLHOR | FLG_JPLHOR_APPROX)) != 0 {
        out &= !(FLG_JPLHOR | FLG_JPLHOR_APPROX);
    }
    out
}

/// Return SWIEPH instead of JPLEPH when the requested JPL filename selected
/// the tier-default image. Only the public flag label changes; the numerical
/// state still comes from its acquired JPL DE kernel.
pub fn degrade_jpl_echo_on_fallback(flags: i32, jpl_fallback_active: bool) -> i32 {
    if (flags & FLG_JPLEPH) != 0 && jpl_fallback_active {
        (flags & !FLG_JPLEPH) | FLG_SWIEPH
    } else {
        flags
    }
}

/// Normalize calculation flags for the common TT and UT1 entries.
///
/// Strip MOSEPH from computation, select JPLEPH over SWIEPH (default), imply
/// SPEED for a SPEED3-only request, and resolve center priorities. Explicit
/// SPEED supersedes SPEED3. Returned public labels are handled separately.
pub fn normalize_calc_flags(flags: i32, jpl_fallback_active: bool) -> i32 {
    let mut out = flags & !FLG_MOSEPH;
    out = degrade_jpl_echo_on_fallback(out, jpl_fallback_active);
    if (out & FLG_JPLEPH) != 0 {
        out &= !FLG_SWIEPH;
    } else if (out & FLG_SWIEPH) == 0 {
        out |= FLG_SWIEPH;
    }
    if (out & FLG_SPEED3) != 0 {
        if (out & FLG_SPEED) != 0 {
            out &= !FLG_SPEED3;
        } else {
            out |= FLG_SPEED;
        }
    }
    resolve_center_flags(out)
}

/// Select one public ephemeris bit, preferring JPLEPH, SWIEPH, then MOSEPH.
/// SWIEPH is the default. This selector does not change rate flags.
pub fn exclusive_ephemeris_bit(flags: i32, jpl_fallback_active: bool) -> i32 {
    let degraded = degrade_jpl_echo_on_fallback(flags, jpl_fallback_active);
    let ephmask = FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH;
    let epheflag = if (degraded & FLG_JPLEPH) != 0 {
        FLG_JPLEPH
    } else if (degraded & FLG_SWIEPH) != 0 {
        FLG_SWIEPH
    } else if (degraded & FLG_MOSEPH) != 0 {
        FLG_MOSEPH
    } else {
        FLG_SWIEPH
    };
    (degraded & !ephmask) | epheflag
}

/// Remove XYZ and RADIANS representation bits before numerical calculation.
/// [`super::finalize_output_flags`] applies the selected representation afterward.
pub fn strip_output_flags(flags: i32) -> i32 {
    flags & !FLG_XYZ & !FLG_RADIANS
}

/// Restore the original request's public rate and ephemeris labels.
/// A SPEED3-only request returns SPEED3; a MOSEPH-only request retains that
/// label. Numerical flags have already been normalized independently.
pub fn echo_request_bits(retflag: i32, raw_flags: i32) -> i32 {
    let mut out = retflag;
    if (raw_flags & FLG_SPEED3) != 0 && (raw_flags & FLG_SPEED) == 0 {
        out = (out & !FLG_SPEED) | FLG_SPEED3;
    }
    if (raw_flags & FLG_MOSEPH) != 0 && (raw_flags & (FLG_JPLEPH | FLG_SWIEPH)) == 0 {
        out = (out & !(FLG_SWIEPH | FLG_JPLEPH)) | FLG_MOSEPH;
    }
    out
}

/// Return flags implied by the selected output and apparent-place policy.
/// J2000/sidereal output implies NONUT; centered or geometric output implies
/// NOGDEFL and NOABERR. This function changes only returned labels.
pub fn implied_retflag_bits(flags: i32) -> i32 {
    let mut extra = 0;
    if (flags & (FLG_J2000 | FLG_SIDEREAL)) != 0 {
        extra |= FLG_NONUT;
    }
    if (flags & (FLG_HELCTR | FLG_BARYCTR | FLG_TRUEPOS)) != 0 {
        extra |= FLG_NOGDEFL | FLG_NOABERR;
    }
    extra
}

/// Return only caller-requested ephemeris bits for TT calculations.
/// Unlike UT1 output, a request without any selection bit returns none.
pub fn calc_tt_epheflag_echo(retflag: i32, raw_flags: i32) -> i32 {
    if (raw_flags & (FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH)) == 0 {
        retflag & !(FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH)
    } else {
        retflag
    }
}

/// Reject an exclusively analytical barycentric request after resolving
/// center and ephemeris priorities. The error is
/// [`Error::UnsupportedFlagCombination`], tagged with the calling entry.
pub fn validate_barycentric_moseph(flags: i32, entry: &'static str) -> Result<()> {
    let resolved = resolve_center_flags(flags);
    if (resolved & FLG_BARYCTR) != 0
        && (resolved & FLG_MOSEPH) != 0
        && (resolved & (FLG_JPLEPH | FLG_SWIEPH)) == 0
    {
        return Err(Error::UnsupportedFlagCombination { entry, flags });
    }
    Ok(())
}
