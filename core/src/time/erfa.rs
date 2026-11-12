//! Native ERFA Earth orientation, precession, nutation, and epoch conversion.
//!
//! Full IAU 2000A/2000 B and IAU 2006 series use bundled coefficient tables.
//! Vondrák 2011 pole series supply long-term
//! precession. The CIO locator, Earth rotation and polar-motion matrices compose
//! the GCRS-to-ITRF transform. ERFA-derived routines and tables retain their
//! licenses and attribution; see `docs/sources.md` and THIRD_PARTY_NOTICES.
//! Regenerate tables with `tools/generators/gen_nutation_tables.py`.
//!
//! Epochs are TT or UT1 Julian days as documented by each function. J2000 day
//! offsets preserve sub-JD-ULP time resolution without reconstructing a scalar
//! date. Fundamental arguments use TT centuries. Angles are radians unless a
//! degree-valued API explicitly says otherwise; matrices are row-major.
//! Numerical evaluation uses binary64, fused Horner polynomials and paired
//! [`f64::sin_cos`] evaluation of trigonometric terms.
//!
//! Names here are Rust interfaces, not interchangeable ERFA signatures:
//! [`obl06`] returns **degrees**, while [`era00`] and the nutation angles use
//! radians. [`gd2gc`] returns geocentric Cartesian metres. A row-major rotation
//! acts on column vectors; consult each entry for its source and target frame.
//!
//! Checked entry points refuse non-finite epochs with
//! [`Error::NonFiniteArgument`]. Calendar conversion preserves its documented
//! negative status codes. Internal offset routines propagate IEEE non-finite
//! arithmetic where the caller owns validation. Free core nutation requires
//! external observations and is not included in these deterministic series.
//!
//! IAU 2000A/2006 models have their highest accuracy within a few centuries of
//! J2000; Vondrák precession extends across ±200 kyr with increasing model
//! uncertainty. Tests cover rotation invariants, analytic derivatives, calendar
//! domains and independent ERFA numerical values.

use crate::errors::{Error, Result};
use crate::nutation_tables::{
    NUT_LS_ROWS, NUT_PL_ROWS, NUT00B_ROWS, S06_S0_ROWS, S06_S1_ROWS, S06_S2_ROWS, S06_S3_ROWS,
    S06_S4_ROWS, S06_SP,
};
use crate::positions::{
    self, FRAME_BIAS_J2000, VONDRAK_PQPER_P_COS, VONDRAK_PQPER_P_SIN, VONDRAK_PQPER_PERIODS,
    VONDRAK_PQPER_Q_COS, VONDRAK_PQPER_Q_SIN, VONDRAK_PQPOL, VONDRAK_XYPER_PERIODS,
    VONDRAK_XYPER_X_COS, VONDRAK_XYPER_X_SIN, VONDRAK_XYPER_Y_COS, VONDRAK_XYPER_Y_SIN,
    VONDRAK_XYPOL, fw2m,
};

// Compatible Rust aliases for selected ERFA-related kernels.
// obl 06 intentionally retains this Rust API's degree output; the angle
// routines and rotation matrices document their own units and rounding
// conventions. Native pfw 06 below supplies fused ERFA polynomial evaluation.
pub use crate::positions::{
    FRAME_BIAS_J2000 as BP06_BIAS_MATRIX, earth_rotation_angle_rad as era00, fw2m as fw2m_matrix,
    geodetic_to_itrf_m as gd2gc, mean_obliquity_iau2006_deg as obl06, nutation_matrix as numat,
};

// ---------------------------------------------------------------------------
// ERFA constants (erfam.h, v 2.0.1)
// ---------------------------------------------------------------------------

/// J2000.0 Julian Day (`ERFA_DJ00`).
pub const DJ00: f64 = 2451545.0;
/// Days per Julian century (`ERFA_DJC`).
pub const DJC: f64 = 36525.0;
/// Arcseconds per full turn (`ERFA_TURNAS`), exact.
pub const TURN_ARCSEC: f64 = 1_296_000.0;
/// Arcsecond-to-radian scale (`ERFA_DAS2R`). The literal keeps the full
/// published spelling; both spellings parse to the same binary64 value
/// under correctly-rounded conversion, so the lint is allowed here
/// rather than by truncating the constant.
#[allow(clippy::excessive_precision)]
pub const DAS2R: f64 = 4.8481368110953599359e-6;
/// Milliarcsecond-to-radian scale (`ERFA_DMAS2R = DAS2R / 1e3`).
pub const DMAS2R: f64 = DAS2R / 1e3;
/// 2 pi (`ERFA_D2PI`).
pub const D2PI: f64 = 2.0 * std::f64::consts::PI;
/// MJD zero point (`ERFA_DJM0`).
pub const DJM0: f64 = 2400000.5;
/// MJD offset of J2000.0 (`ERFA_DJM00`).
pub const DJM00: f64 = 51544.5;
/// Tropical-year length in days (`ERFA_DTY`, Besselian epoch scale; consumed
/// by the `epb2jd` Besselian conversion).
pub const DTY: f64 = 365.242198781;

// ---------------------------------------------------------------------------
// Paired trigonometric evaluation
// ---------------------------------------------------------------------------

/// Sine and cosine of an angle in radians using ordinary binary64 arithmetic.
fn era_sincos(x: f64) -> (f64, f64) {
    x.sin_cos()
}

// ---------------------------------------------------------------------------
// Fundamental arguments (canonical home; `lunar` delegates here)
// ---------------------------------------------------------------------------

/// IERS 2003 (Simon et al. 1994) fundamental-argument polynomials in
/// arcseconds, `(c0, c1, c2, c3, c4)` with `t` in Julian centuries past
/// J2000.0. These coefficients also serve the lunar model.
pub(crate) const FA_D03: (f64, f64, f64, f64, f64) = (
    1072260.703692,
    1602961601.2090,
    -6.3706,
    0.006593,
    -0.00003169,
);
pub(crate) const FA_LP03: (f64, f64, f64, f64, f64) = (
    1287104.793048,
    129596581.0481,
    -0.5532,
    0.000136,
    -0.00001149,
);
pub(crate) const FA_L03: (f64, f64, f64, f64, f64) = (
    485868.249036,
    1717915923.2178,
    31.8792,
    0.051635,
    -0.00024470,
);
pub(crate) const FA_F03: (f64, f64, f64, f64, f64) = (
    335779.526232,
    1739527262.8478,
    -12.7512,
    -0.001037,
    0.00000417,
);
pub(crate) const FA_OM03: (f64, f64, f64, f64, f64) =
    (450160.398036, -6962890.5431, 7.4722, 0.007702, -0.00005939);

/// Evaluate one IERS 2003 fundamental-argument polynomial at `t`
/// centuries, in radians.
///
/// Fused-Horner in arcseconds (one rounding per level via fused
/// multiply-add), truncated remainder against one turn, then
/// the arcsecond-to-radian scale.
pub(crate) fn fa_poly_rad(t: f64, poly: (f64, f64, f64, f64, f64)) -> f64 {
    let (c0, c1, c2, c3, c4) = poly;
    let arcsec = c4
        .mul_add(t, c3)
        .mul_add(t, c2)
        .mul_add(t, c1)
        .mul_add(t, c0);
    (arcsec % TURN_ARCSEC) * DAS2R
}

/// TT Julian Day to Julian centuries since J2000.0.
///
/// Subtract the J2000 pivot before dividing by 36525 days. Finite scalar
/// input is required; its original time resolution cannot be recovered by
/// subsequently splitting the date. The IAU model limits stated at module
/// level apply. A non-finite epoch returns [`Error::NonFiniteArgument`].
pub fn tt_centuries(jd_tt: f64) -> Result<f64> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    Ok((jd_tt - DJ00) / DJC)
}

/// Mean anomaly of the Moon, radians (ERFA `fal03`, IERS 2003).
///
/// `t` is finite TT Julian centuries from J2000. The IERS 2003 polynomial
/// supplies an argument of the IAU series; finite extrapolation does not
/// extend their physical validity. Non-finite `t` returns
/// [`Error::NonFiniteArgument`].
pub fn fal03(t: f64) -> Result<f64> {
    if !t.is_finite() {
        return Err(Error::NonFiniteArgument { value: t });
    }
    Ok(fa_poly_rad(t, FA_L03))
}

/// Mean anomaly of the Sun, radians (ERFA `falp03`, IERS 2003).
///
/// `t` is finite TT Julian centuries from J2000. The IERS 2003 polynomial
/// supplies an argument of the IAU series; finite extrapolation does not
/// extend their physical validity. Non-finite `t` returns
/// [`Error::NonFiniteArgument`].
pub fn falp03(t: f64) -> Result<f64> {
    if !t.is_finite() {
        return Err(Error::NonFiniteArgument { value: t });
    }
    Ok(fa_poly_rad(t, FA_LP03))
}

/// Mean argument of the latitude of the Moon, radians (ERFA `faf03`).
///
/// `t` is finite TT Julian centuries from J2000. The IERS 2003 polynomial
/// supplies an argument of the IAU series; finite extrapolation does not
/// extend their physical validity. Non-finite `t` returns
/// [`Error::NonFiniteArgument`].
pub fn faf03(t: f64) -> Result<f64> {
    if !t.is_finite() {
        return Err(Error::NonFiniteArgument { value: t });
    }
    Ok(fa_poly_rad(t, FA_F03))
}

/// Mean elongation of the Moon from the Sun, radians (ERFA `fad03`).
///
/// `t` is finite TT Julian centuries from J2000. The IERS 2003 polynomial
/// supplies an argument of the IAU series; finite extrapolation does not
/// extend their physical validity. Non-finite `t` returns
/// [`Error::NonFiniteArgument`].
pub fn fad03(t: f64) -> Result<f64> {
    if !t.is_finite() {
        return Err(Error::NonFiniteArgument { value: t });
    }
    Ok(fa_poly_rad(t, FA_D03))
}

/// Mean longitude of the Moon's ascending node, radians (ERFA `faom03`).
///
/// `t` is finite TT Julian centuries from J2000. The IERS 2003 polynomial
/// supplies an argument of the IAU series; finite extrapolation does not
/// extend their physical validity. Non-finite `t` returns
/// [`Error::NonFiniteArgument`].
pub fn faom03(t: f64) -> Result<f64> {
    if !t.is_finite() {
        return Err(Error::NonFiniteArgument { value: t });
    }
    Ok(fa_poly_rad(t, FA_OM03))
}

/// Mean longitude of Mercury, radians (ERFA `fame03`, IERS 2003).
///
/// `t` is finite TT Julian centuries from J2000. The IERS 2003 polynomial
/// supplies an argument of the IAU series; finite extrapolation does not
/// extend their physical validity. Non-finite `t` returns
/// [`Error::NonFiniteArgument`].
pub fn fame03(t: f64) -> Result<f64> {
    if !t.is_finite() {
        return Err(Error::NonFiniteArgument { value: t });
    }
    Ok(2608.7903141574f64.mul_add(t, 4.402608842) % D2PI)
}

/// Mean longitude of Venus, radians (ERFA `fave03`, IERS 2003).
///
/// `t` is finite TT Julian centuries from J2000. The IERS 2003 polynomial
/// supplies an argument of the IAU series; finite extrapolation does not
/// extend their physical validity. Non-finite `t` returns
/// [`Error::NonFiniteArgument`].
pub fn fave03(t: f64) -> Result<f64> {
    if !t.is_finite() {
        return Err(Error::NonFiniteArgument { value: t });
    }
    Ok(1021.3285546211f64.mul_add(t, 3.176146697) % D2PI)
}

/// Mean longitude of Earth, radians (ERFA `fae03`, IERS 2003).
///
/// `t` is finite TT Julian centuries from J2000. The IERS 2003 polynomial
/// supplies an argument of the IAU series; finite extrapolation does not
/// extend their physical validity. Non-finite `t` returns
/// [`Error::NonFiniteArgument`].
pub fn fae03(t: f64) -> Result<f64> {
    if !t.is_finite() {
        return Err(Error::NonFiniteArgument { value: t });
    }
    Ok(628.3075849991f64.mul_add(t, 1.753470314) % D2PI)
}

/// Mean longitude of Mars, radians (ERFA `fama03`, IERS 2003).
///
/// `t` is finite TT Julian centuries from J2000. The IERS 2003 polynomial
/// supplies an argument of the IAU series; finite extrapolation does not
/// extend their physical validity. Non-finite `t` returns
/// [`Error::NonFiniteArgument`].
pub fn fama03(t: f64) -> Result<f64> {
    if !t.is_finite() {
        return Err(Error::NonFiniteArgument { value: t });
    }
    Ok(334.0612426700f64.mul_add(t, 6.203480913) % D2PI)
}

/// Mean longitude of Jupiter, radians (ERFA `faju03`, IERS 2003).
///
/// `t` is finite TT Julian centuries from J2000. The IERS 2003 polynomial
/// supplies an argument of the IAU series; finite extrapolation does not
/// extend their physical validity. Non-finite `t` returns
/// [`Error::NonFiniteArgument`].
pub fn faju03(t: f64) -> Result<f64> {
    if !t.is_finite() {
        return Err(Error::NonFiniteArgument { value: t });
    }
    Ok(52.9690962641f64.mul_add(t, 0.599546497) % D2PI)
}

/// Mean longitude of Saturn, radians (ERFA `fasa03`, IERS 2003).
///
/// `t` is finite TT Julian centuries from J2000. The IERS 2003 polynomial
/// supplies an argument of the IAU series; finite extrapolation does not
/// extend their physical validity. Non-finite `t` returns
/// [`Error::NonFiniteArgument`].
pub fn fasa03(t: f64) -> Result<f64> {
    if !t.is_finite() {
        return Err(Error::NonFiniteArgument { value: t });
    }
    Ok(21.3299104960f64.mul_add(t, 0.874016757) % D2PI)
}

/// Mean longitude of Uranus, radians (ERFA `faur03`, IERS 2003).
///
/// `t` is finite TT Julian centuries from J2000. The IERS 2003 polynomial
/// supplies an argument of the IAU series; finite extrapolation does not
/// extend their physical validity. Non-finite `t` returns
/// [`Error::NonFiniteArgument`].
pub fn faur03(t: f64) -> Result<f64> {
    if !t.is_finite() {
        return Err(Error::NonFiniteArgument { value: t });
    }
    Ok(7.4781598567f64.mul_add(t, 5.481293872) % D2PI)
}

/// General accumulated precession in longitude, radians (ERFA `fapa03`).
///
/// As written: `(0.024381750 + 0.00000538691 * t) * t`, no range reduction.
///
/// `t` is finite TT Julian centuries from J2000. The IERS 2003 polynomial
/// supplies an argument of the IAU series; finite extrapolation does not
/// extend their physical validity. Non-finite `t` returns
/// [`Error::NonFiniteArgument`].
pub fn fapa03(t: f64) -> Result<f64> {
    if !t.is_finite() {
        return Err(Error::NonFiniteArgument { value: t });
    }
    Ok(0.00000538691f64.mul_add(t, 0.024381750) * t)
}

// ---------------------------------------------------------------------------
// MHB 2000 inline arguments (nut 00 a private Delaunay/planetary forms)
// ---------------------------------------------------------------------------

/// Mean anomaly of the Sun, MHB 2000 inline polynomial (radians).
fn mhb_sun_anomaly(t: f64) -> f64 {
    (-0.00001149f64)
        .mul_add(t, 0.000136)
        .mul_add(t, -0.5532)
        .mul_add(t, 129596581.0481)
        .mul_add(t, 1287104.79305)
        % TURN_ARCSEC
        * DAS2R
}

/// Mean elongation of the Moon from the Sun, MHB 2000 inline (radians).
fn mhb_moon_elongation(t: f64) -> f64 {
    (-0.00003169f64)
        .mul_add(t, 0.006593)
        .mul_add(t, -6.3706)
        .mul_add(t, 1602961601.2090)
        .mul_add(t, 1072260.70369)
        % TURN_ARCSEC
        * DAS2R
}

// ---------------------------------------------------------------------------
// Nutation series
// ---------------------------------------------------------------------------

/// 0.1-microarcsecond to radians (`ERFA_DAS2R / 1e7`).
const U2R: f64 = DAS2R / 1e7;

/// Dot product of integer multipliers and fundamental arguments, using fused
/// multiply-add to reduce intermediate rounding.
fn fused_arg_sum(mult: &[i8], fa: &[f64]) -> f64 {
    debug_assert_eq!(mult.len(), fa.len());
    let mut acc = (mult[0] as f64).mul_add(fa[0], mult[1] as f64 * fa[1]);
    for k in 2..mult.len() {
        acc = (mult[k] as f64).mul_add(fa[k], acc);
    }
    acc
}

/// Evaluate independent argument sums without changing any scalar FMA chain.
///
/// The native callers use two or four lanes and five or thirteen columns.
/// Every multiplier, including zero, remains in its original per-row position;
/// compiler scheduling is free to overlap the independent lanes.
fn interleaved_fused_arg_sums<const N: usize, const LANES: usize>(
    multipliers: [&[i8; N]; LANES],
    fa: &[f64; N],
) -> [f64; LANES] {
    let mut sums: [f64; LANES] = std::array::from_fn(|lane| {
        let mult = multipliers[lane];
        (mult[0] as f64).mul_add(fa[0], mult[1] as f64 * fa[1])
    });
    for (k, &argument) in fa.iter().enumerate().skip(2) {
        for (sum, mult) in sums.iter_mut().zip(multipliers.iter()) {
            *sum = (mult[k] as f64).mul_add(argument, *sum);
        }
    }
    sums
}

/// Evaluate two independent argument sums with the shared scalar FMA kernel.
fn paired_fused_arg_sum<const N: usize>(
    first: &[i8; N],
    second: &[i8; N],
    fa: &[f64; N],
) -> [f64; 2] {
    interleaved_fused_arg_sums([first, second], fa)
}

/// Exact binary64 remainder by the stored `D2PI` word for small quotients.
///
/// Below 2^52 turns the truncated rounded quotient is the true integer quotient
/// or one greater. FMA forms the exact residual: both operands are multiples of
/// the divisor's ulp and a residual smaller than the divisor fits in binary64.
/// A negative residual identifies the one-too-large quotient and is corrected
/// by adding the divisor. The final sign preserves negative zero. Larger or
/// non-finite arguments use the ordinary remainder operation.
#[inline]
fn nutation_argument_remainder(angle: f64) -> f64 {
    let magnitude = angle.abs();
    if magnitude < D2PI {
        angle
    } else if magnitude < D2PI * (1_u64 << 52) as f64 {
        let quotient = (magnitude / D2PI).trunc();
        let residual = (-quotient).mul_add(D2PI, magnitude);
        let residual = if residual < 0.0 {
            residual + D2PI
        } else {
            residual
        };
        residual.copysign(angle)
    } else {
        angle % D2PI
    }
}

/// Luni-solar nutation accumulation over one coefficient table.
///
/// Fused argument sums and paired trigonometric evaluation accumulate
/// `dp += (sp + spt*t)*sin + cp*cos`,
/// `de += (ce + cet*t)*cos + se*sin`. Shared by the full MHB 2000
/// series (`NUT_LS_ROWS`) and the 77-term IAU 2000 B series
/// (`NUT00B_ROWS`), which differ only in the table.
fn lunisolar_sum(
    rows: &[([i8; 5], [f64; 6])],
    t: f64,
    el: f64,
    elp: f64,
    f: f64,
    d: f64,
    om: f64,
) -> (f64, f64) {
    let fa = [el, elp, f, d, om];
    let mut dp = 0.0;
    let mut de = 0.0;
    for group in rows.rchunks(4) {
        // Only pure argument sums overlap. Trigonometry and accumulation keep
        // descending row order; short tails evaluate only their actual rows.
        let arguments = match group {
            [a, b, c, d] => interleaved_fused_arg_sums([&d.0, &c.0, &b.0, &a.0], &fa),
            [a, b, c] => {
                let [latest, middle] = paired_fused_arg_sum(&c.0, &b.0, &fa);
                [latest, middle, fused_arg_sum(&a.0, &fa), 0.0]
            }
            [a, b] => {
                let [latest, earliest] = paired_fused_arg_sum(&b.0, &a.0, &fa);
                [latest, earliest, 0.0, 0.0]
            }
            [a] => [fused_arg_sum(&a.0, &fa), 0.0, 0.0, 0.0],
            _ => unreachable!("rchunks(4) yields one to four rows"),
        };
        for ((_, coeff), argument) in group.iter().rev().zip(arguments) {
            let (sarg, carg) = era_sincos(nutation_argument_remainder(argument));
            dp += coeff[1].mul_add(t, coeff[0]).mul_add(sarg, coeff[2] * carg);
            de += coeff[4].mul_add(t, coeff[3]).mul_add(carg, coeff[5] * sarg);
        }
    }
    (dp, de)
}

/// Planetary MHB 2000 accumulation over `NUT_PL_ROWS` (no `t` terms).
// Thirteen scalar arguments mirror the ERFA `xpl` struct field-for-field
// (one per multiplier column); bundling them would obscure the mapping.
#[allow(clippy::too_many_arguments)]
fn planetary_sum(
    al: f64,
    af: f64,
    ad: f64,
    aom: f64,
    alme: f64,
    alve: f64,
    alea: f64,
    alma: f64,
    alju: f64,
    alsa: f64,
    alur: f64,
    alne: f64,
    apa: f64,
) -> (f64, f64) {
    let fa = [
        al, af, ad, aom, alme, alve, alea, alma, alju, alsa, alur, alne, apa,
    ];
    let mut dp = 0.0;
    let mut de = 0.0;
    for group in NUT_PL_ROWS.rchunks(4) {
        let arguments = match group {
            [a, b, c, d] => interleaved_fused_arg_sums([&d.0, &c.0, &b.0, &a.0], &fa),
            [a, b, c] => {
                let [latest, middle] = paired_fused_arg_sum(&c.0, &b.0, &fa);
                [latest, middle, fused_arg_sum(&a.0, &fa), 0.0]
            }
            [a, b] => {
                let [latest, earliest] = paired_fused_arg_sum(&b.0, &a.0, &fa);
                [latest, earliest, 0.0, 0.0]
            }
            [a] => [fused_arg_sum(&a.0, &fa), 0.0, 0.0, 0.0],
            _ => unreachable!("rchunks(4) yields one to four rows"),
        };
        for ((_, coeff), argument) in group.iter().rev().zip(arguments) {
            let (sarg, carg) = era_sincos(nutation_argument_remainder(argument));
            dp += coeff[0].mul_add(sarg, coeff[1] * carg);
            de += coeff[2].mul_add(sarg, coeff[3] * carg);
        }
    }
    (dp, de)
}

/// IAU 2000A nutation (MHB 2000, free core nutation omitted), radians.
/// TT Julian Day in; 1365-term luni-solar plus planetary series;
/// non-finite input refused via `tt_centuries`.
///
/// Native ERFA `nut00a` equivalent: luni-solar series with the IERS
/// 2003Delaunay set plus the MHB 2000 inline Sun/elontation forms,
/// then the planetary series with the MHB 2000 private lunar arguments
/// and IERS 2003 planetary longitudes (the two-argument-set behavior
/// is reproduced faithfully; unifying them changes results below
/// 0.1 microarcsecond).
///
/// TT is a finite Julian day. The IAU series is intended for epochs within
/// a few centuries of J2000 and excludes observed free-core nutation.
/// Non-finite TT returns [`Error::NonFiniteArgument`].
pub fn nut00a(jd_tt: f64) -> Result<(f64, f64)> {
    let t = tt_centuries(jd_tt)?;
    Ok(nut00a_cached_centuries(t))
}

/// Reuse raw IAU 2000A angles only for identical evaluated century words.
/// Checked dates and retained offsets keep their own epoch arithmetic; no
/// scalar JD is reconstructed. The existing native table keeps its capacity
/// and remains separate from the other nutation model. Caller-owned nonfinite
/// input follows the original arithmetic without cache acquisition.
fn nut00a_cached_centuries(t: f64) -> (f64, f64) {
    if !t.is_finite() {
        return nut00a_centuries(t);
    }
    let key = t.to_bits();
    if let Some(angles) = crate::pure_cache::erfa_nutation(key) {
        return angles;
    }
    let angles = nut00a_centuries(t);
    crate::pure_cache::retain_erfa_nutation(key, angles);
    angles
}

/// IAU 2000A series at internally acquired TT centuries.
/// Uses the checked entry's series for an
/// internally acquired epoch whose validation is owned by the caller.
/// `t = (JD_TT - J2000)/36525` is Julian centuries of TT. Returns the
/// `(dpsi, deps)` longitude/obliquity nutation tuple in radians, with the
/// series and accuracy limitations documented on [`nut00a`].
fn nut00a_centuries(t: f64) -> (f64, f64) {
    let el = fa_poly_rad(t, FA_L03);
    let elp = mhb_sun_anomaly(t);
    let f = fa_poly_rad(t, FA_F03);
    let d = mhb_moon_elongation(t);
    let om = fa_poly_rad(t, FA_OM03);
    let (dp_ls, de_ls) = lunisolar_sum(&NUT_LS_ROWS, t, el, elp, f, d, om);

    let al = 8328.6914269554f64.mul_add(t, 2.35555598) % D2PI;
    let af = 8433.466158131f64.mul_add(t, 1.627905234) % D2PI;
    let ad = 7771.3771468121f64.mul_add(t, 5.198466741) % D2PI;
    let aom = (-33.757045f64).mul_add(t, 2.18243920) % D2PI;
    let apa = 0.00000538691f64.mul_add(t, 0.024381750) * t;
    let alme = 2608.7903141574f64.mul_add(t, 4.402608842) % D2PI;
    let alve = 1021.3285546211f64.mul_add(t, 3.176146697) % D2PI;
    let alea = 628.3075849991f64.mul_add(t, 1.753470314) % D2PI;
    let alma = 334.0612426700f64.mul_add(t, 6.203480913) % D2PI;
    let alju = 52.9690962641f64.mul_add(t, 0.599546497) % D2PI;
    let alsa = 21.3299104960f64.mul_add(t, 0.874016757) % D2PI;
    let alur = 7.4781598567f64.mul_add(t, 5.481293872) % D2PI;
    let alne = 3.8127774000f64.mul_add(t, 5.321159000) % D2PI;
    let (dp_pl, de_pl) = planetary_sum(
        al, af, ad, aom, alme, alve, alea, alma, alju, alsa, alur, alne, apa,
    );
    // Convert both accumulated nutation series to radians.
    (dp_ls * U2R + dp_pl * U2R, de_ls * U2R + de_pl * U2R)
}

/// IAU 2000 B abbreviated nutation, radians (native ERFA `nut00b`).
///
/// 77-term luni-solar series with linear-only Delaunay arguments plus
/// the fixed planetary-bias offsets (`-0.135`/`+0.388` mas) that
/// approximate the omitted planetary nutation.
///
/// TT is a finite Julian day. The IAU series is intended for epochs within
/// a few centuries of J2000 and excludes observed free-core nutation.
/// Non-finite TT returns [`Error::NonFiniteArgument`].
pub fn nut00b(jd_tt: f64) -> Result<(f64, f64)> {
    Ok(nut00b_centuries(tt_centuries(jd_tt)?))
}

/// Evaluate the 77-term nutation series at TT Julian centuries from J2000.
/// Both returned angles are radians. The fixed planetary offsets and the
/// linear fundamental arguments are those of [`nut00b`]; nonfinite arguments
/// propagate through the series without another public epoch admission.
fn nut00b_centuries(t: f64) -> (f64, f64) {
    let el = 1717915923.2178f64.mul_add(t, 485868.249036) % TURN_ARCSEC * DAS2R;
    let elp = 129596581.0481f64.mul_add(t, 1287104.79305) % TURN_ARCSEC * DAS2R;
    let f = 1739527262.8478f64.mul_add(t, 335779.526232) % TURN_ARCSEC * DAS2R;
    let d = 1602961601.2090f64.mul_add(t, 1072260.70369) % TURN_ARCSEC * DAS2R;
    let om = (-6962890.5431f64).mul_add(t, 450160.398036) % TURN_ARCSEC * DAS2R;
    let (dp, de) = lunisolar_sum(&NUT00B_ROWS, t, el, elp, f, d, om);
    (dp * U2R - 0.135 * DMAS2R, de * U2R + 0.388 * DMAS2R)
}

/// IAU 2000 B nutation from an already acquired TT day offset from J2000.
///
/// The input is `TT - 2451545.0`, in days. Dividing this offset directly by
/// 36525 preserves sub-Julian-day precision without adding and subtracting
/// J2000 again. The returned longitude and obliquity nutations are radians,
/// with the same 77-term model and accuracy limitations as [`nut00b`].
///
/// Callers provide a finite offset; invalid arithmetic produces non-finite
/// values. Use [`nut00b`] for checked scalar-Julian-day input.
pub fn source_nut00b_from_j2000_offset(offset_days: f64) -> (f64, f64) {
    nut00b_centuries(offset_days / DJC)
}

/// IAU 2006/2000 A nutation, radians (native ERFA `nut06a`).
///
/// The IAU 2000A pair with the P 03 adjustments for the changed
/// precession rates and the J 2 secular variation (Wallace &
/// Capitaine 2006, Eqs. 5):
/// `dpsi += dp * (0.4697e-6 + fj2)`, `deps += de * fj2`.
///
/// TT is a finite Julian day. The IAU series is intended for epochs within
/// a few centuries of J2000 and excludes observed free-core nutation.
/// Non-finite TT returns [`Error::NonFiniteArgument`].
pub fn nut06a(jd_tt: f64) -> Result<(f64, f64)> {
    let t = tt_centuries(jd_tt)?;
    let (dp, de) = nut00a(jd_tt)?;
    let fj2 = -2.7774e-6 * t;
    Ok((dp + dp * (0.4697e-6 + fj2), de + de * fj2))
}

/// Analytic TT derivative of IAU 2006/2000 A obliquity nutation, radians/day.
///
/// Differentiates every luni-solar and planetary term used by [`nut06a`],
/// including the fundamental-argument polynomials, secular amplitudes and P 03
/// J 2 correction. No epoch step or fitted rate is used. TT is a finite Julian
/// day; the value model's validity and missing free-core-nutation limitations
/// apply. Non-finite epochs return [`Error::NonFiniteArgument`].
pub(crate) fn nut06a_obliquity_rate(jd_tt: f64) -> Result<f64> {
    nut06a_rates(jd_tt).map(|(_, obliquity)| obliquity)
}

/// Analytic TT derivative of IAU 2006/2000A longitude nutation, radians/day.
///
/// The derivative includes every luni-solar and planetary term, secular
/// amplitude and P03 product rule in the same value model as [`nut06a`].
/// No finite-difference step, observed EOP or free-core nutation is used.
/// Nonfinite TT is rejected; the value series' physical validity applies.
pub(crate) fn nut06a_longitude_rate(jd_tt: f64) -> Result<f64> {
    nut06a_rates(jd_tt).map(|(longitude, _)| longitude)
}

/// Shared phase evaluation and analytic TT derivatives, in radians/day.
/// Keep the obliquity accumulation order independent of the longitude sum.
fn nut06a_rates(jd_tt: f64) -> Result<(f64, f64)> {
    Ok(nut06a_rates_centuries(tt_centuries(jd_tt)?))
}

/// Longitude-nutation rate at an already acquired scalar TT, radians/day.
/// Invalid internal clocks propagate IEEE arithmetic, like `source_nut06a`;
/// public clock validation and its error precedence belong to the caller.
pub(crate) fn source_nut06a_longitude_rate(jd_tt: f64) -> f64 {
    nut06a_rates_centuries((jd_tt - DJ00) / DJC).0
}

fn nut06a_rates_centuries(t: f64) -> (f64, f64) {
    let polynomial_rate = |(_, c1, c2, c3, c4): (f64, f64, f64, f64, f64)| {
        (4.0 * c4)
            .mul_add(t, 3.0 * c3)
            .mul_add(t, 2.0 * c2)
            .mul_add(t, c1)
            * DAS2R
    };
    let arguments = [
        fa_poly_rad(t, FA_L03),
        mhb_sun_anomaly(t),
        fa_poly_rad(t, FA_F03),
        mhb_moon_elongation(t),
        fa_poly_rad(t, FA_OM03),
    ];
    let rates = [FA_L03, FA_LP03, FA_F03, FA_D03, FA_OM03].map(polynomial_rate);
    let (mut value, mut rate) = (0.0, 0.0);
    let (mut longitude_value, mut longitude_rate) = (0.0, 0.0);
    for (multipliers, coefficients) in NUT_LS_ROWS.iter().rev() {
        let (sine, cosine) = era_sincos(fused_arg_sum(multipliers, &arguments) % D2PI);
        let phase_rate = fused_arg_sum(multipliers, &rates);
        let amplitude = coefficients[4].mul_add(t, coefficients[3]);
        value += amplitude.mul_add(cosine, coefficients[5] * sine);
        rate +=
            coefficients[4] * cosine + (-amplitude * sine + coefficients[5] * cosine) * phase_rate;
        let longitude_amplitude = coefficients[1].mul_add(t, coefficients[0]);
        longitude_value += longitude_amplitude.mul_add(sine, coefficients[2] * cosine);
        longitude_rate += coefficients[1] * sine
            + (longitude_amplitude * cosine - coefficients[2] * sine) * phase_rate;
    }
    let frequencies = [
        8328.6914269554_f64,
        8433.466158131,
        7771.3771468121,
        -33.757045,
        2608.7903141574,
        1021.3285546211,
        628.3075849991,
        334.0612426700,
        52.9690962641,
        21.3299104960,
        7.4781598567,
        3.8127774000,
    ];
    let phases = [
        2.35555598,
        1.627905234,
        5.198466741,
        2.18243920,
        4.402608842,
        3.176146697,
        1.753470314,
        6.203480913,
        0.599546497,
        0.874016757,
        5.481293872,
        5.321159000,
    ];
    let mut planetary_arguments = [0.0; 13];
    let mut planetary_rates = [0.0; 13];
    for i in 0..12 {
        planetary_arguments[i] = frequencies[i].mul_add(t, phases[i]) % D2PI;
        planetary_rates[i] = frequencies[i];
    }
    planetary_arguments[12] = 0.00000538691f64.mul_add(t, 0.024381750) * t;
    planetary_rates[12] = 0.024381750 + 2.0 * 0.00000538691 * t;
    for (multipliers, coefficients) in NUT_PL_ROWS.iter().rev() {
        let (sine, cosine) = era_sincos(fused_arg_sum(multipliers, &planetary_arguments) % D2PI);
        let phase_rate = fused_arg_sum(multipliers, &planetary_rates);
        value += coefficients[2].mul_add(sine, coefficients[3] * cosine);
        rate += (coefficients[2] * cosine - coefficients[3] * sine) * phase_rate;
        longitude_value += coefficients[0].mul_add(sine, coefficients[1] * cosine);
        longitude_rate += (coefficients[0] * cosine - coefficients[1] * sine) * phase_rate;
    }
    // Product rules for the same P03 longitude/obliquity value adjustments.
    let longitude_rate =
        ((1.0 + 0.4697e-6 - 2.7774e-6 * t) * longitude_rate - 2.7774e-6 * longitude_value) * U2R
            / DJC;
    let obliquity_rate = ((1.0 - 2.7774e-6 * t) * rate - 2.7774e-6 * value) * U2R / DJC;
    (longitude_rate, obliquity_rate)
}

/// IAU 2006/2000 A nutation at an already acquired scalar TT.
///
/// Radians are returned with the same series and P 03 adjustments as `nut06a`.
/// Callers provide finite TT; invalid arithmetic produces non-finite values.
/// Use [`nut06a`] for a checked scalar-Julian-day argument. No observational
/// EOP or free-core nutation enters this deterministic series.
pub fn source_nut06a(jd_tt: f64) -> (f64, f64) {
    source_nut06a_from_j2000_offset(jd_tt - DJ00)
}

/// Nutation from an already acquired TT day offset from J2000.
/// The offset is divided by 36525 without reconstructing a scalar Julian day.
/// Outputs are radians, with the same series, limitations and caller-owned
/// finite-input requirement as [`source_nut06a`].
pub fn source_nut06a_from_j2000_offset(offset_days: f64) -> (f64, f64) {
    let t = offset_days / DJC;
    let (dp, de) = nut00a_cached_centuries(t);
    let fj2 = -2.7774e-6 * t;
    (dp + dp * (0.4697e-6 + fj2), de + de * fj2)
}

/// Identity, term count and published accuracy of the native nutation model.
/// Model precision describes the IAU 2006/2000 A series, not a measured bound
/// on a particular floating-point evaluation or observational epoch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NutationModelInfo {
    /// IAU model tag (`"IAU2006_2000A"`).
    pub model: &'static str,
    /// Series rows summed by [`nut06a`].
    pub terms: u32,
    /// Published model precision (`"~0.01-0.05 mas"`).
    pub precision: &'static str,
    /// Evaluation provenance (`"native"`).
    pub source: &'static str,
}

/// Report the nutation model evaluated by [`nut06a`].
///
/// Infallible: the identity is a compile-time property of the shipped
/// tables, independent of epoch, context and provider state.
pub fn nutation_model_info() -> NutationModelInfo {
    NutationModelInfo {
        model: "IAU2006_2000A",
        terms: (NUT_LS_ROWS.len() + NUT_PL_ROWS.len()) as u32,
        precision: "~0.01-0.05 mas",
        source: "native",
    }
}

/// IAU 2006 Fukushima-Williams bias/precession angles at TT Julian Day.
///
/// Return `(gamb, phib, psib, epsa)` in radians: ecliptic-node orientation,
/// inclination, precession in longitude and mean obliquity, respectively.
/// Evaluate the IERS 2010 fifth-degree arcsecond polynomials at Julian
/// centuries from J2000 with ordered fused Horner accumulation, then scale
/// to radians (IERS Conventions 2010, Table 5.1; Fukushima 2003 and
/// Hilton et al. 2006). This realization differs from the separate-product
/// [`crate::positions::iau2006_precession_angles`] under rounding.
/// Non-finite input and overflowing evaluations are refused; model accuracy
/// is limited by IAU 2006 validity, not certified by finite extrapolation.
pub fn pfw06(jd_tt: f64) -> Result<(f64, f64, f64, f64)> {
    let t = tt_centuries(jd_tt)?;
    let (gamb, phib, psib, epsa) = pfw06_centuries(t);
    for value in [gamb, phib, psib, epsa] {
        if !value.is_finite() {
            return Err(Error::MathDomain { value });
        }
    }
    Ok((gamb, phib, psib, epsa))
}

/// The four radian Fukushima-Williams polynomials at acquired TT centuries.
/// Input/output validation belongs to the checked or array-consuming caller.
fn pfw06_centuries(t: f64) -> (f64, f64, f64, f64) {
    let gamb = fw_angle_polynomial(
        t,
        [
            -0.052928,
            10.556378,
            0.4932044,
            -0.00031238,
            -0.000002788,
            0.0000000260,
        ],
    );
    let phib = fw_angle_polynomial(
        t,
        [
            84381.412819,
            -46.811016,
            0.0511268,
            0.00053289,
            -0.000000440,
            -0.0000000176,
        ],
    );
    let psib = fw_angle_polynomial(
        t,
        [
            -0.041775,
            5038.481484,
            1.5584175,
            -0.00018522,
            -0.000026452,
            -0.0000000148,
        ],
    );
    let epsa = fw_angle_polynomial(
        t,
        [
            84381.406,
            -46.836769,
            -0.0001831,
            0.00200340,
            -0.000000576,
            -0.0000000434,
        ],
    );
    (gamb, phib, psib, epsa)
}

/// Fused Horner value for one ascending-order arcsecond polynomial.
/// The dimensionless century argument is validated by the caller; return
/// radians and leave the final arithmetic-finiteness check to that boundary.
fn fw_angle_polynomial(t: f64, coefficients: [f64; 6]) -> f64 {
    coefficients
        .into_iter()
        .rev()
        .fold(0.0_f64, |value, coefficient| value.mul_add(t, coefficient))
        * DAS2R
}

/// IAU 2006/2000 A bias-precession-nutation matrix at TT Julian Day.
///
/// Return a dimensionless row-major map from GCRS column vectors to the
/// true equator/equinox of date, including frame bias. Acquire the radian
/// Fukushima-Williams angles with [`pfw06`] and nutation `(dpsi, deps)`
/// with [`nut06a`], then evaluate
/// `fw2m(gamb, phib, psib + dpsi, epsa + deps)` in ordered pre-rotations.
/// Input and arithmetic errors propagate from those routines; the model
/// domain/accuracy is their IAU 2006/2000 A domain, not arbitrary-epoch
/// astronomical certification. No Earth rotation or polar motion is included.
/// Factoring these same rotations is equivalent in exact arithmetic; this
/// assembly does not imply a universal milliarcsecond accuracy advantage over
/// separate matrix products. Ordered fused evaluation controls rounding only.
pub fn pnm06a(jd_tt: f64) -> Result<[[f64; 3]; 3]> {
    let (gamb, phib, psib, epsa) = pfw06(jd_tt)?;
    let (dpsi, deps) = nut06a(jd_tt)?;
    fw2m(gamb, phib, psib + dpsi, epsa + deps)
}

// ---------------------------------------------------------------------------
// Long-term poles, CIO locator, full frame assembly
// ---------------------------------------------------------------------------

/// Julian epoch (e.g. 2000.0) to TT Julian Day (inverse of
/// `epj = 2000 + (jd - DJ00) / 365.25`).
///
/// This scale definition has no astronomical model interval. Finite input is
/// required; precision is limited by the input's floating-point resolution.
/// Non-finite input returns [`Error::NonFiniteArgument`].
pub fn julian_epoch_to_jd(epj: f64) -> Result<f64> {
    if !epj.is_finite() {
        return Err(Error::NonFiniteArgument { value: epj });
    }
    Ok(DJ00 + (epj - 2000.0) * 365.25)
}

/// TT Julian Day to Julian epoch (inverse of [`julian_epoch_to_jd`]:
/// `epj = 2000 + (jd - DJ00) / 365.25`).
///
/// This scale definition has no astronomical model interval. Finite input is
/// required; precision is limited by the input's floating-point resolution.
/// Non-finite input returns [`Error::NonFiniteArgument`].
pub fn julian_epoch(jd_tt: f64) -> Result<f64> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    Ok(2000.0 + (jd_tt - DJ00) / 365.25)
}

// Long-term ERFA algorithms: Copyright (C) 2013-2023, NumFOCUS Foundation.
// Derived from ERFA v2.0.1; BSD terms/disclaimer in THIRD_PARTY_NOTICES.md.

/// Vondrák 2011/2012 mean ecliptic pole at Julian epoch `epj` (TT).
///
/// Returns a dimensionless unit vector in the J2000 mean equatorial
/// frame. Evaluate eight periodic terms and a cubic polynomial for
/// `P_A, Q_A` at `t = (epj - 2000) / 100`; these are dimensionless pole
/// coordinates after scaling the arcsecond-form coefficients. Define
/// `Z = sqrt(max(0, 1 - P_A² - Q_A²))` and the J2000 obliquity
/// `ε₀ = 84381.406 arcsec`. Return
/// `(P_A, -Q_A cos ε₀ - Z sin ε₀, -Q_A sin ε₀ + Z cos ε₀)`.
/// A scalar-Julian-Day intermediate is deliberately avoided: it loses
/// fine epoch resolution. Overflow is refused
/// before clamping the square-root radicand.
///
/// The published model covers ±200,000 years about J2000, with physical
/// uncertainty growing from modern microarcseconds to tenths of a degree
/// at the extremes (Vondrák, Capitaine & Wallace 2011, A&A 534 A 22;
/// 2012 corrigendum, A&A 541 C 1). Numerical identity is not physical
/// accuracy certification. Non-finite epochs or overflowing evaluations
/// are refused; finite extrapolation does not extend model validity.
pub fn ltpecl(epj: f64) -> Result<[f64; 3]> {
    let t = long_term_centuries(epj)?;
    let w = D2PI * t;
    let (mut p, mut q) = (0.0, 0.0);
    for i in 0..VONDRAK_PQPER_PERIODS.len() {
        let (s, c) = era_sincos(w / VONDRAK_PQPER_PERIODS[i]);
        p += c.mul_add(VONDRAK_PQPER_P_COS[i], s * VONDRAK_PQPER_P_SIN[i]);
        q += c.mul_add(VONDRAK_PQPER_Q_COS[i], s * VONDRAK_PQPER_Q_SIN[i]);
    }
    let mut w = 1.0;
    for (pc, qc) in VONDRAK_PQPOL {
        p = pc.mul_add(w, p);
        q = qc.mul_add(w, q);
        w *= t;
    }
    p *= DAS2R;
    q *= DAS2R;
    let z = pole_z(p, q)?;
    let (s, c) = era_sincos(84381.406 * DAS2R);
    checked_pole([p, (-q).mul_add(c, -z * s), (-q).mul_add(s, z * c)])
}

/// Vondrák 2011/2012 mean equator pole at Julian epoch `epj` (TT).
///
/// Returns the dimensionless J2000-equatorial pole `(X, Y, Z)`, using
/// fourteen periodic terms plus cubic polynomials at the direct
/// Julian-century argument `(epj - 2000) / 100`. `Z` is the non-negative
/// square root of `max(0, 1 - X² - Y²)`. The method, ±200,000-year model
/// validity, growing physical uncertainty and refusals are as in [`ltpecl`];
/// no scalar-JD roundtrip or independently sourced sidereal pole is used.
pub fn ltpequ(epj: f64) -> Result<[f64; 3]> {
    let t = long_term_centuries(epj)?;
    let w = D2PI * t;
    let (mut x, mut y) = (0.0, 0.0);
    for i in 0..VONDRAK_XYPER_PERIODS.len() {
        let (s, c) = era_sincos(w / VONDRAK_XYPER_PERIODS[i]);
        x += c.mul_add(VONDRAK_XYPER_X_COS[i], s * VONDRAK_XYPER_X_SIN[i]);
        y += c.mul_add(VONDRAK_XYPER_Y_COS[i], s * VONDRAK_XYPER_Y_SIN[i]);
    }
    let mut w = 1.0;
    for (xc, yc) in VONDRAK_XYPOL {
        x = xc.mul_add(w, x);
        y = yc.mul_add(w, y);
        w *= t;
    }
    x *= DAS2R;
    y *= DAS2R;
    checked_pole([x, y, pole_z(x, y)?])
}

/// Non-negative pole ordinate from two dimensionless components.
/// Refuse non-finite arithmetic before clamping a negative radicand caused by
/// extrapolation. This does not certify the model outside its published span.
fn pole_z(a: f64, b: f64) -> Result<f64> {
    let radicand = (-b).mul_add(b, (-a).mul_add(a, 1.0));
    if !radicand.is_finite() {
        return Err(Error::MathDomain { value: radicand });
    }
    Ok(radicand.max(0.0).sqrt())
}

/// Direct Julian-century argument for long-term ERFA evaluation.
/// Non-finite input is an argument error; overflow is an arithmetic error.
fn long_term_centuries(epj: f64) -> Result<f64> {
    if !epj.is_finite() {
        return Err(Error::NonFiniteArgument { value: epj });
    }
    let t = (epj - 2000.0) / 100.0;
    if !t.is_finite() {
        return Err(Error::MathDomain { value: t });
    }
    Ok(t)
}

/// Refuse an overflowing pole evaluation instead of returning NaN/Inf.
/// Components are dimensionless and otherwise retain the model's result.
fn checked_pole(pole: [f64; 3]) -> Result<[f64; 3]> {
    if let Some(&value) = pole.iter().find(|v| !v.is_finite()) {
        return Err(Error::MathDomain { value });
    }
    Ok(pole)
}

/// J2000 mean-equatorial to mean-equatorial-of-date precession matrix.
///
/// The row-major matrix acts on a Cartesian column vector. `epj` is a
/// TT Julian epoch; the equinox row is the normalized cross product of
/// the equator and ecliptic poles, the middle row is equator cross
/// equinox, and the final row is the equator pole. Normalize the equinox
/// with its Euclidean norm.
/// Model validity and physical uncertainty follow [`ltpecl`]. An undefined
/// equinox is `SingularMatrix`; arithmetic overflow is `MathDomain`.
pub fn ltp_matrix(epj: f64) -> Result<[[f64; 3]; 3]> {
    let equator = ltpequ(epj)?;
    let ecliptic = ltpecl(epj)?;
    let v = cross_product(equator, ecliptic);
    let norm = v[2].mul_add(v[2], v[0].mul_add(v[0], v[1] * v[1])).sqrt();
    if !norm.is_finite() {
        return Err(Error::MathDomain { value: norm });
    }
    if norm == 0.0 {
        return Err(Error::SingularMatrix);
    }
    let inverse = 1.0 / norm;
    let equinox = checked_pole(v.map(|c| inverse * c))?;
    let middle = checked_pole(cross_product(equator, equinox))?;
    Ok([equinox, middle, equator])
}

/// ICRS to mean-equatorial-of-date matrix including first-order frame bias.
///
/// Apply the IERS 2010 bias angles `dx = -16.617 mas`,
/// `de = -6.8192 mas`, `dr = -14.6 mas` to [`ltp_matrix`] in its
/// explicit ERFA row-update order. For a row `(x, y, z)`, the biased
/// row is `(x - y·dr + z·dx, x·dr + y + z·de, -x·dx - y·de + z)`,
/// with bias angles converted to radians. This is a first-order formulation,
/// not a full three-dimensional bias rotation. TT epoch, output layout,
/// validity and errors follow [`ltp_matrix`]; nutation is not included.
pub(crate) fn ltpb_matrix(epj: f64) -> Result<[[f64; 3]; 3]> {
    apply_long_term_frame_bias(ltp_matrix(epj)?)
}

/// Apply the fixed IERS first-order bias to matrix rows or their TT tangents.
///
/// The ordered linear map, angle realization and finite checks are shared
/// with [`ltpb_matrix`]. Matrices are dimensionless, or per TT day when this
/// maps their derivatives. Nonfinite resulting components return an error.
pub(crate) fn apply_long_term_frame_bias(rp: [[f64; 3]; 3]) -> Result<[[f64; 3]; 3]> {
    let dx = -0.016617 * DAS2R;
    let de = -0.0068192 * DAS2R;
    let dr = -0.0146 * DAS2R;
    let mut rpb = [[0.0; 3]; 3];
    for (out, row) in rpb.iter_mut().zip(rp) {
        *out = checked_pole([
            row[2].mul_add(dx, (-row[1]).mul_add(dr, row[0])),
            row[2].mul_add(de, row[0].mul_add(dr, row[1])),
            (-row[0]).mul_add(dx, -row[1] * de) + row[2],
        ])?;
    }
    Ok(rpb)
}

/// Ordered, fused Cartesian cross product of dimensionless pole vectors.
/// The consuming matrix construction checks its resulting components.
fn cross_product(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1].mul_add(b[2], -a[2] * b[1]),
        a[2].mul_add(b[0], -a[0] * b[2]),
        a[0].mul_add(b[1], -a[1] * b[0]),
    ]
}

/// Vondrák 2011 mean ecliptic frame at scalar TT Julian Day `jd_tt`.
///
/// Convert TT to the Julian epoch `2000 + (jd_tt - 2451545) / 365.25`,
/// then acquire the dimensionless ecliptic and equator poles through
/// [`ltpecl`] and [`ltpequ`] at that epoch. The returned four rows are
/// `(x, y, z, equinox)`: `x` is the normalized cross product of the
/// ecliptic and equator poles, `z` is the ecliptic pole, `y = z × x`,
/// and `equinox = x`. The first three rows map a GCRS column vector
/// to the mean ecliptic of date; nutation is not included.
///
/// Model validity and physical uncertainty follow [`ltpecl`], including
/// its published ±200,000-year span about J2000. Finite extrapolation
/// does not extend that validity. Non-finite TT is `NonFiniteArgument`,
/// overflowing pole evaluation is `MathDomain`, and an undefined
/// equinox is `SingularMatrix`.
pub fn ltp_ecliptic_frame(jd_tt: f64) -> Result<[[f64; 3]; 4]> {
    let epj = julian_epoch(jd_tt)?;
    let pecl = ltpecl(epj)?;
    let pequ = ltpequ(epj)?;
    let mut x = [
        pecl[1] * pequ[2] - pecl[2] * pequ[1],
        pecl[2] * pequ[0] - pecl[0] * pequ[2],
        pecl[0] * pequ[1] - pecl[1] * pequ[0],
    ];
    let n = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
    if n == 0.0 {
        return Err(Error::SingularMatrix);
    }
    x = [x[0] / n, x[1] / n, x[2] / n];
    let y = [
        pecl[1] * x[2] - pecl[2] * x[1],
        pecl[2] * x[0] - pecl[0] * x[2],
        pecl[0] * x[1] - pecl[1] * x[0],
    ];
    Ok([x, y, pecl, x])
}

/// Internal alias for the checked physical frame [`ltp_ecliptic_frame`].
///
/// TT is a finite Julian day, the row layout and units are unchanged, and
/// invalid extrapolation returns the same native arithmetic error.
pub(crate) fn ltp_ecliptic_frame_ieee(jd_tt: f64) -> Result<[[f64; 3]; 4]> {
    ltp_ecliptic_frame(jd_tt)
}

/// CIO locator `s` in radians (native ERFA `s06`): the `s+XY/2`
/// series over `(fal 03, falp 03, faf 03, fad 03, faom 03, fave 03, fae 03,
/// fapa 03)` plus the 6-coefficient polynomial, minus `x*y/2`.
pub(super) fn cio_locator_s(t: f64, x: f64, y: f64) -> f64 {
    let fa = [
        fa_poly_rad(t, FA_L03),
        fa_poly_rad(t, FA_LP03),
        fa_poly_rad(t, FA_F03),
        fa_poly_rad(t, FA_D03),
        fa_poly_rad(t, FA_OM03),
        1021.3285546211f64.mul_add(t, 3.176146697) % D2PI,
        628.3075849991f64.mul_add(t, 1.753470314) % D2PI,
        0.00000538691f64.mul_add(t, 0.024381750) * t,
    ];
    let series = [
        S06_S0_ROWS.as_slice(),
        S06_S1_ROWS.as_slice(),
        S06_S2_ROWS.as_slice(),
        S06_S3_ROWS.as_slice(),
        S06_S4_ROWS.as_slice(),
    ];
    let mut w = [
        S06_SP[0], S06_SP[1], S06_SP[2], S06_SP[3], S06_SP[4], S06_SP[5],
    ];
    for (k, rows) in series.iter().enumerate() {
        // Sum the periodic table from its smallest terms, using fused phases.
        for (mult, coeff) in rows.iter().rev() {
            let mut a = 0.0;
            for j in 0..8 {
                a = (mult[j] as f64).mul_add(fa[j], a);
            }
            let (s_a, c_a) = era_sincos(a);
            w[k] += coeff[0].mul_add(s_a, coeff[1] * c_a);
        }
    }
    w[5].mul_add(t, w[4])
        .mul_add(t, w[3])
        .mul_add(t, w[2])
        .mul_add(t, w[1])
        .mul_add(t, w[0])
        * DAS2R
        - x * y / 2.0
}

/// TIO locator `s'` in radians (native ERFA `sp00`): `-47e-6 * t` arcsec.
fn tio_locator_sp(t: f64) -> f64 {
    -47e-6 * t * DAS2R
}

/// Apply an ERFA `R1` pre-rotation in place (`r <- R1(a) @ r`, with the
/// reference row-operation sequence; the column-index loop updates two
/// rows in lockstep, which the iterator lint cannot express).
#[allow(clippy::needless_range_loop)]
fn apply_rx(a: f64, r: &mut [[f64; 3]; 3]) {
    let (s, c) = era_sincos(a);
    for col in 0..3 {
        let (r1, r2) = (r[1][col], r[2][col]);
        r[1][col] = c * r1 + s * r2;
        r[2][col] = -s * r1 + c * r2;
    }
}

/// Mean-to-true equatorial nutation rotation for angles in radians.
///
/// Starting from the identity, apply `R1(eps)`, `R3(-dpsi)` and
/// `R1(-(eps + deps))` in sequence. Keeping the row updates and their
/// fused products separate preserves cancellation in the small nutation
/// terms; regrouping the three matrices is not numerically equivalent.
/// Inputs and the resulting matrix are checked by the consuming boundary.
/// The paired trigonometric evaluation follows the native ERFA convention.
pub(crate) fn numat_rotation(eps: f64, dpsi: f64, deps: f64) -> [[f64; 3]; 3] {
    let mut r = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    apply_rx(eps, &mut r);
    apply_rz(-dpsi, &mut r);
    let (s, c) = era_sincos(-(eps + deps));
    let (_, rows) = r.split_at_mut(1);
    let (row1, row2) = rows.split_at_mut(1);
    for (out1, out2) in row1[0].iter_mut().zip(row2[0].iter_mut()) {
        let (r1, r2) = (*out1, *out2);
        *out1 = c.mul_add(r1, s * r2);
        *out2 = (-s).mul_add(r1, c * r2);
    }
    r
}

/// Apply an ERFA `R2` pre-rotation in place (reference row-operation sequence).
#[allow(clippy::needless_range_loop)]
fn apply_ry(a: f64, r: &mut [[f64; 3]; 3]) {
    let (s, c) = era_sincos(a);
    for col in 0..3 {
        let (r0, r2) = (r[0][col], r[2][col]);
        r[0][col] = c * r0 - s * r2;
        r[2][col] = s * r0 + c * r2;
    }
}

/// Apply an ERFA `R3` pre-rotation in place (reference row-operation sequence).
#[allow(clippy::needless_range_loop)]
fn apply_rz(a: f64, r: &mut [[f64; 3]; 3]) {
    let (s, c) = era_sincos(a);
    for col in 0..3 {
        let (r0, r1) = (r[0][col], r[1][col]);
        r[0][col] = c * r0 + s * r1;
        r[1][col] = -s * r0 + c * r1;
    }
}

/// Celestial-to-intermediate matrix from dimensionless GCRS CIP `x, y`
/// and radian CIO locator `s`: `Rz(-(e+s)) @ Ry(d) @ Rz(e)`, where
/// `r2 = x²+y²`, `e = atan2(y,x)` (zero when r 2 is zero), and
/// `d = atan(sqrt(r2/(1-r2)))`. Use fused r 2/final mixed-row evaluation.
/// The consuming frame boundary supplies a physical pole with r 2 < 1 and
/// refuses non-finite output; entries are dimensionless, row-major.
fn c2ixys(x: f64, y: f64, s: f64) -> [[f64; 3]; 3] {
    let r2 = x.mul_add(x, y * y);
    let e = if r2 > 0.0 { y.atan2(x) } else { 0.0 };
    let d = (r2 / (1.0 - r2)).sqrt().atan();
    let mut r = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    apply_rz(e, &mut r);
    apply_ry(d, &mut r);
    apply_rz_fused(-(e + s), &mut r);
    r
}

/// Fukushima-Williams bias/precession/nutation rotation from radian angles.
///
/// Starting at identity, pre-rotate by R 3(gamb), R 1(phib), R 3(-psi), then
/// R 1(-eps). Preserve paired trigonometry and fused row updates rather than
/// regrouping products; the last two rotations mix nonzero components and
/// expose cancellation in the small CIP entries. The consuming public
/// boundary owns input/output validation. Matrix entries are dimensionless.
pub(crate) fn fw2m_rotation(gamb: f64, phib: f64, psi: f64, eps: f64) -> [[f64; 3]; 3] {
    let mut r = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    apply_rz(gamb, &mut r);
    apply_rx(phib, &mut r);
    apply_rz_fused(-psi, &mut r);
    let (s, c) = era_sincos(-eps);
    let [_, row1, row2] = r;
    for (col, (r1, r2)) in row1.into_iter().zip(row2).enumerate() {
        r[1][col] = c.mul_add(r1, s * r2);
        r[2][col] = (-s).mul_add(r1, c * r2);
    }
    r
}

/// Fused R 3 pre-rotation of a dimensionless row-major matrix.
/// Preserve the first-product fusion in each of the two mixed rows;
/// separate-product rotation helpers retain their independent convention.
fn apply_rz_fused(angle: f64, matrix: &mut [[f64; 3]; 3]) {
    let (s, c) = era_sincos(angle);
    let [row0, row1, _] = *matrix;
    for (col, (r0, r1)) in row0.into_iter().zip(row1).enumerate() {
        matrix[0][col] = c.mul_add(r0, s * r1);
        matrix[1][col] = (-s).mul_add(r0, c * r1);
    }
}

/// Row-major product of two dimensionless rotations, using fused dot products.
/// A consuming checked boundary validates the resulting components.
fn rotation_product(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            a[i][2].mul_add(b[2][j], a[i][1].mul_add(b[1][j], a[i][0] * b[0][j]))
        })
    })
}

/// Polar-motion matrix (native ERFA `pom00`): `Rx(-yp) @ Ry(-xp) @
/// Rz(sp)` applied to the identity, in that order.
fn pom00(xp: f64, yp: f64, sp: f64) -> [[f64; 3]; 3] {
    let mut r = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    apply_rz(sp, &mut r);
    apply_ry(-xp, &mut r);
    apply_rx(-yp, &mut r);
    r
}

/// IAU 2006/2000 A celestial-to-terrestrial orientation.
///
/// Return a dimensionless row-major GCRS-to-ITRF rotation. TT drives
/// bias/precession/nutation and the CIO locator; UT1 drives Earth rotation.
/// Polar-motion angles `xp` and `yp` are radians. The TIO locator and polar
/// motion complete the IERS orientation chain. J2000 day offsets retain the
/// available scalar-date resolution through the common native model.
///
/// The IAU 2006/2000 A model is most accurate within a few centuries of J2000;
/// the result excludes observed free-core nutation and other unprovided EOP.
/// Non-finite inputs return [`Error::NonFiniteArgument`]; arithmetic overflow
/// returns [`Error::MathDomain`].
pub fn c2t06a(jd_tt: f64, jd_ut1: f64, xp: f64, yp: f64) -> Result<[[f64; 3]; 3]> {
    for value in [jd_tt, jd_ut1, xp, yp] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    let result = c2t06a_from_offsets(jd_tt - DJ00, jd_ut1 - DJ00, xp, yp);
    for &value in result.iter().flatten() {
        if !value.is_finite() {
            return Err(Error::MathDomain { value });
        }
    }
    Ok(result)
}

/// IAU 2006/2000 A orientation from TT and UT1 J2000 offsets, in days.
///
/// This zero-polar-motion entry evaluates the same native model as [`c2t06a`]
/// without constructing a scalar Julian day. The dimensionless row-major
/// result maps GCRS column vectors to ITRF. Callers own finite-input admission;
/// invalid arithmetic produces non-finite components. Model accuracy and
/// omitted EOP are those documented for [`c2t06a`].
pub fn source_c2t06a_from_j2000_offsets(tt_offset: f64, ut_offset: f64) -> [[f64; 3]; 3] {
    c2t06a_from_offsets(tt_offset, ut_offset, 0.0, 0.0)
}

/// Common IAU orientation model for finite scalar dates or retained offsets.
fn c2t06a_from_offsets(tt_offset: f64, ut_offset: f64, xp: f64, yp: f64) -> [[f64; 3]; 3] {
    let t = tt_offset / DJC;
    let (gamb, phib, psib, epsa) = pfw06_centuries(t);
    let (dpsi, deps) = source_nut06a_from_j2000_offset(tt_offset);
    let rbpn = fw2m_rotation(gamb, phib, psib + dpsi, epsa + deps);
    let (rpom, rotation) = c2t06a_parts(t, rbpn, ut_offset, xp, yp);
    rotation_product(rpom, rotation)
}

/// CIO, Earth rotation and polar-motion assembly. All angular units are
/// radians, UT1 is the independent J2000 day offset, and matrices are row-major.
fn c2t06a_parts(
    t: f64,
    rbpn: [[f64; 3]; 3],
    ut_offset: f64,
    xp: f64,
    yp: f64,
) -> ([[f64; 3]; 3], [[f64; 3]; 3]) {
    let (x, y) = (rbpn[2][0], rbpn[2][1]);
    let s = cio_locator_s(t, x, y);
    let rc2i = c2ixys(x, y, s);
    let era = positions::source_earth_rotation_angle_from_j2000_offset(ut_offset);
    let sp = tio_locator_sp(t);
    let rpom = pom00(xp, yp, sp);
    // c 2 tcio: `r = rc2i; Rz(era, r); rc2t = rpom @ r`.
    let mut r = rc2i;
    apply_rz_fused(era, &mut r);
    (rpom, r)
}

// ---------------------------------------------------------------------------
// Epoch converters
// ---------------------------------------------------------------------------

/// Julian epoch to two-part Julian Day (native ERFA `epj2jd`):
/// `(DJM0, DJM00 + (epj - 2000) * 365.25)` as one fused multiply-add
/// (single rounding for the combined operation).
///
/// The returned parts are additive days, not an integer day and a fraction:
/// the first is 2400000.5 and the second is an MJD-sized offset. This defines a
/// Julian epoch coordinate, not a UTC/TT/UT1 conversion or leap-second lookup.
/// Finite input is required; finite extreme inputs can overflow the result.
/// Nonfinite input returns [`Error::NonFiniteArgument`].
pub fn epj2jd(epj: f64) -> Result<(f64, f64)> {
    if !epj.is_finite() {
        return Err(Error::NonFiniteArgument { value: epj });
    }
    Ok((DJM0, (epj - 2000.0).mul_add(365.25, DJM00)))
}

/// Besselian epoch to two-part Julian Day (native ERFA `epb2jd`):
/// `(DJM0, 15019.81352 + (epb - 1900) * DTY)` as one fused
/// multiply-add (separate multiply-then-add is 1 ulp off at some
/// epochs, e.g. epb 1800).
///
/// Returns additive days anchored at 2400000.5; neither part is required to be
/// an integer or a fraction in `[0, 1)`. This is an epoch-coordinate definition,
/// not a time-scale conversion. Finite extreme inputs can overflow the result;
/// nonfinite input returns [`Error::NonFiniteArgument`].
pub fn epb2jd(epb: f64) -> Result<(f64, f64)> {
    if !epb.is_finite() {
        return Err(Error::NonFiniteArgument { value: epb });
    }
    Ok((DJM0, (epb - 1900.0).mul_add(DTY, 15019.81352)))
}

/// Gregorian calendar date to two-part Julian Day (native ERFA
/// `cal2jd`).
///
/// Returns `(djm0, djm, status)` with the ERFA status contract:
/// `0` ok (day out of range still returns a date with `-3`), `-1`
/// year below 4800 BC, `-2` bad month, `-3` bad day for the
/// month/year. Integer division truncates toward zero exactly as the
/// proleptic Gregorian formula. Invalid days still return the corresponding
/// arithmetic date with status `-3`; year and month errors return zero MJD.
/// Dates are at midnight with astronomical year numbering (0 is 1 BCE).
/// The returned additive days are `(2400000.5, MJD)`, not normalized
/// integer-day/fraction parts. No time scale or leap-second conversion is
/// implied; callers must inspect `status` rather than expecting `Result`.
#[allow(clippy::too_many_arguments)]
pub fn cal2jd(iy: i32, im: i32, id: i32) -> (f64, f64, i32) {
    const IYMIN: i32 = -4799;
    const MTAB: [i32; 12] = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    if iy < IYMIN {
        return (DJM0, 0.0, -1);
    }
    if !(1..=12).contains(&im) {
        return (DJM0, 0.0, -2);
    }
    let leap = im == 2 && iy % 4 == 0 && (iy % 100 != 0 || iy % 400 == 0);
    let mut status = 0;
    if id < 1 || id > MTAB[(im - 1) as usize] + i32::from(leap) {
        status = -3;
    }
    // C `my = (im - 14) / 12` truncates toward zero; Rust `/` on i32
    // and the i 64 widening below do the same.
    let my = (im - 14) / 12;
    let iypmy = iy as i64 + my as i64;
    let djm = (1461 * (iypmy + 4800) / 4 + 367 * (im as i64 - 2 - 12 * my as i64) / 12
        - 3 * ((iypmy + 4900) / 100) / 4
        + id as i64
        - 2432076) as f64;
    (DJM0, djm, status)
}

/// Dimensionless GCRS-to-mean-J2000 frame-bias rotation.
///
/// Return the fixed IERS bias angles encoded by [`FRAME_BIAS_J2000`]. This
/// epoch-independent matrix is infallible; it includes neither precession
/// nor nutation.
pub fn bp06_bias() -> [[f64; 3]; 3] {
    FRAME_BIAS_J2000
}

#[cfg(test)]
mod tests {
    #[test]
    fn obliquity_nutation_compatibility_words() {
        // Existing native value/obliquity words are a compatibility contract;
        // the independent series-derivative check below supplies the oracle.
        for (tt, expected) in [
            (
                2415020.5,
                [
                    0x3f16281a684df145_u64,
                    0xbee748df8812e491,
                    0x3e6f07e475c28edd,
                ],
            ),
            (
                2451545.0,
                [0xbf11b4d17b668d56, 0xbefd545ce227381e, 0xbe7ab9d1a6f459f7],
            ),
            (
                2460310.5,
                [0xbefb3e629128cfa7, 0x3f04818633bb7d63, 0x3e8dc04b2af230c5],
            ),
        ] {
            let angles = super::nut06a(tt).unwrap();
            assert_eq!(
                [
                    angles.0.to_bits(),
                    angles.1.to_bits(),
                    super::nut06a_obliquity_rate(tt).unwrap().to_bits()
                ],
                expected
            );
        }
    }

    #[test]
    fn analytic_nutation_rates_match_independent_pyerfa_derivatives() {
        // PyERFA2.0.1.5 nut06a positions, differentiated independently with
        // eight centered samples and exact Decimal weights at each fixed
        // h=1/32,1/64,1/16 TT day. These are series references, not native
        // rates or the former one-day secant used as expected values.
        for (tt, references) in [
            (
                2415020.5,
                [
                    [8.42453930209509e-07, 5.779942456757993e-08],
                    [8.424539300992064e-07, 5.779942455997948e-08],
                    [8.424539304020656e-07, 5.7799424597802135e-08],
                ],
            ),
            (
                2451545.0,
                [
                    [3.5641518015717425e-08, -9.956156942332037e-08],
                    [3.564151801288662e-08, -9.956156942306043e-08],
                    [3.5641518014409156e-08, -9.956156942313122e-08],
                ],
            ),
            (
                2460310.5,
                [
                    [-9.821786275742216e-08, 2.2166331833274136e-07],
                    [-9.821786281456555e-08, 2.2166331832489768e-07],
                    [-9.821786271463721e-08, 2.2166331834374781e-07],
                ],
            ),
        ] {
            let actual = super::nut06a_rates(tt).unwrap();
            for expected in references {
                for (value, oracle) in [(actual.0, expected[0]), (actual.1, expected[1])] {
                    // PLAN's absolute angular-speed floor, degrees/day.
                    assert!(
                        (value - oracle).abs().to_degrees() <= 1e-12,
                        "{tt}: analytic {value}, independent {oracle}"
                    );
                }
            }
            assert_eq!(
                super::source_nut06a_longitude_rate(tt).to_bits(),
                actual.0.to_bits()
            );
        }
    }

    #[test]
    fn analytic_nutation_rates_preserve_checked_and_internal_clock_contracts() {
        for tt in [f64::NAN, f64::NEG_INFINITY, f64::INFINITY] {
            assert!(matches!(super::nut06a_longitude_rate(tt),
                Err(crate::errors::Error::NonFiniteArgument { value })
                if value.to_bits() == tt.to_bits()));
            assert!(matches!(super::nut06a_obliquity_rate(tt),
                Err(crate::errors::Error::NonFiniteArgument { value })
                if value.to_bits() == tt.to_bits()));
            assert!(super::source_nut06a_longitude_rate(tt).is_nan());
        }
    }

    #[test]
    fn analytic_obliquity_nutation_rate_matches_position_derivative() {
        let h = 1.0 / 32.0;
        for tt in [2378496.5, 2451545.0, 2460310.5, 2524593.5] {
            let value = |offset| super::nut06a(tt + offset).unwrap().1;
            let sixth = (45.0 * (value(h) - value(-h)) - 9.0 * (value(2.0 * h) - value(-2.0 * h))
                + value(3.0 * h)
                - value(-3.0 * h))
                / (60.0 * h);
            let analytic = super::nut06a_obliquity_rate(tt).unwrap();
            // Absolute rad/day bound stricter than PLAN's relative daily
            // angular-speed criterion at these nonsingular obliquity rates.
            assert!(
                (analytic - sixth).abs() <= 1e-13,
                "{tt}: {analytic} vs {sixth}"
            );
        }
        assert!(super::nut06a_obliquity_rate(f64::NAN).is_err());
    }

    use super::*;

    // PLAN's unchanged angular comparison limit, expressed in radians.
    const ANGLE_LIMIT_RAD: f64 = 0.001 * DAS2R;

    fn assert_angle_matches(actual: f64, expected: f64) {
        assert!(actual.is_finite() && expected.is_finite());
        assert!(
            (actual - expected).abs() <= ANGLE_LIMIT_RAD,
            "angle {actual} vs reference {expected}"
        );
    }

    fn assert_direction_matches(actual: [f64; 3], expected: [f64; 3]) {
        let cross = [
            actual[1] * expected[2] - actual[2] * expected[1],
            actual[2] * expected[0] - actual[0] * expected[2],
            actual[0] * expected[1] - actual[1] * expected[0],
        ];
        let cross_norm = cross.iter().map(|x| x * x).sum::<f64>().sqrt();
        let dot = actual.iter().zip(expected).map(|(a, b)| a * b).sum::<f64>();
        assert!(
            cross_norm.atan2(dot).abs() <= ANGLE_LIMIT_RAD,
            "direction {actual:?} vs reference {expected:?}"
        );
        let norm2 = actual.iter().map(|x| x * x).sum::<f64>();
        assert!(
            (norm2 - 1.0).abs() <= 256.0 * f64::EPSILON,
            "unit direction norm squared {norm2}"
        );
    }

    fn assert_matrix_matches<const N: usize>(actual: [[f64; 3]; N], expected: [[f64; 3]; N]) {
        for (row, reference) in actual.into_iter().zip(expected) {
            assert_direction_matches(row, reference);
        }
    }

    fn assert_rotation(matrix: [[f64; 3]; 3]) {
        for i in 0..3 {
            for j in 0..3 {
                let dot = matrix[i]
                    .iter()
                    .zip(matrix[j])
                    .map(|(a, b)| a * b)
                    .sum::<f64>();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (dot - expected).abs() <= 256.0 * f64::EPSILON,
                    "rotation row dot {i},{j}: {dot}"
                );
            }
        }
        let handedness = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
            - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
            + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
        assert!((handedness - 1.0).abs() <= 256.0 * f64::EPSILON);
    }

    #[test]
    fn paired_trigonometry_preserves_the_unit_circle() {
        for angle in [
            -1_000_000.0,
            -std::f64::consts::PI,
            -0.25,
            0.0,
            0.25,
            std::f64::consts::FRAC_PI_2,
            1_000_000.0,
        ] {
            let (sine, cosine) = era_sincos(angle);
            assert!((sine * sine + cosine * cosine - 1.0).abs() <= 4.0 * f64::EPSILON);
        }
        let (sine, cosine) = era_sincos(std::f64::consts::FRAC_PI_2);
        assert!((sine - 1.0).abs() <= f64::EPSILON);
        assert!(cosine.abs() <= f64::EPSILON);
    }

    #[test]
    fn retained_offset_orientation_matches_scalar_dates_and_is_a_rotation() {
        for (tt, ut) in [(2451545., 2451544.999), (2460000.5, 2460000.499)] {
            let offsets = source_c2t06a_from_j2000_offsets(tt - DJ00, ut - DJ00);
            let scalar = c2t06a(tt, ut, 0., 0.).unwrap();
            assert_matrix_matches(offsets, scalar);
            assert_rotation(offsets);
            assert_rotation(scalar);
        }
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for values in [
                (bad, DJ00, 0.0, 0.0),
                (DJ00, bad, 0.0, 0.0),
                (DJ00, DJ00, bad, 0.0),
                (DJ00, DJ00, 0.0, bad),
            ] {
                assert!(c2t06a(values.0, values.1, values.2, values.3).is_err());
            }
        }
    }

    const J2000: f64 = 2451545.0;

    #[test]
    fn coefficient_series_are_complete() {
        assert_eq!(NUT_LS_ROWS.len() + NUT_PL_ROWS.len(), 1365);
        assert_eq!(NUT00B_ROWS.len(), 77);
        assert_eq!(
            S06_S0_ROWS.len()
                + S06_S1_ROWS.len()
                + S06_S2_ROWS.len()
                + S06_S3_ROWS.len()
                + S06_S4_ROWS.len(),
            33 + 3 + 25 + 4 + 1
        );
        assert_eq!(S06_SP.len(), 6);
    }

    #[test]
    fn nutation_model_info_matches_shipped_tables() {
        let info = nutation_model_info();
        assert_eq!(info.model, "IAU2006_2000A");
        assert_eq!(info.terms as usize, NUT_LS_ROWS.len() + NUT_PL_ROWS.len());
        assert_eq!(info.terms, 1365);
        assert_eq!(info.precision, "~0.01-0.05 mas");
        assert_eq!(info.source, "native");
    }

    #[test]
    fn nutation_at_j2000_matches_erfa_angles() {
        // Independent ERFA 2.0.1 numerical values.
        let (dpsi, deps) = nut06a(J2000).unwrap();
        assert_angle_matches(dpsi, -6.754425598969513e-05);
        assert_angle_matches(deps, -2.7970831192374137e-05);
        let (dpsi_a, deps_a) = nut00a(J2000).unwrap();
        assert_angle_matches(dpsi_a, -6.754422426417299e-05);
        assert_angle_matches(deps_a, -2.7970831192374137e-05);
        let (dpsi_b, deps_b) = nut00b(J2000).unwrap();
        assert_angle_matches(dpsi_b, -6.754261253992235e-05);
        assert_angle_matches(deps_b, -2.7970923310985653e-05);
    }

    #[test]
    fn general_precession_matches_erfa_angles() {
        // ERFA p 06 e slot 12 is general precession. Independent degree values
        // test the polynomial across separated epochs.
        for (jd, expected_deg) in [
            (2460000.5f64, 0.32339455534432754f64),
            (2400000.5f64, -1.9706927337282851f64),
            (2500000.0f64, 1.8536876027159588f64),
        ] {
            let got = positions::general_precession_deg(jd).unwrap();
            assert_angle_matches(got.to_radians(), expected_deg.to_radians());
        }
    }

    #[test]
    fn fundamental_arguments_at_j2000_follow_iers_constants() {
        assert_angle_matches(fal03(0.0).unwrap(), 485868.249036 * DAS2R);
        assert_angle_matches(falp03(0.0).unwrap(), 1287104.793048 * DAS2R);
        assert_eq!(fapa03(0.0).unwrap(), 0.0);
    }

    #[test]
    fn epoch_converters_match_erfa_identities() {
        // erfa.epj 2 jd(2000) = (DJM 0, DJM 00); epb 2 jd(1900) = (DJM 0, 15019.81352).
        assert_eq!(epj2jd(2000.0).unwrap(), (2400000.5, 51544.5));
        assert_eq!(epb2jd(1900.0).unwrap(), (2400000.5, 15019.81352));
        assert_eq!(julian_epoch_to_jd(2000.0).unwrap(), J2000);
        assert_eq!(julian_epoch(J2000).unwrap(), 2000.0);
        // cal 2 jd: J2000.0 noon TT is 2000-01-01 12:00 -> MJD 51544.0.
        assert_eq!(cal2jd(2000, 1, 1), (2400000.5, 51544.0, 0));
        assert_eq!(cal2jd(-4800, 1, 1).2, -1);
        assert_eq!(cal2jd(2000, 13, 1).2, -2);
        // Day out of range still returns a date with status -3 (ERFA contract).
        let (m0, m, s) = cal2jd(2001, 2, 29);
        assert_eq!((m0, s), (2400000.5, -3));
        assert!(m.is_finite());
        assert!(cal2jd(2000, 2, 29).2 == 0);
    }

    #[test]
    fn long_term_poles_match_erfa_directions() {
        // ERFA 2.0.1: J2000, its immediate
        // epoch neighbor, a modern epoch and both model-span endpoints.
        // Reference directions are retained independently of arithmetic order.
        for (epj, ecliptic, equator) in [
            (
                2000.0,
                [0xbc5455a5b2ff8f9d, 0xbfd9752d885c3387, 0x3fed5c0382e42276],
                [0xbc6455a5b2ff8f9d, 0x3d9552a345878d27, 0x3ff0000000000000],
            ),
            (
                f64::from_bits(2000.0_f64.to_bits() + 1),
                [0xbc545731d40162b4, 0xbfd9752d885c3387, 0x3fed5c0382e42276],
                [0x3c745707eb7b89ff, 0x3d9552a344e7ba85, 0x3ff0000000000000],
            ),
            (
                2065.125,
                [0x3eeca3d69f26defa, 0xbfd972f52975b618, 0x3fed5c7eb20cffd1],
                [0x3f79ea1d0c8a53f6, 0xbf083277d5543a66, 0x3fefffd60655bf21],
            ),
            (
                -198000.0,
                [0x3f72858571246fc2, 0xbfd829cda5bcdc0c, 0x3feda19bb060beea],
                [0x3fd7f3e3a77349ab, 0xbfcb3024763f5879, 0x3fece2ab28c52f5e],
            ),
            (
                202000.0,
                [0x3f93ba4255fd3759, 0xbfda5ea2ebfb94ec, 0x3fed26a710111e68],
                [0xbfd81e1459f855f9, 0xbfd6652088f1b919, 0x3feb71ccc0a32f3b],
            ),
        ] {
            assert_direction_matches(ltpecl(epj).unwrap(), ecliptic.map(f64::from_bits));
            assert_direction_matches(ltpequ(epj).unwrap(), equator.map(f64::from_bits));
        }
    }

    #[test]
    fn long_term_precession_matches_erfa_rotations() {
        // ERFA 2.0.1, including normalization and first-order ICRS bias.
        for (epj, bias, expected) in [
            (
                2000.0,
                false,
                [
                    [0x3ff0000000000000, 0x3c655786a1b6e56a, 0x3c6455a5b2ff1dd8],
                    [0xbc655786a1b67904, 0x3ff0000000000000, 0xbd9552a345878d27],
                    [0xbc6455a5b2ff8f9d, 0x3d9552a345878d27, 0x3ff0000000000000],
                ],
            ),
            (
                2000.0,
                true,
                [
                    [0x3ff0000000000000, 0xbe73002824246eb5, 0x3e75a024b5a30599],
                    [0x3e7300282424e1fd, 0x3ff0000000000000, 0x3e61bf1f581fb721],
                    [0xbe75a024b5a2a04f, 0xbe61bf1f581fb720, 0x3ff0000000000000],
                ],
            ),
            (
                2065.125,
                false,
                [
                    [0x3feffef7955c9d4d, 0xbf8dd3f742337071, 0xbf79ea1d5906c33e],
                    [0x3f8dd3f752d04d1b, 0x3fefff218f06ddf7, 0xbf081df49ad4b05f],
                    [0x3f79ea1d0c8a53f6, 0xbf083277d5543a66, 0x3fefffd60655bf21],
                ],
            ),
            (
                2065.125,
                true,
                [
                    [0x3feffef79514fa99, 0xbf8dd400bac929fb, 0xbf79ea07daab8064],
                    [0x3f8dd400d2c2edd1, 0x3fefff218e796377, 0xbf08195c7670f915],
                    [0x3f79ea076c4883ef, 0xbf0836f7251f3ce2, 0x3fefffd6069b96bf],
                ],
            ),
        ] {
            let actual = if bias {
                ltpb_matrix(epj)
            } else {
                ltp_matrix(epj)
            }
            .unwrap();
            assert_matrix_matches(actual, expected.map(|row| row.map(f64::from_bits)));
            assert_rotation(actual);
        }
    }

    #[test]
    fn fw_and_cio_rotations_match_erfa_orientations() {
        // ERFA 2.0.1 at a fixed TT/UT1 epoch pair. Matrix rows represent
        // physical directions; compare them with the PLAN angular limit.
        let tt = f64::from_bits(0x4142e9aca010624e);
        let ut = f64::from_bits(0x4142e9ac9fd9189c);
        let expected_pnm = [
            [0x3feffea438729d75, 0xbf911aa8f5199574, 0xbf7db7ea866fa25e],
            [0x3f911ab34f851eba, 0x3feffedb6bc2a68f, 0xbf049fb5270666af],
            [0x3f7db78b2ef646e5, 0xbf15759697e95240, 0x3fefffc8cc71d3b8],
        ];
        let expected_c2t = [
            [0x3fefaca3ce7f412d, 0x3fc230d21160c65e, 0xbf7d5e21476ebc00],
            [0xbfc230adc3bd0a96, 0x3fefacda41bea362, 0x3f52386fac367d7c],
            [0x3f7db78b2ef646e5, 0xbf15759697e95241, 0x3fefffc8cc71d3b8],
        ];
        let precession = pnm06a(tt).unwrap();
        let terrestrial = c2t06a(tt, ut, 0.0, 0.0).unwrap();
        assert_matrix_matches(precession, expected_pnm.map(|row| row.map(f64::from_bits)));
        assert_matrix_matches(terrestrial, expected_c2t.map(|row| row.map(f64::from_bits)));
        assert_rotation(precession);
        assert_rotation(terrestrial);
    }

    #[test]
    fn overflowing_long_term_evaluations_are_refused() {
        // ±1 e 60 leaves finite polynomial components but overflows their
        // squared pole radicand; clamping -Inf must not conceal that error.
        for epoch in [f64::MAX, -f64::MAX, 1e60, -1e60] {
            assert!(matches!(ltpecl(epoch), Err(Error::MathDomain { .. })));
            assert!(matches!(ltpequ(epoch), Err(Error::MathDomain { .. })));
            assert!(ltp_matrix(epoch).is_err());
            assert!(ltpb_matrix(epoch).is_err());
            assert!(matches!(
                ltp_ecliptic_frame(epoch),
                Err(Error::MathDomain { .. })
            ));
        }
    }

    #[test]
    fn non_finite_inputs_are_refused() {
        for bad in [f64::NAN, f64::INFINITY] {
            assert!(nut00a(bad).is_err());
            assert!(nut00b(bad).is_err());
            assert!(nut06a(bad).is_err());
            assert!(pnm06a(bad).is_err());
            assert!(ltpecl(bad).is_err());
            assert!(ltpequ(bad).is_err());
            assert!(c2t06a(bad, J2000, 0.0, 0.0).is_err());
            assert!(c2t06a(J2000, bad, 0.0, 0.0).is_err());
            assert!(fal03(bad).is_err());
            assert!(fapa03(bad).is_err());
            assert!(epj2jd(bad).is_err());
            assert!(epb2jd(bad).is_err());
            assert!(julian_epoch(bad).is_err());
            assert!(julian_epoch_to_jd(bad).is_err());
        }
    }

    #[test]
    fn pnm06a_at_j2000_is_bias_with_nutation() {
        // At J2000 the bias-precession part is the frame bias; nutation
        // perturbs it at the ~1 e-4 rad level, never the identity.
        let m = pnm06a(J2000).unwrap();
        let diff = (m[0][0] - FRAME_BIAS_J2000[0][0]).abs();
        assert!(diff > 1e-9 && diff < 1e-3, "diff {diff}");
    }

    #[test]
    fn c2t06a_two_part_date2_grid_is_finer() {
        // The offset representation has a finer grid: a 1 e-10-day
        // nudge after splitting resolves here but vanishes on scalar JD.
        // This cannot recover information already rounded from an input JD.
        let jd = 2458843.7463645935f64;
        let date2 = jd - DJ00;
        let step_single = f64::from_bits(jd.to_bits() + 1) - jd;
        let step_two = f64::from_bits(date2.to_bits() + 1) - date2;
        assert!(
            step_single / step_two >= 256.0,
            "grid ratio {}",
            step_single / step_two
        );
        assert_eq!((jd + 1e-10).to_bits(), jd.to_bits());
        assert_ne!((date2 + 1e-10).to_bits(), date2.to_bits());
    }

    #[test]
    #[allow(clippy::excessive_precision)]
    fn terrestrial_orientation_matches_erfa_reference() {
        // Independent ERFA 2.0.1 orientation at the retained TT/UT1 inputs.
        let m = c2t06a(2458843.7463645935, 2458843.7455658740, 0.0, 0.0).unwrap();
        let expected: [[u64; 3]; 3] = [
            [0xbfeff79136537c7b, 0xbfa7347310d20349, 0x3f5f3bd059ff442a],
            [0x3fa7346f6b7192a7, 0xbfeff79506d22a38, 0xbf1a431a480b1e9e],
            [0x3f5f46a442ede1a3, 0xbeecb0170d833be0, 0x3feffffc2dc3a32b],
        ];
        assert_matrix_matches(m, expected.map(|row| row.map(f64::from_bits)));
        assert_rotation(m);
    }

    #[test]
    fn ltp_frame_rows_are_unit_triad() {
        let frame = ltp_ecliptic_frame(2460000.5).unwrap();
        for row in [frame[0], frame[1], frame[2]] {
            let n = (row[0] * row[0] + row[1] * row[1] + row[2] * row[2]).sqrt();
            assert!((n - 1.0).abs() < 1e-12, "n {n}");
        }
    }

    #[test]
    fn ltp_ecliptic_frame_matches_independent_directions() {
        // Recorded precession-frame values at the Fagan defining epoch and
        // three TT samples around UT1 2460400.25. Expected values are external
        // numerical references, never generated by this implementation.
        for (tt_bits, expected) in [
            (
                0x414290813633e7f9,
                [
                    [0xbfefff64377a4fb2, 0xbf86e4a2d43c2ca7, 0xbf73e6bc5938904c],
                    [0x3f88f65056ebae88, 0xbfed5b15f7c9296f, 0xbfd97665feb012b1],
                    [0xbee4da77d4ba30f1, 0xbfd976e26d73c712, 0x3fed5ba4c9b69ecb],
                    [0xbfefff64377a4fb2, 0xbf86e4a2d43c2ca7, 0xbf73e6bc5938904c],
                ],
            ),
            (
                0x4142c5782019dd6a,
                [
                    [0xbfefffdb5c66aec8, 0x3f7634d9bf02aebe, 0x3f634c07d21d8d44],
                    [0xbf7836438820bf92, 0xbfed5c0fcda8c037, 0xbfd9743c9805e841],
                    [0x3ed4eeead2e1150e, 0xbfd97459da1b4476, 0x3fed5c31652accf7],
                    [0xbfefffdb5c66aec8, 0x3f7634d9bf02aebe, 0x3f634c07d21d8d44],
                ],
            ),
            (
                0x4142c578201a3e81,
                [
                    [0xbfefffdb5c66ad2c, 0x3f7634d9bf7f5958, 0x3f634c07d289def2],
                    [0xbf78364388a8abc0, 0xbfed5c0fcda8bfbf, 0xbfd9743c9805e253],
                    [0x3ed4eeead357cb90, 0xbfd97459da1b3fd2, 0x3fed5c31652acdf8],
                    [0xbfefffdb5c66ad2c, 0x3f7634d9bf7f5958, 0x3f634c07d289def2],
                ],
            ),
            (
                0x4142c578201a9f98,
                [
                    [0xbfefffdb5c66ab91, 0x3f7634d9bffc03ef, 0x3f634c07d2f630b4],
                    [0xbf783643893097f2, 0xbfed5c0fcda8bf47, 0xbfd9743c9805dc67],
                    [0x3ed4eeead3ceaabe, 0xbfd97459da1b3b2e, 0x3fed5c31652acefa],
                    [0xbfefffdb5c66ab91, 0x3f7634d9bffc03ef, 0x3f634c07d2f630b4],
                ],
            ),
        ] {
            let frame = ltp_ecliptic_frame(f64::from_bits(tt_bits)).unwrap();
            assert_matrix_matches(frame, expected.map(|row| row.map(f64::from_bits)));
            assert_rotation([frame[0], frame[1], frame[2]]);
        }
    }

    #[test]
    fn ltp_ecliptic_frame_uses_direct_epoch_poles() {
        for tt in [
            J2000 - 200_000.0 * 365.25,
            2433282.42345905,
            J2000,
            f64::from_bits(J2000.to_bits() + 1),
            2460400.250800908,
            J2000 + 200_000.0 * 365.25,
        ] {
            let frame = ltp_ecliptic_frame(tt).unwrap();
            let epj = julian_epoch(tt).unwrap();
            let equator = ltpequ(epj).unwrap();
            assert_direction_matches(frame[2], ltpecl(epj).unwrap());
            assert_direction_matches(frame[0], frame[3]);
            assert_rotation([frame[0], frame[1], frame[2]]);
            for pole in [frame[2], equator] {
                let dot = frame[0][0] * pole[0] + frame[0][1] * pole[1] + frame[0][2] * pole[2];
                assert!(dot.abs() < 1e-14, "equinox orthogonal to pole at TT {tt}");
            }
        }
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                ltp_ecliptic_frame(bad),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
    }
}

#[cfg(test)]
mod offset_resolution_tests {
    use super::*;

    #[test]
    fn retained_offset_nutation_resolves_below_a_scalar_day_ulp() {
        let offset = 2.0_f64.powi(-33);
        assert_eq!(DJ00 + offset, DJ00);
        let retained = source_nut00b_from_j2000_offset(offset);
        let collapsed = nut00b(DJ00 + offset).unwrap();
        assert!(
            retained.0 != collapsed.0 || retained.1 != collapsed.1,
            "the independently retained epoch must affect the nutation series"
        );
    }
}

#[cfg(test)]
mod exact_nutation_cache_tests {
    use super::*;

    #[test]
    fn native_full_series_keeps_its_own_cached_words_and_input_errors() {
        for epoch in [1000000.25, 2415020.5, 2451545.0, 2488069.5, 5000000.75] {
            let expected = nut00a_centuries(tt_centuries(epoch).unwrap());
            for _ in 0..3 {
                let actual = nut00a(epoch).unwrap();
                assert_eq!(
                    [actual.0.to_bits(), actual.1.to_bits()],
                    [expected.0.to_bits(), expected.1.to_bits()]
                );
            }
        }
        assert!(nut00a(f64::NAN).is_err());
        assert!(nut00a(f64::INFINITY).is_err());
    }

    fn assert_angles(actual: (f64, f64), expected: (f64, f64)) {
        for (actual, expected) in [actual.0, actual.1]
            .into_iter()
            .zip([expected.0, expected.1])
        {
            if expected.is_nan() {
                // Arithmetic NaN payloads are not specified by Rust.
                assert!(actual.is_nan());
            } else {
                assert_eq!(actual.to_bits(), expected.to_bits());
            }
        }
    }

    /// The original uncached offset route, including its P03 operation order.
    fn uncached_offset(offset: f64) -> (f64, f64) {
        let t = offset / DJC;
        let (dp, de) = nut00a_centuries(t);
        let fj2 = -2.7774e-6 * t;
        (dp + dp * (0.4697e-6 + fj2), de + de * fj2)
    }

    #[test]
    fn checked_and_offset_nutation_share_only_evaluated_century_words() {
        for jd in [
            0.0,
            f64::from_bits(1),
            1000000.25,
            2415020.5,
            DJ00,
            2488069.5,
            5000000.75,
        ] {
            let offset = jd - DJ00;
            let t = tt_centuries(jd).unwrap();
            let raw = nut00a_centuries(t);
            let adjusted = uncached_offset(offset);
            for source_first in [false, true] {
                if source_first {
                    assert_angles(source_nut06a_from_j2000_offset(offset), adjusted);
                }
                assert_angles(nut00a(jd).unwrap(), raw);
                assert_angles(nut06a(jd).unwrap(), adjusted);
                assert_angles(source_nut06a_from_j2000_offset(offset), adjusted);
                // Other parallel tests may evict this global entry. Any value
                // still retained at the actual t key must be the raw model.
                if let Some(retained) = crate::pure_cache::erfa_nutation(t.to_bits()) {
                    assert_angles(retained, raw);
                }
            }
        }
        // A mixed JD/t namespace would alias these physically different calls:
        // the old checked key for J2000 equals the source route's t below.
        let offset = DJ00 * DJC;
        assert_eq!((offset / DJC).to_bits(), DJ00.to_bits());
        let distant = uncached_offset(offset);
        for source_first in [false, true] {
            if source_first {
                assert_angles(source_nut06a_from_j2000_offset(offset), distant);
            }
            assert_angles(nut00a(DJ00).unwrap(), nut00a_centuries(0.0));
            assert_angles(source_nut06a_from_j2000_offset(offset), distant);
            assert_angles(nut06a(DJ00).unwrap(), uncached_offset(0.0));
        }
    }

    #[test]
    fn retained_nutation_offsets_keep_subday_and_ieee_key_resolution() {
        let tiny = 2.0_f64.powi(-33);
        assert_eq!(DJ00 + tiny, DJ00);
        assert_eq!(DJ00 - tiny, DJ00);
        let subnormal = f64::from_bits(1);
        assert_eq!((subnormal / DJC).to_bits(), 0.0_f64.to_bits());
        assert_ne!((-0.0_f64 / DJC).to_bits(), (0.0_f64 / DJC).to_bits());
        assert_eq!(tt_centuries(0.0).unwrap(), tt_centuries(subnormal).unwrap());
        for offset in [
            0.0,
            -0.0,
            subnormal,
            -subnormal,
            tiny,
            -tiny,
            DJC.next_down(),
            DJC,
            DJC.next_up(),
        ] {
            let expected = uncached_offset(offset);
            assert_angles(source_nut06a_from_j2000_offset(offset), expected);
            assert_angles(
                source_nut06a(DJ00 + offset),
                uncached_offset((DJ00 + offset) - DJ00),
            );
            assert_angles(source_nut06a_from_j2000_offset(offset), expected);
        }
        let retained = source_nut06a_from_j2000_offset(tiny);
        let collapsed = nut06a(DJ00 + tiny).unwrap();
        assert_ne!(
            [retained.0.to_bits(), retained.1.to_bits()],
            [collapsed.0.to_bits(), collapsed.1.to_bits()]
        );
    }

    #[test]
    fn native_nutation_cache_keeps_original_invalid_and_overflow_semantics() {
        nut00a(DJ00).unwrap();
        source_nut06a_from_j2000_offset(0.0);
        for bad in [
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::from_bits(0x7ff8_0000_0000_1234),
            f64::from_bits(0xfff8_0000_0000_5678),
        ] {
            for result in [nut00a(bad), nut06a(bad)] {
                match result {
                    Err(Error::NonFiniteArgument { value }) => {
                        assert_eq!(value.to_bits(), bad.to_bits());
                    }
                    other => panic!("checked epoch admission changed: {other:?}"),
                }
            }
            assert_angles(source_nut06a_from_j2000_offset(bad), uncached_offset(bad));
            assert_angles(source_nut06a(bad), uncached_offset(bad - DJ00));
            match c2t06a(DJ00, bad, f64::INFINITY, 0.0) {
                Err(Error::NonFiniteArgument { value }) => {
                    assert_eq!(value.to_bits(), bad.to_bits())
                }
                other => panic!("UT1 admission must precede polar motion: {other:?}"),
            }
        }
        for huge in [f64::MAX, -f64::MAX] {
            let t = tt_centuries(huge).unwrap();
            let raw = nut00a_centuries(t);
            assert!(raw.0.is_nan() && raw.1.is_nan());
            for _ in 0..2 {
                assert_angles(nut00a(huge).unwrap(), raw);
                assert_angles(nut06a(huge).unwrap(), uncached_offset(huge - DJ00));
                assert_angles(source_nut06a_from_j2000_offset(huge), uncached_offset(huge));
            }
            assert!(matches!(c2t06a(huge, DJ00, 0.0, 0.0),
                Err(Error::MathDomain { value }) if value.is_nan()));
        }
    }
}

#[cfg(test)]
mod argument_remainder_tests {
    use super::*;

    fn compare(angle: f64) {
        let expected = angle % D2PI;
        let actual = nutation_argument_remainder(angle);
        if expected.is_nan() {
            assert!(actual.is_nan());
        } else {
            assert_eq!(actual.to_bits(), expected.to_bits(), "angle {angle:e}");
        }
    }

    #[test]
    fn integer_turn_boundaries_and_signed_extremes_match_ieee_remainder() {
        for exponent in 0..=54 {
            for turn in [2.0_f64.powi(exponent), 2.0_f64.powi(exponent) - 1.0] {
                let boundary = turn * D2PI;
                for value in [boundary.next_down(), boundary, boundary.next_up()] {
                    compare(value);
                    compare(-value);
                }
            }
        }
        for turn in 0..4096 {
            let boundary = f64::from(turn) * D2PI;
            for value in [boundary.next_down(), boundary, boundary.next_up()] {
                compare(value);
                compare(-value);
            }
        }
        for value in [
            0.0,
            f64::from_bits(1),
            f64::MIN_POSITIVE,
            f64::MAX,
            f64::INFINITY,
            f64::NAN,
        ] {
            compare(value);
            compare(-value);
        }
    }

    #[test]
    fn binary64_samples_cover_fast_quotients_and_large_argument_dispatch() {
        let mut bits = 0x7e3a_950c_219b_8d47_u64;
        for index in 0..262_144 {
            bits ^= bits << 13;
            bits ^= bits >> 7;
            bits ^= bits << 17;
            // Concentrate samples on the optimized quotient domain, while
            // also checking arbitrary exponent/sign/payload combinations.
            let exponent = (1025 + index % 54) as u64;
            let angle = f64::from_bits((bits & ((1_u64 << 52) - 1)) | (exponent << 52));
            // Independent exact-integer oracle: throughout these exponents,
            // both words are integer multiples of 2^-50. The numerator fits
            // u128 and the integer remainder converts to binary64 exactly.
            let mantissa = (1_u64 << 52) | (bits & ((1_u64 << 52) - 1));
            let divisor = (1_u64 << 52) | (D2PI.to_bits() & ((1_u64 << 52) - 1));
            let numerator = u128::from(mantissa) << (exponent - 1025);
            let exact = (numerator % u128::from(divisor)) as f64 * 2.0_f64.powi(-50);
            assert_eq!(
                nutation_argument_remainder(angle).to_bits(),
                exact.to_bits()
            );
            compare(angle);
            compare(-angle);
            compare(f64::from_bits(bits));
        }
    }
}

#[cfg(test)]
mod paired_argument_tests {
    use super::*;

    fn assert_same(actual: f64, expected: f64) {
        if expected.is_nan() {
            // Arithmetic NaN payload/sign is not a portable Rust contract.
            assert!(actual.is_nan());
        } else {
            assert_eq!(actual.to_bits(), expected.to_bits());
        }
    }

    fn compare_table_arguments<const N: usize, const C: usize>(
        rows: &[([i8; N], [f64; C])],
        fa: &[f64; N],
    ) {
        for pair in rows.rchunks(2) {
            if let [earlier, later] = pair {
                let actual = paired_fused_arg_sum(&later.0, &earlier.0, fa);
                assert_same(actual[0], fused_arg_sum(&later.0, fa));
                assert_same(actual[1], fused_arg_sum(&earlier.0, fa));
            }
        }
    }

    // Frozen serial accumulation order, independent of paired loop traversal.
    fn serial_lunisolar(rows: &[([i8; 5], [f64; 6])], t: f64, fa: &[f64; 5]) -> (f64, f64) {
        let mut dp = 0.0;
        let mut de = 0.0;
        for (mult, coeff) in rows.iter().rev() {
            let (sarg, carg) = era_sincos(nutation_argument_remainder(fused_arg_sum(mult, fa)));
            dp += coeff[1].mul_add(t, coeff[0]).mul_add(sarg, coeff[2] * carg);
            de += coeff[4].mul_add(t, coeff[3]).mul_add(carg, coeff[5] * sarg);
        }
        (dp, de)
    }

    fn serial_planetary(fa: &[f64; 13]) -> (f64, f64) {
        let mut dp = 0.0;
        let mut de = 0.0;
        for (mult, coeff) in NUT_PL_ROWS.iter().rev() {
            let (sarg, carg) = era_sincos(nutation_argument_remainder(fused_arg_sum(mult, fa)));
            dp += coeff[0].mul_add(sarg, coeff[1] * carg);
            de += coeff[2].mul_add(sarg, coeff[3] * carg);
        }
        (dp, de)
    }

    #[test]
    fn paired_arguments_preserve_fused_cancellation_in_each_lane() {
        let e = f64::EPSILON;
        for (first, second, fa) in [
            (
                [3, 1, 0, 0, 0],
                [-3, -1, 0, 0, 0],
                [1.0 + e, -(3.0 + 4.0 * e), 0.0, 0.0, 0.0],
            ),
            (
                [0, 1, 3, 0, 0],
                [0, -1, -3, 0, 0],
                [0.0, -(3.0 + 4.0 * e), 1.0 + e, 0.0, 0.0],
            ),
        ] {
            // The exact product is 3 + 3*e. Rounding it before the addition
            // loses this residual; both the initial and later steps must fuse.
            let actual = paired_fused_arg_sum(&first, &second, &fa);
            assert_eq!(actual[0].to_bits(), (-e).to_bits());
            assert_eq!(actual[1].to_bits(), e.to_bits());
            assert_same(actual[0], fused_arg_sum(&first, &fa));
            assert_same(actual[1], fused_arg_sum(&second, &fa));
        }
    }

    #[test]
    fn paired_arguments_match_serial_tables_at_ieee_edges() {
        fn check<const N: usize, const C: usize>(rows: &[([i8; N], [f64; C])]) {
            for value in [
                0.0,
                -0.0,
                f64::from_bits(1),
                -f64::from_bits(1),
                f64::MIN_POSITIVE,
                -f64::MIN_POSITIVE,
                0.125,
                -2.5,
                1.0 + f64::EPSILON,
                f64::MAX,
                -f64::MAX,
            ] {
                compare_table_arguments(rows, &[value; N]);
                let alternating = std::array::from_fn(|k| if k % 2 == 0 { value } else { -value });
                compare_table_arguments(rows, &alternating);
            }
            for value in [
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::from_bits(0xfff8_0000_0000_0123),
            ] {
                for column in 0..N {
                    let mut fa = [0.25; N];
                    fa[column] = value;
                    compare_table_arguments(rows, &fa);
                }
            }
        }
        check(&NUT_LS_ROWS);
        check(&NUT_PL_ROWS);
        check(&NUT00B_ROWS);
    }

    #[test]
    fn paired_series_keep_serial_accumulation_and_odd_tails() {
        fn compare(t: f64, ls: [f64; 5], pl: [f64; 13]) {
            for rows in [&NUT_LS_ROWS[..], &NUT00B_ROWS[..]] {
                let actual = lunisolar_sum(rows, t, ls[0], ls[1], ls[2], ls[3], ls[4]);
                let expected = serial_lunisolar(rows, t, &ls);
                assert_same(actual.0, expected.0);
                assert_same(actual.1, expected.1);
            }
            let actual = planetary_sum(
                pl[0], pl[1], pl[2], pl[3], pl[4], pl[5], pl[6], pl[7], pl[8], pl[9], pl[10],
                pl[11], pl[12],
            );
            let expected = serial_planetary(&pl);
            assert_same(actual.0, expected.0);
            assert_same(actual.1, expected.1);
        }

        for century in [
            -2000.0, -152.0, -60.0, -50.0, -20.0, -10.0, -5.0, -4.0, -3.0, -2.0, -1.0, -0.5, -0.28,
            0.0, 0.01, 0.1, 0.24, 0.26, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 20.0, 50.0, 60.0,
            152.0, 2000.0,
        ] {
            let jd = DJ00 + century * DJC;
            for jd in [jd.next_down(), jd, jd.next_up()] {
                let t = (jd - DJ00) / DJC;
                let ls = [
                    fa_poly_rad(t, FA_L03),
                    mhb_sun_anomaly(t),
                    fa_poly_rad(t, FA_F03),
                    mhb_moon_elongation(t),
                    fa_poly_rad(t, FA_OM03),
                ];
                // Use the full model's distinct planetary arguments. The same
                // loop oracle also exercises every coefficient in the B table.
                let pl = [
                    8328.6914269554_f64.mul_add(t, 2.35555598) % D2PI,
                    8433.466158131_f64.mul_add(t, 1.627905234) % D2PI,
                    7771.3771468121_f64.mul_add(t, 5.198466741) % D2PI,
                    (-33.757045_f64).mul_add(t, 2.18243920) % D2PI,
                    2608.7903141574_f64.mul_add(t, 4.402608842) % D2PI,
                    1021.3285546211_f64.mul_add(t, 3.176146697) % D2PI,
                    628.3075849991_f64.mul_add(t, 1.753470314) % D2PI,
                    334.0612426700_f64.mul_add(t, 6.203480913) % D2PI,
                    52.9690962641_f64.mul_add(t, 0.599546497) % D2PI,
                    21.3299104960_f64.mul_add(t, 0.874016757) % D2PI,
                    7.4781598567_f64.mul_add(t, 5.481293872) % D2PI,
                    3.8127774000_f64.mul_add(t, 5.321159000) % D2PI,
                    0.00000538691_f64.mul_add(t, 0.024381750) * t,
                ];
                compare(t, ls, pl);
            }
        }
        for value in [
            0.0,
            -0.0,
            f64::from_bits(1),
            -f64::from_bits(1),
            f64::MAX,
            -f64::MAX,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NAN,
        ] {
            compare(value, [value; 5], [value; 13]);
        }
        // Empty, one-row and odd prefixes check that no padding row contributes.
        for length in [0, 1, 2, 3, 76, 77, 78, 677, 678] {
            let rows = &NUT_LS_ROWS[..length];
            let fa = [0.125, -0.5, 0.75, 1.25, -2.0];
            let actual = lunisolar_sum(rows, 0.25, fa[0], fa[1], fa[2], fa[3], fa[4]);
            let expected = serial_lunisolar(rows, 0.25, &fa);
            assert_same(actual.0, expected.0);
            assert_same(actual.1, expected.1);
        }
    }

    #[test]
    fn four_argument_chains_preserve_fused_lanes_and_row_tails() {
        let e = f64::EPSILON;
        for (mult, fa) in [
            (
                [
                    [3, 1, 0, 0, 0],
                    [-3, -1, 0, 0, 0],
                    [6, 2, 0, 0, 0],
                    [-6, -2, 0, 0, 0],
                ],
                [1.0 + e, -(3.0 + 4.0 * e), 0.0, 0.0, 0.0],
            ),
            (
                [
                    [0, 1, 3, 0, 0],
                    [0, -1, -3, 0, 0],
                    [0, 2, 6, 0, 0],
                    [0, -2, -6, 0, 0],
                ],
                [0.0, -(3.0 + 4.0 * e), 1.0 + e, 0.0, 0.0],
            ),
        ] {
            // Four distinct residuals detect lane swaps and an unfused product.
            let actual = interleaved_fused_arg_sums([&mult[0], &mult[1], &mult[2], &mult[3]], &fa);
            for (actual, expected) in actual.into_iter().zip([-e, e, -2.0 * e, 2.0 * e]) {
                assert_eq!(actual.to_bits(), expected.to_bits());
            }
        }
        fn compare_four<const N: usize, const C: usize>(rows: &[([i8; N], [f64; C])]) {
            for value in [
                0.0,
                -0.0,
                f64::from_bits(1),
                f64::MIN_POSITIVE,
                1.0 + f64::EPSILON,
                f64::MAX,
                f64::INFINITY,
                f64::NAN,
            ] {
                let fa = std::array::from_fn(|k| if k % 2 == 0 { value } else { -value });
                for group in rows.rchunks(4) {
                    if let [a, b, c, d] = group {
                        let mult = [&d.0, &c.0, &b.0, &a.0];
                        let actual = interleaved_fused_arg_sums(mult, &fa);
                        for (actual, mult) in actual.into_iter().zip(mult) {
                            assert_same(actual, fused_arg_sum(mult, &fa));
                        }
                    }
                }
            }
        }
        compare_four(&NUT_LS_ROWS);
        compare_four(&NUT_PL_ROWS);
        compare_four(&NUT00B_ROWS);

        // Exercise every tail after zero, one and two full groups. The frozen
        // serial oracle ties each argument to its own coefficient row.
        for length in 0..=11 {
            let rows = &NUT_LS_ROWS[..length];
            let fa = [0.125, -0.75, 1.25, -2.5, 4.0];
            let actual = lunisolar_sum(rows, 0.25, fa[0], fa[1], fa[2], fa[3], fa[4]);
            let expected = serial_lunisolar(rows, 0.25, &fa);
            assert_same(actual.0, expected.0);
            assert_same(actual.1, expected.1);
        }
    }
}
