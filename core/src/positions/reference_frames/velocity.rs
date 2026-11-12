//! Analytic TT tangent of the complete vector-Time IAU 2000A nutation series.
//!
//! Values and argument reductions remain with the original position kernels.
//! Phase rates differentiate the underlying periodic model, including at a
//! full-turn reduction, rather than the discontinuous numerical remainder.

use super::{
    ANOMALY_ROWS, Error, FRAME_ASEC2RAD, FRAME_EPOCH_J2000, FUND_ROWS, IAU2000_LUNISOLAR_LON,
    IAU2000_LUNISOLAR_OBL, IAU2000_NALS_T, IAU2000_NAPL_T, IAU2000_NUT_LON, IAU2000_NUT_OBL,
    check_finite, fundamental_arguments_radians, scalar_cos, scalar_sin,
};

const DAYS_PER_CENTURY: f64 = 36525.0;

/// Analytic IAU 2000A longitude and obliquity rates in tenth-microarcseconds
/// per TT day, using the same full series and argument realization as
/// [`super::iau2000a_tenth_microarcsec`].
///
/// `jd_tt` is a finite TT Julian day. All luni-solar secular amplitudes and
/// the quadratic planetary precession argument are differentiated. This is
/// the tangent of the underlying continuous model, not a derivative of its
/// binary64 rounding steps; no P03 adjustment or shortened series is used.
/// The value model's physical limitations remain. Nonfinite input or
/// nonfinite intermediate/output arithmetic returns an error.
pub(crate) fn iau2000a_rate_tenth_microarcsec_per_day(jd_tt: f64) -> Result<(f64, f64), Error> {
    let jd_tt = check_finite(jd_tt)?;
    let t = (jd_tt - FRAME_EPOCH_J2000) / DAYS_PER_CENTURY;
    let arguments = fundamental_arguments_radians(t)?;
    let rates = FUND_ROWS.map(|row| fundamental_rate(row, t));
    for value in arguments.into_iter().chain(rates) {
        check_finite(value)?;
    }
    let (mut dpsi, mut deps) = (0.0, 0.0);
    for ((multipliers, longitude), obliquity) in IAU2000_NALS_T
        .iter()
        .zip(IAU2000_LUNISOLAR_LON.iter())
        .zip(IAU2000_LUNISOLAR_OBL.iter())
    {
        let (phase, phase_rate) = phase_and_rate(multipliers, &arguments, &rates)?;
        let (psi, eps) = luni_solar_term_rate(
            t,
            scalar_sin(phase),
            scalar_cos(phase),
            phase_rate,
            *longitude,
            *obliquity,
        );
        dpsi += check_finite(psi)?;
        deps += check_finite(eps)?;
    }
    let (arguments, rates) = planetary_arguments_and_rates(t);
    for value in arguments.into_iter().chain(rates) {
        check_finite(value)?;
    }
    for ((multipliers, longitude), obliquity) in IAU2000_NAPL_T
        .iter()
        .zip(IAU2000_NUT_LON.iter())
        .zip(IAU2000_NUT_OBL.iter())
    {
        let (phase, phase_rate) = phase_and_rate(multipliers, &arguments, &rates)?;
        let sine = scalar_sin(phase);
        let cosine = scalar_cos(phase);
        dpsi += check_finite(periodic_term_rate(sine, cosine, phase_rate, *longitude))?;
        deps += check_finite(periodic_term_rate(sine, cosine, phase_rate, *obliquity))?;
    }
    Ok((check_finite(dpsi)?, check_finite(deps)?))
}

/// Preserve the original left-to-right dot products, each starting at +0.
fn phase_and_rate<const N: usize>(
    multipliers: &[i32; N],
    arguments: &[f64; N],
    rates: &[f64; N],
) -> Result<(f64, f64), Error> {
    let dot = |weights: &[f64; N]| {
        multipliers
            .iter()
            .zip(weights)
            .fold(0.0, |sum, (&multiplier, &value)| {
                sum + f64::from(multiplier) * value
            })
    };
    Ok((check_finite(dot(arguments))?, check_finite(dot(rates))?))
}

/// Derivative of the unreduced arcsecond polynomial, radians per TT day.
fn fundamental_rate(row: [f64; 5], t: f64) -> f64 {
    (((4.0 * row[4] * t + 3.0 * row[3]) * t + 2.0 * row[2]) * t + row[1]) * FRAME_ASEC2RAD
        / DAYS_PER_CENTURY
}

fn planetary_arguments_and_rates(t: f64) -> ([f64; 14], [f64; 14]) {
    let mut arguments = ANOMALY_ROWS.map(|row| t * row[1] + row[0]);
    let mut rates = ANOMALY_ROWS.map(|row| row[1] / DAYS_PER_CENTURY);
    arguments[13] *= t;
    rates[13] = (2.0 * ANOMALY_ROWS[13][1] * t + ANOMALY_ROWS[13][0]) / DAYS_PER_CENTURY;
    (arguments, rates)
}

fn luni_solar_term_rate(
    t: f64,
    sine: f64,
    cosine: f64,
    phase_rate: f64,
    longitude: [f64; 3],
    obliquity: [f64; 3],
) -> (f64, f64) {
    let psi = longitude[1] / DAYS_PER_CENTURY * sine
        + ((longitude[0] + longitude[1] * t) * cosine - longitude[2] * sine) * phase_rate;
    let eps = obliquity[1] / DAYS_PER_CENTURY * cosine
        + (-(obliquity[0] + obliquity[1] * t) * sine + obliquity[2] * cosine) * phase_rate;
    (psi, eps)
}

fn periodic_term_rate(sine: f64, cosine: f64, phase_rate: f64, coefficients: [f64; 2]) -> f64 {
    (coefficients[0] * cosine - coefficients[1] * sine) * phase_rate
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn polynomial_and_quadratic_arguments_keep_their_complete_tt_tangents() {
        // The derivative of 11 + 2t - 3t² + 4t³ - 5t⁴ at t=1/2 is -1/2.
        assert_eq!(
            fundamental_rate([11.0, 2.0, -3.0, 4.0, -5.0], 0.5),
            -0.5 * FRAME_ASEC2RAD / DAYS_PER_CENTURY
        );
        for t in [-2.0, 0.0, 0.5, 2.0] {
            let (arguments, rates) = planetary_arguments_and_rates(t);
            for (index, row) in ANOMALY_ROWS.iter().enumerate().take(13) {
                assert_eq!(arguments[index], t * row[1] + row[0]);
                assert_eq!(rates[index], row[1] / DAYS_PER_CENTURY);
            }
            let [constant, linear] = ANOMALY_ROWS[13];
            assert_eq!(arguments[13], (linear * t + constant) * t);
            assert_eq!(rates[13], (constant + 2.0 * linear * t) / DAYS_PER_CENTURY);
        }
    }

    #[test]
    fn periodic_tangents_include_phase_and_secular_amplitude_terms() {
        let longitude = [3.0, 5.0, 7.0];
        let obliquity = [11.0, 13.0, 17.0];
        assert_eq!(
            luni_solar_term_rate(2.0, 0.0, 1.0, 0.25, longitude, obliquity),
            (13.0 * 0.25, 13.0 / DAYS_PER_CENTURY + 17.0 * 0.25)
        );
        assert_eq!(
            luni_solar_term_rate(2.0, 1.0, 0.0, 0.25, longitude, obliquity),
            (5.0 / DAYS_PER_CENTURY - 7.0 * 0.25, -37.0 * 0.25)
        );
        assert_eq!(periodic_term_rate(0.0, 1.0, 0.25, [3.0, 7.0]), 0.75);
        assert_eq!(periodic_term_rate(1.0, 0.0, 0.25, [3.0, 7.0]), -1.75);
        assert_eq!(
            luni_solar_term_rate(2.0, 1.0, 0.0, 0.0, longitude, obliquity),
            (5.0 / DAYS_PER_CENTURY, 0.0)
        );
    }

    #[test]
    fn full_series_dimensions_and_invalid_epochs_are_explicit() {
        assert_eq!(IAU2000_NALS_T.len(), 678);
        assert_eq!(IAU2000_LUNISOLAR_LON.len(), IAU2000_NALS_T.len());
        assert_eq!(IAU2000_LUNISOLAR_OBL.len(), IAU2000_NALS_T.len());
        assert_eq!(IAU2000_NAPL_T.len(), 687);
        assert_eq!(IAU2000_NUT_LON.len(), IAU2000_NAPL_T.len());
        assert_eq!(IAU2000_NUT_OBL.len(), IAU2000_NAPL_T.len());
        for jd in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, f64::MAX] {
            assert!(iau2000a_rate_tenth_microarcsec_per_day(jd).is_err());
        }
    }
}
