//! TT derivatives of the actual output-frame transformation.
//!
//! The physical velocity is `R v + dR/dTT r`. A frame's own smooth precession
//! and nutation model determines its derivative; a short body polynomial must
//! not truncate that independent rotation. No position or body state is sampled
//! by this module. All matrices are dimensionless and their derivatives are per
//! TT day. Scalar position callers retain their existing transformations.

use super::*;

type Matrix3 = [[f64; 3]; 3];
type WordMatrix3 = [[super::derivative::AngularWord; 3]; 3];

/// Native JPL output matrix and its TT derivative.
///
/// `jd_tt` is a finite TT Julian day. Flags select the same ICRS/bias,
/// equatorial/ecliptic, J2000 and nutation branches as the native scalar place.
/// IAU 2000A nutation and Vondrák precession are evaluated on their own smooth
/// model, without body/observer/deflector record restrictions. The native
/// NONUT ecliptic-of-date projection keeps its actual precession/obliquity
/// construction. J2000 is static, with an exactly zero matrix derivative.
/// Model, finite-input and differentiation errors propagate.
pub(crate) fn native_output_frame_matrix_rate(
    jd_tt: f64,
    flags: i32,
) -> Result<(Matrix3, Matrix3)> {
    let matrix_at = |epoch| {
        let (dpsi, deps) = reference_frames::iau2000a_tenth_microarcsec(epoch)?;
        let scale = reference_frames::FRAME_TENTH_USEC_2_RAD;
        matrix_with_nutation(epoch, flags, (dpsi * scale, deps * scale), true)
    };
    matrix_rate(jd_tt, flags & FLG_J2000 != 0, &matrix_at)
}

/// Analytic TT tangent of the planet-centered output frame.
///
/// `jd_tt` is TT JD; `nutation` is the actual full vector-Time IAU 2000A
/// `(dpsi,deps)` in radians and `nutation_rate` is its tangent in radians per
/// TT day. The central matrix keeps the native scalar construction. PCTR
/// selects the biased J2000 frame before ICRS, and its NONUT ecliptic frame
/// rotates Julian-epoch precession by the direct-TT Vondrák obliquity.
/// Other branches use the existing analytic output-frame derivatives.
/// Nonfinite operands and frame-model failures propagate. No body state or
/// nutation value is sampled, and scalar position callers remain unchanged.
pub(crate) fn pctr_output_frame_matrix_rate(
    jd_tt: f64,
    flags: i32,
    nutation: (f64, f64),
    nutation_rate: (f64, f64),
) -> Result<(Matrix3, Matrix3)> {
    check_finite_vec3([jd_tt, nutation.0, nutation.1])?;
    check_finite_vec3([nutation_rate.0, nutation_rate.1, 0.0])?;
    let flags = if flags & FLG_J2000 != 0 {
        flags & !FLG_ICRS
    } else {
        flags
    };
    if flags & FLG_J2000 == 0 && flags & FLG_EQUATORIAL != 0 {
        let central = matrix_with_nutation(jd_tt, flags, nutation, true)?;
        let (precession, rate) = pctr_precession_matrix_rate(jd_tt, flags & FLG_ICRS == 0)?;
        let rate = if flags & (FLG_NONUT | FLG_SIDEREAL) != 0 {
            rate
        } else {
            equatorial_nutation_rate(jd_tt, precession, rate, nutation, nutation_rate)?
        };
        return Ok((central, rate));
    }
    if flags & FLG_NONUT != 0 && flags & (FLG_EQUATORIAL | FLG_J2000) == 0 {
        let central = matrix_with_nutation(jd_tt, flags, nutation, true)?;
        let (precession, precession_rate) =
            pctr_precession_matrix_rate(jd_tt, flags & FLG_ICRS == 0)?;
        let (eps, eps_rate) = super::corrections::mean_obliquity_and_rate(jd_tt)?;
        return Ok((
            central,
            non_nutated_ecliptic_rate(precession, precession_rate, eps, eps_rate)?,
        ));
    }
    analytic_output_frame_matrix_rate_with_nutation(jd_tt, flags, nutation, nutation_rate)
}

/// Differentiate the scalar precession's actual angle-built fixed bias.
fn pctr_precession_matrix_rate(jd_tt: f64, bias: bool) -> Result<(Matrix3, Matrix3)> {
    let (precession, rate) = super::corrections::precession_matrix_and_rate(jd_tt, false)?;
    if bias {
        Ok((
            vondrak_precession_matrix(jd_tt, true)?,
            crate::erfa::apply_long_term_frame_bias(rate)?,
        ))
    } else {
        Ok((precession, rate))
    }
}

/// Product rule for Rx(epsilon) P, without regrouping the central position.
fn non_nutated_ecliptic_rate(
    precession: Matrix3,
    precession_rate: Matrix3,
    eps: f64,
    eps_rate: f64,
) -> Result<Matrix3> {
    let (rotation, rotation_rate) = x_rotation_rate(eps, eps_rate);
    let rate = add_matrices(
        mat3_mul(rotation_rate, precession)?,
        mat3_mul(rotation, precession_rate)?,
    );
    for row in rate {
        check_finite_vec3(row)?;
    }
    Ok(rate)
}

/// Standard output matrix with the actual selected nutation value and rate.
///
/// Inputs are TT Julian day, flags, `(dpsi,deps)` radians and their derivatives
/// in radians per TT day. The caller obtains the rate from the same owned
/// nutation polynomial, or differentiates its actual replacement callback on
/// that callback's domain. Advancing this first-order tangent differentiates
/// the central frame by the chain rule without crossing a nutation join.
/// Precession remains its own smooth TT model. The returned central matrix is
/// constructed with the exact supplied values; no changed scalar position is
/// published. J2000 has an exactly zero derivative. Nonfinite inputs and actual
/// frame/differentiation errors propagate.
pub(crate) fn output_frame_matrix_rate_with_nutation(
    jd_tt: f64,
    flags: i32,
    nutation: (f64, f64),
    nutation_rate: (f64, f64),
) -> Result<(Matrix3, Matrix3)> {
    check_finite_vec3([jd_tt, nutation.0, nutation.1])?;
    check_finite_vec3([nutation_rate.0, nutation_rate.1, 0.0])?;
    // The caller has already acquired its actual nutation feed and derivative.
    // Only this pure transformation is reused; no reader/callback is bypassed.
    let key = [
        jd_tt.to_bits(),
        flags as u64,
        nutation.0.to_bits(),
        nutation.1.to_bits(),
        nutation_rate.0.to_bits(),
        nutation_rate.1.to_bits(),
    ];
    if let Some(result) = crate::pure_cache::frame_rate(key) {
        return Ok(result);
    }
    let result = output_frame_matrix_rate_uncached(jd_tt, flags, nutation, nutation_rate)?;
    crate::pure_cache::retain_frame_rate(key, result);
    Ok(result)
}

/// Evaluate the same matrix stencil without retained matrix/rate reuse.
fn output_frame_matrix_rate_uncached(
    jd_tt: f64,
    flags: i32,
    nutation: (f64, f64),
    nutation_rate: (f64, f64),
) -> Result<(Matrix3, Matrix3)> {
    let matrix_at = |epoch| {
        let offset = epoch - jd_tt;
        matrix_with_nutation(
            epoch,
            flags,
            (
                nutation.0 + nutation_rate.0 * offset,
                nutation.1 + nutation_rate.1 * offset,
            ),
            false,
        )
    };
    matrix_rate(jd_tt, flags & FLG_J2000 != 0, &matrix_at)
}

/// Analytic derivative for providers supplying complete apparent kinematics.
/// The central matrix and complete pole/nutation models are unchanged. Older
/// provider paths keep their separately qualified sampled-frame reduction.
pub(crate) fn analytic_output_frame_matrix_rate_with_nutation(
    jd_tt: f64,
    flags: i32,
    nutation: (f64, f64),
    nutation_rate: (f64, f64),
) -> Result<(Matrix3, Matrix3)> {
    check_finite_vec3([jd_tt, nutation.0, nutation.1])?;
    check_finite_vec3([nutation_rate.0, nutation_rate.1, 0.0])?;
    if flags & FLG_J2000 != 0 {
        return Ok((
            matrix_with_nutation(jd_tt, flags, nutation, false)?,
            [[0.0; 3]; 3],
        ));
    }
    if flags & FLG_EQUATORIAL == 0 {
        ecliptic_matrix_rate(jd_tt, flags, nutation, nutation_rate)
    } else {
        equatorial_matrix_rate(jd_tt, flags, nutation, nutation_rate)
    }
}

/// Analytic precession and the full nutation rotation in the equatorial basis.
fn equatorial_matrix_rate(
    jd_tt: f64,
    flags: i32,
    nutation: (f64, f64),
    nutation_rate: (f64, f64),
) -> Result<(Matrix3, Matrix3)> {
    let central = matrix_with_nutation(jd_tt, flags, nutation, false)?;
    let (p, pd) = super::corrections::precession_matrix_and_rate(jd_tt, flags & FLG_ICRS == 0)?;
    if flags & (FLG_SIDEREAL | FLG_NONUT) != 0 {
        return Ok((central, pd));
    }
    Ok((
        central,
        equatorial_nutation_rate(jd_tt, p, pd, nutation, nutation_rate)?,
    ))
}

/// Product rule for the complete nutation rotation times mean precession.
fn equatorial_nutation_rate(
    jd_tt: f64,
    p: Matrix3,
    pd: Matrix3,
    nutation: (f64, f64),
    nutation_rate: (f64, f64),
) -> Result<Matrix3> {
    let (eps, eps_rate) = super::corrections::mean_obliquity_and_rate(jd_tt)?;
    let (a, ad) = x_rotation_rate(eps, eps_rate);
    let (b, bd) = z_rotation_rate(-nutation.0, -nutation_rate.0);
    let (c, cd) = x_rotation_rate(-eps - nutation.1, -eps_rate - nutation_rate.1);
    let ba = mat3_mul(b, a)?;
    let bad = add_matrices(mat3_mul(bd, a)?, mat3_mul(b, ad)?);
    let nd = add_matrices(mat3_mul(cd, ba)?, mat3_mul(c, bad)?);
    let n = nutation_matrix(eps, nutation.0, nutation.1)?;
    Ok(add_matrices(mat3_mul(nd, p)?, mat3_mul(n, pd)?))
}

fn add_matrices(a: Matrix3, b: Matrix3) -> Matrix3 {
    std::array::from_fn(|i| std::array::from_fn(|j| a[i][j] + b[i][j]))
}

/// Passive elementary rotations used by the scalar ERFA nutation sequence.
fn x_rotation_rate(angle: f64, rate: f64) -> (Matrix3, Matrix3) {
    let (s, c) = angle.sin_cos();
    (
        [[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]],
        [
            [0.0; 3],
            [0.0, -s * rate, c * rate],
            [0.0, -c * rate, -s * rate],
        ],
    )
}
fn z_rotation_rate(angle: f64, rate: f64) -> (Matrix3, Matrix3) {
    let (s, c) = angle.sin_cos();
    (
        [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]],
        [
            [-s * rate, c * rate, 0.0],
            [-c * rate, -s * rate, 0.0],
            [0.0; 3],
        ],
    )
}

/// Differentiate the direct ecliptic basis and its longitude-nutation rotation.
/// Obliquity nutation cancels in this physical projection; no sampled frame or
/// shortened nutation model is involved. Bias is a fixed linear transformation.
fn ecliptic_matrix_rate(
    jd_tt: f64,
    flags: i32,
    nutation: (f64, f64),
    nutation_rate: (f64, f64),
) -> Result<(Matrix3, Matrix3)> {
    let central = matrix_with_nutation(jd_tt, flags, nutation, false)?;
    let (mut basis, mut rate) = super::corrections::mean_ecliptic_basis_and_rate(jd_tt)?;
    if flags & FLG_ICRS == 0 {
        basis = mat3_mul(basis, FRAME_BIAS_LTPB)?;
        rate = mat3_mul(rate, FRAME_BIAS_LTPB)?;
    }
    if flags & FLG_NONUT == 0 {
        let (s, c) = nutation.0.sin_cos();
        let rotate = |matrix: Matrix3| {
            [
                std::array::from_fn(|i| c * matrix[0][i] - s * matrix[1][i]),
                std::array::from_fn(|i| s * matrix[0][i] + c * matrix[1][i]),
                matrix[2],
            ]
        };
        let rotated = rotate(basis);
        rate = rotate(rate);
        for column in 0..3 {
            rate[0][column] -= nutation_rate.0 * rotated[1][column];
            rate[1][column] += nutation_rate.0 * rotated[0][column];
        }
    }
    for row in rate {
        check_finite_vec3(row)?;
    }
    Ok((central, rate))
}

/// Evaluate the existing transformation on Cartesian unit columns.
fn matrix_with_nutation(
    jd_tt: f64,
    flags: i32,
    nutation: (f64, f64),
    native_non_nutated: bool,
) -> Result<Matrix3> {
    let frames = frame_matrices_for_date(jd_tt, nutation.0, nutation.1)?;
    let special =
        native_non_nutated && flags & FLG_NONUT != 0 && flags & (FLG_EQUATORIAL | FLG_J2000) == 0;
    let precession = if special {
        Some(vondrak_precession_matrix(jd_tt, flags & FLG_ICRS == 0)?)
    } else {
        None
    };
    let mut matrix = [[0.0; 3]; 3];
    for column in 0..3 {
        let mut unit = [0.0; 3];
        unit[column] = 1.0;
        let value = if let Some(precession) = precession {
            rotate_equatorial_to_ecliptic(
                mat3_vec3(precession, unit)?,
                vondrak_mean_obliquity_rad(jd_tt)?,
            )?
        } else {
            frame_transform(unit, jd_tt, flags, true, frames)?
        };
        for row in 0..3 {
            matrix[row][column] = value[row];
        }
    }
    Ok(matrix)
}

/// Differentiate only a smooth matrix callback, reusing native fourth order.
/// Dyadic TT steps are independent of all geometric provider domains. Each
/// column uses the same adaptive truncation/roundoff rule as a Cartesian place.
fn matrix_rate(
    jd_tt: f64,
    fixed: bool,
    matrix_at: &dyn Fn(f64) -> Result<Matrix3>,
) -> Result<(Matrix3, Matrix3)> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let central = matrix_at(jd_tt)?;
    if fixed {
        return Ok((central, [[0.0; 3]; 3]));
    }
    const HALF_STEP: f64 = 1.0 / 64.0;
    let cache = std::cell::RefCell::new(std::collections::HashMap::new());
    let mut rate = [[0.0; 3]; 3];
    for column in 0..3 {
        let sample = |offset: f64| {
            if offset == 0.0 {
                return Ok(std::array::from_fn(|row| central[row][column]));
            }
            let epoch = jd_tt + offset;
            if epoch == jd_tt || epoch - jd_tt != offset {
                return Err(Error::MathDomain { value: epoch });
            }
            let cached = cache.borrow().get(&offset.to_bits()).copied();
            let matrix = match cached {
                Some(value) => value,
                None => {
                    let value = matrix_at(epoch)?;
                    cache.borrow_mut().insert(offset.to_bits(), value);
                    value
                }
            };
            Ok(std::array::from_fn(|row| matrix[row][column]))
        };
        let derivative = super::derivative::cartesian_rates(
            &sample,
            HALF_STEP,
            (f64::NEG_INFINITY, f64::INFINITY),
        )?;
        for row in 0..3 {
            rate[row][column] = derivative[row];
        }
    }
    Ok((central, rate))
}

#[cfg(test)]
mod pctr_frame_tests {
    use super::*;

    #[test]
    fn central_native_frames_match_and_j2000_bias_precedes_icrs() {
        let nutation = (1.0 / 65536.0, -1.0 / 131072.0);
        let rates = (1.0 / 1048576.0, -1.0 / 2097152.0);
        for tt in [FRAME_J2000, FRAME_J2000 + 12345.125] {
            for bits in 0..16 {
                let choices = [FLG_EQUATORIAL, FLG_J2000, FLG_ICRS, FLG_NONUT];
                let flags = choices.iter().enumerate().fold(0, |flags, (i, bit)| {
                    flags | if bits & (1 << i) != 0 { *bit } else { 0 }
                });
                let normalized = if flags & FLG_J2000 != 0 {
                    flags & !FLG_ICRS
                } else {
                    flags
                };
                let (matrix, rate) =
                    pctr_output_frame_matrix_rate(tt, flags, nutation, rates).unwrap();
                let expected = matrix_with_nutation(tt, normalized, nutation, true).unwrap();
                assert_eq!(
                    matrix.map(|row| row.map(f64::to_bits)),
                    expected.map(|row| row.map(f64::to_bits))
                );
                if flags & FLG_J2000 != 0 {
                    assert_eq!(rate.map(|row| row.map(f64::to_bits)), [[0; 3]; 3]);
                    let (without_icrs, _) =
                        pctr_output_frame_matrix_rate(tt, flags & !FLG_ICRS, nutation, rates)
                            .unwrap();
                    assert_eq!(matrix, without_icrs);
                }
            }
        }
    }

    #[test]
    fn non_nutated_rotation_product_keeps_both_closed_form_tangents() {
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        // At zero angles, precession rotates about z at 1/4 rad/day,
        // while obliquity rotates about x at 1/2 rad/day.
        let precession_rate = [[0.0, 0.25, 0.0], [-0.25, 0.0, 0.0], [0.0; 3]];
        let rate = non_nutated_ecliptic_rate(identity, precession_rate, 0.0, 0.5).unwrap();
        assert_eq!(
            rate,
            [[0.0, 0.25, 0.0], [-0.25, 0.0, 0.5], [0.0, -0.5, 0.0]]
        );
        let position = [3.0, 5.0, 7.0];
        let velocity = [11.0, 13.0, 17.0];
        let rotation = mat3_vec3(rate, position).unwrap();
        let result: [f64; 3] = std::array::from_fn(|i| velocity[i] + rotation[i]);
        assert_eq!(result, [12.25, 15.75, 14.5]);
        // A pure frame rotation cannot contribute a radial tangent.
        assert_eq!(compensated_dot3(position, rotation), 0.0);
    }

    #[test]
    fn sidereal_equatorial_mean_frame_is_independent_of_nutation_motion() {
        for tt in [FRAME_J2000, FRAME_J2000 + 12345.125] {
            for bias in [0, FLG_ICRS] {
                let flags = FLG_EQUATORIAL | FLG_SIDEREAL | bias;
                let first = pctr_output_frame_matrix_rate(
                    tt,
                    flags,
                    (1.0 / 65536.0, -1.0 / 131072.0),
                    (1.0 / 256.0, -1.0 / 512.0),
                )
                .unwrap();
                let second =
                    pctr_output_frame_matrix_rate(tt, flags, (-0.125, 0.25), (0.5, -0.75)).unwrap();
                let mean = pctr_output_frame_matrix_rate(
                    tt,
                    FLG_EQUATORIAL | FLG_NONUT | bias,
                    (0.0, 0.0),
                    (0.0, 0.0),
                )
                .unwrap();
                // SIDEREAL selects the local mean equator. Its matrix and
                // tangent cannot depend on either nutation angle or motion.
                assert_eq!(first, second);
                assert_eq!(first, mean);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn words(pair: (Matrix3, Matrix3)) -> [[u64; 9]; 2] {
        [pair.0, pair.1]
            .map(|matrix| std::array::from_fn(|index| matrix[index / 3][index % 3].to_bits()))
    }

    #[test]
    fn ecliptic_rates_improve_each_high_precision_series_witness() {
        // Independent 100-decimal mpmath differentiation of the full published
        // pole series and normalized cross products, using the exact binary64
        // input coefficients/constants. These are model derivatives, not
        // ephemeris coordinates or an accuracy claim for ancient epochs.
        let cases = [
            (
                1000000.25,
                [
                    [
                        5.382123398617646e-7,
                        -3.4184501426935036e-7,
                        -1.571553919827694e-7,
                    ],
                    [
                        3.7620557931353677e-7,
                        4.939344853342043e-7,
                        2.1385119545109747e-7,
                    ],
                    [
                        -1.5893168177554778e-9,
                        5.8523163558002765e-9,
                        2.603085148793443e-9,
                    ],
                ],
            ),
            (
                2415020.5,
                [
                    [
                        1.6261203029997332e-8,
                        -6.116429113734654e-7,
                        -2.660608449071611e-7,
                    ],
                    [
                        6.670046125955577e-7,
                        1.738525621332096e-8,
                        7.841987690496807e-10,
                    ],
                    [
                        5.057808366902942e-10,
                        5.7123998719095165e-9,
                        2.4781796530105537e-9,
                    ],
                ],
            ),
            (
                2451545.0,
                [
                    [
                        6.95606367042829e-24,
                        -6.121939657207821e-7,
                        -2.6602591422739897e-7,
                    ],
                    [
                        6.674960135078973e-7,
                        2.471566695510153e-9,
                        -5.70072800741611e-9,
                    ],
                    [
                        5.57365702804856e-10,
                        5.700728007416114e-9,
                        2.4715666955101465e-9,
                    ],
                ],
            ),
            (
                2488069.5,
                [
                    [
                        -1.6282651045143106e-8,
                        -6.123811870622751e-7,
                        -2.6583283392322356e-7,
                    ],
                    [
                        6.675907639204606e-7,
                        -1.2464874039737777e-8,
                        -1.2187113674155359e-8,
                    ],
                    [
                        6.087781563807167e-10,
                        5.688671192176997e-9,
                        2.4647934668716766e-9,
                    ],
                ],
            ),
            (
                5000000.75,
                [
                    [
                        -6.760911230268774e-7,
                        9.35546743746177e-8,
                        4.8374617970758846e-8,
                    ],
                    [
                        -1.0498190740731185e-7,
                        -6.231247526483194e-7,
                        -2.624318643404314e-7,
                    ],
                    [
                        3.57142292625423e-9,
                        4.052893851677661e-9,
                        1.6694284658982304e-9,
                    ],
                ],
            ),
        ];
        for (tt, exact) in cases {
            let flags = FLG_NONUT | FLG_ICRS;
            let (_, analytic) =
                analytic_output_frame_matrix_rate_with_nutation(tt, flags, (0.0, 0.0), (0.0, 0.0))
                    .unwrap();
            let (_, numerical) = matrix_rate(tt, false, &|epoch| {
                matrix_with_nutation(epoch, flags, (0.0, 0.0), false)
            })
            .unwrap();
            for row in 0..3 {
                for column in 0..3 {
                    let error = (analytic[row][column] - exact[row][column]).abs();
                    let previous_error = (numerical[row][column] - exact[row][column]).abs();
                    assert!(
                        error <= previous_error,
                        "epoch={tt} row={row} column={column}: {error:e} > {previous_error:e}"
                    );
                    assert!(error < 2e-21);
                }
            }
        }
    }

    #[test]
    fn analytic_ecliptic_rates_preserve_the_central_matrix_and_rotation_metric() {
        for tt in [1000000.25, 2415020.5, 2451545.0, 2488069.5, 5000000.75] {
            for flags in [0, FLG_ICRS, FLG_NONUT, FLG_NONUT | FLG_ICRS] {
                let nutation = (0.00001, -0.00002);
                let nutation_rate = (0.000003, -0.000004);
                let (matrix, analytic) = analytic_output_frame_matrix_rate_with_nutation(
                    tt,
                    flags,
                    nutation,
                    nutation_rate,
                )
                .unwrap();
                let at = |epoch: f64| {
                    matrix_with_nutation(
                        epoch,
                        flags,
                        (
                            nutation.0 + nutation_rate.0 * (epoch - tt),
                            nutation.1 + nutation_rate.1 * (epoch - tt),
                        ),
                        false,
                    )
                };
                let (same, numerical) = matrix_rate(tt, false, &at).unwrap();
                assert_eq!(matrix, same);
                for row in 0..3 {
                    for column in 0..3 {
                        assert!(
                            (analytic[row][column] - numerical[row][column]).abs() < 2e-13,
                            "epoch={tt} flags={flags} row={row} column={column}"
                        );
                        // A rotating basis preserves its constant metric,
                        // including the fixed first-order frame-bias map.
                        let metric_rate: f64 = (0..3)
                            .map(|i| {
                                analytic[i][row] * matrix[i][column]
                                    + matrix[i][row] * analytic[i][column]
                            })
                            .sum();
                        assert!(metric_rate.abs() < 2e-20);
                    }
                }
            }
        }
    }

    #[test]
    fn analytic_equatorial_rates_preserve_the_central_matrix_and_rotation_metric() {
        for tt in [1000000.25, 2415020.5, 2451545.0, 2488069.5, 5000000.75] {
            for option in [
                0,
                FLG_ICRS,
                FLG_NONUT,
                FLG_SIDEREAL,
                FLG_ICRS | FLG_NONUT,
                FLG_ICRS | FLG_SIDEREAL,
            ] {
                let flags = FLG_EQUATORIAL | option;
                let angles = (0.00001, -0.00002);
                let rates = (0.000003, -0.000004);
                let (matrix, analytic) =
                    analytic_output_frame_matrix_rate_with_nutation(tt, flags, angles, rates)
                        .unwrap();
                let (same, numerical) =
                    output_frame_matrix_rate_with_nutation(tt, flags, angles, rates).unwrap();
                assert_eq!(matrix, same);
                for row in 0..3 {
                    for column in 0..3 {
                        assert!(
                            (analytic[row][column] - numerical[row][column]).abs() < 2e-13,
                            "epoch={tt} flags={flags} row={row} column={column}"
                        );
                        let metric_rate: f64 = (0..3)
                            .map(|i| {
                                analytic[i][row] * matrix[i][column]
                                    + matrix[i][row] * analytic[i][column]
                            })
                            .sum();
                        assert!(metric_rate.abs() < 2e-19);
                    }
                }
            }
        }
    }

    #[test]
    fn retained_frame_rates_are_exact_under_replacement_and_eight_threads() {
        let inputs: Vec<_> = (0..40)
            .map(|index| {
                let flags = [
                    0,
                    FLG_EQUATORIAL,
                    FLG_J2000,
                    FLG_NONUT | FLG_EQUATORIAL,
                    FLG_ICRS,
                    FLG_SIDEREAL | FLG_EQUATORIAL,
                ][index % 6];
                let tt = 2415020.5 + index as f64 * 2191.25;
                let angles = (index as f64 * 1e-6, -(index as f64) * 2e-6);
                let rate = (index as f64 * 3e-8, -(index as f64) * 4e-8);
                let expected =
                    words(output_frame_matrix_rate_uncached(tt, flags, angles, rate).unwrap());
                (tt, flags, angles, rate, expected)
            })
            .collect();
        std::thread::scope(|scope| {
            for thread in 0..8 {
                let inputs = &inputs;
                scope.spawn(move || {
                    for index in 0..inputs.len() * 2 {
                        let (tt, flags, angles, rate, expected) =
                            inputs[(index + thread * 3) % inputs.len()];
                        for _ in 0..2 {
                            assert_eq!(
                                words(
                                    output_frame_matrix_rate_with_nutation(tt, flags, angles, rate)
                                        .unwrap()
                                ),
                                expected,
                            );
                        }
                    }
                });
            }
        });
    }

    #[test]
    fn invalid_frame_operands_still_fail_before_retained_reuse() {
        for index in 0..5 {
            for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                let mut operands = [2451545.0, 1e-5, -2e-5, 3e-8, -4e-8];
                operands[index] = invalid;
                assert!(matches!(
                    output_frame_matrix_rate_with_nutation(
                        operands[0],
                        0,
                        (operands[1], operands[2]),
                        (operands[3], operands[4])
                    ),
                    Err(Error::NonFiniteArgument { .. })
                ));
            }
        }
    }

    #[test]
    fn j2000_output_is_static_in_every_retained_basis() {
        for epoch in [2415020.5, 2451545.0, 2488069.5] {
            for flags in [
                FLG_J2000,
                FLG_J2000 | FLG_ICRS,
                FLG_J2000 | FLG_EQUATORIAL,
                FLG_J2000 | FLG_EQUATORIAL | FLG_ICRS,
            ] {
                let (matrix, rate) = native_output_frame_matrix_rate(epoch, flags).unwrap();
                assert_eq!(rate, [[0.0; 3]; 3]);
                if flags & (FLG_EQUATORIAL | FLG_ICRS) == FLG_EQUATORIAL | FLG_ICRS {
                    assert_eq!(matrix, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
                }
            }
        }
    }

    #[test]
    fn a_rotating_frame_obeys_the_closed_form_velocity_and_range_invariant() {
        let epoch = 2451545.0;
        let omega: f64 = 0.03;
        let rotating = |tt: f64| {
            let (sine, cosine) = (0.7 + omega * (tt - epoch)).sin_cos();
            Ok([[cosine, -sine, 0.0], [sine, cosine, 0.0], [0.0, 0.0, 1.0]])
        };
        let (matrix, rate) = matrix_rate(epoch, false, &rotating).unwrap();
        let [x, y, z] = [3.0, -4.0, 12.0];
        let position = mat3_vec3(matrix, [x, y, z]).unwrap();
        let velocity = [0.003, 0.002, -0.001];
        let rotated_velocity = mat3_vec3(matrix, velocity).unwrap();
        let frame_velocity = mat3_vec3(rate, [x, y, z]).unwrap();
        let expected = [-omega * position[1], omega * position[0], 0.0];
        for (actual, exact) in frame_velocity.into_iter().zip(expected) {
            assert!((actual - exact).abs() < 1e-10);
        }
        let output_velocity: [f64; 3] =
            std::array::from_fn(|i| rotated_velocity[i] + frame_velocity[i]);
        let original_dot = x * velocity[0] + y * velocity[1] + z * velocity[2];
        let rotated_dot: f64 = position
            .into_iter()
            .zip(output_velocity)
            .map(|(p, v)| p * v)
            .sum();
        assert!((original_dot - rotated_dot).abs() < 1e-10);
    }

    #[test]
    fn the_selected_nutation_tangent_changes_rates_without_changing_position() {
        let tt = 2460310.5;
        let angles = (0.00001, -0.00002);
        let (plain, zero) =
            output_frame_matrix_rate_with_nutation(tt, 0, angles, (0.0, 0.0)).unwrap();
        let (same, moving) =
            output_frame_matrix_rate_with_nutation(tt, 0, angles, (0.00001, 0.00002)).unwrap();
        assert_eq!(same, plain);
        // In true ecliptic coordinates obliquity nutation cancels, and
        // longitude nutation is precisely a rotation about the ecliptic pole.
        for column in 0..3 {
            let exact = [-0.00001 * same[1][column], 0.00001 * same[0][column], 0.0];
            for row in 0..3 {
                assert!((moving[row][column] - zero[row][column] - exact[row]).abs() < 1e-11);
            }
        }
    }
}

/// Analytic frame construction retaining physical value/rate residuals.
/// The original scalar position matrices are untouched. J2000 follows the same
/// bias/obliquity branches, evaluated before coefficient rounding. Other frame
/// branches retain their currently available binary64 matrix/rate inputs.
pub(crate) fn analytic_output_frame_words_with_nutation(
    jd_tt: f64,
    flags: i32,
    nutation: (f64, f64),
    nutation_rate: (f64, f64),
) -> Result<(WordMatrix3, WordMatrix3)> {
    check_finite_vec3([jd_tt, nutation.0, nutation.1])?;
    check_finite_vec3([nutation_rate.0, nutation_rate.1, 0.0])?;
    // The actual feed and its rate were acquired upstream. These six exact
    // operands completely specify this pure frame; no provider can be skipped.
    let key = [
        jd_tt.to_bits(),
        flags as u64,
        nutation.0.to_bits(),
        nutation.1.to_bits(),
        nutation_rate.0.to_bits(),
        nutation_rate.1.to_bits(),
    ];
    if let Some(frame) = crate::pure_cache::analytic_frame(key) {
        return Ok(analytic_frame_from_storage(frame));
    }
    let result = analytic_output_frame_words_uncached(jd_tt, flags, nutation, nutation_rate)?;
    crate::pure_cache::retain_analytic_frame(key, analytic_frame_to_storage(result));
    Ok(result)
}

/// Store both high/low words directly: no rounding, normalization or arithmetic.
fn analytic_frame_to_storage(
    value: (WordMatrix3, WordMatrix3),
) -> crate::pure_cache::AnalyticFrame {
    crate::pure_cache::AnalyticFrame {
        matrix: value.0.map(|row| row.map(|word| [word.hi, word.lo])),
        rate: value.1.map(|row| row.map(|word| [word.hi, word.lo])),
    }
}

fn analytic_frame_from_storage(
    value: crate::pure_cache::AnalyticFrame,
) -> (WordMatrix3, WordMatrix3) {
    use super::derivative::AngularWord as Word;
    let decode = |matrix: [[[f64; 2]; 3]; 3]| matrix.map(|row| row.map(|[hi, lo]| Word { hi, lo }));
    (decode(value.matrix), decode(value.rate))
}

/// Compute the same analytic frame after the wrapper's ordered finite checks.
/// Misses execute outside the cache lock, including all pure frame failures.
fn analytic_output_frame_words_uncached(
    jd_tt: f64,
    flags: i32,
    nutation: (f64, f64),
    nutation_rate: (f64, f64),
) -> Result<(WordMatrix3, WordMatrix3)> {
    use super::derivative::AngularWord as Word;
    let zero = Word::from(0.0);
    let one = Word::from(1.0);
    if flags & FLG_J2000 != 0 {
        let bias = if flags & FLG_ICRS != 0 {
            [[one, zero, zero], [zero, one, zero], [zero, zero, one]]
        } else {
            FRAME_BIAS_J2000.map(|row| row.map(Word::from))
        };
        if flags & FLG_EQUATORIAL != 0 {
            return Ok((bias, [[zero; 3]; 3]));
        }
        // These are the actual two scalar physical angle definitions, including
        // the distinct existing ICRS branch. No angle or residual is fitted.
        let angle = if flags & FLG_ICRS != 0 {
            vondrak_mean_obliquity_deg(FRAME_J2000)?.to_radians()
        } else {
            OBLIQUITY_J2000_REF_RAD
        };
        let (sine, cosine) = obliquity_sin_cos_words(angle)?;
        let rotation = [
            [one, zero, zero],
            [zero, cosine, sine],
            [zero, -sine, cosine],
        ];
        let matrix = std::array::from_fn(|i| {
            std::array::from_fn(|j| (0..3).fold(zero, |sum, k| sum + rotation[i][k] * bias[k][j]))
        });
        return Ok((matrix, [[zero; 3]; 3]));
    }
    let (_, rate) =
        analytic_output_frame_matrix_rate_with_nutation(jd_tt, flags, nutation, nutation_rate)?;
    let matrix = if flags & FLG_EQUATORIAL != 0 {
        equatorial_frame_words(jd_tt, flags, nutation)?
    } else {
        ecliptic_frame_words(jd_tt, flags, nutation.0)?
    };
    Ok((matrix, rate.map(|row| row.map(Word::from))))
}

/// Direct-TT central ecliptic frame before pole, basis and product rounding.
/// Every original Vondrak term and binary64 constant is retained. The existing
/// century operand and scalar position path are unchanged; only the analytic
/// velocity's central linear map carries its arithmetic residuals onward.
fn ecliptic_frame_words(jd_tt: f64, flags: i32, dpsi: f64) -> Result<WordMatrix3> {
    use super::derivative::AngularWord as Word;
    let zero = Word::from(0.0);
    let one = Word::from(1.0);
    let t = Word::from(super::corrections::julian_centuries(jd_tt));
    let (p, q) = vondrak_pair_words(
        t,
        &VONDRAK_PQPER_PERIODS,
        &VONDRAK_PQPER_P_COS,
        &VONDRAK_PQPER_P_SIN,
        &VONDRAK_PQPER_Q_COS,
        &VONDRAK_PQPER_Q_SIN,
        &VONDRAK_PQPOL,
    )?;
    let (x, y) = vondrak_pair_words(
        t,
        &VONDRAK_XYPER_PERIODS,
        &VONDRAK_XYPER_X_COS,
        &VONDRAK_XYPER_X_SIN,
        &VONDRAK_XYPER_Y_COS,
        &VONDRAK_XYPER_Y_SIN,
        &VONDRAK_XYPOL,
    )?;
    // Preserve the pole model's existing nonpositive-radicand convention.
    let pole_z = |a: Word, b: Word| {
        let squared = one - a * a - b * b;
        if squared.rounded() > 0.0 {
            squared.sqrt()
        } else {
            zero
        }
    };
    let z = pole_z(p, q);
    let (s, c) = frame_sin_cos_words(Word::from(VONDRAK_EPS0_RAD))?;
    let ecliptic = [p, -q * c - z * s, -q * s + z * c];
    let equator = [x, y, pole_z(x, y)];
    let cross = |a: [Word; 3], b: [Word; 3]| {
        std::array::from_fn(|i| {
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            a[j] * b[k] - a[k] * b[j]
        })
    };
    let equinox = cross(equator, ecliptic);
    let length = equinox.into_iter().fold(zero, |sum, v| sum + v * v).sqrt();
    if length.rounded() == 0.0 || !length.rounded().is_finite() {
        return Err(Error::MathDomain {
            value: length.rounded(),
        });
    }
    let equinox = equinox.map(|v| v / length);
    let mut basis = [equinox, cross(ecliptic, equinox), ecliptic];
    if flags & FLG_ICRS == 0 {
        basis = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                (0..3).fold(zero, |sum, k| {
                    sum + basis[i][k] * Word::from(FRAME_BIAS_LTPB[k][j])
                })
            })
        });
    }
    if flags & FLG_NONUT == 0 {
        let (s, c) = frame_sin_cos_words(Word::from(dpsi))?;
        basis = [
            std::array::from_fn(|i| c * basis[0][i] - s * basis[1][i]),
            std::array::from_fn(|i| s * basis[0][i] + c * basis[1][i]),
            basis[2],
        ];
    }
    for row in basis {
        check_finite_vec3(row.map(Word::rounded))?;
    }
    Ok(basis)
}

/// Full periodic/polynomial pole pair in radians, at the original TT centuries.
/// This is the published physical pair, not a shortened series or expression
/// evaluator. Elementary sine/cosine uses the bounded frame-angle helper.
pub(super) fn vondrak_pair_words(
    t: super::derivative::AngularWord,
    periods: &[f64],
    first_cos: &[f64],
    first_sin: &[f64],
    second_cos: &[f64],
    second_sin: &[f64],
    polynomial: &[(f64, f64)],
) -> Result<(
    super::derivative::AngularWord,
    super::derivative::AngularWord,
)> {
    use super::derivative::AngularWord as Word;
    let mut first = Word::from(0.0);
    let mut second = Word::from(0.0);
    for (i, &period) in periods.iter().enumerate() {
        let angle = Word::from(std::f64::consts::TAU) * t / Word::from(period);
        let (s, c) = frame_sin_cos_words(angle)?;
        first = first + c * Word::from(first_cos[i]) + s * Word::from(first_sin[i]);
        second = second + c * Word::from(second_cos[i]) + s * Word::from(second_sin[i]);
    }
    let mut power = Word::from(1.0);
    for &(a, b) in polynomial {
        first = first + Word::from(a) * power;
        second = second + Word::from(b) * power;
        power = power * t;
    }
    let scale = Word::from(ARCSEC_TO_RAD);
    Ok((first * scale, second * scale))
}

/// Central equatorial frame before elementary-rotation/product rounding.
/// Precession retains its actual Julian-epoch realization; the pole-angle
/// obliquity retains its distinct direct-TT realization. The selected feed's
/// nutation inputs are exact binary64 operands. Scalar positions are untouched.
fn equatorial_frame_words(jd_tt: f64, flags: i32, nutation: (f64, f64)) -> Result<WordMatrix3> {
    use super::derivative::AngularWord as Word;
    let p = vondrak_precession_matrix(jd_tt, flags & FLG_ICRS == 0)?.map(|row| row.map(Word::from));
    if flags & (FLG_NONUT | FLG_SIDEREAL) != 0 {
        return Ok(p);
    }
    let zero = Word::from(0.0);
    let one = Word::from(1.0);
    let multiply = |a: WordMatrix3, b: WordMatrix3| -> WordMatrix3 {
        std::array::from_fn(|i| {
            std::array::from_fn(|j| (0..3).fold(zero, |sum, k| sum + a[i][k] * b[k][j]))
        })
    };
    let x_rotation = |angle| -> Result<WordMatrix3> {
        let (s, c) = frame_sin_cos_words(angle)?;
        Ok([[one, zero, zero], [zero, c, s], [zero, -s, c]])
    };
    let eps = Word::from(vondrak_mean_obliquity_rad(jd_tt)?);
    let (s, c) = frame_sin_cos_words(Word::from(-nutation.0))?;
    let b = [[c, s, zero], [-s, c, zero], [zero, zero, one]];
    let a = x_rotation(eps)?;
    let c = x_rotation(-eps - Word::from(nutation.1))?;
    Ok(multiply(multiply(c, multiply(b, a)), p))
}

/// Sine/cosine of the existing obliquity angle, with two binary64 words.
/// Both fixed-frame obliquities satisfy |angle| < 1 radian. Twenty-four updates
/// leave an alternating-series remainder below 1/50! < 4e-65 on this domain,
/// far below the two-word arithmetic error. This evaluates elementary functions;
/// it does not truncate any astronomical precession or nutation model.
fn obliquity_sin_cos_words(
    angle: f64,
) -> Result<(
    super::derivative::AngularWord,
    super::derivative::AngularWord,
)> {
    use super::derivative::AngularWord as Word;
    if !angle.is_finite() || angle.abs() > 1.0 {
        return Err(Error::MathDomain { value: angle });
    }
    frame_sin_cos_words(Word::from(angle))
}

/// Elementary rotations with a two-word angle and bounded Payne-Hanek reduction.
/// The 1280-bit reciprocal covers every finite binary64 exponent with more
/// than 250 guard bits. Only a fixed 21-word integer product is needed; no
/// retained allocation or approximate large-angle fallback is introduced.
pub(super) fn frame_sin_cos_words(
    angle: super::derivative::AngularWord,
) -> Result<(
    super::derivative::AngularWord,
    super::derivative::AngularWord,
)> {
    use super::derivative::AngularWord as Word;
    if !angle.hi.is_finite() || !angle.lo.is_finite() {
        return Err(Error::MathDomain {
            value: angle.rounded(),
        });
    }
    // Preserve the mathematical sign of an exact zero input. The low word
    // has the same sign so that the public rounded sum retains signed zero.
    if angle.hi == 0.0 && angle.lo == 0.0 {
        return Ok((
            Word {
                hi: angle.hi,
                lo: angle.hi,
            },
            Word::from(1.0),
        ));
    }
    let (quadrant_hi, hi) = reduce_frame_angle(angle.hi);
    let (quadrant_lo, lo) = reduce_frame_angle(angle.lo);
    let x = hi + lo;
    let negative_square = -(x * x);
    let mut sine = x;
    let mut cosine = Word::from(1.0);
    let mut sine_term = x;
    let mut cosine_term = Word::from(1.0);
    // The two reduced words sum to at most two radians. The first omitted
    // cosine term is bounded by 2^50/50! < 4e-50, below two-word roundoff.
    for n in 1..=24 {
        let k = 2 * n;
        sine_term = sine_term * negative_square / Word::from((k * (k + 1)) as f64);
        cosine_term = cosine_term * negative_square / Word::from(((k - 1) * k) as f64);
        sine = sine + sine_term;
        cosine = cosine + cosine_term;
    }
    Ok(match (quadrant_hi + quadrant_lo) & 3 {
        0 => (sine, cosine),
        1 => (cosine, -sine),
        2 => (-sine, -cosine),
        _ => (-cosine, sine),
    })
}

/// Reduce one exact binary64 operand to a quarter-turn and a small residual.
fn reduce_frame_angle(angle: f64) -> (u32, super::derivative::AngularWord) {
    use super::derivative::AngularWord as Word;
    if angle.abs() <= 1.0 {
        return (0, Word::from(angle));
    }
    // floor((2/pi) * 2^1280), little-endian 64-bit limbs. This is a
    // mathematical constant, not an ephemeris/model or fitted coefficient.
    const TWO_OVER_PI: [u64; 20] = [
        0xf0cfbc209af4361d,
        0x56033046fc7b6bab,
        0x6bfb5fb11f8d5d08,
        0x3d0739f78a5292ea,
        0x7527bac7ebe5f17b,
        0x4f463f669e5fea2d,
        0x6d367ecf27cb09b7,
        0xef2f118b5a0a6d1f,
        0x1ff897ffde05980f,
        0x9c845f8bbdf9283b,
        0x3991d639835339f4,
        0xe99c7026b45f7e41,
        0xe88235f52ebb4484,
        0xfe1deb1cb129a73e,
        0x06492eea09d1921c,
        0xb7246e3a424dd2e0,
        0xfe5163abdebbc561,
        0xdb6295993c439041,
        0xfc2757d1f534ddc0,
        0xa2f9836e4e441529,
    ];
    let bits = angle.abs().to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32 - 1023;
    let mantissa = (bits & ((1_u64 << 52) - 1)) | (1_u64 << 52);
    let mut product = [0_u64; 21];
    let mut carry = 0_u128;
    for (i, limb) in TWO_OVER_PI.into_iter().enumerate() {
        let value = u128::from(limb) * u128::from(mantissa) + carry;
        product[i] = value as u64;
        carry = value >> 64;
    }
    product[20] = carry as u64;
    let shift = (1280 + 52 - exponent) as usize;
    let bit = |index: usize| (product[index / 64] >> (index % 64)) & 1;
    let extract = |start: usize| -> u64 {
        (0..53).fold(0, |value, offset| value | (bit(start + offset) << offset))
    };
    let round_up = bit(shift - 1);
    let quadrant = ((bit(shift) + 2 * bit(shift + 1) + round_up) & 3) as u32;
    let mut fraction = Word::from(extract(shift - 53) as f64 * 2.0_f64.powi(-53))
        + Word::from(extract(shift - 106) as f64 * 2.0_f64.powi(-106))
        + Word::from(extract(shift - 159) as f64 * 2.0_f64.powi(-159));
    if round_up != 0 {
        fraction = fraction - Word::from(1.0);
    }
    // Mathematical pi/2 split into its nearest binary64 value and residual.
    let half_pi = Word {
        hi: std::f64::consts::FRAC_PI_2,
        lo: 6.123233995736766e-17,
    };
    let reduced = fraction * half_pi;
    if angle.is_sign_negative() {
        ((4 - quadrant) & 3, -reduced)
    } else {
        (quadrant, reduced)
    }
}

#[cfg(test)]
mod compensated_frame_tests {
    use super::super::derivative::AngularWord as Word;
    use super::*;

    #[test]
    fn original_binary64_obliquities_have_independent_two_word_trigonometry() {
        // Mathematical sine/cosine at these exact binary64 inputs, independently
        // evaluated to 100 decimal digits. These values are test references only.
        type Case = (u64, [u64; 2], [u64; 2]);
        let cases: &[Case] = &[
            (
                4601041177052405169,
                [4600837332757394311, 13573441341150731953],
                [4606439164020531830, 4363111936813157835],
            ),
            (
                13824413213907180977,
                [13824209369612170119, 4350069304295956145],
                [4606439164020531830, 4363111936813157835],
            ),
            (0, [0, 0], [4607182418800017408, 0]),
            (
                4602678819172646912,
                [4602308182625945072, 13571467409308189309],
                [4606079780542709072, 13585268877995206696],
            ),
            (
                13826050856027422720,
                [13825680219480720880, 4348095372453413501],
                [4606079780542709072, 13585268877995206696],
            ),
            (
                4607182418800017408,
                [4605754516372524270, 4341579394754863752],
                [4603041830072026764, 13586077967046624444],
            ),
            (
                13830554455654793216,
                [13829126553227300078, 13564951431609639560],
                [4603041830072026764, 13586077967046624444],
            ),
            (
                4427038433705197568,
                [4427038433705197568, 13278112901364012373],
                [4607182418800017408, 13465762885837783040],
            ),
        ];
        let expand = |pair: [u64; 2]| Word {
            hi: f64::from_bits(pair[0]),
            lo: f64::from_bits(pair[1]),
        };
        for &(angle, sine, cosine) in cases {
            let (s, c) = obliquity_sin_cos_words(f64::from_bits(angle)).unwrap();
            for (actual, expected) in [(s, expand(sine)), (c, expand(cosine))] {
                assert_eq!(actual.rounded().to_bits(), expected.hi.to_bits());
                assert!((actual - expected).rounded().abs() <= 2.0_f64.powi(-99));
            }
        }
    }

    #[test]
    fn fixed_frame_products_retain_physical_rotation_residuals_and_zero_rate() {
        // Rx(original obliquity) times the exact stored bias coefficients,
        // independently evaluated to 100 decimal digits, then split into words.
        let expected: [[[u64; 2]; 3]; 3] = [
            [
                [4607182418800017355, 0],
                [13723312729344079124, 0],
                [4500679531294428041, 0],
            ],
            [
                [4495060083190529556, 4255499033251274199],
                [4606439163902084608, 13584090899711104810],
                [4600837333303796994, 13580758646650283867],
            ],
            [
                [13725676689914655412, 13471484629386064509],
                [13824209370158572712, 13581449097177753553],
                [4606439163902084582, 13580818740307725514],
            ],
        ];
        let (matrix, rate) = analytic_output_frame_words_with_nutation(
            FRAME_J2000,
            FLG_J2000,
            (0.0, 0.0),
            (0.0, 0.0),
        )
        .unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let reference = Word {
                    hi: f64::from_bits(expected[i][j][0]),
                    lo: f64::from_bits(expected[i][j][1]),
                };
                assert!((matrix[i][j] - reference).rounded().abs() <= 2.0_f64.powi(-99));
                assert_eq!(rate[i][j].hi, 0.0);
                assert_eq!(rate[i][j].lo, 0.0);
            }
        }
        for flags in [
            FLG_J2000 | FLG_EQUATORIAL,
            FLG_J2000 | FLG_EQUATORIAL | FLG_ICRS,
        ] {
            let (_, rate) = analytic_output_frame_words_with_nutation(
                FRAME_J2000 - 36525.0,
                flags,
                (1e-5, -2e-5),
                (1e-7, 2e-7),
            )
            .unwrap();
            assert!(
                rate.into_iter()
                    .flatten()
                    .all(|word| word.hi == 0.0 && word.lo == 0.0)
            );
        }
    }
}

#[cfg(test)]
mod date_frame_word_tests {
    use super::super::derivative::AngularWord as Word;
    use super::*;
    #[test]
    fn rotation_words_cover_large_finite_angles_and_original_value_residuals() {
        // Independent 1000-decimal mathematical sin/cos of exact input words.
        type Case = (u64, u64, [u64; 2], [u64; 2]);
        let cases: &[Case] = &[
            (
                4609434218613702656,
                0,
                [4607159855645224331, 13578573826839128198],
                [4589761573224315303, 4346251512829300350],
            ),
            (
                4614256656552045848,
                0,
                [4368955796522032135, 13343912218604703677],
                [13830554455654793216, 4126273412488912200],
            ),
            (
                13837628693406821656,
                0,
                [13592327833376807943, 4120540181749927869],
                [13830554455654793216, 4126273412488912200],
            ),
            (
                4621819117588971520,
                0,
                [13826447362944618322, 13584672858912522535],
                [13829104940851424031, 13578440279811922325],
            ),
            (
                13860109328209412096,
                0,
                [4602736155770506680, 4344887152559268091],
                [4605942297449095135, 4362014470795003368],
            ),
            (
                4696837146684686336,
                0,
                [13823348579196524647, 13579026272317008313],
                [4606612732610269996, 4362504945544226432],
            ),
            (
                4906019910204099648,
                0,
                [13827359163295818832, 13586495920377590171],
                [4605056453202808125, 4362038214025942328],
            ),
            (
                6103021453049119613,
                0,
                [13823900616546751321, 4354270060022202799],
                [4606504395019403765, 4355894209148188583],
            ),
            (
                9218868437227405311,
                0,
                [4572370766535020183, 13552030231331168479],
                [13830554344771005301, 13582297432791740657],
            ),
            (
                18442240474082181119,
                0,
                [13795742803389795991, 4328658194476392671],
                [13830554344771005301, 13582297432791740657],
            ),
            (
                4602678819172646912,
                4158027847206421152,
                [4602308182625945072, 13571467409308188169],
                [4606079780542709072, 13585268877995206774],
            ),
            (
                4906019910204099648,
                4661225614328463360,
                [13830311845450013456, 13584618840899164519],
                [4597473793896798811, 13558596909142069089],
            ),
            (
                15326393489903895421,
                5848587129000908155,
                [13830548960268465780, 4363583265611238455],
                [4585194235908298840, 13561509106073662671],
            ),
        ];
        for &(hi, lo, sine, cosine) in cases {
            let angle = Word {
                hi: f64::from_bits(hi),
                lo: f64::from_bits(lo),
            };
            let (s, c) = frame_sin_cos_words(angle).unwrap();
            for (actual, expected) in [(s, sine), (c, cosine)] {
                let expected = Word {
                    hi: f64::from_bits(expected[0]),
                    lo: f64::from_bits(expected[1]),
                };
                assert_eq!(
                    actual.rounded().to_bits(),
                    expected.rounded().to_bits(),
                    "angle=({}, {})",
                    angle.hi,
                    angle.lo
                );
                assert!(
                    (actual - expected).rounded().abs() <= 2.0_f64.powi(-98),
                    "angle=({}, {})",
                    angle.hi,
                    angle.lo
                );
            }
        }
    }

    #[test]
    fn quarter_turn_boundaries_subnormals_and_signed_zero_keep_physical_continuity() {
        // Adjacent binary64 values around reduction and quadrant boundaries,
        // evaluated independently at 1000 decimal digits. The error bound is
        // absolute, since the exact cosine may approach zero at a quarter turn.
        type Case = (u64, [u64; 2], [u64; 2]);
        let cases: &[Case] = &[
            (
                4607182418800017407,
                [4605754516372524269, 4363550326723068176],
                [4603041830072026764, 4362414342797186335],
            ),
            (
                13830554455654793215,
                [13829126553227300077, 13586922363577843984],
                [4603041830072026764, 4362414342797186335],
            ),
            (
                4607182418800017408,
                [4605754516372524270, 4341579394754863752],
                [4603041830072026764, 13586077967046624444],
            ),
            (
                13830554455654793216,
                [13829126553227300078, 13564951431609639560],
                [4603041830072026764, 13586077967046624444],
            ),
            (
                4607182418800017409,
                [4605754516372524271, 4352935053201875278],
                [4603041830072026762, 13577399600525730267],
            ),
            (
                13830554455654793217,
                [13829126553227300079, 13576307090056651086],
                [4603041830072026762, 13577399600525730267],
            ),
            (
                4609753056924675351,
                [4607182418800017408, 13360502858171208664],
                [4374237177576756994, 13353720002776062205],
            ),
            (
                13833125093779451159,
                [13830554455654793216, 4137130821316432856],
                [4374237177576756994, 13353720002776062205],
            ),
            (
                4609753056924675352,
                [4607182418800017408, 13340638250088947016],
                [4364452196894661639, 13339408618977333180],
            ),
            (
                13833125093779451160,
                [13830554455654793216, 4117266213234171208],
                [4364452196894661639, 13339408618977333180],
            ),
            (
                4609753056924675353,
                [4607182418800017408, 13353393706361483170],
                [13593883396231516669, 4128706643053661705],
            ),
            (
                13833125093779451161,
                [13830554455654793216, 4130021669506707362],
                [13593883396231516669, 4128706643053661705],
            ),
            (
                4614256656552045847,
                [4378740777204127490, 13358223602403432705],
                [13830554455654793216, 4146138020571173848],
            ),
            (
                13837628693406821655,
                [13602112814058903298, 4134851565548656897],
                [13830554455654793216, 4146138020571173848],
            ),
            (
                4614256656552045848,
                [4368955796522032135, 13343912218604703677],
                [13830554455654793216, 4126273412488912200],
            ),
            (
                13837628693406821656,
                [13592327833376807943, 4120540181749927869],
                [13830554455654793216, 4126273412488912200],
            ),
            (
                4614256656552045849,
                [13598386995858887165, 4133210242681032203],
                [13830554455654793216, 4139028868761448354],
            ),
            (
                13837628693406821657,
                [4375014959004111357, 13356582279535808011],
                [13830554455654793216, 4139028868761448354],
            ),
            (
                4616991696741409233,
                [13830554455654793216, 4154374413117652938],
                [13606305928836272449, 13362895663598198064],
            ),
            (
                13840363733596185041,
                [4607182418800017408, 13377746449972428746],
                [13606305928836272449, 13362895663598198064],
            ),
            (
                4616991696741409234,
                [13830554455654793216, 4131461981296896241],
                [13594811712176818698, 13349890249418271257],
            ),
            (
                13840363733596185042,
                [4607182418800017408, 13354834018151672049],
                [13594811712176818698, 13349890249418271257],
            ),
            (
                4616991696741409235,
                [13830554455654793216, 4148978333399676298],
                [4380139528331484541, 4137376919173612911],
            ),
            (
                13840363733596185043,
                [4607182418800017408, 13372350370254452106],
                [4380139528331484541, 4137376919173612911],
            ),
            (
                4618760256179416343,
                [13606616413686273794, 4139355165176027410],
                [4607182418800017408, 13378517256680690648],
            ),
            (
                13842132293034192151,
                [4383244376831497986, 13362727202030803218],
                [4607182418800017408, 13378517256680690648],
            ),
            (
                4618760256179416344,
                [13596831433004178439, 4125043781377298367],
                [4607182418800017408, 13358652648598429000],
            ),
            (
                13842132293034192152,
                [4373459396149402631, 13348415818232074175],
                [4607182418800017408, 13358652648598429000],
            ),
            (
                4618760256179416345,
                [4379518558631481853, 13361085879163178513],
                [4607182418800017408, 13371408104870965154],
            ),
            (
                13842132293034192153,
                [13602890595486257661, 4137713842308402705],
                [4607182418800017408, 13371408104870965154],
            ),
            (
                4744861045745790231,
                [13732717203252647593, 13483818939279547000],
                [4607182418800016991, 4361102888008778318],
            ),
            (
                13968233082600566039,
                [4509345166397871785, 4260446902424771192],
                [4607182418800016991, 4361102888008778318],
            ),
            (
                4744861045745790232,
                [13722932222570552323, 13480164395161165322],
                [4607182418800017389, 13586783525316022655],
            ),
            (
                13968233082600566040,
                [4499560185715776515, 4256792358306389514],
                [4607182418800017389, 13586783525316022655],
            ),
            (
                4744861045745790233,
                [4505619348197855708, 4256623339967850478],
                [4607182418800017274, 13583321823941962945],
            ),
            (
                13968233082600566041,
                [13728991385052631516, 13479995376822626286],
                [4607182418800017274, 13583321823941962945],
            ),
            (1, [1, 0], [4607182418800017408, 0]),
            (
                9223372036854775809,
                [9223372036854775809, 0],
                [4607182418800017408, 0],
            ),
            (
                4503599627370495,
                [4503599627370495, 0],
                [4607182418800017408, 0],
            ),
            (
                9227875636482146303,
                [9227875636482146303, 0],
                [4607182418800017408, 0],
            ),
        ];
        for &(angle, sine, cosine) in cases {
            let (s, c) = frame_sin_cos_words(Word::from(f64::from_bits(angle))).unwrap();
            for (actual, expected) in [(s, sine), (c, cosine)] {
                let expected = Word {
                    hi: f64::from_bits(expected[0]),
                    lo: f64::from_bits(expected[1]),
                };
                assert!((actual - expected).rounded().abs() <= 2.0_f64.powi(-98));
                assert!(
                    actual
                        .rounded()
                        .to_bits()
                        .abs_diff(expected.rounded().to_bits())
                        <= 1,
                    "angle={}",
                    f64::from_bits(angle)
                );
            }
        }
        for zero in [0.0_f64, -0.0_f64] {
            let (sine, cosine) = frame_sin_cos_words(Word::from(zero)).unwrap();
            assert_eq!(sine.rounded().to_bits(), zero.to_bits());
            assert_eq!(cosine.rounded(), 1.0);
        }
        for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(frame_sin_cos_words(Word::from(invalid)).is_err());
        }
    }

    #[test]
    fn date_frame_words_preserve_unrotated_precession_and_long_term_domains() {
        for tt in [
            2451545.0,
            2451545.0 - 200000.0 * 365.25,
            2451545.0 + 200000.0 * 365.25,
        ] {
            for flags in [FLG_EQUATORIAL, FLG_EQUATORIAL | FLG_ICRS] {
                let matrix = equatorial_frame_words(tt, flags, (0.0002, -0.0001)).unwrap();
                assert!(
                    matrix
                        .into_iter()
                        .flatten()
                        .all(|v| v.rounded().is_finite())
                );
                for convention in [FLG_NONUT, FLG_SIDEREAL] {
                    let actual =
                        equatorial_frame_words(tt, flags | convention, (0.2, -0.1)).unwrap();
                    let expected = vondrak_precession_matrix(tt, flags & FLG_ICRS == 0).unwrap();
                    for i in 0..3 {
                        for j in 0..3 {
                            assert_eq!(actual[i][j].hi.to_bits(), expected[i][j].to_bits());
                            assert_eq!(actual[i][j].lo, 0.0);
                        }
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod ecliptic_words_tests {
    use super::super::derivative::AngularWord as Word;
    use super::*;

    #[test]
    fn compensated_ecliptic_basis_preserves_physical_orthogonality() {
        for t in [-2000.0, -100.0, -0.456, 0.0, 0.3, 20.0, 2000.0] {
            let tt = FRAME_J2000 + 36525.0 * t;
            for flags in [FLG_ICRS | FLG_NONUT, FLG_ICRS] {
                let matrix = ecliptic_frame_words(tt, flags, 0.0002).unwrap();
                for i in 0..3 {
                    for j in 0..3 {
                        let dot = (0..3)
                            .fold(Word::from(0.0), |sum, k| sum + matrix[i][k] * matrix[j][k]);
                        let target = Word::from(if i == j { 1.0 } else { 0.0 });
                        assert!(
                            (dot - target).rounded().abs() < 2e-28,
                            "t={t}, flags={flags}, rows={i}/{j}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn compensated_central_ecliptic_does_not_change_original_frame_rate() {
        for tt in [2415020.5, 2451545.0, 2460310.5, 7900000.5] {
            for flags in [FLG_SPEED, FLG_SPEED | FLG_NONUT, FLG_SPEED | FLG_ICRS] {
                let nutation = (0.0002, -0.0001);
                let rate = (1e-7, -2e-8);
                let (_, original) =
                    analytic_output_frame_matrix_rate_with_nutation(tt, flags, nutation, rate)
                        .unwrap();
                let (_, expanded) =
                    analytic_output_frame_words_with_nutation(tt, flags, nutation, rate).unwrap();
                for i in 0..3 {
                    for j in 0..3 {
                        assert_eq!(expanded[i][j].hi.to_bits(), original[i][j].to_bits());
                        assert_eq!(expanded[i][j].lo, 0.0);
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod analytic_retention_tests {
    use super::super::derivative::AngularWord as Word;
    use super::*;

    fn bits(value: (WordMatrix3, WordMatrix3)) -> [[[[u64; 2]; 3]; 3]; 2] {
        [value.0, value.1]
            .map(|matrix| matrix.map(|row| row.map(|v| [v.hi.to_bits(), v.lo.to_bits()])))
    }

    #[test]
    fn analytic_retention_roundtrip_keeps_every_high_and_low_bit() {
        let patterns = [
            0,
            1_u64 << 63,
            1,
            (1_u64 << 63) | 1,
            0x0010000000000000,
            0x7fefffffffffffff,
            0xffefffffffffffff,
            0x7ff0000000000000,
            0xfff0000000000000,
            0x7ff0000000000001,
            0x7ff8000000001234,
            0xfff8000000001234,
            1.0_f64.to_bits(),
            (-1.0_f64).to_bits(),
        ];
        for offset in 0..patterns.len() {
            let matrix = std::array::from_fn(|i| {
                std::array::from_fn(|j| Word {
                    hi: f64::from_bits(patterns[(offset + 3 * i + j) % patterns.len()]),
                    lo: f64::from_bits(patterns[(offset + 3 * i + j + 1) % patterns.len()]),
                })
            });
            let rate = matrix.map(|row| row.map(|v| Word { hi: v.lo, lo: v.hi }));
            let original = (matrix, rate);
            assert_eq!(
                bits(analytic_frame_from_storage(analytic_frame_to_storage(
                    original
                ))),
                bits(original)
            );
        }
    }

    #[test]
    fn analytic_retention_hits_match_fresh_words_across_every_frame_branch() {
        let mut nonzero_low_word = false;
        for tt in [
            2451545.0,
            f64::from_bits(2451545.0_f64.to_bits() + 1),
            2460310.5,
            FRAME_J2000 - 200000.0 * 365.25,
            FRAME_J2000 + 200000.0 * 365.25,
        ] {
            for mask in 0..32 {
                let choices = [FLG_EQUATORIAL, FLG_ICRS, FLG_J2000, FLG_NONUT, FLG_SIDEREAL];
                let flags = choices
                    .into_iter()
                    .enumerate()
                    .fold(FLG_SPEED, |flags, (bit, value)| {
                        flags | if mask & (1 << bit) != 0 { value } else { 0 }
                    });
                for (angles, rates) in [
                    ((0.0, -0.0), (-0.0, 0.0)),
                    ((-0.0, 0.0), (0.0, -0.0)),
                    ((0.0002, -0.0001), (1e-7, -2e-8)),
                ] {
                    let fresh =
                        analytic_output_frame_words_uncached(tt, flags, angles, rates).unwrap();
                    nonzero_low_word |= fresh.0.into_iter().flatten().any(|v| v.lo != 0.0);
                    for _ in 0..2 {
                        let actual =
                            analytic_output_frame_words_with_nutation(tt, flags, angles, rates)
                                .unwrap();
                        assert_eq!(bits(actual), bits(fresh));
                    }
                }
            }
        }
        assert!(nonzero_low_word);
    }

    #[test]
    fn analytic_retention_concurrent_replacement_keeps_fresh_matrix_words() {
        std::thread::scope(|scope| {
            for thread in 0..8 {
                scope.spawn(move || {
                    for index in 0..48 {
                        let tt = FRAME_J2000 + (index % 24) as f64;
                        let flags = FLG_SPEED
                            | if (index + thread) % 2 == 0 {
                                FLG_ICRS
                            } else {
                                FLG_EQUATORIAL
                            };
                        let angles = (0.0001, -0.0002);
                        let rate = (1e-7, -2e-8);
                        let fresh =
                            analytic_output_frame_words_uncached(tt, flags, angles, rate).unwrap();
                        let actual =
                            analytic_output_frame_words_with_nutation(tt, flags, angles, rate)
                                .unwrap();
                        assert_eq!(bits(actual), bits(fresh));
                    }
                });
            }
        });
    }

    #[test]
    fn analytic_retention_admission_precedes_even_an_existing_invalid_key() {
        let value = analytic_frame_to_storage(
            analytic_output_frame_words_uncached(FRAME_J2000, FLG_J2000, (0.0, 0.0), (0.0, 0.0))
                .unwrap(),
        );
        for bad in [
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::from_bits(0x7ff800000000abcd),
        ] {
            for index in 0..5 {
                let mut input = [FRAME_J2000, 0.0, 0.0, 0.0, 0.0];
                input[index] = bad;
                let key = [
                    input[0].to_bits(),
                    FLG_J2000 as u64,
                    input[1].to_bits(),
                    input[2].to_bits(),
                    input[3].to_bits(),
                    input[4].to_bits(),
                ];
                crate::pure_cache::retain_analytic_frame(key, value);
                let result = analytic_output_frame_words_with_nutation(
                    input[0],
                    FLG_J2000,
                    (input[1], input[2]),
                    (input[3], input[4]),
                );
                match result {
                    Err(Error::NonFiniteArgument { value }) => {
                        assert_eq!(value.to_bits(), bad.to_bits())
                    }
                    _ => panic!("finite admission must precede retained results"),
                }
            }
        }
        // The first invalid operand keeps its original place in validation.
        assert!(matches!(analytic_output_frame_words_with_nutation(
            f64::INFINITY, FLG_J2000, (f64::NEG_INFINITY, 0.0), (f64::NAN, 0.0)),
            Err(Error::NonFiniteArgument { value }) if value == f64::INFINITY));
        for _ in 0..2 {
            let uncached =
                analytic_output_frame_words_uncached(1e308, 0, (0.0, 0.0), (0.0, 0.0)).map(bits);
            let cached =
                analytic_output_frame_words_with_nutation(1e308, 0, (0.0, 0.0), (0.0, 0.0))
                    .map(bits);
            assert_eq!(format!("{cached:?}"), format!("{uncached:?}"));
        }
    }

    /// Straight-line physical states expose exact velocities and zero
    /// acceleration. Mutable fields inject provider faults, not model changes.
    struct LinearProvider {
        events: std::cell::RefCell<Vec<&'static str>>,
        fail_rate: std::cell::Cell<bool>,
        fail_support: std::cell::Cell<bool>,
        stored_rate: std::cell::Cell<bool>,
    }
    impl StateProvider for LinearProvider {
        fn provider_name(&self) -> &'static str {
            "analytic-frame-retention-test"
        }
        fn has_body(&self, body: i32) -> bool {
            matches!(body, SUN | EARTH | MARS)
        }
        fn body_coverage(&self, body: i32) -> Option<(f64, f64)> {
            self.has_body(body)
                .then_some((FRAME_J2000 - 10.0, FRAME_J2000 + 10.0))
        }
        fn serving_coord_type(&self, body: i32) -> Option<u32> {
            self.has_body(body).then_some(COORD_ICRS_BARY)
        }
        fn serving_entry_snapshot(&self, body: i32) -> Result<EntrySnapshot> {
            if self.has_body(body) {
                Ok(EntrySnapshot {
                    coord_type: COORD_ICRS_BARY,
                })
            } else {
                Err(Error::LebBodyNotFound { body_id: body })
            }
        }
        fn eval_raw(&self, body: i32, tt: f64) -> Result<(Vec<f64>, Vec<f64>)> {
            self.eval_raw_two_part(body, TwoPartEpoch::split(tt)?)
        }
        fn eval_raw_two_part(&self, body: i32, tt: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
            let t = (tt.hi - FRAME_J2000) + tt.lo;
            match body {
                SUN => Ok((vec![0.0; 3], vec![0.0; 3])),
                EARTH => Ok((vec![t * 0.001, 0.0, 0.0], vec![0.001, 0.0, 0.0])),
                MARS => Ok((vec![2.0 + t * 0.003, 0.2, 0.3], vec![0.003, 0.0, 0.0])),
                _ => Err(Error::LebBodyNotFound { body_id: body }),
            }
        }
        fn has_analytic_kinematics(&self) -> bool {
            true
        }
        fn geometric_acceleration_tt(&self, _: i32, _: f64) -> Result<Option<[f64; 3]>> {
            Ok(Some([0.0; 3]))
        }
        fn nutation_rate_tt(&self, _: f64) -> Result<Option<(f64, f64)>> {
            self.events.borrow_mut().push("rate");
            if self.fail_rate.get() {
                return Err(Error::MathDomain { value: 12.0 });
            }
            Ok(self.stored_rate.get().then_some((1e-7, 0.0)))
        }
        fn nutation_smooth_interval_tt(&self, _: f64) -> Result<Option<(f64, f64)>> {
            self.events.borrow_mut().push("support");
            if self.fail_support.get() {
                return Err(Error::MathDomain { value: 13.0 });
            }
            Ok(None)
        }
    }

    #[test]
    fn analytic_retention_preserves_actual_feed_rate_and_fallback_failures() {
        let provider = LinearProvider {
            events: std::cell::RefCell::new(Vec::new()),
            fail_rate: std::cell::Cell::new(false),
            fail_support: std::cell::Cell::new(false),
            stored_rate: std::cell::Cell::new(true),
        };
        let fail_angle = std::cell::Cell::new(false);
        let angle = |tt: f64| {
            provider.events.borrow_mut().push("angle");
            // Reentrant pure-frame use cannot deadlock on a lock held while
            // calling the actual feed. This also exercises a distinct key.
            analytic_output_frame_words_with_nutation(
                FRAME_J2000,
                FLG_J2000 | FLG_EQUATORIAL,
                (0.0, 0.0),
                (0.0, 0.0),
            )?;
            if fail_angle.get() {
                return Err(Error::MathDomain { value: 11.0 });
            }
            Ok((1e-5 + (tt - FRAME_J2000) * 1e-7, -2e-5))
        };
        let serving = SealedServing::new_stored(&provider, &angle);
        let calculate = || serving_compute_tt(&serving, FRAME_J2000, MARS, FLG_SPEED | FLG_NOGDEFL);
        let first = calculate().unwrap();
        let first_events = provider.events.borrow().clone();
        assert!(first_events.contains(&"rate")); // This really enters the analytic stored-rate route.
        provider.events.borrow_mut().clear();
        let second = calculate().unwrap();
        assert_eq!(first.0.map(f64::to_bits), second.0.map(f64::to_bits));
        assert_eq!(first.1, second.1);
        assert_eq!(*provider.events.borrow(), first_events);
        fail_angle.set(true);
        provider.events.borrow_mut().clear();
        assert!(matches!(
            calculate(),
            Err(Error::MathDomain { value: 11.0 })
        ));
        assert_eq!(*provider.events.borrow(), vec!["angle"]);
        fail_angle.set(false);
        provider.fail_rate.set(true);
        provider.events.borrow_mut().clear();
        assert!(matches!(
            calculate(),
            Err(Error::MathDomain { value: 12.0 })
        ));
        assert_eq!(provider.events.borrow().last(), Some(&"rate"));
        provider.fail_rate.set(false);
        provider.stored_rate.set(false);
        let fallback = calculate().unwrap();
        provider.events.borrow_mut().clear();
        let repeated = calculate().unwrap();
        assert_eq!(fallback.0.map(f64::to_bits), repeated.0.map(f64::to_bits));
        assert!(provider.events.borrow().contains(&"support"));
        provider.fail_support.set(true);
        provider.events.borrow_mut().clear();
        assert!(matches!(
            calculate(),
            Err(Error::MathDomain { value: 13.0 })
        ));
        assert_eq!(provider.events.borrow().last(), Some(&"support"));
    }
}
