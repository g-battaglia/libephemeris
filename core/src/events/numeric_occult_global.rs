//! Native finite-disc geometry and shared scalar occultation searches.
use super::*;

/// Target/Moon longitude conjunction at UT1 Julian date `start`.
/// Data-range exhaustion returns None; all other errors propagate.
/// Angles are degrees. The bounded Newton walk uses the shared lunar rate.
pub(super) fn source_conjunction(
    body: &mut impl FnMut(f64) -> Result<(f64, f64)>,
    moon: &mut impl FnMut(f64) -> Result<(f64, f64)>,
    start: f64,
    direction: i32,
) -> Result<Option<(f64, f64, f64)>> {
    source_conjunction_with_star(body, moon, start, direction, None)
}

/// Shared conjunction search with an optional named-star latitude guard.
/// The target is sampled before the Moon. Initial data-range exhaustion is
/// not an event; all later failures propagate. Inputs are degrees and UT1 JD.
pub(super) fn source_conjunction_with_star(
    body: &mut impl FnMut(f64) -> Result<(f64, f64)>,
    moon: &mut impl FnMut(f64) -> Result<(f64, f64)>,
    start: f64,
    direction: i32,
    star_label: Option<&str>,
) -> Result<Option<(f64, f64, f64)>> {
    seek_moon_conjunction(body, moon, start, direction, star_label.is_some())
}

/// Geocentric finite-disc overlap residual in degrees.
/// Barycentric-relative target and Moon vectors and the target radius are AU.
/// Zero distance, non-finite data, overflow and an impossible angular radius
/// return ordinary typed errors. A negative result means overlapping discs.
pub(crate) fn numeric_occult_separation_source(
    body: [f64; 3],
    moon: [f64; 3],
    radius_au: f64,
) -> Result<f64> {
    if !radius_au.is_finite() {
        return Err(Error::NonFiniteArgument { value: radius_au });
    }
    let db = (finite_square(body[0])? + finite_square(body[1])? + finite_square(body[2])?).sqrt();
    let dm = (finite_square(moon[0])? + finite_square(moon[1])? + finite_square(moon[2])?).sqrt();
    let denominator = db * dm;
    if denominator == 0.0 {
        return Err(Error::RationalDivisionByZero);
    }
    let numerator = body[0] * moon[0] + body[1] * moon[1] + body[2] * moon[2];
    let separation = (numerator / denominator)
        .clamp(-1.0, 1.0)
        .acos()
        .to_degrees();
    let lunar_radius = finite_asin(ECL_RMOON_AU / dm)?.to_degrees();
    let body_radius = if radius_au != 0.0 {
        finite_asin(radius_au / db)?.to_degrees()
    } else {
        0.0
    };
    Ok(separation - (lunar_radius + body_radius))
}

fn finite_square(value: f64) -> Result<f64> {
    if !value.is_finite() {
        return Err(Error::NonFiniteArgument { value });
    }
    let square = value * value;
    if !square.is_finite() {
        return Err(Error::NumericalOverflow {
            operation: "cartesian_distance",
            message: "squared Cartesian distance overflowed".into(),
        });
    }
    Ok(square)
}

fn finite_asin(value: f64) -> Result<f64> {
    if !value.is_finite() || !(-1.0..=1.0).contains(&value) {
        return Err(Error::MathDomain { value });
    }
    Ok(value.asin())
}

#[cfg(test)]
#[path = "numeric_occult_global/tests.rs"]
mod tests;
