//! Checked coordinate-rotation layouts over the existing engine kernels.
use crate::errors::Error;

/// Rotate a spherical position. Angles and obliquity are degrees; distance
/// passes through in the caller's length unit. Positive obliquity converts
/// equatorial to ecliptic, negative converts back. Non-finite legs are errors.
pub fn cotrans(
    longitude_or_ra: f64,
    latitude_or_dec: f64,
    distance: f64,
    eps: f64,
) -> Result<[f64; 3], Error> {
    for value in [longitude_or_ra, latitude_or_dec, distance, eps] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    let (longitude, latitude) = crate::positions::cotrans(longitude_or_ra, latitude_or_dec, eps)?;
    Ok([longitude, latitude, distance])
}

/// Rotate a spherical position and its daily rates using the existing analytic
/// derivative. Angles/angle rates are degrees/degrees per day; distance and its
/// rate retain their input length unit. Positive obliquity converts equatorial
/// to ecliptic, negative converts back. Non-finite legs are errors; no numerical
/// difference stencil or new rotation algorithm is introduced.
#[allow(clippy::too_many_arguments)]
pub fn cotrans_sp(
    longitude_or_ra: f64,
    latitude_or_dec: f64,
    distance: f64,
    longitude_or_ra_speed: f64,
    latitude_or_dec_speed: f64,
    distance_speed: f64,
    eps: f64,
) -> Result<([f64; 3], [f64; 3]), Error> {
    crate::utils::cotrans_sp(
        longitude_or_ra,
        latitude_or_dec,
        distance,
        longitude_or_ra_speed,
        latitude_or_dec_speed,
        distance_speed,
        eps,
    )
}
