//! Solar path sampling over independent acquisition channels. Sampling uses
//! inclusive UT1 endpoints and repeated binary64 addition, without grid repair.
use super::*;

/// Staged width acquisition. Main/Time/body lookup failures propagate before
/// the observer arm. Only the observer and apparent-ray arm may return None;
/// altitude, distances and separation are read later outside that catch.
pub(crate) trait SolarPathWidthSource {
    /// Owned apparent Sun/Moon rays and their observer rotation, in AU.
    type Sample;
    /// Public where result, including its consumed local attributes stage.
    fn width_where(&self, jd: f64) -> Result<(i32, f64, f64)>;
    /// Retain a true WGS84 height0 observer and Sun-then-Moon apparent rays.
    fn width_acquire(&self, jd: f64, lon: f64, lat: f64) -> Result<Option<Self::Sample>>;
    /// Unrefracted Sun altitude in degrees at the retained bundled Time.
    fn width_altitude(&self, sample: &Self::Sample) -> Result<f64>;
    /// Observer-relative Sun then Moon distances, in AU.
    fn width_distances(&self, sample: &Self::Sample) -> (f64, f64);
    /// Apparent ICRS angular separation in degrees.
    fn width_separation(&self, sample: &Self::Sample) -> f64;
}

/// Total/annular cone width in kilometres at UT1. Both site coordinates must
/// be present to select the explicit-site branch; either missing coordinate
/// selects public where. Noncentral where results, caught observer failures,
/// a Sun at/below the horizon, separated discs and an exhausted umbra return0.
/// No eager epoch/site or post-derived finite guard is inserted. The cone is
/// projected by sin(altitude), floored at0.05, and capped at1500km. This local
/// vertical-plane model loses accuracy toward grazing sunrise/sunset geometry.
pub(crate) fn solar_path_width_source_via(
    eph: &impl SolarPathWidthSource,
    jd: f64,
    lat: Option<f64>,
    lon: Option<f64>,
) -> Result<f64> {
    Ok(solar_path_width_captured_source_via(eph, jd, lat, lon)?.value)
}

/// Calculate cone width from the apparent Sun and Moon rays.
pub(crate) fn solar_path_width_captured_source_via(
    eph: &impl SolarPathWidthSource,
    jd: f64,
    lat: Option<f64>,
    lon: Option<f64>,
) -> Result<crate::engine_api::CapturedSolarPathWidth> {
    use crate::engine_api::CapturedSolarPathWidth;
    let (lon, lat) = match (lat, lon) {
        (Some(lat), Some(lon)) => (lon, lat),
        _ => {
            let (flags, lon, lat) = eph.width_where(jd)?;
            if flags & ECL_CENTRAL == 0 {
                return Ok(CapturedSolarPathWidth { value: 0.0 });
            }
            (lon, lat)
        }
    };
    let Some(sample) = eph.width_acquire(jd, lon, lat)? else {
        return Ok(CapturedSolarPathWidth { value: 0.0 });
    };
    let altitude = eph.width_altitude(&sample)?;
    if altitude <= 0.0 {
        return Ok(CapturedSolarPathWidth { value: 0.0 });
    }
    let (sun_au, moon_au) = eph.width_distances(&sample);
    solar_path_width_geometry(altitude, sun_au, moon_au, &|| {
        Ok(eph.width_separation(&sample))
    })
}

/// Complete the cone-width calculation from acquired apparent geometry.
/// Altitude and separation are degrees, distances AU and output kilometres.
/// The caller has already applied the horizon gate. Separation is acquired
/// only after the ordered angular-radius divisions; its refusal stays outside
/// the observer/apparent catch. No observer, resource or clock is constructed.
pub(crate) fn solar_path_width_geometry(
    altitude: f64,
    sun_au: f64,
    moon_au: f64,
    separation: &dyn Fn() -> Result<f64>,
) -> Result<crate::engine_api::CapturedSolarPathWidth> {
    if !altitude.is_finite() || !sun_au.is_finite() || !moon_au.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: altitude + sun_au + moon_au,
        });
    }
    if sun_au <= 0.0 || moon_au <= 0.0 {
        return Err(Error::MathDomain {
            value: sun_au.min(moon_au),
        });
    }
    let sun_km = sun_au * ECL_AU_KM;
    let moon_km = moon_au * ECL_AU_KM;
    let sun_radius = (959.63 / 3600.0) / sun_au;
    let moon_radius = (932.56 / 3600.0) * (0.002569 / moon_au);
    let ratio = moon_radius / sun_radius;
    if separation()? > (moon_radius - sun_radius).abs() {
        return Ok(crate::engine_api::CapturedSolarPathWidth { value: 0.0 });
    }
    let alpha = ((PATH_WIDTH_SUN_RADIUS_KM - PATH_WIDTH_MOON_RADIUS_KM) / sun_km).atan();
    let projected = moon_km * alpha.tan();
    let radius = if ratio >= 1.0 {
        let radius = PATH_WIDTH_MOON_RADIUS_KM - projected;
        if radius <= 0.0 {
            return Ok(crate::engine_api::CapturedSolarPathWidth { value: 0.0 });
        }
        radius
    } else {
        (projected - PATH_WIDTH_MOON_RADIUS_KM).abs()
    };
    let sin_alt = altitude.to_radians().sin().max(PATH_WIDTH_SIN_ALT_FLOOR);
    let width = (2.0 * radius / sin_alt).clamp(0.0, PATH_WIDTH_MAX_KM);
    Ok(crate::engine_api::CapturedSolarPathWidth { value: width })
}

/// Limit walk using NaN-preserving bounded angular geometry. The existing
/// checked finite point helper and sealed walks retain their own contracts.
/// Besselian x/y are Earth radii, d/mu degrees and l2 Earth radii. North offsets
/// y by +abs(l2), south by -abs(l2); omitted samples append no components.
pub(crate) fn solar_limit_source_via(
    eph: &impl EclipseLimitEphemeris,
    south: bool,
    start: f64,
    end: f64,
    minutes: f64,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    let name = if south {
        "calc_eclipse_southern_limit"
    } else {
        "calc_eclipse_northern_limit"
    };
    let step = eclipse_sampling_step_days(start, end, minutes, name)?;
    let (mut times, mut lats, mut lons) = (Vec::new(), Vec::new(), Vec::new());
    let mut jd = start;
    while jd <= end {
        let (x, y, d, mu, l2) = eph.limit_besselian(jd)?;
        if let Some((lat, lon)) = solar_limit_point_source(x, y, d, mu, l2, south)? {
            times.push(jd);
            lats.push(lat);
            lons.push(lon);
        }
        jd += step;
    }
    Ok((times, lats, lons))
}

/// One boundary point with spherical shadow-plane geometry followed by the
/// WGS84 geocentric-to-geodetic atan(tan(phi)/(1-f)^2) conversion. Longitude
/// uses floored modulo into [-180,180). Infinite d is a trig-domain failure;
/// NaN declination remains NaN. There is no interpolation between samples.
fn solar_limit_point_source(
    x: f64,
    y: f64,
    d: f64,
    mu: f64,
    l2: f64,
    south: bool,
) -> Result<Option<(f64, f64)>> {
    let radius = l2.abs();
    let gamma = (x * x + y * y).sqrt();
    if !(gamma < 1.0 + PATH_LIMIT_FLATTENING + radius && radius > PATH_LIMIT_UMBRA_FLOOR_ER) {
        return Ok(None);
    }
    let dr = d.to_radians();
    // radians(mu) is consumed even though the converted value is discarded.
    let _ = mu.to_radians();
    let edge = if south { y - radius } else { y + radius };
    let edge_gamma = (x * x + edge * edge).sqrt();
    if edge_gamma.partial_cmp(&(1.0 + PATH_LIMIT_FLATTENING)) != Some(std::cmp::Ordering::Less) {
        return Ok(None);
    }
    let z = if edge_gamma > PATH_LIMIT_SURFACE_GAMMA {
        0.0
    } else {
        let v = 1.0 - edge_gamma * edge_gamma;
        (if v > 0.0 { v } else { 0.0 }).sqrt()
    };
    if dr.is_infinite() {
        return Err(sine_domain(dr));
    }
    let sin_d = dr.sin();
    let cos_d = dr.cos();
    let v = edge * cos_d + z * sin_d;
    let v = v.clamp(-1.0, 1.0);
    let lat = v.asin().to_degrees();
    let lon_offset = if lat.to_radians().cos() > PATH_LIMIT_COS_LAT_FLOOR {
        x.atan2(z * cos_d - edge * sin_d).to_degrees()
    } else {
        0.0
    };
    let lon = -mu + lon_offset;
    let geodetic = (lat.to_radians().tan() / (1.0 - PATH_LIMIT_FLATTENING).powi(2))
        .atan()
        .to_degrees();
    Ok(Some((geodetic, normalize_path_longitude_deg(lon))))
}

/// Retain the infinite radian operand of an actual sine domain failure.
/// The binding acquires genuine interpreter arguments from the held intrinsic;
/// no message, category fold, global error slot or source result is invented.
fn sine_domain(input: f64) -> Error {
    Error::MathDomain { value: input }
}

#[cfg(test)]
#[path = "solar_path/tests.rs"]
mod tests;
